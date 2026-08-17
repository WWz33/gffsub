# Library API

<!-- I18N:START -->

[English](./library-api.md) | **中文**

<!-- I18N:END -->

`libgffsub_core.a` 提供 `gffsub` 二进制使用的全部函数。C++17 链接。

## 链接

```makefile
CXXFLAGS = -std=c++17 -I/path/to/gffsub/src
LDFLAGS = -L/path/to/gffsub -lgffsub_core -pthread
```

## 数据模型

```cpp
// src/record.hpp
struct GffRecord {
    std::string seqid, source, type;
    FeatureClass feat_class = FeatureClass::Unknown;
    int64_t start, end;
    std::optional<double> score;
    std::string score_raw;
    char strand, phase;
    std::string attr_raw;
    std::optional<std::string> id;
    std::optional<std::string> parent_id;
    std::optional<std::string> gene_id;
    std::optional<std::string> transcript_id;
    int line_idx;
    bool kept;
    InputFormat src_fmt;
};

class GffData {
public:
    std::vector<GffRecord> records;
    std::vector<std::string> directives;
    void append(const GffRecord&);
    void append(GffRecord&&);
    size_t size() const;
    iterator begin(), end();
    void clear();
    void reserve(size_t n);
};
```

## 解析

```cpp
// src/parser.hpp
int parse_file(const std::string& filename, GffData& data, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);
InputFormat infer_input_format(const std::string& path);
```

`infer_input_format` 嗅探文件内容：第 9 列含 `=` 和 `;` 判为 GFF3，含 `"` 判为 GTF，3-12 列且 start/end 为整数判为 BED，默认 GFF3。

## AnnotationIndex

```cpp
// src/annotation.hpp
class AnnotationIndex {
public:
    static AnnotationIndex from_file(const std::string& path);
    static AnnotationIndex from_gff3(const std::string& path);
    static AnnotationIndex from_data(GffData data);

    std::optional<GffRecord> find_by_id(std::string_view id) const;
    std::vector<GffRecord> find_all_by_id(std::string_view id) const;
    std::optional<GffRecord> find_gene(std::string_view id) const;
    std::vector<GffRecord> parents_of(std::string_view id) const;
    std::vector<GffRecord> children_of(std::string_view parent_id) const;
    std::vector<GffRecord> descendants_of(std::string_view parent_id) const;
    std::vector<GffRecord> overlap(std::string_view seqid, int64_t start, int64_t end) const;
    std::optional<GffRecord> nearest_gene(std::string_view seqid, int64_t start, int64_t end) const;
    std::vector<GffRecord> with_attribute(std::string_view key, std::string_view value) const;
    std::optional<GeneModel> gene_model(std::string_view id) const;
};

struct GeneModel {
    GffRecord gene;
    std::vector<GffRecord> records;
};
```

`find_gene` 按 ID、gene_id、Name、locus_tag、Alias、Dbxref 查找 gene 类型记录。`find_all_by_id` 返回共享同一 ID 的全部行（多行 CDS）。

## Query

```cpp
// src/query.hpp
struct QueryParams {
    std::vector<std::string> ids;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    std::optional<Region> region;
    std::optional<Region> nearest_region;
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;
    std::string feature_type;
    bool apply_type_filter = true;
};

struct QueryResult {
    GffData records;
};

QueryResult query(const AnnotationIndex& index, const QueryParams& params);
```

## Window

```cpp
// src/window.hpp
struct WindowParams {
    std::string id;
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;
};

GffData window(const AnnotationIndex& index, const WindowParams& params);
```

## Subset

```cpp
// src/subset.hpp
struct SubsetParams {
    std::optional<Region> region;
    std::string bed_file;
    std::string seqid_filter;
    std::string source_filter;
    std::optional<std::optional<double>> score_filter;
    std::optional<char> strand_filter;
    std::optional<char> phase_filter;
    std::string feature;
    bool longest = false;
    size_t threads = 1;
    std::vector<GrepFilter> grep_filters;
    std::vector<ExprNode> include_exprs;
    std::vector<ExprNode> exclude_exprs;
    bool invert_grep = false;
};

void subset(GffData& data, const SubsetParams& params);
```

## 过滤函数

```cpp
// src/filter.hpp
void filter_by_region(GffData& data, const Region& region);
void filter_by_regions_from_file(GffData& data, const std::string& bed_file);
void filter_by_feature(GffData& data, const std::unordered_set<std::string>& types, bool exclude);
void filter_by_seqid(GffData& data, const std::unordered_set<std::string>& seqids, bool exclude);
void filter_by_source(GffData& data, const std::unordered_set<std::string>& sources, bool exclude);
void filter_by_score(GffData& data, std::optional<double> score);
void filter_by_strand(GffData& data, char strand);
void filter_by_phase(GffData& data, char phase);
void filter_longest_isoform(GffData& data, std::string_view feature_type, size_t num_threads = 1);

// src/selector_filter.hpp
void filter_by_grep(GffData& data, const std::vector<GrepFilter>& filters, bool invert);
void filter_by_expr(GffData& data, const std::vector<ExprNode>& filters, bool include);
```

## 表达式过滤

```cpp
// src/expr_parser.hpp
struct ExprFilter {
    std::string field;
    ExprOp op;  // Equal, NotEqual, Regex, NotRegex, Less, LessEqual, Greater, GreaterEqual
    std::string value;
    bool ignore_case = false;
    std::optional<std::regex> compiled;
};

struct ExprNode {
    enum class Kind { Predicate, Not, And, Or };
    Kind kind = Kind::Predicate;
    ExprFilter predicate;
    std::vector<ExprNode> children;
};

bool parse_expr_filters(std::string_view expr, std::vector<ExprNode>& out, std::string& error);
bool compile_filter_regexes(std::vector<GrepFilter>& grep_filters,
                            std::vector<ExprNode>& include_filters,
                            std::vector<ExprNode>& exclude_filters,
                            bool ignore_case, std::string& error);
```

字段名和操作符见 [Expression 过滤](expression-filtering.zh.md)。

## 输出

```cpp
// src/output.hpp
void print_gff3(std::ostream& out, const GffData& data);
void print_gtf3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);  // GTF2 or GTF3
void print_bed(std::ostream& out, const GffData& data);
```

## Summary

```cpp
// src/query_summary.hpp
std::string record_id(const GffRecord& rec);
void print_summary(std::ostream& out, const std::vector<GffRecord>& records);
```

## Region

```cpp
// src/region.hpp
struct Region { std::string seqid; int64_t start, end; };
struct BedRegion { std::string seqid; int64_t start, end; };

std::optional<Region> parse_region(std::string_view region_str);  // "chr1:100-1000"
BedRegion to_bed_region(const GffRecord& rec);      // 1-based GFF3 to 0-based BED
Region from_bed_region(const BedRegion& region);    // 0-based BED to 1-based GFF3
Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware);
```

## GTF 辅助函数

```cpp
// src/gtf_parser.hpp
std::optional<std::string> extract_quoted_value(const std::string& attrs, const std::string& key);
std::string gtf_attrs_to_gff3(const GffRecord& rec);
void apply_gtf_attributes(GffRecord& rec);
```

## 示例

按 gene name 查询并展开子节点，输出 GTF：

```cpp
#include "annotation.hpp"
#include "query.hpp"
#include "output.hpp"

int main() {
    gffsub::AnnotationIndex idx = gffsub::AnnotationIndex::from_gff3("ann.gff3");

    gffsub::QueryParams qp;
    qp.name = "BRCA1";
    qp.include_children = true;
    gffsub::QueryResult res = gffsub::query(idx, qp);

    gffsub::print_gtf(std::cout, res.records, gffsub::OutputFormat::GTF2);
}
```

按 region 过滤后保留最长 isoform：

```cpp
#include "annotation.hpp"
#include "parser.hpp"
#include "subset.hpp"
#include "output.hpp"

int main() {
    gffsub::GffData data;
    gffsub::parse_file("ann.gff3", data, gffsub::InputFormat::GFF3);

    gffsub::SubsetParams sp;
    sp.region = gffsub::parse_region("chr1:1-1000000").value();
    sp.longest = true;
    sp.threads = 4;
    gffsub::subset(data, sp);

    gffsub::print_gff3(std::cout, data);
}
```

feature 周围窗口：

```cpp
#include "annotation.hpp"
#include "window.hpp"
#include "output.hpp"

int main() {
    gffsub::AnnotationIndex idx = gffsub::AnnotationIndex::from_gff3("ann.gff3");

    gffsub::WindowParams wp;
    wp.id = "gene01";
    wp.upstream = 2000;
    wp.downstream = 500;
    wp.strand_aware = true;
    gffsub::GffData res = gffsub::window(idx, wp);

    gffsub::print_gff3(std::cout, res);
}
```

编译链接：

```bash
g++ -std=c++17 -Isrc main.cpp -L. -lgffsub_core -pthread -o example
```
