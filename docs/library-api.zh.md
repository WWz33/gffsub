# Library API

<!-- I18N:START -->

[English](./library-api.md) | **中文**

<!-- I18N:END -->

`libgffsub_core.a` 暴露与 `gffsub` 二进制相同的功能。从 C++17 代码中链接它，可在进程内对 GFF3/GTF 注释做 subset、query 或转换。

## 链接

```makefile
CXXFLAGS = -std=c++17 -I/path/to/gffsub/src
LDFLAGS = -L/path/to/gffsub -lgffsub_core -pthread
```

从 `src/` 引入具体头文件。`src/gff3.hpp` 是 umbrella 头文件，包含常用集合（`annotation`、`filter`、`output`、`parser`、`record`、`region`）。

## 数据模型

```cpp
// src/record.hpp
struct GffRecord {
    std::string seqid, source, type;
    int64_t start, end;
    std::optional<double> score;
    std::string score_raw;
    char strand, phase;
    std::string attr_raw;                 // 原始第 9 列
    std::optional<std::string> id;         // ID= 或从 gene_id/transcript_id 合成
    std::optional<std::string> parent_id; // Parent= 或合成
    std::optional<std::string> gene_id;   // GTF gene_id
    std::optional<std::string> transcript_id; // GTF transcript_id
    int line_idx;
    bool kept;
    InputFormat src_fmt;                   // GFF3、GTF 或 BED
};

class GffData {
public:
    std::vector<GffRecord> records;
    void append(const GffRecord&);
    size_t size() const;
    iterator begin(), end();
    void clear();
    void reserve(size_t n);
};
```

`GffData` 按输入顺序保存记录。

## 解析

```cpp
// src/parser.hpp
int parse_file(const std::string& filename, GffData& data, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);
InputFormat infer_input_format(const std::string& path);
```

`parse_file` 读取 GFF3、GTF 或 BED 文件，填充 `data`，成功返回 0，打开失败返回 -1。`infer_input_format` 把 `.gtf`/`.GTF` 映射为 GTF，`.bed`/`.BED` 映射为 BED，其余映射为 GFF3。`parse_attributes` 把 GFF3 第 9 列字符串拆成 key 到 value 列表的映射，对每个 value 做 URL 解码。逗号是 value 列表分隔符（GFF3 规范），`%2C` 解码为字面逗号。

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

> **迁移说明（IdIndex 已移除）。** 独立的 `IdIndex` 类以及 `parse_file` 的 `IdIndex&` 参数已删除。ID 到 record 的索引现在在 `AnnotationIndex` 内部构建。把 `parse_file(path, data, idx, fmt)` 替换为 `parse_file(path, data, fmt)`，把 `idx.lookup(id)` 替换为 `AnnotationIndex::from_data(std::move(data)).find_by_id(id)`。

`from_file` 根据扩展名推断格式（`.gtf` -> GTF，`.bed` -> BED，其余 GFF3）并构建索引。`from_gff3` 强制按 GFF3 解析。`from_data` 接收已解析的 `GffData`。索引存储 ID 到 record、parent/child 链接、gene name lookup key（ID、gene_id、Name、locus_tag、Alias、Dbxref）以及属性索引。

`find_gene` 按任一命名 key 在 gene 类型记录中查找。`gene_model` 返回 gene 及其所有 transcript 和子节点。`find_all_by_id` 返回共享同一 ID 的所有行（多行 CDS、不连续 feature）。

## Query

```cpp
// src/query.hpp
struct QueryParams {
    std::vector<std::string> ids;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters; // --where KEY=VALUE
    std::optional<Region> region;
    std::optional<Region> nearest_region;
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;
    std::string feature_type;
    bool apply_type_filter = true;
    std::vector<std::string> output_attrs;
    std::string summary_format; // "tsv"、"json" 或空
};

struct QueryResult {
    GffData records;
    std::vector<SummaryRow> summary_rows;
    bool emit_summary = false;
};

QueryResult query(const AnnotationIndex& index, const QueryParams& params);
void print_query_result(std::ostream& out, const QueryResult& result, const QueryParams& params);
```

`query` 解析 selector（ids、name、attr_filters、nearest_region），按 include 标志沿层级扩展，返回匹配记录。设置 `summary_format` 时填充 `summary_rows` 并置 `emit_summary` 为 true。`print_query_result` 根据 `params.summary_format` 输出 GFF3 或摘要。

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

按 ID 查找 feature（回退到 gene name 查找），在其周围构建窗口，返回所有重叠记录。不设 `strand_aware` 时，upstream 向 start 左侧延伸，downstream 向 end 右侧延伸。对 minus strand feature 开启 `strand_aware` 时，延伸方向互换。start 小于 1 时截断为 1。

## Subset

```cpp
// src/subset.hpp
struct SubsetParams {
    std::optional<Region> region;
    std::string bed_file;
    std::string seqid_filter;  // ^ 前缀表示排除
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

依次应用 region、BED、seqid、source、score、strand、phase、grep、expr 和 feature 筛选。原地修改 `data`。`longest` 每 gene 选一个 isoform（见 [Longest Isoform Selection](longest-isoform.zh.md)）。`threads` 按染色体并行执行 longest 选择。

各筛选函数也可单独调用（`src/filter.hpp`、`src/selector_filter.hpp`）：

```cpp
void filter_by_region(GffData& data, const Region& region);
void filter_by_regions_from_file(GffData& data, const std::string& bed_file);
void filter_by_feature(GffData& data, std::string_view feature_type);
void filter_by_seqid(GffData& data, const std::unordered_set<std::string>& seqids, bool exclude);
void filter_by_source(GffData& data, std::string_view source);
void filter_by_score(GffData& data, std::optional<double> score);
void filter_by_strand(GffData& data, char strand);
void filter_by_phase(GffData& data, char phase);
void filter_longest_isoform(GffData& data, std::string_view feature_type, size_t num_threads = 1);

void filter_by_grep(GffData& data, const std::vector<GrepFilter>& filters, bool invert);
void filter_by_expr(GffData& data, const std::vector<ExprNode>& filters, bool include);
```

## 表达式筛选

```cpp
// src/expr_parser.hpp
struct ExprFilter {
    std::string field;
    ExprOp op;      // Equal、NotEqual、Regex、NotRegex、Less、LessEqual、Greater、GreaterEqual
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

把 `length >= 100 && strand == "+"` 之类的表达式字符串解析为 `ExprNode` 树。传入 `filter_by_expr` 前先编译一次正则。字段名和操作符见 [Expression Filtering](expression-filtering.zh.md)。

## 输出

```cpp
// src/output.hpp
void print_gff3(std::ostream& out, const GffData& data);
void print_gtf3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);  // GTF2 或 GTF3
void print_bed(std::ostream& out, const GffData& data);
```

`print_gtf` 配合 `OutputFormat::GTF2` 输出 GTF2.2；配合 `GTF3` 输出 GTF2.2.1 并把 mRNA 改名为 transcript。GFF3 输出原样保留第 9 列，除非记录来自 GTF 输入（`src_fmt == GTF`），此时第 9 列改写为 `tag=value`，合成 `ID=`/`Parent=` 并做 URL 转义。

## 摘要

```cpp
// src/query_summary.hpp
struct SummaryRow {
    std::string query_id, matched_id, matched_by;
    std::string seqid;
    int64_t start, end;
    char strand;
    std::string type, parent_id;
    size_t child_count, transcript_count, exon_count;
    int64_t cds_length;
    std::string status;  // "found" 或 "not_found"
    std::vector<std::string> attrs;
};

SummaryRow make_summary_row(const AnnotationIndex& index,
                            const std::string& query_id,
                            const std::string& matched_by,
                            const GffRecord& rec);
SummaryRow make_not_found_row(const std::string& query_id, const std::string& matched_by);
std::vector<std::string> extract_output_attrs(const GffRecord& rec, const std::vector<std::string>& keys);
void print_summary_tsv(std::ostream& out, const std::vector<SummaryRow>& rows,
                       const std::vector<std::string>& output_attrs);
void print_summary_json(std::ostream& out, const std::vector<SummaryRow>& rows,
                        const std::vector<std::string>& output_attrs);
```

`make_summary_row` 中的计数来自命中 feature 的 `descendants_of`。`extract_output_attrs` 抽取选定的第 9 列属性，优先从 record 字段解析 `gene_id`/`transcript_id`，其次 GFF3 属性，最后 GTF 带引号的值。

## Region 辅助

```cpp
// src/region.hpp
struct Region { std::string seqid; int64_t start; int64_t end; };
struct BedRegion { std::string seqid; int64_t start; int64_t end; };

std::optional<Region> parse_region(std::string_view region_str);  // "chr1:100-1000"
BedRegion to_bed_region(const GffRecord& rec);      // 1-based GFF3 转 0-based BED
Region from_bed_region(const BedRegion& region);    // 0-based BED 转 1-based GFF3
Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware);
```

## GTF 辅助

```cpp
// src/gtf_parser.hpp
std::optional<std::string> extract_quoted_value(const std::string& attrs, const std::string& key);
std::string gtf_attrs_to_gff3(const GffRecord& rec);
void apply_gtf_attributes(GffRecord& rec);
```

`extract_quoted_value` 从 GTF 第 9 列字符串中取出 `key "value";` 键值对，处理转义引号。`gtf_attrs_to_gff3` 把 GTF 第 9 列改写为 GFF3 `tag=value` 对，从 `gene_id`/`transcript_id` 合成 `ID=`/`Parent=` 并对 value 做 URL 转义。`apply_gtf_attributes` 从第 9 列填充 `rec.gene_id` 和 `rec.transcript_id`。

## 示例

构建索引，按 gene name 连同子节点查询，并输出 GTF：

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
    return 0;
}
```

按 region subset 后保留最长 isoform：

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
    return 0;
}
```

feature 周围的窗口：

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
    return 0;
}
```

编译并链接：

```bash
g++ -std=c++17 -Isrc main.cpp -L. -lgffsub_core -pthread -o example
```
