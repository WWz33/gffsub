# Library API

`libgffsub_core.a` exposes the same functionality the `gffsub` binary uses. Link it from C++17 code to subset, query, or convert GFF3/GTF annotations in-process.

## Linking

```makefile
CXXFLAGS = -std=c++17 -I/path/to/gffsub/src
LDFLAGS = -L/path/to/gffsub -lgffsub_core -pthread
```

Include the specific headers from `src/`. `src/gff3.hpp` is an umbrella that pulls in the common set (`annotation`, `filter`, `output`, `parser`, `record`, `region`).

## Data model

```cpp
// src/record.hpp
struct GffRecord {
    std::string seqid, source, type;
    int64_t start, end;
    std::optional<double> score;
    std::string score_raw;
    char strand, phase;
    std::string attr_raw;                 // raw column 9
    std::optional<std::string> id;         // ID= or synthesized from gene_id/transcript_id
    std::optional<std::string> parent_id; // Parent= or synthesized
    std::optional<std::string> gene_id;   // GTF gene_id
    std::optional<std::string> transcript_id; // GTF transcript_id
    int line_idx;
    bool kept;
    InputFormat src_fmt;                   // GFF3, GTF, or BED
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

`GffData` holds records in input order.

## Parsing

```cpp
// src/parser.hpp
int parse_file(const std::string& filename, GffData& data, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);
InputFormat infer_input_format(const std::string& path);
```

`parse_file` reads a GFF3, GTF, or BED file, populates `data`, returns 0 on success or -1 on open failure. `infer_input_format` maps `.gtf`/`.GTF` to GTF, `.bed`/`.BED` to BED, everything else to GFF3. `parse_attributes` splits a GFF3 column-9 string into key to list-of-values, URL-decoding each value. Comma is the value-list separator (GFF3 spec), `%2C` decodes to a literal comma.

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

> **Migration (IdIndex removed).** The standalone `IdIndex` class and the `IdIndex&` parameter of `parse_file` are gone. The ID-to-record index is now built inside `AnnotationIndex`. Replace `parse_file(path, data, idx, fmt)` with `parse_file(path, data, fmt)`, and `idx.lookup(id)` with `AnnotationIndex::from_data(std::move(data)).find_by_id(id)`.

`from_file` infers format from the extension (`.gtf` -> GTF, `.bed` -> BED, else GFF3) and builds the index. `from_gff3` forces GFF3 parsing. `from_data` takes an already-parsed `GffData`. The index stores ID to record, parent/child links, gene-name lookup keys (ID, gene_id, Name, locus_tag, Alias, Dbxref), and attribute indices.

`find_gene` searches gene-type records by any of the naming keys. `gene_model` returns the gene plus all transcripts and their children. `find_all_by_id` returns every line sharing an ID (multi-line CDS, discontinuous features).

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
    std::string summary_format; // "tsv", "json", or empty
};

struct QueryResult {
    GffData records;
    std::vector<SummaryRow> summary_rows;
    bool emit_summary = false;
};

QueryResult query(const AnnotationIndex& index, const QueryParams& params);
void print_query_result(std::ostream& out, const QueryResult& result, const QueryParams& params);
```

`query` resolves the selectors (ids, name, attr_filters, nearest_region), expands along the hierarchy per the include flags, and returns matched records. When `summary_format` is set, `summary_rows` is populated and `emit_summary` is true. `print_query_result` writes GFF3 or summary output based on `params.summary_format`.

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

Finds the feature by ID (falls back to gene-name lookup), builds a window around it, returns all overlapping records. Without `strand_aware`, upstream extends left of the start and downstream extends right of the end. With `strand_aware` on a minus-strand feature, the extensions swap. Start below 1 is clamped to 1.

## Subset

```cpp
// src/subset.hpp
struct SubsetParams {
    std::optional<Region> region;
    std::string bed_file;
    std::string seqid_filter;  // ^-prefixed for exclude
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

Applies region, BED, seqid, source, score, strand, phase, grep, expr, and feature filters in sequence. Modifies `data` in place. `longest` selects one isoform per gene (see [Longest Isoform Selection](longest-isoform.md)). `threads` parallelizes longest selection across chromosomes.

Individual filters are also callable directly (`src/filter.hpp`, `src/selector_filter.hpp`):

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

## Expression filters

```cpp
// src/expr_parser.hpp
struct ExprFilter {
    std::string field;
    ExprOp op;      // Equal, NotEqual, Regex, NotRegex, Less, LessEqual, Greater, GreaterEqual
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

Parse expression strings like `length >= 100 && strand == "+"` into an `ExprNode` tree. Compile regexes once before passing to `filter_by_expr`. See [Expression Filtering](expression-filtering.md) for field names and operators.

## Output

```cpp
// src/output.hpp
void print_gff3(std::ostream& out, const GffData& data);
void print_gtf3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);  // GTF2 or GTF3
void print_bed(std::ostream& out, const GffData& data);
```

`print_gtf` with `OutputFormat::GTF2` emits GTF2.2; with `GTF3` emits GTF2.2.1 with mRNA renamed to transcript. GFF3 output preserves column 9 as-is unless records came from GTF input (`src_fmt == GTF`), in which case column 9 is rewritten as `tag=value` with synthesized `ID=`/`Parent=` and URL escaping.

## Summary

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
    std::string status;  // "found" or "not_found"
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

Counts in `make_summary_row` come from `descendants_of` the matched feature. `extract_output_attrs` pulls selected column-9 attributes, resolving `gene_id`/`transcript_id` from the record fields first, then GFF3 attributes, then GTF quoted values.

## Region helpers

```cpp
// src/region.hpp
struct Region { std::string seqid; int64_t start; int64_t end; };
struct BedRegion { std::string seqid; int64_t start; int64_t end; };

std::optional<Region> parse_region(std::string_view region_str);  // "chr1:100-1000"
BedRegion to_bed_region(const GffRecord& rec);      // 1-based GFF3 to 0-based BED
Region from_bed_region(const BedRegion& region);    // 0-based BED to 1-based GFF3
Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware);
```

## GTF helpers

```cpp
// src/gtf_parser.hpp
std::optional<std::string> extract_quoted_value(const std::string& attrs, const std::string& key);
std::string gtf_attrs_to_gff3(const GffRecord& rec);
void apply_gtf_attributes(GffRecord& rec);
```

`extract_quoted_value` pulls a `key "value";` pair from a GTF column-9 string, handling escaped quotes. `gtf_attrs_to_gff3` rewrites GTF column 9 as GFF3 `tag=value` pairs, synthesizing `ID=`/`Parent=` from `gene_id`/`transcript_id` and URL-escaping values. `apply_gtf_attributes` populates `rec.gene_id` and `rec.transcript_id` from column 9.

## Example

Build an index, query by gene name with children, and print GTF:

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

Subset by region then keep longest isoform:

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

Window around a feature:

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

Compile and link:

```bash
g++ -std=c++17 -Isrc main.cpp -L. -lgffsub_core -pthread -o example
```
