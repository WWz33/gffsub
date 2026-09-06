#ifndef GFFSUB_RECORD_HPP
#define GFFSUB_RECORD_HPP

#include "feature_types.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace gffsub {

// String fields are views into GffData::buffer (the whole input file read
// once). A GffRecord is valid only while its owning GffData lives and keeps
// its buffer; copying a record copies the views, which stay valid as long as
// the buffer does. BED-sourced records have empty attr_raw; print_gff3
// synthesizes "ID=..." at output time.
struct GffRecord {
    std::string_view seqid;
    std::string_view source;
    std::string_view type;
    FeatureClass feat_class = FeatureClass::Unknown;
    int64_t start = 0;
    int64_t end = 0;
    std::optional<double> score;
    std::string_view score_raw;
    char strand = '.';
    char phase = '.';
    std::string_view attr_raw;
    std::optional<std::string> id;
    std::optional<std::string> parent_id;
    std::optional<std::string> gene_id;
    std::optional<std::string> transcript_id;
    int line_idx = 0;
    bool kept = true;
    // Input format of the record's attr_raw. GTF attr_raw uses `key "value";`
    // which is invalid in GFF3 output; print_gff3 rewrites it when src_fmt==GTF.
    InputFormat src_fmt = InputFormat::GFF3;
};

class GffData {
public:
    std::vector<GffRecord> records;
    std::vector<std::string> directives;
    // Backing storage for record string_views. Filled once by parse_file,
    // never grows afterwards, so views stay valid.
    std::string buffer;

    void append(const GffRecord& rec) { records.push_back(rec); }
    void append(GffRecord&& rec) { records.push_back(std::move(rec)); }
    auto size() const { return records.size(); }
    auto begin() { return records.begin(); }
    auto end() { return records.end(); }
    auto begin() const { return records.begin(); }
    auto end() const { return records.end(); }
    void clear() { records.clear(); directives.clear(); buffer.clear(); }
    void reserve(size_t n) { records.reserve(n); }
};

}  // namespace gffsub

#endif  // GFFSUB_RECORD_HPP
