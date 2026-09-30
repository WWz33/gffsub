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
    // Backing storage for record string_views. Either an owned buffer filled
    // once by parse_file (never grows afterwards, so views stay valid) or an
    // mmap'd file exposed through content(); never both.
    std::string buffer;

    // Move-only: copying would leave record views pointing into the source.
    GffData() = default;
    GffData(const GffData&) = delete;
    GffData& operator=(const GffData&) = delete;
    GffData(GffData&& other) noexcept { move_from(std::move(other)); }
    GffData& operator=(GffData&& other) noexcept {
        if (this != &other) {
            release_mapping();
            move_from(std::move(other));
        }
        return *this;
    }
    ~GffData() { release_mapping(); }

    // Byte range the record views point into: the mmap when present, the
    // owned buffer otherwise.
    std::string_view content() const {
        if (mapped_data_ != nullptr) {
            return std::string_view{static_cast<const char*>(mapped_data_), mapped_size_};
        }
        return std::string_view{buffer};
    }

    void set_mapping(void* data, size_t size) {
        release_mapping();
        mapped_data_ = data;
        mapped_size_ = size;
    }

    // Copy the owned storage (records, directives, buffer). Never copies a
    // mapping: record views keep pointing at the source's storage, so the
    // source must outlive the copy — the same borrowed-storage contract as
    // AnnotationIndex::from_data().
    void copy_storage_from(const GffData& other) {
        release_mapping();
        buffer = other.buffer;
        records = other.records;
        directives = other.directives;
    }

    void append(const GffRecord& rec) { records.push_back(rec); }
    void append(GffRecord&& rec) { records.push_back(std::move(rec)); }
    auto size() const { return records.size(); }
    auto begin() { return records.begin(); }
    auto end() { return records.end(); }
    auto begin() const { return records.begin(); }
    auto end() const { return records.end(); }
    void clear() {
        records.clear();
        directives.clear();
        buffer.clear();
        release_mapping();
    }
    void reserve(size_t n) { records.reserve(n); }

private:
    void move_from(GffData&& other) {
        records = std::move(other.records);
        directives = std::move(other.directives);
        buffer = std::move(other.buffer);
        mapped_data_ = other.mapped_data_;
        mapped_size_ = other.mapped_size_;
        other.mapped_data_ = nullptr;
        other.mapped_size_ = 0;
    }

    void release_mapping();

    void* mapped_data_ = nullptr;
    size_t mapped_size_ = 0;
};

}  // namespace gffsub

#endif  // GFFSUB_RECORD_HPP
