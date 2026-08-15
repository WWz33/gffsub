#ifndef GFFSUB_RECORD_HPP
#define GFFSUB_RECORD_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

enum class OutputFormat { GFF3, GTF2, GTF3, BED };
enum class InputFormat { GFF3, GTF, BED };

struct GffRecord {
    std::string seqid;
    std::string source;
    std::string type;
    int64_t start = 0;
    int64_t end = 0;
    std::optional<double> score;
    std::string score_raw;
    char strand = '.';
    char phase = '.';
    std::string attr_raw;
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

    void append(const GffRecord& rec) { records.push_back(rec); }
    auto size() const { return records.size(); }
    auto begin() { return records.begin(); }
    auto end() { return records.end(); }
    auto begin() const { return records.begin(); }
    auto end() const { return records.end(); }
    void clear() { records.clear(); }
    void reserve(size_t n) { records.reserve(n); }
};

class IdIndex {
public:
    std::unordered_map<std::string, std::vector<int>> index;

    void add(const std::string& id, int idx) {
        index[id].push_back(idx);
    }

    std::optional<int> lookup(const std::string& id) const {
        auto it = index.find(id);
        if (it != index.end() && !it->second.empty()) {
            return it->second.front();
        }
        return std::nullopt;
    }

    void clear() { index.clear(); }
};

}  // namespace gffsub

#endif  // GFFSUB_RECORD_HPP
