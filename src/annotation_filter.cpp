#include "filter.hpp"
#include "region.hpp"
#include "record.hpp"
#include "string_utils.hpp"
#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

namespace gffsub {

void filter_by_region(GffData& data, const Region& region) {
    for (auto& rec : data) {
        if (!rec.kept) continue;
        if (rec.seqid != region.seqid || rec.end < region.start || rec.start > region.end) {
            rec.kept = false;
        }
    }
}

void filter_by_region_exclude(GffData& data, const Region& region) {
    for (auto& rec : data) {
        if (!rec.kept) continue;
        if (rec.seqid == region.seqid && rec.end >= region.start && rec.start <= region.end) {
            rec.kept = false;
        }
    }
}

static std::vector<Region> load_regions(const std::string& filename, bool is_bed) {
    std::vector<Region> regions;
    if (is_bed) {
        std::ifstream file(filename);
        // An unopenable file must not masquerade as "zero regions" (which
        // would silently drop every record); surface it as an error.
        if (!file.is_open()) {
            throw std::runtime_error("cannot open BED file: " + filename);
        }

        std::string line;
        while (std::getline(file, line)) {
            // Strip CR from CRLF line endings before validation.
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty() || line[0] == '#') continue;
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            Region r;
            try {
                size_t pos = 0;
                const int64_t bed_start = std::stoll(std::string{cols[1]}, &pos);
                if (pos != cols[1].size()) continue;
                pos = 0;
                r.end = std::stoll(std::string{cols[2]}, &pos);
                if (pos != cols[2].size()) continue;
                if (bed_start < 0 || r.end < bed_start) continue;
                r.seqid = std::string{cols[0]};
                // Zero-length BED [s,s): insertion point, mapped to the
                // one-base 1-based position s+1.
                if (r.end == bed_start) r.end = bed_start + 1;
                r.start = bed_start + 1;  // BED 0-based half-open -> 1-based inclusive
                if (cols.size() > 5 && (cols[5] == "+" || cols[5] == "-")) {
                    r.strand = cols[5][0];
                }
            } catch (const std::exception&) {
                continue;
            }
            regions.push_back(r);
        }
    }
    return regions;
}

// Per-(seqid, strand group) interval index: region starts sorted, with the
// prefix maximum of ends. An overlap with [qs,qe] exists iff some region
// with start <= qe also has end >= qs, answerable in O(log m) by binary
// search + the prefix max.
namespace {

struct IntervalIndex {
    std::vector<int64_t> starts;
    std::vector<int64_t> prefix_max_end;

    void build(std::vector<Region>& regions) {
        std::sort(regions.begin(), regions.end(),
                  [](const Region& a, const Region& b) { return a.start < b.start; });
        starts.reserve(regions.size());
        prefix_max_end.reserve(regions.size());
        int64_t max_end = std::numeric_limits<int64_t>::min();
        for (const auto& r : regions) {
            starts.push_back(r.start);
            max_end = std::max(max_end, r.end);
            prefix_max_end.push_back(max_end);
        }
    }

    bool overlaps(int64_t qs, int64_t qe) const {
        const auto hi = std::upper_bound(starts.begin(), starts.end(), qe);
        if (hi == starts.begin()) return false;
        return prefix_max_end[static_cast<size_t>(hi - starts.begin()) - 1] >= qs;
    }
};

}  // namespace

void filter_by_regions_from_file(GffData& data, const std::string& bed_file,
                                 bool exclude, char strand_mode) {
    auto regions = load_regions(bed_file, true);
    if (strand_mode != 0) {
        // Without a strand column a BED row cannot participate in a
        // same/opposite-strand comparison.
        regions.erase(std::remove_if(regions.begin(), regions.end(),
                                     [](const Region& r) {
                                         return r.strand != '+' && r.strand != '-';
                                     }),
                      regions.end());
    }

    struct SeqIndex {
        IntervalIndex all;    // strand_mode == 0
        IntervalIndex plus;   // strand_mode != 0
        IntervalIndex minus;
    };
    std::unordered_map<std::string, SeqIndex> index;
    {
        std::unordered_map<std::string, std::vector<Region>> all_by_seqid;
        std::unordered_map<std::string, std::vector<Region>> plus_by_seqid;
        std::unordered_map<std::string, std::vector<Region>> minus_by_seqid;
        for (auto& r : regions) {
            if (strand_mode == 0) {
                all_by_seqid[r.seqid].push_back(std::move(r));
            } else if (r.strand == '+') {
                plus_by_seqid[r.seqid].push_back(std::move(r));
            } else {
                minus_by_seqid[r.seqid].push_back(std::move(r));
            }
        }
        for (auto& [seqid, v] : all_by_seqid) {
            index[seqid].all.build(v);
        }
        for (auto& [seqid, v] : plus_by_seqid) {
            index[seqid].plus.build(v);
        }
        for (auto& [seqid, v] : minus_by_seqid) {
            index[seqid].minus.build(v);
        }
    }

    std::string seqid_buf;  // reused lookup key, avoids per-record allocation
    for (auto& rec : data) {
        if (!rec.kept) continue;
        seqid_buf.assign(rec.seqid.data(), rec.seqid.size());
        const auto it = index.find(seqid_buf);
        bool in_region = false;
        if (it != index.end()) {
            const auto& seq_index = it->second;
            if (strand_mode == 0) {
                in_region = seq_index.all.overlaps(rec.start, rec.end);
            } else {
                // '.'/'?' records have no strand to compare and never match.
                const IntervalIndex* group = nullptr;
                if (strand_mode == 's') {
                    if (rec.strand == '+') group = &seq_index.plus;
                    else if (rec.strand == '-') group = &seq_index.minus;
                } else {  // 'o'
                    if (rec.strand == '+') group = &seq_index.minus;
                    else if (rec.strand == '-') group = &seq_index.plus;
                }
                in_region = group && group->overlaps(rec.start, rec.end);
            }
        }
        if (exclude ? in_region : !in_region) {
            rec.kept = false;
        }
    }
}

void filter_by_type(GffData& data, const std::unordered_set<std::string>& types, bool exclude) {
    for (auto& rec : data) {
        if (!rec.kept) continue;
        bool found = types.count(std::string{rec.type}) > 0;
        if (exclude ? found : !found) {
            rec.kept = false;
        }
    }
}

void filter_by_seqid(GffData& data, const std::unordered_set<std::string>& seqids, bool exclude) {
    for (auto& rec : data) {
        if (!rec.kept) continue;
        bool found = seqids.count(std::string{rec.seqid}) > 0;
        if (exclude ? found : !found) {
            rec.kept = false;
        }
    }
}

void filter_by_source(GffData& data, const std::unordered_set<std::string>& sources, bool exclude) {
    for (auto& rec : data) {
        if (!rec.kept) continue;
        bool found = sources.count(std::string{rec.source}) > 0;
        if (exclude ? found : !found) {
            rec.kept = false;
        }
    }
}

void filter_by_score(GffData& data, std::optional<double> score) {
    for (auto& rec : data) {
        if (rec.kept && rec.score != score) {
            rec.kept = false;
        }
    }
}

void filter_by_strand(GffData& data, char strand) {
    for (auto& rec : data) {
        if (rec.kept && rec.strand != strand) {
            rec.kept = false;
        }
    }
}

void filter_by_phase(GffData& data, char phase) {
    for (auto& rec : data) {
        if (rec.kept && rec.phase != phase) {
            rec.kept = false;
        }
    }
}

}  // namespace gffsub
