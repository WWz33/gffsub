#include "region.hpp"

namespace gffsub {

std::optional<Region> parse_region(std::string_view region_str) {
    size_t colon = region_str.find(':');
    if (colon == std::string_view::npos) return std::nullopt;

    std::string seqid(region_str.substr(0, colon));
    auto range_part = region_str.substr(colon + 1);

    size_t dash = range_part.find('-');
    if (dash == std::string_view::npos) return std::nullopt;
    // Reject trailing garbage after the range (e.g. "chr1:1-100:300").
    if (range_part.find('-', dash + 1) != std::string_view::npos) return std::nullopt;

    int64_t start = 0;
    int64_t end = 0;
    try {
        size_t pos = 0;
        const std::string start_str{range_part.substr(0, dash)};
        start = std::stoll(start_str, &pos);
        if (pos != start_str.size()) return std::nullopt;
        const std::string end_str{range_part.substr(dash + 1)};
        end = std::stoll(end_str, &pos);
        if (pos != end_str.size()) return std::nullopt;
    } catch (const std::exception&) {
        return std::nullopt;
    }

    if (start < 1 || end < 1 || start > end) {
        return std::nullopt;
    }

    return Region{seqid, start, end};
}

BedRegion to_bed_region(const GffRecord& rec) {
    return BedRegion{std::string{rec.seqid}, rec.start - 1, rec.end};
}

Region from_bed_region(const BedRegion& region) {
    return Region{region.seqid, region.start + 1, region.end};
}

Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware) {
    int64_t left_extension = upstream;
    int64_t right_extension = downstream;
    if (strand_aware && rec.strand == '-') {
        left_extension = downstream;
        right_extension = upstream;
    }

    int64_t start = rec.start - left_extension;
    if (start < 1) {
        start = 1;
    }
    return Region{std::string{rec.seqid}, start, rec.end + right_extension};
}

}  // namespace gffsub
