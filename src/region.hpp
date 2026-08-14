#ifndef GFFSUB_REGION_HPP
#define GFFSUB_REGION_HPP

#include "record.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace gffsub {

struct Region {
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
};

struct BedRegion {
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
};

std::optional<Region> parse_region(std::string_view region_str);
BedRegion to_bed_region(const GffRecord& rec);
Region from_bed_region(const BedRegion& region);
Region window_region(const GffRecord& rec, int64_t upstream, int64_t downstream, bool strand_aware);

}  // namespace gffsub

#endif  // GFFSUB_REGION_HPP
