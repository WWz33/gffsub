#ifndef GFFSUB_FILTER_HPP
#define GFFSUB_FILTER_HPP

#include "record.hpp"
#include "region.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_set>

namespace gffsub {

void filter_by_region(GffData& data, const Region& region);
void filter_by_region_exclude(GffData& data, const Region& region);
// exclude=true keeps records NOT overlapping any region. strand_mode is 0
// (ignore), 's' (same strand as the BED row), or 'o' (opposite strand).
void filter_by_regions_from_file(GffData& data, const std::string& bed_file,
                                 bool exclude = false, char strand_mode = 0);
void filter_by_type(GffData& data, const std::unordered_set<std::string>& types, bool exclude);
void filter_by_seqid(GffData& data, const std::unordered_set<std::string>& seqids, bool exclude);
void filter_by_source(GffData& data, const std::unordered_set<std::string>& sources, bool exclude);
void filter_by_score(GffData& data, std::optional<double> score);
void filter_by_strand(GffData& data, char strand);
void filter_by_phase(GffData& data, char phase);
void filter_longest_isoform(GffData& data, std::string_view isoform_type, size_t num_threads = 1);

}  // namespace gffsub

#endif  // GFFSUB_FILTER_HPP
