#include "subset.hpp"

#include <string>
#include <string_view>
#include <unordered_set>

namespace gffsub {

void subset(GffData& data, IdIndex& idx, const SubsetParams& params) {
    if (params.region) {
        filter_by_region(data, *params.region);
    }

    if (!params.bed_file.empty()) {
        filter_by_regions_from_file(data, params.bed_file);
    }

    if (!params.seqid_filter.empty()) {
        bool seqid_exclude = false;
        std::string_view sv = params.seqid_filter;
        if (!sv.empty() && sv.front() == '^') {
            seqid_exclude = true;
            sv.remove_prefix(1);
        }
        std::unordered_set<std::string> seqids;
        size_t pos = 0;
        while (pos < sv.size()) {
            auto comma = sv.find(',', pos);
            if (comma == std::string_view::npos) comma = sv.size();
            if (comma > pos) seqids.emplace(sv.substr(pos, comma - pos));
            pos = comma + 1;
        }
        filter_by_seqid(data, seqids, seqid_exclude);
    }

    if (!params.source_filter.empty()) {
        filter_by_source(data, params.source_filter);
    }

    if (params.score_filter) {
        filter_by_score(data, *params.score_filter);
    }

    if (params.strand_filter) {
        filter_by_strand(data, *params.strand_filter);
    }

    if (params.phase_filter) {
        filter_by_phase(data, *params.phase_filter);
    }

    if (!params.grep_filters.empty()) {
        filter_by_grep(data, params.grep_filters, params.invert_grep);
    }

    if (!params.include_exprs.empty()) {
        filter_by_expr(data, params.include_exprs, true);
    }

    if (!params.exclude_exprs.empty()) {
        filter_by_expr(data, params.exclude_exprs, false);
    }

    if (params.longest) {
        filter_longest_isoform(data, idx, params.feature, params.threads);
    } else if (!params.feature.empty()) {
        filter_by_feature(data, params.feature);
    }
}

}  // namespace gffsub
