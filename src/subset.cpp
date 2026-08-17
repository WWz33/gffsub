#include "subset.hpp"

#include "string_utils.hpp"

#include <string>
#include <unordered_set>

namespace gffsub {

void subset(GffData& data, const SubsetParams& params) {
    if (params.region) {
        filter_by_region(data, *params.region);
    }

    if (!params.bed_file.empty()) {
        filter_by_regions_from_file(data, params.bed_file);
    }

    if (!params.seqid_filter.empty()) {
        bool seqid_exclude = false;
        auto seqids = parse_list(params.seqid_filter, seqid_exclude);
        filter_by_seqid(data, seqids, seqid_exclude);
    }

    if (!params.source_filter.empty()) {
        bool source_exclude = false;
        auto sources = parse_list(params.source_filter, source_exclude);
        filter_by_source(data, sources, source_exclude);
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
        filter_longest_isoform(data, params.longest_type, params.threads);
    }
    if (!params.type.empty()) {
        bool type_exclude = false;
        auto types = parse_list(params.type, type_exclude);
        filter_by_type(data, types, type_exclude);
    }
}

}  // namespace gffsub
