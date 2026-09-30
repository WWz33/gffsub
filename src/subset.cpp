#include "subset.hpp"

#include "parser.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace gffsub {

namespace {

// Remove records all of whose Parent references point outside the kept set.
// Runs to a fixpoint because dropping a transcript orphans its exons in
// turn. Records without parents are never orphans.
void drop_orphans(GffData& data) {
    std::vector<std::vector<std::string>> parents(data.records.size());
    for (size_t i = 0; i < data.records.size(); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto it = attrs.find("Parent");
        if (it != attrs.end()) {
            parents[i] = it->second;
        }
        if (rec.parent_id &&
            std::find(parents[i].begin(), parents[i].end(), *rec.parent_id) == parents[i].end()) {
            parents[i].push_back(*rec.parent_id);
        }
    }

    bool changed = true;
    while (changed) {
        changed = false;
        std::unordered_set<std::string> kept_ids;
        for (const auto& rec : data.records) {
            if (rec.kept && rec.id) kept_ids.insert(std::string{*rec.id});
        }
        for (size_t i = 0; i < data.records.size(); ++i) {
            if (!data.records[i].kept || parents[i].empty()) continue;
            const bool has_kept_parent =
                std::any_of(parents[i].begin(), parents[i].end(),
                            [&](const std::string& p) { return kept_ids.count(p) > 0; });
            if (!has_kept_parent) {
                data.records[i].kept = false;
                changed = true;
            }
        }
    }
}

}  // namespace

void subset(GffData& data, const SubsetParams& params) {
    if (params.region) {
        filter_by_region(data, *params.region);
    }

    if (params.exclude_region) {
        filter_by_region_exclude(data, *params.exclude_region);
    }

    if (!params.bed_file.empty()) {
        filter_by_regions_from_file(data, params.bed_file, params.bed_exclude,
                                    params.bed_strand_mode);
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
    if (params.drop_orphans) {
        drop_orphans(data);
    }
}

int sort_records(GffData& data, std::string_view keys, bool reverse) {
    enum class Key { Seqid, NaturalSeqid, Start, End, Length, Type };
    std::vector<Key> order;
    for (const auto k : split_line(keys, ',')) {
        if (k == "seqid") order.push_back(Key::Seqid);
        else if (k == "natural-seqid") order.push_back(Key::NaturalSeqid);
        else if (k == "start") order.push_back(Key::Start);
        else if (k == "end") order.push_back(Key::End);
        else if (k == "length") order.push_back(Key::Length);
        else if (k == "type") order.push_back(Key::Type);
        else return -1;
    }
    if (order.empty()) return -1;

    const auto cmp = [&order](const GffRecord& a, const GffRecord& b) {
        for (const Key k : order) {
            int c = 0;
            switch (k) {
                case Key::Seqid: c = a.seqid.compare(b.seqid); break;
                case Key::NaturalSeqid: c = natural_compare(a.seqid, b.seqid); break;
                case Key::Start: c = a.start < b.start ? -1 : a.start > b.start ? 1 : 0; break;
                case Key::End: c = a.end < b.end ? -1 : a.end > b.end ? 1 : 0; break;
                case Key::Length: {
                    const int64_t la = a.end - a.start, lb = b.end - b.start;
                    c = la < lb ? -1 : la > lb ? 1 : 0;
                    break;
                }
                case Key::Type: {
                    // FeatureClass enum order is the hierarchy order:
                    // Gene < Transcript < Exon < CDS < UTR < Match < Region.
                    const int ra = static_cast<int>(a.feat_class),
                              rb = static_cast<int>(b.feat_class);
                    c = ra < rb ? -1 : ra > rb ? 1 : a.type.compare(b.type);
                    break;
                }
            }
            if (c != 0) return c;
        }
        return 0;
    };

    std::stable_sort(data.records.begin(), data.records.end(),
                     [&](const GffRecord& a, const GffRecord& b) {
                         return reverse ? cmp(b, a) < 0 : cmp(a, b) < 0;
                     });
    return 0;
}

}  // namespace gffsub
