#include "query.hpp"
#include "query_summary.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <string>
#include <unordered_set>

namespace gffsub {

namespace {

bool append_unique(GffData& out, std::unordered_set<int>& seen, const GffRecord& rec) {
    if (!seen.insert(rec.line_idx).second) {
        return false;
    }
    out.append(rec);
    return true;
}

}  // namespace

QueryResult query(const AnnotationIndex& index, const QueryParams& params) {
    QueryResult result;
    std::unordered_set<int> seen;

    // Parse -t value once: comma list with optional ^ exclusion.
    bool type_exclude = false;
    std::unordered_set<std::string> type_set;
    if (!params.type.empty()) {
        type_set = parse_list(params.type, type_exclude);
    }
    const bool use_type_filter = params.apply_type_filter && !params.type.empty();

    auto add_match = [&](const GffRecord& rec) {
        const auto type_ok = [&](const GffRecord& r) {
            if (!use_type_filter) return true;
            const bool found = type_set.count(std::string{r.type}) > 0;
            return type_exclude ? !found : found;
        };
        if (params.include_model) {
            const auto model = index.gene_model(record_id(rec));
            if (model) {
                for (const auto& model_rec : model->records) {
                    if (type_ok(model_rec)) {
                        append_unique(result.records, seen, model_rec);
                    }
                }
                return;
            }
        }
        if (type_ok(rec)) {
            append_unique(result.records, seen, rec);
        }
        if (params.include_parents && rec.id) {
            std::vector<GffRecord> stack = index.parents_of(*rec.id);
            std::unordered_set<int> visited_parents;
            while (!stack.empty()) {
                const auto parent = stack.back();
                stack.pop_back();
                if (!visited_parents.insert(parent.line_idx).second) {
                    continue;
                }
                if (type_ok(parent)) {
                    append_unique(result.records, seen, parent);
                }
                if (parent.id) {
                    for (const auto& grandparent : index.parents_of(*parent.id)) {
                        stack.push_back(grandparent);
                    }
                }
            }
        }
        if (params.include_children && rec.id) {
            for (const auto& child : index.descendants_of(*rec.id)) {
                if (type_ok(child)) {
                    append_unique(result.records, seen, child);
                }
            }
        }
    };

    for (const auto& id : params.ids) {
        const auto recs = index.find_all_by_id(id);
        for (const auto& rec : recs) {
            add_match(rec);
        }
    }

    if (!params.name.empty()) {
        for (const auto& rec : index.find_all_genes(params.name)) {
            add_match(rec);
        }
    }

    for (const auto& [key, value] : params.attr_filters) {
        for (const auto& rec : index.with_attribute(key, value)) {
            add_match(rec);
        }
    }

    if (params.region) {
        for (const auto& rec : index.overlap(params.region->seqid,
                                             params.region->start,
                                             params.region->end)) {
            add_match(rec);
        }
    }

    if (params.nearest_region) {
        const auto rec = index.nearest_gene(params.nearest_region->seqid,
                                            params.nearest_region->start,
                                            params.nearest_region->end);
        if (rec) {
            add_match(*rec);
        }
    }

    std::sort(result.records.records.begin(), result.records.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });

    return result;
}

}  // namespace gffsub
