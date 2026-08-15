#include "query.hpp"

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
    result.emit_summary = !params.summary_format.empty() || !params.output_attrs.empty();
    std::unordered_set<int> seen;

    auto add_summary = [&](const std::string& query_id, const std::string& matched_by, const GffRecord& rec) {
        if (result.emit_summary) {
            auto row = make_summary_row(index, query_id, matched_by, rec);
            if (!params.output_attrs.empty()) {
                row.attrs = extract_output_attrs(rec, params.output_attrs);
            }
            result.summary_rows.push_back(std::move(row));
        }
    };

    auto add_match = [&](const GffRecord& rec, const std::string& query_id, const std::string& matched_by) {
        const auto type_ok = [&](const GffRecord& r) {
            return !params.apply_type_filter || params.feature_type.empty() || r.type == params.feature_type;
        };
        if (params.include_model) {
            const auto model = index.gene_model(record_id(rec));
            if (model) {
                for (const auto& model_rec : model->records) {
                    if (type_ok(model_rec)) {
                        if (append_unique(result.records, seen, model_rec)) {
                            add_summary(query_id, "model", model_rec);
                        }
                    }
                }
                return;
            }
        }
        if (type_ok(rec)) {
            if (append_unique(result.records, seen, rec)) {
                add_summary(query_id, matched_by, rec);
            }
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
                    if (append_unique(result.records, seen, parent)) {
                        add_summary(query_id, "parent", parent);
                    }
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
                    if (append_unique(result.records, seen, child)) {
                        add_summary(query_id, "child", child);
                    }
                }
            }
        }
    };

    for (const auto& id : params.ids) {
        // GFF3 discontinuous features: one ID may span multiple lines.
        const auto recs = index.find_all_by_id(id);
        if (!recs.empty()) {
            for (const auto& rec : recs) {
                add_match(rec, id, "ID");
            }
        } else if (result.emit_summary) {
            auto row = make_not_found_row(id, "ID");
            row.attrs.assign(params.output_attrs.size(), "");
            result.summary_rows.push_back(std::move(row));
        }
    }

    if (!params.name.empty()) {
        const auto rec = index.find_gene(params.name);
        if (rec) {
            add_match(*rec, params.name, infer_gene_match_key(index, params.name, *rec));
        } else if (result.emit_summary) {
            auto row = make_not_found_row(params.name, "name");
            row.attrs.assign(params.output_attrs.size(), "");
            result.summary_rows.push_back(std::move(row));
        }
    }

    for (const auto& [key, value] : params.attr_filters) {
        bool matched = false;
        for (const auto& rec : index.with_attribute(key, value)) {
            matched = true;
            add_match(rec, key + "=" + value, key);
        }
        if (!matched && result.emit_summary) {
            auto row = make_not_found_row(key + "=" + value, key);
            row.attrs.assign(params.output_attrs.size(), "");
            result.summary_rows.push_back(std::move(row));
        }
    }

    if (params.region) {
        const std::string region_str = params.region->seqid + ":" +
            std::to_string(params.region->start) + "-" +
            std::to_string(params.region->end);
        for (const auto& rec : index.overlap(params.region->seqid,
                                             params.region->start,
                                             params.region->end)) {
            add_match(rec, region_str, "region");
        }
    }

    if (params.nearest_region) {
        const auto rec = index.nearest_gene(params.nearest_region->seqid,
                                            params.nearest_region->start,
                                            params.nearest_region->end);
        if (rec) {
            add_match(*rec,
                      params.nearest_region->seqid + ":" +
                          std::to_string(params.nearest_region->start) + "-" +
                          std::to_string(params.nearest_region->end),
                      "nearest");
        } else if (result.emit_summary) {
            auto row = make_not_found_row(
                params.nearest_region->seqid + ":" +
                    std::to_string(params.nearest_region->start) + "-" +
                    std::to_string(params.nearest_region->end),
                "nearest");
            row.attrs.assign(params.output_attrs.size(), "");
            result.summary_rows.push_back(std::move(row));
        }
    }

    std::sort(result.records.records.begin(), result.records.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });

    return result;
}

void print_query_result(std::ostream& out, const QueryResult& result, const QueryParams& params) {
    if (params.summary_format == "json") {
        print_summary_json(out, result.summary_rows, params.output_attrs);
    } else if (result.emit_summary) {
        print_summary_tsv(out, result.summary_rows, params.output_attrs);
    } else {
        print_gff3(out, result.records);
    }
}

}  // namespace gffsub
