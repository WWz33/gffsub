#include "annotation.hpp"
#include "gtf_parser.hpp"
#include "parser.hpp"
#include "region.hpp"

#include <deque>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace gffsub {

AnnotationIndex AnnotationIndex::from_file(const std::string& path) {
    GffData data;
    if (parse_file(path, data, infer_input_format(path)) != 0) {
        throw std::runtime_error("cannot parse file: " + path);
    }
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex AnnotationIndex::from_gff3(const std::string& path) {
    GffData data;
    if (parse_file(path, data, InputFormat::GFF3) != 0) {
        throw std::runtime_error("cannot parse GFF3 file: " + path);
    }
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex AnnotationIndex::from_data(GffData data) {
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex::AnnotationIndex(GffData data) : data_(std::move(data)) {
    for (int i = 0; i < static_cast<int>(data_.records.size()); ++i) {
        const auto& rec = data_.records[i];
        if (rec.id) {
            id_to_record_.emplace(*rec.id, i);
            // GFF3 discontinuous features: the same ID may appear on multiple
            // lines; keep every line index.
            id_to_records_[*rec.id].push_back(i);
        }

        if (rec.type != "gene") {
            continue;
        }

        auto add_gene_key = [&](const std::string& key) {
            if (key.empty()) return;
            auto& vec = gene_lookup_[key];
            if (std::find(vec.begin(), vec.end(), i) == vec.end()) {
                vec.push_back(i);
            }
        };

        if (rec.id) {
            add_gene_key(*rec.id);
        }
        if (rec.gene_id) {
            add_gene_key(*rec.gene_id);
        }

        const auto attrs = parse_attributes(rec.attr_raw);
        for (const char* key : {"Name", "gene_id", "locus_tag", "Alias", "Dbxref"}) {
            const auto it = attrs.find(key);
            if (it == attrs.end()) {
                continue;
            }
            for (const auto& value : it->second) {
                add_gene_key(value);
            }
        }
    }

    for (int i = 0; i < static_cast<int>(data_.records.size()); ++i) {
        const auto& rec = data_.records[i];
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto parent_it = attrs.find("Parent");

        // Collect parent IDs: from the Parent= attribute (GFF3) and from the
        // synthesized rec.parent_id field (GTF, set by the parser from
        // gene_id/transcript_id).
        std::vector<std::string> parent_ids;
        if (parent_it != attrs.end()) {
            parent_ids = parent_it->second;
        }
        if (rec.parent_id) {
            if (std::find(parent_ids.begin(), parent_ids.end(), *rec.parent_id) == parent_ids.end()) {
                parent_ids.push_back(*rec.parent_id);
            }
        }

        for (const auto& parent_id : parent_ids) {
            auto& children = children_by_parent_id_[parent_id];
            if (std::find(children.begin(), children.end(), i) == children.end()) {
                children.push_back(i);
            }
            const auto parent_record = id_to_record_.find(parent_id);
            if (rec.id && parent_record != id_to_record_.end()) {
                auto& parents = parents_by_child_id_[*rec.id];
                if (std::find(parents.begin(), parents.end(), parent_record->second) == parents.end()) {
                    parents.push_back(parent_record->second);
                }
            }
        }
    }
}

std::optional<GffRecord> AnnotationIndex::find_by_id(std::string_view id) const {
    const auto it = id_to_record_.find(std::string{id});
    if (it == id_to_record_.end()) {
        return std::nullopt;
    }
    return data_.records[it->second];
}

std::vector<GffRecord> AnnotationIndex::find_all_by_id(std::string_view id) const {
    std::vector<GffRecord> records;
    const auto it = id_to_records_.find(std::string{id});
    if (it == id_to_records_.end()) {
        return records;
    }
    records.reserve(it->second.size());
    for (const int idx : it->second) {
        records.push_back(data_.records[idx]);
    }
    return records;
}

std::optional<GffRecord> AnnotationIndex::find_gene(std::string_view id) const {
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end() || it->second.empty()) {
        return std::nullopt;
    }
    return data_.records[it->second.front()];
}

std::vector<GffRecord> AnnotationIndex::find_all_genes(std::string_view id) const {
    std::vector<GffRecord> result;
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end()) {
        return result;
    }
    for (int idx : it->second) {
        result.push_back(data_.records[idx]);
    }
    return result;
}

std::vector<GffRecord> AnnotationIndex::parents_of(std::string_view id) const {
    std::vector<GffRecord> parents;
    const auto it = parents_by_child_id_.find(std::string{id});
    if (it == parents_by_child_id_.end()) {
        return parents;
    }
    parents.reserve(it->second.size());
    for (const int idx : it->second) {
        parents.push_back(data_.records[idx]);
    }
    return parents;
}

std::vector<GffRecord> AnnotationIndex::children_of(std::string_view parent_id) const {
    std::vector<GffRecord> children;
    const auto it = children_by_parent_id_.find(std::string{parent_id});
    if (it == children_by_parent_id_.end()) {
        return children;
    }
    children.reserve(it->second.size());
    for (const int idx : it->second) {
        children.push_back(data_.records[idx]);
    }
    return children;
}

std::vector<GffRecord> AnnotationIndex::descendants_of(std::string_view parent_id) const {
    std::vector<GffRecord> descendants;
    std::deque<std::string> pending{std::string{parent_id}};
    std::unordered_set<std::string> visited_ids;
    std::unordered_set<int> visited_indices;

    while (!pending.empty()) {
        const auto current = pending.front();
        pending.pop_front();
        if (!visited_ids.insert(current).second) {
            continue;
        }

        const auto child_it = children_by_parent_id_.find(current);
        if (child_it == children_by_parent_id_.end()) {
            continue;
        }

        for (const int child_idx : child_it->second) {
            if (!visited_indices.insert(child_idx).second) {
                continue;
            }
            const auto& child = data_.records[child_idx];
            descendants.push_back(child);
            if (child.id) {
                pending.push_back(*child.id);
            }
        }
    }

    return descendants;
}

std::vector<GffRecord> AnnotationIndex::overlap(std::string_view seqid, int64_t start, int64_t end) const {
    std::vector<GffRecord> matches;
    for (const auto& rec : data_.records) {
        if (rec.seqid == seqid && rec.end >= start && rec.start <= end) {
            matches.push_back(rec);
        }
    }
    return matches;
}

std::optional<GffRecord> AnnotationIndex::nearest_gene(std::string_view seqid, int64_t start, int64_t end) const {
    std::optional<GffRecord> nearest;
    int64_t nearest_distance = std::numeric_limits<int64_t>::max();

    for (const auto& rec : data_.records) {
        if (rec.type != "gene" || rec.seqid != seqid) {
            continue;
        }

        int64_t distance = 0;
        if (end < rec.start) {
            distance = rec.start - end;
        } else if (start > rec.end) {
            distance = start - rec.end;
        }

        if (distance < nearest_distance) {
            nearest_distance = distance;
            nearest = rec;
        }
    }

    return nearest;
}

std::vector<GffRecord> AnnotationIndex::with_attribute(std::string_view key, std::string_view value) const {
    std::vector<GffRecord> matches;
    std::string skey{key};
    // Accept the attr.KEY form used by --grep/-I for arbitrary column-9 keys.
    if (skey.rfind("attr.", 0) == 0) {
        skey = skey.substr(5);
    }
    const std::string svalue{value};
    for (const auto& rec : data_.records) {
        // Check synthesized rec fields first (needed for GTF where attr_raw
        // uses key "value"; format that parse_attributes cannot parse).
        if (skey == "gene_id" && rec.gene_id && *rec.gene_id == svalue) {
            matches.push_back(rec);
            continue;
        }
        if (skey == "transcript_id" && rec.transcript_id && *rec.transcript_id == svalue) {
            matches.push_back(rec);
            continue;
        }
        if (skey == "ID" && rec.id && *rec.id == svalue) {
            matches.push_back(rec);
            continue;
        }
        if (skey == "Parent" && rec.parent_id && *rec.parent_id == svalue) {
            matches.push_back(rec);
            continue;
        }

        const auto attrs = parse_attributes(rec.attr_raw);
        const auto it = attrs.find(skey);
        if (it != attrs.end()) {
            for (const auto& attr_value : it->second) {
                if (attr_value == svalue) {
                    matches.push_back(rec);
                    break;
                }
            }
            continue;
        }
        // GTF fallback: `key "value";` attributes (same as record_field_value).
        if (extract_quoted_value(rec.attr_raw, skey) == std::optional<std::string>{svalue}) {
            matches.push_back(rec);
        }
    }
    return matches;
}

std::optional<GeneModel> AnnotationIndex::gene_model(std::string_view id) const {
    std::optional<GffRecord> gene = find_gene(id);
    if (!gene) {
        auto rec = find_by_id(id);
        std::unordered_set<std::string> visited;
        while (rec && rec->type != "gene" && rec->id && visited.insert(*rec->id).second) {
            const auto parents = parents_of(*rec->id);
            if (parents.empty()) {
                rec = std::nullopt;
            } else {
                rec = parents.front();
            }
        }
        if (rec && rec->type == "gene") {
            gene = rec;
        }
    }

    if (!gene || !gene->id) {
        return std::nullopt;
    }

    GeneModel model{*gene, {}};
    model.records.push_back(*gene);

    std::unordered_set<int> seen{gene->line_idx};
    for (const auto& rec : descendants_of(*gene->id)) {
        if (seen.insert(rec.line_idx).second) {
            model.records.push_back(rec);
        }
    }

    std::sort(model.records.begin(), model.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });

    return model;
}

}  // namespace gffsub
