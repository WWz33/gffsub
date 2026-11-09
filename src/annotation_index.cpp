#include "annotation.hpp"
#include "parser.hpp"
#include "region.hpp"

#include <deque>
#include <algorithm>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace gffsub {

namespace {

// The spec requires IDs to be unique; a repeated ID on non-matching lines is
// invalid input. This tool does not rewrite IDs, so report it instead.
void warn_duplicate_feature_id(std::string_view id) {
    constexpr int kLimit = 5;
    static int warned = 0;
    if (warned > kLimit) return;
    if (warned == kLimit) {
        std::cerr << "Warning: further duplicate-ID warnings suppressed\n";
    } else {
        std::cerr << "Warning: ID '" << id
                  << "' is shared by lines that are not one feature; the GFF3"
                     " spec requires unique IDs, so lineage lookups may pick"
                     " the wrong line\n";
    }
    ++warned;
}

}  // namespace

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

AnnotationIndex AnnotationIndex::from_owned(GffData data) {
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex AnnotationIndex::from_data(const GffData& data) {
    // Borrowed view: no record copy, maps index into the caller's vector.
    // Valid while the caller's GffData lives and keeps its record order.
    AnnotationIndex index{GffData{}};
    index.recs_ = &data.records;
    index.source_format_ = data.source_format;
    index.directives_ = data.directives;
    index.build_maps(data.records);
    return index;
}

AnnotationIndex::AnnotationIndex(GffData data)
    : data_(std::move(data)), directives_(data_.directives),
      source_format_(data_.source_format) {
    build_maps(data_.records);
}

AnnotationIndex::AnnotationIndex(const AnnotationIndex& other)
    : id_to_record_(other.id_to_record_),
      id_to_records_(other.id_to_records_),
      gene_lookup_(other.gene_lookup_),
      parents_by_child_id_(other.parents_by_child_id_),
      children_by_parent_id_(other.children_by_parent_id_),
      directives_(other.directives_),
      source_format_(other.source_format_) {
    data_.copy_storage_from(other.data_);
    if (other.recs_ == &other.data_.records) {
        recs_ = &data_.records;
    } else {
        recs_ = other.recs_;
    }
}

AnnotationIndex::AnnotationIndex(AnnotationIndex&& other) noexcept
    : data_(std::move(other.data_)),
      id_to_record_(std::move(other.id_to_record_)),
      id_to_records_(std::move(other.id_to_records_)),
      gene_lookup_(std::move(other.gene_lookup_)),
      parents_by_child_id_(std::move(other.parents_by_child_id_)),
      children_by_parent_id_(std::move(other.children_by_parent_id_)),
      directives_(std::move(other.directives_)),
      source_format_(other.source_format_) {
    if (other.recs_ == &other.data_.records) {
        recs_ = &data_.records;
    } else {
        recs_ = other.recs_;
    }
}

AnnotationIndex& AnnotationIndex::operator=(const AnnotationIndex& other) {
    if (this == &other) return *this;
    data_.copy_storage_from(other.data_);
    id_to_record_ = other.id_to_record_;
    id_to_records_ = other.id_to_records_;
    gene_lookup_ = other.gene_lookup_;
    parents_by_child_id_ = other.parents_by_child_id_;
    children_by_parent_id_ = other.children_by_parent_id_;
    directives_ = other.directives_;
    source_format_ = other.source_format_;
    if (other.recs_ == &other.data_.records) {
        recs_ = &data_.records;
    } else {
        recs_ = other.recs_;
    }
    return *this;
}

AnnotationIndex& AnnotationIndex::operator=(AnnotationIndex&& other) noexcept {
    if (this == &other) return *this;
    data_ = std::move(other.data_);
    id_to_record_ = std::move(other.id_to_record_);
    id_to_records_ = std::move(other.id_to_records_);
    gene_lookup_ = std::move(other.gene_lookup_);
    parents_by_child_id_ = std::move(other.parents_by_child_id_);
    children_by_parent_id_ = std::move(other.children_by_parent_id_);
    directives_ = std::move(other.directives_);
    source_format_ = other.source_format_;
    if (other.recs_ == &other.data_.records) {
        recs_ = &data_.records;
    } else {
        recs_ = other.recs_;
    }
    return *this;
}

std::unordered_map<std::string, std::vector<int>> build_children_index(
    const std::vector<GffRecord>& records, bool kept_only) {
    std::unordered_map<std::string, std::vector<int>> children;
    for (int i = 0; i < static_cast<int>(records.size()); ++i) {
        const auto& rec = records[i];
        if (kept_only && !rec.kept) continue;
        // parent_ids() returns each parent once per record, so a plain push
        // cannot duplicate an index under the same parent.
        for (const auto& parent_id : parent_ids(rec)) {
            children[parent_id].push_back(i);
        }
    }
    return children;
}

void AnnotationIndex::build_maps(const std::vector<GffRecord>& records) {
    // Flat GTF fallback: GTF2.2 has no gene feature type, grouping is done
    // purely through the gene_id attribute. When no gene-classified record
    // exists, map each gene_id to every record carrying it so --name returns
    // the gene's lines (exon/CDS records carry no ID of their own, so there
    // is nothing to anchor hierarchy expansion on). window/gene_model on a
    // flat file anchor on the first line's span, not a synthesized gene span.
    bool saw_gene = false;
    std::unordered_set<std::string> parent_record_missing;

    for (int i = 0; i < static_cast<int>(records.size()); ++i) {
        const auto& rec = records[i];
        if (rec.id) {
            id_to_record_.emplace(*rec.id, i);
            // GFF3 discontinuous features: the same ID may appear on multiple
            // lines; keep every line index. Lines that share an ID but differ
            // in type/seqid/strand are not one feature, which the spec forbids
            // (IDs must be unique); warned once here, looked up by first line.
            auto& same_id_lines = id_to_records_[*rec.id];
            if (!same_id_lines.empty()) {
                const auto& first = records[same_id_lines.front()];
                if (first.type != rec.type || first.seqid != rec.seqid ||
                    first.strand != rec.strand) {
                    warn_duplicate_feature_id(*rec.id);
                }
            }
            same_id_lines.push_back(i);
        }

        if (rec.feat_class == FeatureClass::Gene) {
            saw_gene = true;

            // One gene may carry the same value under several keys (ID and
            // gene_id, Name equal to ID); add each key once.
            std::unordered_set<std::string> added;
            auto add_gene_key = [&](const std::string& key) {
                if (key.empty()) return;
                if (added.insert(key).second) {
                    gene_lookup_[key].push_back(i);
                }
            };

            if (rec.id) {
                add_gene_key(*rec.id);
            }
            if (rec.gene_id) {
                add_gene_key(*rec.gene_id);
            }

            const auto attrs = decode_attributes(rec.attr_raw, rec.src_fmt);
            for (const char* key : {"Name", "gene_id", "locus_tag", "Alias", "Dbxref"}) {
                for (const auto& [k, value] : attrs) {
                    if (k == key) {
                        add_gene_key(value);
                    }
                }
            }
        }
    }

    if (!saw_gene) {
        // Flat file: group every record by its gene_id attribute. Only files
        // without gene records pay this extra pass. Each record index is
        // visited once, so no duplicate check is needed.
        for (int i = 0; i < static_cast<int>(records.size()); ++i) {
            const auto& rec = records[i];
            if (rec.gene_id && !rec.gene_id->empty()) {
                gene_lookup_[*rec.gene_id].push_back(i);
            }
        }
    }

    // Hierarchy maps. One builder produces parent -> children; the reverse
    // map and the missing-Parent diagnostics derive from it, so column 9 is
    // scanned once.
    children_by_parent_id_ = build_children_index(records, /*kept_only=*/false);
    for (const auto& [parent_id, children] : children_by_parent_id_) {
        const auto parent_record = id_to_record_.find(parent_id);
        if (parent_record != id_to_record_.end()) {
            for (const int child_idx : children) {
                const auto& child = records[child_idx];
                if (child.id) {
                    parents_by_child_id_[*child.id].push_back(parent_record->second);
                }
            }
            continue;
        }
        // A Parent that names no record in the file: the child links to
        // nothing. Report it once per ID rather than silently emitting a
        // hierarchy with a hole in it. GTF is exempt: a flat GTF has
        // exon/CDS rows whose transcript_id names no transcript row, which
        // is legal there.
        if (!children.empty() &&
            records[children.front()].src_fmt != InputFormat::GTF &&
            !parent_record_missing.count(parent_id)) {
            parent_record_missing.insert(parent_id);
            if (parent_record_missing.size() <= 5) {
                std::cerr << "Warning: Parent '" << parent_id
                          << "' has no record in the input\n";
            } else if (parent_record_missing.size() == 6) {
                std::cerr << "Warning: further missing-Parent warnings suppressed\n";
            }
        }
    }
    for (auto& [child_id, parent_lines] : parents_by_child_id_) {
        std::sort(parent_lines.begin(), parent_lines.end());
        parent_lines.erase(std::unique(parent_lines.begin(), parent_lines.end()),
                           parent_lines.end());
    }

    // GFF3 spec: Parent links that form a cycle are invalid. Every traversal
    // here is visited-guarded, so a cycle cannot hang the tool, but the
    // hierarchy it reports would be nonsense. Walk the parent graph with
    // three-colour marking (unvisited / on stack / done) over every Parent of
    // every line, so a cycle reachable only through a non-first Parent is
    // found too. No depth cutoff: a long chain is still a chain.
    {
        std::unordered_map<std::string, std::vector<std::string>> parents_of_id;
        for (const auto& [parent_id, children] : children_by_parent_id_) {
            for (const int child_idx : children) {
                const auto& child = records[child_idx];
                if (!child.id) continue;
                auto& list = parents_of_id[*child.id];
                if (std::find(list.begin(), list.end(), parent_id) == list.end()) {
                    list.push_back(parent_id);
                }
            }
        }
        // 0 = unvisited, 1 = on the current walk, 2 = fully explored.
        std::unordered_map<std::string, int> state;
        std::unordered_set<std::string> reported_cycles;
        int reported = 0;
        // Explicit stack: a flat file can chain thousands of levels, and the
        // recursion would run out of stack on real input.
        struct Frame {
            std::string id;
            size_t next_parent = 0;
        };
        for (const auto& entry : parents_of_id) {
            if (state.count(entry.first) != 0) continue;
            const std::string root = entry.first;
            std::vector<Frame> stack{{root, 0}};
            state[root] = 1;
            while (!stack.empty()) {
                Frame& frame = stack.back();
                const auto it = parents_of_id.find(frame.id);
                const size_t n = (it == parents_of_id.end()) ? 0 : it->second.size();
                if (frame.next_parent < n) {
                    const std::string& parent = it->second[frame.next_parent++];
                    const auto seen = state.find(parent);
                    if (seen != state.end() && seen->second == 1) {
                        if (reported_cycles.insert(root).second) {
                            if (reported < 5) {
                                std::cerr << "Warning: Parent cycle involving ID '"
                                          << frame.id << "'\n";
                            }
                            ++reported;
                        }
                        continue;
                    }
                    if (seen == state.end()) {
                        state[parent] = 1;
                        stack.push_back(Frame{parent, 0});
                    }
                    continue;
                }
                state[frame.id] = 2;
                stack.pop_back();
            }
        }
        if (reported > 5) {
            std::cerr << "Warning: " << (reported - 5) << " further Parent cycle warning(s) suppressed\n";
        }
    }
}

std::optional<GffRecord> AnnotationIndex::find_by_id(std::string_view id) const {
    const auto it = id_to_record_.find(std::string{id});
    if (it == id_to_record_.end()) {
        return std::nullopt;
    }
    return (*recs_)[it->second];
}

std::vector<GffRecord> AnnotationIndex::find_all_by_id(std::string_view id) const {
    std::vector<GffRecord> records;
    const auto it = id_to_records_.find(std::string{id});
    if (it == id_to_records_.end()) {
        return records;
    }
    records.reserve(it->second.size());
    for (const int idx : it->second) {
        records.push_back((*recs_)[idx]);
    }
    return records;
}

std::optional<GffRecord> AnnotationIndex::find_gene(std::string_view id) const {
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end() || it->second.empty()) {
        return std::nullopt;
    }
    return (*recs_)[it->second.front()];
}

std::vector<GffRecord> AnnotationIndex::find_all_genes(std::string_view id) const {
    std::vector<GffRecord> result;
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end()) {
        return result;
    }
    for (int idx : it->second) {
        result.push_back((*recs_)[idx]);
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
        parents.push_back((*recs_)[idx]);
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
        children.push_back((*recs_)[idx]);
    }
    return children;
}

std::vector<GffRecord> AnnotationIndex::parents_of_record(const GffRecord& rec) const {
    std::vector<GffRecord> parents;
    for (const auto& pid : parent_ids(rec)) {
        const auto it = id_to_record_.find(pid);
        if (it != id_to_record_.end()) {
            parents.push_back((*recs_)[it->second]);
        }
    }
    return parents;
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
            const auto& child = (*recs_)[child_idx];
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
    for (const auto& rec : *recs_) {
        if (rec.seqid == seqid && rec.end >= start && rec.start <= end) {
            matches.push_back(rec);
        }
    }
    return matches;
}

std::optional<GffRecord> AnnotationIndex::nearest_gene(std::string_view seqid, int64_t start, int64_t end) const {
    std::optional<GffRecord> nearest;
    int64_t nearest_distance = std::numeric_limits<int64_t>::max();

    for (const auto& rec : *recs_) {
        if (rec.feat_class != FeatureClass::Gene || rec.seqid != seqid) {
            continue;
        }

        int64_t distance = 0;
        if (end < rec.start) {
            distance = rec.start - end;
        } else if (start > rec.end) {
            distance = start - rec.end;
        }

        // Ties broken deterministically: smaller start, then smaller ID,
        // independent of record order in the file.
        bool better = false;
        if (distance < nearest_distance) {
            better = true;
        } else if (distance == nearest_distance && nearest) {
            if (rec.start != nearest->start) {
                better = rec.start < nearest->start;
            } else {
                const std::string rec_id = rec.id ? *rec.id : "";
                const std::string near_id = nearest->id ? *nearest->id : "";
                better = rec_id < near_id;
            }
        }
        if (better) {
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
    for (const auto& rec : *recs_) {
        // One format-aware value list per key; Parent carries the synthesized
        // fallback and drops a self-link.
        const auto values = (skey == "Parent") ? parent_ids(rec) : attr_values(rec, skey);
        for (const auto& attr_value : values) {
            if (attr_value == svalue) {
                matches.push_back(rec);
                break;
            }
        }
    }
    return matches;
}

std::optional<GeneModel> AnnotationIndex::gene_model(std::string_view id) const {
    std::optional<GffRecord> gene = find_gene(id);
    if (!gene) {
        auto rec = find_by_id(id);
        std::unordered_set<std::string> visited;
        while (rec && rec->feat_class != FeatureClass::Gene && rec->id && visited.insert(*rec->id).second) {
            const auto parents = parents_of(*rec->id);
            if (parents.empty()) {
                rec = std::nullopt;
            } else {
                rec = parents.front();
            }
        }
        if (rec && rec->feat_class == FeatureClass::Gene) {
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
