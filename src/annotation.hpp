#ifndef GFFSUB_ANNOTATION_HPP
#define GFFSUB_ANNOTATION_HPP

#include "record.hpp"

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace gffsub {

struct GeneModel {
    GffRecord gene;
    std::vector<GffRecord> records;
};

class AnnotationIndex {
public:
    static AnnotationIndex from_file(const std::string& path);
    static AnnotationIndex from_gff3(const std::string& path);
    static AnnotationIndex from_data(GffData data);

    std::optional<GffRecord> find_by_id(std::string_view id) const;
    // All lines sharing an ID (GFF3 discontinuous features, e.g. multi-line CDS).
    std::vector<GffRecord> find_all_by_id(std::string_view id) const;
    std::optional<GffRecord> find_gene(std::string_view id) const;
    std::vector<GffRecord> find_all_genes(std::string_view id) const;
    std::vector<GffRecord> parents_of(std::string_view id) const;
    std::vector<GffRecord> children_of(std::string_view parent_id) const;
    std::vector<GffRecord> descendants_of(std::string_view parent_id) const;
    std::vector<GffRecord> overlap(std::string_view seqid, int64_t start, int64_t end) const;
    std::optional<GffRecord> nearest_gene(std::string_view seqid, int64_t start, int64_t end) const;
    std::vector<GffRecord> with_attribute(std::string_view key, std::string_view value) const;
    std::optional<GeneModel> gene_model(std::string_view id) const;

private:
    GffData data_;
    std::unordered_map<std::string, int> id_to_record_;
    std::unordered_map<std::string, std::vector<int>> id_to_records_;
    std::unordered_map<std::string, std::vector<int>> gene_lookup_;
    std::unordered_map<std::string, std::vector<int>> parents_by_child_id_;
    std::unordered_map<std::string, std::vector<int>> children_by_parent_id_;

    explicit AnnotationIndex(GffData data);
};

}  // namespace gffsub

#endif  // GFFSUB_ANNOTATION_HPP
