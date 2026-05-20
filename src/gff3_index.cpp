#include "gff3.hpp"

#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace gffsub {

static std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs) {
    std::unordered_map<std::string, std::vector<std::string>> parsed;
    size_t pos = 0;
    while (pos < attrs.size()) {
        const size_t key_end = attrs.find('=', pos);
        if (key_end == std::string_view::npos) {
            break;
        }

        const std::string key{attrs.substr(pos, key_end - pos)};
        const size_t value_start = key_end + 1;
        size_t value_end = attrs.find(';', value_start);
        if (value_end == std::string_view::npos) {
            value_end = attrs.size();
        }

        const std::string value{attrs.substr(value_start, value_end - value_start)};
        if (!key.empty() && !value.empty()) {
            size_t part_start = 0;
            while (part_start <= value.size()) {
                size_t part_end = value.find(',', part_start);
                if (part_end == std::string::npos) {
                    part_end = value.size();
                }
                const std::string part = value.substr(part_start, part_end - part_start);
                if (!part.empty()) {
                    parsed[key].push_back(part);
                }
                if (part_end == value.size()) {
                    break;
                }
                part_start = part_end + 1;
            }
        }

        pos = (value_end < attrs.size()) ? value_end + 1 : attrs.size();
    }
    return parsed;
}

AnnotationIndex AnnotationIndex::from_gff3(const std::string& path) {
    GffData data;
    IdIndex idx;
    if (parse_file(path, data, idx, InputFormat::GFF3) != 0) {
        throw std::runtime_error("cannot parse GFF3 file: " + path);
    }
    return AnnotationIndex{std::move(data)};
}

AnnotationIndex::AnnotationIndex(GffData data) : data_(std::move(data)) {
    for (int i = 0; i < static_cast<int>(data_.records.size()); ++i) {
        const auto& rec = data_.records[i];
        if (rec.id) {
            id_to_record_.emplace(*rec.id, i);
        }

        if (rec.type != "gene") {
            continue;
        }

        auto add_gene_key = [&](const std::string& key) {
            if (!key.empty()) {
                gene_lookup_[key].push_back(i);
            }
        };

        if (rec.id) {
            add_gene_key(*rec.id);
        }
        if (rec.gene_id) {
            add_gene_key(*rec.gene_id);
        }

        const auto attrs = parse_attributes(rec.attr_raw);
        for (const char* key : {"Name", "gene_id", "locus_tag", "Alias"}) {
            const auto it = attrs.find(key);
            if (it == attrs.end()) {
                continue;
            }
            for (const auto& value : it->second) {
                add_gene_key(value);
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

std::optional<GffRecord> AnnotationIndex::find_gene(std::string_view id) const {
    const auto it = gene_lookup_.find(std::string{id});
    if (it == gene_lookup_.end() || it->second.empty()) {
        return std::nullopt;
    }
    return data_.records[it->second.front()];
}

}  // namespace gffsub
