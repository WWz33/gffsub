#ifndef GFFSUB_GTF_PARSER_HPP
#define GFFSUB_GTF_PARSER_HPP

#include "record.hpp"

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace gffsub {

// Extract a `key "value";` attribute from a GTF column-9 string.
// Handles escaped quotes in the value. Returns nullopt when absent.
std::optional<std::string> extract_quoted_value(std::string_view attrs, std::string_view key);

// Parse all `key "value";` pairs in original order, values unescaped.
// Fragments without a quoted value are skipped.
std::vector<std::pair<std::string, std::string>> parse_gtf_attributes(std::string_view attrs);

// Rewrite a GTF record's column 9 as GFF3 `tag=value` pairs, synthesizing
// ID=/Parent= from the resolved gene_id/transcript_id (AGAT gxf2gxf semantics).
// Values are URL-escaped per the GFF3 spec.
std::string gtf_attrs_to_gff3(const GffRecord& rec);

void apply_gtf_attributes(GffRecord& rec);

}  // namespace gffsub

#endif  // GFFSUB_GTF_PARSER_HPP
