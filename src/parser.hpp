#ifndef GFFSUB_PARSER_HPP
#define GFFSUB_PARSER_HPP

#include "input_source.hpp"
#include "record.hpp"

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace gffsub {

int parse_file(const std::string& filename, GffData& data, InputFormat format);
// Parse records out of data.buffer (already filled).
int parse_content(GffData& data, InputFormat format);
// Read all of stdin into data.buffer, sniff the format from the content,
// then parse. Returns 0 on success; format_out receives the sniffed format.
int parse_stdin(GffData& data, InputFormat& format_out);
// Decode column 9 into ordered (key, decoded value) pairs. This is the one
// column-9 decoder: GFF3 `tag=value` with %XX decoding and the five
// comma-list tags split; GTF `key "value";` with quote-aware splitting,
// backslash unescaping, and Parent split on commas. Repeated keys keep every
// pair, in input order. GTF identity keys (ID, Parent, gene_id,
// transcript_id) are %XX-decoded; other GTF values pass through verbatim.
std::vector<std::pair<std::string, std::string>> decode_attributes(
    std::string_view column9, InputFormat format);
// GFF3 key -> all values. Wrapper over decode_attributes for callers that
// want the map shape.
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);
// Every value of one attribute key, format-aware. The record's synthesized
// fields back the four identity keys when column 9 does not carry them.
std::vector<std::string> attr_values(const GffRecord& rec, std::string_view key);
// All Parent values of a record. GFF3 Parent=a,b and a quoted GTF
// Parent "a,b" both yield a and b; rec.parent_id is included as a fallback.
std::vector<std::string> parent_ids(const GffRecord& rec);

}  // namespace gffsub

#endif  // GFFSUB_PARSER_HPP
