#ifndef GFFSUB_PARSER_HPP
#define GFFSUB_PARSER_HPP

#include "input_source.hpp"
#include "record.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

int parse_file(const std::string& filename, GffData& data, InputFormat format);
// Parse records out of data.buffer (already filled).
int parse_content(GffData& data, InputFormat format);
// Read all of stdin into data.buffer, sniff the format from the content,
// then parse. Returns 0 on success; format_out receives the sniffed format.
int parse_stdin(GffData& data, InputFormat& format_out);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);

}  // namespace gffsub

#endif  // GFFSUB_PARSER_HPP
