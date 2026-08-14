#ifndef GFFSUB_PARSER_HPP
#define GFFSUB_PARSER_HPP

#include "record.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

int parse_file(const std::string& filename, GffData& data, IdIndex& idx, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);

}  // namespace gffsub

#endif  // GFFSUB_PARSER_HPP
