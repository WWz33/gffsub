#ifndef GFFSUB_PARSER_HPP
#define GFFSUB_PARSER_HPP

#include "record.hpp"

#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

int parse_file(const std::string& filename, GffData& data, InputFormat format);
std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs);
InputFormat infer_input_format(const std::string& path);

}  // namespace gffsub

#endif  // GFFSUB_PARSER_HPP
