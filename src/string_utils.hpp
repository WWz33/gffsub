#ifndef GFFSUB_STRING_UTILS_HPP
#define GFFSUB_STRING_UTILS_HPP

#include <string>
#include <string_view>
#include <vector>

namespace gffsub {

// Lowercase a string view, returning a new std::string.
std::string to_lower(std::string_view sv);

// Split a string on a single-character delimiter. Returns all fields
// including empty ones (like split on '\t' for sparse GFF columns).
std::vector<std::string> split_line(const std::string& line, char delimiter);

}  // namespace gffsub

#endif  // GFFSUB_STRING_UTILS_HPP
