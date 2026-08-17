#ifndef GFFSUB_STRING_UTILS_HPP
#define GFFSUB_STRING_UTILS_HPP

#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace gffsub {

// Lowercase a string view, returning a new std::string.
std::string to_lower(std::string_view sv);

// Split a string on a single-character delimiter. Returns all fields
// including empty ones (like split on '\t' for sparse GFF columns).
std::vector<std::string> split_line(const std::string& line, char delimiter);

// Parse comma-separated list with optional ^ prefix for exclusion mode.
// e.g. "exon,CDS" -> {"exon","CDS"}, exclude=false
//      "^gene" -> {"gene"}, exclude=true
std::unordered_set<std::string> parse_list(std::string_view sv, bool& exclude);

}  // namespace gffsub

#endif  // GFFSUB_STRING_UTILS_HPP
