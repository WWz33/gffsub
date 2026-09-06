#ifndef GFFSUB_STRING_UTILS_HPP
#define GFFSUB_STRING_UTILS_HPP

#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

namespace gffsub {

// Lowercase a string view, returning a new std::string.
std::string to_lower(std::string_view sv);

// Split a string view on a single-character delimiter. Returns views into
// the input, including empty fields (like split on '\t' for sparse GFF
// columns). The input must outlive the returned views.
std::vector<std::string_view> split_line(std::string_view line, char delimiter);

// Natural-order compare: embedded digit runs compare by value, so "Gm2"
// sorts before "Gm10". Returns <0/0/>0 like strcmp.
int natural_compare(std::string_view a, std::string_view b);

// URL-decode a GFF3 attribute value per spec: %XX -> byte.
std::string url_decode(std::string_view input);

// Parse comma-separated list with optional ^ prefix for exclusion mode.
// e.g. "exon,CDS" -> {"exon","CDS"}, exclude=false
//      "^gene" -> {"gene"}, exclude=true
std::unordered_set<std::string> parse_list(std::string_view sv, bool& exclude);

}  // namespace gffsub

#endif  // GFFSUB_STRING_UTILS_HPP
