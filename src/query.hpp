#ifndef QUERY_HPP
#define QUERY_HPP

#include "annotation.hpp"
#include "output.hpp"
#include "record.hpp"
#include "region.hpp"

#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace gffsub {

struct QueryParams {
    std::vector<std::string> ids;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    std::optional<Region> region;
    std::optional<Region> nearest_region;
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;
    std::string type;
    bool apply_type_filter = true;
};

struct QueryResult {
    GffData records;
};

// Run a query against an index. Returns matched records.
QueryResult query(const AnnotationIndex& index, const QueryParams& params);

}  // namespace gffsub

#endif  // QUERY_HPP
