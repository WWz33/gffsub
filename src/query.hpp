#ifndef QUERY_HPP
#define QUERY_HPP

#include "annotation.hpp"
#include "output.hpp"
#include "query_summary.hpp"
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
    std::string feature_type;
    bool apply_type_filter = true;
    std::vector<std::string> output_attrs;
    std::string summary_format;
};

struct QueryResult {
    GffData records;
    std::vector<SummaryRow> summary_rows;
    bool emit_summary = false;
};

// Run a query against an index. Returns matched records and optional summary rows.
QueryResult query(const AnnotationIndex& index, const QueryParams& params);

// Print query result to the given stream based on params.summary_format.
void print_query_result(std::ostream& out, const QueryResult& result, const QueryParams& params);

}  // namespace gffsub

#endif  // QUERY_HPP
