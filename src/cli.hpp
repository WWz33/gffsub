#ifndef GFFSUB_CLI_HPP
#define GFFSUB_CLI_HPP

#include "annotation.hpp"
#include "cli_usage.hpp"
#include "expr_parser.hpp"
#include "output.hpp"
#include "parser.hpp"
#include "query.hpp"
#include "record.hpp"
#include "region.hpp"
#include "selector_filter.hpp"
#include "subset.hpp"
#include "version.hpp"
#include "window.hpp"

#include <optional>
#include <string>
#include <vector>

namespace gffsub {

struct CliArgs {
    std::string input_file;

    // query-style selectors
    std::vector<std::string> ids;
    std::string id_list_file;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    std::string nearest_region_str;
    std::string region_str;

    // window shortcut
    std::string upstream_arg;
    std::string downstream_arg;
    bool strand_aware = false;

    // tree expansion
    bool include_children = false;
    bool include_parents = false;
    bool include_model = false;

    // subset filters
    std::string seqid_filter;
    std::string source_filter;
    std::optional<std::optional<double>> score_filter;
    std::optional<char> strand_filter;
    std::optional<char> phase_filter;
    std::string bed_file;
    std::string feature;
    bool do_longest = false;
    size_t num_threads = 6;
    bool threads_set = false;

    // grep / expr
    std::vector<GrepFilter> grep_filters;
    std::string grep_file;
    std::string grep_field;
    bool grep_file_regex = false;
    std::vector<ExprNode> include_expr_filters;
    std::vector<ExprNode> exclude_expr_filters;
    bool invert_grep = false;
    bool ignore_case = false;

    // summary / output
    std::vector<std::string> output_attrs;
    std::string summary_format;
    OutputFormat format = OutputFormat::GFF3;
    std::string output_file;
};

// Parse command-line arguments for the default (non-subcommand) path.
// Returns nullopt on error or help requested; sets help_requested when help was printed.
std::optional<CliArgs> parse_cli_args(int argc, char* argv[], bool& help_requested);

// Helpers used by cli and main
AnnotationIndex load_index(const std::string& path);
std::vector<std::string> split_attr_keys_cli(std::string_view keys);

// Load one ID per non-empty line from a file. Returns nullopt on open failure.
std::optional<std::vector<std::string>> load_id_list_file(const std::string& path);

// Build QueryParams from CliArgs (loads id_list_file, parses regions). Returns nullopt on error.
std::optional<QueryParams> build_query_params(const CliArgs& a);

}  // namespace gffsub

#endif  // GFFSUB_CLI_HPP
