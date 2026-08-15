#include "cli.hpp"
#include "filter.hpp"
#include "output.hpp"
#include "parser.hpp"
#include "query.hpp"
#include "subset.hpp"
#include "window.hpp"

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_set>
#include <vector>

using namespace gffsub;

namespace {

// load_index throws on unreadable files; convert to a clean error + exit.
std::optional<AnnotationIndex> try_load_index(const std::string& path) {
    try {
        return load_index(path);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return std::nullopt;
    }
}

}  // namespace

// --- query subcommand (gffsub query <file> [opts]) ---

static int run_query_subcommand(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        query_usage(prog);
        return 0;
    }
    if (argc < 2) {
        query_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    QueryParams params;
    std::string id_list_file;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            params.ids.push_back(*value);
        } else if (arg == "--name") {
            auto value = require_value("--name");
            if (!value) return 1;
            params.name = *value;
        } else if (arg == "--ids" || arg == "--id-list") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            id_list_file = *value;
        } else if (arg == "--region") {
            auto value = require_value("--region");
            if (!value) return 1;
            params.region = parse_region(*value);
            if (!params.region) {
                std::cerr << "Error: invalid region format " << *value << '\n';
                return 1;
            }
        } else if (arg == "--nearest" || arg == "--nearest-gene") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            params.nearest_region = parse_region(*value);
            if (!params.nearest_region) {
                std::cerr << "Error: invalid nearest region format " << *value << '\n';
                return 1;
            }
        } else if (arg == "--type") {
            auto value = require_value("--type");
            if (!value) return 1;
            params.feature_type = *value;
        } else if (arg == "--where" || arg == "--attr") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            const auto equal_pos = value->find('=');
            if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value->size()) {
                std::cerr << "Error: " << arg << " expects KEY=VALUE\n";
                return 1;
            }
            params.attr_filters.emplace_back(value->substr(0, equal_pos), value->substr(equal_pos + 1));
        } else if (arg == "--output-attrs" || arg == "--out-attrs" || arg == "--attrs") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            const auto keys = split_attr_keys_cli(*value);
            if (keys.empty()) {
                std::cerr << "Error: " << arg << " expects a comma-separated list of keys\n";
                return 1;
            }
            params.output_attrs.insert(params.output_attrs.end(), keys.begin(), keys.end());
        } else if (arg == "--summary" || arg == "--summary-format") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            params.summary_format = *value;
            if (params.summary_format != "tsv" && params.summary_format != "json") {
                std::cerr << "Error: " << arg << " expects tsv or json\n";
                return 1;
            }
        } else if (arg == "-C" || arg == "--children" || arg == "--include-children") {
            params.include_children = true;
        } else if (arg == "--parents" || arg == "--include-parents") {
            params.include_parents = true;
        } else if (arg == "--model" || arg == "--gene-model") {
            params.include_model = true;
        } else if (arg == "-h" || arg == "--help") {
            query_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown query option " << arg << '\n';
            query_usage(prog);
            return 1;
        }
    }

    if ((params.include_children || params.include_parents || params.include_model) &&
        params.ids.empty() && id_list_file.empty() && params.name.empty() &&
        params.attr_filters.empty() && !params.nearest_region) {
        std::cerr << "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest\n";
        return 1;
    }

    if (params.ids.empty() && id_list_file.empty() && params.name.empty() &&
        params.attr_filters.empty() && !params.nearest_region && !params.region) {
        std::cerr << "Error: query requires a selector: --id, --ids, --name, --where, --region, or --nearest\n";
        query_usage(prog);
        return 1;
    }

    if (!id_list_file.empty()) {
        auto ids = load_id_list_file(id_list_file);
        if (!ids) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        params.ids.insert(params.ids.end(), std::make_move_iterator(ids->begin()), std::make_move_iterator(ids->end()));
    }

    const auto index = try_load_index(input_file);
    if (!index) return 1;
    const auto result = query(*index, params);
    print_query_result(std::cout, result, params);
    return 0;
}

// --- window subcommand (gffsub window <file> [opts]) ---

static int run_window_subcommand(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        window_usage(prog);
        return 0;
    }
    if (argc < 2) {
        window_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    WindowParams params;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            if (!params.id.empty()) {
                std::cerr << "Error: window requires exactly one --id\n";
                return 1;
            }
            params.id = *value;
        } else if (arg == "--up" || arg == "--upstream") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            try {
                size_t pos = 0;
                params.upstream = std::stoll(*value, &pos);
                if (pos != value->size()) throw std::invalid_argument{""};
            } catch (const std::exception&) {
                std::cerr << "Error: " << arg << " must be a non-negative integer\n";
                return 1;
            }
            if (params.upstream < 0) {
                std::cerr << "Error: " << arg << " must be non-negative\n";
                return 1;
            }
        } else if (arg == "--down" || arg == "--downstream") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            try {
                size_t pos = 0;
                params.downstream = std::stoll(*value, &pos);
                if (pos != value->size()) throw std::invalid_argument{""};
            } catch (const std::exception&) {
                std::cerr << "Error: " << arg << " must be a non-negative integer\n";
                return 1;
            }
            if (params.downstream < 0) {
                std::cerr << "Error: " << arg << " must be non-negative\n";
                return 1;
            }
        } else if (arg == "--strand-aware") {
            params.strand_aware = true;
        } else if (arg == "-h" || arg == "--help") {
            window_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown window option " << arg << '\n';
            window_usage(prog);
            return 1;
        }
    }

    if (params.id.empty()) {
        std::cerr << "Error: window requires --id\n";
        return 1;
    }

    const auto index = try_load_index(input_file);
    if (!index) return 1;
    try {
        const auto result = window(*index, params);
        print_gff3(std::cout, result);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}

// --- main ---

namespace {

int run_window_from_args(const CliArgs& a) {
    WindowParams wparams;
    wparams.id = a.ids[0];
    auto parse_window_value = [](const std::string& arg, const char* name) -> std::optional<int64_t> {
        try {
            size_t pos = 0;
            const int64_t v = std::stoll(arg, &pos);
            if (pos != arg.size()) {
                std::cerr << "Error: " << name << " must be a non-negative integer\n";
                return std::nullopt;
            }
            if (v < 0) {
                std::cerr << "Error: " << name << " must be non-negative\n";
                return std::nullopt;
            }
            return v;
        } catch (const std::exception&) {
            std::cerr << "Error: " << name << " must be a non-negative integer\n";
            return std::nullopt;
        }
    };
    if (!a.upstream_arg.empty()) {
        auto v = parse_window_value(a.upstream_arg, "--up");
        if (!v) return 1;
        wparams.upstream = *v;
    }
    if (!a.downstream_arg.empty()) {
        auto v = parse_window_value(a.downstream_arg, "--down");
        if (!v) return 1;
        wparams.downstream = *v;
    }
    wparams.strand_aware = a.strand_aware;

    const auto index = try_load_index(a.input_file);
    if (!index) return 1;
    try {
        const auto result = window(*index, wparams);
        print_gff3(std::cout, result);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
    return 0;
}

// True when CliArgs has no subset filters (only query-style selectors are set).
bool no_subset_filters(const CliArgs& a) {
    return a.seqid_filter.empty() && a.source_filter.empty() &&
           !a.score_filter && !a.strand_filter && !a.phase_filter && a.bed_file.empty() &&
           !a.do_longest && !a.threads_set &&
           a.grep_filters.empty() && a.grep_file.empty() && a.grep_field.empty() && !a.grep_file_regex &&
           a.include_expr_filters.empty() && a.exclude_expr_filters.empty() && !a.invert_grep && !a.ignore_case;
}

// Build SubsetParams from CliArgs. Region is parsed separately (already validated).
SubsetParams build_subset_params(const CliArgs& a) {
    SubsetParams s;
    if (!a.region_str.empty()) {
        s.region = parse_region(a.region_str);
    }
    s.bed_file = a.bed_file;
    s.seqid_filter = a.seqid_filter;
    s.source_filter = a.source_filter;
    s.score_filter = a.score_filter;
    s.strand_filter = a.strand_filter;
    s.phase_filter = a.phase_filter;
    s.feature = a.feature;
    s.longest = a.do_longest;
    s.threads = a.num_threads;
    s.grep_filters = a.grep_filters;
    s.include_exprs = a.include_expr_filters;
    s.exclude_exprs = a.exclude_expr_filters;
    s.invert_grep = a.invert_grep;
    return s;
}

}  // namespace

int main(int argc, char* argv[]) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc > 1 && std::string(argv[1]) == "query") {
        return run_query_subcommand(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "window") {
        return run_window_subcommand(argc - 1, argv + 1, argv[0]);
    }

    bool help_requested = false;
    auto args = parse_cli_args(argc, argv, help_requested);
    if (!args) return help_requested ? 0 : 1;
    const auto& a = *args;

    // Window shortcut: --up/--down/--strand-aware with exactly one --id
    if (!a.upstream_arg.empty() || !a.downstream_arg.empty() || a.strand_aware) {
        if (a.ids.size() != 1) {
            std::cerr << "Error: window shortcut requires exactly one --id\n";
            return 1;
        }
        if (!a.id_list_file.empty() || !a.name.empty() || !a.attr_filters.empty() || !a.nearest_region_str.empty() ||
            a.include_children || a.include_parents || a.include_model ||
            !a.output_attrs.empty() || !a.summary_format.empty() || !a.region_str.empty() || !a.bed_file.empty() ||
            !a.seqid_filter.empty() || !a.source_filter.empty() || a.score_filter || a.strand_filter ||
            a.phase_filter || !a.feature.empty() || a.do_longest ||
            !a.grep_filters.empty() || !a.grep_file.empty() || !a.grep_field.empty() || a.grep_file_regex ||
            !a.include_expr_filters.empty() || !a.exclude_expr_filters.empty() || a.invert_grep || a.ignore_case ||
            a.format != OutputFormat::GFF3 || !a.output_file.empty()) {
            std::cerr << "Error: window shortcut only supports --id, --up/--upstream, --down/--downstream, and --strand-aware\n";
            return 1;
        }
        return run_window_from_args(a);
    }

    const bool has_query_style_selector = !a.ids.empty() || !a.id_list_file.empty() ||
        !a.name.empty() || !a.attr_filters.empty() || !a.nearest_region_str.empty();

    // Summary or pure query-style selector (no subset filters) → dispatch to query API
    const bool can_dispatch_to_query = no_subset_filters(a) &&
        a.format == OutputFormat::GFF3 && a.output_file.empty() && a.region_str.empty();

    if (!a.summary_format.empty() || !a.output_attrs.empty()) {
        if (!can_dispatch_to_query) {
            std::cerr << "Error: --summary/--summary-format/--out-attrs only supports query-style selectors; "
                      << "do not combine with --seqid, --source, --score, --strand, --phase, --bed, --longest, --threads, --format/--output-format, or --output\n";
            return 1;
        }
    }

    if (!a.summary_format.empty() || !a.output_attrs.empty() ||
        (has_query_style_selector && can_dispatch_to_query)) {
        auto qparams = build_query_params(a);
        if (!qparams) return 1;

        const auto index = try_load_index(a.input_file);
        if (!index) return 1;
        const auto result = query(*index, *qparams);
        print_query_result(std::cout, result, *qparams);
        return 0;
    }

    // --- main subset path ---

    GffData data;
    IdIndex idx;
    const InputFormat input_fmt = infer_input_format(a.input_file);
    if (parse_file(a.input_file, data, idx, input_fmt) != 0) {
        std::cerr << "Error: cannot parse " << a.input_file << '\n';
        return 1;
    }

    // Query-style selectors: mark selected lines via query() API
    if (has_query_style_selector) {
        auto qparams = build_query_params(a);
        if (!qparams) return 1;
        qparams->apply_type_filter = false;
        qparams->summary_format.clear();
        qparams->output_attrs.clear();
        qparams->region.reset();  // region is handled by subset() filter, not query selector

        const auto index = AnnotationIndex::from_data(GffData{data});
        const auto result = query(index, *qparams);

        std::unordered_set<int> selected_lines;
        for (const auto& rec : result.records.records) {
            selected_lines.insert(rec.line_idx);
        }
        for (auto& rec : data.records) {
            rec.kept = selected_lines.count(rec.line_idx) > 0;
        }
    }

    subset(data, idx, build_subset_params(a));

    // Output
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!a.output_file.empty()) {
        out_file.open(a.output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: cannot open " << a.output_file << '\n';
            return 1;
        }
        out = &out_file;
    }

    switch (a.format) {
        case OutputFormat::GFF3: print_gff3(*out, data); break;
        case OutputFormat::GTF2: print_gtf(*out, data, a.format); break;
        case OutputFormat::GTF3: print_gtf3(*out, data); break;
        case OutputFormat::BED:  print_bed(*out, data); break;
    }

    return 0;
}
