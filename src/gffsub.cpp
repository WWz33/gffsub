#include "cli.hpp"
#include "filter.hpp"
#include "output.hpp"
#include "parser.hpp"
#include "query.hpp"
#include "query_summary.hpp"
#include "string_utils.hpp"
#include "subset.hpp"
#include "window.hpp"

#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_set>
#include <utility>
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

// Describe why a path could not be read, keeping the operating system's
// reason. A missing file, a permission problem, and a corrupt gzip all used
// to print the same sentence.
std::string describe_path_error(const std::string& path) {
    if (path != "-") {
        struct stat st{};
        if (::stat(path.c_str(), &st) != 0) return std::strerror(errno);
        if (S_ISDIR(st.st_mode)) return "is a directory";
        if (::access(path.c_str(), R_OK) != 0) return std::strerror(errno);
    }
    return {};
}

// Split `--flag=value` into the flag and its attached value. A `-`/`--`
// argument without '=' returns the argument unchanged with no attached value.
std::pair<std::string, std::optional<std::string>> split_long_option(const std::string& arg) {
    if (arg.size() > 2 && arg[0] == '-' && arg[1] == '-') {
        const auto eq = arg.find('=');
        if (eq != std::string::npos && eq > 2) {
            return {arg.substr(0, eq), arg.substr(eq + 1)};
        }
    }
    return {arg, std::nullopt};
}

std::optional<OutputFormat> parse_output_format(const std::string& value) {
    if (value == "gff3") return OutputFormat::GFF3;
    if (value == "gtf" || value == "gtf3") return OutputFormat::GTF3;
    if (value == "gtf25") return OutputFormat::GTF25;
    if (value == "gtf22") return OutputFormat::GTF22;
    if (value == "bed") return OutputFormat::BED;
    return std::nullopt;
}

struct OutputOptions {
    std::string input;
    std::string output;  // empty or "-" = stdout
    OutputFormat format = OutputFormat::GFF3;
    std::vector<std::string> out_attrs;
    bool summary = false;
    // Built from the full annotation when the printed data is a subset.
    // Null means build it from the data being printed.
    const GtfModel* gtf_model = nullptr;
};

// Owns the output stream: opens the file, refuses to truncate the input,
// flushes, closes, and reports a failed write. Every command writes through
// this one path.
class OutputSink {
public:
    OutputSink(const std::string& input, const std::string& path) : path_(path) {
        if (path.empty() || path == "-") return;
        // Input gzip is transparent; output gzip is not implemented. Reject
        // the suffix instead of writing plain text into a .gz name.
        if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0) {
            error_ = "gzip output is not supported; write plain text and compress separately: " + path;
            return;
        }
        // Opening the output truncates it. If it is the same file as a plain
        // (mmap'd) input, reading pages past the new EOF faults. Same dev/ino
        // covers symlinks and hard links.
        if (!input.empty() && input != "-") {
            struct stat in_st{};
            struct stat out_st{};
            if (::stat(input.c_str(), &in_st) == 0 &&
                ::stat(path.c_str(), &out_st) == 0 &&
                in_st.st_dev == out_st.st_dev && in_st.st_ino == out_st.st_ino) {
                error_ = "input and output are the same file: " + input;
                return;
            }
        }
        file_.open(path);
        if (!file_.is_open()) {
            error_ = "cannot open " + path + ": " + std::strerror(errno);
            return;
        }
        out_ = &file_;
    }

    bool ok() const { return error_.empty(); }
    const std::string& error() const { return error_; }
    std::ostream& stream() { return *out_; }

    bool finish() {
        if (!ok()) return false;
        out_->flush();
        if (!out_->good()) {
            error_ = "cannot write " + display_path();
            return false;
        }
        if (file_.is_open()) {
            file_.close();
            if (file_.fail()) {
                error_ = "cannot write " + path_;
                return false;
            }
        }
        return true;
    }

private:
    std::string display_path() const {
        return (path_.empty() || path_ == "-") ? std::string{"stdout"} : path_;
    }

    std::ofstream file_;
    std::ostream* out_ = &std::cout;
    std::string path_;
    std::string error_;
};

int write_output(const GffData& data, const OutputOptions& opts) {
    OutputSink sink(opts.input, opts.output);
    if (!sink.ok()) {
        std::cerr << "Error: " << sink.error() << '\n';
        return 1;
    }
    std::ostream& out = sink.stream();
    if (opts.summary) {
        std::vector<GffRecord> kept_records;
        for (const auto& rec : data) {
            if (rec.kept) kept_records.push_back(rec);
        }
        print_summary(out, kept_records);
    } else {
        switch (opts.format) {
            case OutputFormat::GFF3: print_gff3(out, data, opts.out_attrs); break;
            case OutputFormat::GTF3:
            case OutputFormat::GTF25:
            case OutputFormat::GTF22: {
                // Build the gene/transcript/stop-codon model from the whole
                // annotation, then print only the kept rows.
                const GtfModel model = opts.gtf_model ? *opts.gtf_model
                                                      : build_gtf_model(data);
                print_gtf(out, data, model, opts.format, opts.out_attrs);
                break;
            }
            case OutputFormat::BED: print_bed(out, data); break;
        }
    }
    if (!sink.finish()) {
        std::cerr << "Error: " << sink.error() << '\n';
        return 1;
    }
    return 0;
}

}  // namespace

// --- query subcommand (gffsub query <file> [opts]) ---

static int run_query_subcommand(int argc, char* argv[], const char* prog) {
    if (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        query_usage(std::cout, prog);
        return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--version") {
        std::cout << "gffsub " << kVersion << '\n';
        return 0;
    }
    if (argc < 2) {
        query_usage(std::cerr, prog);
        return 1;
    }

    const std::string input_file = argv[1];
    QueryParams params;
    std::vector<std::string> id_list_files;
    bool summary = false;
    OutputOptions out_opts;
    out_opts.input = input_file;

    for (int i = 2; i < argc; ++i) {
        const auto option = split_long_option(argv[i]);
        const std::string& arg = option.first;
        const std::optional<std::string>& attached = option.second;
        auto require_value = [&](const char* name) -> std::optional<std::string> {
            if (attached) {
                if (attached->empty()) {
                    std::cerr << "Error: " << name << " requires a non-empty value\n";
                    return std::nullopt;
                }
                return *attached;
            }
            if (i + 1 >= argc) {
                std::cerr << "Error: " << name << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            // Reject empty selectors like the main parser does: they would
            // otherwise match nothing and exit 0.
            if (argv[i][0] == '\0') {
                std::cerr << "Error: " << name << " requires a non-empty value\n";
                return std::nullopt;
            }
            return std::string{argv[i]};
        };

        if (arg == "-i" || arg == "--id") {
            auto value = require_value("-i");
            if (!value) return 1;
            params.ids.push_back(*value);
        } else if (arg == "-n" || arg == "--name") {
            auto value = require_value("-n");
            if (!value) return 1;
            params.name = *value;
        } else if (arg == "--ids" || arg == "--id-list") {
            auto value = require_value(arg.c_str());
            if (!value) return 1;
            id_list_files.push_back(*value);
        } else if (arg == "-r" || arg == "--region") {
            auto value = require_value("-r");
            if (!value) return 1;
            params.region = parse_region(*value);
            if (!params.region) {
                std::cerr << "Error: invalid region format " << *value << '\n';
                return 1;
            }
        } else if (arg == "-N" || arg == "--nearest" || arg == "--nearest-gene") {
            auto value = require_value("-N");
            if (!value) return 1;
            params.nearest_region = parse_region(*value);
            if (!params.nearest_region) {
                std::cerr << "Error: invalid nearest region format " << *value << '\n';
                return 1;
            }
        } else if (arg == "-t" || arg == "--type") {
            auto value = require_value("--type");
            if (!value) return 1;
            // repeatable: append with comma
            if (params.type.empty()) params.type = *value;
            else params.type += "," + *value;
        } else if (arg == "-w" || arg == "--where" || arg == "--attr") {
            auto value = require_value("-w");
            if (!value) return 1;
            const auto equal_pos = value->find('=');
            if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value->size()) {
                std::cerr << "Error: " << arg << " expects KEY=VALUE\n";
                return 1;
            }
            params.attr_filters.emplace_back(value->substr(0, equal_pos), value->substr(equal_pos + 1));
        } else if (arg == "-s" || arg == "--summary") {
            summary = true;
        } else if (arg == "-C" || arg == "--children" || arg == "--include-children") {
            params.include_children = true;
        } else if (arg == "-p" || arg == "--parents" || arg == "--include-parents") {
            params.include_parents = true;
        } else if (arg == "-m" || arg == "--model" || arg == "--gene-model") {
            params.include_model = true;
        } else if (arg == "-o" || arg == "--output") {
            auto value = require_value("-o");
            if (!value) return 1;
            out_opts.output = *value;
        } else if (arg == "-f" || arg == "--format") {
            auto value = require_value("-f");
            if (!value) return 1;
            const auto fmt = parse_output_format(*value);
            if (!fmt) {
                std::cerr << "Error: unknown output format " << *value << '\n';
                return 1;
            }
            out_opts.format = *fmt;
        } else if (arg == "--out-attrs") {
            auto value = require_value("--out-attrs");
            if (!value) return 1;
            for (const auto tag : split_line(*value, ',')) {
                if (!tag.empty()) out_opts.out_attrs.emplace_back(tag);
            }
        } else if (arg == "-h" || arg == "--help") {
            query_usage(std::cout, prog);
            return 0;
        } else {
            std::cerr << "Error: unknown query option " << arg << '\n';
            query_usage(std::cerr, prog);
            return 1;
        }
    }

    if ((params.include_children || params.include_parents || params.include_model) &&
        params.ids.empty() && id_list_files.empty() && params.name.empty() &&
        params.attr_filters.empty() && !params.nearest_region) {
        std::cerr << "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest\n";
        return 1;
    }

    if (params.ids.empty() && id_list_files.empty() && params.name.empty() &&
        params.attr_filters.empty() && !params.nearest_region && !params.region) {
        std::cerr << "Error: query requires a selector: --id, --ids, --name, --where, --region, or --nearest\n";
        query_usage(std::cerr, prog);
        return 1;
    }

    for (const auto& path : id_list_files) {
        auto ids = load_id_list_file(path);
        if (!ids) {
            std::cerr << "Error: cannot open " << path << '\n';
            return 1;
        }
        if (ids->empty()) {
            std::cerr << "Error: no IDs in " << path << '\n';
            return 1;
        }
        params.ids.insert(params.ids.end(), std::make_move_iterator(ids->begin()), std::make_move_iterator(ids->end()));
    }

    const auto index = try_load_index(input_file);
    if (!index) return 1;
    if (!summary && out_opts.format != OutputFormat::BED && index->source_format() == InputFormat::BED) {
        std::cerr << "Error: BED input is not converted to GFF3; use --summary or convert elsewhere\n";
        return 1;
    }
    const auto result = query(*index, params);
    out_opts.summary = summary;
    GtfModel gtf_model;
    if (!summary && out_opts.format != OutputFormat::GFF3 &&
        out_opts.format != OutputFormat::BED) {
        gtf_model = build_gtf_model(index->records());
        out_opts.gtf_model = &gtf_model;
    }
    return write_output(result.records, out_opts);
}

// --- window subcommand (gffsub window <file> [opts]) ---

static int run_window_subcommand(int argc, char* argv[], const char* prog) {
    if (argc >= 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        window_usage(std::cout, prog);
        return 0;
    }
    if (argc >= 2 && std::string(argv[1]) == "--version") {
        std::cout << "gffsub " << kVersion << '\n';
        return 0;
    }
    if (argc < 2) {
        window_usage(std::cerr, prog);
        return 1;
    }

    const std::string input_file = argv[1];
    WindowParams params;
    OutputOptions out_opts;
    out_opts.input = input_file;

    for (int i = 2; i < argc; ++i) {
        const auto option = split_long_option(argv[i]);
        const std::string& arg = option.first;
        const std::optional<std::string>& attached = option.second;
        auto require_value = [&](const char* name) -> std::optional<std::string> {
            if (attached) {
                if (attached->empty()) {
                    std::cerr << "Error: " << name << " requires a non-empty value\n";
                    return std::nullopt;
                }
                return *attached;
            }
            if (i + 1 >= argc) {
                std::cerr << "Error: " << name << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            // Reject empty selectors like the main parser does: they would
            // otherwise match nothing and exit 0.
            if (argv[i][0] == '\0') {
                std::cerr << "Error: " << name << " requires a non-empty value\n";
                return std::nullopt;
            }
            return std::string{argv[i]};
        };

        if (arg == "--id" || arg == "-i") {
            auto value = require_value("--id");
            if (!value) return 1;
            if (!params.id.empty()) {
                std::cerr << "Error: window requires exactly one --id\n";
                return 1;
            }
            params.id = *value;
        } else if (arg == "--up" || arg == "--upstream" || arg == "-u") {
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
        } else if (arg == "--down" || arg == "--downstream" || arg == "-D") {
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
        } else if (arg == "--strand-aware" || arg == "-a") {
            params.strand_aware = true;
        } else if (arg == "-o" || arg == "--output") {
            auto value = require_value("-o");
            if (!value) return 1;
            out_opts.output = *value;
        } else if (arg == "-f" || arg == "--format") {
            auto value = require_value("-f");
            if (!value) return 1;
            const auto fmt = parse_output_format(*value);
            if (!fmt) {
                std::cerr << "Error: unknown output format " << *value << '\n';
                return 1;
            }
            out_opts.format = *fmt;
        } else if (arg == "--out-attrs") {
            auto value = require_value("--out-attrs");
            if (!value) return 1;
            for (const auto tag : split_line(*value, ',')) {
                if (!tag.empty()) out_opts.out_attrs.emplace_back(tag);
            }
        } else if (arg == "-h" || arg == "--help") {
            window_usage(std::cout, prog);
            return 0;
        } else {
            std::cerr << "Error: unknown window option " << arg << '\n';
            window_usage(std::cerr, prog);
            return 1;
        }
    }

    if (params.id.empty()) {
        std::cerr << "Error: window requires --id\n";
        return 1;
    }

    const auto index = try_load_index(input_file);
    if (!index) return 1;
    if (out_opts.format != OutputFormat::BED && index->source_format() == InputFormat::BED) {
        std::cerr << "Error: BED input is not converted to GFF3; use -f bed\n";
        return 1;
    }
    try {
        const auto result = window(*index, params);
        GtfModel gtf_model;
        if (out_opts.format != OutputFormat::GFF3 && out_opts.format != OutputFormat::BED) {
            gtf_model = build_gtf_model(index->records());
            out_opts.gtf_model = &gtf_model;
        }
        return write_output(result, out_opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
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
    if (a.format != OutputFormat::BED && index->source_format() == InputFormat::BED) {
        std::cerr << "Error: BED input is not converted to GFF3; use -f bed\n";
        return 1;
    }
    try {
        const auto result = window(*index, wparams);
        OutputOptions out_opts;
        out_opts.input = a.input_file;
        out_opts.output = a.output_file;
        out_opts.format = a.format;
        out_opts.out_attrs = a.out_attrs;
        GtfModel gtf_model;
        if (a.format != OutputFormat::GFF3 && a.format != OutputFormat::BED) {
            gtf_model = build_gtf_model(index->records());
            out_opts.gtf_model = &gtf_model;
        }
        return write_output(result, out_opts);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }
}

// True when CliArgs has no subset filters (only query-style selectors are set).
bool no_subset_filters(const CliArgs& a) {
    return a.seqid_filter.empty() && a.source_filter.empty() &&
           !a.score_filter && !a.strand_filter && !a.phase_filter && a.bed_file.empty() &&
           !a.do_longest && !a.threads_set && a.longest_type.empty() &&
           a.grep_filters.empty() && a.grep_file.empty() && a.grep_field.empty() && !a.grep_file_regex &&
           a.include_expr_filters.empty() && a.exclude_expr_filters.empty() && !a.invert_grep && !a.ignore_case;
}

// Build SubsetParams from CliArgs. Region is parsed separately (already validated).
SubsetParams build_subset_params(const CliArgs& a) {
    SubsetParams s;
    if (!a.region_str.empty()) {
        s.region = parse_region(a.region_str);
    }
    if (!a.exclude_region_str.empty()) {
        s.exclude_region = parse_region(a.exclude_region_str);
    }
    s.bed_file = a.bed_file;
    s.bed_exclude = a.bed_exclude;
    s.bed_strand_mode = a.bed_strand_mode;
    s.seqid_filter = a.seqid_filter;
    s.source_filter = a.source_filter;
    s.score_filter = a.score_filter;
    s.strand_filter = a.strand_filter;
    s.phase_filter = a.phase_filter;
    s.type = a.type_filter;
    s.longest_type = a.longest_type;
    s.longest = a.do_longest;
    s.threads = a.num_threads;
    s.drop_orphans = a.drop_orphans;
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

    // Window shortcut: --up/--down/--strand-aware with exactly one --id.
    // The guard below is a blacklist over CliArgs: any new field MUST be
    // added here or window will silently accept (and ignore) its flag.
    if (!a.upstream_arg.empty() || !a.downstream_arg.empty() || a.strand_aware) {
        if (a.ids.size() != 1) {
            std::cerr << "Error: window shortcut requires exactly one --id\n";
            return 1;
        }
        if (!a.id_list_files.empty() || !a.name.empty() || !a.attr_filters.empty() || !a.nearest_region_str.empty() ||
            a.include_children || a.include_parents || a.include_model ||
            a.summary || !a.region_str.empty() || !a.bed_file.empty() ||
            !a.seqid_filter.empty() || !a.source_filter.empty() || a.score_filter || a.strand_filter ||
            a.phase_filter || !a.type_filter.empty() || a.do_longest || !a.longest_type.empty() ||
            !a.grep_filters.empty() || !a.grep_file.empty() || !a.grep_field.empty() || a.grep_file_regex ||
            !a.include_expr_filters.empty() || !a.exclude_expr_filters.empty() || a.invert_grep || a.ignore_case ||
            a.format != OutputFormat::GFF3 || !a.output_file.empty() || !a.sort_keys.empty() || a.sort_reverse ||
            a.threads_set || !a.exclude_region_str.empty() || a.bed_exclude || a.bed_strand_mode != 0 ||
            a.drop_orphans || !a.out_attrs.empty()) {
            std::cerr << "Error: window shortcut only supports --id, --up/--upstream, --down/--downstream, and --strand-aware\n";
            return 1;
        }
        return run_window_from_args(a);
    }

    const bool has_query_style_selector = !a.ids.empty() || !a.id_list_files.empty() ||
        !a.name.empty() || !a.attr_filters.empty() || !a.nearest_region_str.empty();

    // Pure query-style selector (no subset filters) → dispatch to query API
    const bool can_dispatch_to_query = !a.summary && no_subset_filters(a) &&
        a.format == OutputFormat::GFF3 && a.output_file.empty() && a.region_str.empty() &&
        a.sort_keys.empty() && !a.sort_reverse && a.exclude_region_str.empty() &&
        !a.drop_orphans && a.out_attrs.empty() && a.bed_strand_mode == 0 && !a.bed_exclude;

    if (has_query_style_selector && can_dispatch_to_query) {
        auto qparams = build_query_params(a);
        if (!qparams) return 1;

        const auto index = try_load_index(a.input_file);
        if (!index) return 1;
        if (index->source_format() == InputFormat::BED) {
            std::cerr << "Error: BED input is not converted to GFF3; use -f bed or -s\n";
            return 1;
        }
        const auto result = query(*index, *qparams);
        OutputOptions out_opts;
        out_opts.input = a.input_file;
        return write_output(result.records, out_opts);
    }

    // --- main subset path ---

    GffData data;
    InputFormat input_fmt = InputFormat::GFF3;
    if (a.input_file == "-") {
        if (parse_stdin(data, input_fmt) != 0) {
            std::cerr << "Error: cannot parse stdin\n";
            return 1;
        }
    } else {
        input_fmt = infer_input_format(a.input_file);
        if (parse_file(a.input_file, data, input_fmt) != 0) {
            const auto why = describe_path_error(a.input_file);
            std::cerr << "Error: cannot parse " << a.input_file;
            if (!why.empty()) std::cerr << " (" << why << ")";
            std::cerr << '\n';
            return 1;
        }
    }

    // BED is accepted for BED output and summary only; converting it to
    // GFF3/GTF is not supported.
    if (input_fmt == InputFormat::BED && !a.summary && a.format != OutputFormat::BED) {
        std::cerr << "Error: BED input is not converted to GFF3/GTF; use -f bed or -s\n";
        return 1;
    }

    // Query-style selectors: mark selected lines via query() API
    if (has_query_style_selector) {
        auto qparams = build_query_params(a);
        if (!qparams) return 1;
        qparams->apply_type_filter = false;
        qparams->region.reset();  // region is handled by subset() filter, not query selector

        const auto index = AnnotationIndex::from_data(data);
        const auto result = query(index, *qparams);

        std::unordered_set<int> selected_lines;
        for (const auto& rec : result.records.records) {
            selected_lines.insert(rec.line_idx);
        }
        for (auto& rec : data.records) {
            rec.kept = selected_lines.count(rec.line_idx) > 0;
        }
    }

    try {
        subset(data, build_subset_params(a));
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    // Sort after filtering: only order changes, membership does not.
    if (a.sort_keys.empty() && a.sort_reverse) {
        std::cerr << "Error: --reverse requires --sort\n";
        return 1;
    }
    if (!a.sort_keys.empty() && a.summary) {
        // Summary output order is determined by (seqid, type) aggregation,
        // not record order; sorting would have no observable effect.
        std::cerr << "Error: --sort has no effect with --summary\n";
        return 1;
    }
    if (!a.sort_keys.empty()) {
        if (sort_records(data, a.sort_keys, a.sort_reverse) != 0) {
            std::cerr << "Error: unknown sort key in \"" << a.sort_keys
                      << "\" (valid: seqid, natural-seqid, start, end, length, type)\n";
            return 1;
        }
    }

    // Output
    OutputOptions out_opts;
    out_opts.input = a.input_file;
    out_opts.output = a.output_file;
    out_opts.format = a.format;
    out_opts.out_attrs = a.out_attrs;
    out_opts.summary = a.summary;
    return write_output(data, out_opts);
}
