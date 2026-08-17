#include "cli.hpp"
#include "expr_parser.hpp"

#include <cctype>
#include <cstring>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <string_view>

namespace gffsub {

// --- helpers ---

static std::string trim_copy(std::string_view value) {
    size_t start = 0;
    while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start]))) {
        ++start;
    }
    size_t end = value.size();
    while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }
    return std::string{value.substr(start, end - start)};
}

std::vector<std::string> split_attr_keys_cli(std::string_view keys) {
    std::vector<std::string> result;
    size_t pos = 0;
    while (pos <= keys.size()) {
        const size_t comma = keys.find(',', pos);
        const size_t end = (comma == std::string_view::npos) ? keys.size() : comma;
        const auto key = trim_copy(keys.substr(pos, end - pos));
        if (!key.empty()) {
            result.push_back(key);
        }
        if (comma == std::string_view::npos) {
            break;
        }
        pos = comma + 1;
    }
    return result;
}

static std::optional<char> parse_strand_filter(std::string_view value) {
    if (value.size() != 1) return std::nullopt;
    const char strand = value[0];
    if (strand == '+' || strand == '-' || strand == '.' || strand == '?') return strand;
    return std::nullopt;
}

static std::optional<char> parse_phase_filter(std::string_view value) {
    if (value.size() != 1) return std::nullopt;
    const char phase = value[0];
    if (phase == '0' || phase == '1' || phase == '2' || phase == '.') return phase;
    return std::nullopt;
}

static std::optional<std::optional<double>> parse_score_filter(std::string_view value) {
    if (value == ".") return std::optional<double>{};
    if (value.empty()) return std::nullopt;
    try {
        size_t parsed = 0;
        const double score = std::stod(std::string{value}, &parsed);
        if (parsed == value.size() && std::isfinite(score)) {
            return std::optional<double>{score};
        }
    } catch (const std::exception&) {
    }
    return std::nullopt;
}

AnnotationIndex load_index(const std::string& path) {
    return AnnotationIndex::from_file(path);
}

std::optional<std::vector<std::string>> load_id_list_file(const std::string& path) {
    std::ifstream in{path};
    if (!in.is_open()) {
        return std::nullopt;
    }
    std::vector<std::string> ids;
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto id = trim_copy(line);
        if (!id.empty() && id[0] != '#') {
            ids.push_back(id);
        }
    }
    return ids;
}

std::optional<QueryParams> build_query_params(const CliArgs& a) {
    QueryParams q;
    q.ids = a.ids;
    if (!a.id_list_file.empty()) {
        auto ids = load_id_list_file(a.id_list_file);
        if (!ids) {
            std::cerr << "Error: cannot open " << a.id_list_file << '\n';
            return std::nullopt;
        }
        q.ids.insert(q.ids.end(), std::make_move_iterator(ids->begin()), std::make_move_iterator(ids->end()));
    }
    q.name = a.name;
    q.attr_filters = a.attr_filters;
    if (!a.region_str.empty()) {
        q.region = parse_region(a.region_str);
        if (!q.region) {
            std::cerr << "Error: invalid region format " << a.region_str << '\n';
            return std::nullopt;
        }
    }
    if (!a.nearest_region_str.empty()) {
        q.nearest_region = parse_region(a.nearest_region_str);
        if (!q.nearest_region) {
            std::cerr << "Error: invalid nearest region format " << a.nearest_region_str << '\n';
            return std::nullopt;
        }
    }
    q.feature_type = a.feature;
    q.include_children = a.include_children;
    q.include_parents = a.include_parents;
    q.include_model = a.include_model;
    return q;
}

// --- argparse ---

std::optional<CliArgs> parse_cli_args(int argc, char* argv[], bool& help_requested) {
    CliArgs args;
    help_requested = false;

    enum {
        OPT_ID_LIST,
        OPT_SOURCE,
        OPT_SCORE,
        OPT_STRAND_FILTER,
        OPT_PHASE,
        OPT_GREP,
        OPT_GREP_REGEX,
        OPT_GREP_FILE,
        OPT_GREP_FIELD,
        OPT_GREP_FILE_REGEX,
        OPT_INCLUDE_EXPR,
        OPT_EXCLUDE_EXPR,
        OPT_INVERT_MATCH,
        OPT_IGNORE_CASE,
        OPT_VERSION
    };
    static struct option long_options[] = {
        {"id",            required_argument, nullptr, 'i'},
        {"ids",           required_argument, nullptr, OPT_ID_LIST},
        {"id-list",       required_argument, nullptr, OPT_ID_LIST},
        {"name",          required_argument, nullptr, 'n'},
        {"where",         required_argument, nullptr, 'w'},
        {"attr",          required_argument, nullptr, 'w'},
        {"grep",          required_argument, nullptr, OPT_GREP},
        {"grep-regex",    required_argument, nullptr, OPT_GREP_REGEX},
        {"grep-file",     required_argument, nullptr, OPT_GREP_FILE},
        {"grep-field",    required_argument, nullptr, OPT_GREP_FIELD},
        {"grep-file-regex", no_argument,     nullptr, OPT_GREP_FILE_REGEX},
        {"include-expr",  required_argument, nullptr, OPT_INCLUDE_EXPR},
        {"exclude-expr",  required_argument, nullptr, OPT_EXCLUDE_EXPR},
        {"invert-match",  no_argument,       nullptr, OPT_INVERT_MATCH},
        {"ignore-case",   no_argument,       nullptr, OPT_IGNORE_CASE},
        {"summary",       no_argument,       nullptr, 's'},
        {"parents",       no_argument,       nullptr, 'p'},
        {"include-parents", no_argument,      nullptr, 'p'},
        {"model",         no_argument,       nullptr, 'm'},
        {"gene-model",    no_argument,       nullptr, 'm'},
        {"nearest",       required_argument, nullptr, 'N'},
        {"nearest-gene",  required_argument, nullptr, 'N'},
        {"up",            required_argument, nullptr, 'u'},
        {"upstream",      required_argument, nullptr, 'u'},
        {"down",          required_argument, nullptr, 'D'},
        {"downstream",    required_argument, nullptr, 'D'},
        {"strand-aware",  no_argument,       nullptr, 'a'},
        {"seqid",         required_argument, nullptr, 'S'},
        {"source",        required_argument, nullptr, OPT_SOURCE},
        {"score",         required_argument, nullptr, OPT_SCORE},
        {"strand",        required_argument, nullptr, OPT_STRAND_FILTER},
        {"phase",         required_argument, nullptr, OPT_PHASE},
        {"children",      no_argument,       nullptr, 'C'},
        {"include-children", no_argument,     nullptr, 'C'},
        {"region",        required_argument, nullptr, 'r'},
        {"bed",           required_argument, nullptr, 'b'},
        {"feature",       required_argument, nullptr, 'f'},
        {"type",          required_argument, nullptr, 'f'},
        {"longest",       no_argument,       nullptr, 'L'},
        {"threads",       required_argument, nullptr, '@'},
        {"format",        required_argument, nullptr, 't'},
        {"output-format", required_argument, nullptr, 't'},
        {"output",        required_argument, nullptr, 'o'},
        {"help",          no_argument,       nullptr, 'h'},
        {"version",       no_argument,       nullptr, OPT_VERSION},
        {nullptr,         0,                 nullptr, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "r:b:f:CL@:t:o:hI:E:vi:n:w:spmN:u:D:aS:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'i': args.ids.emplace_back(optarg); break;
            case OPT_ID_LIST: args.id_list_file = optarg; break;
            case 'n': args.name = optarg; break;
            case 'w': {
                const std::string value{optarg};
                const auto equal_pos = value.find('=');
                if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value.size()) {
                    std::cerr << "Error: --" << long_options[option_index].name << " expects KEY=VALUE\n";
                    return std::nullopt;
                }
                args.attr_filters.emplace_back(value.substr(0, equal_pos), value.substr(equal_pos + 1));
                break;
            }
            case OPT_GREP:
            case OPT_GREP_REGEX: {
                const auto parsed = parse_field_pattern(optarg);
                if (!parsed) {
                    std::cerr << "Error: --" << long_options[option_index].name << " expects FIELD:PATTERN\n";
                    return std::nullopt;
                }
                GrepFilter filter;
                filter.field = parsed->first;
                filter.pattern = parsed->second;
                filter.use_regex = opt == OPT_GREP_REGEX;
                args.grep_filters.push_back(std::move(filter));
                break;
            }
            case OPT_GREP_FILE:
                args.grep_file = optarg;
                break;
            case OPT_GREP_FIELD:
                args.grep_field = optarg;
                break;
            case OPT_GREP_FILE_REGEX:
                args.grep_file_regex = true;
                break;
            case OPT_INCLUDE_EXPR:
            case 'I': {
                std::string error;
                if (!parse_expr_filters(optarg, args.include_expr_filters, error)) {
                    std::cerr << "Error: invalid include expression: " << error << '\n';
                    return std::nullopt;
                }
                break;
            }
            case OPT_EXCLUDE_EXPR:
            case 'E': {
                std::string error;
                if (!parse_expr_filters(optarg, args.exclude_expr_filters, error)) {
                    std::cerr << "Error: invalid exclude expression: " << error << '\n';
                    return std::nullopt;
                }
                break;
            }
            case OPT_INVERT_MATCH:
            case 'v':
                args.invert_grep = true;
                break;
            case OPT_IGNORE_CASE:
                args.ignore_case = true;
                break;
            case 's':
                args.summary = true;
                break;
            case 'p': args.include_parents = true; break;
            case 'm': args.include_model = true; break;
            case 'N': args.nearest_region_str = optarg; break;
            case 'C': args.include_children = true; break;
            case 'u': args.upstream_arg = optarg; break;
            case 'D': args.downstream_arg = optarg; break;
            case 'a': args.strand_aware = true; break;
            case 'S': args.seqid_filter = optarg; break;
            case OPT_SOURCE: args.source_filter = optarg; break;
            case OPT_SCORE: {
                args.score_filter = parse_score_filter(optarg);
                if (!args.score_filter) {
                    std::cerr << "Error: --score expects a finite floating point number or .\n";
                    return std::nullopt;
                }
                break;
            }
            case OPT_STRAND_FILTER: {
                args.strand_filter = parse_strand_filter(optarg);
                if (!args.strand_filter) {
                    std::cerr << "Error: --strand expects one of +, -, ., ?\n";
                    return std::nullopt;
                }
                break;
            }
            case OPT_PHASE: {
                args.phase_filter = parse_phase_filter(optarg);
                if (!args.phase_filter) {
                    std::cerr << "Error: --phase expects one of 0, 1, 2, .\n";
                    return std::nullopt;
                }
                break;
            }
            case 'r':
                args.region_str = optarg;
                if (!parse_region(optarg)) {
                    std::cerr << "Error: invalid region format " << optarg << '\n';
                    return std::nullopt;
                }
                break;
            case 'b': args.bed_file = optarg; break;
            case 'f': args.feature = optarg; break;
            case 'L': args.do_longest = true; break;
            case '@': {
                args.threads_set = true;
                // stoul accepts a leading '-' and wraps, so reject it explicitly.
                if (optarg[0] == '-') {
                    std::cerr << "Error: --threads must be a non-negative integer\n";
                    return std::nullopt;
                }
                size_t t = 0;
                try {
                    size_t pos = 0;
                    t = std::stoul(optarg, &pos);
                    if (pos != std::strlen(optarg)) throw std::invalid_argument{""};
                } catch (const std::exception&) {
                    std::cerr << "Error: --threads must be a non-negative integer\n";
                    return std::nullopt;
                }
                if (t == 0) t = 1;
                if (t > 256) t = 256;
                args.num_threads = t;
                break;
            }
            case 't': {
                const std::string fmt_str = optarg;
                if (fmt_str == "gff3") args.format = OutputFormat::GFF3;
                else if (fmt_str == "gtf" || fmt_str == "gtf2") args.format = OutputFormat::GTF2;
                else if (fmt_str == "gtf3") args.format = OutputFormat::GTF3;
                else if (fmt_str == "bed") args.format = OutputFormat::BED;
                else {
                    std::cerr << "Error: unknown output format " << fmt_str << '\n'
                              << "Supported formats: gff3, gtf2, gtf3, bed\n";
                    return std::nullopt;
                }
                break;
            }
            case 'o': args.output_file = optarg; break;
            case 'h': usage(argv[0]); help_requested = true; return std::nullopt;
            case OPT_VERSION:
                std::cout << "gffsub " << kVersion << '\n';
                help_requested = true;
                return std::nullopt;
            default: usage(argv[0]); return std::nullopt;
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        return std::nullopt;
    }

    args.input_file = argv[optind];

    // --- post-parse validation ---

    if (!args.grep_file.empty()) {
        if (args.grep_field.empty()) {
            std::cerr << "Error: --grep-file requires --grep-field\n";
            return std::nullopt;
        }
        std::ifstream in{args.grep_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << args.grep_file << '\n';
            return std::nullopt;
        }
        std::string line;
        while (std::getline(in, line)) {
            const auto pattern = trim_copy(line);
            if (pattern.empty()) continue;
            GrepFilter filter;
            filter.field = args.grep_field;
            filter.pattern = pattern;
            filter.use_regex = args.grep_file_regex;
            args.grep_filters.push_back(std::move(filter));
        }
    } else if (!args.grep_field.empty() || args.grep_file_regex) {
        std::cerr << "Error: --grep-field and --grep-file-regex require --grep-file\n";
        return std::nullopt;
    }
    if (args.invert_grep && args.grep_filters.empty()) {
        std::cerr << "Error: --invert-match requires --grep, --grep-regex, or --grep-file\n";
        return std::nullopt;
    }
    {
        std::string error;
        if (!compile_filter_regexes(args.grep_filters, args.include_expr_filters,
                                    args.exclude_expr_filters, args.ignore_case, error)) {
            std::cerr << "Error: " << error << '\n';
            return std::nullopt;
        }
    }

    const bool has_query_style_selector = !args.ids.empty() || !args.id_list_file.empty() ||
        !args.name.empty() || !args.attr_filters.empty() || !args.nearest_region_str.empty();
    if ((args.include_children || args.include_parents || args.include_model) && !has_query_style_selector) {
        std::cerr << "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest\n";
        return std::nullopt;
    }

    return args;
}

}  // namespace gffsub
