// Benchmark: current parse (std::string per field) vs arena parse
// (string_view into one big buffer) on the same file. No changes to src/.
// Usage: ./bench_parse <gff3> [iterations]

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace {

struct RecString {
    std::string seqid, source, type, score_raw, attr_raw;
    int64_t start = 0, end = 0;
    char strand = '.', phase = '.';
    std::optional<std::string> id, parent_id;
};

struct RecView {
    std::string_view seqid, source, type, score_raw, attr_raw;
    int64_t start = 0, end = 0;
    char strand = '.', phase = '.';
    std::optional<std::string_view> id, parent_id;
};

std::vector<std::string_view> split_sv(std::string_view line, char delim) {
    std::vector<std::string_view> cols;
    cols.reserve(9);
    size_t start = 0;
    while (true) {
        auto pos = line.find(delim, start);
        if (pos == std::string_view::npos) { cols.emplace_back(line.substr(start)); break; }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

std::optional<std::string_view> extract_attr_sv(std::string_view attrs, std::string_view key) {
    size_t pos = 0;
    while (pos < attrs.size()) {
        size_t pair_end = attrs.find(';', pos);
        if (pair_end == std::string_view::npos) pair_end = attrs.size();
        const auto pair = attrs.substr(pos, pair_end - pos);
        pos = (pair_end < attrs.size()) ? pair_end + 1 : attrs.size();
        if (pair.empty()) continue;
        const size_t eq = pair.find('=');
        if (eq == std::string_view::npos || eq == 0) continue;
        if (pair.substr(0, eq) == key) {
            const auto value = pair.substr(eq + 1);
            if (value.empty()) return std::nullopt;
            const size_t comma = value.find(',');
            const auto first = (comma == std::string_view::npos) ? value : value.substr(0, comma);
            if (first.empty()) return std::nullopt;
            return first;
        }
    }
    return std::nullopt;
}

int64_t parse_int(std::string_view s) {
    return std::stoll(std::string{s});
}

// --- Current style: std::string fields, copy out of line ---
std::vector<RecString> parse_string(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    std::vector<RecString> recs;
    recs.reserve(1u << 20);
    std::string line;
    line.reserve(256);
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("##FASTA", 0) == 0) break;
        auto cols = split_sv(line, '\t');
        if (cols.size() < 9) continue;
        RecString r;
        r.seqid = cols[0]; r.source = cols[1]; r.type = cols[2];
        r.start = parse_int(cols[3]); r.end = parse_int(cols[4]);
        r.score_raw = cols[5];
        r.strand = cols[6].empty() ? '.' : cols[6][0];
        r.phase = cols[7].empty() ? '.' : cols[7][0];
        r.attr_raw = cols[8];
        if (auto v = extract_attr_sv(r.attr_raw, "ID")) r.id = std::string{*v};
        if (auto v = extract_attr_sv(r.attr_raw, "Parent")) r.parent_id = std::string{*v};
        recs.push_back(std::move(r));
    }
    return recs;
}

// --- Arena style: read whole file, string_view into buffer ---
std::pair<std::string, std::vector<RecView>> parse_arena(const std::string& path) {
    std::ifstream f(path, std::ios::binary | std::ios::ate);
    const auto size = f.tellg();
    f.seekg(0);
    std::string buf(static_cast<size_t>(size), '\0');
    f.read(buf.data(), size);

    std::vector<RecView> recs;
    recs.reserve(1u << 20);
    std::string_view sv{buf};
    size_t pos = 0;
    while (pos < sv.size()) {
        size_t eol = sv.find('\n', pos);
        if (eol == std::string_view::npos) eol = sv.size();
        auto line = sv.substr(pos, eol - pos);
        pos = eol + 1;
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
        if (line.empty() || line[0] == '#') continue;
        if (line.substr(0, 7) == "##FASTA") break;
        auto cols = split_sv(line, '\t');
        if (cols.size() < 9) continue;
        RecView r;
        r.seqid = cols[0]; r.source = cols[1]; r.type = cols[2];
        r.start = parse_int(cols[3]); r.end = parse_int(cols[4]);
        r.score_raw = cols[5];
        r.strand = cols[6].empty() ? '.' : cols[6][0];
        r.phase = cols[7].empty() ? '.' : cols[7][0];
        r.attr_raw = cols[8];
        r.id = extract_attr_sv(r.attr_raw, "ID");
        r.parent_id = extract_attr_sv(r.attr_raw, "Parent");
        recs.push_back(r);
    }
    return {std::move(buf), std::move(recs)};
}

template <typename F>
double time_ms(F&& fn, int iters) {
    double best = 1e18;
    for (int i = 0; i < iters; ++i) {
        const auto t0 = std::chrono::steady_clock::now();
        fn();
        const auto t1 = std::chrono::steady_clock::now();
        best = std::min(best, std::chrono::duration<double, std::milli>(t1 - t0).count());
    }
    return best;
}

}  // namespace

int main(int argc, char* argv[]) {
    if (argc < 2) { std::cerr << "usage: bench_parse <gff3> [iters]\n"; return 1; }
    const std::string path = argv[1];
    const int iters = argc > 2 ? std::atoi(argv[2]) : 3;

    size_t n_string = 0, n_arena = 0;
    const double t_string = time_ms([&] {
        auto recs = parse_string(path);
        n_string = recs.size();
    }, iters);
    const double t_arena = time_ms([&] {
        auto [buf, recs] = parse_arena(path);
        n_arena = recs.size();
    }, iters);

    std::printf("records: string=%zu arena=%zu\n", n_string, n_arena);
    std::printf("string parse: %8.1f ms\n", t_string);
    std::printf("arena  parse: %8.1f ms\n", t_arena);
    std::printf("speedup:      %8.2fx\n", t_string / t_arena);
    return 0;
}
