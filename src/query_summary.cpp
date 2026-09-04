#include "query_summary.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>
#include <vector>

namespace gffsub {

namespace {

std::string tsv_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char ch : s) {
        if (ch == '\t') out += "\\t";
        else if (ch == '\n') out += "\\n";
        else if (ch == '\r') out += "\\r";
        else out.push_back(ch);
    }
    return out;
}

// Integer when integral, one decimal otherwise. Lengths are non-negative.
std::string fmt_len(double v) {
    if (v < 0) v = 0;
    const int64_t whole = static_cast<int64_t>(v);
    if (static_cast<double>(whole) == v) return std::to_string(whole);
    return std::to_string(whole) + "." + std::to_string(static_cast<int64_t>((v - whole) * 10));
}

// Plain median: middle element, or mean of the two middles.
double median_of(const std::vector<int64_t>& sorted, size_t lo, size_t hi) {
    const size_t n = hi - lo;
    if (n == 0) return 0;
    if (n % 2 == 1) return static_cast<double>(sorted[lo + n / 2]);
    return (sorted[lo + n / 2 - 1] + sorted[lo + n / 2]) / 2.0;
}

// Tukey hinges, as used by seqkit stats: Q1/Q3 are medians of the
// lower/upper half (median included in both halves for odd n).
double quantile_tukey(const std::vector<int64_t>& sorted, double p) {
    const size_t n = sorted.size();
    if (n == 0) return 0;
    if (p == 0.5) return median_of(sorted, 0, n);
    const size_t half = (n + 1) / 2;
    if (p < 0.5) return median_of(sorted, 0, half);
    return median_of(sorted, n - half, n);
}

// Union bp: sort intervals by start, merge overlapping or book-ended
// runs (next.start <= cur.end), sum each merged span (bedtools merge).
int64_t union_bp(const std::vector<std::pair<int64_t, int64_t>>& intervals) {
    if (intervals.empty()) return 0;
    std::vector<std::pair<int64_t, int64_t>> sorted{intervals};
    std::sort(sorted.begin(), sorted.end());
    int64_t total = 0;
    int64_t cur_start = sorted.front().first;
    int64_t cur_end = sorted.front().second;
    for (size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i].first <= cur_end + 1) {
            if (sorted[i].second > cur_end) cur_end = sorted[i].second;
        } else {
            total += cur_end - cur_start + 1;
            cur_start = sorted[i].first;
            cur_end = sorted[i].second;
        }
    }
    total += cur_end - cur_start + 1;
    return total;
}

struct Stats {
    int64_t count = 0;
    int64_t sum_len = 0;
    int64_t min_len = 0;
    int64_t max_len = 0;
    // Union bp within this group's seqid. Cross-seqid coverage is the
    // sum of per-seqid unions, never a merge across seqids.
    int64_t coverage = -1;
    std::vector<int64_t> lens;
    std::vector<std::pair<int64_t, int64_t>> intervals;

    void add(int64_t start, int64_t end) {
        const int64_t len = end - start + 1;
        if (count == 0 || len < min_len) min_len = len;
        if (len > max_len) max_len = len;
        ++count;
        sum_len += len;
        lens.push_back(len);
        intervals.push_back({start, end});
    }

    void finalize() {
        if (coverage < 0) coverage = union_bp(intervals);
    }
};

void print_stats_row(std::ostream& out, const std::string& seqid, const std::string& type,
                     Stats& s) {
    std::vector<int64_t> sorted{s.lens};
    std::sort(sorted.begin(), sorted.end());
    s.finalize();
    out << tsv_escape(seqid) << '\t' << tsv_escape(type) << '\t'
        << s.count << '\t' << s.sum_len << '\t'
        << s.min_len << '\t' << fmt_len(static_cast<double>(s.sum_len) / s.count) << '\t'
        << s.max_len << '\t'
        << fmt_len(quantile_tukey(sorted, 0.25)) << '\t'
        << fmt_len(quantile_tukey(sorted, 0.5)) << '\t'
        << fmt_len(quantile_tukey(sorted, 0.75)) << '\t'
        << s.coverage << '\n';
}

}  // namespace

std::string record_id(const GffRecord& rec) {
    if (rec.id) return *rec.id;
    if (rec.gene_id) return *rec.gene_id;
    if (rec.transcript_id) return *rec.transcript_id;
    return "";
}

void print_summary(std::ostream& out, const std::vector<GffRecord>& records) {
    // (seqid, type) -> stats
    std::map<std::pair<std::string, std::string>, Stats> agg;
    for (const auto& rec : records) {
        agg[{rec.seqid, rec.type}].add(rec.start, rec.end);
    }

    std::unordered_set<std::string> seqids;
    for (const auto& [key, _] : agg) seqids.insert(key.first);
    const bool has_all = seqids.size() > 1;

    out << "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len\tQ1\tQ2\tQ3\tcoverage\n";
    for (auto& [key, s] : agg) {
        print_stats_row(out, key.first, key.second, s);
    }

    if (has_all) {
        std::map<std::string, Stats> by_type;
        for (auto& [key, s] : agg) {
            auto& e = by_type[key.second];
            e.count += s.count;
            e.sum_len += s.sum_len;
            if (e.count == s.count || s.min_len < e.min_len) e.min_len = s.min_len;
            if (s.max_len > e.max_len) e.max_len = s.max_len;
            e.lens.insert(e.lens.end(), s.lens.begin(), s.lens.end());
            // cross-seqid coverage = sum of per-seqid unions
            e.coverage = (e.coverage < 0 ? 0 : e.coverage) + std::max<int64_t>(s.coverage, 0);
        }
        for (auto& [type, s] : by_type) {
            print_stats_row(out, "all", type, s);
        }
    }
}

}  // namespace gffsub
