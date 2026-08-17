// Regression smoke tests covering behaviors fixed in rounds 1-5 that had
// no prior coverage. Mirrors the conventions of cli_selector_smoke.cpp:
// standalone main, std::system calls against ./gffsub, expect helpers, and a
// "regression_smoke OK" line on success. argc==2 takes the executable path.

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "test_utils.hpp"

using test_utils::expect_command_failure;
using test_utils::read_file;
using test_utils::require_contains;
using test_utils::require_exit_one_with_error;
using test_utils::require_not_contains;
using test_utils::run_command;

// ---------------------------------------------------------------------------
// Fixture writers
// ---------------------------------------------------------------------------

static bool write_gtf_basic(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";\n"
        << "chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; gene_name \"GENE1\";\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
        << "chr1\tsrc\tCDS\t150\t170\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n";
    return true;
}

static bool write_gtf_unsorted(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    // Exon appears BEFORE its transcript; ID collision must not let the exon
    // masquerade as the transcript when querying --id T1.
    out << "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id \"G1\";\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
        << "chr1\tsrc\ttranscript\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
        << "chr1\tsrc\tCDS\t150\t170\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";\n";
    return true;
}

static bool write_gtf_mrna(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    // Input uses "mRNA"; --format gtf3 must rename it to "transcript".
    out << "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tgene_id \"G1\";\n"
        << "chr1\tsrc\tmRNA\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n";
    return true;
}

static bool write_multi_parent(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t1000\t.\t+\t.\tID=g1\n"
        << "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=t1;Parent=g1\n"
        << "chr1\tsrc\tmRNA\t500\t1000\t.\t+\t.\tID=t2;Parent=g1\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tID=ex_shared;Parent=t1,t2\n"
        << "chr1\tsrc\texon\t210\t260\t.\t+\t.\tID=ex_t1;Parent=t1\n"
        << "chr1\tsrc\texon\t610\t700\t.\t+\t.\tID=ex_t2;Parent=t2\n";
    return true;
}

static bool write_escaped_comma(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1;Note=a%2Cb\n"
        << "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g2;Note=a\n";
    return true;
}

static bool write_bad_coords(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t0\t100\t.\t+\t.\tID=g_start0\n"
        << "chr1\tsrc\tgene\t200\t100\t.\t+\t.\tID=g_start_gt_end\n"
        << "chr1\tsrc\tgene\t10abc\t100\t.\t+\t.\tID=g_trailing\n"
        << "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g_ok\n";
    return true;
}

static bool write_url_encoded(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1;gene_id=Gene%201\n"
        << "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=g2;gene_id=Other\n";
    return true;
}

static bool write_quote_attr(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=g1;Note=has\"quote\n";
    return true;
}

// ---------------------------------------------------------------------------
// Output cleanup
// ---------------------------------------------------------------------------

static void cleanup_outputs() {
    const char* files[] = {
        "regression_gtf_basic.gtf", "regression_gtf_unsorted.gtf",
        "regression_gtf_mrna.gtf", "regression_multi.gff3",
        "regression_esc.gff3", "regression_bad.gff3",
        "regression_url.gff3", "regression_quote.gff3",
        "reg_gtf_children.gff3", "reg_gtf_model.gff3",
        "reg_gtf_unsorted_id.gff3", "reg_gtf_out.gtf",
        "reg_gtf3_rename.gff3",
        "reg_multi_t1_children.gff3", "reg_multi_t2_children.gff3",
        "reg_multi_longest.gff3",
        "reg_esc_match.gff3", "reg_esc_nomatch.gff3",
        "reg_bad_coords.gff3",
        "reg_url_match.gff3",
        "reg_gtf_grep.gff3",
        "reg_gtf_summary.tsv", "reg_gene_summary.tsv",
        "reg_json.json", "reg_json_quote.json",
        "reg_err_missing.err", "reg_err_up.err", "reg_err_threads.err",
        "reg_err_region.err", "reg_err_query_no_selector.err"};
    for (const char* f : files) std::remove(f);
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// Group 1: GTF parent synthesis. --id T1 -C expands to transcript/exon/CDS;
// --model also returns the gene line.
static int test_gtf_parent_synthesis(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " query " + gtf + " --id T1 -C > reg_gtf_children.gff3") != 0 ||
        require_contains("reg_gtf_children.gff3", "transcript\t100\t500") != 0 ||
        require_contains("reg_gtf_children.gff3", "exon\t120\t180") != 0 ||
        require_contains("reg_gtf_children.gff3", "CDS\t150\t170") != 0 ||
        require_not_contains("reg_gtf_children.gff3", "type=\"gene\"") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gtf + " --id T1 --model > reg_gtf_model.gff3") != 0 ||
        require_contains("reg_gtf_model.gff3", "gene\t100\t500") != 0 ||
        require_contains("reg_gtf_model.gff3", "transcript\t100\t500") != 0 ||
        require_contains("reg_gtf_model.gff3", "exon\t120\t180") != 0 ||
        require_contains("reg_gtf_model.gff3", "CDS\t150\t170") != 0) {
        return 1;
    }
    return 0;
}

// Group 2: unsorted GTF (exon before transcript). --id T1 returns the
// transcript line, not the exon that shares transcript_id "T1".
static int test_gtf_no_id_collision(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " query " + gtf + " --id T1 > reg_gtf_unsorted_id.gff3") != 0 ||
        require_contains("reg_gtf_unsorted_id.gff3", "transcript\t100\t500") != 0 ||
        require_not_contains("reg_gtf_unsorted_id.gff3", "exon\t120\t180") != 0) {
        return 1;
    }
    return 0;
}

// Group 3: GTF output. Gene line has NO transcript_id; transcript/exon have
// both gene_id and transcript_id. --format gtf3 renames mRNA -> transcript.
static int test_gtf_output(const std::string& exe, const std::string& gtf_basic,
                           const std::string& gtf_mrna) {
    if (run_command(exe + " " + gtf_basic + " --name G1 --model --output-format gtf > reg_gtf_out.gtf") != 0 ||
        require_contains("reg_gtf_out.gtf", "##gtf-version 2") != 0 ||
        // Gene line: gene_id present, no transcript_id on the same line.
        require_contains("reg_gtf_out.gtf", "gene\t100\t500\t.\t+\t.\tgene_id \"G1\";") != 0 ||
        require_not_contains("reg_gtf_out.gtf", "gene\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id") != 0 ||
        // Transcript + exon lines carry both gene_id and transcript_id.
        require_contains("reg_gtf_out.gtf", "transcript\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";") != 0 ||
        require_contains("reg_gtf_out.gtf", "exon\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf_mrna + " --id T1 -C --output-format gtf3 > reg_gtf3_rename.gff3") != 0 ||
        require_contains("reg_gtf3_rename.gff3", "##gtf-version 2.2.1") != 0 ||
        require_contains("reg_gtf3_rename.gff3", "transcript\t100\t500") != 0 ||
        require_not_contains("reg_gtf3_rename.gff3", "mRNA") != 0) {
        return 1;
    }
    return 0;
}

// Group 4: GFF3 multi-parent. Parent=t1,t2 shared exon is a child of both;
// --longest keeps the shared child of the longest isoform and drops the
// shorter isoform's exclusive children.
static int test_multi_parent(const std::string& exe, const std::string& gff) {
    // t1 children include the shared exon and t1's own exon; not t2's.
    if (run_command(exe + " query " + gff + " --id t1 --children > reg_multi_t1_children.gff3") != 0 ||
        require_contains("reg_multi_t1_children.gff3", "ID=ex_shared") != 0 ||
        require_contains("reg_multi_t1_children.gff3", "ID=ex_t1") != 0 ||
        require_not_contains("reg_multi_t1_children.gff3", "ID=ex_t2") != 0) {
        return 1;
    }
    // t2 children include the shared exon and t2's own exon; not t1's.
    if (run_command(exe + " query " + gff + " --id t2 --children > reg_multi_t2_children.gff3") != 0 ||
        require_contains("reg_multi_t2_children.gff3", "ID=ex_shared") != 0 ||
        require_contains("reg_multi_t2_children.gff3", "ID=ex_t2") != 0 ||
        require_not_contains("reg_multi_t2_children.gff3", "ID=ex_t1") != 0) {
        return 1;
    }
    // --longest -C on the gene: t2 is the longest isoform (exon sum 150 > 110).
    // Shared exon kept, t1 and t1's exclusive exon dropped.
    if (run_command(exe + " " + gff + " --id g1 --longest -C > reg_multi_longest.gff3") != 0 ||
        require_contains("reg_multi_longest.gff3", "ID=g1") != 0 ||
        require_contains("reg_multi_longest.gff3", "ID=t2") != 0 ||
        require_contains("reg_multi_longest.gff3", "ID=ex_shared") != 0 ||
        require_contains("reg_multi_longest.gff3", "ID=ex_t2") != 0 ||
        require_not_contains("reg_multi_longest.gff3", "ID=t1") != 0 ||
        require_not_contains("reg_multi_longest.gff3", "ID=ex_t1") != 0) {
        return 1;
    }
    return 0;
}

// Group 5: escaped comma. Note=a%2Cb parses as ONE value "a,b".
// --where Note=a,b matches; --where Note=a does not match g1.
static int test_escaped_comma(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " --where 'Note=a,b' > reg_esc_match.gff3") != 0 ||
        require_contains("reg_esc_match.gff3", "ID=g1") != 0 ||
        require_not_contains("reg_esc_match.gff3", "ID=g2") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --where Note=a > reg_esc_nomatch.gff3") != 0 ||
        require_contains("reg_esc_nomatch.gff3", "ID=g2") != 0 ||
        require_not_contains("reg_esc_nomatch.gff3", "ID=g1") != 0) {
        return 1;
    }
    return 0;
}

// Group 6: coordinate validation. start=0, start>end, and trailing garbage
// lines are all skipped; only the well-formed line survives.
static int test_coord_validation(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " > reg_bad_coords.gff3") != 0 ||
        require_contains("reg_bad_coords.gff3", "ID=g_ok") != 0 ||
        require_not_contains("reg_bad_coords.gff3", "ID=g_start0") != 0 ||
        require_not_contains("reg_bad_coords.gff3", "ID=g_start_gt_end") != 0 ||
        require_not_contains("reg_bad_coords.gff3", "ID=g_trailing") != 0) {
        return 1;
    }
    return 0;
}

// Group 7: --where URL decode. Attribute gene_id=Gene%201 decodes to "Gene 1";
// --where 'gene_id=Gene 1' matches.
static int test_where_url_decode(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " --where 'gene_id=Gene 1' > reg_url_match.gff3") != 0 ||
        require_contains("reg_url_match.gff3", "ID=g1") != 0 ||
        require_not_contains("reg_url_match.gff3", "ID=g2") != 0) {
        return 1;
    }
    return 0;
}

// Group 8: summary output. -s on query produces TSV with per-record summary
// including child/exon/cds counts.
static int test_summary_scope(const std::string& exe, const std::string& gtf) {
    // Transcript match: -i T1 -s aggregates the transcript record.
    if (run_command(exe + " " + gtf + " -i T1 -s > reg_gtf_summary.tsv") != 0 ||
        require_contains("reg_gtf_summary.tsv", "seqid\ttype\tcount\tsum_len") != 0) {
        return 1;
    }
    {
        const auto text = read_file("reg_gtf_summary.tsv");
        if (text.find("all\t") != std::string::npos) {
            std::cerr << "unexpected all row in single-seqid summary\n";
            return 1;
        }
        if (text.find("chr1\ttranscript\t1\t401\t401\t401\t401") == std::string::npos) {
            std::cerr << "transcript summary row wrong\n";
            return 1;
        }
    }
    // Gene match: -i G1 -s aggregates the gene record.
    if (run_command(exe + " " + gtf + " -i G1 -s > reg_gene_summary.tsv") != 0) {
        return 1;
    }
    {
        const auto text = read_file("reg_gene_summary.tsv");
        if (text.find("chr1\tgene\t1\t401\t401\t401\t401") == std::string::npos) {
            std::cerr << "gene summary row wrong\n";
            return 1;
        }
    }
    return 0;
}

// Group 9: error handling. Each invocation must exit 1 with an "Error:"
// message printed to stderr (no crash).
static int test_error_handling(const std::string& exe, const std::string& gtf) {
    if (require_exit_one_with_error(
            exe + " /nonexistent/file.gff3 > /dev/null 2> reg_err_missing.err",
            "reg_err_missing.err", "cannot parse") != 0) {
        return 1;
    }
    if (require_exit_one_with_error(
            exe + " " + gtf + " --id T1 --up 50abc > /dev/null 2> reg_err_up.err",
            "reg_err_up.err", "non-negative integer") != 0) {
        return 1;
    }
    if (require_exit_one_with_error(
            exe + " " + gtf + " --id T1 --threads -1 > /dev/null 2> reg_err_threads.err",
            "reg_err_threads.err", "non-negative integer") != 0) {
        return 1;
    }
    if (require_exit_one_with_error(
            exe + " " + gtf + " -r chr1:1-100abc > /dev/null 2> reg_err_region.err",
            "reg_err_region.err", "invalid region format") != 0) {
        return 1;
    }
    if (require_exit_one_with_error(
            exe + " query " + gtf + " > /dev/null 2> reg_err_query_no_selector.err",
            "reg_err_query_no_selector.err", "query requires a selector") != 0) {
        return 1;
    }
    return 0;
}

// Group 10: GTF attribute access. --grep attr.gene_name:X matches GTF lines
// carrying that attribute.
static int test_gtf_attr_access(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " " + gtf + " --grep 'attr.gene_name:GENE1' > reg_gtf_grep.gff3") != 0 ||
        require_contains("reg_gtf_grep.gff3", "transcript\t100\t500") != 0 ||
        require_contains("reg_gtf_grep.gff3", "gene_name=GENE1") != 0) {
        return 1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: regression_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string gtf_basic{"regression_gtf_basic.gtf"};
    const std::string gtf_unsorted{"regression_gtf_unsorted.gtf"};
    const std::string gtf_mrna{"regression_gtf_mrna.gtf"};
    const std::string multi{"regression_multi.gff3"};
    const std::string esc{"regression_esc.gff3"};
    const std::string bad{"regression_bad.gff3"};
    const std::string url{"regression_url.gff3"};
    const std::string quote{"regression_quote.gff3"};

    if (!write_gtf_basic(gtf_basic) || !write_gtf_unsorted(gtf_unsorted) ||
        !write_gtf_mrna(gtf_mrna) || !write_multi_parent(multi) ||
        !write_escaped_comma(esc) || !write_bad_coords(bad) ||
        !write_url_encoded(url) || !write_quote_attr(quote)) {
        std::cerr << "cannot write regression fixtures\n";
        cleanup_outputs();
        return 1;
    }

    if (test_gtf_parent_synthesis(exe, gtf_basic) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_no_id_collision(exe, gtf_unsorted) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_output(exe, gtf_basic, gtf_mrna) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_multi_parent(exe, multi) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_escaped_comma(exe, esc) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_coord_validation(exe, bad) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_where_url_decode(exe, url) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_summary_scope(exe, gtf_basic) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_error_handling(exe, gtf_basic) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_attr_access(exe, gtf_basic) != 0) {
        cleanup_outputs();
        return 1;
    }

    cleanup_outputs();
    std::cout << "regression_smoke OK\n";
    return 0;
}
