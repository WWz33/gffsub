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

// One transcript whose CDS carries two distinct IDs (the GFF3 spec's
// frameshift example) against a second isoform with a single longer CDS.
static bool write_cds_variants(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t1000\t.\t+\t.\tID=gene01\n"
        << "chr1\tsrc\tmRNA\t100\t1000\t.\t+\t.\tID=tx01;Parent=gene01\n"
        << "chr1\tsrc\tCDS\t100\t250\t.\t+\t0\tID=cds01;Parent=tx01\n"
        << "chr1\tsrc\tCDS\t500\t750\t.\t+\t2\tID=cds02;Parent=tx01\n"
        << "chr1\tsrc\tmRNA\t100\t800\t.\t+\t.\tID=tx02;Parent=gene01\n"
        << "chr1\tsrc\tCDS\t100\t400\t.\t+\t0\tID=cds03;Parent=tx02\n";
    return true;
}

// Flat GTF: transcript/exon only, no gene rows.
static bool write_flat_gtf_isoforms(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "chr1\tsrc\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TA\";\n"
        << "chr1\tsrc\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TA\";\n"
        << "chr1\tsrc\texon\t200\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TA\";\n"
        << "chr1\tsrc\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TB\";\n"
        << "chr1\tsrc\texon\t1\t50\t.\t+\t.\tgene_id \"G1\"; transcript_id \"TB\";\n";
    return true;
}

// GFF3 with a space after every ';' separator.
static bool write_spaced_separators(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n"
        << "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1; Parent=g1\n"
        << "chr1\tsrc\texon\t1\t50\t.\t+\t.\tID=e1; Parent=t1\n";
    return true;
}

// Multi-value Parent: the record belongs to every listed parent.
static bool write_parent_list(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=t;Parent=g1,g2\n"
        << "chr1\tsrc\tmRNA\t1\t10\t.\t+\t.\tID=u;Parent=g3\n";
    return true;
}

// Smallest useful GFF3 line: 20 bytes, below the 22-byte libc++ SSO limit.
static bool write_tiny_gff3(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "a\tb\tc\t1\t1\t.\t.\t.\tID=x\n";
    return true;
}

// GTF with `; Parent=...` inside a quoted value: the naive `;`/`=` scan must
// not read it as a real Parent key.
static bool write_gtf_quoted_parent(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "chr1\tsrc\texon\t100\t200\t.\t+\t.\t"
           "gene_id \"G1\"; transcript_id \"T1\"; note \"see; Parent=bad\";\n";
    return true;
}

// Transcript with two huge CDS segments whose naive sum would wrap negative
// and lose against a 10 bp rival.
static bool write_cds_overflow(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\ts\tgene\t1\t9223372036854775807\t.\t+\t.\tID=g\n"
        << "chr1\ts\tmRNA\t1\t9223372036854775807\t.\t+\t.\tID=big;Parent=g\n"
        << "chr1\ts\tCDS\t1\t4611686018427387905\t.\t+\t0\tID=c1;Parent=big\n"
        << "chr1\ts\tCDS\t1\t4611686018427387905\t.\t+\t0\tID=c1;Parent=big\n"
        << "chr1\ts\tmRNA\t1\t100\t.\t+\t.\tID=small;Parent=g\n"
        << "chr1\ts\tCDS\t1\t10\t.\t+\t0\tID=cs;Parent=small\n";
    return true;
}

// One discontinuous CDS (same ID, two lines) against a shorter rival.
static bool write_cds_discontinuous(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr1\ts\tgene\t1\t900\t.\t+\t.\tID=g\n"
        << "chr1\ts\tmRNA\t1\t900\t.\t+\t.\tID=t1;Parent=g\n"
        << "chr1\ts\tCDS\t1\t100\t.\t+\t0\tID=c;Parent=t1\n"
        << "chr1\ts\tCDS\t500\t700\t.\t+\t0\tID=c;Parent=t1\n"
        << "chr1\ts\tmRNA\t1\t900\t.\t+\t.\tID=t2;Parent=g\n"
        << "chr1\ts\tCDS\t1\t150\t.\t+\t0\tID=d;Parent=t2\n";
    return true;
}

// One transcript ID used by two different transcripts. The spec requires
// unique IDs; the tool must warn instead of silently mis-assigning gene_id.
static bool write_duplicate_transcript_id(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "##gff-version 3\n"
        << "chr2\tsrc\tgene\t1\t100\t.\t+\t.\tID=G2\n"
        << "chr2\tsrc\tmRNA\t1\t80\t.\t+\t.\tID=T;Parent=G2\n"
        << "chr2\tsrc\texon\t1\t80\t.\t+\t.\tID=E2;Parent=T\n"
        << "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=G1\n"
        << "chr1\tsrc\tmRNA\t1\t40\t.\t+\t.\tID=T;Parent=G1\n"
        << "chr1\tsrc\texon\t1\t40\t.\t+\t.\tID=E1;Parent=T\n";
    return true;
}

// GTF col9 carrying both quoted GTF keys and bare GFF3-style keys.
static bool write_gtf_mixed_keys(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) return false;
    out << "chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tgene_id \"G1\"; ID=g1;\n"
        << "chr1\tsrc\ttranscript\t1\t1000\t.\t+\t.\t"
           "gene_id \"G1\"; transcript_id \"T1\"; ID=t1;Parent=g1;\n"
        << "chr1\tsrc\texon\t100\t200\t.\t+\t.\t"
           "gene_id \"G1\"; transcript_id \"T1\"; ID=e1;Parent=t1,t2;\n";
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
        "regression_cds_variants.gff3", "regression_flat.gtf",
        "regression_spaced.gff3", "regression_parent_list.gff3",
        "regression_tiny.gff3", "regression_tiny.gff3.gz",
        "regression_gtf_quoted.gtf", "regression_cds_overflow.gff3",
        "regression_gtf_mixed.gtf",
        "reg_gtf_children.gff3", "reg_gtf_model.gff3",
        "reg_gtf_unsorted_id.gff3", "reg_gtf_out.gtf",
        "reg_gtf3_rename.gff3",
        "reg_multi_t1_children.gff3", "reg_multi_t2_children.gff3",
        "reg_multi_longest.gff3",
        "reg_esc_match.gff3", "reg_esc_nomatch.gff3",
        "reg_bad_coords.gff3",
        "reg_url_match.gff3",
        "reg_gtf_grep.gff3",
        "reg_cds_sum.gff3", "reg_flat_longest.gff3",
        "reg_spaced.gtf", "reg_parent_eq.gff3", "reg_parent_ne.gff3",
        "reg_parent_grep.gff3", "reg_tiny_plain.gff3", "reg_tiny_gz.gff3",
        "reg_max_down.gff3", "reg_bad_sep.err", "reg_bad_region.err",
        "reg_query_empty.err",
        "reg_gtf_quoted.gff3", "reg_gtf_quoted_expr.gff3",
        "reg_gtf_quoted_grep.gff3", "reg_sat_sum.gff3",
        "reg_mixed_id.gff3", "reg_mixed_children.gff3", "reg_mixed_expr.gff3",
        "reg_cds_disc.gff3", "regression_cds_disc.gff3",
        "regression_dup_id.gff3", "reg_dup.gtf", "reg_dup_gtf.err",
        "reg_dup_id.gff3", "reg_dup_index.err", "reg_dup_disc.gff3", "reg_dup_disc.err",
        "reg_gtf_summary.tsv", "reg_gene_summary.tsv",
        "reg_json.json", "reg_json_quote.json",
        "reg_err_missing.err", "reg_err_up.err", "reg_err_threads.err",
        "reg_err_region.err", "reg_err_query_no_selector.err",
        "reg_gtf22.gtf",
        "reg_spec_quote.gff3", "reg_spec_quote.out",
        "reg_spec_comment.gtf", "reg_spec_comment.out",
        "reg_spec_escape.gtf", "reg_spec_escape.out",
        "reg_spec_comma.gff3", "reg_spec_comma.out",
        "reg_spec_version.gff3", "reg_spec_version.out",
        "reg_spec_phase.gff3", "reg_spec_phase.out", "reg_spec_phase.err",
        "reg_spec_multiparent.gff3", "reg_spec_multiparent.out",
        "reg_spec_cycle.gff3", "reg_spec_cycle.out", "reg_spec_cycle.err",
        "reg_same.gff3", "reg_same.out", "reg_same.err",
        "reg_dir.out", "reg_dir.err",
        "reg_bed_in.bed", "reg_bed_gff3.out", "reg_bed_gff3.err",
        "reg_bed_gtf.out", "reg_bed_gtf.err", "reg_bed_out.bed",
        "reg_bed_sum.tsv", "reg_bed_q.out", "reg_bed_q.err",
        "reg_bed_w.out", "reg_bed_w.err",
        "reg_gtf_multiparent.gtf", "reg_mp_where.gff3",
        "reg_mp_drop.gff3", "reg_mp_out.gtf", "reg_mp_out2.gtf",
        "reg_mp_count1.txt", "reg_mp_count2.txt",
        "reg_idparent.gff3", "reg_idparent.gtf", "reg_idparent.bed",
        "reg_def_gene_link.gtf", "reg_def_gene_link.out", "reg_def_gene_link.gtf.out",
        "reg_def_child_tx.gtf", "reg_def_child_tx.gtf.out",
        "reg_def_stop.gff3", "reg_def_stop.out",
        "reg_def_noid.gff3", "reg_def_noid_p.out", "reg_def_noid_m.out",
        "reg_def_empty.bed", "reg_def_empty.out", "reg_def_empty.err",
        "reg_def_empty_q.out", "reg_def_empty_q.err",
        "reg_def_nostdin.out", "reg_def_nostdin.err",
        "reg_def_extra.gff3", "reg_def_extra.out", "reg_def_extra.err",
        "reg_def_strand.gff3", "reg_def_strand.out", "reg_def_strand.err",
        "reg_def_orphan.gff3", "reg_def_orphan.out", "reg_def_orphan.err",
        "reg_def_big.gff3", "reg_def_big.tsv", "reg_def_big.err",
        "reg_def_types.gff3", "reg_def_types.out",
        "reg_def_operand.out", "reg_def_operand.err",
        "reg_def_attach_q.out", "reg_def_attach_w.out",
        "reg_def_threads.out", "reg_def_threads.err"};
    std::system("rmdir reg_dir_input 2>/dev/null");
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
// both gene_id and transcript_id. --format gtf is the GTF3 vocabulary: mRNA
// is renamed to transcript. --format gtf22 drops gene/transcript rows.
static int test_gtf_output(const std::string& exe, const std::string& gtf_basic,
                           const std::string& gtf_mrna) {
    if (run_command(exe + " " + gtf_basic + " --name G1 --model --output-format gtf > reg_gtf_out.gtf") != 0 ||
        require_contains("reg_gtf_out.gtf", "##gtf-version 3") != 0 ||
        // Gene line: gene_id present, no transcript_id on the same line.
        require_contains("reg_gtf_out.gtf", "gene\t100\t500\t.\t+\t.\tgene_id \"G1\";") != 0 ||
        require_not_contains("reg_gtf_out.gtf", "gene\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id") != 0 ||
        // Transcript + exon lines carry both gene_id and transcript_id.
        require_contains("reg_gtf_out.gtf", "transcript\t100\t500\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";") != 0 ||
        require_contains("reg_gtf_out.gtf", "exon\t120\t180\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf_mrna + " --id T1 -C --output-format gtf3 > reg_gtf3_rename.gff3") != 0 ||
        require_contains("reg_gtf3_rename.gff3", "##gtf-version 3") != 0 ||
        require_contains("reg_gtf3_rename.gff3", "transcript\t100\t500") != 0 ||
        require_not_contains("reg_gtf3_rename.gff3", "mRNA") != 0) {
        return 1;
    }
    // GTF2.2 has no gene/transcript feature: only exon-level rows survive.
    if (run_command(exe + " " + gtf_mrna + " --id T1 -C --output-format gtf22 > reg_gtf22.gtf") != 0 ||
        require_contains("reg_gtf22.gtf", "##gtf-version 2.2") != 0 ||
        require_not_contains("reg_gtf22.gtf", "\tgene\t") != 0 ||
        require_not_contains("reg_gtf22.gtf", "\ttranscript\t") != 0 ||
        require_not_contains("reg_gtf22.gtf", "\tmRNA\t") != 0 ||
        require_contains("reg_gtf22.gtf", "exon\t") != 0) {
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

// Group 11: --longest scores a transcript by its LONGEST CDS. Distinct CDS
// IDs under one transcript are alternative products (the GFF3 spec's
// alternative start codons), so tx01 = max(151, 251) = 251 and tx02 (301)
// wins. Same-ID lines are one discontinuous CDS and sum (checked by the
// docs example fixtures elsewhere).
static int test_longest_cds_max_variant(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " --longest > reg_cds_sum.gff3") != 0 ||
        require_contains("reg_cds_sum.gff3", "ID=tx02") != 0 ||
        require_not_contains("reg_cds_sum.gff3", "ID=tx01") != 0) {
        return 1;
    }
    return 0;
}

// Group 11b: same-ID CDS lines are one discontinuous CDS and are summed:
// t1 (100+201=301) beats t2 (150).
static int test_longest_cds_discontinuous_sum(const std::string& exe,
                                              const std::string& gff) {
    if (run_command(exe + " " + gff + " --longest > reg_cds_disc.gff3") != 0 ||
        require_contains("reg_cds_disc.gff3", "ID=t1") != 0 ||
        require_not_contains("reg_cds_disc.gff3", "ID=t2") != 0) {
        return 1;
    }
    return 0;
}

// Group 12: --longest works without gene rows: isoforms are grouped by their
// parent attribute (flat GTF groups by gene_id).
static int test_longest_without_gene_rows(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " " + gtf + " --longest > reg_flat_longest.gff3") != 0 ||
        require_contains("reg_flat_longest.gff3", "ID=TA") != 0 ||
        require_not_contains("reg_flat_longest.gff3", "ID=TB") != 0) {
        return 1;
    }
    return 0;
}

// Group 13: `; ` separators. Parent resolves, so GTF output carries the whole
// gene_id/transcript_id chain instead of empty values. --format gtf uses the
// GTF3 vocabulary, so the mRNA row is relabelled "transcript".
static int test_spaced_separators(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " --output-format gtf > reg_spaced.gtf") != 0 ||
        require_contains("reg_spaced.gtf",
                         "transcript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";") != 0 ||
        require_contains("reg_spaced.gtf",
                         "exon\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";") != 0) {
        return 1;
    }
    return 0;
}

// Group 14: multi-value attributes. A predicate matches when any value
// matches; != matches only when no value does.
static int test_multi_value_attributes(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " -I 'Parent == g2' > reg_parent_eq.gff3") != 0 ||
        require_contains("reg_parent_eq.gff3", "ID=t") != 0 ||
        require_not_contains("reg_parent_eq.gff3", "ID=u") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " -I 'Parent != g2' > reg_parent_ne.gff3") != 0 ||
        require_contains("reg_parent_ne.gff3", "ID=u") != 0 ||
        require_not_contains("reg_parent_ne.gff3", "ID=t") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --grep Parent:g2 > reg_parent_grep.gff3") != 0 ||
        require_contains("reg_parent_grep.gff3", "ID=t") != 0 ||
        require_not_contains("reg_parent_grep.gff3", "ID=u") != 0) {
        return 1;
    }
    return 0;
}

// Group 15: window arithmetic saturates and selectors reject empty values.
static int test_cli_edge_cases(const std::string& exe, const std::string& tiny) {
    // A tiny gzip input (~20 bytes) must survive the move from the parser's
    // buffer into the index: a short string relocates on move.
    if (run_command("gzip -c " + tiny + " > regression_tiny.gff3.gz") != 0 ||
        run_command(exe + " " + tiny + " > reg_tiny_plain.gff3") != 0 ||
        run_command(exe + " regression_tiny.gff3.gz > reg_tiny_gz.gff3") != 0 ||
        run_command("cmp -s reg_tiny_plain.gff3 reg_tiny_gz.gff3") != 0 ||
        run_command(exe + " regression_tiny.gff3.gz --id x > reg_tiny_gz.gff3") != 0 ||
        require_contains("reg_tiny_gz.gff3", "ID=x") != 0) {
        return 1;
    }
    // --down past INT64_MAX would overflow the window end and match nothing.
    if (run_command(exe + " " + tiny + " --id x --down 9223372036854775807 > reg_max_down.gff3") != 0 ||
        require_contains("reg_max_down.gff3", "ID=x") != 0) {
        return 1;
    }
    // Empty seqid and empty query selector: rejected, not silently empty.
    if (require_exit_one_with_error(exe + " " + tiny + " -r :1-100 > /dev/null 2> reg_bad_region.err",
                                    "reg_bad_region.err", "invalid region") != 0 ||
        require_exit_one_with_error(exe + " query " + tiny + " -i '' > /dev/null 2> reg_query_empty.err",
                                    "reg_query_empty.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + tiny + " -w foo > /dev/null 2> reg_bad_sep.err",
                                    "reg_bad_sep.err", "--where expects KEY=VALUE") != 0) {
        return 1;
    }
    return 0;
}

// Group 16: `; Key=value` inside a GTF quoted value is data, not an
// attribute. The synthesized hierarchy fields must come from the real
// gene_id/transcript_id keys.
static int test_gtf_quoted_value(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " " + gtf + " > reg_gtf_quoted.gff3") != 0 ||
        require_contains("reg_gtf_quoted.gff3", "Parent=T1") != 0 ||
        require_contains("reg_gtf_quoted.gff3", "note=see%3B Parent%3Dbad") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf + " -I 'Parent == T1' > reg_gtf_quoted_expr.gff3") != 0 ||
        require_contains("reg_gtf_quoted_expr.gff3", "exon\t100\t200") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf + " --grep Parent:T1 > reg_gtf_quoted_grep.gff3") != 0 ||
        require_contains("reg_gtf_quoted_grep.gff3", "exon\t100\t200") != 0) {
        return 1;
    }
    return 0;
}

// Group 17: --longest per-ID CDS sums and exon sums saturate instead of
// wrapping negative. The big transcript's two same-ID segments each measure
// 2^62+1 bp; their saturated sum must still beat the 10 bp rival.
static int test_longest_sum_saturation(const std::string& exe, const std::string& gff) {
    if (run_command(exe + " " + gff + " --longest > reg_sat_sum.gff3") != 0 ||
        require_contains("reg_sat_sum.gff3", "ID=big") != 0 ||
        require_not_contains("reg_sat_sum.gff3", "ID=small") != 0) {
        return 1;
    }
    return 0;
}

// Group 18: GTF col9 with both quoted GTF keys and bare GFF3-style keys.
// The bare keys must resolve through the quote-aware parser: IDs select,
// and a bare multi-value Parent matches per value.
static int test_gtf_mixed_keys(const std::string& exe, const std::string& gtf) {
    if (run_command(exe + " " + gtf + " --id e1 > reg_mixed_id.gff3") != 0 ||
        require_contains("reg_mixed_id.gff3", "exon\t100\t200") != 0 ||
        require_contains("reg_mixed_id.gff3", "ID=e1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf + " --id g1 -C > reg_mixed_children.gff3") != 0 ||
        require_contains("reg_mixed_children.gff3", "ID=t1") != 0 ||
        require_contains("reg_mixed_children.gff3", "ID=e1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gtf + " -I 'Parent == t2' > reg_mixed_expr.gff3") != 0 ||
        require_contains("reg_mixed_expr.gff3", "exon\t100\t200") != 0) {
        return 1;
    }
    return 0;
}

// Group 19: a repeated transcript ID is invalid GFF3 (IDs must be unique).
// Both consumers warn: the index for lineage lookups, the GTF writer for the
// gene_id it hands to children. A discontinuous feature (same ID, matching
// type/seqid/strand) stays silent.
static int test_duplicate_id_warning(const std::string& exe, const std::string& gff,
                                     const std::string& discontinuous) {
    if (run_command(exe + " " + gff + " -f gtf > reg_dup.gtf 2> reg_dup_gtf.err") != 0 ||
        require_contains("reg_dup_gtf.err", "Warning: transcript ID") != 0 ||
        require_contains("reg_dup.gtf", "gene_id") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " -i G2 -C > reg_dup_id.gff3 2> reg_dup_index.err") != 0 ||
        require_contains("reg_dup_index.err", "Warning: ID") != 0) {
        return 1;
    }
    if (run_command(exe + " " + discontinuous + " --longest > reg_dup_disc.gff3 2> reg_dup_disc.err") != 0 ||
        require_not_contains("reg_dup_disc.err", "Warning") != 0) {
        return 1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
// Group 16: GFF3 / GTF conformance per the published specs
// ---------------------------------------------------------------------------

// A quoted value on the first data line must not flip a GFF3 file to GTF
// parsing: GFF3 values are `tag=value`, GTF pairs are `tag "value"`.
static int spec_quoted_gff3(const std::string& exe) {
    const std::string f{"reg_spec_quote.gff3"};
    std::ofstream out{f};
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=\"A\"\n"
        << "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n";
    out.close();
    if (run_command(exe + " " + f + " > reg_spec_quote.out") != 0 ||
        require_contains("reg_spec_quote.out", "ID=g1;Name=\"A\"") != 0) {
        return 1;
    }
    return 0;
}

// GTF2.2: `#` outside quotes starts a comment; it must not become an attribute.
static int spec_gtf_trailing_comment(const std::string& exe) {
    const std::string f{"reg_spec_comment.gtf"};
    std::ofstream out{f};
    out << "chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; # tail\n";
    out.close();
    if (run_command(exe + " " + f + " -f gtf > reg_spec_comment.out") != 0 ||
        require_contains("reg_spec_comment.out", "transcript_id \"t1\";\n") != 0 ||
        require_not_contains("reg_spec_comment.out", "tail") != 0) {
        return 1;
    }
    return 0;
}

// GXF2 C-style escapes: `\t` decodes to a tab, and the writer re-escapes it,
// so the round trip keeps the value.
static int spec_gtf_escape_roundtrip(const std::string& exe) {
    const std::string f{"reg_spec_escape.gtf"};
    std::ofstream out{f};
    out << "chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; note \"a\\tb\";\n";
    out.close();
    if (run_command(exe + " " + f + " -f gtf > reg_spec_escape.out") != 0 ||
        require_contains("reg_spec_escape.out", "note \"a\\tb\";") != 0 ||
        require_not_contains("reg_spec_escape.out", "a\\\\tb") != 0) {
        return 1;
    }
    return 0;
}

// Only Parent, Alias, Note, Dbxref and Ontology_term carry comma lists; a
// comma in any other tag is part of the value.
static int spec_comma_scope(const std::string& exe) {
    const std::string f{"reg_spec_comma.gff3"};
    std::ofstream out{f};
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1;Name=A,B;Note=x,y\n";
    out.close();
    if (run_command(exe + " " + f + " -f gtf > reg_spec_comma.out") != 0 ||
        require_contains("reg_spec_comma.out", "Name \"A,B\";") != 0 ||
        require_contains("reg_spec_comma.out", "Note \"x\"; Note \"y\";") != 0) {
        return 1;
    }
    return 0;
}

// Output always declares GFF3 version 3, whatever the input claimed.
static int spec_version_header(const std::string& exe) {
    const std::string f{"reg_spec_version.gff3"};
    std::ofstream out{f};
    out << "##gff-version 2\nchr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=g1\n";
    out.close();
    if (run_command(exe + " " + f + " > reg_spec_version.out") != 0 ||
        require_contains("reg_spec_version.out", "##gff-version 3") != 0 ||
        require_not_contains("reg_spec_version.out", "##gff-version 2") != 0) {
        return 1;
    }
    return 0;
}

// CDS without phase 0/1/2 is kept but reported.
static int spec_cds_phase_warning(const std::string& exe) {
    const std::string f{"reg_spec_phase.gff3"};
    std::ofstream out{f};
    out << "##gff-version 3\n"
        << "chr1\tsrc\tCDS\t1\t100\t.\t+\t.\tID=c1\n";
    out.close();
    if (run_command(exe + " " + f + " > reg_spec_phase.out 2> reg_spec_phase.err") != 0 ||
        require_contains("reg_spec_phase.err", "without phase") != 0 ||
        require_contains("reg_spec_phase.out", "CDS\t1\t100") != 0) {
        return 1;
    }
    return 0;
}

// GTF has no multi-parent feature: a shared exon becomes one row per parent.
static int spec_multi_parent_gtf(const std::string& exe) {
    const std::string f{"reg_spec_multiparent.gff3"};
    std::ofstream out{f};
    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t1\t900\t.\t+\t.\tID=g1\n"
        << "chr1\tsrc\tmRNA\t1\t400\t.\t+\t.\tID=t1;Parent=g1\n"
        << "chr1\tsrc\tmRNA\t500\t900\t.\t+\t.\tID=t2;Parent=g1\n"
        << "chr1\tsrc\texon\t1\t400\t.\t+\t.\tID=e1;Parent=t1,t2\n";
    out.close();
    if (run_command(exe + " " + f + " -f gtf > reg_spec_multiparent.out") != 0 ||
        require_contains("reg_spec_multiparent.out", "transcript_id \"t1\";") != 0 ||
        require_contains("reg_spec_multiparent.out", "transcript_id \"t2\";") != 0) {
        return 1;
    }
    return 0;
}

// Parent links that form a cycle are invalid; the tool warns instead of
// silently producing a wrong hierarchy.
static int spec_parent_cycle(const std::string& exe) {
    const std::string f{"reg_spec_cycle.gff3"};
    std::ofstream out{f};
    out << "##gff-version 3\n"
        << "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=a;Parent=b\n"
        << "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=b;Parent=a\n";
    out.close();
    if (run_command(exe + " " + f + " -i a -p > reg_spec_cycle.out 2> reg_spec_cycle.err") != 0 ||
        require_contains("reg_spec_cycle.err", "Parent cycle") != 0) {
        return 1;
    }
    return 0;
}

// File-vs-file aliasing and non-file inputs: both used to be accepted and
// produce destroyed or empty output.
static int test_io_guards(const std::string& exe) {
    const std::string same{"reg_same.gff3"};
    {
        std::ofstream out{same};
        if (!out.is_open()) return 1;
        out << "##gff-version 3\n"
            << "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=x\n";
    }
    if (expect_command_failure(exe + " " + same + " -o " + same +
                               " > reg_same.out 2> reg_same.err") != 0 ||
        require_contains("reg_same.err", "same file") != 0 ||
        require_contains(same, "ID=x") != 0) {
        return 1;
    }
    if (run_command("mkdir -p reg_dir_input") != 0 ||
        expect_command_failure(exe + " reg_dir_input > reg_dir.out 2> reg_dir.err") != 0 ||
        require_contains("reg_dir.err", "cannot parse") != 0) {
        return 1;
    }
    run_command("rmdir reg_dir_input");
    return 0;
}

// BED is a BED-output and summary input only; converting it to GFF3/GTF is
// rejected on every entry point.
static int test_bed_conversion_removed(const std::string& exe) {
    {
        std::ofstream bed{"reg_bed_in.bed"};
        if (!bed.is_open()) return 1;
        bed << "chr1\t0\t10\tbedrec\t0\t+\n";
    }
    if (expect_command_failure(exe + " reg_bed_in.bed > reg_bed_gff3.out 2> reg_bed_gff3.err") != 0 ||
        require_contains("reg_bed_gff3.err", "not converted") != 0 ||
        expect_command_failure(exe + " reg_bed_in.bed -f gtf > reg_bed_gtf.out 2> reg_bed_gtf.err") != 0 ||
        require_contains("reg_bed_gtf.err", "not converted") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_bed_in.bed -f bed > reg_bed_out.bed") != 0 ||
        require_contains("reg_bed_out.bed", "chr1\t0\t10\tbedrec") != 0 ||
        run_command(exe + " reg_bed_in.bed -s > reg_bed_sum.tsv") != 0 ||
        require_contains("reg_bed_sum.tsv", "chr1\tregion\t1") != 0) {
        return 1;
    }
    if (expect_command_failure(exe + " query reg_bed_in.bed -i bedrec > reg_bed_q.out 2> reg_bed_q.err") != 0 ||
        require_contains("reg_bed_q.err", "not converted") != 0 ||
        expect_command_failure(exe + " window reg_bed_in.bed -i bedrec > reg_bed_w.out 2> reg_bed_w.err") != 0 ||
        require_contains("reg_bed_w.err", "not converted") != 0) {
        return 1;
    }
    return 0;
}

// Quoted GTF Parent "t1,t2" is a two-parent list in every consumer: attribute
// match, orphan retention, and one GTF row per parent.
static int test_gtf_quoted_multi_parent(const std::string& exe) {
    {
        std::ofstream gtf{"reg_gtf_multiparent.gtf"};
        if (!gtf.is_open()) return 1;
        gtf << "chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n"
            << "chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
            << "chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\";\n"
            << "chr1\tx\texon\t1\t50\t.\t+\t.\tgene_id \"g1\"; Parent \"t1,t2\";\n";
    }
    if (run_command(exe + " reg_gtf_multiparent.gtf -w Parent=t2 > reg_mp_where.gff3") != 0 ||
        require_contains("reg_mp_where.gff3", "exon\t1\t50") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_gtf_multiparent.gtf --drop-orphans > reg_mp_drop.gff3") != 0 ||
        require_contains("reg_mp_drop.gff3", "ID=t1") != 0 ||
        require_contains("reg_mp_drop.gff3", "ID=t2") != 0 ||
        require_contains("reg_mp_drop.gff3", "exon\t1\t50") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_gtf_multiparent.gtf -f gtf > reg_mp_out.gtf") != 0 ||
        require_contains("reg_mp_out.gtf", "transcript_id \"t1\";") != 0 ||
        require_contains("reg_mp_out.gtf", "transcript_id \"t2\";") != 0) {
        return 1;
    }
    // One row per parent, each naming its own transcript. Re-reading the
    // output must not expand again: a repeated Parent list doubled the exon
    // row count on every pass.
    if (run_command("grep -c \"exon\" reg_mp_out.gtf > reg_mp_count1.txt") != 0 ||
        require_contains("reg_mp_count1.txt", "2") != 0 ||
        require_not_contains("reg_mp_out.gtf", "Parent \"t1,t2\"") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_mp_out.gtf -f gtf > reg_mp_out2.gtf") != 0 ||
        run_command("grep -c \"exon\" reg_mp_out2.gtf > reg_mp_count2.txt") != 0 ||
        require_contains("reg_mp_count2.txt", "2") != 0) {
        return 1;
    }
    return 0;
}

// GTF output keeps ID/Parent from GFF3 input (AGAT does the same). BED output
// writes an integer score in 0..1000: non-finite, fractional, negative, and
// out-of-range source scores all fall back to 0 or the range edge.
static int test_gtf_attrs_and_bed_score(const std::string& exe) {
    {
        std::ofstream gff{"reg_idparent.gff3"};
        if (!gff.is_open()) return 1;
        gff << "##gff-version 3\n"
            << "chr1\tx\tgene\t1\t100\tinf\t+\t.\tID=g1\n"
            << "chr1\tx\tmRNA\t1\t100\t0.5\t+\t.\tID=t1;Parent=g1\n"
            << "chr1\tx\texon\t1\t100\t-2\t+\t.\tID=e1;Parent=t1\n"
            << "chr1\tx\texon\t1\t100\t1001\t+\t.\tID=e2;Parent=t1\n"
            << "chr1\tx\texon\t1\t100\tabc\t+\t.\tID=e3;Parent=t1\n";
    }
    if (run_command(exe + " reg_idparent.gff3 -f gtf > reg_idparent.gtf") != 0 ||
        require_contains("reg_idparent.gtf", "transcript_id \"t1\"; Parent \"g1\"; ID \"t1\";") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_idparent.gff3 -f bed > reg_idparent.bed") != 0 ||
        require_contains("reg_idparent.bed", "chr1\t0\t100\tg1\t0\t+") != 0 ||
        require_contains("reg_idparent.bed", "chr1\t0\t100\tt1\t0\t+") != 0 ||
        require_contains("reg_idparent.bed", "chr1\t0\t100\te1\t0\t+") != 0 ||
        require_contains("reg_idparent.bed", "chr1\t0\t100\te2\t1000\t+") != 0 ||
        require_contains("reg_idparent.bed", "chr1\t0\t100\te3\t0\t+") != 0) {
        return 1;
    }
    return 0;
}

static int test_spec_conformance(const std::string& exe) {
    return spec_quoted_gff3(exe) || spec_gtf_trailing_comment(exe) ||
           spec_gtf_escape_roundtrip(exe) || spec_comma_scope(exe) ||
           spec_version_header(exe) || spec_cds_phase_warning(exe) ||
           spec_multi_parent_gtf(exe) || spec_parent_cycle(exe);
}

// --- defects recorded in docs/defect-review.md ----------------------------

// A GTF transcript carrying a non-standard Parent attribute still links to
// its gene: --longest must keep the gene and compete the transcripts under it.
// A converter that writes `Parent "gene:g1"` on transcript rows produces this
// shape, so the attribute must not displace the gene_id link.
static int test_gtf_transcript_parent_keeps_gene(const std::string& exe) {
    {
        std::ofstream gtf{"reg_def_gene_link.gtf"};
        if (!gtf.is_open()) return 1;
        gtf << "chr1\tx\tgene\t1\t100\t.\t+\t.\tgene_id \"g1\";\n"
            << "chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\"; Parent \"g1\";\n"
            << "chr1\tx\ttranscript\t1\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\"; Parent \"g1\";\n"
            << "chr1\tx\texon\t1\t50\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
            << "chr1\tx\texon\t1\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t2\";\n";
    }
    if (run_command(exe + " reg_def_gene_link.gtf --longest > reg_def_gene_link.out") != 0) {
        return 1;
    }
    // g1 kept, the longer transcript survives, the shorter one does not.
    if (require_contains("reg_def_gene_link.out", "ID=g1") != 0 ||
        require_contains("reg_def_gene_link.out", "ID=t2;Parent=g1") != 0 ||
        require_not_contains("reg_def_gene_link.out", "ID=t1") != 0) {
        return 1;
    }
    if (run_command(exe + " reg_def_gene_link.gtf -f gtf > reg_def_gene_link.gtf.out") != 0 ||
        require_contains("reg_def_gene_link.gtf.out", "transcript_id \"t1\"; Parent \"g1\";") != 0 ||
        require_contains("reg_def_gene_link.gtf.out", "transcript_id \"t2\"; Parent \"g1\";") != 0) {
        return 1;
    }
    return 0;
}

// A child row that declares its own transcript_id and a Parent list names a
// distinct transcript on each expanded row.
static int test_gtf_child_transcript_id_follows_parent(const std::string& exe) {
    {
        std::ofstream gtf{"reg_def_child_tx.gtf"};
        if (!gtf.is_open()) return 1;
        gtf << "chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"txA\";\n"
            << "chr1\tx\ttranscript\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"txB\";\n"
            << "chr1\tx\texon\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"tx1\"; Parent \"txA,txB\";\n";
    }
    if (run_command(exe + " reg_def_child_tx.gtf -f gtf > reg_def_child_tx.gtf.out") != 0) {
        return 1;
    }
    if (require_contains("reg_def_child_tx.gtf.out", "transcript_id \"txA\";") != 0 ||
        require_contains("reg_def_child_tx.gtf.out", "transcript_id \"txB\";") != 0 ||
        require_not_contains("reg_def_child_tx.gtf.out", "transcript_id \"tx1\";") != 0) {
        return 1;
    }
    return 0;
}

// A stop_codon line shared by two transcripts trims the terminal CDS of both.
static int test_shared_stop_codon(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_stop.gff3"};
        if (!gff.is_open()) return 1;
        gff << "##gff-version 3\n"
            << "chr1\ts\tgene\t1\t200\t.\t+\t.\tID=g\n"
            << "chr1\ts\tmRNA\t1\t103\t.\t+\t.\tID=t1;Parent=g\n"
            << "chr1\ts\tmRNA\t1\t103\t.\t+\t.\tID=t2;Parent=g\n"
            << "chr1\ts\tCDS\t1\t103\t.\t+\t0\tID=c1;Parent=t1\n"
            << "chr1\ts\tCDS\t1\t103\t.\t+\t0\tID=c2;Parent=t2\n"
            << "chr1\ts\tstop_codon\t101\t103\t.\t+\t0\tParent=t1,t2\n";
    }
    if (run_command(exe + " reg_def_stop.gff3 -f gtf > reg_def_stop.out") != 0) {
        return 1;
    }
    if (require_contains("reg_def_stop.out", "CDS\t1\t100\t.\t+\t0\tgene_id \"g\"; transcript_id \"t1\"") != 0 ||
        require_contains("reg_def_stop.out", "CDS\t1\t100\t.\t+\t0\tgene_id \"g\"; transcript_id \"t2\"") != 0) {
        return 1;
    }
    return 0;
}

// An exon without an ID still reaches its ancestors and its gene model because
// its Parent links supply the anchor.
static int test_parent_expansion_without_id(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_noid.gff3"};
        if (!gff.is_open()) return 1;
        gff << "##gff-version 3\n"
            << "chr1\ts\tgene\t1\t900\t.\t+\t.\tID=g\n"
            << "chr1\ts\tmRNA\t1\t400\t.\t+\t.\tID=t1;Parent=g\n"
            << "chr1\ts\tmRNA\t500\t900\t.\t+\t.\tID=t2;Parent=g\n"
            << "chr1\ts\texon\t1\t100\t.\t+\t.\tParent=t1,t2\n";
    }
    if (run_command(exe + " query reg_def_noid.gff3 -w Parent=t2 -p > reg_def_noid_p.out") != 0 ||
        require_contains("reg_def_noid_p.out", "ID=g") != 0 ||
        require_contains("reg_def_noid_p.out", "ID=t1") != 0 ||
        require_contains("reg_def_noid_p.out", "ID=t2") != 0) {
        return 1;
    }
    if (run_command(exe + " query reg_def_noid.gff3 -w Parent=t2 -m > reg_def_noid_m.out") != 0 ||
        require_contains("reg_def_noid_m.out", "ID=g") != 0 ||
        require_contains("reg_def_noid_m.out", "ID=t1") != 0 ||
        require_contains("reg_def_noid_m.out", "ID=t2") != 0 ||
        require_contains("reg_def_noid_m.out", "Parent=t1,t2") != 0) {
        return 1;
    }
    return 0;
}

// An empty BED file is still BED: the conversion guard fires on the format the
// input was read as, not on the presence of records.
static int test_empty_bed_rejected(const std::string& exe) {
    { std::ofstream empty{"reg_def_empty.bed"}; }
    if (expect_command_failure(exe + " reg_def_empty.bed > reg_def_empty.out 2> reg_def_empty.err") != 0 ||
        require_contains("reg_def_empty.err", "not converted") != 0) {
        return 1;
    }
    if (expect_command_failure(exe + " query reg_def_empty.bed -i x > reg_def_empty_q.out 2> reg_def_empty_q.err") != 0 ||
        require_contains("reg_def_empty_q.err", "not converted") != 0) {
        return 1;
    }
    return 0;
}

// A closed stdin is an input error, not an empty annotation.
static int test_closed_stdin_rejected(const std::string& exe) {
    if (expect_command_failure(exe + " - <&- > reg_def_nostdin.out 2> reg_def_nostdin.err") != 0 ||
        require_contains("reg_def_nostdin.err", "cannot parse stdin") != 0) {
        return 1;
    }
    return 0;
}

// A row with more than nine columns is reported and skipped, not silently
// truncated to nine.
static int test_extra_columns_rejected(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_extra.gff3"};
        if (!gff.is_open()) return 1;
        gff << "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\tjunk\n";
    }
    if (run_command(exe + " reg_def_extra.gff3 > reg_def_extra.out 2> reg_def_extra.err") != 0 ||
        require_contains("reg_def_extra.err", "too many columns") != 0 ||
        require_not_contains("reg_def_extra.out", "ID=g1") != 0) {
        return 1;
    }
    return 0;
}

// A multi-character strand or phase is invalid; the row is skipped, and the
// value is not truncated to one character.
static int test_bad_strand_phase_skipped(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_strand.gff3"};
        if (!gff.is_open()) return 1;
        gff << "chr1\ts\tgene\t1\t9\t.\t++\t.\tID=g1\n";
    }
    if (run_command(exe + " reg_def_strand.gff3 > reg_def_strand.out 2> reg_def_strand.err") != 0 ||
        require_contains("reg_def_strand.err", "invalid strand") != 0 ||
        require_not_contains("reg_def_strand.out", "ID=g1") != 0) {
        return 1;
    }
    return 0;
}

// A Parent naming no record in the file is reported once.
static int test_missing_parent_warned(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_orphan.gff3"};
        if (!gff.is_open()) return 1;
        gff << "##gff-version 3\n"
            << "chr1\ts\texon\t10\t20\t.\t+\t.\tParent=missing\n";
    }
    if (run_command(exe + " query reg_def_orphan.gff3 -w Parent=missing -p > reg_def_orphan.out 2> reg_def_orphan.err") != 0 ||
        require_contains("reg_def_orphan.err", "has no record in the input") != 0) {
        return 1;
    }
    return 0;
}

// Summary output saturates on coordinates near the int64 limit instead of
// printing negative lengths.
static int test_summary_no_overflow(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_big.gff3"};
        if (!gff.is_open()) return 1;
        gff << "chr1\tt\texon\t1\t9223372036854775807\t.\t+\t.\tID=r1\n"
            << "chr1\tt\texon\t9223372036854775806\t9223372036854775807\t.\t+\t.\tID=r2\n";
    }
    if (run_command(exe + " reg_def_big.gff3 -s > reg_def_big.tsv 2> reg_def_big.err") != 0 ||
        require_not_contains("reg_def_big.tsv", "-922337203685477") != 0 ||
        require_contains("reg_def_big.tsv", "9223372036854775807") != 0) {
        return 1;
    }
    return 0;
}

// A comma list is split after trimming, so `-t "mRNA, exon"` keeps both types.
static int test_type_list_trims_spaces(const std::string& exe) {
    {
        std::ofstream gff{"reg_def_types.gff3"};
        if (!gff.is_open()) return 1;
        gff << "##gff-version 3\n"
            << "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g\n"
            << "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\n"
            << "chr1\ts\texon\t1\t9\t.\t+\t.\tParent=t\n";
    }
    if (run_command(exe + " reg_def_types.gff3 -t \"mRNA, exon\" > reg_def_types.out") != 0 ||
        require_contains("reg_def_types.out", "mRNA") != 0 ||
        require_contains("reg_def_types.out", "exon") != 0) {
        return 1;
    }
    return 0;
}

// A second positional operand is a usage error, not a silently ignored file.
static int test_extra_operand_rejected(const std::string& exe) {
    if (expect_command_failure(exe + " reg_def_types.gff3 typo.gff3 > reg_def_operand.out 2> reg_def_operand.err") != 0 ||
        require_contains("reg_def_operand.err", "unexpected argument") != 0) {
        return 1;
    }
    return 0;
}

// Subcommands accept --flag=value the same way the default CLI does.
static int test_subcommand_attached_value(const std::string& exe) {
    if (run_command(exe + " query reg_def_types.gff3 --id=g > reg_def_attach_q.out") != 0 ||
        require_contains("reg_def_attach_q.out", "ID=g") != 0) {
        return 1;
    }
    if (run_command(exe + " window reg_def_types.gff3 -i t --up=5 > reg_def_attach_w.out") != 0 ||
        require_contains("reg_def_attach_w.out", "ID=t") != 0) {
        return 1;
    }
    return 0;
}

// A whitespace-prefixed negative thread count is rejected like an unspaced one.
static int test_threads_rejects_space_negative(const std::string& exe) {
    if (expect_command_failure(exe + " reg_def_types.gff3 --longest --threads ' -1' > reg_def_threads.out 2> reg_def_threads.err") != 0 ||
        require_contains("reg_def_threads.err", "non-negative") != 0) {
        return 1;
    }
    return 0;
}

static int test_defect_regressions(const std::string& exe) {
    return test_gtf_transcript_parent_keeps_gene(exe) ||
           test_gtf_child_transcript_id_follows_parent(exe) ||
           test_shared_stop_codon(exe) ||
           test_parent_expansion_without_id(exe) ||
           test_empty_bed_rejected(exe) ||
           test_closed_stdin_rejected(exe) ||
           test_extra_columns_rejected(exe) ||
           test_bad_strand_phase_skipped(exe) ||
           test_missing_parent_warned(exe) ||
           test_summary_no_overflow(exe) ||
           test_type_list_trims_spaces(exe) ||
           test_extra_operand_rejected(exe) ||
           test_subcommand_attached_value(exe) ||
           test_threads_rejects_space_negative(exe);
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
    const std::string cds_variants{"regression_cds_variants.gff3"};
    const std::string flat_gtf{"regression_flat.gtf"};
    const std::string spaced{"regression_spaced.gff3"};
    const std::string parent_list{"regression_parent_list.gff3"};
    const std::string tiny{"regression_tiny.gff3"};
    const std::string gtf_quoted{"regression_gtf_quoted.gtf"};
    const std::string cds_overflow{"regression_cds_overflow.gff3"};
    const std::string cds_disc{"regression_cds_disc.gff3"};
    const std::string dup_id{"regression_dup_id.gff3"};
    const std::string gtf_mixed{"regression_gtf_mixed.gtf"};

    if (!write_gtf_basic(gtf_basic) || !write_gtf_unsorted(gtf_unsorted) ||
        !write_gtf_mrna(gtf_mrna) || !write_multi_parent(multi) ||
        !write_escaped_comma(esc) || !write_bad_coords(bad) ||
        !write_url_encoded(url) || !write_quote_attr(quote) ||
        !write_cds_variants(cds_variants) || !write_flat_gtf_isoforms(flat_gtf) ||
        !write_spaced_separators(spaced) || !write_parent_list(parent_list) ||
        !write_tiny_gff3(tiny) || !write_gtf_quoted_parent(gtf_quoted) ||
        !write_cds_overflow(cds_overflow) || !write_gtf_mixed_keys(gtf_mixed) ||
        !write_cds_discontinuous(cds_disc) ||
        !write_duplicate_transcript_id(dup_id)) {
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
    if (test_longest_cds_max_variant(exe, cds_variants) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_longest_cds_discontinuous_sum(exe, cds_disc) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_longest_without_gene_rows(exe, flat_gtf) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_spaced_separators(exe, spaced) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_multi_value_attributes(exe, parent_list) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_cli_edge_cases(exe, tiny) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_io_guards(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_bed_conversion_removed(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_quoted_multi_parent(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_attrs_and_bed_score(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_quoted_value(exe, gtf_quoted) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_longest_sum_saturation(exe, cds_overflow) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_gtf_mixed_keys(exe, gtf_mixed) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_duplicate_id_warning(exe, dup_id, cds_disc) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_spec_conformance(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    if (test_defect_regressions(exe) != 0) {
        cleanup_outputs();
        return 1;
    }
    cleanup_outputs();
    std::cout << "regression_smoke OK\n";
    return 0;
}
