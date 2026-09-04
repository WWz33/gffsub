#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "test_utils.hpp"

using test_utils::expect_command_failure;
using test_utils::read_file;
using test_utils::require_contains;
using test_utils::require_not_contains;
using test_utils::run_command;

static bool write_test_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1\n"
        << "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tID=exon1;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t150\t170\t.\t+\t0\tID=cds1;Parent=tx1\n"
        << "chr2\tsrc\tgene\t200\t600\t.\t-\t.\tID=gene0002;Name=XYZ1\n";
    return true;
}

static void cleanup_outputs() {
    std::remove("cli_output_attrs_smoke.gff3");
    std::remove("summary_gene.tsv");
    std::remove("summary_exon.tsv");
    std::remove("summary_query.tsv");
    std::remove("summary_help.txt");
    std::remove("summary_bad.out");
    std::remove("summary_bad.err");
    std::remove("summary_gene_out.tsv");
    std::remove("summary_overlap.gff3");
    std::remove("summary_overlap.tsv");
    std::remove("summary_tbed.tsv");
    std::remove("summary_gtf_attrs.gff3");
    std::remove("summary_gtf.gtf");
    std::remove("summary_borrowed.gff3");
    std::remove("summary_borrowed_out.gff3");
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cli_output_attrs_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string gff{"cli_output_attrs_smoke.gff3"};
    if (!write_test_annotation(gff)) {
        std::cerr << "cannot write test annotation\n";
        return 1;
    }

    // --help shows -s
    if (run_command(exe + " --help > summary_help.txt 2>&1") != 0 ||
        require_contains("summary_help.txt", "-s, --summary") != 0) {
        return 1;
    }

    // -t gene -s: aggregated TSV by seqid × type
    if (run_command(exe + " " + gff + " -t gene -s > summary_gene.tsv") != 0 ||
        require_contains("summary_gene.tsv", "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len") != 0 ||
        require_contains("summary_gene.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0 ||
        require_contains("summary_gene.tsv", "chr2\tgene\t1\t401\t401\t401\t401") != 0) {
        return 1;
    }

    // multi seqid: all row per type with full column check
    // (Q1=301 Q2=351 Q3=401; coverage = 301 + 401, intervals on distinct seqids)
    if (require_contains("summary_gene.tsv", "all\tgene\t2\t702\t301\t351\t401\t301\t351\t401\t702") != 0) {
        return 1;
    }

    // overlapping exons on two seqids: coverage dedups overlaps within a
    // seqid; all-row coverage sums the per-seqid unions.
    {
        std::ofstream ov{"summary_overlap.gff3"};
        if (!ov.is_open()) return 1;
        ov << "##gff-version 3\n"
           << "chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=e1\n"
           << "chr1\tsrc\texon\t150\t250\t.\t+\t.\tID=e2\n"
           << "chr2\tsrc\texon\t500\t600\t.\t+\t.\tID=e3\n";
    }
    // chr1: [100,250] union = 151; chr2: 101; all: 252
    // lens 101,101,101 → Q1=Q2=Q3=101
    if (run_command(exe + " summary_overlap.gff3 -t exon -s > summary_overlap.tsv") != 0 ||
        require_contains("summary_overlap.tsv", "chr1\texon\t2\t202\t101\t101\t101\t101\t101\t101\t151") != 0 ||
        require_contains("summary_overlap.tsv", "chr2\texon\t1\t101\t101\t101\t101\t101\t101\t101\t101") != 0 ||
        require_contains("summary_overlap.tsv", "all\texon\t3\t303\t101\t101\t101\t101\t101\t101\t252") != 0) {
        return 1;
    }

    // -t bed is a type filter, not an output format: no record has type
    // "bed", so output is empty (header only for -s)
    if (run_command(exe + " " + gff + " -t bed -s > summary_tbed.tsv") != 0 ||
        require_contains("summary_tbed.tsv", "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len") != 0 ||
        require_not_contains("summary_tbed.tsv", "chr1") != 0) {
        return 1;
    }

    // GTF output preserves non-gene_id/transcript_id attributes
    // (AGAT "attribute conserved: All"), URL-decoded and quoted
    {
        std::ofstream ga{"summary_gtf_attrs.gff3"};
        if (!ga.is_open()) return 1;
        ga << "##gff-version 3\n"
           << "chr1\tsrc\tgene\t100\t300\t.\t+\t.\tID=g1;Name=MyGene\n"
           << "chr1\tsrc\tmRNA\t100\t300\t.\t+\t.\tID=t1;Parent=g1;biotype=protein_coding\n"
           << "chr1\tsrc\texon\t100\t200\t.\t+\t.\tParent=t1;Note=a%2Cb\n";
    }
    if (run_command(exe + " summary_gtf_attrs.gff3 --format gtf -o summary_gtf.gtf") != 0 ||
        require_contains("summary_gtf.gtf", "gene_id \"g1\"; Name \"MyGene\";") != 0 ||
        require_contains("summary_gtf.gtf", "gene_id \"g1\"; transcript_id \"t1\"; biotype \"protein_coding\";") != 0 ||
        require_contains("summary_gtf.gtf", "Note \"a,b\";") != 0) {
        return 1;
    }

    // borrowed index (from_data const&): nearest + seqid subset combination
    // (regression guard for recs_ view vs data_ copy)
    {
        std::ofstream bi{"summary_borrowed.gff3"};
        if (!bi.is_open()) return 1;
        bi << "##gff-version 3\n"
           << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=ga\n"
           << "chr1\tsrc\texon\t100\t150\t.\t+\t.\tID=ea;Parent=ga\n"
           << "chr2\tsrc\tgene\t400\t500\t.\t+\t.\tID=gb\n";
    }
    if (run_command(exe + " summary_borrowed.gff3 --nearest chr1:300-320 --seqid chr1 > summary_borrowed_out.gff3") != 0 ||
        require_contains("summary_borrowed_out.gff3", "ID=ga") != 0 ||
        require_not_contains("summary_borrowed_out.gff3", "ID=gb") != 0) {
        return 1;
    }

    // -t exon -s: no all row when single seqid
    if (run_command(exe + " " + gff + " -t exon -s > summary_exon.tsv") != 0 ||
        require_contains("summary_exon.tsv", "chr1\texon\t1\t61\t61\t61\t61") != 0) {
        return 1;
    }
    // single seqid → no all row
    const auto exon_text = read_file("summary_exon.tsv");
    if (exon_text.find("all\t") != std::string::npos) {
        std::cerr << "unexpected all row in single-seqid summary\n";
        return 1;
    }

    // query -i gene0001 -s
    if (run_command(exe + " query " + gff + " -i gene0001 -s > summary_query.tsv") != 0 ||
        require_contains("summary_query.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }

    // -s with -o writes to file
    if (run_command(exe + " " + gff + " -t gene -s -o summary_gene_out.tsv") != 0 ||
        require_contains("summary_gene_out.tsv", "chr1\tgene\t1\t301") != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_output_attrs_smoke OK\n";
    return 0;
}
