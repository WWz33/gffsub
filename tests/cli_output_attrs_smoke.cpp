#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "test_utils.hpp"

using test_utils::expect_command_failure;
using test_utils::read_file;
using test_utils::require_contains;
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

    // multi seqid: all row per type with summed counts
    if (require_contains("summary_gene.tsv", "all\tgene\t2\t702\t301\t351\t401") != 0) {
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
