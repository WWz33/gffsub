#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

static bool write_test_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1;gene_id=G1;locus_tag=Locus1;Alias=ABC-1,LegacyABC;Dbxref=GeneID:123\n"
        << "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tID=exon1;Parent=tx1\n"
        << "chr1\tsrc\tgene\t600\t700\t.\t+\t.\tID=gene0002;Name=XYZ1\n";
    return true;
}

static std::string read_file(const std::string& path) {
    std::ifstream in{path};
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

static int run_command(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "command failed: " << command << '\n';
    }
    return status;
}

static bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

static int require_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (!contains(text, needle)) {
        std::cerr << "missing '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

static int require_not_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (contains(text, needle)) {
        std::cerr << "unexpected '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

static int compare_files(const std::string& lhs_path, const std::string& rhs_path) {
    const auto lhs = read_file(lhs_path);
    const auto rhs = read_file(rhs_path);
    if (lhs != rhs) {
        std::cerr << "output mismatch: " << lhs_path << " vs " << rhs_path << '\n';
        return 1;
    }
    return 0;
}

static void cleanup_outputs() {
    std::remove("cli_selector_smoke.gff3");
    std::remove("selector_id.gff3");
    std::remove("selector_query_id.gff3");
    std::remove("selector_name.gff3");
    std::remove("selector_alias.gff3");
    std::remove("selector_dbxref.gff3");
    std::remove("selector_name_summary.tsv");
    std::remove("selector_gene_id_summary.tsv");
    std::remove("selector_locus_tag_summary.tsv");
    std::remove("selector_alias_summary.tsv");
    std::remove("selector_dbxref_summary.tsv");
    std::remove("selector_parent.gff3");
    std::remove("selector_children.gff3");
    std::remove("selector_window_top.gff3");
    std::remove("selector_window_command.gff3");
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cli_selector_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string gff{"cli_selector_smoke.gff3"};
    if (!write_test_annotation(gff)) {
        std::cerr << "cannot write test annotation\n";
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 > selector_id.gff3") != 0 ||
        require_contains("selector_id.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_id.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " query " + gff + " --id gene0001 > selector_query_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_query_id.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name ABC1 > selector_name.gff3") != 0 ||
        require_contains("selector_name.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_name.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC > selector_alias.gff3") != 0 ||
        require_contains("selector_alias.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_alias.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 > selector_dbxref.gff3") != 0 ||
        require_contains("selector_dbxref.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_dbxref.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name ABC1 --summary-format tsv > selector_name_summary.tsv") != 0 ||
        require_contains("selector_name_summary.tsv", "ABC1\tgene0001\tName") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name G1 --summary-format tsv > selector_gene_id_summary.tsv") != 0 ||
        require_contains("selector_gene_id_summary.tsv", "G1\tgene0001\tgene_id") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name Locus1 --summary-format tsv > selector_locus_tag_summary.tsv") != 0 ||
        require_contains("selector_locus_tag_summary.tsv", "Locus1\tgene0001\tlocus_tag") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC --summary-format tsv > selector_alias_summary.tsv") != 0 ||
        require_contains("selector_alias_summary.tsv", "LegacyABC\tgene0001\tAlias") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 --summary-format tsv > selector_dbxref_summary.tsv") != 0 ||
        require_contains("selector_dbxref_summary.tsv", "GeneID:123\tgene0001\tDbxref") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --attr Parent=tx1 > selector_parent.gff3") != 0 ||
        require_contains("selector_parent.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_parent.gff3", "ID=gene0001") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --include-children > selector_children.gff3") != 0 ||
        require_contains("selector_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_top.gff3") != 0 ||
        run_command(exe + " window " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_command.gff3") != 0 ||
        compare_files("selector_window_top.gff3", "selector_window_command.gff3") != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_selector_smoke OK\n";
    return 0;
}
