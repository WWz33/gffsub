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

static bool write_qc_annotation(const std::string& path) {
    std::ofstream out{path};
    if (!out.is_open()) {
        return false;
    }

    out << "##gff-version 3\n"
        << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=dup_gene\n"
        << "chr1\tsrc\tgene\t300\t400\t.\t+\t.\tID=dup_gene\n"
        << "chr1\tsrc\tmRNA\t500\t600\t.\t+\t.\tID=orphan_tx;Parent=missing_gene\n";
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
    std::remove("cli_selector_qc.gff3");
    std::remove("cli_selector_ids.txt");
    std::remove("selector_id.gff3");
    std::remove("selector_query_id.gff3");
    std::remove("selector_attr_id.gff3");
    std::remove("selector_where_id.gff3");
    std::remove("selector_id_list.gff3");
    std::remove("selector_id_list_verbose.gff3");
    std::remove("selector_query_id_list.gff3");
    std::remove("selector_query_id_list_verbose.gff3");
    std::remove("selector_id_list_children.gff3");
    std::remove("selector_id_list_children_verbose.gff3");
    std::remove("selector_query_id_list_children.gff3");
    std::remove("selector_query_id_list_children_verbose.gff3");
    std::remove("selector_name.gff3");
    std::remove("selector_query_name.gff3");
    std::remove("selector_alias.gff3");
    std::remove("selector_dbxref.gff3");
    std::remove("selector_name_summary.tsv");
    std::remove("selector_name_summary_verbose.tsv");
    std::remove("selector_gene_id_summary.tsv");
    std::remove("selector_locus_tag_summary.tsv");
    std::remove("selector_alias_summary.tsv");
    std::remove("selector_dbxref_summary.tsv");
    std::remove("selector_parent.gff3");
    std::remove("selector_parent_attr.gff3");
    std::remove("selector_query_parent.gff3");
    std::remove("selector_query_parent_attr.gff3");
    std::remove("selector_children.gff3");
    std::remove("selector_query_children.gff3");
    std::remove("selector_query_children_short.gff3");
    std::remove("selector_children_short.gff3");
    std::remove("selector_children_type.gff3");
    std::remove("selector_query_children_type.gff3");
    std::remove("selector_region_intersection.gff3");
    std::remove("selector_help.txt");
    std::remove("selector_bed_short.bed");
    std::remove("selector_bed_format.bed");
    std::remove("selector_bed_output_format.bed");
    std::remove("selector_window_top.gff3");
    std::remove("selector_window_top_short.gff3");
    std::remove("selector_window_command.gff3");
    std::remove("selector_window_command_short.gff3");
    std::remove("selector_qc_top.tsv");
    std::remove("selector_qc_command.tsv");
    std::remove("selector_query_help.txt");
    std::remove("selector_window_help.txt");
    std::remove("selector_qc_help.txt");
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

    if (run_command(exe + " query --help > selector_query_help.txt 2>&1") != 0 ||
        require_contains("selector_query_help.txt", "Most workflows can use the top-level form") != 0 ||
        require_contains("selector_query_help.txt", "--ids FILE") != 0 ||
        require_contains("selector_query_help.txt", "Verbose alias for --ids") != 0 ||
        require_contains("selector_query_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_query_help.txt", "--summary FMT") != 0 ||
        require_contains("selector_query_help.txt", "Verbose alias for --summary") != 0) {
        return 1;
    }
    if (run_command(exe + " window --help > selector_window_help.txt 2>&1") != 0 ||
        require_contains("selector_window_help.txt", "Top-level equivalent") != 0 ||
        require_contains("selector_window_help.txt", "--up N") != 0 ||
        require_contains("selector_window_help.txt", "--down N") != 0) {
        return 1;
    }
    if (run_command(exe + " qc --help > selector_qc_help.txt 2>&1") != 0 ||
        require_contains("selector_qc_help.txt", "Top-level equivalent") != 0) {
        return 1;
    }
    if (run_command(exe + " --help > selector_help.txt 2>&1") != 0 ||
        require_contains("selector_help.txt", "--format FMT") != 0 ||
        require_contains("selector_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_help.txt", "--output-format remains as a verbose alias") != 0) {
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

    if (run_command(exe + " " + gff + " --attr ID=gene0001 > selector_attr_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_attr_id.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --where ID=gene0001 > selector_where_id.gff3") != 0 ||
        compare_files("selector_id.gff3", "selector_where_id.gff3") != 0) {
        return 1;
    }

    const std::string id_list{"cli_selector_ids.txt"};
    {
        std::ofstream out{id_list};
        if (!out.is_open()) {
            std::cerr << "cannot write ID list\n";
            return 1;
        }
        out << "gene0001\n"
            << "gene0002\n";
    }
    if (run_command(exe + " " + gff + " --ids " + id_list + " > selector_id_list.gff3") != 0 ||
        run_command(exe + " " + gff + " --id-list " + id_list + " > selector_id_list_verbose.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " > selector_query_id_list.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id-list " + id_list + " > selector_query_id_list_verbose.gff3") != 0 ||
        compare_files("selector_id_list.gff3", "selector_id_list_verbose.gff3") != 0 ||
        compare_files("selector_id_list.gff3", "selector_query_id_list.gff3") != 0 ||
        compare_files("selector_query_id_list.gff3", "selector_query_id_list_verbose.gff3") != 0 ||
        require_contains("selector_id_list.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_id_list.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --ids " + id_list + " --include-children > selector_id_list_children.gff3") != 0 ||
        run_command(exe + " " + gff + " --id-list " + id_list + " --include-children > selector_id_list_children_verbose.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " --include-children > selector_query_id_list_children.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id-list " + id_list + " --include-children > selector_query_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_query_id_list_children.gff3") != 0 ||
        compare_files("selector_query_id_list_children.gff3", "selector_query_id_list_children_verbose.gff3") != 0 ||
        require_contains("selector_id_list_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_id_list_children.gff3", "ID=exon1") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name ABC1 > selector_name.gff3") != 0 ||
        require_contains("selector_name.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_name.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --name ABC1 > selector_query_name.gff3") != 0 ||
        compare_files("selector_name.gff3", "selector_query_name.gff3") != 0) {
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

    if (run_command(exe + " " + gff + " --name ABC1 --summary tsv > selector_name_summary.tsv") != 0 ||
        run_command(exe + " " + gff + " --name ABC1 --summary-format tsv > selector_name_summary_verbose.tsv") != 0 ||
        compare_files("selector_name_summary.tsv", "selector_name_summary_verbose.tsv") != 0 ||
        require_contains("selector_name_summary.tsv", "ABC1\tgene0001\tName") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name G1 --summary tsv > selector_gene_id_summary.tsv") != 0 ||
        require_contains("selector_gene_id_summary.tsv", "G1\tgene0001\tgene_id") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name Locus1 --summary tsv > selector_locus_tag_summary.tsv") != 0 ||
        require_contains("selector_locus_tag_summary.tsv", "Locus1\tgene0001\tlocus_tag") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC --summary tsv > selector_alias_summary.tsv") != 0 ||
        require_contains("selector_alias_summary.tsv", "LegacyABC\tgene0001\tAlias") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 --summary tsv > selector_dbxref_summary.tsv") != 0 ||
        require_contains("selector_dbxref_summary.tsv", "GeneID:123\tgene0001\tDbxref") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --where Parent=tx1 > selector_parent.gff3") != 0 ||
        require_contains("selector_parent.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_parent.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --attr Parent=tx1 > selector_parent_attr.gff3") != 0 ||
        run_command(exe + " query " + gff + " --where Parent=tx1 > selector_query_parent.gff3") != 0 ||
        run_command(exe + " query " + gff + " --attr Parent=tx1 > selector_query_parent_attr.gff3") != 0 ||
        compare_files("selector_parent.gff3", "selector_parent_attr.gff3") != 0 ||
        compare_files("selector_parent.gff3", "selector_query_parent.gff3") != 0 ||
        compare_files("selector_query_parent.gff3", "selector_query_parent_attr.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --include-children > selector_children.gff3") != 0 ||
        require_contains("selector_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --include-children > selector_query_children.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 -C > selector_query_children_short.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children_short.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 -C > selector_children_short.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_children_short.gff3") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --include-children -f mRNA > selector_children_type.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id gene0001 --include-children --type mRNA > selector_query_children_type.gff3") != 0 ||
        compare_files("selector_children_type.gff3", "selector_query_children_type.gff3") != 0 ||
        require_contains("selector_children_type.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_children_type.gff3", "ID=gene0001") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --region chr1:130-140 > selector_region_intersection.gff3") != 0 ||
        require_contains("selector_region_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " -r chr1:100-400 -t bed > selector_bed_short.bed") != 0 ||
        run_command(exe + " " + gff + " -r chr1:100-400 --format bed > selector_bed_format.bed") != 0 ||
        run_command(exe + " " + gff + " -r chr1:100-400 --output-format bed > selector_bed_output_format.bed") != 0 ||
        compare_files("selector_bed_short.bed", "selector_bed_format.bed") != 0 ||
        compare_files("selector_bed_short.bed", "selector_bed_output_format.bed") != 0 ||
        require_contains("selector_bed_format.bed", "chr1\t99\t400\tgene0001") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_top.gff3") != 0 ||
        run_command(exe + " " + gff + " --id gene0001 --up 50 --down 10 > selector_window_top_short.gff3") != 0 ||
        run_command(exe + " window " + gff + " --id gene0001 --upstream 50 --downstream 10 > selector_window_command.gff3") != 0 ||
        run_command(exe + " window " + gff + " --id gene0001 --up 50 --down 10 > selector_window_command_short.gff3") != 0 ||
        compare_files("selector_window_top.gff3", "selector_window_top_short.gff3") != 0 ||
        compare_files("selector_window_top.gff3", "selector_window_command.gff3") != 0 ||
        compare_files("selector_window_command.gff3", "selector_window_command_short.gff3") != 0) {
        return 1;
    }

    const std::string qc_gff{"cli_selector_qc.gff3"};
    if (!write_qc_annotation(qc_gff)) {
        std::cerr << "cannot write QC test annotation\n";
        return 1;
    }
    if (run_command(exe + " " + qc_gff + " --qc > selector_qc_top.tsv") != 0 ||
        run_command(exe + " qc " + qc_gff + " > selector_qc_command.tsv") != 0 ||
        compare_files("selector_qc_top.tsv", "selector_qc_command.tsv") != 0 ||
        require_contains("selector_qc_top.tsv", "duplicate_id") != 0 ||
        require_contains("selector_qc_top.tsv", "missing_parent") != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_selector_smoke OK\n";
    return 0;
}
