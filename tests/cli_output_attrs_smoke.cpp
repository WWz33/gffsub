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
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123\n"
        << "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n";
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

static int compare_files(const std::string& lhs_path, const std::string& rhs_path) {
    const auto lhs = read_file(lhs_path);
    const auto rhs = read_file(rhs_path);
    if (lhs != rhs) {
        std::cerr << "output mismatch: " << lhs_path << " vs " << rhs_path << '\n';
        return 1;
    }
    return 0;
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

    const std::string keys{"ID,Name,Alias,Dbxref"};
    if (run_command(exe + " " + gff + " --id gene0001 --output-attrs " + keys + " > output_attrs_main.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --out-attrs " + keys + " > output_attrs_short.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --attrs " + keys + " > output_attrs_old.tsv") != 0) {
        return 1;
    }
    if (compare_files("output_attrs_main.tsv", "output_attrs_short.tsv") != 0 ||
        compare_files("output_attrs_main.tsv", "output_attrs_old.tsv") != 0) {
        return 1;
    }

    if (run_command(exe + " query " + gff + " --id gene0001 --output-attrs " + keys + " > query_output_attrs_main.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --out-attrs " + keys + " > query_output_attrs_short.tsv") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --attrs " + keys + " > query_output_attrs_old.tsv") != 0) {
        return 1;
    }
    if (compare_files("query_output_attrs_main.tsv", "query_output_attrs_short.tsv") != 0 ||
        compare_files("query_output_attrs_main.tsv", "query_output_attrs_old.tsv") != 0) {
        return 1;
    }

    std::cout << "cli_output_attrs_smoke OK\n";
    return 0;
}
