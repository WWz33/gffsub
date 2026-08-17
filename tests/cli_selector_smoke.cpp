#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "test_utils.hpp"

using test_utils::compare_files;
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
        << "chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=gene0001;Name=ABC1;gene_id=G1;locus_tag=Locus1;Alias=ABC-1,LegacyABC;Dbxref=GeneID:123\n"
        << "chr1\tsrc\tmRNA\t100\t400\t42.5\t+\t.\tID=tx1;Parent=gene0001;Name=ABC1.1\n"
        << "chr1\tsrc\texon\t120\t180\t.\t+\t.\tID=exon1;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t150\t170\t.\t+\t0\tID=cds1;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t200\t220\t.\t+\t1\tID=cds2;Parent=tx1\n"
        << "chr1\tsrc\tCDS\t230\t250\t.\t+\t2\tID=cds3;Parent=tx1\n"
        << "chr1\tsrc\tgene\t600\t700\t.\t-\t.\tID=gene0002;Name=XYZ1;biotype=protein_coding;Note=transposon-like\n"
        << "chr2\tother\tgene\t100\t200\t.\t+\t.\tID=gene0003;Name=CHR2\n"
        << "chr2\tother\texon\t250\t280\t.\t.\t.\tID=exon2\n"
        << "chr2\tother\tmRNA\t300\t380\t.\t+\t.\tID=orphan_tx\n"
        << "chr2\tother\texon\t320\t360\t.\t+\t.\tID=orphan_exon;Parent=orphan_tx\n";
    return true;
}

static void cleanup_outputs() {
    std::remove("cli_selector_smoke.gff3");
    std::remove("cli_selector_ids.txt");
    std::remove("cli_selector_patterns.txt");
    std::remove("selector_id.gff3");
    std::remove("selector_query_id.gff3");
    std::remove("selector_attr_id.gff3");
    std::remove("selector_where_id.gff3");
    std::remove("selector_grep_id.gff3");
    std::remove("selector_grep_regex_id.gff3");
    std::remove("selector_grep_regex_seqid.gff3");
    std::remove("selector_grep_file.gff3");
    std::remove("selector_grep_invert.gff3");
    std::remove("selector_grep_ignore_case.gff3");
    std::remove("selector_include_expr_biotype.gff3");
    std::remove("selector_include_expr_logic.gff3");
    std::remove("selector_include_expr_quoted_regex.gff3");
    std::remove("selector_include_expr_numeric.gff3");
    std::remove("selector_include_expr_score.gff3");
    std::remove("selector_exclude_expr_note.gff3");
    std::remove("selector_expr_bad.out");
    std::remove("selector_expr_bad.err");
    std::remove("selector_grep_bad.out");
    std::remove("selector_grep_bad.err");
    std::remove("selector_invert_bad.out");
    std::remove("selector_invert_bad.err");
    std::remove("selector_id_list.gff3");
    std::remove("selector_id_list_verbose.gff3");
    std::remove("selector_query_id_list.gff3");
    std::remove("selector_query_id_list_verbose.gff3");
    std::remove("selector_id_list_children.gff3");
    std::remove("selector_id_list_children_alias.gff3");
    std::remove("selector_id_list_children_verbose.gff3");
    std::remove("selector_query_id_list_children.gff3");
    std::remove("selector_query_id_list_children_alias.gff3");
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
    std::remove("selector_children_alias.gff3");
    std::remove("selector_query_children.gff3");
    std::remove("selector_query_children_alias.gff3");
    std::remove("selector_query_children_short.gff3");
    std::remove("selector_children_short.gff3");
    std::remove("selector_children_type.gff3");
    std::remove("selector_query_children_type.gff3");
    std::remove("selector_parents.gff3");
    std::remove("selector_parents_alias.gff3");
    std::remove("selector_query_parents.gff3");
    std::remove("selector_query_parents_alias.gff3");
    std::remove("selector_parents_gene.gff3");
    std::remove("selector_query_parents_gene.gff3");
    std::remove("selector_parents_children.gff3");
    std::remove("selector_parents_bad.out");
    std::remove("selector_parents_bad.err");
    std::remove("selector_query_parents_bad.out");
    std::remove("selector_query_parents_bad.err");
    std::remove("selector_model.gff3");
    std::remove("selector_gene_model_alias.gff3");
    std::remove("selector_query_model.gff3");
    std::remove("selector_query_gene_model_alias.gff3");
    std::remove("selector_model_gene.gff3");
    std::remove("selector_query_model_gene.gff3");
    std::remove("selector_model_orphan_children.gff3");
    std::remove("selector_model_bad.out");
    std::remove("selector_model_bad.err");
    std::remove("selector_query_model_bad.out");
    std::remove("selector_query_model_bad.err");
    std::remove("selector_nearest.gff3");
    std::remove("selector_nearest_alias.gff3");
    std::remove("selector_query_nearest.gff3");
    std::remove("selector_query_nearest_alias.gff3");
    std::remove("selector_nearest_overlap.gff3");
    std::remove("selector_nearest_children.gff3");
    std::remove("selector_nearest_seqid_keep.gff3");
    std::remove("selector_nearest_seqid_drop.gff3");
    std::remove("selector_nearest_summary.tsv");
    std::remove("selector_nearest_not_found.tsv");
    std::remove("selector_nearest_bad.out");
    std::remove("selector_nearest_bad.err");
    std::remove("selector_query_nearest_bad.out");
    std::remove("selector_query_nearest_bad.err");
    std::remove("selector_region_intersection.gff3");
    std::remove("selector_seqid_chr2.gff3");
    std::remove("selector_seqid_gene.gff3");
    std::remove("selector_seqid_id_intersection.gff3");
    std::remove("selector_source_other.gff3");
    std::remove("selector_source_gene.gff3");
    std::remove("selector_seqid_source_intersection.gff3");
    std::remove("selector_score_value.gff3");
    std::remove("selector_score_missing.gff3");
    std::remove("selector_score_feature.gff3");
    std::remove("selector_score_bad.out");
    std::remove("selector_score_bad.err");
    std::remove("selector_score_nan.out");
    std::remove("selector_score_nan.err");
    std::remove("selector_strand_minus.gff3");
    std::remove("selector_strand_dot.gff3");
    std::remove("selector_strand_gene.gff3");
    std::remove("selector_strand_bad.out");
    std::remove("selector_strand_bad.err");
    std::remove("selector_phase_zero.gff3");
    std::remove("selector_phase_dot.gff3");
    std::remove("selector_phase_cds.gff3");
    std::remove("selector_phase_bad.out");
    std::remove("selector_phase_bad.err");
    std::remove("selector_help.txt");
    std::remove("selector_bed_short.bed");
    std::remove("selector_bed_format.bed");
    std::remove("selector_bed_output_format.bed");
    std::remove("selector_window_top.gff3");
    std::remove("selector_window_top_short.gff3");
    std::remove("selector_window_command.gff3");
    std::remove("selector_window_command_short.gff3");
    std::remove("selector_query_help.txt");
    std::remove("selector_window_help.txt");
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: cli_selector_smoke <gffsub-executable>\n";
        return 2;
    }

    const std::string exe = std::string{"\""} + argv[1] + "\"";
    const std::string exe_raw = argv[1];
    const std::string gff{"cli_selector_smoke.gff3"};
    if (!write_test_annotation(gff)) {
        std::cerr << "cannot write test annotation\n";
        return 1;
    }

    if (run_command(exe + " query --help > selector_query_help.txt 2>&1") != 0 ||
        require_contains("selector_query_help.txt", "About:") != 0 ||
        require_contains("selector_query_help.txt", "--ids FILE") != 0 ||
        require_contains("selector_query_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_query_help.txt", "--children") != 0 ||
        require_contains("selector_query_help.txt", "--parents") != 0 ||
        require_contains("selector_query_help.txt", "--model") != 0 ||
        require_contains("selector_query_help.txt", "--nearest REGION") != 0 ||
        require_contains("selector_query_help.txt", "--summary FMT") != 0) {
        return 1;
    }
    if (run_command(exe + " window --help > selector_window_help.txt 2>&1") != 0 ||
        require_contains("selector_window_help.txt", "About:") != 0 ||
        require_contains("selector_window_help.txt", "--up N") != 0 ||
        require_contains("selector_window_help.txt", "--down N") != 0) {
        return 1;
    }
    if (run_command(exe + " --help > selector_help.txt 2>&1") != 0 ||
        require_contains("selector_help.txt", "Program: gffsub") != 0 ||
        require_contains("selector_help.txt", "--format FMT") != 0 ||
        require_contains("selector_help.txt", "--where KEY=VALUE") != 0 ||
        require_contains("selector_help.txt", "--seqid LIST") != 0 ||
        require_contains("selector_help.txt", "--source SOURCE") != 0 ||
        require_contains("selector_help.txt", "--score SCORE") != 0 ||
        require_contains("selector_help.txt", "--strand STRAND") != 0 ||
        require_contains("selector_help.txt", "--phase PHASE") != 0 ||
        require_contains("selector_help.txt", "--grep FIELD:PATTERN") != 0 ||
        require_contains("selector_help.txt", "--grep-regex FIELD:REGEX") != 0 ||
        require_contains("selector_help.txt", "--grep-file FILE") != 0 ||
        require_contains("selector_help.txt", "--include-expr EXPR") != 0 ||
        require_contains("selector_help.txt", "--exclude-expr EXPR") != 0 ||
        require_contains("selector_help.txt", "--invert-match") != 0 ||
        require_contains("selector_help.txt", "--model") != 0 ||
        require_contains("selector_help.txt", "--nearest REGION") != 0 ||
        require_contains("selector_help.txt", "--output-attrs") != 0) {
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

    if (run_command(exe_raw + " " + gff + " --grep ID:gene000 > selector_grep_id.gff3") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_grep_id.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_id.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex ID:gene000[12] > selector_grep_regex_id.gff3") != 0 ||
        require_contains("selector_grep_regex_id.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_regex_id.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_grep_regex_id.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_regex_id.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex seqid:chr[0-9]+ -f gene > selector_grep_regex_seqid.gff3") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_grep_regex_seqid.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_regex_seqid.gff3", "ID=tx1") != 0) {
        return 1;
    }
    {
        std::ofstream patterns{"cli_selector_patterns.txt"};
        if (!patterns.is_open()) {
            std::cerr << "cannot write grep pattern list\n";
            return 1;
        }
        patterns << "gene0001\n"
                 << "gene0003\n";
    }
    if (run_command(exe_raw + " " + gff + " --grep-file cli_selector_patterns.txt --grep-field ID > selector_grep_file.gff3") != 0 ||
        require_contains("selector_grep_file.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_file.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_grep_file.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_grep_file.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep ID:gene000 -v > selector_grep_invert.gff3") != 0 ||
        require_contains("selector_grep_invert.gff3", "ID=tx1") != 0 ||
        require_contains("selector_grep_invert.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_grep_invert.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep name:abc1 --ignore-case > selector_grep_ignore_case.gff3") != 0 ||
        require_contains("selector_grep_ignore_case.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_grep_ignore_case.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_grep_ignore_case.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"type==gene && attr.biotype==protein_coding\" > selector_include_expr_biotype.gff3") != 0 ||
        require_contains("selector_include_expr_biotype.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_include_expr_biotype.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_include_expr_biotype.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"(type==gene && attr.biotype==protein_coding) || !seqid==chr1\" > selector_include_expr_logic.gff3") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=gene0003") != 0 ||
        require_contains("selector_include_expr_logic.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_include_expr_logic.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_include_expr_logic.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"attr.ID~\\\"gene000[13]\\\"\" > selector_include_expr_quoted_regex.gff3") != 0 ||
        require_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_include_expr_quoted_regex.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_include_expr_quoted_regex.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"type==gene && length>=101\" > selector_include_expr_numeric.gff3") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_include_expr_numeric.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_include_expr_numeric.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I \"score==42.5\" > selector_include_expr_score.gff3") != 0 ||
        require_contains("selector_include_expr_score.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_include_expr_score.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -E \"attr.Note~transposon|retroelement\" -f gene > selector_exclude_expr_note.gff3") != 0 ||
        require_contains("selector_exclude_expr_note.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_exclude_expr_note.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_exclude_expr_note.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -I type > selector_expr_bad.out 2> selector_expr_bad.err") == 0 ||
        require_contains("selector_expr_bad.err", "Error: invalid include expression") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " --grep-regex ID:[ > selector_grep_bad.out 2> selector_grep_bad.err") == 0 ||
        require_contains("selector_grep_bad.err", "Error: invalid regex in --grep-regex") != 0) {
        return 1;
    }
    if (run_command(exe_raw + " " + gff + " -v > selector_invert_bad.out 2> selector_invert_bad.err") == 0 ||
        require_contains("selector_invert_bad.err", "Error: --invert-match requires --grep, --grep-regex, or --grep-file") != 0) {
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
    if (run_command(exe + " " + gff + " --ids " + id_list + " --children > selector_id_list_children.gff3") != 0 ||
        run_command(exe + " " + gff + " --ids " + id_list + " --include-children > selector_id_list_children_alias.gff3") != 0 ||
        run_command(exe + " " + gff + " --id-list " + id_list + " --include-children > selector_id_list_children_verbose.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " --children > selector_query_id_list_children.gff3") != 0 ||
        run_command(exe + " query " + gff + " --ids " + id_list + " --include-children > selector_query_id_list_children_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id-list " + id_list + " --include-children > selector_query_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_id_list_children_alias.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_id_list_children_verbose.gff3") != 0 ||
        compare_files("selector_id_list_children.gff3", "selector_query_id_list_children.gff3") != 0 ||
        compare_files("selector_query_id_list_children.gff3", "selector_query_id_list_children_alias.gff3") != 0 ||
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

    if (run_command(exe + " " + gff + " --id gene0001 --children > selector_children.gff3") != 0 ||
        require_contains("selector_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --include-children > selector_children_alias.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_children_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --children > selector_query_children.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --include-children > selector_query_children_alias.gff3") != 0 ||
        compare_files("selector_children.gff3", "selector_query_children_alias.gff3") != 0) {
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
    if (run_command(exe + " " + gff + " --id exon1 --parents > selector_parents.gff3") != 0 ||
        require_contains("selector_parents.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_parents.gff3", "ID=tx1") != 0 ||
        require_contains("selector_parents.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_parents.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --include-parents > selector_parents_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --parents > selector_query_parents.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --include-parents > selector_query_parents_alias.gff3") != 0 ||
        compare_files("selector_parents.gff3", "selector_parents_alias.gff3") != 0 ||
        compare_files("selector_parents.gff3", "selector_query_parents.gff3") != 0 ||
        compare_files("selector_query_parents.gff3", "selector_query_parents_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --parents -f gene > selector_parents_gene.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --parents --type gene > selector_query_parents_gene.gff3") != 0 ||
        compare_files("selector_parents_gene.gff3", "selector_query_parents_gene.gff3") != 0 ||
        require_contains("selector_parents_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_parents_gene.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_parents_gene.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id tx1 --parents -C > selector_parents_children.gff3") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=exon1") != 0 ||
        require_contains("selector_parents_children.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_parents_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --parents > selector_parents_bad.out 2> selector_parents_bad.err") == 0 ||
        require_contains("selector_parents_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --parents > selector_query_parents_bad.out 2> selector_query_parents_bad.err") == 0 ||
        require_contains("selector_query_parents_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --model > selector_model.gff3") != 0 ||
        require_contains("selector_model.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_model.gff3", "ID=tx1") != 0 ||
        require_contains("selector_model.gff3", "ID=exon1") != 0 ||
        require_contains("selector_model.gff3", "ID=cds1") != 0 ||
        require_contains("selector_model.gff3", "ID=cds2") != 0 ||
        require_contains("selector_model.gff3", "ID=cds3") != 0 ||
        require_not_contains("selector_model.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --gene-model > selector_gene_model_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --model > selector_query_model.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --gene-model > selector_query_gene_model_alias.gff3") != 0 ||
        compare_files("selector_model.gff3", "selector_gene_model_alias.gff3") != 0 ||
        compare_files("selector_model.gff3", "selector_query_model.gff3") != 0 ||
        compare_files("selector_query_model.gff3", "selector_query_gene_model_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id exon1 --model -f CDS > selector_model_gene.gff3") != 0 ||
        run_command(exe + " query " + gff + " --id exon1 --model --type CDS > selector_query_model_gene.gff3") != 0 ||
        compare_files("selector_model_gene.gff3", "selector_query_model_gene.gff3") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds1") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds2") != 0 ||
        require_contains("selector_model_gene.gff3", "ID=cds3") != 0 ||
        require_not_contains("selector_model_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_model_gene.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id orphan_tx --model -C --seqid chr2 > selector_model_orphan_children.gff3") != 0 ||
        require_contains("selector_model_orphan_children.gff3", "ID=orphan_tx") != 0 ||
        require_contains("selector_model_orphan_children.gff3", "ID=orphan_exon") != 0 ||
        require_not_contains("selector_model_orphan_children.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --model > selector_model_bad.out 2> selector_model_bad.err") == 0 ||
        require_contains("selector_model_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --model > selector_query_model_bad.out 2> selector_query_model_bad.err") == 0 ||
        require_contains("selector_query_model_bad.err", "Error: --children/--parents/--model require --id, --ids, --name, --where, or --nearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 > selector_nearest.gff3") != 0 ||
        require_contains("selector_nearest.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest-gene chr1:450-500 > selector_nearest_alias.gff3") != 0 ||
        run_command(exe + " query " + gff + " --nearest chr1:450-500 > selector_query_nearest.gff3") != 0 ||
        run_command(exe + " query " + gff + " --nearest-gene chr1:450-500 > selector_query_nearest_alias.gff3") != 0 ||
        compare_files("selector_nearest.gff3", "selector_nearest_alias.gff3") != 0 ||
        compare_files("selector_nearest.gff3", "selector_query_nearest.gff3") != 0 ||
        compare_files("selector_query_nearest.gff3", "selector_query_nearest_alias.gff3") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:610-620 > selector_nearest_overlap.gff3") != 0 ||
        require_contains("selector_nearest_overlap.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_nearest_overlap.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 -C > selector_nearest_children.gff3") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=tx1") != 0 ||
        require_contains("selector_nearest_children.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_nearest_children.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --seqid chr1 > selector_nearest_seqid_keep.gff3") != 0 ||
        require_contains("selector_nearest_seqid_keep.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest_seqid_keep.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --seqid chr2 > selector_nearest_seqid_drop.gff3") != 0 ||
        require_not_contains("selector_nearest_seqid_drop.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_nearest_seqid_drop.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 --summary tsv > selector_nearest_summary.tsv") != 0 ||
        require_contains("selector_nearest_summary.tsv", "chr1:450-500\tgene0001\tnearest") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr9:1-100 --summary tsv > selector_nearest_not_found.tsv") != 0 ||
        require_contains("selector_nearest_not_found.tsv", "chr9:1-100\t\tnearest") != 0 ||
        require_contains("selector_nearest_not_found.tsv", "not_found") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr1-450-500 > selector_nearest_bad.out 2> selector_nearest_bad.err") == 0 ||
        require_contains("selector_nearest_bad.err", "Error: invalid nearest region format chr1-450-500") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --nearest chr1-450-500 > selector_query_nearest_bad.out 2> selector_query_nearest_bad.err") == 0 ||
        require_contains("selector_query_nearest_bad.err", "Error: invalid nearest region format chr1-450-500") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --id gene0001 --region chr1:130-140 > selector_region_intersection.gff3") != 0 ||
        require_contains("selector_region_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_region_intersection.gff3", "ID=gene0002") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr2 > selector_seqid_chr2.gff3") != 0 ||
        require_contains("selector_seqid_chr2.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_seqid_chr2.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr1 -f gene > selector_seqid_gene.gff3") != 0 ||
        require_contains("selector_seqid_gene.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_seqid_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_seqid_gene.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_seqid_gene.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --id gene0001 --seqid chr2 > selector_seqid_id_intersection.gff3") != 0 ||
        require_not_contains("selector_seqid_id_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_seqid_id_intersection.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --source other > selector_source_other.gff3") != 0 ||
        require_contains("selector_source_other.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_source_other.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --source src -f gene > selector_source_gene.gff3") != 0 ||
        require_contains("selector_source_gene.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_source_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_source_gene.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_source_gene.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --seqid chr2 --source src > selector_seqid_source_intersection.gff3") != 0 ||
        require_not_contains("selector_seqid_source_intersection.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_seqid_source_intersection.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score 42.5 > selector_score_value.gff3") != 0 ||
        require_contains("selector_score_value.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_score_value.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score . > selector_score_missing.gff3") != 0 ||
        require_contains("selector_score_missing.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_score_missing.gff3", "ID=tx1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score . -f gene > selector_score_feature.gff3") != 0 ||
        require_contains("selector_score_feature.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_score_feature.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_score_feature.gff3", "ID=exon1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score abc > selector_score_bad.out 2> selector_score_bad.err") == 0 ||
        require_contains("selector_score_bad.err", "Error: --score expects a finite floating point number or .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --score nan > selector_score_nan.out 2> selector_score_nan.err") == 0 ||
        require_contains("selector_score_nan.err", "Error: --score expects a finite floating point number or .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand - > selector_strand_minus.gff3") != 0 ||
        require_contains("selector_strand_minus.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_strand_minus.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand . > selector_strand_dot.gff3") != 0 ||
        require_contains("selector_strand_dot.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_strand_dot.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand - -f gene > selector_strand_gene.gff3") != 0 ||
        require_contains("selector_strand_gene.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_strand_gene.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_strand_gene.gff3", "ID=exon2") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --strand forward > selector_strand_bad.out 2> selector_strand_bad.err") == 0 ||
        require_contains("selector_strand_bad.err", "Error: --strand expects one of +, -, ., ?") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 0 > selector_phase_zero.gff3") != 0 ||
        require_contains("selector_phase_zero.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_phase_zero.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_phase_zero.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase . > selector_phase_dot.gff3") != 0 ||
        require_contains("selector_phase_dot.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_phase_dot.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_phase_dot.gff3", "ID=cds1") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 1 -f CDS > selector_phase_cds.gff3") != 0 ||
        require_contains("selector_phase_cds.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 3 > selector_phase_bad.out 2> selector_phase_bad.err") == 0 ||
        require_contains("selector_phase_bad.err", "Error: --phase expects one of 0, 1, 2, .") != 0) {
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

    cleanup_outputs();
    std::cout << "cli_selector_smoke OK\n";
    return 0;
}
