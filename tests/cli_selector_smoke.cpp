#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "test_utils.hpp"

using test_utils::compare_files;
using test_utils::expect_command_failure;
using test_utils::read_file;
using test_utils::require_contains;
using test_utils::require_exit_one_with_error;
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
    std::remove("selector_overlap.bed");
    std::remove("selector_strand2.bed");
    std::remove("selector_window_patterns.txt");
    std::remove("selector_odd.gff3");
    std::remove("selector_odd.gff3.gz");
    std::remove("selector_odd_fifo");
    std::remove("odd_plain.gff3");
    std::remove("odd_gz.gff3");
    std::remove("odd_stdin.gff3");
    std::remove("odd_fifo.gff3");
    std::remove("pad_where.gff3");
    std::remove("pad_attrs.gff3");
    std::remove("selector_bed_overlap.gff3");
    std::remove("selector_bed_exclude.gff3");
    std::remove("selector_bed_same.gff3");
    std::remove("selector_bed_opp.gff3");
    std::remove("selector_strand_bad.out");
    std::remove("selector_strand_bad.err");
    std::remove("selector_strand2_any.gff3");
    std::remove("selector_strand2_same.gff3");
    std::remove("selector_strand2_opp.gff3");
    std::remove("selector_excl_region.gff3");
    std::remove("selector_excl_touch.gff3");
    std::remove("selector_mRNA_nodrop.gff3");
    std::remove("selector_mRNA_drop.gff3");
    std::remove("selector_exon_drop.gff3");
    std::remove("selector_exon_mRNA_drop.gff3");
    std::remove("selector_out_attrs.gff3");
    std::remove("selector_out_attrs_note.gff3");
    std::remove("selector_window.gff3");
    std::remove("win_plain.gff3");
    std::remove("win_strand.gff3");
    std::remove("win_mrna.gff3");
    std::remove("selector_window_plus.gff3");
    std::remove("win_plus_plain.gff3");
    std::remove("win_plus_a.gff3");
    std::remove("selector_gz_probe.gz");
    std::remove("selector_gz_fifo");
    std::remove("gz_ref.gff3");
    std::remove("gz_file.gff3");
    std::remove("gz_stdin.gff3");
    std::remove("gz_fifo.gff3");
    std::remove("gz_trunc.out");
    std::remove("gz_trunc.err");
    std::remove("window_bad.out");
    std::remove("window_bad.err");
    std::remove("selector_bad.out");
    std::remove("selector_bad.err");
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
    std::remove("cli_selector_tie.gff3");
    std::remove("cli_selector_tie_rev.gff3");
    std::remove("selector_tie_a.gff3");
    std::remove("selector_tie_b.gff3");
    std::remove("selector_type_repeat.gff3");
    std::remove("selector_type_list.gff3");
    std::remove("selector_query_type_list.gff3");
    std::remove("selector_query_type_exclude.gff3");
    std::remove("cli_selector_nc.gff3");
    std::remove("selector_longest_nc_auto.gff3");
    std::remove("selector_lt_requires.out");
    std::remove("selector_lt_requires.err");
    std::remove("selector_longest_nc_type.gff3");
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
        require_contains("selector_query_help.txt", "-s, --summary") != 0) {
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
        require_contains("selector_help.txt", "-s, --summary") != 0) {
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
    if (run_command(exe_raw + " " + gff + " --grep-regex seqid:chr[0-9]+ -t gene > selector_grep_regex_seqid.gff3") != 0 ||
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
    if (run_command(exe_raw + " " + gff + " -E \"attr.Note~transposon|retroelement\" -t gene > selector_exclude_expr_note.gff3") != 0 ||
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

    if (run_command(exe + " " + gff + " --name ABC1 -s > selector_name_summary.tsv") != 0 ||
        run_command(exe + " " + gff + " --name ABC1 -s > selector_name_summary_verbose.tsv") != 0 ||
        compare_files("selector_name_summary.tsv", "selector_name_summary_verbose.tsv") != 0 ||
        require_contains("selector_name_summary.tsv", "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len") != 0 ||
        require_contains("selector_name_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name G1 -s > selector_gene_id_summary.tsv") != 0 ||
        require_contains("selector_gene_id_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name Locus1 -s > selector_locus_tag_summary.tsv") != 0 ||
        require_contains("selector_locus_tag_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name LegacyABC -s > selector_alias_summary.tsv") != 0 ||
        require_contains("selector_alias_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }

    if (run_command(exe + " " + gff + " --name GeneID:123 -s > selector_dbxref_summary.tsv") != 0 ||
        require_contains("selector_dbxref_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
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

    if (run_command(exe + " " + gff + " --id gene0001 --include-children -t mRNA > selector_children_type.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --id exon1 --parents -t gene > selector_parents_gene.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --id exon1 --model -t CDS > selector_model_gene.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --nearest chr1:450-500 -s > selector_nearest_summary.tsv") != 0 ||
        require_contains("selector_nearest_summary.tsv", "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len") != 0 ||
        require_contains("selector_nearest_summary.tsv", "chr1\tgene\t1\t301\t301\t301\t301") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --nearest chr9:1-100 -s > selector_nearest_not_found.tsv") != 0 ||
        require_contains("selector_nearest_not_found.tsv", "seqid\ttype\tcount\tsum_len\tmin_len\tavg_len\tmax_len") != 0 ||
        require_not_contains("selector_nearest_not_found.tsv", "chr9") != 0) {
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

    // nearest tie: two genes equidistant from the query region must resolve
    // deterministically to the one with the smaller start, independent of
    // record order in the file.
    {
        std::ofstream tie{"cli_selector_tie.gff3"};
        if (!tie.is_open()) return 1;
        tie << "##gff-version 3\n"
            << "chr1\tsrc\tgene\t400\t500\t.\t+\t.\tID=gb\n"
            << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=ga\n";
    }
    {
        std::ofstream tie_rev{"cli_selector_tie_rev.gff3"};
        if (!tie_rev.is_open()) return 1;
        tie_rev << "##gff-version 3\n"
                << "chr1\tsrc\tgene\t100\t200\t.\t+\t.\tID=ga\n"
                << "chr1\tsrc\tgene\t400\t500\t.\t+\t.\tID=gb\n";
    }
    if (run_command(exe + " cli_selector_tie.gff3 --nearest chr1:280-320 > selector_tie_a.gff3") != 0 ||
        run_command(exe + " cli_selector_tie_rev.gff3 --nearest chr1:280-320 > selector_tie_b.gff3") != 0 ||
        compare_files("selector_tie_a.gff3", "selector_tie_b.gff3") != 0 ||
        require_contains("selector_tie_a.gff3", "ID=ga") != 0 ||
        require_not_contains("selector_tie_a.gff3", "ID=gb") != 0) {
        return 1;
    }

    // -t repeatable: -t A -t B equals -t A,B
    if (run_command(exe + " " + gff + " -t gene -t mRNA > selector_type_repeat.gff3") != 0 ||
        run_command(exe + " " + gff + " -t gene,mRNA > selector_type_list.gff3") != 0 ||
        compare_files("selector_type_repeat.gff3", "selector_type_list.gff3") != 0 ||
        require_contains("selector_type_repeat.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_type_repeat.gff3", "ID=tx1") != 0 ||
        require_not_contains("selector_type_repeat.gff3", "ID=exon1") != 0) {
        return 1;
    }

    // query -t supports comma list and ^ exclusion (same semantics as main path)
    if (run_command(exe + " query " + gff + " --id gene0001 --children -t mRNA,exon > selector_query_type_list.gff3") != 0 ||
        require_contains("selector_query_type_list.gff3", "ID=tx1") != 0 ||
        require_contains("selector_query_type_list.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_query_type_list.gff3", "ID=cds1") != 0) {
        return 1;
    }
    if (run_command(exe + " query " + gff + " --id gene0001 --children --type ^mRNA > selector_query_type_exclude.gff3") != 0 ||
        require_contains("selector_query_type_exclude.gff3", "ID=exon1") != 0 ||
        require_contains("selector_query_type_exclude.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_query_type_exclude.gff3", "ID=tx1") != 0) {
        return 1;
    }

    // --longest-type: isoform selection independent of output type filter
    {
        std::ofstream nc{"cli_selector_nc.gff3"};
        if (!nc.is_open()) return 1;
        nc << "##gff-version 3\n"
           << "chr1\tsrc\tgene\t100\t500\t.\t+\t.\tID=gnc\n"
           << "chr1\tsrc\tncRNA\t100\t300\t.\t+\t.\tID=n1;Parent=gnc\n"
           << "chr1\tsrc\tncRNA\t200\t500\t.\t+\t.\tID=n2;Parent=gnc\n"
           << "chr1\tsrc\texon\t100\t300\t.\t+\t.\tID=ne1;Parent=n1\n"
           << "chr1\tsrc\texon\t200\t500\t.\t+\t.\tID=ne2;Parent=n2\n";
    }
    // auto-detect now recognizes ncRNA via transcript class
    if (run_command(exe + " cli_selector_nc.gff3 --longest > selector_longest_nc_auto.gff3") != 0 ||
        require_contains("selector_longest_nc_auto.gff3", "ID=n2") != 0 ||
        require_not_contains("selector_longest_nc_auto.gff3", "ID=n1") != 0 ||
        require_contains("selector_longest_nc_auto.gff3", "ID=ne2") != 0) {
        return 1;
    }
    // --longest-type requires --longest
    if (run_command(exe + " cli_selector_nc.gff3 --longest-type ncRNA > selector_lt_requires.out 2> selector_lt_requires.err") == 0 ||
        require_contains("selector_lt_requires.err", "Error: --longest-type requires --longest") != 0) {
        return 1;
    }
    // --longest --longest-type + -t output filter compose independently
    if (run_command(exe + " cli_selector_nc.gff3 --longest --longest-type ncRNA -t ncRNA > selector_longest_nc_type.gff3") != 0 ||
        require_contains("selector_longest_nc_type.gff3", "ID=n2") != 0 ||
        require_not_contains("selector_longest_nc_type.gff3", "ID=gnc") != 0 ||
        require_not_contains("selector_longest_nc_type.gff3", "ID=ne2") != 0) {
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
    if (run_command(exe + " " + gff + " --seqid chr1 -t gene > selector_seqid_gene.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --source src -t gene > selector_source_gene.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --score . -t gene > selector_score_feature.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --strand - -t gene > selector_strand_gene.gff3") != 0 ||
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
    if (run_command(exe + " " + gff + " --phase 1 -t CDS > selector_phase_cds.gff3") != 0 ||
        require_contains("selector_phase_cds.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=cds1") != 0 ||
        require_not_contains("selector_phase_cds.gff3", "ID=gene0001") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " --phase 3 > selector_phase_bad.out 2> selector_phase_bad.err") == 0 ||
        require_contains("selector_phase_bad.err", "Error: --phase expects one of 0, 1, 2, .") != 0) {
        return 1;
    }
    if (run_command(exe + " " + gff + " -r chr1:100-400 --format bed > selector_bed_short.bed") != 0 ||
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

    // --- complement and strand-aware BED overlap ---

    // A small BED set: one chr1 interval on +, one chr1 interval on -,
    // one chr2 interval. Written by hand here because the tests below
    // assert on its exact contents.
    {
        std::ofstream bed{"selector_overlap.bed"};
        bed << "chr1\t99\t150\tp\t0\t+\n"
            << "chr1\t599\t800\tm\t0\t-\n"
            << "chr2\t99\t150\tc2\t0\t+\n";
    }
    // -b keeps overlapping records on both strands.
    if (run_command(exe + " " + gff + " -b selector_overlap.bed > selector_bed_overlap.gff3") != 0 ||
        require_contains("selector_bed_overlap.gff3", "ID=gene0001") != 0 ||   // chr1 100-400 hits [99,150)
        require_contains("selector_bed_overlap.gff3", "ID=gene0002") != 0 ||   // chr1 600-700 hits [599,800)
        require_contains("selector_bed_overlap.gff3", "ID=gene0003") != 0 ||   // chr2 100-200 hits [99,150)
        require_not_contains("selector_bed_overlap.gff3", "ID=orphan_tx") != 0) {  // chr2 300-380 misses
        return 1;
    }
    // -b ^FILE drops exactly what -b keeps.
    if (run_command(exe + " " + gff + " -b '^selector_overlap.bed' > selector_bed_exclude.gff3") != 0 ||
        require_not_contains("selector_bed_exclude.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_bed_exclude.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_bed_exclude.gff3", "ID=gene0003") != 0 ||
        require_contains("selector_bed_exclude.gff3", "ID=orphan_tx") != 0) {
        return 1;
    }
    // --same-strand: gene0001 (+ over the + interval) and gene0002
    // (- over the - interval) both in; the chr2 interval is + and
    // gene0003 is +, so it stays too.
    if (run_command(exe + " " + gff + " -b selector_overlap.bed --same-strand > selector_bed_same.gff3") != 0 ||
        require_contains("selector_bed_same.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_bed_same.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_bed_same.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    // --opposite-strand: no record overlaps an interval of the opposite
    // strand (the two chr1 candidates match same-strand intervals, chr2
    // is + over +), so nothing survives but the header.
    if (run_command(exe + " " + gff + " -b selector_overlap.bed --opposite-strand > selector_bed_opp.gff3") != 0 ||
        require_not_contains("selector_bed_opp.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_bed_opp.gff3", "ID=gene0002") != 0 ||
        require_not_contains("selector_bed_opp.gff3", "ID=gene0003") != 0 ||
        require_not_contains("selector_bed_opp.gff3", "ID=orphan_tx") != 0) {
        return 1;
    }
    // Strand modes require -b.
    if (run_command(exe + " " + gff + " --same-strand > selector_strand_bad.out 2> selector_strand_bad.err") == 0 ||
        require_contains("selector_strand_bad.err", "Error:") != 0) {
        return 1;
    }

    // Second strand fixture: exercises a '.'-strand record (exon2, overlaps
    // the chr2 interval but has no strand) and gives --opposite-strand a
    // positive control (gene0001/mRNA/cds2 are + and overlap the - interval).
    {
        std::ofstream bed{"selector_strand2.bed"};
        bed << "chr1\t199\t260\topp\t0\t-\n"
            << "chr2\t249\t281\td\t0\t+\n";
    }
    // Pattern file for the window-guard matrix below (valid so the command
    // passes cross-validation and actually reaches the guard).
    {
        std::ofstream pat{"selector_window_patterns.txt"};
        pat << "gene\n";
    }
    // No strand mode: exon2 overlaps and survives, proving the exclusion in
    // the next two cases is due to its strand, not its position.
    if (run_command(exe + " " + gff + " -b selector_strand2.bed > selector_strand2_any.gff3") != 0 ||
        require_contains("selector_strand2_any.gff3", "ID=exon2") != 0 ||
        require_contains("selector_strand2_any.gff3", "ID=cds2") != 0) {
        return 1;
    }
    // --same-strand: nothing matches (the + records overlap only the -
    // interval; exon2's '.' never matches).
    if (run_command(exe + " " + gff + " -b selector_strand2.bed --same-strand > selector_strand2_same.gff3") != 0 ||
        require_not_contains("selector_strand2_same.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_strand2_same.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_strand2_same.gff3", "ID=exon2") != 0) {
        return 1;
    }
    // --opposite-strand: the + records match the - interval (positive
    // control); exon2 still drops out, so the dot-strand exclusion is
    // exercised by a record that demonstrably overlaps.
    if (run_command(exe + " " + gff + " -b selector_strand2.bed --opposite-strand > selector_strand2_opp.gff3") != 0 ||
        require_contains("selector_strand2_opp.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_strand2_opp.gff3", "ID=cds2") != 0 ||
        require_not_contains("selector_strand2_opp.gff3", "ID=exon2") != 0) {
        return 1;
    }

    // --- --exclude-region ---

    if (run_command(exe + " " + gff + " --exclude-region chr1:100-400 > selector_excl_region.gff3") != 0 ||
        require_not_contains("selector_excl_region.gff3", "ID=gene0001") != 0 ||
        require_not_contains("selector_excl_region.gff3", "ID=tx1") != 0 ||
        require_contains("selector_excl_region.gff3", "ID=gene0002") != 0 ||
        require_contains("selector_excl_region.gff3", "ID=gene0003") != 0) {
        return 1;
    }
    // A touching-but-not-overlapping region keeps the feature (1-based
    // inclusive): 401-500 does not overlap gene0001's 100-400.
    if (run_command(exe + " " + gff + " --exclude-region chr1:401-500 > selector_excl_touch.gff3") != 0 ||
        require_contains("selector_excl_touch.gff3", "ID=gene0001") != 0 ||
        require_contains("selector_excl_touch.gff3", "ID=gene0002") != 0) {
        return 1;
    }

    // --- --drop-orphans ---

    // Without --drop-orphans the type filter is the only filter.
    if (run_command(exe + " " + gff + " -t mRNA > selector_mRNA_nodrop.gff3") != 0 ||
        require_contains("selector_mRNA_nodrop.gff3", "ID=tx1") != 0 ||
        require_contains("selector_mRNA_nodrop.gff3", "ID=orphan_tx") != 0) {
        return 1;
    }
    // With it, tx1 is dropped: its Parent=gene0001 left the kept set when
    // the type filter removed the gene. orphan_tx has no Parent, so it
    // is never an orphan and survives.
    if (run_command(exe + " " + gff + " -t mRNA --drop-orphans > selector_mRNA_drop.gff3") != 0 ||
        require_not_contains("selector_mRNA_drop.gff3", "ID=tx1") != 0 ||
        require_contains("selector_mRNA_drop.gff3", "ID=orphan_tx") != 0) {
        return 1;
    }
    // Fixpoint: -t exon drops every parented exon (exon1's tx1 is gone,
    // orphan_exon's orphan_tx is gone); only the parentless exon2 stays.
    if (run_command(exe + " " + gff + " -t exon --drop-orphans > selector_exon_drop.gff3") != 0 ||
        require_contains("selector_exon_drop.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_exon_drop.gff3", "ID=exon1") != 0 ||
        require_not_contains("selector_exon_drop.gff3", "ID=orphan_exon") != 0) {
        return 1;
    }
    // Keeping mRNA too lets orphan_tx (no Parent) survive, and its child
    // orphan_exon with it: one extra fixpoint iteration.
    if (run_command(exe + " " + gff + " -t exon,mRNA --drop-orphans > selector_exon_mRNA_drop.gff3") != 0 ||
        require_contains("selector_exon_mRNA_drop.gff3", "ID=orphan_tx") != 0 ||
        require_contains("selector_exon_mRNA_drop.gff3", "ID=orphan_exon") != 0 ||
        require_contains("selector_exon_mRNA_drop.gff3", "ID=exon2") != 0 ||
        require_not_contains("selector_exon_mRNA_drop.gff3", "ID=tx1") != 0) {
        return 1;
    }

    // --- --out-attrs ---

    if (run_command(exe + " " + gff + " --out-attrs Name > selector_out_attrs.gff3") != 0 ||
        require_contains("selector_out_attrs.gff3", "ID=gene0001;Name=ABC1") != 0 ||
        require_contains("selector_out_attrs.gff3", "ID=tx1;Parent=gene0001;Name=ABC1.1") != 0 ||
        require_not_contains("selector_out_attrs.gff3", "locus_tag") != 0 ||
        require_not_contains("selector_out_attrs.gff3", "biotype") != 0) {
        return 1;
    }
    // A tag not present on a record simply drops out; ID/Parent always survive.
    if (run_command(exe + " " + gff + " --out-attrs Note,nonexistent > selector_out_attrs_note.gff3") != 0 ||
        require_contains("selector_out_attrs_note.gff3", "ID=gene0002;Note=transposon-like") != 0 ||
        require_contains("selector_out_attrs_note.gff3", "ID=exon1;Parent=tx1") != 0 ||
        require_not_contains("selector_out_attrs_note.gff3", "Name=") != 0) {
        return 1;
    }

    // --- window-shortcut guard rejects every non-window flag ---

    // The window path (triggered by --up/--down/--strand-aware with one -i)
    // only accepts -i/-u/-D/-a. Each flag below must be rejected; this is a
    // blacklist guard that needs every future CliArgs field added to it, so
    // we cover the surface explicitly.
    {
        const std::string wbase = exe + " " + gff + " -i gene0001 -u 50 ";
        const std::vector<std::string> reject = {
            "--ids nope.txt",          "-n ABC1",           "-w ID=gene0001",
            "--grep type:gene",        "--grep-regex type:.+",
            "--grep-file selector_window_patterns.txt --grep-field type",
            "--grep-file selector_window_patterns.txt --grep-field type --grep-file-regex",
            "-I type==gene",           "-E type==gene",
            "-v --grep type:gene",     "-y",
            "-C",                      "-p",                "-m",
            "-N chr1:1-100",           "-s",
            "-r chr1:1-100",           "-b selector_overlap.bed",
            "--exclude-region chr1:1-100",
            "-b selector_overlap.bed --same-strand",
            "-b selector_overlap.bed --opposite-strand",
            "--drop-orphans",          "--out-attrs Name",
            "-S chr1",                 "--source src",
            "-c 0",                    "--strand +",        "--phase 0",
            "-t gene",                 "-L",                "-L --longest-type mRNA",
            "-@ 2",                    "-k start",          "-R",
            "-f gtf",                  "-o /dev/null",
        };
        for (const auto& extra : reject) {
            if (expect_command_failure(wbase + extra + " > window_bad.out 2> window_bad.err") != 0 ||
                // Assert on the guard's own message, not just "Error:": a row
                // that fails earlier (e.g. cross-validation) would otherwise
                // pass without exercising the guard clause.
                require_contains("window_bad.err", "window shortcut only supports") != 0) {
                std::cerr << "window guard failed to reject: " << extra << '\n';
                return 1;
            }
        }
    }

    // --- window subcommand functionality ---

    {
        std::ofstream wf{"selector_window.gff3"};
        wf << "##gff-version 3\n"
           << "chr1\t.\tgene\t5000\t6000\t.\t-\t.\tID=wg\n"
           << "chr1\t.\tmRNA\t5000\t6000\t.\t-\t.\tID=wt;Parent=wg\n"
           << "chr1\t.\texon\t4930\t4980\t.\t-\t.\tID=wleft;Parent=wt\n"
           << "chr1\t.\texon\t6050\t6080\t.\t-\t.\tID=wright;Parent=wt\n";
    }
    // Plain window around wg: [4900,6010] keeps the left exon, not the right.
    if (run_command(exe + " window selector_window.gff3 -i wg -u 100 -D 10 > win_plain.gff3") != 0 ||
        require_contains("win_plain.gff3", "ID=wleft") != 0 ||
        require_not_contains("win_plain.gff3", "ID=wright") != 0) {
        return 1;
    }
    // Strand-aware on the minus strand: upstream extends toward larger
    // coordinates, so the window is [4990,6100] and the right exon replaces
    // the left one.
    if (run_command(exe + " window selector_window.gff3 -i wg -u 100 -D 10 -a > win_strand.gff3") != 0 ||
        require_contains("win_strand.gff3", "ID=wright") != 0 ||
        require_not_contains("win_strand.gff3", "ID=wleft") != 0) {
        return 1;
    }
    // The selector takes any record ID, not only a gene name.
    if (run_command(exe + " window selector_window.gff3 -i wt -u 10 -D 10 > win_mrna.gff3") != 0 ||
        require_contains("win_mrna.gff3", "ID=wt") != 0 ||
        require_not_contains("win_mrna.gff3", "ID=wleft") != 0) {
        return 1;
    }
    // Plus strand: -a must NOT swap the extensions (only minus-strand genes
    // extend toward larger coordinates), so -a and plain agree.
    {
        std::ofstream wf{"selector_window_plus.gff3"};
        wf << "##gff-version 3\n"
           << "chr1\t.\tgene\t8000\t9000\t.\t+\t.\tID=pg\n"
           << "chr1\t.\texon\t7930\t7980\t.\t+\t.\tID=pleft;Parent=pg\n"
           << "chr1\t.\texon\t9050\t9080\t.\t+\t.\tID=pright;Parent=pg\n";
    }
    if (run_command(exe + " window selector_window_plus.gff3 -i pg -u 100 -D 10 > win_plus_plain.gff3") != 0 ||
        run_command(exe + " window selector_window_plus.gff3 -i pg -u 100 -D 10 -a > win_plus_a.gff3") != 0 ||
        compare_files("win_plus_plain.gff3", "win_plus_a.gff3") != 0 ||
        require_contains("win_plus_a.gff3", "ID=pleft") != 0 ||
        require_not_contains("win_plus_a.gff3", "ID=pright") != 0) {
        return 1;
    }

    // --- gzip: file, stdin, and FIFO all agree ---

    const bool have_gz_tools =
        std::system("command -v gzip > /dev/null 2>&1") == 0 &&
        std::system("command -v mkfifo > /dev/null 2>&1") == 0;
    if (!have_gz_tools) {
        std::cerr << "skipping gzip tests: gzip or mkfifo not available\n";
    }
    if (have_gz_tools) {
        if (run_command(exe + " " + gff + " > gz_ref.gff3") != 0 ||
            run_command("gzip -c " + gff + " > selector_gz_probe.gz") != 0 ||
            run_command(exe + " selector_gz_probe.gz > gz_file.gff3") != 0 ||
            run_command("cat selector_gz_probe.gz | " + exe + " - > gz_stdin.gff3") != 0 ||
            compare_files("gz_ref.gff3", "gz_file.gff3") != 0 ||
            compare_files("gz_ref.gff3", "gz_stdin.gff3") != 0) {
            return 1;
        }
        // FIFO carrying gzip data: the writer starts first and blocks until
        // the reader opens the pipe.
        if (run_command("rm -f selector_gz_fifo") != 0 ||
            run_command("mkfifo selector_gz_fifo") != 0 ||
            run_command("(cat selector_gz_probe.gz > selector_gz_fifo &) ; " + exe +
                        " selector_gz_fifo > gz_fifo.gff3") != 0 ||
            compare_files("gz_ref.gff3", "gz_fifo.gff3") != 0) {
            return 1;
        }
        // Comment line with tabs (a pre-fix mis-sniff trigger), a padded
        // key after '; ', and a final line ending in CR without LF: all
        // routes must agree and keep both records.
        {
            std::ofstream odd{"selector_odd.gff3"};
            odd << "#comment\twith\ttab\n"
                << "##gff-version 3\n"
                << "chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g1; Name=n1\n"
                << "chr1\t.\texon\t2\t3\t.\t+\t.\tID=e1;Parent=g1\r";
        }
        if (run_command("gzip -c selector_odd.gff3 > selector_odd.gff3.gz") != 0 ||
            run_command(exe + " selector_odd.gff3 > odd_plain.gff3") != 0 ||
            run_command(exe + " selector_odd.gff3.gz > odd_gz.gff3") != 0 ||
            run_command("cat selector_odd.gff3.gz | " + exe + " - > odd_stdin.gff3") != 0 ||
            run_command("(cat selector_odd.gff3.gz > selector_odd_fifo &) ; " + exe +
                        " selector_odd_fifo > odd_fifo.gff3") != 0 ||
            compare_files("odd_plain.gff3", "odd_gz.gff3") != 0 ||
            compare_files("odd_plain.gff3", "odd_stdin.gff3") != 0 ||
            compare_files("odd_plain.gff3", "odd_fifo.gff3") != 0 ||
            require_contains("odd_gz.gff3", "ID=e1") != 0) {
            return 1;
        }
        // `; Key=value` spacing: the attribute indexer and --out-attrs must
        // agree on the key.
        if (run_command(exe + " selector_odd.gff3 --where Name=n1 > pad_where.gff3") != 0 ||
            run_command(exe + " selector_odd.gff3 --out-attrs Name > pad_attrs.gff3") != 0 ||
            require_contains("pad_where.gff3", "ID=g1") != 0 ||
            require_contains("pad_attrs.gff3", "Name=n1") != 0) {
            return 1;
        }
        // Truncated gzip on stdin must fail, not silently emit nothing.
        if (expect_command_failure("gzip -c " + gff + " | head -c 40 | " + exe +
                                   " - > gz_trunc.out 2> gz_trunc.err") != 0 ||
            require_contains("gz_trunc.err", "Error:") != 0) {
            return 1;
        }
    }

    // --- empty-value rejection on the guarded flags ---

    if (require_exit_one_with_error(exe + " " + gff + " --nearest '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -N '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -k '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -o '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " --grep-file '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " --grep-field '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " --seqid '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -n '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -i '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -t '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -r '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " --exclude-region '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " --ids '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0 ||
        require_exit_one_with_error(exe + " " + gff + " -L --longest-type '' > selector_bad.out 2> selector_bad.err",
                                    "selector_bad.err", "non-empty") != 0) {
        return 1;
    }

    cleanup_outputs();
    std::cout << "cli_selector_smoke OK\n";
    return 0;
}
