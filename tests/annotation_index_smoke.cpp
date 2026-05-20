#include "gff3.hpp"

#include <fstream>
#include <iostream>
#include <string>

static int check_soybean_annotation(const std::string& path) {
    const std::string gene_id{"SoyL04_01G000000"};
    const auto index = gffsub::AnnotationIndex::from_gff3(path);

    const auto gene = index.find_by_id(gene_id);
    if (!gene || gene->type != "gene" || gene->seqid != "Chr01") {
        std::cerr << "find_by_id failed for " << gene_id << '\n';
        return 1;
    }

    const auto gene_by_lookup = index.find_gene(gene_id);
    if (!gene_by_lookup || gene_by_lookup->id != gene->id) {
        std::cerr << "find_gene failed for " << gene_id << '\n';
        return 1;
    }

    const auto children = index.children_of(gene_id);
    if (children.size() != 1 || children.front().type != "mRNA") {
        std::cerr << "children_of failed for " << gene_id << '\n';
        return 1;
    }

    const auto parents = index.parents_of("SoyL04_01G000000.m1");
    if (parents.size() != 1 || parents.front().id != gene->id) {
        std::cerr << "parents_of failed for SoyL04_01G000000.m1\n";
        return 1;
    }

    const auto descendants = index.descendants_of(gene_id);
    if (descendants.size() != 7) {
        std::cerr << "descendants_of failed for " << gene_id << '\n';
        return 1;
    }

    return 0;
}

static int check_multi_parent_graph() {
    const std::string path{"/tmp/gffsub_annotation_index_multiparent.gff3"};
    std::ofstream out{path};
    out << "chr1\t.\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=GeneOne\n"
        << "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1\n"
        << "chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=tx2;Parent=gene1\n"
        << "chr1\t.\texon\t10\t20\t.\t+\t.\tID=exon_shared;Parent=tx1,tx2\n";
    out.close();

    const auto index = gffsub::AnnotationIndex::from_gff3(path);
    const auto gene = index.find_gene("GeneOne");
    if (!gene || gene->id != "gene1") {
        std::cerr << "find_gene failed for Name=GeneOne\n";
        return 1;
    }

    const auto parents = index.parents_of("exon_shared");
    if (parents.size() != 2) {
        std::cerr << "multi-parent parents_of failed\n";
        return 1;
    }

    const auto tx1_children = index.children_of("tx1");
    const auto tx2_children = index.children_of("tx2");
    if (tx1_children.size() != 1 || tx2_children.size() != 1) {
        std::cerr << "multi-parent children_of failed\n";
        return 1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: annotation_index_smoke <annotation.gff3>\n";
        return 2;
    }

    if (check_soybean_annotation(argv[1]) != 0) {
        return 1;
    }
    if (check_multi_parent_graph() != 0) {
        return 1;
    }

    std::cout << "annotation_index_smoke OK\n";
    return 0;
}
