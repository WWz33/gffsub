#include "gff3.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 2) {
        std::cerr << "usage: annotation_index_smoke <annotation.gff3>\n";
        return 2;
    }

    const std::string gene_id{"SoyL04_01G000000"};
    const auto index = gffsub::AnnotationIndex::from_gff3(argv[1]);

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

    std::cout << "annotation_index_smoke OK\n";
    return 0;
}
