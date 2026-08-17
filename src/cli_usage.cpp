#include "cli_usage.hpp"
#include "version.hpp"

#include <iostream>

namespace gffsub {

void usage(const char* prog) {
    std::cerr
        << "Program: gffsub (subset GFF3/GTF annotations)\n"
        << "Version: " << kVersion << "\n"
        << "\n"
        << "Usage:   " << prog << " <input> [options]\n"
        << "         " << prog << " query <input> [options]\n"
        << "         " << prog << " window <input> -i ID [options]\n"
        << "\n"
        << "Options:\n"
        << "  -i, --id ID              feature ID (repeatable)\n"
        << "      --ids FILE           one ID per line\n"
        << "  -n, --name NAME          gene by ID/Name/gene_id/locus_tag/Alias/Dbxref\n"
        << "  -w, --where KEY=VALUE    exact column-9 attribute (repeatable)\n"
        << "      --grep FIELD:PATTERN substring match (repeatable)\n"
        << "      --grep-regex FIELD:REGEX  regex match (repeatable)\n"
        << "      --grep-file FILE     patterns from file (one per line)\n"
        << "      --grep-field FIELD   field for --grep-file\n"
        << "      --grep-file-regex    treat --grep-file lines as regex\n"
        << "  -I, --include-expr EXPR  keep matching records (repeatable)\n"
        << "  -E, --exclude-expr EXPR  drop matching records (repeatable)\n"
        << "  -v, --invert-match       invert grep matches\n"
        << "      --ignore-case        case-insensitive grep/expression\n"
        << "  -C, --children           include descendants\n"
        << "  -p, --parents            include ancestors\n"
        << "  -m, --model              include full gene model\n"
        << "  -N, --nearest REGION     nearest gene to CHR:START-END\n"
        << "  -s, --summary            output TSV summary of filtered records\n"
        << "  -u, --up N               upstream window (bp) [0]\n"
        << "  -D, --down N             downstream window (bp) [0]\n"
        << "  -a, --strand-aware       window follows feature strand\n"
        << "  -r, --region REGION      overlap CHR:START-END (1-based)\n"
        << "  -b, --bed FILE           overlap BED intervals (0-based)\n"
        << "  -S, --seqid LIST         keep seqids; ^LIST excludes\n"
        << "      --source SOURCE      exact source column\n"
        << "      --score SCORE        score column (. for missing)\n"
        << "      --strand STRAND      strand column: +|-|.|?\n"
        << "      --phase PHASE        phase column: 0|1|2|.\n"
        << "  -f, --feature TYPE       type column (repeatable)\n"
        << "  -L, --longest            one transcript per gene\n"
        << "  -@, --threads N          threads for --longest [1]\n"
        << "  -t, --format FMT         output: gff3|gtf|gtf2|gtf3|bed [gff3]\n"
        << "  -o, --output FILE        write to FILE [stdout]\n"
        << "  -h, --help               show this help\n"
        << "      --version            print version and exit\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " ann.gff3 -r chr1:1-100000 -f gene\n"
        << "  " << prog << " ann.gff3 -i GeneA -C\n"
        << "  " << prog << " ann.gff3 -i GeneA -m -t gtf\n"
        << "  " << prog << " ann.gff3 --longest -@ 4\n"
        << "  " << prog << " query ann.gff3 -i GeneA -s tsv\n"
        << "  " << prog << " window ann.gff3 -i GeneA -u 2000 -D 500 -a\n";
}

void query_usage(const char* prog) {
    std::cerr
        << "About:   query features by ID, name, region, or attributes\n"
        << "Usage: " << prog << " query <input> [options]\n"
        << "\n"
        << "Options:\n"
        << "  -i, --id ID              feature ID\n"
        << "      --ids FILE           one ID per line\n"
        << "  -n, --name NAME          gene naming keys\n"
        << "  -r, --region REGION      overlap (1-based)\n"
        << "      --type TYPE          restrict by feature type\n"
        << "  -w, --where KEY=VALUE    exact attribute\n"
        << "  -C, --children           include descendants\n"
        << "  -p, --parents            include ancestors\n"
        << "  -m, --model              full gene model\n"
        << "  -N, --nearest REGION     nearest gene\n"
        << "  -s, --summary            output TSV summary\n"
        << "  -h, --help               show this help\n"
        << "\n"
        << "Example:\n"
        << "  " << prog << " query ann.gff3 -i GeneA -s\n";
}

void window_usage(const char* prog) {
    std::cerr
        << "About:   records overlapping a window around one feature\n"
        << "Usage: " << prog << " window <input> -i ID [options]\n"
        << "\n"
        << "Options:\n"
        << "  -i, --id ID              target feature ID or gene name\n"
        << "  -u, --up N               bases upstream [0]\n"
        << "  -D, --down N             bases downstream [0]\n"
        << "  -a, --strand-aware       window follows feature strand\n"
        << "  -h, --help               show this help\n"
        << "\n"
        << "Example:\n"
        << "  " << prog << " window ann.gff3 -i GeneA -u 2000 -D 500\n";
}

}  // namespace gffsub
