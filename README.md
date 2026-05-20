# gffsub

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

`gffsub` is a small C++17 tool for extracting, querying, summarizing, and checking GFF3/GTF-style genome annotation files. It is built for day-to-day annotation work: subset a region, find a gene by multiple identifiers, pull a gene window, keep the longest isoform, or run quick QC before a downstream pipeline.

## Start With Your Task

| I want to... | Use this |
|-------------|----------|
| Extract features in a genomic interval | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| Use BED regions as input | `gffsub annotation.gff3 -b regions.bed -f exon` |
| Find one gene or feature by ID/name/attribute | `gffsub query annotation.gff3 --id GeneA` |
| Query many IDs and include their children | `gffsub query annotation.gff3 --id-list genes.txt --include-children` |
| Produce a pipeline-friendly summary | `gffsub query annotation.gff3 --id GeneA --summary-format tsv` |
| Extract upstream/downstream context | `gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware` |
| Keep the longest transcript per gene | `gffsub annotation.gff3 --longest` |
| Check annotation graph problems | `gffsub qc annotation.gff3` |

## Install

### Requirements

- C++17 compiler (`g++` or `clang++`)
- CMake is optional; the repository also includes a `Makefile`

### Build With Make

```bash
cd gffsub
make
```

### Build And Test With CMake

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

## Scenario: Extract Features From Regions

Use the classic mode when you already know the interval or have a BED file.

```bash
# Extract genes in a 1-based inclusive GFF region
./gffsub annotation.gff3 -r chr1:1-100000 -f gene

# Extract exons overlapping BED intervals
./gffsub annotation.gff3 -b regions.bed -f exon

# Convert overlapping features to BED
./gffsub annotation.gff3 -r chr1:1-100000 -t bed
```

Coordinate rules are explicit:

| Input or output | Coordinate system |
|-----------------|-------------------|
| GFF3/GTF records | 1-based inclusive |
| `CHR:START-END` regions | 1-based inclusive |
| BED input/output | 0-based half-open |

## Scenario: Find Genes And Gene Models

Use `query` when the question starts from an identifier, name, attribute, or region.

```bash
# Exact feature ID
./gffsub query annotation.gff3 --id Glyma.01G000100

# Gene lookup by ID, Name, gene_id, locus_tag, Alias, or full Dbxref value
./gffsub query annotation.gff3 --name ABC1

# Attribute selector
./gffsub query annotation.gff3 --attr biotype=protein_coding

# Include descendants such as transcript, exon, CDS, and UTR records
./gffsub query annotation.gff3 --id Glyma.01G000100 --include-children
```

For batch pipelines, ask for summaries instead of raw GFF3:

```bash
./gffsub query annotation.gff3 --id-list genes.txt --summary-format tsv
./gffsub query annotation.gff3 --id GeneA --summary-format json
```

Summary fields include query ID, matched ID, matched field, coordinates, strand, feature type, parent ID, child count, transcript count, exon count, CDS length, and status.

## Scenario: Extract Upstream Or Downstream Windows

Use `window` when you want a local annotation context around a gene or feature.

```bash
# Genomic expansion: left/right on the reference sequence
./gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500

# Biological upstream/downstream: strand-aware
./gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500 --strand-aware
```

Without `--strand-aware`, upstream means lower genomic coordinates and downstream means higher genomic coordinates. With `--strand-aware`, the interpretation follows the feature strand.

## Scenario: Keep One Transcript Per Gene

Use `--longest` when a downstream tool expects one representative transcript per gene.

```bash
./gffsub annotation.gff3 --longest
./gffsub annotation.gff3 --longest -@ 6
```

The longest isoform logic follows the existing AGAT-style rule in this project: if a gene has CDS isoforms, compare CDS length; otherwise compare exon length.

## Scenario: Check Annotation Quality

Run `qc` before feeding annotations into a graph-aware workflow.

```bash
./gffsub qc annotation.gff3
```

Current checks report duplicate IDs, invalid ranges, missing parents, and child features outside parent coordinates. Output is TSV so it can be filtered with standard command-line tools.

## Output Formats

| Format | Header | Coordinate |
|--------|--------|------------|
| gff3 | `##gff-version 3` | 1-based inclusive |
| gtf2 | `##gtf-version 2` | 1-based inclusive |
| gtf3 | `##gtf-version 2.2.1` | 1-based inclusive |
| bed | no header | 0-based half-open |

## C++ API

The public C++ API centers on `AnnotationIndex`:

```cpp
auto index = gffsub::AnnotationIndex::from_gff3("annotation.gff3");
auto gene = index.find_gene("GeneA");
auto model = index.gene_model("GeneA");
auto hits = index.overlap("chr1", 1000, 2000);
```

Available queries include `find_by_id`, `find_gene`, `parents_of`, `children_of`, `descendants_of`, `gene_model`, `overlap`, `nearest_gene`, and `with_attribute`.

## Distribution

`gffsub` builds as a single command-line binary. Copy the binary to a glibc-based Linux x86_64 machine and run it with your annotation files.

## License

MIT License
