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

## Command Shortcuts

These are the shortest forms for the common workflows:

```bash
# Region subset
gffsub annotation.gff3 -r chr1:1-100000

# Region subset, one feature type
gffsub annotation.gff3 -r chr1:1-100000 -f gene

# BED-driven subset
gffsub annotation.gff3 -b regions.bed

# Convert to BED
gffsub annotation.gff3 -r chr1:1-100000 -t bed

# Write output to a file
gffsub annotation.gff3 -r chr1:1-100000 -o subset.gff3

# Longest transcript per gene
gffsub annotation.gff3 -L

# Query one ID
gffsub query annotation.gff3 --id GeneA

# Query one name or gene lookup key
gffsub query annotation.gff3 --name ABC1

# Query a batch of IDs
gffsub query annotation.gff3 --id-list genes.txt

# Query a region
gffsub query annotation.gff3 --region chr1:1-100000

# Query and summarize
gffsub query annotation.gff3 --id GeneA --summary-format tsv

# Extract an upstream/downstream window
gffsub window annotation.gff3 --id GeneA --upstream 2000 --downstream 500

# QC table
gffsub qc annotation.gff3
```

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

## CLI Parameters

`gffsub` has four command surfaces: the classic subset mode, `query`, `window`, and `qc`.

### Classic Subset Mode

```bash
gffsub <input.gff3> [options]
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input GFF3/GTF-style annotation file. |
| `-r`, `--region` | `CHR:START-END` | Keep features overlapping a 1-based inclusive region. |
| `-b`, `--bed` | file | Keep features overlapping BED intervals; BED is read as 0-based half-open. |
| `-f`, `--feature` | type | Keep only records whose third column matches the feature type, such as `gene`, `mRNA`, `transcript`, `exon`, or `CDS`. |
| `-L`, `--longest` | flag | Keep the longest transcript isoform per gene. CDS length is used when CDS isoforms exist; otherwise exon length is used. |
| `-@`, `--threads` | integer | Set worker threads for `--longest`; values above 256 are capped. Set this explicitly when reproducibility of resource use matters. |
| `-t`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | Select output format. `gtf` is accepted as `gtf2`. Default output is `gff3`. |
| `-o`, `--output` | file | Write output to a file instead of stdout. |
| `-h`, `--help` | flag | Show help for classic subset mode. |

Options can be combined. For example, `-r chr1:1-100000 -f gene -t bed` first filters by region, then feature type, then prints BED coordinates.

### Query Mode

```bash
gffsub query <input.gff3> [options]
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input annotation file. |
| `--id` | ID | Query an exact feature `ID`. This option can be repeated. |
| `--name` | key | Gene lookup by `ID`, `Name`, `gene_id`, `locus_tag`, `Alias`, or full `Dbxref` value. |
| `--id-list` | file | Read one exact feature ID per non-empty line. |
| `--region` | `CHR:START-END` | Query features overlapping a 1-based inclusive region. |
| `--type` | type | Restrict query output to one feature type. |
| `--attr` | `KEY=VALUE` | Query exact attribute matches. This option can be repeated. |
| `--include-children` | flag | Include descendants of matched records, such as transcript, exon, CDS, and UTR features. |
| `--summary-format` | `tsv`, `json` | Print summary rows instead of GFF3 records. |
| `-h`, `--help` | flag | Show help for query mode. |

When `--summary-format` is used, the summary columns are `query_id`, `matched_id`, `matched_by`, `seqid`, `start`, `end`, `strand`, `type`, `parent_id`, `child_count`, `transcript_count`, `exon_count`, `cds_length`, and `status`.

### Window Mode

```bash
gffsub window <input.gff3> --id ID [options]
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input annotation file. |
| `--id` | ID or gene key | Required target. The command first tries exact `ID`, then gene lookup. |
| `--upstream` | integer | Bases to add upstream of the target; must be non-negative. Default is `0`. |
| `--downstream` | integer | Bases to add downstream of the target; must be non-negative. Default is `0`. |
| `--strand-aware` | flag | Interpret upstream/downstream biologically by feature strand. Without it, upstream means lower genomic coordinates and downstream means higher genomic coordinates. |
| `-h`, `--help` | flag | Show help for window mode. |

The output is GFF3 records overlapping the expanded window.

### QC Mode

```bash
gffsub qc <input.gff3>
```

`qc` takes only the input file. It writes a TSV table with `severity`, `code`, `line_idx`, `id`, and `message`. Current check codes are `duplicate_id`, `invalid_range`, `missing_parent`, and `child_outside_parent`.

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
