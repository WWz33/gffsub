# gffsub

<!-- README-I18N:START -->

**English** | [汉语](./README.zh.md)

<!-- README-I18N:END -->

`gffsub` is a small C++17 tool for extracting, querying, summarizing, and checking GFF3/GTF-style genome annotation files. It is built for day-to-day annotation work: subset a region, find a gene by multiple identifiers, pull a gene window, keep the longest isoform, or run quick QC before a downstream pipeline.

## Start With Your Task

| I want to... | Use this |
|-------------|----------|
| Extract features in a genomic interval | `gffsub annotation.gff3 -r chr1:1-100000 -f gene` |
| Extract one chromosome or contig | `gffsub annotation.gff3 --seqid chr1` |
| Extract records from one annotation source | `gffsub annotation.gff3 --source Gnomon` |
| Extract records with one score | `gffsub annotation.gff3 --score 42.5` |
| Extract records on one strand | `gffsub annotation.gff3 --strand -` |
| Extract records with one phase | `gffsub annotation.gff3 --phase 0 -f CDS` |
| Use BED regions as input | `gffsub annotation.gff3 -b regions.bed -f exon` |
| Find one feature by exact ID | `gffsub annotation.gff3 --id GeneA` |
| Extract many exact IDs | `gffsub annotation.gff3 --ids genes.txt` |
| Find a gene or feature by name | `gffsub annotation.gff3 --name GeneA` |
| Find features by attribute value | `gffsub annotation.gff3 --where biotype=protein_coding` |
| Find the nearest gene to a region | `gffsub annotation.gff3 --nearest chr1:1000-2000` |
| Extract many IDs and include their children | `gffsub annotation.gff3 --ids genes.txt -C` |
| Extract a feature and its parents | `gffsub annotation.gff3 --id ExonA --parents` |
| Extract the full gene model from any feature | `gffsub annotation.gff3 --id ExonA --model` |
| Produce a pipeline-friendly summary | `gffsub annotation.gff3 --id GeneA --summary tsv` |
| Extract selected attribute values | `gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent` |
| Extract upstream/downstream context | `gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware` |
| Keep the longest transcript per gene | `gffsub annotation.gff3 --longest` |
| Check annotation graph problems | `gffsub annotation.gff3 --qc` |

## Command Shortcuts

These are the shortest forms for the common workflows:

```bash
# Region subset
gffsub annotation.gff3 -r chr1:1-100000

# Whole seqid subset
gffsub annotation.gff3 --seqid chr1

# Source column subset
gffsub annotation.gff3 --source Gnomon

# Score column subset
gffsub annotation.gff3 --score 42.5

# Strand column subset
gffsub annotation.gff3 --strand -

# Phase column subset
gffsub annotation.gff3 --phase 0

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

# Extract one exact ID
gffsub annotation.gff3 --id GeneA

# Extract one name or gene lookup key
gffsub annotation.gff3 --name ABC1

# Extract exact attribute matches
gffsub annotation.gff3 --where biotype=protein_coding

# Find the nearest gene to a 1-based inclusive region
gffsub annotation.gff3 --nearest chr1:1000-2000

# Extract a batch of IDs
gffsub annotation.gff3 --ids genes.txt

# Extract an ID and its descendants
gffsub annotation.gff3 --id GeneA -C

# Extract an ID and its ancestors
gffsub annotation.gff3 --id ExonA --parents

# Extract the full gene model containing an ID
gffsub annotation.gff3 --id ExonA --model

# Extract a region
gffsub annotation.gff3 --region chr1:1-100000

# Query and summarize
gffsub annotation.gff3 --id GeneA --summary tsv

# Query and extract selected attributes
gffsub annotation.gff3 --id GeneA --out-attrs ID,Name,Parent

# Extract an upstream/downstream window
gffsub annotation.gff3 --id GeneA --up 2000 --down 500

# QC table
gffsub annotation.gff3 --qc
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

# Extract all records from one GFF3 seqid column value
./gffsub annotation.gff3 --seqid chr1

# Extract all records from one GFF3 source column value
./gffsub annotation.gff3 --source Gnomon

# Extract all records from one GFF3 score column value
./gffsub annotation.gff3 --score 42.5

# Extract all records from one GFF3 strand column value
./gffsub annotation.gff3 --strand -

# Extract all records from one GFF3 phase column value
./gffsub annotation.gff3 --phase 0

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

Use the top-level selector form when the question starts from an identifier, name, attribute, or region.

```bash
# Exact feature ID
./gffsub annotation.gff3 --id Glyma.01G000100

# Gene lookup by ID, Name, gene_id, locus_tag, Alias, or full Dbxref value
./gffsub annotation.gff3 --name ABC1

# Attribute selector
./gffsub annotation.gff3 --where biotype=protein_coding

# Nearest gene to a 1-based inclusive region
./gffsub annotation.gff3 --nearest chr1:1000-2000

# Include descendants such as transcript, exon, CDS, and UTR records
./gffsub annotation.gff3 --id Glyma.01G000100 -C

# Include ancestors linked through Parent
./gffsub annotation.gff3 --id ExonA --parents

# Return the full gene model containing this feature
./gffsub annotation.gff3 --id ExonA --model
```

For batch pipelines, ask for summaries instead of raw GFF3:

```bash
./gffsub annotation.gff3 --ids genes.txt --summary tsv
./gffsub annotation.gff3 --id GeneA --summary json
```

Summary fields include query ID, matched ID, matched field, coordinates, strand, feature type, parent ID, child count, transcript count, exon count, CDS length, and status.

### Attribute Keys In Search And Output

GFF3 stores record attributes in column 9 as semicolon-separated `KEY=VALUE` pairs:

```gff3
chr1	src	gene	100	400	.	+	.	ID=gene0001;Name=ABC1;Alias=ABC-1;Dbxref=GeneID:123
```

Use `--id` for exact `ID` lookup, `--name` for gene lookup across common naming keys, and `--where KEY=VALUE` for any exact attribute-value filter:

```bash
./gffsub annotation.gff3 --id gene0001
./gffsub annotation.gff3 --name ABC1
./gffsub annotation.gff3 --where Alias=ABC-1
./gffsub annotation.gff3 --where Dbxref=GeneID:123
```

Non-region selectors use these column-9 keys:

| Task | Command | Keys used |
|------|---------|-----------|
| Exact feature lookup | `--id gene0001` | `ID` |
| Batch exact feature lookup | `--ids genes.txt` | `ID` values, one per line |
| Gene lookup | `--name ABC1` | gene records by `ID`, `gene_id`, `Name`, `locus_tag`, `Alias`, or full `Dbxref` value |
| Any exact attribute filter | `--where Parent=gene0001` | any column-9 `KEY=VALUE`, including `ID`, `Name`, `Alias`, `Parent`, `Dbxref`, `Accession`, or `Parent_Accession` |
| Nearest gene lookup | `--nearest chr1:1000-2000` | same-seqid gene with the shortest distance to a 1-based inclusive region; overlapping genes have distance 0; ties use input file order |
| Include matched descendants | `-C`, `--children` | child records linked by `Parent`; starts from records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`; `--include-children` is a verbose alias |
| Include matched ancestors | `--parents` | parent records reached by walking `Parent` links upward; starts from records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`; `--include-parents` is a verbose alias |
| Extract full gene model | `--model`, `--gene-model` | starts from records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`, then returns the containing gene and its descendants |
| Print selected attributes | `--out-attrs ID,Name,Parent` | prints selected column-9 keys after records are matched |

Use `--out-attrs` when the records are already selected and you want selected column-9 attributes added to TSV/JSON summary output:

```bash
./gffsub annotation.gff3 --id gene0001 --out-attrs ID,Name,Alias,Dbxref
```

`--attr KEY=VALUE` is a compatibility alias for `--where KEY=VALUE`. `--output-attrs` is a verbose alias for `--out-attrs`. `--attrs` remains as a deprecated compatibility alias.

## Scenario: Extract Upstream Or Downstream Windows

Use the top-level window options when you want a local annotation context around a gene or feature.

```bash
# Genomic expansion: left/right on the reference sequence
./gffsub annotation.gff3 --id GeneA --upstream 2000 --downstream 500
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500

# Biological upstream/downstream: strand-aware
./gffsub annotation.gff3 --id GeneA --up 2000 --down 500 --strand-aware
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

Run `--qc` before feeding annotations into a graph-aware workflow.

```bash
./gffsub annotation.gff3 --qc
```

Current checks report invalid `##gff-version` directives, invalid column counts, invalid attribute syntax, duplicate attribute tags, invalid comma-separated attribute values, invalid seqids, invalid feature types, duplicate IDs, invalid or duplicate `##sequence-region` directives, non-positive coordinates, invalid ranges, invalid score values, features outside `##sequence-region` bounds, invalid strand or phase column values, invalid CDS phases, duplicate or missing parents, and child features outside parent coordinates. Output is TSV so it can be filtered with standard command-line tools.

## CLI Parameters

`gffsub` is top-level first: common GFF3 work starts as `gffsub <input.gff3> [options]`. The `query`, `window`, and `qc` subcommands remain as compatible advanced entry points and share the same output semantics.

### Top-Level Mode

```bash
gffsub <input.gff3> [options]
```

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input GFF3/GTF-style annotation file. |
| `--id` | ID | Keep the exact feature `ID`. This option can be repeated. This is the top-level selector for exact-ID extraction. |
| `--ids`, `--id-list` | file | Read one exact feature ID per non-empty line. This is the top-level selector for batch exact-ID extraction. `--id-list` is a verbose alias. |
| `--name` | key | Keep one gene found by `ID`, `Name`, `gene_id`, `locus_tag`, `Alias`, or full `Dbxref` value. This is the top-level selector for common gene naming keys. |
| `--where`, `--attr` | `KEY=VALUE` | Keep features with an exact GFF3 attribute value. This option can be repeated. This is the top-level selector for exact column-9 `KEY=VALUE` matches. `--attr` is a compatibility alias. |
| `-C`, `--children`, `--include-children` | flag | Include descendants of records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`. `--include-children` is a verbose alias. |
| `--parents`, `--include-parents` | flag | Include ancestors of records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`. `--include-parents` is a verbose alias. |
| `--model`, `--gene-model` | flag | Include the full gene model containing records matched by `--id`, `--ids`, `--name`, `--where`, or `--nearest`; this returns the gene plus transcript/exon/CDS/UTR descendants. `--gene-model` is a verbose alias. |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | Keep the nearest gene on the same seqid as a 1-based inclusive region; ties use input file order. `--nearest-gene` is a verbose alias. |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | Print selected column-9 attributes as extra TSV/JSON fields. Combine only with query-style selectors. |
| `--attrs` | `KEY1,KEY2,...` | Deprecated compatibility alias for `--out-attrs`. |
| `--summary`, `--summary-format` | `tsv`, `json` | Print summary rows instead of GFF3 records. Combine only with query-style selectors. `--summary-format` is a verbose alias. |
| `--up`, `--upstream` | integer | With `--id`, extract records overlapping the upstream-expanded target window. |
| `--down`, `--downstream` | integer | With `--id`, extract records overlapping the downstream-expanded target window. |
| `--strand-aware` | flag | With window extraction, interpret upstream/downstream biologically by feature strand. |
| `--qc` | flag | Run annotation QC. |
| `--seqid` | seqid | Keep records whose first GFF3 column exactly matches the value. |
| `--source` | source | Keep records whose second GFF3 column exactly matches the value. |
| `--score` | number, `.` | Keep records whose sixth GFF3 column matches the numeric score, or whose score column is `.`. |
| `--strand` | `+`, `-`, `.`, `?` | Keep records whose seventh GFF3 column exactly matches the value. This is a filter, unlike `--strand-aware`, which only changes window interpretation. |
| `--phase` | `0`, `1`, `2`, `.` | Keep records whose eighth GFF3 column exactly matches the value. For CDS records, GFF3 phase is normally `0`, `1`, or `2`; `.` matches any record whose phase column is `.`. |
| `-r`, `--region` | `CHR:START-END` | Keep features overlapping a 1-based inclusive region. |
| `-b`, `--bed` | file | Keep features overlapping BED intervals; BED is read as 0-based half-open. |
| `-f`, `--feature`, `--type` | type | Keep only records whose third column matches the feature type, such as `gene`, `mRNA`, `transcript`, `exon`, or `CDS`. |
| `-L`, `--longest` | flag | Keep the longest transcript isoform per gene. CDS length is used when CDS isoforms exist; otherwise exon length is used. |
| `-@`, `--threads` | integer | Set worker threads for `--longest`; values above 256 are capped. Set this explicitly when reproducibility of resource use matters. |
| `-t`, `--format`, `--output-format` | `gff3`, `gtf`, `gtf2`, `gtf3`, `bed` | Select output format. `gtf` is accepted as `gtf2`. Default output is `gff3`. `--output-format` is a verbose alias. |
| `-o`, `--output` | file | Write output to a file instead of stdout. |
| `-h`, `--help` | flag | Show help for top-level mode. |

Options can be combined. For example, `-r chr1:1-100000 -f gene -t bed` first filters by region, then feature type, then prints BED coordinates.

### Query Compatibility Mode

```bash
gffsub query <input.gff3> [options]
```

Most default GFF3 and summary workflows can be written without the `query` subcommand. This compatibility mode is kept for existing scripts and explicit query-style command lines.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input annotation file. |
| `--id` | ID | Query an exact feature `ID`. This option can be repeated. |
| `--name` | key | Gene lookup by `ID`, `Name`, `gene_id`, `locus_tag`, `Alias`, or full `Dbxref` value. |
| `--ids`, `--id-list` | file | Read one exact feature ID per non-empty line. `--id-list` is a verbose alias. |
| `--region` | `CHR:START-END` | Query features overlapping a 1-based inclusive region. |
| `--nearest`, `--nearest-gene` | `CHR:START-END` | Query the nearest gene on the same seqid as a 1-based inclusive region; ties use input file order. `--nearest-gene` is a verbose alias. |
| `--type` | type | Restrict query output to one feature type. |
| `--where`, `--attr` | `KEY=VALUE` | Query exact attribute matches. This option can be repeated. `--attr` is a compatibility alias. |
| `--out-attrs`, `--output-attrs` | `KEY1,KEY2,...` | Append selected attribute values as extra TSV/JSON fields. |
| `--attrs` | `KEY1,KEY2,...` | Deprecated compatibility alias for `--out-attrs`. |
| `-C`, `--children`, `--include-children` | flag | Include descendants of matched records, such as transcript, exon, CDS, and UTR features. `--include-children` is a verbose alias. |
| `--parents`, `--include-parents` | flag | Include ancestors of matched records by walking GFF3 `Parent` links upward. `--include-parents` is a verbose alias. |
| `--model`, `--gene-model` | flag | Include the full gene model containing matched records. `--gene-model` is a verbose alias. |
| `--summary`, `--summary-format` | `tsv`, `json` | Print summary rows instead of GFF3 records. `--summary-format` is a verbose alias. |
| `-h`, `--help` | flag | Show help for query mode. |

When `--summary` is used, the summary columns are `query_id`, `matched_id`, `matched_by`, `seqid`, `start`, `end`, `strand`, `type`, `parent_id`, `child_count`, `transcript_count`, `exon_count`, `cds_length`, and `status`. If `--out-attrs` is present, those keys are appended as extra TSV columns or as an `attrs` object in JSON. Without `--summary`, `--out-attrs` prints TSV.

### Window Compatibility Mode

```bash
gffsub window <input.gff3> --id ID [options]
```

The same workflow can be written at top level with `gffsub <input.gff3> --id ID --up N --down N`.

| Parameter | Value | Meaning |
|-----------|-------|---------|
| `<input.gff3>` | file | Input annotation file. |
| `--id` | ID or gene key | Required target. The command first tries exact `ID`, then gene lookup. |
| `--up`, `--upstream` | integer | Bases to add upstream of the target; must be non-negative. Default is `0`. |
| `--down`, `--downstream` | integer | Bases to add downstream of the target; must be non-negative. Default is `0`. |
| `--strand-aware` | flag | Interpret upstream/downstream biologically by feature strand. Without it, upstream means lower genomic coordinates and downstream means higher genomic coordinates. |
| `-h`, `--help` | flag | Show help for window mode. |

The output is GFF3 records overlapping the expanded window.

### QC Compatibility Mode

```bash
gffsub qc <input.gff3>
```

The same workflow can be written at top level with `gffsub <input.gff3> --qc`. QC writes a TSV table with `severity`, `code`, `line_idx`, `id`, and `message`. Current check codes are `invalid_gff_version`, `invalid_column_count`, `invalid_attribute_syntax`, `duplicate_attribute_tag`, `invalid_attribute_multivalue`, `invalid_seqid`, `invalid_feature_type`, `duplicate_id`, `invalid_sequence_region`, `duplicate_sequence_region`, `invalid_coordinate`, `invalid_range`, `invalid_score`, `outside_sequence_region`, `invalid_strand`, `invalid_phase`, `invalid_cds_phase`, `invalid_gap`, `invalid_is_circular`, `invalid_target`, `duplicate_parent`, `parent_cycle`, `missing_derives_from`, `missing_parent`, and `child_outside_parent`. In strict GFF3 QC, comma-separated values are accepted only for `Parent`, `Alias`, `Note`, `Dbxref`, and `Ontology_term`; literal commas in other values should be URL-escaped. A `region` feature marked `Is_circular=true` may make features on that seqid wrap past the end of their `##sequence-region`.

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
