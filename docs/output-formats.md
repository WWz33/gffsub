# Output Formats

<!-- I18N:START -->

**English** | [中文](./output-formats.zh.md)

<!-- I18N:END -->

## -t / --format flag

Syntax: `-t FMT` or `--format FMT` (alias: `--output-format FMT`)

| Format | Value | Description |
|--------|-------|-------------|
| GFF3 | `gff3` | GFF3 with original column 9 attributes (default) |
| GTF | `gtf` | GTF2.2 with gene_id and transcript_id |
| GTF2 | `gtf2` | alias for gtf |
| GTF3 | `gtf3` | GTF2.2.1 with type renaming (mRNA to transcript) |
| BED | `bed` | 6-column BED (chrom, start, end, name, score, strand) |

## GFF3 output

- Header: `##gff-version 3`
- Column 9 preserved as-is from input
- Records output in input order

## GTF output

- Header: `##gtf-version 2`
- Every feature line has gene_id (required by GTF2.2)
- transcript_id on non-gene features only (GTF2.2: gene_id on all lines, transcript_id on all except gene)
- gene_id resolved from Parent chain: child features (exon, CDS) look up their mRNA's gene via an internal mRNA-to-gene map
- When gene_id cannot be resolved, an empty value is emitted (`gene_id "";`) rather than dropping the feature
- mRNA type stays as "mRNA" in gtf/gtf2, renamed to "transcript" in gtf3
- Attribute values are escaped (backslash and quote)

## GTF3 output

- Header: `##gtf-version 2.2.1`
- Only GTF3-compatible feature types pass through: gene, transcript, exon, CDS, start_codon, stop_codon, five_prime_utr, three_prime_utr, Selenocysteine, mRNA
- mRNA renamed to transcript

## BED output

- 6 columns: chrom, start (0-based), end, name, score, strand
- start = record start - 1 (GFF3 is 1-based, BED is 0-based)
- score from column 6, or `0` when missing or `.`
- name from ID attribute, or the feature type when ID is missing

## Format conversion examples

Sample input (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
```

GFF3 (default):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
...
```

GTF:
```
##gtf-version 2
chr1	src	gene	100	1000	.	+	.	gene_id "gene01";
chr1	src	mRNA	100	1000	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

BED:
```
chr1	99	1000	gene01	0	+
chr1	99	1000	tx01	0	+
chr1	99	250	ex01	0	+
chr1	99	250	cds01	0	+
```

## Commands

```bash
# GFF3 to GTF
./gffsub demo.gff3 -t gtf

# GFF3 to GTF3 (rename mRNA to transcript)
./gffsub demo.gff3 -t gtf3

# GFF3 to BED
./gffsub demo.gff3 -t bed

# GTF to GFF3
./gffsub input.gtf -t gff3

# output to file
./gffsub demo.gff3 -t gtf -o output.gtf
```
