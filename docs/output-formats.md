# Output Formats

<!-- I18N:START -->

**English** | [中文](./output-formats.zh.md)

<!-- I18N:END -->

## -t / --format

Syntax: `-t FMT` or `--format FMT` (alias `--output-format FMT`).

| Value   | Output |
|---------|--------|
| `gff3`  | GFF3, original column 9 (default) |
| `gtf`   | GTF2.2, gene_id + transcript_id |
| `gtf2`  | alias for `gtf` |
| `gtf3`  | GTF2.2.1, mRNA renamed to transcript, type whitelist |
| `bed`   | 6-column BED: chrom, start, end, name, score, strand |

## GFF3

- Header: `##gff-version 3`
- Column 9 preserved from input.
- GTF input: column 9 rewritten as `tag=value` with `ID=`/`Parent=` and URL escaping.

## GTF

- Header: `##gtf-version 2`
- `gene_id` on every line.
- `transcript_id` on non-gene lines only.
- `mRNA` stays as `mRNA` in `gtf`/`gtf2`.
- Values escaped: backslash, quote, tab, newline, CR.
- When `gene_id` cannot be resolved, an empty value is emitted (`gene_id "";`).

## GTF3

- Header: `##gtf-version 2.2.1`
- Type whitelist: `gene`, `transcript`, `exon`, `CDS`, `start_codon`, `stop_codon`, `five_prime_utr`, `three_prime_utr`, `Selenocysteine`, `mRNA`, `5UTR`, `3UTR`, `inter`, `inter_CNS`, `intron_CNS`.
- `mRNA` renamed to `transcript`.
- UTR normalized to `five_prime_utr`/`three_prime_utr`.

## BED

- 6 columns: chrom, start (0-based), end, name, score, strand.
- `start` = record start - 1.
- `score` from column 6, or `0` when missing or `.`.
- `name` from `ID`, or the type when `ID` is missing.

## Example

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

```bash
./gffsub demo.gff3 --format gtf
./gffsub demo.gff3 --format gtf3
./gffsub demo.gff3 --format bed
./gffsub input.gtf --format gff3
./gffsub demo.gff3 --format gtf -o output.gtf
```
