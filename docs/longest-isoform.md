# Longest Isoform Selection

<!-- I18N:START -->

**English** | [中文](./longest-isoform.zh.md)

<!-- I18N:END -->

## Overview

The `-L` / `--longest` flag keeps one isoform per gene. Genes with at least two isoforms of the chosen type are reduced to the longest. The longest is scored by CDS span when the gene has any CDS, otherwise by exon span. The non-selected isoforms and their child features (exons, CDS) are dropped. Genes with one isoform are left unchanged. Genes with no isoform of the chosen type are dropped.

## Algorithm

1. Group isoform-type records (mRNA or transcript) by their parent gene.
2. Group genes by chromosome.
3. For each gene:
   - If the gene has no isoform children of the chosen type, drop the gene.
   - If the gene has exactly one isoform, keep it and its children unchanged.
   - Otherwise compute a per-gene flag: does any isoform of this gene have a CDS child.
   - For each isoform, compute a span:
     - When the gene has CDS: CDS length. Distinct CDS IDs under one transcript are alternative protein variants; the isoform is scored by its LONGEST variant (lines sharing a CDS ID are one discontinuous CDS and are summed). Isoforms without CDS children are skipped, so they cannot be selected.
     - When the gene has no CDS: sum of exon lengths. Isoforms without exon children are skipped, so they cannot be selected.
   - Keep the isoform with the maximum span. Ties keep the first encountered (strict greater-than comparison).
   - Drop non-kept isoforms and their children. Children of the kept isoform stay.
4. Threads (`-@` / `--threads`) parallelize across chromosomes. Each chromosome is processed independently.

## Isoform type auto-detection

When `-f` / `--feature` is not specified:
- If the data has `mRNA` records, use `mRNA`.
- Otherwise, if the data has `transcript` records, use `transcript`.
- If neither exists, default to `mRNA`. All genes are then dropped because no isoform matches.

Covers both GFF3 (typically mRNA) and GTF (typically transcript).

When `-f TYPE` is specified, that type is the isoform type, and the run additionally drops every record whose type is not TYPE (genes, CDS, exons, and others). Longest selection runs first, then the feature filter keeps only TYPE records.

## Sample data

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	CDS	100	400	.	+	0	ID=cds03;Parent=tx02
chr2	src	gene	200	600	.	-	.	ID=gene02
chr2	src	mRNA	200	600	.	-	.	ID=tx03;Parent=gene02
chr2	src	exon	200	400	.	-	.	ID=ex04;Parent=tx03
chr2	src	exon	450	600	.	-	.	ID=ex05;Parent=tx03
```

tx01 CDS span: 151 + 251 = 402 bp.
tx02 CDS span: 301 bp.
tx03 exon span: 201 + 151 = 352 bp (no CDS).

## Commands

```bash
# keep one isoform per gene (tx01 for gene01, tx03 for gene02)
./gffsub demo.gff3 --longest

# use 4 threads
./gffsub demo.gff3 --longest -@ 4

# restrict to a specific isoform type
./gffsub demo.gff3 --longest -f transcript
```

## Result for the sample

- gene01: tx01 kept (CDS span 402), tx02 and cds03 dropped.
- gene02: tx03 kept and unchanged (only one isoform).
- Gene records are kept when they have at least one isoform of the chosen type.

## Combining with other filters

Region, feature type, and attribute filters apply before longest selection. To longest-select within a region:

```bash
./gffsub demo.gff3 -r chr1:1-1000 --longest
```
