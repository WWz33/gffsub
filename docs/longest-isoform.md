# Longest Isoform Selection

<!-- I18N:START -->

**English** | [中文](./longest-isoform.zh.md)

<!-- I18N:END -->

## -L / --longest

Keep one transcript per gene. Length metric:

- CDS length (sum of CDS segments) if the gene has CDS, otherwise exon length (sum of exon segments).
- Discontinuous CDS (same ID, multiple lines): segments summed.
- Alternative CDS variants (distinct IDs under one transcript): longest variant used, not summed.
- Ties keep the first encountered.
- Auto-detects isoform type: the most frequent transcript-class type in the file (`mRNA`, `transcript`, `ncRNA`, `tRNA`, ...). Ties prefer `mRNA`, then `transcript`.
- Genes with one isoform are left unchanged.

Sample data (demo.gff3):

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

- tx01 CDS length: 151 + 251 = 402 bp.
- tx02 CDS length: 301 bp.
- tx03 exon length: 201 + 151 = 352 bp (no CDS).

```bash
# keep tx01 for gene01 (402 bp), tx03 for gene02 (only isoform)
./gffsub demo.gff3 --longest

# keep only the longest transcript rows (drop gene/exon)
./gffsub demo.gff3 --longest -t mRNA
```

## --longest-type TYPE

Isoform type for `--longest` selection. Independent of `-t`: `--longest-type` picks which records compete as isoforms, `-t` filters output rows. Requires `--longest`.

```bash
# ncRNA file: select longest ncRNA per gene, keep hierarchy
./gffsub lnc.gff3 --longest --longest-type ncRNA

# same, but output only the ncRNA rows
./gffsub lnc.gff3 --longest --longest-type ncRNA -t ncRNA
```

## -@ / --threads

Parallelize by chromosome. Default 1, max 256.

```bash
./gffsub demo.gff3 --longest --threads 4
```
