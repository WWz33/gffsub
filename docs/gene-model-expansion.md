# Gene Model Expansion

<!-- I18N:START -->

**English** | [中文](./gene-model-expansion.zh.md)

<!-- I18N:END -->

gffsub follows `Parent=` links to expand selector hits. GTF uses synthesized links from `gene_id` / `transcript_id`.

Sample data (demo.gff3):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex03;Parent=tx02
chr2	src	gene	200	600	.	-	.	ID=gene02;Name=XYZ1
chr2	src	mRNA	200	600	.	-	.	ID=tx03;Parent=gene02
chr2	src	exon	200	400	.	-	.	ID=ex04;Parent=tx03
```

## -C / --children

Include the selector hit and all its descendants.

```bash
./gffsub demo.gff3 -i gene01 -C
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex03;Parent=tx02
```

Selecting a transcript stops below it:

```bash
./gffsub demo.gff3 -i tx01 -C
```

```
##gff-version 3
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## -p / --parents

Include the selector hit and all its ancestors up to the gene.

```bash
./gffsub demo.gff3 -i ex01 -p
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
```

## -m / --model

Include the full gene model: the gene, the selector hit, all sibling transcripts, and all their children.

```bash
./gffsub demo.gff3 -i tx01 -m
```

The output is the same as `-i gene01 -C` above, because tx01 belongs to gene01 and gene01 has no other children.

## Comparison

| Flag | Direction | Includes |
|------|-----------|----------|
| `-C` / `--children` | downward | selector hit + all descendants |
| `-p` | upward | selector hit + all ancestors to gene |
| `-m` | both | full gene model |

## --drop-orphans

Remove records all of whose Parent references point outside the kept set. Applied after every other filter and repeated until a fixpoint, so dropping a transcript also drops its exons. Records without parents are never dropped.

```bash
./gffsub demo.gff3 -t exon --drop-orphans
```

```
##gff-version 3
```

`-t exon` removes the transcript rows, so every exon loses its parent and the fixpoint drops the whole set.
