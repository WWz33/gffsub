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
```

## -C / --children

Include the selector hit and all its descendants.

```bash
# gene01 + all mRNAs, exons, CDS
./gffsub demo.gff3 -i gene01 -C

# tx01 + its exons and CDS
./gffsub demo.gff3 -i tx01 -C
```

`-i gene01 -C` returns: gene01, tx01, ex01, cds01, ex02, cds02, tx02, ex03.

`-i tx01 -C` returns: tx01, ex01, cds01, ex02, cds02.

## -p / -p

Include the selector hit and all its ancestors up to the gene.

```bash
./gffsub demo.gff3 -i ex01 -p
```

Returns: ex01, tx01, gene01.

## -m / -m

Include the full gene model: the gene, the selector hit, all sibling transcripts, and all their children.

```bash
./gffsub demo.gff3 -i tx01 -m
```

Returns: gene01, tx01, ex01, cds01, ex02, cds02, tx02, ex03.

## Comparison

| Flag | Direction | Includes |
|------|-----------|----------|
| `-C` / `-C` | downward | selector hit + all descendants |
| `-p` | upward | selector hit + all ancestors to gene |
| `-m` | both | full gene model |
