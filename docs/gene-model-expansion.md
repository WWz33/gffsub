# Gene Model Expansion

gffsub builds a parent/child index from the ID and Parent attributes. Three flags expand selector hits along the hierarchy.

## Hierarchy

```
gene
  └── mRNA / transcript
        ├── exon
        └── CDS
```

The index follows Parent= links. GTF input uses synthesized parent/child links from gene_id and transcript_id.

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

## --children / -C: include descendants

Recursively includes all children of the selector hit.

```bash
# gene01 plus all its mRNAs, exons, CDS
./gffsub demo.gff3 --id gene01 -C

# tx01 plus its exons and CDS
./gffsub demo.gff3 --id tx01 -C
```
--id tx01 -C returns: tx01, ex01, cds01, ex02, cds02.
--id gene01 -C returns: gene01, tx01, ex01, cds01, ex02, cds02, tx02, ex03.

## --parents: include ancestors

Walks up the Parent chain to the root (gene).

```bash
# exon and its parent mRNA and grandparent gene
./gffsub demo.gff3 --id ex01 --parents
```
Returns: ex01, tx01, gene01.

## --model: full gene model

Includes the complete gene model: the gene, the matching selector hit, all siblings (other transcripts of the same gene), and all their children.

```bash
# full model containing tx01
./gffsub demo.gff3 --id tx01 --model
```
Returns: gene01, tx01, ex01, cds01, ex02, cds02, tx02, ex03.

## Comparison

| Flag | Direction | Includes |
|------|-----------|----------|
| `-C` / `--children` | downward | selector hit + all descendants |
| `--parents` | upward | selector hit + all ancestors to gene |
| `--model` | both | full gene model (gene + all transcripts + all children) |

## With multiple selectors

When combined with --name or multiple --id:
- Each selector hit is expanded independently
- Results are deduplicated by line index

```bash
# full models for two genes
./gffsub demo.gff3 --id gene01 --id gene02 --model
```

## With region filter

Region, seqid, and feature-type filters apply first. Gene model expansion happens after filtering, so it can pull in records that the region filter excluded.
