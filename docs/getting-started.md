# Getting Started

<!-- I18N:START -->

**English** | [中文](./getting-started.zh.md)

<!-- I18N:END -->

gffsub subsets GFF3 and GTF annotation files by region, feature ID, attribute, or gene model.

## Build

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

Produces `./gffsub`. Requires a C++17 compiler (g++ 9+, clang 10+) and zlib (`-lz`; ships with macOS, install `zlib1g-dev` on Debian/Ubuntu).

## Inputs

- GFF3, GTF, or BED, detected from content (not the extension)
- Plain or gzip-compressed files read directly (`ann.gtf.gz`)
- `-` reads stdin, plain or gzip-compressed

## Example data

Save as `demo.gff3`:

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## Examples

Subset by region:

```bash
./gffsub demo.gff3 -r chr1:200-600
```

Select by ID with children:

```bash
./gffsub demo.gff3 -i tx01 -C
```

Filter by type:

```bash
./gffsub demo.gff3 -t exon
```

Convert to GTF:

```bash
./gffsub demo.gff3 --format gtf
```

Longest isoform per gene:

```bash
./gffsub demo.gff3 --longest
```

## Longer commands

```bash
./gffsub ann.gff3 -w biotype=protein_coding -L -m -f gtf -o coding.gtf
```

```bash
./gffsub ann.gff3 -S chr1 --exclude-region chr1:500000-800000 --drop-orphans
```

```bash
./gffsub ann.gff3 --ids genes.txt -m --out-attrs Name,biotype
```

```bash
./gffsub ann.gff3 -N chr1:1000000-1000500 -m
```

```bash
./gffsub ann.gff3.gz -k seqid,natural-seqid,start -t gene -f bed
```

```bash
./gffsub ann.gff3 -w biotype=protein_coding -C | ./gffsub - -L -t mRNA
```

## Command structure

```
gffsub <input> [options]
gffsub query <input> [options]
gffsub window <input> [options]
```

## Help

```bash
./gffsub -h
./gffsub query -h
./gffsub window -h
```
