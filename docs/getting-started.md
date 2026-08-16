# Getting Started

gffsub subsets GFF3 and GTF annotation files by region, feature ID, attribute, or gene model.

## Installation

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

Produces `./gffsub`. Needs a C++17 compiler (g++ 9+, clang 10+), no external dependencies.

## Quick Start

Sample data (save as `demo.gff3`):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

### Subset by region

Keep records overlapping a genomic interval.

```bash
./gffsub demo.gff3 -r chr1:200-600
```

### Select by ID with children

Keep one feature and all its descendants.

```bash
./gffsub demo.gff3 --id tx01 -C
```

### Filter by feature type

```bash
./gffsub demo.gff3 -f exon
```

### Convert to GTF

```bash
./gffsub demo.gff3 --id gene01 -C -t gtf
```

### Longest isoform

Keep the longest mRNA per gene.

```bash
./gffsub demo.gff3 --longest
```

## Command Structure

```
gffsub <input.gff3> [options]
gffsub query <input.gff3> [options]
gffsub window <input.gff3> [options]
```

- Default mode: subset and filter records, write annotation output.
- `query`: look up features by ID, name, or attribute, with summary output.
- `window`: records overlapping a window around one feature.

## Help

```bash
./gffsub -h          # full flag list
./gffsub query -h    # query-specific options
./gffsub window -h   # window-specific options
```
