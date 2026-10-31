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

- GFF3 or GTF, detected from content (not the extension)
- BED files are read for BED output (`-f bed`) and summaries (`-s`) only; BED is not converted to GFF3 or GTF
- Plain or gzip-compressed files read directly (`ann.gtf.gz`)
- `-` reads stdin, plain or gzip-compressed

## Parsing warnings

Malformed input does not stop the run. Stderr reports, once per file:

- lines skipped for having too few columns
- lines skipped for invalid coordinates
- CDS lines whose phase is not 0, 1 or 2
- Parent links that form a cycle

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

```bash
# records overlapping chr1:200-600
./gffsub demo.gff3 -r chr1:200-600

# tx01 plus its exons and CDS
./gffsub demo.gff3 -i tx01 -C

# exon rows only
./gffsub demo.gff3 -t exon

# GTF3 output
./gffsub demo.gff3 --format gtf

# one transcript per gene
./gffsub demo.gff3 --longest
```

Filters combine in one invocation. Against a full annotation file:

```bash
# protein-coding genes, longest isoform, full gene model, GTF output
./gffsub ann.gff3 -w biotype=protein_coding -L -m -f gtf -o coding.gtf

# chr1 outside chr1:500000-800000, no dangling Parent references
./gffsub ann.gff3 -S chr1 --exclude-region chr1:500000-800000 --drop-orphans

# IDs from a file, full gene models, column 9 reduced to Name and biotype
./gffsub ann.gff3 --ids genes.txt -m --out-attrs Name,biotype

# gene nearest to a 500 bp region, with its full model
./gffsub ann.gff3 -N chr1:1000000-1000500 -m

# gzipped input, sorted by seqid then start, gene rows as BED
./gffsub ann.gff3.gz -k seqid,natural-seqid,start -t gene -f bed

# second pass reads the first pass from stdin
./gffsub ann.gff3 -w biotype=protein_coding -C | ./gffsub - -L -t mRNA
```

## Selector combination

In the default command, selectors (`-i`, `--ids`, `-n`, `-w`, `-N`) choose candidate records and filters (`-r`, `-b`, `-t`, ...) narrow them. Both sides must match.

```bash
# ex01 spans chr1:100-250, so it does not overlap chr1:600-650
./gffsub demo.gff3 -i ex01 -r chr1:600-650
```

```
##gff-version 3
```

The `query` subcommand adds records from every selector instead. The same flags return ex01 by ID plus the records overlapping the region:

```bash
./gffsub query demo.gff3 -i ex01 -r chr1:600-650
```

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
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
