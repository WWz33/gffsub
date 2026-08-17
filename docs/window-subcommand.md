# Window Subcommand

<!-- I18N:START -->

**English** | [中文](./window-subcommand.zh.md)

<!-- I18N:END -->

## Usage

```
gffsub window <input> --id ID [--up N] [--down N] [--strand-aware]
```

Shortcut form (default mode):

```
gffsub <input> --id ID --up N --down N
```

## --id ID

Target feature ID or gene name. An ID match takes priority; if none matches, gene-name lookup runs next (Name, gene_id, locus_tag, Alias, Dbxref). Required.

## --up N / --upstream N

Bases upstream of the feature start. Default 0. Must be a non-negative integer; trailing garbage (e.g. `50abc`) is rejected.

## --down N / --downstream N

Bases downstream of the feature end. Default 0. Must be a non-negative integer.

## --strand-aware

- Plus strand: upstream extends left of start, downstream extends right of end.
- Minus strand: upstream extends right of end, downstream extends left of start.
- Without the flag: upstream always left, downstream always right.
- Window start below 1 is clamped to 1.

## Discontinuous features

Multi-line CDS sharing one ID: window covers the full span across all segments.

## Example

Sample input (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

```bash
# 200 bp upstream, 500 bp downstream of gene01
# window chr1:1-1500 (start 100-200=-100, clamped to 1; end 1000+500=1500)
./gffsub window demo.gff3 --id gene01 --up 200 --down 500

# strand-aware window for minus-strand gene02
# window chr2:1-800 (downstream 500 left of start: 200-500=-300, clamped to 1; upstream 200 right of end: 600+200=800)
./gffsub window demo.gff3 --id gene02 --up 200 --down 500 --strand-aware

# lookup by gene name
./gffsub window demo.gff3 --id BRCA1 --up 1000 --down 1000

# shortcut form in default mode
./gffsub demo.gff3 --id gene01 --up 200 --down 500
```

Shortcut requires exactly one `--id` plus any of `--up`, `--down`, `--strand-aware`. No other filter flags are accepted; passing one aborts with an error naming the allowed set.
