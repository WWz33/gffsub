# Window Subcommand

## Syntax

```
gffsub window <input.gff3> --id ID [--up N] [--down N] [--strand-aware]
```

Returns all records overlapping the genomic window around one feature.

## Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--id ID` | required | target feature by ID or gene name |
| `--up N` | 0 | bases upstream of the feature start |
| `--down N` | 0 | bases downstream of the feature end |
| `--strand-aware` | off | window follows feature strand |

The `--up` flag has the long alias `--upstream`; `--down` has `--downstream`.

## Window construction

Without --strand-aware (default):
- Window start = feature start - upstream
- Window end = feature end + downstream
- Direction-independent: window extends equally left and right

With --strand-aware:
- Plus strand: upstream extends left of the start, downstream extends right of the end
- Minus strand: the two extensions swap. Upstream extends right of the end, downstream extends left of the start

Any computed start below 1 is clamped to 1 before overlap lookup.

## --id lookup

--id accepts both feature IDs and gene-name keys, using the same lookup as --name. An ID match takes priority. If no ID matches, the gene-name lookup runs next. The indexed gene-name attributes are Name, gene_id, locus_tag, Alias, and Dbxref.

## Sample data

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr2	src	gene	200	600	.	-	.	ID=gene02
```

## Commands

```bash
# 200 bp upstream, 500 bp downstream of gene01
./gffsub window demo.gff3 --id gene01 --up 200 --down 500
# window: chr1:1 to 1500 (start 100 - 200 = -100, clamped to 1; end 1000 + 500 = 1500)

# strand-aware window for a minus-strand gene
./gffsub window demo.gff3 --id gene02 --up 200 --down 500 --strand-aware
# window: chr2:1 to 800 (downstream 500 extends left of start: 200 - 500 = -300, clamped to 1; upstream 200 extends right of end: 600 + 200 = 800)

# lookup by gene name
./gffsub window demo.gff3 --id BRCA1 --up 1000 --down 1000
```

## Window shortcut in default mode

In default (non-subcommand) mode, --up and --down with exactly one --id trigger a window query internally:

```bash
./gffsub demo.gff3 --id gene01 --up 200 --down 500
```

This path requires exactly one --id, plus any of --up, --down, or --strand-aware. No other filter flags are accepted. Passing another filter flag aborts with an error naming the allowed set: --id, --up/--upstream, --down/--downstream, and --strand-aware.

## Validation

- --up and --down must be non-negative integers
- Trailing garbage after the integer (e.g. "50abc") is rejected
- --id is required
