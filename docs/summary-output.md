# Summary Output

<!-- I18N:START -->

**English** | [中文](./summary-output.zh.md)

<!-- I18N:END -->

## -s, --summary

Output an aggregated TSV summary of the filtered records instead of feature records. Works on the main subset path and the `query` subcommand.

## TSV columns

`seqid`, `type`, `count`, `sum_len`, `min_len`, `avg_len`, `max_len`, `Q1`, `Q2`, `Q3`, `coverage`.

One row per distinct (seqid, type). `length` = `end - start + 1`. `avg_len` is truncated to one decimal. `Q1`/`Q2`/`Q3` are Tukey hinges (lower/upper half medians). `coverage` is the union bp: intervals of the group are merged when overlapping or book-ended (bedtools merge semantics), so bases shared by overlapping records (e.g. one exon in two isoforms) count once. Rows are sorted by seqid, then type.

## all rows

When the filtered records span more than one distinct `seqid`, one `all` row per type is appended: `count`/`sum_len` are totals, `min_len`/`max_len` are global extremes, `avg_len` is the overall mean, `Q1`/`Q2`/`Q3` are computed over all lengths, `coverage` is the sum of per-seqid unions (intervals on different seqids never merge). If all records share one `seqid`, `all` rows are omitted.

## Example

Sample input (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr2	src	gene	500	600	.	-	.	ID=gene02;Name=XYZ1
```

```bash
./gffsub demo.gff3 -s
```

Output:
```tsv
seqid	type	count	sum_len	min_len	avg_len	max_len	Q1	Q2	Q3	coverage
chr1	CDS	2	402	151	201	251	151	201	251	402
chr1	exon	2	402	151	201	251	151	201	251	402
chr1	gene	1	901	901	901	901	901	901	901	901
chr1	mRNA	1	901	901	901	901	901	901	901	901
chr2	gene	1	101	101	101	101	101	101	101	101
all	CDS	2	402	151	201	251	151	201	251	402
all	exon	2	402	151	201	251	151	201	251	402
all	gene	2	1002	101	501	901	101	501	901	1002
all	mRNA	1	901	901	901	901	901	901	901	901
```
