# Summary Output

<!-- I18N:START -->

**English** | [中文](./summary-output.zh.md)

<!-- I18N:END -->

## -s, --summary

Output a TSV summary of the filtered records instead of feature records. Works on the main subset path and the `query` subcommand. Counts (`child_count`, `transcript_count`, `exon_count`, `cds_length`) are computed from the full annotation index, not the filtered subset, so parent-child relationships remain intact.

## TSV columns

`seqid`, `start`, `end`, `strand`, `type`, `length`, `child_count`, `transcript_count`, `exon_count`, `cds_length`.

- `length` = `end - start + 1`
- `child_count` = number of direct children (records with `Parent` pointing to this record's ID)
- `transcript_count` = number of mRNA/transcript descendants
- `exon_count` = number of exon descendants
- `cds_length` = sum of CDS descendant lengths

Tab, newline, and CR in values escaped as `\t`, `\n`, `\r`.

## all row

When the filtered records span more than one distinct `seqid`, a final `all` row is emitted with summed numeric columns. Its `start`, `end`, `strand`, and `type` are `NA` to distinguish it from per-record GFF rows. If all records share one `seqid`, the `all` row is omitted.

## Example

Sample input (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
chr2	src	gene	500	600	.	-	.	ID=gene02;Name=XYZ1
```

```bash
./gffsub demo.gff3 -f gene -s
```

Output:
```tsv
seqid	start	end	strand	type	length	child_count	transcript_count	exon_count	cds_length
chr1	100	1000	+	gene	901	1	1	2	402
chr2	500	600	-	gene	101	0	0	0	0
all	NA	NA	NA	NA	1002	1	1	2	402
```
