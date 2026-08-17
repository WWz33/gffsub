# Summary Output

<!-- I18N:START -->

**English** | [中文](./summary-output.zh.md)

<!-- I18N:END -->

## --summary FMT

Values: `tsv`, `json`. Outputs tabular rows instead of feature records. Only valid with query-style selectors (`--id`, `--ids`, `--name`, `--where`, `--nearest`, query-subcommand `--region`). Cannot be combined with `--grep`, `-I`, `-E`, `--seqid`, `--source`, `--score`, `--strand`, `--phase`, `--bed`, `--longest`, `--threads`, `--format`/`--output-format`, or `--output`.

## --out-attrs KEYS

Comma-separated attribute keys appended as extra columns (aliases `--output-attrs`, `--attrs`). Resolves `gene_id`/`transcript_id` from record fields first, then GFF3 attributes, then GTF quoted values.

## TSV columns

`query_id`, `matched_id`, `matched_by`, `seqid`, `start`, `end`, `strand`, `type`, `parent_id`, `child_count`, `transcript_count`, `exon_count`, `cds_length`, `status`, plus one column per `--out-attrs` key. Tab, newline, and CR in values escaped as `\t`, `\n`, `\r`.

## JSON format

Array of objects with the same fields plus an `attrs` object holding the `--out-attrs` keys.

## status field

`found` or `not_found`. On a miss, `matched_by` keeps the lookup key, numeric counts are 0, other fields empty.

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
```

```bash
./gffsub demo.gff3 --id gene01 --summary tsv
./gffsub demo.gff3 --id gene01 --summary json --out-attrs Name,biotype
```

TSV output:
```tsv
query_id	matched_id	matched_by	seqid	start	end	strand	type	parent_id	child_count	transcript_count	exon_count	cds_length	status
gene01	gene01	ID	chr1	100	1000	+	gene		1	1	2	402	found
```

JSON output:
```json
[
  {"query_id":"gene01","matched_id":"gene01","matched_by":"ID","seqid":"chr1","start":100,"end":1000,"strand":"+","type":"gene","parent_id":"","child_count":1,"transcript_count":1,"exon_count":2,"cds_length":402,"status":"found","attrs":{"Name":"BRCA1","biotype":"protein_coding"}}
]
```
