# Summary Output

## Overview

Instead of printing annotation records, gffsub can output a tabular summary of selector hits. This is useful for inspecting gene structure details.

## --summary flag

Syntax: `--summary FMT`

| Format | Value |
|--------|-------|
| TSV | `tsv` |
| JSON | `json` |

Default: GFF3 annotation output (no summary).

## --out-attrs: selected attributes

Syntax: `--out-attrs KEYS` (comma-separated)

Aliases: `--output-attrs`, `--attrs`

Adds selected column-9 attributes as extra columns in the summary output. KEYS is a comma-separated list of attribute names.

## TSV columns

| Column | Description |
|--------|-------------|
| query_id | the ID or name queried |
| matched_id | the matched feature ID |
| matched_by | lookup key that matched (ID, Name, gene_id, etc.) |
| seqid | column 1 |
| start | column 4 |
| end | column 5 |
| strand | column 7 |
| type | column 9 |
| parent_id | feature ID of parent |
| child_count | direct children count |
| transcript_count | count of child transcripts |
| exon_count | count of child exons |
| cds_length | total CDS length |
| status | found or not_found |
| (extra) | one per --out-attrs key |

## JSON output

Array of objects. Each object has the same fields as TSV columns, plus an `attrs` object containing the requested attribute keys.

Control characters in field values are escaped per RFC 8259 (short escapes for backspace, formfeed, newline, carriage return, tab; \uXXXX for other chars below 0x20).

## Sample data

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	exon	500	750	.	+	.	ID=ex02;Parent=tx01
chr1	src	CDS	500	750	.	+	2	ID=cds02;Parent=tx01
```

## Commands

```bash
# TSV summary of one gene
./gffsub query demo.gff3 --id gene01 --summary tsv

# JSON summary with selected attributes
./gffsub query demo.gff3 --id gene01 --summary json --out-attrs ID,Name,biotype

# TSV with attributes
./gffsub query demo.gff3 --id gene01 --summary tsv --out-attrs Name,biotype
```

## TSV output example

```tsv
query_id	matched_id	matched_by	seqid	start	end	strand	type	parent_id	child_count	transcript_count	exon_count	cds_length	status
gene01	gene01	ID	chr1	100	1000	+	gene		1	1	2	402	found
```

## JSON output example

```json
[
  {"query_id":"gene01","matched_id":"gene01","matched_by":"ID","seqid":"chr1","start":100,"end":1000,"strand":"+","type":"gene","parent_id":"","child_count":1,"transcript_count":1,"exon_count":2,"cds_length":402,"status":"found","attrs":{"Name":"BRCA1","biotype":"protein_coding"}}
]
```

## not_found status

When a query ID does not match any feature, the summary includes a row with status `not_found` and empty fields. This allows batch queries to report missing IDs.

## Restrictions

--summary and --out-attrs are only valid with query-style selectors (--id, --ids, --name, --where, --grep). They cannot be combined with --seqid, --source, --score, --strand, --phase, --bed, --longest, --threads, --format/--output-format, or --output.
