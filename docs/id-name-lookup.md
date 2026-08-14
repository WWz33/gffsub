# Feature ID and Gene Name Lookup

## By ID: --id

Syntax: `--id ID` (repeatable)

- Exact match on the feature's ID attribute (column 9 `ID=` in GFF3, synthesized from gene_id/transcript_id in GTF)
- Multiple --id flags are OR'd: records matching any ID are kept
- Without -C/--children/--parents/--model, only the exact ID match is returned

Sample data (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;Alias=BRCA-1;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01;Name=BRCA1.1
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex02;Parent=tx02
```

Commands:
1. `./gffsub demo.gff3 --id gene01` — the gene record only
2. `./gffsub demo.gff3 --id tx01 --id tx02` — both mRNA records
3. `./gffsub demo.gff3 --id gene01 -C` — gene plus all descendants (see [Gene Model Expansion](gene-model-expansion.md))

## By ID file: --ids

Syntax: `--ids FILE`

- One ID per non-empty line
- Equivalent to multiple --id flags
- Blank lines and lines starting with # are skipped

Example file (ids.txt):
```
gene01
tx01
```
Command: `./gffsub demo.gff3 --ids ids.txt -C`

## By gene name: --name

Syntax: `--name NAME`

- Searches gene records by multiple naming keys (not just the ID attribute)
- Lookup keys (checked in order):
  1. ID
  2. Name
  3. gene_id
  4. locus_tag
  5. Alias
  6. Dbxref
- Only matches gene-type features (type == "gene")
- Returns the gene record; use -C or --model to expand

Commands:
1. `./gffsub demo.gff3 --name BRCA1` — gene with Name=BRCA1
2. `./gffsub demo.gff3 --name BRCA-1` — gene with Alias=BRCA-1
3. `./gffsub demo.gff3 --name GeneID:672` — gene with Dbxref=GeneID:672
4. `./gffsub demo.gff3 --name gene01 -C` — gene by ID, with children

## query subcommand

The `query` subcommand is the preferred interface for ID/name lookup. It adds summary output support.

```
./gffsub query demo.gff3 --id gene01 -C
./gffsub query demo.gff3 --name BRCA1 --summary tsv
./gffsub query demo.gff3 --ids ids.txt --summary json
```

See also: [Gene Model Expansion](gene-model-expansion.md), [Summary Output](summary-output.md).
