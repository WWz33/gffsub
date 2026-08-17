# Feature ID and Gene Name Lookup

<!-- I18N:START -->

**English** | [中文](./id-name-lookup.zh.md)

<!-- I18N:END -->

Sample data (`demo.gff3`):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;Alias=BRCA-1;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	mRNA	100	800	.	+	.	ID=tx02;Parent=gene01
chr1	src	exon	100	800	.	+	.	ID=ex02;Parent=tx02
```

## --id

Syntax: `--id ID` (repeatable)

- Exact feature ID match
- Multiple `--id` flags OR together

```bash
./gffsub demo.gff3 --id gene01
./gffsub demo.gff3 --id tx01 --id tx02
./gffsub demo.gff3 --id gene01 -C
```

## --ids

Syntax: `--ids FILE`

- One ID per line
- Blank lines and lines starting with `#` skipped

Example file (`ids.txt`):

```
gene01
tx01
```

```bash
./gffsub demo.gff3 --ids ids.txt -C
```

## --name

Syntax: `--name NAME`

- Gene lookup by multiple keys: `ID`, `gene_id`, `Name`, `locus_tag`, `Alias`, `Dbxref`
- Only matches gene-type features

```bash
./gffsub demo.gff3 --name BRCA1
./gffsub demo.gff3 --name BRCA-1
./gffsub demo.gff3 --name GeneID:672
./gffsub demo.gff3 --name gene01 -C
```

## --nearest REGION

Syntax: `--nearest CHR:START-END`

- Finds nearest gene on the same seqid to a 1-based region

```bash
./gffsub demo.gff3 --nearest chr1:1500-2000
```

See also: [Gene Model Expansion](gene-model-expansion.md), [Summary Output](summary-output.md).
