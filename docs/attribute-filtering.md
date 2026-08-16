# Attribute Filtering

## --where: exact attribute match

Syntax: `--where KEY=VALUE` (repeatable)

- Matches GFF3 column-9 attributes exactly
- For GTF input, gene_id and transcript_id are accessible (synthesized by the parser)
- Multiple --where flags are OR'd (a record matching any condition is kept)
- Attribute values are URL-decoded before comparison (GFF3 spec): a file value `Note=a%2Cb` is matched by the decoded query form `--where Note=a,b`

Sample data (demo.gff3):
```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1;biotype=protein_coding;Dbxref=GeneID:672
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01;biotype=protein_coding
chr2	src	gene	200	600	.	-	.	ID=gene02;Name=TP53;biotype=tumor_suppressor
```

Commands:
1. `./gffsub demo.gff3 --where biotype=protein_coding -f gene` — gene with protein_coding biotype
2. `./gffsub demo.gff3 --where biotype=tumor_suppressor` — all records with that biotype
3. `./gffsub demo.gff3 --where biotype=protein_coding --where Name=BRCA1` — either condition may match

## --grep: substring search

Syntax: `--grep FIELD:PATTERN` (repeatable)

- Searches for substring PATTERN within the named FIELD
- Fields: column fields (seqid, source, type, start, end, score, strand, phase) and column-9 attribute keys (ID, Name, Parent, gene_id, transcript_id, biotype, etc.)
- Keys not in the predefined list need the `attr.` prefix
- Multiple --grep flags are OR'd (any match keeps the record)
- `-v` / `--invert-match` inverts: keep records that do NOT match

Commands:
1. `./gffsub demo.gff3 --grep type:gene` — genes (substring "gene" in type column)
2. `./gffsub demo.gff3 --grep Name:BRCA` — records with "BRCA" in Name attribute
3. `./gffsub demo.gff3 --grep type:mRNA -v` — everything except mRNA records

## --grep-regex: regex search

Syntax: `--grep-regex FIELD:REGEX` (repeatable)

- ECMAScript regex syntax (std::regex)
- Same field list as --grep
- Multiple flags are OR'd

Commands:
1. `./gffsub demo.gff3 --grep-regex ID:ex0.` — records whose ID matches the regex
2. `./gffsub demo.gff3 --grep-regex Name:^BRCA` — names starting with BRCA

## --grep-file: batch patterns

Syntax: `--grep-field FIELD --grep-file FILE [--grep-file-regex]`

- Reads patterns from FILE (one per line)
- Matches against the single FIELD specified by --grep-field
- Without --grep-file-regex: substring matching
- With --grep-file-regex: regex matching

Example file (patterns.txt):
```
BRCA
TP53
```
Command: `./gffsub demo.gff3 --grep-field Name --grep-file patterns.txt`

## --ignore-case

Syntax: `--ignore-case`

- Makes --grep, --grep-regex, --grep-file, -I, and -E case-insensitive

Command: `./gffsub demo.gff3 --grep type:mrna --ignore-case`

## Field reference

| Field | Source |
|-------|--------|
| seqid | column 1 |
| source | column 2 |
| type | column 3 |
| start | column 4 |
| end | column 5 |
| score | column 6 |
| strand | column 7 |
| phase | column 8 |
| ID | column 9 attribute |
| Name | column 9 attribute |
| Parent | column 9 attribute |
| Alias | column 9 attribute |
| Dbxref | column 9 attribute |
| Note | column 9 attribute |
| gene_id | column 9 attribute |
| transcript_id | column 9 attribute |
| biotype | column 9 attribute |
| locus_tag | column 9 attribute |
| attr.KEY | column 9 attribute (any key) |

See also: [Expression Filtering](expression-filtering.md) for numeric comparison, logical operators, and complex predicates.
