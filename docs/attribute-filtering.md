# Attribute Filtering

<!-- I18N:START -->

**English** | [中文](./attribute-filtering.zh.md)

<!-- I18N:END -->

## --where KEY=VALUE

Exact match on a column-9 attribute. Repeatable; multiple flags are OR'd.

```bash
./gffsub demo.gff3 --where biotype=protein_coding
```

## --grep FIELD:PATTERN

Substring match on a field or attribute. Repeatable; multiple flags are OR'd.

Column fields: `seqid`, `source`, `type`, `start`, `end`, `length`, `score`, `strand`, `phase`, `attrs`.

Attribute fields: `ID`, `Name`, `Parent`, `Alias`, `Dbxref`, `Note`, `biotype`, `gene_id`, `transcript_id`, `locus_tag`. Any other attribute uses the `attr.KEY` form.

```bash
./gffsub demo.gff3 --grep Name:ABC
```

## --grep-regex FIELD:REGEX

Regex match (ECMAScript). Repeatable; multiple flags are OR'd.

```bash
./gffsub demo.gff3 --grep-regex ID:exon.*
```

## --grep-file FILE

Read patterns from FILE, one per line. Requires `--grep-field FIELD` to set the target field.

- Without `--grep-file-regex`: substring matching.
- With `--grep-file-regex`: regex matching.

Patterns file (patterns.txt):

```
BRCA
TP53
```

```bash
./gffsub demo.gff3 --grep-field Name --grep-file patterns.txt
```

## -v / --invert-match

Invert grep matches: keep records that do NOT match.

```bash
./gffsub demo.gff3 --grep type:mRNA -v
```

## --ignore-case

Case-insensitive matching for `--grep`, `--grep-regex`, `--grep-file`, `-I`, and `-E`.

```bash
./gffsub demo.gff3 --grep type:mrna --ignore-case
```

See also: [Expression Filtering](expression-filtering.md) for numeric comparison and logical operators.
