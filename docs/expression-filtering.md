# Expression Filtering

<!-- I18N:START -->

**English** | [中文](./expression-filtering.zh.md)

<!-- I18N:END -->

gffsub filters records by field values, numeric comparisons, and attribute matching using `-I` / `--include-expr` (keep matching) and `-E` / `--exclude-expr` (drop matching).

## Syntax

```
-I 'FIELD OP VALUE'
```

Multiple `-I` and `-E` flags combine with AND logic. Each expression supports `&&` (and), `||` (or), `!` (not), and parentheses `()`.

### Operators

| Operator | Meaning | Example |
|----------|---------|---------|
| `==` | string equal | `-I 'type == "exon"'` |
| `!=` | string not equal | `-I 'type != "gene"'` |
| `~` | regex match | `-I 'ID ~ "exon.*"'` |
| `!~` | regex not match | `-I 'Name !~ "test.*"'` |
| `<` | numeric less than | `-I 'length < 100'` |
| `<=` | numeric less or equal | `-I 'length <= 200'` |
| `>` | numeric greater than | `-I 'length > 500'` |
| `>=` | numeric greater or equal | `-I 'length >= 1000'` |

Quote values containing spaces or special characters with double quotes.

### Logical operators

Combine conditions within a single expression:

```bash
# exon length >= 100 on plus strand
-I 'length >= 100 && strand == "+"'

# gene or transcript
-I 'type == "gene" || type == "mRNA"'

# not a gene
-I '!(type == "gene")'
```

## Fields

### Built-in fields

| Field | Type | Description |
|-------|------|-------------|
| `seqid` | string | column 1 (sequence name) |
| `source` | string | column 2 (source) |
| `type` | string | column 3 (feature type) |
| `start` | numeric | column 4 (1-based start) |
| `end` | numeric | column 5 (1-based end) |
| `length` | numeric | `end - start + 1` |
| `score` | numeric | column 6 |
| `strand` | string | column 7 (`+ - . ?`) |
| `phase` | string | column 8 (`0 1 2 .`) |
| `attrs` | string | raw column 9 |

### Attribute fields

GFF3 column-9 attributes are accessible by key name. Keys not in the built-in list need the `attr.` prefix.

| Field | GFF3 key | Notes |
|-------|----------|-------|
| `ID` | `ID` | feature identifier |
| `Name` | `Name` | display name |
| `Parent` | `Parent` | parent feature ID |
| `Alias` | `Alias` | alias |
| `Dbxref` | `Dbxref` | database cross-reference |
| `Note` | `Note` | free-text note |
| `gene_id` | `gene_id` | gene ID (GTF / GFF3) |
| `transcript_id` | `transcript_id` | transcript ID (GTF / GFF3) |
| `biotype` | `biotype` | biotype |
| `locus_tag` | `locus_tag` | locus tag |
| `attr.KEY` | `KEY` | any column-9 attribute |

GTF input: `gene_id`, `transcript_id`, `ID`, and `Parent` resolve from the parser-synthesized fields. Other GTF attributes use the `attr.` prefix with the GTF key name.

## Numeric comparisons

Fields `start`, `end`, `length`, and `score` are numeric for `<`, `<=`, `>`, `>=`, `==`, and `!=`. Non-numeric values fail the comparison (no crash, record excluded).

## Missing attributes

When a record lacks the referenced attribute:

- `==`, `~`, `<`, `<=`, `>`, `>=` do not match (record excluded by `-I`).
- `!=` and `!~` do match (record kept by `-I`).

So `-I 'attr.biotype != "lncRNA"'` keeps records with no `biotype` attribute, treating missing as "not equal".

## Length filtering examples

Filter features by genomic span:

```bash
# exons at least 100 bp
gffsub ann.gff3 -f exon -I 'length >= 100'

# CDS shorter than 200 bp
gffsub ann.gff3 -f CDS -I 'length < 200'

# mRNAs longer than 1 kb
gffsub ann.gff3 -f mRNA -I 'length > 1000'

# exons between 50 and 200 bp
gffsub ann.gff3 -f exon -I 'length >= 50 && length <= 200'

# drop features shorter than 10 bp
gffsub ann.gff3 -E 'length < 10'
```

## Combining with feature type

`-f` / `--feature` restricts by column 3 type before expression evaluation:

```bash
# long exons on plus strand only
gffsub ann.gff3 -f exon -I 'length >= 200 && strand == "+"'
```
