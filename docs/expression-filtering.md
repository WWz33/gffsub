# Expression Filtering

<!-- I18N:START -->

**English** | [中文](./expression-filtering.zh.md)

<!-- I18N:END -->

## Syntax

```
-I 'FIELD OP VALUE'   # keep matching records
-E 'FIELD OP VALUE'   # drop matching records
```

Multiple `-I` and `-E` flags combine with AND logic. Within one expression: `&&` (and), `||` (or), `!` (not), `()`.

Quote values with spaces or special characters using double quotes.

## Operators

| Operator | Meaning |
|----------|---------|
| `==` | string equal |
| `!=` | string not equal |
| `~` | regex match |
| `!~` | regex not match |
| `<` | numeric less than |
| `<=` | numeric less or equal |
| `>` | numeric greater than |
| `>=` | numeric greater or equal |

## Logic examples

```bash
-I 'length >= 100 && strand == "+"'
-I 'type == "gene" || type == "mRNA"'
-I '!(type == "gene")'
```

## Built-in fields

| Field | Source |
|-------|--------|
| `seqid` | column 1 |
| `source` | column 2 |
| `type` | column 3 |
| `start` | column 4 |
| `end` | column 5 |
| `length` | `end - start + 1` |
| `score` | column 6 |
| `strand` | column 7 |
| `phase` | column 8 |
| `attrs` | column 9 (raw) |

## Attribute fields

| Field | GFF3 key |
|-------|----------|
| `ID` | `ID` |
| `Name` | `Name` |
| `Parent` | `Parent` |
| `Alias` | `Alias` |
| `Dbxref` | `Dbxref` |
| `Note` | `Note` |
| `biotype` | `biotype` |
| `gene_id` | `gene_id` |
| `transcript_id` | `transcript_id` |
| `locus_tag` | `locus_tag` |
| `attr.KEY` | any column-9 attribute |

GTF input: `gene_id` / `transcript_id` resolve from record fields. Other attributes use the `attr.` prefix.

## Numeric comparison

`start`, `end`, `length`, and `score` are numeric for `<`, `<=`, `>`, `>=`, `==`, `!=`.

## Missing attributes

- `==`, `~`, `<`, `<=`, `>`, `>=` do not match (record excluded by `-I`).
- `!=` and `!~` match (record kept by `-I`).

## Examples

```bash
gffsub ann.gff3 -f exon -I 'length >= 100'
gffsub ann.gff3 -E 'length < 10'
gffsub ann.gff3 -f exon -I 'length >= 200 && strand == "+"'
```
