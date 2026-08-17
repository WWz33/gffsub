# GTF Input Handling

<!-- I18N:START -->

**English** | [中文](./gtf-input.zh.md)

<!-- I18N:END -->

## Format detection

gffsub infers input format from the file extension:
- `.gtf` or `.GTF` -> GTF parsing
- `.bed` or `.BED` -> BED parsing
- everything else -> GFF3 parsing

## GTF attribute model

GTF uses a different attribute format from GFF3:
- GFF3: `key=value;` (key and value separated by `=`, pairs separated by `;`)
- GTF: `key "value";` (key and quoted value separated by space, pairs separated by `;`)

gffsub's `parse_attributes` function parses GFF3-style attributes. For GTF, the parser extracts `gene_id` and `transcript_id` into dedicated record fields using a GTF-aware extractor that handles quoted values and escaped quotes.

## Parent/child synthesis

GTF has no explicit `ID=` and `Parent=` attributes. gffsub synthesizes these from `gene_id` and `transcript_id`:

| Feature type | Synthesized ID | Synthesized Parent |
|--------------|---------------|-------------------|
| gene | gene_id | none (no parent) |
| transcript / mRNA | transcript_id | gene_id |
| exon, CDS, etc. | none | transcript_id |

The synthesis builds the parent/child index so `--children`, `--parents`, and `--model` work on GTF input. Exon and CDS records get no synthesized ID; they link to their transcript via the parent field only. This avoids ID collisions when multiple exons share the same transcript_id.

## GTF-specific parsing

- Quoted values: `gene_id "ENSG000001";` extracts the string between quotes.
- Escaped quotes: `gene_id "test\"name";` handles backslash-escaped quotes within the value.
- Empty values: `gene_id "";` is treated as missing (the record gets no gene_id).

## Attribute access for GTF

For `--where`, `--grep`, and `-I` expressions:
- `gene_id` and `transcript_id` resolve from the synthesized record fields, not from parsing column 9.
- Other GTF attributes (like `gene_name`, `gene_biotype`) are accessible via the `attr.` prefix, but they use GTF parsing which extracts the quoted value directly.

Filter by gene_id in GTF:
```bash
./gffsub input.gtf --where gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## --longest on GTF

The isoform type is auto-detected. GTF typically uses `transcript` rather than `mRNA`, so gffsub checks which type exists in the data:
- If `mRNA` records exist, use `mRNA`.
- Otherwise if `transcript` records exist, use `transcript`.

`--longest` works on GTF without needing `-f transcript`.

## GTF output

When converting GTF to GTF, gffsub reads the gene_id and transcript_id, rebuilds the hierarchy, and outputs GTF2.2-compliant attributes:
- gene_id on every line (required by GTF2.2).
- transcript_id on non-gene lines only.
- Gene lines do not get transcript_id (per GTF2.2 and GENCODE convention).
- Attribute values are escaped (backslash and quote).

See [Output Formats](output-formats.md) for details.

## Example

Sample GTF (demo.gtf):
```
chr1	src	gene	100	1000	.	+	.	gene_id "g1"; gene_name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	250	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	500	750	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	250	.	+	0	gene_id "g1"; transcript_id "t1";
```

```bash
# select by gene_id
./gffsub demo.gtf --where gene_id=g1

# children of a transcript
./gffsub demo.gtf --id t1 -C

# longest isoform (auto-detects transcript type)
./gffsub demo.gtf --longest

# convert GTF to GFF3: column 9 is rewritten as tag=value with synthesized
# ID=/Parent= (gene lines get ID=<gene_id>, transcripts ID=<transcript_id>
# +Parent=<gene_id>, children Parent=<transcript_id>); remaining attributes
# are converted to key=value form with URL escaping
./gffsub demo.gtf --id g1 -C -t gff3
```
