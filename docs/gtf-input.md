# GTF Input Handling

<!-- I18N:START -->

**English** | [中文](./gtf-input.zh.md)

<!-- I18N:END -->

## Format detection

gffsub sniffs the first data line: the first column-9 fragment decides. `tag=value` is GFF3, `tag "value"` is GTF, 3-12 columns with integer start/end is BED, default GFF3.

The first fragment decides, so punctuation inside a value does not flip the verdict: `ID=g1;Name="A"` stays GFF3, `note "a=b"` stays GTF.

## GTF vs GFF3 attributes

- GFF3: `key=value;` pairs.
- GTF: `key "value";` pairs.
- `gene_id` and `transcript_id` extracted to record fields.

GTF column 9 details:

- `#` outside a quoted value starts a comment that runs to the end of the line; a `#` inside quotes is data.
- `\"`, `\\`, `\t`, `\n`, `\r` are decoded.
- Unquoted values are accepted.
- A bare `ID=x;Parent=y` fragment (GFF3 style pasted into GTF) is accepted; `key="value"` keeps `key` as the tag name.
- A quoted multi-value `Parent "t1,t2"` is split on commas. The record links to both transcripts, `--where Parent=t1` matches, and GTF output writes one row per parent.

## Parent/child synthesis

GTF has no `ID=`/`Parent=`. gffsub synthesizes:

| Feature type        | Synthesized ID | Synthesized Parent |
|---------------------|----------------|--------------------|
| gene                | gene_id        | none               |
| transcript / mRNA   | transcript_id  | gene_id            |
| exon, CDS, etc.     | none           | transcript_id      |

Enables `-C`, `-p`, `-m` on GTF.

## Attribute access

For `--where`, `--grep`, and `-I`: `gene_id` and `transcript_id` resolve from record fields. Other GTF attributes use the `attr.` prefix (e.g. `attr.gene_name`).

```bash
./gffsub input.gtf -w gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## -L / --longest on GTF

Auto-detects the isoform type: the most frequent transcript-class type (`transcript` for GTF). No need for `-t transcript`.

## GTF to GFF3 conversion

Column 9 rewritten as `tag=value`. `ID=`/`Parent=` synthesized from `gene_id`/`transcript_id`. Other attributes URL-escaped.

```bash
./gffsub input.gtf --format gff3
```

## GTF to GTF

`gene_id` on every line. `transcript_id` on non-gene lines only. Values escaped.

Sample input (demo.gtf):
```
chr1	src	gene	100	1000	.	+	.	gene_id "g1"; gene_name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	100	250	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	exon	500	750	.	+	.	gene_id "g1"; transcript_id "t1";
chr1	src	CDS	100	250	.	+	0	gene_id "g1"; transcript_id "t1";
```

```bash
./gffsub demo.gtf -w gene_id=g1
./gffsub demo.gtf -i t1 -C
./gffsub demo.gtf --longest
./gffsub demo.gtf -i g1 -C --format gff3
```
