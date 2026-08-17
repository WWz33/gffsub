# GTF Input Handling

<!-- I18N:START -->

**English** | [中文](./gtf-input.zh.md)

<!-- I18N:END -->

## Format detection

gffsub sniffs file content: column 9 with `=` and `;` is GFF3; column 9 with `"` is GTF; 3-12 columns with integer start/end is BED; default GFF3.

## GTF vs GFF3 attributes

- GFF3: `key=value;` pairs.
- GTF: `key "value";` pairs.
- `gene_id` and `transcript_id` extracted to record fields.

## Parent/child synthesis

GTF has no `ID=`/`Parent=`. gffsub synthesizes:

| Feature type        | Synthesized ID | Synthesized Parent |
|---------------------|----------------|--------------------|
| gene                | gene_id        | none               |
| transcript / mRNA   | transcript_id  | gene_id            |
| exon, CDS, etc.     | none           | transcript_id      |

Enables `--children`, `--parents`, `--model` on GTF.

## Attribute access

For `--where`, `--grep`, and `-I`: `gene_id` and `transcript_id` resolve from record fields. Other GTF attributes use the `attr.` prefix (e.g. `attr.gene_name`).

```bash
./gffsub input.gtf --where gene_id=ENSG000001
./gffsub input.gtf --grep gene_id:ENSG
./gffsub input.gtf -I 'gene_id == "ENSG000001"'
```

## --longest on GTF

Auto-detects the isoform type: uses `transcript` if no `mRNA` records exist. No need for `-f transcript`.

## GTF to GFF3 conversion

Column 9 rewritten as `tag=value`. `ID=`/`Parent=` synthesized from `gene_id`/`transcript_id`. Other attributes URL-escaped.

```bash
./gffsub input.gtf -t gff3
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
./gffsub demo.gtf --where gene_id=g1
./gffsub demo.gtf --id t1 -C
./gffsub demo.gtf --longest
./gffsub demo.gtf --id g1 -C -t gff3
```
