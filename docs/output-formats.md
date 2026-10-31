# Output Formats

<!-- I18N:START -->

**English** | [中文](./output-formats.zh.md)

<!-- I18N:END -->

## -f / --format

Syntax: `-f FMT` or `--format FMT` (alias `--output-format FMT`).

| Value | Output |
|-------|--------|
| `gff3` | GFF3, original column 9 (default) |
| `gtf`, `gtf3` | GTF3 vocabulary: `gene`, `transcript`, `exon`, `CDS`, `start_codon`, `stop_codon`, `five_prime_utr`, `three_prime_utr`, `Selenocysteine` |
| `gtf25` | GTF2.5 vocabulary: as GTF3, but both UTR types collapse to `UTR` |
| `gtf22` | GTF2.2 vocabulary: `CDS`, `exon`, `start_codon`, `stop_codon`, `5UTR`, `3UTR`, `inter`, `inter_CNS`, `intron_CNS`; no gene or transcript rows |
| `bed` | 6-column BED |

GTF version history and who reads which:

| Version | Types | First shipped by |
|---------|-------|------------------|
| GTF2.2 | CDS, exon, start_codon, stop_codon, 5UTR, 3UTR, inter, inter_CNS, intron_CNS | published spec, 2007 |
| GTF2.5 | gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, UTR | Ensembl release 75 (unofficial) |
| GTF3 | gene, transcript, exon, CDS, Selenocysteine, start_codon, stop_codon, five_prime_utr, three_prime_utr | Ensembl release 82 (unofficial) |

Ensembl and GENCODE distribute GTF3 today, and StringTie, featureCounts and ballgown read it. GTF2.2 is the only version with a published specification.

## GTF output

All three GTF flavors share these rules:

- `gene_id` on every line, `transcript_id` on every non-gene line, both first in column 9.
- `gene_id` and `transcript_id` are taken from the column-9 attribute of that name when the record has one, otherwise from `ID` / `Parent`. Ensembl GFF3 labels a gene `gene:gene-X` and carries `gene_id=gene-X`; the GTF output uses `gene-X`.
- Remaining attributes follow the mandatory pair, URL-decoded from GFF3 and quoted.
- Values are escaped: `"`, `\`, tab, newline, CR.
- The terminal CDS ends before its `stop_codon` when the codon abuts the CDS 3' end of the same transcript.
- Types outside the selected vocabulary are dropped.
- `mRNA` is renamed to `transcript` for GTF3 and GTF2.5.

Header written per flavor: `##gtf-version 3`, `##gtf-version 2.5`, `##gtf-version 2.2`.

## GFF3

- Header: `##gff-version 3`
- Column 9 preserved from input.
- GTF input: column 9 rewritten as `tag=value` with `ID=`/`Parent=` and URL escaping.

## BED

- 6 columns: chrom, start (0-based), end, name, score, strand.
- `start` = record start - 1.
- `score` from column 6, or `0` when missing or `.`.
- `name` from `ID`, or the type when `ID` is missing.

## Example

Sample input (demo.gff3):

```
##gff-version 3
chr1	src	gene	100	1000	.	+	.	ID=gene01;Name=BRCA1
chr1	src	mRNA	100	1000	.	+	.	ID=tx01;Parent=gene01
chr1	src	exon	100	250	.	+	.	ID=ex01;Parent=tx01
chr1	src	CDS	100	250	.	+	0	ID=cds01;Parent=tx01
chr1	src	stop_codon	248	250	.	+	0	ID=stop01;Parent=tx01
```

GTF3 (`-f gtf`):

```
##gtf-version 3
chr1	src	gene	100	1000	.	+	.	gene_id "gene01"; Name "BRCA1";
chr1	src	transcript	100	1000	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	247	.	+	0	gene_id "gene01"; transcript_id "tx01";
chr1	src	stop_codon	248	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

GTF2.2 (`-f gtf22`) drops the gene and transcript rows. The CDS is still trimmed:

```
##gtf-version 2.2
chr1	src	exon	100	250	.	+	.	gene_id "gene01"; transcript_id "tx01";
chr1	src	CDS	100	247	.	+	0	gene_id "gene01"; transcript_id "tx01";
chr1	src	stop_codon	248	250	.	+	0	gene_id "gene01"; transcript_id "tx01";
```

BED:

```
chr1	99	1000	gene01	0	+
chr1	99	1000	tx01	0	+
chr1	99	250	ex01	0	+
chr1	99	250	cds01	0	+
chr1	247	250	stop01	0	+
```

```bash
./gffsub demo.gff3 --format gtf
./gffsub demo.gff3 --format gtf25
./gffsub demo.gff3 --format gtf22
./gffsub demo.gff3 --format bed
./gffsub input.gtf --format gff3
./gffsub demo.gff3 --format gtf -o output.gtf
```

## --out-attrs

Emit only the listed column-9 tags (comma-separated). `ID` and `Parent` are always kept so the feature tree stays consistent.

```bash
./gffsub demo.gff3 --out-attrs Name
./gffsub demo.gff3 --out-attrs Name,gene_id --format gtf
```

Records whose column 9 is `.` are unaffected.
