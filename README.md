# gffsub

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

Subset GFF3 and GTF annotation files by region, feature ID, attribute, or gene model. Pick longest isoforms, switch formats, write summaries.

## Getting Started

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j

# minimal GFF3
printf '%s\n' '##gff-version 3' \
  $'chr1\t.\tgene\t1\t100\t.\t+\t.\tID=GeneA;Name=GeneA' \
  $'chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=GeneA' \
  $'chr1\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1' > data/smoke.gff3

./gffsub data/smoke.gff3 -r chr1:1-100 -f gene
./gffsub data/smoke.gff3 --id GeneA -C
```

## Usage

```text
gffsub <input.gff3> [options]
gffsub query|window <input.gff3> [options]
```

See `./gffsub -h` for the full flag list.

### Common options

| Flag | Default | Description |
|------|---------|-------------|
| `-r, --region` | — | 1-based inclusive `CHR:START-END` |
| `-b, --bed` | — | BED intervals (0-based half-open) |
| `--id` / `--ids` | — | exact feature ID(s) |
| `--name` | — | gene by common naming keys |
| `--where KEY=VALUE` | — | exact column-9 attribute |
| `-C, --children` | off | include descendants of selector hits |
| `--parents` | off | include ancestors |
| `--model` | off | full gene model |
| `--nearest REGION` | — | nearest gene on same seqid |
| `--up` / `--down` | 0 | window (bp) around `--id` |
| `--strand-aware` | off | window follows feature strand |
| `-f, --feature` | — | type column filter |
| `-L, --longest` | off | one transcript/gene: CDS length if any CDS, else exon |
| `-@, --threads` | 1 | threads for `--longest` |
| `-t, --format` | gff3 | `gff3\|gtf\|gtf2\|gtf3\|bed` |
| `-o, --output` | stdout | output file |

## Input / Output

| Kind | Notes |
|------|--------|
| Annotation | GFF3/GTF-style feature records |
| GFF/GTF / `-r` | 1-based inclusive |
| BED | 0-based half-open |
| `--summary` | `tsv` or `json` rows instead of features |

Annotation records only. FASTA and VCF are not supported.

## Documentation

Full guides in [`docs/`](docs/Home.md):

- [Getting Started](docs/getting-started.md)
- [Region and BED Filtering](docs/region-bed-filtering.md)
- [Feature ID and Gene Name Lookup](docs/id-name-lookup.md)
- [Gene Model Expansion](docs/gene-model-expansion.md)
- [Attribute Filtering](docs/attribute-filtering.md)
- [Expression Filtering](docs/expression-filtering.md)
- [Longest Isoform Selection](docs/longest-isoform.md)
- [Window Subcommand](docs/window-subcommand.md)
- [Output Formats](docs/output-formats.md)
- [Summary Output](docs/summary-output.md)
- [GTF Input Handling](docs/gtf-input.md)

## License

MIT License
