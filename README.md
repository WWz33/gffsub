# gffsub

<!-- README-I18N:START -->

**English** | [中文](./README.zh.md)

<!-- README-I18N:END -->

Subset GFF3 and GTF annotation files by region, feature ID, attribute, or gene model. Pick longest isoforms, switch formats, write summaries.

## Getting Started

```bash
git clone https://github.com/WWz33/gffsub.git
cd gffsub && make -j
```

Minimal GFF3:

```bash
printf '%s\n' '##gff-version 3' \
  $'chr1\t.\tgene\t1\t100\t.\t+\t.\tID=GeneA;Name=GeneA' \
  $'chr1\t.\tmRNA\t1\t100\t.\t+\t.\tID=mRNA1;Parent=GeneA' \
  $'chr1\t.\texon\t1\t100\t.\t+\t.\tID=exon1;Parent=mRNA1' > data/smoke.gff3

./gffsub data/smoke.gff3 -r chr1:1-100 -t gene
./gffsub data/smoke.gff3 -i GeneA -C
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
| `-i, --id` / `--ids` | — | exact feature ID(s) |
| `-n, --name` | — | gene by common naming keys |
| `-w, --where KEY=VALUE` | — | exact column-9 attribute |
| `-C, --children` | off | include descendants of selector hits |
| `-p, --parents` | off | include ancestors |
| `-m, --model` | off | full gene model |
| `-N, --nearest REGION` | — | nearest gene on same seqid |
| `-u, --up` / `-D, --down` | 0 | window (bp) around `-i` |
| `-a, --strand-aware` | off | window follows feature strand |
| `-t, --type` | — | type column filter; comma list, `^LIST` excludes; repeatable |
| `-L, --longest` | off | one transcript/gene: CDS length if any CDS, else exon |
| `--longest-type` | auto | isoform type for `-L`; auto-detects transcript class |
| `-@, --threads` | 1 | threads for `-L` |
| `--format` | gff3 | `gff3\|gtf\|gtf2\|gtf3\|bed` |
| `-o, --output` | stdout | output file |
| `-s, --summary` | off | TSV stats by seqid × type (like `seqkit stats`) |
| `-S, --seqid` | — | keep seqids; `^LIST` excludes |

## Input / Output

| Kind | Notes |
|------|--------|
| Annotation | GFF3/GTF-style feature records |
| GFF/GTF / `-r` | 1-based inclusive |
| BED | 0-based half-open |

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
- [Library API](docs/library-api.md)

## License

MIT License
