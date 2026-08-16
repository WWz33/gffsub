# gffsub Documentation

Subset GFF3 and GTF annotation files by region, feature ID, attribute, or gene model.

## Getting Started

- [Installation and Quick Start](getting-started.md)

## Filtering

1. [Region and BED Filtering](region-bed-filtering.md): `-r`, `-b`, `--seqid`
2. [Feature ID and Gene Name Lookup](id-name-lookup.md): `--id`, `--ids`, `--name`
3. [Gene Model Expansion](gene-model-expansion.md): `-C`, `--parents`, `--model`
4. [Attribute Filtering](attribute-filtering.md): `--where`, `--grep`, `--grep-regex`, `--grep-file`
5. [Expression Filtering](expression-filtering.md): `-I`, `-E`

## Selection

- [Longest Isoform Selection](longest-isoform.md): `--longest`, `-@`
- [Window Subcommand](window-subcommand.md): `window`, `--up`, `--down`, `--strand-aware`

## Output

- [Output Formats](output-formats.md): GFF3, GTF, GTF2, GTF3, BED (`-t`)
- [Summary Output](summary-output.md): TSV and JSON tabular summaries (`--summary`, `--out-attrs`)

## Input

- [GTF Input Handling](gtf-input.md): GTF attribute parsing, Parent/ID synthesis, isoform auto-detection

## Library

- [Library API](library-api.md): link `libgffsub_core.a` and call `query()`, `window()`, `subset()` from C++17
