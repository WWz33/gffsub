# gffsub Documentation

Subset GFF3 and GTF annotation files by region, feature ID, attribute, or gene model.

## Getting Started

- [Installation and Quick Start](getting-started.md)

## Filtering

Progressive from simple to advanced:

1. [Region and BED Filtering](region-bed-filtering.md): subset by genomic coordinates (`-r`, `-b`, `--seqid`)
2. [Feature ID and Gene Name Lookup](id-name-lookup.md): select by ID or name (`--id`, `--ids`, `--name`)
3. [Gene Model Expansion](gene-model-expansion.md): include children, parents, or full gene model (`-C`, `--parents`, `--model`)
4. [Attribute Filtering](attribute-filtering.md): match column-9 attributes (`--where`, `--grep`, `--grep-regex`, `--grep-file`)
5. [Expression Filtering](expression-filtering.md): numeric comparison and logical operators (`-I`, `-E`)

## Selection

- [Longest Isoform Selection](longest-isoform.md): keep one isoform per gene (`--longest`, `-@`)
- [Window Subcommand](window-subcommand.md): records around a feature (`window`, `--up`, `--down`, `--strand-aware`)

## Output

- [Output Formats](output-formats.md): GFF3, GTF, GTF2, GTF3, BED (`-t`)
- [Summary Output](summary-output.md): TSV and JSON tabular summaries (`--summary`, `--out-attrs`)

## Input

- [GTF Input Handling](gtf-input.md): GTF attribute parsing, Parent/ID synthesis, and isoform auto-detection
