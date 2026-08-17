# Changelog

All notable changes to gffsub are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2026-08-16

First tracked release. gffsub is a C++17 command-line tool for subsetting GFF3
and GTF annotation files by region, feature ID, attribute, or gene model.

### Added
- Region filtering (`-r`, `--region`) with 1-based inclusive coordinates
- BED interval overlap filtering (`-b`, `--bed`) using 0-based half-open semantics
- Feature ID lookup (`--id`, `--ids`) with exact match
- Gene name lookup (`--name`) across ID, Name, gene_id, locus_tag, Alias, Dbxref
- Attribute filtering (`--where KEY=VALUE`)
- Grep filtering (`--grep`, `--grep-regex`, `--grep-file`) with case-insensitive mode
- Boolean expression filtering (`-I`, `--include-expr`, `-E`, `--exclude-expr`)
- Gene model expansion (`--children`, `--parents`, `--model`)
- Nearest gene query (`--nearest`)
- Window subcommand and shortcut (`--up`, `--down`, `--strand-aware`)
- Feature type filter (`-f`, `--feature`)
- Column filters: `--seqid`, `--source`, `--score`, `--strand`, `--phase`
- Longest isoform selection (`-L`, `--longest`) with multi-threading (`-@`, `--threads`)
- Output format conversion: GFF3, GTF2, GTF3 (Ensembl), BED (`-t`, `--format`)
- Summary output in TSV or JSON (`--summary`)
- Selected attribute projection (`--out-attrs`)
- `query` and `window` subcommands for query-style API access
- Input format auto-detection by content sniffing (GFF3, GTF, BED)
- `--version` flag
- GTF3 type whitelist with case-insensitive matching
- GTF type label normalization (mRNA to transcript, UTR spelling)
- GTF parent ID synthesis from gene_id / transcript_id
- URL decoding for GFF3 attribute values
- Coordinate validation (start >= 1, start <= end)
- `##directive` preservation in output headers
- Control character escaping in GTF output
- Feature type classification module (`feature_types`) with `FeatureClass` enum
- Shared string utilities module (`string_utils`)
- Version header (`version.hpp`) for single source of version truth
- MIT License

### Changed
- Split monolithic `gffsub.cpp` into library + CLI architecture (lib + thin main)
- Extracted type classification from `record.hpp` into independent `feature_types` module
- Moved GTF type label mapping and GTF3 whitelist from `annotation_output.cpp`
- Consolidated duplicate `to_lower` / `lowercase_copy` into shared `string_utils`
- Consolidated duplicate `split_line` into shared `string_utils`
- Consolidated test helpers from 3 test files into `tests/test_utils.hpp`
- Unified version number across CMakeLists.txt, help text, and `--version` flag
- README simplified to behavior-only documentation

### Fixed
- GFF3 to GTF conversion preserves `##directive` headers
- GTF output escapes control characters in attribute values
- TSV summary escaping for special characters
- Window subcommand handles discontinuous features correctly
- GTF to GFF3 conversion applies URL escaping to comma-containing values
- GTF3 whitelist matches case-insensitively (GFF3 `five_prime_UTR` no longer dropped)
- UTR normalization runs before GTF3 emittable filter
- Unsorted GTF input does not cause ID collision between exon and transcript
- Multi-parent exon handling in longest isoform selection
- Coordinate validation rejects start=0 and start > end
- Input format inference recognizes `.bed` extension
- GTF format inference handles malformed lines gracefully

[0.1.0]: https://github.com/WWz33/gffsub/releases/tag/v0.1.0
