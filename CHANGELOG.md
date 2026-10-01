# Changelog

All notable changes to gffsub are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Region complement (`--exclude-region`) and BED complement (`-b ^FILE`)
- Strand-aware BED overlap (`--same-strand`, `--opposite-strand`)
- Orphan removal (`--drop-orphans`): drops records all of whose Parent
  references leave the kept set, to a fixpoint
- Column-9 projection (`--out-attrs LIST`); ID and Parent are always kept
- Transparent gzip input for files, stdin, FIFOs/devices, and format
  sniffing; concatenated gzip members are read in full
- mmap-backed parsing for regular files
- Per-seqid interval index (sorted starts + prefix max end) for BED overlap
- Flat GTF `--name` lookup falls back to the gene_id attribute when the file
  has no gene feature records
- GTF attribute values may be unquoted; semicolons inside quoted values are
  preserved; GFF3 attributes separated by `; ` are recognized
- BED9/BED12 input detection; zero-length BED intervals map to a 1-based
  point
- Region arguments split on the last colon, so seqids containing `:` work
- `##sequence-region` directives are pruned to surviving seqids;
  `##gff-version` is always the first output line; `##gtf-version` is
  filtered out of GFF3 output

### Fixed
- GTF col9 mixing quoted keys with bare GFF3-style `ID=x;Parent=y` keys
  resolves the bare keys again: the quote-aware GTF parser accepts the
  `key=value` bare form, and ID/Parent are read through it; a bare
  multi-value `Parent=t1,t2` matches per value in `-I`/`-E`/`--grep`
- GTF: `; Key=value` inside a quoted attribute value is no longer read as a
  real attribute; the hierarchy fields (ID/Parent/gene_id/transcript_id) of
  GTF records come from the quote-aware GTF parsers, and `-I`/`-E`/`--grep`
  prefer the record fields on GTF input
- `--longest` saturates the CDS/exon segment sum so transcripts near the
  coordinate limit cannot wrap negative and lose to a shorter isoform
- GFF3 attributes after `; ` keep their `Parent`/`gene_id`/`transcript_id`
  values, so `--format gtf` no longer emits empty `gene_id ""`
- `--longest` sums every CDS segment of a transcript; distinct CDS IDs under
  one transcript are segments of it, not competing variants
- `--longest` applies to files without gene records (flat GTF, GFF3 with
  transcript rows only), grouping isoforms by their parent attribute
- `-I`/`-E`/`--grep` on a multi-value attribute (`Parent=g1,g2`) match any
  value; `!=`/`!~` match only when no value does
- Records keep valid string views when a parsed buffer shorter than the SSO
  limit is moved into the index
- `--up`/`--down` saturate instead of overflowing the window coordinates
- `-r :100-200` (empty seqid) and empty selectors on the `query`/`window`
  subcommands are rejected instead of matching nothing
- `-w VALUE` without `=` reports `--where`, not an unrelated option name
- `-b/--bed` with an unopenable file exits 1 instead of silently dropping
  every record
- `--longest` is deterministic for multi-parent records: winner-priority
  kept assignment, no cross-thread record writes
- `-b/--bed` no longer mis-detects BED9/BED12 as GFF3
- window subcommand accepts the documented `-i/-u/-D/-a` short forms
- `--ids`/`--id-list` are repeatable and accumulate
- Empty values for `--seqid`, `--source`, `--type`, `--name`, `--ids`,
  `--region`, `--bed`, `--longest-type`, `--up`, `--down`, `--sort`,
  `--output`, `--grep-file`, `--grep-field`, `--nearest` are rejected
- Truncated gzip input is rejected on every route (a cut member used to
  parse as a partial file on the regular-file path)
- gzip format sniffing skips comment and blank lines and strips a trailing
  CR, matching the plain-file path
- `--out-attrs` matches keys after the same trim the attribute indexer
  applies
- `--threads` help text and README state the actual default (6)

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
