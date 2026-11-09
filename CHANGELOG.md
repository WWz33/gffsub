# Changelog

All notable changes to gffsub are documented in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `-f gtf25` (GTF2.5 vocabulary) and `-f gtf22` (GTF2.2
  vocabulary: CDS, exon, start_codon, stop_codon, 5UTR, 3UTR, inter,
  inter_CNS, intron_CNS; no gene or transcript rows)
- `--format gtf` writes the GTF3 vocabulary (gene, transcript, exon, CDS,
  start_codon, stop_codon, five_prime_utr, three_prime_utr, Selenocysteine);
  `gtf3` is an alias for it. This is what Ensembl, GENCODE, StringTie,
  featureCounts and ballgown exchange
- GTF headers name the vocabulary: `##gtf-version 3`, `##gtf-version 2.5`,
  `##gtf-version 2.2`
- GTF `gene_id` / `transcript_id` prefer the column-9 attribute of that name
  and fall back to `ID` / `Parent`, so Ensembl GFF3 `gene:` / `transcript:`
  prefixes no longer reach the GTF output
- GTF output trims the terminal CDS at its `stop_codon` when the codon abuts
  the CDS 3' end of the same transcript
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

### Changed
- BED input requires `-f bed` or `-s`; other output formats exit 1
- GTF output keeps `ID` and `Parent` as extra quoted attributes from GFF3 input
- A quoted GTF `Parent "t1,t2"` links to both parents in the attribute index,
  `--drop-orphans`, `--longest` and GTF output; `--where Parent=t1` matches the
  record and GTF output writes one row per parent
- BED output writes `0` for non-finite scores (`inf`, `nan`)
- `-o` naming the input file is rejected; Linux faulted on the mmap'd pages
  after the output open truncated the same file
- A directory path is rejected instead of parsing as an empty annotation
- Bytes after the last gzip member are ignored on stdin and FIFO, matching the
  file path
- Missing `<cmath>` and `<vector>` includes are added so GCC/libstdc++ builds

### Removed
- BED to GFF3 and BED to GTF conversion
- JSON summary output (`--summary-format json`); the TSV form remains
- `CMakeLists.txt`; the Makefile is the only build

### Fixed

#### Architecture review 2026-11-10
- One column-9 decoder (`decode_attributes`) replaces six scanners. `-w`,
  `-I`, `--grep`, `--name` and the GTF writer now agree on a value:
  `Parent=t1, t2`, repeated GTF keys (`tag "a"; tag "b";`), and `%2C`
  identities resolve the same way everywhere
- A flat GTF keeps `gene_id` when converted to GFF3 (child rows carry
  `gene_id=` instead of losing the gene)
- `--name` reads a GTF `Name` attribute
- A GFF3 value containing a comma outside the five list tags is no longer
  truncated on GTF conversion
- `query()` and `window()` results print every matched record: a `kept` flag
  left by an earlier `subset()` no longer hides them
- Parsing into a non-empty `GffData` clears the previous file first; record
  views no longer point at the old mapping
- `AnnotationIndex::from_data(GffData&&)` is deleted, so a temporary cannot
  leave a dangling index
- `copy_storage_from` copies `source_format`; `clear()` resets it
- `-o -` writes to stdout instead of creating a file named `-`; `-o x.gz`
  is rejected
- `query` and `window` accept `-o`, `-f`, and `--out-attrs`
- `query` output keeps the input's `##` directives, matching the subset path
- `--help` writes to stdout; `query --version` and `window --version` work
- `-i=ID`, `-t=type`, and `-o=FILE` are accepted like their `--long=value`
  forms
- An empty `--ids`/`--grep-file` file, or a directory passed to either,
  is an error instead of silently disabling the filter
- A 3-column BED with `--same-strand`/`--opposite-strand` is an error
  instead of an empty or unfiltered result
- BED output replaces tab, newline, and CR in the name with a space
- A multi-parent GTF expansion writes a unique ID per row (`<id>.<parent>`)
- `stop_codon` trimming requires the same seqid and strand
- GTF output falls back to `transcript_id` when a transcript has no gene,
  instead of writing `gene_id "";`
- GFF3 output maps `5UTR`/`3UTR` to `five_prime_UTR`/`three_prime_UTR`
- `--out-attrs` keeps `ID` and `Parent` in both output formats
- GTF output uses canonical type capitalization (`CDS`, `Selenocysteine`)
- `--longest` uses one shared child index; the index build is linear again
  (80k flat records: 0.90s to 0.16s)
- `print_gtf` takes a `GtfModel` built once from the full annotation, so a
  library caller printing a query subset gets the same `gene_id` and CDS
  trim as the CLI
- Format detection reads the first column-9 fragment instead of searching the
  whole column: `ID=g1;Name="A"` stays GFF3 (`tag=value`) and
  `note "a=b"` stays GTF (`tag "value"`), so a quoted GFF3 value no longer
  sends the file through the GTF parser
- GTF: a `#` outside a quoted value starts a comment instead of becoming an
  attribute; `\t`, `\n` and `\r` are decoded so a value round-trips
  unchanged; `key="value"` keeps `key` as the tag name
- GFF3: the comma list is split only for Parent, Alias, Note, Dbxref and
  Ontology_term; a comma in any other tag is part of the value
- GFF3 output always declares `##gff-version 3`; an input line claiming
  another version is not reproduced
- GTF output writes one row per parent for a feature with several parents
- Parsing reports what it skipped or could not check, once per file: lines
  with too few columns, lines with invalid coordinates, CDS lines whose phase
  is not 0/1/2, and Parent links that form a cycle
- GTF col9 mixing quoted keys with bare GFF3-style `ID=x;Parent=y` keys
  resolves the bare keys again: the quote-aware GTF parser accepts the
  `key=value` bare form, and ID/Parent are read through it; a bare
  multi-value `Parent=t1,t2` matches per value in `-I`/`-E`/`--grep`
- GTF: `; Key=value` inside a quoted attribute value is no longer read as a
  real attribute; the hierarchy fields (ID/Parent/gene_id/transcript_id) of
  GTF records come from the quote-aware GTF parsers, and `-I`/`-E`/`--grep`
  prefer the record fields on GTF input
- `--longest` saturates the per-CDS and per-exon segment sums so transcripts
  near the coordinate limit cannot wrap negative and lose to a shorter isoform
- `--longest` scores a transcript by its longest CDS. Lines sharing an ID are
  one discontinuous CDS and are summed; distinct CDS IDs under one transcript
  are alternative products (the spec's alternative start codons), not
  segments, per the GFF3 specification's canonical gene example
- `--longest` applies to files without gene records (flat GTF, GFF3 with
  transcript rows only), grouping isoforms by their parent attribute
- GFF3 attributes after `; ` keep their `Parent`/`gene_id`/`transcript_id`
  values, so `--format gtf` no longer emits empty `gene_id ""`
- A feature ID reused by lines that differ in type, seqid or strand warns:
  the spec requires unique IDs, and lineage lookups or GTF `gene_id`
  assignment may otherwise pick the wrong line
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
- `--threads` help text states the actual default (6)

### Fixed (defect review 2026-11-09)

Details, commands, and before/after output for each entry are in
`docs/defect-review.md`.

- A GTF transcript carrying a `Parent` attribute keeps its `gene_id` link, so
  `--longest` no longer drops the gene as childless; a transcript is never
  listed as its own parent
- GTF row expansion for a multi-parent child writes each parent's transcript
  on its own row: an existing `transcript_id` no longer overwrites the row's
  parent, and the repeated Parent list is dropped, so re-reading the output
  keeps the row count instead of doubling it
- A `stop_codon` shared by several transcripts trims each parent's terminal
  CDS, not only the first
- `--parents` and `--model` expand an ID-less child through its own Parent
  links instead of returning the child alone
- A quoted GTF `Parent` list is URL-decoded, matching the GFF3 path
- GTF to GFF3 output writes the full `Parent` list, not just the first value
- The Parent cycle warning walks every Parent of every line, with no depth
  cutoff: a cycle reachable only through a non-first parent is reported
- A Parent that names no record in an input GFF3 file is reported once
- An empty BED file is still recognized as BED, so the conversion guard fires;
  an empty file takes its format from the `.bed` / `.gtf` suffix
- BED `score` is written as an integer clamped to 0..1000; fractional,
  negative, out-of-range, unparseable and non-finite values write `0`
- A write failure after a successful open (full disk, closed stdout) exits 1
  instead of 0 with a truncated file
- A closed stdin exits 1 instead of printing an empty annotation
- A UTF-8 BOM is stripped before parsing and format sniffing, so it no longer
  becomes part of the first seqid or hides a `##` directive
- GFF3/GTF rows with more than nine columns (BED: more than 12) are reported
  and skipped instead of silently dropping the extra fields
- Multi-character strand and phase values are rejected; the row is skipped
  instead of truncating `++` to `+`
- `parse_list` trims each item, so `-t "mRNA, exon"` keeps both types
- Summary lengths, sums, coverage and the median saturate at the int64 limit
  instead of printing negative values
- `-t`/`-@` validation: a whitespace-prefixed negative thread count is
  rejected like an unspaced one
- A second positional operand is a usage error instead of an ignored file
- `query` and `window` accept `--flag=value` like the default CLI
- Input errors keep the operating system's reason (`No such file or
  directory`, `Permission denied`, `is a directory`)
- `--longest` dedups isoform and child indices with a set: a gene with tens of
  thousands of isoforms no longer costs quadratic time
- `docs/library-api.md`: field names match the headers (`type`, `filter_by_type`,
  `std::string_view` record fields), `infer_input_format` is attributed to
  `input_source.hpp`, `from_owned` / `find_all_genes` / `parents_of_record` are
  documented, the borrowed-index lifetime requirement is stated, and the link
  command includes `-lz`
- `docs/output-formats.md`: the GTF examples show the `ID` / `Parent`
  attributes the writer emits

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
- Summary output in TSV (`--summary`); the JSON formatter was removed before the first release
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
