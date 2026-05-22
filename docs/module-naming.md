# Module Naming Guide

This document defines the naming direction for `gffsub` modules before any file split or rename. The goal is to avoid a half-migrated codebase where strict GFF3 semantics, GTF3 output, annotation graph code, and command names all mean different things.

## Source References

- AGAT GxF guide: https://github.com/NBISweden/AGAT/blob/master/docs/gxf.md
- Sequence Ontology GFF3 specification: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md

Project scope: focus on strict GFF3 semantics and GTF3 output. Use AGAT's GxF guide only to name and separate format flavors clearly, not to expand the project into a full GFF/GTF family toolkit. Use the Sequence Ontology GFF3 specification when changing strict GFF3 semantics, especially column 9 attributes, coordinates, phase, and parent-child feature rules.

## Naming Principles

| Term | Meaning | Use for |
|------|---------|---------|
| `gff3` | Strict GFF3 flavor and input semantics | Parsing, attributes, parent-child graph semantics, coordinates, phase, and user-facing GFF3 behavior. |
| `gtf3` | GTF3 output flavor | Output format selection and GTF3-specific output rules. |
| `annotation` | Biological annotation data model | Records, indexes, feature graph, gene models, lookup APIs. |
| `attribute` | GFF3 column 9 key/value semantics | Strict GFF3 `tag=value` attributes. Keep GTF-style quoted attributes separate if needed for compatibility. |
| `region` | Coordinates and intervals | Region strings, GFF3 coordinate conversion, windows, and any compatibility conversion helpers. |
| `filter` | Keep/drop operations on records | Region/type filters and other selection filters. |
| `isoform` | Transcript selection rules | Longest transcript and future representative transcript logic. |
| `query` | User-facing lookup commands | CLI command orchestration and result assembly. |
| `summary` | Structured non-GFF output | TSV/JSON summaries and selected attribute projection. |

Compatibility terms:

- `bed` may remain as a format flavor for region input/output compatibility, but it should not drive module names.
- `gtf2` may remain as a legacy output option if the CLI keeps supporting it, but new naming should not center it.
- Avoid broad `gxf_*` names unless a module genuinely owns shared flavor-neutral formatting code for GFF3 and GTF3.

## Compatibility Policy

Compatibility flavors are allowed, but they do not define new module boundaries.

| Surface | Policy | Naming impact |
|---------|--------|---------------|
| `InputFormat::GFF3` | Core input semantics. | May drive parser and attribute module names. |
| `InputFormat::GTF` | Compatibility shim only. | Do not create broad `gtf_parser` modules unless GTF input becomes a core scope decision. |
| `InputFormat::BED` | Compatibility for region-style input. | Keep BED naming local to region conversion or I/O helpers. |
| `OutputFormat::GFF3` | Core output flavor. | May drive explicit GFF3 output helpers. |
| `OutputFormat::GTF3` | Core output flavor. | May drive explicit GTF3 output helpers. |
| `OutputFormat::GTF2` | Legacy compatibility output. | Do not use as a module prefix. |
| `OutputFormat::BED` | Compatibility output. | Do not use as a module prefix. |

If a compatibility flavor starts requiring substantial independent behavior, revisit project scope before adding new modules.

## GTF3 Output Policy

GTF3 is a core output flavor, not a core input semantics layer. The explicit `print_gtf3()` API owns the core GTF3 output path, while `print_gtf()` remains a compatibility dispatch entry point for callers that still pass `OutputFormat`.

GTF2 support, if kept, should remain a legacy branch inside output compatibility code unless the project explicitly expands scope.

## Output Split Rule

Use a single output module only while GFF3 and GTF3 output share most implementation. Split output modules when any of these become true:

- GFF3 and GTF3 need different attribute builders beyond small flavor branches.
- GTF3 feature-type filtering or transcript naming grows beyond local output logic.
- Compatibility output such as GTF2 or BED starts obscuring GFF3/GTF3 core behavior.

Preferred split, if needed:

```text
gff3_output.cpp       # strict GFF3 writer
gtf3_output.cpp       # GTF3 writer
compat_output.cpp     # optional legacy GTF2/BED output helpers, if still supported
```

## Public CLI Philosophy

The public CLI should follow the `bcftools`-style `gffsub <file> [selectors] [modifiers] [output]`
shape where possible. Concepts such as ID lookup, attribute selection, region selection, and
upstream/downstream windows are common GFF annotation semantics, not user-facing module names.

Internal source modules may use names such as `region`, `query`, `summary`, or `annotation_filter`
to keep the code auditable. Those names do not automatically justify public subcommands.

Public subcommands should be kept only when the output contract is genuinely different from
record subsetting. `qc` may remain a candidate because it emits diagnostic rows, not annotation
records. `query` and `window` should be treated as compatibility or transition surfaces until
top-level selectors cover the same semantics.

## Query And Summary Extraction Policy

`query` and `summary` are first-class internal concepts, but they can remain inside `gffsub.cpp`
while public CLI naming is still being normalized. Do not extract them into new source modules
until the user-facing selector model is settled.

Suggested threshold for extraction:

- Add `query_command.cpp` when query option parsing, matching, and result assembly are no longer easy to audit together.
- Add `summary_output.cpp` when TSV/JSON summary formatting or `--attrs` projection grows beyond simple row printing.
- Keep `gffsub.cpp` as the command dispatcher and top-level option parser.

## Current Naming Audit

| Current name | Status | Issue | Preferred direction |
|--------------|--------|-------|---------------------|
| `src/annotation.hpp` | Aligned | Public annotation data model and `AnnotationIndex` API. | Keep as the annotation-facing public API surface. |
| `src/gff3.hpp` | Acceptable for now | Compatibility umbrella header for existing users plus strict GFF3 helpers and format utilities. | Keep including `annotation.hpp`; split narrower format headers only when public API migration is planned. |
| `src/attributes.cpp` | Aligned | Owns strict GFF3 `tag=value` attribute projection helper. | Keep GTF quoted attribute parsing separate from this module. |
| `src/gff3_parser.cpp` | Aligned | Owns file and line parsing for strict GFF3 input plus compatibility GTF/BED parsing paths. | Keep parser-only behavior here. |
| `src/region.cpp` | Aligned | Owns region strings, GFF/BED coordinate conversion, and strand-aware windows. | Keep coordinate and window helpers here. |
| `src/annotation_filter.cpp` | Aligned | Owns keep/drop filters by region, region file, and feature type. | Keep selection filters here. |
| `src/annotation_index.cpp` | Aligned | Implements `AnnotationIndex`, feature graph, lookup, overlap, nearest gene, and gene model logic. | Keep as the graph/query core. |
| `src/isoform_filter.cpp` | Aligned | Implements longest isoform selection. | Keep for transcript selection rules; add representative isoform logic here only if it gets a distinct rule. |
| `src/annotation_output.cpp` | Acceptable for now | Writes GFF3, GTF2, GTF3, and BED from annotation records while GFF3 and GTF3 output still share most implementation. | Keep as one module until the output split rule above is triggered. |
| `src/gffsub.cpp` | Acceptable for now | CLI entry point also contains query, summary, window, and QC command bodies. | Keep initially; extract `query_command.cpp` or `summary_output.cpp` only when thresholds above are met. |
| `tests/annotation_index_smoke.cpp` | Aligned | Test name matches API under test. | Keep. |
| `CMakeLists.txt` / `Makefile` source lists | Mechanically tied to old names | They will need updates during file renames. | Update only in the same commits as source file renames. |

## Current API Naming Audit

| Current API | Status | Notes |
|-------------|--------|-------|
| `AnnotationIndex` | Aligned | Good domain name. It should remain independent of file format flavor. |
| `AnnotationIndex::from_gff3()` | Aligned | Keep. The project focuses on GFF3 input semantics. |
| `GffRecord` | Acceptable | It represents parsed GFF3 annotation records. Consider `AnnotationRecord` only if public API naming is broadened later. |
| `GffData` | Partially aligned | Same issue as `GffRecord`; acceptable until public API naming is revised. |
| `InputFormat::{GFF3,GTF,BED}` | Partially aligned | `GFF3` is core. `GTF` and `BED` are compatibility names under the compatibility policy. |
| `OutputFormat::{GFF3,GTF2,GTF3,BED}` | Partially aligned | `GFF3` and `GTF3` are core. `GTF2` and `BED` are compatibility output flavors under the compatibility policy. |
| `parse_attributes()` | Aligned | Shared strict GFF3 attribute helper now lives in `attributes.cpp`. GTF quoted attributes stay separate. |
| `parse_region()` | Aligned | Region naming is format-neutral. |
| `to_bed_region()` / `from_bed_region()` | Aligned | Names correctly describe coordinate conversion. |
| `window_region()` | Aligned | Name fits region/window semantics. |
| `filter_by_region()` / `filter_by_feature()` | Aligned | Filter names are clear; they should live in a filter module. |
| `filter_longest_isoform()` | Aligned | Explicit API for longest transcript isoform selection. |
| `filter_longest()` | Compatibility | Legacy wrapper kept to avoid breaking existing C++ callers. |
| `print_gff3()` / `print_gtf3()` / `print_bed()` | Aligned | Core output helpers are explicit for GFF3, GTF3, and BED compatibility output. |
| `print_gtf()` | Compatibility | Legacy dispatch entry point kept for GTF2/GTF3 compatibility callers. |

## Recommended Target Layout

Do not apply this all at once. This is the target vocabulary for staged migration.

```text
src/
  annotation.hpp          # public annotation data model and index API, if the public header is split later
  gff3.hpp                # strict GFF3 flavor helpers, if format helpers are split later
  gff3_parser.cpp         # file and line parsing for strict GFF3 input
  attributes.cpp          # strict GFF3 tag=value attribute helpers
  region.cpp              # region strings, coordinate conversion, windows
  annotation_filter.cpp   # region/type filters
  isoform_filter.cpp      # longest/representative transcript rules
  annotation_index.cpp    # AnnotationIndex and graph/query API
  annotation_output.cpp   # GFF3 and GTF3 writers, plus any compatibility output
  gffsub.cpp              # CLI entry point and command dispatch
```

## Migration Plan

1. **Document vocabulary first**
   - Keep this document as the source of truth for names.
   - Do not rename files until the target names are agreed.

2. **Rename without moving behavior**
   - Rename one module at a time.
   - Update `CMakeLists.txt`, `Makefile`, and includes in the same commit.
   - Run the smoke test after every rename.

3. **Finish naming before more splits**
   - Do not split new CLI helper modules while public command semantics are still being normalized.
   - Keep `query`, `window`, and summary helpers in `gffsub.cpp` until top-level selectors are settled.
   - Prefer compatibility aliases over behavior changes when introducing clearer names.

4. **Split internal boundaries only after names are stable**
   - `region.cpp` now owns `parse_region()`, BED/GFF conversion, and `window_region()`.
   - `annotation_filter.cpp` now owns `filter_by_region()`, `filter_by_regions_from_file()`, and `filter_by_feature()`.
   - `attributes.cpp` now owns shared strict GFF3 attribute parsing.
   - Keep GTF quoted attribute parsing separate from strict GFF3 attributes, and only keep it if compatibility requires it.

5. **Clean public CLI surfaces last**
   - Make top-level selectors the canonical user-facing path for common GFF operations.
   - Keep or remove subcommands based on output contract, not internal source module names.
   - Before removing a public subcommand, provide an equivalent top-level command and test output equivalence.

6. **Delay public API renames**
   - Avoid renaming `GffRecord`, `GffData`, or `from_gff3()` until file/module naming has settled.
   - If public API names change, provide a clear migration note.

## Naming Checkpoint

The first naming pass is considered complete for internal modules:

- `annotation.hpp` owns the public annotation model and `AnnotationIndex`.
- `attributes.cpp` owns strict GFF3 `tag=value` attribute helpers.
- `gff3_parser.cpp` owns file and line parsing.
- `region.cpp` owns region strings, coordinate conversion, and window regions.
- `annotation_filter.cpp` owns keep/drop filters.
- `isoform_filter.cpp` owns longest isoform selection.
- `annotation_index.cpp` owns graph/query data structures.
- `annotation_output.cpp` owns output helpers while GFF3/GTF3 output remains shared.

The remaining non-aligned names are intentionally deferred:

- `GffRecord` and `GffData` remain stable public API names for now.
- `InputFormat::{GTF,BED}` and `OutputFormat::{GTF2,BED}` remain compatibility enum values.
- `query` and `window` remain transition CLI surfaces until top-level selectors cover the same semantics.
- `print_gtf()` and `filter_longest()` remain compatibility wrappers for existing C++ callers.

Do not start further source boundary splits until this checkpoint is still accurate after the next CLI
normalization pass.

## Decisions For Future Work

- Treat `GxF` as an AGAT reference term, not as the default project module prefix.
- Treat `GFF3` and `GTF3` as the project focus.
- Treat `GTF2` and `BED` as compatibility flavors, not broad module names.
- Treat `annotation` as the project term for graph/index/query biological semantics.
- Keep `gffsub` as the CLI binary name. It expresses the project philosophy even if internal modules use broader GxF naming.
