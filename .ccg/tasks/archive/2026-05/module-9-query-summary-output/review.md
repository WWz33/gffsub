# Review

## Scope

- Added `gffsub query --summary-format tsv|json`.
- Kept default `gffsub query` output as GFF3.
- Added found/not_found summary rows for ID/name/attribute selectors.
- Added per-record child/transcript/exon/CDS summary fields.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: summary output shares query result de-duplication with GFF3 output, so multi-parent descendants are not emitted twice in summaries.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub query ... --summary-format tsv`
- `gffsub query ... --summary-format json` validated with `python -m json.tool`
- Default query output remained GFF3.
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
