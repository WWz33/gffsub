# Review

## Scope

- Added exact `Dbxref` value support to `AnnotationIndex::find_gene()`.
- Added `Dbxref` to query summary `matched_by` inference.
- Added smoke coverage for a full `Dbxref=GeneID:12345` lookup value.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: lookup is exact on the full Dbxref value. It does not split `GeneID:12345` into `12345`, avoiding ambiguous cross-database matches.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub query ...` default GFF3 regression
- `gffsub query ... --summary-format json` validated with `python -m json.tool`
- `gffsub window ... --strand-aware`
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
