# Review

## Scope

- Added `window_region()` for genomic and strand-aware upstream/downstream windows.
- Added `gffsub window <input.gff3> --id ID --upstream N --downstream N [--strand-aware]`.
- Added smoke coverage for plus- and minus-strand window coordinate behavior.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: this module only outputs GFF3 for `window`; existing CLI and `query` output semantics were not changed.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub window ... --upstream 100 --downstream 50`
- `gffsub window ... --upstream 100 --downstream 50 --strand-aware`
- `gffsub query ...` default GFF3 regression
- `gffsub query ... --summary-format json` validated with `python -m json.tool`
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
