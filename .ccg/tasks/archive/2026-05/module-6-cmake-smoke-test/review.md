# Review

## Scope

- Added a CTest target for `annotation_index_smoke`.
- Made the smoke test self-contained when no annotation path is provided.
- Kept optional SoyL04 real-data mode with `annotation_index_smoke <annotation.gff3>`.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in both Windows and WSL checks, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: the self-contained descendant assertion follows the current graph-edge behavior, where a multi-parent child can be visited through each parent path.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
