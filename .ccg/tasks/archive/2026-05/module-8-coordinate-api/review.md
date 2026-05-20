# Review

## Scope

- Added explicit `BedRegion` type for 0-based half-open BED intervals.
- Added `to_bed_region()` and `from_bed_region()` conversion APIs.
- Added smoke coverage for GFF 1-based inclusive to BED 0-based half-open round-trip conversion.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: existing CLI BED output behavior was not changed.

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
