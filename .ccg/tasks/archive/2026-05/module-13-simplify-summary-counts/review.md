# Review

## Scope

- Simplified `make_summary_row()` feature counting.
- Removed the redundant direct-child counting pass when a gene model is available.
- Did not touch AGAT attribute parsing or CLI output formatting.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: summary counts remain sourced from the gene model when available, and fall back to direct children otherwise.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub query ... --summary-format tsv`
- `gffsub query ... --summary-format json` validated with `python -m json.tool`
- `gffsub query ...` default GFF3 regression
- `gffsub window ... --strand-aware`
- `gffsub qc /home/ww/gffsub/SoyL04.gene.gff`
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
