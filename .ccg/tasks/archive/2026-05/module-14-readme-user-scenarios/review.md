# Review

## Scope

- Rewrote `README.md` around user-driven scenarios.
- Rewrote `README.zh.md` with the same structure and language selector.
- Covered current commands: classic region extraction, `query`, summary output, `window`, `--longest`, `qc`, output formats, and C++ API.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: README variants each keep exactly one `README-I18N` selector block and have matching code fence counts.

## Verification

- `README.md` selector count: 1 start / 1 end.
- `README.zh.md` selector count: 1 start / 1 end.
- Both README variants have 18 fence markers and no same-file anchor links.
- WSL with `/home/ww/miniforge3/bin/mamba`:
  - `cmake -S /mnt/f/codex/gffsub -B "$build"`
  - `cmake --build "$build" -j2`
  - `ctest --test-dir "$build" --output-on-failure`
  - `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
  - `gffsub query ... --summary-format json` validated with `python -m json.tool`
  - `gffsub window ... --strand-aware`
  - `gffsub qc /home/ww/gffsub/SoyL04.gene.gff`
  - Existing CLI regression on SoyL04:
    - region gene output: 2 lines
    - region gtf3 output: 9 lines
    - longest output: 782042 lines
