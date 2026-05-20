# Review

## Scope

- Added read-only `gffsub qc <input.gff3>` command.
- Reports duplicate IDs, invalid ranges, missing parents, and child coordinates outside parent coordinates.
- Outputs TSV to stdout and returns non-zero only for argument or parse errors.

## Findings

- Critical: none found.
- Warning: external CCG Gemini/Claude wrapper was unavailable in this environment, so this review records local inspection plus WSL verification instead of dual-model output.
- Info: a Windows-created CRLF test fixture showed carriage returns in raw attribute-derived values. Parser newline normalization is a broader concern and was not changed in this QC-only module.

## Verification

WSL with `/home/ww/miniforge3/bin/mamba`:

- `cmake -S /mnt/f/codex/gffsub -B "$build"`
- `cmake --build "$build" -j2`
- `ctest --test-dir "$build" --output-on-failure`
- `"$build/annotation_index_smoke" /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub qc /home/ww/gffsub/SoyL04.gene.gff`
- `gffsub qc` against a synthetic bad GFF3 covering duplicate ID, invalid range, missing Parent, and child outside Parent.
- `gffsub query ... --summary-format json` validated with `python -m json.tool`
- `gffsub window ... --strand-aware`
- Existing CLI regression on SoyL04:
  - region gene output: 2 lines
  - region gtf3 output: 9 lines
  - longest output: 782042 lines
