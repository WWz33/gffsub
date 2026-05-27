# Review

- Scope: `src/gffsub.cpp` only.
- Change: extracted repeated percent-encoding QC output into `print_percent_encoding_qc`.
- Risk: low; emitted error code and message text remain unchanged.
- Spec update: none. This did not change GFF3 semantics or add a reusable project convention.

Verification:
- `cmake --build build` could not run because `cmake` is not on PATH.
- `mingw32-make test` compiled successfully and passed `annotation_index_smoke`, then failed at `cli_output_attrs_smoke` because the Makefile invokes `./gffsub`, which Windows `cmd` rejects.
- Manual Windows-style smoke tests passed:
  - `.\annotation_index_smoke.exe`
  - `.\cli_output_attrs_smoke.exe .\gffsub.exe`
  - `.\cli_selector_smoke.exe .\gffsub.exe`
