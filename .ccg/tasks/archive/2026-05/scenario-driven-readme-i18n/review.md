# Review

## Scope

- Rewrote `README.md` as the English source-of-truth for a scenario-driven bioinformatics CLI README.
- Rewrote `README.zh.md` as the Chinese sibling with matching section hierarchy and preserved command/code mechanics.
- Did not touch the unrelated `src/gffsub.cpp` working-tree change.

## External review

- Gemini analysis/review attempts failed because Gemini CLI requires interactive OAuth in this environment.
- Claude analysis completed and recommended the scenario-first rewrite, supported input/output section, scannable QC matrix, and reduced duplicated command shortcut material.
- Claude review completed and identified issues that were fixed: final newline, Chinese `flag` cells, `--summary`/`--out-attrs` query-style constraint, compatibility subcommand notes, and `gtf` alias visibility.
- A focused review agent found one additional usability warning: a quick-start row mixed full commands and option fragments. Fixed by splitting score, strand, and phase into separate copy-pasteable rows.

## Verification

- `git diff --check -- README.md README.zh.md` passed.
- Each README contains exactly one `README-I18N` selector block.
- Both READMEs contain 24 code-fence delimiters / 12 fenced code blocks.
- Fenced code block contents match byte-for-byte between `README.md` and `README.zh.md`.
- Heading levels match between `README.md` and `README.zh.md`.
- Same-file anchor count is zero in both READMEs, so no localized anchor rewrites were required.
- Both README files end with LF.

## Spec evolution

No `.ccg/spec` update. This task changed README organization and localization only; it did not add a new reusable project coding or GFF3 semantics convention.