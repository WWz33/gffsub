# Requirements

- Refactor the repository README around scenario-driven bioinformatics workflows.
- Keep `README.md` as the English source-of-truth.
- Keep the existing localized sibling `README.zh.md` and update it consistently.
- Preserve GitHub Markdown mechanics: selector block, command literals, code fences, tables, file links, and option names.
- Keep project identifiers, CLI flags, paths, and commands unchanged except where examples are intentionally reordered or grouped.
- Do not touch the existing unstaged C++ cleanup in `src/gffsub.cpp`.

## Notes

- Existing target language is Chinese via `README.zh.md`.
- The configured CCG `codeagent-wrapper` is not present in this environment, so dual-model analysis/review will use the available OMX team runner with `gemini` and `claude` workers.
