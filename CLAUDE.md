# pycircuit

## graphify knowledge graph

This project has a knowledge graph at `graphify-out/` (gitignored) with god nodes,
community structure, and cross-file relationships. Rebuild it with `graphify update .`.

Rules:
- For codebase questions, first run `graphify query "<question>"` when
  `graphify-out/graph.json` exists. Use `graphify path "<A>" "<B>"` for relationships,
  `graphify explain "<concept>"` for focused concepts, `graphify affected "<X>"` to find
  what a change impacts, and `graphify god-nodes --top 10` for architectural hubs.
  These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- Use the results to decide which files to read.
- Read `graphify-out/GRAPH_REPORT.md` only for broad architecture review, or when
  query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost).
  Community names are rebuilt as hub names by `update`; restore them with
  `graphify label . --backend claude-cli`.
