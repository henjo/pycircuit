# pycircuit

## Knowledge graph (graphify)

`graphify-out/graph.json` holds a knowledge graph of this repo (code-only AST
extraction, ~9.9k nodes). It is gitignored — rebuild with `graphify update .`.

Before grepping or reading files to answer a question about the architecture,
what calls what, or where something lives, **query the graph first**:

```bash
graphify query "<question>"        # scoped subgraph for a plain-language question
graphify explain "Transient"       # one node and its neighbors
graphify path "DC" "SubCircuit"    # shortest path between two concepts
graphify affected "Branch"         # what breaks if this changes
graphify god-nodes --top 10        # architectural hubs
```

Run these from the repo root. Use the results to pick which files to read.
