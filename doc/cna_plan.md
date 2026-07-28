# CNA revival — plan

Reasoning and scope: `doc/cna_conclusions.md`. Architecture context: P4 in
`doc/architecture.md`.

## Architect-lens review, done before writing code

*"What does this force to change outside itself, and how do I avoid bloating
those parts?"*

The new element is purely additive — a new class, no change to `circuit.py`,
`toolkit.py`, or the existing `Nullor`/`VCCS`/`CCVS`. It needs no new toolkit
capability: it's a topology choice (which unknowns an element introduces), not
an arithmetic or solve-representation choice, so it sits at the same layer as
any other element in `elements.py` and reaches the solver through the same
`G`/`branches` contract every other linear element uses. The only things that
change outside the new element itself are: one new test file, one new/extended
doc page (following the existing `exec-rst` live-measurement pattern), and the
deletion of two files nothing imports. Nothing about `AC`/`DC`/`Noise`/`DDD`
analysis code changes — this is the same lever noted for DDD in
architecture.md §6: prefer adding at the element/toolkit-primitive layer over
touching analysis code. Passes the lens: the blast radius is one new file plus
deletions, not a design that ripples outward.

## Stages and gates

Each stage states its success condition before running it, per the project's
gate discipline. A negative result at Stage 2 is kept and recorded, not
discarded.

### Stage 1 — the CNA `Nullor`, correctness against the paper's own example

Add a new element (working name `NullorCNA`; naming open for review) to
`elements.py`: `terminals = ('inp', 'inn', 'outp', 'outn')`, no `branches`,
`G` returns an all-zero `(n, n)` matrix — i.e. the direct port of
`elements_cna.py`'s `Nullor.G`, against the current `Circuit` base and
`self.toolkit.zeros` rather than the old file's hardcoded numeric-only zeros.

**Gate, stated before running it — two parts, both against the paper itself
(`doc/cna_conclusions.md`'s "open question" section), not an invented test
circuit:**

1. **Correctness.** Build the paper's own inverting-amplifier example
   (Fig. 8: `Ri` from input to `v-`, `Rf` feedback, `NullorCNA` across
   `(v-, gnd)`–`(vout, gnd)`), driven by pycircuit's `IS` at the input node
   rather than the paper's voltage-source-to-nullor transform (unnecessary
   here — `IS` is already NA-compatible). Solve under `numeric`, and check
   the resulting `Vo/Vin` against the published closed form `-Rf/Ri` (eq. 3).
2. **Compactness.** Check whether the assembled system's dimension for that
   circuit already matches what the paper's Compacted-PNA would give it
   (order reduced by 1 relative to what the mainline branch-based `Nullor`
   needs for the same circuit), *without* implementing the paper's explicit
   row/column merge algorithm — on the hypothesis that pycircuit's existing
   NA-compatible elements make that algorithm unnecessary (see
   `cna_conclusions.md`).

**Success:** both hold — the element is correct and the reduction is free.
**Partial failure (1 holds, 2 doesn't):** the element is still a valid,
correct `Nullor` alternative and worth keeping, but the "compact" framing
needs qualifying — record as a partial negative result and reconsider whether
Stage 2/3 are worth pursuing given the actual size of the win.
**Failure (1 doesn't hold):** the zero-`G` no-branch formulation is not
electrically equivalent in this topology — stop, do not proceed to Stage 2,
and re-derive from the nullor's definition rather than pattern the old code
blindly; record as a full negative result next to the other five in
`[[ddd-implementation]]`'s style.

### Stage 2 — measure the compactness claim

Rebuild `example10`'s MFB filter (or Stage 1's gate circuit, if it's the
better fit) with `NullorCNA` in place of `Nullor`, run it through
`ddd_toolkit` and/or `ginac_toolkit`/`symbolic_poly`, and compare `|DDD|` /
symbolic expression size / solve time against the existing `Nullor` version.
Add this as a live `exec-rst` block in a doc page (new page, or a new section
of `doc/src/circuit/ddd.rst`) so the numbers regenerate at build time rather
than being pasted.

**Gate — declared now, before running it:** a reduction is a win worth
documenting as such; a wash or a regression is a sixth negative result,
recorded next to the other five in `[[ddd-implementation]]`'s style, and
Stage 3 does not proceed without a new reason to.

### Stage 3 (conditional on Stage 2) — `CCVS` via `NullorCNA` + `R`

Only attempted if Stage 2 is positive. Add a `CCVS`-via-nullor composition
(`doc/cna_conclusions.md`'s table: 2 branches → 0 unknowns, the other clear
win from the survey). Same gate shape as Stage 1: correctness against the
mainline `CCVS`, plus a symbolic-size comparison for a circuit that uses one.

### Stage 4 — delete the old files

`circuit_cna.py` and `elements_cna.py` are removed regardless of Stage 2/3's
outcome — nothing in this plan depends on repairing them, and
`doc/cna_conclusions.md` records why they were beyond repair rather than
merely unused. Update `doc/architecture.md`'s P4 entry to `*(fixed)*` with a
pointer to this plan and its outcome (including if Stage 2/3 turned out
negative — say so plainly).

## Testing and doc conventions to follow (already established in this repo)

- Full suite via the scratchpad venv (`[[running-tests]]`); doc build via the
  repo's own `.venv` with sphinx.
- New tests go in `pycircuit/circuit/tests/`, following
  `test_element_jacobians.py`'s pattern of building inside a `SubCircuit` with
  an explicit toolkit rather than depending on `default_toolkit`.
- The explanation for whichever element(s) ship goes in the same commit as the
  code (architecture.md's own rule, and CLAUDE.md's).
