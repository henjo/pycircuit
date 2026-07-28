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

### Stage 1 — the CNA `Nullor`, correctness only

Add a new element (working name `NullorCNA`; naming open for review) to
`elements.py`: `terminals = ('inp', 'inn', 'outp', 'outn')`, no `branches`,
`G` returns an all-zero `(n, n)` matrix — i.e. the direct port of
`elements_cna.py`'s `Nullor.G`, against the current `Circuit` base and
`self.toolkit.zeros` rather than the old file's hardcoded numeric-only zeros.

**Gate:** build a small feedback circuit two ways — once with the mainline
`Nullor`, once with `NullorCNA` — solve both under `DC` and `AC` with the
`numeric` toolkit, and assert node voltages / transfer function agree to
numerical precision. A unity-gain buffer or the existing MFB topology with a
resistive load both work; pick whichever composes more simply with what
`example10.py` already builds, so Stage 2 can reuse the same circuit.
**Success:** matches. **Failure:** the zero-`G` formulation is not
electrically equivalent for some topology — stop and re-derive from the
nullor's definition (nullator/norator, i.e. `v=0` at the input port, arbitrary
current; arbitrary voltage at the output port, `i=0`) rather than pattern the
old code blindly.

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
