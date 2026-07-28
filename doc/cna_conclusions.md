# Reviving Compact Nodal Analysis (CNA) — reasoning

*Written before the plan, per the project's working method: this records why,
the plan records what and in what order.*

## What CNA is

`elements_cna.py`'s `Nullor` docstring names the source:

> Esteban Tlelo-Cuautle / Arturo Sarmiento-Reyes, *A Pure Nodal-Analysis
> Method Suitable for Analog Circuits Using Nullors*, Journal of Applied
> Research and Technology, vol. 1, no. 3.

**The paper itself is not available locally** — nothing in `~/pycircuit_agy/papers/`
or elsewhere in this checkout. Per the extraction rule, no number or formula
from it is used here; everything below is re-derived from circuit theory and
cross-checked against the *existing, tested* MNA elements, which is a stronger
standard for this purpose than trusting an unseen citation would be. If the
paper turns up later, treat this as independently re-derived, not "confirmed by
the paper."

The idea, confirmed by reading both the old CNA code and the current MNA code
side by side: MNA gives every element that defines a *voltage* (a source, or a
nullor's zero-input constraint) an extra unknown — a branch current — plus a
row and column enforcing that constraint. CNA's nullor instead stamps an
**all-zero G with no branch at all**: `elements_cna.py`'s `Nullor.G` returns a
zero matrix and declares no `branches`, versus the mainline
`elements.py:722` `Nullor`, which adds one branch and four stamped entries.
Other elements (inductors, voltage sources, controlled sources) are then
*composed* from nullors, gyrators and resistors rather than given their own
branch equations.

## What was measured, and what it changes about the plan

I read the composed forms in `elements_cna.py` and counted unknowns against
the current `elements.py` equivalents:

| element | mainline MNA | CNA composition |
|---|---|---|
| `Nullor` | +1 branch | +0 |
| `CCVS` | +2 branches | +0 (`Nullor` + `R`, no internal nodes) |
| `VCCS` | +0 (already a direct 4×4 stamp) | +2 internal nodes — **worse** |
| `VCVS`, `CCCS` | +1 branch / composed | composed via 2 internal nodes each — no clear win |

This is the reason the plan below is scoped to the bare `Nullor` (and,
optionally, `CCVS`) rather than the whole element zoo `elements_cna.py`
reimplements. CNA is not a universal improvement — the mainline's `VCCS` is
already optimal, and building it through a nullor composition would spend
unknowns to no benefit. The paper's title says "suitable for circuits *using*
nullors" — i.e. circuits that already contain a nullor as a primitive (op-amp
macromodels, nullor-based feedback topologies), not a recipe for rebuilding
every controlled source from scratch. Scoping to that reading is a deliberate
choice, recorded here so it can be reconsidered if wrong.

I also confirmed neither `circuit_cna.py` nor `elements_cna.py` currently
imports under Python 3 (`from constants import *` / `from circuit import *`,
both Python-2 implicit relative imports) — so despite being touched by the
Python-3 port, the NumPy-2.0 fixes and the toolkit refactor (confirmed by
`git log`), neither has been exercised in years. `circuit_cna.py` (1377 of the
1528 lines) is a near-total duplicate of the general `Circuit`/`SubCircuit`
engine — grepped for "nullor" and "transform" (the docstring's own words) and
found nothing CNA-specific in it at all. **Revival means a clean-room
reimplementation of the ~150 lines of real content in `elements_cna.py`
against the current `circuit.py`/toolkit, not a repair of either file.** Both
old files are deleted outright as part of this work, independent of how the
revival itself goes.

## Why this is worth doing at all

Reducing system dimension for nullor-containing circuits is not an abstract
win — it plugs directly into the thing this codebase's symbolic backends
already struggle with. `[[ddd-implementation]]` and `[[cna-jax-vectorization-state]]`
both record that a fully-symbolic determinant is exponential in the number of
unknowns, and that fraction-free/DDD methods only win with *few* generators.
One fewer unknown per nullor is one row and column off every determinant a
nullor-containing circuit produces under `symbolic`, `symbolic_poly`, `ginac`
or `ddd_toolkit`. `example10.py` already builds exactly such a circuit — an
MFB filter with one `Nullor` — for GiNaC/DDD benchmarking, so the plan reuses
it rather than inventing a synthetic test circuit.

This is a **hypothesis to measure, not a given.** One unknown fewer is a small
change for a 2-3 node filter; whether it's large enough to show up as a
meaningful `|DDD|` or symbolic-size reduction on a realistic circuit is an
open question the plan's Stage 2 gate answers empirically, following rule #6
(never type a measured number into prose — an `exec-rst` block regenerates it
at doc-build time, the pattern already used for `soe_symbolic.rst` and
`ginac_native.rst`).

## Scope

**In scope:**
- A new, additive `Nullor` variant using the CNA (zero-`G`, no-branch)
  formulation, built against the current `circuit.py` and toolkit — not
  replacing the existing branch-based `Nullor`, which stays exactly as is
  (it's load-bearing: `example10.py`, the DDD calibration circuits, and
  presumably other examples depend on its current behaviour and matrix shape).
- Correctness validation: circuits built with the new element must match the
  mainline `Nullor`'s DC/AC behaviour to numerical precision.
- A live, build-time measurement of whether the new element actually shrinks
  `|DDD|`/symbolic expression size on `example10`'s MFB filter (or an
  equivalent), reusing the existing `ddd_toolkit`/`ginac_toolkit`.
- Deletion of `circuit_cna.py` and `elements_cna.py`.

**Out of scope, with reconsider-if:**
- `VCCS` built via nullor composition. *Reconsider if:* a future circuit
  specifically needs a `VCCS` that shares internal state with an existing
  nullor network in a way the direct 4×4 stamp cannot express — no such case
  exists today, and the composition is strictly worse by the unknown count.
- `VCVS`, `CCCS`, inductor/voltage-source-via-gyrator composed forms.
  *Reconsider if:* Stage 2's measurement shows the nullor saving is large
  enough on realistic circuits to justify extending the pattern, and a
  concrete circuit needs one of these specifically in CNA form (e.g. a larger
  op-amp macromodel with several nullors and controlled sources together).
- Reviving `circuit_cna.py`/`elements_cna.py` as files, or trying to make them
  importable. *Reconsider if:* never — they duplicate the core engine with no
  CNA-specific content in the duplicated 1377 lines, so repairing them would
  just resurrect a parallel `Circuit` implementation to maintain forever.
- Making the new element the *default* `Nullor`, or switching existing
  examples/tests to it. *Reconsider if:* Stage 2's measurement is strongly
  positive and the existing `Nullor`'s branch-based matrix shape turns out not
  to be depended on anywhere it matters — that would need its own audit, not
  assumed here.
