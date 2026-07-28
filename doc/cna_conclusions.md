# Reviving Compact Nodal Analysis (CNA) — reasoning

*Written before the plan, per the project's working method: this records why,
the plan records what and in what order.*

## What CNA is

`elements_cna.py`'s `Nullor` docstring names the source, and the paper is now
downloaded and read (not text-extracted — rendered pages via the `Read` tool,
per the project's PDF-verification rule): Esteban Tlelo-Cuautle & Arturo
Sarmiento-Reyes, *A Pure Nodal-Analysis Method Suitable for Analog Circuits
Using Nullors*, Journal of Applied Research and Technology, vol. 1 no. 3,
2003, pp. 235-247 — `~/pycircuit_agy/papers/cna/2003_JART_TleloCuautle_pure-nodal-nullors.pdf`.
A related follow-up by the same group is also saved:
`2005_JART_TleloCuautle_enhancing-symbolic-analysis.pdf` ("Compact system of
equations," avoiding multiplications by zero in determinant evaluation — not
yet read in full; noted for later if Stage 2/3 need it).

**The mechanism is subtler than "the element stamps zero," and matters for
scope.** A nullor is a nullator (input port, `V=0, I=0`) plus a norator
(output port, arbitrary `V, I`) (§2.1, Fig. 1). The paper's base method,
"Pure Nodal Analysis" (PNA, §3, eq. 1), transforms *every* non-NA-compatible
element — including plain independent voltage sources (§2.2, Fig. 2) — into a
nullor-based equivalent, so the assembled matrix contains only R, C, nullors
and independent current sources. This transform-everything approach adds
nodes (the worked inverting-amplifier example, Fig. 8→9, goes from 3 physical
nodes to a 4-node nullor circuit, partly because the input voltage source
itself is turned into a nullor + current source per Fig. 2). "Compacted PNA"
(CPNA, §4, eq. 2) then removes that self-inflicted bloat by four row/column
rules on the assembled `Y_PNA` matrix: merge the two columns of a floating
nullator's terminals (its voltage is 0, so they're the same unknown), merge
the two rows of a floating norator's terminals (conservation makes its
current bookkeeping redundant), or delete the column/row outright when either
is grounded. §6: "the order of the matrix is reduced in one by each nullor."
Two worked, published examples anchor this claim with real numbers: the
CCII-based GIC (Fig. 12/13) goes from a 26×26 `Y_PNA` to a 12×12 `Y_CPNA` with
14 nullors (eq. 9); the lossy Gm-C integrator (Fig. 14/15) goes from 7×7 to
2×2 with 5 nullors, with the closed-form result `Vo/Vin = 1/(1 + sC/gm)`
(eq. 14) — both `26 − 14 = 12` and `7 − 5 = 2` match the headline rule
exactly.

**Open question this plan resolves empirically, not by assumption.**
pycircuit already has proper NA-compatible primitives — `IS` is a genuine
independent current source, so nothing needs the paper's voltage-source→nullor
transform (Fig. 2) that inflates their node count in the first place. A single
4-terminal `Circuit` element stamping an all-zero `G` and declaring no
`branches`, wired directly onto a circuit's *existing* terminal nodes (no
extra internal nodes synthesized), might therefore land at the same reduced
order "for free" — without ever needing the explicit row/column merge
algorithm, simply because it never introduces the artifact nodes that
algorithm exists to remove. That is a hypothesis, not a derivation: hand-
tracing the paper's own intermediate matrix (eq. 4) from a photographed page
risks exactly the transcription error this project's own working method warns
against (rule #10). The sound test is behavioural, matching how this codebase
already validates elements elsewhere (`test_element_jacobians.py`; DDD's
calibration against published `|DDD|` values rather than re-deriving the
recursion): build the paper's own inverting-amplifier example (Fig. 8, using
pycircuit's `IS`/`R` so no source-nullor transform is needed) with the
candidate zero-`G` no-branch element, solve it with pycircuit's existing,
tested numeric machinery, and check the transfer function against the
published closed form `Vo/Vi = -Rf/Ri` (eq. 3). If it matches, the direct
implementation is validated end-to-end without needing the paper's merge
algorithm at all. If it does not, that means the naive element is not
electrically equivalent to a real nullor in this topology, and Stage 1 must
implement something closer to the literal nullator/norator decomposition
instead — a materially bigger feature. Stage 1 in the plan is this test,
stated as a gate before it's run, exactly so a negative result here is
recorded rather than quietly worked around.

The old `elements_cna.py`'s `Nullor.G` (all-zero, no branches, four terminals
mapped directly) is the direct-terminal form, not the paper's literal
nullator+norator decomposition — consistent with treating it as the
hypothesis to test rather than an already-working reference implementation
(it doesn't import under Python 3 either way, see below).

**Hypothesis REFUTED (2026-07-28) — see `cna_plan.md`'s Stage 1 result.** The
direct-terminal, zero-`G`, no-branch element produces a genuinely singular
system for the paper's own example (`det() == 0`, verified symbolically, not
inferred from an exception). It fails to enforce the nullator's constraint
(the old, dead `elements_cna.py` code never actually worked, even on its own
terms) and fails to give the norator's node a current-balancing degree of
freedom. The paper's row/column merge-and-delete rules are load-bearing, not
an artifact of their own methodology that pycircuit's better element set
lets us route around. Getting the real "compact" result requires operating
on node indices at circuit-assembly time — classifying each nullor's ports
as floating or grounded and aliasing/eliminating accordingly before the
matrix is built — which is a materially larger, more architecturally
invasive feature than adding one new element class, and changes the
cost/benefit calculus for continuing. See `cna_plan.md` for the decision
point this raises.

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
`ginac_native.rst`). The paper's own two examples (26→12 for 14 nullors,
7→2 for 5 nullors) are a real published upper bound on what's possible for a
*nullor-dense* circuit — `example10`'s single-nullor MFB filter should show
much less, proportionally, than either.

One more thing from the paper worth carrying forward honestly: its own
conclusion (§6) states the method should be used "only if low accuracy is
required" — nullors are ideal, zero-order models. That's not a new limitation
CNA introduces; the mainline `elements.py` `Nullor` is exactly as ideal
already, so nothing changes about what accuracy trade-off a user is making by
using a nullor at all. It only bears on whether *this specific reduction
technique* is being oversold as a general-purpose accuracy improvement, which
it isn't and doesn't need to be — the target is smaller symbolic determinants
for circuits that already use idealized nullors.

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
