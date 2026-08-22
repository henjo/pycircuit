# idtmod: research and implementation plan

Research report on implementing the Verilog-A `idtmod` (circular / modulo / periodic
integrator) operator properly in pycircuit. Covers the LRM semantics, a survey of how
ngspice, Xyce, VACASK/OpenVAF and gnucap handle (or fail to handle) the operator, how
general DAE/ODE tooling treats wrapped integrated states, an audit of pycircuit's current
`Idtmod`, and a concrete three-phase implementation recommendation validated against the
pycircuit solver code.

*Status, 2026-08-22: research complete; **Phases 1 AND 2 (sections 5.1, 5.2) are
implemented** on `cna-jax-vectorization` (per the owner's direction). Phase 1:
`floored_wrap`, the shared `_IdtBase`, pure `i()`, `eval_*_pure` batching, `ic` with
DC pinning via `epar.analysis_kind`, the `IC_KIND='state'` uic route. Phase 2: the
`periodic_states()` hook (Circuit/SubCircuit/Idtmod), the post-acceptance gauge shift
on both CPU paths (`_apply_periodic_shifts`, with the `_q_cache` invalidation) and
branchless inside the JAX `do_accept` (including post-shift `sig_max`/`ref_running`,
so the sigglobal reference no longer inflates with run length); `solve_batched` drops
the shift for elements whose modulus/offset are swept per-lane (Phase-1 fallback,
still correct). **Measured payoff** (`test_long_run_precision_payoff`): phase error
2.2e-16 -- one ulp, flat from tend=500 to 2000 -- with the shift, vs 1.4e-12 growing
to 2.2e-11 without: five orders of magnitude, confirming sec. 3.3 row 3. Footnotes
from landing: the float wrap can land one ulp outside the half-open range at exact
multiples (documented on `floored_wrap`); the ic-without-uic guard exempts
`IC_KIND='state'`; and the bounded state is decoupled from the time grid's
accumulated rounding, so a sample landing exactly ON a wrap can report either limit
-- sawtooth comparisons in tests are congruence-based (two-sided). **Phase 3 (sec.
5.3) is also implemented**, with one measured deviation from the plan:
`Idtmod.accept_step` caches the gauge-invariant `(wrapped phase, rate)` pair -- the
rate read exactly off the element's own ODE right-hand side (`v_ip − v_in`), not
finite-differenced, so the Phase-2 shift ordering is irrelevant -- and `next_event`
predicts the crossing (whole-period advance when sitting on the corner; inf when
stalled, un-seeded, or in idt-degradation mode). Measured: adaptive runs land steps
on crossings to ~1e-15; 22 wraps under BOTH controller paths complete with <=3
rejections and congruence at machine precision; a VSin-driven (varying-rate) VCO
shows only ordinary LTE error. **The kink-gate extension was therefore NOT applied**:
landed corners mean no history ever straddles a wrap, and the collapse it was
reserved for does not occur -- the gate stays TLine-only, per the branch's
gates-are-measurements convention. The Euler pin is retired in effect: a Gear-2/
trap/Euler consistency test asserts congruence under all three. One repair landed
alongside: `JAXTransient` now calls `reset_state` before its static breakpoint
enumeration (Stage 8(d)'s "every analysis" contract; previously a stale element
cache from a prior CPU run could leak into `collect_breakpoints`). JAX dynamic wrap
breakpoints remain architecturally out (static `t_breaks_array`), as sec. 5.3
documents. **Section 7 adds `IdtmodCircular`** (requested as "idtmod_circular"): the
two-state phasor-pair variant with Baumgarte-style orbit correction — researched
(Baumgarte 1972, Ascher et al. 1994, feedback-integrator gain theory, quaternion
normalization feedback, DSP quadrature-oscillator AGC; no circuit-simulation
precedent found) and implemented with a dimensionless per-radian gain `gamma`
(default 1) whose effective per-step strength self-tracks the admissible bracket
through the LTE controller. Measured: Gear-2, 50 cycles, γ=0 → radius decays to
0.48; γ=1 → `max|r−1| = 1.1e-2`; trapezoidal at constant ω conserves the radius
exactly after a single first-step order-drop injection of ½(ωh)², confirming the
Cayley-transform prediction. The scalar `Idtmod` remains the primary
recommendation (ulp-exact phase); the phasor form trades accumulating LTE-level
phase error for a globally smooth state.*

---

## 1. What idtmod is

### 1.1 LRM semantics

Verilog-AMS LRM 2.4.0, §4.5.5 "Circular integrator operator" (identical in the 2023
Accellera edition). Grammar (§A.8.2):

```
idt    ( expr [ , ic [ , assert [ , abstol | nature ] ] ] )
idtmod ( expr [ , ic [ , modulus [ , offset [ , abstol | nature ] ] ] ] )
```

Note the asymmetry: `idt`'s third argument is `assert` (a reset flag), `idtmod`'s third is
`modulus`. `idtmod` has no reset form.

Defining relation (Table 4-19): `idtmod(expr, ic, modulus, offset)` returns **k** where

```
offset <= k < offset + modulus
∫[t0..t] expr dτ + ic  =  n·modulus + k,     n ∈ ℤ
```

Normative prose, verbatim:

> "The output of the idtmod() function shall remain in the range
> `offset <= idtmod < offset+modulus`. The modulus shall be an expression which evaluates
> to a positive value. If the modulus is not specified, then idtmod() shall behave like
> idt() and not limit the output of the integrator. The default for offset shall be
> zero (0)."

and the invariant tying the two operators together:

> "If `y = idt(expr, ic); z = idtmod(expr, ic, modulus, offset);` then
> `y = n * modulus + z` where `offset ≤ z < modulus + offset`."

DC/IC semantics (§4.5.4/4.5.5):

> "The initial condition is optional. If the initial condition is not specified, it
> defaults to zero (0). Regardless, the initial condition shall force the DC solution to
> the system."

i.e. at DC/OP the operator returns `ic`; without an `ic` the operator must sit inside a
feedback loop that forces `expr` to zero, otherwise its DC output is undefined.

§4.5.2 explicitly sanctions the implementation strategy every serious simulator uses:

> "Occasionally, analog operators require new equations and new unknowns be introduced by
> the simulator to convert a module description into a set of first-order differential
> equations."

**Standard hole:** the LRM says the modulus *shall* be positive and that an *absent*
modulus degrades to `idt`, but defines no behaviour for a supplied non-positive modulus.
The sane, recommended choice (matching what a user would expect from "not limiting"):
treat `modulus <= 0` the same as absent, i.e. plain `idt`. This is an implementation
decision and should be documented as such.

**The modulo must be floored, not truncated.** The wrap is
`k = y − m·floor((y − offset)/m) ` (+ offset folded in), which stays in `[offset,
offset+m)` for negative arguments too. C-style truncating `%` is wrong for negative
integrals — and this is not hypothetical: Xyce shipped exactly this bug for the Verilog-A
`%` operator (Xyce Release Notes 7.4/7.5 known-defects table, Gitlab issue 41: Xyce/ADMS
"simply emits the equivalent C++ expression using the same operator").

### 1.2 The canonical use case (Kundert)

Ken Kundert, "Modeling Jitter in PLL-based Frequency Synthesizers"
(designers-guide.org/analysis/PLLjitter.pdf), §7.3, on the VCO model:

> "The phase is computed with idtmod, a function that provides integration followed by a
> modulus operation. **This serves to keep the phase bounded, which prevents a loss of
> numerical precision that would otherwise occur when the phase became large after a long
> period of time.**"

The canonical code (also verbatim in his published `vco.va`,
designers-guide.org/verilog-ams/functional-blocks/vco/vco.va):

```verilog
// phase is the integral of the freq modulo 2π
phase = 2*`M_PI*idtmod(freq, 0.0, 1.0, -0.5);
@(cross(phase + `M_PI/2, +1, ttol) or cross(phase - `M_PI/2, +1, ttol)) begin
    ...
```

Two design lessons hide in these four lines:

1. **The motivation is floating-point precision, not solver convenience.** Over 10⁶ RF
   cycles an unwrapped phase reaches ~10⁷ rad; the double-precision ulp there is
   ~10⁻⁹ rad, so per-cycle phase resolution — the very quantity a jitter simulation
   measures — is quantized away. An implementation that keeps the internal state unbounded
   and wraps only the output does *not* deliver the property the operator exists for.
2. **`idtmod` interacts with event detection.** The modulus is 1 with offset −0.5: the
   range is symmetric so that the `cross()` thresholds at ±π/2 sit away from the wrap
   point. The wrap itself is a sawtooth discontinuity in the returned value; a naive
   implementation produces a spurious "crossing" once per cycle. Any implementation must
   reason about how the wrap interacts with breakpoints/timestep control.

The LRM's own example (§4.5.5) is the same VCO:
`phase = idtmod(fc + gain*V(in), 0, 1, 0); V(OUT) <+ sin(2*`M_PI*phase);`.

---

## 2. Survey: how existing simulators handle it

Summary table (repos examined under `~/source/`, August 2026):

| Tool | `idt` | `idtmod` | Mechanism |
|---|---|---|---|
| ngspice (ADMS era + 47/OSDI) | no (dead vendor code only) | **no** | — |
| ngspice XSPICE `int` model | yes (plain, saturating limits) | no | companion-model inversion, state in `CKTstates` |
| Xyce (ADMS + expression lib) | `SDT()` in B-sources only | **no** | local trapezoid history |
| gnucap (modelgen-verilog) | **yes**, incl. assert/reset | **no** | per-call-site storage element |
| VACASK + OpenVAF-Reloaded | yes | **yes** | extra DAE unknown; wrap on output only |

`idtmod` does not exist in ngspice or Xyce at all — zero grep hits, no parser token, no
error message, not even an entry on the Xyce/ADMS `TODO`. Searches under the alternative
names *circular integrator*, *periodic integrator*, *wrapping/modulo integrator*, *phase
accumulator* also return nothing relevant in any surveyed tree. Whatever we build here is
genuinely absent from the classical open-source SPICE world; the only working prior art is
OpenVAF-Reloaded (§2.4).

### 2.1 Xyce

- **ADMS front-end** (`utils/ADMS/adms.implicit.xml:275,297`): `idt` is parsed and marked
  bias-dependent, then listed among "simulator-handled builtins" with an *empty* handler.
  The back end (`utils/ADMS/xyceBasicTemplates_nosac.xml`) implements `ddt` by *erasing*
  it — the argument is routed into the Q (reactive) vector and Xyce's integrator
  differentiates it — but has **no `idt` branch at all**: `idt(...)` falls through to the
  generic emitter, producing a call to a nonexistent C++ function plus a literal
  `FIXME/*derivatives not implemented*/`, i.e. the generated device does not compile.
  The only `idt` in the tree is the PSP102/103 NQS spline
  (`V(SPLINE1) <+ vnorm_inv * idt(-Tnorm*fk1, Qp1_0);`) behind an `ifdef NQSmodel` that
  Xyce never defines — shipped dead code confirming `idt` is unsupported.
- **Expression library `SDT(x)`** (`src/UtilityPKG/ExpressionSrc/ast.h:6092-6142`) is the
  closest working thing: a per-call-site trapezoid,
  `integral = integral_old + 0.5·(v1+v2)·Δt`, with history rotated by an external
  `processSuccessfulTimeStep()` sweep the author's own comment calls "a kludge"
  (`newExpression.C:2398-2467`). Documented defects: unusable derivatives (known defect
  1006-SON: "this can result in a lack of robustness for circuits that contain SDT-Bsrc
  devices") and broken re-entrancy when called inside a `.func` used twice (patched with a
  per-call-site state stash keyed by a string). **A cautionary tale**: a purely local,
  fixed-order history integrator bolted onto the expression tree, invisible to the
  Jacobian and to LTE control, accumulates exactly the robustness debt their release notes
  describe.

### 2.2 ngspice

- ADMS templates (`src/spicelib/devices/adms/admst/ngspiceMODULE.hxx.xml:79`):
  `#define _DDT(q) q` — same erase-into-charge strategy; no `idt` branch; `idt` falls into
  a generic "emit call by name" that produces nonexistent `idt()`/`d_idt()` C calls. Same
  dead PSP102 NQS vendor code. ADMS was deprecated in ngspice-43 and removed since.
- B-source expressions: `ddt` only (added ngspice-40), with **its Jacobian entry hard-zeroed**
  (`inpptree.c:573-576`). No `idt`/`sdt`.
- **XSPICE `int` code model — the best classical pattern found.** The model
  (`src/xspice/icm/analog/int/cfunc.mod`) calls `cm_analog_integrate()`
  (`src/xspice/cm/cm.c:216-309` → `cm_static_integrate`, `cm.c:536`, "a modified version
  of the function NIintegrate()"), which **inverts the companion model**: SPICE's
  integrator normally computes `dq/dt = ag0·q + ceq` from stored charge history; here,
  given the desired *derivative* (the integrand), it solves for the state:

  ```
  integral = (integrand − ceq) / ag0        ∂integral/∂integrand = 1/ag0
  ```

  The integral lives in a registered slot of the multi-level SPICE state vector
  `CKTstates[0..CKTorder]`, so step rejection/rollback, variable order (Gear 1–6, trap)
  and LTE control (`MIFtrunc`/`MIFterr`, `src/xspice/mif/miftrunc.c:62-180` — the state
  contributes to timestep selection exactly like a capacitor charge) all come for free.
  One warning from the same file: the model's output *limits* are implemented by clamping
  the state and zeroing the partial — correct for a *saturating* integrator, and exactly
  what a *wrapping* integrator must never do (§3).

### 2.3 gnucap (modelgen-verilog)

gnucap is the only surveyed tool with a serious, tested `idt` implementation — including
the LRM's 3-arg assert/reset form — yet `idtmod` appears exactly once in the whole suite:
as a commented-out line of LRM grammar marked "not implemented"
(`src/mg_in_analog.cc:79-86`).

Architecture (two layers):

- Compile time (`mgvams/mg_filt_xdt.cc`): each `idt()` call site becomes an internal
  branch `_b_idt_<n>` with an attached element (`class XDT`, `class IDT` at :281,
  registered `d_idt` at :374; `max_args()==4` — the `idt` list, not `idtmod`'s five). The
  generated code emits per-call-site `tr_advance/tr_eval/tr_regress/tr_accept/tr_review`
  hooks; the reset argument is delivered to the runtime device via a magic parameter index
  and acted on at **acceptance/event time, never inside the residual**.
- Run time (`mgsim/d_va_filter.cc:257`): `class DEV_IDT : public DEV_CPOLY_CAP` reuses
  gnucap's generic `STORAGE`/`integrate()` linear-multistep machinery — the same `_y[]` /
  `_i[]` history arrays a capacitor uses — with LTE in `tr_review()` and analytic
  `1/jω` for AC. `DEV_IDT` is also the state element gnucap's Laplace filters are built
  from: `idt` is their canonical integrator building block.

Relevance: gnucap stores the integrator state *in the element's own LMS history*. Adding
`idtmod` there is precisely the "wrap the state consistently in the integrator history"
problem — you cannot wrap `_y[0]` without also wrapping `_y[1]`, `_y[2]` used by
`integrate()`. They never did it. Their VCO test (`tests/mg_vco.0.gc:22`) is written the
unbounded way — `cos(idt(6.28*(vc*kvco + fc), 0.0))` — i.e. exactly the case `idtmod`
exists to fix.

### 2.4 VACASK + OpenVAF-Reloaded — the only working idtmod

VACASK itself has no `idtmod` anywhere (its built-in behavioral-source language supports
only `idt(x)` and `idt(x, ic)`, `lib/rpnexprva.cpp:114-159`); it delegates all Verilog-A
to the bundled OpenVAF-Reloaded compiler (23.5.0), whose builtin table **does** include
`idtmod` and whose "unsupported" list does not.

**The lowering** (OpenVAF `openvaf/hir_lower/src/expr.rs`): every `idt`/`idtmod` becomes
an implicit unknown `val` plus one extra DAE row, expressed in resistive/reactive residual
parts:

```
transient:  resist = −expr,      react = val         (⇒ d(val)/dt = expr)
DC/IC:      resist = val − ic,   react = ic          (⇒ pins val = ic)
```

switched by a simulator-supplied `EnableIntegration` parameter. At the OSDI ABI level the
extra unknown is indistinguishable from an ordinary internal node — which is why the OSDI
header has no integrator concept at all, and why simulators can get it wrong silently:
VACASK needed a 2026-08 fix (`lib/coretran.cpp:955-975`, commit 381c7f7, replicated in the
PSS engine) because the flag combination it passed left every `idt` frozen.

**The wrap — a documented failure and fix.** The original OpenVAF `idtmod` wrapped the
state *inside the residual* (and had two further bugs: offset read from `args[2]`, the
modulus slot; and a copy-pasted underflow branch). The OpenVAF-Reloaded fix (commit
`d98d73bc`, 2026-07-04) is the money quote for this whole document — from the source
comments:

> "Always integrate the DAE state unbounded; for `idtmod` the modulo wrap is applied to
> the *returned value* (below), not the state. **Wrapping the state inside the residual
> makes the reactive residual jump by `modulus` at each wrap, so the transient
> integrator's d/dt term (based on the previous charge, ~modulus) diverges at the wrap.**"

The returned value is wrapped with the floored modulo
`offset + (val−offset) − m·floor((val−offset)/m)`.

**The accepted caveat:** OpenVAF-Reloaded delivers the LRM's *bounded return value* but
not the *bounded state* — the internal DAE unknown still grows without bound, so
Kundert's precision motivation (§1.2) is not addressed. They chose output-wrapping
because state-wrapping (done naively) breaks the implicit integrator. §3 shows there is a
third option they did not take.

### 2.5 Others

qucs-s: nothing (its `.va` files are digital; no ADMS scripts in tree). Trilinos: nothing
(below the Verilog-A layer, as expected). One instructive near-miss in VACASK
(`docs/dev-builtin-src.md:160`): periodic *waveform sources* wrap their **time argument**
into `[0, period)` before interpolation — the safe "wrap the input of an algebraic map"
pattern, with no integrator history involved. Useful contrast: the danger is never the
modulo itself, it is modulo applied to a *state with history*.

---

## 3. The wrap problem in DAE/ODE terms

### 3.1 Two distinct discontinuities

The idtmod value is a sawtooth: it ramps with slope `expr` and drops by exactly `modulus`
at each wrap. Two problems hide in that, and conflating them is what produced the original
OpenVAF bug:

1. **The integrator's problem.** Implicit multistep formulas (trapezoidal, Gear/BDF) and
   divided-difference LTE estimators assume the state history lies on a smooth polynomial.
   A `modulus`-sized jump in the history makes `dq/dt` estimates and error estimates
   diverge: step rejection to `h_min`, trapezoidal ringing, Newton failure.
2. **The circuit's problem.** Whatever *consumes* the wrapped value sees a genuine
   discontinuity at the wrap instant. This is inherent to the operator's semantics — no
   implementation removes it; it can only be scheduled cleanly (§3.4).

### 3.2 The key insight: the wrap is a gauge shift, not a dynamic event

The jump is always by an exact constant, `n·modulus`. Subtracting a constant from a state:

- does not change its derivative,
- does not change any divided difference or linear-multistep formula, **provided the same
  constant is subtracted from every stored history value of that state** (LMS formulas are
  translation-invariant).

The SUNDIALS documentation states the underlying doctrine precisely:

> "SUNDIALS solvers implement multistep methods and use past solution history stored
> internally to advance the solution in time. Once a discontinuity occurs, the solutions
> at previous steps are not suitable anymore and must be discarded."

Their remedy is rootfinding + `CVodeReInit` — stop at the discontinuity, discard history,
restart at order 1. Correct, but for a VCO the wrap happens **once per RF cycle**: paying
an integrator restart every cycle is prohibitive. Modelica is the same story: `wrapAngle`
is an output-only algebraic block; wrapping a *state* requires `reinit()` inside a
`when`-clause, i.e. a zero-crossing event with an integrator stop/cold-restart.

But for idtmod the general remedy is unnecessary, because the jump is not an arbitrary
discontinuity — it is an exact translation. Shift the state **and its entire stored
history** by the same `n·modulus` between accepted steps, and the integrator provably
never sees anything: every difference it forms is unchanged. No event, no restart, no
order drop, no approximation.

### 3.3 The four options

| # | Approach | Integrator history | Accuracy/convergence | Bounded state? | Who does it |
|---|---|---|---|---|---|
| 1 | Wrap **output** only | untouched, smooth | clean; returned value discontinuous but only *consumed* | **no** — precision still degrades | OpenVAF-Reloaded, Modelica `wrapAngle`, current pycircuit (partly) |
| 2 | Wrap state **inside the residual** | silently inconsistent | `dq/dt` sees a modulus-sized charge jump → diverges; Newton residual discontinuous | yes, but broken | original OpenVAF (reverted), current pycircuit (the `%` in `i()`) |
| 3 | Wrap state **and shift the whole history** at an accepted timepoint | consistent by construction | exact — LMS/LTE invariant under uniform translation | **yes** | **nobody surveyed** |
| 4 | Event + reinit (Modelica `reinit`, CVODE root + `CVodeReInit`) | discarded, rebuilt | correct but order-1 restart every wrap (once per RF cycle) | yes | general ODE tools |

Row 3 is correct, cheap, and unimplemented anywhere we looked — it is pycircuit's
opportunity to be genuinely ahead of the state of the art on this operator.

**Soundness caveat.** The gauge argument is legal *only because nothing in the system may
couple to the unwrapped value*: the idtmod contract guarantees the model only ever sees
the wrapped `k`, and the state's own equation depends only on `dx/dt`. A plain `idt`
state must never be wrapped this way — some consumer may depend on its absolute value.
This is also why the wrap must be exempted for `modulus <= 0` (idt-degradation mode).

### 3.4 The consumer discontinuity

Even with row 3, the *output* is still a sawtooth. Three cases:

- **Consumer periodic with period = modulus** (the typical case: `sin(2π·k/T)`): the wrap
  is invisible; nothing to do.
- **Consumer non-periodic in k** (threshold comparators, square-wave generation): the
  circuit genuinely carries a sawtooth waveform. Best handling is a **predicted breakpoint**:
  after an accepted step, linearly extrapolate the crossing time
  `t_cross ≈ t + (bound − s)/ṡ` and make the solver land a step boundary there — the
  SPICE analogue of a Modelica when-event, but *without* the integrator restart, because
  the history shift keeps the history valid at full order. Prediction only needs to
  bracket the corner; the step-controller's order-drop machinery does the rest.
- **Post-processing**: never interpolate a saved waveform across a wrap point.

---

## 4. Audit of the current pycircuit implementation

`Idt` (`pycircuit/circuit/elements.py:1203`) and `Idtmod` (`elements.py:1245`). Both add a
private node `idt_node`; the state obeys row `idt_node`: `(v_ip − v_in) + ds/dt = 0`, so
`s = −∫(v_ip − v_in)dt`, and the output branch row defines `v_out` from `s`. The entire
modulo behaviour is `Idtmod.i()` (`elements.py:1316-1327`):

```python
_i[branchindex] -= self._G[branchindex, self._idt_index] * x[self._idt_index]
v_mod = (-x[self._idt_index] % self.modulus) + self.offset
_i[branchindex] += v_mod
```

Deficiencies, in decreasing order of severity:

1. **Output-map-only modulo over an unbounded state** — option 1 *and* option 2 of §3.3
   at once, with the drawbacks of both: the state grows forever (precision loss ∝ run
   length; with `relref='sigglobal'` the admissible absolute error on `s` also grows with
   `|s|`, so the wrapped output degrades linearly in run time), *and* the `%` sits inside
   the residual, so `i()` is discontinuous within a Newton solve.
2. **Offset congruence bug.** `(-s % m) + offset` violates the LRM invariant
   `integral + ic ≡ k (mod m)` whenever `offset % m ≠ 0` (the result range is right but
   the *value* is wrong — the wrap boundary must be at `offset + n·m`, not at `n·m`).
   Kundert's canonical `idtmod(freq, 0, 1.0, -0.5)` hits this immediately.
3. **Jacobian vs. residual.** `G[branch, idt] = −1` is the exact slope a.e., but the
   residual jumps by `modulus` at the wrap surface; Newton is handed a linearization valid
   only on the current sawtooth branch. Harmless while the output row is feed-forward;
   a genuine non-smooth-Newton hazard the moment the wrapped output feeds back into the
   integrand (a real PLL, a ΔΣ modulator).
4. **In-place mutation breaks the JAX toolkit.** `_i[branchindex] -= ...` on an immutable
   `jnp` array raises. Untested: the P24 backend-parity probe
   (`tests/test_backend_parity_phaseA.py:791`) only calls `G()`.
5. **No `ic`, no DC story.** Neither element accepts an initial condition; at DC the
   `idt_node` row degenerates to `v_ip − v_in = 0`, which (accidentally) implements the
   LRM's *no-ic* branch — the OP is singular unless the input is forced to zero, so every
   user must know to pass `uic=True` (`doc/time_stepping.rst:154` documents the symptom).
   The missing half is the `ic` branch: DC must pin the output to `ic`. The internal
   `idt_node` is also unreachable from analysis-level `ic={}` (element-private node), so
   starting at a nonzero phase requires hand-building `x0`.
6. **The two step-control paths disagree about the wrap.** The default
   `IntegralController` measures LTE on charge; `q[idt_node] = s` is smooth, so the wrap
   is *invisible* to step control — the doctest's clean sawtooth samples are a
   uniform-grid artefact, not accuracy. The `SolutionLTEController`/coupled path and
   `max_dv_step` *do* see the output jump and will drive the step toward `minstep` at
   every wrap. Related: the wrap sample's left/right-limit convention is
   integrator-dependent — `test_Idtmod_modulo` is pinned to `EulerIntegrator()` with a
   comment that Gear-2 lands the other limit (`tests/test_elements.py:355-360`).
7. **`modulus`/`offset` frozen at `__init__`** (`elements.py:1307-1308`), with no
   `update()` observer hook — a modulus given as an expression over a parent subcircuit
   parameter is silently ignored (every other parameterised element follows the
   `R`/`C`/`L` update pattern).
8. **No `eval_i_pure`/`eval_q_pure`** — no autodiff Jacobian, no vectorised stamping, and
   `solve_batched` refuses any sweep touching `modulus`/`offset`
   (`jaxtransient.py:2261-2280`).
9. Smaller items: the modulo is entirely absent from AC/small-signal (transfer is a clean
   `1/s` — *correct* in principle, since `∂wrap/∂s = 1` a.e., but currently undocumented
   and untested as a deliberate choice); `Idtmod` is not in the finite-difference Jacobian
   harness nor the PCNR inventory; `Idt`/`Idtmod` duplicate ~30 lines of identical stamps;
   negative modulus silently yields the wrong range with Python's `%`.

---

## 5. Recommended implementation

Three phases. Phase 1 is self-contained in the element and reaches parity with
OpenVAF-Reloaded (the current state of the art). Phase 2 adds the §3.3-row-3 gauge shift —
the part nobody else has — and touches the two solver loops. Phase 3 is optional polish.

### 5.1 Phase 1 — a correct element (OpenVAF-Reloaded parity, no solver changes)

Files: `pycircuit/circuit/elements.py`, `pycircuit/circuit/transient.py` (ic plumbing),
`pycircuit/circuit/dcanalysis.py` (flag).

1. **Shared `_IdtBase(Circuit)`** holding the common topology (private `idt_node`, output
   branch), with `G`/`C` built in an `update(subject)` observer hook (the `R`/`C`/`L`
   pattern) — fixes the frozen-parameter defect. `Idt = _IdtBase` with identity wrap;
   `Idtmod` adds `modulus`, `offset`, and both gain a new `ic` parameter (default 0).
   Element-local x layout is fixed by `update_node_map` (terminals, then private nodes,
   then branches): `iplus=0, iminus=1, oplus=2, ominus=3, idt_node=4, branch=5`.
2. **Floored-wrap helper**, module level, toolkit-generic:

   ```python
   def floored_wrap(y, m, o, toolkit):
       # LRM idtmod wrap; m <= 0 handled by caller (degrade to idt)
       return y - m * toolkit.floor((y - o) / m)
   ```

   `numpy.floor` and `jnp.floor` agree; the symbolic toolkit cannot evaluate `floor`, and
   that is acceptable — AC/symbolic analyses use only `G`/`C` (see below). Fixes the
   offset congruence bug and the negative-integral case. `modulus <= 0` (or absent) →
   identity, per the §1.1 decision.
3. **Pure functional `i()`** — no in-place mutation: precompute `_G0` (= `_G` with the
   `(branch, idt)` entry zeroed) and a one-hot branch vector `_e_branch`; then
   `i(x) = dot(_G0, x) + _e_branch * (−floored_wrap(−x[4], m, o))`. Pure array
   arithmetic → JAX-traceable, symbolic-safe for the parts that matter.
4. **`eval_i_pure` / `eval_q_pure` staticmethods** (signature per `R.eval_i_pure`):
   `eval_q_pure → [0,0,0,0, x[4], 0]`;
   `eval_i_pure → [0, 0, x[5], −x[5], x[0]−x[1], −(x[2]−x[3]) − wrap(−x[4])]`.
   Admits `Idtmod` to `evaluation_groups`/`batched_contributions`
   (`toolkit.py:349-457`), which lifts the `solve_batched` refusal — `modulus`/`offset`
   sweeps become possible. Under `jax.jacfwd`, `floor`'s zero gradient makes
   `∂wrap/∂s = 1` automatically, so the autodiffed stamp matches the stored `G`.
5. **DC/IC semantics.** Mechanism: an **epar-scoped analysis flag** (precedent:
   `epar.bypasstol`, `analysis.py:189`), set by `DC.solve` on entry and restored in a
   `finally`. Under the flag the element's `G`/`i`/`u` switch row 4 to the pin
   `s = −ic` (G row: `[0,0,0,0,1,0]`, `u[4] = ic`), so the output is `wrap(ic)` and the
   OP is non-singular — `uic=True` stops being mandatory. This is exactly OpenVAF's
   `EnableIntegration` split (§2.4), including its lesson: the flag *must* be scoped
   correctly, because `Transient._solve_operating_point` forwards the same epar object
   into its inner DC — a leaked flag freezes every integrator for the whole run (VACASK
   shipped that bug).
   For the `uic=True` path, add a third kind `IC_KIND='state'` to `_apply_element_ics`
   (`transient.py:886`): the element declares `state_ic() → [(local_row, value)]`
   (for `Idtmod`: `[(4, −floored_wrap(ic, m, o))]`), mapped through
   `cir.elementnodemap[name]`.
6. **AC**: unchanged and now *documented*: transfer stays `1/s`; the wrap is invisible to
   small-signal because its derivative is 1 a.e. (state the assumption: the linearization
   point is not exactly on the wrap surface).

### 5.2 Phase 2 — bounded state via the post-acceptance gauge shift

The differentiator (§3.3 row 3). Files: `pycircuit/circuit/circuit.py`,
`pycircuit/circuit/transient.py`; JAX part later in `pycircuit/circuit/jaxtransient.py`.

1. **Registration API.** The element *cannot* do the shift itself: `accept_step` hands
   elements a fancy-indexed **copy** of `x` (`circuit.py:1476-1483`), and the history
   rings live in the analysis. So the element declares and the analysis acts:

   ```python
   class Circuit:
       def periodic_states(self):
           """[(local_x_index, modulus, offset)] for states defined only up to
           n*modulus. Contract: q[row] == x[row] (unit C diagonal), so shifting
           x and the q history by the same n*modulus is an exact translation."""
           return []

   class SubCircuit:
       def periodic_states(self):
           out = []
           for name, el in self.elements.items():
               for r, m, o in el.periodic_states():
                   out.append((self.elementnodemap[name][r], m, o))
           return out
   ```

   `Idtmod` returns one entry for row 4, with the (modulus, offset) expressed in the
   internal `s = −phase` convention so the state is kept within one modulus of the range
   the wrap in `i()` expects. Polled by the analysis once per solve, after
   `update_iparv` (parameters are static per run); assert the unit-C-diagonal contract at
   poll time so a non-conforming element fails loudly. `modulus <= 0` entries are dropped
   (idt-degradation mode — never wrap a plain integral, §3.3 caveat).

2. **CPU accept-time shift** (`transient.py`; two call sites — after the history push in
   `_solve`, and the analogous point in `_solve_coupled`):

   ```python
   def _apply_periodic_shifts(self, X):
       for row, m, o in self._periodic_rows:
           n = floor((X[-1][row] - o) / m)
           if n != 0:
               d = n * m
               for k in (-1, -2, -3):        # bounded by len(X)
                   X[k][row] -= d
               self._qlast[:, row] -= d      # whole ring
               self._q_cache = None
   ```

   The shift set was verified complete against every consumer of history in the tree:
   - Integrators (`integrator.py`): `compute_derivatives` reads `q_curr, q_last[0..1],
     iq_last[0]`; `compute_lte` reads `q_curr, q_last[0..2], iq_last[0..1]` — all
     invariant under a uniform shift of the q column; `iq` (a derivative) is invariant
     outright, so `_iqlast` and `_dt_last*` are untouched.
   - `IntegralController`/`PIController`: consume fresh `q(x)` + `_qlast` + `_iqlast` +
     `J` — consistent. Their running references (`_ref_running` etc.) are running maxima;
     with the state now bounded they simply stop inflating — no correction needed.
   - `SolutionLTEController` and `fang_timestep`: read exactly the last 3 solution
     vectors — precisely the entries shifted.
   - PCNR: per-timepoint, no persistent state — untouched.
   - The next step's Newton starts from `X[-1]`, i.e. in-range — which is the point.

   **One hazard found in review:** `_q_at` (`transient.py:803-822`) memoises `(x, q)` with
   an identity-first guard; shifting the accepted `x` in place would keep identity and
   serve a **stale q**, silently corrupting the LTE estimate. The `_q_cache = None`
   invalidation above is mandatory and deserves its own regression test.

   Order at the accept site: push history rings → shift → *then* call `accept_step`, so a
   Phase-3 element caches the post-shift state (`TLine` reads only port rows; unaffected).
   The shift is `n·modulus` with integer `n` — an exact translation, applied only to
   accepted state, never inside a Newton solve.

3. **JAX traced path** (`jaxtransient.py` — future work on its own branch; explicitly
   *not* part of `cna-jax-vectorization`). The shift is attractively branchless: pre-trace,
   poll `periodic_states()` into static arrays `(p_rows, p_mods, p_offs)`; inside
   `do_accept`, before the results buffer is written:

   ```python
   n = jnp.floor((x_curr[p_rows] - p_offs) / p_mods)
   d = n * p_mods
   x_curr     = x_curr.at[p_rows].add(-d)
   x_hist_new = x_hist_new.at[:, p_rows].add(-d)
   q_hist_new = q_hist_new.at[:, p_rows].add(-d)
   ```

   `iq_history`, `h_history` and the running references are invariant. Empty `p_rows`
   reduces to no-ops — untraced circuits compile identically. The candidate step's LTE was
   computed pre-shift in a single consistent gauge, which is correct. Fully
   `vmap`-compatible, no event detection in the trace — consistent with the branch's
   existing policy of keeping rescue/decision logic out of the traced step.

4. **Recorded-waveform caveat** (document as an invariant): the shift rewrites the private
   `idt_node` row of the last ≤3 *recorded* points, so that row's stored waveform is
   gauge-dependent and not user-meaningful. The observable output — the branch-defined
   wrapped voltage at `oplus` — is untouched and correct at every timepoint, including
   within-step values before a rewrap (the wrap in `i()` covers any excursion up to the
   next acceptance).

### 5.3 Phase 3 — wrap breakpoints and controller interaction (optional)

- **CPU wrap breakpoints.** `next_event(t)` is polled live each loop iteration but takes
  only `t`; so `Idtmod.accept_step(t, subx, epar)` caches `(s, ṡ)` and `next_event(t)`
  returns the linearly predicted crossing time (strictly `> t` per the contract, `inf` if
  the rate is ~0). Prediction only needs to *bracket* the corner; the existing
  break-step → first-step order-drop machinery handles the rest. This resolves the
  consumer discontinuity (§3.4) for non-periodic consumers and makes the wrap sample's
  left/right-limit convention integrator-independent (the Euler-pinned test can be
  un-pinned).
- **Controller interaction.** Phase 2 does not remove the *output* sawtooth; the
  `SolutionLTEController`/coupled path will treat each wrap like a pulse edge (shrink
  toward minstep) while the charge-based default controller sails through. The remedy is
  the existing kink discipline (`transient.py:3010-3040`, currently gated on
  `_tline_tds`): extend the gate to "circuit has periodic states AND this accepted step
  wrapped" (detectable as `n != 0` in the shift helper).
- **JAX breakpoints are impossible today** — `t_breaks_array` is precomputed statically
  before tracing, and `max_timestep()` is polled once per run; a state-dependent event
  cannot enter either. Document honestly as CPU-only; a branchless in-trace `dt` cap from
  the predicted crossing is listed as future work, not promised.

### 5.4 Test plan

1. **FD Jacobian**: add `Idtmod` to the finite-difference harness pattern
   (`tests/test_element_jacobians.py`): `G` vs central differences of `i` away from the
   wrap surface; plus an explicit ±ε-across-the-surface test (branch-row jump of exactly
   `m`, slope unchanged).
2. **`floored_wrap` unit test** against the LRM congruence, including negative integrals
   and `offset` not a multiple of `modulus` (catches the live offset bug, §4 item 2).
3. **DC with ic**: OP returns `v(out) == wrap(ic)`; transient without `uic` starts there;
   `uic=True` + element `ic` seeds `x0` via the new `IC_KIND='state'` route.
4. **Long-run precision** — the Phase-2 payoff, measured: constant input into
   `modulus=1` for ~10⁴ wraps; assert the phase error vs analytic stays flat with the
   gauge shift, and demonstrate the `eps·|state|` drift of the unbounded (Phase-1-only)
   variant. This is the claim the recommendation stands on.
5. **Integrator consistency at wraps**: the `test_Idtmod_modulo` circuit under
   Euler/Trapezoidal/Gear2 must agree at wrap samples (un-pin Euler once Phase 3 lands).
6. **JAX parity**: `Idtmod.i()` under the JAX toolkit (regression for the in-place
   mutation crash), stamp parity in the P24 harness, and a JAX-vs-CPU transient waveform
   comparison spanning ≥2 wraps.
7. **`solve_batched`** sweep of `modulus`/`offset` lanes (exercises the `eval_*_pure`
   batching path end to end).
8. **Feedback/PLL-ish loop**: wrapped output driving a comparator feeding back into the
   integrand — convergence and no minstep collapse at wraps on both the default and
   `coupled_lte=True` paths, and controller-path agreement per the `_compare_methods`
   idiom in `tests/test_analysis_transient_stress.py`.

### 5.5 Risks and open questions

- **epar flag scoping** is the weakest joint of Phase 1: Transient shares its epar with
  its inner DC solve; a missed `finally` restore silently pins every integrator for the
  whole run (VACASK's actual shipped bug). Alternative if this worries review: thread
  `analysis_kind` as an explicit kwarg through `G`/`i` — far more invasive.
- **Recorded gauge history**: rewriting the last ≤3 result rows of a private node is safe
  today, but becomes a landmine if `idt_node` rows are ever exposed in results
  post-processing. Name the invariant: *only the branch output is observable*.
- **`_q_cache` identity guard**: the one place an in-place shift corrupts silently;
  explicit invalidation + regression test.
- **Coupled path at wraps without Phase 3**: expect step collapse at each wrap under
  `coupled_lte=True` (same class as the known TLine kink findings) — Phase 3's gate
  extension is the fix; until then, document the limitation.
- **Symbolic toolkit** cannot evaluate `floor`: symbolic transient/DC of `Idtmod` is out
  of scope (AC unaffected).
- **`hdl.py` (`Behavioural` DSL)** has only `ddt`; adding `idt`/`idtmod` sympy functions
  that lower onto `_IdtBase` stamps is a natural later extension — noted, not designed.

---

## 7. IdtmodCircular — the two-state circular integrator with Baumgarte orbit correction

Section 3.3 listed the phasor-pair state as an alternative "to document, not recommend."
This section is its full treatment: research into the drift problem and its cures, and
the implemented element (`IdtmodCircular` in `pycircuit/circuit/elements.py`; named per
the repo's CamelCase element convention, requested as "idtmod_circular").

### 7.1 The idea

Instead of a scalar phase state that must be folded, integrate the phase's unit phasor
`(c, s) = (cos θ, sin θ)`, `θ = 2π·(k − offset-fold)/modulus`:

```
dc/dt = −ω·s − γ·|ω|·c·(r² − 1)
ds/dt = +ω·c − γ·|ω|·s·(r² − 1),     ω = 2π·(v_ip − v_in)/modulus,  r² = c² + s²
k     = floored_wrap(modulus·atan2(s, c)/(2π), modulus, offset)
```

The **state is smooth for all time** — no wrap, no Phase-2 gauge shift, no discontinuity
anywhere in the trajectory. The `atan2` branch cut is exactly cancelled by the output
wrap (both are mod-`modulus` in `k`), leaving only the one sawtooth corner idtmod
semantics require — handled by the same `_WrapEvents` breakpoint prediction as `Idtmod`
(the output rate is `v_ip − v_in` for both).

### 7.2 Why orbit correction is needed: the drift problem

Discretized rotation does not conserve the circle invariant `r² = 1`. The
per-method behaviour, from the literature and confirmed by measurement:

| method | `r²` behaviour (constant ω) | mechanism |
|---|---|---|
| forward Euler | grows ×(1+(ωh)²) per step | amplification factor > 1 on the imaginary axis |
| backward Euler | decays ×1/(1+(ωh)²) | over-stability / numerical damping |
| Gear2/BDF2 | decays (slower than BE) | numerical damping — Kundert's classic LC-tank observation: "if either Gear2 or backward Euler are used, the amplitude of the oscillation asymptotically decays" |
| trapezoidal | **exactly conserved** | for linear constant-ω the trap update equals implicit midpoint = the Cayley transform of a skew matrix, which is orthogonal |
| trapezoidal, varying ω(t) | O(h²·ω̇) drift per step | endpoint matrices differ; no single Cayley transform |

Two important corrections to naive expectations, both verified:

- **"Trapezoidal preserves quadratic invariants" is false in general** — trap is
  Lobatto IIIA, symmetric but NOT symplectic (Hairer/Lubich/Wanner, *Geometric
  Numerical Integration*, ch. IV/VI); only for the linear constant-coefficient case
  does it coincide with the (invariant-conserving) implicit midpoint rule.
- **Even trapezoidal runs drift in practice here**, because the Phase-3 wrap
  breakpoints order-drop to Euler at every corner, and the run's first step is an
  Euler step too. Measured: pure trap at constant ω with no wrap crossed shows exactly
  ONE injection of ½(ωh)² = 7.9e-5 from the first-step order drop, then stays flat;
  a 50-cycle Gear-2 run with γ=0 decays the radius to **0.48**. VACASK's own
  oscillator demo (`demo/trannoise/osc.sim:52`) documents the same bind in
  production: Gear chosen to avoid trap ringing, then "a small timestep to reduce
  Gear's damping" — if users run Gear for ringing reasons, a phasor state WILL
  spiral inward unless corrected.

### 7.3 The Baumgarte correction

Baumgarte's original method (J. Baumgarte, "Stabilization of constraints and integrals
of motion in dynamical systems", *CMAME* 1:1–16, 1972) replaces a position constraint
`g = 0` by the critically-damped `g̈ + 2α·ġ + β²·g = 0`. For a *first-order invariant*
like ours the analogue is `ḣ + γ_eff·h = 0`, obtained by appending `−γ_eff·∇V` to the
RHS with `V = ¼(r²−1)²`, `∇V = (c,s)·(r²−1)` — exactly the terms above. The radial
error then obeys

```
d(r²)/dt = −2·γ·|ω|·r²·(r²−1)
```

so the unit circle is exponentially attracting with rate `2γ|ω|`. Two properties worth
stating precisely:

- **The correction is purely radial**: its contribution to the phase flow is
  `(c·ṡ − s·ċ)_corr = c·(corr·s) − s·(corr·c) = 0`. It cannot bias the phase — only
  rescue the radius. (The identical structure is standard for quaternion attitude
  kinematics: `q̇ = ½Ω·q + η·(1−‖q‖²)·q`, e.g. arXiv:1107.1119; named as Baumgarte
  stabilization over SO(3) by Gros et al., IEEE CDC 2015.)
- **`(c,s) = (0,0)` is a degenerate equilibrium** the correction cannot leave
  (`d(r²)/dt = 0` at r = 0) — which is why the element always pins/seeds on the
  circle (sec. 7.5).

### 7.4 Gain selection — the classic critique, and why γ is per-radian here

The single most-cited weakness of Baumgarte stabilization is that parameter choice is
ad hoc: "a popular method, but the choice of parameters to make it robust is unclear
in practice" (Ascher, Chin, Petzold, Reich, "Stabilization of DAEs and invariant
manifolds", *Numer. Math.* 67:131–149, 1994), and the admissible region depends on the
integrator AND the step size (Flores et al., *JCND* 6(1), 2011) — γ is a property of
(model, integrator, h), not of the model alone. The cleanest quantitative theory found
("Feedback Integrators", arXiv:2512.01528) gives, for Euler discretization with
per-step feedback strength `β = h·γ_eff`:

```
admissible:  β ∈ (0, 2/L),   optimal:  β* = 1/L,   L = sup‖∇²V‖ ≈ 2 near r = 1
```

i.e. `γ_eff·h ≲ 1` with `γ_eff ≈ 1/(2h)` near-optimal (consistent with Ascher's
`γ ≈ 1/h` recommendation up to the factor-2 convention).

`IdtmodCircular` cannot know `h` — elements never see the step size — so its gain is
made **dimensionless per radian of phase travel**: the RHS term is `γ·|ω|·x·(r²−1)`,
giving `β = γ·(|ω|·h) = γ·(phase per step)`. The LTE controller keeps phase-per-step
roughly constant (a fraction of a radian), so `β` automatically lands inside the
admissible bracket at any oscillation frequency, and the default `γ = 1` is safe
without tuning. Measured (50 Gear-2 cycles at ωh ≈ 0.13): γ=0 → radius 0.48;
γ=1 → `max|r−1| = 1.1e-2`, final `|r−1| < 3e-4`.

### 7.5 Alternatives considered and rejected

- **Post-step orthogonal projection** (`(c,s) ← (c,s)/r`, or the one-Newton-step
  `×(3−r²)/2` renormalization from the quaternion literature): preferred in the ODE
  literature — parameter-free, exact per step (Ascher et al.'s bottom line is "use
  orthogonal projections whenever possible"). Rejected HERE because projection is not
  a linear-multistep invariant: unlike Phase 2's translation gauge (which LMS
  formulas provably cannot see), radial rescaling of the accepted state is
  inconsistent with the stored q/x history, and rewriting history radially is
  inexact. The Baumgarte term is pure RHS + Jacobian content — invisible to the
  integrator machinery — which matches both Ascher's advice to "discretize the
  stabilizing term simply, differently from the ODE" and VACASK's in-repo idiom of
  filtering the predictor while leaving the solution untouched.
- **Periodic renormalization** (INS practice: normalize every N steps at low rate):
  same history-consistency objection, in milder form; noted, not needed.
- **CORDIC-style precomputed gain correction**: inapplicable — CORDIC's magnitude
  error is known a priori and constant; ours is solution- and history-dependent.
- **Low-pass-filtered Baumgarte correction**: practiced where the correction is an
  EXPLICIT, sampled feedback loop around a noisy violation signal — real-time
  physics/robotics engines (gain-limited, slop-deadzoned, compliance-relaxed
  corrections; ODE's CFM is a leaky first-order-filtered constraint), co-simulation
  force coupling (filtering the exchanged force avoids energy injection at the
  sampling rate), and decimated INS quaternion renormalization. Filtering trades
  correction bandwidth for noise rejection at the cost of phase lag. Not adopted
  here: our violation `r²−1` is deterministic (no measurement noise, no contact
  chatter), and the correction is co-solved IMPLICITLY inside Newton under A-stable
  integrators — no sampling loop to alias against, which is the structural
  difference from the settings that need the filter. One scenario could reopen
  this: trapezoidal ringing appearing in `r²−1` and being fed back with
  alternating sign. The idiomatic fix would then be VACASK's pattern (filter the
  signal used by the correction, leave the solution untouched) — at the cost of
  one extra filter state per element. Reconsider on evidence: ringing-coupled
  radius oscillation, or step-size collapse traced to the γ term.
  **The trigger is now a live sentinel** (`test_circular_trap_ringing_sentinel`):
  the mechanism is exact — trap on the radial mode `λ = 2γ|ω|` has amplification
  `(1−λh/2)/(1+λh/2)`, negative for `λh > 2`, so `γ|ω|h > 1` rings the radius at
  the step Nyquist. Measured: alternation fraction switches 0.00 → 1.00 cleanly
  at the predicted boundary (γ=1,5 vs γ=20,40 at ωh=0.0628); the default γ=1
  sits ~15× below it; and even in the fully-ringing regime the OUTPUT error is
  unchanged (1.6e-3 vs 1.7e-3) — the ringing is confined to the radial mode by
  the correction's phase-orthogonality, so the failure the filter would fix
  cannot reach the waveform. The sentinel fails if a future change pushes the
  default into the ringing regime, and this bullet is its remedy note.

**Negative finding:** no precedent was found for a phasor-pair idtmod state in any
circuit simulator — Verilog-A VCO models integrate a scalar phase and wrap it; the
quadrature VCOs in the literature are circuit topologies, not numerical techniques.
`IdtmodCircular` is a transplant from multibody/aerospace/DSP into circuit
simulation. The DSP world does carry the exact same correction: coupled-form
quadrature oscillators drift off the unit circle and US patent 6,642,796 feeds
`(sin²θ + cos²θ)` back into the loop gain — literally the `−γ(r²−1)` term, discrete.

### 7.6 Implementation, semantics, and the honest trade-off

Element `IdtmodCircular` (`elements.py`): 4 terminals as `Idtmod`, two private nodes
(`cos_node`, `sin_node`) + output branch; single-source-of-truth `eval_i_pure`
(autodiff Jacobian under JAX, central differences under numeric — the Semiconductor
pattern); constant C (unit diagonal on the phasor rows); `eval_q_pure` for the
batched path. Deliberate semantic differences from `Idtmod`:

- `modulus` must be positive — a phasor cannot represent an unbounded plain integral,
  so there is NO idt-degradation mode; `modulus <= 0` raises at construction.
- `ic` defaults to 0 and ALWAYS defines the DC solution (pin folded into
  `eval_i_pure`'s dc branch, so FD/autodiff of one function yields the DC stamp too)
  and the `uic` seed (`IC_KIND='state'`, both phasor rows) — the state must start on
  the circle, since the centre is a degenerate equilibrium.
- Under `JAXTransient` use `uic=True` or `x0`: its internal DC solve does not carry
  the `epar.analysis_kind` pin and would settle on the degenerate centre.
- No `periodic_states()` — there is nothing to shift; the smooth state is the point.
- Not supported on the symbolic toolkit (`atan2`/`floor`).

**The honest trade-off, measured:** the scalar `Idtmod` remains the primary
recommendation. Head-to-head (1 V input, modulus 1, fixed dt = 0.01, ωh ≈ 0.063;
congruence error at end of run):

| integrator | element | 2 cycles | 20 cycles | 200 cycles |
|---|---|---|---|---|
| Euler | Idtmod | 2.2e-16 | 3.1e-13 | 3.7e-11 |
| Euler | IdtmodCircular | 1.2e-3 | 1.3e-2 | 1.3e-1 |
| Trapezoidal | Idtmod | 6.7e-16 | 3.2e-13 | 3.7e-11 |
| Trapezoidal | IdtmodCircular | 6.7e-4 | 6.7e-3 | 6.7e-2 |
| Gear2 | Idtmod | 1.4e-15 | 3.2e-13 | 3.7e-11 |
| Gear2 | IdtmodCircular | 2.6e-3 | 2.6e-2 | 2.6e-1 |

Two different error regimes. `Idtmod` has NO method error for constant input (a
linear phase is integrated exactly by every LMS formula — note the integrator-
independence); the residue is float noise of accumulation, ~N·ulp. `IdtmodCircular`
pays a **carrier-rate tax**: per-cycle phase lag of the method at the oscillation
frequency — Euler ≈ (ωh)²/6, trap ≈ (ωh)²/12, Gear2 ≈ (ωh)²/3 per cycle, matching
the table exactly — accumulating linearly with cycle count, present even for
constant input, reducible quadratically in h but never absent. At 200 cycles the
gap is ~9 orders of magnitude. For accumulated-phase measurements (jitter, long
PLL runs) the scalar form is the only right tool; `IdtmodCircular` is for when a
globally smooth STATE matters more than ultimate phase accuracy — no discontinuity
anywhere in the trajectory, and a quadrature (cos, sin) pair available internally.

### 7.7 IdtmodQuadrature — the two states as outputs, no atan2

The terminal-level answer to "when would the phasor form actually win": as
implemented, `IdtmodCircular`'s output is the same sawtooth as `Idtmod`'s, so at
the terminals it never wins. `IdtmodQuadrature` removes the sawtooth instead of
re-deriving it: six terminals, the same Baumgarte-corrected phasor dynamics, and
the two output branches carry the states directly —
`v = amplitude·cos(2π·phase/modulus)` and `amplitude·sin(...)` — with **no wrapped
output and no `atan2` anywhere**. Consequences:

- Every node voltage and every state is smooth for all time: no wrap corner, no
  `_WrapEvents` breakpoints, no events — measured, zero breakpoints fire where the
  scalar and circular forms land one per cycle.
- **Periodic-steady-state readiness** (the reason it exists as a separate
  element): over one output period the state vector returns exactly to itself,
  which shooting/PSS methods and their monodromy matrices are built on. The
  scalar `Idtmod` phase structurally cannot do this — unwrapped it advances by
  one modulus per period; wrapped it jumps. When a PSS/shooting analysis is
  implemented, this is the phase-accumulator representation it can consume.
- Natural macromodel for sinusoidal VCOs, I/Q sources, and mixer LOs (the
  quadrature pair is the deliverable, not the phase).
- `offset` is meaningless without a wrapped output and is dropped; `ic` is the
  initial phase (starting angle `2π·ic/modulus`), always pinning DC and seeding
  `uic`; `modulus`, `gamma`, and all `IdtmodCircular` semantics otherwise carry
  over — including the accuracy regime (per-cycle carrier-rate lag, sec. 7.6):
  when accumulated phase is the measurement, `Idtmod` remains the tool.

### 7.8 References for this section

- Baumgarte, *CMAME* 1:1–16 (1972) — https://www.sciencedirect.com/science/article/abs/pii/0045782572900187
- Ascher, Chin, Petzold, Reich, *Numer. Math.* 67:131–149 (1994) —
  https://link.springer.com/article/10.1007/s002110050020 (TR-92-17:
  https://www.cs.ubc.ca/sites/default/files/tr/1992/TR-92-17.pdf ; note: the
  "γ ≈ 1/h near-optimal" statement was confirmed at abstract level only)
- Flores et al., "A Parametric Study on the Baumgarte Stabilization Method",
  *JCND* 6(1) 011019 (2011); and the Eng. w/ Computers 2025 integrator-dependence study
- "Feedback Integrators: gain selection under Euler discretization" — arXiv:2512.01528
- Gros et al., "Baumgarte Stabilisation over the SO(3) Rotation Group for Control",
  IEEE CDC 2015 — https://cdn.syscop.de/publications/Gros2015.pdf (equations not
  independently verified — PDF undecodable during research)
- Quaternion normalization feedback: arXiv:1107.1119; Deutschmann et al. 1993
- Hairer, Lubich, Wanner, *Geometric Numerical Integration*, 2nd ed., ch. IV/VI
  (https://www.unige.ch/~hairer/poly_geoint/intro.pdf)
- Kundert, *The Designer's Guide to SPICE and Spectre*, App. A (Gear/BE damping vs
  trapezoidal ringing; passage surfaced by search, page unverified)
- Lyons, "A New Contender in the Quadrature Oscillator Race" —
  https://www.dsprelated.com/showarticle/1467.php ; US patent 6,642,796
- VACASK: `lib/coretran.cpp:1275` (trap ringing filter, "solution left untouched"),
  `demo/trannoise/osc.sim:52` (Gear-vs-damping oscillator note)
- Repo greps: **zero** Baumgarte hits in trilinos/xyce/vacask/gnucap/ngspice;
  Trilinos "quaternion" hits are tensor type names only.

---

## 8. Roadmap — follow-up work (2026-08-22)

The idtmod plan proper (sections 5 and 7) is complete and pushed. What remains,
in priority order agreed with the owner:

**Near-term fixes (small, semantic-gap closers):**

1. **`solve_batched` instance-name quirk** — **DONE 2026-08-22.** The tree is
   class-keyed *by design* (per-class vmap groups); an instance key is now
   remapped when unambiguous and refused with the class spelling otherwise.
   The investigation also caught a latent Phase-2 bug: the periodic-shift
   exclusion matched by instance name against class-keyed trees — now matches
   by class (sharp test: swept modulus 0.7 against trace-time 1.0).
2. **JAX DC pin** — **DONE 2026-08-22.** `solve`'s guard now exempts
   state-kind ics (its DC already went through `dcanalysis.DC`);
   `solve_batched`'s traced DC runs inside the shared `analysis_kind`
   scope (new helper in analysis.py). Fixing it surfaced that the batched
   path assembles through `eval_i_pure`, so the Idt/Idtmod pin moved to a
   single fold-into-i convention shared by both paths (mixed ic/no-ic
   groups fall back loudly to the singular no-ic behaviour, never a
   silent pin).
3. **Per-lane modulus in the traced gauge shift** — a swept `modulus` currently
   drops the Phase-2 shift for that element (correct fallback, unbounded
   state). Make `p_mods`/`p_offs` per-lane arrays fed from the override tree
   so swept lanes keep the bounded-state property.
4. **Doc-build health** — pre-existing dead `exec-rst` blocks (ginac import on
   two pages, a `ValueError` on distortion_limits) warn "figures NOT live" on
   every build. Fix or consciously retire; a live-doc system that ships dead
   blocks erodes exactly the trust it exists to build.

**Larger arcs (need scoping):**

5. **PSS/shooting for phase circuits** — the payoff `IdtmodQuadrature` was
   built for: make the shooting analysis consume phase-accumulator circuits
   (oscillator steady state; phase noise later). Its state vector is exactly
   periodic over a cycle, which the scalar phase structurally is not.
6. **`hdl.py` `idt`/`idtmod`** — the Behavioural sympy DSL knows only `ddt`;
   add `idt`/`idtmod` functions lowering onto `_IdtBase`/`Idtmod` stamps, a
   step toward a Verilog-A-shaped frontend (sec. 5.5).
7. **JAX in-trace wrap-crossing dt cap** — dynamic breakpoints are
   architecturally out of the traced loop (static `t_breaks_array`); a
   branchless per-step `dt` clamp from the predicted crossing is the listed
   future work from Phase 3 (sec. 5.3).

Also standing: the trap-ringing sentinel's reopening trigger (sec. 7.5) and
the transient-rescue-in-trace trigger parked in the branch's own docs.

## 9. References (sections 1–5)

**Standards & canonical usage**
- Verilog-AMS LRM 2.4.0, §4.5.2/4.5.4/4.5.5 —
  https://designers-guide.org/verilog-ams/VlogAMS-2.4.0.pdf
  (Accellera mirror: https://www.accellera.org/images/downloads/standards/v-ams/VAMS-LRM-2-4.pdf)
- Kundert, *Modeling Jitter in PLL-based Frequency Synthesizers* —
  https://designers-guide.org/analysis/PLLjitter.pdf
- Kundert, `vco.va` functional block —
  https://designers-guide.org/verilog-ams/functional-blocks/vco/vco.va

**OpenVAF / VACASK**
- OpenVAF issue #55 "Verilog-A idt operator supported?" —
  https://github.com/pascalkuthe/OpenVAF/issues/55
- OpenVAF PR #57 "fix idt operator" — https://github.com/pascalkuthe/OpenVAF/pull/57
- OpenVAF-Reloaded — https://github.com/OpenVAF/OpenVAF-Reloaded
  (commit `d98d73bc`, 2026-07-04: "idtmod wrap/offset" fix; `openvaf/hir_lower/src/expr.rs`)
- VACASK: `lib/coretran.cpp:955-975` (EnableIntegration fix, commit 381c7f7),
  `lib/rpnexprva.cpp:114-159`, `docs/dev-builtin-behavioral.md`

**Classical SPICE prior art**
- ngspice XSPICE: `src/xspice/cm/cm.c:216-309,536-632` (`cm_analog_integrate` /
  `cm_static_integrate`), `src/xspice/mif/miftrunc.c:62-180`,
  `src/xspice/icm/analog/int/cfunc.mod`
- Xyce: `src/UtilityPKG/ExpressionSrc/ast.h:6092-6210` (`sdtOp`),
  `newExpression.C:2398-2467`, `utils/ADMS/adms.implicit.xml:275,297`,
  `utils/ADMS/xyceBasicTemplates_nosac.xml:5339-5370`, Release Notes 7.4/7.5
  known-defects (1006-SON, Gitlab issue 41)
- gnucap: `gnucap-modelgen-verilog/mgvams/mg_filt_xdt.cc` (`XDT`/`IDT`),
  `mgsim/d_va_filter.cc:257` (`DEV_IDT`), `src/mg_in_analog.cc:79-86`,
  `tests/mg_idt.*.gc`, `tests/mg_vco.0.gc`

**DAE/ODE tooling**
- SUNDIALS discontinuity handling — https://computing.llnl.gov/projects/sundials/usage-notes ,
  https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html
- Modelica: `Modelica.Math.wrapAngle` — https://doc.modelica.org/om/Modelica.Math.wrapAngle.html ;
  `reinit`/events — ModelicaReference Operators;
  *Event Handling in the OpenModelica Compiler and Runtime System* —
  https://liu.diva-portal.org/smash/get/diva2:67/FULLTEXT01.pdf

**pycircuit (audited at fca1640, branch cna-jax-vectorization)**
- `pycircuit/circuit/elements.py:1203-1327` (`Idt`, `Idtmod`)
- `pycircuit/circuit/transient.py` (residual assembly ~1159, `get_diff`:741,
  `_q_at`:803-822, `_apply_element_ics`:886, kink discipline ~3010)
- `pycircuit/circuit/integrator.py`, `pycircuit/circuit/stepcontroller.py`
- `pycircuit/circuit/jaxtransient.py` (`do_accept` ~1589, `solve_batched` ~2248)
- `pycircuit/circuit/circuit.py` (`accept_step`/`next_event` hooks:453-487,606;
  `update_node_map`:1281)
