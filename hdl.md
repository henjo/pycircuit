# hdl.py — the Behavioural element DSL

Research, design, and reference for `pycircuit/circuit/hdl.py`: a
Verilog-A-shaped language for writing circuit elements as *equations*,
compiled symbolically into MNA stamps.

*Status, 2026-08-22: implemented and in use. `pycircuit/circuit/
elements_hdl.py` re-expresses ten hand-written elements in the DSL and
`pycircuit/circuit/tests/test_elements_hdl.py` (28 tests) proves each
equivalent to its `elements.py` counterpart; `benchmarks/hdl_overhead.py`
measures the cost. Every number below is measured, not estimated.*

---

## 1. What it is

An element's terminal behaviour is written as *contribution statements* —
the DSL form of Verilog-A's `<+` — inside an `analog()` method whose
argument names declare the terminals:

```python
class MyResistor(Behavioural):
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm',
                            default=1e3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, b.V / r)
```

That is the whole element. The metaclass compiles it, at class-definition
time, into everything a pycircuit element must supply: `i`/`q`/`u`
vectors, exact `G`/`C` Jacobians, `dudt`, the noise matrix `CY`, DC
variants for pinned states, `eval_i_pure`/`eval_q_pure` for the JAX
batched path, `state_ic` seeding, and `periodic_states` declarations.

**Why a DSL at all**, when `elements.py` exists: an element written by
hand states its physics *twice* — once as `i(x)` and again as the `G(x)`
you differentiated yourself — and nothing checks that the two agree. That
is not hypothetical; it is exactly the defect `test_element_jacobians.py`
was written for after `VCVS_limited` shipped a `G` that dropped its gain
entirely. In the DSL the equation is stated once and the Jacobian is
*derived*, so the two cannot disagree.

---

## 2. Theory: from contributions to stamps

### 2.1 The residual

pycircuit's transient residual (`transient.py`) is

```
f(x, t) = i(x) + d/dt q(x) + u(t) = 0,     J = G(x) + companion(C(x))
```

so an element must decompose its behaviour into a static current `i`, a
charge `q` whose time derivative the *simulator's integrator* takes, and
an independent excitation `u(t)`. Everything below is about performing
that decomposition automatically and provably.

### 2.2 The unknown vector

The compiler builds the element's local unknowns in exactly the order
`Circuit.update_node_map` assigns them:

```
x = [ terminal nodes | internal nodes | state nodes | branch currents ]
```

Terminals come from `analog()`'s argument names. Internal nodes are any
`Node` mentioned in the equations that is not a terminal (sorted by name,
so the layout is deterministic — an unsorted set made it depend on the
hash seed). State nodes are generated, one per `idt`/`idtmod`
application. Branch currents are generated, one per branch that is either
V-contributed or current-probed.

### 2.3 Contribution lowering

**Current contribution**, `Contribution(b.I, expr)` — Verilog-A's
`I(a,b) <+ expr`. The expression is added to the KCL row of `b.plus` and
subtracted from `b.minus`:

```
i[p] += expr        i[m] -= expr
```

**Voltage contribution**, `Contribution(b.V, expr)` — Verilog-A's
`V(a,b) <+ expr`. A branch-current unknown `i_b` is created; the KCL rows
carry it, and the branch row states the constitutive relation:

```
i[p] += i_b     i[m] -= i_b     i[b] += -(v_p - v_m) + expr
```

This is the standard MNA extra-equation trick — the same shape as
`elements.VS` and the `Idtmod` output row. It is what makes voltage
sources, inductors and controlled voltage sources expressible at all.

### 2.4 Term routing — the heart of the compiler

Each contribution's right-hand side is expanded and routed term by term:

| term contains | goes to | why |
|---|---|---|
| `ddt(arg)`, or `k*ddt(arg)` with `k` x-independent | `q` | the simulator integrates it, at *its* order |
| nothing x-dependent, but contains `TIME` | `u(t)` | a genuine independent excitation |
| nothing x-dependent, no `TIME` | `i` | part of the static characteristic |
| `white_noise`/`flicker_noise` | `CY` | small-signal PSD only (LRM: zero in DC/transient) |
| anything else | `i` | static current |

Three of those rows encode decisions worth defending:

**`ddt` goes to the charge vector, not to a difference quotient.** The
element never computes a time derivative; it declares a charge and the
transient engine differentiates it with whatever order and step-size
control it is using. Xyce's expression-level `ddt` computes the
derivative itself and is *hardwired to backward Euler* regardless of the
integrator in use — with the author's own comment in the source saying so
(`ast.h:6314`, `// for now, hardwire to backward Euler`, and at `:6320`
`// ERK. this isn't right`). LTspice's `ddt` silently degenerates to the
identity in `.DC`. ngspice has no `ddt` at all and the manual tells you to
build dynamics from a real capacitor plus an F-source. Routing to `q` is
the only formulation that inherits the simulator's integration order for
free.

**A state-dependent scaling of `ddt` is refused, loudly.** `g(v)*ddt(v)`
is *not* the time derivative of any charge — there is no `Q` with
`dQ/dt = g(v)·dv/dt` unless `g` is constant. Accepting it would stamp a
capacitance that conserves nothing. The compiler raises and tells you to
move the factor inside the `ddt`. No other tool surveyed performs this
check; they integrate the wrong thing quietly.

**A constant belongs to `i`, not `u`.** A diode's `-IS` is part of its
I–V characteristic, not an independent source. Putting it in `u` — as the
DSL's first version did, and as the "constant ⇒ source" rule in several
ABM implementations does — means it also appears in the *AC excitation
vector*, injecting a spurious AC current source at every bias point.
`u` here means exactly one thing: time-varying independent excitation.
(`test_hdl_element_ac_has_no_bias_leak` pins it.)

### 2.5 States: `idt` and `idtmod`

Verilog-A's integral operators cannot be expressed as algebra over the
existing unknowns, so each *application* is lowered to a new state
unknown `s` with its own equation:

```
ds/dt = arg      →      q[s] += s,   i[s] += -arg
```

and the application itself evaluates to `s` (for `idt`) or to the floored
wrap `s - m·floor((s-o)/m)` (for `idtmod`, see `idtmod.md` §1.1). This is
the same construction OpenVAF uses (an implicit equation with resistive
and reactive residual parts) and the same one `elements._IdtBase`
implements by hand.

Because the state is a real system unknown, everything the solver knows
how to do applies to it:

* **DC pinning.** With an `ic` argument, the DC variant of the residual
  replaces the state row with `s - ic`, so the operating point satisfies
  the LRM's "the initial condition shall force the DC solution". Without
  an `ic`, the DC row stays `ds/dt = arg`, which is singular for a driven
  integrator — the LRM's other branch, and the reason `uic=True` exists.
* **`uic` seeding** through the generated `state_ic()` and
  `IC_KIND='state'`.
* **Bounded states.** An `idtmod` state is declared through
  `periodic_states()`, so the transient engine's gauge shift keeps it
  inside one modulus forever (`idtmod.md` §5.2) — the generated element
  gets the ulp-exact phase property with no extra work.

`IdtmodHdl` in `elements_hdl.py` is one line of `analog()` and receives
all of the above; `test_idtmod_full_treatment` proves the DC pin, the
congruence-correct wrapped output, and the bounded state.

### 2.6 Jacobians

`G = ∂i/∂x` and `C = ∂q/∂x` are sympy Jacobians of the assembled vectors,
compiled with `lambdify(..., cse=True)`. Two consequences:

* they are **exact and closed-form** — no finite differencing, no
  per-node AD tape rebuilt every iteration;
* common subexpressions are shared **across the whole element**, because
  sympy sees `i` and `q` as complete vectors rather than as independent
  scalar expressions.

Where a stamp turns out to be x-independent (every linear element), it is
computed once and returned by reference, and the `update()` observer drops
the cache whenever parameters change. That is not an optimization
invented for its own sake: `benchmarks/hdl_overhead.py` measured a
resistor's `G` at **28× the hand-written cost** before it, because the
hand-written element returns a stored matrix while the generated one was
recomputing per Newton iteration.

---

## 3. The language

### 3.1 Statements and access functions

| DSL | Verilog-A | notes |
|---|---|---|
| `Contribution(b.I, expr)` | `I(a,b) <+ expr` | accumulates |
| `Contribution(b.V, expr)` | `V(a,b) <+ expr` | creates a branch-current unknown |
| `b.V`, `node.V` | `V(a,b)`, `V(a)` | potential probe |
| `b.I` | `I(b)` | flow probe (branch must be V-contributed) |
| `Branch(a, b)` | `branch (a,b)` | |
| `Node('mid')` | internal node declaration | discovered from use |

`analog()` returns one `Contribution` or a tuple of them.

### 3.2 Operators and functions

| DSL | Verilog-A | status |
|---|---|---|
| `ddt(x)` | `ddt` | ✅ charge routing |
| `idt(x[, ic])` | `idt` | ✅ generated state, DC pin, uic seed |
| `idtmod(x, ic, m, o)` | `idtmod` | ✅ + floored wrap + gauge shift |
| `ddx(expr, probe)` | `ddx` | ✅ exact symbolic partial |
| `limexp(x)` | `limexp` | ✅ C¹ linear continuation above 80 |
| `white_noise(p)` | `white_noise` | ✅ → `CY` |
| `flicker_noise(p, e)` | `flicker_noise` | ✅ → `CY`, `p/f^e` |
| `vt()`, `TEMP` | `$vt`, `$temperature` | ✅ bound from `epar.T` |
| `TIME` | `$abstime` | ✅ in source terms |
| any sympy function | Table 9-11 math | ✅ `exp`, `log`, `tanh`, `sqrt`, `Piecewise`, … |

`ddx` deserves a note: gnucap needs ~145 lines of AD-tape plumbing for it
and still marks flow probes `incomplete()`; here it is
`expr.subs(probe, d).diff(d)` — one line, both probe kinds, because
`Quantity` is a sympy atom. It is also, by usage count, one of the most
wanted operators in modern compact models (≈100 call sites each in the
vacask and gnucap-models device libraries).

### 3.2a Convergence aids: `limexp`, and how it relates to PCNR

`limexp(x)` is `exp` continued **linearly** above `x0 = 80`, C¹ at the
seam. It does not change the answer in the physical region and it is not
a smoothing approximation of the device: it exists so that a Newton
iterate which wanders far forward produces a *finite* residual and
Jacobian instead of `inf`. Measured on a generated diode
(`test_limexp_keeps_residual_finite_where_exp_overflows`):

| v | plain `exp`: i, G | `limexp`: i, G |
|---|---|---|
| 0.7 V | 0.0580, 2.25 | 0.0580, 2.25 (identical) |
| 3 V | 2.6e37, 1.0e39 | 2.1e23, 2.1e23 |
| 50 V | **inf, inf** | 1.0e25, 2.1e23 (finite) |

**When PCNR is enabled, it does the limiting — including for generated
elements.** pycircuit's PCNR layer (Aadithya, Keiter & Mei, *"Predictor/
Corrector Newton-Raphson"*) makes each limited quantity an explicit
unknown with its own residual row `v_Dk − (e_a − e_b) = 0`, so every
device linearises at **its own** copy and no two devices can fight over a
shared point. A device joins by supplying three things — which terminal
pair the quantity spans, its terminal currents as a function of that
quantity alone, and their derivative — because the layer *removes* the
device's ordinary `i`/`G` contribution and re-stamps it at the limited
voltage.

The DSL generates all three, plus the limiter, whenever the element
qualifies:

* exactly one current-contributed branch, no voltage contributions;
* no charge (`q ≡ 0`) — the layer refuses charge on a participant,
  since the charge term would have to move to the limited unknown too;
* no `idt`/`idtmod` states;
* the whole current is a function of that one branch voltage;
* and the current is recognisably **exponential**, with the scale read
  off the expression rather than assumed: `VT` from the argument of the
  `exp`, and `IS` from `VT·f'(0)` — which is exactly `IS` for the
  textbook `IS·(exp(v/VT) − 1)`.

Measured (`test_generated_pcnr_protocol_matches_hand_written_diode`): a
generated diode's `pcnr_i`, `pcnr_didv` and `pcnr_limit` equal the
hand-written `elements.Diode`'s to the last digit, and under PCNR the
`limexp` and raw-`exp` versions of the same model converge to *the same
answer* — the limiting is PCNR's, not the model's.

One change to the PCNR layer made this possible, and it is not an
invention — it is what the paper specifies. Aadithya, Keiter & Mei
describe the correction phase as *"requesting each device/sub-circuit to
limit the solution variables it owns"*, so that *"limiting is explicitly
invoked by the simulator, rather than being implicitly done by the
devices without the simulator's knowledge"*. pycircuit's `pcnr.refine`
nonetheless hardcoded `pnjlim` and read an `IS` parameter **by name**,
which meant only a device shaped exactly like `elements.Diode` could
ever participate. It now calls the device's `pcnr_limit` when one is
declared and falls back to that same `pnjlim` otherwise (`Diode` is
unchanged).

Two further points from the paper are worth carrying here, because they
tell you when to switch PCNR on:

* The device API PCNR wants is **stateless** — the paper's whole
  complaint about traditional limiting is that each device must track an
  internal junction voltage `v̂` across iterations, so `g` and its
  Jacobian "are no longer functions of just `x`; they also depend on the
  history of `x`". PCNR's device instead "can simply do all its
  calculations at the current iteration's `v_D`". A *generated* device is
  stateless by construction, which is why the DSL fits this interface
  better than it fits classical limiting: `elements.Diode` needs its
  `_vlim` bookkeeping and `limit()` write-back; a generated element needs
  none.
* The cost is real and bounded: the paper measures **the same number of
  NR iterations** as traditional limiting (and proves PCNR never takes
  more), at **≈1.7× the runtime**, the extra being the two phases plus
  the augmented system that the Schur complement only partly offsets.
  That is why `pcnr=False` remains pycircuit's default and why `limexp`
  still matters.

One divergence found while reconciling, now **fixed to match the paper**:
`_limiting._pnjlim` lacked Listing 1's escape, which leaves a step
unlimited when `|v_new − v_old| ≤ 2·V_T`, and compressed every step taken
above `vc` however small. Since `log(1+ε) ≈ ε` that moved no operating
point — but "nearly the identity" is exactly what the escape exists to
avoid: the paper's form makes limiting a **strict no-op** near the
solution, which is what lets a simulator use "was limiting applied?" as a
convergence signal (its footnote 3) rather than watching an O(ε²)
perturbation. Both implementations were updated together — the scalar one
and the branchless one the traced JAX/PCNR path uses — and
`tests/test_limiting.py` now checks each against a verbatim transcription
of Listing 1 and against the other.

One departure from the listing is **kept, deliberately**: our `v_new > 0`
guard. With a large `IS`, `vc` goes negative, so a negative `v_new` can
sit above it, and the listing's `v_old ≤ 0` branch then evaluates
`log(v_new/V_T)` at a negative argument — which raises. The test asserts
both halves of that: the transcribed listing raises where our guarded
version returns the step unchanged.

**Where `limexp` still earns its place:** elements that do *not* qualify
(a charge-storing junction, a state-carrying model, a polynomial or
multi-branch nonlinearity), and every run with `pcnr=False`, which is the
default. The two compose: with PCNR on, `limexp`'s clamp is simply never
reached, because PCNR keeps the iterate in range.

### 3.3 Parameters

`instparams` uses the ordinary `pycircuit.utilities.param.Parameter`;
inside `analog()` the parameter names are bound as sympy symbols, and
values are read from the *resolved* `iparv` at call time — so hierarchical
string expressions (`r='2*rbase'`) reach the generated code.

> `hdl.Parameter` used to mix in `sympy.Symbol`. sympy **interns symbols
> by name**, so two elements each declaring a parameter named `c0` shared
> one object and whichever class was defined last silently set the default
> for both. The mixin bought nothing (the compiler builds its own symbols)
> and is gone; `hdl.Parameter` is now an alias for the plain
> `Parameter`.

---

## 4. Where this sits in the landscape

Behavioral modelling across simulators falls into four archetypes:

| archetype | examples | dynamics | Jacobian |
|---|---|---|---|
| static expression source | ngspice B/ASRC | **none** (use an external C) | symbolic diff of the parse tree, cached at parse |
| rich expression source | Xyce B/E/G, HSPICE/PSpice E/G/F/H, LTspice B | `ddt`/`sdt` (Xyce, BE-hardwired); `ddt`/`idt`/`idtmod` (LTspice); `LAPLACE`/`FREQ` (HSPICE) | forward AD over the AST (Xyce); analytic, undocumented (others) |
| per-branch I(V)+Q(V) device | **Qucs EDD** | charge-based `Q(V,I)` | symbolic diff of the equation AST at setup, ≤20 branches |
| compiled HDL | Verilog-A via Spectre/OpenVAF/VACASK; XSPICE code models; gnucap `bm_*` | full operator set | compiler AD, or hand-written partials |

`Behavioural` is an **archetype-2/3 hybrid with an archetype-4 surface
syntax**. Its closest relatives are Qucs EDD structurally (both split
static current from charge and differentiate symbolically at setup) and
Xyce ABM on Jacobian quality — but with these differences:

1. The `i`/`q` split is **derived** from `ddt` in a single contribution,
   not hand-separated by the user into `I<n>` and `Q<n>` equations.
2. `ddt` inherits the **simulator's** integration order (§2.4).
3. A non-conservative `ddt` scaling is **rejected** rather than silently
   mis-integrated.
4. `idt`/`idtmod` become **real DAE states**, not hidden accumulators —
   which is why they work at DC and under `uic`, where LTspice's `idt`
   famously does not.
5. There is a **JAX autodiff path**: generated `eval_i_pure`/
   `eval_q_pure` admit the element to vectorised, batched, per-lane
   parameter sweeps. No simulator in the survey vectorises behavioral
   evaluation across devices.
6. The DSL is **Python and sympy**, so parameters are symbols and the
   model can be differentiated with respect to *parameters*, not just
   states.

**And where compiled Verilog-A still wins** (§6): events, `$limit`,
`laplace_*`, `transition`/`slew`, and the general caveat every archetype
1–3 tool shares — a contribution expression must be valid for every value
the solver might try, not merely for physical ones, and an infinite-slope
transition will produce timestep-too-small failures with no built-in
smoothing escape hatch.

---

## 5. Cost — measured

`benchmarks/hdl_overhead.py`, this machine, 2026-08-22:

| call | hand-written | HDL | ratio |
|---|---|---|---|
| `R.i` | 0.53 µs | 1.30 µs | 2.5× |
| `R.G` | 0.07 µs | 0.55 µs | 8.2× |
| `C.q` | 0.76 µs | 0.78 µs | 1.0× |
| `C.C` | 0.26 µs | 0.11 µs | **0.4×** |
| `Diode.i` | 1.88 µs | 1.78 µs | **0.95×** |
| `Diode.G` | 1.61 µs | 2.44 µs | 1.5× |

Compile cost: **6–9 ms per element, once**, at class definition;
importing all ten of `elements_hdl` costs **136 ms** with the package
already loaded. That is the price of the exactness — symbolic assembly,
two Jacobians, `lambdify` with CSE — and it is why `elements.py` stays
the default catalogue while `elements_hdl.py` is imported on demand.

End-to-end, an 8-stage RC ladder transient (1000 steps, identical
waveforms to 1.1e-15): **0.62 s hand-written vs 0.70 s HDL — 1.14×**.

Read it this way: for a *nonlinear* element the generated code is at
parity or better, because both sides evaluate the same transcendental and
the generated one is CSE'd. The remaining gap is on *trivial linear*
stamps, where the hand-written element does almost nothing at all
(`return self._G`) and any function-call overhead dominates a 2×2 matrix.
The 8× on `R.G` is 0.5 µs of Python call machinery, not expression
walking.

---

## 6. Capability map against Verilog-A

Priorities below come from call-site counts across the vacask device
library (59 models) and gnucap-models (19), plus gnucap's own
`vams/*.vams` — SPICE primitives written *as* Verilog-A elements, which is
precisely this DSL's target.

### Implemented

`I <+`, `V <+`, accumulation, `V`/`I` probes, internal nodes, named
branches, `ddt`, `idt`, `idtmod`, `ddx`, `limexp`, `white_noise`,
`flicker_noise`, `$temperature`/`$vt`, `$abstime`, parameters (with
hierarchical expressions), analog helper functions (any Python function
returning a sympy expression — the LRM's restrictions make them
inlinable by construction), `Piecewise` conditionals, exact `G`/`C`,
`dudt`, `CY`, DC-variant stamps, `state_ic`/`IC_KIND`, `periodic_states`,
`eval_*_pure`.

### Not implemented, ranked by measured demand

| capability | class | usage (vacask/gnucap-models) | note |
|---|---|---|---|
| `$limit` | HARD | 273 / 0 | needs previous-iterate memory **and** a "not converged" channel. **PCNR participation is generated instead** (§3.2a) for qualifying elements, which is the better answer where it applies |
| `$param_given` | TRIVIAL | 1871 / 263 | most-called system function in real models; needs a "was it defaulted" flag on `Parameter` |
| unconditional node collapse (`V(a,b) <+ 0`) | TRIVIAL | common idiom | merge the node symbols before building `x` |
| `@(initial_step)` | TRIVIAL | 2 / 0 | a phase flag |
| `@cross`, `@timer` | EVENT | 0 / 0 | `elements._WrapEvents` already implements the contract (`accept_step` + `next_event`) — promote when a comparator/oscillator macromodel needs it |
| `absdelay` | HARD | 4 / 0 | history buffer; `elements.TLine` covers the real use |
| `transition`, `slew` | EVENT | 0 / 0 | gnucap's `slew` is self-described as a stub |
| `laplace_*`, `zi_*` | STATE | 0 / 0 | gnucap ships only 2 of 4 `zi_*`, with two open bugs |
| `last_crossing`, `noise_table` | EVENT/— | 0 / 0 | gnucap does not implement either |
| conditional node collapse / switch branches | HARD | 1 idiom | changes sparsity per iteration; diagnose and reject rather than half-implement |
| AC excitation from HDL | small | — | `u` is zeroed for AC; a behavioral AC source needs an `ac`-variant vector |
| symbolic-toolkit transient | small | — | compilation targets numpy/jax; AC via `G`/`C` works |

The shape of that table is the point: **eleven of gnucap's fourteen
Verilog-A-written SPICE primitives need nothing outside what is
implemented.** The two that do — `vpulse` (`transition`, `@timer`) and
`tline` (`absdelay`) — already exist as hand-written pycircuit elements.

---

## 7. Usage

### Writing an element

```python
import sympy
from pycircuit.utilities.param import Parameter
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   ddt, idtmod, limexp, white_noise,
                                   vt, TIME, TEMP)

class MyDiode(Behavioural):
    instparams = [Parameter(name='IS', desc='Sat. current', unit='A',
                            default=1e-13)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return (Contribution(b.I, IS * (limexp(b.V / vt()) - 1)),
                Contribution(b.I, white_noise(2 * 1.38e-23 * TEMP * IS)))
```

Then use it like any element: `c['D'] = MyDiode('a', 'b', IS=2e-14)`.

### A VCO in four lines

```python
class VCO(Behavioural):
    instparams = [Parameter(name='kvco', desc='Gain', unit='Hz/V',
                            default=1e6),
                  Parameter(name='fc', desc='Centre', unit='Hz',
                            default=1e6)]

    @staticmethod
    def analog(ctrl, gnd_, outp, outn):
        bc, bo = Branch(ctrl, gnd_), Branch(outp, outn)
        phase = idtmod(fc + kvco * bc.V, 0, 1, -0.5)
        return Contribution(bo.V, sympy.sin(2 * sympy.pi * phase))
```

The bounded phase, the DC pin, the wrap breakpoints and the gauge shift
all come from the `idtmod` lowering.

### Reading the generated code

```python
el = MyDiode('a', 'b'); el.update_iparv()
el.i(x)        # generated numpy function
el.G(x)        # exact symbolic Jacobian, compiled
el.linear      # False -- inferred from x-dependence of G
```

`elements_hdl.py` is the worked catalogue: ten elements, each exercising
one capability, each proven against its hand-written twin.

---

## 8. Testing and provenance

* `pycircuit/circuit/tests/test_elements_hdl.py` — 28 tests. Stamp
  equivalence (R, G, C, VCCS, Diode) at arbitrary operating points;
  waveform equivalence in real circuits (RC lowpass, RL, VCVS DC, HDL
  source, `Idt`, full `Idtmod` treatment); the generated Jacobian against
  central differences; noise against `elements.R.CY`; `ddx` and `limexp`;
  the two historic compiler crashers; the `ddt`-scaling refusal; the
  constant-stamp cache following parameters; the AC bias-leak property;
  JAX pure-form parity.
* `pycircuit/circuit/tests/test_hdl.py` — the original two tests, still
  green, unchanged.
* `benchmarks/hdl_overhead.py` — §5's numbers, regenerated on demand;
  it asserts the two circuits take the same steps and agree to 1e-12, so
  the timings compare equal work.

### Defects found and fixed while building this

Each was caught by a test that now guards it:

1. **`Contribution(b.V, ...)` returned `None`** and the caller unpacked
   it — voltage sources, inductors and controlled voltage sources were
   structurally impossible.
2. **`c*ddt(v)` was misclassified** as static current, producing an
   unevaluated `Derivative` that killed `lambdify`.
3. **`ddt(c*tanh(v))` crashed at class creation** — sympy's printer
   dispatches on the class *name*, hit `_print_Quantity` (written for
   `sympy.physics.units.Quantity`) and demanded a `.name`.
4. **Method signatures did not match the element protocol**, so no
   `Behavioural` element could be placed in a `SubCircuit` at all.
5. **Parameters were read from `ipar`, not `iparv`** — hierarchical
   expressions never resolved.
6. **`hdl.Parameter` interned across elements** via `sympy.Symbol`
   (§3.3).
7. **Internal nodes were compiled but never registered**, so the
   generated vectors were longer than `element.n`.
8. **`CY` was a hardcoded zero matrix** — every generated element was
   silently noiseless.
9. **Bias constants leaked into `u`**, i.e. into the AC excitation.
10. **`linear` stayed `True`** for nonlinear models, hiding them from
    distortion analysis.
11. **The DSL's atoms carried no `real` assumption**, so sympy refused to
    build any inequality over a voltage — every `Piecewise` died at class
    creation with "Invalid comparison of non-real". That made `limexp`
    and *all* conditional models unusable inside an element. It survived
    the first round of tests because those exercised `limexp` only as a
    standalone sympy function, with no `Quantity` in the expression; it
    surfaced when the operator was first put into a circuit alongside
    PCNR.
12. **A subclass changing only `instparams`** inherited its parent's
    compiled code — which reads a different parameter list — and would
    have answered confidently with the wrong numbers. Now refused with
    an explanatory `TypeError`.


---

## 9. Roadmap

Ordered by (measured demand) × (what it unblocks) ÷ (cost). "Demand" is
call sites across the vacask device library (59 models) and
gnucap-models (19), counted in §6's survey; a capability with zero call
sites in either is deferred on purpose and says so.

Every item states what it touches, how it would be proven, and its risk.
Two conventions from §8 apply to all of them and are not repeated per
item: **a new operator is not done until a test exercises it inside a
circuit** (defect 11 passed a standalone-function test and was still
fatal in every element), and **the compiler must refuse what it cannot
do, per element, rather than emit something plausible** (§3.2a's
qualification rules are the pattern).

### Phase A — convergence: finish what PCNR started

**A1. Widen PCNR qualification.** *This is the most valuable item in the
list and the one most likely to be mis-ranked, because it looks like
infrastructure rather than a feature.* §3.2a generates the PCNR protocol
only for a single-branch, charge-free, state-free, exponential element.
Every real junction device fails at least one of those: a diode with
junction capacitance has charge, a BJT has two junctions sharing a base,
a MOSFET has four terminals. The blockers are in the PCNR layer, not the
DSL:

* **charge.** `pcnr.augmented_system` raises when a participating device
  has any `C` entry, and the refusal is correct as written: the charge
  term is evaluated at the node voltage while the current is evaluated
  at `v_lim`, which is exactly the inconsistency PCNR exists to remove.
  The fix is to move the charge to the limited unknown too — `q(v_lim)`
  and `dq/dv_lim` alongside `pcnr_i`/`pcnr_didv` — so the companion
  terms are formed from the same quantity. The paper is explicit that
  PCNR "works for differential-algebraic equations as well" but
  (footnote 1) develops only the algebraic case, so the DAE bookkeeping
  is ours to derive and must be tested against a transient, not just a
  DC point.
* **multiple junctions per device.** The layer already carries a
  sequence of `(anode, cathode)` pairs per device and `limit_junctions`
  already handles the shared-terminal case with its `move` index (a
  BJT's two junctions share the base, so limiting both by adjusting the
  anode has the second undo the first). What is missing is generating
  those pairs from an expression with two or more exponentials in
  different branch voltages.
* **non-exponential nonlinearity.** `pcnr_limit` currently derives
  `(VT, IS)` by reading the `exp` argument. For anything else there is
  no principled scale, and inventing one would be the kind of unvalidated
  heuristic §7 rejected. Options, in order of honesty: leave it
  unqualified (today's behaviour); let the *user* declare a limiter in
  the element; derive a local scale `f/f'` and prove it on a benchmark
  before shipping it.

  Proof: extend `test_elements_hdl.py`'s protocol-equality test to a
  charge-storing diode against a hand-written reference, plus a transient
  (not just DC) where PCNR and non-PCNR runs agree; a two-junction
  element against `elements`' BJT-shaped limiting. Risk: medium-high —
  it changes shared solver code that the hand-written devices also use,
  so the existing PCNR tests are the regression gate.

**A2. `$limit`.** 273 call sites, the highest raw demand in the survey,
and the fallback for every model A1 cannot reach. Needs two things the
DSL has no route to today: memory of the previous iterate, and a channel
to tell the solver "limiting fired, so this iteration is not converged"
(the paper's footnote 3). `elements.Diode._vlim` and `VCVS_limited` show
the hand-written shape; `SubCircuit.limit` already dispatches. The
awkwardness is that a generated element is otherwise *stateless*, which
is precisely what makes it fit PCNR — so `$limit` should be implemented
as the non-PCNR path and documented as the lesser option, not as the
default. Proof: a model that converges with `$limit` and demonstrably
does not without it (the measurement that `limexp` alone could not
produce here — see §3.2a, where pycircuit's continuation ladders rescued
plain `exp` anyway). Risk: medium; the honest failure mode is shipping a
capability whose benefit cannot be measured on this simulator's circuits.

### Phase B — model surface: cheap, high-frequency, boring

**B1. `$param_given` + parameter ranges + `aliasparam`.** By call count
the most-used system function in the entire survey (**1871** vacask /
**263** gnucap-models), and trivial: `Parameter` already carries
`desc`/`unit`/`default`, so `$param_given` is a "was this defaulted"
flag set when `ipar.set()` names it, ranges are a validator on
`update_iparv`, and `aliasparam` is a name map. Touches
`pycircuit/utilities/param.py` and the DSL's parameter binding. Proof:
a model that branches on `$param_given` compiling to two different
stamps; an out-of-range parameter refused with the range in the message.
Risk: low — but it touches `param.py`, which every element uses, so the
regression gate is the whole suite.

**B2. Unconditional node collapse.** `V(a,b) <+ 0` should merge the two
nodes and delete the unknown, which is the standard idiom for an
optional series resistance (`if (rs) I(a,ia) <+ ...; else V(a,ia) <+ 0;`).
Adopt gnucap's restriction to *unconditionally executed* contributions
(`mg_in_analog.cc:2065-2073`): a conditional collapse changes the
sparsity pattern per iteration and is Phase D. Implementation is a
substitution before the unknown vector is built — merge the node symbols,
drop the row. Proof: an element with `rs=0` producing exactly the stamp
of the element without the branch, and node count reduced by one. Risk:
low, and it is pure DSL work.

**B3. AC excitation.** `u` is currently zeroed for `analysis='ac'` so a
device's bias constants cannot leak into the small-signal drive (§2.4).
The missing half is a *deliberate* AC source: an `ac`-variant vector, so
a behavioural element can be an AC stimulus rather than merely
transparent to one. Follows `VS.u`/`IS.u`'s analysis switch. Proof: an
HDL source driving an AC analysis to the analytic transfer function.
Risk: low.

### Phase C — analyses the DSL cannot currently reach

**C1. Events: `@cross`, `@timer`.** Zero call sites in either corpus, so
this is *not* ranked on model demand — it is ranked on what it unblocks
here: comparator, switch and oscillator macromodels, which are exactly
the behavioural models a DSL is for. The machinery exists:
`elements._WrapEvents` already implements the `accept_step` +
`next_event` contract with the same first-order crossing prediction
gnucap uses, and `SubCircuit.next_event` already polls. The work is
generating a predictor for an arbitrary expression rather than the
idtmod wrap. Proof: a comparator whose output edge lands on the crossing
to solver tolerance, and a step-count comparison against the same model
without the event. Risk: medium — event prediction interacts with the
step controller, and §5.3 of `idtmod.md` records how that went for the
wrap breakpoints.

**C2. Symbolic-toolkit transient.** Compilation targets numpy/jax, so a
symbolic transient of a generated element is out (AC works, since it
uses `G`/`C` only). The fix is to keep the sympy expressions rather than
lambdifying when the toolkit is symbolic — the expressions already exist,
so this is plumbing, not derivation. Proof: the existing symbolic tests
extended to a generated element. Risk: low. Value: mostly for the
distortion/DDD analyses, which is where this repo's symbolic work lives.

### Phase D — deferred, with the reason

Not "someday" — these have a stated trigger, so the deferral can be
argued with:

| capability | why deferred | what would reopen it |
|---|---|---|
| `laplace_*`, `zi_*` | 0 call sites in both corpora; gnucap ships only 2 of 4 `zi_*` and documents two open bugs in them | a user needing a filter macromodel; an explicit state-space element is probably the better answer even then |
| `absdelay` | 4 sites, all transmission lines, and `elements.TLine` already covers that | a delay model that is not a T-line |
| `transition`, `slew` | 0 real sites; gnucap's own `slew` is a self-described stub; `VPulse` exists | a digital-to-analog boundary model |
| `last_crossing`, `noise_table` | 0 sites, and gnucap implements neither | demand |
| conditional node collapse / switch branches | changes sparsity per iteration; needs gnucap's whole per-iteration re-stamping machine | a model that genuinely needs V and I on one branch under a condition — diagnose and refuse clearly until then |
| `$discontinuity` | 41 sites, but it is advisory and gnucap's implementation is an empty body | nothing; **parse and ignore** is the cheapest correct answer and unblocks those 41 sites today |

### Sequencing

B1 → B2 → B3 first: they are cheap, they are the highest-frequency
things real models use, and none of them touch the solver. Then A1,
which is the item that decides whether this DSL can express *production*
device models or only behavioural ones — and which should be attempted
before A2, because if PCNR can be widened far enough, `$limit` matters
much less. C1 and C2 follow demand. Phase D stays deferred until its
trigger fires.

## 10. References

**Local repositories** (`~/source`, evidence for §4 and §6)

- ngspice: `src/spicelib/parser/inpptree.c:231-234` (derivative trees
  built at parse), `:255` (`PTdifferentiate`), `:473-553` (naive
  derivatives for step/round/min/max), `doc/ngspice.texi:3021-3086`
  (function list; the "integrate a current into a capacitor" idiom);
  XSPICE `src/xspice/cmpp/mod_yacc.y:463-468` (`PARTIAL()`).
- Xyce: `src/UtilityPKG/ExpressionSrc/ast.h:241` (`dx`), `:6092`
  (`sdtOp`), `:6268-6320` (`ddtOp`, hardwired BE and the author's note),
  `doc/Users_Guide/Xyce_UG_ABM.tex:344-380` (expressions must be valid for
  non-physical iterates), `:427-440` (infinite-slope failures),
  `doc/Reference_Guide/B_Device.tex:80-93` (`smoothbsrc`).
- Qucs-S EDD: `qucsator_rf/src/components/devices/eqndefined.cpp:196-265`
  (G and C by symbolic differentiation, incl. the `dQ/dI·dI/dV` chain),
  `:470-511` (the 20-branch cap), `src/equation.cpp:988-996` (silent
  fallback when no derivative rule exists).
- gnucap: `gnucap-modelgen-verilog/mgvams/mg_filt_xdt.cc:34-63`
  (`ddt`/`idt` as filter elements with state), `src/mg_func_evt.cc:134`
  (`@cross`), `src/mg_task_limit.cc:36-117` (`$limit` and its
  `set_converged(false)` channel), `mgvams/mg_func_math.cc:491`
  (`limexp` as ordinary math), `mgvams/mg_filt_ddx.cc:70-145` (`ddx` via
  the AD tape; flow probes `incomplete()`), `src/mg_in_analog.cc:1020`,
  `:2065-2073` (node collapse, restricted to unconditional
  contributions), `bm/bm_tanh.cc:101-118` (hand-written value+slope).
- VACASK: `lib/rpnexprva.cpp:114-159` (b-sources translated to Verilog-A
  and handed to OpenVAF).

**Standards and literature**

- Verilog-AMS LRM 2.4.0 —
  <https://www.accellera.org/images/downloads/standards/v-ams/VAMS-LRM-2-4.pdf>
- Aadithya, K.V.; Keiter, E.R.; Mei, T., *"Predictor/Corrector
  Newton-Raphson (PCNR): A Simple, Flexible, Scalable, Modular, and
  Consistent Replacement for Limiting in Circuit Simulation"*, Sandia
  National Laboratories; Springer 2020,
  <https://doi.org/10.1007/978-3-030-44101-2_19> (OSTI 1523781) — the
  method pycircuit's `pcnr.py` implements and that §3.2a's generated
  device protocol targets. Its Fig. 2 gives the predict/correct flow and
  the Schur-complement elimination; Table 1 the cost (same NR iteration
  count, ≈1.7× runtime).
- Kundert & Zinke, *The Designer's Guide to Verilog-AMS* —
  <https://designers-guide.org/verilog-ams/dg-vams/index.html>
- HSPICE behavioral modelling ch. 22 —
  <https://class.ece.uw.edu/cadta/hspice/chapter_22.pdf>
- PSpice ABM and convergence —
  <https://www.pspice.com/resources/application-notes/analog-behavioral-modeling>,
  <https://resources.pcb.cadence.com/pspiceuserguide/b-convergence-and-time-step-too-small-errors>
- LTspice B-source reference —
  <https://ltwiki.org/?title=B_sources_(complete_reference)>
- Spectre reference (bsource, AHDL/CMI) —
  <https://ee.kpi.ua/~yv/edu/ok/book/spectre_refManual.pdf>

**In this repository**

- `idtmod.md` — the `idt`/`idtmod` research record: LRM semantics, the
  simulator survey, the gauge shift, Baumgarte orbit correction.
- `doc/src/circuit/idtmod.rst` — the user-facing circular-integrator page.
