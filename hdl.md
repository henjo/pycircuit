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
| `Branch(a, b)` | `branch (a,b)` | anonymous |
| `Branch(a, b, 'ch')` | `branch (a,b) ch;` | a **distinct** branch across the same pair, with its own current |
| `simparam(name, default)` | `$simparam(name, default)` | resolved at compile time |
| `Collapse(br, when)` | `if (R>0) I<+…; else V<+0;` | node collapse on a **parameter** condition |
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

**Where `limexp` still earns its place** — and an earlier version of this
section got this wrong, so the correction is the interesting part. It said
that with PCNR on, `limexp`'s clamp "is simply never reached, because PCNR
keeps the iterate in range". Measurement says the opposite, and the
mechanism matters: `pcnr.augmented_system` assembles `cir.i(x)` — which
*includes* the participant — and then subtracts that device's own `i(sub)`
again, since its current is about to be re-stamped at `v_lim`. The
cancellation is exact in arithmetic and worthless in floating point once
the term is `inf`: **`inf − inf = nan` poisons the whole system.** PCNR
bounds the *limited quantity*; it does not bound the node voltage at which
the device's own `i()` is evaluated during assembly.

Measured on a 20 V, 1 Ω forward drive: the raw-`exp` generated diode
gives `i = inf` at the trial point and PCNR fails outright, while the
`limexp` version stays finite (3.9e24) and converges to the same 0.849861
the non-PCNR path finds. So the guidance is the reverse of what was
written: **use `limexp` in a model you intend to run under PCNR**, not
merely in one you don't. `test_limexp_is_what_makes_a_pcnr_participant_robust`
pins it. `limexp` is also still the aid for elements that do not qualify
at all (a state-carrying model, a polynomial or multi-branch
nonlinearity) and for every `pcnr=False` run, which is the default.

### 3.2b The let-chain: `var()`, and why a compact model needs it

Everything in §3.2 assembles one expression per contribution and hands it
to `Matrix.jacobian`. That is the right thing for a diode, and it is
structurally incapable of compiling a surface-potential model.

The reason is not model size. It is **reuse**. A compact model computes a
quantity once and mentions it several times:

```
u1 = sqrt(1 + u0*u0) + u0*tanh(u0)/2
```

`u0` appears three times. Substitute `u0`'s own definition in and you get
three copies of it; substitute again and you get nine. After *d* levels
the expression **tree** holds 3^d copies of the innermost subexpression,
while the **DAG** has only *d* nodes. sympy works on the tree.

Measured (`benchmarks/hdl_chain_scaling.py`, table 1) on exactly that
shape:

| depth | eager | growth vs previous | let-chain |
|---|---|---|---|
| 1 | 0.04 s | — | 0.05 s |
| 2 | 0.08 s | 1.9× | 0.02 s |
| 3 | 0.18 s | 2.4× | 0.02 s |
| 4 | 0.58 s | 3.2× | 0.02 s |
| 5 | 16.2 s | **27.9×** | 0.03 s |
| 16 | *not run* | | 0.07 s |
| 48 | *not run* | | 0.22 s |

The growth ratio is itself growing, because sympy's operations are
superlinear in tree size — so the wall arrives sooner than 3^d predicts.
PSP103 carries roughly 1400 intermediates. There is no version of the
eager path that reaches it.

**`var(expr, name)`** names an intermediate and returns a symbol standing
for it:

```python
from pycircuit.circuit.hdl import var

def analog(d, g, s, b):
    bds, bgs = Branch(d, s), Branch(g, s)
    vg   = var(bgs.V / vt(), 'vg')
    surf = var(sympy.sqrt(1e-9 + vg*vg), 'surf')
    chan = var(surf * sympy.tanh(bds.V / (1 + surf)), 'chan')
    return Contribution(bds.I, IS * chan)
```

Nothing is substituted. The compiler emits a topological **let-chain** —
one assignment per intermediate, in dependency order — and obtains the
Jacobian by **forward accumulation** over the chain: each definition's
gradient is written in terms of the gradients of the definitions it
actually mentions, so the work is proportional to the number of DAG
*edges*. Definitions no output reaches are pruned, so `i`, `q`, `u` and
the AC stamp each compile only their own part of the chain.

At PSP scale (table 2 of the same benchmark — 4 terminals, 3
contributions, and the sqrt/exp/ln/tanh kernel):

| intermediates | compile | emitted lines | `i` eval | `G` eval | ‖J − FD‖ |
|---|---|---|---|---|---|
| 200 | 3.3 s | 934 | 59 µs | 0.64 ms | 1.0e-16 |
| 600 | 11.6 s | 2932 | 180 µs | 2.2 ms | 2.9e-17 |
| 1400 | 34.0 s | 6932 | 400 µs | 5.2 ms | 1.0e-19 |

Compile time is paid once per model class. Every Jacobian is checked
against finite differences.

**Two things the classification has to get right.** `_split_terms` sorts
a contribution's terms into `i`, `q` and `u` by asking which symbols they
touch — and an intermediate is an *opaque symbol*. So the compiler tracks,
transitively, which intermediates depend on the solution and which carry
`TIME`, and feeds both sets into the split. Without the first, a
conductance written through `var()` is filed as a constant source and the
device has no conductance at all; without the second, a time-varying
source behind a `var()` is filed as a static characteristic and the source
is dead. `test_hdl_chain.py::TestTermClassification` pins both.

**What the chain gives up, it gives up loudly.** A chained model reports
`pure_spec = None` (no JAX `eval_*_pure`), `sym_spec = None` (no symbolic
toolkit), `const_G = const_C = False` (no constant-stamp caching) and
`pcnr_funcs = None`. That last one matters: PCNR's shape detector reads
the exponential straight out of the expression to recover `IS` and `VT`,
and behind a `var()` it would find none — which is indistinguishable from
a linear device. The compiler declines PCNR explicitly rather than
silently registering the device as having nothing to limit.
`test_hdl_chain.py::TestCapabilityRefusals` pins each of these.

The path is chosen per model: declare no intermediates and nothing
changes. `_hdl_info['chained']` says which path a class took.

### 3.2c The math kernel, and both-arms-safe conditionals

A compiled conditional evaluates **both arms** — sympy's `Piecewise`
becomes numpy's `select`, which picks afterwards. So the obvious guard

```python
sympy.Piecewise((sympy.sqrt(x), x > 0), (0.0, True))
```

still takes the square root of the negative number. The *value* survives
(the bad arm is discarded), but two things do not:

- the floating-point flag is raised on **every** evaluation, and under
  `np.seterr(all='raise')` — which anyone debugging a convergence failure
  will set — the model becomes unusable at ordinary biases;
- the **derivative** frequently does not survive. sympy differentiates
  the clamped form `sqrt(max(x,0))` to `d(max)/dx / (2*sqrt(max(x,0)))`,
  which is `0/0` at every `x ≤ 0`. Measured: 2013 non-finite derivatives
  over a 4001-point sweep. One NaN in the Jacobian loses the whole Newton
  step, not one entry.

PSP103 has roughly 70 bias-dependent conditionals of exactly this kind.
`hdl.py` therefore ships the kernel PSP is written in, with every arm
valid for every input:

| function | what it is | note |
|---|---|---|
| `expl(x)` | two-sided clamped `exp` | PSP's: `exp` inside ±`se05` = 230.2585, its own 3rd-order Taylor outside, **C-3** at both seams |
| `expl_low(x)` | lower tail only | stays **strictly positive**; plain `exp` underflows to exactly 0 below −745 |
| `expl_high(x)` | upper tail only | the upper half of the same |
| `hypsmooth(x, eps)` | smooth `max(x, 0)` | no branch at all — the right first thing to reach for |
| `safe_sqrt(x)` | `sqrt(hypsmooth(x))` | finite derivative below zero |
| `safe_ln(x)` | `ln(hypsmooth(x))` | no `−inf` at zero |
| `safe_div(a, b)` | `a·b / (b² + ε²)` | finite at `b = 0` |

Three findings from building these are worth stating, because each one
broke a version that looked correct:

**Floating point destroys `hypsmooth`'s whole point if you write it
literally.** `0.5*(x + sqrt(x*x + 4ε²))` is strictly positive in exact
arithmetic. At `x = −100`, `ε = 1e-12`, `sqrt(x*x + 4e-24)` rounds to
exactly `100.0`, the sum cancels to exactly `0.0`, and every downstream
consumer promised a positive number gets a zero — this is what turned
`safe_ln` into 2012 infinities. The negative side is therefore evaluated
through the conjugate form `2ε²/(sqrt(x*x + 4ε²) − x)`, algebraically
identical but a sum of positives with nothing to cancel. Each arm is
additionally clamped to *its own* side, or the conjugate arm evaluated at
large positive `x` has `r − x` cancel to zero and divides by it — the
mirror image of the bug it was written to avoid.

**`sign` is not usable in a guard.** The natural regularisation
`a/(b + ε·sign(b))` differentiates to a `DiracDelta`, which no numeric
backend prints — so the model fails to *compile*, not to converge. Hence
the `a·b/(b² + ε²)` form.

**`limexp` is deliberately left unsafe.** Clamping its argument turns
`exp(v/vt)` into `exp(min(v/vt, 80))`, and PCNR's shape detector no longer
recognises a junction — the device silently loses its limiting. Since
PCNR is exactly what bounds the argument, the overflow in the discarded
arm is the better half of that trade.
`test_hdl_kernel.py::TestLimexpIsDeliberatelyNotSafe` pins it.

**`limexp` and `expl_high` are different functions, not two spellings.**
`limexp` is the LRM's: continue *linearly* from 80. PSP's `expl_high`
continues with exp's own third-order Taylor from `se05` = 230.2585. They
agree below 80 and diverge enormously above: at `x = 100`, `limexp`
gives `exp(80)·21` where `expl_high` gives `exp(100)` exactly — seven
orders of magnitude apart. Use `limexp` for the LRM's function and the
`expl` family when translating a PSP-family model.

> **Correction (2026-08-22).** The `expl` family originally shipped here
> with a threshold of 80 and linear/hyperbolic continuations, documented
> as "PSP's `expl`". That was wrong on both counts. Reading the shipped
> sources — `Common103_macrodefs.include:68-89`, and identically in the
> PSP-based `mosvar` — PSP uses `se05 = ln(1e100) = 230.2585` and
> continues with `P3(u) = 1 + u + u²/2 + u³/6`, exp's own Taylor series
> about the seam, which is what makes the family **C-3** continuous
> rather than merely C¹. The 80 appears to have been carried over from
> `limexp`'s conventional threshold. Now corrected to match the source;
> the both-arms-safe clamping is kept, since the vendor's own arms do
> overflow. PSP103 and JUNCAP200 call the family **31** times between
> them.

Underflow is not guarded anywhere: `0.0` is the correctly-rounded answer
for `exp(-1000)`, and it is the one benign IEEE condition. Overflow,
invalid and divide-by-zero all mean an arm produced garbage, and those are
what the kernel exists to stop. Above `|x| ≈ 1e154` the radicand of
`hypsmooth` overflows; circuit quantities do not reach there, and buying
that range back would cost a comparison on every evaluation of the most
frequently called function in the kernel.

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

* **charge — DONE 2026-08-22, and the plan above was wrong about why.**
  `pcnr.augmented_system` raised whenever a participant had any `C`
  entry, on the stated grounds that leaving the charge in the MNA block
  reintroduces "the exact inconsistency PCNR exists to remove". That
  reasoning conflated two different things, and the refusal cost every
  junction device with capacitance — which is to say every real one.
  What PCNR removes is a **clash between devices** over a shared
  linearisation point; a charge is limited by nobody, so it has no clash
  and no owner to fight over. And the Newton system is **exactly
  consistent** as it stands: `g_MNA`'s charge part is a function of
  `x_MNA` with its derivative `Geq` in `J_MNA/MNA`, while the device's
  current is a function of `v_lim` alone with its derivative in
  `J_MNA/lim` — `J` is `dg/dx` for both blocks. The refusal is gone; the
  claim is finite-differenced in `tests/test_pcnr_charge.py`, which also
  checks a charge-storing participant's DC and transient against
  `pcnr=False`. What is genuinely given up is that the reactive part is
  not limited along the iteration path, which costs nothing at the
  answer: convergence already requires `|g_lim|` below tolerance, i.e.
  `v_lim` equals the node voltage the charge was evaluated at. Moving
  the charge to the limited unknown as well (`q(v_lim)`, `dq/dv_lim`)
  remains possible but is now an optimisation of the iteration path, not
  a correctness fix — and the paper does not derive it either (footnote
  1: PCNR "works for differential-algebraic equations as well, but for
  simplicity, we only consider algebraic equations").
* **multiple junctions per device — DONE 2026-08-22.** The blocker was
  in the protocol, not the generator: `pcnr_i(v, params, epar, toolkit)`
  had no way to say *which* junction was being asked about, so a device
  owning two of them could not be written even by hand. The layer now
  passes a per-device junction index `jn` (computed by counting each
  instance's pairs in declaration order; single-junction devices never
  see anything but 0 and may ignore it), and the DSL detects junctions
  **per contribution** rather than by reverse-engineering the assembled
  `ivec` — each `I`-contribution that is a function of its own branch
  voltage alone and carries a single exponential scale becomes one
  limited quantity. Measured on a two-junction element sharing a base
  (the BJT shape, and the clash case PCNR exists for): both junctions
  declared, per-junction currents differing exactly as their saturation
  currents do, and DC identical with PCNR on and off.

  **A silent wrong answer found while doing it, and then the cause
  removed.** The traced JAX PCNR path did not call the device at all — it
  rebuilt every junction as `IS*(exp(v/VT)−1)` with one global VT,
  reading `IS` *by name*. A device whose saturation current is not called
  `IS` (the two-junction element has `ISE`/`ISC`) got `IS = 0`: a
  junction carrying **no current**, and a confident wrong answer.

  The first fix was to verify and refuse. The **second was to stop
  guessing altogether**: the traced loop now asks each device for its own
  junction current and conductance, through the same `pcnr_i`/`pcnr_didv`
  the CPU path calls, with the generated elements emitting a jax-printed
  twin. The junction set is static at trace time, so the per-junction
  loop unrolls under `jit`. That removed **all three** of the traced
  backend's restrictions at once — arbitrary junction shape, several
  junctions per device, and charge storage (the assembly subtracts the
  junction term at the node voltage and adds it at `v_lim`; it never
  touched the charge) — closing the asymmetry where the CPU accepted
  models the traced backend refused. Measured: a two-junction element
  agrees with the CPU to six digits, a charge-storing one to 0.3% on a
  2.6 V swing. The verify-and-refuse logic survives as the *fallback*,
  for a device that supplies no traceable evaluation.
* **non-exponential nonlinearity.** `pcnr_limit` currently derives
  `(VT, IS)` by reading the `exp` argument. For anything else there is
  no principled scale, and inventing one would be the kind of unvalidated
  heuristic §7 rejected. Options, in order of honesty: leave it
  unqualified (today's behaviour); let the *user* declare a limiter in
  the element; derive a local scale `f/f'` and prove it on a benchmark
  before shipping it.

  Proof, for the multi-junction half still outstanding: a two-junction
  element against `elements`' BJT-shaped limiting. Risk: medium-high —
  it changes shared solver code that the hand-written devices also use,
  so the existing PCNR tests are the regression gate. (The charge half
  shipped with `tests/test_pcnr_charge.py`, 5 tests, and the DSL's
  qualification rule dropped its no-charge condition accordingly.)

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

**B1. `$param_given` + parameter ranges + `aliasparam` — DONE
2026-08-22.** One design point was not obvious in advance: `$param_given`
had to become a **runtime** value, not a compile-time one. The element is
compiled once per *class* while givenness is a property of each
*instance*, so the flag is bound per instance and used in a `Piecewise`
exactly as Verilog-A uses it in an `if`. Note what the operator actually
asks: not "does this differ from the default" but "did the user write
it" — a parameter given its own default value is still *given*, and
`test_param_given_selects_a_formulation` pins that distinction.
`ParameterDict` now records which names were set (`is_given`), validates
declared `minval`/`maxval` on both assignment and after string-expression
resolution, and `update_values` deliberately does **not** mark the names
it copies, since a resolved `iparv` has a value for everything.
`aliasparams = {alias: canonical}` on the element class completes it.
Original scoping: By call count
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

**B2. Unconditional node collapse — DONE 2026-08-22.** `V(a,b) <+ 0`
now merges the nodes and deletes both the internal node and the branch
unknown, so the optional-series-resistance idiom costs nothing when the
resistance is absent (measured: `n` = 2 instead of 3, and the stamp is
the plain junction). Two restrictions, both deliberate: chains are
chased to a fixed point, and **only an internal node can be absorbed** —
terminals belong to the parent circuit's node map, so a `V <+ 0` between
two terminals stays an honest zero-volt source rather than being
silently mis-collapsed. Original scoping: `V(a,b) <+ 0` should merge the two
nodes and delete the unknown, which is the standard idiom for an
optional series resistance (`if (rs) I(a,ia) <+ ...; else V(a,ia) <+ 0;`).
Adopt gnucap's restriction to *unconditionally executed* contributions
(`mg_in_analog.cc:2065-2073`): a conditional collapse changes the
sparsity pattern per iteration and is Phase D. Implementation is a
substitution before the unknown vector is built — merge the node symbols,
drop the row. Proof: an element with `rs=0` producing exactly the stamp
of the element without the branch, and node count reduced by one. Risk:
low, and it is pure DSL work.

**B3. AC excitation — DONE 2026-08-22.** `ac_stim(mag, phase)` routes to
an AC-only source vector: live in AC analysis, identically zero in DC and
transient, and the bias-leak guard of §2.4 is unchanged (only `ac_stim`
terms can drive a small-signal analysis). Original scoping: `u` is currently zeroed for `analysis='ac'` so a
device's bias constants cannot leak into the small-signal drive (§2.4).
The missing half is a *deliberate* AC source: an `ac`-variant vector, so
a behavioural element can be an AC stimulus rather than merely
transparent to one. Follows `VS.u`/`IS.u`'s analysis switch. Proof: an
HDL source driving an AC analysis to the analytic transfer function.
Risk: low.

### Phase C — analyses the DSL cannot currently reach

**C1. Events: `@cross` — DONE 2026-08-22** (`@timer` not done; it is a
source-side concept with no call sites in either corpus). `Cross(expr,
direction)` is returned alongside contributions; the element caches two
accepted points, and `next_event` extrapolates the zero crossing
linearly, honouring direction and the strictly-future contract.
Measured on a comparator across a 1 kHz sine: without it the nearest
timepoint to an edge is 5.5e-6–1.4e-5 s away, with it 5.4e-7–1.9e-6 s —
about an order of magnitude closer, for **one** extra accepted step and
no extra rejections. Original scoping: Zero call sites in either corpus, so
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

**C2. Symbolic toolkit — DONE 2026-08-22.** Elements now keep their
symbolic `i`/`q`/`G`/`C` alongside the compiled ones and substitute into
them exactly when the toolkit is symbolic. This mattered more than
"plumbing" suggests: the numpy lambda survived symbols by duck typing for
plain arithmetic and failed the moment an expression contained a `floor`
or a `Piecewise` — the wrap of an `idtmod`, the branch of a `limexp` —
which is precisely the interesting case. Original scoping: Compilation targets numpy/jax, so a
symbolic transient of a generated element is out (AC works, since it
uses `G`/`C` only). The fix is to keep the sympy expressions rather than
lambdifying when the toolkit is symbolic — the expressions already exist,
so this is plumbing, not derivation. Proof: the existing symbolic tests
extended to a generated element. Risk: low. Value: mostly for the
distortion/DDD analyses, which is where this repo's symbolic work lives.

### Phase D — revisited 2026-08-22

Phase D was the *deferred* table: items with a stated trigger rather than
a date. Revisiting it deliberately, three turned out to be worth doing
now and the rest keep their deferral with the reasoning sharpened.

**Done:**

* **`$discontinuity` — parsed and ignored.** It is advisory: it lets a
  simulator drop order rather than discover a corner by rejection, and
  this one does that from breakpoints (`Cross`) instead. gnucap accepts
  it and generates an empty body for the same reason. Accepting it costs
  nothing and lets a model written elsewhere compile unchanged — 41 call
  sites — where rejecting it would fail those models over a hint they
  can do without.
* **`laplace_nd` and `laplace_zp`.** The deferral said "0 call sites;
  gnucap ships only 2 of 4 `zi_*` with two open bugs; an explicit
  state-space element is probably the better answer even then". The last
  clause is what changed: a rational transfer function **is** a
  state-space element here, and the compiler already emits states for
  `idt`, so the realisation is a generalisation of machinery that
  existed rather than new solver work. An order-*N* denominator becomes
  *N* unknowns in controllable canonical form, integrated by the
  simulator's own method with an exact Jacobian through them — which a
  convolution-based implementation cannot offer. Verified against
  |H(jω)| across three decades for first and second order, and
  `laplace_zp` expands to the same filter. Improper transfer functions
  are refused (that is `ddt`'s job).
* **Switch branches — now *refused*, which was the point.** The plan said
  "diagnose and reject clearly rather than half-implement", and the
  compiler was doing neither: V and I contributed to one branch were
  silently accepted and built a voltage source with a conductance in
  parallel — defined, but not what the model says. Named now.

**A sympy trap found while doing it**, worth recording because it is
invisible: `sympy.Float(0.0) == 0` is **False**. Sympy's `==` is
structural, not numeric, so a real pole written as `(-1/tau, 0.0)` took
the complex-conjugate branch and produced a filter of *twice* the
intended order — no error, just the wrong answer. Use `.is_zero`.

**Still deferred, with the trigger:**

| capability | why | what would reopen it |
|---|---|---|
| `zi_*` (Z-domain) | 0 call sites; needs a sample clock, i.e. an event-driven state the analog solver does not have | a mixed-signal model with a real sample rate |
| `absdelay` | 4 sites, all transmission lines, and `elements.TLine` already covers that | a delay model that is not a T-line |
| `transition`, `slew` | 0 real call sites; gnucap's own `slew` is a self-described stub; `VPulse` exists. **Note the deferral is now weaker than it was**: `Cross` plus `accept_step` give exactly the event-and-history machinery these need, so the cost has dropped from "new infrastructure" to "a few days of care" | a digital-to-analog boundary model, or an op-amp macromodel wanting a real slew rate |
| `last_crossing` | 0 sites, and gnucap does not implement it either. `Cross` now makes it cheap — it is the same prediction, returned as a value instead of a breakpoint | demand |
| `noise_table` | 0 sites; unimplemented in gnucap | demand |
| conditional node collapse | changes the sparsity pattern per iteration; needs gnucap's whole per-iteration re-stamping machine | a model that genuinely needs it — until then the switch-branch refusal names it |
| `@timer` | 0 sites; a source-side concept, and `VPulse`/`VSin` cover the waveforms | a model that schedules its own events |

### Phase E — the PSP target

Everything above was aimed at *behavioural* models. The target that
decides whether this DSL is a production compact-model language is
**PSP103** — the industry-standard surface-potential MOSFET, whose full
Verilog-A source (7400 lines across 13 files, version 103.8.2) and
359-parameter model cards ship in the IHP Open PDK at
`~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/verilog-a/psp103/`. That PDK
is the test basis: real extracted parameters, three simulator flavours,
and no binning (PSP's own geometry scaling does the job), so one model
card exercises the whole geometry range.

The research ranked three blockers. Two are now cleared:

**E1. Expression-graph blowup — DONE 2026-08-22.** The `#1 practical
risk, ahead of any language feature`. Eager substitution is exponential
in intermediate reuse; `var()` plus a forward-accumulated let-chain is
linear. §3.2b has the theory and the measured tables;
`benchmarks/hdl_chain_scaling.py` and
`pycircuit/circuit/tests/test_hdl_chain.py` are the evidence. 1400
intermediates across 4 terminals compile in 21 s with the Jacobian
matching finite differences to 1e-19.

**E2. Both-arms-safe conditionals — DONE 2026-08-22.** PSP's ~70
bias-dependent guards, whose untaken arms overflow, divide by zero or
take roots of negatives. Measured first: numpy's `select` discards the
bad *value* correctly on both compile paths, so there is no wrong answer
— but it raises a floating-point flag on every evaluation, and the
clamped-`sqrt` *derivative* is genuinely `0/0`. §3.2c has the kernel
(`expl` family, `hypsmooth`, `safe_*`), the three traps that broke
earlier versions, and the `limexp`/PCNR trade;
`pycircuit/circuit/tests/test_hdl_kernel.py` pins all of it.

**E3. Language surface.** Found by reading the shipped PSP source rather
than the LRM, which is why two of the three turned out different from the
original estimate:

**Named branches — DONE 2026-08-22.** `Branch(a, b, 'ch')` is
Verilog-A's `branch (a,b) ch;`. The compiler keyed branches on the node
pair alone, so `branch (a,b) br1, br2;` merged into one unknown and the
second constitutive relation vanished. Worse, `Quantity`'s sympy
identity dropped the name too, so `I(br1)` and `I(br2)` were the *same
atom* and both resolved to whichever current was substituted last — a
wrong stamp rather than a missing one. Both fixed; an unnamed branch
keeps its old identity, so nothing that worked before changes.

**`$simparam` — DONE 2026-08-22.** Resolved at compile time, since every
parameter pycircuit can answer is fixed for the run. `gmin` answers
`0.0`, which is the truth and not a placeholder: pycircuit's gmin is a
*continuation schedule* inside the DC solver, ramped away before the
answer is returned, not a standing conductance models should shunt
themselves with — and `0.0` is exactly what PSP passes as its default.
An unknown name with no default raises rather than returning zero.

**`$mfactor` — not needed, measured.** PSP 103.8.2 references it at
exactly two lines (`PSP103_module.include:2183-2184`), both inside the
noise *operating-point output* block, and its own Changelog line 3 reads
"Remove $mfactor from OP output variables calculation". Multiplicity in
PSP flows through its own `MULT_i` model parameter, which is an ordinary
parameter the DSL already handles. Auto-scaling every contribution by a
multiplicity would change the meaning of every existing element to serve
a construct the model does not use.

What is left:

**Optional series resistances — the workaround measured, the feature
sized.** PSP's seven parasitics are each

```
if ((R) > 0.0) I(N1,N2) <+ (G) * V(N1,N2);
else           V(N1,N2) <+ 0.0;
```

Read as topology that is a collapse decided by a parameter *value*, and
the element's unknown count depends on its parameters. Read as algebra
both arms are one relation, `V(N1,N2) <+ I(N1,N2)*R`, which the DSL
already expresses: the MNA row is `-(v1 - v2) + i_br*R = 0`, the resistor
for `R > 0` and an exact short for `R = 0`. Measured bit-exact at every
resistance from 0 to 1e9 (`test_hdl_collapse.py`), where the conductance
form divides by zero — which is why the model branches at all.

It is a workaround, not the answer, and `benchmarks/collapse_cost.py`
says how much of one. A chain of 100 devices each carrying the seven
parasitics:

| devices | unknowns collapsed | unknowns as branches | DC solve | slowdown |
|---|---|---|---|---|
| 5 | 8 | 78 | 0.57 → 0.95 ms | 1.7× |
| 20 | 23 | 303 | 1.11 → 5.80 ms | 5.2× |
| 50 | 53 | 753 | 3.61 → 37.9 ms | 10.5× |
| 100 | 103 | 1503 | 7.02 → 100.6 ms | **14.3×** |

Identical answers; purely a matrix-size argument, and a decisive one.
In the IHP card most of the seven *are* zero for the ordinary device
(`rsh`, `rshd`, `rvpoly` are 0; `rgo` and `rbulko` are gated on
`rfmode`), so this is the common case rather than a corner.

**A collapse bug found on the way, which would have hit PSP directly.**
The collapse test was `sympify(rhs) == 0`, and sympy's `==` is
*structural*: `Float(0.0) == 0` is False. So `V(a,b) <+ 0` collapsed
while `V(a,b) <+ 0.0` silently did not — and `0.0` is exactly how PSP's
own macro spells it. Now `.is_zero`. A second bug in the same code
dropped any branch merely *touching* a collapsed node, deleting the
constitutive relation of whatever hung off the absorbed node; only the
branch whose two ends meet is gone now. A third dropped branch *names*
during the rewrite, merging two named branches that a collapse brought
onto one node pair.

**Parameter-driven collapse — DONE 2026-08-22.** `Collapse(branch, when)`
alongside the resistor contribution:

```python
br = Branch(n1, n2, 'rd')
return (Contribution(br.I, br.V / rd),
        Collapse(br, rd <= 0))
```

When the condition holds for an instance's parameters, the branch's nodes
merge and every contribution to it is dropped — so the `V/rd` that would
be infinite is never compiled, which is the whole reason the model
branches. The condition may mention **parameters only**; one that
depended on the operating point would move the sparsity every Newton
iteration, and is refused by name.

Implemented by compiling a **subclass per mask** and retargeting the
instance's `__class__`, rather than by teaching every generated method to
look up a per-instance function table. The metaclass already compiles
everything a class needs from its `analog`, so re-entering it with the
mask baked in gets correct stamps, branches, PCNR participation, state
metadata and pure forms with no second code path to keep in step.
Variants are built on demand and cached on the base class, so N instances
sharing a parameter set share one compilation. A gating parameter changed
*after* construction raises: it would change the element's size behind a
node map the circuit has already fixed.

Measured on the same benchmark, the collapsed variant tracks the
hand-written topology exactly:

| devices | unknowns hand / branch / `Collapse` | DC solve hand / branch / `Collapse` |
|---|---|---|
| 20 | 23 / 303 / 23 | 1.13 / 5.87 / 1.72 ms |
| 100 | 103 / 1503 / 103 | 6.99 / 99.3 / 8.10 ms |

**First production compact model — DONE 2026-08-22.** Before committing
to PSP's 7417 lines, the whole pipeline is now proved end to end on a
smaller real model. The PDK ships a ladder of them, each with its own
compiled binary and model cards: `cap_cmomf` (166 lines), `cap_cmomi`
(314), `r3_cmc` (1645), `mosvar` (1380), `psp103` (7417).

`pycircuit/circuit/compact.py::CapCmomi` is the 314-line interdigitated
MoM capacitor, translated from the shipped Verilog-A. It exercises,
together and for the first time, essentially everything Phase E added:
eight **named branches** (three pairs sharing a node pair — the skin arm
is `L || R`, the substrate termination is `R || C`), `ddt` of a branch
**current** as well as of a voltage, parameter-only conditionals
selecting fitted coefficient sets by layer count / feed style /
bottom-metal index, `floor` on a parameter expression, and clamps.

Checked against IHP's own OSDI build of the same model, swept in ngspice
over 1 MHz – 100 GHz for ten configurations reaching every coefficient
branch:

**worst relative impedance error 2.0×10⁻⁸**, over all cases and all
frequencies — at the reference's own printed precision (ngspice `wrdata`
writes 9 significant digits).

The low-frequency capacitance is *additionally* checked against the
vendor's tiling arithmetic done by hand, so a wrong reference cannot make
a wrong model pass. Note what the model costs: the eight named branches
need only **two** branch-current unknowns, because only the two inductors
are voltage-defined.

Notable: `Contribution(b.V, ddt(L * b.I))` — an inductor as a
voltage-defined branch whose own current is inside the `ddt` — already
worked (`LHdl`), and the two inductors here are the first use of it in a
production model.

**The surface-potential kernel — DONE 2026-08-22.**
`pycircuit/circuit/psp_kernel.py` translates PSP103's `sp_s`, `sigma`,
`sigma2` and `MINA` — the *explicit* approximation to the root of the
surface-potential equation, which is the heart of the model and the piece
everything else rests on. PSP does not iterate for it: it computes a
starting point and applies two closed-form corrections, so the whole
thing is straight-line code the simulator can differentiate.

It is checked **without a vendor binary or a model card**, two
independent ways:

* the surface-potential equation is evaluated at the returned root. That
  residual is not quoted from a paper — it is the expression PSP's own
  final correction computes as its `qC` term, so the test and the code
  come from the same source. Worst over the physical envelope
  (`Gf` 0.3–5, `xn` 15–45, `xg` −100…300): **≤ 6.7×10⁻⁸, typically
  2×10⁻⁹**.
* the DSL's forward-accumulated `d(sp)/d(xg)` is compared against the
  **exact implicit derivative** `−F_xg/F_x` of that same equation —
  agreeing to ~10⁻⁷. Translation and Jacobian, checked by different
  routes.

Two limits are pinned rather than avoided: the approximation degrades to
~10⁻⁵ as `Gf → 0` (the equation becomes degenerate — the root turns into
a double root) and at `xn = 0` (no inversion barrier). Both are far
outside any real device.

**Note the PDK ships two versions of `sp_s`.** The PSP-based `mosvar`
varactor carries a simplified one without the ξ₀/ξ₁/ξ₂ terms, so it
solves a slightly different equation. This implements PSP103's. Internal
consistency was checked rather than assumed: the macro's `xi1` is exactly
`dξ₀/dx = 4x/(2+x²)²`, which is what it must be for `pC` to be `−F′(x)`
and the final step to be the Halley correction it is written as.

**Three findings from making it safe**, each of which cost a debugging
session:

1. **`safe_div`'s `eps` is squared, so it has a floor.** The kernel first
   used `eps = 1e-300`; `1e-300²` underflows to *exactly zero*, the
   regularisation silently vanished, and flat band gave `0/0`. `safe_div`
   now refuses an `eps` whose square underflows.
2. **Both-arms-safe has to be applied at the REGIME level, not
   operation by operation.** With every individual operation guarded and
   every intermediate finite, the inversion arm at `xg = −77` still
   produced `ei ~ 10¹¹⁴` — finite, but large enough that products in
   `sigma2`'s *derivative* overflowed to `inf`, and `inf − inf = NaN`.
   Chasing operations does not converge. Clamping each arm's **input to
   its own validity domain** fixes the whole class at once, and is exact
   where the arm is selected.
3. **The chain rule leaks NaN out of unselected arms.** A `Piecewise`'s
   partial derivative w.r.t. an intermediate used only in a discarded arm
   is *zero* — and `0 × NaN` is `NaN`. So it is not enough for the
   selected arm to be finite; every arm must be. This is the concrete
   form of the `0 * inf` hazard the original research flagged.

A fourth, smaller: `sp_s` **cannot** be built as one flattened
expression — it does not finish in ten minutes, where the let-chain takes
0.13 s. `hdl.compile_chain` was added so a kernel routine can be compiled
and tested outside a model without that fallback existing at all.

**A surface-potential MOSFET — DONE 2026-08-22.**
`psp_kernel.ids_long_channel` assembles PSP's intrinsic drain current
from the surface potentials at the two channel ends, and
`compact.PspMosLongChannel` wraps it as an element. PSP builds the
current as `Ids = BET*(FdL*qim1*dps*Gvsatinv)`
(`PSP103_module.include:1178`); with channel-length modulation and
velocity saturation at unity — the long-channel intrinsic device — what
remains is `Ids = β·qim1·Δψ`, the symmetric-linearisation charge-sheet
current.

There is **no vendor reference** for this one, deliberately: it is PSP's
core without the short-channel, mobility, series-resistance and geometry
layers, so it is not the same device as any model card and comparing I–V
curves would compare two different things. What is checked instead is
what the formulation guarantees *by construction* — a stronger statement
than a curve fit, because a construction error breaks these exactly
while a fitting error does not:

* **source and drain are exactly interchangeable** — swapping them
  negates the current bit for bit. Threshold-voltage models are famously
  asymmetric about `Vds = 0`, which shows up as spurious distortion in
  passing gates and mixers; here symmetry falls out of evaluating the
  inversion charge at the midpoint potential;
* **exactly zero current at zero drain bias**, at any gate bias, with no
  epsilon;
* **one expression spans accumulation to strong inversion** — no
  regional equations, no smoothing functions pasted between them, and
  monotone over every decade.

The body effect appears without being put in by hand: a larger `Gf`
degrades the subthreshold slope, in the right order.

**The charge model** (`psp_kernel.charges_long_channel`) comes with the
same character. `Qg`, `Qb` and `Qd` are contributed on branches referred
to the source, so `Qs` is whatever is left and **charge conservation is
structural** rather than numerical — measured exact to rounding, and the
capacitance matrix's rows and columns both sum to zero, so no
common-mode shift of all four terminals injects current.

The drain/source split is PSP's `SWQPART = 0` linear distribution, and it
reproduces the textbook Ward–Dutton result without being fitted to it:

| Vds | Qd/(Qd+Qs) |
|---|---|
| 0 V | **0.5000** (exactly — symmetric channel) |
| 0.5 V | 0.4304 |
| ≥ 0.8 V | **0.395** (the 40/60 partition) |

Swapping source and drain swaps `Qd` and `Qs` to **six ulp** — not
bit-exact like the current, because the two are algebraically mirror
images but not literally the same expression; the gate charge, symmetric
in `dps`, does come back bit-identical. And below threshold the gate
charge is terminated by the bulk to within 3×10⁻¹², the residue being
the subthreshold inversion charge, which had better not be zero since it
*is* subthreshold conduction.

Two measured limits are pinned rather than smoothed away. The current is
strictly monotone above **1e-18 A**; below that it is a difference of two
surface potentials agreeing to nine digits, and the cancellation leaves a
few percent of wobble — one attoamp is a thousandth of any real
subthreshold measurement floor, and the only cure would be more precision
than a double has. And `psp_kernel.XN_FLOOR` bounds the normalised
quasi-Fermi splitting at zero: a junction forward-biased past its
built-in potential would drive `delta = exp(-xn)` to 1e152 and overflow,
so Newton gets a bounded wrong answer there instead of a NaN, exactly as
PSP clamps its own junction voltages.

**Mobility reduction and velocity saturation — DONE 2026-08-22.** The
first layer on top of the ideal core, and the biggest determinant of
strong-inversion accuracy. Defaults are the IHP SG13G2 n-channel card's
values (`mue = 0.779`, `themu = 2.05`, `cs = 0.316`, `thecs = 1.18`,
`thesat = 0.398`); zeroing `mue`, `cs` and `thesat` recovers the ideal
core **bit-exactly**, so the layer is genuinely additive.

The mobility ratio against the ideal device falls monotonically with
gate bias — 0.73 at Vg = 0.4 V to 0.55 at 1.8 V — which is the whole
physical content of the term. And the construction survives being built
on: both factors depend only on midpoint quantities and on `dps`
*squared*, so both are **even** under the source/drain exchange, the lone
odd factor stays `dps`, and the antisymmetry of the current is still
`==` rather than `approx`.

> **A real bug here, worth recording.** PSP writes Coulomb scattering as
> `cs*exp(0.5*thecs*ln(Pm/(Pm+Dm)))`, because Verilog-A's `pow` with a
> variable exponent is awkward. Compiled literally that nests two
> `Piecewise`s (`expl` inside `safe_ln`'s `hypsmooth`), and
> differentiating the nest emitted conditions built from
> `logical_and.reduce` over scalars that numpy's `select` **rejects
> outright** — the Jacobian raised rather than losing precision.
> Algebraically the whole thing is just `cs * ratio**(0.5*thecs)`, and
> the power form is the one to compile. The ratio is floored because
> `thecs = 1.18` puts the exponent below 1, so the derivative diverges as
> the ratio goes to zero — which it does exactly at flat band.

### Compiled-code speed

Profiling the new device's Jacobian exposed two costs in the chain
printer that every chained model was paying:

* sympy prints `Piecewise` as `numpy.select([conds], [values])`, which
  is built for arrays and spends its time broadcasting and allocating
  for the SCALAR arguments a device model is called with — **950 calls
  and 60% of the runtime** per Jacobian. `numpy.where` has identical
  both-arms-evaluated semantics, which the kernel's safety work depends
  on, and none of the machinery;
* `Min`/`Max` print as `functools.reduce(numpy.minimum, [...])`, which
  allocates a list and enters `reduce` for what is nearly always a
  two-argument call — another 18%.

Both are fixed in `_ChainPrinter`. Jacobian evaluation on the MOSFET went
from **11.6 ms to 6.3 ms**, and the chain's own arithmetic now dominates,
which is where it should be.

**The gap, measured — 2026-08-22.** Everything above is validated by
construction properties and internal consistency. Those say nothing
about how close the device is to the transistor a foundry ships, so
`pycircuit/circuit/psp_scaling.py` maps a real card through PSP's
geometry scaling into the element, and `benchmarks/psp_gap.py` compares
the result against the committed PSP103 reference. Nothing is tuned: the
parameters come from the card.

| sweep | W | L | first | +QM | +Rs | +`CT` | +CLM | +GPE | +poly | +`FdL` |
|---|---|---|---|---|---|---|---|---|---|---|
| `nmos_long_idvd` | 10 µm | 1 µm | 1.29 | 1.095 | 1.041 | 1.008 | 1.010 | 1.009 | **0.998** | 1.035 |
| `nmos_long_idvg` | 10 µm | 1 µm | — | — | — | 1.10 | 1.10 | 1.10 | — | **1.00** |
| `nmos_idvd_vg1p2` | 1 µm | 0.13 µm | — | 1.221 | 1.115 | 0.985 | 1.021 | 0.959 | 0.950 | **0.996** |
| `nmos_idvd_vg0p6` | 1 µm | 0.13 µm | — | — | — | — | — | 0.735 | 0.732 | **1.111** |
| `nmos_idvg_vd1p2` | 1 µm | 0.13 µm | — | — | 1.028 | 0.868 | 0.974 | 0.898 | 0.892 | **1.032** |
| `nmos_idvg_vd0p05` | 1 µm | 0.13 µm | — | 1.703 | 1.288 | 1.115 | 1.117 | 1.019 | **1.010** | 1.011 |
| `nmos_idvg_vb_m1` | 1 µm | 0.13 µm | — | — | — | — | — | 1.023 | **1.011** | 1.014 |

(median ratio, ours / PSP103.) Every sweep is now within a few percent,
from a card, with no fitting anywhere. Which geometry is *most* accurate
has inverted twice along the way, which is worth more than any single
number: the long device led until `GPE`, the short device leads now.

Each column is a term the *shape* of the previous residual named:

* the first pass was 23–33% high with a ratio that **grew with current**
  → the quantum-mechanical correction (below);
* what remained was a nearly **flat** ~10%, the signature of something
  multiplicative → **series resistance**, which PSP folds into the
  mobility as an extra `Gmob` term (`THER = 2·BET·RS`) rather than adding
  a network element. That keeps the device four-terminal *and* symmetric,
  since `qim` is a midpoint quantity. Worth 1.095 → 1.041;
* what then remained was a flat offset again, and the diagnostic for it
  had to be *built*: a **long-device transfer curve** was added to the
  reference, because a gain error gives a ratio flat in Vg while a drive
  error gives one that varies. The measured ratio fell from 1.175 near
  threshold to 1.022 at Vg = 1.5 — a drive error → **the effective
  thermal voltage**. PSP does not normalise by `phit` but by
  `phit1 = phit·(1 + CT)`, and the gate drive, the quasi-Fermi levels and
  the charges are all in units of *that*
  (`PSP103_macrodefs.include:503`). The card sets `CTO = 0.0546`, worth
  ~7% here. 1.041 → **1.008**, and the short device improved with it
  (1.288 → 1.115). The body factor is deliberately not rescaled: PSP
  builds `G_0` on the plain `phit` and only moves it under `SWFIX`, which
  defaults to 0;
* **channel-length modulation, on the second attempt.** With the flat
  offset gone the residual was centred rather than one-sided, and the
  same CLM that had made every sweep worse now helps in four of seven —
  dramatically on the short device (median |ratio−1| 0.098 → 0.033 on
  `nmos_idvd_vg1p2`, 0.132 → 0.026 on `nmos_idvg_vd1p2`). The long
  device barely moves, which is right: CLM is a short-channel effect;
* what remains sits on the short device at low drain bias — and it is
  **not** what it looked like. See below.

> **The threshold block: measured, and mostly already there —
> 2026-08-23.** The residual pointed at the short-channel threshold, so
> the obvious next layer was DIBL. Measuring first said otherwise, twice
> over.
>
> `CF`, PSP's DIBL coefficient, scales to 1.1 × 10⁻⁷ on the long device
> and contributes **under a millivolt** at Vd = 0.05 — which is exactly
> where two of the three worst sweeps sit. DIBL cannot be the
> explanation for a 12% error at a bias where it does nothing.
>
> So a better instrument was built. ngspice will print an OSDI model's
> operating-point outputs, and PSP exposes its own **`vth`** — which
> turns "is our scaling right?" from an inference about currents into a
> direct comparison. A current can be right for compensating reasons; a
> threshold cannot. The reference now records it at four geometries.
>
> | geometry | PSP shift vs long | our shift | error |
> |---|---|---|---|
> | L = 0.5 µm | +69.2 mV | +65.5 mV | −3.7 mV |
> | L = 0.13 µm | +237.3 mV | +214.3 mV | −22.9 mV |
> | L = 0.13 µm, W = 10 µm | +198.9 mV | +179.4 mV | −19.5 mV |
>
> The *absolute* values differ by a nearly constant 70–93 mV across a
> 7.7× range of channel length, which makes that a definitional
> difference — PSP's `vth` is its own extraction, not
> `vfb + phib + γ√phib`. The *shift* is the comparable quantity, and we
> already have **90% of the reverse-short-channel effect**, out of the
> pocket-implant doping formula plus the length terms in `VFB` and
> `DPHIB`.
>
> The missing 23 mV is worth about 3% of drain current at Vg = 1.2 —
> **a quarter** of the ~12% the short device is off by. So the threshold
> block is largely present, DIBL is not the next layer, and what remains
> on the short device is predominantly a **gain** error.

> **The gain error, found term by term — 2026-08-23.** The `vth` trick
> generalises: PSP exposes its whole post-scaling parameter set as
> `lp_*` operating-point outputs, so the scaling layer can be checked
> **term by term** instead of inferred from currents. The reference now
> records fourteen of them at four geometries.
>
> Every term matched to five digits — `vfb`, `tox`, `ct`, `mue`,
> `themu`, `cs`, `thecs`, `rs`, and `neff` with the quantum correction
> off — **except `betn`, which was 12% high**. Back-solving gave PSP's
> `GPE` as 1.4234 where ours was 1.271, and the reason is in the source:
> `FBET1e = FBET1·(1 + FBET1W·iWE)` and `LP1e = LP1·max(1 + LP1W·iWE,
> 1e-3)` (`PSP103_scaling.include:284-285`). The trailing `e` is not
> decoration — both are **width-adjusted** before use, and the raw card
> values were being used instead.
>
> Since `BETN = UO·WE·GWE/(GPE·LE)`, that was a 12% gain error, and it
> is most of what the short device was off by:
>
> | sweep | before | after |
> |---|---|---|
> | `nmos_idvg_vd0p05` | 1.117 | **1.019** |
> | `nmos_idvg_vb_m1` | 1.128 | **1.023** |
> | `nmos_long_idvd` | 1.010 | 1.009 |
>
> The long device barely moves, which is the check: `GPE` is 1.02 there,
> so the bug was invisible on the geometry that had been driving the
> whole investigation. No amount of staring at I-V curves would have
> said *which* term; one comparison did.
>
> The scaling is now guarded term by term, which is the strongest form
> this check can take.
>
> Two things the comparison also surfaced: `lp_np` is 4.6 × 10²⁶, so
> **polysilicon depletion is active** and the core's ideal-gate
> assumption (`eta_p = 1`) is a real omission; and the short device is
> now *low* at high drain bias (0.90–0.96), where before it was high.
> DIBL raises the drive, so its direction is now right — the term ruled
> out earlier is worth revisiting, for the same reason CLM was.

> **The channel-shortening factor `FdL` — 2026-08-23.** The residual left
> the short device 27% low at Vg = 0.6, with the reference current
> climbing **2.4× through saturation** where ours was flat. DIBL was the
> obvious suspect and the measurement killed it: PSP's own `vth` moves
> **3.5 mV over 1.35 V of drain bias** — 2.6 mV/V. (Stated more exactly
> than it was: we do not model DIBL at all — `CF`, `CFB` and `CFD` are
> absent from the scaling layer — so the measurement is what licenses
> the omission, not agreement with a term we have. PSP's own shift is
> 3.2 mV at Vds = 1.2 V here, worth about 1% of the current at a
> strong-inversion `gm/Id`.) So the climb happens at essentially constant
> threshold
> and is in the current formula.
>
> It is `FdL = (1 + dL1 + dL1²)·GdL`
> (`PSP103_module.include:1137`), which I had left out on the reading
> that `ALP1` and `ALP2` are zero — they would make `dL1 = dL` and `FdL`
> exactly 1. **They are not zero**: PSP's own `lp_alp2` is **4.5** on a
> 0.13 µm device. That is not visible in the card, because `ALP2` has no
> geometry-independent coefficient and its value comes out of `ALP2L1`
> through a length power; only the term-by-term comparison showed it.
>
> `ALP2` multiplies the *bulk* charge, so it acts near threshold —
> exactly where the deficit was. The worst sweep goes **0.732 → 1.111**,
> and the summed median error across all seven sweeps **more than
> halves**, 0.480 → 0.216.
>
> It is mixed by *count* — better in three sweeps, worse in four — and it
> costs the long device (0.998 → 1.035). That is the price of our `dL`
> being an approximation to PSP's: `FdL` amplifies whatever error `dL`
> carries. Which makes the `Vdsat`/`Vdse` machinery, skipped twice
> already, the next thing genuinely worth building.

> **Polysilicon depletion, and a symmetry bug it exposed —
> 2026-08-23.** The gate is not a perfect conductor: a depletion layer
> forms on its silicon side and takes a share of the applied voltage.
> PSP folds this in as `eta_p = 1/sqrt(1 + kp·xgm)` on the charge slope
> *and* as a correction that shifts the midpoint potential itself
> (`PSP103_macrodefs.include:697-724`). Implemented whole rather than as
> the `eta_p` factor alone — half of a coupled correction is not a
> smaller version of it, which is the CLM lesson applied in advance.
>
> Worth about 1% here, and it took the **long device to 0.998**
> (0.979–1.025) with `nmos_idvg_vd0p05` to 1.010 and `nmos_idvg_vb_m1`
> to 1.011.
>
> **It also broke source/drain symmetry — and revealed that two earlier
> layers already had.** The intrinsic core is exactly antisymmetric
> because every quantity it uses is either a midpoint one (even) or
> `dps` (odd). The layers on top are not automatically so: series
> resistance reads a source–bulk voltage, which *changes* under the
> exchange, and channel-length modulation reads a signed `Vds` through a
> function that is not odd. Both had been breaking symmetry since they
> went in, **invisibly — because `rs` and `alp` default to zero and
> every symmetry test built its element with defaults.** Only a test
> using real card parameters could see it.
>
> Fixed the way PSP fixes it: `Vsbx = Vsbstar + 0.5·(Vds − Vdsx)`
> (`macrodefs:472`), which evaluates to the *lower* of the two junction
> voltages under either polarity, plus a smoothed `|Vds|`. Both are even,
> so the antisymmetry survives — now 2.8 × 10⁻¹⁶ rather than bit-exact,
> the smoothing costing an ulp. The guard now asserts the card actually
> switches `rs`, `alp` and `kp` on, so it cannot silently test nothing
> again.

> **CLM was tried, measured, rejected, and then vindicated —
> 2026-08-23.** The rejection stands as a record of the reasoning, and
> the sequel is the point: with `CT` in place the same code helps. Read
> together they are the methodological lesson of this whole exercise.
>
> **The rejection.** The saturation
> tail is CLM's signature, and ruling out the alternative was easy with
> data already committed: impact ionisation flows to the *bulk*, and the
> reference's bulk current is 8×10⁻¹¹ A against a drain rise of
> 8.6×10⁻⁵ A — **0.0%** of it.
>
> But adding CLM made the fit **worse in all six sweeps and better in
> none** (long device median |ratio−1| 0.041 → 0.045; `nmos_idvd_vg0p6`
> 0.158 → 0.185). Two reasons, both worth keeping:
>
> * PSP's `s1` is a *ratio* of logarithms involving `Vdse`, a smoothly
>   saturation-limited drain voltage that it also uses for the drain
>   surface potential (`PSP103_macrodefs.include:628-632`). This core
>   computes that potential at the true drain bias, so `dps` saturates on
>   its own and the second logarithm was dropped — which recovers the
>   classic `ln(1 + (Vds − Vdsat)/VP)` form but captured only about a
>   third of the measured rise. A faithful CLM needs PSP's `Vdsat`
>   machinery (`Phi_0`, `Phi_2`, `asat`), which is a block of its own.
> * More importantly: the device is already ~4% HIGH, and CLM only
>   pushes it higher. **A term that is individually correct can worsen
>   the fit when a compensating term is missing.** So the flat +4% is
>   not missing CLM, and chasing CLM before finding it is the wrong
>   order.
>
> Ruled out for the +4% by measurement: `Rxcor` (unity at Vsb = 0, which
> every sweep uses), poly depletion (`NP` absent from the card), the edge
> transistor (`SWEDGE = 0`), gate current (1e-6 of the drain current),
> and temperature (the card's `TR = 27 °C` against the 300 K used, worth
> 0.05% on the thermal voltage).
>
> **The sequel.** The +4% turned out to be the effective thermal voltage
> (above). With that in, CLM was re-applied unchanged and now helps in
> four sweeps of seven. Nothing about the CLM code changed between the
> two measurements — only what it was sitting on top of.
>
> The lesson, and it is worth carrying: **a residual's shape tells you
> which term is missing, but not which to add first.** When two terms are
> missing and they pull opposite ways, adding either alone looks like a
> regression. The order that works is to fix what is flat before what is
> shaped, because a flat error contaminates every diagnosis made in its
> presence.

> **The quantum-mechanical correction, found by measuring.** The first
> comparison ran **23–33% high** with a ratio that grew with current.
> The card sets `QMC = 1`: carriers in the inversion layer occupy
> quantised states below the surface, so the surface potential at
> threshold is higher and the body factor larger
> (`PSP103_macrodefs.include:322-327`). It was the only term of that size
> the core was missing, and adding it took the long device from 1.23–1.33
> to **1.05–1.12**. A test pins that switching it off reopens the gap, so
> it cannot be removed as cosmetic.

What is left for a running PSP103 card:

**Model-card ingestion — DONE 2026-08-22.**
`pycircuit/utilities/spicecard.py`. A foundry card is not a list of
numbers: PSP103's is 359 parameters, most of them quoted expressions
referencing `.param` multipliers a *corner section* defines,
`.include`d from inside a `.subckt` so the card can also read instance
parameters. None of it resolves without following that whole chain.

```python
deck = spicecard.read('cornerMOSlv.lib', section='mos_tt')
p = deck.model_params('sg13g2_lv_nmos_psp', w=1e-6, l=0.13e-6, ng=1)
```

Measured on the real PDK: **371 parameters, all resolving to floats**,
with the corner multipliers applied (`dphibo = −0.25737 × 0.9915` in
`mos_tt`, and three distinct values across tt/ss/ff) and instance
parameters reaching the card (`dlq` reads `pre_layout`, `cfrw` divides by
`ng`). Every model family in the PDK reads: MOS lv/hv, the MoM
capacitors, and the `r3_cmc` resistors.

Not a netlist parser and not trying to be one — it reads the declarative
part (`.lib` sections, `.include`, `.param`, `.subckt`, `.model`,
`.if`) and skips the rest. Three decisions worth stating:

* **`.LIB` sections are opt-in.** A corner file defines the same names
  differently per section, so reading them all would silently give
  whichever came last. Without a `section=` argument every block is
  skipped.
* **Statistical functions return their nominal value.** `agauss`,
  `gauss`, `aunif`, `unif` describe a distribution; without a draw the
  centre is the only defensible answer, and it is what a nominal corner
  uses. Monte-Carlo is not implemented rather than faked.
* **The expression namespace is closed.** A vendor file is *data*;
  evaluating it over an open namespace would make reading a card
  equivalent to running it. A test pins that `__import__` cannot be
  reached.

Two bugs the tests caught, both of which silently truncate rather than
fail: a comment line *between* `+` continuations ended the logical line
(real cards put commentary there), and circular `.include` detection was
built but never actually tested against.

| item | size | note |
|---|---|---|
| geometry scaling layer | medium | `PSP103_scaling.include`, 849 lines, pure parameter arithmetic — no solver involvement, so it is bulk rather than difficulty |
| the rest of the surface-potential core | large | `PSP103_module.include`, 2371 lines. The `sp_s` kernel it is built on is done and validated; what remains is the current, charge and geometry layers around it |

> **The saturation voltage, a floor hidden in the scaling, and two
> things it broke — 2026-08-23.** PSP does not evaluate the drain
> surface potential at the applied drain bias. It computes a saturation
> voltage `Vdsat` from **source-side quantities alone**, smoothly limits
> `Vds` to it, and uses the limited `Vdse` for the drain quasi-Fermi
> level: `xn_d = xn_s + Vdse/φt` (`PSP103_macrodefs.include:596-632`).
> This was the missing block behind two earlier compromises — the
> dropped second logarithm in CLM's `s1`, and `FdL` amplifying an
> approximate `dL`.
>
> Built whole: `xgs`, `qis`, `qbs` and a source-end `Gmob`, then `Φ₀`,
> `Φ₂`, `asat`, `Φ_sat` and `Vdsat = Φ_sat − φt·ln(1 + …)`. Result
> across the six n-channel sweeps, nothing tuned — five improved, one
> regressed, summed median error **0.207 → 0.138**, and the long device
> went **1.035 → 1.010** with its range tightened to 0.996–1.014.
>
> **First, though, it made things worse — by 14%, on the sweeps where it
> should have changed nothing.** `Vd = 0.05` is far below saturation;
> a limiter has no business acting there. It acted because the exponent
> `AX` governing it is **floored at 2** in the scaling layer
> (`PSP103_scaling.include:743`), and the card gives no hint: on a
> 0.13 µm device `AXO/(1 + AXL·iLE)` evaluates to **0.88**, and an
> exponent below one makes `(1 + (Vds/Vdsat)^AX)^(−1/AX)` soft enough to
> bite everywhere. The floor is now verified against PSP's own `lp_ax`
> at four geometries — exact at all four, and binding at both short
> ones. `lp_ax`, `lp_thesat`, `lp_thesatb`, `lp_thesatg` and `lp_alp`
> were added to the recorded reference for it; the regeneration was
> purely additive, every existing number reproducing bit-identically.
>
> **It broke source/drain symmetry, and the fix is structural.**
> `Vdsat` is built from the source end alone, so it is not an odd
> function of `Vds` — the reverse current came out 16% larger than the
> forward. There is no cleverer formula available, and PSP does not look
> for one: it **orders the terminals**, computes the device forward from
> the lower junction, and applies the sign on the way out. Same here,
> reusing the `vsbx`/`vdsx` pair that already existed:
> `xn_s = (φ_B + Vsbx)/φt`, `xn_d = xn_s + Vdsx/φt`, current
> `sgn·Ids` with `sgn = Vds/Vdsx`, and the channel charge interpolated
> back onto the real terminals. The antisymmetry is now a property of
> the **topology** rather than of the algebra — which is the only form
> of it that survives adding non-odd physics, and it had already failed
> twice in the algebraic form. It costs the bit-exactness: the two
> polarities reach the same number by different roundings, so the
> agreement is 3 × 10⁻¹⁶ rather than the last bit.
>
> **And it exposed a hard clamp that had been fine only by luck.** The
> kernel floors `xn` at zero — the quasi-Fermi level reaching the
> built-in potential. A hard floor is fine for the value and poison for
> the Jacobian: at exactly `Vd = −φ_B` the analytic derivative and a
> finite difference disagreed by **60%**, and below it every conductance
> froze bit-identically, so a solver that overshot had nothing telling
> it how to return. It had gone unnoticed because the floor only ever
> bit the *drain* end, where the device is off and the sensitivity is
> negligible; `Vdsat` reads the *source* end, which multiplies
> everything. Replaced with PSP's own conditioning
> (`macrodefs:330-334, 1104-1105`): take the lower junction, clip it at
> `−0.95·φ_B` through the smooth `MINA`, lift `Vsb` by what the clip
> removed. Sub-millivolt at ordinary bias, so it applies
> unconditionally. Jacobian error at the old kink: 0.6 → 3 × 10⁻⁷.
>
> **A saturating model needs a Newton limiter — that is new
> infrastructure, not a workaround.** Saturation is the point of
> `Vdse`, and its consequence is that `dIds/dVds` falls to 10⁻¹¹ by
> 500 V and 10⁻²⁸ by 10⁷ V. A solver that lands out there is not slow,
> it is stuck: the row goes numerically empty and the matrix is reported
> singular. Two devices in a stack driving their own internal node were
> enough to walk out there and stay — and the *old* model reached the
> same absurd biases and escaped only because its current kept growing,
> so the divergence was always latent. `PspMosLongChannel.limit` now
> plays the part SPICE gives `DEVfetlim`, in the state-free convention
> the tree prefers (return a limited copy, so residual and Jacobian are
> never taken at different points). The source is left where the solver
> put it — it is some other device's drain, and limiting it here would
> have two elements undo each other.
>
> **One general lesson worth carrying, about the compiler rather than
> the physics.** `sympy.Max(f, c)` must be applied to an **atom**, never
> to an expression: it is differentiated by expanding `f` into every
> branch condition, and when `f` contains a `hypsmooth`, sympy
> rationalises it there as `2ε²/(√(z² + 4ε²) − z)` — whose denominator
> cancels to *exactly* zero for any ordinary `z`, because `4ε²` sits far
> below `z`'s last bit. The value was finite; the Jacobian divided by
> zero. Naming the argument with `var()` first keeps it opaque and the
> derivative stays a one-line `where`. The same shape of trap as
> compiling Coulomb scattering as a power rather than nested
> exponentials, and as `rat**ax` needing its base floored strictly above
> zero (sympy writes that derivative as `ax·rat**ax/rat`).
>
> Still absent from the velocity-saturation term: `THESATB` and
> `THESATG`, the body- and gate-bias modulations of `THESAT`
> (`macrodefs:596-607`). Both are nonzero on this card — 0.082 and
> 0.115 — and both are now recorded in the reference so the gap is
> written down rather than remembered.

> **The channel-length modulation, finally translated rather than
> approximated — 2026-08-23.** Immediately after the saturation voltage
> went in, and the largest single accuracy step this core has taken.
>
> PSP writes the pinch-off length as a RATIO of logarithms
> (`PSP103_macrodefs.include:753`):
>
> ```
> s1 = ln((1 + (Vds − dps)/VP) / (1 + (Vdse − dps)/VP))
> ```
>
> The denominator had been dropped, and honestly labelled as "an
> APPROXIMATION to PSP, not a translation of it", because without a
> saturation-limited drain voltage there was nothing to put in it. What
> was left recovered the classic `ln(1 + (Vds − Vdsat)/VP)` with `dps`
> standing in for the saturation voltage. `Vdse` exists now, so this is
> just the real expression.
>
> **Result across the six n-channel sweeps: every one improved, and
> every one is now within 0.4%.**
>
> | sweep | before | after |
> |---|---|---|
> | `nmos_long_idvd` | 1.010 | **1.002** |
> | `nmos_idvd_vg1p2` | 1.019 | **1.002** |
> | `nmos_idvd_vg0p6` | 1.071 | **1.003** |
> | `nmos_idvg_vd0p05` | 1.004 | **1.001** |
> | `nmos_idvg_vd1p2` | 1.028 | **0.996** |
> | `nmos_idvg_vb_m1` | 1.006 | **1.004** |
>
> Summed median error **0.138 → 0.016**, a factor of 8.6, with nothing
> tuned and no parameter touched — every value still comes off the card
> through the scaling layer.
>
> Why the denominator matters more than it looks: both logarithms grow
> together in the linear region, where `Vdse` tracks `Vds`, so their
> ratio stays near 1 and `dL` stays near zero. That is the physical
> statement that there is no pinch-off region to shorten yet. The
> single-logarithm form had no way to say it and leaned on `dps` to say
> it instead, which is why its error showed up as output conductance in
> the wrong places.
>
> **This also closes the CLM story recorded above.** CLM was tried,
> measured, rejected, then vindicated unchanged once `CT` went in — and
> the rejection note said a faithful version "needs PSP's `Vdsat`
> machinery (`Phi_0`, `Phi_2`, `asat`), which is a block of its own".
> That was correct, and the ordering it implied was correct: the block
> had to exist before the term that consumes it could be written
> properly. The whole sequence — reject on measurement, keep the
> reasoning, build the prerequisite, come back — is the methodological
> point of this exercise, and this is where it pays.

> **Two more correct terms, both measured worse, both kept and switched
> off — 2026-08-23.** DIBL (`CF`, `CFB`, `CFD`,
> `PSP103_macrodefs.include:473-476`) and the body- and gate-bias
> modulation of the velocity-saturation parameter (`THESATB`, `THESATG`,
> `macrodefs:596-607`). Both implemented, both faithful, both off by
> default behind `to_long_channel(..., all_terms=True)`.
>
> Summed median error over the six n-channel sweeps: **0.016** without
> them, 0.034 with the `THESAT` modulation alone, **0.041** with both.
>
> What the evidence rules out is more useful than the numbers. **All
> twenty-one scaled parameters the reference now records match PSP's own
> `lp_*` outputs exactly at four geometries** — including `cf`, `cfb`,
> `thesatb`, `thesatg`, `alp`, `alp1` and `alp2`. So the scaling layer is
> not the problem, and neither are `FdL`'s coefficients. `s1`, `s2`,
> `r1`, `r2` and `FdL` were re-read against the vendor source line by
> line and are transcribed correctly. The terms are right.
>
> **Where it points.** Our `delVg` comes to **3.6 mV at Vds = 1.35**, and
> PSP's own `vth` was measured moving **3.5 mV** over exactly that range.
> At this device's 85 mV/decade, 3.6 mV is a **9% current change in weak
> inversion**. So DIBL is a real part of the 2.4× climb at Vg = 0.6 —
> the climb that `FdL` was accepted for explaining *in full*, on a
> residual that already contained this omission. And the near-threshold
> sweep is precisely where enabling DIBL now overshoots, 1.003 → 1.023,
> while the two high-Vd sweeps it was predicted to help both improved
> (0.011 → 0.006 and 0.011 → 0.003). `FdL` is the term to re-examine,
> with these enabled.
>
> This is the third time on this branch that an individually correct term
> has made the fit worse, and the second time it has been resolved by
> leaving it out and recording why. The first was channel-length
> modulation — which came back as the single largest accuracy gain the
> model has taken. **Discarding correct physics because it exposes an
> error elsewhere is how that error gets preserved**; switching it off
> with the measurement attached is not.
>
> One numerical fix came out of the same work. `Vds*(1 + rat^ax)^(−1/ax)`
> overflows by `rat = 1e46` for the long device's `ax = 6.5`, and then
> `huge * (1 + inf)^(−1/ax)` is `huge * 0` — **NaN, not a large number**,
> and a diverging Newton step reaches it. Rewritten as the two forms that
> are each stable on their own side of `rat = 1`, the second being the
> algebraically identical `Vdsat*(1 + rat^−ax)^(−1/ax)`. Each arm's input
> is clamped to its own domain rather than the arms merely being selected
> between, for the reason recorded earlier: a Piecewise's derivative
> w.r.t. something used only in the discarded arm is zero, and zero times
> NaN is NaN. Finite and correctly saturated out to `Vds = 1e60` now,
> and still antisymmetric there.

> **The p-channel device — 2026-08-23.** `PspPmosLongChannel`, and the
> striking thing is how little of it is p-channel-specific.
>
> PSP is written once and run for both polarities. It converts the
> terminal voltages to an internal n-like convention on the way in
> (`PSP103_module.include:1035-1048`), writes every equation for that
> convention alone, and restores the polarity on every contribution on
> the way out (`:1697-1795`). Its **849-line geometry-scaling layer
> contains not one reference to the channel type**, and IHP's p-channel
> card is not a mirrored set of numbers — every parameter is a positive
> magnitude carrying the same signs as the n-channel one, `vfbo`
> negative on both. The polarity lives in a single `type = -1` on the
> `.model` line.
>
> Exactly **five** things genuinely differ, four of them in the kernel:
>
> * the effective-field weighting `eta_mu`, ½ for electrons and ⅓ for
>   holes, at both the source end and the midpoint (`module:736-743`);
> * `ysat = ysat/sqrt(1 + ysat)` in the saturation-voltage block
>   (`macrodefs:609-613`);
> * `zsat = zsat/(1 + thesat1·dps)` in velocity saturation
>   (`macrodefs:766-771`) — these two being one statement, that holes
>   follow `v ~ E/(1 + E/Ec)` where electrons follow
>   `v ~ E/sqrt(1 + (E/Ec)²)`;
> * `QMN → QMP` in the quantum correction (`module:727-734`), a ratio of
>   1.2515 that moves `phib` and the body factor together.
>
> Everything else is a sign. `charges_long_channel`, `sp_s`,
> `surface_state` and the Newton limiter needed no change at all.
>
> **Result on the first measurement, nothing p-channel-specific tuned:**
> long device 1.029 and 1.019, short device 1.028 / 0.939 / 1.075. For
> comparison the n-channel's first measurement was **1.29**, and took
> eight terms to bring inside a percent.
>
> And every p-channel scaled parameter matches PSP's own `lp_*` exactly
> at four geometries — `vfb`, `tox`, `ct`, `mue`, `themu`, `cs`,
> `thecs`, `rs`, `ax`, `thesat`, `alp`, `alp1`, `alp2`, `cf`, `cfb`, and
> `neff` once the quantum correction is switched off (we fold it into
> the doping where PSP applies it to `phib`; the p-channel inflation
> factor is 1.80–1.87 against the n-channel's 1.60–1.67, which is
> `QMP/QMN` doing its job). The scaling layer needed no p-channel work
> whatsoever, exactly as the vendor's does not.
>
> **How the polarity is carried matters.** `T` and `pmos` are Python
> values closed over at class-creation time, not symbols, so each
> channel type compiles to an expression exactly the size it would be if
> the other did not exist. A symbolic polarity would double every branch
> in the model and branch on something that cannot vary with bias. The
> one trap: the DSL keys compilation on `analog` appearing in the CLASS
> BODY, so a subclass that merely inherited it would silently reuse the
> n-channel expression — and assigning `analog` after the class
> statement produces a 0x0 element rather than an error.
>
> **Two overflows fixed on the way, neither cosmetic.** `safe_div`
> regularises as `b/(b² + eps²)`, which squares its argument — the right
> thing for a quantity that can change sign, and half the exponent range
> thrown away for one that cannot. `expl` is allowed to return 1e100 by
> construction, so `safe_div(1, expl(x))` inside the surface-potential
> solver overflowed at biases a diverging Newton step reaches; it is a
> floored reciprocal now. And the mobility's effective field feeds a
> power with an exponent near 2, so it is clamped above as well as
> below. An overflow to `inf` is not self-limiting: it becomes
> `inf − inf` downstream, which is a NaN in a Jacobian row. Both are
> now covered by tests that run with warnings as errors, at bias points
> lifted from a real stacked-pair solve.

Deferred, unchanged from the original research verdict:

| item | why |
|---|---|
| NQS block (18 `V<+`, 9 `idt`, branch-current readback) | genuinely inexpressible today; **gnucap ships PSP without it**, and the IHP Xyce flavour has no NQS model cards at all |
| self-heating (`psp103t.va`, thermal node + `Temp()`/`Pwr()` discipline) | needs a thermal discipline; **no IHP model card uses it** |

### Sequencing

**Phases A–D are complete** (2026-08-22): phases A, B, C and the actionable
part of D are all done,
in the order given except that A1's two halves landed either side of B.
Phase E carries the work forward to a production compact model: E1 and
E2 are done, E3 lists what a PSP103 card still needs.
What else remains is phase D, which is deferred by design, plus `@timer`
(zero call sites) and lifting the traced backend's PCNR restrictions —
it rebuilds junctions itself rather than calling the device, so
multi-junction and charge-storing participants are refused there and the
refusals are now checked rather than assumed. Original sequencing: B
first, then A1,
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
