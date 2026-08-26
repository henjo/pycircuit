# hdl.py roadmap — supporting production models, and building a library

**2026-08-24.** A review of `pycircuit/circuit/hdl.py` after the first
production compact model was built in it, with a plan for what to change
and what to build next.

## 0. How this was arrived at

Five independent passes, because the question ("how should the DSL
change?") has at least five different right answers depending on who is
asking:

1. **the compact-model author** — an audit of `compact.py`,
   `psp_kernel.py` and `juncap.py` for every place the author worked
   *around* the DSL rather than expressing physics;
2. **the circuit designer** — an empirical probe of what a newcomer
   writing a 20-line element actually hits, by making the mistakes;
3. **the librarian** — an inventory of `elements.py` against what the
   DSL can express today, and what a serious library is missing;
4. **the compiler** — a structural recon of `hdl.py`'s pipeline, its
   caches and its extension points;
5. **the stopwatch** — measured compile and evaluation cost, below.

The external target is
[designers-guide.org's Verilog-A library](https://designers-guide.org/verilog-ams/index.html),
which is the de-facto statement of what a designer expects to find.

---

## 1. Measured cost

Measured 2026-08-24 on this machine, not quoted from the record.

### Compile (once, at class definition)

| element | compile |
|---|---|
| `R` — one term, linear | 12–27 ms |
| `RC` — two terms with `ddt` | 12 ms |
| `Diode` — one exponential | 78 ms |
| `elements_hdl` — all ten, cold import | **139 ms** |
| `elements_hdl` — all **26**, cold import (2026-08-25) | **6.0 s** |
| `compact.py` — PSP n+p and the MoM capacitor | **93.5 s** |

### Evaluation

| call | hand-written | HDL | ratio |
|---|---|---|---|
| `R.i` | 0.53 µs | 1.58 µs | 3.0× |
| `R.G` | 0.07 µs | 0.62 µs | 8.8× |
| `Diode.i` | 1.90 µs | 1.91 µs | **1.00×** |
| `Diode.G` | 1.62 µs | 2.63 µs | 1.6× |

| PSP MOSFET (153 parameters, 2 internal nodes) | per call |
|---|---|
| `i(x)` | **1.22 ms** |
| `G(x)` | **23.2 ms** |
| `q(x)` | 0.58 ms |
| `C(x)` | 12.8 ms |
| `CY(x, f)` | 0.59 ms |

### Reading these

**For small elements the DSL is already the right default.** The ratios
above 1 are on operations costing well under a microsecond, where the
absolute difference is noise against a solve; for the one *nonlinear*
element in the set it is at parity on `i` and 1.6× on `G`, and it
derived the Jacobian that the hand-written twin got wrong twice
(`elements.py` records both).

**For a production compact model the Jacobian is the problem.** 23 ms
per `G` means a 100-device DC solve spends ~23 s per Newton iteration in
Jacobian evaluation alone. A compiled C compact model is ~10 µs. That
single number decides whether this DSL is usable for production models
at circuit scale, and it is the target of §2 and §3.

**93.5 s of compile is paid by every process that imports the model.**
It is cached nowhere across runs.

> **2026-08-25: and so is `elements_hdl`'s, which is no longer small.**
> The library was ten elements at 139 ms when this table was written and
> is now 26 at **6.0 s** (measured three times, 5.9-6.1 s; 3.5 s at the
> commit before the fourth batch). Every `import
> pycircuit.circuit.elements_hdl` pays it, including every test session
> and every `--collect-only`. Per class, with a warm sympy cache:
> `MosLevel1Hdl` 327 ms, its p-channel twin 302 ms, `OpAmpHdl` 127 ms,
> `GummelPoonNpnHdl` 496 ms. This is the same problem §2 records for
> `compact.py` at a smaller scale, and it has the same answer -- the
> persistent on-disk compile cache §3's S2 names as its companion item.
> Nothing is broken; it is a number that has grown 43x without anyone
> re-reading it, which is exactly the shape §11 keeps recording.

---

## 2. The structural finding: the chained-path cliff

`generate_code` forks twice on `if chain_defs:` —
i.e. on whether the model used `var()`. Every production model does.
The chained arm silently loses **five** capabilities:

| lost on the chained path | line | consequence |
|---|---|---|
| `G`/`C` as symbolic matrices | 1702 | set to `zeros` |
| `pure_spec` | 1780 | **no JAX form** |
| `sym_spec` | 1789 | no symbolic toolkit |
| `pcnr_spec` | 1832 | no PCNR participation |
| `const_G`/`const_C` | 1906 | no constant-stamp caching |

Verified on the shipped model: `PspMosLongChannel._hdl_info['chained']`
is `True` and `pure_spec` is `None`.

**So the branch's headline feature does not apply to its flagship
model.** `solve_batched` — the vmapped, jit-compiled path measured at
22.5× on a 512-lane sweep — cannot take a PSP transistor. A corner sweep
over the one device anybody would want to sweep runs on the CPU loop.

The fork exists for a real reason, recorded in the comment above the
first `if chain_defs:` in `generate_code`:
`Matrix.jacobian` differentiates the assembled expression *tree*, which
for a model that reuses intermediates is exponential in nesting depth
even though the DAG is linear. The chain path avoids that by
differentiating the let-chain forward instead. The reason is sound; the
consequence — five capabilities gated on one flag — is not.

---

## 3. Five structural changes, ranked

### S1 — give the chained path a pure form  *(unlocks JAX for real models)*

A let-chain **is** a straight-line program, and a JAX-traceable function
is exactly a straight-line program. The chain compiler already emits
Python source (`_chain_compile`, `fn._src`); emitting the
same chain against `jnp` instead of `numpy` is a printer change, not a
new derivation. Jacobians then come from `toolkit.jacobian` (autodiff),
which is what the pure path already does for non-chained elements —
and autodiff over a straight-line program of ~700 intermediates is
exactly the case JAX is good at.

Expected payoff: batched corner sweeps over real transistors, and a
plausible route to the 23 ms `G` problem, since XLA fuses what numpy
cannot.

*Risk:* `_recip2`/`_rdiv` are mapped to numpy implementations
(the `pure_spec` assignment in `generate_code`) — pure arithmetic, so
they trace, but every future
primitive must be checked. `Piecewise` → `jnp.where` is already handled.

### S2 — generalise the variant mechanism to elide inert blocks

> ⛔ **REFUSED ON MEASUREMENT, 2026-08-26 (§29).** The pitch below is
> kept as written, because the *reason* it was wrong is the useful part.
> Eliding the blocks a default PSP card makes inert — 15% of the
> generated Jacobian, exponentials included — measures **1.00x on numpy
> and 1.00x in C**, bit-identical. The blocks are cheap: they contain
> **zero** smoothing-primitive calls, and a compact model's cost is its
> **regularisers**, not its physics.
>
> The paragraph beginning *"It is also the most direct attack on
> evaluation cost"* is the specific claim that failed. What *was* the
> direct attack is **§30–31's card-constant folding** (`hdl.fold_card`,
> 1.54x numpy / 1.48x C) — which subsumes this idea anyway, since with
> the parameters numeric a kernel's own build-time guards fire on their
> own.
>
> The **correctness** half below still stands: a block whose parameter is
> zero is still compiled and still evaluated, and `0 × inf = NaN`.

`_collapse_mask_of` / `_collapse_variant` already
recompile-and-cache a class per parameter mask and retarget instances by
assigning `__class__`. The mask is an **opaque tuple**; nothing in the
machinery knows it means "collapse".

Generalising it to "elide any block whose parameters make it inert"
needs: a new `Statement` subclass (`Collapse` is the template, including
its parameters-only validator in `Collapse.__init__`), a partition and
an application step beside the ones `generate_code` runs for
`collapse_conds`, a widened mask tuple, and a smarter freeze check —
because a guard that removes *terms* but not *nodes* does not change the
unknown count, so it may legally change after build.

*(Anchors were `hdl.py:971`, `:1172`, `:1189` when this was written on
2026-08-24. All three are long stale — `Collapse` is now near line 2481
— and are replaced with symbol names per `record-upkeep`.)*

This fixes a trap that cost real time: **a block whose parameter is zero
is still compiled and still evaluated**, because the element compiles
from symbolic parameters, so `if param != 0` is true at build time
whatever the value — and `0 × inf = NaN`. The model currently guards
these blocks by hand with build-time tests of the form
`if not isinstance(p, sympy.Expr) and p == 0.0:`, which only fire for
the *default* build; an instance whose card sets the parameter to zero
still evaluates the whole block.

*(Recorded as "twelve blocks" on 2026-08-24; counted again 2026-08-26 it
is **nine**, in `psp_kernel.py`. The shape is unchanged.)*

It is also the most direct attack on evaluation cost: a card that
switches off gate current, avalanche and half the parasitics should not
pay for them on every Newton iteration.

*Risk:* variant explosion. Mitigate by hashing the *reduced statement
list* rather than the mask, so masks that produce identical models share
one variant, and by keeping compilation lazy (it already is). A
first-instance stall of tens of seconds for a PSP-scale model is
user-visible latency, not a background cost — which makes a persistent
on-disk compile cache (§6) a companion item rather than a nicety.

### S3 — automatic both-arms domain clamping

Both arms of a compiled conditional are evaluated, so an unselected arm
must still be finite — *and* its derivative must be, because the chain
rule multiplies a zero partial by that arm's derivative and `0 × NaN` is
`NaN`. Today the author does the compiler's domain propagation by hand,
at **fifteen-plus sites**, each one a clamped copy of an input fed to an
arm that never selects it.

A `hdl.select((cond, expr), ...)` that substitutes the branch condition
into each arm — clamping the free symbols the condition constrains —
would absorb the whole class. Even a partial version handling the common
shapes (`x > c`, `x < c`, `x*x < c`) would remove most of them.

This is the single largest source of *silent* bugs in the model: three
separate overflows in the last block written, none visible at ordinary
bias.

### S4 — hold every intermediate automatically

`var()`'s docstring still says naming is *"optional and only affects
generated-code readability"*. The model's experience contradicts it
flatly: bare sub-expressions get substituted into one another before
printing, so the compiled form evaluates products the written form never
creates, and the n-channel device lost its Jacobian at `Vd = 1e26`
*purely from evaluation order*. There are **342** hand-tagged `_v()`
calls in `psp_kernel.py`.

If the compiler held every non-atomic sub-expression automatically
(structural hashing at build time), `var()` would go back to being the
naming convenience its docstring claims. Until then the docstring is
wrong and should say so.

### S5 — a parameter namespace instead of `__globals__` mutation

`analogfunc.__globals__.update(...)` (now `_analog_function`) writes parameter
symbols into the *defining module's* globals — permanently, because
`copy.copy` of a function is a no-op. Three measured consequences: 63
`# noqa: F821` in `compact.py` alone; a parameter named `vt` shadows the
`hdl.vt()` helper for every later class in the module; and an undeclared
or misspelled parameter silently resolves to a *previous element's*
symbol, failing much later as a `NameError` from inside
`<lambdifygenerated-12>`.

Passing a namespace — `def analog(d, g, s, b, p)` with `p.vfb` — removes
all three, makes the DSL legible to a linter and an IDE, and gives the
model a supported way to ask for its own parameters (which it currently
does by reading the compiler's injection back out of `globals()`).

### Kernel additions — small, and each removes a rule the author must remember

| addition | replaces |
|---|---|
| `maxc`/`minc` as `Function`s with their own `fdiff` | the undocumented "`Max` only on an atom" rule, and ~50 sites obeying it by construction |
| `safe_pow(b, e)` with a documented base range | six hand-floored, hand-capped power bases |
| `sign` with `fdiff → 0`; real-domain `Abs` | `sign` failing to *compile* (`DiracDelta`), and conditions rewritten algebraically to dodge `Abs` |
| `softplus`, `mne`/`mxe`, `bias_mod`, `p3` | generic helpers currently hand-built in `psp_kernel.py`; `p3` is outright duplicated |
| `limit_delta(probe, vmax)` | 35 hand-written lines of FET Newton limiting per model |
| **range contracts** on every `safe_*`, and a build-time "this guard can never fire" check | the exponent-range audit the author now does in prose, per call site |

---

## 4. The ergonomic surface

The pattern, stated in one line: **the compiler defends its semantic
invariants superbly and its syntactic surface not at all.** The errors
that are already right — switch branches, state-scaled `ddt`, collapse
conditions, instparams drift — share a shape (what you wrote, why it
cannot work, what to write instead) and should be the template. The
errors a newcomer actually reaches have no message at all:

| mistake | what happens today |
|---|---|
| forgetting `return` (Verilog-A has no `return`) | `TypeError: 'NoneType' object is not iterable`, from inside the compiler |
| `def analog(self, plus, minus)` | **silently** creates a terminal named `self`; every connection shifts by one |
| declaring `terminals = (...)` in the class body | **silently discarded**; pin order comes from the `analog()` signature |
| an internal `Node('out')` in an element with an `out` terminal | **silently** merges into the terminal |
| `if v > 0` on a circuit quantity | raw sympy `TypeError: cannot determine truth value of Relational` |
| `math.exp` / `numpy.exp` | `TypeError: Cannot convert expression to float` |
| wrong number of connection nodes | accepted silently |
| unknown parameter name at instantiation | `KeyError: 'parameter R not in parameter dictionary'` — no class, no valid names |

Six checks in `generate_code`, a non-mutating globals
injection, and a public `explain()` (terminal list, `x`-vector layout,
symbolic `i`/`q`/`G`/`C`, generated source — all of which exist
privately today) would remove essentially every silent or opaque failure
a designer can hit.

Add to that a public `check_jacobians(element, x)` — `G` against `i` and
`C` against `q` by finite differences, plus a finiteness scan. The
documentation warns beautifully about NaN Jacobians and then hands the
reader no instrument.

---

## 5. The library

### First: correct the capability record

`hdl.md` §6's "Not implemented" table and the matching paragraph in
`doc/src/circuit/hdl.rst` both predate Phase B/C/D and **understate the
DSL by six capabilities**. Verified present and tested:

| §6 says missing | actually at | test |
|---|---|---|
| `$limit` | `hdl.py`, `limit_pnj` | `test_limit_pnj_compresses_a_wild_step` |
| `laplace_nd`/`laplace_zp` | `:621` | `test_laplace_nd_matches_the_analytic_response` |
| `@cross` | `:891` | `test_cross_lands_timepoints_on_the_crossing` |
| `$param_given` | `:774` | `test_param_given_selects_a_formulation` |
| node collapse | `:931` | `test_node_collapse_removes_the_internal_node` |
| AC excitation | `:791` | `test_ac_stim_drives_a_small_signal_analysis` |

Anyone planning migration work off that table would conclude the DSL
cannot do things it does. Fix it before planning anything else.

**Genuinely still missing:** `absdelay`, `transition`/`slew`, `zi_*`,
`last_crossing`, `noise_table`, `@timer`, and conditional node collapse
(refused by design).

### What the behavioural half needs

Reading the designers-guide models directly, the phase-frequency
detector — the archetype — uses `@cross` (have it), `@initial_step`
(trivial), `transition()` (missing), and **`integer state` held across
evaluations** (missing). That last one is the real gap: `idt`/`idtmod`
are *continuous* states; a PFD, flip-flop, counter, divider and
quantizer all need a *latched discrete* state updated on an event.

**Three additions — `transition()`, event-latched discrete state, and
`@initial_step` — unlock the sampled-data half of that library**:
PFD, charge pump, divider, DFF, latch, counter, quantizer, S/H,
comparator, ADC/DAC. That is the highest-leverage feature cluster for
the library, and it is well-scoped: `Cross` already provides the event
and cost ~50 lines.

### The fifteen worth building first

Ranked by (value × already-buildable) ÷ size. All are writable with
today's DSL unless noted.

| # | model | family | ~lines | why now |
|---|---|---|---|---|
| 1 | **Self-heating thermal node** + Foster/Cauer ladder | electrothermal | 60 | ~10 lines per device, no new operator: the thermal node's "voltage" is ΔT, the device contributes `v·i` into it and reads `TEMP + V(th)` back. Retroactively upgrades every device, and closes a named gap in `PspMosLongChannel`. **DONE**; `GummelPoonNpnThermalHdl` added 2026-08-25, and the ladder is measured -- the thermal pin is a PARALLEL path to the device's own `rth`, not a series one (§12.6d) |
| 2 | Full SPICE diode (N, RS, TT, CJ0/VJ/M, BV, shot+flicker noise) | semiconductor | 80 | today's `Diode` is a bare exponential; `juncap.py` already holds the kernels |
| 3 | SPICE Gummel-Poon BJT | semiconductor | 300 | the most conspicuous absence in the repo |
| 4 | Opamp macromodel (Boyle-class) | behavioural | 150 | deletes `macromodels.OpAmp`'s runtime-`import jax` lambdas; the most-used block in any testbench. **DONE 2026-08-25** as `OpAmpHdl`; the `import jax` claim is TRUE and still live, but it is the smaller half -- see §11 |
| 5 | Quartz crystal / BVD resonator | passive | 40 | nothing in the repo can simulate an oscillator's frequency-setting element |
| 6 | VCO with real phase noise | behavioural | 40 | noise injected into the *frequency*, so jitter accumulates as a random walk — a capability the repo's own survey found in no other simulator |
| 7 | Comparator with hysteresis | behavioural | 60 | the **first production user of `Cross`**, which Phase C1 built and nothing consumes |
| 8 | Skin-effect R + ferrite bead | passive | 100 | the **first production users of `laplace_zp`**, likewise unconsumed |
| 9 | EKV 2.6 MOSFET | semiconductor | 150 | one tenth of PSP's size, all-region, and an independent cross-check on the surface-potential kernel |
| 10 | MOS Level 1 + Level 3 | semiconductor | 330 | legacy netlists and teaching. **Level 1 DONE 2026-08-25** as `MosLevel1Hdl`/`MosLevel1PmosHdl`, and it is `limit_together`'s first production user; Level 3 not started |
| 11 | Memristor (HP + Biolek window) | passive | 70 | the clearest small demonstration that `idt` is a real DAE state |
| 12 | Mixer + PLL loop filter + charge pump + continuous ÷N | RF/behavioural | 140 | with #6, a complete usable PLL minus the tri-state PFD, from three existing capabilities that have no library user |
| 13 | Photodiode / solar cell / LED | semiconductor | 210 | falls out of #1 and #2; the optical node is free in MNA |
| 14 | Shot noise across the existing junction devices | noise | 40 | hours of work, and it makes `CY` mean something beyond `R` |
| 15 | Statz/Curtice MESFET + Angelov HEMT | semiconductor | 320 | the RF gap; both are `tanh`-shaped and map almost line for line |

### Migration of `elements.py`

Most of it can move, and several elements get *better* by moving —
`VCVS_limited` and the switches gain a derived Jacobian; `BSource`'s
Python callables become a real expression with exact derivatives and
noise. Four should not move, and the reasons should be recorded so
nobody rediscovers them:

| element | why not |
|---|---|
| `TLine` | needs `absdelay` (history at `t − TD`); its interpolation and monotone limiter are solver-coupled work no contribution expresses |
| `VPWL`/`IPWL` | a list-valued parameter changes the expression, i.e. per-instance topology |
| `Nullor` | expressible, but as an unproven idiom that reads far worse than the four-entry stamp |
| `BSource` | callable parameters — nothing symbolic to compile (and superseded in purpose) |

---

## 6. Sequence

**Phase 1 — make the record true** (hours). Correct §6 and `hdl.rst`.
Everything below is planned against what the DSL *actually* does.

**Phase 2 — the designer's surface** (days, low risk, no architecture).
The six syntactic checks, the non-mutating globals injection,
`explain()`, `check_jacobians()`, and the two documentation sections a
designer needs and cannot find: "writing your first element" and "when
it goes wrong". This is the cheapest work with the widest audience.

**Phase 3 — the kernel additions** (days, low risk, contained).
`maxc`/`minc`, `safe_pow`, `sign`, the promoted generic helpers, range
contracts. Each removes a rule the author must remember, and the
`differentiable-numerics` skill records why each exists.

**Phase 4 — S1, the pure form for chained models** (weeks, high value).
The single change that would let the branch's headline feature reach its
flagship model, and the most promising attack on the 23 ms Jacobian.
Gate it on a measurement: batched PSP sweep against the CPU loop.

**Phase 5 — S2, generalised variants** (weeks, medium risk).
Needs the compile cache to land with it, or first-instance latency
becomes the new complaint.

**Phase 6 — the library**, continuously, starting with the thermal node
and the diode. Each model is also a test of Phases 2–3: if writing it
still needs a workaround, that workaround is the next item.

**S3 and S4 are research.** Automatic domain clamping and automatic
intermediate holding are the two changes that would most improve the
compact-model author's life, and both are genuinely hard. They should be
prototyped against `psp_kernel.py`, which is now a 1671-line regression
test for exactly these properties.

---

## 7. What this plan does not claim

The 23 ms Jacobian is *diagnosed*, not solved. S1 and S2 are the two
plausible attacks and both are unproven; it is equally possible that the
answer is a different code generator (`_chain_compile` and
`_ChainPrinter` are a self-contained replacement point) targeting C or LLVM rather than Python source. That should
be measured before it is chosen — the same discipline the model itself
was built under.

---

## 8. Addendum, 2026-08-24 (same day): where the 23 ms actually goes

§1 measured whole elements and §2 blamed the chained code generator.
Neither is enough to act on, because neither says what the evaluation
is *spending its time on*. `benchmarks/hdl_model_cost.py` (added with
this section) takes the ladder apart, and the answer changes what
Phase 4 should be.

### False start, recorded because it is instructive

The first pass measured primitives through `sympy.lambdify(modules=
['numpy'])` and found that sympy lowers `Piecewise`, `Min`, `Max` and
`Heaviside` all to `numpy.select` — at 6.03 us per scalar call against
`numpy.minimum`'s 0.60 us. It concluded that `select` dispatch was the
cost centre of production models, and measured a 184x speedup from
`modules=['math']` on `d(expl)`.

**That conclusion was wrong, and `hdl.py` already contained the
refutation.** `_ChainPrinter` exists precisely to
replace `select` with nested `numpy.where`, and its docstring records
its own profiling: *950 select calls per evaluation, 60% of the total
runtime* on this same MOSFET, plus 18% more from `functools.reduce` on
`Min`/`Max`. Every production model takes the chained path, so none of
them were paying the cost being "discovered". Dumping the generated
source confirms it: the chained diode's `G` emits **zero** `select`.

The lesson is the one this project keeps relearning: a measurement of
the wrong code path is not evidence about this one. `lambdify(modules=
['numpy'])` is not what compiles a production model, and the difference
was one `_src` dump away.

### The flagship, measured directly

PSP103 n-channel, saturation bias, this machine:

    compile+instantiate                        105.9 s
    i                                           1414 us
    G                                          26523 us     (18.8x)

Generated source, and the primitives it actually emits per call:

    i:   676 lines     940 numpy calls
         where 187, minimum 343, maximum 359, exp 51, _rdiv 209
    G:  2675 lines   14030 numpy calls
         where 7096, minimum 3558, maximum 2905, exp 303,
         reduce 168, _recip2 1132, _rdiv 1775

Priced at the measured per-call scalar costs (`where` 0.94 us,
`minimum`/`maximum` 0.60 us, `exp` 0.10 us):

    i:  predicted dispatch   602 us of  1414 us measured   (43%)
    G:  predicted dispatch 10578 us of 26523 us measured   (40%)

**So numpy scalar dispatch is ~40% of the Jacobian, and the arithmetic
spread over 2675 lines of Python is the other ~60%.**

### What follows

**Scalar code generation is worth ~1.7x, not 184x.** Emitting straight-
line scalar Python (`if`/`else`, `math.*`) removes the 40%, taking `G`
from ~26.5 ms to ~16 ms. That is real and cheap, but it does not make
the model usable at circuit scale, and the 184x figure above must not
be quoted for it.

**The number that matters is 14,030.** One device's Jacobian issues
fourteen thousand dispatched operations. No amount of printer tuning
changes the count — only the price per item. At ~1-5 ns per operation a
compiled backend would put this at 15-70 us, which is the 300-1000x
that would actually change what the DSL can be used for. §7 listed a C
or LLVM code generator as one of three candidates; on these numbers it
is the only one that reaches the goal, and `_chain_compile` /
`_ChainPrinter` remain the self-contained replacement point.

**A live defect, found by this measurement.** 168 `reduce(` calls
survive in `G`. `_ChainPrinter._minmax` exists to eliminate exactly
those and its docstring claims 18% of runtime for them, so some
`Min`/`Max` is being printed by a path the custom printer does not
intercept — a Min with more than two arguments, or an expression
printed by the stock `NumPyPrinter` rather than `_P`. Worth finding: it
is a bug in a fix that was already made, and it is invisible to every
test that checks values rather than emitted source.

### The instrument

`benchmarks/hdl_model_cost.py:emitted()` and the `_src` dump are the
tools here. A primitive is only as cheap as what the printer turns it
into, and that is invisible from the sympy side — so a kernel primitive
needs a test asserting **what it emits**, not only what it evaluates
to. Neither of the two errors above would have survived one such test.

**And a caution against the tempting version of all this.** Do not
"fix" the cost by removing regularisers. They cost what they cost
because both arms of a conditional are always evaluated, and that is
exactly what keeps the Jacobian finite. The target is the dispatch and
the operation count, never the safety.

---

## 9. The chained path costs TWO capabilities, not one (2026-08-25)

§2 recorded that a model taking the chained path gets no `pure_spec`,
so `solve_batched` cannot take it. Measured while correcting a stale
comment in `elements_hdl.py`, there is a second one with the same root
cause. Three elements differing by exactly one thing:

    flat, no charge         chained=False   pcnr_i=True
    flat, WITH charge       chained=False   pcnr_i=True
    chained, with charge    chained=True    pcnr_i=False

**`var()` disqualifies a model from PCNR.** The shape detector recovers
a junction's `IS` and `VT` by reading `exp(arg)` out of the expression,
and a let-chain hides it behind an intermediate symbol.

The mechanism is the same one as `pure_spec`: the chained code
generator produces something no downstream consumer can pattern-match
against. So the ledger for `var()` is:

    loses   pure_spec   ->  no JAX batching (solve_batched)
    loses   pcnr_i      ->  no PCNR limiting
    costs   ~3x         ->  on the Jacobian, vs the eager path

against a gain of forward-accumulated derivatives that is what makes a
model the size of PSP compile at all. Every production model takes it.

**This raises the value of S1** (§3): a pure form for the chained path
would recover both capabilities, not just batching. It does not change
the §8 conclusion about where the 23 ms goes -- that is dispatch and
operation count, and S1 addresses neither.

> **SUPERSEDED IN PART, 2026-08-25 (same day).** The `pcnr_i` half of
> that ledger is fixed -- see §10.2 as implemented. The detector now
> walks the chain instead of reading a flattened expression, so a
> chained model qualifies:
>
>     chained, with charge    chained=True    pcnr_i=True
>
> The ledger for `var()` is now:
>
>     loses   pure_spec   ->  no JAX batching (solve_batched)   STILL TRUE
>     loses   pcnr_i      ->  FIXED
>     costs   ~3x         ->  on the Jacobian                   STILL TRUE
>
> S1's value is correspondingly back to what §2 said: batching. The
> paragraph above is kept because the *reasoning* -- that the chained
> generator produces something downstream consumers cannot
> pattern-match -- is what predicted the other four bugs of this shape
> found the same week, and it is still the right thing to check first.

**And a correction it replaces.** The `elements_hdl.py` comment said
the disqualifier was the CHARGE, because the layer refuses a
participant that stores it. That was true when written and is not true
now: `pcnr.py`'s "CHARGE STAYS IN THE MNA BLOCK" comment in
`augmented_system` records the refusal being removed, precisely because
it "cost every junction device with capacitance, which is to say every
real one". The conclusion in that comment -- use `expl`, not `limexp`
-- was right throughout; only its stated reason had rotted.

That is the fourth claim this campaign has found right in conclusion
and wrong in reason, after the `Max`-on-an-atom rule, `limexp`'s
"is it a junction" test, and `_ChainPrinter`'s `reduce` premise. The
pattern is worth naming: **a reason decays faster than the conclusion
it supports, and nothing fails when it does.** Recording the
measurement next to the rule is what makes the decay visible later.

---

## 10. PCNR: what it would take to admit chained models, and FETs

Two questions this plan did not cover, worked out against the code on
2026-08-25. Both are about the same layer, and they have very
different costs, so they are separated here rather than bundled as
"improve PCNR".

### 10.1 First, a premise to discard: charge is not the blocker

It is easy to reach for the charge as the reason a real diode cannot
join PCNR, and `elements_hdl.py` said exactly that until this section
was written. `pcnr.py`'s "CHARGE STAYS IN THE MNA BLOCK" comment in
`augmented_system` records the refusal being removed, on the
grounds that it "cost every junction device with capacitance, which is
to say every real one", and `test_pcnr_charge.py` proves the charge
case by finite differences.

Measured, three elements differing by one thing each:

    flat, no charge         chained=False   pcnr_i=True
    flat, WITH charge       chained=False   pcnr_i=True
    chained, with charge    chained=True    pcnr_i=False

**`var()` is the blocker.** See section 9.

### 10.2 Admitting chained models

The gate is explicit, in `generate_code`:

    if (not chain_defs and not states and not vbranches
            and len(terminalnodes) >= 2):

and the reasoning above it is sound: the detector reads
`f.atoms(sympy.exp)` off the assembled expression, so behind a `var()`
the exponential is an opaque symbol and it would find none — silently
declaring the device un-limitable. Refusing loudly beats that.

**The wrong fix is to inline the chain before detecting.** Inlining is
the thing the let-chain exists to prevent, and on a chain the size of
PSP's it does not terminate in useful time. It is also unnecessary: the
detector never needs the flattened expression, only two facts about it,
and both are available by walking the chain in linear time.

**(a) Which branch voltages does the contribution reach?**
`_chain_compile` already prunes a chain to the
definitions an output actually reaches — `defmap`, `wanted`, `_leaves`.
Reuse that pruning, then take the union of `xsyms` over the reachable
definitions. Exactly `{kp, km}` means "a function of one branch voltage
alone", which is the condition the current free-symbol test expresses.

**(b) The exponential scale.** For each `exp(arg)` in the reachable
sub-chain, the scale needs `d(arg)/dv_branch`. That is forward
accumulation over the chain — the machinery
`_chain_compile(want_jacobian_of=...)` already implements, linear in
the number of definitions rather than exponential in nesting depth. If
every such derivative is free of `xsyms` and they all agree, that is
the single `VT` the detector is looking for.

Then compile `pcnr_i`/`pcnr_didv` through `_chain_compile` on the
pruned chain with the branch voltage as the sole `x` argument, in place
of the `sympy.lambdify` that builds `pcnr_funcs`.

Contained: roughly 60-100 lines, no new concepts, and the three-element
table in 10.1 is already its regression test.

**What it does NOT deliver, stated plainly.** It does not admit
`DiodeSpiceHdl`. That model also carries an internal node for `rs`, so
its current spans two branch voltages and the free-symbol test still
refuses it. Only the `rs = 0` collapsed variant would qualify. A
worked example for this change should therefore be a chained diode
*without* series resistance, and the limitation should be said out loud
rather than discovered by whoever tries it on the real model.

> **BOTH SENTENCES ABOVE ARE WRONG, corrected 2026-08-25 by building
> it.** They are kept because how they are wrong is the useful part.
>
> *The claim.* `DiodeSpiceHdl`'s **default card qualifies**. `rs`
> defaults to 0, and the collapse variant removes the series
> contribution entirely -- so the shipping default IS the collapsed
> variant. "Only the `rs = 0` variant would qualify" was the accurate
> half; "does not admit `DiodeSpiceHdl`" was not.
>
> *The reason.* Measured: **each** of the two contributions is a
> function of its own branch voltage alone -- that is exactly what the
> internal node buys, and the free-symbol test passes. What refuses
> `rs = 2` is that the series resistor's current is not exponential,
> and the rule is *every* current, not *some*.
>
> *And one more.* "Exactly `{kp, km}` is the chained equivalent of the
> free-symbol test" is not equivalent -- it would wrongly accept
> `f(V(a) + V(b))`. What was implemented substitutes the branch voltage
> into each definition separately and requires no x-symbol to survive
> in any of them, which *is* equivalent.
>
> Sixth instance in this campaign of a claim right in conclusion and
> wrong in its stated reason, and the third of those was mine.

**`PspMosLongChannel` still does not qualify -- and not for the reason
this section gives.** It is refused at its **first contribution**,
`I(g,gi) <+ V(g,gi)/rg`: a plain gate resistor with no exponential, so
the single-scale test fires and the detector never reaches the drain
current. The structural obstacle §10.3 describes is real, but it is not
the refusal that fires, and anyone reasoning from §10.3 about what to
change would have optimised the wrong thing.

### 10.3 FET limiting

Here the obstacle is structural rather than a matter of degree, and it
is worth stating exactly.

**PCNR's unknown is one scalar per two-terminal branch.**
`solve_dc` builds `v_lim` as one entry per junction, `x[ra] - x[rb]`;
`augmented_system` reads `v = float(v_lim[idx])`; and
`pcnr_i(v, ...)` takes that scalar. A MOSFET's drain current is
irreducibly `f(vgs, vds, vbs)`. No scalar `v_lim` expresses it, so
`f.free_symbols & set(xsyms)` fails **by construction** — the detector
is not missing a case, it is correctly reporting that the shape does
not fit.

So there are two different deliverables here.

**(a) FET limiting on the ordinary Newton path.** This is the cheap one
and probably the one worth having first. That path is `$limit` /
`limit_spec`, which is per-PROBE and carries none of PCNR's
scalar-branch assumption. `limit_pnj(probe, IS, VT)` is
the existing declaration; `limit_fet(probe, vto)` and
`limit_vds(probe)` belong beside it, implementing SPICE's `fetlim` and
`limvds`.

One real edit is needed rather than a pure addition: `limit_spec`
hardcodes an `(IS, VT)` pair per probe, so it has to
carry a limiter KIND plus that kind's own parameters instead. Days of
work, and it gives MOSFET limiting in ordinary runs — which is where it
matters most, since `pcnr=False` is the default.

**(b) FET limiting inside PCNR.** This is a genuine extension, not a
feature. PCNR has to become vector-valued per device:

    junctions      -> per-device, carrying a LIST of terminal pairs
    v_lim          -> a block per device, not one entry
    pcnr_i/didv    -> take a vector; return an (n_rows x n_lim) block
                      for J_MNA/lim
    refine()       -> call a DEVICE-level limiter that receives all of
                      that device's voltages at once

The last point is not cosmetic. SPICE's sequence is order-dependent:
`fetlim(vgs)` runs first, then `limvds(vds)` uses the ALREADY-LIMITED
vgs. A per-branch limiter cannot express that, so the vector form is
forced by the physics of the limiter, not merely convenient.

(b) also makes 10.2 a prerequisite, since any real FET is chained.

> **10.3(b) SPLITS IN TWO, and the first half landed 2026-08-25.**
> The list above bundles a *device-level limiter* with *vector PCNR*,
> and they are not the same project. The device-level limiter is
> `limit_together`, on the ORDINARY Newton path, and it cost a day:
>
>     vgs, vds = limit_together(limit_fet(bgs.V, VTO), limit_vds(bds.V))
>
> Only the PCNR half — vector `v_lim`, block-valued `pcnr_i/didv`, a
> `refine()` that hands a device its whole block — is still a project.
>
> **Two claims above are wrong, and both are wrong in the same
> direction: they name SPICE's ORDERING as the thing that forces the
> device-level form.**
>
> *First.* "The vector form is forced by the physics of the limiter,
> not merely convenient." Measured, and it is not. SPICE's ordering is
> now implementable (`limit_together(..., sequential=True)`) and it is
> NOT better: over 48 operating points of the cascode in
> `test_device_limiter.py`, 1029 Jacobian evaluations against 927 for
> the order-independent form, and equal at the reference point. It is
> a difference, offered for SPICE comparability; it is not physics.
>
> *Second, and this is the one that matters.* What actually forces the
> device-level form is the WRITE-BACK, on a device with more than two
> probes. SPICE's own MOSFET declaration is four probes — `fetlim` on
> `vgs`, `limvds` on `vds`, `pnjlim` on each bulk junction — over four
> terminals, and the per-probe write-back gives each probe a terminal
> of its own, so the fourth finds both of its terminals claimed and
> `generate_code` REFUSES TO COMPILE THE MODEL. The four-terminal
> MOSFET is not expressible per-probe at all. Grouped it compiles, the
> triangle `(b,s)-(b,d)-(d,s)` is resolved at run time by dropping the
> smallest constraint, and on a body-biased cascode it takes 896
> Jacobian evaluations against 3478 unlimited, with 6 of the 48 points
> not converging at all without it.
>
> **What the device-level form does NOT buy is iterations.** On the
> two-probe cascode it ties the per-probe form at 34 of 48 points and
> costs exactly one iteration at the other 14 — 927 against 909. Said
> plainly because §12.1's expectation was the opposite, and because the
> number that expectation rested on was stale (see 12.1).
>
> **The design, in one paragraph.** A group's probes are read in a
> canonical order (largest correction first, ties by row index), each
> reading the shift the earlier ones left, and then written back as a
> whole: the probes are a graph over the device's terminals, a maximum
> spanning forest by correction size keeps every constraint that does
> not close a cycle, and in each component the LEAST DRIFTED node is
> held while every other node is derived from it. If no probe bit,
> nothing is written. **Reading the probes independently instead —
> every one from the unlimited vector — looks like the natural
> device-level thing and OVER-CORRECTS**: probes share terminals, so
> one write often satisfies the next probe outright, and enforcing a
> correction that is no longer needed is what pushes the write onto a
> source-pinned node. Measured, that mistake costs 1040 against 927.

### 10.4 Sequence, and the gate on the expensive one

1. **10.3(a)**, FET limiting on the ordinary path. Independent of
   everything else here, and the only item that helps PSP today.
2. **10.2**, the chained detector. Contained, well-tested by
   construction, and it removes one of the two capabilities `var()`
   costs (section 9).
3. **10.3(b)**, vector PCNR. Its own project.

**2026-08-25: items 1 and 2 are done, and so is the device-level
limiter that was bundled into item 3.** What remains of item 3 is PCNR
alone. The gate below is unchanged and is now the only thing standing
between here and it — and the limiting to measure against is
`limit_together`, not the per-probe form, since a four-terminal MOSFET
can only be declared with it.

**Gate item 3 on a measurement, the way section 8 gates the Jacobian
work.** PCNR's value for a FET is real in principle — several
transistors sharing a node is exactly the clash over a shared
linearisation point that PCNR exists to remove — but "real in
principle" is what this campaign has repeatedly found to be the
weakest kind of claim. Build a convergence case that classical
limiting handles badly and PCNR should handle well, measure it with
the limiting from item 1, and only then decide whether the
architecture is worth rebuilding around vector unknowns.

---

## 11. Status, 2026-08-25

This plan was written on 2026-08-24 and parts of it were implemented
the next day, so several sections above describe as PENDING things that
now exist. Rather than rewrite them — the reasoning is worth keeping in
the form it was argued — here is what actually landed.

    section  item                                    status
    3        S5  parameter namespace                 DONE, both halves: the leak
                                                     (_analog_function) and the
                                                     surface (params_as, sec. 17)
    3        kernel: maxc/minc                       DONE
    3        kernel: safe_pow                        DONE
    3        kernel: sign, real-domain Abs           DONE
    3        kernel: softplus/mne/mxe/p3 promoted    DONE
    3        kernel: range contracts                 DONE
    3        kernel: limit_delta                     not started
    3        S1  pure form for chained path          REFUSED on measurement
                                                     (sec. 22/23): the two fast
                                                     paths are disjoint (22
                                                     chained / 15 pure_spec /
                                                     0 both), so S1 buys only
                                                     batching on top of a C
                                                     backend already at 212x
    3        C backend for the EAGER path            REFUSED on measurement
                                                     (sec. 23): the elements'
                                                     own arithmetic is 18.2%
                                                     of an eager transient,
                                                     43.3% at 121 devices --
                                                     ceiling 1.76x against the
                                                     C backend's measured
                                                     212x, for the same build
    30/31    card-constant folding                   DONE (sec. 31):
                                                     `hdl.fold_card`,
                                                     1.54x numpy / 1.48x C
                                                     on PSP, compiles in
                                                     17 s against 28.  NOT
                                                     bit-identical -- a
                                                     validated model must be
                                                     re-validated folded
    3        S2  generalised variants                REFUSED on measurement
                                                     (sec. 29): eliding the
                                                     blocks a default card
                                                     makes inert -- 15% of the
                                                     generated Jacobian,
                                                     exponentials included --
                                                     is 1.00x on numpy and
                                                     1.00x in C, bit-identical.
                                                     The cost is the
                                                     REGULARISERS (minc/maxc,
                                                     4.5 of 17.2 ms), not the
                                                     physics
    3        S3  automatic domain clamping           DONE (sec. 21): hdl.select,
                                                     6 of 20 sites converted,
                                                     bit-identical
    3        S4  automatic intermediate holding      research
    4        the six syntactic checks                DONE (nine of them)
    4        explain()                               DONE
    4        check_jacobians()                       DONE
    5        thermal node, full SPICE diode          DONE
    6        Phase 1 make the record true            DONE
    6        Phase 2 designer's surface              DONE
    6        Phase 3 kernel additions                DONE
    8        printer/scalar codegen                  DEAD -- see §13/§14
                                                     (~1.2x now, overtaken
                                                      by the kernel work)
    16       compile cache                           DONE: import 5.1 s -> 0.14 s,
                                                     PSP 41 s -> 0.6 s
    14       backend: the spike                      DONE (sec. 19): C at -O2 is
                                                     219-229x, 0 ulp on 4985 pts
    14       backend: BUILT                          DONE (sec. 20): 20-270x,
                                                     bitwise, default numpy
    10.2     PCNR admits chained models              DONE
    10.3(a)  FET limiting, ordinary Newton path      DONE
    10.3(b)  device-level limiter (limit_together)   DONE
    15       PCNR eligibility measured, explain()    DONE (2 of 26 qualify)
             names the refusal
    15       Stage 0: inter-device clash              MEASURED: reproduces
    15       circuit-level write-back forest         MEASURED: dead
    15       vector PCNR Stage 1 (BJT, CPU)          DONE
    15       vector PCNR Stage 2 (FETs, CPU)         DONE -- diff pair
                                                     [14,FAIL] -> [22,22]
    15       vector PCNR Stage 3 (JAX)               PARKED by the user;
                                                     traced path raises
                                                     NotImplementedError
    15       the unlimited-node failure              DONE 2026-08-26 (sec.
                                                     27).  It was the SEED,
                                                     not the node: v_lim_init
                                                     used raw branch
                                                     voltages.  Mirror and
                                                     cascode now converge
                                                     from every start, both
                                                     orders.  A gmin damper
                                                     ships OFF (it trades
                                                     circuits); the EKV pair
                                                     stays on the fallback
    12.1     limiter write-back moves the wrong node DONE
    12.2     check_jacobians "not resolvable"        DONE 2026-08-25
    12.3     no way to declare a device needs gmin   DONE 2026-08-25
    12.4     limit_pnj parameters cannot use var()   DONE
    12.6(c)  limiter params may read the solution   DONE (at last iterate)
    12.4     limit_fet's parameter-only vto         DONE (von, body-biased)
    10.3(b)  BJT adopting limit_together            MEASURED: NO (§12.7) 2026-08-25
    12.6     fourth library batch: findings          2026-08-25 (evening)
    5        Gummel-Poon BJT, EKV core, crystal,
             ferrite, skin-effect R, comparator,
             memristor                               DONE
    5        item 4  opamp macromodel                DONE 2026-08-25
    5        item 10 MOS level 1 (n and p)           DONE 2026-08-25
    5        item 1  self-heating BJT                DONE 2026-08-25
    5        the remaining five: item 6 (VCO), 10
             (Level 3), 12, 13, 15                   DONE 2026-08-26 (sec. 18):
                                                     all fifteen built, 37 classes

Suite: 1874 tests at `4d359e5` before this work, **2351 passed / 7
skipped / 3 xfailed / 0 failed** as of `cb477a8`. After 12.2, 12.3 and
12.4: **2220 passed / 7 skipped / 3 xfailed** in 695 s, plus the two
`test_elements_hdl_library3.py` assertions 12.3 deliberately overturned
and left for that file's owner. `sphinx -W` clean. After
`limit_together` (§10.3): **2239 passed / 7 skipped / 3 xfailed / 0
failed**, the 17 new tests all in `test_device_limiter.py`.

**2026-08-25 (evening), after the fourth library batch (§12.6):
2296 passed / 7 skipped / 3 xfailed / 0 failed**, the 57 new tests all
in `test_elements_hdl_library4.py`.  Run in four chunks (682 + 698 +
831 in `pycircuit/circuit/tests`, 85 elsewhere) because the whole suite
in one process gets killed part-way on this machine.  `sphinx -W`
clean.

**Corrections this campaign made to its own record.** Ten claims were
found right in conclusion and wrong in their stated reason:

- the `Max`-on-an-atom rule (the divide-by-zero symptom does not
  reproduce on sympy 1.14; a `lambdify` that never terminates does);
- `limexp`'s "is it a junction" test (the real condition is PCNR
  eligibility, and the disqualifier was `var()`, not charge — and that
  correction has itself since expired, because §10.2 made `var()` stop
  disqualifying anything. Two revisions in two days on one docstring);
- `_ChainPrinter`'s `reduce` premise (`reduce` measures FASTER than the
  bare call it was said to cost);
- this plan's own section 8 first draft, which measured the eager path
  and generalised it to production models;
- §10.2's "does not admit `DiodeSpiceHdl`" and the reason given for it
  (the default card qualifies; and what refuses `rs = 2` is that a
  resistor's current is not exponential, not a two-branch dependence);
- §12.1's `fetlim` attribution (the write-back moved the wrong node; the
  recorded explanation about a volt being forty thermal voltages was
  plausible and not the cause);
- **§12.2's own proposed fix**, `max(atol, eps·|value|/h)` — ten decades
  too small on the very bias that motivated it, because the quantum is
  set by an internal cancellation the outputs do not show. The
  *diagnosis* on the line above it was exactly right and named that
  cancellation; the formula on the next line used the output anyway;
- **§12.4's stated failure mode**, "not in scope ... a bare `NameError`
  from inside the lambdified parameter function". Measured: the class
  compiles, the parameter function returns the sympy Symbol, and the
  first Newton iteration raises `TypeError: Cannot convert expression to
  float`. `hdl.py` had the correct version of this the whole time, in
  the comment on the source vectors;
- **§12.1's `both = 30`** (2026-08-25). The conclusion — two probes
  competing for a terminal costs iterations — holds at that operating
  point. The number does not: `both` was 13 at the very commit that
  recorded 30, because 30 was copied down from the before-table and
  never re-run, and across 48 operating points `both` is the cheapest
  variant there is. **The cheapest reason to check is a number, and it
  is the easiest to believe** — this one justified building a feature,
  and the feature turned out to be worth building for something else
  entirely (§10.3).

**The pattern is worth naming: a reason decays faster than the
conclusion it supports, and nothing fails when it does.** Every one was
caught by taking a measurement rather than by reading, and most were
one `_src` dump or one `hasattr` away the whole time. Most were
written by whoever was correcting the previous one -- so this is not a
failure of care, it is what happens to
a stated reason when the code under it keeps moving. Line numbers in this document decayed the same way, twice: twelve were
stale within a day of being written because `hdl.py` grew ~1400 lines,
were re-verified, and were all stale again the next day after another
~400. **They have been replaced with symbol names throughout.** A
symbol survives an edit above it; a line number is a measurement of
something nobody is holding still, and re-verifying it is a chore that
comes due again every time the file is touched.

---

## 12. Open items found by building, 2026-08-25

Everything here was found by writing a model against the DSL rather
than by reading it, which is §5's thesis working: eleven library models
across three batches produced more defects than any amount of review
had.

### 12.1 The FET limiter moved nodes a source pins — FIXED 2026-08-25

> **Resolved.** The write-back now chooses its terminal at RUNTIME:
> whichever end of the branch has drifted further from the last accepted
> point. Same circuit, same start:
>
>     limited, compile-time write-back    225 Jacobian evaluations
>     limited, runtime write-back           8
>     unlimited                            25
>
> and the worst displacement of a pinned gate went from **5e48 V** to
> tens of volts. Limiting now helps where it previously hurt 9×.
>
> **The attribution flipped completely**, which is the part worth
> keeping. Before: `fetlim` alone FAILED to converge, and the recorded
> explanation was that SPICE's `fetlim` bounds a step to about a volt
> and a volt is forty thermal voltages to a subthreshold exponential —
> plausible, self-consistent, and not the reason. It failed because its
> write-back always moved the gate. After: `fetlim` alone is the *best*
> of the four (12 iterations against `limvds`'s 57). A seventh instance
> of right-conclusion-wrong-reason, and this one had a measurement
> attached to the wrong cause.
>
> Two smaller findings fell out of it. A limiter that did not bite must
> write NOTHING — `a - (a - b)` is not `b`, so the write alone broke
> "a step of at most 2·VT passes through exactly", which is what lets
> "did limiting fire?" be a convergence signal. And when two probes want
> the same shared terminal, the one applying the LARGER correction
> should get it: a BJT's two junctions both hang off the base and see
> identical drift, so the tie-break decides real behaviour, and ordering
> by row index handed the base to the base-collector probe.
>
> Still open, now measured: `both` (30 iterations) is worse than `fet`
> alone (12), because two probes may not move the same terminal and the
> second is pushed onto a node it would not have chosen. That is what
> §10.3(b)'s device-level limiter removes.
>
> > **THE 30 IS WRONG. It is 13, and it always was** (corrected
> > 2026-08-25 by re-running the table at the commit that introduced
> > it, `de78f42`). 30 is the number from the row ABOVE — the
> > before-the-fix table — copied down into the after-table without
> > being re-measured. `limit()`, `StandardNewton` and the test have
> > not changed since, so the gap being paid for was never 2.5x; it is
> > one iteration.
> >
> > The conclusion survives at this point and only at this point.
> > Across the 48-operating-point grid of `test_device_limiter.py`,
> > `both` is the CHEAPEST of the four per-probe variants — 909
> > Jacobian evaluations against `fet`'s 1222 (which also fails to
> > converge at one point) and `vds`'s 1404. "Adding a second, correct
> > limiter costs iterations" is false as a general statement.
> >
> > And the device-level limiter does not remove the one iteration: it
> > measures 927 over the same grid, never better at a point, one
> > iteration worse at 14 of them. What it removes is a REFUSAL — see
> > §10.3 — not a count. **Ninth instance in this campaign of a claim
> > right in conclusion and wrong in its stated reason; this one's
> > reason was a number nobody re-ran, which is the cheapest kind to
> > check and the easiest to believe.**
>
> And a side effect worth knowing: a stacked pair at (5 V, 2.5 V) that
> used to raise `SingularMatrix` now converges, because the better
> write-back keeps the iterate out of the cutoff region where the lower
> device's conductance underflows to exactly zero. The hazard is NOT
> gone — (40 V, 0.2 V) still reaches it — only harder to reach. §12.3
> stands.
>
> **2026-08-25:** §12.3 landed, and (40 V, 0.2 V) now converges. The
> hazard is still real — the row still goes bit-exactly empty at iterate
> 0 — but the solver now has an anchor for it, and the answer it returns
> is a `gmin = 0` solution, not an anchored one.

#### The original finding, kept for the mechanism

`limit_fet(V(g,s))` writes back as `x[g] = x[s] + vlim`, so it moves
the **gate** — and a gate is almost always driven. `vlim` is bounded;
`x[s]` is not.

On a 40 V cascode started from a uniform 10 V:

    limited      225 Jacobian evaluations
    unlimited     25

and instrumenting `SubCircuit.limit` shows **203 of those 225
iterations move a gate node pinned at 1.0 V by a source**, by as much
as **5e48 V**. Two devices sharing that gate is worse: each writes it
in turn and the second overwrites the first.

This is a defect in §10.3(a) as shipped, not a limitation of it. The
fix is either §10.3(b)'s device-level limiter or, far cheaper, refusing
to move a node that a source pins — and the numbers above are recorded
so a fix can be measured rather than argued.

**It is a caveat, not a verdict.** Measured on the same machinery: a
milliamp forced into a diode-connected EKV converges in 7 Jacobians
with the limiters and raises `SingularMatrix` without them; the 40 V
cascode from the origin is 9 with and 25 without, same answer. The
limiters earn their place. The write-back rule does not.

### 12.2 `check_jacobians` needs a "not resolvable" verdict — FIXED 2026-08-25

> **Resolved, and THE PROPOSED FORMULA WAS WRONG BY TEN DECADES.**
>
> `check_jacobians` now has three per-entry verdicts —
> `JAC_PASS` / `JAC_UNRESOLVED` / `JAC_FAIL` — surfaced as
> `res.verdict` (`'ok'` / `'UNRESOLVED'` / `'FAILED'` /
> `'NOT COMPARABLE'`), as `res.resolved` beside the existing `res.ok`,
> and as `res.unresolved` / `res.failures` lists of `JacEntry`
> records carrying the analytic value, the difference, the error, the
> floor and the reason. `res.ok` stays true when entries are merely
> unresolved, so no existing assertion changes meaning; `res.resolved`
> is the stronger claim and has to be asked for.
>
> **What the item got wrong.** It prescribed a floor of
> `max(atol, eps·|value|/h)`, and both models' comments proposed the
> same. Measured on EKV at `vgs = -3 V`, the bias that motivated it:
>
>     eps*|q|/h        (the proposal)      9.4e-36
>     eps*max|q|/h     (vector form)       9.0e-26
>     the real quantum / 2h                2.1e-24
>     the entry it must cover              1.6e-25
>
> The proposal is **ten decades too small** and even the vector form
> misses by 2×, because the cancellation that sets `q`'s representable
> step happens on an **internal magnitude of 1.9e-15** that no output
> shows. A floor derived from the outputs cannot see it. So the
> roundoff floor is **measured**: when a value comes back bitwise
> identical at all four probe points, the step is widened by a decade
> at a time until the value clears its own quantisation, and the floor
> is that quantum over the widened step. On the EKV point that turns a
> difference of *exactly 0.0* into `-1.26e-25` against an analytic
> `-1.28e-25` — 1.5% corroboration where there had been none.
>
> **Widening is a discriminator, not a whitewash.** A value genuinely
> independent of `x[k]` stays frozen at every rung, and a non-zero
> analytic entry against it still FAILS. That case has its own test.
>
> The other two floors are measured too. **Truncation** comes from
> Richardson — the same column differenced at `h` and `2h`, whose
> disagreement is three times the `h²` term. **A kink is detected from
> the VALUE alone**: the one-sided disagreement
> `|f(x+h) - 2f(x) + f(x-h)| / h` halves with `h` for a smooth
> function and stays put at a corner. That detector never looks at the
> Jacobian, which is what makes it a fact about the model.
>
> **Is a kink distinguishable from a real error? Partly, and the limit
> is measured.** The non-smoothness is distinguishable — it is a
> property of `f`. What is *not* distinguishable is an error smaller
> than the jump, because at a corner the derivative is genuinely
> two-valued. The floor is 0.75 of the jump, so on `i = g·max(v,0)` at
> `v = 0` an analytic entry anywhere in `[-0.25g, 1.25g]` passes and
> anything outside fails. That band is asserted from both sides rather
> than described.
>
> **A fourth mechanism turned up that the item did not list**: a step
> so large that the Taylor expansion is not converging at all. On the
> memristor's `ron = 1, roff = 1e9` card at `x = 1`, `R` moves by 100×
> across one 1e-7 step and the difference reports 3.5e6 against a
> correct 7.0e8. No finite floor covers that honestly, so when
> Richardson's term is a quarter of the difference itself the entry is
> UNRESOLVED with an unbounded floor.
>
> **The cost of that rule, measured, as a count.** Run the memristor's
> whole 162-point sweep against a deliberately HALVED `G`: **156
> caught, 6 missed**, and all six are the stiff card at `x = 1`. Not
> because the rule is loose — because a 1e-7 central difference at that
> point carries no information about the derivative, and because at a
> corner whose forward arm is flat a halved Jacobian coincides exactly
> with the average the difference returns. The alternative is a FALSE
> FAILED on a correct model, which is what 12.2 exists to remove.
>
> **A coincidence worth knowing before writing a test here.** At a
> corner between a slope of `g` and a slope of `0`, the central
> difference sits at `g/2` — so a *halved* analytic Jacobian lands
> exactly on it. `assert 0.5*g != g/2` is not a test anybody can pass.
> The same coincidence appears at the memristor's stiff clamp corner.
> Corrupt by three, not by a half, when the point is a kink.
>
> **Workarounds removed.** EKV's two hand-written `atol = 1e-24`
> overrides are gone. The memristor's finite-difference sweep gained
> the `ron = 1, roff = 1e9` card (4 real FAILURES before, 0 after), the
> clamp corners `x = 0` and `x = 1`, and `|v| = 1e4` — 162 points where
> there were 56, of which **122 are fully resolved** and the other 40
> say which mechanism stopped them.
>
> **The cost**, since it is a diagnostic that library tests call
> hundreds of times: four evaluations per column instead of two, plus
> the ladder where it fires. Measured — BJT (n = 6) 3.4 ms → **6.4 ms**
> per call, EKV (n = 4) 0.8 ms → **2.1 ms**, EKV in deep cutoff where
> the ladder runs **2.7 ms**. The memristor test now covers 162 points
> where it covered 56 and the file still runs in ~20 s.
>
> **The floors do not widen a resolvable entry**, which is the property
> that keeps the instrument an instrument. On a diode at 0.42 V the
> band is 4.4e-12 and the floors are 2.5e-15 (roundoff) and 7.1e-15
> (truncation) — three decades under. `test_check_jacobians.py`
> asserts that, and runs every deliberately-built mechanism a second
> time with a corrupted derivative.

#### The original finding, kept for the mechanism

Two independent models, three distinct mechanisms, all reporting FAILED
where the model is right:

- a genuine **kink** (a window seam, a clamp corner) — central
  differencing across a jump, and no `h` helps because a jump has no
  scale;
- an entry whose own signal is **below one ulp of the value** it
  differentiates, so the FD returns exactly `0.0` (EKV's normalised
  charges at deep cutoff: `C` entries ~1e-25 F against a band of
  1e-32);
- **truncation** on a stiff card where `dR/dx = 1e9`.

There is no single `h` that serves all three — measured: `h = 1e-9`
fixes the third and breaks the default card.

**The fix, now well evidenced: report a per-entry noise floor
`max(atol, eps·|value|/h)` and mark entries below it UNRESOLVED rather
than FAILED.** Both models currently carry hand-written `atol`
overrides with the diagnosis in a comment beside them, which is the
workaround this would remove.

> The conclusion held; **the formula in that last paragraph is wrong by
> ten decades** and is corrected above. An eighth
> right-in-conclusion-wrong-in-its-reason, and this one was written by
> whoever had just diagnosed the bias correctly — the diagnosis
> ("computed as a difference of order-0.5 quantities that cancels") is
> exactly right and names an internal magnitude, and the formula
> proposed on the very next line uses the OUTPUT.

### 12.3 No way to say "this device needs a GMIN anchor" — FIXED 2026-08-25

> **Resolved, and the item's own framing was wrong in two places.**
>
> `DC` now takes a `gmin` parameter (default **1e-12 S**, SPICE's
> `GMIN`, `gmin=0` to disable) and `nrsolver.GminAnchorNewton` sits
> **outermost** in the DC chain. What it engages on is the correction
> that mattered: every layer below it engages on `NoConvergenceError`,
> and a singularity passes straight through all four of them by
> §stage-6(a) design, which is exactly why an empty row had no rescue at
> all. The anchor takes the singularity, and only that.
>
> **The item said "structurally singular". It is not.** The stacked
> pair is singular *at iterate 0*, where both devices are in cutoff and
> the `nm` column is bit-exactly zero, and perfectly well posed *at its
> answer*, where the same conductance is **3.8e-8 S** — four decades
> above the anchor that got it there. That distinction is the whole
> item: **gmin must rescue a numerically-empty row and must not disguise
> a structurally-empty one**, and the two are the same matrix shape.
>
> Two gates tell them apart, both asked about the `gmin = 0` system
> because that is the only system the caller asked about:
>
> 1. `_structural_singularity` on the pure Jacobian **at the converged
>    point**. Still-empty column → reject and name it. Catches the
>    capacitor-only node and the current-source-only node.
> 2. **the answer must not depend on the anchor** — re-solve at
>    `gmin/10` and require agreement to the solver's own `conv_x`. A
>    quantity gmin is holding up moves a decade when gmin does.
>
> **And the answer normally owes gmin nothing at all**, because the
> ladder finishes with a plain `gmin = 0` solve from the anchored seed.
> On every circuit measured for this item that solve converges, so gmin
> is a homotopy path and never a term:
>
>     gmin        v(nm), 40 V stacked pair    anchor retained?
>     1e-13       0.01805186588               no
>     1e-12       0.01805186588               no
>     1e-11       0.01805186588               no
>     1e-10       0.01805186588               no
>
> Four decades of gmin, 4e-14 V of movement. And it is the *right*
> answer, checked by taking the anchor away in the limit:
>
>     explicit leak from nm to gnd   v(nm)             rel. to anchored
>     1e9  Ohm                       0.017588037157    2.57e-2
>     1e12 Ohm                       0.018051385528    2.66e-5
>     1e14 Ohm                       0.018051861024    2.69e-7
>     1e16 Ohm                       0.018051865779    5.60e-9
>
> — which also says the **1 GOhm testbench leak the library tests use is
> 2.6% off**, i.e. the anchor is the more faithful of the two.
>
> **The model's refusal was right, and now has a number.** A *permanent*
> 1 pS to ground — what a naive always-on gmin would be — moves the EKV
> drain current by **15% at `vgs = 0`**, 4.4% at 0.05 V, 1.3% at 0.10 V.
> The shipped anchor moves it by **nothing**: `np.array_equal` on the
> whole solution vector, because on a solve that succeeds the anchor is
> not in the matrix at any iterate.
>
> **Cost.** 24 Jacobian assemblies for the rescued 40 V pair against 5
> before the unanchored solve gives up. Zero on any circuit that already
> converged.
>
> **Gate 2 has teeth on a real circuit, measured by sabotage.** Delete
> the final `gmin = 0` solve so the anchor becomes permanent, and gate 2
> *rejects the 40 V stacked pair itself* — the `vdd` branch current
> moves from -5.688e-10 A to -4.968e-10 A when gmin drops a decade,
> **68x** its own convergence tolerance. At that bias the supply current
> is leakage of the same order as the anchor, so "1 pS is negligible" is
> false there and the gate says so. Two other sabotages were run and
> caught: anchoring branch rows as well as node rows (1 test), and
> setting the default `gmin` to 0 (6 of 12 tests).
>
> **One incidental bug, found by the gate's own message.**
> `analysis.reduced_row_names` built branch names with
> `getattr(b, 'name', i)`. The attribute EXISTS on a `Branch` and is
> usually `None`, so the default never applied and every branch row in
> every diagnostic in this tree was called `'branch None'`. Fixed to
> `getattr(b, 'name', None) or i`.
>
> **A gate that reads well and is fooled, recorded because it was the
> first design.** Gate 2 was going to be "does the anchored answer
> satisfy the unanchored KCL to `reltol * I_scale + abstol`?" —
> `conv_f` on the pure system. `I_scale` is `|J|.|x|`, so an anchor that
> drives a node to 5e8 V inflates the tolerance by the same factor it
> inflated the answer: the 1 mA island below misses by 5e-4 A against a
> tolerance of **100 A** and sails through. *A tolerance computed from
> the corrupted answer cannot detect the corruption.* Pinned as a test.
>
> **Where the distinction is imperfect, stated plainly.** The verdict is
> a property of the ANSWER, not of the netlist. A node whose only DC
> path runs through a device that is still exactly off at the answer is
> reported structural — defensible ("nothing determines it") but it
> means the same netlist can be rescued at one bias and refused at
> another. And one shape is accepted with the anchor **retained**: an
> island with no injection has infinitely many solutions, every common
> mode satisfies KCL exactly, and gmin picks the minimum-norm one — what
> every SPICE does. `DC.gmin_anchor_retained` and a logged warning are
> how a caller finds out; it is the only circuit found in this campaign
> that sets it.
>
> **Scope: DC only.** The transient's inner Newton has no anchor, and
> `doc/transient_review.md`'s "no gmin stamping in the transient loop"
> row still stands.
>
> **Two behaviour changes in a file this work does not own**, both of
> them the recorded claim this item exists to overturn. Left alone, as
> instructed; each needs its singular assertion inverting or pinning at
> `gmin=0` by whoever owns `test_elements_hdl_library3.py`.
>
> - `test_a_stacked_pair_without_a_dc_path_is_structurally_singular`
>   asserts `pytest.raises(SingularMatrix)` on `DC(_pair(40, 0.2))`.
>   That call now converges to `v(nm) = 0.0180519`. The rest of the test
>   — the 5 V bias, the leak-resistor insensitivity — is unaffected.
> - `test_fet_limiting_rescues_a_solve_that_otherwise_goes_singular`
>   asserts `isinstance(res2, SingularMatrix)` for the milliamp-forced
>   diode-connected EKV **with its `$limit` declarations removed**. The
>   anchor rescues it — and this is the more interesting of the two,
>   because it lands on **the limited model's own answer**:
>
>       limited     v(nd) = 0.3069549236
>       unlimited   v(nd) = 0.3069549236   (anchor not retained)
>       relative difference                 5.7e-11
>
>   So the limiters buy a shorter path and not a different answer, and
>   the anchor is an independent route to the same point. That
>   measurement is *not* pinned by a test in this work: reproducing it
>   needs `_ekv_no_limit`, which strips `$limit` from the analog block
>   and lives in the file this work must not edit. `test_gmin_anchor.py`
>   deliberately does not copy it — a first attempt did, using a
>   `limit_spec = ()` subclass, and that attribute does not exist on
>   `EkvNmosHdl`, so the "unlimited" model was the limited one and the
>   test compared a number with itself and passed.

#### The original finding, kept for the mechanism

EKV's channel conductance underflows to *exactly* zero in deep cutoff —
`softplus(-800)` is `0.0`, not a denormal — so a stacked pair with no
other DC path to the intermediate node is structurally singular.
`compact.PspMosLongChannel` carries a private `GLEAK = 1e-12` for
precisely this, and the library model deliberately does **not**, because
1 pS is a picoamp at a volt and would sit on top of the weak-inversion
currents the model exists to measure.

Both facts are tests. This belongs at the simulator level as `gmin` (or
`gmin` stepping), not in every model that has an off region.

### 12.4 Smaller, each with its cost known

- ~~**`limit_pnj`'s parameters cannot be written with `var()`.**~~
  **FIXED 2026-08-25 — the first route, resolution.** A limiter
  parameter that mentions a `var()` symbol is now compiled through
  `_chain_compile` against the let-chain, the same path the PCNR
  detector's `VT` and `IS` take, so the chain is **read, not inlined**
  — a parameter written over a 40-deep chain costs what the chain
  costs, and there is a test that says so. A parameter that mentions
  nothing but the card still goes through `sympy.lambdify` exactly as
  before, so nothing on the common path moved.

  Resolution turned out to be cleaner than refusal *because it forces
  the question refusal would have dodged*: what a limiter parameter may
  legitimately be. It is evaluated from the card, once, BEFORE the
  device is evaluated — so a chain that reaches the iterate is now
  refused by name at compile time ("`limit_pnj`'s parameters ... cannot
  depend on the solution; this one reaches plus (through
  var('isbad'))"), which is a class of wrong answer that had no error
  attached to it at all before. Time-dependent chains are refused the
  same way.

  **And the failure mode this item states is wrong** — a ninth
  right-in-conclusion-wrong-in-its-reason. It says the symbol "is not
  in scope", and `elements_hdl.py`'s comment beside the duplicate says
  "the compile fails with a bare `NameError` from inside the lambdified
  parameter function". Measured against the unfixed code: the class
  compiles, the instance constructs, and the parameter function
  *returns the sympy Symbol* — `[_v0_isT, 0.02584269662921349]`. The
  failure is `TypeError: Cannot convert expression to float` from
  inside `limit()` on the first Newton iteration. `hdl.py` already had
  this right, a hundred lines away, in the comment on the source
  vectors: "lambdify handed one it cannot resolve produces a function
  that quietly evaluates to nothing useful rather than complaining".

  **`elements_hdl.py`'s Gummel-Poon BJT can now drop `_isr`** and write
  `limit_pnj(bbe.V, isT, nf * _vt(T))`. Verified by rebuilding the
  model with the substitution applied: the limiter parameters are
  bit-identical (9.772967086868698e-17, 0.02584269662921349) and
  `limit()` agrees to 0.0 over 200 random (x, x0) pairs. Not done here
  — `elements_hdl.py` was not this change's to edit.

  This also **expired the stated reason in `limit_fet`'s docstring**,
  which forbade a bias-dependent `vto` because "the parameter
  expressions are lambdified over `paramsyms + [TEMP]` and nothing
  else". The rule is unchanged and now enforced; its reason is *when*
  the limiter runs, not what its namespace contains.
- **`limit_fet`'s parameter-only `vto` bites, and the size is known**
  *(CLOSED 2026-08-25, §12.7 -- `von` is now bias-dependent and the
  565 mV is gone)*:
  for the EKV card at 2 V of body bias the true turn-on is 1.06 V
  against `vto = 0.50`, so every clamp sits **565 mV** low. Looser, not
  wrong — the step is still bounded and the no-op band intact.
- **`0 × d(√0)` is reachable from a legitimate CARD, not just a wild
  bias.** EKV's pinch-off root with `gamma = 0` — a perfectly good "no
  body effect" card, and the one every textbook-asymptote test uses —
  sits exactly at zero with an infinite derivative. This is the
  `differentiable-numerics` rule in a form that skill does not spell
  out: the dangerous value was a *parameter*, not a bias.
- **A per-probe write-back can UNDO another probe's branch, and the
  size is known (2026-08-25).** "No probe is undone" was verified at
  the level of NODES — no row is written twice — and that is not the
  same claim. Two probes sharing a terminal: the one with the larger
  correction takes the shared node, and the other probe's branch, which
  hangs off it, follows and lands somewhere its own law never chose.
  Measured over 813 random steps in which both probes of the
  `(vgs, vds)` star bite: **27 of them (3.3%) leave a branch its own
  law would still move, the worst by 36 V.** `limit_together` has no
  such case (0 of 813) because it solves branch CONSTRAINTS rather than
  applying node DISPLACEMENTS. Not fixed in the per-probe form: it
  would mean re-solving after each write, which is the device-level
  form with extra steps, and the per-probe form measures fine on the
  grid regardless (909 against 927). Recorded so the next person does
  not read the test name as the stronger claim.
- **`hdl.var` vs SPICE's `VAR`.** `elements_hdl.py` aliases the
  let-chain binder as `_var`, and a comment now says that is "a
  readability choice, not a workaround". It is a workaround again: the
  Gummel-Poon BJT could not otherwise use SPICE's own parameter name
  for the reverse Early voltage. §3's S5 fixed the leak, not the
  collision.

### 12.5 What the kernel work bought, measured across three batches

    batch 1   three hand-rolled clamps, each with a floor-placement
              argument and a smoothing-width scan
    batch 2   zero hand-rolled clamps
    batch 3   zero; `safe_pow` replaced the diode's three-line floored
              power with one call, and `softplus` is the whole of EKV's
              all-region behaviour (the literal form overflows at
              x = 710 and returns 710)

That is §3's kernel additions paying for themselves in the only
currency that counts here: what an author has to think about that is
not physics.


### 12.6 The fourth library batch, 2026-08-25 (evening)

Three models -- `MosLevel1Hdl`/`MosLevel1PmosHdl`, `OpAmpHdl`,
`GummelPoonNpnThermalHdl` -- and 57 tests in
`test_elements_hdl_library4.py`. As before, the models found things
review had not.

#### What worked, said first because it is the unusual half

Four pieces of machinery met a real model for the first time in this
batch and **all four behaved**, which has not been true of any previous
batch:

- **`limit_together` is exactly what its docstring says it is.**
  `MosLevel1Hdl` declares SPICE's own four probes -- `fetlim(vgs)`,
  `limvds(vds)`, `pnjlim(vbs)`, `pnjlim(vbd)` -- and the per-probe form
  refuses the production model with the message that names the fix:
  *"the $limit probes over-determine this device: both terminals of
  branch (b,di) have already been moved by earlier probes ... Declare
  them with limit_together()"*. Nothing had to be worked around.
- **`param_given` with `$limit` AND `Collapse` in one model.** First
  time all three have been in one element. Works. (The §12.4 fix to
  `limit_pnj`'s parameters, and `test_limit_with_param_given.py`, are
  what made the first two combinable; the third had never been tried
  with either.)
- **`Collapse` takes a CONJUNCTION.** `Collapse(brd, And(rd <= 0,
  rsh*nrd <= 0))` is what SPICE's "`rd` if given, else `rsh` times the
  squares" needs, and it compiles and caches per variant like any other
  condition.
- **`check_jacobians`' UNRESOLVED verdict fires where it should and is
  silent where it should not.** Ten points on the MOSFET: seven fully
  resolved, three not, and both mechanisms are properties of the model
  rather than of the checker --
  *roundoff* on a reverse-biased junction conductance of ~1e-46 S
  against a measured floor of 4.4e-23, and *truncation* at `vgs = vth`
  exactly, where `(beta/2)*max(vgs - vth, 0)^2` is C1 but not C2 so a
  central difference returns `beta*h/4` (4.4e-11 S against an analytic
  3.9e-13 S). The control -- the same element with `G` and with `C`
  multiplied by **three**, not by a half, per §12.2's corner
  coincidence -- FAILS in both cases.
- **the gmin anchor neither helped nor surprised.** It changed no
  answer here (four decades of `gmin` move the forced-junction solve by
  literally zero), and it twice REFUSED to manufacture one: for an
  output current demanded beyond the opamp's `isc`, and for a
  thermally-runaway BJT above its critical `rth`. Both refusals name
  the node. §12.3's gate is doing the job it was built for.
- **`explain()` prints the group** (`fetlim on (g,s) [group 1]`, four
  times) and the compile-time messages were right every time they
  fired, including the `NameError` for an un-imported `limit_together`,
  which listed the class's declared parameters.
- **`param_given` inside `var()`** works, which is not obvious given
  that state operators inside `var()` were unreachable until Phase D
  and `$limit` markers were being stripped from `var()` definitions as
  recently as this campaign. Four of them are in `MosLevel1Hdl`'s
  chain.
- **the Gummel-Poon body moved out of its closure without moving a
  number.** `_gummel_poon`'s nested `analog` became a module-level
  `_gp_core` so the self-heating variant could share it, and
  `test_elements_hdl_library3.py`'s 67 tests -- which include an
  independent numpy transcription of the whole model and an RK4
  integration of a common-emitter stage -- passed unchanged on the
  first run. That is what a test file written against the physics
  rather than the implementation buys.

#### The friction, in the order it was met

**(a) A helper shared by two classes cannot see its parameters, and the *(CLOSED 2026-08-26, §17: `params_as` hands the helper a namespace; `_with_params` is gone.)*
BJT's card has forty-two of them.** `hdl._analog_function` injects the
parameter symbols into `analog()`'s OWN globals, so a module-level
helper it calls has a different `__globals__` and the bare names do not
resolve. `_spice_diode` pays for this by taking all nineteen of its
parameters as explicit arguments; forty-two was not tolerable, so
`elements_hdl` now has

    def _with_params(func, ns):
        return types.FunctionType(func.__code__, ns, func.__name__,
                                  func.__defaults__, func.__closure__)

    ## at the call site, inside analog():
    _with_params(_gp_core, globals())(heat.T, +1, c, b, e)

which is `_analog_function`'s own constructor applied by hand. **What
should have been written is `def analog(c, b, e, th, tha, p)` with
`p.bf`** -- roadmap S5's parameter namespace, which §11 records as DONE
because the *leak* was fixed. The leak was; the ergonomic half was not,
and this is what it costs. The workaround also does not compose: it
only works from inside the defining module, because the namespace it
copies has to be `elements_hdl`'s.

**(b) Three SPICE MOSFET card names cannot be Python identifiers.** *(CLOSED 2026-08-26, §17: `p['lambda']` in a `params_as` model, and `aliasparams` on Level 1.)*
`LAMBDA` and `AS` are keywords and `IS` is one letter from being one.
They ship here as `lambd`, `asrc` and `IS`. `globals()['lambda']`
would in fact resolve inside `analog()` and `cls(**{'lambda': 0.02})`
would set it, so the name is reachable -- but neither spelling is one
anybody should write. A parameter namespace (item (a)) fixes this too:
`p['lambda']` is legal Python.

**(c) A self-heating device's limiter parameters reach the solution,
and the DSL is right to refuse them.** *(CLOSED 2026-08-25 -- see
§12.7: the refusal was right about the order and wrong about the
conclusion, and is gone.)* `limit_pnj(bbe.V, isT, ...)` is
the spelling §12.4 went to some trouble to make legal, and with a
thermal node `isT` depends on `V(th,tha)`. `generate_code` refuses by
name -- *"cannot depend on the solution; this one reaches th (through
var('dT'))"* -- which is correct and was a five-minute fix rather than
a bisection. The model now builds a SECOND, ambient-temperature copy of
the saturation current for the limiter alone, under a build-time
`if tlim is None`, so the isothermal element is untouched. The cost is
five duplicated lines and a limiter placed against a junction that may
be 100 K hotter than it thinks -- looser, not wrong, exactly as
`limit_fet`'s parameter-only `vto` is. **What would remove it:** a
limiter parameter evaluated at the LAST ACCEPTED iterate rather than
from the card. That is information a limiter already has (it is handed
`x0`), and it is the natural home for every "SPICE recomputes `von`
each iteration" case, including `limit_fet`'s.

**(d) `SelfHeating`'s thermal pin is a PARALLEL path, not a series
one.** The docstring says the pin "can carry an external Foster or
Cauer ladder", and it can -- but the device's own `rth` stays across
the same `(th, tha)` branch, so the ladder is in parallel with it and
the obvious way to say "no internal thermal resistance, use my ladder"
is exactly the trap the same docstring warns about: `rth = 0`
collapses the branch to a zero-volt source and SHORTS the ladder out.
The only route today is a large `rth`, whose value then leaks into the
answer -- measured, `rth = 1e6` against a 5000 K/W ladder still puts
**0.5%** of the heat through the device's own path (the parallel value
is 4975 K/W, not 5000). **The fix is small and is in
the model rather than the DSL:** a `SelfHeating(..., external=True)`
that simply does not emit the RC contribution, chosen by a build-time
Python `if`, which is legal. Not done here; the test
`test_the_thermal_pin_takes_an_external_foster_ladder` asserts the
parallel behaviour against the exact parallel value so that the
behaviour is pinned either way.

**(e) `VS`'s `vac` defaults to 1.** Not an HDL matter, and it cost an
hour: every node in an opamp AC testbench came back at exactly 1.0 and
looked like a singular matrix, because the supplies and the
non-inverting input were all driving the AC analysis too. Every AC
testbench in `test_elements_hdl_library4.py` now sets `vac` explicitly
on every source and says why.

**(f) `expl` is not finite for every double.** `expl`'s continuation is
a third-order Taylor, so its cubic term overflows: bisected on the
compiled function, the last finite argument is **4.76e69**, where it
returns 1.798e308. On the MOSFET that is a bulk-node voltage between
1e68 V (finite) and 1e69 V (not), so the gap is academic -- but
`expl`'s own docstring and `differentiable-numerics` both said "finite
for every double", and a guarantee that is nearly true is the kind that
gets relied on at the one bias where it is not.

> **CLOSED 2026-08-25, and this entry had the same defect it reports.**
> Both places were corrected in the commit that landed this batch. And
> the number above was first written here as "about 2.6e70", which was
> an ESTIMATE from reading the Taylor form and was wrong by 5.5x; the
> 4.76e69 is bisected on the function. An unmeasured number in a
> friction log about unmeasured numbers -- eleventh instance in this
> campaign of the shape §11 catalogues, and the cheapest kind to check.

(The opamp's own boundary is `hypsmooth`'s: the Jacobian goes at
1e107 V on a supply pin and the value at 1e160, both from a radicand
that squares a quantity 1e4 larger than the argument.)

#### The mutation check, and the two that got through first time

Twenty-two deliberate corruptions of `elements_hdl.py`, each applied,
run against the batch's own tests, and reverted. **Twenty were caught
immediately. Two were not, and both were failures of the TEST rather
than of the corruption:**

- **the drain overlap capacitance moved to the source branch** --
  `cgdo` replaced by `cgso` -- and the card had `cgso == cgdo`, so the
  "corrupted" model was bit-identical to the correct one. This is
  §12.2's warning in its other form: *check that a deliberately
  corrupted fixture is actually corrupt at the point you test it.* Fixed
  by giving the three overlap densities three different values, which
  is better validation design anyway (`validation-design`: "a zero
  coefficient does not test a scaling" -- an EQUAL coefficient does not
  test a routing);
- **the `-IBC` term deleted from the self-heating BJT's dissipation.**
  Every bias in the test was forward active, where the reverse
  transport current is 1e-16 A and the term is 1e-13 of the answer.
  Fixed by adding two SATURATED biases (`Vce = 0.2 V` and `0.05 V`),
  where the measured beta falls to a twentieth of `bf` and the term is
  17% of the collector current.

A third, added after the first pass, also got through and is recorded
for the same reason: **removing the factor of 2 from SPICE's
forward-bulk threshold continuation** changed nothing, because every
bias point in the reference comparison had `Vbs <= 0` and that arm is
never selected there. Two forward-bulk points were added and it is a
6.4% error.

All three are the same shape and it is the shape `validation-design`
names: **a bias no sweep visits tests nothing, and neither does a
coefficient no card distinguishes.** Corruptions were by a factor of
three rather than a half wherever a corner was involved, per §12.2.

#### The measurement that was worth taking

**On a level-1 MOSFET, `fetlim` and `limvds` earn almost nothing and
`pnjlim` earns everything.** Over a 48-point cascode grid with `gmin`
off and no rescue ladder, limited and unlimited cost **552 and 554**
Jacobian evaluations -- a difference of two. The reason is structural
and worth stating because it does not generalise: **a level-1 channel
is a POLYNOMIAL.** There is no subthreshold exponential for `fetlim` to
compress, which is why `EkvNmosHdl` -- whose channel is `softplus^2`,
exponential below threshold -- is where those two limiters were
measured to be worth having.

What is exponential in a level-1 MOSFET is the two bulk junctions, and
there the four-probe group is decisive: forcing a current into the bulk
node, over four decades from 1 uA to 1 A, **plain Newton without the
limiters does not converge at any of them**, and with them it converges
in **at most three** Jacobian evaluations, to the same answer the full
`DC` rescue chain finds for the unlimited element (1e-13 relative).

So §10.3's conclusion holds and its reason is confirmed on a production
model: what `limit_together` buys is a CAPABILITY -- the four-probe
declaration exists at all -- and the iteration count it buys depends
entirely on which probes the device's physics actually needs.

#### What the opamp found about the claim that motivated it

§5 item 4 says an HDL opamp "deletes `macromodels.OpAmp`'s
runtime-`import jax` lambdas". **Verified, and the claim is true, still
live, and the smaller half of the difference.**

- `macromodels.OpAmp` still carries `import jax.numpy as jnp` INSIDE
  the closure it hands to `BSource`, so it runs on every evaluation and
  on every finite-difference probe of one. The *element* half of that
  pattern was fixed by `toolkit.derivative`; the *user callable* half
  was not, because it is user code. 0.43 us per call, measured.
- `BSource`'s Jacobian is `toolkit.derivative`, a central difference at
  `eps = 1e-6`. Accurate at this bias (1.2e-10 relative on a `tanh`),
  so this is a structural point -- two extra evaluations per entry, and
  nothing a tracing backend can differentiate -- not an accuracy one.
- **The accuracy difference is elsewhere, and it is large.**
  `macromodels.OpAmp` limits its output with `Vspan*tanh(v/Vspan)`,
  whose gain error is `(v/Vspan)^2/3`: **0.15% at 1 V out of a 15 V
  swing and 12.6% at 8.7 V**, which is why its own test asserts the DC
  gain only to 2e-3. `OpAmpHdl`'s `hypsmooth` pair, with a smoothing
  width written relative to the swing, is **2e-8** at the same points.
- **It is not a drop-in.** `macromodels.OpAmp` is a four-terminal
  `SubCircuit` with a differential output and no supply pins;
  `OpAmpHdl` has five terminals, two of which are the supplies it takes
  its rails and its PSRR from. Anything wanting a floating output
  reference should keep the old one, and `macromodels.OpAmp` is
  therefore left in place.

`OpAmpHdl` also carries what the old one has no expression for: slew
-rate limiting, output current limiting, input offset voltage and bias
and offset currents, CMRR and PSRR, and rails taken from real supply
pins rather than from parameters.

**One structural consequence of having no ground pin, recorded because
a user will meet it.** The output is referred to the supply MIDPOINT --
an element can only reach the nodes in its own signature, so "referred
to the global ground" is not expressible from inside `analog()`, which
is the same limitation `SelfHeating` records for its ambient pin. So
moving one rail alone moves the reference and the output follows it: a
naive PSRR measurement to ground reads 1.5 V per volt of `vcc` where
the differential-supply measurement reads 2.0. Both numbers are in the
test.

---

## 13. Re-measured 2026-08-25 (evening): §8's numbers moved, and so did its conclusion

§8 was measured before the kernel additions of §3 landed. Re-measured on
the same machine, same card, same bias:

    quantity                     2026-08-24      now      change
    PSP compile+instantiate         67.4 s      66.4 s     --
    i                               1414 us     1282 us    -9%
    G                              26523 us    17304 us    -35%
    G generated lines               2675        2675       --
    G dispatch calls               14030        3471       -75%
    G/i                             18.8x       13.5x

**The kernel work paid for itself in a currency it was not aimed at.**
`maxc`/`minc` were added for *ergonomics* -- to remove the undocumented
"`Max` only on an atom" rule -- and rebuilding `expl`/`hypsmooth` on
them cut the Jacobian's dispatched calls by three quarters and its
runtime by a third. Nothing in §3 predicted that, and §8's own
"op counts are the wrong proxy" caution is why: the win came from
emitting FEWER CALLS, not from a smaller expression. The line count did
not move at all.

### What this does to the backend decision

Priced at the measured scalar costs, dispatch is now roughly

    2642 where x 0.94 + 454 min/max x 0.60 + 303 exp x 0.10  =  ~2.8 ms

of 17.3 ms — **about 16%, where §8 measured 40%.** So:

- **scalar code generation is now worth ~1.2x, not ~1.7x.** It has
  been overtaken by a change made for other reasons, and on these
  numbers it is no longer worth doing on its own.
- **the remaining ~84% is the Python interpreter**, not numpy: 2675
  lines of straight-line float arithmetic plus 2775 `_recip2`/`_rdiv`
  Python-level calls. That is not a dispatch problem any printer change
  can reach.

**This strengthens the compiled-backend case and weakens every
alternative to it.** The cost is now interpretation itself, and only a
backend that stops interpreting removes it. `_chain_compile` and
`_ChainPrinter` remain the replacement point.

**And it is a warning about §8's own remaining numbers.** The 300-1000x
estimate for a compiled backend was extrapolated from 14030 operations
at 1-5 ns each. The operation count is now 3471 dispatched calls over
2675 lines, so that extrapolation needs redoing before anyone commits
to it. **Do not quote §8's ratios without re-running
`benchmarks/hdl_model_cost.py`** -- this section exists because they
were quoted once already, four weeks stale, in the space of one day.

---

## 14. The backend decision, and the work that answers it

§8 diagnosed the Jacobian's cost and §13 re-measured it. Neither
decides anything. This section states the decision, and specifies the
work that would turn it from an argument into a measurement.

### 14.1 What is being decided

**Whether this DSL is a research instrument or something a real circuit
can be put through.**

PSP's Jacobian costs **17.3 ms**. The arithmetic that matters is not
that number but what it becomes:

    100 devices, one Newton iteration            1.7 s of Jacobian
    1000-step transient, 5 iterations per step   ~2.4 hours
    a compiled C compact model, same job         ~10 us per device

Everything else in this plan improves the DSL for the person *writing*
a model. This is the only item that changes who can *use* one.

### 14.2 Where the cost is, measured today

    G = 17304 us over 2675 generated lines
      dispatched numpy calls   3471   ~2.8 ms   ~16%
      everything else                 ~14.5 ms  ~84%

The 84% is **the Python interpreter**: straight-line float arithmetic
over 2675 lines, plus 2775 `_recip2`/`_rdiv` Python-level calls. No
printer change reaches it.

### 14.3 The options, priced

| option | worth | effort | verdict |
|---|---|---|---|
| do nothing | — | — | legitimate if 17 ms blocks nothing you run |
| scalar codegen (`math`, `if`/`else`) | **~1.2x** | days | **dead as a standalone item** — §8 priced it at 1.7x and the kernel work took the rest |
| compiled backend (C / LLVM / numba) | unknown | weeks+ | the only option that reaches parity |

**The cheap option died while nobody was looking**, and that is the most
useful thing §13 found. It was overtaken by a change made for
ergonomics. It survives only as a component of the third option.

### 14.4 The number this plan is NOT entitled to quote

§8 said a compiled backend would be "300-1000x", extrapolated from
**14030 operations at 1-5 ns each**. That count is now **3471**. The
extrapolation has not been redone and **must not be repeated until it
is**. This campaign has produced ten claims that were right in
conclusion and wrong in their stated reason, and several were exactly
this: a figure that was true when taken, quoted after the ground moved.

### 14.5 The work: a spike, not a build

**Do not choose between C, LLVM and numba by argument.** Build the
smallest thing that measures each, on the real model, and let the
numbers choose. Estimated one day.

**Input.** `PspMosLongChannel`'s compiled `G` chain — 2675 lines of
straight-line scalar float code with `where`/`min`/`max`/`exp` and two
helper calls. It is already in the shape every candidate wants;
`_chain_compile` and `_ChainPrinter` are the replacement point.

**Candidates, cheapest first:**

1. **numba** on the existing generated source. No new code generator at
   all — the chained path already emits its ideal input. Measure
   compile time as well as run time: a 2675-line function may cost more
   to JIT than a device model can amortise.
2. **C via cffi**, emitted from the same chain. Most control, most
   work, no runtime dependency beyond a compiler.
3. **LLVM via llvmlite**. No external toolchain; between the two on
   effort.

**Measure, for each:**

    per-call G time at the same bias        (against 17304 us)
    compile/JIT time, once per class        (against 66.4 s)
    agreement with the numpy path           bitwise, or to 1 ulp,
                                            over a bias sweep INCLUDING
                                            the extremes the safety
                                            primitives exist for
    behaviour at +-1e30 and beyond          both arms still finite?

**Acceptance.** A candidate is worth building out if it reaches **at
least 50x** on `G` while agreeing with the numpy path over the full
sweep. Below 10x it is not worth the surface area. Between the two,
report and decide.

### 14.6 What it must not break, and this is the whole risk

The DSL's value is exact symbolic derivatives with **provable finiteness
under any Newton iterate**. A backend that is fast and loses either is
worth nothing here.

- **both arms of every conditional are still evaluated.** A backend
  that "optimises" a branch away changes the safety semantics the
  kernel is built on; `differentiable-numerics` records what that
  costs.
- **the range contracts still hold.** `safe_*` primitives spend
  exponent headroom deliberately; a different float path can move where
  they overflow.
- **the extremes are the test, not the typical bias.** Agreement at
  1.2 V proves nothing. The sweep must include the biases that motivated
  `expl` and `hypsmooth` in the first place.

### 14.7 What would make the answer "don't"

Recorded so that "no" stays a real option rather than a failure:

- the spike shows under 10x — the interpreter is not the bottleneck it
  looks like;
- JIT or compile cost per class exceeds what a simulation amortises;
- agreement at the extremes cannot be had without giving up the
  finiteness properties;
- **or nobody is actually blocked.** If the circuits being run are ten
  devices, 17 ms is a non-problem and the eight remaining library
  models are worth more. That is a question about use, not about the
  code, and it should be asked first.


---

## 12.7 Limiter parameters at the last accepted iterate (2026-08-25)

Closes §12.6(c) and the `vto` item of §12.4, and answers the open
question of whether the BJT should adopt `limit_together`.

### The rule that was right about the order and wrong about the conclusion

`_limit_par_fn` refused any limiter parameter that reached the solution:
*"evaluated from the parameters alone, and BEFORE the device is, so
they cannot depend on the solution."* The order is true. The conclusion
does not follow: a limiter is handed `x0`, the last accepted iterate,
precisely so it can measure the proposed step against it, and a
parameter evaluated **there** is as well-defined as `vold` is. It is
also exactly what SPICE does -- `mos1load.c` passes `fetlim` a `von`
recomputed from the *previous* iterate's bulk bias, not a card
constant.

So a parameter that reads the solution is now compiled with `x` as a
leading argument (through `_chain_compile`, chain read not inlined; an
`unpack` for the flat path) and the generated `limit()` calls it with
`x0`. The callable carries `_wants_x`; `explain()` prints
`[params at last iterate]`. Time is still refused.

Verified discriminatingly on both code generators: an `IS` that falls
thirteen decades per volt gives `pnjlim` a critical voltage of ~0.73 V
at `x0 = 0` and ~1.4 V at the proposed 0.9 V, so whether the clamp
fires says which point it read. It fires.

### What it retired

- **The thermal BJT's second saturation current.** `_gp_core` took a
  `tlim` argument so the self-heating class could hand the limiter an
  ambient-temperature copy of `isT` -- five duplicated lines and a
  critical voltage placed against a junction up to 100 K hotter than
  it thought. Gone; the limiter reads the heated `isT`/`vtT` directly,
  and the test that asserted the limiter was *independent* of the
  thermal node now asserts that it moves with it.
- **EKV's parameter-only `vto`.** Now `von = vtoT + gamma*(sqrt(phiT +
  vsb) - sqrt(phiT))`, evaluated at `x0`. The 100 V jump from an off
  device at 2 V of body bias lands at `von + 0.5 = 1.56 V` where it
  used to land at `vto + 0.5 = 1.0 V`; at zero body bias nothing
  changed. The test that pinned the 565 mV looseness now pins its
  absence.
- **Level 1's `vtoT`.** Same, with SPICE's forward-bias continuation
  built on the raw `V(b,s)`.

### Should the BJT adopt `limit_together`? Measured: no.

Batch 4 left this as "worth adopting, but measure first". Measured on
the shipped card, `gmin = 0`, plain `DC`, Jacobian evaluations:

    circuit          start    per-probe   grouped   same answer
    common-emitter    0 V        11         11        yes
    common-emitter   -5 V        11         11        yes
    common-emitter   20 V         7         92        yes
    current mirror    0 V         9          9        yes
    current mirror   -5 V         9          9        yes
    current mirror   20 V        29         29        yes
    (both, uniform 5 V start)  FAIL       FAIL     -- the ladder itself
                                                     gives up; not a
                                                     limiter question

Grouped is never better and is once **13x worse**. The per-probe form
stays. This is the third measurement pointing the same way (§10.3(b):
ties or +1 on the 48-point grid; batch 4: 552 vs 554 on Level 1):
**`limit_together` is a capability -- it makes four probes over four
terminals expressible at all -- and not a speed-up**, and a two-probe
device should not declare it.

### A crash seen once and not isolated, recorded so it is not lost

While writing the discriminating test: `IS * expl(-30.0 * b.V)` inside
`analog()` crashed in **sympy's `Piecewise` constructor** --
`KeyError: 'pop from an empty set'` from `_eval_rewrite_as_ITE` --
before the compiler was reached. `b.V.free_symbols` is `set()` (a
`Quantity` is atomic and reports none), and routing the same expression
through `vt()` (which carries the `TEMP` symbol) works. But the minimal
form `Piecewise((1.0, b.V < 1.0), (0.0, True))` does *not* reproduce
it, so the trigger is narrower than "a condition on a bare Quantity".
Not chased; the tests use `b.V / vt()` and this note is the only record.


---

## 15. "Is PCNR possible with all models?" -- measured 2026-08-25

**No: 2 of 26 shipped models qualify, and for transistors the obstacle
is structural.** PCNR has one scalar limited unknown per two-terminal
branch; a MOSFET's current is `f(vgs, vds, vbs)` and a BJT's transport
current `f(vbe, vbc)`. Every model, default card, first refusal:

    VBRANCH  (a branch-current unknown)   9   L, VSin, VCVS, RThermal,
                                              DiodeSpiceThermal,
                                              GummelPoonNpnThermal, Xtal,
                                              Comparator, OpAmp
    STATES   (idt / laplace_*)            5   Idt, Idtmod, FerriteBead,
                                              RSkin, Memristor
    MULTI-V  (reads other node voltages)  6   GummelPoon x2, Ekv x2,
                                              MosLevel1 x2
    linear   (no exponential scale)       3   R, G, and PSP at its FIRST
                                              contribution I(g,gi)/rg
    none     (no resistive current)       1   C
    QUALIFY                               2   DiodeHdl; DiodeSpiceHdl at
                                              rs = 0 only

"Partial participation" (qualifying junctions join, the rest stays in
MNA) was checked as a cheap shortcut and **does not** rescue the thermal
models: at `rth > 0` their junction current reads the thermal node
through `isT`, so it is MULTI-V in its own right. It would admit the
`rs > 0` diode and little else.

`explain()` now names the rule that refused (it used to say only "does
not qualify"; PSP's gate-resistor refusal took a re-measurement to
learn). `hdl.md` §3.2a's qualification bullets, which said "exactly one
branch" and "no charge" for months, are corrected.

### Vector PCNR: three premises corrected before anything was built

1. **The O(k) Schur trick survives vector unknowns.** `J_ml J_lm` is a
   sum over *unknowns*, not devices; a device with m probes over t
   terminals is m rank-one updates touching 2t entries.
   `doc/pcnr_native_design.md` §6(a) already said so.
2. **`refine` does not need sequential reading.** That was forced on the
   ordinary path by node write-back; PCNR has none. `limit_together` is
   not a consumer of vector PCNR -- under PCNR a four-probe MOSFET is
   expressible per-probe trivially.
3. **Chained models do have a traced PCNR twin** (`_pcnr_compiled`).
   The traced path's obstacle is that `jaxtransient`'s junction
   machinery is hard-wired to the two-terminal `[i, -i]` pnjlim shape.

And one trap for whoever builds it: `pcnr_junctions` pairs are consumed
as **gmin targets** by `dcanalysis._jrows` and
`jaxtransient._gmin_junction_rows`, on the ordinary path too. A `fet`
or `vds` probe must never reach them -- a gmin across `vgs` is a gate
leak.

### The gate, and where it stands

PCNR removes exactly one thing: two devices limiting the SAME branch
voltage fighting over a shared node. On transistors this was looked for
once before -- `doc/transient_work_plan.md` ("A hypothesis I tested and
had to withdraw"), `doc/pcnr_native_design.md` §8 -- and **did not
reproduce**: the second limiter reads the already-limited value and
composes with it. Every PCNR benchmark and test in the tree is diodes.

The full design, staged with the measurement first and an exit
condition written in advance, is in the plan file this section was
written from; the measurement's result follows.

### Stage 0 result: the clash REPRODUCES on FETs (2026-08-25)

Two independent measurements, both plain Newton (`StandardNewton`, no
ladder, no gmin anchor -- a `DC()` solve rescues every case below and
would have hidden all of it), Jacobian counts, each circuit built in
BOTH instance orders. Same equations, same x layout; only
`SubCircuit.limit`'s visit order differs.

**Measurement 1** -- 28 cases (mirror with a body-bias mismatch, diff
pair with a resistive tail, cascode, BJT mirror; EKV, Level 1, Gummel-
Poon; starts 0 / -5 / +20 V / one wild node): **25 order-independent, 3
signatures**:

    diff pair   EKV       start 0      5 vs 7        order-dependent count
    cascode     Level 1   start +20    FAIL vs 12    order-dependent failure
    BJT mirror  GP        start +20    FAIL vs 7     order-dependent failure

Traced, iteration by iteration (the 1e45/1e109 figures that appear in
such traces are the size of the Newton PROPOSALS being corrected, not
limiter outputs -- checked before reading anything into them):

- *cascode, Level 1, +20 V, M2 first*: M2's bulk-junction probes push
  the shared `bulk` node up by +9.5, +19.5, +49.5, +140 V on successive
  iterations and M1's push it back down by almost the same each time.
  A geometric tug-of-war over one node; diverges. M1 first: settles.
- *BJT mirror, +20 V, q1 first*: q1 limits the shared base; q2's
  base-collector probe is then anchored off the ALREADY-LIMITED base
  and lands the collector 2.8 V forward of it -- deep saturation,
  overflow. q2 first: the same probe anchors off the unlimited base,
  lands sane, converges in 7. Which point a derived node is anchored
  from depends on who wrote the shared node first.

**Measurement 2** (the design agent's, ideal 200 uA tail source,
5 kOhm loads, inputs 2.5 +- vin, `[M1 first, M2 first]`):

    vin          -1.0          -0.3        0        0.3         1.0
    both       [45, 14]      [26, 14]  [15, 15]  [14, 26]  [14, NoConv]
    group      [45, 14]      [26, 14]  [15, 15]  [14, 26]  [14, NoConv]
    Level 1    6 / 6 everywhere (polynomial channel, nothing to fall off)

Antisymmetric in `vin`, as a clash must be. Mechanism at iteration 0,
`vin = 1.0`: Newton proposes tail ~ -1.45e8 V; **M1's `fetlim` alone
wants tail = 2.5 V, M2's wants 0.5 V**, and dictionary order decides
which the tail gets. *(Corrected the same night by a trace: it is NOT
"last writer wins". M2, visited first, writes tail = 0.5; M1 then sees
the tail drifted only 0.5 V against its gate's 3.5 V from the zero
start, and its own move-the-wilder-node rule moves the SOURCE-PINNED
GATE `inp` to 1.5 V. Newton restores the gate and M1 is left at
vgs = 3 V from there on. Visited first, M1 writes 2.5 and M2 simply
does not bite -- composition, as the earlier BJT search found. So the
write-back rule of §12.1 is what turns a conflict into a wrong move;
the conflict itself is the paper's.)* In the losing order M1 is left at `vgs = 3 V` on a
subthreshold exponential and the tail creeps up by `N*VT` per
iteration until Newton falls over. `limit_together` cannot help -- it
is per DEVICE and this is BETWEEN devices. `DC()`'s ladder rescues
both orders to the same 2.818678 V.

**This is the paper's Figure 1 with FETs in place of diodes, and the
§10.4 gate is passed.** The one prior search for it (`transient_work_
plan.md`, two parallel BJTs) looked at a shared-GATE case, which
composes; the clash lives where two devices' limiters want DIFFERENT
values on one shared node -- a tail, a bulk, a derived collector.

**Before Stage 1: price the cheaper competitor.** `device_writeback`
already resolves competing probes WITHIN a device by a spanning
forest. Lifting it to `SubCircuit.limit` -- one forest over every
element's targets -- is classical limiting with no new unknowns. If it
makes the diff pair order-independent at ~14, vector PCNR's case is
capability again. If "two constraints on one node with different
values" forces the forest to drop one (it is a conflict, not a cycle),
that is the number that commissions Stage 1. Spike in progress; result
recorded below when it lands. Vector PCNR itself is NOT started, per
the approved scope.

### The cheaper competitor, priced: a circuit-level forest is order-independent and WORSE

`_limiting.CIRCUIT_LEVEL` (default off) + `element_targets()` +
`SubCircuit.limit` collecting every element's probe targets into ONE
`device_writeback` over the whole circuit. Plain Newton, gmin = 0:

    diff pair `both`   -1.0      -0.3       0       0.3      1.0
    baseline         [45,14]   [26,14]  [15,15]  [14,26]  [14,F]
    forest           [F, F]    [43,43]  [15,15]  [43,43]  [F, F]
    forest+pinned    [F, F]    [32,32]  [15,15]  [32,32]  [F, F]

    cascode L1 +20 V   (20,2,.8)  (40,2,.8)  (5,2,.8)  (20,2,.8,vb=-2)
    baseline           [17,10]    [F,10]     [F,12]    [F,21]
    forest             [13,13]    [12,12]    [19,19]   [21,21]

    48-point grid      both     group    mos4-group
    baseline           909/0    927/0    896/0
    forest            1878/8   1842/8   1866/9
    forest+pinned     1052/0   1051/0   1175/2

**Order-independent: yes. Better: no -- it fails in BOTH orders where
the baseline failed in one.** Two devices asking one node for 2.5 V and
0.5 V is a conflict, not a cycle: `(inp,tail)` and `(inn,tail)` are two
edges of a star, so the forest keeps both, anchors on the least-drifted
node and derives `tail = 0.5` then `inp = 1.5` -- a source-pinned gate
moved 2 V, the exact §12.1 failure, now reached identically in both
orders. Pinning VS nodes turns the conflict into a dropped edge, but
WHICH edge is decided by anchor choice, and that hands the tail to M2,
the device that is OFF. **Nothing in `|x - x0|` or in the correction
sizes knows which device is on.** The forest can only pick a loser
without the information to pick the right one.

Not adopted; left in the tree behind the switch, because the test that
records this (`test_a_circuit_level_forest_is_order_independent_and_
worse`) measures the intervention against its absence, per
`validation-design`.

**This is the number that commissions Stage 1.** The clash is real, the
classical fix cannot resolve it in principle, and PCNR removes it by
construction -- each device limits its OWN `vgs` unknown and the tail
is solved from both currents at their own limited voltages. Vector PCNR
is still NOT started: the approved scope was Stage 0 only, and this
result goes back to the user.

### Stage 1 delivered (2026-08-25): the BJT participates, and three things the design got wrong

`pcnr.py` now has per-device records (`PcnrDevice`: rows, probes,
kinds, offset), `augmented_system`/`schur_reduce`/`predict`/`refine`
over device blocks, and an **adapter** that presents every existing
scalar participant as m = 1 -- measured bit-identical on `g_lim`,
`J_ml`, `J_lm`, `f_eff`/`J_eff` and `refine` against the old formulas,
same iteration counts on the paper's Fig. 1 and `TwoJunction`. The
declared-probe route in `generate_code` admits a device whose every
current reaches the solution only through its `$limit` probes (pnj-only
in this stage); `explain()` prints `vector route, 2 probes over 3 rows
(pnjlim on (b,e), pnjlim on (b,c))` for the NPN. The O(k) Schur trick
survives as predicted (verified equal to the dense product). Finite-
difference `J == dg/d[x;v]` on the full NPN card with charge, 2e-4. Same
answers as the ordinary path to 1e-9. `dcanalysis` untouched: the pair
view is filtered to `pnj` probes, so no gmin ever lands across a `vgs`.
+29 tests; suite 2335 / 7 / 3 / 0.

**1. "Assemble, then subtract the participant" is inexact, not merely
fragile.** hdl.md 3.2a had recorded `inf - inf = nan` and concluded
"use `limexp` under PCNR". With the participant kept finite, its
current at the unlimited node reached 1e72 and the subtraction left
~1e56 of ORDER-DEPENDENT noise: the same mirror took 11 iterations in
one instance order and 179 in the other. The participant is now
**excluded from assembly**, never evaluated at the node voltages. The
raw-`exp` diode converges under PCNR; the `limexp` test was inverted;
the charge diode's 14 DC iterations became 8, the 14 having been a
defect of the cancellation. Twelfth right-conclusion-wrong-reason: the
`nan` was real and the lesson drawn from it ("harden the participant")
was the wrong one.

**2. Order-independence is delivered; convergence from +20 V is not.**
The Stage 0 acceptance case, full NPN card mirror, `[q1 first, q2 first]`:

    start     plain Newton     PCNR
    0 V        [9, 9]          [9, 9]
    -5 V       [9, 9]          [9, 9]
    +5 V       [FAIL, FAIL]    [164, 164]     <- rescued where the ordinary
                                                 ladder itself gave up
    +10 V      [7, 7]          [6, 6]
    +20 V      [FAIL, 7]       [FAIL, FAIL]   <- order-independent, and fails

The structural claim holds everywhere -- PCNR never depends on instance
order. But from +20 V it fails in BOTH orders, and the trace says why:
SPICE `pnjlim`'s `vold <= 0` branch maps a 9e45 V proposal to
`vbc = 2.83 V`, 1e33 S, singular. Plain Newton's one-order "win" was a
stale-anchor accident (2.83 V written against the not-yet-limited base,
which then moved, leaving `vbc = -12.8 V`). Behind that sits a second
obstacle: PCNR carries `-9e45` in the accepted iterate and the
correction's ulp is 1e30. **So the +20 V case is a limiter-LAW question
(what `pnjlim` does with an astronomically wild proposal), not an
architecture one**, and it is the same question on both paths.

**3. JAX + `pcnr=True` + a vector device now raises `NotImplementedError`
at setup** where the BJT previously, silently, did not participate.
Stage 3.

Also fixed: a degenerate probe on one row (diode-connected `(nb,nb)`)
must accumulate into dense `J_lm`, not assign. `transient.py`'s
hand-written Schur loop and `fang_timestep`'s `dx_lim` copy are gone in
favour of the shared functions. Per-component convergence replaces the
`max|v|`-scaled one (the old criterion accepted `[2e-5, 0]` against
`[0.7, 40]`).

**Stage 2 (FETs) and Stage 3 (JAX) remain the user's call.** Stage 2 is
where the diff-pair clash -- the case that commissioned all this -- gets
its number.

### Stage 2 delivered (2026-08-26): the diff-pair clash is gone, and the old solver had been declaring fake convergences

`fetlim`/`limvds` probes are PCNR unknowns. A cycle in the declaration
(`(b,s)`, `(b,d)`, `(d,s)`) is handled by taking a spanning TREE of
probes as unknowns and applying the redundant probe's law over the tree
(Kruskal on |correction|, `pcnr.limit_block`); `g_lim` consistency and
the implied `vbd` are asserted to 1e-12 at convergence. Measured: with
the redundant law made a no-op, Level 1 fails **48 of 48** grid points
from a 20 V start -- `vbd = vbs - vds` is left unbounded forward and
the bulk-drain diode overflows -- so the law over the tree is
load-bearing, and `limit_together`'s coupling is NOT "trivially
per-probe under PCNR" when the declaration has a cycle (a premise of
the design, corrected). EKV is refused as predicted (`var(vsb) reads
b, g, which no $limit probe limits`); the minimal admission is an
identity probe kind (`limit_identity(b.V)`: an unknown with no law),
proposed and not built. +15 tests; suite 2351 / 7 / 3 / 0.

**The number this was commissioned for** -- `_fet('both')` diff pair,
ideal tail, `[M1 first, M2 first]`:

    vin        -1.0       -0.3        0        0.3       1.0
    plain    [45,14]    [26,14]  [15,15]   [14,26]  [14,FAIL]
    PCNR     [22,22]    [23,23]  [15,15]   [23,23]  [22,22]

Order-independent at every point AND converging everywhere, including
`vin = +1.0` where plain Newton failed in one order; tails equal to
`DC()`'s in both orders. Price: 22 against the lucky order's 14.

    48-point cascode grid    plain      PCNR (budget 200)
    _fet('both')             909 / 0    1281 / 0
    _mos4('group')           896 / 0    1430 / 0

~40% dearer, every answer equal on the nodes; the slowest point (152)
is a `limvds` limit cycle riding Newton's `N*VT`-per-step creep down a
1e25 A exponential -- slow, not divergent.

    Level-1 cascode, 20 V start     plain      PCNR
    (20, 2, .8)                     [13, 8]    [12, 12]
    (40, 2, .8)                     [17, 8]    [12, 12]
    (5, 2, .8)                      [31, 9]    [19, 19]
    (20, 2, .8, vb = -2)            [13, 8]    [FAIL, FAIL]

Order-independent at all four, converging at three; the fourth is
Stage 1's obstacle again -- from an all-off start `mid` is held by
nothing, PCNR limits no NODE, and the iterate reaches -3.4e117 --
architecture, the same on FETs as on BJTs. (The plain column no longer
matches Stage 0's `[17,10] [F,10] [F,12] [F,21]`: `von` now follows the
bulk bias, §12.6(c), which Stage 0 predates.)

**Two defects in the EXISTING solver, found on the way:**

1. **`solve_dc` was declaring fake convergences.** Its criterion was
   `max|dx| < reltol * max|x_new|` over the whole vector, source branch
   currents included. On the `fet` cascode at (40, 1, 1.2) the VDD
   current hit 7.5e27 A, the node tolerance became 7e23 V, and it
   "converged" in 5 iterations at `mid = -2.12 V` with a KCL residual
   of 7.5e27 A -- **18 of 48 grid points, up to 5e9 V off**, all
   reported as converged. Diodes never exposed it. Replaced by
   `StandardNewton`'s per-component step test plus a residual test.
   Two existing counts moved because of it (charge diode 8 -> 9, BJT
   +10 V 6 -> 7): each had stopped with a residual 4-5x over tolerance
   that the old criterion never examined.
2. **`DC(pcnr=True)` gated PARTICIPATION on the pnj-only pair view**,
   which exists for the gmin ladders and is EMPTY for a circuit of pure
   `fetlim`/`limvds` devices -- so a MOSFET differential pair asked for
   PCNR silently got the ordinary solver. Fixed in `dcanalysis` and
   `transient` (gate on `pcnr_devices()`), with a test that counts the
   `solve_dc` calls. Every Stage 2 table was taken through
   `pcnr.solve_dc` directly, which is how the fall-through was noticed.

**Bottom line for the question that started §15.** Vector PCNR does
what it was commissioned for: the inter-device clash is gone in both
orders, converging where plain Newton failed, at a ~40-60% iteration
premium. It does not remove the unlimited-node failure (a node no probe
touches, driven to 1e117), which is architectural and shared by every
path. **Stage 3 (JAX) is the user's call**; a vector device under
`pcnr=True` on the traced path currently raises `NotImplementedError`
rather than silently not participating.

### The "unlimited-node failure", re-measured (2026-08-26): not what was recorded, and the real gap fixed

Listed after Stage 2 as "the one failure vector PCNR could not remove,
architectural on every path". Re-measured at every ENTRY POINT rather
than under the plain-Newton harness the Stage tables used:

    Level-1 cascode (20,2,.8, vb=-2)   uniform 20 V start   zero start
      DC() default (ladder + anchor)        52 (0.7796 V)    136 (0.7797 V)
      DC(pcnr=True)                          24 (0.7796 V)     17 (0.7796 V)
      pcnr.solve_dc direct                   30 (0.7797 V)     22 (0.7797 V)

**It converges at every entry point, from both starts.** Stage 2's
`[FAIL, FAIL]` at this point does not reproduce with this harness; the
row is left in the Stage 2 table with this note beside it rather than
deleted, because the two harnesses have not been reconciled.

    BJT mirror, uniform 20 V start
      DC() default                          146 (3.996 V)
      DC(pcnr=True)                         FAIL LinAlgError
      pcnr.solve_dc direct                  FAIL LinAlgError

**This one is real, and it is an integration gap, not an architectural
limit**: `DC(pcnr=True)` called `solve_dc` and returned. A PCNR failure
was a hard failure, on a circuit the default path solves, because the
PCNR path had NO rescue chain while the ordinary path has four ladders
plus the anchor.

Three native rescues were built as adaptive ladders around `solve_dc`
(the driver `_adaptive_conductance_ladder` is generic and
`augmented_system` already exposes the hooks) and **all three fail**
on this circuit, in both instance orders:

    gmin shunt to ground (e0 = -3)     FAIL after 673 / 727 evaluations
    source stepping (lambda from 1e-3) FAIL at the FIRST rung
    Jacobian damping (e0 = -3, -5)     FAIL after 673 / 727, 300 / 300

The failure is PCNR's UNDAMPED first `predict` step from a wild start
-- 9e45 V proposed, which `pnjlim`'s `vold <= 0` branch maps to a
2.83 V forward junction, 1e33 S, singular -- and it is the same on
every rung, because every rung takes that step. (A first damping
attempt was inconclusive rather than negative: `solve_dc` raises
`RuntimeError`, which the ladder driver does not catch, and 1 S of
damping against 1 kOhm loads contracts by 0.999 per step. Wired
properly, it is negative.)

**Fix shipped:** `DC(pcnr=True)` and the transient's PCNR step fall
back to the ordinary chain when PCNR raises, log it, and set
`DC.pcnr_fell_back`. Order-independence is lost for THAT solve; the
answer never is. Pinned by two tests -- the fallback fires on the
mirror and matches `DC()`'s answer, and it does NOT fire on a solve
PCNR handles (a fallback that always fires would pass the first test
alone). The experiment hooks were removed from `solve_dc` after the
measurement; they have no user.

What remains genuinely open is narrower than the heading said: **a
native rescue for PCNR from a wild start**, which would need a damped
`predict` (not a damped rung), and has no measurement demanding it yet.


---

## 16. The compile cache (2026-08-26): profile first, then cache -- and a harness that had drifted

### Profile, before any change (`elements_hdl` cold import, 26 classes, 6.0 s)

    _chain_compile (chained path)        4.15 s   69%   3.58 s of it Jacobian
                                                        compiles; inside, sympy
                                                        diff ~ printing ~ half each
    resolve()  (Quantity -> x subs)      ~1.0 s
    _run_analog                          0.36 s
    sympy.lambdify (eager path)          0.18 s
    _pcnr_declared_route                 0.16 s
    Matrix.jacobian (eager path)         0.11 s
    _chain_prune / forward_d / simplify  0.03 / 0.00 / 0.00 s

The cost is diffuse sympy work (assumptions engine, `Mul.flatten`, the
printer), not one hot spot -- so a cache IS the right fix. The brief's
suspicion of the PCNR detector's `simplify` was wrong: two calls, nil.

**One waste fixed first, bit-identical:** `i_dc`/`G_dc`/`u_dc` were
compiled separately although identical to `i`/`G`/`u` for every class
without a DC-pinned state -- all 26 -- and `G_dc` on the chained path is
a full forward-mode Jacobian compile. Shared when `dc_pins` is empty:
**6.0 -> 5.1 s (-15%)**. `BehaviouralMeta` also recomputed
`Matrix.jacobian(ivec)` to decide `linear`; it now reads the Jacobian
it already has.

### The cache (`_hdl_cache.py`)

Key: format number; source hashes of `hdl.py` and `_hdl_cache.py`;
sympy/numpy/Python versions; module + qualname; `inspect.getsource
(analog)`; every instparam declaration; the collapse mask; class-body
`terminals`; the class's module source and that of every pycircuit
module a named global lives in; **and the body's bindings** -- closure
cells, defaults, referenced globals (sympy by `srepr`, arrays by bytes,
functions by source + their bindings). That last item was not in the
brief and was REQUIRED: the suite builds classes with identical analog
text whose behaviour is decided by an enclosing variable, and with
source-only keys 28 tests failed on hits from the wrong twin.

Refused (cold compile, status says why): unrecoverable source; a body
that writes into a non-local; a global with no stable repr; a rebuilt
namespace that would bind a loaded global to a different object.
Cached: the whole `_hdl_info`, functions as generated source
(`_src`) re-executed on load; `sym_spec`/`pure_spec` as pickled sympy.
Not cached: the JAX twins, built lazily at first traced use.
`$PYCIRCUIT_HDL_CACHE=0` disables; `_DEBUG=1` explains a miss;
`cls._hdl_cache_status` always says which.

Bit-identity: all 26 elements + 5 test models, `i,G,q,C,CY` bytes at
four points, `explain()` text, `limit()` output, collapse variants --
identical cold vs hit, in-process and across two interpreters.
Invalidation: one character of `analog`, an instparam default, the
format, the compiler hash, a helper module, a closure value -- all
miss. Damage: garbage / truncated / foreign / unrebuildable entries
miss and repair; two concurrent processes on one key leave a valid
entry.

### Numbers (fresh interpreter per column)

                                       cache off   cold miss   warm hit
    import elements_hdl (26 classes)     5.11 s      5.17 s     0.14 s
    GummelPoonNpnThermalHdl              1.08        1.07       0.018
    GummelPoonNpnHdl                     0.80        0.78       0.009
    PspMosLongChannel (compact import)     --        41.4       0.59

§13's "PSP compiles in ~66 s" included instantiation and parameter
resolution; the compile alone is 41 s. Either way it is now 0.6 s
warm. `benchmarks/hdl_model_cost.py --cache` reproduces the table.

Suite: 2390 passed with the cache ENABLED and with it DISABLED (36 new
tests; 2354 vs the 2353 baseline is one environment-dependent skip
that passes here).

### A harness that had drifted, found by this work

The chunked suite recipe -- `sed -n '1,33p' / '34,66p' / '67,99p'` --
was written when there were 97 test files. There are 101. Since
`f029d6c` (vector PCNR Stage 1) `test_vss_gear2.py` has been outside
every chunk, and every "verified 2335 / 2351 / 2353" in §15 was taken
without it -- including across four commits that touched
`transient.py`, which it tests. Re-run: it passes, so nothing was
broken, only unverified. The lesson is in `validation-design`: a
chunked run is a coverage claim, and it must be reconciled against
`--collect-only` EVERY time, not once. The recipe now derives its
ranges from the file count.


---

## 17. The parameter namespace (2026-08-26): S5's second half

`params_as = 'p'` on a `Behavioural` class makes `analog()`'s FIRST
argument a `ParamNamespace`: `p.bf`, `p['lambda']`, `p.given('rb')`,
`p.names`. Its attributes ARE the parameter symbols, so generated code,
`_hdl_info` and `explain()` are byte-identical to the bare-name
spelling -- pinned for an eager and a chained model, with a mutation
control. Bare names stay the default and every other model is
untouched; the two styles coexist across classes, not within one, and
each style's mistakes get a message naming the class and the declared
names.

**Adopted where it removed a workaround:** `_gp_core(p, ...)` and
`_spice_diode(p, a, c, T)` -- the latter down from 19 arguments -- with
`_with_params` deleted. Bit-identity of the five affected classes
(three BJTs, two SPICE diodes) is pinned as SHA-256 digests of
`i/G/q/C` at three random points on a full card plus `explain()`,
taken from the pre-change source and re-taken after: 20 of 20 equal.
`MosLevel1Hdl` now declares SPICE's own `lambda` and `as` through
`Behavioural.aliasparams` (which existed; `Parameter` itself has no
aliases), canonical names unchanged.

**Three things the brief got wrong**, all found by building it:

1. *"The first argument named `p`"* as the opt-in would have broken
   shipped tests -- `analog(p, n)` is a PIN spelling in ten existing
   models. Hence the explicit class attribute.
2. A keyword-named parameter needs SYMBOL MANGLING: the let-chain
   printer emits `def _f(x, lambda, ...)` verbatim -- a SyntaxError --
   and four JAX sites rebuild symbols from names. `_param_symbol()`
   maps unsafe names to `_hdl_kw_<name>` at all five sites.
3. A bare class and a `params_as` class with the SAME analog text would
   have shared a cache key (different pin count). `params_as` is in the
   key now, with a test.

The claim that attribute access is visible to linters was not machine-
checked: no linter is installed in the venv.


---

## 18. Library batch 5 and `limit_identity` (2026-08-26): the list is built, and a DC-pin bug on the chained path

Eleven new classes (26 -> 37): `VcoHdl` (phase noise injected into the
FREQUENCY through an `idtmod` phase, so jitter random-walks -- output
PSD `(2 pi va)^2 (sf + kff/f) / w^2` to 1e-9 over three decades),
`DividerHdl` / `MixerHdl` / `ChargePumpHdl`, `PhotodiodeHdl` / `LedHdl`
(single-diode closed forms, implicit I-V with `rs`+`rsh` to 1e-8, LED
linear above threshold, optocoupler CTR), `MesfetCurticeHdl` /
`MesfetStatzHdl` / `HemtAngelovHdl` (Curtice 1980, `mes1.va`, Angelov
1992), `MosLevel3Hdl` / `Pmos` (a numpy transcription of `mos3load.c`
with every knob on, 14 biases, 1e-9; = Level 1 at gamma = 0). All
written with `params_as = 'p'` -- "the first thing reached for, zero
workarounds"; `lambda` and `as` declared under SPICE's own names.
+64 tests, 20 mutations (18 caught, one test added for the survivor,
one mutation inexpressible). Suite 2482 / 7 / 3 / 0.

**`limit_identity(probe)`** -- kind `'id'`, no law, `apply_limit`
returns `vnew` itself, never bites on the ordinary path (verified
bitwise: 100 wild steps write nothing), refused by `limit_together`
(a grouped probe that did not bite is held as a constraint, and an
identity has nothing to hold). It declares a PCNR unknown, and EKV now
qualifies: `PCNR: vector route, 3 probes over 4 rows (identity (no
law) on (s,b), fetlim on (g,s), limvds on (d,s))`. Measured on the
diff pair, EKV `w = 20 um`, `[M1 first, M2 first]`:

    zero start   vin   -1.0    -0.3      0      0.3    1.0
    plain            [7,5]   [7,5]   [7,7]   [5,7]  [5,7]    Stage 0's 5 vs 7
    PCNR             [7,7]   [7,7]   [7,7]   [7,7]  [7,7]
    +20 V start
    plain            [6,6]   [6,6]   [6,6]   [6,6]  [6,6]
    PCNR             [F,F]   [F,F]   [F,F]   [F,F]  [F,F]    at the FIRST predict

Order-independence delivered; convergence equal from zero and WORSE
from a uniform wild start: both channels cut off, the tail's whole
conductance is 2e-19 S, plain Newton's LU pivots on it and the limiter
repairs the step, PCNR's Schur complement rounds it to a zero pivot.
That is the unlimited-node case of §15 again, on a FET -- and
`DC(pcnr=True)` falls back and gets the answer, which is what the
fallback of §15 was for.

**A DSL bug, fixed outside the batch's stated scope:** on the chained
(`var()`) path `i()`/`G()` returned the TRANSIENT stamps
unconditionally, so an `idt`/`idtmod` with an `ic` was never DC-pinned
and the operating point was singular. `IdtmodHdl` is flat and never
showed it; the VCO did. Eight lines plus
`test_a_chained_idtmod_is_pinned_at_dc`; reverting breaks six tests.
Sixth instance this week of machinery written against the eager path
and silently wrong for the chained one.

**`Cross`, second production user, exonerated.** Batch 2's ten-decade
accuracy scatter and 3x step cost were properties of its LATCHING
comparator, not of `Cross`: on the divider it is uniformly first-order
(worst edge 1.7e-5 -> 8.2e-6 s) at 1.03x the steps.

Friction that remains: `idtmod(expr)` refuses noise in `expr` and
V-contributions refuse noise, so frequency noise needs an internal node
and a 1 S contribution; ONE LAW PER BRANCH -- SPICE's MESFET runs
`pnjlim` and `fetlim` on the same `vgs`, which is not expressible;
Level 3's source/drain exchange has no antisymmetric form, so it is a
`Piecewise` with two channel chains (1.2-1.4 s per class, cold); a
helper's internal node is unreachable from its caller (`_spice_diode`
grew a `junction=` hook, bit-identical when unused); Level 1's `vtoT`
omits `mos1temp.c`'s `(Eg(tnom) - Eg(T))/2` term (pinned; Level 3 has
it). And the limiter measurement on a bounded tanh channel: the
`fetlim`/`limvds` group COSTS iterations (Curtice `[5,5] [4,3] [6,3]`
limited vs not); what earns its keep on a MESFET is `pnjlim` on the
Schottky gates.

Cache: cold import of 37 classes 9.3 s (Level 3 x2 is 2.2 s of it),
warm 0.91 s, 37/37 hit, none refuse.


---

## 19. The backend spike, measured (2026-08-26): C at -O2 is 219-229x and bitwise exact

§14.5's spike, run as specified (harness `benchmarks/backend_spike.py`,
full report `doc/backend_spike_260826.md`). Sweep = 4096 extreme-grid
combinations of ±1e30/±100/±1/0/0.7/1.2 on all four terminals plus 889
reference-sweep biases = 4985 points; the orchestrator's independent
re-run used 1189.

    candidate                     G us     speedup   compile   agreement   finiteness   verdict
    numpy (today)               17 356        1x         --
    numba, default              JIT never finished (> 600 s on G)                       don't (14.7)
    numba, NUMBA_OPT=0             622       28x       68 s     2.4e-10      kept        between; don't
    C gcc -O2, bitwise-faithful   75-79    219-229x  22-32 s   0 ulp, all   kept        WORTH BUILDING OUT
    C gcc -O2, x*x for **2          14     1235x     19 s      2.4e-10      kept        if 2.4e-10 accepted
    llvmlite                      skipped: numba shows LLVM's optimiser is what does not finish

`i`: 1.27 ms -> 8 us (159x), same fidelity.

**Both-arms semantics survive** because they are a property of the
CALL: arms are computed as arguments before `where` picks, and both
candidates keep them as arguments (`_sel(c, a, b)` in C, an `@njit
_where` in numba); neither uses inline `if/else` or `a*c + b*(1-c)`.
The one ulp is instructive: every one of ~740 intermediates matches
bitwise except `**2` -- numpy's scalar `**2` is glibc `pow`, not
correctly rounded, and gcc folds `pow(x, 2)` to `x*x`, so the compiled
path was MORE accurate and a cancellation in the derivative chain
amplified 1 ulp to 2.4e-10. Bitwise fidelity costs 5.4x
(`-fno-builtin-pow`, 5618 real `pow` calls) and is still 229x.
§14.6's overflow-migration risk never occurred; the numpy path itself
is NaN at 96 extreme points (Vb and Vs/Vd both at ±1e30), which is
model-level and unchanged by any backend.

**§13's census was wrong in reason.** "3471 calls, 16% dispatch, 84%
interpreter" counted numpy calls; the real Python-level census of `G`
is ~21 400 calls (`minc` 3833, `maxc` 2911, `_step` 2046, `_rdiv`/
`_recip2` 2775, ~2800 comparisons ...) and the kernel helpers alone are
~41% of the time. §14.5's "numba needs no new code generator" was
wrong: three rewrites and it still does not JIT. §14.4's forbidden
"300-1000x" was, by accident, about right (1235x non-faithful, 229x
faithful). Fifteenth right-conclusion-wrong-reason.

**What building it takes (days, not weeks):** a `_CChainPrinter`
beside `_ChainPrinter`, printing the same sympy lists `_chain_compile`
already walks (NOT re-parsing `_src`); a `_KERNEL_C` prelude beside
`_KERNEL_NUMPY`; `backend=` on `_chain_compile`, keeping `_src` as the
numpy reference for `explain()`, the cache and the tests; the `.so`
keyed by C-source SHA + compiler in `_hdl_cache` with numpy fallback;
an `(x, p, out)` convention with parameters packed once per
`update_iparv`. cffi releases the GIL; ~0.6-1 MB per class. Tests: the
bitwise sweep, a both-arms overflow model, a `pow` sentinel, a cache
round-trip.

**Whether to build it is the user's decision**, per §14.7 -- and its
first question is unchanged: whether 17 ms per Jacobian blocks
anything actually run.

**A record incident from the same day, kept here because it is the
kind that leaves no trace:** the commit before this one deleted
§12-§15 of this file -- 1300 lines -- by patching a status-table row
with a regex that searched for "the next numbered row" and found one
in §15's tables. No test can see a document shrink; the spike agent
noticed because §14, the section it was working from, was gone. The
file was restored from the previous commit and the row re-patched with
the search scoped to the table.


---

## 20. The C backend, built (2026-08-26)

§19's spike said "worth building out"; this is it. Default is still
**numpy** -- `PYCIRCUIT_HDL_BACKEND=c`, `hdl.set_backend(which, cls)`
or a per-class `hdl_backend` turns it on, `explain()` prints the
backend, and `cls._hdl_backend_status` always says what actually
happened (`'c'`, or `'numpy (<why not>)'`). C requested and
unavailable warns and falls back; it never silently degrades.

`_CChainPrinter` derives from the SAME `NumPyPrinter` subclass as
`_ChainPrinter` (refactored to a classmethod), so Add/Mul structure,
parenthesisation, Piecewise nesting and Min/Max binarisation are
literally the numpy path's code -- only leaf syntax is overridden.
`_KERNEL_C` sits beside `_KERNEL_NUMPY` with NaN semantics written as
numpy defines them. `_hdl_cbackend.py` does compiler discovery, strict
flags (`-O2 -fno-fast-math -ffp-contract=off -fno-builtin-pow/exp/log`)
and the `.so` store; the cache carries the C text (format 2).

**Measured (independently re-run by the orchestrator):**

    class                numpy        C        speedup   bitwise
    MosLevel3Hdl        1861 us     6.9 us      270x      yes
    GummelPoonNpnHdl     340 us     7.6 us       45x      yes
    EkvNmosHdl            93 us     4.6 us       20x      yes
    PSP103 G           17 508 us    82.8 us      212x      yes
    PSP103 i            1 282 us    10.7 us      120x      yes
    PSP103 q / C                                 89x / 215x

Suite **2547 / 6 / 3 / 0 in BOTH backend states**, chunks derived from
N=105 and reconciled against `--collect-only` (2553).

**Bit-identity, with its exceptions named rather than hidden.** PSP's
`i`, `G`, `q` are byte-identical over the spike's full 4985-point
sweep; `C` is value-identical with 384 sign-of-exact-zero byte
differences. All 22 chained library classes x {i,G,q,C} over 56 points
including +-1e30 are value-identical, and byte-identical except:

1. **tanh** -- numpy ships its own, differing from libm by 1 ulp on
   ~30% of arguments. Every tanh class is bitwise-identical to its own
   source run with libm tanh, which converts the tolerance into a
   cause;
2. **the sign of exact zeros** -- numpy computes some subchains in
   int64 (`where(c,1,0)`), where `-1 * 0 = +0`; C is all doubles, where
   `-1.0 * 0.0 = -0.0`. Values compare equal.

**Deviations from §19, each forced by a measurement:** the C text is
emitted ALWAYS (+19% on a cold compile) so a warm cache switches
backend with zero sympy -- keying the cache by backend would have made
the switch a 41 s recompile; the compiler's identity lives in the `.so`
FILENAME, not the key, so a warm store serves a machine with no
compiler at all (tested); `q` and `C` are compiled too (89x, 215x);
`T` is a slot in the packed parameters; and `set_backend` must chase
collapse variants -- without that MosLevel3 ran at exactly 1x, caught
by the benchmark and not by any test.

**Nine brief claims wrong**, the sharpest being that `-fno-builtin-pow`
is necessary but NOT sufficient: numpy's `**2` on the 0-d array a
`where` returns is `square`, i.e. `x*x` -- the opposite correction --
so the printer emits `_sq` for Piecewise-based bases at exponent 2 and
real `pow` elsewhere. Sixteenth right-conclusion-wrong-reason.

**A claim from the backend report, checked and WITHDRAWN
(2026-08-26).** It said six chained classes "raise
`ZeroDivisionError` from `G` at DEFAULT parameters on the numpy path".
Measured across all 37 library classes at 28 bias points each --
zero, four methods, defaults: **no class raises anything.** What does
raise is the UNCOLLAPSED BASE function (`cls._hdl_info['funcs']['G']`
called directly), and no instance ever runs it: `Collapse` retargets
every instance to a variant with the dividing branch removed, which is
exactly what `rc = 0` is declared to do. The `KW` overrides in
`test_hdl_cbackend.py` are therefore a harness convenience -- they let
the sweep exercise the parasitic branches -- and not a defect list.

The behavioural difference the report noticed is real and narrower:
where the numpy path would raise, C returns `inf`. That matters only
for code that calls a base function directly, which is test-harness
code. Seventeenth right-conclusion-wrong-reason of this campaign: the
`ZeroDivisionError` was observed, and the conclusion drawn about who
suffers it was wrong.

A c-mode suite run leaves ~440 small `.so` files (~157 MB with the
pickles) in the user cache; `_hdl_cache.clear()` removes them.


---

## 21. S3 -- automatic both-arms domain clamping (2026-08-26)

`hdl.select((expr, cond), ..., (default, True))`: a `Piecewise` that
substitutes, inside each arm, a clamped copy of every symbol the arm's
own selection region bounds. The hand-written

    Piecewise((0.0, p <= 0.0), (1.0 / maxc(p, 1e-30), True))

becomes `select((0.0, p <= 0.0), (1.0 / p, True), margin=1e-30)`, and
the `maxc` is inserted by the compiler -- which is domain propagation
the compiler always had the information to do.

**Why it is exact.** `minc(v, c)` returns `v` bit for bit wherever
`v <= c`, and differentiates to `_step(c, v) = 1` there with *exactly
zero weight on `c`*. So a clamped arm agrees with a hand-written one to
the last bit of the residual AND of the Jacobian -- which is the half
that matters, since the chain rule multiplies a zero partial by the
discarded arm's derivative and `0 x NaN` is `NaN`.

**Supported, and nothing else:** `v < c` / `<=` / `>` / `>=` and their
mirrors; `u < v` between two atoms (both clamped simultaneously);
`|v| < c` and `v*v < c`; `|v| > c` (pushed out of the hole); `And` as
the intersection of boxes; `Or` in a *negated* position by De Morgan;
`True` as the complement of the earlier arms. `c` may mention `v`
itself (`vds < vdsat`) because the substitution is simultaneous, so the
bound is evaluated at the unclamped `v` and no cycle appears -- and it
costs nothing in finiteness, because the condition was evaluating that
bound anyway.

**Everything else is REFUSED at build time** with `SelectRefused`
naming the condition and the sympy class that made it unsupported --
an `Or` in a positive position (a union of boxes is not a box), an
equality, a relational whose sides are both compound. That refusal is
the design's load-bearing half: **a `select` that silently failed to
clamp would be worse than `Piecewise`**, because the author would trust
it.

**Adoption: 6 of the 20 surveyed sites**, all in `elements_hdl.py`
(the four `_recip` guards in the Gummel-Poon, and two more).
`psp_kernel.py` was NOT touched -- its four sites are inside a model
validated to 1e-6 against a vendor, and the conversion was not worth
the risk for this change.

**Verified independently by the orchestrator:** all **37** library
classes digested over 30 random biases x {i, G, q, C} with the
adoption in place and again with `elements_hdl.py` reverted --
**37 identical, 0 differ**. Suite **2622 / 6 / 3 / 0 in BOTH backend
states**; `sphinx -W` clean.

**What it does not absorb:** 14 of the 20 sites remain hand-clamped.
The mechanism covers their shapes, so the residue is adoption work,
not a limitation -- with `psp_kernel`'s four deliberately excluded
until someone wants to re-validate PSP against the vendor afterwards.


---

## 22. S1, measured before building it (2026-08-26): the two fast paths are DISJOINT

S1 -- a pure form for the chained path -- has one remaining
justification: `solve_batched`. Measured before committing days to it.

**The 22.5x reproduces.** `benchmarks/batched_sweep.py`, unchanged, on
this machine's GPU:

    lanes    cpu loop    jax warm   speedup
        8      819 ms     1726 ms      0.5x
       32     3377        1932         1.7x
      128    13451        2060         6.5x
      512    53462        2577        20.7x

**But the fast paths do not overlap at all.** Every `Behavioural` class
in the library, partitioned:

    chained  -> the C backend (20-270x)        22
    pure_spec -> solve_batched (up to 20.7x)   15
    BOTH                                        0
    NEITHER                                     0

A model is chained *because* it calls `var()`, which is what makes it
too big to differentiate as a tree -- and that is exactly the property
that denies it a `pure_spec`. So the 22 models that carry real physics
(PSP, the BJTs, EKV, both MOS levels, the MESFETs) get C and cannot
batch; the 15 that can batch (R, C, L, the ideal diode, the passives)
get no C at all. Verified directly: running the batched-sweep
rectifier with HDL elements and `set_backend('c')` changes nothing and
reports `numpy (eager path: the C backend serves let-chain models
only)` -- the backend refuses honestly rather than pretending.

### What this does to S1

S1's value is no longer "batching for the flagship model" in the
abstract; it is **whatever batching wins ON TOP OF a C backend that is
already 212x on PSP's Jacobian**. That is not measurable until one of
the two exists for the same model, and it is plausibly small: the C
path removes the per-evaluation cost that batching was amortising.

### And a gap the roadmap never named

The partition is symmetric, and only one half was ever written down.
**The eager path has no C backend** -- 15 shipped classes, including
every passive, evaluate through `lambdify` at numpy speed. Emitting C
for the eager path looks *cheaper* than S1, not dearer: there is no
let-chain to walk, the expressions are already flat, and
`_CChainPrinter` exists. Neither has been measured against the other.

**So S1 is NOT the next thing to build.** The honest ranking is:

1. measure what fraction of an eager transient is device evaluation --
   the rectifier CPU loop is 14.2 s for 128 lanes with HDL elements
   against 13.5 s with hand-written ones, which suggests the SOLVER
   dominates and that C for the eager path may buy little;
2. if it does dominate, neither C-for-eager nor S1 is worth building,
   and the batching win at 512 lanes is already the answer for the
   models that can use it;
3. if device evaluation is significant, C-for-eager is the cheaper of
   the two and should go first.

Recorded rather than acted on: this is the roadmap's own rule for the
backend work -- gate the expensive item on a measurement -- applied to
the item that was next in line.

## 23. §22's measurement, taken (2026-08-26): the device is not the cost, and the marshalling is

§22 ranked three items and put a measurement first: *what fraction of an
eager transient is device evaluation?* Taken. It answers both open
backend items, and the answer to each is **no**.

### The 5% inference in §22 was wrong

§22 reasoned from "14.2 s with HDL elements against 13.5 s with
hand-written ones" that device evaluation was ~5% and the solver
dominated. That comparison **confounds two things**: both sides pay the
same circuit-level assembly cost, so the difference measures only the
gap between two element implementations, not the size of the phase they
sit in. Profiled directly on the batched-sweep rectifier's single lane
(4 elements, `tend=2ms`, `timestep=1e-5`):

| phase | ms | share |
|---|---|---|
| circuit assembly (`cir.i/G/q/C/u`) | 77.9 | 59.4% |
| — of which the **elements' own** `i/G/q/C` | **23.9** | **18.2%** |
| — of which the **stamping loop** | 54.0 | 41.2% |
| solver / LU / control | 53.3 | 40.6% |
| total | 131.2 | |

The phase is 59%, not 5%. But only the element half of it is what a C
backend can touch, and that is **18.2%** -- a ceiling of **1.22x** even
with device evaluation driven to zero.

### The share grows with circuit size, and still does not clear 2x

Four elements could understate it, so it was scaled -- n cascaded
RC-diode stages, all eager classes:

| devices | nodes | total ms | element ms | element share | ceiling |
|---|---|---|---|---|---|
| 4 | 2 | 61.2 | 10.8 | 17.6% | 1.21x |
| 13 | 5 | 190.0 | 55.3 | 29.1% | 1.41x |
| 49 | 17 | 606.9 | 239.3 | 39.4% | 1.65x |
| 121 | 41 | 1338.0 | 579.2 | 43.3% | **1.76x** |

It rises and flattens toward ~1.8x. Against the C backend's measured
**212x** on PSP's Jacobian, an eager C backend is a different kind of
item: same build cost, two orders of magnitude less return. **Do not
build it**, and do not build S1 either -- S1's value was only ever
"batching on top of C", and there is no model that can have both (§22).

### What the profile actually points at

The top cost centres in a 49-element run are not arithmetic:

```
2259 calls  0.174s  circuit.py:1674(_add_element_subvectors)
794299 calls 0.147s  {built-in getattr}
1781 calls  0.143s  circuit.py:1589(_add_element_submatrices)
421670 calls 0.104s  param.py:218(ParameterDict.__getattr__)
125029 calls 0.091s  hdl.py:5092(_args_of)
198058 calls 0.083s  {ndarray.flatten}
125029 calls 0.058s  hdl.py:5017(_params_of)
480149 calls 0.046s  {built-in hasattr}
```

`_args_of` runs **125 029 times** and rebuilds the same list every time:
one `getattr(self.iparv, name)` per parameter plus one `is_given` per
given-flag, on every `i`/`G`/`q`/`C`/`u` call of every element of every
Newton iteration of every timestep. **Parameters do not change during a
transient.** This is marshalling, not physics, and it is pure repeat
work.

### Measured ceiling of caching it

A per-instance cache keyed on `T` (no invalidation -- a ceiling probe,
not an implementation):

| devices | base ms | cached ms | speedup |
|---|---|---|---|
| 13 | 172.5 | 150.1 | 1.15x |
| 49 | 533.4 | 445.0 | 1.20x |
| 121 | 1149.5 | 955.7 | **1.20x** |

Answers **bit-identical** at n=4 and n=16, first and last node,
`max|diff| = 0.0`.

So roughly **15 lines buys 1.20x** -- more than an entire eager C
backend would return on a small circuit (1.21x), and about two-thirds of
what one would return on a large one (1.76x), at a fraction of the cost
and with no second code path to keep in step.

### The one risk, and the hook that answers it

A naive cache goes stale the moment a parameter is reassigned -- and
sweeps, `.alter`, and temperature loops all do exactly that.
`ParameterDict` already has the invalidation point: **`notify`**, called
from both `__setattr__` (`param.py:230`) and `update_values`
(`param.py:189`). The real version must clear the cache from there, not
guess at a key. `T` still belongs in the key, since `epar.T` varies
within a run without touching `ipar`.

Not built -- measured and recorded, per the same rule §22 applied.

### Ranking after this measurement

1. **`_args_of` caching** -- 1.20x, ~15 lines, needs the `notify` hook.
2. **The stamping loop** (`_add_element_subvectors` /
   `_add_element_submatrices`, 41% of the transient and the single
   largest phase) -- unexamined; 198 k `flatten` calls and 480 k
   `hasattr` calls in it suggest the same repeat-work shape.
3. ~~C backend for the eager path~~ -- **refused on measurement**,
   ceiling 1.76x.
4. ~~S1~~ -- **refused on measurement** (§22), no model can hold both
   fast paths.

Two of the four items the roadmap was carrying are now closed by
measurement rather than by building them.

## 24. `_args_of` memoised (2026-08-26): §23's first item, built -- 1.16x

§23 ranked the marshalling cache first and did not build it. Built now.

### What it is

`hdl._args_of` builds the trailing argument list every compiled function
is called with -- parameter values, `T`, givenness flags. It runs once
per `i`/`G`/`q`/`C`/`u` call, so once per element per Newton iteration
per timestep: **125 029 times** in a 49-element transient, rebuilding an
identical list every time. It is now memoised on the instance.

### Invalidation is the whole design

A value cache is worth exactly its invalidation, and this one reads
**two** parameter dicts -- values from `iparv`, givenness flags from
`ipar`. It hangs off the `iparv` observer (`update`), the same one that
already drops the constant stamps `_hdl_Gc`/`_hdl_Cc` and the C
backend's packed vector `_hdl_cp`. That covers both dicts, because the
routes converge:

```
iparv written  -> notify -> update
ipar  written  -> _ipar_changed -> update_iparv
                             -> iparv.update_values -> notify -> update
```

So a sweep, an `.alter`, an `ipar.set()`, an explicit `update_iparv()`
and a late-resolved parameter expression all invalidate. **`T` is in the
KEY, not the cached list** -- `epar.T` varies within a run without
touching either dict, so it must select between entries rather than be
stored beside them. A non-scalar `T` (symbolic, array) is not cached at
all: the key comparison would be ambiguous, and no production path does
it.

**The key is `(value, type)`, not just value** -- found while auditing
rather than by a failure. `defaultepar` carries `T` as an **`int`**
(300, not 300.0), and `300 == 300.0` is True in Python, so a key matched
on value alone would hand a float-`T` caller the list built for int-`T`.
Numerically harmless in every expression this tree evaluates, but it is
a divergence from the uncached path, and the C backend work already
recorded one bug where an integer subchain changed the sign of an exact
zero. Comparing `type()` as well costs one `is` per hit against
rebuilding the whole list. `test_an_int_temperature_does_not_answer_a_float_one`
pins it and fails without the type check.

Two further exposures checked and closed: **no caller mutates the
returned list** (all nine sites unpack it or slice a copy -- audited, and the
docstring now says a mutating caller must copy first); and **class
retargeting for node collapse cannot outlive an entry**, because it runs
in `__init__` before any evaluation and the variant carries identical
`instparams`.

### Measured

`benchmarks/args_cache.py`, shipped implementation against an uncached
twin, one process, interleaved, three repeats, minimum of each:

| devices | uncached ms | cached ms | speedup |
|---|---|---|---|
| 13 | 171.5 | 155.3 | 1.10x |
| 49 | 531.1 | 464.5 | 1.14x |
| 121 | 1159.0 | 1002.5 | **1.16x** |

Answers **bit-identical** at every size and repeat. The element phase
fell from 43.3% to 36.3% of a 121-device transient -- the win came out
of exactly the phase §23 predicted.

⚠ **THIS NUMBER WAS PUBLISHED WRONG TWICE BEFORE LANDING AT 1.16x**,
which is worth more than the number itself:

| figure | where it came from | why it was too high |
|---|---|---|
| 1.20x | §23's ceiling probe | cached with NO invalidation -- cheaper than the real thing |
| 1.18x | first interleaved A/B | taken before the key was made type-exact |
| **1.16x** | final, machine idle | -- |

Each was caught by re-running rather than by carrying the earlier figure
forward. That is the rule from §11's `both = 30` incident applied three
times in one afternoon, and it earned its keep every time: **every
intermediate measurement in this work was optimistic, and none of them
was wrong by enough to look wrong.**

📌 **AND 1.16x DID NOT STAY 1.16x EITHER -- the code did not change.**
Re-run after §26, the same benchmark on the same unchanged cache reads
**1.16x / 1.22x / 1.24x**. §26 removed a 425 ns `epar` lookup from the
same method, and taking overhead out of the denominator raises the share
this cache accounts for. **A speedup figure is a property of the whole
call, not of the change alone**, so it expires whenever anything else in
that call moves -- a fourth way for a number to go stale, and the only
one where re-running the *same* benchmark on the *same* code gives a
different answer. `benchmarks/args_cache.py` now says to re-run rather
than quote.

### Tests

`pycircuit/circuit/tests/test_hdl_argcache.py`, nine tests. The one that
gives the others meaning is `test_the_cache_is_actually_used` -- every
other test in the file passes with the memoisation deleted, since they
would then only be re-testing that `_args_of` reads live parameters.

**The invalidation tests were verified to BITE**: with the
`_hdl_args` pop removed from `update()`, five of the nine fail,
including a parameter sweep frozen at `[0.01, 0.01, 0.01, 0.01]` --
a resistor that ignores every value after the first.

### It also caught a stale number in the user-facing docs

`hdl.rst` claimed the DSL's end-to-end overhead was **1.14x** the
hand-written elements on an 8-stage RC ladder. Re-measured on the same
machine at HEAD, before this change, it was **1.22x** -- the figure had
gone stale at some point in the campaign and nobody re-ran it. With the
cache the same benchmark reads **1.15x** -- most of the drift recovered,
and close to the figure the docs had been claiming all along.

| 8-stage RC ladder | ratio, 4 runs | median |
|---|---|---|
| at HEAD (no cache) | 1.22 / 1.22 / 1.19 / 1.23 | **1.22x** |
| with the cache | 1.16 / 1.15 / 1.15 / 1.16 | **1.15x** |

⚠ Single runs of this benchmark said **1.25x** and **1.12x** -- the two
extremes -- and quoting them would have overstated the change in both
directions. Four runs each is what the numbers above rest on.

`hdl.rst` now says 1.15x and carries a note recording that it read 1.14x
until 2026-08-26 and had drifted to 1.22x. **The published figure was
wrong in the direction that flatters the DSL**, which is the direction
nobody checks -- the same failure mode as roadmap §11's `both = 30`, and
the second time in this campaign a number went stale in prose while the
test around it stayed green.

Worth noting what this says about the cache's reach: the ladder is 17
elements of ONE parameter each, the least favourable case there is, and
it still gains ~1.10x in a direct A/B. The win is not in the arithmetic
the parameters feed -- it is in `ParameterDict.__getattr__` and the
attribute protocol, which cost the same whether a model has one
parameter or forty.

**A timing assertion was considered and NOT added**, and the reason is
recorded rather than left implicit. A bound loose enough not to flake on
a shared machine (say "under 1.6x") would not have caught the drift this
work found -- 1.22x sails under it -- and a bound tight enough to catch
it would fail on a busy runner. The guard is STRUCTURAL instead:
`test_the_cache_is_actually_used` asserts the memoisation is live, which
is the property whose loss caused the regression, and it cannot flake.
The `record-upkeep` skill's new rule ("re-run the published figures that
measure a path you touched") covers the rest, because that is a
discipline question and not a test-shaped one.

### Ranking after this

1. **The stamping loop** -- `_add_element_subvectors` /
   `_add_element_submatrices`, still the largest single phase (41%),
   still unexamined, 198 k `flatten` and 480 k `hasattr` calls.
2. ~~`_args_of` caching~~ -- **DONE**, 1.16x.
3. ~~C backend for the eager path~~ -- refused on measurement (§23).
4. ~~S1~~ -- refused on measurement (§22).

## 25. The stamping loop, measured and reduced (2026-08-26): `flatten` was copying for nothing

§23 left the stamping loop as the largest unexamined phase (41%). Opened.

### Where a stamp call actually goes

One `cir.i()` on 49 elements, phases timed by replicating each in
isolation (accounted 95.6% of the call):

| phase | us | share |
|---|---|---|
| element loop TOTAL | 109.93 | 87.2% |
| — of which the element's own work | 77.99 | 61.8% |
| — of which **loop overhead** | **31.93** | **25.3%** |
| `_scatter_1d` | 9.20 | 7.3% |
| prologue (`zeros` + probes) | 1.42 | 1.1% |
| `batched_contributions` (no groups) | 0.05 | 0.0% |

So the loop is 25.3% overhead at 651 ns per element, against 1637 ns of
real work. The scatter -- the part STAGE 2b rewrote -- is already down to
7.3% and is not where the remaining time is.

### The 651 ns, itemised

| operation | ns | per |
|---|---|---|
| `np.asarray(rhs).flatten()` | 397 | element |
| `x[nodemap]` | 135 | element |
| `idxmap.get(instance)` | 57 | element |
| `getattr(element, meth)` | 41 | element |

`flatten` is **61% of the overhead**, and it is waste: **`flatten` always
copies**, while a generated 2-D stamp is already contiguous. Measured on
the real G stamp, `np.asarray(rhs).flatten()` is **412 ns** against
`ravel`'s **100 ns** -- a 4x difference, because `ravel` returns a view.

The copy buys nothing, because the only consumer is the
`np.concatenate` inside `_scatter_1d`/`_scatter_2d`, which allocates its
own buffer regardless. Three sites changed (the two pending appends and
the `build_sparse` append); the two `has_add_at` sites are the traced
path and were left alone.

### Measured

| devices | flatten | ravel | speedup |
|---|---|---|---|
| 13 | 154.3 ms | 143.5 ms | 1.075x |
| 49 | 463.3 ms | 421.2 ms | 1.100x |
| 121 | 997.6 ms | 895.8 ms | **1.114x** |

**Cumulative with §24's `_args_of` cache**, against 62b6a3b:

| devices | 62b6a3b | now | speedup |
|---|---|---|---|
| 13 | 168.9 ms | 144.2 ms | 1.17x |
| 49 | 524.8 ms | 423.4 ms | 1.24x |
| 121 | 1145.2 ms | 903.0 ms | **1.27x** |

Answers **bit-identical to 62b6a3b** at every size, `max|diff| = 0.0`.

### This one is NOT a DSL win -- it is a simulator-wide one

Worth stating plainly, because the two changes in §24 and §25 are
different in kind. `_args_of` is DSL machinery: only generated elements
pay it, so caching it narrowed the DSL's overhead RATIO (the 8-stage
ladder went 1.22x -> 1.15x). The stamping loop is shared: **hand-written
elements go through the identical assembly**, so `ravel` speeds them up
too and the ratio does not move.

On the `hdl_overhead.py` ladder the hand-written side went from ~0.63-0.69 s
to ~0.58 s while the generated side went from ~0.72 s to ~0.67 s -- both
about 10% faster, ratio unchanged at ~1.15x. **A ratio that stays put is
not evidence that nothing happened**, and reading it that way would have
buried this result: the ratio is blind to any change that helps both
sides equally, which is exactly what a shared-path fix is.

### What was NOT taken, and why

`hasattr(toolkit, 'add_at')` costs **1024 ns** -- it misses, so it enters
`Toolkit.__getattr__`, which formats an error message and raises.
Counted per transient: **4624 failed lookups on 13 devices** (~4.6 ms,
3.2% of that run) and 5823 on 49 devices (1.4%). `add_at` and
`build_sparse` are the whole of it.

The obvious fix is to memoise failures alongside successes -- and
`toolkit.py` **records that decision already, and decided against it**:
a negative cache "would also make an attribute that appears later
permanently invisible", and the hot probes were hoisted out of the
per-element loop in STAGE 2b instead. That hoist is why the count is now
per stamp CALL rather than per element.

**Left alone.** The residual is 1.4-3.2%, it does not grow with circuit
size (it is per call, not per element), and reversing a recorded design
decision for that is the user's call, not a silent refactor. Recorded
here with the numbers so the decision can be revisited on evidence
rather than rediscovered.

### Tests

`pycircuit/circuit/tests/test_stamp_assembly.py`, 7 tests. `ravel`
trades a copy for a view, so the property they pin is that **nothing
observable aliases**: the assembled matrix shares memory with no
element stamp, writing to the result cannot corrupt an element, a cached
constant stamp survives repeated assembly, and the assembled `G` equals
a hand-built sum (not a golden number).

Honest note on what bites: the aliasing tests pass under BOTH `flatten`
and `ravel` -- they are regression guards for a future change that
removes the copy `concatenate` provides, not proof of this one. The
hand-sum test is the one that discriminates, and it was verified to fail
on a deliberately transposed scatter.

### Ranking after this

1. **Element evaluation itself** -- now 61.8% of a stamp call and the
   largest remaining phase by far. The eager path returns a Python
   **list**, so `np.asarray` converts it on every call; the chained path
   already returns an array. Unmeasured.
2. `hasattr` on the toolkit -- 1.4-3.2%, blocked on the recorded
   decision above.
3. ~~stamping loop~~ -- **DONE**, 1.11x.
4. ~~`_args_of` caching~~ -- **DONE**, 1.16x.
5. ~~C backend for the eager path~~ / ~~S1~~ -- refused on measurement.

## 26. Element evaluation (2026-08-26): two thirds of it was never the model

§25 left element evaluation as the largest phase (61.8% of a stamp
call). Opened, and the answer is that most of it was not evaluation.

### The generated function is a fifth of what the method costs

| element | `el.i()` | the generated `f(x, *args)` | wrapper |
|---|---|---|---|
| `RHdl` (eager) | 1236 ns | 277 ns | 959 ns |
| `DiodeHdl` (eager) | 1439 ns | 473 ns | 966 ns |
| `EkvNmosHdl` (chained) | 14776 ns | 13058 ns | 1718 ns |

For a small eager element **the arithmetic is 22% of the call**. The
rest is the method around it.

⚠ **A confident wrong guess, recorded because it was nearly acted on.**
959 ns is almost exactly the 1024 ns that a *missing* toolkit attribute
costs (§25), and `i()` contains `getattr(self.toolkit, 'symbolic',
False)` -- so the obvious reading was that `symbolic` misses and every
element pays `Toolkit.__getattr__`'s raise. Measured before touching it:
`symbolic` **resolves normally, in 39 ns**. The shape of a number
matching a known mechanism is not evidence that it IS that mechanism.

### Where it actually went

| operation | ns | note |
|---|---|---|
| `_dc(epar)` | 425 | `getattr(epar, 'analysis_kind', ...)` on a `ParameterDict` |
| `_args_of` (cache HIT) | 314 | the hit path itself, not the rebuild |
| the generated `f` | 279 | the work |
| `getattr(toolkit,'symbolic')` | 39 | |
| `f.__dict__.get('_hdl_c')` | 38 | |
| `info['chained']` | 28 | |

### Fix 1 -- short-circuit on the cheap operand (`has_dc_pins`)

Every site read `_dc(epar) and state_meta['dc_pins']`. Python evaluates
the **left** operand first, so all 37 classes paid the 425 ns lookup and
then discarded it -- and **only 3 of the 37 have DC pins at all**
(`IdtmodHdl`, `MemristorHdl`, `VcoHdl`). Both operands are used for
truthiness only at all five sites, so the order is free: hoisting
`has_dc_pins = bool(state_meta['dc_pins'])` and testing it first means
34 of 37 classes never call `_dc`.

`el.i()`: **1236 -> 826 ns.**

### Fix 2 -- identity before value in the `_args_of` key

The §24 cache still cost 314 ns *on a hit*: an `isinstance` against a
4-tuple, an `==`, and a `type()`. But `epar.T` is read from the same
`ParameterDict` slot every time -- instrumented over a real transient,
**one `T` object across 60 924 calls**. So an `is` settles the common
case in a single pointer compare, placed *in front of* the value test
rather than replacing it: a caller that rebuilds `T` merely falls
through and pays what it used to. Strictly safer, never weaker.

`_args_of` **314 -> 208 ns**; `el.i()` **826 -> 648 ns**, against 279 ns
of arithmetic.

### Measured, cumulative

| devices | 62b6a3b | now | speedup |
|---|---|---|---|
| 13 | 168.9 ms | 136.1 ms | 1.24x |
| 49 | 524.8 ms | 388.0 ms | 1.35x |
| 121 | 1145.2 ms | 816.4 ms | **1.40x** |

Answers **bit-identical to 62b6a3b**, `max|diff| = 0.0`. `el.i()` itself
is **1.91x** (1236 -> 648 ns).

### The DSL's end-to-end overhead is now PARITY

The figure `hdl.rst` exists to publish -- generated elements against
hand-written ones on the 8-stage ladder -- has closed:

| | ratio |
|---|---|
| documented (stale) | 1.14x |
| actually, at HEAD | 1.22x |
| after §24 (`_args_of`) | 1.15x |
| after §25 (`ravel`) | 1.15x (shared path -- helps both sides) |
| **after §26** | **1.01x** (0.99-1.03 over five idle runs) |

Per call the generated elements are now **faster than hand-written
wherever there is real arithmetic** -- the diode's `i` is **0.47x**, its
`G` 0.93x, the capacitor's `C` 0.62x. What is left above 1.0 is the
trivial linear stamp (`R G` at 2.55x), where the hand-written element
does almost nothing and any call overhead dominates a 2x2 matrix --
which is what `hdl.rst` already said, and it is still the right
explanation.

⚠ Measured on an **idle** machine. The same benchmark run while the test
suite was going gave 0.97-1.05x, a spread wide enough to read as
"faster than hand-written" if a single run were quoted.

### A test that did not discriminate, caught before it was believed

The first version of `test_the_identity_fast_path_is_used` asserted that
two calls return the same list. **It passed with the fast path deleted**
-- the value+type path returns the cached object too, so the assertion
was blind to the thing it named.

Rewritten on **NaN**, which separates the two exactly: `nan is nan` is
True for one object while `nan == nan` is False, so the entry can only be
reused through the identity check. Verified to fail when the path is
removed.

This is the fourth test in this campaign that had to be re-aimed after
being written, and the pattern is always the same: **the assertion was
true of the code, but would have been true without the feature too.**
Deleting the feature and re-running is the only check that catches it.

### Ranking after this

1. **`_dc` for the 3 DC-pin classes** and the `_args_of` hit's residual
   208 ns -- both now small; no measurement demands them.
2. `hasattr` on the toolkit -- 1.4-3.2%, blocked on the recorded
   decision in `toolkit.py` (§25).
3. The eager path returns a Python **list**, so `np.asarray` converts it
   in the stamping loop every call. Moving the conversion into the
   generated code does not remove it -- unmeasured whether anything can.
4. ~~element evaluation~~ / ~~stamping loop~~ / ~~`_args_of`~~ -- DONE,
   **1.40x cumulative**.
5. ~~C backend for the eager path~~ / ~~S1~~ -- refused on measurement.

## 27. The damped `predict`, built (2026-08-26): and the failure was mostly the SEED

§15 left "a native rescue for PCNR from a wild start, which would need a
damped `predict`". Built. The damping works, is off by default, and is
**not** what fixed most of the problem.

### Instrumenting the first step, before designing anything

BJT mirror, uniform 20 V start, iteration by iteration:

```
start x       : [20, 20, 20, 0, 20]
v_lim init    : [20, 0, 20, 0]        <-- vbe seeded at 20 V
it 0  cond(S)=4.6e+94  max|g_mna|=1.6e+92  max|dx_mna|=4.5e+45
      after refine: v = [0.849, 0, 0.849, 2.811]
```

**`v_lim_init` seeded the unknowns with the RAW branch voltages.** So the
very first `augmented_system` evaluated the transport current at
`exp(20/vt)` and the first step was computed from a Jacobian with
condition number 1e94. `refine` then pulled `v` back to a sane 0.849 V --
one step too late, with `x` already at -4.5e45.

That is also **why no ladder could ever have worked** (§25): every rung
began by building the same Jacobian at the same unlimited seed.

### Fix 1 -- limit the seed (this is the one that mattered)

`v_lim_init(..., limit=True)` passes the seed through each device's own
law, exactly as `refine` does for every later iterate.

| start | mirror before | after |
|---|---|---|
| 0 V | 10 | 10 |
| −5 V | 10 | 10 |
| +5 V | 165 | **9** |
| +20 V | **LinAlgError** | **8** |

It is **inert where the start is sane** -- seeds up to ~0.8 V come back
unchanged on the Gummel-Poon card; 1 V → 0.094, 5 V → 0.136, 20 V →
0.172. And it is scoped: **`solve_dc` clamps, the transient does not.**
The transient seeds from the previous time point's *accepted solution*,
which is a real operating point and the best information available;
clamping it moved the rectifier's waveform 1.7e-8 V and its time points
6.5e-12 s against the classic-limiting path it is required to track.

### Fix 2 -- `gmin` on the Schur diagonal, OFF by default

Two cheaper dampers were measured first and **both are wrong**:

| damper | EKV diff pair | Level-1 cascode | grid (2,2,2) |
|---|---|---|---|
| truncate the step (`cap/max|dx|`) | fixed, 10 V…1e6 V | **broken at every cap** | -- |
| backtrack on the residual norm | **not fixed** | 9 → 5 | -- |
| `gmin` on the Schur diagonal | fixed | fixed | **broken, 1e-12…1e-6** |

A magnitude cap is not a safe damper: **a huge step can be better than a
truncated one.** A line search controls the step but cannot condition the
matrix -- the step it accepts is still one whose *next* Jacobian is
singular. Only regularising the matrix addresses the actual defect, which
is a **rank deficiency**: the EKV pair from a wild start measures the
Schur complement at **rank 8 of 9**, both channels cut off, 2e-19 S to
the tail. A rank-deficient solve does not return a large step, it returns
a meaningless one (3e48 V).

The anchor goes on the **Jacobian only, never the residual**, so the
converged point solves the untouched system -- matched to `DC()` at
1e-11 on all three circuits.

**It still defaults to off**, because it trades circuits: it rescues the
EKV pair and breaks cascode grid point (2,2,2) at *every* value swept.
"Nothing that converged stops converging" outranks rescuing a start the
fallback already handles. It is also non-monotone -- the Level-1 cascode
fails at 1e-14/1e-12/1e-11, converges from 1e-10 up, **and converges at
exactly zero** -- so an iteration count from a wild start is not a stable
acceptance number.

### Measured, everything

Every PCNR case, both instance orders, answers matched against `DC()`:

```
BJT mirror        0 V [10,10]   +5 V [9,9]   +20 V [8,8]   −5 V [10,10]
EKV diff pair     0 V [8,8]     +20 V [9,9]   (with gmin=1e-9)
Level-1 cascode   0 V [7,7]     +20 V [9,9]   all four grid points
```

Counts equal in both orders everywhere: **order-independence preserved**,
which is what the fallback costs and this does not. And at the +5 V start
**PCNR now converges in 9 where the ordinary chain -- four ladders plus
the anchor -- fails outright.**

### ⚠ The recorded diagnosis of the cascode was wrong

`test_the_level_1_cascode_from_20_v_under_pcnr` explained its `[FAIL,
FAIL]` row as structural: *"the middle node is held by nothing... PCNR
limits no node... which is Stage 1's second obstacle on a FET."* Every
observation in that is real. The conclusion is not.

**The runaway was seeded, not structural.** `vgs = 20 V` went into the
first Jacobian, the 4.9e9 first step came out of it, and once the seed is
limited the node never leaves its basin -- with no node write involved.
PCNR still writes no node; that was never what made this fail.
**Eighteenth right-conclusion-wrong-reason**, and the docstring keeps the
wrong version beside the right one.  (Written as "sixteenth" when §27
landed; §20 had already used sixteenth and seventeenth, so the count had
gone backwards -- corrected on the next read of the document.)

### Three tests inverted, one re-aimed

- the mirror's **claim 3** (*"PCNR fails identically in both orders at
  +20 V"*) -- written as an equality *specifically* so a later fix would
  pass it unchanged, which it did. Now also asserts both CONVERGE, so a
  regression to `[FAIL, FAIL]` cannot pass.
- the cascode's fourth row -- was asserted `is None`.
- **the fallback test had to move circuits.** Its precondition was that
  PCNR fails on the mirror; that stopped being true. A precondition that
  stops holding makes a test **vacuous**, and this one failed loudly
  (`DID NOT RAISE`) rather than passing quietly -- the good outcome. It
  now runs on the EKV pair, which still fails, and for a different
  reason.

## 28. The `select` residue, surveyed (2026-08-26): §21's "adoption work" was wrong

§21 converted 6 of 20 surveyed sites and recorded that the remaining 14
were **"adoption work, not a limitation -- the mechanism covers their
shapes."** Went to do that adoption. **The claim does not hold**, and the
reason is structural rather than incidental.

### What `select` can substitute, and what this codebase writes

`select` replaces a **symbol appearing inside an arm expression** with a
clamped copy derived from that arm's condition. So it needs the arm to
*contain* the symbol.

This codebase almost never writes them that way. Its house style binds
every intermediate with `var()` first, so the arm is a **bare symbol**
and any guard lives in the referenced definition, evaluated
unconditionally somewhere above:

```python
qlo = _var(cj * vjx * (1.0 - ub) / (1.0 - mc), 'qlo')   # guard is in `ub`
qhi = _var(cj * (f1 + ...), 'qhi')
return _var(sympy.Piecewise((qlo, v < fcvj), (qhi, True)), 'qj')
```

There is nothing in `(qlo, v < fcvj)` for `select` to clamp. The same
shape covers the depletion charges, `ids`, `isbd`/`isbs`, `cbdb`/`cbsb`,
`rdx`/`rsx`, `rbmx`, `vp`, `vgtn`, `vdsn`, `kp0` and the Biolek `stp`.

Two further shapes defeat it for different reasons:

| shape | example | why |
|---|---|---|
| clamp **hoisted** into its own `var` | `xjf = maxc(p.xj, 1e-12)`, used by `fshort`; `vdsc = leff*maxc(p.vmax,1e-30)/us`, used by `fdrain` | the clamp is not in the arm and is shared by other users |
| clamp on a **different symbol** than the condition bounds | `wfact`: condition is `vgs < von`, clamp is `maxc(xn, 1.0)` | the condition says nothing about `xn` |
| **compound** argument | `_recip(p.area * p.ikf)` | already recorded at the site: the product is not a node, so the substitution silently does nothing |

**Exactly one of the residue was convertible: `alpha` in MOS level 3**,
where the clamp `maxc(nsm3, 1.0)` sits inside the arm and the condition
is `nsm3 > 0.0`. Converted, with `margin=1.0` reproducing the hand clamp
(`nsm3` is a doping density near 1e21, so a 1.0 floor is far below
anything selected). `elements_hdl.py` is now 7 `select` against 32
`Piecewise`; `psp_kernel.py`'s 24 remain deliberately untouched.

**The honest status is therefore: `select` is done, and its adoption is
done.** The remaining `Piecewise` sites are not waiting on anyone --
they are a different shape, and converting them would mean first
un-binding intermediates that `var()` exists to bind. That trade
(§21 already noted it for `_recip`: "it also changes the emitted chain")
is not obviously worth making, and no measurement asks for it.

### ⚠ And the verification method §21 accepted is BLIND to this class of change

§21's evidence was: *"all 37 library classes digested over 30 random
biases x {i, G, q, C} with the adoption in place and again reverted --
37 identical, 0 differ."* The same recipe was run here, and it reported
**37 identical** for the `alpha` conversion too.

Then the margin was deliberately set **wrong** -- `1e-30` instead of
`1.0`, which is not the hand-written clamp at all -- and the digest still
came back **bit-identical across all 37 classes**.

The cause: `nsm3 = maxc(nsub, 0) * 1e6` and **`nsub` defaults to 0**, so
at default cards the `nsm3 > 0` arm is never selected and the margin
cannot matter. A digest over default parameters cannot see a change in
an arm that default parameters do not select. `validation-design`
already carries the rule -- *a bias no sweep visits tests nothing* -- and
the digest walked into it.

Validated instead on a card that reaches the arm (`nsub = 1e16`, plus
the geometry the level-3 model needs): 1600 samples of `i`/`G`/`q`/`C`
over 40 random biases, **bit-identical** against the pre-change source.
Pinned by `test_select_with_a_margin_reproduces_the_hand_clamped_piecewise`,
which compiles both forms and compares value AND derivative across the
selected region, the boundary and outside it -- and by
`test_the_mos3_alpha_arm_is_reached_only_with_nsub_set`, which exists to
record why the default card is not enough.

**And it is worse than one blind spot: the evidence for §21's SIX
conversions was vacuous too.** Tested directly rather than left as a
suspicion -- `select` was reduced to a plain `sympy.Piecewise`, so that
**no clamping happened anywhere in the library**, and the digest was
re-run:

```
select reduced to a plain Piecewise (NO clamping anywhere)
  -> NO DIFFERENCE across all 37 classes
```

So "37 identical, 0 differ" would have been reported whether the
clamping worked or not. It was never evidence for the adoption.

**Why**, and this is the transferable part: a clamp is **insurance
against the tail**. It only changes a number where an arm would
otherwise be non-finite or wildly out of range -- and a digest over 30
random biases in a normal window samples the **bulk**, where every arm
is finite already and the clamp is the identity by construction. The
method could not have worked.

**This does NOT mean the conversions are wrong.** `select` itself is
well covered: `test_hdl_select.py` exercises the mechanism directly, and
neutralising the margin fails **19** of its tests. What was vacuous is
the *adoption* evidence -- the step that claimed these particular six
sites still compute what they used to. Re-running it properly means, for
each converted arm, a card and a bias that **select that arm and put it
near its boundary**, which is where a clamp can differ from the
expression it replaced.

## 29. S2 measured before building it (2026-08-26): the inert blocks are CHEAP

§3 proposed S2 -- generalise the per-parameter-mask class-variant
mechanism from "collapse nodes" to "elide any block whose parameters
make it inert" -- and called it *"the most direct attack on evaluation
cost: a card that switches off gate current, avalanche and half the
parasitics should not pay for them on every Newton iteration."*

Weeks of work, medium risk. Measured first, as §22 and §23 were.

### The measurement is available without building anything

`psp_kernel` **already contains the elision.** Its build-time guards fire
when a parameter arrives as a NUMBER rather than a symbol -- and today it
never is, which is precisely the trap S2 describes. Handing the kernel a
literal `0.0` for the switches a card sets to zero produces exactly the
code S2 would emit, so the two can be compiled side by side and compared.

(The compile cache must be off: the variant's `analog` source and
`instparams` are identical to the base's, so it hashes to the same key.)

**And the default card is already the interesting case.** `alp`, `alp1`,
`alp2`, `kp`, `rs`, `rsg`, `rsb`, `gc3`, `gco`, `iginv`, `igov` **all
default to 0.0** -- channel-length modulation, series resistance and gate
tunnelling are all switched off on a default build and all fully
compiled.

### Result: nothing, on either backend

| elided | G source | numpy | C backend |
|---|---|---|---|
| the `mob` switches (`alp`, `kp`, `rs`, …) | 2675 → 2469 lines | **0.996x** | **1.000x** |
| the whole **gate-current** block | 2675 → 2275 lines | **1.000x** | -- |

Four interleaved runs each, minimum taken; spreads of 17.13–17.47 ms.
Every comparison **bit-identical** in `i` and `G` (`max|diff| = 0.0`),
which also confirms the kernel's guards are correct.

Removing **15% of the generated Jacobian source**, exponentials included,
changes the evaluation time by nothing measurable.

### Why -- and the mechanism is the interesting part

Executed primitive calls, counted per evaluation, base against
gate-block-removed:

```
base      total/eval 27659   maxc=2911 minc=3833 rdiv=1605 recip2=1084 step=2046
no gate   total/eval 27659   maxc=2911 minc=3833 rdiv=1605 recip2=1084 step=2046
```

**Identical to the last call.** The 400 removed lines contain **zero**
calls to the smoothing primitives, which are what the cost is made of:

| primitive | us/call | calls/eval | ms/eval |
|---|---|---|---|
| `minc` | 0.673 | 3833 | 2.578 |
| `maxc` | 0.679 | 2911 | 1.977 |
| `step` / `rdiv` / `recip2` | ~0.045 | 4735 | 0.205 |
| **total** | | **11 479** | **4.76 = 28% of the 17.2 ms** |

The optional physics blocks are plain array arithmetic -- C-level ufunc
work that the call counter does not even see. The 845 emitted numpy
calls the gate block costs are worth roughly **0.25 ms of 17.2 ms, 1.5%**,
which is inside the run-to-run spread. So the null is not a measurement
failure; it is the right answer.

**The cost of a compact model is its REGULARISERS, not its physics.**
`minc` and `maxc` alone -- the domain clamping that keeps both arms of
every conditional finite -- are 4.5 ms of the 17.2 ms Jacobian, and they
are spread through every block, including the ones no card switches off.

### Verdict

**S2 is REFUSED on measurement.** Its stated justification was evaluation
cost, and eliding the blocks it targets buys 1.00x on numpy and 1.00x in
C, on the card where the most blocks are inert.

What remains true is the **correctness** half of §3's argument -- *a
block whose parameter is zero is still compiled and still evaluated, and
`0 x inf = NaN`*. That is a real trap and it is why the nine hand guards
in `psp_kernel.py` exist. But hand guards are nine lines, and S2 is weeks
with a variant-explosion risk. If the trap bites again, the cheap answer
is another hand guard, or a checked helper -- not a compiler feature.

Fifth item closed by measuring instead of building (§22 S1, §23
C-for-eager, §25 toolkit memoisation, §28 the `select` residue, this).

### What the measurement DOES point at

If PSP's Jacobian is ever worth attacking again on the numpy path, the
target is now named and quantified: **`minc`/`maxc`, 6744 calls and
4.5 ms per evaluation**. They are Python-level functions called once per
clamped subexpression. The C backend already removes them entirely
(212x), which is why this is recorded rather than pursued.

## 30. Card parameters as compile-time constants (2026-08-26): 1.54x, and it is what S2 was reaching for

Asked directly: *can evaluation be faster if every model-card parameter
— everything not set at instantiation — is treated as a fixed constant?*
Measured. **Yes: 1.54x on PSP's Jacobian**, to machine precision in the
regime the model is validated in.

This is the first item this session that measured **well**.

### Two different mechanisms, and only one of them works

| | what it removes | measured |
|---|---|---|
| **§29 S2** — elide blocks a card makes inert | whole blocks, 15% of the source | **1.00x** |
| **§30 folding** — make every card parameter a literal | parameter-only arithmetic *throughout* | **1.54x** |

§29 established that a compact model's cost is its **regularisers**, not
its physics. Folding wins precisely because it attacks those: the
`minc`/`maxc`/`_step`/`_rdiv` calls in the emitted source fall from
**11 565 to 7 251, −37%**, because a clamp whose bound is now a literal
often folds away entirely. Eliding blocks removed **zero** of them.

**And folding subsumes S2**: with the parameters numeric, `psp_kernel`'s
build-time guards fire on their own — which is exactly the trick §29 used
to simulate S2 in the first place.

### Measured

Real IHP sg13g2 card (`cornerMOSlv.lib`, `mos_tt`, w = 10 µm, l = 1 µm),
68 of 153 parameters supplied, compile cache off:

```
                  G source     primitive calls    G evaluation
symbolic params   2675 lines        11 565        17.90 ms
folded params     1971 lines         7 251        11.62 ms   1.54x
```

Four interleaved runs; spreads 17.90–18.12 and 11.61–11.64 ms.

### Accuracy — and why the first number was misleading

Folding reassociates arithmetic, so this is **not** bit-identical. What
matters is where the difference lands:

| observable | max relative difference |
|---|---|
| `G`, at the measured bias | 1.4e-13 |
| drain current, strong inversion (>1e-6 A), 120 points | **8.5e-16** |
| drain current, validated window (1e-9…1e-6 A) | **2.7e-15** |
| `i` at 40 **random** biases on all six nodes | **1.6e-6** |

⚠ The random-bias figure was measured **first**, and taken at face value
it would have killed the idea: 1.6e-6 is the size of PSP's entire
vendor-validation budget (1.3e-6 at the worst of twelve sweeps). Banding
the same comparison by current magnitude on a real gate sweep shows the
validated regime agrees to **machine precision**, and the 1.6e-6 comes
from bias combinations no sweep visits — random values on internal nodes
included. **A tolerance is meaningless without the window it applies
in**, which is `validation-design`'s own rule, and it nearly cost a 1.54x
result here.

The random-bias regime is **not characterised** and would need to be
before adopting this — it is recorded as an open question, not waved
away.

### What it would cost to build

A compiled variant **per distinct model card**, which is the same
variant-explosion risk §3 named for S2 — but with a better ratio: one
card typically serves every transistor of that flavour in a design, so
the variant count is the number of cards, not the number of masks.

The two things that made it expensive when §3 was written are now in
place: the **on-disk compile cache** (§16) makes the ~30 s build a
one-time cost per card, and `_collapse_variant` already does
compile-a-class-per-key and retarget-by-`__class__`, so the machinery
exists and would be keyed on the card instead of a mask.

**Unmeasured:** whether this stacks on the **C backend**. C already
removes the Python-level primitives entirely (212x), so folding's −37%
of them may buy much less there — and the C path is where a large model
would actually run.

### Status

**Not built.** It is the strongest remaining performance item by a wide
margin, it has a real accuracy question attached, and the C-backend
interaction should be measured before anyone commits days to it.
`benchmarks/card_constant_folding.py` reproduces all of the above.

## 31. `hdl.fold_card` BUILT (2026-08-26): 1.54x numpy, 1.48x C

§30 measured card-constant folding and left it unbuilt pending one
question -- *does it stack on the C backend, which already removes the
Python primitives entirely?* Measured first: **it does, 1.44x**. Built.

### The API, and the design error the measurement caught

```python
Nmos = hdl.fold_card(PspMosLongChannel, instance=('w', 'l'), **card)
m1 = Nmos('d', 'g', 's', 'b', w=10e-6, l=1e-6)
```

**Everything not named in `instance` is folded** -- at its `card` value
where one is given, and at its **declared default** otherwise.

The first version took only `**card` and folded exactly what was passed.
That left **87 of PSP's 153 parameters symbolic**, almost nothing folded,
and measured **0.96x -- a slowdown**. A parameter sitting at its default
is just as constant as one the card sets. The API had to name what stays
*variable*, not what is fixed, and that is what the user's own framing
said: *"all parameters that are not part of the instantiation."*

Pinned by `test_everything_outside_instance_is_folded`.

### Measured, shipped API

```
                  G source     primitives    compile
symbolic          2675 lines     11 565      ~28 s
folded            1971 lines      7 251      17.2 s

numpy   symbolic 17.73 ms   folded 11.55 ms   1.535x
C       symbolic  83.5 us   folded  56.5 us   1.478x
        max relative dG 1.44e-13
```

It **compiles faster too** -- 17.2 s against ~28 s -- because there is
less code to print and lambdify. The C build halves for the same reason.

### Three edges, each with a test

**It is not bit-identical, and cannot be** (§30): folding reassociates
arithmetic. Machine precision in PSP's validated window, 1.6e-6 at
arbitrary internal-node biases. A model validated against a reference
must be **re-validated with its card folded** -- said in the docstring,
in `hdl.rst`, and in the test module's header.

**A folded parameter cannot be changed.** Its value is in the compiled
code, not in `iparv`, so setting one on an instance would give an
element that *reports* the new value and *evaluates* the old. It raises,
naming `fold_card`. The folded parameters' **defaults move to the card's
values** so an instance reports what it actually evaluates -- which is
also what makes the compile-cache keys differ (see below).

**Folding at bare defaults is REFUSED.** A model guards `1/p` for the
`p` it expects to be zero; a folded card evaluates that quotient at
build time, and a parameter left at `0.0` makes it `ComplexInfinity` --
which surfaces from inside sympy's printer as `KeyError:
'ComplexInfinity'`, saying nothing. Caught and renamed. **Defaults are
not a physical card**, and both PSP and the Gummel-Poon demonstrate it.

### Two of my own tests did not discriminate, and the bite check caught both

- *"the instance parameter stays variable"* first evaluated at two values
  of `bf` and required the currents to differ. They did not -- at a
  saturated operating point the terminal currents barely depend on the
  forward beta. **The test was testing my choice of bias.** Rewritten
  structurally: `bf` must still be a symbol in the emitted body.
- *"two cards do not share a cache entry"* passed with the `fold=` key
  entry **neutralised**. What actually separates the keys is the
  rewritten **defaults**, which `_param_decl` already hashes. The test
  now pins that mechanism, and says plainly that `fold=` is defence in
  depth.

Both found by deleting the feature and re-running -- `validation-design`'s
rule, earning its place for the third time today.

### Status

Suite **2657 passed / 7 skipped / 3 xfailed / 0 failed**, reconciled as
2663 collected + 4 module-level `importorskip`. `sphinx -W` clean.
`benchmarks/card_constant_folding.py` reproduces the numbers;
`hdl.rst` documents it with the warning attached.
