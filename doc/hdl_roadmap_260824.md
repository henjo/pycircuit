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

`_collapse_mask_of` / `_collapse_variant` already
recompile-and-cache a class per parameter mask and retarget instances by
assigning `__class__`. The mask is an **opaque tuple**; nothing in the
machinery knows it means "collapse".

Generalising it to "elide any block whose parameters make it inert"
needs: a new `Statement` subclass (`Collapse` is the template, including
its parameters-only validator at `:971`), a partition at `:1172`, an
application step beside `:1189`, a widened mask tuple, and a smarter
freeze check — because a guard that removes *terms* but not *nodes*
does not change the unknown count, so it may legally change after build.

This fixes a trap that cost real time: **a block whose parameter is zero
is still compiled and still evaluated**, because the element compiles
from symbolic parameters, so `if param != 0` is true at build time
whatever the value — and `0 × inf = NaN`. The model currently guards
twelve blocks by hand with `if isinstance(p, sympy.Expr) or p != 0.0:`,
which only fires for the *default* build; an instance whose card sets
the parameter to zero still evaluates the whole block.

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
    3        S1  pure form for chained path          not started
    3        S2  generalised variants                not started
    3        S3  automatic domain clamping           research
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
    14       backend: the spike that decides it      NOT STARTED, and
                                                     needs a decision
                                                     before it starts
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
    15       the unlimited-node failure              OPEN, architectural
                                                     on every path
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
