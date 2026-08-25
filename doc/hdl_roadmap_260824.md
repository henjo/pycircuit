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
| 1 | **Self-heating thermal node** + Foster/Cauer ladder | electrothermal | 60 | ~10 lines per device, no new operator: the thermal node's "voltage" is ΔT, the device contributes `v·i` into it and reads `TEMP + V(th)` back. Retroactively upgrades every device, and closes a named gap in `PspMosLongChannel` |
| 2 | Full SPICE diode (N, RS, TT, CJ0/VJ/M, BV, shot+flicker noise) | semiconductor | 80 | today's `Diode` is a bare exponential; `juncap.py` already holds the kernels |
| 3 | SPICE Gummel-Poon BJT | semiconductor | 300 | the most conspicuous absence in the repo |
| 4 | Opamp macromodel (Boyle-class) | behavioural | 150 | deletes `macromodels.OpAmp`'s runtime-`import jax` lambdas; the most-used block in any testbench |
| 5 | Quartz crystal / BVD resonator | passive | 40 | nothing in the repo can simulate an oscillator's frequency-setting element |
| 6 | VCO with real phase noise | behavioural | 40 | noise injected into the *frequency*, so jitter accumulates as a random walk — a capability the repo's own survey found in no other simulator |
| 7 | Comparator with hysteresis | behavioural | 60 | the **first production user of `Cross`**, which Phase C1 built and nothing consumes |
| 8 | Skin-effect R + ferrite bead | passive | 100 | the **first production users of `laplace_zp`**, likewise unconsumed |
| 9 | EKV 2.6 MOSFET | semiconductor | 150 | one tenth of PSP's size, all-region, and an independent cross-check on the surface-potential kernel |
| 10 | MOS Level 1 + Level 3 | semiconductor | 330 | legacy netlists and teaching |
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

### 10.4 Sequence, and the gate on the expensive one

1. **10.3(a)**, FET limiting on the ordinary path. Independent of
   everything else here, and the only item that helps PSP today.
2. **10.2**, the chained detector. Contained, well-tested by
   construction, and it removes one of the two capabilities `var()`
   costs (section 9).
3. **10.3(b)**, vector PCNR. Its own project.

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
    3        S5  parameter namespace                 DONE (_analog_function)
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
    5        thermal node                            DONE
    5        full SPICE diode                        DONE
    5        the other thirteen models               not started
    6        Phase 1 make the record true            DONE
    6        Phase 2 designer's surface              DONE
    6        Phase 3 kernel additions                DONE
    8        printer/scalar codegen                  measured, not started
    10.2     PCNR admits chained models              DONE
    10.3(a)  FET limiting, ordinary Newton path      DONE
    10.3(b)  FET limiting inside PCNR (vector)       not started
    5        Gummel-Poon BJT, EKV core, crystal,
             ferrite, skin-effect R, comparator,
             memristor                               DONE
    5        the remaining eight models              not started

Suite: 1874 tests at `4d359e5` before this work, **2184 passed / 7
skipped / 3 xfailed / 0 failed** as of `cb477a8`. `sphinx -W` clean.

**Corrections this campaign made to its own record.** Four claims were
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

**The pattern is worth naming: a reason decays faster than the
conclusion it supports, and nothing fails when it does.** Every one was
caught by taking a measurement rather than by reading, and most were
one `_src` dump or one `hasattr` away the whole time. Three of the five
were written by whoever was correcting the previous one -- including
two of mine -- so this is not a failure of care, it is what happens to
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
> And a side effect worth knowing: a stacked pair at (5 V, 2.5 V) that
> used to raise `SingularMatrix` now converges, because the better
> write-back keeps the iterate out of the cutoff region where the lower
> device's conductance underflows to exactly zero. The hazard is NOT
> gone — (40 V, 0.2 V) still reaches it — only harder to reach. §12.3
> stands.

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

### 12.2 `check_jacobians` needs a "not resolvable" verdict — third sighting

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

### 12.3 No way to say "this device needs a GMIN anchor"

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

- **`limit_pnj`'s parameters cannot be written with `var()`.** They are
  lambdified over the parameter namespace, so an `_IntermediateSymbol`
  is not in scope. The BJT spells its temperature-scaled `IS` twice as
  a result. Resolving `var()` symbols against the chain before
  lambdifying removes the duplicate; refusing loudly at compile time
  would at least remove the guessing.
- **`limit_fet`'s parameter-only `vto` bites, and the size is known:**
  for the EKV card at 2 V of body bias the true turn-on is 1.06 V
  against `vto = 0.50`, so every clamp sits **565 mV** low. Looser, not
  wrong — the step is still bounded and the no-op band intact.
- **`0 × d(√0)` is reachable from a legitimate CARD, not just a wild
  bias.** EKV's pinch-off root with `gamma = 0` — a perfectly good "no
  body effect" card, and the one every textbook-asymptote test uses —
  sits exactly at zero with an infinite derivative. This is the
  `differentiable-numerics` rule in a form that skill does not spell
  out: the dangerous value was a *parameter*, not a bias.
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
