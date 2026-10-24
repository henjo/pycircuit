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

> ⚠ **THAT ANNOTATION WAS HALF WRONG, corrected 2026-08-27.** It read:
> *"the docstring complaint below is stale: `var()` now opens 'This is
> what makes a real compact model compilable'"*. The **opening** line had
> indeed been fixed. The **closing** line had not — the docstring still
> ended with the exact sentence S4 complains about, *"Naming is optional
> and only affects generated-code readability"*, contradicting its own
> first paragraph. Reading the top of a docstring is not reading the
> docstring. Fixed now, with the wrong sentence kept as a quotation.
>
> 📌 **And the finiteness half of the argument below does NOT reproduce.**
> Measured 2026-08-27 by replacing `_v` with the identity for selected
> holds and rebuilding PSP against a real IHP card:
>
> * **compilability — confirmed, dramatically.** All holds dropped: does
>   not finish in **10 minutes**, against **62 s** held.
> * **finiteness — not reproduced.** Dropping only the `sg_*` holds, the
>   ones `_sigma_body`'s docstring says lost the Jacobian at `Vd = 1e26`
>   *"purely from evaluation ORDER"*, leaves the model finite to
>   **1e36** — unchanged from held. Something since (most likely
>   `hdl._rdiv`'s regrouping, which moved the bound from 1e27/1e24 to
>   1e33 on both channel types) has made those particular holds no longer
>   load-bearing for finiteness. The docstring is left as written and
>   annotated, not deleted: it records what was true when the holds went
>   in.
>
> Two further measurements that shape the work: **52% of the holds are
> ≤3 ops** (median 3, max 54) — so "hold everything" mostly holds naming;
> and **434 of 1184 kernel-regulariser calls take a bare non-atomic
> argument**, which is where the risk actually concentrates. See §36.

`var()`'s docstring used to say naming is *"optional and only affects
generated-code readability"*. The model's experience contradicts that
flatly: bare sub-expressions get substituted into one another before
printing, so the compiled form evaluates products the written form never
creates, and the n-channel device lost its Jacobian at `Vd = 1e26`
*purely from evaluation order*. There are **342** hand-tagged `_v()`
calls in `psp_kernel.py`.

If the compiler held every non-atomic sub-expression automatically
(structural hashing at build time), `var()` would go back to being the
naming convenience its docstring once claimed. ~~Until then the
docstring is wrong and should say so.~~ — done; the docstring now
states the opposite, and correctly.

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
| `limit_delta(probe, vmax)` | 35 hand-written lines of FET Newton limiting per model — **and see the note below: it now has a named consumer** |
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
| 10 | MOS Level 1 + Level 3 | semiconductor | 330 | legacy netlists and teaching. **Level 1 DONE 2026-08-25** as `MosLevel1Hdl`/`MosLevel1PmosHdl`, and it is `limit_together`'s first production user; **Level 3 DONE too** as `MosLevel3Hdl`/`MosLevel3PmosHdl` (`elements_hdl.py`) -- this cell read "Level 3 not started" until the 2026-08-31 audit, long after both classes shipped and were listed PCNR-eligible in the reviewer's guide |
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
    3        kernel: limit_delta                     DONE 2026-08-26 (sec.
                                                     34): the family is
                                                     complete.  EKV adopted it
                                                     on its bulk at 10 V (sec.
                                                     35), which retired the
                                                     last PCNR fallback
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
    32       re-validate PSP folded, vs the vendor   DONE 2026-08-26: worst
                                                     point 1.305e-06 BOTH
                                                     folded and direct over
                                                     all 19 sweeps; the two
                                                     agree to 6.7e-14.  The
                                                     figure is now PINNED by
                                                     a test -- it never was.
                                                     Was raised as: OPEN --
                                                     the item sec. 31
                                                     created.  Folding is
                                                     not bit-identical and
                                                     the 1.3e-6 vendor
                                                     figure has NO test
                                                     guarding it (the
                                                     suite asserts 1e-4)
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
    3        S4  automatic intermediate holding      DONE, cut down (sec.
                                                     36): the regularisers
                                                     hold their own args,
                                                     safe_div included.
                                                     PSP G 1.22x faster.
                                                     safe_div was excluded
                                                     for a day on a
                                                     misreading -- sec. 37.
                                                     Full S4 refused
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
    15       vector PCNR Stage 3 (JAX)               DONE 2026-09-01.  All
                                                     five gates; DC to
                                                     4.4e-16 and transient
                                                     to ~1e-5 vs the CPU,
                                                     premium +4-13% wall
                                                     and 0% steps.  Sec. 49
                                                     carries the work; 49.7
                                                     gave the CHAINED path
                                                     a jax backend on the
                                                     way, which is a
                                                     capability and not a
                                                     speed-up (~18x slower
                                                     than the CPU, 49.8)
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
    38       published perf figures are unguarded    DONE 2026-08-27 (sec. 39):
             by any test                             sec. 38 OVERSTATED -- fold_
                                                     card and _autohold already
                                                     had DIRECTION guards.  A
                                                     magnitude bound on those was
                                                     REFUSED (4.4% margin = noise).
                                                     The overhead figure, guarded by
                                                     neither, now has test_perf_
                                                     guards.py: interleaved ratio,
                                                     confirm-before-failing, bite-
                                                     checked with a calibrated slow-
                                                     down

    40-44    PCNR on the cards a PDK ships           DONE 2026-08-27: 16 of 37
                                                     models eligible, against ~4
                                                     on real cards before.  Affine
                                                     remainder for constant series
                                                     R; probe declarations for GP
                                                     rbb, SelfHeating dT and the
                                                     diode series branch; and the
                                                     V-branch gate relaxed (it
                                                     refused on PRESENCE).  All
                                                     gated on OPERATING-POINT
                                                     agreement, 0.0e+00 to 9.0e-17.
                                                     Only PhotodiodeHdl still
                                                     refused with an exponential:
                                                     iph reads optical INPUTS

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

> ✅ **Measured before building, §31: it stacks — 1.44x on C**
> (82.2 → 57.2 µs), and the C build halves (27.1 → 14.8 s) because
> there is less code to compile. That is what justified building it.

### Status

> ✅ **BUILT, same day — see §31.** `hdl.fold_card`, 1.54x numpy /
> 1.48x C on PSP, and it compiles faster than the symbolic form. The
> accuracy question below is answered and carried into the API's
> docstring, `hdl.rst` and the test module: **not bit-identical, so a
> model validated against a reference must be re-validated folded.**
>
> §31 also records the design error this section did not anticipate:
> folding only the parameters a card *mentions* leaves the rest
> symbolic and measures **0.96x**. The API has to name what stays
> VARIABLE.

Written before the build, and all three held:

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

## 32. DONE: re-validate PSP against the vendor with its card folded

> ✅ **Done 2026-08-26, same day it was raised.** Folding costs
> **nothing**: worst point against IHP's compiled PSP103 is **1.305e-06
> folded and 1.305e-06 direct**, over all 19 recorded sweeps, with the
> two builds agreeing to **6.7e-14** point by point. And the number is
> now **pinned by a test**, which was the real point. Results at the
> bottom of this section.

*(Written before the work:)* **This is the one item `fold_card` created,
and it is the one that matters most, because that model's whole value is
a number that folding is not guaranteed to preserve.**

### Why it exists

`hdl.fold_card` (§31) is **not bit-identical** and cannot be — folding
reassociates arithmetic. §30 measured the drain current banded by
magnitude and found **8.5e-16** relative in strong inversion and
**2.7e-15** in the 1e-9…1e-6 A window, i.e. machine precision where it
counts. That is strong evidence and it is **not the same thing** as
having re-run the validation.

`compact.PspMosLongChannel`'s claim is: *"agrees with IHP's own compiled
PSP103 to 1.3e-6 at the worst point of twelve sweeps, two geometries and
both channel types — the reference's own printed precision."* Nobody has
checked that this survives folding.

### Why the existing evidence does not settle it

Three reasons, and they are the ones this campaign keeps rediscovering:

1. **The measured window was not the validated one.** §30's banding used
   a gate sweep at one drain bias. The vendor comparison is *twelve*
   sweeps across two geometries, both channel types and a body-bias
   ladder — including the `alp2` and subthreshold cases that took the
   longest to close, and the `1e-9…1e-6 A` window has a **junction
   leakage floor** beneath it that this core does not model.
2. **The tolerance is 43x tighter than the test that guards it.**
   `test_psp_gap.py` asserts a part in a thousand
   (`0.9999 < ratio < 1.0001`); the 1.3e-6 is the *measured* worst point
   and lives in prose. So **the suite would stay green through a
   regression of three orders of magnitude** — the exact shape of §25's
   drifted `hdl.rst` figure.
3. **Folding at the vendor's card is not folding at defaults.** §31
   showed defaults are not a physical card (they raise). The real card
   supplies 66 of 153 parameters; the other 87 fold at *their defaults*,
   and no measurement has looked at what that does to the terms those
   defaults switch off.

### What the work is

Re-run `test_psp_gap.py`'s twelve-sweep comparison against
`tests/data/psp103_ihp_iv.json` with the device built through
`hdl.fold_card(..., instance=('w', 'l'), **card)` instead of directly,
and report the worst point per sweep beside the unfolded number.

Then one of:

* **it holds at 1.3e-6** — record it, and add the folded build as a
  second parametrisation of the existing tests so it cannot drift; or
* **it degrades** — report by how much and in which sweep, and either
  document folding as unsuitable for a vendor-validated model, or find
  the term whose reassociation costs the precision.

⚠ **Pin the number this time.** Whatever it comes to, assert it — the
1.3e-6 currently has no test, which is why this item can even be asked.

### Size

Half a day. The harness, the reference data and the card loader all
exist; `benchmarks/card_constant_folding.py` already builds the folded
device from the same card.

### RESULT (2026-08-26)

`_compare` gained a `fold` switch rather than a second code path, so
both builds run the identical comparison. Worst |ours/ref − 1| per
sweep, above the 1e-6 floor:

| | worst point |
|---|---|
| direct | **1.305e-06** (`nmos_long_idvg`) |
| folded | **1.305e-06** (`nmos_long_idvg`) |
| folded vs direct, point by point | **6.7e-14** worst of 19 sweeps |

**Folding is free here.** The concern in reason 3 above -- that 87 of
153 parameters fold at their *defaults* -- turned out not to bite: the
worst sweep is the same one, at the same value.

**Two tests, both verified to bite:**

* `test_the_worst_point_over_every_sweep_is_pinned[direct|folded]`
  bounds the worst point at 1.5e-6 over **every** sweep. A 1e-5
  perturbation of the model -- eight times the bound, and *under* the
  1e-4 the surrounding tests assert -- fails it.
* `test_folding_the_card_does_not_move_the_vendor_agreement` compares
  the two curves **point by point**, because two builds could each sit
  1.3e-6 from the reference and 2.6e-6 from each other. A 1e-11
  perturbation of the folded path alone fails it.

**And reason 2 was worse than written.** The 1.3e-6 had no test at all,
and the surrounding tests assert a part in a thousand -- so the suite
would have stayed green through a **three-order-of-magnitude**
regression, for the whole life of the model. That is now closed.

⚠ **A stale count found on the way:** the model's docstring, this
document and eight other places say *"twelve sweeps"*. The reference
holds **19** -- it is regenerated additively, so twelve was true when
written. Corrected in `compact.py` (the model's own public claim, and
the one a user reads) with the old value kept visible; the dated
narrative in `hdl.md` is left as the record of what was true then.


## 33. `limit_delta` has acquired a consumer (2026-08-26), and S4 has acquired a blocker

Prioritising the two remaining kernel items against each other, after
`fold_card` landed. They are not close.

### `limit_delta` — half a day, and something is waiting for it

The limiter family shipped complete except this one: `limit_pnj`,
`limit_fet`, `limit_vds`, `limit_identity` and `limit_together` all
exist. §3 justified `limit_delta` generically ("35 hand-written lines
per model"). It now has a **specific** consumer, which §3 could not have
known:

`EkvNmosHdl` declares `limit_identity(bsb.V)` on its bulk, and the
source says what that is: *"NOT a limiter -- it never moves the branch
-- but the declaration that admits this model to vector PCNR."* The
device owns all three branch unknowns, so the tail-node clash cannot
form; but the bulk probe has **no law**.

That is exactly why the EKV differential pair is the one circuit still
riding the PCNR fallback (§27, §31): the wild-start fix was to limit the
SEED through each device's own law, and an identity law clamps nothing.
`limit_delta(bsb.V, vmax)` would give that probe a real law -- bound the
change per iteration -- with no device physics required.

**Ready-made acceptance case**: the EKV pair from a uniform 20 V start,
which currently raises `LinAlgError` and falls back
(`test_dc_pcnr.py::test_dc_pcnr_true_falls_back_...` runs on it *because*
it still fails). Convergence there, order-independently, would close the
last PCNR gap and retire a fallback.

Risk is low: additive, with four templates to follow.

### S4 — weeks, and blocked behind §32

S4 would change what the compiler emits for **every** model. §31 just
demonstrated what that means: `fold_card` changes emission and is
machine-precision but **not bit-identical**. S4 deserves the same
scrutiny.

Which is the blocker: **PSP's 1.3e-6 vendor agreement has no test**
(§32). Without one, there is no way to tell whether S4 preserved it —
the suite asserts a part in a thousand and would stay green through a
three-order-of-magnitude regression. **Doing S4 before §32 means
changing every model's arithmetic with the flagship model's only real
acceptance number unguarded.**

Nothing is waiting on S4, either. 342 hand-tagged `_v()` calls is
friction for whoever writes the next big model; it is not a defect, and
no measurement asks for it.

### Verdict

**`limit_delta` first**, by a wide margin — days against weeks, a named
consumer against none, a ready acceptance test against none, and low
blast radius against total. **S4 is not merely lower priority; it is
sequenced behind §32**, and §32 is half a day.

## 34. `limit_delta` BUILT (2026-08-26): the limiter family is complete

§33 ranked it first: half a day, and something waiting for it. Built.

```python
vsb = limit_delta(bsb.V, 1.0)       # do not move this branch >1 V per step
```

`pnjlim` compresses a junction logarithmically, `fetlim` works about a
threshold, `limvds` has hard-coded breakpoints — each needs a model that
*has* such a quantity. `deltalim` needs nothing, which makes it the law
for a branch whose model offers no better one. The kind table, the
dispatch and the group validator now derive from one dictionary, so a
sixth kind touches three lines rather than eight sites.

**Not invented here.** `compact.PspMosLongChannel.limit` already bounds
d, g and b alike by a symmetric `|delta| <= 1 V` about the source, by
hand — and `_limvds`'s docstring already called that out as *"a third,
blunter law"*. This is that law, named and shared, which is exactly what
§3 asked for.

### The acceptance case: it fixes the EKV diff pair completely

`EkvNmosHdl` declares `limit_identity` on its bulk — a probe with **no
law**, present only to admit the model to vector PCNR. §27's wild-start
fix limits the *seed* through each device's own law, and an identity law
clamps nothing, which is why this is the one circuit still falling back.

Diff pair, both instance orders, seven `vin` values, `.` = converged:

```
             start 0 V                start +20 V
identity     .. .. .. .. .. .. ..     XX XX XX XX XX XX XX
delta  1     .. .. .. .. .. .. ..     .X .. .. .. .. .. X.
delta  5     .. .. .. .. .. .. ..     .. .. .. .. .. .. ..
delta 10     .. .. .. .. .. .. ..     .. .. .. .. .. .. ..
```

**14 of 14 failures become 14 of 14 convergences at `vmax >= 5`**, order-
independently, with every zero-start result unchanged and every tail
equal to `DC()`'s.

⚠ **The bound is NOT monotone**, and 1 V — the value `compact.py` uses —
is the wrong choice here: it leaves the two `|vin| = 1` cases
**order-dependent**, which is the one property PCNR exists to provide.
Swept at `vin = +1`: 0.1, 0.25, 0.5 fail; 1.0 half-works; 2.0 fails; 5.0
and 10.0 work. Same chaotic signature as §27's `gmin` sweep, and the same
lesson: **an iteration count from a wild start is not a stable
acceptance number**, so pick the bound from a grid, never from one point.

### EKV is deliberately NOT changed

The measurement says it should be — and that is a separate decision:

* it changes numerics for anyone using `EkvNmosHdl`;
* `test_limit_identity.py` is built around that declaration. Its claims
  2 and 3 assert EKV qualifies *via identity* and measure this very diff
  pair with it, so the change invalidates a test module rather than
  breaking a line of it;
* the right `vmax` needs its own justification — 1 V is what the only
  in-tree precedent uses and 1 V is the value that does not work here.

So: the primitive ships, the measurement is recorded, and the model
change is the user's call with the numbers in hand.

### Tests

`test_limit_delta.py`, 18 tests, both mechanisms verified to bite
(neutralising the law fails 6; removing the positivity check fails 3).
Covers the law itself (bites symmetrically; **returns `vnew` itself**
when it does not, which is how both write-back paths detect "did
limiting fire?"), the declaration through a real model, `explain()`,
vector-PCNR qualification, a parameter-valued bound, refusal of a
non-positive bound, and grouping — which a delta probe may join and an
identity probe may not.

## 35. EKV switched to `limit_delta` (2026-08-26): the PCNR fallback has no vehicle left

§34 built the primitive, measured that it fixes the EKV diff pair, and
**deliberately did not change the model** — that being a separate call.
Made now.

```python
vsb = _var(limit_delta(bsb.V, 5.0), 'vsb')   # was limit_identity(bsb.V)
```

### 10 V is chosen from TWO grids — and 5 V, chosen from one, was wrong

The bound is not monotone, so it was swept on the full grid — seven
`vin` values × both instance orders × a uniform 20 V start:

```
identity     XX XX XX XX XX XX XX
delta  2     XX .. .. .. .. .. XX
delta  3     .. .. .. .. .. .. ..
delta  5     .. .. .. .. .. .. ..      <- chosen
delta 10     .. .. .. .. .. .. ..
delta 20     .. .. .. .. .. .. ..
delta 50     XX XX XX XX XX XX XX      <- too loose to limit at all
```

**The working plateau is 3–20 V.** Below 3 it is partial *and
order-dependent at |vin| = 1*, which is the one property PCNR exists to
provide. **1 V — the value `compact.PspMosLongChannel.limit` uses by
hand, and the obvious guess — does not work here.**

5 V was picked from that grid, and the suite caught it. The 40 V cascode
benchmark (`test_fet_limiting_cuts_the_iteration_count_on_a_hard_solve`)
dropped from a ratio of 1.9 to 1.79 and failed its bound. Swept against
*that* circuit, Jacobian evaluations:

| bulk law | identity | δ1 | δ2 | δ3 | **δ5** | δ10 | δ20 |
|---|---|---|---|---|---|---|---|
| cascode | 13 | 20 | 16 | 13 | **14** | 13 | 12 |

**5 V is the one value inside the plateau that costs that solve an
iteration.** 10 V is inside both grids: fourteen of fourteen on the diff
pair, and it ties the identity probe on the cascode. Shipped at 10.

⚠ Two grids, two different answers about which end of the range is good.
A bound picked from one circuit is a bound picked from a point, however
many values of it were swept.

Zero-start results and every converged tail are unchanged.

### ⚠ Nothing makes PCNR fail any more

Probed after the change: **72 configurations** — EKV diff pairs, the BJT
mirror, Level-1 cascodes, `_fet` cascodes at three geometries — starts
from **−40 V to +100 V**, both instance orders. `pcnr.solve_dc`
**converged on every one**.

That is a good result and an awkward one, because two tests used "PCNR
fails here" as their **precondition**:

* **the `DC(pcnr=True)` fallback test** has now run on three circuits.
  It began on the BJT mirror (fixed by §27's seed), moved to the EKV
  pair (fixed here), and has run out. It now **forces** the failure by
  raising from `solve_dc` — a weaker test of *when* the fallback fires
  and a better test of *what it does*, and one that cannot rot again.
  A companion test asserts both former vehicles now go through PCNR.
* **the `gmin` damper test** (§31) asserted "raises undamped, converges
  with it". The first half is false now.

**The fallback stays.** "Nothing I can find fails" is not "nothing
fails", and it costs nothing when it never fires.

### The gmin damper has lost its case

§31 built it for exactly this circuit and shipped it **off** because it
traded the pair for cascode grid point (2, 2, 2). The pair no longer
needs it, so it now **rescues no circuit anyone can name** while still
costing the grid point when switched on. Kept — it is the only tool for
a genuinely rank-deficient start — but `solve_dc`'s docstring now says
plainly that it should not be reached for without a failing case in
hand.

### `limit_identity` now has no user in the library

It remains supported and tested (section 1 of `test_limit_identity.py`,
on synthetic devices), and it is still the right declaration for a
branch whose model genuinely has no opinion. It is simply no longer what
any shipped model uses — because §27 made the *seed* the thing that
matters, and a probe with no law has nothing to apply to it.

### Five tests re-aimed, none patched

Each kept its property and lost its incidental:

| test | was | now |
|---|---|---|
| `..._qualifies_for_vector_pcnr_...` | qualifies *through the identity probe* | qualifies through **a bulk probe**; the kind is incidental |
| `..._refusal_still_fires_...` | without *the identity probe* | without **a bulk probe of any kind** |
| `..._gmin_pair_view_stays_pnj_only` | kinds `('id', 'fet', 'vds')` | `('delta', ...)`; the load-bearing `pairs == []` untouched |
| `..._diff_pair_from_twenty_volts` | recorded `[F, F]`, deliberately not asserted | **asserts convergence** — otherwise a regression to `[F, F]` passes |
| `..._gmin_damper_is_off_by_default_...` | "raises undamped, converges damped" | off by default, and **cannot move the answer** |
| `..._declares_fetlim_on_the_gate_...` | "the bulk is never written" | the bulk **is** a probe terminal now, and is written 147/300 |

**And one of those was a behavioural change, not a rename.** *"The bulk
is not a probe terminal and is never written"* was true only because
`limit_identity` never bit. A probe with a law has to write back, so the
bulk now moves — and three probes contend for the shared source node
where two did before. Measured: **3 of 600 branch-bound checks (0.5%)
end with a branch wilder than it went in**, which is the same per-probe
write-back interference `limit_vds`'s docstring already records at 3.3%
over 813 two-probe steps. `limit_together` is what satisfies all the
constraints at once, and EKV does not use it. Recorded in the test with
the number, not smoothed away by loosening the bound.

The fourth is worth its own line: its docstring said *"a future damped
`predict` that converges here should make this docstring stale, not the
test red"* — and it did exactly that. **It was right about the outcome
and wrong about the cause.** What fixed it was not a damped `predict`
but a seed with a law to apply.

## 36. S4, cut down and measured (2026-08-27): the regularisers hold their own arguments

S4 — *"hold every non-atomic sub-expression automatically"* — has sat as
**research, weeks** since §3. Investigated properly, the full form is
still not worth it, but a narrow slice is, and the investigation
corrected the record S4 is argued from.

### What the investigation found, before any code was written

| question | measured |
|---|---|
| is `var()` load-bearing for **compilability**? | **yes, overwhelmingly** — with `_v` the identity, PSP does not finish in **10 min**, against 62 s held |
| for **finiteness**, as recorded? | **the recorded case does not reproduce** — dropping the `sg_*` holds, whose docstring blames evaluation order for losing the Jacobian at `Vd = 1e26`, leaves the model finite to **1e36** |
| how big are the holds? | median **3 ops**, max 54 — **52% are ≤3 ops**, i.e. naming |
| where is the risk? | **434 of 1184** regulariser calls take a bare non-atomic argument, `safe_div` one of **44 ops** |

So "hold everything" mostly holds naming, and the risk is concentrated
at the regularisers. That is the slice.

### Phase 1 — `benchmarks/hold_value.py`

Drops holds by prefix and reports compile time, emitted size and
finiteness reach. **342 assumed-necessary holds become a measured list:**

```
dropped              compile   G lines    finite to
nothing (baseline)    132.0s      2860        1e+36
sg_*                  133.0s      2800        1e+36
ids_*                  >300s        --   DID NOT BUILD
avl_*                 123.8s      2846        1e+36
EVERYTHING             >300s        --   DID NOT BUILD
```

**The `ids_*` holds are load-bearing for compilability; `sg_*` and
`avl_*` are not, for either property.** `">300s" is a result, not a
hang` — and the harness had to learn that: the first survey sat for
thirty minutes on one prefix with nothing printed, so each build now
runs in a bounded subprocess.

**A compile-time budget warning** ships with it (`hdl.COMPILE_WARN_SECONDS
= 180`). Compile time is the *only* usable signal for a forgotten hold:
the emitted line count goes **down** when a hold is dropped, the
expressions are small either way, and the finiteness reach did not move.

### Phase 2 — `_autohold` in the regularisers

`safe_div`, `expl`, `hypsmooth`, `safe_ln`, `safe_pow` and `softplus`
hold their own arguments, so the author does not have to. **PSP `G` is
1.22x faster** — 25 paired samples on a loaded machine, per-pair ratio
median **1.216x** (p10 1.168, p90 1.284) — with the build unchanged and
the finiteness reach still 1e36. Holding a repeated sub-expression
computes it once instead of substituting it twice.

(Without `safe_div`, which was the shipped state for a day, it is
1.12x: 17.680 ms unheld, 15.784 without, 14.121 with.)

### ⚠⚠ The defect the four planned gates did not catch

`var()` is what makes a model **chained**, and a chained model has **no
`eval_i_pure`** — so it loses `solve_batched` (20.7x at 512 lanes) and
the JAX path. §22 measured the two fast paths as **disjoint**;
auto-holding unconditionally moves any eager model that divides or
exponentiates out of the batchable set.

**Caught by `KernelAllEager` in `test_hdl_kernel.py`, which lost
`eval_i_pure` outright.** Not by the vendor gate, not by the finiteness
gate, not by the 37-class digest — and not by a flip-check I ran across
all 37 library classes, which reported *no flips*, because the 15
batchable classes are passives that never call a regulariser. **The
library could not have shown me this; only a test model did.**

Fixed by holding **only once the body has declared an intermediate
itself** (`_VAR_STACK[-1]` non-empty): it is chained anyway by then, so
holding more is free.

### ⚠ `safe_div` was excluded for a day, on a misreading — see §37

The site where the value is largest — the biggest bare arguments, 100 of
the 434. It was excluded because holding there made
`TestPspBitIdentity` report `{'equal': 989, 'zero-sign': 16, 'value':
0}` on `G`, read as *16 signed zeros*; a signed zero is observable
(`1/-0.0` is `-inf`), so refusing to erode the C backend's bitwise
guarantee was right **on that reading**.

**The reading was wrong.** Traced in §37: they are 240 **NaN sign bits**
at 16 bias points, all where `G` is already NaN in both paths. IEEE-754
does not define that bit. `safe_div` now holds, and the shipped figure
is the one in the table below.

### Accuracy

Not bit-identical, and cannot be: holding stops sympy flattening across
the boundary, so the arithmetic reassociates. Measured on cards that
**reach** the changed arms (§28's lesson — `nsub` defaults to 0, so a
default card never gets there):

| class | max relative | median |
|---|---|---|
| MOS level 3 (n and p) | **2.535e-16** | 1.669e-16 |
| `DiodeSpiceThermalHdl` | **1.213e-15** | 1.568e-16 |
| `DiodeSpiceHdl` | **unchanged** — only its `explain()` text moved | |

One ulp. Five recorded digests re-recorded with the measurement beside
them; the vendor figure (1.305e-06) and the finiteness ladder (1e33) are
untouched.

**`test_hdl_select.py`'s recorded digests caught what my own 37-class
digest missed**, because they carry cards that reach the arms. It is the
better instrument.

### Status

Phases 0–2 done, about two days as estimated. **Phase 3 is not
authorised and Phase 4 is refused** — the full form would intercept
sympy construction for every model, and 52% of what it would hold is
≤3 ops.

## 37. The "16 signed zeros" were 240 NaN sign bits (2026-08-27)

§36 excluded `safe_div` from `_autohold` — the single most valuable site,
100 of the 434 bare arguments and the largest of them — because enabling
it made `TestPspBitIdentity` report:

```
('G', {'equal': 989, 'zero-sign': 16, 'value': 0})
```

read as *"16 signed zeros in PSP's `G`"*. The reasoning was sound on that
reading: a signed zero is an **observable value** (`1/-0.0` is `-inf`),
the C backend's shipped guarantee is bitwise identity, and the tree had
already drawn a deliberate line — `C` may carry signed zeros (*"the
integer-lattice exception — PSP's charge Jacobian has integer-valued zero
cells"*), `i`, `G` and `q` may not. Crossing a line someone drew on
purpose, to buy speed, was refused.

### What they actually are

Traced entry by entry:

```
240 byte-differing entries, 16 bias points
by KIND : [('nan-payload', 240)]
  numpy fff8000000000000   C 7ff8000000000000
  bias (0.0, 1.2, 1e+30, -1e+30) -> G: 20 of 36 cells NaN
```

**Not one zero, signed or otherwise. Not one differing value.** They are
NaN *sign bits*, at extreme biases where `G` is already NaN in both
paths.

**IEEE-754 does not define the sign of a NaN produced by arithmetic**,
nothing preserves it, and no consumer can observe it — every use of a
NaN yields a NaN. The exclusion was protecting a bit that is not a
computation, and it cost 15.784 → 14.121 ms.

### The instrument was ambiguous, and that is the transferable part

`_compare` returned `'zero-sign'` whenever *values compared equal and
bytes did not* — which is true of a signed zero **and** of two NaNs with
different sign bits. One label, two findings of very different weight,
and the label named only the serious one. So the tally read as the worse
of the two possibilities every time.

Split now: `'nan-bits'` and `'zero-sign'` are separate tallies. The
strict assertion still requires `value == 0` and `zero-sign == 0` — a
signed zero remains banned, because it is observable — while `nan-bits`
is **recorded rather than ignored**: `G <= 32`, `i == 0`, `q == 0`, so a
change that starts producing NaNs somewhere new still fails here.

**A diagnostic that collapses two findings into the name of the graver
one will make you refuse the right thing for the wrong reason.** This
cost a day and the largest available win. The tell was there to be
noticed: the tally said `zero-sign` on a *Jacobian* at ±1e30, where zeros
are not what one should expect to find.

## 38. Writing the branch up for a reviewer (2026-08-27): five wrong claims in my own draft, in thirty minutes

The roadmap has no open items left except the one the user parked. The
remaining bottleneck is not engineering: **551 commits, 442 files,
+160,343 / −3,904, none of it reviewed.** So the work became a document
that lets the branch be reviewed by subsystem —
`doc/branch_review_260827.md`.

Summarising is not the same activity as building, and it has a failure
mode of its own: **a summary is written from memory, and memory is
exactly where the rot is.** Verifying each claim before writing it down
caught five errors in half an hour, every one of them mine, and every
one of them something I would have asserted to the repo owner with
confidence.

| I was about to write | Actually |
|---|---|
| "36 models" | **37.** `RSkinHdl` is built by the `skin_effect_resistor()` factory, so `grep '^class .*(Behavioural)'` cannot see it. Reflection finds it; text search does not. |
| "toolkit memoisation — refused on measurement" | **Built**, and it is §24's 1.16x. I had it filed under the wrong outcome. |
| "none of the root scratch files is referenced" | `example10.py` is a **documented example** — two `.rst` files include it. Deleting it would have broken the docs build. |
| "`MesfetCurticeHdl`/`HemtAngelovHdl` are the under-evidenced models" | They are checked against their **published equations**, by finite-difference Jacobians, with limiters bite-checked. |
| "library equivalence across all models, via recorded digests" | Digests pin **five** models. The library-wide guarantee is a different mechanism — the reflection-discovered C-backend parity sweep. |

### The transferable part: a coverage proxy misleads in BOTH directions

I reached twice for a cheap proxy for "how well is this model tested",
and got two different wrong answers.

- **Test files per model** said `MesfetCurticeHdl` and `HemtAngelovHdl`
  were the weak spot. They have one test file, which happens to be a
  *dedicated* one with published-equation checks in it.
- **Test functions per model** corrected that, and is better — but it
  still cannot see a model covered only through a parametrized sweep
  over a reflected class list, which is precisely how every model gets
  its floor here.

A count is evidence about the count. Reading four test names settled in
seconds what two proxies had got wrong. **Where a number is a stand-in
for a judgement, print the number and then go and look.**

### The finding worth acting on

**No performance figure on this branch is asserted by any test.** The
cumulative 1.40x, the DSL overhead at parity, `fold_card`'s 1.54x,
`_autohold`'s 1.22x — all measured once, on one machine, recorded in
prose.

> ⚠ **CORRECTED 2026-08-27, same day (§39).** That sentence is
> overstated, and the correction is the useful part. Two of the four
> optimisations already have a **direction** guard: a test asserts
> `fold_card` emits fewer primitives than the symbolic class, and
> another asserts the regularisers really hold. Neither can silently
> stop working. What none of the four has is a **magnitude** guard, and
> the end-to-end DSL overhead had no guard of *either* kind. The
> distinction matters because it changes what is worth building: three
> of the figures needed nothing, and the fourth needed a real test. `hdl.rst`'s overhead figure has already drifted once under
exactly these conditions (1.14x published, 1.25x when re-measured, with
a green benchmark in the tree the whole time, because it asserted the
waveforms agreed and never the ratio).

This is a standing open item, not a closed one: a bound-assertion around
each published figure would convert silent decay into a failing run.

### The nine root scratch files, deleted (2026-08-27)

A bare `pytest` at the repo root aborted during collection —
`Interrupted: 5 errors` — because of five `test_*.py` files at the root
with **no test functions** and broken imports, swept in by a
`git add -A` on 2026-07-31 in a commit about `Toolkit.__getattr__`.
Nine such files in all. They are gone; `pytest` from the root now
collects **2,699 tests with no errors**.

⚠ The first check — "is anything referencing them?" — was the wrong
question, and asking only it would have been an expensive mistake twice
over.

- It nearly deleted **`example10.py`**, which *is* referenced: two
  `.rst` files include it as a documented example. Grep for references
  across `.rst` and `.md`, not just `.py`.
- It would not have noticed that two of the nine were **prototypes of
  code that shipped**: `scratch_coupled.py` is the origin of
  `Transient._coupled_band`'s simultaneous `(dx, dh)` Schur solve, and
  `temp_nl.py` is the origin of `jaxtransient.py`'s TLine ring buffer
  and `interp_tlines`. Nothing was lost, but that was established by
  **reading all nine and tracing each idea into the tree**, not by
  grepping for their names.

"Nothing references it" proves a deletion will not break a build. It
says nothing about whether the content is worth keeping. Those are two
different questions and only the second one is about value.

The files had an active cost too: five of them set
`circuit.default_toolkit = jax_toolkit` **at module scope**, and pytest
imports every collected module into one process. `test_sym_dc.py` fails
alone with `SymbolicToolkit has no attribute 'insert'` and fails under
collection with a JAX dtype error, because a sibling had already
replaced the global toolkit. **A broken file that reports the wrong
reason costs more than a broken file.**

## 39. A bound around the overhead figure (2026-08-27): three of the four needed nothing

⚠ **The guard this section builds could not report a regression, and
that was not found until §48 (2026-08-31).** Everything below about
*which* figures deserve a guard stands; what did not stand is the claim
that the one guard built here would fail on decay. Read §48 before
trusting this section's closing paragraphs.

§38 proposed guarding every published performance figure. Measuring
first cut that to one, and the cutting is the useful half.

**What already existed.** `fold_card` and `_autohold` each have a
DIRECTION guard — `test_everything_outside_instance_is_folded` asserts
the folded class emits fewer primitives than the symbolic one, and
`test_each_regulariser_holds_its_argument` asserts the regularisers
hold. Neither optimisation can silently stop working, which is the
failure that actually matters. §38 said otherwise and was wrong.

**What was refused.** Strengthening those from `<` to `< 0.75 *`, to
pin magnitude as well as direction. On the Gummel-Poon card that test
uses, folding moves the primitive count **275 → 263** — a 4.4% margin.
A bound there pins noise. The 38% reduction that produces the published
1.54x lives on PSP, whose folded compile is **19.5 s**, too slow to
spend on every run. Recorded rather than built.

**What was built.** `test_perf_guards.py`: the end-to-end DSL overhead,
the one figure with no guard of either kind, because it is *call*
overhead and has no emitted-work proxy to assert on instead.

### What the measurement design cost, and why it is not a `time.time()` assert

| finding | number |
|---|---|
| Absolute wall time is unusable | the machine ranged **load 7 to 17** within the hour |
| An in-run ratio is not load-immune either | the benchmark reported **0.86x** at load 7 against 1.01x quiet — it times all of A then all of B, so drift lands entirely on B. Fixed by **interleaving** A B A B and taking the min of each |
| Shrinking the run is not a free saving | 1000 steps → min-ratio 0.966; 300 → 1.080; 150 → 1.131. Fewer steps is noisier **and biased upward**, because setup is amortised over less work |
| A single reading flakes | **1 failure in 27 runs**, ~4%. Enough to fire every dozen full-suite runs and earn the test a reputation for crying wolf |

The bound cannot simply be widened to absorb that: **1.22x is the state
being guarded against**, so there is nowhere above 1.15x to go without
making the test unable to fail for its own reason. Instead a breach is
**confirmed, not believed** — re-measured once, and the less damning
reading kept. Two independent excursions square the false-alarm rate; a
real regression is a level shift, not an excursion, and breaches both
times. The second measurement is only paid for when the first says
something interesting.

### The bite-check found the bug, as usual

The injected slowdown was too small — **1.096x against a 1.15x bound**
— so the harness could not have detected the regression it claimed to
guard. Calibrating it produced a finding worth keeping:

    1 dot of 2000 elements  ->  1.089x
    3 dots of 400 elements  ->  1.167x
    5 dots of 400 elements  ->  1.398x

**Twenty-five times less arithmetic, and it costs more.** The price is
the per-call numpy entry, not the flops — which is precisely why the
DSL's overhead is a call-count question, and precisely why it has no
emitted-work proxy and needed a timing test after all.

Both assertions are load-bearing together: equality without the ratio is
exactly what let `hdl.rst` drift from 1.14x to 1.25x with a green
benchmark in the tree, and the ratio without equality could be "won" by
quietly doing less work.

## 40. Lifting PCNR's internal-node constraint (2026-08-27): the dependence is affine, and PCNR already has the hook

### The problem, stated properly

Vector PCNR refuses any current that still reads a node voltage once the
declared `$limit` probes have been substituted: *"every resistive
current must reach the solution only through the declared probes"*.

**That rule is why PCNR does not work on the cards a PDK ships.**
Measured, on model cards rather than defaults:

| model | with defaults | with a parasitic |
|---|---|---|
| `DiodeSpiceHdl` | scalar | `rs=2` → **refused** |
| `GummelPoonNpnHdl` | vector | `rb=100` → **refused** |
| `MosLevel1Hdl` | vector | `rd=5` → **refused** |
| `MosLevel3Hdl` | vector | `rs,rd=5` → **refused** |
| `MesfetStatzHdl` | vector | `rs=2` → **refused** |

⚠ This was known and pinned: `test_pcnr_vector.py::test_the_refusals_
name_their_rule` asserts the `rb=10.0` refusal, and the PCNR tests set
the parasitics to zero **in their docstrings** -- *"at `rb = re = rc =
0` so that it qualifies"*. What was NOT recorded anywhere is the
consequence: eligibility counts quoted from default cards (§15's "2 of
26", and "13 of 37" measured this session) both overstate it. **A
default card zeroes the parasitic that does the refusing** -- the same
blindness §28 and §33 already recorded, in a third guise.

### The measurement that says it is liftable

The rule is right for a NONLINEAR node dependence and too strong for a
linear one: a series resistance contributes `g * (x_a - x_b)`, an
ordinary conductance that Newton solves exactly and that limiting would
do nothing for.

Instrumenting the detector and capturing every surviving expression:

    I(c)     = _pcnr_v1*area/rc + _x0*area/rc - _x3*area/rc
    var(dT)  = _x3 - _x4
    I(d)     = (_pcnr_v0 - _pcnr_v1 + _x0 - _x1)/_v91_rdx

**Nine survivors across the five refusing models, and every one has an
identically zero second derivative in every node symbol it reads.** Not
one coefficient depends on a probe. The shapes are conductances.

### The hook that already exists

`augmented_system` does not merely ignore a participating device -- it
**shadows its `i`/`G` with `_zero_vector`/`_zero_matrix`** while it
re-stamps the device at `v_lim` (`pcnr.py:428-431`). So the affine
remainder has a natural home and needs no new assembly path: **shadow
with the remainder instead of with zeros.** The linear part is then
stamped by the ordinary MNA machinery exactly as it is on the
`pcnr=False` path, and PCNR carries only the probe-dependent part.

### Built so far

`hdl._affine_in_nodes` -- the split, with the two refusals that keep it
honest:

* `d/dx` still reading a node ⇒ not affine ⇒ refuse (`x*x`, `exp(x)`);
* a coefficient reading a PROBE ⇒ bilinear `v_lim * x` ⇒ refuse. That
  is a conductance that MOVES with the limited unknown, and the ordinary
  assembly is built before `v_lim` is known. No measured case has one.

The refusal still stands and now names the class
(`[affine remainder available]`), so the five models above are
identified without any behaviour change. `test_pcnr_affine.py`, 15
tests, bite-checked: with the detection disabled, **9 fail**.

### What remains

1. Compile the remainder as a stampable `i`/`G` pair.
2. Shadow with it in `augmented_system` instead of zeros.
3. Accept the device: return the spec rather than the refusal.
4. Validate: `GummelPoonNpnHdl(rb=100)` under PCNR must reach the same
   operating point as `pcnr=False`, and the diff-pair clash case must
   stay fixed. **This is where a mistake would be silent**, so the gate
   is agreement with the unlimited path, not convergence alone.

⚠ Two of the seven non-qualifying nonlinear models are NOT this class
and must not be swept in: `PhotodiodeHdl` reads optical drive nodes
(`optn`, `optp`), and `LedHdl`/`DiodeSpiceThermalHdl` carry a
branch-current unknown. The bite-check test asserts they are not
claimed.

### ⚠ THE "CONSUMPTION IS WRONG" VERDICT BELOW WAS ITSELF WRONG

**Corrected same day. The lift works and is ON.** All three liftable
models reach the unlimited path's operating point to **1.5e-18,
1.1e-16 and 2.8e-17** relative. `PCNR_LIFT_AFFINE = True`.

The cause of the wrong verdict is worth more than the verdict was.
`_pcnr_declared_route` returns a spec dict; `generate_code` then
**REBUILDS** it as `pcnr_vec = dict(probes=..., kinds=..., ...)` rather
than passing it through. A field added to the route's return and not
added to the rebuild is dropped in silence -- so the lift ran, computed
`lin_terms`, and the class saw `None`. PCNR then zeroed the resistor and
never re-stamped it, which is the 0.5 V error.

I diagnosed that as "the lift is wrong" and wrote it into this document.
It was one missing line in a dict literal. **A field that must be added
in two places will eventually be added in one**; the rebuild now carries
a comment saying so.

Three further things the working version needed, each found only by
running it:

* the coefficients must be **back-substituted through the chain** before
  compiling. They are not parameter expressions: the measured ones are
  `1/_v91_rdx` and `area/rc`, the first a chain variable. Safe because
  `_affine_in_nodes` has already proved no coefficient reaches a probe.
* `paramsyms` is a `generate_code` local and does not exist in the
  metaclass `__init__`, so the compile has to happen at the former.
* the coefficients must be evaluated with **numpy scalars, not python
  floats**. A compiled coefficient is a `select` over both arms and the
  unselected arm is still evaluated: MOS level 1 emits
  `select([rd > 0, True], [1/rd, 1/(nrd*rsh)])` and `nrd` defaults to 0.
  Python raises `ZeroDivisionError`; numpy gives `inf`, which `select`
  discards. The same both-arms hazard `hdl.select` exists for (sec. 21),
  met in a new place.

And one more, found only by running the FULL suite after enabling it:

* **a collapsed variant must not inherit the conductance block.**
  Card-collapse forks a subclass, so setting `pcnr_lin_G` only on the
  classes that have a remainder leaves the others inheriting one sized
  for the parent's geometry -- a `(6,6)` block dotted into a four-node
  variant's local voltages. Three tests failed at once with
  `shapes (6,6) and (4,) not aligned`, including a cascode convergence
  case that looked like a solver regression and was not. It is now set
  on every class that reaches the vector route, `None` included, and
  pinned by a test.

  ⚠ The narrow validation is the lesson. Three single-device DC circuits
  agreed to 1e-18 and I read that as "the lift works". It worked for
  three geometries out of the ones the library builds. **An agreement
  gate proves the cases it runs, and card-collapse silently multiplies
  the cases.**

Final: **2633 passed / 3 skipped / 3 xfailed / 0 failed**, lift ON.
`MosLevel1Hdl`, `MosLevel3Hdl` and `MesfetStatzHdl` now take PCNR with
their series resistances declared -- cards a PDK would actually ship.

The original (wrong) entry is kept below, because the SYMPTOM it records
is the valuable part and is exactly what the gate is for.

### The lift, first attempt (2026-08-27): the gate caught a silent wrong answer

Scaffolded behind `hdl.PCNR_LIFT_AFFINE` (default **False**, so nothing
in the library changes). With it on, the three liftable models stop
being refused -- **and the answer is wrong.**

    vd=2.0, vg=1.5, MosLevel1Hdl(rd=5, w=10u, l=1u)

      ordinary   internal drain node  1.998874596 V
      PCNR lift  internal drain node  1.498874596 V

Out by **exactly the gate voltage**: the drain resistor is zeroed out of
the ordinary assembly and never re-stamped, because `lin_terms` arrives
`None` and no `pcnr_lin_G` is attached. Detection (`_affine_in_nodes`,
independently tested) is correct; consumption is not.

⚠ **The device's terminal CURRENT still agreed to seven digits**
(-2.2508079e-04 against -2.2508072e-04) **and the solver converged
without complaint.** Only the operating-point comparison saw it. That is
the entire argument for §40's gate being *agreement with the unlimited
path* rather than convergence -- a limiter bug does not announce itself
by failing to converge, it announces itself by converging somewhere
else. Pinned as a `strict=True` xfail so it cannot be forgotten and
cannot silently start passing.

### Two corrections this work forced

**1. "Not one coefficient depends on a probe" was wrong.** That was
measured with a surface check on the coefficient's own symbols. A
coefficient is typically a CHAIN VARIABLE -- `1/_v68_rbb` -- whose
definition may itself read a probe. Checked transitively,
`GummelPoonNpnHdl`'s base resistance **is** probe-dependent: Gummel-Poon
models current crowding, so `rbb` moves with bias. `g*(x_a - x_b)` with
`g` moving per iteration is bilinear, not a constant conductance, and it
is correctly refused. Liftable set: **3 models, not 5** -- `MosLevel1`,
`MosLevel3`, `MesfetStatz`.

**2. A chain-variable survivor is never liftable**, even when affine.
`GummelPoonNpnThermalHdl` gives `var(dT) = _x3 - _x4`: affine in the
thermal nodes, and `dT` then drives every temperature-dependent quantity
exponentially. Splitting a DEFINITION linearly says nothing about how it
reaches the current. Only a survivor in a CURRENT is a conductance.

### ⚠ And a real defect in the compile cache, found by being bitten

`PCNR_LIFT_AFFINE` is a RUNTIME flag that changes what is compiled. The
cache key covers `hdl.py`'s SOURCE (via the dependency hashes), so an
edit invalidates -- but flipping a flag in a running process leaves
every file byte-identical. **The two builds shared a key.**

One probe with the lift on wrote **924 entries** that then served the
lift-off suite, and the symptom was `MosLevel1Hdl(rd=5)` no longer being
refused, with an *empty* refusal string, in a process where the flag was
False. It read exactly like a bug in the detector and cost a confusing
half hour.

Fixed: `'pcnr_lift=%r'` is now part of `_hdl_cache.key_for`, with a test
asserting the two keys differ and that the off-key is stable. **The
general rule, which this cache did not previously state: anything that
changes emission and is not a file must be in the key.**

## 41. Gummel-Poon, addressed (2026-08-27): `var()` was hiding the cancellation

§40 refused `GummelPoonNpnHdl` because its base-resistance coefficient
`1/_v68_rbb` reads a probe. That was measured with a check that follows
chain REFERENCES, and it over-refuses. SPICE's base resistance is

    rbmx = var(Piecewise((rb, rbm < 0), (rbm, True)))
    rbb  = var((rbmx + (rb - rbmx) / qb) / area)

where `rbm < 0` is SPICE's sentinel for *not given*. When it is not
given, `rbmx = rb`, `(rb - rbmx)` is an **exact zero**, the whole `qb`
term vanishes, and **`rbb = rb/area` is a constant**.

⚠ **sympy cannot see that while `rbmx` is a held `var()`.** It reads
`25.0 - _v67_rbmx`, which is not zero. This is the same "`var()` is
arithmetic" trap the DSL notes elsewhere, arriving as a *false refusal*
rather than a wrong number: the compiler could not prove a conductance
was constant because the model had helpfully named it.

**Fix: resolve the chain and ask what survives**, rather than following
references. `_reads_a_probe` now substitutes each chain definition
(reverse order, one pass) before testing. The two answers are genuinely
different physics, not two spellings:

| card | `rbb` | verdict |
|---|---|---|
| `rbm` folded to its default | `rb/area`, constant | **liftable** |
| `rbm` left symbolic | carries `qb` | correctly **refused** |

So `fold_card` is what makes a BJT eligible -- which is what anyone
shipping a production model does anyway, and is exactly what §30 built
folding for.

**Measured:** a folded Gummel-Poon common-emitter stage under PCNR
reaches the ordinary path's operating point to **2.220e-17**. The
unfolded card is still refused, and a test holds both.

Bite-checked: reverting `_reads_a_probe` to reference-following fails
the two new tests and passes the other 24, so they discriminate on
exactly the change that matters.

### Liftable set, current

| model | needs |
|---|---|
| `MosLevel1Hdl`, `MosLevel3Hdl`, `MesfetStatzHdl` | nothing -- `rd`/`rs` cards work directly |
| `GummelPoonNpnHdl` | `fold_card` (so `rbm` resolves) |
| `GummelPoonNpnThermalHdl` | still refused: chain-variable survivor, `dT` drives physics exponentially |
| `PhotodiodeHdl`, `LedHdl` | still refused: optical drive nodes / branch-current unknown |

### `rbm > 0`: the refusal is on the arithmetic, not the sentinel (2026-08-27)

§41 said folding makes Gummel-Poon liftable "because `rbm < 0` means
unset". That is the common case, and stating it that way would have
been a rule about a *sentinel*. The check is on whether the arithmetic
cancels, which is a stronger and more useful thing. Measured, card
`rb = 25`:

| `rbm` | `rb - rbmx` | verdict | why |
|---|---|---|---|
| unset (−1) | **0** | liftable | SPICE's "not given" ⇒ `rbmx = rb` |
| `25.0` | **0** | **liftable** | numerically equal to `rb` -- cancels just as exactly |
| `10.0` | 15 | refused | a real high-injection card; `qb` survives |
| `0.0` | 25 | refused | `qb` survives |

**A sentinel-based check would get row 2 wrong** and quietly refuse a
card it could take. Resolving the chain and asking what survives gets
all four.

⚠ **And the refusal is guarding something.** If `rbm = 10` gave the same
answer as unset, refusing it would prove nothing and the liftable tests
would pass for a lift that was simply broken. It does not: the operating
point moves by **1.07e-3**. Meanwhile `rbm == rb` is **bit-identical**
to unset -- 0.0 difference, not merely small -- which is what the whole
lift rests on. Both are asserted, because "refused" and "liftable" are
only meaningful if the two cards actually differ and the two identical
ones actually agree.

The refused cards lose nothing but PCNR: the ordinary path solves them
normally, and `pcnr.solve_dc` raises a plain "no device in this circuit
declares a PCNR junction" rather than degrading silently.

## 42. High injection: solved by a probe, not by more solver (2026-08-27)

§41 left `GummelPoonNpnHdl` liftable only on cards where `rbb` happens
to be constant. On a real high-injection card (`rbm < rb`) it is not,
and the term is **bilinear** -- `g(v_lim) * (x_a - x_b)`, a conductance
that changes every iteration, which the ordinary MNA assembly cannot
carry because it is built before `v_lim` is known.

### The design that was NOT built

Extending the lift to bilinear terms is tractable and was scoped: the
shadow closes over `g` evaluated at the current `v_lim` (that part is
nearly free, and gives `J_mm` correctly), and `stamp` gains the extra
`dg/dv * (x_a - x_b)` contribution to `J_ml`, which needs `x` at stamp
time -- available, since `augmented_system` already reads `x[ra]`.
Perhaps a day, plus the convergence risk of a new coupling.

### The one-line change that was

**Declare the base-resistance branch as a `limit_identity` probe.**

    Contribution(brb.I, limit_identity(brb.V) / rbb)

`brb.V` was a raw node difference, which is exactly what PCNR's rule
refuses. As a declared probe it is an unknown like any other, so the
current reaches the solution legitimately -- and `rbb`'s dependence on
`qb` stops mattering, because `qb` is itself a function of the two
junction probes. **The device qualifies on ANY card, folded or not,
whatever `rbm` says.** Cost: one PCNR unknown per device.

`limit_identity` carries no law, so it must not move a number. Measured:
**bit-identical `i`/`G`/`q`/`C` across 25 random biases for all three
Gummel-Poon classes -- max relative difference exactly 0.000e+00.** The
recorded-digest test corroborates independently: its POINT digests are
unchanged and only the EXPLAIN text moved, which is the separation that
test exists to make.

Operating-point agreement under PCNR: **2.5e-17** (`rbm` unset) and
**1.78e-16** (`rbm = 10`, the high-injection card the affine remainder
could not take at all).

### The two mechanisms compose

Gummel-Poon now carries an identity probe AND a conductance block: the
probe handles `rbb`, which is bias-dependent, and the affine remainder
handles `rc` and `re`, which are plain constants. Each takes the
resistance the other cannot.

⚠ **The lesson is about reaching for the solver.** The bilinear
extension was a real design and I had it costed before asking whether
the model could simply declare what it already knew. A limiting
framework's rule -- *every resistive current must reach the solution
through a declared probe* -- is more usefully read as an instruction to
the MODEL than as a constraint on the solver.

Records updated: five tests asserted the old refusal and are inverted
rather than deleted, with the superseded assertion kept in a comment;
two probe-list assertions gained the entry; three explain digests
re-recorded with the bit-identical evidence beside them.

## 43. Self-heating and the series branch take PCNR too (2026-08-27)

§42's one-line trick generalises. Three more declarations, each turning
a raw node difference into a declared unknown:

| where | was | now |
|---|---|---|
| `SelfHeating.__init__` | `dT = var(b.V)` | `dT = var(limit_identity(b.V))` |
| `_spice_diode` | `vd = var(bd.V)` | `vd = var(limit_pnj(bd.V, isT, n*vtT))` |
| `_spice_diode` | `Contribution(brs.I, brs.V * area/rs)` | `limit_identity(brs.V)` |

**`dT` is the deepest dependence there was.** It is not a parasitic on
the side: `dT` sets `Tj`, and `Tj` enters *every* current through
`exp(v/(n*k*Tj/q))`. Self-heating therefore disqualified a device for
the temperature the whole model is evaluated at.

**The SPICE diode's `rs` was refused for a reason that is not the
resistor.** `pdiss = Branch(a, c).V * idio` -- the TERMINAL voltage,
spanning both the series resistance and the junction. Declaring the
series branch lets the spanning tree resolve `x_a - x_c` into probe
symbols. The junction itself now takes `limit_pnj` rather than being
detected by the scalar route, which is the better mechanism anyway
(§15 built the vector route to fix the inter-device clash).

### Measured

**Bit-identical `i`/`G`/`q`/`C` on all nine affected model/card
combinations** -- `DiodeSpiceHdl`, `DiodeSpiceThermalHdl`,
`GummelPoonNpnThermalHdl`, `RThermalHdl`, at several cards each. The
recorded-digest test corroborates: three EXPLAIN digests moved, every
POINT digest held.

Operating-point agreement under PCNR, against the unlimited path:

    DiodeSpiceHdl rs=2            0.000e+00
    DiodeSpiceHdl plain           0.000e+00
    DiodeSpiceThermal rth=200 rs=1 0.000e+00
    GummelPoonNpnThermal rth=100  0.000e+00

### Still refused, and correctly

At the default `rth = 0` the thermal branch **collapses to a zero-volt
source**, leaving a genuine branch-current unknown. That is the
isothermal limit, where a thermal unknown means nothing. Pinned by a
test so a later relaxation of the branch-unknown rule has to think
about the case rather than sweep it in.

A thermal diode at 0.8 V with `rth = 200` and no series resistance does
not converge on EITHER path -- an actual thermal runaway, not a solver
failure. PCNR does not rescue it and should not be claimed to.

### Records

Seven tests recorded the old refusals and are inverted rather than
deleted, each keeping its superseded assertion in a comment: two
`explain` refusal strings, a probe-kind list, a `limit_spec` list, and
`test_the_full_spice_diode_qualifies_only_without_rs` -- whose NAME was
the claim, and which had already expired once before under a different
reason. Three EXPLAIN digests re-recorded with the bit-identical
evidence beside them.

## 44. The branch-unknown gate was too strong (2026-08-27)

`_pcnr_declared_route` opened with

    if states or vbranches:
        return None, '... the device carries ...'

-- a refusal on the PRESENCE of a V-contributed branch, not on whether
it interacts with anything. The reasoning behind it is sound-sounding: a
V-branch adds a current unknown to `x`, and that row's equation is a
voltage residual rather than KCL, so stamping a current into it would be
wrong.

**Measured, it is not wrong.** `pcnr_i` is compiled from the same
`i_exprs` the ordinary assembly uses, so it reproduces the voltage
residual on that row exactly as the element would have. Relaxed:

| device | agreement with the unlimited path |
|---|---|
| self-heating BJT at `rth = 0` (collapsed thermal branch) | **9.025e-17** |
| `LedHdl` (carries a V-branch outright) | **0.000e+00** |

Full suite with the gate relaxed: **2644 passed, and not one numerical
failure** -- the only three were tests asserting the refusal itself.

### ⚠ One of those three was mine, written an hour earlier

`test_the_isothermal_limit_is_still_refused_and_should_be` said the
`rth = 0` refusal was *"correct rather than incidental: it is the
isothermal limit, where a thermal unknown means nothing"*. A fine
argument, and not a measurement. Solving the circuit gives 9.025e-17.

The collapsed branch does pin its probe to zero, so it is one unknown
that cannot move -- which makes it **useless**, which is what I
reasoned, and not **harmful**, which is what I asserted. Only the second
justifies a refusal. The test is inverted, with that written in it.

### What was NOT relaxed

`states` -- an `idt`/`idtmod`/`laplace` state is device MEMORY, and
nothing in PCNR has been measured against one. Pinned separately so the
next relaxation has to argue for itself.

⚠ And the evidence covers the shapes this library has: a collapsed
thermal branch, and an optical-drive LED. A V-branch doing something
else is untested, not blessed.

### Where PCNR stands after §40-44

**16 of 37 models eligible at defaults**, against 13 this morning and
about 4 on cards a PDK would actually ship. Exactly one model with an
exponential is still refused: `PhotodiodeHdl`, whose `iph` reads its
optical drive nodes -- an INPUT rather than a parasitic, and the one
remaining case that is about physics rather than plumbing.

The arc from §40: the constraint was never really "PCNR cannot do this".
It was three separate places where a model handed the solver a raw node
difference and the solver, correctly, refused to guess. Two lines of
model declaration and one over-broad gate.

## 45. The photodiode, and a gate that nearly hid its own point (2026-08-27)

The last model with an exponential that PCNR would not take.
`iph = var(resp * bopt.V)` reads the optical port -- a raw node
difference -- so the device was refused. Declared:

    iph = _var(p.resp * limit_identity(bopt.V), 'iph')

**16 of 37 models eligible becomes 17, and every model with an
exponential now qualifies.**

### ⚠ I investigated this and reported the wrong conclusion

Asked for both sides of the decision, I wrote that there was **no
demonstrated convergence benefit**: a lit cell at 1, 10 and 50 suns into
10 ohm to 10 M, and 4- and 12-cell series strings in both instance
orders, all converge on the ordinary path and agree with PCNR. On that
basis I recommended *not now*.

The recommendation was answered "add it", and adding it produced the
measurement I had not made. My agreement gate -- the same 1e-9 used for
every other model -- **failed** on three of the five circuits, differing
by about 3e-9. The first reading was "PCNR is wrong here".

It is the other way round:

| circuit | ordinary \|f\| | PCNR \|f\| |
|---|---|---|
| 1 sun, 100k (near Voc) | 2.933e-09 | **5.424e-16** |
| 50 suns, 10M | 3.113e-09 | **6.545e-15** |
| 4-cell string | 4.231e-09 | **7.932e-16** |
| 12-cell string | 3.007e-09 | **5.940e-15** |

Five to six orders of magnitude, on a Jacobian conditioned at **41** --
so the residual maps almost one-to-one onto voltage error. The ordinary
path stops at its tolerance; PCNR drives the residual to machine
precision. **The "disagreement" was the REFERENCE being loose.**

The investigation asked *did both converge* and *did they agree*, and
never asked **which one was right**. `reference-differential` has the
rule already -- *a residual is not an error, and root error is
`residual / |dF/dx|`* -- and I applied it to neither side. An agreement
gate silently assumes one party is correct.

So this device is gated on its **residual**, not on agreement with the
ordinary path, which is the honest comparison when neither side is
known-correct in advance.

### The caveat that still stands

This is the one declaration that is a boundary CONDITION rather than a
device quantity. Gummel-Poon's base resistance, `SelfHeating`'s junction
temperature and the diode's series branch are things the device HAS;
`V(optp,optn)` is an input someone drives. Recorded in the model, because
"declare it" stops discriminating if it is the answer to every refusal.

Bit-identical `i`/`G`/`q`/`C` on both cards tested, as an identity probe
must be.

## 46. Which models need limiting -- the criterion, measured (2026-08-27)

§45 left 20 models PCNR-refused and the obvious question: is that
correct, or merely unexamined? The `states` gate in particular refused
on the same blanket-presence basis the V-branch gate did, and that one
had already been wrong.

**It is correct.** But the two criteria that look obviously right are
both wrong, and that is the part worth keeping.

### ✗ Does the current overflow?

`DiodeHdl` leaves float range at **18.35 V** -- the classic `709*V_T`
signature, and well inside what a Newton iterate reaches. But
`DiodeSpiceHdl` and `GummelPoonNpnHdl` **never overflow at any
voltage**, because the DSL's `expl` already bounds them, and they
plainly still need limiting.

Overflow measures whether the DSL has a regulariser, not whether the
device needs a limiter.

### ✗ How many decades does `|G|` span?

A **false positive**. `ChargePumpHdl` spans **8.1 decades** over a
+/-1 V sweep -- more than `EkvNmosHdl` (6.8) or `PhotodiodeHdl` (4.7),
both eligible. Its nonlinearity is a `tanh` with a narrow switching
width, so the range comes from the TAILS going tiny: its max `|G|`
actually **falls**, 1e-3 to 8e-12.

Dynamic range conflates *G becomes huge* with *G becomes small*. Only
the first is a problem for Newton.

### ✓ The GROWTH LAW of `max|G|`

Limiting exists because Newton steps into a region where the
linearisation is worthless, and that is when the conductance grows
**super-polynomially**. Measured, sweeping every row, `max|G|(10 V) /
max|G|(1 V)`:

| population | ratio |
|---|---|
| 13 eligible models (`Diode`, both `GummelPoon`, both `MosLevel`, `MesfetStatz`, `Led`, `Photodiode`, the thermals) | **1.95e+87** (`DiodeHdl` 1.56e+151) |
| every refused model | **7.92e+03 and below** -- `Memristor` 7.9e3, `Mixer` 1.0e1, the rest at 1.0 |

**Eighty-four orders of magnitude between the populations.** The
threshold is not a tuned parameter; anything in six decades of slack
gives the same answer.

Built: `benchmarks/limiting_need.py` (the table, with both wrong
criteria recorded) and
`test_no_model_that_needs_limiting_is_refused_by_pcnr`. Bite-checked by
removing the photodiode's probe: the guard names the model, its ratio
**1.9468e+87**, and the refusal text.

### ⚠ The criterion holds in ONE direction only

A low ratio is **not** proof a device is safe. `EkvNmosHdl` reads
**3.11e+01** here -- bounded, like the refused models -- and it
genuinely needs limiting: a probe from zero bias never reaches its
exponential region. `MesfetCurticeHdl` (1.6e1) and `HemtAngelovHdl`
(1.0) are the same.

So the guard asserts *everything that visibly needs limiting is
eligible*, which is the half that can fail usefully. The conclusion
that the 20 refusals are correct rests additionally on their
nonlinearities being **structurally bounded** -- `tanh`, `maxc`/`minc`,
products -- rather than merely un-probed.

⚠ And the first version of the probe swept only row 0, which made
`GummelPoonNpnHdl` and `MosLevel1Hdl` read **flat**, because their
exponential lives on an internal junction. That nearly went into a
report as data. Both the benchmark and the test sweep every row, and
say why.

## 47. `pcnr=True` is a request; `pcnr_status` is what happened (2026-08-27)

`DC(pcnr=True)` could quietly not use PCNR, in two different ways, and
the caller had no way to tell.

**Defect 1 -- the attribute did not exist when it was not interesting.**
`self.pcnr_fell_back = False` was assigned *inside* `if self.par.pcnr`,
so a caller checking `dc.pcnr_fell_back` after an ordinary solve got an
`AttributeError`.

**Defect 2 -- the flag could not separate "ran" from "did nothing".**
`pcnr_fell_back` reads **False** both when PCNR solved the circuit and
when `pcnr=True` was asked for and *no device declared a probe*.
Measured:

    has a diode    participants=1   pcnr_fell_back=False
    linear only    participants=0   pcnr_fell_back=False

⚠ That is the **silent fallthrough this module's own comment records
fixing** for the MOSFET differential pair (Stage 2, §15): participation
used to be gated on the pnj-only pair view, which is empty for pure
`fetlim`/`limvds` devices, so `pcnr=True` fell through without saying
so. The GATE was fixed then. The REPORTING was not, and the same
failure came back through the other door.

**Fixed:** `pcnr_status`, always set, one of

    'off'              the parameter was False
    'used'             PCNR solved this
    'no-participants'  asked for, no device declares a probe
    'fell-back'        PCNR raised; ordinary Newton solved it instead

and the no-participants case now `logging.warning`s. `pcnr_fell_back`
is kept -- three tests read it -- and still means what it meant.

Bite-checked: reverting the reporting fails three tests, one of them on
the original `AttributeError`.

### `Transient` -- the same gap, closed with a counter

`transient.py:2062` carried the identical fallback and recorded
nothing. It differs in a way that decides the shape of the fix: PCNR is
attempted **per timestep** and can fall through on any of them, so
there is no single state to report. A transient can lose PCNR on some
steps and leave no trace but log lines.

So it counts:

    pcnr_solves      timesteps PCNR carried
    pcnr_fallbacks   timesteps that fell through to the ordinary solver
    pcnr_status      'off' | 'used' | 'no-participants'
                     | 'partial' | 'fell-back'

**`'partial'` is the case the counter exists for**: a run that loses
PCNR on some steps is neither 'used' nor 'fell-back', and before this
it looked exactly like a clean PCNR run from outside.

⚠ The counts are **solver invocations, not accepted steps**. A rejected
step is solved and then thrown away, and it is still a step PCNR did or
did not carry -- conflating the two would make the number disagree with
`statistics.accepted_steps` for a reason nobody could reconstruct.

**SETTLED 2026-08-31.** This was raised as an open trade-off and put to
the branch author, who ruled: *keep counting solver invocations, that
is the honest number.* So the consequence is intended, not an
inconsistency awaiting repair -- `pcnr_solves + pcnr_fallbacks` will
generally **exceed** `statistics.accepted_steps`, and anyone reconciling
the two is removing information rather than adding it. The alternative
(counting accepted steps only) is **refused**: it hides solver work that
actually happened, which is the whole reason the counter exists.

Measured on a diode transient: `pcnr=True` gives `used`, **58 solves, 0
fallbacks**; the same circuit with `pcnr=False` gives `off, 0, 0`; an
all-linear circuit gives `no-participants`. Bite-checked by injecting a
failure on every third PCNR step -- status `partial`, both counters
non-zero -- and by disabling the counters, which fails three tests.

## 48. The one magnitude guard could not report a regression (2026-08-31)

Section 39 built `test_perf_guards.py` to convert silent performance
decay into a failing run -- the branch's ONLY magnitude guard, since
pinning `fold_card` and `_autohold` was measured and refused as
noise-pinning.  Re-measuring the suite (sec. 9 of the reviewer's guide)
showed it skipping on some runs and not others, which was the thread
worth pulling.

**It could not fail.** `_require_a_still_machine` was consulted only
inside the `if breached:` arm of `_confirmed_ratio`, so the sequence for
a real regression was breach -> gate -> abstain: on any machine too busy
to measure, the one outcome the file exists to report was the one it
could not report.  Demonstrated by driving the guard's own path with the
file's calibrated `_SlowR` -- `Skipped` where a failure was due.

⚠ **The bite-check could not have caught this, and looked like it
should.** It drives `want='above'`, where a successful bite means
breached=False and returns BEFORE any gate is consulted.  It exercises a
different path from the one a real regression takes, so "the guard has
been shown to bite" was true and did not mean what it appeared to.

**And the comfortable side was never checked at all**, on the recorded
argument that "a reading that is already comfortable needs no
defending".  That argument is wrong: the bias runs both ways.  Measured
at the then-current 3 samples a side under load, a true ~1.01x parity
read anywhere in **0.916x to 1.054x** -- so a 1.20x regression reads
1.09x and passes, silently, exactly when the machine is least able to
see it.

### Two candidate fixes, refused with numbers

**CPU time as a load-immune instrument -- refused.** Wall clock and
`time.process_time` recorded simultaneously on identical work under
four concurrent openEMS runs: drift `2.371 / 1.149 / 3.353` wall against
`2.371 / 1.149 / 3.353` CPU.  Identical to three decimals.  Contention
here does not deschedule the process, so CPU time is not a different
instrument.

**Fewer samples for the coarse bite-check -- refused, same day, by its
own first run.** The bite-check resolves a ~1.40x injected slowdown
rather than a 15% effect, so it was given 3 samples on the reasoning
that a coarse effect needs fewer.  It skipped immediately: two
measurements disagreed by **1.147x** against a 1.05x limit.  The
reasoning conflated separating an effect with reproducing an estimate.
The gate constrains ESTIMATOR PRECISION, which at 3 samples is +/-14%
however large the effect is.  Both callers use the same sample count.

### What the fix is

**Gate on reproducibility, not on stillness.** `drift` -- the spread of
raw timings inside one measurement -- measures the MACHINE, and the
question is whether the ESTIMATE reproduces.  They disagree completely:
raw spread read 1.86x-3.35x in trials whose min-of-12 estimates agreed
across trials to 1.5%.  Gating on drift would have abstained on every
one of them.

**Readability is decided first and for every outcome**, so a comfortable
reading taken on an unreadable machine is an unread instrument rather
than a pass.

**Retry rather than abstain, bounded at four measurements.** Measured
disagreement between consecutive N-sample estimates under load 17:

| N | consecutive-pair disagreement | within 1.05 |
|---|---|---|
| 3 | 1.098 1.020 1.163 1.045 1.104 1.049 1.118 1.020 ... | 8 of 16 |
| 6 | 1.131 1.084 1.025 1.054 1.009 1.058 1.013 1.022 | 4 of 8 |
| 12 | 1.034 1.049 1.029 1.004 | 4 of 4 |
| 24 | 1.081 1.032 | 1 of 2 |

Agreement improves with N but never becomes reliable -- dependable 1.05
would cost ~48 samples a side, minutes per test.  A busy machine has
quiet patches, so a second or third attempt usually lands in one, and
the cost is paid only when the first pair disagrees.  `PAIRS` rose 3 ->
12 on the same table.

⚠ Agreement is judged on the **two most recent** estimates, never the
closest two of the set.  Shopping a growing set for its friendliest pair
would find agreement most easily when both readings are wrong the same
way.

**Abstention is loud.** It warns before it skips, so an unchecked figure
lands in the warnings summary of a plain run rather than only behind
`-rs`, and the reason names both numbers: the disagreement is why the
reading was dropped, the drift says whether the machine was why.

### The part that is testable without a machine

The decision is now a pure function of `(ratio, disagreement)`, so the
paths a real regression takes are covered on any machine at any load, in
under a second: **16 tests**, against 2 before.  That matters beyond
convenience -- the reason the original hole survived review is that
every existing test needed a quiet machine to say anything, and this
machine has not been quiet for hours.

Measured after the fix, three consecutive runs at load 13.8-15.7: **18
passed, 18 passed, 18 passed**, no abstentions, at 92 s / 170 s / 85 s.
The 170 s run is the retry earning its place rather than a cost.

⚠ **The suite pays for this.** The file went from 2 tests to 18 and from
roughly 30-45 s to 85-170 s.  `PAIRS` and `MAX_MEASUREMENTS` are the
knobs if that is judged too expensive; lowering `MAX_DISAGREEMENT` is
not, since 1.05 is what makes an estimate precise enough to assert a
1.15x bound against a true 1.01x.

## 49. Vector PCNR Stage 3 (JAX): opened 2026-08-31, and scoped by probing rather than by reading

Stage 3 was parked by the branch author after Stages 1 and 2 landed on the
CPU. Opened at their request. This section is the plan, and the scope below
is **measured**, because the one-line description in the status table turned
out to be true but far too coarse to plan from.

### What the traced path ALREADY has (do not rebuild it)

P19 reversed the parity review's "one-sided by design" verdict and ported
PCNR to JAX. Present today in `jaxtransient.py`:

* `PcnrState` (`x`, `v_lim`, `iters`, `converged`) and `pcnr_inner_loop` — the
  limited quantities live in the traced state, which is the insight that made
  the port possible: classic limiting mutates per-device Python `_vlim`, which
  a `lax.while_loop` cannot host, while PCNR's junction unknowns are just more
  state vector.
* `pcnr_controller_jacobian`, the rank-k Schur update as four scatter-adds per
  junction, and the CORRECT phase as a branchless `pnjlim` over `v_lim`.
* `_junction_arrays`: static trace-time metadata, and a setup-time refusal of
  charge-storing junction devices, mirroring the CPU's `augmented_system`.

### What is missing — four gaps, each measured on 2026-08-31

**1. The unit is the wrong one.** `_pcnr_setup` gates on `pcnr_junctions`,
which is an ALIAS for `pcnr_junction_pairs` — the pnj-only pair view. Vector
PCNR's unit is the DEVICE (`pcnr_devices`, each carrying `m` probes). This is
exactly the distinction Stage 2 turned on: `transient.py` records that the pair
view "is empty for a circuit of pure `fetlim`/`limvds` devices, so `pcnr=True`
on a MOSFET differential pair used to fall through to the ordinary solver
SILENTLY". Measured on a `MosLevel1Hdl`: `pcnr_devices` → 1 device with
**m = 3**; `pcnr_junction_pairs` → 2 pairs. Different objects, different laws.

**2. A signature mismatch is what actually raises today.** `_junction_arrays`
builds per-junction scalar callables by calling `pi(v, params, epar, tk, jn=jn)`
and taking `[0]`. The HDL-generated device exposes
`pcnr_i(v, params, epar, toolkit)` — **no `jn` keyword** — which takes the
device's whole probe vector and returns the element's local row vector
(shape `(6,)` for `MosLevel1Hdl`). So the current refusal reads

    junction 0 of 'M' could not be evaluated for the traced backend
    (pcnr_i() got an unexpected keyword argument ...)

⚠ That message names a keyword-argument error, which invites someone to
"fix the signature". **Do not**: the mismatch is a symptom of gap 1. The
per-junction scalar callable is the wrong shape to ask a vector device for.

**3. Mixed limiting kinds.** `MosLevel1Hdl.pcnr_probes` is
`(((1,5),'fet'), ((4,5),'vds'), ((3,5),'pnj'))` — three probes, three
DIFFERENT laws. The traced CORRECT phase applies `pnjlim` only. `fetlim` and
`limvds` must join it, branchlessly (no host control flow, so `jnp.where` over
a per-probe kind mask assembled at trace time).

**4. The affine remainder (§40).** On the CPU the surviving affine dependence
is stamped by the ordinary MNA assembly — `pcnr.py` shadows a participating
device's `i`/`G` so the shadow IS the remainder — and PCNR carries only the
probe-dependent part. The traced assembly needs the same split, or every
device with a series resistance loses PCNR here exactly as it did on the CPU
before §40.

### Gates, declared before the work

* **G1 — it does what Stage 2 did.** The MOSFET differential pair that took
  Stage 2 from `[14, FAIL]` to `[22, 22]` on the CPU, run under
  `JAXTransient(pcnr=True)`: converges from both orders where plain Newton
  fails. Same circuit, same measure.
* **G2 — agreement, not convergence.** §40's lesson: "a convergent PCNR answer
  is not a correct one" — the first affine lift converged happily with an
  internal node out by a whole gate voltage and terminal current still matching
  to seven digits. Gate on agreement with the UNLIMITED path (CPU measured
  0.0e+00 to 9.0e-17), not on the solve succeeding.
* **G3 — no silent fall-through.** §47's defect, which this path is structurally
  primed to repeat: a vector device under `pcnr=True` must never quietly not
  participate. `pcnr_status`/`pcnr_solves`/`pcnr_fallbacks` parity with the CPU.
* **G4 — no regression for what already works.** P19's scalar port stays
  bit-identical where it applies.
* **G5 — cost, measured not assumed.** The CPU carries a 40–60% iteration
  premium for vector PCNR. Measure the traced equivalent; a materially worse
  number is a finding, not a footnote.

### What could refuse this, and should be allowed to

* **`solve_batched`.** Per-device probe blocks whose `m` differs across lanes
  are a shape problem, not a value problem. If per-lane participation cannot be
  made static, batched vector PCNR is refused and says so — the ordinary
  `solve` path is unaffected.
* **Charge-storing participants stay refused**, as they are on both paths
  today (§46: `states` was deliberately NOT relaxed with `vbranches`).
* **`limit_together` (§10.3(b))** is device-level Python state. PCNR REPLACES
  limiting rather than composing with it, so it is expected to be irrelevant
  here — but that is an assumption, and the first measurement should confirm it
  rather than inherit it.

⚠ **The first task is not code.** It is to decide whether `_junction_arrays`
grows a device-shaped sibling or is replaced by one, because gap 1 makes the
existing per-junction metadata the wrong spine for this work — and the
tempting local fix (gap 2's keyword) would build on that wrong spine.

### 49.1 Feasibility, measured the same day — two gaps shrank, one grew, and the fork moved

Probed before writing solver code. Three findings change the plan above.

**The maths transfers, exactly.** Shadowing a participant's `i`/`G` works under
the JAX toolkit, and the resulting excluded assembly matches the CPU's
`augmented_system` **bit for bit** — `|J_mm - J_mm_jax| = 0.000e+00` on a
`MosLevel1Hdl` at a fixed `(x, v_lim)`. And `pcnr_i`/`pcnr_didv` evaluate under
JAX unchanged, returning `(6,)` and `(6, 3)`. So gap 2 really is only a calling
convention, and the assembly half of Stage 3 is mechanical.

⚠ **The exclusion is STATIC here, which the CPU's design obscures.** `pcnr.py`
shadows `i`/`G` per Newton call because it must; the traced loop assembles once
at trace time, so the shadow is applied while the body is traced and never
again. That is cheaper than the CPU's arrangement, not merely equal to it.

**Gap 3 is bigger than this section said.** It claimed `fetlim` and `limvds`
must join the branchless `pnjlim`. Measured, by jitting each law:

    id     TRACEABLE
    pnj    TracerBoolConversionError
    fet    TracerBoolConversionError
    vds    TracerBoolConversionError
    delta  TracerBoolConversionError

**Four of five laws branch on values**, `pnj` included -- the traced path's
existing "branchless pnjlim" is a SECOND implementation of that law, not the
shared one. And the device's own `pcnr_limit` cannot be called at all: it dies
at `hdl.py:6481`, `float(f(x_old_sub, *args))`, coercing a limiter parameter
(SPICE's `von`, read off the last accepted sub-vector) to a Python scalar.

**The redundant-probe rule is NOT on the critical path, which reprieves G1.**
`pcnr.limit_block` resolves over-determined laws by taking probes in order of
decreasing correction and dropping the one that closes a cycle -- a
data-dependent sort, the hardest thing here to make branchless. But it only
runs for a device with a REDUNDANT probe, and measured:

    MosLevel1Hdl      3 probes + 1 redundant   -> needs the rule
    EkvNmosHdl        3 probes, none           -> per-probe loop, exactly
    GummelPoonNpnHdl  3 probes, none           -> per-probe loop, exactly

G1's circuit is the **EKV** differential pair (`test_pcnr_vector.py` ->
`test_limit_identity._diffpair`), which has no redundant probe. So the headline
gate is reachable without the greedy rule, and `MosLevel1Hdl` can be refused
explicitly in a first increment rather than blocking it.

### 49.2 The fork, now worth deciding on evidence

**(a) Reimplement four laws branchlessly inside `jaxtransient`.** Contained,
no CPU blast radius -- and makes a second copy of each law, which the traced
`pnjlim` already shows is a real pattern here rather than a hypothetical one.
Two implementations of a limiting law will diverge; the only question is when.

**(b) Make `_limiting.apply_limit` branchless by toolkit dispatch**, replacing
each `if` with `toolkit.where`, and drop the `float()` at `hdl.py:6481`. Then
ONE law serves both paths, the traced duplicate `pnjlim` can be retired, and
the device's own `pcnr_limit` becomes callable from the trace -- which is the
paper's modularity claim and this repo's stated single-source-of-truth rule.

⚠ **CORRECTION, same day: (b) is right AND cheap, and the sentence that stood
here was wrong.** It said (b) "edits the CPU's limiting laws, which sit under
digest-pinned models and the whole PCNR suite". It does not, because the
pattern already exists: `_limiting.py` carries `_pnjlim` and
`_pnjlim_branchless` SIDE BY SIDE -- "same law, arithmetic selects instead of
Python `if`s, so it runs on numpy scalars/arrays and traced jax arrays alike"
-- with `jaxtransient` importing the twin and the CPU keeping the original. I
wrote the blast-radius claim from the call sites before reading the module.

So (b) is ADDITIVE: three more twins beside their originals. Nothing on the
CPU path moves, and the "bit-identity proof" is a unit test against the
original rather than a suite-wide risk.

**DONE 2026-08-31.** `_fetlim_branchless`, `_limvds_branchless` and
`_deltalim_branchless` added, with `_sel`/`_minimum`/`_maximum` helpers.
Bit-identical to their originals over every branch boundary (~1,900 `fetlim`
comparisons across four `vto`, plus `limvds` and `deltalim` sweeps; zero
mismatches), and all four laws now `jit`, `vmap` AND `grad` -- vmap is what
`solve_batched` needs, grad what PCNR's Jacobian needs, and a law that traced
but differentiated to a non-finite value would pass a jit check and still be
useless. Bite-checked: removing `limvds`'s reverse-mode fold, and a 2.0 -> 2.5
slip in `fetlim`, each fail their own test and only their own.

⚠ Two things worth carrying. `_limvds` RECURSES once for `vold < 0`
(`mos1load.c`'s reverse-mode call); a recursion cannot be selected over, so
the twin unrolls it into a sign -- depth is exactly one, so the unrolling is
complete rather than truncated. And the arithmetic select `a*c + b*(1-c)` is
exact only where BOTH arms are finite (`0*inf = nan`), which is precisely why
`_pnjlim_branchless` floors its `log` arguments before taking them instead of
selecting afterwards; these three are pure `abs`/`min`/`max` on finite values
and have no such arm.

**Recommended: (b), gated on a bit-identity proof taken FIRST** — every law,
over a sweep that includes each branch boundary, `where`-form against
`if`-form, before any of it is wired to a solver. If that proof fails on any
law, that law falls back to (a) and the reason is recorded.

**Increment 1 is then:** tree-probe devices only (EKV, Gummel-Poon), a
device-shaped `_device_arrays` sibling to `_junction_arrays` (the pair view
keeps its own consumer -- gmin targets on the `pcnr=False` path), the traced
augmented assembly proven against the CPU's, and `MosLevel1Hdl` refused by
name until the redundant rule is traced.

### 49.3 The traced limiter, built — and the one blocker left in the compiler

**`pcnr_limit_branchless` is generated beside `pcnr_limit`** for every vector
device, same twin pattern as the laws themselves. It drops the `float()` that
raised (a limiter parameter reading the solution is traced too, so coercing it
cannot work) and dispatches the branchless laws through
`apply_limit_branchless`.

Measured: it reproduces `pcnr_limit` **exactly** on `EkvNmosHdl` and
`GummelPoonNpnHdl`, and Gummel-Poon now **traces** -- `jit` equals eager and
the Jacobian is finite, which is what PCNR's own solve needs.

⚠ **A REDUNDANT-PROBE DEVICE IS REFUSED BY NAME.** `limit_block`'s
over-determined rule is a data-dependent sort, so `MosLevel1Hdl` (4 probes, 3
unknowns) raises with that reason rather than silently receiving a different
law -- which is §47's failure mode, and this path is where it would recur.

⚠ **EKV DOES NOT TRACE YET, and the cause is one line in the compiler.** Its
`fet` law takes a parameter that READS THE SOLUTION (SPICE's `von`,
`_wants_x`), and those chains are compiled against the numpy kernel
unconditionally:

    mods = dict(_KERNEL_NUMPY, _wrapfloor=np.floor)      # hdl.py

so the function cannot accept a tracer. Gummel-Poon's limiter parameters are
constants, which is the whole reason it traces and EKV does not.

**The fix is the treatment `_pcnr_vec_compiled` already gives `pcnr_i` and
`pcnr_didv`**: keep the sympy expression, compile a jax twin lazily under
`_kernel_jax`.

**BUILT 2026-08-31 as `_limit_par_for`, and it works on a COLD BUILD** --
EkvNmosHdl then traces, `jit` equals eager, the Jacobian is finite and the
answer is bit-equal to the CPU limiter.

⚠ **AND IT IS DEFEATED BY THE COMPILE CACHE, which is the finding.** The
ingredients are attached by `_limit_par_fn` at build time; the cache's
rehydration path reconstructs these functions WITHOUT running it, so a warm
process has the generated code and not the ingredients. Measured, same class,
same cache dir:

    first run into an empty cache   ingredients present   traces
    second run (cache hit)          ingredients absent    does not

Bumping `CACHE_FORMAT` 2 -> 3 was necessary -- this is an emission change and
the roadmap's own PCNR_LIFT_AFFINE trap says emission changes must move the
key -- but NOT sufficient: the bump invalidates old entries, and the newly
written ones drop the ingredients just the same. **The payload has to carry
them.**

**DONE 2026-08-31 (`CACHE_FORMAT` 4).** The frozen function record already
carried one attribute -- `_wants_x`, described there as "the one attribute the
compiler hangs on a function that its callers read" -- and there are now two.
`_hdl_limit_par` travels beside it; the COMPILED twin deliberately does not,
being derived, the same reason the `_jax` key is skipped. The ingredients
pickle at 314 and 1429 bytes for EKV's two laws. Measured after: the same
class traces on a cold build AND on a cache hit, with `x_old_sub` genuinely
traced.

⚠ **Invalidating a cache and fixing what it stores are different repairs**,
and the first looks like the second for exactly one run. The 2 -> 3 bump made
the next run work -- because it was then a cold build -- and the run after that
failed again.

⚠ **AND THE OBVIOUS TESTS ONLY BITE WHEN THE CACHE IS WARM.** Asserting the
class has its ingredients passes on any first run, because the class was built
rather than restored, so CI with a fresh cache would sail past this regression.
`test_the_freeze_thaw_round_trip_keeps_the_recompile_ingredients` freezes and
thaws directly, which is the only form of the assertion that holds in both
conditions -- verified by breaking the thaw and watching it fail on a COLD
cache, where the other two pass.

⚠ **The failure is now LOUD.** `_limit_par_for` refuses with the cause and the
remedy (clear the cache, or `PYCIRCUIT_HDL_CACHE=0`) rather than returning a
numpy function that dies several frames deep inside a compiled chain naming
nothing. §47's rule, applied to a path that had every opportunity to repeat it.

⚠ **And the first version of that refusal was too broad, which the first full
suite run caught.** It fired whenever the ingredients were missing, rather than
when a twin was actually NEEDED -- grounding every device whose limiter
parameters are constants, Gummel-Poon among them, and failing two of this
work's own tests. A parameter that does not read the solution never receives a
tracer, so it needs no twin; only `_wants_x` does. Keying the refusal on that
is the fix. Worth recording because a loud refusal is only an improvement on a
silent fall-through while it refuses the right thing.

An EXPIRING test asserts whichever state the process is in -- twin available on
a cold build, refusal on a warm one -- and says to replace it with the positive
form when the payload carries the ingredients.

⚠ **G1's circuit is the EKV pair**, so this blocker is on the critical path
after all -- §49.1 said the redundant-probe rule was off it, and that is still
true, but this takes its place. The order is now: fix the parameter-chain
compilation, then `_device_arrays`, then the traced assembly and loop.

### 49.4 The traced solver: G1 and G2 met on the circuit Stage 2 was built for

`_device_arrays`, `pcnr_vector_blocks` and `pcnr_vector_inner_loop` are in
`jaxtransient.py`. The unit is the DEVICE, and `_junction_arrays` keeps its own
consumer -- the gmin ladders on the `pcnr=False` path, where a conductance
across a FET's `vgs` would be a gate leak.

**The assembly is `pcnr.augmented_system`'s twin, and agrees with it exactly.**
All five blocks -- `g_mna`, `g_lim`, `J_mm`, `J_ml`, `J_lm` -- match to
**0.000e+00** at a RANDOM `(x, v_lim)`, on EKV and on Gummel-Poon. Asserted at
a random point rather than at a solution, where a wrong block can still look
right.

⚠ **The exclusion is static here.** The CPU shadows a participant's `i`/`G` per
Newton call because it must; the traced loop assembles once, so the shadow is
applied while the body is traced and never again -- with §40's affine remainder
shadowed IN rather than zeroed, so a series resistance does not cost the device
PCNR.

**G1/G2, the EKV differential pair** -- two devices, three probes each, the
inter-device clash vector PCNR was commissioned for:

    CPU  pcnr.solve_dc      9 iterations
    JAX  vector PCNR        9 iterations, max|x - x_cpu| = 4.441e-16

and from three starts (0 V, 1 V, 20 V uniform) it reaches the ordinary DC
answer to 0.000e+00. **Gated on agreement, not on convergence** -- §40's
recorded failure is a PCNR solve that converged with an internal node out by a
whole gate voltage while the terminal current matched to seven digits.

⚠ **One inherited criterion was wrong for this path and is fixed here.**
`pcnr.lim_converged` scores the limited residual PER COMPONENT, and says why:
"with vector devices a 40 V `vds` in the same vector would loosen a `vbe`
component forty-fold". The scalar traced loop uses the global
`max|g_lim| < reltol*max(|v|)` form, which is right for a circuit of
one-branch diodes -- where the two coincide -- and wrong the moment a device
owns several unknowns of different scale.

⚠ **NOT YET WIRED.** `_pcnr_setup` still gates on `_junction_arrays`, so
`JAXTransient(pcnr=True)` on a vector device raises as before. The machinery is
tested directly; routing it, and the transient's per-timestep use of it, is the
next step -- with G3 (no silent fall-through, `pcnr_status` parity) and G5
(cost against the CPU's 40-60% iteration premium) still to take.

Two starts cost me a debugging pass each and are worth recording. The source
vector: `augmented_system` takes `u_extra`, and without it the loop solves
`i(x) = 0`, converges in one iteration and returns the start unchanged -- false
convergence that looks like success. And the reference row is PINNED, not
solved, so a start with a voltage on it keeps that voltage for ever; that is
the solver being right and the caller being wrong.

### 49.5 Wired — and the charge, not the current, is what limits the transient

`_pcnr_setup` now routes on the DEVICE view first, so a circuit whose
participants declare `pcnr_probes` reaches the vector solver instead of reading
as empty in the pair view and falling through to the ordinary Newton in
silence. Both `pcnr_meta` consumers were taught the new shape: the per-timestep
Newton and the step controller's Jacobian, the latter through
`pcnr_vector_controller_jacobian` -- `schur_reduce`'s docstring records that
open-coding `circuit.G(x)` there once produced a step count **6.6x too large**.

⚠ **AND THE TRANSIENT CANNOT RUN YET, for a reason that is not PCNR's.**
Measured: **no class declaring `pcnr_probes` has a traceable `q`** -- twelve
classes, none with `eval_q_pure`, and `circuit.q` raises
`TracerArrayConversionError` for every one.

**The asymmetry is the interesting part, and it explains why DC works.** PCNR
SHADOWS a participant's `i`/`G` out of the ordinary assembly and re-stamps its
current at `v_lim` through `pcnr_i`, which does trace -- so the device's own
`i` never runs, and `EkvNmosHdl`, whose `i` does NOT trace, nevertheless solves
to machine precision at DC (§49.4). Its `q` is not shadowed: **charge stays in
the MNA block at the node voltages**, the CPU's documented trade, so a
transient calls it and stops there.

So vector PCNR on the traced backend is, today:

    DC          works, agrees with the CPU to 4.4e-16 on the Stage 2 diff pair
    transient   refused at setup, by reason

⚠ **The refusal had to be added, because wiring made things worse before it
made them better.** Before this work, `JAXTransient(pcnr=True)` on a vector
device raised `NotImplementedError` from the junction path. After routing, it
reached the vector solver and died as a `TracerArrayConversionError` several
frames inside a compiled chain. A setup-time probe (`jax.eval_shape` over the
charge assembly, once) restores a named refusal. **Routing a path before its
failure mode is as good as the one it replaces is a regression, even when the
new path is better.**

**What unblocks the transient is a traceable charge on the participating
devices** -- `eval_q_pure`, the same treatment the batchable classes already
have -- which is DSL work, not PCNR work, and is where G3 and G5 wait.

### 49.6 What actually blocks the traced transient: the CHAINED path has no jax backend

Chasing §49.5's charge blocker to its root, and it is much larger than PCNR.

**25 of 38 library classes are CHAINED**, which is every real compact model --
`EkvNmos/Pmos`, both `GummelPoon`, `DiodeSpice`, `DiodeSpiceThermal`, all four
`MosLevel`, PSP. The other 13 are eager: passives and simple behavioural
blocks. And the chained assembly reads

    if info['chained']:
        ck = funcs['q'].__dict__.get('_hdl_c')      # the C backend hook
        if ck is not None:
            return ck(self, x, epar)
        return np.asarray(funcs['q'](x, *_args_of(self, epar)), dtype=float)

-- a hook for the **C** backend, and for jax a `np.asarray` that cannot take a
tracer. So **the JAX backend cannot run any chained model at all**, with or
without PCNR. That is why `EkvNmosHdl` fails under a plain `JAXTransient` too,
which is worth stating plainly because it means §49.5's refusal is not a PCNR
limitation wearing a PCNR message.

⚠ **§22's disjointness does NOT refuse repairing this, and the distinction is
the point.** §22 measured that a chained model has no `eval_i_pure` and so
loses `solve_batched` -- that finding is about BATCHING, under `vmap`, where
every lane must share one program. TRACING is a different capability: a
`lax.while_loop` needs the assembly to be expressible in jax primitives, not to
be uniform across lanes. Vector PCNR at DC already proves the difference --
`EkvNmosHdl` solves there to 4.4e-16 (§49.4) precisely because PCNR shadows its
untraceable `i` out and re-stamps through `pcnr_i`, which is jax-compiled.

**The work is a `_hdl_jax` hook beside `_hdl_c`**, compiling the chain under
`_kernel_jax` on demand and caching it -- the same pattern applied three times
already today (`_pcnr_vec_compiled` for `pcnr_i`/`pcnr_didv`, `_limit_par_for`
for limiter parameters, and the frozen-record carriage that makes both survive
the compile cache). The machinery exists: `_chain_compile` already accepts
`modules_map=dict(_kernel_jax(jnp), numpy=jnp)`.

**What it would unlock is not only Stage 3.** It is every compact model on the
JAX transient -- which is the backend's whole reason to exist -- and it turns
§49.5's "transient refused" into G3 and G5 being takeable. It is also the
largest single capability gap between the two backends, and it has been
invisible because the batchable 13 are exactly the classes whose tests exercise
the JAX path.

### 49.7 The chained path gets a jax backend — and the traced transient runs

`_chained_jax` beside `_hdl_c`, and `_chained_eval` so the rule that picks
between them is written once for `i`, `G`, `q` and `C` rather than four times.

**It was cheap because every piece already existed.** The compiled function
keeps its generated source in `_src` -- that is how the on-disk cache rebuilds
it -- and `_chain_namespace` is the single place that says what namespace a
chain runs in. Handing that same source a `_kernel_jax` namespace produces a
traceable twin, compiled on demand and cached on the function.

Measured on `EkvNmosHdl`, a chained model the JAX backend could not evaluate at
all an hour earlier:

    circuit.i   TRACES   max|jax - numpy| = 0.000e+00
    circuit.G   TRACES   max|jax - numpy| = 3.388e-21
    circuit.q   TRACES   max|jax - numpy| = 0.000e+00
    circuit.C   TRACES   max|jax - numpy| = 0.000e+00

**And the transient §49.5 had to refuse now runs.** `JAXTransient(pcnr=True)`
on that compact model, against a CPU reference at `reltol` 1e-9:

    JAX reltol 1e-6    177 points   max|v(d) - cpu| = 2.663e-04
    JAX reltol 1e-8    944 points                     2.000e-05
    JAX reltol 1e-10  5077 points                     2.202e-05

falling with tolerance and plateauing at the reference's own residual -- the
two backends converging to the same waveform, ~1e-5 relative on a 1.79 V
signal.

⚠ **Compare converged against converged.** The first reading of this was
4.5e-03 and did NOT improve with tolerance, which looked like a defect in the
new path. It was the CPU reference running at its default tolerance: 85 points
against JAX's 5077. Tightening the reference moved the disagreement by two
orders. The same trap as T4's output node, met again in one day.

⚠ **A shape error caught a scoring error.** The limited-unknown convergence
test compared `g_lim` (length k) against `abstol`, which here is the per-ROW
MNA vector (length n, `iabstol` on node rows and `vabstol` on branch rows).
That is a shape error when the two differ -- which is how it surfaced -- and it
would have been a SILENT mis-scoring had they happened to match. The limited
unknowns are branch voltages and take a scalar floor.

**What this unlocks is bigger than Stage 3**: every chained compact model is
now evaluable on the JAX transient, which is 25 of 38 library classes and the
reason that backend exists. §22's disjointness is untouched -- these classes
still have no `eval_i_pure` and still cannot join a `vmap` group, so
`solve_batched` is unchanged. Tracing was never what that measurement refused.

### 49.8 G5, the cost gate — and an honest note on what the traced path costs

**The PCNR premium is small.** Circuit built once so the trace is paid once,
minimum of three repeats after a warm call:

    easy MOS    CPU  pcnr=False 0.502 s -> pcnr=True 0.532 s    +6.1%
    easy MOS    JAX  pcnr=False 9.171 s -> pcnr=True 9.561 s    +4.3%
    diff pair   JAX  pcnr=False 10.168 s -> pcnr=True 11.499 s  +13.1%

with the accepted step count **unchanged in every case**. Against the CPU's
recorded **40-60% iteration premium** for vector PCNR that is comfortable --
with the caveat that these are WALL and STEP figures and the recorded one is
ITERATIONS, which the traced path does not expose, so it is not the same metric
and should not be quoted as if it were.

⚠ **AND THE TRACED CHAINED PATH IS ~18x SLOWER THAN THE CPU on these
circuits.** 9.2 s against 0.5 s, and it is not compile: successive solves on
one instance run 11.9, 9.5, 9.1, 9.3 s, so ~2.7 s is tracing and ~9.2 s is
execution, every time. A compact model's chain is hundreds of primitives, and
running it inside a `lax.while_loop` pays per-op dispatch on all of them.

**So §49.7 bought a CAPABILITY, not a speed-up, and the distinction should
survive into whatever reads this next.** What it is worth:

* the transient parity Stage 3 needed, which is now measurable rather than
  refused;
* running the same compact model on both backends, which is how a backend
  disagreement becomes findable at all;
* a floor to improve from -- 18x slower is a number to attack, where "cannot
  evaluate" was not.

**What it is NOT worth yet** is a reason to run a compact model on JAX for
speed. `solve_batched` is where that backend pays (22.5x at 512 lanes), and
§22's disjointness still keeps chained models out of it: they have no
`eval_i_pure`, so they cannot join a `vmap` group. A chained model on the JAX
transient today is slower than the CPU and correct, which is the right order to
achieve those two in.

## 2026-09-15 — PSP flicker noise kept a (0.2 mV)² floor at Vds = 0 (peer report) — FIXED

**Reported** by a peer session (pycircuit-only evidence, plus an external-reference agreement not recorded here) and
**reproduced**: IHP `sg13g2_lv_nmos_psp`, W/L 10/1 µm, `swign = 0`, Vg = 0.8 V — flicker part of `CY[d,d]` at 1 Hz
(CY(1 Hz) − CY(1 MHz)) = **2.5686e-23 A²/Hz at Vds = 0**; a fit `K (Vds² + v0²)` gives **v0 = 0.2000 mV**; flicker/Vds² =
3.21e-15 / 7.45e-16 / 6.68e-16 / 6.42e-16 at 0.1 / 0.5 / 1 / 5 mV.

**Mechanism.**  PSP orders the terminals exactly and computes `Sfl = Sfl_prefac · Ids · Gvsatinv · (…)/N1` with
`(…) ∝ Delta_N1 ∝ dps` (`PSP103_module.include:1824-1827`), so `Sfl ∝ Vds²` and is 0 at the origin
(contributed as `flicker_noise(sigVds · MULT_i · Sfl)`, `:1950`).  `compact.py` carries the ordering arithmetically:
the core runs on `vdsa = 2·hypsmooth(vds, 1e-4) − vds = sqrt(Vds² + (2e-4)²)`, so `Ids_core` and `dps` both scale with
`vdsa` and the density with `vdsa²`.  The terminal current is multiplied back by `sgn = vds/vdsa` (0 at the origin); the
density, a power, was not.

**Fix:** `flicker_noise(mult · n_sfl · sgn², ef)` — `(|Vds|/vdsa)²` removes the floor to first order and is `1 − 4e-8/Vds²`
away from it.  After: flicker at Vds = 0 exactly 0; flicker/Vds² 6.4212 / 6.4211 / 6.4205 / 6.4197 / 6.4127e-16 at 0.1 /
0.2 / 0.5 / 1 / 5 mV (p-channel card at negative bias 2.749 / 2.747 / 2.744 / 2.719e-16, flicker(0) = 0); white
density unchanged (3.780778e-23 vs 3.780780e-23).  The 34 PSP noise tests pass.

**Audit of other powers built on `vdsa`.**  The white channel density is EVEN and finite at the origin: its smoothing
offset measured −1.9e-4 relative at Vds = 0 (slope vs `vdsa` × 2e-4) — recorded, not corrected (no clean form).  The
induced-gate pair's coupling already carries `sgn` (an amplitude).  ⚠ **Found absent, not fixed:** PSP's gate-tunnelling
shot noise (`shot_igs`, `shot_igd`, `module:1951-1952`), junction shot noise (`jnoise_s`, `jnoise_d`, `:1953-1954`) and
the edge-transistor flicker/thermal (`Sfledge`, `sqidedge`, `:1955-1956`) have no contribution in the element — the
bulk row of `CY` is exactly 0.  Small for this card's leakage, but structurally missing.

**Pinned:** `test_psp_gap.py::test_the_flicker_density_goes_as_vds_squared_through_the_origin` (n and p: flicker(0) ≤ 1e-12
× flicker(5 mV); flicker/Vds² at 0.1 mV within 1 % of 1 mV).

## 2026-09-15 — PSP gate-tunnelling, junction and avalanche SHOT noise added (found absent during the flicker fix)

PSP contributes `2q|I|` of currents the element already computes (`PSP103_module.include:1886-1906, 1951-1954`), and
none of it was in `compact.py` — `CY`'s bulk row was exactly zero.  Added as white sources carrying `mult` (no `T`):

- gate–source / gate–drain: `2q(gsel|Igcs| + (1 − gsel)|Igcd| + |Igsov|)` and mirror — the CHANNEL tunnelling part
  follows the ordering (as PSP's forward/reverse `shot_igs`), the OVERLAP part does not;
- bulk–source / bulk–drain: each junction's `2q|Ijun|` plus the avalanche shot `2q (mavl + 1) |Iimpact|` on the HIGH
  terminal (`(1 ± sgn)/2`).
- The edge transistor (`Sfledge`, `sqidedge`) stays out: SWEDGE = 0 on this PDK's cards and the element has no edge path.

On the IHP n-channel 10/1 µm card the gate current is 4.2e-9 A at (Vd 0.05, Vg 1.2) — shot 1.3e-27 A²/Hz, ~1e-7 of the gate
row; the card sets no junction saturation currents; avalanche 4.7e-13 A at Vds 1.2 V.  Small, but PSP-faithful, and it
dominates the gate row at weak inversion (1e-30 against induced-gate densities of 1e-35..1e-41).

**Pinned:** `test_the_gate_and_avalanche_shot_noise_are_2q_times_their_currents` — identities against the element's OWN
terminal currents (no reference number): at Vds = 0, `−CY[g,s] = 2q|I_s|`, `−CY[g,d] = 2q|I_d|` (rel 1e-9); at Vds = 1.2
with gate leakage off, `−CY[b,d] = 2q (mavl + 1)|I_b|`, `mavl = |I_b|/(|I_d| − |I_b|)` (rel 1e-6); n and p.  Fails on HEAD
(the entries were 0).

**Tests whose MEASUREMENT had to change (not their tolerances):**
- the two all-zero tests (`test_it_is_off_without_a_card`, `test_zero_resistance_collapses_the_branch_and_its_noise`)
  and the induced-gate `swign`/`fnt` zero tests and the `f²` slope test now switch leakage off (`NO_LEAKAGE`) — shot noise
  is not switched by FNT/SWIGN, exactly as in PSP;
- `test_the_gate_density_matches_psp` compares PSP's `sig` (induced only) with `sig(swign on) − sig(swign off)`, which
  removes the uncorrelated gate shot exactly with the card unchanged;
- `_split` (PSP `sid`/`sfl` comparisons) reads the drain diagonal with the shot taken out, `CY[d,d] + CY[gi,d] + CY[b,d]`
  (each shot source is the only one between its nodes).  ⚠ A first attempt used `−CY[d,s]` and was WRONG — the correlated
  induced-gate source also sits on noi–s, and its cross term lands on [d,s]: `sid` read +51 % (long) / ×14.5 (short).
All 423 PSP gap / C-backend / compile-budget tests pass.

## 2026-10-02 — 50. The generated chain's repeated subexpressions, computed once (bit-identical)

Andreas (2026-10-02, after the speed plan): "fused compiled hdl evaluation 1
Create a deep explanation for me. 2 do a deep investigation of the topic to
understand the topic and limitations. 3 make the plan".  Plan approved
(`plans/twinkly-bouncing-finch.md`): F0 measure, F1a CSE, F1b the C
backend's libm, F2 fusion of i/q/G/C, F3 evaluation sessions, F4 'auto'.

**The waste, counted in PSP's cache entry.**  `G` repeats `i`'s 719 value
definitions and adds 2148 derivative statements, each re-printing the local
partials `diff(expr, parent)` inline -- once per unknown (6): exp 303, sqrt
1747, where 2313 calls against `i`'s 51 / 210 / 187.  `C` repeats `q`'s 388
and adds 1362 (1222 text-identical to `G`'s).  Of `G`'s ~24.7k calls ~3.75k
are distinct.  `_chain_compile` does no CSE (the eager path's lambdify does,
per function).

**Why not in sympy.**  Substituting a symbol for a subexpression changes
sympy's canonical ordering and flattening of `Add`/`Mul` -- `Mul(2,a,b,dj)`
against `Mul(P,dj)` -- so association and summation order, hence bits,
move; and every compact-model result with them (as `fold_card`, sec. 30, and
`_autohold`, sec. 36, did).  ON THE GENERATED PYTHON'S SYNTAX TREE it is exact:
a subtree evaluates as a unit in the same order inlined or hoisted, the
function is SSA, the primitives are pure, and the printer emits nothing lazy
(the design review scanned the 300 newest cache entries: no IfExp/
BoolOp/lambda/comprehension; `numpy.where` evaluates both arms).  Constants keyed by type
and repr (`0` / `0.0` / `-0.0`), never hoisted alone.

**F0 (measured, no source change).**  A prototype on the cached `_src` of
PSP and the 22 chained library classes: every one of 88 functions
byte-identical over its sweep (PSP: the C backend's 1005-point spike sweep),
raise-mode behaviour unchanged.  PSP `G` 13.9 -> 3.2 ms (4.31x), `C` 7.2 ->
1.6 ms (4.54x), `i` 1.32x, `q` 1.31x; library `G` 1.22-5.17x, median 3.0x.
Fusion by statement name on top (F0.3): all four at one state 23.4 ms ->
6.3 ms CSE'd separately -> 3.4 ms fused (0.54x; library 0.46-0.70x),
byte-identical.  The C backend with the `-fno-builtin-*` libm declared
`__attribute__((const))` (F0.4): `G` 50.0 -> 22.5 us, `C` 26.4 -> 15.6 us,
`i`/`q` identical; `G`/`C` differ from the plain kernels in NaN SIGN BITS at
16 of 1005 points, every one already NaN on every path (the backend's own
contract counts the same 16 against numpy); `pow` calls in `G` 5503 ->
2843.  End to end with the functions swapped in-process (F0.5): the PSP
stage's PSS gear 15.5 -> 4.5 s, radau 24.4 -> 6.9 s, waveforms identical.

**F1a (built).**  `_hdl_cse.py`: an iterative walk and hash-consing,
DAG reference counts, hoisting each repeated non-trivial subtree before
its first use (reusing a variable whose whole value it is), text spliced
from the source spans (no `ast.unparse`: a ~300-deep statement exhausts its
recursion), and `_verify` -- every variable and the return value expanded
into the inputs must be the same expression.  Run after the compile cache
(`BehaviouralMeta`, before the C backend attaches), stored in `cse/` beside
it (keyed on the module's own source, the Python version and the input
text), `PYCIRCUIT_HDL_CSE=0` off.  `_src` stays the reference (explain, the
JAX twin, the C kernels' print), the twin in `_src_cse`, the reference
function in `_hdl_ref`; `compiled_jacobian_size` reads `_hdl_codelen`, the
reference size (PSP's optimised `G` is a third of it, and 18 KB on C is what
the 'auto' threshold of 10 KB is calibrated against).

**Measured traps.**  The benchmark circuit's PSP is a COLLAPSE VARIANT
(`PspMosLongChannel_collapse10`) with its own functions -- a swap on the
base class reached nothing and read as "no gain" until the instance's class
was used; a test reading the base class's functions for an element built as
a variant fails with an IndexError (fewer unknowns).

**Gates.**  G75 (the first build): every transient, PSS and PAC call
identical to G74; five failures, `test_hdl_params`' adopter records, which
digest `explain()` -- the CSE note had been a new line; it rides on the
`backend:` line now, which those records leave out (an evaluation fact, not
the compile).  G76: 3565 passed, every call identical.  One warnings record
moved: `TestARealPcnrSolve` logged two "overflow encountered in exp" before
and none after -- the twin warns from other locations (`<hdl-chain-cse>`),
and Python reports a location once; under `-W always` the same overflows
occur (3 emissions without CSE, 2 with: the duplicates removed).  The PSP
benchmark (`benchmarks/speed_analyses.py --psp`): gear 16.2 -> 5.4 s, radau
24.4 -> 6.8 s.  The cost: each `hdl.py` edit recompiles every model once
(~70 s; the compile cache keys on its source), the pass itself ~2 s for
PSP on a store miss.

## 2026-10-02 — 51. Fused compiled evaluation, F2 and F3: one fused pass per state, and the stage paths' evaluations recorded

**F2 (`_hdl_cse.fused`)**: several of a chained class's i/q/G/C at one
state computed by one function -- their statements merged by name (the same
name is the same text in every function `_chain_compile` printed for a
class: asserted, `Refused` otherwise), each function's order kept, the
union optimised and stored like the F1a twins.  Built lazily for the exact
set a session names (any subset of two or more), so a site never computes a
chain it does not read.  Every fused set byte-identical to its separate
functions (PSP's 1005-point sweep; MosLevel3, thermal Gummel-Poon, EKV,
SPICE diode in the tests).

**F3 (`pycircuit/circuit/_evalhint.py`, `_hdl_cse.take`)**: an evaluation
session names the methods a site will request at one state; inside it a
chained element's first request runs the fused pass and keeps the other
outputs for the calls that follow (keyed on the session, the state's bytes
and the `_args_of` list -- a new list after `update()` or a temperature
change; each output handed out once; at most 8 states per element).  The
analyses still make every call -- the circuit passes, the tests' Jacobian
counters and PCNR's instance shadows see the same calls.  Bypassed for a
JAX or symbolic toolkit, `params_tree`, DC pins, a method shadowed on the
instance, a C-bound function, a class whose reference `G` is under 5 KB of
bytecode; a fused pass that raises falls back to the separate call.
`PYCIRCUIT_HDL_FUSE=0` turns sessions off.  The sites: the multistep step's
`_residual_and_jacobian` ({q, i, G}, + C unless cached -- `_C_lookup` split
out of `_C_at_state`), the chord's `residual_only` ({q, i}), the converged
point as ONE session object the branch screen's `C` read and
`jacobian_only`'s `q`, `G` both re-enter; the dense Radau stage ({q, i, G,
C}), the transform's iterations ({q, i}) and start ({C, G}), the DIRK
stage's `func_i`, `_finish_stage_step` (what the memo lacks).

**The stage paths' device memo, completed (the F0 side finding).**  The
PSP radau PSS evaluated `C` about 9 times a step at 3 distinct states and
`G` 5 times at 3: the transform's start reads `x_n` = the previous step's
final stage, which `_finish_stage_step` had evaluated and not recorded; the
shooting's stage step re-read `C` at the stages the branch screen had just
read, and `G` at the final stage.  Memo records are now partial and merged
(`_memo_put`); `_finish_stage_step`, the transform's start and
`_C_at_state` record where the memo rolls (`_memo_ok`: the coupled stage
steps -- the transform rolls its own generation now -- no stateful limiter,
no bypass); the transform's start, `_coupled_stage_context`, PSS `_C_at`
and `_G_at` read it.  Recorded values are the same function's at the same
state: bit for bit.

The PSP stage (`_cs_amp`, 40 points; waveforms byte-identical to the
reference before F1a at every row):

| | before F1a | F1a (CSE) | + memo | + fusion |
|---|---|---|---|---|
| gear PSS | 15.5 s | 4.5 s | -- | 2.9 s |
| radau PSS + factored period | 24.4 s | 6.9 s | 4.6 s | 4.2 s |

Tests (`test_hdl_cse.py`, 12 more): the session API, fusion and its
refusal, fused sets byte-identical, one fused pass per state per session
with distinct writable outputs and a recomputation after a parameter
change, every bypass, MosLevel3 inverter transients (gear, radau, trbdf2)
byte-identical with sessions and the memo off (gear: the same circuit `G`
call count).

Gate G77: 3576 passed; every transient, PSS and PAC call identical to G76
(the new tests' own calls aside); no warning differs.

## 2026-10-02 — 52. Fused compiled evaluation, F1b and F4: the C kernels' libm calls declared const; 'auto' re-calibrated

Andreas: "Do both F1b and F4".

**F1b.**  `_KERNEL_C` declares the sixteen libm functions `CFLAGS` keeps
real calls (`-fno-builtin-<f>`) `__attribute__((const))`: gcc computes an
identical call once instead of at every repeat (a chained `G` prints each
local partial once per unknown) -- still glibc's call, so still its value.
PSP `G` 50 -> 22 us, `C` 26 -> 16 us, `pow` calls 5503 -> 2843; `i`/`q`
byte-identical; `G`/`C` differ from the plain kernels only in NaN sign
bits, at the 16 spike points where the plain kernels already differ from
numpy so (the backend's contract).  `source_key` hashes the prelude, so
every `.so` rebuilds on its own.  The PSP stage's PSS on C: gear 0.40 ->
0.32 s, radau 1.06 -> 0.85 s (full Newton).

**F4.**  The 2026-10-01 'auto' table re-measured with the twins, the fused
passes and the const libm running (`doc/transient_history.md`, "'auto'
re-calibrated"): the gains shrank (numpy chord 5-32 %, transform up to 68
%), no crossover moved -- `AUTO_JACOBIAN_CODE`, `C_KERNEL_SHARE` and the
reference-size proxy stand (the reference size does not move with
`PYCIRCUIT_HDL_CSE`, which keeps that switch bit-identical end to end).
The option descriptions' percentages corrected.

## 2026-10-02 — 53. Speed round 2, stage A: the C backend's defects closed before it becomes the default

Andreas: "Make plan do implement 1, 2, 3 and 4, if 4 is worth it" -- the C
backend as the default (every chained model, built at a class's first
instance), exact scalar fast paths on the numpy path, a compiled-code store,
fused C kernels if they pay (`plans/twinkly-bouncing-finch.md`).  Stage A,
nothing changed on the numpy default:

* DEFECT: numpy 2.5 on x86 returns the SECOND operand of `maximum`/`minimum`
  on a +-0 tie (it follows `maxsd`), the prelude's `_npmax`/`_npmin` the
  first -- a zero's sign and `1/x`'s infinity flipped against numpy.  The rule
  is measured at import (`_numpy_tie_rule`) and printed into the prelude.
* DEFECT: `_chained_eval` tried the C kernel before the jax twin -- a C-bound
  class failed under the JAX backend on its tracer.
* DEFECT: `CKernel` raised on an array temperature (or a parameter that is not
  one number); it now steps aside for the numpy function, which broadcasts.
* DEFECT (F2's): `_hdl_cache.freeze` of a live `info` after a session raised
  `Uncacheable` on the fused `_f` -- `_RUNTIME_KEYS` are not frozen.
* cffi declared (and checked: `'numpy (cffi not installed)'`), `info['_c_bound']`
  first in `take`, an exponent that is a `numpy.where` value refused by the C
  printer, a `flock` build lock per key.

## 2026-10-02 — 54. Speed round 2, stage B: compiled code from a content-addressed store

Each process compiled every generated function it loaded (0.78 s of a 1.80 s
library import; 0.27 s of PSP's 0.28 s class creation).
`_hdl_cache.code_for(text, filename)` serves the marshalled code object from
`cache_dir()/code/` (keyed on the interpreter's bytecode identity, the
filename and the text), compiling and storing on a miss; used by the chain
compile, the cache thaw, the CSE twins, the fused functions and the JAX
twins.  Import 1.80 -> 1.06 s, PSP's class 0.28 -> 0.02 s; no cache-format
bump; every result bit-identical (gate G80).

## 2026-10-02 — 55. Speed round 2, stage C: exact scalar fast paths for the numpy twins

`_hdl_fast` + `_hdl_cse._fast_rewrite`: the twins call fast paths for
`numpy.where`, `_step`, `maxc`/`minc`, `numpy.maximum`/`minimum` and the
comparisons, each answering exactly what numpy answers on the scalars it
decides and handing every other case (ties, NaN, arrays, int64, sympy) to
numpy.  The two measured places a numpy scalar differs from numpy's 0-d
array are closed: a `where` whose value reaches `**` stays `numpy.where`
(numpy squares a 0-d array, a scalar goes through `pow`), and a fast `where`
is only ever finite (which NaN two NaNs give follows the operand order).
PSP `G` 3.16 -> 1.90 ms, library median `G` 2.5x, the PSP stage PSS gear
2.75 -> 1.71 s, radau 3.30 -> 1.97 s; bit-identical (gate G81).
`PYCIRCUIT_HDL_FAST=0` turns it off.

## 2026-10-02 — 56. Speed round 2, stage D: the C backend is the default

`'auto'`: every chained model runs C where a compiler (or its stored
objects), cffi and the compile cache serve it, resolved at a class's first
instance (`_hdl_cbackend.ensure` from `Behavioural.__init__`); numpy
quietly elsewhere, the reason in `_hdl_backend_status`.  Cold builds: the
library 0.04-0.66 s per class, PSP 16.7 s, once per machine.  On the suite,
every moved result is either the 'auto' Newton options reading a C kernel at
a hundredth (verified: with the share at 1 all but the `tanh` models are
identical to numpy) or the `tanh` ulp (comparator, op-amp, divider); C
raises no numpy warnings.  `PYCIRCUIT_HDL_BACKEND=numpy` restores the old
default.  A vacuous F1a test (a swap on a collapse variant's base) fixed.

## 2026-10-02 — 57. Speed round 2, stage E: fused C kernels measured and refused

A multi-output C kernel per session subset (PSP's `_csrc` merged by
statement name) computes all four at a state in 28 us against 66 us
separately, byte-identical -- but end to end the PSP stage's PSS on C moves
-10 % under gear and +1 % under radau (whose sessions are mostly `(i, q)`,
where the session bookkeeping eats the saving), against a gate of >= 10 %
under both; and four more PSP-sized objects would double to quadruple the
cold-build CPU.  Not built; the numbers in the log.  What is left in a C run
is the per-step Python (the kernels are 28-36 % of it).

## 2026-10-02 — 58. Speed round 3, stage L1: the limiter's parameter values once per call

Measured for the plan: on a 20-MosLevel1 chain the hdl `$limit` was 62 % of a
gear step -- each probe's parameter functions (chain-compiled lambdas over
`maxc`/`minc`) evaluated twice per probe per Newton iteration, ten
evaluations where one is needed.  The generated `limit()` now builds the
values once per call: the `x0` readers at the previous iterate, the pure
ones kept on the instance against the argument list (`_hdl_lp`, by the
function object, dropped by `update()`).  20 MosLevel1: 0.382 -> 0.202 s,
20 Gummel-Poon 0.447 -> 0.218 s, PSP untouched; bit-identical (gate G84).
`LIMIT_PAR_CACHE = False` is the old behaviour.

## 2026-10-02 — 59. Speed round 3, stage W: the leaner kernel call

`CKernel.__call__` cost 4.3 us around a mid-sized model's 0.3 us of C: two
`ffi.cast` of `arr.ctypes.data` and `np.ndim(T)`.  Now `ffi.from_buffer` on
the state and the output with a pointer type resolved once, a type test for
the usual float `T`, the dtype singleton; the fresh output, the fallbacks
and the ValueError kept; and a C-bound class's element methods skip the
session call.  MosLevel1 4.3 -> 1.5 us a call, PSP `G` 23.2 -> 20.3; gear
transients: 20 MosLevel1 0.202 -> 0.150 s, 20 Gummel-Poon 0.218 -> 0.155 s,
20 PSP 0.372 -> 0.309 s; bit-identical (gate G85).

## 2026-10-02 — 60. Speed round 3, stages O' and L2: one pass per state refused; the limiter chains as twins

O' (the plan's item 1, one evaluation per element per state inside the
stamp plan) was gated on a ceiling measurement: one hand-rolled loop against
the four circuit passes saves 42 us per state on a 20-MosLevel1 chain (18 %
of the passes) -- ~5 % of a gear step with its ~1.1 Newton iterations --
and nothing on 20 PSP, where alternating four large kernels per element
thrashes the instruction cache (+3 to -10 %); the gate was 8 %.  Not built.
After stage W the per-element pass overhead is ~1.5 us; what remains in a
step is the kernels and the per-step solver machinery.

L2: `_hdl_cse.optimise_limit_pars` wraps each chain-compiled limiter
parameter over the CSE + fast-path twin of its inner (`von` 10.6 -> 5.6 us);
20 MosLevel1 0.150 -> 0.144 s, bit-identical (gate G86).

Round 3 in sum (gear transients, 100 fixed steps, 20 devices): MosLevel1
0.382 -> 0.144 s (2.65x), Gummel-Poon 0.447 -> 0.154 s (2.9x), PSP 0.383 ->
0.310 s (1.24x); the PSP stage PSS 0.301 -> 0.266 s.

## 2026-10-02 — 61. Speed round 3, after L2: `-O3` refused, the limiter body streamlined (L3)

`-O3` for the C kernels (fast-math and contraction still off): byte-identical
over 307 points and no faster (PSP `G` 19.1 -> 19.3 us) -- straight-line
scalar code gives the vectoriser nothing.  Not adopted.

L3: the generated `limit()` takes its branch voltages from Python lists
and reuses the ranking's limit where the inputs are the same bits (a
signed-zero guard): `limit()` 20.5 -> 18.3 us a MosLevel1 call, the chain
0.144 -> 0.138 s; a differential test against the old body over 2000
random cases.  Round 3 closes here: 20 MosLevel1 0.382 -> 0.138 s (2.8x),
20 Gummel-Poon 0.447 -> 0.149 s (3.0x), 20 PSP 0.383 -> 0.310 s (1.24x),
the PSP stage PSS 0.301 -> 0.266 s.  What is left: the kernels' own C work
on PSP-class models, and the per-step solver machinery on small circuits.

## 2026-10-03 — 62. Speed round 5: the hdl limiter in C (`_hdl_climit`; LC1, the kernel proven offline)

Andreas, after speed round 4: "Plan for the hdl limiter C kernel".  The
round's in-run tree had left `cir.limit` at 33 % of a 20-MosLevel1 gear
step -- 18.8 us of generated Python a call after round 3's trimming (the
body 17.6: lists, two rankings, sorts, the write-back; the parameter
chain `von` 6.2; the laws 2.6) against 1.7 us for a C kernel call.  So
the whole closure is printed, from the spec and the parameter functions'
kept sympy statements (`_hdl_limit_par`), through the chain printer's
own symbol map (`hdl._c_symmap`, factored out of `_render_c` with its
text unchanged); the laws and `device_writeback` are transliterated in
Python's own forms in a prelude of the kernel's own (`_LIMIT_C`, never
`_KERNEL_C`'s: that would rebuild every object).

Three things the design review changed.  Python's `sorted` and `min` on
keys holding a NaN order by timsort's and the set's accidents, which no
other sort reproduces (measured: a stable insertion sort disagrees in 29
% of random NaN-keyed lists) -- so the kernel returns an `int` and
DECLINES such a call, leaving the closure to answer it with its own
order and its own warnings (`_hdl_cbackend._dlopen` and `load_kernel`
take a `cdef` for that; `kernel_for` is the wrapper the chain functions
keep).  MosLevel1/3's `von` evaluates both arms of a `where` and warns
on `sqrt` of a negative at a forward-biased body where C is silent, so
the gate's warning counters are to be censused before the kernel is
bound (LC2).  Binding is per collapse variant, with `nx` the `i`
kernel's, never the terminal count.

LC1 proves the kernel offline: `test_hdl_climit.py` sweeps the five laws
through a probe kernel against `apply_limit` (20 000 random points and
every pair of special values), every library class with `$limit` against
the closure over 2000 random cases, a chained sequential group, the
inputs the closure accepts, the temperatures, the declines, and shows
the sweep can fail.  Standalone the kernel is 1.7 us a call where the
closure was 3.8 (the SPICE diode) to 20 (the thermal BJT).

### LC2 (2026-10-03): bound and taken

`_hdl_cbackend._resolve` renders the limiter beside the chain functions
and binds it after them; `detach` unbinds; the generated `limit()` tries
the kernel first (not on a JAX or symbolic toolkit) and keeps the closure
for what the kernel hands back as None; `explain()`'s `backend:` line
says `(limit: c)`.  The tests that count the closure's parameter
functions pin numpy, and the reference transliteration checks both
paths.  Measured against LC1 (d8f74c42): 20 MosLevel1 1256.7 -> 801.5 us
a step (-36.2 %), 20 Gummel-Poon -12.7 %, the PSP cases unchanged (no
`$limit`); `el.limit` 2.05-2.31 us from 3.8-17.2; bytes and statistics
identical.  Since round 3's close the 20-MosLevel1 chain is -42 % a
step, since round 1's start 4.8x.  What a mid-sized model's step is now:
the kernels' call floor (~1.7 us a call, 1.2 of it wrapper) and the
assembly's numpy floor (round 4's refusal); nothing of the limiter.

## 2026-10-03 — 63. Speed round 6: one C call per class per pass (`_hdl_batch`)

Andreas, after round 5 closed: "Continue".  What a mid-sized model's
pass was made of after rounds 2 and 5 was the CALL around each kernel:
on the 20-MosLevel1 chain `cir.G` cost 56 us, of which the twenty
kernels' C work was 6.7 -- the other 50 us twenty wrapped calls (`el.G`
-> the generated method -> `_chained_eval` -> `CKernel.__call__`, 2.0
us each) through the stamp plan's per-element loop (0.8-1.1 us each);
in the run the four passes were 560 of a 971 us step.  The driver
prototype (a loop in C over the class's elements calling the kernel's
own function pointer on `x[NM]` and each element's own pack) gave the
pass in 7.5 us, bit for bit today's, before anything was built.

`_hdl_batch.py`: `split` at the plan's build gathers the eligible
elements of each C-bound chained class into a `Batch` -- no DC pins
(the three classes with them keep the per-element path's `G_dc`/`i_dc`
choice), the method the metaclass generated by CODE identity against
`BehaviouralMeta.__init__`'s nested definitions (a spy, a `wraps` spy
and a Python override each carry their own code object; `getattr(cls,
m) is K.__dict__[m]` would pass a class-patched spy), the circuit's
toolkit, a bound kernel of the element's size and the slots' shape, two
or more of them; batches group by class, so a plain subclass has its
own.  `Batch.run` per pass: `ENABLED`, the binding, the kernel's
identity (re-cast when it moved, refused when its layout did), the
method's code, the temperature one number (`CKernel`'s rule); per
element an instance shadow (PCNR installs `el.i`/`el.G` on the
participating devices DURING a solve) or a failed pack hands THAT
element back and the rest of its class stays one call.  The driver
writes the temperature into each element's OWN pack as `CKernel.
__call__` does (one consistency surface), and the batch MIRRORS the
tuple the element holds at the pass, never repacking an existing one:
`state_restore` (the branch check's speculative Newton) puts an old
`__dict__` back with no epoch move, and both paths must read that pack.
The driver is one small C source of its own exporting the chain
functions' entry name with its own signature, built and loaded once per
process through `load_kernel` under its own key (the limiter's
precedent).

Measured (parent dc2214da, 5 interleaved rounds, bytes and statistics
the same on every case): 20 MosLevel1 802.3 -> 445.9 us a step (-44.4
%), 20 Gummel-Poon 1116.2 -> 667.9 (-40.2 %), 20 PSP -16.5 %, the PSP
stage and its PSS unchanged; the gate >= 25 %.  Since round 1's start
the 20-MosLevel1 chain is 3.82 -> 0.45 ms a step (8.6x).  ⚠ The in-run
tree wrapped every element's methods on the instance -- an instance
shadow, which is exactly what sends an element to the per-element path
-- so its first run showed no batch; it wraps the hand-written elements
only now and times `Batch.run`.

### Commit 2 (2026-10-03): the zero source term skipped

After the batches the largest circuit piece of a 20-MosLevel1 step was
`cir.u` (105 of 544 us with timers): forty-four element `u` calls, every
one returning `[0, 0, ...]` -- 22 of 22 chained library classes have no
source term (their compiled `u`, `u_dc` and `dudt` end in a `return` of
literal zeros) -- each allocating a vector to add exact +0.0 to a bin
that starts at +0.0.  `_add_element_subvectors` skips such an element as
it skips `Circuit.u`'s default (`_hdl_batch.zero_source`): decided once
per class by an `ast` read of the compiled function's source, kept as
`info['_u_zero']` (a runtime key of the compile cache, never pickled),
checked per call against the method's code identity (a subclass that
overrides `u` shares its parent's info and must be called), on the
default skip's conditions plus no `dtype` and not the 'ac' analysis (the
complex `uac` path, where a class with `ac_stim` terms answers).
`PYCIRCUIT_HDL_ZERO_U=0` keeps the calls.  Measured (parent 40a976d6, 5
interleaved rounds, bytes and statistics the same): 20 MosLevel1 446.9
-> 412.7 us a step (-7.6 %), 20 Gummel-Poon -7.1 %, 20 PSP -2.9 %, the
PSP stage and its PSS unchanged; the gate >= 5 %.

### Commit 3 (2026-10-03): the limiter walk; the round closed

After commit 2 the largest circuit piece left in a 20-MosLevel1 step was
`cir.limit` (78 of 544 us with timers): twenty wrapped `el.limit` calls
at 2.3 us around twenty limit kernels of 0.34, plus the loop's Python.
The loop is sequential -- the next element reads the earlier write-back
(every chain device shares gnd, a cascode shares a node), `x0` often IS
`x`, duplicate rows (source and bulk on gnd) take the last value -- so
the passes' batch is not its shape.  `_hdl_climit.limit_walk` walks
every run of C-kernel elements in ONE C call on the live `x`: the loop's
order, its gathers before the call and its scatter after, the kernel's
own pointer on the element's own pack with the temperature written as
`CLimitKernel.__call__` writes it; the driver returns at the first
element Python must answer -- no kernel (a hand-written limiter, a nested
SubCircuit, a patched or shadowed `limit`, a failed pack) or a call its
kernel declines (a NaN among the keys) -- and the caller runs the loop's
own statement for it before the walk resumes behind it.  Where the walk
does not serve (a state that is not one writeable float64 vector of the
circuit's size, a strided or complex `x0`, a temperature that is not one
number, a JAX or symbolic toolkit, the circuit-level resolution) nothing
is touched and the loop runs.  A rebound kernel is taken per call; the
walk's cache is keyed as the stamp plan's.  `PYCIRCUIT_HDL_LIMIT_WALK=0`
keeps the loop.

⚠ Two faults of the first cuts, both found by measurement or the suite:
the walk asked the cache "is anything here kernel-capable?" LAST, after
the state, temperature and driver checks, so a circuit whose limiters
are all hand-written Python (PSP's own, `compact.py`) paid 2-3 us a call
for nothing -- the three PSP cases read +0.5..+1.6 % twice, and a count
showed `cir.limit` called 2.15 times a step there where I had claimed
it never was (every hdl class has a `limit`; PSP its own); the cache is
asked first now.  And a capable class with no kernel bound at both the
build and the call dereferenced None (`test_device_limiter`'s PCNR grids,
whose classes have a `$limit` spec and no C source).

Measured (parent 1aff4af1, 5 interleaved rounds, bytes and statistics
the same on every case): 20 MosLevel1 409.7 -> 352.1 us a step (-14.1
%), 20 Gummel-Poon 626.3 -> 510.5 (-18.5 %), 20 PSP -0.3 %, the PSP
stage +0.2 %, its PSS +0.1 %; the gate >= 5 %.  The round closed: since
dc2214da (round 5's close) the 20-MosLevel1 chain is 802 -> 352 us a
step (-56 %), 20 Gummel-Poon 1116 -> 511 (-54 %), 20 PSP 2853 -> 2329
(-18 %); since round 1's start 20-MosLevel1 is 3.82 -> 0.35 ms a step
(10.8x).  What a 20-MosLevel1 step is now (the in-run tree, timers on,
392 us from 1048 before the round): `cir.C/G/q/i` 39/30/31/31 us, `cir.u` 34, `cir.limit`
19; the rest the Newton's own bookkeeping and linear algebra, the
companion assembly and the predictor.

## 2026-10-03 — 64. Speed round 7, stage 0: the Newton iterate in C -- the ceiling measured, the build refused

The plan (the iterate out of Python: three C primitives behind a slim
loop, bit-identical through numpy's own `scipy_dgesv_64_` and
`scipy_cblas_dgemv64_`) put a two-day ceiling first, with a gate of
1.4x on the 20-MosLevel1 step and a refusal below 1.3x.  The evaluate
half was built in scratch, bit for bit the Python path on every trial,
and measured warm in the run: 20-MosLevel1 357 -> 311 us a step (-12.9
%), 20-Gummel-Poon -6.9 %, the PSP cases -2..-4 %.  1.15x: refused.  The
premise -- the evaluation at 143 of a 392 us step -- came from the
in-run tree's INCLUSIVE timers (fifteen nested inside the residual, ~1
us each); the warm standalone cost is 66 of 355, and the rest of the
iterate is 36 us.  The full build would have reached ~1.2x.

What the profile found: the Python around the existing C -- the four
per-pass `Batch.run` check loops (~17 % of a step; one check serves all
four passes), `zero_source`'s per-element `is_generated` (~8 %), the
assembly wrappers, the walk's per-call loop, the predictor.  An
estimated 15-20 % for a few days, no new C.  Lesson, recorded in the log
and the memory: a tree of inclusive timers ranks pieces, it does not
size them; size a piece standalone and warm before a plan rests on it.

## 2026-10-03 — 65. Speed round 7, commit 1: the evaluate core (`_tran_core`)

After stage 0's refusal (sec. 64) the Python trims the profile had
suggested were sized standalone at ~3 %, and Andreas chose to build the
prototype for real.  `_tran_core.py`: one C call evaluates the passes a
site asks for (the batches through the kernels' own pointers, the
constant vector groups through numpy's own `cblas_dgemv`, the templates,
the bincounts from +0.0 in array order), then the companion in the
integrator's own operation order -- gear2, euler, trap, theta -- and the
residual; the three evaluation sites of the LMM step take it
(`_residual_and_jacobian`; `jacobian_only` at the converged point; the
chord's `residual_only` with the conductance from the held `C`), each
leaving the state the Python path leaves, Python still deciding the
active integrator and asking `_C_lookup` first.  The tables sit in one
C struct filled once per stamp plan; a lone C-bound element is a batch
of one; the readiness check runs once per call over the unique
elements.  Everything the core cannot serve declines before touching
state.  Pinned against its reference (`test_pinned_pairs`: the pair
fired once during the build, when a lint fix moved the twin -- the
practice working).

Measured (parent 82819d61, 5 interleaved rounds, bytes and statistics
the same on every case): 20 MosLevel1 357.8 -> 277.3 us a step (-22.5
%), 20 Gummel-Poon -18.9 %, the PSP stage -7.6 %, its PSS -8.0 %, 20 PSP
-3.2 %; the gate >= 12 % and >= 6 %.  Since round 5's close the
20-MosLevel1 chain is 2.9x, since round 1's start 13.8x.  What is left
of a step: the Newton's own bookkeeping and the step machinery, in
pieces of a few microseconds each; the C work is ~15 us of ~277.

## 2026-10-06 — 66. Speed round 12, stage 1: the C library against numpy, by CPU

Section 19's contract -- the C kernels return the numpy path's bytes --
rests on numpy's float64 transcendentals being the C library's.  They are
on this box but for `tanh` (exception 1 above); numpy's AVX-512 loops
bring their own `exp`, `log`, ... (a peer machine failed 50 kernel tests).
`_hdl_cbackend.libm_check()` now measures it once a process (10.5 ms: each
function of the prelude, in the forms the numpy path calls it, against
the C library through ctypes), and a class whose C calls a function that
differs KEEPS C and names it -- `cls._hdl_backend_ulp`, the backend line of
`explain`, one `simwarnings.PlatformWarning` a process beyond tanh
(Andreas's decision, kernels and limiters alike).  The exceptions above
become: (1) every LISTED function, tanh on every CPU measured -- the class
bitwise its source run with the C library's functions; (2) unchanged.  The
tests apply the rule wherever C meets numpy, and a pytest plugin
(`scripts/plant_libm.py`) simulates such a CPU on this one.  Record:
`doc/pss_log_260902.md`, "speed round 12, stage 1".


## 2026-10-07 — 67. Speed round 12, stage 4: the pack's parameter values at dict speed

The pack (`_params_of`: every parameter of an element, each time its
packed array is rebuilt) read each value by `getattr`, ~2.3 k
instructions through `ParameterDict.__getattr__` -- the failed lookup and
the hook, which no leaner body avoids.  It now reads the values dict
exactly where `getattr` would reach it (an exact `ParameterDict`, the
class's own `__getattr__` and the default lookup, the name a parameter, on
neither the instance nor the class), `getattr` otherwise; switch
`PYCIRCUIT_PARAM_DIRECT`, a pinned pair with the whole `ParameterDict`
class.  `BSource` reads its function once a call.  Measured: the element
cases -1.0..-3.3 % in instructions, -0.9..-1.4 % in time (12/12).
Record: `doc/pss_log_260902.md`, "speed round 12, stage 4".


## 2026-10-07 — 68. Speed round 12, stage 4 reverted

Section 67's direct read and `BSource`'s single read are reverted:
the stage landed under its plan's refusal line on mos1 (-1.25 % against
-1.5 %), and Andreas kept the line.  `_params_of` reads every value by
`getattr` again.  Record: `doc/pss_log_260902.md`, "speed round 12, stage
4 REFUSED".


## 2026-10-07 — 69. `BSourceHdl`: the behavioural source in the HDL

Andreas, on the open question of a C path for circuits with Python
elements (the stamp plan's passes are 36 % of the van der Pol PSS): "Can
we build the BSource from the HDL?"  `elements.BSource` holds arbitrary
Python callables, so the evaluate core refuses any circuit that has one
(`_tran_core`: a non-constant element that is not a C-bound chained HDL
class) and every pass runs element by element in Python and numpy.

`elements_hdl.BSourceHdl(inp, inn, outp, outn, i_func=f, q_func=g)` takes
the same terminals and callables, TRACES each callable once with a
symbolic control voltage, and compiles ``I(outp,outn) <+ i(u) +
ddt(q(u))`` with ``u = var(V(inp,inn))`` -- one class for each pair of
expressions, made at the first instance and kept, its compile cached by
the expressions (the closure's `srepr`), the instance returned is of it.
The `var` is not decoration: a class without an intermediate is not
CHAINED, and only chained classes run in the core -- prototyped without
it, the van der Pol PSS was 6 % SLOWER than with `BSource`; with it the
core takes the circuit whole and the C Newton every step: 0.470 -> 0.161 s
(2.9x).  What traces: Python arithmetic, powers, sympy's functions
(`sympy.Piecewise` for a branch); numpy's functions, an `if` on the
voltage, a free symbol or a complex value are refused, naming
`elements.BSource` for them.

NOT `BSource`'s numbers, by construction: `BSource` takes its derivative
as a central difference (`Toolkit.derivative`, step 1e-6), this one is
exact; and the expression is sympy's (``u**3/3.0`` becomes ``0.333... *
u**3``).  On the van der Pol PSS the waveforms agree to 2.7e-14 and the
monodromy matrices to 1.1e-10 -- the difference's own error.  So it is an
element to choose, not a replacement: `BSource` is unchanged, `vdp_pss`
stays on it, and `vdp_pss_hdl` is the same PSS on the HDL source
(`benchmarks/step_machinery.py`; a case that needs what a tree may lack
is left out of a comparison whose tree lacks it -- `NEEDS`).  Pinned pair:
the behavioural source and its HDL form.  Record:
`doc/pss_log_260902.md`, "BSourceHdl".

## 2026-10-08 — 70. The MOS gate charge: charge-conserving, following Meyer (the SPICE benchmark plan's stage 6)

`elements_hdl._mos_gate_charge` and four variants,
`MosLevel{1,1Pmos,3,3Pmos}GateChargeHdl` (`_mos1_analog`/`_mos3_analog`
with `gate_charge=True`); `MosLevel{1,3}[Pmos]Hdl` are unchanged.  The
importer builds the variants where SPICE adds Meyer's capacitances
(level 3 always; level 1 when TOX is given and not zero).  Andreas
(2026-10-07): a charge model whose capacitances follow Meyer's, not
SPICE's integration of Meyer.

The charge: inversion `(2/3)*Cox*(a^2 + ab + b^2)/(a + b)` over the gate
drives `a`, `b` at source and drain (its gradient IS Meyer's triode and
saturation capacitances at `vdsat = vgs - von`), and below threshold two
functions of `max(us, ud)` whose slopes are Meyer's depletion/accumulation
and near-threshold capacitances.  Checked against a transcription of
ngspice's `DEVqmeyer` (`test_mos_gate_charge.py`): the gate's total
capacitance Meyer's to 1e-15 of Cox in accumulation, depletion,
transition and inversion; the source/drain/bulk split exact in inversion;
inside Meyer's 25 mV sliver his total at most 1.47x; conservation (the
4x4 columns and rows sum to zero, both polarities, body effect on); round
a closed terminal loop Meyer's integration leaves 0.06 Cox, this one the
quadrature's 1e-4.  Named deviations (the docstring): below threshold the
charge follows `vgs`; the sliver; body-effect transcapacitances; level 3's
triode (its `vdsat` is not `vgs - von`: at most 0.25 Cox with `vmax`);
level 3's threshold without `eta`'s DIBL.

Found while building:

* **DIBL made a negative capacitance.**  A charge following a `von` that
  falls with `vds` has the drain's own capacitance negative
  (`-(4/15)*Cox*etal`), and these decks give the drain no other
  (`ad = pd = 0`, `cgdo = 0`): an inverter chain on mike2's cards grew a
  mode at 3.9e13/s and the step collapsed to 6e-18 s.  The charge now sees
  `von` at `vds = 0`.
* **The hard `max(us, ud)`** moved the near-threshold capacitance from
  source to drain in one step as `vds` changed sign; gm3 failed at the
  time-step floor.  Rounded over `|vds| < 25 mV` (C1) it runs.
* **The channel split decides the timing.**  An inverter chain on mike2's
  cards against ngspice's Meyer, 50 % crossings: Ward-Dutton +20 ps a
  stage, every channel charge to the source -30 ps, the split that
  follows Meyer's terminal capacitances (`Qd = -(Q(a,b) - Q(a,0) +
  Q(b,0))/2`: nothing in saturation, half at `vds = 0`) +10 ps.  Andreas:
  both, selectable -- `qpart` (0 Ward-Dutton, the default; 1 following
  Meyer; a blend between), the importer's `mos_qpart`, the suite's
  `--qpart`.  Exact reciprocity throughout needs a `d ln d` term (an
  infinite slope at `vds = 0`).  The ablations `gamma = 0`, `vmax = 0`,
  `eta = 0` (on both simulators) moved it little.

The gold decks (Xyce's metric, 1 or below passes; stage 4 had no gate
charge; ngspice on our own written deck beside them):

| deck | stage 4 | qpart 0 | qpart 1 | ngspice |
|---|---|---|---|---|
| gm1 | 39.7 | 5.81 | 4.29 | 0.511 |
| gm2 | 20.3 | 2.37 | 1.71 | 0.612 |
| gm3 | 39.6 | 2.16 | 1.34 | 0.651 |
| gm17 | 3.2e9 | 96.9 | fails (Newton stalls at 5.75 ns) | 12.9 |
| mike2 | 32.7 | 1.92 | 1.86 | 0.681 |
| todd3 | 4.36 | 1.06 | **0.73** | 3.18 |
| 4049osc (no TOX: unchanged) | 0.385 | 0.385 | 0.385 | 0.406 |

The remaining gap is the model, not the numerics (gm2/mike2 at reltol
1e-5 score the same), and the metric swings with small changes (the
rounded `max` moved gm1 2.74 -> 5.81, mike2 2.07 -> 1.92).  The gate
charge adds real dynamics: steps roughly double, at ngspice's own count on
the inverter chain (430 against its 419 points); the force-accepts (20-40
on these decks, at the sources' corners) and the rejections are open.
`qpart = 1` doubles the force-accepts and fails gm17.  The plan's lift of
the junctions' 0.9 grading clamp was dropped: no deck has `mj` near it.

## 2026-10-08 — 71. Gummel-Poon: the forward diffusion charge is `TF*IF/qb` (a defect, fixed)

SPICE's forward diffusion charge is `TF_eff*IF/qb` (ngspice `bjtload.c`:
`cbe = cbe*(1 + argtf)/qb`, then `qbe = tf*cbe`); `_gummel_poon` had
`tff*ifwd`, without the base charge `qb`.  Found by the SPICE benchmark
plan's stage 7: the rca wideband amplifier (CircuitSim90; VAF = 50, so
`qb` ~ 0.89 at its bias) ran 8 % of its output swing away from ngspice on
the same written deck, with both simulators converged (ours to 0.05 ns
steps, ngspice's `.tran` capped at 0.05 ns) and the operating points equal
to 3 uV -- a difference in the stored charge.  With `/qb` the two agree to
1e-4 of the swing.

The test reference had been transcribed from the model, not from SPICE
(`_gp_reference`'s `qbe = tff*ifw`), so it could not see it.  Now from
`bjtload.c`; the RK4 transient's ODE carries the full charge Jacobian
(`Qbe` sees `vbc` through `qb`); the closed-form capacitance writes
`d(TF*IF/qb)/dVbe` out.  Reverting `/qb` fails all three.  Every class
built on `_gummel_poon` moves: `i` and `G` bit-identical at 25 biases a
class, `q` and `C` moved; the recorded point digests re-recorded (old ones
in the notes beside them).  ngspice divides only where `vbe > 0`; here
always, so the charge is continuous (`TF*IS*(1/qb - 1)` apart below zero
bias).

rca's gold metric: 2.76 -> 2.31 with the fix at the deck's own stepping,
0.800 (passing) when our step is also capped at SPICE's default `TMAX =
min(TSTEP, TSTOP/50)` -- the cap `write_ngspice` gives ngspice.  Whether
the importer should apply that default to our own runs is open: the gold
is Xyce's, which has no such cap, and the cap costs 4049osc 5x the steps
(0.385 -> 0.389).

## 2026-10-08 — 72. The bipolar substrate junction: `GummelPoon{Npn,Pnp}4Hdl` (the SPICE benchmark plan's stage 7)

Terminals `(c, b, e, s)`, built on `_gp_core` with a substrate (the
3-terminal classes and their recorded digests unchanged): SPICE's
substrate charge (`bjtload.c`: depletion below zero bias, the capacitance
linearly extended above, no current), `cjs`/`vjs`/`mjs` with the
junctions' temperature path, attached by `subs` -- 1 vertical, on the
internal collector (Xyce's, whose gold the suite scores against), -1
lateral, on the internal base (ngspice's default for a p-n-p), polarity
`type*subs` as ngspice's.  The importer builds them where a card gives
CJS: on the line's substrate node, or ground on a 3-node line; vertical in
Xyce's dialect, a p-n-p lateral in ngspice's; `write_ngspice` writes SUBS
explicitly.  `test_bjt_substrate.py`: the charge is ngspice's on both
sides of zero bias, on the collector (or the base) and nowhere else; a
zero junction is the 3-terminal device bit for bit.

The stage's census (the decks' cards): no deck needs `ptf`; `irb` only the
HB gilbert cell (its card gives one -- the first census saw only the MCNC
decks), which now stops there; the substrate junction latch and opampal.

| deck | before | now | ngspice on our deck |
|---|---|---|---|
| latch | refused (CJS) | **0.107** | 0.307 |
| opampal | refused (CJS) | DC operating point fails | finds it |
| rca | 2.76 | 2.31 (0.800 with SPICE's TMAX; section 71) | 0.588 |
| schmitecl | 0.108 | 0.126 | 1.54 |

OPEN: opampal's DC.  ngspice finds the operating point of our written
deck; ours fails (pseudo-transient exhausted, a residual 1e5 over at the
-35 V supply), and fails too when started FROM ngspice's point (internal
nodes at their terminals): residual 1e4 over at `xb:15`, updates of 2.8 V
at `xd:qt.ci`.  The cards are unusual (a p-n-p with BR = 25990, IKR =
12.8 uA, MJE = MJC = 0.99); not yet localised.

## 2026-10-08 — 73. opampal's DC: source stepping falls back to adaptive steps (section 72's open item, closed)

Localised in four measurements, each ruling a cause in or out:

1. **The models agree.**  Each of the deck's two cards alone, against
   ngspice on the same written deck, 25 biases (forward, saturated,
   reverse, `vce` -1..30 V): equal to ~1e-5 wherever current flows.  The
   only differences are at 1e-14..1e-10 A, and they are SPICE's GMIN (1e-12
   S across each junction, kept permanently by ngspice; ours has none).
   (A first version of this comparison read ngspice's output wrongly and
   its NaNs compared as "agree": `NaN > tol` is False.  The flag is now
   `not (rel <= tol)`.)
2. **GMIN is not the cause.**  1 TOhm across every junction: the same
   failure, to the same residual.
3. **Our equations hold ngspice's point.**  Started from ngspice's
   operating point (its node voltages, each internal node at its
   terminal), our DC converges and agrees with it to 3.5e-5 V.  (A first
   attempt left the internal nodes at 0 V and read as a singular Jacobian
   at ngspice's point -- an artefact of the starting vector.)
4. **The continuation was the failure.**  `SourceSteppingNewton` stepped
   the sources 0, 0.01, 0.1, 1 and nothing between: on +-35 V supplies and
   516 unknowns the last jump failed and so did every later rung.  Steps
   that grow on success and shrink on failure, as SPICE steps its sources,
   reach the operating point ngspice finds (to 2 uV) in 16 steps.

The fix: the fixed ladder unchanged (a circuit it solves is solved exactly
as before), the adaptive steps only where it fails, from the last factor
it solved (`SourceSteppingNewton._adaptive`; tested with a stub solver
that converges only for small jumps).  Its reach is wider than one deck: a
circuit the fixed ladder gave up on may now converge -- gate 6-2's test,
which caps Newton at 2 iterations to force a failure and check its
message, converged (steps that small need 2), and now switches the
fallback off (`ADAPT_MIN`) to keep testing the message.  opampal now passes its gold:
**0.013** (ngspice on our deck 0.019).

Cost, measured and left: its DC takes ~200 s (the transient 101 s in
all), ngspice's whole run 0.2 s.  176 s of it is 156 FAILED leaf Newton
solves (100 iterations at ~11 ms on 519 unknowns), because each source
step's inner solve is the whole gmin ladder, which runs to exhaustion
before a step is cut.  A faster DC (SPICE's order: gmin stepping, then
source stepping with plain Newton steps) is its own piece of work.

## 2026-10-08 — 74. Gummel-Poon `irb`: SPICE's current-crowding base resistance (stage 7, finished)

`irb > 0`: ngspice's law on the base current `IB` (`bjtload.c`): `x =
max(IB/(IRB*area), 1e-9)`, `z = (sqrt(1 + 14.59025 x) - 1)/(2.4317
sqrt(x))`, `rbb = RBM + 3(RB - RBM)(tan z - z)/(z tan^2 z)`, its series
`1 - 4z^2/15 - 4z^4/105 - 8z^6/1575` below `z^2 = 1e-4` (where `tan z - z`
cancels; ~1e-11 apart at the seam).  `irb = 0` keeps the `qb` law: the
point digests did not move; the explain digests did (a parameter and new
intermediates), and one `select` adopter digest moved by nine signed zeros
in one Jacobian entry (`G[12]`, +0.0 <-> -0.0, values equal), each
recorded beside its digest.  The importer reads IRB.

Against ngspice on a crowding card (RB 100, RBM 10, IRB 100 uA; the base
current crosses IRB in the sweep): 1e-5 at every bias; committed as
ngspice's numbers in `test_bjt_substrate.py`, with `irb = 0` shown to be
another device there (>1 %).  ⚠ The first reference was ngspice's `.dc`
sweep at its default tolerances and read as a 9e-4 disagreement: a sweep
continues each point from the last at `reltol = 1e-3`, so it carries
~1e-3 of its own error.  Taken at `reltol = 1e-9` it agrees with ours.

The stage's verification, where it lives: the substrate charge against a
transcription of `bjtload.c` (section 72); the transient against ngspice
on rca (1e-4 of its swing, section 71); `irb` above.

**PTF (excess phase) stays refused.**  Only the HB gilbert cell gives it
(14.6 degrees, with TF = 0.53 ns): its delay is ~PTF[rad] * TF = 0.135 ns,
about 1e-6 of the deck's 10 kHz period -- below anything its HB gold can
resolve.  Whether to import such a card with a note (excess phase left
out) or to build the plan's Bessel network is the stage-8 decision, when
gilbert's HB runs.

## 2026-10-08 — 75. PTF (excess phase), and SPICE's default TMAX for imported decks

**PTF** (Andreas: support it).  Weil's approximation as ngspice's
`bjtload.c`: the forward transport current `IF/qb` through `td^2 y'' + 3
td y' + 3 y = 3 IF/qb`, `td = PTF[rad]*TF`, and `y` in its place in the
collector current.  ngspice integrates it inside the device by backward
Euler on its own step; here it is two internal states (`td y' = w`, `td
w' = 3(IF/qb - y - w)`, scaled to volts by `_XP_R`) the simulator
integrates with the circuit, on nodes referenced to the EMITTER TERMINAL
and collapsed away at `PTF <= 0` -- referenced to the internal emitter the
collapse could not chain (it left zero-volt branches: 7 unknowns where
today's device has 3).  The states are raw node voltages, so a `PTF > 0`
device is declined by vector PCNR; identity probes would have kept it, at
two PCNR unknowns on EVERY device, and broke the thermal class's affine
remainder -- measured, refused.  At `PTF = 0` the variants an instance
runs are today's; the recorded digests moved by last-bit reassociation
only (worst 3.3e-15, measured on every recorded point and sweep), noted
beside them.  A PCNR test helper read the BASE class, which no instance
runs and which now carries the states: it reads `type(instance)`.

Against ngspice (a common-emitter stage at 100 MHz, TF = 1 ns, PTF = 60):
the collector within 0.12 mV of ngspice's (committed), where PTF = 0 is
~80 mV away; the operating point does not see PTF.  The importer reads
PTF; the HB gilbert cell imports (no deck in the census is refused now
but for MOS level 2).

**TMAX** (Andreas: use SPICE's default).  `Imported.tmax()`: the `.tran`'s
TMAX or `min(TSTEP, (TSTOP - TSTART)/50)`, for our runs and `write_ngspice`
alike (which had used TSTOP alone).

The gold suite with both (Xyce's metric; ngspice on our deck beside):
4049osc 0.389 (0.406), toronto 0.462 (0.242), rca **0.800** (0.588),
schmitecl 0.126 (1.537), latch 0.101 (0.307), opampal 0.013 (0.019), gm1
5.81 (0.511), gm2 2.38 (0.612), gm3 2.16 (0.651), gm17 96.9 (12.9), mike2
1.88 (0.681), rich3 29.5 (95.3), todd3 1.06 (3.18); no gold: arom 4.0, gm19
72, jge 50 against ngspice.

⚠ **slowlatch FAILS (minstep at 10.0 ns) -- a stage-6 regression**, with or
without the cap: without the gate charge it runs (138, as in stage 4).
No growing mode (unlike the DIBL case); Newton stalls at node 99 (residual
55x over with updates already at tolerance): a kink.  The gate charge
still has capacitance JUMPS -- the gate's total at threshold for vds ~ 0
(2/3 Cox -> Cox, Meyer's sliver) and the source/drain rows at threshold
for any vds (the near-threshold image split 50/50, Ward-Dutton's 40/60
above) -- and at a small step `C/h` turns a jump into a Newton stall.  A
gate charge with continuous capacitances is the fix; next, before stage 8.

## 2026-10-09 — 76. slowlatch's step collapse: the error test on source currents (`lte_algebraic`)

**Section 75's explanation was wrong.**  The continuous-capacitance gate
charge it called the fix was built (a 25 mV rounding band; every
capacitance continuous) and slowlatch failed at the same 10.0167 ns;
without it, once the step no longer collapses, the run is the same (666
steps, metric 136.2 either way).  Not built.

What collapses the step, each link measured:

* the error estimate is dominated by `br4` -- VSEU2's current (-2.6 mA,
  tolerance 7.7 uA, estimate 5-7 uA), mapped there through `J^-1` from
  node 12's charge (the gate it drives);
* it does not fall with the step: ~5e-6 from h = 1e-12 down to 1e-16 s,
  so it is not truncation error; not Newton noise (vabstol 1e-6, 1e-8,
  1e-10: the same estimate); no growing mode;
* node 12's charge history at the collapse gives the source current
  2.00, 2.24, 2.47 mA over three steps of ~0.9 fs: real femtosecond
  content, from the internal source/drain nodes behind the 40 ohm sheet
  resistances coupled through the gate charge -- `rsh = 0` and it runs.

ngspice runs it because it does not judge source currents at all (a
source declares no truncation; each device judges its own charges in
charge units).  VACASK judges every unknown, but on predictor minus
corrector of its values; with the current moving ~10 % a femtosecond it
would most likely hold the step down too -- not measured.

**`Transient.lte_algebraic`**: 'judge' (the default, today's) or 'skip',
which leaves out of the error test every branch current no charge depends
on -- `C`'s column and row zero at a generic state, read once a run: a
voltage source's current is left out, an inductor's (its flux) judged,
node voltages always.  The Python chain and the C error test carry the
same mask (`StepLTEInputs.lte_mask`; the C sets the left-out entries to
+0.0 before the maximum): slowlatch and gm2 under 'skip' give the same
waveform bytes on either path; the pinned pair re-recorded after that.

**Measured as the default and refused there.**  With 'skip' everywhere:
18 tests failed, six of them on accuracy -- an HDL Laplace element (3 %),
a VCO's phase, the driven fold on radau (0.36 % vs 0.05 %), a BJT's
turn-off no longer lengthening with TF, the state-event period, a diode's
forced steps -- and the two ladder timing cases changed bytes (a third
fewer error tests).  A current no charge depends on can still be the
quantity a circuit is about.  slowlatch, gm2 and gm3 gained (gm2 20 -> 0
forced steps, 753 -> 362 steps; gm3 29 -> 0, 579 -> 177) with their gold
scores unchanged (2.34, 2.09).  Andreas: an opt-in, not even the
importer's default -- set only where a case needs it: the suite runner's
`OVERRIDES` gives slowlatch 'skip'.

slowlatch now completes (665 steps, 6 forced) and scores 136 against its
gold -- stage 4's 138 (before the gate charge), where ngspice on our deck
scored 0.494: a separate discrepancy, open.

## 2026-10-09 — 77. Stage 8: PSS on the SPICE decks (record in `doc/pss_log_260902.md`)

The 4049 oscillator's limit cycle (autonomous, gear on a one-period
`lte_grid`): period 196.98892 us against the gold's settled 196.9905,
every node within 3 mV.  Xyce's two HB decks by driven shooting from the
operating point (from zeros the common-emitter stage diverges):
common-emitter xyce_verify 0.094 / 0.808 (the 0.808 is the gold's own
50-harmonic truncation), harmonics within 2e-4 of the fundamental; the
Gilbert cell 0.016 / 0.009, h1/h3/h5 to 4 digits.  Three gold tests, and
`spice_ce_pss` in step_machinery.  Open: radau 50-450x slower than gear
on these decks, unsized.

## 2026-10-09 — 78. MosLevel3's channel noise in weak inversion (reported from another session's Si2302CDS model)

The channel noise was the Klaassen-Prins integral
`(vgtn^2 - vgtn vdsn + vdsn^2/3) / max(vgtn - vdsn/2, 1e-9)`, `vgtn` from
`von`, `vdsn = min(|vds|, vdsat)`.  With `nfs`, `vdsat` comes from the
clamped `vgsx = max(vgs, von)` and is measured from `vth`: below threshold
it stays ~`xn Vt` while `vgtn` is 0, so the term was `(vdsn^2/3)/1e-9` --
1e4..1e5 S, 1e8..1e10 x `2 q Id` (a Si2302CDS-like card, VTO 1.17 V, W
0.273 m: the floor 1e-13 A^2/Hz flat below 1.3 V).  Above threshold the
same mismatch made the denominator small near `von` (S/(8/3 kT gm) 2.39
at Vgs 1.35) and the vds -> 0 value fall short of `4 k T g_ds` (0.28 at
1.3 V).

Now: `vdsn` at most `vgtn` (the integral falls to 0 with `vgtn`, as
`2 vgtn/3`; unchanged where `vdsat <= vgtn`), plus the weak-inversion
diffusion current's shot noise `2 q |Id| coth(|vds|/2 Vt)` -- `2 q Id`
saturated, `4 k T Id/vds` as `vds -> 0` -- weighted by
`exp(-vgtn/(xn Vt))` (`_mos3_channel` returns `xn` now).  Measured on
that card: weak inversion saturated `S = 2 q Id` exactly; `vds` = 1 mV,
`S/(4 k T g_ds)` 1.01 in weak inversion (0.83-0.96 near threshold, the
strong part's `von`/`vth` offset; 1.04 at 3 V); saturated, `S/((8/3) k T
gm)` monotone from 2.71 (= 3 xn/4, the shot value) through 1.59, 1.44,
1.18, 1.05 to 1.00 at 1.8 V -- where the old form spiked and then dipped
to 0.84.  Test: `test_noise_devices.py`, the level-3 anchors.

EKV's noise is charge-based (`4kT beta n Ut (qs+qd)`: `2 q I` in weak
inversion) and level 1 has no subthreshold current: neither has the
defect.  There is no level 2 yet (stage 9).

**The reporter's second item, flicker against ngspice (5e4 x on an
AD8606 input pair), is a CONVENTION, found in ngspice-47's source**: its
levels 1, 2 and 3 default to `nlev = 2` (`mos{1,2,3}set.c`), the form
`KF gm^2 / (f^AF W Leff Cox)`, not `KF Id^AF / (f Cox Leff^2)` (its
`nlev = 0`, ours).  The ratio, ours over ngspice's, is `W Vov^2 / (4 L Id)`
for a square law: 5.0e4 for that pair (W/L 1600, KP 10u, Id 0.6 mA).  Not
changed: a choice for Andreas (an `nlev` parameter, and the importer's
default).

## 2026-10-09 — 79. `nlev`: ngspice's MOSFET flicker forms, and `ddx` through `var` intermediates (Andreas: "Add nlev option")

`MosLevel1Hdl`/`MosLevel3Hdl` (and their PMOS and gate-charge classes)
take `nlev`, ngspice's flicker selector, `cox` per area, `Leff = l - 2 ld`:

| nlev | PSD | who |
|---|---|---|
| 0 (default) | `kf |Id|^af / (f cox Leff^2)` | SPICE2; ngspice's compatibility mode; pycircuit until now |
| 1 | `kf |Id|^af / (f cox W Leff)` | the gate area -- as a commercial simulator normalises (relayed, pss_log 2026-09: not measured here) |
| 2, 3 | `kf gm^2 / (f^af cox W Leff)` | ngspice-47's default for levels 1-3 (`mos{1,2,3}set.c`) |

`gm` is the device's own `dId/dVgs`, drain, source and bulk held -- the
limited `vgs` the current reads, so SPICE's `gm` -- taken by `ddx`, which
now sees through `var` intermediates as Verilog-A's does through local
variables: forward accumulation over the let-chain, each intermediate
between the probe and the expression given its derivative as an
intermediate of its own, linear in the chain; the probe may be an
intermediate (`ddx(ids, vgs)`).  Measured: every form against its formula
to 1e-6 on both models, two frequencies at `af = 1.2`; `ddx` through a
chain against the analytic derivative and the compiled Jacobian to 1e-14.
The importer reads `NLEV` (it refused the card before); without it the
device keeps 0.  The default stays 0: changing it moves every flicker
result (an AD8606 pair: 5e4 between 0 and 2).

## 2026-10-09 — 80. Stage 9: MOS level 2 (`MosLevel2Hdl`), and `hdl.quartic_root`

**The model.**  `MosLevel2Hdl` / `MosLevel2PmosHdl` and their gate-charge
classes, transcribed from ngspice-47's `mos2load.c`/`mos2temp.c` into
`_mos2_channel` (one arm, level 3's `Piecewise` on `vds` for the two
modes) and `_mos2_analog` (level 3's junctions, overlaps, `rd`/`rs`,
limiting group, noise, `nlev` flicker and gate charge).  Everything in the
channel: Grove-Frohman bulk charge with the 3/2-power body term, `xj`
short-channel gamma, `delta`, `nfs` subthreshold, `ucrit`/`uexp`
mobility, `vmax` velocity saturation, `lambda` given or computed from
`nsub` (and with `vmax`, `neff`), the punch-through limit; `vto` from
`nsub` with `tpg`/`nss` where not given; `phi` 0.6 without `nsub`.
Deliberate departures, each in the docstring: no `lvds <= 1e-10` band
(below); `xlamda*lvds` formed without dividing by `lvds`; a NaN quartic
candidate invalid; ngspice's ueff/xv read the model's (untempered) `uo`,
transcribed as such.

**`hdl.quartic_root(a1, b1, c1, d1)`**: SPICE's velocity-saturation root
(Ferrari through the resolvent cubic, the smallest positive root with
`|P| <= 1e-6`, 0.0 for none), its derivative implicit (`-(dP/da_k)/P'`);
scalar in `math` (the libm the C prelude calls), per element on numpy
arrays (bit for bit the scalar's), masked on jax; the C prelude twin.
SPICE's closed form gives its root to ~1e-8 relative -- a finite
difference needs a 1e-4 step to see the derivative.

**Verified**, `tests/test_mos_level2.py`: an independent transcription
(`tests/_mos2_reference.py`) against ngspice-47's own currents (13 cards,
every feature, 1560 points recorded in `tests/data/mos2_ngspice_iv.json`:
worst 3.8e-9 of a card's peak, ngspice's solve tolerance); the model
against it on pycircuit's constants (worst 4e-11, rounding); gm/gds/gmb
from the Jacobian against its differences (worst 4e-7, the differences'
own); numpy against C on every value of i/G/q/C (bit for bit, 12 cards x
32760 values); PMOS the NMOS mirrored; the importer (LEVEL=2 to the
gate-charge class, VTO from NSUB with NSS/TPG).  Compile ~4 s.

**g1310's DC point found the band.**  SPICE zeroes the current for
`lvds <= 1e-10` and keeps `gds` for its Jacobian; transcribed, the exact
derivative of that zero is zero, so a switch held on at vds ~ 0 (x6:mn,
vds = 6 nV) had no conductance inside the band and Newton bounced across
it -- from ngspice's own operating point too.  Without the band the
expressions reach 0 with vds and their slope there is SPICE's `gds`.

**The gold decks** (Xyce's metric; ngspice on our deck beside):
ab_integ **0.004** (0.003), e1480 **0.008** (ngspice aborts: 'Timestep
too small' at its first point -- the suite now records that), gm6 4.14
(3.88), mux8 1.05 (0.34), g1310 4.52 (1.49), nand 2.40 (0.55), mosrect
7.34 (0.94), ring 72.6 (24.4); hussamp refused (an `.ic` on a node a
source holds -- the importer's, not level 2's).  The DC channel is
ngspice's (above, and on mosrect's card at |vds| to 11 V, forward bulk to
0.6 V: 1e-7), so the transient gaps are in the charges.  mosrect is
traced: the bridge devices' bulk is another bridge node, and our gate
charge's named deviation (1) -- below threshold Meyer's gate-BULK
capacitance appears against the SOURCE -- loads node 4 where Meyer loads
node 2 (without our intrinsic charge the gap falls 0.19 -> 0.07 V while
ngspice keeps Meyer's).  nand's stacked NMOS have the same topology
(bulk at ground, sources inside the stack); not separated there.
Open, for Andreas: a bulk-referred accumulation/depletion charge.
Not done yet from the plan: voter25, the ring's PSS, mem_plus,
`spice_ring` in step_machinery.

## 2026-10-10 — 81. SPICE's junction GMIN on the Gummel-Poon (Andreas: "add GMIN", then "default 0")

`gmin` (default **0**, the card's physics; SPICE's 1e-12 S for parity)
across the base-emitter, base-collector and substrate junctions, where
`bjtload.c` adds it: in the base currents (`cben += gmin*vbe`, `cbcn +=
gmin*vbc`, so irb's crowding law reads it as ngspice's does) and on the
substrate junction (ngspice's `gdsub` without ISS).  Default 0 because with
1e-12 a dark high-beta phototransistor's collector current is GMIN's
(5 V x 1e-12 x beta 1200 ~ 6 nA against ~0.1 pA of junction leakage).

WHY IT WAS WANTED: another project's optocoupler (nested), its
phototransistor's base fed by nothing but a dark photocurrent, failed DC
with the base "in no equation" -- bisected to 4ec03518, whose last-bit
reassociation changed the Newton path.  From the old answer the new code
converges (base column 8e-14 S there): a PATH problem.  `gmin=1e-12` on
that device is SPICE's remedy.  PCNR does not reach it (off by default in
DC; it limits forward overshoot).

ngspice's REVERSE junction limiting (`DEVpnjlim`'s `vnew < 0` branch,
Gillespie, ngspice NEWS "Ng-spice-rework-13", 2000; absent from SPICE3f5
and the PCNR paper's listing) was ported and MEASURED, and reverted: it
cost the high-voltage starts -- a level-1 cascode from 40 V under PCNR 77
iterations (budget 30), a PCNR grid point that plain Newton solves no
longer converging, the grouped MOS limiter's advantage gone -- while the
mechanism it addresses was not the one at work.  The reporter's
iterate at the failure (2026-10-10): base -2.3e83 V, the internal
collector +2.4e74 V, and a node held by a 3.3 V SOURCE at -2.2e58 V -- the
linear solve itself had returned garbage, the near-singular Jacobian (the
base held by ~1e-18 S beside ~1e-3 S) blowing a step up in EVERY unknown;
the empty column (every exponential exactly 0 at 1e83 V) is the symptom.
Junction limiting cannot cure that; GMIN, by restoring the conditioning,
does (`gmin=1e-12` on that device, now the reporter's fix).  OPEN, a
solver item: DC Newton should reject a step that violates the linear
source rows that grossly, as a failed linear solve, rather than carry it
to the structural check.

With `gmin=0` every recorded point digest of the three classes is the
previous record bit for bit (measured); the explain digests moved with
the parameter.  The physics tests pin `gmin=0` explicitly; the reference
transcription carries GMIN as `bjtload.c` does.  Test:
`test_spice_gmin_holds_every_junction`.

## 2026-10-10 — 82. The gate charge's depletion is bulk-referred (Andreas: "fix it")

`_mos_gate_charge`'s named deviation (1) is gone: below threshold the gate
charge followed `m = vgs - von(vbs)`, so Meyer's gate-BULK capacitance sat
against the SOURCE -- invisible where the bulk is the source, wrong where
it is not.  No charge that is a function of `vgs - von(vbs)` can move it
(any such function responds to the source), so the depletion and
accumulation charge now takes its own argument, `mdep = vgb - von0`,
`von0` the threshold at zero bulk bias (level 3 also without DIBL, as its
gate charge already was) -- built by a new `hdl.with_values`, which
re-evaluates an intermediate with inputs substituted (the dependent part
of the chain rebuilt, as `ddx`'s forward pass).  `m` still sets the
near-threshold transition and the inversion charge, which have a channel.
At zero bulk bias in forward mode `mdep = m`; apart where `vbs != 0` or a
device runs reversed -- the regions there are the gate-bulk bias's (the
flat-band picture), Meyer's `vgs - von(vbs)` the named deviation now.

Measured: mosrect (each bridge device's bulk another bridge node) 0.19 V
from ngspice -> 0.002 V; gold 7.34 -> **0.917** (ngspice 0.939); mux8 1.05
-> **0.611** (0.342); every other gold deck unchanged (their devices'
bulk is their source).  nand stays 2.40 (0.55): another cause, open.
Test: `test_below_threshold_the_gate_capacitance_is_to_the_bulk` (levels
1-3, bulk 1 V from the source: deep accumulation all to the bulk, the
depletion ramp `-mdep*Cox/phi`; fails on the old charge).

(A measurement note: a gold run during which `elements_hdl.py` was
briefly swapped and restored read mosrect as 7.34 -- a class compiled
from the swapped file.  The rerun on an untouched tree is the record.)

## 2026-10-10 — 83. Meyer's channel split, smoothed, is the importer's default; gm17 runs (Andreas: "do as you suggest")

After the bulk-referred depletion (82) the two channel splits were scored
on every MOS gold deck.  Meyer's (`qpart = 1`) won on ten -- toronto 0.47
-> 0.32, todd3 1.06 -> **0.68 (passes)**, mux8 0.61 -> 0.40, nand 2.40 ->
1.32, gm1 6.0 -> 3.6, gm2 2.38 -> 1.35, gm3 2.19 -> 1.49, g1310 4.5 -> 3.0,
rich3 19.6 -> 16.7, ring 75 -> 66 -- lost on gm6 (4.14 -> 5.39), the rest
unchanged.  The importer's `mos_qpart` defaults to 1 now (the decks'
simulators use Meyer); the classes keep Ward-Dutton (0).

Meyer's split had a DEFECT, found on the way: its drain charge `-(Q(a,b)
- Q(a,0) + Q(b,0))/2` has slope `-coxt/3` in `b` at the saturation edge,
so `Cdd` jumped by Cox/3 there (measured 1.15e-14 -> 0 on a 10u/2u
device).  Now `qd = -Q(a,b) b^2/(a^2 + b^2)`: 0 in saturation, half at
vds = 0, symmetric, conserving, and its slope zero at the edge (3.5e-18
-> 0).  Its interior capacitances are not Meyer's -- named.

gm17 under Meyer's split collapsed its step at 5.75 ns -- not the jump
(it stayed with the smooth split), but slowlatch's mechanism: the error
test on its supply's current (with `lte_algebraic='skip'` it completes in
94 steps, no forced step).  The suite's `OVERRIDES` give gm17 the skip, as
slowlatch.  Its v(5) column stays ~96 under either split: open.
Passing the gold now: 4049osc, toronto, ab_integ, e1480, mosrect, mux8,
rca, schmitecl, latch, opampal, todd3 (11 of 23 scored).

## 2026-10-10 — 84. Level 2's gate charge saturates at its own vdsat; the MOS junctions' MJ unclamped; hussamp (Andreas: "1 first")

**The ring ran 4.8 % slow** (11.066 ns against the gold's 10.56; ngspice on
our deck 10.68).  Two charge defects, measured:

1. The square-law gate charge saturated where `vgs - von` says, Meyer's
   (SPICE's) at the model's `vdsat` -- which the ring card's VMAX and
   body effect put well below it, so past `vdsat` our drain side kept a
   Miller capacitance.  `_mos_gate_charge(sat=(vds, vdsat_f, vdsat_r))`
   scales the drain drive to reach 0 at `vdsat` (`b = a(1 - vds/vdsat)`
   forward, the mirror reverse), the two blended over |vds| < 25 mV by a
   C1 smoothstep (Meyer's MAGIC_VDS band): capacitances continuous through
   vds = 0 (max jump 6e-20 F over 2e-7 V).  LEVEL 2 ONLY: on level 3 the
   charge inherits its `vdsat`'s kinks and gm3 and jge collapsed their
   step -- level 3 needs a kink-free `vdsat` first (open).
2. The MOS bulk junctions clamped MJ at 0.9 (`_pn_depletion_charge`, the
   bipolar's rule); ngspice's MOS take MJ as given, and these cards give
   1.067 -- the junction capacitance was 14-39 % high at reverse bias.
   `mos=True` unclamps it (m = 1 itself kept off); now exact against
   SPICE's formula.  (The bipolar's 0.9 is 0.999 in `bjttemp.c`: open.)

Ring 11.066 -> 10.872 ns (ngspice 10.68: the rest, ~1.8 %, open -- Meyer's
interior capacitances and SPICE's non-conserving integration of them are
the candidates).  Gold, level 2 (vs before): ring 66 -> **48**, mux8 0.40 ->
**0.26** (ngspice 0.34), gm6 5.39 -> 5.24, g1310 3.02 -> 3.12, nand 1.32 ->
1.38, mosrect 0.917 -> 0.938 (passes); level 3 unchanged.

**hussamp passes (0.767)**: its `.NODESET V(1)=5` is on its supply, which
DC refused to hold; the importer's existing `.ic` rule (a source holds the
node: no effect in SPICE, left out and said) now covers `.nodeset` too.

**The supply-current error test, now on three decks** (slowlatch, gm17,
the ring: with `lte_algebraic='skip'` the ring 7885 steps / 729 forced ->
2877 / 0 at the same accuracy): still an opt-in per deck (Andreas's
choice); the ring's PSS waits on it -- its `lte_grid` drew 25,563 points
from the forced-step run and the shooting ran out of memory.

## 2026-10-10 — 85. Stage 9 continued: the ring's PSS, voter25 runs (level 2's VMAX vdsat floored), `spice_ring`

**The ring's PSS** converges: 10.877 / 10.856 / 10.852 ns at 320 / 1013 /
2864 grid points (gear), against ngspice's 10.68 on our deck -- the gate
charge's open 1.6 %.  The key was `lte_grid`'s reference: its per-unknown
default (`relref='pointlocal'`, kept deliberately for switched-noise
grids) at reltol 1e-5 drew 24,336 points through the nine fast edges and
the shooting ran out of memory; the GLOBAL reference at 1e-4 draws 320
and refines to the same period.  `PSS(lte_algebraic=...)` now reaches
its transients.  Test: `test_the_ring_oscillator_s_periodic_steady_state`
(5 s).  `spice_ring` joins step_machinery (the ring's transient to 20
ns, 0.43 s).

**voter25 stalled at 32.6 ns -- and ngspice stalls on it too (Xyce runs,
8.9 s).**  An off PMOS (|vgs| 0.5 V, vds 22 uV) whose gate sat at `vbin`:
with VMAX the quartic gives `vdsat = x^2 - (phi - vbs)`, NEGATIVE below
threshold (-0.33 V there), which ngspice does not floor (its Grove-Frohman
branch does), and the saturation arm then carries a spurious 1.7 uA,
switched on abruptly at `vgs = vbin` -- a jump Newton cannot cross.  The
VMAX `vdsat` is floored at 0 now (Andreas: "make level 2 continuous"):
the current is 0 below threshold and grows continuously.  A named
departure from ngspice; the reference transcription gains `floor=` (the
model's comparisons use it, the ngspice one does not).  voter25 then runs:
36.5 s, 14463 steps / 890 forced; with its supply currents out of the
error test 10.4 s, 3989 / 2 (the FOURTH deck: slowlatch, gm17, the ring,
voter25).  Against Xyce's waveforms (no gold): the logic levels agree at
every sampled time; Xyce's metric with a 1 V offset 4.6-6.2 -- edge
timing, the same open charge gap.

mem_plus (17,789 unknowns) imports in 19 s at 468 MB; its transient is
not run yet.

**mem_plus (2026-10-10)**: imported (21,743 elements, 17,789 unknowns, 19 s,
467 MB); its operating point did not finish in 25 minutes (a 10 GB cap,
~3.1 GB reached) -- the size of ONE dense 17,789^2 matrix (2.5 GB) with an
LU of ~1e12 operations an iteration.  The large decks need the sparse
path: the evidence for stage 10's decision.

## 2026-10-10 — 86. The channel's split starts at 1/2 where inversion starts; gm17 runs under the full error test (Andreas: "how can we solve this without loosening the tolerance", "start gm17")

**gm17 was a model corner, not the error test.**  Judged (supply currents
in the error test) it stopped at 5.75 ns, the step at 1e-18 s, Newton
failing.  The drain capacitance of one PMOS (`x21:mp1`, reverse, |vds|
28.7 mV, at threshold) jumped by 4.86e-14 F under qpart 1 and 0.98e-14
under qpart 0 within 0.1 mV -- the ratio 5 that the mechanism predicts:
where a device crosses its threshold with ONE side, below it the
transition charge `qtr`'s image is halved between source and drain, above
it the channel's split gives the inverting side all of the channel
(Meyer; Ward-Dutton 3/5), so `dqd` jumped by coxt/3 (coxt/15) -- forward
as well as reverse.  Now the split starts at 1/2 there and is the chosen
split from a drive of 25 mV on: `qdch = -qinv/2 + h (qdch_split +
qinv/2)`, `h = 1 - (1 - s)^2`, `s = min((a^2 + b^2)/(25 mV)^2, 1)` -- a
function of a^2 and b^2, flat at both edges, exact Meyer / Ward-Dutton
beyond 25 mV (the existing split tests hold to 1e-12 COX).  Named
deviation (6) in `_mos_gate_charge`.  Test:
`test_the_capacitances_are_continuous_where_one_side_crosses_threshold`
(L1 and L3, both splits, both directions) -- on the old code it fails
with jumps of exactly 0.0667 and 0.3333 COX.

**The first fix was wrong, and the eigenvalues said so.**  Moving `qtr`'s
image with vds instead (Meyer's 0 forward / 1 reverse, blended over the
band, tapered by `(1 + 2m/phi)^2`) removed the jump but put a
drain-source transcapacitance into near-threshold devices: gm17 then
crawled at 5.77 ns with every step ACCEPTED (0.83 of tolerance) and
shrinking, `br0`'s estimate GROWING as h fell -- a growing mode.  At that
state `eig(-G, C)` had lambda = +1.17e13/s on `x21:mp1`'s internal
nodes, its C indefinite (internal self-capacitances 54/82 fF -> 23/46
fF); with the halved image at the SAME state: no unstable mode.

**Judged, after the fix**: gm17 completes (517 steps, 26 forced, metric
96.06 against the skip's 96.11); slowlatch 1025 / 22 (136.23 both). The
remaining forced steps sit on no capacitance or conductance jump (every
device scanned +-2 mV at each forced state): the supply current's
estimate falls only ~linearly with h at 0.2-60 ps steps (stiff, h >> the
~0.3 ps internal-node time constants) while the controller sizes its
retries for h^3, so four retries end at 1.0-1.25x tolerance and the step
is forced -- the step controller's order assumption, the next item.  The
ring's forced steps (7824 / 657) are that, unchanged by this fix.

**The gold decks** (suite overrides unchanged): gm3 and jge, which
collapsed at the time-step floor, now RUN (gm3 metric 1.49, ngspice
0.65); gm1 3.87 -> 3.62, rich3 17.2 -> 16.9, toronto 0.338 -> 0.325;
the rest unchanged or within noise (g1310 3.12 -> 3.16).  voter25 with
the skip 3995 / 1 (was 3989 / 2).

**Open: voter25 judged now crawls at 40.04 ns** (it ran before, 36.7 s,
890 forced).  Not the error test (1e-5 of tolerance): Newton, in a
period-3 cycle of 1-2 mV on `xfblk0#b:xload:mp5.si` / `mp6.di` at h =
6.4e-17 s, converging at half that.  No C or G jump along the cycle's
path, C and G agree with finite differences to 1e-9, no physical
unstable mode.  Not yet explained; a suspect (from reading, unverified):
level 2's reverse drive `b0 (1 + vds/vdsatr)` blended into the 25 mV
band, with `vdsatr` near threshold tiny.

## 2026-10-10 — 87. The step controller's second retry shrinks at the order its tries measured (Andreas: "start on the step controller retry order", "build it, take on the JAX side later")

The forced steps that stayed after 86 were the retry law, not the model:
a supply's branch current falls as ~h^1.2-1.5 between retries at one
point (q ~ 0.85-1.2 in the sequences that ended forced) where the law
assumed h^3.  From the second retry the exponent is the one the two tries
measured, clamped to [1, p]: the ring judged 657 -> 0 forced, mux8 624 ->
131, gm1 59 -> 16, gm17 judged 26 -> 2, slowlatch 22 -> 10, the time the
same or less; native circuits unchanged.  The measurements and the
design: `doc/transient_history.md` (`RETRY_ORDER_MIN`).  The JAX transient
has it since (Andreas: "do it"), and both take `retry_order` --
'assumed' / 'measured' (default) / 'carried' (Andreas: "implement but
controlled by an option"; 'carried' trades 20-27 % of the rejections for
8-17 % more steps and no fewer forced, so it is not the default).
The one suite failure, van der Pol's LTE-grid PSS test, was the test's
window (short of a period by up to one coarse step; it passed or failed
by where `tend` fell, on the old code too); built exactly, 12 of 12
converge.  `PSS.lte_grid` is not affected (`doc/transient_history.md`).

## 2026-10-10 — 88. Level 2's gate-charge saturation floor is Meyer's 25 mV, not 1 uV -- voter25 runs judged (Andreas: "investigate and fix")

**voter25 judged crawled at 40.04 ns** (2-100 as steps, every try
accepted by the error test at 1e-5 of tolerance): the step's Newton
cycled.  ⚠ The first diagnosis (85, then 86's open item) was VOID: the
saved Newton iterates were the solver's REDUCED vector (ground row 13
removed), read against the full vector's names -- every unknown past
ground was its neighbour, and "no C/G jump, the Jacobian equals FD" was
measured on misaligned vectors.  Redone aligned, then with the solver's
own residual function captured: J equals FD of the residual at every
iterate, yet the residual along the cycle's segment is linear with one
STEEP RISE -- the node-3 (clock) charge gains 1.06e-17 C over 2.9 uV of
one PMOS's (`xcblk0:xload:mp6`) gate drive, a picofarad on a 6 fF device.
Seen at the deck's 125 C (`tr.epar`), invisible at `defaultepar`'s 27 C
(the first per-element scans used that).  A guessed cause was REFUTED
(the VMAX vdsat's no-root branch: changing it moved nothing, both
backends, each tree importing its own source).  The cause: the gate
charge's saturation block divides by the model's vdsat floored at 1 uV;
in reverse mode the forward arm `a0 (1 + |vds|/vdsat)` (weight 0.34 at
-5 mV) sweeps through thousands of `a0` as the forward vdsat leaves 0
near threshold.  Floored at Meyer's MAGIC_VDS, 25 mV, rounded C1 (25 mV
+ vdsat^2/0.1 V below 50 mV, vdsat above), the step is gone (largest
charge step = the median one).

**voter25 judged completes**: 16614 steps, 107 forced, 36.7 s (it crawled
indefinitely; before 86 it ran, 14463 / 890 forced).  The level-2 gold
decks: unchanged except mux8 0.259 -> 0.263 (131 -> 145 forced) and nand
29 -> 33 forced.  Test: `test_level2_s_saturation_floor_is_meyer_s_not_a_
microvolt` (the two states, 2.9 uV apart: 1.06e-17 C on the old code).

**Open, found on the way: the band's near-threshold gap.**  At the same
state the gate's own capacitance still JUMPS ~0.4 coxt (2.3 fF) within
5 uV (HEAD: 1.4 fF): inside |vds| < 25 mV the rounded `m` ends `qtr` a few
mV before inversion starts -- the docstring's deviation (2), where Meyer
spreads the transition x1.5 near vds = 0.  A capacitance jump, which
Newton crosses (voter25 runs); the candidate fix is Meyer's sliver plus
drives measured from `m`.

## 2026-10-10 — 89. Stage 10: the decision's numbers (the rule as the plan wrote it before measuring)

**The rule** (the plan, stage 10): C stepping loop next if the median
orchestration share on the gold cases is >= 40 % AND the circuits wanted
fit the dense limit (expected gain ~1/(1 - 0.8 share)); sparse path next
if the wanted circuits (mem_plus, ram2k, the IBM grids) are refused by the
guard, or assembly + linear algebra is >= 50 % of a step at their size;
vectorised sources first if the per-source pass dominates the mesh steps.

**The split, measured** -- instructions in the timed solve only
(cachegrind switched on and off around it, `_bench.CountRegion`), every
function attributed to the shared object that defines it; OpenBLAS pinned
to one thread; 26 gold decks with n <= 429 (the long ones capped at ~2 s
native, hundreds of steps each).  Attribution checked by callgrind on two
decks: the unnamed code is the interpreter binary's own; memcpy/memset is
half numpy's copies, half the in-memory C core's matrix traffic; libm is
the kernels' math.  Unattributed at most 1.2 %.

| share of a step (26 decks) | median | quartiles |
|---|---|---|
| the Python interpreter alone | 36.1 % | 13.4 / 43.9 % |
| orchestration: interpreter + numpy + libc + half the memory traffic | 49.1 % | 26.4 / 59.0 % |
| assembly + linear algebra: BLAS/LAPACK + half the memory traffic | 25.0 % | 12.0 / 34.9 % |
| C kernels (HDL kernels, C Newton/LTE cores, libm) | 25.4 % | 13.7 / 38.0 % |

n <= 100 (16 decks): orchestration 54.9 %, assembly + LA 15.6 %; n > 100
(10): orchestration 31.8 %, assembly + LA 34.2 %.  Dense LU already
dominates some 100-300-unknown decks: e1480 (n 205) 80 % BLAS, ab_integ
(95) 68 %, g1310 (98) 65 %.  Instruction shares understate the
interpreter's TIME share (its IPC is the lowest of the four).

**The other inputs** (stage 5, stage 9): the dense wall ~10 k unknowns on
this 30 GB box; on the RC mesh assembly + LA is 52 % of a step at n 1190,
77 % at 9226; the per-source pass 0-8 % (not dominant: vectorised sources
are not first); every IBM grid above 1/16 refused (whole, merged: ~61 GB);
mem_plus (17,789 unknowns) has no operating point dense in 25 min; ram2k
(13,880 MOSFETs) and pc_frame (14,265) are past the wall.

**What the rule says.**  Orchestration median 49 % (36 % counting the
interpreter alone) -- over the C loop's 40 % line by the broad reading,
under it by the narrow one; its gain ~1.65x (1.4x narrow), on circuits
that fit.  The circuits the plan names as wanted -- mem_plus, ram2k, the
IBM grids -- are all refused dense: if they are wanted, the sparse path is
next (capability before a constant factor).  Presented to Andreas.

## 2026-10-10 — 90. The C stepping loop plan, stage 0: sized in cycles -- cleared its line (Andreas: "starts with the C stepping loop", "write the plan", "go ahead")

The plan: `~/.claude/plans/c-stepping-loop.md` -- `Transient`'s multistep loop
in C, sparse-ready (the matrix through an interface with dense, chord and KLU
backends), stage 0 a measurement with its line written first: build only if the
cycle-weighted conservative ceiling is >= 1.5x in median over the 26 gold decks
and >= 1.8x over those with n <= 100.

**The instrument** (`benchmarks/loop_split.py`): each deck in its own process,
pinned, warm; the CPU's instruction and cycle counters read at the boundaries of
`_SteppingLoop`'s five phases and of every COMPUTE call (the C Newton, the C error
test and the HDL batch drivers' functions, a `CKernel` call, the leaf linear
solvers, dense LAPACK through `numpy.linalg` / `scipy.linalg`); the readings' own
cost calibrated in the process and taken off.  REMOVABLE = in the loop and in no
compute call; CEILING = 1 / (1 - removable).  Every counted run's fast-path
counts are compared with an unwrapped run's -- and three wraps were CAUGHT moving
the path, then removed: the C core's Python callers (the C Newton declines when
they are not genuine, `newton_c:patched`; they are wrapped only on a deck the C
Newton never serves), and `numpy.linalg.solve` (the C error test twins it,
`lte_c:patched`).  Proved by a planted cost: a busy loop added to one phase came
back in it to -0.09 % (accept, gm6) and -0.40 % (take, mosrect).

**The result**: ceiling median **1.70x in cycles** over 26 decks (quartiles 1.34 /
3.06), **3.00x over the 16 with n <= 100**; in instructions 1.28x / 2.00x -- the
interpreter's cycles per instruction are the highest of the four.  CLEARED on
both.  Per deck (cycles): 4049osc 9.35x, gm2 4.25x, todd3 3.38x, gm6 3.32x ...
ring 1.68x, gm1 1.72x, voter25 1.44x, rich3 1.35x, gm17 1.17x (the operating
point is 64 % of it), gm19 1.06x (76 %).

**Found on the way**: the branch check's per-step SVD screen (`_branch_screen_C`
-> `numpy.linalg.svd`) is 51 % of ab_integ's run and dominates e1480 and g1310 --
compute, its own item.

**Census** (the loop's scope -- why it would NOT run): the 26 gold decks all in
scope; the C Newton serves on 23, the three n >= 250 decks (rich3, arom, gm19) on
KLU (`newton_c:linsolver`), stage 2's.  The suite (a pytest plugin wrapping
`_SteppingLoop.__init__`; its one failure, the pinned pair "the stepping loop and
its unread Jacobian", is the plugin's own patch -- 58/58 without it): 1187
stepping loops in 544 tests, 632 (53 %) in scope; declines: the stage family 334,
stateful limiters (Diode) 125, `accept_step` 59, the coupled family 52,
`provided_function` 41, tests' own controllers 26, PI 9, state events 7, the
band 2, a source 1, the toolkit 1.  Tests: 255 with every loop in scope, 26 with
some, 263 with none.

## 2026-10-10 — 91. The C loop plan amended after a measurement and four reviews; stage R1a, the C-against-C guard test (Andreas: "measure and deep review first and amend the plan ... not just speed ... also to make the code understandable")

**Measured** (`insitu.py lines`, instructions per call of `_tran_newton_c.solve`): on
spice_4049 (n 15, fixed step) 407 k a call -- the C call 114 k (28 %), the branch check's
screen after it 87 k (21 %), the source pass with a 21-key counter snapshot for rollback 65 k
(16 %), eligibility / type checks, the temperature probe, the bypass and limiter decisions ~75 k
(18 %), marshalling (history rows, 5-7 copies in, 6 out of which 3 n x n, `insert_row`) ~55 k
(13 %).  On spice_ring (n 88, adaptive): 3761 k, the C 3335 k (89 %).  The error test's wrapper:
735 k on the ring, 658 k in C.

**Reviewed** (read-only, four areas): eligibility re-decided every call (28 decline points in
the C Newton's wrapper, 23 in the error test's, each rebuilding a tuple of the method chain to
compare); the loop's state on ~45 fields; the companion setup and its 8-attribute write-back
three times; couplings a C twin must copy (the retry's "same point" by object identity, the
excursion veto after a recorded accept, the order-drop signal's four hops); two SCOPE GAPS in the
first plan (the branch check screens after every served Newton; the C error test declines every
Euler step, i.e. after every breakpoint); a WRONG FACT in it (n 250..999 uses SuperLU, not KLU --
checked: rich3, arom, gm19 `solver.auto:SuperLUSolver`; and today's KLU pattern is
value-derived); the C side: 14 units in Python strings, 7 structs written twice with nothing
comparing them, the core's 21-argument signature re-typed in 5 places, `dgetrs_t` LP64 in two
units and ILP64 in a third, magic status / pass / formula codes, the HDL prelude compiled into
every runtime object, the `.c` deleted after the build, an 11 GB object store never pruned; and
the gate counts fast paths, which a C loop moves.

**The plan amended** (`~/.claude/plans/c-stepping-loop.md`): two clean-up stages before the
sources -- R1, the C side (one home, real files, one header declaring every interface once; no
answer moves) and R2, the Python side (one run context checked once per run, one declared list
of the loop's state, one companion setup, the couplings explicit; line: the per-call Python
around the C Newton -35 %) -- then sources, the linear system on a structural pattern, the error
test's Euler formula and the branch proxy in C, the fixed-step loop, the adaptive loop.  Code
clarity is an acceptance criterion of every stage.  OPEN for Andreas before stage 2: the n >= 250
runs decline the C loop, or every sparse run moves to KLU on the structural pattern (one
deliberate re-baseline).

**R1a** (`tests/test_c_abi.py`): each of the 7 structs parsed by cffi's own parser twice -- from
its cdef and from the typedef in its C text -- and compared field by field (name, offset, size,
kind; 301 fields, every size equal); each caller's copy of the core's (21 arguments) and the
limiter walk's (11) signature compared with the callee's own definition and declaration by ABI.
All agree today; alive on a planted missing field, an extra field, a narrowed argument and two
swapped arguments.

## 2026-10-10 — 92. Stage R1b: the C build and load machinery in its own module, `_cbuild.py` (no object rebuilt)

The flags, the sanitizer, compiler discovery, the object store, the key, the lock, the
build, `dlopen` and `load_kernel` moved out of `_hdl_cbackend` (an HDL module) into
`_cbuild.py`, CUT from the file byte for byte, under a docstring that states its contract.
`_hdl_cbackend` re-exports the names for its callers and reaches them as `_cbuild.X` itself.
Proof: the 14 recorded keys (the 12 runtime C sources, the flags, the default entry cdef)
identical before and after -- every stored object is served as it was.  Tests that PATCH a
moved name (`CFLAGS` x2, `_loaded` x7, `_compiler` x3) now patch it on `_cbuild`, where the
code reads it; a patch left on the old name fails `test_every_key_ingredient_changes_the_key`
(planted, then restored).  The pinned pair "the printed pass kernels and their fused kernel"
moved on its twin (`bind_fused`: `load_kernel` / `CompileError` named through `_cbuild`, the
same objects) -- checked, re-recorded (twin digest only).  Full suite: 5027 passed and that
pair.

## 2026-10-10 — 93. Stage R1c: the runtime units' C in real files, `pycircuit/circuit/csrc/*.c` (no object rebuilt)

The twelve C texts the runtime builds -- `newton.c`, `lte.c`, `core.c`, `radau.c`,
`radau_tc.c`, `radau_fold.c`, `sens.c`, `pass.c`, `walk.c`, `limit.c` (the limiter's
template), `psp_limit.c` (PSP's, %-formatted per value as before), `watch.c` -- moved out
of their Python string literals into files, BYTE FOR BYTE (each literal asserted equal to
its value before the move, raw strings only, no continuation); each module reads its file
through `_cbuild.csrc(name)` (no newline translation) into the same attribute.  About 1350
lines of C leave the Python.  Proof: the 14 recorded keys identical, and the pinned pairs
and the ABI test pass WITHOUT re-recording (they digest the strings' values).  The files
ship with the package (`package-data` in pyproject.toml and setup.py; a built wheel checked:
all 12 inside).  Full suite 5028 passed.

## 2026-10-10 — 94. Stage R1d, step 1: one interface text per unit, the shared types once (machine code identical)

Each runtime unit `X.c` now has a header `csrc/X.h` -- its struct (cut from the `.c` byte for
byte) and its `hdl_fn` declaration (cut from the module's cdef).  The module's cdef IS that
header's text (`_cbuild.csrc('X.h')`) and `X.c` includes it, so cffi lays the struct out from
the very text the compiler compiles, and the compiler checks the definition against the
declaration cffi reads.  Before, each struct was written twice (C and cdef) and R1a's test
compared them; that comparison would now compare a text with itself, so the test asserts the
construction instead (the cdef is the header, the unit includes it once, no `.c` declares a
type, no header includes another).

The function types the units passed each other by hand-written copies -- `core_fn_t` and
`walk_fn_t` (3 copies each), `dgesv_t` (4), `dgetrf_t` and `dgemv_t` (2 each), and `dgetrs_t`
(ONE name for TWO signatures: 32-bit integers in newton.c and sens.c, 64-bit in radau_tc.c) --
are declared once in `csrc/pcx_abi.h`, under names that say what they are: `pass_fn_t` (core.c's
`hdl_fn_t` and pass.c's `hdl_batch_fn_t`, the same contract), `fused_fn_t`, `lim_fn_t`,
`dgesv64_t` / `dgetrs64_t` / `dgemv64_t` (numpy's ILP64 OpenBLAS), `dgetrf32_t` / `dgetrs32_t`
(SciPy's LP64), `klu_refactor_t` / `klu_solve_t`.  `_cbuild.csrc` expands `#include "x.h"` in
place, so a unit's key and its pinned pair cover every line compiled.

Proof: every unit's MACHINE CODE identical -- the `.text` and `.rodata` bytes and the exported
symbols of all ten units, compiled as the store compiles them, before and after (an A/A of the
instrument identical first); each new cdef declares what the old one did (normalised text and
cffi's layout); no unit gains a `-Wall -Wextra` warning.  The 12 pinned pairs that digest a
unit's text moved on the twin only (re-recorded; no reference moved).  `limit.c` and
`psp_limit.c` untouched (no model limiter rebuilt).  Full suite 5028 passed; the sanitizer's
fast tier clean (4722 passed).

## 2026-10-10 — 95. Stage R1d, step 2: the codes named once, in the headers (machine code identical)

The numbers the C units return and are called with were bare in the C and mirrored by hand in
Python: each unit's results (`BAIL` tables in five modules, `status != 1`, `status == 8`), the
core's passes (`want`: 1 G, 2 C, 4 i; the fused mask's q bit 8; `_PASS_BITS`, `WANT`), its
formulas (`FORMULA`: 0 gear2 .. 3 theta, 4 the passes alone), the error test's reference modes
(`_RELREF`) and CBLAS's `101, 111`.  A real trap among them: "the limiter walk stopped" is 8 in
newton.c and radau.c and 10 in radau_tc.c -- the Python tested each by its own number.

Now each unit's header names its results (`enum newton_status { NEWTON_CONVERGED = 1,
NEWTON_BAIL_X = 2, ... }`, a line saying what each means; lte.h also `enum lte_relref`), and
pcx_abi.h names the core's (`enum pcx_pass`, `enum pcx_formula`, `enum pcx_cblas`).  The C says
the names; Python reads its tables from the header text (`_cbuild.constants`, `_cbuild.bails`):
`BAIL = _cbuild.bails(_CODES, 'NEWTON_BAIL_')`, `status != _CODES['NEWTON_CONVERGED']`,
`FORMULA` / `WANT` / `_PASS_BITS` from pcx_abi.h -- the strings the Python reports (`'walkstop'`,
`'getrf'` ...) are the names' suffixes, unchanged.

Proof: the machine code of all ten units identical again; every Python table equal to its old
literal (the eight compared as data); cffi reads each header's enums as `constants` does.
test_c_abi gains: the two readers agree on every header; pcx_pass's bits are `1 << k` of the
pass order Python builds the fused masks from (`_tran_core.PASSES`, `_hdl_cbackend.FUSE_PASSES`);
no `.c` says a code by its number (the check finds 80 in the previous commit's C, 0 now); an
enumerator spelling no value refuses.  11 pinned pairs moved on the twin only (pass.c and walk.c
among them: they include pcx_abi.h, which gained the enums -- the watcher's, which does not,
stayed).  Full suite 5040 passed; the sanitizer's fast tier clean (4776 passed).

## 2026-10-10 — 96. Stage R1d, step 3: the runtime units built without the HDL kernels' prelude

Every object was built from `hdl._KERNEL_C` + its source: the prelude of the printed kernels
(`<math.h>`, the libm functions re-declared `__attribute__((const))`, `_sel`, `_npmax`,
`quartic_root` ...).  The runtime units used none of it knowingly, but one depended on it:
`sens.c` calls `isfinite` and never included `<math.h>` -- the prelude did.  Now
`_cbuild.load_unit(csrc, cdef)` builds a unit as the complete C file it is (`prelude=False`
through `source_key`, `_build`, `load_kernel`; the key hashes an empty prelude part, so no HDL
kernel's key moved), `sens.c` includes `<math.h>`, and the ten unit drivers no longer import the
HDL backend to build (`_cbuild.load_unit`, `_cbuild.CompileError`).

Proof: nine units' machine code identical (the objects the store built, hashed).  `radau_tc`'s
differs, by one thing: the prelude declared `hypot` `const`, and it is not -- it raises the
floating-point flags this unit tests.  Its four `hypot` calls sit between a `feclearexcept` and
the `fetestexcept` that judges them, each result used before that test, so data dependencies
already fix their place: the same four call sites either way, 26 instructions fewer without the
attribute; the values and the flags cannot move.  The 102 Radau tests (their C-against-Python
bit-identity among them) pass.  test_c_abi gains: every unit compiles ON ITS OWN under `-Wall
-Wextra -Werror` (alive: `sens.c` without its include is refused, implicit `isfinite`).  The two
tests that build a planted variant of a unit (core.c, pass.c) build it as the unit is built.
One pinned pair moved on its twin (sens.c's include).  Full suite 5050 passed; the sanitizer's
fast tier clean (4758 passed).

## 2026-10-10 — 97. Stage R1d, step 4: `#line` -- a report names the repository file and line (R1 done)

`_cbuild.csrc` marks each header it expands: `#line 1 "pycircuit/circuit/csrc/x.h"` before it,
`#line N "pycircuit/circuit/csrc/unit.c"` after it (the path from the source tree's root, the
same wherever the package lives -- it is part of the key).  A unit is compiled with nothing before
it, so its own lines need no other marker; a file including no header (limit.c, psp_limit.c) is
its text exactly -- no model limiter's key moved.  A failed build or a sanitizer report (the
sanitized builds carry `-g`) named a temporary `.c` that `_build` had already deleted; it now
names the file and line in the tree.

The plan's other half -- the `.c` kept beside its `.so` -- is not built, by reason: a unit's key
is the hash of its exact text, so an object served for the current source was compiled from
exactly the files its markers name; a kept copy would add a file per object (4566 in an 11 GB
store) for nothing a report needs.  The HDL kernels (generated text, no file in the tree) are where
a kept `.c` would help; left until a report in one needs it.

Proof: the machine code of all ten units identical (no `-g`, `#line` emits nothing).
test_c_abi gains: an `#error` planted before the last line of newton.c, of newton.h and of
pcx_abi.h is reported at that file and line (alive: with the markers stripped the report names
`unit.c:193`; a first plant at a file's middle sat in a block comment, where `#error` is inert).
12 pinned pairs moved on the twin only.  Full suite 5053 passed; the sanitizer's fast tier clean
(4762 passed).

R1 (stages R1a-R1d) is done: the C units are files, one interface text each, shared types and
codes named once, complete without the HDL prelude, and reported at their own lines -- every step
proved by unchanged machine code (radau_tc's one difference: the prelude's false `const` on
`hypot`).

## 2026-10-10 — 98. R1's line held, and the one cost it hid: radau_tc's residual loops reloaded after every `hypot`

R1's line (written before measuring): bit-identical everywhere, `step_machinery --check` no SLOWER,
against the tree before R1 (02235bcc).  Held: all 18 cases "no change", bytes, statistics and
fast-path counts SAME.  But the child leaned slower in most rounds (+0.1..+1.1 %, some intervals
clear of zero: ladder_gear +1.1 % [+0.5, +1.7]) -- under the check's band (+3 % flags every case,
+1 % none), so counted rather than shrugged off: `--count` (instructions, perf) flagged ONE case,
mos1_radau's step +1.00 % (2.759 M -> 2.787 M), with ladder_radau +0.33 % (spread 0.01 %); every
other case within its tree's own spread.

The cause is step 3's one machine-code change.  The prelude declared `hypot` `const`; without that
lie, the compiler must assume each `hypot` call may write memory, so after every call in the
residual check's loops it reloaded the struct's array pointers (`s->Ax`, `s->Ai`, `s->y1`,
`s->ab`, the bound `s->Ap[j + 1]`).  Fixed by saying it honestly: the arrays in local pointers and
the inner bound read once (the arithmetic and its order untouched).  Recounted against the same
parent: mos1_radau's step -0.10 %, ladder_radau -0.07 %, pss_radau -0.03 % -- no MORE anywhere.
Lesson: a time verdict of "no change" with a one-sided lean is worth one `--count`; it named the
case and the cause in a minute.
Two pinned pairs moved on the twin (radau_tc.c) and were re-recorded; the 103 Radau tests pass.  Full
suite 5053 passed; the sanitizer's fast tier clean (4746 passed).

## 2026-10-10 — 99. Stage R2a: the stepping loop's state, measured and declared (`transient._RUN_STATE`)

What does one step hand the next?  The C loop must marshal exactly that, and no reading of the
code had listed it.  Measured instead, two ways, with `sys.monitoring` on `_SteppingLoop.execute`
and `accept` (nothing on the loop replaced: wrapping a loop method makes `_paths.genuine` false
and changes the very path observed):

* a WRITE census over the 20 benchmark cases (13 with a stepping loop; PSS shoots outside it):
  every attribute reachable from the transient fingerprinted at each accepted step (identity,
  array bytes) -- 39 transient attributes written between steps, 11 on the loop, the rest in
  the core's buffers, the limiters, the linear solver, the integrator and the controller;
* a FIRST-ACCESS trace (class-level `__getattribute__` / `__setattr__` on `Transient` while a
  loop runs): per step, is an attribute's first access a read (it was CARRIED in) or a write
  (it is STEP scratch)?  ~80 attributes touched a step: 14 read before written, 16 written
  before read, the rest the run's context (mostly methods).

Crossed (a "carried" one written once a run is a lazily memoised RUN constant), they give three
roles, now declared once in `transient._RUN_STATE` (49 names, each with its role and a line of
meaning):

* carried (15): the history (`_qlast`, `_iqlast`, `_dt_last`, `_dt_last2`, `_is_first_step`,
  `_no_history`), the predictor's (`_pred_hist`, `_pred_fast`), `_C_cache`, the device memo,
  the branch check's bound, the Nordsieck GLM's four;
* run (21): the C paths' records and decline memos, the lazily decided verdicts
  (`_has_limiters`, `_jacobian_expensive`), the caches, the linear solver, the loop's mark;
* step (16): the companion block (`_q_cache`, `_iq`, `_Geq`, `_Cmat`, `_companion_coeffs`,
  `_effective_method`, `active_integrator`) -- written in FOUR places, not the three the plan
  counted (the Python reference `_tran_companion`, `_tran_core`, `_tran_newton_c`,
  `_tran_stages`) -- with `_dt`, the source memo, the stages' values.

Found on the way: `statistics_glm_startups` is a counter kept ON the transient, not in its
`statistics` (R2d's list).  `test_run_state.py`: 14 runs covering every step family (each
integrator adaptive on a pulsed RC, the coupled family, a fixed grid, an HDL MOS chain fixed /
gear / Radau through the C paths) -- any attribute whose binding changes while a loop runs must
be declared; every declared name has a role and is assigned in the package; the observer proved
alive on a planted write.  No code path changed (the manifest is data; the tests observe).  Full suite 5069 passed.

## 2026-10-10 — 100. Stage R2b measured and REFUTED: one run context per call does not pay

R2b's premise: the C Newton's ~28 per-call checks become ONE run context, made once a run and
re-validated by ONE stamp test per call.  Sized first (`insitu.py lines spice_4049
_tran_newton_c:solve`, 408.5 k a call): the C 113.7 k, the branch screen 87.8 k (stage 3's), the
source pass 49.2 k (stage 1's), the run-level re-checks ~95 k.  R2's line was RE-SIZED before R2b
was built (the plan, 2026-10-10): the R2-addressable ~158 k down >= 35 %.  A throwaway ceiling
prototype (every run-level check skipped, unsafely) put the most R2b could save at 46.5 k a call.

Built exactly: a `_Run` context after a call whose every check passed, at the counter's present
value; its test -- the counter, plus ONE pass over what no dict watcher sees (the transient's own
attributes the checks read, the machinery's class and module attributes, the switches,
`ParameterDict`'s epoch, the linear solver's choice, values in unwatched dicts, each element's
`__dict__` object, each pass's code, the circuit's size).  First run: built but held 0 of 400
calls -- the step's source memo (`_RUN_STATE` 'step') is a new dict every step, and the context
compared it by identity where the check asks only that one exists (fixed: held 397 of 400).
Then COUNTED against its parent (`step_machinery --count`): spice_4049 +1.24 %, mos1 and stage
+1.0 % a step -- MORE.  Bytes and statistics the same.

Why: the stamps already skip what a dict watcher can see; what is left on every call is exactly
what NO watcher sees, item by item -- and the fused test must compare every one of those too (~70
identities through mapped tuples), which costs what today's inline checks cost, plus the
context's own reads.  The work a context could truly skip (the limiter closure built to test `is
not None`, the walk's lookup, the core record's call) is a few k.  Reverted, by the plan's rule (a
stage showing MORE is not kept, even for clarity).  The design is not lost: validated ONCE PER
C-LOOP CALL -- many steps a call -- it is the C loop's entry gate (stage 4), where it pays; the
diff is kept with the plan's working notes.

## 2026-10-10 — 101. Stage R2c, batch 1: per-call work the C Newton's entry did for nothing (-2.6 % a step)

Exact changes, counted against their parent (`step_machinery --count`, instructions; bytes and
statistics the same everywhere):

* the source-count snapshot taken on EVERY call only to roll back a bail -- a dict built by a
  comprehension of `_PC.get(k, 0)` over 21 keys (15.9 k a call in situ) -- now the values,
  `tuple(map(_PC.get, keys, repeat(0)))`, the dict made only on the rare rollback (the form
  `_tran_radau_tc` already used); in `_tran_newton_c` and `_tran_radau_c`;
* `x.dtype == np.float64` on the per-call checks (numpy's dtype coercion, ~50 ns a test) as
  `x.dtype is _F64 or x.dtype == np.float64` (~19 ns; the same verdict for every array): the
  seed, the source vector, the cached state, `_tran_core._row`, `evaluate`'s and `passes`'
  inputs;
* the converged `x` is the C's own full-width `xc` (the reference row's +0.0 in place), copied,
  where `insert_row` rebuilt it from the reduced solution (7.6 k a call) -- the same array byte
  for byte (checked on 99 served solves before the change).

spice_4049 -3.49 %, mos1 -2.63 % a step, stage -2.72 % a step, mos1_adaptive -1.52 %; the Radau
cases no change (another path); no MORE.  On spice_4049, ~24.5 k a call of R2's re-sized
~55 k line.  Five pinned pairs moved on the twin (the Newton's, the core's, the coupled stage's)
and were re-recorded; one, "the step's end passes and the C's", moved on its REFERENCE
(`_tran_core.passes`, its input guard's dtype test): checked -- the C twin makes the calls
`_stage_end_passes` makes, which a guard's spelling cannot change -- and re-recorded.  Full
suite 5064 passed and those five.

## 2026-10-10 — 102. CORRECTION to 100 and 101: `--count` has a 1-1.7 % offset between identical trees -- R2b was neutral, not MORE

An A/A (`step_machinery --count` with the parent tree at HEAD and this tree identical to it,
checked file by file): spice_4049 +1.66 % / +1.69 %, mos1 +1.20 % / +1.22 % (its step +0.89 /
+0.88 %), stage +1.17 % / +1.19 % (its step +1.22 / +1.12 %), mos1_adaptive +1.41 / +1.33 % --
the child MORE, the same within 0.1 % run to run: a SYSTEMATIC offset, not noise.  Roles
swapped (this tree the parent, the worktree the child): -0.2 .. -0.3 % -- so it follows the
role and the tree both (the address layout each process gets with randomisation off: object
hashes, allocation sizes), and moves when either tree's contents change.

So a count difference below ~2 % between two trees is not attributable without an A/A taken in
the same configuration.  Consequences:

* 100 (R2b) said "counted MORE: +1.24 %, +1.0 % a step".  Against this offset that is NEUTRAL
  (~0, within the band).  R2b stays reverted -- no gain was measured, and it adds a context to
  keep exact -- but it was not shown to cost.  Its design still goes to stage 4's entry.
* 101 (R2c batch 1, -3.49 % spice_4049, -2.6 / -2.7 % a step) was counted in the same default
  configuration, whose A/A offset is positive: its true effect is at least that (about -5 % if
  the offset then was today's).  It stands.
* Every later count of the plan names its A/A, taken in the same configuration beside it.

## 2026-10-10 — 103. Stage R2c, batch 2: the limiter test without its closure, one capability check (-0.6 .. -1.0 % against an A/A)

* `has_limiter = tr._newton_limiter() is not None` built the limiting closure on every call to
  test it.  The chain check has already held `_newton_limiter` to be the genuine method, which
  answers None exactly where the transient's `_has_limiters` memo is falsy: read the memo (the
  method only where it is unset).
* A pre-check before the record looked the walk up (`_walk_for`) on every call to decline a
  limiting element with no C limiter -- which `_walk_ready` checks itself on its slow path, with
  the same kept decline (`walk_ok` False -> `_keep('limit')`); while the walk's stamp stands the
  answer cannot change.  Removed; the decline keys per call are the same (a circuit with the walk
  switched off declines 'newton_c:limit' every call either way).

Counted WITH ITS A/A (the change stashed, the same configuration, beside it -- roadmap 102):
change vs A/A, spice_4049 +0.705 % vs +1.689 % -> -0.98 %; mos1's step +0.11 vs +0.92 -> -0.81 %;
stage's step +0.53 vs +1.14 -> -0.61 %; mos1_adaptive +0.952 vs +1.354 -> -0.40 %.  ~7 k a call
on spice_4049, as sized; R2c now ~31 k of R2's ~55 k line.  The Newton's pinned pair moved on
its twin, re-recorded.  Full suite 5068 passed and that pair.

## 2026-10-10 — 104. Stage R2 closed: what stands, what was reverted for speed, the couplings written down

* ONE COMPANION SETUP AND WRITE-BACK (R2c's clarity item), built and COUNTED: the step's
  companion set up (the step, the one before it, the integrator after an order drop, its
  formula, history rows, coefficients) is the same fifteen lines in `_tran_core.evaluate` and
  `_tran_newton_c.solve`, and the six-attribute write-back the same in both (the stage family's
  differs in kind: no active integrator, `Geq = a h G`).  One shared function each (an 11-tuple
  out) counted against an A/A: spice_4049 +0.44 %, mos1's step +0.43 %, stage's step +0.29 %,
  mos1_adaptive +0.15 % -- ~3 k a call for two calls, a tuple and a test.  Speed is ranked first:
  REVERTED; each copy now names its twin in a comment ("a change here is a change there"), and
  the pinned pairs hold each to its reference.
* R2d, THE COUPLINGS A TWIN MUST COPY, written into `_SteppingLoop`'s docstring: who decided a
  step's size (`clamped` / `imposed` / `hold`, each flag's phase), `landing` and `order_drop` as
  one signal, the order drop's path (`where` -> `_is_first_step` -> `check_order_drop` ->
  `_effective_method` -> `order_drops`), the veto after an accept the controller has recorded,
  the retry's "same point" by object identity, `X` both output and history.  Two plan items not
  built, by reason: the "same point" as an INDEX belongs to the C loop (stages 4-5) -- the
  periodic shift edits `X`'s last entries in place, so a plain count of accepted points is not
  the same test, and the index must move there too; `_is_first_step` NOT renamed -- private, its
  meaning stated in `_RUN_STATE` and in `where`, and the alias the plan asked for would be a
  property on every read.

R2 against its line (re-sized before R2b: the R2-addressable per-call Python on spice_4049,
>= 55 k a call): ~31 k (batch 1 ~24.5 k, batch 2 ~7 k) -- NOT MET.  By the plan's own rule
what paid stays (both batches counted FEWER, beyond their A/A), what did not is reverted (R2b,
neutral; the shared companion, +0.4 %), and R2a's manifest and R2d's couplings stay (no cost).
The per-call Python left on spice_4049 is dominated by the branch screen (88 k, stage 3) and
the source pass (49 k, stage 1) -- the C loop's own stages.  Full suite 5069 passed.

## 2026-10-10 — 105. Stage 1: the sources and their breakpoints in C (`csrc/sources.c`, `_tran_sources.Table`)

Sized first: on spice_4049 the source pass is 48.5 k a call, served by the direct pass
(`_stamp_plan.source_direct`) -- its validity check 22.2 k (what no watcher sees: the elements'
dicts and types, the classes' `u` and `function`, `ParameterDict`'s lookup), the evaluation
16.8 k.  So the C sources save little per PYTHON step; their place is the C loop (stage 4), where
they run every step and their validity once a loop call.  A prerequisite, as the plan says.

`csrc/sources.c` (+ `sources.h`): SEVEN waveforms -- the plan said five; DC (the base
`TimeFunction`, `f` the int 0), Sin, Pulse, PWL, Exp, AM, SFFM -- in the Python's own operations
and order: Python's float `%` and `//` (CPython's `float_rem` and `_float_div_mod`) for the
pulse's period, the `where` chain's last true branch with each branch's arithmetic as written,
`_sin_scalar` and its guards, Exp's three branches, the AM / SFFM products as Python groups
them; the direct pass's per-source `v = a + f(t)` and its scatter in element order from +0.0
(a voltage source's `-v` into its branch row, a current source's `+v`, `-v`); each
`next_event` (the pulse's candidates from the period index, strictly after `t`).  DECLINE,
NEVER DIFFER: anything the Python makes by another route bails -- a non-finite intermediate,
a value past the fast add's 1e300 guard, a sine past `_sin_scalar`'s guards, a PWL time that
finds no segment (the Python's None raises).

`_tran_sources.table_for(sd)`: the table of a direct pass, its sources in the pass's own order
and rows, built only where every value comes the fast float way (the waveform exactly one of
the seven, the numeric toolkit's `exp` / `sin` / `where` / `pi`, every parameter a float or a
small int, the DC value and the waveform not both ints -- the Python adds those as ints, an
exponential's taus nonzero -- the Python raises -- numpy's exp / sin the C library's on this CPU
by the HDL backend's probe, `libm_check`).  Per call: the toolkit's functions, the sine's
switch, numpy's underflow state.  `Table.holds()`: every input as read (each source's
`function`, its class, its attributes, a PWL's points, the DC value).

Tests (`test_tran_sources.py`): every waveform V and I on a 3001-point grid and its corners
(+-1 ulp): `u` byte for byte and the next event, none declined; a hypothesis twin per waveform
(signed zeros, subnormals, 1e299 .. 1e301, infinities, NaN, ints among the values; times
negative and past the period) -- wherever the C answers it is the Python's (185 vectors
compared, 53 bails, 63 declines at build in the gate profile); each decline; a table whose input
changed does not hold; a parameter a ulp off is caught.  The unit in `test_c_abi` (complete,
`-Werror`, its header its cdef); a pinned pair "the sources and their C" (the waveforms' `f` /
`next_event`, `_sin_scalar`, `_SourceDirect.run` against the C and its table).

NOT wired into the Python step's source pass: the benchmark decks hold one DC source each
(spice_4049, spice_ring), whose pass is plumbing, not waveforms -- wired only if it counts FEWER
against an A/A, which a single DC source will not show.  Stage 4 calls the table.
Full suite 5084 passed; the sanitizer's fast tier clean (4816 passed).

## 2026-10-10 — 106. Stage 2, step a: every sparse run factors on the circuit's STRUCTURAL pattern with KLU (Andreas's decision; one re-baseline)

Measured first.  The three sparse gold decks (rich3 269 unknowns, arom 295, gm19 429) ran SuperLU
in Python (`AutoSolver`: SuperLU below 1000, the C Newton declining `newton_c:linsolver`); a
profile of rich3 found per step a dense n x n assembly, the error test's dense solve, SuperLU
re-ordering at every factorisation and `nonzero` building the CSC from the dense J.  A
throwaway prototype -- KLU on a FIXED pattern, the values taken with one `take` -- first replaced
`solve` only and barely ran (the sparse Newton factors once a solve and reuses the factors:
`factor`), then `factor` too: rich3 -25 %, arom -25 %, gm19 -5 % at 15 / 15 / 5 % of the run.

Built: `analysis.jacobian_structure(cir, iref)` -- every position an element can stamp (the full
block over its node map, G's and C's alike) and the whole diagonal (gmin, shunt, pseudo-transient
steps), the reference row and column out and renumbered as `_reduce_ndarray`, CSC order;
`linearsolver.StructuralKLU` -- ONE analysis a run, each call's values at the pattern's positions,
`solve` refactoring with `KLUSolver`'s residual validation (a fresh factor where it fails, the old
numeric freed), `factor` a numeric of its own with the transposed solve; a matrix the pattern does
not cover (another shape or dtype, a nonzero off it: a monodromy, a block system) goes to
`fallback()` -- today's choice for its size -- and leaves the factorisation alone;
`AutoSolver(structure=...)` takes it for a sparse matrix exactly where an analysis hands a
structure (`Analysis._get_linearsolver`; the analysis held weakly -- a closure over it made a
cycle, and a solved PSS waited for the collector: `test_a_solved_pss_is_freed...` caught it).
Without a structure (a solver built directly), today's choice, unchanged.

THE RE-BASELINE.  In the suite, one test runs a sparse analysis
(`test_linearsolver.py::test_the_analyses_default_to_one_autosolver`, 262 unknowns: now on
`StructuralKLU`, passing); every other run is dense and unchanged (a census plugin over the whole
suite: 9 tests made a sparse choice -- that one, three that build `AutoSolver` directly, keeping
today's choice, and the new tests).  The sparse runs are the benchmark decks, full length, both
versions the same day (`benchmarks/spice_suite.py`): rich3 15.02 -> 11.32 s (-25 %), arom 14.77 ->
11.11 s (-25 %), gm19 44.54 -> 37.12 s (-17 %); steps 3639 -> 3620, 3151 -> 3129, 3901 -> 3882;
rich3's gold metric 16.813 both (arom and gm19 have no gold file).  The answers move by more than
their last bits: a last bit moves a step decision, and the run takes another valid grid -- on the
old grid, rich3's largest |new - old| is 3.2e-2 V at a fast edge (on a 4.25 V swing), 6.5e-4 at
tend.  THAT IT IS THE INTEGRATION ERROR, measured: at reltol 1e-4 (the deck's) max |new - old|
3.23e-2, median over unknowns 2.0e-3; at 3e-5, 7.2e-3 and 9.6e-4 -- it shrinks with the
tolerance.  (At reltol 1e-3 the OLD SuperLU path raises `SingularMatrix` on rich3; the structural
KLU runs.  At 1e-5 both fail the same step.)

Tests (`test_structural_klu.py`): the structure covers every nonzero of G, C and G + a C at drawn
states of a 300-section ladder, in CSC order without duplicates; the solver equals the dense
solve (one and several right-hand sides, complex, the factor's solve and transposed solve); a
dense matrix goes to the fallback and leaves the numeric alone; a refactor failing its residual
is factored again; `AutoSolver` takes it exactly where a structure is given; a sparse transient
runs on it and agrees with the dense solver.  Two pinned pairs moved on their REFERENCE
(`Analysis._get_linearsolver`, `AutoSolver`): the readiness stamps read the same solver identity
(one per analysis, the same attribute), and the frozen transform takes the analysis's solver for
its real block where the choice is not dense -- both twins checked, re-recorded.

Not yet: the C side (the core assembling on the pattern, a C linsys with KLU, the C Newton serving
sparse runs) -- stage 2's next steps; `_tran_radau_tc`'s complex KLU still takes its pattern from
the values.  Found on the way: the editable install's import finder serves THIS tree's package
whatever `sys.path` says, so a parent worktree cannot be A/B'd by path -- stash in the same tree.

## 2026-10-10 — 107. Stage 2, step b: a sparse run's error test solves on its structural factorisation (-15 .. -23 % more)

After step a, rich3's profile put the error test's dense LU first among what was left: the C
error test (`_tran_lte_c`) maps the charge LTE by `J^-1` with numpy's own `dgesv` on the m x m
reduced Jacobian, 0.52 ms a step at m = 268 against KLU's 63 us factor of the same matrix.  The
reference (`StepController._charge_lte`) solved with the toolkit's dense solve.  For a sparse run
-- the analysis's `AutoSolver` having chosen `StructuralKLU` -- the transient now hands the step
controller that solver (`StepLTEInputs.linsolver`, from `_lte_inputs`' `_sparse_linsolver()`;
None everywhere else) and `_charge_lte` solves `J^-1 Eg` with it, as the Newton's solves are;
the C error test declines there (`lte_c:sparse`: its `dgesv` is no longer the reference).  Dense
runs: unchanged (None: the toolkit's solve, the C serving).  The same decision as step a (a
sparse run's linear algebra on the structural pattern); its answers move once more.

Measured, full length, both versions: rich3 11.71 -> 9.48 s, arom 11.45 -> 9.79 s, gm19 39.49 ->
30.49 s; with step a, against SuperLU and the dense error test: rich3 15.0 -> 9.5 s (-37 %),
arom 14.8 -> 9.8 s (-34 %), gm19 44.5 -> 30.5 s (-31 %).  The answers move at the integration
error's level as in step a (on the old grid max |new - old| 3.4e-2 / 2.5e-2 / 4.5e-2, at tend
3.7e-3 / 3.3e-3 / 7.3e-4).  `test_structural_klu.py` gains: a sparse run's error test goes to its
`StructuralKLU` (the C declining `lte_c:sparse`), a dense run's does not.  `StepLTEInputs`' field
list grew by `linsolver` (its contract test's input set with it).  The pinned pair "the error
test and its C" moved on both sides -- the reference's sparse branch and the twin's decline for
exactly it -- checked, re-recorded.  Full suite 5090 passed.

## 2026-10-10 — 108. Stage 2, step c: the evaluation core writes only the circuit's pattern (byte for byte the same)

After steps a and b, the core's `evaluate` was the largest piece of a sparse deck's step (rich3:
0.66 s of 2.37 s profiled, ~0.2 ms a call): each call zeroed and filled dense n x n buffers -- the
G and C bins, the companion `Geq`, the Jacobian `J` -- about 2 MB of memory traffic at n = 269,
where the circuit's pattern holds 1816 entries.  Now a caller whose matrices are zero off the
pattern asks the core to write only the pattern's entries (`PCX_WANT_PATTERN` in pcx_abi.h; the
core's `npat` / `pat`: every entry a matrix pass can write, and the diagonal): the bins zeroed
there, `Geq` and `J` formed there.  The same arithmetic on every entry written; the entries left
are the +0.0 the full loops write there (`a0 * 0.0`, `0.0 / h`, `0.0 + 0.0`: every formula's
coefficient is positive).  Opt-in per call: `evaluate` and `passes` (their fresh matrices
`np.zeros` where the core has a pattern) and the C Newton (its buffers zeroed once, `s.pattern`);
the Radau units keep every entry.  Used where the pattern is at most a quarter of n x n
(`PATTERN_MAX_FILL`; env `PYCIRCUIT_CORE_PATTERN=0` turns it off): mos1 (24.8 %), gp, psp and the
sparse decks; not the 7-unknown stage.

Proved: `test_core_pattern.py` -- a 30-stage chain, fixed gear (the C Newton) and adaptive
trapezoid, and a 250-stage chain (a sparse run: `evaluate` and KLU), each the same bytes and
statistics with the pattern on and off; a pattern missing the first drain's diagonal entry makes
another run (a singular J).  The sparse decks full length, on and off: byte for byte the same,
rich3 9.71 -> 9.39 s, arom 9.16 -> 8.86 s, gm19 28.95 -> 27.99 s (-3.3 % each; the fresh zeroed
matrices a call keeps give back part of the saving).  Counted against HEAD with an A/A beside
it: gp -2.9 % a step; mos1, psp, mos1_adaptive -0.2 .. -0.6 %; stage and spice_4049 (no pattern)
+0.3 / +0.5 % first, then recounted twice: +0.01 .. +0.26 % against an A/A of +0.25 .. +0.51 % --
noise, not a cost.  Ten pinned pairs moved, nine on the twin (the core's and the Newton's C, and
every unit including pcx_abi.h, whose enum gained the bit); "the step's end passes and the C's"
also on its reference (`passes`' zeroed C and the bit) -- its twin (radau_tc's end passes)
takes the dense route to the same bytes -- checked, re-recorded.  Full suite 5095 passed; the sanitizer's fast
tier clean (4824 passed).

## 2026-10-10 — 109. Stage 2, step d: the C Newton serves sparse runs -- the chord on KLU from C, byte for byte the Python's

Sized first (rich3, a quarter of the run): the Python Newton was 1.78 s of 2.60 s, of which the C
core 0.46 s and KLU 0.10 s -- 1.22 s (47 % of the run) of Python around them, ~320 us an
iteration against the core's 121 us.  The C Newton declined every sparse run
(`newton_c:linsolver`: the AutoSolver's choice not the dense one).  Now, where the choice is the
run's `StructuralKLU` and the Newton is the chord, `newton.c` factors and solves with KLU itself:
at the first iteration the pattern's values gathered straight from the core's full J (`kat`: each
CSC entry's place there -- no reduced copy, no coverage count) and `klu_factor` on StructuralKLU's
own analysis and common -- the numeric its `factor` makes; each iteration `klu_solve` on -F, as
`_KLUNumeric.solve`; |J| written on the pattern only (`kabs`; zero off it, so the convergence
test's `dgemv` reads the matrix numpy's `dot` reads); the numeric freed on every return
(`hdl_fn` wraps `iterate`).  `_tran_newton_c._klu_bind` points the C at the solver once a run,
declining where it is not analysed yet or where the core's pattern does not lie within KLU's (a
nonzero there would go to its fallback in the Python); a full Newton on KLU declines
(`newton_c:sparse_full`: its solve refactors the run's numeric, validated by the residual).  A
served call counts the numeric (`numerics`), as the Python's `factor` would.

Proved: `test_newton_c.py` -- a 250-stage MOS chain with capacitors, gear / trapezoid / Euler
adaptive and gear fixed: the same bytes, steps, iterations and StructuralKLU counts with the C on
and off; a bail on every call leaves no trace; a core pattern given an entry KLU's lacks declines
(`newton_c:klu`) with the Python's answer; alive: the first drain's diagonal read from the
second's place makes another run.  The sparse decks full length, C Newton on and off (off = the
path before): byte for byte the same; rich3 9.16 -> 6.52 s (-29 %), arom 8.66 -> 6.24 s (-28 %),
gm19 27.54 -> 24.00 s (-13 %; every Newton call served -- its remaining cost is elsewhere, below).
Counted on the dense cases against HEAD with an A/A: no change on every case.  Ten pinned pairs'
twins moved (pcx_abi.h gained `klu_factor_t` / `klu_free_numeric_t`; newton.c/h): checked; "the
Newton solve and its C" now names `StructuralKLU.factor`, `._factor` and `_KLUNumeric._run` as
reference and `_klu_bind` as twin -- re-recorded.

Found while sizing gm19 (not changed): the DC operating point solves DENSE whatever its size --
`DCAnalysis`'s two `solve_system` calls pass no `linsolver`, so every rung of the gmin / gshunt
continuation takes `toolkit.linearsolver` (numpy's dgesv): on gm19 (429 unknowns) 3227 dense
solves at a tenth of the run.  Routing it to the analysis's AutoSolver would move the sparse
circuits' operating points in their last bits: a decision, put to Andreas.

Full suite 5102 passed; the sanitizer's fast tier clean (4734 passed).

## 2026-10-10 — 110. The operating point solves on the analysis's linear solver (Andreas's decision) -- gm19's DC op on structural KLU

Found sizing step d (109): `DC`'s two `solve_system` calls passed no `linsolver`, so the operating
point solved dense at any size -- every rung of the gmin / gshunt / source-stepping chain took
`toolkit.linearsolver` (numpy's dgesv) -- while the transient after it ran on the circuit's
structural KLU; a user's `linearsolver` never reached a DC either.  Now both calls pass
`self._get_linearsolver()`: below 250 unknowns the AutoSolver's choice is `DenseSolver`, whose
`solve` IS `toolkit.linearsolver` -- byte for byte the call before (the default `StandardNewton`
and the continuation wrappers solve through `linsolver.solve`); from 250 the circuit's
`StructuralKLU` (or today's value-patterned sparse choice where the DC's system is reduced
further: rich3 and arom hold nodes -- SuperLU, their DC 0.02 s either way).

Measured (a scratch wrapper on `DC._chain` first, then the change): gm19 (429 unknowns) DC op
17.8 -> 6.9 s on a loaded box (3227 dense solves at a tenth of the run were the DC's); its
transient answer moves by at most 1.3e-23 absolute (max |x| 10.2); rich3 and arom identical.
Tests (`test_structural_klu.py`): a 30-section diode ladder's operating point is byte for byte
the call with no solver; a 300-section one factors on the structure (nothing foreign) and agrees
with the dense solve to 1e-12.  Full suite 5104 passed, nothing else moved.

## 2026-10-10 — 111. Stage 2 closed: the reduced system's dense plumbing in one header (`csrc/pcx_dense.h`)

The reference row out and back and LAPACK's column-major copy were written out by hand at each
site: seven loops in `newton.c` (the seed reduced; each evaluation's x, the walk's two vectors
and the converged point put back; F reduced; the walk's result reduced; J copied with |J|) and
three in `lte.c` (J copied with its finiteness bail, the right-hand side reduced, the result put
back).  Now four inline helpers, named for the Python they stand for: `pcx_drop_ref`
(`remove_row_col` on a vector), `pcx_put_ref` (`insert_row`: +0.0 at the reference),
`pcx_reduce_cm` and `pcx_reduce_cm_abs` (the reduced J column-major -- numpy's / SciPy's copy for
LAPACK -- and |J| row-major for the convergence test; each returns whether every entry copied was
finite).  The Newton's chord reads that flag instead of scanning A again; the error test bails on
it as before.  Radau's 3 x 3 block layouts keep their own loops (another layout, not this one).

Not built, by choice: one polymorphic "linsys" over LAPACK and KLU.  The bits depend on which
library a solve calls (numpy's ILP64 dgesv for the full Newton and the error test, SciPy's LP64
getrf/getrs for the chord, KLU's factor/solve for the sparse chord), and a struct that hides the
choice hides the one thing a reader of a twin must see; the three solves stay where they are
called, each a few lines, each naming the Python it stands for.

Proved: the Newton and error-test byte-identity suites (test_newton_c, test_lte_c) unchanged and
passing; `-Wall -Wextra -Werror` clean; the header added to the `#line` plant test (an error in it
reported at pcx_dense.h:line).  Two pinned pairs moved on the twin only ("the Newton solve and
its C", "the error test and its C"): checked, re-recorded.  Counted against HEAD with an A/A beside it:
no change on every case (-0.34 .. +0.30 %, the A/A +0.04 .. +0.57 %).  Full suite 5105 passed; the sanitizer's
fast tier clean (4790 passed).

## 2026-10-11 — 112. Stage 3, step a: Euler's error estimate in C -- the step after a breakpoint no longer hands the error test back

Sized first (a census of 26 gold decks at a fifth of their runs): the C error test served every
gear step and declined two kinds -- every sparse run (`lte_c:sparse`) and EULER, the step after
each breakpoint where the order drops (`lte_c:integrator`: jge 1037 of 3198 error tests, mux8
244, 1..47 on the others; confirmed Euler by type, every one).  Now `lte.c` has Euler's estimate
(`LTE_FORMULA_EULER`: `coef * ((q - q1) / h - iq1)`, `iq1` the last companion current, checked
finite), the glue its coefficient in Python floats -- `-0.5 * (2 h / (h + h_last) if h_last else
1.0)`, `EulerIntegrator.compute_lte`'s variable-step correction -- and its order (2.0); the J
solve, the reference and the tolerance are the gear and trapezoid's.  A run's first step (no
history: zeros) still declines (`lte_c:formula`, 1-2 a run).

Proved: test_lte_c -- the pulse's order drops now served (`lte_c:euler` > 0, no `integrator`
decline) with every attempt's verdict, error and running reference the Python's; the drawn
bit-for-bit twin draws Euler too (its own companion currents, a last step missing or zero): 3000
examples, 759 Euler error tests served, every one the chain's bits.  jge and mux8 at half
length: byte for byte the same.  Speed: neutral, as the stage's line asks -- jge 23.0 s both
ways (instructions +0.41 %: the C's reduced, column-major copy of J -- an element-wise strided
transpose -- costs ~1 M instructions at n = 244 where numpy's copies cost less, and the gear steps
already pay it: a known cost, for the loop in C to remove); mux8 instructions the same; the
step_machinery cases (no Euler step among them) no change against an A/A.  "The error test and
its C" now names `EulerIntegrator.compute_lte` among its references: re-recorded.  Full suite 5105 passed; the
sanitizer's fast tier clean (4798 passed).

## 2026-10-11 — 113. Stage 3, step b: a sparse run's error test in C -- `StructuralKLU.solve` twinned; rich3 -4.8 %, arom -6.1 %

The census (112) found every sparse run's error test in Python (`lte_c:sparse`: 15 % of rich3
and arom, 8 % of gm19).  Now `lte.c`'s `sparse_solve` is `StructuralKLU.solve` on the reduced J:
the reduced J row-major (`pcx_reduce_rm`, remove_row_col's copy), the pattern's values by
`take`, `_ours`'s check (as many nonzeros as the matrix -- `pcx_count_nonzero`; another matrix
is handed back, its fallback solver is the Python's), `klu_refactor` on the run's numeric,
`klu_solve`, and the refactor's validation (max|R x - b| / max(max|b|, 1e-300) against the
solver's `REFACTOR_RESIDUAL_TOL`, R x numpy's `dgemv` -- `arr @ x` is that call); a failed
validation is handed back too, so the Python alone ever replaces the numeric, and a run's first
error test (no numeric yet) declines (`lte_c:sparse_first`).  A served call counts the refactor
(`refactors`), as `solve` would.  The glue binds once a run (`_klu_bind`: the solver's CSC
pattern, `take`, analysis, common, the library's entries, numpy's `dgemv` from the Newton's
driver).  And J is read where it stands when C-contiguous (the C only reads it): an n x n copy
a call saved, on the dense path too.

First built with the C's own element-by-element loops (a finiteness test per entry in the
reduce, a plain nonzero count): gm19 (m = 428) read +3.5 % instructions, +1.2 % wall -- numpy's
vectorised copies and counts beat them at that size.  Then: the reduce is two `memcpy` a row and
tests nothing (a non-finite entry off the pattern is a nonzero the count hands back; one on it
makes the residual non-finite, handed back too), the count runs in blocks of eight, and J is not
copied: gm19 -0.28 % instructions.

Proved: test_lte_c -- a 250-stage MOS chain, gear and trapezoid: every attempt's verdict, error
and running reference the Python's; a residual bound below zero hands every solve back with the
Python's run; a drawn twin at the `StepLTEInputs` level -- banded patterns, values, a nonzero
off the pattern now and then, a failing residual bound, gear / trapezoid / Euler -- compares the
bits, the running reference, the warnings and the solver's own counts, each side on a solver
primed alike: 3000 examples, 1601 sparse solves served, 797 foreign and 345 failed-residual hand-
backs.  The decks full length against HEAD, twice each: byte for byte; rich3 6.50 -> 6.19 s
(-4.8 %), arom 6.24 -> 5.86 s (-6.1 %), gm19 no change (17.12 .. 17.73 s both ways).  The
step_machinery cases no change against an A/A.  "The error test and its C" names
`AutoSolver.solve`, `StructuralKLU.solve` and `._ours` among its references and `_klu_bind`
among its twin; "the Newton solve and its C" moved on its twin (pcx_dense.h): re-recorded.  Full suite 5109 passed; the
sanitizer's fast tier clean (4801 passed).

## 2026-10-11 — 114. Stage 3, step c: the C reads the Newton's and the error test's n x n inputs where they stand -- gm19 -4.9 %, rich3 -2.4 %

What 3c named -- the error test reading the Newton's buffers, no n x n copy of J -- landed with
3b for J (read in place).  The C Newton made the same copy the other way: every call copied the
cached `C` at the seed (`_C_cache`, which the C only reads -- the first evaluation's companion,
the chord's later ones) into its own buffer, 580 KB a call at n = 269, 1.5 MB at n = 429.  Now
it is read where it stands.  Measured both ways first: in place everywhere read +0.1 .. +0.25
points of a step on the small cases (n <= 25: the cffi buffer object and re-pointing the struct
cost more than copying a few hundred doubles), neutral at n = 45 -- so both inputs (this `C` and
the error test's J) go in place from `_cbuild.INPLACE_MIN_ENTRIES` (4096 entries, n >= 64) and
are copied below it.

Measured: rich3 6.14 -> 5.99 s (-2.4 %), gm19 17.54 -> 16.68 s (-4.9 %) full length, byte for
byte; the step_machinery cases no change against an A/A (within 0.12 points).  The Newton and
error-test byte-identity suites pass both routes (the 250-stage chains in place, the small
circuits copied).  Two pinned pairs moved on the twin only (`_tran_newton_c`'s `solve` and
`_Ctx`, `_tran_lte_c`'s): re-recorded.  Full suite 5109
passed (no C source changed).

## 2026-10-11 — 115. Stage 3, step d: the branch check's cheap screen in C -- the SVD stays Python's; rich3 -3.3 %, stage 3 closed

Sized first: after every converged step the branch check screens `C` at the converged point
(`_branch_screen_C`): |C|'s largest entry, the run's running maximum, the diagonal's largest and
its smallest positive entry -- six numpy reductions over n x n, 6 us a call at n = 15 (4.4 % of
4049osc), 27 us at 171, 60 us at 269 (4.7 % of rich3).  (Found on the way, not changed: on
ab_integ the screen fires on EVERY step and its SVD clears it -- 823 us a step, 60 % of the run;
e1480 and g1310 the same: the plan's "after this plan" item, now the largest single cost left on
those decks.)  Now `csrc/branch.c` makes the cheap level -- the same comparisons on the same
values, Python's `max(cmax, nrm)` included -- and hands back QUIET, COLLAPSED or SVD; the SVD is
the Python's, split out as `_BranchCheck._branch_svd` (one copy for both routes); a non-finite
entry hands the screen back.  THE PATTERN: where `C` is one the core wrote on its pattern -- the
C Newton's converged point, `_tran_core.evaluate`'s fresh C -- the writer records it with the
pattern (`_C_pattern`, declared in `_RUN_STATE`), and the screen reads only those entries of
that object (an identity test: any other C is read whole up to 4096 entries, numpy's above).
Declines on a patched or shadowed `_branch_screen_C` (a test counts its calls).

Proved: test_branch_c -- a drawn twin (random patterns with the diagonal, entries at the
thresholds -- 1e-9, 1e-6, subnormals, signed zeros, 1e300 -- and a drawn running maximum)
compares verdict, largest entry and running maximum bit for bit by both routes: 3000 examples,
quiet 450, collapsed 1060, SVD 4490; non-finite entries hand back; the declines; alive: a pattern
missing the largest entry moves the largest entry; four whole runs (gear, trapezoid, fixed, a
250-stage sparse run on the pattern route) the same bytes, statistics, warnings and running
maximum.  Decks against HEAD, twice each: byte for byte; rich3 5.98 -> 5.79 s (-3.3 %), mux8
(half) 14.33 -> 14.18 s (-1.1 %), 4049osc ~ -2 % (noisy).  Counted against an A/A: stage
-1.55 points, spice_4049 -2.5, the rest no change.  New pinned pair "the branch screen and its
C"; "the transient evaluation and the core" and "the Newton solve and its C" moved on the twin
(the marker's writes): re-recorded.  Full suite 5119 passed; the sanitizer's fast tier clean (4830
passed).

## 2026-10-11 — 116. Stage 4 sized, its line re-sized for psp; step a: the stage predictor's fast path in C (not yet wired)

SIZED (scratchpad s4/size.py, size4049.py: the marginal fixed step at 100 vs 300 steps, median of
5, the C entries timed inside it): stage 157 us a step of which 54 in C (ceiling 2.9x), mos1 107 /
20.5 (5.2x), gp 131 / 38 (3.4x), psp 901 / 761 (1.18x), spice_4049 91 / 7.8 (11.7x).  mos1's
Python by the profile: the C Newton's glue ~40 %, the predictor ~15 %, the accept and the rings
~15 %, the source pass ~12 %, the breakpoint poll ~8 %.  Inside a C loop no Python runs between
steps, so every check the glue makes per call is made once at the loop's entry.  THE LINE,
RE-SIZED BEFORE BUILDING: stage, mos1, gp, spice_4049 >= 1.3x each (unchanged); psp >= 1.05x (its
ceiling is 1.18x: 761 of its 901 us are PSP's evaluation).  Steps: 4a the predictor in C, 4b the
loop, 4c the write-back proven complete, 4d the line.

STEP 4a: `csrc/pred.c` / `_tran_pred_c` -- `_StagePredictor._predict_fast` with
`_pred_fast_nodes` and `_fit_fast`: the nodes newest first (the 1e-13 dedupe, the smallest gap,
the tie test), the fit's times in Python's float order, `np.vander`'s running products, numpy's
own `dgesv` for the weights (kept for the last `tau`, as the memo keeps them) and its `dgemv` for
`w @ X` -- `cblas_dgemv(RowMajor, Trans, ...)` IS that product (3000 of 3000 drawn; a plain loop
1608: the kernel's own order) -- and the clamp in numpy's `clip` order (`m = x > lo ? x : lo; m <
hi ? m : hi`: every chunk length, signed zeros included, matches the ufunc's vector path).  A
flag numpy would report, and every case the fast path sends to the general one, hand back.  Not
yet wired (as stage 1's sources): the loop (4b) calls it.  `pcx_abi.h` gains `PCX_CBLAS_TRANS`
(every unit's twin moved with it).

Proved: test_pred_c -- a drawn twin (histories uniform, graded, duplicated inside the dedupe,
tied, out of order, at or past the target; states of every magnitude; degrees 1-4): 3000 drawn,
437 served and 525 "no prediction" byte for byte the Python's, 1761 to the general path both
ways, 277 handed back on a flag (extreme magnitudes); the general path's cases; a linear ramp
predicted exactly; alive: the clamp's factor a ulp off moves a prediction the clamp bites.  New
pinned pair "the predictor and its C"; every unit's twin re-recorded for pcx_abi.h.  Full suite 5125 passed; the
sanitizer's fast tier clean (4817 passed).

## 2026-10-11 — 117. Stage 4, steps b and c: the fixed-step multistep loop in C, its write-back proven complete

`csrc/fixed.c` / `_tran_fixed_c` run `_SteppingLoop`'s fixed multistep steps themselves --
where (the grid; the next breakpoint as `_next_breakpoint` makes it: numpy's maximum of t and
the source table's earliest event, past one closer than `minbreak`; the order drop a landing
asks for; the grid's end), take (the predictor in C, 4a; `check_order_drop` and `bdf2_alphas`
in the glue's float order; u(t) from the source table, stage 1; the C Newton through its own
struct, set up as its glue sets it up, its converged C and Geq alternating between two buffers
so a step handed back after its Newton leaves the last accepted step's; the branch screen in C,
3d) and accept (the rings, `_pred_note`'s ordered, deduplicated nodes, the outputs, the
statistics) -- driving the units through their entries and structs, as the Newton drives the
core.  It stops BEFORE a step it cannot take as the Python would -- a Newton that fails, a
screen that fires (the SVD is the Python's), a predictor or a table that hands back -- and the
Python loop takes that step; the next enters the C again.

THE RUN CONTEXT: inside the C loop no Python runs between steps, so the C Newton glue's
per-call checks hold for every step of an entry when its last call SERVED and the watch counter
has not moved since (`_tran_newton_c` marks `served` with `_watch.now()`, cleared on any other
outcome; watching off, the loop does not run).  Checked once an entry: every piece the steps
call is its module's own and not shadowed (the loop's, the predictor's, the history's, the
screen's, `_source_at`, `SubCircuit.u` / `next_event` / `accept_step`, the integrators'
`check_order_drop`, the family's, `_tran_newton_c.solve`), the circuit's methods not shadowed
on the instance (a source perturbed there, as the shooting tests make one -- found by
test_source_memo, which the first version failed), no `accept_step` hook, every breakpoint
the table's (the sources' `next_event`, HDL devices' inf), the table holding, no event rows or
periodic states or stateful limiters, the predictor's settings, the cached C at the last point.
State copied in (the rings, the predictor's nodes, the cached C, the step record) and written
back as the Python steps leave it (`X` and `timelist` a fresh array and float a step, the
statistics, the rings, `_q_cache` / `_C_cache` / `_C_pattern`, `_iq`, `_Geq`, `_Cmat`, the
active integrator and its coefficients, the predictor's history, the step record, the screen's
maximum, the StructuralKLU's `numerics`) with the `_paths` tallies the Python steps make, under
their keys.

PROVED (test_fixed_c): the run with the C loop on is the run with it off -- every output byte,
time, statistic, warning, the `_paths` tallies (the loop's own build counts aside) and every
`_RUN_STATE` attribute (the predictor's two memos, the core's check counters and an id-keyed
cache aside: 4c, 49 of 53 compared byte for byte), and the linear solver's counts -- on MOS
level 1, Gummel-Poon on the chord, the PSP stage behind the screen, a 250-stage sparse run, a
pulse's landings and order drops, the trapezoid, Euler, a grid that does not divide `tend`,
entries that fill their outputs (seven rows); through its hand-backs (a screen made to fire, a
Newton held to two iterations); its declines (an `accept_step` hook, a patched piece, a
shadowed method, a perturbed source); alive: the predictor's clamp a ulp off is another run.
The watch tests that count the Python steps' per-call checks run with the C loop off (it makes
them once an entry), and the watched/unwatched comparison sets its own counts aside as it sets
`once:` aside.  New pinned pair "the fixed-step loop and its C"; "the Newton solve and its C"
and "the stepping loop and its unread Jacobian" moved on their twins (the served mark; the
loop's call): re-recorded.

The marginal fixed step, C loop on / off (alternated, median of 7): stage 38.6 / 153.8 us
(3.98x), mos1 19.4 / 105.8 us (5.46x), gp 36.2 / 130.5 us (3.60x), psp 735.6 / 900.3 us
(1.22x) -- above the sizing's ceilings on stage (2.9x): the Newton called from C costs less
than the cffi call the sizing timed.

STEP 4d, THE LINE (`step_machinery --check` against the parent tree, 15-min load < 3, six paired
rounds, the whole timed case -- the operating point and the Python's first steps in it): stage
-60.6 % [-60.9, -59.6] (2.54x), mos1 -58.3 % [-59.8, -57.8] (2.40x), gp -43.5 % [-44.1, -43.3]
(1.77x), psp -15.6 % [-16.3, -15.3] (1.18x), spice_4049 -82.9 % [-83.4, -82.7] (5.85x); child
faster 6 of 6 on every case; bytes and statistics the same across the trees, the paths moved by
the loop's own counts only.  THE LINE HOLDS: stage, mos1, gp, spice_4049 >= 1.3x with the
interval above it; psp >= 1.05x.  Stage 4 closed.  Full suite 5140 passed; the sanitizer's fast
tier clean (4866 passed).

## 2026-10-11 — 118. Stage 5 sized; step a: the step controller's law in C (not yet wired)

SIZED (scratchpad s5/size.py, 26 gold decks full length, the C calls timed inside the run): the C
Newton serves 95-100 % of every deck's attempts; the ceiling -- the run over its operating point
plus its time in C -- median 3.08x: 3-10x on the small and medium decks, 1.25-1.6x on jge, mux8,
gm19, voter25, rich3 and arom (the core's C dominates them).  ab_integ, e1480 and g1310 will gain
little whatever the loop: their screen fires every step and the SVD is the Python's.  The line
stands (>= 1.4x median over the decks the loop covers).  What an adaptive attempt adds to a
fixed step: the judge (the C error test of stage 3, the controller law), a rejection retried at
the same point within the budget, breakpoint truncation, the cap; a rejected attempt leaves its
Newton write-backs and the error test's running reference behind, so the loop keeps them per
attempt.  Steps: 5a the law, 5b the adaptive path in the one loop unit, 5c the write-back, 5d
the line.

STEP 5a: `csrc/ctrl.c` / `_tran_ctrl_c` -- `IntegralController._evaluate` past its error test,
the default band: `0.9 ** p` (libm's pow, as Python's float `**`), the shrink with
`_retry_exponent` under every `retry_order` (the last rejection at the SAME point -- the object
the controller kept, by `is`, an index here -- larger and above the target: the order the two
tries measured, clamped, and carried), the grow, the cap, and `_judged`'s retry state; Python's
`min` and `max` as its builtins pick (the first argument unless the second is strictly beyond).
Not yet wired: the loop (5b) calls it.

Proved: test_ctrl_c -- a drawn twin against `IntegralController._evaluate` (its error test
stubbed): 3000 drawn, 1702 accepted (981 at the cap) and 1298 rejected, 127 retries that measured
their order and 237 that carried one, every verdict, next step, `_retry` and `_q_carried` bit for
bit; a second retry at one point measures its order, at another it does not.  New pinned pair
"the controller law and its C".  Full suite 5144 passed; the sanitizer's fast tier clean (4818
passed).

## 2026-10-11 — 119. A FIX to 106/110: a structural KLU refactor that meets a zero pivot factors again -- opampal's operating point had failed since the DC routing

Found sizing stage 5: opampal (519 unknowns) had dropped out of two censuses (26 of 27 decks) and
nobody asked why.  Its operating point failed: "singular Jacobian: 'xa:qna.bi' appears in no
equation".  The cause is 106's `StructuralKLU.solve`, reached by the DC since 110 (212a5b06): a
refactor reuses the first factor's pivot ORDER, and where a pivot of that order went to zero it
raised "singular" -- without the fresh factor a failed residual already gets.  An operating
point's gmin and source ladders move the matrix far between solves; opampal's met it.  Measured,
its operating point under each solver: dense 176 s, the structural KLU as it was FAILED (24.5 s),
KLU factored fresh every solve 40.6 s, SuperLU 74 s.  Now a refactor that fails is followed by a
fresh factor (`pivot_fallbacks` counts them; only the fresh factor's failure is the matrix's):
opampal's operating point converges on the structural KLU in 36 s (34 pivot fallbacks), within
2.2e-16 of the dense one (max |x| 35).  The transients never reach the branch (their factors
stay one: rich3 and gm19 byte for byte as before).  `KLUSolver` (the value-patterned one, n >=
1000 without a structure) keeps the old behaviour -- the same shape, not reached here; noted.

Proved: test_structural_klu -- a refactor on a permutation the first factor's pivots cannot
take factors again and solves it (one pivot fallback), a rank-one matrix still raises; the test
fails on the code before ("Singular matrix").  "The error test and its C" moved on its reference
(`StructuralKLU.solve`): its C hands a failed refactor back, which the Python now factors again
-- checked, re-recorded.  The lesson kept: a solver or DC change runs every gold deck, its
operating point included, before it is committed -- the suite runs none of them.  Full suite
5145 passed.

## 2026-10-11 — 120. Stage 5, steps b and c: the adaptive multistep loop in C -- one loop unit for both modes

`csrc/fixed.c` became `csrc/loop.c` (`_tran_fixed_c` -> `_tran_loop_c`, `PYCIRCUIT_LOOP_C`,
`loop_c:` counts; the tests `test_loop_c.py`): one stepping loop in C, fixed or adaptive, as the
Python has one.  An adaptive attempt adds to a fixed step: `where`'s truncation onto the next
breakpoint (a landing) and the cap; the judge -- the C error test (stage 3) driven through its
own struct (the estimate's coefficient and step sums as `_tran_lte_c` makes them in Python
floats; its sparse solve bound to the run's StructuralKLU as its glue binds it; an `lte_algebraic`
mask copied in once a run), then the controller's law (5a), the "same point" the index of the
accepted point the attempt left from; a rejection retried at the point within the family's
budget (`max_retries`, `max_reject`, `minstep`) -- the attempt's Newton write-backs, the screen's
maximum, the error test's running reference, the law's retry record and `last_err` kept PER
ATTEMPT (the converged C and Geq alternating between two buffers per attempt, so an attempt
handed back leaves the last one's); a force-accept handed back (its warning's record is the
Python's), the law's record restored first; the accept as the fixed loop's, the next step
`min(h_next, max_step)`.  A rejected attempt keeps `where`'s landing: its retry drops the order,
as the Python's `where` reads the landing flag the rejected attempt left.  THE ENTRY: the error
test's glue marks its last serve as the Newton's does (`_tran_lte_c`'s `served`), and an adaptive
entry needs both; checked once an entry besides the fixed loop's: the controller exactly the
`IntegralController` with the default band, a `retry_order` the law knows, no excursion bound, the
judge's pieces genuine (the loop's and the family's `judge`, `_lte_inputs`, the controller's
methods, `_tran_lte_c.max_error`) and unshadowed, the run's tolerance and caps of the types the C
reads.  THE BACKOFF: an entry that makes no attempt (a screen that fires every step, the SVD the
Python's -- mosrect, ab_integ, e1480, g1310) doubles the Python attempts before the next, up to
256: mosrect's 307 hand-backs became 8 and its run no slower.

Proved (test_loop_c, 26 tests): the fixed cases as before, and adaptive -- gear on a hard-driven
chain, a pulse's landings, order drops and rejections retried in C, the trapezoid, Gummel-Poon on
the chord, the PSP stage behind the screen, a 250-stage sparse run (the error test on KLU), relref
alllocal and pointlocal, retry_order assumed and carried, an LTE mask: the same bytes, times,
statistics, warnings, `_paths` tallies, every `_RUN_STATE` attribute (memos and counters aside)
and the controller's state (its running reference, retry record, carried order, last error); a
controller subclass declines (the error test's glue declines it first), an excursion bound is
refused by the judge's check; alive: the law's shrink floor a ulp off is another run.  Found on the
way, NOT this stage's: an excursion bound (`max_dv_step=5.0`) on the 8-stage MOS chain makes the
Python loop crawl (40 ns of a run not done in 60 s, at HEAD too) -- noted for its own look.

Every gold deck (27, at 0.3 of each run, the C loop on against off): the same x, t, statistics,
warnings and controller state.  The C makes 95-100 % of the attempts on 22 of them; ab_integ,
e1480, g1310 and mosrect none (the screen fires: handed back, backed off), todd3 81 %, gm17 83 %.
Wall time at 0.3 (one run each, not the line): 4049osc 0.43 -> 0.08 s, gm2 0.14 -> 0.04, gm6
0.11 -> 0.04, mike2 0.10 -> 0.04, latch 0.11 -> 0.05, nand 0.47 -> 0.24, ring 0.62 -> 0.33, gm1
0.44 -> 0.25, rich3 1.85 -> 1.12, arom 1.76 -> 1.10, mux8 8.75 -> 6.30, gm19 6.69 -> 5.36, jge
13.94 -> 11.66, voter25 2.65 -> 2.22, opampal 32.6 -> 31.8 (its operating point dominates); the
screen's four unchanged.  The line (5d) is measured full length, paired against the parent tree.

A FIX to 117, found by the readability review (below): `loop_c:steps` (then `fixed_c:steps`)
counted both the steps made in C and a decline (the run's opening steps, no history yet) -- the
tests read it as steps made.  The declines are now `loop_c:history` and `loop_c:step_bounds`.
Stale lines the rename left (the "fixed-step loop" in transient.py, `_tran_newton_c`, the pred and
ctrl docstrings, test_watch) corrected.  Pinned: "the stepping loop and its C" (renamed from "the
fixed-step loop and its C", its reference list grown by the judge's pieces: `_LMMSteps.judge`,
`_lte_inputs`, `IntegralController._evaluate`, `_tran_lte_c.max_error`, the event breakpoint),
re-recorded with "the error test and its C", "the Newton solve and its C" and "the stepping loop
and its unread Jacobian" (twin side only).  Full suite 5158 passed; sanitizer fast tier clean.

## 2026-10-11 — 121. Stage 5, step d: THE LINE MET -- the median gain over the 23 gold decks the loop covers is 1.53x (line >= 1.4x); stage 5 closed

Every gold deck full length, the parent tree (8fc455a8: the fixed-step loop in C, the adaptive
in Python) and this one (47d782b7) alternated per round and per deck, five rounds, one process
each, no address randomisation, the 15-min load below 3 at the start (scratchpad s5/line.py).
The paired ratio's median per deck; "covered" = the C makes more than half of the run's attempts.

| deck | parent | child | gain | C makes |
|---|---|---|---|---|
| 4049osc | 1.37 s | 0.15 s | 9.08x | 100 % |
| gm2 | 0.31 | 0.11 | 2.97x | 99 % |
| gm6 | 0.21 | 0.09 | 2.34x | 99 % |
| mike2 | 0.22 | 0.10 | 2.25x | 99 % |
| latch | 0.33 | 0.15 | 2.14x | 95 % |
| nand | 1.05 | 0.53 | 1.98x | 98 % |
| ring | 2.06 | 1.11 | 1.85x | 100 % |
| gm1 | 1.11 | 0.65 | 1.72x | 99 % |
| arom | 5.50 | 3.38 | 1.63x | 100 % |
| rich3 | 5.85 | 3.62 | 1.62x | 99 % |
| toronto | 1.17 | 0.73 | 1.59x | 98 % |
| slowlatch | 0.36 | 0.23 | 1.53x | 94 % |
| gm3 | 0.96 | 0.66 | 1.45x | 99 % |
| rca | 0.09 | 0.06 | 1.45x | 99 % |
| gm19 | 15.42 | 10.96 | 1.41x | 100 % |
| mux8 | 28.17 | 20.01 | 1.41x | 98 % |
| hussamp | 0.09 | 0.07 | 1.35x | 97 % |
| todd3 | 0.40 | 0.31 | 1.29x | 86 % |
| jge | 45.21 | 35.68 | 1.26x | 99 % |
| voter25 | 9.14 | 7.56 | 1.21x | 92 % |
| schmitecl | 0.07 | 0.06 | 1.19x | 97 % |
| gm17 | 0.32 | 0.29 | 1.09x | 95 % |
| opampal | 36.22 | 36.73 | 0.99x | 100 % |
| ab_integ | 5.53 | 5.59 | 0.99x | 0 % |
| e1480 | 1.07 | 1.08 | 1.00x | 0 % |
| g1310 | 1.19 | 1.20 | 0.99x | 0 % |
| mosrect | 0.26 | 0.27 | 0.94x | 0 % |

Median over the 23 covered: 1.53x; over all 27: 1.45x.  The rounds agree within ~2 % on every
deck but opampal (0.77-1.16: its KLU operating point, ~36 s, is the run) and ab_integ.  No A/A
beside it: the margin (9 %) is several times the 1-2 % a tree-and-role offset reads.  Where the
gain is small the run is elsewhere: the core's C on jge, voter25, mux8, gm19 (the sizing's
ceiling there was 1.25-1.6x), the operating point on opampal and gm17.  NOT a gain: the four
decks the screen hands back pay a cost -- mosrect 0.94x in all five rounds (~15 ms on 0.26 s),
the other three 0.99-1.00x; not yet located (the unit's load and bind, the entries before the
backoff: to be measured, not guessed).

## 2026-10-11 — 122. Stage 5e: the units' load precompiled -- a C unit loads in ~0.2 ms, not 2-3; the loop's cost on a hand-back deck 13 -> 4 ms; the smallest decks 1.3-1.7x

LOCATED (121's open item, mosrect: the C loop on cost 15.0 ms mean, 13.0 median over five
in-tree pairs, on a 0.26 s run whose every step the screen hands back): the first real entry
12.1 ms once, 998 backoff skips 2.1 ms (2.1 us each), nine entries whose Newton the screen then
discards ~1 ms.  The first entry was cffi's: a unit's load is ~90 % `ffi.cdef` (pycparser, 1.7-2.7
ms a header, ~6 ms more on the parser's first use in a process) and the first use of each type
string in `ffi.cast` / `from_buffer` / `new` ~0.2 ms (72 in the loop's `_bind`); `dlopen` is
nothing.  Every unit pays it: 16.5 ms of loads in mosrect's solve (11 units).

BUILT: `_cbuild._dlopen` takes its FFI from `_ffi_for(cdef)` -- one per declaration text a process
(the HDL kernels share the chain functions' entry), from `_precompiled`: cffi's out-of-line ABI
mode, the declaration's type tables as a generated module in the object store (`pcx_ffi_<sha>.py`,
keyed by the text and cffi's version, written to a temporary file and renamed into place),
parsed as before when it cannot be written or read; every library opened is held for the process
(an out-of-line FFI's library closes its object when collected).  `_tran_loop_c.run` checks the
backoff's skip first (the `_watch` import hoisted): a skip is now ~0.6 us.

Proved: for every unit (15) the precompiled FFI declares what parsing declares -- every typedef,
every struct's size, alignment and fields (offset, type, bit width: newton 87, loop 77, radau_tc
72, lte 54, ...), and with the unit's object opened through each, every function's type and
every enumerator's value; the kernels' entry signature; alive: one field's type changed and one
enumerator moved are each caught (the plant first caught the check itself comparing typedef
NAMES, not layouts -- fixed before it counted).  Twelve processes racing on one new module all
got it precompiled, one file, no temporary left; a truncated module falls back to parsing.  Full
suite 5187 passed; sanitizer fast tier clean; every gold deck's final x and t one digest, the
same as 8b5d2cca's.

THE LINE (written before building): loads during mosrect's solve 16.5 -> 2.4 ms (<= 3: met);
mosrect's loop on-minus-off 13.0 -> 4.05 ms median over 15 pairs, 4.94 mean (<= 5: met; the first
build alone read 7.1 -- the skip-first check closed it); no gold deck slower in all five paired
rounds (met).  Paired against 8b5d2cca, full length, five rounds: rca 1.66x, schmitecl 1.64x,
hussamp 1.48x, gm6 1.33x, mike2 1.32x, gm2 1.28x, latch 1.18x, 4049osc 1.17x, mosrect 1.14x,
slowlatch 1.12x, todd3 1.10x, gm17 1.09x; the long decks 1.00-1.05x; median 1.04x over all 27.
gm19 read 0.96x (rounds 0.997, 0.963, 0.912, 0.963, 1.003): its stepping is 0.03 s FASTER in
every one of eight further pairs (6.77-6.81 s against 6.82-6.84 s); its operating point takes
3.6, 4.1-4.2 or 4.4-4.6 s by process in BOTH trees, the same answer -- a matter of its own.

A MEASUREMENT LESSON: the first probe walked the stack (`traceback.extract_stack`, which reads
source lines) INSIDE the timed region: each load read ~0.3 ms high, the baseline 21 ms instead of
16.5 and the build's 5.4 ms instead of 2.4.  Take the time, then describe what was timed.

## 2026-10-11 — 123. Readability stage R3, steps a-d: the map, the comments, the tables, the names, the typed wiring -- no answer moved, no case slower

Andreas: "Start R3" (a readability stage after stage 5; "it does not need to be small"). Two
read-only reviews sized it (the plan's R3 section); seven steps, a-d done:

* R3a (a8709599, a9705bf6, 5a83f1f8): `doc/transient_runtime.md`, the map of the C runtime, with a
  test that it names every C file, glue module and switch; every C struct field documented by
  role, every entry's contract, the why of every bit-identity form (comments only: all 26 files
  token-identical, all 15 units' machine code identical); the glue's docstrings say what each
  module stands for and owns; the history (81 lines of Python, the C's plan tags) moved to
  `doc/transient_history.md`; tests keep history out of the C and the glue.
* Fixes R3a's review found (4fa1530f, 68079113): the C loop ignored `PYCIRCUIT_PRED_C=0` and
  raised where the predictor could not load; a Pulse's zero edge set after `__init__` was served
  where the Python raises; the core packed an unpacked element with the FIRST batch's kernel
  (another class's layout); `_RUN_STATE` lacked three records. Each with a test that fails before.
* R3b (beb97d51, 33f00b94, f56cb48b): the C loop's 57 pieces one table (`_PIECES` -> `_own`,
  `_pieces`, `_SHADOWS`; the copies had drifted -- `_pred_fast_nodes` unshadowable, now declined);
  every decline reason the run's or the attempt's in one table (`_KEEPS`; three reasons that meant
  two things split); the write-back's counts in `_count`, each with its Python site, and a moved
  count named by the A/B. `--check` against 68079113: 18 of 18 no change.
* R3c (2a9e0998, 02510917): the glue says the codes the C reads by the headers' names (a test
  forbids the numbers); one flags macro and THE FLAG RULE in `pcx_fenv.h`; the predictor's, the
  core's and the walk's capacities and the core's ztab columns as header enums both sides read;
  the shrink guard named on both sides, a test holding them (machine code identical).
* R3d (d2b682b6, e009ac19): every unit's entry `pcx_<unit>` (loop.c's six `#define hdl_fn`
  renames gone; its calls typed by `__typeof__` of the callee's own declaration -- .text
  identical, the exports renamed alone); newton_t's dead fields out, its locals local.

Proved at R3d's end against f56cb48b (before R3c): `--check` 18 of 18 no change (spice_ce_pss
+0.9 % [+0.7, +1.3] under the +1 % line, pss -0.7 %, ladder_radau -0.5 %); all 27 gold decks one
digest, the same in both trees. Found on the way and recorded in the plan: the gm19 operating
point's per-process modes; hand-backs on underflow (numpy ignores it); a few latent cases in the
Radau glue to test. Steps e-g (one job per function in C and Python, shared helpers) remain.

## 2026-10-11 — 124. Readability stage R3, steps e-g: one job per function, the shared helpers -- R3 closed; no answer moved, no case slower; one criterion not met

Andreas: "Start" (continue R3 after the a-d checkpoint).

* R3e (9113cec3): newton.c's `iterate` and loop.c's `pcx_loop` as their named phases (the Newton:
  finish, cache hit, solve, limit, converged, fold, chord breaks; the loop: where, take, judge,
  accept, one attempt's values in `loop_attempt_t`), so the C loop reads as
  `_SteppingLoop.execute`. `--check` against e009ac19: 18 of 18 no change; gold decks the same
  digests.
* R3f (3c74b517): `_tran_loop_c.run` as its jobs (why not the run, the sources, the read, the
  bounds, the load, the call, the counts); the per-step `_tran_newton_c.solve` and
  `_tran_lte_c.max_error` sectioned and numbered rather than split (each a call a step); the C
  loop's state on the Python loop one object (`LoopState`, `attach`).
* R3g (a36646b6): `csrc/pcx_minmax.h` -- Python's `max`/`min` and numpy's `maximum`/`clip` by name,
  the 21 hand-written ternaries of five units saying which rule each follows (all 15 units'
  machine code identical; pred.c keeps its own written out: the helpers moved its instruction
  schedule, and a unit that tests flags keeps the code it was proven with);
  `_cbuild.openblas` / `address`, one way to find numpy's and SciPy's own routines (five copies
  of the lookup, three spellings of the cast); one `_undo_source`.

Proved at R3's end against 9113cec3 (R3e): `--check` 18 of 18 no change (SLOWER none, FASTER
none); all 27 gold decks one digest, the same in both trees; full suite 5221 and the sanitizer's
fast tier clean at each commit. opampal's one-round 1.38x was a slow parent round: over three
rounds 0.99x (0.99, 0.98, 1.01).

Decided, not shared (the plan says why for each): the eight `driver()`s (each states its own
reasons), the buffer binders and walk records (each struct's), `cmul` and the bail helpers
(different by design: FMA order, flag clears), radau / radau_tc's common lines.

The readability line, item by item: the map names every C file, glue module and switch; every
struct field documented; no history in the C or the glue; no literal code across the boundary --
MET, each held by a test. "No function over ~80 lines without a measured cost in its docstring"
-- NOT MET: `_tran_newton_c.solve` 276 lines, `_tran_lte_c.max_error` 204, `_Core.__init__` 192,
`evaluate` 116, `probe` 83, `_tran_loop_c._bind` 104, `run` 85; C `pcx_core` 148, `pcx_lte` 123,
`pcx_pred` 110, `judge` 97, `iterate` 85. Left for a follow-up: split the run-once ones (`_bind`,
`_Core.__init__`, `judge`); measure a split's cost on the per-step ones and write it where they
stay whole.

## 2026-10-12 — 125. Readability stage R3h: the long functions split -- every function of the loop units and their glue at most 80 lines, a test holds it; no answer moved; the per-step split's cost measured, inside the line

Andreas: "Split now" (R3's one unmet criterion: fourteen functions over 80 lines). The line,
written first: every answer byte-identical; `--check` no case SLOWER and the four gold decks
the C loop does not cover >= 0.99x; every function in scope <= 80 lines or `KEPT WHOLE:` with
its measured cost, a test holding it.

* (a) 3d7473be, the run-once Python: `_tran_loop_c._bind` (its buffers, the units, the screen,
  the judge; `_plug` one way to point the loop at a unit) and `run`'s checks (`_eligible`; the
  gate stays: a cheap return made once a Python attempt), `_Core.__init__` (`_Flat`,
  `_lone_batches`, `_struct`, `_index_elements`) and `probe`'s re-check, `_hdl_climit.render`.
  Shown unchanged: every core build's digest on 29 circuits (A/A-stable; two planted faults
  moved it), every limiter's C text on 32 classes.
* (b) 52b02a42, the C: loop.c `judge` (`error_test`, `retry_at_point`), newton.c `iterate`
  (`evaluate`), lte.c `pcx_lte` (five helpers, each flag-tested region whole in one), core.c
  `pcx_core` (`fused_groups`, `one_pass`), pred.c `pcx_pred` (`fast_nodes`, `fit`,
  `clamped_product`) -- moved along anchors by script; the five units' machine code moved
  (contraction off: no result can), the two flag-testing units run against their old selves
  on 400 000 / ~360 000 drawn calls, identical, each harness alive to a planted flag test.
  `test_c_abi`'s bare-code check widened to a returned code.
* (c) f8a95566, the per-step Python: `_tran_newton_c.solve` (`_admit`, `_companion`, `_record`,
  `_source`, `_load`, `_serve`), `_tran_lte_c.max_error` (`_admit`, `_scalars`, `_arrays`,
  `_load`), `_tran_core.evaluate` (`_companion`, `_c_inputs`, `_leave`), `_hdl_climit.limit_walk`
  (`_take_elements`, `_run_walk`). Timed standalone against 52b02a42: `solve` +1.8 us a call,
  `max_error` +0.37, `evaluate` +0.06, `limit_walk` +0.13 where it serves -- about 2.3 us of a
  ~240 us mosrect step; the four uncovered decks over seven rounds mosrect 0.9954x, ab_integ
  0.9977x, e1480 0.9991x, g1310 0.9990x. `test_runtime_map` holds the 80 lines;
  `test_paths` accepts a decline passed on (`# paths: counted by <helper>`) only from a
  helper it checks itself.

Proved at R3h's end against 445d5a2f: `--check` 18 of 18 no change (SLOWER none, FASTER none);
all 27 gold decks one digest, the same in both trees; full suite 5233 (+12 tests). One round's
outliers re-measured: nand 0.9996x, hussamp 0.9956x, voter25 0.9986x over five rounds.
opampal is unresolvable: its run takes 31.5..43 s from process to process in either tree
(a fixed hash seed does not remove it), nine paired rounds median 1.02x, the answer the same
-- a per-process mode like gm19's operating point, an open item.

## 2026-10-12 — 126. The DC operating point, D0-D3: the ladders' dense temporaries out, the Jacobian assembled reduced -- opampal's operating point 1.52x pinned, 1.22x unpinned, the per-process modes gone; every answer and warning the same

The C loop plan done, its "after this plan" list sized (D0): the operating point was a third of
the 27 gold decks' time, nearly all opampal's -- ~19 900 Newton iterations in the gmin /
gshunt ladders, each dense 519 x 519 numpy work (the gshunt `J + eye(n) * g`, junction gmin's
copy and Python loop, the dense G, its reduced copy, KLU's whole-matrix nonzero count), the
KLU solve under 5 %. Its 31..43 s spread by process was glibc's dynamic mmap threshold.
Plan `dc-operating-point`; D1 the instrument (an operating-point census, x0 and warnings
digests, a timer pinned and unpinned -- unpinned A/A reads up to 7.8 %).

D2 (`_plus_diagonal`, `_junction_gmin`: bytes and warnings the expressions') read opampal 1.35x
pinned but 0.93x unpinned: the in-place diagonal, allocating less, put glibc into its
trim-and-regrow mode (12-17 M minor faults against 4-6 M). Andreas: fold it into D3. D3
(`assemble_matrix(drop=)`: each kept bin summed from its slots in order -- one (n-1)^2 array
for an n^2 and an (n-1)^2) took the pair of large arrays a call away, and the mode with it.
f7d0be3a, against 658e6756: opampal 1.524x pinned, 1.22x unpinned (26.6 s every round); gm19
1.31x / 1.67x; jge 1.13x; gold decks full length: opampal 1.50x, gm19 1.08x, none below 0.99x,
every digest the same. Next: D4 (the evaluations on the C core), sized from what remains.

## 2026-10-12 — 127. The DC operating point, D4: the base Newton solve in C -- opampal's operating point 4.0x pinned, 5.5x unpinned over D3, every base solve of every gold deck the Python's, byte for byte

Sized first (no code changed): opampal's operating point makes 448 base Newton solves (19 737
iterations), 186 failing at the iteration limit -- so the C must run failing solves to their
end, not hand them back; StructuralKLU never foreign on the DC, the J on its pattern by
construction. Built (0389f9b0): `csrc/dcnewton.c`, one `StandardNewton.solve_system` on the
compositions the ladders make -- the DC residual, source stepping's factor, gshunt, junction
gmin -- read from tags the wrappers carry; numpy's dgesv or StructuralKLU's own refactor /
fresh factor / residual sequence; the limiter walk; the convergence test with numpy's dgemv.
A flag, a non-finite value, a failed factor or a walk stop hands the solve back with KLU's
numeric as it was, the Python's rerun replaying it exactly; the iteration limit is served and
raised by the Python's own `_standard_failure`.

The census: 717 base solves of the 27 gold decks, each solution's bytes and iterations (or
exception text), KLU's counters, x0 and the warnings -- C on against off, 27 of 27. Line, against
67bda752: opampal's operating point 4.03x pinned / 5.54x unpinned (>= 3x), gm19 2.45x / 2.43x
(>= 2x), gm3 2.44x, mike2 1.28x; full gold runs pinned: opampal 3.54x, gm3 1.57x, gm19 1.19x,
none below 0.99x.

The first gold line failed on seven tiny decks (0.97..0.99x), none of them served. Two causes,
both outside the solve: a `PCX_SKIP_Q` bit added to the core's pass loop (the DC never reads
q) cost the transient's core -- 4049osc, no DC solve at all, 0.988x with it and 0.999x without,
so the bit was dropped and core.c is untouched; and the glue's import, 0.2 ms (plus 0.35 ms of
C source reads at import), landing inside a small deck's first operating point because
nrsolver itself is imported there: rca's operating point read 0.90x, hussamp's 0.88x, though
rca's one solve was declined. The fix is a gate in nrsolver ahead of a lazy import: below
`DC_C_MIN_N` = 32 unknowns (gm6, 20: 0.97x served) and the run's first solve (AutoSolver's
choice not made) are declined there, counted; the C sources are read at the first load. rca,
schmitecl, hussamp's operating points then 1.00x, 1.00x, 0.99x. Not served yet: held solves
(`.ic` without UIC: gm17, jge's 1.2 s).

## 2026-10-12 — 128. PCNR's own ladders: 24 of the 26 gold operating points solved by PCNR, from 17

The road here. KLU below 250 (entry 127's sequel) "failed" gm3's transient and e1480's operating
point; traced, both were non-unique or non-contracting Newton solves -- seven gold operating
points are SINGULAR today (a MOS2/3 series-stack node with no DC path), and the ordinary DC
Newton's failing ladders come from the state-free limiter's write-back (a later device moves a
shared or source-held node and undoes an earlier device's junction limit; gm3 fails 21 of its 27
solves that way). A standing junction gmin (roadmap Ag) makes the seven unique and matches
ngspice-47 (e1480 5.0000 V against 4.99992); a SPICE-style limiter prototype (each device
evaluated at its own limited point) converged all 26 with a ladder. Andreas chose the PCNR route
-- PCNR gives each limited quantity its own unknown, so no write-back exists -- and asked first
for PCNR's own ladders, in Python. Full record: the plan `dc-final-gmin`.

fc4999a5: `pcnr.solve_dc` gains `gshunt`, `source`, `v_lim0` (inert at their defaults);
`pcnr.solve_dc_continued` runs the plain solve, then a gshunt ladder (1e-2 down by decades,
geometric retries), then source stepping; `DC(pcnr=True)` reports `pcnr_path` and falls back to
the ordinary chain only where all three fail. Ladders had been tried once (2026-08-26) and failed
because every rung began from an unlimited seed -- since fixed in `v_lim_init`.

Line met. With the standing gmin (scratch elements): PCNR solves 24/26 without falling back
(plain 17, gshunt 6, source 1 -- opampal); gm1 and jge fall back. Without it 23/26. Bytes
unchanged wherever PCNR solved plainly; 18 decks equal the ordinary chain to 1e-6, the other
four non-unique points (both answers satisfy KCL to 1.1e-11 A). PCNR is off by default: no
default answer moved; suite 5253. Speed is Python's: e1480 200 s, opampal 161 s, gm19 104 s,
jge 812 s before falling back. Measured for the next step: batching `pcnr_didv` per model class
is bit-identical and 3x; evaluating devices through the compiled C kernels at their own points
is the real lever, moving PCNR's answers at rounding level.

## 2026-10-13 — 129. `lte_algebraic='charge'`, the default: a current that closes a capacitive loop judged by the charge it moves

**Why now.**  The standing gmin (the DC plan's final stage) broke jge's
transient at 2.186 ns.  Traced, the error test's worst component was the
clock source VPh1H's current -- section 76's mechanism, not the gmin's:
with the gmin the step had to fall below ~4e-18 s, where the step matrix
breaks; without it the run survived by 143 force-accepted steps.  And
'judge' already gave up on the gold decks: 384 steps force-accepted over
tolerance, on 15 of the 27.

**What those currents are.**  A voltage source whose two ends are joined
by capacitors and other voltage-defined branches -- a supply or a clock
driving gate capacitance -- is where MNA's index is 2: its current is the
slope of capacitor voltages, determined by them, not integrated, and its
truncation estimate does not fall with the step (76: ~5e-6 A from h =
1e-12 to 1e-16 s).  RADAU5 scales such components' error by h; ngspice
judges device charges in charge units with a charge floor, never a
source's current.

**`lte_algebraic='charge'`, the default**: those currents judged by the
charge they move in a step -- the absolute tolerance max(lte_iabstol,
lte_chgtol / h), `lte_chgtol` 1e-15 C (`stepcontroller.charge_floor`),
never below 'judge''s -- every other unknown as under 'judge'.  Which
currents: structural, from C at the error test's generic state, once a run
(`_RunHistory._lte_charge`; `_tran_history.capacitive_loop_currents`: the
charge-free branch currents whose ends the capacitive couplings and the
other such currents join, a bridge test on the groups the couplings make,
checked against the definition on 300 drawn circuits).  The multistep
error test and the stage family's take the floor; the coupled path already
leaves every charge-free row out (P22).  The C error test floors the same
rows the same way (lte.c: `pcx_np_maximum(abstol, chgtol / h)`), the C
loop binds them once a run; both pinned pairs re-recorded after the twins
(the drawn bit-for-bit test with drawn rows, a pulse driving a capacitor
in the C loop: the same run, and another under 'judge').  The suite's four
'skip' OVERRIDES are gone.

**Sized before it was built** (a prototype, the Python loop; the DC plan's
"ALGEBRAIC ERROR CONTROL"):

| the 27 gold decks | steps | rejected | forced | time | gold columns passing |
|---|---|---|---|---|---|
| 'judge' (+ the four 'skip' overrides) | 73 342 | 17 685 | 384 | 133 s | 69 / 99 |
| 'skip' everywhere | 51 020 | 7 155 | 1 | 81 s | 69 / 99 |
| **'charge'** | 53 228 | 8 224 | **1** | 87 s | 69 / 99 |
| a looser band, gamma_max 3 / 8 | 77 920 / 71 077 | 7 324 / 3 054 | 173 / 38 | 123 / 107 s | 69 / 99 |

The one forced step is 4049osc's, in every mode.  'charge' moves 44
columns closer to Xyce's gold, 29 further, 26 not at all; none crosses the
pass line.  jge: 11 036 steps / 143 forced / 44.5 s -> 4 490 / 0 / 17.7 s.
A looser band (VACASK accepts ~55x the raw tolerance before a rejection,
we 7x) is no remedy: jge still forces 26 steps at gamma_max 3, and without
the overrides voter25 goes 4 006 -> 11 962 steps.  chgtol: 1e-16 takes
8 % more steps, 1e-14 (ngspice's) 3 % fewer but forces a step on two
decks.  On radau: 227 147 forced (todd3 224 646) -> 20, steps 486 220 ->
39 280, the same 69 columns passing.

**The suite.**  'skip' as the default fails 10 tests, 76's six accuracy
ones among them; 'charge' none -- those circuits' currents close no
capacitive loop, so 'charge' judges them as 'judge' does.  Of 655 tests
that reach the error test, 46 see such currents and 11 have a step
decided otherwise; all pass.

**Not decision 0.3d's charge guard** (`benchmarks/transient_decisions.py`):
that was a second acceptance test on every unknown, and it switched step
control off; this floors only the index-2 rows, and the decks keep 8 224
rejections.

**Gate.**  The full suite: 5306 passed (11 of them new), none failed.  Every
gold deck on the C paths -- the C loop served 47 330 of the 53 228 steps --
gives the prototype's numbers exactly, deck by deck: steps, rejections,
forced steps and every gold metric.

## 2026-10-13 — 130. 'charge''s floor 1e-18 C, and six rejections where it floors a current: those currents as accurate as under 'judge', no step forced

**What 129's gate could not see.**  Its score was the gold columns, node
voltages but one floored current of 99.  Measured on the currents
themselves -- the 21 gold decks with a loop-closing source, each against a
fixed-step Gear-2 reference of 1e5 steps; the waveform's RMS error over the
reference's, from 1 % of the run on (a start from zeros puts a spike in the
first step whose height is the step's), sources carrying only rounding left
out -- the worst current per deck was 15.1 % at 1e-15 C (the median over
the decks), 3.2x 'judge''s (up to 15.8x) and near 'skip''s 20.2 %; the
charge it moved within 2.5 %; node voltages within tolerance in every mode.
The floor's size was the cause: TRTOL chgtol / h is 0.7 mA at a 10 ps
step, so those currents went almost unjudged at ordinary steps.  (A
reference of 4e5 steps moves the larger errors by less than 0.3 points.)

**The fix.**  A smaller floor restores them -- 1e-17 C: 6.1 % (1.06x
'judge''s), 1e-18 C: 4.2 % (1.00x, the worst deck 1.6x) -- but forced
steps return, and traced they are the multistep family's cap: at a corner
the floored current's error falls 714, 227, 44, 4.8 times its tolerance
over four tries, and the third rejection forced the step.  So six
rejections where 'charge' floors a current, and none forced.  Not where it
does not: without the floor such an estimate may never fall -- six let
gm3's step under 'judge' shrink to minstep and fail where three
force-accept it -- and the comparator oscillator, with no floored current,
took six as a change in its event timing.  A floor on retries only (the
first try judged in full) bought nothing (gm3: 6.3 % against 7.1 %).

| the 27 gold decks | steps | rejected | forced | gold columns passing |
|---|---|---|---|---|
| 'judge', the default until 129 | 73 342 | 17 685 | 384 | 69 / 99 |
| 129: 1e-15 C, three rejections | 53 228 | 8 224 | 1 | 69 / 99 |
| **1e-18 C, six where it floors** | 78 883 | 21 867 | **0** | 69 / 99 |

jge with the standing gmin completes at 1e-18 C: 10 667 steps, none forced.

**Radau** keeps its own budget of twelve: at 1e-18 C it forces near
'judge''s counts at the corners (jge 1 007, mux8 633, arom 92); a budget of
40 removes them but lets arom's step reach minstep and fail.  Its step
controller -- the RMS norm over every unknown, a fixed-order retry -- is
the limit; sized next (the DC plan's R).

**Gate.**  The full suite: 5306 passed, none failed (the first run caught
the two force-accept tests and the comparator's jitter while six applied
to every 'charge' run).  Every gold deck on the C paths gives the
candidate's numbers exactly, deck by deck -- steps, rejections, forced
steps and every gold metric; the C loop served 72 825 of the 78 883
steps.  jge with the standing gmin: 10 667 steps, none forced.

## 2026-10-13 — 131. The stage family's retry: the worst unknown at the order its tries measured -- radau's forced steps 2 190 -> 2, gm19 completes

**What 130 left.**  Radau, TR-BDF2 and ESDIRK43 judge a step by their
embedded estimate's RMS over the unknowns and, until now, sized a retry by
the accepted step's law: that RMS at the fixed order p, the shrink bounded
by one half.  At a loop-closing source current's corners one unknown is
over its tolerance and the rest are not, so the RMS diluted a 7x excess to
~1.1, each retry shrank ~12 %, and the twelfth rejection force-accepted
the step: radau at 1e-18 C forced 2 190 steps on the 27 gold decks
(jge 1 007, mux8 633, rich3 132, arom 92) and failed gm19 on a singular
matrix.  A budget of 40 removed the forced steps but let arom's step reach
minstep.

**The fix.**  The multistep family's retry (2026-10-10,
`stepcontroller.RETRY_ORDER_MIN`) on the WORST unknown: from a point's
second try, the order its two tries measured, clamped to [1, p]; the
shrink bounded by `MIN_SHRINK_RATIO` (0.2).  Acceptance is unchanged --
the RMS <= 1 -- so nothing is loosened: only the step after a rejection
moves.  `retry_order` now means the same in both families ('assumed' keeps
p, 'carried' carries the last measured order).  Three rules were sized on
radau, every gold deck (the DC plan's R):

| radau, 27 gold decks | steps | rejected | forced | failed | solve time |
|---|---|---|---|---|---|
| the old retry | 63 537 | 45 733 | 2 190 | gm19 | 599 s |
| RMS, measured order, bound 0.5 | 68 771 | 34 036 | 438 | -- | 627 s |
| RMS, measured order, bound 0.2 | 70 340 | 23 970 | 56 | -- | 554 s |
| **worst unknown, measured order, bound 0.2** | 72 949 | 17 879 | **2** | -- | 515 s |

The 2 are 4049osc's start (the all-zero state jumping to the sources' in
the first step), forced under every rule.  Without gm19: 10 % more steps,
27 % less time.  Every gold column passes or fails as before (69 / 99);
each column's metric against the old retry's has median 1.00, range
0.70-1.45.  Against the budget-40 runs (9 decks): 2-13 % more steps,
faster on all nine (jge 212 -> 118 s).

The same rule on the family's other members (the old retry against the
new, 27 decks each):

| | steps | rejected | forced | solve time | gold |
|---|---|---|---|---|---|
| esdirk43 | 30 887 -> 32 502 | 3 868 -> 2 273 | 35 -> 1 | 353 -> 305 s | same |
| trbdf2 | 84 495 -> 96 383 | 57 714 -> 25 306 | 1 583 -> 2 | 938 -> 719 s | same |

No run fails under either, the 2 and the 1 are 4049osc's start, and each
column's metric has median 1.00 (TR-BDF2 0.69-1.41, ESDIRK43 0.63-1.65);
the counts were taken in a worktree at 130's commit with the new rule
patched in, the times beside a running suite (both rules interleaved deck
by deck).  The GLMs step through the same judge; the suite covers them, no
gold deck was run on them.

One case differs from the multistep rule, as measured: a worst error that
did not fall since the previous try is retried at p (the multistep rule
takes 1 there).  Under the old 'judge' mode (not the default) radau on
gm3 now takes 1 213 steps (39 forced) against 593 (40): the supply's
unfloored current cannot be satisfied, and the stronger shrink spends
steps on it.

**Gate.**  The full suite: 5 309 passed, none failed.  Every gold deck on
radau with the built code gives the prototype's numbers exactly, deck by
deck -- steps, rejections, forced steps and every gold metric.

## 2026-10-14 — 132. The early abort: an attempt whose error estimate at its first Newton iterate is already over twice the tolerance is rejected there -- 9 % fewer Newton iterations, 6 % less solve time

**Where the rejections are.**  Andreas asked how to reduce the rejections
(21 867 of 100 750 attempts on the 27 gold decks after 130).  Classified
on every deck (the Python loop, which made the C census exactly): all but
two are the error test's; 84 % have a loop-closing source current as the
worst unknown (node voltages alone would pass 79 %); rejected attempts
take 28.5 % of the Newton iterations.  They are the price of judging those
currents at 1e-18 C (130): where node voltages set the step our rejection
rate is ngspice-47's (arom 24 vs 23 %, mux8 24 vs 20), where the currents
do it is not (jge 26 vs 0.4 %) -- ngspice's error test reads charges only.
A first retry sized for the currents' measured order (three variants) cut
the rejections 20 % and the Newton work 1-3 %, the larger shrink paid back
in steps: refused.  `lte_chgtol` 1e-17 would save 17 % at 1.06x the
currents' error (the worst deck 4.0x): left as it is.

**The fix: stop paying for a rejection's Newton.**  The error estimate at
the attempt's first iterate -- the error test's own chain on x1 and its
charge, J^-1 by the factorisation of the seed's Jacobian the Newton already
holds (dgesv's LU, the chord's factors, KLU's numeric) -- told the
rejections apart: over 2, 18 063 of the 21 863 and 9 accepted steps.
Abandoned there, the attempt is rejected at that estimate (the controller's
law, its measured-order record; `IntegralController.reject_early`).  Sized
first in Python (the retry sized from the estimate: -9.4 % Newton
iterations, none forced; from a quarter of it, or halving the step, the
point's budget ran out and steps were forced -- refused), then in C in a
worktree against a refusal line of 3 % (-6.5 % unpinned, -7.1 % pinned,
A/A within 0.6 / 2.9 %).

**Built** (`Transient(early_abort=2.0)`, 0 off): armed by the stepping
loop for an adaptive multistep attempt with history while a rejection
would leave a full attempt in the point's budget (`_SteppingLoop.
early_arm`); the check in `_StepNewton._newton` at the solve's second
evaluation (`_early_estimate`, under `np.errstate(all='ignore')`), its C
twins in newton.c (`early`, on lte.c's factored solve) for the Python
loop's C Newton (`_tran_lte_c.early_context`) and the C loop (loop.c's
`take` and `early_reject`).  Not on the full Newton of a sparse matrix (no
factors held), a 2nd-order step with no third point, or a non-finite
estimate.  `TransientStatistics.early_aborts`.

| the 27 gold decks | steps | rejected | Newton iterations | forced | solve time |
|---|---|---|---|---|---|
| before (130 + 131) | 78 883 | 21 867 | 240 476 | 0 | 96.6 s |
| **early abort at 2** | 78 989 | 21 901 (9 302 abandoned) | **217 833** | 0 | 90.6 s |

Every gold column passes or fails as before.  Solve time, the built code, `early_abort` 0 and 2 in one
process per deck, [off, on, off] x 3: pinned -6.2 % (sum of minima; A/A
0.0 %), -5.7 % (medians; A/A +0.1 %); unpinned -6.9 % (A/A -0.6 %), -6.4 %
(A/A -0.1 %); per deck median 0.955 -- todd3 0.49, schmitecl 0.60, latch
0.70, voter25 0.90.  ⚠ THE PYTHON LOOP PAYS FOR THE ARM: the four decks the
C loop does not serve (mosrect, ab_integ, e1480, g1310) prepare the error
test's context in Python at every armed attempt (`_tran_lte_c.
early_context`, ~22 us), so a deck with no rejection loses: mosrect 1.12x
(388 checks, no abort, 0.23 s), ab_integ 1.03x; on the C loop a check costs
~1 us (4049osc 1.02x on 3 676 checks).  The context kept across a run, its
step fields alone refreshed, is the follow-up.  ⚠ EQUAL CONVERGED POINTS
NO LONGER MEAN EQUAL STEPS: whether an attempt is abandoned is decided at
its first Newton iterate, so two Newton routes to the same points take
different sequences (PCNR 186 / limiting 197 steps on gate 13-6's
rectifier); three tests whose subject was measured on the old stepping run
with it off (gate 13-6, the comparator's jitter, PSS's closing column).

**Gate.**  The three paths -- the Python loop and Newton, the Python loop
on the C Newton, the C loop -- give the same steps, rejections, Newton
iterations, aborts, checks and every gold metric on all 27 decks; with
`early_abort=0` every deck is the census before.  The full suite: 5 318 passed (9 new), none failed; the
pairs re-recorded after their twins (the controller law, the stepping
loop, the Newton solve, the error test, the loop's unread Jacobian).

## 2026-10-14 — 133. The branch check's certificate: the rank proved before the SVD -- four gold decks run in the C loop, 6.8 % less solve time on the 27

**What the screen cost.**  Sizing the early abort's setup on the four gold
decks that ran outside the C loop (mosrect, ab_integ, e1480, g1310) found
why they were there: the branch screen's cheap level handed every
converged attempt to the SVD -- a constant small capacitance at 1e-8 to
4e-7 of the running maximum, under its `1e-6 * ref` test -- and the SVD
never fired (0 of 5 834; the rank at its structural value, the smallest
singular value above tol 1.6x to 191x it).  The C loop bails on such a
screen.  The SVD was most of their time (ab_integ 3.29 of 5.54 s, e1480
0.67 of 0.99), and with the screen costing nothing they ran 5-8x faster in
the C loop (the four 7.76 -> 1.40 s).  Sized and refused on the way: the
early abort's context kept per run (a third of its setup back), that
refresh in C (slower than the Python's: the glue is the cost).  Sized: the
SVD in C (full 5.60 s, values only 3.32 s -- the LAPACK time stays).

**The certificate.**  The screen asks whether sigma_r0(C) <= tol.  For any
r0 x r0 submatrix B of C, sigma_r0(C) >= sigma_min(B) (interlacing); a
block-diagonal B has its blocks' singular values; sigma_min(B_k) >=
1 / ||B_k^-1||_F.  So every block's 1 / ||B_k^-1||_F above tol (1 + delta_k),
delta_k = 10 k eps ||B_k||_F ||B_k^-1||_F (the inverse's rounding at the
block's own conditioning), proves the rank without the SVD.  B: C's
structurally nonzero rows and columns without the reference (where they
are more than r0, an r0 subset by pivoted QR), its blocks the generic
pattern's pieces -- made once, from the structural rank's generic C, at the
first screen the cheap level fires on.  It never fires: what it cannot
prove (a bound under tol, an entry between two blocks, a singular block)
the SVD decides, as before.  Prototyped in C in a worktree first:
identical runs on all 27 decks, 5 834 of 5 834 proved, the four 7.83 ->
1.61 s.

**Built** in the reference (`_BranchCheck._branch_screen_C`: the cheap
level, the certificate, the SVD; `_branch_certify`, `_branch_cert_data`,
`_certify`, `_cert_data`) and its C twin (branch.c's `certify`, on numpy's
own dgesv -- numpy's inv is dgesv on the identity -- with the sums of
squares sequential as numpy's cumsum makes them, bit for bit), bound by
`_tran_branch_c` on the Python loop's screen and by `_tran_loop_c` on the C
loop's own once the Python has the run's data.  `Transient.BRANCH_CERT`
(True; False: every fired cheap level to the SVD), `BRANCH_CERT_GROWTH`
(10), `BRANCH_CERT_MAX_N` (1000: a larger generic C is not kept for it).
`_paths` `branch_c:cert`.

**Gate.**  Suite 5 328 passed.  The census on all 27 gold decks, three
paths (every C off, the Python loop, the C loop): every count and gold
metric as at 132 (78 989 steps, 21 901 rejections, 217 833 Newton
iterations, 0 forced, 9 302 aborts); every screen the cheap level fired on
was proved -- 5 837 on the pure Python path, all on the four decks, each
one the SVD also found quiet; the C loop's steps 73 108 -> 78 829.  Timing
(one process a deck, `BRANCH_CERT` off / on / off interleaved, three
rounds, pinned and unpinned malloc): the 27's min-sum 90.78 -> 84.62 s
pinned (-6.8 %, A/A +0.1 %), 90.89 -> 84.73 s unpinned (-6.8 %, A/A
-0.1 %); mosrect 0.13x, e1480 0.22x, g1310 0.32x, ab_integ 0.19x; the
other 23 +0.03 % / +0.07 %, none outside 3 %.  e1480's tightest proof
stands 1.4 % above its bound: a step the certificate cannot prove there is
one SVD, as before.

## 2026-10-14 — 134. The early abort's threshold 1.25 (was 2): 5.8 % fewer Newton iterations, 3.3-3.7 % less solve time; seeding a retry from the rejected attempt refused

**Where the rejections' work was.**  Andreas asked what more could cut the
rejections or their cost.  The census on the current code (after 132 and
133; the pure Python path, the same run): 100 889 attempts, 78 989
accepted, 9 302 abandoned at their first iterate, 12 598 judged and
rejected -- those took 16.7 % of the Newton iterations (2.9 each), the
abandoned 4.3 %.  Every judged rejection had estimated 2 or less at its
first iterate (5 507 in (1, 1.25]); of the accepted attempts only 74 were
in (1, 1.25] and 7 above.  Retries cost 2.46 Newton iterations against
2.07 for first tries, 29 % of them rejected again; first rejections came
mostly at a constant step with the error 2-5x what the last step
predicted -- an edge arriving, which a controller looking back cannot
foresee; 79 % were judged on a floored loop-closing current (the accuracy
target, 129-130).

**Sized against lines written first** (the Python loop and Newton, 27
decks): (1) a lower threshold alone needs >= 2 % of the Newton
iterations with no forced step, no gold column across its pass and no
deck's Newton up over 5 %; (2) seeding a retry needs >= 3 points more than
that threshold alone, likewise; (3) then >= 2 % solve on the C paths.
SEEDING -- the retry predicted with the rejected attempt's converged point
(an abandoned one's first iterate) as an extra node, `_predict_state`'s
`extra`, an interpolation inside the rejected step -- brought a retry to a
first try's cost: alone at 2, -4.07 %.  Thresholds alone: 1.5 -2.88 %, 1.25
-5.79 %, 1.1 -7.37 % (1.5 and 1.1 each slowed a deck -- todd3, slowlatch --
by over 5 %); seeding on top of 1.25 added 2.56 points and slowed
slowlatch 14.9 %: REFUSED by (2).  1.25 passed (1); on the C paths (one
process a deck, 2 / 1.25 / 2 interleaved, three rounds): min-sum
84.76 -> 81.66 s pinned (-3.7 %, A/A 0.0 %), 84.98 -> 82.15 s unpinned
(-3.3 %, A/A -0.6 %), the C paths' counts the Python sizing's on every
deck, no deck more than 3 % slower, hussamp 0.67x, slowlatch 0.85x, todd3
0.87x.

**Built**: `Transient.early_abort`'s default 1.25.  Gate: suite 5 329 passed.  The census on all
27 gold decks, three paths: the same run on each (counts and gold
metrics), every deck as the Python sizing measured it -- 78 925 steps,
21 935 rejections, 205 231 Newton iterations (-5.79 %), 0 forced, 16 258
attempts abandoned (were 9 302); no gold column across its pass.  It found
the C loop counting one early check more than the Python on rich3: an
attempt whose check was made and whose screen then fired was counted in
C and again by the Python attempt that made it again -- loop.c now counts
an attempt's checks where it is made or abandoned (`attempt_made`,
`early_reject`; `test_loop_c`'s adaptive run with the screen made to
fire).  One test cut with the abort off (as three were at 132): the
folds' walk (`test_shooting_events`), its vdP grids measured before the
abort -- at 1.25 a retry leaves one accepted step 2.1x under the one
before, and 286 / 287 points where there were 273 / 272.

## 2026-10-15 — 135. SPICE's noncon for the stateful limiter: a Newton iteration on which `elements.Diode` limited is not converged -- a diode clamp with a parallel path no longer returns the bare divider

**The defect** (found 2026-10-14 on the way to the gmin build; plan
`dc-final-gmin.md`, D-b).  5 V - 1 k - b, a diode b to ground and 1 k
from b to ground: `DC()` returned v(b) = 2.500000, the bare divider, where
the clamp holds 0.6297 V (PCNR and the SPICE importer's state-free diode);
a fresh evaluation at 2.5 V reads a KCL residual of 9.5e28 A.
`elements.Diode` is the one STATEFUL limiter: `limit` keeps the limited
junction voltage `_vlim` on the instance and returns None, and `i` / `G`
read the device as the tangent at `_vlim`.  `StandardNewton` accepted on
|dx| and |F|, and `F` was the circuit with that tangent at the clamped
voltage (0.24 V when it returned, iteration 2), which carries nothing: the
divider's first iterate does not move, both tests hold.  Present at every
commit tried back to 64d7ac19 (2026-08-21), C and Python paths alike; the
plain clamp escaped only because its source current is ~0 (the relative
update test failed on it).  Reach: of 36 clamps (V1 0.5..20 V x R2 none,
100 ohm..1 Mohm) 14 returned the bare divider (V1 10 V, R2 100 ohm:
0.909 V); a fixed-step transient through the clamp passed false points at
the source's edge on every path -- Newton, the chord (`chord_jacobian`),
Fang's coupled step -- 0.25 .. 11.8 V off, one chord run ENDING at
9.84 V.  The radau stage Newton already asked its limiters
(`_limiters_at_rest`); its docstring's claim that the ordinary Newton was
saved by its residual test was wrong and is corrected.

**Sized before building** (a pytest plugin patching the solvers' own
source, nothing in the tree): the full suite with the veto, 107 failures
-- the SAME 107 in an A/A control (the same patched source, no veto), so
every one the patching's own (pinned pairs reading the source, the C paths
declining a patched solver); no test result moved.  The veto met 15
acceptances in 12 tests, all `StandardNewton`'s: every one FALSE, a fresh
residual of 254 A .. 3.4e5 A (9e66 A on the loop-gain test, NaN on the
minstep test).  Under the veto each of those solves converged true -- the
acceptance after it within 2e-12 .. 2e-9 A, at most 4e-5 of the solver's
tolerance -- or failed where its test expects a failure (the minstep
test; gate 1.3, an absurd bypass tolerance that must raise, now raising
from the inner DC its docstring names).  So those 12 tests had passed on
false operating points: none reads the operating point's KCL.  The chord
never met the veto in the suite; the gold decks have no `elements.Diode`.

**Built.**  `_limiting.NONCON`, SPICE's `CKTnoncon`: `Diode.limit` adds one
where `_pnjlim` moved the junction (`_pnjlim` hands `vnew` itself back
wherever it does not limit, so a new object is a limited junction; not
under an autodiff toolkit, where `i` / `G` read the node voltage).
`StandardNewton`, `DampedNewton`, `ChordNewton`, `SchurCoupledNewton` and
Fang's coupled step (its non-PCNR branch; PCNR's has `lim_converged`) read
it around their limiter call (`limited_since`) and do not accept an
iteration that moved it.  `NONCON_CHECK` switches it off, for the tests
that prove the check alive.  No C twin: the transient's C Newton declines
every circuit with a stateful limiter (`newton_c:stateful`), the DC's
(`dcnewton_c:unservable`) -- the two pinned pairs say so; their references
re-recorded, the twins unchanged.  State-free limiters never count: they
move the node, and the update the test reads is the limited one.

**Tried, not built**: counting also an evaluation OFF the node voltage
(`Diode.i` at a `_vlim` the previous limiter call left -- SPICE's own
reading of noncon, which counts the load's evaluation): it did not move
the one discrepancy left.

**Known gap, open**: the chord on this diode reads 1.8 .. 2.4 mV off the
state-free diode on the three clamp transients, where the chord on the
state-free diode is within 5e-6 V of the full Newton and the full Newton
on this diode within 2e-8 V; `test_the_chord_through_the_clamp_...`
bounds it at 1e-2 V, which separates it from the defect.

**Gate**: suite 5 393 passed, 6 skipped, 3 xfailed (5 336 + the 57 new).  The gold census on all 27 decks identical to f2a500d5's on every deck (steps, rejections, Newton iterations, forced steps, the operating point's and the waveform's digests, the gold metrics: 78 925 steps, 205 231 Newton iterations, 0 forced).  `test_noncon.py` (57): the clamp
DC by both Newtons, the 36-clamp grid by a fresh residual, the three
transients on Newton / Fang (2e-8 V of the state-free diode) and chord,
each solver refusing a limited iteration on a bumped count, the count's
owners (a state-free limiter and an autodiff toolkit never count), and
the check switched off bringing every face of the defect back.

## 2026-10-15 — 136. The chord's settled test: a chord iterate whose increments contract monotonically is not accepted until its estimated error is under the tolerance -- the diode clamp's 2 mV gone, slowlatch's 21-tolerance early stops gone, every oscillating chord as before

**The gap 135 left.**  After noncon the chord (`chord_jacobian`, 'auto' --
on by default where compiled device Jacobians are expensive: 22 of the 27
gold decks) still read 1.8 .. 2.4 mV off the state-free diode on the clamp
transients, where the full Newton on the same diode was within 2e-8 V and
the chord on the state-free diode within 5e-6 V of the full Newton.  The
cause: a few solves at the source's edge whose increments contracted at
0.97 a step -- the held Jacobian taken with the junction at a clamped
voltage, its conductance far under the root's -- crawled 80 .. 90
iterations until an increment passed the update test (|dx| under its
tolerance), and stopped ~30 tolerances from the root: a linearly
convergent iteration's error after an increment dx is lambda/(1-lambda)
|dx|.  Counting also the diode's reads off its node voltage (SPICE's own
reading of noncon) did not move it.

**Sized against lines written first** (plan `dc-final-gmin.md`, 'A1
SIZING' and 'A1 REFINEMENT'; a source patch of `ChordNewton`, the C paths
off so the Python chord saw every solve).  A1 -- refuse where |dx| x
theta/(1-theta) reaches the tolerance, theta the ratio of the increments'
max norms -- refused 20 acceptances on 7 gold decks: no gold crossing, no
forced step, Newton -0.03 %; but the truth (each step re-solved from the
accepted point by the full Newton to 1e-3 of the tolerances; the chord
itself continued lands on that same root, 0.00 tolerances away, every case)
put only 8 of them early (slowlatch 2.3 .. 21 tolerances, todd3 1.85, arom
1.01) and 12 within 0.26 .. 0.41 tolerances where A1 estimated 1.04 ..
2.9.  Over 56 689 multi-increment acceptances with their truth, no
magnitude refinement separated them (weighted / RMS norms, two-ratio
means; a theta threshold 0.7 did, cutting at slowlatch's and todd3's 0.707
.. 0.710 -- a fitted constant).  THE SIGN did, as predicted before it was
measured: every one of the 12 oscillates (lambda -0.55 .. -0.78, error
|lambda|/(1+|lambda|) |dx|, which reproduces each truth to 2 decimals),
every early stop is monotone (+0.57 .. +0.97, lambda/(1-lambda) |dx| = the
truth).  A2 -- the signed ratio at the last increment's dominant weighted
component -- refuses the 11 monotone ones and none of the 12; gold census
in its mode: 24 of 27 decks bit-identical, arom +1.18 %, slowlatch -4.17 %,
todd3 0 Newton iterations, no gold crossing, no forced step, +0.008 % in
total.  Beyond any ratio rule: one rich3 acceptance 1.81 tolerances off at
theta 0.05.

**Built.**  `nrsolver._chord_settled`: lambda = this increment over the
last at the last one's dominant weighted component (`np.argmax`, the
first); settled under 1/2 (the oscillating estimate is then under half a
tolerance at any accepted increment), never at 1 or above, else where max
|dx|/tol x lambda/(1-lambda) is under 1.  `ChordNewton` asks it of an
iterate its update and residual tests accept, after its second increment.
Its twin is newton.c's `chord_settled` / `chord_dominant` (the same
operations in the same order); `_tran_newton_c`'s chain holds
`_chord_settled`, so a patched rule takes the Python Newton.  The pinned
pair 'the Newton solve and its C' carries it, both sides re-recorded.

**Cost seen**: a refused monotone acceptance keeps the chord crawling until
it settles (the clamp 16 more iterations; a falling-edge clamp of
state-free diodes +9 % Newton iterations, 437 vs 400) -- today's chord
already crawls 80-90 iterations on those solves before its early stop;
handing a slow chord to the full Newton is sized next.

**Gate**: Suite 5 404 passed, 6 skipped, 3 xfailed (5 393 + the 11 new).  The gold census on all 27 decks on the three paths (Python loop and Newton, Python loop and C Newton, the C loop) each identical to the A2 sizing on every deck (steps, rejections, Newton iterations, forced steps, the waveform's digest) and to 135's but on arom, slowlatch and todd3; no gold crossing; 0 forced; 205 247 Newton iterations on each path (135: 205 231).  `test_chord_settled.py` (11): a
scalar chord with a held slope 10 (monotone 0.9: the update test alone
stops ~9 tolerances short, the settled one within one) and 0.625
(oscillating -0.6: the same iterate as the update test), the rule's
signed-ratio table, and on a falling-edge clamp of state-free diodes the C
chord byte for byte the Python's, the rule acting there (off, the answer
moves), a patched rule taking the Python Newton.  The clamp chord test of
135 tightened from 1e-2 to 1e-6 V against the state-free diode's chord.

## 2026-10-15 — 137. No rounding-level final step: the run's end merges a source corner within `minbreak`, and a step ending within it short of the end ends at it

**Found** diagnosing the gmin route's slow small decks (plan `dc-final-gmin.md`, '(2) DIAGNOSED'). With the standing
gmin inside the junction devices, hussamp, nand and schmitecl each ran one transient Newton to 100 iterations per run
before a rescue: their LAST step.  hussamp's `.tran` stop and its pulse's corner round one ulp apart (the end
1.0000000000000001e-07, the corner 1e-07 -- a stop of 100n read as 100 x 1e-9 against the card's 100E-9; nand the
same): the loop landed on the corner and then stepped h = 1.3e-23 s to the end; schmitecl's 5 ns capped steps summed
to 4.999999999999999e-07 against an end of 5.000000000000001e-07 (500 x 1e-9) and left a 2.1e-22 s step.  At
that step the companion conductances are ~1e23 and J is numerically singular (singular values 7.8e11 .. 4.7e-11,
the near-null modes the sources' branch currents): the solve returned noise (|J dx + F| = |F|), the Newton
repeated a useless step.  Without the gmin the same step happened to converge.  `_next_breakpoint` merged two
corners within `minbreak` of each other, but never a corner with the run's end; and an adaptive step had no residue
rule (the fixed grid's: a residue under 1e-9 of its step ends the run).

**Built.**  SPICE keeps its stop time in the breakpoint table, so a corner that close merges into it: here
`_StepEvents._next_breakpoint(t, tend)` returns `tend` for a corner within `minbreak` of it (relative to `tend`), the
three step families passing `run.tend`; loop.c's `next_breakpoint` the same with `s->tend` (its twin, the same
operations; the pinned pair 'the stepping loop and its C' re-recorded).  The step that would land on the corner
lands on the end, and the run's last time is `tend` exactly.  And in `_SteppingLoop.where` (loop.c's `where`, its
twin), adaptive runs only: a step ending within `minbreak` short of the end ends AT it (`h = tend - t`, exact
there: Sterbenz), and a residue under `minbreak` left at the end ends the run -- the adaptive twin of the grid's
rule.  A fixed grid keeps its own (its landing flag reads the merged time).  The pairs 'the stepping loop and its C'
and '... its unread Jacobian' (reference; a fixed run's pieces unchanged) re-recorded.

**Gate**: Suite 5 412 passed, 6 skipped, 3 xfailed (5 404 + the 8 new).  The gold census on the three paths, each against the A2 commit's: hussamp and nand (the corner merge), schmitecl and mosrect (the near-miss rule) one step fewer each (Newton -6 / -8 / -4 / -1), their gold metrics within 3e-16, the other 23 decks identical, every path the same, no crossing, 0 forced.  `test_final_step.py` (8): an end one ulp past a corner landed on in one
step with no rounding-level step (red before), an end one ulp before a corner, the C loop byte for byte the
Python's with the C serving, the merge's reach (only within `minbreak`, either side; without `tend` as before), and
`where` on a stub: a step two ulps short of the end ends at it, a one-ulp residue ends the run (both red on the old
loop), a step 1e-13 s short and a fixed grid as before.  (A natural reproducer of the accumulation was not found
outside schmitecl: the controller's ramp and a step cap tied to `tend` break the uniform sum.)

## 2026-10-15 — 138. The lift's remainder joined to the nodes it physically joins: a series resistance no longer measured from the gate

**Found** sizing the gmin route's ladder (plan `dc-final-gmin.md`, 'gm1 DIAGNOSED'). gm1's operating point never
converged on PCNR: every iteration passed the residual and limited tests and failed the step test, |dx| 2e-5 to 6e-5 V
against ~1e-6 V, on nodes of a latch whose supply is 0 V at DC -- a cluster held only by cut-off channels and junction
leakage (the reduced Jacobian's smallest singular values 5e-13 to 6.7e-13).  The step scaled as 1/gmin (from the same
start, 2e-6 V at 1e-11, 1.7e-7 V at 1e-10): a fixed 2e-17 A of noise along the cluster's mode.  The noise was in the
raw KCL sum: on the internal drain and source rows of the ON transistors (gate at 5 V) the view's remainder carried
+-1.05 A and the devices' probe currents -+1.05 A, cancelling to what flows and leaving ulp(1.05) = 2.2e-16 A.  The
remainder of `MosLevel3GateChargeHdl` stamped rd between the drain and the GATE column: `_affine_in_nodes` splits a
current after the probes are substituted, when x_di is the probe tree's root plus probe voltages (x_g - vgs + vds,
the root the first probe's plus terminal), so the lift kept G (x_d - x_g) and the probe current G (vgs - vds) -- each
0.21 S x 5 V, exact in sum (every answer right), ulp(1 A) of KCL lost.  `_kernel_map` followed it, mirroring d and s
onto the gate.  (The same root bit the first version of the lift, sec. 40: the internal drain out by exactly the gate
voltage.)  The element-gmin prototype's 'solve' of gm1 was a lucky step: the same jitter, one step 2.4e-8 V at
iteration 15.

**Built.**  `hdl._paired_lift`, called after the lift in `_pcnr_declared_route`: the remainder completed from the
columns of the unprobed nodes, which carry the physical entries (column q as lifted, row q its transpose, each
reached row's diagonal its zero row sum), and what leaves it folded into the rest, `rest += (L - K) pot` -- probe-only,
since the difference's rows sum to zero and the root cancels.  Node unknowns only: a V-contributed branch's stamp
(`LedHdl`'s) is no conductance and is kept as lifted; where a rest would still read a node voltage the lift's own
split is kept.  `pcnr_lin_G` and `pcnr_i`/`pcnr_didv` change together; the device's total is unchanged by
construction; the compile cache keys on hdl.py's source.  Re-paired: MOS levels 1, 2, 3 (NMOS, PMOS, gate charge) and
the Statz MESFET -- rd d-di, rs s-si, `_kernel_map` mirroring d onto di and s onto si.  PCNR is off by default in DC
and Transient, so no gold deck's production run calls either.

**Gate**: Suite 5 442 passed, 6 skipped, 3 xfailed (5 412 + the 30 new).  The PCNR operating point of all 27 gold decks (DC pcnr=True, the continued solve, no gmin;
HEAD's tree against the fix, Python loop): PCNR iterations 4 343 -> 2 658 (-39 %) -- gm1 1 288 and a fallback -> 7
PLAIN, mux8 200 (gshunt) -> 66 plain, slowlatch 126 (gshunt) -> 48 plain, g1310 126 -> 35, e1480 564 -> 468, opampal
277 -> 233; up on toronto 60 -> 73, voter25 466 -> 495; fallbacks 2 -> 1 (opampal); answers within 1 mV, most within
1e-9 V, but slowlatch (17 nodes, bistable) and e1480's one floating group (78 nodes moving together 0.63 -> 8.89 V:
both points solve the circuit, max|f| 4e-14 A, 99.93 % of the difference along a mode of sigma 3.6e-14 -- without a
gmin nothing fixes it).  The gold census on the three paths against HEAD's: identical -- no deck moved, no crossing, 0 forced (PCNR off in production, as it must be).  `test_pcnr_affine.py` +4: the lift serves the series-resistance models
(the control); every lifted class's remainder joins each terminal only to its own internal node, rows summing to zero
(red before on 13 classes); remainder + probe current is the element at drawn points; at a 5 V gate the remainder
carries what flows (red before: +-1.000004 A).  (A one-transistor cluster held by 1e-12 S converges either way --
each device's +-1 A lands as an equal and opposite pair on rows the stiff cluster's mode weights alike -- so the stall
itself is gm1's, rows shared by several devices summing different roundings.)

## 2026-10-16 — 139. PCNR's junction-gmin term and the continued solve's route: true gmin stepping, nested in source stepping where a step fails

**Why** (plan `dc-final-gmin.md`, from 'C2's NO-REBUILD KNOB' to 'ARM N MEETS ALL FOUR LINES').  The gmin build needs
the SPICE decks' operating points on PCNR, and PCNR's continued solve -- plain, a gshunt ladder, source steps of plain
solves -- left two gold decks to the ordinary chain: jge (with the standing gmin, its gshunt ladder's last rungs fail
on modes the gmin holds) and opampal (148 Gummel-Poon BJTs from a 35 V supply: PCNR's equations right -- at the
chain's point its residual 2e-14 A and one iteration from there -- but from a cold start every ladder's first rung
ran away, junctions predicted at 34-736 V and cut back by volts each step).  The chain itself solves opampal only by
source stepping with node-gmin AND junction-gmin ladders nested in each step.  Sized on all 27 decks, in order: true
gmin stepping by device-parameter writes rebuilt PCNR's record every rung (~55 ms of a 63 ms rung); as a term in the
equations instead (the 'knob') FIRST, it costs g1310 2.09x (a hard stretch of g1310's junction-gmin path, 2e-4 to
6e-4 S: failed rungs, each a solver call) -- refused, and 100x rungs did not help; AFTER gshunt (C1) it missed its time
line by 0.19 %, all opampal's failing ladder; nesting the junction ladder in the source steps (the chain's shape) PCNR
solves opampal itself (1.99 s against the chain's 5.41) and that route met every line: answers as before, every gold
operating point on PCNR, the operating points 52 % faster in total (measured with device gmin and the held-node
prototype, in a worktree).

**Built.**  `pcnr.solve_dc(jgmin=g)`: `g` across every junction pair (`pcnr_junction_pairs`) IN the equations --
`F[a] += g vj`, `F[b] -= g vj`, J's four entries, after gshunt's rows -- in the Python loop and in `pcnrdc.c`'s
`assemble`, its twin (the pairs' rows set once per record; the pinned pair 'the PCNR DC solve and its C' re-recorded,
its reason naming the term).  `solve_dc_continued`: plain -> the gshunt ladder and its pure solve -> `_jgmin_ladder`
(1e-2 down by decades to 1e-11, a failed rung retried as gshunt's, then the pure solve) -> `_source_steps`, a step
that fails retried by `_source_step`: the gshunt ladder and then the junction-gmin ladder at its source level.  A
gshunt ladder that came down with its pure solve failing no longer ends the route: the junction ladder deforms the
junctions, not the nodes.  Behind DC(pcnr=True): no production default moves (held nodes in PCNR are the next commit).

**Gate**: Suite 5 448 passed, 6 skipped, 3 xfailed (5 442 + the 6 new).  The 27 gold decks' PCNR operating points against 8e8d33c0 (no gmin, the Python loop): every deck but opampal bit-identical, answer and iterations; opampal from the chain's fallback to PCNR's own nested source steps (2 271 its), within 2.8e-14 V of the chain's answer; fallbacks 1 -> 0.  The gold census on the three paths identical (PCNR off in production).  `test_pcnr_jgmin.py` (6): the term in C bit for bit the Python's
(two conductances, the C serving); the circuit with explicit conductances across the junction pairs has the term's
solution as a fixed point (one iteration, 1e-12 relative); the route's order on a stand-in solve -- the junction
ladder after gshunt, by decades from 1e-2 to 1e-11 with its pure solve; a failed source step retried by both ladders
at its level; the plain solve alone.  `test_dc_pcnr.py`'s two forced-failure tests now pin the new route (the gshunt
ladder's eleven decades, then the junction ladder's ten; every solve failing: each of the five source levels tried
plain and by both ladders' first rungs).

## 2026-10-16 — 140. Held nodes in PCNR, eliminated as the chain eliminates them: DC(pcnr=True) with .ic / .nodeset

**Why.**  `DC(pcnr=True)` with `pin` or `nodeset` raised ("not supported"), so a SPICE deck's `.ic` or `.nodeset`
kept its operating point off PCNR -- the gmin build's sizing held such nodes through temporary voltage sources in a
harness (plan `dc-final-gmin.md`, S1).  The ordinary chain eliminates a held node (`_Held`, roadmap's stage 3 of the
SPICE benchmark plan): not an unknown, its KCL row and column out, its value `lam` times its volts.

**Built.**  `pcnr.solve_dc(held={row: volts})`: `x[held] = source * volts`, the held rows left out of the reduced
Schur solve with the reference's (`predict`'s `held`: their step 0) and out of the step and residual tests (the hold
supplies their current).  In C a row map replaces the single reference (`pcnrdc.h`'s `m` and `red`, set per call by
the glue: `solve_dense`, `solve_klu`, `solve` and `converged` read it; without a hold it is the old arithmetic's
positions exactly) and the KLU choice is made on the same reduction.  `solve_dc_continued(held=, v_lim0=)`.  `DC`:
with `pcnr`, the chain's refusal first (a node a voltage source or an inductor holds), then PCNR -- a nodeset held
and released from that answer (but for the pins), a failed held solve leaving the hint unused as the chain's does;
the RELEASE starts from the held solve's limited unknowns: re-seeded from x (clamped, as a cold start is), the
cross-coupled latch fell from its picked state to its metastable point (q = qb = 2.56 V).  A PCNR failure falls
back to the chain as before.  Pinned pairs re-recorded: 'the PCNR DC solve and its C' (both sides; the reduction in
its reason) and 'the DC base solve and its C' (reference only: `DC._solve_dc`'s PCNR branch; its C twin's tests pass).

**Gate**: Suite 5 455 passed, 6 skipped, 3 xfailed (5 448 + 7).  Every gold operating point without a hold, production DC(pcnr=True) on the C path, bit-identical to 99af7e86's (answer, path, the C's iterations).  The six gold decks with .ic / .nodeset (arom, gm17, hussamp, jge, rich3, ring), which raised before, solved by DC(pcnr=True) on PCNR, no fallback; against the chain's held operating points within 1.3e-6 V but jge -- its held solves agree off its 18 floating modes (both solve the circuit to 2e-17 A; 99.93 % of a 44 mV difference along them) and the release lands in another valid state of a multistable deck (neither answer keeps every nodeset).  The gold census on the three paths identical (PCNR off in production).  Tests: `test_dc_held.py` -- the refusal test replaced
by PCNR's holds against the chain's on the latch (a pin, each nodeset state, a pin with a nodeset: the chain's answer
to 0.0 relative) and the same source-hold refusal; `test_pcnr_jgmin.py` -- a held solve in C bit for bit the
Python's, the C serving (full sources, half -- the hold scales -- and with the junction term).

## 2026-10-16 — 141. SPICE's junction gmin inside the device models: `gmin` on MOS levels 1-3 and the SPICE diode

**Why** (plan `dc-final-gmin.md`, S2).  SPICE puts `CKTgmin` (1e-12 S) across every junction inside the device
loads -- `mos*load.c`'s `CKTgmin*vbs` and `CKTgmin*vbd` on the bulk junctions, `dioload.c`'s `CKTgmin*vd` -- and
our models carried none, so a node held only by junctions floats (e1480's group of 78 nodes, gm1's unpowered latch).
Sized in the gmin plan, the standing gmin inside the devices costs nothing per step where elements across the pairs
cost 1-4 % (MOS3 decks).  The Gummel-Poon BJTs had it already.

**Built.**  A `gmin` parameter, SPICE's semantics, default 0 (the card's physics; the importer will set SPICE's
1e-12): MOS levels 1, 2 and 3, NMOS and PMOS, plain and gate charge (12 classes), `gmin*vbs` and `gmin*vbd` on the two
bulk junctions; the SPICE diode (with its thermal, LED and photodiode variants) a contribution of its own across the
junction, `gmin*vd` at the limited voltage, not area-scaled -- NOT inside the junction current `id`, which also feeds
the shot noise, the dissipation and the LED's light (inside it, the LED emitted light from the gmin leakage: the
first test run caught it).  The importer is untouched: no gold answer moves.

**Gate**: Suite 5 474 passed, 6 skipped, 3 xfailed (5 455 + the 19 new).  The gold census on the three paths identical.  `test_hdl_params.py`'s records: the diodes'
explain digests re-recorded (the parameter and the contribution added), their numeric POINT digests unmoved -- the
compile at gmin = 0 bit for bit as before.  `test_device_gmin.py` (19): the 16 carriers enumerated (the control);
`G(g) - G(0)` exactly `g` across each junction pair (`pcnr_junction_pairs`) and `i(g) - i(0)` that matrix times the
potentials, at drawn points, every carrier; the class C kernels with the gmin on bit for bit the Python's (PCNR's
solve on a latch, the C serving, at 1e-12 and 1e-6).

## 2026-10-16 — 142. The operating point by PCNR by default, and SPICE's junction gmin on every imported junction

**Why** (plan `dc-final-gmin.md`, S2 commit B).  With roadmap 139-141 PCNR solves every gold deck's operating point on
its own (true gmin stepping, nested source steps, held nodes) and the device models carry SPICE's junction gmin; this
commit makes them the default: `DC` takes PCNR wherever a device declares a probe, and the importer builds every
junction element with SPICE's `CKTgmin`.

**Built.**  `DC`'s `pcnr` default is None -- AUTO: PCNR where a device takes part, the ordinary chain silently where none
does; True still asks for it (and reports 'no-participants'), False is the chain.  The transient's operating point is a
`DC`, so it follows.  The importer: `Imported.gmin`, SPICE's 1e-12 unless `.options gmin` (ngspice's ungrouped, Xyce's
`.options device gmin`) sets it, given to every element that takes `gmin` (MOS 1-3, the SPICE diode family, the
Gummel-Poon BJTs); `import_netlist(gmin=)` overrides both (0, the cards' physics alone); `write_ngspice` writes a
non-default.  `pcnr.solve_dc` builds its view FIRST: adding the elements to it notifies their parameters, which moves the
epoch, and a transient's stamp plan, built earlier in its operating point, was rebuilt once after it (two plan tests).

**The suite's moved tests, each with its cause.**  Fifteen tests OF THE CHAIN now took PCNR's operating point: a conftest
fixture `ordinary_dc` (DC's default False for the test) keeps them on it -- the DC Newton C twin and the gmin anchor
(module-wide), noncon's check-off, gate 1.3's bypasstol, the DC sweep's continuation count, two FET-limiting counts,
gate 13.7 (its operating point assumed the classic path) -- and explicit `pcnr=False` on `test_pcnr_affine`'s five
ordinary-path references and the status test.  The DC twin's decks import with `gmin=0`: with SPICE's gmin mike2 and gm3
converge on the plain solve at once and leave no ladder for the C to serve.  The level-2 channel test imports with
`gmin=0` (the junction gmin puts 1e-12 S x vbd = 2e-12 A on the drain).  The BJT excess-phase test's operating-point
agreement 1e-12 -> 1e-10 (PCNR solves the two circuits -- PTF adds unknowns -- to 8.4e-12 of each other).  The SPICE
diode's rs test asks the solve for reltol 1e-12 and asserts 1e-12: its old 1e-8 "logarithmic shift" was solver
tolerance (PCNR at the default reltol 1.3e-7 off, the chain 6.7e-10, asked 2.2e-15) -- DC's default answer is accurate
to its tolerances, no longer to the chain's quadratic luck.  The imported-deck test builds its hand circuit with gmin.

**Not built: a bound on PCNR's node steps.**  From a cold start PCNR's iterates swing far (slowlatch 1.7e5 V, opampal
1.3e4) where the chain's stay near +-50 V: the chain limits the node vector itself, PCNR only each device's unknowns.
Three bounds were sized (ngspice's node damping, the update scaled uniformly to the largest source, each node clipped to
+-max(U, |x|)): each cut the gold decks' swings and work (the last -17 % iterations, four decks newly plain) and each
broke what PCNR's hard-start tests pin -- an EKV diff pair's convergence (unbounded it returns from -5.7e5 V in a step:
the devices see their limited unknowns, so the excursion is harmless; bounded it cycled), order independence, parity with
plain Newton, and an LED's optical node (nodes are untyped: 3.8e18 in its fixture).  Recorded in the plan.

**Gate**: Suite 5 478 passed, 6 skipped, 3 xfailed.  Every gold operating point on PCNR (production `DC`), 4049osc failing there as before.  The gold
census on the three paths against 97089a2d, identical to each other: every deck moves (gmin, PCNR's operating point),
0 forced; e1480's v(8) 0.0084 -> 0.0001 (its floating group pinned by gmin), gm17 v(5) 96 -> 28, toronto and mux8
better, gm1 marginally worse on columns it already failed, jge +8 % steps; two latches with no .ic / .nodeset land in
their OTHER exact state (both points solve the circuit to 1e-17 A; a latch's state is the solver's early dynamics --
same first iterate, then swings): todd3 fails two columns (i(vseu) 0.68 -> 3.21, i(vamp) 0.25 -> 2.74), slowlatch passes
two it failed (v(4) 32 -> 0.73, v(5) 43 -> 0.62).  The gold decks' transient solve time (production, the C path, pinned malloc, median of runs 2-4) against 97089a2d: 81.52 -> 79.42 s in total (-2.6 %), no deck over max(1.10x, +0.5 s) -- opampal 5.80 -> 2.83 s, g1310 0.32 -> 0.19, gm3 0.31 -> 0.16, gm17 0.41 -> 0.26, gm19 6.87 -> 6.22; jge 27.6 -> 29.3 (its +8 % steps), e1480 0.21 -> 0.36 and todd3 0.08 -> 0.17 (within +0.5 s).  New tests: the importer's gmin (default, ngspice,
Xyce), DC's default (PCNR where a device takes part, silently the chain where none, reported when asked).

## 2026-10-16 — 143. PCNR's source steps: the junction-gmin ladder first in a failed step, and the first step 0.05

**Why** (plan `dc-final-gmin.md`, 'PCNR SPEED').  Where PCNR's operating-point work goes, every rung of every gold deck's
operating point logged (the Python loop): of 5 461 iterations, 3 668 (67 %) are spent in rungs that fail and 3 030 (55 %)
in rungs that run to the 100-iteration limit.  Stopping a failing rung sooner is not safe: successful rungs run to 98
iterations (nand's plain solve), and wander first -- up to 81 iterations without a new minimum of max|f| before
converging (gm19's gshunt rung at 1e-4 goes from 9e-3 to 4e22 and converges at 61) -- so every 'no new minimum in N'
rule below ~82 kills successful rungs (N 20: 17 of 149; N 50: 8, jge's and mux8's plain solves among them).  The waste
that can go is structural and nearly all opampal's (2 571 of 5 461), the one gold deck whose operating point steps its
sources: its first step 0.1 failed on every try (~300 iterations), and at 0.05 the gshunt ladder spent ~900 in a dead end
(its 1e-9 rung and eight retries at the limit) before the junction-gmin ladder succeeded there in ~150.

**Built.**  `_source_step`: a failed source step tries the junction-gmin ladder before the gshunt ladder (the route's own
order at full sources is unchanged: plain, gshunt, junction gmin); `SOURCE_STEP` 0.1 -> 0.05 (the steps before giving up
0.05, 0.025, 0.0125).  Sized first, with its lines written before measuring (both alone and together): opampal's PCNR
iterations 2 522 -> 1 563 (order) / 2 225 (start) / 1 266 (both); its answer within 3.7e-12 V; the 25 other decks that
solve identical (path, iterations, the operating point's bytes); its operating point on the C path 1.898 -> 0.992 s
(median of rounds 2-4, interleaved).  Tests: the nested-step route test now asks for the junction ladder first and, where
it fails, the gshunt ladder after it at the same level; the give-up test for three steps from 0.05.

**Gate**: suite 5 478 passed, 6 skipped, 3 xfailed (as commit B).  Every gold operating point by production `DC`
identical to commit B's but opampal's (C-path PCNR iterations 2 552 -> 1 286, the answer within 3.7e-12 V); 26 of 27 on
PCNR, 4049osc failing there as before.  The gold census on the three paths against commit B's: only opampal moves -- the
same steps, Newton iterations and 0 forced, its metric by 1.6e-14 -- no crossing.  The 27 decks' runs against commit B,
quiet (1-min load <= 2.7), alternating trees per deck: opampal 2.828 -> 1.925 s (0.68x; the operating point's 0.9 s),
the total 79.30 -> 78.35 s (-1.2 %), every other deck within +-3 % on identical steps (code they do not run).

## 2026-10-16 — 144. The full Newton's sparse solve in C, and the transient's KLU from 40 unknowns

**Why** (plan `dc-final-gmin.md`, 'TRANSIENT SPEED', 'SPARSE FULL NEWTON IN C').  Where the big decks' time goes (perf,
one BLAS thread, leaf samples): jge (244 unknowns), mux8 (171) and voter25 (200) spend 34-43 % in dense LAPACK -- below
`MIN_N_FOR_SPARSE` (250), set 2026-10-02 on the Python loop with SuperLU, before the C stepping loop and KLU on the
circuit's structure (2026-10-10).  On today's C loop KLU on the structure makes every gold deck's transient from 40
unknowns faster (mike2 (40) 0.97x .. voter25 0.53x, jge 0.64x; three rounds interleaved), the 27 decks' transients 77.0
-> ~52.6 s, no gold column crossing.  But the C Newton served the sparse solve on the chord alone: a run on the FULL
Newton (chord 'auto' is off for cheap devices) was declined (`newton_c:sparse_full`) and its whole transient ran the
Python loop -- the small decks that run it on KLU 4-10x slower, and every full-Newton circuit of 250+ unknowns on the
Python loop today.  So the floor could not come down for every transient until the full Newton's sparse solve was in C.
Sized first (its lines in the plan): on each deck's own Jacobians, warm, the full Newton's solve 6.8 -> 2.8 us at 40
unknowns, 38.5 -> 6.5 at 97, 692 -> 25 at 243 (dense LU goes threaded from ~100 unknowns: 5x its cubic trend); the
forced-full-Newton transients predicted 0.93 .. 0.32x, the floor 40.

**Built.**  newton.c's `sparse_full`: each iteration's solve is `StructuralKLU.solve` on its refactor path -- klu_refactor
of the run's numeric, klu_solve, the validation in numpy's operations and order (`arr @ x` is numpy's matmul's dgemv, the
NaN-propagating maximum, `max(max |b|, 1e-300)`, `not (resid <= tol)`) -- as lte.c's twin of the same path.  Every other
course is handed back and the Python's repeats it from the same seed: no numeric yet, a non-finite matrix or right-hand
side, a zero pivot, a failed solve or validation, a nonzero off the pattern -- exact, because a refactor writes the
numeric's values, never its pivot order.  A core with no pattern reduces J whole and makes `_ours`'s count.  The glue
gives each call the run's numeric (`_record`; the C loop before its call) and counts a refactor an iteration on a serve.
The floor: `linearsolver.MIN_N_FOR_STRUCTURAL_KLU = 40`, `AutoSolver(structural_min_n=)`; the transient asks for it on
the multistep methods the C loop serves (Gear-2 -- the default -- Euler, trapezoidal, theta); Radau, the stage methods
and the GLMs keep 250 (not measured), as do DC and every other analysis, a matrix the structure does not cover and a run
without one.

**Checked.**  The C against the Python (C loop and C Newton off), the full Newton forced, KLU from 0: 21 decks byte for
byte -- waveform, times, steps, rejections, Newton iterations, the StructuralKLU's counts (voter25 through a maxiter
hand-back and two loop hand-backs).  Tests: the full Newton on a sparse ladder byte for byte the Python's and served in
C; every validation refused (tolerance 0) and every refactor refused in C (a stand-in), each handed back with the run the
Python's; a core with no pattern.  The byte test is alive: a one-ulp perturbation of the accepted step fails it (a
skipped refactor is caught by the validation and handed back -- the bytes hold).  `test_linearsolver`'s historical-path
test sized below the new floor.  Pinned pairs: the Newton pair gains `StructuralKLU.solve` and `_ours`, `_record` and
`_serve`; five re-recorded.

**Gate**: suite 5 482 passed, 6 skipped, 3 xfailed (5 478 + the four new).  Every gold operating point identical to
HEAD's (DC's floor unchanged; 26 of 27 on PCNR, 4049osc as before).  The gold census on the three paths, identical to
each other, against HEAD: the 16 decks of 40-249 unknowns move (KLU's last bits: mux8 17 000 -> 17 356 steps, voter25
9 971 -> 9 924, jge 10 745 -> 10 692, the largest metric change toronto's 0.018), the 11 others identical; no column
crosses, 0 forced.  The 27 decks' runs against HEAD, quiet, alternating trees per deck: 78.27 -> 55.12 s (0.70x) --
voter25 0.56x, gm1 0.60x, mux8 0.62x, gm17 0.62x, toronto 0.65x, jge 0.66x, latch 0.69x, ab_integ 0.74x -- the decks
under 40 unknowns and from 250 within +-1.2 %.  Forced onto the full Newton (three rounds interleaved, medians), sparse
against dense in the build: every deck from 41 unknowns 0.50 .. 0.93x (mike2 0.93, voter25 0.58, e1480 0.50); a
full-Newton run from 250 unknowns, the C against HEAD's Python loop: rich3 11.56 -> 4.68 s, arom 12.79 -> 4.28, gm19
24.96 -> 9.60, opampal 4.21 -> 1.28.

## 2026-10-16 — 145. A held operating point's solve on the structure without its held rows

**Why** (plan `dc-final-gmin.md`, 'THE HELD-NODE DEFECT').  Found while sizing the transient's sparse floor: a DC
operating point with nodes held (`.ic`, `.nodeset`: PCNR eliminates them, roadmap 140) on a circuit at or above DC's
sparse floor (250 unknowns) handed DC's one `AutoSolver` the held solve's smaller matrix first.  It did not match the
circuit's structure, so the solver chose SuperLU on it -- and kept that choice for the release too -- and PCNR's C
declined every solve (`pcnr_c:linsolver`): the whole operating point ran PCNR's Python loop.  Two gold decks today,
rich3 (269 unknowns, `.ic` on 2 nodes) and arom (295, 3), both converging plain, so it cost little there; with DC's
floor forced to 0, jge's operating point (its ladders) took 31.9 s against 0.68 dense.

**Built.**  `AutoSolver.held(drop)`: the solver of this one's matrices with the rows and columns `drop` removed -- its
floors, and its structure without them, renumbered (the reduced matrix's own pattern: `predict`'s keep, the C's row
map) -- one per set of rows, kept on the parent with its choice and factorisation.  `pcnr.solve_dc` hands a held solve
that solver, to the C and the Python alike; DC's own then chooses on the release's full matrix.

**Checked.**  rich3 and arom: the held solve on StructuralKLU, served by PCNR's C, byte for byte the Python PCNR's on the
same solvers; their operating points move 5.0e-18 and 6.5e-19 V.  A test: a 300-section SPICE-diode ladder with a node
held, the held and the release solves both StructuralKLU (301 and 302 unknowns), PCNR's C serving both, the answer the
Python's byte for byte and the dense solve's to 1e-12 -- on the code before, its solver is SuperLU.  jge with DC's floor
forced to 0: 31.9 -> 0.16 s (the dense solve 0.69).  Pinned pairs: the PCNR DC solve's and the frozen transform's
re-recorded (`AutoSolver`, `pcnr.solve_dc` on their reference side; their twins take the same solver object).

**Gate**: suite 5 483 passed, 6 skipped, 3 xfailed (19:31 under another session's simulations: the slow tests
uniformly 1.4-1.8x).  Every gold operating point identical but rich3's and arom's (now served in C).  The gold census on
the three paths, identical to each other: only rich3 (3 618 -> 3 625 steps) and arom (3 119 -> 3 115) move, their
metrics by at most 1.1e-9; no crossing, 0 forced.

## 2026-10-17 — 146. A square is one multiplication: `_sq(x)`, `x * x`, in the numpy chain and the C alike

**Why** (plan `dc-final-gmin.md`, 'POW IN THE DEVICE MODELS'; Andreas: squares only).  libm `pow` was 13-26 % of the
big MOS decks' time (perf).  ~90 % of the MOS models' `pow` calls take an INTEGER exponent from their symbolic
derivatives, most of them 2 (MOS3's Jacobian pass: 194 of 218).  They were `pow` because the C backend's contract is
bit-identity with the numpy path, whose scalar `x ** 2` is glibc `pow` -- which differs from `x * x` in 0.07 % of
arguments.  So both printers change together.  Sized first (a timing prototype rewriting the kernels' C at load): the
four big decks' transients -26 % on squares alone, -29 % with cubes and fourth powers too -- but `x*x*x` differs from
`pow(x, 3)` on 25 % of arguments (two roundings, where `x * x` is correctly rounded), so squares only.

**Built.**  The numpy chain printer prints a power whose exponent is the number 2 (`2`, `2.0`: sympy's Float(2.0) is not
Integer(2)) as `_sq(base)` -- `_sq_numpy`, `x * x`, registered in `_KERNEL_NUMPY`, numpy's `square` on an array -- and
the C printer, derived from it, inherits the rule (the prelude's `_sq`, `a * a`); its array-square branch is gone.  Every
HDL model's value moves in its last bit on those 0.07 % of squares, toward the correctly rounded one.  Tests: the `pow`
sentinel (`-fno-builtin-pow`) moved to `u ** -1.0` -- gcc folds `pow(x, -1.0)` to `1.0/x`, which glibc's differs from on
44 in 50 000 arguments, and its mutation arm still fires (`u ** -1` itself prints as a division); a square at a double
where `pow` and `x * x` differ is `x * x` in both paths, byte for byte; the fast path's `where` guard's test moved to a
`-1.0` exponent (numpy reciprocates a 0-d array where a scalar goes through `pow`; a square no longer reaches it); four
adopter records' `explain` digests re-recorded (the explanation shows the chain source, now `_sq(..)`; their value
digests did not move).

**Gate**: suite 5 479 passed + the 5 above, re-made (137 of 137 in the three changed files); the C-against-numpy suites
green, no pinned pair moved.  The gold census on the three paths, identical to each other: 25 decks move (the last bits),
no column crosses, 0 forced, the largest metric change ring's 0.011.  The 27 decks against HEAD (quiet, alternating
trees): 54.92 -> 42.02 s (0.77x) -- jge 0.68x, gm19 0.81x, mux8 0.79x, voter25 0.80x, rich3 0.76x, arom 0.77x; the four
big decks 42.27 -> 31.51 s (-25.5 %).  One deck the other way: opampal 1.91 -> 2.06 s -- its transient faster (0.98 ->
0.92 s), its operating point slower (1.00 -> 1.21 s): the same route (source steps, the junction-gmin ladder) taking
another path through its failing rungs (PCNR iterations 1 286 -> 1 632) -- a hard operating point's sensitivity to its
last bits, as a latch's state is.

## 2026-10-17 — 147. Level 3 without `vmax`: a finite Jacobian (the reverse-mode ratio clamped below)

**Found** diagnosing the gold decks' MOS failures (plan `dc-final-gmin.md`, 'THE OPEN TOPICS'): an inverter chain on
gm1's level-3 cards with `vmax=0` -- SPICE's default, every card that does not give it -- had no operating point (PCNR's
iterate not finite at iteration 2, the chain exhausted).  The device's current was finite and its Jacobian NaN in
saturation and near threshold.  Cause: `ratio = minc(vds/vdsatf, 1)` was clamped above only.  The reverse-mode copy of
the channel runs with `vds < 0` at every bias, where `vdsat` floors at 1e-9: the ratio reached -2.5e9, its fourth power
in `dl_below` 4e37, `delxl` rounded to `leff` and `xlfact = 1/(1 - delxl/leff)` divided by zero.  With `vmax > 0` that
branch is the constant 0, so no gold card (all give `vmax`) met it, and the finite-difference test's no-`vmax` points
on the library's card did not either.

**Fixed.**  The ratio is clamped below at 0 as well: where its arm is selected it is in [0, 1] already, the clamp the
identity there -- every selected value and derivative bit-identical.

**Checked.**  gm1's card without `vmax`: the current and Jacobian finite on a 1 050-point grid (forward and reverse, cut
off to saturated); the inverter chain's operating point solves; its DC current agrees with ngspice's on the same deck
text (48 points, within 0.5 %, one 1 nA subthreshold point 3 %).  Test: `test_level3_without_vmax_is_finite_at_every_
bias` (the library's card and gm1's; gm1's fails on the code before).

**Gate**: suite 5 486 passed, 6 skipped, 3 xfailed (5 484 + the two cases); every gold operating point identical; the gold
census on the three paths: no deck moves.

## 2026-10-17 — 148. An opt-in SPICE Meyer mode for MOS levels 2 and 3 (`import_netlist(mos_charge='meyer')`)

**Why** (plan `dc-final-gmin.md`, 'THE OPEN TOPICS', 'AN OPT-IN SPICE MEYER MODE'; Andreas: "Opt-in SPICE Meyer mode").
The gold decks' failing columns, classified by the three-way check (ngspice on the same imported circuit against the same
gold): 14 of 30 fail in ngspice too; 16, in 8 decks -- all MOS level 2 / 3 -- fail in ours alone.  Measured, not
assumed: not integration (tolerances x0.1 move nothing), not DC (the drain current is ngspice's on a 48-point grid, both
devices, the same deck text read by both), but our edges LATE -- 20..30 ps a gate stage on gm1's cards, g1310's chain
accumulating +60 a stage -- the gate's total capacitance Meyer's, its split not.  On a 5-stage chain switched off one at
a time on both sides: eta's DIBL left out of the charge's threshold ~30 %, the junctions a little, the body effect and
velocity saturation nothing; the rest is a charge-conserving model against Meyer's capacitances, which no charge has --
and which Xyce's gold and ngspice integrate without conserving anything (the choice declined at design time, 2026-10-07).

**Built.**  `_meyer_caps`: Meyer's capacitances as ngspice's `DEVqmeyer` computes them (its 25 mV `MAGIC_VDS` floor,
its halves doubled), each mode's own `von` (eta's DIBL in, as SPICE's) and `vdsat`, SPICE's mode swap.  `_meyer_currents`:
each as a capacitance-defined current `i = C(v) dv/dt` through a scratch node (`V(y) - V(minus) = tau dv/dt`, the
memristor's idiom; the branch current `(C/tau) V(y)`; tau = 1 ps), compiled, on the core and the C loop, nothing at DC.
`MosLevel{2,3}[Pmos]MeyerHdl`; `import_netlist(mos_charge='meyer')` ('conserving', today's, the default; level 1 refused);
`spice_suite --meyer`.  Found building: (a) the compile took over 30 minutes -- vector PCNR's route analysis resolved the
Meyer coefficient through the whole chain (sympy `subs` on nested `Piecewise`) only to refuse it: `Behavioural.
pcnr_declined`, a class's own reason, skips that analysis, and the Meyer classes set it (their operating points the
ordinary chain's, as SPICE's); (b) with tau = 1 ns a 5-stage chain collapsed its step at 1e-14 s -- the scratch value
carries v's Newton jitter times tau/h.

**Measured.**  M1: the capacitances against ngspice's own (`@m[cgs/cgd/cgb]`, gm1's card, 24 biases): 2.7e-6 in inversion
and accumulation (the oxide constant), 6e-4 .. 4e-3 near threshold -- the formula fed ngspice's own `von` and `vdsat`
leaves the same residual: the DC model's `phi` / `von` constants there, not Meyer's.  M2: the 5-stage chain against
ngspice, per node, -2 .. +1 ps (the conserving model +134 / +147 ps by the fifth).  M3 (the gold decks under the mode):
the 16 ours-only columns pass but slowlatch's (not run), the ring passes both (48.6 / 39.8 by default), no column passing
by default fails; failing columns 30 -> 7, each one ngspice fails too and ours at or below its (rich3 2.1-2.7 against
~90).  M4: the default path bit-identical -- the gold census on the three paths moves no deck, every operating point
identical.  Cost: gm19 12.1 s against 5.0, mux8 11.6 against 7.7.  NOT RUN under the mode: jge (its operating point,
PCNR declined, the chain alone fails it), slowlatch and todd3 (the step collapses mid-run: their 1.5e-12 A tolerance
against the scratch coupling C/tau ~ 10 S).  Tests: `test_mos_meyer.py` -- DEVqmeyer in both modes on levels 2 and 3 to
1e-9 Cox, the PMOS mirror, the DC current the charge classes' byte for byte (gmin on), PCNR declined by the reason, the
importer, an imported chain on the C core switching sooner than the conserving model's; `test_hdl_cache`'s reclass keeps
`pcnr_declined` (dropped, its copies ran the route analysis for hours); `test_device_gmin` and the PCNR coverage guard
take a declared decline.

**Gate**: suite 5 542 passed, 6 skipped, 3 xfailed; the gold census (default) on the three paths: no deck moves.

## 2026-10-17 — 149. The Meyer mode runs every gold deck: its scratch relation potential-contributed

**Why** (Andreas: "Make the Meyer mode run everywhere").  Under roadmap 148's mode three decks did not run: jge had no
operating point and slowlatch and todd3 collapsed their step.  One cause: each scratch node's relation `V(y) - V(minus)
= tau dv/dt` was a node equation, and the conductances a solver adds to every node -- the gshunt ladder, pseudo-transient
continuation, the step rescues -- moved it; through `C/tau ~ 10 S` that was a fake gate current (jge's operating point:
a residual of 1.2e4 A at node 20, the continuation failing at g = 1e6).

**Built.**  The relation is POTENTIAL-contributed (`Contribution(by.V, tau*ddt(br.V))`): a branch-current unknown
enforces it exactly, and no added conductance moves it -- at DC the Meyer currents are exactly 0 whatever a ladder
does.  One unknown more a capacitance.  The prototype (a constant capacitance this way against a linear capacitor) is
unchanged: 1.9e-13 V, the same 274 steps.

**Measured.**  All 27 gold decks run under the mode (jge 207 s).  Its failing columns: 12, every one but slowlatch's a
column ngspice fails too and ours at or below its (todd3 3.16 / 2.73 against 3.18 / 2.73 -- the same latch state);
slowlatch's three (v(3) 136, v(4) 32, v(5) 43) are its OTHER latch state -- no .ic or .nodeset, so the operating point's
solver picks it, and under the mode that is the ordinary chain (the Meyer classes decline PCNR), which picks the other
one, as it did by default before roadmap 142.  The chain against ngspice: -2 .. 0 ps a node.

**Gate**: suite 5 542 passed, 6 skipped, 3 xfailed; the gold census (default) on the three paths: no deck moves.

## 2026-10-17 — 150. Meyer devices take part in PCNR: their operating points are the default mode's; a redundant probe's row spans the whole block

**Why** (Andreas: "Do the follow-up now" -- roadmap 149's: let the Meyer devices take part in PCNR, since their currents
vanish at DC).  The Meyer classes declined vector PCNR (`Behavioural.pcnr_declined`, roadmap 148): the route analysis ran
over 30 minutes resolving the Meyer coefficient through the chain (`_reads_a_probe`) only to refuse, so the ordinary chain
solved the mode's operating points.

**Built.**  `_meyer_currents` reads each scratch voltage through an identity probe (`limit_identity(by.V)`): the Meyer
current is then a function of the device's probes, the route analysis takes 2.6 s and qualifies the device -- the charge
class's three unknowns and three more, the scratch voltages, which the branch row pins at 0 at DC.  `pcnr_declined`
removed (no class declines any more) with its test exemptions; `test_device_gmin` sweeps the four Meyer classes too
(16 -> 20).  **A latent PCNR defect, found by it**: a redundant probe's coefficients were taken over the unknowns that
existed when it closed its cycle, so a probe declared AFTER it (the Meyer MOSFETs' identity probes, after the group's
redundant `(b,d)`) left the row short -- `limit_block` raised on the ragged matrix and every solve fell back.  The row is
padded with zeros to the block's width (`_pcnr_declared_route`); no library device had declared a probe after its
redundant one.  `test_pcnr_vector::test_a_probe_declared_after_the_redundant_one_widens_its_row` fails without the fix.

**Measured.**  The operating points under the mode against the default's: gm1, rich3, todd3 within 1e-14 V, PCNR used
in both.  jge: the two modes' PCNR iterates on the shared nodes are bit-identical for two iterations; the third differs
by 1.4e-7 -- the default's own held solve is heading into a near-singular step (the next iterates reach 7e5 V and 5e10 V)
that amplifies rounding ~1e9x -- and the two take different ladders (the default's held solve fails, the Meyer's does
not) to the same point within 3e-6 V.  Seeding the scratch nodes consistently changes nothing (the inconsistent start
was refuted).  **slowlatch's operating point is a coin flip in BOTH modes**: from starts nudged by 1e-12 V, DC lands in
each latch state about half the time (default 8 / 8, Meyer 7 / 9 of 16) -- a symmetric latch with no .ic or .nodeset.
The default matched ngspice by luck; roadmap 149's "the chain picks the other one" attributed a lottery to a solver.
Cost (the operating point, PCNR against the chain under the mode, interleaved, minimum of 3, on a loaded machine): gm19
3.6 s against 14.3, gm17 0.05 against 0.58, gm3 0.01 against 0.20; e1480 2.75 against 0.06 -- five failing rungs of the
same gshunt ladder the default climbs, which the C core hands back on floating-point flags (`pcnr_c:bail:flags`) and the
Python loop reruns, 0.3 .. 0.8 s each, where the default's fail in C as non-finite; jge 13.3 against 3.9 (2 332 unknowns
against the default's 244, 489 C iterations).

**Gate**: suite 5 557 passed, 6 skipped, 3 xfailed; the gold census (default) on the three paths: no deck moves; all 27
decks under the mode, 12 failing columns as in 149 (gm17 v(5) 10.0 -> 10.9, ngspice 12.9; slowlatch's three its other
latch state); the chain against ngspice -2 .. 0 ps a node.

## 2026-10-19 — 151. A NaN goes through the C PCNR core's KLU sequence: the failing rungs it handed back are served

**Why** (roadmap 150's open cost; Andreas: "Fix the C core's flag bail").  Under the Meyer mode e1480's operating point
took 2.64 s against the chain's 0.06: five rungs of its gshunt ladder the C core handed back (`pcnr_c:bail:flags`) and
the Python's loop reran, 0.3 .. 0.8 s each.  The default circuit (205 unknowns) solves dense and the same rungs fail in
C; the Meyer circuit (373) takes StructuralKLU, and `solve_klu` handed back on any entry or rhs that was not finite --
"the sequence's guards are twinned on finite values only".

**Measured first.**  The five reruns emit no numpy warning; in each the first non-finite KLU input is the rung's LAST
solve, NaN only (8 .. 54 entries, 2 .. 14 rhs values, no inf), and it ends as the C would have said: two non-finite
iterates (`NONFINITE`), three failed factors (`Singular matrix`).  `pcx_kluseq.h`'s residual guard already propagates
a NaN as numpy's `max` does (`nan_r`): every decision the sequence makes on a NaN is the Python's.

**Built.**  `solve_klu` lets a NaN through; an infinite entry or rhs still hands back (unmeasured: the residual's `arr @
x - rhs` may be numpy's warning).  `test_pcnr_c`: a NaN on the rhs, in the matrix (`_Root`, a square-root conductance at
a negative probe -- a `var` so that the class is chained, which the C core requires) and on a held numeric (the refactor,
then the residual's fresh factor) served as the Python's error with StructuralKLU's counters the same, each failing on
the old C; an infinite rhs still handed back.  The hand-back test it replaces pinned the old behaviour.  The pinned pair
"the PCNR DC solve and its C" re-recorded (its reference unchanged).  `dcnewton.c` keeps its own, broader finite checks
(StandardNewton's semantics; not measured here).

**Measured.**  e1480 under the mode: all 22 rungs served in C (2 non-finite, 3 singular), each rung's outcome and the
operating point bit-identical C on and off; the operating point 2.64 -> 0.26 s (back to back; the chain 0.06).  jge
11.5 s against the chain's 3.9 (its 2 332 unknowns), gm19 2.8 against 6.1.

**Gate**: suite 5 559 passed, 6 skipped, 3 xfailed (`test_c_abi`'s no-history rule then moved the measurement out of the
C comment, here); the gold census (default) on the three paths: no deck moves; the 27 decks under the mode: no metric and
no step count moves; the chain -2 .. 0 ps a node.

## 2026-10-19 — 152. The C PCNR core on the Jacobian's pattern: no dense pass, every verdict the Python's -- jge's Meyer operating point 11.7 -> 1.1 s, opampal's 1.37 -> 0.36 s

**Why** (Andreas: "jge's Meyer operating point", then "Build the sparse PCNR core").  Sized first (plan 'JGE'S MEYER
OPERATING POINT, SIZED'): the Meyer mode's jge took as many PCNR iterations as the default's (489 against 579) at 23.9 ms
each against 1.0 -- perf: 67 % the C's own dense n x n passes (2 332 unknowns), 12 % two dense dgemv's, KLU 0.7 %.

**Design: the Python reference unchanged, bit for bit.**  Everything dense but two products is data movement and adds
at known positions; the two products -- `|J_eff| |x_new|` for the convergence test, StructuralKLU's refactor residual
`arr @ x - rhs` -- feed only VERDICTS.  (a) J stays a dense buffer, zero off a pattern P (the view's core pattern -- binG
is zero off it -- each Schur position, the junction pairs' entries, the diagonal): allocated zero and written only on P,
so `J[P] = binG[P]` is the dense copy, KLU's entries are gathered where they stand and `_ours`' foreign count is a test
of P's entries off KLU's pattern.  (b) A sum rounds only where two operands holding a nonzero term meet (a zero product
joins exactly, fused or not), so a row with k nonzeros summed in order is within 2 gamma_(k+1) T of numpy's dgemv in any
summation tree (T = sum |a| |x|); with d = 4 (k + 2) eps both tests -- monotone in the product, rounding monotone -- are
decided where they hold or fail at both ends of the band, and the dense product decides where they do not (also any
non-finite value, any row reaching 1e300).  (c) `pcx_kluseq.h`'s sequence takes the sparse verdict behind a flag
(`a_built`); `dcnewton.c` keeps the dense test.  The premise -- numpy's BLAS rounds a row once a nonzero term -- is a
test against the installed BLAS: `test_numpys_dgemv_rounds_a_row_once_a_nonzero_term` (sizes 2 .. 3001, cancelling and
nonnegative rows, one thread and the pool; alive: 29 % of rows differ from the in-order sum, the worst 2.8 eps T, 0.07
of the bound).

**Measured.**  jge under the Meyer mode: every rung bit-identical C on and off, the operating point 11.7 -> 1.08 s
(line: 3.0 s); e1480 under the mode 0.26 -> 0.15 s.  The convergence test never needed the dense product; the residual
test did where the residual sits at its rounding floor (rhs tiny against |A| |x|: T / scale 1e6 .. 1e11; jge 31 times,
e1480 57), where the reference's own verdict is set by its dgemv's rounding and no bound can decide.  **The default
decks gain too** (my sizing said they would not): their operating points, old C against new interleaved twice, the best
of 12, a loaded box with free cores -- opampal 1 373 -> 361 ms, gm19 344 -> 136, rich3 26.9 -> 12.4, jge 689 -> 499,
mux8 30.7 -> 25.4, voter25 315 -> 265, e1480 224 -> 199, toronto 57.3 -> 53.2, arom 13.5 -> 13.2, gm1 5.34 -> 5.19;
no deck slower (line: +2 %).

**Gate**: suite 5 562 passed, 6 skipped, 3 xfailed; the gold census (default) on the three paths: no deck moves; the 27
decks under the Meyer mode: no metric, no step count moves; the chain -2 .. 0 ps a node.  Pinned pairs "the PCNR DC
solve and its C" and "the DC base solve and its C" re-recorded (their references unchanged).

## 2026-10-19 — 153. KLU below 250 for DC, Radau and the stage methods: measured, refused at 40; the C PCNR core's solver check fixed

**Why** (Andreas: "Start Radau's and DC's sparse thresholds").  StructuralKLU starts at 40 unknowns for the multistep
methods (measured 2026-10-16) and at 250 for everything else, `MIN_N_FOR_SPARSE`, measured on the Python loop with
SuperLU.  The question: do DC, Radau (RadauIIA3), TRBDF2 and ESDIRK43 gain from KLU in the 16 gold decks with 40 <= n <
250?  NB the transient's initial operating point is DC's own solve: DC's floor moves every band deck's operating point.

**Lines** (plan, before measuring): floor 40 for an analysis/method only if the band's summed time falls >= 5 %, no
deck is slower than +3 % (best of 3 interleaved, fresh processes) and no gold column flips pass/fail.

**Measured** (scratch g12/thr/).  All four REFUSED by the per-deck line; no gold column flips anywhere.  DC: the band's
operating points 1019 -> 379 ms with the C serving, but mike2 +31, gm3 +15, ring +14, gm1 +11.5, latch +5.5, slowlatch
+4.9 %; every deck from 149 faster (-33 .. -83 %); gm17's operating point flips latch state (a rounding lottery).
Radau: band +2.0 % (ring +49, mike2 +21 ...; from 130 gm1 -28 %, voter25 and jge flat).  TRBDF2: band -12.4 % (mike2
+13, todd3 +10, gm3 +9, slowlatch +7 %; from 130 all -5 .. -19 %).  ESDIRK43: band -13.8 % (todd3 +11, mike2 +9 %; from
130 all -6 .. -20 %).  Below ~100 unknowns KLU's fixed cost loses; from ~150 it wins for DC and the stage methods.
**Radau stays flat on the big decks** because its C transformed stage solve (`radau_tc`) twins only the DENSE real
factor -- numpy's LU, with the complex factor on complex KLU always -- so with a KLU real factor every iteration runs
in Python: `radau_tc:served` 0 on every deck at 40.  Where the C served most iterations (jge: 17 139 of 20 332 frozen
transforms) losing it cancels the cheaper solve (110.97 -> 111.82 s); where it rarely served -- it hands back on
floating-point flags, gm1 833 bails against 175 served -- the cheaper solve shows (gm1 -28 %).

**Found, fixed: the C PCNR core's solver check** (`_pcnr_c._solver`) took the AutoSolver's first choice below `min_n`
as dense and made it on a zero matrix.  With a structural floor below `min_n` that chose an unanalysed StructuralKLU,
the C declined ('linsolver': the band's operating points +861 % in the first measurement) and the Python then solved on
a choice its own first matrix might not have made.  The check now uses `_select`'s own floor (the lower one where a
structure exists) and lets the first Schur matrix choose above it.  Unreachable today (DC's two floors are both 250);
`test_a_structural_floor_below_the_sparse_one_is_served` fails on the old code.  `dcanalysis._blas_held`'s same
assumption is kept, commented: holding BLAS to one thread is the safe guess there (harmless to KLU, needed by dense).
**Found, open: Radau crashes** on mux8 at floor 40 -- its dense coupled stage Newton (`_stage_newton`) calls
`np.linalg.solve` bare, and a singular matrix's LinAlgError escapes the transient past its own rescues.

**Gate**: suite 5 562 passed (the no-history rule then moved a date out of `_pcnr_c`'s comment); the gold census on the
three paths: no deck moves.  The pair "the PCNR DC solve and its C" re-recorded (its reference unchanged).

## 2026-10-19 — 154. A singular coupled Radau system fails the attempt, not the run

**Why** (roadmap 153's open finding; Andreas: "handle the two findings").  Radau's mux8 with KLU's floor at 40 crashed:
on that trajectory the transform fell back to the dense coupled stage Newton (`_tran_radau._stage_newton`), its
3m x 3m solve met a singular matrix, and numpy's LinAlgError left the transient -- past the step's own rescues, which
catch `NoConvergenceError` only: the rejection and smaller retry, and at `minstep` the gshunt ladder, whose
conductance makes exactly this matrix regular.  The DC solvers wrap a singular factorisation as `NoConvergenceError`;
the coupled PCNR path already falls back on any exception (Andreas's rule there: "a singular Jacobian included -- must
not end the run"); this was the one bare solve on a Radau step path.  Its C twin (`_tran_radau_c`) hands a singular
system back to the Python, which therefore always answered it.

**Built.**  `_stage_newton`'s solve raises `NoConvergenceError('Radau IIA(3) coupled stage Newton met a singular
Jacobian' [at gshunt=...])` from the LinAlgError.  `test_radau_c`: a zero Jacobian planted once in the Python Newton's
fifth assembly -- the run finishes, one rejection more, its end within 1e-6 V of the run without it (fails on the old
code with the LinAlgError); the C's singular hand-back now answers that error.  Pinned pairs "the coupled stage
Newton and its C" and "the stage paths' passes and the core's" re-recorded (their references moved, the twins not).

**Measured.**  mux8 under Radau with the floor at 40 now runs: 13 726 steps (13 722 at 250), 4 268 rejections (4 234),
every gold column the same.

**Gate**: suite 5 564 passed, 6 skipped, 3 xfailed; the gold census on the three paths: no deck moves.

## 2026-10-20 — 155. KLU on the structure from 150 unknowns for DC (every operating point) and the stage methods

**Why** (roadmap 153's finding; Andreas: "handle the two findings").  Floor 40 was refused for DC, TRBDF2 and
ESDIRK43 by its per-deck line, but every deck from ~150 up had gained.  A floor of 150 was a NEW line, written before a
fresh measurement (plan 'A FLOOR OF 150 FOR DC AND THE STAGE METHODS'): over the five decks it moves (toronto 153,
mux8 171, voter25 200, e1480 205, jge 244; gm17 at 149, a rounding-lottery operating point, stays dense), each of
DC's operating point, TRBDF2's and ESDIRK43's transient: summed time -5 % or better, no deck slower than +3 %, no gold
column flipping, DC's operating points within 1e-8 V.  Radau keeps 250: its C transformed solve twins only a dense
real factor, and losing it cancelled the gain (roadmap 153).

**Measured** (fresh processes, three rounds interleaved, best of 3; both floors lowered together, as built): DC's
operating points 927 -> 328 ms (-64.6 %: toronto -40, mux8 -35, voter25 -57, e1480 -28, jge -83 %), within 2e-9 V;
TRBDF2 308.7 -> 262.2 s (-15.1 %, each deck -3.0 .. -18 %); ESDIRK43 106.0 -> 86.0 s (-18.9 %, -11 .. -23 %); no gold
column flipped.  The line met by all three.

**Built.**  `linearsolver.MIN_N_FOR_STRUCTURAL_KLU_OP_AND_STAGES = 150`; `DC._STRUCTURAL_KLU_MIN_N` (so every
analysis's operating point, a transient's included); `Transient._STRUCTURAL_KLU_MIN_N` per method -- 40 on the C
loop's multistep methods, 150 on TRBDF2 and ESDIRK43, `min_n` on Radau and the GLMs.  The C PCNR core serves the band
since roadmap 153's solver-check fix.  `test_structural_klu`: an operating point at 140 sections the dense call byte
for byte and at 160 on KLU within rounding; each method's floor from both sides.

**Gate**: suite 5 571 passed, 6 skipped, 3 xfailed.  The gold census on the three paths moves EXACTLY the five decks
-- their operating points now on KLU, their transients' last bits and step counts with them (jge 10 746 -> 10 734,
mux8 17 000 -> 16 921, toronto 766 -> 774, voter25 9 959 -> 9 988, e1480 unchanged) -- metrics within 0.006, no
pass/fail crossing: the re-baseline this line named.

## 2026-10-20 — 156. Radau's transform C: an underflow no longer hands the loop back, and maxiter is served

**Why** (Andreas: "Address the two points related to Radau now").  By default `radau_tc` handed the transform loop
back on floating-point flags on most calls (gm1 833 of 1 205, voter25 5 983 of 9 670), and the Python loop ran them.
**Measured first** (the C's flag tests instrumented, eight decks): about 23 000 hand-backs, every one FE_UNDERFLOW (the
residual's products, a few in the right-hand sides) but 6 FE_INVALID on voter25, every call under numpy's default
error state, which ignores underflow.  The C tested underflow too (pcx_fenv.h: "a needless hand-back at worst");
`dcnewton.c` and `pcnrdc.c` test only what numpy reports and their glue declines where underflow is not ignored.
A flag never changes a value.

**Built.**  `radau_tc.c` tests overflow, invalid and divide (`TC_FLAGS`); `_tran_radau_tc.solve` declines where numpy's
error state does not ignore underflow (`radau_tc:underflow`, `func._underflow_ignored`).  Found before measuring: with
underflow no longer handing back, a loop that does not converge ran all 100 iterations in C, handed back, and the
Python repeated them before the dense fallback (latch 1.88 -> 2.08 s).  The C's maxiter is the Python loop's own False
(its iterates, bit for bit): served as such (`RADAU_TC_MAXITER`, `radau_tc:noconv`; the source memo's hits counted as
the loop makes them), and `_rk_step_transformed` raises as before.  `test_radau_tc`: an underflow in the residual (the
tableau scaled by 1e-300) served bit for bit, and under `errstate(under='warn')` declined with the Python's warning;
maxiter the Python's False at 1 and 3 iterations; the drawn-seed liveness test now asks for both served outcomes --
each failing on the old code.  Pinned pairs "the transform Newton and its C", "the step's end passes and the C's" and
"the stage paths' passes and the core's" re-recorded.

**Measured** (Radau's default transient on the 16 decks from 40 to 249 unknowns and the four from 250; today's tree
against this one, fresh processes, three rounds interleaved, best of 3; the line written first): every waveform
bit-identical; 454.6 -> 421.6 s (-7.3 %): latch -61, slowlatch -47, todd3 -37, gm17 -19, gm1 -18, mux8 -18, nand -17,
toronto -14, voter25 -12, e1480 -10, mike2 -6, jge -3.8 %; the decks the C never serves within -0.5 .. +0.8 %.

**Refused on the way** (plan 'RADAU'S TWO C GAPS').  (1) The same rule in `radau.c`, the coupled stage Newton's C:
built with the above and measured with it -- arom +5.4 % (each round, both orders): its coupled C then served 39 solves
it had handed back, and at 3m = 882 a served solve costs more than numpy's (unmeasured why).  Dropped; the line
re-measured without it (above).  (2) The transform's real factor on StructuralKLU in C (`pcx_kluseq.h`'s sequence, the
C serving 99 % of the steps of rich3, arom, gm19 and opampal, every waveform bit-identical): -3.5 % summed against a
-5 % line (rich3 -7.7, arom -6.7, gm19 -0.9, opampal -0.3 %).  Not committed.  Why flat: on opampal (single-threaded
profile, 15.7 s) the whole transform loop is 1.7 s; Radau's error estimate is 6.6 s, 3.6 s of it a dense
`numpy.linalg.solve` of `(gamma_r/h) C + G` every step -- the transform's real factor in principle, built again from
`GAMMA_REAL` (2 ulp from the transform's `eig` value) and factored afresh, dense even on a KLU run.

**Gate**: suite 5 573 passed, 6 skipped, 3 xfailed.  The gold census does not run Radau (no deck sets the integrator): not run.

## 2026-10-20 — 157. Radau's error estimate solves on the transform's real factor

**Why** (roadmap 156's profile; Andreas: "Do the next lever").  On opampal (518 unknowns; cProfile, single-threaded)
Radau's error estimate took 6.6 of 15.7 s, 3.6 s of it a fresh dense `numpy.linalg.solve` of `(GAMMA_REAL / h) C(x_n)
+ G(x_n)` every step -- the toolkit's solver, dense even on a KLU run.  On the transform path that is the step's real
factor in principle: the same C and G at x_n, already factored once a step (`_radau_frozen`); radau5 uses the one
factor for both.  Its scalar is the transform's `lam_0` (`eig`'s), 5 ulp above GAMMA_REAL (3.637834252744498 against
...496, the exact 3 + 3^(2/3) - 3^(1/3)).

**Built.**  `_rk_step_transformed` hands `_finish_radau` its real solve with the C and G at x_n it is of:
`_FrozenTransform.real_solve` (the frozen solve's own real solve, extracted: the kept LU for a finite right-hand side,
the analysis solver otherwise), or, where no frozen factors were made, the per-iteration solve's own
(`_radau_transform_solve`'s factor and solver: the same bytes, so the frozen transform stays the per-iteration one's).
`_radau_error_estimate` solves on it where its own C and G at x_n are those (compared: a stateful limiter's sync may
move them) on the numeric toolkit (`radau.est:transform`); elsewhere -- the coupled paths, a mismatch -- as before
(`radau.est:fresh`).  A RE-BASELINE of Radau's adaptive answers: the estimate moves by the 5 ulp (and by KLU's
rounding from 250 unknowns), the next step size with it.  `test_radau_frozen`: the estimate on the kept LU is numpy's
solve of `(lam_0 / h) Cr + Gr` bit for bit, on StructuralKLU its solve; within rounding of the old one; another C ->
the fresh solve, the old bytes; a KLU run the same with the frozen factors off.  Pinned pairs "the stage paths'
passes and the core's", "the transform solve and its frozen form" and "the transform Newton and its C" re-recorded.

**Measured** (line written first: plan 'THE ERROR ESTIMATE ON THE STEP'S REAL FACTOR'; Radau's default transient on
the 20 decks, fresh processes, three rounds interleaved, best of 3 -- the disk filled during round 3, the 8 largest
decks have 2): summed 425.2 -> 372.9 s (-12.3 %): gm19 -34, opampal -34, rich3 -13, gm3 -13, arom -12, jge -9.3, mux8
-7.6 %; no gold column flips.  THE LINE FAILED on two criteria: gm1 +19.1 % and slowlatch +13.9 % (3 rounds each), and
accepted steps beyond +-2 % on gm3 -10.9, g1310 -7.0, gm17 -3.9, arom -3.4 %.  A lottery check on today's code --
GAMMA_REAL alone moved by -2 / +2 ulp -- moves those step counts as much (gm3 -7.0 / -12.2, g1310 -3.6 / +9.5, gm17
-0.8 / +3.1 %): the controller's sensitivity to the last bits of its estimate, which no reuse of the factor can avoid
(its matrix must be the transform's).  Andreas's decision, offered refusal, a new lottery-aware line, or acceptance:
accepted despite the line.

**Gate**: suite 5 576 passed, 6 skipped, 3 xfailed.  The gold census does not run Radau: not run.

## 2026-10-21 — 158. glibc's allocator mode pinned on request; Radau's frozen factors made in place, refused

**Why** (Andreas: "Do 1 and 2", the open topics after 157 -- the cost of `fold`, Radau's frozen factors, and the
coupled stage Newton's C `radau.c` at large n; then, on what item 2 found, "Do as you suggest").

**Item 2, sized: the allocator's, not the C's.**  `radau.c` under 156's underflow rule made arom +5.4 % and was dropped.
It is not slower than the Python: arom on that tree, C on against off -- unpinned 28.37 / 26.60 s; a cache-friendly
tiled assembly 28.72 s (a strided-write hypothesis refuted); MALLOC_MMAP_THRESHOLD_ 128 KB 40.84 / 41.88 s; 32 MB
30.61 / 30.60 s; 32 MB with MALLOC_TRIM_THRESHOLD_ 256 MB 26.07 / 26.22 s; bits identical everywhere.  The
single-threaded profile puts its 39 solves at +0.3 s net and +1.6 s on functions every step runs.  glibc serves an
allocation above `M_MMAP_THRESHOLD` (128 KB at start) with fresh mmap'ed pages and raises the threshold (and the trim
threshold to twice it) when such a chunk is freed: a process settles into a MODE by its allocation history, which the
C's changes.  Nothing built for `radau.c`.

**The allocator's mode: a decision rule written first** (plan 'GLIBC'S ALLOCATOR MODE'): pin both thresholds BY
DEFAULT if, pinned, the 27 gold decks' summed time falls >= 3 %, no deck is slower than +3 % (the 27 and Radau's 20),
every waveform is bit-identical and no deck's peak RSS grows more than 10 %; else on request, which the timing
harnesses make.  **Measured** (the 27 decks' default transient, unpinned against MALLOC_MMAP_THRESHOLD_ 32 MB +
MALLOC_TRIM_THRESHOLD_ 256 MB, fresh processes, three rounds interleaved, best of 3): every waveform bit-identical;
summed 40.61 -> 40.62 s (+0.0 %); the decks under 0.2 s -15.8 .. +9.0 % (noise at 30 ms), the rest -0.9 .. +0.8 %;
peak RSS at most +4.3 % (mux8, 501 -> 522 MB), but a pinned process keeps more at its end (mux8 459 -> 508, jge 354 ->
402, voter25 336 -> 366 MB).  Radau on its 20 decks (the fold A/B's today's-tree rows): 366.24 -> 362.69 s (-1.0 %),
opampal -5.1 %, worst +1.1 %, bits identical.  The rule says ON REQUEST.

**Built.**  `pycircuit.circuit._malloc`, imported first by `pycircuit.circuit`: with env `PYCIRCUIT_MALLOC_PIN=1`, or
on `_malloc.pin()`, `mallopt` sets `M_MMAP_THRESHOLD` 32 MB (glibc's ceiling for it on 64 bits) and `M_TRIM_THRESHOLD`
256 MB, which also switches the adjustment off; glibc only; a threshold the user set (the `MALLOC_*_` variables, a
`glibc.malloc.` tunable) wins; `STATUS` says which.  Pinned at import is the variables' pinning at process start
(opampal 9.41 / 9.49 s, arom 23.24 / 23.15; unpinned 10.00, 23.56; three rounds, bits identical).
`benchmarks/_bench.pin_malloc()` beside `pin_threads()` in `micro`, `insitu` and `loop_split`, and the gold census
(`spice_suite`) sets the variable for its decks' processes unless it is set: an A/B there measures a change, not the
mode its allocations flip.  `test_malloc_pin`: by glibc's own count of mmap'ed chunks (`mallinfo2`), a 24 MB array
comes from the heap once pinned (on import and by `pin()`) and from mmap under a threshold the user set.

**Item 1, built and refused.**  With the error estimate on the transform's factor, `fold` -- `_radau_frozen`'s two
factors and the complex one's packed values in one C call -- was the largest single piece of opampal's step (518
unknowns, single-threaded profile: 3.3 ms a step, 28 % of the run).  Sized on opampal's own factors (warm, a call):
2.28 ms -- the inputs copied into the C's buffers 0.18, the C 1.53 (the gather of the packed values column by column
through the m x m factor 0.66), the factors copied out 0.27.  Built, byte for byte as before: the C read `Cr` and `Gr` in
place and wrote the two factors into the arrays `fold` returns; the gather took the last record's own CSC and proved
the pattern exact by count (every position nonzero, rows strictly ascending, the spans in order, as many nonzeros in
the factor as the pattern holds); 2.28 -> 1.04 ms a call.  Measured (line written first; Radau's default transient on
the 20 decks, three rounds interleaved, best of 3, unpinned and pinned): -0.9 / -0.7 % after 1.5 rounds -- `fold` ran
only where the step's `P` was cached, so an adaptive run made its factors in numpy almost every step (latch: 31 folds
in 402 steps).  A new stage with its own line (>= 3 % in each mode): `P` first under its own `errstate`, `fold` on every
step whose scalars fold.  Every waveform bit-identical; summed 366.24 -> 359.73 s unpinned, 362.69 -> 356.00 pinned
(-1.8 % each); opampal -22.8 / -23.2 %, the rest -3.4 .. +0.4 %.  Below its line: REFUSED, not committed (the patch
kept outside the repository).  The gain is a fixed-step dense deck's alone.

**Gate**: suite 5 579 passed, 6 skipped, 3 xfailed.  Nothing changes by default; the 27 gold decks pinned are bit-identical
to unpinned (above).

## 2026-10-21 — 159. Radau's early abort: an attempt over its tolerance at the transform's first iterate is abandoned there

**Why** (Andreas: "Start with 1 and 4", then "build prototype", then "accept").  158 left item 4 -- `radau.c`'s
underflow rule re-measured with the allocator pinned -- and item 1, an early abort for Radau, which the multistep
family's (132, 134) does not reach.  **Item 4: refused** (line: the summed time falls, no deck +3 %, bits identical;
measured +0.3 %, bits identical, every deck -2.7 .. +1.5 %): the C now served the solves it handed back (voter25 323
against 47), worth nothing measurable -- the coupled path is rare -- and it settles that 156's +5.4 % on arom was the
allocator (pinned, arom -1.2 %).  Nothing committed for it.

**Item 1, sized first** (a time census in the default configuration and a separation census with the transform loop
in Python, the same bytes, attempt by attempt; 394 s over Radau's 20 decks).  Rejected attempts are 45.9 % of Radau's
time.  736 attempts ran the transform's simplified Newton to `maxiter` (100; the loop has no rate test) and then the
dense coupled Newton: 31.0 % of the time, 614 of them REJECTED after (voter25 73 % of its time, todd3 63 %).  The
embedded estimate at the transform's first iterate, normalised as the step is judged, is the final one to a few %
(median ratio 1.00, p10 0.99, p90 1.04).  Over 2 there: 27.8 % of the time saved gross, 81 accepted attempts wrongly
abandoned.  radau5's rate exit to the fallback was sized too and is NOT built: it sends slow convergers to the dense
fallback (gm19 -36 % at NIT 20).

**Built.**  `Transient.radau_early_abort` (default 2; 0 off; else >= 1).  The stepping loop arms an adaptive Radau
attempt whose estimate the step makes (`_SteppingLoop.radau_arm` -> `tr._radau_arm`, a `_RadauArm`) while a rejection
would leave a full attempt in the point's budget, above `minstep`, as the multistep abort is armed.
`_rk_step_transformed` runs the transform's first iteration alone -- the C's (`_tran_radau_tc.solve(..., probe=True)`:
its False the loop's own, the iterate read in its buffer) or the Python loop's -- and, not converged there,
`_radau_early_check`: `_radau_error_estimate`'s formula on the step's frozen C at x_n and its real solve, made in the
reduced space, `f0 = -(i(x_n) + u(t_n))` once a point (its retries share it), normalised by `_StageSteps.normalised`
(`judge`'s, extracted).  Over the threshold the attempt is abandoned (`_RadauAbort`) and `take` rejects it by the
family's retry law (`_StageSteps.reject_early` -> `_retry_after`, `judge`'s, extracted: the worst unknown, the order
two tries measured) -- no iterations after the first, no step end, no dense fallback; counted in `rejected_steps`
and `early_aborts`.  Else the loop goes on from the iterate: `_tran_radau_tc.resume` (no readiness checks and no
source call -- nothing the C reads moved; the iterates one call makes), or, handed back, the Python loop from the
iterate.  `test_radau_early_abort`: abandoned attempts on the latch (none forced, fewer fallbacks and iterations);
a check that never fires is the run without it byte for byte; the C and the Python loop the same run; `judge`'s retry
law; the arm's budget; a threshold under 1 refused.  Pinned pairs re-recorded after checking their twins (`radau_tc.c`
unchanged; loop.c runs the multistep family alone): "the stepping loop and its C", "... and its unread Jacobian",
"the transform Newton and its C", "the stage paths' passes and the core's", "the step's end passes and the C's".

**Measured** (line written first: plan 'THE EARLY ABORT FOR RADAU, PROTOTYPED'; Radau's default transient on the 20
decks, pinned, today's tree against this one, fresh processes, three rounds interleaved, best of 3).  (a) With the
threshold at 1e300 (every check made, none firing) every waveform bit-identical.  (b) Summed 363.99 -> 245.80 s
(-32.5 %): voter25 -71.5, todd3 -58.9, rich3 -36.4, gm19 -33.7, slowlatch -31.6, latch -30.5, gm1 -29.0, mux8 -24.1,
arom -16.2, jge -6.9, g1310 -6.1 %; the transform's dense fallbacks 736 -> 220.  (d) No gold column flips pass/fail,
no forced step, no run fails.  (c) FAILED: mike2 +6.8 % (0.57 -> 0.61 s, steps 614 -> 656, 31 aborts); every other deck
within +1.4 %.  The line's escape for (c), a lottery of GAMMA_REAL +-2 ulp on today's tree, was a dead instrument --
since 157 the transform path's estimate does not read GAMMA_REAL (mike2's digest identical at -2, 0, +2 ulp).  Read for
the decision: mike2's runs part at its point 55, the abort's retry 1.2e-5 apart from the judged one (the estimate at
the first iterate, not the last); today's tree under reltol * (1 +- 1e-5 .. 1e-4) spans 614 .. 651 steps, 0.575 ..
0.639 s, where reltol +-2 ulp leaves it at 614 every time: 614 is its low end.  The first build's check cost ~110 us
and its two C calls ~26 us (ring +7 %); rebuilt before the A/B (`resume`; the lean check): ring +0.7 %.  Andreas's
decision, offered refusal or acceptance: accepted despite the line.  Every gold deck's default transient (Gear2) and
operating point bit-identical to before (54 runs).

**Gate**: suite 5 585 passed, 6 skipped, 3 xfailed.

## 2026-10-21 — 160. DC from a caller's x0: PCNR seeds its limited unknowns unclamped first

**Why** (a bug report from a user project, bisected to 327df7e3, PCNR by default).  `DC(circuit).solve(x0=x)`
from a settled transient's last point -- a current-sense amplifier driving a P-FET into an LED load at 0.25 A -- returned
vout -189 V on a 3.3 V circuit, load -1e-20 A, as a success (path 'jgmin'); the chain (`pcnr=False`) finds 1.655 V.
`DC` passed the caller's `x0` with no limited unknowns, so `pcnr.solve_dc` seeded them through `v_lim_init(limit=True)`,
which limits each junction voltage as a step from 0 V -- built to tame a wild start (a 20 V uniform start overflowing
the first Jacobian), and documented as never for a converged point.  It moved two of the 34 seeds, both in the load:
the LED 1.155 -> 0.098 V and the opto-transistor's base-emitter 0.935 -> 0.093 V.  The load switched off, the plain
solve and the gshunt ladder failed, and the junction-gmin ladder came down to a far root.  A DC sweep's continuation,
which hands each point the previous one's answer, went through the same clamp.

**Built.**  `pcnr.solve_dc_continued(raw_seed=True)`: where the clamp would move `x0`'s seed, the plain solve first
from the UNCLAMPED seed (`v_lim_init(limit=False)`); failing, the plain solve from the clamped seed and the ladders as
before.  `DC` asks for it wherever the caller gave `x0` (its plain solve and a nodeset's held solve).  `x0=None` (zeros:
the clamp is inert) is unchanged, and so is every operating point the transient and the gold census make.  A wild
start's unclamped solve fails at once (the 20 V BJT mirror: singular at its first iterate, 2 ms) and the clamped one
follows, its answer as before.  `test_dc_pcnr`: from its own operating point a diode at 1 A converges at the first
iterate on the plain path (the clamped seed, 0.774 -> 0.088 V, takes three), `DC(x0=op)` returns it; the 20 V mirror
through `DC` still solves, both instance orders.  Pinned pair "the DC base solve and its C" re-recorded (`DC._solve_dc`
moved in its PCNR branch alone; the chain's C, its twin, unchanged).  The reporter's circuit is not in this repository
(its macromodels are the project's); the fix is to be confirmed there.  Not built: rejecting a PCNR answer far outside
the supplies (the -189 V root) -- a boost converter or an inductive kick is outside them legitimately; for Andreas.

**Gate**: every gold deck's default transient (Gear2) and operating point bit-identical (54 runs); suite 5 587 passed, 6 skipped, 3 xfailed.

## 2026-10-21 — 161. DC's acceptance check: an answer must solve the circuit's equations, not only pass the Newton test

**Why** (160's report, item 2; Andreas: "Then investigate item 2", then "YEs, do it").  The user project's diagnostic
on its real circuit (a current-sense amplifier regulating a P-FET into an LED load) at 780882fc: PCNR's ladders
returned a point with the whole supply side near 1.37e10 V on a 3.3 V circuit as a success -- `DC()` from zeros returned
1.18e10 V (gshunt ladder) where the chain finds 1.655 V.  It is NOT A ROOT: its KCL residual, the devices at its node
voltages, is 5.6e79 A.  It passed because the Newton tests are per row and relative, ``|F| < reltol (|J||x| + |F|) +
abstol``, and at 1.4e10 V ``|J||x|`` is ~2.4e91 (relative residual 2.3e-12).  The test is shared: the chain, started
there, accepts it too.  The gshunt rungs show the jump: 3.07 V at gshunt 5.6e-4 and 3.2e-4, then 1.36e10 V at 3.2e-5,
kept down to the pure solve.

**Built.**  `pcnr.dc_residual`: the largest ``|i(x) + u_dc|`` at an answer, the devices at its node voltages (a stateful
limiter synced for the read and put back; the sources read first), over the rows a solve tests;
`pcnr.DC_ACCEPT_RESIDUAL` = 1 (A on a node row, V on a branch row: a true root at reltol 1e-6 exceeds it only on a row
carrying over 1e6 A).  `pcnr.solve_dc_continued` checks every answer it would return -- the plain solve's, each
ladder's pure solve -- and one that does not solve counts as that route failing: the next is tried, and where all fail
`DC` falls back to the chain as before.  `DC._accept` checks the chain's answer (and a pinned solve's, its held rows not
tested): one that does not solve raises `NoConvergenceError` naming the row.  `test_dc_pcnr`: a planted non-solution
from PCNR's plain solve is refused and the next route gives the operating point (`solve_dc_continued` and `DC()`); one
from the chain raises.  It also refuses what the noncon check exists for -- a limiter holding a device short of its
node voltages, `_limiting.NONCON_CHECK` off: a 1e28 A residual (`test_noncon`: the old test runs with the acceptance
check off too, a sibling pins the refusal).  The other tests it touched: the ladder-order tests' fake circuits stub
it (`test_pcnr_jgmin`); the coupled-path PCNR gate does not count its read of the devices as limiting (`test_pcnr`); the
source stamp's arming bound measures its subject -- a few arms against ~105 served calls, `<= 4` -- where `<= 2` held
only in the suite's order (alone it arms 3 times at 333fe000; the check's read after PCNR's solve makes 4).  Pinned pair
"the DC base solve and its C" re-recorded (`DC._solve_dc`'s tail, after the base solve; its C twin unchanged).

**Refused: the far-voltage warning** (a node beyond 1e3 times the largest voltage source, warned and kept; in the same
line, which asked the suite to pass).  The suite failed 13 tests on it: the internal state nodes of behavioural models
are not voltages -- a VCO's phase integral (1.05e4 against a 0.5 V source), an op-amp model's gain node (1.59e4
against 15 V), a macromodel's mid node (-1 V against a 1e-5 V input) -- and an ideal open-loop amplifier (3.8e18 V) and
a gmin-anchored floating node (5e8 V) reach it by design.  Not built; the acceptance check is the fix.

**Measured** (line written first: plan 'PCNR'S FAR ROOT').  Every gold deck's default transient (Gear2) and
standalone operating point bit-identical to 780882fc (54 runs); no gold deck's operating point fails the check.  On
the reporter's circuit -- not in this repository, its macromodels are the project's -- the confirmation (csa_pfet from
zeros at 1.655 V) is the reporter's to make.

**Gate**: suite 5 590 passed, 6 skipped, 3 xfailed.

## 2026-10-21 — 162. A lottery that is alive: each deck's own noise, for the lines that cite it

**Why** (Andreas: "Do as you suggest", then "Commit"; item 2 of the list after 158).  159's line had an escape for a
deck slower than +3 %: a lottery of GAMMA_REAL moved by +-2 ulp on today's tree, slower than its slowest arm.  It was
DEAD -- since 157 the transform path's estimate does not read GAMMA_REAL; mike2's digest was the same at -2, 0, +2 ulp
-- so it excused nothing, and the line failed on mike2 +6.8 %.  And a live last-bits lottery can still be flat: mike2
under reltol +-2 ulp took 614 steps every time, every digest different, while a 1e-5 nudge -- the scale 159's abort
moved a retry -- spans 614 .. 651 steps and 0.575 .. 0.639 s.

**Built.**  `benchmarks/lottery.py`: per deck the base run and four arms with `reltol` scaled by 1 -+ 1e-5 and
1 -+ 3e-5; an arm counts only where it is ALIVE (its waveform's digest differs from the base's), a deck with none is
DEAD; each run a process of its own (one BLAS thread, the allocator pinned, as the A/B lines run), rounds interleaved
with the order reversed every other round, best of the rounds, under the benchmark lock, the conditions stamped.
`table` reads the runs (raising where a base or an arm gave two digests: no lottery in a run that is not
deterministic); `bound` is the rule a line cites -- a deck fails only if slower than its slowest live arm by more than
3 % (a dead deck: than its base).  `test_benchmark_lottery`: the liveness, the best of the rounds, the span over the
live arms, the bound, a dead deck, a run that is not deterministic, a failed run left out.

**Measured** (a740330e, Radau's 20 decks, 2 rounds, pinned; the table rebuilt by `table` from the same 200 runs
equals the one read when they were made): 19 decks alive at all four arms, AB_INTEG DEAD (its digest unmoved by
reltol at 3e-5).  Over the base and the live arms: arom -4.0 .. +18.3 %, todd3 -4.4 .. +11.3, gm19 -7.8 .. +8.4, gm3
-1.3 .. +5.4, gm17 -9.6 .. +3.9, latch 0 .. +3.0, gm1 -31.9 .. 0 (its base its slowest), the others within -6.5 ..
+2.2 %; the 20-deck sum -3.4 .. +2.9 %.  So a per-deck "+3 %" fails on noise alone on five decks, and a summed
criterion means something only beyond 3 %.  The table describes the tree it was measured on: measured again where a
change moves the step controller.

**Gate**: suite 5 594 passed, 6 skipped, 3 xfailed.

## 2026-10-21 — 163. PCNR's last route: pseudo-transient continuation; opampal's operating point, its path measured

**Why** (Andreas: "Can we make PCNR more robust to handle these cases without falling back", then "Do 1 and 2 now",
"Finish 2", "Finish 1").  On a user's converter (161) PCNR's three routes all fail from zeros and `DC` falls back to the
chain; the chain ends with pseudo-transient continuation (`nrsolver.PseudoTransientNewton`), PCNR had nothing like it.

**Built.**  `pcnr._ptc_ladder`, `solve_dc_continued`'s LAST route (path 'ptc'): the chain's schedule
(`nrsolver._adaptive_conductance_ladder`, g from 1 down to 1e-12 by halving steps in exponent space, escalating to
1e6 before the first rung converges), each rung a gshunt rung ANCHORED at its seed, ``g (x - x_seed)``, the limited
unknowns carried from rung to rung, then the pure solve and the acceptance check (161).  `solve_dc(anchor=)` carries
the anchor in the sources the loop is handed (``-g anchor`` on the node rows), so the C loop takes it as the Python
loop does -- no C change.  `test_dc_pcnr`: with the plain solve and the three ladders failing, the route alone solves
the diode at 1 A from zeros and the BJT mirror from its 20 V start to the chain's answers (and `DC()` reports 'ptc');
anchored at the solution a rung solves to it, anchored at 0 V it is the plain gshunt rung byte for byte.  The ladder
test that forces every solve to fail now sees the route's rungs after source stepping.  Pinned pairs "the PCNR DC
solve and its C" and "the DC base solve and its C" re-recorded (the anchor rides in the sources the C is handed; a
comment in `DC._solve_dc`).

**Measured, and decided** (plan 'PCNR'S MISSING ROUTE'; line written first).  Sized ALONE on the 27 gold decks (the plain
solve and the ladders skipped): 24 solved, often in fewer rungs (voter25 13 of today's 18), but only 5 of the 7 decks
today's ladders solve (not e1480, opampal), and slowlatch and todd3 to the OTHER state of their latches (it follows the
pseudo-time trajectory from the zeros seed): its line's clause "solves at least what the ladders solve" FAILED.  As the
LAST route no gold deck reaches it: every gold deck's default transient and operating point bit-identical (54 runs).
Andreas's decision: build it last.  Whether it rescues the user's cases is the user's to measure.

**opampal's operating point** (the open item "sensitive to last-bit changes").  Its ANSWER is robust: one operating
point from six random starts (+-1 V; PCNR the base bit for bit, the chain within 3.5e-12 V), a resistor x (1 + 1e-12)
moving it 2.5e-12 V.  Its PATH is not: under 12 live nudges (each of 12 resistors x (1 + 1e-12)) 916 .. 1 841 Newton
iterations, 0.28 .. 0.41 s, 29 .. 46 rungs, 9 failure patterns -- the cost is 11 .. 16 rungs failing AT THE ITERATION
LIMIT, which ones moving with the last bits.  No exit rule tells them from rungs that converge (557 solves over the 27
decks and the nudges): converged rungs take up to 93 iterations (no cap is free); a residual plateau of 30 iterations
cuts 13 that converge; an exact cycle catches 1 of 141; a rate prediction cuts 20 .. 42 -- PCNR rungs converge by a
sudden drop after a long plateau or rise (the limiter's slew: a thermal voltage a step).  Closed without a fix; its
time is covered by 162's lottery rule.

**Gate**: suite 5 596 passed, 6 skipped, 3 xfailed; every gold deck bit-identical (54 runs).

## 2026-10-21 — 164. Radau's coupled stage C: an underflow no longer hands the solve back

**Why** (Andreas: "Fix Radau's C solve flag bails").  `radau.c`, the C twin of the coupled stage Newton
(`_stage_newton`), still tested PCX_FP_FLAGS -- underflow included -- after 156 took underflow out of `radau_tc.c`'s
test: an underflow numpy's default error state ignores handed the solve back and the Python Newton repeated it.  The
same rule was refused twice on time, not on its answers: in 156 (arom +5.4 %, the allocator's mode) and in 159's item 4
(pinned, three rounds: +0.3 % summed, every deck -2.7 .. +1.5 %, its line asked the sum to fall).  Built now on
Andreas's decision; `radau.c` now tests what numpy reports, as `radau_tc.c` (156), `dcnewton.c` and `pcnrdc.c` do.

**Built** (the patch measured in 159, unchanged on today's tree).  `radau.c` tests overflow, invalid and divide
(`RC_FLAGS`) where it tested PCX_FP_FLAGS; `_tran_radau_c.solve` declines where numpy's error state does not ignore
underflow (`radau_c:underflow`, `func._underflow_ignored`).  `test_radau_c`: a 1e300-ohm resistor across a load (a
1e-300 S entry, the blocks' products underflow) is served, the Python Newton's stages and memos bit for bit, and under
`errstate(under='warn')` declined with the Python's warning -- failing on the old C (`radau_c:bail:flags`).  Pinned pair
"the coupled stage Newton and its C" re-recorded (the twin; the reference did not move).

**Measured** (Radau's default transient on its 20 decks, today's tree against this one, one round, pinned): every
waveform and step count bit-identical; every flag hand-back of the coupled C now served (voter25 107 solves against 17,
gm19 20 / 1, arom 20 / 0, slowlatch 18 / 2, latch 19 / 8, rich3 16 / 1, gm1 2 / 1, todd3 14 / 13, jge 1 / 0), none
left.  Time: 159's three pinned rounds stand (+0.3 %, within noise); this round's jge 167 -> 87 s is that deck's own
spread (87 .. 144 s in 159's rounds, either tree), one coupled solve apart.

**Gate**: suite 5 597 passed, 6 skipped, 3 xfailed.  The gold census does not run Radau: not run.

## 2026-10-21 — 165. PCNR's pseudo-transient route powered up and stepped by progress

**Why** (Andreas: "Can the Pseudo-transient be improved?", "push the size items 1, 2 and 3", "Do it", then "take
power-up combined with adaptive steps despite its line").  163's route ran the chain's schedule at full sources:
alone on the gold decks it solved 24 of 27, not e1480 and opampal, and took slowlatch and todd3 to the other state of
their latches.

**Measured first** (plan 'IMPROVING THE PSEUDO-TRANSIENT ROUTE', lines written first; each variant ALONE on the 27 gold
decks, PCNR's other routes skipped, against today's PCNR answer).  Three changes sized alone: the circuit's own
capacitance as the pseudo-term (left 4 unsolved, gm17 to its other state); the sources POWERED UP with the anchored
conductance (2 unsolved, no state changed); the steps BY PROGRESS, switched evolution relaxation (1 unsolved, slowlatch
to its other state).  Power-up and progress together (both ride in `solve_dc`'s inputs: no C change): 26 of 27 (jge
falls back to the chain; progress alone solved it), every answer within DC's tolerance of today's PCNR (worst 0.85 of
it), NO deck in another state, e1480 and opampal solved; 2 062 Newton iterations against the old route's 1 976 on the
decks both solve (+4.4 %); 73.2 against 79.7 s.  Its line asked all 27 and no more iterations: FAILED on both.
Andreas's decision: build it despite the line.

**Built.**  `pcnr._ptc_ladder`, the same anchored gshunt rungs (`solve_dc(anchor=)`), now: each rung's sources at
`_ptc_source(g)` -- `PTC_SOURCE_MIN` (1e-3) from g = 1 up, in full from `PTC_SOURCE_FULL` (1e-6) down, the logarithm
linear in g's between; after a converged rung the next g is its g times the ratio of `dc_residual` at its answer (full
sources) to the one before, clipped to [`PTC_RATIO_MIN`, `PTC_RATIO_MAX`] = [1e-3, 0.5]; a failed rung retried halfway
in exponent to the last converged, giving up within 5 % of it; before any rung converges g escalates by 100 up to 1e6;
from `PTC_G_END` (1e-12) the pure solve; 60 rungs at most.  The sizing re-run with the built route reproduced the
prototype's on all 27 decks (answer, rungs, iterations, failures).  `test_dc_pcnr`: the schedule on the 1 A diode (each
rung's sources, each next g from the residuals, the pure solve after `PTC_G_END`); the forced-failure sequence's
route is now its four escalating rungs at a thousandth of the sources; the route alone still solves the diode and
the 20 V BJT mirror to the chain's answers.

**Gate**: suite 5 598 passed, 6 skipped, 3 xfailed; every gold deck bit-identical as the last route (54 runs).

## 2026-10-21 — 166. The pseudo-transient route's second pass at full sources; its steps during the power-up, refused

**Why** (Andreas: "What more can we improve in the pseudo transient", "Do as you suggest", then "Build it").  Alone on
the gold decks 165's route takes 23 of 27 decks in exactly 24 solves -- its schedule's minimum: the residual its steps
read is at full sources and cannot fall while the sources come up, so the ramp is all halvings -- and leaves jge to
the chain.

**Measured first** (plan 'THE ROUTE'S STEPS DURING ITS POWER-UP, AND JGE', lines written first; each variant alone on the
27 decks against 165's route and today's PCNR answers; the instrument -- a scratch copy of the route -- reproducing
165's on 27 of 27 decks).  (1) The steps' progress read at the rung's OWN source level, two readings: both residuals at
the new level -- 285 rungs against 664, but 3 586 Newton iterations against 2 768 (+30 %) and gm17 in its other
state; each at its own level -- hussamp unsolved, gm17 in its other state.  REFUSED: the halving ramp is what keeps
the states.  (2) jge traced: every rung converges, in both of DC's passes; the route ends at its pure solve, failing
at the iteration limit.  A guess that the ramp's last two steps (each the largest, 1e-3: the residual before them
read at full sources from an answer made at 0.72 of them) caused it is REFUTED: a halving there and honest steps after
still leave jge unsolved (its iterations -8.9 % on 26 decks, every deck's path moved; refused on its line).  A SECOND
PASS at full sources where the powered-up one fails -- the steps-by-progress variant that solved jge alone -- meets its
line: jge solved in 0.5 s (72 s of failing and the chain's fallback before), within DC's tolerance of today's PCNR and
in its state; the other 26 decks' answers, rungs and iterations identical to 165's.

**Built.**  `pcnr._ptc_ladder(powerup=False)`: every rung at full sources; `solve_dc_continued` calls it where the
powered-up pass returns nothing or a point the acceptance check refuses (path 'ptc' either way).  Second, not first:
at full sources alone it takes slowlatch to its other state and does not solve opampal.  The sizing re-run with the
built route reproduces the variant's on all 27 decks.  `test_dc_pcnr`: the second pass at full sources where the first
returns nothing and where its point is refused (the 1 A diode's operating point); the forced-failure sequence ends in
both passes' four escalating rungs, at a thousandth of the sources, then full.

**Gate**: suite 5 600 passed, 6 skipped, 3 xfailed; every gold deck bit-identical as the last route (54 runs).

## 2026-10-21 — 167. A PCNR rung that jumps far fails as a rung: a user's converter solved by PCNR's own routes

**Why** (Andreas: "Can we make PCNR more robust to handle these cases without falling back", "Ask the peer to send the
circuit", then on the cases "Do as you suggest" and "Do it").  A user's three far-root circuits arrived as a generated,
self-contained script (scratch only: vendor model numbers).  On it the CSA + P-FET converter from zeros raises with PCNR
on and off: PCNR's routes all end refused by the acceptance check (161) or failing, and the chain converges to a
non-solution (residual 2.5e62).  Traced (every rung's deformation, iterations and largest node voltage): in the gshunt
ladder, in source stepping's gshunt rungs and in the pseudo-transient route's full-source pass alike, ONE rung near
1e-5 S goes from ~3 V to 6.8..12 kV after 46..98 iterations and "converges" there by the relative Newton test; every
later rung and the pure solve stay at that far point, refused only at the end.  The pseudo-transient route (163, 165,
166) rescues none of the three cases.

**Sized first** (plan 'PCNR ROBUST ON ITS OWN', its line written before the cases arrived: (a) PCNR's own routes solve
the three cases; (b) every gold deck's standalone operating point the same, none falling back, no gold column flipping,
the 27 decks' summed DC time no worse than +3 %; (c) the suite).  Rung continuity -- a rung moving more than JUMP x
(its start's largest node voltage + the largest DC source) counts as failed -- sized on the gold decks at JUMP 1e3
before the cases (bit-identical; the largest in-ladder move 0.50 of that scale, jge).  On the cases at 1e3 the
converter still falls back (retried rungs jump just under it); at 100 and at 10 PCNR SOLVES IT by source stepping, the
expected 1.65518481 V, in 0.29 s against 1.93 s of failing; the other cases unchanged; the gold decks' operating points
bit-identical at both.  100 chosen after 1e3 failed (the line's clauses unchanged): 200x above the gold decks' largest
move, under the converter's jumps.

**Built.**  `pcnr._rung`: handed a previous rung's limited unknowns (inside a ladder, not its first rung), a rung whose
answer moves more than `RUNG_JUMP` (100) x (max |x| of its start + the largest DC source) over the nodes is a failed
rung (`_jumps`): the ladder retries smaller.  Every route's rungs go through it (the gshunt and junction ladders,
source stepping, both pseudo-transient passes, the pure solves after a ladder).  The sources are read only where the
move exceeds `RUNG_JUMP` x max |x| -- never on the gold decks.  No C change.  `test_dc_pcnr`: the jump's scale (from 0 V
the source alone, from 5 V both; the sources unread below it); a jumping rung is None, a first rung is not judged; the
gshunt ladder on the 1 A diode, its 1e-3 rung planted to jump 1e9 V, retried at the geometric mean and coming down to
the diode's operating point -- each failing without the check.

**Gate**: (a) met -- the converter solved by PCNR's source stepping (1.65518481 V), the other cases unchanged;
(b) the 54 gold runs (each deck's default transient and standalone operating point) bit-identical, none falling
back; the 27 decks' summed standalone DC time +0.34 % (0.8553 -> 0.8581 s, best of 3 interleaved rounds, pinned,
a quiet box; single decks -0.8 .. +2.2 %, every answer and path counter the same); (c) the suite, 5 603 passed.

## 2026-10-22 — 168. Eager nonlinear HDL classes and DC-pinned ones on the evaluate core: one class without `var()` no longer takes a circuit off C

**Why** (Andreas: "Size the auto-chaining for unchained Behavioural classes", then "Try to get meristor setved by the C
core", then "Build both once the suite passes").  The transient's evaluate core (`_tran_core`) serves a non-constant
HDL element only if its class is chained -- compiled through the let-chain because its `analog()` declares a `var()`
-- and C-bound; one element that is not makes the whole circuit Unservable, and every step runs in Python.  A user's
buck converter: one comparator output class without `var()`, 20.2 s; with `va = var(ba.V, 'va')` 2.9 s (the user's
session measured it).  Nothing said so.  In the library three classes were eager and non-constant -- `DiodeHdl`,
`MixerHdl`, `MemristorHdl` -- and the memristor was refused a second time for its DC pins.  No gold deck holds one (all
27 are served today).

**Sized first** (plan 'AUTO-CHAINING EAGER BEHAVIOURAL CLASSES', the line written before measuring: (a) one small
transient each for the three classes, the user's comparator and a 6-stage Cockcroft-Walton multiplier of 12
`DiodeHdl`, served by the core, summed at least 2x faster, none slower; (b) the moved classes keep `eval_i_pure`,
`solve_batched`, the JAX path and the symbolic forms, every other library class compiled the same, the 54 gold runs
bit-identical; (c) each moved circuit moved no more than a 1e-12 V nudge of its start moves it on today's code; (d)
cold compile within max(+1 s, 2x), warm import unchanged, the suite).  ⚠ hdl.py's own rule: NEVER turn an eager model
into a chained one -- a chained model has no `eval_i_pure` and loses `solve_batched` (20.7x at 512 lanes) and the JAX
path.  So the design is dual: the numeric functions through the chain, the eager pure and symbolic forms kept.

Measured (best of 3 interleaved rounds): rectifier 871.9 -> 10.2 ms, multiplier 1 527 -> 30.0, mixer 33.6 -> 2.8,
comparator 356.5 -> 6.5, memristor 262 -> 193 (not served: its DC pins), summed 3.05 -> 0.24 s; 51 of the library's 54
classes compile to the same sources, the three keep `eval_i_pure` where they had it; the 54 gold runs bit-identical;
the rectifier and mixer bit-identical, the comparator 8.0e-12 V against its nudges' 8.2e-12 and 1.3e-11, the
memristor 3.2e-15 against 2.0e-12; the multiplier NOT EVALUABLE as written -- under both +-1e-12 V nudges today's code
fails (SingularMatrix) -- on the prototype a nudge moves it 8.33 mV and the change 8.38 mV; cold compile + first C build
+0.04..0.11 s a class, the warm import 1.119 -> 1.100 s.  Then the memristor: the core takes a class with DC pins and
declines a DC epar; served, bit-identical to the Python path, 193 -> 9.5 ms.

**Built.**  `hdl.generate_code`: a model with no `var()` whose G or C depends on the unknowns compiles its numeric
functions through the let-chain (no definitions) -- the C kernel, the core -- with `pure_spec` and `sym_spec` still the
eager forms; `_hdl_info['auto_chain']`; a model with constant G and C stays eager.  Its methods take the chain except on
a symbolic toolkit, which keeps the eager symbolic forms.  `_tran_core._lone_batches` no longer refuses a class with DC
pins -- it evaluates its pinned functions only where the epar's `analysis_kind` is 'dc', and the core serves a
transient's evaluations; the core records `pins`, and `probe` -- and its twin `passes_held`, for Radau's C --
declines a DC epar ('pins_dc').  `_hdl_batch.split` still keeps such a class out of its batches (they serve DC too).  `test_tran_core`: `DiodeHdl` chained and auto-chained,
keeping `eval_i_pure` and its symbolic forms, `RHdl` eager; a `DiodeHdl` rectifier served, bit for bit with the core
off; a memristor served, bit for bit, and a DC epar declined -- each failing on the old code.  Six tests that told the eager
code generator by `chained is False`, or took a flat class for one without C kernels, now read `auto_chain` or pin
a class to numpy; the sequential-limiter test switches the C limiter off (built from the spec as the class bound,
it cannot see a spec reversed in place).  Pinned pairs "the transient evaluation and the core" and "the step's end
passes and the C's" re-recorded.  Side findings, not
fixed: `DiodeHdl` has no limiting, and a transient at a 5 V or 2 V drive (or 1 V with the default first step) goes
numerically SINGULAR and ABORTS rather than retrying the step smaller, on the old code as on the new; nine hand-written
nonlinear elements (`Diode`, `VSwitch`, `ISwitch`, `BSource`, `NonLinearVCCS`, `VCVS_limited`, three `Idtmod`) still
take a circuit off the core.

**A correction to 167.**  The user's far-root script was not faithful to their circuits (the user's session found
it: the nodes created in another order, which steers the solver's paths from zeros, and every parameter written as
given, which moves a MOSFET model that reads givenness -- the rebuilt P-FET differed by up to ~1 A).  On the
corrected script HEAD's PCNR solves the converter from zeros by the JUNCTION-GMIN ladder (1.65518481 V), where
a740330e fails it with PCNR and without; 167's result holds, its route was jgmin, not source stepping.  The user's
other circuits (two op-amps, an inverting and a referenced current-sense, a low-side sink, a current monitor, the
buck) solve from zeros with PCNR, and all thirteen are served by the evaluate core.

**Gate**: (a) met but for the memristor in the first stage (then served); (b) met; (c) met where evaluable (the
multiplier's as above); (d) met; the 54 gold runs on the built tree bit-identical; the suite, 5 626 passed.

## 2026-10-22 — 169. A C Newton's iteration limit raised from the C: a failing step no longer repeats its iterations in Python

**Why** (Andreas: "Continue improving speed", and the user's session's circuits).  The user's buck converter (20 us
of its own testbench), served by the evaluate core and the C stepping loop, profiled: 4.9 s of 6.3 s were 66 PYTHON
Newton solves of ~98 iterations each -- the steps whose Newton fails.  By design (`_tran_newton_c`, TRANSACTIONAL)
any unusual outcome, the iteration limit included, rolled the C attempt back and `_newton` repeated the same
arithmetic in Python, for the same exception, warnings and statistics: a failing step was solved three times (the C
loop's attempt, the C Newton's, the Python's), the last ~10x the cost of a C one.  The DC's C Newton already served
its limit (`_dc_newton_c`: `nrsolver._standard_failure` on its last iteration).

**Sized first** (plan 'A C NEWTON'S MAXITER FAILURE', the line written before measuring: (a) the 27 gold transients
and the buck bit-identical, no warning lost; (b) the buck at least 2x faster; (c) no gold deck slower than its
lottery bound, the summed time no worse; (d) the suite, and a test failing on the old code).  Measured: (a) the 27
gold transients and the buck bit-identical (digest, steps, rejections, Newton count), no warning on either tree; the
limit is never raised on a gold deck (their remaining maxiter hand-backs, slowlatch 2 and todd3 10, are chord solves,
which fall back to the full Newton); (b) the buck's 20 us 1.044 -> 0.337 s (x3.1, best of 5 in a process, quiet box; x3.06 interleaved under load);
(c) the 27 gold transients' summed time -0.03 % (one solve a process, 3 rounds interleaved, under a falling load);
the two it read over +3 % -- mosrect +8.7 %, schmitecl +11.5 %, a cold single solve of 25..50 ms, mosrect's lottery
dead -- re-timed warm, best of 10 in a process, against an A/A tree: mosrect 29.03 / 29.08 (A/A) / 29.00 ms,
schmitecl 8.73 / 8.73 / 8.74 ms.

**Built.**  `_tran_newton_c.solve`: where the full Newton stops at `maxiter` with every iterate finite, and no damped
retry or rescue is armed, StandardNewton's exception is raised from the C's last iteration (`_unconverged`,
`_maxiter`: `nrsolver._standard_failure` on its residual, limited step and iterate, read as `_newton` reads it -- a
message naming a singularity promoted to SingularMatrix; `newton_c:maxiter`), the attempt's trace kept as the repeat
leaves it and the early-abort arm popped as `_newton` pops it.  A chord's limit still hands back
(the Python chord falls back to the full Newton).  `test_newton_c`: a transient whose `maxiter` is 1 (its operating
point given) raises 31 limits from the C, its answer, statistics and warnings the Python's -- failing on the old
code; the two bail tests that faked code 11 after a converged solve fake a non-finite iterate (still a hand-back).
Pinned pairs "the Newton solve and its C" (its twin now names `_unconverged` and `_maxiter`) and "the stepping loop
and its C" re-recorded.
What remains on the buck (6.3 -> 1.36 s under cProfile): the C loop 0.69 s and the C Newton 0.42 s, ~0.27 s of it
the C Newton re-solving each step the C loop hands back on its Newton failure -- a next stage, the loop rejecting
such a step itself.

**Gate**: (a) met; (b) met (x3.1); (c) met (summed -0.03 %, the two flagged decks flat against their A/A); (d)
met; the suite, 5 627 passed.

## 2026-10-22 — 170. The C stepping loop rejects its own failed Newton step: a failing step no longer leaves the C

**Why** (Andreas: "Build the C loop and try to improve", after the sizing).  On the user's buck converter (20 us of
its testbench) after 169, every one of the C loop's 88 hand-backs was the C Newton at its iteration limit (100
iterations): the loop stopped, `run` returned and `execute` entered again at once -- 37 times the same attempt failed
again in C, 0 attempts -- and then the Python attempt solved it a third time to raise the limit from the C and reject
it (`take`: a quarter of the step).  Failing steps cost 15 300 C Newton iterations against 7 975 converged.

**Sized first** (plan 'THE C LOOP REJECTS ITS OWN FAILED NEWTON STEP', the line written before measuring: (a) the 27
gold transients and the buck bit-identical, every count but the loop's own the same; (b) the buck at least 1.3x
faster; (c) no gold deck slower than its lottery bound, the summed time no worse; (d) the suite, and a test failing on
the old code).  Measured: (a) the 27 gold transients bit-identical (digest, steps, rejections, Newton count) and their
loop counts too -- the rejection never fires there (their maxiter hand-backs are chord solves); the buck bit-identical
(061ac43ff24c, 1 730 steps, 593 rejected, 7 976 Newton), every count but `loop_c:*` the same, no warning on either
tree; (b) the buck's 20 us 0.338 -> 0.187 s (x1.80;
best of 5 in a process, 5 rounds rotated, quiet box; an A/A tree of the same commit read -0.4 %); (c) every gold deck
-0.5..+0.8 % (the A/A tree -0.7..+0.7 %), the 27 decks' summed time -0.05 % (A/A -0.20 %), every digest the same.

**Built.**  loop.c's `take`: where the C Newton stops at its limit (`NEWTON_BAIL_MAXITER`) and the run allows it
(`newton_reject`: adaptive, the full Newton, no damped retry or rescue armed -- where `_tran_newton_c._unconverged`
would raise the limit) and a quarter of the step is not under `minstep`, the attempt is rejected in C
(`LOOP_JUDGE_NEWTON`, `newton_reject`) as `_SteppingLoop.take` rejects the exception: the step a quarter, `imposed`
where the step was held (`clamped`, now set by `where`), the point's flags kept as a judged rejection keeps them, a
rejection counted, the seed's `pred:fast`, nothing of the attempt recorded as made, its iterations not counted.
`imposed` is carried in and written back (judged rejections, early aborts and accepts clear it, as the Python's do).
The glue counts what the Python's failed attempt counted (`newton_c:maxiter`, the source pass's counts) and
`loop_c:newton_rejects`.  Below `minstep` (the rescue), on a fixed grid (the fallback's warning record), a chord's
limit (the full Newton after it) and every other status still hand back.  `test_loop_c`: an adaptive run at
`maxiter` 2 rejects its failed steps in C (`loop_c:bail:newton` 0), its bytes, statistics, warnings, counts and state
the Python loop's -- failing on the old code.  Pinned pair "the stepping loop and its C" re-recorded (the twin
`_SteppingLoop.take` unchanged; its reason names the rejection).

**Gate**: the suite, 5 628 passed (6 skipped, 3 xfailed); the 27 gold transients bit-identical on the build
(digest, steps, rejections, Newton count, every count), the buck bit-identical (`loop_c:newton_rejects` 65,
`loop_c:bail:newton` 0).

## 2026-10-22 — 171. The sparse path's verdicts made sparse: the error test and the Newton stop forming dense products whose only use is a yes or no

**Why** (Andreas: "Build the C loop and try to improve").  Profiled with DWARF call graphs (perf), the C stepping loop
is 96 % of the large gold decks' solves (62 % of opampal's), and inside it, on the sparse (KLU) path, four O(m^2)
passes over matrices that are zero off the circuit's pattern each fed only a VERDICT: the error test's sparse solve
reduced the whole J (a memmove), counted its nonzeros (`_ours`'s foreign test) and validated the refactor with a dense
dgemv; the Newton's convergence test formed |J| |x_next| with a dense dgemv every iteration, and the sparse full Newton
validated each refactor with another.  On opampal they were ~70 % of the C loop, on gm19 ~25 %, on jge ~9 % (its
devices dominate).  152 had made the same verdicts sparse in PCNR's DC (`pcx_klu_resid_sparse`).

**Sized first** (plan 'THE SPARSE PATH'S DENSE VERDICTS', the line written before measuring: (a) the 27 gold
transients and the buck bit-identical, every count; (b) the summed time at least 5 % lower, no deck slower than its
lottery bound against an A/A tree; (c) the fallback alive, the sparse verdict never contradicting the dense one; (d)
the suite, tests failing on the old code, pinned pairs re-recorded).  Measured: (a) the 27 gold transients
bit-identical (digest, steps, rejections, Newton count, every count), the buck bit-identical; (b) the 27 decks' summed
transient time -11.4 % (39.67 -> 35.16 s; the A/A tree -0.01 %), best of 5 in a process (best of 1 after a warm-up over
3 s), 3 rounds rotated, quiet box: opampal -25.3 %, gm19 -22.8, arom -17.9, rich3 -16.6, voter25 -11.1, gm1 -9.7,
mux8 -9.6, toronto -8.9, gm17 -7.7, jge -6.5, the rest -6.0..+0.2 % (the A/A tree -0.8..+0.6 %); the buck -7.6 %;
(c) `test_verdict_c`: rows drawn at the threshold, each verdict seen holding, failing and undecided, never against
numpy's; opampal's run left 2 of ~1 000 validations to the dense test, its Newton none.

**Built.**  `csrc/pcx_verdict.h` (new): `pcx_resid_sparse` (152's residual verdict, moved from pcx_kluseq.h, which now
calls it; dcnewton.c and pcnrdc.c include the header) and `pcx_scale_sparse` (nrsolver's residual test: every term
non-negative, so T_i = S_i and the dgemv's row lies in [max(S_i - e_i, 0), S_i + e_i]; one row failing decides).
newton.c: the convergence test makes the step's test first and the residual's only where it holds (its product feeds
nothing else, no flag is read); on a KLU path its verdict by `pcx_scale_sparse`, numpy's dgemv where undecided; the
sparse full Newton's validation by `pcx_resid_sparse`, the dense dgemv where undecided (`vdense` counts them).  lte.c:
where the loop's error test reads its Newton's J -- the core's, on its pattern alone, which lies within the run's KLU
pattern, both one StructuralKLU's (`_tran_loop_c._on_pattern`) -- the sparse solve gathers the values from J by `kat`,
no reduce and no count, and takes the validation's verdict from the bound where it decides and its sums raised no
flag, else the dense test on R scattered from the values (`loop_c:lte_on_pattern`, `loop_c:lte_dense`).  Where the
bound does not decide nothing changes; a decided verdict is numpy's (the premise: test_pcnr_c's
`test_numpys_dgemv_rounds_a_row_once_a_nonzero_term`).  Tests: `test_verdict_c` (both verdicts against numpy's at the
threshold, failing on the old code: no header); `test_loop_c` (the sparse adaptive run's error test on the pattern).
Pinned pairs "the stepping loop and its C", "the Newton solve and its C", "the error test and its C", "the DC base
solve and its C", "the PCNR DC solve and its C" re-recorded (each reference unchanged; their reasons name the
verdicts).

**Gate**: the suite, 5 638 passed (6 skipped, 3 xfailed); the 27 gold transients bit-identical on the build
(digest, steps, rejections, Newton count, every count but the new `loop_c:lte_on_pattern`), the buck bit-identical.

## 2026-10-22 — 172. A sine source's phase and its AC phase made two parameters

**Why** (Andreas: "fix sine bugs", after the exporter's sizing found it).  `VSin` and `ISin` declared `phase` twice:
VS's (IS's) AC analysis phase and the sine's phase in degrees -- one value in the ParameterDict, so a sine's phase
turned its AC phasor too (`VSin(vac=1, phase=90)`: AC source -1j), and an AC phase moved the transient's sine.  The
PSS log (`doc/pss_log_260902.md`, 2026-09-22) met it as an instrument trap ("divide by the perturbed source's AC
phase"); the SPICE importer refused a SIN phase other than the AC phase.

**Built.**  `VSin` / `ISin`: `phase` is the sine's (SPICE's SIN phase, as every caller in the tree passes it); the
AC phasor's phase is `acphase` (VS's `phase` dropped from their parameters), read through `_ac_phase` (VS and IS:
`phase`).  `spice_import`: `AC mag phase` on a SIN source maps to `acphase`, the refusal removed; the deck written
back (`_source_spec`) carries it.  A sine source given `phase` alone now has an AC phase of 0 -- an AC or PAC result
that read the sine's phase as its AC phase changes by that rotation.  Where the code meant the two to be one it now
says so: `ProbeShooting._build` gives its probes `acphase` = `phase` (`_pac_response`'s image sideband is right only
with the excitation in the probe's own phase -- its PAC Jacobian 100 % off otherwise, which its validation caught), as
do three tests whose reference divides by the source's AC phase (`test_shooting`, the staged PAC's three methods).
Tests: `test_elements` (each phase turns only
its own: VSin and ISin), `test_spice_import` (a SIN with AC phase 60 and SIN phase 30 maps both; written back and read
again, both kept) -- failing on the old code.  Pinned pair "the sources and their direct pass" re-recorded (the
reference `VS.u` / `IS.u` changed in the AC branch alone; the twin, the time-domain add, unchanged).

**Gate**: the suite, 5 642 passed (6 skipped, 3 xfailed) -- its first run failed six tests, the probe's PAC
Jacobian and three references that divided by the sine's phase as the AC phase, each now giving `acphase` = `phase`.

## 2026-10-22 — 173. A circuit written out as standalone Python: `export.to_python`, checked against the original

**Why** (the Howland session's proposal, 2026-10-22; Andreas: "build exporter", after its sizing met the line).
Every reproduction a user sent needed a flat circuit written out as Python, and each went wrong on a pycircuit
internal a hand copy cannot keep in step: the node order (the unknowns' order steers the solver from zeros), the
parameters' givenness (a level 3 MOS computes PHI from NSUB only where PHI is not given), a keyword-named parameter
(`as`), a collapse variant (`X_collapseNN`), a `BSourceHdl_<hash>` class, a user class.

**Sized first** (plan 'A CIRCUIT-TO-PYTHON EXPORTER', the line written before measuring; prototype in scratch): (a)
all 27 gold decks round-trip exactly (node and branch order, i/G/C/q at random x, u, the DC answer bit for bit, the 22
transients under 5 s identical); (b) 90 library classes, 180 of 180 cases exact or refused by name; (c) the user
session's 13 circuits exact, its buck's transient identical; (d) the check alive (givenness dropped, first-use node
order, a replaced source function: each caught); (e) jge exported and checked in 0.74 s.  It also found 172's
doubled `phase`.

**Built.**  `pycircuit/circuit/export.py`: `to_python(cir, name, check=True, allow_source=True)` -- the source of a
module whose `name()` rebuilds the flat circuit: the nodes in `cir.nodes` order first (gnd in its slot), the elements
in order, each its class (a collapse variant as its base; a BSourceHdl from `bsource_exprs` by `sympy.srepr`; a class
outside the library from its source) and only its GIVEN parameters, raw from `ipar`, a keyword through `**{}`, values
by `repr`; refused by name (`ExportError`): hierarchy, another toolkit, a callable or object parameter, a class local
to a function.  The check rebuilds each element beside the original and compares i/q/G/C at random x, u (dc, tran,
ac), `next_event` and CY bit for bit, the stateful limiters' analysis state cleared on both sides and put back --
state the parameters do not hold raises, never a silently different circuit.  `selfcheck(cir, src)` compares the
whole rebuilt circuit the same way (None, or the first difference); `rebuild(src, name)`; `x0_to_python(x)`.
Measured on the build: the gold decks rca, mosrect, gm1, opampal, jge exact with the DC bit for bit (jge's 568
elements: export with its check 0.64 s, the selfcheck 0.52 s, best of 3 under a suite's load -- the noise matrix
compared at one x, an HDL class's CY being Python's); the user session's 13 circuits exact, the buck's transient identical
(scratch only).  Tests: `test_export.py` -- an imported deck (level 3 MOS with `as`, a collapsed diode, a pulse, a
sine with its two phases) round-trips with its DC and transient bit for bit; a copy with every parameter given, and
one with its nodes in first-use order, caught; a BSourceHdl and a user class; the refusals; a replaced source
function; a diode's limiter state compared from none and put back; a sine's AC phase seen by the check; a resistance
an ulp off caught; a start vector exact; another toolkit refused.  `test_radau_fold`'s counted A/B runs now hold the
collector off (test_watch's convention): in this build's first full run one xdist worker freed a watched dict in the
middle of a run -- the watch counter moved and one fused hand-off was dropped (`radau.fuse:unhanded`, a count, never a
value), which the test's exact count comparison read as a difference.

**Gate**: the suite, 5 653 passed (6 skipped, 3 xfailed); the exporter's module doctest; the user session's 13
circuits exact on the build (scratch only).

## 2026-10-22 — 174. A waveform source's parameters written after construction reach its waveform

**Why** (Andreas: "Start with S3", the item-2 plan's correctness stage, found while mapping which parameters the
per-circuit caches hold).  Every source whose time function is made from its parameters -- `VSin`, `ISin`, `VPulse`,
`IPulse`, `VPWL`, `IPWL`, `VExp`, `IExp`, `VAM`, `IAM`, `VSFFM`, `ISFFM` -- built `self.function` once, in
`__init__`.  A frequency, a delay, a pulse width or a PWL point written afterwards reached `iparv` and nothing else:
`u`, `dudt` and `next_event` kept the first waveform, so a sweep of a sine's frequency returned the unswept answer
(checked: a 1 kHz `VSin` written to 2 kHz, `update_iparv` called, still drove 1 kHz).  A waveform parameter given as an
expression of a parent's parameter, unresolved at construction, kept its DEFAULT forever (a `freq='fs'` sine ran at
0 Hz).  Two more elements had the same shape: `SVCVS` built its stamps in `__init__` (a coefficient written later
never seen) and `VCVS_limited`'s `G` differentiated a limiter built in `__init__` while its `i` read `level` itself
(a `level` written later left the two limiting at different levels).

**Built.**  `_FunctionOfParameters` (elements): the twelve sources build `function` in `update()` -- which runs on
every parameter write (the `iparv` observer) and once at construction -- from a per-class `_time_function`; the IS
family, which had no `update`, gets one through it (its `iparv` is now observed).  `SVCVS`: the checks and stamps
moved to `_check_coefficients` / `_stamps`, run at construction and in `update()`; the state nodes stay its
topology, added at construction, and a write that changes the denominator's ORDER is refused (ValueError).
`VCVS_limited`: its limiter built in `update()`.  Tests (`test_elements`, each failing on the old code, 28 cases):
every waveform class written after construction (by attribute and by `set`) equals one built with the new values,
`u` / `dudt` / `next_event` to the bit; a sine's frequency from a parent's expression follows the parent's parameter;
a sine's frequency written into a circuit gives the transient (its C source table included) of a circuit built with
it, to the bit; `VCVS_limited`'s `G` and `i` after a `level` write; `SVCVS` after a coefficient write, and the order
change refused.

**Not changed, found on the way** (the next entry): a PCNR view's build re-evaluates every shared element's parameters
against the view's own (empty) parameters, so a parameter given as an expression of a parent's resets to its default
-- the circuit's own element, and the DC answer, with it.

**Gate**: the suite, 5 681 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone DC, against
3c7451a2: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 175. The PCNR view no longer re-resolves the circuit's own elements: a parameter from a parent's expression survives a DC

**Why** (found on the way to 174, reproduced before the fix).  `pcnr_view` builds a circuit that HOLDS the circuit's
own element objects, and it added each with `add_instance` -- whose last step resolves the instance's parameters
against its new parent's, here the view's, which has none.  A parameter given as an expression of a parent's (a
`SubCircuit` subclass's `vdd`, a source written `v='vdd'`) failed to resolve, under `ignore_errors`, to the parameter's
DEFAULT -- in the circuit's own element.  Any DC that built the view (a circuit with a PCNR participant: a diode, a
BJT, a MOSFET) returned the answer for the defaults (0 V at the source for 2 V) and left the circuit so; the second
solve repeated it.  A source inside a subcircuit instance whose own parameter is that expression went the same way
(the instance re-resolved, its children with it).  SPICE decks were not reached: the importer resolves `.param` to
numbers.

**Built.**  `SubCircuit._place_instance`: `add_instance` without the parameter resolution (`add_instance` is it
followed by `update_iparv`).  `pcnr_view` places the circuit's own elements and adds only its stubs.  Building the view
therefore no longer runs every shared element's `update()` (an R's or a VS's stamp rebuilt, an HDL element's packs
dropped) -- only the view's and its stubs' construction move the parameter epoch.  Tests (`test_pcnr_c`, failing on
the old code at 0 V): a top-level source and a source in a subcircuit instance, each from the parent's `vdd`, beside a
diode: the DC equals the numeric circuit's to the bit, and the element keeps its value.

**Gate**: the suite, 5 682 passed (6 skipped, 3 xfailed) and one failure, the pinned pair "the PCNR DC solve and its
C" (its reference holds `pcnr_view`; the twin, the C solve, runs on the same view, unchanged) -- re-recorded,
its file rerun (72 passed); the 27 gold decks, transient and standalone DC, against 174: every digest, step,
rejection, Newton count and path count identical.

## 2026-10-22 — 176. A source's value written between two solves keeps the circuit's set-up: live parameters

**Why** (Andreas: "HO can we speed up the Howland circuits even more", then the item-2 plan's S1).  A DC sweep point on
the Howland session's converters (n 85..164) cost 6.7..12.7 ms for ~1 ms of Newton: every per-circuit cache -- the
stamp plan, the PCNR view and its C record, the cores, the limiter walk -- is keyed on ONE global sentinel,
`ParameterDict._epoch`, replaced on any parameter write, so writing the DAC's `v` threw all of it away although none
of it holds a source's value (`u` is read on every pass, and the C records copy it in per solve).  One `ipar.v` write
moved the epoch twice; `update_iparv` on a two-element circuit five times (every element's `update()` runs twice).
pycircuit's own `DCSweep` writes and calls `update_iparv` at every point: 20..24 ms a point.  The ceiling, the epoch
put back after the write (scratch g19/ceiling.py): 2.4..3.5x a point, every answer bit-identical.

**Built.**
* `Circuit._live_params`: the parameters a class only reads at evaluation time -- VS `v`, `vac`, `phase`, `noisePSD`;
  IS `i`, `iac`, `phase` and its three noise parameters; every parameter of the twelve waveform sources (read, or made
  into their waveform by `update()`, 174).  Valid only while the class's `update` is the declaring class's own
  (`_live_names`): a subclass that overrides it may bake them and is not live.  Set on an element's `ipar` and
  `iparv` (`ParameterDict._live`) once it is built.
* `ParameterDict.notify(names=...)`: the epoch moves unless every written name is live; observers are notified as
  before.  `update_values` names only the values it CHANGED (`_unchanged`: the same object, or an equal number or
  string of the same type; NaN a change) and a given-flag set that changes counts as a change (a compact model
  branches on `$param_given`) -- so an `update_iparv` that changes nothing moves nothing.
* No-op rewrites removed from the elements' dicts, which the dict-watch stamps read (`_watch`): `update_values`
  replaces `_given` only when it differs, VS builds its topology-only `_G` once per toolkit, and a waveform source
  rebuilds `function` only when an argument or the toolkit changed.  Over a 25-point sweep the source pass's stamp now
  holds (2..3 arms, 50..66 checks served; every point re-armed before).
* `ParameterDict.copy` keeps the live names (a copied circuit's sources keep the set-up).
* `DC`'s docstring says what a solve keeps for the next; `_stamp_plan`'s rebuild rule names the exception.

Tests (`test_live_parameters`, 13; 8 failing on the old code): a source's value written between two DCs (VS, VSin,
VPulse; an IS) keeps the plan, the view and the C record and gives a fresh circuit's answer to the bit; a written
`v`, `vac`, `phase`, `noisePSD` reach the transient, AC and noise as a fresh circuit's; a pulse's width written after a
DC keeps the set-up and gives the fresh pulse's transient; a resistor's `r` and a participant's `IS` still rebuild; a
subclass overriding `update` is not live; an `update_iparv` that changes nothing keeps everything; a given flag
changed alone moves the epoch; a copied element keeps its live names; each class's live names are its own.

**Measured** (scratch g20/s1_bench.py, the Howland session's seven 25-point sweeps as they run them -- `ipar.v`
written, a new `DC` per point from the last answer -- and `DCSweep` on three; every answer bit-identical to 175):
  | converter | n | 175 ms/pt | S1 ms/pt | x |
  |---|---|---|---|---|
  | two_opamps | 153 | 10.54 | 3.82 | 2.76 |
  | csa_inverting | 118 | 10.58 | 4.56 | 2.32 |
  | low_side_sink | 85 | 5.72 | 1.62 | 3.54 |
  | csa_ref | 119 | 8.76 | 2.65 | 3.31 |
  | current_monitor | 90 | 6.14 | 1.66 | 3.71 |
  | csa_pfet | 111 | 7.53 | 2.10 | 3.58 |
  | parallel_howland | 164 | 10.83 | 3.30 | 3.28 |

  Median 3.31x, lowest 2.32x (the line: median >= 2x, none < 1.5x).  `DCSweep` (26 points, continuation): two_opamps
  21.81 -> 12.19 ms a point (1.79x), csa_ref 18.79 -> 11.10 (1.69x), parallel_howland 18.45 -> 10.97 (1.68x) (the
  line: >= 1.5x each).  Quiet box, AC power, warm best of 5 per process, 3 processes per tree interleaved, the
  process minima spread < 1 %.

**Gate**: the suite, 5 695 passed (6 skipped, 3 xfailed) before `ParameterDict.copy` kept the live names, 5 699 with it and the next
entry; the 27 gold decks, transient and standalone DC, against 175: every digest, step, rejection, Newton count and path
count identical.  The line's last clause -- no gold transient slower than its lottery band + 3 % -- is
NOT YET MEASURED (the box was not quiet); it is measured after this commit.

## 2026-10-22 — 177. An element's own parameter write resolves against its parent: expressions of a parent's parameter survive it

**Why** (found while planning item 2, confirmed by running it; Andreas: "Build A", of three options).  An element's
own write (`elem.ipar.x = ...`) re-resolved ALL its parameters through `_ipar_changed` -> `update_iparv()` with no
parent, so every other parameter given as an expression of a parent's failed to resolve -- under `ignore_errors` --
and fell back to its default: `R(r='10+x')` in a circuit with x = 20 read 30 ohm, and 1 kohm after its `noisy` was
written.  A new expression written onto a leaf never resolved at all.  The two other options -- keep the last value
where an expression cannot resolve; re-resolve only the written names -- leave a new parent expression written onto
a leaf silently stale.

**Built.**  `Circuit.update_iparv` REMEMBERS the context it resolved in (`_resolved_in`: the parent's resolved
parameters and the globals; a missing part keeps the one remembered, so the parent's observer path, which passes no
globals, does not drop them), stored only when it changes -- the dict-watch stamps read the element's dict.
`_ipar_changed` resolves against it, so a nested instance's own write resolves against its parent and re-resolves
its children.  `Circuit.__copy__` keeps the context; `SubCircuit.__copy__` points a copied element whose context was
the original circuit at the copy.  The PCNR view places elements without resolving them (175), so it never becomes
their context.  Tests (`test_circuit`, 3, failing on the old code): a leaf's own write keeps another parameter's
parent expression, and a new parent expression written onto it resolves, a nested instance's too; a copied circuit's
element resolves against the copy; a re-resolution in the same parent leaves the context object as it was.

**Found, not changed**: copying a circuit that holds a subcircuit instance raises `KeyError` in `update_node_map`
(before and after this change).

**Gate**: the suite, 5 699 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone DC, against 176: every
digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 178. A parameter written by attribute is given and range-checked, as one written by set -- the MOS levels' `gamma` and `phi` written after construction are the device's (an attribute write was ignored: they stayed derived from nsub and tox)

**Why**  `pd.x = v` stored the value and nothing else: the name was not GIVEN and its declared range was not checked,
both of which `set` does.  A model that branches on $param_given -- SPICE's MOS levels 1-3 take `gamma` and `phi` as
given or derive them from `nsub` and `tox` -- ignored a value written the usual way, `m.ipar.gamma = 0.5`: 21
class/parameter pairs in the library (a sweep of every element class's parameters written after construction against a
fresh instance built with them).

**Built.**  param: `__setattr__` checks the range and gives the name; `eval_expressions` writes the resolved values
without giving them and carries what the source was given; `forget(name)` un-gives a name; a change of given flags
that only live names make keeps the epoch. dcanalysis: DCSweep restores whether the swept parameter was given.

**Tests.**  tests (6, failing on the old code): test_param -- an attribute write gives and range-checks, forget
un-gives, resolution carries what was given; test_elements_hdl -- MOS level 1 and 3 `gamma`, level 1 gate charge
`phi`, written by attribute, give the stamps of a device built with them, to the bit; test_dcsweep -- a swept
parameter that was not given is left not given.

**Gate**: the suite, 5 705 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone DC, against 177:
every digest, step, rejection, Newton count and path count identical; the library sweep of parameters written after
construction: 34 pairs differing before, 13 after (the 21 given-or-derived MOS ones gone).

## 2026-10-22 — 179. AC and noise without a given operating point solve it at their own temperature (at the default 300.15 K until now: a diode stage at 400 K read a gain of 0.95 for 0.12)

**Why**  `dc_steady_state` solved its operating point with `DC(cir)` -- the default temperature, whatever the
analysis's `epar` said -- and then linearised at the analysis's temperature around it.  FeedbackLoopAnalysis and
TwoPortAnalysis, which run AC, inherited it.

**Built.**  analysis_ss: `DC(cir, epar=epar)`.  tests: test_analysis_numeric -- AC and noise at 400 K without dcx
equal the ones given the 400 K operating point, to the bit (failing on the old code).

**Gate**: the suite, run once over this commit and the three after it (179-182): 5 713 passed (6 skipped, 3 xfailed)
with one failure -- a test that reads the source of OpAmp for its limiter's per-call JAX import, which 180 moved from
__init__ to _functions; pointed there, it passes; the 27 gold decks, transient and standalone DC, over the four
against 178: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 180. MOS's and OpAmp's parameters written after construction reach their child elements (the children kept the values of the construction)

**Why**  MOS and OpAmp are built from child elements, each given a parameter's value at construction; OpAmp's two
behavioural sources captured theirs in closures.  A later write of `gm`, `gmb`, `gds`, `Cgs`, `Cgd`, `Cdb` (MOS) or
`Aol`, `GBW`, `Rin`, `Rout`, `Vdd`, `Vss` (OpAmp) changed nothing (12 of the library sweep's pairs).

**Built.**  mos: the children take expressions of MOS's parameters (resolved against it, so a write reaches them); the
noise term rebuilt in `update`.  macromodels: the resistors and the pole capacitor take expressions (the capacitor's
in the order its arithmetic always ran, so its value is the same to the bit); the sources' functions are rebuilt by
`update` only when their parameters changed, so a re-resolution that changes nothing keeps the circuit's set-up.
Stamps at construction unchanged, to the bit.

**Tests.**  tests: test_elements -- each parameter of both written after construction gives the stamps of one built
with it, to the bit; the OpAmp's sources kept through a no-op re-resolution and replaced by a write of `Aol` (failing
on the old code).

**Gate**: the suite, run once over this commit and the three after it (179-182): 5 713 passed (6 skipped, 3 xfailed)
with one failure -- a test that reads the source of OpAmp for its limiter's per-call JAX import, which 180 moved from
__init__ to _functions; pointed there, it passes; the 27 gold decks, transient and standalone DC, over the four
against 178: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 181. IdtmodCircular refuses a `gamma_tau` write that would add or remove its filter state (it was silently ignored)

**Why**  `gamma_tau` > 0 adds the circular integrator's filter state -- a node -- at construction, by design; a write
after construction from 0 to a positive value, or back, was silently ignored.  elements: `update` raises ValueError
when the written value would add or remove the state; one positive value over another is a value, as before.  tests:
test_idtmod (failing on the old code).

**Gate**: the suite, run once over this commit and the three after it (179-182): 5 713 passed (6 skipped, 3 xfailed)
with one failure -- a test that reads the source of OpAmp for its limiter's per-call JAX import, which 180 moved from
__init__ to _functions; pointed there, it passes; the 27 gold decks, transient and standalone DC, over the four
against 178: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 182. A subcircuit instance's later structure reaches the circuits holding it, and a deletion rebuilds the branch rows after it (an instance given its elements after it was added raised KeyError in every analysis and copy)

**Why**  A subcircuit instance added to its parent and given its elements afterwards left the parent without them: its
nodes and branches stayed as at the addition, and DC (or a copy) raised KeyError on the instance's internal node.
Deleting an instance removed its branches by equality -- ambiguous for parallel elements -- and left every later
instance's branch span unshifted, the deleted one's kept: after `del c['V0']` an inductor's row pointed one past the
end.

**Built.**  circuit: an instance added with `add_instance` remembers its parents (weakly); a change to its own
structure -- an element or node added, an element deleted -- refreshes them (`_refresh_instance`: its internal nodes
follow, the branches are rebuilt), and theirs in turn.  The branch list and the spans are rebuilt from the instances
in order (`_rebuild_branches`), on a deletion too.  A copied circuit's subcircuit instances remember the copy.
Circuits built elements first are unchanged.

**Tests.**  tests: test_circuit -- an instance given its elements after it was added: the elements-first circuit's
nodes, branches and DC to the bit, and its copy; two levels, an element added to the inner instance reaches the top
and its deletion leaves it; after a deletion each later instance's branch rows are its own (failing on the old code).

**Gate**: the suite, run once over this commit and the three after it (179-182): 5 713 passed (6 skipped, 3 xfailed)
with one failure -- a test that reads the source of OpAmp for its limiter's per-call JAX import, which 180 moved from
__init__ to _functions; pointed there, it passes; the 27 gold decks, transient and standalone DC, over the four
against 178: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 184. Copying a SubCircuit no longer updates the fresh construction's elements before replacing them -- a symbolic circuit holding a MOS copies again (raised TypeError since MOS's children took expressions, 03a5413f)

**Why**  `Circuit.__copy__` builds a fresh instance, copies the parameters into it and calls its `update`;
`SubCircuit.__copy__` then replaces the fresh instance's elements with copies of the original's.  The update ran on
the fresh construction's elements -- on the default toolkit, not the copied circuit's -- and once MOS's children took
expressions of its parameters (03a5413f) it resolved a symbolic `gm` into a numeric VCCS: `Cannot convert expression
to float`, so TwoPortAnalysis (which copies its circuit) on the regulated-cascode example failed.  The gates of
03a5413f missed it: they ran in a worktree, and the example runs in a fresh interpreter that imports the main tree's
installed package.

**Built.**  circuit: `Circuit.__copy__` skips the update for a SubCircuit, whose elements are replaced next (the final
copy is what it was).

**Tests.**  tests: test_elements -- a symbolic circuit holding a MOS, copied with the default toolkit numeric: the
copy's children are the original's (failing on the old code with that TypeError).

**Gate**: the suite IN THE MAIN TREE (the examples' interpreter imports it), 5 715 passed (6 skipped, 3 xfailed), the
regulated-cascode example among them; the 27 gold decks, transient and standalone DC, against 182: every digest, step,
rejection, Newton count and path count identical.

## 2026-10-22 — 176b. Live parameters (176): the line's last clause measured -- no gold transient slower

176 (84dc36d7) went in with its line's last clause -- no gold transient slower than its lottery band + 3 % -- NOT YET
MEASURED (the box was busy).  Measured on a quiet box (AC power, each deck run only when the 1-minute load was under
1.5): the 27 gold transients by `run.py` (one cold solve a process), 9baca22f against 84dc36d7 interleaved, twice
each: every digest identical; summed 35.927 s against 35.932 s.  Seven decks read over 3 % slower there -- mosrect
1.20x (a load spike to 4.5 during its runs), schmitecl 1.14x, latch 1.10x, mike2 1.08x, gm2, gm6, toronto 1.03..1.04x,
all of 30..250 ms, where one cold solve is known to read 9..12 % as noise -- so those seven were re-timed WARM
(best of 10 in a process, three processes per tree, the base run twice as an A/A, rotated): 84dc36d7 / base 0.993 ..
0.999 on every one of them, the A/A 0.993 .. 1.005, digests identical.  The clause is met: 176's line is met in full.

## 2026-10-22 — 183. Under the JAX toolkit a parameter written after the evaluation groups were built is evaluated, and a JAX transient solves its operating point at its own temperature (the write was ignored; a run at 400 K started from the 300.15 K point)

**Why**  Two defects of the kinds the previous commits fixed, on the JAX path. The JAX toolkit evaluates a circuit's
elements in groups whose parameters are copied when the groups are built -- by `update_iparv` alone: a parameter
written afterwards was not evaluated (a diode's `IS` written, the circuit's `i` and `G` unchanged).  And the JAX
transient, given no start, solved its operating point with `DC(self.cir)`, the default temperature: at 400 K it
started 0.17 V off.

**Built.**  circuit: `SubCircuit._current_groups` rebuilds the groups where the parameter epoch moved since they were
built -- except inside a JAX transformation (an argument a tracer), where arrays built would leak out of it: the
groups as built serve, and the next evaluation outside rebuilds them (a JAX transient evaluates the circuit outside
its compiled step first).  jaxtransient: `DC(self.cir, epar=self.epar)`.

**Tests.**  tests: test_jaxtransient, 3, failing on the old code -- a parameter written after the groups were built
gives a fresh circuit's `i` and `G`, to the bit; a JAX transient given its start after such a write runs as the fresh
circuit; a JAX transient at 400 K starts from the 400 K operating point.

**Gate**: the suite in the main tree, 5 718 passed (6 skipped, 3 xfailed) and 1 error: the leak detector at
test_hdl_diagnostics' module teardown (DiodeHdl's 'auto' class resolved while a test hid the compiler) -- NOT this
change: the module run alone gives the same 2 errors in the main tree at 9baca22f, before today's commits, and none in
a clean worktree at any of today's commits (an environment- and order-dependent leak of the main tree, open); pinned
pairs 'the assembly loops and the plan' and 'the sources and their direct pass' re-recorded (their references, the
assembly loops, now fetch the groups through _current_groups; the twins -- the stamp plan, the direct source pass --
serve the numeric toolkit, where the groups are empty: unchanged); the 27 gold decks, transient and standalone DC,
against 184: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 185. An 'auto' class's backend resolved through a subclass is noted on the class that owns it -- the leak detector's intermittent report at test_hdl_diagnostics' teardown, fixed at its cause

**Why**  A plain Python subclass of a Behavioural class (one that overrides a method, as test_hdl_diagnostics'
HalfG(DiodeHdl) does) shares its base's `_hdl_info`.  Built before the base was, its first instance resolved the
SHARED info's 'auto' backend -- C bound, pending cleared -- but `_resolve` noted the status, the ULP list and the
limiter on the subclass: the base's own status still read 'auto (resolved at the first instance)'.  The leak detector
compares the base's state before and after each test and reported "an 'auto' class resolved ... (not C: a compiler,
cffi or the compile cache was hidden from it)" -- whenever no earlier test in that worker had built the base itself,
so in the main tree only, and in the full suite only for some orders (it surfaced on the previous commit's gate; the
module run alone reproduced it in the main tree at 9baca22f, before today's commits).

**Built.**  _hdl_cbackend: `ensure` resolves against `_owner(cls, info)`, the class in the MRO whose `_hdl_info` the
info is.

**Tests.**  tests: test_hdl_cbackend -- a fresh 'auto' class, a subclass of it built first: the base's info resolved
and the base's status 'c', the subclass reading it through (failing on the old code: the base's status still the
pending note).  test_hdl_diagnostics run alone in the main tree: 34 passed, the two leak errors gone.

**Gate**: the suite in the main tree, 5 719 passed (6 skipped, 3 xfailed), no failure, no error, no leak (a first run
failed only test_runtime_map's no-history-in-C-modules rule on this change's own docstring, fixed);
test_hdl_diagnostics alone in the main tree 34 passed; the pinned pairs 72 passed; the 27 gold decks, transient and
standalone DC, against 183: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 186. The PCNR device list kept on the circuit and the seed's limiting through each class's C refine kernel -- the Howland session's seven DAC sweeps 1.28..2.08x a point again (median 1.68x), every answer bit-identical

**Why**  Item 2's S2a.  After the set-up was kept across source writes, a sweep point still redid two things on every
solve: `pcnr_devices` built new device records up to three times a solve (13..15 % of a point, profiled), and the
seed's limiting -- `v_lim_init` -> `refine` -- ran each participant's generated Python law (37..48 %), though each
class has a C refine kernel, the one PCNR's C solve runs and pinned to `PcnrDevice.limit` bit for bit.

**Built.**  pcnr: `pcnr_devices` keeps its records on the circuit under the view's key (the parameter epoch, the
elements dict, the node map, the size), a new list each call; `refine` runs `_pcnr_c.refine` per device where the
class has a kernel and the C is on, the device's Python law where it declines or the C is off.  _pcnr_c: `_packed`
declines a class whose `i` has no C binding (it raised AttributeError, asked for a refine alone). Pinned pair "the
PCNR DC solve and its C" re-recorded (its reference holds `refine`; the twin's digest unchanged).

**Tests.**  Measured (scratch bench, the sweeps as the Howland session runs them -- `ipar.v` written, a new DC per
point -- quiet box, warm best of 5, two processes per tree interleaved), b06c28fe -> this, ms a point: two_opamps 3.86
-> 2.33, csa_inverting 4.56 -> 3.57, low_side_sink 1.63 -> 0.78, csa_ref 2.66 -> 1.71, current_monitor 1.67 -> 0.83,
csa_pfet 2.11 -> 1.26, parallel_howland 3.31 -> 1.80 (median 1.68x, the line >= 1.25x); DCSweep 1.12..1.46x; every
answer bit-identical.

tests: test_pcnr_c, 2, failing on the old code -- the device records kept through a source write and renewed by a
participant's parameter or an element added; the seed's refine with the C on equals it with the C off, to the bit, and
the kernel served.

**Gate**: the suite in the main tree, 5 719 passed (6 skipped, 3 xfailed), no failure, no error, no leak (before the
two tests were added; test_pcnr_c with them, 48 passed); the pinned pairs 72 passed; the 27 gold decks, transient and
standalone DC, against 185: every digest, step, rejection, Newton count and path count identical.

## 2026-10-22 — 187. A DC's default linear solver kept on the circuit -- one DC per sweep point keeps the solver's choice, its KLU analysis and its factorisation, as one DC object across the points always did (answers move at the rounding level where KLU reuses its pivots)

**Why**  Item 2's S2b (Andreas: "Build s2a and s2b", the bit moves accepted). The usual way to write a DC sweep makes
a new DC per point, and each new DC made a new AutoSolver: on a circuit of 150 unknowns or more it chose StructuralKLU
again, built the circuit's structure and ran KLU's symbolic analysis, every point.  One DC object kept across the
points (and DCSweep) never paid that.

**Built.**  dcanalysis: `DC._get_linearsolver` takes, where no solver is given, the circuit's kept AutoSolver
(`_circuit_solver`: per element node map object, size, reference node and structural floor; its structural pattern
read from the circuit, held weakly).  An explicit solver is the DC's own.  Where the solver is KLU a point's
factorisation reuses the last point's pivots: against a fresh solver the answers move at the rounding level (the
Howland sweeps: two_opamps 2.2e-16 V at most, parallel_howland 8.8e-14 V; the five dense ones bit-identical).

**Tests.**  Measured (scratch bench, the sweeps as the Howland session runs them, the box idle, two processes per tree
interleaved), 1ea29982 -> this, ms a point: two_opamps 2.37 -> 1.29 (1.84x) and parallel_howland 1.80 -> 0.70 (2.57x)
-- the two that factor with KLU; the five dense ones and DCSweep unchanged (0.99..1.02x).

tests: test_dcsweep -- a new DC per point keeps the circuit's solver (KLU on a 160-section ladder), an explicit one is
the DC's, a topology change brings a new one (failing on the old code); a point on the kept solver equals a fresh
solver's to 1e-12.

**Gate**: the suite in the main tree, 5 723 passed (6 skipped, 3 xfailed), no failure, no error, no leak
(test_dcsweep's helper then given f-strings for the lint; its file rerun, 11 passed); the pinned pairs 72 passed; the
27 gold decks, transient and standalone DC, against 186: every digest, step, rejection, Newton count and path count
identical (one DC a process: nothing kept to reuse).

## 2026-10-22 — 188. A hierarchical circuit is evaluated through its flat view -- the C paths serve circuits of subcircuit instances (they declined them: every element ran the Python paths, 7..40x slower; a user's Howland testbench, its first DC 40x)

**Why**  The C paths -- the evaluation core, the C Newton and stepping loop, PCNR's C solve -- take a flat list of
leaf elements, each a constant stamp or a compiled HDL kernel.  A subcircuit instance is neither: a container its own
Python assembly sums.  The core declined a circuit that held one ('unservable') and every element of it ran the Python
paths; PCNR, which collects its participants from the top level, did not see the devices inside the instances at all.
The SPICE importer flattens .subckt instances, so no gold deck ever met it.  Measured: a Howland current pump's
testbench, hierarchical (its op-amp, transistors and optocoupler as subcircuit instances), took 1.667 s for its first
DC against 0.042 s flattened; with only its op-amp wrapped in one instance a DC from zeros took 231 ms against 21 ms,
a 20 us transient 205 ms against 29 ms.

**Built.**  circuit: `SubCircuit._flat_view` -- the circuit itself where it holds no subcircuit instance, else its
leaves (the same objects) placed on one level, each on its own nodes as seen from the top (an instance's terminals its
connections, its internal nodes `inst.node`), the circuit's nodes and branches in its own order, so a solution vector
means the same in both; kept per topology; an instance whose class overrides an evaluation method stays whole; a view
whose nodes or branches would differ is not used (counted).  The leaves are placed, not added, so their parameters
stay resolved in their own parents and their toolkit as it was (a circuit may mix them); the view's `update_iparv` is
the circuit's.  An element whose own name has dots (a view's leaf) is found by that name first.  analysis: an analysis
keeps the circuit given (`circuit`) and evaluates its view (`cir`), re-synced at a DC's or transient's solve where the
topology changed since (`_sync_cir`); `CircuitResult` reports against the circuit.  dcanalysis: DCSweep writes the
swept parameter into the circuit given (a subcircuit instance's own parameter, which the view's leaves resolve).

**Tests.**  Measured, the IMD testbench with its op-amp in one instance: a DC from zeros 231 -> 17 ms, a 20 us
transient 205 -> 15 ms, the answers within 1e-10 V; a flat circuit gets no view and runs as before.

tests: test_initial_conditions -- an ic inside a subcircuit, refused until now (a nested instance owned one span over
its children's branches), is placed on its inductor's own branch row: the test that pinned the refusal now pins the
placement.  test_flat_view, 7 (4 failing on the old code) -- a hierarchical circuit served by the PCNR C solve and the
C Newton, its answers the hierarchical evaluation's; the view's unknowns the circuit's and its results reported
against it; a subcircuit's parameter swept by DCSweep; a topology change after the analysis was made in its solve; a
leaf's own toolkit kept; a flat circuit without a view and a nested one complete.

**Gate**: the suite in the main tree, 5 729 passed (6 skipped, 3 xfailed) and one failure: test_initial_conditions'
refusal of an ic inside a subcircuit, which the view now places -- the test rewritten to pin the placement, its file
rerun, 24 passed; the pinned pairs and the runtime map 82 passed; the 27 gold decks (flat: no view), transient and
standalone DC, against 187: every digest, step, rejection, Newton count and path count identical.  Timing, the box
idle, three processes per tree interleaved: the seven flat Howland sweeps 0.984..1.012x against 187 (median 0.994; the
view's own cost measured 0.4 us a DC point, 0.05 %, the rest within the trees' +-1 % spread), bits identical; the
hierarchical IMD testbench 13..14x (DC from zeros 231 -> 17 ms, a 20 us transient 205 -> 15 ms).

## 2026-10-22 — 189. The dense multistep walk carries its monodromy on the live columns -- the columns a charge reaches, after the opening (a Howland testbench's gear PSS 1.48x, 11.84 -> 8.01 s, the orbit bit-identical)

**Why**  The shooting Newton's monodromy rides along every step of a multistep walk (`_lmm_recursion`): ``S = sum a_k
C_{n-k} P_{n-k}``, ``P_n = -Jf^-1 S``, ``Pq_n = a_0 C_n P_n + S`` -- every column of the unknown, a dense product and
a dense solve on each.  Every step reads the state through a capacitance, so after the opening a column survives only
where a charge reaches it; the rest are exactly zero for good (the recursion is linear in the state and the
monodromy's carries no forcing).  Measured on a user's Howland current pump testbench (90 unknowns, its two-tone gear
PSS: 5 walks of 3 200 steps): the recursion was 6.9 s of 11.8, on gear's pair 180 columns of which 42 are ever nonzero
-- twice the 21 columns of its 90 x 90 capacitance.

**Built.**  _steps: `_LiveMonodromy` carries ``(P_n, P_{n-1}, Pq_n)`` and drops a column once it is zero in all three;
`full`/`ring` put the zeros back. Each kept column is computed as it was.  A step whose result is not finite runs
again at full width from the full state, and the walk keeps every column from there on: a dropped column multiplied
into a non-finite `C` or factor is NaN at full width (``0 * inf``), not zero (`pss.live:nonfinite`).  `LIVE_COLUMNS`
(True) switches it.  _pss_walks: `_walk_lmm` carries its dense monodromy through it; a capture and the walk's `P` read
it at full width.

**Tests.**  Measured, the box idle, three processes per tree interleaved: the testbench's gear PSS 11.84 -> 8.01 s
(1.48x), the orbit's digest identical in all six.

tests: test_shooting_live_columns, 9 (all failing on the old code: no live columns) -- the walk's monodromy on the
live columns equal to every column's (an RLC and a diode clipper, trap, euler and gear's pair; zero pattern identical,
values to 1e-13; the steps past the opening solve the columns a charge reaches); a PSS answering as on every column
(trap, gear); a non-finite step run at full width, its NaN where every column puts it.  test_shooting_pss: the
explicit-inverse test's width assertion -- one matrix solve a step, the opening's at full width, the rest on the live
columns.

**Gate**: the suite in the main tree, 5 739 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 188: every digest identical; the shooting files 645 + 9 passed before it (the explicit-inverse
test's width assertion updated).  Timing, the box idle (1-min load 1.5..1.8), three processes per tree interleaved:
the Howland testbench's gear PSS 11.84 -> 8.01 s (1.48x; the line 1.3x), the orbit's digest identical in all six.

## 2026-10-22 — 190. The stage walk carries its monodromy on the live columns too, and the multistep walk's products read the capacitance's nonzero columns -- a Howland testbench's Radau PSS 1.14x (62.4 -> 54.7 s), its gear PSS a further 1.04x (7.94 -> 7.61 s), the orbits bit-identical

**Why**  The stage family's dense walk (Radau IIA, TR-BDF2, ESDIRK; `_walk_stage`) propagates its `m x m` monodromy
through every step as the multistep walk did: ``C_n P`` stacked per stage and solved against the stage block (the C
step map, `_sens_c`, or the Python `_StageStep.solve`).  After the first step a column survives only where a charge
reaches it -- 21 of the testbench's 90 -- and the block's solve ran on all 90.

**Built.**  _steps: `_LiveColumns` carries a stage walk's `P` on its live columns, as `_LiveMonodromy` (189) carries
the multistep walk's; a step whose result is not finite runs again at full width through the Python step, and the walk
keeps every column from there on.  _pss_walks: `_walk_stage` carries its dense monodromy through it -- the C step map
and the Python step alike; a capture and the walk's `P` read it at full width.

**Tests.**  _steps: `_LiveMonodromy` reads each capacitance through its nonzero columns (`_ColumnsOf`: ``C[:, k] @
P[k]``; a NaN or inf is nonzero, so it is read) -- the terms dropped are exact zeros, so a product is the dense one's
to the rounding, its sums grouped otherwise -- and a step whose result is not finite runs at full width on the dense
products. It looks for newly dead columns on the first four steps and every 64th after, not on every step: a column
dies with the opening's identity leaving the ring, or by a coincidence.

Measured, the box idle, three processes per tree interleaved, against 189: the testbench's Radau PSS 62.42 -> 54.68 s
(1.14x), its gear PSS 7.94 -> 7.61 s (1.04x); the orbit's digest identical in all twelve. OpenBLAS's threads do not
enter: the gear PSS 7.62 s with the default, 7.57 s on one.

tests: test_shooting_live_columns, 7 more (failing on the old code: no live columns in a stage walk) -- a stage walk's
monodromy on the live columns equal to every column's (an RLC and a diode clipper; Radau through the C step map on the
RLC, the Python step where a junction declines it; TR-BDF2, ESDIRK43); a Radau PSS answering as on every column.  The
non-finite test compares the finite entries to the rounding.

**Gate**: the suite in the main tree, 5 746 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 189: every digest identical; the pinned pairs 72 passed (the C step map's twin unchanged: it is
called on the live columns).  Timing, the box idle (1-min load < 1.6), three processes per tree interleaved, against
189: the Howland testbench's Radau PSS 62.42 -> 54.68 s (1.14x), its gear PSS 7.94 -> 7.61 s (1.04x); the orbit's
digest identical in all twelve.

## 2026-10-22 — 191. A failed rescue is a TransientStepError whichever way the ladder's last rung failed -- one ending on a singular Jacobian raised a bare SingularMatrix, past a caller's `except NoConvergenceError` (a user's step Monte Carlo)

**Why**  `Transient._rescue_step` re-solves a step stuck at `minstep` with the continuation chain armed (gmin, gshunt,
pseudo-transient) and, where that fails, raises `TransientStepError` -- a `NoConvergenceError` -- naming the point.
It caught only `NoConvergenceError` from the chain: a ladder whose last rung ends on a singular Jacobian raises
`SingularMatrix`, which is not one, and escaped bare.  Found by a user's Howland current pump testbench, whose step
Monte Carlo retries a failed transient at a looser tolerance on `NoConvergenceError`: on their machine the
pseudo-transient rung ended "Singular Jacobian" and the retry never ran (on ours the same step's ladder ends
non-converged, and the error was the TransientStepError).  The step itself failed on its tolerance -- reltol 1e-9 with
iabstol 1e-15, below the Newton update's noise on near-zero source currents at 0.04..1 ohm nodes -- not on a path.

**Built.**  transient: `_rescue_step` catches `SingularMatrix` with `NoConvergenceError`; the error is the
`TransientStepError` either way, the singular matrix its cause.

**Tests.**  tests: test_analysis_transient -- a rescue ending singular is the step's TransientStepError, a
NoConvergenceError with the SingularMatrix as its cause, whether or not the path carries a ladder (failing on the old
code: the SingularMatrix escaped).

**Gate**: the suite in the main tree, 5 747 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone
DC, against 190: every digest identical; the new test failing on the old code, passing on the new.

## 2026-10-22 — 192. A fully implicit stage method's coupled stage block factored with KLU on its structural pattern -- a Howland testbench's Radau PSS 1.28x (54.88 -> 42.86 s), the orbit bit-identical

**Why**  The stage walk factors every step's coupled stage block (Radau's `3 m`: ``J[i][j] = delta_ij C(Y_i) + h A_ij
G(Y_j)``) through the analysis's solver, and that solver decides dense or sparse once, on the first matrix it sees: on
the testbench (90 unknowns) a 90-unknown Jacobian, so every 270-unknown block was a dense LU -- 540 us a step, of a
block with 2 580 nonzeros.  On a circuit whose first matrix was a block, the same choice is SuperLU on the values'
pattern.

**Built.**  _pss_inner: `_InnerTransient._step_klu(n)` -- for a coupled stage block (``n = s m``, ``s >= 2``) from
`SPARSE_BLOCK_MIN_N` unknowns (`linearsolver.MIN_N_FOR_SPARSE`), on the numeric toolkit and where no `linearsolver=`
was given, a `StructuralKLU` on the block's structural pattern (the circuit's, `analysis.jacobian_structure`, in every
block), one per size and topology; `_factorise` factors such a block there, each factor a numeric of its own (a
function of its block alone), a matrix off the pattern going on to the caller's solver.  Measured on the testbench's
block: KLU's factor and 21-column solve 227 + 65 us against dense LAPACK's 540 + ~150.  `PYCIRCUIT_SHOOT_KLU=0`
(`SHOOT_KLU`) switches it.  A step's own Jacobian is not routed: on the testbench's gear PSS (90 unknowns) KLU's
per-call cost lost 3 % (7.57 -> 7.78 s).

**Tests.**  _sens_c: the C step map declines a block the Python step factors on its structure (`pss.sens:sparse`,
before its solver check) -- its dense LU is not that step's twin; `_step_klu` joins the walk names an instance shadow
of declines.

tests: test_shooting_step_klu, 3 (failing on the old code: no `_step_klu`) -- a large block factored on its structure,
the C step map declining it, the orbit the old solver's to 1e-9; a step's Jacobian, a small block and a caller's
solver left alone; a factored period's block factors and their transposed solves (PAC's adjoint row and sideband
response) answering as without.  test_pinned_pairs: 'the stage step map and its C' re-recorded -- the Python step's
`_factorise` routes a large block to KLU, the C declines exactly those (its twin unchanged where it serves: the
caller's dense solver); `_step_klu` added to its reference.

**Gate**: the suite in the main tree, 5 750 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 191: every digest identical; the shooting files, pinned pairs and runtime map 674 passed before
it (the stage step map's pair re-recorded).  Timing, the box idle (1-min load < 1.4), three processes per tree
interleaved, against 191: the Howland testbench's Radau PSS 54.88 -> 42.86 s (1.28x), its gear PSS 7.61 -> 7.57 s
(unchanged); the orbit's digest identical in all twelve.  KLU on the step Jacobian, measured and refused: gear 7.57 ->
7.78 s (pinned), 7.70 -> 7.84 (unpinned).

## 2026-10-22 — 193. The dense multistep walk's monodromy step in C -- bit for bit the Python step, a Howland testbench's gear PSS 1.11x (7.57 -> 6.83 s; built under a 1.12x line it missed by a hair, committed on Andreas's decision)

**Why**  The multistep walk's monodromy step on its live columns (`_LiveMonodromy.step`: `_lmm_recursion`, each
capacitance read through its nonzero columns by `_ColumnsOf`, the step's solve through the caller's dense solver --
numpy's `linalg.solve`) took 165 us a step in place on the testbench's gear PSS (90 unknowns, 38..42 live columns),
around ~80 us of numpy's own products and solve.

**Built.**  csrc/lmmsens.c, lmmsens.h: the step in one call, `pcx_lmmsens` -- numpy's own `cblas_dgemm` (row-major,
neither transposed) on each capacitance's nonzero columns and the monodromy's rows they read, gathered C-contiguous as
numpy gathers them (numpy's matmul of those arrays is that call: checked equal on 336 shapes, threaded and single),
numpy's own loop for a one-column factor (``0 + a b``) and zeros for none (numpy's special cases); numpy's own `dgesv`
on the column-major copies `linalg.solve` makes; every elementwise operation in `_lmm_recursion`'s order.  It hands
back where numpy's solve raises (a singular `Jf`) and where the result is not finite (the Python step's full-width
fallback).  The capacitances' column lists are kept from call to call: the ring moves each array on one place a step.
`dgemm64_t` is declared in its own header, the one unit that multiplies matrices.

**Tests.**  shooting/_lmmsens_c: the glue -- the switch `PYCIRCUIT_LMMSENS_C` (`LMMSENS_C`), numpy's OpenBLAS by
address, the pieces the C stands in for checked genuine (the step, the recursion, `_ColumnsOf`, the solvers' `solve`
and choice, the toolkit's `linearsolver`, numpy's `linalg.solve`: an array-function dispatcher, told by its type), the
toolkit and solver checked again where they change, the arrays' kind and layout (an array the last call made or read
known); `pss.lmm:` counts each decline and hand-back and `pss.lmm:served`.  _steps: `_LiveMonodromy.advance` -- the
state one step on, the step's own bookkeeping, shared by the Python step and the C's.  _pss_walks: a dense walk
keeping no step asks the C first.

Measured, the box idle, three rounds interleaved in one process: the testbench's gear PSS 7.57..7.63 -> 6.83 s (1.108x
best, 1.116x median), every walk's monodromy and the orbit bit-identical in every round, 15 974 of 16 000 steps served
(25 at each walk's opening: the Python step's results, numpy's solve's, are F-ordered until the first prune copies
them; one before the solver chose).  The C 108 us a step in place against ~74 us of LAPACK and BLAS; the glue ~20.

tests: test_lmmsens_c, 34 -- the C step against the Python recursion bit for bit over every product branch (no, one,
several and all columns), with and without the companion's `b`, two- and three-term companions; kept column lists and
pointers answering as fresh ones; a singular Jacobian and a non-finite result handed back; a walk's monodromy the
Python step's (RLC, a diode clipper, a one-capacitor clipper; gear's pair, trap, euler), a PSS's orbit; the declines
(every column carried, a caller's solver, a patched piece) and a kept walk left alone. test_pinned_pairs: 'the
multistep step map and its C' recorded. test_c_abi, test_paths: the unit and its entry listed.  _testing/leaks:
`LMMSENS_C`, `MAXM` -- and `SHOOT_KLU` (192), `LIVE_COLUMNS` (189), the switches those stages added unlisted.

**Gate**: the suite in the main tree, 5 792 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 192: every digest identical; the sanitized fast tier (scripts/sanitize_suite.sh:
AddressSanitizer and UndefinedBehaviorSanitizer on every C object) 5 511 passed, no sanitizer report, and the two
allocator-pin tests failing as they fail on 192 under the sanitizer (ASan's allocator refuses glibc's mallopt; the
script does not deselect them); the pinned pairs, ABI, paths and runtime-map checks 231 passed.

## 2026-10-23 — 194. A shooting Newton's full step that settles is not halved -- a Howland testbench's Radau PSS 2.10x (42.27 -> 20.14 s), its gear PSS 1.52x (6.83 -> 4.48 s), the orbit within 6.3e-14 V

**Why**  Profiled (Andreas: "do it"), the testbench's default Radau PSS spent 8 traversals of the period where 4
would do.  The first shooting step put the residual at the arithmetic's floor (1e-14); the next full step (6e-15)
did not lower that rounding, so `fsolve`'s line search halved it four times -- four more traversals with
sensitivities, 5.5 s each -- and the convergence test, run after the search, then took a step it would have taken
from the full one: a halved step passes the step test whenever the full one does, so the halvings could change only
the point returned, never the verdict.  `radau_transform=True` had looked 1.31x faster on this case for the same
reason: its rounding fell on the first halving.

**Built.**  `analysis.fsolve(line_search=True)`, without a `limiter`: a full step that passes the step test, to a
point whose OWN residual passes the residual test, is taken -- not halved.  `_ShootingNewton._matrix_free_newton` the
same (the point's residual scale, a matvec, made only where the step test holds and reused by the convergence test).
⚠ The trial stays.  Built first without it -- the convergence test run before the trial, the undamped loop's verdict
at its point, 2.90x and 2.04x -- the suite refused it: the convergence test reads the step and the PREVIOUS point's
residual, and only the trial sees the point returned.  On a period map with a jump (a state reset landing on a grid
point, `test_a_state_reset_needs_no_saltation_but_grid_alignment_is_a_cliff`) a 3.4e-10 step crossed it to a residual
of 8.27e-3, and the solve returned that point as converged; the line search had evaluated it and halved back to
1.7e-10.  With the trial kept that step does not settle and is halved as before.

**Tests.**  test_transient_branch_and_dae: a full step at a growing rounding floor settles -- three residuals, the
undamped point bit for bit (seven before); a step across a jump in the map is still halved to the side the solve
stood on (it fails under the version without the trial).  test_shooting_pss: the matrix-free Newton's floor the same
(seven builds before, three now).

Measured, separate processes, HEAD (a worktree) and this tree interleaved: Radau 42.93 / 42.69 -> 20.22 / 20.33 s,
and again the box idle (load < 1.5) 43.50 / 42.27 -> 20.14 / 20.23 s, 8 -> 4 traversals; gear 6.83 / 6.90 -> 4.48 /
4.54 s, 5 -> 3; every one of the 60 nodes' orbits within 4.7e-14 V
(Radau) and 6.3e-14 V (gear) of HEAD's on a 3.3 V scale, IM3 -63.57360 and -63.57383 dBc both, to five decimals.

**Gate**: the suite in the main tree, 5 795 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 193: every digest identical.  No C changed.

## 2026-10-23 — 195. Radau's PSS on the walk's KLU: the C stage Newton solves its block on it, the walk's step factors the C's last assembly, the converged replay is the last shooting walk, the memo takes the C's slots -- a Howland testbench's Radau PSS 2.33x (20.32 -> 8.73 s), the orbit within 5.2e-14 V

**Why**  Re-profiled after 194 (Andreas: "continue with next", then "implement all 4 items"), the testbench's
default Radau PSS (90 unknowns, a 270-unknown stage block, 3200 steps a period) took 20.1 s: the C stage
Newton 12.0 s (its dense `dgesv` on the 270 block ~470 us a step, its glue ~200), the sensitivities' steps
4.6 s (the block made again, 1.3 s, and factored on KLU, 2.5 s), the converged replay 3.6 s -- a walk of the
period the solve had just made.

**Built.**  (1) `_InnerTransient._new_transient` gives a frozen-grid, fully implicit traversal the walk's
StructuralKLU (`_step_klu(3 m)`) as `tr._stage_klu`; the coupled stage Newton solves its block on it
(`StructuralKLU.solve`: a refactor on the numeric it holds, a fresh factor on a zero pivot, the residual's
check), and the C stage Newton (radau.c) makes that sequence (pcx_kluseq.h, as the DC Newton's) with
StructuralKLU's numeric handed in and out (pcx_klu.h: a hand-back leaves the entry's).  The core's pattern
lies within the KLU's in each of the nine blocks (checked when it is bound), so the block is zero off it:
with the KLU bound an assembly makes the pattern's entries alone, in its order, no dense block, scan or
gather -- and, for their flags, the core's entries on the reference row and column, which the Python makes
before it drops them.  (3) The last assembly, at the converged stages, is the block the walk's step makes
from them, the same bits: a live stage step takes its values (`_tran_radau_c.final_block`, the very arrays
the C returned, the same `h` and tableau) and factors them (`StructuralKLU.factor_values`, `factor`'s own
second half).  (2) The shooting's period map records a stage walk's states and stage states as the
converged replay would; the replay takes the record where the solve returned the very point that walk
started from, at the same period and grid (`_stage_walk_reused`), and walks otherwise.  (4) The device
memo records the C's output slots themselves, views of a pooled set no earlier record holds a view of (its
arrays' reference counts as made; at most 3 sets and 48 MB), copies where none is free.

**Tests.**  test_radau_c: on a stage KLU the C Newton is the Python's (solution, statistics, warnings, memo,
the KLU's counters and numeric) on the four fixtures; a KLU without the core's pattern declines; a
hand-back after the C's refactors (the assembly cap at two) leaves the Python's sequence; the memo holds
the output slots and never one still held, copies where no set is free.  test_shooting_step_klu: on the
walk's KLU the orbit and the monodromy are the bytes of the Python Newton's and `_stage_block`'s; the
converged replay is the last walk (the waveform and PAC as walked, to the rounding); a record at another
point is not taken.  test_pinned_pairs: 'the coupled stage Newton and its C', 'the stage paths' passes and
the core's', 'the stage block and its one pass', 'the stage step map and its C', 'the Newton solve and its
C' re-recorded; 'the walk's stage block and the C stage Newton's last assembly' and 'the converged stage
replay and the shooting's last walk' recorded.  test_c_abi: the unit embeds pcx_klu.h.

Measured, the box idle, separate processes, HEAD (a worktree) and this tree interleaved, two rounds:
Radau 20.40 / 20.32 -> 8.73 / 8.77 s (3 walks of the period where 4 were: the replay taken); gear 4.49 /
4.56 -> 4.55 / 4.60 s, its orbit bit-identical (none of it is gear's); every one of the 60 nodes' Radau
orbits within 5.2e-14 V of HEAD's on a 3.3 V scale (1.6e-14 relative), IM3 -63.57360 dBc both.  The steps
on the way, one process: KLU in the C Newton 20.1 -> 15.5 s, the converged block 14.3, the replay 11.5,
the pattern-only assembly 9.35, the memo's views 8.74.

**Gate**: the suite in the main tree, 5 809 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 194: every digest identical; the sanitized fast tier (scripts/sanitize_suite.sh) 5 557
passed, no sanitizer report: clean.

## 2026-10-23 — 196. A Newton that fails at its iteration limit reports the solver's own tests -- the residual against `reltol (|J| |x| + |f|) + abstol` and a stateful limiter that held it (a user's 4 mV step reported a residual 1e8x over a tolerance it met 1e5x under; the update had failed)

**Why**  A peer's step test (a 4 mV step into an amplifier macromodel, tolerances reltol 1e-9, iabstol 1e-13)
fell to the minimum step at the step's corner with "StandardNewton failed to converge after 100 iterations;
residual worst at vgate: |f| = 1.319e-06 against a tolerance of 1.013e-13".  The report judged the residual
against `reltol |f| + abstol`; the solver's test is `reltol (|J| |x_next| + |f|) + abstol`, and at h = 1e-18 s
a 354 pF capacitor's companion makes `|J| |x|` ~ 1e8: the residual met its test ~1e5x under, and what had
failed -- the update of a 0 V source's current inside the macromodel, changing by 7.1e-8 A each iteration, about one
unit of rounding in the voltages times C/h -- came second, unnamed.  (Why the update cannot settle at that h is
a separate matter, not this entry's.)

**Built.**  `_worst_row_report` takes the residual's scale `I_scale` (the solver's `|J| |x_next| + |f|`) and
judges against it; without one it says the tolerance is "at least" `reltol |f| + abstol`.  A stateful limiter
that limited the last step (the solver's noncon) is named: both tests can hold at a point that is not a
solution.  StandardNewton and DampedNewton pass their own (DampedNewton the step it took, `alpha dx`).  The C
Newton (newton.c) makes the last iteration's `|J| |x_next|` at the limit -- its test makes it only where the
update held -- and `_tran_newton_c` and `_dc_newton_c` pass it with `|F|`, so the C's message is the Python's.

**Tests.**  test_nrsolver_variants: a Jacobian ~1e8 with a weak direction and a residual of 1e-6 -- inside its
own test, the step along the weak direction over the update's -- names the update and not the residual
(standard and damped; the old report named the residual 1e6x over); a limiter that limits every step names
itself.  test_pinned_pairs: 'the Newton solve and its C', 'the DC base solve and its C', 'the stepping loop
and its C' re-recorded.

**Gate**: the suite in the main tree, 5 812 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 195: every digest identical (a report changed, no iterate); the sanitized fast tier
(scripts/sanitize_suite.sh) 5 573 passed, no sanitizer report: clean.

## 2026-10-23 — 197. A Newton update that has stopped contracting is judged against the rounding its own solve carries -- a user's 4 mV supply step, which fell to the minimum step at the step's corner, completes in 0.07 s (the 2 mV step's waveform within 1.6e-11 V)

**Why**  The step test of 196 (a 4 mV step into an op-amp macromodel's supply; reltol 1e-9, iabstol 1e-13) failed
because one update could not settle: a 0 V source's current in the macromodel, set by the KCL of a node a 354 pF
capacitor loads.  At a small step the capacitor's companion makes that KCL a difference of terms ~C/h |v|, and
the residual is evaluated to eps of them, so the current moves by ~eps C |v| / h every iteration -- 7e-8 A at
h = 1e-18 s against a tolerance of 1e-13 A -- while the residual's own test holds.  It crosses the tolerance at
h ~ 1e-12 s, where the failures began; each retry halved h and doubled the noise, down to `minstep`.  The current
is not a capacitive loop's (the node's other side is a state through 1 ohm), so the error test's 'charge' mode is
right to leave it alone: what failed is the Newton's update test, which no tolerance setting reaches.  The 2 mV
step converged only because each such Newton happened to draw a noise sample under the tolerance.

**Built.**  `StandardNewton`'s fourth condition (`nrsolver._within_rounding`): where the residual's test holds,
no stateful limiter limited, and the update fails and has stopped contracting (its largest entry not under half
the last one's: `_step_size`), every unknown failing the update test passes if it moved by less than its
tolerance plus ``eps (|J^-1| I_scale)_i`` -- the residual's own scale (``|J| |x_next| + |F|``) through row i of
J^-1.  The rows by transposed solves with the failing unknowns' unit vectors, `FLOOR_ROWS` (32) at a time, on the
factorisation the step was just solved with (`LinearSolver.solve_transposed`, new: numpy's solve of J^T on a
dense J, `klu_tsolve` on the numeric `StructuralKLU` refactored; numpy's solve of J^T where a scaler
equilibrated the step), each floor ``eps * np.dot(|row|, I_scale)``.  A Newton that converges by (1)-(3) is
untouched; the floor is consulted only on a stall.  Accepted, the point is as good as the arithmetic allows: the
residual's test, which holds, already bounds its error by reltol / eps (~4.5e6) times this floor.  The C Newtons
make it bit for bit: newton.c (`step_size`, `floor_holds`; the dense rows by dgesv on J reduced row-major --
J^T's column-major copy -- the KLU rows by klu_tsolve on the run's numeric; the floors by numpy's own ddot) and
dcnewton.c (the same, its flags read before the verdict); `_tran_newton_c._chain` and `_dc_newton_c` take the
Python where either helper is patched.

**Measured** on the peer's case (scratch only; its macromodels are vendor files), every failing Newton before
this change: over 217 000 stalled iterations the observed moves were 0.023 (median) .. 0.27 (max) of the floor;
the moves above it were genuine first steps from the predictor, which the floor refuses (it needs a stall).
After: the 4 mV step passes (309 steps, 15 rejected, 666 Newton iterations; 0.07 s with the C, the Python's
bytes); the 2 mV step's three node waveforms move by at most 1.6e-11 V (tolerance ~2e-9 V), 515 Newton
iterations against 525.

**Tests.**  test_nrsolver_variants: an update at its residual's rounding converges at the second iteration, one
11x over the floor is refused; a contracting update never consults the floor; one block or two, the same answer.
test_structural_klu: the transposed solve is on the numeric the solve left (numpy's to the rounding, nothing
factored or counted again; a foreign matrix answered by the fallback).  test_newton_c: a VCVS-and-354 pF chain at
femtosecond steps, dense and on KLU -- the C's bytes, statistics and warnings the Python's, whose Newtons the
floor converged.  test_dc_newton_c: the stalled step as an operating point (a behavioural g V for the companion, a
current source for its history), dense and on KLU (130 cells, five blocks) -- the C's bytes the Python's, and
neither converges without the floor.  test_pinned_pairs: 'the Newton solve and its C' and 'the DC base solve and
its C' take the floor's pieces and are re-recorded, 'the stepping loop and its C' (newton.h) and 'the transform
solve and its frozen form' (AutoSolver) too, and thirteen more whose C includes pcx_abi.h (a ddot typedef: each
as recorded with the header as it was).

**Gate**: the suite in the main tree, 5 820 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 196: every digest, step, rejection and Newton count identical -- the floor fires on none
of them; the sanitized fast tier (scripts/sanitize_suite.sh) 5 516 passed, clean (the tier deselects by the last
full run's durations: fewer than at 196), and the four C-twin cases of the floor under `--deep`, clean.

## 2026-10-23 — 198. A Radau PSS's walks refactor one KLU numeric, take a stage time's sources once a solve, and seed the three stages at once -- a Howland testbench's Radau PSS 1.27x (8.88 -> 6.99 s), its gear PSS 1.09x (4.52 -> 4.14 s), every orbit bit-identical

**Why**  Re-profiled after 197 (Andreas: "start as is"), the testbench's default Radau PSS took 8.8 s, three walks of
3 200 steps; the C stage Newton is two assemblies a step already.  Around it, per step: the sensitivity factored
the C's last assembly afresh (klu_factor, 167 us) and freed it at once; the sources at the three stage times were
assembled again every walk (the step's memo lives one step: 28 800 vectors for 9 600 times); the three stage
seeds each gathered the predictor's nodes, the last step's motion and the nodes' values anew; the coupled solver
walked the circuit for its stateful limiters every step (45 us) though the run holds them; the accept re-sorted
the predictor's whole node list.

**Built.**  (C) `StructuralKLU.walk`: a `_KLUWalk` holds one numeric for a walk -- the first block factored, every
later one refactored on it (klu_refactor, 57 us) and checked as `StructuralKLU.solve` checks its refactor, by a
residual (of a probe: the block times ones) against REFACTOR_RESIDUAL_TOL, a fresh factor where it fails or meets
a zero pivot; a factorisation it handed out refuses to answer after the next block's (`_KLUBorrowed`).
`_walk_stage` gives it to every coupled step of a walk that keeps none (the shooting Newton's monodromy:
`_final_lu`, `_walk_factor`); a kept walk -- PAC's and pnoise's factored period -- factors each step on its own.
(B) `_stamp_plan.source_memo(cir)`, opened by `PSS.solve` for the whole solve: `source_direct` keeps each vector it
makes by time and analysis and hands a copy for that time asked again, before its checks -- the opener vouches
that nothing a source reads moves within the solve (a sweep writes between solves, which drop the memo).  (A
memo on the program itself, renewed when its checks re-arm, never hit -- the program re-arms at every walk's
start -- and could not have seen a time function replaced in place.)  (A) `_coupled_stage_solver` reads the run's
`_stateful_lims`, as `_newton` does.  The extension (Andreas: "extend first"): `_pred_stages`, the three seeds at
once -- per target `_pred_or`'s bytes and counts, the nodes, motion and values made once (`_fit_fast` takes the
values); `_pred_note` merges new nodes after the history without sorting it again (the same list).  Not built: the
branch screen of the stages in C (`branch_c:large`: no pattern mark on the Radau C's matrices) -- served whole,
it measured no faster than numpy's reductions; the time is the Python around it.

**Tests.**  test_structural_klu: a walk refactors one numeric and checks it (a zero pivot factored afresh, a stale
factorisation refused); the transposed solve as before.  test_shooting_step_klu: a shooting walk factors once and
refactors 19 times in 20 steps, no check failing, the orbit and the monodromy within the rounding of a fresh factor
every step.  test_source_direct: a solve's memo hands a copy of the vector made at a time (nested, scoped, the
counts); a PSS by radau and by gear makes each time's sources once, its orbit the bytes of one without the memo;
the PSP stage's direct-pass count split into served and memo.  test_radau_c: the coupled solver reads the run's
stateful limiters.  test_predictor_fast: the seeds at once are the per-target ones, bytes and counts (hypothesis,
targets behind and ahead of the history); a patched `_pred_or` sees every stage; the accept's merge is the whole
sort's list (hypothesis, new nodes before, among and after the history, within the dedupe's tolerance).
test_pinned_pairs: nine re-recorded after checking their twins -- 'the walk's stage block and the C stage
Newton's last assembly', 'the sources and their direct pass' and 'the stage predictor and its multistep fast
path' take the new pieces; 'the stepping loop and its C' (`_pred_note`: loop.c's history is the same list), 'the
predictor and its C' (`_fit_fast`'s arithmetic as before), 'the coupled stage Newton and its C', 'the stage paths'
passes and the core's', 'the stage block and its one pass', 'the stage step map and its C' (its C declines the
KLU blocks).

Measured, HEAD (a worktree) and this tree interleaved, separate processes, load ~2: Radau 8.84 / 8.92 -> 6.91 /
7.06 s, gear 4.52 / 4.51 -> 4.12 / 4.15 s; the whole orbit (`tpss.x`) the same bytes, both methods.  The pieces,
one at a time against 8.8 s: the walk's numeric -1.2 s, the solve's sources -0.5 s, the stateful limiters ~0.1 s,
the seeds at once, the accept's merge and the memo before the checks -0.35 s.

**Gate**: the suite in the main tree, 5 829 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 197: every digest, step, rejection, Newton count and path count identical.  No C changed.

## 2026-10-23 — 199. A shooting walk's step on its KLU in one C call -- the refactor, its check and the stage system's solve (`csrc/sensklu.c`) -- a Howland testbench's Radau PSS 1.13x (6.94 -> 6.15 s), every orbit bit-identical

**Why**  The plan after 198 (g24/plan_stage_walk.md; Andreas: "start with the next lever"): the Radau PSS's step was
~730 us, ~300 of it arithmetic (the C stage Newton 179, KLU, BLAS) and ~400 Python around it.  Its first stage, the
sensitivity: the plain walk's step on the walk's KLU (198's `StructuralKLU.walk`) took 243 us -- `_stage_step` 152
(the refactor 57, its probe's check 41 in eight small numpy calls, `_C_at` 15) and `_StageStep.solve` 91 (the
solve of the stacked `C_n P`, 19 live columns, 62) -- and the dense C step map (`_sens_c`) declines KLU blocks.

**Built.**  `csrc/sensklu.c` (`pcx_sensklu`) and its glue `shooting/_sensklu_c`: `block_step(walk, Ax, base, s)`
makes `_KLUWalk.factor_values` of the C stage Newton's last assembly (`_tran_radau_c.final_block`'s values) --
the numeric refactored, the probe in `np.bincount`'s order with numpy's max and Python's `max`, a fresh factor
where the refactor or the check fails, the first block factored -- then the stage system's right-hand side,
`C_n P` (numpy's product, made in the glue as `_StageStep.solve` makes it) once a stage, laid column-major as
`_KLUNumeric._run` copies it, and klu_solve of its columns at once; the new `P` is the last stage block of the
transpose of a C-ordered buffer (the Fortran strides `_run` returns).  The walk's numeric, generation and counts
and the solver's `numerics` are left as `factor_values` leaves them; a fresh factor that fails raises its
LinAlgError, a failed solve `_run`'s; a flag numpy would warn of in the check hands the step back (after a
refactor in place, which the Python's repeat makes again).  `_walk_stage` asks it for a KLU walk's plain step
(the dense map keeps the dense walks); `step` returns, with `P`, the step on a full-width `P` on this step's
factorisation (`_LiveColumns.advance`'s non-finite case).

**Tests.**  test_structural_klu: a walk's blocks through the C and through `factor_values` and the stacked solve --
the first factored, refactors, a zero pivot factored afresh, a singular block raised with the numeric gone, one
after it, and three stages -- the solution's bytes and strides, the exceptions and every count the Python's.
test_shooting_step_klu: a PSS with the C step map and without -- the orbit and the monodromy the same bytes,
every walk's counts and the solver's `numerics` the same; the large-block test reads the KLU map's count.
test_c_abi: the unit.  test_pinned_pairs: 'the KLU step map and its C' recorded.  doc/transient_runtime.md: its
row; doc/shooting_history.md: its entry.

Measured, HEAD (a worktree) and this tree interleaved, separate processes, load ~2.5: Radau 6.91 / 6.97 -> 6.17 /
6.13 s (1.13x, over the stage's 1.10x line), gear 4.13 / 4.09 -> 4.08 / 4.08 s (none of it gear's); the whole orbit
the same bytes, both methods.  9 599 of 9 600 steps served (the run's first: no C last assembly yet).

**Gate**: the suite in the main tree, 5 838 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 198: every digest, step, rejection, Newton count and path count identical; the sanitized
fast tier (scripts/sanitize_suite.sh) 5 553 passed, clean, and the KLU step map's ten tests under `--deep`, clean.

## 2026-10-23 — 200. A shooting walk's plain stage steps in C, a chunk a call -- the seeds, the C stage Newton, the screen, the finish and roll, the step map on the walk's KLU (`csrc/stagewalk.c`) -- a Howland testbench's Radau PSS 1.92x (6.10 -> 3.17 s), every orbit bit-identical

**Why**  S2 of the stage-walk plan (g24/plan_s2_stagewalk.md; Andreas: "Do both in the order you think is best"):
S2b, the whole plain walk in one C call, first -- S2a's Python trims proved smaller (the guards stamp-based
already, the memo puts twin-compared).  After 199 a plain Radau walk's step was the C stage Newton, the C step
map and the predictor's numpy fit with the Python between them -- the inner transient's step chain, the stage
context and sources, three screens, the finish and roll, the walk's loop, `C_n P`, `advance` -- about as much
again as the C.

**Built.**  `csrc/stagewalk.c` (`pcx_stagewalk`, contract `stagewalk.h`) and its glue `shooting/_stagewalk_c`:
a chunk of steps a call (16, doubling to 256), each made as the inner transient's step and the walk's step map
make it -- pred.c's three seeds (`_fit_fast`'s twin), radau.c as `_tran_radau_c.solve` calls it (its readiness
checks factored out as `_ready`), branch.c on each stage's `C` (the core's pattern, the run's certificate),
sensklu.c on the walk's numeric with `C_n P` numpy's own dgemm / dgemv, the live columns' drop; the finish and
the roll (the rings, `_pred_note`'s merge in Python's stable-sort order with its dedupe and cap) in C.  The
stage times and sources for a chunk are made first, in the Python's order.  The state is read from the Python
walk and the inner transient before a call and written back after it as the Python leaves it: the rings,
`_pred_hist`, `_q_cache`, `_rk_Y`, `_iq`, `_Cmat`, `_Geq`, `_C_cache`, `_branch_cmax`, `_dt` and its history,
the device memo's last step's records, both numerics with their counts and generation, the live columns, each
step's records.
  A stop is exact three ways: before the step's Newton nothing of it is made; after a converged Newton its
result is held (`_tran_radau_c.hold`) and the next `solve` serves it once where its inputs are the same bytes
and the stage KLU and limiter walk as held -- without it, a Newton that factored afresh is repeated as a
refactor of what it left (forced: refactors / factors 170 / 1 against the Python's 168 / 3); after the step map
the step is made and the glue does what the Python does next (`advance` at full width, the LinAlgError).  The
step map writes to the one of two buffers the carried `P` is not in; a point whose reference entry is not +0.0
stops (the next step starts from the inserted +0.0); a constant `C` is screened (the Python's kept verdict is
the same screen).  `_walk_stage` resumes the C after a step it stopped before, and spaces out calls that made
none (1, 2, 4 ... 64 Python steps).  `_sensklu_c.full_step` factored out of `step` for the full-width step.
Serves what the C stage Newton serves on the walk's KLU: no stateful limiter, no PCNR, the cost transform off.
No function of the unit or its glue over 80 lines (both now in test_runtime_map's lists): the C walk is its jobs
(the node order, the seeds, the screen, the step map, the commit, the live drop), the glue its own (`_eligible`,
`_screen_ready`, `_state_ready`, `_bind`, `_table`, `_chunks`, `_settle`, `_write_back`, `_after`) -- the split
timed against the whole: 3.17 against 3.16 s.

**Tests.**  test_stagewalk_c (10): a PSS with the C walk and without -- the orbit, the monodromy, the replay's
records, the inner transient's state, every walk's and the stage KLU's counts the same -- on a linear ladder and
a MOSFET chain; Python steps between calls and chunks of one to three; the C Newton's assembly cap handing steps
back mid-walk; a stop after a fresh-factor Newton served from the hold, and with the hold spoiled or dropped the
counts parting; a NaN in a carried column (`advance` at full width); a singular step map raised with the same
state; the declines (switch, a patched piece, the transform); one ulp planted in the C's points showing in the
orbit.  test_shooting_step_klu: its three unit tests switch the stage walk off.  test_c_abi: the unit.
test_runtime_map: the unit and its glue in the long-function lists; test_paths: `walk` among the declining.
test_pinned_pairs: 'the stage walk and its C' recorded; 'the coupled stage Newton and its C' and 'the KLU step
map and its C' re-recorded on their twins (`_ready`, the held serve; `full_step`), the references unmoved.
doc/transient_runtime.md: its rows; doc/shooting_history.md: its entry.

Measured, HEAD (a worktree) and this tree interleaved, separate processes, load ~1.5-3: Radau 6.10 / 6.09 / 6.17
-> 3.19 / 3.15 / 3.17 s (1.92x, over the stage's 1.4x line), gear 4.06 -> 4.09 s (none of it gear's); the whole
orbit the same bytes.  9 599 of 9 600 steps made in C (the run's first: the stage KLU not yet analysed).

**Gate**: the suite in the main tree, 5 856 passed (6 skipped, 3 xfailed); the 27 gold decks (no PSS), transient and
standalone DC, against 199: every digest, step, rejection, Newton count and path count identical (DC but for its
timing); the sanitized fast tier (scripts/sanitize_suite.sh) 5 562 passed, clean, and under `--deep` the stage
walk's, the step KLU's and the radau C's tests (57, the stage walk's ten among them), clean.

## 2026-10-23 — 201. A plain Radau transient's Python trimmed -- the estimate from the memo and the run's limiters, the seeds from pred.c, the screens in C on the core's pattern, the sources from the C table -- a Howland testbench's Radau transient 1.095x (0.877 -> 0.801 s), every answer bit-identical

**Why**  Andreas, on S2a after 200 ("Both, trims first"): S2a as planned (the PSS step's Python) was moot -- the C walk
replaced it -- and a plain Radau transient, which the C walk does not serve, spent 2.22 ms an attempt on the Howland
testbench (91 unknowns, adaptive, from DC): 1.47 ms in the C stage Newton's dense LU, the rest Python -- the embedded
estimate 224 us, the sources 126, the seeds 103, the screens 67 (in situ wrappers, every path count unmoved).  Plan
g24/plan_s2a_trims.md: four bit-identical trims, ~1.14x expected, LINE 1.08x.  The stage block on KLU (~2x, a
rounding-level change) follows as 202.

**Built.**  T1, `_RadauStages._radau_error_estimate`: the run's stateful limiters (`_stateful_lims`, as 198 made the
stage solver read them) instead of a walk of the circuit's elements every attempt, and x_n's C, G and i from the device
memo, which the step that reached x_n recorded (no stateful limiter: the memo's premise -- the core's values), else the
core's passes as before (`radau.est:memo`, 302 of 399 attempts).  T2, `_StagePredictor._pred_stages`: the three seeds
from pred.c (`_tran_pred_c.stages`: `_predict_fast`'s pinned twin a target, on a state kept on the transient, `_pred_c`;
the count `pred:fast` the Python's); pred.c hands back an underflow only where numpy's error state reports one
(`pred_t.under_ok`, set by `stages` and the stage walk; the C loop's and `predict`'s unchanged) -- the Howland's tiny
states raised it on 323 of 399 calls.  T3, `_tran_radau_c.solve`: the three stage C's the screen reads next (the memo's,
every entry the core's) marked with the core's pattern (`_C_pattern` a tuple of them; `_tran_branch_c.screen` takes one
of them), so the C screen serves where a large dense C declined.  T4, `_stamp_plan.source_direct`: the vector from the C
sources' table kept on the direct pass's program (`_direct_c`) where the program's stamp holds -- `_SourceDirect.arm`
now watches each source's values and time function, the table checked once a stamp (`holds`), a PWL (whose points no
watcher sees) the Python's; the table made and called quietly (`quiet`), so a watched run counts what an unwatched one
does, `u:called` the pass's own.  `_RUN_STATE` declares `_pred_c`, and `_stagewalk_c` (200's, set only on a shooting
walk's inner transient, which the registry's loops never observe).
  ⚠ FOUND ON THE WAY, A WRONG ANSWER IN THE C LOOP: `_tran_sources.table_for` made a time function's CLASS waveform
where the function carried its own `f` (a user's injected noise: `IS.function.f = ...`) -- the C stepping loop ran such
a source as its class's (a DC 0: a noisy RC's adaptive gear run 74 steps where the Python takes 124, another waveform),
since the C loop's sources (2026-10-10).  T4's table met it on test_shooting_sampled's noisy-transient jitter (the
scatter 2.5e-18 s where 3.3e-11).  `table_for` now declines a function with its own `f` or `next_event`
(`sources_c:own`), for the loop and the direct pass alike.

**Tests.**  test_radau_trims (8): each trim switched off -- pred.c, the C screen (a 68-unknown chain, past the dense
screen's size), the C sources -- the same run (waveform, sweep, statistics, warnings, the device memo, the counts the
Python makes); the estimate's memo read made unavailable, the same run; the run's limiters read once; a source's value
and its time function's attribute written between runs, each run the C-free one's; a PWL source the Python's; pred.c
through an underflow, the Python's seeds, and handed back under `np.errstate(under='warn')`.  test_tran_sources: a
time function of its own declined, and run through the C loop and the direct pass, the Python's.  test_pinned_pairs: the
moved pairs re-recorded after checking (the coupled stage Newton, the branch screen, the predictor, the stage predictor's
fast path, the stage paths' passes, the sources and their C, the sources and their direct pass, the stepping loop (pred.c
in its twin), the stage walk); `stages` and `_direct_c` joined their pairs' twins.

Measured, HEAD (a worktree, 85190140) and this tree interleaved, separate processes, load ~2.6-3.8: the plain Radau
transient (4 periods) 0.883 / 0.877 / 0.875 -> 0.801 / 0.813 / 0.800 s (1.095x, over the 1.08x line); the Howland's
Radau PSS 3.19 -> 3.05 s, gear 4.12 -> 4.07 s; every waveform and statistic the same bytes.

**Gate**: the suite in the main tree, 5 865 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone
DC, against 200: every digest, step, rejection, Newton count and path count identical (DC but for its timing); the
sanitized fast tier 5 535 passed, clean, and under `--deep` the trims', the predictor's, the sources', the screen's and
the stage walk's tests (88), clean.

## 2026-10-23 — 202. A plain run's coupled stage block on a StructuralKLU from 250 unknowns (default on) -- a Howland testbench's Radau transient 2.09x (0.801 -> 0.383 s), the answers ~1e-8 V from the dense LU's

**Why**  The second half of Andreas's "Both, trims first": after 201 a plain Radau transient's attempt on the Howland
testbench spent ~1.46 of ~2.0 ms in the C stage Newton's dense LU of the 270-unknown stage block, which the shooting
walks already solve on a StructuralKLU (`_InnerTransient._step_klu`).  Prototyped (g24/klu_proto.py): 2.1-2.5x -- and
the answers NOT to the rounding, as Andreas had been told: a fixed grid's waveform 1.13e-8 V apart on 3.3 V over four
periods.  Asked again with those numbers, Andreas: "Build, default on".  Plan g24/plan_202_stage_klu.md, LINE 1.5x.

**Built.**  `_RadauStages._run_stage_klu`, called by `_open_run` (a plain run's start, never a PSS walk's): for a fully
implicit tableau on the numeric toolkit with no `linearsolver=` given, from `MIN_N_FOR_SPARSE` (250) block unknowns,
a StructuralKLU on the circuit's structural pattern tiled over the stage blocks (the shooting walk's), made once per
topology and size (`_stage_klu_run`), set as `_stage_klu` -- which the C stage Newton and the Python's already solve on;
a caller's own `_stage_klu` (a PSS's) never replaced; `PYCIRCUIT_RADAU_KLU=0` keeps the dense solve.  `_RUN_STATE`
declares `_stage_klu_run` and says who sets `_stage_klu`.

**Tests.**  test_radau_stage_klu (3): an 84-stage MOSFET chain (its block past 250) on its KLU, the C serving on it,
the fixed-grid waveform within 1e-9 relative of the dense run's and an adaptive run's steps the same; the C and the
Python stage Newton the same bits on the run's KLU; none for a small block, a multistep method, a given
`linearsolver=`, the switch off; a caller's own kept; one per topology, kept from run to run.

Measured, HEAD (a worktree, aad6c5f3) and this tree interleaved, separate processes, load ~1-3: the Howland plain
Radau transient (4 periods, adaptive) 0.808 / 0.797 / 0.801 -> 0.383 / 0.383 / 0.394 s (2.09x, over the 1.5x line), the
same 302 steps, `max_step` 1e-9 relative apart, the waveform within 9.6e-10 V at the dense run's times; on a fixed grid
(12 800 steps) 1.13e-8 V apart.

**Gate**: the suite in the main tree, 5 868 passed (6 skipped, 3 xfailed); the 27 gold decks (gear: none takes the
change), transient and standalone DC, against 201: every digest, step, rejection, Newton count and path count
identical (DC but for its timing); the sanitized fast tier 5 614 passed, clean, and under `--deep` the stage KLU's, the
radau C's, the trims' and the stage walk's tests (59), clean.

## 2026-10-23 — 203. Every C path against a method set on an element's instance -- one more wrong answer: the PCNR operating point's C made a participant's class kernels where the instance held its own `i` or `G`

**Why**  Andreas, "to 1, 2 and 3 then continue to 4" -- 1, the audit 201's sources bug called for: a C path that
picks its work by an object's CLASS while the Python calls a method set on the INSTANCE gives another answer, silently
(201's: a time function's own `f`, the C loop's sources making its class's waveform for two weeks).

**Measured.**  A sweep (g24/ovr_sweep.py, ovr_drive.sh): six circuits -- HDL MOS level 1 and gate-charge chains, a
built-in diode with its stateful limiter, a charged SPICE diode, Gummel-Poon chains plain and charged -- x DC, gear and
Radau transients, Radau and gear PSS, with an instance override planted on a device or a source (`i`, `q`, `G`, `C`,
`u`, `dudt`, `limit`, `next_event`, `accept_step`; values moved by 1e-3, limiting removed, calls counted, breakpoints
added), each run as is against the run with every twin claimed bit-identical switched off (their env switches; the dict
watcher and the algorithmic choices as they are): 820 runs.  30 differed, every one an HDL device's `i` or `G` in DC
or in a transient opened on its DC point -- bisected to the PCNR operating point's C (`_pcnr_c`), whose `_class_ok`
checked the class's methods but not the participant's instance.  Every other override (charges, sources, limiters,
events, the step hook) agreed, each shown to take effect in its analyses.

**Built.**  `_pcnr_c._record` declines a participant whose instance holds its own `i`, `G` or `pcnr_limit` (`_OWN`,
`pcnr_c:own`): the Python PCNR reads them (`PcnrDevice._kernel_stamp`'s `el.i`, `el.G`; `PcnrDevice.limit`'s
`pcnr_limit`).  The sweep again: 820 runs, none differing.

**Tests.**  test_instance_overrides (12): a charged SPICE diode and a gate-charge MOSFET chain x DC, gear and Radau
transients, Radau and gear PSS x eight plants, every twin on against every twin off, in-process; alive: the decline
emptied, exactly the device's `i` and `G` differ; the decline counted and the operating point moved.
test_pinned_pairs: 'the PCNR DC solve and its C' re-recorded (the decline).

**Gate**: the suite in the main tree, 5 880 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone
DC, against 202: every digest, step, rejection, Newton count and path count identical (DC but for its timing).  No C
changed (the sanitizer tier not run).

## 2026-10-23 — 204. A shooting walk's plain multistep steps in C, a chunk a call -- the C stepping loop's own attempt on the walk's grid and the step map after it (`csrc/loop.c`'s walk) -- a Howland testbench's gear PSS 1.86x (4.05 -> 2.18 s), every orbit and count bit-identical

**Why**  Andreas, "to 1, 2 and 3 then continue to 4" -- 2, a C walk for gear PSS.  After 193 and 198 a gear walk's
step was the Python's inner-transient step around the C Newton, and one C call for its step map: on the Howland
testbench the walks were 3.24 of the PSS's 4.03 s, the C Newton 0.55 and the map 1.03 of those -- the rest Python.

**Built.**  A walk mode of the C stepping loop (`csrc/loop.c`: `walk`, `walk_map`, `walk_dies`; loop.h's THE WALK) --
the loop already made the multistep attempt bit for bit, so no second unit: on the walk's grid each step is `take`
(the end time the grid's, `have_tn`) and `accept`, every step accepted, then the step map (lmmsens.c through its own
struct).  Its glue `shooting/_lmmwalk_c`: a chunk of steps a call (16, doubling to 256), entered where the C Newton's
last call served at the watch counter as it stands (as the C loop enters), its checks the C loop's for the run and the
walk's own (a plain dense walk, its pieces read off the walk object's class, the step map's toolkit and solver, the
monodromy's arrays as the step map takes them).  Each chunk's sources are the Python's own: `_source_at` at each step's
end in a step memo of its own (`_tran_newton_c._source`'s pass), so the solve's source memo is read, written and
counted as the Python's steps do; the rows of steps not made are taken back (`_undo`).  The state is read before a
call and written back after it as the Python leaves it -- the C loop's write-back with the step lengths the grid's own
objects (numpy floats, as `_dt`, its history and the coefficients hold them), the walk object's `_Jf`, `_Geq`, `_C`,
`_coeffs`, the monodromy, each step's record.  A prune that drops a column is the Python's (`advance`: its views,
which the Python's next maps decline as 'layout' and take in numpy); a map handed back (singular, not finite) is made
as `_walk_lmm` makes it.  `_walk_lmm` asks the C from where a Python step's Newton served; a call that made none is
asked again after 1, 2, 4 ... 64 Python steps.  `_tran_loop_c._bind` takes no sources table (the walk's, `wu`).

**Found on the way** (not changed here): the C sources table (`_tran_sources`, the C loop's and the direct pass's)
makes a waveform's class sine where the class's `f` is patched -- a plain transient with `func.Sin.f` patched gives
another answer with the table on than off.  205's.

**Tests.**  test_lmmwalk_c (13): a PSS with the C walk and without -- the orbit, the monodromy, the inner transient's
and the walk object's state (the step lengths' and coefficients' types with their values), the Newton KLU's counts,
every count but the walk's own -- under Gear-2 (dense and on the KLU), the trapezoid and Euler; a call cut short every
third step and chunks of one or two; the C Newton handing back steps keyed by their seed; the step map handing back
(singular, not finite) steps keyed by their J; every plain walk recording, the records the same; a drive whose own `f`
raises, the same error behind the same state; the declines (switch, a patched piece, `_record_cj`); one ulp planted
in a call's last point showing.  test_paths: `walk` among the declining; test_runtime_map: the glue among the long-
function list; test_instance_overrides: its switch in OFF.  test_pinned_pairs: 'the multistep walk and its C'
recorded; 'the stepping loop and its C' re-recorded on its twin (the walk mode, `_bind`'s sources) and 'the multistep
step map and its C' on its reference (`_walk_lmm`'s call), the other sides unmoved.  doc/transient_runtime.md: its
rows; doc/shooting_history.md: its entry.

Measured, HEAD (a worktree) and this tree interleaved, separate processes, load ~1.3: gear 4.07 / 4.05 / 4.00 ->
2.18 / 2.20 / 2.14 s (1.86x, over the 1.4x line), Radau 3.04 -> 3.01 s (none of it Radau's); the orbit the same bytes.
9 579 of the walks' 9 599 steps made in C (the rest: each walk's first served step and the steps a prune's views
leave to the Python).

**Gate**: the suite in the main tree, 5 896 passed (6 skipped, 3 xfailed); the 27 gold decks, transient and standalone
DC, against 203: every digest, step, rejection, Newton count and path count identical (DC but for its timing); the
sanitized fast tier (scripts/sanitize_suite.sh) 5 653 passed, clean, and under `--deep` the walk's and the C loop's
tests (49), clean.

## 2026-10-23 — 205. The C sources' table against a waveform patched on its CLASS -- one more wrong answer: the table made the genuine sine where `func.Sin.f` was patched

**Why**  Found building 204: a twin test patched `func.Sin.f` on the class (a drive that raises) and the patch never ran --
the C sources' table, which picks its arithmetic by a source's waveform class, declined an instance's own `f` or
`next_event` (201) but not the class's.  Checked: a plain gear transient of a MOSFET chain with `func.Sin.f` patched
gives another waveform with the table on (the C loop's sources, the direct pass's `_direct_c`) than off.  The same
family as 203's (a C path that does not see what the Python calls).

**Built.**  `_tran_sources`: the seven waveforms' own `f` and `next_event`, with `Sin.f`'s scalar path
`func._sin_scalar`, read once as `func` defines them (`_genuine`, `_paths.genuine`); `table_for` declines a source whose
class holds another (`_patched`, `sources_c:patched`); a table records the classes it was built under and is neither
held (`holds`, the C loop's entry and the direct pass's stamp) nor called (`_call`, every direct pass) once one of them
is patched (`Table.genuine`: a few `is` tests on its own classes).

**Tests.**  test_tran_sources (+5): `Sin.f`, `Sin.next_event`, `_sin_scalar`, `Pulse.f` patched on the class -- no table
under the patch, a held one refused (`sources_c:patched` twice), a gear transient (the C loop) and a Radau one (the
direct pass) the runs without the table, each moved by its patch; alive: the check emptied, the run differs again.
test_pinned_pairs: 'the sources and their C' re-recorded on its twin, `_genuine` and `_patched` added to it, the
reference unmoved.  doc/transient_history.md: its entry.

**Gate** (one run, 206's change beside it): the suite in the main tree, 5 902 passed (6 skipped, 3 xfailed); the 27
gold decks, transient and standalone DC, against 204: every digest, step, rejection, Newton count and path count
identical (DC but for its timing).  No C changed (the sanitizer tier not run).

## 2026-10-23 — 206. A PCNR rung that runs away fails unwarned -- numpy's "invalid value encountered in multiply" at the Schur sum was an overflowed diode derivative times a consistent probe's zero

**Why**  A peer session's report (its user approved sending the case; it carries vendor models, kept in scratch): a laser
driver testbench's DC printed numpy's "invalid value encountered in multiply" at `schur_reduce` (`col * g_lim`), and under
`np.seterr(invalid='raise')` the FloatingPointError was swallowed by `solve_dc_continued`'s ladder; the answer the same.
Its capture: by iteration 53 the rung's iterate had run away (to 1e174, a junction's limited unknown to -1.6e164), a
diode's `didv` column `[0, -inf, inf]` met its consistent probe's `g_lim` of 0.  It reproduces only on the peer's AVX-512
Zen 4, whose numpy `exp` differs from the C library's (the C kernels decline, the Newton takes another path); on this
box's i7-12800HX the same DC never runs away (f326cb88 and e5521385, the C PCNR on and off).

**Built.**  `pcnr.schur_reduce`'s sum under `np.errstate(invalid='ignore')`: an invalid operation there (`inf * 0`,
`inf - inf`) comes only from an overflowed derivative at a run-away iterate; its NaN is kept -- the same bytes -- and fails
the step at the solve or at `solve_dc`'s non-finite test as before; no warning.  One `errstate` a call, not a probe (0.7 us
each).  A finite check before `predict` was weighed and refused: it would change the C twin's served NONFINITE path (its
KLU counts) and needs a C change for no other gain.  The C twin still hands such a step back on its invalid flag.

**Tests.**  test_pcnr_vector (+1): the reported column on a four-row system -- the bytes numpy makes unwarned, no warning
and no error under `invalid='raise'` (fails on the old sum: `FloatingPointError`).  test_pinned_pairs: 'the PCNR DC solve
and its C' re-recorded on its reference (`schur_reduce`), the twin unmoved.  doc/transient_history.md: its entry.

**Gate** (one run, 205's change beside it): the suite in the main tree, 5 902 passed (6 skipped, 3 xfailed); the 27
gold decks, transient and standalone DC, against 204: every digest, step, rejection, Newton count and path count
identical (DC but for its timing).  No C changed (the sanitizer tier not run).

## 2026-10-24 — 207. A plain Radau transient's attempts in C -- the C stepping loop's stage family: a Howland testbench's Radau transient 1.96x (0.384 -> 0.196 s), every byte the same

**Why**  Item 3 of Andreas's "1, 2 and 3 then continue to 4": the default analysis, a plain Radau transient, spent about
half of each attempt in Python around the C stage Newton.  Measured at ba0bcd9a (a Howland testbench, 4 periods, 399
attempts): ~954 us an attempt -- the C stage Newton 439, the estimate 160 (its dense solve 36), the sources 70, the
screens 42, the seeds 26, the breakpoints 17, the attempt's and the loop's bookkeeping ~140.  Planned ~1.7x; the line
1.4x, interleaved against HEAD in separate processes.

**Built.**  `csrc/loop.c`'s stage family (`family`: loop.h's `loop_family`), beside the multistep one: `stage_take` (the
stage times from the step's end, the step's nodes in merge order, the start's record -- `stage_find` in the last made
attempt's set, the entry's memo, or `stage_pass`: the core's passes as `_tran_core.passes` makes them, counted --
pred.c's seeds, `stage_u`'s sources with the step memo's dedupe, radau.c into the set the last made attempt's is not
in, the screens, the end's sources, `stage_estimate`: F1, numpy's dgemv, the real factor made over the full width and
reduced, numpy's dgesv on it), `stage_judge` (lte.c's new LTE_FORMULA_STAGE: `reference_stage`, the quotients,
`stage_norm` -- the kept unknowns, the mask, the squares, numpy's pairwise sum, the RMS and the worst; then the law,
`stage_retry` with the retry order, the budget), `stage_made` (the generations, the counts), `stage_accept` (the
rings, the note); `csrc/pcx_stagestep.h` (the nodes' order, the seeds, the screen, the note) shared with stagewalk.c,
whose walk is bit-identical on it.  `_tran_stageloop_c`: the entry (on the C stage Newton's last serve at the watch
counter as it stands), the checks (the C loop's of the run; the family, the pieces genuine and unshadowed, the
transform and PCNR off, numpy's underflow ignored, the core's passes as stamped, numpy's sum the C's), the bind (the
C loop's buffers and units, the radau record, an lte.c context, the stage buffers; two output sets: the radau
record's own, which no record views, and one of its own), the load, the hold after a converged Newton, and the
write-back: the loop's and the statistics' as the C loop's, `_finish_stage_step`'s outputs, `_rk_est`, the pending
nodes of a rejected attempt, the pattern mark, the screen's `_C_cache`, the memo's two generations (copies; the older
rebuilt by the core calls its Newton made, or the entry's memo after one attempt), the family's `running`, `_retry`,
`_q_carried`, the stage KLU's numeric and counts, and every count the Python's calls make.  `_tran_loop_c.attach` and
`run` hand a stage-family loop to it; `_tran_radau_c`: `served` (the stamp of the last serve), `hold`'s `outs`.

**Found.**  (1) lte.c's flags test took an underflow numpy ignores as a hand-back: 9 of a pulsed chain's 46 attempts
handed back for nothing, now `under_ok` as pred.c's.  (2) An entry right after a stop that followed made attempts made
the stopped attempt again: on the stage KLU from the numeric its first Newton moved, so a held result could have been
another's; the stage loop now leaves that attempt to the Python before its next entry.  The multistep loop's entry
does the same (an attempt's work spent twice at such a stop; its Newton holds nothing, so no answer moves) -- not
changed here.  (3) Two output sets of the loop's own, each `MAXA` 3 (2 n^2 + 3 n) doubles, would have added ~770 MB
at n = 1000: the loop borrows the radau record's own set and keeps its one within `_tran_radau_c.OUT_BYTES` (48 MiB,
n <= 362; `size` above it, where the estimate's dense O(n^3) solve outweighs what the loop saves).

**Measured** (interleaved against ba0bcd9a, separate processes, 7 rounds and 3): the Howland testbench's plain Radau
transient 0.384 -> 0.196 s (1.96x) over 4 periods, 1.835 -> 0.875 s (2.10x) over 20; the waveform's and the
statistics' digests the same, 397 of 399 attempts in C.

**Tests.**  test_stageloop_c (new, 29): 16 runs the Python's to the byte (the state, the family, the stage KLU's
counts); an entry filling its outputs; the stops (a screen that fires, the estimate's solve and the sources failing
in the C -- before and after the Newton -- the budget spent: each after-Newton stop held and taken; a Newton handed
back); the declines; the decline table; alive: the law's safety factor and the estimate's eigenvalue a ulp off.
test_radau_c: the serve's memo-views test with the stage loop off (its write-back copies).  test_instance_overrides:
the stage loop among the twins switched off.  test_runtime_map: the glue among the loop units' (no function over 80
lines).  test_pinned_pairs: 'the stage loop and its C' (new); re-recorded on their twins, each twin's tests run
first: the stepping loop, the multistep walk (loop.c), the stage walk (pcx_stagestep.h, `hold`), the error test
(lte.c's STAGE judgement), the Newton solve (lte.h's `under_ok` in NEWTON_C), the coupled stage Newton (`served`,
`hold`).  doc/transient_runtime.md: its rows; doc/transient_history.md: its entry.

**Gate** (the main tree): the suite, 5 932 passed (6 skipped, 3 xfailed); the 27 gold decks, gear2 transient and
standalone DC, against 205's: every digest, step, rejection, Newton count and path count identical (DC but for its
timing); the 27 under Radau against HEAD (ba0bcd9a, run from its worktree): every digest, step, rejection and path
count but the loop's own identical.  A first HEAD baseline, run from the main tree, had imported the main tree --
`benchmarks/spice_suite.py` puts its own root first on `sys.path` -- and was caught by its `stageloop_c:` counts.
Under Radau the loop serves 5 decks (one run each, not interleaved: 4049osc 3.29 -> 0.33 s, mosrect 0.415 -> 0.098,
gm2 0.533 -> 0.148, rca 0.158 -> 0.073, schmitecl 0.091 -> 0.041); the other 22 run the cost transform ('auto': their
compiled Jacobians expensive), which it leaves to the Python.  The sanitizer's fast tier, 5 703 passed, and its deep
tier on the loop tests, 111, clean.

## 2026-10-24 — 208. A plain Radau transient's attempts on the cost transform in C -- the C stepping loop's transform mode: 18 of the 22 transform gold decks served, 2.0x over them, every byte the same

**Why**  Item 4 measured (Andreas chose B: keep 'auto' and its answers, build the loop): 'auto' sends 22 of the 27 gold
decks to the cost transform (a compiled Jacobian of 10 KB or more; one SPICE diode is 20 KB), where 207's loop does not
serve; on them the C transform Newton was 8-39 % of an attempt (gm6 43 of ~510 us, ring 312 of ~1 080, mux8 575 of
~1 880, voter25 1 005 of ~2 610, jge 1 259 of ~4 000), the rest Python and numpy.  Planned ~1.85x over the 18 dense
decks; the line 1.4x on mux8 AND voter25, against HEAD in separate processes, default settings.

**Built.**  `csrc/loop.c`'s stage family in a transform mode (`tmode`): `tc_take` (`tc_start`: x_n's charge before the
memo's roll, the roll, x_n's C and G; `tc_frozen`: `tc_scalars` -- P's rows from numpy's M = diag(lam) Tinv by numpy's
complex division by h, Smith's form (20 000 samples the same), s0 and c1 -- radau_fold.c, the LU buffer column-major;
`tc_repattern` where the complex nonzeros leave the pattern: the CSC gathered as `_csc_of_dense`, `klu_analyze`,
`klu_z_factor`; the seeds, the stage sources; `tc_newton`: radau_tc.c on the loop's own radau_tf_t -- armed, the
probe, `tc_check` (the early abort's estimate, f0 once a point) and the resume; on a new pattern the Python loop's
sequence, a fresh factor then a refactor; `tc_end`; `tc_screen`, C alone by a core pass at the first two stages;
`tc_estimate`, dgetrs on the Newton's LU), `tc_abort` (`reject_early`'s law, the attempt's counts), `tc_undo` (the
memo's roll, f0, a new pattern's analysis taken back from an attempt not made), the memo by state and passes
(`tg_*`, three slots).  `_tran_tcloop_c`: the entry (`_tran_radau_tc`'s record now stamps `served` when its checks
ran, converged or not), the checks (the run, the family with the transform on, its pieces genuine and unshadowed --
the frozen factors' chain, the C Newton's, the step's --, the switches, the analysis solver dense, the complex solver
analysed on the last pattern, the early abort's threshold one the C arms alike), the bind (the stage loop's units by
`_tran_stageloop_c.bind_common`; its own radau_tf_t, radau_fold_t and buffers), the load, the settle (the complex
solver's refactors, a new pattern's analyses, factors, handles and record), the write-back (the stage outputs, the memo
made again, `_radau_Pc`, `early_f0`) and the counts.  `_tran_loop_c.run` sends a stage family whose transform served
to it; `LoopState.tkept`: its own kept decline (207's 'transform' decline had kept it out of nearly every entry).
pcx_stagestep.h: a seed with no prediction counted as `_pred_stages` counts it (with two nodes or more).

**Found.**  (1) A twin run's watch epoch moved 28 times inside one attempt where the garbage collector freed the
earlier runs' watched dicts (`linearsolver.py:214`, a ctypes call's allocations): each twin run collects first.  (2)
`ComplexKLUSolver` never frees a KLU analysis or factor: each new pattern leaks the old two (~1 150 each a voter25
run) -- the C mirrors it; a fix of its own (the Python and the C together), no answer moves.  (3) 207's shared kept
decline blocked the transform's entries.

**Measured** (the 27 gold decks under Radau, HEAD 8bfd4ccf from its worktree and this tree, separate processes,
interleaved per deck, one run each, default settings): over the 18 served 164.6 -> 81.2 s (2.03x); mux8 35.88 ->
16.89 s (2.12x), voter25 25.88 -> 17.02 (1.52x) -- over the line on both --, jge 79.94 -> 38.95 (2.05x), ring 6.82 ->
2.39 (2.86x), ab_integ 6.17 -> 1.64 (3.76x), gm1 2.38 -> 1.43, toronto 1.26 -> 0.55, nand 1.19 -> 0.38, gm6 0.30 ->
0.06, mike2 0.51 -> 0.10; the 27, 236.6 -> 151.8 s.  What the C leaves to the Python on voter25: 253 of its 9 742
attempts -- the 121 it stops at (102 where the transform does not converge, its dense fallback the Python's: the
plan's third step, not built; 19 the C Newton's other stops), those after them until the C Newton serves again, the
run's first --; outside the C and its glue ~3 s of the 17, the setup and the operating point among them.

**Tests.**  test_tcloop_c (new, 23): eleven runs the Python's to the byte (the state, `_radau_Pc`, the family, the
complex solver's counts); the early abort and new patterns made in C; an entry filling its outputs; the stops (a
screen that fires, the real solves and the sources failing in the C, a transform not converging, a stop after a new
pattern, the budget spent); the declines; the decline table; alive: P's first entry and the real factor's eigenvalue
a ulp off.  test_radau_tc_end, test_radau_fold: their switch-off comparisons leave out the loop's own counts (its
declines and sources table); test_instance_overrides: the transform loop among the twins switched off;
test_runtime_map: its glue among the loop units'.  test_pinned_pairs: 'the transform loop and its C' (new); the
stepping loop, the stage loop, the stage and multistep walks, the transform Newton and the end passes re-recorded on
their twins, each twin's tests run first.  doc/transient_runtime.md: its rows; doc/transient_history.md: its entry.

**Gate** (the main tree): the suite, 5 957 passed (6 skipped, 3 xfailed); the 27 gold decks, gear2 transient and
standalone DC, against 207's: every digest, step, rejection, Newton count and path count identical (DC but for its
timing); the 27 under Radau against HEAD (8bfd4ccf, run from its worktree; a tell: HEAD made no `tcloop_c:` count):
every digest, step, rejection and path count but the loops' own identical, and the A/B again, interleaved per deck:
the 18 served 159.7 -> 80.1 s (2.00x), mux8 34.06 -> 16.70 s (2.04x), voter25 25.18 -> 16.94 (1.49x).  The
sanitizer's fast tier, 5 718 passed, and its deep tier on the loop tests, 169, clean.  A first gate's suite had failed
two rules -- an argument never read (`_tc_ready`'s), a type declared in loop.c (the analysis's pointer type, now
loop.h's: in pcx_abi.h it would have moved 21 pinned pairs) -- and ran again.

## 2026-10-24 — 209. The complex KLU solver frees what it replaces -- `ComplexKLUSolver` and the C stepping loop's transform mode: voter25's peak memory 448 -> 390 MB, every byte the same

**Why**  Found building 208 (Andreas: fix it): `ComplexKLUSolver.solve_prepared` overwrote its analysis and factor at a
new pattern, and its factor at a residual fallback, without freeing either -- a voter25 Radau run left 1 152 of each
(a new pattern every eight attempts), mux8 912 patterns and 94 fallbacks -- and the C stepping loop's transform mode
(208) mirrored it.

**Built.**  `ComplexKLUSolver._free_numeric` / `_free_symbolic` (`klu_z_free_numeric`, `klu_free_symbolic`, declared
with its other entries): at a new key the factor and the analysis it holds are freed before the new analysis; at a
residual fallback the factor before the fresh one; at the solver's end (`__del__`) both.  `analyses_freed` and
`factors_freed` count the first two.  loop.c's `tc_kept`: a made attempt's new pattern (an attempt judged, or abandoned
at the early abort) frees the analysis and factor it replaced; an attempt handed back frees its own and puts the old
back (`tc_undo`), as before.  Its counts (`tc_symfreed`, `tc_numfreed`) are settled into the solver by
`_tran_tcloop_c._settle`, which now runs however the call ends.  The holders of the handles compare them with the
solver's before each use (`_tran_radau_tc`'s key, the loop's binding key): a freed handle is never used, and a new
object at a freed address is the solver's own.

**Not changed (found).**  `KLUSolver.solve` (real) overwrites its analysis and factor the same way at a new pattern
and a residual fallback.  Its handles are borrowed by four C units (`_tran_newton_c`, `_dc_newton_c`, `_pcnr_c`,
`_tran_lte_c`), each its own audit before a free can be added: an item of its own.

**Measured.**  voter25 under Radau, each tree from its own directory: peak RSS 448 -> 390 MB (the run's growth 154 ->
95 MB), 1 152 analyses and 1 152 factors freed, one of each left; 17.0 s both.

**Tests.**  test_linearsolver: the complex solver frees the handle it held at a new pattern and at a residual
fallback, one of each live, both at its end (a library wrapper noting each free).  test_tcloop_c: the twins compare
the frees as well (the C frees what the Python frees); a run with new patterns made in C leaves one analysis and one
factor.  test_pinned_pairs: 'the transform Newton' on its reference, 'the transform loop' on both sides, 'the transform
solve and its frozen form', the stepping loop, the stage loop and the multistep walk on their twin, re-recorded, each
twin's tests run first.  doc/transient_runtime.md: the transform loop's row; doc/transient_history.md:
its entry.

**Gate** (the main tree): the suite, 5 958 passed (6 skipped, 3 xfailed); the 27 gold decks, gear2 transient and
standalone DC, against 208's, and under Radau against HEAD (28bf80ef, run from its worktree, which lacks the free: a
tell): every digest, step, rejection, Newton count and path count identical, the time the same (the 18 transform
decks 80.0 -> 79.9 s).  The sanitizer's fast tier, 5 723 passed, and its deep tier on the loop and solver tests, 254,
clean -- the frees and the loop's twins under the address sanitizer.

## 2026-10-24 — 210. The transform's dense fallback in C -- `csrc/loop.c`'s transform mode continues into `_rk_step_coupled`'s fallback in the same attempt

**Why**  Andreas: "do 2" (208's third step, sized at ~2 s on voter25).  208's loop stopped at each attempt whose
transform Newton did not converge (voter25: 121 stops, 102 of them that; slowlatch 13): the Python repeated the
attempt -- its transform Newton again (~15 ms of C), then the dense fallback -- and the attempt after it too, the
fallback having moved the watch (each stop two Python attempts: 245 of voter25's 253).  LINE: voter25 at least 1.07x
against 209 in the gold A/B, no deck slower, every byte the same.

**Built.**  loop.c: `tc_newton` marks a Newton at `maxiter` (`fb`; `noconv` where the C Newton's own, `_outcome`'s
count, not the Python loop's on a new pattern, `pyloop`); `tc_body` then runs `tc_fallback` -- `tc_fb_newton`: the
memo's second roll (`tg_roll2`, the fourth slot: the attempt's start's two stand, a stop restores both), x_n's charge
as the coupled context reads it (the memo's, else a 'q' pass, not put), the seeds again (`_pred_stages`, counted
again), the stage sources (the step's memo), radau.c on the transient's record into the coupled set the last made
fallback's is not in (`gen0`), its assemblies' records put; then the screen on the last assembly's C's, the sources,
x_n's record as the estimate reads it and `stage_estimate` (207's dense estimate).  A made fallback: `tc_fb_made`
(the coupled set, its stages for `_C_pattern`, the stage KLU's counts), x_{n+1}'s record in the coupled set (`fbx`,
`fbv`).  `_tran_tcloop_c`: `_fb_ready` (the transient's `_radau_c` made by the Python's first fallback, on this core,
size and reference row, `_tran_radau_c._ready` passing -- nothing remade -- the slots and the loop's set within
OUT_BYTES), `_bind_fb` (the radau plug, the two coupled sets, the dense estimate's work), `_load_fb` (the record's run
settings as `solve` sets them, the stage KLU's numeric), `_fb_settle` (207's: the numeric and counts, the hold on a
stop after the Newton), `_fb_back` (`_radau_transform_fallbacks`, `_C_pattern`: the memo's records' C where it holds
them), `_count_fb` (the fallback, the C Newton's serves and assemblies, the marked screens, the fresh estimates).

**Tests.**  test_tcloop_c (+4, 27): fallbacks made in C on a fixed grid and on an adaptive one with the early abort
armed (the pulsed chain, three Newton iterations), the Python's to the byte with `_radau_transform_fallbacks`; a stop
after a fallback's converged Newton (the dense estimate's dgesv said failed) held and served; the transform and the
dense Newton both not converging (two iterations), the Python's; alive: the fallback record's A[0][0] a ulp off.  The
scratch twins on the gold decks: gm6 1, slowlatch 17 fallbacks made in C, every state and count the same.
test_pinned_pairs: 'the transform loop and its C' gains the coupled step's references and the fallback's glue; it, the
stepping loop, the stage loop and the multistep walk re-recorded, each twin's tests run first.
doc/transient_runtime.md: the transform loop's rows; doc/transient_history.md: its entry.

**Measured** (the 27 gold decks under Radau, HEAD 9eb458c2 from its worktree and this tree, separate processes,
interleaved per deck): voter25 16.84 -> 14.62 s (1.15x; the line 1.07x), slowlatch 0.30 -> 0.23, latch 0.23 -> 0.17,
todd3 0.19 -> 0.14 (about 1.3x each), the 18 served 79.7 -> 77.3 s, the rest unchanged (gm3, mike2 and g1310 read
0.92-0.98x on one run, 1.04-1.05x on five interleaved runs each).  Fallbacks made in C: voter25 101, latch 18,
slowlatch 17, todd3 13, jge 1 -- each deck's first the Python's (it makes the record).

**Gate** (the main tree): the suite, 5 962 passed (6 skipped, 3 xfailed); the 27 gold decks, gear2 transient and
standalone DC, against 209's, and under Radau against HEAD (a tell: HEAD made no `tcloop_c:fallbacks`): every digest,
step, rejection, Newton count and path count but the loops' own identical.  The sanitizer's fast tier, 5 728 passed,
and its deep tier on the loop, Radau and solver tests, 285, clean -- the fallback's tests under the address sanitizer.

## 2026-10-24 — 211. The transform's fallback in C left the dense Newton's record on the loop's buffers -- a fallback the Python made after the C's read the C's last attempt's charge, sources and seeds

**Found**  210's fallback twin, run alone (the first transient of a fresh process), counted otherwise: one fallback's
dense Newton not converging where the Python-only run's did; the waveform the same, the Python's own Newton then
converging.  210's glue let the C point the transient's `_radau_c` record's inputs and outputs (`qn`, `u`, `seed`,
`Ya` ...) at the loop's buffers and left them there; `_tran_radau_c.solve` writes its own buffers, set once at the
record's making, and the C read the loop's.  A Python fallback after the C's -- here where a cold process split the
attempts otherwise -- read the C's last attempt's charge, sources and seeds; a Newton converging on them would have
handed back another step's stages.  207's glue saves and restores them around its call; 210's did not.  The gold decks
never met it: past a deck's first fallback every one was the C's.

**Fixed**  `_tran_tcloop_c._call` saves the record's pointers before the C and restores them after, as
`_tran_stageloop_c._call` does.

**Tests**  test_tcloop_c (+2, 29): after the C's fallbacks the record's `qn`, `u` and `seed` its own buffers; the
fixed-grid fallback twin in a fresh process (the C's run its first transient) the Python's -- both failing on 210's
glue.  test_pinned_pairs: 'the transform loop and its C' re-recorded on its twin.

**Gate** (one for 211, 212 and 213; the main tree): the suite, 5 969 passed (6 skipped, 3 xfailed); the 27 gold
decks, gear2 transient and standalone DC, against 210's, and under Radau against HEAD (06c34cd0, run from its
worktree): every digest, step, rejection, Newton count and path count identical, the time the same.  The sanitizer's
fast tier, 5 723 passed, and its deep tier on the loop, Radau, solver and BLAS tests, 297, clean -- the new tests
among them.

## 2026-10-24 — 212. The small-signal analyses hold BLAS at one thread -- AC, Noise and Transimpedance (TwoPort and FeedbackLoop through them): 1.19-1.48x from 171 to 295 unknowns

**Why**  Andreas chose to hold AC to one thread as transient, PSS and DC are held.  `_blas.py` said threads never
lost on AC (0.80-0.98x); measured 2026-10-22 on three converters (119-164 unknowns) one thread was 1.4-1.7x, and a
peer's parallel workers took ~1 s a frequency threaded against 0.35-4.5 ms on one.

**Built**  `_blas.blas_held`, a decorator holding `_single_threaded_blas` for the call, on `AC.solve`, `Noise.solve`
and `TransimpedanceAnalysis.solve`; `TwoPortAnalysis` and `FeedbackLoopAnalysis` reach it through them.  The periodic
small-signal analyses after a PSS (`PAC`, `pnoise`, the spectra) are not held: not measured past the threading size.

**Measured** (AC over 61 frequencies at eight gold decks' operating points, the first voltage source driven, one
process alternating three rounds, a quiet box): ring (88 unknowns) 1.00x, the answers bit for bit; mux8 (171) 1.30x,
voter25 (200) 1.32x, jge (244) 1.19x, rich3 (269) 1.32x, arom (295) 1.48x, gm19 (429) 1.00x, opampal (519) 0.95x.
The answers move at rounding level -- 4e-16 to 5e-11 relative, rich3's 5e-9 (an ill-conditioned system).  The hold
is unconditional above 400 unknowns too: the answers the same whatever the core count, and no oversubscription in a
parallel sweep, against 5 % at 519.

**Tests**  test_blas_limit (+3): AC, Noise and Transimpedance read one BLAS thread inside, the pool opened to four
outside -- each failing without the hold.

**Gate** (one for 211, 212 and 213; the main tree): the suite, 5 969 passed (6 skipped, 3 xfailed); the 27 gold
decks, gear2 transient and standalone DC, against 210's, and under Radau against HEAD (06c34cd0, run from its
worktree): every digest, step, rejection, Newton count and path count identical, the time the same.  The sanitizer's
fast tier, 5 723 passed, and its deep tier on the loop, Radau, solver and BLAS tests, 297, clean -- the new tests
among them.

## 2026-10-24 — 213. The real KLU solvers free what they replace and what they hold at their end -- `KLUSolver` at a new pattern and a residual fallback, every one at its end: 2 000 pattern changes +115 -> +9 MB

**Why**  Andreas: fix the real KLU solver's leak (found at 209).  `KLUSolver.solve` overwrote its analysis and factor
at a new pattern, and its factor at a residual fallback, without freeing them; and no real KLU solver freed anything
at its end -- `solve`'s analysis and factor, `factor`'s analyses -- so each run's `StructuralKLU` (one to three a gold
deck) left its analysis and its last numeric.

**Audited**  The C units that borrow a `StructuralKLU`'s handles (`_tran_newton_c`, `_dc_newton_c`, `_pcnr_c`,
`_tran_lte_c`, `_tran_radau_c`) hold the solver itself and read its numeric before each call; pcx_kluseq.h frees the
numeric a sequence replaces, on a serve and on a hand-back, and newton.c frees its chord's; `StructuralKLU.solve`
frees what its fallbacks replace.  The leaks: the plain solver's two overwrites, and every solver's end.

**Built**  `KLUSolver._free_numeric` / `_free_symbolic` (`klu_free_symbolic` declared beside `klu_free_numeric`),
`analyses_freed`, `factors_freed`: `solve` frees the factor and the analysis before a new pattern's analysis, the
factor before a residual fallback's; `__del__` frees `solve`'s factor and analysis and `factor`'s analyses -- none of
its factorisations is left then, each `_KLUNumeric` holding its solver.  `StructuralKLU` inherits the end.  `factor`'s
analyses stay one a pattern while the solver lives (its cache).

**Measured** (400 unknowns, each tree from its own directory): the plain solver over 2 000 alternating patterns, peak
RSS +115 -> +9 MB (1 999 analyses freed); 300 structural solvers made and dropped, +20 -> +4 MB.

**Tests**  test_linearsolver (+2): the plain solver frees the handle it held at a new pattern and at a residual
fallback, at its end the last of each and `factor`'s analysis; a structural solver frees its analysis and its last
numeric at its end.

**Gate** (one for 211, 212 and 213; the main tree): the suite, 5 969 passed (6 skipped, 3 xfailed); the 27 gold
decks, gear2 transient and standalone DC, against 210's, and under Radau against HEAD (06c34cd0, run from its
worktree): every digest, step, rejection, Newton count and path count identical, the time the same.  The sanitizer's
fast tier, 5 723 passed, and its deep tier on the loop, Radau, solver and BLAS tests, 297, clean -- the new tests
among them.

## 2026-10-24 — 214. The BLAS hold sized for AC and extended to PAC -- AC's family holds below 440 unknowns (from there the threads pay: 0.46x on one thread at 1 503), every `PAC` surface holds (1.1-1.8x)

**Why**  Andreas: take on the two items 212 left open -- AC's hold above ~400 unknowns (one thread 5 % slower at
opampal's 519), and PAC and pnoise, not measured past the threading size.

**Measured**  AC over 61 frequencies on an RC ladder (the dense solve's cost the circuit's size, whatever its kind),
one thread against OpenBLAS's 24, three rounds each, alternating: 303 unknowns 1.23x, 403 1.27x, 423 1.11-1.13x;
443 0.97x, 463 0.87-0.90x, 483 0.94x, 503 0.80-0.85x, 603 0.72x, 803 0.60x, 1 003 0.58x, 1 503 0.46x; the gold
decks gm19 (429) 1.00x, opampal (519) 0.95x.  With the box busy (a desktop's bursts) two threaded runs, at 443 and
483, took 2.4-4x the one-thread time.  PAC on a driven noisy RC ladder with a diode (PSS of 60 points; `solve` at
three frequencies, `pnoise` at one; five rounds): under gear, 53 to 803 unknowns, `solve` 0.92-1.76x and `pnoise`
1.10-1.28x on one thread (a threaded `pnoise` at 103 once 35x slower, the box busy); under Radau, 103 to 503, `solve`
0.97-1.09x and `pnoise` 1.23-1.57x.  PAC's work is many circuit-size operations through the period, not one large
dense solve a frequency: its threads never paid.

**Built**  `_blas.AC_THREADS_FROM` (440) and `blas_held_below(n)`: AC, Noise and Transimpedance hold below it and
leave BLAS its threads from it.  `_blas.blas_held_methods`, a class decorator holding every public method of a class,
its own and its bases', on `PAC`: 30 surfaces (`solve`, `pnoise`, the spectra, the covariances, the adjoint rows,
the sampled noise).  Source introspection is unchanged (`inspect` follows `__wrapped__`): no pinned pair moved.

**Tests**  test_blas_limit (+5, 13): AC, Noise and Transimpedance hold below the size and leave the pool alone on a
450-section ladder; PAC's `solve` and `pnoise` read one thread inside -- each failing without the hold.

**Gate** (the main tree): the suite, 5 974 passed (6 skipped, 3 xfailed).  No C, and no transient, DC or solver path
changed (the gold decks run neither AC nor PAC): the gold decks and the sanitizer not run.
