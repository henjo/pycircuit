# Can the DDD work be reused for symbolic distortion? — feasibility

**Verdict: yes, and the fit is closer than it first appears — but not for the
reason that motivated the question, and the largest available win may not need
DDD at all.**

Written before any code, as the reasoning behind
`doc/distortion_ddd_plan.md`. Everything numeric here was measured, not
argued; the spikes are reproducible from the snippets quoted.

## 1. Where this came from

`doc/distortion_mimo_plan.md` §8.3 records an open question. Stage E measured
symbolic expression **size** and found it polynomial in truncation order, with
multi-node costing a small constant factor. It measured nothing about **cost**,
and `sympy.cancel` failed to finish in 900 s at `U^7` on a two-node system.
Those are different claims and only one was established.

The suggestion to look at the DDD work is a good one because the two problems
share a structural property that is unusually specific.

## 2. The structural match

The distortion method's central claim, from the theory page: **every step is a
solve against the same linear operator `Y(s) = G + sC`, evaluated at a
different frequency.** Nothing needs a nonlinear solve or a new matrix.

The DDD work's central objects are the answer to almost exactly that:

| distortion needs | DDD provides |
|---|---|
| `det(G + sC)` at many harmonic frequencies | `s_expand` — determinant split by powers of `s`, **coefficients shared as diagrams**; one construction serves every harmonic |
| one matrix, very many right-hand sides (one per graded `(harmonic, power)` key) | `DDDFamily` — numerators expressed as **cofactors of the original matrix**, so a new right-hand side is a *recombination*, not a rebuild |
| many evaluations of overlapping structure | `eval_roots` — one pass over shared structure, memoised across roots |

The `DDDFamily` docstring states its own motivation as noise analysis, where
"one transfer function per noise source costs little more than a single
transfer function: the transimpedances share almost everything." **A graded
perturbation solve has the same shape**: dozens of right-hand sides against one
matrix, at a handful of distinct frequencies.

## 3. What was measured

Three spikes on the gm-C biquad pencil `[[-g2 + s*C1, -g4], [-g3, s*C2]]`,
symbolic entries, evaluated against `numpy.linalg.solve`.

**Complex `s` works.** This was the first thing that could have killed it: DDD
is built for a symbolic frequency variable, and harmonics need `s = j m w`.

```
m=1  max rel err 2.50e-15
m=2  max rel err 2.13e-16
m=3  max rel err 1.90e-16
```

**`s_expand` reproduces the determinant across harmonics** from one
construction: degree 2, **5 shared vertices**, `det` correct to 5e-15 at
`m = 1` and `m = 3`.

**The right-hand side really is free, but only on the fast path.** This is the
important measurement:

| route | cost per right-hand side |
|---|---|
| `DDDCombination.eval(env)` per solve | **3.05 ms** |
| `eval_roots` memo once per frequency, then recombine cofactors | **3.07 us** |

**A factor of 1000, for the same answer to 2.5e-15.** The DDD memory records
this trap as a 200x `subs` cost fixed by a `_resolve` fast path and keeping
factors unmultiplied; here the same trap is worse. **Any prototype that reaches
for `eval` in a loop will measure DDD as hopeless and be wrong.**

## 4. What DDD does *not* solve — the honest part

**The blow-up in §8.3 is not in the determinant.** After `p` passes of the
graded recurrence, a coefficient is a sum of products of up to `p` cofactor and
determinant factors *at different frequencies*. DDD compacts each factor and
shares them; it does not by itself stop the product tree growing. Anyone
expecting "use DDD and the symbolic cost problem goes away" will be
disappointed.

**DDD is not a numeric accelerator, and must not be sold as one.** Even the
fast path is ~3 us per solve against numpy's sub-microsecond dense solve on a
2x2. The existing numeric path must stay numpy. Nothing in this work should
make the current tests faster, and a claim that it does would be false.

**The 2x2 biquad proves nothing about scale.** Five vertices is trivial.
Whether sharing pays on a matrix where it matters is unmeasured — and that is
the whole question.

## 5. So what would actually help

Ordered by expected value per unit of work, cheapest first. **This ordering is
the main output of this document.**

1. **Stop calling `cancel`, and keep coefficients factored.** The single
   measured fact in §3 is that a memo plus recombination beats substitution by
   1000x. That lesson is *transferable without DDD*: represent each graded
   coefficient as a sum of products of unevaluated factors, and evaluate
   through one memo. **It is entirely possible this alone answers §8.3.**
2. **DDD for the linear solves.** Replace the `solve` callable with a
   DDD-backed one, so every harmonic draws on one shared construction.
3. **Canonical sharing across the whole computation.** DDD's distinctive
   contribution over (1) is that identical minors recurring across harmonics
   *and* across perturbation orders collapse to one vertex. sympy has no
   canonical merge. This is the part that could beat a plain factored
   representation, and the part nobody has evidence for.

**Reconsider-if for the whole exercise:** if stage A of the plan shows that a
plain factored graded ring already makes `U^7` tractable symbolically, then
DDD is an optimisation of an already-solved problem and this work should be
stopped and recorded as such. That outcome is a success, not a failure.

## 6. Why this is a good place to try it

The DDD work's own calibration boundary is a real limitation: for named
op-amps "a published `|DDD|` pins an order of magnitude, not a number", because
formulation and device model differ between papers. Its *exact* checks were
confined to full matrices.

The distortion implementation is a much better oracle. It is gated against
**four published closed forms** across two circuits (eqs. 41, 42, 48, 52), an
**independent ODE integration** (magnitude, phase to 0.0000 deg, and full
waveform to 2e-10), and a **2-D numerical Fourier extraction** for the
two-tone exponential. A DDD-backed path can be required to reproduce all of it
to floating point.

**That is the strongest argument for doing this at all**: not that DDD will be
faster, but that for once there is an oracle good enough to prove a symbolic
representation correct rather than merely plausible.

## 7. Constraints and traps carried over

- **Use `eval_roots`, never `eval`, in any loop.** §3. This is not a
  micro-optimisation; it is the difference between feasible and not.
- **`DDDFamily` defaults to `'min-degree'`, not `'auto'`**, deliberately:
  `'auto'` reorders the matrix and the class addresses cofactors by *original*
  row and column. Do not "improve" this to `'auto'`.
- **The exact-rational blow-up that killed GiNaC cannot arise here** (DDD
  negative result P3): a diagram never multiplies entries symbolically. With
  complex float entries we are further from it still.
- **Keep factors unmultiplied as tuples**, as `DDDCombination` does — a
  premultiplied product is a sympy expression and costs a full substitution.
- **Separate files and a separate label throughout**, per the instruction that
  prompted this: the existing implementation must remain the reference, not a
  thing to be refactored around.

## 8. Which diagram for which regime — the framing that matters

Added after the first draft, prompted by two things: the question of whether a
*determinant* diagram is even the right shape, and the goal of running
distortion analysis on **bigger circuits**. The second inverts part of §4.

The cost has two independent sources, and they dominate in opposite regimes:

| | small circuit, high truncation order | big circuit, low order |
|---|---|---|
| dominant cost | the **product tree** of graded coefficients | the **determinant and cofactors** |
| right tool | a canonical *polynomial* diagram, or just factored evaluation | **DDD, exactly as built** |
| evidence | tree/DAG ratio 1.6 → 4.7 over `U^3`..`U^9`; `cancel` unfinished at `U^7` on 2 nodes | µA741 flat `\|DDD\|` 1040, hierarchy 1040 → 156, noise 11 088 → 26 vertices |

**§4 said "the blow-up is not in the determinant". That is true of the 2x2
biquad and false of a 25x25 op-amp.** On the biquad the determinant has five
vertices and the products dominate; on a real amplifier the determinant is the
millions-of-terms object DDD was built for, and the graded products sit *on
top of* it. For the stated goal — bigger circuits — DDD is not the wrong
shape. It is the necessary substrate, and the polynomial growth rides above it.

The two compose rather than compete: a graded coefficient on a large circuit
is a product of cofactors, each of which DDD represents compactly.

## 9. Other diagram variants, and an honest note on sourcing

**These are directions to evaluate, not results.** The project rule applies:
*only use external data if both the thing measured and the measurement are
understood*. None of the literature below has been re-read for this document
and **no number from it should be quoted until the paper is in hand**. What
follows is a shape-matching argument, which is all it claims to be.

### 9.1 The natural extension of what already exists — a *graded* diagram

`SExpandedDDD` is "one root per power of `s`, sharing subgraphs, memo key
`(rows, cols, power)`". The distortion work grades by `(harmonic, power of the
drive)`. So the direct analogue is **one root per graded key, memo key
`(rows, cols, harmonic, power)`** — the same multiroot construction Shi & Tan
use for `s`, applied to the grading this analysis already carries.

This is the cheapest variant to try because the machinery, the sharing
argument and the trimming of zero coefficients all exist; only the expansion
variable changes. **It is also the one most likely to pay**, because the
grading is exactly where the repeated structure lives.

### 9.2 Diagrams built for polynomials rather than determinants

A determinant is one specific polynomial — a signed sum over permutations —
and DDD's variable order and zero-suppression are tuned to it. The graded
numerator is a *general* multivariate polynomial. Two families were built for
that case:

- **Taylor Expansion Diagrams (TED)** — canonical form for multivariate
  polynomials by Taylor decomposition, `f = f(0) + x f'(0) + x^2 f''(0)/2 +
  ...`, each child a diagram. Developed for datapath verification, where BDDs
  blow up on arithmetic. The decomposition is the same one this analysis
  already performs on its devices, which is at least suggestive.
- **Multiplicative binary moment diagrams (`*BMD`, `K*BMD`)** — moment
  decomposition `f = f_0 + x f_1` with **multiplicative edge weights**. Those
  weights are the canonical version of the trick that gave 1000x in §3:
  keeping a common factor out of the representation instead of distributing
  it.

**Reconsider-if, and it is a real one:** both are canonical forms for
polynomials over a *fixed* variable set. Here each harmonic contributes its
own set of transfer-function values, so the variable count grows with the
truncation order. Whether either stays compact under that is exactly the
question, and it is unmeasured.

### 9.3 Considered and set aside

- **ZDD over factor multisets.** Terms are *multisets* — a factor can appear
  squared — and ZDDs are canonical for *sets*. Encoding multiplicity as
  distinct variables is possible but reintroduces the growth it was meant to
  remove. **Reconsider if** a formulation appears where every factor is
  square-free.
- **ADD / MTBDD.** Terminals carrying values help when many coefficients
  repeat exactly. Here they are distinct rational functions.

### 9.4 The lever for *bigger circuits*, which needs no new diagram at all

`HierarchicalDDD` and `suppression_order` already exist and already work:
1040 → 156 vertices on the µA741, 11 088 → 26 and 86 s → 0.48 s on noise.
Suppressing internal nodes before expanding is a property of **the matrix**,
and the distortion path solves the same matrix at every harmonic.

**So the most promising route to larger circuits may not be a new diagram
variant but the existing hierarchy applied to the existing solve** — cheaper
than anything in §9.2 and already validated on a circuit far larger than
anything the distortion work has run. It goes in the plan ahead of the exotic
options for that reason.
