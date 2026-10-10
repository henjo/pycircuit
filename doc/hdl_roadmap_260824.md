
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
