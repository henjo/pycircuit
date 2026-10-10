
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
