# The transient's C runtime: a map

The transient steps in Python (`transient.py`: `Transient`, `_SteppingLoop` and the step
families). Where it can, a C unit makes the same arithmetic, **bit for bit**, and the Python
continues as if it had made it itself. This page says which unit stands for which Python, what
each one's glue module owns, how a C path decides it may run, and how each one is proven. Read it
before any `csrc/*.c` or `_tran_*_c.py` file. History (why a unit exists, what it measured when it
was built) is in `doc/transient_history.md` and the roadmap, not here.

## The one rule

**Exact or not at all.** Every C path reproduces its Python reference to the last bit -- the same
floating-point operations in the same order, numpy's own BLAS and LAPACK routines (their
addresses, not a copy), the same libm, Python's `min`/`max` tie rules -- or it does not run.
The Python is both the reference and the fallback:

* a C path that cannot reproduce the case **declines** before it starts: the Python runs, and the
  decline is counted in `_paths` as `<prefix>:<reason>`;
* a C path that meets something it cannot reproduce part-way (a non-finite value, a
  floating-point flag numpy would report, a case its Python sends elsewhere) **hands back**: it
  writes nothing, and the Python redoes the work from the same state.

Every C unit can be switched off with its environment variable set to `0` (below). With all of
them off the run is the pure Python one; with any of them on, the answer is the same bytes.

## The layers

From the outside in. "Called from" is where the Python reaches the glue.

| Layer | Python reference | C unit | Glue module | Switch | `_paths` | Called from |
|---|---|---|---|---|---|---|
| Stepping loop | `_SteppingLoop.where / take / judge / accept` with a multistep family | `csrc/loop.c` | `_tran_loop_c` | `PYCIRCUIT_LOOP_C` | `loop_c:` | `_SteppingLoop.execute` |
| Controller law | `IntegralController._evaluate` (past its error test) | `csrc/ctrl.c` | `_tran_ctrl_c` | `PYCIRCUIT_CTRL_C` | -- | `loop.c` |
| Error test | `StepController._charge_lte`, `_normalised`, the maximum | `csrc/lte.c` | `_tran_lte_c` | `PYCIRCUIT_LTE_C` | `lte_c:` | `StepController._max_error`, `loop.c` |
| Stage predictor | `_StagePredictor._predict_fast` | `csrc/pred.c` | `_tran_pred_c` | `PYCIRCUIT_PRED_C` | `pred_c:` | `loop.c`, `_StagePredictor._pred_stages` (`stages`), `stagewalk.c` |
| Newton solve | `StandardNewton`, `ChordNewton`, `jacobian_only` | `csrc/newton.c` | `_tran_newton_c` | `PYCIRCUIT_NEWTON_C` | `newton_c:` | `Transient._solve_timestep`, `loop.c` |
| DC base solve | `StandardNewton.solve_system` on the DC operating point's compositions (`dcanalysis`'s residual, source stepping, gshunt, junction gmin: `_dc_newton_c.composition`) | `csrc/dcnewton.c` | `_dc_newton_c` | `PYCIRCUIT_DCNEWTON_C` | `dcnewton_c:` | `StandardNewton.solve_system` (the DC's ladders) |
| PCNR DC solve | `pcnr.solve_dc`'s loop from its start: the view's assembly (`pcnr_view`), each participant at its own unknowns (`PcnrDevice._kernel_stamp`), the Schur step, the solve, `refine` (each class's refine kernel: its law parameters, `limit.c`'s laws, `limit_block`'s Kruskal), the three tests; a NaN carried as numpy carries it to the non-finite iterate's RuntimeError (`pcnr.NONFINITE`) | `csrc/pcnrdc.c`; per class `csrc/pcnr_refine.c` after `csrc/limit.c` | `_pcnr_c` | `PYCIRCUIT_PCNR_C` | `pcnr_c:` | `pcnr.solve_dc` (`DC(pcnr=True)`'s ladders, `solve_dc_continued`) |
| Branch screen | `_BranchCheck._branch_screen_C` (the SVD stays Python) | `csrc/branch.c` | `_tran_branch_c` | `PYCIRCUIT_BRANCH_C` | `branch_c:` | `_tran_branch._branch_screen`, `loop.c` |
| Assembly | `_residual_and_jacobian`, `jacobian_only`, `residual_only` | `csrc/core.c` | `_tran_core` | `PYCIRCUIT_TRAN_CORE` (`PYCIRCUIT_CORE_PASSES`, `PYCIRCUIT_CORE_PATTERN`, `PYCIRCUIT_SKIP_UNREAD_J`) | `core.<want>:`, `core.ready:` | the transient's evaluation sites, `newton.c` |
| Sources | `_stamp_plan.source_direct` | `csrc/sources.c` | `_tran_sources` | `PYCIRCUIT_SOURCES_C` | `sources_c:` | `loop.c` (its table), `_stamp_plan.source_direct` (`_direct_c`: a table on the direct pass's program) |
| Limiter walk | the Newton's limiter walk; the hand-written limiters | `csrc/walk.c`, `csrc/limit.c`, `csrc/psp_limit.c` | `_hdl_climit` | `PYCIRCUIT_HDL_LIMIT_WALK`, `PYCIRCUIT_HDL_CLIMIT`, `PYCIRCUIT_PSP_LIMIT_C` | `walk:`, `climit:` | the Newtons, `newton.c` |
| Device passes | the circuit passes `G`, `C`, `i`, `q` of the HDL elements | generated kernels; `csrc/pass.c` (a batch: one kernel over its elements) | `_hdl_cbackend`, `_hdl_batch` | `PYCIRCUIT_HDL_BACKEND`, `PYCIRCUIT_HDL_BATCH` (`PYCIRCUIT_HDL_CFUSE`: the fused kernel a call's passes share; `PYCIRCUIT_HDL_ZERO_U`: the source passes' skip) | `batch.<pass>:`, `once:backend:` | the stamp plans' batches; `core.c` calls the kernels itself |
| Invalidation | -- (a dict watcher) | `csrc/watch.c` | `_watch` | `PYCIRCUIT_WATCH` | -- | every readiness check |
| Build and load | -- | -- | `_cbuild` | `PYCIRCUIT_C_SANITIZE`, `PYCIRCUIT_HDL_CC`, `PYCIRCUIT_HDL_CACHE` (the store; `PYCIRCUIT_HDL_CACHE_DIR` its place) | -- | every glue module |

Outside the multistep loop, the same rules: the Radau stage solve (`csrc/radau.c`, `_tran_radau_c`,
`PYCIRCUIT_RADAU_C`, `radau_c:`; on a StructuralKLU of the stage block -- a shooting walk's, or a plain run's from
250 block unknowns, `_RadauStages._run_stage_klu`, `PYCIRCUIT_RADAU_KLU`), its tensor-product form and fold (`csrc/radau_tc.c`,
`csrc/radau_fold.c`, `_tran_radau_tc`, `PYCIRCUIT_RADAU_TC`, `PYCIRCUIT_RADAU_FOLD`,
`PYCIRCUIT_RADAU_TC_END`, `radau_tc:`, the fold's `radau.frozen:`), and the shooting sensitivities
(`csrc/sens.c`, `shooting/_sens_c`, `PYCIRCUIT_SENS_C`, `pss.sens:`; on the walk's KLU
`csrc/sensklu.c`, `shooting/_sensklu_c`, `PYCIRCUIT_SENSKLU_C`, `pss.sensklu:`; the multistep walk's
`csrc/lmmsens.c`, `shooting/_lmmsens_c`, `PYCIRCUIT_LMMSENS_C`, `pss.lmm:`), and the plain stage
walk that drives the stage solve, the predictor, the branch screen and the step map on the walk's
KLU a chunk of steps a call (`csrc/stagewalk.c`, `shooting/_stagewalk_c`, `PYCIRCUIT_STAGEWALK_C`,
`pss.stagewalk:`; a stop after a converged Newton held for the Python's repeat of the step,
`_tran_radau_c.hold`, `radau_c:held`), and the plain multistep walk on the C stepping loop's own walk
mode, a chunk of steps a call -- each step `take` and `accept` on the walk's grid, its sources the Python's
in a table, then the step map (`csrc/loop.c`, `shooting/_lmmwalk_c`, `PYCIRCUIT_LMMWALK_C`, `pss.lmmwalk:`;
a column's prune and a map handed back the Python's). A few glue modules also make a
count the Python would have made on the same path: `umemo:hit` (the source memo, by the Newton and
Radau glue), `pss.sync:skipped` (the sensitivities), `chord:chosen` (the Newton).

**The linear system** has no unit of its own. Dense solves are numpy's and SciPy's LAPACK
(`dgesv`, `dgetrf`/`dgetrs`), called by address from `newton.c`, `lte.c` and `pred.c` -- the very
libraries the packages mapped, found one way (`_cbuild.openblas`, `_cbuild.address`). Sparse
solves are `linearsolver.StructuralKLU` (KLU on the circuit's structure): `newton.c`, `lte.c`,
`dcnewton.c` and `pcnrdc.c` call KLU through the run's own symbolic and numeric objects, and the glue adds to the solver
object's counts (`factors`, `refactors`, `numerics`) what the Python would; a refactor that meets
a zero pivot bails to the Python, which factors afresh and counts it (`pivot_fallbacks`). The
Newtons' rounding floor (`nrsolver._within_rounding`, in `newton.c` and `dcnewton.c`) solves with
J^T -- numpy's `dgesv` of it, or `klu_tsolve` on the numeric the step refactored -- and makes its
floors with numpy's own `ddot`.

Every unit has the same shape: a header `csrc/X.h` that is both what the C compiles and what cffi
reads (its structs, its entry `pcx_X(X_t *)`, named for the unit, its codes as named enums); a `csrc/X.c` that
includes it; and a glue module that checks the case, fills the struct, calls the entry, reads the
status, and writes back. `csrc/pcx_abi.h` holds what units share (the core's and the walk's
signatures, the LAPACK, BLAS and KLU typedefs, the pass bits and formula codes);
`csrc/pcx_dense.h` the small dense helpers;
`csrc/pcx_klu.h` the run's StructuralKLU as a unit holds it (a struct `dcnewton.h` and `pcnrdc.h`
embed: their cdef is its text and then their own) and `csrc/pcx_kluseq.h` StructuralKLU.solve's
own sequence, written once for `dcnewton.c` and `pcnrdc.c`;
`csrc/pcx_fenv.h` the floating-point flags a unit tests and THE FLAG RULE every such unit keeps
(each operation a flag test covers stores its result before the test).
`csrc/pcx_minmax.h` min and max as Python and numpy pick them (which operand a tie or a NaN
keeps is part of the bits: `pcx_py_max`, `pcx_py_min`, `pcx_np_maximum`, `pcx_np_clip`).

## Who owns what

What each glue module keeps on objects it does not define. What a run keeps on the `Transient` is
declared in `transient._RUN_STATE`, with its owner and lifetime
(`tests/test_run_state.py` holds the list to the code).

| Glue | Keeps | Writes back on a serve |
|---|---|---|
| `_tran_loop_c` | on the `_SteppingLoop`: `c`, a `LoopState` (`attach`; None where the loop never tries) -- `run` (the bound run: structs, buffers, the units' pointers), `tab` (the source table), `kept` (the watch counter of the last kept decline), `skip` / `backoff` (the spacing of entries that made no attempt) | the loop's position (`t`, `h`, `h_next`, `landing`, `order_drop`, `point_retries`, `t_break`, `J`, `imposed`), the transient's run state and statistics, the controller's `last_err`, `_q_carried`, `_ref_running`, `_retry` |
| `_tran_newton_c` | on the transient: `_newton_c` (its record: buffers, pointers, the walk's and KLU's bindings, `served`), `_newton_c_no` (a kept decline); on the limiter walk its tables (`F`, `PR`, `mirror`); on each limiting element `_hdl_cp` | the step's state: `_q_cache`, `_C_cache`, `_Cmat`, `_Geq`, `_iq`, `_companion_coeffs`, `_C_pattern`, `active_integrator`, `_effective_method`; the statistics' `newton_iterations`; the `StructuralKLU`'s `numerics` |
| `_dc_newton_c` | on the DC analysis: `_dc_newton_c` (its record: buffers, the core, the KLU and walk bindings) | the `StructuralKLU`'s `_numeric`, `refactors`, `factors`, `pivot_fallbacks`, `residual_fallbacks` |
| `_pcnr_c` | on the circuit: `_pcnr_c` (its record: the view's core, the participants' tables, buffers, the KLU binding); on each participating class's `_hdl_info`: `_c_pcnr_refine` (its refine kernel), `_c_pcnr_refine_status`, `_c_pcnr_kernels` | the `StructuralKLU`'s `_numeric` (None after a failed factor), `refactors`, `factors`, `pivot_fallbacks`, `residual_fallbacks` |
| `_tran_lte_c` | on the step controller: `_lte_c` (its context: buffers, the KLU binding, `served`) | the controller's `_ref_running` |
| `_tran_core` | on the transient: `_tran_core` (the plan, the core, their stamp), `_passes_no` (a kept decline); on each HDL element: `_hdl_cp` (its packed parameters); on the stamp plan, the `methods` it lacks | the same step state as the Newton |
| `_tran_radau_c` | on the transient: `_radau_c` (its record: buffers, the output sets the device memo takes as views, the stage KLU's binding -- a shooting walk's `_stage_klu`, or a plain run's (`_run_stage_klu`) -- the last serve's converged block, and `pre`: a converged result held for the next call, `hold`) | the stage state; on the stage KLU its `_numeric`, `refactors`, `factors`, `pivot_fallbacks`, `residual_fallbacks` |
| `_tran_radau_tc`, `_sens_c` | on the transient: `_radau_tc`, `_tc_end`; `_sens_c` | their stage state |
| `_lmmsens_c` | on the walk object (the PSS): `_lmmsens_c` (its record: the struct, buffers, the last call's arrays and the column lists they keep) | -- |
| `_sensklu_c` | on the walk's StructuralKLU: `_sensklu_c` (its context: the struct, the probe's buffers, the binding) | the walk's `_KLUWalk`: its numeric, generation, `factors`, `refactors`, `fallbacks`; the solver's `numerics` |
| `_stagewalk_c` | on the inner transient: `_stagewalk_c` (its context: the struct, the state's and a chunk's buffers, the predictor state, the screen's struct) | the stage state (`_qlast`, `_iqlast`, `_pred_hist`, `_q_cache`, `_rk_Y`, `_iq`, `_Cmat`, `_Geq`, `_C_cache`, `_branch_cmax`, `_dt` and its history), the device memo's records of the last step's stages; the stage KLU's `_numeric` and counts; the walk's `_KLUWalk` and the solver's `numerics` as `_sensklu_c`; the live columns; on a stop after a converged Newton, `_radau_c.pre` |
| `_lmmwalk_c` | on the walk object (the PSS): `_lmmwalk_c` (its context: the C loop's run bound for the walk without a sources unit, the walk's buffers -- the reduced C ring, the monodromy's, a chunk's lengths, ends and sources -- its step map's struct) | the multistep state as `_tran_loop_c` writes a stepping loop's (`_qlast`, `_iqlast`, `_pred_hist`, `_q_cache`, `_C_cache`, `_C_pattern`, `_iq`, `_Cmat`, `_Geq`, the active integrator and `_companion_coeffs`, `_dt` and its history -- the grid's own objects -- `_branch_cmax`); the C Newton's KLU counts; the walk object's `_Jf`, `_Geq`, `_C`, `_coeffs`; the live monodromy (`P`, `Pq`, `n`); the solve's source memo and its counts, by the Python's own passes |
| `_tran_pred_c` | on the transient: `_pred_c` (`stages`' struct and history buffers) | -- (the seeds, returned) |
| `_tran_ctrl_c`, `_tran_branch_c`, `_tran_sources` | nothing: their state lives in the caller's struct or table (the direct pass keeps its table on its program, `_SourceDirect.tab`) | -- |

## When a C path may run: the run context

A C path rests on facts that cannot change while nothing changes: the elements are the ones the
plan was built for, no method was patched or shadowed on an instance, the parameters are the ones
packed. Checking them costs far more than the work, so they are checked once and then watched:

* `_watch` keeps ONE counter, bumped by a C callback on any change to any watched dict (an item
  set, added or deleted; the dict cleared or freed). A checker that passed **arms** the dicts it
  read and keeps the counter's value as its stamp; while the counter stands, the check stands.
  What a dict watcher cannot see (a reassigned `__dict__`, a class attribute, a list's length) the
  checker compares on every call.
* **Served marks.** `_tran_newton_c`'s record and `_tran_lte_c`'s context set `served` to the
  counter at a call they served, and clear it on every other outcome. The C loop enters only when
  the Newton's mark (and, adaptive, the error test's) equals the counter now: the last Python
  attempt passed every check, and nothing watched has moved since. With watching off the counter
  reads -1 and the loop never enters.
* **Kept declines.** A decline whose reason is the run's, not the attempt's, is kept and not
  re-checked: the C loop's with the watch counter (`c.kept`, until something watched moves), the
  Newton's and the core's with the circuit's stamp plan (`_newton_c_no`, `_passes_no`, until the
  plan is another).
* **The backoff.** An entry that made no attempt (the screen fired on the first one) doubles the
  Python attempts before the next entry, up to `_tran_loop_c.BACKOFF_MAX`; an entry that makes any
  resets it.
* **Genuine pieces.** Each glue checks that the Python it stands for is the genuine one (not
  patched, subclassed or shadowed): `_tran_loop_c._own`, `_tran_newton_c._chain`,
  `_tran_lte_c._globals`. A test that patches one of them gets the Python, every time.

## The write-back

A served C call leaves exactly the state the Python would have left: every `_RUN_STATE`
attribute, the statistics, the controller's state.

The `_paths` counts say WHICH path ran, not what its arithmetic was: a unit counts its own serve
(`newton_c:served`), and the counts the Python would have made inside the work the C did -- the
core's passes, the limiter walk, the memos, inside a C Newton -- are not made at all. The C loop
is the one exception: it stands in for the Python loop calling the same C units, so
`_tran_loop_c._write_back` re-creates the counts their glue would have made (`newton_c:`,
`lte_c:`, `branch_c:`, the sources' and the predictor's). The A/B tests compare, between a unit on
and off: the bytes of every output, the times, the statistics, the warnings, every `_RUN_STATE`
attribute and the controller's state -- and, for the loop, every count but its own `loop_c:`.

## How each unit is proven

* **Pinned pairs** (`tests/test_pinned_pairs.py`): each Python reference and its C twin as a
  digest pair; a change to one side fails naming the other. Re-record only after checking the
  twin, in the same commit.
* **A/B tests** per unit (`test_newton_c.py`, `test_lte_c.py`, `test_loop_c.py`, `test_tran_core.py`,
  `test_tran_sources.py`, ...): byte equality on built
  circuits and drawn cases, every decline and hand-back reached, and an **alive** test that plants
  a one-ulp change in the C path and sees another run.
* **The ABI test** (`tests/test_c_abi.py`): each header is the cdef its unit is loaded through,
  each unit compiles alone without a warning, each code reads as cffi reads it, and each
  precompiled FFI declares exactly what parsing the header declares.
* **The sanitizer tier** (`scripts/sanitize_suite.sh`): the suite with every unit built under
  AddressSanitizer and UndefinedBehaviorSanitizer.
* **The gold decks**: every deck's answer the same bytes with the C on and off (any change to a
  solver or the operating point runs them all).

## Building and loading

`_cbuild.load_unit(csrc, cdef, entry)` returns the unit's object from the store (`~/.cache/pycircuit/hdl`
or a per-process directory with the cache off), keyed by the exact source text and the flags,
built under the key's lock when absent. `CFLAGS` are strict IEEE (no fast-math, no FMA
contraction, no `-march`, every transcendental a real libm call): the bit identity rests on them.
The FFI a unit is called through comes from `_cbuild._ffi_for`: one per header text a process,
precompiled into the store as cffi's out-of-line module (no header or type string is parsed at
run time), parsed as before when that cannot be written or read. Every opened library is held for
the process.

## Adding a C unit

1. `csrc/X.h` (structs, the entry `pcx_X`, the codes as named enums; no `#include`) and `csrc/X.c`
   (includes `X.h` once); shared types in `pcx_abi.h`, not in a unit's header.
2. The glue `_tran_X_c.py`: its switch, `driver()`, the eligibility checks and their decline
   reasons, the call, the write-back; its docstring says what it stands for, what it keeps on
   others' objects, its `_paths` prefix, what it declines and what it hands back.
3. `tests/test_c_abi.py` `UNITS`; `tests/test_X_c.py` (A/B, drawn twin, every decline, alive);
   a pinned pair in `tests/test_pinned_pairs.py`; any new transient attribute in `_RUN_STATE`.
4. If the C loop drives it: `loop.h` fields and a bail code, `loop.c`'s call, and in
   `_tran_loop_c` the genuine pieces, the bind and the write-back, with the counts it re-creates.
5. This page: its row in the layers, what it owns, its switch.
6. No function of the unit or its glue over 80 lines (`tests/test_runtime_map.py`, which lists
   the files): a long one is split into its jobs, or says `KEPT WHOLE:` with the measured cost
   of the split it refused.  A split of a per-step function adds a call to every step -- time it
   standalone against its parent before the check.  A glue helper that declines counts its own
   reason; its caller passes the None on with `# paths: counted by <helper>`, the helper itself
   in `tests/test_paths.py`'s `DECLINING`.
