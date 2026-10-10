
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
