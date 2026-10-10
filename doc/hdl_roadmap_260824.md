
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
