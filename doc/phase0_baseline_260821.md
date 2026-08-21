# Phase-0 baseline — 2026-08-21

Recorded per `doc/transient_review_260820.md`, Phase 0 item 6, after the Phase-0 fixes
(F2(a), F6(a), F8, F9(a), F18) landed. These are the yardsticks the Phase-2/3 fixes are
judged against: each behavior-changing fix records its own delta as it lands; suite
expectations are pinned once, at the Phase-3 exit gate.

Machine: ThinkPad P16 Gen 1, i7-12800HX, RTX A1000 (4 GB), venv Python 3.14.4,
jax 0.11.1/CUDA 13, `XLA_PYTHON_CLIENT_PREALLOCATE=false`.

## Performance (benchmarks/phase0_baseline.py — one warm-up, one timed run)

rc-vsin (VSin 1 kHz → R 1k → C 100n), `tend=2e-3`, `timestep=1e-5`:

| configuration | wall | accepted | rejected |
|---|---|---|---|
| CPU `Transient`, Euler (shipped default) | 150.1 ms | 412 | 8 |
| CPU `Transient`, Gear2 | 78.7 ms | 211 | 0 |
| JAX `JAXTransient`, gear (shipped default) | 612.5 ms | 210 | 0 |

Three facts worth naming, all relevant to the review:

- **At matched order the backends agree almost exactly** (211 vs 210 accepted steps) —
  the step-count parity F19(b) says is only visible once the integrator is pinned. The
  shipped defaults compare 412 vs 210, i.e. different methods, which is F19(b)'s point.
- **The JAX path is ~7.8× slower than CPU-Gear2 on this single small circuit.** Not a
  defect — n≈3 matrices cannot feed a GPU, and the JAX path's value proposition is
  `solve_batched` — but until now that trade-off was an impression, not a number. Any
  Phase-3 fix that moves this ratio materially in either direction should say so.
- **CPU-Euler's 8 rejections vs Gear2's 0** on the same circuit is the F17-adjacent
  population to watch when the JAX safety factor lands.

## Coverage (coverage 7.13, `--source=pycircuit/circuit`, the eleven transient-related
test files at `-m ""`; 204 tests, 72 s)

| module | stmts | miss | cover |
|---|---|---|---|
| `_lte_kernels.py` | 107 | 9 | 92% |
| `integrator.py` | 97 | 10 | 90% |
| `transient.py` | 735 | 110 | 85% |
| `jaxtransient.py` | 477 | 73 | 85% |
| `stepcontroller.py` | 174 | 49 | 72% |
| `nrsolver.py` | 241 | 72 | 70% |

Dark territory worth knowing before trusting a green suite (from `--show-missing`):

- **`transient.py:805–857` — `_apply_voltage_ics`, the spanning-tree capacitor-IC solve,
  is essentially untested**, along with most of the element-IC machinery (686–704,
  728–752). The review praised its design from reading; nothing currently guards it.
- `transient.py:536–549` — the `_newton` exception-classification branches
  (SingularMatrix promotion) are unexercised.
- `nrsolver.py` at 70% — the DampedNewton backtracking paths (the F15 territory) and
  `JAXNewtonSolver` are largely dark, consistent with F15's "delete-or-port; grep for
  users first".
- `stepcontroller.py` at 72% — the PI controller and parts of the band logic (the F10
  territory).

Reproduce: `XLA_PYTHON_CLIENT_PREALLOCATE=false .venv/bin/python -m coverage run
--source=pycircuit/circuit -m pytest -q -m "" <the eleven files>` then
`coverage report --show-missing`.
