# The batched-sweep payoff — 2026-08-21

The number this branch exists for.  The Phase-0 baseline recorded the
single-circuit cost of the JAX path (7.8× slower than CPU-Gear2 at n≈3 —
expected: tiny matrices cannot feed a GPU) and noted that the value
proposition, `solve_batched`, was "an impression, not a number".  This is
the number.

**Benchmark:** `benchmarks/batched_sweep.py` — an N-lane resistor corner
sweep of the half-wave rectifier (VSin 5 V / 1 kHz → Diode → R ∥ 100 nF,
r log-spaced 300 Ω … 30 kΩ), `tend = 2 ms`, `timestep = 1e-5`, shipped
defaults (Gear-2 on both backends, per P22).  The JAX side is ONE
`solve_batched` call — per-lane DC bias through the P18/P25 continuation
ladders included, junction lanes and all.  The CPU side is the honest
alternative a user writes today: a Python loop over `Transient`,
rebuilding the circuit per lane.  Every row is fully measured — the CPU
loop really runs all N lanes at every N; nothing is extrapolated.

Machine: ThinkPad P16 Gen 1, i7-12800HX, RTX A1000 (4 GB, CUDA 13,
jax 0.11.1), venv Python 3.14.4, `XLA_PYTHON_CLIENT_PREALLOCATE=false`.

| lanes | CPU loop | JAX cold (incl. jit) | JAX warm | speedup (warm) | max final-sample Δv |
|---|---|---|---|---|---|
| 8   | 0.87 s  | 3.66 s | 1.93 s | 0.5× | 7.5e-4 V |
| 32  | 3.60 s  | 3.39 s | 2.15 s | 1.7× | 7.6e-4 V |
| 128 | 13.88 s | 3.31 s | 2.01 s | 6.9× | 7.3e-4 V |
| 512 | 54.71 s | 4.19 s | **2.43 s** | **22.5×** | 9.7e-4 V |

Facts worth naming:

- **The crossover is ~16 lanes** on this machine, jit compile amortized.
  Counting the compile against a single cold run moves it to ~32 lanes.
  Below that, the CPU loop wins and is the honest recommendation.
- **The JAX wall is nearly flat in N**: 1.93 s → 2.43 s across a 64×
  lane increase (~1 ms marginal cost per lane at 512).  At this circuit
  size the vmapped loop is latency-bound (kernel-launch and
  chunk-boundary Python), not throughput-bound — the GPU has capacity
  far beyond 512 lanes of an n=3 circuit, so the speedup keeps growing
  roughly linearly with N until that changes.
- **All lanes march at the slowest lane's step count** (a vmapped
  `while_loop` runs until every lane's condition clears).  The r sweep
  spans a 100× RC spread, so this cost is IN the measured numbers, not
  hidden by a homogeneous sweep.
- **Cross-backend agreement ≤ 1 mV on a 5 V drive** at the final sample,
  per lane, at matched order (both backends default Gear-2) — final
  points land on different adaptive grids, so this is a loose functional
  check, not the conformance harness's bit-level one.
- **Deferred-item triggers, checked and NOT fired:** the flat wall means
  finished-lane dispatch (batch compaction, F1's rejected strategy (b))
  is invisible at these sizes, and this benchmark does not profile
  breakpoint cost, so `Sin.next_event` quarter-period coalescing stays
  parked on its own trigger.

Reproduce: `XLA_PYTHON_CLIENT_PREALLOCATE=false .venv/bin/python
benchmarks/batched_sweep.py`.
