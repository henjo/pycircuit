# Evidence for `doc/transient_review.md`

Probe scripts from the four-lens transient review, 2026-07-30, preserved because the
review's numbers are only worth what their reproduction is worth. They were written in a
session scratchpad that does not survive, which is why they are here.

These are **probes, not benchmarks**: written to answer one question each, not maintained,
not wired into the suite, and several hard-code paths or assume the tree state at
`a675455` (i.e. BEFORE the GBW compensation raised the leapfrog to 136 unknowns). Re-run
them expecting to fix a path, not expecting them to pass a gate.

| script | what it establishes |
|---|---|
| `prof_tran.py` | cProfile: assembly 96%, linear solve 2.1%, `Toolkit.__getattr__` 547k calls |
| `proto_assembly.py` | hoisted-`hasattr` + `bincount` assembly, correctness `0.00e+00`, **2.0x** |
| `compound.py` | the **10.5x** compound A/B, 8.80 s -> 0.84 s, waveform drift `0.00e+00` |
| `threadab.py`, `densecheck.py` | BLAS threads: 0.238 ms vs 4.462 ms; 2.3-2.5x on a full transient |
| `csc_reuse.py` | dense vs `spsolve` vs fixed-pattern `splu` vs resolve, real Jacobians n=29..542 |
| `scaling.py` | the same to n=5000, with LU fill and dense memory (200 MB at n=5000) |
| `probe_sparse.py` | `sparse_numeric` crashes when the *circuit* uses it (`AxisError`, `ValueError`) |
| `diag.py` | diagnostics: the singular-matrix message, parameter surfaces, absent run statistics |
| `silentdc.py` | the silent-DC-failure demo, and per-step assembly counts (G/C/i/u 3.17, q 5.22) |

Run with `PYTHONPATH=<repo> MPLBACKEND=Agg python3 -u <script>`.

## Stage 0 gate evidence (2026-07-30)

These three are **not** probes in the sense above. They were written to answer a numbered
question in `doc/transient_work_plan.md`, they run against the tree as it is, and their
output is quoted in that document's `OUTCOME:` lines. Keep them runnable.

| script | question | result |
|---|---|---|
| `stage0_2a_integrator_choice.py` | 0.2a — is the IM3 harness's integrator choice an artefact of breakpoint seeding? | No. Suppressing `Sin.next_event` moves the step count 1.8%, and the ranking is 1.199 vs 1.206. The recorded 10x stands. |
| `stage0_2b_lte_solve_cost.py` | 0.2b — what does the LTE's `J^-1` solve cost per step? | Under the 10% decision threshold. Also: there is no `G` evaluation to avoid, so half the charge-domain path's claimed advantage never existed. |
| `stage0_1d_coupled_livelock.py` | 0.1d — what does `_solve_coupled` do when Newton fails? | Livelocks. Each outer iteration costs 10 Newton solves and buys `h*0.25^10`; measured rate matches that prediction to a ratio of 1.0000. |
| `stage0_1c_pac_memory.py` | 0.1c — how big is PAC's dense operator? | 420 GiB at N=137, M=1000. Exists because the review's "~150 TB" was a typed number and wrong by 365x; this one is generated. |

The first two need `benchmarks/` on `PYTHONPATH` as well as the repo root, because they
import the IM3 harness:

```
PYTHONPATH=<repo>:<repo>/benchmarks MPLBACKEND=Agg python3 -u <script>
```

`stage0_2a` and `stage0_2b` each run several full transients on the leapfrog and take
5-10 minutes; `stage0_1d` is seconds.

## PCNR design-study evidence (`doc/pcnr_native_design.md`)

Like the stage 0 scripts, these run against the tree as it is, and the design study quotes
their output.

| script | question | result |
|---|---|---|
| `companion_consistency.py` | Is `Diode.G` the exact derivative of `Diode.i`? | Yes, to 5.0e-10 by central difference. The documented claim that they are taken at different points was false. |
| `pcnr_narrow_fix.py` (+ `.log`) | Can `Diode` stop reading a stale `_vlim` without touching PCNR assembly? | No. On a 1 A diode the PCNR node reaches 2.58e11 V, an honest diode overflows `exp`, and 199 of 200 assemblies are nan. The narrow fix needs the skip set first. |
| `stage13_10_11_gates.py` | Gates 13-10 (`_pnjlim`'s 2 VT condition) and 13-11 (skip on the JAX batched path); also the entry evidence for 12-6. Modes `table`, `snapshot`/`compare` (25 DC and transient cases), `mismatch`, `jax` (needs the repo `.venv`). | Steps <= 2 VT no longer limited; 13-9c's 1518 mismatches -> 0; CE stage +1 iteration (traced: SPICE's full step overshoots further); PCNR on JAX matches numeric to 1.1e-17. |
| `stage13_8_9_gates.py` | Gates 13-8/13-9: the skip set, then the classic-limiting scope. `snapshot`/`compare` check bit-identity over 18 cases; `probes` measures call counts, assembly cost, Jacobian honesty and PSS. | Classic paths bit-identical except one explained 3.2e-14; PCNR iterations 29% cheaper than limiting; `cir.G` exact outside limiting; PSS on a rectifier 0.000 V → 4.2481 V (reference 4.2485). |

The numerical lens's own scripts (`e1`-`e15`: error constants, the trapezoidal `(-1)^n`
contamination, PI-controller stability, BDF-2 step-ratio limits, DAE index) were not
recovered before the scratchpad was reused. Their results are recorded in
`doc/transient_review.md` section 4 with enough detail to reconstruct them, but the
scripts themselves are gone -- worth knowing before citing those numbers as reproducible.
