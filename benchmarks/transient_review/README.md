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

The numerical lens's own scripts (`e1`-`e15`: error constants, the trapezoidal `(-1)^n`
contamination, PI-controller stability, BDF-2 step-ratio limits, DAE index) were not
recovered before the scratchpad was reused. Their results are recorded in
`doc/transient_review.md` section 4 with enough detail to reconstruct them, but the
scripts themselves are gone -- worth knowing before citing those numbers as reproducible.
