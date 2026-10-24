"""BLAS thread control for the analyses that make many small solves.

Circuit matrices are small (n ~ 10^2), so a threaded LAPACK spends more time
spawning and synchronising threads than doing the work: the penalty is
per-CALL thread-pool overhead, crossing over around 50-100 solves, not a
property of problem size.  An analysis that makes many small assemblies and
solves holds BLAS at one thread for its whole run:

* `Transient.solve` -- 1.72x on the leapfrog on 4 cores, 14-20x on 24;
* `PSS.solve` -- its walks drive the inner transient's steps directly;
* `DC.solve`, where the AutoSolver's dense LU is past the size OpenBLAS
  threads (`THREADED_FROM` unknowns) and below its sparse floor -- the C
  Newton's ladders make hundreds to thousands of solves: over the 26 gold
  operating points 17.9 -> 11.1 s, voter25 30x, gm17 3.3x, no deck slower,
  and the 10-30x run-to-run swings of the threaded solve gone.  Smaller
  systems are not threaded, and the sparse ones factor with KLU, which the
  threads do not reach: there the limit would only cost its own call.

* the small-signal analyses -- `AC.solve`, `Noise.solve`,
  `TransimpedanceAnalysis.solve` (`blas_held`) and what is built on them
  (`TwoPortAnalysis`, `FeedbackLoopAnalysis`): one dense solve or more at
  each frequency.  Over 61 frequencies at eight gold decks' operating
  points, 1.19-1.48x from 171 to 295 unknowns, the same at 88 and 429,
  0.95x at 519 (where a dense complex LU starts to pay for its threads);
  the answers move at rounding level (5e-9 relative at worst, an
  ill-conditioned system's).  And a sweep's parallel workers no longer
  oversubscribe the cores: one measured ~1 s a point threaded against
  0.35-4.5 ms on one thread.

The periodic small-signal analyses after a PSS (`PAC`, `pnoise`) are not
held: not measured past the threading size.

`threadpoolctl` is a dependency.  The import stays guarded anyway, and
`blas_single_thread_available()` stays in the API, so a stripped environment
degrades to threaded BLAS rather than failing to import a circuit simulator.

History (the measurements, the 2026-08-31 refusal for DC that DC's later
measurement overturned, and AC's, overturned the same way):
`doc/transient_history.md`, `_threadpool_limits`.
"""
import contextlib
import functools
import sys

try:
    from threadpoolctl import threadpool_limits as _threadpool_limits
except ImportError:
    _threadpool_limits = None

#: the size from which OpenBLAS threads a dense LU (dgetrf), in unknowns of the
#: matrix solved: below it the limit changes nothing (g1310's 97 did not move,
#: gm1's 129 ran 6.6x faster)
THREADED_FROM = 100


def blas_single_thread_available():
    """True if the BLAS thread limit can be applied from inside the process."""
    return _threadpool_limits is not None


## THE CONTROLLER, KEPT: building one rescans every shared library the
## process has loaded -- 3.4 M instructions, two steps of a 20-MosLevel1
## transient, on every `solve` (`threadpool_limits` builds one per call) --
## for a set that changes when a package is imported: a BLAS arrives with its
## package.  So one, rebuilt when `len(sys.modules)` moved; a library loaded
## through ctypes alone is not seen until then.  The limit itself is the same
## call on it.  In a container: the leak detector reads module-level scalars.
_BLAS_CONTROLLER = {}


def _single_threaded_blas():
    """Limit BLAS to one thread for the duration, if that is possible here.

    NOTE ON REPRODUCIBILITY: a threaded LAPACK reduces in a different order, so
    results differ in the last bits between one and many threads.  Measured on the
    leapfrog: identical step count, `max|dv| = 2.16e-15 V`.  That is BLAS
    non-associativity, not a change of method -- but it does mean a run is only
    bit-reproducible against another run with the same thread count.
    """
    if _threadpool_limits is None:
        return contextlib.nullcontext()
    n = len(sys.modules)
    rec = _BLAS_CONTROLLER.get('blas')
    if rec is None or rec[0] != n:
        from threadpoolctl import ThreadpoolController
        rec = _BLAS_CONTROLLER['blas'] = (n, ThreadpoolController())
    return rec[1].limit(limits=1, user_api='blas')


def blas_held(solve):
    """`solve` with BLAS held at one thread for its whole run
    (`_single_threaded_blas`): the small-signal analyses' (the module's
    list)."""
    @functools.wraps(solve)
    def held(*args, **kwargs):
        with _single_threaded_blas():
            return solve(*args, **kwargs)
    return held
