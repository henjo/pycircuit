"""BLAS thread control for the analyses that make many small solves.

Circuit matrices are small (n ~ 10^2), so a threaded LAPACK spends more time
spawning and synchronising threads than doing the work: the penalty is
per-CALL thread-pool overhead, crossing over around 50-100 solves, not a
property of problem size.  An analysis that makes many small assemblies and
solves holds BLAS at one thread for its whole run:

* `Transient.solve` -- 1.72x on the leapfrog on 4 cores, 14-20x on 24;
* `PSS.solve` -- its walks drive the inner transient's steps directly;
* `DC.solve`, where the AutoSolver's system is past the size OpenBLAS
  threads (`THREADED_FROM` unknowns) -- below its sparse floor the C
  Newton's ladders make hundreds to thousands of dense solves: over the 26
  gold operating points 17.9 -> 11.1 s, voter25 30x, gm17 3.3x, no deck
  slower, and the 10-30x run-to-run swings of the threaded solve gone; on
  KLU above it the solve's other dense calls are threaded still -- an
  active chain's DC at 4 000 stages 1.34x, a PSP chain's at 1 200 1.16x,
  meshes and a quiet chain even.  Smaller systems are not threaded.

* the small-signal analyses -- `AC.solve`, `Noise.solve`,
  `TransimpedanceAnalysis.solve` (`blas_held_below`) and what is built on
  them (`TwoPortAnalysis`, `FeedbackLoopAnalysis`) -- below
  `AC_THREADS_FROM` unknowns: one dense complex solve or more at each
  frequency.  Over 61 frequencies at eight gold decks' operating points,
  1.19-1.48x from 171 to 295 unknowns; on an RC ladder 1.11-1.27x from 303
  to 423, and from 443 up the threads pay (0.84-0.97x to 503, 0.46x at
  1503).  The answers move at rounding level (5e-9 relative at worst, an
  ill-conditioned system's), and a sweep's parallel workers no longer
  oversubscribe the cores below the size (one measured ~1 s a point
  threaded against 0.35-4.5 ms on one thread);
* the periodic small-signal analyses after a PSS -- every public method
  of `PAC` (`blas_held_methods`: `solve`, `pnoise`, the spectra, the
  covariances, the adjoint rows): many circuit-size operations through the
  period, where the threads never paid -- one thread 0.97-1.76x for
  `solve`, 1.11-1.57x for `pnoise`, from 53 to 803 unknowns under gear and
  Radau.

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
## A library its loader opens itself (`dlopen`) comes with no import: that
## loader says so (`library_loaded`) -- CHOLMOD brings the system's OpenBLAS.
_BLAS_CONTROLLER = {}

## THE HOLDS OPEN NOW, outermost first.  A library reported while one is open
## arrived inside it (CHOLMOD loads when a transient's solver is built, in the
## transient's own hold), so the limit the hold set on entry missed it: it is
## held from then to the outermost hold's end.  Unheld, ibmpg3t's CHOLMOD ran
## on 24 threads, its refactor 1.42x slower than on one, and the run that
## loaded it computed other bits than the runs after it.
_OPEN_HOLDS = []


def library_loaded():
    """A shared library opened by its loader rather than with a package (a
    BLAS can come with it): the kept controller forgotten, so the next hold
    builds one that sees it, and inside a hold the library held at one thread
    to the outermost hold's end (its own count back then)."""
    _BLAS_CONTROLLER.clear()
    if _OPEN_HOLDS and _threadpool_limits is not None:
        from threadpoolctl import ThreadpoolController
        _OPEN_HOLDS[0].later.append(ThreadpoolController().limit(limits=1, user_api='blas'))


class _Hold:
    """One hold: the kept controller's one-thread limit, set on entry and
    restored on exit -- after the limits `library_loaded` added while it was
    the outermost hold open, restored last in first out."""
    __slots__ = ('controller', 'later', 'limit')

    def __init__(self, controller):
        self.controller = controller
        self.later = []
        self.limit = None

    def __enter__(self):
        self.limit = self.controller.limit(limits=1, user_api='blas')
        _OPEN_HOLDS.append(self)
        return self

    def __exit__(self, *exc):
        _OPEN_HOLDS.remove(self)
        while self.later:
            self.later.pop().restore_original_limits()
        self.limit.restore_original_limits()
        return False


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
    return _Hold(rec[1])


#: the size from which the small-signal analyses' dense complex solve, one a
#: frequency, pays for OpenBLAS's threads: on an RC ladder (61 frequencies)
#: one thread 1.11-1.27x from 303 to 423 unknowns, 0.84-0.97x from 443 to
#: 503, 0.46x at 1503; gm19 (429) even, opampal (519) 0.95x
AC_THREADS_FROM = 440


def blas_held(solve):
    """`solve` with BLAS held at one thread for its whole run
    (`_single_threaded_blas`)."""
    @functools.wraps(solve)
    def held(*args, **kwargs):
        with _single_threaded_blas():
            return solve(*args, **kwargs)
    return held


def blas_held_below(threads_from):
    """`blas_held` for an analysis's method where its circuit has fewer
    than `threads_from` unknowns, BLAS's own threads from there (the
    small-signal analyses', `AC_THREADS_FROM`)."""
    def sized(solve):
        held = blas_held(solve)

        @functools.wraps(solve)
        def run(self, *args, **kwargs):
            if self.cir.n >= threads_from:
                return solve(self, *args, **kwargs)
            return held(self, *args, **kwargs)
        return run
    return sized


def blas_held_methods(cls):
    """`cls` with each public method, its own and its bases', `blas_held`
    (`PAC`'s: the module's list)."""
    import inspect
    for name in dir(cls):
        fn = inspect.getattr_static(cls, name)
        if not name.startswith('_') and inspect.isfunction(fn):
            setattr(cls, name, blas_held(fn))
    return cls
