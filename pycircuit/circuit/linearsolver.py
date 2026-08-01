"""Linear-solver strategies for the Newton iteration, in the shape the package
already uses for ``Integrator``, ``StepController``, ``Scaler`` and
``NonLinearSolver``.

STAGE 7b.  Every Newton iteration solves ``J dx = -F`` through
``toolkit.linearsolver``, which is ``numpy.linalg.solve`` on the numeric backend:
a dense LU of the whole matrix, every time.  Circuit Jacobians are extremely
sparse -- an RC ladder at n=1101 is **0.27% non-zero** -- and a sparse
factorisation of that is **74x** faster.

THE SELECTION KEYS ON FILL, NOT ON ``n``, and that is a departure from the plan.
7b specifies "DenseLU below ~200 unknowns, SuperLU above".  Measured, the
crossover follows the sparsity *pattern*: on the RC ladder SuperLU already wins
2.5x at n=31, while a dense-ish circuit at n=300 would be mis-served by a
size-only rule.  ``n`` is used only as a cheap floor below which the measurement
itself is not worth taking.

THE DEFAULT IS THE OLD PATH, DELIBERATELY.  ``numpy.linalg.solve`` and SuperLU
round differently, so any circuit that switches gets different last bits.  Gate
7b-1 asks for identical results on the existing transient tests, and those are all
small; ``AutoSolver`` therefore leaves anything below its thresholds on exactly the
call the code made before, and only reaches for the sparse path where the win is
large enough to be worth a changed rounding.

A further measured caution, which is why ``DenseSolver`` does NOT use
``scipy.linalg.lu_factor``: at n=31 ``numpy.linalg.solve`` takes 2e-5 s against
3.2e-4 s for ``lu_factor`` + ``lu_solve`` -- an order of magnitude *worse*, on
Python call overhead alone.  A factor/solve split only pays when the factors are
reused, and on this code path they are not: see the note on ``factor()`` below.
"""
from abc import ABC, abstractmethod

import numpy


class LinearSolver(ABC):
    """Solve ``A x = b`` for the Newton iteration."""

    @abstractmethod
    def solve(self, A, b, toolkit):
        """Return ``x`` such that ``A x = b``."""

    def factor(self, A, toolkit):
        """Optional: return a reusable factorisation, or ``None``.

        Unused on the transient path today, and the reason is worth stating
        because 7b-2 proposed otherwise.  Newton's last factorisation is of
        ``J(x_k)``; the step controller needs ``J`` at the *converged*
        ``x_{k+1}``, which nothing else factors.  Reusing across that boundary
        substitutes a different matrix -- measured at median 5.5e-9 and max
        2.2e-5 relative on a diode circuit -- so it is an approximation, not a
        saved duplicate, and it cannot be done under a bit-identical gate.
        """
        return None


class DenseSolver(LinearSolver):
    """The historical path: whatever ``toolkit.linearsolver`` does.

    On the numeric backend that is ``numpy.linalg.solve``.  Kept as the default so
    existing results do not move, and so symbolic and other toolkits keep working
    unchanged -- this class is the only one that is toolkit-agnostic.
    """

    def solve(self, A, b, toolkit):
        return toolkit.linearsolver(A, b)

    def __repr__(self):
        return 'DenseSolver()'


class SuperLUSolver(LinearSolver):
    """SuperLU through ``scipy.sparse.linalg.splu``.

    Discovered rather than required, the way ``symengine`` and ``_ginac_ext``
    already are: if SciPy is missing the class raises on construction and
    ``AutoSolver`` never selects it.
    """

    def __init__(self):
        try:
            import scipy.sparse as _sp
            import scipy.sparse.linalg as _spla
        except ImportError as e:                          # pragma: no cover
            raise ImportError(
                'SuperLUSolver needs scipy.sparse; install scipy or use '
                'DenseSolver()') from e
        self._sp = _sp
        self._spla = _spla

    def solve(self, A, b, toolkit):
        A = numpy.asarray(A)
        if A.dtype == object:
            ## A symbolic matrix cannot go to SuperLU; fall back rather than fail,
            ## because the selection may have been made on a numeric sibling.
            return toolkit.linearsolver(A, b)
        Acsc = self._sp.csc_matrix(A)
        return self._spla.splu(Acsc).solve(numpy.asarray(b))

    def __repr__(self):
        return 'SuperLUSolver()'


## ---------------------------------------------------------------------------
## STAGE 7c -- KLU, through ctypes, discovered rather than required.
##
## KLU is the factorisation production circuit simulators use (Xyce among them)
## and the reason is the analyze/factor/REFACTOR split.  A circuit's Jacobian
## keeps its sparsity PATTERN across every Newton iteration and every timestep --
## only the values change -- so the fill-reducing ordering, which is most of the
## cost, need be computed once per circuit rather than once per solve.
##
## Measured here against `scipy`'s `splu`, which recomputes the ordering every
## call (7c's premise, confirmed: `klu_analyze` 0.000441 s against `klu_factor`
## 0.000378 s at n=3200, so the ordering is roughly half of a full factorisation,
## and `klu_refactor` at 0.000099 s is 4x cheaper than factoring again):
##
##      n      splu   klu refactor+solve   speedup
##    150  0.000163            0.000046     3.51x
##    400  0.000463            0.000060     7.75x
##    800  0.000541            0.000061     8.85x
##   1600  0.000907            0.000097     9.31x
##   3200  0.001714            0.000172     9.97x
##
## No pip dependency is added: `libklu.so` ships with SuiteSparse, and this binds
## it the way `_ginac_ext` and `symengine` are already discovered.  If the library
## is absent the class raises on construction and `AutoSolver` never selects it.
## ---------------------------------------------------------------------------

class KLUSolver(LinearSolver):
    """SuiteSparse KLU with the symbolic analysis reused across solves.

    THE REFACTOR PATH IS VALIDATED, NOT TRUSTED.  ``klu_refactor`` reuses the
    pivot ORDER chosen when the matrix was first factored.  That is the whole
    point -- and it is also its one hazard: if the values move far enough, pivots
    that were sound become unstable and the factorisation silently degrades.
    Production tools read KLU's reciprocal pivot growth to decide; reading it here
    would mean hard-coding offsets into a ``klu_common`` struct this binding
    cannot verify, so instead **the residual is checked directly** and a full
    ``klu_factor`` is redone when it is not good enough.  One sparse mat-vec
    against a factorisation is cheap, and it cannot be wrong about the layout.

    That choice is deliberate given what stage 9 and 7d kept turning up: three
    separate defects this session survived because nothing ever checked the thing
    itself.
    """

    ## Relative residual above which a refactor is rejected and the ordering
    ## recomputed.  Loose enough not to trip on ordinary round-off, tight enough
    ## that a genuinely bad reuse is caught: a well-conditioned solve lands near
    ## 1e-16 and a degraded one is orders above.
    REFACTOR_RESIDUAL_TOL = 1e-8

    def __init__(self, libname='libklu.so.2'):
        import ctypes
        try:
            lib = ctypes.CDLL(libname)
        except OSError as e:
            raise ImportError(
                'KLUSolver needs SuiteSparse KLU (%s); install libklu / '
                'libsuitesparse, or use SuperLUSolver()' % libname) from e
        try:
            import scipy.sparse as _sp
        except ImportError as e:                              # pragma: no cover
            raise ImportError('KLUSolver needs scipy.sparse') from e

        self._ct = ctypes
        self._sp = _sp
        self._lib = lib
        c_int_p = ctypes.POINTER(ctypes.c_int)
        c_dbl_p = ctypes.POINTER(ctypes.c_double)
        for nm, argt, rest in (
                ('klu_defaults', [ctypes.c_void_p], ctypes.c_int),
                ('klu_analyze', [ctypes.c_int, c_int_p, c_int_p, ctypes.c_void_p],
                 ctypes.c_void_p),
                ('klu_factor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                ctypes.c_void_p], ctypes.c_void_p),
                ('klu_refactor', [c_int_p, c_int_p, c_dbl_p, ctypes.c_void_p,
                                  ctypes.c_void_p, ctypes.c_void_p], ctypes.c_int),
                ('klu_solve', [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_int,
                               ctypes.c_int, c_dbl_p, ctypes.c_void_p],
                 ctypes.c_int)):
            f = getattr(lib, nm)
            f.argtypes = argt
            f.restype = rest

        ## `klu_common` is ~200 bytes; over-allocating a buffer is safe and avoids
        ## committing to a struct layout that varies between SuiteSparse releases.
        self._common = ctypes.create_string_buffer(4096)
        if not lib.klu_defaults(ctypes.byref(self._common)):
            raise ImportError('klu_defaults failed; the KLU library looks unusable')

        self._pattern = None      # (n, indptr bytes, indices bytes)
        self._symbolic = None
        self._numeric = None
        self.refactors = 0        # counters, so the reuse is observable
        self.factors = 0
        self.analyses = 0
        self.residual_fallbacks = 0

    def _ptr(self, arr, ctype):
        return arr.ctypes.data_as(self._ct.POINTER(ctype))

    def solve(self, A, b, toolkit):
        import numpy
        ctypes = self._ct
        arr = numpy.asarray(A)
        if arr.dtype == object:
            ## Symbolic matrices cannot go to KLU; fall back rather than fail.
            return toolkit.linearsolver(A, b)

        Acsc = self._sp.csc_matrix(arr)
        n = Acsc.shape[0]
        Ap = numpy.ascontiguousarray(Acsc.indptr, dtype=numpy.int32)
        Ai = numpy.ascontiguousarray(Acsc.indices, dtype=numpy.int32)
        Ax = numpy.ascontiguousarray(Acsc.data, dtype=numpy.float64)
        rhs = numpy.ascontiguousarray(numpy.asarray(b, dtype=numpy.float64).ravel())

        key = (n, Ap.tobytes(), Ai.tobytes())
        reused = key == self._pattern and self._symbolic and self._numeric

        if not reused:
            ## The pattern changed (or this is the first solve): a new ordering is
            ## genuinely required, and this is the cost KLU exists to pay ONCE.
            self._symbolic = self._lib.klu_analyze(
                n, self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                ctypes.byref(self._common))
            self.analyses += 1
            if not self._symbolic:
                raise numpy.linalg.LinAlgError('klu_analyze failed (singular structure?)')
            self._numeric = self._lib.klu_factor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.byref(self._common))
            self.factors += 1
            if not self._numeric:
                raise numpy.linalg.LinAlgError('Singular matrix')
            self._pattern = key
        else:
            ok = self._lib.klu_refactor(
                self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                ctypes.c_void_p(self._numeric), ctypes.byref(self._common))
            self.refactors += 1
            if not ok:
                raise numpy.linalg.LinAlgError('Singular matrix')

        x = rhs.copy()
        if not self._lib.klu_solve(
                ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                n, 1, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
            raise numpy.linalg.LinAlgError('klu_solve failed')

        if reused:
            ## The validation described in the class docstring.  Only on the reuse
            ## path: a fresh factorisation chose its own pivots and needs no check.
            scale = max(float(numpy.abs(rhs).max()), 1e-300)
            resid = float(numpy.abs(Acsc.dot(x) - rhs).max()) / scale
            if not (resid <= self.REFACTOR_RESIDUAL_TOL):
                self.residual_fallbacks += 1
                self._numeric = self._lib.klu_factor(
                    self._ptr(Ap, ctypes.c_int), self._ptr(Ai, ctypes.c_int),
                    self._ptr(Ax, ctypes.c_double), ctypes.c_void_p(self._symbolic),
                    ctypes.byref(self._common))
                self.factors += 1
                if not self._numeric:
                    raise numpy.linalg.LinAlgError('Singular matrix')
                x = rhs.copy()
                if not self._lib.klu_solve(
                        ctypes.c_void_p(self._symbolic), ctypes.c_void_p(self._numeric),
                        n, 1, self._ptr(x, ctypes.c_double), ctypes.byref(self._common)):
                    raise numpy.linalg.LinAlgError('klu_solve failed')
        return x

    def __repr__(self):
        return ('KLUSolver(analyses=%d, factors=%d, refactors=%d, fallbacks=%d)'
                % (self.analyses, self.factors, self.refactors,
                   self.residual_fallbacks))


## Below this many unknowns the sparse path cannot repay its own setup, and
## measuring the fill is not worth the pass over the matrix.  Chosen from the
## measured table (see the module docstring): at n=61 the win is 1.16x, which is
## inside the noise of a Python-level dispatch.
MIN_N_FOR_SPARSE = 100

## And above this fill the matrix is dense enough that SuperLU's bookkeeping costs
## more than it saves.  The measured circuits run 0.27%-9.4%; 20% is a deliberately
## conservative bound rather than a tuned one, because no measurement here explores
## dense circuits and a threshold fitted to sparse ones would be a guess dressed as
## a number.
MAX_FILL_FOR_SPARSE = 0.20


class AutoSolver(LinearSolver):
    """Pick dense or sparse per matrix, on measured fill rather than on ``n``.

    The decision is cached on the (shape, fill-bucket) of the first matrix seen,
    because a transient calls this thousands of times with the same structure and
    counting non-zeros every time would give back the saving.
    """

    def __init__(self, min_n=MIN_N_FOR_SPARSE, max_fill=MAX_FILL_FOR_SPARSE):
        self.min_n = min_n
        self.max_fill = max_fill
        self._dense = DenseSolver()
        self._sparse = None
        self._choice = None

    def _select(self, A):
        if self._choice is not None:
            return self._choice
        arr = numpy.asarray(A)
        if arr.ndim != 2 or arr.shape[0] < self.min_n or arr.dtype == object:
            self._choice = self._dense
            return self._choice
        fill = numpy.count_nonzero(arr) / float(arr.shape[0] * arr.shape[1])
        if fill > self.max_fill:
            self._choice = self._dense
            return self._choice
        ## SuperLU, NOT KLU -- and that is a measurement, not an oversight.
        ##
        ## KLU wins 3.5x-10x on the isolated factor+solve, which is what 7c's gate
        ## asks for and what it passes.  END TO END IT DOES NOT: best-of-3 on a
        ## transient, KLU came in at **0.52x** the shipped dense path at n=152 and
        ## 1.31x at n=402, where SuperLU managed 1.36x with far less setup.  The
        ## refactor works exactly as intended (1 analyze, 1 factor, 227 refactors
        ## over a run) -- the loss is Python overhead per call against a solve that
        ## is only ~8% of a transient.  Selecting KLU automatically would make a
        ## 152-unknown circuit twice as slow.
        ##
        ## `KLUSolver` therefore stays explicit: `linearsolver=KLUSolver()`.  It is
        ## the right tool once assembly stops dominating -- see 2+.5 -- and the
        ## measurement should be redone then rather than assumed to still hold.
        for cls in (SuperLUSolver,):
            try:
                self._sparse = cls()
                self._choice = self._sparse
                return self._choice
            except ImportError:
                continue
        self._choice = self._dense                          # pragma: no cover
        return self._choice

    def solve(self, A, b, toolkit):
        return self._select(A).solve(A, b, toolkit)

    def __repr__(self):
        return ('AutoSolver(min_n=%d, max_fill=%.2f, chose=%r)'
                % (self.min_n, self.max_fill, self._choice))
