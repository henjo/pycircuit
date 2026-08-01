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
        try:
            self._sparse = SuperLUSolver()
            self._choice = self._sparse
        except ImportError:                                # pragma: no cover
            self._choice = self._dense
        return self._choice

    def solve(self, A, b, toolkit):
        return self._select(A).solve(A, b, toolkit)

    def __repr__(self):
        return ('AutoSolver(min_n=%d, max_fill=%.2f, chose=%r)'
                % (self.min_n, self.max_fill, self._choice))
