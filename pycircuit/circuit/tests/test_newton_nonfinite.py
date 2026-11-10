"""`StandardNewton` at an iterate that is not finite (`STANDARD_NONFINITE`):
the point is evaluated once and its step solved, and there the iteration
ends -- a singular solve raising "Singular Jacobian" as before, which the DC
and `GminAnchorNewton` read as singular (`_is_singular_failure`), any other
the non-finite text, which none reads so.  A compact model's chain with no
limiter turns an iterate infinite on its DC ladders' rungs; past that point
each evaluation was at a point not finite and the iteration went on to
`maxiter` for nothing."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import circuit
from pycircuit.circuit.nrsolver import (
    STANDARD_NONFINITE,
    NoConvergenceError,
    StandardNewton,
    _is_singular_failure,
)

numeric = circuit.numeric


def _run(singular_past_it):
    """`x - 1 = 0` in two unknowns from 0, the first's derivative 1e-310:
    the first step, 1 / 1e-310, overflows to an infinite iterate.  There
    the Jacobian is as before, or [[1, 1], [1, 1]] (`singular_past_it`:
    singular in its values, every row and column holding an entry, so no
    structural test reads it).  The exception and the evaluations made."""
    calls = []

    def eval_FJ(x):
        calls.append(np.array(x, dtype=float))
        finite = bool(np.isfinite(x).all())
        J = (np.ones((2, 2)) if (singular_past_it and not finite)
             else np.array([[1e-310, 0.0], [0.0, 1.0]]))
        return np.asarray(x, dtype=float) - 1.0, J
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        with pytest.raises(NoConvergenceError) as exc:
            StandardNewton().solve_system(np.zeros(2), eval_FJ, numeric, 1e-6, 1e-12, 1e-12, 20)
    return exc.value, calls


def test_a_non_finite_iterate_ends_the_iteration_after_its_step():
    """The step at the infinite iterate solves: the iteration ends there,
    at its second evaluation, with the non-finite text -- not read as
    singular -- where it went on to the iteration limit."""
    e, calls = _run(singular_past_it=False)
    assert len(calls) == 2 and not np.isfinite(calls[1]).all()
    assert str(e) == STANDARD_NONFINITE % 2
    assert not _is_singular_failure(e)


def test_a_non_finite_iterate_whose_step_is_singular_says_so():
    """The Jacobian at the infinite iterate is singular: the solve there
    raises, and the iteration says "Singular Jacobian" -- read as singular,
    so the DC's anchor rescue still sees it -- as it did."""
    e, calls = _run(singular_past_it=True)
    assert len(calls) == 2
    assert str(e).startswith('Singular Jacobian') and _is_singular_failure(e)
