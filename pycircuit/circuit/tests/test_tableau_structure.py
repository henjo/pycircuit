"""A stage method's tableau classified once (speed round 11): the stage step
asks `is_stiffly_accurate` and `stage_structure` on every step
(`_solve_timestep_rk`), and each answer is kept with the tableau `butcher`
returned for it.  The answers as the formulas give them for every
structure; a run that does not classify its tableau again; a tableau made
anew classified anew."""
import numpy as np
import pytest

from pycircuit.circuit import circuit, integrator
from pycircuit.circuit.integrator import (
    ESDIRK43Integrator,
    RadauIIA3Integrator,
    RungeKuttaIntegrator,
    TRBDF2Integrator,
)
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.transient import Transient


def _answers(A, B, C):
    """The two answers by the formulas, as they stood."""
    stiff = bool(np.allclose(B, A[-1]) and abs(C[-1] - 1.0) < 1e-14)
    if not np.allclose(np.triu(A, 1), 0.0):
        return stiff, 'full'
    diag = np.diag(A)
    expl0 = abs(diag[0]) < 1e-14
    impl = diag[1:] if expl0 else diag
    equal = impl.size > 0 and np.allclose(impl, impl[0]) and impl[0] != 0.0
    return stiff, ('esdirk' if expl0 and equal else 'sdirk' if equal else 'dirk')


def _custom(A, B, C):
    return type('Custom', (RungeKuttaIntegrator,), {'A': A, 'B': B, 'C': C})()


CASES = {
    'radau': RadauIIA3Integrator,
    'trbdf2': TRBDF2Integrator,
    'esdirk43': ESDIRK43Integrator,
    'sdirk': lambda: _custom([[0.5, 0.0], [0.5, 0.5]], [0.5, 0.5], [0.5, 1.0]),
    'dirk': lambda: _custom([[0.25, 0.0], [0.5, 0.5]], [0.5, 0.5], [0.25, 1.0]),
    'not stiffly accurate': lambda: _custom([[0.5, 0.0], [0.5, 0.5]], [0.25, 0.75],
                                            [0.5, 1.0]),
    'explicit first, unequal': lambda: _custom([[0.0, 0.0, 0.0], [0.2, 0.3, 0.0],
                                                [0.3, 0.3, 0.4]], [0.3, 0.3, 0.4],
                                               [0.0, 0.5, 1.0]),
}


@pytest.mark.parametrize('case', list(CASES))
def test_the_answers_are_the_formulas(case):
    integ = CASES[case]()
    want = _answers(*integ.butcher())
    for _ in range(3):
        assert (integ.is_stiffly_accurate(), integ.stage_structure()) == want
    assert type(integ.is_stiffly_accurate()) is bool
    assert integ.is_fully_implicit() == (want[1] == 'full')


class _Counting:
    """`numpy`, its `allclose` and `triu` counted (the integrator module's
    `np`)."""

    def __init__(self):
        self.calls = 0

    def __getattr__(self, name):
        f = getattr(np, name)
        if name in ('allclose', 'triu'):
            def counted(*a, **k):
                self.calls += 1
                return f(*a, **k)
            return counted
        return f


@pytest.mark.parametrize('make', [RadauIIA3Integrator, TRBDF2Integrator],
                         ids=['radau', 'trbdf2'])
def test_a_run_does_not_classify_its_tableau_again(make, monkeypatch):
    """Twenty fixed steps: the tableau classified at the first (two
    `allclose` for a fully implicit one, three and a `triu` for an ESDIRK),
    not at every step (the parent: those again at each)."""
    counting = _Counting()
    monkeypatch.setattr(integrator, 'np', counting)
    tr = Transient(_mos_chain(), toolkit=circuit.numeric, integrator=make())
    tr.solve(tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
    assert 0 < counting.calls <= 4, counting.calls


def test_a_new_tableau_is_classified_anew():
    integ = RadauIIA3Integrator()
    assert integ.stage_structure() == 'full' and integ.is_stiffly_accurate()
    ## (the instance's own tableau, `butcher`'s cache dropped: a new object)
    integ.A = [[0.5, 0.0], [0.5, 0.5]]
    integ.B = [0.25, 0.75]
    integ.C = [0.5, 1.0]
    assert integ.stage_structure() == 'full', 'the cached tableau stands'
    integ._butcher_cache = None
    assert integ.stage_structure() == 'sdirk'
    assert not integ.is_stiffly_accurate()
