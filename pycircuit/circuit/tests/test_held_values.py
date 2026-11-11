"""THE SAME VALUES REUSE THE NUMERIC (`linearsolver.StructuralKLU`,
csrc/pcx_held.h): a solve whose values are, bit for bit, those the run
numeric was last refactored on makes no refactor -- the same solution --
and counts it in `reused`; a value one ulp off, or a zero of the other sign,
refactors; a fresh factor holds none.  A checkpoint handed back leaves the
held values, the count and the numeric as it found them; one whose numeric
was replaced restores nothing.  A linear circuit's run, where the step's
coefficients repeat, reuses the numeric in C as in Python: the same answers
and counts -- Gear-2 (the stepping loop, its Newton and its error test) and
Radau (the stage loop's stage system)."""
import importlib
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit
from pycircuit.circuit.analysis import jacobian_structure
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.linearsolver import (
    HELD_OK,
    StructuralKLU,
    held_back,
    held_commit,
    held_open,
)
from pycircuit.circuit.tests.test_klu_order import _layers
from pycircuit.circuit.tests.test_structural_klu import _reduced
from pycircuit.circuit.transient import Transient

numeric = circuit.numeric


def _system(k=4):
    """Two meshes joined by 0 V vias (`test_klu_order._layers`), linear:
    the pattern, and the reduced G and C."""
    cir = _layers(k)
    iref = cir.get_node_index(circuit.gnd)
    st = jacobian_structure(cir, iref)
    x = np.zeros(cir.n)
    return st, _reduced(cir.G(x, defaultepar), iref), _reduced(cir.C(x, defaultepar), iref)


def _plain(sk, A, b):
    """`sk`'s solve with nothing held: a refactor, as without the rule."""
    sk._held_st[HELD_OK] = 0
    return sk.solve(A, b, numeric)


def _held(sk):
    return sk._held_v[:sk.nnz].copy()


def test_the_same_values_reuse_the_numeric_bit_for_bit():
    """The run numeric refactored on J, J again: no refactor, the solution
    the refactor's to the bit; one ulp off, or -0.0 where 0.0 stood (an
    always-zero place of the structure: a via's branch diagonal), a
    refactor; held, those values reused in turn."""
    st, G, C = _system()
    sk, plain = StructuralKLU(*st, fallback=None), StructuralKLU(*st, fallback=None)
    rng = np.random.default_rng(3)
    J = G + C / 1e-11
    for _ in range(2):
        b = rng.normal(size=st[0])
        assert sk.solve(J, b, numeric).tobytes() == _plain(plain, J, b).tobytes()
    assert (sk.factors, sk.refactors, sk.reused) == (1, 1, 0)
    for _ in range(3):
        b = rng.normal(size=st[0])
        assert sk.solve(J, b, numeric).tobytes() == _plain(plain, J, b).tobytes()
    assert (sk.factors, sk.refactors, sk.reused) == (1, 1, 3)
    J1 = J.copy()
    i, j = np.argwhere(J1 != 0.0)[0]
    J1[i, j] = np.nextafter(J1[i, j], np.inf)
    b = rng.normal(size=st[0])
    assert sk.solve(J1, b, numeric).tobytes() == _plain(plain, J1, b).tobytes()
    assert (sk.refactors, sk.reused) == (2, 3)
    J2 = J1.copy()
    zeros = [(r, c) for r, c in zip(sk.Ai, np.repeat(np.arange(sk.m), np.diff(sk.Ap)))
             if J2[r, c] == 0.0]
    assert zeros, 'no always-zero place on the pattern'
    J2[zeros[0]] = -0.0
    assert sk.solve(J2, b, numeric).tobytes() == _plain(plain, J2, b).tobytes()
    assert (sk.refactors, sk.reused) == (3, 3)
    assert sk.solve(J2, b, numeric).tobytes() == _plain(plain, J2, b).tobytes()
    assert (sk.refactors, sk.reused) == (3, 4)


def test_a_fresh_factor_holds_no_values():
    """A reuse whose residual check fails falls back to a fresh factor, as
    a refactor's does -- whose numeric is no refactor's: the same values
    next are refactored."""
    st, G, C = _system()
    sk = StructuralKLU(*st, fallback=None)
    J, b = G + C / 1e-11, np.ones(st[0])
    sk.solve(J, b, numeric)
    sk.solve(J, b, numeric)
    assert sk._held_st[HELD_OK] == 1
    sk.REFACTOR_RESIDUAL_TOL = -1.0
    sk.solve(J, b, numeric)
    assert (sk.reused, sk.residual_fallbacks, sk._held_st[HELD_OK]) == (1, 1, 0)
    del sk.REFACTOR_RESIDUAL_TOL
    r = sk.refactors
    sk.solve(J, b, numeric)
    assert (sk.refactors, sk.reused, sk._held_st[HELD_OK]) == (r + 1, 1, 1)


def test_a_hand_back_leaves_the_held_values_and_the_numeric_as_found():
    """A checkpoint opened on J held, a refactor of J' made inside it, then
    handed back: J held again, the count as it was, and the numeric J's --
    J solved next reuses it, the very bytes it gave before; a commit keeps
    J'; one opened inside an open one is not its owner; a numeric replaced
    inside (a fallback's fresh factor) is not restored, and nothing is
    held."""
    st, G, C = _system()
    sk = StructuralKLU(*st, fallback=None)
    J, J1 = G + C / 1e-11, G + C / 2e-11
    b = np.linspace(-1.0, 1.0, st[0])
    sk.solve(J, b, numeric)
    sk.solve(J, b, numeric)
    xa = sk.solve(J, b, numeric)
    held, r, u = _held(sk), sk.refactors, sk.reused
    own = held_open(sk)
    assert own == 1 and held_open(sk) == 0
    sk.solve(J1, b, numeric)
    assert sk.refactors == r + 1 and _held(sk).tobytes() != held.tobytes()
    held_back(sk, own)
    assert sk._held_st[HELD_OK] == 1 and _held(sk).tobytes() == held.tobytes()
    assert sk.reused == u
    ## (the numeric J's again: its reuse passes its residual check -- J1's
    ## would fail it, and the fallback's fresh factor can give the same bytes)
    fb = (sk.factors, sk.residual_fallbacks)
    xb = sk.solve(J, b, numeric)
    assert (sk.refactors, sk.reused) == (r + 1, u + 1)
    assert (sk.factors, sk.residual_fallbacks) == fb
    assert xb.tobytes() == xa.tobytes()
    own = held_open(sk)
    sk.solve(J1, b, numeric)
    held_commit(sk, own)
    assert _held(sk).tobytes() != held.tobytes() and sk._held_st[HELD_OK] == 1
    own = held_open(sk)
    sk.REFACTOR_RESIDUAL_TOL = -1.0
    sk.solve(J, b, numeric)
    del sk.REFACTOR_RESIDUAL_TOL
    f = sk.factors
    held_back(sk, own)
    assert sk._held_st[HELD_OK] == 0 and sk.factors == f


_C_MODULES = ('_tran_newton_c', '_tran_loop_c', '_tran_lte_c', '_dc_newton_c', '_tran_branch_c',
              '_pcnr_c', '_tran_radau_c', '_tran_radau_tc', '_tran_stageloop_c', '_tran_tcloop_c',
              '_tran_ctrl_c', '_tran_pred_c')


def _run(on, monkeypatch, make):
    """`_layers(8)`'s transient, every C path on or off: its answer, its
    statistics and every StructuralKLU's counts -- and the run's counts."""
    for name in _C_MODULES:
        monkeypatch.setattr(importlib.import_module('pycircuit.circuit.' + name), 'ENABLED', on)
    sks, init = [], StructuralKLU.__init__

    def klu_init(self, *a, **k):
        init(self, *a, **k)
        sks.append(self)
    monkeypatch.setattr(StructuralKLU, '__init__', klu_init)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_layers(8), toolkit=numeric, **make)
        res = tr.solve(tend=6e-10, timestep=1e-11)
    st = tr.statistics
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations,
           [(s.factors, s.refactors, s.reused, s.pivot_fallbacks, s.residual_fallbacks)
            for s in sks])
    return run, _paths.since(before)


@pytest.mark.parametrize('method, make, served', [
    ('gear2', {}, 'loop_c:steps'),
    ('radau', {'integrator': RadauIIA3Integrator()}, 'stageloop_c:steps'),
])
def test_a_linear_circuits_run_reuses_the_numeric_in_c_as_in_python(method, make, served,
                                                                     monkeypatch):
    """A linear circuit's run: the same answer, statistics and counts --
    `reused` among them -- with every C path on as with none, the numeric
    reused, and the steps made in C."""
    py, _ = _run(False, monkeypatch, make)
    got, counts = _run(True, monkeypatch, make)
    assert got == py
    assert sum(k[2] for k in got[5]) > 0, got[5]
    assert sum(v for k, v in counts.items() if k.startswith(served)) > 0, counts
