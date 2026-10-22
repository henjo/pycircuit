"""A fully implicit stage method's coupled stage block factored with KLU
on its structural pattern (`_InnerTransient._step_klu`) from
`SPARSE_BLOCK_MIN_N` unknowns; the answers the dense LU's to the
rounding."""
import numpy as np

from pycircuit.circuit import _paths, circuit, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import C, Diode, R, VSin
from pycircuit.circuit.linearsolver import DenseSolver, StructuralKLU
from pycircuit.circuit.shooting import PAC, PSS, _pss_inner
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet

#: sections enough for Radau's `3 m` block to reach the floor (`m =
#: sections + 2`)
SECTIONS = _pss_inner.SPARSE_BLOCK_MIN_N // 3


def _ladder(sections, diode=True):
    """A sine into an RC ladder of `sections`, a diode at its end:
    ``sections + 2`` unknowns."""
    c = SubCircuit()
    c['vs'] = VSin('n0', gnd, vac=1.0, va=1.0, freq=1e3)
    for k in range(sections):
        c[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=100.0)
        c[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-8)
    if diode:
        c['D'] = Diode(f'n{sections}', gnd)
    return c


def _pss(cir, klu, monkeypatch, npts=20, **kw):
    circuit.default_toolkit = circuit.numeric
    monkeypatch.setattr(_pss_inner, 'SHOOT_KLU', klu)
    pss = PSS(cir, method='radau', reltol=1e-10, **kw)
    with quiet(AccuracyWarning):
        res = pss.solve(period=1e-3, timestep=1e-3 / npts)
    assert pss.converged
    out = 'n{}'.format(len([e for e in cir.elements if e.startswith('R')]))
    return pss, np.asarray(res['tpss'].v(out).get_y(), dtype=float)


def test_a_large_stage_block_factors_on_its_structure_and_the_c_step_declines(monkeypatch):
    """Radau's `3 m` block from `SPARSE_BLOCK_MIN_N`, on a linear ladder
    the C step map (`_sens_c`) would otherwise take: it declines the block
    ('pss.sens:sparse') -- the Python step factors it on its structure, and
    the C's dense LU is not that step's twin.  (Without it the analysis's
    solver chose on the first matrix it saw, here SuperLU on the block.)"""
    out = {}
    for klu in (False, True):
        before = _paths.snapshot()
        pss, y = _pss(_ladder(SECTIONS, diode=False), klu, monkeypatch)
        out[klu] = (y, _paths.since(before), pss)
    (y0, d0, _), (y1, d1, pss1) = out[False], out[True]
    n = 3 * (pss1.cir.n - 1)
    assert n >= _pss_inner.SPARSE_BLOCK_MIN_N
    assert not d0.get('pss.sens:sparse')
    assert d1.get('pss.sens:sparse', 0) > 0 and not d1.get('pss.sens:served')
    klu = pss1._step_klus[n][1]
    assert type(klu) is StructuralKLU and klu.numerics > 0 and klu.foreign == 0
    assert np.allclose(y1, y0, rtol=1e-9, atol=1e-12)


def test_a_step_jacobian_and_a_small_block_stay_with_the_caller_s_solver(monkeypatch):
    monkeypatch.setattr(_pss_inner, 'SHOOT_KLU', True)
    big = PSS(_ladder(SECTIONS), method='radau')
    m = big.cir.n - 1
    assert big._step_klu(m) is None and big._step_klu(3 * m) is not None
    small = PSS(_ladder(SECTIONS - 3), method='radau')
    ms = small.cir.n - 1
    assert 3 * ms < _pss_inner.SPARSE_BLOCK_MIN_N and small._step_klu(3 * ms) is None
    own = PSS(_ladder(SECTIONS), method='radau', linearsolver=DenseSolver())
    assert own._step_klu(3 * m) is None


def test_a_factored_period_on_the_structure_answers_as_dense(monkeypatch):
    """The kept walk's block factors (`_factorise`) and their transposed
    solves -- PAC's adjoint row and its sideband response -- on the
    structure answer as without it."""
    got = {}
    for klu in (False, True):
        pss, _y = _pss(_ladder(SECTIONS), klu, monkeypatch)
        pac = PAC(pss.cir, toolkit=circuit.numeric)
        with quiet(AccuracyWarning):
            row, _info = pac.adjoint_transfer_row(pss, 300.0, f'n{SECTIONS}')
            res = pac.solve(pss, [300.0])
        got[klu] = (np.asarray(row), np.asarray(res.v(f'n{SECTIONS}', gnd)))
        if klu:
            assert pss._step_klus[3 * (pss.cir.n - 1)][1].numerics > 0
    for a, b in zip(got[True], got[False]):
        assert np.allclose(a, b, rtol=1e-8, atol=1e-14 * np.max(np.abs(b)))
