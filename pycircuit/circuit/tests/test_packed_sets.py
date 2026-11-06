"""THE PACKED SETS (newton.h, `_tran_newton_c.SETS`; plan 261's L4d2,
roadmap 270): from `SETS_FROM` unknowns, where the Newton solves on the
StructuralKLU and the core holds a pattern, the C Newton's C, Geq and J sets,
the C loop's cached C and Geq slots, and what the error test and the branch
screen read of them hold the pattern's values alone (npat each) -- the core
asked for them packed (`PCX_WANT_SETS`) -- where n x n sets, written at a few
entries a row, were resident whole under huge pages.  Every operation on
them was already the pattern's: the run is the n x n run's bit for bit, and
what it leaves the Python (the write-backs: `_C_cache`, `_Cmat`, `_Geq`)
n x n, scattered on pattern pages."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_newton_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import Gear2Integrator, TrapezoidalIntegrator
from pycircuit.circuit.transient import Transient


def _chain(n):
    """A pulsed MOS chain: a capacitor to ground, a gate-drain one at every
    stage, one three stages back at every third, a junction's depletion
    charge at every fifth node (C moving with the state)."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=2e-8, per=5e-8)
    for k in range(n):
        g = 'g0' if k == 0 else f'd{k - 1}'
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'cm{k}'] = C(f'd{k}', g, c=2e-15)
        if k >= 3 and k % 3 == 0:
            c[f'cl{k}'] = C(f'd{k}', f'd{k - 3}', c=1e-15)
        if k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    return c


def _b(a):
    return np.asarray(a, dtype=float).tobytes()


def _run(monkeypatch, sets, integ=Gear2Integrator, n=300, **kw):
    """The chain's run with THE PACKED SETS on from 100 unknowns, or off:
    its bytes and statistics, the solver's counts, the state it leaves the
    Python, the record's sets and their size, its paths."""
    monkeypatch.setattr(_tran_newton_c, 'SETS', sets)
    monkeypatch.setattr(_tran_newton_c, 'SETS_FROM', 100)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(n), toolkit=circuit.numeric, integrator=integ())
        res = tr.solve(**(kw or {'tend': 6e-8, 'timestep': 1e-9}))
    st, ch, rec = tr.statistics, tr._get_linearsolver()._choice, tr.__dict__['_newton_c']
    run = (_b(res.x), _b(res.sweep_values), st.accepted_steps, st.rejected_steps,
           st.newton_iterations, st.chord_fallbacks, st.early_aborts)
    klu = tuple(getattr(ch, k, None) for k in ('analyses', 'factors', 'refactors', 'chord_factors',
                                          'chord_refactors', 'residual_fallbacks', 'pivot_fallbacks'))
    left = (_b(tr._Cmat), tr._Cmat.shape, _b(tr._Geq), tr._Geq.shape, _b(tr._C_cache[1]),
            _b(tr._q_cache[1]), _b(tr._iq))
    ## (the process's once-only counts left out: the first run makes them)
    paths = {k: v for k, v in _paths.since(before).items() if not k.startswith('once:')}
    return run, klu, left, (int(rec.s.sets), rec.buf['Cj'].size, rec.core.npat, rec.n), paths


@pytest.fixture(autouse=True)
def _warm(monkeypatch):
    """A first run of each kind for the process's builds made once (the
    sources' unit among them: a first run's loop declines until it is)."""
    for kw in ({'tend': 4e-9, 'timestep': 1e-9},
               {'tend': 4e-9, 'timestep': 2e-10, 'fixed_timestep': True}):
        with monkeypatch.context() as m:
            _run(m, False, n=20, **kw)


@pytest.mark.parametrize('integ, kw', [
    (Gear2Integrator, {}), (TrapezoidalIntegrator, {}),
    (Gear2Integrator, {'tend': 3e-8, 'timestep': 2e-10, 'fixed_timestep': True})],
    ids=['gear2', 'trap', 'gear2-fixed'])
def test_a_run_on_the_packed_sets_is_the_n_x_n_run(monkeypatch, integ, kw):
    """A 305-unknown chain (the StructuralKLU, the core's pattern): with the
    sets packed the run, its statistics, the solver's counts and every path
    count are the n x n run's; the state left the Python n x n and the
    same bytes; the record's sets npat values each, the C loop serving."""
    a = _run(monkeypatch, False, integ, **kw)
    b = _run(monkeypatch, True, integ, **kw)
    assert a[0] == b[0] and a[1] == b[1] and a[2] == b[2]
    assert a[4] == b[4]
    assert a[3][0] == 0 and a[3][1] == a[3][3] ** 2
    assert b[3][0] == 1 and b[3][1] == b[3][2] > 0
    assert b[4].get('loop_c:attempts', 0) > 20 and b[4].get('newton_c:served', 0) > 20


def test_under_the_size_the_sets_stay_n_x_n(monkeypatch):
    """Below `SETS_FROM` (its default 1 024) the record's sets are n x n."""
    monkeypatch.setattr(_tran_newton_c, 'SETS', True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(300), toolkit=circuit.numeric)
        tr.solve(tend=4e-9, timestep=1e-9)
    rec = tr.__dict__['_newton_c']
    assert rec.s.sets == 0 and rec.buf['Cj'].shape == (rec.n, rec.n) and rec.n > 300


def test_a_shooting_walk_on_packed_sets_is_the_pythons(monkeypatch):
    """A PSS whose transient's sets are packed: the C walk declines
    (`pss.lmmwalk:sets` -- its map reads J n x n) and the PSS is the one
    the sets n x n make."""
    from pycircuit.circuit.tests import test_lmmwalk_c as W
    monkeypatch.setattr(_tran_newton_c, 'SETS_FROM', 10)
    monkeypatch.setattr(_tran_newton_c, 'SETS', False)
    want, _c = W._run(monkeypatch, 'gear', True, n=40)
    monkeypatch.setattr(_tran_newton_c, 'SETS', True)
    got, counts = W._run(monkeypatch, 'gear', True, n=40)
    assert counts.get('pss.lmmwalk:sets', 0) > 0, counts
    assert [k for k in want if want[k] != got[k]] == []
