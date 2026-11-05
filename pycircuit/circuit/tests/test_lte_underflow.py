"""The multistep error test's underflow rule (lte.c's `lte_flags`, the
struct's `under_ok`): J^-1 Eg decays down a long chain -- an inverter at a
rail passes ~1e-3 a stage -- to subnormal entries, and its products, its
quotients and the refactor's validation raise UNDERFLOW.  Where numpy's
error state ignores underflow (its default) that is no warning and the
values are the IEEE ones, so the C test -- the C loop's, and the one served
from Python -- keeps the attempt where it handed every such one back; where
numpy would warn it still hands back.  Either way the run is the Python's."""
import warnings

import numpy as np

from pycircuit.circuit import _paths, _tran_loop_c, _tran_lte_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.transient import Transient


def _chain(n=400):
    """A pulsed MOS chain with capacitors: at 400 stages its error test's
    J^-1 Eg reaches the subnormal range on a few attempts (`StructuralKLU`,
    gear2)."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def _ladder(n=60):
    """A linear ladder that attenuates 1e-6 a stage (1 MOhm in series, 1 Ohm
    and 1 pF to ground): past ~50 stages its node charges are subnormal, so
    the charge's divided difference underflows as well as J^-1 Eg."""
    c = SubCircuit()
    c.add_node('n0')
    c['vin'] = VPulse('n0', gnd, v1=0.0, v2=1.0, td=5e-9, tr=2e-9, tf=2e-9, pw=4e-8, per=1e-7)
    for k in range(n):
        c.add_node(f'n{k + 1}')
        c[f'rs{k}'] = R(f'n{k}', f'n{k + 1}', r=1e6)
        c[f'rp{k}'] = R(f'n{k + 1}', gnd, r=1.0)
        c[f'c{k}'] = C(f'n{k + 1}', gnd, c=1e-12)
    return c


def _run(monkeypatch, loop=True, lte=True, under='ignore', make=_chain):
    """The run's bytes and statistics, and the paths it took."""
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte)
    before = _paths.snapshot()
    with warnings.catch_warnings(), np.errstate(under=under):
        warnings.simplefilter('ignore')
        tr = Transient(make(), toolkit=circuit.numeric)
        res = tr.solve(tend=1e-7, timestep=1e-9)
    st = tr.statistics
    key = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations)
    return key, _paths.since(before)


def test_the_loop_keeps_an_attempt_whose_error_test_underflows(monkeypatch):
    """Under numpy's default the C loop makes every attempt but those after
    its Newton's hand-backs -- no error test handed back, none sent to the
    dense validation -- and the run is the Python's (no C loop, no C error
    test)."""
    py, _d = _run(monkeypatch, loop=False, lte=False)
    got, d = _run(monkeypatch)
    assert got == py
    assert d.get('loop_c:bail:lte', 0) == 0 and d.get('loop_c:lte_dense', 0) == 0
    assert d['loop_c:attempts'] >= 90


def test_where_numpy_would_warn_the_loop_hands_it_back(monkeypatch):
    """Under `np.errstate(under='warn')` the same attempts are the Python's
    -- the loop's error test hands them back (an underflow its flag) -- and
    the run is the same."""
    py, _d = _run(monkeypatch, loop=False, lte=False)
    got, d = _run(monkeypatch, under='warn')
    assert got == py
    assert d.get('loop_c:bail:lte', 0) >= 1


def test_the_error_test_served_from_python_keeps_it_too(monkeypatch):
    """Without the C loop every attempt's error test is served from Python
    (`_tran_lte_c`): under numpy's default none is handed back on a flag,
    under a warning state some are -- each run the Python's."""
    py, _d = _run(monkeypatch, loop=False, lte=False)
    got, d = _run(monkeypatch, loop=False)
    assert got == py
    assert d.get('lte_c:bail:flags', 0) == 0 and d['lte_c:served'] >= 100
    got, d = _run(monkeypatch, loop=False, under='warn')
    assert got == py
    assert d.get('lte_c:bail:flags', 0) >= 1


def test_a_charge_in_the_subnormal_range(monkeypatch):
    """The attenuating ladder, its error test served from Python: under
    numpy's default no attempt is handed back on a flag -- its charges'
    divided differences underflow too -- under a warning state most are;
    each run the Python's."""
    py, _d = _run(monkeypatch, loop=False, lte=False, make=_ladder)
    got, d = _run(monkeypatch, loop=False, make=_ladder)
    assert got == py
    assert d.get('lte_c:bail:flags', 0) == 0 and d['lte_c:served'] >= 60
    got, d = _run(monkeypatch, loop=False, under='warn', make=_ladder)
    assert got == py
    assert d.get('lte_c:bail:flags', 0) >= 20
