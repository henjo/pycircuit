"""The evaluation core on the circuit's pattern (the C stepping loop plan's
stage 2, 2026-10-10): where a caller's matrices are zero off the pattern --
every entry a matrix pass can write, and the diagonal -- the core zeroes and
forms only the pattern's entries (`PCX_WANT_PATTERN`): G's and C's bins, the
companion `Geq` and `J`.  The same arithmetic on every entry it writes, and
the entries it leaves are the +0.0 the full loops write there, so a run is
the same bytes either way.  Held on the fixed-step chains (the C Newton) and
on a sparse chain (`evaluate` and KLU), adaptive and fixed; and alive: a
pattern missing one entry it should hold makes another run.
"""
import numpy as np
import pytest

from pycircuit.circuit import _tran_core, gnd, numeric
from pycircuit.circuit import elements as E
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.integrator import Gear2Integrator, TrapezoidalIntegrator
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.transient import Transient


def _chain(n=8):
    c = SubCircuit()
    c['vdd'] = E.VS('vdd', gnd, v=1.8)
    c['vg'] = E.VPulse('g', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for k in range(n):
        c[f'R{k}'] = E.R('vdd', f'd{k}', r=5e3)
        c[f'C{k}'] = E.C(f'd{k}', gnd, c=2e-14)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


RUNS = {
    'chain-fixed-gear': (lambda: _chain(30), Gear2Integrator,
                         {'tend': 4e-8, 'timestep': 2e-10, 'fixed_timestep': True}),
    'chain-adaptive-trap': (lambda: _chain(30), TrapezoidalIntegrator,
                            {'tend': 4e-8, 'timestep': 2e-10}),
    ## (a sparse run -- 250 unknowns and more: `evaluate` and KLU, the C
    ## Newton declining)
    'chain-sparse-gear': (lambda: _chain(250), Gear2Integrator, {'tend': 1e-8, 'timestep': 2e-10}),
}


def _run(make, integ, kw, monkeypatch, pattern):
    monkeypatch.setattr(_tran_core, 'PATTERN', pattern)
    tr = Transient(make(), toolkit=numeric, integrator=integ(), reltol=1e-5)
    with quiet(AccuracyWarning):
        res = tr.solve(**kw)
    core = _tran_core.core_for(tr)
    assert core is not None, 'the core does not serve this circuit'
    st = tr.statistics
    return (np.asarray(res.x).tobytes(), st.accepted_steps, st.rejected_steps,
            st.newton_iterations, core.npat)


@pytest.mark.parametrize('name', sorted(RUNS))
def test_a_run_on_the_pattern_is_the_run_on_every_entry(name, monkeypatch):
    make, integ, kw = RUNS[name]
    on = _run(make, integ, kw, monkeypatch, True)
    off = _run(make, integ, kw, monkeypatch, False)
    assert on[4] > 0 and off[4] == 0, 'the pattern was not used, or not switched off'
    assert on[:4] == off[:4]


def test_a_pattern_missing_an_entry_makes_another_run(monkeypatch):
    """Alive: drop the first drain's diagonal entry from the pattern -- its
    `J` entry is then never formed (the zero it started at) -- and the run
    is another (here: a singular J, raised)."""
    make, integ, kw = RUNS['chain-fixed-gear']
    ref = _run(make, integ, kw, monkeypatch, True)
    orig = _tran_core._Core.__init__

    def init(self, cir, *a, **k):
        orig(self, cir, *a, **k)
        pat, n = self.arrays['pat'], self.n
        d = cir.get_node_index('d0')
        at = int(np.nonzero(pat == d * (n + 1))[0][0])
        pat[at] = pat[0] if at else pat[1]             # gone (a duplicate instead)
    monkeypatch.setattr(_tran_core._Core, '__init__', init)
    try:
        got = _run(make, integ, kw, monkeypatch, True)
    except Exception as exc:                            # noqa: BLE001
        got = exc                                       # (a singular J: another run, loudly)
    assert type(got) is not tuple or (got[4] > 0 and got[:4] != ref[:4])
