"""A plain run's coupled stage block on a StructuralKLU (`_run_stage_klu`,
roadmap 202): a fully implicit tableau from `MIN_N_FOR_SPARSE` block
unknowns solves on the shooting walk's tiled structural pattern -- its
answers the dense LU's to the conditioning-amplified rounding, the C and
the Python stage Newton the same bits on it -- and nothing else does: a
small circuit, a given `linearsolver=`, the switch off, a caller's own."""
import warnings

import numpy as np

from pycircuit.circuit import _paths, _tran_radau, _tran_radau_c, circuit
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.linearsolver import MIN_N_FOR_SPARSE, DenseSolver, StructuralKLU
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.transient import Transient

#: a chain whose 3 m stage block is past `MIN_N_FOR_SPARSE`; its compiled
#: Jacobian would turn 'auto''s cost transform on
STAGES = 84
KW = {'tend': 1e-7, 'timestep': 2e-9}


def _tr(build=lambda: _mos_chain(n=STAGES, cap=1e-13), **make):
    return Transient(build(), toolkit=circuit.numeric, integrator=make.pop('integ', None)
                     or RadauIIA3Integrator(), radau_transform=False, **make)


def _run(tr, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        before = _paths.snapshot()
        res = tr.solve(**dict(KW, **kw))
        return res, _paths.since(before)


def test_a_large_block_solves_on_its_klu_near_the_dense_lu(monkeypatch):
    """The chain's 3 m block past `MIN_N_FOR_SPARSE` is solved on the run's
    StructuralKLU by the C stage Newton; on a fixed grid the waveform is the
    dense run's within the stage solves' amplified rounding, and an
    adaptive run takes the same steps."""
    tr = _tr()
    assert 3 * (tr.cir.n - 1) >= MIN_N_FOR_SPARSE
    out = {}
    for on in (True, False):
        monkeypatch.setattr(_tran_radau, 'RADAU_KLU', on)
        t = _tr()
        res, d = _run(t, fixed_timestep=True)
        out[on] = (np.asarray(res.x), d, t)
        a, _da = _run(_tr())
        out[on] += (a,)
    (xk, dk, tk, ak), (xd, dd, td, ad) = out[True], out[False]
    assert type(tk._stage_klu) is StructuralKLU and '_stage_klu' not in td.__dict__
    assert dk.get('radau_c:on_klu', 0) > 0 and not dd.get('radau_c:on_klu'), dk
    assert xk.shape == xd.shape
    assert np.max(np.abs(xk - xd)) < 1e-9 * max(1.0, np.max(np.abs(xd)))
    assert len(ak.sweep_values) == len(ad.sweep_values)


def test_the_c_and_the_python_stage_newton_are_one_on_the_runs_klu(monkeypatch):
    """The run's StructuralKLU is the C stage Newton's and the Python's: the
    same bits with the C off."""
    out = {}
    for on in (True, False):
        monkeypatch.setattr(_tran_radau_c, 'ENABLED', on)
        res, d = _run(_tr())
        out[on] = (np.asarray(res.x).tobytes(), np.asarray(res.sweep_values).tobytes(), d)
    assert out[True][:2] == out[False][:2]
    assert out[True][2].get('radau_c:on_klu', 0) > 0 and not out[False][2].get('radau_c:on_klu')


def test_only_a_large_fully_implicit_plain_run_takes_it(monkeypatch):
    """None for a small block, a multistep method, a given `linearsolver=`,
    the switch off; a caller's own `_stage_klu` kept; one per topology
    and size, kept from run to run."""
    small = _tr(lambda: _mos_chain(cap=1e-13))
    _run(small)
    assert '_stage_klu' not in small.__dict__
    gear = _tr(integ=Gear2Integrator())
    _run(gear)
    assert '_stage_klu' not in gear.__dict__
    given = _tr(linearsolver=DenseSolver())
    _run(given)
    assert '_stage_klu' not in given.__dict__
    tr = _tr()
    own = StructuralKLU(3, np.arange(3), np.arange(3), fallback=DenseSolver)
    tr._stage_klu = own
    tr._run_stage_klu()
    assert tr._stage_klu is own
    tr = _tr()
    _run(tr)
    first = tr._stage_klu
    _run(tr)
    assert tr._stage_klu is first
    monkeypatch.setattr(_tran_radau, 'RADAU_KLU', False)
    _run(tr)
    assert '_stage_klu' not in tr.__dict__
