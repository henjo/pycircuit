"""A plain Radau transient's Python around the C stage Newton, trimmed
(roadmap 201): the embedded estimate reads x_n's C, G and i from the device
memo and the run's stateful limiters; the stage seeds come from pred.c; the
three screens run in C on the core's pattern; the source vectors come from
the C sources' table where the direct pass's stamp holds it.  Each the same
run, bit for bit -- waveform, statistics, warnings, the device memo -- as
without it."""
import sys
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_branch_c,
    _tran_pred_c,
    _tran_radau,
    _tran_radau_c,
    _tran_sources,
    circuit,
)
from pycircuit.circuit.elements import VPWL, C, R
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_radau_c import _memo
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_radau_c.driver() is None,
                                reason='the radau C did not build')

#: an adaptive run: the embedded estimate made every attempt
KW = {'tend': 4e-7, 'timestep': 2e-8}


def _chain():
    return _mos_chain(cap=1e-13)


def _long_chain():
    """A chain whose `C` the C screen takes only on the core's pattern: 68
    unknowns (past `_cbuild.INPLACE_MIN_ENTRIES` dense), its pattern sparse."""
    return _mos_chain(n=64, cap=1e-13)


def _run(build, patch=None, cir=None, make=None, **kw):
    """`(the run's answer, the paths it took)`, `patch(mp)` applied first;
    `cir` a circuit to run instead of a built one, `make` the transient's
    options."""
    with pytest.MonkeyPatch.context() as mp:
        if patch is not None:
            patch(mp)
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build() if cir is None else cir, toolkit=circuit.numeric,
                           integrator=RadauIIA3Integrator(), **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**dict(KW, **kw))
            d = _paths.since(before)
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _memo(tr)), d


def _direct(cir):
    """The circuit's direct source pass program (`_stamp_plan._SourceDirect`)."""
    return cir.__dict__['_stamp_plan'].methods['src:u'].direct


def _chain_run_table():
    """A chain run with the C sources on: its circuit."""
    cir = _chain()
    _run(None, cir=cir)
    return cir


def test_the_estimate_reads_x_n_from_the_memo():
    """Where the step that reached x_n recorded it, the estimate takes its
    C, G and i from the device memo -- the core's values, as its passes
    make them again."""
    def unread(mp):
        get = Transient._memo_get

        def memo_get(self, x):
            if sys._getframe(1).f_code.co_name == '_radau_error_estimate':
                return None
            return get(self, x)
        mp.setattr(Transient, '_memo_get', memo_get)
    a, da = _run(_chain)
    b, db = _run(_chain, unread)
    assert a == b
    assert da.get('radau.est:memo', 0) > 0 and not db.get('radau.est:memo'), (da, db)


def test_the_estimate_reads_the_runs_stateful_limiters(monkeypatch):
    """The estimate takes the run's stateful limiters (`_stateful_lims`,
    `_begin_run`'s) as the stage solver does, rather than walking the
    circuit's elements every attempt."""
    a, _da = _run(_chain)
    calls = []
    orig = _tran_radau.stateful_limiters

    def counted(cir):
        calls.append(1)
        return orig(cir)
    monkeypatch.setattr(_tran_radau, 'stateful_limiters', counted)
    b, _db = _run(_chain)
    assert calls == [] and a == b


@pytest.mark.parametrize('unit', ['pred', 'screen', 'sources'])
def test_each_trim_is_the_pythons_run(unit):
    """The seeds from pred.c, the screens in C on the core's pattern, the
    source vectors from the C table: switched off, the same run."""
    mod = {'pred': _tran_pred_c, 'screen': _tran_branch_c, 'sources': _tran_sources}[unit]
    ## (the long chain's compiled Jacobian turns 'auto''s cost transform on)
    build, kw = ((_long_chain, {'tend': 1e-7, 'make': {'radau_transform': False}})
                 if unit == 'screen' else (_chain, {}))
    a, da = _run(build, **kw)
    b, db = _run(build, lambda mp: mp.setattr(mod, 'ENABLED', False), **kw)
    assert a == b
    assert da.get('pred:fast') == db.get('pred:fast') and da['u:called'] == db['u:called']
    if unit == 'pred':
        assert da.get('pred:fast', 0) > 0 and not any(k.startswith('pred_c:') for k in da), da
    elif unit == 'screen':
        att = da['radau_c:served']
        assert da.get('branch_c:pattern', 0) == 3 * att and not db.get('branch_c:pattern'), da
    else:
        ## (made and called quietly: the program holds it)
        assert type(_direct(_chain_run_table()).tab) is _tran_sources.Table


def test_a_source_written_between_runs_is_read_again():
    """The C table kept on the direct pass's program is checked again when
    a source's value or its time function's attribute is written: each run
    the C-free run's."""
    def runs(enabled):
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(_tran_sources, 'ENABLED', enabled)
            cir = _chain()
            out = [_run(None, cir=cir)[0]]
            cir['vdd'].ipar.v = 1.7
            out.append(_run(None, cir=cir)[0])
            cir['vg'].function.amplitude = 0.4
            out.append(_run(None, cir=cir)[0])
            return out
    a, b = runs(True), runs(False)
    assert a == b and a[0] != a[1] != a[2]


def test_a_pwl_source_is_the_pythons():
    """A PWL's points, which no dict watcher sees, keep the Python's pass."""
    def build():
        c = circuit.SubCircuit()
        c['vin'] = VPWL('a', circuit.gnd, tvpairs=[0.0, 0.0, 1e-7, 1.0, 3e-7, 0.2])
        c['r'] = R('a', 'b', r=1e3)
        c['c'] = C('b', circuit.gnd, c=1e-10)
        return c
    cir = build()
    a, _da = _run(None, cir=cir)
    b, _db = _run(build, lambda mp: mp.setattr(_tran_sources, 'ENABLED', False))
    assert a == b and _direct(cir).tab is False


def test_pred_c_seeds_through_an_underflow():
    """Seeds whose fit underflows: where numpy ignores underflow pred.c makes
    them (no hand-back) and they are the Python's; where numpy would report
    it, pred.c hands them back."""
    from pycircuit.circuit.tests.test_pred_c import _P
    p = _P()
    p._pred_hist = [(k * 1e-9, np.array([1e-308 * (1 + 0.3 * k), 1.0 + 0.1 * k]))
                    for k in range(4)]
    x = p._pred_hist[-1][1]
    ts = [3.4e-9, 3.7e-9, 4e-9]
    py = [p._predict_fast(t, 3) for t in ts]
    before = _paths.snapshot()
    got = _tran_pred_c.stages(p, p._pred_hist, x, ts, 3)
    assert got is not None, _paths.since(before)
    assert [g.tobytes() for g in got] == [np.asarray(q).tobytes() for q in py]
    with np.errstate(under='warn'):
        before = _paths.snapshot()
        assert _tran_pred_c.stages(p, p._pred_hist, x, ts, 3) is None
        assert _paths.since(before).get('pred_c:bail:flags', 0) == 1
