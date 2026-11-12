"""The transform loop's packed sets (`_tran_tcloop_c.PACKED`, loop.h's
`tc_packed`): on the structural pattern with the real factor on KLU and the
core's pattern, every C and G set the loop holds is the pattern's values in
its order, the fold reads and writes packed, and one dense scratch takes the
C a dense reader needs -- so the loop binds where its dense buffers would
not fit, its memory growing with the pattern, and makes the Python's run."""
import warnings
from types import SimpleNamespace

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_radau_c,
    _tran_radau_tc,
    _tran_tcloop_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.transient import Transient


def _chain(n=250):
    """A pulsed MOS chain with capacitors: at 250 stages 254 unknowns --
    `StructuralKLU`, the core's pattern -- every fifth node with a junction's
    depletion charge, so C moves with the state (C_n not C_n+1)."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-8, tr=1e-8, tf=1e-8, pw=2.4e-7, per=5e-7)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        if k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    return c


def _run(monkeypatch, loops=True, packed='size', budget=None, n=250, tend=3e-7, binds=None,
         no_memo=None):
    """The run's bytes, statistics and both factors' counts, and its paths;
    `budget` the C stage Newton's output budget (`OUT_BYTES`, which the
    dense loop's buffers must fit)."""
    for mod in (_tran_radau_tc, _tran_radau_c, _tran_tcloop_c, _tran_loop_c):
        monkeypatch.setattr(mod, 'ENABLED', loops)
    monkeypatch.setattr(_tran_tcloop_c, 'PACKED', packed)
    if budget is not None:
        monkeypatch.setattr(_tran_radau_c, 'OUT_BYTES', budget)
    if binds is not None:
        bind = _tran_tcloop_c._bind

        def caught(e):
            run_ = bind(e)
            binds.append((e, run_))
            return run_
        monkeypatch.setattr(_tran_tcloop_c, '_bind', caught)
    if no_memo is not None:
        load = _tran_tcloop_c._load_x

        def blind(e, run_, xl):
            ## (the entry finds no record of x_n: the C makes its passes)
            no_memo.append(1)
            held = e.tr.__dict__.pop('_dev_memo', None)
            try:
                return load(e, run_, xl)
            finally:
                if held is not None:
                    e.tr._dev_memo = held
        monkeypatch.setattr(_tran_tcloop_c, '_load_x', blind)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_chain(n), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
        res = tr.solve(tend=tend, timestep=1e-8)
    st = tr.statistics
    zs = tr.__dict__.get('_radau_zsolver')
    sk = tr._get_linearsolver()._choice
    key = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations,
           None if zs is None else (zs.analyses, zs.factors, zs.refactors, zs.pivot_fallbacks),
           *(getattr(sk, k, None) for k in ('factors', 'refactors', 'pivot_fallbacks',
                                            'residual_fallbacks')))
    return key, _paths.since(before)


def test_where_the_dense_buffers_do_not_fit_the_packed_sets_serve(monkeypatch):
    """With the budget below the dense buffers' (`_fits` false), 'size'
    binds the packed sets and the loop makes the Python's run; 'never'
    declines for the size, as before."""
    py, _d = _run(monkeypatch, loops=False)
    got, d = _run(monkeypatch, budget=1)
    assert got == py
    assert d['tcloop_c:packed'] >= 1 and d['tcloop_c:attempts'] > 20
    got, d = _run(monkeypatch, budget=1, packed='never')
    assert got == py
    assert d.get('tcloop_c:packed', 0) == 0 and d['tcloop_c:size'] >= 1


def test_always_packed_is_the_dense_loops_run(monkeypatch):
    """'always' runs the packed sets where the dense ones fit too: the run
    the dense loop's and the Python's, each loop serving its share."""
    py, _d = _run(monkeypatch, loops=False)
    dense, dd = _run(monkeypatch)
    packed, dp = _run(monkeypatch, packed='always')
    assert dense == packed == py
    assert dd.get('tcloop_c:packed', 0) == 0 and dp['tcloop_c:packed'] >= 1
    assert dd['tcloop_c:attempts'] > 20 and dp['tcloop_c:attempts'] > 20


def test_the_loops_own_passes_at_the_start(monkeypatch):
    """Where an entry finds no record of x_n in the memo, the loop makes x_n's
    passes itself -- C made in the scratch and gathered, G gathered from the
    bin -- and the run is the Python's."""
    py, _d = _run(monkeypatch, loops=False)
    seen = []
    got, d = _run(monkeypatch, packed='always', no_memo=seen)
    assert got == py
    assert seen and d['tcloop_c:packed'] >= 1 and d['tcloop_c:attempts'] > 20


def _bytes(run_):
    """The loop's bound buffers' bytes, the scratch and the fresh pages'
    work (never or rarely written) aside."""
    return sum(a.nbytes for k, a in run_.tbuf.items()
               if isinstance(a, np.ndarray) and k != 'cs')


def test_the_packed_buffers_grow_with_the_pattern(monkeypatch):
    """Chains of 250 and 500 stages (`StructuralKLU`) on the packed sets:
    the loop's buffers about double (the dense ones grow fourfold), and the
    end sets and x_n's hold the pattern's values alone; the scratch is the
    one n x n."""
    sizes = []
    for n in (250, 500):
        binds = []
        _run(monkeypatch, packed='always', n=n, tend=1e-7, binds=binds)
        e, run_ = [b for b in binds if b[0].packed][-1]
        npat = int(e.rec.core.npat)
        assert run_.tbuf['xr_C'].shape == (npat,) and run_.tbuf['te02'].shape == (npat,)
        assert run_.tbuf['cs'].shape == (e.n * e.n,)
        assert run_.tbuf['fold_rf'].shape == (int(run_.fo.pn) + 1,)
        sizes.append(_bytes(run_))
    assert sizes[1] < 2.3 * sizes[0]


def _entry(memo_C, n=4):
    """A packed entry's `_load_x` inputs: a 4 x 4 core pattern (the
    diagonal and (0, 1)), the memo's record at x_n."""
    pat = np.array([0, 1, 5, 10, 15], dtype=np.int64)
    xl = np.arange(n, dtype=float)
    rec = {'q': np.ones(n), 'i': np.full(n, 2.0), 'C': memo_C,
           'G': np.diag(np.arange(1.0, n + 1.0))}
    e = SimpleNamespace(packed=True, rec=SimpleNamespace(core=SimpleNamespace(arrays={'pat': pat})),
                        tr=SimpleNamespace(_dev_memo=({xl.tobytes(): rec}, {})))
    tb = {'xr_q': np.zeros(n), 'xr_i': np.zeros(n), 'xr_C': np.zeros(5), 'xr_G': np.zeros(5)}
    return e, SimpleNamespace(s=SimpleNamespace(), tbuf=tb), xl


def test_the_memos_sets_are_loaded_at_the_pattern():
    """x_n's C and G from the memo: the pattern's values in its order; a
    memo matrix with a nonzero off the pattern is no packed set (False:
    the entry declined)."""
    Cm = np.zeros((4, 4))
    Cm[0, 0], Cm[0, 1], Cm[3, 3] = 1.0, -2.0, 3.0
    e, run_, xl = _entry(Cm)
    assert _tran_tcloop_c._load_x(e, run_, xl) is True
    assert np.array_equal(run_.tbuf['xr_C'], [1.0, -2.0, 0.0, 0.0, 3.0])
    assert np.array_equal(run_.tbuf['xr_G'], [1.0, 0.0, 2.0, 3.0, 4.0])
    assert run_.s.te_x == -1 and run_.s.xfields
    off = Cm.copy()
    off[2, 0] = 1e-30
    e, run_, xl = _entry(off)
    assert _tran_tcloop_c._load_x(e, run_, xl) is False


def test_the_real_factors_take_maps_onto_the_fold():
    """The analysis solver's positions in the reduced factor become their
    places among the fold's; one off the core's pattern the zero slot."""
    run_ = SimpleNamespace(tbuf={'fold_rmaj': np.array([4, 0, 7, 2], dtype=np.int64)})
    e = SimpleNamespace(packed=True)
    sk = SimpleNamespace(take=np.array([0, 2, 3, 4, 7], dtype=np.int64))
    assert np.array_equal(_tran_tcloop_c._rk_take(run_, e, sk), [1, 3, 4, 0, 2])
    e.packed = False
    assert np.array_equal(_tran_tcloop_c._rk_take(run_, e, sk), sk.take)


def test_zero_pages_are_zero_and_writable():
    a = _tran_radau_tc.zero_pages(300_000)
    assert a.shape == (300_000,) and a.dtype == np.float64 and not a.any()
    a[12345] = 1.5
    assert a[12345] == 1.5 and np.count_nonzero(a) == 1
    b = _tran_radau_tc.zero_pages(10, np.int32)
    assert b.dtype == np.int32 and b.shape == (10,) and not b.any()


@pytest.mark.parametrize('no_memo', [False, True], ids=['memo', 'no-memo'])
def test_the_packed_sets_are_the_packed_loops_run(monkeypatch, no_memo):
    """THE PACKED SETS (`SETS`, from `SETS_FROM` unknowns: here from any) on
    the packed transform: its core calls write C packed and read G from the
    core's packed bin, its screen reads the packed C -- the run the bin's
    and the Python's, every count the bin's, the loop serving; with x_n's
    passes the loop's own (no memo record) too."""
    from pycircuit.circuit import _tran_radau
    monkeypatch.setattr(_tran_tcloop_c, 'SETS_FROM', 1)
    ## (the estimate on C's pattern, as from 1 000 unknowns: no dense reader of C)
    monkeypatch.setattr(_tran_radau, 'PATTERN_FROM', 1)
    py, _d = _run(monkeypatch, loops=False)
    with monkeypatch.context() as m:
        m.setattr(_tran_tcloop_c, 'SETS', False)
        binds_off = []
        off, doff = _run(m, packed='always', binds=binds_off, no_memo=[] if no_memo else None)
    binds_on = []
    on, don = _run(monkeypatch, packed='always', binds=binds_on, no_memo=[] if no_memo else None)
    assert on == off == py
    assert don == doff
    assert don['tcloop_c:packed'] >= 1 and don['tcloop_c:attempts'] > 20
    assert all(int(r.s.tc_sets) == 1 for _e, r in binds_on)
    assert all(int(r.s.tc_sets) == 0 for _e, r in binds_off)
