"""The branch check's cheap screen in C (`_tran_branch_c`, the C stepping
loop plan's stage 3, step d; 2026-10-11): bit for bit
`_BranchCheck._branch_screen_C` up to its SVD -- the verdict, |C|'s largest
entry and the running maximum -- on drawn matrices, read whole or on the
pattern the core wrote; whole runs the same with it on and off; its
declines; and alive: a pattern missing the largest entry is another screen.
"""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_branch_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit._tran_branch import _BranchCheck
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import Gear2Integrator, TrapezoidalIntegrator
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_branch_c.driver() is None,
                                reason='the branch screen unit did not build')

_SVD = 'svd'


class _Screen(_BranchCheck):
    """The screen alone: its tolerance and running maximum, the SVD stubbed
    (the C hands that level back to the Python's own)."""

    def _branch_svd(self, C, r0, ref):
        return _SVD


def _python(C, cmax, g=None):
    t = _Screen()
    if cmax is not None:
        t._branch_cmax = cmax
    if g is not None:
        t.irefnode = g
    verdict, nrm = t._branch_screen_C(C, C.shape[0])
    return verdict, nrm, t._branch_cmax


def _c(C, cmax, pat=None, g=None):
    t = _Screen()
    if cmax is not None:
        t._branch_cmax = cmax
    if g is not None:
        t.irefnode = g
    if pat is not None:
        t._C_pattern = (C, pat)
    got = _tran_branch_c.screen(t, C)
    if got is None:
        return None
    code, nrm, ref = got
    verdict = ((False, None) if code == _tran_branch_c.QUIET
               else (True, None) if code == _tran_branch_c.COLLAPSED else _SVD)
    return verdict, nrm, ref


## (magnitudes about the screen's thresholds: its tolerance 1e-9 and the
## positive minimum's 1e-6 of the running maximum, subnormals, signed zeros)
_entry = st.one_of(st.floats(-2.0, 2.0), st.sampled_from(
    [0.0, -0.0, 5e-324, 1e-300, 1e-12, 1e-9, 1.0000000001e-9, 1e-6, 9.99e-7, 1e-3, 1e300]))
_cmax = st.one_of(st.none(), st.sampled_from([0.0, 1e-12, 1.0, 1e3, 1e300]),
                  st.floats(0.0, 10.0))


@settings(deadline=None)
@given(data=st.data(), n=st.integers(1, 8), cmax=_cmax)
def test_the_screen_is_the_python_one_bit_for_bit(data, n, cmax):
    """Drawn `C` (a random pattern, the diagonal in it, the rest zero) and a
    drawn running maximum: the same verdict, |C|'s largest entry and new
    running maximum, read whole and on the pattern."""
    mask = np.array(data.draw(st.lists(st.booleans(), min_size=n * n, max_size=n * n)))
    mask = mask.reshape(n, n) | np.eye(n, dtype=bool)
    vals = np.array(data.draw(st.lists(_entry, min_size=n * n, max_size=n * n))).reshape(n, n)
    C = np.where(mask, vals, 0.0)
    pat = np.flatnonzero(mask).astype(np.int64)
    py = _python(C, cmax)
    whole, on_pattern = _c(C, cmax), _c(C, cmax, pat)
    assert whole is not None and on_pattern is not None
    for got in (whole, on_pattern):
        assert got[0] == py[0]
        assert float(got[1]).hex() == float(py[1]).hex()
        assert float(got[2]).hex() == float(py[2]).hex()


def _packed_of(C, mask):
    """`C` on the pattern `mask` as a `_packed.PatternMatrix` (its values in
    the pattern's own order)."""
    from pycircuit.circuit import _packed
    rows, cols = np.nonzero(mask)
    pt = _packed.Pattern(C.shape[0], rows.astype(np.int64), cols.astype(np.int64))
    return _packed.PatternMatrix(pt, C[pt.rows, pt.cols].copy())


def _same(got, py):
    assert got is not None
    assert got[0] == py[0]
    assert float(got[1]).hex() == float(py[1]).hex()
    assert float(got[2]).hex() == float(py[2]).hex()


@settings(deadline=None)
@given(data=st.data(), n=st.integers(1, 8), cmax=_cmax)
def test_past_dense_max_the_screen_without_the_reference_is_the_python_one(data, n, cmax):
    """Past `_packed.DENSE_MAX` (here from one unknown) on the packed C, a
    reference node drawn (`irefnode`): its diagonal out of the scale and of
    the diagonal's tests on both sides (`_tran_branch._off_reference`,
    branch.h's `iref1`) -- the same verdict, |C|'s largest entry and running
    maximum (where the C serves: the packed sets hand a fired cheap level
    back); a dense C there, and the packed C below DENSE_MAX, leave nothing
    out, the same as without a reference."""
    from pycircuit.circuit import _packed
    mask = np.array(data.draw(st.lists(st.booleans(), min_size=n * n, max_size=n * n)))
    mask = mask.reshape(n, n) | np.eye(n, dtype=bool)
    vals = np.array(data.draw(st.lists(_entry, min_size=n * n, max_size=n * n))).reshape(n, n)
    g = data.draw(st.integers(0, n - 1))
    C = np.where(mask, vals, 0.0)
    pm = _packed_of(C, mask)
    pat = np.flatnonzero(mask).astype(np.int64)
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(_packed, 'DENSE_MAX', 1)
        py = _python(pm, cmax, g)
        got = _c(pm, cmax, g=g)
        if not (got is None and py[0] == _SVD):
            _same(got, py)
        ## (a dense C past DENSE_MAX: nothing left out)
        dense = _python(C, cmax, g)
        assert dense[0] == _python(C, cmax)[0] and float(dense[1]).hex() == float(_python(C, cmax)[1]).hex()
        _same(_c(C, cmax, g=g), dense)
        _same(_c(C, cmax, pat, g=g), dense)
    ## (below DENSE_MAX: nothing left out, packed or dense)
    below = _python(pm, cmax, g)
    plain = _python(pm, cmax)
    assert below[0] == plain[0] and float(below[1]).hex() == float(plain[1]).hex()
    got = _c(pm, cmax, g=g)
    if not (got is None and below[0] == _SVD):
        _same(got, below)


def test_past_dense_max_the_reference_diagonal_is_out_of_the_scale(monkeypatch):
    """The reference's diagonal the largest entry by far (a grid's: every
    capacitance to ground summed) and a capacitor under 1e-6 of it but not
    of the rest: past DENSE_MAX, packed, quiet at the cheap level with |C|'s
    largest entry the rest's, on both sides; below it, or without a
    reference, the cheap level fires on it (the packed sets handing the SVD
    back), as it did."""
    from pycircuit.circuit import _packed
    n, g = 6, 2
    C = np.diag([1e-12, 2e-12, 0.0, 1e-12, 1e-12, 3e-17])
    C[g, g] = 5e-11
    pm = _packed_of(C, np.ones((n, n), dtype=bool))
    assert _python(pm, None, g)[0] == _SVD and _c(pm, None, g=g) is None
    monkeypatch.setattr(_packed, 'DENSE_MAX', 1)
    assert _python(pm, None, g) == ((False, None), 2e-12, 2e-12)
    assert _c(pm, None, g=g) == ((False, None), 2e-12, 2e-12)
    assert _python(pm, None)[0] == _SVD and _c(pm, None) is None
    assert _python(C, None, g)[0] == _SVD and _c(C, None, g=g)[0] == _SVD
    ## (branch.c's own guards: an index past n, or a C not packed, leaves
    ## nothing out)
    assert _c(pm, None, g=g) is not None
    drv = _tran_branch_c.driver()
    s, cfn = drv[2], drv[1]
    assert s.sets == 1 and s.iref1 == g + 1
    for iref1, sets in ((n + 1, 1), (g + 1, 0)):
        s.iref1, s.sets, s.cmax = iref1, sets, 0.0
        if not sets:
            s.npat = 0
            buf = np.ascontiguousarray(C)
            s.C = drv[0].from_buffer('double *', buf)
        assert cfn(s) == _tran_branch_c.SVD and s.nrm == 5e-11, (iref1, sets, s.nrm)


def test_a_non_finite_entry_hands_the_screen_back():
    for bad in (np.nan, np.inf, -np.inf):
        C = np.eye(3)
        C[1, 2] = bad
        before = _paths.snapshot()
        assert _c(C, 1.0) is None
        assert _paths.since(before).get('branch_c:bail:nonfinite', 0) == 1


def test_the_declines():
    """A large `C` without the core's pattern (numpy's reductions), a
    patched `_branch_screen_C`, a `C` that is not a C-contiguous float64
    square array."""
    before = _paths.snapshot()
    assert _c(np.eye(65), 1.0) is None
    assert _c(np.eye(65), 1.0, np.arange(65, dtype=np.int64) * 66) is not None
    assert _c(np.asfortranarray(np.arange(4.0).reshape(2, 2)), 1.0) is None
    assert _c(np.eye(3, dtype=np.float32), 1.0) is None
    t = _Screen()
    t._branch_screen_C = lambda C, r0: None
    assert _tran_branch_c.screen(t, np.eye(3)) is None
    d = _paths.since(before)
    assert d.get('branch_c:large', 0) == 1 and d.get('branch_c:C', 0) == 2
    assert d.get('branch_c:patched', 0) == 1


def test_a_pattern_missing_the_largest_entry_is_another_screen():
    """Alive: the pattern route reads only the pattern -- leave out the
    entry holding |C|'s largest and the screen's largest entry moves."""
    C = np.eye(4) * 1e-3
    C[0, 3] = 5.0
    pat = np.flatnonzero(C).astype(np.int64)
    full = _c(C, 0.0, pat)
    short = _c(C, 0.0, pat[pat != 3])
    assert full[1] == 5.0 and short[1] != full[1]


## -- whole runs ---------------------------------------------------------------

def _mos_chain(n=6, cap=1e-13):
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for k in range(n):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'c{k}'] = C(f'd{k}', gnd, c=cap)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


def _run(build, on, monkeypatch, **kw):
    monkeypatch.setattr(_tran_branch_c, 'ENABLED', on)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
        res = tr.solve(**kw)
    st_ = tr.statistics
    stats = {k: getattr(st_, k) for k in st_.__slots__ if 'seconds' not in k} \
        if hasattr(st_, '__slots__') else {k: v for k, v in vars(st_).items() if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), stats, sorted(str(w.message) for w in W),
            getattr(tr, '_branch_cmax', None), _paths.since(before))


@pytest.mark.parametrize('name, build, kw', [
    ('chain', _mos_chain, {'tend': 3e-8, 'timestep': 2e-10}),
    ('chain-trap', _mos_chain, {'tend': 3e-8, 'timestep': 2e-10,
                                'make': {'integrator': TrapezoidalIntegrator()}}),
    ('chain-fixed', _mos_chain, {'tend': 2e-8, 'timestep': 2e-10, 'fixed_timestep': True,
                                 'make': {'integrator': Gear2Integrator()}}),
    ('sparse', lambda: _mos_chain(250), {'tend': 8e-9, 'timestep': 2e-10}),
])
def test_a_run_is_the_same_with_the_screen_in_c(name, build, kw, monkeypatch):
    off = _run(build, False, monkeypatch, **dict(kw))
    on = _run(build, True, monkeypatch, **dict(kw))
    assert on[:4] == off[:4]
    assert on[4].get('branch_c:served', 0) > 5
    if name == 'sparse':
        assert on[4].get('branch_c:pattern', 0) > 5
