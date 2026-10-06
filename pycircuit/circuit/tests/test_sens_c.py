"""The shooting's coupled stage step map in C (speed round 11, item 4,
`shooting/_sens_c.py`): the plain dense walk's monodromy step -- the stage
block from the stages' full `C` and `G`, SciPy's `getrf`, `C_n P` stacked,
SciPy's `getrs` -- against the Python step on the same inputs (the
readers' reduction, `_stage_block` or the loop, `lu_factor`, `solve`'s
`lu_solve`): drawn and enumerated blocks -- sizes across OpenBLAS'
thresholds, every reference row, signed zeros, subnormals and overflow
(the assembly's flags), non-finite values, singular blocks -- the C
serving the same new `P` (bytes and strides) and handing back exactly
where the Python step declines, warns or raises; whole radau PSSes the
same with it off; every decline counted."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, circuit
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.linearsolver import DenseSolver
from pycircuit.circuit.shooting import _numerics, _pss_walks, _sens_c
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage

#: Radau IIA(3)'s Butcher matrix
_A = np.array([[0.19681547722366, -0.06553542585020, 0.02377097434822],
               [0.39442431473909, 0.29207341166523, -0.04154875212600],
               [0.37640306270047, 0.51248582618842, 0.11111111111111]])


def _drv():
    drv = _sens_c.driver()
    if drv is None:
        pytest.skip(f'the C is off: {_sens_c.STATUS}')
    return drv


def _python_map(Cfs, Gfs, iref, h, A, base):
    """The Python step on the same inputs, as `_stage_step` and
    `_StageStep.solve` make it: ('ok', P) | ('raise', type) | ('warn', P),
    and `_stage_block`'s counts."""
    tk = circuit.numeric
    Cs = [np.asarray(remove_row_col((C,), iref, tk)[0]) for C in Cfs]
    Gs = [np.asarray(remove_row_col((G,), iref, tk)[0]) for G in Gfs]
    s, m = A.shape[0], Cs[0].shape[0]
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        try:
            Jb = _pss_walks._stage_block(Cs, Gs, h, A)
            if Jb is None:
                Jb = np.zeros((s * m, s * m))
                for i in range(s):
                    for jj in range(s):
                        blk = h * A[i, jj] * Gs[jj]
                        if i == jj:
                            blk = Cs[i] + blk
                        Jb[i * m:(i + 1) * m, jj * m:(jj + 1) * m] = blk
            fac = DenseSolver().factor(Jb, tk)
            Z = _numerics._lu_solve_split(fac, np.vstack([base] * s))
            out = ('ok', Z[(s - 1) * m:s * m])
        except ValueError:
            out = ('raise', ValueError)
    d = {k: v for k, v in _paths.since(before).items() if k.startswith('pss.block:')}
    if out[0] == 'ok' and any('Singular' in str(w.message) for w in W):
        out = ('warn', out[1])
    return out, d


def _check(Cfs, Gfs, iref, h, A, base):
    """The C against the Python step: the outcome each makes."""
    c = _sens_c._map({}, _drv(), A, Cfs, Gfs, iref, h, base)
    p, blk = _python_map(Cfs, Gfs, iref, h, A, base)
    if type(c) is not int:
        assert p[0] == 'ok' and blk == {'pss.block:served': 1}, (p[0], blk)
        P = p[1]
        assert c.dtype == P.dtype and c.shape == P.shape and c.strides == P.strides
        assert c.tobytes() == P.tobytes()
    elif c == 2:                    # a non-finite input: declined, then raised
        assert blk == {'pss.block:nonfinite': 1} and p == ('raise', ValueError), (blk, p[0])
    elif c == 3:                    # a flag in the assembly: `_stage_block`'s errstate
        assert blk == {'pss.block:fp': 1}, blk
    elif c == 4:                    # a non-finite `C_n P`: `lu_solve` refuses it --
        ## after the block assembled, or after it declined on a flag in its
        ## assembly: the C tests `C_n P` before it assembles, the Python after
        ## (a drawn case, found 2026-10-06 when other draws hit it)
        assert blk in ({'pss.block:served': 1}, {'pss.block:fp': 1}) \
            and p == ('raise', ValueError), (blk, p[0])
    elif c == 5:                    # a singular block: `lu_factor` warns
        assert blk == {'pss.block:served': 1} and p[0] == 'warn', (blk, p[0])
    else:
        raise AssertionError(f'status {c}')
    return c if type(c) is int else 1


def _bits(b):
    return float(np.frombuffer(np.uint64(b).tobytes(), dtype=np.float64)[0])


SPECIAL = (0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e300, -1e300, np.inf,
           -np.inf, _bits(0x7FF8000000000000), _bits(0xFFF8000000000123), 1.0, -2.5)


def _system(rng, s, nf, scale=1.0):
    """Stages' full `C` and `G` of a well-posed block: `C` near a small
    identity, `G` a diagonal-heavy conductance."""
    Cfs = [1e-12 * (np.eye(nf) + 0.1 * rng.standard_normal((nf, nf))) for _ in range(s)]
    Gfs = [scale * (np.eye(nf) * 1e-2 + 1e-3 * rng.standard_normal((nf, nf))) for _ in range(s)]
    return Cfs, Gfs


@settings(deadline=None, max_examples=200)
@given(data=st.data())
def test_a_drawn_step_map_is_the_python_one(data):
    s = data.draw(st.sampled_from((3, 1, 2, 4)))
    A = _A if s == 3 else np.array(data.draw(st.lists(st.floats(-1, 1), min_size=s * s,
                                                      max_size=s * s))).reshape(s, s)
    nf = data.draw(st.integers(2, 12))
    iref = data.draw(st.integers(0, nf - 1))
    k = data.draw(st.sampled_from((nf - 1, 1, 3)))
    rng = np.random.default_rng(data.draw(st.integers(0, 2 ** 32 - 1)))
    Cfs, Gfs = _system(rng, s, nf, data.draw(st.sampled_from((1.0, 1e-200, 1e200))))
    base = rng.standard_normal((nf - 1, k))
    ## special values planted anywhere -- the reference row too
    for _ in range(data.draw(st.integers(0, 3))):
        which = data.draw(st.sampled_from(('C', 'G', 'base')))
        v = data.draw(st.sampled_from(SPECIAL))
        M = (Cfs if which == 'C' else Gfs)[data.draw(st.integers(0, s - 1))] if which != 'base' else base
        M.flat[data.draw(st.integers(0, M.size - 1))] = v
    h = data.draw(st.sampled_from((2.5e-8, 1e-9, 1.0, 1e300, 5e-324, 0.0, -0.0, np.inf)))
    with np.errstate(all='ignore'):
        _check(Cfs, Gfs, iref, h, A, base)


@pytest.mark.parametrize('m', [1, 6, 20, 33, 40, 50])
def test_sizes_across_the_library_thresholds(m):
    """Blocks of 3 to 150 unknowns (OpenBLAS' `getrf` blocks and threads
    past ~100): every reference row of the smaller, three of the larger."""
    rng = np.random.default_rng(m)
    nf = m + 1
    rows = range(nf) if nf <= 7 else (0, nf // 2, nf - 1)
    for iref in rows:
        Cfs, Gfs = _system(rng, 3, nf)
        assert _check(Cfs, Gfs, iref, 2.5e-8, _A, rng.standard_normal((m, m))) == 1


def test_every_hand_back():
    rng = np.random.default_rng(7)
    nf = 5
    Cfs, Gfs = _system(rng, 3, nf)
    base = rng.standard_normal((nf - 1, nf - 1))
    ## a non-finite input in the reduced part -- and in the reference row,
    ## which the readers' reduction drops: served
    bad = [c.copy() for c in Cfs]
    bad[1][2, 3] = np.nan
    assert _check(bad, Gfs, 0, 2.5e-8, _A, base) == 2
    bad = [c.copy() for c in Cfs]
    bad[1][0, 3] = np.nan
    assert _check(bad, Gfs, 0, 2.5e-8, _A, base) == 1
    ## the assembly's flags: an overflowing product (finite factors) and
    ## an underflowing one
    _C, G_big = _system(rng, 3, nf, 1e100)
    assert _check(Cfs, G_big, 0, 1e300, _A, base) == 3
    assert _check(Cfs, Gfs, 0, 5e-324, _A, base) == 3
    ## a non-finite `C_n P`
    b2 = base.copy()
    b2[1, 1] = np.inf
    assert _check(Cfs, Gfs, 0, 2.5e-8, _A, b2) == 4
    ## ... with a flag in the assembly as well: the C sees `C_n P` first
    assert _check(Cfs, G_big, 0, 1e300, _A, b2) == 4
    ## a singular block
    Z = [np.zeros((nf, nf)) for _ in range(3)]
    assert _check(Z, Z, 0, 2.5e-8, _A, base) == 5


## -- whole PSSes ----------------------------------------------------------------------

def _pss(build, on, monkeypatch, **kw):
    from pycircuit.circuit.shooting import PSS
    monkeypatch.setattr(_sens_c, 'SENS_C', on)
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        p = PSS(build(), method='radau', reltol=1e-8, **kw)
        before = _paths.snapshot()
        p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
        d = _paths.since(before)
    monkeypatch.undo()
    wf = p.waveform
    return ((np.asarray(wf[0], float).tobytes(), np.asarray(wf[1], float).tobytes(),
             np.asarray(p._monodromy, float).tobytes(), sorted(str(w.message) for w in W)), d)


@pytest.mark.parametrize('build', [_stage, lambda: _mos_chain(cap=1e-13)], ids=['psp', 'mos'])
def test_a_radau_pss_is_the_same_with_the_c_off(build, monkeypatch):
    a, da = _pss(build, False, monkeypatch)
    b, db = _pss(build, True, monkeypatch)
    assert a == b
    ## the counts: the readers' sync skips the same; the C's step maps in
    ## place of the stage reads' and the block's
    assert da.get('pss.sync:skipped') == db.get('pss.sync:skipped')
    served = db.get('pss.sens:served', 0)
    assert da.get('pss.reads:served', 0) - db.get('pss.reads:served', 0) == served
    assert da.get('pss.block:served', 0) - db.get('pss.block:served', 0) == served


def test_the_c_serves_a_radau_pss(monkeypatch):
    """The PSP stage's radau PSS: every dense step but the walks' first maps
    in C.  (Apart from the comparison above: with the other fast paths off
    the step's memo holds no stage matrices, and this fails.)"""
    _drv()
    _b, d = _pss(_stage, True, monkeypatch)
    assert d.get('pss.sens:served', 0) >= 100, d


def _declines(monkeypatch, build=_stage, **kw):
    _b, d = _pss(build, True, monkeypatch, **kw)
    return {k: v for k, v in d.items() if k.startswith('pss.sens:') and k != 'pss.sens:served'}


def test_the_declines(monkeypatch):
    """Each reason, on a whole PSS: the switch, a piece patched (on a
    class, in a module), a solver that is not the dense one, the trims
    off (the readers' sync would run), a system past `MAXN` unknowns, a
    tableau past `MAXS`."""
    _drv()
    with monkeypatch.context() as mp:
        mp.setattr(_sens_c, 'SENS_C', False)
        _b, d = _pss(_stage, False, mp)
    assert d.get('pss.sens:off', 0) >= 100 and 'pss.sens:served' not in d
    import scipy.linalg
    real = scipy.linalg.lu_factor
    with monkeypatch.context() as mp:
        mp.setattr(scipy.linalg, 'lu_factor', lambda *a, **k: real(*a, **k))
        assert set(_declines(mp)) == {'pss.sens:patched'}
    from pycircuit.circuit.shooting import _steps
    real_solve = _steps._StageStep.solve
    with monkeypatch.context() as mp:
        mp.setattr(_steps._StageStep, 'solve', lambda self, *a, **k: real_solve(self, *a, **k))
        assert set(_declines(mp)) == {'pss.sens:patched'}
    from pycircuit.circuit.linearsolver import SuperLUSolver
    with monkeypatch.context() as mp:
        assert 'pss.sens:solver' in _declines(mp, linearsolver=SuperLUSolver())
    from pycircuit.circuit.shooting import _pss_inner
    with monkeypatch.context() as mp:
        mp.setattr(_pss_inner, 'SHOOT_TRIM', False)
        assert set(_declines(mp)) == {'pss.sens:limits'}
    with monkeypatch.context() as mp:
        mp.setattr(_sens_c, 'MAXN', 17)
        d = _declines(mp)
    assert d.get('pss.sens:size', 0) >= 100 and d.get('pss.sens:solver') == 1, d
    with monkeypatch.context() as mp:
        mp.setattr(_sens_c, 'MAXS', 2)
        d = _declines(mp)
    ## (and, before it, the first step's: the solver not chosen yet --
    ## `AutoSolver._select` chooses and counts that, in the Python step)
    assert d.get('pss.sens:tableau', 0) >= 100 and d.get('pss.sens:solver') == 1, d
