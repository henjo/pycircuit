"""THE FROZEN FACTORS IN C (speed round 12, stage 5b, `_tran_radau_tc.fold`):
`_radau_frozen`'s two factors -- ``(lam0/h) Cr + Gr`` and ``(lam1/h) Cr +
Gr`` -- with the complex one's packed values, from one C call, every
decision left where it was.  The factors numpy's bytes and its errstate's
flags, enumerated and drawn; the packed values `_csc_of_dense`'s where the
pattern is the last record's, none elsewhere; a non-finite entry numpy's;
`prepare_values` `prepare`'s record; the frozen transform, a radau PSS and
the transform transients the same with the fold off."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, circuit, linearsolver
from pycircuit.circuit import _tran_radau as R
from pycircuit.circuit import _tran_radau_tc as TC
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.tests.test_radau_tc import _on
from pycircuit.circuit.transient import Transient


def _fold_or_skip():
    _on()
    if TC._fold_driver() is None:
        pytest.skip(f"the fold is off: {TC._FOLD.get('status')}")


def _bits(b):
    return float(np.frombuffer(np.uint64(b).tobytes(), dtype=np.float64)[0])


SPECIAL = (0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e-200, -1e-200, 1e200,
           -1e200, 1.7e308, np.inf, -np.inf, _bits(0x7FF8000000000000),
           _bits(0xFFF8000000000123), 1.0, -2.5)


def _transient():
    circuit.default_toolkit = circuit.numeric
    return Transient(_stage(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                     radau_transform=True)


_LAM = []


def _lam():
    """The transform's eigenvalues, as `_radau_frozen` reads them (of a
    transient that has stepped: the tableau is its integrator's)."""
    if not _LAM:
        tr = _transient()
        tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
        _LAM.append(tr._radau_transform_matrices()[0])
    return _LAM[0]


def _scalars(h):
    lam = _lam()
    return lam[0].real / h, lam[1] / h


def _numpy(Cr, Gr, s0, c1):
    """`_radau_frozen`'s factors, its way: None where its errstate raises."""
    try:
        with np.errstate(all='raise'):
            return s0 * Cr + Gr, c1 * Cr + Gr
    except FloatingPointError:
        return None


def _check(Cr, Gr, s0, c1, last=None):
    """`fold` against numpy -- the kind of its answer."""
    got = TC.fold(Cr, Gr, s0, c1, last)
    if not (np.isfinite(Cr).all() and np.isfinite(Gr).all()):
        assert got is None
        return 'numpy'
    ref = _numpy(Cr, Gr, s0, c1)
    if ref is None:
        assert got == (True, None, None, None)
        return 'flagged'
    flagged, rf, cf, data = got
    assert flagged is False
    for a, b in ((rf, ref[0]), (cf, ref[1])):
        assert a.dtype == b.dtype and a.shape == b.shape and a.flags.c_contiguous
        assert a.tobytes() == b.tobytes()
    if last is None:
        assert data is None
        return 'made'
    rec, pat = linearsolver._csc_of_dense(ref[1], linearsolver._csc_matvec(), last)
    if pat[1] is not last[1] or not last[2].shape[0]:
        ## (another pattern, or one without nonzeros: `prepare`'s to make)
        assert data is None
        return 'other'
    assert data is not None and data.dtype == np.complex128 and data.flags.c_contiguous
    assert data.tobytes() == rec[4].tobytes()
    return 'same'


def _pattern(Cr, Gr, s0, c1):
    with np.errstate(all='ignore'):
        cf = c1 * Cr + Gr
    return linearsolver._csc_of_dense(cf, linearsolver._csc_matvec(), None)[1]


@pytest.mark.parametrize('where', ['C', 'G'])
@pytest.mark.parametrize('v', SPECIAL, ids=[repr(v) for v in SPECIAL])
def test_a_special_value_anywhere_is_numpys(v, where):
    """Each value at every position of either input, at two steps: numpy's
    factors, its flags (overflow, underflow) as the flagged answer, a
    non-finite entry left to numpy; with the pattern of the plain inputs."""
    _fold_or_skip()
    rng = np.random.default_rng(5)
    m = 4
    for h in (2.5e-8, 1e-12):
        s0, c1 = _scalars(h)
        for k in range(m * m):
            Cr = 1e-12 * rng.standard_normal((m, m))
            Gr = 1e-3 * rng.standard_normal((m, m))
            last = _pattern(Cr, Gr, s0, c1)
            (Cr if where == 'C' else Gr).flat[k] = v
            _check(Cr, Gr, s0, c1, last)


def test_flags_and_patterns_are_reached():
    """The enumeration's kinds all occur: made, flagged (an overflow, an
    underflow), left to numpy, the same pattern and another."""
    _fold_or_skip()
    s0, c1 = _scalars(2.5e-8)
    rng = np.random.default_rng(9)
    Cr, Gr = 1e-12 * rng.standard_normal((3, 3)), 1e-3 * rng.standard_normal((3, 3))
    last = _pattern(Cr, Gr, s0, c1)
    assert _check(Cr, Gr, s0, c1) == 'made'
    assert _check(Cr * 1.5, Gr * 0.5, s0, c1, last) == 'same'
    big, tiny, nan = Cr.copy(), Cr.copy(), Gr.copy()
    ## (an overflow past ~1.8e308, a subnormal product: with s0 ~ 1.5e8)
    big[1, 1], tiny[1, 1], nan[2, 0] = 1e305, 3e-318, np.nan
    assert _check(big, Gr, s0, c1, last) == 'flagged'
    assert _check(tiny, Gr, s0, c1, last) == 'flagged'
    assert _check(Cr, nan, s0, c1, last) == 'numpy'
    Cz, Gz = Cr.copy(), Gr.copy()
    Cz[0, 2], Gz[0, 2] = 0.0, -0.0
    assert _check(Cz, Gz, s0, c1, last) == 'other'


def test_the_packed_values_fall_in_the_last_pattern_or_none_are_made():
    """Structural zeros (of both inputs) kept: the values, in its order; a
    zero appearing or filled: none; a real part cancelling exactly is no
    zero of the complex factor (its imaginary part is not)."""
    _fold_or_skip()
    s0, c1 = _scalars(2.5e-8)
    rng = np.random.default_rng(11)
    m = 5
    Cr, Gr = 1e-12 * rng.standard_normal((m, m)), 1e-3 * rng.standard_normal((m, m))
    for i, j in ((0, 1), (2, 3), (4, 0)):
        Cr[i, j] = Gr[i, j] = 0.0
    last = _pattern(Cr, Gr, s0, c1)
    assert _check(Cr * 1.1, Gr * 0.9, s0, c1, last) == 'same'
    C2, G2 = Cr.copy(), Gr.copy()
    C2[1, 1], G2[1, 1] = 0.0, 0.0
    assert _check(C2, G2, s0, c1, last) == 'other'
    C3 = Cr.copy()
    C3[0, 1] = 1e-12
    assert _check(C3, Gr, s0, c1, last) == 'other'
    G4 = Gr.copy()
    G4[3, 3] = -(s0 * Cr[3, 3])
    with np.errstate(all='ignore'):
        G4[3, 3] = -float((c1 * Cr[3, 3]).real)
    assert _check(Cr, G4, s0, c1, last) == 'same'
    ## (a pattern laid out as `A.T` as well as the usual `A`)
    lastc = (np.ascontiguousarray(last[0]), *last[1:])
    assert _check(Cr * 1.1, Gr * 0.9, s0, c1, lastc) == 'same'
    assert _check(C2, G2, s0, c1, lastc) == 'other'


@settings(deadline=None, max_examples=300)
@given(data=st.data())
def test_drawn_factors_are_numpys(data):
    _fold_or_skip()
    m = data.draw(st.integers(1, 9))
    rng = np.random.default_rng(data.draw(st.integers(0, 2 ** 32 - 1)))
    sc, sg = data.draw(st.sampled_from(((1e-12, 1e-3), (1.0, 1.0), (1e-160, 1e-160),
                                        (1e150, 1e150), (1e-300, 1.0))))
    Cr, Gr = sc * rng.standard_normal((m, m)), sg * rng.standard_normal((m, m))
    for _ in range(data.draw(st.integers(0, 3))):
        Z = Cr if data.draw(st.booleans()) else Gr
        Z.flat[data.draw(st.integers(0, m * m - 1))] = 0.0
    h = data.draw(st.sampled_from((2.5e-8, 1e-9, 1e-12, 1.0, 1e-280, 1e280)))
    s0, c1 = _scalars(h)
    last = _pattern(Cr, Gr, s0, c1) if data.draw(st.booleans()) else None
    for _ in range(data.draw(st.integers(0, 3))):
        Z = Cr if data.draw(st.booleans()) else Gr
        Z.flat[data.draw(st.integers(0, m * m - 1))] = data.draw(st.sampled_from(SPECIAL))
    _check(Cr, Gr, s0, c1, last)


def test_inputs_of_another_kind_are_numpys(monkeypatch):
    """A float32 input, a non-square one, two sizes, the switch off: no
    factors (counted), the caller's numpy -- a strided input is copied in."""
    _fold_or_skip()
    s0, c1 = _scalars(2.5e-8)
    rng = np.random.default_rng(3)
    Cr, Gr = 1e-12 * rng.standard_normal((4, 4)), 1e-3 * rng.standard_normal((4, 4))
    before = _paths.snapshot()
    assert TC.fold(Cr.astype(np.float32), Gr, s0, c1, None) is None
    assert TC.fold(Cr[:, :3], Gr[:, :3], s0, c1, None) is None
    assert TC.fold(Cr, Gr[:3, :3], s0, c1, None) is None
    monkeypatch.setattr(TC, 'FOLD', False)
    assert TC.fold(Cr, Gr, s0, c1, None) is None
    monkeypatch.undo()
    assert _paths.since(before) == {'radau.frozen:fold.input': 3, 'radau.frozen:fold.off': 1}
    big = 1e-12 * rng.standard_normal((8, 8))
    assert _check(big[::2, ::2], Gr, s0, c1) == 'made'


def test_prepare_values_is_prepares_record(monkeypatch):
    """`prepare_values(data, A)` where `A`'s nonzeros fall in the last
    pattern: `prepare(A)`'s record and pattern; no pattern, or no numpy
    marshalling: `prepare(A)` itself."""
    try:
        z1, z2 = linearsolver.ComplexKLUSolver(), linearsolver.ComplexKLUSolver()
    except ImportError as e:
        pytest.skip(f'libklu not available: {e}')
    rng = np.random.default_rng(13)
    A = rng.standard_normal((5, 5)) + 1j * rng.standard_normal((5, 5))
    A[0, 3] = A[4, 1] = 0.0
    A2 = 1.5 * A
    assert z2.prepare_values(A2.T[A2.T != 0], A2)[5] == z1.prepare(A2)[5]
    z1.prepare(A)
    z2.prepare(A)
    r1, r2 = z1.prepare(A2), z2.prepare_values(A2.T[A2.T != 0], A2)
    assert r1[1] == r2[1] and r1[5] == r2[5]
    for k in (2, 3, 4):
        assert r1[k].dtype == r2[k].dtype and r1[k].tobytes() == r2[k].tobytes()
    assert r2[2] is z2._csc_last[1] and r2[3] is z2._csc_last[2]
    x = rng.standard_normal(5) + 1j * rng.standard_normal(5)
    assert r1[0](x).tobytes() == r2[0](x).tobytes()
    assert z1._csc_last[0].tobytes() == z2._csc_last[0].tobytes()
    monkeypatch.setattr(linearsolver, '_csc_matvec', lambda: None)
    r3 = z2.prepare_values(A2.T[A2.T != 0], A2)
    assert r3[2].tobytes() == r1[2].tobytes() and r3[4].tobytes() == r1[4].tobytes()


def _frozen_snaps(on, monkeypatch):
    """Every frozen transform of a radau PSS of the PSP stage, as
    `_radau_frozen` returns it (its factors, kept LU, record and the
    solver's pattern), and the run's counts."""
    snaps = []
    orig = R._RadauStages._radau_frozen

    def spy(self, Cr, Gr, h):
        fz = orig(self, Cr, Gr, h)
        if fz is None:
            snaps.append(None)
            return fz
        lu, pr = fz.lu, fz.prep
        last = self._radau_zsolver._csc_last
        snaps.append((fz.rf.tobytes(), fz.cf.tobytes(), fz.m, fz.p, fz.v,
                      None if lu is None else (lu._a[0].tobytes(), lu._info, lu.n),
                      None if pr is None else (pr[0].func, pr[0].args[1], pr[0].args[4].tobytes(),
                                               pr[1], pr[2].tobytes(), pr[3].tobytes(),
                                               pr[4].tobytes(), pr[5]),
                      None if last is None else (last[0].tobytes(), last[1].tobytes(),
                                                 last[2].tobytes(), last[3])))
        return fz
    monkeypatch.setattr(R._RadauStages, '_radau_frozen', spy)
    out, counts = _pss(on)
    monkeypatch.undo()
    return snaps, out, counts


def _pss(on):
    from pycircuit.circuit.shooting import PSS
    old = TC.FOLD
    TC.FOLD = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            p = PSS(_stage(), method='radau', reltol=1e-8)
            before = _paths.snapshot()
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
            d = _paths.since(before)
    finally:
        TC.FOLD = old
    out = (np.asarray(p.waveform[1], float).tobytes(), np.asarray(p._monodromy, float).tobytes(),
           sorted(str(w.message) for w in W))
    return out, {k: v for k, v in d.items() if not k.startswith('once:')}


def _transient_run(on, **kw):
    old = TC.FOLD
    TC.FOLD = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = _transient()
            before = _paths.snapshot()
            res = tr.solve(tend=6e-7, timestep=2e-8, **kw)
            d = _paths.since(before)
    finally:
        TC.FOLD = old
    out = (np.asarray(res.x, float).tobytes(), sorted(str(w.message) for w in W))
    return out, {k: v for k, v in d.items() if not k.startswith('once:')}


def _moved_but_own(da, db):
    keys = {k for k in set(da) | set(db) if not k.startswith('radau.frozen:fold')}
    return {k: (da.get(k, 0), db.get(k, 0)) for k in sorted(keys) if da.get(k, 0) != db.get(k, 0)}


def test_the_frozen_transform_is_the_same_with_the_fold_off(monkeypatch):
    _fold_or_skip()
    sa, a, da = _frozen_snaps(False, monkeypatch)
    sb, b, db = _frozen_snaps(True, monkeypatch)
    assert len(sa) == len(sb) > 100
    for k, (x, y) in enumerate(zip(sa, sb, strict=True)):
        assert x == y, k
    assert a == b and _moved_but_own(da, db) == {}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
def test_a_transform_transient_is_the_same_with_the_fold_off(fixed):
    _fold_or_skip()
    kw = {'fixed_timestep': True} if fixed else {}
    a, da = _transient_run(False, **kw)
    b, db = _transient_run(True, **kw)
    assert a == b
    assert _moved_but_own(da, db) == {}


def test_the_factors_are_folded():
    """Served: on the PSP stage's radau PSS every step's factors from the C
    but where the step's `P` is new, each in the last record's pattern (the
    parent: none)."""
    _fold_or_skip()
    _b, db = _pss(True)
    assert db.get('radau.frozen:fold', 0) > 100, db
    assert db.get('radau.frozen:fold.same', 0) == db['radau.frozen:fold'], db
