"""Radau's transform Newton in C (`_tran_radau_tc`, speed round 10, B3.3):
the cost transform's simplified Newton loop (`_transform_loop`) in one C
call.  Its stages, the source memo, the kept LU's and the complex solver's
state and the warnings are the Python loop's -- on whole transients and a
PSS, and on one captured step driven from drawn seeds, step sizes and
entering charges; every hand-back (a non-finite stage, a floating-point
exception, a singular real factor, a failed refactor, the residual check,
the walk stopping, maxiter) leaves the Python loop to answer as before,
with no trace; numpy's complex product read, not assumed; the inputs
it keeps from step to step set again wherever they change (B3.5); and
the declines."""
import ctypes
import math
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _hdl_climit, _paths, circuit
from pycircuit.circuit import _tran_radau_tc as TC
from pycircuit.circuit import linearsolver as LS
from pycircuit.circuit._tran_radau import _RadauStages
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_newton_c import _mos_chain
from pycircuit.circuit.tests.test_psp_limit_c import _stage
from pycircuit.circuit.transient import Transient

DBL_MAX = float(np.finfo(float).max)


def _on():
    try:
        LS.ComplexKLUSolver()
    except ImportError as e:
        pytest.skip(f'libklu not available: {e}')
    if TC.driver() is None:
        pytest.skip(f'the transform C is off: {TC.STATUS}')


def _run(build, on, make=None, **kw):
    old = TC.ENABLED
    TC.ENABLED = on
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(),
                           **(make or {}))
            before = _paths.snapshot()
            res = tr.solve(**kw)
            d = _paths.since(before)
    finally:
        TC.ENABLED = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W)), d


CASES = {
    'psp': (_stage, {}),
    'psp-transform': (_stage, {'radau_transform': True}),
    'mos-transform': (_mos_chain, {'radau_transform': True}),
    'mos-caps-transform': (lambda: _mos_chain(cap=1e-13), {'radau_transform': True}),
}


@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
@pytest.mark.parametrize('case', list(CASES))
def test_a_radau_transient_is_the_same_with_the_c_off(case, fixed):
    _on()
    build, make = CASES[case]
    kw = {'tend': 6e-7, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, _da = _run(build, False, make, **dict(kw))
    b, db = _run(build, True, make, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert db.get('radau_tc:served', 0) > 0, db


def test_a_radau_pss_is_the_same_with_the_c_off():
    """PSS's default method on a PSP circuit: the transform under 'auto'."""
    _on()
    from pycircuit.circuit.shooting import PSS

    def run(on):
        old = TC.ENABLED
        TC.ENABLED = on
        try:
            with warnings.catch_warnings(record=True) as W:
                warnings.simplefilter('always')
                warnings.simplefilter('ignore', ResourceWarning)
                p = PSS(_stage(), method='radau', reltol=1e-8)
                before = _paths.snapshot()
                p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
                d = _paths.since(before)
        finally:
            TC.ENABLED = old
        wf = p.waveform
        return (np.asarray(wf[0], float).tobytes(), np.asarray(wf[1], float).tobytes(),
                sorted(str(w.message) for w in W)), d
    a, _da = run(False)
    b, db = run(True)
    assert a == b
    assert db.get('radau_tc:served', 0) > 0, db


## -- one step's loop, captured: drawn seeds, steps and charges -------------------------

def _captured(build=_stage, steps=4, **make):
    """A radau transient on the transform path whose last step's loop was
    kept: the C's arguments, the step's frozen inputs `(Cr, Gr, h)` and a
    copy of the source memo as the loop found it."""
    tr = Transient(build(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator(), **make)
    got = {}
    real_fz = tr._radau_frozen
    real_tc = TC.solve

    def spy_fz(Cr, Gr, h):
        got.update(Cr=np.array(Cr), Gr=np.array(Gr), h=h)
        return real_fz(Cr, Gr, h)

    def spy_tc(tr_, ctx, fz, seed, src, pf, lims, nobypass, reltol, abstol, maxit):
        got.update(ctx=ctx, seed=[np.array(y, dtype=float) for y in seed], src=src,
                   u=dict(tr_.__dict__['_u_memo']), tols=(reltol, abstol, maxit))
        return real_tc(tr_, ctx, fz, seed, src, pf, lims, nobypass, reltol, abstol, maxit)
    tr._radau_frozen = spy_fz
    TC.solve = spy_tc
    try:
        tr.solve(tend=steps * 2e-8, timestep=2e-8, fixed_timestep=True)
    finally:
        TC.solve = real_tc
        del tr._radau_frozen
    assert 'ctx' in got, 'the transform path did not run'
    return tr, got


#: (only this module's captured-step tests read it: a captured transient is
#: recorded against the test that made it)
_CAP = {}


def _cap():
    if 'psp' not in _CAP:
        _CAP['psp'] = _captured()
    return _CAP['psp']


def _seed(dx=None):
    Y = [np.array(y, dtype=float) for y in _cap()[1]['seed']]
    if dx is not None:
        for j, d in dx.items():
            Y[j] = Y[j] + d
    return Y


#: a KLU entry that fails (returns 0), for the C's pointers alone
_FAIL6 = ctypes.CFUNCTYPE(ctypes.c_int, *([ctypes.c_void_p] * 6))(lambda *a: 0)


def _attempt(seed, on, h=None, qn=None, Cr=None, Gr=None, maxit=None, tol=None,
             zero_record=False, fail=None):
    """The captured step's loop from `seed`, as `_rk_step_transformed` runs
    it -- the C, and the Python loop where the C declines or hands back --
    the C on or off, from one state: the source memo as the loop found it,
    a fresh kept LU and a fresh complex solver primed as the step before
    leaves it (the record analysed and factored; its next call refactors).
    The stages' bytes (or the exception, or 'noconv'), the source memo and
    its counts, both factorisations' state, the warnings with their lines;
    and every count."""
    tr, got = _cap()
    ctx = got['ctx']
    reltol, abstol, maxit0 = got['tols']
    h0, qn0 = ctx.h, ctx.qn
    h = h0 if h is None else h
    tr.__dict__['_u_memo'] = dict(got['u'])
    zs = tr._radau_zsolver = LS.ComplexKLUSolver()
    old = TC.ENABLED
    before = _paths.snapshot()
    try:
        ctx.h = h
        if qn is not None:
            ctx.qn = qn
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            warnings.simplefilter('ignore', ResourceWarning)
            fz = tr._radau_frozen(got['Cr'] if Cr is None else Cr,
                                  got['Gr'] if Gr is None else Gr, h)
            if fz is not None and fz.prep is not None:
                zs.solve_prepared(fz.prep, np.ones(fz.m, dtype=complex))
                if zero_record:
                    fz.prep[4][:] = 0.0
            if tol is not None:
                zs.REFACTOR_RESIDUAL_TOL = tol
            if fail is not None:
                ## (the C's refactor or solve planted as failing; the Python's
                ## own calls are the library's)
                ref = ctypes.cast(zs._lib.klu_z_refactor, ctypes.c_void_p).value
                sol = ctypes.cast(zs._lib.klu_z_solve, ctypes.c_void_p).value
                bad = ctypes.cast(_FAIL6, ctypes.c_void_p).value
                TC._MOD['klu'] = (zs._lib, bad if fail == 'refactor' else ref,
                                  bad if fail == 'solve' else sol)
            Y = [np.array(y, dtype=float) for y in seed]
            TC.ENABLED = on
            try:
                Yc = TC.solve(tr, ctx, fz, Y, got['src'], None, [], True, reltol, abstol,
                              maxit0 if maxit is None else maxit)
                if Yc is not None:
                    out = [y.tobytes() for y in Yc]
                elif tr._transform_loop(ctx, fz, Y, [None] * 3, [], None, None, reltol,
                                        abstol, maxit0 if maxit is None else maxit):
                    out = [y.tobytes() for y in Y]
                else:
                    out = 'noconv'
            except Exception as e:                             # noqa: BLE001
                out = (type(e).__name__, str(e))
            finally:
                TC.ENABLED = old
                if fail is not None:
                    TC._MOD.pop('klu', None)
    finally:
        ctx.h, ctx.qn = h0, qn0
    d = _paths.since(before)
    um = sorted((repr(k), np.asarray(v).tobytes())
                for k, v in tr.__dict__['_u_memo'].items())
    src = {k: v for k, v in d.items() if k.startswith(('umemo:', 'u:', 'src.u:'))}
    fac = None
    if fz is not None:
        fac = (None if fz.lu is None else (fz.lu._info, fz.lu._a[0].tobytes()),
               zs.analyses, zs.factors, zs.refactors, zs.residual_fallbacks,
               zs._fresh is fz.prep)
    return (out, um, src, fac,
            [(w.category.__name__, str(w.message), w.filename, w.lineno) for w in W]), d


def _both(seed, **kw):
    """The captured loop from `seed` with the C off and on: the two
    outcomes and the C run's counts.  A hand-back leaves no trace: every
    count the Python's but the C's own and the once-only ones."""
    _on()
    a, da = _attempt(seed, False, **kw)
    b, db = _attempt(seed, True, **kw)
    if not db.get('radau_tc:served'):
        own = lambda d: {k: v for k, v in d.items()
                         if not k.startswith(('radau_tc:', 'once:'))}
        assert own(db) == own(da), (da, db)
    return a, b, db


def _counts(db):
    return {k: v for k, v in db.items() if k.startswith('radau_tc:')}


def _bail(db, why):
    c = _counts(db)
    assert c.get('radau_tc:bail:' + why) == 1 and 'radau_tc:served' not in c, c


def test_the_captured_step_is_served_and_the_pythons():
    a, b, db = _both(_seed())
    assert a == b
    assert _counts(db).get('radau_tc:served') == 1, db
    assert a[3][-1] and a[3][3] >= 1, a[3]          # (the record refactored, and kept)


def test_every_hand_back_answers_as_the_python():
    """Each way the C hands the loop back, planted: the Python loop's
    answer, exception, warnings, memo and factorisations, and the
    hand-back counted."""
    tr, got = _cap()
    n = len(_seed()[0])
    m = n - 1
    ## a non-finite stage: `passes` declines, the Python evaluates
    bad = _seed()
    bad[1][2] = np.nan
    a, b, db = _both(bad)
    assert a == b
    _bail(db, 'nonfinite')
    ## an overflow in the residual: numpy warns
    qn = np.where(np.arange(n) % 2, DBL_MAX, -DBL_MAX)
    a, b, db = _both(_seed(), qn=qn)
    assert a == b and any('overflow' in w[1] or 'invalid' in w[1] for w in a[4]), a[4]
    _bail(db, 'flags')
    ## a singular real factor (its first diagonal entry cancels exactly:
    ## `(gamma/h) 1 + (-(gamma/h))`), the complex one not: numpy's error,
    ## from the Python's first solve
    gam = tr._radau_transform_matrices()[0][0].real / got['ctx'].h
    Gs = np.eye(m)
    Gs[0, 0] = -gam
    a, b, db = _both(_seed(), Cr=np.eye(m), Gr=Gs)
    assert a == b and a[0] == ('LinAlgError', 'Singular matrix'), a[0]
    _bail(db, 'lu')
    ## the record's values zeroed after its factor: KLU's refactor of this
    ## block structure succeeds, its solve is not finite, the residual check
    ## hands back -- and the Python's fallback factor raises
    a, b, db = _both(_seed(), zero_record=True)
    assert a == b and a[0] == ('LinAlgError', 'Singular matrix'), a[0]
    _bail(db, 'residual')
    ## KLU's refactor or solve failing in the C (planted in its pointers
    ## alone): handed back with nothing left behind, the Python's answer
    for fail in ('refactor', 'solve'):
        a, b, db = _both(_seed(), fail=fail)
        assert a == b
        _bail(db, 'refactor' if fail == 'refactor' else 'klusolve')
    ## the residual check failing: the Python's fallback factors afresh
    a, b, db = _both(_seed(), tol=-1.0)
    assert a == b and a[3][4] >= 1, a[3]
    _bail(db, 'residual')
    ## maxiter: the Python loop does not converge either
    far = _seed({j: 0.3 for j in range(3)})
    a, b, db = _both(far, maxit=1)
    assert a == b and a[0] == 'noconv', a[0]
    _bail(db, 'maxiter')


def test_a_stopped_walk_hands_back_and_is_taken_again():
    """The walk's first entry planted as a stop (its kernel's address 0, so
    every setup takes it): the C hands back (`walkstop`), the Python's walk
    runs the loop's own statement for that element (the same law), and once
    the address is back the next C call takes the tables again and is
    served."""
    _on()
    tr = _cap()[0]
    a, _da = _attempt(_seed(), False)
    w = _hdl_climit._walk_for(tr.cir)
    addr0 = w.addr[0]
    w.addr[0] = 0
    w.F[0] = 0
    try:
        b, db = _attempt(_seed(), True)
    finally:
        w.addr[0] = addr0
        w.F[0] = addr0
    assert a == b
    _bail(db, 'walkstop')
    rec = tr.__dict__['_radau_tc']
    assert rec.walk_ok is None
    c, dc = _attempt(_seed(), True)
    assert c == a and _counts(dc).get('radau_tc:served') == 1, dc


def test_what_the_c_keeps_follows_each_step():
    """The C's kept inputs -- `P` and `V`, the kept LU's addresses, the
    pattern's arrays, KLU's handles -- set again wherever the step's
    change (speed round 10, B3.5): on one context, loops whose step size
    leaves and returns, whose LU is a fresh one while the last one's
    reader lives (a stale address would solve another factor), and whose
    complex factor loses an entry (a stale pattern would mismatch its
    values) -- each the Python loop's, served, and what is kept the
    step's own."""
    _on()
    tr, got = _cap()
    h0 = got['ctx'].h

    def kept():
        rec, zs = tr.__dict__['_radau_tc'], tr._radau_zsolver
        assert rec.p_of is tr._radau_Pc[1] and rec.lu_addr == tr._radau_lu.addr
        assert rec.Ap_obj is zs._csc_last[1] and rec.klu_kk[2:4] == (zs._symbolic, zs._numeric)
        return rec
    marks = []
    for h in (h0, h0 * 0.5, h0 * 0.5, h0):
        a, b, db = _both(_seed(), h=h)
        assert a == b and _counts(db).get('radau_tc:served') == 1, (h, db)
        marks.append(kept().p_of)
    assert marks[1] is not marks[0] and marks[2] is marks[1] and marks[3] is not marks[2]
    ## the last LU's reader alive: the C's LU is another one, at its own address
    addr = kept().lu_addr
    keep = tr._radau_frozen(got['Cr'] * 2.0, got['Gr'], h0)
    assert keep.lu.addr == addr
    a, b, db = _both(_seed())
    assert a == b and _counts(db).get('radau_tc:served') == 1, db
    assert kept().lu_addr != addr and keep.lu._a[0].tobytes() == keep.rf.tobytes()
    ## the complex factor's smallest off-diagonal entry removed: a new pattern
    Cr, Gr = np.array(got['Cr']), np.array(got['Gr'])
    m = Cr.shape[0]
    off = [(abs(Cr[i, j]) + abs(Gr[i, j]), i, j) for i in range(m) for j in range(m)
           if i != j and (Cr[i, j] != 0 or Gr[i, j] != 0)]
    _s, i, j = min(off)
    nnz = int(kept().Ap_obj[-1])
    Cr[i, j] = Gr[i, j] = 0.0
    a, b, db = _both(_seed(), Cr=Cr, Gr=Gr)
    assert a == b and _counts(db).get('radau_tc:served') == 1, db
    assert int(kept().Ap_obj[-1]) == nnz - 1


SCALES = (0.0, 1e-12, 1e-6, 1e-3, 0.05, 0.4, 3.0, 40.0, 1e3, 1e30, 1e150, 1e300)


@settings(deadline=None)
@given(data=st.data())
def test_drawn_seeds_steps_and_charges_are_the_pythons(data):
    _on()
    seed = _seed()
    n = len(seed[0])
    for j in range(3):
        sc = data.draw(st.sampled_from(SCALES))
        if sc:
            for i in data.draw(st.lists(st.integers(0, n - 1), min_size=1, max_size=4)):
                seed[j][i] += sc * data.draw(st.sampled_from((1.0, -1.0, 0.37, -2.5)))
        if data.draw(st.integers(0, 40)) == 0:
            seed[j][data.draw(st.integers(0, n - 1))] = data.draw(
                st.sampled_from((np.nan, np.inf, -np.inf, -0.0)))
    h0 = _cap()[1]['ctx'].h
    h = h0 * data.draw(st.sampled_from((1.0, 1.0, 1.0, 1e-6, 1e-3, 30.0, 1e6, 1e300)))
    qn = None
    if data.draw(st.integers(0, 5)) == 0:
        qn = np.array(_cap()[1]['ctx'].qn, dtype=float)
        qn[data.draw(st.integers(0, n - 1))] = data.draw(st.sampled_from(
            (1e-12, -3.0, 1e300, -DBL_MAX, 1e305)))
    a, b, _db = _both(seed, h=h, qn=qn)
    assert a == b


def test_drawn_seeds_reach_the_c_and_its_hand_backs():
    """Not vacuous: on a walk of drawn seeds the C serves most loops, over
    one iteration and several, and hands some back."""
    _on()
    rng = np.random.default_rng(5)
    served = handed = 0
    iters = set()
    for k in range(48):
        sc = (1e-9, 1e-3, 0.05, 0.5, 2.0, 1e3)[k % 6]
        seed = _seed({j: sc * rng.standard_normal(len(_seed()[0])) for j in range(3)})
        a, b, db = _both(seed)
        assert a == b
        c = _counts(db)
        if c.get('radau_tc:served'):
            served += 1
            iters.add(c.get('radau_tc:iterations'))
        else:
            handed += 1
    assert served >= 24 and handed >= 1 and len(iters) >= 2, (served, handed, iters)


## -- numpy's complex product, read ----------------------------------------------------

def test_numpys_complex_product_is_read_not_assumed():
    """`cmul_mode` names the form numpy's own product takes here, and that
    form is numpy's on fresh draws (scalar by vector, every position of
    short vectors, and a real vector cast); the samples it reads tell the
    two forms apart."""
    _on()
    mode = TC.cmul_mode()
    assert mode in (0, 1)
    a, b, _za, _tiny = TC._cmul_samples()
    fused = [math.fma(x.real, y.real, -(x.imag * y.imag)) for x, y in zip(a, b)]
    sep = [x.real * y.real - x.imag * y.imag for x, y in zip(a, b)]
    assert all(f != s for f, s in zip(fused, sep))
    rng = np.random.default_rng(99)
    for L in (1, 2, 3, 5, 6, 7, 9, 16, 17):
        for _ in range(40):
            x = complex(*(rng.standard_normal(2) * 10.0 ** rng.integers(-15, 15, size=2)))
            ys = (rng.standard_normal(L) + 1j * rng.standard_normal(L)) * 10.0 ** rng.integers(
                -15, 15, size=L)
            got = np.complex128(x) * ys
            for t in range(L):
                y = complex(ys[t])
                if mode:
                    want = (math.fma(x.real, y.real, -(x.imag * y.imag)),
                            math.fma(x.real, y.imag, x.imag * y.real))
                else:
                    want = (x.real * y.real - x.imag * y.imag, x.real * y.imag + x.imag * y.real)
                assert TC._same((float(got[t].real), float(got[t].imag)), want), (x, y)


def test_a_product_of_neither_form_turns_the_path_off(monkeypatch):
    _on()
    monkeypatch.setattr(TC, 'STATUS', TC.STATUS)
    monkeypatch.setitem(TC._MOD, 'cmul', None)
    monkeypatch.setattr(TC, '_driver', None)
    assert TC.driver() is None and 'neither form' in TC.STATUS
    monkeypatch.undo()
    assert TC.driver() is not None


## -- the declines ---------------------------------------------------------------------

def _declined(why, **kw):
    a, b, db = _both(_seed(), **kw)
    assert a == b
    return _counts(db).get('radau_tc:' + why, 0)


def test_the_declines(monkeypatch):
    """Before anything is touched: a patched piece of the loop, an instance
    shadow, a provided function, a stateful limiter, a bypass, a solver not
    analysed for the record, a singular kept LU -- each counted, the Python
    loop's answer."""
    _on()
    tr, got = _cap()
    ctx = got['ctx']
    reltol, abstol, maxit = got['tols']
    ## (the source memo as the loop found it: a finished run drops its own)
    tr.__dict__['_u_memo'] = dict(got['u'])
    fz = tr._radau_frozen(got['Cr'], got['Gr'], ctx.h)
    Y = _seed()

    def no(why, *a, **kw):
        before = _paths.snapshot()
        assert TC.solve(*a, **kw) is None
        assert _paths.since(before).get('radau_tc:' + why) == 1, _paths.since(before)
    args = (tr, ctx, fz, Y, got['src'])
    no('frozen', tr, ctx, None, Y, got['src'], None, [], True, reltol, abstol, maxit)
    no('pf', *args, lambda t: 0.0, [], True, reltol, abstol, maxit)
    no('stateful', *args, None, [object()], True, reltol, abstol, maxit)
    no('bypass', *args, None, [], False, reltol, abstol, maxit)
    ## a solver whose numeric is not of the record's pattern (fresh)
    tr._radau_zsolver = LS.ComplexKLUSolver()
    fz2 = tr._radau_frozen(got['Cr'], got['Gr'], ctx.h)
    no('klu', tr, ctx, fz2, Y, got['src'], None, [], True, reltol, abstol, maxit)
    ## a kept LU found singular
    fz2.zs.solve_prepared(fz2.prep, np.ones(fz2.m, dtype=complex))
    fz2.lu._info = 2
    no('lu', tr, ctx, fz2, Y, got['src'], None, [], True, reltol, abstol, maxit)
    fz2.lu._info = None
    ## a patched piece of the loop, on the class or on the instance
    real = _RadauStages._transform_loop
    monkeypatch.setattr(_RadauStages, '_transform_loop', lambda self, *a: real(self, *a))
    no('patched', tr, ctx, fz2, Y, got['src'], None, [], True, reltol, abstol, maxit)
    monkeypatch.undo()
    tr._transform_loop = tr._transform_loop
    try:
        no('shadow_tr', tr, ctx, fz2, Y, got['src'], None, [], True, reltol, abstol, maxit)
    finally:
        del tr._transform_loop
    monkeypatch.setattr(TC, 'ENABLED', False)
    no('off', tr, ctx, fz2, Y, got['src'], None, [], True, reltol, abstol, maxit)


def test_a_subclass_declines():
    _on()

    class T2(Transient):
        pass
    tr = T2(_stage(), toolkit=circuit.numeric, integrator=RadauIIA3Integrator())
    before = _paths.snapshot()
    tr.solve(tend=1e-7, timestep=2e-8, fixed_timestep=True)
    d = _paths.since(before)
    assert d.get('radau_tc:class', 0) > 0 and 'radau_tc:served' not in d, d
