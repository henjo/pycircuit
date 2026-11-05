"""The error test on the seed's factors (`Transient(lte_jacobian='seed')`):
an adaptive Gear-2 attempt whose sparse chord (KLU) converged within
`_tran_newton.LTE_SEED_ITERS` iterations takes J^-1 Eg by the chord's own
numeric of the seed's Jacobian -- its converged point evaluated for C and q
alone -- where 'converged' (the default) evaluates and factors J there.  Off
by default; the Python step and chain, the C error test, the C Newton
served from Python and the C loop alike."""
import sys
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _paths,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton,
    _tran_newton_c,
    circuit,
)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.integrator import RadauIIA3Integrator, TrapezoidalIntegrator
from pycircuit.circuit.stepcontroller import StepController
from pycircuit.circuit.transient import Transient

#: the paths a run may take, as (the C Newton, the C loop, the C error test)
PATHS = {'python': (False, False, False), 'lte_c': (False, False, True),
         'served': (True, False, True), 'loop': (True, True, True)}


def _chain(n=250, fast=False):
    """A pulsed MOS chain with capacitors: at 250 stages 254 unknowns, KLU
    and the chord ('auto').  `fast`: edges of half a nanosecond into
    heavier stages -- error tests that raise floating-point flags (handed
    back), solves past every bound."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    if fast:
        c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=1e-9, tr=5e-10, tf=5e-10, pw=4e-9, per=1e-8)
    else:
        c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-8, tr=1e-8, tf=1e-8, pw=2.4e-7, per=5e-7)
    rl, cap = (5e3, 1e-13) if fast else (2e4, 2e-14)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=rl)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=cap)
    return c


def _run(monkeypatch, lte='seed', path='loop', tend=3e-7, step=1e-8, cir=None, rt=None, solve=None,
         errs=None, **make):
    newton_c, loop_c, lte_c = PATHS[path]
    monkeypatch.setattr(_tran_newton_c, 'ENABLED', newton_c)
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', loop_c)
    monkeypatch.setattr(_tran_lte_c, 'ENABLED', lte_c)
    before = _paths.snapshot()
    with warnings.catch_warnings(), np.errstate(**(errs or {})):
        warnings.simplefilter('ignore')
        tr = Transient(cir if cir is not None else _chain(), toolkit=circuit.numeric,
                       lte_jacobian=lte, **make)
        if rt is not None:
            tr.par.reltol = rt
        res = tr.solve(tend=tend, timestep=step, **(solve or {}))
    return tr, res, _paths.since(before)


def _key(tr, res):
    """The run's bytes, statistics and KLU's counts."""
    st = tr.statistics
    ls = tr._get_linearsolver()._choice
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st.accepted_steps, st.rejected_steps, st.newton_iterations, type(ls).__name__,
            *(getattr(ls, k, None) for k in ('factors', 'refactors', 'chord_refactors',
                                             'chord_factors', 'chord_pivot_fallbacks',
                                             'chord_residual_fallbacks')))


def _bound(monkeypatch, k):
    monkeypatch.setattr(_tran_newton, 'LTE_SEED_ITERS', k)
    monkeypatch.setitem(_tran_newton_c._mods(), 'lte_seed', k)


def test_off_by_default_and_checked(monkeypatch):
    """'converged' by default: nothing held, the converged point's G made
    for every judged attempt; another value raises at the solve."""
    tr, _res, d = _run(monkeypatch, lte='converged', path='python', tend=2e-7)
    assert Transient(_chain(2), toolkit=circuit.numeric).par.lte_jacobian == 'converged'
    assert d.get('lte:seed', 0) == 0 and d.get('core.c:served', 0) == 0 and d['core.j:served'] > 50
    tr.par.lte_jacobian = 'stale'
    with pytest.raises(ValueError, match='lte_jacobian'):
        tr.solve(tend=1e-8, timestep=1e-9)


@pytest.mark.parametrize('bound', [2, 3])
def test_every_path_is_the_python_one(monkeypatch, bound):
    """The Python step and chain, the C error test, the C Newton served from
    Python and the C loop make the same run -- bytes, statistics, KLU's
    counts, the held tests and the error tests handed back on a flag -- at
    a bound holding some of the solves (2) and most (3); each C path
    serving its share, the loop its own full tests after held ones (2), and
    entering before the run's KLU numeric exists.  Under
    `np.errstate(under='warn')`: the fast chain's flags are underflows, which
    under numpy's default are no flag (test_lte_underflow)."""
    _bound(monkeypatch, bound)
    bound_runs, bound_ = [], _tran_loop_c._bound

    def caught(*a, **k):
        r = bound_(*a, **k)
        bound_runs.append(r)
        return r
    monkeypatch.setattr(_tran_loop_c, '_bound', caught)
    runs = {p: _run(monkeypatch, path=p, cir=_chain(fast=True), tend=2.4e-8, step=2e-10,
                    errs={'under': 'warn'})
            for p in PATHS}
    keys = {p: _key(tr, res) for p, (tr, res, _d) in runs.items()}
    assert all(keys[p] == keys['python'] for p in keys)
    d = {p: r[2] for p, r in runs.items()}
    held = d['python']['lte:seed']
    tried = sum(d['python'].get(k, 0) for k in ('core.j:served', 'core.c:served'))
    assert 0 < held < tried
    assert all(d[p]['lte:seed'] == held for p in d)
    assert d['lte_c']['lte_c:held'] == d['served']['lte_c:held'] == d['loop']['lte_c:held'] > 0
    assert (d['lte_c']['lte_c:bail:flags'] == d['served']['lte_c:bail:flags']
            == d['loop']['lte_c:bail:flags'] > 0)
    assert d['served']['newton_c:lte_seed'] > 0 and d['loop']['loop_c:attempts'] > 0
    assert bound == 3 or d['loop']['loop_c:lte_on_pattern'] > 0
    ## (the loop never declined for the run's KLU numeric)
    assert d['loop'].get('loop_c:sparse', 0) == 0
    ## (the factors in hand are a call's alone: every error test's context
    ## left without them -- the Python's C error test's, the C loop's)
    assert all(r[0].step_controller._lte_c.s.lu_mode == 0 for p, r in runs.items() if p != 'python')
    assert bound_runs and all(r.lte.s.lu_mode == 0 for r in bound_runs if r is not None)


def test_the_error_test_takes_the_chords_factors_of_the_seed(monkeypatch):
    """Each held error test solves with the Newton's own factors --
    `_Factored` over the chord's `fac` -- and no J made (the early abort's
    estimates on the same factors, with the seed's J, as before); every
    other test makes its J; the converged point of each held attempt
    evaluated without G."""
    seen = []
    charge_lte = StepController._charge_lte

    def spy(self, s):
        seen.append((type(s.linsolver), s.J is None, type(getattr(s.linsolver, 'fac', None)),
                     sys._getframe(1).f_code.co_name))
        return charge_lte(self, s)
    monkeypatch.setattr(StepController, '_charge_lte', spy)
    _tr, _res, d = _run(monkeypatch, path='python')
    tests = [x for x in seen if x[3] == '_max_error']
    held = [x for x in tests if x[0] is _tran_newton._Factored]
    assert 0 < len(held) <= d['lte:seed']
    assert all(no_j and fac.__name__ == '_Factored' for _c, no_j, fac, _w in held)
    assert all(not no_j for c, no_j, _f, _w in tests if c is not _tran_newton._Factored)
    assert all(not no_j for _c, no_j, _f, w in seen if w == '_early_estimate')
    assert d['core.c:served'] == d['lte:seed']


def test_a_full_test_without_the_run_numeric_is_handed_back(monkeypatch):
    """The C loop meeting a full error test with no run numeric to refactor
    -- presented here at every entry -- hands that attempt back by its own
    test (`numeric`; KLU would refuse the NULL numeric too, as an error
    test's bail), and the Python's repeat makes the run the plain one.  At
    bound 2, where the loop makes full tests of its own."""
    _bound(monkeypatch, 2)
    plain = _run(monkeypatch, path='loop', cir=_chain(fast=True), tend=2.4e-8, step=2e-10)
    load = _tran_loop_c._load

    def numeric_missing(loop, tr, td, run_, *a, **k):
        r = load(loop, tr, td, run_, *a, **k)
        if run_.sk is not None:
            run_.lte.s.klu_num = run_.lte.ffi.NULL
        return r
    monkeypatch.setattr(_tran_loop_c, '_load', numeric_missing)
    got = _run(monkeypatch, path='loop', cir=_chain(fast=True), tend=2.4e-8, step=2e-10)
    assert _key(*got[:2]) == _key(*plain[:2])
    assert got[2]['loop_c:bail:numeric'] > 0 and got[2]['loop_c:attempts'] > 20


class _Stale:
    """Held factors no solve made: solving with them fails loudly."""

    def solve(self, _J, _b, _toolkit):
        raise AssertionError('a stale factorisation reached the error test')


@pytest.mark.parametrize('path', ['python', 'served'])
def test_a_solve_drops_held_factors_it_did_not_make(monkeypatch, path):
    """Every Newton solve -- `_newton`, the C Newton -- starts from no held
    factors: one planted before each solve never reaches an error test, and
    the run is the plain one."""
    _bound(monkeypatch, 2)
    plain = _key(*_run(monkeypatch, path=path)[:2])
    if path == 'python':
        newton = Transient._newton

        def planted(self, *a, **k):
            self._lte_fac = _Stale()
            return newton(self, *a, **k)
        monkeypatch.setattr(Transient, '_newton', planted)
    else:
        solve = _tran_newton_c.solve

        def planted(tr, *a, **k):
            tr._lte_fac = _Stale()
            return solve(tr, *a, **k)
        monkeypatch.setattr(_tran_newton_c, 'solve', planted)
    assert _key(*_run(monkeypatch, path=path)[:2]) == plain


@pytest.mark.parametrize('case', ['chord_off', 'dense', 'fixed', 'trapezoid', 'radau'])
def test_inert_where_nothing_holds(monkeypatch, case):
    """'seed' makes the run 'converged' makes -- bytes and counts, nothing
    held -- with the chord off, on a dense chord (a small chain), on a fixed
    grid, and for the trapezoid and Radau (Gear-2 runs alone)."""
    kw = {'chord_off': {'chord_jacobian': False},
          'dense': {'cir': _chain(8), 'chord_jacobian': True},
          'fixed': {'solve': {'fixed_timestep': True}},
          'trapezoid': {'integrator': TrapezoidalIntegrator()},
          'radau': {'integrator': RadauIIA3Integrator(), 'tend': 1e-7}}[case]

    def kws():
        out = dict(kw)
        if 'integrator' in out:
            out['integrator'] = type(out['integrator'])()
        if 'cir' in out:
            out['cir'] = _chain(8)
        return out
    a = _run(monkeypatch, lte='converged', **kws())
    b = _run(monkeypatch, lte='seed', **kws())
    assert _key(*a[:2]) == _key(*b[:2])
    assert b[2].get('lte:seed', 0) == 0


def test_the_loop_holds_its_record_to_the_runs_mark(monkeypatch):
    """The C loop enters only where the C Newton's record was armed as the
    run's mark asks (`lte_seed_for`): a record armed otherwise -- here the
    served solves arming nothing -- keeps it out ('newton')."""
    real = _tran_newton_c.lte_seed_for

    def skewed(td, M, chord):
        if sys._getframe(1).f_code.co_name == '_why_not_the_run':
            return real(td, M, chord)
        return 0
    monkeypatch.setattr(_tran_newton_c, 'lte_seed_for', skewed)
    _tr, _res, d = _run(monkeypatch, path='loop')
    assert d.get('loop_c:attempts', 0) == 0 and d.get('loop_c:newton', 0) > 0


def _err(ref, res, nodes=('d0', 'd1', 'd2', 'd3', 'd4', 'd5')):
    out = []
    for nd in nodes:
        wr, w = ref.v(nd), res.v(nd)
        tr_, yr = np.asarray(wr.x[0], float), np.asarray(wr.y, float)
        y = np.interp(tr_, np.asarray(w.x[0], float), np.asarray(w.y, float))
        out.append(np.max(np.abs(y - yr)) / (yr.max() - yr.min()))
    return max(out)


def test_a_seed_run_is_within_its_tolerance(monkeypatch):
    """Off and on against a reltol/100 reference over a microsecond: the
    answers move, and on the error of the first six stages (their swing
    full) is at most 1.25x the off run's + reltol."""
    _t, ref, _d = _run(monkeypatch, lte='converged', rt=1e-6, tend=1e-6)
    t0, off, _d = _run(monkeypatch, lte='converged', tend=1e-6)
    t1, on, d = _run(monkeypatch, lte='seed', tend=1e-6)
    assert d['lte:seed'] > 50 and _key(t0, off) != _key(t1, on)
    assert _err(ref, on) <= 1.25 * _err(ref, off) + t1.par.reltol
