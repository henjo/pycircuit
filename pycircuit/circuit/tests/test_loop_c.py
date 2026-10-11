"""The stepping loop in C (`_tran_loop_c`, the C stepping loop plan's
stages 4 and 5; 2026-10-11): a fixed run with the C loop on is the run with
it off -- every output byte, time, statistic and warning, the `_paths`
tallies the Python steps make, and every attribute the steps leave on the
transient (`transient._RUN_STATE`, its memos and counters aside) -- on the
chains of the speed cases (MOS level 1, Gummel-Poon on the chord, the PSP
stage behind the branch screen), a sparse run, a pulse's landings and
order drops, the trapezoid and Euler, a grid that does not divide `tend`,
an entry that fills its outputs; through its hand-backs (a screen that
fires, a Newton that fails); its declines; and alive: a seed a hair off is
another run.
"""
import math
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_loop_c, circuit, compact
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.elements import VS, C, R, SubCircuit, VPulse, VSin, gnd
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    TrapezoidalIntegrator,
)
from pycircuit.circuit.transient import Transient

pytestmark = pytest.mark.skipif(_tran_loop_c.driver() is None,
                                reason='the fixed-step loop unit did not build')


def _chain(make, n=8, vdd=1.8, vg=0.9, va=2e-2, cap=None):
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vg'] = VSin('g0', gnd, v=vg, va=va, freq=1e6)
    for k in range(n):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = make(f'd{k}', 'g0' if k == 0 else f'd{k - 1}')
        if cap:
            c[f'c{k}'] = C(f'd{k}', gnd, c=cap)
    return c


def _mos(n=8, **kw):
    return _chain(lambda d, g: eh.MosLevel1Hdl(d, g, gnd, gnd), n=n, **kw)


def _gp(n=6):
    c = SubCircuit()
    c['vcc'] = VS('vcc', gnd, v=3.0)
    c['vin'] = VSin('in', gnd, v=0.75, va=2e-2, freq=1e6)
    prev = 'in'
    for k in range(n):
        c[f'rb{k}'] = R(prev, f'b{k}', r=1e4)
        c[f'rc{k}'] = R('vcc', f'c{k}', r=1e3)
        c[f'Q{k}'] = eh.GummelPoonNpnHdl(f'c{k}', f'b{k}', gnd)
        prev = f'c{k}'
    return c


def _psp_stage():
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd, fnt=1.0)
    return c


def _pulsed(n=6):
    """A MOS chain behind a pulse: breakpoints, landings, order drops."""
    c = _mos(n, cap=1e-13)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=3e-9, tr=1e-9, tf=1e-9, pw=6e-9, per=2e-8)
    return c


## (what the steps leave that is not state they compute from: the
## predictor's memos, the core's check counters, a cache keyed by ids)
_NOT_STATE = {'_pred_fast', '_pred_wmemo', '_tran_core', '_row_names_cache'}
## (counts the C loop's own machinery makes: its build, its table's)
_OWN_COUNTS = ('once:', 'sources_c:', 'loop_c:', 'pred_c:', 'branch_c:bail')


def _norm(v, depth=0):
    if isinstance(v, np.ndarray):
        return ('nd', v.dtype.str, v.shape, v.tobytes())
    if isinstance(v, (list, tuple)) and depth < 4:
        return (type(v).__name__, tuple(_norm(e, depth + 1) for e in v))
    if isinstance(v, dict) and depth < 4:
        return ('dict', tuple(sorted((str(k), _norm(e, depth + 1)) for k, e in v.items())))
    if isinstance(v, (int, float, str, bool, type(None), np.floating, np.integer)):
        return (type(v).__name__, repr(v))
    if type(v).__module__ == 'pycircuit.circuit.integrator':
        return ('obj', type(v).__name__, _norm(vars(v), depth + 1))
    return ('obj', type(v).__name__)


def _run(build, on, monkeypatch, make=None, setup=None, **kw):
    monkeypatch.setattr(_tran_loop_c, 'ENABLED', on)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **(make or {}))
        if setup:
            setup(tr)
        kw.setdefault('fixed_timestep', True)
        res = tr.solve(**kw)
    st = tr.statistics
    stats = {k: getattr(st, k) for k in st.__slots__ if 'seconds' not in k}
    counts = {k: v for k, v in _paths.since(before).items() if not k.startswith(_OWN_COUNTS)}
    state = {k: _norm(tr.__dict__.get(k, '<absent>')) for k in trmod._RUN_STATE
             if k not in _NOT_STATE}
    ## (the linear solver's own counts: a sparse run's factors and numerics)
    state['linear solver'] = repr(tr.__dict__.get('_default_linearsolver'))
    ## (the controller's own: its running reference, its retry record, its
    ## carried order, its last error -- its C error test's context aside)
    ctrl = getattr(tr, 'step_controller', None)
    if ctrl is not None:
        state['controller'] = _norm({k: v for k, v in vars(ctrl).items() if k != '_lte_c'})
    return {'x': np.asarray(res.x, float).tobytes(),
            't': np.asarray(res.sweep_values, float).tobytes(), 'stats': stats,
            'warnings': sorted(str(w.message) for w in W), 'counts': counts,
            'state': state}, _paths.since(before)


def _same(build, monkeypatch, steps_min=1, **kw):
    a, ca = _run(build, True, monkeypatch, **dict(kw))
    b, _cb = _run(build, False, monkeypatch, **dict(kw))
    for k in ('x', 't', 'stats', 'warnings', 'counts', 'state'):
        if k == 'state':
            diff = [n for n in a[k] if a[k][n] != b[k][n]]
            assert not diff, ('the state moved', diff)
        elif k == 'counts':
            diff = {n: (a[k].get(n), b[k].get(n)) for n in set(a[k]) | set(b[k])
                    if a[k].get(n) != b[k].get(n)}
            assert not diff, (('the counts the C loop makes for the Python attempts\' calls '
                               '(`_tran_loop_c._count`) moved -- (C, Python):'), diff)
        else:
            assert a[k] == b[k], (k, a[k] if k != 'x' else '', b[k] if k != 'x' else '')
    assert ca.get('loop_c:steps', 0) >= steps_min, ca
    return ca


@pytest.mark.parametrize('name, build, make, kw', [
    ('mos1', _mos, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('gp-chord', _gp, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('psp-screen', _psp_stage, None, {'tend': 60 * 2e-8, 'timestep': 2e-8}),
    ('sparse', lambda: _mos(250, cap=1e-13), None, {'tend': 20 * 2e-10, 'timestep': 2e-10}),
    ('pulse', _pulsed, None, {'tend': 4e-8, 'timestep': 2.5e-10}),
    ('trap', _pulsed, {'integrator': TrapezoidalIntegrator()},
     {'tend': 3e-8, 'timestep': 2.5e-10}),
    ('euler', _pulsed, {'integrator': EulerIntegrator()}, {'tend': 3e-8, 'timestep': 2.5e-10}),
    ('ragged-end', _mos, None, {'tend': 61.37 * 2e-8, 'timestep': 2e-8}),
])
def test_a_fixed_run_is_the_python_one(name, build, make, kw, monkeypatch):
    ca = _same(build, monkeypatch, steps_min=10, make=make, **kw)
    if name == 'pulse':
        assert ca.get('loop_c:steps', 0) > 100


def test_an_entry_that_fills_its_outputs_comes_back(monkeypatch):
    """Seven rows an entry: the run in many entries, the same run."""
    monkeypatch.setattr(_tran_loop_c, 'CAP_MAX', 7)
    _same(_mos, monkeypatch, steps_min=40, tend=60 * 2e-8, timestep=2e-8)


def test_a_screen_that_fires_hands_the_step_back(monkeypatch):
    """The PSP stage driven hard (its C moves) and the screen's tolerance
    raised on the instance, so it fires where C dips below the run's
    largest: those steps go to the Python (its SVD, its confirmation), the
    rest stay in C -- the same run, the same `branch_screens`."""
    def build():
        c = _psp_stage()
        c['vg'] = VSin('g', gnd, v=0.7, va=0.5, freq=1e6)
        return c

    def setup(tr):
        tr.BRANCH_SCREEN_TOL = 0.9
    ca = _same(build, monkeypatch, steps_min=5, setup=setup, tend=80 * 2e-8, timestep=2e-8)
    assert ca.get('loop_c:bail:screen', 0) > 0


def test_a_newton_that_fails_hands_the_step_back(monkeypatch):
    """Two Newton iterations at most: the steps that need more fail in C and
    are the Python's (its failure, its smaller retry) -- the same run."""
    def setup(tr):
        tr.par.maxiter = 2
    ca = _same(_pulsed, monkeypatch, steps_min=5, setup=setup, tend=2e-8, timestep=2.5e-10)
    assert ca.get('loop_c:bail:newton', 0) > 0


def test_the_declines(monkeypatch):
    """An element with an `accept_step` hook, a patched piece of the step,
    a shadowed method, a source perturbed on the circuit: the Python loop,
    every step."""
    class Hooked(R):
        def accept_step(self, t, x, epar):
            pass

    def hooked():
        c = _mos(4)
        c['rh'] = Hooked('d0', gnd, r=1e6)
        return c
    _, ca = _run(hooked, True, monkeypatch, tend=20 * 2e-8, timestep=2e-8)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:circuit', 0) > 0
    real = Transient._pred_note
    monkeypatch.setattr(Transient, '_pred_note', lambda self, *a, **k: real(self, *a, **k))
    _, ca = _run(_mos, True, monkeypatch, tend=20 * 2e-8, timestep=2e-8)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:patched', 0) > 0
    monkeypatch.undo()

    def shadow(tr):
        tr._q_at = lambda x: Transient._q_at(tr, x)
    _, ca = _run(_mos, True, monkeypatch, setup=shadow, tend=20 * 2e-8, timestep=2e-8)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:shadow', 0) > 0

    def perturbed():
        ## (a source perturbed on the instance, as the shooting tests make
        ## one: the Python steps call it, so the table cannot stand in)
        c = _mos(4)
        real = c.u
        c.u = lambda t, epar=None, analysis=None, **kw: np.asarray(
            real(t, epar, analysis=analysis, **kw), float) * 1.01
        return c
    _, ca = _run(perturbed, True, monkeypatch, tend=20 * 2e-8, timestep=2e-8)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:circuit', 0) > 0


def test_a_seed_a_hair_off_is_another_run(monkeypatch):
    """Alive: the predictor's clamp a ulp off in the C (its seeds move, and
    the converged points with them) -- the comparison sees another run."""
    a, _ = _run(_pulsed, False, monkeypatch, tend=2e-8, timestep=2.5e-10)
    bind = _tran_loop_c._bind

    def planted(*args, **kw):
        run = bind(*args, **kw)
        run.pr[0].clamp = math.nextafter(run.pr[0].clamp, 0.0)
        return run
    monkeypatch.setattr(_tran_loop_c, '_bind', planted)
    b, cb = _run(_pulsed, True, monkeypatch, tend=2e-8, timestep=2.5e-10)
    assert cb.get('loop_c:steps', 0) > 10 and a['x'] != b['x']


## -- adaptive (stage 5) -------------------------------------------------------

def _hard_mos(n=8):
    """A MOS chain with capacitors driven hard: the controller rejects."""
    return _mos(n, va=0.6, cap=1e-13)


@pytest.mark.parametrize('name, build, make, kw', [
    ('gear', _hard_mos, None, {'tend': 4e-6, 'timestep': 2e-8}),
    ('pulse', _pulsed, None, {'tend': 6e-8, 'timestep': 1e-9}),
    ('trap', _pulsed, {'integrator': TrapezoidalIntegrator()}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('gp-chord', _gp, None, {'tend': 3e-6, 'timestep': 2e-8}),
    ('psp-screen', _psp_stage, None, {'tend': 3e-6, 'timestep': 2e-8}),
    ('sparse', lambda: _mos(250, cap=1e-13), {'reltol': 1e-5}, {'tend': 1e-8, 'timestep': 2e-10}),
    ('alllocal', _pulsed, {'relref': 'alllocal'}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('pointlocal', _pulsed, {'relref': 'pointlocal'}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('assumed', _pulsed, {'retry_order': 'assumed'}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('carried', _pulsed, {'retry_order': 'carried'}, {'tend': 4e-8, 'timestep': 1e-9}),
    ('mask', _pulsed, {'lte_algebraic': 'skip'}, {'tend': 4e-8, 'timestep': 1e-9}),
])
def test_an_adaptive_run_is_the_python_one(name, build, make, kw, monkeypatch):
    """Every attempt -- accepted, rejected at a point and retried, a landing,
    an order drop -- the same bytes, statistics, warnings, tallies, state and
    controller (its running reference, retry record, carried order)."""
    kw = dict(kw, fixed_timestep=False)
    ca = _same(build, monkeypatch, steps_min=10, make=make, **kw)
    if build is _pulsed:
        ## (the pulse's corners reject: retried at the point, in C)
        assert ca.get('loop_c:attempts', 0) > ca.get('loop_c:steps', 0), ca


def test_a_piece_shadowed_on_the_transient_is_declined(monkeypatch):
    """A piece the attempts call, shadowed on the transient instance -- here
    the fast predictor's node choice, which the Python's attempts then call
    -- declines (`loop_c:shadow`): the C would not call it."""
    def shadow(tr):
        tr._pred_fast_nodes = lambda hist: None
    ca = _same(_mos, monkeypatch, steps_min=0, setup=shadow, tend=1e-6, timestep=2e-8)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:shadow', 0) > 0, ca


def test_a_predictor_off_or_unbuilt_is_declined(monkeypatch):
    """The predictor's switch off, or its unit not built: the loop declines
    (`loop_c:pred`), as for the screen and the law -- the Python attempts,
    nothing raised."""
    from pycircuit.circuit import _tran_pred_c
    kw = {'tend': 1e-6, 'timestep': 2e-8}
    monkeypatch.setattr(_tran_pred_c, 'ENABLED', False)
    _, ca = _run(_mos, True, monkeypatch, **dict(kw))
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:pred', 0) > 0
    monkeypatch.setattr(_tran_pred_c, 'ENABLED', True)
    monkeypatch.setattr(_tran_pred_c, 'driver', lambda: None)
    _, ca = _run(_mos, True, monkeypatch, **dict(kw))
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('loop_c:pred', 0) > 0


def test_the_adaptive_declines(monkeypatch):
    """Another controller (a subclass: its law may be another) -- the
    Python attempts, every one -- and an excursion bound, which the judge's
    check refuses."""
    from pycircuit.circuit import stepcontroller as sc

    class Mine(sc.IntegralController):
        pass

    def mine(tr):
        tr.step_controller = Mine()
    _, ca = _run(_hard_mos, True, monkeypatch, setup=mine, tend=1e-6, timestep=2e-8,
                 fixed_timestep=False)
    ## (declined where the C error test declines it: no serve to enter on)
    assert ca.get('loop_c:steps', 0) == 0 and ca.get('lte_c:controller', 0) > 0
    ## (an excursion bound, asked of the judge's check directly: a run with
    ## one set crawls in the Python loop on this chain -- a matter of its
    ## own, 2026-10-11)
    from types import SimpleNamespace
    tr = Transient(_hard_mos(), toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(tend=2e-7, timestep=2e-8)
    fam = SimpleNamespace(max_retries=10, max_reject=3)
    loop = SimpleNamespace(run=SimpleNamespace(abstol=1e-6, max_step=1e-8), minstep=1e-18,
                           family=fam)
    O = _tran_loop_c._own()
    assert type(_tran_loop_c._judge(loop, tr, O, tr.cir.n)) is tuple
    tr.par.max_dv_step = 5.0
    assert _tran_loop_c._judge(loop, tr, O, tr.cir.n) == 'veto'


def test_a_law_a_hair_off_is_another_adaptive_run(monkeypatch):
    """Alive: the shrink's floor a ulp off in the C's law (its rejections
    retry at another step) -- the comparison sees another run."""
    kw = {'tend': 6e-8, 'timestep': 1e-9, 'fixed_timestep': False}
    a, _ = _run(_pulsed, False, monkeypatch, **dict(kw))
    bind = _tran_loop_c._bind

    def planted(*args, **k):
        run = bind(*args, **k)
        run.ct.min_shrink = math.nextafter(run.ct.min_shrink, 1.0)
        run.ct.order_min = 1.5
        return run
    monkeypatch.setattr(_tran_loop_c, '_bind', planted)
    b, cb = _run(_pulsed, True, monkeypatch, **dict(kw))
    assert cb.get('loop_c:attempts', 0) > 10 and (a['x'] != b['x'] or a['t'] != b['t'])


#: (the functions whose returned strings are decline reasons)
_REASONERS = ('_judge', '_why_not_the_run', '_sources', '_read')


def _reasons(text):
    """Every reason the loop's text declines with: `_declined(loop, '<x>')`
    and the strings the reasoners return (`_REASONERS`)."""
    import ast
    tree = ast.parse(text)
    out = set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Call) and getattr(node.func, 'id', None) == '_declined'
                and len(node.args) > 1 and isinstance(node.args[1], ast.Constant)):
            out.add(node.args[1].value)
        if isinstance(node, ast.FunctionDef) and node.name in _REASONERS:
            out |= {r.value.value for r in ast.walk(node) if isinstance(r, ast.Return)
                    and isinstance(r.value, ast.Constant) and type(r.value.value) is str}
    return out


def test_every_decline_says_whether_it_holds_for_the_run():
    """`_KEEPS` names every reason the loop declines with (a reason missing
    would raise at the decline) and no other."""
    import inspect
    found = _reasons(inspect.getsource(_tran_loop_c))
    assert found == set(_tran_loop_c._KEEPS), (found ^ set(_tran_loop_c._KEEPS))
    ## alive: a reason the table does not name, and a judge's
    planted = "def run(loop):\n    return _declined(loop, 'nowhere')\n" \
              "def _judge(loop):\n    return 'elsewhere'\n"
    assert _reasons(planted) == {'nowhere', 'elsewhere'}


def test_the_shrink_guard_is_one_value():
    """loop.c's `SHRINK_GUARD_RATIO` is `integrator.SHRINK_GUARD_RATIO` (a
    float: not an enum the glue could read, so held here)."""
    import re

    from pycircuit.circuit import _cbuild, integrator
    m = re.search(r'^#define SHRINK_GUARD_RATIO (\S+)$', _cbuild.csrc('loop.c'), re.MULTILINE)
    assert m and float(m.group(1)) == integrator.SHRINK_GUARD_RATIO
