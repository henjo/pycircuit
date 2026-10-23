"""The multistep Newton solve in C (`_tran_newton_c`, speed round 8, stage
4): bit for bit `_newton` and `jacobian_only`, its declines before anything
is touched, its bails without a trace."""
import warnings

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pycircuit.circuit import _paths, _tran_newton_c, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import (
    VS,
    C,
    Diode,
    R,
    SubCircuit,
    VSin,
    gnd,
)
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    Gear2Integrator,
    TrapezoidalIntegrator,
)
from pycircuit.circuit.transient import Transient


def _mos_chain(n=4, va=0.6, rl=5e3, cap=None):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=0.9, va=va, freq=1e6)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=rl)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', f'g{k}' if k == 0 else f'd{k - 1}', gnd, gnd)
        if cap:
            c[f'c{k}'] = C(f'd{k}', gnd, c=cap)
    return c


def _gp_chain(n=3, va=0.05):
    c = SubCircuit()
    c.add_node('vcc')
    c['vcc'] = VS('vcc', gnd, v=3.0)
    c.add_node('in')
    c['vin'] = VSin('in', gnd, v=0.75, va=va, freq=1e6)
    prev = 'in'
    for k in range(n):
        c.add_node(f'b{k}')
        c.add_node(f'c{k}')
        c[f'rb{k}'] = R(prev, f'b{k}', r=1e4)
        c[f'rc{k}'] = R('vcc', f'c{k}', r=1e3)
        c[f'Q{k}'] = eh.GummelPoonNpnHdl(f'c{k}', f'b{k}', gnd)
        prev = f'c{k}'
    return c


def _run(build, on, **kw):
    """`(x bytes, time bytes, statistics, warnings, path counts)` of one
    transient with the C solve on or off."""
    old = _tran_newton_c.ENABLED
    _tran_newton_c.ENABLED = on
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings(record=True) as W:
            warnings.simplefilter('always')
            ## (a finalizer's warning is not the run's: a file another test left in
            ## a reference cycle closes whenever the collector runs -- inside either
            ## run, 2026-10-04)
            warnings.simplefilter('ignore', ResourceWarning)
            tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
            res = tr.solve(**kw)
    finally:
        _tran_newton_c.ENABLED = old
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__
           if 'seconds' not in k} if hasattr(tr.statistics, '__slots__') else \
        {k: v for k, v in vars(tr.statistics).items() if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _paths.since(before))


def _same(build, served_min=1, **kw):
    a = _run(build, False, **dict(kw))
    b = _run(build, True, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert b[4].get('newton_c:served', 0) >= served_min, b[4]
    return a, b


@pytest.mark.parametrize('fixed', [True, False])
def test_the_plain_newton_is_the_python_one(fixed):
    """MOS level 1 (no chord): fixed steps and the default adaptive gear."""
    kw = {'tend': 2e-6, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    _same(_mos_chain, **kw)


@pytest.mark.parametrize('fixed', [True, False])
def test_the_chord_is_the_python_one(fixed):
    """Gummel-Poon with the chord (SciPy's LU, called from C)."""
    kw = {'tend': 2e-6, 'timestep': 2e-8, 'make': {'chord_jacobian': True}}
    if fixed:
        kw['fixed_timestep'] = True
    a, b = _same(_gp_chain, **kw)
    assert b[4].get('chord:chosen', 0) == a[4].get('chord:chosen', 0) > 0


@pytest.mark.parametrize('integ', [Gear2Integrator, EulerIntegrator, TrapezoidalIntegrator])
def test_every_companion_and_a_capacitive_chain(integ):
    """The core's companion enum, on a chain with capacitors (the branch
    screen reads the converged point's `C` from the C solve)."""
    _same(lambda: _mos_chain(cap=1e-13), tend=1.5e-6, timestep=2e-8,
          make={'integrator': integ()})


@settings(max_examples=12, deadline=None)
@given(n=st.integers(1, 5), va=st.floats(0.05, 0.9), rl=st.floats(1e3, 2e4),
       cap=st.sampled_from([None, 1e-14, 1e-12]))
def test_drawn_chains_are_the_python_ones(n, va, rl, cap):
    _same(lambda: _mos_chain(n, va, rl, cap), served_min=0, tend=1e-6, timestep=2e-8)


def test_a_bail_leaves_no_trace(monkeypatch):
    """Every C solve runs to its end and then hands the step back (its
    entry wrapped to report a non-finite iterate): the Python Newton
    repeats each one -- the same answer, statistics, warnings, and every
    count but the bail's own (the source memo's entry and counts rolled
    back).  (Not `maxiter`: that one is raised from the C.)"""
    drv = _tran_newton_c.driver()
    if drv is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    ffi, cfn, lapack = drv
    ## (once first, unmeasured: a process's first C build counts its
    ## `once:` keys in whichever run comes first)
    _run(_mos_chain, False, tend=6e-7, timestep=2e-8)
    a = _run(_mos_chain, False, tend=6e-7, timestep=2e-8)

    def bailing(s):
        cfn(s)
        return 2
    monkeypatch.setattr(_tran_newton_c, '_driver', (ffi, bailing, lapack))
    b = _run(_mos_chain, True, tend=6e-7, timestep=2e-8)
    assert a[:4] == b[:4]
    assert b[4].get('newton_c:served', 0) == 0 and b[4].get('newton_c:bail:x', 0) > 10
    rest = {k: v for k, v in b[4].items() if not k.startswith('newton_c:')}
    assert rest == {k: v for k, v in a[4].items() if not k.startswith('newton_c:')}


def test_the_iteration_limit_is_raised_from_the_c():
    """A full Newton stopped at `maxiter` with every iterate finite raises
    StandardNewton's exception from the C's last iteration -- no Python
    repeat -- and the run is the Python's: the same answer, statistics and
    warnings, the C counting `newton_c:maxiter` and handing nothing back.
    A `maxiter` of 1 for the transient alone (its operating point given):
    steps fail, and the controller retries them smaller."""
    from pycircuit.circuit.dcanalysis import DC
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    x0 = DC(_mos_chain()).solve().x
    kw = {'tend': 6e-7, 'timestep': 2e-8, 'x0': x0, 'make': {'maxiter': 1}}
    a = _run(_mos_chain, False, **dict(kw))
    b = _run(_mos_chain, True, **dict(kw))
    assert a[:4] == b[:4]
    assert b[4].get('newton_c:maxiter', 0) > 0, b[4]
    assert b[4].get('newton_c:bail:maxiter', 0) == 0, b[4]


def test_the_declines():
    """A stateful limiter (`Diode`) and a provided function: today's
    Newton, counted -- the circuit's reason once, kept with its plan."""
    def diode():
        c = _mos_chain(2)
        c['D'] = Diode('d0', gnd)
        return c
    _, b = _same(diode, served_min=0, tend=4e-7, timestep=2e-8)
    assert b[4].get('newton_c:stateful', 0) > 0 and b[4].get('newton_c:served', 0) == 0
    before = _paths.snapshot()
    Transient(_mos_chain(2), toolkit=circuit.numeric, uic=True).solve(
        tend=2e-7, timestep=2e-8, provided_function=lambda t: 0.0)
    d = _paths.since(before)
    assert d.get('newton_c:pf', 0) > 0 and d.get('newton_c:served', 0) == 0


def test_an_instance_shadow_of_the_newton_is_honoured():
    """A caller wrapping `_newton` on the instance expects its calls."""
    tr = Transient(_mos_chain(2), toolkit=circuit.numeric)
    calls = []
    real = tr._newton

    def wrapped(*a, **k):
        calls.append(1)
        return real(*a, **k)
    tr._newton = wrapped
    tr.solve(tend=2e-7, timestep=2e-8)
    assert len(calls) >= 10


def test_a_class_patch_and_a_subclass_take_the_python_newton(monkeypatch):
    """`_newton` patched on the class, or overridden by a subclass: today's
    Newton runs, every step, counted (the C stood in for both until
    2026-10-04)."""
    calls = []
    real = Transient._newton

    def counting(self, *a, **k):
        calls.append(1)
        return real(self, *a, **k)
    monkeypatch.setattr(Transient, '_newton', counting)
    before = _paths.snapshot()
    Transient(_mos_chain(2), toolkit=circuit.numeric).solve(tend=4e-7, timestep=2e-8,
                                                            fixed_timestep=True)
    d = _paths.since(before)
    assert len(calls) >= 15 and d.get('newton_c:served', 0) == 0, d
    assert d.get('newton_c:patched', 0) >= 15
    monkeypatch.undo()

    class Mine(Transient):
        def _newton(self, *a, **k):
            calls.append(2)
            return super()._newton(*a, **k)
    before = _paths.snapshot()
    Mine(_mos_chain(2), toolkit=circuit.numeric).solve(tend=4e-7, timestep=2e-8,
                                                       fixed_timestep=True)
    d = _paths.since(before)
    assert calls.count(2) >= 15 and d.get('newton_c:class', 0) >= 15, d


def test_the_c_solve_serves_the_multistep_step():
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    before = _paths.snapshot()
    Transient(_mos_chain(), toolkit=circuit.numeric).solve(tend=1e-6, timestep=2e-8,
                                                            fixed_timestep=True)
    d = _paths.since(before)
    assert d.get('newton_c:served', 0) >= 40
    assert d.get('newton_c:jfold', 0) == d.get('newton_c:served', 0)


## -- the sparse chord (the C loop plan's stage 2, step d) ------------------

def _sparse_chain():
    """A 250-stage MOS chain with capacitors: 254 unknowns, a sparse run --
    `StructuralKLU`, the chord."""
    return _mos_chain(250, cap=1e-13)


def _sparse_run(on, **kw):
    old = _tran_newton_c.ENABLED
    _tran_newton_c.ENABLED = on
    before = _paths.snapshot()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr = Transient(_sparse_chain(), toolkit=circuit.numeric, **kw.pop('make', {}))
            res = tr.solve(**kw)
    finally:
        _tran_newton_c.ENABLED = old
    sk = tr._get_linearsolver()._choice
    st_ = tr.statistics
    return ((np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
             st_.accepted_steps, st_.rejected_steps, st_.newton_iterations),
            (type(sk).__name__, sk.analyses, sk.factors, sk.refactors, sk.numerics, sk.foreign),
            _paths.since(before))


@pytest.mark.parametrize('integ, fixed', [(Gear2Integrator, False), (TrapezoidalIntegrator, False),
                                          (EulerIntegrator, False), (Gear2Integrator, True)])
def test_the_sparse_chord_is_the_python_one(integ, fixed):
    """KLU from C on `StructuralKLU`'s own analysis: the same bytes, steps
    and iterations, and the solver's counts the same (the numeric the C
    makes counted as the Python's `factor`)."""
    kw = {'tend': 8e-9, 'timestep': 2e-10, 'make': {'integrator': integ()}}
    if fixed:
        kw['fixed_timestep'] = True
    a = _sparse_run(False, **dict(kw, make=dict(kw['make'])))
    kw['make'] = {'integrator': integ()}
    b = _sparse_run(True, **kw)
    assert a[0] == b[0]
    assert a[1] == b[1] and a[1][0] == 'StructuralKLU'
    assert b[2].get('newton_c:sparse', 0) > 10 and a[2].get('newton_c:sparse', 0) == 0


def test_a_sparse_bail_leaves_no_trace(monkeypatch):
    """Every sparse C solve runs to its end -- its numeric made and freed
    -- and hands the step back: the Python's answer and counts."""
    drv = _tran_newton_c.driver()
    if drv is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    ffi, cfn, lapack = drv
    kw = {'tend': 4e-9, 'timestep': 2e-10}
    a = _sparse_run(False, **dict(kw))

    def bailing(s):
        cfn(s)
        return 2
    monkeypatch.setattr(_tran_newton_c, '_driver', (ffi, bailing, lapack))
    b = _sparse_run(True, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1]
    assert b[2].get('newton_c:served', 0) == 0 and b[2].get('newton_c:bail:x', 0) > 10


def test_a_core_pattern_outside_the_klu_pattern_declines(monkeypatch):
    """The core's pattern given an entry KLU's does not hold (two drains
    no element joins; the core writes its zero there): the C declines
    (`newton_c:klu`) -- in the Python a nonzero there would go to the
    solver's fallback -- and the run is the Python's."""
    from pycircuit.circuit import _tran_core
    orig = _tran_core._Core.__init__

    def init(self, cir, *a, **k):
        orig(self, cir, *a, **k)
        n = self.n
        extra = cir.get_node_index('d10') * n + cir.get_node_index('d200')
        pat = np.ascontiguousarray(np.append(self.arrays['pat'], extra), dtype=np.int64)
        self.arrays['pat'] = pat
        self.handles.append(self.ffi.from_buffer('int64_t *', pat))
        self.cs.pat = self.handles[-1]
        self.npat += 1
        self.cs.npat = self.npat
    monkeypatch.setattr(_tran_core._Core, '__init__', init)
    kw = {'tend': 4e-9, 'timestep': 2e-10}
    a = _sparse_run(False, **dict(kw))
    b = _sparse_run(True, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1]
    assert b[2].get('newton_c:klu', 0) > 0 and b[2].get('newton_c:sparse', 0) == 0


def test_the_sparse_comparison_sees_a_planted_difference(monkeypatch):
    """Alive: the first drain's diagonal read from the second's place in J
    (the C's map planted) -- another run."""
    bind = _tran_newton_c._klu_bind

    def planted(rec, sk):
        if not bind(rec, sk):
            return False
        kat, n = rec.klu_keep[1], rec.n
        cir = planted.cir
        d0, d1 = (cir.get_node_index(k) * (n + 1) for k in ('d0', 'd1'))
        kat[np.nonzero(kat == d0)[0][0]] = d1
        return True
    kw = {'tend': 4e-9, 'timestep': 2e-10}
    a = _sparse_run(True, **dict(kw))
    planted.cir = _sparse_chain()
    monkeypatch.setattr(_tran_newton_c, '_klu_bind', planted)
    b = _sparse_run(True, **dict(kw))
    assert a[0] != b[0]


## -- the full Newton's rounding floor (roadmap 197) -------------------------

def _companion_chain(cells):
    """Per cell a VCVS setting `b{k}` to 0.7 V(s{k}) -- `s{k}` a divider's
    node, fed by the cell before -- and 354 pF on `b{k}`: at a femtosecond
    step the capacitor's companion makes the VCVS's current a ~2e-4 A
    difference of ~6e5 A terms, and its update meets the rounding of its
    residual (`nrsolver._within_rounding`).  One cell: a dense run; 16:
    StructuralKLU."""
    from pycircuit.circuit.elements import VCVS
    c = SubCircuit()
    c['v'] = VSin('in', gnd, v=1.65, va=0.1, freq=1e6)
    prev = 'in'
    for k in range(cells):
        c[f'r1{k}'] = R(prev, f's{k}', r=1e3)
        c[f'r2{k}'] = R(f's{k}', gnd, r=3e3)
        c[f'e{k}'] = VCVS(f's{k}', gnd, f'b{k}', gnd, g=0.7)
        c[f'c{k}'] = C(f'b{k}', gnd, c=354e-12)
        prev = f'b{k}'
    return c


@pytest.mark.parametrize('cells', [1, 16])
def test_the_rounding_floor_is_the_python_one(cells, monkeypatch):
    """The full Newton's rounding floor in C (newton.c's `floor_holds`: the
    rows of J^-1 by numpy's dgesv of J^T, or by klu_tsolve on the step's
    numeric): twenty femtosecond steps, the same bytes, statistics and
    warnings as the Python's -- whose Newtons the floor converged."""
    from pycircuit.circuit import nrsolver
    if _tran_newton_c.driver() is None:
        pytest.skip(f'the C solve is off: {_tran_newton_c.STATUS}')
    kw = {'tend': 2e-14, 'timestep': 1e-15, 'fixed_timestep': True,
          'make': {'reltol': 1e-9, 'iabstol': 1e-13, 'vabstol': 1e-12}}
    seen, orig = [], nrsolver._within_rounding

    def within(*a):
        seen.append(orig(*a))
        return seen[-1]
    monkeypatch.setattr(nrsolver, '_within_rounding', within)
    a = _run(lambda: _companion_chain(cells), False, **dict(kw, make=dict(kw['make'])))
    monkeypatch.setattr(nrsolver, '_within_rounding', orig)
    b = _run(lambda: _companion_chain(cells), True, **kw)
    assert a[0] == b[0] and a[1] == b[1] and a[2] == b[2] and a[3] == b[3]
    assert sum(seen) >= 10, seen
    assert b[4].get('newton_c:served', 0) == 19, b[4]
    assert (b[4].get('newton_c:sparse', 0) > 0) == (cells > 1), b[4]
