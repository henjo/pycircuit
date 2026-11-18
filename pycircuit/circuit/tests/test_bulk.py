"""THE BULK LINEAR ELEMENTS (`_bulk`, plan 261 L21 B2): resistors and
capacitors held as arrays in a bank beside a circuit's elements stamp as
the same elements would -- added after every other element, the same
circuit gives the same bytes: its matrices and vectors through the stamp
plan and the legacy loops, its Jacobian's structure, its operating point
and transient on the dense and the packed layouts, through the C paths."""
import copy
import hashlib
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _packed, _paths, _pcnr_c, _stamp_plan, numeric, pcnr
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.analysis import jacobian_structure
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import IS, VS, C, IPulse, R
from pycircuit.circuit.transient import Transient


def _grid(nx, ny, banked, extra=None):
    """An nx x ny RC grid (each node to ground through a capacitor, the
    supply's 0 F: a stamp with no entries), a supply at one corner, a
    pulsed and a constant current source: its resistors and capacitors
    elements or two banks, after every other element either way
    (`extra(c)` first), its nodes added first in one order."""
    c = SubCircuit()
    nm = lambda i, j: f'n{i}_{j}'
    names = [nm(i, j) for i in range(nx) for j in range(ny)]
    c.add_nodes(*names)
    c['vdd'] = VS('n0_0', gnd, v=1.0)
    c['ip'] = IPulse(nm(nx - 1, ny - 1), gnd, i1=0.0, i2=1e-3, td=1e-10, tr=1e-10, tf=1e-10,
                     pw=5e-10, per=2e-9)
    c['is'] = IS(nm(nx // 2, ny // 2), gnd, i=2e-4)
    if extra is not None:
        extra(c)
    rs, cs = [], []
    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                rs.append((nm(i, j), nm(i + 1, j), 1.0 + 0.1 * ((7 * i + 3 * j) % 5)))
            if j + 1 < ny:
                rs.append((nm(i, j), nm(i, j + 1), 2.0 + 0.1 * ((5 * i + j) % 4)))
            cs.append((nm(i, j), gnd, 1e-12 * (1 + (i + j) % 3) if i + j else 0.0))
    if banked:
        nodes = [*names, gnd]
        ix = {n: k for k, n in enumerate(names)}
        ix[gnd] = len(names)
        c.add_bank('rb', 'R', nodes, [ix[a] for a, _b, _v in rs], [ix[b] for _a, b, _v in rs],
                   [v for _a, _b, v in rs])
        c.add_bank('cb', 'C', nodes, [ix[a] for a, _b, _v in cs], [ix[b] for _a, b, _v in cs],
                   [v for _a, _b, v in cs])
    else:
        for k, (a, b, v) in enumerate(rs):
            c[f'r{k}'] = R(a, b, r=v)
        for k, (a, b, v) in enumerate(cs):
            c[f'c{k}'] = C(a, b, c=v)
    return c


def _h(a):
    return hashlib.sha1(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()[:12]


@pytest.mark.parametrize('size', [(3, 4), (12, 15)])
@pytest.mark.parametrize('plan', [True, False])
def test_a_bank_stamps_as_its_elements_would(size, plan, monkeypatch):
    """G, C, i and q at a random x -- through the stamp plan (its constant
    groups: a bank a run after every element) and through the legacy loops
    -- and the Jacobian's structure: the element-built circuit's bytes."""
    monkeypatch.setattr(_stamp_plan, 'ENABLED', plan)
    e, b = _grid(*size, False), _grid(*size, True)
    assert [n.name for n in e.nodes] == [n.name for n in b.nodes] and len(b.elements) == 3
    x = np.random.default_rng(1).standard_normal(e.n)
    for m in ('G', 'C', 'i', 'q'):
        assert getattr(e, m)(x).tobytes() == getattr(b, m)(x).tobytes(), m
    je, jb = jacobian_structure(e, 0), jacobian_structure(b, 0)
    assert je[0] == jb[0] and np.array_equal(je[1], jb[1]) and np.array_equal(je[2], jb[2])


def test_the_plans_fallback_takes_the_banks():
    """An element answering G and i with an object array sends the stamp
    plan to its fallback (`_legacy_matrix` / `_legacy_vector`): the banks'
    entries there too, the values the element-built circuit's."""
    def odd(c):
        r = R('n1_1', 'n2_1', r=3.0)
        G0, i0 = r.G, r.i
        r.G = lambda x, *a: np.asarray(G0(x, *a), dtype=object)
        r.i = lambda x, *a: np.asarray(i0(x, *a), dtype=object)
        c['rx'] = r
    got = []
    for banked in (False, True):
        c = _grid(4, 5, banked, odd)
        x = np.random.default_rng(4).standard_normal(c.n)
        before = _paths.snapshot()
        got.append([np.asarray(getattr(c, m)(x), float).tobytes() for m in ('G', 'i')])
        d = _paths.since(before)
        assert d.get('plan.G:legacy') == 1 and d.get('plan.i:legacy') == 1, d
    assert got[0] == got[1]


@pytest.mark.parametrize('packed', [False, True])
def test_a_banked_run_is_the_element_built_one(packed, monkeypatch):
    """The operating point and an adaptive Gear-2 run, dense and on THE
    PACKED LAYOUT (from any size): the same bytes and steps both ways, the
    C paths serving the banked circuit as they serve the element-built."""
    monkeypatch.setattr(_packed, 'PACKED_FROM', 0 if packed else 1 << 30)
    got = []
    for banked in (False, True):
        c = _grid(12, 15, banked)
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(c, toolkit=numeric).solve()
            tr = Transient(c, toolkit=numeric)
            res = tr.solve(tend=4e-9, timestep=1e-11)
        d = _paths.since(before)
        got.append(((_h(dc.x), _h(res.x), tr.statistics.accepted_steps),
                    {k: d.get(k, 0) for k in ('newton_c:served', 'lte_c:served')}))
    assert got[0] == got[1], got
    assert got[1][1]['newton_c:served'] > 10 and got[1][1]['lte_c:served'] > 10, got[1][1]
    ## (the banked circuit with every C path off: the Python's bytes the C's)
    from pycircuit.circuit import (
        _dc_newton_c,
        _tran_loop_c,
        _tran_lte_c,
        _tran_newton_c,
    )
    for mod in (_dc_newton_c, _tran_loop_c, _tran_lte_c, _tran_newton_c):
        monkeypatch.setattr(mod, 'ENABLED', False)
    c = _grid(12, 15, True)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dc = DC(c, toolkit=numeric).solve()
        tr = Transient(c, toolkit=numeric)
        res = tr.solve(tend=4e-9, timestep=1e-11)
    d = _paths.since(before)
    assert (_h(dc.x), _h(res.x), tr.statistics.accepted_steps) == got[1][0]
    assert not any(d.get(k) for k in ('dcnewton_c:served', 'newton_c:served', 'lte_c:served',
                                      'loop_c:steps')), d
    assert (d.get('core:packed', 0) > 0) == packed, d


def test_noise_refuses_a_resistor_bank():
    """`CY` (noise) where a resistor bank is held: refused (`elements.R` is
    noisy by default, a bank has no noise model); a capacitor bank alone
    adds no noise, as a capacitor does not."""
    b = _grid(3, 4, True)
    x = np.zeros(b.n)
    with pytest.raises(NotImplementedError, match='no noise model'):
        b.CY(x, 1.0)
    cap = SubCircuit()
    cap.add_nodes('a')
    cap['r'] = R('a', gnd, r=1.0)
    cap.add_bank('cb', 'C', ['a', gnd], [0], [1], [1e-12])
    alone = SubCircuit()
    alone.add_nodes('a')
    alone['r'] = R('a', gnd, r=1.0)
    alone['c'] = C('a', gnd, c=1e-12)
    assert np.array_equal(cap.CY(np.zeros(cap.n), 1.0), alone.CY(np.zeros(alone.n), 1.0))


def test_a_copy_and_the_views_carry_the_banks(monkeypatch):
    """A copy of a banked circuit, the flat view of one holding a
    subcircuit instance (its transient by the C paths), and PCNR's view of one holding a diode (its
    operating point by PCNR, its C and with the C off its Python): each the
    element-built circuit's bytes."""
    e, b = _grid(4, 5, False), _grid(4, 5, True)
    x = np.random.default_rng(2).standard_normal(e.n)
    assert copy.copy(b).G(x).tobytes() == e.G(x).tobytes()

    class Stage(SubCircuit):
        terminals = ('a', 'b')

        def __init__(self, *args, **kw):
            super().__init__(*args, **kw)
            self['R1'] = R('a', 'm', r=10.0)
            self['C1'] = C('m', 'b', c=1e-12)

    def stage(c):
        c['X1'] = Stage('n0_1', 'n1_1')
    e, b = _grid(4, 5, False, stage), _grid(4, 5, True, stage)
    x = np.random.default_rng(3).standard_normal(e.n)
    assert b._flat_view() is not b and getattr(b._flat_view(), '_banks', None)
    for m in ('G', 'C', 'i', 'q'):
        assert getattr(e, m)(x).tobytes() == getattr(b, m)(x).tobytes(), m
    ## (the transient's C takes the hierarchy through its flat view)
    got = []
    for banked in (False, True):
        c = _grid(4, 5, banked, stage)
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr = Transient(c, toolkit=numeric)
            got.append((_h(tr.solve(tend=2e-9, timestep=1e-11).x), tr.statistics.accepted_steps))
        d = _paths.since(before)
        assert d.get('flat:view', 0) > 0 and d.get('newton_c:served', 0) > 0, d
    assert got[0] == got[1]

    def diode(c):
        c['d1'] = eh.DiodeSpiceHdl('n2_2', gnd, IS=1e-14, toolkit=numeric)
    views, view = [], pcnr.pcnr_view

    def counted(circuit, *a):
        v = view(circuit, *a)
        views.append(len(getattr(v, '_banks', ())))
        return v
    monkeypatch.setattr(pcnr, 'pcnr_view', counted)
    got = []
    for banked, c_on in ((False, True), (True, True), (True, False)):
        monkeypatch.setattr(_pcnr_c, 'ENABLED', c_on)
        c, views[:] = _grid(4, 5, banked, diode), []
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            got.append(_h(DC(c, toolkit=numeric).solve().x))
        d = _paths.since(before)
        assert (d.get('pcnr_c:served', 0) > 0) == c_on, d
        assert views and set(views) == {2 if banked else 0}, views
    assert got[0] == got[1] == got[2]


def test_a_bank_refuses_what_it_cannot_be():
    """A kind other than R, C or I, ends outside its nodes, arrays of
    unequal length, a current source bank without a pulse for each source,
    a waveform on a resistor bank, and a name an element or a bank holds."""
    c = SubCircuit()
    c['r0'] = R('a', gnd, r=1.0)
    with pytest.raises(ValueError, match='R, C or I'):
        c.add_bank('lb', 'L', ['a', gnd], [0], [1], [1.0])
    with pytest.raises(ValueError, match='a pulse'):
        c.add_bank('ib', 'I', ['a', gnd], [0, 0], [1, 1], [0.0, 0.0], [[0, 1, 0, 1, 1, 1, 0]])
    with pytest.raises(ValueError, match='no waveform'):
        c.add_bank('rb', 'R', ['a', gnd], [0], [1], [1.0], [[0, 1, 0, 1, 1, 1, 0]])
    with pytest.raises(ValueError, match='outside'):
        c.add_bank('rb', 'R', ['a', gnd], [0], [2], [1.0])
    with pytest.raises(ValueError, match='one length'):
        c.add_bank('rb', 'R', ['a', gnd], [0, 1], [1], [1.0])
    with pytest.raises(ValueError, match='of that name'):
        c.add_bank('r0', 'R', ['a', gnd], [0], [1], [1.0])
    c.add_bank('rb', 'R', ['a', gnd], [0], [1], [1.0])
    with pytest.raises(ValueError, match='of that name'):
        c.add_bank('rb', 'R', ['a', gnd], [0], [1], [1.0])


_DECK = """* bulk elements
v1 1 0 1.0
r1 1 2 10
r2 2 3 20
r3 3 0 30
r4 3 4 40
c1 2 0 1p
c2 3 0 2p
c3 4 0 1p ic=0.25
i1 0 4 pulse(0 1m 1n 0.1n 0.1n 1n 4n)
i2 0 3 dc 1u ac 1
i3 0 2 pulse(0 0.5m 0.5n 0.2n 0.1n 1n 3n) ac 1
i4 2 0 pulse(0.1m 0.3m 0.7n 0.1n 0.3n 0.5n 2n)
.tran 10p 5n uic
.end
"""


def test_the_importer_banks_the_plain_resistors_and_capacitors(tmp_path):
    """`import_netlist(..., bulk=True)`: the resistors and the capacitors
    whose one parameter is a number, and the pulse current sources without
    an AC value, held as three banks (a capacitor with an initial condition,
    a DC current source and a pulse source with an AC value stay elements),
    the operating point and the transient by node name the element-built
    import's within rounding (the deck's resistors and sources sit among
    each other, a bank's after them), the same steps."""
    from pycircuit.circuit.spice_import import import_netlist
    p = tmp_path / 'bulk.cir'
    p.write_text(_DECK)
    got = {}
    for bulk in (False, True):
        imp = import_netlist(str(p), bulk=bulk)
        cir = imp.circuit
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(cir, epar=imp.epar()).solve()
            tr, kw = imp.transient()
            res = tr.solve(**kw)
        got[bulk] = (sorted(cir.elements), {n: float(dc.v(n)) for n in ('1', '2', '3', '4')},
                     {n: np.asarray(res.v(n).y, float) for n in ('2', '3', '4')}, tr.statistics.accepted_steps)
    assert got[True][0] == ['c3', 'i2', 'i3', 'v1']
    assert [bk.kind for bk in imp.circuit._banks] == ['R', 'C', 'I']
    assert [len(bk) for bk in imp.circuit._banks] == [4, 2, 2]
    assert got[True][3] == got[False][3]
    for n in ('1', '2', '3', '4'):
        assert abs(got[True][1][n] - got[False][1][n]) <= 1e-14, n
    for n in ('2', '3', '4'):
        assert np.max(np.abs(got[True][2][n] - got[False][2][n])) <= 1e-12, n


## -- THE CURRENT SOURCE BANKS (`_bulk`, kind 'I'; plan 261 L21 B3) ------------

def _sources(m=14, seed=3):
    """`m` pulse current sources' nodes and values: random levels, edges
    and periods -- no period, zero rise and fall (the floor), a zero delay,
    a delay longer than the period (its next event from the period before
    the one `t` falls in) among them -- and a DC value on some."""
    rng = np.random.default_rng(seed)
    out = []
    for k in range(m):
        a, b = f'n{int(rng.integers(1, 9))}', ('0' if k % 3 else f'n{int(rng.integers(1, 9))}')
        i1, i2 = (float(x) for x in rng.uniform(-1e-3, 1e-3, 2))
        td, tr, tf, pw = (float(x) for x in rng.uniform(0.0, 2e-9, 4))
        if k == 1:
            tr = tf = 0.0
        if k == 2:
            td = 0.0
        per = (0.0, 5e-9, 7.3e-9, 3.1e-9)[k % 4]
        if k == 5:
            ## (between one period and two: the period before `t`'s gives
            ## the next event just after each period's start)
            td, per = 2.55e-9, 1.7e-9
        elif per:
            per = max(per, td + tr + pw + tf + 1e-10)
        out.append((a, b if b != a else '0', 0.0 if k % 2 else float(rng.uniform(-1e-4, 1e-4)),
                    (i1, i2, td, tr, tf, pw, per)))
    return out


def _pulsed(banked, srcs=None, int_source=False):
    """An RC chain with a supply and a pulsed voltage source (`int_source`:
    and a constant current source of an int value, which the direct pass
    adds as `u` does, in ints), then the pulse current sources `srcs`
    (`_sources`) LAST: elements, or one current source bank -- the same
    sources in the same order either way."""
    srcs = _sources() if srcs is None else srcs
    c = SubCircuit()
    if int_source:
        c['iint'] = IS('n2', gnd, i=1)
    c['vdd'] = VS('n0', gnd, v=1.0)
    for k in range(1, 9):
        c[f'r{k}'] = R(f'n{k}', f'n{k - 1}', r=50.0 + 7 * k)
        c[f'c{k}'] = C(f'n{k}', gnd, c=2e-13 * k)
    from pycircuit.circuit.elements import VPulse
    c['vp'] = VPulse('n9', gnd, v1=0.0, v2=0.5, td=0.3e-9, tr=0.2e-9, tf=0.2e-9, pw=1e-9, per=4e-9)
    c['r9'] = R('n9', 'n4', r=500.0)
    node = lambda n: gnd if n == '0' else n
    if banked:
        names, code = [], {}
        for a, b, _i, _w in srcs:
            for n in (a, b):
                if n not in code:
                    code[n] = len(names)
                    names.append(node(n))
        c.add_bank('ib', 'I', names, [code[a] for a, _b, _i, _w in srcs],
                   [code[b] for _a, b, _i, _w in srcs], [i for _a, _b, i, _w in srcs],
                   [w for _a, _b, _i, w in srcs])
    else:
        for k, (a, b, i, w) in enumerate(srcs):
            c[f'i{k}'] = IPulse(node(a), node(b), i=i, i1=w[0], i2=w[1], td=w[2], tr=w[3],
                                tf=w[4], pw=w[5], per=w[6])
    return c


def _corner_times(srcs):
    """Times on, just before and just after every source's edges in its
    first periods (as `Pulse` sums them), and some between."""
    ts = {0.0, 1e-12, 2.5e-9, 1.37e-8, 4.1e-8}
    for _a, _b, _i, (_i1, _i2, td, tr, tf, pw, per) in srcs:
        tr, tf = max(tr, 1e-18), max(tf, 1e-18)
        edges = (td, td + tr, td + tr + pw, td + tr + pw + tf)
        for j in range(3 if per else 1):
            for e in edges + ((per,) if per else ()):
                t = j * per + e
                ts.update((t, np.nextafter(t, -1.0), np.nextafter(t, 1.0), t + 3e-11))
    return sorted(float(t) for t in ts if t >= 0.0)


@pytest.mark.parametrize('path', ['direct', 'direct-int', 'plan', 'loop'])
def test_a_source_bank_is_its_sources(path, monkeypatch):
    """u, dudt and the next breakpoint of a current source bank: the same
    sources as elements give the same bytes at every edge of every source,
    a hair either side of it and between -- through the direct source pass
    (its Python floats, and its arrays where a source's value is an int),
    the plan's source loop and the legacy loop -- and nothing outside the
    time domain."""
    if path in ('plan', 'loop'):
        monkeypatch.setattr(_stamp_plan, 'SOURCE_DIRECT', False)
    if path == 'loop':
        monkeypatch.setattr(_stamp_plan, 'SOURCE_PLAN', False)
    arrays = []
    real = _stamp_plan._SourceDirect._arrays
    monkeypatch.setattr(_stamp_plan._SourceDirect, '_arrays',
                        lambda self, *a: arrays.append(1) or real(self, *a))
    from pycircuit.circuit import _tran_sources
    monkeypatch.setattr(_tran_sources, 'ENABLED', False)
    srcs = _sources()
    e, b = (_pulsed(k, srcs, int_source=path == 'direct-int') for k in (False, True))
    assert e.n == b.n and [n.name for n in e.nodes] == [n.name for n in b.nodes]
    ts = _corner_times(srcs)
    assert len(ts) > 300
    before = _paths.snapshot()
    for t in ts:
        for an in ('tran', 'dc'):
            assert e.u(t, analysis=an).tobytes() == b.u(t, analysis=an).tobytes(), (t, an)
        assert e.dudt(t, analysis='tran').tobytes() == b.dudt(t, analysis='tran').tobytes(), t
        assert e.next_event(t) == b.next_event(t), t
    d = _paths.since(before)
    key = {'direct': 'src.direct:served', 'direct-int': 'src.direct:served', 'plan': 'src.u:served',
           'loop': 'u:called'}[path]
    assert d.get(key, 0) >= 2 * len(ts), d
    assert (len(arrays) >= 2 * len(ts)) == (path == 'direct-int'), len(arrays)
    for an in ('ac', None):
        assert np.array_equal(e.u(0.0, analysis=an), b.u(0.0, analysis=an)), an
    ## (and the breakpoints walked as a run feeds them back: every corner)
    t, seen = 0.0, 0
    while t < 2e-8:
        nxt = e.next_event(t)
        assert nxt == b.next_event(t) and nxt > t
        t, seen = nxt, seen + 1
    assert seen > 60


def test_the_source_table_holds_a_bank():
    """The C sources' table (`_tran_sources`) built from the direct pass's
    program holds the bank's sources after the elements': its `u` and its
    earliest event the Python pass's bytes; a bank alone (no element source)
    makes one too; a bank's arrays replaced, it no longer holds."""
    from pycircuit.circuit import _tran_sources
    srcs = _sources()
    for c in (_pulsed(True, srcs), _bank_alone(srcs)):
        c.u(0.0, analysis='tran')
        sd = c._stamp_plan.methods['src:u'].direct
        assert len(sd.banks) == 1
        tab = _tran_sources.table_for(sd)
        assert tab is not None and tab.ns == len(sd.ent) + len(srcs) and tab.holds()
        for t in _corner_times(srcs)[::7]:
            assert tab.u(t).tobytes() == sd.run(c, (t, None, 'tran')).tobytes(), t
            assert tab.next_event(t) == _bank_next(c, t), t
        bk = c._banks[0]
        bk.wave = bk.wave.copy()
        assert not tab.holds()


def _bank_alone(srcs):
    """The current source bank with resistors to ground and no source element."""
    c = SubCircuit()
    names = sorted({n for a, b, _i, _w in srcs for n in (a, b)} - {'0'})
    for k, n in enumerate(names):
        c[f'r{k}'] = R(n, gnd, r=100.0 + k)
    full = names + [gnd]
    ix = {n: k for k, n in enumerate(names)}
    ix['0'] = len(names)
    c.add_bank('ib', 'I', full, [ix[a] for a, _b, _i, _w in srcs], [ix[b] for _a, b, _i, _w in srcs],
               [i for _a, _b, i, _w in srcs], [w for _a, _b, _i, w in srcs])
    return c


def _bank_next(c, t):
    """The circuit's earliest event before `SubCircuit.next_event`'s
    `maximum` with `t`: its elements' and its banks'."""
    from pycircuit.circuit import _bulk
    events = [el.next_event(t) for _i, el in c._hook_elements('next_event')]
    return min(events + [_bulk.next_event(c, t)])


@pytest.mark.parametrize('method', ['gear2', 'stage', 'transform'])
def test_a_source_bank_run_is_its_sources_run(method, monkeypatch):
    """The operating point and an adaptive run (Gear-2: the multistep C
    loop; Radau: the stage loop, the transform loop) with the sources as
    elements and as a bank: the same bytes and steps, the C loop and the C
    table serving both; with the loops' C, then the table's, switched off
    the banked run's bytes the same."""
    from pycircuit.circuit import (
        _tran_loop_c,
        _tran_sources,
        _tran_stageloop_c,
        _tran_tcloop_c,
    )
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    srcs = _sources()

    def run(banked):
        c = _pulsed(banked, srcs)
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dc = DC(c, toolkit=numeric).solve()
            kw = ({} if method == 'gear2' else
                  {'integrator': RadauIIA3Integrator(), 'radau_transform': method == 'transform'})
            tr = Transient(c, toolkit=numeric, **kw)
            res = tr.solve(tend=1.2e-8, timestep=1e-10)
        d = _paths.since(before)
        st = tr.statistics
        return (_h(dc.x), _h(res.x), st.accepted_steps, st.rejected_steps), d
    (e, de), (b, db) = run(False), run(True)
    assert e == b, (e, b)
    loop = {'gear2': 'loop_c:steps', 'stage': 'stageloop_c:steps',
            'transform': 'tcloop_c:steps'}[method]
    for d in (de, db):
        assert d.get('sources_c:built', 0) >= 1 and d.get(loop, 0) > 0.8 * e[2], (loop, d)
    for mod in (_tran_loop_c, _tran_stageloop_c, _tran_tcloop_c):
        monkeypatch.setattr(mod, 'ENABLED', False)
    off, d = run(True)
    assert off == b and not d.get(loop) and d.get('src.direct:served', 0) > e[2], d
    monkeypatch.setattr(_tran_sources, 'ENABLED', False)
    off, d = run(True)
    assert off == b and not d.get('sources_c:built'), d


@pytest.mark.parametrize('i1, i2', [(1e-3, 0.4e-3), (0.4e-3, 1e-3), (-0.5e-3, 0.6e-3)])
def test_the_auto_excursion_scale_reads_a_source_bank(i1, i2):
    """`max_di_step='auto'`'s current scale (`_source_signal_scales`) the
    sources' largest swing as elements and as a bank -- the first level,
    the second and their difference each the largest once."""
    srcs = [(a, b, i, (i1, i2) + w[2:]) for a, b, i, w in _sources()[:3]]
    got = [Transient(_pulsed(b, srcs), toolkit=numeric)._source_signal_scales() for b in (False, True)]
    assert got[0] == got[1] and got[1][1] == max(abs(i1), abs(i2), abs(i2 - i1)), got


def test_jax_refuses_a_bank():
    """`JAXTransient` on a circuit holding a bank: refused, naming it."""
    pytest.importorskip('jax')
    from pycircuit.circuit.jaxtransient import JAXTransient
    with pytest.raises(NotImplementedError, match=r"bulk elements.*\['ib'\]"):
        JAXTransient(_pulsed(True))
