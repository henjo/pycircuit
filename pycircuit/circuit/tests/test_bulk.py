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
    """A kind other than R or C, ends outside its nodes, arrays of unequal
    length, and a name an element or a bank holds."""
    c = SubCircuit()
    c['r0'] = R('a', gnd, r=1.0)
    with pytest.raises(ValueError, match='R or C'):
        c.add_bank('lb', 'L', ['a', gnd], [0], [1], [1.0])
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
.tran 10p 5n uic
.end
"""


def test_the_importer_banks_the_plain_resistors_and_capacitors(tmp_path):
    """`import_netlist(..., bulk=True)`: the resistors and the capacitors
    whose one parameter is a number held as two banks (a capacitor with an
    initial condition stays an element), the operating point and the
    transient by node name the element-built import's within rounding (the
    deck's resistors sit among its sources, a bank's after them: node 4's
    voltage one ulp apart), the same steps."""
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
    assert got[True][0] == ['c3', 'i1', 'v1']
    assert [bk.kind for bk in imp.circuit._banks] == ['R', 'C']
    assert [len(bk) for bk in imp.circuit._banks] == [4, 2]
    assert got[True][3] == got[False][3]
    for n in ('1', '2', '3', '4'):
        assert abs(got[True][1][n] - got[False][1][n]) <= 1e-14, n
    for n in ('2', '3', '4'):
        assert np.max(np.abs(got[True][2][n] - got[False][2][n])) <= 1e-12, n
