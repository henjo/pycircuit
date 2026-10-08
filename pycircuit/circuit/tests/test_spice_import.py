"""`spice_import.import_netlist`: a SPICE netlist mapped onto pycircuit as
SPICE defines each element (the SPICE benchmark plan's stage 2).  An
imported deck is the circuit built by hand -- the same DC and transient
bytes; each mapping rule pinned; subcircuits flattened with `:` names;
what cannot be mapped refused, every occurrence with its file and line.
Synthetic netlists only, but for the census over the fetched benchmark
decks, which skips without them."""
import os
import textwrap
import warnings

import numpy as np
import pytest

from pycircuit._testing import benchdata
from pycircuit.circuit import circuit, elements, elements_hdl, integrator
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.spice_import import SpiceImportError, import_netlist
from pycircuit.circuit.transient import Transient


def _write(tmp_path, text, name='a.cir'):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text).lstrip('\n'))
    return str(p)


def _import(tmp_path, text, **kw):
    circuit.default_toolkit = circuit.numeric
    return import_netlist(_write(tmp_path, text), **kw)


DECK = """
    a deck of every kind
    V1 in 0 DC 1 AC 1
    Vp clk 0 PULSE(0 5 1n 1n 1n 10n 20n)
    R1 in a 1k
    C1 a 0 1p
    L1 a b 1u
    L2 b 0 2u
    K1 L1 L2 0.5
    R2 b out {2*rv}
    .param rv = 500
    G1 out 0 in 0 1m
    D1 out 0 dmod area=2
    .model dmod d is=1e-14 rs=10 cjo=1p
    Q1 c clk 0 qmod
    Rc vcc c 2k
    Vcc vcc 0 5
    .model qmod npn bf=100 is=1e-16 vaf=50
    X1 c y inv
    Rl y 0 10k
    .subckt inv in out
    M1 out in 0 0 nch w=2u l=1u
    M2 out in vdd vdd pch w=4u l=1u
    Vdd vdd 0 5
    .model nch nmos level=1 vto=0.7 kp=1e-4
    .model pch pmos level=3 vto=-0.8 uo=200 tox=2e-8
    .ends
    .tran 1n 40n
    .end
    """


def _by_hand():
    """`DECK` built by hand, in the importer's order (the coupled pair
    last, where the coupling replaced its inductors)."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['v1'] = elements.VS('in', gnd, v=1.0, vac=1.0)
    c['vp'] = elements.VPulse('clk', gnd, v=0.0, vac=0.0, v1=0.0, v2=5.0, td=1e-9, tr=1e-9,
                              tf=1e-9, pw=10e-9, per=20e-9)
    c['r1'] = elements.R('in', 'a', r=1e3)
    c['c1'] = elements.C('a', gnd, c=1e-12)
    c['r2'] = elements.R('b', 'out', r=1000.0)
    c['g1'] = elements.VCCS('in', gnd, 'out', gnd, gm=1e-3)
    c['d1'] = elements_hdl.DiodeSpiceHdl('out', gnd, IS=1e-14, rs=10.0, cjo=1e-12, area=2.0)
    c['q1'] = elements_hdl.GummelPoonNpnHdl('c', 'clk', gnd, bf=100.0, IS=1e-16, vaf=50.0)
    c['rc'] = elements.R('vcc', 'c', r=2e3)
    c['vcc'] = elements.VS('vcc', gnd, v=5.0, vac=0.0)
    c['x1:m1'] = elements_hdl.MosLevel1Hdl('y', 'c', gnd, gnd, vto=0.7, kp=1e-4, phi=0.6,
                                           w=2e-6, l=1e-6)
    c['x1:m2'] = elements_hdl.MosLevel3PmosGateChargeHdl('y', 'c', 'x1:vdd', 'x1:vdd', vto=0.8,
                                               u0=200.0, tox=2e-8, phi=0.6, w=4e-6, l=1e-6)
    c['x1:vdd'] = elements.VS('x1:vdd', gnd, v=5.0, vac=0.0)
    c['rl'] = elements.R('y', gnd, r=10e3)
    c['k1'] = elements.CoupledInductors('a', 'b', 'b', gnd, L1=1e-6, L2=2e-6, K=0.5)
    return c


def test_an_imported_deck_is_the_circuit_built_by_hand(tmp_path):
    """Every kind of element the importer maps, a subcircuit, a scoped
    model, an expression: the operating point and 20 fixed transient
    steps the same bytes as the circuit built by hand."""
    imp = _import(tmp_path, DECK)
    assert not [r for r in imp.report if 'not mapped' not in r], imp.report
    got, ref = imp.circuit, _by_hand()
    assert list(got.elements) == list(ref.elements)
    assert [str(n) for n in got.nodes] == [str(n) for n in ref.nodes]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        xa, xb = DC(got).solve().x, DC(ref).solve().x
        assert np.array_equal(np.asarray(xa), np.asarray(xb))
        ta, kw = imp.transient()
        assert kw == {'tend': 40e-9, 'timestep': 1e-9}
        ra = ta.solve(tend=20e-9, timestep=1e-9, fixed_timestep=True)
        rb = Transient(ref, epar=imp.epar()).solve(tend=20e-9, timestep=1e-9,
                                                   fixed_timestep=True)
    assert np.array_equal(np.asarray(ra.x), np.asarray(rb.x))
    assert np.asarray(ra.x).shape[1] >= 21


#: `.tran 1n 100n`'s TSTOP as SPICE reads `100n`: the digits times the scale
#: factor (not the literal 100e-9, another double).
TSTOP = 100 * 1e-9


def test_the_ngspice_deck_written_back_reads_as_the_same_circuit(tmp_path):
    """`write_ngspice` (the stage 4 references' decks): every element's
    class, parameters (to the bit -- `repr` round-trips a float) and nodes
    (under the written names); the coupling, the subcircuit's internal
    nodes, a PMOS threshold's sign, `.tran` and its print."""
    imp = _import(tmp_path, DECK)
    path = str(tmp_path / 'w.cir')
    names = imp.write_ngspice(path, probes=['out', 'x1:vdd', ('i', 'vcc'), ('i', 'x1:vdd')])
    assert names['x1:vdd'] == 'x1_vdd' and names['0'] == '0'
    ## (a source's current: printed as i(name), ngspice's column name#branch)
    assert names[('i', 'vcc')] == 'vcc#branch' and names[('i', 'x1:vdd')] == 'v_x1_vdd#branch'
    with open(path) as fh:
        text = fh.read()
    assert ('.print tran\n+ v(out) v(x1_vdd) i(vcc) i(v_x1_vdd)\n' in text
            and text.endswith('.end\n'))
    again = import_netlist(path, dialect='ngspice')
    assert len(again.elements) == len(imp.elements) == 15
    for a, b in zip(imp.elements, again.elements, strict=True):
        assert (a.cls, [names[n] for n in a.nodes], a.params) == (b.cls, b.nodes, b.params)
    assert [again.tran[k] for k in ('tstep', 'tstop', 'start')] == [1e-9, 40 * 1e-9, 'op']


def _one(tmp_path, line, extra='', tran='.tran 1n 100n'):
    imp = _import(tmp_path, f'title\n{line}\n{extra}\n{tran}\n')
    (m,) = [m for m in imp.elements if m.name == line.split()[0].lower()]
    return m


@pytest.mark.parametrize('line, cls, nodes, params', [
    ('R1 a b 2k', elements.R, ['a', 'b'], {'r': 2e3}),
    ('R1 a b r=2k', elements.R, ['a', 'b'], {'r': 2e3}),
    ('R1 a b 0', elements.VS, ['a', 'b'], {'v': 0.0, 'vac': 0.0}),
    ('C1 a 0 10pF', elements.C, ['a', '0'], {'c': 10 * 1e-12}),
    ('L1 a 0 1uH', elements.L, ['a', '0'], {'L': 1e-6}),
    ('G1 p n cp cn 2m', elements.VCCS, ['cp', 'cn', 'p', 'n'], {'gm': 2e-3}),
    ('V1 a 0 3', elements.VS, ['a', '0'], {'v': 3.0, 'vac': 0.0}),
    ('V1 a 0 DC 3 AC', elements.VS, ['a', '0'], {'v': 3.0, 'vac': 1.0}),
    ('V1 a 0 AC 2 45', elements.VS, ['a', '0'], {'v': 0.0, 'vac': 2.0, 'phase': 45.0}),
    ('I1 a 0 1mA', elements.IS, ['a', '0'], {'i': 1e-3, 'iac': 0.0}),
    ('V1 a 0 DC 9 PULSE(0 1)', elements.VPulse, ['a', '0'],
     {'v': 0.0, 'vac': 0.0, 'v1': 0.0, 'v2': 1.0, 'td': 0.0, 'tr': 1e-9, 'tf': 1e-9,
      'pw': TSTOP, 'per': TSTOP}),
    ('V1 a 0 PULSE 0 1 2n 0 0 5n 0', elements.VPulse, ['a', '0'],
     {'v': 0.0, 'vac': 0.0, 'v1': 0.0, 'v2': 1.0, 'td': 2e-9, 'tr': 1e-9, 'tf': 1e-9,
      'pw': 5e-9, 'per': TSTOP}),
    ('V1 a 0 sin 0 1', elements.VSin, ['a', '0'],
     {'v': 0.0, 'vac': 0.0, 'vo': 0.0, 'va': 1.0, 'freq': 1 / TSTOP, 'td': 0.0,
      'theta': 0.0, 'phase': 0.0}),
    ('I1 a 0 SIN(1m 2m 1meg 1n 0 30)', elements.ISin, ['a', '0'],
     {'i': 0.0, 'iac': 0.0, 'io': 1e-3, 'ia': 2e-3, 'freq': 1e6, 'td': 1e-9, 'theta': 0.0,
      'phase': 30.0}),
    ('V1 a 0 EXP(0 1 2n)', elements.VExp, ['a', '0'],
     {'v': 0.0, 'vac': 0.0, 'v1': 0.0, 'v2': 1.0, 'td1': 2e-9, 'tau1': 1e-9,
      'td2': 2e-9 + 1e-9, 'tau2': 1e-9}),
    ('V1 a 0 PWL(0 0, 1n 1, 2n 0)', elements.VPWL, ['a', '0'],
     {'v': 0.0, 'vac': 0.0, 'tvpairs': [0.0, 0.0, 1e-9, 1.0, 2e-9, 0.0]}),
])
def test_each_two_terminal_maps_as_spice_defines_it(tmp_path, line, cls, nodes, params):
    """Values as SPICE reads them; `vac` 0 unless AC; SPICE's waveform
    defaults from `.tran 1n 100n` (TR, TF: TSTEP where absent or 0; PW,
    PER: TSTOP; SIN's FREQ 1/TSTOP; EXP's TAUs TSTEP, TD2 TD1 + TSTEP);
    the DC value given beside a waveform is not the transient's; G's
    control pins first."""
    m = _one(tmp_path, line)
    assert (m.cls, m.nodes, m.params) == (cls, nodes, params)


@pytest.mark.parametrize('model, inst, cls, params', [
    ('.model n1 nmos vto=0.7 kp=1e-4 lambda=0.02', '', elements_hdl.MosLevel1Hdl,
     {'vto': 0.7, 'kp': 1e-4, 'lambd': 0.02, 'phi': 0.6}),
    ('.model n1 pmos level=1 vto=-0.7 uo=300 tox=1e-8', 'w=2u l=1u as=1p',
     elements_hdl.MosLevel1PmosGateChargeHdl,
     {'vto': 0.7, 'tox': 1e-8, 'kp': 300.0 * 1e-4 * (3.9 * 8.854214871e-12) / 1e-8,
      'phi': 0.6, 'w': 2e-6, 'l': 1e-6, 'asrc': 1e-12}),
    ('.model n1 nmos level=3 vto=0.9 uo=600 nsub=1e16 tpg=1 nss=0', 'as=2p',
     elements_hdl.MosLevel3GateChargeHdl, {'vto': 0.9, 'u0': 600.0, 'nsub': 1e16, 'as': 2e-12}),
    ('.model n1 nmos (level=3 vto=0.9 phi=0.7 l=3u)', '', elements_hdl.MosLevel3GateChargeHdl,
     {'vto': 0.9, 'phi': 0.7, 'l': 3e-6}),
    ('.model n1 nmos level=1 vto=0.7 kp=1e-4 tox=0', '', elements_hdl.MosLevel1Hdl,
     {'vto': 0.7, 'kp': 1e-4, 'tox': 0.0, 'phi': 0.6}),
])
def test_mosfets_map_as_spice_reads_their_cards(tmp_path, model, inst, cls, params):
    """LEVEL picks the class -- with the gate charge where SPICE adds
    Meyer's (level 3; level 1 with a TOX that is not zero); a PMOS threshold its magnitude; level 1's KP
    from UO and TOX where KP is not given (mos1temp.c); PHI 0.6 where
    neither PHI nor NSUB is given; TPG and NSS without effect where VTO is
    given; a model's L and W the instance's default; `as` level 1's
    `asrc`."""
    m = _one(tmp_path, f'M1 d g s b n1 {inst}', model)
    assert (m.cls, m.nodes, m.params) == (cls, ['d', 'g', 's', 'b'], params)


def test_diodes_and_bipolars(tmp_path):
    m = _one(tmp_path, 'D1 a k dm 3', '.model dm d (is=2e-14 cj0=1p pb=0.7 mj=0.4 bv=5)')
    assert (m.cls, m.nodes) == (elements_hdl.DiodeSpiceHdl, ['a', 'k'])
    assert m.params == {'IS': 2e-14, 'cjo': 1e-12, 'vj': 0.7, 'm': 0.4, 'bv': 5.0, 'area': 3.0}
    q = _one(tmp_path, 'Q1 c b e qp area=2', '.model qp pnp (bf=50 va=40 ik=0.1 cjs=0)')
    assert (q.cls, q.nodes) == (elements_hdl.GummelPoonPnpHdl, ['c', 'b', 'e'])
    assert q.params == {'bf': 50.0, 'vaf': 40.0, 'ikf': 0.1, 'area': 2.0}
    q4 = _one(tmp_path, 'Q1 c b e s qn', '.model qn npn bf=80')
    assert (q4.cls, q4.nodes, q4.params) == (elements_hdl.GummelPoonNpnHdl, ['c', 'b', 'e'],
                                             {'bf': 80.0})


def test_flattening_names_and_nodes(tmp_path):
    """Instance paths joined by `:` (Xyce's), internal nodes likewise,
    ports mapped, `0` global; a subcircuit used before its definition and
    one defined inside another; instance parameters reaching the
    expressions and the models inside."""
    imp = _import(tmp_path, """
        title
        X1 in out buf g=3
        .subckt buf a y g=1
        X2 a mid stage
        R1 mid y {g*1k}
        .subckt stage p q
        R1 p q 1
        C1 q 0 1p
        .ends stage
        .ends buf
        """)
    got = [(m.name, m.nodes, m.params) for m in imp.elements]
    assert got == [('x1:x2:r1', ['in', 'x1:mid'], {'r': 1.0}),
                   ('x1:x2:c1', ['x1:mid', '0'], {'c': 1e-12}),
                   ('x1:r1', ['x1:mid', 'out'], {'r': 3000.0})]
    assert 'x1:x2:c1' in imp.circuit.elements and imp.circuit['x1:r1'].ipar.r == 3000.0


def test_initial_conditions_only_under_uic(tmp_path):
    imp = _import(tmp_path, """
        title
        C1 a 0 1p IC=2
        L1 a 0 1u IC=1m
        R1 a 0 1k
        .ic v(a)=2
        .tran 1n 10n UIC
        """)
    by = {m.name: m.params for m in imp.elements}
    assert by['c1'] == {'c': 1e-12, 'ic': 2.0} and by['l1'] == {'L': 1e-6, 'ic': 1e-3}
    tr, _kw = imp.transient()
    assert tr.par.uic and dict(tr.par.ic) == {'a': 2.0}
    imp = _import(tmp_path, 'title\nC1 a 0 1p IC=2\nR1 a 0 1k\n.tran 1n 10n\n')
    assert {m.name: m.params for m in imp.elements}['c1'] == {'c': 1e-12}
    assert any('IC= ignored without UIC' in r for r in imp.report)


def test_an_ic_without_uic_holds_the_operating_point(tmp_path):
    """`.ic` without UIC: the transient's operating point with those nodes
    held, released at t = 0 (the plan's stage 3); `.nodeset` the point it
    is solved from."""
    imp = _import(tmp_path, """
        title
        R1 a 0 1k
        C1 a 0 1n
        R2 b 0 1k
        .ic v(a)=1
        .nodeset v(b)=0.5
        .tran 10n 1u
        """)
    tr, kw = imp.transient()
    assert not tr.par.uic and dict(tr.par.ic) == {'a': 1.0}
    assert dict(tr.par.nodeset) == {'b': 0.5}
    res = tr.solve(**kw)
    assert res.v('a', gnd).y[0] == 1.0


def test_merge_shorts_makes_one_node_of_what_a_0v_source_joins(tmp_path):
    """A power grid's pads (ibmpg1t's `vb9 _Y_n2 0 0`): the 0 V sources no
    `.print` reads left out, the nodes they join one -- ground the node a
    class becomes where ground is in it; a read one, and a 1.8 V one, kept;
    the same answer at the nodes left; a merge that would short a source
    left in refused."""
    text = """
        title
        V1 vdd 0 1.8
        R1 vdd a 1k
        Vs a b 0
        R2 b y 1k
        Vg y 0 0
        Vm b c 0
        R3 c 0 2k
        .print tran i(vm)
        .tran 1n 10n
        """
    plain = _import(tmp_path, text)
    imp = import_netlist(_write(tmp_path, text, 'b.cir'), merge_shorts=True)
    assert [m.name for m in imp.elements] == ['v1', 'r1', 'r2', 'vm', 'r3']
    assert imp.merged == {'b': 'a', 'y': '0'}
    assert {m.name: m.nodes for m in imp.elements}['r2'] == ['a', '0']
    assert imp.circuit.n == plain.circuit.n - 4
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        va = DC(plain.circuit).solve().v('a')
        assert DC(imp.circuit).solve().v('a') == pytest.approx(va, rel=1e-12, abs=1e-15)
    assert any('merge_shorts: 2 0 V sources left out, 2 nodes merged' in r for r in imp.report)
    with pytest.raises(SpiceImportError, match='merge_shorts would short this source'):
        import_netlist(_write(tmp_path, 'title\nV1 a 0 1\nVz a 0 0\nR1 a 0 1\n', 'c.cir'),
                       merge_shorts=True)


def test_an_ic_a_source_holds_is_left_out_as_spice_does(tmp_path):
    """SPICE lets the source win (CircuitSim90's gm17 holds such a node):
    the `.ic` left out and said, the others kept."""
    imp = _import(tmp_path, """
        title
        V1 s 0 5
        R1 s a 1k
        C1 a 0 1p
        .ic v(s)=1 v(a)=2
        .tran 1n 10n
        """)
    assert imp.ic == {'a': 2.0}
    assert any('.ic v(s): a voltage source or an inductor holds the node' in r
               for r in imp.report)


def test_the_transient_follows_the_netlist(tmp_path):
    """TMAX, NOOP (zeros), the options' method and reltol, the
    temperature; an override replaces any."""
    imp = _import(tmp_path, """
        title
        R1 a 0 1k
        V1 a 0 1
        .tran 1n 100n 0 2n NOOP
        .options timeint method=trap reltol=1e-3 abstol=1e-9
        .options device temp=50
        """)
    tr, kw = imp.transient()
    assert kw == {'tend': TSTOP, 'timestep': 1e-9} and imp.tran['start'] == 'noop'
    assert tr.par.uic and tr.par.timestep_max == 2e-9 and tr.par.reltol == 1e-3
    assert isinstance(tr.par.integrator, integrator.TrapezoidalIntegrator)
    assert tr.par.epar.T == 273.15 + 50
    assert any('abstol=1e-9: not mapped' in r for r in imp.report)
    assert imp.transient(reltol=1e-6)[0].par.reltol == 1e-6
    ng = import_netlist(_write(tmp_path, 'title\nR1 a 0 1\n.options temp=75 reltol=1e-2\n'
                                         '.tran 1 2\n', 'b.cir'), dialect='ngspice')
    assert ng.temp == 75.0 and ng.options == {'reltol': 1e-2}


def test_every_refusal_is_listed_with_its_line(tmp_path):
    """One error naming every element and statement that cannot be
    mapped; `strict=False` reports them and leaves those elements out."""
    text = """
        title
        E1 a 0 b 0 2
        M1 d g 0 0 n2
        .model n2 nmos level=2 vto=1
        Q1 c b 0 qs
        .model qs npn cjs=1p
        X1 a b nosuch
        R1 a b 1k tc1=0.01
        V1 a 0 SFFM(0 1 1k 5 1)
        .ic v(a)=1
        .tran 1n 10n
        """
    with pytest.raises(SpiceImportError) as e:
        _import(tmp_path, text)
    msg = str(e.value)
    for line, what in ((2, 'e1: a voltage-controlled voltage source is not supported'),
                       (3, 'MOS LEVEL 2 is not supported'),
                       (5, 'a substrate junction (CJS) is not supported'),
                       (7, 'x1: no subcircuit nosuch'),
                       (8, "resistor parameters ['tc1'] are not supported"),
                       (9, 'a SFFM waveform is not supported')):
        assert f'a.cir:{line}: ' in msg and what in msg, (line, what, msg)
    assert '.ic' not in msg
    imp = _import(tmp_path, text, strict=False)
    assert [m.name for m in imp.elements] == []
    assert sum('not supported' in r for r in imp.report) >= 5


def test_a_second_element_of_one_name_is_renamed_and_reported(tmp_path):
    imp = _import(tmp_path, 'title\nR1 a 0 1k\nR1 a 0 2k\nV1 a 0 1\n')
    assert [(m.name, m.params) for m in imp.elements] == [
        ('r1', {'r': 1e3}), ('r1#1', {'r': 2e3}), ('v1', {'v': 1.0, 'vac': 0.0})]
    assert any('a second element of that name, renamed r1#1' in r for r in imp.report)
    assert imp.circuit['r1#1'].ipar.r == 2e3


@pytest.mark.parametrize('line, says', [
    ('M1 d g s b nm', 'VTO from NSUB is not supported'),
    ('M1 d g s b ok m=2', 'instance parameters'),
    ('R1 a', 'two nodes needed'),
    ('C1 a b 1p 2p', 'one value expected'),
    ('K1 L1 L9 0.5', 'no inductor'),
    ('V1 a 0 PULSE(0)', 'PULSE takes 2 to 7 values'),
    ('R1 a.b 0 1', 'hierarchy separator'),
    ('D1 a 0 dd', 'cjo given twice'),
], ids=['vto-nsub', 'm', 'nodes', 'values', 'coupling', 'pulse', 'dot', 'alias'])
def test_refusals_say_why(tmp_path, line, says):
    text = (f'title\n{line}\nL1 a 0 1u\n.model nm nmos level=1 nsub=1e16\n'
            '.model ok nmos level=1 vto=1\n.model dd d cjo=1p cj0=2p\n.tran 1n 10n\n')
    with pytest.raises(SpiceImportError, match=says):
        _import(tmp_path, text)


#: The fetched decks' known gaps: what each import refuses today, by the
#: plan's stages (7: the bipolar substrate junction; 9: MOS level 2).  A
#: deck absent here imports.
CENSUS_GAPS = {
    'latch.cir': ('CJS',), 'opampal.cir': ('CJS',), 'gilbert_cell_hb.cir': ('CJS',),
}
_LEVEL2 = ('ab_ac', 'ab_integ', 'ab_opamp', 'cram', 'e1480', 'g1310', 'gm6', 'hussamp',
           'mosrect', 'mux8', 'nand', 'pump', 'ring', 'schmitfast', 'schmitslow')
CENSUS_GAPS.update({f'{n}.cir': ('MOS LEVEL 2',) for n in _LEVEL2})
#: The large decks import in seconds to tens of seconds: their census is a
#: benchmark's (`benchmarks/large_circuits.py`), not the suite's.
_LARGE = 'Netlists/CircuitSim90/MOS2_LARGE/'


def _census_decks():
    return sorted(f['path'] for f in benchdata.manifest()['files']
                  if f['source'] == 'xyce_regression' and f['path'].startswith('Netlists/')
                  and f['path'].endswith(('.cir', '.cir_NORUN'))
                  and not f['path'].startswith(_LARGE))


@pytest.mark.parametrize('path', _census_decks())
def test_each_fetched_deck_imports_or_names_exactly_its_known_gaps(path):
    p = benchdata.spice_data(path)
    if p is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    circuit.default_toolkit = circuit.numeric
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        imp = import_netlist(p, strict=False)
    refused = [r for r in imp.report
               if not any(s in r for s in ('not mapped', 'not read', 'ignored', 'a 0 V source',
                                           'left unconnected', 'starts at 0',
                                           'ignored under UIC', 'no effect'))]
    gaps = CENSUS_GAPS.get(os.path.basename(path), ())
    for r in refused:
        assert any(g in r for g in gaps), (path, r)
    for g in gaps:
        assert any(g in r for r in refused), (path, g)
