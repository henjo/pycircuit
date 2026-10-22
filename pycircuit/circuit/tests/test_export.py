"""The circuit-to-Python exporter (`export`): what a rebuild must keep -- the
node order, givenness, the class (a collapse variant's base, a BSourceHdl's
expression, a user class's source), keyword names, exact values -- each
shown kept, and each shown CAUGHT by the check where a copy loses it; what
it refuses, by name."""
import textwrap
import warnings

import numpy as np
import pytest

from pycircuit.circuit import circuit, elements, elements_hdl, export, func, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
from pycircuit.circuit.spice_import import import_netlist
from pycircuit.circuit.transient import Transient
from pycircuit.utilities.param import Parameter

pytestmark = pytest.mark.usefixtures('restore_default_toolkit')

#: a deck of the kinds the importer makes: a level 3 MOS whose PHI comes from
#: NSUB (given: not PHI) and whose `as` is a Python keyword, a diode whose RS
#: is 0 (a collapse variant), a pulse, a sine with an AC phase of its own
DECK = """
    exporter deck
    Vdd vdd 0 DC 3
    Vin in 0 AC 1 60 SIN(1.2 0.3 1meg 0 0 30)
    Vck ck 0 PULSE(0 3 2n 1n 1n 20n 50n)
    M1 out in 0 0 n3 w=4u l=1u as=2p
    R1 vdd out 10k
    C1 out 0 50f
    D1 0 out d1
    R2 ck out 100k
    .model n3 nmos level=3 vto=0.7 uo=600 nsub=1e16 tox=2e-8
    .model d1 d is=1e-14
    .tran 1n 200n
"""


def _deck(tmp_path):
    p = tmp_path / 'deck.cir'
    p.write_text(textwrap.dedent(DECK).lstrip('\n'))
    circuit.default_toolkit = circuit.numeric
    return import_netlist(str(p))


def _digest(res):
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes())


class Pinch(Behavioural):
    """A user class (module level: written from its source)."""
    instparams = [Parameter(name='g', desc='G', unit='S', default=1e-3)]  # noqa: RUF012

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        return Contribution(b.I, g * b.V ** 3)                  # noqa: F821


def test_an_imported_deck_round_trips_to_the_bit(tmp_path):
    """The rebuilt deck is the original as the solvers see it: the selfcheck,
    the DC answer and the transient, bit for bit -- its source holding the
    level 3 MOS's `as` through `**{}`, the diode as its base class, the sine
    source's two phases."""
    imp = _deck(tmp_path)
    c = imp.circuit
    src = export.to_python(c, 'deck')
    assert "**{'as': 2e-12}" in src and 'elements_hdl.DiodeSpiceHdl(' in src
    assert 'acphase=60.0' in src and 'phase=30.0' in src
    assert type(c['d1']).__name__.startswith('DiodeSpiceHdl_collapse')
    assert export.selfcheck(c, src, name='deck') is None
    new = export.rebuild(src, 'deck')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        xa, xb = (np.asarray(DC(k, epar=imp.epar()).solve().x) for k in (c, new))
        ra, rb = (Transient(k, epar=imp.epar()).solve(tend=60e-9, timestep=1e-9) for k in (c, new))
    assert xa.tobytes() == xb.tobytes()
    assert _digest(ra) == _digest(rb)


def test_a_copy_with_every_parameter_given_is_caught(tmp_path):
    """Givenness: the level 3 MOS computes PHI from NSUB only where PHI is not
    given -- a copy writing every parameter (at its value) is another
    device, and the selfcheck says where."""
    c = _deck(tmp_path).circuit
    every = export._to_python(c, 'deck', True, given_only=False)
    assert export.selfcheck(c, every, name='deck') is not None
    assert export.selfcheck(c, export._to_python(c, 'deck', True), name='deck') is None


def test_a_copy_with_its_nodes_in_first_use_order_is_caught():
    """The node order: gnd not last and a node added before its first use --
    a copy that lets the elements add them is another unknowns' order."""
    c = SubCircuit()
    c.add_node('b')
    c['V1'] = elements.VS('a', gnd, v=1.0)
    c['R1'] = elements.R('a', 'b', r=1e3)
    c['R2'] = elements.R('b', gnd, r=2e3)
    loose = export._to_python(c, 'c', True, node_order=False)
    assert export.selfcheck(c, loose, name='c') == 'node order'
    assert export.selfcheck(c, export.to_python(c, 'c'), name='c') is None


def test_a_bsource_and_a_user_class_round_trip():
    """A BSourceHdl from its traced expression (srepr: exact, the same
    class) and a class outside the library from its source."""
    c = SubCircuit()
    c['V1'] = elements.VS('a', gnd, v=0.7)
    c['B1'] = elements_hdl.BSourceHdl('a', gnd, 'a', 'b',
                                      i_func=lambda v: 1e-3 * v ** 2 + 1e-6 * v,
                                      q_func=lambda v: 1e-15 * v)
    c['P1'] = Pinch('b', gnd, g=2e-3)
    src = export.to_python(c, 'c')
    assert 'BSourceHdl(' in src and 'class Pinch(Behavioural)' in src
    new = export.rebuild(src, 'c')
    assert type(new['B1']) is type(c['B1'])
    assert export.selfcheck(c, src, name='c') is None


def test_what_it_cannot_write_it_refuses_by_name():
    """A callable parameter, a subcircuit instance, a class local to a
    function: ExportError naming the element."""
    c = SubCircuit()
    c['B'] = elements.BSource('a', gnd, 'a', gnd, i_func=lambda v: v)
    with pytest.raises(export.ExportError, match="element 'B'"):
        export.to_python(c)

    class Local(Pinch):
        pass
    c = SubCircuit()
    c['L'] = Local('a', gnd)
    with pytest.raises(export.ExportError, match="element 'L'.*local to a function"):
        export.to_python(c)
    inner = SubCircuit()
    inner['R'] = elements.R('a', gnd)
    c = SubCircuit()
    c['X'] = inner
    with pytest.raises(export.ExportError, match="element 'X'.*hierarchy"):
        export.to_python(c)


def test_state_the_parameters_do_not_hold_is_caught():
    """A source whose time function was replaced after construction: its
    parameters rebuild another source, and the export says so."""
    c = SubCircuit()
    c['V1'] = elements.VS('a', gnd, v=1.0)
    c['R1'] = elements.R('a', gnd, r=1e3)
    c['V1'].function = func.Sin(0.0, 1.0, 1e3, toolkit=circuit.numeric)
    with pytest.raises(export.ExportError, match="element 'V1'.*(u\\(|next_event)"):
        export.to_python(c)


def test_a_limiter_s_state_is_compared_from_none_and_put_back():
    """elements.Diode keeps its limiter's last voltage (`_vlim`) -- the
    solve's state, not the circuit's: both sides compared from none, the
    original's put back."""
    c = SubCircuit()
    c['V1'] = elements.VS('a', gnd, v=0.8)
    c['R1'] = elements.R('a', 'b', r=100.0)
    c['D1'] = elements.Diode('b', gnd)
    DC(c).solve()
    before = dict(vars(c['D1']))
    src = export.to_python(c, 'c')
    assert export.selfcheck(c, src, name='c') is None
    after = vars(c['D1'])
    assert {k: repr(v) for k, v in after.items()} == {k: repr(v) for k, v in before.items()}


def test_a_sine_source_s_two_phases_survive():
    """`phase` (the sine's) and `acphase` (the AC phasor's) are written and
    compared apart: the check's AC source vector sees the AC phase."""
    c = SubCircuit()
    c['V1'] = elements.VSin('a', gnd, vac=1.0, va=1.0, freq=1e3, phase=30.0, acphase=60.0)
    c['I1'] = elements.ISin('a', gnd, iac=1e-3, ia=1e-3, freq=1e3, phase=-45.0)
    c['R1'] = elements.R('a', gnd, r=1e3)
    src = export.to_python(c, 'c')
    new = export.rebuild(src, 'c')
    for k in ('V1', 'I1'):
        assert (new[k].iparv.phase, new[k].iparv.acphase) == (c[k].iparv.phase,
                                                              c[k].iparv.acphase)
    assert np.array_equal(new.u(0.0, analysis='ac'), c.u(0.0, analysis='ac'))
    shifted = src.replace('acphase=60.0', 'acphase=60.000000000000014')
    assert export.selfcheck(c, shifted, name='c').startswith('u(ac)')


def test_a_parameter_a_hair_off_is_caught():
    """The selfcheck is alive on values: one resistance an ulp off reads
    as a different i or G at the first x."""
    c = SubCircuit()
    c['V1'] = elements.VS('a', gnd, v=1.0)
    c['R1'] = elements.R('a', gnd, r=2e3)
    src = export.to_python(c, 'c')
    assert 'r=2000.0' in src
    assert export.selfcheck(c, src.replace('r=2000.0', 'r=2000.0000000000002'),
                            name='c') in ('i at x #0', 'G at x #0')


def test_a_start_vector_round_trips_exactly():
    """`x0_to_python`: every entry by repr -- a NaN, an infinity, a signed
    zero, a subnormal kept."""
    x = np.array([1.0 / 3.0, -0.0, np.nan, np.inf, -np.inf, 5e-324, 1.65518481e10])
    ns = {'np': np}
    exec(export.x0_to_python(x, 'x0'), ns)  # noqa: S102
    y = ns['x0']
    assert y.dtype == np.float64 and y.tobytes() == x.tobytes()
    with pytest.raises(export.ExportError):
        export.x0_to_python(np.zeros((2, 2)))


def test_only_a_default_toolkit_circuit_is_exported():
    """The rebuild takes the default toolkit: a circuit on another one is
    refused."""
    c = SubCircuit(toolkit=circuit.symbolic)
    c['R1'] = elements.R('a', gnd, r=1e3, toolkit=circuit.symbolic)
    with pytest.raises(export.ExportError, match='toolkit'):
        export.to_python(c)
