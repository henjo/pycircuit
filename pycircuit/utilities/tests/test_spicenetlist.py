"""`spicenetlist.read`: a SPICE netlist as SPICE wrote it -- element lines
with their files and lines, subcircuits and their scoped models and
parameters, analyses, options, prints, `.ic`, `.nodeset`, `*COMP` -- on
synthetic netlists (no benchmark data needed).  The SPICE benchmark
plan's stage 2."""
import os
import textwrap

import pytest

from pycircuit.utilities import spicenetlist
from pycircuit.utilities.spicenetlist import SpiceNetlistError, Where


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text).lstrip('\n'))
    return str(p)


def _read(tmp_path, text, name='a.cir'):
    return spicenetlist.read(_write(tmp_path, name, text))


def test_the_first_line_is_the_title_whatever_it_looks_like(tmp_path):
    net = _read(tmp_path, """
        R1 1 0 1k
        R2 1 0 2k
        """)
    assert net.title == 'R1 1 0 1k' and [c.name for c in net.cards] == ['r2']
    assert _read(tmp_path, '* a comment as the title\nr1 1 0 1\n').title == '* a comment as the title'


def test_an_element_line_its_words_and_parameters(tmp_path):
    """Parentheses and commas separate words; `name=value` pairs (spaces
    around `=` or not) are parameters; an expression is one word; names
    and nodes in lower case, expressions as written."""
    net = _read(tmp_path, """
        title
        V1 IN 0 DC 5 AC 1 PULSE(0 1 1n 1n 1n 5n 10n)
        vin 3 6 pulse (0v 10v 10ns)
        I1 1 0 PWL(0 0, 1n 1)
        m$615 30 57 57 57 pn2 w=37u l=2u ad=185p as=185p
        R2 1 2 {RV*2}
        X1 a b sub w = 2u L='3*u0'
        """)
    by = {c.name: c for c in net.cards}
    assert by['v1'].words == ['in', '0', 'dc', '5', 'ac', '1', 'pulse', '0', '1', '1n', '1n',
                              '1n', '5n', '10n'] and by['v1'].params == {}
    assert by['vin'].words == ['3', '6', 'pulse', '0v', '10v', '10ns']
    assert by['i1'].words == ['1', '0', 'pwl', '0', '0', '1n', '1']
    assert by['m$615'].words == ['30', '57', '57', '57', 'pn2']
    assert by['m$615'].params == {'w': '37u', 'l': '2u', 'ad': '185p', 'as': '185p'}
    assert by['r2'].words == ['1', '2', '{RV*2}']
    assert by['x1'].words == ['a', 'b', 'sub'] and by['x1'].params == {'w': '2u', 'l': "'3*u0'"}
    assert by['m$615'].kind == 'm' and by['x1'].where.line == 7


def test_where_is_a_continued_lines_first_line(tmp_path):
    net = _read(tmp_path, """
        title
        * a comment
        R1 1
        + 0
        +
        + 1k
        C1 1 0 1p
        """)
    assert [(c.name, c.where.line, c.words) for c in net.cards] == [
        ('r1', 3, ['1', '0', '1k']), ('c1', 7, ['1', '0', '1p'])]
    assert str(net.cards[0].where) == 'a.cir:3'


def test_subcircuits_their_ports_parameters_and_models(tmp_path):
    """A subcircuit used before it is defined; its parameters' defaults
    in its scope; a model inside it is its own (found from inside, not from
    the top level), and a nested definition sees its definer's."""
    net = _read(tmp_path, """
        title
        X1 1 2 inv w=2u
        .subckt INV in out params: w=1u l=1u
        .param k = 'w/l'
        M1 out in 0 0 nch w={w} l={l}
        .model nch nmos level=1 vto=0.7
        .subckt inner a
        R1 a 0 {k*1k}
        .ends inner
        .ends INV
        .model top npn bf=100
        """)
    inv = net.subckts['inv']
    assert inv.ports == ['in', 'out'] and [c.name for c in inv.cards] == ['m1']
    assert net.subckts['inner'].parent is inv and [c.name for c in net.subckts['inner'].cards] == ['r1']
    assert net.model('nch') is None and net.model('nch', inv).type == 'nmos'
    assert net.model('NCH', net.subckts['inner']) is net.model('nch', inv)
    assert net.model('top', net.subckts['inner']).type == 'npn'
    assert net.deck.evaluate('k', inv.scope) == 1.0
    assert net.deck.evaluate('k', inv.scope, w=4e-6) == 4.0
    assert net.deck.values(net.model('nch', inv)) == {'level': 1.0, 'vto': 0.7}


def test_analyses_options_prints_temperatures(tmp_path):
    net = _read(tmp_path, """
        title
        R1 1 0 1
        .tran .05U 500U NOOP
        .op
        .options timeint method=gear reltol=1e-3
        .OPTIONS device temp=125
        .print tran precision=10 width=19 {V(8)+4} v(5) V(a, b) i(VSRC)
        .temp 27 50
        """)
    assert [(d.name, d.words) for d in net.analyses] == [('tran', ['.05u', '500u', 'noop']),
                                                         ('op', [])]
    assert [(d.words, d.params) for d in net.options] == [
        (['timeint'], {'method': 'gear', 'reltol': '1e-3'}), (['device'], {'temp': '125'})]
    (p,) = net.prints
    assert p.words == ['tran', '{V(8)+4}', 'v(5)', 'V(a, b)', 'i(VSRC)']
    assert p.params == {'precision': '10', 'width': '19'}
    assert net.temps[0].words == ['27', '50']


def test_initial_conditions_nodesets_and_comps(tmp_path):
    net = _read(tmp_path, """
        title
        R1 7 0 1
        .ic v(7)=0.0 V(OUT) = 2
        .nodeset v(9)=1.5
        *COMP V(8) reltol=0.02
        """)
    assert net.ic == {'7': '0.0', 'out': '2'} and net.nodeset == {'9': '1.5'}
    assert [c for c, _w in net.comps] == ['V(8) reltol=0.02']


def test_include_and_library_sections(tmp_path):
    (tmp_path / 'sub').mkdir()
    _write(tmp_path, 'sub/parts.inc', """
        R9 9 0 9
        .model d1 d is=1e-14
        """)
    _write(tmp_path, 'lib.lib', """
        .lib tt
        .param corner = 1
        .model qx npn bf='100*corner'
        .endl tt
        .lib ss
        .param corner = 0.5
        .model qx npn bf='100*corner'
        .endl ss
        """)
    net = _read(tmp_path, """
        title
        .include sub/parts.inc
        .lib lib.lib ss
        R1 1 0 1
        """)
    assert [c.name for c in net.cards] == ['r9', 'r1']
    assert net.cards[0].where.file.endswith('parts.inc') and net.cards[0].where.line == 1
    assert net.deck.values(net.model('d1')) == {'is': 1e-14}
    assert net.deck.values(net.model('qx')) == {'bf': 50.0}
    _write(tmp_path, 'loop.inc', '.include b.cir\n')
    with pytest.raises(SpiceNetlistError, match='circular'):
        _read(tmp_path, 'title\n.include loop.inc\n', name='b.cir')


def test_end_control_blocks_and_what_is_not_read(tmp_path):
    net = _read(tmp_path, """
        title
        R1 1 0 1
        .control
        run
        .endc
        .measure tran t1 when v(1)=0.5
        .step r1 1 2 1
        .end
        R2 1 0 2
        """)
    assert [c.name for c in net.cards] == ['r1']
    assert [(u, w.line) for u, w in net.unsupported] == [('.measure', 6), ('.step', 7)]


@pytest.mark.parametrize('text, says', [
    ('title\n.ends\n', '.ends without .subckt'),
    ('title\n.subckt a x\nR1 x 0 1\n', 'without .ends'),
    ('title\n.subckt a x\n.ends\n.subckt A y\n.ends\n', 'defined again'),
    ('title\nR1 1 0 r=\n', 'malformed assignment'),
    ('title\n.control\nrun\n', 'without .endc'),
    ('title\n.ic 7=0\n', 'cannot read .ic'),
    ("title\nR1 1 0 '1k\n", 'cannot read'),
], ids=['ends', 'subckt', 'duplicate', 'assignment', 'control', 'ic', 'quote'])
def test_what_cannot_be_read_is_refused_saying_where(tmp_path, text, says):
    with pytest.raises(SpiceNetlistError, match=says):
        _read(tmp_path, text)


def test_where_names_the_file_and_line():
    assert str(Where(os.path.join('x', 'y.cir'), 12)) == 'y.cir:12'
