"""Reading real foundry model cards.

Two layers.  The first uses small decks written inline, so every rule of
the format is pinned independently of any PDK being installed.  The
second reads the actual IHP Open PDK if it is present and skips if not --
that one is the reason the module exists, but it cannot be a hard
dependency of the test suite.
"""
import ast
import os
import re
import textwrap

import pytest

from pycircuit.utilities import spicecard
from pycircuit.utilities.spicecard import SpiceCardError

PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
needs_pdk = pytest.mark.skipif(not os.path.isdir(PDK),
                               reason='IHP Open PDK not installed')


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(textwrap.dedent(text).lstrip())
    return str(p)


class TestParameters(object):

    def test_plain_values_and_references(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param a = 2.5
            .param b = 'a * 4'
            .model m1 nmos vth='b + 1' w=3
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['vth'] == pytest.approx(11.0, rel=1e-12, abs=0.0)
        assert p['w'] == pytest.approx(3.0, rel=1e-12, abs=0.0)

    def test_a_wrapped_comma_separated_card_reads_like_a_plain_one(
            self, tmp_path):
        """The form vendor macromodels use: `TYPE(a=1, b=2)`, the
        parenthesis glued to the type or not, over continuation lines.
        Commas inside an expression are the expression's."""
        f = _write(tmp_path, 'a.sp', """
            .model dx D(IS=1E-14,RS=5)
            .model pox PMOS (LEVEL=2,KP=10E-6,VTO=-0.328)
            .model np NPN(Bf=1200 Vaf=140
            + Ikf=100m)
            .model m1 nmos a='max(1, 2)' b={min(3, 4)}, c=max(5,6)
            """)
        d = spicecard.read(f)
        assert d.models['dx'].type == 'd'
        assert d.model_params('dx') == {'is': 1e-14, 'rs': 5.0}
        assert d.models['pox'].type == 'pmos'
        assert d.model_params('pox') == pytest.approx(
            {'level': 2.0, 'kp': 10e-6, 'vto': -0.328}, rel=1e-12, abs=0.0)
        assert d.model_params('np') == pytest.approx(
            {'bf': 1200.0, 'vaf': 140.0, 'ikf': 0.1}, rel=1e-12, abs=0.0)
        assert d.model_params('m1') == {'a': 2.0, 'b': 3.0, 'c': 6.0}

    def test_engineering_suffixes(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos a=1u b=2n c=3p d=1meg e=4k g=5m h=6f
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['a'] == pytest.approx(1e-6, rel=1e-12, abs=0.0)
        assert p['b'] == pytest.approx(2e-9, rel=1e-12, abs=0.0)
        assert p['c'] == pytest.approx(3e-12, rel=1e-12, abs=0.0)
        assert p['d'] == pytest.approx(1e6, rel=1e-12, abs=0.0), 'meg must beat m'
        assert p['e'] == pytest.approx(4e3, rel=1e-12, abs=0.0)
        assert p['g'] == pytest.approx(5e-3, rel=1e-12, abs=0.0)
        assert p['h'] == pytest.approx(6e-15, rel=1e-12, abs=0.0)

    def test_exponent_notation_is_not_eaten_by_the_suffix_rule(self,
                                                              tmp_path):
        """`1e-6` must stay 1e-6, not become 1 with an `e` suffix."""
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos a=1e-6 b=2.5E+3 c=1.44e-15
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['a'] == pytest.approx(1e-6, rel=1e-12, abs=0.0)
        assert p['b'] == pytest.approx(2500.0, rel=1e-12, abs=0.0)
        assert p['c'] == pytest.approx(1.44e-15, rel=1e-12, abs=0.0)

    def test_functions(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos
            + a='max(3, 7)' b='min(3, 7)' c='sqrt(16)'
            + d='pow(2, 10)' e='abs(-4)' g='log(100)' h='ln(1)'
            """)
        p = spicecard.read(f).model_params('m1')
        assert (p['a'], p['b'], p['c']) == (7.0, 3.0, 4.0)
        assert (p['d'], p['e'], p['g'], p['h']) == (1024.0, 4.0, 2.0, 0.0)

    def test_statistical_functions_give_the_nominal_value(self, tmp_path):
        """Without a random draw the centre is the only honest answer."""
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos a='agauss(5, 1, 3)' b='gauss(2, 0.1, 1)'
            + c='aunif(7, 2)' d='unif(9, 0.5)'
            """)
        p = spicecard.read(f).model_params('m1')
        assert (p['a'], p['b'], p['c'], p['d']) == (5.0, 2.0, 7.0, 9.0)

    def test_comments_and_continuations(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            * a full-line comment
            .param a = 1 ; trailing comment
            .param b = 2 $ another style
            .model m1 nmos
            + x='a+b'
            * comment between continuations is fine
            + y=10
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['x'] == pytest.approx(3.0, rel=1e-12, abs=0.0)
        assert p['y'] == pytest.approx(10.0, rel=1e-12, abs=0.0)

    def test_a_semicolon_inside_a_quoted_expression_survives(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param a = 'max(1, 2)'
            .model m1 nmos x='a * 3'
            """)
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(6.0, rel=1e-12, abs=0.0)


class TestLibSections(object):

    def test_a_section_is_opt_in(self, tmp_path):
        """Unrequested `.LIB` blocks are skipped, not concatenated.

        A corner file defines the SAME names differently per section, so
        reading them all would silently give whichever came last.
        """
        f = _write(tmp_path, 'c.lib', """
            .LIB tt
            .param k = 1.0
            .ENDL tt
            .LIB ss
            .param k = 0.5
            .ENDL ss
            .model m1 nmos x='k'
            """)
        assert spicecard.read(f, section='tt').model_params('m1')['x'] == 1.0
        assert spicecard.read(f, section='ss').model_params('m1')['x'] == 0.5
        with pytest.raises(SpiceCardError, match='undefined parameter'):
            spicecard.read(f).model_params('m1')

    def test_lib_with_a_file_argument_is_an_include(self, tmp_path):
        _write(tmp_path, 'inner.lib', """
            .LIB tt
            .param k = 3.0
            .ENDL tt
            """)
        f = _write(tmp_path, 'top.sp', """
            .lib inner.lib tt
            .model m1 nmos x='k * 2'
            """)
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(6.0, rel=1e-12, abs=0.0)

    def test_include_paths_are_relative_to_the_including_file(self, tmp_path):
        sub = tmp_path / 'sub'
        sub.mkdir()
        (sub / 'p.lib').write_text('.param k = 4.0\n')
        f = _write(tmp_path, 'top.sp', """
            .include sub/p.lib
            .model m1 nmos x='k'
            """)
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(4.0, rel=1e-12, abs=0.0)


class TestScoping(object):

    def test_a_subckt_parameter_shadows_the_global_one(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param w = 1.0
            .subckt cell a b w=2.0
            .model inner nmos x='w'
            .ends
            .model outer nmos x='w'
            """)
        d = spicecard.read(f)
        assert d.model_params('inner')['x'] == pytest.approx(2.0, rel=1e-12, abs=0.0)
        assert d.model_params('outer')['x'] == pytest.approx(1.0, rel=1e-12, abs=0.0)

    def test_an_override_beats_the_subckt_default(self, tmp_path):
        """That is what an instance parameter IS.

        `.subckt cell ... w=0.5u` is a default the instantiation
        replaces; the card downstream must see the replacement.
        """
        f = _write(tmp_path, 'a.sp', """
            .subckt cell a b w=0.5u ng=1
            .param area = 'w * 10'
            .model inner nmos x='area / ng'
            .ends
            """)
        d = spicecard.read(f)
        assert d.model_params('inner')['x'] == pytest.approx(5e-6, rel=1e-12, abs=0.0)
        assert d.model_params('inner', w=2e-6)['x'] == pytest.approx(2e-5, rel=1e-12, abs=0.0)
        assert d.model_params('inner', w=2e-6, ng=4)['x'] == pytest.approx(5e-6, rel=1e-12, abs=0.0)

    def test_overrides_are_case_insensitive(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .subckt cell a b W=1.0
            .model inner nmos x='w'
            .ends
            """)
        d = spicecard.read(f)
        assert d.model_params('inner', W=7.0)['x'] == pytest.approx(7.0, rel=1e-12, abs=0.0)
        assert d.model_params('INNER', w=8.0)['x'] == pytest.approx(8.0, rel=1e-12, abs=0.0)


class TestRefusals(object):
    """Every failure names what is wrong, rather than yielding a number."""

    def test_unknown_model(self, tmp_path):
        f = _write(tmp_path, 'a.sp', '.model m1 nmos x=1\n')
        with pytest.raises(SpiceCardError, match='no model'):
            spicecard.read(f).model_params('nope')

    def test_undefined_parameter(self, tmp_path):
        f = _write(tmp_path, 'a.sp', ".model m1 nmos x='missing * 2'\n")
        with pytest.raises(SpiceCardError, match='undefined parameter'):
            spicecard.read(f).model_params('m1')

    def test_circular_definition(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param a = 'b + 1'
            .param b = 'a + 1'
            .model m1 nmos x='a'
            """)
        with pytest.raises(SpiceCardError, match='circular'):
            spicecard.read(f).model_params('m1')

    def test_circular_include(self, tmp_path):
        _write(tmp_path, 'b.sp', '.include a.sp\n')
        f = _write(tmp_path, 'a.sp', '.include b.sp\n')
        with pytest.raises(SpiceCardError, match='circular'):
            spicecard.read(f)

    def test_missing_file(self, tmp_path):
        f = _write(tmp_path, 'a.sp', '.include nope.lib\n')
        with pytest.raises(SpiceCardError, match='no such file'):
            spicecard.read(f)

    def test_division_by_zero_is_reported_as_such(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param z = 0
            .model m1 nmos x='1 / z'
            """)
        with pytest.raises(SpiceCardError, match='division by zero'):
            spicecard.read(f).model_params('m1')

    def test_an_expression_cannot_reach_arbitrary_python(self, tmp_path):
        """A vendor file is data.  Reading it must not run it."""
        f = _write(tmp_path, 'a.sp',
                   ".model m1 nmos x='__import__(\"os\").system(\"true\")'\n")
        with pytest.raises(SpiceCardError):
            spicecard.read(f).model_params('m1')


class TestSpiceLexicalRules:
    """SPICE's lexical rules (2026-10-07, the SPICE benchmark plan's
    stage 1) -- each failed on the parent."""

    @pytest.mark.parametrize('text, digits, mult', [
        ('1MEG', 1, 1e6), ('2.5Meg', 2.5, 1e6), ('1meg', 1, 1e6),
        ('1M', 1, 1e-3), ('1Mohm', 1, 1e-3), ('1MEGohm', 1, 1e6),
        ('2mil', 2, 25.4e-6), ('2MILS', 2, 25.4e-6),
        ('1.5P', 1.5, 1e-12), ('10pF', 10, 1e-12), ('10PF', 10, 1e-12),
        ('7.0F', 7.0, 1e-15), ('1.8mA', 1.8, 1e-3), ('10ns', 10, 1e-9),
        ('10NS', 10, 1e-9), ('37U', 37, 1e-6), ('4.7u', 4.7, 1e-6),
        ('3k', 3, 1000), ('2K', 2, 1000), ('1T', 1, 1e12), ('3G', 3, 1e9),
        ('2a', 2, 1e-18), ('5V', 5, None), ('.05V', 0.05, None),
        ('0v', 0, None), ('5Hz', 5, None), ('-2.5e-3', -2.5e-3, None),
        ('+3k', 3, 1000), ('-4u', -4, 1e-6), ('1e-6u', 1e-6, 1e-6),
        ('2.5E+3', 2500.0, None), ('007', 7, None),
    ])
    def test_a_number_is_its_digits_times_its_scale_factor(self, text, digits, mult):
        """A scale factor in any case (`meg` and `mil` before `m`, which
        is milli -- so `1Mohm` is a milliohm, as in SPICE); the letters
        after it a unit, ignored (`7.0F` is seven femto)."""
        want = float(digits if mult is None else digits * mult)
        assert spicecard.number(text) == want

    @pytest.mark.parametrize('text', ['', 'abc', 'u1', '1.2.3', '1 2', '--1', '1u5'])
    def test_what_is_not_a_number_is_refused(self, text):
        with pytest.raises(SpiceCardError, match='not a number'):
            spicecard.number(text)

    def test_a_card_reads_them_as_number_does(self, tmp_path):
        vals = ['1MEG', '1.5P', '10pF', '1.8mA', '5V', '7.0F', '10NS', '2MILS', '1Mohm']
        cards = ' '.join(f'p{k}={v}' for k, v in enumerate(vals))
        f = _write(tmp_path, 'a.sp', f".model m1 nmos {cards}\n.param e = '2*1MEG + 10pF'\n")
        d = spicecard.read(f)
        p = d.model_params('m1')
        assert [p[f'p{k}'] for k in range(len(vals))] == [spicecard.number(v) for v in vals]
        assert d.param('e') == 2 * 1e6 + 10 * 1e-12

    def test_model_parameters_in_parentheses(self, tmp_path):
        """Glued to the type or not; the closing one on a line of its own;
        an empty continuation line."""
        f = _write(tmp_path, 'a.sp', """
            .model q1 NPN(BF=100 IS=1e-16)
            .model q2 npn ( bf = 50
            + is=2e-16 )
            .MODEL d1 D(
            + IS=14.34f
            +
            + RS=10
            + )
            """)
        d = spicecard.read(f)
        assert [d.models[m].type for m in ('q1', 'q2', 'd1')] == ['npn', 'npn', 'd']
        assert d.model_params('q1') == {'bf': 100.0, 'is': 1e-16}
        assert d.model_params('q2') == {'bf': 50.0, 'is': 2e-16}
        assert d.model_params('d1') == {'is': 14.34 * 1e-15, 'rs': 10.0}

    def test_subckt_ports_end_at_params_or_the_first_assignment(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .subckt inv in out params: w=1u l=2u
            .model m1 nmos x='w/l'
            .ends
            .SUBCKT cell a b c W = 2u
            .ends
            """)
        d = spicecard.read(f)
        assert d.subckt_ports == {'inv': ['in', 'out'], 'cell': ['a', 'b', 'c']}
        assert d.model_params('m1')['x'] == 1e-6 / 2e-6

    def test_a_dollar_starts_a_comment_only_where_a_token_starts(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .model q$1 npn bf=50 $ a comment
            .model q2 npn bf=60 ,$ after a comma too
            $ a whole line
            .param a = 2 ; and a semicolon
            .param b = '3'$ c=4
            """)
        d = spicecard.read(f)
        assert sorted(d.models) == ['q$1', 'q2']
        assert d.model_params('q$1') == {'bf': 50.0} and d.model_params('q2') == {'bf': 60.0}
        assert d.param('a') == 2.0
        assert (d.param('b'), d.param('c')) == (3.0, 4.0)    # after a quote: no comment


class TestExpressions:
    """Expressions are parsed (2026-10-07): SPICE's operators with C's
    precedence where Python has none -- each failed on the parent -- and
    Python's where it has, as the regex translation had them."""

    @pytest.mark.parametrize('expr, value', [
        ('1 > 0 ? 2 : 3', 2.0),
        ('0 ? 2 : 3', 3.0),
        ('0 ? 1 : 0 ? 2 : 3', 3.0),           # right-associative
        ('1 ? 0 ? 4 : 5 : 6', 5.0),           # a conditional in the middle
        ('1 + 1 > 2 ? 7 : 8', 8.0),           # the loosest operator
        ('1 || 0 && 0', 1.0),                 # && binds tighter than ||
        ('0 && 1 || 1', 1.0),
        ('!0', 1.0), ('!3', 0.0),
        ('!0 + 1', 2.0),                      # C's unary !, not Python's loose `not`
        ('2 <> 3', 1.0), ('2 <> 2', 0.0),
        ('if(1 > 2, 4, 5)', 5.0),
        ('2*1MEG + 10pF', 2 * 1e6 + 10 * 1e-12),
        ('agauss(1, 0.1, (1 != 1 ? 0 : 1))', 1.0),   # the IHP mismatch cards' form
    ])
    def test_spice_operators(self, tmp_path, expr, value):
        f = _write(tmp_path, 'a.sp', f".model m1 nmos x='{expr}'\n")
        assert spicecard.read(f).model_params('m1')['x'] == value

    @pytest.mark.parametrize('expr, value', [
        ('-2^2', -4.0), ('2^3^2', 512.0), ('2**-1', 0.5), ('-2**2', -4.0),
        ('7 % 4', 3.0), ('1 + 2 * 3', 7.0), ('(1 + 2) * 3', 9.0),
        ('IF(1, 4, 5)', 4.0), ('1 < 2 and not 3 < 2', 1.0),
    ])
    def test_python_precedence_where_python_has_the_operator(self, tmp_path, expr, value):
        f = _write(tmp_path, 'a.sp', f".model m1 nmos x='{expr}'\n")
        assert spicecard.read(f).model_params('m1')['x'] == value

    def test_an_if_chain_selects_its_branch(self, tmp_path):
        for k, want in ((2, 10.0), (3, 20.0), (4, 30.0)):
            f = _write(tmp_path, f'a{k}.sp', f"""
                .param k = {k}
                .if (k == 2)
                .param a = 10
                .elseif (k == 3 && 1)
                .param a = 20
                .else
                .param a = 30
                .endif
                .model m1 nmos x='a'
                """)
            assert spicecard.read(f).model_params('m1')['x'] == want, k

    def test_a_parameter_named_like_a_python_word(self, tmp_path):
        """`as`, `is`, `lambda` -- a source area, a saturation current, a
        channel-length modulation -- read as parameters (the IHP cards'
        `.if (as <= 1e-50)` could not be read)."""
        f = _write(tmp_path, 'a.sp', """
            .subckt dio a b as=2p is=1e-14 lambda=0.02
            .if (as <= 1e-50)
            .param k = 0
            .else
            .param k = 1
            .endif
            .model d1 d x='as*k + is' y='lambda*2'
            .ends
            """)
        d = spicecard.read(f)
        assert d.model_params('d1') == {'x': 2 * 1e-12 * 1.0 + 1e-14, 'y': 0.04}
        assert d.model_params('d1', **{'as': 0.0})['x'] == 1e-14
        with pytest.raises(SpiceCardError, match="undefined parameter 'is'"):
            spicecard.read(_write(tmp_path, 'b.sp', ".model d2 d x='is'\n")).model_params('d2')

    @pytest.mark.parametrize('expr', [
        '(lambda: 1)()', '[1][0]', '1 .real', '().__class__', "'a'", '1 if 1 else 0',
    ])
    def test_only_numbers_names_calls_and_operators_parse(self, tmp_path, expr):
        """A card is data: on the parent each of these reached `eval` (a
        closed namespace stops a name, not a lambda or an attribute walk
        from a literal)."""
        f = _write(tmp_path, 'a.sp', f'.model m1 nmos x={{{expr}}}\n')
        with pytest.raises(SpiceCardError, match='cannot parse'):
            spicecard.read(f).model_params('m1')


#: The translation the parser replaced (`Deck._pythonise` at 0c1b2174,
#: verbatim): wherever it wrote Python, the parser's text must parse to
#: the same tree -- an expression it read means what it meant.
_OLD_SUFFIX = [('meg', 1e6), ('mil', 25.4e-6), ('t', 1e12), ('g', 1e9),
               ('k', 1e3), ('m', 1e-3), ('u', 1e-6), ('n', 1e-9),
               ('p', 1e-12), ('f', 1e-15), ('a', 1e-18)]
_OLD_NUM = re.compile(
    r'\b(\d+\.?\d*(?:[eE][-+]?\d+)?|\.\d+(?:[eE][-+]?\d+)?)'
    r'(meg|mil|t|g|k|m|u|n|p|f|a)?([a-zA-Z_]*)\b')


def _old_pythonise(expr):
    e = expr.strip()
    if e[:1] in "'\"" and e[-1:] == e[:1]:  # noqa: SIM114 (verbatim)
        e = e[1:-1]
    elif e.startswith('{') and e.endswith('}'):
        e = e[1:-1]
    e = e.replace('^', '**')

    def num(m):
        base, suf, tail = m.group(1), m.group(2), m.group(3)
        if tail:
            return m.group(0)
        if not suf:
            return base
        for name, mult in _OLD_SUFFIX:
            if suf.lower() == name:
                return '(%s*%g)' % (base, mult)  # noqa: UP031 (verbatim)
        return m.group(0)

    return _OLD_NUM.sub(num, e)


def _tree(text):
    try:
        return ast.dump(ast.parse(text, mode='eval'))
    except SyntaxError:
        return None


@needs_pdk
def test_every_pdk_expression_the_regex_read_parses_to_its_tree():
    """Every expression of the IHP cards -- `.param`, `.model`, instance
    and `.subckt` values, `.if` conditions: where the regex translation
    wrote Python the parser's parses to the same tree (measured
    2026-10-07: 1849 the same, none different); the parser reads the
    rest (the 38 with `?:`, an upper-case `1G`, and `as`)."""
    exprs = set()
    for fn in sorted(os.listdir(PDK)):
        for text, _here in spicecard._logical_lines(os.path.join(PDK, fn)):
            head = text.split(None, 1)[0].lower()
            if head in ('.if', '.elseif'):
                exprs.add(text[text.find('(') + 1:text.rfind(')')])
            body = spicecard._model_card(text)[2] if head == '.model' else text
            exprs.update(raw for _name, raw in spicecard._assignments(body))
    same, new, differ = 0, 0, []
    for raw in sorted(exprs):
        old, now = _tree(_old_pythonise(raw)), _tree(spicecard.Deck._pythonise(raw))
        assert now is not None, raw
        if old is None:
            new += 1
        elif old == now:
            same += 1
        else:
            differ.append(raw)
    assert not differ, differ[:5]
    assert same > 1800 and new >= 40, (same, new)


@needs_pdk
class TestTheRealPDK(object):
    """The card this module exists for: PSP103, 359 parameters."""

    CORNER = os.path.join(PDK, 'cornerMOSlv.lib')
    INST = dict(w=1e-6, l=0.13e-6, ng=1, m=1, pre_layout=1)

    def test_the_whole_psp103_card_resolves_to_numbers(self):
        d = spicecard.read(self.CORNER, section='mos_tt')
        p = d.model_params('sg13g2_lv_nmos_psp', **self.INST)
        assert len(p) > 350
        assert all(isinstance(v, float) for v in p.values())
        assert p['level'] == pytest.approx(103.6, rel=1e-12, abs=0.0)
        assert p['type'] == pytest.approx(1.0, rel=1e-12, abs=0.0)

    def test_corner_multipliers_are_actually_applied(self):
        """The card holds `'-0.25737*sg13g2_lv_nmos_dphibo'`.

        If the corner section were not followed, this would either fail
        to resolve or come back as the bare coefficient.  It comes back
        multiplied, and differently per corner.
        """
        got = {}
        for sec in ('mos_tt', 'mos_ss', 'mos_ff'):
            d = spicecard.read(self.CORNER, section=sec)
            got[sec] = d.model_params('sg13g2_lv_nmos_psp', **self.INST)
        assert got['mos_tt']['dphibo'] == pytest.approx(-0.25737 * 0.9915, rel=1e-12, abs=0.0)
        assert len({round(g['dphibo'], 9) for g in got.values()}) == 3
        assert len({round(g['rsw1'], 9) for g in got.values()}) == 3

    def test_instance_parameters_reach_the_card(self):
        """`dlq` reads `pre_layout`; `cfrw` divides by `ng`."""
        d = spicecard.read(self.CORNER, section='mos_tt')
        pre1 = d.model_params('sg13g2_lv_nmos_psp',
                              **dict(self.INST, pre_layout=1))
        pre0 = d.model_params('sg13g2_lv_nmos_psp',
                              **dict(self.INST, pre_layout=0))
        assert pre0['dlq'] == pytest.approx(pre1['dlq'] - 2e-8, rel=1e-12, abs=0.0)
        for ng in (1, 2, 4):
            p = d.model_params('sg13g2_lv_nmos_psp',
                               **dict(self.INST, ng=ng))
            assert p['cfrw'] == pytest.approx(2e-16 / ng, rel=1e-12, abs=0.0)

    def test_the_rf_cards_read(self):
        """`sg13g2_lv_nmos_psp_rf`'s `dlq` holds `(ng<3 ? 4e-08 : 0)`:
        SPICE's conditional, which the parent could not read."""
        d = spicecard.read(self.CORNER, section='mos_tt')
        for ng, extra in ((1, 4e-08), (4, 0)):
            p = d.model_params('sg13g2_lv_nmos_psp_rf', **dict(self.INST, ng=ng, rfmode=1))
            assert p['dlq'] == -1.3721e-08 - ((1 - 1) * 2e-08) + 1 * (-1.5368e-08 + extra)

    @pytest.mark.parametrize('lib,section,model', [
        ('cornerMOSlv.lib', 'mos_tt', 'sg13g2_lv_pmos_psp'),
        ('cornerMOShv.lib', 'mos_tt', 'sg13g2_hv_nmos_psp'),
        ('cornerCAP.lib', 'cap_typ', 'cap_cmomi_mod'),
        ('cornerRES.lib', 'res_typ', 'rmod_rsil'),
    ])
    def test_every_model_family_in_the_pdk_reads(self, lib, section, model):
        path = os.path.join(PDK, lib)
        if not os.path.exists(path):
            pytest.skip('%s not in this PDK checkout' % lib)
        d = spicecard.read(path, section=section)
        assert model in d.models
        p = d.model_params(model, **self.INST)
        assert all(isinstance(v, float) for v in p.values())

    def test_the_psp_card_values_are_physically_sane(self):
        """A spot check that the numbers mean something.

        Oxide thickness of a 130 nm node is a couple of nanometres, and
        the flat-band voltage is order -1 V.  Wrong scoping tends to
        produce values that are off by the multiplier, which this catches.
        """
        d = spicecard.read(self.CORNER, section='mos_tt')
        p = d.model_params('sg13g2_lv_nmos_psp', **self.INST)
        assert 1e-9 < p['toxo'] < 5e-9
        assert -2.0 < p['vfbo'] < 0.0
        assert 1e22 < p['nsubo'] < 1e24
