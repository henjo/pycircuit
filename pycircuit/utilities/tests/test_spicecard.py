"""Reading real foundry model cards.

Two layers.  The first uses small decks written inline, so every rule of
the format is pinned independently of any PDK being installed.  The
second reads the actual IHP Open PDK if it is present and skips if not --
that one is the reason the module exists, but it cannot be a hard
dependency of the test suite.
"""
import os
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
        assert p['vth'] == pytest.approx(11.0)
        assert p['w'] == pytest.approx(3.0)

    def test_engineering_suffixes(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos a=1u b=2n c=3p d=1meg e=4k g=5m h=6f
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['a'] == pytest.approx(1e-6)
        assert p['b'] == pytest.approx(2e-9)
        assert p['c'] == pytest.approx(3e-12)
        assert p['d'] == pytest.approx(1e6), 'meg must beat m'
        assert p['e'] == pytest.approx(4e3)
        assert p['g'] == pytest.approx(5e-3)
        assert p['h'] == pytest.approx(6e-15)

    def test_exponent_notation_is_not_eaten_by_the_suffix_rule(self,
                                                              tmp_path):
        """`1e-6` must stay 1e-6, not become 1 with an `e` suffix."""
        f = _write(tmp_path, 'a.sp', """
            .model m1 nmos a=1e-6 b=2.5E+3 c=1.44e-15
            """)
        p = spicecard.read(f).model_params('m1')
        assert p['a'] == pytest.approx(1e-6)
        assert p['b'] == pytest.approx(2500.0)
        assert p['c'] == pytest.approx(1.44e-15)

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
        assert p['x'] == pytest.approx(3.0)
        assert p['y'] == pytest.approx(10.0)

    def test_a_semicolon_inside_a_quoted_expression_survives(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .param a = 'max(1, 2)'
            .model m1 nmos x='a * 3'
            """)
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(6.0)


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
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(6.0)

    def test_include_paths_are_relative_to_the_including_file(self, tmp_path):
        sub = tmp_path / 'sub'
        sub.mkdir()
        (sub / 'p.lib').write_text('.param k = 4.0\n')
        f = _write(tmp_path, 'top.sp', """
            .include sub/p.lib
            .model m1 nmos x='k'
            """)
        assert spicecard.read(f).model_params('m1')['x'] == pytest.approx(4.0)


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
        assert d.model_params('inner')['x'] == pytest.approx(2.0)
        assert d.model_params('outer')['x'] == pytest.approx(1.0)

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
        assert d.model_params('inner')['x'] == pytest.approx(5e-6)
        assert d.model_params('inner', w=2e-6)['x'] == pytest.approx(2e-5)
        assert d.model_params('inner', w=2e-6, ng=4)['x'] == pytest.approx(5e-6)

    def test_overrides_are_case_insensitive(self, tmp_path):
        f = _write(tmp_path, 'a.sp', """
            .subckt cell a b W=1.0
            .model inner nmos x='w'
            .ends
            """)
        d = spicecard.read(f)
        assert d.model_params('inner', W=7.0)['x'] == pytest.approx(7.0)
        assert d.model_params('INNER', w=8.0)['x'] == pytest.approx(8.0)


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
        assert p['level'] == pytest.approx(103.6)
        assert p['type'] == pytest.approx(1.0)

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
        assert got['mos_tt']['dphibo'] == pytest.approx(-0.25737 * 0.9915)
        assert len({round(g['dphibo'], 9) for g in got.values()}) == 3
        assert len({round(g['rsw1'], 9) for g in got.values()}) == 3

    def test_instance_parameters_reach_the_card(self):
        """`dlq` reads `pre_layout`; `cfrw` divides by `ng`."""
        d = spicecard.read(self.CORNER, section='mos_tt')
        pre1 = d.model_params('sg13g2_lv_nmos_psp',
                              **dict(self.INST, pre_layout=1))
        pre0 = d.model_params('sg13g2_lv_nmos_psp',
                              **dict(self.INST, pre_layout=0))
        assert pre0['dlq'] == pytest.approx(pre1['dlq'] - 2e-8)
        for ng in (1, 2, 4):
            p = d.model_params('sg13g2_lv_nmos_psp',
                               **dict(self.INST, ng=ng))
            assert p['cfrw'] == pytest.approx(2e-16 / ng)

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
