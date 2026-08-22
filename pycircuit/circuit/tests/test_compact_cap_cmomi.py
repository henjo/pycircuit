"""`CapCmomi` against IHP's own compiled build of the same model.

`pycircuit/circuit/compact.py` translates the PDK's 314-line Verilog-A
interdigitated MoM capacitor into the Behavioural DSL.  These tests check
it against `data/cap_cmomi_ihp_ac.json` -- impedance from 1 MHz to
100 GHz, produced by sweeping the vendor's OSDI binary in ngspice
(`benchmarks/cap_cmomi_reference.py`).  The tests never invoke ngspice.

This is the first COMPLETE production compact model expressed in the DSL,
so it is also the first end-to-end check of the features added for the
PSP target: eight named branches (three pairs sharing a node pair),
`ddt` of a branch current as well as of a voltage, parameter-only
conditionals selecting fitted coefficient sets, `floor` on a parameter
expression, and clamps.

Two things are checked independently of the reference, so that a wrong
reference cannot make a wrong model pass: the low-frequency capacitance
is computed from the vendor's published tiling arithmetic by hand, and
the coefficient sets are checked to actually differ between layer counts.
"""
import json
import math
import os

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import IS, SubCircuit
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.compact import CapCmomi
from pycircuit.circuit.hdl import Node


DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    'data', 'cap_cmomi_ihp_ac.json')

#: Measured worst case over all ten configurations and all 51 frequencies
#: is 2.0e-8, and the reference itself carries only ~9 significant digits
#: (ngspice `wrdata`), so this is at the data's own precision.  The bound
#: is set an order above that and is still four orders below any
#: modelling difference that would matter.
ZTOL = 1e-6


@pytest.fixture(scope='module')
def ref():
    with open(DATA) as fh:
        return json.load(fh)


def _impedance(case, freqs):
    """Two-terminal Z(f) of `CapCmomi`, as the reference stores it.

    A 1 A AC source between the node and ground: SPICE's `I n+ n-`
    convention drives current from `n+` THROUGH the source to `n-`, so
    it draws from the node and the node voltage is `-Z`.  The reference
    generator applies the same negation, so both sides here are genuine
    impedances -- capacitive means negative imaginary part.
    """
    c = SubCircuit()
    n = c.add_node('PLUS')
    c['iac'] = IS(n, gnd, iac=1.0)
    c['X'] = CapCmomi(n, gnd, gnd, w=case['w'], l=case['l'],
                      mmin=case['mmin'], mmax=case['mmax'],
                      feed=case['feed'])
    c.update_iparv()
    res = AC(c, toolkit=numeric).solve(freqs=np.asarray(freqs, float))
    return -np.asarray(res.v(n, gnd))


CASES = ['w5_l5_N5_double', 'w5_l5_N5_same', 'w5_l5_N5_none',
         'w10_l20_N5_double', 'w30_l60_N5_double', 'w5_l5_N3_double',
         'w5_l5_N4_double', 'w5_l5_N2_same', 'w2_l2_N5_double',
         'w60_l3_N5_same']


class TestAgainstTheVendorBinary(object):

    def test_the_reference_is_what_it_claims(self, ref):
        assert 'cap_cmomi' in ref['source']
        assert 'ngspice' in ref['simulator']
        assert set(ref['cases']) == set(CASES)

    @pytest.mark.parametrize('name', CASES)
    def test_impedance_matches_over_the_whole_band(self, ref, name):
        case = ref['cases'][name]
        f = np.asarray(case['f'], float)
        zref = (np.asarray(case['z_real'], float)
                + 1j * np.asarray(case['z_imag'], float))
        z = _impedance(case, f)
        assert z.shape == zref.shape
        rel = np.abs(z - zref) / np.abs(zref)
        k = int(np.argmax(rel))
        assert rel.max() < ZTOL, (
            '%s: worst %.2e at %.4g Hz (got %r, want %r)'
            % (name, rel.max(), f[k], z[k], zref[k]))

    @pytest.mark.parametrize('name', CASES)
    def test_it_is_capacitive_at_low_frequency(self, ref, name):
        """A sign check the impedance comparison alone would not catch.

        Both sides could be negated together and still agree; this pins
        that the shared convention is the physical one.
        """
        case = ref['cases'][name]
        assert case['z_imag'][0] < 0
        assert _impedance(case, case['f'][:1])[0].imag < 0


class TestAgainstHandArithmetic(object):
    """Independent of the reference, so a wrong reference cannot hide."""

    @staticmethod
    def _expected_cap_fF(w, l, mmin, mmax, feed):
        """The vendor's published tiling arithmetic, done by hand.

        Unit cell 0.84 x 0.89 um; one y-row is spent on the feed.  Area
        density and per-micron feed capacitance are keyed by layer count.
        """
        lu, wu = l * 1e6, w * 1e6
        nlay = max(mmax - mmin + 1, 1)
        ax = max(math.floor(lu / 0.84 + 1e-6), 1)
        ay = max(math.floor(wu / 0.89 + 1e-6) - 1, 1)
        area = ax * ay * 0.84 * 0.89
        density = {2: 0.55, 3: 0.82, 4: 1.09}.get(min(nlay, 5), 1.36)
        if nlay <= 2:
            density = 0.55
        cfeed_pu = {3: 0.97, 4: 1.28}.get(nlay, 1.46 if nlay >= 5 else 0.70)
        feed_width = ay * 0.89 + 0.64
        return density * area + (cfeed_pu * feed_width if feed == 1 else 0.0)

    @pytest.mark.parametrize('name', CASES)
    def test_low_frequency_capacitance(self, ref, name):
        """C from the model's own impedance vs the tiling arithmetic.

        At 1 MHz the series R and L are utterly negligible against a
        multi-megohm reactance, so `Im(Z) = -1/(wC)` recovers `Cmain`
        directly.
        """
        case = ref['cases'][name]
        f0 = case['f'][0]
        z0 = _impedance(case, [f0])[0]
        got_fF = -1.0 / (2 * math.pi * f0 * z0.imag) * 1e15
        want_fF = self._expected_cap_fF(case['w'], case['l'], case['mmin'],
                                        case['mmax'], case['feed'])
        assert got_fF == pytest.approx(want_fF, rel=2e-3), name

    def test_the_layer_count_actually_changes_the_density(self, ref):
        """N=3/4/5 must give 0.82/1.09/1.36 fF/um^2 on one geometry."""
        caps = {}
        for name, nlay in (('w5_l5_N3_double', 3), ('w5_l5_N4_double', 4),
                           ('w5_l5_N5_double', 5)):
            case = ref['cases'][name]
            z0 = _impedance(case, [case['f'][0]])[0]
            caps[nlay] = -1.0 / (2 * math.pi * case['f'][0] * z0.imag)
        ## 5x5 um -> ax=5, ay=4 -> area = 14.952 um^2
        for nlay, density in ((3, 0.82), (4, 1.09), (5, 1.36)):
            assert caps[nlay] * 1e15 == pytest.approx(density * 14.952,
                                                      rel=1e-6)

    def test_the_same_side_feed_adds_its_capacitance(self, ref):
        """feed=same adds cfeed_pu * feed_width; double and none do not."""
        z_d = _impedance(ref['cases']['w5_l5_N5_double'], [1e6])[0]
        z_n = _impedance(ref['cases']['w5_l5_N5_none'], [1e6])[0]
        z_s = _impedance(ref['cases']['w5_l5_N5_same'], [1e6])[0]
        cap = lambda z: -1.0 / (2 * math.pi * 1e6 * z.imag)
        assert cap(z_d) == pytest.approx(cap(z_n), rel=1e-9)
        ## feed_width = 4*0.89 + 0.64 = 4.2 um, cfeed_pu(N=5) = 1.46
        assert (cap(z_s) - cap(z_d)) * 1e15 == pytest.approx(1.46 * 4.2,
                                                             rel=1e-6)


class TestStructure(object):
    """The topology the DSL built, which is where named branches matter."""

    def _elem(self, **kw):
        kw.setdefault('w', 5e-6)
        kw.setdefault('l', 5e-6)
        e = CapCmomi(Node('P'), Node('M'), Node('S'), **kw)
        e.update_iparv()
        return e

    def test_four_internal_nodes_and_two_branch_currents(self):
        """Only the two inductors are voltage-defined.

        Everything else is a current contribution and needs no unknown --
        so the eight named branches cost exactly two rows, not eight.
        """
        e = self._elem()
        assert len(e.nodes) == 7          # 3 terminals + nsk/ncap/nl/nox
        assert len(e.branches) == 2
        assert {b.name for b in e.branches} == {'lskin', 'lcore'}
        assert e.n == 9

    def test_the_parallel_pairs_did_not_merge(self):
        """`Lskin || Rskin` and `Rsub || Csub` share node pairs.

        Keyed on the node pair alone, each pair would collapse to one
        branch and lose a constitutive relation -- silently, since the
        remaining one still conducts.  The conductance matrix is the
        witness: the skin resistor and the substrate resistor must both
        appear.
        """
        e = self._elem(mmin=1, mmax=5, feed=2)
        G = np.asarray(e.G(np.zeros(e.n)), float)
        C = np.asarray(e.C(np.zeros(e.n)), float)
        ## Rskin (PLUS-nsk), Rseries (nl-MINUS), Rsub (nox-SUB) and the two
        ## GLEAK anchors are all conductances; the four capacitors and two
        ## inductors are all reactances.
        assert np.count_nonzero(G) >= 16
        assert np.count_nonzero(C) >= 10

    def test_dc_is_well_posed_and_the_anchor_is_the_only_path(self):
        """The GLEAK anchors exist so a floating island is not singular.

        Trace the DC network with `MINUS` and `SUB` both grounded: the
        two inductors are shorts, so `nsk` sits on `PLUS` and `ncap`/`nl`
        sit on `MINUS`; the four capacitors are open, so the PLUS/nsk
        island touches ground ONLY through `br_lkp`.  The other anchor,
        `br_lkm`, spans MINUS-SUB and is shorted out here.

        So the answer is `1 nA / 1 pS`, through one anchor -- and it is
        negative because SPICE's `I n+ n-` draws current out of `n+`.
        Getting a finite number at all is the point; getting exactly this
        one shows the anchors are where the model puts them.
        """
        c = SubCircuit()
        n = c.add_node('P')
        c['iac'] = IS(n, gnd, i=1e-9)
        c['X'] = CapCmomi(n, gnd, gnd, w=5e-6, l=5e-6)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        v = float(res.v(n, gnd))
        assert np.isfinite(v)
        assert v == pytest.approx(-1e-9 / 1e-12, rel=1e-6)

    def test_subblock_has_no_electrical_effect(self):
        """Accepted for interface parity, as the vendor documents."""
        a = self._elem(subblock=0)
        b = self._elem(subblock=1)
        x = np.zeros(a.n)
        assert np.allclose(np.asarray(a.G(x), float),
                           np.asarray(b.G(x), float))
        assert np.allclose(np.asarray(a.C(x), float),
                           np.asarray(b.C(x), float))


class TestParameterBranches(object):
    """Every conditional in the model, reached and shown to do something."""

    def _cap_fF(self, **kw):
        c = SubCircuit()
        n = c.add_node('P')
        c['iac'] = IS(n, gnd, iac=1.0)
        c['X'] = CapCmomi(n, gnd, gnd, **kw)
        c.update_iparv()
        z = -np.asarray(AC(c, toolkit=numeric).solve(freqs=np.array([1e6]))
                        .v(n, gnd))[0]
        return -1.0 / (2 * math.pi * 1e6 * z.imag) * 1e15

    def test_every_density_branch_is_distinct(self):
        caps = [self._cap_fF(w=5e-6, l=5e-6, mmin=mmin, mmax=5, feed=2)
                for mmin in (4, 3, 2, 1)]        # N = 2, 3, 4, 5
        assert caps == sorted(caps), caps
        assert len(set(round(c, 6) for c in caps)) == 4

    def test_the_floor_tiling_is_a_staircase(self):
        """Capacitance steps at the 0.84 um pitch, it does not ramp."""
        caps = [self._cap_fF(w=5e-6, l=x * 1e-6, mmin=1, mmax=5, feed=2)
                for x in (4.20, 4.60, 5.03, 5.05, 5.90)]
        ## 4.20/0.84 = 5.0 exactly; 4.60 and 5.03 still floor to 5;
        ## 5.05/0.84 = 6.01 -> 6; 5.90 -> 7
        assert caps[0] == pytest.approx(caps[1], rel=1e-9)
        assert caps[1] == pytest.approx(caps[2], rel=1e-9)
        assert caps[3] > caps[2] * 1.15
        assert caps[4] > caps[3] * 1.15

    def test_an_exact_pitch_multiple_lands_on_the_upper_step(self):
        """The vendor's epsilon before `floor`, which is load-bearing.

        `l = 21 um` is exactly 25 pitches, and 21e-6*1e6/0.84 lands one
        ulp low without the nudge -- floor would then give 24 and the
        model would disagree with the PCell by a whole row.
        """
        got = self._cap_fF(w=5e-6, l=21e-6, mmin=1, mmax=5, feed=2)
        ## 25 * 4 cells * 0.84 * 0.89 * 1.36 fF/um^2
        assert got == pytest.approx(25 * 4 * 0.84 * 0.89 * 1.36, rel=1e-9)

    def test_the_feed_code_switches_the_core_inductance(self):
        """double takes the loop fit, same/none the single-side one.

        The difference is 15-90x, so it shows up as a self-resonance far
        below the single-side one on a large device.
        """
        f = np.array([3e10])
        base = dict(w=30e-6, l=60e-6, mmin=1, mmax=5)
        zs = {}
        for feed in (0, 1, 2):
            c = SubCircuit()
            n = c.add_node('P')
            c['iac'] = IS(n, gnd, iac=1.0)
            c['X'] = CapCmomi(n, gnd, gnd, feed=feed, **base)
            c.update_iparv()
            zs[feed] = -np.asarray(AC(c, toolkit=numeric).solve(freqs=f)
                                   .v(n, gnd))[0]
        ## feed=double has resonated (inductive); the others have not
        assert zs[2].imag > 0, zs[2]
        assert zs[0].imag < 0, zs[0]
        assert zs[1].imag < 0, zs[1]
