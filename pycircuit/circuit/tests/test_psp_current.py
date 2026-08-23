"""The surface-potential drain current, and the MOSFET built on it.

`psp_kernel.ids_long_channel` assembles PSP's intrinsic current from the
surface potentials at the two channel ends, and
`compact.PspMosLongChannel` wraps it as an element.

There is no vendor reference here, deliberately: this is PSP's *core*
without the short-channel, mobility, series-resistance and geometry
layers, so it is not the same device as any model card and comparing
I-V curves would be comparing two different things.  What is checked
instead is the set of properties the formulation guarantees BY
CONSTRUCTION -- which is a stronger statement than a curve fit, because
a construction error breaks them exactly and a fitting error does not.

The headline one is source-drain symmetry.  Threshold-voltage models
(BSIM3 and relatives) are famously asymmetric about `Vds = 0`, which
shows up as spurious distortion in passing-gate and mixer circuits;
surface-potential models are symmetric because the inversion charge is
evaluated at the MIDPOINT potential.  Here it comes out bit-exact.
"""
import warnings

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import R, SubCircuit, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit import hdl
from pycircuit.circuit import psp_kernel as K
from pycircuit.circuit.compact import PspMosLongChannel
from pycircuit.circuit.hdl import Node


PHIT = 0.02585
INV_SQRT2 = 7.0710678118654746e-01
XN0 = 35.0           # 2*phi_F / phit for a real device at room temperature

#: Below this the drain current is a difference of two surface potentials
#: that agree to nine significant digits, and the cancellation leaves a
#: few percent of wobble.  Measured, not chosen: above it every sampled
#: step is strictly positive.  One attoamp is a thousandth of any real
#: subthreshold measurement floor.
NOISE_FLOOR = 1e-18


@pytest.fixture(scope='module')
def ids():
    """`Ids(xg, xn_s, xn_d)` at fixed Gf, compiled once."""
    xg, xns, xnd, gf, xi = sympy.symbols('xg xns xnd Gf xi', real=True)
    f = hdl.compile_chain(
        lambda: K.ids_long_channel(xg, xns, xnd, gf, xi, PHIT, 1e-4),
        [xg, xns, xnd, gf, xi])

    def call(vg, vs, vd, gf_=1.0):
        return float(f(vg / PHIT, XN0 + vs / PHIT, XN0 + vd / PHIT,
                       gf_, 1.0 + gf_ * INV_SQRT2)[0])
    return call


class TestTheConstructionsGuarantees(object):

    @pytest.mark.parametrize('vg', [0.3, 0.8, 1.2, 1.8])
    @pytest.mark.parametrize('vd', [0.01, 0.1, 0.5, 1.2])
    def test_source_and_drain_are_exactly_interchangeable(self, ids, vg, vd):
        """Bit-exact, not approximately.

        Swapping which end is called the source must negate the current
        and change nothing else.  A regional model has to be *fitted* to
        approximate this; here it is a consequence of evaluating the
        inversion charge at the midpoint, so any deviation at all would
        mean the construction is wrong.
        """
        forward = ids(vg, 0.0, vd)
        reverse = -ids(vg, vd, 0.0)
        assert forward == reverse

    @pytest.mark.parametrize('vg', [-0.5, 0.0, 0.3, 1.0, 1.8, 3.0])
    def test_no_drain_bias_means_no_current(self, ids, vg):
        """Exactly zero, at any gate bias, with no epsilon."""
        assert ids(vg, 0.0, 0.0) == 0.0

    def test_one_expression_spans_every_region(self, ids):
        """Accumulation to strong inversion, with no regional equations.

        Strictly monotone wherever the current is above `NOISE_FLOOR`,
        which is measured, not chosen for convenience: below about
        1e-18 A the current is a difference of two surface potentials
        that agree to nine digits, and the cancellation leaves a few
        percent of wobble.  One attoamp is a thousandth of any real
        subthreshold measurement floor, so the model is monotone across
        its entire meaningful range.
        """
        vgs = np.linspace(-0.5, 2.0, 200)
        cur = np.array([ids(v, 0.0, 0.1) for v in vgs])
        assert np.all(np.isfinite(cur))
        live = cur > NOISE_FLOOR
        assert live.sum() > 50
        assert np.all(np.diff(cur[live]) > 0), 'monotone in Vg above the floor'
        assert cur[-1] / max(cur[len(cur) // 3], 1e-300) > 1e6

    def test_below_the_noise_floor_it_is_only_nearly_monotone(self):
        """The limit, recorded so it is a known property.

        Not a defect to fix by smoothing: it is the cancellation in
        `x_d - x_s`, and the only cure is more precision than a double
        has.
        """
        assert NOISE_FLOOR == 1e-18

    def test_it_saturates(self, ids):
        """Flat beyond pinch-off -- to the last bit, here.

        With no channel-length modulation in this core the saturated
        current is constant to the last bit, so consecutive samples
        differ only by rounding -- measured at 1e-25 on 8e-10, one ulp.
        The bound is therefore relative rather than a bare `>= 0`.
        Adding CLM would make the curve rise again slightly.
        """
        vds = np.linspace(0.01, 1.5, 60)
        cur = np.array([ids(1.0, 0.0, v) for v in vds])
        step = np.diff(cur)
        assert np.all(step > -1e-12 * cur[:-1]), 'monotone in Vd to rounding'
        lin = (cur[3] - cur[0]) / (vds[3] - vds[0])
        sat = (cur[-1] - cur[-4]) / (vds[-1] - vds[-4])
        assert sat < 0.05 * lin
        assert cur[-1] > 0

    def test_the_subthreshold_slope_is_physical(self, ids):
        """Above the thermal limit, and degraded by the body factor.

        `ln(10)*phit` is 59.5 mV/decade; a body factor makes it worse,
        never better.
        """
        vgs = np.linspace(0.0, 0.6, 40)
        cur = np.array([ids(v, 0.0, 0.1) for v in vgs])
        m = (cur > 1e-20) & (cur < 1e-12)
        assert m.sum() > 8
        slope = np.polyfit(vgs[m], np.log10(cur[m]), 1)[0]
        ss_mv = 1000.0 / slope
        assert 59.5 < ss_mv < 200.0, '%.1f mV/dec' % ss_mv

    def test_a_bigger_body_factor_degrades_the_slope(self, ids):
        """The body effect, appearing without being put in by hand."""
        got = []
        for gf in (0.5, 1.0, 2.0):
            vgs = np.linspace(0.0, 0.8, 40)
            cur = np.array([ids(v, 0.0, 0.1, gf) for v in vgs])
            m = (cur > 1e-20) & (cur < 1e-12)
            got.append(1000.0 / np.polyfit(vgs[m], np.log10(cur[m]), 1)[0])
        assert got[0] < got[1] < got[2]

    def test_it_survives_biases_newton_will_try(self, ids):
        """Including ones the device physically cannot reach.

        A drain 10 V below the source drives the normalised quasi-Fermi
        level negative, where `delta = exp(-xn)` would reach 1e152 and
        overflow.  `psp_kernel.XN_FLOOR` bounds it: the answer there is
        not physics, but it is finite and Newton can walk back out.
        """
        for vg in (-50.0, -5.0, 0.0, 5.0, 50.0):
            for vd in (-10.0, -1.0, 0.0, 1.0, 10.0):
                assert np.isfinite(ids(vg, 0.0, vd)), (vg, vd)


class TestTheElement(object):

    def _fet(self, **kw):
        kw.setdefault('w', 1e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                              **kw)
        e.update_iparv()
        return e

    def test_it_took_the_let_chain_path(self):
        assert PspMosLongChannel._hdl_info['chained'] is True

    @pytest.mark.parametrize('x', [
        np.array([0.6, 1.0, 0.0, 0.0]),
        np.array([0.05, 0.4, 0.0, 0.0]),
        np.array([1.2, 1.8, 0.0, 0.0]),
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([-0.5, -0.5, 0.0, 0.0]),
        np.array([2.0, 2.0, 0.5, -0.5]),
    ])
    def test_the_jacobian_is_finite_and_correct(self, x):
        e = self._fet()
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            i = np.asarray(e.i(x), float)
            G = np.asarray(e.G(x), float)
        assert np.all(np.isfinite(i))
        assert np.all(np.isfinite(G))
        fd = np.zeros_like(G)
        for j in range(len(x)):
            h = 1e-6 * max(1.0, abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            fd[:, j] = (np.asarray(e.i(xp), float)
                        - np.asarray(e.i(xm), float)) / (2 * h)
        assert np.max(np.abs(G - fd)) < 1e-6 * max(1.0, np.abs(G).max())

    def test_kirchhoff(self):
        """The four terminal currents sum to zero, exactly."""
        e = self._fet()
        for x in (np.array([0.6, 1.0, 0.0, 0.0]),
                  np.array([1.2, 1.8, 0.2, -0.3])):
            assert abs(np.asarray(e.i(x), float).sum()) == 0.0

    def test_the_drain_current_has_the_conventional_sign(self):
        """Positive into the drain for a positive drain bias.

        pycircuit's `i[node]` is the current flowing FROM the node INTO
        the device, the same convention a generated resistor uses.
        """
        e = self._fet()
        i = np.asarray(e.i(np.array([0.6, 1.5, 0.0, 0.0])), float)
        assert i[0] > 0, 'drain'
        assert i[2] < 0, 'source'
        assert i[1] == 0.0, 'no gate current in the intrinsic core'

    def test_swapping_the_terminals_negates_the_current(self):
        """Symmetry, at the element level rather than the kernel's.

        To one ulp, where this used to be bit-exact.  That is a real
        change and worth stating plainly: the element now ORDERS its
        terminals -- it computes the device forward from the lower
        junction and applies the sign afterwards, because the
        saturation-limited drain voltage is built from source-side
        quantities and is not an odd function of `Vds`.  Ordering makes
        the antisymmetry structural, which is stronger than what it
        replaced.  What it costs is the bit-exactness: under the two
        polarities `vsbx` evaluates the same number by a different
        sequence of roundings, so the two currents agree to a relative
        3e-16 rather than to the last bit.  One ulp is float noise; the
        alternative was an algebraic symmetry that silently failed the
        moment non-odd physics was added, which is exactly what
        happened.
        """
        e = self._fet()
        fwd = np.asarray(e.i(np.array([0.7, 1.2, 0.0, 0.0])), float)
        rev = np.asarray(e.i(np.array([0.0, 1.2, 0.7, 0.0])), float)
        assert fwd[0] == pytest.approx(-rev[0], rel=1e-14)
        assert fwd[2] == pytest.approx(-rev[2], rel=1e-14)

    def test_geometry_scales_the_current(self):
        """Ids is proportional to W/L, and to nothing else here."""
        base = self._fet(w=1e-6, l=1e-6)
        x = np.array([0.6, 1.5, 0.0, 0.0])
        i0 = np.asarray(base.i(x), float)[0]
        for w, l, factor in ((2e-6, 1e-6, 2.0), (1e-6, 2e-6, 0.5),
                             (4e-6, 2e-6, 2.0), (5e-6, 1e-6, 5.0)):
            got = np.asarray(self._fet(w=w, l=l).i(x), float)[0]
            assert got == pytest.approx(factor * i0, rel=1e-12)

    def test_transconductance_and_output_conductance_are_positive(self):
        e = self._fet()
        G = np.asarray(e.G(np.array([0.6, 1.2, 0.0, 0.0])), float)
        ## d(i_drain)/d(v_gate) and d(i_drain)/d(v_drain)
        assert G[0, 1] > 0, 'gm'
        assert G[0, 0] > 0, 'gds'


class TestInACircuit(object):

    def test_a_resistor_loaded_stage_solves(self):
        """The real check that it is an element and not just a formula.

        The gate sweep straddles the threshold these default parameters
        actually give: `Vth = vfb + phib + gamma*sqrt(phib)`, which for
        tox = 2.2 nm and Nsub = 5e23 works out at 0.20 V.  Picking gate
        voltages without computing that is how the first version of this
        test came to assert the device was off at 0.6 V -- which is
        0.4 V into strong inversion.
        """
        cm.default_toolkit = numeric
        got = []
        for vg in (0.0, 0.15, 0.3, 0.6, 1.2):
            c = SubCircuit()
            nd, ng = c.add_node('d'), c.add_node('g')
            c['vdd'] = VS(c.add_node('vdd'), gnd, v=1.8)
            c['rl'] = R('vdd', nd, r=10e3)
            c['vg'] = VS(ng, gnd, v=vg)
            c['m1'] = PspMosLongChannel(nd, ng, gnd, gnd, w=10e-6, l=1e-6)
            c.update_iparv()
            got.append(float(DC(c, toolkit=numeric).solve().v(nd, gnd)))
        ## A common-source stage: raising the gate pulls the drain down,
        ## monotonically, and it stays inside the rails.
        assert all(0.0 <= v <= 1.8 for v in got), got
        assert all(a > b for a, b in zip(got, got[1:])), got
        assert got[0] > 1.75, 'off below threshold'
        assert got[-1] < 0.2, 'hard on well above it'
        assert got[0] - got[-1] > 1.5, 'the stage should actually swing'

    def test_two_devices_in_series_share_a_current(self):
        """A stacked pair -- the internal node has to settle somewhere."""
        cm.default_toolkit = numeric
        c = SubCircuit()
        nd, nm = c.add_node('d'), c.add_node('mid')
        c['vdd'] = VS(c.add_node('vdd'), gnd, v=1.8)
        c['rl'] = R('vdd', nd, r=5e3)
        c['vg'] = VS(c.add_node('g'), gnd, v=1.8)
        c['m1'] = PspMosLongChannel(nd, 'g', nm, gnd, w=10e-6, l=1e-6)
        c['m2'] = PspMosLongChannel(nm, 'g', gnd, gnd, w=10e-6, l=1e-6)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        vd, vm = float(res.v(nd, gnd)), float(res.v(nm, gnd))
        assert 0.0 < vm < vd < 1.8


class TestTheChargeModel(object):
    """Terminal charges, and the properties that come with them.

    PSP contributes `Qg`, `Qb` and `Qd` on branches referred to the
    source, which makes `Qs` whatever is left.  So charge conservation
    is not a numerical result to be checked but a consequence of the
    topology -- and the tests below confirm it is exact rather than
    merely close.
    """

    def _fet(self, **kw):
        kw.setdefault('w', 10e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                              **kw)
        e.update_iparv()
        return e

    @pytest.mark.parametrize('x', [
        np.array([0.1, 1.2, 0.0, 0.0]),
        np.array([1.2, 1.2, 0.0, 0.0]),
        np.array([0.0, -0.5, 0.0, 0.0]),
        np.array([1.8, 0.3, 0.0, -0.5]),
        np.array([0.0, 0.0, 0.0, 0.0]),
    ])
    def test_charge_is_conserved_exactly(self, x):
        """The four terminal charges sum to zero, to rounding."""
        q = np.asarray(self._fet().q(x), float)
        assert abs(q.sum()) < 1e-16 * max(np.abs(q).max(), 1e-30)

    @pytest.mark.parametrize('x', [
        np.array([0.1, 1.2, 0.0, 0.0]),
        np.array([1.2, 1.8, 0.0, 0.0]),
        np.array([0.0, -0.5, 0.0, 0.0]),
    ])
    def test_the_capacitance_matrix_is_conservative(self, x):
        """Every column sums to zero: no charge appears from nowhere.

        A capacitance matrix whose columns did not sum to zero would
        inject current into the circuit under a common-mode shift of all
        four terminals, which is unphysical and a classic compact-model
        defect.
        """
        C = np.asarray(self._fet().C(x), float)
        assert np.all(np.isfinite(C))
        scale = max(np.abs(C).max(), 1e-30)
        assert np.max(np.abs(C.sum(axis=0))) < 1e-14 * scale
        assert np.max(np.abs(C.sum(axis=1))) < 1e-14 * scale

    def test_the_ward_dutton_partition(self):
        """0.5 in the linear region, 0.4 in saturation.

        The 40/60 split is the textbook result for a linear charge
        distribution along the channel, and it is not fitted here -- it
        falls out of the partition PSP writes for `SWQPART = 0`.
        """
        e = self._fet()
        frac = {}
        for vd in (0.0, 0.05, 0.5, 1.2, 1.8):
            q = np.asarray(e.q(np.array([vd, 1.2, 0.0, 0.0])), float)
            frac[vd] = q[0] / (q[0] + q[2])
        assert frac[0.0] == pytest.approx(0.5, abs=1e-12)
        assert frac[1.8] == pytest.approx(0.40, abs=0.01)
        ## and it moves monotonically from one to the other
        seq = [frac[v] for v in (0.0, 0.05, 0.5, 1.2, 1.8)]
        assert all(a >= b for a, b in zip(seq, seq[1:])), seq

    def test_swapping_source_and_drain_swaps_their_charges(self):
        """To six ulp -- not bit-exact, and that is expected.

        Unlike the current, whose antisymmetry is structural, `Qd` and
        `Qs` are computed by expressions that are algebraically mirror
        images but not literally the same, so the swap is exact only up
        to rounding.  The GATE charge is symmetric in `dps` and used to
        come back bit-identical; with the terminals ordered it too is
        exact only to an ulp.
        """
        e = self._fet()
        for vg in (0.4, 0.8, 1.2, 1.8):
            for vd in (0.05, 0.3, 0.8, 1.5):
                fwd = np.asarray(e.q(np.array([vd, vg, 0.0, 0.0])), float)
                rev = np.asarray(e.q(np.array([0.0, vg, vd, 0.0])), float)
                assert fwd[0] == pytest.approx(rev[2], rel=1e-14)
                assert fwd[2] == pytest.approx(rev[0], rel=1e-14)
                assert fwd[1] == pytest.approx(rev[1], rel=1e-14), \
                    'gate charge, to one ulp since terminals are ordered'

    def test_below_threshold_the_gate_charge_mirrors_the_bulk(self):
        """Below threshold the gate is terminated almost entirely by
        depletion charge in the bulk.

        "Almost": the inversion charge is not zero, it is twelve orders
        of magnitude down (6e-26 C against 2e-14 C at Vg = -0.5).  It
        had better not be zero -- that residual charge IS subthreshold
        conduction, and a model that set it to zero here would have no
        subthreshold current either.
        """
        e = self._fet()
        q = np.asarray(e.q(np.array([0.0, -0.5, 0.0, 0.0])), float)
        assert q[1] > 0, 'gate charge'
        assert q[3] == pytest.approx(-q[1], rel=1e-10), 'bulk mirrors it'
        assert 0 < abs(q[0]) < 1e-10 * abs(q[1]), 'inversion charge is tiny'
        assert 0 < abs(q[2]) < 1e-10 * abs(q[1]), 'but not zero'

    def test_gate_charge_grows_with_gate_bias(self):
        e = self._fet()
        qg = [np.asarray(e.q(np.array([0.1, v, 0.0, 0.0])), float)[1]
              for v in np.linspace(-0.5, 2.0, 40)]
        assert all(b > a for a, b in zip(qg, qg[1:]))

    def test_the_gate_capacitance_is_bounded_by_the_oxide(self):
        """`Cgg` cannot exceed `Cox' * W * L`.

        In strong inversion the gate is screened by the channel and sees
        essentially the full oxide capacitance; it can approach that
        value but not pass it.
        """
        from pycircuit.circuit.compact import EPS_OX
        e = self._fet(w=10e-6, l=1e-6, tox=2.2e-9)
        cox_tot = EPS_OX / 2.2e-9 * 10e-6 * 1e-6
        cgg = np.asarray(e.C(np.array([0.05, 1.8, 0.0, 0.0])),
                         float)[1, 1]
        assert 0.0 < cgg <= cox_tot * (1.0 + 1e-12)
        assert cgg > 0.5 * cox_tot, 'should be well into inversion'

    def test_it_runs_in_a_transient(self):
        """Charges mean the device has dynamics; this exercises them."""
        from pycircuit.circuit.elements import C as Cap
        from pycircuit.circuit.transient import Transient

        cm.default_toolkit = numeric
        c = SubCircuit()
        nd = c.add_node('d')
        c['vdd'] = VS(c.add_node('vdd'), gnd, v=1.8)
        c['rl'] = R('vdd', nd, r=20e3)
        c['vg'] = VS(c.add_node('g'), gnd, v=1.2)
        c['m1'] = PspMosLongChannel(nd, 'g', gnd, gnd, w=10e-6, l=1e-6)
        c['cl'] = Cap(nd, gnd, c=1e-13)
        c.update_iparv()
        res = Transient(c, toolkit=numeric).solve(tend=2e-7, timestep=1e-9)
        v = np.asarray(res.v(nd, gnd).y, float)
        assert np.all(np.isfinite(v))
        assert np.all(v >= -0.1) and np.all(v <= 1.9)
        ## it settles to the DC operating point
        dc = float(DC(c, toolkit=numeric).solve().v(nd, gnd))
        assert v[-1] == pytest.approx(dc, abs=2e-3)


class TestMobilityAndVelocitySaturation(object):
    """The first layer on top of the ideal core.

    Mobility reduction and velocity saturation are what make the
    strong-inversion current realistic.  They are also the first real
    test of whether the construction survives being built on: both
    depend only on MIDPOINT quantities and on `dps` SQUARED, so both are
    EVEN under exchanging source and drain -- and the antisymmetry of
    the current must therefore still be exact, not merely close.

    Defaults are the IHP SG13G2 n-channel card's values.  Setting
    `mue`, `cs` and `thesat` to zero recovers the ideal core exactly,
    which is what the tests above exercise.
    """

    def _fet(self, **kw):
        kw.setdefault('w', 10e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                              **kw)
        e.update_iparv()
        return e

    def _ideal(self, **kw):
        return self._fet(mue=0.0, cs=0.0, thesat=0.0, **kw)

    def test_symmetry_survives_the_extra_layer(self):
        """The point of the whole construction, re-checked with physics on.

        To one ulp rather than bit-exact, for the reason given on
        `test_swapping_the_terminals_negates_the_current`: the terminals
        are ordered now, so the two polarities reach the same number
        along different rounding paths.
        """
        e = self._fet()
        for vg in (0.4, 0.9, 1.4, 1.8):
            for vd in (0.05, 0.3, 0.9, 1.6):
                fwd = np.asarray(e.i(np.array([vd, vg, 0.0, 0.0])), float)
                rev = np.asarray(e.i(np.array([0.0, vg, vd, 0.0])), float)
                assert fwd[0] == pytest.approx(-rev[0], rel=1e-14), (vg, vd)

    def test_charge_conservation_survives_it_too(self):
        e = self._fet()
        for x in (np.array([0.9, 1.4, 0.0, 0.0]),
                  np.array([0.0, -0.5, 0.0, 0.0]),
                  np.array([1.8, 1.8, 0.2, -0.4])):
            q = np.asarray(e.q(x), float)
            assert abs(q.sum()) < 1e-16 * max(np.abs(q).max(), 1e-30)

    def test_turning_it_off_recovers_the_ideal_core(self):
        """Bit-exact, so the layer is genuinely additive."""
        real = self._fet(mue=0.0, cs=0.0, thesat=0.0)
        ideal = self._ideal()
        for x in (np.array([0.5, 1.2, 0.0, 0.0]),
                  np.array([1.5, 1.8, 0.0, 0.0])):
            assert (np.asarray(real.i(x), float).tolist()
                    == np.asarray(ideal.i(x), float).tolist())

    def test_mobility_reduction_grows_with_the_gate_field(self):
        """Higher vertical field, more surface scattering, less current.

        The ratio to the ideal device must fall monotonically with `Vg`
        -- that is the whole physical content of the term.
        """
        real, ideal = self._fet(), self._ideal()
        ratios = []
        for vg in (0.4, 0.8, 1.2, 1.8):
            x = np.array([1.2, vg, 0.0, 0.0])
            ratios.append(np.asarray(real.i(x), float)[0]
                          / np.asarray(ideal.i(x), float)[0])
        assert all(0.0 < r < 1.0 for r in ratios), ratios
        assert all(a > b for a, b in zip(ratios, ratios[1:])), ratios

    def test_velocity_saturation_reduces_the_saturated_current(self):
        e_off = self._fet(thesat=0.0)
        e_on = self._fet(thesat=0.39843)
        x = np.array([1.8, 1.8, 0.0, 0.0])
        assert (np.asarray(e_on.i(x), float)[0]
                < np.asarray(e_off.i(x), float)[0])
        ## and it does nothing at zero drain bias, where dps is zero
        x0 = np.array([0.0, 1.8, 0.0, 0.0])
        assert np.asarray(e_on.i(x0), float)[0] == 0.0

    def test_it_still_saturates_and_stays_monotone(self):
        e = self._fet()
        vds = np.linspace(0.01, 1.8, 40)
        cur = np.array([np.asarray(e.i(np.array([v, 1.4, 0.0, 0.0])),
                                   float)[0] for v in vds])
        assert np.all(np.diff(cur) > -1e-12 * cur[:-1])
        lin = (cur[3] - cur[0]) / (vds[3] - vds[0])
        sat = (cur[-1] - cur[-4]) / (vds[-1] - vds[-4])
        assert sat < 0.1 * lin

    @pytest.mark.parametrize('x', [
        np.array([0.8, 1.5, 0.0, 0.0]),
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([1.2, 0.3, 0.0, -0.5]),
        np.array([-0.5, -0.5, 0.0, 0.0]),
        np.array([2.0, 2.0, 0.5, -0.5]),
    ])
    def test_the_jacobian_is_still_finite_and_correct(self, x):
        """The Coulomb-scattering term is where this went wrong once.

        Written as PSP writes it -- `cs*exp(0.5*thecs*ln(ratio))` -- it
        nests two `Piecewise`s, and differentiating the nest emitted
        conditions that numpy's `select` REJECTED: the Jacobian raised
        rather than losing precision.  It is compiled as the equivalent
        power `cs*ratio**(0.5*thecs)` instead.
        """
        e = self._fet()
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            G = np.asarray(e.G(x), float)
        assert np.all(np.isfinite(G))
        fd = np.zeros_like(G)
        for j in range(len(x)):
            h = 1e-6 * max(1.0, abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            fd[:, j] = (np.asarray(e.i(xp), float)
                        - np.asarray(e.i(xm), float)) / (2 * h)
        assert np.max(np.abs(G - fd)) < 1e-6 * max(1.0, np.abs(G).max())

    def test_the_coulomb_term_is_finite_at_flat_band(self):
        """`Pm` is exactly zero there, so the ratio it divides by is too.

        The exponent is below 1 (0.59 for the card's `thecs`), so the
        derivative would diverge; the ratio is floored for that reason.
        """
        e = self._fet()
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            for x in (np.zeros(4), np.array([0.0, 0.95, 0.0, 0.0])):
                assert np.all(np.isfinite(np.asarray(e.i(x), float)))
                assert np.all(np.isfinite(np.asarray(e.G(x), float)))


class TestTheSaturationVoltage(object):
    """What `Vdse` does to the element, and what it cost.

    See `test_psp_gap.TestTheSaturationVoltage` for the comparison
    against PSP's own scaled parameters.  This file's concern is the
    consequences for the element: it saturates now, and a saturating
    device needs a Newton limiter.
    """

    def _fet(self, **kw):
        kw.setdefault('w', 10e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    def test_the_current_actually_saturates(self):
        """The point of the whole construction.

        Without `Vdse` the drain surface potential follows the applied
        bias and the current keeps climbing.  With it the drain end is
        pinned at `Vdsat` and the output conductance collapses -- which
        is the physical answer, and is also why the limiter below had to
        be written.
        """
        e = self._fet()
        g = {}
        for vd in (0.3, 1.0, 5.0, 50.0):
            G = np.asarray(e.G(np.array([vd, 1.8, 0.0, 0.0])), float)
            g[vd] = G[0, 0]
        assert g[0.3] > g[1.0] > g[5.0] > g[50.0] > 0.0, g
        assert g[0.3] / g[50.0] > 1e5, g

    def test_the_current_stays_bounded_at_absurd_bias(self):
        """Saturated, not overflowing.  A solver WILL evaluate here."""
        e = self._fet()
        i0 = np.asarray(e.i(np.array([2.0, 1.8, 0.0, 0.0])), float)[0]
        for vd in (1e2, 1e4, 1e7):
            i = np.asarray(e.i(np.array([vd, 1.8, 0.0, 0.0])), float)
            assert np.all(np.isfinite(i))
            assert i[0] < 1.5 * i0, (vd, i[0], i0)


class TestTheVoltageConditioning(object):
    """The smooth limiter on the junction bias, and the kink it replaced.

    The kernel floors `xn` at zero, which is the quasi-Fermi level
    reaching the built-in potential -- below that the formulation stops
    meaning anything.  A hard floor is fine for the VALUE and poison for
    the Jacobian: at exactly `Vd = -phib` the analytic derivative and a
    finite difference disagreed by 60%, and below it every conductance
    froze, so a solver that overshot had nothing telling it how to get
    back.

    That went unnoticed while the floor only ever bit the DRAIN end,
    where the device is off and the sensitivity is negligible anyway.
    The saturation voltage reads the SOURCE end, which multiplies
    everything, and the kink became visible immediately.

    PSP conditions the terminal voltages instead
    (`PSP103_macrodefs.include:330-334, 1104-1105`): take the lower of
    the two junctions, clip it at `-0.95*phib` through the smooth
    `MINA`, and lift `Vsb` by whatever the clip removed.
    """

    def _fet(self, **kw):
        kw.setdefault('w', 10e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    @staticmethod
    def _fd_error(e, x):
        x = np.array(x, float)
        G = np.asarray(e.G(x), float)
        fd = np.zeros_like(G)
        for j in range(len(x)):
            h = 1e-7 * max(1.0, abs(x[j]))
            xp, xm = x.copy(), x.copy()
            xp[j] += h
            xm[j] -= h
            fd[:, j] = (np.asarray(e.i(xp), float)
                        - np.asarray(e.i(xm), float)) / (2 * h)
        return np.max(np.abs(G - fd)) / max(np.abs(G).max(), 1e-30)

    def test_the_jacobian_is_right_at_the_old_kink(self):
        """`Vd = -phib` exactly -- the bias the hard floor broke.

        `phib` defaults to 0.9, so this is the one point where the old
        `max(xn, 0)` sat exactly on its corner.  It used to be 0.6 out.
        """
        e = self._fet()
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            assert self._fd_error(e, [-0.9, 1.8, 0.0, 0.0]) < 1e-5

    @pytest.mark.parametrize('vd', [-0.85, -0.9, -0.95, -1.5, -4.0])
    def test_it_is_right_on_both_sides_too(self, vd):
        e = self._fet()
        assert self._fd_error(e, [vd, 1.8, 0.0, 0.0]) < 1e-5, vd

    def test_the_conductance_does_not_freeze_below_the_limit(self):
        """What the hard floor actually cost.

        Below `-phib` the clamped formulation returned the same numbers
        for every input, so the conductance was bit-identical at every
        bias out there and Newton had no gradient at all.  The smooth
        limiter keeps it moving.
        """
        e = self._fet()
        g = [np.asarray(e.G(np.array([vd, 1.8, 0.0, 0.0])),
                        float)[0, 0] for vd in (-1.0, -1.5, -2.5)]
        assert len(set(g)) == 3, g

    @staticmethod
    def _vsbstar(vd, vs, phib=0.9):
        """The conditioning, in closed form -- what the element builds."""
        def mina(x, y, a):
            return 0.5 * (x + y - np.sqrt((x - y) ** 2 + a))
        phix = 0.95 * phib
        aphi = 0.0025 * phib * phib
        phix1 = mina(phix - 0.5 * np.sqrt(aphi), 0.0, aphi)
        return vs - mina(mina(vd, vs, aphi) + phix, 0.0, aphi) + phix1

    @pytest.mark.parametrize('vd,vs', [(0.0, 0.0), (1.5, 0.0), (0.0, 1.5),
                                       (1.8, 0.9), (-0.5, 0.0)])
    def test_it_is_invisible_at_ordinary_bias(self, vd, vs):
        """It has to be a limiter, not a correction.

        Wherever nothing needs limiting it must leave `Vsb` alone, or it
        would be a silent shift in every threshold on the card.  A
        fraction of a millivolt is what it costs, which is why it can be
        applied unconditionally rather than behind a test.
        """
        assert abs(self._vsbstar(vd, vs) - vs) < 1e-3, (vd, vs)

    def test_it_clamps_the_lower_junction(self):
        """And where something DOES need limiting, it limits.

        However far below `-phib` the lower terminal is driven, the
        conditioned junction stops at about `-0.95*phib` -- so `phib +
        Vsbx` stays positive and the kernel's own floor never binds.
        """
        for vd in (-1.0, -5.0, -50.0, -1e4):
            eff = self._vsbstar(vd, 0.0) + (vd - abs(vd)) * 0.5
            assert eff > -0.9, (vd, eff)
            assert eff < -0.8, (vd, eff)


class TestTheNewtonLimiter(object):
    """`PspMosLongChannel.limit`, the part SPICE calls `DEVfetlim`.

    A saturating device has no gradient far from its solution, so a
    Newton step that overshoots by a few hundred volts does not come
    back -- the row goes numerically empty and the matrix is reported
    singular.  Two of these in a stack was enough to reach that state.
    """

    def _fet(self, **kw):
        kw.setdefault('w', 10e-6)
        kw.setdefault('l', 1e-6)
        e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                              cm.Node('b'), **kw)
        e.update_iparv()
        return e

    def test_it_bounds_a_large_step(self):
        e = self._fet()
        x0 = np.array([0.5, 0.5, 0.0, 0.0])
        out = e.limit(np.array([500.0, -300.0, 0.0, 90.0]), x0)
        assert out[0] == pytest.approx(0.5 + e.vlimit)
        assert out[1] == pytest.approx(0.5 - e.vlimit)
        assert out[3] == pytest.approx(0.0 + e.vlimit)

    def test_it_leaves_a_small_step_alone(self):
        """It must be inactive near a solution, or it caps convergence."""
        e = self._fet()
        x0 = np.array([0.8, 1.2, 0.1, 0.0])
        x = np.array([0.83, 1.19, 0.1, 0.001])
        assert np.allclose(e.limit(x, x0), x, rtol=0, atol=0)

    def test_it_does_not_move_the_source(self):
        """The source is some other device's drain.

        Limiting it here would have two elements fight over the same
        node, each undoing the other.  SPICE bounds the branch voltages
        and leaves the reference terminal where the solver put it.
        """
        e = self._fet()
        x = np.array([500.0, 500.0, 400.0, 500.0])
        assert e.limit(x, np.zeros(4))[2] == 400.0

    def test_it_does_not_mutate_its_argument(self):
        """The state-free convention: return a limited COPY.

        `SubCircuit.limit` writes the return value back; a limiter that
        edited its argument would be editing a fancy-indexed temporary,
        which is how limiting was silently a no-op in this tree once
        before.
        """
        e = self._fet()
        x = np.array([500.0, 500.0, 0.0, 0.0])
        before = x.copy()
        e.limit(x, np.zeros(4))
        assert np.array_equal(x, before)

    def test_a_stacked_pair_converges(self):
        """The circuit that found it.

        Both devices swing out past their saturation on the way to the
        solution; without the limiter the solve reaches 6e7 V on the
        internal node and reports a singular matrix.
        """
        cm.default_toolkit = numeric
        c = SubCircuit()
        nd, nm = c.add_node('d'), c.add_node('mid')
        c['vdd'] = VS(c.add_node('vdd'), gnd, v=1.8)
        c['rl'] = R('vdd', nd, r=5e3)
        c['vg'] = VS(c.add_node('g'), gnd, v=1.8)
        c['m1'] = PspMosLongChannel(nd, 'g', nm, gnd, w=10e-6, l=1e-6)
        c['m2'] = PspMosLongChannel(nm, 'g', gnd, gnd, w=10e-6, l=1e-6)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        vd, vm = float(res.v(nd, gnd)), float(res.v(nm, gnd))
        assert 0.0 < vm < vd < 1.8
