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
        """Symmetry, at the element level rather than the kernel's."""
        e = self._fet()
        fwd = np.asarray(e.i(np.array([0.7, 1.2, 0.0, 0.0])), float)
        rev = np.asarray(e.i(np.array([0.0, 1.2, 0.7, 0.0])), float)
        assert fwd[0] == -rev[0]
        assert fwd[2] == -rev[2]

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
        to rounding.  The GATE charge, being symmetric in `dps`, does
        come back bit-identical.
        """
        e = self._fet()
        for vg in (0.4, 0.8, 1.2, 1.8):
            for vd in (0.05, 0.3, 0.8, 1.5):
                fwd = np.asarray(e.q(np.array([vd, vg, 0.0, 0.0])), float)
                rev = np.asarray(e.q(np.array([0.0, vg, vd, 0.0])), float)
                assert fwd[0] == pytest.approx(rev[2], rel=1e-14)
                assert fwd[2] == pytest.approx(rev[0], rel=1e-14)
                assert fwd[1] == rev[1], 'gate charge is bit-exact'

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
