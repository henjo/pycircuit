"""The charge-conserving MOS gate charge (`elements_hdl._mos_gate_charge`,
the SPICE benchmark plan's stage 6): its capacitances against SPICE's
Meyer model where the docstring claims them, its deviations where it
names them, and conservation.

The reference is a numpy transcription of ngspice's `DEVqmeyer`
(`devsup.c`), doubled -- SPICE returns half of each capacitance and adds
the previous step's half -- with its `MAGIC_VDS` 25 mV floor on `vdsat`.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.tests.test_hdl_cbackend import needs_cc

EPSOX = 3.9 * 8.854187817e-12
VTO, PHI, TOX, W, L = 0.7, 0.7, 2e-8, 10e-6, 2e-6
COX = EPSOX / TOX * W * L
D, G, S, B = range(4)


def meyer(vgs, vgd, von, vdsat, phi, cox):
    """`DEVqmeyer` doubled: ``(cgs, cgd, cgb)``, forward mode (``vds >= 0``)."""
    vgst = vgs - von
    vdsat = max(vdsat, 0.025)
    if vgst <= -phi:
        return 0.0, 0.0, cox
    if vgst <= -phi / 2:
        return 0.0, 0.0, -vgst * cox / phi
    vds = vgs - vgd
    if vgst <= 0:
        cgb = -vgst * cox / phi
        cgs = 2 * (vgst * cox / (1.5 * phi) + cox / 3)
        if vds >= vdsat:
            return cgs, 0.0, cgb
        d2 = (2 * vdsat - vds) ** 2
        return cgs * (1 - (vdsat - vds) ** 2 / d2), cgs * (1 - vdsat ** 2 / d2), cgb
    if vdsat <= vds:
        return 2 * cox / 3, 0.0, 0.0
    d2 = (2 * vdsat - vds) ** 2
    return (2 * cox / 3 * (1 - (vdsat - vds) ** 2 / d2),
            2 * cox / 3 * (1 - vdsat ** 2 / d2), 0.0)


def _device(cls, **kw):
    card = {"vto": VTO, "kp": 5e-5, "tox": TOX, "w": W, "l": L, "phi": PHI}
    card.update(kw)
    e = cls(*[Node(n) for n in 'dgsb'], **card)
    e.update_iparv()
    return e


def _C(e, vd, vg, vs, vb):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(e.C(np.array([vd, vg, vs, vb], float), defaultepar), float)


def _gate_row(C):
    """Meyer's three, from the gate row: ``C_gx = -dQg/dVx``."""
    return -C[G, S], -C[G, D], -C[G, B]


def _sliver(vgst, vds):
    """Where the two differ by design below 25 mV of `vds`: Meyer's 25 mV
    floor on `vdsat` (to 25 mV above threshold) and this model's rounded
    `max(us, ud)` (every gate drive below that)."""
    return vds < 0.025 and vgst < 0.025


def test_level1_gate_capacitances_are_meyers_where_claimed():
    """Level 1, bulk at the source: the gate's TOTAL capacitance is
    Meyer's to rounding in accumulation, depletion, the transition and
    inversion; source/drain/bulk split exactly in inversion; below
    threshold the charge follows `vgs`, so Meyer's gate-bulk capacitance
    reads against the source (the total and `Cgd` unchanged); and in
    Meyer's 25 mV sliver his total is at most 1.5x this one."""
    e = _device(eh.MosLevel1GateChargeHdl)
    worst_ratio, n_sliver, n_exact = 1.0, 0, 0
    for vgs in np.linspace(-1.5, 3.0, 181):
        for vds in (0.0, 0.01, 0.02, 0.03, 0.2, 1.0, 2.5):
            vgst = vgs - VTO
            cgs, cgd, cgb = _gate_row(_C(e, vds, vgs, 0.0, 0.0))
            mgs, mgd, mgb = meyer(vgs, vgs - vds, VTO, max(vgst, 0.0), PHI, COX)
            if _sliver(vgst, vds):
                n_sliver += 1
                tot, mtot = cgs + cgd + cgb, mgs + mgd + mgb
                worst_ratio = max(worst_ratio, mtot / tot, tot / mtot)
                continue
            n_exact += 1
            assert cgs + cgd + cgb == pytest.approx(mgs + mgd + mgb, rel=0, abs=1e-12 * COX)
            if vgst > 0:
                assert (cgs, cgd, cgb) == pytest.approx((mgs, mgd, mgb), rel=0,
                                                        abs=1e-12 * COX), (vgs, vds)
            else:
                assert cgd == pytest.approx(mgd, rel=0, abs=1e-12 * COX)
                assert cgs + cgb == pytest.approx(mgs + mgb, rel=0, abs=1e-12 * COX)
    assert n_sliver > 20 and n_exact > 900
    assert worst_ratio <= 1.5 + 1e-9, worst_ratio
    print(f'Meyer sliver: {n_sliver} points, his total at most {worst_ratio:.3f}x')


def _level3_vdsat(vgst, theta, vmax, u0=600.0):
    """`mos3load.c`'s `vdsat` at `gamma = delta = eta = nfs = 0`."""
    us = u0 * 1e-4 / (1.0 + theta * vgst)
    if vmax <= 0:
        return vgst
    vdsc = (L * vmax) / us
    return vgst + vdsc - np.sqrt(vgst * vgst + vdsc * vdsc)


@pytest.mark.parametrize('card', [{}, {'vmax': 1e5, 'theta': 0.1}], ids=['plain', 'vmax'])
def test_level3_is_meyers_at_vds_zero_and_in_saturation_and_bounded_between(card):
    """Level 3 against Meyer fed level 3's own `vdsat`: exact at `vds = 0`
    and wherever both saturate (`vds >= vgs - von`); between, where
    `vmax` and `theta` make level 3's `vdsat` smaller than `vgs - von`,
    the named triode approximation -- printed, and bounded by the 2/3 Cox
    the channel can move."""
    e = _device(eh.MosLevel3GateChargeHdl, **card)
    worst = 0.0
    for vgs in np.linspace(0.75, 3.0, 46):
        vgst = vgs - VTO
        vdsat = _level3_vdsat(vgst, card.get('theta', 0.0), card.get('vmax', 0.0))
        for vds in (0.0, 0.3 * vdsat, vdsat, 0.5 * (vdsat + vgst), vgst, 1.5 * vgst + 0.1):
            got = _gate_row(_C(e, vds, vgs, 0.0, 0.0))
            want = meyer(vgs, vgs - vds, VTO, vdsat, PHI, COX)
            if vgst < 0.025:
                continue
            if vds == 0.0 or vds >= vgst:
                assert got == pytest.approx(want, rel=0, abs=1e-12 * COX), (vgs, vds)
            else:
                worst = max(worst, max(abs(a - b) for a, b in zip(got, want)) / COX)
    if not card:
        assert worst < 1e-12, worst
    assert worst <= 2.0 / 3.0
    print(f'level 3 {card}: triode deviation at most {worst:.3f} Cox')


@pytest.mark.parametrize('cls', ['MosLevel1GateChargeHdl', 'MosLevel3GateChargeHdl',
                                 'MosLevel1PmosGateChargeHdl', 'MosLevel3PmosGateChargeHdl'])
def test_the_charge_is_conserved_and_sees_only_differences(cls):
    """Every column of the 4x4 capacitance matrix sums to zero (the four
    charges sum to zero: conservation) and every row does (a common shift
    of the four potentials moves no charge), at biases in every region,
    with the body effect on and the bulk away from the source."""
    extra = {'lambd': 0.02} if 'Level1' in cls else {'eta': 0.5, 'nfs': 1e11}
    e = _device(getattr(eh, cls), gamma=0.5, qpart=0.3, **extra)
    rng = np.random.default_rng(20261008)
    sign = -1.0 if 'Pmos' in cls else 1.0
    for _ in range(200):
        vd, vg, vs = rng.uniform(-0.5, 3.0, 3)
        vb = min(vd, vs) - rng.uniform(0.0, 1.5)
        C = _C(e, *(sign * np.array([vd, vg, vs, vb])))
        scale = np.abs(C).max()
        assert np.abs(C.sum(axis=0)).max() <= 1e-12 * scale
        assert np.abs(C.sum(axis=1)).max() <= 1e-12 * scale


def test_the_pmos_is_the_nmos_mirrored():
    """A p-channel's capacitances at `-v` are the n-channel's at `v`."""
    n = _device(eh.MosLevel1GateChargeHdl, gamma=0.5)
    p = _device(eh.MosLevel1PmosGateChargeHdl, gamma=0.5)
    rng = np.random.default_rng(7)
    for _ in range(50):
        v = rng.uniform(-1.0, 3.0, 4)
        np.testing.assert_allclose(_C(p, *(-v)), _C(n, *v), rtol=0, atol=1e-12 * COX)


def _loop(e, n):
    """The gate current integrated round a closed loop of the drain and
    gate potentials (source and bulk at 0) crossing every region, by the
    midpoint rule on `n` intervals: ``(ours, Meyer's)``, in Cox."""
    t = np.linspace(0.0, 1.0, n + 1)
    V = [1.0 + np.sin(2 * np.pi * t + 1.3), 1.2 + 1.8 * np.sin(2 * np.pi * t), 0 * t, 0 * t]
    ours = meyers = 0.0
    for k in range(n):
        mid = [0.5 * (a[k] + a[k + 1]) for a in V]
        dv = np.array([a[k + 1] - a[k] for a in V])
        ours += _C(e, *mid)[G] @ dv
        vgs, vgd = mid[G] - mid[S], mid[G] - mid[D]
        if vgs >= vgd:
            mgs, mgd, mgb = meyer(vgs, vgd, VTO, max(vgs - VTO, 0.0), PHI, COX)
        else:
            mgd, mgs, mgb = meyer(vgd, vgs, VTO, max(vgd - VTO, 0.0), PHI, COX)
        meyers += mgs * (dv[G] - dv[S]) + mgd * (dv[G] - dv[D]) + mgb * (dv[G] - dv[B])
    return ours / COX, meyers / COX


def test_a_closed_loop_returns_no_charge_where_meyer_does():
    """Drive the terminals round a closed loop crossing every region and
    integrate the gate current, `sum C_gx dV_x`: this model's is the
    quadrature's error alone -- it shrinks with the step (the capacitances
    jump at threshold, so first order) -- where Meyer's converges to a net
    charge, the reason for the stage."""
    e = _device(eh.MosLevel1GateChargeHdl)
    coarse, fine = _loop(e, 1000), _loop(e, 16000)
    assert abs(fine[0]) <= abs(coarse[0]) / 5, (coarse, fine)
    assert abs(fine[1] - coarse[1]) <= 0.1 * abs(fine[1])      # Meyer's has converged
    assert abs(fine[1]) >= 50 * abs(fine[0]), fine
    print(f'net gate charge round the loop: ours {coarse[0]:.1e} -> {fine[0]:.1e} Cox, '
          f'Meyer {coarse[1]:.2e} -> {fine[1]:.2e} Cox')


@pytest.mark.parametrize('cls', ['MosLevel1GateChargeHdl', 'MosLevel3GateChargeHdl'])
def test_the_two_channel_splits_at_their_defining_points(cls):
    """`qpart = 0`, Ward-Dutton: the drain takes 2/5 of the channel's
    response to the gate in saturation, half at `vds = 0`.  `qpart = 1`,
    following Meyer: none in saturation (`Cdg = 0`, as his `Cgd`) and half
    at `vds = 0` (`Cdg = Cgd = Cox/2`).  The gate row is the same for
    both, and in triode `qpart = 1` keeps the drain's own capacitance
    positive."""
    wd = _device(getattr(eh, cls))
    ml = _device(getattr(eh, cls), qpart=1.0)
    vgs = 2.0
    for e, sat_share in ((wd, 0.4), (ml, 0.0)):
        C = _C(e, 3.0, vgs, 0.0, 0.0)                     # saturated
        cgg = C[G, G]
        assert -C[D, G] == pytest.approx(sat_share * cgg, rel=1e-9, abs=1e-12 * COX)
        C0 = _C(e, 0.0, vgs, 0.0, 0.0)                    # vds = 0
        assert -C0[D, G] == pytest.approx(0.5 * C0[G, G], rel=1e-9)
        assert -C0[D, G] == pytest.approx(COX / 2, rel=1e-9)
    for vds in (0.05, 0.3, 0.8, 1.2):
        a, b = _C(wd, vds, vgs, 0.0, 0.0), _C(ml, vds, vgs, 0.0, 0.0)
        np.testing.assert_allclose(a[G], b[G], rtol=0, atol=1e-12 * COX)
        assert b[D, D] > 0.0


@needs_cc
def test_an_imported_inverter_chain_runs_on_the_c_core(tmp_path):
    """A level-3 deck imports to the gate-charge variants (SPICE adds
    Meyer to every level-3 device), and its transient's Newton is the C
    core's (`newton_c:served`): the variants are chained."""
    from pycircuit.circuit import _paths, circuit, spice_import
    deck = tmp_path / 'chain.cir'
    deck.write_text("""inverter chain
vdd vdd 0 dc 5
vin in 0 pulse 0 5 1n 1n 1n 10n 22n
.model n nmos level=3 vto=0.8 uo=600 tox=3e-8 nsub=1e16 vmax=1.5e5 theta=0.05
.model p pmos level=3 vto=-0.8 uo=250 tox=3e-8 nsub=1e16 vmax=1.5e5 theta=0.05
mp1 a in vdd vdd p w=20u l=2u as=140p ad=140p
mn1 a in 0 0 n w=10u l=2u as=70p ad=70p
mp2 b a vdd vdd p w=20u l=2u as=140p ad=140p
mn2 b a 0 0 n w=10u l=2u as=70p ad=70p
cl b 0 50f
.tran 0.1n 10n
.end
""")
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(str(deck))
    mos = [e for e in imp.circuit.elements.values() if 'Mos' in type(e).__name__]
    assert len(mos) == 4 and all(isinstance(e, (eh.MosLevel3GateChargeHdl,
                                                eh.MosLevel3PmosGateChargeHdl)) for e in mos)
    tr, kw = imp.transient()
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tr.solve(**kw)
    assert _paths.since(before).get('newton_c:served', 0) > 50
    vb = np.asarray(res.v('b').y, float)
    assert vb.min() < 0.5 and vb.max() > 4.5             # it switched both ways
