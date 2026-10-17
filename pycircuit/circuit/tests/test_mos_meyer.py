"""SPICE's Meyer integration as an opt-in MOS gate capacitance
(`elements_hdl._meyer_currents`, `import_netlist(mos_charge='meyer')`;
Andreas, 2026-10-17): Meyer's capacitances as ngspice computes them,
integrated as capacitance-defined currents ``i = C(v) dv/dt`` -- which no
charge has, so they conserve nothing, as SPICE's do not.  Each current runs
through a scratch branch ``y``, potential-contributed ``V(y) - V(minus) = tau
dv/dt``.

The reference is `test_mos_gate_charge.meyer`: ngspice's `DEVqmeyer`
(devsup.c) transcribed and doubled, its `MAGIC_VDS` floor on `vdsat`.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.tests.test_mos_gate_charge import COX, PHI, VTO, _device, meyer

_SCRATCH = ('ygs', 'ygd', 'ygb')


def _state(e, vd, vg, vs, vb):
    """The element's state at its terminals' voltages, each scratch node at
    its branch's minus node (`V(y) - V(minus) = 0`: no time derivative)."""
    lay = [lab for _k, lab, _kind in hdl.x_layout(e)]
    x = np.zeros(len(lay))
    v = {'d': vd, 'g': vg, 's': vs, 'b': vb}
    for k, lab in enumerate(lay):
        if lab in v:
            x[k] = v[lab]
    ## (the scratch nodes reference their branch's minus: the source for gs
    ## and gd -- the internal source/drain collapse with rd = rs = 0 -- the
    ## bulk for gb; n-channel orientation)
    ref = {'ygs': vs, 'ygd': vd, 'ygb': vb}
    for k, lab in enumerate(lay):
        if lab in ref:
            x[k] = ref[lab]
    return x, lay


def _meyer_row(e, vd, vg, vs, vb):
    """The gate's effective capacitances ``(cgs, cgd, cgb)``: the scratch
    nodes eliminated from ``G + jw C`` at the state, the gate row's
    imaginary part over ``w``."""
    x, lay = _state(e, vd, vg, vs, vb)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Gm = np.asarray(e.G(x, defaultepar), float)
        Cm = np.asarray(e.C(x, defaultepar), float)
    ## (every unknown but the four terminals: the scratch nodes and their
    ## branch currents)
    A = [k for k, lab in enumerate(lay) if lab not in ('d', 'g', 's', 'b')]
    K = [k for k in range(len(lay)) if k not in A]
    w = 1e3
    Y = Gm + 1j * w * Cm
    Yr = Y[np.ix_(K, K)] - Y[np.ix_(K, A)] @ np.linalg.solve(Y[np.ix_(A, A)], Y[np.ix_(A, K)])
    pos = {lay[k]: i for i, k in enumerate(K)}
    g = pos['g']
    return tuple(-Yr[g, pos[t]].imag / w for t in ('s', 'd', 'b'))


@pytest.mark.parametrize('cls', [eh.MosLevel3MeyerHdl, eh.MosLevel2MeyerHdl])
def test_the_capacitances_are_devqmeyers_in_both_modes(cls):
    """Accumulation, depletion, the near-threshold strip, triode and
    saturation, MAGIC_VDS's sliver included -- ngspice's DEVqmeyer, doubled,
    with the model's own `von` (`vto` on this card) and `vdsat` (`vgs -
    von`); in reverse (``vds < 0``) SPICE's swap: the drain as the source."""
    e = _device(cls)
    n = 0
    for vgs in np.linspace(-1.5, 3.0, 31):
        for vds in (0.0, 0.01, 0.03, 0.2, 1.0, 2.5):
            vdsat = max(vgs - VTO, 0.0)
            got = _meyer_row(e, vds, vgs, 0.0, 0.0)
            want = meyer(vgs, vgs - vds, VTO, vdsat, PHI, COX)
            assert got == pytest.approx(want, rel=0, abs=1e-9 * COX), (vgs, vds, got, want)
            ## reverse: the drain at -vds, so the drain is the source
            vgd = vgs + vds
            got_r = _meyer_row(e, -vds, vgs, 0.0, 0.0)
            ws, wd, wb = meyer(vgd, vgs, VTO, max(vgd - VTO, 0.0), PHI, COX)
            if vds > 0.0:
                assert got_r == pytest.approx((wd, ws, wb), rel=0, abs=1e-9 * COX), (vgs, vds)
            n += 1
    assert n == 186


def test_the_pmos_is_the_nmos_mirrored():
    n, p = _device(eh.MosLevel3MeyerHdl), _device(eh.MosLevel3PmosMeyerHdl)
    for vgs, vds in ((2.0, 0.1), (1.0, 2.0), (0.6, 0.5), (-0.5, 1.0)):
        assert _meyer_row(p, -vds, -vgs, 0.0, 0.0) == pytest.approx(
            _meyer_row(n, vds, vgs, 0.0, 0.0), rel=1e-12, abs=1e-12 * COX)


@pytest.mark.parametrize('pair', [(eh.MosLevel3MeyerHdl, eh.MosLevel3GateChargeHdl),
                                  (eh.MosLevel2MeyerHdl, eh.MosLevel2GateChargeHdl)],
                         ids=['level3', 'level2'])
def test_the_dc_current_is_the_charge_classes(pair):
    """Meyer adds nothing at DC: at a state with no time derivative every
    terminal current is the gate-charge class's, byte for byte -- SPICE's
    junction gmin on (`test_device_gmin` leaves these classes out: they have
    no PCNR device)."""
    m, c = _device(pair[0], gmin=1e-3), _device(pair[1], gmin=1e-3)
    for vd, vg, vb in ((0.1, 2.0, 0.0), (2.5, 1.5, -1.0), (-1.0, 2.0, 0.0), (1.0, 0.3, 0.0)):
        xm, lm = _state(m, vd, vg, 0.0, vb)
        xc, lc = _state(c, vd, vg, 0.0, vb)
        im, ic = np.asarray(m.i(xm), float), np.asarray(c.i(xc), float)
        for t in 'dgsb':
            assert im[lm.index(t)] == ic[lc.index(t)], (vd, vg, vb, t)


def test_the_meyer_classes_decline_pcnr_by_their_own_reason():
    for cls in (eh.MosLevel3MeyerHdl, eh.MosLevel2PmosMeyerHdl):
        _device(cls)
        assert cls._hdl_info['pcnr_refusal'] == cls.pcnr_declined


def test_the_importer_builds_them_on_request_and_only_then(tmp_path):
    """`mos_charge='meyer'` maps levels 2 and 3 to the Meyer classes, the
    default keeps the charge-conserving ones, and a level-1 card with TOX
    (no Meyer class yet) is refused, as an unknown value is."""
    from pycircuit.circuit import spice_import
    deck = tmp_path / 'm.cir'
    deck.write_text('m\n.model n3 nmos level=3 vto=0.8 tox=3e-8 nsub=1e16\n'
                    '.model p2 pmos level=2 vto=-0.8 tox=3e-8 nsub=1e16\n'
                    'vd d 0 dc 1\nm1 d d 0 0 n3 w=10u l=2u\nm2 d d 0 0 p2 w=10u l=2u\n.op\n.end\n')
    def mos(imp):
        return [e for e in imp.circuit.elements.values() if 'Mos' in type(e).__name__]
    ## (an imported element is its class's collapse specialisation: isinstance)
    m3, p2 = mos(spice_import.import_netlist(str(deck), mos_charge='meyer'))
    assert isinstance(m3, eh.MosLevel3MeyerHdl) and isinstance(p2, eh.MosLevel2PmosMeyerHdl)
    m3, p2 = mos(spice_import.import_netlist(str(deck)))
    assert isinstance(m3, eh.MosLevel3GateChargeHdl) and isinstance(p2, eh.MosLevel2PmosGateChargeHdl)
    deck.write_text('m\n.model n1 nmos level=1 vto=0.8 tox=3e-8\nvd d 0 dc 1\n'
                    'm1 d d 0 0 n1 w=10u l=2u\n.op\n.end\n')
    with pytest.raises(spice_import.SpiceImportError):
        spice_import.import_netlist(str(deck), mos_charge='meyer')
    with pytest.raises(ValueError):
        spice_import.import_netlist(str(deck), mos_charge='spice')


def test_an_imported_inverter_chain_runs_meyer_on_the_c_core_and_switches_sooner(tmp_path):
    """The inverter chain of `test_mos_gate_charge` in the Meyer mode: its
    Newton the C core's, and its output switches sooner than the charge-
    conserving model's (whose gate-drain capacitance in triode is larger:
    the gold decks' late edges, 2026-10-17)."""
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
    crossing = {}
    for mode in ('meyer', 'conserving'):
        imp = spice_import.import_netlist(str(deck), mos_charge=mode)
        tr, kw = imp.transient()
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tr.solve(**kw)
        if mode == 'meyer':
            assert _paths.since(before).get('newton_c:served', 0) > 50
        t, vb = np.asarray(res.v('b').x[0], float), np.asarray(res.v('b').y, float)
        assert vb.min() < 0.5 and vb.max() > 4.5
        k = int(np.argmax(vb > 2.5))
        crossing[mode] = t[k - 1] + (t[k] - t[k - 1]) * (2.5 - vb[k - 1]) / (vb[k] - vb[k - 1])
    assert crossing['meyer'] < crossing['conserving'], crossing
