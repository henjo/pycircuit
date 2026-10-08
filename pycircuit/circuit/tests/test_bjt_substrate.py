"""The Gummel-Poon substrate junction (`GummelPoon{Npn,Pnp}4Hdl`, the SPICE
benchmark plan's stage 7): its charge against a transcription of ngspice's
(`bjtload.c`), where it attaches, and that a zero junction is the
3-terminal device.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node, defaultepar

CARD = {'IS': 1e-16, 'bf': 100.0, 'vaf': 50.0, 'cje': 1e-12, 'cjc': 5e-13, 'tf': 1e-10,
        'rb': 50.0, 'rc': 10.0, 're': 1.0}
SUB = {'cjs': 2e-12, 'vjs': 0.6, 'mjs': 0.4}


def ngspice_qsub(v, cj, vj, m):
    """`bjtload.c`'s substrate charge: depletion below zero bias, the
    capacitance linearly extended above."""
    if v < 0:
        arg = 1 - v / vj
        return vj * cj * (1 - arg * arg ** -m) / (1 - m)
    return v * cj * (1 + m * v / (2 * vj))


def _device(cls, **kw):
    e = cls(*[Node(n) for n in 'cbes'[:len(cls.terminals)]], **kw)
    e.update_iparv()
    return e


def _qC(e, x):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return (np.asarray(e.q(x, defaultepar), float),
                np.asarray(e.C(x, defaultepar), float))


def _x(e, v):
    """The state with the terminals at `v` (c, b, e, s) and each internal
    node at its terminal's potential."""
    names = _names(e)
    out = np.zeros(len(names))
    for k, nm in enumerate(names):
        t = nm[0]
        out[k] = v['cbes'.index(t)]
    return out


def _names(e):
    """The state's names in order (`hdl.x_layout`: (index, name, kind))."""
    return [name for _k, name, _kind in hdl.x_layout(type(e))]


@pytest.mark.parametrize('cls, sign', [('GummelPoonNpn4Hdl', 1.0), ('GummelPoonPnp4Hdl', -1.0)])
def test_the_substrate_charge_is_ngspices_on_the_collector(cls, sign):
    """Vertical (`subs = 1`, Xyce's): the substrate's own charge is
    ngspice's `qsub` of `type*(V(s) - V(c))`, below zero bias and above,
    and it is the collector's image -- nothing on the base or emitter."""
    e = _device(getattr(eh, cls), **CARD, **SUB)
    zero = _device(getattr(eh, cls), **CARD)
    names = _names(e)
    s, ci = names.index('s'), names.index('ci')
    for vsc in (-5.0, -1.0, -0.3, 0.0, 0.2, 0.5):
        v = sign * np.array([0.0, 0.7, 0.0, vsc])          # c b e s
        x = _x(e, v)
        q, _C = _qC(e, x)
        q0, _ = _qC(zero, x)
        want = sign * ngspice_qsub(vsc, SUB['cjs'], SUB['vjs'], SUB['mjs'])
        assert q[s] - q0[s] == pytest.approx(want, rel=1e-12, abs=1e-27), vsc
        dq = q - q0
        assert dq[ci] == pytest.approx(-want, rel=1e-12, abs=1e-27)
        others = [k for k in range(len(q)) if k not in (s, ci)]
        assert np.abs(dq[others]).max() <= 1e-27


def test_lateral_puts_it_on_the_base():
    """`subs = -1` (lateral, ngspice's p-n-p default): the junction is
    between the substrate and the INTERNAL base, polarity `type*subs`."""
    e = _device(eh.GummelPoonPnp4Hdl, **CARD, **SUB, subs=-1.0)
    zero = _device(eh.GummelPoonPnp4Hdl, **CARD)
    names = _names(e)
    s, bi = names.index('s'), names.index('bi')
    for vsb in (-2.0, -0.4, 0.3):
        v = np.array([0.0, -0.7, 0.0, -0.7 + vsb])         # V(s) - V(b) = vsb
        x = _x(e, v)
        dq = _qC(e, x)[0] - _qC(zero, x)[0]
        want = ngspice_qsub(vsb, SUB['cjs'], SUB['vjs'], SUB['mjs'])   # ttype = (-1)(-1)
        assert dq[s] == pytest.approx(want, rel=1e-12, abs=1e-27)
        assert dq[bi] == pytest.approx(-want, rel=1e-12, abs=1e-27)


@pytest.mark.parametrize('n3, n4', [('GummelPoonNpnHdl', 'GummelPoonNpn4Hdl'),
                                    ('GummelPoonPnpHdl', 'GummelPoonPnp4Hdl')])
def test_a_zero_substrate_junction_is_the_three_terminal_device(n3, n4):
    """With `cjs = 0` the 4-terminal device is the 3-terminal one, bit for
    bit on the shared unknowns, the substrate carrying nothing."""
    e3 = _device(getattr(eh, n3), **CARD)
    e4 = _device(getattr(eh, n4), **CARD)
    l3, l4 = _names(e3), _names(e4)
    pick = [l4.index(n) for n in l3]
    rng = np.random.default_rng(5)
    for _ in range(20):
        x4 = rng.uniform(-1, 1, len(l4))
        x3 = x4[pick]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            for m in ('i', 'q'):
                a = np.asarray(getattr(e3, m)(x3, defaultepar), float)
                b = np.asarray(getattr(e4, m)(x4, defaultepar), float)
                assert np.array_equal(a, b[pick]), m
                rest = [k for k in range(len(l4)) if k not in pick]
                assert np.all(b[rest] == 0.0), m


def test_the_importer_builds_the_substrate_where_the_card_asks(tmp_path):
    """A card with CJS gives the 4-terminal class: on the line's fourth
    node, or ground on a 3-node line (SPICE's default); vertical in
    Xyce's dialect, a p-n-p lateral in ngspice's; written back with SUBS
    explicit, so ngspice reads the same junction."""
    from pycircuit.circuit import spice_import

    deck = tmp_path / 'q.cir'
    deck.write_text("""title
q1 c b e qn
q2 c b e sub qp
.model qn npn bf=80 cjs=1p vjs=0.6 mjs=0.3
.model qp pnp bf=40 cjs=2p
vc c 0 5
vb b 0 0.7
.tran 1n 10n
.end
""")
    for dialect, psubs in (('xyce', 1.0), ('ngspice', -1.0)):
        imp = spice_import.import_netlist(str(deck), dialect=dialect)
        q = {m.name: m for m in imp.elements if m.name in ('q1', 'q2')}
        assert q['q1'].cls is eh.GummelPoonNpn4Hdl and q['q1'].nodes[3] == '0'
        assert q['q1'].params['cjs'] == 1e-12 and q['q1'].params['mjs'] == 0.3
        assert q['q1'].params['subs'] == 1.0
        assert q['q2'].cls is eh.GummelPoonPnp4Hdl and q['q2'].nodes[3] == 'sub'
        assert q['q2'].params['subs'] == psubs
        written = tmp_path / f'w_{dialect}.cir'
        imp.write_ngspice(str(written), probes=['c'])
        cards = [ln.lower() for ln in written.read_text().splitlines()
                 if ln.lower().startswith('.model')]
        assert any('cjs=' in c and f'subs={psubs!r}' in c for c in cards if ' pnp' in c), cards


#: ngspice-47 on this deck (`.dc vb 0.6 1.0 0.1`, `.print dc i(vc) i(vb)`,
#: `.options reltol=1e-9 abstol=1e-18 vntol=1e-12`), 2026-10-08: the base
#: current crosses IRB = 100 uA in the sweep, so the base resistance moves
#: from RB = 100 toward RBM = 10 along SPICE's law.  ⚠ At ngspice's default
#: tolerances a SWEEP carries ~1e-3 of its own error (each point continued
#: from the last): the first reference taken that way read as a 9e-4
#: disagreement that was ngspice's.
IRB_DECK = """irb reference
vc c 0 5
vb b 0 {vbe}
q1 c b 0 QI
.model QI NPN IS=1e-15 BF=100 BR=2 VAF=50 RB=100 RBM=10 IRB={irb} RE=1 RC=5 ISE=1e-14 NE=1.5
.end
"""
IRB_NGSPICE = {0.6: (-1.29017e-05, -1.70579e-07), 0.7: (-5.89129e-04, -6.09017e-06),
               0.8: (-1.34061e-02, -1.29162e-04), 0.9: (-6.07214e-02, -5.78865e-04),
               1.0: (-1.26997e-01, -1.21333e-03)}


def _irb_currents(tmp_path, vbe, irb):
    from pycircuit.circuit import spice_import
    from pycircuit.circuit.dcanalysis import DC
    deck = tmp_path / f'irb_{vbe}_{irb}.cir'
    deck.write_text(IRB_DECK.format(vbe=vbe, irb=irb))
    imp = spice_import.import_netlist(str(deck), strict=False)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = DC(imp.circuit).solve()
    return float(r.i('vc.plus')), float(r.i('vb.plus'))


def test_irb_is_spices_base_resistance_law(tmp_path):
    """`irb > 0`: the current-crowding base resistance (ngspice
    `bjtload.c`), against ngspice's own sweep of the same card to 1e-4 --
    across the base current's crossing of IRB; and `irb = 0` (the `qb`
    law) is measurably another device there, so the agreement is the
    law's and not a coincidence."""
    worst, gap = 0.0, 0.0
    for vbe, (ic_ng, ib_ng) in IRB_NGSPICE.items():
        ic, ib = _irb_currents(tmp_path, vbe, 1e-4)
        worst = max(worst, abs(ic - ic_ng) / abs(ic_ng), abs(ib - ib_ng) / abs(ib_ng))
        ic0, ib0 = _irb_currents(tmp_path, vbe, 0)
        gap = max(gap, abs(ib0 - ib_ng) / abs(ib_ng))
    assert worst < 1e-4, worst
    assert gap > 1e-2, gap


PTF_DECK = """ce excess phase
vcc vcc 0 10
vin in 0 sin(0.75 0.02 100meg 0 0)
rb in b 200
rc vcc c 1k
q1 c b 0 qx
.model qx npn is=1e-16 bf=100 vaf=50 tf=1n ptf={ptf} cje=0.2p cjc=0.1p rb=10
.tran 0.01n 40n
.end
"""
#: ngspice-47's collector on PTF_DECK at PTF = 60 (`.options reltol=1e-6`,
#: `.tran 0.005n 40n 0 0.005n`), 2026-10-08, interpolated at these times.
PTF_NGSPICE = {30e-9: 9.661455, 32.5e-9: 9.574411, 35e-9: 9.426587, 37.5e-9: 9.512939}


def _ptf_collector(tmp_path, ptf):
    from pycircuit.circuit import spice_import
    deck = tmp_path / f'ptf_{ptf}.cir'
    deck.write_text(PTF_DECK.format(ptf=ptf))
    imp = spice_import.import_netlist(str(deck))
    tr, kw = imp.transient()
    tr.par.timestep_max = 0.01e-9
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = tr.solve(**kw)
    v = r.v('c')
    return np.asarray(v.x[0], float), np.asarray(v.y, float)


def test_ptf_is_spices_excess_phase(tmp_path):
    """PTF (Weil's excess phase, ngspice `bjtload.c`): a common-emitter
    stage at 100 MHz with TF = 1 ns, against ngspice's collector to 0.5 mV
    (0.2 % of its 0.28 V swing, measured 0.12; ngspice integrates the delay inside the
    device by backward Euler, here the simulator integrates its two states)
    -- where PTF = 0 is ~80 mV away, so the agreement is the delay's.  The
    operating point does not see PTF (the delay's DC gain is one)."""
    t, v = _ptf_collector(tmp_path, 60)
    t0, v0 = _ptf_collector(tmp_path, 0)
    for tt, want in PTF_NGSPICE.items():
        assert abs(np.interp(tt, t, v) - want) < 5e-4, (tt, np.interp(tt, t, v), want)
    g = np.linspace(20e-9, 40e-9, 801)
    assert np.abs(np.interp(g, t, v) - np.interp(g, t0, v0)).max() > 0.05
    assert v[0] == pytest.approx(v0[0], rel=1e-12)
