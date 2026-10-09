"""MOS level 2 (`elements_hdl.MosLevel2Hdl`, the SPICE benchmark plan's
stage 9, 2026-10-09) against an independent transcription of ngspice-47's
`mos2load.c` (`_mos2_reference`), which is itself held to ngspice's own
currents (recorded in `data/mos2_ngspice_iv.json`: 13 cards, every level-2
feature, 1560 points).  And `hdl.quartic_root`, the velocity-saturation
root the model needs."""
import json
import math
import os

import numpy as np
import pytest
import sympy

from pycircuit.circuit import circuit, constants, hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.tests import _mos2_reference as R

with open(os.path.join(os.path.dirname(__file__), 'data', 'mos2_ngspice_iv.json')) as _fh:
    DATA = json.load(_fh)
#: pycircuit's constants (SI-2019 q and k, eps0 8.854187817e-12): the
#: model's, against ngspice's 3e-6 away -- the reference takes either
PYCIRCUIT = {'CHARGE': constants.qelectron, 'CONSTboltz': constants.kboltzmann,
             'EPS0': 8.854187817e-12}


@pytest.fixture
def model_constants():
    R.use_constants(**PYCIRCUIT)
    yield
    R.use_constants(**R.NGSPICE)


def _model(card, cls=eh.MosLevel2Hdl, **extra):
    kw = {('lambda' if k == 'lambd' else k): v for k, v in card.items()}
    return cls('d', 'g', 's', 'b', IS=1e-30, **kw, **extra)


def _channel_points(points):
    """The recorded points whose bulk junctions stay off (the reference
    has none): neither junction above 0.45 V forward."""
    return [p for p in points if max(p[2], p[2] - p[1]) <= 0.45]


def test_the_reference_is_ngspice_s_level_2():
    """The transcription against ngspice-47's own currents, every card:
    within ngspice's solve tolerance (measured worst 3.8e-9 of a card's
    peak current; pointwise 1e-7 on the mosrect card at |vds| 6 V)."""
    R.use_constants(**R.NGSPICE)
    for name, rec in DATA['cards'].items():
        c = R.card(**rec['card'])
        pts = np.array(rec['points'])
        ref = np.array([R.cdrain(c, vg, vd, vb) for vg, vd, vb, _i in pts])
        scale = max(np.max(np.abs(pts[:, 3])), 1e-12)
        assert np.max(np.abs(ref - pts[:, 3])) / scale < 2e-8, name
        big = np.abs(pts[:, 3]) > 1e-9
        assert np.max(np.abs(ref[big] / pts[big, 3] - 1.0)) < 1e-6, name


def test_the_level2_channel_is_the_reference(model_constants):
    """The model's drain current against the transcription at the recorded
    points (junctions off), on the model's own constants: rounding only
    (measured worst 4e-11 relative, on a 1.7e-8 A subthreshold point)."""
    circuit.default_toolkit = circuit.numeric
    for name, rec in DATA['cards'].items():
        m = _model(rec['card'])
        c = R.card(**rec['card'])
        for vg, vd, vb, _i in _channel_points(rec['points']):
            x = np.zeros(m.n)
            x[0], x[1], x[3] = vd, vg, vb
            got = float(np.asarray(m.i(x), dtype=float)[0])
            ref = R.cdrain(c, vg, vd, vb)
            assert abs(got - ref) <= 1e-9 * abs(ref) + 1e-21, (name, vg, vd, vb, got, ref)


@pytest.mark.parametrize('name', ['nsubxj', 'nfs', 'vmax', 'vmaxall', 'short', 'mosrect'])
def test_the_level2_jacobian_is_the_reference_s_derivative(model_constants, name):
    """gm, gds and gmb from the model's Jacobian (exact derivatives of its
    current, the quartic root's implicit) against central differences of
    the transcription, off the seams (vbs = 0 is the body term's, `von`
    the subthreshold clamp's): measured worst 4e-7 -- the differences'
    own truncation."""
    circuit.default_toolkit = circuit.numeric
    card = DATA['cards'][name]['card']
    m = _model(card)
    c = R.card(**card)
    h = 1e-6
    for vbs in (-0.5, -1.5):
        for vgs in (1.0, 1.5, 2.5):
            for vds in (-0.35, 0.3, 1.0, 2.5):
                x = np.zeros(m.n)
                x[0], x[1], x[3] = vds, vgs, vbs
                G = np.asarray(m.G(x), dtype=float)
                for col, dv in ((0, (0, 1, 0)), (1, (1, 0, 0)), (3, (0, 0, 1))):
                    def f(s, dv=dv, vgs=vgs, vds=vds, vbs=vbs):
                        return R.cdrain(c, vgs + s * dv[0], vds + s * dv[1], vbs + s * dv[2])
                    fd = (f(h) - f(-h)) / (2 * h)
                    scale = max(abs(fd), abs(G[0, col]), 1e-9)
                    assert abs(G[0, col] - fd) / scale < 2e-6, (name, vgs, vds, vbs, col, G[0, col], fd)


def test_the_level2_c_kernel_is_numpy_s_bit_for_bit():
    """The compiled C chain (with `quartic_root`'s C twin) against the numpy
    chain, every value of `i`, `G`, `q`, `C` on a grid through every region
    of the all-features card: bit for bit (measured: all 12 cards, 32760
    values each)."""
    circuit.default_toolkit = circuit.numeric
    card = DATA['cards']['vmaxall']['card']
    cls = eh.MosLevel2Hdl
    pts = [np.array([vd, vg, 0.0, vb]) for vb in (0.0, -1.0, 0.3)
           for vg in np.arange(-0.5, 3.01, 0.5) for vd in np.arange(-1.0, 3.01, 0.4)]

    def values():
        m = _model(card)
        return [np.concatenate([np.asarray(m.i(x), dtype=float),
                                np.asarray(m.G(x), dtype=float).ravel(),
                                np.asarray(m.q(x), dtype=float),
                                np.asarray(m.C(x), dtype=float).ravel()]) for x in pts]
    hdl.set_backend('numpy', cls)
    try:
        ref = values()
    finally:
        hdl.set_backend(None, cls)
    hdl.set_backend('c', cls)
    try:
        got = values()
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
    finally:
        hdl.set_backend(None, cls)
    for a, b in zip(ref, got):
        assert np.array_equal(a, b), np.max(np.abs(a - b))


def test_a_level2_pmos_is_the_nmos_mirrored():
    """`MosLevel2PmosHdl` with the same (positive) card at mirrored terminal
    voltages carries the mirrored current, both modes."""
    circuit.default_toolkit = circuit.numeric
    card = DATA['cards']['vmaxall']['card']
    n, p = _model(card), _model(card, cls=eh.MosLevel2PmosHdl)
    for vg, vd, vb in ((1.5, 2.0, -1.0), (2.5, -0.7, -0.5), (0.9, 1.0, 0.0)):
        x = np.array([vd, vg, 0.0, vb])
        assert np.asarray(p.i(-x), dtype=float)[0] == pytest.approx(
            -float(np.asarray(n.i(x), dtype=float)[0]), rel=1e-12)


def _draw_quartics(count, seed):
    """Level-2 saturation quartics drawn over SPICE's own coefficient
    construction (`mos2load.c`), roots present and absent."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(count):
        gammad = 0.0 if rng.random() < 0.2 else rng.uniform(0.0, 1.5)
        v2 = max(rng.uniform(0.4, 1.2) - rng.uniform(-3.0, 0.5), 0.05)
        v1 = rng.uniform(-0.5, 5.0) + v2
        xv = rng.uniform(0.01, 20.0)
        out.append((gammad / 0.75, -2.0 * (v1 + xv), -2.0 * gammad * xv,
                    2.0 * v1 * (v2 + xv) - v2 * v2 - 4.0 / 3.0 * gammad * math.sqrt(v2) ** 3))
    return out


def test_quartic_root_is_spice_s_root():
    """`hdl.quartic_root`'s value: the transcription's `_quartic` bit for bit
    where it finds a root, 0.0 where it finds none (measured 11316 roots,
    8684 none, of 20000 drawn); the array path is the scalar per element,
    bit for bit."""
    found = none = 0
    qs = _draw_quartics(4000, 7)
    for args in qs:
        try:
            ref = R._quartic(*args)
        except ValueError:          # (a negative radicand: NaN in the C)
            ref = None
        got = hdl._quartic_root_scalar(*args)
        if ref is None:
            assert got == 0.0, args
            none += 1
        else:
            assert got == ref, args
            found += 1
    assert found > 1000 and none > 500, (found, none)
    arr = np.array(qs).T
    vec = hdl._quartic_root_numpy(*arr)
    assert np.array_equal(vec, [hdl._quartic_root_scalar(*a) for a in qs])


def test_quartic_root_s_derivative_is_implicit():
    """`quartic_root`'s derivative, ``-(dP/da_k)/P'(x)``, against central
    differences of its value where a root is found, every coefficient.
    The step is 1e-4, not 1e-7: SPICE's closed form (Ferrari through a
    resolvent cubic) gives its root to ~1e-8 relative, and a 1e-7 step
    turned that into a 0.26 error on a derivative of 0.11 (checked by hand:
    ``-x^3/P'(x)`` at x = 1.154, P' = -13.79)."""
    a = sympy.symbols('a1 b1 c1 d1')
    expr = hdl.quartic_root(*a)
    grads = [sympy.lambdify(a, expr.diff(s), modules=[{'quartic_root': hdl._quartic_root_numpy},
                                                       'numpy']) for s in a]
    checked = 0
    for args in _draw_quartics(400, 11):
        x = hdl._quartic_root_scalar(*args)
        if x == 0.0:
            continue
        for k in range(4):
            h = 1e-4 * max(abs(args[k]), 1.0)
            up = list(args); dn = list(args)
            up[k] += h; dn[k] -= h
            xu, xd = hdl._quartic_root_scalar(*up), hdl._quartic_root_scalar(*dn)
            if xu == 0.0 or xd == 0.0 or abs(xu - xd) > 1e-3 * x:
                continue                     # (the selected root changed)
            fd = (xu - xd) / (2 * h)
            an = float(grads[k](*args))
            assert abs(an - fd) <= 2e-3 * max(abs(fd), 1e-3), (args, k, an, fd)
            checked += 1
    assert checked > 400, checked


def test_a_level2_card_imports_with_vto_from_nsub(model_constants, tmp_path):
    """The importer maps LEVEL=2 to the gate-charge class, and a card with
    NSUB and no VTO keeps NSS/TPG for the model to derive VTO from, as
    SPICE does (levels 1 and 3 still refuse it)."""
    from pycircuit.circuit import spice_import
    circuit.default_toolkit = circuit.numeric
    card = DATA['cards']['vtonsub']['card']
    deck = tmp_path / 'l2.cir'
    deck.write_text('* l2\nM1 d g 0 b NX L=2u W=20u\nVD d 0 1\nVG g 0 1.5\nVB b 0 -1\n'
                    '.MODEL NX NMOS (LEVEL=2 IS=1e-30 NSUB=2e16 TOX=4e-8 UO=500 NSS=1e10)\n'
                    '.op\n.end\n')
    imp = spice_import.import_netlist(str(deck))
    m1 = next(e for k, e in imp.circuit.elements.items() if k.lower() == 'm1')
    ## (the importer builds the class's collapse variant: a subclass)
    assert isinstance(m1, eh.MosLevel2GateChargeHdl), type(m1)
    x = np.zeros(m1.n)
    x[0], x[1], x[3] = 1.0, 1.5, -1.0
    ref = R.cdrain(R.card(**card), 1.5, 1.0, -1.0)
    assert float(np.asarray(m1.i(x), dtype=float)[0]) == pytest.approx(ref, rel=1e-9)
