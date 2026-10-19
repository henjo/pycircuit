"""SPICE's junction GMIN inside the device models: `gmin` on MOS levels 1, 2
and 3 (NMOS, PMOS, plain and gate charge) and on the SPICE diode's junction
(its thermal, LED and photodiode variants share it) -- `mos*load.c`'s
`CKTgmin*vbs` / `CKTgmin*vbd` and `dioload.c`'s `CKTgmin*vd` (not area-scaled).
The default 0 is the card's physics; the importer sets SPICE's 1e-12."""
import numpy as np
import pytest

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import pcnr
from pycircuit.circuit.analysis import defaultepar
from pycircuit.circuit.elements import SubCircuit
from pycircuit.circuit.tests.test_dc_held import _latch
from pycircuit.circuit.tests.test_pcnr_c import _twin

CARRIERS = sorted(n for n in dir(eh)
                  if (n.startswith('MosLevel') or n in ('DiodeSpiceHdl', 'DiodeSpiceThermalHdl',
                                                        'LedHdl', 'PhotodiodeHdl'))
                  and isinstance(getattr(eh, n), type)
                  and any(p.name == 'gmin' for p in getattr(eh, n).instparams))


def test_every_series_carries_it():
    """The control for the test below: what it enumerates."""
    assert len(CARRIERS) == 20, CARRIERS


def _device(name, g):
    cls = getattr(eh, name)
    c = SubCircuit()
    pins = [f'n{k}' for k in range(len(cls.terminals))]
    for nd in pins:
        c.add_node(nd)
    kw = {'gmin': g}
    if 'w' in {p.name for p in cls.instparams}:
        kw.update(w=10e-6, l=1e-6)
    c['m1'] = cls(*pins, **kw)
    dev, = pcnr.pcnr_devices(c)
    return c, dev


@pytest.mark.parametrize('name', CARRIERS)
def test_gmin_is_a_conductance_across_each_junction(name):
    """`G(g) - G(0)` is `g` stamped across each junction pair -- a MOSFET's
    bulk-drain and bulk-source, a diode's junction -- and `i(g) - i(0)` is
    that matrix times the potentials, at drawn points."""
    g = 1e-3
    c0, d0 = _device(name, 0.0)
    _cg, dg = _device(name, g)
    t = d0.t
    local = {r: k for k, r in enumerate(d0.rows)}
    K = np.zeros((t, t))
    pairs = pcnr.pcnr_junction_pairs(c0)
    assert pairs, name
    for _inst, _dev, ra, rb in pairs:
        a, b = local[ra], local[rb]
        K[a, a] += g
        K[a, b] -= g
        K[b, a] -= g
        K[b, b] += g
    rng = np.random.default_rng(3)
    for _ in range(4):
        pot = rng.uniform(-0.6, 0.6, t)
        dG = (np.asarray(dg.element.G(pot, defaultepar), dtype=float).reshape(t, t)
              - np.asarray(d0.element.G(pot, defaultepar), dtype=float).reshape(t, t))
        di = (np.asarray(dg.element.i(pot, defaultepar), dtype=float).reshape(t)
              - np.asarray(d0.element.i(pot, defaultepar), dtype=float).reshape(t))
        assert np.allclose(dG, K, rtol=1e-9, atol=1e-12 * g), (name, dG, K)
        assert np.allclose(di, K @ pot, rtol=1e-9, atol=1e-12 * g), (name, di, K @ pot)


@pytest.mark.parametrize('g', [1e-12, 1e-6])
def test_the_c_kernels_carry_it(g, monkeypatch):
    """The latch's MOSFETs with the junction gmin on: PCNR's solve in C (the
    class kernels) bit for bit the Python's, the C serving."""
    c = _latch()
    for nm in ('m1', 'm2'):
        c[nm].ipar.gmin = g
    got, _counts = _twin(c, monkeypatch)
    assert got[0] == 'ok', got
