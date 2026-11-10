"""PSP's limiter in C.  `compact.PspMosLongChannel` declares SPICE's MOSFET
group -- `fetlim` on gate-source, `limvds` on drain-source, `pnjlim` on the
two bulk junctions, limited as one -- and an identity probe on its noise
node, so its `limit` is the generated one, its C kernel the generated
kernel (`_hdl_climit`) and its PCNR refine the generated refine
(`_pcnr_c`).  Here: the declaration on both channel types and cards; the
walk and the C Newton serving PSP circuits with the Python path's answers;
PCNR's refine kernel and DC solve in C against their Python on PSP.
`_stage` and `_chain` are the PSP fixtures other C-twin tests share."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _hdl_batch,
    _hdl_climit,
    _paths,
    _pcnr_c,
    circuit,
    compact,
    pcnr,
)
from pycircuit.circuit.elements import VS, R, SubCircuit, VSin, gnd
from pycircuit.circuit.tests.test_pcnr_c import _start, _twin
from pycircuit.circuit.transient import Transient

PSP = compact.PspMosLongChannel


def _stage(pmos=False):
    c = SubCircuit()
    for n in ('g', 'd', 'vdd'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=1.2)
    if pmos:
        c['vg'] = VSin('g', gnd, v=0.5, va=2e-2, freq=1e6)
        c['rl'] = R('d', gnd, r=5e3)
        c['M'] = compact.PspPmosLongChannel('d', 'g', 'vdd', 'vdd', fnt=1.0)
    else:
        c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
        c['rl'] = R('vdd', 'd', r=5e3)
        c['M'] = PSP('d', 'g', gnd, gnd, fnt=1.0)
    return c


def _chain(n=3, va=2e-2):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=0.7, va=va, freq=1e6)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = PSP(f'd{k}', f'g{k}' if k == 0 else f'd{k - 1}', gnd, gnd, fnt=1.0)
    return c


def _cmos(n=6, vin=0.3, **kw):
    """`n` PSP inverters in a row, each p over n, the input held at `vin`."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vin'] = VS('d_in', gnd, v=vin)
    prev = 'd_in'
    for k in range(n):
        c[f'MP{k}'] = compact.PspPmosLongChannel(f'd{k}', prev, 'vdd', 'vdd', w=20e-6, l=1e-6, **kw)
        c[f'MN{k}'] = PSP(f'd{k}', prev, gnd, gnd, w=10e-6, l=1e-6, **kw)
        prev = f'd{k}'
    return c


## ------------------------------------------------------------------------
## The declaration.

@pytest.mark.parametrize('card', [{}, {'rg': 10.0, 'fnt': 1.0}], ids=['default', 'rg-fnt'])
@pytest.mark.parametrize('cls', [PSP, compact.PspPmosLongChannel], ids=['nmos', 'pmos'])
def test_psp_declares_spices_group_and_pcnr_takes_it(cls, card):
    """Both channel types, on the default card (gate resistor and noise
    node collapsed) and with both kept: SPICE's four probes over the four
    terminals, the noise node's identity probe where it exists, the
    generated `limit` with its C kernel, and vector PCNR's route --
    which, without the declaration, refuses the channel current for
    reading three node voltages."""
    el = cls('d', 'g', 's', 'b', w=10e-6, l=1e-6, **card)
    el.update_iparv()
    info = type(el)._hdl_info
    kinds = sorted(k for _rows, k, _m, _p in info['limit_spec'])
    assert kinds == sorted(['fet', 'vds', 'pnj', 'pnj'] + (['id'] if card else [])), kinds
    assert _hdl_batch.is_generated(type(el), 'limit')
    if info.get('_c_limit_status') != 'c':
        pytest.skip(f"no C limit kernel here: {info.get('_c_limit_status')}")
    assert info['pcnr_refusal'] is None
    assert sorted(info['pcnr_vector']['kinds']) == sorted(['fet', 'vds', 'pnj'] + (['id'] if card else []))


## ------------------------------------------------------------------------
## The walk and the C Newton.

def _run(build, on, monkeypatch, **kw):
    monkeypatch.setattr(_hdl_climit, 'ENABLED', on)
    before = _paths.snapshot()
    with warnings.catch_warnings(record=True) as W:
        warnings.simplefilter('always')
        ## (a finalizer's warning is not the run's: a file another test left in
        ## a reference cycle closes whenever the collector runs)
        warnings.simplefilter('ignore', ResourceWarning)
        tr = Transient(build(), toolkit=circuit.numeric, **kw.pop('make', {}))
        res = tr.solve(**kw)
    st_ = {k: getattr(tr.statistics, k) for k in tr.statistics.__slots__ if 'seconds' not in k}
    return (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
            st_, sorted(str(w.message) for w in W), _paths.since(before))


def _same(build, monkeypatch, served_min=1, **kw):
    a = _run(build, False, monkeypatch, **dict(kw))
    b = _run(build, True, monkeypatch, **dict(kw))
    assert a[0] == b[0] and a[1] == b[1], 'the solution moved'
    assert a[2] == b[2], ('the statistics moved', a[2], b[2])
    assert a[3] == b[3], 'the warnings moved'
    assert b[4].get('newton_c:served', 0) >= served_min, b[4]
    return a, b


@pytest.mark.parametrize('pmos', [False, True], ids=['nmos', 'pmos'])
@pytest.mark.parametrize('fixed', [True, False], ids=['fixed', 'adaptive'])
def test_the_c_newton_serves_psp_with_the_python_answers(pmos, fixed, monkeypatch):
    kw = {'tend': 1e-6, 'timestep': 2e-8}
    if fixed:
        kw['fixed_timestep'] = True
    a, _b = _same(lambda: _stage(pmos), monkeypatch, **kw)
    assert a[4].get('newton_c:served', 0) == 0


def test_a_chain_driven_hard_enough_to_limit(monkeypatch):
    """A chain whose gates swing by 6 V a step: the law changes iterates
    (counted on the Python path), and the walk's kernel changes the same
    rows by the same amounts."""
    _same(lambda: _chain(3, va=6.0), monkeypatch, tend=1e-6, timestep=1e-7, fixed_timestep=True)
    cls = type(_chain(1)['M0'])
    real, bites = cls.limit, []

    def counting(self, x, *a, **k):
        out = real(self, x, *a, **k)
        bites.append(not np.array_equal(np.asarray(out), np.asarray(x)))
        return out
    monkeypatch.setattr(cls, 'limit', counting)
    _run(lambda: _chain(3, va=6.0), True, monkeypatch, tend=1e-6, timestep=1e-7,
         fixed_timestep=True)
    assert sum(bites) >= 5, (sum(bites), len(bites))


## ------------------------------------------------------------------------
## PCNR in C.

@pytest.mark.parametrize('card', [{}, {'fnt': 1.0}], ids=['default', 'noise-node'])
def test_the_refine_kernel_is_the_devices_limit(card):
    """PSP's refine kernel (`_pcnr_c.render_refine`) gives
    `PcnrDevice.limit` bit for bit -- `von` at the last accepted
    sub-vector, the laws, the redundant `(b,d)` probe's Kruskal and
    potentials, the noise node's identity -- on draws where the laws
    bite."""
    c = _cmos(1, vin=0.6, **card)
    for dev in pcnr.pcnr_devices(c):
        params, epar = dev.params(), circuit.defaultepar
        rng = np.random.default_rng(7)
        scale = np.array([0.8 if k == 'pnj' else 5.0 for k in dev.kinds])
        bit = 0
        for _ in range(300):
            vo = rng.uniform(-1, 1, dev.m) * scale
            vn = vo + rng.standard_normal(dev.m) * scale * 10.0 ** rng.uniform(-3, 1.5)
            xs = rng.uniform(-3, 3, dev.t)
            ref = dev.limit(vn.copy(), vo.copy(), params, epar, xs.copy())
            got = _pcnr_c.refine(dev, vn, vo, epar, xs)
            assert got is not None
            assert np.array_equal(np.asarray(ref).view(np.int64), got.view(np.int64))
            bit += not np.array_equal(ref, vn)
        assert bit > 50, bit


@pytest.mark.parametrize('card', [{}, {'fnt': 1.0}], ids=['default', 'noise-node'])
@pytest.mark.parametrize('start', [0.0, 1.2, 5.0])
def test_a_cmos_chain_solve_in_c_is_the_pythons(card, start, monkeypatch):
    """PCNR's DC on a PSP CMOS chain, C against Python bit for bit, served
    by the C, from three starts -- and converged."""
    c = _cmos(6, **card)
    got, _counts = _twin(c, monkeypatch, x0=_start(c, start))
    assert got[0] == 'ok', got
