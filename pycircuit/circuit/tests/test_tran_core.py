"""The evaluate core (`pycircuit/circuit/_tran_core.py`, speed round 7's
first commit, 2026-10-03): the transient's three evaluation sites served by
one C call must give the Python path's bytes and leave the Python path's
state, decline to it wherever they cannot serve, and step a transient and a
PSS bit for bit."""
import contextlib
import os
import subprocess
import sys
import types
import warnings

import numpy as np
import pytest

from pycircuit.circuit import PSS, _cbuild, _tran_core, simwarnings
from pycircuit.circuit import circuit as cm
from pycircuit.circuit.elements import VS, Diode, R, SubCircuit, VSin, gnd
from pycircuit.circuit.integrator import (
    EulerIntegrator,
    Gear2Integrator,
    ThetaIntegrator,
    TrapezoidalIntegrator,
)
from pycircuit.circuit.tests.test_hdl_batch import mos_chain
from pycircuit.circuit.transient import Transient

METHODS = ('gear', 'euler', 'trap', 'theta')
INTEGRATOR = {'gear': Gear2Integrator, 'euler': EulerIntegrator,
              'trap': TrapezoidalIntegrator, 'theta': ThetaIntegrator}


@contextlib.contextmanager
def core(on):
    was = _tran_core.CORE
    _tran_core.CORE = on
    try:
        yield
    finally:
        _tran_core.CORE = was


def _stage():
    """The harness's PSP stage: one C-bound device, a lone element."""
    from pycircuit.circuit import compact
    cm.default_toolkit = cm.numeric
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd, fnt=1.0)
    return c


def _tr(cir, method='gear', steps=6):
    """A transient a few steps in, with the step state a residual needs,
    and its last accepted state."""
    tr = Transient(cir, toolkit=cm.numeric, integrator=INTEGRATOR[method]())
    res = tr.solve(tend=steps * 2e-8, timestep=2e-8, fixed_timestep=True)
    if _tran_core.STATUS != 'c':
        pytest.skip(f'the core is off: {_tran_core.STATUS}')
    if _tran_core.core_for(tr) is None:
        from pycircuit.circuit.tests.test_hdl_batch import _why
        facts = _why(cir, 'G')
        if any(v['bound'] for v in facts.values()):
            pytest.fail(f'a bound class the core does not serve: {facts}')
        pytest.skip('the chain is not served (a class not C-bound)')
    x = np.ascontiguousarray(np.asarray(res.x, float)[:, -1])
    tr._u_memo = {}
    return tr, x


def _state(tr):
    out = [tr._iq.tobytes(), np.asarray(tr._Geq).tobytes(), np.asarray(tr._Cmat).tobytes(),
           repr(tr._companion_coeffs), tr._effective_method, type(tr.active_integrator).__name__]
    for name in ('_q_cache', '_C_cache'):
        c = getattr(tr, name, None)
        out.append(None if c is None else (c[0].tobytes(), np.asarray(c[1]).tobytes()))
    return out


def _python_j(tr, x):
    """`jacobian_only`'s body (transient.py), the Python path."""
    _iq, Geq = tr._companion_at(x)
    return tr.cir.G(x, tr.epar) + Geq


def _python_f(tr, x, t):
    """`residual_only`'s body, the Python path."""
    q = tr.cir.q(x, tr.epar)
    iq, _geq = tr.get_diff(q, tr._Cmat)
    u = tr._source_at(t, None)
    return tr.cir.i(x, tr.epar) + iq + u


def _states(x, k, seed=0):
    rng = np.random.default_rng(seed)
    out = [x]
    for _ in range(k):
        out.append(x + rng.normal(0.0, 0.3, x.size))
    out.append(np.full(x.size, 0.7))
    out.append(np.zeros(x.size))
    return out


@pytest.mark.parametrize('method', METHODS)
def test_the_residual_is_the_python_paths_bytes_and_state(method):
    tr, x = _tr(mos_chain(8), method)
    t = 1.1e-7
    for j, xt in enumerate(_states(x, 12)):
        tr._is_first_step = (j % 5 == 4)        # the order drop, where the method has one
        with core(True):
            a = tr._residual_and_jacobian(xt, t)
            sa = _state(tr)
        with core(False):
            b = tr._residual_and_jacobian(xt, t)
            sb = _state(tr)
        tr._is_first_step = False
        assert a[0].dtype == b[0].dtype and a[1].shape == b[1].shape
        assert a[0].tobytes() == b[0].tobytes(), (method, j)
        assert a[1].tobytes() == b[1].tobytes(), (method, j)
        assert sa == sb, (method, j)


@pytest.mark.parametrize('make', [lambda: mos_chain(8), _stage], ids=['chain', 'lone'])
def test_the_converged_point_and_the_chord_residual(make):
    tr, x = _tr(make())
    t = 1.1e-7
    for xt in _states(x, 6, seed=2):
        with core(True):
            none, J = _tran_core.evaluate(tr, xt, t, None, 'j')
            sa = _state(tr)
        assert none is None
        with core(False):
            Jp = _python_j(tr, xt)
            sb = _state(tr)
        assert J.tobytes() == Jp.tobytes() and sa == sb
        ## the chord's residual: the conductance from the held `_Cmat`
        with core(True):
            f = _tran_core.evaluate(tr, xt, t, None, 'f')
            sa = _state(tr)
        with core(False):
            fp = _python_f(tr, xt, t)
            sb = _state(tr)
        assert f.tobytes() == fp.tobytes() and sa == sb


def _patterned():
    """A chain with the pattern on (20 stages: 161 entries of 27 x 27): `(tr,
    x, core)`."""
    tr, x = _tr(mos_chain(20))
    core_ = _tran_core.core_for(tr)
    assert core_.npat, 'the pattern is off: the test says nothing'
    return tr, x, core_


def _lone_negative_zeros(tr, x, which):
    """The entries of `which` ('G' or 'C') whose one stamp at `x` is -0.0."""
    contrib = {}
    for name, el in tr.cir.elements.items():
        idx = np.asarray(tr.cir.elementnodemap[name], dtype=np.int64)
        M = np.asarray(getattr(el, which)(x[idx], tr.epar), float)
        for a in range(M.shape[0]):
            for b in range(M.shape[1]):
                contrib.setdefault((int(idx[a]), int(idx[b])), []).append(M[a, b])
    return [rc for rc, v in contrib.items() if len(v) == 1 and v[0] == 0.0 and np.signbit(v[0])]


def test_a_lone_negative_zero_slot_is_bincounts_positive_zero():
    """On the pattern every entry is bincount's -- from +0.0, its slots added
    in their order: a capacitance entry whose one stamp is -0.0 (the chain's
    gate couplings, `-C` of a zero C) is +0.0, and every byte of the four
    passes the Python's.  (The pattern's first test with a lone -0.0 slot:
    the chains of the core's other tests are too small for the pattern.)"""
    tr, x, _core = _patterned()
    lone = _lone_negative_zeros(tr, x, 'C')
    assert len(lone) >= 10, lone
    got = _tran_core.passes(tr, x, ('G', 'C', 'i', 'q'))
    want = {'G': tr.cir.G(x, tr.epar), 'C': tr.cir.C(x, tr.epar),
            'i': tr.cir.i(x, tr.epar), 'q': tr.cir.q(x, tr.epar)}
    for m in ('G', 'C', 'i', 'q'):
        assert np.asarray(got[m]).tobytes() == np.asarray(want[m], float).tobytes(), m
    for r, c in lone:
        assert got['C'][r, c] == 0.0 and not np.signbit(got['C'][r, c]), (r, c)


def test_a_cached_capacitance_is_used_as_python_uses_it():
    """`_C_at_state` serves `C` from the memo or the converged-point cache
    before evaluating; the core asks the same lookup and skips the pass."""
    tr, x = _tr(mos_chain(6))
    fake = np.full((tr.cir.n, tr.cir.n), 1e-15)
    for on in (True, False):
        tr._C_cache = (x, fake)
        with core(on):
            _f, J = tr._residual_and_jacobian(x, 1.1e-7)
        assert tr._Cmat is fake
        if on:
            Jc = J
    assert Jc.tobytes() == J.tobytes()


def test_what_the_core_declines_takes_the_python_path():
    tr, x = _tr(mos_chain(6))
    t = 1.1e-7
    def ok():
        return _tran_core.evaluate(tr, x, t, None, 'fj') is not None
    assert ok()
    with core(False):
        assert not ok()
    assert _tran_core.evaluate(tr, x.astype(np.float32), t, None, 'fj') is None
    nan = x.copy()
    nan[2] = np.nan
    assert _tran_core.evaluate(tr, nan, t, None, 'fj') is None
    assert _tran_core.evaluate(tr, x[:-1], t, None, 'fj') is None
    ## a temperature that is not one number
    tr.epar = types.SimpleNamespace(T=np.array([300.0, 310.0]), bypasstol=-1.0)
    assert not ok()
    tr.epar = types.SimpleNamespace(T=np.array(310.0), bypasstol=-1.0)
    assert ok()                                 # (a 0-d array is one number)
    ## an integrator outside the enum
    base = tr.base_integrator
    tr.base_integrator = types.SimpleNamespace(
        check_order_drop=lambda h, hl, first: types.SimpleNamespace())
    assert not ok()
    tr.base_integrator = base
    ## an instance shadow of one of the four methods, then its removal
    el = tr.cir['M2']
    el.__dict__['q'] = lambda xx, epar=None, params_tree=None: np.zeros(el.n)
    assert not ok()
    del el.__dict__['q']
    assert ok()
    ## the switch off, and the bytes of the path it falls to
    with core(True):
        a = tr._residual_and_jacobian(x, t)
    with core(False):
        b = tr._residual_and_jacobian(x, t)
    assert a[1].tobytes() == b[1].tobytes()


def test_an_element_the_core_cannot_evaluate_leaves_the_circuit_to_python():
    cir = mos_chain(6)
    cir['D'] = Diode('d3', gnd)
    tr = Transient(cir, toolkit=cm.numeric)
    tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    assert _tran_core.core_for(tr) is None
    ## and the record is kept until the plan rebuilds
    assert tr.__dict__['_tran_core'][1] is None
    got = {}
    for on in (True, False):
        c2 = mos_chain(6)
        c2['D'] = Diode('d3', gnd)
        with core(on):
            res = Transient(c2, toolkit=cm.numeric).solve(
                tend=20 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = np.asarray(res.x, float)
    assert got[True].tobytes() == got[False].tobytes()


def test_a_constant_only_circuit_is_served_with_no_kernel_at_all():
    """A circuit of hand-written constant elements has no kernel to call
    and still goes through the core: the templates, the constant
    products and the companion in C -- the parity guard's RC ladder runs
    twice as fast (401 -> 214 ms for its pairs), bit for bit."""
    from pycircuit.circuit.elements import C
    got = {}
    for on in (True, False):
        cm.default_toolkit = cm.numeric
        cir = SubCircuit()
        cir['vs'] = VSin('n0', gnd, va=2.0, freq=1e5)
        for k in range(8):
            cir[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
            cir[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-9)
        with core(on):
            tr = Transient(cir, toolkit=cm.numeric)
            res = tr.solve(tend=20e-6, timestep=1e-6)
        got[on] = np.asarray(res.x, float).tobytes()
        if on and _tran_core.STATUS == 'c':
            c_ = _tran_core.core_for(tr)
            assert c_ is not None and c_.nb == 0 and c_.ng > 0
    assert got[True] == got[False]


def _served_and_python(build, **solve):
    """`build()`'s transient with the core on and off: the two results'
    bytes, and whether the core served it (None: no C driver here)."""
    got, served = {}, None
    for on in (True, False):
        cm.default_toolkit = cm.numeric
        cir, kw = build()
        with core(on):
            tr = Transient(cir, toolkit=cm.numeric, **kw)
            res = tr.solve(**solve)
        got[on] = np.asarray(res.x, float).tobytes()
        if on and _tran_core.STATUS == 'c':
            served = _tran_core.core_for(tr)
    return got, served


def test_an_eager_nonlinear_class_is_served_and_keeps_its_eager_forms():
    """Roadmap 168: an HDL class with no `var()` whose G depends on the
    unknowns (`DiodeHdl`) compiles through the let-chain too, so the core
    serves its circuit -- before, one such element took the WHOLE circuit
    to Python -- bit for bit with the Python path; it keeps its eager forms
    (`eval_i_pure` for `solve_batched` and the JAX path, the symbolic
    ones), and a class with constant G and C (`RHdl`) stays eager."""
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.elements import C
    info = eh.DiodeHdl._hdl_info
    assert info['chained'] and info['auto_chain']
    assert info['pure_spec'] is not None and info['sym_spec'] is not None
    assert hasattr(eh.DiodeHdl, 'eval_i_pure')
    assert not eh.RHdl._hdl_info['chained']
    assert not eh.RHdl._hdl_info['auto_chain']

    def build():
        cir = SubCircuit()
        cir['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
        cir['r'] = R('a', 'b', r=1e3)
        cir['d'] = eh.DiodeHdl('b', 'o')
        cir['c'] = C('o', gnd, c=1e-6)
        cir['rl'] = R('o', gnd, r=1e3)
        return cir, {}
    got, served = _served_and_python(build, tend=25 * 4e-5, timestep=4e-5,
                                     fixed_timestep=True)
    if _tran_core.STATUS == 'c':
        assert served is not None
    assert got[True] == got[False]


def test_a_class_with_dc_pins_is_served_and_declines_a_dc_epar():
    """Roadmap 168: a class with DC pins (`MemristorHdl`'s state) evaluates
    its pinned functions only in a DC analysis, so the core takes it --
    it refused the whole circuit before -- bit for bit with the Python
    path, and declines a DC analysis's epar ('pins_dc')."""
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.analysis import analysis_kind
    assert eh.MemristorHdl._hdl_info['state_meta']['dc_pins']

    def build():
        cir = SubCircuit()
        cir['vs'] = VSin('n1', gnd, va=1.0, freq=1.0)
        cir['m'] = eh.MemristorHdl('n1', gnd, ron=100.0, roff=16e3,
                                   muv=1e-14, d=10e-9, ic=0.1, p=1.0)
        return cir, {'uic': True}
    got, served = _served_and_python(build, tend=0.25, timestep=5e-3,
                                     fixed_timestep=True)
    assert got[True] == got[False]
    if _tran_core.STATUS != 'c':
        return
    assert served is not None and served.pins
    ep = cm.defaultepar.copy()
    assert type(served.probe(ep)) is float
    with analysis_kind(ep, 'dc'):
        assert served.probe(ep) == 'pins_dc'


@pytest.mark.parametrize('method', METHODS)
def test_a_transient_steps_bit_for_bit(method):
    got = {}
    for on in (True, False):
        cir = mos_chain(6)
        with core(on):
            tr = Transient(cir, toolkit=cm.numeric, integrator=INTEGRATOR[method]())
            res = tr.solve(tend=25 * 2e-8, timestep=2e-8, fixed_timestep=True)
        got[on] = (np.asarray(res.x, float).tobytes(), tr.statistics.newton_iterations,
                   _tran_core.core_for(tr) is not None, cir)
    assert got[True][0] == got[False][0] and got[True][1] == got[False][1]
    if _tran_core.STATUS == 'c':
        from pycircuit.circuit.tests.test_hdl_batch import _why
        assert got[True][2], _why(got[True][3], 'G')


def test_the_stage_pss_is_bit_for_bit():
    got = {}
    for on in (True, False):
        ## (the harness's PSS: its grid is coarse for the accuracy it asks,
        ## and says so; the comparison is the bytes, so the warning is quiet)
        with core(on), warnings.catch_warnings():
            warnings.simplefilter('ignore', simwarnings.AccuracyWarning)
            p = PSS(_stage(), method='gear', reltol=1e-8)
            p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
        got[on] = np.asarray(p.waveform[1], float).tobytes()
    assert got[True] == got[False]


def test_a_core_with_the_wrong_accumulation_order_is_caught(monkeypatch):
    tr, x = _tr(mos_chain(6))
    ## (every bin of three slots or more summed last to first, where the
    ## pull sums them in slot order, bincount's: two slots in either order
    ## are one sum)
    src = _tran_core.CORE_C
    for a, b in (('bin[wb[q]] = ((0.0 + V[ix[0]]) + V[ix[1]]) + V[ix[2]];',
                  'bin[wb[q]] = ((0.0 + V[ix[2]]) + V[ix[1]]) + V[ix[0]];'),
                 ('bin[wb[q]] = (((0.0 + V[ix[0]]) + V[ix[1]]) + V[ix[2]]) + V[ix[3]];',
                  'bin[wb[q]] = (((0.0 + V[ix[3]]) + V[ix[2]]) + V[ix[1]]) + V[ix[0]];'),
                 ('for (k = 0; k < cnt; k++) t += V[ix[k]];',
                  'for (k = cnt - 1; k >= 0; k--) t += V[ix[k]];')):
        assert src.count(a) == 1, a
        src = src.replace(a, b)
    assert src != _tran_core.CORE_C
    ffi, cfn, key, _cold, _secs = _cbuild.load_unit(src, _tran_core.CORE_CDEF, 'pcx_core')
    assert key != _cbuild.source_key(_tran_core.CORE_C, prelude=False)
    drv = _tran_core.driver()
    monkeypatch.setattr(_tran_core, '_driver', (ffi, cfn, drv[2]))
    tr.__dict__.pop('_tran_core', None)
    with core(True):
        a = tr._residual_and_jacobian(x, 1.1e-7)
    with core(False):
        b = tr._residual_and_jacobian(x, 1.1e-7)
    assert a[1].tobytes() != b[1].tobytes()


def test_the_switch_is_read_from_the_environment():
    env = dict(os.environ, PYCIRCUIT_TRAN_CORE='0')
    code = 'from pycircuit.circuit import _tran_core; print(_tran_core.CORE)'
    out = subprocess.run([sys.executable, '-c', code], env=env,
                         capture_output=True, text=True, check=True)
    assert out.stdout.strip() == 'False'


def test_an_element_is_packed_by_its_own_class_s_kernel(monkeypatch):
    """`_Core.probe` packs an element whose pack was dropped (a parameter
    change drops it) with the kernel of ITS class, as its batch does: two
    classes whose layouts put the temperature at another slot would
    otherwise have one's parameters packed in the other's places, silently."""
    from pycircuit.circuit import _hdl_cbackend as cb
    from pycircuit.circuit import elements_hdl as eh
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VSin('g', gnd, v=0.9, va=0.3, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['M'] = eh.MosLevel1Hdl('d', 'g', gnd, gnd)
    c['D'] = eh.DiodeSpiceHdl('d', gnd)
    tr, _x = _tr(c)
    core = tr.__dict__['_tran_core'][1]
    kern_of = {}
    for bt in core.batches:
        kern_of.setdefault(bt.cls, set()).add(bt.kern)
    assert len(kern_of) == 2
    calls = []
    pack = cb.CKernel.pack

    def spy(self, el):
        calls.append((self, el.__class__))
        return pack(self, el)
    monkeypatch.setattr(cb.CKernel, 'pack', spy)
    for el in core.uniq_els:
        el.__dict__.pop('_hdl_cp', None)
    assert _tran_core.core_for(tr) is not None and core.probe(tr.epar) is not None
    assert {cls for _k, cls in calls} == set(kern_of)
    assert all(k in kern_of[cls] for k, cls in calls), calls


def test_a_fused_group_s_row_is_where_its_column_names_say():
    """core.h's `core_ztab` names the columns `_Core.__init__` writes and
    core.c reads: a fused group's elements, inputs, passes and each pass's
    output length, read by name from the row a real core built."""
    tr, _x = _tr(_stage())
    core = tr.__dict__['_tran_core'][1]
    if not core.zgroups:
        pytest.skip('no fused group (the fused kernel is off)')
    Z = _cbuild.constants(_tran_core.CORE_CDEF)
    rows = core.arrays['ztab'].reshape(-1, Z['ZT_WIDTH'])
    assert len(rows) == len(core.zgroups)
    for row, (_info, fk, els, mask, _znm) in zip(rows, core.zgroups):
        assert row[Z['ZT_NE']] == len(els) and row[Z['ZT_K']] == fk.nx
        assert row[Z['ZT_PASSES']] == mask
        for m, p in enumerate(_tran_core.PASSES):
            assert row[Z['ZT_G_LEN'] + 2 * m] == fk.sizes.get(p, 0)


## ---------------------------------------------------------------------------
## THE SMALL PRODUCTS (core.c's `small_gemv`): a constant group's 2 x 2 and
## 3 x 3 products by OpenBLAS's dgemv's own operation order, where the check
## made once a process (`_small_gemv_ok`) found this process's dgemv to make
## exactly that.

def _capacitive_chain():
    """`mos_chain(8)` with a capacitor at each drain: constant groups of 2 x 2
    (the capacitors' q, the loads) and 3 x 3 (the sources)."""
    from pycircuit.circuit.elements import C
    cir = mos_chain(8)
    for k in range(8):
        cir[f'c{k}'] = C(f'd{k}', gnd, c=1e-13 * (k + 1))
    return cir


def test_the_small_products_are_checked_once_and_set_on_the_core():
    if not (_tran_core._small_gemv_ok(2) and _tran_core._small_gemv_ok(3)):
        pytest.skip("this process's dgemv makes a small product otherwise: the core calls it")
    assert _tran_core._small_gemv_ok(2) is _tran_core._small_gemv_ok(2)
    tr, _x = _tr(_capacitive_chain())
    cs = _tran_core.core_for(tr).cs
    ks = {int(cs.g_k[g]) for g in range(cs.ng)}
    assert {2, 3} <= ks, ks
    assert cs.small_gemv == (1 << 2) | (1 << 3)


def test_the_core_s_small_products_are_dgemv_s_bits():
    """The core with its small products made in C and with each a dgemv
    call: the residual and the Jacobian bit for bit, on states with signed
    zeros, huge and tiny values.  (A product's lone -0.0 reaches no output --
    the bins add each slot onto +0.0, as bincount does -- so the formula's
    final `0.0 +` is held by the check against dgemv, `_small_gemv_ok`, not
    here.)"""
    tr, x = _tr(_capacitive_chain())
    cs = _tran_core.core_for(tr).cs
    if not cs.small_gemv:
        pytest.skip("this process's dgemv makes a small product otherwise")
    states = _states(x, 8, seed=234) + [np.full(x.size, -0.0), -x, x * 1e150, x * 1e-160]
    bits = cs.small_gemv
    for xt in states:
        xt = np.ascontiguousarray(xt, dtype=float)
        got = []
        for on in (bits, 0):
            cs.small_gemv = on
            with np.errstate(all='ignore'):
                F, J = _tran_core.evaluate(tr, xt.copy(), 1.1e-7, None, 'fj')
            got.append((np.asarray(F).tobytes(), np.asarray(J).tobytes()))
        cs.small_gemv = bits
        assert got[0] == got[1]


@pytest.mark.parametrize('k', [2, 3])
def test_another_order_fails_the_check(monkeypatch, k):
    """The check sees an order of the terms other than dgemv's on its random
    products, and a lost `0.0 +` (a -0.0 kept) on its edges."""
    real = _tran_core._small_formula
    swapped = (lambda row, x: real(row[::-1], x[::-1]))

    def bare(row, x):
        t = float(np.float64(row[1]) * np.float64(x[1]))
        t = _tran_core._fma(row[0], x[0], t)
        if len(x) == 3:
            t = _tran_core._fma(row[2], x[2], t)
        return t
    for wrong in (swapped, bare):
        monkeypatch.setattr(_tran_core, '_SMALL_GEMV', {})
        monkeypatch.setattr(_tran_core, '_small_formula', wrong)
        assert _tran_core._small_gemv_ok(k) is False


def _core_once(core, x, T, bits, formula, Cin, Geq, J):
    """One call of the core entry as `evaluate` makes it, on fixed inputs:
    its C, q, iq and F (the Geq and J passed in written in place)."""
    ffi, cfn, _dgemv = _tran_core.driver()
    fb, d = ffi.from_buffer, core.dptr
    n = core.n
    rng = np.random.default_rng(20261104)
    q1, q2, iq1, u = (rng.standard_normal(n) for _ in range(4))
    C, q, iq, F = np.zeros((n, n)), np.empty(n), np.empty(n), np.empty(n)
    cfn(core.cs, fb(d, x), T, bits, formula, 1.5e9, -2.0e9, 0.5e9, 1e-9, 0.6,
        fb(d, q1), fb(d, q2), fb(d, iq1), fb(d, u), fb(d, Cin), fb(d, C), fb(d, q), fb(d, iq),
        fb(d, Geq), fb(d, F), fb(d, J))
    return C, q, iq, F


@pytest.mark.parametrize('formula', ['GEAR2', 'EULER', 'TRAPEZOIDAL', 'THETA'])
def test_the_packed_j_is_the_dense_one_on_the_pattern(formula):
    """PCX_WANT_PACKED: J's pattern entries in the pattern's order, each the
    dense call's bits, for every formula and both C sources (made, held);
    Geq not touched (a planted value stays); C, q, iq and F the dense
    call's -- and a residual-only call leaves Geq and J alone."""
    from pycircuit.circuit.tests.test_newton_c import _sparse_chain
    tr = Transient(_sparse_chain(), toolkit=cm.numeric, integrator=Gear2Integrator())
    tr.solve(tend=2e-9, timestep=1e-9)
    core = _tran_core.core_for(tr)
    T = core.ready(tr.epar)
    assert core.npat and T is not None
    n, pat = core.n, core.arrays['pat']
    abi = _tran_core._ABI
    f = abi['PCX_FORMULA_' + formula]
    x = np.random.default_rng(5).normal(0.0, 0.3, n)
    held = np.zeros((n, n))
    held.ravel()[pat] = np.random.default_rng(6).standard_normal(len(pat))
    pt, pk = abi['PCX_WANT_PATTERN'], abi['PCX_WANT_PACKED']
    for bits, Cin in ((abi['PCX_WANT_ALL'], held), (abi['PCX_PASS_G'] | abi['PCX_PASS_I'], held)):
        Gd, Jd = np.zeros((n, n)), np.zeros((n, n))
        dense = _core_once(core, x, T, bits | pt, f, Cin, Gd, Jd)
        Gp, Jp = np.full((n, n), 7.0), np.empty(len(pat))
        packed = _core_once(core, x, T, bits | pt | pk, f, Cin, Gp, Jp)
        assert all(a.tobytes() == b.tobytes() for a, b in zip(dense, packed))
        assert Jp.tobytes() == Jd.ravel()[pat].tobytes()
        assert (Gp == 7.0).all()
    Gr, Jr = np.full((n, n), 7.0), np.full(len(pat), 7.0)
    _core_once(core, x, T, abi['PCX_PASS_I'] | pt | pk, f, held, Gr, Jr)
    assert (Gr == 7.0).all() and (Jr == 7.0).all()
