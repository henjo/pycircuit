"""Every C twin against its Python reference on DRAWN inputs (hypothesis;
testing for development, stage 5, 2026-10-04).

The pinned pairs keep each C twin's source beside its reference's, and the
bit-identity tests compare them at fixed inputs (56 points per class, the
limiter's 2000 cases, the fast-path tests' states).  Here hypothesis draws
the inputs: states built from special values (signed zeros, subnormals,
huge, infinite, NaN), from the class's own parameter values give or take
an ulp (where its regions switch) and from ordinary values; parameters
perturbed; the temperature given the three ways a caller gives it.  Held to
their contracts (`test_hdl_cbackend._compare`): every chained library
device's C kernels against the numpy functions they were printed from (and
those against the uncompressed reference), the PSP model's, the limiter
kernel against the Python closure, and the evaluate core and the limiter
at the widths their C buffers hold and one past (`_wide_elements`).

The suite runs the `gate` profile (registered in the root conftest:
derandomized -- the same examples every run -- and no example database);
`--hypothesis-profile deep` draws thousands, fresh each run, keeping
failures in `~/.cache/pycircuit/hypothesis`.  A failure prints the
shrunk example.
"""
import math
import types
import warnings

import numpy as np
import pytest

hypothesis = pytest.importorskip('hypothesis')
from hypothesis import assume, given
from hypothesis import strategies as st

from pycircuit.circuit import _hdl_climit as cl
from pycircuit.circuit import _paths, _tran_core, hdl
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import R, SubCircuit, VSin, gnd
from pycircuit.circuit.tests import _wide_elements as wide
from pycircuit.circuit.tests.test_hdl_cbackend import (
    CHAINED,
    KW,
    _c_funcs,
    _compare,
    _instance,
    _libm_twin,
    _twin_exact,
    _ulp_names,
    c_backend,
    needs_cc,
)
from pycircuit.circuit.tests.test_hdl_climit import (
    _limiting_classes,
    _py_limit,
)

cm.default_toolkit = cm.numeric

SPECIALS = (0.0, -0.0, 5e-324, -5e-324, 2.2250738585072014e-308, 1e-300, -1e-300,
            1e-12, -1e-12, 0.5, -0.5, 1.0, -1.0, 40.0, -40.0, 700.0, -700.0,
            1e30, -1e30, 1e300, -1e300, math.inf, -math.inf, math.nan)
EPARS = (defaultepar, types.SimpleNamespace(T=350.0), types.SimpleNamespace(T=300),
         types.SimpleNamespace(T=np.float64(280.5)))


def _ulps(v, k):
    for _ in range(abs(k)):
        v = float(np.nextafter(v, math.inf if k > 0 else -math.inf))
    return v


def coordinate(params):
    """One state coordinate: ordinary, special, or a parameter value of the
    class give or take two ulps (a region boundary is often one)."""
    out = [st.floats(-3.0, 3.0), st.sampled_from(SPECIALS)]
    finite = sorted({float(p) for p in params if math.isfinite(float(p))})
    if finite:
        out.append(st.tuples(st.sampled_from(finite), st.integers(-2, 2))
                   .map(lambda t: _ulps(t[0], t[1])))
        out.append(st.sampled_from(finite).map(lambda v: -v))
    return st.one_of(out)


def states(n, params):
    return st.lists(coordinate(params), min_size=n, max_size=n).map(
        lambda v: np.ascontiguousarray(v, dtype=float))


## -- the device kernels ------------------------------------------------------------------

def _numeric(e):
    return [nm for nm in e.ipar._paramnames
            if isinstance(getattr(e.ipar, nm), (int, float))
            and not isinstance(getattr(e.ipar, nm), bool)]


def _ref(f, x, args):
    """`f(x, *args)` or None where Python raises (a parameter set the
    contract excludes: there the kernel makes inf or NaN)."""
    try:
        with np.errstate(all='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore')
            return np.asarray(f(x, *args), float)
    except (ArithmeticError, ValueError):
        return None


def _check_kernels(e, found, x, epar, args, ulp, strict=()):
    """Every C function of `found` at `x` against its numpy function (and
    that against the raw reference it was optimised from) -- run with the C
    library's functions for those of `ulp` (`_ulp_names`: numpy's own differ
    on this CPU, tanh here); False, checking nothing more, where Python
    raises (outside the contract)."""
    for fname, f in found:
        ref = _ref(f, x, args)
        if ref is None:
            return False
        raw = f.__dict__.get('_hdl_ref')
        if raw is not None:
            r0 = _ref(raw, x, args)
            ## the printed function is the optimised twin of the raw one:
            ## the same bytes, a NaN's (undefined) sign bit aside
            assert r0 is not None and _compare(r0, ref) in ('equal', 'nan-bits'), \
                (fname, x.tolist())
        out = f.__dict__['_hdl_c'](e, x, epar)
        assert out is not None, (fname, 'the kernel declined', epar)
        if ulp:
            ref = _ref(_libm_twin(f, ulp), x, args)
        got = _compare(ref, out)
        assert got != 'value', (fname, x.tolist(), epar)
        if fname in strict:
            assert got != 'zero-sign', (fname, x.tolist())
    return True


def one_at_a_time(n, params, bases):
    """Every special value and every parameter value (and one ulp either
    side, and its negative) at every coordinate of every base state.
    DRAWING finds a fault tied to one value at one coordinate rarely (a
    NaN masked at `x[0]`, a coordinate equal to a parameter: ~1 % an
    example; 60 examples missed both planted ones, 2026-10-04); this
    walks them all."""
    vals = list(SPECIALS)
    finite = sorted({float(q) for q in params if math.isfinite(float(q))})
    vals += sorted({_ulps(q, k) for q in finite for k in (-1, 0, 1)} | {-q for q in finite if q})
    for base in bases:
        for j in range(n):
            for v in vals:
                x = np.array(base, dtype=float)
                x[j] = v
                yield x


def _bases(n):
    return [np.zeros(n), np.full(n, 0.7), np.linspace(-0.9, 1.3, n)]


#: THE FUNCTIONS AN INSTANCE RUNS.  A MOSFET at its defaults (no `rd`,
#: `rs`) is a collapse variant (`..._collapse11`): its base class's
#: functions keep the internal nodes and divide by the zero resistances --
#: NaN in the drain current almost everywhere, so a comparison there sees
#: nothing (a planted fault in `i` passed every test, 2026-10-04).  So every
#: case tests `type(instance)`; the four MOSFETs also with series
#: resistances, where the uncollapsed functions run.
SERIES = {n: {'rd': 10.0, 'rs': 10.0} for n in ('MosLevel1Hdl', 'MosLevel1PmosHdl',
                                               'MosLevel3Hdl', 'MosLevel3PmosHdl',
                                               'MosLevel1GateChargeHdl',
                                               'MosLevel1PmosGateChargeHdl',
                                               'MosLevel3GateChargeHdl',
                                               'MosLevel3PmosGateChargeHdl')}
CASES = [(n, 'defaults') for n in CHAINED] + [(n, 'series') for n in SERIES]


def _device(name, how):
    kw = dict(KW.get(name, {}), **(SERIES[name] if how == 'series' else {}))
    e = _instance(getattr(eh, name), **kw)
    return e, type(e)


@needs_cc
@pytest.mark.parametrize('name, how', CASES, ids=[f'{n}-{h}' for n, h in CASES])
def test_every_c_kernel_takes_every_special_and_parameter_value_at_every_coordinate(name, how):
    e, cls = _device(name, how)
    n = len(hdl.x_layout(cls))
    args = [float(v) for v in hdl._args_of(e, defaultepar)]
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        found = _c_funcs(cls)
        assert found, 'no C kernels attached'
        ulp = _ulp_names(found)
        if not _twin_exact(ulp):
            pytest.skip("numpy's power is not the C library's pow on this CPU, and `**` "
                        'cannot be swapped in the twin')
        checked = sum(_check_kernels(e, found, x, defaultepar, args, ulp)
                      for x in one_at_a_time(n, args, _bases(n)))
        ## (not vacuous: the reference is a number in most entries -- NaN
        ## compares equal to NaN, and a fault there would pass unseen)
        for fname, f in found:
            share = np.mean([np.isfinite(r).mean() for x in one_at_a_time(n, args, _bases(n))
                             if (r := _ref(f, x, args)) is not None])
            assert share > 0.5, (fname, share)
    assert checked > 0


@needs_cc
@pytest.mark.parametrize('name, how', CASES, ids=[f'{n}-{h}' for n, h in CASES])
def test_every_c_kernel_answers_its_numpy_function_on_drawn_states(name, how):
    e, cls = _device(name, how)
    n = len(hdl.x_layout(cls))
    pnames = _numeric(e)
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        found = _c_funcs(cls)
        assert found, 'no C kernels attached'
        ulp = _ulp_names(found)
        if not _twin_exact(ulp):
            pytest.skip("numpy's power is not the C library's pow on this CPU, and `**` "
                        'cannot be swapped in the twin')

        @given(data=st.data())
        def drawn(data):
            epar = data.draw(st.sampled_from(EPARS), label='epar')
            ## a parameter perturbed (restored): its value moves the regions
            pick = data.draw(st.sampled_from([None] + pnames), label='parameter')
            old = getattr(e.ipar, pick) if pick is not None else None
            if pick is not None:
                factor = data.draw(st.sampled_from((0.5, 0.9, 1.0 + 2 ** -40, 1.1, 2.0)),
                                   label='factor')
                setattr(e.ipar, pick, old * factor)
                e.update_iparv()
            try:
                args = [float(v) for v in hdl._args_of(e, epar)]
                x = data.draw(states(n, args), label='x')
                assume(_check_kernels(e, found, x, epar, args, ulp))
            finally:
                if pick is not None:
                    setattr(e.ipar, pick, old)
                    e.update_iparv()
        drawn()


@needs_cc
def test_the_psp_kernels_answer_their_numpy_functions_on_drawn_biases():
    """PSP's stricter contract: no value difference anywhere, and on `i`,
    `G`, `q` not even a zero's sign (`TestPspBitIdentity`)."""
    from pycircuit.circuit import compact
    e = compact.PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                                  fnt=1.0)
    e.update_iparv()
    cls = type(e)
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        found = _c_funcs(cls)
        assert found, 'no C kernels attached'
        ulp = _ulp_names(found)
        if not _twin_exact(ulp):
            pytest.skip("numpy's power is not the C library's pow on this CPU, and `**` "
                        'cannot be swapped in the twin')

        @given(data=st.data())
        def drawn(data):
            bias = [data.draw(st.one_of(st.floats(lo, hi), st.sampled_from(SPECIALS)), label=nm)
                    for nm, lo, hi in (('vd', -0.3, 1.5), ('vg', -0.5, 1.8),
                                       ('vs', -0.2, 0.2), ('vb', -0.8, 0.2))]
            with np.errstate(all='ignore'):
                x = np.ascontiguousarray(e.bias(*bias), dtype=float)
            epar = data.draw(st.sampled_from(EPARS), label='epar')
            args = [float(v) for v in hdl._args_of(e, epar)]
            assume(_check_kernels(e, found, x, epar, args, ulp, strict=('i', 'G', 'q')))
        drawn()


## -- the fused kernel -----------------------------------------------------------------------

def _check_fused(e, fk, x, epar):
    """The class's fused kernel at `x` (`_hdl_cbackend.fuse_csrc`), every
    subset of its passes, against the passes' own kernels: the same bytes
    (a NaN's sign aside), and nothing written for a pass outside the
    subset.  False where a pass kernel declines the call."""
    from pycircuit.circuit import _hdl_cbackend as cb
    funcs = type(e)._hdl_info['funcs']
    ref = {}
    for m in fk.parts:
        r = funcs[m].__dict__['_hdl_c'](e, x, epar)
        if r is None:
            return False
        ref[m] = np.asarray(r, float).reshape(-1)
    _p, pcast = e.__dict__['_hdl_cp']
    ffi = fk.ffi
    dptr = ffi.typeof('double *')
    xa = np.ascontiguousarray(x, dtype=float)
    for want in range(1, 16):
        if want & ~fk.mask:
            continue
        outs = [np.full(max(fk.sizes.get(m, 0), 1), 1.25) for m in cb.FUSE_PASSES]
        fk.cfn(ffi.from_buffer(dptr, xa), pcast, want, *(ffi.from_buffer(dptr, o) for o in outs))
        for b, m in enumerate(cb.FUSE_PASSES):
            if m not in fk.parts:
                continue
            got = outs[b][:fk.sizes[m]]
            if want & (1 << b):
                assert _compare(ref[m], got) in ('equal', 'nan-bits'), (m, want, x.tolist())
            else:
                assert (got == 1.25).all(), (m, want, 'written outside its passes')
    return True


def _fused(cls):
    fk = cls._hdl_info.get('_c_fused')
    assert fk is not None, cls._hdl_info.get('_c_fused_status')
    return fk


#: AT ITS DEFAULTS A DEVICE CARRIES NO CHARGE: every capacitance zero, so a
#: fault in a charge statement multiplies out to the same zero (a planted one
#: passed, 2026-10-05).  These cases load the charges and the other zeros.
LOADED = {
    'MosLevel1Hdl': {'vto': 0.5, 'gamma': 0.4, 'lambd': 0.02, 'cgso': 2e-10, 'cgdo': 2e-10,
                     'cgbo': 1e-10, 'cbd': 1e-14, 'cbs': 1e-14, 'cj': 1e-4, 'cjsw': 1e-10,
                     'ad': 1e-12, 'asrc': 1e-12, 'pd': 4e-6, 'ps': 4e-6},
    ## the intrinsic gate charge (stage 6): TOX gives it its oxide, and a
    ## threshold and body effect put its regions inside the drawn states
    'MosLevel1GateChargeHdl': {'vto': 0.5, 'gamma': 0.4, 'tox': 2e-8, 'w': 1e-5, 'l': 1e-6,
                               'cbd': 1e-14, 'cbs': 1e-14},
    'MosLevel3GateChargeHdl': {'vto': 0.5, 'gamma': 0.4, 'tox': 2e-8, 'w': 1e-5, 'l': 1e-6,
                               'eta': 0.3, 'nfs': 1e11, 'vmax': 1e5, 'theta': 0.1},
    ## the substrate junction (stage 7), vertical and lateral
    'GummelPoonNpn4Hdl': {'rb': 100.0, 'rc': 2.0, 're': 1.0, 'cjs': 1e-12, 'vjs': 0.6,
                          'mjs': 0.4, 'cje': 1e-12, 'tf': 1e-10, 'vaf': 50.0,
                          'rbm': 10.0, 'irb': 1e-4, 'ptf': 30.0},
    'GummelPoonPnp4Hdl': {'rb': 100.0, 'rc': 2.0, 're': 1.0, 'cjs': 1e-12, 'vjs': 0.6,
                          'mjs': 0.4, 'subs': -1.0},
    'GummelPoonNpnHdl': dict(KW.get('GummelPoonNpnHdl', {}), vaf=50.0, ikf=0.1, ise=1e-15,
                             cje=1e-12, cjc=5e-13, tf=1e-10, tr=1e-8, xtf=1.0, vtf=2.0,
                             itf=0.1),
}
FUSED_CASES = CASES + [(n, 'loaded') for n in LOADED]


def _fused_device(name, how):
    if how != 'loaded':
        return _device(name, how)
    e = _instance(getattr(eh, name), **LOADED[name])
    return e, type(e)


@needs_cc
@pytest.mark.parametrize('name, how', FUSED_CASES, ids=[f'{n}-{h}' for n, h in FUSED_CASES])
def test_the_fused_kernel_takes_every_special_and_parameter_value_at_every_coordinate(name, how):
    from pycircuit.circuit import _hdl_cbackend as cb
    if not cb.FUSE:
        pytest.skip('the fused kernels are off (PYCIRCUIT_HDL_CFUSE=0)')
    e, cls = _fused_device(name, how)
    n = len(hdl.x_layout(cls))
    args = [float(v) for v in hdl._args_of(e, defaultepar)]
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        fk = _fused(cls)
        checked = sum(_check_fused(e, fk, x, defaultepar)
                      for x in one_at_a_time(n, args, _bases(n)))
    assert checked > 0


@needs_cc
@pytest.mark.parametrize('name, how', FUSED_CASES, ids=[f'{n}-{h}' for n, h in FUSED_CASES])
def test_the_fused_kernel_answers_the_passes_on_drawn_states(name, how):
    from pycircuit.circuit import _hdl_cbackend as cb
    if not cb.FUSE:
        pytest.skip('the fused kernels are off (PYCIRCUIT_HDL_CFUSE=0)')
    e, cls = _fused_device(name, how)
    n = len(hdl.x_layout(cls))
    pnames = _numeric(e)
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        fk = _fused(cls)

        @given(data=st.data())
        def drawn(data):
            epar = data.draw(st.sampled_from(EPARS), label='epar')
            pick = data.draw(st.sampled_from([None] + pnames), label='parameter')
            old = getattr(e.ipar, pick) if pick is not None else None
            if pick is not None:
                factor = data.draw(st.sampled_from((0.5, 0.9, 1.0 + 2 ** -40, 1.1, 2.0)),
                                   label='factor')
                setattr(e.ipar, pick, old * factor)
                e.update_iparv()
            try:
                args = [float(v) for v in hdl._args_of(e, epar)]
                x = data.draw(states(n, args), label='x')
                assume(_check_fused(e, fk, x, epar))
            finally:
                if pick is not None:
                    setattr(e.ipar, pick, old)
                    e.update_iparv()
        drawn()


@needs_cc
def test_the_psp_fused_kernel_answers_its_passes_on_drawn_biases():
    from pycircuit.circuit import _hdl_cbackend as cb
    from pycircuit.circuit import compact
    if not cb.FUSE:
        pytest.skip('the fused kernels are off (PYCIRCUIT_HDL_CFUSE=0)')
    e = compact.PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                                  fnt=1.0)
    e.update_iparv()
    cls = type(e)
    with c_backend(cls):
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        fk = _fused(cls)
        assert fk.stats['union'] < fk.stats['printed'] * 0.6, fk.stats

        @given(data=st.data())
        def drawn(data):
            bias = [data.draw(st.one_of(st.floats(lo, hi), st.sampled_from(SPECIALS)), label=nm)
                    for nm, lo, hi in (('vd', -0.3, 1.5), ('vg', -0.5, 1.8),
                                       ('vs', -0.2, 0.2), ('vb', -0.8, 0.2))]
            with np.errstate(all='ignore'):
                x = np.ascontiguousarray(e.bias(*bias), dtype=float)
            epar = data.draw(st.sampled_from(EPARS), label='epar')
            assume(_check_fused(e, fk, x, epar))
        drawn()


## -- the limiter kernel ---------------------------------------------------------------------

def _limiter_matches_on(e, kern, n, epar_strategy, params):
    @given(data=st.data())
    def drawn(data):
        epar = data.draw(epar_strategy, label='epar')
        x0 = data.draw(states(n, params), label='x0')
        x = data.draw(states(n, params), label='x')
        got = kern(e, x, x0, epar)
        if got is None:
            ## the kernel declines only where a sort key is not a number
            assert not (np.isfinite(x).all() and np.isfinite(x0).all()), (x, x0)
            return
        try:
            ref = _py_limit(e, x, x0, epar)
        except (ArithmeticError, ValueError):
            assume(False)
        assert got.tobytes() == ref.tobytes(), (x.tolist(), x0.tolist(), epar)
    drawn()


@needs_cc
@pytest.mark.parametrize('name', _limiting_classes())
def test_the_limiter_kernel_takes_every_special_value_at_every_coordinate(name):
    cls = getattr(eh, name)
    e = _instance(cls)
    info = type(e)._hdl_info
    if not info.get('_c_bound'):
        pytest.fail(f'{name} is not bound to C: {type(e)._hdl_backend_status}')
    kern = cl.bind(type(e), info)
    assert kern is not None, info.get('_c_limit_status')
    n = e.n
    params = [float(v) for v in hdl._args_of(e, defaultepar)]
    rng = np.random.default_rng(n)
    for x in one_at_a_time(n, params, _bases(n)):
        x0 = x + rng.normal(0.0, 0.5, n)
        for a, b in ((x, x0), (x0, x)):          # the special in the step and in the old state
            got = kern(e, a, b, defaultepar)
            if got is None:
                assert not (np.isfinite(a).all() and np.isfinite(b).all()), (a, b)
                continue
            assert got.tobytes() == _py_limit(e, a, b).tobytes(), (a.tolist(), b.tolist())


@needs_cc
@pytest.mark.parametrize('name', _limiting_classes())
def test_the_limiter_kernel_answers_the_closure_on_drawn_steps(name):
    cls = getattr(eh, name)
    e = _instance(cls)
    info = type(e)._hdl_info
    if not info.get('_c_bound'):
        pytest.fail(f'{name} is not bound to C: {type(e)._hdl_backend_status}')
    kern = cl.bind(type(e), info)
    assert kern is not None, info.get('_c_limit_status')
    params = [float(v) for v in hdl._args_of(e, defaultepar)]
    _limiter_matches_on(e, kern, e.n, st.sampled_from(EPARS), params)


## -- the C buffers' limits -----------------------------------------------------------------

def _ring_circuit(cls, k):
    """`cls` (k terminals) on k nodes, each through 1 k to ground, one
    driven by a sine."""
    c = SubCircuit()
    c['vs'] = VSin('n0', gnd, va=0.2, freq=1e6)
    c['W'] = cls(*[f'n{j}' for j in range(k)])
    for j in range(1, k):
        c[f'R{j}'] = R(f'n{j}', gnd, r=1e3)
    return c


def _core_tr(cir):
    from pycircuit.circuit.tests.test_tran_core import _tr
    return _tr(cir)


@needs_cc
def test_an_element_as_wide_as_the_core_holds_is_served_bit_for_bit():
    from pycircuit.circuit.tests.test_tran_core import _state, core
    cir = _ring_circuit(wide.Wide64, 64)
    tr, x = _core_tr(cir)
    assert _tran_core.core_for(tr) is not None
    t = 1.1e-7
    rng = np.random.default_rng(64)
    for j in range(6):
        xt = x if j == 0 else x + rng.normal(0.0, 0.05, x.size)
        before = _paths.snapshot()
        with core(True):
            a = tr._residual_and_jacobian(xt, t)
            sa = _state(tr)
        assert _paths.since(before).get('core.fj:served', 0) == 1
        with core(False):
            b = tr._residual_and_jacobian(xt, t)
            sb = _state(tr)
        assert a[0].tobytes() == b[0].tobytes() and a[1].tobytes() == b[1].tobytes(), j
        assert sa == sb, j


@needs_cc
def test_an_element_one_wider_than_the_core_holds_is_left_to_python():
    """65 unknowns: the core is not built (`once:core.build:unservable`) and
    every call declines as unservable -- the transient runs the Python
    path, the same bytes."""
    cir = _ring_circuit(wide.Wide65, 65)
    before = _paths.snapshot()
    from pycircuit.circuit.transient import Transient
    tr = Transient(cir, toolkit=cm.numeric)
    res = tr.solve(tend=4 * 2e-8, timestep=2e-8, fixed_timestep=True)
    d = _paths.since(before)
    assert _tran_core.core_for(tr) is None
    assert d.get('once:core.build:unservable', 0) >= 1 and d.get('core.fj:served', 0) == 0, d
    assert d.get('core.fj:unservable', 0) > 0, d
    from pycircuit.circuit.tests.test_tran_core import core
    with core(False):
        ref = Transient(_ring_circuit(wide.Wide65, 65), toolkit=cm.numeric).solve(
            tend=4 * 2e-8, timestep=2e-8, fixed_timestep=True)
    assert np.asarray(res.x).tobytes() == np.asarray(ref.x).tobytes()


@needs_cc
def test_a_limiter_as_wide_as_the_write_back_holds_is_the_kernels():
    """32 limited unknowns: the class binds its limiter kernel, which
    answers the closure on drawn steps; 33: the kernel refuses itself (the
    class stays on C) and the closure answers."""
    e32 = wide.LimWide32(*[f'n{j}' for j in range(32)])
    e32.update_iparv()
    info = type(e32)._hdl_info
    assert info.get('_c_bound'), type(e32)._hdl_backend_status
    assert info.get('_c_limit') is not None, info.get('_c_limit_status')
    _limiter_matches_on(e32, info['_c_limit'], 32, st.just(defaultepar), [1.0])
    e33 = wide.LimWide33(*[f'n{j}' for j in range(33)])
    e33.update_iparv()
    i33 = type(e33)._hdl_info
    assert i33.get('_c_bound'), type(e33)._hdl_backend_status
    assert i33.get('_c_limit') is None, i33.get('_c_limit_status')
