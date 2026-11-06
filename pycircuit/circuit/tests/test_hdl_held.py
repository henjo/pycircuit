"""The held statements (`_hdl_cbackend`, The held statements): a class's
kernel makes the statements no `x` reaches once a parameter state into the
pack's tail and reads them as loads -- made again where the temperature
moves, never where it does not; a new parameter state is a new pack; the
flags their prep raised reach every call wanting their pass, as numpy's
warnings do, and no other.  (The kernel's bytes against numpy on drawn and
special states, every pass set: `test_twins_random.py`.)"""
import ctypes
import re
import types

import numpy as np

from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.tests.test_core_fused import _src
from pycircuit.circuit.tests.test_hdl_cbackend import _instance, c_backend, needs_cc

_LIBM = ctypes.CDLL('libm.so.6')
_LIBM.fetestexcept.argtypes = _LIBM.feclearexcept.argtypes = [ctypes.c_int]
#: (glibc x86-64: invalid 1, divide 4, overflow 8, underflow 16; inexact 32 not
#: among them)
FE_INVALID, FE_DIVBYZERO, FE_OVERFLOW, FE_UNDERFLOW = 1, 4, 8, 16
FLAGS = FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW | FE_UNDERFLOW
_NP = {'invalid': FE_INVALID, 'divide': FE_DIVBYZERO, 'overflow': FE_OVERFLOW,
       'underflow': FE_UNDERFLOW}


def test_the_held_statements_are_those_no_state_reaches():
    """`L_a` and `L_c` read parameters alone (`L_c` through `L_a`): held,
    made by the prep from the tail, loads in the kernel; `L_b` and `L_d`
    read `x`.  Without a layout nothing is held."""
    srcs = {'G': _src('  const double L_a = p[0]*2.0;', '  const double L_b = x[0]*L_a;',
                      '  const double L_c = L_a + p[1];', '  out[0] = L_b + L_c;'),
            'q': _src('  const double L_a = p[0]*2.0;', '  const double L_d = L_a*x[1];',
                      '  out[0] = L_d;')}
    text, stats = cb.fuse_csrc(srcs, (2, None))
    assert stats['held'] == 2 and stats['tail'] == 2 + cb.HELD_TAIL
    prep = text[text.index('void hdl_prep('):text.index('void hdl_fn(')]
    assert '  p[2] = p[0]*2.0;' in prep and '  p[3] = p[2] + p[1];' in prep
    body = text[text.index('void hdl_fn('):]
    assert '    L_a = p[2];' in body and '    L_c = p[3];' in body
    assert '    L_b = x[0]*L_a;' in body and '    L_d = L_a*x[1];' in body
    assert '  if (!hdl_held(p)) hdl_prep(p);' in body
    plain, pstats = cb.fuse_csrc(srcs)
    assert pstats['held'] == 0 and 'hdl_prep' not in plain


def test_the_held_flags_are_the_runtime_s():
    """The prep keeps PCX_FP_FLAGS -- the flags every unit tests -- by its
    own spelling of them."""
    from pycircuit.circuit import _cbuild
    fenv = _cbuild.csrc('pcx_fenv.h')
    unit = re.search(r'#define PCX_FP_FLAGS (\(.*\))', fenv).group(1)
    names = lambda s: sorted(re.findall(r'FE_\w+', s))
    assert names(unit) == names(cb.HELD_FLAGS)


def _mos(**kw):
    e = _instance(eh.MosLevel1Hdl, vto=0.5, gamma=0.4, lambd=0.02, cgso=2e-10, **kw)
    return e, type(e)


def _x(cls):
    x = np.full(cls._hdl_info['_c_fused'].nx, 0.3)
    x[0] = 1.0
    return x


def _tail(e, cls):
    """The element's pack, its held values and its stamp."""
    fk = cls._hdl_info['_c_fused']
    p = e.__dict__['_hdl_cp'][0]
    n_p, held = fk.layout[0], fk.stats['held']
    return p, p[n_p:n_p + held], p[n_p + held]


@needs_cc
def test_the_held_values_are_made_once_a_parameter_state():
    e, cls = _mos()
    with c_backend(cls):
        fk = cls._hdl_info['_c_fused']
        assert fk.stats['held'] > 0 and fk.layout[1] is not None, fk.stats
        x, epar = _x(cls), types.SimpleNamespace(T=300.0)
        ref = e.G(x, epar).copy()
        p, held, stamp = _tail(e, cls)
        assert p.size == fk.layout[0] + fk.tail and stamp == 1.0
        made = held.copy()
        ## the values are READ: a planted one is what the next call sees...
        held[:] = np.nan
        assert np.isnan(e.G(x, epar)).any()
        ## ...until the temperature moves (made again, at its bits)
        e.G(x, types.SimpleNamespace(T=301.0))
        assert e.G(x, epar).tobytes() == ref.tobytes()
        assert _tail(e, cls)[1].tobytes() == made.tobytes()
        ## a parameter set: a new pack, its stamp clear until a call
        e.ipar.lambd = 0.03
        e.update_iparv()
        assert '_hdl_cp' not in e.__dict__
        e.G(x, epar)
        assert _tail(e, cls)[0] is not p and _tail(e, cls)[2] == 1.0


@needs_cc
def test_a_pack_carries_its_class_kernel_s_tail():
    """Whichever kernel packs an element -- a pass's or the limiter's (same
    layout) -- the pack has its class's tail (the kernel's own, then the
    limiter's region: test_climit_held), the kernel's stamp clear."""
    e, cls = _mos()
    with c_backend(cls):
        fk = cls._hdl_info['_c_fused']
        n = fk.layout[0] + fk.tail
        for kern in [cls._hdl_info['funcs'][m].__dict__['_hdl_c'] for m in ('G', 'q')] + \
                [cls._hdl_info.get('_c_limit')]:
            if kern is None:
                continue
            p, _pc = kern.pack(e)
            assert p.size == n and p[fk.layout[0] + fk.stats['tail'] - cb.HELD_TAIL] == 0.0


def _numpy_flags(f, x, args):
    got = []
    old = np.seterrcall(lambda kind, flag: got.append(kind))
    try:
        with np.errstate(all='call'):
            f(x, *args)
    finally:
        np.seterrcall(old)
    return sum(v for k, v in _NP.items() if any(k in g for g in got))


@needs_cc
def test_a_held_statement_s_flag_reaches_every_call_wanting_its_pass():
    """`kp w` overflowing (a statement of the parameters alone, `G`'s and
    `i`'s): every call wanting `G` or `i` -- the one that makes the held
    values and each after it -- raises what numpy's pass raises, and a call
    wanting `C` or `q` alone what theirs do (nothing)."""
    e, cls = _mos(kp=1e300, w=1e300)
    with c_backend(cls):
        info = cls._hdl_info
        fk = info['_c_fused']
        x = _x(cls)
        args = hdl._args_of(e, defaultepar)
        want_np = {m: _numpy_flags(info['funcs'][m], x, args) for m in fk.sizes}
        assert want_np['G'] & FE_OVERFLOW and not want_np['C'], want_np
        ## (a fresh pack, its stamp clear: the first call makes the values)
        p, pcast = info['funcs']['G'].__dict__['_hdl_c'].pack(e)
        e.__dict__['_hdl_cp'] = (p, pcast)
        if fk.layout[1] is not None:
            p[fk.layout[1]] = 300.0
        ffi = fk.ffi
        xa = ffi.from_buffer('double *', np.ascontiguousarray(x))
        for rep in range(3):
            for b, m in enumerate(cb.FUSE_PASSES):
                if m not in fk.sizes:
                    continue
                out = np.empty(fk.sizes[m])
                o = ffi.from_buffer('double *', out)
                _LIBM.feclearexcept(FLAGS)
                fk.cfn(xa, pcast, 1 << b, o, o, o, o)
                got = _LIBM.fetestexcept(FLAGS)
                _LIBM.feclearexcept(FLAGS)
                assert got == want_np[m], (rep, m, got, want_np[m])
        assert _tail(e, cls)[2] == 1.0


def _short(v, e):
    """`e`'s pack as a kernel without the tail makes it (a stale or foreign
    one): the parameters alone."""
    k = cb.CKernel(v.ffi, v.cfn, v.shape, v.nx, (v.n_p, v.t_index), v.key, v.built_s, bit=v.bit)
    return k.pack(e)


@needs_cc
def test_a_short_pack_is_declined_by_the_per_element_call():
    """A pack shorter than the class's kernel writes its held statements
    into: `CKernel` declines it (numpy answers) and writes nothing."""
    e, cls = _mos()
    with c_backend(cls):
        v = cls._hdl_info['funcs']['G'].__dict__['_hdl_c']
        assert v.tail > 0
        p, pc = _short(v, e)
        guard = p.copy()
        e.__dict__['_hdl_cp'] = (p, pc)
        assert v(e, _x(cls), defaultepar) is None
        assert p.tobytes() == guard.tobytes()


@needs_cc
def test_a_short_pack_is_declined_by_the_core():
    from pycircuit.circuit.tests.test_watch import _on, _served
    _on()
    tr, core = _served()
    el = core.uniq_els[0]
    v = type(el)._hdl_info['funcs']['G'].__dict__['_hdl_c']
    assert v.tail > 0
    old = el.__dict__['_hdl_cp']
    el.__dict__['_hdl_cp'] = _short(v, el)
    try:
        assert core.probe(tr.epar) == 'pack'
    finally:
        el.__dict__['_hdl_cp'] = old
    assert core.probe(tr.epar).__class__ is not str
