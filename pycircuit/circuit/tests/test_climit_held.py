"""THE LIMITER'S HELD STATEMENTS (`_hdl_climit`, THE HELD STATEMENTS; plan
261's L12b): a class limiter's parameter-chain statements no state reaches
are made once a parameter state and temperature (`lim_prep`) into the
limiter's own region of the element's pack -- after the class kernel's tail
-- and a call reads them as loads.  The region is the pack's: the class's
kernel, every view and the limiter carry the whole tail, a pack is made with
all of it clear, and a pack shorter than it is declined where the limiter
adopts it."""
import types

import numpy as np
import pytest

from pycircuit.circuit import _hdl_climit as cl
from pycircuit.circuit import _paths
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.tests.test_hdl_batch import mos_chain
from pycircuit.circuit.tests.test_hdl_climit import _bound, _inst, _py_limit
from pycircuit.circuit.tests.test_hdl_limit_walk import _both


def _layout(e):
    """`(kern, fk, n_p, base, nh)`: the class's limiter and kernel, the
    parameters' length, where the limiter's region starts, its held
    values' count."""
    info = _bound(e)
    kern, fk = info['_c_limit'], info['_c_fused']
    assert isinstance(kern, cl.CLimitKernel) and fk is not None
    n_p = kern.n_p
    ltail = kern.tail - fk.stats['tail']
    assert ltail > 2, 'a MOSFET limiter holds statements'
    return kern, fk, n_p, n_p + fk.stats['tail'], ltail - 2


def _limited_pairs(e, n):
    """States a Newton step from the accepted one whose limiter's answer
    moves with the temperature (a junction's and the fet law's parameters
    hold it)."""
    rng = np.random.default_rng(12)
    out = []
    while len(out) < 6:
        x0 = rng.uniform(-1.0, 2.0, n)
        x = x0 + rng.normal(0.0, 3.0, n)
        a = _py_limit(e, x, x0, types.SimpleNamespace(T=300.0))
        b = _py_limit(e, x, x0, types.SimpleNamespace(T=380.0))
        if a.tobytes() != b.tobytes():
            out.append((x, x0))
    return out


def test_the_class_pack_carries_the_limiter_s_region_clear():
    """Whichever kernel packs an element -- a pass's view or the limiter --
    the pack is the parameters, the class kernel's tail and the limiter's
    region, all clear (every stamp with it); the class kernel's `tail` is
    the whole, its own in its stats."""
    e = _inst(eh.MosLevel1Hdl)
    kern, fk, n_p, _base, nh = _layout(e)
    assert fk.tail == kern.tail == fk.stats['tail'] + nh + 2
    info = type(e)._hdl_info
    for k in [info['funcs'][m].__dict__['_hdl_c'] for m in ('G', 'q', 'i')] + [kern]:
        assert k.tail == fk.tail
        p, _pc = k.pack(e)
        assert p.size == n_p + fk.tail and not p[n_p:].any()


def test_the_held_values_are_read_and_made_again_at_another_temperature():
    """A call makes the region once (its stamp set, the temperature it was
    made at kept); the next calls READ it -- the region a call at 380 K made,
    planted under 300 K's bits, gives 380 K's answer at 300 K -- until the
    temperature's bits move, when it is made again and the answer is the
    closure's."""
    e = _inst(eh.MosLevel1Hdl)
    kern, _fk, _n_p, base, nh = _layout(e)
    t300, t380 = types.SimpleNamespace(T=300.0), types.SimpleNamespace(T=380.0)
    for x, x0 in _limited_pairs(e, e.n):
        ref300 = kern(e, x, x0, t300)
        p = e.__dict__['_hdl_cp'][0]
        assert p[base + nh] == 1.0 and p[base + nh + 1] == 300.0
        assert ref300.tobytes() == _py_limit(e, x, x0, t300).tobytes()
        ref380 = kern(e, x, x0, t380)
        assert p[base + nh + 1] == 380.0 and ref380.tobytes() != ref300.tobytes()
        made380 = p[base:base + nh].copy()
        kern(e, x, x0, t300)
        assert p[base:base + nh].tobytes() != made380.tobytes()
        ## (380 K's values under 300 K's bits: read, not made)
        p[base:base + nh] = made380
        assert kern(e, x, x0, t300).tobytes() == ref380.tobytes()
        ## (another temperature's bits: made again)
        kern(e, x, x0, types.SimpleNamespace(T=301.0))
        assert kern(e, x, x0, t300).tobytes() == ref300.tobytes()


def test_a_parameter_set_makes_a_new_pack_and_new_values():
    """A parameter set drops the pack; the next call packs anew (the region
    clear until it is made) and answers as the closure does with the new
    parameter."""
    e = _inst(eh.MosLevel1Hdl)
    kern, _fk, _n_p, base, nh = _layout(e)
    (x, x0), = _limited_pairs(e, e.n)[:1]
    kern(e, x, x0, defaultepar)
    old = e.__dict__['_hdl_cp'][0]
    before = old[base:base + nh].copy()
    e.ipar.vto = e.ipar.vto + 0.2
    e.update_iparv()
    assert '_hdl_cp' not in e.__dict__
    got = kern(e, x, x0, defaultepar)
    p = e.__dict__['_hdl_cp'][0]
    assert p is not old and p[base + nh] == 1.0
    assert p[base:base + nh].tobytes() != before.tobytes()
    assert got.tobytes() == _py_limit(e, x, x0, defaultepar).tobytes()


def _short_pack(e, fk, n_p):
    """The element's pack as a binding without the limiter's region made
    it: the class kernel's own tail alone -- a view at the head of a longer
    buffer whose rest (the guard, returned) a write past the pack would
    change, where it cannot corrupt the heap."""
    full, _pc = e.__dict__['_hdl_cp']
    m = n_p + fk.stats['tail']
    buf = np.full(full.size + 16, 7.0)
    buf[:m] = full[:m]
    p = buf[:m]
    e.__dict__['_hdl_cp'] = (p, fk.ffi.cast('double *', p.ctypes.data))
    return buf[m:]


def test_a_short_pack_is_declined_by_the_kernel():
    """A pack shorter than the limiter's tail -- made by another binding --
    is never written past: the kernel declines it (`climit:pack`) and the
    closure answers."""
    e = _inst(eh.MosLevel1Hdl)
    kern, fk, n_p, _base, _nh = _layout(e)
    (x, x0), = _limited_pairs(e, e.n)[:1]
    kern(e, x, x0, defaultepar)
    guard = _short_pack(e, fk, n_p)
    before = _paths.COUNTS['climit:pack']
    assert kern(e, x, x0, defaultepar) is None
    assert _paths.COUNTS['climit:pack'] == before + 1
    assert (guard == 7.0).all()


def test_a_short_pack_is_left_to_the_loop_by_the_walk():
    """On the walk an element whose pack is short is left to the loop
    (`walk:pack`), its answer the loop's -- the walk's and the loop's bytes
    the same."""
    cir = mos_chain(8)
    m3 = cir['M3']
    if not isinstance(type(m3)._hdl_info.get('_c_limit'), cl.CLimitKernel):
        pytest.skip('no limit kernel here')
    _kern, fk, n_p, _base, _nh = _layout(m3)
    rng = np.random.default_rng(3)
    x0 = rng.uniform(-1.0, 2.0, cir.n)
    x = x0 + rng.normal(0.0, 3.0, cir.n)
    _both(cir, x, x0)
    guard = _short_pack(m3, fk, n_p)
    before = _paths.COUNTS['walk:pack']
    _both(cir, x, x0)
    assert _paths.COUNTS['walk:pack'] > before and (guard == 7.0).all()


def test_the_c_newton_s_walk_declines_a_short_pack():
    """The C Newton's walk setup (`_tran_newton_c._walk_full`: the Radau
    C's and the transform loop's walk too) declines a short pack as
    `_take_elements` does -- behind the core, which declines it first, a
    second line: the setup is False and the walk's table never points at
    it."""
    from pycircuit.circuit import _tran_newton_c
    cir = mos_chain(8)
    m3 = cir['M3']
    if not isinstance(type(m3)._hdl_info.get('_c_limit'), cl.CLimitKernel):
        pytest.skip('no limit kernel here')
    drv = cl._walk_drv()
    if drv is None:
        pytest.skip('no walk driver here')
    _kern, fk, n_p, _base, _nh = _layout(m3)
    _tran_newton_c._mods()
    w = cl._walk_for(cir)
    assert _tran_newton_c._walk_full(w, drv[0])
    e = next(k for k, (inst, _el, _nm) in enumerate(w.entries) if inst == 'M3')
    _short_pack(m3, fk, n_p)
    short = m3.__dict__['_hdl_cp'][0]
    assert not _tran_newton_c._walk_full(w, drv[0])
    assert w.PR[e] != short.ctypes.data and w.mirror[e] is not m3.__dict__['_hdl_cp']
