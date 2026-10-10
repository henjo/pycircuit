"""The C units' interfaces agree with each other (the C stepping loop plan's
stage R1a, 2026-10-10).

The C code is compiled from Python strings and loaded through cffi's ABI
mode, which lays a struct out from the cdef TEXT it is handed and never
sees the compiled struct: a field added to a unit's C but not to its cdef
(or a type changed in one of the two) loads without complaint and puts
every later field at the wrong offset.  Likewise a unit that calls another
through a function pointer re-types the callee's signature by hand.  No
test compared either until this one.

Each struct is parsed twice by cffi's own parser -- from the cdef, and from
the typedef in the C text -- and compared field by field (name, offset,
size, kind) and by its size.  Each re-typed function-pointer copy is
compared with the callee's own definition by its ABI: the number of
arguments, each one's kind (pointer, integer, floating) and size, and the
result's.  The comparison is proved alive on planted differences.
"""
import re

import cffi
import pytest

from pycircuit.circuit import (
    _hdl_climit,
    _tran_core,
    _tran_lte_c,
    _tran_newton_c,
    _tran_radau_c,
    _tran_radau_tc,
)
from pycircuit.circuit.shooting import _sens_c

#: (the unit, its C text, its cdef, the struct both declare)
STRUCTS = [
    ('_tran_newton_c', _tran_newton_c.NEWTON_C, _tran_newton_c.NEWTON_CDEF, 'newton_t'),
    ('_tran_lte_c', _tran_lte_c.LTE_C, _tran_lte_c.LTE_CDEF, 'lte_t'),
    ('_tran_core', _tran_core.CORE_C, _tran_core.CORE_CDEF, 'hdl_core_t'),
    ('_tran_radau_c', _tran_radau_c.RADAU_C, _tran_radau_c.RADAU_CDEF, 'radau_t'),
    ('_tran_radau_tc', _tran_radau_tc.RADAU_TC, _tran_radau_tc.RADAU_TC_CDEF, 'radau_tf_t'),
    ('_tran_radau_tc.fold', _tran_radau_tc.FOLD_C, _tran_radau_tc.FOLD_CDEF, 'radau_fold_t'),
    ('_sens_c', _sens_c.SENS_C_SRC, _sens_c.SENS_CDEF, 'sens_t'),
]

#: the units that call the evaluation core and the limiter walk through a
#: pointer, each with its own copy of their signatures
CALLERS = [('_tran_newton_c', _tran_newton_c.NEWTON_C),
           ('_tran_radau_c', _tran_radau_c.RADAU_C),
           ('_tran_radau_tc', _tran_radau_tc.RADAU_TC)]


def _strip(c):
    """The C text without comments and preprocessor lines."""
    c = re.sub(r'/\*.*?\*/', ' ', c, flags=re.DOTALL)
    c = re.sub(r'//[^\n]*', ' ', c)
    return '\n'.join(ln for ln in c.split('\n') if not ln.lstrip().startswith('#'))


def _typedefs(c):
    """Every top-level `typedef ...;` statement of a C text, in order."""
    c = _strip(c)
    out, depth, start = [], 0, None
    for m in re.finditer(r'\btypedef\b|[{};]', c):
        tok = m.group()
        if tok == 'typedef' and depth == 0 and start is None:
            start = m.start()
        elif tok == '{':
            depth += 1
        elif tok == '}':
            depth -= 1
        elif tok == ';' and depth == 0 and start is not None:
            out.append(c[start:m.end()])
            start = None
    return out


def _definition(c, ret):
    """The `ret hdl_fn(...)` definition of a C text, as a typedef of a
    pointer to it named `_def_t`."""
    c = _strip(c)
    m = re.search(r'\b' + ret + r'\s+hdl_fn\s*\(([^)]*)\)\s*\{', c)
    assert m, 'no definition of hdl_fn in the C text'
    return f'typedef {ret} (*_def_t)({m.group(1)});'


def _kind(ffi, ct):
    if ct.kind == 'pointer':
        return ('ptr', 8)
    if ct.kind == 'primitive':
        return ('float' if ct.cname in ('double', 'float') else 'int', ffi.sizeof(ct))
    if ct.kind == 'array':
        return ('array', ct.length, _kind(ffi, ct.item))
    return (ct.kind, ffi.sizeof(ct))


def _layout(text, name):
    """`name`'s fields as cffi lays them out from `text`: (name, offset,
    size, kind) each, then the struct's size."""
    ffi = cffi.FFI()
    ffi.cdef(text)
    t = ffi.typeof(name)
    return [(f, fld.offset, ffi.sizeof(fld.type), _kind(ffi, fld.type)) for f, fld in t.fields] \
        + [('<sizeof>', ffi.sizeof(t))]


def _signature(text, name):
    """A function-pointer type's ABI: (result kind, (argument kinds))."""
    ffi = cffi.FFI()
    ffi.cdef(text)
    t = ffi.typeof(name)
    t = t.item if t.kind == 'pointer' else t
    return _kind(ffi, t.result), tuple(_kind(ffi, a) for a in t.args), t.ellipsis


@pytest.mark.parametrize('unit,c,cdef,name', STRUCTS, ids=[s[0] for s in STRUCTS])
def test_each_struct_s_cdef_is_its_c_struct(unit, c, cdef, name):
    """The cdef cffi lays the struct out from is the struct the C compiles."""
    assert _layout(cdef, name) == _layout('\n'.join(_typedefs(c)), name), unit


@pytest.mark.parametrize('unit,c', CALLERS, ids=[u for u, _ in CALLERS])
def test_each_copy_of_the_core_s_signature_is_the_core_s(unit, c):
    """`core_fn_t` in each caller is the evaluation core's own `hdl_fn`, as
    defined in `CORE_C` and declared in `CORE_CDEF` (its first argument a
    `const hdl_core_t *` there, a `const void *` in the copies: both a
    pointer)."""
    defined = _signature(_tran_core.CORE_CDEF + _definition(_tran_core.CORE_C, 'long'), '_def_t')
    declared = _signature(_tran_core.CORE_CDEF.replace(' hdl_fn(', ' (*_decl_t)(')
                          .replace('long (*_decl_t)', 'typedef long (*_decl_t)'), '_decl_t')
    assert defined == declared
    assert _signature('\n'.join(_typedefs(c)), 'core_fn_t') == defined, unit


@pytest.mark.parametrize('unit,c', CALLERS, ids=[u for u, _ in CALLERS])
def test_each_copy_of_the_walk_s_signature_is_the_walk_s(unit, c):
    """`walk_fn_t` in each caller is the limiter walk's own `hdl_fn`, as
    defined in `_WALK_C` and declared in `_WALK_CDEF`."""
    defined = _signature(_definition(_hdl_climit._WALK_C, 'long'), '_def_t')
    declared = _signature('typedef ' + _hdl_climit._WALK_CDEF.replace(' hdl_fn(', ' (*_decl_t)('),
                          '_decl_t')
    assert defined == declared
    assert _signature('\n'.join(_typedefs(c)), 'walk_fn_t') == defined, unit


def test_the_comparisons_see_a_planted_difference():
    """Alive: a field the cdef lacks, two fields swapped, a width changed
    in one argument -- each is caught."""
    c = '\n'.join(_typedefs(_tran_lte_c.LTE_C))
    base = _layout(c, 'lte_t')
    cdef = _tran_lte_c.LTE_CDEF
    first = re.search(r'\{\s*([^;]+;)', cdef).group(1)
    assert _layout(cdef.replace(first, '', 1), 'lte_t') != base          # one field fewer
    assert _layout(cdef.replace('{', '{ long planted;', 1), 'lte_t') != base
    copy = '\n'.join(_typedefs(_tran_newton_c.NEWTON_C))
    good = _signature(copy, 'core_fn_t')
    assert _signature(copy.replace('long want', 'int want', 1), 'core_fn_t') != good
    assert _signature(copy.replace('double T, long want', 'long want, double T', 1),
                      'core_fn_t') != good
