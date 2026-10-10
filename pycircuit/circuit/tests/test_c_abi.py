"""The C units' interfaces agree with each other (the C stepping loop plan's
stages R1a and R1d, 2026-10-10).

The C units are loaded through cffi's ABI mode, which lays a struct out from
the cdef TEXT it is handed and never sees the compiled struct: a field in
one but not the other loads without complaint and puts every later field at
the wrong offset.  R1a compared the two texts; R1d made them ONE: each unit
`X.c` has a header `csrc/X.h` (its struct and its `hdl_fn` declaration),
the module's cdef IS that header's text, and `X.c` includes it, so the
compiler checks the definition against the very declaration cffi reads.
The function types the units share are declared once, in `pcx_abi.h`.
These tests hold that construction: each cdef is its header, each unit
includes its own, no C file declares a type, no header includes another.

The codes the C returns and is called with (a unit's results, the core's
passes and formulas) are named once, in the headers, and Python reads its
tables from them (`_cbuild.constants`): tested here against cffi's own
reading, against the pass order Python builds the fused masks from, and
for any bare number left where a name belongs.

What a header cannot check: `core_fn_t` and `walk_fn_t` (pcx_abi.h) re-type
the core's and the walk's entries for their callers, and the core's first
argument is a `const void *` there.  Those two are compared with their
definitions by ABI -- the number of arguments, each one's kind (pointer,
integer, floating) and size, and the result's -- and the comparison is
proved alive on planted differences.
"""
import glob
import os
import re

import cffi
import pytest

from pycircuit.circuit import (
    _cbuild,
    _hdl_batch,
    _hdl_cbackend,
    _hdl_climit,
    _tran_core,
    _tran_lte_c,
    _tran_newton_c,
    _tran_radau_c,
    _tran_radau_tc,
    _watch,
)
from pycircuit.circuit.shooting import _sens_c

#: (the unit, the cdef its module loads it through)
UNITS = [
    ('newton', _tran_newton_c.NEWTON_CDEF),
    ('lte', _tran_lte_c.LTE_CDEF),
    ('core', _tran_core.CORE_CDEF),
    ('radau', _tran_radau_c.RADAU_CDEF),
    ('radau_tc', _tran_radau_tc.RADAU_TC_CDEF),
    ('radau_fold', _tran_radau_tc.FOLD_CDEF),
    ('sens', _sens_c.SENS_CDEF),
    ('pass', _hdl_batch.PASS_CDEF),
    ('walk', _hdl_climit._WALK_CDEF),
    ('watch', _watch.WATCH_CDEF),
]


def _raw(name):
    with open(os.path.join(_cbuild.CSRC_DIR, name), encoding='utf-8', newline='') as fh:
        return fh.read()


def _strip(c):
    """The C text without comments and preprocessor lines."""
    c = re.sub(r'/\*.*?\*/', ' ', c, flags=re.DOTALL)
    c = re.sub(r'//[^\n]*', ' ', c)
    return '\n'.join(ln for ln in c.split('\n') if not ln.lstrip().startswith('#'))


def _declared_types(c):
    """The names every `typedef` of a C text declares."""
    c = _strip(c)
    return re.findall(r'\btypedef\b[^;{]*?(?:\{.*?\}\s*(\w+)|\(\*\s*(\w+)\s*\)[^;]*|(\w+))\s*;',
                      c, flags=re.DOTALL)


def _function_type(c, name):
    """The `typedef ... (*name)(...);` statement of a C text."""
    m = re.search(r'typedef[^;{]*?\(\*\s*' + name + r'\s*\)\s*\([^;]*\)\s*;', _strip(c))
    assert m, f'no function type {name}'
    return m.group(0)


def _definition(c, ret):
    """The `ret hdl_fn(...)` definition of a C text, as a typedef of a
    pointer to it named `_def_t`."""
    c = _strip(c)
    m = re.search(r'\b' + ret + r'\s+hdl_fn\s*\(([^)]*)\)\s*\{', c)
    assert m, 'no definition of hdl_fn in the C text'
    return f'typedef {ret} (*_def_t)({m.group(1)});'


def _declaration(h, ret):
    """The `ret hdl_fn(...);` declaration of a header, as `_decl_t`."""
    m = re.search(r'\b' + ret + r'\s+hdl_fn\s*\(([^)]*)\)\s*;', _strip(h))
    assert m, 'no declaration of hdl_fn in the header'
    return f'typedef {ret} (*_decl_t)({m.group(1)});'


def _kind(ffi, ct):
    if ct.kind == 'pointer':
        return ('ptr', 8)
    if ct.kind == 'primitive':
        return ('float' if ct.cname in ('double', 'float') else 'int', ffi.sizeof(ct))
    if ct.kind == 'array':
        return ('array', ct.length, _kind(ffi, ct.item))
    return (ct.kind, ffi.sizeof(ct))


def _signature(text, name):
    """A function-pointer type's ABI: (result kind, (argument kinds))."""
    ffi = cffi.FFI()
    ffi.cdef(text)
    t = ffi.typeof(name)
    t = t.item if t.kind == 'pointer' else t
    return _kind(ffi, t.result), tuple(_kind(ffi, a) for a in t.args), t.ellipsis


@pytest.mark.parametrize('unit,cdef', UNITS, ids=[u for u, _ in UNITS])
def test_each_unit_is_loaded_through_the_header_it_compiles_with(unit, cdef):
    """The cdef is the header's text; the unit includes the header once, and
    its compiled text holds it."""
    header = _raw(unit + '.h')
    assert cdef == header
    assert _raw(unit + '.c').count(f'#include "{unit}.h"\n') == 1
    assert _cbuild.csrc(unit + '.c').count(header) == 1


def test_no_c_file_declares_a_type_and_no_header_includes_another():
    """Every type lives in a header -- a struct in its unit's, a shared
    function type in pcx_abi.h, each declared once; one level of local
    includes, so `csrc` expands them all."""
    for path in glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.c')):
        assert not _declared_types(_raw(os.path.basename(path))), path
    seen = {}
    for path in glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.h')):
        h = _raw(os.path.basename(path))
        assert not re.search(r'^#include "', h, flags=re.MULTILINE), path
        for names in _declared_types(h):
            name = next(n for n in names if n)
            assert name not in seen, (name, path, seen.get(name))
            seen[name] = path
    for path in glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.c')):
        assert not re.search(r'^#include "', _cbuild.csrc(os.path.basename(path)),
                             flags=re.MULTILINE), path


def test_core_fn_t_is_the_core_s_entry():
    """`core_fn_t` (pcx_abi.h) is the core's own `hdl_fn`, as core.c defines
    it and core.h declares it (its first argument a `const hdl_core_t *`
    there, a `const void *` in the type: both a pointer)."""
    core_h = _raw('core.h')
    defined = _signature(core_h + _definition(_raw('core.c'), 'long'), '_def_t')
    declared = _signature(core_h + _declaration(core_h, 'long'), '_decl_t')
    assert defined == declared
    assert _signature(_function_type(_raw('pcx_abi.h'), 'core_fn_t'), 'core_fn_t') == defined


def test_walk_fn_t_is_the_walk_s_entry():
    """`walk_fn_t` (pcx_abi.h) is the limiter walk's own `hdl_fn`, as walk.c
    defines it and walk.h declares it."""
    defined = _signature(_definition(_raw('walk.c'), 'long'), '_def_t')
    declared = _signature(_declaration(_raw('walk.h'), 'long'), '_decl_t')
    assert defined == declared
    assert _signature(_function_type(_raw('pcx_abi.h'), 'walk_fn_t'), 'walk_fn_t') == defined


#: a code said by its number where the header names it: a bail, a status
#: set, a formula tested
_BARE_CODE = re.compile(r'_bail\(\w+, \d+\)|status = \d+;|formula [=!]= \d+')


@pytest.mark.parametrize('unit,cdef', UNITS, ids=[u for u, _ in UNITS])
def test_each_header_s_codes_read_as_cffi_reads_them(unit, cdef):
    """`_cbuild.constants` (how the modules read their tables) gives what
    cffi's own parser gives for every enum of the header."""
    ffi = cffi.FFI()
    ffi.cdef(cdef)
    by_cffi = {}
    for name in re.findall(r'\benum\s+(\w+)\s*\{', _strip(cdef)):
        by_cffi.update(ffi.typeof('enum ' + name).relements)
    assert _cbuild.constants(cdef) == by_cffi


def test_the_pass_bits_are_the_order_python_builds_the_masks_from():
    """pcx_pass's G, C, i, q bits are `1 << k` of the pass order the core's
    and the fused kernel's masks are built from in Python, and a want of
    every pass is G | C | i (q is always made)."""
    abi = _cbuild.constants(_raw('pcx_abi.h'))
    for order in (_tran_core.PASSES, _hdl_cbackend.FUSE_PASSES):
        assert [abi['PCX_PASS_' + p.upper()] for p in order] == [1 << k for k in range(4)]
    assert abi['PCX_WANT_ALL'] == abi['PCX_PASS_G'] | abi['PCX_PASS_C'] | abi['PCX_PASS_I']
    assert _tran_core.WANT['fj'] == abi['PCX_WANT_ALL']


def test_no_c_file_says_a_code_by_its_number():
    """A bail, a status set and a formula test say the header's name."""
    for path in glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.c')):
        assert not _BARE_CODE.findall(_strip(_raw(os.path.basename(path)))), path


def test_the_checks_see_a_planted_difference():
    """Alive: a width changed in one argument, two arguments swapped, a type
    declared in a C file, a header including another, a code said by its
    number, an enumerator spelling no value -- each is caught."""
    copy = _function_type(_raw('pcx_abi.h'), 'core_fn_t')
    good = _signature(copy, 'core_fn_t')
    assert _signature(copy.replace('long want', 'int want', 1), 'core_fn_t') != good
    assert _signature(copy.replace('double T, long want', 'long want, double T', 1),
                      'core_fn_t') != good
    assert _declared_types(_raw('newton.c') + 'typedef long (*x_t)(long);\n')
    assert _declared_types('typedef struct {\n    long a;\n} y_t;\n')
    assert re.search(r'^#include "', '#include "pcx_abi.h"\n', flags=re.MULTILINE)
    assert _BARE_CODE.findall(_raw('lte.c').replace('LTE_BAIL_DGESV', '4', 1))
    assert _BARE_CODE.findall('if (formula == 2) {')
    with pytest.raises(ValueError):
        _cbuild.constants('enum e { A = 1, B };')
