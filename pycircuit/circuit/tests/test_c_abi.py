"""The C units' interfaces agree with each other (the C stepping loop plan's
stages R1a and R1d, 2026-10-10).

The C units are loaded through cffi's ABI mode, which lays a struct out from
the cdef TEXT it is handed and never sees the compiled struct: a field in
one but not the other loads without complaint and puts every later field at
the wrong offset.  R1a compared the two texts; R1d made them ONE: each unit
`X.c` has a header `csrc/X.h` (its struct and its entry `pcx_X`'s declaration),
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
import shutil
import subprocess

import cffi
import pytest

from pycircuit.circuit import (
    _cbuild,
    _hdl_batch,
    _hdl_cbackend,
    _hdl_climit,
    _tran_branch_c,
    _tran_core,
    _tran_ctrl_c,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton_c,
    _tran_pred_c,
    _tran_radau_c,
    _tran_radau_tc,
    _tran_sources,
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
    ('sources', _tran_sources.SOURCES_CDEF),
    ('branch', _tran_branch_c.BRANCH_CDEF),
    ('pred', _tran_pred_c.PRED_CDEF),
    ('ctrl', _tran_ctrl_c.CTRL_CDEF),
    ('loop', _tran_loop_c.LOOP_CDEF),
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


def _definition(c, ret, entry):
    """The `ret entry(...)` definition of a C text, as a typedef of a
    pointer to it named `_def_t`."""
    c = _strip(c)
    m = re.search(r'\b' + ret + r'\s+' + entry + r'\s*\(([^)]*)\)\s*\{', c)
    assert m, f'no definition of {entry} in the C text'
    return f'typedef {ret} (*_def_t)({m.group(1)});'


def _declaration(h, ret, entry):
    """The `ret entry(...);` declaration of a header, as `_decl_t`."""
    m = re.search(r'\b' + ret + r'\s+' + entry + r'\s*\(([^)]*)\)\s*;', _strip(h))
    assert m, f'no declaration of {entry} in the header'
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


@pytest.mark.parametrize('unit,cdef', UNITS, ids=[u for u, _ in UNITS])
def test_each_unit_compiles_on_its_own_without_a_warning(unit, cdef, tmp_path):
    """A unit is a complete C file (R1d step 3: built without the HDL
    kernels' prelude): what it includes and declares is all it needs, and
    the compiler finds nothing to warn about."""
    cc = _cbuild.find_compiler()[0]
    if cc is None:
        pytest.skip('no C compiler')
    path = tmp_path / (unit + '.c')
    path.write_text(_cbuild.csrc(unit + '.c'))
    proc = subprocess.run([cc, *_cbuild.CFLAGS, '-Wall', '-Wextra', '-Werror', '-fsyntax-only',
                           str(path)], capture_output=True, text=True, check=False)
    assert proc.returncode == 0, proc.stderr


@pytest.mark.parametrize('planted', ['newton.c', 'newton.h', 'pcx_abi.h', 'pcx_dense.h'])
def test_a_compiler_message_names_the_file_and_line_it_came_from(planted, tmp_path,
                                                                monkeypatch):
    """`csrc` marks each expanded header with `#line` (R1d step 4): an error
    planted at a line of a unit, or of a header it includes, is reported
    at that file and line -- and so is a sanitizer's report (`-g`)."""
    cc = _cbuild.find_compiler()[0]
    if cc is None:
        pytest.skip('no C compiler')
    d = tmp_path / 'csrc'
    shutil.copytree(_cbuild.CSRC_DIR, d)
    lines = (d / planted).read_text().split('\n')
    k = len(lines) - 2           # before the last line: code, not a comment
    lines.insert(k, '#error planted')
    (d / planted).write_text('\n'.join(lines))
    monkeypatch.setattr(_cbuild, 'CSRC_DIR', str(d))
    path = tmp_path / 'unit.c'
    path.write_text(_cbuild.csrc('newton.c'))
    proc = subprocess.run([cc, *_cbuild.CFLAGS, '-fsyntax-only', str(path)],
                          capture_output=True, text=True, check=False)
    assert f'{_cbuild.CSRC_NAME}{planted}:{k + 1}:' in proc.stderr, proc.stderr


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
    """`core_fn_t` (pcx_abi.h) is the core's own entry, as core.c defines
    it and core.h declares it (its first argument a `const hdl_core_t *`
    there, a `const void *` in the type: both a pointer)."""
    core_h = _raw('core.h')
    defined = _signature(core_h + _definition(_raw('core.c'), 'long', 'pcx_core'), '_def_t')
    declared = _signature(core_h + _declaration(core_h, 'long', 'pcx_core'), '_decl_t')
    assert defined == declared
    assert _signature(_function_type(_raw('pcx_abi.h'), 'core_fn_t'), 'core_fn_t') == defined


def test_walk_fn_t_is_the_walk_s_entry():
    """`walk_fn_t` (pcx_abi.h) is the limiter walk's own entry, as walk.c
    defines it and walk.h declares it."""
    defined = _signature(_definition(_raw('walk.c'), 'long', 'pcx_walk'), '_def_t')
    declared = _signature(_declaration(_raw('walk.h'), 'long', 'pcx_walk'), '_decl_t')
    assert defined == declared
    assert _signature(_function_type(_raw('pcx_abi.h'), 'walk_fn_t'), 'walk_fn_t') == defined


#: a code said by its number where the header names it: a bail, a status
#: set, a formula tested, a code returned (a helper's: 0 goes on, 1 is a
#: truth; every bail and status a header names here is 2 or more)
_BARE_CODE = re.compile(r'_bail\(\w+, \d+\)|status = \d+;|formula [=!]= \d+'
                        r'|return (?:[2-9]|\d{2,});')


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
    assert _BARE_CODE.findall(_raw('lte.c').replace('return LTE_BAIL_DGESV;', 'return 4;', 1))
    assert _BARE_CODE.findall(_raw('lte.c').replace('lte_bail(s, LTE_BAIL_NONFINITE)',
                                                    'lte_bail(s, 2)', 1))
    assert _BARE_CODE.findall('if (formula == 2) {')
    with pytest.raises(ValueError):
        _cbuild.constants('enum e { A = 1, B };')


def _declared(ffi):
    """Every type `ffi` declares: a typedef's C name and kind, and a
    struct's or union's -- named, or behind a typedef, as every header
    declares them -- size, alignment and fields (offset, type, bit
    width)."""
    typedefs, structs, unions = ffi.list_types()
    types = {f'struct {n}': ffi.typeof(f'struct {n}') for n in structs}
    types.update({f'union {n}': ffi.typeof(f'union {n}') for n in unions})
    types.update({n: ffi.typeof(n) for n in typedefs})
    out = {}
    for name, t in types.items():
        out[name] = (t.cname, t.kind)
        if t.kind in ('struct', 'union'):
            out[name] += (ffi.sizeof(t), ffi.alignof(t), tuple(
                (f, ft.offset, ft.type.cname, ft.bitsize) for f, ft in t.fields))
    return out


def _opened(ffi, path):
    """The unit's object opened through `ffi`: each name it declares --
    a function's C type, an enumerator's value."""
    lib = ffi.dlopen(path)
    out = {}
    for name in dir(lib):
        v = getattr(lib, name)
        out[name] = v if type(v) is int else ffi.typeof(v).cname
    return out


@pytest.mark.parametrize('unit,cdef', UNITS, ids=[u for u, _ in UNITS])
def test_each_unit_s_precompiled_ffi_declares_what_parsing_declares(unit, cdef):
    """The FFI a unit is loaded through (`_cbuild._ffi_for`: cffi's
    out-of-line module, read, not parsed) declares what parsing the header
    declares -- every type's layout, every function's signature, every
    code -- so which one a process got moves no byte."""
    pre = _cbuild._precompiled(cdef)
    assert pre is not None and type(pre).__module__ == '_cffi_backend'
    assert _cbuild._ffi_for(cdef) is pre or type(_cbuild._ffi_for(cdef)) is type(pre)
    parsed = cffi.FFI()
    parsed.cdef(cdef)
    assert _declared(pre) == _declared(parsed)
    path = _cbuild.load_unit(_cbuild.csrc(unit + '.c'), cdef, 'pcx_' + unit)[2]
    path = _cbuild.find_so(path)
    assert _opened(pre, path) == _opened(parsed, path)


def test_the_kernels_entry_s_precompiled_ffi_is_the_parsed_one():
    """The chain functions' entry, shared by every HDL kernel: its one
    signature as parsing declares it."""
    from pycircuit.circuit import hdl
    cdef = _cbuild._fn_cdef()
    pre = _cbuild._precompiled(cdef)
    parsed = cffi.FFI()
    parsed.cdef(cdef)
    assert f' {hdl._C_ENTRY}(' in cdef
    sig = 'void(*)(const double *, const double *, double *)'
    ## (cffi's C names drop `const`: the same name from both is the check)
    assert pre.typeof(sig).cname == parsed.typeof(sig).cname


def test_the_precompiled_check_sees_a_planted_difference():
    """Alive: one field's type changed, one enumerator's value moved --
    each is another declaration."""
    cdef = _tran_newton_c.NEWTON_CDEF
    pre = _cbuild._precompiled(cdef)
    planted = cffi.FFI()
    planted.cdef(re.sub(r'\bdouble (\w+);', r'float \1;', cdef, count=1))
    assert _declared(pre) != _declared(planted)
    m = re.search(r'(NEWTON_\w+) = (\d+)', cdef)
    assert m, 'no enumerator to plant on'
    moved = cffi.FFI()
    moved.cdef(cdef.replace(m.group(0), f'{m.group(1)} = {int(m.group(2)) + 100}', 1))
    path = _cbuild.find_so(_cbuild.load_unit(_cbuild.csrc('newton.c'), cdef, 'pcx_newton')[2])
    assert _opened(pre, path) != _opened(moved, path)


#: (history in a comment: a plan's stage or round, a date, the logs -- the
#: history is `doc/transient_history.md`'s and the roadmap's)
_HISTORY = re.compile(r'\b(?:stage [0-9R]|stages [0-9]|speed round|round [0-9]+\b|'
                      r'20[0-9]{2}-[0-9]{2}-[0-9]{2}|roadmap|pss_log)', re.IGNORECASE)


def _undocumented(text):
    """`(struct, line)` for every field line of `text`'s structs that no
    comment covers: none on its line, and no comment line heading its group
    (a blank line ends a group)."""
    out = []
    for m in re.finditer(r'typedef struct(?:\s+\w+)?\s*\{(.*?)\}\s*(\w+)\s*;', text, re.DOTALL):
        in_comment = covered = False
        for raw in m.group(1).split('\n'):
            line = raw.strip()
            if in_comment:
                if '*/' not in line:
                    continue
                in_comment = False
                line = line.split('*/', 1)[1].strip()
                if not line:
                    covered = True
                    continue
            if not line:
                covered = False
            elif line.startswith('/*'):
                in_comment = '*/' not in line
                covered = True
            elif ';' in line:
                has = '/*' in line
                in_comment = has and '*/' not in line.split('/*', 1)[1]
                if not (has or covered):
                    out.append((m.group(2), line))
    return out


def test_every_field_of_a_unit_s_struct_says_what_it_holds():
    """A field's comment, or its group's, says what it holds, who sets it
    and when it is valid."""
    for path in sorted(glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.h'))):
        assert not _undocumented(_raw(os.path.basename(path))), path


def test_no_c_file_tells_history():
    """Comments state contracts; why a unit exists and what it measured when
    it was built is the history's (`doc/transient_history.md`)."""
    for path in sorted(glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.[ch]'))):
        found = _HISTORY.findall(_raw(os.path.basename(path)))
        assert not found, (path, found)


def test_the_comment_checks_see_a_planted_gap():
    """Alive: a field with no comment, a comment telling a stage and a date."""
    assert _undocumented('typedef struct {\n    /* a */\n    long a;\n\n    long b;\n} t_t;\n') \
        == [('t_t', 'long b;')]
    assert not _undocumented('typedef struct {\n    long a;   /* a */\n    /* b,\n       c */\n'
                             '    long b;\n    long c;\n} t_t;\n')
    assert _HISTORY.findall('/* the plan\'s stage 4, 2026-10-11 */') == ['stage 4', '2026-10-11']
    assert not _HISTORY.findall('/* the next stage, the stage solve, a round-off */')


#: (a function a header declares: a line from column 0, not a typedef)
_FUNCTION = re.compile(r'^(?!typedef)[a-z][\w\s]*?\**\s*\b(\w+)\s*\(', re.MULTILINE)


@pytest.mark.parametrize('unit,cdef', UNITS, ids=[u for u, _ in UNITS])
def test_each_unit_s_header_declares_its_entry_and_no_other_function(unit, cdef):
    """A unit's entry is named for it, `pcx_<unit>`: the one function its
    header declares, what `load_unit` loads it by, what a backtrace or a
    sanitizer report names -- and what another unit calls it by (loop.c)."""
    assert set(_FUNCTION.findall(_strip(cdef))) == {'pcx_' + unit}


def test_the_entry_check_sees_a_planted_function():
    text = 'typedef long (*x_t)(long);\nlong pcx_a(a_t *s);\nvoid *pcx_b(int w);\n'
    assert set(_FUNCTION.findall(text)) == {'pcx_a', 'pcx_b'}
