"""`doc/transient_runtime.md`, the map of the transient's C runtime, names
every C file, every module that builds or loads C, and every switch those
modules read -- so a unit added without its row fails here."""
import glob
import os
import re

import pycircuit
from pycircuit.circuit import _cbuild

_PKG = os.path.dirname(os.path.abspath(pycircuit.__file__))
_MAP = os.path.join(os.path.dirname(_PKG), 'doc', 'transient_runtime.md')

#: (a module builds or loads C when it says one of these)
_LOADS = re.compile(r'_cbuild\.(?:csrc|load_unit|load_kernel)\(|\bload_(?:unit|kernel)\(')


def _map():
    with open(_MAP, encoding='utf-8') as fh:
        return fh.read()


def _glue():
    """Every module of the package (tests aside) that builds or loads C:
    `{name relative to pycircuit/circuit: its text}`."""
    out = {}
    for path in glob.glob(os.path.join(_PKG, '**', '*.py'), recursive=True):
        if os.sep + 'tests' + os.sep in path:
            continue
        with open(path, encoding='utf-8') as fh:
            text = fh.read()
        if _LOADS.search(text):
            out[os.path.splitext(os.path.relpath(path, os.path.join(_PKG, 'circuit')))[0]] = text
    return out


def _missing_files(text):
    files = sorted(os.path.basename(p) for p in glob.glob(os.path.join(_cbuild.CSRC_DIR, '*.c')))
    assert len(files) >= 15
    return [f for f in files if f'csrc/{f}' not in text]


def _missing_modules(text):
    glue = _glue()
    assert {'_tran_loop_c', '_tran_newton_c', '_cbuild'} <= set(glue)
    return sorted(m for m in glue if f'`{os.path.basename(m)}`' not in text and f'`{m}`' not in text)


def _missing_switches(text):
    switches = set()
    for src in _glue().values():
        switches |= set(re.findall(r'\bPYCIRCUIT_[A-Z0-9_]+', src))
    switches.discard('PYCIRCUIT_BENCH_PIN')   # (the benchmarks', not a unit's)
    assert len(switches) >= 15
    return sorted(s for s in switches if f'`{s}`' not in text)


def test_the_map_names_every_c_file():
    assert not _missing_files(_map()), 'not in doc/transient_runtime.md'


def test_the_map_names_every_module_that_builds_or_loads_c():
    assert not _missing_modules(_map()), 'not in doc/transient_runtime.md'


def test_the_map_names_every_switch_those_modules_read():
    assert not _missing_switches(_map()), 'not in doc/transient_runtime.md'


def test_the_map_s_checks_see_a_missing_row():
    """Alive: the map without the controller's row misses its C file, its
    module and its switch."""
    text = _map()
    for name in ('csrc/ctrl.c', '`_tran_ctrl_c`', '`PYCIRCUIT_CTRL_C`'):
        text = text.replace(name, '')
    assert _missing_files(text) == ['ctrl.c']
    assert _missing_modules(text) == ['_tran_ctrl_c']
    assert _missing_switches(text) == ['PYCIRCUIT_CTRL_C']


#: (history in a docstring or comment: a plan's stage or round, a date, the
#: logs -- `doc/transient_history.md` and the roadmap keep it)
_HISTORY = re.compile(r'\b(?:stage [0-9R]|stages [0-9]|speed round|round [0-9]+\b|'
                      r'20[0-9]{2}-[0-9]{2}-[0-9]{2}|roadmap|pss_log)', re.IGNORECASE)


def _history(glue):
    return {m: found for m, text in glue.items() if (found := _HISTORY.findall(text))}


def test_no_module_that_builds_or_loads_c_tells_history():
    """Their docstrings and comments state contracts: why a module exists
    and what it measured when it was built is the history's."""
    assert not _history(_glue())


def test_the_history_check_sees_a_planted_line():
    """Alive: a stage and a date in one module's comment, and nothing in a
    line naming the next stage of a solve."""
    glue = {'a': '## (the plan\'s stage 4; 2026-10-11)\n', 'b': '## (the next stage solve)\n'}
    assert _history(glue) == {'a': ['stage 4', '2026-10-11']}


#: (a code the C reads, said in the glue by its number: a companion, a pass
#: mask, the converged point's passes, the law's point mark)
_BARE_CODE = re.compile(r'formula (?:==|!=|in) [0-9(]|bits (?:&|&=|\|=) ~?[0-9]|'
                        r'fold_j (?:=|==|!=) \(?[0-9]|retry_point (?:=|==|!=) -?[0-9]|'
                        r'form\[0\] == [0-9]')


def test_no_module_that_builds_or_loads_c_says_a_code_by_its_number():
    """The codes are the headers' names (`_cbuild.constants`), on both sides
    of the boundary."""
    found = {m: _BARE_CODE.findall(text) for m, text in _glue().items()}
    assert not {m: f for m, f in found.items() if f}


def test_the_code_check_sees_a_planted_number():
    for text in ('if formula == 0:', 'J if bits & 1 else', 's.fold_j = 3 if j else 2',
                 'ct.retry_point != -2', "x if form[0] == 0 else y"):
        assert _BARE_CODE.findall(text), text
    assert not _BARE_CODE.findall('if formula == _GEAR2 and bits & _C:')
