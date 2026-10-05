"""The fast-path counters (`_paths`; testing for development, stage 2,
2026-10-04).

Every fast path steps aside without changing an answer, so a path that
stopped serving passes every bit-identity check.  Each decline is counted
with its reason, each served call at the coarse sites, and the gate records
the counts per test (family `paths`).  Here: every `return None` in the
fast paths' entry functions IS a count (an AST check, with written
exemptions), the counts follow the paths on a real transient, and the gate's
comparison names a fast path that stopped serving.
"""
import ast
import inspect
import os
import textwrap

import pytest

from pycircuit.circuit import (
    _hdl_batch,
    _hdl_climit,
    _hdl_cse,
    _paths,
    _stamp_plan,
    _tran_companion,
    _tran_core,
    _tran_lte_c,
    _tran_newton_c,
    _tran_radau,
    _tran_radau_c,
    _tran_radau_tc,
)
from pycircuit.circuit import circuit as cm
from pycircuit.circuit.transient import Transient

pytest_plugins = ['pytester']

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

#: the fast paths' entry functions: every decline in them returns `_paths.no(...)`
DECLINING = [
    (_tran_core, 'evaluate'),
    (_tran_core, '_Core.ready'),
    (_tran_core, 'passes'),
    (_tran_newton_c, 'solve'),
    (_tran_radau_c, 'solve'),
    (_tran_radau, '_RadauStages._radau_frozen'),
    (_tran_radau_tc, 'solve'),
    (_tran_lte_c, 'max_error'),
    (_hdl_batch, 'Batch.run'),
    (_hdl_climit, 'CLimitKernel.__call__'),
    (_hdl_climit, 'limit_walk'),
    (_stamp_plan, 'assemble_matrix'),
    (_stamp_plan, 'assemble_vector'),
    (_stamp_plan, 'assemble_source'),
    (_hdl_cse, 'take'),
    (_tran_companion, '_CompanionModel._C_lookup'),
]
#: a `return None` that is not a decline says so on its line, and why
EXEMPT = '# paths: not a decline'


def _returns(fn):
    """The `return` statements of `fn` itself (not of a nested def)."""
    stack = list(fn.body)
    while stack:
        n = stack.pop()
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.Lambda, ast.ClassDef)):
            continue
        if isinstance(n, ast.Return):
            yield n
        stack.extend(ast.iter_child_nodes(n))


def uncounted(src):
    """The lines of a function's source whose `return` gives None without a
    count (and without the exemption)."""
    src = textwrap.dedent(src)
    lines = src.splitlines()
    fn = ast.parse(src).body[0]
    bad = []
    for r in _returns(fn):
        v = r.value
        if (v is None or (isinstance(v, ast.Constant) and v.value is None)) \
                and EXEMPT not in lines[r.lineno - 1]:
            bad.append(lines[r.lineno - 1].strip())
    return sorted(bad)


@pytest.mark.parametrize('mod, qual', DECLINING, ids=[q for _m, q in DECLINING])
def test_every_decline_is_counted(mod, qual):
    obj = mod
    for part in qual.split('.'):
        obj = getattr(obj, part)
    assert uncounted(inspect.getsource(obj)) == []


def test_the_check_sees_an_uncounted_decline():
    src = ("def f(x):\n"
           "    if x > 1:\n"
           "        return None\n"
           "    if x > 2:\n"
           "        return\n"
           "    if x > 3:\n"
           "        return None    # paths: not a decline (a lookup's miss)\n"
           "    def g():\n"
           "        return None\n"
           "    return _paths.no('k')\n")
    assert uncounted(src) == ['return', 'return None']


def test_no_counts_and_returns_none():
    before = _paths.snapshot()
    assert _paths.no('test.paths:probe') is None
    assert _paths.since(before) == {'test.paths:probe': 1}


## -- the counts on a real transient ---------------------------------------------------

def _core_or_skip():
    _tran_core.driver()                         # (loaded at first use)
    if _tran_core.STATUS != 'c':
        pytest.skip(f'the core is off: {_tran_core.STATUS}')


def _counted(cir):
    _core_or_skip()
    before = _paths.snapshot()
    Transient(cir, toolkit=cm.numeric).solve(tend=6 * 2e-8, timestep=2e-8,
                                             fixed_timestep=True)
    return _paths.since(before)


def test_a_served_transient_counts_its_core_and_no_decline():
    from pycircuit.circuit.tests.test_hdl_batch import mos_chain
    d = _counted(mos_chain(6))
    if d.get('once:core.build:built', 0) == 0 and d.get('core.fj:unservable'):
        pytest.skip('the chain is not served here (a class not C-bound)')
    assert d.get('core.fj:served', 0) > 0, d
    declined = {k: v for k, v in d.items()
                if k.startswith('core.') and not k.endswith(':served')}
    assert declined == {}, declined


def test_a_shadow_on_the_circuit_is_counted_where_the_core_stepped_aside():
    """An instance attribute shadowing `cir.G` (a test counting the calls, a
    harness timing them): the core declines every call -- the same answer,
    on the Python path -- and the counts say so."""
    from pycircuit.circuit.tests.test_hdl_batch import mos_chain
    cir = mos_chain(6)
    cir.G = cir.G
    d = _counted(cir)
    assert d.get('core.fj:served', 0) == 0
    assert d.get('core.fj:shadow_cir', 0) > 0, d


def test_a_switch_off_is_counted_by_its_reason(monkeypatch):
    from pycircuit.circuit.tests.test_hdl_batch import mos_chain
    monkeypatch.setattr(_tran_core, 'CORE', False)
    monkeypatch.setattr(_hdl_batch, 'ENABLED', False)
    d = _counted(mos_chain(6))
    assert d.get('core.fj:off', 0) > 0 and d.get('core.fj:served', 0) == 0, d
    assert d.get('batch.G:off', 0) > 0 and d.get('batch.G:served', 0) == 0, d


## -- the gate names it --------------------------------------------------------------------

def test_the_gate_comparison_names_a_fast_path_that_stopped_serving(pytester, monkeypatch, capsys):
    """End to end: the recorder's `paths` family on a run with and without a
    planted shadow; the comparison exits 4 (only the counts differ) and names
    the test and the counters that moved."""
    _core_or_skip()
    for var in ('PYCIRCUIT_TRAN_CORE', 'PYCIRCUIT_HDL_BATCH', 'PYCIRCUIT_STAMP_PLAN',
                'PYCIRCUIT_HDL_LIMIT_WALK', 'PYCIRCUIT_HDL_CLIMIT', 'PYCIRCUIT_HDL_FUSE',
                'PYCIRCUIT_HDL_ZERO_U', 'TRANREC_COVERAGE'):
        monkeypatch.delenv(var, raising=False)
    pytester.makepyfile(test_one="""
        import os
        from pycircuit.circuit import circuit as cm
        from pycircuit.circuit.transient import Transient
        from pycircuit.circuit.tests.test_hdl_batch import mos_chain

        def test_chain():
            cir = mos_chain(4)
            if os.environ.get('PLANT_SHADOW'):
                cir.G = cir.G
            Transient(cir, toolkit=cm.numeric).solve(tend=4e-8, timestep=2e-8,
                                                     fixed_timestep=True)
    """)
    monkeypatch.setenv('PYTHONPATH', os.pathsep.join(
        [ROOT, os.path.join(ROOT, 'benchmarks', 'tranrec'), os.environ.get('PYTHONPATH', '')]))
    monkeypatch.setenv('TRANREC_FAMILIES', 'paths')
    out = {}
    for label, plant in (('a', None), ('b', '1')):
        out[label] = str(pytester.path / label)
        monkeypatch.setenv('TRANREC_OUT', out[label])
        if plant:
            monkeypatch.setenv('PLANT_SHADOW', plant)
        r = pytester.runpytest_subprocess('-q', '-p', 'no:randomly', '-p', 'tran_recorder')
        r.assert_outcomes(passed=1)
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        '_tranrec_compare_paths', os.path.join(ROOT, 'benchmarks', 'tranrec', 'compare.py'))
    cmp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cmp)
    capsys.readouterr()
    assert cmp.main([out['a'], out['a']]) == 0
    assert cmp.main([out['a'], out['b']]) == 4
    said = capsys.readouterr().out
    assert "DIFF ('test_one.py::test_chain', 'paths', 0)" in said
    assert 'core.fj:shadow_cir 0 ->' in said and 'core.fj:served' in said
