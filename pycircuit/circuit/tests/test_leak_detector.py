"""The leak detector (`pycircuit/_testing/leaks.py`, robust testing stage 2,
2026-10-03): every kind of planted leak is caught and named at the right
checkpoint, a restored change and a normal first-instance resolution are
not, and `fail` mode turns a leak into a failure.  Each scenario runs in a
pytester SUBPROCESS with the repository's root conftest (an in-process inner
run would share this worker's state, and would load pytest-randomly)."""
import json
import os

import pytest

pytest_plugins = ['pytester']

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


@pytest.fixture
def project(pytester, monkeypatch, tmp_path):
    ## (an inner run must not inherit the outer run's switches: a gate's
    ## state dump, a replay's timing switch, an explicit order)
    for var in ('PYCIRCUIT_TEST_TIMINGS', 'PYCIRCUIT_TEST_ORDER', 'PYCIRCUIT_STATE_DUMP',
                'PYCIRCUIT_LEAKS', 'PYCIRCUIT_LEAKS_REPORT', 'PYCIRCUIT_LEAKS_RUN'):
        monkeypatch.delenv(var, raising=False)
    with open(os.path.join(ROOT, 'conftest.py')) as f:
        pytester.makeconftest(f.read())
    report = tmp_path / 'leaks'
    monkeypatch.setenv('PYCIRCUIT_LEAKS_REPORT', str(report))
    monkeypatch.setenv('PYCIRCUIT_TEST_TIMINGS', '0')
    pytester.report = report
    return pytester


def _reports(project):
    out = []
    if project.report.is_dir():
        for p in sorted(project.report.iterdir()):
            out += [json.loads(ln) for ln in p.read_text().splitlines() if ln.strip()]
    return out


def _run(project, src, mode='warn', *args):
    was = os.environ.get('PYCIRCUIT_LEAKS')
    os.environ['PYCIRCUIT_LEAKS'] = mode
    try:
        project.makepyfile(test_planted=src)
        return project.runpytest_subprocess('-q', '-p', 'no:randomly', '-n', '0', *args)
    finally:
        if was is None:
            del os.environ['PYCIRCUIT_LEAKS']
        else:
            os.environ['PYCIRCUIT_LEAKS'] = was


def _changes(project, kind=None):
    return [c for r in _reports(project) if kind is None or r['kind'] == kind for c in r['changes']]


def test_a_switch_left_off_is_named_and_restored(project):
    _run(project, """
        from pycircuit.circuit import _tran_core
        def test_leaks():
            _tran_core.CORE = False
        def test_sees_it_restored():
            assert _tran_core.CORE is True
    """, 'fail').assert_outcomes(passed=2, errors=1)
    ch = _changes(project, 'test')
    assert any('pycircuit.circuit._tran_core.CORE: True -> False' in c for c in ch), ch


def test_a_restored_change_is_not_a_leak(project):
    _run(project, """
        from pycircuit.circuit import _tran_core
        def test_restores(monkeypatch):
            monkeypatch.setattr(_tran_core, 'CORE', False)
        def test_try_finally():
            was = _tran_core.CORE
            _tran_core.CORE = False
            try:
                pass
            finally:
                _tran_core.CORE = was
    """).assert_outcomes(passed=2)
    assert _changes(project) == []


def test_the_toolkit_env_path_and_numpy_state_are_caught_and_restored_in_fail(project):
    _run(project, """
        import os, sys
        import numpy as np
        from pycircuit.circuit import circuit
        def test_toolkit():
            circuit.default_toolkit = circuit.symbolic
        def test_env():
            os.environ['PYCIRCUIT_SOMETHING'] = '1'
        def test_path():
            sys.path.insert(0, '/nonexistent')
        def test_numpy():
            np.seterr(all='raise')
        def test_all_restored():
            assert circuit.default_toolkit is circuit.numeric
            assert 'PYCIRCUIT_SOMETHING' not in os.environ
            assert '/nonexistent' not in sys.path
            assert np.geterr()['divide'] != 'raise'
    """, 'fail').assert_outcomes(passed=5, errors=4)
    ch = ' | '.join(_changes(project, 'test'))
    for needle in ('toolkit:', 'env PYCIRCUIT_SOMETHING', 'sys.path:', 'np.geterr:'):
        assert needle in ch, (needle, ch)


def test_a_replaced_method_is_caught(project):
    _run(project, """
        from pycircuit.circuit import _tran_core
        def test_replaces():
            _tran_core.evaluate = lambda *a: None
    """).assert_outcomes(passed=1)
    assert any('_tran_core.evaluate replaced' in c for c in _changes(project, 'test'))


def test_a_module_fixture_leak_lands_at_the_module_boundary(project):
    _run(project, """
        import pytest
        from pycircuit.circuit import _stamp_plan
        @pytest.fixture(scope='module', autouse=True)
        def pins():
            _stamp_plan.ENABLED = False
            yield
        def test_one():
            pass
        def test_two():
            pass
    """).assert_outcomes(passed=2)
    assert _changes(project, 'test') == []
    assert any('_stamp_plan.ENABLED: True -> False' in c for c in _changes(project, 'module'))


def test_an_import_time_mutation_lands_at_collection(project):
    project.makepyfile(test_importer="""
        from pycircuit.circuit import _hdl_batch
        _hdl_batch.ENABLED = False
        def test_nothing():
            pass
    """)
    _run(project, """
        def test_other():
            pass
    """).assert_outcomes(passed=2)
    assert any('_hdl_batch.ENABLED: True -> False' in c for c in _changes(project, 'session'))


def test_hdl_backend_leaks_and_the_invariant(project):
    """A pin left in place, a stripped kernel, a swapped chain function --
    and a first instance resolving to C, which is normal."""
    _run(project, """
        import numpy as np
        from pycircuit.circuit import hdl
        from pycircuit.circuit import elements_hdl as eh
        from pycircuit.circuit.elements import gnd
        def test_first_instance_is_normal():
            eh.MosLevel1Hdl('d', 'g', gnd, gnd)
        def test_pin_left():
            hdl.set_backend('numpy', eh.DiodeSpiceHdl)
        def test_strip_a_kernel():
            v = type(eh.MosLevel1Hdl('d', 'g', gnd, gnd))
            v._hdl_info['funcs']['G'].__dict__.pop('_hdl_c', None)
        def test_swap_a_function():
            v = type(eh.MosLevel1Hdl('d', 'g', gnd, gnd))
            f = v._hdl_info['funcs']['q']
            v._hdl_info['funcs']['q'] = lambda *a: f(*a)
    """).assert_outcomes(passed=4)
    by = {r['where'].split('::')[-1]: ' | '.join(r['changes']) for r in _reports(project)}
    assert 'test_first_instance_is_normal' not in by, by
    assert 'pin' in by['test_pin_left'] or 'status' in by['test_pin_left'], by
    assert 'INVARIANT' in by['test_strip_a_kernel'] and 'carry no kernel' in by['test_strip_a_kernel']
    assert 'entries replaced' in by['test_swap_a_function']


def test_fail_mode_fails_the_leaking_test_only(project):
    r = _run(project, """
        from pycircuit.circuit import _tran_core
        def test_leaks():
            _tran_core.CORE = False
        def test_clean():
            assert _tran_core.CORE is True
    """, 'fail')
    r.assert_outcomes(passed=2, errors=1)
    assert 'LEAK' in r.stdout.str()


def test_warn_mode_observes_only_and_fail_mode_restores_values(project):
    """The first warn-mode gate's lesson: a lazily registered parameter is
    not a leak, a changed value is; warn leaves it changed, fail puts the
    value back -- and never the parameter table's structure."""
    _run(project, """
        from pycircuit.circuit import circuit
        from pycircuit.circuit.analysis import analysis_kind
        def test_registers():
            with analysis_kind(circuit.defaultepar, 'tran'):
                pass
        def test_moves_t():
            circuit.defaultepar.T = 350.0
        def test_still_moved_in_warn():
            assert circuit.defaultepar.T == 350.0
            assert circuit.defaultepar.analysis_kind is None
    """).assert_outcomes(passed=3)
    by = {r['where'].split('::')[-1]: r['changes'] for r in _reports(project)}
    assert 'test_registers' not in by, by
    assert any('defaultepar.T: 300.15 -> 350.0' in c for c in by['test_moves_t'])
    r = _run(project, """
        from pycircuit.circuit import circuit
        from pycircuit.circuit.analysis import analysis_kind
        def test_registers():
            with analysis_kind(circuit.defaultepar, 'tran'):
                pass
        def test_moves_t():
            circuit.defaultepar.T = 350.0
        def test_restored_in_fail():
            assert circuit.defaultepar.T == 300.15
            assert circuit.defaultepar.analysis_kind is None
    """, 'fail')
    r.assert_outcomes(passed=3, errors=1)


def test_off_mode_still_restores_the_toolkit(project):
    _run(project, """
        from pycircuit.circuit import circuit
        def test_toolkit():
            circuit.default_toolkit = circuit.symbolic
        def test_sees_numeric():
            assert circuit.default_toolkit is circuit.numeric
    """, 'off').assert_outcomes(passed=2)
    assert _changes(project) == []


def test_a_new_class_at_a_freed_address_is_a_new_class():
    """A class that dies frees its address, and a new class can be allocated
    there: keyed by `id`, the detector read a fixture's fresh test-local
    class as the old one with "functions replaced" (a shuffled run,
    2026-10-03).  The key is a serial held weakly per live class."""
    import gc

    from pycircuit._testing import state

    class Fresh:
        pass
    k1 = state.class_key(Fresh, True)
    addr = id(Fresh)
    del Fresh
    gc.collect()
    for _ in range(5000):
        class Fresh:
            pass
        if id(Fresh) == addr:
            break
    assert state.class_key(Fresh, True) != k1
