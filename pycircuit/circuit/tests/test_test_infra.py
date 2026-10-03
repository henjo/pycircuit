"""The test infrastructure itself (2026-10-03, robust testing stage 1): the
root conftest's order modes and run log, the per-process backend state dump,
and the replay and polluter scripts.  Every inner run is a SUBPROCESS
(`runpytest_subprocess`) with `-p no:randomly` unless the test is about
randomly: an in-process inner run would read no `pytest.ini`, load randomly
and reseed this worker's global RNG."""
import json
import os
import subprocess
import sys

import pytest

pytest_plugins = ['pytester']

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
SCRIPTS = os.path.join(ROOT, 'scripts')
sys.path.insert(0, SCRIPTS) if SCRIPTS not in sys.path else None
import test_replay  # scripts/test_replay.py

sys.path.remove(SCRIPTS)


@pytest.fixture
def project(pytester, monkeypatch):
    """A tiny project with the repository's root conftest, three tests and a
    full-suite timing record that orders them b, c, a."""
    ## (an inner run must not inherit the outer run's switches: a gate's
    ## state dump, a replay's timing switch, an explicit order)
    for var in ('PYCIRCUIT_TEST_TIMINGS', 'PYCIRCUIT_TEST_ORDER', 'PYCIRCUIT_STATE_DUMP',
                'PYCIRCUIT_LEAKS', 'PYCIRCUIT_LEAKS_REPORT', 'PYCIRCUIT_LEAKS_RUN'):
        monkeypatch.delenv(var, raising=False)
    with open(os.path.join(ROOT, 'conftest.py')) as f:
        pytester.makeconftest(f.read())
    pytester.makepyfile(test_three="""
        def test_a():
            pass
        def test_b():
            pass
        def test_c():
            pass
    """)
    rec = pytester.path / 'test_timings'
    rec.mkdir()
    (rec / '20000101T000000Z_x.json').write_text(json.dumps(
        {'tests': 1000, 'durations': [['test_three.py::test_b', 3.0],
                                      ['test_three.py::test_c', 2.0],
                                      ['test_three.py::test_a', 1.0]]}))
    return pytester


def _order(result):
    return [ln.split('::')[-1] for ln in result.outlines if ln.startswith('test_three.py::')]


def test_the_sort_is_longest_first_and_steps_aside_on_request(project, monkeypatch):
    r = project.runpytest_subprocess('--co', '-q', '-p', 'no:randomly')
    assert _order(r) == ['test_b', 'test_c', 'test_a']
    monkeypatch.setenv('PYCIRCUIT_TEST_ORDER', 'collected')
    r = project.runpytest_subprocess('--co', '-q', '-p', 'no:randomly')
    assert _order(r) == ['test_a', 'test_b', 'test_c']


def test_the_sort_yields_to_a_shuffled_run(project):
    """randomly shuffles before every plain hook; the sort must not undo it:
    over a few seeds the shuffled order differs from the sorted one."""
    orders = set()
    for seed in (1, 2, 3, 4, 5, 6):
        r = project.runpytest_subprocess('--co', '-q', '-p', 'randomly', f'--randomly-seed={seed}')
        orders.add(tuple(_order(r)))
    assert len(orders) > 1, orders


def test_the_run_log_is_written_by_an_ordinary_run_only(project, monkeypatch):
    log = project.path / 'test_timings' / 'runs.csv'
    project.runpytest_subprocess('-q', '-p', 'no:randomly').assert_outcomes(passed=3)
    rows = log.read_text().splitlines()
    assert rows[0].split(',')[:6] == ['timestamp', 'commit', 'scope', 'workers', 'load',
                                      'wall_seconds']
    assert len(rows) == 2
    project.runpytest_subprocess('-q', '-p', 'randomly', '--randomly-seed=3')
    monkeypatch.setenv('PYCIRCUIT_TEST_TIMINGS', '0')
    project.runpytest_subprocess('-q', '-p', 'no:randomly')
    assert len(log.read_text().splitlines()) == 2


def test_every_process_dumps_its_backend_state_when_asked(project, monkeypatch, tmp_path):
    out = tmp_path / 'state'
    monkeypatch.setenv('PYCIRCUIT_STATE_DUMP', str(out))
    project.makepyfile(test_hdl="""
        from pycircuit.circuit import elements_hdl as eh
        from pycircuit.circuit.elements import gnd
        def test_one_device():
            eh.MosLevel1Hdl('d', 'g', gnd, gnd)
    """)
    project.runpytest_subprocess('-q', '-p', 'no:randomly', '-n', '2').assert_outcomes(passed=4)
    names = sorted(p.name for p in out.iterdir())
    assert names == ['state-gw0.json', 'state-gw1.json', 'state-main.json']
    dumps = [json.loads((out / n).read_text()) for n in names]
    assert all(d['invariant_breaks'] == {} for d in dumps)
    assert any(any(k.endswith('MosLevel1Hdl_collapse11') for k in d['state']) for d in dumps)


def test_a_replay_file_is_read_in_start_order_and_cut_at_the_victim(tmp_path):
    p = tmp_path / '.pytest-replay-gw3.txt'
    recs = []
    for nid in ('m.py::a', 'm.py::b', 'm.py::c'):
        recs += [{'nodeid': nid, 'start': 1.0}, {'nodeid': nid, 'start': 1.0, 'finish': 2.0}]
    p.write_text('\n'.join(json.dumps(r) for r in recs) + '\n')
    order, lines = test_replay.read(str(p), 'm.py::b')
    assert order == ['m.py::a', 'm.py::b']
    assert len(lines) == 4
    assert test_replay.replay_file(str(tmp_path), 'gw3') == str(p)
    with pytest.raises(SystemExit):
        test_replay.read(str(p), 'm.py::zzz')


def test_the_polluter_is_found_from_a_recorded_order(pytester, monkeypatch):
    """End to end: a test that leaves a module global set, a victim that
    fails on it, a recorded order with innocents around them -- the script
    names the polluter."""
    pytester.makepyfile(shared="STATE = []\n")
    pytester.makepyfile(test_pol="""
        import shared
        def test_innocent_1():
            pass
        def test_polluter():
            shared.STATE.append(1)
        def test_innocent_2():
            pass
        def test_victim():
            assert shared.STATE == []
    """)
    order = ['test_pol.py::test_innocent_1', 'test_pol.py::test_polluter',
             'test_pol.py::test_innocent_2', 'test_pol.py::test_victim']
    rp = pytester.path / '.pytest-replay-gw0.txt'
    rp.write_text('\n'.join(json.dumps({'nodeid': n}) for n in order) + '\n')
    ## (an inner run must not inherit the outer run's switches: a gate's
    ## state dump, a replay's timing switch, an explicit order)
    for var in ('PYCIRCUIT_TEST_TIMINGS', 'PYCIRCUIT_TEST_ORDER', 'PYCIRCUIT_STATE_DUMP',
                'PYCIRCUIT_LEAKS', 'PYCIRCUIT_LEAKS_REPORT', 'PYCIRCUIT_LEAKS_RUN'):
        monkeypatch.delenv(var, raising=False)
    env = dict(os.environ, PYTHONPATH=str(pytester.path) + os.pathsep + os.environ.get('PYTHONPATH', ''))
    out = subprocess.run([sys.executable, os.path.join(SCRIPTS, 'find_polluter.py'),
                          'test_pol.py::test_victim', str(rp)],
                         cwd=pytester.path, env=env, capture_output=True, text=True, timeout=300,
                         check=False)
    assert 'test_pol.py::test_polluter' in out.stdout, out.stdout + out.stderr
