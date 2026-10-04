"""Session-wide pytest configuration.

Kept free of ``jax`` imports at module level: the environment variable below
only takes effect if it is set before JAX initialises its backend, and this
file is imported before any test module.
"""

import os

# BLAS THREADS: one per worker.  The suite ran at 24-25 min wall against a
# 4200 s serial sum -- an ideal of ~7 min at -n 10 -- and reordering the
# collection (longest first, interleaved) did not move it: the gap was
# CONTENTION, not scheduling.  Measured 2026-09-08: the dense-linear-algebra
# tests ran 6x slower inside the suite than alone (B16: 25 s alone, 180 s in
# the suite) while a Newton-ladder test ran only 1.4x slower.  scipy-openblas
# opens a thread per core (24 here) in EVERY worker, so 10 workers fight over
# 240 threads.  Set before numpy loads; `pytest_configure` below also limits
# any pool already open (threadpoolctl) for the controller process.
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

# JAX preallocates ~75% of the GPU's memory in each process the first time a
# device is used.  Under pytest-xdist every worker is a separate process, so on
# a GPU with modest VRAM all but the first worker fails to obtain memory --
# and the resulting failures surface as ordinary assertion errors rather than
# an out-of-memory message, which makes them easy to misread as numerical bugs
# in the solver.  Allocate on demand instead.
#
# Measured on an RTX A1000 (4 GB) with ``-n 8``: 15 failures with
# preallocation, 0 without.  Harmless on CPU-only installations.
#
# ``setdefault`` so an explicit setting in the environment still wins.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

# Some example tests call ``pylab.show()``.  On a developer machine with DISPLAY
# set, matplotlib picks an interactive backend (tkagg) and show() opens a window
# and blocks until something closes it -- so the suite appears to hang, for
# minutes at a time, and nondeterministically depending on the window manager.
# test_my_elements.py::test_transient_plot took 148 s this way against 0.96 s
# headless.  Tests should never open a window; force the non-interactive
# backend.  Again setdefault, so running with MPLBACKEND=tkagg to eyeball a
# plot still works.
os.environ.setdefault("MPLBACKEND", "Agg")


# ---------------------------------------------------------------------------
# SUITE TIMING RECORD (Andreas, 2026-09-08: "start recording suite execution
# time; some tests are extending the suite's duration a lot").
#
# Every run writes ``test_timings/<utc-timestamp>_<commit>.json`` -- the call
# duration of every test, sorted slowest first -- and appends one line to
# ``test_timings/runs.csv`` (timestamp, commit, scope, workers, load, wall
# seconds, test count, the slowest test and its seconds).  Both are local and
# gitignored since 2026-10-03 (`history.csv`, the tracked record before, lost
# rows whenever it was restored after a targeted run).  Under pytest-xdist
# the controller receives every worker's reports, so the record is complete
# and the hooks are skipped inside workers.  The one recorded fact that
# motivated this: a single injection-locking gate held the whole suite at
# 99 % for 14 minutes (30 CPU-minutes), found only with a root py-spy dump.
# ---------------------------------------------------------------------------
import json as _json
import subprocess as _subprocess
import time as _time
from datetime import datetime as _datetime, timezone as _timezone

## THE LEAK DETECTOR (2026-10-03): a snapshot of all process-global state
## around every test, every module and the collection; see its module note.
pytest_plugins = ['pycircuit._testing.leaks']

_TIMINGS = {}


def _is_xdist_worker(config):
    return hasattr(config, 'workerinput')


def pytest_configure(config):
    """Limit BLAS pools already open in this process to one thread (see the
    environment pin at the top: that covers pools opened AFTER it)."""
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:
        pass
    _bench_lock_shared(config)


# ---------------------------------------------------------------------------
# THE BENCHMARK LOCK, SHARED (2026-10-03; `benchmarks/_bench.py`): a test
# session of this repository holds `~/.cache/pycircuit/bench.lock` shared,
# benchmarks take it exclusively -- so a benchmark never measures beside this
# suite, and a suite started during a benchmark waits for it (saying so).
# The controller holds it for the session; workers do not need to.
# ---------------------------------------------------------------------------
def _bench_lock_shared(config):
    if _is_xdist_worker(config) or os.environ.get('PYCIRCUIT_BENCH_LOCK') == '0':
        return
    try:
        import fcntl
    except ImportError:
        return
    path = os.path.join(os.path.expanduser('~'), '.cache', 'pycircuit', 'bench.lock')
    try:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        fd = os.open(path, os.O_RDWR | os.O_CREAT, 0o644)
    except OSError:
        return
    try:
        fcntl.flock(fd, fcntl.LOCK_SH | fcntl.LOCK_NB)
    except OSError:
        print('waiting for a running benchmark (the benchmark lock, benchmarks/_bench.py)',
              flush=True)
        fcntl.flock(fd, fcntl.LOCK_SH)
    config._pycircuit_bench_lock = fd


def pytest_unconfigure(config):
    fd = getattr(config, '_pycircuit_bench_lock', None)
    if fd is not None:
        try:
            os.close(fd)
        except OSError:
            pass


# ---------------------------------------------------------------------------
# THE FAST TIER (2026-10-04, testing for development): `--tier fast` deselects
# every test whose last COMPLETE record took >= `PYCIRCUIT_SLOW_SECONDS`
# (default 5 s) -- in the 2026-10-03 record 263 tests, 82 % of the suite's
# test-seconds; the rest runs in ~3 min at -n 8.  For iterating; the full
# suite stays the gate.  A test the record does not know runs (it is new).
# The `slow` marker is a different thing (transient/ODE tests, pytest.ini)
# and is left alone.  Deselection happens where collection does: in every
# xdist worker, from the same record, so the workers agree.
# ---------------------------------------------------------------------------
def pytest_addoption(parser):
    parser.addoption(
        '--tier', choices=('all', 'fast'), default='all',
        help='fast: deselect the tests whose last complete timing record took '
             '>= PYCIRCUIT_SLOW_SECONDS (default 5) seconds -- for iterating; '
             'the full suite before a commit')


def _slow_seconds():
    return float(os.environ.get('PYCIRCUIT_SLOW_SECONDS', '5'))


def pytest_report_header(config):
    if config.getoption('tier', 'all') != 'fast':
        return None
    rec = _last_full_record(os.path.dirname(os.path.abspath(__file__)))
    if not rec:
        return 'tier fast: no complete timing record yet -- every test runs'
    limit = _slow_seconds()
    n = sum(1 for d in rec.values() if d >= limit)
    return (f'tier fast: the {n} tests that took >= {limit:g} s in the last complete '
            f'record ({len(rec)} tests) are deselected')


def pytest_sessionstart(session):
    if _is_xdist_worker(session.config):
        return
    session.config._suite_t0 = _time.perf_counter()


def pytest_runtest_logreport(report):
    ## ``call`` only: setup/teardown are not what makes a test slow here.
    if report.when == 'call':
        _TIMINGS[report.nodeid] = float(report.duration)


def _randomly_active(config):
    """pytest-randomly loaded (the ini blocks it; `-p randomly` unblocks)
    and reordering."""
    return (config.pluginmanager.hasplugin('randomly')
            and bool(config.getoption('randomly_reorganize', False)))


def _replaying(config):
    return bool(config.getoption('replay_files', None))


def _dump_state(config):
    """THE PER-WORKER BACKEND STATE (2026-10-03): with `PYCIRCUIT_STATE_DUMP`
    set, every process -- each xdist worker and the controller -- writes the
    hdl backend state it ended with (`pycircuit._testing.state.hdl_state`),
    so two workers can be diffed after a one-worker anomaly (gates G97/G98:
    a class bound to C with no kernels, in one worker of eight)."""
    out = os.environ.get('PYCIRCUIT_STATE_DUMP')
    if not out:
        return
    try:
        from pycircuit._testing import state as _state
        os.makedirs(out, exist_ok=True)
        who = getattr(config, 'workerinput', {}).get('workerid', 'main')
        with open(os.path.join(out, f'state-{who}.json'), 'w') as f:
            _json.dump({'worker': who, 'state': _state.hdl_state(),
                        'invariant_breaks': _state.invariant_breaks()}, f, indent=1,
                       sort_keys=True)
    except Exception as e:  # noqa: BLE001 -- a diagnostic must never fail a run
        print(f'state dump failed: {e!r}')


def pytest_sessionfinish(session, exitstatus):
    config = session.config
    _dump_state(config)
    if _is_xdist_worker(config) or not _TIMINGS:
        return
    ## (not the runs that are not the suite's timing: a shuffled order, a
    ## replay, a bisection's many sub-runs, or anyone who says so)
    if (os.environ.get('PYCIRCUIT_TEST_TIMINGS', '1') == '0'
            or _randomly_active(config) or _replaying(config)):
        return
    root = os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.join(root, 'test_timings')
    os.makedirs(outdir, exist_ok=True)
    try:
        commit = _subprocess.check_output(
            ['git', 'rev-parse', '--short', 'HEAD'], cwd=root,
            stderr=_subprocess.DEVNULL).decode().strip()
    except Exception:
        commit = 'unknown'
    stamp = _datetime.now(_timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    wall = _time.perf_counter() - getattr(config, '_suite_t0', _time.perf_counter())
    items = sorted(_TIMINGS.items(), key=lambda kv: -kv[1])
    with open(os.path.join(outdir, '%s_%s.json' % (stamp, commit)), 'w') as f:
        _json.dump({'timestamp': stamp, 'commit': commit, 'wall_seconds': wall,
                    'tests': len(items), 'complete': _complete_run(session),
                    'args': list(getattr(config, 'invocation_params').args),
                    'durations': items}, f, indent=1)
    ## THE RUN LOG (2026-10-03): `runs.csv`, local and untracked.  Until then
    ## `history.csv` was tracked in git, and restoring it after targeted runs
    ## discarded the full runs' rows too (29 lost in one day); it is kept as
    ## the old record.  Each row now says how the suite ran (workers, the
    ## load average) so two rows are comparable, and appends under a lock.
    hist = os.path.join(outdir, 'runs.csv')
    ## `scope`: the positional arguments, so a subset run is told apart from
    ## the full suite when reading the history (`pycircuit` = everything);
    ## the option values that follow a flag are not part of it.
    args = list(config.invocation_params.args)
    scope = ' '.join(a for j, a in enumerate(args) if not a.startswith('-')
                     and not (j and args[j - 1] in ('-n', '-p', '-k', '-m', '-c', '-o'))) or '.'
    workers = config.getoption('numprocesses', None)
    try:
        load = f'{os.getloadavg()[0]:.1f}'
    except OSError:
        load = ''
    line = (f"{stamp},{commit},{scope.replace(',', ';')},{workers},{load},{wall:.1f},"
            f"{len(items)},{items[0][0]},{items[0][1]:.1f}\n")
    with open(hist, 'a') as f:
        try:
            import fcntl
            fcntl.flock(f, fcntl.LOCK_EX)
        except (ImportError, OSError):
            pass
        if f.tell() == 0:
            f.write('timestamp,commit,scope,workers,load,wall_seconds,tests,'
                    'slowest_test,slowest_seconds\n')
        f.write(line)


def _complete_run(session):
    """THE RECORD GUARD (2026-10-04): a run is the suite's timing only when it
    ran everything -- the package (or the root) and no selection: no `-m`,
    `-k`, `--deselect`, last-failed / failed-first / stepwise, `--tier
    fast`, and not stopped early (`-x`, `--maxfail`).  Until then any run of
    1000 tests or more counted, and a `-m "not slow"` run (~3830 tests)
    became THE record the sort and the fast tier read, with the deselected
    tests missing from it.  Decided from the options, in the controller:
    under xdist only the workers collect and deselect."""
    config = session.config
    if config.getoption('tier', 'all') != 'all':
        return False
    if config.getoption('markexpr', '') or config.getoption('keyword', ''):
        return False
    if config.getoption('deselect', None):
        return False
    for opt in ('lf', 'failedfirst', 'stepwise', 'stepwise_skip'):
        if config.getoption(opt, False):
            return False
    if config.getoption('maxfail', 0) or session.shouldstop or session.shouldfail:
        return False
    root = os.path.dirname(os.path.abspath(__file__))
    ## pytest's own parse of the positional paths (an option's value is
    ## never one), the invocation directory when none were given
    for a in config.args:
        rel = os.path.relpath(os.path.abspath(os.path.join(
            str(config.invocation_params.dir), a.split('::', 1)[0])), root)
        if rel not in ('.', 'pycircuit'):
            return False
    return True


# ---------------------------------------------------------------------------
# LONGEST-FIRST COLLECTION ORDER, from the timing record (2026-09-08).
#
# The first full record said where the suite's time goes: the serial sum of
# all call durations was 4201 s, so ten workers could finish in ~7 min, yet
# the run took 24.5 min -- fourteen tests over a minute, the longest 193 s,
# landed wherever collection order put them and the tail waited on them.
# Sorting collected items by their last recorded duration, longest first,
# is the classic longest-processing-time heuristic: xdist's load scheduler
# hands out items in collection order, so every worker gets a long test
# early and the short ones fill the gaps.  Tests with no record run last
# (they are new, and new tests are usually short).  No test changes; a
# missing or unreadable record leaves the order untouched.
# ---------------------------------------------------------------------------
def _last_full_record(root):
    outdir = os.path.join(root, 'test_timings')
    if not os.path.isdir(outdir):
        return {}
    best = None
    for name in sorted(os.listdir(outdir)):
        if not name.endswith('.json'):
            continue
        try:
            with open(os.path.join(outdir, name)) as f:
                d = _json.load(f)
        except Exception:
            continue
        ## the most recent record of a COMPLETE run (`_complete_run`); a
        ## record from before 2026-10-04 has no flag and counts when it
        ## covered 1000 tests or more, as it did then
        if d.get('complete', d.get('tests', 0) >= 1000) is True:
            best = d
    return dict(best['durations']) if best else {}


def _recorded_duration(rec):
    """`nodeid -> seconds or None` from a record.  ⚠ BY NAME WHEN THE NODEID
    IS NEW (2026-09-27): a test MOVED to another file keeps its name, and
    without this every moved test ran LAST as "unknown" -- measured by
    simulation for the split of test_analysis_shooting.py, 1304 s -> ~1370 s
    until the next full run.  Only names unique in the record are used."""
    by_name, seen = {}, set()
    for nid, dur in rec.items():
        nm = nid.split('::', 1)[-1]
        if nm in seen:
            by_name.pop(nm, None)
        else:
            by_name[nm] = dur
            seen.add(nm)

    def dur(nodeid):
        d = rec.get(nodeid)
        return by_name.get(nodeid.split('::', 1)[-1]) if d is None else d
    return dur


def _select_tier(config, items, rec, dur):
    if config.getoption('tier', 'all') != 'fast' or not rec:
        return
    limit = _slow_seconds()
    keep, slow = [], []
    for it in items:
        d = dur(it.nodeid)
        ## (`@pytest.mark.fast_tier`: the cheapest test of a feature no
        ## faster test reaches, kept -- chosen by coverage, see pytest.ini)
        if d is not None and d >= limit and it.get_closest_marker('fast_tier') is None:
            slow.append(it)
        else:
            keep.append(it)
    if slow:
        config.hook.pytest_deselected(items=slow)
        items[:] = keep


import pytest as _pytest


@_pytest.hookimpl(tryfirst=True)
def pytest_collection_modifyitems(session, config, items):
    """Longest first, from the last full-suite record, IN EVERY PROCESS.

    THE ORDER MODES (2026-10-03).  `tryfirst`, so an order another plugin
    imposes afterwards wins: pytest-replay's (`--replay`) and
    detect-test-pollution's.  And it steps aside entirely for a shuffled run
    (pytest-randomly shuffles in a tryfirst WRAPPER before every plain hook,
    so sorting after it would undo the shuffle), for a replay, and under
    `PYCIRCUIT_TEST_ORDER=collected` (the collection's own order).

    ⚠⚠ MEASURED WRONG TWICE, 2026-09-08: the first version returned early
    on xdist workers -- and under xdist ONLY the workers collect, so the
    reorder never ran; the "plain longest-first made no difference" (24:13)
    and "interleaved made no difference" (25:28) readings were both of the
    UNMODIFIED order.  Every worker sorts the same list from the same file,
    so they agree (xdist requires identical collections).

    Pairs with `--maxschedchunk=1` in pytest.ini: xdist's load scheduler
    otherwise hands each worker a consecutive quarter-share to open and
    then half-shares of the remainder as it drains, none of which it can
    take back, so a long test late in a chunk leaves the other workers
    idle -- the 10-minute tail at 96 % that the record shows.  One test at
    a time from a longest-first list is LPT list scheduling.
    """
    rec = _last_full_record(os.path.dirname(os.path.abspath(__file__)))
    dur = _recorded_duration(rec)
    ## the fast tier selects in EVERY order mode (a shuffled or replayed
    ## fast run is still the fast tier); only the sort steps aside
    _select_tier(config, items, rec, dur)
    if (_randomly_active(config) or _replaying(config)
            or os.environ.get('PYCIRCUIT_TEST_ORDER') == 'collected'):
        return
    if not rec:
        return

    def key(it):
        d = dur(it.nodeid)
        return -(-1.0 if d is None else d)
    items.sort(key=key)
