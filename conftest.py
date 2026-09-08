"""Session-wide pytest configuration.

Kept free of ``jax`` imports at module level: the environment variable below
only takes effect if it is set before JAX initialises its backend, and this
file is imported before any test module.
"""

import os

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
# ``test_timings/history.csv`` (timestamp, commit, wall seconds, test count,
# the slowest test and its seconds).  The per-run files are gitignored; the
# history is meant to be committed so the trend survives.  Under pytest-xdist
# the controller receives every worker's reports, so the record is complete
# and the hooks are skipped inside workers.  The one recorded fact that
# motivated this: a single injection-locking gate held the whole suite at
# 99 % for 14 minutes (30 CPU-minutes), found only with a root py-spy dump.
# ---------------------------------------------------------------------------
import json as _json
import subprocess as _subprocess
import time as _time
from datetime import datetime as _datetime, timezone as _timezone

_TIMINGS = {}


def _is_xdist_worker(config):
    return hasattr(config, 'workerinput')


def pytest_sessionstart(session):
    if _is_xdist_worker(session.config):
        return
    session.config._suite_t0 = _time.perf_counter()


def pytest_runtest_logreport(report):
    ## ``call`` only: setup/teardown are not what makes a test slow here.
    if report.when == 'call':
        _TIMINGS[report.nodeid] = float(report.duration)


def pytest_sessionfinish(session, exitstatus):
    config = session.config
    if _is_xdist_worker(config) or not _TIMINGS:
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
                    'tests': len(items),
                    'args': list(getattr(config, 'invocation_params').args),
                    'durations': items}, f, indent=1)
    hist = os.path.join(outdir, 'history.csv')
    new = not os.path.exists(hist)
    ## `scope`: the positional arguments, so a subset run is told apart from
    ## the full suite when reading the history (`pycircuit` = everything).
    scope = ' '.join(a for a in config.invocation_params.args if not a.startswith('-')) or '.'
    with open(hist, 'a') as f:
        if new:
            f.write('timestamp,commit,scope,wall_seconds,tests,slowest_test,slowest_seconds\n')
        f.write('%s,%s,%s,%.1f,%d,%s,%.1f\n' % (stamp, commit, scope.replace(',', ';'), wall,
                                               len(items), items[0][0], items[0][1]))


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
        ## the most recent record that covered the whole package
        if d.get('tests', 0) >= 1000:
            best = d
    return dict(best['durations']) if best else {}


def pytest_collection_modifyitems(session, config, items):
    if _is_xdist_worker(config):
        return
    rec = _last_full_record(os.path.dirname(os.path.abspath(__file__)))
    if not rec:
        return
    items.sort(key=lambda it: -rec.get(it.nodeid, -1.0))
    ## ⚠ MEASURED 2026-09-08: plain longest-first made NO difference (24:13
    ## against baselines of 24:33 and 24:40).  xdist's load scheduler opens
    ## by sending each worker a CONSECUTIVE slice of ~len/(4 n) items, so a
    ## longest-first list puts every long test into the FIRST worker's
    ## opening slice -- the opposite of balance -- and only the dispatch
    ## after that is dynamic.  Interleave: deal the sorted list round-robin
    ## into n bins and concatenate, so every opening slice carries its
    ## share of long tests in descending order.  Serial runs are untouched.
    try:
        n = int(getattr(config.option, 'numprocesses', 0) or 0)
    except (TypeError, ValueError):
        n = 0
    if n > 1:
        bins = [items[b::n] for b in range(n)]
        items[:] = [it for b in bins for it in b]
