"""THE GATE'S TRANSIENT RECORDER (2026-09-23; in the repo since 2026-09-27):
a pytest plugin that records the OUTPUT of every `Transient.solve` call in
the suite, keyed by (test nodeid, call index within the test), so that two
gate runs can be compared call by call with `compare.py`.

A passing suite says every test met its own tolerance; this says every
transient in it produced the SAME numbers -- a change that moves a
transient by 1e-9 passes the tests and shows here.  It has gated every
commit since the events refactor (706 calls in the full suite).

Recorded per call: the time points (`sweep_values`), the solution `x`, the
step statistics without the timing entries, `event_times`, the PCNR
counters, or the exception if the solve raised (re-raised unchanged).

Usage (the plugin is opt-in, nothing in `pycircuit` imports it):

    PYTHONPATH=benchmarks/tranrec TRANREC_OUT=/tmp/rec_A \
        pytest pycircuit -q -p no:cacheprovider -p tran_recorder
    ... change the code, record again into /tmp/rec_B ...
    python benchmarks/tranrec/compare.py /tmp/rec_A /tmp/rec_B

Each xdist worker writes `<TRANREC_OUT>/rec_<worker>.pkl` at session end,
~60 MB for the full suite: keep recordings OUT of the repo.
"""
import os, pickle
import numpy as np

_REC = []
_CUR = {'id': None, 'k': 0}
_TIMING = ('solve_seconds', 'total_seconds')


def pytest_configure(config):
    from pycircuit.circuit import transient as T
    orig = T.Transient.solve
    if getattr(orig, '_tranrec', False):
        return

    def solve(self, *a, **kw):
        k = _CUR['k']; _CUR['k'] += 1
        rec = {'test': _CUR['id'], 'k': k, 'cls': type(self).__name__}
        try:
            res = orig(self, *a, **kw)
        except BaseException as e:           # record, then re-raise unchanged
            rec['exc'] = '%s: %s' % (type(e).__name__, str(e)[:400])
            _REC.append(rec)
            raise
        try:
            rec['t'] = np.array(np.asarray(res.sweep_values), dtype=float)
            rec['x'] = np.array(np.asarray(res.x), dtype=float)
        except Exception as e:                # noqa: BLE001
            rec['t_err'] = repr(e)[:200]
        st = getattr(self, 'statistics', None)
        if st is not None and hasattr(st, 'as_dict'):
            rec['stats'] = {kk: v for kk, v in st.as_dict().items() if kk not in _TIMING}
        et = getattr(self, 'event_times', None)
        if et is not None:
            rec['event_times'] = np.array(np.asarray(et, dtype=float))
        for a_ in ('pcnr_solves', 'pcnr_fallbacks', 'pcnr_status'):
            if hasattr(self, a_):
                rec[a_] = getattr(self, a_)
        _REC.append(rec)
        return res
    solve._tranrec = True
    T.Transient.solve = solve


def pytest_runtest_setup(item):
    _CUR['id'] = item.nodeid
    _CUR['k'] = 0


def pytest_sessionfinish(session, exitstatus):
    out = os.environ.get('TRANREC_OUT')
    if not out:
        return
    os.makedirs(out, exist_ok=True)
    wid = os.environ.get('PYTEST_XDIST_WORKER', 'main')
    with open(os.path.join(out, 'rec_%s.pkl' % wid), 'wb') as f:
        pickle.dump(_REC, f)
