"""Compare two recordings of `tran_recorder.py`, call by call:

    python benchmarks/tranrec/compare.py DIR_A DIR_B [TOL]

Each call is IDENTICAL, DIFF (the field, and its relative difference over
TOL, default 0: bit-for-bit), MISSING on one side, or CACHED: missing on
one side, but a test that shares a per-process cached helper made the
identical call there (see `CACHING_TESTS`).  The last line counts them.
"""
import sys, glob, pickle, os
import numpy as np
def load(d):
    out = {}
    for fn in glob.glob(os.path.join(d, 'rec_*.pkl')):
        for r in pickle.load(open(fn, 'rb')):
            out[(r['test'], r['k'])] = r
    return out
A, B = load(sys.argv[1]), load(sys.argv[2]); tol = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0


def sig(r):
    ## the call's CONTENT, whatever test made it: a helper that caches its
    ## transient per process (test_distortion's _MEASUREMENT_CACHE) runs it
    ## under whichever test asks first on that worker, so the same call can
    ## move between tests, or be made once where the other run made it twice
    parts = [repr(r.get('exc')), repr(sorted((r.get('stats') or {}).items()))]
    for f in ('t', 'x', 'event_times'):
        v = r.get(f)
        parts.append('-' if v is None else np.ascontiguousarray(v).tobytes().hex())
    return hash(tuple(parts))


## ⚠ CONTENT ALONE IS TOO WEAK: 201 of 706 calls are byte-identical to a call
## in ANOTHER test (shared fixtures), so a match anywhere would hide a
## genuinely dropped call ~30 % of the time.  A call legitimately moves only
## between tests that share a per-process cached helper -- today only
## test_distortion_vs_transient._measure (_MEASUREMENT_CACHE); keep this in
## step with the tests (grep for per-process caches around Transient).
CACHING_TESTS = (
    'test_distortion_vs_transient.py::',
    'test_distortion.py::test_higher_truncation_improves_accuracy',
    'test_distortion.py::test_the_improvement_is_large_enough_to_be_worth_it',
)
caching = lambda k: any(c in str(k[0]) for c in CACHING_TESTS)
sigA = {sig(r) for k, r in A.items() if caching(k)}
sigB = {sig(r) for k, r in B.items() if caching(k)}
keys = sorted(set(A) | set(B), key=lambda k: (str(k[0]), k[1]))
bad = 0; same = 0; cached = 0; worst = (0.0, None)
for k in keys:
    if k not in A or k not in B:
        side, other = ('A', sigA) if k not in A else ('B', sigB)
        rec = B[k] if k not in A else A[k]
        if caching(k) and sig(rec) in other:
            ## an identical call exists on the other side under another test
            print('CACHED', side, k); cached += 1; continue
        print('MISSING', side, k); bad += 1; continue
    a, b = A[k], B[k]
    msgs = []
    if a.get('exc') != b.get('exc'):
        msgs.append('exc %r vs %r' % (a.get('exc'), b.get('exc')))
    for f in ('t', 'x', 'event_times'):
        if (f in a) != (f in b):
            msgs.append('%s present in one only' % f); continue
        if f not in a: continue
        x, y = a[f], b[f]
        if x.shape != y.shape:
            msgs.append('%s shape %s vs %s' % (f, x.shape, y.shape)); continue
        if x.size == 0: continue
        sc = max(float(np.max(np.abs(x))), 1e-300)
        e = float(np.max(np.abs(x - y))) / sc
        if e > worst[0]: worst = (e, (k, f))
        if e > tol: msgs.append('%s rel %.2e' % (f, e))
    sa, sb = a.get('stats'), b.get('stats')
    if sa != sb:
        d = {kk: (sa.get(kk) if sa else None, sb.get(kk) if sb else None)
             for kk in set(sa or {}) | set(sb or {}) if (sa or {}).get(kk) != (sb or {}).get(kk)}
        msgs.append('stats %s' % d)
    for f in ('pcnr_solves', 'pcnr_fallbacks', 'pcnr_status'):
        if a.get(f) != b.get(f): msgs.append('%s %r vs %r' % (f, a.get(f), b.get(f)))
    if msgs:
        bad += 1; print('DIFF', k, '; '.join(msgs)[:600])
    else:
        same += 1
print('compared %d calls: %d identical, %d cached elsewhere, %d differ; worst array diff %.2e at %s' % (len(keys), same, cached, bad, worst[0], worst[1]))
