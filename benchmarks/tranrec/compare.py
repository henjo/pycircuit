"""Compare two recordings of `tran_recorder.py`, call by call:

    python benchmarks/tranrec/compare.py DIR_A DIR_B [TOL] [--by-name]
                                         [--family F] [--quiet]

Each call is IDENTICAL, DIFF (what differs, with the relative difference
of each array over TOL, default 0: BIT-FOR-BIT), MISSING on one side, or
CACHED: missing on one side, but a test that shares a per-process cached
helper made the identical call there (see `CACHING_TESTS`).  Per test the
warnings are compared by (category, message); a warning whose text is the
same but whose reported location moved is counted, not failed.  The last
lines count each family.

BIT MODE (TOL 0) compares the bytes: a NaN is equal only to the same NaN,
-0.0 differs from 0.0.  TOLERANCE MODE requires the non-finite entries to
sit in the same places with the same values and takes the relative
difference over the finite ones.  (Before 2026-09-27 a NaN difference read
as IDENTICAL here: `nan > tol` is False.)

`--only-common` compares only the tests that RAN on both sides (and the
families recorded on both) -- a subset run against a full recording, while
iterating.  A test that ran on both sides with calls on one only is
MISSING (before 2026-10-04 "common" meant "made a recorded call on both
sides", and such a test dropped out).  A recording older than the run list
falls back to the tests with calls.

FAMILIES are read from each recording's header: a family recorded on one
side only is NAMED and not compared (a gate recorded with fewer families,
or a recording older than a new family), instead of every call reading
MISSING.

`--passed-only` compares only the tests that PASSED on both sides, and
names how many it left out (the switches-off check: a test asserting that
a fast path served fails by design when it is off).

`--ignore-added-state` ignores a `state` leaf present in B only: across the
recorder's 2026-09-27 fix (it compared ids, and a reassigned attribute whose
new value was born at the freed address read as unchanged), a recording made
before it misses some reassignments at random.

`--by-name` keys on the test NAME (`file.py::name[params]` -> `name[params]`)
instead of the nodeid, for a split test file; a name that occurs in two
files of one recording keeps its nodeid (counted on the last line).
"""
import sys, glob, pickle, os, argparse
from collections import Counter
import numpy as np


## ⚠ CONTENT ALONE IS TOO WEAK: 201 of 706 calls are byte-identical to a call
## in ANOTHER test (shared fixtures), so a match anywhere would hide a
## genuinely dropped call ~30 % of the time.  A call legitimately moves only
## between tests that share a per-process cached helper -- today only
## test_distortion_vs_transient._measure (_MEASUREMENT_CACHE); keep this in
## step with the tests (grep for per-process caches around Transient, PSS and
## PAC -- none caches a PSS or PAC result across tests as of 2026-09-27).
CACHING_TESTS = (
    'test_distortion_vs_transient.py::',
    'test_distortion.py::test_higher_truncation_improves_accuracy',
    'test_distortion.py::test_the_improvement_is_large_enough_to_be_worth_it',
)


def _key_name(nodeid):
    return nodeid.split('::', 1)[1] if '::' in nodeid else nodeid


def load(d, by_name, meta=None):
    """`(calls, warnings, collisions)`; `meta`, a dict when given, receives
    the families header (`families`), the tests that ran (`ran`) and their
    outcomes (`outcomes`) -- each None when a recording predates it."""
    calls, warns = [], []
    fams, ran, outcomes = set(), set(), {}
    old_fams = old_ran = False
    for fn in sorted(glob.glob(os.path.join(d, 'rec_*.pkl'))):
        with open(fn, 'rb') as f:
            data = pickle.load(f)
        if isinstance(data, list):               # the transient-only format
            calls.extend(data)
            fams.add('transient')
            old_ran = True
        else:
            calls.extend(data['calls'])
            warns.extend(data['warnings'])
            if 'families' in data:
                fams.update(data['families'])
            else:
                old_fams = True
            if 'ran' in data:
                ran.update(data['ran'])
                outcomes.update(data.get('outcomes', {}))
            else:
                old_ran = True
    names = {}
    if by_name:
        per = {}
        for r in calls + warns:
            per.setdefault(_key_name(str(r['test'])), set()).add(r['test'])
        names = {nm: ids for nm, ids in per.items()}
    collisions = {nm for nm, ids in names.items() if len(ids) > 1}

    def tkey(t):
        ## (a cache-sharing test keeps its nodeid: `CACHING_TESTS` names files)
        if not by_name or any(c in str(t) for c in CACHING_TESTS):
            return t
        nm = _key_name(str(t))
        return t if nm in collisions else nm
    out = {}
    for r in calls:
        out[(tkey(r['test']), r.get('fam', 'transient'), r['k'])] = r
    w = {}
    for r in warns:
        w.setdefault(tkey(r['test']), []).append(r)
    if meta is not None:
        meta['families'] = None if old_fams else fams
        meta['ran'] = None if old_ran else {tkey(t) for t in ran}
        meta['outcomes'] = None if old_ran else {tkey(t): v for t, v in outcomes.items()}
    return out, w, len(collisions)


def rel(x, y):
    """The relative difference of two finite-masked arrays, inf when the
    non-finite entries differ in place or value."""
    fx, fy = np.isfinite(x), np.isfinite(y)
    if not np.array_equal(fx, fy):
        return float('inf')
    if not np.array_equal(x[~fx], y[~fy], equal_nan=True):
        return float('inf')
    if not fx.any():
        return 0.0
    xf, yf = x[fx], y[fy]
    sc = max(float(np.max(np.abs(xf))), 1e-300)
    return float(np.max(np.abs(xf - yf))) / sc


def leaf_diff(x, y, tol):
    """None when equal (bit-for-bit at tol 0, within tol otherwise), else
    (message, relative difference)."""
    if isinstance(x, tuple) and x and x[0] == '__big__':
        if not (isinstance(y, tuple) and y and y[0] == '__big__'):
            return 'big vs not', float('inf')
        if x[1] == y[1]:
            return None
        if x[2:4] != y[2:4]:
            return 'big %s vs %s' % (x[2:4], y[2:4]), float('inf')
        e = rel(x[4], y[4])
        return (None if tol and e <= tol else ('big digest, sample rel %.2e' % e, e))
    xa, ya = isinstance(x, np.ndarray), isinstance(y, np.ndarray)
    if xa != ya:
        return '%s vs %s' % (type(x).__name__, type(y).__name__), float('inf')
    if not xa:
        if isinstance(x, (float, complex)) and isinstance(y, (float, complex)):
            if repr(x) == repr(y):
                return None
            e = rel(np.atleast_1d(np.asarray(x)), np.atleast_1d(np.asarray(y)))
            return None if tol and e <= tol else ('%r vs %r' % (x, y), e)
        if type(x) is not type(y) or x != y:
            return ('%r vs %r' % (x, y))[:200], float('inf')
        return None
    if x.dtype != y.dtype:
        return 'dtype %s vs %s' % (x.dtype, y.dtype), float('inf')
    if x.shape != y.shape:
        return 'shape %s vs %s' % (x.shape, y.shape), float('inf')
    if x.dtype.kind not in 'fc':
        return None if np.array_equal(x, y) else ('values differ', float('inf'))
    if x.size == 0 or x.tobytes() == y.tobytes():
        return None
    e = rel(x, y)
    if tol and e <= tol:
        return None
    return 'rel %.2e' % e, e


IGNORE_ADDED_STATE = False


def leaves_diff(A, B, what, tol, msgs, worst, k):
    for p in sorted(set(A) | set(B)):
        if IGNORE_ADDED_STATE and what == 'state' and (
                p not in A or (p == '' and A[p] == '<empty dict>')):
            continue
        if p not in A or p not in B:
            msgs.append('%s%s present in one only' % (what, p))
            continue
        d = leaf_diff(A[p], B[p], tol)
        if d is None:
            continue
        m, e = d
        if e > worst[0] and np.isfinite(e):
            worst[:] = [e, (k, what + p)]
        msgs.append('%s%s %s' % (what, p, m))


def compare_call(a, b, tol, worst, k):
    msgs = []
    if a.get('exc') != b.get('exc'):
        msgs.append('exc %r vs %r' % (a.get('exc'), b.get('exc')))
    ## (a transient-only recording, before 2026-09-27, carries no name)
    if a.get('name') and b.get('name') and a['name'] != b['name']:
        msgs.append('call %r vs %r' % (a.get('name'), b.get('name')))
    for f in ('t', 'x', 'event_times'):
        if (f in a) != (f in b):
            msgs.append('%s present in one only' % f)
        elif f in a:
            leaves_diff({'': a[f]}, {'': b[f]}, f, tol, msgs, worst, k)
    for f in ('inp', 'out', 'state'):
        if (f in a) != (f in b):
            msgs.append('%s present in one only' % f)
        elif f in a:
            leaves_diff(a[f], b[f], f, tol, msgs, worst, k)
    sa, sb = a.get('stats'), b.get('stats')
    if sa != sb:
        d = {kk: (sa.get(kk) if sa else None, sb.get(kk) if sb else None)
             for kk in set(sa or {}) | set(sb or {}) if (sa or {}).get(kk) != (sb or {}).get(kk)}
        msgs.append('stats %s' % d)
    for f in ('pcnr_solves', 'pcnr_fallbacks', 'pcnr_status'):
        if a.get(f) != b.get(f):
            msgs.append('%s %r vs %r' % (f, a.get(f), b.get(f)))
    return msgs


def sig(r):
    ## the call's CONTENT, whatever test made it: a helper that caches its
    ## transient per process (test_distortion's _MEASUREMENT_CACHE) runs it
    ## under whichever test asks first on that worker, so the same call can
    ## move between tests, or be made once where the other run made it twice
    ## (the arrays and the exception, not the statistics: a counter added to
    ## the statistics must not hide a moved call as MISSING)
    parts = [repr(r.get('exc'))]
    for f in ('t', 'x', 'event_times'):
        v = r.get(f)
        parts.append('-' if v is None else np.ascontiguousarray(v).tobytes().hex())
    return hash(tuple(parts))




def main(argv):
    ap = argparse.ArgumentParser()
    ap.add_argument('a'); ap.add_argument('b')
    ap.add_argument('tol', nargs='?', type=float, default=0.0)
    ap.add_argument('--by-name', action='store_true')
    ap.add_argument('--family', default=None)
    ap.add_argument('--quiet', action='store_true')
    ap.add_argument('--only-common', action='store_true')
    ap.add_argument('--passed-only', action='store_true')
    ap.add_argument('--ignore-added-state', action='store_true')
    o = ap.parse_args(argv)
    global IGNORE_ADDED_STATE
    IGNORE_ADDED_STATE = o.ignore_added_state
    MA, MB = {}, {}
    A, WA, ca = load(o.a, o.by_name, MA)
    B, WB, cb = load(o.b, o.by_name, MB)
    if o.family:
        A = {k: v for k, v in A.items() if k[1] == o.family}
        B = {k: v for k, v in B.items() if k[1] == o.family}
    say = (lambda *s: None) if o.quiet else print
    caching = lambda k: k[1] == 'transient' and any(c in str(k[0]) for c in CACHING_TESTS)
    ## (the cached-elsewhere matches from the WHOLE recording, before any
    ## filter: a cached call the full run made in a test the subset did not
    ## run is still the same call)
    sigA = {sig(r) for k, r in A.items() if caching(k)}
    sigB = {sig(r) for k, r in B.items() if caching(k)}
    ## a family one side did not record is not compared, and said so
    fa, fb = MA['families'], MB['families']
    if fa is not None and fb is not None and fa != fb:
        for fam in sorted(fa ^ fb):
            print(f"NOTICE family {fam!r} recorded in {'A' if fam in fa else 'B'} only: "
                  'not compared')
        A = {k: v for k, v in A.items() if k[1] in fb}
        B = {k: v for k, v in B.items() if k[1] in fa}
    keep = None
    if o.only_common:
        ra = MA['ran'] if MA['ran'] is not None else {k[0] for k in A}
        rb = MB['ran'] if MB['ran'] is not None else {k[0] for k in B}
        keep = ra & rb
    if o.passed_only:
        if MA['outcomes'] is None or MB['outcomes'] is None:
            sys.exit('--passed-only needs recordings with outcomes (2026-10-04 or later)')
        passed = {t for t, v in MA['outcomes'].items()
                  if v == 'passed' and MB['outcomes'].get(t) == 'passed'}
        ran = (set(MA['outcomes']) | set(MB['outcomes'])) if keep is None else keep
        left = sorted((ran - passed), key=str)
        names = (': ' + ', '.join(map(str, left[:8])) + (' ...' if len(left) > 8 else '')
                 if left else '')
        print(f'NOTICE {len(left)} tests not compared (not passed on both sides){names}')
        keep = passed if keep is None else keep & passed
    if keep is not None:
        A = {k: v for k, v in A.items() if k[0] in keep}
        B = {k: v for k, v in B.items() if k[0] in keep}
        WA = {t: v for t, v in WA.items() if t in keep}
        WB = {t: v for t, v in WB.items() if t in keep}
    keys = sorted(set(A) | set(B), key=lambda k: (str(k[0]), k[1], k[2]))
    count = {}
    worst = {}
    for k in keys:
        c = count.setdefault(k[1], Counter())
        w = worst.setdefault(k[1], [0.0, None])
        if k not in A or k not in B:
            side, other = ('A', sigA) if k not in A else ('B', sigB)
            rec = B[k] if k not in A else A[k]
            if caching(k) and sig(rec) in other:
                say('CACHED', side, k); c['cached'] += 1; continue
            print('MISSING', side, k, rec.get('name')); c['bad'] += 1; continue
        msgs = compare_call(A[k], B[k], o.tol, w, k)
        if msgs:
            c['bad'] += 1
            print('DIFF', k, A[k].get('name'), '; '.join(msgs)[:800])
        else:
            c['same'] += 1
    ## ⚠ WARNINGS PER TEST ARE NOT ALL DETERMINISTIC under xdist (measured,
    ## two runs of one commit: 146 tests differed).  pytest's OWN warnings
    ## (a class-scoped fixture's deprecation, once per class per worker)
    ## land on whichever test of the class a worker runs first -- they are
    ## about the suite, not the code, and are dropped; and a warning from a
    ## cache-sharing helper lands on whichever of those tests runs it first,
    ## so theirs are pooled into one bucket, as their calls are matched.
    def _norm(W):
        out = {}
        for t, rs in W.items():
            ## a warning whose stacklevel lands on the recorder's own wrapper
            ## is located there by line; the line moves when the recorder is
            ## edited, which is not the code under test moving
            for r in rs:
                if r['where'].startswith('tran_recorder.py:'):
                    r['where'] = 'tran_recorder.py'
            key = '<the tests sharing a cache>' if any(c in str(t) for c in CACHING_TESTS) else t
            out.setdefault(key, []).extend(r for r in rs if not r['cat'].startswith('Pytest'))
        ## the pooled bucket as a SET: the helper runs once per worker that
        ## asks, and how many workers ask varies (4 emissions vs 5, measured)
        pool = out.get('<the tests sharing a cache>')
        if pool:
            uniq = {(r['cat'], r['msg'], r['where']): r for r in pool}
            out['<the tests sharing a cache>'] = list(uniq.values())
        return {t: rs for t, rs in out.items() if rs}
    WA, WB = _norm(WA), _norm(WB)
    wbad = wmoved = 0
    if WA or WB:
        for t in sorted(set(WA) | set(WB), key=str):
            ra, rb = WA.get(t, []), WB.get(t, [])
            ma = Counter((r['cat'], r['msg']) for r in ra)
            mb = Counter((r['cat'], r['msg']) for r in rb)
            if ma != mb:
                wbad += 1
                print('WARN', t, 'only A: %s' % list((ma - mb).elements())[:3],
                      'only B: %s' % list((mb - ma).elements())[:3])
            elif (Counter((r['cat'], r['msg'], r['where']) for r in ra)
                  != Counter((r['cat'], r['msg'], r['where']) for r in rb)):
                wmoved += 1
    for fam in sorted(count):
        c = count[fam]
        print('%-9s compared %d calls: %d identical, %d cached elsewhere, %d differ; '
              'worst array diff %.2e at %s' % (
                  fam, sum(c.values()), c['same'], c['cached'], c['bad'],
                  worst[fam][0], worst[fam][1]))
    print('warnings: %d tests differ, %d tests with a moved location only; '
          'name collisions kept as nodeids: %d / %d' % (wbad, wmoved, ca, cb))
    return 1 if (wbad or any(c['bad'] for c in count.values())) else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
