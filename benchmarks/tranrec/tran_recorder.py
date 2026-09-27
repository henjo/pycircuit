"""THE GATE'S OUTPUT RECORDER (2026-09-23; in the repo since 2026-09-27;
families since 2026-09-27): a pytest plugin that records the OUTPUT of the
analyses the suite runs, keyed by (test, family, call index within the test
and family), so that two gate runs can be compared call by call with
`compare.py`.

A passing suite says every test met its own tolerance; this says every
analysis in it produced the SAME numbers -- a change that moves a result by
1e-9 passes the tests and shows here.  The transient family has gated every
commit since the events refactor (706 calls in the full suite).

FAMILIES (`TRANREC_FAMILIES`, comma separated or `all`; default `transient`,
so a recording made with the old command compares with a new one):

- `transient`: every `Transient.solve` call (nested ones too, as always).
  Recorded: the time points (`sweep_values`), the solution `x`, the step
  statistics without the timing entries, `event_times`, the PCNR counters.
- `pss`: `PSS.solve` and the other public PSS entry points in `PSS_METHODS`.
- `pac`: every public method of `PAC` (found through its MRO, so a split of
  `pac.py` into mixins keeps them).
  For both: the arguments (`inp`), the return value (`out`) and the public
  attributes the call (re)assigned on its instance (`state`), each FLATTENED
  to `{path: leaf}` -- arrays keep their dtype, pycircuit objects are never
  pickled (a moved class would not load), a solved PSS argument reads as
  `<PSS solve#k>` (the pss record that solved it) and other objects as a
  type marker.  ONLY THE OUTERMOST pss/pac call is recorded (one depth
  counter shared by both): public methods call each other (`pnoise` calls
  `adjoint_sideband_row`, `phase_psd` calls `coloured_diffusion_resolved`,
  PAC calls `pss.ppv` ...), and a refactor that reroutes those must not read
  as calls gone MISSING.

For every family the exception, if the call raised (re-raised unchanged),
and per test the warnings pytest recorded (category, message, location).

Usage (the plugin is opt-in, nothing in `pycircuit` imports it):

    PYTHONPATH=benchmarks/tranrec TRANREC_OUT=/tmp/rec_A \
    TRANREC_FAMILIES=all pytest pycircuit -q -p no:cacheprovider -p tran_recorder
    ... change the code, record again into /tmp/rec_B ...
    python benchmarks/tranrec/compare.py /tmp/rec_A /tmp/rec_B

`TRANREC_COVERAGE=1` also measures branch coverage of
`pycircuit.circuit.shooting` in every worker (data files in TRANREC_OUT;
`coverage combine` then `coverage report` there).

Each xdist worker writes `<TRANREC_OUT>/rec_<worker>.pkl` at session end,
~60 MB for the transient family alone: keep recordings OUT of the repo.
"""
import os, pickle, hashlib, inspect, functools, types
import numpy as np

_REC = []
_WARN = []
_CUR = {'id': None, 'k': {}, 'depth': 0, 'orbits': {}}
_TIMING = ('solve_seconds', 'total_seconds')
_COV = []

## the public PSS entry points recorded (`solve_timestep` and the
## `factored_period*` builders are plumbing, called per step / per replay)
PSS_METHODS = ('solve', 'find_initial_solution', 'event_grid', 'lte_grid',
               'ppv', 'frequency_aware_ppv', 'floquet_modes', 'grid_error',
               'warping_estimate')
BIG = 1_000_000          # arrays larger than this are stored as a digest
MAX_LEAVES = 20000       # per flattened value
MAX_DEPTH = 6


def _families():
    f = os.environ.get('TRANREC_FAMILIES', 'transient')
    if f.strip() == 'all':
        return {'transient', 'pss', 'pac'}
    return {x.strip() for x in f.split(',') if x.strip()}


def _next_k(fam):
    k = _CUR['k'].get(fam, 0)
    _CUR['k'][fam] = k + 1
    return k


def _sha(a):
    return hashlib.sha1(np.ascontiguousarray(a).tobytes()).hexdigest()


def _marker(v):
    """A string standing for an object that is not recorded by value, or
    None when `v` is to be flattened."""
    from pycircuit.circuit.analysis import Analysis
    from pycircuit.circuit.circuit import Circuit
    from pycircuit.circuit.shooting import PSS
    if isinstance(v, PSS):
        return '<PSS %s>' % _CUR['orbits'].get(id(v), 'unsolved-or-unrecorded')
    if isinstance(v, (Analysis, Circuit)):
        return '<%s>' % type(v).__name__
    if isinstance(v, (types.ModuleType, type, types.FunctionType,
                      types.BuiltinFunctionType, types.MethodType,
                      functools.partial)):
        return '<callable %s>' % getattr(v, '__qualname__', type(v).__name__)
    return None


def _flat(v, path='', out=None, depth=0, seen=None):
    """`v` as `{path: leaf}`; leaves are numpy arrays (copied, dtype kept),
    Python scalars, strings, None, or marker strings."""
    if out is None:
        out, seen = {}, set()
    if len(out) >= MAX_LEAVES:
        out['<truncated>'] = '<more than %d leaves>' % MAX_LEAVES
        return out
    if v is None or isinstance(v, (bool, int, float, complex, str, bytes)):
        out[path] = v
    elif isinstance(v, np.generic):
        out[path] = np.array(v)
    elif isinstance(v, np.ndarray):
        if v.dtype == object:
            for i, e in enumerate(v.ravel()[:MAX_LEAVES]):
                _flat(e, '%s[%d]' % (path, i), out, depth + 1, seen)
            out[path + '.shape'] = str(v.shape)
        elif v.size > BIG:
            fl = v.ravel()
            out[path] = ('__big__', _sha(v), v.shape, str(v.dtype),
                         np.array(fl[::max(1, fl.size // 4096)][:4096]))
        else:
            out[path] = np.array(v, copy=True)
    elif depth >= MAX_DEPTH:
        out[path] = '<depth %s>' % type(v).__name__
    elif isinstance(v, dict):
        for kk in sorted(v, key=repr):
            _flat(v[kk], '%s[%r]' % (path, kk), out, depth + 1, seen)
        if not v:
            out[path] = '<empty dict>'
    elif isinstance(v, (list, tuple)):
        if v and all(isinstance(e, (int, float, complex, np.number))
                     and not isinstance(e, bool) for e in v):
            out[path] = np.asarray(v)
            out[path + '.type'] = type(v).__name__
        else:
            for i, e in enumerate(v):
                _flat(e, '%s[%d]' % (path, i), out, depth + 1, seen)
            if not v:
                out[path] = '<empty %s>' % type(v).__name__
    else:
        m = _marker(v)
        if m is not None:
            out[path] = m
        elif id(v) in seen:
            out[path] = '<seen %s>' % type(v).__name__
        else:
            seen.add(id(v))
            ## `object.__getstate__`, not `vars`: `ResultDict` defines a
            ## METHOD named `__dict__`, so `vars()` of a PSS result is a method
            st = object.__getstate__(v)
            fields = {}
            for part in (st if isinstance(st, tuple) else (st,)):
                if isinstance(part, dict):
                    fields.update(part)
            fields = {kk: x for kk, x in fields.items()
                      if not kk.startswith('__')}
            if not fields:
                out[path] = '<%s>' % type(v).__name__
            for kk in sorted(fields):
                _flat(fields[kk], '%s.%s' % (path, kk), out, depth + 1, seen)
    return out


def _safe_flat(v):
    """`_flat`, which must never fail the test it records: a failure is
    recorded as its (deterministic) message."""
    try:
        return _flat(v)
    except Exception as e:                    # noqa: BLE001
        return {'<flatten error>': '%s: %s' % (type(e).__name__, str(e)[:200])}


def _instance_dict(obj):
    st = object.__getstate__(obj)
    for part in (st if isinstance(st, tuple) else (st,)):
        if isinstance(part, dict):
            return part
    return {}


def _public_ids(obj):
    return {n: id(x) for n, x in _instance_dict(obj).items() if not n.startswith('_')}


def _wrap(owner, name, fam):
    static = inspect.getattr_static(owner, name)
    kind = type(static)
    fn = static.__func__ if isinstance(static, (staticmethod, classmethod)) else static
    if getattr(fn, '_tranrec', False):
        return
    qual = '%s.%s' % (owner.__name__, name)
    bound = kind is not staticmethod

    @functools.wraps(fn)
    def wrapper(*a, **kw):
        if _CUR['id'] is None or _CUR['depth'] > 0:
            _CUR['depth'] += 1
            try:
                return fn(*a, **kw)
            finally:
                _CUR['depth'] -= 1
        k = _next_k(fam)
        self_ = a[0] if bound and a else None
        rec = {'test': _CUR['id'], 'fam': fam, 'name': qual, 'k': k,
               'inp': _safe_flat((a[1:] if bound else a, kw))}
        before = (_public_ids(self_) if self_ is not None
                  and kind is not classmethod else None)
        _CUR['depth'] += 1
        try:
            res = fn(*a, **kw)
        except BaseException as e:           # record, then re-raise unchanged
            rec['exc'] = '%s: %s' % (type(e).__name__, str(e)[:400])
            _REC.append(rec)
            raise
        finally:
            _CUR['depth'] -= 1
        rec['out'] = _safe_flat(res)
        if before is not None:
            changed = {n: x for n, x in _instance_dict(self_).items()
                       if not n.startswith('_') and before.get(n) != id(x)}
            rec['state'] = _safe_flat(changed)
        if fam == 'pss' and name == 'solve':
            _CUR['orbits'][id(self_)] = 'solve#%d' % k
        _REC.append(rec)
        return res

    wrapper._tranrec = True
    if kind is staticmethod:
        wrapper = staticmethod(wrapper)
    elif kind is classmethod:
        wrapper = classmethod(wrapper)
    setattr(owner, name, wrapper)


def _public_methods(cls, package):
    names = []
    for c in cls.__mro__:
        if not c.__module__.startswith(package):
            continue
        for n, v in vars(c).items():
            if n.startswith('_') or n in names:
                continue
            if isinstance(v, (staticmethod, classmethod)) or inspect.isfunction(v):
                names.append(n)
    return names


def _record_transient(orig):
    def solve(self, *a, **kw):
        k = _next_k('transient')
        rec = {'test': _CUR['id'], 'fam': 'transient', 'name': 'Transient.solve',
               'k': k, 'cls': type(self).__name__}
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
    return solve


def pytest_configure(config):
    fams = _families()
    if 'transient' in fams:
        from pycircuit.circuit import transient as T
        if not getattr(T.Transient.solve, '_tranrec', False):
            T.Transient.solve = _record_transient(T.Transient.solve)
    if 'pss' in fams:
        from pycircuit.circuit.shooting import PSS
        for n in PSS_METHODS:
            _wrap(PSS, n, 'pss')
    if 'pac' in fams:
        from pycircuit.circuit.shooting import PAC
        for n in _public_methods(PAC, 'pycircuit.circuit.shooting'):
            _wrap(PAC, n, 'pac')
    ## ⚠ Coverage starts AFTER pycircuit (and with it numpy / scipy) is
    ## imported: numpy.fft's C extension, first imported under coverage,
    ## fails with "cannot load module more than once per process" -- which a
    ## run on one test file never showed (its conftest had imported
    ## everything), and `pytest pycircuit` did.  So import-time lines (`def`,
    ## class bodies) read as MISSED; function bodies and branches are measured.
    if os.environ.get('TRANREC_COVERAGE') and os.environ.get('TRANREC_OUT'):
        import pycircuit.circuit.shooting               # noqa: F401
        import coverage
        out = os.environ['TRANREC_OUT']
        os.makedirs(out, exist_ok=True)
        _COV.append(coverage.Coverage(
            data_file=os.path.join(out, '.coverage'), data_suffix=True,
            branch=True, source=['pycircuit.circuit.shooting']))
        _COV[0].start()


def pytest_runtest_setup(item):
    _CUR['id'] = item.nodeid
    _CUR['k'] = {}
    _CUR['depth'] = 0
    _CUR['orbits'] = {}


def pytest_warning_recorded(warning_message, when, nodeid, location):
    ## only in the process that ran the test (the xdist controller sees the
    ## forwarded copies with `_CUR['id']` None)
    if when != 'runtest' or _CUR['id'] is None or nodeid != _CUR['id']:
        return
    w = warning_message
    _WARN.append({'test': nodeid, 'cat': w.category.__name__,
                  'msg': str(w.message)[:500],
                  'where': '%s:%s' % (os.path.basename(str(w.filename)), w.lineno)})


def pytest_sessionfinish(session, exitstatus):
    for cov in _COV:
        cov.stop()
        cov.save()
    out = os.environ.get('TRANREC_OUT')
    if not out:
        return
    os.makedirs(out, exist_ok=True)
    wid = os.environ.get('PYTEST_XDIST_WORKER', 'main')
    with open(os.path.join(out, 'rec_%s.pkl' % wid), 'wb') as f:
        pickle.dump({'version': 2, 'families': sorted(_families()),
                     'calls': _REC, 'warnings': _WARN}, f)
