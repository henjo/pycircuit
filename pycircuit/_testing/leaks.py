"""THE LEAK DETECTOR (2026-10-03; robust testing, stage 2).

A test that changes process-global state and does not put it back makes a
LATER test -- often in another xdist worker's order, never in a single
process -- see a different world: a class bound to numpy where it should be
C, a switch off, the default toolkit symbolic.  The cases that motivated
this: gates G97/G98 found a MosLevel1 collapse variant bound to C with no
kernels in one worker of eight (numpy ran; one solve moved 1.4e-12), and the
suite's only reset (`reset_global_toolkit`) silently repaired 28 unrestored
toolkit writes.

Three checkpoints, so every kind of leak lands on the code that made it:

* PER TEST -- the `_pycircuit_leak_check` autouse fixture, function-scoped
  and the first such fixture to set up (it is a plugin's), so it snapshots
  after the module's fixtures and compares after the test's own teardown;
* PER MODULE -- after the last test of a module (its module- and class-
  scoped fixtures torn down), against the snapshot taken before that
  module's first test: a module fixture that does not restore, a module
  that mutates at import;
* SESSION -- the state at `pytest_sessionstart` against the first module's
  start: what collection (every test module's import) changed.

The snapshot is plain data: every scalar (bool/int/float/str/None) module
attribute of every loaded non-test `pycircuit` module -- a switch added
tomorrow is covered without registering it; the identity of every function
and class in those modules' dicts and in their classes' dicts (a replaced
method is a leak); the default toolkit; the default parameter set; the
environment (minus pytest's own); `sys.path`; the working directory;
numpy's error and print state; mpmath's precision; jax's x64 flag once jax
is loaded; the BLAS thread count; the hdl compiler's evaluation stacks and
session; and every hdl class's backend state (`pycircuit._testing.state`):
pins, status, bound flag, pending flag, the kernel of every chain function.

ALLOWED, because normal use does it: a lazy driver, status or memo going
from unloaded to loaded; an 'auto' class resolving to C at its first
instance (pending -> `'c'`, kernels None -> loaded); new modules, classes,
collapse variants; jax's x64 turning on when the jax toolkit is first
imported.  REPORTED although normal-looking: an 'auto' class resolving to
NUMPY -- on this box, with a compiler and cffi, that happens only when a
test has hidden one of them (the compile cache off, a fake compiler, cffi
missing), and the class then runs numpy for the rest of that process.

THE INVARIANT, checked at every checkpoint (`state.invariant_breaks`): a
class bound to C carries its kernels, one not bound carries none.  Each
broken class is reported once, at the first checkpoint that sees it.

`PYCIRCUIT_LEAKS` = `fail` (the default since stage 3: a leak fails the
test at teardown, or the module's last test), `warn` (report to a file; a
pure observer, the run is the run without it), `off`.  In `fail` mode the detector RESTORES only plain
data it owns the meaning of -- the known switches, the environment, the
toolkit, the VALUES of the default parameters, `sys.path`, the working
directory, numpy's and mpmath's state -- so one leak does not cascade;
everything else is reported, never touched.  (Restoring an object's
internal fields is not safe: the first warn-mode gate restored the default
`ParameterDict`'s name list and value dict but not its parameter table,
and 190 later tests failed on a parameter with no value.)  The toolkit is
restored even when `off`.  Reports go to
`$PYCIRCUIT_LEAKS_REPORT/leaks-<worker>.jsonl` (default
`test_timings/leaks/`).
"""
import json
import math
import os
import sys
import types

import pytest

MODE = os.environ.get('PYCIRCUIT_LEAKS', 'fail').strip().lower() or 'fail'

#: lazy loads and memos: may change from an unloaded value to anything
LAZY = {'STATUS', 'WALK_STATUS', '_driver', '_walk_driver', '_compiler', '_tmp_dir', '_GEN',
        '_PD', '_DEFAULT_EPAR', '_OWN'}
UNLOADED = (None, 'not loaded', False)
#: counters and identities that move on every use
IGNORED = {('pycircuit.circuit._hdl_cache', '_counter'),
           ## the PSF grammar's parser counts names as it reads them
           ('pycircuit.post.cds.yapps.runtime', 'in_name')}

#: the switches the detector restores (module, name): the inventory of
#: 2026-10-03 -- every other scalar is reported, never written
SWITCHES = (
    ('pycircuit.circuit._stamp_plan', 'ENABLED'),
    ('pycircuit.circuit._hdl_batch', 'ENABLED'),
    ('pycircuit.circuit._hdl_batch', 'SKIP_ZERO_SOURCE'),
    ('pycircuit.circuit._hdl_climit', 'ENABLED'),
    ('pycircuit.circuit._hdl_climit', 'WALK'),
    ('pycircuit.circuit._tran_core', 'CORE'),
    ('pycircuit.circuit._hdl_cse', 'ENABLED'),
    ('pycircuit.circuit._hdl_cse', 'FAST_ENABLED'),
    ('pycircuit.circuit._hdl_cse', 'FUSE_ENABLED'),
    ('pycircuit.circuit._hdl_cse', 'FUSE_MIN_CODE'),
    ('pycircuit.circuit._hdl_cse', 'MEMO_ENTRIES'),
    ('pycircuit.circuit._hdl_cache', 'ENABLED'),
    ('pycircuit.circuit.hdl', 'EMIT_C_SOURCE'),
    ('pycircuit.circuit.hdl', 'BACKEND'),
    ('pycircuit.circuit.hdl', 'LIMIT_PAR_CACHE'),
    ('pycircuit.circuit.hdl', 'PCNR_LIFT_AFFINE'),
    ('pycircuit.circuit.hdl', 'AUTOHOLD_MIN_OPS'),
    ('pycircuit.circuit.hdl', 'COMPILE_WARN_SECONDS'),
    ('pycircuit.circuit._limiting', 'CIRCUIT_LEVEL'),
    ('pycircuit.circuit._tran_companion', 'U_MEMO'),
    ('pycircuit.circuit._tran_predictor', 'PRED_WEIGHT_MEMO'),
    ('pycircuit.circuit._tran_predictor', 'PRED_WEIGHT_MEMO_SIZE'),
    ('pycircuit.circuit._tran_predictor', 'PRED_FAST'),
    ('pycircuit.circuit._tran_newton', 'C_KERNEL_SHARE'),
    ('pycircuit.circuit.linearsolver', 'MIN_N_FOR_KLU'),
    ('pycircuit.circuit.semiconductors', 'EXP_ARG_MAX'),
    ('pycircuit.circuit._ginac', 'MAX_COMPILE_CHARS'),
)

_MISSING = object()


def _test_module(name):
    return '.tests.' in name or name.endswith('.tests') or '.test.' in name or name.rsplit('.', 1)[-1].startswith('test_')


def _modules():
    return [(n, m) for n, m in list(sys.modules.items())
            if m is not None and (n == 'pycircuit' or n.startswith('pycircuit.'))
            and not _test_module(n) and not n.startswith('pycircuit._testing')]


def _scalar(v):
    return v is None or type(v) in (bool, int, float, str)


def _same(a, b):
    if a is b:
        return True
    if type(a) is float and type(b) is float and math.isnan(a) and math.isnan(b):
        return True
    try:
        r = a == b
    except Exception:  # noqa: BLE001 -- an object without a usable ==
        r = None
    return r if isinstance(r, bool) else repr(a) == repr(b)


_BLAS = None


def _blas_threads():
    global _BLAS
    if _BLAS is None:
        try:
            from threadpoolctl import ThreadpoolController
            _BLAS = ThreadpoolController()
        except Exception:  # noqa: BLE001
            _BLAS = False
    if not _BLAS:
        return None
    try:
        return tuple(sorted((c.internal_api, c.num_threads) for c in _BLAS.lib_controllers))
    except Exception:  # noqa: BLE001
        return None


def snapshot():
    """The state the detector compares, as plain data and identities."""
    from pycircuit._testing import state as _state
    scal, calls = {}, {}
    for n, m in _modules():
        d = getattr(m, '__dict__', None)
        if d is None:
            continue
        for k, v in list(d.items()):
            if k.startswith('__'):
                continue
            if _scalar(v):
                scal[(n, k)] = v
            elif isinstance(v, (types.FunctionType, type)):
                ## (every one, whoever defined it: a test that puts a lambda
                ## in place of a library function is the leak to catch)
                calls[(n, k)] = id(v)
                if isinstance(v, type) and v.__module__ == n:
                    for a, w in list(v.__dict__.items()):
                        if isinstance(w, (types.FunctionType, staticmethod, classmethod, property)):
                            calls[(n, k, a)] = id(w)
    snap = {'scalars': scal, 'callables': calls}
    cm = sys.modules.get('pycircuit.circuit.circuit')
    if cm is not None:
        tk = getattr(cm, 'default_toolkit', None)
        snap['toolkit'] = (id(tk), type(tk).__name__)
        ep = getattr(cm, 'defaultepar', None)
        if ep is not None:
            ## the VALUES of its parameters (a parameter registered on first
            ## use -- `analysis.analysis_kind` -- is the library's design,
            ## so names may appear; values may not move)
            snap['defaultepar'] = dict(vars(ep).get('_values', {}))
    snap['environ'] = {k: v for k, v in os.environ.items() if not k.startswith('PYTEST_')}
    ## the path as the import system resolves it: entries normalised, each
    ## once (a test that prepends the repository root under another
    ## spelling changes no import)
    seen, path = set(), []
    for p in sys.path:
        r = os.path.realpath(p) if isinstance(p, str) and p else p
        if r not in seen:
            seen.add(r)
            path.append(r)
    snap['sys.path'] = path
    snap['cwd'] = os.getcwd()
    np = sys.modules.get('numpy')
    if np is not None:
        snap['np.geterr'] = dict(np.geterr())
        snap['np.printoptions'] = {k: v for k, v in np.get_printoptions().items()}
    mp = sys.modules.get('mpmath')
    if mp is not None:
        snap['mpmath.dps'] = mp.mp.dps
    snap['jaxtoolkit'] = 'pycircuit.circuit._jaxtoolkit' in sys.modules
    jax = sys.modules.get('jax')
    if jax is not None:
        try:
            snap['jax.x64'] = bool(jax.config.jax_enable_x64)
        except AttributeError:
            snap['jax.x64'] = None
    snap['blas'] = _blas_threads()
    hdl = sys.modules.get('pycircuit.circuit.hdl')
    if hdl is not None:
        snap['hdl.stacks'] = (len(getattr(hdl, '_VAR_STACK', ())), len(getattr(hdl, '_NODE_TRACE', ())))
    eh = sys.modules.get('pycircuit.circuit._evalhint')
    if eh is not None:
        snap['evalhint'] = eh.current() is None
    snap['hdl'] = _state.hdl_state(ids=True)
    snap['hdl_funcs'] = {_state.class_key(c, True):
                         tuple(sorted((k, id(f)) for k, f in (c.__dict__['_hdl_info'].get('funcs') or {}).items()))
                         for c in _state.hdl_classes()}
    return snap


def _hdl_diff(a, b):
    out = []
    for name, sa in a['hdl'].items():
        sb = b['hdl'].get(name)
        if sb is None:
            continue
        resolved = (sa['pending'] and not sb['pending'])
        for key in ('pin', 'status', 'c_bound', 'pending', 'limit_status'):
            if sa[key] == sb[key]:
                continue
            if resolved and key in ('status', 'c_bound', 'pending', 'limit_status'):
                continue
            out.append(f'hdl {name}.{key}: {sa[key]!r} -> {sb[key]!r}')
        if resolved and sb['library'] and sb['status'] != 'c':
            out.append(f"hdl {name}: an 'auto' class resolved to {sb['status']!r} (not C: "
                       f'a compiler, cffi or the compile cache was hidden from it)')
        for f, fa in sa['funcs'].items():
            fb = sb['funcs'].get(f)
            if fb is None:
                out.append(f'hdl {name}.funcs[{f!r}] removed')
                continue
            if fa['kernel'] != fb['kernel'] and not (resolved and fa['kernel'] is None):
                out.append(f'hdl {name}.funcs[{f!r}] kernel {fa["kernel"]!r:.20} -> {fb["kernel"]!r:.20}')
        fa_ids, fb_ids = a['hdl_funcs'].get(name), b['hdl_funcs'].get(name)
        if fa_ids is not None and fb_ids is not None and fa_ids != fb_ids:
            changed = sorted({k for k, _ in set(fa_ids) ^ set(fb_ids)})
            out.append(f'hdl {name}: info["funcs"] entries replaced: {changed}')
    for name, sb in b['hdl'].items():
        if (name not in a['hdl'] and sb['library'] and sb['chained']
                and str(sb['status']).startswith('numpy (auto')):
            out.append(f"hdl {name} (new): an 'auto' class resolved to {sb['status']!r}")
    return out


def diff(a, b, skip=()):
    """The changes from snapshot `a` to `b` that are not allowed, as text;
    `skip` names whole sections to leave out."""
    out = []
    for key, va in a['scalars'].items():
        if key in IGNORED or key not in b['scalars']:
            continue
        vb = b['scalars'][key]
        if _same(va, vb):
            continue
        if key[1] in LAZY and va in UNLOADED:
            continue
        out.append(f'{key[0]}.{key[1]}: {va!r} -> {vb!r}')
    for key, ia in a['callables'].items():
        ib = b['callables'].get(key)
        if ib is not None and ib != ia and key[-1] not in LAZY:
            out.append(f'{".".join(key)} replaced')
    da, db = a.get('defaultepar', {}), b.get('defaultepar', {})
    for k in sorted(set(da) & set(db)):
        if not _same(da[k], db[k]):
            out.append(f'defaultepar.{k}: {_short(da[k])} -> {_short(db[k])}')
    for key in ('toolkit', 'sys.path', 'cwd', 'np.geterr', 'np.printoptions',
                'mpmath.dps', 'blas', 'hdl.stacks', 'evalhint'):
        if key in skip:
            continue
        if key in a and key in b and not _same(a[key], b[key]):
            out.append(f'{key}: {_short(a[key])} -> {_short(b[key])}')
    first_jax = b.get('jax.x64') and not a.get('jaxtoolkit') and b.get('jaxtoolkit')
    if ('jax.x64' in a and 'jax.x64' in b and a['jax.x64'] != b['jax.x64']
            and not first_jax):
        out.append(f"jax.x64: {a['jax.x64']} -> {b['jax.x64']}")
    ea, eb = a['environ'], b['environ']
    for k in sorted(set(ea) | set(eb)):
        if ea.get(k) != eb.get(k):
            out.append(f'env {k}: {ea.get(k)!r} -> {eb.get(k)!r}')
    out += _hdl_diff(a, b)
    return out


def _short(v, n=120):
    r = repr(v)
    return r if len(r) <= n else r[:n] + '...'


def restore(a):
    """Put back the plain data the detector owns the meaning of."""
    for mod, name in SWITCHES:
        m = sys.modules.get(mod)
        if m is not None and (mod, name) in a['scalars']:
            setattr(m, name, a['scalars'][(mod, name)])
    cm = sys.modules.get('pycircuit.circuit.circuit')
    if cm is not None and 'toolkit' in a:
        restore_toolkit(a)
        ep = getattr(cm, 'defaultepar', None)
        vals = vars(ep).get('_values') if ep is not None else None
        if vals is not None and 'defaultepar' in a:
            ## values only, of parameters present in both: never the dict's
            ## structure (a half-restored ParameterDict -- names without
            ## values -- broke 190 tests in the first warn-mode gate)
            for k, v in a['defaultepar'].items():
                if k in vals and not _same(vals[k], v):
                    vals[k] = v
    env = a['environ']
    for k in [k for k in os.environ if not k.startswith('PYTEST_') and k not in env]:
        del os.environ[k]
    for k, v in env.items():
        if os.environ.get(k) != v:
            os.environ[k] = v
    if [os.path.realpath(p) if isinstance(p, str) and p else p for p in sys.path] \
            != a['sys.path']:
        keep = set(a['sys.path'])
        sys.path[:] = [p for p in sys.path
                       if (os.path.realpath(p) if isinstance(p, str) and p else p) in keep]
    if os.getcwd() != a['cwd']:
        os.chdir(a['cwd'])
    np = sys.modules.get('numpy')
    if np is not None and 'np.geterr' in a:
        np.seterr(**a['np.geterr'])
        np.set_printoptions(**a['np.printoptions'])
    mp = sys.modules.get('mpmath')
    if mp is not None and 'mpmath.dps' in a:
        mp.mp.dps = a['mpmath.dps']


_TOOLKIT_OBJ = {}


def restore_toolkit(a):
    cm = sys.modules.get('pycircuit.circuit.circuit')
    obj = _TOOLKIT_OBJ.get(a.get('toolkit', (None,))[0])
    if cm is not None and obj is not None and cm.default_toolkit is not obj:
        cm.default_toolkit = obj


def _remember_toolkit():
    cm = sys.modules.get('pycircuit.circuit.circuit')
    if cm is not None:
        tk = cm.default_toolkit
        _TOOLKIT_OBJ[id(tk)] = tk


# ---------------------------------------------------------------------------
# reporting

_REPORTED_BREAKS = set()


def _report_dir(config):
    d = os.environ.get('PYCIRCUIT_LEAKS_REPORT')
    if d:
        return d
    return os.path.join(str(config.rootpath), 'test_timings', 'leaks')


def _worker(config):
    return getattr(config, 'workerinput', {}).get('workerid', 'main')


def _run_id():
    return os.environ.get('PYCIRCUIT_LEAKS_RUN', 'run')


def pytest_configure(config):
    ## one id per run, set by the controller before the workers start (they
    ## inherit the environment), so a run's report files are told apart
    if not hasattr(config, 'workerinput') and 'PYCIRCUIT_LEAKS_RUN' not in os.environ:
        import time
        os.environ['PYCIRCUIT_LEAKS_RUN'] = time.strftime('%Y%m%dT%H%M%S') + f'-{os.getpid()}'


def _record(config, where, kind, changes):
    rec = {'where': where, 'kind': kind, 'changes': changes, 'worker': _worker(config)}
    try:
        d = _report_dir(config)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f'leaks-{_run_id()}-{_worker(config)}.jsonl'), 'a') as f:
            f.write(json.dumps(rec) + '\n')
    except OSError:
        pass
    config._pycircuit_leaks = getattr(config, '_pycircuit_leaks', 0) + 1


def note_toolkit_reset(config, nodeid, old):
    """`reset_global_toolkit` (tests/conftest.py) repairs a toolkit a test
    left behind BEFORE the detector looks (it is torn down first): it reports
    what it repaired here, until stage 3 removes it."""
    if MODE != 'off':
        msg = f'toolkit: left as {type(old).__name__} (repaired by reset_global_toolkit)'
        _record(config, nodeid, 'test', [msg])


def _breaks(where):
    from pycircuit._testing import state as _state
    new = {k: v for k, v in _state.invariant_breaks().items() if k not in _REPORTED_BREAKS}
    _REPORTED_BREAKS.update(new)
    return [f'INVARIANT {k}: {v} (first seen {where})' for k, v in new.items()]


# ---------------------------------------------------------------------------
# the pytest side

def pytest_sessionstart(session):
    if MODE == 'off':
        return
    ## THE MODULES THAT HOLD SWITCHES, imported now: a module first imported
    ## while the tests are collected has no value from before, so a test
    ## module that changes one at import could not be seen
    import importlib
    for mod in sorted({m for m, _ in SWITCHES} | {'pycircuit.circuit.circuit'}):
        try:
            importlib.import_module(mod)
        except ImportError:
            continue                                 # (an optional backend: ginac)
    _remember_toolkit()
    session.config._pycircuit_leak_session = snapshot()


@pytest.fixture(autouse=True)
def _pycircuit_leak_check(request):
    if MODE == 'off':
        cm = sys.modules.get('pycircuit.circuit.circuit')
        tk = cm.default_toolkit if cm is not None else None
        yield
        if cm is not None and cm.default_toolkit is not tk:
            cm.default_toolkit = tk
        return
    _remember_toolkit()
    before = snapshot()
    yield
    after = snapshot()
    changes = diff(before, after) + _breaks(request.node.nodeid)
    if changes:
        if MODE == 'fail':
            restore(before)
        _record(request.config, request.node.nodeid, 'test', changes)
        if MODE == 'fail':
            pytest.fail('LEAK: this test changed global state and did not restore it:\n  '
                        + '\n  '.join(changes), pytrace=False)


_MODULE = {}
_SESSION_CHECKED = []


@pytest.hookimpl(hookwrapper=True, tryfirst=True)
def pytest_runtest_setup(item):
    if MODE != 'off':
        mod = item.nodeid.split('::', 1)[0]
        if _MODULE.get('name') != mod:
            _MODULE.update(name=mod, snap=snapshot())
            sess = getattr(item.config, '_pycircuit_leak_session', None)
            if not _SESSION_CHECKED and sess is not None:
                _SESSION_CHECKED.append(True)
                ## (pytest itself puts the test directories on sys.path while
                ## it collects)
                changes = diff(sess, _MODULE['snap'], skip=('sys.path',)) + _breaks('collection')
                if changes:
                    _record(item.config, 'collection (test module imports)', 'session', changes)
    yield


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_teardown(item, nextitem):
    yield
    if MODE == 'off':
        return
    mod = item.nodeid.split('::', 1)[0]
    if nextitem is None or nextitem.nodeid.split('::', 1)[0] != mod:
        before = _MODULE.get('snap')
        if before is not None and _MODULE.get('name') == mod:
            changes = diff(before, snapshot()) + _breaks(f'end of {mod}')
            if changes:
                if MODE == 'fail':
                    restore(before)
                _record(item.config, mod, 'module', changes)
                item._pycircuit_module_leak = changes
        _MODULE.clear()


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    if MODE != 'fail' or call.when != 'teardown':
        return
    changes = getattr(item, '_pycircuit_module_leak', None)
    if changes:
        rep = outcome.get_result()
        rep.outcome = 'failed'
        rep.longrepr = ('LEAK: module ' + item.nodeid.split('::', 1)[0]
                        + ' (its module/class fixtures or its import) changed global state:\n  '
                        + '\n  '.join(changes))


def pytest_terminal_summary(terminalreporter, config):
    if MODE == 'off' or hasattr(config, 'workerinput'):
        return
    d = _report_dir(config)
    if not os.path.isdir(d):
        return
    n = 0
    for name in os.listdir(d):
        if name.startswith(f'leaks-{_run_id()}-') and name.endswith('.jsonl'):
            with open(os.path.join(d, name)) as f:
                n += sum(1 for _ in f)
    if n:
        terminalreporter.write_line(f'pycircuit leak detector ({MODE}): {n} reports in {d}')
