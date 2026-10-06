"""What one function costs inside a real run, in instructions (speed round
12, stage 0, 2026-10-06; from speed round 11's scratch tools).

Four modes, each in ONE fresh process per target -- a wrap is the only
change to the run, and the run is a `step_machinery` case:

    python benchmarks/insitu.py rank  CASE [--top N]
        the case's functions ranked by self time and by calls (cProfile),
        each named as `module:Qual.name` -- the candidates to size
    python benchmarks/insitu.py calls CASE TARGET [TARGET ...]
        per call: the target's inclusive instructions, its calls, its share
        of the run -- the wrapper's own cost calibrated and taken off
    python benchmarks/insitu.py lines CASE TARGET
        per statement: the target's top-level statements each preceded by a
        checkpoint (an AST rewrite), the checkpoints' own cost calibrated
        and taken off; a compound statement carries everything inside it
    python benchmarks/insitu.py callers CASE TARGET
        who calls the target in the timed region, and how often (cProfile;
        nothing wrapped) -- the sites behind a call count

CASE: any `step_machinery` case (`stage`, `mos1`, `gp`, `psp`, `mos1_radau`,
`pss`, `pss_radau`, `vdp_pss`, `ladder_gear`, `ladder_radau`, `pnoise`, ...).
TARGET: `module:Qual.name` or a short name of `SHORT`.  `--tree DIR` runs
another tree (its package and its `step_machinery`).

WHY IN SITU.  A tree of inclusive timers ranks pieces but does not size them
(2026-10-03: a plan rested on one), and a standalone timing misses the cold
caches between a step's calls; the CPU's instruction counter
(`_bench.InstrCounter`), read while it runs, sizes a call where it is made.
Speed round 11 found its two largest wins this way (a tableau classified at
every step, `np.allclose` 131 k a call; small reductions by slice copies) --
neither was where anyone looked.

WHAT A WRAP CAN CHANGE, AND IS SHOWN.  A fast path that checks the
genuineness of the function it stands in for (`_paths.genuine`) declines
when that function is wrapped: the run takes another path.  Every count of
`_paths` is compared between an unwrapped run and the counted one, and any
difference is printed -- a sizing of a path the run no longer takes is not a
sizing.  A target bound by name elsewhere (`from m import f`) is wrapped
there too (`aliases`); one kept in a container (a record, a dict) is not,
and its calls go uncounted -- a call count of 0 says so.

THE RUNS.  The case once to warm it (builds, compiles, caches), once more
unwrapped (its total and its fast-path counts), then wrapped: once warm,
once counted.  What is counted is the run's TIMED call, the region
`step_machinery --count` counts (`_Region`): a case's setup -- its circuit
built, the PSS that pnoise reads -- is not (a first version counted the
whole case and read the circuit's construction as 12 % of mos1's solve).
Counts are instructions of the CPU this process is pinned to
(`_bench.pin_cpu`); a reading the counter did not take whole (multiplexed)
is rejected, not scaled.  Run on a quiet box: other processes do not move an
instruction count, but they move the heap layout that numpy's small calls
are sensitive to (~0.3 %).
"""
import argparse
import ast
import importlib
import inspect
import os
import struct
import subprocess
import sys
import textwrap

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

_R = 'pycircuit.circuit._tran_radau:_RadauStages.'
#: Short names for the targets speed rounds 11 and 12 sized
SHORT = {
    'pss_step': 'pycircuit.circuit.shooting._pss_inner:_InnerTransient.solve_timestep',
    'tr_step': 'pycircuit.circuit.transient:Transient.solve_timestep',
    'solve_ts_': 'pycircuit.circuit.transient:Transient._solve_timestep',
    'ts_rk': 'pycircuit.circuit._tran_stages:_SequentialStages._solve_timestep_rk',
    'rk_coupled': _R + '_rk_step_coupled',
    'rk_transformed': _R + '_rk_step_transformed',
    'context': _R + '_coupled_stage_context',
    'frozen': _R + '_radau_frozen',
    'finish': _R + '_finish_radau',
    'end_passes': _R + '_stage_end_passes',
    'tc_solve': 'pycircuit.circuit._tran_radau_tc:solve',
    'finish_stage': 'pycircuit.circuit._tran_stages:_SequentialStages._finish_stage_step',
    'pred_or': 'pycircuit.circuit._tran_predictor:_StagePredictor._pred_or',
    'branch': 'pycircuit.circuit._tran_branch:_BranchCheck._branch_after_coupled',
    'screen': 'pycircuit.circuit._tran_branch:_BranchCheck._branch_screen',
    'memo_get': 'pycircuit.circuit._tran_companion:_CompanionModel._memo_get',
    'memo_put': 'pycircuit.circuit._tran_companion:_CompanionModel._memo_put',
    'C_lookup': 'pycircuit.circuit._tran_companion:_CompanionModel._C_lookup',
    'source_at': 'pycircuit.circuit._tran_companion:_CompanionModel._source_at',
    'passes': 'pycircuit.circuit._tran_core:passes',
    'core_for': 'pycircuit.circuit._tran_core:core_for',
    'probe': 'pycircuit.circuit._tran_core:_Core.probe',
    'plan_for': 'pycircuit.circuit._stamp_plan:_plan_for',
    'ineligible': 'pycircuit.circuit._stamp_plan:_ineligible',
    'assemble_matrix': 'pycircuit.circuit._stamp_plan:assemble_matrix',
    'assemble_vector': 'pycircuit.circuit._stamp_plan:assemble_vector',
    'pd_getattr': 'pycircuit.utilities.param:ParameterDict.__getattr__',
    'newton': 'pycircuit.circuit._tran_newton:_StepNewton._newton',
    'newton_tol': 'pycircuit.circuit._tran_newton:_StepNewton._newton_tolerances',
    'roll_history': 'pycircuit.circuit._tran_history:_RunHistory._roll_history',
    'push_history': 'pycircuit.circuit._tran_history:_RunHistory._push_history',
    'sens_step': 'pycircuit.circuit.shooting._sens_c:step',
    'sens_map': 'pycircuit.circuit.shooting._sens_c:_map',
    'walk_stage': 'pycircuit.circuit.shooting._pss_walks:_PeriodWalks._walk_stage',
    'stage_reads': 'pycircuit.circuit.shooting._pss_walks:_stage_reads',
}


## -- the child: one target in one process ----------------------------------

def _setup(tree):
    """`_bench` and `step_machinery` of `tree`, this process pinned."""
    sys.path.insert(0, tree)
    sys.path.insert(1, os.path.join(tree, 'benchmarks'))
    import _bench
    _bench.pin_threads()
    _bench.pin_cpu()
    import step_machinery
    return _bench, step_machinery


def _resolve(target):
    """`(owner, name, the static attribute, the full name)`."""
    full = SHORT.get(target, target)
    mod, sep, qual = full.partition(':')
    if not sep or not qual:
        raise SystemExit(f'insitu: a target is module:Qual.name, not {target!r}')
    owner = importlib.import_module(mod)
    parts = qual.split('.')
    for p in parts[:-1]:
        owner = getattr(owner, p)
    return owner, parts[-1], inspect.getattr_static(owner, parts[-1]), full


def _function_of(static):
    if isinstance(static, property):
        return static.fget
    if isinstance(static, (staticmethod, classmethod)):
        return static.__func__
    return static


def _rebind(owner, name, static, fn):
    """Put `fn` where `static` was, as the same kind of attribute."""
    if isinstance(static, property):
        setattr(owner, name, property(fn, static.fset, static.fdel, static.__doc__))
    elif isinstance(static, staticmethod):
        setattr(owner, name, staticmethod(fn))
    elif isinstance(static, classmethod):
        setattr(owner, name, classmethod(fn))
    else:
        setattr(owner, name, fn)


def _aliases(fn, new):
    """Every module-level name in a pycircuit module bound to `fn` (a `from m
    import f`), rebound to `new`: its names."""
    out = []
    for mname, m in list(sys.modules.items()):
        if m is None or not mname.startswith('pycircuit'):
            continue
        for k, v in list(vars(m).items()):
            if v is fn:
                setattr(m, k, new)
                out.append(f'{mname}.{k}')
    return out


class _Region:
    """The run's timed call -- the region `step_machinery --count` counts --
    with hooks for entering and leaving it: a case's setup (its circuit
    built, the PSS pnoise reads solved once) is outside it.  `step_machinery`
    calls `_timed` by its module's name, so the patch reaches every case."""

    def __init__(self, SM, fd):
        self.armed = False
        self.enter = self.leave = None
        self.total = 0
        self.whole = True
        self.fd = fd
        orig = SM._timed

        def timed(fn):
            def region():
                if not self.armed:
                    return fn()
                v0 = struct.unpack('QQQ', os.read(self.fd, 24))
                if self.enter:
                    self.enter()
                try:
                    return fn()
                finally:
                    if self.leave:
                        self.leave()
                    v1 = struct.unpack('QQQ', os.read(self.fd, 24))
                    self.total += v1[0] - v0[0]
                    ## (the counter ran the whole region, or the reading is rejected)
                    self.whole = self.whole and (v1[1] - v0[1]) == (v1[2] - v0[2])
            return orig(region)
        SM._timed = timed

    def run(self, SM, case, enter=None, leave=None):
        """One counted run: the region's instructions (None: rejected)."""
        self.enter, self.leave, self.total, self.whole = enter, leave, 0, True
        self.armed = True
        try:
            SM.run_case(case)
        finally:
            self.armed = False
        return self.total if self.whole else None


def _reference(SM, case, _paths, region):
    """The case's warm run, then its reference: `(the region's instructions
    or None, the run's fast-path counts)`."""
    SM.run_case(case)
    b0 = _paths.snapshot()
    total = region.run(SM, case)
    return total, _paths.since(b0)


def _path_changes(ref, got):
    keys = sorted(set(ref) | set(got))
    return [(k, ref.get(k, 0), got.get(k, 0)) for k in keys
            if not k.startswith('once:') and ref.get(k, 0) != got.get(k, 0)]


def _report_paths(ref, got):
    ch = _path_changes(ref, got)
    if not ch:
        print('  fast paths: as the unwrapped run (every count)')
        return
    print(f'  ⚠ fast paths CHANGED by the wrap ({len(ch)} counts) -- this sizes another path:')
    for k, a, b in ch[:12]:
        print(f'      {k}: {a} -> {b}')


def _of_run(x, total):
    return f'{x / 1e6:.2f} M' + (f' of {total / 1e6:.1f} M ({100 * x / total:.2f} %)'
                                 if total else ' (the run\'s total rejected)')


def child_calls(tree, case, target):
    _bench, SM = _setup(tree)
    from pycircuit.circuit import _paths
    ic = _bench.InstrCounter()
    ic.start()
    region = _Region(SM, ic.fd)
    total, ref = _reference(SM, case, _paths, region)
    owner, name, static, full = _resolve(target)
    fn = _function_of(static)
    rd, up, fd = os.read, struct.unpack, ic.fd
    ## [counting, calls, instructions, inside the target]
    box = [False, 0, 0, False]

    def make(f):
        def counted(*a, **k):
            if not box[0] or box[3]:
                return f(*a, **k)
            box[3] = True
            v0 = up('QQQ', rd(fd, 24))[0]
            try:
                return f(*a, **k)
            finally:
                v1 = up('QQQ', rd(fd, 24))[0]
                box[3] = False
                box[1] += 1
                box[2] += v1 - v0
        counted.__wrapped__ = f
        counted.__name__ = getattr(f, '__name__', 'counted')
        counted.__qualname__ = getattr(f, '__qualname__', 'counted')
        return counted

    ## the wrapper's own cost: the same wrapper around a function that does
    ## nothing, read the same way
    noop = make(lambda: None)
    box[0] = True
    for _ in range(5000):
        noop()
    floor = box[2] / box[1]
    box[:] = [False, 0, 0, False]
    wrapped = make(fn)
    _rebind(owner, name, static, wrapped)
    al = _aliases(fn, wrapped)
    SM.run_case(case)                                   # warm, wrapped
    b0 = _paths.snapshot()

    def on():
        box[0] = True

    def off():
        box[0] = False
    ok = region.run(SM, case, on, off) is not None
    got = _paths.since(b0)
    calls, instr = box[1], box[2] - floor * box[1]
    print(f'{case}  {full}')
    if not ok:
        print('  ⚠ the counter was multiplexed during the run: rejected')
    print(f'  calls {calls}   per call {instr / max(calls, 1) / 1e3:.2f} k   in the run '
          f'{_of_run(instr, total)}   (wrapper floor {floor:.0f})')
    steps = SM._RUN_STEPS[0]
    if steps:
        print(f'  per step ({steps} steps, the solve\'s own setup included): '
              f'{instr / steps / 1e3:.1f} k, {calls / steps:.2f} calls')
    if al:
        print('  aliases wrapped too: ' + ', '.join(al))
    if not calls:
        print('  ⚠ no call seen in the timed region: the target may be held in a container '
              '(a record, a dict), or run only in the case\'s setup')
    _report_paths(ref, got)


def _rewrite(fn, mod, owner):
    """`fn` recompiled with a checkpoint before each top-level statement and
    the open interval closed in a `finally`: `(new function, {label: text})`,
    or why it cannot be.  A method is compiled inside a class of its owner's
    name, so its private names (`self.__x`) mangle as they did."""
    code = fn.__code__
    if code.co_freevars:
        return f'it closes over {code.co_freevars} (a closure or zero-argument super())'
    if code.co_flags & (inspect.CO_GENERATOR | inspect.CO_COROUTINE | inspect.CO_ASYNC_GENERATOR):
        return 'a generator or coroutine'
    src = textwrap.dedent(inspect.getsource(fn))
    tree = ast.parse(src)
    fdef = tree.body[0]
    if not isinstance(fdef, ast.FunctionDef):
        return 'not a plain function definition'
    fdef.decorator_list = []
    lines = src.splitlines()
    text, body = {}, []
    stmts = list(fdef.body)
    if (stmts and isinstance(stmts[0], ast.Expr) and isinstance(stmts[0].value, ast.Constant)
            and isinstance(stmts[0].value.value, str)):
        body.append(stmts.pop(0))                    # the docstring
    inner = []
    for k, st in enumerate(stmts):
        L = st.lineno
        text[L] = lines[L - 1].strip()[:72]
        inner.append(ast.parse(f'{"__ck0__" if k == 0 else "__ck__"}({L})').body[0])
        inner.append(st)
    body.append(ast.Try(body=inner, handlers=[], orelse=[],
                        finalbody=[ast.parse('__end__()').body[0]]))
    fdef.body = body
    if isinstance(owner, type):
        cls = ast.parse(f'class {owner.__name__}:\n    pass').body[0]
        cls.body = [fdef]
        tree.body = [cls]
    ast.fix_missing_locations(tree)
    ## (the source's line 1 is the definition's first line in its file)
    ast.increment_lineno(tree, code.co_firstlineno - 1)
    ns = {}
    ## (the rewritten definition, run in its own module's namespace)
    exec(compile(tree, inspect.getsourcefile(fn), 'exec'), mod.__dict__, ns)  # noqa: S102
    new = ns[owner.__name__].__dict__[fdef.name] if isinstance(owner, type) else ns[fdef.name]
    return new, text


def child_lines(tree, case, target):
    _bench, SM = _setup(tree)
    from pycircuit.circuit import _paths
    ic = _bench.InstrCounter()
    ic.start()
    region = _Region(SM, ic.fd)
    total, ref = _reference(SM, case, _paths, region)
    owner, name, static, full = _resolve(target)
    fn = _function_of(static)
    mod = sys.modules[fn.__module__]
    rd, up, fd = os.read, struct.unpack, ic.fd
    st = {'on': False, 'last': None, 'v': 0}
    acc, hits = {}, {}

    def __ck0__(label):
        ## (the first statement: an interval an exception left open in an
        ## earlier call is dropped, never attributed)
        st['last'] = label
        st['v'] = up('QQQ', rd(fd, 24))[0]

    def __ck__(label):
        v = up('QQQ', rd(fd, 24))[0]
        last = st['last']
        if st['on'] and last is not None:
            acc[last] = acc.get(last, 0) + (v - st['v'])
            hits[last] = hits.get(last, 0) + 1
        st['last'] = label
        st['v'] = up('QQQ', rd(fd, 24))[0]

    def __end__():
        v = up('QQQ', rd(fd, 24))[0]
        last = st['last']
        if st['on'] and last is not None:
            acc[last] = acc.get(last, 0) + (v - st['v'])
            hits[last] = hits.get(last, 0) + 1
        st['last'] = None

    mod.__dict__['__ck0__'], mod.__dict__['__ck__'], mod.__dict__['__end__'] = __ck0__, __ck__, __end__
    out = _rewrite(fn, mod, owner)
    if isinstance(out, str):
        raise SystemExit(f'insitu: {full} cannot be rewritten: {out}')
    new, text = out
    new.__wrapped__ = fn
    ## the checkpoints' own cost: two back to back, the interval between them
    st['on'] = True
    for _ in range(3000):
        __ck__(-1)
        __ck__(-2)
    __end__()
    cal = (acc.pop(-1) + acc.pop(-2, 0)) / (hits.pop(-1) + hits.pop(-2, 0))
    st['on'] = False
    _rebind(owner, name, static, new)
    al = _aliases(fn, new)
    SM.run_case(case)                                   # warm, rewritten
    acc.clear()
    hits.clear()
    b0 = _paths.snapshot()

    def on():
        st['on'] = True

    def off():
        st['on'] = False
    region.run(SM, case, on, off)
    got = _paths.since(b0)
    first = fn.__code__.co_firstlineno
    calls = max(hits.values()) if hits else 0
    print(f'{case}  {full}   (checkpoint {cal:.0f} instructions, taken off each statement)')
    rows, tot = [], 0.0
    for L in sorted(acc):
        net = acc[L] - cal * hits[L]
        tot += net
        rows.append((L, hits[L], net))
    for L, h, net in rows:
        print(f'  {L + first - 1:5d} {h:7d}x {net / max(calls, 1) / 1e3:9.2f} k a call  '
              f'{text.get(L, "")}')
    print(f'  sum {tot / max(calls, 1) / 1e3:.1f} k a call over {calls} calls; in the run '
          f'{_of_run(tot, total)}')
    if al:
        print('  aliases rewritten too: ' + ', '.join(al))
    _report_paths(ref, got)


def child_rank(tree, case, top):
    import cProfile
    import pstats
    _bench, SM = _setup(tree)
    ic = _bench.InstrCounter()
    ic.start()
    region = _Region(SM, ic.fd)
    SM.run_case(case)
    pr = cProfile.Profile()
    region.run(SM, case, pr.enable, pr.disable)
    stats = pstats.Stats(pr).stats
    names = {}

    def qualname(path, line, func):
        if 'site-packages' in path:
            rel = path.split('site-packages' + os.sep, 1)[1]
            return f'{rel[:-3].replace(os.sep, ".")}:{func}'
        if not path.startswith(tree):
            return f'{os.path.basename(path)}:{func}'
        if path not in names:
            m = {}
            rel = os.path.relpath(path, tree)
            modname = rel[:-3].replace(os.sep, '.').removesuffix('.__init__')
            try:
                with open(path) as fh:
                    t = ast.parse(fh.read())
            except (OSError, SyntaxError):
                t = None

            def visit(node, prefix):
                for ch in ast.iter_child_nodes(node):
                    if isinstance(ch, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                        q = prefix + ch.name
                        lines = [ch.lineno] + [d.lineno for d in ch.decorator_list]
                        if not isinstance(ch, ast.ClassDef):
                            for ln in lines:
                                m[ln] = q
                        visit(ch, q + '.')
            if t is not None:
                visit(t, '')
            names[path] = (modname, m)
        modname, m = names[path]
        return f'{modname}:{m.get(line, func)}'

    rows = []
    for (path, line, func), (cc, nc, tt, ct, _callers) in stats.items():
        rows.append((tt, ct, nc, qualname(path, line, func)))
    tot = sum(r[0] for r in rows) or 1.0
    print(f'{case}: the top {top} by self time (cProfile: ranks, does not size)')
    for tt, ct, nc, q in sorted(rows, reverse=True)[:top]:
        print(f'  {100 * tt / tot:5.1f} %  self {tt * 1e3:8.1f} ms  incl {ct * 1e3:8.1f} ms  '
              f'{nc:8d} calls  {q}')
    print(f'{case}: the top {top} by calls')
    for tt, ct, nc, q in sorted(rows, key=lambda r: -r[2])[:top]:
        print(f'  {nc:8d} calls  self {tt * 1e3:8.1f} ms  {q}')


def child_callers(tree, case, target, top):
    import cProfile
    import pstats
    _bench, SM = _setup(tree)
    ic = _bench.InstrCounter()
    ic.start()
    region = _Region(SM, ic.fd)
    SM.run_case(case)
    _owner, name, static, full = _resolve(target)
    fn = _function_of(static)
    ## (numpy's functions are C dispatchers around the Python one: its code
    ## is the `__wrapped__`'s, else matched by its module's file and name)
    code = getattr(fn, '__code__', None) or getattr(getattr(fn, '__wrapped__', None),
                                                    '__code__', None)
    pr = cProfile.Profile()
    region.run(SM, case, pr.enable, pr.disable)
    stats = pstats.Stats(pr).stats
    found = None
    if code is not None:
        found = stats.get((code.co_filename, code.co_firstlineno, code.co_name))
    if found is None:
        src = getattr(sys.modules.get(getattr(fn, '__module__', ''), None), '__file__', None)
        for (path, _line, func), v in stats.items():
            if func == name and src and os.path.exists(path) and os.path.samefile(path, src):
                found = v
                break
    if found is None:
        raise SystemExit(f'insitu: {full} was not called in the timed region')
    _cc, nc, _tt, _ct, callers = found
    print(f'{case}  {full}: {nc} calls in the timed region; its callers:')
    rows = sorted(((v[1] if isinstance(v, tuple) else v, k) for k, v in callers.items()),
                  reverse=True)
    for n, (cf, cl, cfn) in rows[:top]:
        where = os.path.relpath(cf, tree) if cf.startswith(tree) else cf
        print(f'  {n:8d}  {where}:{cl} {cfn}')


## -- the parent: one child per target --------------------------------------

def _spawn(tree, args):
    env = dict(os.environ, PYTHONHASHSEED='0')
    r = subprocess.run([sys.executable, os.path.abspath(__file__), *args, '--tree', tree,
                        '--child'], env=env, check=False)
    return r.returncode


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('mode', choices=('rank', 'calls', 'lines', 'callers'))
    ap.add_argument('case')
    ap.add_argument('targets', nargs='*')
    ap.add_argument('--tree', default=ROOT)
    ap.add_argument('--top', type=int, default=30)
    ap.add_argument('--child', action='store_true', help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    tree = os.path.abspath(a.tree)
    if a.child:
        if a.mode == 'rank':
            child_rank(tree, a.case, a.top)
        elif a.mode == 'callers':
            child_callers(tree, a.case, a.targets[0], a.top)
        elif a.mode == 'calls':
            child_calls(tree, a.case, a.targets[0])
        else:
            child_lines(tree, a.case, a.targets[0])
        sys.stdout.flush()
        os._exit(0)
    if a.mode == 'rank':
        return _spawn(tree, ['rank', a.case, '--top', str(a.top)])
    if not a.targets:
        ap.error(f'{a.mode} needs a target')
    if a.mode in ('lines', 'callers') and len(a.targets) != 1:
        ap.error(f'{a.mode} takes one target')
    if a.mode == 'callers':
        return _spawn(tree, ['callers', a.case, a.targets[0], '--top', str(a.top)])
    rc = 0
    for t in a.targets:
        rc = _spawn(tree, [a.mode, a.case, t]) or rc
    return rc


if __name__ == '__main__':
    sys.exit(main())
