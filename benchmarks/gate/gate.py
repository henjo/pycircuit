"""THE GATE: what a change to the simulator must show, run and read in one place.

    python benchmarks/gate/gate.py all --out DIR --base WORKTREE [--since REV] [--reach KEY] [--line 'mesh 63x63 gear2 1e-6 >= 1.25']

`pairs`    timing lines: interleaved pairs against the base tree, medians, on a settled box (uninstrumented; run alone).
`suite`    the test suite (`--coverage`: under the coverage build, for the 'fixture' class).
`tree`     the gold decks and the circuits at size, this tree against the base: results the same, or where they moved.
`matrix`   the gold decks under the forced settings (`settings.MATRIX`), the C paths against the Python's, on the coverage build.
`size`     the same at size (`settings.SIZE_TWINS`), and the C paths further (`settings.SIZE_RUNS`).
`report`   the verdict (`verdict.py`): SAME, NOT REACHED or DIFFERS, with what each comparison reached; `--since REV` classes every
           changed C line and arm by what compared it; `--reach KEY` asks a path counter to have moved.
`all`      in that order, `report` last.

A run happens in the tree under test (cwd and PYTHONPATH) and says which package it imported: a gate that ran the wrong tree is a failed gate,
not a pass (`spice_suite` puts its own root first on `sys.path`, whatever PYTHONPATH says).  The coverage build has a store of its own under
DIR, its entries links to the default store's, so the HDL kernels stay warm and removing DIR removes what the build made.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import re
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import gcov_read
import settings
import verdict as V

PY = os.path.join(ROOT, '.venv', 'bin', 'python')
SIZE_CAP = '11G'


def _default_store():
    return os.environ.get('PYCIRCUIT_HDL_CACHE_DIR') or os.path.join(
        os.environ.get('XDG_CACHE_HOME') or os.path.expanduser('~/.cache'), 'pycircuit', 'hdl')


def seed_store(store):
    """A store of its own, its entries links to the default store's (the kernels' objects, the FFIs); what is built in it stays in it."""
    os.makedirs(store, exist_ok=True)
    src = _default_store()
    have = set(os.listdir(store))
    for name in os.listdir(src):
        if name != 'locks' and name not in have:
            os.symlink(os.path.join(src, name), os.path.join(store, name))
    return store


def _cov_env(out, scope):
    notes = os.path.join(out, 'notes')
    os.makedirs(notes, exist_ok=True)
    return {'PYCIRCUIT_C_COVERAGE': notes, 'PYCIRCUIT_HDL_CACHE_DIR': os.path.join(out, 'store'),
            'GCOV_PREFIX': os.path.join(out, 'scopes', scope), 'GCOV_PREFIX_STRIP': str(notes.count('/'))}


def _run(cmd, tree, env, timeout):
    """One runner: its JSON row, or a failed row saying why; a row from another tree than asked for is failed."""
    e = dict(os.environ, OPENBLAS_NUM_THREADS='1', PYTHONPATH=tree, **env)
    try:
        p = subprocess.run(cmd, cwd=tree, env=e, capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired:
        return {'failed': True, 'why': f'over {timeout} s'}
    last = p.stdout.strip().rsplit('\n', 1)[-1] if p.stdout.strip() else ''
    if not last.startswith('{'):
        return {'failed': True, 'why': (p.stderr or '').strip()[-300:] or f'exit {p.returncode}'}
    row = json.loads(last)
    if not os.path.realpath(row['file']).startswith(os.path.realpath(tree) + os.sep):
        return {'failed': True, 'why': f"ran {row['file']}, not {tree}"}
    return row


def _deck_cmd(deck, method, setting):
    return [PY, os.path.join(HERE, 'run_deck.py'), deck, method, setting]


def _size_cmd(spec):
    kind, size, method, tend = spec
    return ['timeout', '1500', 'systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={SIZE_CAP}', '-p', 'MemorySwapMax=0',
            PY, os.path.join(HERE, 'run_size.py'), kind, str(size), method, repr(tend)]


def _path_env(path, method):
    return {k: '0' for k in settings.PATHS[path].get(method, ())}


class _Log:
    """Rows appended as they come, so a stopped gate leaves what it made."""

    def __init__(self, out, name):
        os.makedirs(out, exist_ok=True)
        self.path = os.path.join(out, name)
        self.lock = threading.Lock()
        with open(self.path, 'w'):
            pass

    def add(self, **row):
        with self.lock, open(self.path, 'a') as fh:
            fh.write(json.dumps(row) + '\n')


def _pool(jobs, tasks):
    with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
        for f in concurrent.futures.as_completed([ex.submit(t) for t in tasks]):
            f.result()


def matrix(a):
    """The decks under each forced setting, on each path."""
    if a.coverage:
        seed_store(os.path.join(a.out, 'store'))
    log = _Log(a.out, 'matrix.jsonl')
    pairs = [tuple(s.split('/')) for s in a.settings.split(',')] if a.settings else settings.MATRIX
    decks = a.decks.split(',') if a.decks else settings.DECKS
    paths = a.paths.split(',')

    def one(deck, method, setting, path):
        def run():
            env = _path_env(path, method)
            if a.coverage:
                env.update(_cov_env(a.out, f'{method}_{setting}_{path}'))
            row = _run(_deck_cmd(deck, method, setting), a.tree, env, a.timeout)
            log.add(group=f'{method}/{setting}', deck=deck, path=path, row=row)
        return run
    tasks = [one(d, m, s, p) for d in decks for m, s in pairs for p in paths if p == 'c' or (m, s) not in settings.NO_TWIN]
    if a.coverage and tasks:
        tasks[0]()          # (the units built once, not eight times over the locks)
        tasks = tasks[1:]
    _pool(a.jobs, tasks)


def size(a):
    """The circuits at size: C against Python where the Python's side is affordable, the C's alone further."""
    if a.coverage:
        seed_store(os.path.join(a.out, 'store'))
    log = _Log(a.out, 'size.jsonl')

    def one(spec, path):
        def run():
            env = _path_env(path, spec[2])
            if a.coverage:
                env.update(_cov_env(a.out, f'size_{path}'))
            log.add(group=' '.join(map(str, spec)), deck=' '.join(map(str, spec[:2])), path=path, row=_run(_size_cmd(spec), a.tree, env, 1600))
        return run
    _pool(a.size_jobs, [one(s, p) for s in settings.SIZE_TWINS for p in ('c', 'py')] + [one(s, 'c') for s in settings.SIZE_RUNS])


def tree(a):
    """This tree against the base, uninstrumented: the decks by default, the circuits at size on the C paths."""
    log = _Log(a.out, 'tree.jsonl')

    def deck(d, method, setting, side, where):
        def run():
            log.add(group=method if setting == 'default' else f'{method}/{setting}', deck=d, path=side,
                    row=_run(_deck_cmd(d, method, setting), where, {}, a.timeout))
        return run

    def circ(spec, side, where):
        def run():
            log.add(group='size', deck=' '.join(map(str, spec)), path=side, row=_run(_size_cmd(spec), where, {}, 1600))
        return run
    sides = (('new', a.tree), ('base', a.base))
    groups = [('gear2', 'default'), ('radau', 'default'), ('dc', 'default'), *sorted(settings.NO_TWIN)]
    _pool(a.jobs, [deck(d, m, st, s, w) for d in settings.DECKS for m, st in groups for s, w in sides])
    _pool(a.size_jobs, [circ(c, s, w) for c in settings.SIZE_TWINS + settings.SIZE_RUNS for s, w in sides])


def _settle(limit=1.5, tries=90):
    for _ in range(tries):
        load = os.getloadavg()[0]
        if load < limit:
            return load
        time.sleep(10)
    return os.getloadavg()[0]


def pairs(a):
    """Interleaved pairs against the base for each timing line, base first: run alone, on a settled box."""
    log = _Log(a.out, 'pairs.jsonl')
    _settle()
    for line in a.line or ():
        kind, sz, method, tend, _op, _bound = line.split()
        spec = (kind, sz, method, float(tend))
        for _ in range(a.pairs):
            for side, where in (('base', a.base), ('new', a.tree)):
                log.add(line=line, side=side, row=_run(_size_cmd(spec), where, {}, 1600))


def suite(a):
    env = dict(os.environ, PYTHONPATH=a.tree)
    if a.coverage:
        seed_store(os.path.join(a.out, 'store'))
        env.update(_cov_env(a.out, 'suite'))
    os.makedirs(a.out, exist_ok=True)
    cmd = [PY, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', '-n', '8', '-rfEs', '--tb=short'] + list(a.pytest or ())
    p = subprocess.run(cmd, cwd=a.tree, env=env, capture_output=True, text=True, check=False)
    with open(os.path.join(a.out, 'suite.txt'), 'w') as fh:
        fh.write(p.stdout[-20000:])
    print(p.stdout.strip().rsplit('\n', 1)[-1])
    return p.returncode


def _rows(out, name):
    p = os.path.join(out, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return [json.loads(ln) for ln in fh]


def _coverage(a):
    """What each scope took, and each class: `None` where the gate ran without the coverage build."""
    notes, scopes = os.path.join(a.out, 'notes'), os.path.join(a.out, 'scopes')
    if not os.path.isdir(scopes):
        return None
    csrc = set(os.listdir(os.path.join(a.tree, 'pycircuit', 'circuit', 'csrc')))
    U_lines, U_arms = gcov_read.universe(notes, os.path.join(scopes, '_universe'), csrc)
    L, A = {}, {}
    for s in os.listdir(scopes):
        if s != '_universe':
            ls, ar = gcov_read.read(os.path.join(scopes, s), notes, csrc)
            L[s], A[s] = gcov_read.taken(ls), gcov_read.taken(ar)

    def union(table, want):
        return set().union(*[v for k, v in table.items() if want(k)]) if table else set()
    default = ('gear2_default_c', 'radau_default_c', 'dc_default_c')
    is_forced = lambda k: k.endswith('_c') and k not in default and not k.startswith('size')
    cl = {'lines': [union(L, lambda k: k in default), union(L, is_forced), L.get('size_c', set()), L.get('suite', set())],
          'arms': [union(A, lambda k: k in default), union(A, is_forced), A.get('size_c', set()), A.get('suite', set())]}
    return {'U_lines': U_lines, 'U_arms': U_arms, 'cl': cl, 'A': A, 'L': L}


def _class_table(cov):
    """Lines and arms by class, each counted once, in `V.CLASSES`' order; and the arms' table by file."""
    out = {}
    for what, U in (('lines', cov['U_lines']), ('arms', cov['U_arms'])):
        c = collections.Counter(V.arm_class(k, *cov['cl'][what]) for k in U)
        out[what] = {k: c[k] for k in V.CLASSES}
    byf = collections.defaultdict(collections.Counter)
    for k in cov['U_arms']:
        byf[k[0]][V.arm_class(k, *cov['cl']['arms'])] += 1
    return out, byf


def _changed(a, cov):
    diff = subprocess.run(['git', 'diff', '-U0', a.since, '--', 'pycircuit/circuit/csrc'], cwd=a.tree, capture_output=True, text=True, check=False).stdout
    arms_at = collections.defaultdict(list)
    for k in cov['U_arms']:
        arms_at[k[:2]].append(k)
    return V.changed_report(V.changed_lines(diff), cov['U_lines'],
                            lambda k: V.arm_class(k, *cov['cl']['lines']), lambda k: sorted(arms_at[k]),
                            lambda k: V.arm_class(k, *cov['cl']['arms']))


def _print_comparisons(title, comps):
    print(f'\n{title}')
    for g, c in sorted(comps.items(), key=lambda kv: str(kv[0])):
        print(f"  {g!s:34s} {c['runs']:3d} runs, {len(c['differ'])} differing")
        for d, what in c['differ'][:6]:
            print(f'      {d}: {str(what)[:200]}')


def _print_lines(a):
    rows = _rows(a.out, 'pairs.jsonl')
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        by[r['line']][r['side']].append(r['row'])
    missed = []
    if by:
        print('\nTiming lines (interleaved pairs, base first):')
    for line, sides in by.items():
        *_spec, op, bound = line.split()
        pairs_ = [(b, n) for b, n in zip(sides['base'], sides['new']) if not b.get('failed') and not n.get('failed')]
        rat = [b['solve'] / n['solve'] for b, n in pairs_]
        m, met = V.line_met(rat, op, float(bound))
        print(f"  {line}: {[round(x, 3) for x in rat]}, median {m and round(m, 3)}, peak {pairs_ and pairs_[0][0]['peak']} -> {pairs_ and pairs_[0][1]['peak']} MB: "
              f"{'MET' if met else 'MISSED'}")
        if not met:
            missed.append(line)
    return missed


def _instrument_check(matrix_rows, tree_rows):
    """The coverage build moves no answer: the matrix's default C rows against the uninstrumented `new` rows of `tree`."""
    new = {(r['group'], r['deck']): r['row'] for r in tree_rows if r['path'] == 'new' and r['group'] != 'size'}
    bad, n = [], 0
    for r in matrix_rows:
        m, s = r['group'].split('/')
        if s == 'default' and r['path'] == 'c' and (m, r['deck']) in new:
            n += 1
            if V.differing(r['row'], new[(m, r['deck'])]):
                bad.append((m, r['deck']))
    return n, bad


def report(a):
    mx, sz, tr = _rows(a.out, 'matrix.jsonl'), _rows(a.out, 'size.jsonl'), _rows(a.out, 'tree.jsonl')
    comps = {}
    twins = V.compare(V.with_side(mx, 'py'))
    _print_comparisons('C against Python, the decks:', twins)
    comps.update({f'twin {k}': v for k, v in twins.items()})
    off = V.compare(V.with_side(mx, 'off'), right='off')
    if off:
        _print_comparisons('C against every C path off:', off)
        comps.update({f'off {k}': v for k, v in off.items()})
    at_size = V.compare(V.with_side(sz, 'py'))
    _print_comparisons('C against Python, at size:', at_size)
    comps.update({f'size {k}': v for k, v in at_size.items()})
    trees = V.compare(tr, left='new', right='base')
    _print_comparisons('This tree against the base:', trees)
    moves = V.compare_counts(tr, allow=re.compile(a.moves) if a.moves else None)
    _print_comparisons('Its path counters' + (f' (bar {a.moves})' if a.moves else '') + ':', moves)
    if a.answers_move:
        print('  (--answers-move: a difference in these is reported, not a failure)')
    else:
        comps.update({f'tree {k}': v for k, v in trees.items()})
        comps.update({f'counters {k}': v for k, v in moves.items()})
    n, bad = _instrument_check(mx, tr)
    if n:
        print(f'\nThe instrument: {n} decks on the coverage build against without it, {len(bad)} differing {bad[:8]}')
    reached = V.reach(mx + sz, a.reach or ())
    if reached:
        print('\nReached on the compared C runs:', reached)
    named = [{'name': f"{r['row'].get('kind')} {r['row'].get('size')} {r['row'].get('method')}", 'n': r['row'].get('n'), 'swing': r['row'].get('swing')}
             for r in sz if r['path'] == 'c' and not r['row'].get('failed')]
    q = V.quiet(named)
    if named:
        print('\nSwing over 100 mV at size:', {x['name']: f"{x['swing']} of {x['n']}" for x in named}, 'QUIET:', q)
    changed = None
    cov = _coverage(a)
    if cov:
        tab, byf = _class_table(cov)
        print('\nC lines by what compared them, each counted once:', tab['lines'])
        print('C branch arms the same way:', tab['arms'], '(compiler-made arms included: "none" is an upper bound)')
        print(f"  {'file':16s}{'arms':>6s}" + ''.join(f'{c:>9s}' for c in V.CLASSES))
        for f, c in sorted(byf.items(), key=lambda kv: -sum(kv[1].values())):
            print(f"  {f:16s}{sum(c.values()):6d}" + ''.join(f'{c[k]:9d}' for k in V.CLASSES))
        if a.since:
            changed = _changed(a, cov)
            print(f'\nChanged since {a.since}: {len(changed["lines"])} executable lines, {len(changed["arms"])} arms:',
                  collections.Counter(c for _k, c in changed['arms']))
    missed = _print_lines(a)
    out, why = V.verdict(comps, reached, changed)
    print(f'\nVERDICT: {out}')
    for w in why[:40]:
        print(f'  {w}')
    if missed:
        print('TIMING LINES MISSED:', missed)
    return 0 if out == 'SAME' and not missed else 1


def all_(a):
    t = {}
    for name, step in (('pairs', pairs), ('suite', suite), ('tree', tree), ('matrix', matrix), ('size', size)):
        if name == 'pairs' and not a.line:
            continue
        if name == 'suite':
            cov, a.coverage = a.coverage, False
        t0 = time.perf_counter()
        step(a)
        t[name] = round(time.perf_counter() - t0, 1)
        if name == 'suite':
            a.coverage = cov
        print(f'{name} done, {t[name]} s', flush=True)
    with open(os.path.join(a.out, 'times.json'), 'w') as fh:
        json.dump(t, fh)
    print('times (s):', t)
    return report(a)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=('all', 'pairs', 'suite', 'tree', 'matrix', 'size', 'report'))
    ap.add_argument('--out', required=True)
    ap.add_argument('--tree', default=ROOT, help='the tree under test')
    ap.add_argument('--base', help='the tree to compare with (a worktree)')
    ap.add_argument('--decks', help='comma-separated: instead of every deck')
    ap.add_argument('--settings', help='comma-separated method/setting: instead of the matrix')
    ap.add_argument('--paths', default='c,py', help='comma-separated: c, py, off')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--size-jobs', type=int, default=3)
    ap.add_argument('--timeout', type=int, default=1200)
    ap.add_argument('--no-coverage', dest='coverage', action='store_false', default=True)
    ap.add_argument('--coverage-suite', action='store_true', help='`suite` under the coverage build')
    ap.add_argument('--line', action='append', help="a timing line: 'kind size method tend >= bound'")
    ap.add_argument('--pairs', type=int, default=2)
    ap.add_argument('--since', help='class the C lines and arms changed since this revision')
    ap.add_argument('--reach', action='append', help='a path counter the compared C runs must have moved')
    ap.add_argument('--moves', help='a regular expression: the path counters this change moves on purpose')
    ap.add_argument('--answers-move', action='store_true', help='the change moves answers on purpose: tree against tree is reported, not judged')
    ap.add_argument('pytest', nargs='*', help='with `suite`: arguments for pytest')
    a = ap.parse_args()
    a.out = os.path.abspath(a.out)
    a.tree = os.path.abspath(a.tree)
    a.base = a.base and os.path.abspath(a.base)
    if a.step == 'suite':
        a.coverage = a.coverage_suite
    sys.exit({'all': all_, 'pairs': pairs, 'suite': suite, 'tree': tree, 'matrix': matrix, 'size': size, 'report': report}[a.step](a) or 0)


if __name__ == '__main__':
    main()
