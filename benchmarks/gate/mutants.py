"""HOW MUCH WOULD WE NOTICE?  A mutation sample of the transient's C units, or named changes, each run against a gate's own record:

    python benchmarks/gate/mutants.py sample --out DIR --ref GATE_DIR --wt WORKTREE [--n 60] [--seed 288]
    python benchmarks/gate/mutants.py spec FILE --out DIR --ref GATE_DIR --wt WORKTREE [--reach KEY]

A line changed in a worktree at this commit (the tree itself untouched), built in a store of its own, and run on the gate's matrix against the C rows of `GATE_DIR`: the Python's side
does not move, so it need not run again.  Decks fastest first, stopping at the first difference.  A status per mutant:
    'invalid'          does not compile: not a measurement, replaced in a sample;
    'killed:default'   a result differs on a deck by default (or it crashed, or hung, or ran away with memory);
    'killed:forced'    only under a forced setting;
    'killed:counts'    every result the same and a path counter moved: what tree against tree reads;
    'killed:tests'     the decks saw nothing and the unit's own tests fail;
    'killed:static'    no behaviour noticed, and the compile's own check does (a dropped statement leaving a variable unused);
    'survived'         nothing noticed: every one to be READ: an equivalent change, or a gap.

⚠ A TEST THAT READS THE C'S TEXT CATCHES EVERYTHING AND MEANS NOTHING.  A pinned pair fails on any edit, by design, so with it in a unit's tests every mutant reaching them read as caught (the first
pilot: 7 of 13).  They are left out (`TEXT_TESTS`), and the control proves it: a change with no effect, in every unit, must pass that unit's tests.
"""
import argparse
import concurrent.futures
import json
import os
import re
import subprocess
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gate as G
import gcov_read
import mutate
import settings
import verdict as V

#: the units sampled, and the modules that load each: their tests are what a survivor meets next
UNITS = {'loop.c': ('_tran_loop_c', '_tran_tcloop_c', '_tran_stageloop_c'), 'newton.c': ('_tran_newton_c',), 'lte.c': ('_tran_lte_c',),
         'radau_tc.c': ('_tran_radau_tc',), 'core.c': ('_tran_core',), 'radau.c': ('_tran_radau_c',), 'pcnrdc.c': ('_pcnr_c',), 'dcnewton.c': ('_dc_newton_c',)}
#: a mutant's run: a tight memory cap (a deck takes well under this), and a time from the record
DECK_CAP = '4G'
TEST_TIMEOUT = 900
#: what reads the C's text or the doc's, not what the C does: no catch at all
TEXT_TESTS = ('test_pinned_pairs.py', 'test_runtime_map.py')
#: what compiles each unit and reads the compiler's warnings: a catch of its own kind, asked last
STATIC_TESTS = ('test_c_abi.py',)


def csrc_dir(tree):
    return os.path.join(tree, 'pycircuit', 'circuit', 'csrc')


def _tests_of(tree, unit):
    d = os.path.join(tree, 'pycircuit', 'circuit', 'tests')
    pat = re.compile('|'.join(UNITS[unit]))
    def names_it(f):
        with open(os.path.join(d, f)) as fh:
            return pat.search(fh.read())
    return sorted(os.path.join('pycircuit/circuit/tests', f) for f in os.listdir(d)
                  if f.startswith('test_') and f.endswith('.py') and f not in TEXT_TESTS + STATIC_TESTS and names_it(f))


def _order(unit, ref):
    """The record's runs, a default setting first, each fastest first; Radau's units need no Gear-2 run and no DC."""
    groups = [f'{m}/{s}' for m, s in settings.MATRIX]
    if unit in ('radau.c', 'radau_tc.c'):
        groups = [g for g in groups if g.startswith('radau/')]
    groups.sort(key=lambda g: not g.endswith('/default'))
    return [(g, sorted((r for r in ref if r['group'] == g), key=lambda r: r['row'].get('s') or 0)) for g in groups]


def _deck_run(r, wt, store):
    method, setting = r['group'].split('/')
    s = r['row'].get('s') or 1.0
    cmd = ['timeout', str(int(max(90, 8 * s))), 'systemd-run', '--user', '--scope', '-q', '-p', f'MemoryMax={DECK_CAP}', '-p', 'MemorySwapMax=0',
           G.PY, os.path.join(HERE, 'run_deck.py'), r['deck'], method, setting]
    return G._run(cmd, wt, {'PYCIRCUIT_HDL_CACHE_DIR': store}, int(max(120, 10 * s)))


def counts_moved(a, b):
    return V.count_moves(a.get('counts', {}), b.get('counts', {}))


def _run_group(rows, g, wt, store, jobs, seen):
    """A group's decks, stopping at the first difference: `(differences, counters moved)`."""
    stop = threading.Event()
    first, moved = [], []

    def one(r):
        if stop.is_set():
            return
        row = _deck_run(r, wt, store)
        what = V.differing(row, r['row'])
        if what:
            first.append((r['deck'], what, row.get('why')))
            stop.set()
            return
        for k in seen:
            seen[k] += row.get('counts', {}).get(k, 0)
        c = counts_moved(row, r['row'])
        if c:
            moved.append((r['deck'], g, c))
    with concurrent.futures.ThreadPoolExecutor(jobs) as ex:
        list(ex.map(one, rows))
    return first, moved


def run_decks(unit, ref, wt, store, jobs, reach=()):
    """`(status, what)`: the first difference, group by group; or what the counters did; 'none' where nothing moved."""
    moved = None
    seen = dict.fromkeys(reach, 0)
    for g, rows in _order(unit, ref):
        first, mv = _run_group(rows, g, wt, store, jobs, seen)
        if first:
            status = 'killed:default' if g.endswith('/default') else 'killed:forced'
            return status, {'group': g, 'deck': first[0][0], 'what': str(first[0][1])[:240], 'why': first[0][2]}
        if mv and moved is None:
            moved = mv[0]
    lost = [k for k, v in seen.items() if not v]
    if moved or lost:
        return 'killed:counts', {'moved': str(moved)[:300], 'reach_lost': lost}
    return 'none', None


def run_tests(unit, wt, store, limit=TEST_TIMEOUT, tests=None):
    """The unit's own tests: None where they pass, else what failed.  Over `limit` seconds (three times their time with nothing changed): a hang the change made."""
    tests = _tests_of(wt, unit) if tests is None else tests
    if not tests:
        return None
    env = dict(os.environ, PYTHONPATH=wt, PYCIRCUIT_HDL_CACHE_DIR=store, OPENBLAS_NUM_THREADS='1')
    p = subprocess.run(['timeout', str(int(limit)), G.PY, '-m', 'pytest', '-q', '-x', '-p', 'no:cacheprovider', '-n', '8', '--tb=no', '-rf', *tests],
                       cwd=wt, env=env, capture_output=True, text=True, check=False)
    if p.returncode == 0:
        return None
    if p.returncode == 124:
        return f'over {int(limit)} s'
    names = re.findall(r'^FAILED \S+::(\S+)', p.stdout, re.MULTILINE)
    return names[:3] or f'exit {p.returncode}'


def compiles(unit, wt, store):
    """The mutant built, in the store, under the key its glue will ask for: None, or why not."""
    code = ('import sys\nfrom pycircuit.circuit import _cbuild\n'
            f'src = _cbuild.csrc({unit!r})\n'
            'key = _cbuild.source_key(src, False)\n'
            '_cbuild.find_so(key) or _cbuild._build(src, _cbuild.so_path(key), False)\n')
    env = dict(os.environ, PYTHONPATH=wt, PYCIRCUIT_HDL_CACHE_DIR=store)
    p = subprocess.run([G.PY, '-c', code], cwd=wt, env=env, capture_output=True, text=True, timeout=300, check=False)
    return None if p.returncode == 0 else (p.stderr.strip().rsplit('\n', 1)[-1])[:160]


def _clean_store(store):
    for name in os.listdir(store):
        p = os.path.join(store, name)
        if name.startswith('c-') and not os.path.islink(p):
            os.unlink(p)


class Applied:
    """A file of the worktree with one change, restored whatever happens."""

    def __init__(self, path, change):
        self.path, self.change = path, change

    def __enter__(self):
        with open(self.path) as fh:
            self.text = fh.read()
        with open(self.path, 'w') as fh:
            fh.write(self.change(self.text))

    def __exit__(self, *exc):
        with open(self.path, 'w') as fh:
            fh.write(self.text)


def one_mutant(m, change, a, ref, store, log):
    unit, t0 = m['unit'], time.perf_counter()
    path = os.path.join(a.wt, m.get('dir', 'pycircuit/circuit/csrc'), m.get('file', unit))
    with Applied(path, change):
        why = compiles(unit, a.wt, store) if path.endswith('.c') else None
        if why:
            status, what = 'invalid', why
        else:
            status, what = run_decks(unit, ref, a.wt, store, a.jobs, a.reach or ())
            if status == 'none':
                failed = run_tests(unit, a.wt, store, a.test_limit.get(unit, TEST_TIMEOUT))
                static = None if failed else run_tests(unit, a.wt, store, tests=[os.path.join('pycircuit/circuit/tests', t) for t in STATIC_TESTS])
                status, what = ('killed:tests', failed) if failed else ('killed:static', static) if static else ('survived', None)
    _clean_store(store)
    rec = dict(m, status=status, what=what, s=round(time.perf_counter() - t0, 1))
    log.add(**rec)
    print(f"{unit}:{m.get('line', m.get('name'))} {m.get('kind', '')} -> {status} ({rec['s']} s) {str(what)[:110] if what else ''}", flush=True)
    return status


def _lines(a):
    """`{unit: {line: text}}`: the lines the compiler made code for, from the gate's notes."""
    notes = os.path.join(a.ref, 'notes')
    csrc = set(UNITS)
    U_lines, _arms = gcov_read.universe(notes, os.path.join(a.out, '_universe'), csrc)
    out = {}
    for unit in UNITS:
        with open(os.path.join(csrc_dir(a.wt), unit)) as fh:
            text = fh.read().split('\n')
        out[unit] = {ln: text[ln - 1] for u, ln in U_lines if u == unit and ln - 1 < len(text)}
    return out


def sample(a, ref, store, log):
    tally = {}
    for d in mutate.draw(_lines(a), a.n, a.seed).values():
        done = 0
        for m in d['drawn']:
            if done >= d['want']:
                break
            old, new, line = m['old'], m['new'], m['line']

            def change(text, line=line, old=old, new=new):
                ls = text.split('\n')
                assert ls[line - 1] == old
                ls[line - 1] = new
                return '\n'.join(ls)
            st = one_mutant(m, change, a, ref, store, log)
            tally[st] = tally.get(st, 0) + 1
            done += st != 'invalid'
    return tally


def spec(a, ref, store, log):
    tally = {}
    for it in a.spec_items:
        def change(text, it=it):
            assert text.count(it['old']) == 1, it['name']
            return text.replace(it['old'], it['new'])
        st = one_mutant(dict(it, unit=it['unit']), change, a, ref, store, log)
        tally[st] = tally.get(st, 0) + 1
    return tally


def control(ref, a, store):
    """The worktree unchanged, on the record: nothing may move, or the record cannot be read against and every status would mean nothing.  The time it takes is
    a survivor's."""
    t0 = time.perf_counter()
    status, what = run_decks('newton.c', ref, a.wt, store, a.jobs)
    print(f'control, nothing changed: {status} ({round(time.perf_counter() - t0, 1)} s)', flush=True)
    if status != 'none':
        raise SystemExit(f'the record does not reproduce here: {status} {what}')
    a.test_limit = {}
    for unit in sorted({it['unit'] for it in a.spec_items} if a.spec_items else UNITS):
        t0 = time.perf_counter()
        failed = run_tests(unit, a.wt, store)
        took = time.perf_counter() - t0
        print(f'control, {unit} tests: {len(_tests_of(a.wt, unit))} files, {round(took, 1)} s, {failed or "pass"}', flush=True)
        if failed:
            raise SystemExit(f'the tests of {unit} fail with nothing changed: {failed}')
        a.test_limit[unit] = max(120, 3 * took)
        ## (a change with no effect: an empty statement opening the unit's first function.  Whatever fails on it reads the text, and would catch every mutant)
        with Applied(os.path.join(csrc_dir(a.wt), unit), lambda text: text.replace('\n{\n', '\n{\n    ;\n', 1)):
            why = compiles(unit, a.wt, store)
            failed = why or run_tests(unit, a.wt, store)
        _clean_store(store)
        print(f'control, {unit}, a change with no effect: {failed or "pass"}', flush=True)
        if failed:
            raise SystemExit(f'a change with no effect fails in {unit}: {failed}: this catches every mutant and means nothing')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('how', choices=('sample', 'spec'))
    ap.add_argument('file', nargs='?')
    ap.add_argument('--out', required=True)
    ap.add_argument('--ref', required=True, help="a gate's directory: its matrix.jsonl, its notes")
    ap.add_argument('--wt', required=True, help='a worktree at this commit, clean')
    ap.add_argument('--n', type=int, default=60)
    ap.add_argument('--seed', type=int, default=288)
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--reach', action='append')
    a = ap.parse_args()
    a.out, a.ref, a.wt = (os.path.abspath(x) for x in (a.out, a.ref, a.wt))
    dirty = subprocess.run(['git', 'status', '--short', '--', 'pycircuit'], cwd=a.wt, capture_output=True, text=True, check=False).stdout.strip()
    if dirty:
        raise SystemExit(f'the worktree is not clean: {dirty[:200]}')
    ref = [r for r in G._rows(a.ref, 'matrix.jsonl') if r['path'] == 'c' and not r['row'].get('failed')]
    store = G.seed_store(os.path.join(a.out, 'store'))
    log = G._Log(a.out, 'mutants.jsonl')
    a.spec_items = None
    if a.how == 'spec':
        with open(a.file) as fh:
            a.spec_items = json.load(fh)
    control(ref, a, store)
    tally = (sample if a.how == 'sample' else spec)(a, ref, store, log)
    print('tally:', tally)


if __name__ == '__main__':
    main()
