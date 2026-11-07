"""EVERY HAND-BACK, FIRED.  `python benchmarks/gate/handback.py --out DIR --ref GATE_DIR --wt WORKTREE`

A hand-back's whole contract is that the run goes on as if the C had never tried: nothing of the attempt left, the Python repeating it.  Most of the one-line checks that make one never fire by themselves (a census
of 866a0300: up to 109 of 154 in any test at all), so nothing has ever tested that contract there.  Here each is forced true at a chosen asking, in a worktree (`mutate.force_guards`; the tree itself untouched) on a
deck's run, and the run must equal the gate's record of the same deck and setting.  For each check, on the settings whose gcov record shows it asked:
    'same'        it fired, and every result and count is the record's;
    'counts'      it fired, the results the record's, a KLU count not;
    'differs'     it fired and a result moved, or the run crashed or hung: read it.  A status a C unit acts on itself is not a hand-back, and forcing it takes another course, not a wrong one: the error test's
                  solve under the early check, where a failed solve means no abort (`early`, newton.c) as the Python's own failing solve does; forcing one side alone says something the two cannot disagree on;
    'not fired'   nothing here reached it (named: the shooting's units, and what no deck meets).
The control: nothing forced, on every deck by default, must be the record.
"""
import argparse
import collections
import concurrent.futures
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gate as G
import gcov_read
import mutate
import settings
import verdict as V

#: units a gold deck can reach; the shooting's (`stagewalk.c`, `sens*.c`, `lmmsens.c`) are their tests' to fire
UNITS = ('branch.c', 'core.c', 'ctrl.c', 'dcnewton.c', 'lte.c', 'loop.c', 'newton.c', 'pcnrdc.c', 'pred.c', 'radau.c', 'radau_fold.c', 'radau_tc.c', 'sources.c', 'walk.c')
RESULTS = ('x', 'err', 'steps', 'rej', 'newton', 'fb')
#: how many decks a check is fired on, at each of these askings (not one fixed gap: plants at a fixed gap resonate)
DECKS_A_CHECK = 2
ASKINGS = (1, 3, 17)


def _env_name(unit):
    return 'PCX_FAULT_' + unit[:-2].upper()


def _instrument(wt):
    """Every unit's checks forced, in the worktree: `{unit: {line: text}}`, `{unit: original text}`."""
    sites, originals = {}, {}
    for unit in UNITS:
        p = os.path.join(wt, 'pycircuit', 'circuit', 'csrc', unit)
        with open(p) as fh:
            text = fh.read()
        new, got = mutate.force_guards(text, _env_name(unit))
        if got:
            originals[unit] = text
            sites[unit] = got
            with open(p, 'w') as fh:
                fh.write(new)
    return sites, originals


def _restore(wt, originals):
    for unit, text in originals.items():
        with open(os.path.join(wt, 'pycircuit', 'circuit', 'csrc', unit), 'w') as fh:
            fh.write(text)


def _run(wt, store, r, fault):
    """One deck of the record `r`, the check `fault` (`(unit, line, asking)`, or None) armed: `(row, fired)`."""
    method, setting = r['group'].split('/')
    s = r['row'].get('s') or 1.0
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', PYTHONPATH=wt, PYCIRCUIT_HDL_CACHE_DIR=store)
    if fault:
        env[_env_name(fault[0])] = f'{fault[1]}:{fault[2]}'
    limit = int(max(120, 8 * s))
    cmd = ['timeout', str(limit), 'systemd-run', '--user', '--scope', '-q', '-p', 'MemoryMax=4G', '-p', 'MemorySwapMax=0',
           G.PY, os.path.join(HERE, 'run_deck.py'), r['deck'], method, setting]
    try:
        p = subprocess.run(cmd, cwd=wt, env=env, capture_output=True, text=True, timeout=limit + 30, check=False)
    except subprocess.TimeoutExpired:
        return {'failed': True, 'why': 'over time'}, True
    out = p.stdout.strip().rsplit('\n', 1)[-1] if p.stdout.strip() else ''
    fired = 'PCX_FAULT_FIRED' in p.stderr
    if not out.startswith('{'):
        return {'failed': True, 'why': p.stderr.strip()[-200:] or f'exit {p.returncode}'}, fired
    return json.loads(out), fired


def _asked_in(a):
    """`{(unit, line): [scope, ...]}`: the settings whose coverage record shows a line asked."""
    notes = os.path.join(a.ref, 'notes')
    got = collections.defaultdict(list)
    for scope in sorted(os.listdir(os.path.join(a.ref, 'scopes'))):
        if not scope.endswith('_c') or scope.startswith('size'):
            continue
        method, setting = scope[:-2].split('_', 1)
        if (method, setting) in settings.NO_TWIN:
            continue            # (its state follows the calls: a repeat of an attempt is not as if nothing happened, by design)
        lines, _arms = gcov_read.read(os.path.join(a.ref, 'scopes', scope), notes, set(UNITS))
        for (u, ln), c in lines.items():
            if c:
                got[(u, ln)].append(f'{method}/{setting}')
    return got


def one_check(unit, line, text, groups, ref_by, a, store):
    """A check on up to `DECKS_A_CHECK` decks of the first settings that ask it, at each asking."""
    runs, fired_decks, worst, why = 0, set(), 'not fired', None
    rank = {'not fired': 0, 'same': 1, 'counts': 2, 'differs': 3}
    for g in groups:
        for r in ref_by.get(g, ()):
            if len(fired_decks) >= DECKS_A_CHECK:
                break
            for asking in ASKINGS:
                row, fired = _run(a.wt, store, r, (unit, line, asking))
                runs += 1
                if not fired:
                    if asking == ASKINGS[0] and not row.get('failed'):
                        break           # (this deck never asks it)
                    if not row.get('failed'):
                        continue        # (asked fewer times than that: nothing fired, nothing to read)
                fired_decks.add(r['deck'])
                d = V.differing(row, r['row'])
                res = {k: v for k, v in d.items() if k in RESULTS or k == 'failed'}
                st = 'differs' if res else 'counts' if d else 'same'
                if rank[st] > rank[worst]:
                    worst, why = st, {'group': g, 'deck': r['deck'], 'asking': asking, 'what': str(res or d)[:300], 'how': row.get('why')}
        if len(fired_decks) >= DECKS_A_CHECK:
            break
    return {'unit': unit, 'line': line, 'check': text.strip(), 'status': worst, 'decks': sorted(fired_decks), 'runs': runs, 'what': why}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', required=True)
    ap.add_argument('--ref', required=True, help="a gate's directory: its matrix.jsonl, its coverage")
    ap.add_argument('--wt', required=True, help='a clean worktree at this commit')
    ap.add_argument('--jobs', type=int, default=8)
    ap.add_argument('--only', help='comma-separated units')
    a = ap.parse_args()
    a.out, a.ref, a.wt = (os.path.abspath(x) for x in (a.out, a.ref, a.wt))
    if subprocess.run(['git', 'status', '--short', '--', 'pycircuit'], cwd=a.wt, capture_output=True, text=True, check=False).stdout.strip():
        raise SystemExit('the worktree is not clean')
    ref = [r for r in G._rows(a.ref, 'matrix.jsonl') if r['path'] == 'c' and not r['row'].get('failed')]
    ref_by = collections.defaultdict(list)
    for r in sorted(ref, key=lambda r: r['row'].get('s') or 0):
        ref_by[r['group']].append(r)
    store = G.seed_store(os.path.join(a.out, 'store'))
    asked = _asked_in(a)
    log = G._Log(a.out, 'handback.jsonl')
    sites, originals = _instrument(a.wt)
    try:
        ## the control: nothing armed, every deck by default, on the instrumented build
        defaults = [r for g in ('gear2/default', 'radau/default', 'dc/default') for r in ref_by[g]]
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            bad = [(r['group'], r['deck'], V.differing(row, r['row']), row.get('why')) for r, (row, _f) in zip(defaults, ex.map(lambda r: _run(a.wt, store, r, None), defaults))
                   if V.differing(row, r['row'])]
        print(f'control, nothing armed: {len(defaults)} runs, {len(bad)} differing {str(bad)[:900]}', flush=True)
        if bad:
            raise SystemExit('the instrumented build is not the record')
        units = a.only.split(',') if a.only else sorted(sites)
        tasks = [(u, ln, t) for u in units for ln, t in sites.get(u, {}).items()]
        print(f'{len(tasks)} checks in {len(units)} units', flush=True)

        def do(task):
            u, ln, t = task
            rec = one_check(u, ln, t, asked.get((u, ln), []), ref_by, a, store)
            log.add(**rec)
            print(f"{u}:{ln} {rec['status']} {rec['decks']} ({rec['runs']} runs) {str(rec['what'])[:170] if rec['what'] else ''}", flush=True)
            return rec['status']
        with concurrent.futures.ThreadPoolExecutor(a.jobs) as ex:
            tally = collections.Counter(ex.map(do, tasks))
        print('tally:', dict(tally))
    finally:
        _restore(a.wt, originals)


if __name__ == '__main__':
    main()
