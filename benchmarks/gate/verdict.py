"""What a gate concludes, as pure functions of what it measured: nothing here runs a circuit, so every branch is tested in the suite
(tests/test_gate_verdict.py) whatever the machine is doing.

A VERDICT HAS THREE OUTCOMES.  'DIFFERS': a comparison differs, or a run the gate asked for is missing or failed (a deck missing from a
census is a failure until shown otherwise).  'NOT REACHED': every comparison agrees, and something the gate was asked about was not
compared: nothing was compared at all, a named counter stayed at zero, or a changed C line or arm was not taken on a compared run.
'SAME': they agree and all of it was reached.  "0 differing" with the code unreached reads as a pass otherwise: it did, twice on 2026-11-07.
"""
import collections
import re

#: what a C path and its Python must agree on, and two trees on the same path
KEYS = ('x', 'err', 'steps', 'rej', 'newton', 'fb', 'klu')
#: an arm or line's class: taken on a real deck by default, only under a forcing, only at size, only by the suite's tests, by nothing
CLASSES = ('decks', 'forced', 'size', 'fixture', 'none')
COMPARED = ('decks', 'forced', 'size')
#: the fewest unknowns of a benchmark that should move before it says much about a long run
QUIET_SHARE = 0.01


def differing(a, b, keys=KEYS):
    """The keys on which rows `a` and `b` differ, `{key: (a's, b's)}`.  A failed row differs from everything, another failed row too."""
    if a is None or b is None or a.get('failed') or b.get('failed'):
        return {'failed': (bool(a is None or a.get('failed')), bool(b is None or b.get('failed')))}
    return {k: (a.get(k), b.get(k)) for k in keys if a.get(k) != b.get(k)}


def compare(rows, left='c', right='py', keys=KEYS):
    """Rows `{'group': hashable, 'deck': name, 'path': p, 'row': dict}` paired by group and deck: `{group: {'runs': n, 'differ': [(deck,
    what)]}}`.  A deck with one side only differs."""
    by = collections.defaultdict(dict)
    for r in rows:
        by[(r['group'], r['deck'])][r['path']] = r['row']
    out = collections.defaultdict(lambda: {'runs': 0, 'differ': []})
    for (g, d), v in sorted(by.items(), key=lambda kv: str(kv[0])):
        if left not in v and right not in v:
            continue
        out[g]['runs'] += 1
        what = differing(v.get(left), v.get(right), keys)
        if what:
            out[g]['differ'].append((d, what))
    return dict(out)


def with_side(rows, side):
    """The rows of the (group, deck) pairs that hold a row on path `side`: only those are compared with it."""
    have = {(r['group'], r['deck']) for r in rows if r['path'] == side}
    return [r for r in rows if (r['group'], r['deck']) in have]


def reach(rows, keys, path='c'):
    """`{key: total count}` over the compared `path` rows, every asked key present: 0 where nothing counted it."""
    got = {k: 0 for k in keys}
    for r in rows:
        if r['path'] == path and not r['row'].get('failed'):
            c = r['row'].get('counts', {})
            for k in keys:
                got[k] += c.get(k, 0)
    return got


def count_moves(a, b, allow=None):
    """The path counters that differ between two runs, `{key: (a's, b's)}`: results the same with a counter moved is a fast path that stopped serving, or one
    that started.  Not a `once:` key (a fact decided once a process, landing on whichever run came first), nor one matching `allow`, a stage's own moves."""
    keys = set(a) | set(b)
    return {k: (a.get(k, 0), b.get(k, 0)) for k in sorted(keys)
            if a.get(k, 0) != b.get(k, 0) and not k.startswith('once:') and not (allow and allow.search(k))}


def compare_counts(rows, left='new', right='base', allow=None):
    """`compare`'s shape on the path counters, the rows that failed left to `compare`."""
    by = collections.defaultdict(dict)
    for r in rows:
        by[(r['group'], r['deck'])][r['path']] = r['row']
    out = collections.defaultdict(lambda: {'runs': 0, 'differ': []})
    for (g, d), v in sorted(by.items(), key=lambda kv: str(kv[0])):
        a, b = v.get(left), v.get(right)
        if a is None or b is None or a.get('failed') or b.get('failed'):
            continue
        out[g]['runs'] += 1
        moved = count_moves(a.get('counts', {}), b.get('counts', {}), allow)
        if moved:
            out[g]['differ'].append((d, moved))
    return dict(out)


def arm_class(key, decks, forced, size, fixture):
    """The class of a line or arm, by the sets that took it (the first holding it, in `CLASSES`' order)."""
    for name, taken in zip(CLASSES, (decks, forced, size, fixture)):
        if key in taken:
            return name
    return 'none'


def changed_lines(diff):
    """A unified diff with no context (`git diff -U0`): `{file base name: {line numbers on the new side}}`."""
    out = collections.defaultdict(set)
    name = None
    for ln in diff.splitlines():
        if ln.startswith('+++ '):
            name = None if ln[4:].strip() == '/dev/null' else ln[4:].strip().rsplit('/', 1)[-1]
        elif ln.startswith('@@') and name:
            m = re.match(r'@@ -\S+ \+(\d+)(?:,(\d+))? @@', ln)
            start, count = int(m.group(1)), int(m.group(2) or 1)
            out[name].update(range(start, start + count))
    return dict(out)


def changed_report(changed, universe_lines, line_class, arms_of, arm_class_of):
    """Each changed line a compiler made code for, classed.  `universe_lines` holds `(file, line)`; `line_class(key)` and `arm_class_of(key)`
    class a line and an arm; `arms_of(key)` lists a line's arms `(file, line, i)`.  `{'lines': [(key, class)], 'arms': [(arm, class)]}`:
    non-executable lines (a comment, a declaration) are not in it."""
    lines, arms = [], []
    for f in sorted(changed):
        for ln in sorted(changed[f]):
            key = (f, ln)
            if key in universe_lines:
                lines.append((key, line_class(key)))
                arms.extend((a, arm_class_of(a)) for a in arms_of(key))
    return {'lines': lines, 'arms': arms}


def verdict(comparisons, reached=None, changed=None):
    """`(outcome, reasons)`.  `comparisons`: `{name: {'runs': n, 'differ': [...]}}` (`compare`'s, any number merged); `reached`: `reach`'s;
    `changed`: `changed_report`'s."""
    reasons = []
    differ = [(n, d) for n, c in comparisons.items() for d in c['differ']]
    if differ:
        return 'DIFFERS', [f'{n}: {d[0]} {d[1]}' for n, d in differ]
    if not sum(c['runs'] for c in comparisons.values()):
        reasons.append('nothing was compared')
    reasons += [f'{k} reached 0 times on the compared runs' for k, v in (reached or {}).items() if not v]
    for key, cls in (changed or {}).get('lines', ()):
        if cls not in COMPARED:
            reasons.append(f'changed line {key[0]}:{key[1]} not compared ({cls})')
    for key, cls in (changed or {}).get('arms', ()):
        if cls not in COMPARED:
            reasons.append(f'changed arm {key[0]}:{key[1]}#{key[2]} not compared ({cls})')
    return ('NOT REACHED', reasons) if reasons else ('SAME', [])


def quiet(rows):
    """The benchmarks of `rows` (`{'name', 'n', 'swing'}`) where under `QUIET_SHARE` of the unknowns move over 100 mV."""
    return [r['name'] for r in rows if r.get('swing') is not None and r['swing'] < QUIET_SHARE * r['n']]


def line_met(ratios, op, bound):
    """A timing line on the median of interleaved pairs' ratios: `(median, met)`.  No pairs is not met."""
    if not ratios:
        return None, False
    s = sorted(ratios)
    m = s[len(s) // 2] if len(s) % 2 else (s[len(s) // 2 - 1] + s[len(s) // 2]) / 2
    return m, (m >= bound if op == '>=' else m <= bound)
