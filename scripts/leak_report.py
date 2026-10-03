#!/usr/bin/env python
"""Summarise the leak detector's reports (2026-10-03).

    python scripts/leak_report.py <dir with leaks-*.jsonl> [--by module|change|test]

Groups every report line -- per test, per module boundary, at collection --
by the change reported (default), by the test module it came from, or lists
every test.  Changes are normalised (object addresses and long values cut)
so the same leak from several tests counts once per kind.
"""
import collections
import json
import os
import re
import sys


def load(d):
    out = []
    for name in sorted(os.listdir(d)):
        if name.startswith('leaks-') and name.endswith('.jsonl'):
            with open(os.path.join(d, name)) as f:
                out += [json.loads(ln) for ln in f if ln.strip()]
    return out


def norm(change):
    c = re.sub(r'@[0-9a-f]{6,}', '@…', change)
    c = re.sub(r"'[0-9a-f]{20,}'", "'<key>'", c)
    c = re.sub(r'\(first seen [^)]*\)', '', c)
    return c[:160]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        print(__doc__)
        return 2
    by = 'change'
    if '--by' in argv:
        by = argv[argv.index('--by') + 1]
    recs = load(argv[0])
    print(f'{len(recs)} reports')
    if by == 'test':
        for r in recs:
            print(f"[{r['kind']}] {r['where']}")
            for c in r['changes']:
                print(f'    {c[:200]}')
        return 0
    counts = collections.Counter()
    where = collections.defaultdict(set)
    for r in recs:
        mod = r['where'].split('::', 1)[0]
        for c in r['changes']:
            key = norm(c) if by == 'change' else mod
            counts[key] += 1
            where[key].add(mod if by == 'change' else norm(c))
    for key, n in counts.most_common():
        print(f'{n:5d}  {key}')
        for w in sorted(where[key])[:6]:
            print(f'         {w}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
