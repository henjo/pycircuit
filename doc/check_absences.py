#!/usr/bin/env python3
"""Assert that each literature gap the PSS roadmap records as ABSENT is still absent -- so an acquisition
makes this FAIL and names the file, instead of a "recorded so nobody searches again" note rotting.

Two of this repository's absence markers rotted within a day this week (Lamour-Marz-Tischendorf, the
Sickenberger Part I report), and the docs session found four of its ten stale in one sweep.  A marker that
forbids the retry cannot be maintained by remembering to retry; a check that fails when the gap fills can.

Every ABSENT pattern is paired with PRESENT controls searched by the SAME mechanism, so a search that fails
open (the docs session's first Lamour grep used BRE alternation inside an ERE and matched a literal) is
caught: a control that does not hit means the search is broken, not that the corpus is empty.

Usage:  python doc/check_absences.py [corpus-root]      (default ~/docs; reads .corpus/manifest.tsv and
        walks the tree)   exit 0 = every absence still absent and every control present; 1 otherwise.
"""
import os
import re
import sys

ROOT = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/docs')

## (label, regex over the file path, where the roadmap records it)
ABSENT = [
    ('Lamour-Marz-Tischendorf Ch. 5 (GLM sections 5.2.3/5.3.3/5.5.3)', r'Lamour-2013.*\(Ch5\b|Lamour-2013.*Ch ?5\b',
     'pss_roadmap_260902.md, absence markers re-audited 2026-09-09'),
    ('Lamour-Marz-Tischendorf Ch. 10 (index sections 10.2.2.2-3)', r'Lamour-2013.*Ch ?10\b',
     'pss_roadmap_260902.md, absence markers re-audited 2026-09-09'),
    ("Hegazi & Abidi, the Leeson chapter", r'Hegazi|Abidi',
     'pss_roadmap_260902.md, absence markers re-audited 2026-09-09'),
    ('SSP for GLMs (Spijker; Ferracina & Spijker; Higueras)', r'Spijker|Ferracina|Higueras',
     'pss_roadmap_260902.md, absence markers re-audited 2026-09-09'),
]
## controls: known-present files the same search must hit (regexes chosen so they cannot match a sibling:
## "Part I " with the space and a following word, not "Part I" which is a prefix of "Part II")
PRESENT = [
    ('Sickenberger 2007 Part I (ODEs and DAEs)', r'Sickenberger-2007.*Part I ODEs'),
    ('Hairer, Lubich & Roche 1989 (RK for DAEs)', r'Hairer-1989-Runge-Kutta methods for differential-algebraic'),
    ('Lamour-Marz-Tischendorf Ch. 3', r'Lamour-2013.*Ch3'),
    ('Wright 2002 thesis', r'Wright-2002-General Linear Methods with Inherent'),
    ('Voigtmann thesis', r'Voigtmann-2006-General Linear Methods'),
]


def paths():
    seen = set()
    man = os.path.join(ROOT, '.corpus', 'manifest.tsv')
    if os.path.exists(man):
        for line in open(man, encoding='utf-8', errors='replace'):
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 5:
                seen.add(parts[4].lstrip('./'))
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in dirnames if not d.startswith('.')]
        for f in filenames:
            seen.add(os.path.relpath(os.path.join(dirpath, f), ROOT))
    return sorted(seen)


def main():
    allp = paths()
    if not allp:
        print('no corpus at %s -- nothing checked' % ROOT)
        return 2
    ok = True
    for label, rx in PRESENT:
        hits = [p for p in allp if re.search(rx, p, re.I)]
        if hits:
            print('control present : %-48s -> %s' % (label, hits[0]))
        else:
            ok = False
            print('CONTROL MISSING : %-48s (pattern %r hit nothing: the SEARCH is broken)' % (label, rx))
    for label, rx, where in ABSENT:
        hits = [p for p in allp if re.search(rx, p, re.I)]
        if hits:
            ok = False
            print('FILLED          : %-48s -> %s   (update %s)' % (label, hits[0], where))
        else:
            print('still absent    : %-48s' % label)
    print('%d files scanned; %s' % (len(allp), 'all markers current' if ok else 'ACTION NEEDED'))
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
