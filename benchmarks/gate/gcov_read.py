"""The coverage build's counts read through gcov (`_cbuild`'s `PYCIRCUIT_C_COVERAGE`): which C lines and taken branch arms a scope reached.

A scope is one directory of counts (`GCOV_PREFIX`).  Lines and arms are keyed `(file, line)` and `(file, line, i)` on `csrc/`'s own names and
lines (the build's `#line` markers).  ⚠ A LINE IS NOT A BRANCH: a one-line `if (x) return BAIL;` counts as reached when it never fires, so arms are what to
read; and at -O2 some arms are the compiler's own and cannot be taken, so a list of untaken arms is an upper bound, to be read."""
import collections
import glob
import json
import os
import subprocess


def _gcov(cwd, stem, csrc):
    out = subprocess.run(['gcov', '--json-format', '--stdout', '--branch-probabilities', stem], cwd=cwd, capture_output=True, text=True, check=False).stdout
    for doc in out.splitlines():
        if not doc.startswith('{'):
            continue
        for f in json.loads(doc)['files']:
            b = os.path.basename(f['file'])
            if b in csrc:
                yield b, f['lines']


def read(scope_dir, notes_dir, csrc):
    """`(lines, arms)`: Counters of the scope's counts, every line and arm its units hold (at 0 where not reached)."""
    lines, arms = collections.Counter(), collections.Counter()
    for g in glob.glob(os.path.join(scope_dir, '*.gcda')):
        stem = os.path.basename(g)[:-5]
        note = os.path.join(notes_dir, stem + '.gcno')
        if not os.path.exists(note):
            continue
        link = os.path.join(scope_dir, stem + '.gcno')
        if not os.path.exists(link):
            os.symlink(note, link)
        for b, ls in _gcov(scope_dir, stem, csrc):
            for ln in ls:
                lines[(b, ln['line_number'])] += ln['count']
                for i, br in enumerate(ln.get('branches', ())):
                    if not br.get('throw'):
                        arms[(b, ln['line_number'], i)] += br['count']
    return lines, arms


def universe(notes_dir, work_dir, csrc):
    """Every line and arm any built unit could reach: `(set of lines, set of arms)`, from the notes alone."""
    os.makedirs(work_dir, exist_ok=True)
    lines, arms = set(), set()
    for note in glob.glob(os.path.join(notes_dir, '*.gcno')):
        stem = os.path.basename(note)[:-5]
        link = os.path.join(work_dir, stem + '.gcno')
        if not os.path.exists(link):
            os.symlink(note, link)
        for b, ls in _gcov(work_dir, stem, csrc):
            for ln in ls:
                lines.add((b, ln['line_number']))
                arms.update((b, ln['line_number'], i) for i, br in enumerate(ln.get('branches', ())) if not br.get('throw'))
    return lines, arms


def taken(counter):
    return {k for k, c in counter.items() if c}
