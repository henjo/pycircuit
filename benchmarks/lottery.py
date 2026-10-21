"""A deck's own noise, from a lottery that is alive: what a change that only
nudges the step sequence moves its time by, per deck.

    python benchmarks/lottery.py                       # Radau's 20 decks
    python benchmarks/lottery.py mike2 arom            # by name
    python benchmarks/lottery.py --method gear2        # the default method
    python benchmarks/lottery.py --json table.json     # the table, for a line

Why: a re-baseline -- a change that moves the step sequence by a little --
moves an adaptive run's time and step count by the step controller's own
sensitivity, not by the change's cost, and a per-deck line ("no deck slower
than +3 %") then fails on noise.  Roadmap 159's line had an escape for that,
a lottery of GAMMA_REAL moved by +-2 ulp, and it was DEAD: since 157 the
transform path's estimate does not read GAMMA_REAL, mike2's digest was the
same at -2, 0 and +2 ulp, and the escape excused nothing.  And a live
last-bits lottery can still be flat (mike2 under reltol +-2 ulp: 614 steps
every time, every digest different) while a 1e-5 nudge -- the scale 159's
abort moved a retry -- spans 614 .. 651 steps, 0.575 .. 0.639 s.

So: per deck, the base run and arms with `reltol` scaled by (1 - 3e-5),
(1 - 1e-5), (1 + 1e-5), (1 + 3e-5) (`ARMS`); an arm counts only where it is
ALIVE -- its waveform's digest differs from the base's -- and a deck with no
live arm is reported dead (its lottery excuses nothing).  Each run in a
process of its own, BLAS on one thread and glibc's allocator pinned, as the
A/B lines are measured (the CPU is not pinned: the table describes their
conditions); rounds interleaved, the arm order reversed every other round;
best of the rounds per arm.  Under the benchmark lock (`_bench.lock`).

THE RULE A LINE CITES (`bound`): a deck fails only if slower than its
slowest live arm by more than `SLACK` (3 %); a dead deck against its base.
A summed criterion means something only beyond the sum's own spread.
Measured 2026-10-21 (a740330e, Radau): every deck but ab_integ alive, single
decks -32 % (gm1, its base its slowest) .. +18 % (arom), the 20-deck sum
-3.4 .. +2.9 %.  The table describes the tree it was measured on: measure
again where a change moves the step controller itself.
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

#: Radau's 20 decks (the stage methods' A/B set)
RADAU_DECKS = ('mike2', 'slowlatch', 'todd3', 'latch', 'nand', 'gm3', 'ring', 'ab_integ', 'g1310', 'gm1',
               'gm17', 'toronto', 'mux8', 'voter25', 'e1480', 'jge', 'rich3', 'arom', 'gm19', 'opampal')
#: the arms: `reltol` scaled by these
ARMS = (0.99997, 0.99999, 1.00001, 1.00003)
#: a deck's allowance over its slowest live arm (`bound`)
SLACK = 0.03
METHODS = ('radau', 'gear2')


def child(deck, method, factor):
    """One run: the deck's default transient under `method`, `reltol` scaled by
    `factor` -- its time, steps, rejections and waveform digest."""
    sys.path.insert(0, HERE)
    import _bench
    _bench.pin_threads()
    _bench.pin_malloc()
    sys.path.insert(0, ROOT)
    import time
    import warnings

    import numpy as np
    import spice_suite

    from pycircuit.circuit import circuit, spice_import
    from pycircuit.circuit import integrator as I
    warnings.simplefilter('ignore')
    circuit.default_toolkit = circuit.numeric
    imp = spice_import.import_netlist(spice_suite.benchdata.spice_data(spice_suite.CASES[deck][0]))
    tr, kw = imp.transient(**spice_suite.OVERRIDES.get(deck, {}))
    if method == 'radau':
        tr.par.integrator = I.RadauIIA3Integrator()
    if factor != 1.0:
        tr.par.reltol = float(tr.par.reltol) * factor
    t0 = time.perf_counter()
    res = tr.solve(**kw)
    s = time.perf_counter() - t0
    st = tr.statistics
    x = hashlib.sha1(np.asarray(res.x, float).tobytes()
                     + np.asarray(res.sweep_values, float).tobytes()).hexdigest()[:12]
    return {'deck': deck, 'method': method, 'f': factor, 's': s, 'steps': st.accepted_steps,
            'rej': st.rejected_steps, 'x': x}


def run_child(deck, method, factor):
    out = subprocess.run([sys.executable, os.path.abspath(__file__), '--child', deck, method,
                          repr(factor)], capture_output=True, text=True, cwd=ROOT, check=False)
    line = [ln for ln in out.stdout.splitlines() if ln.startswith('{')]
    if not line:
        return {'deck': deck, 'method': method, 'f': factor, 'error': out.stderr.strip()[-300:]}
    return json.loads(line[-1])


def table(rows):
    """{deck: record} from the runs `rows` (`child`'s records, several rounds):
    the base's best time and steps; the arms ALIVE (their digest not the
    base's) and their best times; the time and step span over the base and the
    live arms; `dead` where no arm is alive.  Raises where a deck's base or
    one arm gave two digests -- a run that is not deterministic has no lottery
    to read."""
    by = {}
    for r in rows:
        if 'error' in r:
            continue
        by.setdefault(r['deck'], {}).setdefault(r['f'], []).append(r)
    out = {}
    for deck, arms in by.items():
        for f, runs in arms.items():
            if len({r['x'] for r in runs}) != 1:
                raise ValueError(f'{deck}: the runs at {f} gave different waveforms -- not deterministic')
        base = arms.get(1.0)
        if not base:
            continue
        bx = base[0]['x']
        live = sorted(f for f in arms if f != 1.0 and arms[f][0]['x'] != bx)
        best = {f: min(r['s'] for r in arms[f]) for f in [1.0] + live}
        steps = {f: arms[f][0]['steps'] for f in [1.0] + live}
        out[deck] = {'base_s': best[1.0], 'min_s': min(best.values()), 'max_s': max(best.values()),
                     'base_steps': steps[1.0], 'min_steps': min(steps.values()),
                     'max_steps': max(steps.values()), 'live': len(live),
                     'dead': not live, 'arms': len([f for f in arms if f != 1.0]),
                     'rounds': len(base)}
    return out


def bound(tab, deck, slack=SLACK):
    """The time a deck may take under a line citing the table `tab`: its
    slowest live arm's plus `slack` -- its base's where it is dead."""
    rec = tab[deck]
    return (rec['base_s'] if rec['dead'] else rec['max_s']) * (1.0 + slack)


def show(tab):
    s0 = s1 = s2 = 0.0
    for deck, t in tab.items():
        b = t['base_s']
        s0 += b
        s1 += t['min_s']
        s2 += t['max_s']
        live = 'DEAD' if t['dead'] else f"live {t['live']}/{t['arms']}"
        print(f"{deck:10s} {live:9s} base {b:8.2f} s  {100 * (t['min_s'] / b - 1):+6.1f} .. "
              f"{100 * (t['max_s'] / b - 1):+6.1f} %  steps {t['min_steps']} .. {t['max_steps']} "
              f"(base {t['base_steps']})  bound {bound(tab, deck):8.2f} s")
    if s0:
        print(f"summed base {s0:.2f} s, the arms {100 * (s1 / s0 - 1):+.1f} .. {100 * (s2 / s0 - 1):+.1f} %")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('decks', nargs='*')
    ap.add_argument('--method', choices=METHODS, default='radau')
    ap.add_argument('--rounds', type=int, default=2)
    ap.add_argument('--json')
    ap.add_argument('--child', nargs=3, metavar=('DECK', 'METHOD', 'FACTOR'), help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.child:
        print(json.dumps(child(a.child[0], a.child[1], float(a.child[2]))))
        return 0
    sys.path.insert(0, HERE)
    import _bench
    decks = a.decks or list(RADAU_DECKS)
    rows = []
    with _bench.lock(what='lottery'):
        cond = _bench.stamp()
        for rnd in range(a.rounds):
            for deck in decks:
                order = (1.0,) + ARMS
                for f in (order if rnd % 2 == 0 else order[::-1]):
                    r = run_child(deck, a.method, f)
                    r['round'] = rnd + 1
                    rows.append(r)
                    if 'error' in r:
                        print(f"{deck} at {f}: FAILED {r['error'][-120:]}", file=sys.stderr)
    tab = table(rows)
    show(tab)
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump({'conditions': cond, 'method': a.method, 'arms': ARMS, 'slack': SLACK,
                       'table': tab, 'runs': rows}, fh, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
