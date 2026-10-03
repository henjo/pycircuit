"""Speed round 4's harness: what one gear step costs around the device
evaluations (2026-10-02, the per-step solver machinery plan).

Five cases: the 7-node PSP stage (one device; the PSS target's circuit),
a 20-MosLevel1 chain, a 20-Gummel-Poon chain and a 20-PSP chain (45-61
unknowns), each a DC
operating point then 100 fixed gear steps, and the PSP stage's PSS (gear,
40 points, reltol 1e-8).  A case is timed over its `solve` and reported
per step; the solution's bytes and the run's statistics (every slot but
the seconds) are printed so two trees can be checked against each other
as well as timed -- bit-identity is the contract of every stage.

    python benchmarks/step_machinery.py                  # this tree, every case
    python benchmarks/step_machinery.py stage psp        # by name
    python benchmarks/step_machinery.py --compare DIR    # the tree at DIR (the
                                                         # parent) against this one:
                                                         # subprocesses, interleaved
    python benchmarks/step_machinery.py --tree stage     # the in-run timer tree
                                                         # of one case (per step)

`--compare` runs each tree in its own interpreter from its own directory,
alternating parent / child per round (default 5), and reports the MIN and
the MEDIAN per-step time of each, the change of each, and whether the
bytes and the statistics agree.  Read both numbers: the minimum is the
quiet-box figure, the median says whether the box was quiet (`uptime`
first; another session's benchmark here moves the median first).

Round-4 baseline (4af78ddd, this box, min of 5): stage 0.71 ms/step, 20
MosLevel1 1.40 ms/step, 20 PSP 3.11 ms/step, the stage PSS 0.266 s.  The
in-run tree (`--tree`) attributes a step to its pieces with inclusive
`perf_counter` timers; the timers and the cold caches between the PSP
kernel's calls inflate every Python piece 1.5-3x over its standalone
time, so attribute a saving to the piece it was predicted from, never
compare a tree's number with a standalone one.  The hdl elements' i/q/G/C
are not wrapped on the instance (that would be an instance shadow, which
sends them to the per-element path): their batch is timed as one piece
(`Batch.run`, speed round 6).
"""
import functools
import hashlib
import json
import os
import statistics
import subprocess
import sys
import time
import warnings

warnings.simplefilter('ignore')
import numpy as np

from pycircuit.circuit import PSS, circuit, compact
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, R, SubCircuit, VSin, gnd
from pycircuit.circuit.transient import Transient

circuit.default_toolkit = circuit.numeric

STEPS = 100


def chain(make, ndev=20, vdd=1.8, vg=0.9):
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c.add_node('g0')
    c['vg'] = VSin('g0', gnd, v=vg, va=2e-2, freq=1e6)
    for k in range(ndev):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'M{k}'] = make(f'd{k}', f'g{k}' if k == 0 else f'd{k-1}')
    return c


def gp_chain(ndev=20):
    """Common-emitter stages in a chain: each base through 10 k from the
    previous collector (the limiter's other library case: single probes,
    no parameter reading the solution)."""
    c = SubCircuit()
    c.add_node('vcc')
    c['vcc'] = VS('vcc', gnd, v=3.0)
    c.add_node('in')
    c['vin'] = VSin('in', gnd, v=0.75, va=2e-2, freq=1e6)
    prev = 'in'
    for k in range(ndev):
        c.add_node(f'b{k}')
        c.add_node(f'c{k}')
        c[f'rb{k}'] = R(prev, f'b{k}', r=1e4)
        c[f'rc{k}'] = R('vcc', f'c{k}', r=1e3)
        c[f'Q{k}'] = eh.GummelPoonNpnHdl(f'c{k}', f'b{k}', gnd)
        prev = f'c{k}'
    return c


def stage():
    c = SubCircuit()
    for n in ('g', 'd', 'vdd'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=1.2)
    c['vg'] = VSin('g', gnd, v=0.7, va=2e-2, freq=1e6)
    c['rl'] = R('vdd', 'd', r=5e3)
    c['M'] = compact.PspMosLongChannel('d', 'g', gnd, gnd, fnt=1.0)
    return c


BUILD = {
    'stage': stage,
    'mos1': lambda: chain(lambda d, g: eh.MosLevel1Hdl(d, g, gnd, gnd)),
    'psp': lambda: chain(lambda d, g: compact.PspMosLongChannel(
        d, g, gnd, gnd, fnt=1.0), vdd=1.2, vg=0.7),
    'gp': gp_chain,
}
CASES = ('stage', 'mos1', 'gp', 'psp', 'pss')


def _stats(tr):
    s = getattr(tr, 'statistics', None)
    d = getattr(s, '__dict__', None)
    if d is None:
        d = {k: getattr(s, k) for k in getattr(s, '__slots__', ())}
    return {k: v for k, v in sorted(d.items()) if 'seconds' not in k}


def run_case(name):
    """One run: `(seconds, per_step_us or None, sha of x, stats)`."""
    if name == 'pss':
        c = stage()
        p = PSS(c, method='gear', reltol=1e-8)
        t0 = time.perf_counter()
        p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
        dt = time.perf_counter() - t0
        x = np.asarray(p.waveform[1], float)
        return dt, None, hashlib.sha256(x.tobytes()).hexdigest()[:12], None
    c = BUILD[name]()
    tr = Transient(c, toolkit=circuit.numeric)
    t0 = time.perf_counter()
    res = tr.solve(tend=STEPS * 2e-8, timestep=2e-8, fixed_timestep=True)
    dt = time.perf_counter() - t0
    x = np.asarray(res.x, float)
    return (dt, 1e6 * dt / STEPS, hashlib.sha256(x.tobytes()).hexdigest()[:12],
            _stats(tr))


def measure(cases, rounds):
    """Every case warmed once, then `rounds` interleaved runs; a dict per
    case with the times, the digest and the statistics."""
    out = {k: {'times': [], 'per_step': []} for k in cases}
    for k in cases:
        run_case(k)
    for _ in range(rounds):
        for k in cases:
            dt, us, sha, st = run_case(k)
            out[k]['times'].append(dt)
            out[k]['per_step'].append(us)
            out[k]['sha'] = sha
            out[k]['stats'] = st
    return out


def _fmt(k, r):
    if k == 'pss':
        return f'{min(r["times"]):.3f} s (median {statistics.median(r["times"]):.3f})'
    return (f'{min(r["per_step"]):7.1f} us/step (median '
            f'{statistics.median(r["per_step"]):7.1f})')


def compare(parent, cases, rounds):
    """The tree at `parent` against this one, each in its own interpreter
    from its own directory, alternating per round."""
    child = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    trees = {'parent': os.path.abspath(parent), 'child': child}
    got = {lab: {k: {'times': [], 'per_step': []} for k in cases} for lab in trees}
    for _ in range(rounds):
        for lab, tree in trees.items():
            ## THIS script in both trees (the parent may predate it): the
            ## package comes from PYTHONPATH, which precedes site-packages
            env = dict(os.environ, PYTHONPATH=tree)
            cmd = [sys.executable, os.path.abspath(__file__),
                   '--json', '--rounds', '1', *cases]
            res = subprocess.run(cmd, cwd=tree, env=env, capture_output=True, text=True,
                                 check=True)
            for line in res.stdout.splitlines():
                if not line.startswith('{'):
                    continue
                rec = json.loads(line)
                g = got[lab][rec['case']]
                g['times'] += rec['times']
                g['per_step'] += rec['per_step']
                g['sha'] = rec['sha']
                g['stats'] = rec['stats']
    for k in cases:
        p, c = got['parent'][k], got['child'][k]
        key = 'times' if k == 'pss' else 'per_step'
        dmin = 100.0 * (min(c[key]) / min(p[key]) - 1.0)
        dmed = 100.0 * (statistics.median(c[key]) / statistics.median(p[key]) - 1.0)
        same = ('bytes ' + ('SAME' if p['sha'] == c['sha'] else 'DIFFER')
                + ', stats ' + ('SAME' if p['stats'] == c['stats'] else 'DIFFER'))
        print(f'{k:5s} parent {_fmt(k, p)} | child {_fmt(k, c)} | '
              f'min {dmin:+.1f} % median {dmed:+.1f} % | {same}', flush=True)


## -- the in-run timer tree -----------------------------------------------------

def tree(name):
    """Inclusive `perf_counter` timers on the step's pieces over one run of
    `name`, printed per step, biggest first."""
    from pycircuit.circuit import _tran_newton, nrsolver
    from pycircuit.circuit import transient as TR
    acc = {}

    def timed(f, lab):
        acc.setdefault(lab, [0, 0])

        @functools.wraps(f)
        def w(*a, **k):
            t0 = time.perf_counter_ns()
            try:
                return f(*a, **k)
            finally:
                acc[lab][0] += time.perf_counter_ns() - t0
                acc[lab][1] += 1
        return w

    def wrap(obj, attr, lab=None):
        setattr(obj, attr, timed(getattr(obj, attr), lab or attr))

    orig_rr = _tran_newton.refnode_removed
    _tran_newton.refnode_removed = lambda f, i, tk: timed(
        orig_rr(f, i, tk), 'eval_FJ closure (reinsert + residual + reduce)')
    for cls in (nrsolver.StandardNewton, nrsolver.ChordNewton):
        wrap(cls, 'solve_system', f'{cls.__name__}.solve_system')
    for nm in ('where', 'take', 'judge', 'event', 'accept'):
        wrap(TR._SteppingLoop, nm, f'loop.{nm}')
    c = BUILD[name]()
    tr = Transient(c, toolkit=circuit.numeric)
    tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    tr = Transient(c, toolkit=circuit.numeric)
    for nm in ('solve_timestep', '_newton', '_residual_and_jacobian', '_predict_state',
               '_push_history', '_branch_after_solve', '_branch_screen',
               '_newton_abstol_vector', '_newton_xtol_vector', '_source_at', 'get_diff',
               '_C_lookup', '_companion_at', '_C_at_state'):
        wrap(tr, nm)
    orig_nl = tr._newton_limiter

    def newton_limiter():
        f = orig_nl()
        return None if f is None else timed(
            f, 'limiter_func (reinsert x2 + cir.limit + remove)')
    tr._newton_limiter = newton_limiter
    for nm in ('i', 'q', 'G', 'C', 'u', 'limit', 'accept_step', 'next_event'):
        wrap(c, nm, f'cir.{nm}')
    ## (an instance wrap IS an instance shadow, and a shadow sends an hdl
    ## element to the per-element path -- so the hdl elements' i/q/G/C are
    ## not wrapped, and their batch (`_hdl_batch`, one C call per class per
    ## pass) is timed as one piece; `u` is not batched and stays counted)
    from pycircuit.circuit import _hdl_batch
    for el in c.elements.values():
        hdl_el = getattr(type(el), '_hdl_info', None) is not None
        for nm in ('i', 'q', 'G', 'C', 'u'):
            if hasattr(el, nm) and not (hdl_el and nm != 'u'):
                wrap(el, nm, f'element.{nm} calls')
    wrap(_hdl_batch.Batch, 'run', 'Batch.run (one C call per class per pass)')
    ls = tr._get_linearsolver()
    wrap(ls, 'solve', 'linsolver.solve')
    wrap(ls, 'factor', 'linsolver.factor')
    t0 = time.perf_counter_ns()
    tr.solve(tend=STEPS * 2e-8, timestep=2e-8, fixed_timestep=True)
    total = (time.perf_counter_ns() - t0) / STEPS
    print(f'== {name}: n={c.n} {total / 1e3:7.1f} us/step (solve wall / {STEPS})')
    for lab, (ns, cnt) in sorted(acc.items(), key=lambda kv: -kv[1][0]):
        if cnt:
            print(f'   {lab:52s} {ns / STEPS / 1e3:7.1f} us/step {100 * ns / STEPS / total:5.1f} %'
                  f'  ({cnt / STEPS:6.2f} calls/step, {ns / cnt / 1e3:7.2f} us each)')


def main(argv):
    rounds = 5
    as_json = False
    cmp_dir = None
    tree_case = None
    cases = []
    it = iter(argv)
    for a in it:
        if a == '--rounds':
            rounds = int(next(it))
        elif a == '--json':
            as_json = True
        elif a == '--compare':
            cmp_dir = next(it)
        elif a == '--tree':
            tree_case = next(it)
        else:
            cases.append(a)
    cases = cases or list(CASES)
    if tree_case:
        tree(tree_case)
        return
    if cmp_dir:
        compare(cmp_dir, cases, rounds)
        return
    got = measure(cases, rounds)
    for k in cases:
        r = got[k]
        if as_json:
            print(json.dumps({'case': k, **r}), flush=True)
        else:
            print(f'{k:5s} {_fmt(k, r)}  (digest {r["sha"]})', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
