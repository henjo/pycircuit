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
                                                         # of one case (RANKS pieces)
    python benchmarks/step_machinery.py --sample mos1    # a sampled profile (py-spy)
    python benchmarks/step_machinery.py --check          # the per-commit check against
                                                         # the parent tree (verdicts)

`--check` (2026-10-04, testing for development, stage 3): every case of
`CASES` -- the four fixed-step gear steps, the stage's PSS, and the analysis
cases (`ANALYSES`: adaptive gear and Radau transients on a diode ladder, a
van der Pol PSS with its factored period, pnoise, the PPV) -- against the
tree `scripts/parent_tree.sh` keeps at HEAD (`--parent DIR` for another,
`--child DIR` for another child), with a verdict each (`verdict`): SLOWER,
FASTER, UNSURE (more rounds while `--budget`, 300 s, lasts) or "no change".
Exit 3 on SLOWER, 2 when bytes or statistics differ.  The fast-path counts
of every timed call are compared and reported.  `--plant P` slows the
child by P %: the check's own proof (+3 % gives SLOWER on every case, +1 %
on none).

`--compare` (2026-10-03, `_bench.py`): under the benchmark lock, each child
pinned to one P-core with one BLAS thread, the order alternating per round,
rounds with the core's sibling busy discarded, and the result the median of
PAIRED per-round ratios with a bootstrap 95 % interval; every run's bytes
and statistics are checked and any difference exits 2.  An A/A run
(`--compare .`) shows the noise floor: its interval covers 0.

Each side of a round is the minimum of five warm runs in one process; the
line under each case lists the per-round ratios, so an outlier round is
visible rather than averaged in.  The printed min and median per side are
over the rounds: the minimum is the quiet-box figure, the median says
whether the box was quiet.

Round-4 baseline (4af78ddd, this box, min of 5): stage 0.71 ms/step, 20
MosLevel1 1.40 ms/step, 20 PSP 3.11 ms/step, the stage PSS 0.266 s.  The
in-run tree (`--tree`) attributes a step to its pieces with inclusive
`perf_counter` timers; the timers and the cold caches between the PSP
kernel's calls inflate every Python piece 1.5-3x over its standalone
time, so attribute a saving to the piece it was predicted from, never
compare a tree's number with a standalone one.  The hdl elements are not
wrapped on the instance (that would be an instance shadow, which sends
them to the per-element path), and no element's `u` is (a shadow sends a
zero `u` to a call): their batch is timed as one piece (`Batch.run`,
speed round 6), `cir.u` as one pass.
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
## THE CONDITIONS FIRST (robust timing, 2026-10-03; `_bench.py`): one BLAS
## thread BEFORE numpy opens its pool (the `pss` case ran on 24 threads:
## PSS never enters `Transient.solve`'s single-thread scope), and, in a
## `--compare` child, one P-core
import _bench  # benchmarks/, this script's directory

_bench.pin_threads()
if os.environ.get('PYCIRCUIT_BENCH_PIN'):
    _bench.pin_cpu()
import numpy as np

from pycircuit.circuit import PSS, circuit, compact
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import (
    IS,
    VS,
    BSource,
    Diode,
    L,
    R,
    SubCircuit,
    VSin,
    gnd,
)
from pycircuit.circuit.elements import C as Cap
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
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


## -- the analysis cases (2026-10-04, testing for development, stage 3) -----------
##
## The performance check (`--check`) covers more than the fixed-step gear
## step: adaptive transients (gear, Radau), a PSS with its factored period,
## pnoise and the PPV.  Their builders are COPIED here, not imported from
## the tests or `speed_analyses.py`: under `--compare` the parent runs this
## script against its own package, and a fixture it imported from its own
## tests could be another circuit.  The pnoise and PPV cases time the
## analysis on a PSS solved once per process (untimed): the solve is the
## `vdp_pss` case's work.

def _ladder(n=20):
    """An RC ladder driven by a sine, a diode to ground at its end
    (`speed_analyses.ladder`)."""
    c = SubCircuit()
    c['vs'] = VSin('n0', gnd, va=2.0, freq=1e3)
    for k in range(n):
        c[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
        c[f'C{k}'] = Cap(f'n{k + 1}', gnd, c=1e-8)
    c['D'] = Diode(f'n{n}', gnd)
    return c


def _vdp(mu=1.0, psd=1e-6):
    """van der Pol with a white current source at its core
    (`tests/_shooting_fixtures._vdp_with_noise`) and its PSS seed."""
    c = SubCircuit()
    c.add_node('v')
    c['C'] = Cap('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: mu * (u - u ** 3 / 3.0))
    c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    x0 = np.zeros(c.n - 1)
    x0[0] = 2.0
    return c, x0


def _sha(a):
    return hashlib.sha256(np.ascontiguousarray(np.asarray(a, float)).tobytes()).hexdigest()[:12]


def _paths_now():
    try:
        from pycircuit.circuit import _paths
    except ImportError:            # a tree from before the fast-path counts
        return None
    return _paths.snapshot()


def _paths_since(before):
    if before is None:
        return None
    from pycircuit.circuit import _paths
    return {k: v for k, v in sorted(_paths.since(before).items()) if not k.startswith('once:')}


def _timed(fn):
    """`(seconds, fn(), the fast-path counts it added)`.  With
    `PYCIRCUIT_BENCH_PLANT=p` (`--plant`, the check's own proof) the call is
    followed by a busy wait of `p` times its time, inside the timing: a
    calibrated slowdown of exactly `p`."""
    p0 = _paths_now()
    t0 = time.perf_counter()
    out = fn()
    dt = time.perf_counter() - t0
    paths = _paths_since(p0)
    plant = float(os.environ.get('PYCIRCUIT_BENCH_PLANT', '0') or 0)
    if plant:
        end = time.perf_counter() + plant * dt
        while time.perf_counter() < end:
            pass
        dt = time.perf_counter() - t0
    return dt, out, paths


_SETUP = {}


def _case_ladder(method):
    integ = Gear2Integrator() if method == 'gear' else RadauIIA3Integrator()
    tr = Transient(_ladder(), integrator=integ, reltol=1e-5)
    dt, res, paths = _timed(lambda: tr.solve(tend=3e-3, timestep=1e-5))
    return dt, None, _sha(res.x), _stats(tr), paths


def _case_vdp_pss():
    cir, x0 = _vdp()
    p = PSS(cir, method='gear', reltol=1e-12)

    def go():
        p.solve(period=6.6634, timestep=6.6634 / 200, x0=x0, maxiterations=60)
        p.factored_period()
    dt, _out, paths = _timed(go)
    return dt, None, _sha(p.waveform[1]), None, paths


def _case_pnoise():
    from pycircuit.circuit.shooting import PAC
    if 'pnoise' not in _SETUP:
        cir, x0 = _vdp()
        p = PSS(cir, method='gear', reltol=1e-12)
        p.solve(period=6.6634, timestep=6.6634 / 200, x0=x0, maxiterations=60)
        _SETUP['pnoise'] = (cir, p, np.logspace(-3, -0.5, 20) / float(p.period))
    cir, p, fs = _SETUP['pnoise']
    pac = PAC(cir, toolkit=circuit.numeric)
    dt, out, paths = _timed(lambda: pac.pnoise(p, fs, 0, maxsidebands=8))
    return dt, None, _sha(out[0]), None, paths


def _case_ppv():
    if 'ppv' not in _SETUP:
        cir, x0 = _vdp()
        p = PSS(cir, method='radau', reltol=1e-10)
        p.solve(period=6.6634, timestep=6.6634 / 200, x0=x0, maxiterations=60)
        _SETUP['ppv'] = p
    p = _SETUP['ppv']
    dt, out, paths = _timed(p.ppv)
    v = out[0] if isinstance(out, tuple) else out
    return dt, None, _sha(getattr(v, 'y', v)), None, paths


ANALYSES = {
    'ladder_gear': lambda: _case_ladder('gear'),
    'ladder_radau': lambda: _case_ladder('radau'),
    'vdp_pss': _case_vdp_pss,
    'pnoise': _case_pnoise,
    'ppv': _case_ppv,
}
CASES = CASES + tuple(ANALYSES)


def _stats(tr):
    s = getattr(tr, 'statistics', None)
    d = getattr(s, '__dict__', None)
    if d is None:
        d = {k: getattr(s, k) for k in getattr(s, '__slots__', ())}
    return {k: v for k, v in sorted(d.items()) if 'seconds' not in k}


def run_case(name):
    """One run: `(seconds, per_step_us or None, sha of the result, stats,
    the fast-path counts of the timed call or None)`."""
    if name in ANALYSES:
        return ANALYSES[name]()
    if name == 'pss':
        c = stage()
        p = PSS(c, method='gear', reltol=1e-8)
        dt, _out, paths = _timed(lambda: p.solve(period=1e-6, timestep=1e-6 / 40,
                                                 maxiterations=60))
        return dt, None, _sha(p.waveform[1]), None, paths
    c = BUILD[name]()
    tr = Transient(c, toolkit=circuit.numeric)
    dt, res, paths = _timed(lambda: tr.solve(tend=STEPS * 2e-8, timestep=2e-8,
                                             fixed_timestep=True))
    return dt, 1e6 * dt / STEPS, _sha(res.x), _stats(tr), paths


def measure(cases, rounds):
    """Every case warmed once, then `rounds` interleaved runs; a dict per
    case with the times, the digest, the statistics and the path counts."""
    out = {k: {'times': [], 'per_step': [], 'shas': [], 'stats_all': [], 'paths_all': []}
           for k in cases}
    for k in cases:
        run_case(k)
    for _ in range(rounds):
        for k in cases:
            dt, us, sha, st, paths = run_case(k)
            r = out[k]
            r['times'].append(dt)
            r['per_step'].append(us)
            r['shas'].append(sha)
            r['stats_all'].append(st)
            r['paths_all'].append(paths)
    for k in cases:
        ## EVERY run's bytes and statistics, not the last one's: a run that
        ## differs from its own siblings is a finding before any timing
        r = out[k]
        same = (len(set(r['shas'])) == 1 and all(s == r['stats_all'][0] for s in r['stats_all'])
                and all(p == r['paths_all'][0] for p in r['paths_all']))
        r['sha'] = r['shas'][0] if same else 'INCONSISTENT:' + ','.join(sorted(set(r['shas'])))
        r['stats'] = r['stats_all'][0] if same else None
        r['paths'] = r['paths_all'][0] if same else None
        del r['shas'], r['stats_all'], r['paths_all']
    return out


def _key(r):
    return 'per_step' if r['per_step'] and r['per_step'][0] is not None else 'times'


def _fmt(k, r):
    if _key(r) == 'times':
        return f'{min(r["times"]):.4f} s (median {statistics.median(r["times"]):.4f})'
    return (f'{min(r["per_step"]):7.1f} us/step (median '
            f'{statistics.median(r["per_step"]):7.1f})')


## -- comparing two trees ------------------------------------------------------------

def _run_rounds(trees, cases, n, got, max_busy, repeats, plant, attempts_left):
    """`n` more kept rounds into `got`; `(kept, discarded, attempts used)`."""
    done = discarded = used = 0
    while done < n and used < attempts_left:
        used += 1
        r = len(got['parent'][cases[0]]['times']) + done
        order = list(trees) if r % 2 == 0 else list(trees)[::-1]
        this, busy = {}, 0.0
        for lab in order:
            tree_dir = trees[lab]
            ## THIS script in both trees (the parent may predate it): the
            ## package comes from PYTHONPATH, which precedes site-packages
            env = dict(os.environ, PYTHONPATH=tree_dir, PYCIRCUIT_BENCH_PIN='1')
            env.pop('PYCIRCUIT_BENCH_PLANT', None)
            if plant and lab == 'child':
                env['PYCIRCUIT_BENCH_PLANT'] = str(plant)
            ## (one memory layout and one hash order for every process: no
            ## address randomisation -- `setarch -R`, no root needed -- and a
            ## fixed hash seed; per-process layouts made the per-round
            ## ratios BIMODAL, -2.5 % or +1 %, and an A/A median of -1.7 %)
            env['PYTHONHASHSEED'] = '0'
            cmd = _bench.no_aslr() + [sys.executable, os.path.abspath(__file__),
                                      '--json', '--rounds', str(repeats), *cases]
            with _bench.Idle() as idle:
                res = subprocess.run(cmd, cwd=tree_dir, env=env, capture_output=True,
                                     text=True, check=True)
            busy = max(busy, idle.sibling_busy)
            this[lab] = [json.loads(ln) for ln in res.stdout.splitlines()
                         if ln.startswith('{')]
        if busy > max_busy:
            discarded += 1
            print(f'  round discarded: CPU {_bench.CPU}\'s sibling was {100 * busy:.0f} % '
                  'busy', flush=True)
            continue
        done += 1
        for lab, recs in this.items():
            for rec in recs:
                g = got[lab][rec['case']]
                ## each side of a round: the MIN of its process's warm runs
                g['times'].append(min(rec['times']))
                if rec['per_step'][0] is not None:
                    g['per_step'].append(min(rec['per_step']))
                g['sha'].add(rec['sha'])
                g['stats'].append(rec['stats'])
                g['paths'].append(rec.get('paths'))
    return done, discarded, used


def _summary(p, c):
    """`(key, same_bytes, same_stats, paths note, (med, lo, hi, wins) or None)`."""
    key = _key(p)
    same_bytes = len(p['sha'] | c['sha']) == 1
    same_stats = all(s == p['stats'][0] for s in p['stats'] + c['stats'])
    pp = [x for x in p['paths'] if x is not None]
    cp = [x for x in c['paths'] if x is not None]
    if not cp:
        note = 'paths: not counted'
    elif not pp:
        note = 'paths: the parent has no counts'
    elif all(x == pp[0] for x in pp) and all(x == cp[0] for x in cp) and pp[0] == cp[0]:
        note = 'paths SAME'
    else:
        a, b = pp[0], cp[0]
        moved = [f'{k} {a.get(k, 0)}->{b.get(k, 0)}' for k in sorted(set(a) | set(b))
                 if a.get(k, 0) != b.get(k, 0)]
        note = 'paths MOVED: ' + ('; '.join(moved[:6]) + (' ...' if len(moved) > 6 else '')
                                  if moved else 'between runs of one side')
    return key, same_bytes, same_stats, note, _bench.paired_summary(p[key], c[key])


#: `--check`'s verdict: SLOWER only when the whole interval sits above +1 %
#: and the median above +2 % (FASTER mirrored); UNSURE when the median is
#: beyond +-1.5 % without that -- a 2 % change may be there, and more rounds
#: are run while the budget lasts; "no change" otherwise (no change of 2 %
#: or more seen).  Measured 2026-10-04: five A/A checks, no verdict (the
#: short analysis cases' intervals reach +-2-4 %, their medians stay within
#: +-1.5 %); a planted +3 % made all ten cases SLOWER, a planted +1 % none.
VERDICT_MEDIAN = 0.02
VERDICT_EDGE = 0.01
VERDICT_UNSURE = 0.015


def verdict(med, lo, hi):
    if lo > 1 + VERDICT_EDGE and med > 1 + VERDICT_MEDIAN:
        return 'SLOWER'
    if hi < 1 - VERDICT_EDGE and med < 1 - VERDICT_MEDIAN:
        return 'FASTER'
    if abs(med - 1) > VERDICT_UNSURE:
        return 'UNSURE'
    return 'no change'


def compare(parent, cases, rounds, max_busy=0.25, repeats=5, child=None, check=False,
            budget=300.0, plant=0.0):
    """The tree at `parent` against `child` (this script's tree), each in its
    own interpreter from its own directory, pinned to one P-core, under the
    benchmark lock (`_bench`).  The order alternates per round (parent
    first, then child first), a round in which the timed core's hyperthread
    sibling was busy more than `max_busy` is discarded and re-run, and the
    result is the median of the per-round PAIRED ratios child/parent with a
    bootstrap 95 % interval -- drift during the run lands on both sides of a
    pair.  Each side of a round is the MIN of `repeats` warm runs in its
    process (timing noise only adds: one slow run made a whole round an
    outlier, and six such rounds an A/A interval of [-3, +20] %).  The bytes
    and statistics of every run are checked; any difference exits 2.  The
    fast-path counts of the timed calls are compared and reported (a speed
    change usually MOVES them; they never block a verdict).

    `check` (`--check`): a verdict per case (`verdict`); rounds are added,
    two at a time, while a case is UNSURE and `budget` seconds last.  Exit 3
    when a case is SLOWER."""
    child = child or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    trees = {'parent': os.path.abspath(parent), 'child': os.path.abspath(child)}
    got = {lab: {k: {'times': [], 'per_step': [], 'sha': set(), 'stats': [], 'paths': []}
                 for k in cases} for lab in trees}
    t_start = time.time()
    with _bench.lock(exclusive=True, what='step_machinery --check' if check else
                     'step_machinery --compare'):
        print('conditions:', _bench.stamp(), flush=True)
        if plant:
            print(f'PLANTED: the child is slowed by {100 * plant:.1f} % (the check\'s proof)',
                  flush=True)
        kept, discarded, used = _run_rounds(trees, cases, rounds, got, max_busy, repeats,
                                            plant, 3 * rounds)
        while check and time.time() - t_start < budget:
            pending = []
            for k in cases:
                s = _summary(got['parent'][k], got['child'][k])
                if s[4] is not None and verdict(*s[4][:3]) == 'UNSURE':
                    pending.append(k)
            if not pending:
                break
            per_round = (time.time() - t_start) / max(kept + discarded, 1)
            if time.time() - t_start + 2 * per_round > budget:
                break
            print(f'  unsure {", ".join(pending)}: two more rounds', flush=True)
            k2, d2, u2 = _run_rounds(trees, cases, 2, got, max_busy, repeats, plant, 6)
            kept, discarded, used = kept + k2, discarded + d2, used + u2
    print(f'{kept} rounds kept, {discarded} discarded; each side of a round the min of '
          f'{repeats} warm runs; {time.time() - t_start:.0f} s', flush=True)
    differ, slower, faster = False, [], []
    for k in cases:
        p, c = got['parent'][k], got['child'][k]
        key, same_bytes, same_stats, note, summ = _summary(p, c)
        differ = differ or not (same_bytes and same_stats)
        if summ is None:
            print(f'{k:12s} no rounds kept', flush=True)
            continue
        med, lo, hi, wins = summ
        v = ''
        if check:
            v = 'NOT-COMPARABLE' if not (same_bytes and same_stats) else verdict(med, lo, hi)
            slower += [k] if v == 'SLOWER' else []
            faster += [k] if v == 'FASTER' else []
            v = f'{v:14s} '
        print(f'{k:12s} {v}parent {_fmt(k, p)} | child {_fmt(k, c)} | paired '
              f'{100 * (med - 1):+.1f} % [{100 * (lo - 1):+.1f}, {100 * (hi - 1):+.1f}] '
              f'child faster {wins}/{len(p[key])} | bytes {"SAME" if same_bytes else "DIFFER"}'
              f', stats {"SAME" if same_stats else "DIFFER"}, {note}', flush=True)
        print('      per-round ratios: ' + ' '.join(
            f'{100 * (cv / pv - 1):+.1f}' for pv, cv in zip(p[key], c[key])), flush=True)
    if differ:
        print('BYTES OR STATISTICS DIFFER between or within the trees: not comparable',
              flush=True)
        sys.exit(2)
    if check:
        print(f'check: SLOWER {slower or "none"}; FASTER {faster or "none"}', flush=True)
        if slower:
            sys.exit(3)


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
    ## NOTHING IS WRAPPED ON AN INSTANCE (2026-10-03).  An instance wrap is
    ## an instance shadow, and the fast paths decline on shadows by design:
    ## on `cir` or `tr` the evaluate core declines and the Python path runs,
    ## on an hdl element the batch does, on a hand-written element its
    ## constant stamp is no longer recognised.  The tree once timed the
    ## Python path because of it.  So the wraps are on CLASSES (the core
    ## reads instance dicts only) and the elements are not wrapped at all;
    ## and the run is checked against the untimed one below.
    ref = run_case(name)
    from pycircuit.circuit.circuit import SubCircuit as _Sub
    for nm in ('solve_timestep', '_newton', '_residual_and_jacobian', '_predict_state',
               '_push_history', '_branch_after_solve', '_branch_screen',
               '_newton_abstol_vector', '_newton_xtol_vector', '_source_at', 'get_diff',
               '_C_lookup', '_companion_at', '_C_at_state'):
        wrap(Transient, nm)
    orig_nl = Transient._newton_limiter

    def newton_limiter(self):
        f = orig_nl(self)
        return None if f is None else timed(
            f, 'limiter_func (reinsert x2 + cir.limit + remove)')
    Transient._newton_limiter = newton_limiter
    for nm in ('i', 'q', 'G', 'C', 'u', 'limit', 'accept_step', 'next_event'):
        wrap(_Sub, nm, f'SubCircuit.{nm} (Python passes)')
    from pycircuit.circuit import _hdl_batch, _tran_core
    wrap(_hdl_batch.Batch, 'run', 'Batch.run (one C call per class per pass)')
    served = [0, 0]
    orig_eval = _tran_core.evaluate

    def evaluate(*a, **k):
        r = orig_eval(*a, **k)
        served[r is None] += 1
        return r
    _tran_core.evaluate = timed(evaluate, '_tran_core.evaluate (the passes + companion in C)')
    c = BUILD[name]()
    tr = Transient(c, toolkit=circuit.numeric)
    tr.solve(tend=4e-8, timestep=2e-8, fixed_timestep=True)
    for v in acc.values():
        v[0] = v[1] = 0
    served[:] = [0, 0]
    c = BUILD[name]()
    tr = Transient(c, toolkit=circuit.numeric)
    ls = tr._get_linearsolver()
    wrap(ls, 'solve', 'linsolver.solve')
    wrap(ls, 'factor', 'linsolver.factor')
    t0 = time.perf_counter_ns()
    res = tr.solve(tend=STEPS * 2e-8, timestep=2e-8, fixed_timestep=True)
    total = (time.perf_counter_ns() - t0) / STEPS
    sha = hashlib.sha256(np.asarray(res.x, float).tobytes()).hexdigest()[:12]
    same = sha == ref[2] and _stats(tr) == ref[3]
    print(f'== {name}: n={c.n} {total / 1e3:7.1f} us/step (solve wall / {STEPS}) -- '
          f'the timers RANK pieces, they do not size them (each adds ~1 us and '
          f'nests); size a piece with benchmarks/micro.py')
    print(f'   the timed run is the untimed run: {"YES" if same else "NO -- the timers changed the path"}'
          f'; the evaluate core served {served[0]} calls, declined {served[1]}')
    for lab, (ns, cnt) in sorted(acc.items(), key=lambda kv: -kv[1][0]):
        if cnt:
            print(f'   {lab:52s} {ns / STEPS / 1e3:7.1f} us/step {100 * ns / STEPS / total:5.1f} %'
                  f'  ({cnt / STEPS:6.2f} calls/step, {ns / cnt / 1e3:7.2f} us each)')


def _sampled_round(name):
    ## the frame `sample` keeps: only stacks through it are counted, so the
    ## imports, the compile-cache loads and the warm-up run drop out
    return run_case(name)


def _sample_child(name, rounds):
    run_case(name)
    for _ in range(rounds):
        _sampled_round(name)


def sample(name, rounds=20):
    """A sampled profile of `rounds` warm runs of `name` (py-spy in LAUNCH
    mode: the target is its child, which `ptrace_scope=1` permits), no
    wrapper in the code: the self time of each function, largest first,
    over the samples taken inside the measured runs."""
    import tempfile
    spy = os.path.join(os.path.dirname(sys.executable), 'py-spy')
    out = tempfile.mktemp(suffix='.txt')
    cmd = [spy, 'record', '--format', 'raw', '--rate', '997', '-o', out, '--',
           sys.executable, os.path.abspath(__file__), '--sample-child', str(rounds), name]
    ## (py-spy 0.4.2 writes the profile, then exits 1 with "No child
    ## process" when it reaps the child here: the file decides, not the code)
    subprocess.call(cmd, stdout=subprocess.DEVNULL)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        print('py-spy could not sample (ptrace refused?).  To allow it: '
              '`sudo setcap cap_sys_ptrace+ep $(readlink -f ' + spy + ')` or run '
              'under a session that permits ptrace of children.')
        return
    self_t, total, outside = {}, 0, 0
    with open(out) as f:
        for ln in f:
            stack, _, n = ln.rstrip().rpartition(' ')
            if not stack:
                continue
            n = int(n)
            ## inside a measured run AND inside its timed call (everything
            ## under `_timed`; the circuit's construction is not timed)
            if '_sampled_round (' not in stack or ';_timed (' not in stack:
                outside += n
                continue
            total += n
            leaf = stack.rsplit(';', 1)[-1]
            leaf = leaf.rsplit(':', 1)[0] + ')' if leaf.endswith(')') else leaf
            self_t[leaf] = self_t.get(leaf, 0) + n
    os.unlink(out)
    print(f'== {name}: {total} samples in the timed calls of {rounds} warm runs '
          f'({outside} elsewhere dropped), self time by function (sampled, no wrappers)')
    for fn, n in sorted(self_t.items(), key=lambda kv: -kv[1])[:30]:
        print(f'   {100 * n / total:5.1f} %  {fn}')


PARENT_TREE = os.environ.get('PYCIRCUIT_PARENT_TREE') or os.path.join(
    os.path.expanduser('~'), '.cache', 'pycircuit', 'wt_parent')


def main(argv):
    rounds = None
    as_json = False
    cmp_dir = None
    child = None
    check = False
    budget = 300.0
    plant = 0.0
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
        elif a == '--check':
            check = True
        elif a == '--parent':
            cmp_dir = next(it)
        elif a == '--child':
            child = next(it)
        elif a == '--budget':
            budget = float(next(it))
        elif a == '--plant':
            plant = float(next(it)) / 100.0
        elif a == '--tree':
            tree_case = next(it)
        elif a == '--sample':
            sample(next(it), rounds=rounds or 20)
            return
        elif a == '--sample-child':
            n = int(next(it))
            _sample_child(next(it), n)
            return
        else:
            cases.append(a)
    cases = cases or list(CASES)
    if tree_case:
        tree(tree_case)
        return
    if check:
        cmp_dir = cmp_dir or PARENT_TREE
        if not os.path.isdir(cmp_dir):
            sys.exit(f'no parent tree at {cmp_dir}: make it with scripts/parent_tree.sh')
    if cmp_dir:
        with warnings.catch_warnings():
            compare(cmp_dir, cases, rounds or (6 if check else 5), child=child, check=check,
                    budget=budget, plant=plant)
        return
    got = measure(cases, rounds or 5)
    for k in cases:
        r = got[k]
        if as_json:
            print(json.dumps({'case': k, **r}), flush=True)
        else:
            print(f'{k:5s} {_fmt(k, r)}  (digest {r["sha"]})', flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
