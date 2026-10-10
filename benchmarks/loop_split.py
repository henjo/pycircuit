"""What a C stepping loop could save on the SPICE gold decks, in
instructions AND cycles (the C stepping loop plan's stage 0,
`~/.claude/plans/c-stepping-loop.md`, 2026-10-10).

    python benchmarks/loop_split.py                    # every deck of DECKS
    python benchmarks/loop_split.py gm6 ring           # by name
    python benchmarks/loop_split.py --json out.json    # the records too
    python benchmarks/loop_split.py --plant accept 20000 gm6   # the proof

Each deck runs in a process of its own, pinned to one P-core with one BLAS
thread, warm (the deck solved once before the timed solve).  In the timed
solve the CPU's counters -- instructions and cycles, user space -- are read
at the boundaries of `_SteppingLoop`'s five phases and of every COMPUTE
call: each compiled entry point's cached `_driver` function (the C Newton,
the C error test, the C evaluation core, the batched HDL kernels), a
`CKernel` call, the leaf linear solvers (dense, chord, SuperLU, KLU) and
dense LAPACK reached through `numpy.linalg` / `scipy.linalg`.
Compute calls nest (a solver inside a kernel's caller): only the outermost
is counted.

    outside      the solve outside the loop's phases (the operating point,
                 the result) -- not a C loop's to remove
    compute      inside the loop, inside a compute call -- what C costs
                 and keeps costing
    removable    inside the loop, in no compute call: the Python (and the
                 numpy it drives) a C loop replaces
    ceiling      the run with `removable` gone, 1 / (1 - removable)

THE READINGS COST WHAT THEY MEASURE, so they are calibrated in the same
process before the timed solve -- a counted call of an empty function,
nested as the real ones are -- and taken off: each counted call's cost to
the region around it and to its own region.  The whole method is proved by
`--plant PHASE N`: a busy loop of known cost added to one phase must come
back in that phase, to 2 % in instructions (the deck run with and without
it, the difference against the loop's standalone count times its calls).
"""
import argparse
import json
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

#: the gold decks (n 15 .. 429), and the share of each one's `.tran` run
#: (the long ones capped near 2 s native, hundreds of steps each)
DECKS = {'4049osc': 1, 'toronto': 1, 'ab_integ': 0.35, 'e1480': 1, 'g1310': 1, 'gm6': 1,
         'hussamp': 1, 'mosrect': 1, 'mux8': 0.07, 'nand': 1, 'ring': 0.4, 'slowlatch': 1,
         'rca': 1, 'schmitecl': 1, 'latch': 1, 'gm1': 1, 'gm2': 1, 'gm3': 1, 'gm17': 1,
         'mike2': 1, 'rich3': 0.15, 'todd3': 1, 'arom': 0.15, 'voter25': 0.2, 'jge': 0.05,
         'gm19': 0.05}
PHASES = ('where', 'take', 'judge', 'event', 'accept')
#: the compiled entry points, each module's cached `(ffi, cfn, ...)` driver
#: -- not `_tran_core`'s, whose function the C Newton takes as a POINTER
#: (a Python stand-in breaks that cast): its Python callers instead,
#: `CORE_CALLERS`, their packing counted as compute (conservative) -- and only
#: on a deck the C Newton never serves (wrapped, they are not genuine and the
#: C Newton declines: `newton_c:patched`); where it serves, the core runs
#: inside it and is counted with it
DRIVERS = ('_tran_newton_c', '_tran_lte_c', '_hdl_batch')
CORE_CALLERS = ('evaluate', 'passes')
#: dense LAPACK reached from Python by module attribute (the branch check's
#: per-step SVD screen was 51 % of ab_integ's run, uncounted, until these)
#: -- not `numpy.linalg.solve`: the C error test twins it and declines when
#: it is not genuine (`lte_c:patched`); dense solves are counted through the
#: solver classes instead
LAPACK = {'numpy.linalg': ('svd', 'inv', 'eig', 'eigh', 'eigvals', 'eigvalsh', 'lstsq',
                           'qr', 'matrix_rank', 'det', 'slogdet', 'pinv'),
          'scipy.linalg': ('lu_factor', 'lu_solve', 'svd', 'solve', 'qr', 'eig', 'eigvals',
                           'lstsq', 'svdvals', 'null_space')}
#: the leaf linear solvers (`AutoSolver` only delegates to these)
SOLVERS = (('DenseSolver', ('solve', 'factor')), ('SuperLUSolver', ('solve', 'factor')),
           ('KLUSolver', ('solve', 'factor')), ('NumpyLU', ('solve',)),
           ('_Factored', ('solve', 'solve_transposed')),
           ('_KLUNumeric', ('solve', 'solve_transposed')))


#: the independent sources the C loop's source table holds (DC, PULSE, PWL,
#: SIN, EXP; V and I) -- any other `VS`/`IS` declines the run
SPICE_SOURCES = {'VS', 'IS', 'VSin', 'ISin', 'VPulse', 'IPulse', 'VPWL', 'IPWL', 'VExp', 'IExp'}
LMM_INTEGRATORS = {'Gear2Integrator', 'TrapezoidalIntegrator', 'EulerIntegrator', 'ThetaIntegrator'}


def _leaves(cir):
    for el in cir.elements.values():
        if getattr(el, 'elements', None):
            yield from _leaves(el)
        else:
            yield el


def census(tr, family, provided_function):
    """Why the C loop would NOT run this transient (the plan's scope), as a
    sorted list of reasons -- empty when it would run."""
    from pycircuit.circuit import circuit as _c
    from pycircuit.circuit._limiting import stateful_limiters
    from pycircuit.circuit.elements import IS, VS
    why = []
    if type(family).__name__ != '_LMMSteps':
        why.append('family:' + type(family).__name__)
    if type(tr.base_integrator).__name__ not in LMM_INTEGRATORS and not why:
        why.append('integrator:' + type(tr.base_integrator).__name__)
    ctrl = getattr(tr, 'step_controller', None)
    if ctrl is not None and type(ctrl).__name__ != 'IntegralController':
        why.append('controller:' + type(ctrl).__name__)
    elif ctrl is not None and (ctrl.lte_gamma_min, ctrl.lte_gamma_max, ctrl.lte_eta) != (0.0, 1.0, None):
        why.append('band')
    if getattr(tr, '_ev_rows', None) is not None:
        why.append('state_events')
    if provided_function is not None:
        why.append('provided_function')
    if stateful_limiters(tr.cir):
        why.append('stateful_limiter')
    hooks = sorted({type(el).__name__ for _i, el in tr.cir._hook_elements('accept_step')})
    if hooks:
        why.append('accept_step:' + '+'.join(hooks))
    odd = sorted({type(el).__name__ for el in _leaves(tr.cir)
                  if isinstance(el, (VS, IS)) and type(el).__name__ not in SPICE_SOURCES})
    if odd:
        why.append('source:' + '+'.join(odd))
    if tr.toolkit is not _c.numeric:
        why.append('toolkit')
    return sorted(why)


def census_hook(record):
    """Wrap `_SteppingLoop.__init__` so every stepping loop the process runs
    is censused: `record(tr, reasons)`."""
    from pycircuit.circuit import transient as T
    init = T._SteppingLoop.__init__

    def counted_init(self, tr, family, run, X, tend, timestep, provided_function):
        init(self, tr, family, run, X, tend, timestep, provided_function)
        try:
            record(tr, census(tr, family, provided_function))
        except Exception as e:  # noqa: BLE001 -- a census must never break a run
            record(tr, ['census-error:' + type(e).__name__])
    T._SteppingLoop.__init__ = counted_init


class Counters:
    """Instructions and cycles of this process, user space, on the P-core's
    PMU (`_bench.InstrCounter`'s attributes, a second event for cycles).
    `read()` is (instructions, cycles); `whole` turns False if either was
    multiplexed -- such a run is rejected, not scaled."""

    def __init__(self):
        import ctypes

        import _bench
        libc = ctypes.CDLL(None, use_errno=True)
        libc.syscall.restype = ctypes.c_long
        self._libc, self.fds, self.whole = libc, [], True
        pmu = _bench._core_pmu_type()
        for event in (_bench.PERF_COUNT_HW_INSTRUCTIONS, 0):   # 0: CPU cycles
            config = event | (pmu << 32 if pmu is not None else 0)
            buf = ctypes.create_string_buffer(_bench._perf_attr(config), 64)
            fd = libc.syscall(ctypes.c_long(_bench._NR_PERF_EVENT_OPEN), buf, ctypes.c_int(0),
                              ctypes.c_int(-1), ctypes.c_int(-1), ctypes.c_ulong(8))
            if fd < 0:
                raise OSError(ctypes.get_errno(), 'perf_event_open')
            self.fds.append(fd)

    def start(self):
        import _bench
        for fd in self.fds:
            self._libc.ioctl(fd, _bench.PERF_EVENT_IOC_RESET, 0)
            self._libc.ioctl(fd, _bench.PERF_EVENT_IOC_ENABLE, 0)

    def read(self):
        out = []
        for fd in self.fds:
            v, en, run = struct.unpack('QQQ', os.read(fd, 24))
            if run != en:
                self.whole = False
            out.append(v)
        return out[0], out[1]


class Split:
    """The counted regions of one run: per phase, per compute label, and
    the call counts the calibration is taken off by."""

    def __init__(self, ctr):
        self.ctr = ctr
        self.acc = {}           # label -> [instr, cycles, calls]
        self.loop_depth = 0
        self.compute_depth = 0
        self.plant = None       # (phase, function)

    def add(self, label, d0, d1):
        a = self.acc.setdefault(label, [0, 0, 0])
        a[0] += d0
        a[1] += d1
        a[2] += 1

    def phase(self, f, name):
        def counted(*a, **k):
            i0, c0 = self.ctr.read()
            self.loop_depth += 1
            try:
                if self.plant is not None and self.plant[0] == name:
                    self.plant[1]()
                return f(*a, **k)
            finally:
                self.loop_depth -= 1
                i1, c1 = self.ctr.read()
                self.add('phase.' + name, i1 - i0, c1 - c0)
        return counted

    def compute(self, f, label):
        def counted(*a, **k):
            if self.compute_depth:
                return f(*a, **k)
            i0, c0 = self.ctr.read()
            self.compute_depth = 1
            try:
                return f(*a, **k)
            finally:
                self.compute_depth = 0
                i1, c1 = self.ctr.read()
                self.add(('compute.' if self.loop_depth else 'outside-compute.') + label,
                         i1 - i0, c1 - c0)
        return counted


def install(split, core_callers):
    """The phase and compute wraps -- on the classes and module drivers the
    run uses, before it starts."""
    import importlib

    from pycircuit.circuit import _hdl_cbackend, linearsolver
    from pycircuit.circuit import transient as T
    for ph in PHASES:
        setattr(T._SteppingLoop, ph, split.phase(getattr(T._SteppingLoop, ph), ph))
    for modname in DRIVERS:
        mod = importlib.import_module('pycircuit.circuit.' + modname)
        drv = mod._driver or mod.driver()
        if drv:
            mod._driver = (drv[0], split.compute(drv[1], modname)) + tuple(drv[2:])
    from pycircuit.circuit import _tran_core
    for name in CORE_CALLERS if core_callers else ():
        setattr(_tran_core, name, split.compute(getattr(_tran_core, name), '_tran_core.' + name))
    _hdl_cbackend.CKernel.__call__ = split.compute(_hdl_cbackend.CKernel.__call__, 'CKernel')
    for modname, names in LAPACK.items():
        mod = importlib.import_module(modname)
        for name in names:
            if hasattr(mod, name):
                setattr(mod, name, split.compute(getattr(mod, name), modname + '.' + name))
    for cls, methods in SOLVERS:
        klass = getattr(linearsolver, cls, None)
        for m in methods if klass is not None else ():
            if m in vars(klass):
                setattr(klass, m, split.compute(vars(klass)[m], 'solver.' + cls))


def calibrate(split, rounds=2000):
    """The readings' own cost, as the run pays it: (a) what one counted
    compute call adds to the region AROUND it, (b) what it records INSIDE
    itself, (c) what one counted phase adds around it -- each per call, from
    an empty function, best of three."""
    ctr = split.ctr

    def empty():
        return None
    inner = split.compute(empty, '_cal')
    outer = split.phase(lambda: None, '_cal')
    best = None
    for _ in range(3):
        split.acc.clear()
        split.loop_depth = 1
        i0, c0 = ctr.read()
        for _ in range(rounds):
            empty()
        i1, c1 = ctr.read()
        for _ in range(rounds):
            inner()
        i2, c2 = ctr.read()
        for _ in range(rounds):
            outer()
        i3, c3 = ctr.read()
        rec = split.acc.get('compute._cal', [0, 0, 0])
        a = ((i2 - i1) - (i1 - i0) - rec[0]) / rounds, ((c2 - c1) - (c1 - c0) - rec[1]) / rounds
        b = rec[0] / rounds, rec[1] / rounds
        ph = split.acc.get('phase._cal', [0, 0, 0])
        c = ((i3 - i2) - (i1 - i0) - ph[0]) / rounds, ((c3 - c2) - (c1 - c0) - ph[1]) / rounds
        cand = (a, b, c)
        if best is None or a[0] + b[0] + c[0] < sum(x[0] for x in best):
            best = cand
    split.acc.clear()
    split.loop_depth = 0
    return best


def child(deck, fraction, plant):
    sys.path.insert(0, HERE)
    import _bench
    _bench.pin_threads()
    _bench.pin_cpu()
    sys.path.insert(0, ROOT)
    import warnings

    import spice_suite

    from pycircuit._testing import benchdata
    from pycircuit.circuit import circuit, spice_import
    circuit.default_toolkit = circuit.numeric
    path = spice_suite.CASES[deck][0]
    over = spice_suite.OVERRIDES.get(deck, {})

    def build():
        imp = spice_import.import_netlist(benchdata.spice_data(path))
        tr, kw = imp.transient(**over)
        kw['tend'] = kw['tend'] * fraction
        return tr, kw
    from pycircuit.circuit import _paths

    def paths_of_a_run():
        tr, kw = build()
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr.solve(**kw)
        return {k: v for k, v in _paths.since(before).items() if not k.startswith('once:')}
    paths_of_a_run()                            # warm: compiles, caches
    plain = paths_of_a_run()                    # the fast paths an unwrapped run takes
    ctr = Counters()
    split = Split(ctr)
    install(split, core_callers=not plain.get('newton_c:served', 0))
    ctr.start()
    cal = calibrate(split)
    plant_cost = None
    if plant is not None:
        n = int(plant[1])

        def busy():
            s = 0
            for k in range(n):
                s += k
            return s
        costs = []
        for _ in range(5):
            i0, _c0 = ctr.read()
            busy()
            i1, _c1 = ctr.read()
            costs.append(i1 - i0)
        plant_cost = sorted(costs)[2]
        split.plant = (plant[0], busy)
    tr, kw = build()
    split.acc.clear()
    before = _paths.snapshot()
    i0, c0 = ctr.read()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(**kw)
    i1, c1 = ctr.read()
    counted = {k: v for k, v in _paths.since(before).items() if not k.startswith('once:')}
    ## A WRAP CAN CHANGE THE PATH (a fast path checks the genuineness of what
    ## it stands in for): every fast-path count compared with the unwrapped run
    moved = {k: (plain.get(k, 0), counted.get(k, 0)) for k in set(plain) | set(counted)
             if plain.get(k, 0) != counted.get(k, 0)}
    ## (the C Newton declining for another REASON on a deck it never serves is
    ## the same path: `newton_c:patched` -- the core callers wrapped -- in
    ## place of the run's own reason, the same count)
    declines = {k: v for k, v in moved.items() if k.startswith('newton_c:')}
    if (declines == moved and not plain.get('newton_c:served', 0)
            and sum(v[0] for v in declines.values()) == sum(v[1] for v in declines.values())):
        moved = {}
    if not ctr.whole:
        return {'deck': deck, 'rejected': 'multiplexed'}
    acc = split.acc
    (ao_i, ao_c), (ai_i, ai_c), (ap_i, ap_c) = cal
    n_comp_loop = sum(v[2] for k, v in acc.items() if k.startswith('compute.'))
    n_comp_out = sum(v[2] for k, v in acc.items() if k.startswith('outside-compute.'))
    n_phase = sum(v[2] for k, v in acc.items() if k.startswith('phase.'))
    rec = {'deck': deck, 'n': int(tr.cir.n), 'steps': int(tr.statistics.accepted_steps),
           'paths_moved': moved,
           'attempts': int(tr.statistics.accepted_steps + tr.statistics.rejected_steps),
           'calibration': {'around': [ao_i, ao_c], 'inside': [ai_i, ai_c], 'phase': [ap_i, ap_c]},
           'calls': {'compute_in_loop': n_comp_loop, 'compute_outside': n_comp_out, 'phases': n_phase}}
    for j, unit in ((0, 'instr'), (1, 'cycles')):
        corr = (ao_i, ao_c)[j], (ai_i, ai_c)[j], (ap_i, ap_c)[j]
        total = (i1 - i0, c1 - c0)[j] - n_phase * corr[2] - (n_comp_loop + n_comp_out) * (corr[0] + corr[1])
        loop = (sum(v[j] for k, v in acc.items() if k.startswith('phase.'))
                - n_comp_loop * (corr[0] + corr[1]))
        compute = sum(v[j] for k, v in acc.items() if k.startswith('compute.')) - n_comp_loop * corr[1]
        removable = loop - compute
        rec[unit] = {'total': total, 'outside': (total - loop) / total, 'compute': compute / total,
                     'removable': removable / total, 'ceiling': total / (total - removable),
                     'phases': {p: (acc.get('phase.' + p, [0, 0, 0])[j]) / total for p in PHASES},
                     'compute_by': {k.split('.', 1)[1]: v[j] / total for k, v in acc.items()
                                    if k.startswith('compute.')}}
    if plant is not None:
        rec['plant'] = {'phase': plant[0], 'calls': acc.get('phase.' + plant[0], [0, 0, 0])[2],
                        'cost': plant_cost, 'phase_instr': acc.get('phase.' + plant[0], [0, 0, 0])[0]}
    return rec


def census_child(deck, fraction):
    sys.path.insert(0, HERE)
    sys.path.insert(0, ROOT)
    import warnings

    import spice_suite

    from pycircuit._testing import benchdata
    from pycircuit.circuit import _paths, circuit, spice_import
    circuit.default_toolkit = circuit.numeric
    got = []
    census_hook(lambda tr, why: got.append(why))
    imp = spice_import.import_netlist(benchdata.spice_data(spice_suite.CASES[deck][0]))
    tr, kw = imp.transient(**spice_suite.OVERRIDES.get(deck, {}))
    kw['tend'] = kw['tend'] * fraction
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr.solve(**kw)
    p = _paths.since(before)
    newton = {k: v for k, v in p.items() if k.startswith('newton_c:') and k not in
              ('newton_c:iterations', 'newton_c:jfold')}
    return {'deck': deck, 'n': int(tr.cir.n), 'why': got[0] if got else ['no-loop'],
            'newton_c': newton, 'chord': p.get('chord:chosen', 0)}


def run_child(deck, fraction, plant=None):
    cmd = [sys.executable, os.path.abspath(__file__), '--child', deck, str(fraction)]
    if plant is not None:
        cmd += ['--plant'] + [str(p) for p in plant]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT, check=False)
    line = [ln for ln in out.stdout.splitlines() if ln.startswith('{')]
    if not line:
        return {'deck': deck, 'rejected': 'failed: ' + out.stderr.strip()[-300:]}
    return json.loads(line[-1])


def show(rec):
    if 'rejected' in rec:
        print(f"{rec['deck']:10s} REJECTED: {rec['rejected']}")
        return
    i, c = rec['instr'], rec['cycles']
    if rec.get('paths_moved'):
        print(f"{rec['deck']:10s} ⚠ the counted run took other paths: {rec['paths_moved']}")
    print(f"{rec['deck']:10s} n {rec['n']:3d}  steps {rec['steps']:5d}  "
          f"outside {100 * c['outside']:5.1f}%  compute {100 * c['compute']:5.1f}%  "
          f"removable {100 * c['removable']:5.1f}% (instr {100 * i['removable']:5.1f}%)  "
          f"ceiling {c['ceiling']:4.2f}x (instr {i['ceiling']:4.2f}x)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('decks', nargs='*')
    ap.add_argument('--json')
    ap.add_argument('--plant', nargs=2, metavar=('PHASE', 'N'))
    ap.add_argument('--child', nargs=2, metavar=('DECK', 'FRACTION'), help=argparse.SUPPRESS)
    ap.add_argument('--census', action='store_true',
                    help='why the C loop would not run each deck (the plan\'s scope)')
    a = ap.parse_args(argv)
    if a.child and a.census:
        print(json.dumps(census_child(a.child[0], float(a.child[1]))))
        return 0
    if a.child:
        print(json.dumps(child(a.child[0], float(a.child[1]), a.plant)))
        return 0
    if a.census:
        for d in a.decks or list(DECKS):
            out = subprocess.run([sys.executable, os.path.abspath(__file__), '--census', '--child',
                                  d, str(DECKS.get(d, 1))], capture_output=True, text=True,
                                 cwd=ROOT, check=False)
            line = [ln for ln in out.stdout.splitlines() if ln.startswith('{')]
            r = json.loads(line[-1]) if line else {'deck': d, 'why': ['failed'], 'n': 0}
            print(f"{r['deck']:10s} n {r['n']:3d}  {'RUNS' if not r['why'] else 'declines: ' + ', '.join(r['why'])}"
                  f"  | C Newton {r.get('newton_c')}  chord {r.get('chord')}")
        return 0
    import numpy as np
    decks = a.decks or list(DECKS)
    if a.plant:
        deck = decks[0]
        base = run_child(deck, DECKS[deck])
        plant = run_child(deck, DECKS[deck], a.plant)
        p = plant['plant']
        got = plant['plant']['phase_instr'] - base['instr']['phases'][a.plant[0]] * base['instr']['total']
        want = p['calls'] * p['cost']
        print(f"plant {a.plant[0]} x{p['calls']} calls at {p['cost']} instr: expected {want:.4g}, "
              f"found {got:.4g} ({100 * (got / want - 1):+.2f} %)")
        return 0 if abs(got / want - 1) < 0.02 else 1
    recs = [run_child(d, DECKS.get(d, 1)) for d in decks]
    for r in recs:
        show(r)
    ok = [r for r in recs if 'rejected' not in r]
    if ok:
        for unit in ('cycles', 'instr'):
            c = np.array([r[unit]['ceiling'] for r in ok])
            small = [r[unit]['ceiling'] for r in ok if r['n'] <= 100]
            print(f"{unit:6s}: ceiling median {np.median(c):.2f}x over {len(ok)} decks "
                  f"(quartiles {np.percentile(c, 25):.2f} / {np.percentile(c, 75):.2f}); "
                  f"n <= 100: {np.median(small):.2f}x over {len(small)}")
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(recs, fh, indent=1)
    return 0


if __name__ == '__main__':
    sys.exit(main())
