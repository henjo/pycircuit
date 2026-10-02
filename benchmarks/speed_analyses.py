"""Where the analyses spend their time (2026-10-01, the speed plan's P0).

The cases the plan was profiled on, each timed on its own and, with
`--profile`, counted: the circuit assembly passes (`SubCircuit.
_add_element_submatrices` / `_add_element_subvectors`), the reverse-replay
steps (`*.adjoint` of the factored period's steps) and the period walks.
Attribute a speed change by those counts and a profile, never by wall
clock alone -- this box is shared (`uptime` first).

    python benchmarks/speed_analyses.py                 # every case but PSP
    python benchmarks/speed_analyses.py ladder20 pnoise # by name
    python benchmarks/speed_analyses.py --psp           # + the PSP stage (minutes)
    python benchmarks/speed_analyses.py --profile       # + counts and a top list

Baseline 2026-10-01 (13729471, quiet box, no profiler): ladder20 gear 0.56
s, radau 0.50 s; ladder300 gear 2.82 s; vdP gear PSS + factored period 1.38
s; pnoise 20 x 8 sidebands 0.80 s; radau vdP first ppv 0.33 s; PSP gear PSS
16.3 s, radau 27.8 s (its factored period 6.1 s).  Under cProfile the
assembly (`_add_element_*`) was 73-86 % of the transient cases -- the
profiler inflates call-heavy code ~2-3x, but the per-pass timings agree:
~75 % unprofiled.
"""
import cProfile
import io
import pstats
import sys
import time
import warnings

import numpy as np

from pycircuit.circuit import circuit
from pycircuit.circuit.elements import C, Diode, R, SubCircuit, VSin, gnd
from pycircuit.circuit.integrator import Gear2Integrator, RadauIIA3Integrator
from pycircuit.circuit.shooting import PAC, PSS
from pycircuit.circuit.transient import Transient

circuit.default_toolkit = circuit.numeric


def ladder(n, diode_every=None):
    """An RC ladder driven by a sine, with diodes to ground: one at the
    end (`diode_every=None`) or one every `diode_every` sections."""
    c = SubCircuit()
    c['vs'] = VSin('n0', gnd, va=2.0, freq=1e3)
    for k in range(n):
        c[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
        c[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-8)
    if diode_every is None:
        c['D'] = Diode(f'n{n}', gnd)
    else:
        for k in range(0, n, diode_every):
            c[f'D{k}'] = Diode(f'n{k + 1}', gnd)
    return c


def vdp():
    from pycircuit.circuit.tests._shooting_fixtures import _vdp_with_noise
    cir = _vdp_with_noise(1e-6)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    return cir, x0


def case_ladder20(method):
    integ = Gear2Integrator() if method == 'gear' else RadauIIA3Integrator()
    tr = Transient(ladder(20), integrator=integ, reltol=1e-5)
    tr.solve(tend=3e-3, timestep=1e-5)
    return f'{tr.statistics.accepted_steps} steps'


def case_ladder300():
    tr = Transient(ladder(300, diode_every=10), integrator=Gear2Integrator(),
                   reltol=1e-5)
    tr.solve(tend=1e-3, timestep=1e-5)
    return f'{tr.statistics.accepted_steps} steps'


def case_vdp_pss_gear():
    cir, x0 = vdp()
    p = PSS(cir, method='gear', reltol=1e-12)
    p.solve(period=6.6634, timestep=6.6634 / 400, x0=x0, maxiterations=60)
    p.factored_period()
    return f'converged {p.converged}'


def case_pnoise():
    cir, x0 = vdp()
    p = PSS(cir, method='gear', reltol=1e-12)
    p.solve(period=6.6634, timestep=6.6634 / 400, x0=x0, maxiterations=60)
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / p.period
    fs = np.logspace(-3, -0.5, 20) * f0
    t0 = time.perf_counter()
    for f in fs:
        pac.pnoise(p, f, 0, maxsidebands=8)
    t1 = time.perf_counter()
    pac.pnoise(p, fs, 0, maxsidebands=8)        # one call (P4)
    t2 = time.perf_counter()
    return f'pnoise alone: a loop {t1 - t0:.2f} s, one array call {t2 - t1:.2f} s'


def case_radau_ppv():
    cir, x0 = vdp()
    p = PSS(cir, method='radau', reltol=1e-10)
    p.solve(period=6.6634, timestep=6.6634 / 200, x0=x0, maxiterations=60)
    t0 = time.perf_counter()
    p.ppv()
    return 'first ppv alone %.2f s' % (time.perf_counter() - t0)


def case_psp(method):
    from pycircuit.circuit.tests.test_shooting_pnoise import _cs_amp
    p = PSS(_cs_amp(2e-2, fnt=1.0), method=method, reltol=1e-8)
    p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=60)
    t0 = time.perf_counter()
    p.factored_period()
    return 'factored period alone %.2f s' % (time.perf_counter() - t0)


CASES = [('ladder20_gear', lambda: case_ladder20('gear')),
         ('ladder20_radau', lambda: case_ladder20('radau')),
         ('ladder300_gear', case_ladder300),
         ('vdp_pss_gear', case_vdp_pss_gear),
         ('pnoise', case_pnoise),
         ('radau_ppv', case_radau_ppv)]
PSP_CASES = [('psp_gear', lambda: case_psp('gear')),
             ('psp_radau', lambda: case_psp('radau'))]

COUNTED = ('_add_element_submatrices', '_add_element_subvectors', 'adjoint',
           '_walk_lmm', '_walk_stage', '_walk_glm', 'solve_timestep')


def run(name, fn, profile):
    if profile:
        pr = cProfile.Profile()
        t0 = time.perf_counter()
        pr.enable()
        note = fn()
        pr.disable()
    else:
        t0 = time.perf_counter()
        note = fn()
    dt = time.perf_counter() - t0
    print(f'{name:<16s} {dt:7.2f} s   {note}', flush=True)
    if not profile:
        return
    st = pstats.Stats(pr)
    counts = {}
    for (_f, _l, fn_), (_cc, nc, _tt, ct, _c) in st.stats.items():
        if fn_ in COUNTED:
            n0, c0 = counts.get(fn_, (0, 0.0))
            counts[fn_] = (n0 + nc, max(c0, ct))
    for k in COUNTED:
        if k in counts:
            calls, cum = counts[k]
            print(f'    {k:<26s} calls {calls:8d}  cum {cum:6.2f} s')
    s = io.StringIO()
    pstats.Stats(pr, stream=s).sort_stats('tottime').print_stats(12)
    for line in s.getvalue().split('\n'):
        if '/' in line or 'built-in' in line or 'method' in line:
            print('    ' + line.strip()[:140])


def main(argv):
    warnings.simplefilter('ignore')
    profile = '--profile' in argv
    names = [a for a in argv if not a.startswith('--')]
    cases = CASES + (PSP_CASES if '--psp' in argv else [])
    for name, fn in cases:
        if names and not any(n in name for n in names):
            continue
        run(name, fn, profile)


if __name__ == '__main__':
    main(sys.argv[1:])
