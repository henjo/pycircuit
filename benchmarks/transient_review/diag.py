"""What does the user see when a run goes wrong, and what does DC vabstol cost?"""
import sys, time, traceback
import numpy as np
sys.path.insert(0, '/home/andreas/sources/pycircuit')

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.elements import R, C, VS, VSin, IS, Diode, VPulse
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
import pycircuit.circuit.nrsolver as nrs

print('=== 1. non-convergent DC: what is reported? ===')
c = SubCircuit()
c['vs'] = VS('a', gnd, v=5.0)
c['R'] = R('a', 'b', r=1e-3)
c['D1'] = Diode('b', gnd)   # forward-biased hard through a tiny resistor
try:
    DC(c).solve()
    print('  converged (no failure to show)')
except Exception as e:
    print('  %s: %s' % (type(e).__name__, e))

print('\n=== 2. singular matrix (floating node): what is reported? ===')
c = SubCircuit()
c['vs'] = VS('a', gnd, v=1.0)
c['R'] = R('a', 'b', r=1e3)
c['C'] = C('c', gnd, c=1e-9)     # node c totally floating
try:
    DC(c).solve()
    print('  converged?!')
except Exception as e:
    print('  %s: %s' % (type(e).__name__, e))

print('\n=== 3. Newton iteration count vs vabstol on the leapfrog DC ===')
def build(stages=5):
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        build_leapfrog_network(cir, stages=stages, tones=((1e-3, 8e3),))
        return cir
    finally:
        circuit_module.default_toolkit = saved

for vab in (1e-12, 1e-6):
    cir = build()
    counts = []
    orig = nrs.StandardNewton.solve_system
    def spy(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None):
        r = orig(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler)
        counts.append(r[1])
        return r
    nrs.StandardNewton.solve_system = spy
    t0 = time.perf_counter()
    try:
        res = DC(cir, vabstol=vab).solve()
        ok = 'ok'
    except Exception as e:
        ok = 'FAIL %s' % e
    dt = time.perf_counter() - t0
    nrs.StandardNewton.solve_system = orig
    print('  vabstol=%.0e: %s, newton iters %s, %.3f s' % (vab, ok, counts, dt))

print('\n=== 4. is there any run statistic on the result object? ===')
c = SubCircuit()
c['IS'] = IS(gnd, 'n1', i=1e-3)
c['R'] = R('n1', gnd, r=1e3)
c['C'] = C('n1', gnd, c=1e-9)
tr = Transient(c)
res = tr.solve(tend=1e-5, timestep=1e-7)
print('  result attrs:', [a for a in dir(res) if not a.startswith('_')])
print('  analysis attrs:', sorted(a for a in vars(tr) if not a.startswith('__')))
print('  n accepted steps:', len(res.sweep_values))

print('\n=== 5. Transient parameter surface ===')
for p in tr.parameters:
    print('   %-12s default=%-24r %s' % (p.name, p.default, p.desc))

print('\n=== 6. DC parameter surface ===')
d = DC(c)
for p in d.parameters:
    print('   %-12s default=%-24r %s' % (p.name, p.default, p.desc))
