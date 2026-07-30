"""Compound effect: fast assembly + drop the redundant post-Newton func(x)."""
import sys, time
import numpy as np
sys.path.insert(0,'/home/andreas/sources/pycircuit')
sys.path.insert(0,'/tmp/claude-1000/-home-andreas-pycircuit-claude/e98a2fae-be82-4481-bb9a-0deaaae9b6e6/scratchpad')
import pycircuit.circuit.circuit as cm
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit as SC
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.transient import Transient
import pycircuit.circuit.nrsolver as nrs
from proto_assembly import fast_subvectors, fast_submatrices  # reuse verified impls

def build():
    saved=cm.default_toolkit; cm.default_toolkit=numeric
    try:
        cir=SubCircuit(toolkit=numeric)
        _,amp,_=build_leapfrog_network(cir,stages=5,tones=((1e-3,8e3),))
        cir['nl']=BSource(amp[0]['e1'],gnd,amp[0]['e1'],gnd,i_func=lambda v:1e-3*v**3)
        return cir
    finally: cm.default_toolkit=saved

def run(tag):
    c=build(); t0=time.perf_counter()
    r=Transient(c,toolkit=numeric).solve(refnode=gnd,tend=20e-6,timestep=1.0/(8e3*128))
    dt=time.perf_counter()-t0
    print('%-34s %6.2f s  (%d steps)' % (tag, dt, len(r.sweep_values)), flush=True)
    return dt, np.asarray(r.x)

base, xb = run('baseline')

orig_v, orig_m = SC._add_element_subvectors, SC._add_element_submatrices
SC._add_element_subvectors, SC._add_element_submatrices = fast_subvectors, fast_submatrices
a, xa = run('+ fast assembly')

# also drop the redundant post-Newton re-evaluation: make solve_system hand back F,J
orig_ss = nrs.StandardNewton.solve_system
def ss(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter=None, scaler=None):
    r = orig_ss(self, x0, eval_FJ, toolkit, reltol, abstol, xtol, maxiter, limiter, scaler)
    self.last_FJ = eval_FJ  # placeholder; real fix returns cached F,J
    return r
orig_st = Transient.solve_timestep
def solve_timestep(self, x0, t, refnode=gnd, provided_function=None):
    cache = {}
    dt = self._dt
    def func(x):
        C = self.cir.C(x); q = self.cir.q(x)
        iq, Geq = self.get_diff(q, C)
        f = self.cir.i(x) + iq + self.cir.u(t, analysis=self.par.analysis)
        J = self.cir.G(x) + Geq
        f = self.toolkit.array(f, dtype=float); J = self.toolkit.array(J, dtype=float)
        cache['f'], cache['J'], cache['C'], cache['q'] = f, J, C, q
        return f, J
    x = self._newton(func, x0)
    f, J = cache['f'], cache['J']     # last evaluated, not re-evaluated
    return x, None, J, f
Transient.solve_timestep = solve_timestep
b, xc = run('+ no redundant func(x)')
Transient.solve_timestep = orig_st
SC._add_element_subvectors, SC._add_element_submatrices = orig_v, orig_m

print('\nspeedups vs baseline: assembly %.2fx, +no-redundant %.2fx' % (base/a, base/b))
print('waveform drift vs baseline: fast=%.2e  fast+noredundant=%.2e (rel to %.3g)'
      % (np.max(np.abs(xa-xb)), np.max(np.abs(xc-xb)), np.max(np.abs(xb))))
