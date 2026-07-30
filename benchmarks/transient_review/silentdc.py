import sys, numpy as np
sys.path.insert(0,'/home/andreas/sources/pycircuit')
from pycircuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import R, C, VS, VSin, IS
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.dcanalysis import DC

# circuit whose DC op fails (floating node) but whose transient "works"
c = SubCircuit()
c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
c['R']  = R('a','b', r=1e3)
c['C']  = C('b', gnd, c=1e-7)
c['Cfl']= C('float', gnd, c=1e-9)     # node 'float' has no DC path -> singular
try:
    DC(c).solve(); print('DC converged')
except Exception as e:
    print('DC FAILS: %s: %s' % (type(e).__name__, e))

t = Transient(c)
res = t.solve(tend=2e-3, timestep=1e-5)
print('Transient returned normally, %d points, v(b)[-1]=%.4g' %
      (len(res.sweep_values), res.v('b')[-1]))
print('  -> no warning, no exception, x0 silently became all-zeros')

# count full assemblies per accepted step
import pycircuit.circuit.circuit as cm
n_calls = {'G':0,'C':0,'i':0,'q':0,'u':0}
SC = cm.SubCircuit
om, ov = SC._add_element_submatrices, SC._add_element_subvectors
def wm(self, name, *a, **k):
    n_calls[name]+=1; return om(self, name, *a, **k)
def wv(self, name, *a, **k):
    n_calls[name]+=1; return ov(self, name, *a, **k)
SC._add_element_submatrices, SC._add_element_subvectors = wm, wv
c2 = SubCircuit()
c2['vs']=VSin('a',gnd,va=1.0,freq=1e3); c2['R']=R('a','b',r=1e3); c2['C']=C('b',gnd,c=1e-7)
r2 = Transient(c2).solve(tend=2e-3, timestep=1e-5)
SC._add_element_submatrices, SC._add_element_subvectors = om, ov
ns = len(r2.sweep_values)
print('\n%d accepted steps; assemblies: %s' % (ns, n_calls))
print('per step: ' + ', '.join('%s=%.2f'%(k,v/ns) for k,v in n_calls.items()))
