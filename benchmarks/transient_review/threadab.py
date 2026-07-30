import sys, time
import numpy as np
sys.path.insert(0,'/home/andreas/sources/pycircuit')
import pycircuit.circuit.circuit as cm
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.transient import Transient
saved=cm.default_toolkit; cm.default_toolkit=numeric
cir=SubCircuit(toolkit=numeric)
_,amp,_=build_leapfrog_network(cir,stages=5,tones=((1e-3,8e3),))
cir['nl']=BSource(amp[0]['e1'],gnd,amp[0]['e1'],gnd,i_func=lambda v:1e-3*v**3)
cm.default_toolkit=saved
t0=time.perf_counter()
r=Transient(cir,toolkit=numeric).solve(refnode=gnd,tend=20e-6,timestep=1.0/(8e3*128))
print('steps=%d wall=%.2f s  %.4f s/step' % (len(r.sweep_values), time.perf_counter()-t0,
      (time.perf_counter()-t0)/len(r.sweep_values)))
