import pycircuit.circuit._jaxtoolkit as jax_toolkit
from pycircuit.circuit import circuit, R, VS, gnd
circuit.default_toolkit = jax_toolkit
from pycircuit.circuit.elements import SubCircuit, TLine
from pycircuit.circuit.dcanalysis import DC
import numpy as np

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')
cir['V1'] = VS('in', gnd, v=1.0)
cir['T1'] = TLine('in', gnd, 'out', gnd, Z0=50.0, TD=1e-9)
cir['R1'] = R('out', gnd, r=50.0)

x0 = np.zeros(cir.n)
print("u=", cir.u(0, analysis='dc'))
print("G=", cir.G(x0))
