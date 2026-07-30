import pycircuit.circuit._jaxtoolkit as jax_toolkit
from pycircuit.circuit import circuit
circuit.default_toolkit = jax_toolkit

from pycircuit.circuit import R, VS, gnd
from pycircuit.circuit.elements import SubCircuit, TLine
from pycircuit.circuit.func import Pulse
from pycircuit.circuit.jaxtransient import JAXTransient

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')
cir['V1'] = VS('in', gnd, v=0)
cir['V1'].function = Pulse(0, 1, 1e-9, 1e-12, 1e-12, 10e-9, 20e-9, toolkit=jax_toolkit)
cir['T1'] = TLine('in', gnd, 'out', gnd, Z0=50.0, TD=1e-9)
cir['R1'] = R('out', gnd, r=50.0)

try:
    res = JAXTransient(cir).solve(gnd, tend=5e-9, timestep=1e-10)
    print("TLine worked in JAXTransient!")
except Exception as e:
    import traceback
    traceback.print_exc()
