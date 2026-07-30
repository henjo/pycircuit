import pycircuit.circuit._jaxtoolkit as jax_toolkit
from pycircuit.circuit import circuit
circuit.default_toolkit = jax_toolkit
from pycircuit.circuit.elements import SubCircuit, TLine

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')
tline = TLine('in', '0', 'out', '0', Z0=50.0, TD=1e-9)
tline.history = []
tline.iparv = type('IPARV', (), {'Z0': 50.0, 'TD': 1e-9})
tline.toolkit = jax_toolkit

G = tline.G([0]*6)
print(G)
