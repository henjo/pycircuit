from pycircuit.circuit.elements import SubCircuit, TLine
tline = TLine('p1', 'm1', 'p2', 'm2')
tline.history = []
import pycircuit.circuit._jaxtoolkit as toolkit
tline.toolkit = toolkit
from pycircuit.circuit import ParameterValue
tline.iparv = ParameterValue({'Z0': 50.0, 'TD': 1e-9})
G = tline.G([0,0,0,0,0,0])
print(G)
