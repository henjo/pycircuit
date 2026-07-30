from pycircuit.circuit import circuit, R, VS, gnd
from pycircuit.circuit.elements import SubCircuit, TLine

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')
cir['V1'] = VS('in', gnd, v=1.0)
cir['T1'] = TLine('in', gnd, 'out', gnd, Z0=50.0, TD=1e-9)
cir['R1'] = R('out', gnd, r=50.0)

print("Nodes:", cir.nodes)
print("Branches:", cir.branches)
print("V1 indices:", cir.elementnodemap['V1'])
print("T1 indices:", cir.elementnodemap['T1'])
