from pycircuit.circuit import Circuit, R, VS, gnd, circuit
from pycircuit.circuit.elements import SubCircuit

cir = SubCircuit()
cir.add_node('in')
cir.add_node('out')
cir['V1'] = VS('in', gnd, v=1.0)
print(cir.get_node_index(gnd))
print(cir.get_node_index('in'))
print(cir.get_node_index('out'))
