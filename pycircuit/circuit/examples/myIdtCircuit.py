from pycircuit.circuit import numeric, SubCircuit, gnd
from pycircuit.circuit.elements import VCVS, VCCS, C

c = SubCircuit()
iplus,iminus,int0,oplus,ominus = c.add_nodes('1', '2', '3', '4', '5')
c['vccs'] = VCCS(iplus,iminus,int0,gnd,gm=1.)
c['vcvs'] = VCVS(int0,gnd,oplus,ominus,g=-1.)
c['c'] = C(int0,gnd, c=1.)


