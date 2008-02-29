from scipy import *
from circuit import *
from analysis import *

cir=SubCircuit()
net1 = cir.addNode("net1")
net2 = cir.addNode("net2")

cir['R1'] = R(net1, net2, 1e3)
cir['C1'] = C(net2, gnd, 1e-12)
cir['VS'] = VS(net1, gnd, 1.0)

ac = AC(cir)
ac.run(logspace(6,9))

res = ac.getResult()
vnet1 = db20(res.getSignal('net1'))
print vnet1
