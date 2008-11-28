# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from scipy import *
from circuit import *
from analysis import *
from pylab import *

cir=SubCircuit()
net1 = cir.add_node("net1")
net2 = cir.add_node("net2")

cir['R1'] = R(net1, net2, 1e3)
cir['C1'] = C(net2, gnd, 1e-12)
cir['VS'] = VS(net1, gnd, 1.0)

ac = AC(cir)
ac.run(freqs=logspace(6,9))

res = ac.getResult()
vnet1 = res['net2'].db20()

print vnet1
vnet1.semilogx()
grid(True)
show()
