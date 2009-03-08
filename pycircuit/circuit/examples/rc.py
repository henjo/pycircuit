# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from scipy import *
from pycircuit.circuit import *
from pylab import *
from pycircuit.post.functions import *

cir=SubCircuit()
net1 = cir.add_node("net1")
net2 = cir.add_node("net2")

cir['R1'] = R(net1, net2, r=1e3)
cir['C1'] = C(net2, gnd, c=1e-12)
cir['VS'] = VS(net1, gnd, vac=1.0)

ac = AC(cir)
res = ac.solve(freqs=logspace(6,9))

vnet1 = db20(res['net2'])

print vnet1
vnet1.semilogx()
grid(True)
show()
