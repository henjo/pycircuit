from pycircuit.sim import ngspice
from pycircuit.post import Waveform, db20
from pycircuit.post.testing import *
from pylab import plot, show

import numpy as np

#Single transistor netlist with model
simple_netlist = """
.model cmosn nmos (level=2 ld=0.265073u tox=418.0e-10
+ nsub=1.53142e+16 vto=0.844345 kp=4.15964e-05 gamma=0.863074
+ phi=0.6 uo=503.521 uexp=0.163917 ucrit=161166
+ delta=1e-06 vmax=55903.5 xj=0.400000u lambda=0.01
+ nfs=3.5934e+12 neff=1.001 nss=1e+12 tpg=1.000000
+ rsh=29.3 cgdo=2.18971e-10 cgso=2.18971e-10
+ cj=0.0003844 mj=0.488400 cjsw=5.272e-10 mjsw=0.300200 pb=0.700000)
MN 2 1 0 0 cmosn  l=0.80u w=300u
VDD 3 0 2.5
VIN 1 0 0.5
RL 3 2 2e2
"""

#Create a circuit from the netlist
cir = ngspice.Circuit(simple_netlist)

#Create a simulation
sim = ngspice.Simulation(cir, direct=False)
sim.command('') # flush

sim.command('op')
sim.command('') # flush

res = sim.command('show mn')
res = sim.command('') #flush

mn_d = dict([l.split() for l in res.split('\n')])

print res

#Set up DC-sweep for id and gm
# sim.command('print dc id(MN) gm(MN)')
# res_id_gm=sim.command('dc VIN 0 2.5 1m',parse_result=True)

# #Set up DC-sweep for gds
# sim.command('print dc gds(MN)')
# res_gds=sim.command('dc VDD 0 2.5 1m',parse_result=True)

# gmvec=res_id_gm['gm(MN)']
# idvec=res_id_gm['id(MN)']
# gdsvec=res_gds['gds(MN)']

# plot(gmvec/idvec)
# show()
