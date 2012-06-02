#This code finds the gain of a simple common-source amplifier stage, using symgolic toolkit.

from sympy import *
from pycircuit.circuit import *
from pycircuit.circuit import mos
c=SubCircuit(toolkit=symbolic)
inp=c.add_node('inp')
out=c.add_node('out')
vdd=c.add_node('vdd')
var('R1,gm1,s')
c['VDD']=VS(vdd,gnd,v=5,vac=0)
c['R1']=R(vdd,out,r=R1)
c['Vin']=VS(inp,gnd,v=1,vac=1)
c['M1']=mos.MOS(inp,out,gnd,gnd,gm=gm1,gds=0,gmb=0,toolkit=symbolic)
ac=AC(c)
res=ac.solve(s,complexfreq=True)
gain=sympy.simplify(res.v('out')/res.v('inp'))
print "The gain of the CS stage is:"
sympy.pprint(gain)

