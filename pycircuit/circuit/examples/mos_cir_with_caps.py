#This code finds the gain of a CS stage, now considering the intrinsic MOS capacitantes. The results is the same as in the book Design Analog CMOS Integrated Circuits, by Behzad Razavi - pg 174.

from sympy import *
from pycircuit.circuit import *
from pycircuit.circuit import mos
c=SubCircuit(toolkit=symbolic)
inp=c.add_node('inp')
inp1=c.add_node('inp1')
out=c.add_node('out')
vdd=c.add_node('vdd')
var('R_L,R_S,gm1,gmb1,ro1,Cgs1,Cgd1,Cdb1,s')
c['VDD']=VS(vdd,gnd,v=5,vac=0)
c['R_L']=R(vdd,out,r=R_L)
c['R_S']=R(inp,inp1,r=R_S)
c['Vin']=VS(inp,gnd,v=1,vac=1)
c['M1']=mos.MOS(inp1,out,gnd,gnd,gm=gm1,gds=0,gmb=0,Cgs=Cgs1,Cgd=Cgd1,Cdb=Cdb1,toolkit=symbolic)
ac=AC(c)
res=ac.solve(s,complexfreq=True)
gain=simplify(res.v('out')/res.v('inp'))
print "\nThe transfer function of the CS stage is:"
sympy.pprint(gain)
print "\nShowing the denominator as polynomial:"
sympy.pprint(denom(gain).as_poly(s))
