from pycircuit.circuit.elements import VS, R, SubCircuit, gnd
from pycircuit.circuit import DC
from pycircuit.circuit import symbolic
import sympy

c = SubCircuit(toolkit=symbolic)
V0 = sympy.Symbol('V0')
R1 = sympy.Symbol('R1')
c['V0'] = VS(1, gnd, v=V0, toolkit=symbolic)
c['R1'] = R(1, gnd, r=R1)

dc = DC(c, toolkit=symbolic)
res = dc.solve()
print("Success:", res)
