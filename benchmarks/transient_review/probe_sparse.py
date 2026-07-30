"""Does the sparse_numeric toolkit actually stay sparse through Transient/DC?"""
import sys, time, traceback
import numpy as np
sys.path.insert(0, '/home/andreas/sources/pycircuit')

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.toolkit import sparse_numeric
from pycircuit.circuit.elements import R, C, IS, VS, Diode
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.dcanalysis import DC


def rc_chain(N, toolkit):
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = toolkit
    try:
        cir = SubCircuit(toolkit=toolkit)
        cir['IS'] = IS(gnd, 'n0', i=1e-3)
        for k in range(N):
            a = 'n%d' % k
            b = 'n%d' % (k + 1)
            cir['R%d' % k] = R(a, b, r=1e3)
            cir['C%d' % k] = C(a, gnd, c=1e-9)
        cir['Rend'] = R('n%d' % N, gnd, r=1e3)
        return cir
    finally:
        circuit_module.default_toolkit = saved


print("=== what type does G() return under sparse_numeric? ===")
cir = rc_chain(4, sparse_numeric)
x0 = np.zeros(cir.n)
G = cir.G(x0)
print('type(G) =', type(G))
print('has build_sparse:', hasattr(sparse_numeric, 'build_sparse'))

print("\n=== DC with sparse_numeric on RC chain ===")
try:
    r = DC(rc_chain(4, sparse_numeric), toolkit=sparse_numeric).solve()
    print('DC ok, v(n0) =', r.v('n0'))
except Exception:
    traceback.print_exc()

print("\n=== Transient with sparse_numeric on RC chain ===")
try:
    r = Transient(rc_chain(4, sparse_numeric), toolkit=sparse_numeric).solve(
        tend=1e-6, timestep=1e-7)
    print('TRAN ok, last v(n0) =', r.v('n0')[-1])
except Exception:
    traceback.print_exc()

print("\n=== is the reduced J passed to linearsolver dense or sparse? ===")
import pycircuit.circuit._sparse_numeric as sn
orig = sn.linearsolver
seen = []
def spy(*a, **k):
    seen.append((type(a[0]).__name__, getattr(a[0], 'shape', None)))
    return orig(*a, **k)
sn.linearsolver = spy
try:
    DC(rc_chain(4, sparse_numeric), toolkit=sparse_numeric).solve()
except Exception as e:
    print('err', e)
sn.linearsolver = orig
print('linearsolver saw:', seen[:5])
