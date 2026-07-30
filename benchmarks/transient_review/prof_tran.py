"""Profile the default Transient path on the numeric leapfrog replica."""
import cProfile, pstats, io, sys, time
import numpy as np

sys.path.insert(0, '/home/andreas/sources/pycircuit')
sys.path.insert(0, '/home/andreas/sources/pycircuit/benchmarks')

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.transient import Transient

F1 = 8e3
G3 = 1e-3


def build(stages=5):
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        node, amp, _ = build_leapfrog_network(cir, stages=stages,
                                              tones=((1e-3, F1),))
        nl = amp[0]['e1']
        cir['nl'] = BSource(nl, gnd, nl, gnd, i_func=lambda v: G3 * v ** 3)
        return cir, amp
    finally:
        circuit_module.default_toolkit = saved


cir, amp = build()
print('n unknowns =', cir.n, 'nodes =', len(cir.nodes), 'branches =', len(cir.branches))
print('n elements (top level) =', len(cir.elements))

# matrix sparsity at x=0
x0 = np.zeros(cir.n)
t0 = time.time()
G = cir.G(x0)
t_G = time.time() - t0
t0 = time.time()
C = cir.C(x0)
t_C = time.time() - t0
t0 = time.time()
i = cir.i(x0)
t_i = time.time() - t0
t0 = time.time()
q = cir.q(x0)
t_q = time.time() - t0
t0 = time.time()
u = cir.u(0.0, analysis='tran')
t_u = time.time() - t0

J = G + C / 1e-9
nnz = int(np.count_nonzero(J))
print('J shape %s nnz %d density %.3f%%' % (J.shape, nnz, 100.0 * nnz / J.size))
print('assembly times: G=%.4f C=%.4f i=%.4f q=%.4f u=%.4f  total=%.4f s'
      % (t_G, t_C, t_i, t_q, t_u, t_G + t_C + t_i + t_q + t_u))

# solve cost, dense
Jr = np.delete(np.delete(J, 0, 0), 0, 1)
b = np.random.rand(Jr.shape[0])
N = 50
t0 = time.time()
for _ in range(N):
    np.linalg.solve(Jr, b)
t_dense = (time.time() - t0) / N
print('dense np.linalg.solve: %.5f s' % t_dense)

import scipy.sparse as sp
import scipy.sparse.linalg as spl
Jsp = sp.csc_matrix(Jr)
t0 = time.time()
for _ in range(N):
    spl.spsolve(sp.csc_matrix(Jr), b)
t_spsolve = (time.time() - t0) / N
print('scipy spsolve (rebuild csc each time): %.5f s' % t_spsolve)

t0 = time.time()
for _ in range(N):
    lu = spl.splu(Jsp)
    lu.solve(b)
t_splu = (time.time() - t0) / N
print('scipy splu factor+solve (csc prebuilt): %.5f s' % t_splu)

lu = spl.splu(Jsp)
t0 = time.time()
for _ in range(N):
    lu.solve(b)
t_reuse = (time.time() - t0) / N
print('scipy splu solve only (reuse factors): %.5f s' % t_reuse)

# fill-in of LU
print('splu L nnz %d U nnz %d (total %d vs A nnz %d)'
      % (lu.L.nnz, lu.U.nnz, lu.L.nnz + lu.U.nnz, Jsp.nnz))

# ---- profile a real (short) transient
tran = Transient(cir, toolkit=numeric)
pr = cProfile.Profile()
t0 = time.time()
pr.enable()
res = tran.solve(refnode=gnd, tend=20e-6, timestep=1.0 / (F1 * 128))
pr.disable()
wall = time.time() - t0
nsteps = len(res.sweep_values)
print('\ntransient: %d steps in %.2f s = %.4f s/step' % (nsteps, wall, wall / nsteps))
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('cumulative').print_stats(45)
print(s.getvalue())
s = io.StringIO()
pstats.Stats(pr, stream=s).sort_stats('tottime').print_stats(35)
print(s.getvalue())
