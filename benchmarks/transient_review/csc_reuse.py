"""Cost of the four candidate linear-solve strategies, on the real leapfrog J.

(a) dense np.linalg.solve                     -- what pycircuit does now
(b) csc_matrix(dense) + spsolve               -- what sparse_numeric does now
(c) fixed-pattern CSC (data write) + splu     -- symbolic pattern reused
(d) fixed-pattern CSC + lu.solve only         -- factors reused (bypass/latency)
"""
import sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, '/home/andreas/sources/pycircuit')

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource

def build(stages):
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        _, amp, _ = build_leapfrog_network(cir, stages=stages, tones=((1e-3, 8e3),))
        cir['nl'] = BSource(amp[0]['e1'], gnd, amp[0]['e1'], gnd,
                            i_func=lambda v: 1e-3 * v ** 3)
        return cir
    finally:
        circuit_module.default_toolkit = saved

def timeit(f, reps=200):
    f()
    t0 = time.perf_counter()
    for _ in range(reps):
        f()
    return (time.perf_counter() - t0) / reps

print('%7s %6s %7s %8s %11s %11s %11s %11s %7s' %
      ('stages', 'n', 'nnz', 'dens%', '(a)dense', '(b)spsolve', '(c)splu', '(d)resolve', 'fill'))
for stages in (1, 5, 12, 20):
    print("building %d..." % stages, flush=True)
    cir = build(stages)
    n = cir.n
    x0 = np.zeros(n)
    J = cir.G(x0) + cir.C(x0) / 1e-9
    Jr = np.delete(np.delete(J, 0, 0), 0, 1)
    b = np.random.rand(n - 1)
    A = sp.csc_matrix(Jr)
    nnz = A.nnz
    # fixed pattern: keep indices/indptr, only overwrite .data
    pat_i, pat_p = A.indices.copy(), A.indptr.copy()
    data = A.data.copy()
    def c():
        M = sp.csc_matrix((data, pat_i, pat_p), shape=A.shape)
        lu = spl.splu(M)
        return lu.solve(b)
    lu0 = spl.splu(A)
    reps = 200 if n < 800 else 20
    ta = timeit(lambda: np.linalg.solve(Jr, b), reps)
    tb = timeit(lambda: spl.spsolve(sp.csc_matrix(Jr), b), reps)
    tc = timeit(c, reps)
    td = timeit(lambda: lu0.solve(b), max(reps, 200))
    print('%7d %6d %7d %8.3f %11.6f %11.6f %11.6f %11.6f %7.2f' %
          (stages, n, nnz, 100.0*nnz/(n-1)**2, ta, tb, tc, td,
           (lu0.L.nnz + lu0.U.nnz)/nnz))
