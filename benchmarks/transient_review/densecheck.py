import time, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
n=136
rng=np.random.default_rng(0)
A=np.eye(n)*4+rng.standard_normal((n,n))*0.01
b=rng.standard_normal(n)
def t(f,r=500):
    f(); t0=time.perf_counter()
    for _ in range(r): f()
    return (time.perf_counter()-t0)/r
print('dense solve n=136: %.6f s' % t(lambda: np.linalg.solve(A,b)))
print('dense lu_factor+solve: ', end='')
import scipy.linalg as sl
print('%.6f s' % t(lambda: sl.lu_solve(sl.lu_factor(A), b)))
