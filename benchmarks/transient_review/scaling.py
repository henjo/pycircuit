"""Dense LU vs SuperLU (factor / refactor-reuse) on MNA-like sparsity."""
import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

def mna_like(n, seed=0):
    """Tridiagonal-ish ladder plus a few long-range couplings: ~5 nnz/row."""
    rng = np.random.default_rng(seed)
    rows, cols, vals = [], [], []
    for i in range(n):
        rows.append(i); cols.append(i); vals.append(4.0 + rng.random())
        if i + 1 < n:
            for (r, c) in ((i, i+1), (i+1, i)):
                rows.append(r); cols.append(c); vals.append(-1.0 - 0.1*rng.random())
    # a few random couplings (feedback / controlled sources)
    for _ in range(n // 4):
        i = int(rng.integers(0, n)); j = int(rng.integers(0, n))
        rows.append(i); cols.append(j); vals.append(0.3)
    A = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsc()
    return A

def timeit(f, reps):
    t0 = time.perf_counter()
    for _ in range(reps):
        f()
    return (time.perf_counter() - t0) / reps

print('%6s %8s %9s %10s %10s %10s %10s %8s' %
      ('n', 'nnz', 'dens%', 'dense_s', 'splu_s', 'resolve_s', 'dense_MB', 'fill'))
for n in (137, 300, 600, 1200, 2500, 5000):
    A = mna_like(n)
    b = np.random.rand(n)
    Ad = A.toarray()
    reps = max(2, min(30, int(2e7 / (n**2))))
    t_dense = timeit(lambda: np.linalg.solve(Ad, b), reps)
    t_splu = timeit(lambda: spl.splu(A).solve(b), reps)
    lu = spl.splu(A)
    t_res = timeit(lambda: lu.solve(b), max(reps, 20))
    fill = (lu.L.nnz + lu.U.nnz) / A.nnz
    print('%6d %8d %9.3f %10.5f %10.5f %10.5f %10.1f %8.2f' %
          (n, A.nnz, 100.0*A.nnz/(n*n), t_dense, t_splu, t_res,
           n*n*8/1e6, fill))
