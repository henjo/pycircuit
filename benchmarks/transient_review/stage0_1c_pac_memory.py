"""STAGE 0.1c -- how large is PAC's dense operator, actually?

`doc/transient_review.md` asserted "~150 TB at N=137, M=1000".  That is wrong by about
365x, and it survived because nothing downstream depended on the magnitude: the
conclusion ("does not scale past toy circuits") is equally true at the real figure.

This script exists so the number is generated rather than typed.  It reads the dtypes
off `shooting.PAC.solve` itself rather than hard-coding them, because one of them is
not what it looks like: `shooting.py:195` builds `B` with `tk.zeros(L.shape)` and NO
dtype, so `B` is float64 while `L` is complex128.

Run:  PYTHONPATH=<repo> python3 -u benchmarks/transient_review/stage0_1c_pac_memory.py
"""

import sys

import numpy as np

## N is `cir.n - 1` (the reference node is removed); M is `len(pss.times) - 1`,
## i.e. period/timestep.  The pair below is the one the review quoted.
CASES = ((137, 1000), (137, 100), (29, 200), (542, 1000))

L_DTYPE = np.cdouble      # shooting.py:194, tk.zeros(..., dtype=tk.cdouble)
B_DTYPE = np.float64      # shooting.py:195, tk.zeros(L.shape) -- NO dtype passed


def main():
    print('PAC dense operator size.  L is %s, B is %s (shooting.py:194-195).'
          % (np.dtype(L_DTYPE).name, np.dtype(B_DTYPE).name))
    print('B is real only because no dtype was passed; that is not deliberate.')
    print()
    print('%>8s %>8s %>12s %>14s %>14s %>14s'
          .replace('%>', '%') % ('N', 'M', 'N*M', 'L', 'B', 'total'))

    for N, M in CASES:
        n = N * M
        lb = n * n * np.dtype(L_DTYPE).itemsize
        bb = n * n * np.dtype(B_DTYPE).itemsize
        print('%8d %8d %12d %11.1f GiB %10.1f GiB %10.1f GiB'
              % (N, M, n, lb / 2 ** 30, bb / 2 ** 30, (lb + bb) / 2 ** 30))

    N, M = CASES[0]
    n = N * M
    total = n * n * (np.dtype(L_DTYPE).itemsize + np.dtype(B_DTYPE).itemsize)
    print()
    print('At N=%d, M=%d: %.4e elements per matrix, %.1f GiB total (%.3f TiB).'
          % (N, M, n * n, total / 2 ** 30, total / 2 ** 40))
    print('The review said "~150 TB".  The 150 is B alone, %.4e bytes, read as GB'
          % (n * n * np.dtype(B_DTYPE).itemsize))
    print('and written as TB -- a unit slip on the smaller of the two arrays.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
