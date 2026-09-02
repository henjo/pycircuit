"""Does the matrix-free gate survive a SPARSE baseline?  UNMEASURED.

RAISED BY EXTERNAL REVIEW (2026-09-02).  Every recorded matrix-free speedup
-- 1.36x/1.51x/2.13x on the driven solved-history path and the rest -- was
measured with BOTH sides on dense LAPACK, because `_traverse_factored` and
`_step_sensitivity` reached past the caller's `linearsolver` to
`scipy.linalg` and `toolkit.linearsolver` respectively.  A circuit Jacobian
at m=1002 is very sparse and a dense LU there is ~3e8 flops per step, `N`
times over, so that baseline is not what a real simulator pays.

BOTH PATHS NOW GO THROUGH THE CALLER'S SOLVER (`PSS._factorise` and
`_step_sensitivity`), so the comparison can be run twice and the m~250 gate
re-asked.  This driver does that.

⚠ AND IT HAS NOT BEEN RUN ON A QUIET BOX, so nothing here is recorded as a
result.  The attempt on 2026-09-02 read load average 31-43 with a SECOND
session benchmarking the same machine, which is exactly the condition the
campaign's own record says makes a wall-clock ratio worthless -- and the
readings taken then disagreed with earlier ones by 25-30% on the dense/dense
pair alone.  They are not quoted here in either direction.

Re-run this on a quiet box.  The questions it settles, in order:

  1. Does matrix-free still win once the baseline is sparse?  Its own
     per-step solves get sparse too, so both sides move and the sign is not
     obvious from either side alone.
  2. Where is the crossover, if `m~250` moves?  The gate was set from a
     propagation share measured entirely on dense solves.
  3. Is `splu/dense` ever simply the best answer -- i.e. is "use a sparse
     solver" a better recommendation than "go matrix-free" for some range?

⚠ AND THE LADDER FLATTERS SPARSITY.  `rc_ladder` is tridiagonal, about the
most favourable structure a sparse solver can be handed; a real circuit's
Jacobian is sparse but not this sparse.  Whatever this measures is an upper
bound on the sparse baseline's advantage, and the conclusion should say so.
"""
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pss_matrix_free_ceiling import rc_ladder, NPTS

from pycircuit import circuit
from pycircuit.circuit import gnd
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit.linearsolver import DenseSolver, SuperLUSolver

ITERS = 3
SIZES = (108, 240, 500, 1000)


def main():
    circuit.default_toolkit = circuit.numeric
    print('driven solved-history (gear), %d Newton iterations.' % ITERS)
    print('⚠ Run this SINGLE-THREADED (OMP_NUM_THREADS=1) on a QUIET box, '
          'and check the load average before quoting anything.\n')
    print('%6s | %-21s | %-21s | %s'
          % ('m', 'DenseSolver', 'SuperLUSolver', 'verdict'))
    print('%6s | %-10s %-10s | %-10s %-10s |'
          % ('', 'dense', 'mfree', 'dense', 'mfree'))
    print('-' * 88)
    for sections in SIZES:
        row = []
        for solv in (DenseSolver, SuperLUSolver):
            ts = []
            for mf in (False, True):
                pss = PSS(rc_ladder(sections), method='gear', reltol=1e-6,
                          linearsolver=solv())
                t0 = time.perf_counter()
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    pss.solve(refnode=gnd, period=1e-3,
                              timestep=1e-3 / NPTS,
                              maxiterations=ITERS, matrix_free=mf)
                ts.append(time.perf_counter() - t0)
            row.append(ts)
        (dd, dm), (sd, sm) = row
        names = ['dense/dense', 'dense/mfree', 'splu/dense', 'splu/mfree']
        best = int(np.argmin([dd, dm, sd, sm]))
        print('%6d | %10.3f %10.3f | %10.3f %10.3f | mf %4.2fx dense, '
              '%4.2fx splu; best %s'
              % (sections + 2, dd, dm, sd, sm, dd / dm, sd / sm, names[best]))


if __name__ == '__main__':
    main()
