"""Does the matrix-free gate survive a SPARSE baseline?  YES -- MEASURED.

RESULT (2026-09-02, quiet box at load 1.7, single-threaded, 3 Newton
iterations each).  Wall seconds:

      m | DenseSolver         | SuperLUSolver       | matrix-free gain
        | dense      mfree    | dense      mfree    | vs dense  vs splu
    ----+---------------------+---------------------+------------------
    110 |  0.713     0.646    |  0.739     0.682    |  1.10x     1.08x
    242 |  2.166     1.596    |  2.045     1.531    |  1.36x     1.34x
    502 |  9.301     6.614    |  7.897     5.298    |  1.41x     1.49x
   1002 | 51.083    26.049    | 40.804    17.674    |  1.96x     2.31x

THE GATE HOLDS.  The worry that prompted this -- that every recorded
speedup had a dense baseline on both sides, so a sparse one might erase the
win -- is NOT borne out.  Matrix-free gains the SAME or MORE against
SuperLU than against dense LAPACK at every size, and m~250 remains where it
crosses into being worth asking for.

⚠ AND THE TWO ARE COMPLEMENTARY, not alternatives.  `splu/mfree` is fastest
at every m >= 242, and at m=1002 it is 17.674 s against the 51.083 s of the
dense/dense path this analysis shipped with -- 2.89x, most of which neither
change delivers alone.  A sparse solver alone buys 1.25x there;
matrix-free alone buys 1.96x.

⚠ SuperLU is a LOSS at m=110 (0.739 against 0.713): sparse bookkeeping
costs more than it saves on a small matrix, the usual crossover.

⚠ AND `rc_ladder` IS TRIDIAGONAL -- about the most favourable structure a
sparse solver can be handed, so this OVERSTATES the sparse baseline\'s
advantage.  That is the direction that makes the conclusion safe: a sparse
baseline did not erode matrix-free\'s gain even where sparsity is most
flattered.  A denser real circuit moves `splu/dense` back toward
`dense/dense` and can only widen the matrix-free margin.

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

⚠ THE FIRST ATTEMPT WAS THROWN AWAY, and why is worth keeping.  It read
load average 31-43 with a SECOND session benchmarking the same machine, and
its readings disagreed with earlier ones by 25-30% on the dense/dense pair
alone.  Under this campaign's own rule that is not a number.  Check
`ps -eo pid,pcpu,args --sort=-pcpu` and `uptime` before quoting anything.

The questions it was written to settle, and what it found:

  1. Does matrix-free still win once the baseline is sparse?  YES, by the
     same margin or more -- both sides get sparse solves and scale together.
  2. Does the `m~250` crossover move?  NO.
  3. Is "use a sparse solver" ever the better recommendation on its own?
     NO -- `splu/dense` never wins a row.  Above m~250 the answer is BOTH,
     below it neither.

⚠ AND THE LADDER FLATTERS SPARSITY.  `rc_ladder` is tridiagonal, about the
most favourable structure a sparse solver can be handed; a real circuit's
Jacobian is sparse but not this sparse.  Whatever this measures is an upper
bound on the sparse baseline's advantage, and the conclusion should say so.
"""
import builtins
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
    ## ⚠ FLUSHED, because this one runs for minutes and its output goes to a
    ## file more often than to a terminal.  Without it the rows sit in the
    ## buffer and the run looks hung, which is how the first attempt got
    ## killed and restarted for no reason.
    import functools
    print = functools.partial(builtins.print, flush=True)
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
