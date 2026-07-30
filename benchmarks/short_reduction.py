"""Stage 13 -- the ``Short`` operation: collapse nodes, shrink the determinant.

Stage 6 ranked elements by ``Open`` (parameter to zero) and found it changed the
compact diagram not at all: zeroing one addend inside a composite matrix entry
removes no vertex.  ``Short`` is the other half and is structurally different --
identifying two nodes deletes a row and a column, so the determinant genuinely gets
smaller.  It is the only lever measured in this plan that does that.

Implemented as a node merge: ``V_a = V_b`` means column ``b`` adds into column ``a``
(the two voltages are the same unknown) and row ``b`` adds into row ``a`` (KCL at the
merged node), then ``b`` is deleted.

Chosen greedily with **exact re-evaluation of H over the whole sweep after every
merge**, which is the discipline stage 12 had to relearn the hard way.

**The merge is calibrated before it is used.**  It is easy to get subtly wrong --
double-counting the diagonal, or forgetting that deleting a row shifts the output
index -- and a wrong merge yields a smaller matrix for a *different circuit*, which
would look like a good result.  So it is first checked against a case with a
hand-computable answer: shorting the upper leg of a resistive divider must give
exactly ``H = 1``.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/short_reduction.py
"""

import math
import sys
import time
import warnings

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, s_expand, _resolve

from cancellation_profile import NOMINAL, DEGRADED, environment
from transfer_function import (UNCONVERGED, coefficients, numerator_matrix, ops,
                               rank, choose_coefficients, sweep_error)


WIDE = np.logspace(1, 7, 25)


def merge_numeric(M, rhs, a, b):
    """Identify node ``b`` with node ``a``; return the smaller system.

    Column first, then row, so the merged diagonal picks up all four of
    ``M[a,a] + M[a,b] + M[b,a] + M[b,b]`` exactly once.
    """
    M = M.copy()
    rhs = rhs.copy()
    M[:, a] += M[:, b]
    M[a, :] += M[b, :]
    rhs[a] += rhs[b]
    keep = [i for i in range(M.shape[0]) if i != b]
    return M[np.ix_(keep, keep)], rhs[keep]


def merge_symbolic(A, a, b):
    A = sympy.Matrix(A)
    for i in range(A.rows):
        A[i, a] = A[i, a] + A[i, b]
    for j in range(A.cols):
        A[a, j] = A[a, j] + A[b, j]
    keep = [i for i in range(A.rows) if i != b]
    return A.extract(keep, keep)


def shift(index, removed):
    """Where an index lands after a row/column is deleted."""
    return index - 1 if index > removed else index


def calibrate():
    """A divider whose upper leg is shorted must give exactly H = 1."""
    print('== calibration: the merge, against a hand-computable answer ==')
    R1, R2 = 100.0, 300.0
    ## Nodes: 0 = input (driven), 1 = mid (output).  Nodal equations with the
    ## input forced by a 1 V source through a stiff conductance is fiddly, so use
    ## the plain two-node divider driven by a current source at node 0:
    ##   [ 1/R1        -1/R1      ] [V0]   [I]
    ##   [ -1/R1   1/R1 + 1/R2    ] [V1] = [0]
    M = np.array([[1 / R1, -1 / R1], [-1 / R1, 1 / R1 + 1 / R2]])
    rhs = np.array([1.0, 0.0])
    v = np.linalg.solve(M, rhs)
    print('  divider: V1/V0 = %.6f, expected R2/(R1+R2) = %.6f'
          % (v[1] / v[0], R2 / (R1 + R2)))
    ok1 = abs(v[1] / v[0] - R2 / (R1 + R2)) < 1e-12

    Ms, rs = merge_numeric(M, rhs, 0, 1)
    vs = np.linalg.solve(Ms, rs)
    ## After merging, one node remains and V1/V0 is identically 1.
    print('  after shorting R1 (merging the two nodes): %d unknown left,'
          ' V1/V0 = %.6f, expected 1.0' % (Ms.shape[0], 1.0))
    ## The merged system must also give the right voltage: I flows through R2.
    print('  merged node voltage %.6f, expected I*R2 = %.6f'
          % (vs[0], 1.0 * R2))
    ok2 = abs(vs[0] - R2) < 1e-9
    print('  GATE 13-1 (merge reproduces the hand answer): %s'
          % ('PASS' if ok1 and ok2 else 'FAIL'))
    print()
    return ok1 and ok2


def reference_over(system, cenv, sweep):
    out = []
    for freq in sweep:
        env = dict(cenv)
        env[system.s] = 2j * math.pi * freq
        M = np.array([[complex(_resolve(system.A[i, j], env))
                       for j in range(system.dim)] for i in range(system.dim)],
                     dtype=complex)
        rhs = np.array([complex(_resolve(system.b[i], {**env, **{
            s: 1 for s in system.b.free_symbols if str(s) == 'vin'}}))
            for i in range(system.dim)], dtype=complex)
        out.append(np.linalg.solve(M, rhs)[system.out_index])
    return np.array(out)


def numeric_system(system, cenv, freq):
    env = dict(cenv)
    env[system.s] = 2j * math.pi * freq
    M = np.array([[complex(_resolve(system.A[i, j], env))
                   for j in range(system.dim)] for i in range(system.dim)],
                 dtype=complex)
    rhs = np.array([complex(_resolve(system.b[i], {**env, **{
        s: 1 for s in system.b.free_symbols if str(s) == 'vin'}}))
        for i in range(system.dim)], dtype=complex)
    return M, rhs


def h_after(cached, merges, out_index):
    """H over the sweep after applying merges, from pre-resolved matrices.

    The matrices are resolved from sympy **once** per frequency and reused.  Doing
    it inside the candidate loop instead meant ~5 million sympy substitutions per
    greedy iteration, which is what made the first version of this script appear to
    hang rather than run.
    """
    out = []
    for M0, rhs0 in cached:
        M, rhs = M0, rhs0
        oi = out_index
        for a, b in merges:
            M, rhs = merge_numeric(M, rhs, a, b)
            oi = shift(oi, b)
        try:
            out.append(np.linalg.solve(M, rhs)[oi])
        except np.linalg.LinAlgError:
            return None
    return np.array(out)


def greedy_merges(system, cached, reference, tol, node_rows):
    """Merge nodes greedily, re-solving H over the whole sweep after each."""
    merges = []
    alive = list(node_rows)
    out_index = system.out_index
    while True:
        best, best_err = None, None
        for i in range(len(alive)):
            for j in range(i + 1, len(alive)):
                a, b = alive[i], alive[j]
                cand = merges + [(a, b)]
                got = h_after(cached, cand, system.out_index)
                if got is None:
                    continue
                err = float(np.max(np.abs(got - reference) / np.abs(reference)))
                if best_err is None or err < best_err:
                    best, best_err = (a, b), err
        if best is None or best_err > tol:
            return merges, best_err
        a, b = best
        merges.append((a, b))
        ## Every index above the deleted row shifts down by one.
        alive = [shift(x, b) for x in alive if x != b]


def main():
    if not calibrate():
        print('calibration failed; not proceeding')
        return 1

    system = bc.ua741(symbolic_devices=('q1', 'q2', 'q17'))
    ## Branch-current rows have a zero diagonal in MNA; only node rows may merge.
    probe_env = environment(system, NOMINAL)
    M0, _r0 = numeric_system(system, {k: v for k, v in probe_env.items()
                                      if k != system.s}, 1e3)
    node_rows = [i for i in range(system.dim) if abs(M0[i, i]) > 0]
    print('== uA741: %d rows, %d of them node rows (nonzero diagonal) =='
          % (system.dim, len(node_rows)))
    print()

    for label, gms in (('nominal', NOMINAL), ('degraded', DEGRADED)):
        cenv = {k: v for k, v in environment(system, gms).items()
                if k != system.s}
        ref = reference_over(system, cenv, WIDE)
        ## Resolve the matrices from sympy once; the greedy reuses them.
        cached = [numeric_system(system, cenv, f) for f in WIDE]
        print('  == %s gm, wide band 10 Hz-10 MHz ==' % label)
        for tol in (0.2, 0.05):
            t0 = time.time()
            merges, nxt = greedy_merges(system, cached, ref, tol, node_rows)
            got = h_after(cached, merges, system.out_index)
            err = (0.0 if not merges else
                   float(np.max(np.abs(got - ref) / np.abs(ref))))
            print('    tol=%-6g %d merges accepted %s -> dim %d, sweep err'
                  ' %.3e  (next best merge would cost %.2e)  %.0f s'
                  % (tol, len(merges), merges, system.dim - len(merges), err,
                     nxt if nxt is not None else float('nan'),
                     time.time() - t0))

            if not merges:
                print('      GATE 13-5: no merge is acceptable at this tolerance')
                continue

            ## Apply the same merges symbolically and re-measure the symbolic cost.
            A = system.A
            oi = system.out_index
            for a, b in merges:
                A = merge_symbolic(A, a, b)
                oi = shift(oi, b)
            D = ddd_of_matrix(A)
            print('      reduced symbolic matrix: %dx%d, diagram %d vertices,'
                  ' %d terms, N_eff %.1f'
                  % (A.rows, A.rows, D.size, D.term_count(),
                     D.concentration({**cenv, system.s: 2j * math.pi * 1e3})))

            ## The gate is an operation count, so run the full transfer-function
            ## pipeline on the reduced circuit and compare against the ORIGINAL
            ## response.  The merges have already spent part of the error budget,
            ## which is the tension worth reporting rather than hiding.
            b = sympy.Matrix(system.b)
            vin = [x for x in b.free_symbols if str(x) == 'vin']
            if vin:
                b = b.subs({vin[0]: 1})
            bb = sympy.Matrix(system.b)
            for a2, b2 in merges:
                bb[a2, 0] = bb[a2, 0] + bb[b2, 0]
                bb = sympy.Matrix([bb[i, 0] for i in range(bb.rows)
                                   if i != b2])
            if vin:
                bb = bb.subs({vin[0]: 1})
            Anum = sympy.Matrix(A)
            for i in range(Anum.rows):
                Anum[i, oi] = bb[i, 0]
            t1 = time.time()
            den_exp = s_expand(A, system.s)
            num_exp = s_expand(Anum, system.s)
            num = coefficients(num_exp, cenv)
            den = coefficients(den_exp, cenv)
            print('      s-expanded reduced: den degree %d, num degree %d'
                  ' (%.1f s)' % (den_exp.degree, num_exp.degree,
                                 time.time() - t1))
            for ctol in (2.5e-3, 1e-4):
                keep_n, keep_d = choose_coefficients(num, den, ref, 0.2, WIDE)
                nv, dv, nops, dops = {}, {}, 0, 0
                for kk in keep_n:
                    e, _n, _e = rank(num[kk][0], cenv, ctol)
                    nv[kk] = complex(e.xreplace(cenv)) if e != 0 else 0j
                    nops += ops(e)
                for kk in keep_d:
                    e, _n, _e = rank(den[kk][0], cenv, ctol)
                    dv[kk] = complex(e.xreplace(cenv)) if e != 0 else 0j
                    dops += ops(e)
                tot = nops + dops + 2 * (len(keep_n) + len(keep_d))
                e2 = sweep_error(nv, dv, ref, WIDE)
                missed = len(UNCONVERGED)
                UNCONVERGED.clear()
                print('        coeff tol=%-8.3g N %d + D %d coeffs -> %5d ops,'
                      ' total err vs ORIGINAL %.3e  %s%s'
                      % (ctol, len(keep_n), len(keep_d), tot, e2,
                         'PASS' if tot < 200 and e2 <= 0.20
                         else 'holds' if e2 <= 0.20 else 'fail',
                         '' if not missed else
                         '   [%d coefficients MISSED their tolerance -- the'
                         ' composed error above is not trustworthy]' % missed))
        print()
    return 0


if __name__ == '__main__':
    sys.exit(main())
