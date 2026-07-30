"""Stage 6 of `doc/cancellation_ranking_plan.md` -- rank elements, not terms.

Both literature sweeps (`doc/cancellation_ranking_conclusions.md` §14c, §15d)
converge on one lesson: **stop accounting error per term; evaluate the candidate
answer.**  Hu, Shi, Tai & Lee (ISCAS 2015) rank circuit *elements* by the exact
error their removal causes and get expressions matching textbook hand-derived
formulas; Guerra et al. (DATE 2000) score candidate approximations against a
global numerical reference for the same reason.

`kappa` cannot enter such a ranking -- nothing ever sums term magnitudes -- so
the wall the rest of this plan has been climbing is not in the path.  It is also
*simplification before generation*: fewer symbols means a smaller diagram for
everything downstream, which is the reason it is measured here rather than only
as a circuit-reduction tool.

Two deliberate limitations, stated rather than skipped:

* their Short/Open pair are GPDD **edge** operations; on an MNA matrix "Open"
  (parameter to zero) is immediate while "Short" (parameter to infinity, i.e.
  node merging) is **not implemented**.  This is the weaker half of their method.
* the error is measured over a **sampled** frequency sweep.  The literature is
  unanimous that this guarantees nothing between samples
  (Rodriguez-Garcia et al., DATE 1999: *"the error specs are met at the sampling
  frequencies, but exceeded at intermediate ones"*).

Selection is greedy with **exact joint re-evaluation** after every accepted
removal.  Individual removal errors are never added -- §11 of the reasoning
document measured what that does.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/element_ranking.py
"""

import math
import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.ddd import ddd_of_matrix, _resolve


DEVICES = tuple('q%d' % k for k in range(1, 25))
## Decades either side of the 1 kHz the rest of this work used, so a symbol that
## only matters out of band cannot be dropped for free.
FREQS = np.logspace(1, 6, 21)


def numeric_matrices(system, symbols):
    """Split A into a constant part and one matrix per symbol.

    Every entry of the uA741's matrix is linear in each ``gm`` (each device
    contributes its transconductance to a few entries), so
    ``A = A0 + sum_i gm_i * Ai`` exactly.  Splitting once means an evaluation is
    a few array adds instead of a sympy substitution, which is what makes an
    O(n^2) greedy search cheap.  Verified against direct substitution below.
    """
    A = system.A
    n = A.rows
    s = system.s
    base_env = {}
    for sym in A.free_symbols:
        if sym is s or sym in symbols:
            continue
        base_env[sym] = complex(system.params[sym])

    def at(env):
        return np.array([[complex(_resolve(A[i, j], env)) for j in range(n)]
                         for i in range(n)], dtype=complex)

    out = {}
    for freq in FREQS:
        env = dict(base_env)
        env[s] = 2j * math.pi * freq
        for sym in symbols:
            env[sym] = 0.0
        A0 = at(env)
        mats = {}
        for sym in symbols:
            e = dict(env)
            e[sym] = 1.0
            mats[sym] = at(e) - A0
        out[freq] = (A0, mats)
    return out, base_env


def transfer(split, freq, values, out_index, in_index, rhs):
    """Transfer function at one frequency with the given symbol values."""
    A0, mats = split[freq]
    A = A0.copy()
    for sym, v in values.items():
        if v:
            A = A + v * mats[sym]
    try:
        x = np.linalg.solve(A, rhs)
    except np.linalg.LinAlgError:
        return None
    return x[out_index]


def response(split, values, out_index, in_index, rhs):
    out = []
    for freq in FREQS:
        v = transfer(split, freq, values, out_index, in_index, rhs)
        if v is None:
            return None
        out.append(v)
    return np.array(out)


def worst_error(ref, got):
    """Worst relative error over the sweep; inf if the circuit went singular."""
    if got is None:
        return math.inf
    scale = np.abs(ref)
    ok = scale > 0
    if not ok.any():
        return math.inf
    return float(np.max(np.abs(got[ok] - ref[ok]) / scale[ok]))


def greedy(split, nominal, ref, out_index, in_index, rhs, tol):
    """Drop symbols one at a time, re-evaluating the joint error exactly.

    Returns ``(order, errors)``: the removal order and the *joint* error after
    each removal.  Nothing is ever summed.
    """
    values = dict(nominal)
    alive = [s for s in nominal]
    order, errors = [], []
    while alive:
        best, best_err = None, None
        for sym in alive:
            trial = dict(values)
            trial[sym] = 0.0
            err = worst_error(ref, response(split, trial, out_index, in_index,
                                            rhs))
            if best_err is None or err < best_err:
                best, best_err = sym, err
        if best_err > tol:
            break
        values[best] = 0.0
        alive.remove(best)
        order.append(best)
        errors.append(best_err)
    return order, errors, values


def main():
    t0 = time.time()
    system = bc.ua741(symbolic_devices=DEVICES)
    symbols = sorted((s for s in system.A.free_symbols
                      if str(s).startswith('gm_')), key=str)
    print('uA741 with every transistor gm symbolic: %dx%d, %d gm symbols'
          % (system.dim, system.dim, len(symbols)))
    print('frequency sweep: %d points, %.3g Hz to %.3g Hz (SAMPLED -- nothing is'
          ' guaranteed between samples)' % (len(FREQS), FREQS[0], FREQS[-1]))

    split, base_env = numeric_matrices(system, set(symbols))
    nominal = {s: complex(system.params[s]) for s in symbols}
    ## The right-hand side carries `vin`, which never appears in A, so it is
    ## absent from base_env -- resolve it from the fixture's own parameters.
    rhs_env = dict(base_env)
    rhs_env.update({k: 0.0 for k in symbols})
    for sym in system.b.free_symbols:
        if sym not in rhs_env:
            rhs_env[sym] = complex(system.params[sym])
    rhs = np.array([complex(_resolve(system.b[i], rhs_env))
                    for i in range(system.dim)], dtype=complex)
    print('split A = A0 + sum gm_i*Ai and prepared %d frequencies in %.1f s'
          % (len(FREQS), time.time() - t0))

    ## Check the split reproduces sympy's own substitution, so a wrong
    ## factorisation cannot masquerade as a circuit result.
    probe = FREQS[len(FREQS) // 2]
    env = dict(base_env)
    env[system.s] = 2j * math.pi * probe
    env.update(nominal)
    direct = np.array([[complex(_resolve(system.A[i, j], env))
                        for j in range(system.dim)]
                       for i in range(system.dim)], dtype=complex)
    A0, mats = split[probe]
    rebuilt = A0 + sum(nominal[s] * mats[s] for s in symbols)
    print('split check vs direct substitution: max abs diff %.2e'
          % np.max(np.abs(direct - rebuilt)))

    ref = response(split, nominal, system.out_index, system.in_index, rhs)
    print('reference |H| at %g Hz = %.6e' % (probe, abs(ref[len(FREQS) // 2])))
    print()

    for tol in (0.01, 0.05):
        t0 = time.time()
        order, errors, values = greedy(split, nominal, ref, system.out_index,
                                       system.in_index, rhs, tol)
        secs = time.time() - t0
        survivors = sorted((str(s) for s, v in values.items() if v), key=str)
        print('== joint error tolerance %.0f%% ==' % (100 * tol))
        print('  dropped %d of %d symbols in %.1f s (%d exact solves)'
              % (len(order), len(symbols), secs, len(FREQS) *
                 sum(range(len(symbols) - len(order) + 1, len(symbols) + 1))))
        print('  removal order, with the JOINT error after each step:')
        for sym, err in zip(order, errors):
            print('    %-12s %.3e' % (str(sym), err))
        print('  monotone in the greedy order: %s'
              % all(b >= a for a, b in zip(errors, errors[1:])))
        print('  SURVIVORS (%d): %s' % (len(survivors), ', '.join(survivors)))

        ## Now the point of doing this here: does pruning help the term ranking?
        keep = {s for s in symbols if values[s]}
        sub = {s: sympy.Integer(0) for s in symbols if s not in keep}
        reduced = system.A.subs(sub)
        env2 = dict(base_env)
        env2[system.s] = 2j * math.pi * 1e3
        for s in keep:
            env2[s] = nominal[s]
        full_env = dict(env2)
        for s in symbols:
            full_env.setdefault(s, nominal[s])

        t1 = time.time()
        D_full = ddd_of_matrix(system.A)
        D_red = ddd_of_matrix(reduced)
        k_full = D_full.cancellation(full_env)
        k_red = D_red.cancellation(env2)
        print('  diagram of the FULL circuit:      %d vertices, %d terms,'
              ' kappa %.3e' % (D_full.size, D_full.term_count(), k_full))
        print('  diagram of the SIMPLIFIED circuit: %d vertices, %d terms,'
              ' kappa %.3e' % (D_red.size, D_red.term_count(), k_red))
        print('  kappa change: %.2fx   (%.1f s)'
              % (k_full / k_red if k_red else float('inf'), time.time() - t1))
        print('  GATE 3 (kappa falls by >= 10x): %s'
              % ('PASS' if k_red and k_full / k_red >= 10 else 'FAIL'))
        print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
