"""The LTE-chosen grid's prize, measured -- and now reachable through PSS.

RECORDED SCOPE ITEM 5: a transient adapts because it cannot see the future;
PSS re-solves the SAME interval repeatedly, so it can be handed a grid that
was chosen well once and then frozen.  The grid still never moves inside a
solve, so the shooting Newton stays exact.

This driver keeps its OWN `(x_0, T)` unknowns and a FINITE-DIFFERENCE
Jacobian, because that is what made it a control:

    method  grid              steps  outcome        period       err ppm
    trap    LTE-chosen         1105  converged     162.834707      -47.3
    trap    uniform, same N    1105  NoConvergence      -              -
    trap    uniform (20000)   20000  converged     162.832543      -60.6
    gear    LTE-chosen         1105  NoConvergence      -              -

18x fewer points AND better accuracy than the uniform grid that does
converge -- and at 1105 points a uniform grid does not converge at all.

⚠ THE REAL ANALYSIS NOW REPRODUCES IT (2026-09-02): `PSS.solve(grid=...)`
converges on the same grid at -73.8 ppm (trap) and -100.6 ppm (gear), the
latter being a case NO driver here ever solved.

⚠ AND THE REASON RECORDED HERE FOR WHY IT DID NOT WAS WRONG.  This header
read: "the finite-difference Jacobian here is effectively exact; the
shipped plain path ... measures ~30% off ... a stiff relaxation oscillator
is simply the first circuit that cannot tolerate that 30%."  Measured
since: giving the real analysis an exact finite-difference Jacobian does
NOT fix it, and after the real fix both Jacobians agree to six digits.
The blocker was the MANUFACTURED OPENING STEP -- `h[0] = 1.4845` against a
median of 4.62e-04, an order-dropped Euler step the inner Newton could not
solve.  This driver never met it because its unknown is `x_0` itself.
`_period_grid` now opens a coarse-opening grid on its own finest step.

Kept as a control: it is the independent implementation the analysis is
checked against, and its -47.3 ppm is still the better answer.
"""
import warnings, sys
import numpy as np
from copy import copy
sys.path.insert(0, 'benchmarks')
from pss_stiff_autonomous import van_der_pol, VDP_PERIOD
from pycircuit import circuit
from pycircuit.circuit import gnd
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit import analysis as an
from pycircuit.circuit.analysis import remove_row_col

circuit.default_toolkit = circuit.numeric
T = VDP_PERIOD


def lte_grid(reltol=1e-7):
    """One period of ACCEPTED steps from an adaptive transient, as fractions.

    This is the whole idea of item 5: a transient adapts because it cannot
    see the future; PSS re-solves the SAME interval over and over, so it
    can be handed a grid that was chosen well once.
    """
    cir = van_der_pol()
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=reltol).solve(refnode=gnd, tend=1200.0 + T,
                                                  timestep=0.05, x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    xs = np.asarray(res.x, dtype=float)
    ## a settled window of exactly one period
    j0 = int(np.searchsorted(t, t[-1] - T))
    win_t, win_x = t[j0:], xs[:, j0:]
    hs = np.diff(win_t)
    hs = hs / hs.sum()                       # fractions of the period
    iref = cir.get_node_index(gnd)
    seed = np.concatenate((win_x[:iref, 0], win_x[iref + 1:, 0]))
    return hs, seed


def shoot(method, fracs, seed, reltol=1e-7, maxiter=25):
    """Shooting on a FROZEN, possibly non-uniform grid; FD Jacobian (m=2)."""
    ## No setup solve: `irefnode` is set in __init__ and `_transient()` is
    ## lazy, so the only thing to supply is the autonomy verdict that
    ## `_solves_history` consults.  (A uniform-grid setup solve would fail
    ## on this circuit for exactly the reason under test.)
    pss = PSS(van_der_pol(), method=method, reltol=reltol)
    pss.autonomous = True
    m = pss.cir.n - 1
    calls = [0]

    def phi(z):
        x0, per = np.asarray(z[:m], float), float(z[-1])
        hs = fracs * per
        ## trap has `b != 0` so `_install_history` refuses it -- it belongs
        ## on the plain path, opened with a flat seed and an Euler step.
        if pss._companion_reach() >= 2:
            pss._install_history(x0, x0, hs[0])
        else:
            pss._begin_period(x0)
        x, t = copy(x0), 0.0
        for h in hs:
            t += h
            x = copy(pss.solve_timestep(x, t, h))
        return np.asarray(x, float)

    k = int(np.argmax(np.abs(seed))) if np.max(np.abs(seed)) else 0
    pin = float(seed[k])

    def F(z):
        calls[0] += 1
        out = np.zeros(m + 1)
        out[:m] = np.asarray(z[:m], float) - phi(z)
        out[m] = float(z[k]) - pin
        return out

    def func(z):
        z = np.asarray(z, float); f0 = F(z)
        J = np.zeros((m + 1, m + 1))
        for j in range(m + 1):
            s = 1e-7 * max(abs(z[j]), 1.0 if j < m else T)
            zp = z.copy(); zp[j] += s
            J[:, j] = (F(zp) - f0) / s
        return f0, J

    tol = an.newton_tolerance_vectors(len(pss.cir.nodes),
                                      len(pss.cir.branches), pss.par.iabstol,
                                      pss.par.vabstol, pss.toolkit)[1]
    (tol,) = remove_row_col((tol,), pss.irefnode, pss.toolkit)
    try:
        z, _i, ier, _mg = an.fsolve(
            func, np.concatenate((seed, [T])), maxiter=maxiter, reltol=reltol,
            abstol=np.concatenate((tol, [tol[k]])),
            xtol=np.concatenate((tol, [1e-15 * T])),
            toolkit=pss.toolkit, full_output=True)
        return ('converged' if ier == 1 else 'NOT converged'), float(z[-1]), calls[0]
    except Exception as exc:
        return type(exc).__name__, None, calls[0]


fr, seed = lte_grid()
print('LTE-chosen grid: %d steps over one period' % len(fr))
print('  step spread: min %.3e  max %.3e  ratio %.0f'
      % (fr.min() * T, fr.max() * T, fr.max() / fr.min()))
print('  uniform equivalent would need h = %.4f to match the finest'
      % (fr.min() * T))
print('  -> a uniform grid at that step is %d points\n'
      % int(round(1.0 / fr.min())))

print('%-6s %-18s %7s %-16s %12s %10s'
      % ('method', 'grid', 'steps', 'outcome', 'period', 'err ppm'))
for m_ in ('trap', 'gear'):
    for label, g in (('LTE-chosen', fr),
                     ('uniform, same N', np.ones(len(fr)) / len(fr))):
        out, per, nc = shoot(m_, g, seed)
        err = ('%10.1f' % (1e6 * (per - T) / T)) if per else '%10s' % '-'
        print('%-6s %-18s %7d %-16s %12s %s'
              % (m_, label, len(g), out,
                 ('%12.6f' % per) if per else '%12s' % '-', err))
