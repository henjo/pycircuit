"""GATE 4-D, first half: does the failure mode the guard exists to prevent occur?

Decision 0.3d kept a charge-domain check as a "step-collapse guard for the
near-singular case": if `J` is near-singular, `J^-1` amplifies the charge
residual unboundedly and the solution-domain LTE becomes noise.  A floating or
weakly-grounded node is the case to try.

Measured here: how often the `linearsolver` in `IntegralController` actually
fails (the silent fallback fires), and how ill-conditioned J_reduced gets.
"""
import warnings
import numpy as np
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import stepcontroller as sc
from pycircuit.circuit.analysis import remove_row_col


def healthy():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3); c['C'] = C('b', gnd, c=1e-7)
    return c, dict(tend=2e-4, timestep=1e-5)

def series_caps():
    """Midpoint 'm' has no DC path to ground -- SPICE's classic case."""
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C1'] = C('b', 'm', c=1e-7); c['C2'] = C('m', gnd, c=1e-7)
    return c, dict(tend=2e-4, timestep=1e-5)

def dangling_cap():
    """'f' hangs off 'b' through a capacitor and touches nothing else."""
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3); c['C'] = C('b', gnd, c=1e-7)
    c['Cf'] = C('b', 'f', c=1e-12)
    return c, dict(tend=2e-4, timestep=1e-5)

def weak_ground():
    """'m' reaches ground only through 1 Tohm."""
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C1'] = C('b', 'm', c=1e-7); c['Rleak'] = R('m', gnd, r=1e12)
    return c, dict(tend=2e-4, timestep=1e-5)

CASES = [('healthy', healthy), ('series-caps', series_caps),
         ('dangling-cap', dangling_cap), ('weak-ground', weak_ground)]

print('%-14s %8s %9s %12s %12s  %s'
      % ('circuit', 'steps', 'fallback', 'max cond(J)', 'max |lte|', 'outcome'))
for name, build in CASES:
    stats = {'fallback': 0, 'cond': 0.0, 'lte': 0.0, 'calls': 0}
    real = sc.IntegralController.evaluate_step
    orig_solve = numeric.linearsolver
    def spy_solve(A, b, _o=orig_solve, _s=stats):
        _s['calls'] += 1
        try:
            c = np.linalg.cond(np.asarray(A, dtype=float))
            if np.isfinite(c):
                _s['cond'] = max(_s['cond'], c)
        except Exception:
            pass
        try:
            out = _o(A, b)
        except Exception:
            _s['fallback'] += 1
            raise
        _s['lte'] = max(_s['lte'], float(np.max(np.abs(out))) if np.size(out) else 0.0)
        return out
    numeric.linearsolver = spy_solve
    try:
        cir, kw = build()
        tran = Transient(cir, toolkit=numeric, reltol=1e-5)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tran.solve(coupled_lte=False, **kw)
            out = 'completed'
            steps = tran.statistics.accepted_steps
        except Exception as e:
            out = '%s: %s' % (type(e).__name__, str(e)[:28])
            steps = getattr(getattr(tran, 'statistics', None), 'accepted_steps', -1)
    finally:
        numeric.linearsolver = orig_solve
    print('%-14s %8d %9d %12.3e %12.3e  %s'
          % (name, steps, stats['fallback'], stats['cond'], stats['lte'], out))
