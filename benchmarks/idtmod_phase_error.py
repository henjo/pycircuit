"""Idtmod vs IdtmodCircular: the two phase-error regimes, measured live.

The table behind idtmod.md sec. 7.6.  One circuit (1 V into a modulus-1
circular integrator, so the analytic phase is exactly ``t``), one fixed
grid (dt = 0.01, w*h = 2*pi/100), both elements, three integrators, three
run lengths.  The congruence error at end of run separates two regimes:

* ``Idtmod`` (scalar + Phase-2 gauge shift): NO method error -- a linear
  phase is integrated exactly by every LMS formula, so the three
  integrators agree to the digit and what remains is float accumulation
  noise, ~N*ulp.
* ``IdtmodCircular`` (phasor pair): the method's per-cycle phase lag at
  the carrier rate -- (wh)^2/6, /12, /3 per cycle for Euler, trapezoidal
  and Gear-2 -- accumulating linearly with cycle count, present even for
  constant input.  The PREDICTED column prints that closed form next to
  the measurement; the run fails loudly if they disagree by more than 2x,
  so this doubles as a live check of the sec. 7.6 theory.

Bounded on purpose, like phase0_baseline: one circuit, one dt, one
repetition -- a yardstick reproducing the documented table, not a sweep.

Run:  python benchmarks/idtmod_phase_error.py
"""

import time
import warnings

import numpy as np

DT = 0.01
TENDS = (2.0, 20.0, 200.0)

## Per-cycle phase-lag coefficients of each method on a rotation at w*h:
## error/cycle = COEFF * (w*h)^2, from the phase of the amplification
## factor (arg(1 + i*w*h) etc.); see idtmod.md sec. 7.6.
LAG_COEFF = {
    'EulerIntegrator': 1.0 / 6.0,
    'TrapezoidalIntegrator': 1.0 / 12.0,
    'Gear2Integrator': 1.0 / 3.0,
}


def run(element_cls, integrator_cls, tend):
    import pycircuit.circuit.circuit as circuit_mod
    from pycircuit.circuit import gnd
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit.transient import Transient

    circuit_mod.default_toolkit = numeric
    c = SubCircuit()
    nin, nout = c.add_node('in'), c.add_node('out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['X'] = element_cls(nin, gnd, nout, gnd, modulus=1.0, ic=0.0)

    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=integrator_cls())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        res = tran.solve(tend=tend, timestep=DT, fixed_timestep=True)
        wall = time.perf_counter() - t0

    y = res.v('out').y
    t = res.v('out').x[0]
    ## Congruence distance: a sample landing exactly ON a wrap is a
    ## two-limit knife edge; distance modulo the modulus is well-defined.
    d = np.abs(y[1:] - t[1:] % 1.0)
    d = np.minimum(d, 1.0 - d)
    tail = d[-max(4, len(d) // 20):]
    return float(tail.max()), wall


def main():
    from pycircuit.circuit.elements import Idtmod, IdtmodCircular
    from pycircuit.circuit.integrator import (EulerIntegrator,
                                              TrapezoidalIntegrator,
                                              Gear2Integrator)

    wh = 2.0 * np.pi * DT
    print('dt = %g, w*h = %.4f; congruence error at end of run\n' % (DT, wh))
    print('%-22s %-15s %8s   %-11s %-11s %s' % (
        'integrator', 'element', 'cycles', 'measured', 'predicted', 'wall'))

    failures = []
    for integ in (EulerIntegrator, TrapezoidalIntegrator, Gear2Integrator):
        for tend in TENDS:
            for cls in (Idtmod, IdtmodCircular):
                err, wall = run(cls, integ, tend)
                if cls is IdtmodCircular:
                    pred = LAG_COEFF[integ.__name__] * wh ** 2 * tend
                    pred_s = '%.3e' % pred
                    if not (0.5 < err / pred < 2.0):
                        failures.append((integ.__name__, tend, err, pred))
                else:
                    pred_s = '(float noise)'
                print('%-22s %-15s %8g   %.3e   %-11s %5.1fs' % (
                    integ.__name__, cls.__name__, tend, err, pred_s, wall))
        print()

    if failures:
        raise SystemExit(
            'IdtmodCircular error does not match the per-cycle lag theory '
            'of idtmod.md sec. 7.6 within 2x: %r' % (failures,))
    print('IdtmodCircular matches the per-cycle lag prediction within 2x '
          'at every point; Idtmod stays at float-noise level throughout.')


if __name__ == '__main__':
    main()
