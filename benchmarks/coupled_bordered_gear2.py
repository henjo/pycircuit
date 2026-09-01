"""B3a: is the coupled 'bordered' branch sound under Gear-2?

`doc/backend_parity_260821.md` records that landing Gear-2 on the coupled path
"mistunes `bordered`'s grow-back (1181 points where 'approx' took 350 on the
pulsed RC)" and lists "re-derive `bordered`'s denominator pairing" as scope
that was never executed.  The JAX coupled path defaults to Gear-2, so porting
`bordered` there means porting it into the one configuration it is recorded as
misbehaving in.  This measures that before any JAX code is written.

Gate (set before running, from the plan): if bordered under Gear-2 still
spends ~3.4x the time points of 'approx', it is not sound and the port does
not happen.
"""
import argparse, time, warnings

import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VPulse, VSin
from pycircuit.circuit.integrator import EulerIntegrator, Gear2Integrator
from pycircuit.circuit.transient import Transient

TD, TR, PW, TF, PER = 1e-5, 1e-7, 2e-5, 1e-7, 5e-5


def pulsed_rc():
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=TD, tr=TR, tf=TF, pw=PW,
                     per=PER)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def smooth_rc():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


CIRCUITS = {'pulsed': (pulsed_rc, 1.2e-4, 1e-6),
            'smooth': (smooth_rc, 5e-3, 1e-4)}


def run(build, tend, step, method, integrator, reltol):
    tran = Transient(build(), toolkit=numeric, reltol=reltol,
                     integrator=integrator, coupled_method=method)
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=step, coupled_lte=True)
    wall = time.perf_counter() - t0
    t = np.asarray(res.v('b').x, float).ravel()
    v = np.asarray(res.v('b').y, float).ravel()
    st = tran.statistics
    return dict(points=len(t), newton=st.newton_iterations,
                rejected=st.rejected_steps, wall=wall, t=t, v=v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reltol', type=float, default=1e-5)
    ap.add_argument('--circuit', default='both',
                    choices=['pulsed', 'smooth', 'both'])
    args = ap.parse_args()

    names = ['pulsed', 'smooth'] if args.circuit == 'both' else [args.circuit]
    print('reltol %g\n' % args.reltol)
    print('%-8s %-8s %-9s %8s %9s %9s %8s' %
          ('circuit', 'method', 'integ', 'points', 'newton', 'rejected', 'wall s'))
    for name in names:
        build, tend, step = CIRCUITS[name]
        ref = {}
        for integ_name, integ in (('euler', EulerIntegrator),
                                  ('gear2', Gear2Integrator)):
            for method in ('approx', 'bordered'):
                r = run(build, tend, step, method, integ(), args.reltol)
                ref[(integ_name, method)] = r
                print('%-8s %-8s %-9s %8d %9d %9d %8.2f' %
                      (name, method, integ_name, r['points'], r['newton'],
                       r['rejected'], r['wall']))
        for integ_name in ('euler', 'gear2'):
            a = ref[(integ_name, 'approx')]
            b = ref[(integ_name, 'bordered')]
            n = min(len(a['t']), len(b['t']))
            dev = float(np.max(np.abs(
                b['v'] - np.interp(b['t'], a['t'], a['v']))))
            print('   %-6s %-6s bordered/approx: points %.2fx  newton %.2fx'
                  '   max|b-a| %.3e' %
                  (name, integ_name, b['points'] / a['points'],
                   b['newton'] / max(a['newton'], 1), dev))
        print()


if __name__ == '__main__':
    main()
