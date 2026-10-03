"""SIZE one piece of a step, standalone and warm (2026-10-03; robust
timing).  The speed round that planned the Newton iterate in C sized its
pieces from a timer tree of nested inclusive wrappers, which RANKS pieces
but inflates each by ~1 us per nested timer: the evaluation it planned
around read 143 of 392 us there and was 66 of 355 measured standalone.  This
is the standalone measurement, through pyperf (calibrated loops, several
worker processes, warm-ups, the spread reported), on the harness's own
circuits (`step_machinery.BUILD`) a few steps into a transient.

    python benchmarks/micro.py --list
    python benchmarks/micro.py residual mos1        # one piece, one case
    python benchmarks/micro.py all mos1 -o out.json # every piece; pyperf options pass through

pyperf runs this script again in each worker; `--copy-env` keeps the thread
variables and `PYTHONPATH`, `--affinity` pins the workers (CPU 8 by default,
`_bench.CPU`).  Piece names map to callables at a state (`PIECES`).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _bench

_bench.pin_threads()
import warnings

warnings.simplefilter('ignore')
import numpy as np
import pyperf
import step_machinery as sm

from pycircuit.circuit import circuit
from pycircuit.circuit.transient import Transient


def _state(case):
    c = sm.BUILD[case]()
    tr = Transient(c, toolkit=circuit.numeric)
    res = tr.solve(tend=6 * 2e-8, timestep=2e-8, fixed_timestep=True)
    x = np.ascontiguousarray(np.asarray(res.x, float)[:, -1])
    tr._u_memo = {}
    return c, tr, x


def _pieces(case):
    c, tr, x = _state(case)
    t = 1.1e-7
    ep = tr.epar
    ir = tr.irefnode
    xr = np.ascontiguousarray(np.concatenate((x[:ir], x[ir + 1:])))
    f, J = tr._residual_and_jacobian(x, t)
    from pycircuit.circuit.analysis import _reduce_ndarray
    Jr, Fr = _reduce_ndarray(J, ir), _reduce_ndarray(f, ir)
    lim = tr._newton_limiter()
    return {
        'residual': lambda: tr._residual_and_jacobian(x, t),
        'pass_G': lambda: c.G(x, ep),
        'pass_C': lambda: c.C(x, ep),
        'pass_i': lambda: c.i(x, ep),
        'pass_q': lambda: c.q(x, ep),
        'sources': lambda: c.u(t, ep, 'tran'),
        'limit': (lambda: lim(xr.copy(), xr)) if lim is not None else None,
        'reduce': lambda: (_reduce_ndarray(J, ir), _reduce_ndarray(f, ir)),
        'solve': lambda: np.linalg.solve(Jr, -Fr),
        'predict': lambda: tr._predict_state(t),
    }


def main():
    argv = sys.argv[1:]
    if '--worker' not in argv:
        ## the master: the piece and the case travel to pyperf's workers in
        ## the environment (`--copy-env`), its own options in `sys.argv`
        if '--list' in argv:
            print('pieces:', ' '.join(_pieces('stage')))
            print('cases :', ' '.join(sm.BUILD))
            return
        if len(argv) < 2 or argv[0].startswith('-'):
            print(__doc__)
            return
        os.environ['PYCIRCUIT_MICRO'] = f'{argv[0]}:{argv[1]}'
        rest = argv[2:]
        if not any(a.startswith('--affinity') for a in rest):
            rest.append(f'--affinity={_bench.CPU}')
        if '--copy-env' not in rest:
            rest.append('--copy-env')
        sys.argv = [sys.argv[0]] + rest
    piece, case = os.environ['PYCIRCUIT_MICRO'].split(':')
    runner = pyperf.Runner()
    runner.metadata['case'] = case
    pieces = _pieces(case)
    names = list(pieces) if piece == 'all' else [piece]
    for nm in names:
        fn = pieces[nm]
        if fn is not None:
            runner.bench_func(f'{case}.{nm}', fn)


if __name__ == '__main__':
    main()
