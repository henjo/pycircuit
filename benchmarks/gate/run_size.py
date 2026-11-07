"""One circuit at size: `run_size.py kind size method tend` -> one JSON line: the digest, the statistics, the KLU's counts, the swing and
the phases a timing pair reads (the operating point, the C loops' calls and glue, the Python's rest), the peak.  Run from the tree under
test (cwd and PYTHONPATH): `file` says which package it imported."""
import functools
import hashlib
import json
import os
import resource
import sys
import time
import warnings

warnings.simplefilter('ignore')
import circuits
import numpy as np

import pycircuit
from pycircuit.circuit import (
    _paths,
    _tran_initial,
    _tran_loop_c,
    _tran_tcloop_c,
    circuit,
)
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.transient import Transient

KLU = ('factors', 'refactors', 'pivot_fallbacks', 'residual_fallbacks', 'chord_factors', 'chord_refactors')
T = {'op': 0.0, 'loop_call': 0.0, 'tc_call': 0.0, 'run': 0.0}


def _timed(key, f):
    @functools.wraps(f)
    def w(*a, **k):
        t = time.perf_counter()
        try:
            return f(*a, **k)
        finally:
            T[key] += time.perf_counter() - t
    return w


def main():
    kind, size, method, tend = sys.argv[1], sys.argv[2], sys.argv[3], float(sys.argv[4])
    _tran_initial._InitialState._solve_operating_point = _timed('op', _tran_initial._InitialState._solve_operating_point)
    _tran_loop_c._call = _timed('loop_call', _tran_loop_c._call)
    _tran_tcloop_c._call = _timed('tc_call', _tran_tcloop_c._call)
    _tran_loop_c.run = _timed('run', _tran_loop_c.run)
    c = circuits.make(kind, size)
    tr = Transient(c, toolkit=circuit.numeric, **({'integrator': RadauIIA3Integrator()} if method == 'radau' else {}))
    b = _paths.snapshot()
    t0 = time.perf_counter()
    res = tr.solve(tend=tend, timestep=1e-9)
    tot = time.perf_counter() - t0
    x = np.asarray(res.x, float)
    t = np.asarray(res.sweep_values, float)
    h = hashlib.sha1()
    for row in x:
        h.update(row.tobytes())
    h.update(t.tobytes())
    ax = 0 if x.shape[0] == t.size else -1
    st = tr.statistics
    ch = tr._get_linearsolver()._choice
    print(json.dumps({'kind': kind, 'size': size, 'method': method, 'tend': tend, 'n': tr.cir.n, 'x': h.hexdigest()[:12], 'file': pycircuit.__file__,
                          'solve': round(tot, 3), 'op': round(T['op'], 3), 'c_call': round(T['loop_call'] + T['tc_call'], 3),
                          'python': round(tot - T['op'] - T['run'], 3), 'steps': st.accepted_steps, 'rej': st.rejected_steps,
                          'newton': st.newton_iterations, 'fb': getattr(st, 'chord_fallbacks', 0),
                          'klu': {k: getattr(ch, k, None) for k in KLU}, 'swing': int(((x.max(axis=ax) - x.min(axis=ax)) > 0.1).sum()),
                          'peak': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, 'load': os.getloadavg()[0],
                          'counts': _paths.since(b)}))


if __name__ == '__main__':
    main()
