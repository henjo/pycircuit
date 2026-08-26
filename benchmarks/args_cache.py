# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""What memoising `hdl._args_of` buys (roadmap sec. 23).

`_args_of` builds the trailing argument list every compiled function is
called with, so it runs once per ``i``/``G``/``q``/``C``/``u`` call --
once per element per Newton iteration per timestep.  Profiling a
49-element transient counted 125 029 calls, each rebuilding an identical
list.  Parameters do not change during a transient, so it is memoised on
the instance and dropped by the ``iparv`` observer.

This script measures the shipped implementation against an uncached twin
IN ONE PROCESS, interleaved and repeated, and asserts the answers stay
bit-identical.  Run it from the repository root::

    python benchmarks/args_cache.py

Recorded 2026-08-26 (4-core Linux, numpy backend), re-run after the
element-evaluation work of roadmap sec. 26::

      devs uncached_ms   cached_ms   speedup
        13       156.5       135.4     1.16x
        49       463.3       380.2     1.22x
       121      1001.5       809.4     1.24x

Quote **1.24x**, and re-run this rather than quoting it: what the cache
is worth **depends on what else the call does**.  It measured 1.16x when
it landed and 1.24x an hour later, unchanged, because sec. 26 removed a
425 ns `epar` lookup from the same method -- taking overhead out of the
denominator raises the share this cache accounts for.  A speedup figure
is a property of the whole call, not of the change alone.

Two earlier figures were higher and both were wrong to publish: a
ceiling probe taken before the implementation said 1.20x (it cached
without invalidation, which is cheaper than the real thing), and a first
interleaved run said 1.18x (taken before the cache key was made
type-exact).
"""

import time
import warnings

import numpy as np

warnings.simplefilter('ignore')

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd, elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import hdl as H

cm.default_toolkit = numeric

_real = H._args_of


def _uncached(self, epar):
    """`_args_of` as it was before the cache -- the A side."""
    n = len(self._hdl_paramnames)
    v = H._params_of(self)
    return v[:n] + [H._epar_T(epar)] + v[n:]


def build(n):
    """`n` cascaded diode-RC stages: 3n+1 elements, n+1 nodes."""
    c = SubCircuit()
    c['vs'] = eh.VSinHdl('n0', gnd, va=5.0, freq=1e3)
    for k in range(n):
        c['D%d' % k] = eh.DiodeHdl('n%d' % k, 'n%d' % (k + 1))
        c['R%d' % k] = eh.RHdl('n%d' % (k + 1), gnd, r=1e3)
        c['C%d' % k] = eh.CHdl('n%d' % (k + 1), gnd, c=1e-7)
    return c


def timed(n, fn):
    H._args_of = fn
    try:
        c = build(n)
        t = time.perf_counter()
        r = Transient(c).solve(refnode=gnd, tend=5e-4, timestep=1e-5)
        return time.perf_counter() - t, np.asarray(r.v('n%d' % n, gnd))
    finally:
        H._args_of = _real


def main(sizes=(4, 16, 40), repeats=3):
    print('%6s %11s %11s %9s'
          % ('devs', 'uncached_ms', 'cached_ms', 'speedup'))
    for n in sizes:
        us, cs = [], []
        for _ in range(repeats):
            u, vu = timed(n, _uncached)
            c, vc = timed(n, _real)
            ## The cache may not move a single bit of the answer.
            assert np.array_equal(vu, vc), 'ANSWERS DIFFER at n=%d' % n
            us.append(u)
            cs.append(c)
        ## Minimum, not mean: the fastest run is the one least polluted
        ## by whatever else the machine was doing.
        u, c = min(us), min(cs)
        print('%6d %11.1f %11.1f %8.2fx' % (3 * n + 1, u * 1e3, c * 1e3,
                                            u / c))
    print()
    print('answers bit-identical at every size and repeat.')


if __name__ == '__main__':
    main()
