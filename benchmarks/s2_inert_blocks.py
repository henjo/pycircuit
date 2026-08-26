# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""What would S2 buy on PSP?  (roadmap sec. 29 -- answer: nothing.)

Run it from the repository root, with the compile cache OFF::

    PYCIRCUIT_HDL_CACHE=0 python benchmarks/s2_inert_blocks.py

Recorded 2026-08-26: eliding the `mob` switches takes the generated `G`
from 2675 to 2469 lines and measures **0.996x** on numpy, **1.000x** on
the C backend, bit-identical.  Eliding the whole gate-current block
(2275 lines) is **1.000x**.  The cost is the regularisers, not the
physics: `minc`/`maxc` are 4.5 ms of the 17.2 ms Jacobian.


S2 would compile, per parameter mask, a variant with every inert block
elided.  `psp_kernel` already contains the elision -- its build-time
guards fire when a parameter arrives as a NUMBER rather than a symbol,
and today it never is.  So handing the kernel literal 0.0 for the
switches a card sets to zero produces exactly what S2 would emit, and
the two can be compared directly.

The compile cache MUST be off: the variant's `analog` source and
instparams are identical to the base's, so it would hash to the same key.
"""
import os, time, warnings, numpy as np
assert os.environ.get('PYCIRCUIT_HDL_CACHE') == '0', 'run with PYCIRCUIT_HDL_CACHE=0'
warnings.simplefilter('ignore')
import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
cm.default_toolkit = numeric
import pycircuit.circuit.psp_kernel as pk
from pycircuit.circuit.compact import PspMosLongChannel
from pycircuit.circuit.circuit import defaultepar

OFF = ('alp', 'alp1', 'alp2', 'kp', 'rs', 'rsg', 'rsb')

t0 = time.perf_counter()
class PspBase(PspMosLongChannel):
    analog = staticmethod(PspMosLongChannel.analog)
    instparams = list(PspMosLongChannel.instparams)
t_base = time.perf_counter() - t0

_orig = pk.intrinsic
def _elided(*a, **kw):
    mob = kw.get('mob')
    if mob:
        for k in OFF:
            if k in mob:
                mob[k] = 0.0
    return _orig(*a, **kw)

pk.intrinsic = _elided
try:
    t0 = time.perf_counter()
    class PspElided(PspMosLongChannel):
        analog = staticmethod(PspMosLongChannel.analog)
        instparams = list(PspMosLongChannel.instparams)
    t_elid = time.perf_counter() - t0
finally:
    pk.intrinsic = _orig

def stats(cls):
    info = cls._hdl_info
    src = info['funcs']['G'].__dict__.get('_src', '')
    return len(src.splitlines()), src.count('\n')

for tag, cls, tc in (('base (all blocks)', PspBase, t_base),
                     ('elided (S2 variant)', PspElided, t_elid)):
    n, _ = stats(cls)
    print('%-22s G source %5d lines   compile %6.1f s' % (tag, n, tc))

def bench(cls, k=120):
    el = cls('d', 'g', 's', 'b')
    el.update_iparv()
    x = np.zeros(el.n); x[0] = 1.0; x[1] = 0.8
    for _ in range(5):
        el.G(x, defaultepar)
    t = time.perf_counter()
    for _ in range(k):
        el.G(x, defaultepar)
    g = np.asarray(el.G(x, defaultepar), dtype=float)
    i = np.asarray(el.i(x, defaultepar), dtype=float)
    return (time.perf_counter() - t) / k * 1e3, g, i

## interleaved and repeated -- a single run of this has misled me before
bs, es = [], []
for _ in range(4):
    tb, gb, ib = bench(PspBase)
    te, ge, ie = bench(PspElided)
    bs.append(tb); es.append(te)
tb, te = min(bs), min(es)
print()
print('G evaluation (min of 4, interleaved)')
print('   base   %7.3f ms   %s' % (tb, ' '.join('%.2f' % v for v in bs)))
print('   elided %7.3f ms   %s' % (te, ' '.join('%.2f' % v for v in es)))
print('   speedup %.3fx' % (tb / te))
print('i/G identical: %s / %s   max|dG| %.3e  max|di| %.3e'
      % (np.array_equal(ib, ie), np.array_equal(gb, ge),
         float(np.max(np.abs(gb - ge))), float(np.max(np.abs(ib - ie)))))
