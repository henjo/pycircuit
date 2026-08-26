# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`hdl.fold_card` on PSP with a real IHP sg13g2 card (roadmap sec. 30-31).

Run from the repository root with the compile cache OFF -- the folded
variant shares its `analog` source with the base, and the point of the
run is to time the compile too::

    PYCIRCUIT_HDL_CACHE=0 python benchmarks/card_constant_folding.py

Recorded 2026-08-26, w=10um l=1um, 66 card parameters folded, w/l left
per-instance::

    symbolic  G  2675 lines   primitives  11565
    folded    G  1971 lines   primitives   7251
    fold_card compile 17.2 s   (base builds in ~28 s)
    numpy   symbolic  17.73 ms   folded  11.55 ms   1.535x
    C       symbolic   83.5 us   folded   56.5 us   1.478x
            max relative dG 1.44e-13

⚠ Fold EVERYTHING outside `instance`.  Folding only the parameters the
card happens to MENTION leaves the rest symbolic -- on PSP that was 87
of 153 -- and measures **0.96x**, a slowdown.
"""
import os, time, warnings, numpy as np
warnings.simplefilter('ignore')
import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
cm.default_toolkit = numeric
from pycircuit.circuit import hdl, psp_scaling
from pycircuit.utilities import spicecard
from pycircuit.circuit.compact import PspMosLongChannel
from pycircuit.circuit.circuit import defaultepar

PDK = os.path.expanduser('~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'), section='mos_tt')
w, l = 10e-6, 1e-6
CARD = psp_scaling.to_long_channel(
    deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1, pre_layout=1),
    w=w, l=l, T=273.15 + 27.0)
names = {p.name for p in PspMosLongChannel.instparams}
KW = {k: v for k, v in CARD.items() if k in names and isinstance(v, (int, float))}
## geometry stays an INSTANCE parameter -- that is the whole point of the split
INST = {k: KW.pop(k) for k in ('w', 'l') if k in KW}
print('card folded: %d parameters   instance: %s' % (len(KW), sorted(INST)))

class PspSym(PspMosLongChannel):
    analog = staticmethod(PspMosLongChannel.analog)
    instparams = list(PspMosLongChannel.instparams)

t = time.perf_counter()
PspFold = hdl.fold_card(PspSym, instance=('w', 'l'), **KW)
PspFold._hdl_info  # force compile
print('fold_card compile: %.1f s' % (time.perf_counter() - t))

for tag, cls in (('symbolic', PspSym), ('folded', PspFold)):
    src = cls._hdl_info['funcs']['G'].__dict__.get('_src', '')
    prims = sum(src.count(k+'(') for k in ('minc','maxc','_step','_rdiv','_recip2'))
    print('%-9s G %5d lines   primitives %6d' % (tag, len(src.splitlines()), prims))

def mk(cls):
    kw = dict(INST) if cls is PspFold else dict(KW, **INST)
    el = cls('d','g','s','b', **kw); el.update_iparv(); return el

def bench(cls, k, us):
    el = mk(cls)
    x = np.zeros(el.n); x[0]=1.0; x[1]=0.8
    for _ in range(20 if us else 5): el.G(x, defaultepar)
    t=time.perf_counter()
    for _ in range(k): el.G(x, defaultepar)
    return ((time.perf_counter()-t)/k*(1e6 if us else 1e3),
            np.asarray(el.G(x, defaultepar), float))

a,b = [],[]
for _ in range(4):
    ta,ga = bench(PspSym,120,False); tb,gb = bench(PspFold,120,False)
    a.append(ta); b.append(tb)
print()
print('numpy   symbolic %7.2f ms   folded %7.2f ms   %.3fx' % (min(a), min(b), min(a)/min(b)))

for cls in (PspSym, PspFold):
    hdl.set_backend('c', cls)
a,b = [],[]
for _ in range(4):
    ta,ga = bench(PspSym,3000,True); tb,gb = bench(PspFold,3000,True)
    a.append(ta); b.append(tb)
print('C       symbolic %7.1f us   folded %7.1f us   %.3fx' % (min(a), min(b), min(a)/min(b)))
nz = np.abs(ga) > 0
print('        max relative dG %.2e' % float(np.max(np.abs(ga[nz]-gb[nz])/np.abs(ga[nz]))))
