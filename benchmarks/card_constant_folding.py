# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Constant-fold every CARD parameter: does evaluation get cheaper?

Roadmap sec. 30.  Answer: **yes, 1.54x on PSP**, and to machine
precision in the regime the model is validated in.

Run from the repository root with the compile cache OFF (the folded
variant has the same `analog` source and instparams as the base, so it
would hash to the same key)::

    PYCIRCUIT_HDL_CACHE=0 python benchmarks/card_constant_folding.py

Recorded 2026-08-26, real IHP sg13g2 card, w=10um l=1um:
  G source        2675 -> 1971 lines   (-26%)
  primitive calls 11565 -> 7251        (-37%)
  G evaluation    17.90 -> 11.62 ms    (1.54x)
  drain current   8.5e-16 relative in strong inversion,
                  2.7e-15 in the validated 1e-9..1e-6 A window


`analog()` receives the parameters as sympy Symbols.  Handing it the
card's NUMBERS instead makes sympy fold every parameter-only
subexpression at build time -- `a*b*x` becomes `c*x`.  The compiled
function still takes the (now unused) parameter arguments, so nothing
downstream changes.
"""
import os, time, warnings, numpy as np, sympy
assert os.environ.get('PYCIRCUIT_HDL_CACHE') == '0'
warnings.simplefilter('ignore')
import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
cm.default_toolkit = numeric
from pycircuit.circuit import hdl
from pycircuit.circuit.compact import PspMosLongChannel
from pycircuit.circuit.circuit import defaultepar

import os as _os
from pycircuit.circuit import psp_scaling
from pycircuit.utilities import spicecard
PDK = _os.path.expanduser('~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')
_deck = spicecard.read(_os.path.join(PDK, 'cornerMOSlv.lib'), section='mos_tt')
_w, _l = 10e-6, 1e-6
CARD = psp_scaling.to_long_channel(
    _deck.model_params('sg13g2_lv_nmos_psp', w=_w, l=_l, ng=1, m=1,
                       pre_layout=1), w=_w, l=_l, T=273.15 + 27.0)
DEF = {p.name: p.default for p in PspMosLongChannel.instparams}
DEF.update({k: v for k, v in CARD.items() if isinstance(v, (int, float))})
print('real card: %d of %d parameters supplied' % (len(CARD), len(DEF)))

class PspSym(PspMosLongChannel):
    analog = staticmethod(PspMosLongChannel.analog)
    instparams = list(PspMosLongChannel.instparams)

_orig = hdl._analog_function
def _folded(func, paramnames, paramsyms):
    vals = [sympy.Float(DEF[n]) if isinstance(DEF.get(n), (int, float))
            else paramsyms[i] for i, n in enumerate(paramnames)]
    return _orig(func, paramnames, vals)

hdl._analog_function = _folded
try:
    class PspFold(PspMosLongChannel):
        analog = staticmethod(PspMosLongChannel.analog)
        instparams = list(PspMosLongChannel.instparams)
    ok = True
except Exception as e:
    ok = False
    print('BUILD FAILED with folded parameters: %s: %s' % (type(e).__name__, str(e)[:160]))
finally:
    hdl._analog_function = _orig

for tag, cls in [('symbolic params', PspSym)] + ([('folded params', PspFold)] if ok else []):
    src = cls._hdl_info['funcs']['G'].__dict__.get('_src', '')
    import re
    prims = sum(src.count(k+'(') for k in ('minc','maxc','_step','_rdiv','_recip2'))
    print('%-18s G lines %5d   primitive calls in src %6d'
          % (tag, len(src.splitlines()), prims))

if ok:
    names = {p.name for p in PspSym.instparams}
    KW = {k: v for k, v in CARD.items() if k in names}

    def bench(cls, k=120):
        el = cls('d', 'g', 's', 'b', **KW)
        el.update_iparv()
        x = np.zeros(el.n); x[0] = 1.0; x[1] = 0.8
        for _ in range(5):
            el.G(x, defaultepar)
        t = time.perf_counter()
        for _ in range(k):
            el.G(x, defaultepar)
        return ((time.perf_counter() - t) / k * 1e3,
                np.asarray(el.G(x, defaultepar), float))
    a, b = [], []
    for _ in range(4):
        ta, ga = bench(PspSym); tb, gb = bench(PspFold)
        a.append(ta); b.append(tb)
    print()
    print('G  symbolic %7.3f ms  %s' % (min(a), ' '.join('%.2f' % v for v in a)))
    print('   folded   %7.3f ms  %s' % (min(b), ' '.join('%.2f' % v for v in b)))
    print('   speedup %.3fx' % (min(a)/min(b)))
    nz = np.abs(ga) > 0
    rel = np.max(np.abs(ga[nz]-gb[nz])/np.abs(ga[nz])) if nz.any() else 0.0
    print('   G: max|abs| %.3e   max|rel| %.3e   (|G| up to %.3e)'
          % (float(np.max(np.abs(ga-gb))), float(rel), float(np.max(np.abs(ga)))))
    ## The drain current on a proper gate sweep, banded by magnitude.
    ## PSP is validated against the vendor to 1.3e-6 in the window
    ## 1e-9..1e-6 A -- so that band is the one that decides this.
    e1 = PspSym('d','g','s','b', **KW); e1.update_iparv()
    e2 = PspFold('d','g','s','b', **KW); e2.update_iparv()
    bands = {'>1e-6 A (strong)': [], '1e-9..1e-6 A (validated)': [],
             '<1e-9 A (below the reference floor)': []}
    for vg in np.linspace(0.0, 1.2, 121):
        xx = np.zeros(e1.n); xx[0] = 1.0; xx[1] = vg
        i1 = np.asarray(e1.i(xx, defaultepar), float)
        i2 = np.asarray(e2.i(xx, defaultepar), float)
        d = float(np.max(np.abs(i1)))
        if d == 0.0:
            continue
        r = float(np.max(np.abs(i1 - i2)) / d)
        k = ('>1e-6 A (strong)' if d > 1e-6 else
             '1e-9..1e-6 A (validated)' if d > 1e-9 else
             '<1e-9 A (below the reference floor)')
        bands[k].append(r)
    print()
    print('   drain current, gate sweep 0..1.2 V, by magnitude band:')
    for k, v in bands.items():
        if v:
            print('     %-38s n=%3d  max rel %.3e' % (k, len(v), max(v)))
