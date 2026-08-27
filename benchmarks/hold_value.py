# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Which of a compact model's `var()` holds are LOAD-BEARING?

`psp_kernel.py` hand-tags **342** intermediates with `_v()`.  Every one
is assumed necessary and none is measured, which is the state roadmap
S4's argument is made from.  This answers it mechanically: replace `_v`
with the identity -- globally, or only for names matching a prefix --
rebuild the model, and report what changed.

Three properties are reported, because a hold can be load-bearing for
any of them independently:

* **compile time** -- without holds an expression that mentions a
  previous result twice doubles the tree, so a model nested `n` deep has
  ``2**n`` occurrences and every sympy traversal walks all of them;
* **emitted size** -- definitions inlined rather than named;
* **finiteness reach** -- the largest bias at which `i` and `G` are
  still finite, against the same ladder `test_psp_gap.py` uses.

Run from the repository root, with the compile cache OFF (the variants
share their `analog` source, and the compile is what is being timed)::

    PYCIRCUIT_HDL_CACHE=0 python benchmarks/hold_value.py           # survey
    PYCIRCUIT_HDL_CACHE=0 python benchmarks/hold_value.py sg_ ids_  # named

Recorded 2026-08-27, IHP sg13g2 nmos, w=10um l=1um:

===================  ==========  ===========  ==============
dropped              compile      G lines      finite out to
===================  ==========  ===========  ==============
nothing (baseline)      62 s         2675          1e36
``sg_*``                65 s         2575          1e36
everything            **>600 s**      --            --
===================  ==========  ===========  ==============

So `var()` is load-bearing for **compilability** and the case is
overwhelming.  The `sg_*` holds are **not** load-bearing for finiteness
today, though `_sigma_body`'s docstring records that they once were --
see roadmap sec. 36.
"""

import os
import sys
import time
import warnings

import numpy as np

warnings.simplefilter('ignore')

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric

cm.default_toolkit = numeric

#: The same ladder `test_psp_gap.py::test_the_jacobian_stays_finite_with_a_card`
#: walks, extended one decade past it so an IMPROVEMENT is visible too.
REACH = [1.0e2, 1.0e3, 1.0e4, 1.0e7, 1.0e12, 1.0e20, 1.0e26, 1.0e30,
         1.0e33, 1.0e36]

PDK = os.path.expanduser(
    '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models')


def _card():
    from pycircuit.circuit import psp_scaling
    from pycircuit.utilities import spicecard
    deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                          section='mos_tt')
    w, l = 10e-6, 1e-6
    return psp_scaling.to_long_channel(
        deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1, m=1,
                          pre_layout=1), w=w, l=l, T=273.15 + 27.0)


def measure(drop=None):
    """Build PSP with some holds dropped; return (seconds, lines, reach).

    `drop` is None or `()` (keep every hold), ``'*'`` (drop them all) or
    a tuple of name prefixes.  A survey runs each of these through
    `_one`, which bounds it in a subprocess.
    """
    from pycircuit.circuit import psp_kernel
    real = psp_kernel._v
    if drop == '*':
        psp_kernel._v = lambda expr, name=None: expr
    elif drop:
        def sel(expr, name=None, _p=tuple(drop)):
            return expr if (name and name.startswith(_p)) else real(expr, name)
        psp_kernel._v = sel
    try:
        t0 = time.perf_counter()
        from pycircuit.circuit.compact import PspMosLongChannel as P
        ## A fresh subclass so the metaclass recompiles under the patch.
        cls = type('PspHoldProbe', (P,),
                   dict(analog=staticmethod(P.analog),
                        instparams=list(P.instparams), __module__=P.__module__))
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                **_card())
        e.update_iparv()
        secs = time.perf_counter() - t0
    finally:
        psp_kernel._v = real

    src = e._hdl_info['funcs']['G'].__dict__.get('_src', '')
    reach = None
    for v in REACH:
        x = e.bias(v, 1.2, 0.0, 0.0)
        with np.errstate(all='ignore'):
            i = np.asarray(e.i(x), float)
            g = np.asarray(e.G(x), float)
        if not (np.all(np.isfinite(i)) and np.all(np.isfinite(g))):
            break
        reach = v
    return secs, len(src.splitlines()), reach


#: Per-build ceiling for a survey, in seconds.  **A build that exceeds it
#: is a RESULT, not a hang**: the point of dropping a hold is that the
#: tree may blow up.  Learned the hard way -- the first survey sat for
#: thirty minutes on one prefix with nothing printed, because `tail`
#: buffers and the blow-up is inside sympy where it cannot be
#: interrupted cooperatively.
BUILD_BUDGET = 300.0


def _one(drop, budget):
    """One measurement, in a bounded SUBPROCESS.

    A subprocess because the blow-up cannot be interrupted from inside,
    and because a half-built class would poison the next measurement in
    the same interpreter.
    """
    import subprocess
    arg = '*' if drop == '*' else ','.join(drop or ())
    try:
        out = subprocess.run(
            [sys.executable, os.path.abspath(__file__), '--one', arg],
            capture_output=True, text=True, timeout=budget,
            env=dict(os.environ, PYCIRCUIT_HDL_CACHE='0'))
    except subprocess.TimeoutExpired:
        return None
    for line in out.stdout.splitlines():
        if line.startswith('RESULT '):
            _, secs, lines, reach = line.split()
            return (float(secs), int(lines),
                    None if reach == 'none' else float(reach))
    return None


def main(prefixes, budget=BUILD_BUDGET):
    print('%-22s %10s %9s %14s' % ('dropped', 'compile', 'G lines',
                                   'finite to'))
    rows = [('nothing (baseline)', ())]
    rows += [('%s*' % p, (p,)) for p in prefixes]
    rows += [('EVERYTHING', '*')]
    for tag, drop in rows:
        r = _one(drop, budget)
        if r is None:
            print('%-22s %10s %9s %14s'
                  % (tag, '>%.0fs' % budget, '--', 'DID NOT BUILD'))
        else:
            secs, lines, reach = r
            print('%-22s %9.1fs %9d %14s'
                  % (tag, secs, lines, ('%.0e' % reach) if reach else 'FAILED'))
    print()
    print('">%.0fs" is a RESULT, not a failure of the run: those holds are'
          % budget)
    print('load-bearing for COMPILABILITY, which is the property that '
          'matters most.')


if __name__ == '__main__':
    if len(sys.argv) > 2 and sys.argv[1] == '--one':
        a = sys.argv[2]
        drop = '*' if a == '*' else (tuple(a.split(',')) if a else None)
        _s, _l, _r = measure(drop)
        print('RESULT %.3f %d %s'
              % (_s, _l, 'none' if _r is None else repr(_r)))
    else:
        main(sys.argv[1:] or ['sg_', 'ids_', 'avl_', 'ov_', 'gc_'])
