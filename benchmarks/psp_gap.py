"""How far is our surface-potential core from the real PSP103?

Everything in `psp_kernel` and `compact.PspMosLongChannel` is validated
by construction properties -- exact source/drain symmetry, structural
charge conservation, the surface-potential equation satisfied to 1e-9.
Those are strong, and they say nothing about how close the device is to
the transistor a foundry ships.

This measures that, end to end and with nothing hand-tuned:

    IHP model card  ->  spicecard  ->  psp_scaling  ->  the element
                                                          |
                    committed PSP103 reference  <---- compared

Run:  python benchmarks/psp_gap.py [--pdk PATH]

The point is not the headline number but the SHAPE of the discrepancy,
which says what to build next.  Reading it:

  * READ THE SPREAD BEFORE THE MEDIAN.  The median says whether the
    GAIN is right; the spread across a sweep says whether the bias
    dependence is.  A flat ratio at 1.09 is one missing gain term -- a
    well-posed question with one answer.  A ratio sweeping from 1.03 to
    0.44 is a broken drive, and no median makes that the better model.
    Ranking terms by their effect on the median once kept a correct
    block switched off for a commit;
  * a ratio that is flat in Vds and near 1 on the long device means the
    threshold and gain factor are right;
  * a ratio that GROWS with current is series resistance, which PSP
    folds into the mobility as an extra `Gmob` term;
  * a ratio that differs between the two GEOMETRIES is a scaling term,
    and the fastest way to find which one is not to reason about the
    shape at all -- it is to compare every scaled parameter against
    PSP's own `lp_*` operating-point outputs, which is what
    `test_psp_gap.py` does at four geometries on both channel types.
    Three terms have been found that way, each of them invisible in the
    card because its coefficient was zero on the device being measured;
  * a ratio that differs between the two CHANNEL TYPES is either a term
    with a genuinely different form for holes -- there are exactly four
    of those, all in `psp_kernel` -- or, more likely on the evidence so
    far, a scaling that one of the two cards switches off.

What the core deliberately lacks, and therefore what this measures the
cost of: `PSCE` (all-zero on this card anyway), gate and junction
leakage, impact ionisation, overlap and fringe capacitances, the
non-quasi-static block, self-heating, and every temperature coefficient.

Three of those cannot appear here at ALL, which is worth stating rather
than leaving for someone to rediscover: gate and junction leakage are
four to six orders of magnitude below the 1e-6 A floor below, and
overlap and fringe capacitances contribute identically zero DC current.
They matter to a CV comparison, which this is not.
"""

import argparse
import json
import os
import sys
import warnings

import numpy as np


REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..',
                   'pycircuit', 'circuit', 'tests', 'data',
                   'psp103_ihp_iv.json')

#: n-channel sweeps from the committed reference, long device first.
SWEEPS = ['nmos_long_idvd', 'nmos_idvd_vg1p2', 'nmos_idvd_vg0p6',
          'nmos_idvg_vd0p05', 'nmos_idvg_vd1p2', 'nmos_idvg_vb_m1',
          'nmos_long_idvg_vb_m1',
          'pmos_long_idvd', 'pmos_long_idvg', 'pmos_idvd_vg1p2',
          'pmos_idvg_vd1p2', 'pmos_idvg_vd0p05']

#: Currents below this are not compared: the reference carries leakage
#: and junction terms the core does not model at all, and they dominate
#: there.
FLOOR = 1e-6


def compare(deck, sweep):
    """Ratio of our drain current to PSP103's, over one sweep."""
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.compact import (PspMosLongChannel,
                                           PspPmosLongChannel)
    from pycircuit.circuit import psp_scaling

    w, l = sweep['w'], sweep['l']
    ## The sweep names its own device, so the p-channel curves pick the
    ## p-channel card and the p-channel class without a second code path.
    pmos = 'pmos' in sweep['device']
    model = 'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp'
    cls = PspPmosLongChannel if pmos else PspMosLongChannel
    card = deck.model_params(model, w=w, l=l, ng=1, m=1, pre_layout=1)
    kw = psp_scaling.to_long_channel(card, w=w, l=l)
    e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'), **kw)
    e.update_iparv()

    v = np.asarray(sweep['v'], float)
    ref = np.abs(np.asarray(sweep['i_d'], float))
    b = sweep['bias']
    if sweep['sweep'] == 'Vd':
        got = np.array([np.asarray(e.i(np.array([x, b['Vg'], b['Vs'],
                                                 b['Vb']])), float)[0]
                        for x in v])
    else:
        got = np.array([np.asarray(e.i(np.array([b['Vd'], x, b['Vs'],
                                                 b['Vb']])), float)[0]
                        for x in v])
    m = ref > FLOOR
    return v[m], ref[m], np.abs(got[m]), kw


def main(argv=None):
    from pycircuit.utilities import spicecard
    import ngspice_ref

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pdk', default=ngspice_ref.DEFAULT_PDK)
    args = ap.parse_args(argv)
    if not os.path.isdir(args.pdk):
        ap.error('PDK not found at %s' % args.pdk)

    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    cm.default_toolkit = numeric

    with open(os.path.normpath(REF)) as fh:
        ref = json.load(fh)['sweeps']
    deck = spicecard.read(os.path.join(args.pdk, 'models',
                                       'cornerMOSlv.lib'),
                          section='mos_tt')

    print()
    print('Our surface-potential core vs IHP\'s own PSP103, drain current.')
    print('Nothing tuned: parameters come from the card through the')
    print('scaling layer.  Ratio = ours / PSP103.')
    print()
    print('%-20s %7s %7s %9s %9s %9s %9s' %
          ('sweep', 'W (um)', 'L (um)', 'median', 'min', 'max', 'spread'))
    for name in SWEEPS:
        if name not in ref:
            continue
        _, r, g, _ = compare(deck, ref[name])
        if not len(r):
            continue
        ratio = g / r
        print('%-20s %7.1f %7.2f %9.3f %9.3f %9.3f %9.3f'
              % (name, ref[name]['w'] * 1e6, ref[name]['l'] * 1e6,
                 np.median(ratio), ratio.min(), ratio.max(),
                 ratio.max() / ratio.min() - 1.0))
        sys.stdout.flush()

    print()
    print('The long device used to be the one to read, on the grounds')
    print('that short-channel physics was omitted and mattered least')
    print('there.  That is no longer why: the short-channel block is')
    print('largely present now, and the two geometries agree, which is')
    print('worth more than either number.  What separates the sweeps now')
    print('is bias condition rather than geometry -- and the SPREAD')
    print('column is the one to read, because it is the bias dependence')
    print('that a missing term breaks.')
    print()

    ## The detail of the long device, where the comparison is meaningful.
    v, r, g, kw = compare(deck, ref['nmos_long_idvd'])
    print('nmos_long_idvd in detail (W = 10 um, L = 1 um, Vg = 1.2 V):')
    print('%8s %14s %14s %8s' % ('Vd (V)', 'PSP103 (A)', 'ours (A)',
                                 'ratio'))
    for k in range(0, len(v), max(1, len(v) // 10)):
        print('%8.2f %14.5e %14.5e %8.3f' % (v[k], r[k], g[k], g[k] / r[k]))
    print()
    print('scaled parameters used:')
    for k in sorted(kw):
        print('   %-8s = %.6g' % (k, kw[k]))


if __name__ == '__main__':
    warnings.simplefilter('ignore')
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
