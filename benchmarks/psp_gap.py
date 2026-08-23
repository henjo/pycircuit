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

  * a ratio that is flat in Vds and near 1 on the long device means the
    threshold and gain factor are right;
  * a ratio that GROWS with current is series resistance -- the core has
    none, and PSP folds it into the mobility as an extra `Gmob` term;
  * a large ratio on the SHORT device only is the short-channel block:
    DIBL, channel-length modulation, and the much larger relative
    contribution of series resistance at 0.13 um.

What the core deliberately lacks, and therefore what this measures the
cost of: series resistance, channel-length modulation, DIBL and the rest
of the short-channel parameters, gate and junction leakage, impact
ionisation, overlap and fringe capacitances, and every temperature
coefficient.
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
    print('%-20s %7s %7s %9s %9s %9s' %
          ('sweep', 'W (um)', 'L (um)', 'median', 'min', 'max'))
    for name in SWEEPS:
        if name not in ref:
            continue
        _, r, g, _ = compare(deck, ref[name])
        if not len(r):
            continue
        ratio = g / r
        print('%-20s %7.1f %7.2f %9.3f %9.3f %9.3f'
              % (name, ref[name]['w'] * 1e6, ref[name]['l'] * 1e6,
                 np.median(ratio), ratio.min(), ratio.max()))
        sys.stdout.flush()

    print()
    print('The long device is the one to read: at L = 1 um the')
    print('short-channel physics the core omits matters least, so what is')
    print('left there is mostly the missing series resistance.  At')
    print('L = 0.13 um the omissions dominate and the ratio is not a')
    print('measure of the core at all.')
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
