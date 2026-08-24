"""Characterising a foundry MOSFET from its own model card.

Reads IHP's PSP103 card for the SG13G2 process, builds the device, and
measures it the way a bench would: a transfer curve, an output curve,
and the small-signal operating point.  Nothing is fitted -- every
parameter comes off the card through the scaling layer.
"""
import numpy as np

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd, psp_scaling, DCSweep, SubCircuit, VS
from pycircuit.circuit import defaultepar
from pycircuit.circuit.compact import PspMosLongChannel
from pycircuit.circuit.toolkit import numeric
from pycircuit.utilities import spicecard

PDK = '~/source/IHP-Open-PDK/ihp-sg13g2/libs.tech/ngspice/models'
W, L = 1e-6, 0.13e-6

## The card is specified at its own reference temperature, and the
## scaling layer must be told the SAME temperature the element will be
## evaluated at.  They are not linked, and disagreeing about it is a
## quiet way to be wrong.
cm.default_toolkit = numeric
defaultepar.T = 273.15 + 27.0


def device(**over):
    """One transistor, scaled from the card for this geometry."""
    import os
    deck = spicecard.read(os.path.join(os.path.expanduser(PDK),
                                       'cornerMOSlv.lib'),
                          section='mos_tt')
    card = deck.model_params('sg13g2_lv_nmos_psp', w=W, l=L,
                             ng=1, m=1, pre_layout=1)
    kw = psp_scaling.to_long_channel(card, w=W, l=L, T=defaultepar.T)
    kw.update(over)
    return PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                             cm.Node('b'), **kw), kw


## ---------------------------------------------------------------- 1.
## A transfer curve, through a real DC solve.
def transfer(vg, vd=1.2):
    _, kw = device()
    c = SubCircuit()
    d, g = c.add_node('d'), c.add_node('g')
    c['Vd'] = VS(d, gnd, v=vd)
    c['Vg'] = VS(g, gnd, v=0.0)
    c['M1'] = PspMosLongChannel(d, g, gnd, gnd, **kw)
    res = DCSweep(c).solve('Vg', 'v', vg)
    return -np.asarray(res.i('Vd.plus'), float)


## ---------------------------------------------------------------- 2.
## An output curve, evaluated directly on the element.  `e.bias()`
## builds the unknown vector, `e.i()` returns the terminal currents in
## the order the terminals were given: d, g, s, b.
def output(vd, vg=1.2):
    e, _ = device()
    e.update_iparv()
    return np.array([np.asarray(e.i(e.bias(v, vg, 0.0, 0.0)),
                                float)[0] for v in vd])


## ---------------------------------------------------------------- 3.
## The small-signal operating point.  `G = di/dx` is the exact Jacobian
## the DSL generated, so the transconductances are read straight off it
## rather than differenced.
##
## LOOK THE NODES UP BY NAME.  The element owns internal nodes -- a gate
## node behind the gate resistance, and the auxiliary node the induced
## gate noise lives on -- so the unknown vector is longer than the four
## terminals and `gm` is the derivative with respect to the INTERNAL
## gate.  Indexing `G[0, 1]` returns exactly zero, because the drain
## current does not depend on the external gate at all: `rg` separates
## them.
def operating_point(vd=1.2, vg=1.2):
    e, _ = device()
    e.update_iparv()
    n = {str(node): k for k, node in enumerate(e.nodes)}
    x = e.bias(vd, vg, 0.0, 0.0)
    i = np.asarray(e.i(x), float)
    G = np.asarray(e.G(x), float)
    d, gi, b = n['d'], n.get('gi', n['g']), n['b']
    return dict(ids=i[d], gm=G[d, gi], gds=G[d, d], gmb=G[d, b])


if __name__ == '__main__':
    vg = np.linspace(0.4, 1.2, 5)
    print('transfer  Vd=1.2:', ' '.join('%.4e' % v for v in transfer(vg)))
    vd = np.linspace(0.2, 1.2, 5)
    print('output    Vg=1.2:', ' '.join('%.4e' % v for v in output(vd)))
    op = operating_point()
    print('op        Ids=%.5e  gm=%.5e  gds=%.5e  gmb=%.5e'
          % (op['ids'], op['gm'], op['gds'], op['gmb']))

    ## The committed reference records PSP's own operating point at this
    ## bias, so the example checks itself against the vendor.
    import json, os
    import pycircuit.circuit
    ref = json.load(open(os.path.join(
        os.path.dirname(pycircuit.circuit.__file__), 'tests', 'data',
        'psp103_ihp_iv.json')))['op']['short']
    pt = next(p for p in ref['points']
              if (p['vg'], p['vd'], p['vb']) == (1.2, 1.2, 0.0))
    print('psp103    Ids=%.5e  gm=%.5e  gds=%.5e  gmb=%.5e'
          % (pt['ids'], pt['gm'], pt['gds'], pt['gmb']))
