"""Generate a golden PSP103 I-V reference from the IHP Open PDK.

pycircuit's PSP work needs something to be right *against*.  The IHP
SG13G2 open PDK ships the full PSP 103.8.2 Verilog-A source, compiled
OSDI binaries, and extracted 359-parameter model cards, so a real
industry-standard surface-potential MOSFET can be swept in ngspice and
the curves kept.

This script writes ``pycircuit/circuit/tests/data/psp103_ihp_iv.json``.
The tests read that file; they never invoke ngspice.  That split is
deliberate -- the reference is a *fixed artifact* checked into the repo,
so a test failure means pycircuit changed, not that someone's PDK
checkout moved.  Re-run this only when the reference itself should
change, and say so in the commit.

Requires: ngspice with OSDI support (>= 42) and the IHP PDK.  Point
``--pdk`` at ``.../ihp-sg13g2/libs.tech/ngspice`` if it is not in the
default place.

Run:  python benchmarks/psp_reference.py [--pdk PATH] [--out PATH]
"""

import argparse
import json
import os
import sys

import ngspice_ref


DEFAULT_OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', 'pycircuit', 'circuit', 'tests', 'data', 'psp103_ihp_iv.json')


#: The sweeps.  Chosen to cover the three regimes a compact model has to
#: get right, at two geometries so the geometry-scaling layer is
#: exercised rather than one binned corner.
SWEEPS = [
    dict(name='nmos_idvd_vg1p2', device='sg13_lv_nmos', w=1e-6, l=0.13e-6,
         sweep='Vd', start=0.0, stop=1.5, step=0.05,
         bias=dict(Vg=1.2, Vs=0.0, Vb=0.0),
         note='output characteristic, strong inversion'),
    dict(name='nmos_idvd_vg0p6', device='sg13_lv_nmos', w=1e-6, l=0.13e-6,
         sweep='Vd', start=0.0, stop=1.5, step=0.05,
         bias=dict(Vg=0.6, Vs=0.0, Vb=0.0),
         note='output characteristic, near threshold'),
    dict(name='nmos_idvg_vd0p05', device='sg13_lv_nmos', w=1e-6, l=0.13e-6,
         sweep='Vg', start=0.0, stop=1.5, step=0.025,
         bias=dict(Vd=0.05, Vs=0.0, Vb=0.0),
         note='transfer characteristic, linear region -- subthreshold '
              'slope lives here'),
    dict(name='nmos_idvg_vd1p2', device='sg13_lv_nmos', w=1e-6, l=0.13e-6,
         sweep='Vg', start=0.0, stop=1.5, step=0.025,
         bias=dict(Vd=1.2, Vs=0.0, Vb=0.0),
         note='transfer characteristic, saturation -- DIBL lives here'),
    dict(name='nmos_idvg_vb_m1', device='sg13_lv_nmos', w=1e-6, l=0.13e-6,
         sweep='Vg', start=0.0, stop=1.5, step=0.025,
         bias=dict(Vd=0.05, Vs=0.0, Vb=-1.0),
         note='body effect'),
    dict(name='nmos_long_idvd', device='sg13_lv_nmos', w=10e-6, l=1e-6,
         sweep='Vd', start=0.0, stop=1.5, step=0.05,
         bias=dict(Vg=1.2, Vs=0.0, Vb=0.0),
         note='long/wide device -- geometry scaling, little short-channel'),
    dict(name='nmos_long_idvg', device='sg13_lv_nmos', w=10e-6, l=1e-6,
         sweep='Vg', start=0.0, stop=1.5, step=0.025,
         bias=dict(Vd=0.05, Vs=0.0, Vb=0.0),
         note='transfer characteristic of the LONG device -- the one '
              'diagnostic that separates a gain error (ratio flat in Vg) '
              'from a threshold or charge error (ratio varying with Vg), '
              'on the geometry where short-channel physics does not '
              'confound it'),
    ## THE BODY-BIASED LONG DEVICE.  Added because `XCOR` -- a body-bias
    ## mobility correction -- is exactly 1 at Vsb = 0 and therefore
    ## untouched by every other sweep here, while scaling to zero on the
    ## short geometry the one existing body-biased sweep uses.  It was a
    ## term no measurement could see.  This is the measurement.
    dict(name='nmos_long_idvg_vb_m1', device='sg13_lv_nmos', w=10e-6,
         l=1e-6, sweep='Vg', start=0.0, stop=1.5, step=0.025,
         bias=dict(Vd=0.05, Vs=0.0, Vb=-1.0),
         note='long device with reverse body bias -- the only sweep that '
              'exercises the body-bias mobility correction, which is '
              'identically 1 at Vsb = 0 and scales to zero at 0.13 um'),
    ## THE BODY-BIAS LADDER.  Two body biases cannot separate the two
    ## things that live in a body-bias difference: a body-factor error
    ## grows as `sqrt(phib + Vsb)` and the `XCOR` mobility correction as
    ## the saturating `(1 + 0.2*XCOR*Vsb)/(1 + XCOR*Vsb)`.  Both fit two
    ## points equally well and neither fits five.  Added because the
    ## n-channel subthreshold threshold offset grows with body bias and
    ## the existing pair of sweeps could not say which term it is.
] + [
    dict(name='nmos_long_idvg_vb%s' % ('%.2f' % vb).replace('-', 'm')
                                                  .replace('.', 'p'),
         device='sg13_lv_nmos', w=10e-6, l=1e-6, sweep='Vg',
         ## Starts at ZERO, and that is a limit rather than a choice.
         ## Below Vg = 0 this device's reference current stops being
         ## diffusion subthreshold and becomes junction and GIDL
         ## LEAKAGE, which PSP models and this element does not: the
         ## curve goes flat, `d ln I / dVg` goes to zero, and the
         ## implied-threshold measurement divides by it.  Tried at
         ## -0.4 V and the numbers were nonsense.  So the clean window
         ## is Vg >= 0, which on this geometry is about three decades.
         start=0.0, stop=1.0, step=0.025,
         bias=dict(Vd=0.05, Vs=0.0, Vb=vb),
         note='body-bias ladder: separates a sqrt(phib+Vsb) body-factor '
              'error from the saturating XCOR mobility correction')
    for vb in (0.0, -0.2, -0.4, -0.7, -1.0, -1.5)
] + [

    dict(name='pmos_idvd_vg1p2', device='sg13_lv_pmos', w=1e-6, l=0.13e-6,
         sweep='Vd', start=-1.5, stop=0.0, step=0.05,
         bias=dict(Vg=-1.2, Vs=0.0, Vb=0.0),
         note='p-channel output characteristic'),

    ## The p-channel LONG device.  The one 130 nm p-channel curve above
    ## was the least favourable geometry to judge a long-channel core on
    ## -- it is where the omitted short-channel physics dominates, so a
    ## ratio taken there says little about the core.  These mirror the
    ## n-channel set the whole measurement programme is built on.
    dict(name='pmos_long_idvd', device='sg13_lv_pmos', w=10e-6, l=1e-6,
         sweep='Vd', start=-1.5, stop=0.0, step=0.05,
         bias=dict(Vg=-1.2, Vs=0.0, Vb=0.0),
         note='p-channel long/wide device -- geometry scaling, little '
              'short-channel'),
    dict(name='pmos_long_idvg', device='sg13_lv_pmos', w=10e-6, l=1e-6,
         sweep='Vg', start=-1.5, stop=0.0, step=0.025,
         bias=dict(Vd=-0.05, Vs=0.0, Vb=0.0),
         note='p-channel transfer characteristic of the LONG device -- '
              'separates a gain error from a threshold error on the '
              'geometry where short-channel physics does not confound '
              'it'),
    dict(name='pmos_idvg_vd1p2', device='sg13_lv_pmos', w=1e-6, l=0.13e-6,
         sweep='Vg', start=-1.5, stop=0.0, step=0.025,
         bias=dict(Vd=-1.2, Vs=0.0, Vb=0.0),
         note='p-channel transfer characteristic, saturation'),
    dict(name='pmos_idvg_vd0p05', device='sg13_lv_pmos', w=1e-6,
         l=0.13e-6, sweep='Vg', start=-1.5, stop=0.0, step=0.025,
         bias=dict(Vd=-0.05, Vs=0.0, Vb=0.0),
         note='p-channel transfer characteristic, linear region -- '
              'subthreshold slope, and where the hole-specific '
              'effective-field weighting shows'),
]


def _deck(pdk, spec, out):
    bias = '\n'.join('V%s %s 0 dc %g' % (k[1:], k[1:], v)
                     for k, v in sorted(spec['bias'].items()))
    sweep_node = spec['sweep'][1:]
    return """* PSP103 reference: {name}
.lib {pdk}/models/cornerMOSlv.lib mos_tt

* THE SOLVER TOLERANCES ARE PART OF THE REFERENCE, AND THE DEFAULTS
* ARE NOT GOOD ENOUGH.
*
* A `dc` sweep seeds each point from the previous one, so Newton needs
* few iterations and stops as soon as `reltol*|i| + abstol` is met.  At
* ngspice's defaults that leaves up to 9.6e-4 of relative error on
* currents of order 1e-5 A -- point to point, so it reads as a KINK in
* the curve rather than as an offset, and it was the dominant error in
* this file at every bias above the leakage floor.
*
* WHICH KNOB MATTERS DEPENDS ON THE SWEEP, and measuring it on one
* sweep type and generalising gave the wrong answer once already:
*
*   sweep   default   abstol only   reltol only   both
*   idvg    6.41e-4      4.30e-8       4.30e-8    4.30e-8
*   idvd    7.97e-4      7.97e-4       6.61e-9    6.61e-9
*
* `reltol` is the one that matters and it fixes BOTH.  `abstol` fixes
* the transfer sweeps only: on an `idvd` sweep in saturation the current
* barely moves from point to point, so `reltol*|i|` is the binding term
* and no `abstol` can help.  `vntol` and `gmin` do nothing either way.
*
* `abstol` is kept as well because it costs nothing and guards the
* low-current end, where `reltol*|i|` vanishes.
*
* Both are far inside convergence.  `reltol` is flat from 1e-5 to 1e-10
* (and ngspice fails to converge on the `idvd` sweep by 1e-11, so 1e-6
* is deliberately mid-plateau rather than as tight as possible);
* `abstol` is bit-identical from 1e-14 to 1e-18.
*
* The `op` decks below do NOT need this and do not have it: a single
* operating point iterates to convergence with margin, and all 271 of
* its outputs are bit-identical either way.  It is specifically the
* SWEEP that stops early -- which is why regenerating moves every sweep
* and not one bit of the `lp_*`, threshold, capacitance or noise data.
*
* This was invisible for as long as the model was further out than
* 1e-3.  It became the limit on what could be measured once the model
* reached 1e-4, at which point the reference's own convergence was the
* floor rather than its physics.
.options reltol=1e-6 abstol=1e-15

{bias}
V{sn} {sn} 0 dc 0
X1 d g s b {device} w={w:g} l={l:g} ng=1 m=1

.control
osdi {pdk}/osdi/psp103.osdi
dc {sweep} {start:g} {stop:g} {step:g}
wrdata {out} v({sn}) i(Vd) i(Vg) i(Vb)
.endc
.end
""".format(name=spec['name'], pdk=pdk, bias=bias, sn=sweep_node,
           device=spec['device'], w=spec['w'], l=spec['l'],
           sweep=spec['sweep'], start=spec['start'], stop=spec['stop'],
           step=spec['step'], out=out)


def _run(pdk, spec):
    """Run one sweep and return the terminal currents.

    The deck records the DRAIN current always, never "the current in the
    source being swept": for a Vg sweep that would be gate leakage, four
    orders of magnitude below the quantity of interest and entirely
    plausible-looking until plotted.
    """
    ## Four vectors: the swept node voltage, then the three currents.
    ## `_deck` is handed a literal '{out}' so the placeholder survives its
    ## own `.format` and the harness fills in the data file.
    _sweep, (v, i_d, i_g, i_b) = ngspice_ref.run(
        _deck(pdk, spec, '{out}'), nvectors=4, label=spec['name'])
    ## The forced-voltage source current is the NEGATIVE of the current
    ## into the device terminal, so flip it -- what the reference stores
    ## is terminal current, the convention every consumer wants.
    return dict(v=v,
                i_d=[-v for v in i_d],
                i_g=[-v for v in i_g],
                i_b=[-v for v in i_b])


#: Geometries at which to record PSP's own internal threshold voltage.
#:
#: `vth` is an operating-point output of the model itself, so recording
#: it turns "is our scaling right?" from an inference about currents into
#: a direct comparison.  Note it is PSP's OWN definition of threshold and
#: not `vfb + phib + gamma*sqrt(phib)`; absolute values are not
#: comparable across models, but the SHIFT between geometries is.
VTH_GEOMETRIES = [
    dict(name='long', w=10e-6, l=1e-6),
    dict(name='mid', w=1e-6, l=0.5e-6),
    dict(name='short', w=1e-6, l=0.13e-6),
    dict(name='wide_short', w=10e-6, l=0.13e-6),
]

VTH_DECK = """* PSP103 internal vth: {name}
.lib {pdk}/models/cornerMOSlv.lib mos_tt
.temp {tempc:g}
Vd d 0 dc {vd:g}
Vg g 0 dc {vg:g}
Vs s 0 dc 0
Vb b 0 dc {vb:g}
X1 d g s b {device} w={w:g} l={l:g} ng=1 m=1
.control
osdi {pdk}/osdi/psp103.osdi
op
show all
.endc
.end
"""


#: Small-signal quantities and the charge model, at a grid of bias
#: points.  Recording these opens a dimension the currents cannot reach.
#:
#: `gm`/`gds`/`gmb` are DERIVATIVES, and a derivative localises an error
#: that an integrated current smears out: a drain current can be right
#: at every point of a sweep while its slope is wrong in the middle,
#: and the residual left on the n-channel has exactly that shape.
#:
#: The capacitances matter more.  Our charge model has been validated by
#: CONSTRUCTION -- conservation exact, source/drain swap exact,
#: Ward-Dutton partition reproducing the textbook 0.5 -> 0.4 -- and
#: never once against the vendor.  Those properties say the bookkeeping
#: is right; they say nothing about whether the charge itself is.  A
#: model whose charges are wrong is unusable in transient and AC however
#: good its DC is.
#:
#: One caveat governs the comparison: PSP's terminal capacitances
#: include overlap and fringe contributions this core does not model.
#: On the 0.13 um device that is about a fifth of `cgg`; on the 10 um /
#: 1 um device it is a few percent, which is why the long geometry is
#: the one to read here.
OP_OUTPUTS = ('ids', 'gm', 'gds', 'gmb', 'vth',
              ## Channel noise: the white drain density, the flicker
              ## density at 1 Hz, and the corner between them.  All three
              ## are bias-dependent, so they belong on the grid rather
              ## than in the scaled-parameter list -- and `fknee` is
              ## worth recording separately from the two densities
              ## because it is their RATIO, so it checks the two against
              ## each other rather than both against the current they
              ## share.
              'sid', 'sfl', 'fknee',
              ## INDUCED GATE NOISE.  `sig` is the gate density at 1 kHz
              ## -- a frequency, not a 1 Hz reference, because the
              ## density grows as f^2 and is meaningless extrapolated to
              ## DC.  `cigid` is the correlation with the drain, which
              ## is the whole reason the term is not just another
              ## independent source: an implementation can get `sig`
              ## right with `cigid` wrong, and the pair is what pins the
              ## coupling rather than only its magnitude.
              'sig', 'cigid',
              ## The JUNCTION, per component.  `cjs`/`cjd` are the
              ## totals and `cjsbot`/`cjssti`/`cjsgat` the three parts
              ## `SWJUNCAP = 3` switches on -- recording them separately
              ## is what turns a junction implementation from one number
              ## to check into three, which is the difference between
              ## finding a wrong component and knowing only that
              ## something is wrong.  `ijs`/`ijd` are the currents, and
              ## are recorded so their irrelevance stays measured rather
              ## than assumed.
              'cjs', 'cjd', 'cjsbot', 'cjssti', 'cjsgat',
              'ijs', 'ijd',
              ## `cgsol`/`cgdol` are the OVERLAP capacitances, and PSP
              ## reports them SEPARATELY from `cgg`.  Recording them is
              ## what settles a question that is otherwise easy to get
              ## backwards: `cgg` is the INTRINSIC capacitance, so a
              ## comparison against an intrinsic-only model is
              ## like-for-like and a residual in it is NOT overlap.
              'cgsol', 'cgdol',
              'cgg', 'cgs', 'cgd', 'cgb',
              'cdg', 'cdd', 'cds', 'cdb',
              'csg', 'csd', 'css', 'csb',
              'cbg', 'cbd', 'cbs', 'cbb')

#: Bias grid.  Chosen to span the regimes the residual moves through --
#: near threshold, moderate and strong inversion, at linear and
#: saturated drain bias -- plus one reverse-body point.
OP_BIASES = [
    dict(vg=0.4, vd=0.05, vb=0.0), dict(vg=0.4, vd=1.2, vb=0.0),
    dict(vg=0.8, vd=0.05, vb=0.0), dict(vg=0.8, vd=0.6, vb=0.0),
    dict(vg=0.8, vd=1.2, vb=0.0),
    dict(vg=1.2, vd=0.05, vb=0.0), dict(vg=1.2, vd=0.6, vb=0.0),
    dict(vg=1.2, vd=1.2, vb=0.0),
    dict(vg=1.2, vd=0.6, vb=-1.0),
]

OP_GEOMETRIES = [
    dict(name='long', w=10e-6, l=1e-6),
    dict(name='short', w=1e-6, l=0.13e-6),
]


#: PSP's own SCALED parameters, exposed as `lp_*` operating-point
#: outputs.  Recording these turns the geometry-scaling layer from
#: something inferred out of currents into something checkable term by
#: term -- a current can be right for compensating reasons, a scaled
#: parameter cannot.  This is how the `GPE` width-adjustment bug was
#: found: `lp_betn` was 12% off while every other term matched exactly.
LP_PARAMS = ('vfb', 'tox', 'neff', 'dphib', 'ct', 'cf', 'cfb', 'betn',
             'mue', 'themu', 'cs', 'thecs', 'rs', 'np',
             ## The saturation-voltage group.  `ax` is here because its
             ## scaling has a FLOOR (`PSP103_scaling.include:743` clips
             ## it at 2) that the card gives no hint of: `AXO/(1+AXL*iLE)`
             ## comes to 0.88 at 0.13 um, and an exponent below one makes
             ## the drain-voltage limiter bite far below saturation.  The
             ## only way to see the floor is to ask PSP what it used.
             ## `thesatb` and `thesatg` are recorded because they are
             ## NOT zero on this card and we do not model them yet -- a
             ## known gap is better written down than remembered.
             'ax', 'thesat', 'thesatb', 'thesatg', 'alp', 'alp1', 'alp2',
             ## Everything else the element takes and PSP will tell us.
             ## The point of this list is to leave NOTHING the model uses
             ## unchecked: three separate terms have now been found wrong
             ## or missing because a coefficient was zero on the one card
             ## they were tested against, and a parameter that is never
             ## compared is a parameter that is being assumed.
             ## `xcor` earns its place specially -- it is a body-bias
             ## mobility correction that is EXACTLY 1 at Vsb = 0, so it
             ## is invisible on every sweep but the body-biased one.
             'vp', 'rsg', 'rsb', 'dnsub', 'vnsub', 'nslp', 'cfd',
             'feta', 'xcor',
             ## `cox` is the TOTAL oxide capacitance PSP builds the
             ## charge model on, and it is not the geometric one: it
             ## uses the CV effective dimensions (`scaling:38-39, 359`),
             ## which the card shifts by `DLQ` and `DWQ`.  Recording it
             ## turns "our capacitances are 24% high" into two separate
             ## checkable numbers.
             'cox')


def _op_outputs(pdk, spec, names):
    """Named operating-point outputs of the OSDI model, via `show`."""
    import re
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.sp', delete=False) as fh:
        fh.write(VTH_DECK.format(pdk=pdk, **spec))
        path = fh.name
    try:
        out = subprocess.run(['ngspice', '-b', path], capture_output=True,
                             text=True, timeout=300).stdout
    finally:
        os.unlink(path)
    got = {}
    for n in names:
        m = re.search(r'^\s*%s\s+([-\d.eE+]+)\s*$' % re.escape(n),
                      out, re.M)
        if m:
            got[n] = float(m.group(1))
    return got


def _op_grid(pdk, device, sign):
    """Small-signal outputs and charges over `OP_BIASES`, per geometry.

    One `op` per bias point.  Slow -- a couple of dozen ngspice runs --
    and it is recorded data, so it is paid once.
    """
    out = {}
    for geom in OP_GEOMETRIES:
        pts = []
        for b in OP_BIASES:
            spec = dict(geom, device=device, tempc=27.0,
                        vd=sign * b['vd'], vg=sign * b['vg'],
                        vb=sign * b['vb'])
            got = _op_outputs(pdk, spec, OP_OUTPUTS)
            if not got:
                continue
            pts.append(dict(got, vg=sign * b['vg'], vd=sign * b['vd'],
                            vb=sign * b['vb']))
        out[geom['name']] = dict(w=geom['w'], l=geom['l'], points=pts)
    return out


def _vth(pdk, spec):
    """PSP's own `vth` operating-point output, via `show`."""
    import re
    import subprocess
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.sp', delete=False) as fh:
        fh.write(VTH_DECK.format(pdk=pdk, **spec))
        path = fh.name
    try:
        out = subprocess.run(['ngspice', '-b', path], capture_output=True,
                             text=True, timeout=300).stdout
    finally:
        os.unlink(path)
    m = re.search(r'\bvth\s+([-\d.eE+]+)', out)
    if not m:
        raise RuntimeError('no vth in ngspice output for %s' % spec['name'])
    return float(m.group(1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pdk', default=ngspice_ref.DEFAULT_PDK)
    ap.add_argument('--out', default=os.path.normpath(DEFAULT_OUT))
    args = ap.parse_args(argv)

    if not os.path.isdir(args.pdk):
        ap.error('PDK not found at %s' % args.pdk)

    ver = ngspice_ref.version()

    data = dict(
        source='IHP-Open-PDK SG13G2, PSP 103.8.2 via OSDI',
        simulator=ver,
        corner='mos_tt',
        note=('Terminal currents, sign convention "into the terminal". '
              'Generated by benchmarks/psp_reference.py; do not edit.'),
        sweeps={},
    )
    for spec in SWEEPS:
        sys.stdout.write('%-22s ' % spec['name'])
        sys.stdout.flush()
        res = _run(args.pdk, spec)
        res.update({k: spec[k] for k in
                    ('device', 'w', 'l', 'sweep', 'bias', 'note')})
        data['sweeps'][spec['name']] = res
        print('%d points, max |Id| = %.4g A'
              % (len(res['v']), max(abs(v) for v in res['i_d'])))

    ## PSP's own threshold voltage, per geometry.
    data['vth'] = {}
    data['scaled'] = {}
    ## The p-channel set goes under its OWN top-level keys rather than
    ## into `scaled`/`vth`.  A test asserts the exact key set of
    ## `scaled`, which is what makes it a completeness check rather than
    ## a lookup; mixing channel types in would quietly destroy that.
    data['vth_pmos'] = {}
    data['scaled_pmos'] = {}
    for dev, bias, vkey, skey in (
            ('sg13_lv_nmos', (0.05, 1.2), 'vth', 'scaled'),
            ('sg13_lv_pmos', (-0.05, -1.2), 'vth_pmos', 'scaled_pmos')):
        for geom in VTH_GEOMETRIES:
            spec = dict(geom, device=dev, vd=bias[0], vg=bias[1],
                        vb=0.0, tempc=27.0)
            v = _vth(args.pdk, spec)
            data[vkey][spec['name']] = dict(w=spec['w'], l=spec['l'],
                                            vth=v)
            lp = _op_outputs(args.pdk, spec,
                             ['lp_' + n for n in LP_PARAMS])
            data[skey][spec['name']] = dict(
                w=spec['w'], l=spec['l'],
                **{k[3:]: val for k, val in lp.items()})
            print('vth %-6s %-12s W=%-6.4g L=%-8.4g %9.6f V  (%d params)'
                  % (dev[-4:], spec['name'], spec['w'], spec['l'], v,
                     len(lp)))

    ## THE SCALED PARAMETERS AT A SECOND TEMPERATURE.
    ##
    ## The card specifies everything at `TR = 27 C` and PSP scales it
    ## from there; 100 C is 73 K away, far enough that a wrong exponent
    ## or a flipped sign cannot hide.  Recording it turns the whole
    ## temperature layer from something argued to something checked --
    ## the same move the `lp_*` comparison made for geometry.
    data['scaled_hot'] = {}
    data['scaled_hot_pmos'] = {}
    for dev, key in (('sg13_lv_nmos', 'scaled_hot'),
                     ('sg13_lv_pmos', 'scaled_hot_pmos')):
        for geom in VTH_GEOMETRIES:
            spec = dict(geom, device=dev, vd=0.05, vg=1.2, vb=0.0,
                        tempc=100.0)
            lp = _op_outputs(args.pdk, spec,
                             ['lp_' + n for n in LP_PARAMS])
            data[key][geom['name']] = dict(
                w=geom['w'], l=geom['l'], tempc=100.0,
                **{k[3:]: val for k, val in lp.items()})
        print('%-16s %d geometries at 100 C' % (key, len(data[key])))

    ## Small-signal and charge outputs over a bias grid.
    for key, dev, sign in (('op', 'sg13_lv_nmos', 1.0),
                           ('op_pmos', 'sg13_lv_pmos', -1.0)):
        data[key] = _op_grid(args.pdk, dev, sign)
        n = sum(len(g['points']) for g in data[key].values())
        print('%-8s %d bias points over %d geometries'
              % (key, n, len(data[key])))

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    print('\nwrote %s (%d sweeps)' % (args.out, len(data['sweeps'])))


if __name__ == '__main__':
    main()
