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
    dict(name='pmos_idvd_vg1p2', device='sg13_lv_pmos', w=1e-6, l=0.13e-6,
         sweep='Vd', start=-1.5, stop=0.0, step=0.05,
         bias=dict(Vg=-1.2, Vs=0.0, Vb=0.0),
         note='p-channel output characteristic'),
]


def _deck(pdk, spec, out):
    bias = '\n'.join('V%s %s 0 dc %g' % (k[1:], k[1:], v)
                     for k, v in sorted(spec['bias'].items()))
    sweep_node = spec['sweep'][1:]
    return """* PSP103 reference: {name}
.lib {pdk}/models/cornerMOSlv.lib mos_tt

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

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as fh:
        json.dump(data, fh, indent=1, sort_keys=True)
    print('\nwrote %s (%d sweeps)' % (args.out, len(data['sweeps'])))


if __name__ == '__main__':
    main()
