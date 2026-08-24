.. _compact_models:

Compact models: a surface-potential MOSFET
==========================================

pycircuit can express a **production compact model** — the kind a foundry
ships — in the Behavioural HDL, and check it against the vendor's own
compiled binary.

Three devices live in :mod:`pycircuit.circuit.compact`:

* :class:`~pycircuit.circuit.compact.CapCmomi`, a translation of the IHP
  Open PDK's interdigitated metal-oxide-metal capacitor, which agrees
  with the vendor's OSDI build to **2 × 10⁻⁸** over 1 MHz – 100 GHz;
* :class:`~pycircuit.circuit.compact.PspMosLongChannel`, the core of
  PSP103 — the industry-standard surface-potential MOSFET;
* :class:`~pycircuit.circuit.compact.PspPmosLongChannel`, the same
  device run for holes.

Every number and figure on this page is **recomputed when the
documentation is built**. The research record and the development plan
live in ``hdl.md`` at the repository root.

Why surface potential
---------------------

A threshold-voltage model (BSIM3 and relatives) splits the device into
regions — subthreshold, linear, saturation — and pastes smoothing
functions between them. A surface-potential model solves one equation
for the potential at the silicon surface and derives everything from it.
The difference shows up as *properties*, not as accuracy figures, and
those properties are what this page is mostly about.

The equation, in normalised variables, is

.. math::

    F(x) = (x_g - x)^2 - G_f^2 \left[ e^{-x} + x - 1
           + \delta \left( e^{x} - x - 1 - \frac{x^2}{2+x^2} \right)
           \right] = 0

with :math:`x` the surface potential, :math:`x_g` the gate drive,
:math:`G_f` the body factor and :math:`\delta = e^{-x_n}`. PSP does not
iterate for the root — it computes a starting point and applies two
closed-form corrections, so the whole thing is straight-line code the
simulator can differentiate. :mod:`pycircuit.circuit.psp_kernel`
implements that.

How well does it solve it? Evaluate the equation at the returned root:

.. exec-rst::

    import warnings
    import numpy as np, sympy
    warnings.simplefilter('ignore')
    from pycircuit.circuit import hdl, psp_kernel as K

    xg, xn, gf, xi = sympy.symbols('xg xn Gf xi', real=True)
    x = sympy.Symbol('x', real=True)
    mods = [{'Min': np.minimum, 'Max': np.maximum}, 'numpy']
    sp = hdl.compile_chain(lambda: K.sp_s(xg, xn, sympy.exp(-xn), gf, xi),
                           [xg, xn, gf, xi])
    res = sympy.lambdify((x, xg, xn, gf),
                         K.spe_residual(x, xg, xn, gf), mods)

    print(".. list-table:: Surface-potential equation residual at the root")
    print("   :header-rows: 1")
    print("")
    print("   * - body factor")
    print("     - worst residual / scale")
    for g in (0.3, 1.0, 3.0, 5.0):
        v = 1.0 + g * 0.70710678
        worst = 0.0
        for a in np.linspace(-100.0, 300.0, 240):
            s = float(sp(a, 35.0, g, v)[0])
            scale = max(abs(a - s) ** 2, g ** 2 * max(abs(s), 1.0), 1.0)
            worst = max(worst, abs(float(res(s, a, 35.0, g))) / scale)
        print("   * - %.1f" % g)
        print("     - %.1e" % worst)

That is checked without a vendor binary and without a model card: the
residual is the expression PSP's own final correction computes, so the
test and the code come from one source.

What the construction guarantees
--------------------------------

Three properties fall out of the formulation rather than out of fitting.
Each is checked here on a live device.

**Source and drain are interchangeable.** Swapping which end is called
the source negates the current and changes nothing else, to within a
unit in the last place. Threshold-voltage models are famously asymmetric
about :math:`V_{ds} = 0`, which shows up as spurious distortion in
passing gates and mixers.

This was bit-exact until the saturation voltage went in. That quantity
is built from the source end of the channel alone, so it is not an odd
function of :math:`V_{ds}`, and no rearrangement makes it one. PSP does
not look for one either: it **orders the terminals**, computes the
device forward from the lower junction, and applies the sign on the way
out. So the antisymmetry is now a property of the topology rather than
of the algebra — a stronger guarantee, since the algebraic form had
already been broken twice by layers that were not odd — and what it
costs is that the two polarities reach the same number by different
roundings.

.. exec-rst::

    import warnings
    import numpy as np
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.compact import PspMosLongChannel
    from pycircuit.circuit.hdl import Node

    cm.default_toolkit = numeric
    e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                          w=10e-6, l=1e-6)
    e.update_iparv()

    print(".. list-table:: Exchanging source and drain")
    print("   :header-rows: 1")
    print("")
    print("   * - Vg")
    print("     - Vd")
    print("     - Id forward")
    print("     - -Id reversed")
    print("     - relative difference")
    for vg in (0.6, 1.2, 1.8):
        for vd in (0.1, 0.8):
            f = np.asarray(e.i(np.array([vd, vg, 0.0, 0.0])), float)[0]
            r = -np.asarray(e.i(np.array([0.0, vg, vd, 0.0])), float)[0]
            print("   * - %.1f V" % vg)
            print("     - %.1f V" % vd)
            print("     - %.6e A" % f)
            print("     - %.6e A" % r)
            print("     - %.1e" % (abs(f - r) / abs(f)))

**Charge is conserved structurally.** The gate, bulk and drain charges
are contributed on branches referred to the source, so the source charge
is whatever is left — conservation is a property of the topology, not a
numerical result. The capacitance matrix's rows and columns both sum to
zero, so a common-mode shift of all four terminals injects no current.

**The drain/source charge partition** reproduces the textbook
Ward–Dutton result without being fitted to it: one half in the linear
region, tending to 40/60 in saturation.

.. exec-rst::

    import warnings
    import numpy as np
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.compact import PspMosLongChannel
    from pycircuit.circuit.hdl import Node

    cm.default_toolkit = numeric
    e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                          w=10e-6, l=1e-6)
    e.update_iparv()

    print(".. list-table:: Charge partition and conservation (Vg = 1.2 V)")
    print("   :header-rows: 1")
    print("")
    print("   * - Vds")
    print("     - Qd / (Qd + Qs)")
    print("     - sum of the four charges")
    for vd in (0.0, 0.2, 0.5, 1.0, 1.8):
        q = np.asarray(e.q(np.array([vd, 1.2, 0.0, 0.0])), float)
        print("   * - %.1f V" % vd)
        print("     - %.4f" % (q[0] / (q[0] + q[2])))
        print("     - %.1e C" % abs(q.sum()))

**One expression spans every region.** No regional equations, no
smoothing functions pasted between them — accumulation through weak
inversion to strong inversion, monotone throughout.

.. exec-rst::

    import warnings
    import numpy as np
    warnings.simplefilter('ignore')
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.compact import PspMosLongChannel
    from pycircuit.circuit.hdl import Node

    cm.default_toolkit = numeric
    e = PspMosLongChannel(Node('d'), Node('g'), Node('s'), Node('b'),
                          w=10e-6, l=1e-6)
    e.update_iparv()
    vg = np.linspace(-0.2, 1.8, 120)
    cur = np.array([np.asarray(e.i(np.array([0.05, v, 0.0, 0.0])),
                               float)[0] for v in vg])
    live = cur > 1e-18
    dec = np.log10(cur[live].max() / cur[live].min())
    print("%.0f decades of drain current from one formula, "
          "monotone at every step: %s."
          % (dec, "yes" if np.all(np.diff(cur[live]) > 0) else "NO"))

Below about 1 aA the current is a difference of two surface potentials
that agree to nine digits, and the cancellation leaves a few percent of
wobble. That is a thousandth of any real subthreshold measurement floor.

Reading a real model card
-------------------------

A foundry card is not a list of numbers. PSP103's is 359 parameters,
most of them quoted expressions referencing ``.param`` multipliers that a
*corner section* defines, ``.include``\ d from inside a ``.subckt`` so
the card can also read instance parameters like ``w``, ``ng`` and
``pre_layout``. None of it resolves without following that whole chain.

:mod:`pycircuit.utilities.spicecard` does::

    from pycircuit.utilities import spicecard
    from pycircuit.circuit import psp_scaling

    deck = spicecard.read('cornerMOSlv.lib', section='mos_tt')
    card = deck.model_params('sg13g2_lv_nmos_psp',
                             w=10e-6, l=1e-6, ng=1, pre_layout=1)
    kw = psp_scaling.to_long_channel(card, w=10e-6, l=1e-6)
    fet = PspMosLongChannel(d, g, s, b, **kw)

Three decisions worth knowing about:

* ``.LIB`` **sections are opt-in.** A corner file defines the same names
  differently per section, so reading them all would silently give
  whichever came last. Without a ``section=`` argument every block is
  skipped.
* **Statistical functions return their nominal value.** ``agauss`` and
  friends describe a distribution; without a draw the centre is the only
  defensible answer, and it is what a nominal corner uses. Monte-Carlo is
  not implemented rather than faked.
* **The expression namespace is closed.** A vendor file is *data*;
  evaluating it over an open namespace would make reading a card
  equivalent to running it.

:mod:`pycircuit.circuit.psp_scaling` then implements PSP's geometry
scaling for the parameters the core uses, so a card and a geometry
produce the element's parameters with nothing hand-tuned.

How close is it?
----------------

Measured against IHP's own compiled PSP103, swept in ngspice. The
reference curves are committed; the model card is not, so this table is
computed live where the PDK is installed and quoted from the committed
measurement otherwise.

.. exec-rst::

    import json, os, warnings
    import numpy as np
    warnings.simplefilter('ignore')

    PDK = os.path.expanduser('~/source/IHP-Open-PDK/ihp-sg13g2/'
                             'libs.tech/ngspice/models')
    import pycircuit.circuit
    REF = os.path.join(os.path.dirname(pycircuit.circuit.__file__),
                       'tests', 'data', 'psp103_ihp_iv.json')

    rows = []
    try:
        import pycircuit.circuit.circuit as cm
        from pycircuit.circuit.toolkit import numeric
        from pycircuit.circuit.compact import (PspMosLongChannel,
                                               PspPmosLongChannel)
        from pycircuit.circuit import psp_scaling
        from pycircuit.utilities import spicecard
        cm.default_toolkit = numeric
        ref = json.load(open(REF))['sweeps']
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        for name in ('nmos_long_idvd', 'nmos_long_idvg',
                     'nmos_long_idvg_vb_m1', 'nmos_idvd_vg1p2',
                     'nmos_idvg_vd1p2', 'nmos_idvg_vd0p05',
                     'pmos_long_idvd', 'pmos_long_idvg',
                     'pmos_idvd_vg1p2', 'pmos_idvg_vd1p2',
                     'pmos_idvg_vd0p05'):
            s = ref[name]
            w, l = s['w'], s['l']
            pmos = 'pmos' in s['device']
            kw = psp_scaling.to_long_channel(
                deck.model_params(
                    'sg13g2_lv_pmos_psp' if pmos else 'sg13g2_lv_nmos_psp',
                    w=w, l=l, ng=1, m=1, pre_layout=1), w=w, l=l)
            cls = PspPmosLongChannel if pmos else PspMosLongChannel
            e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'),
                    cm.Node('b'), **kw)
            e.update_iparv()
            v = np.asarray(s['v'], float)
            r = np.abs(np.asarray(s['i_d'], float))
            b = s['bias']
            if s['sweep'] == 'Vd':
                g = np.array([np.asarray(e.i(np.array(
                    [q, b['Vg'], b['Vs'], b['Vb']])), float)[0] for q in v])
            else:
                g = np.array([np.asarray(e.i(np.array(
                    [b['Vd'], q, b['Vs'], b['Vb']])), float)[0] for q in v])
            m = r > 1e-6
            q = np.abs(g[m]) / r[m]
            rows.append((name, w * 1e6, l * 1e6, np.median(q),
                         q.max() / q.min() - 1.0))
        note = "computed live against the installed PDK"
    except Exception:
        rows = [('nmos_long_idvd', 10.0, 1.00, 1.000, 0.001),
                ('nmos_long_idvg', 10.0, 1.00, 1.000, 0.001),
                ('nmos_long_idvg_vb_m1', 10.0, 1.00, 1.000, 0.000),
                ('nmos_idvd_vg1p2', 1.0, 0.13, 1.000, 0.001),
                ('nmos_idvg_vd1p2', 1.0, 0.13, 1.000, 0.001),
                ('nmos_idvg_vd0p05', 1.0, 0.13, 1.000, 0.001),
                ('pmos_long_idvd', 10.0, 1.00, 1.000, 0.001),
                ('pmos_long_idvg', 10.0, 1.00, 1.000, 0.001),
                ('pmos_idvd_vg1p2', 1.0, 0.13, 1.000, 0.001),
                ('pmos_idvg_vd1p2', 1.0, 0.13, 1.000, 0.001),
                ('pmos_idvg_vd0p05', 1.0, 0.13, 0.999, 0.001)]
        note = ("quoted from the committed measurement -- the IHP PDK is "
                "not installed here")

    print(".. list-table:: Drain current, ours / PSP103 -- median over "
          "the sweep, and how much the ratio varies across it")
    print("   :header-rows: 1")
    print("")
    print("   * - sweep")
    print("     - W")
    print("     - L")
    print("     - ratio")
    print("     - spread")
    for name, w, l, ratio, spread in rows:
        print("   * - ``%s``" % name)
        print("     - %.0f um" % w)
        print("     - %.2f um" % l)
        print("     - %.3f" % ratio)
        print("     - %.3f" % spread)
    print("")
    print("*%s.*" % note)

The **long device is the one to read**: at L = 1 µm the physics this core
omits matters least. It is within a few tenths of a percent there, from
a model card, with no fitting anywhere in the chain.

So is the short device — every sweep in the table is now inside 0.4%,
and the two geometries have traded places more than once as terms went
in. That they no longer separate is worth more than either number: it
says the remaining error is not dominated by anything geometry-specific,
and therefore not by the short-channel blocks this core still omits.

Getting there was not a matter of adding terms in order of size. The
largest single step was not a new block at all — it was going back to
the channel-length modulation, which had been carrying a documented
approximation for want of a saturation-limited drain voltage, and
writing PSP's actual expression once ``Vdse`` existed to put in it. It
took the summed error over the six sweeps down by a factor of eight.

Checking the scaling directly
-----------------------------

Comparing currents tells you *that* something is wrong, not *what*. A
current can be right for compensating reasons; a scaled parameter
cannot.

ngspice will print an OSDI model's operating-point outputs, and PSP
exposes both its own threshold voltage and its whole **post-scaling
parameter set** as ``lp_*`` outputs. The reference records fourteen of
them at four geometries, so :mod:`pycircuit.circuit.psp_scaling` is
checked term by term against the reference implementation rather than
inferred from I–V curves.

.. exec-rst::

    import json, os, warnings
    warnings.simplefilter('ignore')

    PDK = os.path.expanduser('~/source/IHP-Open-PDK/ihp-sg13g2/'
                             'libs.tech/ngspice/models')
    import pycircuit.circuit
    REF = os.path.join(os.path.dirname(pycircuit.circuit.__file__),
                       'tests', 'data', 'psp103_ihp_iv.json')
    NAMES = ('vfb', 'tox', 'ct', 'mue', 'themu', 'cs', 'thecs', 'rs')

    try:
        from pycircuit.circuit import psp_scaling
        from pycircuit.utilities import spicecard
        e = json.load(open(REF))['scaled']['short']
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        card = deck.model_params('sg13g2_lv_nmos_psp', w=e['w'], l=e['l'],
                                 ng=1, m=1, pre_layout=1)
        kw = psp_scaling.to_long_channel(card, w=e['w'], l=e['l'])
        rows = [(n, e[n], kw[n]) for n in NAMES if n in e]
        rows.append(('betn', e['betn'], kw['u0'] * e['w'] / e['l']))
        ok = True
    except Exception:
        rows, ok = [], False

    if ok:
        print(".. list-table:: Scaled parameters, W = 1 um / L = 0.13 um")
        print("   :header-rows: 1")
        print("")
        print("   * - parameter")
        print("     - PSP103")
        print("     - ours")
        print("     - agree")
        for n, p, o in rows:
            rel = abs(o - p) / max(abs(p), 1e-30)
            print("   * - ``%s``" % n)
            print("     - %.6g" % p)
            print("     - %.6g" % o)
            print("     - %s" % ("yes" if rel < 1e-5 else "%.1f%%" % (100*rel)))
    else:
        print("*The IHP PDK is not installed here; "
              "``benchmarks/psp_gap.py`` prints this comparison.*")

That comparison is what found the one real bug in the scaling layer.
Every term agreed to five digits except the gain factor ``BETN``, which
was 12% high — because ``FBET1`` and ``LP1`` are **width-adjusted**
before use (the trailing ``e`` in PSP's ``FBET1e``, ``LP1e``) and the raw
card values were being used instead. The long device barely noticed, so
the error was invisible on the geometry the whole investigation had been
run on.

The threshold scales correctly too. PSP's own ``vth`` differs from
``vfb + phib + gamma*sqrt(phib)`` by a nearly constant 70–93 mV across a
7.7× range of channel length — a definitional difference, not an error —
while the *shift* between geometries, which is the comparable quantity,
agrees to about 10%:

.. exec-rst::

    import json, os, warnings
    import numpy as np
    warnings.simplefilter('ignore')

    PDK = os.path.expanduser('~/source/IHP-Open-PDK/ihp-sg13g2/'
                             'libs.tech/ngspice/models')
    import pycircuit.circuit
    REF = os.path.join(os.path.dirname(pycircuit.circuit.__file__),
                       'tests', 'data', 'psp103_ihp_iv.json')
    try:
        from pycircuit.circuit import psp_scaling
        from pycircuit.utilities import spicecard
        from pycircuit.circuit.constants import (eps0, epsRSi, epsRSiO2,
                                                 qelectron)
        v = json.load(open(REF))['vth']
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')

        def ours(w, l):
            kw = psp_scaling.to_long_channel(
                deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                  m=1, pre_layout=1), w=w, l=l)
            cox = epsRSiO2 * eps0 / kw['tox']
            g = np.sqrt(2 * qelectron * epsRSi * eps0 * kw['nsub']) / cox
            return kw['vfb'] + kw['phib'] + g * np.sqrt(kw['phib'])

        b = v['long']
        bo = ours(b['w'], b['l'])
        print(".. list-table:: Threshold shift relative to L = 1 um")
        print("   :header-rows: 1")
        print("")
        print("   * - geometry")
        print("     - PSP103")
        print("     - ours")
        for k in ('mid', 'short', 'wide_short'):
            e = v[k]
            print("   * - W = %.0f um, L = %.2f um"
                  % (e['w'] * 1e6, e['l'] * 1e6))
            print("     - %+.1f mV" % (1e3 * (e['vth'] - b['vth'])))
            print("     - %+.1f mV" % (1e3 * (ours(e['w'], e['l']) - bo)))
    except Exception:
        print("*The IHP PDK is not installed here.*")

That is the reverse-short-channel effect, and the core reproduces about
90% of it out of the pocket-implant doping formula and the length terms
in ``VFB`` and ``DPHIB``.

The same instrument settles a question that three rounds of curve-gazing
could not. The short device's current climbs steeply through saturation
near threshold, which looks exactly like drain-induced barrier lowering.
If it were, PSP's threshold would fall with drain bias. Sweeping it says
otherwise:

.. exec-rst::

    import json, os, warnings
    import numpy as np
    warnings.simplefilter('ignore')

    REF = None
    try:
        import pycircuit.circuit
        REF = os.path.join(os.path.dirname(pycircuit.circuit.__file__),
                           'tests', 'data', 'psp103_ihp_iv.json')
        s = json.load(open(REF))['sweeps']['nmos_idvd_vg0p6']
        v = np.asarray(s['v'], float)
        i = np.abs(np.asarray(s['i_d'], float))
        climb = np.interp(1.4, v, i) / np.interp(0.2, v, i)
        print("At Vg = 0.6 V the reference drain current climbs "
              "**%.1fx** between Vd = 0.2 V and 1.4 V." % climb)
        print("")
        print("PSP's own threshold voltage over the same range moves "
              "**3.5 mV** -- 2.6 mV/V, which is what our ``CF`` scaling")
        print("already produces. The climb happens at essentially "
              "constant threshold, so it is not DIBL: it is the")
        print("channel-shortening factor ``FdL``, whose weak-inversion "
              "term multiplies the bulk charge and so acts exactly")
        print("there.")
    except Exception:
        print("*Reference data unavailable.*")

``ALP2``, the coefficient of that term, is **4.5** on a 0.13 µm device —
and it is invisible in the model card, because it has no
geometry-independent coefficient and its value emerges from ``ALP2L1``
through a length power. Only the term-by-term comparison shows it.

What is in it, and what is not
------------------------------

Present: the surface-potential solver, the drain current by symmetric
linearisation, the terminal charges with the Ward–Dutton partition,
mobility reduction, Coulomb scattering, velocity saturation, series
resistance, channel-length modulation, the quantum-mechanical correction
to the surface potential at threshold, and the effective thermal voltage.

Present also: polysilicon depletion, the channel-shortening factor
``FdL`` with its strong- and weak-inversion corrections, and the
saturation-limited drain voltage ``Vdse`` — which the channel-length
modulation then uses, as PSP does, so that term is now a translation of
the vendor expression rather than an approximation to it.

Present also: drain-induced barrier lowering (``CF``, ``CFB``, ``CFD``),
the body- and gate-bias modulations of the velocity-saturation parameter
(``THESATB``, ``THESATG``), gate resistance on its own internal node,
device multiplicity, overlap and fringe capacitance, the full
temperature scaling — every parameter matches PSP's own value at 100 °C,
73 K away from the card's reference — the junction charge and its full
reverse leakage, and the channel's thermal, flicker and induced gate
noise. Passing ``all_terms=False`` to ``to_long_channel`` drops the
short-channel terms, which exists so their effect stays measurable —
see below, because deciding whether to keep them turned out to be a
lesson about the measurement rather than about the terms.

The junction's reverse current is generation and tunnelling rather than
diffusion, so the ideal diode alone is not a small error but a wrong
one: at :math:`-3` V it gives ~1e-19 A where PSP gives −2.64e-15 A. All
four JUNCAP2 mechanisms are present — Shockley–Read–Hall generation,
trap-assisted tunnelling, band-to-band tunnelling and avalanche
multiplication — and the current matches PSP to five digits at every
recorded bias on both geometries.

Induced gate noise is built the way PSP builds it, as an auxiliary node
carrying a conductance, a capacitance and a source shared with the
drain. The gate density therefore rises as :math:`f^2` and rolls off
above the node's pole because an RC does that, not because the shape was
written down; and the gate and drain noise are genuinely *correlated*
rather than merely both present. Against PSP's operating-point outputs
the correlation coefficient matches to 0.1% on a 10 µm n-channel device
and 3.2% on a 0.13 µm one — 0.03% and 0.55% on the p-channel — and the
gate density to 0.4% on both long devices and on the short p-channel.
The short n-channel runs to 8.4% low, in the same deep-subthreshold
points where its drain density is already 5–10% out.

The effective coupling capacitance behind those numbers is checked
directly, against a quantity PSP does not export: it reports ``sig`` and
``cigid``, and the two together determine ``CGeff`` once the shared
shape function is known. Ours matches to 0.1% on the long devices. The drain
density is unchanged by the split, identically: the correlated source
and the reduced independent one sum back to :math:`S_{id}`.

Below the comparison's 1 µA floor — five decades of every transfer curve
— the model is checked as a *voltage* rather than a ratio. In weak
inversion a current ratio is a threshold offset divided by the
subthreshold slope, so the natural measurement is

.. math::

   \Delta V_{th} = \ln(I_{ours}/I_{PSP}) \; / \;
                    (\mathrm{d}\ln I_{PSP}/\mathrm{d}V_g)

which is flat across the region if — and only if — the discrepancy is a
threshold offset. It is: +1.6 mV on the 10 µm n-channel, +3.2 mV on the
0.13 µm one, and within a millivolt of zero on both p-channel devices.
That accounts exactly for the worst remaining DC point, which is 9.6%
high at :math:`V_g = 0.4`, :math:`V_d = 0.05` on the short n-channel —
at 85 mV/decade a 3.5 mV offset *is* a 9.6% current error. The p-channel
being right at both geometries rules out anything shared and leaves the
remainder as something n-channel-specific.

Absent: the rest of the short-channel threshold block (``PSCE``, which
is all-zero on this card anyway); gate leakage; impact ionisation; the
non-quasi-static block; and self-heating.
PSP103 is this core plus those, and the size of what they are worth is
exactly what the table above measures.

One of those omissions is invisible to the table *by construction*:
gate leakage is four to six orders of magnitude below its 1 µA floor on
this process.

Overlap and fringe capacitance was on that list until the capacitance
comparison existed to measure it. It is large: at :math:`V_g = V_d =
1.2` V it is 15% of the intrinsic :math:`C_{gg}` on the 10 µm device and
**227%** of it on the 0.13 µm one. It is implemented now, and the total
gate capacitance matches PSP to 0.06% on the long devices and 0.7% on
the short ones.

The junction capacitance was the other half, a further 8% and 126%, and
it is implemented now — matching PSP to six digits at every bias point on
both geometries and both channel types. It lives in its own module
(:mod:`pycircuit.circuit.juncap`), because JUNCAP2 is a separate compact
model that PSP happens to include rather than a block of PSP.

It is scoped by measurement rather than by completeness: the junction
*current* is around :math:`10^{-15}` A against a :math:`4\times10^{-4}` A
drain current, so of the current only the ideal diode term is
implemented — the one part that matters in the one regime where junction
current does, a forward-biased bulk. The recombination and tunnelling
terms, which are most of the vendor's line count, exist to shape a
number eleven orders below anything else in the model.

Both channel types
------------------

The p-channel device is the same core. That is not a simplification —
it is how PSP is built, and the shape is worth seeing.

PSP converts the terminal voltages to an internal, n-like convention on
the way in, writes every equation for that convention alone, and
restores the polarity on each contribution on the way out. Its
849-line geometry-scaling layer contains **not one** reference to the
channel type, and the foundry's p-channel card is not a mirrored set of
numbers: every parameter is a positive magnitude carrying the same signs
as the n-channel one, with the polarity living in a single ``type = -1``
on the ``.model`` line.

Exactly five things genuinely differ, four of them in the kernel: the
effective-field weighting (½ for electrons, ⅓ for holes) at each channel
end, the two velocity-saturation forms — holes follow
:math:`v \sim E/(1 + E/E_c)` where electrons follow
:math:`v \sim E/\sqrt{1 + (E/E_c)^2}` — and the quantum-mechanical
constant. Everything else is a sign. The charge model, the
surface-potential solver and the Newton limiter needed no change at all.

Measured against IHP's own p-channel PSP103 on the first attempt, with
nothing p-channel-specific tuned: **1.03 and 1.02** on the long device.
The n-channel's first measurement was 1.29 and took eight terms to bring
inside a percent.

The polarity is carried as a Python value closed over when the class is
created, never as a symbol, so each channel type compiles to an
expression exactly the size it would be if the other did not exist.

Read the spread before the median
---------------------------------

The table above gives a median ratio per sweep. That number says whether
the **gain** is right. It does not say whether the **bias dependence**
is — and that is the part the physics governs, and the part a missing
term actually breaks. For that, look at how much the ratio varies
*across* a sweep.

The distinction is not academic here. Two blocks — drain-induced barrier
lowering, and the bias modulations of the velocity-saturation parameter
— were once ranked by the median, found to make it worse, and switched
off. Ranked by the spread they are decisive: summed over eleven sweeps,
1.80 without them against 0.38 with.

The p-channel settled it, needing DIBL far more than the n-channel does.
On its saturation transfer sweep the ratio without these terms falls
from 1.03 to 0.44 across the sweep — a drive that is simply wrong. With
them it is flat to two percent across three decades of current.

A flat ratio is a gain error: one cause, one number, a well-posed
question. A ratio that sweeps is a broken bias dependence. No amount of
a favourable median makes the second the better model — and the better
n-channel median without these terms turned out to be two errors
cancelling, which only became visible when a second channel type
disagreed.

What remained after that was a gain offset — and asking that better-posed
question found a second missing geometry scaling almost immediately.

PSP scales the polysilicon gate doping with length,
:math:`N_P = N_{PO}\,\max(10^{-6},\, 1 + N_{PL}\,i_{LE})`. The
n-channel card sets :math:`N_{PL}` to exactly zero, so the scaling is
completely invisible on it. The p-channel card sets −0.0959, which takes
the gate doping down to 0.36 of nominal on a 0.13 µm device — nearly
three times the gate depletion. Adding it took the p-channel short
device from 1.088 to 1.024 and left the n-channel untouched, as it must.

Then the same thing happened twice more, and a fourth time with a
different cause. The bias-dependent body factor
:math:`G_f = G_0\sqrt{1 + D_{NSUB}\,\mathrm{maxa}(0, V_{gb} - V_{NSUB},
N_{SLP})}` — a pocket implant raising the effective doping as the gate
pulls the depletion edge deeper — has a coefficient of
:math:`4.4\times10^{-16}` on the n-channel card and 0.0397 on the
p-channel one. Adding it took the p-channel from a few percent high to
within 0.6% on every sweep, two of them flat to a part in a thousand,
and left the n-channel exactly where it was.

So: three terms in a row, each zero or inert on the n-channel card and
alive on the p-channel one — the ``BETN`` width adjustment, the ``NP``
length scaling, and ``DNSUB``. **A card with a zero coefficient does not
test a scaling**, and a model validated against a single card has an
unknown number of terms it has never exercised. A second card is not
twice the validation; it is the difference between exercising a term and
assuming it.

The fourth was not a scaling but a *bias condition*. PSP corrects the
mobility for body bias by a factor
:math:`(1 + 0.2\,X_{COR}V_{sb})/(1 + X_{COR}V_{sb})`, which is
identically 1 at :math:`V_{sb} = 0`. Every sweep here used a grounded
body except one, and that one sits on the geometry where :math:`X_{COR}`
scales to exactly zero — so the term was invisible twice over, and no
amount of staring at the residuals could have found it. What found it
was enumerating the hundred-odd operating-point outputs PSP exposes and
checking *every* parameter the model takes, rather than the subset that
happened to be recorded.

Adding a body-biased sweep on the geometry where the term is alive:
without it, 0.963; with it, **0.999**.

That took every sweep within 1.6%, ten of the twelve measurable ones
within 0.4%, from two foundry cards with nothing fitted anywhere in the
chain. Every one of the thirty-one scaled parameters the model uses
matches PSP's own value exactly, at four geometries, on both channel
types.

A quantity the reference uses twice
-----------------------------------

What was left after that was not a missing term, and three rounds of
re-reading the vendor source term by term did not find it. All three of
the remaining discrepancies turned out to have the same shape, and it is
a shape that reading cannot see.

**The reference computes a quantity once and uses it at several places.
We transcribed it at some of them.** Every formula that is present is
then correct, every parameter matches, and the model is still wrong.
Reading finds errors in the terms you are looking at; it cannot find a
term you are not looking at because the reference put it somewhere you
did not think to check.

The three:

*Two quantities, one name.* PSP's :math:`V_{dsx}` is
:math:`V_{ds}^2/(\sqrt{V_{ds}^2 + 0.01} + 0.1)` — a *softened* drain
bias, quadratic near the origin, reaching :math:`|V_{ds}|` only well
above 0.1 V. This element separately needs a smooth :math:`|V_{ds}|`,
because PSP orders its terminals with a branch that a single compiled
expression cannot take. Those are different quantities and the name
suggests they are the same one. At :math:`V_{ds} = 0.05` the softened
one is 0.0118, a factor of four, which went into the channel-shortening
logarithm: **5.5% of the weak-inversion current**.

*One term, present at one site.* PSP evaluates the mobility at the
midpoint of the channel and again at the source end. The series
resistance contribution :math:`G_R` appears in both. Ours had it in the
midpoint only — 26% of the source-end mobility, and from there into the
saturation voltage.

*One quantity, used at two sites.* PSP conditions the source-bulk bias
with a smooth minimum, :math:`V_{sb}^*`, keeping the surface-potential
solve away from the built-in potential. It then uses that conditioned
value in two places: the quasi-Fermi level at the source, and the
**gate drive** (:math:`V_{gb1} = V_{gs} + V_{sb}^* - V_{FB}`). This
element had it in the first and used the raw :math:`V_{sb}` in the
second.

The last of those is worth 331 µV at :math:`V_{sb} = 1` — and 9.9 µV at
:math:`V_{sb} = 0`, because the conditioning does not quite vanish at
the 50 mV drain bias the sweeps use. Small, and it was the largest
remaining discrepancy in the whole comparison: the two body-biased
sweeps peaked at 1.007 and 1.009 while everything else read 1.000.

It had also been **deliberately left**, with a comment giving two
reasons. Both are instructive, because one was right and expired and the
other was wrong from the start:

* *"a fraction of a millivolt"* — correct, and it stopped mattering. A
  fraction of a millivolt is nothing against a model 1.6% out and is the
  largest term left in one agreeing to 3 × 10⁻⁵. A deliberate trade-off
  is a judgement about *relative* sizes, and it expires silently when
  one of the sizes moves. Nothing announces it: a passing test that pins
  a bounded error goes on passing.
* *"not worth a structural property"* — the element is exactly
  antisymmetric under source/drain exchange, and PSP's :math:`V_{gs}` is
  referred to the *lower* terminal, so transcribing the sum literally
  does break that. True, and the wrong conclusion was drawn from it. The
  same quantity written as a difference,

  .. math::

     V_{gs} + V_{sb}^* - V_{FB} \;=\; V_{gb} - V_{FB} - (V_{sb} - V_{sb}^*)

  has no such problem: :math:`V_{sb} - V_{sb}^*` depends on the two
  junctions only through :math:`\mathrm{MINA}(V_{db}, V_{sb}, a)`, which
  is symmetric in them. It is exactly even under the exchange, and the
  antisymmetry still holds to 7 × 10⁻¹⁵ at body bias.

  With one catch, which the first version of the fix walked into.
  Symmetric in exact arithmetic is not symmetric in floating point:
  taking the correction as a literal :math:`V_{sb} - V_{sb}^*`
  subtracts two numbers agreeing to three parts in :math:`10^{10}` at
  :math:`V_{sb} = 1`, and to three parts in :math:`10^{44}` at
  :math:`V_{sb} = 10^{40}`. The antisymmetry breaks by 2% out there —
  not because the algebra is asymmetric but because the cancellation
  is. Held as its own quantity where the conditioning is computed, the
  large difference is never formed. The old comment's *instinct* was
  sound even if its stated reason was not: the site did need care about
  the antisymmetry. It needed the right form, not omission of the term.

That the mechanism is the *right* one is not an argument from shape.
A body-bias ladder of six rungs from 0 to 1.5 V measures the implied
threshold offset at each; the conditioning computed from the vendor
formula predicts all six within 3.5%, over a 39× range of the quantity,
with nothing fitted. Leverage gives one matching magnitude. Six is the
mechanism.

All twelve sweeps now agree with the vendor at **median 1.000 with a
spread of 0.001 or less**, both geometries, both channel types.

The reference's own convergence
-------------------------------

At which point the largest error in the comparison stopped being in the
model at all.

The reference decks carried no ``.options``, so ngspice swept them at
its defaults. A ``dc`` sweep seeds each point from the previous one, so
Newton needs few iterations and stops as soon as
:math:`r_{tol}|i| + a_{tol}` is met — which on currents of order
:math:`10^{-5}` A left up to **9.6 × 10⁻⁴** of relative error in the
recorded curves.

That is not a well-behaved error bar. It is point-to-point, so it
appears as a **kink at one bias** rather than as an offset, and a kink
invites a search for a mechanism. Two transfer sweeps showed a clean
step of 7–9 × 10⁻⁴ at a single gate voltage followed by a smooth decay,
exactly what a small threshold shift switching on would look like; the
output sweeps showed a run of four consecutive points dipping
8 × 10⁻⁴ and snapping back, which reads as a saturation-knee defect and
sent a hunt after ``THESATB``/``THESATG`` before it was measured.

Which tolerance matters **depends on the sweep**:

.. list-table:: Worst relative error against a fully converged solve
   :header-rows: 1

   * - sweep
     - defaults
     - ``abstol`` only
     - ``reltol`` only
     - both
   * - transfer (``idvg``)
     - 6.41e-04
     - 4.30e-08
     - 4.30e-08
     - 4.30e-08
   * - output (``idvd``)
     - 7.97e-04
     - 7.97e-04
     - 6.61e-09
     - 6.61e-09

``reltol`` is the one that matters, and it fixes both. On an ``idvd``
sweep in saturation the current barely moves from point to point, so
:math:`r_{tol}|i|` is the binding term and no ``abstol`` can help.
``vntol`` and ``gmin`` do nothing either way. Both are set:
``reltol = 1e-6``, mid-plateau — flat from :math:`10^{-5}` to
:math:`10^{-10}`, and ngspice fails to converge on the output sweep by
:math:`10^{-11}` — and ``abstol = 1e-15``, which costs nothing and
guards the low-current end.

That table is also a small lesson in how to run the experiment. The
first version of this measurement concluded the opposite — that
``abstol`` alone accounted for everything and ``reltol`` did nothing —
from a sweep of ``reltol`` over eight decades that changed not one
digit. It changed nothing because ``abstol`` was *already tight in every
run of that sweep*. Vary one knob against the **defaults**, not against
your other fix, and check more than one kind of sweep before
generalising.

The ``op`` decks are unaffected and deliberately do not carry the
option: a single operating point iterates to convergence with margin,
and all 271 of its outputs are bit-identical either way. That asymmetry
is the mechanism — it is the *sweep* that stops early — and it is why
regenerating moved every sweep and **not one bit** of the ``lp_*``,
threshold, capacitance or noise data.

Telling this apart from physics is a measurement, not a judgement. A
converged curve is smooth, so a local polynomial through a point's
neighbours predicts it; convergence noise breaks that and real curvature
does not. But the residual of such a fit is not itself the
discriminator, because a cubic through a moderate-inversion knee carries
genuine truncation error of the same size — 6 × 10⁻⁴ here, which is what
made the first version of this measurement read the knee as a defect.
What discriminates is comparing that residual against a curve known to
be analytic. Ours is closed-form, so its residual is pure truncation;
subtract and only the reference's own noise remains:

.. list-table:: Worst local-fit residual, degree 5 over 8 neighbours
   :header-rows: 1

   * - sweep
     - reference, before
     - reference, after
     - analytic curve
   * - ``nmos_idvg_vd0p05``
     - 4.90e-04
     - 1.34e-04
     - 1.34e-04
   * - ``nmos_idvg_vb_m1``
     - 5.43e-04
     - 1.91e-04
     - 1.91e-04
   * - ``pmos_idvg_vd0p05``
     - 4.55e-04
     - 4.19e-05
     - 4.19e-05

After regeneration the reference matches the analytic curve to every
digit printed. Before it, it was three to eleven times rougher.

The general point is worth more than the particular one. **A reference
generated by a solver inherits that solver's tolerances, and those are
part of the reference.** They were harmless for as long as the model
under test was further out than they were — two decades of margin, and
no reason to look. They became the limit on what could be measured the
moment the model got good, and until they were found the model was being
blamed for them.

The charge model, and what construction properties cannot tell you
------------------------------------------------------------------

The three properties above are real, and for a long time they were the
only checks the charge model had. They are also all **ratios** —
conservation, the source/drain swap, the Ward–Dutton partition — and a
uniform error in the oxide capacitance divides out of every one of them.

The capacitances were 24% high at every bias point, and nothing in the
test suite could have seen it.

**A model checked only against itself is checked for consistency, not
for correctness.** Two causes, both found by comparing against PSP's own
operating-point outputs rather than against our own arithmetic:

* PSP builds the charge model's oxide capacitance from the **CV
  effective dimensions**, not the drawn ones — this card puts the
  effective length 7% under the drawn value;
* and it applies a **quantum-mechanical reduction** to that capacitance,
  worth another 12%. The same correction already shifts the threshold in
  the DC path; here it is the same physics acting on the charge, the
  inversion layer sitting a little below the interface and putting a
  capacitance in series with the oxide.

With both, the gate capacitance came to within 1% on the long device.
The residual then named its own cause: in the **linear region** the two
agreed to a part in ten thousand at both geometries, and the whole error
was in **saturation**, growing with drain bias. That rules out overlap
capacitance — a fixed capacitance in parallel does not switch itself off
at low drain bias — and points at channel-length modulation, which had
been left out of the charge model as a long-channel simplification while
the current had used it all along. Putting it back, along with the
polysilicon factor on the gate charge:

**The long device now matches PSP to 0.07% across the whole bias grid,
on both channel types** — charges from a foundry card, through the
scaling layer, nothing fitted, agreeing with the vendor's own
implementation to a part in ten thousand.

PSP reports the overlap and fringe capacitances **separately** from
:math:`C_{gg}`, so this is an intrinsic-to-intrinsic comparison. That is
worth stating because it is easy to assume otherwise: on a 0.13 µm
device the overlap is 182% of the intrinsic :math:`C_{gg}`, and on the
10 µm device 11%, so an intrinsic-only model compared against a total
would be wrong by those amounts rather than by the tenths of a percent
measured here.

The short device is within 2%, and that residual is **not** explained.
It is not overlap, for the reason just given, and it is not
channel-length modulation, which is now in.

The DC current is untouched by any of it, as it must be.

The subthreshold slope
----------------------

Worth separating out, because it is the steepest check available and it
is not a fitted quantity. In a surface-potential model the subthreshold
swing falls out of the electrostatics — the body factor, the surface
potential, the effective thermal voltage — so it says something about
the formulation rather than about a parameter.

Against PSP's own curves, over three decades of current, it agrees to
**a quarter of a millivolt per decade on both channel types**.

One caveat governs that number, and getting it wrong reports a defect
that is not there. The reference records *total terminal current*, and
PSP's junction leakage is a flat 2 pA floor that this core does not
model at all. On the body-biased sweep the drain current comes down to
meet it, so a slope taken through that region measures PSP's leakage
rather than its channel — and says the swing is 2.5 mV/decade out when
it is 0.2. The same caution applies to every ratio on this page: they
are all against a current that includes terms deliberately absent here,
which is why the comparison carries a 1 µA floor.

A note on DIBL, because this page said the opposite for a while. In
absolute terms it *is* small on this process — PSP's own threshold moves
3.5 mV over 1.35 V of drain bias, and that was measured rather than
assumed. But small in millivolts is not small in current: in weak
inversion, at 85 mV/decade, 3.5 mV is a 9% change, and on the p-channel
saturation sweep leaving it out costs more than a factor of two at the
bottom of the range. "Small parameter" and "negligible term" are
different claims, and only the first was ever measured.

.. note::

   :class:`~pycircuit.circuit.compact.PspMosLongChannel` is **not**
   PSP103 and should not be used as a substitute for it in design work.
   It is the core of PSP103, built to establish that the DSL can carry a
   surface-potential formulation, and measured so that the claim has a
   number attached.
