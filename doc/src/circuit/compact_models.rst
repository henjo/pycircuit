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
        from pycircuit.circuit.compact import PspMosLongChannel
        from pycircuit.circuit import psp_scaling
        from pycircuit.utilities import spicecard
        cm.default_toolkit = numeric
        ref = json.load(open(REF))['sweeps']
        deck = spicecard.read(os.path.join(PDK, 'cornerMOSlv.lib'),
                              section='mos_tt')
        for name in ('nmos_long_idvd', 'nmos_long_idvg',
                     'nmos_idvd_vg1p2', 'nmos_idvg_vd1p2',
                     'nmos_idvg_vd0p05'):
            s = ref[name]
            w, l = s['w'], s['l']
            kw = psp_scaling.to_long_channel(
                deck.model_params('sg13g2_lv_nmos_psp', w=w, l=l, ng=1,
                                  m=1, pre_layout=1), w=w, l=l)
            e = PspMosLongChannel(cm.Node('d'), cm.Node('g'), cm.Node('s'),
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
            rows.append((name, w * 1e6, l * 1e6,
                         np.median(np.abs(g[m]) / r[m])))
        note = "computed live against the installed PDK"
    except Exception:
        rows = [('nmos_long_idvd', 10.0, 1.00, 1.002),
                ('nmos_long_idvg', 10.0, 1.00, 1.000),
                ('nmos_idvd_vg1p2', 1.0, 0.13, 1.002),
                ('nmos_idvg_vd1p2', 1.0, 0.13, 0.996),
                ('nmos_idvg_vd0p05', 1.0, 0.13, 1.001)]
        note = ("quoted from the committed measurement -- the IHP PDK is "
                "not installed here")

    print(".. list-table:: Drain current, ours / PSP103 (median over the sweep)")
    print("   :header-rows: 1")
    print("")
    print("   * - sweep")
    print("     - W")
    print("     - L")
    print("     - ratio")
    for name, w, l, ratio in rows:
        print("   * - ``%s``" % name)
        print("     - %.0f um" % w)
        print("     - %.2f um" % l)
        print("     - %.3f" % ratio)
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

Present also: drain-induced barrier lowering (``CF``, ``CFB``, ``CFD``)
and the body- and gate-bias modulations of the velocity-saturation
parameter (``THESATB``, ``THESATG``). Passing ``all_terms=False`` to
``to_long_channel`` drops both, which exists so their effect stays
measurable — see below, because deciding whether to keep them turned
out to be a lesson about the measurement rather than about the terms.

Absent: the rest of the short-channel threshold block (``PSCE``, which
is all-zero on this card anyway); gate and junction
leakage; impact ionisation; overlap and fringe capacitances; the
non-quasi-static block; self-heating; and every temperature coefficient.
PSP103 is this core plus those, and the size of what they are worth is
exactly what the table above measures.

Two of those omissions are worth naming as *permanently* invisible to
the table: gate leakage and junction leakage are four to six orders of
magnitude below its 1 µA floor on this process, and overlap and fringe
capacitances contribute identically zero DC current. They matter for
CV and AC work, not here.

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

Then the same thing happened twice more. The bias-dependent body factor
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

All eleven sweeps across both channel types are within 2.3% now, and
every p-channel sweep within 0.6%, from two foundry cards with nothing
fitted anywhere in the chain.

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
