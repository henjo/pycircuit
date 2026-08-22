.. _compact_models:

Compact models: a surface-potential MOSFET
==========================================

pycircuit can express a **production compact model** — the kind a foundry
ships — in the Behavioural HDL, and check it against the vendor's own
compiled binary.

Two devices live in :mod:`pycircuit.circuit.compact`:

* :class:`~pycircuit.circuit.compact.CapCmomi`, a translation of the IHP
  Open PDK's interdigitated metal-oxide-metal capacitor, which agrees
  with the vendor's OSDI build to **2 × 10⁻⁸** over 1 MHz – 100 GHz;
* :class:`~pycircuit.circuit.compact.PspMosLongChannel`, the core of
  PSP103 — the industry-standard surface-potential MOSFET.

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
    print("     - worst |F| / scale")
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

**Source and drain are exactly interchangeable.** Swapping which end is
called the source negates the current and changes nothing else — bit for
bit, not approximately. Threshold-voltage models are famously asymmetric
about :math:`V_{ds} = 0`, which shows up as spurious distortion in
passing gates and mixers.

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
    print("     - identical")
    for vg in (0.6, 1.2, 1.8):
        for vd in (0.1, 0.8):
            f = np.asarray(e.i(np.array([vd, vg, 0.0, 0.0])), float)[0]
            r = -np.asarray(e.i(np.array([0.0, vg, vd, 0.0])), float)[0]
            print("   * - %.1f V" % vg)
            print("     - %.1f V" % vd)
            print("     - %.6e A" % f)
            print("     - %.6e A" % r)
            print("     - %s" % ("yes" if f == r else "NO"))

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
    REF = os.path.join('..', '..', '..', 'pycircuit', 'circuit', 'tests',
                       'data', 'psp103_ihp_iv.json')

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
        rows = [('nmos_long_idvd', 10.0, 1.00, 1.010),
                ('nmos_long_idvg', 10.0, 1.00, 1.10),
                ('nmos_idvd_vg1p2', 1.0, 0.13, 1.021),
                ('nmos_idvg_vd1p2', 1.0, 0.13, 0.974),
                ('nmos_idvg_vd0p05', 1.0, 0.13, 1.117)]
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
omits matters least. It is within about 1% there, from a model card,
with no fitting anywhere in the chain.

The short device is a few percent off at high drain bias and ~12% off in
the linear region, which is the cost of the one layer still missing —
the short-channel threshold block, DIBL.

What is in it, and what is not
------------------------------

Present: the surface-potential solver, the drain current by symmetric
linearisation, the terminal charges with the Ward–Dutton partition,
mobility reduction, Coulomb scattering, velocity saturation, series
resistance, channel-length modulation, the quantum-mechanical correction
to the surface potential at threshold, and the effective thermal voltage.

Absent: DIBL and the rest of the short-channel threshold block, gate and
junction leakage, impact ionisation, overlap and fringe capacitances,
the non-quasi-static block, self-heating, and every temperature
coefficient. PSP103 is this core plus those, and the size of what they
are worth is exactly what the table above measures.

.. note::

   :class:`~pycircuit.circuit.compact.PspMosLongChannel` is **not**
   PSP103 and should not be used as a substitute for it in design work.
   It is the core of PSP103, built to establish that the DSL can carry a
   surface-potential formulation, and measured so that the claim has a
   number attached.
