.. _idtmod:

Circular integrators
====================

This page documents pycircuit's phase-accumulator elements: the plain
integrator :class:`~pycircuit.circuit.elements.Idt`, the circular (modulo)
integrator :class:`~pycircuit.circuit.elements.Idtmod` -- pycircuit's
implementation of the Verilog-A ``idtmod`` operator -- and the two
phasor-pair variants :class:`~pycircuit.circuit.elements.IdtmodCircular`
and :class:`~pycircuit.circuit.elements.IdtmodQuadrature`.

Every number on this page is **recomputed when the documentation is
built**: the tables and measurements below are produced by running the
elements, not pasted in.  The full research record behind the design --
simulator survey, literature, and the measured design decisions -- is in
``idtmod.md`` at the repository root.

The elements at a glance
------------------------

.. list-table::
   :header-rows: 1
   :widths: 18 42 40

   * - element
     - what it computes
     - use it when
   * - ``Idt``
     - :math:`v_{out} = \int (v_{ip} - v_{in})\,dt`
     - you need the plain, unbounded integral
   * - ``Idtmod``
     - the integral folded into ``[offset, offset+modulus)``
     - phase accumulators: VCO/PLL phase, anything measured in
       accumulated phase.  **The default choice.**
   * - ``IdtmodCircular``
     - the same folded value, computed from a smooth
       :math:`(\cos, \sin)` state pair
     - a research/PSS stepping stone; at the terminals ``Idtmod`` is
       strictly more accurate
   * - ``IdtmodQuadrature``
     - :math:`A\cos` and :math:`A\sin` of the phase on two output pairs
     - sinusoidal VCO / I/Q / LO macromodels; periodic-steady-state
       work -- no sawtooth anywhere in the circuit

A first example -- the wrapped ramp of :doc:`example 6
<examples/example6>` with Kundert's canonical arguments (modulus 1,
offset -0.5, so the phase lives in the symmetric range
:math:`[-\tfrac12, \tfrac12)`):

.. plot::
    :include-source: True
    :width: 10cm

    import pylab
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit import SubCircuit, VS, R, Idtmod, gnd

    c = SubCircuit()
    nin, nout = c.add_nodes('in', 'out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0, offset=-0.5)

    tran = Transient(c, toolkit=numeric, uic=True)
    result = tran.solve(tend=3.0, timestep=0.01)

    t = result.v(nout).x[0]
    pylab.plot(t, result.v(nout).y)
    pylab.xlabel('t [s]'); pylab.ylabel('v(out) [V]')
    pylab.title('idtmod(1V, 0, 1, -0.5): phase folded into [-0.5, 0.5)')
    pylab.grid(True)

The contract
------------

``Idtmod`` implements the Verilog-AMS ``idtmod(expr, ic, modulus,
offset)`` operator.  Writing :math:`\phi(t) = \int_0^t expr\,d\tau + ic`
for the running integral, the output is the unique :math:`k` with

.. math::

   \phi(t) = n \cdot modulus + k, \qquad
   offset \le k < offset + modulus, \qquad n \in \mathbb{Z}

realised by the *floored* wrap
:math:`k = y - m\lfloor (y-o)/m \rfloor`.  The distinction from a naive
remainder matters whenever ``offset`` is not a multiple of the modulus --
the congruence with the integral must hold, not just the range:

.. exec-rst::

    from pycircuit.circuit.elements import floored_wrap
    from pycircuit.circuit.toolkit import numeric

    y, m, o = 1.7, 1.0, -0.5
    k = floored_wrap(y, m, o, numeric)
    naive = (y % m) + o
    print("For an integral of %.1f with modulus %.0f and offset %.1f:" % (y, m, o))
    print("")
    print("* ``floored_wrap`` returns **%.1f** -- and %.1f - (%.1f) = %.0f" % (k, y, k, y - k))
    print("  is a whole number of moduli, as the contract requires;")
    print("* the naive ``(y %% m) + offset`` returns %.1f, which differs from" % naive)
    print("  the integral by %.1f -- NOT a whole number of moduli.  (This was" % (y - naive))
    print("  a live bug in the pre-rework element.)")

Two boundary rules complete the contract, both taken from the LRM:

* ``modulus <= 0`` (or unset) degrades ``Idtmod`` to a plain ``Idt`` --
  the operator then "does not limit the output".
* The result may land one float ulp outside the half-open range when the
  integral sits exactly on a boundary; consumers of a wrapped value are
  periodic in it, so this is harmless and deliberately not clamped.

**Why the operator exists at all** is floating-point precision, not
convenience.  A VCO phase integrated without folding grows without bound,
and a double has only so many digits:

.. exec-rst::

    import numpy as np
    for cycles in (1e3, 1e6, 1e9):
        phase = 2 * np.pi * cycles
        print("* after %.0e cycles the unwrapped phase is ~%.1e rad and its ulp" % (cycles, phase))
        print("  is %.1e rad -- %.1e of a cycle is simply gone." % (np.spacing(phase), np.spacing(phase) / (2 * np.pi)))

A bounded state keeps the ulp at the size of one modulus forever.  How
pycircuit achieves that *without* upsetting the integrator is the subject
of `Exactness: the gauge shift`_ below.

Parameters
----------

.. parametertable:: pycircuit.circuit.elements.Idtmod.instparams

``ic`` deserves its own words, because it controls the DC behaviour:

* **With** ``ic`` **given**, the LRM applies -- "the initial condition
  shall force the DC solution" -- and the operating point solves with the
  output pinned to the wrapped ``ic``.  ``uic=True`` is then optional.
* **Without** ``ic`` (``None``), the LRM's other branch applies: the
  integrator has a DC solution only if feedback forces its input to zero.
  A driven integrator's bias solve is then singular *by specification*,
  and the transient must start from explicit initial conditions
  (``Transient(..., uic=True)``).

Both behaviours, measured at build time:

.. exec-rst::

    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, Idtmod, gnd
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.analysis import NoConvergenceError, SingularMatrix

    def build(**kw):
        c = SubCircuit()
        nin, nout = c.add_nodes('in', 'out')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['R'] = R(nout, gnd, r=1e3)
        c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, **kw)
        return c

    res = DC(build(modulus=1.0, offset=-0.5, ic=1.7), toolkit=numeric).solve()
    print("* ``ic=1.7`` with modulus 1, offset -0.5: the DC output is pinned to")
    print("  wrap(1.7) = **%.1f V** (measured: %.12g V)." % (-0.3, float(res.v('out'))))
    try:
        DC(build(modulus=1.0), toolkit=numeric).solve()
        print("* without ic: unexpectedly converged -- THIS IS A REGRESSION")
    except (NoConvergenceError, SingularMatrix):
        print("* without ``ic``: the bias solve fails loudly (as the LRM's no-ic")
        print("  branch requires for a driven integrator) -- use ``uic=True``.")

Under ``uic=True`` an element ``ic`` seeds the internal state directly
(wrapped into range), so the run starts at the requested phase without a
bias solve.

Exactness: the gauge shift
--------------------------

A folded phase presents a genuine numerical dilemma.  The wrap is a jump
of exactly one modulus, and multistep integrators (trapezoidal, Gear)
assume the state history lies on a smooth polynomial -- feed them a
jump and the step control collapses.  But *not* folding the state loses
precision, as computed above.  General-purpose ODE tools resolve this
with an event: stop the integrator at the wrap, apply the jump, restart
from scratch -- once per cycle, which is prohibitive for a VCO.

pycircuit uses the structure of the problem instead.  The jump is always
by an **exact multiple of the modulus**, and linear-multistep formulas
are *translation invariant*: they consume differences of history values,
so subtracting the same constant from the current state **and from every
stored history entry** changes nothing they can see.  After each accepted
step the solver therefore re-folds the state by :math:`n \cdot modulus`
across the whole live history window -- an exact operation, no event, no
restart, no order loss.  (Elements declare such states through
``periodic_states()``; the machinery is generic.)

The measurable consequence: for a constant input the folded phase is
**exact to the ulp, indefinitely** -- there is no method error left, only
float noise.  Measured at build time with the per-step increment made
deliberately inexact (input 0.2 V, dt 0.5, so each step adds exactly
1/10 modulus and the analytic reference is the exact cyclic sequence
``k/10``):

.. exec-rst::

    import warnings
    import numpy as np
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, Idtmod, gnd
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import EulerIntegrator

    c = SubCircuit()
    nin, nout = c.add_nodes('in', 'out')
    c['vin'] = VS(nin, gnd, v=0.2)
    c['R'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)
    tran = Transient(c, toolkit=numeric, uic=True, integrator=EulerIntegrator())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=500.0, timestep=0.5, fixed_timestep=True)
    y = res.v('out').y
    k = np.arange(len(y))
    d = np.abs(y[1:] - ((k[1:] % 10) / 10.0))
    d = np.minimum(d, 1.0 - d)
    print("After **%d wraps** (%d steps), the largest phase error against the" % (100, len(y) - 1))
    print("exact analytic value is **%.1e** -- about one ulp, and flat in run" % d.max())
    print("length.  Without the bounded state the same run drifts at the")
    print("eps-times-integral level (measured 1.4e-12 here, growing linearly;")
    print("see ``test_long_run_precision_payoff``).")

Wrap corners and breakpoints
----------------------------

The folded *output* is a sawtooth -- that discontinuity is what the
operator means, and no implementation removes it.  What the solver can do
is refuse to be surprised by it.  ``Idtmod`` predicts the time of its
next wrap from the last accepted point (the phase rate is just the input
voltage, so the prediction is a one-liner) and publishes it through the
element-event mechanism; the adaptive controller then **lands a step
boundary on the corner** and drops the integration order across it, so no
polynomial is ever fitted through the jump.  Measured at build time:

.. exec-rst::

    import warnings
    import numpy as np
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, Idtmod, gnd
    from pycircuit.circuit.transient import Transient

    c = SubCircuit()
    nin, nout = c.add_nodes('in', 'out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0, offset=-0.5)
    tran = Transient(c, toolkit=numeric, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=2.2, timestep=1e-2)
    t = res.v('out').x[0]
    s = res.statistics
    print("Wrap crossings at t = 0.5 and 1.5 s; the run landed timepoints within")
    print("**%.1e s** and **%.1e s** of them (%d breakpoints hit, %d steps, %d rejections)." % (
        np.min(np.abs(t - 0.5)), np.min(np.abs(t - 1.5)),
        s.breakpoints_hit, s.accepted_steps, s.rejected_steps))

One consequence worth knowing when you compare waveforms: a sample
landing *exactly on* a wrap is a knife edge -- the sawtooth is
double-valued there, and which limit the solver reports is decided by
sub-ulp rounding.  Compare wrapped waveforms **congruence-style**
(distance modulo the modulus), not sample-by-sample.

The phasor variants
-------------------

Theory: a circle instead of a sawtooth
``````````````````````````````````````

The classical alternative to folding a scalar phase is to integrate its
unit phasor.  With :math:`\omega = 2\pi(v_{ip}-v_{in})/modulus`,

.. math::

   \dot c = -\omega s, \qquad \dot s = +\omega c

rotates the point :math:`(c, s)` around the unit circle; the folded phase
is recovered as :math:`k = \mathrm{wrap}(m\,\mathrm{atan2}(s,c)/2\pi)`.
The **state** is now smooth for all time -- nothing ever jumps.  The
price is twofold.  First, accuracy: the state carries the oscillation's
curvature, so every integrator leaves its characteristic per-cycle phase
lag (quantified below).  Second, the **radius drifts**: numerical
integration does not conserve :math:`r^2 = c^2 + s^2`.  Gear-type
methods damp it (the classic "Gear damping" that shrinks LC-tank
oscillations in SPICE), forward rules grow it, and even trapezoidal --
exactly radius-preserving for constant :math:`\omega`, because its update
is then a Cayley transform of a skew-symmetric matrix, i.e. a pure
rotation -- drifts in real runs through order-dropped steps and varying
input.

Baumgarte orbit correction
``````````````````````````

The cure implemented here is a *Baumgarte-style* constraint feedback,
the same device used to stop quaternion attitude states drifting off
unit norm.  The dynamics gain a term pointing back to the circle:

.. math::

   \dot c = -\omega s - \gamma|\omega|\,c\,(r^2-1), \qquad
   \dot s = +\omega c - \gamma|\omega|\,s\,(r^2-1)
   \quad\Rightarrow\quad
   \frac{d(r^2)}{dt} = -2\gamma|\omega|\,r^2(r^2-1)

so the unit circle is exponentially attracting.  Two properties matter
to a user:

* **The correction cannot bias the phase.**  Its contribution to the
  angular rate is :math:`(c\dot s - s\dot c)_{corr} =
  c(\mathrm{corr}\cdot s) - s(\mathrm{corr}\cdot c) = 0` -- it is purely
  radial.
* **The gain needs no tuning.**  The classic complaint about Baumgarte
  stabilization is that its gain depends on the integrator and the step
  size.  Here :math:`\gamma` is *per radian of phase travel*: the
  effective per-step strength is :math:`\gamma\,|\omega| h`, and since
  the step controller keeps phase-per-step roughly constant, the default
  :math:`\gamma = 1` stays inside the theoretically admissible range at
  any frequency.  (Push :math:`\gamma|\omega|h` past 1 and trapezoidal
  will ring the *radius* at the step rate -- a sentinel test guards that
  boundary.)

The drift and the cure, measured at build time (Gear-2, 50 cycles):

.. plot::
    :include-source: True
    :width: 10cm

    import warnings
    import numpy as np
    import pylab
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, gnd
    from pycircuit.circuit.elements import IdtmodCircular
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import Gear2Integrator

    def radius(gamma):
        c = SubCircuit()
        nin, nout = c.add_nodes('in', 'out')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['R'] = R(nout, gnd, r=1e3)
        c['X'] = IdtmodCircular(nin, gnd, nout, gnd, modulus=1.0,
                                gamma=gamma, ic=0.0)
        rows = [i for i, nm in enumerate(map(str, c.nodes))
                if 'cos_node' in nm or 'sin_node' in nm]
        tran = Transient(c, toolkit=numeric, uic=True,
                         integrator=Gear2Integrator())
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=50.0, timestep=2e-2, fixed_timestep=True)
        X = np.asarray(res.x, float)
        return np.asarray(res.v('out').x[0], float), np.hypot(X[rows[0]], X[rows[1]])

    for gamma, label in ((0.0, 'gamma = 0 (no correction)'),
                         (1.0, 'gamma = 1 (default)')):
        t, r = radius(gamma)
        pylab.plot(t, r, label=label)
    pylab.xlabel('t [s]'); pylab.ylabel('phasor radius r')
    pylab.title('Gear-2 damps the phasor; Baumgarte correction holds the circle')
    pylab.legend(); pylab.grid(True)

Accuracy: the two error regimes
```````````````````````````````

The scalar and phasor forms are not interchangeable in accuracy.  The
scalar ``Idtmod`` integrates a *linear* phase exactly (any LMS formula
does), so its only residue is float noise.  The phasor pays each method's
per-cycle phase lag at the carrier rate -- :math:`(\omega h)^2/6`,
:math:`/12`, :math:`/3` per cycle for Euler, trapezoidal and Gear-2 --
accumulating linearly with cycle count.  Regenerated at build time
(constant 1 V input, dt = 0.01), with the closed-form prediction beside
each phasor measurement:

.. exec-rst::

    import warnings
    import numpy as np
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, gnd
    from pycircuit.circuit.elements import Idtmod, IdtmodCircular
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import (EulerIntegrator,
                                              TrapezoidalIntegrator,
                                              Gear2Integrator)

    DT = 0.01
    WH = 2 * np.pi * DT
    COEFF = {'EulerIntegrator': 1 / 6.0, 'TrapezoidalIntegrator': 1 / 12.0,
             'Gear2Integrator': 1 / 3.0}

    def err(cls, integ, tend):
        c = SubCircuit()
        nin, nout = c.add_nodes('in', 'out')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['R'] = R(nout, gnd, r=1e3)
        c['X'] = cls(nin, gnd, nout, gnd, modulus=1.0, ic=0.0)
        tran = Transient(c, toolkit=numeric, uic=True, integrator=integ())
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=tend, timestep=DT, fixed_timestep=True)
        y = res.v('out').y
        t = res.v('out').x[0]
        d = np.abs(y[1:] - t[1:] % 1.0)
        d = np.minimum(d, 1.0 - d)
        return d[-max(4, len(d) // 20):].max()

    print(".. list-table:: Phase error at end of run (measured at build time)")
    print("   :header-rows: 1")
    print("")
    print("   * - integrator")
    print("     - cycles")
    print("     - Idtmod")
    print("     - IdtmodCircular")
    print("     - predicted lag")
    for integ in (EulerIntegrator, TrapezoidalIntegrator, Gear2Integrator):
        for tend in (2.0, 20.0):
            pred = COEFF[integ.__name__] * WH ** 2 * tend
            print("   * - %s" % integ.__name__.replace('Integrator', ''))
            print("     - %g" % tend)
            print("     - %.1e" % err(Idtmod, integ, tend))
            print("     - %.1e" % err(IdtmodCircular, integ, tend))
            print("     - %.1e" % pred)
    print("")
    print("At 200 cycles the gap reaches ~9 orders of magnitude; the full table,")
    print("with a loud consistency check against the closed form, is")
    print("``benchmarks/idtmod_phase_error.py``.")

The practical rule follows directly: **when accumulated phase is the
measurement** -- jitter, long PLL transients, cycle counting -- **use**
``Idtmod``.

IdtmodQuadrature: quadrature outputs, no sawtooth
`````````````````````````````````````````````````

When what you actually need is the *sinusoid* (a sinusoidal VCO, an I/Q
source, a mixer LO), fold nothing: ``IdtmodQuadrature`` has six
terminals and drives its two output pairs with
:math:`A\cos(2\pi\phi/m)` and :math:`A\sin(2\pi\phi/m)` directly.  There
is no wrapped output, no ``atan2``, and therefore **no discontinuity
anywhere in the circuit** -- no wrap events fire at all.  It is also the
representation suited to periodic-steady-state / shooting analysis when
one is implemented: over one period this state vector returns exactly to
itself, which a scalar phase (advancing by one modulus per period)
structurally cannot do.

.. parametertable:: pycircuit.circuit.elements.IdtmodQuadrature.instparams

.. exec-rst::

    import warnings
    import numpy as np
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import SubCircuit, VS, R, gnd
    from pycircuit.circuit.elements import IdtmodQuadrature
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.integrator import TrapezoidalIntegrator

    def build(**kw):
        c = SubCircuit()
        nin = c.add_node('in')
        nc, ns = c.add_nodes('outc', 'outs')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['Rc'] = R(nc, gnd, r=1e3)
        c['Rs'] = R(ns, gnd, r=1e3)
        c['X'] = IdtmodQuadrature(nin, gnd, nc, gnd, ns, gnd, **kw)
        return c

    res = DC(build(ic=0.25, modulus=1.0, amplitude=2.0), toolkit=numeric).solve()
    print("* DC with ``ic=0.25`` (a quarter turn) and ``amplitude=2``: the outputs")
    print("  pin to (cos, sin) = (%.3g, %.3g) V -- the phasor starts ON the circle." % (
        float(res.v('outc')), float(res.v('outs'))))

    c = build(ic=0.0, modulus=1.0)
    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=TrapezoidalIntegrator())
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        r = tran.solve(tend=2.0, timestep=1e-2, fixed_timestep=True)
    t = r.v('outc').x[0]
    ec = np.abs(r.v('outc').y[1:] - np.cos(2 * np.pi * t[1:])).max()
    print("* two cycles of oscillation: max deviation from the analytic cosine")
    print("  %.1e V, and **%d wrap breakpoints fired** -- there is no sawtooth" % (
        ec, r.statistics.breakpoints_hit))
    print("  for the solver to negotiate (Idtmod lands one per cycle here).")

Choosing an element
-------------------

.. list-table::
   :header-rows: 1
   :widths: 45 25

   * - you want
     - use
   * - a phase accumulator; accumulated phase is the measurement
       (jitter, PLL transients, cycle counts)
     - ``Idtmod``
   * - a plain unbounded integral
     - ``Idt`` (or ``Idtmod`` with ``modulus <= 0``)
   * - sinusoidal/quadrature outputs; no discontinuities anywhere;
       future periodic-steady-state work
     - ``IdtmodQuadrature``
   * - the folded value from a smooth state (research, PSS stepping
       stone)
     - ``IdtmodCircular``

Limitations, briefly: the phasor elements need ``modulus > 0`` and always
pin their DC point from ``ic`` (their state must start on the circle);
under ``JAXTransient`` start them with ``uic=True``.  None of the four
support the symbolic toolkit's transient (``Idt``/``Idtmod`` AC analysis
works and gives the expected :math:`1/s`).

Further reading
---------------

* ``idtmod.md`` (repository root) -- the full research record: LRM
  semantics, how ngspice/Xyce/VACASK/gnucap handle (or don't handle) the
  operator, the gauge-shift theory, Baumgarte gain selection, and every
  measured design decision with references.
* ``benchmarks/idtmod_phase_error.py`` -- the complete accuracy
  comparison, self-checking against the closed-form lag predictions.
* ``pycircuit/circuit/tests/test_idtmod.py`` -- the executable contract:
  26 tests covering everything claimed on this page.
