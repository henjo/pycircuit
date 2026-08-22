.. _hdl:

Behavioural elements: the HDL
=============================

pycircuit lets you write a circuit element as **equations** rather than
as matrix stamps. The language is shaped like Verilog-A's analog block —
you write contribution statements — and the
:class:`~pycircuit.circuit.hdl.Behavioural` metaclass compiles them,
symbolically, into everything an element must supply: the current and
charge vectors, exact Jacobians, sources, noise, initial conditions, and
the hooks the JAX and PCNR machinery need.

Every number and figure on this page is **recomputed when the
documentation is built**. The research record, the full capability map
against Verilog-A, and the development plan live in ``hdl.md`` at the
repository root.

A first element
---------------

The argument names of ``analog()`` declare the terminals; the body
returns contribution statements:

.. exec-rst::

    import inspect
    from pycircuit.circuit import elements_hdl
    src = inspect.getsource(elements_hdl.RHdl)
    print(".. code-block:: python")
    print("")
    for line in src.splitlines():
        print("   " + line)

That is the entire element. Used like any other:

.. exec-rst::

    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.dcanalysis import DC

    cm.default_toolkit = numeric
    c = SubCircuit()
    na, nb = c.add_nodes('a', 'b')
    c['vs'] = VS(na, gnd, v=3.0)
    c['R1'] = eh.RHdl(na, nb, r=2e3)
    c['R2'] = eh.RHdl(nb, gnd, r=1e3)
    res = DC(c, toolkit=numeric).solve()
    print("A 2k/1k divider from two HDL resistors gives "
          "``v(b) = %.6f V`` (3 x 1/3 = 1.0)." % float(res.v('b')))

**Why write elements this way?** A hand-written element states its
physics twice — once as ``i(x)``, once as the ``G(x)`` you differentiated
yourself — and nothing checks that the two agree. That is not a
hypothetical failure: ``test_element_jacobians.py`` exists because
``VCVS_limited`` shipped a ``G`` that had dropped its gain entirely. In
the HDL the equation is written once and the Jacobian is *derived*, so
they cannot disagree:

.. exec-rst::

    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import elements_hdl as eh

    cm.default_toolkit = numeric
    el = eh.DiodeHdl('p', 'n', IS=1e-13)
    el.update_iparv()
    x = np.array([0.42, 0.1])
    eps = 1e-7
    J = np.zeros((2, 2))
    for k in range(2):
        d = np.zeros(2); d[k] = eps
        J[:, k] = (np.asarray(el.i(x + d), float)
                   - np.asarray(el.i(x - d), float)) / (2 * eps)
    err = np.max(np.abs(np.asarray(el.G(x), float) - J))
    print("For the generated diode at v = 0.32 V, the symbolic ``G`` and a")
    print("central difference of ``i`` agree to **%.1e** (relative to a" % err)
    print("conductance of %.3g S)." % float(np.asarray(el.G(x), float)[0, 0]))

The language
------------

Statements and probes
`````````````````````

.. list-table::
   :header-rows: 1
   :widths: 34 30 36

   * - DSL
     - Verilog-A
     - notes
   * - ``Contribution(b.I, expr)``
     - ``I(a,b) <+ expr``
     - accumulates; the usual form
   * - ``Contribution(b.V, expr)``
     - ``V(a,b) <+ expr``
     - creates a branch-current unknown
   * - ``b.V``, ``node.V``
     - ``V(a,b)``, ``V(a)``
     - potential probe
   * - ``b.I``
     - ``I(b)``
     - flow probe (V-contributed branches)
   * - ``Node('mid')``
     - internal node
     - discovered from use

``analog()`` returns one ``Contribution`` or a tuple of them.

Operators
`````````

``ddt``, ``idt``, ``idtmod``, ``ddx``, ``limexp``, ``white_noise``,
``flicker_noise``, ``vt()``/``TEMP`` (``$temperature``/``$vt``), ``TIME``
(``$abstime``), and any sympy function — ``exp``, ``log``, ``tanh``,
``sqrt``, ``Piecewise`` for conditionals.

Two of these carry decisions worth knowing about:

**ddt declares a charge; it never differentiates.** The term goes into
the charge vector and *the simulator's integrator* differentiates it, at
whatever order and step size it is using. A state-dependent scaling is
refused rather than silently mis-integrated, because ``g(v)·ddt(v)`` is
not the time derivative of any charge:

.. exec-rst::

    import sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ddt)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric
    try:
        class Bad(Behavioural):
            instparams = []

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, b.V * ddt(b.V))
        print("no error -- THIS IS A REGRESSION")
    except NotImplementedError as exc:
        print("Writing ``b.V * ddt(b.V)`` is refused at class creation::")
        print("")
        print("   " + str(exc).split('. ')[0] + ".")

**idt and idtmod become real state unknowns**, with their own equation
added to the system — so they work at DC (where an ``ic`` pins the
solution), under ``uic``, and with the bounded-state gauge shift
described in :ref:`idtmod`. A whole circular integrator is one line:

.. exec-rst::

    import inspect
    from pycircuit.circuit import elements_hdl
    src = inspect.getsource(elements_hdl.IdtmodHdl.analog)
    print(".. code-block:: python")
    print("")
    for line in src.splitlines():
        print("   " + line.rstrip())

Parameters
``````````

Declare them with the ordinary
:class:`~pycircuit.utilities.param.Parameter`; inside ``analog()`` the
names are bound as symbols, and values are read from the *resolved*
parameters at call time, so hierarchical expressions (``r='2*rbase'``)
reach the generated code.

What you get for free
---------------------

Exact Jacobians, and JAX
````````````````````````

``G = ∂i/∂x`` and ``C = ∂q/∂x`` are sympy Jacobians of the whole element,
compiled with common-subexpression elimination. Generated elements also
carry the pure forms the JAX backend needs, so they can be evaluated
batched and swept per lane.

Noise
`````

A ``white_noise`` or ``flicker_noise`` contribution populates the noise
correlation matrix — something a SPICE B-source cannot express at all:

.. exec-rst::

    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.elements import R
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       white_noise, KBOLTZMANN, TEMP)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class RNoisy(Behavioural):
        instparams = [Parameter(name='r', desc='R', unit='ohm', default=1e3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, b.V / r),
                    Contribution(b.I, white_noise(4 * KBOLTZMANN * TEMP / r)))

    a = RNoisy('p', 'n', r=2e3); a.update_iparv()
    b = R('p', 'n', r=2e3); b.update_iparv()
    x = np.zeros(2)
    ga = float(np.asarray(a.CY(x, 1.0), float)[0, 0])
    gb = float(np.asarray(b.CY(x, 1.0), float)[0, 0])
    print("The generated thermal-noise stamp is **%.4g A^2/Hz**, matching" % ga)
    print("the hand-written ``elements.R`` exactly (%.4g)." % gb)

Newton limiting via PCNR
````````````````````````

When ``pcnr=True``, limiting is done by the simulator's PCNR layer
(Predictor/Corrector Newton-Raphson) rather than by each device behind
the simulator's back. A generated element joins in automatically when
its current is an exponential of a single branch voltage: the compiler
reads the exponential scale off the expression and emits the whole PCNR
device protocol.

.. exec-rst::

    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, vt
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class MyDiode(Behavioural):
        instparams = [Parameter(name='IS', desc='Is', unit='A',
                                default=1e-13)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * (sympy.exp(b.V / vt()) - 1))

    got = np.asarray(MyDiode.pcnr_i(0.6, {'IS': 1e-13}, defaultepar,
                                    numeric), float)
    ref = np.asarray(Diode.pcnr_i(0.6, {'IS': 1e-13}, defaultepar,
                                  numeric), float)
    print("``MyDiode`` declares the junction %s and its generated PCNR"
          % (MyDiode.pcnr_junctions,))
    print("current at 0.6 V is **%.6g A** -- the hand-written ``Diode``" % got[0])
    print("gives %.6g A, i.e. identical." % ref[0])

An element that does *not* qualify — a polynomial nonlinearity, several
branches, a generated state — simply declares nothing and is solved the
ordinary way. The compiler never claims a capability it cannot honour.

.. admonition:: Use ``limexp`` in models you intend to run under PCNR

   PCNR bounds the *limited quantity*, but during assembly the layer
   still evaluates the device's own ``i()`` at the raw node voltage --
   it adds the whole circuit's current and then subtracts the
   participant's, because that current is about to be re-stamped at
   ``v_lim``.  The cancellation is exact in arithmetic and worthless in
   floating point once the term is ``inf``: ``inf - inf`` is ``nan`` and
   poisons the system.  Measured on a 20 V, 1 ohm forward drive, a
   generated diode written with a raw ``exp`` fails under PCNR while the
   ``limexp`` version converges to the same answer as the ordinary path.

What it costs
-------------

The generated code is a compiled numpy function, not an interpreter.
Regenerated at build time:

.. exec-rst::

    import time, warnings
    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import SubCircuit, VSin, R, C, Diode
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.transient import Transient

    cm.default_toolkit = numeric

    def bench(fn, n=4000):
        for _ in range(3):
            fn()
        t0 = time.perf_counter()
        for _ in range(n):
            fn()
        return (time.perf_counter() - t0) / n * 1e6

    rows = []
    for label, ref, hdl_el, meth, x in (
            ('R.i', R('p', 'n', r=1e3), eh.RHdl('p', 'n', r=1e3), 'i',
             np.array([0.3, 0.0])),
            ('C.C', C('p', 'n', c=1e-9), eh.CHdl('p', 'n', c=1e-9), 'C',
             np.array([0.3, 0.0])),
            ('Diode.i', Diode('p', 'n', IS=1e-13),
             eh.DiodeHdl('p', 'n', IS=1e-13), 'i', np.array([0.4, 0.0]))):
        ref.update_iparv(); hdl_el.update_iparv()
        a = bench(lambda: getattr(ref, meth)(x))
        b = bench(lambda: getattr(hdl_el, meth)(x))
        rows.append((label, a, b, b / a))

    print(".. list-table:: Per-call cost, hand-written vs generated")
    print("   :header-rows: 1")
    print("")
    print("   * - call")
    print("     - hand-written")
    print("     - generated")
    print("     - ratio")
    for label, a, b, r in rows:
        print("   * - ``%s``" % label)
        print("     - %.2f us" % a)
        print("     - %.2f us" % b)
        print("     - %.2fx" % r)
    print("")
    print("Compilation happens once, at class definition (~6-9 ms per")
    print("element).  For a *nonlinear* element the generated code is at")
    print("parity or better -- both sides evaluate the same transcendental,")
    print("and the generated one shares subexpressions.  The gap is on")
    print("trivial linear stamps, where the hand-written element does almost")
    print("nothing at all and any call overhead dominates a 2x2 matrix.")

End to end, an 8-stage RC ladder transient built from generated elements
runs at about 1.14x the hand-written one and produces identical waveforms
(``benchmarks/hdl_overhead.py`` measures both and asserts the agreement).

Where this sits
---------------

Behavioural modelling across simulators falls into four families:

.. list-table::
   :header-rows: 1
   :widths: 26 30 22 22

   * - family
     - examples
     - dynamics
     - Jacobian
   * - static expression source
     - ngspice B-source
     - none
     - symbolic diff of the parse tree
   * - rich expression source
     - Xyce B, HSPICE/PSpice E/G, LTspice B
     - ``ddt``/``sdt``, ``idt`` (LTspice)
     - forward AD over the AST
   * - per-branch I(V)+Q(V)
     - Qucs EDD
     - charge-based
     - symbolic, ≤20 branches
   * - compiled HDL
     - Verilog-A (Spectre, OpenVAF), XSPICE
     - full operator set
     - compiler AD

pycircuit's ``Behavioural`` is a hybrid of the middle two with the
surface syntax of the last. Its distinguishing properties, in decreasing
order of how much they matter:

* ``ddt`` inherits the simulator's integration order, because it declares
  a charge instead of computing a derivative. Xyce's expression-level
  ``ddt`` is hardwired to backward Euler regardless of the integrator in
  use, with a comment in its source saying so; LTspice's silently
  degenerates to the identity in ``.DC``; ngspice has none.
* A non-conservative ``ddt`` scaling is **rejected**, not integrated
  anyway.
* ``idt``/``idtmod`` become real DAE states rather than hidden
  accumulators — which is why they behave at DC and under ``uic``.
* Noise contributions are supported; expression sources are noiseless.
* Exact closed-form Jacobians over the whole element, computed once.
* There is a vectorised (JAX) evaluation path — no simulator in the
  survey vectorises behavioural evaluation across devices.

And the honest other side: no ``$limit``, no ``laplace_*``, no
``transition``/``slew``, no events yet — and, in common with every
expression-based approach, **a contribution must be valid for every value
the solver might try, not only for physical ones**. An expression with a
square root or a division that a non-physical iterate can reach will
produce ``nan`` and a "timestep too small" failure rather than a helpful
message.

Extending the language
----------------------

The development plan lives in ``hdl.md`` section 9, in four phases: widen
PCNR participation (charge-storing junctions are done; multi-junction
devices are next), ``$limit``, then the cheap high-frequency model
surface (``$param_given``, parameter ranges, node collapse, AC
excitation), then events and symbolic-toolkit support. Deferred items
each name the trigger that would reopen them.

Two conventions bind any addition, both learned by getting them wrong:

#. **An operator is not finished until a test exercises it inside a
   circuit.** ``limexp`` once passed a standalone test of the sympy
   function while being unusable in *every* element — the atoms carried
   no ``real`` assumption, so sympy refused to build the comparison in
   its ``Piecewise``.
#. **Refuse per element what the compiler cannot honour.** Qualification
   rules that emit nothing are better than a plausible-looking stamp; a
   silent wrong answer is the expensive failure.

Reference
---------

.. autoclass:: pycircuit.circuit.hdl.Behavioural

.. autoclass:: pycircuit.circuit.hdl.Contribution

The worked catalogue is ``pycircuit/circuit/elements_hdl.py``: ten
elements, each exercising one capability, each proven equivalent to its
hand-written counterpart in
``pycircuit/circuit/tests/test_elements_hdl.py``.

Further reading
---------------

* ``hdl.md`` — research record, capability map against the Verilog-AMS
  LRM with call-site counts from 78 real device models, measured cost,
  the development plan, and the defects found while building it.
* :ref:`idtmod` — the circular integrators the ``idt``/``idtmod``
  operators lower onto.
* ``benchmarks/hdl_overhead.py`` — the cost measurement.
