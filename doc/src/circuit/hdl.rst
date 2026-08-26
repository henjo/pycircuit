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
(``$abstime``), ``param_given`` (``$param_given``), ``ac_stim``,
``limit_pnj``/``limit_fet``/``limit_vds`` (``$limit``),
``limit_together`` (those probes limited as one device), and any sympy function — ``exp``, ``log``,
``tanh``, ``sqrt``, ``Piecewise`` for conditionals.  Two statements
besides ``Contribution``: ``Cross(expr, direction)`` for ``@cross``
events, and ``Contribution(b.V, 0)`` which **collapses** a branch,
merging an internal node away.

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

The math kernel
```````````````

``exp``, ``log`` and ``sqrt`` are enough for a behavioural element and
not enough for a compact model, for two reasons that only show up once
the element is differentiated.

**Both arms of every conditional are evaluated.** A compiled
``Piecewise`` becomes a ``where``, which selects *afterwards*.  The
discarded arm still takes the square root of its negative number, still
divides by its zero, and still contributes its NaN to anything that
combines the two.

**Clamping fixes the value and breaks the derivative.**
``sqrt(Max(x, 0))`` differentiates to ``d(Max)/dx / (2*sqrt(Max(x,0)))``
— which is ``0/0`` for every ``x <= 0``.  The current comes out right and
the Jacobian comes out NaN, which is the worse failure: nothing looks
wrong and Newton is silently poisoned.

:mod:`pycircuit.circuit.hdl` therefore ships a small kernel of *smoothed*
replacements, all of which are finite in value **and** in derivative
everywhere:

.. exec-rst::

    import warnings
    import numpy as np, sympy
    warnings.simplefilter('ignore')
    from pycircuit.circuit import hdl

    x = sympy.Symbol('x')
    rows = [('sqrt(Max(x,0))', sympy.sqrt(sympy.Max(x, 0))),
            ('safe_sqrt(x)',   hdl.safe_sqrt(x)),
            ('log(x)',         sympy.log(x)),
            ('safe_ln(x)',     hdl.safe_ln(x)),
            ('1/x',            1 / x),
            ('safe_div(1, x)', hdl.safe_div(1.0, x))]
    mods = [{'Max': np.maximum, 'Min': np.minimum,
             'Heaviside': lambda a, *r: np.heaviside(a, 0.5)}, 'numpy']

    def ev(g):
        try:
            v = float(g(0.0))
            return '%.4g' % v if np.isfinite(v) else '**not finite**'
        except Exception:
            return '**raises**'

    print(".. list-table:: Evaluated at x = 0, where every one of these "
          "is used")
    print("   :header-rows: 1")
    print("")
    print("   * - expression")
    print("     - value")
    print("     - derivative")
    for name, e in rows:
        f = sympy.lambdify(x, e, mods)
        d = sympy.lambdify(x, sympy.diff(e, x), mods)
        print("   * - ``%s``" % name)
        print("     - %s" % ev(f))
        print("     - %s" % ev(d))

The family is ``expl`` / ``expl_low`` / ``expl_high`` (bounded
exponentials, PSP's, continued by their own third-order Taylor and
``C-3`` at the seams), ``hypsmooth`` (a smooth ``max(x, 0)``, and the
right first thing to reach for), ``safe_sqrt``, ``safe_ln``,
``safe_div``, ``safe_abs``, and ``limexp`` (the LRM's, deliberately left
unclamped so PCNR can recognise it).

**They are not free, and the cost is exponent range.**  Each regulariser
squares something: ``safe_div`` squares its ``eps``, ``hypsmooth`` and
``safe_abs`` square theirs, and ``safe_ln`` inherits ``hypsmooth``\ 's.
Squaring halves the range the function can carry, so

.. exec-rst::

    import warnings
    import numpy as np, sympy
    warnings.simplefilter('ignore')
    from pycircuit.circuit import hdl

    x = sympy.Symbol('x')
    f = sympy.lambdify(x, hdl.safe_ln(x), 'numpy')
    g = sympy.lambdify(x, sympy.log(x), 'numpy')
    print(".. list-table:: ``safe_ln`` against plain ``log``, large "
          "argument")
    print("   :header-rows: 1")
    print("")
    print("   * - x")
    print("     - ``safe_ln(x)``")
    print("     - ``log(x)``")
    for v in (1e150, 1e160, 4.4e222):
        a, b = float(f(v)), float(g(v))
        fmt = lambda q: ('%.6g' % q) if np.isfinite(q) else '**inf**'
        print("   * - %.1e" % v)
        print("     - %s" % fmt(a))
        print("     - %s" % fmt(b))

So **ask what a guard is protecting against, and whether that case is
reachable.**  If the argument is provably positive — ``1 + expl(...)``
is at least 1 by construction — the plain function is both correct and
wider, and reaching for the safe one buys nothing while spending range
you may need.  This is not hypothetical: it cost a compact model here a
drain current of ``-inf`` at a bias a diverging Newton step reaches,
from an expression whose every ingredient was finite.

Two related traps, both from the same model:

* a **product** of two bounded exponentials is not bounded.  ``expl``
  continues *polynomially* past its seam rather than saturating, so two
  factors of ~1e186 overflow.  Where only the logarithm of the product
  is needed, carry the exponents instead: ``log(1 + exp(z))`` is bounded
  by ``|z| + 0.7``;
* a smoothing **constant** needs a scale.  ``hypsmooth(x, 1e-3)`` whose
  result is divided by a model parameter depends on the ratio, and that
  parameter varied by 44000× between two cards of the same process.
  Write the width relative to the scale it is compared against.

``hdl.md`` §3.2c has the full kernel table and the measurements behind
these.

Parameters
``````````

Declare them with the ordinary
:class:`~pycircuit.utilities.param.Parameter`; inside ``analog()`` the
names are bound as symbols, and values are read from the *resolved*
parameters at call time, so hierarchical expressions (``r='2*rbase'``)
reach the generated code.  Parameters may declare a range
(``minval``/``maxval``, Verilog-A's ``from [lo:hi]``), which is checked
both when the value is set and after a string expression resolves, and a
class may declare ``aliasparams = {alias: canonical}`` for the
alternative spellings compact models carry by the dozen.

``param_given('rs')`` answers Verilog-A's ``$param_given``: did the user
supply this parameter, as opposed to letting it default?  It is a
*runtime* value -- the element compiles once per class, while givenness
belongs to each instance -- so use it in a ``Piecewise``:

.. exec-rst::

    import sympy
    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       param_given)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class Rsel(Behavioural):
        instparams = [Parameter(name='g1', desc='g', unit='S', default=1e-3),
                      Parameter(name='g2', desc='g', unit='S', default=5e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            g = sympy.Piecewise((g2, param_given('g2') > 0.5), (g1, True))
            return Contribution(b.I, g * b.V)

    x = np.array([1.0, 0.0])
    a = Rsel('p', 'n'); a.update_iparv()
    b = Rsel('p', 'n', g2=5e-3); b.update_iparv()
    print("Defaulted, the element conducts %.0e S; given ``g2=5e-3`` -- the"
          % float(np.asarray(a.i(x), float)[0]))
    print("same value it would have defaulted to -- it conducts %.0e S."
          % float(np.asarray(b.i(x), float)[0]))
    print("Givenness is not \"differs from the default\", which is exactly")
    print("why the operator exists.")

Parameters as a namespace
`````````````````````````

The bare-name convention above has three limits, all measured on the
shipped library.  A **helper shared by several classes** cannot read the
parameters -- the compiler binds them into a private copy of
``analog()``'s globals, and a helper has its own -- so the SPICE diode
body took nineteen explicit arguments and the Gummel-Poon body was
callable only through a hand-applied rebind of its globals.  A **Python
keyword** cannot be a bare name, so SPICE's ``LAMBDA`` and ``AS`` ship
on `MosLevel1Hdl` as ``lambd`` and ``asrc``.  And a linter sees every
bare parameter as an undefined name (63 ``# noqa: F821`` in one file).

``params_as = 'p'`` on the class hands ``analog()`` a
:class:`~pycircuit.circuit.hdl.ParamNamespace` as its **first argument**
instead of binding bare names; the terminals are the arguments after it.
``p.g`` and ``p['g']`` are the *same* symbol the bare style binds, so
the generated code, the compilation record and ``explain()`` are
identical between the two spellings (``test_hdl_params.py`` pins that,
and pins the five library adopters to hashes recorded before they were
converted).  What changes is scope and spelling:

.. exec-rst::

    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    def channel(p, d, s):
        """Shared by both classes below; reads whichever card called it."""
        b = Branch(d, s)
        clm = var(1 + p['lambda'] * b.V, 'clm')       # a keyword-named parameter
        return Contribution(b.I, p.kp * clm * b.V)

    class Nch(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='kp', desc='k', unit='A/V^2', default=2e-5),
                      Parameter(name='lambda', desc='clm', unit='1/V', default=0.0)]

        @staticmethod
        def analog(p, d, s):
            return channel(p, d, s)

    class Pch(Behavioural):
        params_as = 'p'
        instparams = [Parameter(name='kp', desc='k', unit='A/V^2', default=8e-6),
                      Parameter(name='lambda', desc='clm', unit='1/V', default=0.0)]

        @staticmethod
        def analog(p, s, d):
            return channel(p, d, s)

    x = np.array([1.0, 0.0])
    n = Nch('d', 's', **{'lambda': 0.5}); n.update_iparv()
    q = Pch('s', 'd'); q.update_iparv()
    print("Terminals are the arguments after ``p``: Nch %s, Pch %s."
          % (list(Nch.terminals), list(Pch.terminals)))
    print("One helper, two cards: Nch conducts %.1e A at 1 V with lambda = 0.5,"
          % float(np.asarray(n.i(x), float)[0]))
    print("Pch %.1e A with its own kp and lambda = 0." % float(np.asarray(q.i(x), float)[0]))

``p.given('rs')`` is `param_given` with the name checked against the
declaration, and ``p.names`` is the declared order.  A parameter that
is not declared is refused **at compile time**, naming the class and
the declared names, exactly as an undeclared bare name is; a bare name
inside a ``params_as`` class gets a message saying where the parameters
went.  A keyword-named parameter can only be declared by a
``params_as`` class -- the bare style refuses it with the fix in the
message -- and in generated source it appears under a mangled name
(``_hdl_kw_lambda``), since no Python function can take an argument
called ``lambda``; ``explain()`` lists it by its declared name.

The two styles coexist across classes, not inside one: every model that
reads bare names keeps working unchanged, and a ``params_as`` class
binds no bare names at all.  `MosLevel1Hdl` still reads bare names and
keeps ``lambd``/``asrc`` as the canonical spellings; SPICE's own are
accepted on the instance through ``aliasparams``.

Conditionals: both arms are evaluated
-------------------------------------

A ``Piecewise`` compiles to ``where``, which **picks afterwards**.  Both
arms run at every bias, so an arm that is not selected still overflows,
still raises floating-point flags, and still has its derivative taken --
and the derivative is the half that matters, because the chain rule
multiplies a zero partial by that arm's derivative and ``0 * nan`` is
``nan``.

The model obeys that rule by hand, with a clamped copy of the input fed
to the arm that does not want it::

    sympy.Piecewise((0.0, p <= 0.0), (1.0 / maxc(p, 1e-30), True))

That ``maxc`` is domain propagation a compiler can do.  The second arm is
selected exactly where ``p > 0``, so inside it ``p`` may be replaced by
anything that equals ``p`` there.  :func:`~pycircuit.circuit.hdl.select`
makes the substitution::

    hdl.select((0.0, p <= 0.0), (1.0 / p, True), margin=1e-30)

For each arm, the region it is selected in -- its own condition **and**
the complement of every earlier one -- is turned into a box of bounds on
individual symbols, and each bounded symbol is replaced inside that arm
by a ``minc``/``maxc`` clamped copy.

**The clamp is the identity where the arm is selected, in value and in
derivative, exactly.**  ``minc(v, c)`` returns ``v`` bit for bit
wherever ``v <= c``, and differentiates to ``_step(c, v) = 1`` there,
with exactly zero weight on ``c`` -- so a converted arm returns the same
bits from ``i``, ``G``, ``q`` and ``C`` as the hand-clamped one.  The
tie ``v == c`` is included on purpose: it is a selected point whenever
the condition is non-strict.

Supported conditions, and nothing else:

===============================  =====================================
``v < c``, ``v <= c``            ``v -> minc(v, c)``
``v > c``, ``v >= c``            ``v -> maxc(v, c)``
``a < v``                        the mirror of the above
``u < v`` (both atoms)           both, simultaneously
``|v| < c``, ``v*v < c``         ``v -> minc(maxc(v, -c), c)``
``|v| > c``                      pushed out of the hole
``And(...)``                     the intersection of the boxes
``Or(...)``, negated             ditto, by De Morgan
``True``                         the complement of the earlier arms
===============================  =====================================

``c`` may be any expression, including one that mentions ``v`` itself
(``vds < vdsat``): the substitution is simultaneous, so each bound is
evaluated at the unclamped symbol and no cycle is created.  It costs
nothing in finiteness either -- the bound was already being evaluated,
by the condition.

Everything else is **refused**, loudly, with ``SelectRefused`` naming
the condition and the reason: an ``Or`` in a positive position (a union
of boxes is not a box), an equality, a relational whose sides are both
compound.  A ``select`` that quietly did not clamp would be worse than
the ``Piecewise`` it replaced, because the author would trust it.  The
refusal fires only when the arm actually *uses* one of the symbols the
condition failed to bound -- a refusal over a condition the arm does not
read is noise, and noise is how a strict mode gets turned off.
``strict=False`` accepts any shape unclamped, and ``hdl.unclamped(expr)``
marks one arm as deliberately raw.

``margin`` (default ``0.0``) is subtracted from every bound that came
from a **strict** inequality.  It exists because clamping to the
boundary is not always enough: ``p > 0`` clamps to ``maxc(p, 0)``, and
an arm that divides by ``p`` still divides by zero.  The price is that
the clamp is no longer the identity on the sliver ``0 < p < margin`` --
which is why the default is zero, so that a selected value can never
change by accident.  A margin is never applied to a non-strict bound:
there the boundary is *inside* the selected region, where moving it
would change an answer.

.. exec-rst::

    import numpy as np
    import sympy
    from pycircuit.circuit import hdl

    v = sympy.Symbol('v', real=True)
    e = hdl.select((sympy.log(v - 1.0), v > 1.0), (0.0, True), margin=1e-6)
    print("The arm the model wrote is ``log(v - 1)``; what is compiled is::\n")
    print("    %s\n" % e.args[0].expr)
    f = sympy.lambdify(v, e.args[0].expr, 'numpy')
    g = sympy.lambdify(v, sympy.diff(e.args[0].expr, v), 'numpy')
    x = np.array([-1e30, 0.5, 2.0])
    print("Discarded at ``v = -1e30`` and ``v = 0.5``, where the unclamped")
    print("arm is ``nan``, it is %s, with derivative %s;" % (f(x)[:2], g(x)[:2]))
    print("selected at ``v = 2`` it is exactly ``log(1) = %s``." % f(x)[2])

What it cannot do
`````````````````

A clamp can only reach a symbol the arm **mentions**.  If the arm uses
an intermediate that ``var()`` bound earlier from the same symbol --
``vdsc = leff*vmax/us``, used in an arm conditioned on ``vmax > 0`` --
substituting ``vmax`` inside the arm does not change ``vdsc``, and the
clamp is a no-op.  ``select`` detects exactly that case, and refuses.

That is the binding limit, not a corner case.  Surveying the 33
``Piecewise`` sites in ``elements_hdl.py`` and the 24 in
``psp_kernel.py``: most arms are built from ``var()`` intermediates and
never name the symbol their own condition constrains, so the clamp has
nothing to substitute into.  Automatic domain propagation and automatic
intermediate holding are the same problem seen twice, and a version of
``select`` that reached inside the let-chain would have to clamp a
*definition*, not an arm.

Two smaller limits, both refusals rather than silence:

* **a compound side.**  ``area*ikf <= 0`` constrains neither factor, and
  clamping the product would need a substitution for the compound --
  which sympy does not keep: ``1.0/(area*ikf)`` is stored as
  ``1.0*area**-1*ikf**-1``, in which ``area*ikf`` is not a node, so the
  substitution would silently do nothing.  Binding the product with
  ``var()`` first gives the condition an atom;
* **a singular arm at a selected point.**  ``select((sqrt(v), v >= 0),
  ...)`` clamps to ``sqrt(maxc(v, 0))``, whose derivative is ``0/0``
  below zero -- and no margin can fix it, because ``v = 0`` is a point
  the arm is *selected* at.  That is not a clamping problem: smooth
  instead of clip, with `safe_sqrt`.

Writing your first element
--------------------------

The whole loop on a real element -- a transconductor that saturates and
carries its own output capacitance -- from the equations to a solved
circuit.  Four rules cover almost everything:

#. ``analog()`` is a **staticmethod**, and its argument names *are* the
   terminals, in pin order;
#. it **returns** its contribution statements (Verilog-A has no
   ``return``; this is Python and you need one);
#. every expression is **symbolic** -- ``sympy.tanh``, not
   ``math.tanh`` -- because the compiler differentiates it;
#. parameters are declared with
   :class:`~pycircuit.utilities.param.Parameter` and used by bare name.

.. exec-rst::

    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ddt, check_jacobians)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class SoftGm(Behavioural):
        instparams = [
            Parameter(name='gm', desc='Transconductance', unit='S',
                      default=1e-3),
            Parameter(name='vsat', desc='Saturation voltage', unit='V',
                      default=0.2),
            Parameter(name='cout', desc='Output capacitance', unit='F',
                      default=1e-12)]

        @staticmethod
        def analog(inp, inn, outp, outn):
            bin_ = Branch(inp, inn)
            bout = Branch(outp, outn)
            return (Contribution(bout.I,
                                 gm * vsat * sympy.tanh(bin_.V / vsat)),
                    Contribution(bout.I, ddt(cout * bout.V)))

    el = SoftGm('ip', 'in', 'op', 'on', gm=2e-3, vsat=0.2)
    el.update_iparv()
    print(".. code-block:: text")
    print("")
    for line in el.explain(source=False, symbolic=False).splitlines():
        print("   " + line)

``explain()`` is the first thing to call on a new element, because the
one thing you cannot read off the source is the **x-vector layout**: the
order of the unknowns that ``i``, ``q``, ``G`` and ``C`` use.  Terminals
come first, then internal nodes, then ``idt``/``idtmod`` states -- which
are not in your source at all -- then branch currents.  Reconstructing
that order by hand from ``el.nodes`` and ``el.branches`` is how a model
written this way earned an ``IndexError: index 5 is out of bounds for
axis 0 with size 5``.

The second thing to call is ``check_jacobians()``, which differentiates
``i`` and ``q`` numerically and compares them against the ``G`` and
``C`` the compiler derived, and scans everything for non-finite entries:

.. exec-rst::

    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ddt, check_jacobians)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class SoftGm(Behavioural):
        instparams = [
            Parameter(name='gm', desc='Transconductance', unit='S',
                      default=1e-3),
            Parameter(name='vsat', desc='Saturation voltage', unit='V',
                      default=0.2),
            Parameter(name='cout', desc='Output capacitance', unit='F',
                      default=1e-12)]

        @staticmethod
        def analog(inp, inn, outp, outn):
            bin_ = Branch(inp, inn)
            bout = Branch(outp, outn)
            return (Contribution(bout.I,
                                 gm * vsat * sympy.tanh(bin_.V / vsat)),
                    Contribution(bout.I, ddt(cout * bout.V)))

    el = SoftGm('ip', 'in', 'op', 'on', gm=2e-3, vsat=0.2)
    el.update_iparv()
    print(".. code-block:: text")
    print("")
    for line in repr(check_jacobians(el, [0.05, 0.0, 1.0, 0.0])).splitlines():
        print("   " + line)

Then use it like any other element.  A 0.5 V input is 2.5 saturation
voltages in, so the output current is well into ``tanh``'s flat region --
which is the point of the model and the reason to check it against the
closed form rather than against itself:

.. exec-rst::

    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, ddt)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class SoftGm(Behavioural):
        instparams = [
            Parameter(name='gm', desc='Transconductance', unit='S',
                      default=1e-3),
            Parameter(name='vsat', desc='Saturation voltage', unit='V',
                      default=0.2),
            Parameter(name='cout', desc='Output capacitance', unit='F',
                      default=1e-12)]

        @staticmethod
        def analog(inp, inn, outp, outn):
            bin_ = Branch(inp, inn)
            bout = Branch(outp, outn)
            return (Contribution(bout.I,
                                 gm * vsat * sympy.tanh(bin_.V / vsat)),
                    Contribution(bout.I, ddt(cout * bout.V)))

    c = SubCircuit()
    nin, nout = c.add_nodes('in', 'out')
    c['vs'] = VS(nin, gnd, v=0.5)
    c['gm'] = SoftGm(nin, gnd, gnd, nout, gm=2e-3, vsat=0.2)
    c['rl'] = R(nout, gnd, r=1e3)
    v = float(DC(c, toolkit=numeric).solve().v('out'))
    closed = 2e-3 * 0.2 * np.tanh(0.5 / 0.2) * 1e3
    print("Driven at 0.5 V into a 1 k load the element gives ``v(out) =")
    print("%.6f V``; the closed form ``gm*vsat*tanh(vin/vsat)*RL`` is" % v)
    print("%.6f V. The compiler derived the Jacobian of that ``tanh``," % closed)
    print("and the DC solve converged on it.")

When it goes wrong
------------------

The compiler defends its *semantic* invariants closely -- a switch
branch, a state-scaled ``ddt``, an operating-point-dependent collapse
and an ``instparams`` drift are all refused with a message saying what
you wrote, why it cannot work, and what to write instead.  Its
*syntactic* surface used to have no messages at all, and three of the
mistakes below were **silent**: they produced a working element that was
wired differently from the one you wrote.

Every message in this table is produced by actually making the mistake,
here, when this page is built -- so it cannot go stale.  Only the first
sentence is shown; each message continues with the fix.

.. exec-rst::

    import math
    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Node,
                                       Contribution)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric
    P = [Parameter(name='g', desc='g', unit='S', default=1e-3)]

    def build(name, analog, **extra):
        d = dict(instparams=list(P), analog=staticmethod(analog))
        d.update(extra)
        return type(name, (Behavioural,), d)

    def no_return(plus, minus):
        Contribution(Branch(plus, minus).I, g * Branch(plus, minus).V)

    def with_self(self, plus, minus):
        return Contribution(Branch(plus, minus).I, g * Branch(plus, minus).V)

    def good(plus, minus):
        return Contribution(Branch(plus, minus).I, g * Branch(plus, minus).V)

    def collide(plus, out):
        n = Node('out')
        return (Contribution(Branch(plus, n).I, g * Branch(plus, n).V),
                Contribution(Branch(n, out).I, g * Branch(n, out).V))

    def python_if(plus, minus):
        b = Branch(plus, minus)
        if b.V > 0:
            return Contribution(b.I, g * b.V)
        return Contribution(b.I, 0)

    def float_math(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * math.exp(b.V / 0.026))

    def undeclared(plus, minus):
        return Contribution(Branch(plus, minus).I,
                            gain * Branch(plus, minus).V)

    cases = [
        ("no ``return`` in ``analog()``", "compiler traceback",
         lambda: build('NoReturn', no_return)),
        ("``def analog(self, plus, minus)``", "**silent**",
         lambda: build('WithSelf', with_self)),
        ("``terminals`` in the class body", "**silent**",
         lambda: build('Decl', good, terminals=('minus', 'plus'))),
        ("``Node('out')`` beside an ``out`` terminal", "**silent**",
         lambda: build('Collide', collide)),
        ("``if b.V > 0``", "raw sympy",
         lambda: build('IfOnV', python_if)),
        ("``math.exp(b.V)``", "raw sympy",
         lambda: build('FloatMath', float_math)),
        ("a parameter no class declared", "``NameError`` at first call",
         lambda: build('Undeclared', undeclared)),
        ("three nodes for a two-pin element", "**silent**",
         lambda: build('Ok1', good)('a', 'b', 'c')),
        ("``R=1e3`` where the parameter is ``g``", "bare ``KeyError``",
         lambda: build('Ok2', good)('a', 'b', R=1e3)),
    ]

    def first_sentence(text):
        text = ' '.join(text.split())
        for k in range(30, len(text)):
            if text[k] == '.' and (k + 1 == len(text) or text[k + 1] == ' '):
                return text[:k + 1]
        return text

    print(".. list-table::")
    print("   :header-rows: 1")
    print("   :widths: 26 16 58")
    print("")
    print("   * - what you wrote")
    print("     - used to be")
    print("     - what it says now")
    for label, before, fn in cases:
        try:
            fn()
        except Exception as e:
            msg = first_sentence(str(e))
        else:
            msg = "NO ERROR -- this table is out of date"
        print("   * - %s" % label)
        print("     - %s" % before)
        print("     - %s" % msg)

Two instruments
```````````````

``explain(element)`` prints the compilation record: terminals, the
x-vector layout, the parameters *as the generated code will read them*
(so a hierarchical expression nothing has resolved is visible rather
than silently defaulted), which compilation path the model took, what
the compiler found in it -- states, collapses, crossings, ``$limit``
probes, PCNR junctions, constant stamps, JAX pure forms -- the symbolic
``i``/``q``/``G``/``C``, and the generated source.

``check_jacobians(element, x)`` answers the question the NaN warnings
above raise and never answer.  This is the classic one -- a depletion
charge written the way the textbook writes it, which is finite over the
whole physical range and ``inf`` at its edge:

.. exec-rst::

    import numpy as np, sympy
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ddt, check_jacobians)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class Varicap(Behavioural):
        instparams = [Parameter(name='cj0', desc='C', unit='F',
                                default=1e-12),
                      Parameter(name='vj', desc='V', unit='V', default=0.8)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            q = 2 * cj0 * vj * (1 - sympy.sqrt(1 - b.V / vj))
            return Contribution(b.I, ddt(q))

    el = Varicap('p', 'n')
    el.update_iparv()
    print(".. code-block:: text")
    print("")
    for v in (0.0, 0.5, 0.79, 0.8, 1.2):
        r = check_jacobians(el, [v, 0.0])
        print("   v = %5.2f   ok = %-5s   C[0,0] = %.6g"
              % (v, r.ok, r.results['C']['ana'][0, 0]))
    print("")
    print("and at the first bias that fails:")
    print("")
    print(".. code-block:: text")
    print("")
    for line in repr(check_jacobians(el, [0.8, 0.0])).splitlines():
        print("   " + line)

Nothing about that element looks wrong, it solves happily at every bias
a physical circuit visits, and a Newton iterate that overshoots
``vj`` turns the whole matrix to ``nan``.  ``safe_sqrt``, or the
junction's own ``fc``-limited formulation, is the fix; the instrument is
what tells you which bias to look at.

The finite-difference step is ``max(1e-7, 1e-7*|x[k]|)`` per column, and
the tolerance is **one band for the whole matrix** -- ``atol`` defaults
to ``1e-7`` of its largest entry.  That is deliberate: a per-entry
relative tolerance passes an entry that is small *because it is wrong*.
Before trusting a pass, tighten ``rtol`` and confirm a deliberately
broken model fails.

A third verdict: UNRESOLVED
```````````````````````````

A finite difference is not an oracle, and three mechanisms make it
report a large discrepancy on a model that is right: a **kink**, where
central differencing returns the average of two one-sided slopes and no
``h`` helps because a jump has no scale; **roundoff**, where the entry's
signal is below the representable step of the value it differentiates
and the difference comes back exactly ``0.0``; and **truncation** on a
stiff card.  Each was hit by a real model here before it was understood.

Entries like that come back ``UNRESOLVED`` rather than ``FAILED``.  The
floors are measured per entry, not assumed -- truncation by differencing
the same column at ``h`` and ``2h``, a kink by whether the one-sided
disagreement shrinks with ``h``, roundoff by widening the step until the
value clears its own quantisation.  So:

* ``res.ok`` is unchanged: no failing entry and nothing non-finite.  An
  unresolved entry does not make it false;
* ``res.resolved`` is the stronger claim -- every entry got a real
  verdict;
* ``res.verdict`` is ``'ok'``, ``'UNRESOLVED'``, ``'FAILED'`` or
  ``'NOT COMPARABLE'``;
* ``res.unresolved`` and ``res.failures`` list the entries, each with
  its analytic value, the difference, the error, the floor and the
  reason.

.. exec-rst::

    import numpy as np
    import pycircuit.circuit.circuit as cm
    from pycircuit.circuit.toolkit import numeric
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       maxc, check_jacobians)
    from pycircuit.utilities.param import Parameter

    cm.default_toolkit = numeric

    class Rectifier(Behavioural):
        instparams = [Parameter(name='g', desc='S', unit='S', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * maxc(b.V, 0.0))

    el = Rectifier('p', 'n')
    el.update_iparv()
    print(".. code-block:: text")
    print("")
    for line in repr(check_jacobians(el, [0.0, 0.0])).splitlines():
        print("   " + line)

**What that verdict does and does not claim.**  The non-smoothness is a
fact about the value, found without looking at the Jacobian at all.  At
a corner the derivative is genuinely two-valued, so an error *smaller
than the jump* cannot be separated from the ambiguity and is not
reported; an error larger than it is.  ``UNRESOLVED`` is the honest
statement of that, and it is why a model author should ask for
``res.resolved`` at biases where the model ought to be smooth.

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

A device may own **several** junctions -- two sharing a base, the BJT
shape -- and may store charge; both work on the CPU and traced backends,
because the solver asks the device for its junction rather than assuming
a textbook diode.  An element that does *not* qualify — a polynomial
nonlinearity, several branches, a generated state — simply declares
nothing and is solved the ordinary way. The compiler never claims a capability it cannot honour.

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
runs at about 1.15x the hand-written one and produces identical waveforms
(``benchmarks/hdl_overhead.py`` measures both and asserts the agreement).

.. note::

   This figure read ``1.14x`` until 2026-08-26 and had quietly gone
   stale: measured again on the same machine it was **1.22x** (median of
   four runs), because the argument list every generated call is made
   with was being rebuilt on every call -- 125 029 times in a 49-element
   transient.  Memoising it, invalidated by the same ``iparv`` observer
   that drops the constant stamps, brought the ladder back to 1.15x and
   gained 1.16x on a 121-device transient, with bit-identical answers.
   See ``benchmarks/args_cache.py`` and ``doc/hdl_roadmap_260824.md``
   sec. 24.

   Quote these from a median of several runs: single runs of the ladder
   benchmark ranged from 1.12x to 1.25x across the change.

The C backend
-------------

A model that uses ``var()`` compiles to straight-line scalar Python --
for the surface-potential MOSFET's Jacobian, some 21 000 interpreted
calls per evaluation.  The same statements, printed to C and compiled
by the system compiler, run two orders of magnitude faster and return
**the same bytes**.  That second backend ships with the DSL and is off
by default::

    PYCIRCUIT_HDL_BACKEND=c python my_simulation.py

or, from code, process-wide or per class::

    from pycircuit.circuit import hdl
    hdl.set_backend('c')                  # classes compiled from now on
    hdl.set_backend('c', MyModel)         # this class, immediately

``explain()`` prints which backend a class actually runs, and
``cls._hdl_backend_status`` says the same thing programmatically:
``'c'``, or ``'numpy (<why not>)'`` -- an eager-path class (no
``var()``), a missing compiler, a failed build.  A request for C that
cannot be served **falls back to numpy and warns**; it never silently
degrades and never breaks the class.

What it buys (measured 2026-08-26, gcc 15.2, one core;
``benchmarks/hdl_model_cost.py --backend`` reproduces the table):

.. list-table::
   :header-rows: 1

   * - model
     - ``i`` numpy
     - ``i`` C
     - ``G`` numpy
     - ``G`` C
   * - Gummel-Poon BJT
     - 112 us
     - 4.2 us (27x)
     - 347 us
     - 4.5 us (77x)
   * - MOS level 3
     - 205 us
     - 8.1 us (25x)
     - 2011 us
     - 8.1 us (248x)
   * - PSP103 (359 parameters)
     - 1288 us
     - 10.6 us (121x)
     - 17 349 us
     - 81.8 us (212x)

What it costs: the C text is printed at class-compile time whatever the
backend (+19% on a cold compile, then held by the compile cache), and
the shared objects are built by the system ``cc`` once per (source,
flags) -- for PSP about 90 s, in parallel across ``i``/``G``/``q``/``C``,
stored beside the compile cache and thereafter ``dlopen``-ed in
milliseconds.  A warm store even serves a machine with **no compiler at
all**.

The fidelity contract is bit-identity with the numpy path, and it is a
test, not an aspiration: every chained library element and the PSP
MOSFET are compared byte-for-byte over sweeps that include ``+-1e30``
(``test_hdl_cbackend.py``).  The compiler is told to keep every
``pow``, ``exp``, ``log`` a real libm call (``-fno-builtin-*``, no
fast-math, no FMA contraction), because gcc would otherwise fold
``pow(x, 2.0)`` into the *more accurate* ``x*x`` and drift an ulp from
numpy's glibc ``pow``.  Two measured exceptions, both named in the
tests: numpy ships its own vectorised ``tanh`` (an ulp from libm's in
~30% of arguments -- a ``tanh``-using model agrees bitwise with its
source run against libm's ``tanh``, and to that ulp with numpy's), and
the sign of exact zeros (numpy computes integer-typed subchains, whose
zeros carry no sign; values compare equal).  Where the numpy path
returns NaN, C returns the same NaN; where a degenerate parameter set
makes the numpy path raise ``ZeroDivisionError`` from a pure-parameter
division, C returns ``inf`` instead -- the one behavioural difference.

The kernel's both-arms guarantee survives the translation because
selection stays a function call: ``numpy.where(c, a, b)`` becomes
``_sel(c, a, b)``, whose arguments C evaluates -- both of them --
before the pick, exactly as Python does.

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

And the honest other side. ``$limit``, ``laplace_*`` and ``@cross`` are
now implemented — this sentence said otherwise until 2026-08-24, having
been written before they landed, and it is corrected here rather than
quietly deleted because a stale capability claim is worse than none: it
is trusted like a measurement.

What is genuinely absent is ``absdelay``, ``transition``/``slew``,
``zi_*``, ``@timer`` and an event-latched discrete state — the last two
of which are what a phase-frequency detector, a divider or a
sample-and-hold need. ``doc/hdl_roadmap_260824.md`` in the repository
plans them, with measured compile and evaluation costs.

And, in common with every expression-based approach, **a contribution
must be valid for every value the solver might try, not only for
physical ones**. An expression with a
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

.. autoclass:: pycircuit.circuit.hdl.ParamNamespace

.. autofunction:: pycircuit.circuit.hdl.select

.. autofunction:: pycircuit.circuit.hdl.unclamped

.. autofunction:: pycircuit.circuit.hdl.explain

.. autofunction:: pycircuit.circuit.hdl.x_layout

.. autofunction:: pycircuit.circuit.hdl.check_jacobians

The worked catalogue is ``pycircuit/circuit/elements_hdl.py``:
**twenty-six** elements as of 2026-08-25.  The first ten each exercise
one language capability and each is proven equivalent to its
hand-written counterpart in
``pycircuit/circuit/tests/test_elements_hdl.py``; the rest are the
Phase-6 library -- junction diode, quartz resonator, ferrite bead,
skin-effect resistor, comparator, memristor, Gummel-Poon bipolar (both
polarities, and a self-heating n-p-n), EKV 2.6 MOSFET, SPICE level-1
MOSFET and a Boyle-class opamp macromodel -- tested in
``test_elements_hdl_library.py`` through ``..._library4.py``.

Further reading
---------------

* :ref:`compact_models` — what the language is ultimately for: a
  production surface-potential MOSFET and a foundry MoM capacitor,
  written in this DSL and measured against the vendor's own compiled
  binaries.
* ``hdl.md`` — research record, capability map against the Verilog-AMS
  LRM with call-site counts from 78 real device models, measured cost,
  the development plan, and the defects found while building it.
* :ref:`idtmod` — the circular integrators the ``idt``/``idtmod``
  operators lower onto.
* ``benchmarks/hdl_overhead.py`` — the cost measurement.
