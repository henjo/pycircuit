The symbolic_poly toolkit
=========================

:mod:`pycircuit.circuit` selects its math backend through a *toolkit*: the same
analysis code runs numerically or symbolically depending on the toolkit passed
to it.  ``symbolic_poly`` is an experimental symbolic toolkit that works in
sympy's polynomial/rational *domains* instead of on free-form expressions.

It behaves like the stock ``symbolic`` toolkit but

* solves the MNA system **fraction-free** (``DomainMatrix.solve_den``) so
  intermediate results stay polynomials instead of the nested fractions
  ``Matrix.LUsolve`` produces — which lets it handle larger circuits, and
* exposes the result as a rational **transfer function** ``N(s)/D(s)`` with
  poles, zeros, gain and a fast numeric evaluator.

It is opt-in; the ``symbolic`` toolkit is unchanged.

.. note::

   Declare the frequency variable as ``sympy.Symbol('s', imaginary=True)`` (it
   represents :math:`j\omega`).  This lets conjugation resolve cleanly in the
   noise analysis and keeps results expressed in ``s``.

.. note::

   A *fully* symbolic circuit (every component a distinct symbol) has a network
   determinant with exponentially many terms — no method yields a compact
   closed form.  ``symbolic_poly`` wins when there are **few symbolic
   generators**; numeric component values with only ``s`` symbolic is the sweet
   spot and scales to circuits with many nodes.


Quick start
-----------

.. code-block:: python

    import sympy
    from pycircuit.circuit import symbolic_poly, use_toolkit, SubCircuit, R, C, VS, gnd
    from pycircuit.circuit.analysis_ss import AC

    s = sympy.Symbol('s', imaginary=True)
    R1, C1 = sympy.symbols('R C', positive=True)

    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        cir['R']  = R('in', 'out', r=R1)
        cir['C']  = C('out', gnd, c=C1)

    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)

    H = res.tf('out', gnd)          # the transfer function to v(out)
    H.canonical()                   # -> 1/(C*R*s + 1)
    H.poles()                       # -> {-1/(C*R): 1}
    H.dcgain()                      # -> 1


Transfer functions
------------------

An AC solve with ``symbolic_poly`` returns a ``CircuitResultACPoly`` that keeps
the numerator vector ``N`` and the shared denominator ``D`` (the network
determinant).  From it you get transfer functions without a swelling
``cancel``:

* ``res.tf(plus, minus=None)`` — the voltage transfer function to
  ``v(plus[, minus])``.
* ``res.tf_i(branch_or_term)`` — the current transfer function into a terminal
  or branch (transimpedance / current gain).
* ``res.poles()`` — the circuit poles, computed **once** from the shared
  denominator for the whole circuit.

Each ``res.tf*`` returns a :class:`~pycircuit.circuit.transferfunction.TransferFunction`:

.. code-block:: python

    H = res.tf('out', gnd)
    H.canonical()          # cancelled num/den expression
    H.num, H.den           # raw numerator / denominator
    H.poles()              # {root: multiplicity}
    H.zeros()
    H.dcgain()             # H(s=0)
    H.bode()               # a fast lambdified callable f(s)
    H.frequencyresponse([1e3, 1e6, 1e9])   # complex response at those Hz

Current transfer functions and input impedance:

.. code-block:: python

    # source current of the RC low-pass above
    res.tf_i('R.plus').canonical()     # -> C*s/(C*R*s + 1)
    res.tf_i('R.plus').zeros()         # -> {0: 1}   (no current at DC)

    # input impedance: drive with a unit current source, then Zin = v(node)
    from pycircuit.circuit import IS
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['IS'] = IS(gnd, 'in', iac=1)
        cir['R']  = R('in', gnd, r=R1)
        cir['C']  = C('in', gnd, c=C1)          # R || C
    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    res.tf('in', gnd).canonical()      # -> R/(C*R*s + 1)


Poles and zeros at scale
------------------------

Symbolic root finding has no closed form for degree :math:`\ge 5`.  When the
transfer function has **numeric** coefficients, pass ``numeric=True`` to find
the roots with ``numpy.roots`` — fast and reliable regardless of degree:

.. code-block:: python

    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['VS'] = VS('in', gnd, vac=1)
        prev = 'in'
        for i in range(6):                       # a 6-section RC ladder
            cir['R%d' % i] = R(prev, 'n%d' % i, r=1e3)
            cir['C%d' % i] = C('n%d' % i, gnd, c=1e-9)
            prev = 'n%d' % i

    res = AC(cir, toolkit=symbolic_poly).solve(s, complexfreq=True)
    res.poles(numeric=True)          # 6 complex poles, all on the negative real axis

``numeric=True`` raises ``ValueError`` if any coefficient is still symbolic —
substitute parameter values first.


Noise
-----

The noise analysis works with ``symbolic_poly`` unchanged.  Internally the
output noise power ``zm**T CY conj(zm)`` is kept as a single shared-denominator
rational ``N**T CY conj(N) / (D conj(D))`` instead of an
:math:`O(n^2)` sum of divided rationals — the same value, a far smaller
expression:

.. code-block:: python

    from pycircuit.circuit import Noise
    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        cir['vs'] = VS(1, gnd, vac=1)
        cir['R']  = R(1, 2, r=R1)
        cir['C1'] = C(2, gnd, c=C1)
        noise = Noise(cir, inputsrc='vs', outputnodes=('2', gnd), toolkit=symbolic_poly)
    res = noise.solve(s, complexfreq=True)
    sympy.simplify(res['Svnout'])    # -> -4*R*T*k/(C**2*R**2*s**2 - 1)


Selecting the toolkit
---------------------

Pass ``toolkit=`` to an analysis, or scope the construction-time default with
the ``use_toolkit`` context manager (the recommended replacement for assigning
the deprecated ``circuit.default_toolkit`` global):

.. code-block:: python

    from pycircuit.circuit import use_toolkit, symbolic_poly

    with use_toolkit(symbolic_poly):
        cir = SubCircuit()
        ...                     # circuits built here use symbolic_poly
    # default restored on exit

Why fraction-free: a live benchmark
-----------------------------------

The table below is **generated when this page is built**.  It runs an AC solve of
an :math:`N`-section RC ladder with *numeric* component values and a *symbolic*
complex frequency :math:`s` -- the regime ``symbolic_poly`` targets -- under the
stock ``symbolic`` toolkit (``LUsolve``) and ``symbolic_poly``
(``DomainMatrix.solve_den``), reporting the solve time and the size of the
resulting node-voltage expression (``sympy.count_ops``).

Fraction-free elimination keeps intermediate results as a single
:math:`N(s)/D(s)` instead of the nested fractions ``LUsolve`` accumulates, so
the expression stays small and the solve is faster as :math:`N` grows.

.. exec-rst::

    import time, sympy
    from pycircuit.circuit import SubCircuit, R, C, VS, gnd, symbolic, symbolic_poly
    from pycircuit.circuit.analysis_ss import AC

    def run(N, tk):
        s = sympy.Symbol('s', imaginary=True)
        cir = SubCircuit(toolkit=tk)
        cir['VS'] = VS('n0', gnd, vac=1)
        for i in range(N):
            cir['R%d' % i] = R('n%d' % i, 'n%d' % (i + 1), r=100.0 * (i + 1))
            cir['C%d' % i] = C('n%d' % (i + 1), gnd, c=1e-9 * (i + 1))
        t = time.time()
        res = AC(cir, toolkit=tk).solve(s, complexfreq=True)
        v = res.v('n%d' % N, gnd)
        return time.time() - t, sympy.count_ops(v)

    print(".. list-table:: AC solve of an N-section RC ladder"
          " (numeric R,C + symbolic s)")
    print("   :header-rows: 1")
    print("   :widths: 8 16 14 18 14")
    print("")
    print("   * - N")
    print("     - symbolic time")
    print("     - symbolic ops")
    print("     - symbolic_poly time")
    print("     - symbolic_poly ops")
    for N in (2, 4, 6, 8):
        ta, oa = run(N, symbolic)
        tb, ob = run(N, symbolic_poly)
        print("   * - %d" % N)
        print("     - %.3f s" % ta)
        print("     - %d" % oa)
        print("     - %.3f s" % tb)
        print("     - %d" % ob)

.. note::

   The win is specific to **few symbolic generators** (numeric components with a
   symbolic ``s``, the usual AC/noise case).  For *fully* symbolic circuits (all
   ``R``/``C`` symbolic) the fraction-free determinant is a large multivariate
   polynomial and ``symbolic_poly`` offers no advantage -- use the stock
   ``symbolic`` toolkit there.
