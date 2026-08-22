Example 10: MFB filter -- symbolic post-processing with GiNaC
-------------------------------------------------------------

.. image:: mfb.*

This example revisits the Multiple-FeedBack (MFB) active filter of
:doc:`Example 2 <example2>`, which derives the transfer function, DC gain,
poles, natural frequency :math:`\omega_0` and quality factor :math:`Q` by hand
with sympy (the classic "Symby" analysis). Here we obtain **the same results**
through the native GiNaC backend: the linear solve runs in ``GinacToolkit`` and
the post-processing -- collecting the transfer function into ``N(s)/D(s)`` and
reading off the filter parameters -- uses :func:`pycircuit.circuit._ginac.poly_coeffs`,
which does the cancellation and collection in GiNaC and returns only the compact
coefficient vectors.

.. literalinclude:: example10.py

Because the collection happens in GiNaC, the post-processing is a handful of
lines: extract the ``N(s)/D(s)`` coefficients, then read off the DC gain
:math:`N(0)/D(0)`, and -- for the second-order denominator
:math:`D(s) = d_0 + d_1 s + d_2 s^2` -- the standard lowpass parameters
:math:`\omega_0 = \sqrt{d_0/d_2}` and :math:`Q = \sqrt{d_0 d_2}/d_1`.

Live results
````````````

The expressions below are **generated when this page is built** by running the
GiNaC analysis, and the final line verifies they equal the pure-sympy result.

.. exec-rst::

    try:
        from pycircuit.circuit import _ginac
    except ImportError:
        ## Optional extension absent: a DECISION, rendered as a note --
        ## not an accident for the loud fallback to shout about.
        _ginac = None
        print(".. note::")
        print("")
        print("   The optional GiNaC extension is not built in this")
        print("   documentation environment")
        print("   (``pycircuit/circuit/build_ginac_ext.sh``), so the live")
        print("   expressions on this page were not regenerated")
        print("   for this build.")
    if _ginac is not None:
        import sympy
        from sympy import Symbol, symbols
        from pycircuit.circuit import SubCircuit, R, C, IS, Nullor, gnd
        from pycircuit.circuit.analysis_ss import AC
        from pycircuit.circuit.toolkit import symbolic, ginac_toolkit
        from pycircuit.circuit import _ginac

        R1, R2, R3, C1, C2, i_s = symbols('R1 R2 R3 C1 C2 i_s', real=True, positive=True)
        s = Symbol('s', complex=True)

        cir = SubCircuit(toolkit=symbolic)
        cir['R1'] = R(1, 3, r=R1); cir['R2'] = R(1, 2, r=R2); cir['R3'] = R(1, gnd, r=R3)
        cir['C1'] = C(1, gnd, c=C1); cir['C2'] = C(2, 3, c=C2)
        cir['Nullor'] = Nullor(2, gnd, 3, gnd); cir['ISource'] = IS(1, gnd, iac=i_s)

        res = AC(cir, toolkit=ginac_toolkit).solve(s, complexfreq=True)
        H = res.v(3, gnd) / i_s
        num, den = _ginac.poly_coeffs(H, s)          # N(s)/D(s), GiNaC collect
        N = sum(c * s**k for k, c in enumerate(num))
        D = sum(c * s**k for k, c in enumerate(den))

        dc = sympy.simplify(num[0] / den[0])
        d0, d1, d2 = den
        omega0 = sympy.simplify(sympy.sqrt(d0 / d2))
        Q = sympy.simplify(sympy.sqrt(d0 * d2) / d1)
        poles = sympy.solve(D, s)

        def show(label, expr):
            print(".. math::")
            print("")
            print("   " + label + " = " + sympy.latex(expr))
            print("")

        show(r"H(s)", sympy.simplify(H))
        show(r"H_{\mathrm{DC}}", dc)
        show(r"\omega_0", omega0)
        show(r"Q", Q)
        show(r"s_{1,2}", sympy.Matrix(poles).T)

        H_ref = sympy.simplify(
            AC(cir, toolkit=symbolic).solve(s, complexfreq=True).v(3, gnd) / i_s)
        ok = sympy.simplify(N / D - H_ref) == 0
        print("Matches the pure-sympy (Symby) transfer function: **%s**." % ok)
