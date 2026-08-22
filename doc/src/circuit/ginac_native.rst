=====================================================
Native GiNaC result object (fast symbolic, no bridge)
=====================================================

.. note::

   Experimental / research backend (see ``doc/ginac_fully_symbolic.md``).
   Requires the optional compiled ``_ginac_ext`` extension
   (``pycircuit/circuit/build_ginac_ext.sh``); when it is not built the live
   tables below are replaced by a note saying so.

Solving a symbolic circuit in GiNaC returns the numerators and shared
denominator as a ratio of network determinants -- expressions that are
determinant-sized for a large circuit. Converting that *whole* result back to
sympy dominates the end-to-end time, yet callers almost always want only a
**compact** piece: a transfer function (the shared denominator cancels) or the
denominator polynomial (for poles).

:func:`pycircuit.circuit._ginac.solve_native` keeps ``num[]``/``den`` as native
GiNaC ``ex`` inside a light :class:`~pycircuit.circuit._ginac.GinacResult`
wrapper and does the cancellation *in GiNaC*, so only the small result crosses
back to sympy::

    from pycircuit.circuit.toolkit import ginac_toolkit

    res = ginac_toolkit.solve_native(A, b)
    H   = res.tf(out, inp)            # x_out/x_in, cancelled natively -> compact
    num, den = res.tf_coeffs(out, inp, s)   # canonical N(s)/D(s) coefficients
    poly = res.denominator_coeffs(s)        # characteristic polynomial (poles)

Extracting a transfer function
==============================

The table below is **generated when this page is built**. It solves a fully
symbolic ``N``-section RC ladder and extracts the transfer function
``v_{N-1}/v_0`` two ways -- the eager path (sympify the full result, then
``cancel`` in sympy) and the native path (``solve_native`` + ``tf``, cancelled
in GiNaC) -- timing each and checking the two agree.

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
        print("   tables on this page were not regenerated")
        print("   for this build.")
    if _ginac is not None:
        import time
        import numpy as np
        import sympy
        from pycircuit.circuit import _ginac

        def rc_ladder(n):
            R = sympy.symbols('R1:%d' % (n + 1), positive=True)
            C = sympy.symbols('C1:%d' % (n + 1), positive=True)
            s = sympy.Symbol('s', imaginary=True)
            A = sympy.zeros(n, n)
            for k in range(n):
                g = 1 / R[k]
                gnext = 1 / R[k + 1] if k + 1 < n else 0
                A[k, k] = g + gnext + s * C[k]
                if k + 1 < n:
                    A[k, k + 1] = -gnext
                    A[k + 1, k] = -gnext
            b = np.zeros(n, dtype=object)
            b[0] = 1 / R[0]
            return np.array(A.tolist(), dtype=object), b

        print(".. list-table:: Transfer function v_{N-1}/v_0 -- eager vs native")
        print("   :header-rows: 1")
        print("   :widths: 6 16 16 10 12")
        print("")
        print("   * - N")
        print("     - eager (s)")
        print("     - native (s)")
        print("     - speed-up")
        print("     - agree")
        for N in (4, 5, 6):
            A, b = rc_ladder(N)

            t = time.time()
            num, den = _ginac.linearsolver_num_den(A, b)
            tf_eager = sympy.cancel(num[N - 1] / num[0])
            t_eager = time.time() - t

            t = time.time()
            res = _ginac.solve_native(A, b)
            tf_native = res.tf(N - 1, 0)
            t_native = time.time() - t

            ok = sympy.simplify(tf_eager - tf_native) == 0
            print("   * - %d" % N)
            print("     - %.3f" % t_eager)
            print("     - %.3f" % t_native)
            print("     - %.1fx" % (t_eager / t_native))
            print("     - %s" % ok)

The native path avoids the determinant-sized sympify round-trip; the advantage
widens with ``N`` (at ``N=8`` it is ~2.8x, 10.5 s to 3.8 s).

Coefficient extraction (canonical N/D)
======================================

``tf_coeffs`` collects the transfer function into polynomial-in-``s``
coefficient vectors -- the canonical ``N(s)/D(s)`` form for poles/zeros or a
``scipy.signal`` hand-off. The collection runs natively in GiNaC. The table
compares it against extracting the same coefficients with sympy ``Poly`` (both
starting from the native ``tf``), and reconstructs ``N/D`` to confirm agreement.

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
        print("   tables on this page were not regenerated")
        print("   for this build.")
    if _ginac is not None:
        import time
        import numpy as np
        import sympy
        from pycircuit.circuit import _ginac

        def rc_ladder(n):
            R = sympy.symbols('R1:%d' % (n + 1), positive=True)
            C = sympy.symbols('C1:%d' % (n + 1), positive=True)
            s = sympy.Symbol('s', imaginary=True)
            A = sympy.zeros(n, n)
            for k in range(n):
                g = 1 / R[k]
                gnext = 1 / R[k + 1] if k + 1 < n else 0
                A[k, k] = g + gnext + s * C[k]
                if k + 1 < n:
                    A[k, k + 1] = -gnext
                    A[k + 1, k] = -gnext
            b = np.zeros(n, dtype=object)
            b[0] = 1 / R[0]
            return np.array(A.tolist(), dtype=object), b, s

        print(".. list-table:: N(s)/D(s) coefficients -- native tf_coeffs vs sympy Poly")
        print("   :header-rows: 1")
        print("   :widths: 6 16 16 10 12")
        print("")
        print("   * - N")
        print("     - sympy Poly (s)")
        print("     - native (s)")
        print("     - speed-up")
        print("     - reconstructs")
        for N in (4, 5, 6):
            A, b, s = rc_ladder(N)
            res = _ginac.solve_native(A, b)
            tf = res.tf(N - 1, 0)

            t = time.time()
            num, den = res.tf_coeffs(N - 1, 0, s)
            t_native = time.time() - t

            t = time.time()
            tfc = sympy.cancel(tf)
            Pn = sympy.Poly(sympy.numer(tfc), s).all_coeffs()
            Pd = sympy.Poly(sympy.denom(tfc), s).all_coeffs()
            t_sympy = time.time() - t

            Nn = sum(c * s**k for k, c in enumerate(num))
            Dd = sum(c * s**k for k, c in enumerate(den))
            ok = sympy.simplify(Nn / Dd - tf) == 0
            print("   * - %d" % N)
            print("     - %.4f" % t_sympy)
            print("     - %.4f" % t_native)
            print("     - %.1fx" % (t_sympy / t_native))
            print("     - %s" % ok)

Both operate on the same native ``ex`` and return only the compact structured
result, so the "manipulate in GiNaC, convert the small piece" win compounds as
the circuit grows.
