==========================================================
Compact symbolic solving: Sequence of Expressions (SoE)
==========================================================

.. note::

   Experimental / research prototype (see ``doc/ginac_fully_symbolic.md``), not
   a shipped analysis. It illustrates a representation for *fully symbolic*
   circuits, where the exploded transfer function ``N(s)/D(s)`` is exponentially
   large.

A fully symbolic transfer function is a ratio of network determinants, and a
symbolic determinant has exponentially many terms. Instead of expanding it,
solve ``A x = b`` as a **sequence of expressions**: symbolic Gaussian
elimination whose intermediate results are kept as fresh symbols. The solution
is then a small graph of assignments (``O(n^3)``; linear for a tridiagonal
ladder), directly evaluable, rather than one exploded expression.

The solver:

.. literalinclude:: ../../prototypes/soe_symbolic.py
   :pyobject: solve_soe

Operating on SoE results
========================

Because the SoE is a *shared* computation graph, results compose. Any node
voltage ``x[i]`` is a small expression over the shared assignment list, so a
transfer function is just a division of two of them::

    assigns, x = solve_soe(A, b)
    H = x[out] / x[in]          # one operation over the shared assignments

``H`` stays compact -- it references the same assignments -- and is evaluated by
resolving the assignment list numerically and then the ratio, never expanding to
the exponential ``N(s)/D(s)``.

Live example
============

The table below is **generated when this page is built**. It solves a fully
symbolic ``N``-section RC ladder with the SoE solver, forms the transfer
function ``v_out / v_in`` by dividing two solved node voltages, reports the
op-count of the compact SoE form, and checks the value against a direct numeric
solve.

.. exec-rst::

    import numpy as np
    import sympy
    from soe_symbolic import solve_soe, _ladder

    print(".. list-table:: Transfer function v_out/v_in via the SoE"
          " (fully symbolic RC ladder)")
    print("   :header-rows: 1")
    print("   :widths: 8 22 24")
    print("")
    print("   * - N")
    print("     - SoE-form op count")
    print("     - matches numeric solve")
    for N in (4, 8, 12, 16):
        A, b, s, R, C = _ladder(N)
        assigns, x = solve_soe(A, b)
        H = x[N - 1] / x[0]                       # divide two node voltages

        ops = sum(sympy.count_ops(e) for _, e in assigns) + sympy.count_ops(H)

        env = {R[i]: 100.0 * (i + 1) for i in range(N)}
        env.update({C[i]: 1e-9 * (i + 1) for i in range(N)})
        env[s] = 1j * 2 * np.pi * 1e6

        # evaluate H via the SoE (resolve assignments in order, then the ratio)
        e = dict(env)
        for t, expr in assigns:
            e[t] = complex(expr.subs(e))
        H_soe = complex(H.subs(e))

        # reference: substitute values first, then solve numerically
        An = np.array([[complex(A[i, j].subs(env)) for j in range(N)]
                       for i in range(N)])
        bn = np.array([complex(b[i].subs(env)) for i in range(N)])
        xn = np.linalg.solve(An, bn)
        ok = abs(H_soe - xn[N - 1] / xn[0]) < 1e-9 * max(1, abs(xn[N - 1] / xn[0]))

        print("   * - %d" % N)
        print("     - %d" % ops)
        print("     - %s" % ok)

The op-count grows only linearly in ``N`` while the exploded ``N(s)/D(s)`` grows
exponentially -- yet the divided transfer function evaluates to the same value.
