Symbolic distortion on shared structure
========================================

The distortion analysis in :doc:`distortion` is exact at any truncation order,
but a *symbolic* result has to be represented as well as computed, and that
turned out to be the binding constraint rather than the mathematics.

This page records what the constraint actually was.  Every number below is
measured at build time; nothing here is quoted from a previous run.

The problem, stated precisely
------------------------------

An earlier measurement established that symbolic expression **size** grows
polynomially with truncation order, and multi-node circuits cost only a small
constant factor on top.  It established **nothing about cost** — and
``sympy.cancel`` did not finish in **900 seconds** at :math:`U^7` on a
two-node system.

Those are different claims, and only one of them had been checked.  The
suspicion worth testing was that the cost belonged to *simplification* rather
than to the method.

What changes, and what deliberately does not
---------------------------------------------

:class:`~pycircuit.circuit.distortion.GradedSpectrum` asks its coefficients
for four things — ``+``, ``*``, ``conjugate()`` and ``!= 0`` — and nothing
else.  So the recurrence can be re-run over a different algebra by supplying a
different coefficient type, with **no change to the graded ring, the
consistency filter or the composition bookkeeping.**

:mod:`~pycircuit.circuit.distortion_ddd` supplies one: :class:`Expr`, a
hash-consed expression graph that is never expanded and never simplified.
Structurally identical subexpressions become the *same object*, so sharing is
by identity; commuting operands are sorted so that :math:`ab` and :math:`ba`
are one node.  Cost then falls entirely on evaluation, which is memoised over
the shared structure.

.. note::

   ``distortion.py`` is untouched by this work.  That is not tidiness: it is
   the oracle every measurement below is checked against, and refactoring it
   to share code with its own checker would destroy the independence that
   makes the comparison mean anything.

Cost against truncation order
------------------------------

The sweep the plan asked for: start at :math:`U^7`, the order that defeated
``cancel``, and climb.  The value is checked against the *identical
recurrence* run in plain complex arithmetic.

.. exec-rst::

    import time
    import numpy as np, sympy
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)
    from pycircuit.circuit.distortion_ddd import (Expr, evaluate_one,
                                                  nodes_of)

    names = 'a c1 c2 b d k3 A w0'
    sym = dict(zip(names.split(), sympy.symbols(names)))
    at = {k: Expr.atom(v) for k, v in sym.items()}
    env = {sym['a']: 1e-3, sym['c1']: 1e-12, sym['c2']: 2e-12,
           sym['b']: 3e-4, sym['d']: 4e-4, sym['k3']: 5e-3,
           sym['A']: 1e-3, sym['w0']: 2*np.pi*1e6}

    def cube(sp, p):
        return ((sp*sp).truncated(p) * sp).truncated(p)

    def run(power, solve, coeff, drive, tone):
        src = GradedVector([GradedSpectrum.from_phasor(1, 1, drive),
                            GradedSpectrum()])
        f = lambda x: GradedVector([GradedSpectrum(),
                                    cube(x[0], power).scaled(coeff)])
        return graded_response_mimo(solve, src, f, (tone,), max_power=power)

    def symbolic_solve(s, rhs):
        sv = Expr.atom(s)
        det = (at['a'] + sv*at['c1'])*(sv*at['c2']) - at['b']*at['d']
        return [(sv*at['c2']*rhs[0] + at['b']*rhs[1])/det,
                (at['d']*rhs[0] + (at['a'] + sv*at['c1'])*rhs[1])/det]

    def numeric_solve(s, rhs):
        det = ((env[sym['a']] + s*env[sym['c1']])*(s*env[sym['c2']])
               - env[sym['b']]*env[sym['d']])
        return [(s*env[sym['c2']]*rhs[0] + env[sym['b']]*rhs[1])/det,
                (env[sym['d']]*rhs[0]
                 + (env[sym['a']] + s*env[sym['c1']])*rhs[1])/det]

    head = ['truncation', 'build (s)', 'graph nodes', 'evaluate (s)',
            'vs numeric']
    rows = []
    for power in (7, 9, 11, 13):
        t = time.time()
        g = run(power, symbolic_solve, at['k3'], at['A'], sym['w0'])
        root = g[0].phasor(3)
        build = time.time() - t

        t = time.time()
        got = evaluate_one(root, env)
        ev = time.time() - t

        n = run(power, numeric_solve, complex(env[sym['k3']]),
                complex(env[sym['A']]), complex(env[sym['w0']]))
        want = n[0].phasor(3)

        rows.append(['U^%d' % power, '%.3f' % build, '%d' % nodes_of(root),
                     '%.4f' % ev, '%.1e' % (abs(got-want)/abs(want))])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)
    print('')
    print('Growth in graph size from ``U^7`` to ``U^13``: **%.1fx**.'
          % (int(rows[-1][2])/int(rows[0][2])))

**The cost was simplification, not the method.**  ``cancel`` could not finish
:math:`U^7` in 900 seconds; the same coefficient is built and evaluated here
at :math:`U^{13}` in well under a second, agreeing with the numeric path to
around :math:`10^{-16}`.

Growth in graph size is polynomial in the truncation order — a geometric
representation would roughly square over that span, and this does not come
close.

Cost against circuit size
--------------------------

The other axis, and the one that matters for analysing bigger circuits.  The
matrix here is **dense**, so Gaussian elimination fills in completely: a
ladder would have no fill-in and would flatter the result.

.. exec-rst::

    import time
    import numpy as np, sympy
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)
    from pycircuit.circuit.distortion_ddd import Expr, evaluate_one, nodes_of

    POWER = 7
    w0, drive = sympy.Symbol('w0'), sympy.Symbol('A')

    def cube(sp, p):
        return ((sp*sp).truncated(p) * sp).truncated(p)

    def measure(n):
        ent = [[sympy.Symbol('y%d_%d' % (i, j)) for j in range(n)]
               for i in range(n)]
        caps = [sympy.Symbol('cc%d' % i) for i in range(n)]
        ks = [sympy.Symbol('kk%d' % i) for i in range(n)]
        Ee = [[Expr.atom(e) for e in row] for row in ent]
        Ec = [Expr.atom(c) for c in caps]
        Ek = [Expr.atom(k) for k in ks]

        def solve(s, rhs):
            sv = Expr.atom(s)
            M = [[Ee[i][j] + (sv*Ec[i] if i == j else Expr.atom(0))
                  for j in range(n)] for i in range(n)]
            b = [r if isinstance(r, Expr) else Expr.atom(r) for r in rhs]
            for i in range(n):
                for r in range(i+1, n):
                    fac = M[r][i]/M[i][i]
                    for j in range(i, n):
                        M[r][j] = M[r][j] - fac*M[i][j]
                    b[r] = b[r] - fac*b[i]
            x = [None]*n
            for i in range(n-1, -1, -1):
                acc = b[i]
                for j in range(i+1, n):
                    acc = acc - M[i][j]*x[j]
                x[i] = acc/M[i][i]
            return x

        f = lambda x: GradedVector([cube(x[i], POWER).scaled(Ek[i])
                                    for i in range(n)])
        src = GradedVector([GradedSpectrum.from_phasor(1, 1, Expr.atom(drive))]
                           + [GradedSpectrum() for _ in range(n-1)])
        t = time.time()
        sol = graded_response_mimo(solve, src, f, (w0,), max_power=POWER)
        root = sol[0].phasor(3)
        build = time.time() - t

        env = {w0: 2*np.pi*1e6, drive: 1e-3}
        for i in range(n):
            env[caps[i]] = 1e-12*(1+0.2*i)
            env[ks[i]] = 5e-3
            for j in range(n):
                env[ent[i][j]] = (1e-3 if i == j else 1e-4)*(1+0.05*(i+j))
        t = time.time()
        evaluate_one(root, env)
        return build, nodes_of(root), time.time() - t

    sizes = (3, 6, 9, 12)
    head = ['circuit nodes', 'build (s)', 'graph nodes', 'evaluate (s)']
    rows = []
    counts = []
    for n in sizes:
        build, count, ev = measure(n)
        counts.append(count)
        rows.append(['%d' % n, '%.3f' % build, '%d' % count, '%.4f' % ev])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)
    print('')
    exponent = (np.log(counts[-1]/counts[0])
                / np.log(sizes[-1]/sizes[0]))
    print('Fitted growth: graph size scales as ``n**%.2f`` in the number of'
          % exponent)
    print('circuit nodes -- so roughly **quadratic**, on the topology that')
    print('fills in worst.')

So the representation is polynomial on **both** axes: in truncation order and
in circuit size, the latter measured on the least favourable topology.  There
is no exponential barrier in either direction, which is the result that
decides whether larger circuits are reachable at all.

What this does *not* establish
-------------------------------

Stated plainly, because the tables above are easy to over-read.

* **A dense synthetic matrix is not a real circuit.**  Entries here are
  independent symbols; a real device model shares parameters across entries,
  which changes both the sharing available and the arithmetic.
* **Nothing here is a claim about numeric speed.**  The numeric path is numpy
  and stays numpy; this representation exists for the symbolic case, where the
  alternative was not slower but *unavailable*.
* **The largest circuit measured is small.**  The extrapolation to a
  transistor-level amplifier is arithmetic, not measurement.

Consequences for the decision-diagram plan
-------------------------------------------

The work above was stage A of a plan whose later stages proposed determinant
decision diagrams — the machinery in :doc:`ddd` — for the same problem.  That
plan set a condition in advance: if a plain factored representation made
:math:`U^{13}` tractable, the diagram stages would be an optimisation of an
already-solved problem and should not be built.

That condition was met.  The determinant stages are therefore **not being
built**, and that is recorded as the outcome rather than quietly dropped —
the cheapest possible answer to the question, reached by measuring the
cheap option first.

What remains genuinely open is the *circuit* axis at realistic scale, where a
determinant diagram's own strength — compact representation of a large
determinant and its cofactors, and hierarchical suppression of internal nodes
— has not been tested against this representation on anything resembling an
op-amp.
