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

A real circuit
---------------

The measurements above use synthetic matrices, which is the right way to
isolate scaling but says nothing about whether a real circuit behaves the
same.  The µA741 below is transcribed from Tan & Shi's TCAD 2000 paper — the
circuit their published diagram sizes are measured on — and its matrix is
built through the *same* path :class:`~pycircuit.circuit.analysis_ss.AC` uses,
so it is the analysis's own matrix rather than a hand-written one.

It is 26x26 with 103 nonzeros: sparse, and mostly numeric with a handful of
symbolic device parameters, which is what a real analysis looks like and is
quite unlike a dense matrix of distinct symbols.

.. exec-rst::

    import time
    import numpy as np, sympy
    from pycircuit.circuit import benchmark_circuits as bc
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)
    from pycircuit.circuit.distortion_ddd import Expr, evaluate_one, nodes_of

    def cube(sp, p):
        return ((sp*sp).truncated(p) * sp).truncated(p)

    def expr_solve(A, s_sym, n):
        rows = [[A[i, j] for j in range(n)] for i in range(n)]
        def solve(s, rhs):
            M = [[Expr.atom(e.subs(s_sym, s)
                            if getattr(e, 'free_symbols', set()) else e)
                  if e != 0 else Expr.atom(0) for e in row] for row in rows]
            b = [r if isinstance(r, Expr) else Expr.atom(r) for r in rhs]
            for i in range(n):
                for r in range(i+1, n):
                    if M[r][i]._is_zero():
                        continue
                    fac = M[r][i]/M[i][i]
                    for j in range(i, n):
                        if not M[i][j]._is_zero():
                            M[r][j] = M[r][j] - fac*M[i][j]
                    b[r] = b[r] - fac*b[i]
            x = [None]*n
            for i in range(n-1, -1, -1):
                acc = b[i]
                for j in range(i+1, n):
                    if not M[i][j]._is_zero():
                        acc = acc - M[i][j]*x[j]
                x[i] = acc/M[i][i]
            return x
        return solve

    def num_solve(A, s_sym, n, env):
        rows = sympy.Matrix(A)
        def solve(s, rhs):
            loc = dict(env); loc[s_sym] = s
            M = np.array([[complex(rows[i, j].subs(loc)) for j in range(n)]
                          for i in range(n)], dtype=complex)
            return list(np.linalg.solve(M, np.asarray(rhs, dtype=complex)))
        return solve

    sysm = bc.ua741(symbolic_devices=('q1', 'q2', 'q3', 'q4'))
    n, s_sym = sysm.dim, sysm.s
    IN = [i for i in range(n) if sysm.b[i] != 0][0]
    OUT = sysm.out_index
    NL = [17, 18, 22]
    ks = [sympy.Symbol('kk%d' % i) for i in NL]
    drive = sympy.Symbol('Adrv')
    env = {k: 1e-9 for k in ks}
    env[drive] = 1e-3
    for d in sysm.A.free_symbols:
        if str(d).startswith('gm_'):
            env[d] = 1e-3

    def make_f(coeffs, power):
        def f(x):
            out = [GradedSpectrum() for _ in range(n)]
            for idx, node in enumerate(NL):
                out[node] = cube(x[node], power).scaled(coeffs[idx])
            return GradedVector(out)
        return f

    def source(value):
        vec = GradedVector([GradedSpectrum() for _ in range(n)])
        vec.components[IN].terms.update(
            GradedSpectrum.from_phasor(1, 1, value).terms)
        return vec

    tone = 2*np.pi*1e3
    head = ['truncation', 'build (s)', 'graph nodes', 'evaluate (s)',
            'vs numeric']
    rows = []
    for power in (3, 5, 7):
        t = time.time()
        g = graded_response_mimo(expr_solve(sysm.A, s_sym, n),
                                 source(Expr.atom(drive)),
                                 make_f([Expr.atom(k) for k in ks], power),
                                 (tone,), max_power=power)
        root = g[OUT].phasor(3)
        build = time.time() - t
        t = time.time(); got = evaluate_one(root, env); ev = time.time() - t

        nn = graded_response_mimo(num_solve(sysm.A, s_sym, n, env),
                                  source(complex(env[drive])),
                                  make_f([complex(env[k]) for k in ks], power),
                                  (tone,), max_power=power)
        want = nn[OUT].phasor(3)
        rows.append(['U^%d' % power, '%.2f' % build, '%d' % nodes_of(root),
                     '%.4f' % ev, '%.1e' % (abs(got-want)/abs(want))])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)

**A transistor-level operational amplifier is comfortably in reach**, with the
symbolic result agreeing with the numeric one to around :math:`10^{-14}`.

Two things about that table are worth knowing rather than inferring.

The graph is barely sensitive to how many device parameters are left
symbolic — going from two symbolic transconductances to eight changed the node
count by under 1%, because symbols are *leaves* on structure that is shared
regardless.

And on cascaded amplifiers the growth is **linear** in circuit size when the
nonlinearity sits far from the output.  Placing it *near* the output instead
makes the size constant, which looks better and measures the wrong thing: it
reports the locality of the nonlinear subcircuit rather than the scale of the
circuit.  The linear figure is the honest one.

.. warning::

   A size measurement says nothing about whether the analysis did anything.
   An earlier version of this page's measurement drove node 0 and read node 0,
   producing entirely plausible build times and node counts against a third
   harmonic that was **identically zero**.  It was caught only by comparing
   against the numeric path, which is why every table here carries that
   column.

Consequences for hierarchy
---------------------------

The plan behind this work proposed hierarchical node suppression — the
machinery in :doc:`ddd`, which reduces the µA741's determinant diagram from
1040 vertices to 156 — as the route to larger circuits.

It is not needed here.  Suppression exists to rescue a representation that is
blowing up, and this one is linear in circuit size and builds a transistor-level
amplifier in under a second.  That stage is therefore **not being built**, and
this paragraph records why rather than leaving the plan's proposal dangling.

What kind of object this is — and what it is not
-------------------------------------------------

.. warning::

   **This page originally invited a comparison it should not have.**  The
   graph sizes above are much smaller than the diagram sizes in :doc:`ddd`,
   and it is natural to read that as the representation being better.  It is
   not: **the two represent different objects**, and only one of them is what
   symbolic analysis usually means.

Measured on the same circuit — the µA741 with all 24 transconductances and
:math:`s` symbolic:

.. exec-rst::

    import time
    from pycircuit.circuit import benchmark_circuits as bc
    from pycircuit.circuit.ddd import ddd_of_matrix
    from pycircuit.circuit.distortion_ddd import Expr, nodes_of

    sysm = bc.ua741(symbolic_devices=tuple('q%d' % i for i in range(1, 25)))
    n = sysm.dim

    D = ddd_of_matrix(sysm.A, order='min-degree')

    rows = [[sysm.A[i, j] for j in range(n)] for i in range(n)]
    M = [[Expr.atom(e) if e != 0 else Expr.atom(0) for e in row]
         for row in rows]
    det = Expr.atom(1)
    for i in range(n):
        det = det * M[i][i]
        for r in range(i+1, n):
            if M[r][i]._is_zero():
                continue
            fac = M[r][i]/M[i][i]
            for j in range(i, n):
                if not M[i][j]._is_zero():
                    M[r][j] = M[r][j] - fac*M[i][j]

    print('The **same determinant**, two representations:')
    print('')
    print('* a determinant decision diagram -- ``%d`` vertices, standing for'
          % D.size)
    print('  ``%d`` expanded product terms;' % D.term_count())
    print('* this page\'s expression graph -- ``%d`` nodes.' % nodes_of(det))

The graph is smaller because it is **a record of the elimination's
arithmetic** — a straight-line program of :math:`O(n^3)` operations — and it
never expands, so the size of the expanded polynomial never costs it
anything.  The diagram is a *canonical form of that expanded polynomial*.

The consequence is a capability difference, not a size one:

.. code-block:: text

    expression graph        decision diagram
    ----------------        ----------------
    evaluate fast, at       evaluate
    many parameter values
                            coefficients of s**k  (poles, zeros)
                            rank terms by magnitude
                            term-ranked approximation
                            sensitivity to any parameter

**So the honest summary of this page is narrower than it first appears.**  The
cost recorded in :doc:`distortion` section 8.3 was removed by *not producing
the expensive object*: no normal form is ever built, so nothing has to be
simplified.  That is exactly the right trade when the goal is a **number** —
distortion at given parameter values, swept over frequency or over a
parameter, which is what this analysis is usually for.  It is the wrong trade
when the goal is an **expression** a designer reads, and for that the
diagrams in :doc:`ddd` remain the tool.

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

The circuit axis, at op-amp scale
---------------------------------

The paragraph above records the circuit axis as open: this representation had never
been tried on anything resembling an op-amp.  It has now, on the 5th-order leapfrog
filter of :func:`~pycircuit.circuit.benchmark_circuits.leapfrog_5th_order` — five
µA741s, **127 unknowns** — with a cubic ``i = kk*v**3`` on the first amplifier's input
differential pair.

Two things make the measurement worth having. The **drive amplitude and the
nonlinearity strength both stay symbolic**, so one build per order serves an entire
sweep: the expensive step happens once and every amplitude is a walk over the same
graph. And the question asked is *self-convergence* — at which truncation order does
the third harmonic stop moving — which needs no external oracle.

.. exec-rst::

    import numpy as np
    from order_convergence import (AMPLITUDES, F0, G3_VALUE, NODE_NAME,
                                   build_symbolic, turning_point)
    from pycircuit.circuit import benchmark_circuits as bc
    from pycircuit.circuit.benchmark_circuits import _UA741_NODE_NAMES
    from pycircuit.circuit.distortion_ddd import evaluate_one

    ## Three well-separated orders rather than every one: the build cost is
    ## cumulative and this page's timeout has to clear the worst case, not the
    ## typical one.  The full sweep to U^17 is tabulated below.
    orders = (3, 7, 11)

    system = bc.leapfrog_5th_order()
    names = ['in'] + ['s%d_%s' % (k, b) for k in range(5)
                      for b in _UA741_NODE_NAMES]
    node_row = names.index(NODE_NAME)
    drive_row = [i for i in range(system.dim) if system.b[i] != 0][0]
    base = {s: complex(system.params[s]) for s in system.A.free_symbols
            if s is not system.s}
    vturn, gnode = turning_point(system, node_row, G3_VALUE, F0)

    built = {}
    for order in orders:
        b = build_symbolic(system, node_row, system.out_index, drive_row, order)
        b['vals'] = {}
        for amp in AMPLITUDES:
            env = dict(base)
            env[b['drive']] = complex(amp)
            env[b['kk']] = complex(G3_VALUE)
            b['vals'][amp] = (evaluate_one(b['fund'], env),
                              evaluate_one(b['third'], env),
                              evaluate_one(b['nodev'], env))
        built[order] = b

    print(".. list-table:: Leapfrog, %d unknowns: third harmonic by truncation order"
          % system.dim)
    print("   :header-rows: 1")
    print("   :widths: 12 14 10 14 14 14")
    print("")
    print("   * - drive (V)")
    ## NOT "|v|": reST reads |...| as a substitution reference and errors, and it
    ## reports the error at the top of the page rather than at the block.
    print("     - node volts")
    ## A plain string, not a %-format, so "%%" would render literally as "%%" --
    ## which is exactly what the built page showed until 2026-07-30.  Neighbouring
    ## lines here ARE %-formatted, which is how the escape got copied in.
    print("     - percent of v_turn")
    print("     - HD3")
    for prev, order in zip(orders, orders[1:]):
        print("     - change at U^%d" % order)
    for amp in AMPLITUDES:
        top = built[orders[-1]]['vals'][amp]
        hd3 = abs(top[1]) / abs(top[0])
        print("   * - %g" % amp)
        print("     - %.3e" % abs(top[2]))
        print("     - %.0f%%" % (100.0 * abs(top[2]) / vturn))
        print("     - %.3e" % hd3)
        for prev, order in zip(orders, orders[1:]):
            a = built[prev]['vals'][amp][1]
            c = built[order]['vals'][amp][1]
            print("     - %.1e" % (abs(c - a) / abs(c) if abs(c) else float('nan')))
    print("")
    print("Built live: %d graph nodes at U^%d, and the cubic's validity limit is"
          " v_turn = %.3e V (g = %.3e S at %s)."
          % (built[orders[-1]]['nodes'], orders[-1], vturn, gnode, NODE_NAME))

**The order needed tracks how close the node is to the cubic's turning point**, not the
drive level as such.  That limit matters and is easy to forget: a truncated cubic
``i = g(v + a*v**3)`` has negative differential conductance beyond
:math:`v_{\mathrm{turn}} = 1/\sqrt{3|a|}` and **is not a physical device** there, so an
order-convergence study only means something inside it.  An earlier version of this
table reported a non-convergence at 366% of :math:`v_{\mathrm{turn}}` as if it were a
property of the series; it was a property of a model that had stopped being a circuit.

Frequency matters as much as amplitude, for a reason worth stating: at 1 kHz the
amplifier's loop gain holds that input pair at a virtual ground — 2.5 µV for a 10 mV
drive — and every amplitude converges at :math:`U^3` with HD3 near
:math:`10^{-14}`.  Distortion rises where loop gain falls.  The table above is taken at
100 kHz, where this µA741's open-loop gain is down to roughly 8 dB and the same drive
gives 6.3e-05 at that node, **25× larger**.

Extending to :math:`U^{17}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The numbers below are **not** built on this page.  The full sweep costs about 330 s of
symbolic builds, and this page's ``exec_rst_timeout`` is 300 s — a block that overruns
renders its own source instead of a table, so the live example above stops at
:math:`U^{11}`.  Reproduce these with
``PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/order_convergence.py``.

.. list-table:: Build cost per order (127 unknowns, one build serves every amplitude)
   :header-rows: 1
   :widths: 12 16 14 22

   * - order
     - graph nodes
     - build
     - evaluate, 6 amplitudes
   * - :math:`U^3`
     - 6 415
     - 1.7 s
     - 0.13 s
   * - :math:`U^9`
     - 30 075
     - 20.1 s
     - 1.09 s
   * - :math:`U^{13}`
     - 54 079
     - 53.0 s
     - 2.09 s
   * - :math:`U^{15}`
     - 70 495
     - 79.7 s
     - 2.77 s
   * - :math:`U^{17}`
     - 84 427
     - 123.5 s
     - 3.50 s

Successive corrections to the third harmonic, at 12% and 34% of
:math:`v_{\mathrm{turn}}`:

.. list-table:: Convergence to :math:`U^{17}`
   :header-rows: 1
   :widths: 14 12 12 12 12 12 12 12

   * - drive
     - :math:`U^5`
     - :math:`U^7`
     - :math:`U^9`
     - :math:`U^{11}`
     - :math:`U^{13}`
     - :math:`U^{15}`
     - :math:`U^{17}`
   * - 1 V (12%)
     - 1.6e-02
     - 2.9e-04
     - 5.8e-06
     - 1.2e-07
     - 2.6e-09
     - 5.8e-11
     - 1.3e-12
   * - 3 V (34%)
     - 1.7e-01
     - 2.7e-02
     - 4.7e-03
     - 8.8e-04
     - 1.7e-04
     - 3.4e-05
     - 7.1e-06

**Both rows are converging geometrically; the 3 V row simply needs more orders than
17.**  Its successive ratios are 0.16, 0.17, 0.19, 0.19, 0.20, 0.21 — a convergent
series creeping toward its radius, not a divergent one.  At 12% the ratio is a steady
0.022.

That gives a usable rule.  The per-two-order ratio scales close to the square of the
distance to the turning point: :math:`0.022` at 12% and :math:`0.20` at 34% are
:math:`1.5(v/v_{\mathrm{turn}})^2` and :math:`1.7(v/v_{\mathrm{turn}})^2`.  So the
order needed for a target accuracy can be estimated from the node voltage alone, before
any expansion is built — and at a third of the turning point, six-figure accuracy in the
third harmonic wants roughly :math:`U^{25}`.
