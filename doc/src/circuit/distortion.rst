Distortion analysis by perturbation
===================================

A weakly nonlinear circuit is one designed to be linear and failing slightly.
Its response to a sinusoid is *almost* a sinusoid: the deviation shows up as
harmonics, and their size relative to the fundamental is what a designer
calls distortion.

This page describes an analysis that predicts those harmonics in closed form,
by treating the circuit as a **perturbation of its own linearised model**.
The method is Buonomo and Lo Schiavo's; the primary reference is *Predicting
nonlinear distortion in multistage amplifiers and gm-C filters*, Analog
Integrated Circuits and Signal Processing, vol. 77, pp. 483–493 (2013), and
the relation to Volterra series proved in IEEE Trans. Circuits and Systems I,
vol. 52, no. 8, pp. 1620–1631 (2005).

.. note::

   This covers stages 1–4 of the plan in ``doc/distortion_plan.md``: any
   number of nonlinear devices; single-tone harmonic distortion for cubic or
   exponential nonlinearities, and two-tone intermodulation for cubic ones.

The idea
--------

Write the circuit in modified nodal form, split into its linear and strictly
nonlinear parts, and solve for the node voltages :math:`x`:

.. math::

   Y(s)\, x = u - f(x), \qquad Y(s) = G + sC

where :math:`u` is the source vector and :math:`f` collects the nonlinear
device currents.  :math:`f` is *strictly* nonlinear — any linear part of a
device has already been absorbed into :math:`G`, so :math:`f(0) = 0` and
:math:`f'(0) = 0`.

Inverting the linear part gives the form the method actually uses:

.. math::

   x = \mathcal{G}u - \mathcal{G}f(x), \qquad \mathcal{G} = Y(s)^{-1}

Now introduce a bookkeeping parameter :math:`\varepsilon` on the nonlinear
term and expand the solution in powers of it:

.. math::

   x = \mathcal{G}u - \varepsilon\,\mathcal{G}f(x), \qquad
   x = x_0 + \varepsilon x_1 + \varepsilon^2 x_2 + \dots

Matching powers of :math:`\varepsilon` gives the recurrence

.. math::

   x_0 &= \mathcal{G}u \\
   x_1 &= -\mathcal{G}f(x_0) \\
   x_2 &= -\mathcal{G}\,\frac{\partial f}{\partial x}\bigg|_{x_0} x_1

and setting :math:`\varepsilon = 1` at the end recovers the original problem.

Two things about :math:`\varepsilon` are worth being explicit about, because
they are easy to misread.  It is **fictitious** — it has no physical meaning
and no dimension, and no circuit quantity is equal to it.  It exists only to
order the terms.  What actually has to be small is :math:`\mathcal{G}f(x)`
compared with :math:`\mathcal{G}u`, and **none of the papers in this line
gives a quantitative bound on how small.**  That gap is real and is discussed
under `Validity`_ below.

Why every step is a linear solve
--------------------------------

Look at what the recurrence asks for.  :math:`x_0` is the ordinary linear
solution.  :math:`x_1` evaluates the nonlinearity **on a solution that is
already known** and pushes the result back through the same
:math:`\mathcal{G}`.  So does :math:`x_2`.  At no point is there a nonlinear
system to solve — no Newton iteration, no convergence criterion, no initial
guess.

That is the property the implementation is built around, and it is why this
analysis needs no changes to the circuit or element layers: :math:`Y(s)` is
exactly the matrix AC analysis already assembles.

Where harmonic balance enters
-----------------------------

The recurrence above is in the time domain, and :math:`f(x_0)` is a periodic
waveform, not a phasor.  The step that makes the method practical is to
expand it in a Fourier series — and this is the only place harmonic balance
appears.

It is worth being precise here, because the phrase "perturbation plus
harmonic balance" suggests the two methods alternate.  They do not.
**Harmonic balance closes each perturbation order; it is not used to find the
fundamental.**  The fundamental comes from an ordinary linear AC solve.

Concretely, for a single tone :math:`u(t) = \mathrm{Re}[U e^{j\omega_0 t}]`,
the linear solution is :math:`X_1 = \mathcal{G}(j\omega_0)U`.  Squaring it,

.. math::

   \left(\mathrm{Re}\!\left[X_1 e^{j\theta}\right]\right)^2
     = \mathrm{Re}\!\left[\tfrac{1}{2}X_1^2 e^{2j\theta}\right]
       + \tfrac{1}{2}|X_1|^2

so a quadratic nonlinearity :math:`b\,v^2` produces a second harmonic of
amplitude :math:`\tfrac{1}{2}b X_1^2` and a DC shift.  Each of those
components is then passed back through :math:`\mathcal{G}` **evaluated at its
own frequency** — the second harmonic through
:math:`\mathcal{G}(j2\omega_0)`, not :math:`\mathcal{G}(j\omega_0)`.  That
per-harmonic re-evaluation is the whole content of the step.

For a cubic nonlinearity :math:`f(v) = b v^2 + c v^3` the result to second
order is

.. math::

   \hat{X}_1 &= \mathcal{G}(j\omega_0)\,U \\
   \hat{X}_2 &= -\tfrac{1}{2}\,\mathcal{G}(j2\omega_0)\, b\,\hat{X}_1^2 \\
   \hat{X}_3 &= -\tfrac{1}{4}\,\mathcal{G}(j3\omega_0)
                \left[c\,\hat{X}_1^3
                      - 2b^2\,\hat{X}_1^3\,\mathcal{G}(j2\omega_0)\right]

and the distortion figures are :math:`\mathrm{HD}_2 = |\hat{X}_2|/|\hat{X}_1|`
and :math:`\mathrm{HD}_3 = |\hat{X}_3|/|\hat{X}_1|`.

The second term in :math:`\hat{X}_3` deserves attention.  It contains no cubic
coefficient at all — it is the *quadratic* nonlinearity acting twice, once to
make a second harmonic and once to mix that back up against the fundamental.
A circuit whose dominant nonlinearity is quadratic can therefore have its
third-harmonic distortion dominated by :math:`b`, and an analysis that stops
at first perturbation order will miss it entirely.

What the truncation is actually truncating
------------------------------------------

Everything in this method is sized by powers of the drive amplitude
:math:`U`:

.. code-block:: text

    fundamental   ~ U
    2nd harmonic  ~ U^2
    3rd harmonic  ~ U^3

So computing the third harmonic means finding **every route that produces
something of size** :math:`U^3`, and ignoring anything smaller.  There are
exactly two:

.. code-block:: text

    ROUTE 1 -- one big step
        the cubic coefficient acts on the fundamental
            c * (U)^3  =  U^3                     <- first perturbation order

    ROUTE 2 -- two small steps
        the quadratic coefficient makes a 2nd harmonic
            b * (U)^2  =  U^2
        that 2nd harmonic then beats against the fundamental
            b * U^2 * U  =  U^3                   <- second perturbation order

Both arrive at :math:`U^3`, so the method keeps both.  Everything else is
smaller and is dropped — the second harmonic beating against itself gives
:math:`U^4`, a DC shift against the third harmonic gives :math:`U^5`, and so
on.

Here is the part that catches people out.  **"Perturbation order" and "power
of** :math:`U` **" are not the same thing.**  Route 2 is *second*-order
perturbation yet produces :math:`U^3`.  So the rule cannot be "keep first
order, drop second order"; it has to be "keep :math:`U^3`, drop
:math:`U^4`" — which happens to take exactly one term out of the second-order
convolution and leave the rest behind.

That is why the recurrence looks lopsided.  It is not keeping "first order
plus a bit"; it is keeping *all of* :math:`U^3` and *none of* :math:`U^4`.

The scaling is checkable rather than asserted.  Both surviving terms must
hold a constant ratio to :math:`U^3` as the drive changes, and they do —
regenerated on every build:

.. exec-rst::

    import numpy as np
    from pycircuit.circuit import numeric
    from pycircuit.circuit.distortion import ExponentialNonlinearity, Harmonic

    IS, VT = 1e-6, 25e-3
    alpha = IS/VT
    w0 = 2*np.pi*1e4
    Y = lambda s: 1e-3 + s*1e-9 + alpha
    aG = lambda h, r: [np.asarray(r, dtype=complex)[0]/Y(1j*h.frequency((w0,)))]
    nl = ExponentialNonlinearity(IS, VT, port=0, n=1)

    rows = []
    for drive in (1e-6, 2e-6, 4e-6, 8e-6):
        X1 = np.asarray(aG(Harmonic((1,)), [drive]), dtype=complex)
        F2, F3, B1 = nl.harmonic_sources(X1, numeric)
        X2 = np.asarray(aG(Harmonic((2,)), [-F2[0]]), dtype=complex)
        mixing = B1(X2)[0]/2
        U = abs(X1[0])
        rows.append((U, abs(F3[0])/U**3, abs(mixing)/U**3))

    head = ['X1 (V)', 'route 1 / U^3', 'route 2 / U^3']
    w = [12, 15, 15]
    sep = ' '.join('='*x for x in w)
    print(sep)
    print(' '.join('%-*s' % (x, h) for x, h in zip(w, head)))
    print(sep)
    for U, a, b in rows:
        print(' '.join('%-*s' % (x, c) for x, c in
                       zip(w, ['%.3e' % U, '%.4e' % a, '%.4e' % b])))
    print(sep)

Both columns stay flat across a factor of eight in drive.  If either route
were secretly a different order in :math:`U`, its column would slope.

Several nonlinear devices at once
---------------------------------

A real circuit has more than one nonlinear device, and they interact: a
distortion product generated in one stage is amplified by the next, and may
partially cancel one generated there.  Nothing in the derivation above
assumed a single device — promoting :math:`b` and :math:`c` from scalars to
matrices covers the general case.

Entry :math:`b_{sk}` is the quadratic coefficient of the nonlinear current
injected at node :math:`s` and controlled by node :math:`k`, so the
nonlinearity is

.. math::

   f_s = \sum_k \left( b_{sk}\,x_k^2 + c_{sk}\,x_k^3 \right)

and the recurrence is unchanged — the products above simply become
matrix–vector products with elementwise powers.

One restriction comes with this form and is worth stating plainly: the sum
is over **self-terms only**.  A device whose current depends on the *product*
of two different controlling voltages, :math:`x_j x_k`, is outside the
formulation.

The worked example in the reference — a three-stage amplifier with reversed
nested Miller compensation, carrying three nonlinearities across two nodes,
one of them a nonlinear output conductance with a negative coefficient — is
reproduced in ``test_distortion``.

The fundamental is not corrected
--------------------------------

Note that :math:`\hat{X}_1` above carries no perturbation term.  At this
truncation the method gives the fundamental no nonlinear correction, so
**gain compression is not predicted**.  All four papers in this line share the
property; two of them derive the correction and then drop it at the final
assembly step.  It is a consequence of truncating at the third power of the
input, not an oversight.

Relation to the Volterra series
-------------------------------

The Volterra series is the established method for this problem, and the
obvious question is whether this one gives a different answer.  It does not.
Buonomo and Lo Schiavo prove that the two agree term by term at each order,
and that for a cubic nonlinearity the resulting HD2 and HD3 are *identical*,
not merely asymptotically close.

That is a useful fact for a simulator, because it turns into a test that
does not depend on the method's own authors being right.  ``test_distortion``
derives HD2 and HD3 a second time from Volterra kernels computed from
scratch, and requires symbolic agreement.  The kernels for
:math:`Y(s)v = u - (bv^2 + cv^3)` are

.. math::

   H_1(s) &= 1/Y(s) \\
   H_2(s_1,s_2) &= -H_1(s_1{+}s_2)\, b\, H_1(s_1)H_1(s_2) \\
   H_3(s_1,s_2,s_3) &= -H_1(s_1{+}s_2{+}s_3)
     \left[c\,H_1(s_1)H_1(s_2)H_1(s_3)
           + 2b\,\mathrm{sym}\{H_1(s_i)H_2(s_j,s_k)\}\right]

which is a genuinely different derivation — a kernel recursion rather than an
order-by-order harmonic balance — and it reaches the same expressions.

The practical advantage claimed for the perturbation form is that it needs
only transfer functions and Fourier coefficients, never explicit kernels, and
that it extends to non-polynomial nonlinearities (the exponential of a diode
or bipolar junction) for which the kernels become unwieldy.

Where the Taylor coefficients come from
---------------------------------------

The method needs :math:`b` and :math:`c` — the second and third Taylor
coefficients of each device characteristic about the operating point.  Every
other analysis in pycircuit needs only the first, which elements supply as
their ``G()`` stamp.

These are **derived from the element's own model function**, never declared
separately:

.. math::

   b = \tfrac{1}{2} i''(x_0), \qquad c = \tfrac{1}{6} i'''(x_0)

taken from the same ``eval_i_pure`` the element already exposes, by repeated
differentiation through the toolkit.  Under a symbolic toolkit that
differentiation is exact at every order.

The alternative — giving each element a method that states its Taylor
coefficients — was rejected deliberately.  It would state each device model
twice, and the two statements could silently disagree.  That is not
hypothetical: ``VCVS_limited`` carried exactly that defect for years, its
hand-written ``G()`` disagreeing with its ``i()`` in three separate ways while
every simulation still converged, just to a slightly wrong answer.  Deriving
the coefficients means they cannot drift, because they *are* the model.

Validity
--------

The method assumes :math:`\mathcal{G}f(x)` is small compared with
:math:`\mathcal{G}u`.  **No paper in this line quantifies that**, and only the
2005 one acknowledges that the convergence radius of the series is hard to
predict at all.

This matters more than it might appear.  The analysis will return a
confident, plausible, entirely wrong number for a circuit outside its range,
with nothing in the result to indicate it.  A later stage of the plan
therefore reports the ratio :math:`\|\mathcal{G}f(x_0)\| / \|\mathcal{G}u\|`
alongside every result, so that the assumption is visible rather than
implicit.

Validation against published closed forms
-----------------------------------------

The multi-node path,
:func:`~pycircuit.circuit.distortion.graded_response_mimo`, is checked against
**closed-form expressions printed by Buonomo & Lo Schiavo** (*Analog Integr.
Circ. Sig. Process.* 77:483–493, 2013) for both of that paper's worked
examples.  These are not comparisons against our own arithmetic, and nothing
is fitted: the reference is an independently-derived algebraic formula, and
the circuits were chosen by the authors, not by us.

Every number below is recomputed on each build.

**Worked example 2 — Tow–Thomas gm-C biquad** (their eq. 46), two nodes and
four cubic transconductor nonlinearities, against their eq. (48) for the third
harmonic at the output node:

.. exec-rst::

    import numpy as np
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)

    g1 = g2 = 31.26e-6; g3 = 625.2e-6; g4 = -625.2e-6
    C1 = C2 = 9.3054e-12
    alpha = -0.0535                      # every cubic tied to its linear
    g1c, g2c, g3c, g4c = (alpha*g for g in (g1, g2, g3, g4))
    q = lambda w: g3*g4 + C1*C2*w**2 + g2*C2*1j*w

    def solve(s, rhs):
        M = np.array([[-g2 + s*C1, -g4], [-g3, s*C2]], dtype=complex)
        return np.linalg.solve(M, np.asarray(rhs, dtype=complex))

    def cube(sp, n):
        return ((sp*sp).truncated(n) * sp).truncated(n)

    def ours(Xin, w, n=3):
        src = (GradedSpectrum.from_phasor(1, 1, g1*Xin)
               + cube(GradedSpectrum.from_phasor(1, 1, Xin), n).scaled(g1c))
        f = lambda x: GradedVector([
            cube(x[0], n).scaled(-g2c) + cube(x[1], n).scaled(-g4c),
            cube(x[0], n).scaled(-g3c)])
        return graded_response_mimo(
            solve, GradedVector([src, GradedSpectrum()]), f, (w,), n)

    def eq48(w, Xin):
        first = 3*(g4c*g1**3*g3**3 - g1c*q(w)**3)
        second = -g1**3*C2**2*w**2*(g3c*g4 + 3*g2c*C2*1j*w)
        return Xin**3*C2*1j*w*(first + second)/(4*q(3*w)*q(w)**3)

    Xin = 1e-3
    head = ['f (Hz)', 'ours', 'their eq. (48)', 'rel. difference']
    rows = []
    for f in (1e5, 1e6, 5e6, 1.06931e7, 2e7):
        w = 2*np.pi*f
        a, b = abs(ours(Xin, w)[0].phasor(3)), abs(eq48(w, Xin))
        rows.append(['%.4g' % f, '%.6e' % a, '%.6e' % b,
                     '%.1e' % (abs(a-b)/b)])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)

**Worked example 1 — three-stage RNMC amplifier** (their eq. 45), three nodes,
with **quadratic as well as cubic** terms, against their eqs. (41) and (42)
for :math:`HD_2` and :math:`HD_3`.  Their formulas are restricted to the
output transconductor's nonlinearities, so the comparison switches the others
off to match:

.. exec-rst::

    import numpy as np
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)

    gm1, g01 = 245e-6, 1/98e3
    gm2, g02 = 200e-6, 1/107e3
    gm3, g03 = 1e-3, 1/23e3
    gm3q, gm3c = 2e-3, 3e-3
    CL, C1, C2 = 10e-12, 2e-12, 1e-12

    p = lambda w: g03 + 1j*w*CL
    r = lambda w: g01*g02 + 1j*w*(C1*g02 + C2*(g02+gm2))
    q = lambda w: (g01*g02*g03
                   + 1j*w*(CL*g01*g02 + C1*(g02*g03 + gm2*gm3)
                           + C2*g03*(g02+gm2))
                   - w**2*((C1+C2)*CL*g02 + C2*CL*gm2))

    def solve(s, rhs):
        ## Open-loop limit of their eq. (45): k_in -> 1, k_3 -> 0,
        ## g_03e -> g_03.  Checked against their eq. (39) in the test suite.
        M = np.array([[g01 + (C1+C2)*s, -C2*s, -C1*s],
                      [gm2, g02, 0.0],
                      [0.0, gm3, -g03 - CL*s]], dtype=complex)
        return np.linalg.solve(M, np.asarray(rhs, dtype=complex))

    def ours(Xin, w, n=3):
        def f(x):
            sq = (x[1]*x[1]).truncated(n)
            return GradedVector([GradedSpectrum(), GradedSpectrum(),
                                 sq.scaled(gm3q)
                                 + (sq*x[1]).truncated(n).scaled(gm3c)])
        src = GradedVector([GradedSpectrum.from_phasor(1, 1, gm1*Xin),
                            GradedSpectrum(), GradedSpectrum()])
        return graded_response_mimo(solve, src, f, (w,), n)

    Xin = 1e-4
    head = ['f (Hz)', 'HD2 ours', 'HD2 eq. (41)', 'rel.',
            'HD3 ours', 'HD3 eq. (42)', 'rel.']
    rows = []
    for fr in (1e2, 1e3, 1e4, 1e5, 1e6, 1e7):
        w = 2*np.pi*fr
        g = ours(Xin, w)
        ## Their X_{1,3} is the *linearised* fundamental; the graded form
        ## additionally carries the fundamental's own U**3 correction.
        lin = -gm1*gm2*gm3*Xin/q(w)
        h2, h3 = abs(g[2].phasor(2)/lin), abs(g[2].phasor(3)/lin)
        e41 = abs(gm3q*gm1*gm2*Xin*p(w)**2*r(2*w)/(2*gm3*q(w)*q(2*w)))
        e42 = (abs(gm1**2*gm2**2*Xin**2*p(w)**3*r(3*w)
                   / (gm3*q(w)**2*q(3*w)))
               * abs(gm3c/4 - 1j*w*C1*gm3q**2*gm2/q(2*w)))
        rows.append(['%.4g' % fr, '%.4e' % h2, '%.4e' % e41,
                     '%.0e' % (abs(h2-e41)/e41), '%.4e' % h3,
                     '%.4e' % e42, '%.0e' % (abs(h3-e42)/e42)])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)

**Two tones on the same amplifier**, against their eq. (43) for the
:math:`2\omega_1-\omega_2` intermodulation product, at the paper's own tone
ratio :math:`f_2 = 0.9 f_1`:

.. exec-rst::

    import numpy as np
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)

    gm1, g01 = 245e-6, 1/98e3
    gm2, g02 = 200e-6, 1/107e3
    gm3, g03 = 1e-3, 1/23e3
    gm3q, gm3c = 2e-3, 3e-3
    CL, C1, C2 = 10e-12, 2e-12, 1e-12

    p = lambda w: g03 + 1j*w*CL
    r = lambda w: g01*g02 + 1j*w*(C1*g02 + C2*(g02+gm2))
    q = lambda w: (g01*g02*g03
                   + 1j*w*(CL*g01*g02 + C1*(g02*g03 + gm2*gm3)
                           + C2*g03*(g02+gm2))
                   - w**2*((C1+C2)*CL*g02 + C2*CL*gm2))

    def solve(s, rhs):
        M = np.array([[g01 + (C1+C2)*s, -C2*s, -C1*s],
                      [gm2, g02, 0.0],
                      [0.0, gm3, -g03 - CL*s]], dtype=complex)
        return np.linalg.solve(M, np.asarray(rhs, dtype=complex))

    def ours(Xin, w1, w2, n=3):
        src = (GradedSpectrum.from_phasor((1, 0), 1, gm1*Xin)
               + GradedSpectrum.from_phasor((0, 1), 1, gm1*Xin))
        def f(x):
            sq = (x[1]*x[1]).truncated(n)
            return GradedVector([GradedSpectrum(), GradedSpectrum(),
                                 sq.scaled(gm3q)
                                 + (sq*x[1]).truncated(n).scaled(gm3c)])
        return graded_response_mimo(
            solve, GradedVector([src, GradedSpectrum(), GradedSpectrum()]),
            f, (w1, w2), n)

    Xin = 1e-4
    head = ['f1 (Hz)', 'IM3 ours', 'their eq. (43)', 'rel. difference']
    rows = []
    for f1 in (1e2, 1e3, 1e4, 1e5, 1e6, 1e7):
        w1 = 2*np.pi*f1
        w2 = 0.9*w1
        got = (abs(ours(Xin, w1, w2)[2].phasor((2, -1)))
               / abs(-gm1*gm2*gm3*Xin/q(w1)))
        want = abs(gm1**2*gm2**2*Xin**2*p(w1)**2*np.conj(p(w2))*r(2*w1-w2)
                   / (gm3*q(w1)*np.conj(q(w2))*q(2*w1-w2))) \
            * abs(3*gm3c/4 - C1*gm3q**2*gm2*(1j*w1/q(2*w1)
                                             + 1j*(w1-w2)/q(w1-w2)))
        rows.append(['%.4g' % f1, '%.6e' % got, '%.6e' % want,
                     '%.1e' % (abs(got-want)/want)])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)

The two-tone index is a pair :math:`(m, n)` meaning :math:`m\omega_1 +
n\omega_2`, and because harmonics are stored two-sided the product
:math:`(2,-1)` arises as :math:`(1,0)+(1,0)+(0,-1)` — a componentwise sum,
with no case analysis over which combinations of sum and difference
frequencies land where.

Three points about what these tables do and do not establish.

**They agree to floating point, not merely closely.**  That is the expected
outcome and the reason the gate was set there: the published forms *are*
consistent :math:`U^3` truncations of the same expansion, so anything short of
round-off would mean one of the two is wrong.  A comparison that agreed to
"within a few percent" would have been a failure, not a success.

**The amplifier's matrix is reconstructed, and that is checked separately.**
The paper prints only the feedback configuration (its eq. 45); the open-loop
matrix used here is its limit.  A reconstruction is an assumption, so it has
its own gate against their eq. (39), which gives the linearised solution at
all three nodes — otherwise a wrong matrix would surface as a distortion
discrepancy and be misdiagnosed as one.

**The ratios use the linearised fundamental**, which is what the papers'
:math:`\hat X_{1,3}` denotes.  Dividing by the graded fundamental instead
disagrees by about half a percent at low frequency — not an error, but the
fundamental's own :math:`U^3` correction, which the published form drops and
this one carries.  The distinguishing evidence is that the gap falls by
:math:`100\times` per decade of drive, the signature of an :math:`O(U^3)` term
against an :math:`O(U)` fundamental; see `The fundamental is not corrected`_.

.. warning::

   **These closed forms are moduli of products, so they are blind to a whole
   class of error.**  Every one of them is written as
   :math:`\left|\,\cdot\,\right|` around a product and ratio of transfer
   functions, and conjugating any factor leaves a modulus unchanged.

   Two consequences, both observed here rather than supposed.  The conjugates
   the paper writes as :math:`\tilde p`, :math:`\tilde q` make **no
   difference** to the comparison — dropping them changes nothing.  And the
   pencil printed as eq. (46) has right-half-plane poles without any of these
   comparisons noticing.

   So agreement to :math:`10^{-15}` against these formulas is strong evidence
   about **magnitudes** and no evidence at all about **phase**.  Anything that
   depends on phase — cancellation between contributions, stability, a
   time-domain waveform — needs a different check.

   That check exists and is described under `Phase`_ below, against a
   time-domain integration.  The gap is closed; this warning stands because
   it explains *why* a second kind of check was needed.

Phase
-----

A waveform is not a modulus, so integrating the circuit's own differential
equations tests what none of the published closed forms can.  Regenerated on
every build:

.. exec-rst::

    import numpy as np
    from scipy.integrate import solve_ivp
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)

    g1 = g2 = 31.26e-6; g3 = 625.2e-6; g4 = -625.2e-6
    C1 = C2 = 9.3054e-12
    alpha = -0.0535
    g1c, g2c, g3c, g4c = (alpha*g for g in (g1, g2, g3, g4))
    w0 = np.sqrt(-g3*g4/(C1*C2))

    def solve(s, rhs):
        M = np.array([[g2 + s*C1, -g4], [-g3, s*C2]], dtype=complex)
        return np.linalg.solve(M, np.asarray(rhs, dtype=complex))

    def cube(sp, n):
        return ((sp*sp).truncated(n) * sp).truncated(n)

    Xin = 0.10
    src = (GradedSpectrum.from_phasor(1, 1, g1*Xin)
           + cube(GradedSpectrum.from_phasor(1, 1, Xin), 11).scaled(g1c))
    f = lambda x: GradedVector([cube(x[0], 11).scaled(-g2c)
                                + cube(x[1], 11).scaled(-g4c),
                                cube(x[0], 11).scaled(-g3c)])
    got = graded_response_mimo(solve, GradedVector([src, GradedSpectrum()]),
                               f, (w0,), 11)

    def rhs(t, y):
        x1, x2 = y
        u = Xin*np.cos(w0*t)
        return [(-g2*x1 + g4*x2 + g1*u + g1c*u**3
                 + g2c*x1**3 + g4c*x2**3)/C1,
                (g3*x1 + g3c*x1**3)/C2]

    T = 2*np.pi/w0
    tend, keep = 200*T, 32
    ## Window starts at a whole number of drive periods, so the FFT's phase
    ## reference is cos(w0 t) -- the same one the phasors use.
    ts = np.linspace(tend - keep*T, tend, keep*256, endpoint=False)
    sol = solve_ivp(rhs, (0, tend), [0.0, 0.0], t_eval=ts, method='DOP853',
                    rtol=1e-12, atol=1e-18, max_step=T/60)
    S = np.fft.rfft(sol.y[0])/len(sol.y[0])

    head = ['harmonic', 'magnitude error', 'phase error (deg)']
    rows = []
    for m in (1, 3, 5):
        a, b = got[0].phasor(m), 2*S[m*keep]
        rows.append(['%d' % m, '%.1e' % (abs(abs(a)-abs(b))/abs(b)),
                     '%.4f' % np.angle(a/b, deg=True)])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)
    print('')

    rebuilt = np.full_like(ts, float(np.real(got[0].phasor(0))))
    for m in range(1, 12):
        X = got[0].phasor(m)
        if X != 0:
            rebuilt = rebuilt + np.real(X*np.exp(1j*m*w0*ts))
    err = np.max(np.abs(rebuilt - sol.y[0]))/np.max(np.abs(sol.y[0]))
    print('Waveform rebuilt from the phasors, against the integrated one:')
    print('peak error ``%.2e`` of peak amplitude -- every harmonic\'s'
          % err)
    print('magnitude *and* phase at once.')
    print('')

    ## The criterion has to be shown capable of failing.  These leave |X|
    ## bit-identical, so a magnitude-only comparison accepts all of them.
    ref = got[0].phasor(3)
    b = 2*S[3*keep]
    muts = [('sign flip', -ref), ('90 degree rotation', 1j*ref),
            ('conjugation', np.conj(ref))]
    print('Mutations that preserve ``|X|`` exactly:')
    print('')
    for name, z in muts:
        print('* %s -- magnitude error ``%.1e``, phase error ``%+.1f`` deg'
              % (name, abs(abs(z)-abs(b))/abs(b), np.angle(z/b, deg=True)))

Every mutation above leaves the magnitude error unchanged to the digit shown,
so a magnitude-only comparison accepts all three.  That is the point of
including them: it distinguishes "the phases agree" from "the phases agree up
to an offset nobody checked".

The test is ``test_biquad_phase_and_waveform_match_the_ode_integration``.

Higher order on a multi-node circuit
------------------------------------

The published forms above are :math:`U^3` truncations.  Raising the truncation
on a *multi-node* circuit is checked separately, against a direct integration
of the biquad's own differential equations — a reference that shares the
circuit with the analysis and nothing else.  Errors in the third harmonic at
the resonant node:

.. exec-rst::

    import numpy as np
    from scipy.integrate import solve_ivp
    from pycircuit.circuit.distortion import (graded_response_mimo,
                                              GradedSpectrum, GradedVector)

    g1 = g2 = 31.26e-6; g3 = 625.2e-6; g4 = -625.2e-6
    C1 = C2 = 9.3054e-12
    alpha = -0.0535
    g1c, g2c, g3c, g4c = (alpha*g for g in (g1, g2, g3, g4))
    w0 = np.sqrt(-g3*g4/(C1*C2))

    def solve(s, rhs):
        ## Damping sign flipped from eq. (46): as printed the pencil has
        ## right-half-plane poles, so it has no periodic steady state to
        ## integrate toward.  See the note below the table.
        M = np.array([[g2 + s*C1, -g4], [-g3, s*C2]], dtype=complex)
        return np.linalg.solve(M, np.asarray(rhs, dtype=complex))

    def cube(sp, n):
        return ((sp*sp).truncated(n) * sp).truncated(n)

    def ours(Xin, n):
        src = (GradedSpectrum.from_phasor(1, 1, g1*Xin)
               + cube(GradedSpectrum.from_phasor(1, 1, Xin), n).scaled(g1c))
        f = lambda x: GradedVector([cube(x[0], n).scaled(-g2c)
                                    + cube(x[1], n).scaled(-g4c),
                                    cube(x[0], n).scaled(-g3c)])
        return graded_response_mimo(
            solve, GradedVector([src, GradedSpectrum()]), f, (w0,), n)

    def reference(Xin, cycles=200, keep=32):
        def rhs(t, y):
            x1, x2 = y
            u = Xin*np.cos(w0*t)
            return [(-g2*x1 + g4*x2 + g1*u + g1c*u**3
                     + g2c*x1**3 + g4c*x2**3)/C1,
                    (g3*x1 + g3c*x1**3)/C2]
        T = 2*np.pi/w0
        tend = cycles*T
        ts = np.linspace(tend - keep*T, tend, keep*256, endpoint=False)
        s = solve_ivp(rhs, (0, tend), [0.0, 0.0], t_eval=ts, method='DOP853',
                      rtol=1e-11, atol=1e-18, max_step=T/40)
        S = np.fft.rfft(s.y[0])/len(s.y[0])
        return 2*abs(S[3*keep])

    ORDERS = (3, 5, 7, 9, 11)
    head = ['drive X_in'] + ['U^%d' % n for n in ORDERS] + ['monotone?']
    rows = []
    for Xin in (0.05, 0.10, 0.20, 0.30):
        want = reference(Xin)
        errs = [abs(abs(ours(Xin, n)[0].phasor(3)) - want)/want
                for n in ORDERS]
        ok = all(b < a for a, b in zip(errs, errs[1:]))
        rows.append(['%.2f' % Xin] + ['%.2e' % e for e in errs]
                    + ['yes' if ok else 'NO'])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)

**Below the bound every added order helps, and by orders of magnitude.**  The
last row is above the bound, and there the sequence stops being monotone: it
improves overall but wanders on the way, which is the same asymptotic
behaviour the single-node case shows (see :doc:`distortion_limits`).  Monotone
improvement is a property of being inside the radius of convergence, not of
the method.

One caveat on reading the smallest numbers in the table.  At the lowest drive
and highest truncation the discrepancy is down at the level of **the
reference's own accuracy**, not the analysis's: tightening the integrator
(more cycles, smaller tolerance) moves that entry by an order of magnitude
while leaving every other entry unchanged.  Those cells bound how well the
comparison can be made, not how well the method does.  The columns to trust
for the method are the ones where the error is still comfortably above the
integrator's floor.

.. note::

   **The pencil printed as eq. (46) has right-half-plane poles**,
   :math:`\mathrm{Re}\,s = +1.68\times10^{6}`.  The table above therefore
   flips the sign of the damping term, because an unstable system has no
   periodic steady state for the reference to integrate toward.

   This is invisible in the source: every published result for this circuit
   is a modulus, and :math:`|q(j\omega)|` is unchanged by conjugating the
   pencil — which is why the comparisons against eqs. (47) and (48) above
   pass using the matrix exactly as printed.

   It does **not** follow that the sign never affects a magnitude.  It
   affects this one by 1–2%: the graded computation *adds* complex quantities
   across harmonics, and addition is not invariant under conjugation.  Only
   products and ratios of :math:`q` are, which is precisely the form the
   published closed forms take.

Limitations of the present implementation
-----------------------------------------

* Single tone.  :func:`~pycircuit.circuit.distortion.intermodulation_response`
  covers two tones, but only on the scalar path.
* The fundamental is uncorrected in the published second-order form, so gain
  compression and AM-AM effects are outside *that* model — though the graded
  forms do carry the correction; see `The fundamental is not corrected`_.
* **No longer a limitation:** several nonlinear elements across several nodes,
  and cross-terms :math:`x_j x_k` between different controlling variables.
  The source formulation assumes self-terms only and the scalar
  :func:`~pycircuit.circuit.distortion.graded_response` inherited that, but
  :func:`~pycircuit.circuit.distortion.graded_response_mimo` takes the
  nonlinearity as a callable in the graded ring, where a cross term costs
  nothing extra.  No published example in this line uses one, so the capability
  is **untested against an external reference** — the representation permits
  it, which is not the same as it having been checked.

Devices that are not polynomials
--------------------------------

A diode or bipolar junction is exponential, not cubic:

.. math::

   f(v) = I_S\left(e^{v/V_T} - 1 - v/V_T\right)

with the constant and linear terms subtracted so that the nonlinearity is
strictly nonlinear, its linear part having been absorbed into the circuit's
admittance as the junction conductance :math:`I_S/V_T`.

The usual move is to fit a cubic and reuse the machinery above.  That is what
the Volterra literature generally does, and it is **not adequate here**: the
2005 reference measures a cubic-truncated exponential converging to a second
harmonic materially below the true value — around 20% on its Fig. 3 — while
still looking entirely reasonable.

It is also unnecessary, because the exponential case has an exact closed form
at *every* harmonic order.  The Fourier coefficients of
:math:`\exp(a\cos\theta)` are modified Bessel functions of the first kind,
so

.. math::

   F_m = 2 I_S\, I_m\!\left(|X_1|/V_T\right), \qquad m > 1

with no truncation and no fitting.  The recurrence is untouched — only the
supplier of :math:`F_2`, :math:`F_3` and :math:`B_1` changes, which is why
those three quantities are the interface between the nonlinearity and the
method.

Two consequences worth stating.  At small drive the exponential and cubic
treatments must converge, since the cubic fit *is* the exponential's Taylor
expansion; ``test_distortion`` pins that, and it is the check that would catch
a wrong argument scaling in the Bessel functions, which a working-amplitude
comparison alone would not.  At working drive they must diverge, in the
direction of the exact treatment predicting *more* distortion.

Node quantities are not output quantities
------------------------------------------

The harmonics this analysis produces are node voltages.  Distortion measured
at the node controlling a nonlinearity is generally **not** the distortion at
the circuit's output: the mapping from one to the other carries a feedforward
path and a frequency-dependent factor of its own.

For the bipolar amplifier of the 2005 reference the output current is

.. math::

   y = \mathcal{E}u + \mathcal{N}x, \qquad
   \mathcal{E} = -\beta_F/R, \qquad
   \mathcal{N} = (\beta_F/R)(1 + j\omega RC)

and referring HD to :math:`y` rather than to :math:`v_{be}` changes it by a
factor of ten for the published component values.  Ten is a plausible number
to see and no error is raised, which is why ``test_distortion`` pins the
distinction explicitly.

Two tones, and intermodulation
------------------------------

Harmonic distortion tells only part of the story.  Drive a circuit with two
closely spaced tones :math:`\omega_1` and :math:`\omega_2` and the
nonlinearity generates products at :math:`2\omega_1 - \omega_2` and
:math:`2\omega_2 - \omega_1` — which land *inside* the passband and cannot be
filtered away.  That is usually the figure a designer actually cares about.

Nothing about the method changes.  The recurrence is the same: evaluate the
nonlinearity on a known solution, push each component back through
:math:`\mathcal{G}` at its own frequency.  Only the frequency index becomes
richer — a pair :math:`(m, n)` denoting :math:`m\omega_1 + n\omega_2` in place
of a single harmonic number.  That is why :class:`Harmonic` has been a tuple
since the first commit of this module.

Which coefficient reaches the product, and how
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the part worth understanding, because it is not symmetric.

**Only the cubic coefficient reaches** :math:`2\omega_1 - \omega_2`
**directly.**  Squaring two tones produces components at :math:`2\omega_1`,
:math:`2\omega_2`, :math:`\omega_1 \pm \omega_2` and DC — never at
:math:`2\omega_1 - \omega_2`.  So at first perturbation order the quadratic
coefficient contributes nothing at all.

It arrives at *second* order instead, by mixing: the quadratic nonlinearity
creates products at :math:`2\omega_1` and :math:`\omega_1 - \omega_2`, and
those mix back up against a tone.  There are exactly two routes,

.. math::

   2\omega_1 - \omega_2 = (2\omega_1) + (-\omega_2)
                        = (\omega_1) + (\omega_1 - \omega_2)

and the implementation enumerates both.  The 2005 reference writes the same
thing as a double convolution over all index pairs; enumerating is equivalent
once the sum is truncated to third order in the input.

The two contributions can very nearly cancel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the reference amplifier they do.  The first- and second-order
contributions differ by 0.74 dB and partially cancel, leaving a total some
22 dB below either one.

Two consequences follow, and both matter more than the headline number.  An
analysis truncated at first perturbation order would report IM3 roughly 21 dB
*too high* here — not a refinement, a wrong answer.  And any implementation
error in either contribution is amplified in the total, which is why
``test_distortion`` checks the two contributions separately as well as their
sum: checking only the total would let compensating errors through, and
checking only the terms would miss a sign error in how they combine.

Not implemented: two tones on an exponential device
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-tone derivation in the reference covers cubic polynomials only, and
:class:`~pycircuit.circuit.distortion.ExponentialNonlinearity` raises rather
than returning a number from a derivation nobody has checked.

The obstacle is *not* the mathematics, which is easy.  The exponential
factorises over a sum of tones,

.. math::

   e^{a_1\cos\theta_1 + a_2\cos\theta_2}
     = \sum_{m,n} I_m(a_1)\, I_n(a_2)\, e^{j(m\theta_1 + n\theta_2)}

so the coefficient at :math:`(m,n)` is an ordinary *product* of Bessel
functions — no special two-argument function is involved.

What is missing is phase, and the shape of that gap is worth understanding
because it is a good example of a defect that hides.

The Jacobi–Anger expansion is stated for a real cosine drive, so a drive
:math:`A e^{j\phi}` gives harmonic :math:`m` a factor :math:`e^{jm\phi}`.
Omit it and every magnitude is still exactly right while every argument is
wrong.  For **one tone and one device that is unobservable**: the same factor
multiplies :math:`F_m` and the second-order mixing term, so it becomes a
common multiplier on the whole harmonic — precisely a choice of time origin —
and cancels in every ratio.  No single-tone answer was ever wrong because of
it, and no single-device test could have detected it.

It stops being free as soon as the phases cannot all be absorbed into one time
origin.  Two cases do that: two tones (only one phase can be absorbed), and
**two nonlinear devices whose controlling nodes sit at different phases**.
The second is reachable today with
:class:`~pycircuit.circuit.distortion.CompositeNonlinearity`, and
``test_distortion`` uses it — an exponential device and a cubic one on two
capacitively coupled nodes about 65° apart — as the regression test that the
single-device configuration could not provide.

So extending to two tones means carrying complex amplitudes through the Bessel
path rather than magnitudes: a real piece of work, but a bounded one, and not
a mathematical obstacle.

.. _when-more-terms-stop-helping:





Why this method suits a *symbolic* tool
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Everything above weighs accuracy.  For a symbolic circuit analyser there is a
prior consideration that settles the matter more firmly: **what the answer
looks like as an expression.**

Take the scalar circuit :math:`Yx + bx^2 = u` and carry both constructions
through with every quantity symbolic.  The perturbation terms come out as
single monomials:

.. code-block:: text

    x(0) =      u/Y         1 term,   1 op
    x(1) =  -  b u^2/Y^3    1 term,   5 ops
    x(2) =    2 b^2 u^3/Y^5 1 term,   6 ops
    x(3) =  - 5 b^3 u^4/Y^7 1 term,   7 ops
    x(4) =   14 b^4 u^5/Y^9 1 term,   6 ops

An iterative alternative -- substituting the current estimate back into the
circuit equation and repeating -- doubles in length every pass:

.. code-block:: text

    iterate 1:   2 terms,    6 ops
    iterate 2:   4 terms,   19 ops
    iterate 3:   8 terms,   47 ops
    iterate 4:  16 terms,  103 ops
    iterate 5:  32 terms,  215 ops

Two things follow.

**The perturbation terms are readable.**  :math:`-5b^3u^4/Y^7` states directly
how the third-order contribution scales with the nonlinearity, the drive and
the admittance — which is the entire point of computing it symbolically rather
than numerically.  A 32-term iterate is a number generator wearing algebra.
(The coefficients 1, 1, 2, 5, 14 are the Catalan numbers, the expected
signature of a quadratic fixed point; the structure is clean, not accidental.)

**Iteration cannot be carried symbolically in any case.**  Each pass
substitutes the whole previous expression into :math:`f`, so a cubic
nonlinearity cubes the expression size every time — the doubling above is the
mild quadratic case.  Iterative schemes that work in the time domain, by
sampling and transforming, are not symbolic operations at all.

(:doc:`distortion_limits` measures what such an approach buys numerically.
The short answer is that it converges to the exact steady state where it
converges at all, and that this is not available symbolically.)

So the truncation is not merely a choice about accuracy or consistency.  The
closed-form Fourier coefficients of `Where the Taylor coefficients come from`_
keep every order at bounded degree, and **that boundedness is what makes
symbolic distortion analysis possible at all**.  An iterative method, however
accurate, forfeits it.
