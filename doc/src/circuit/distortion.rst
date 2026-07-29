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

Limitations of the present implementation
-----------------------------------------

* Cubic polynomial nonlinearity, single tone, one nonlinear element.
* Cross-terms :math:`x_j x_k` between different controlling variables are not
  covered — the source formulation assumes self-terms only.
* The fundamental is uncorrected, so gain compression and AM-AM effects are
  outside the model.

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

When more terms stop helping
----------------------------

The natural response to a prediction that is off by 20% is to take another
term.  For this method that works up to a point and then stops working
entirely, and the boundary is sharp enough to state.

The perturbation series is a fixed-point iteration:
:math:`x \leftarrow \mathcal{G}(u - f(x))`.  Its terms shrink only if the
iteration is a contraction, i.e. if
:math:`|\mathcal{G}\,f'(x)| < 1`.  For a junction whose own conductance
dominates the node, :math:`\mathcal{G} \approx 1/(I_S/V_T)` and
:math:`f'(v) = (I_S/V_T)(e^{v/V_T}-1)`, so the contraction factor at a drive
amplitude :math:`a = |X_1|/V_T` is simply

.. math::

   |\mathcal{G}f'| \approx e^{a} - 1

which crosses 1 at :math:`a = \ln 2 \approx 0.693`.

Running the iteration and watching successive terms confirms it — magnitudes
relative to the linear solution, order by order, regenerated on every build:

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from pycircuit.circuit import numeric

    IS_D = 1e-13
    VT = numeric.kboltzmann * 300 / numeric.qelectron
    Ibias, Rload, Cload, f0 = 1e-3, 1e3, 1e-9, 1e4
    v0 = brentq(lambda v: v/Rload + IS_D*np.expm1(v/VT) - Ibias, 0, 1)
    ISeff = IS_D*np.exp(v0/VT); alpha = ISeff/VT
    w0 = 2*np.pi*f0
    NH, N = 48, 512
    Gm = np.array([1.0/(1/Rload + 1j*m*w0*Cload + alpha) for m in range(NH+1)])
    f_nl = lambda v: ISeff*(np.expm1(v/VT) - v/VT)

    rows = []
    for Isig in (2e-5, 1e-4, 3e-4, 6e-4):
        U = np.zeros(NH+1, dtype=complex); U[1] = Isig
        Xk = Gm*U; x0 = abs(Xk[1]); a = x0/VT
        mags = []
        for _ in range(5):
            pad = np.concatenate([Xk, np.zeros(N//2+1-len(Xk))])
            xt = np.real(np.fft.irfft(pad*N/2))
            F = np.fft.rfft(f_nl(xt))/(N/2); F[0] /= 2
            Xn = Gm*(U - F[:NH+1])
            mags.append((abs(Xn[1]-Xk[1])+abs(Xn[2]-Xk[2])+abs(Xn[3]-Xk[3]))/x0)
            Xk = Xn
        rows.append((a, np.exp(a)-1, mags))

    w = [9, 13, 46]
    sep = ' '.join('=' * x for x in w)
    print(sep)
    print('%-9s %-13s %s' % ('a', 'contraction', 'successive term magnitudes'))
    print(sep)
    for a, contraction, mags in rows:
        print('%-9.3f %-13.3f %s' % (a, contraction,
              ', '.join('%.1e' % m for m in mags)))
    print(sep)

Below :math:`\ln 2` the terms fall geometrically and extra orders buy a great
deal — at :math:`a = 0.044` each term is some twenty times smaller than the
last.  Approaching :math:`\ln 2` they barely fall and extra orders buy little.
Above it **the terms grow**: the series diverges, and no number of additional
terms recovers anything.  A prediction that is 50% wrong there cannot be
repaired by working harder at it.

For a junction at room temperature :math:`a < \ln 2` means a signal swing
below roughly 17 mV.

Two practical consequences.  Increasing the *harmonic* order is a different
axis — it is ordinary Fourier truncation and always converges — but it is
subordinate to this bound, since more harmonics of a diverging series are
still diverging.  And when a prediction is poor, the question to ask first is
which truncation is responsible: approximating an exponential device by a
cubic is a separate error from truncating the perturbation series, and
:class:`~pycircuit.circuit.distortion.ExponentialNonlinearity` removes the
former exactly while leaving the latter untouched.

.. note::

   None of the source papers gives a quantitative validity bound; the 2005
   one states only that the convergence radius is hard to predict.  The
   threshold above is specific to an exponential device on a
   junction-dominated node, derived and then measured here.  It is not a
   general result, but it is a concrete one, and
   :meth:`~pycircuit.circuit.distortion.DistortionSolution.perturbation_ratio`
   reports a computable proxy for it on every solve.

Past the bound, a better model can give a worse answer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Comparing the polynomial treatments against the exact one, all three referred
to a transient simulation of the same diode, makes the point sharply.
Relative error, regenerated on every build:

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from pycircuit.circuit import SubCircuit, gnd, R, C, ISin, numeric
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.distortion import (harmonic_response,
                                              PolynomialNonlinearity,
                                              ExponentialNonlinearity)

    IS_D = 1e-13
    VT = numeric.kboltzmann * 300 / numeric.qelectron
    Ibias, Rload, Cload, f0 = 1e-3, 1e3, 1e-9, 1e4
    cm.default_toolkit = numeric
    v0 = brentq(lambda v: v/Rload + IS_D*np.expm1(v/VT) - Ibias, 0, 1)
    ISeff = IS_D*np.exp(v0/VT); alpha = ISeff/VT
    w0 = 2*np.pi*f0
    Y = lambda s: 1/Rload + s*Cload + alpha
    apply_G = lambda h, r: [np.asarray(r, dtype=complex)[0]
                            / Y(1j*h.frequency((w0,)))]

    models = [('cubic',   PolynomialNonlinearity.taylor_of_exponential(ISeff, VT, 3)),
              ('quintic', PolynomialNonlinearity.taylor_of_exponential(ISeff, VT, 5)),
              ('exact',   ExponentialNonlinearity(ISeff, VT, port=0, n=1))]

    def measure(Isig, per=16, pts=256, keep=8):
        c = SubCircuit(toolkit=numeric); n = c.add_node('n1')
        c['I'] = ISin(gnd, n, io=Ibias, ia=Isig, freq=f0)
        c['R'] = R(n, gnd, r=Rload); c['C'] = C(n, gnd, c=Cload)
        c['D'] = Diode(n, gnd, IS=IS_D)
        res = Transient(c, toolkit=numeric).solve(
            refnode=gnd, tend=per/f0, timestep=1.0/(f0*pts))
        w = res.v(n); t = np.asarray(w.x[0]); v = np.asarray(w.y)
        g = np.linspace(t[-1]-keep/f0, t[-1], keep*pts, endpoint=False)
        S = np.fft.rfft(np.interp(g, t, v))/(keep*pts)
        return abs(S[2*keep])/abs(S[keep]), abs(S[3*keep])/abs(S[keep])

    rows = []
    for Isig in (2e-5, 1e-4, 4e-4, 1e-3):
        m2, m3 = measure(Isig)
        errs = []
        for _, nl in models:
            sol = harmonic_response(apply_G, [Isig], nl, tones=(w0,),
                                    toolkit=numeric)
            errs += [abs(m2/sol.HD2(0)-1)*100, abs(m3/sol.HD3(0)-1)*100]
            ratio = sol.perturbation_ratio(0)
        a = abs(apply_G(__import__('pycircuit.circuit.distortion',
                fromlist=['Harmonic']).Harmonic((1,)), [Isig])[0])/VT
        rows.append((a, ratio, errs))

    head = ['a', 'ratio', 'HD2 cub', 'HD3 cub', 'HD2 qui', 'HD3 qui',
            'HD2 exact', 'HD3 exact']
    order = [0, 1, 2, 3, 4, 5]
    w = [7, 7, 9, 9, 9, 9, 10, 10]
    sep = ' '.join('='*x for x in w)
    print(sep)
    print(' '.join('%-*s' % (x, h) for x, h in zip(w, head)))
    print(sep)
    for a, ratio, e in rows:
        cells = ['%.3f' % a, '%.3f' % ratio] + ['%.2f%%' % e[i] for i in order]
        print(' '.join('%-*s' % (x, c) for x, c in zip(w, cells)))
    print(sep)

Two things stand out, and neither is what one would guess.

**The quintic fit essentially *is* the exact treatment.** Below the bound its
errors match the Bessel result to the second decimal.  Fifth order captures
practically all of the device-model truncation, which means the residual
error at :math:`a \approx 0.88` — still around 13% — is **not** the device
model at all.  It is the perturbation series being truncated at second order,
and no amount of further polynomial refinement touches it.  That is the clean
separation of the two truncations: raise the polynomial order to remove the
first, and only then is what remains attributable to the second.

**Past the bound the ordering inverts.** At the largest drive the crude cubic
is the *closest* to the simulation and the exact model the furthest.  Nothing
is wrong with the exact model: it feeds a diverging iteration exactly, while
the cubic under-represents the device in a direction that partly offsets the
divergence.  Quintic and exact bracket the cubic monotonically in the wrong
direction, which is what confirms the cubic's apparent accuracy is
coincidence rather than merit — it would not survive a change of operating
point.

(At the smallest drive all three sit near a tenth of a percent, which is the
transient simulation's own numerical floor rather than a property of any of
them; the block above deliberately uses a cheap transient so the page builds
quickly.)

The lesson is worth stating because it inverts the usual intuition: outside
the convergence bound, improving an input can degrade the output, and
apparent agreement carries no information.  This is exactly why
:meth:`~pycircuit.circuit.distortion.DistortionSolution.perturbation_ratio`
is reported alongside every result rather than left for the caller to
wonder about.

What the truncation is actually truncating
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Keeping more terms makes it worse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second-order step produces the third harmonic through a convolution whose
terms include :math:`B_0X_3`, :math:`B_1X_2`, :math:`B_2X_1`, :math:`B_3X_0`
and their conjugates.  The recurrence above keeps **only** :math:`B_1X_2`,
the sources having shown every other term to be :math:`O(U^4)` where that one
is :math:`O(U^3)`.

That looks like an approximation worth removing, and
:func:`~pycircuit.circuit.distortion.harmonic_response_spectral` removes it —
it forms the first-order correction as a waveform, multiplies by the
derivative as a waveform and transforms back, so *every* convolution term up
to a chosen cutoff survives.  The cutoff is a free parameter, so the natural
experiment is to double it.

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from pycircuit.circuit import SubCircuit, gnd, R, C, ISin, numeric
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.distortion import (harmonic_response_spectral,
                                              harmonic_response, Harmonic,
                                              PolynomialNonlinearity,
                                              ExponentialNonlinearity)

    IS_D = 1e-13
    VT = numeric.kboltzmann*300/numeric.qelectron
    Ib, Rl, Cl, f0 = 1e-3, 1e3, 1e-9, 1e4
    cm.default_toolkit = numeric
    v0 = brentq(lambda v: v/Rl + IS_D*np.expm1(v/VT) - Ib, 0, 1)
    ISe = IS_D*np.exp(v0/VT); al = ISe/VT
    w0 = 2*np.pi*f0
    Y = lambda s: 1/Rl + s*Cl + al
    aG = lambda h, r: [np.asarray(r, dtype=complex)[0]/Y(1j*h.frequency((w0,)))]
    fe  = lambda v: ISe*(np.expm1(v/VT) - v/VT)
    fpe = lambda v: (ISe/VT)*np.expm1(v/VT)

    def meas(Isig, per=16, pts=256, keep=8):
        c = SubCircuit(toolkit=numeric); n = c.add_node('n1')
        c['I'] = ISin(gnd, n, io=Ib, ia=Isig, freq=f0)
        c['R'] = R(n, gnd, r=Rl); c['C'] = C(n, gnd, c=Cl)
        c['D'] = Diode(n, gnd, IS=IS_D)
        res = Transient(c, toolkit=numeric).solve(
            refnode=gnd, tend=per/f0, timestep=1.0/(f0*pts))
        w = res.v(n); t = np.asarray(w.x[0]); v = np.asarray(w.y)
        g = np.linspace(t[-1]-keep/f0, t[-1], keep*pts, endpoint=False)
        S = np.fft.rfft(np.interp(g, t, v))/(keep*pts)
        return abs(S[2*keep])/abs(S[keep]), abs(S[3*keep])/abs(S[keep])

    rows = []
    for Isig in (2e-5, 1e-4, 4e-4):
        m2, m3 = meas(Isig)
        trunc = harmonic_response(aG, [Isig],
                    ExponentialNonlinearity(ISe, VT, port=0, n=1),
                    tones=(w0,), toolkit=numeric)
        cells = [abs(m2/trunc.HD2(0)-1)*100, abs(m3/trunc.HD3(0)-1)*100]
        for NH in (3, 6):
            s = harmonic_response_spectral(aG, [Isig], fe, fpe, (w0,), numeric,
                                           n_harmonics=NH)
            cells += [abs(m2/s.HD2(0)-1)*100, abs(m3/s.HD3(0)-1)*100]
        a = abs(aG(Harmonic((1,)), [Isig])[0])/VT
        rows.append((a, cells))

    head = ['a', 'HD2 trunc', 'HD3 trunc', 'HD2 N=3', 'HD3 N=3',
            'HD2 N=6', 'HD3 N=6']
    w = [7, 11, 11, 10, 10, 10, 10]
    sep = ' '.join('='*x for x in w)
    print(sep)
    print(' '.join('%-*s' % (x, h) for x, h in zip(w, head)))
    print(sep)
    for a, cells in rows:
        c = ['%.3f' % a] + ['%.2f%%' % v for v in cells]
        print(' '.join('%-*s' % (x, v) for x, v in zip(w, c)))
    print(sep)

Two results, and both cut against the instinct that prompted the experiment.

**Doubling the cutoff buys almost nothing.**  Harmonics above the third
barely feed back at second order, so ``N=6`` and ``N=3`` agree to within a
fraction of a percent.  The harmonic cutoff is simply not where the error
lives.

**Keeping the discarded convolution terms makes the answer markedly worse**,
and at moderate drive catastrophically so — the second-harmonic error grows
by more than an order of magnitude.  This is not a defect in either
implementation: the two agree as :math:`U \to 0`, exactly as they must, and
the discrepancy scales as :math:`U^2` relative to the leading term, which is
the signature of the :math:`O(U^4)` contributions being restored.

Why remains **open**, and an earlier version of this page asserted a reason it
should not have.  The tempting explanation is that the restored contributions
are incomplete — they are the :math:`O(U^4)` terms from second-order
perturbation, while third-order perturbation also produces :math:`O(U^4)`
terms that are still missing — so the expansion is left unbalanced.  That
story is plausible and may even be right here, but it is contradicted as a
*general* principle by the scalar comparison in
`Picard iteration, and why it is not simply "higher order"`_, where carrying
an incomplete higher-order fragment makes the answer **better** at every
order.

So what is established is the measurement, not the mechanism: on this circuit,
at this drive, restoring those terms degrades agreement with the transient.
Whether that is the inconsistency, the harmonic truncation in the spectral
form, or something specific to how HD is referred to a node has not been
determined.

So the truncation in the published method is principled rather than
convenient, and "keep more terms" is the wrong instinct here.  The way to
improve accuracy is to take the *next whole order*, not part of it — and
:ref:`the convergence bound <when-more-terms-stop-helping>` says when even
that will fail.

.. _picard-iteration:

Picard iteration, and why it is not simply "higher order"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Raising the perturbation order needs a way to generate the terms.  The
mechanical route is *Picard iteration*: the circuit equation
:math:`Y x + f(x) = u` rearranges into a fixed point,

.. math::

   x = \mathcal{G}\bigl(u - f(x)\bigr), \qquad \mathcal{G} = Y^{-1}

and one simply substitutes the current guess back in:

.. code-block:: text

    x_0 = G u                    <- the linear solution
    x_1 = G u - G f(x_0)
    x_2 = G u - G f(x_1)
    ...

Expanding the second iterate shows the connection to the perturbation series:

.. math::

   x_2 &= \mathcal{G}u - \mathcal{G}f\bigl(x^{(0)} + x^{(1)}\bigr) \\
       &= \mathcal{G}u - \mathcal{G}\Bigl[f(x^{(0)})
          + f'(x^{(0)})x^{(1)}
          + \tfrac{1}{2}f''(x^{(0)})\bigl(x^{(1)}\bigr)^2 + \dots\Bigr] \\
       &= x^{(0)} + x^{(1)} + x^{(2)} + \bigl(\text{a tail of higher-order terms}\bigr)

The :math:`n`-th iterate therefore contains **all** perturbation terms up to
order :math:`n` — this is the 2005 reference's Theorem 2 — **plus a fragment
of the orders beyond it**.

.. important::

   Picard iteration is **not** the same thing as raising the perturbation
   order, and the two must not be conflated.  They agree at :math:`n = 0` and
   :math:`n = 1` and differ from :math:`n = 2` onward.  For the scalar problem
   :math:`Yx + bx^2 = u` the difference at the second iterate is exactly
   :math:`-\mathcal{G}bx_1^2`, which is one of the two terms of
   :math:`x^{(3)} = -\mathcal{G}b\bigl(2x_0x_2 + x_1^2\bigr)` — a fragment
   of the next order, not the whole of it.

   Genuinely raising the perturbation order means constructing
   :math:`x^{(3)}, x^{(4)}, \dots` explicitly from the composition formula
   (Faà di Bruno, which the 1997 reference uses for exactly this purpose).
   **That is not implemented here.**  The ``order`` parameter of
   :func:`~pycircuit.circuit.distortion.harmonic_response_spectral` counts
   Picard iterations, and its name is a convenience rather than a claim.

Which of the two converges faster is **not settled**, and the little evidence
available points against the intuition that a clean truncation must win.  On
that same scalar problem, comparing successive partial sums of the true
perturbation series against successive Picard iterates:

.. code-block:: text

    n    true perturbation      Picard
    0        6.594e-02        6.594e-02
    1        2.156e-02        2.156e-02
    2        9.067e-03        6.387e-03
    3        4.331e-03        1.955e-03
    4        2.234e-03        5.927e-04

Both are monotone, and the *inconsistent* one — Picard, carrying its partial
higher-order fragment — is the more accurate at every order past the first.

Convergence is the familiar condition.  Picard converges precisely when the
map is a contraction, :math:`|\mathcal{G}f'(x)| < 1` — the same
:math:`a < \ln 2` threshold derived earlier.  The two regimes are worth
separating clearly:

**Below the bound**, iterating *to convergence* stops approximating anything.
The fixed point of the iteration is the exact periodic steady state, so this
becomes harmonic balance solved by successive substitution — better than any
finite perturbation order, because it is not a truncation at all.

**Above the bound**, it diverges.  A large change between one iterate and the
next there is not refinement; it is the iteration failing to settle, and
running more passes makes it worse.

This is why
:func:`~pycircuit.circuit.distortion.harmonic_response_spectral` takes
``order`` as an experimental knob rather than as an improvement over
:func:`~pycircuit.circuit.distortion.harmonic_response`.  It exists to
measure what raising the order and the harmonic cutoff actually buy —
which, as the tables above show, is much less than one would expect, and
below the bound is very nearly nothing, because second order has already
converged.

Raising the order: the answer is "yes, but only all the way"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Putting the two knobs together settles the question of whether more
perturbation orders help.  Below the convergence bound, iterating **to
convergence** does better than any truncation — but stopping partway can be
worse than not having started.  Regenerated on every build:

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from pycircuit.circuit import numeric
    from pycircuit.circuit.distortion import (harmonic_response,
        harmonic_response_spectral as spectral, ExponentialNonlinearity,
        Harmonic)

    IS_D = 1e-13
    VT = numeric.kboltzmann*300/numeric.qelectron
    Ib, Rl, Cl, f0 = 1e-3, 1e3, 1e-9, 1e4
    v0 = brentq(lambda v: v/Rl + IS_D*np.expm1(v/VT) - Ib, 0, 1)
    ISe = IS_D*np.exp(v0/VT); al = ISe/VT
    w0 = 2*np.pi*f0
    Y = lambda s: 1/Rl + s*Cl + al
    aG = lambda h, r: [np.asarray(r, dtype=complex)[0]/Y(1j*h.frequency((w0,)))]
    fe = lambda v: ISe*(np.expm1(v/VT) - v/VT)
    fpe = lambda v: (ISe/VT)*np.expm1(v/VT)
    exact = ExponentialNonlinearity(ISe, VT, port=0, n=1)

    drive = 1e-4                        # a ~ 0.22, well inside the bound
    a = abs(aG(Harmonic((1,)), [drive])[0])/VT

    rows = [('published, order 2',
             harmonic_response(aG, [drive], exact, tones=(w0,),
                               toolkit=numeric).HD3(0))]
    for order in (2, 5, 12):
        rows.append(('Picard, order %d' % order,
                     spectral(aG, [drive], fe, fpe, (w0,), numeric,
                              n_harmonics=6, order=order).HD3(0)))

    w = [22, 14]
    sep = ' '.join('='*x for x in w)
    print('Drive a = X1/VT = %.3f, harmonic cutoff 6.' % a)
    print('')
    print(sep)
    print(' '.join('%-*s' % (x, h) for x, h in zip(w, ['method', 'HD3'])))
    print(sep)
    for name, hd3 in rows:
        print(' '.join('%-*s' % (x, c) for x, c in
                       zip(w, [name, '%.4e' % hd3])))
    print(sep)

Against the transient measurement of this circuit the published second-order
truncation lands within a fraction of a percent, the Picard iterates at order
2 and 5 are several times further away, and by order 12 the iteration has
converged and is closer than any of them.

The shape of that needs a caveat, and it is an important one.  These are
*Picard iterates*, not perturbation-order truncations — see the box in
`Picard iteration, and why it is not simply "higher order"`_.  The
non-monotonicity visible here (order 5 worse than order 2) is therefore a
property of **this iteration under a harmonic cutoff**, and specifically not
evidence that the perturbation series behaves that way.  On a scalar problem
where the true series can be written down term by term, it converges
monotonically.

Whether a genuine third- or fifth-order perturbation truncation would improve
these numbers is **unanswered**, because it has not been built.

The harmonic cutoff meanwhile **saturates**.  Raising it from 6 to 12 to 24
changes nothing at all — every digit is identical.  Six harmonics is simply
enough to represent this waveform, and no further Fourier resolution is
available to be had.

So the honest summary for anyone asking whether to raise the order:

* **Below the bound, iterate to convergence or do not bother.**  Converged,
  the iteration is exact and beats every truncation; stopped partway it can
  be worse than the cheap second-order form.
* **Raising the harmonic cutoff is not the lever.**  It saturates quickly and
  cheaply.
* **Above the bound none of it applies** — the iteration diverges, and at
  large drive it overflows outright rather than returning a poor answer,
  which is at least an honest failure.
* **The published second-order form remains the right default.**  It is
  consistent at :math:`U^3`, costs three linear solves, and is within a
  fraction of a percent wherever the method is valid at all.
