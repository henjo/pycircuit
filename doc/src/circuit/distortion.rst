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

What is missing is phase.  For a single tone the drive can be taken real by
absorbing its phase into the time origin, which is what the implementation
does.  With two tones only one phase can be absorbed; the relative phase is
physical, and since the two contributions to IM3 nearly cancel, the total
depends on it.  Supporting exponential devices with two tones means carrying
complex amplitudes through the Bessel path rather than magnitudes — a real
piece of work, but a bounded one.
