Limits of the perturbation method
=================================

The method in :doc:`distortion` is a second-order perturbation truncated at
the third power of the drive.  Both of those are choices, and the obvious
question is what relaxing them would buy.

This page answers that with measurements rather than argument.  It is separate
from the theory page deliberately: none of it is needed to *use* the analysis,
and much of it documents experiments whose conclusion was "this does not
help".  Those are worth recording — the alternatives are the ones a reader
would otherwise try — but they are not the method.

.. note::

   The apparatus used here,
   :func:`~pycircuit.circuit.distortion.harmonic_response_spectral`, is a
   **measuring instrument, not an analysis path**.  It is FFT-based, so it
   cannot be used symbolically at all, and nothing in the production path
   calls it.  See :doc:`distortion` under "Why this method suits a symbolic
   tool" for why that matters.

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
