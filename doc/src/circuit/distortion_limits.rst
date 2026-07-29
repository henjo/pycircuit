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
last.  Approaching :math:`\ln 2` they barely fall.  Above it **the terms
grow** and the series diverges.

.. important::

   **That statement is about the infinite series, and an earlier version of
   this page wrongly extended it to the truncations.**  It said that above the
   bound "no number of additional terms recovers anything".  That is false,
   and measurably so: at :math:`a = 0.88`, where the contraction factor is
   1.34 and the series is genuinely divergent, the error of the *truncated*
   result still falls from 67.7% at :math:`U^3` to 0.86% at :math:`U^{11}`.

   This is ordinary asymptotic-series behaviour.  A divergent series can have
   truncations that improve up to an optimal order and worsen only beyond it,
   and the optimal order gets earlier as the drive rises.  So the practical
   rule is not a hard cutoff:

   * **Below the contraction bound** — the series converges and higher order
     converges with it.
   * **Above it** — higher order still helps, up to a drive-dependent optimum,
     after which it degrades.  Measured by watching successive truncations of
     the third harmonic: at :math:`a = 1.33` (contraction 2.6) the change
     between orders falls to 4.7e-2 and then *rises* to 8.9e-2, which is the
     optimum being passed.  The optimum arrives earlier as the drive rises —
     around :math:`U^{19}` at :math:`a = 0.88` against :math:`U^{17}` at
     :math:`a = 1.33` — exactly as an asymptotic series behaves.

   The bound remains the right thing to *report*, because it marks where
   accuracy stops being guaranteed.  It is not the point where higher order
   becomes useless.

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

.. warning::

   **The section below measures Picard iterations, not perturbation order,
   and its conclusion has been superseded.**  Raising the *genuine*
   perturbation order — building the terms by the composition formula and
   truncating consistently by power of the drive — improves accuracy
   monotonically and substantially: against the transient cross-check on a
   biased diode, HD3 error at a drive ratio of 0.05 falls from 2.85% at
   :math:`U^3` to under 0.01% at :math:`U^5`.  See
   ``doc/distortion_higher_order_plan.md`` (stages C and D) and
   ``test_distortion.py``.

   What follows remains a correct account of what *Picard iteration under a
   harmonic cutoff* does, which is a different and less useful thing.  It is
   kept because the distinction is exactly what an earlier version of this
   work got wrong.

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

Harmonic amplitudes against the transient, order by order
---------------------------------------------------------

The tables below give each harmonic's amplitude at four drive levels and five
truncation orders, with the transient simulation as reference and the
difference from it in parentheses.  Regenerated on every build.

``a = |X1|/V_T`` is the drive measured against the thermal voltage; the
contraction bound sits near ``a = 0.7``.  A dash means the harmonic does not
exist at that truncation: harmonic ``m`` first appears at exactly ``U**m``, so
the fifth harmonic is genuinely absent below ``U**5`` rather than zero.

.. exec-rst::

    import numpy as np, math
    from scipy.optimize import brentq
    from pycircuit.circuit import SubCircuit, gnd, R, C, ISin, numeric
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.distortion import graded_response, Harmonic

    IS_D = 1e-13
    VT = numeric.kboltzmann*300/numeric.qelectron
    Ib, Rl, Cl, f0 = 1e-3, 1e3, 1e-9, 1e4
    cm.default_toolkit = numeric
    v0 = brentq(lambda v: v/Rl + IS_D*np.expm1(v/VT) - Ib, 0, 1)
    ISe = IS_D*np.exp(v0/VT); al = ISe/VT
    w0 = 2*np.pi*f0
    Y = lambda s: 1/Rl + s*Cl + al
    resp = lambda s: 1.0/Y(s)
    aG = lambda h, r: [np.asarray(r, dtype=complex)[0]/Y(1j*h.frequency((w0,)))]
    taylor = [ISe/(math.factorial(n)*VT**n) for n in range(2, 14)]

    ORDERS = (3, 5, 7, 9, 11)
    HARMONICS = (1, 2, 3, 5)
    DRIVES = (1e-4, 2e-4, 4e-4, 6e-4)

    def measure(A, per=16, pts=256, keep=8):
        c = SubCircuit(toolkit=numeric); n = c.add_node('n1')
        c['I'] = ISin(gnd, n, io=Ib, ia=A, freq=f0)
        c['R'] = R(n, gnd, r=Rl); c['C'] = C(n, gnd, c=Cl)
        c['D'] = Diode(n, gnd, IS=IS_D)
        res = Transient(c, toolkit=numeric).solve(
            refnode=gnd, tend=per/f0, timestep=1.0/(f0*pts))
        w = res.v(n); t = np.asarray(w.x[0]); v = np.asarray(w.y)
        g = np.linspace(t[-1]-keep/f0, t[-1], keep*pts, endpoint=False)
        S = np.fft.rfft(np.interp(g, t, v))/(keep*pts)
        return {m: 2*abs(S[m*keep]) for m in HARMONICS}

    ref = {A: measure(A) for A in DRIVES}
    got = {}
    for A in DRIVES:
        for P in ORDERS:
            g = graded_response(resp, A, taylor[:P-1], (w0,), max_power=P)
            for m in HARMONICS:
                got[(A, P, m)] = abs(g.phasor(m))

    for m in HARMONICS:
        head = ['a'] + ['U^%d' % P for P in ORDERS] + ['transient']
        rows = []
        for A in DRIVES:
            a = abs(aG(Harmonic((1,)), [A])[0])/VT
            r = ref[A][m]
            cells = ['%.3f' % a]
            for P in ORDERS:
                if m > P:
                    cells.append('--')
                else:
                    v = got[(A, P, m)]
                    cells.append('%.4g (%+.1f%%)' % (v, (v/r - 1)*100))
            cells.append('%.4g' % r)
            rows.append(cells)

        ## Column widths from the content: a fixed guess overflows the moment
        ## a percentage needs an extra digit, and reST rejects the table
        ## rather than wrapping it.
        widths = [max(len(r[i]) for r in rows + [head])
                  for i in range(len(head))]
        sep = ' '.join('=' * w for w in widths)

        print('')
        print('**Harmonic %d** (volts)' % m)
        print('')
        print(sep)
        print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
        print(sep)
        for cells in rows:
            print(' '.join('%-*s' % (w, c) for w, c in zip(widths, cells)))
        print(sep)

Three things to read out of these.

**Convergence below the bound is excellent and orderly.**  At ``a = 0.22``
every harmonic settles by ``U**5``–``U**7`` to within a few hundredths of a
percent of the simulation, and further orders change nothing.

**Above the bound the truncations overshoot rather than merely lag.**  At
``a = 1.33`` the successive orders sail past the reference — the third
harmonic reaching about +7% and the fifth about +41% at ``U**11`` — which is
the asymptotic series turning around, not slow convergence.  This is visible
here in a way the error-only tables elsewhere on this page obscure, because
the sign of the discrepancy is what distinguishes overshoot from lag.

**The fifth harmonic cannot be validated well at either end.**  At low drive
it is around 1e-7 V, where the transient's own resolution limits the
comparison to a few percent; at high drive the truncation has begun to
diverge.  Its agreement in the middle rows is the meaningful part, and the
outer rows say more about the reference than about the method.

The 1 dB point, order by order
------------------------------

Sweeping the input current amplitude ``i_in`` and watching the fundamental
depart from its small-signal extrapolation gives a single scalar figure of
merit, and so a sharper test of truncation order than any individual harmonic:
each order must predict *where* on the drive axis a threshold is crossed, not
merely an amplitude at a drive chosen for it.

.. note::

   **On this circuit the 1 dB point is an expansion point, not a compression
   point.**  The fundamental grows *faster* than the drive, so the departure
   is :math:`+1` dB.  This is a property of the topology rather than a sign
   convention: the drive is a *current* into a shunt diode, so on the negative
   half-cycle the diode turns off and the node sees the load resistor alone —
   a small-signal conductance of ``1/R`` instead of ``1/R + g_d``, here a
   factor of 17 more transimpedance.  The positive half-cycle compresses, but
   only logarithmically, and the negative half wins.  Both the perturbation
   result and the independent transient simulation agree on the sign.

   A configuration with ``1/R >> g_d`` does not compress either; it simply
   becomes linear, the diode having stopped mattering.  Getting genuine
   compression out of this cell requires a different topology (a series
   element with a clamped load), not a different drive level.

.. exec-rst::

    import numpy as np, math
    from scipy.optimize import brentq
    from pycircuit.circuit import SubCircuit, gnd, R, C, ISin, numeric
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.distortion import graded_response

    IS_D = 1e-13
    VT = numeric.kboltzmann*300/numeric.qelectron
    Ib, Rl, Cl, f0 = 1e-3, 1e3, 1e-9, 1e4
    cm.default_toolkit = numeric
    v0 = brentq(lambda v: v/Rl + IS_D*np.expm1(v/VT) - Ib, 0, 1)
    ISe = IS_D*np.exp(v0/VT); al = ISe/VT
    w0 = 2*np.pi*f0
    Y = lambda s: 1/Rl + s*Cl + al
    resp = lambda s: 1.0/Y(s)
    taylor = [ISe/(math.factorial(n)*VT**n) for n in range(2, 14)]
    g0 = abs(resp(1j*w0))               # small-signal transimpedance, V/A

    ORDERS = (3, 5, 7, 9, 11)
    TARGET = 10**(1/20.)                # +1 dB in amplitude

    def transient_H1(A, per=16, pts=256, keep=8):
        c = SubCircuit(toolkit=numeric); n = c.add_node('n1')
        c['I'] = ISin(gnd, n, io=Ib, ia=A, freq=f0)
        c['R'] = R(n, gnd, r=Rl); c['C'] = C(n, gnd, c=Cl)
        c['D'] = Diode(n, gnd, IS=IS_D)
        res = Transient(c, toolkit=numeric).solve(
            refnode=gnd, tend=per/f0, timestep=1.0/(f0*pts))
        w = res.v(n); t = np.asarray(w.x[0]); v = np.asarray(w.y)
        g = np.linspace(t[-1]-keep/f0, t[-1], keep*pts, endpoint=False)
        S = np.fft.rfft(np.interp(g, t, v))/(keep*pts)
        return 2*abs(S[keep])

    ## Bracket kept tight: each transient evaluation is a full simulation, and
    ## a wide bracket costs build time without buying accuracy.
    A_tr = brentq(lambda A: transient_H1(A)/(g0*A) - TARGET,
                  2.0e-4, 5.0e-4, xtol=1e-8)

    def dev(A, P):
        g = graded_response(resp, A, taylor[:P-1], (w0,), max_power=P)
        return abs(g.phasor(1))/(g0*A) - TARGET

    ## No pipes in generated text: reST reads |X1| as a substitution
    ## reference and errors, and reports it at line 1 of the page because
    ## generated reST carries no source-line provenance.
    head = ['truncation', 'i_in at +1 dB', 'a = X1/VT', 'error vs transient']
    rows = []
    for P in ORDERS:
        Ap = brentq(lambda A: dev(A, P), 1e-5, 3e-3, xtol=1e-10)
        rows.append(['U^%d' % P, '%.2f uA' % (Ap*1e6),
                     '%.4f' % (abs(resp(1j*w0)*Ap)/VT),
                     '%+.2f%%' % ((Ap/A_tr - 1)*100)])
    rows.append(['transient', '%.2f uA' % (A_tr*1e6),
                 '%.4f' % (abs(resp(1j*w0)*A_tr)/VT), '--'])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for cells in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, cells)))
    print(sep)
    print('')
    print('Contraction bound ``a = ln 2 = %.4f``.' % math.log(2))

**The threshold estimate converges monotonically and fast** — the second-order
form misplaces the 1 dB point by about 11%, and each further order removes
roughly a factor of five, reaching the transient's own resolution by
:math:`U^9`.

**The 1 dB point sits essentially at the convergence bound**, ``a = 0.68``
against ``ln 2 = 0.693``.  That near-coincidence is not luck: both are set by
the same ratio of nonlinear to linear response, one asking when the series
stops contracting and the other asking when the fundamental has been visibly
pulled off its linear extrapolation.  It has a practical consequence — this
figure of merit is computable by the method almost exactly up to the point
where the method stops being valid, and not beyond.  A circuit whose 1 dB
point lay well above the bound could not be characterised this way at any
order.

A genuine 1 dB compression point
---------------------------------

The diode cell above *expands*, and neither published circuit in
:doc:`distortion` reaches 1 dB of compression anywhere the series converges.
A cell that does is easy to build once the requirement is stated properly:
the nonlinearity must **load the node** through :math:`Y(s)` — otherwise the
perturbation machinery does no work — and the output must be taken as the
**transconductor current** rather than a node voltage.

.. code-block:: text

    v_in --Rs-- node v --+-- C to ground
                         +-- i(v) = g*(v + alpha*v**3) to ground

with :math:`\alpha = -0.0535\,\mathrm{V}^{-2}` from the 2013 paper.  The
topology is ours rather than a reference's, so what follows **demonstrates**
the method; the validation is the time-domain comparison in the last column.
Regenerated on every build:

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from scipy.integrate import solve_ivp
    from pycircuit.circuit.distortion import graded_response

    g, alpha = 625.2e-6, -0.0535
    Rs, C, f0 = 1.0/g, 1e-12, 1e6
    w0 = 2*np.pi*f0
    resp = lambda s: 1.0/(1/Rs + g + s*C)
    small = g*abs(resp(1j*w0))/Rs

    def node(Vin, n):
        return graded_response(resp, Vin/Rs, [0.0, g*alpha], (w0,),
                               max_power=n)

    def dev_db(Vin, n):
        v = node(Vin, n)
        cube = ((v*v).truncated(n)*v).truncated(n)
        iout = (v + cube.scaled(alpha)).scaled(g)
        return 20*np.log10(abs(iout.phasor(1))/(small*Vin))

    def ode_db(Vin, cycles=120, ramp=50, keep=16):
        ## Ramped drive: a step start overshoots the cubic's turning point
        ## and the model runs away.  See the note below.
        T = 2*np.pi/w0
        def rhs(t, y):
            v = y[0]
            u = Vin*min(1.0, t/(ramp*T))*np.cos(w0*t)
            return [((u - v)/Rs - g*(v + alpha*v**3))/C]
        tend = cycles*T
        ts = np.linspace(tend-keep*T, tend, keep*512, endpoint=False)
        s = solve_ivp(rhs, (0, tend), [0.], t_eval=ts, method='DOP853',
                      rtol=1e-11, atol=1e-16, max_step=T/50)
        i = g*(s.y[0] + alpha*s.y[0]**3)
        S = np.fft.rfft(i)/len(i)
        return 20*np.log10(2*abs(S[keep])/(small*Vin))

    ORDERS = (3, 5, 7, 9, 11)
    head = ['order', 'P1dB (V)', 'shift from previous']
    rows = []
    previous = None
    for n in ORDERS:
        p = brentq(lambda V: dev_db(V, n) + 1.0, 1.0, 5.0, xtol=1e-8)
        rows.append(['U^%d' % n, '%.4f' % p,
                     '--' if previous is None else '%+.4f' % (p - previous)])
        previous = p

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)
    print('')

    p1 = previous
    v1 = abs(node(p1, 11).phasor(1))
    v_turn = np.sqrt(1.0/(3*abs(alpha)))
    print('At ``Vin = %.4f V`` the node fundamental is ``%.4f V`` and the'
          % (p1, v1))
    print('contraction factor is ``%.4f`` against a bound of 1, so the series'
          % abs(3*g*alpha*v1**2*resp(1j*w0)))
    print('is still converging comfortably where the cell compresses by 1 dB.')
    print('')
    print('Against the integration: ``%+.4f dB`` at ``Vin = 2 V`` and'
          % (dev_db(2.0, 11) - ode_db(2.0)))
    print('``%+.4f dB`` at the compression point.'
          % (dev_db(p1, 11) - ode_db(p1)))
    print('')
    print('The cubic turns over at ``%.4f V``, so the node sits at ``%.0f%%``'
          % (v_turn, 100*v1/v_turn))
    print('of the amplitude where the device model stops being physical.')

**A cubic model can only just represent 1 dB compression, and the margin is a
constant of the model.**  A device :math:`i = g(v + \alpha v^3)` with
:math:`\alpha < 0` turns over at :math:`v_{\mathrm{turn}} = 1/\sqrt{3|\alpha|}`;
past that its differential conductance is negative and it is not a physical
device.  Setting the fundamental to :math:`-1` dB gives

.. math::

   \frac{v_{\mathrm{turn}}}{v_{1\mathrm{dB}}}
     = \sqrt{\frac{1/3}{\left(1 - 10^{-1/20}\right)/0.75}}

.. exec-rst::

    import numpy as np

    ## Evaluated here rather than typed: this page has already shipped one
    ## pasted measurement that a later change falsified.
    ratio = np.sqrt((1/3.)/((1 - 10**(-1/20.))/0.75))
    print('which is ``%.6f``, **independent of** ``alpha`` -- so 1 dB always'
          % ratio)
    print('falls at **%.1f%%** of the amplitude where the model breaks down.'
          % (100/ratio))

That is the real reason none of the published circuits reaches 1 dB: not the
perturbation method, and not the choice of circuit, but that **a cubic is a
weakly nonlinear model and 1 dB compression is about where weak nonlinearity
stops describing the device.**

The loaded cell is worse than the isolated device, sitting nearer that limit
than the ratio above suggests, because the same negative :math:`\alpha` that
compresses
the output current also reduces the node's loading — so the node voltage
*expands* toward the breakdown even as the output compresses.  It is also why
the reference integration ramps its drive: started as a step, the transient
overshoots :math:`v_{\mathrm{turn}}` and the cubic diverges to overflow.

A saturating model: ``tanh``
-----------------------------

The cubic's difficulty above is not the perturbation method's — it is that a
cubic stops being a device before it finishes compressing.  A saturating model
has neither problem.
:class:`~pycircuit.circuit.distortion.TanhNonlinearity` implements
:math:`i(v) = I_0 \tanh(v/V_x)`, whose differential conductance
:math:`(I_0/V_x)\,\mathrm{sech}^2(v/V_x)` is positive everywhere: it saturates
instead of turning over, and stays physical at any amplitude.

It needed no new machinery.  :func:`~pycircuit.circuit.distortion.graded_response`
already accepts an arbitrary coefficient list, so a ``tanh`` device is a
supplier of Taylor coefficients rather than a new code path — and those
coefficients are *generated* from :math:`y' = 1 - y^2` rather than transcribed.

.. exec-rst::

    import numpy as np
    from scipy.optimize import brentq
    from scipy.integrate import solve_ivp
    from pycircuit.circuit.distortion import (graded_response,
                                              TanhNonlinearity)

    I_0, V_x = 625.2e-6, 1.0
    g = I_0/V_x
    Rs, Cl, f0 = 1.0/g, 1e-12, 1e6
    w = 2*np.pi*f0
    resp = lambda s: 1.0/(1/Rs + g + s*Cl)
    device = TanhNonlinearity(I_0, V_x)
    series = TanhNonlinearity.series_coefficients(13)
    small = g*abs(resp(1j*w))/Rs

    def dev_db(Vin, power=13):
        v = graded_response(resp, Vin/Rs, device.coefficients(power),
                            (w,), max_power=power)
        iout = v.scaled(g)
        term = v
        for k in range(2, power + 1):
            term = (term*v).truncated(power)
            if series[k]:
                iout = iout + term.scaled(I_0*series[k]/V_x**k)
        return 20*np.log10(abs(iout.phasor(1))/(small*Vin))

    def ode(Vin, cycles=8, keep=4):
        T = 2*np.pi/w
        def rhs(t, y):
            return [((Vin*np.cos(w*t) - y[0])/Rs - I_0*np.tanh(y[0]/V_x))/Cl]
        tend = cycles*T
        ts = np.linspace(tend-keep*T, tend, keep*512, endpoint=False)
        s = solve_ivp(rhs, (0, tend), [0.], t_eval=ts, method='DOP853',
                      rtol=1e-10, atol=1e-16, max_step=T/40)
        S = np.fft.rfft(I_0*np.tanh(s.y[0]/V_x))/len(s.y[0])
        return (20*np.log10(2*abs(S[keep])/(small*Vin)),
                np.max(np.abs(s.y[0]))/V_x)

    head = ['Vin (V)', 'U^5', 'U^9', 'U^13', 'integration', 'peak v/Vx']
    rows = []
    for Vin in (0.5, 1.0, 1.5, 2.0, 2.5):
        ref, peak = ode(Vin)
        rows.append(['%.2f' % Vin] +
                    ['%.4f' % dev_db(Vin, n) for n in (5, 9, 13)] +
                    ['%.4f' % ref, '%.4f' % peak])

    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    sep = ' '.join('=' * w for w in widths)
    print(sep)
    print(' '.join('%-*s' % (w, h) for w, h in zip(widths, head)))
    print(sep)
    for r in rows:
        print(' '.join('%-*s' % (w, c) for w, c in zip(widths, r)))
    print(sep)
    print('')

    p1 = brentq(lambda V: dev_db(V) + 1.0, 0.5, 3.0, xtol=1e-9)
    v1 = abs(graded_response(resp, p1/Rs, device.coefficients(13), (w,),
                             max_power=13).phasor(1))
    print('1 dB compression at ``Vin = %.4f V``, node fundamental'
          ' ``%.4f V`` --' % (p1, v1))
    print('``%.0f%%`` of the Taylor radius ``pi/2 V_x = %.4f V``.'
          % (100*v1/device.radius_of_convergence(),
             device.radius_of_convergence()))

**Every entry compresses; none expands.**  That is the property the cubic
cell could not offer, and it holds by construction rather than by luck.

Both models place 1 dB at roughly two thirds of their limit — the exact
fractions are computed above and in the cubic section — but **the
two limits are not the same kind of thing.**  The cubic's is where the
*device* stops being physical; past it a time-domain solve runs away.
:math:`\tanh`'s is only where the *series* stops converging, and the last row
above is past it: the node peak exceeds :math:`\pi/2` while :math:`U^{13}` is
still within a few thousandths of a dB of the integration, in the usual
asymptotic way.

So a saturating model can be pushed into compression and beyond without ever
leaving the region where it describes a device.  For fitting there are two
independent handles — :math:`I_0` is the saturation level and
:math:`V_x = I_0/g_m` the knee — so a measured small-signal gain and a
measured output limit can both be matched;
:meth:`~pycircuit.circuit.distortion.TanhNonlinearity.from_transconductance`
takes exactly that pair.  A cubic also has two parameters, but its second one
only sets how quickly the model runs out of validity.
