# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Harmonic distortion of weakly nonlinear driven circuits, by perturbation.

A weakly nonlinear circuit is analysed as a *perturbation of its linearised
model*: the linear solution is found first, the nonlinearity is evaluated on
it, and the resulting harmonics are pushed back through the same linear
circuit.  Repeating that order by order gives the harmonics, and their ratios
to the fundamental are the distortion figures.

The method is Buonomo & Lo Schiavo's; the primary reference is *Predicting
nonlinear distortion in multistage amplifiers and gm-C filters*, Analog
Integrated Circuits and Signal Processing 77:483-493 (2013), whose equation
numbers are cited throughout.  ``doc/distortion_plan.md`` has the staged plan
and its gates; ``doc/src/circuit/distortion.rst`` has the theory.

**The property the whole design rests on:** every order is a solve against
the *same* linear operator ``Y(s) = G + sC``, evaluated at a different
frequency.  Nothing here needs a nonlinear solve, a new matrix structure, or
a change to ``circuit.py`` -- ``Y`` is what AC analysis already assembles.

Scope of this module today (stages 1-4 of the plan): any number of nonlinear
devices; single-tone harmonic distortion with either cubic-polynomial or
exponential (diode/bipolar) nonlinearities, and two-tone intermodulation for
cubic ones.  Two-tone for an exponential device is not implemented because no
reference derives it -- see
:meth:`ExponentialNonlinearity.intermodulation_sources`.

**The harmonics returned are node quantities.**  Distortion at the
nonlinearity's controlling node is *not* distortion at the circuit output:
mapping to an output generally carries both a feedforward path and a
frequency-dependent factor.  For the bipolar example in ``test_distortion``
the two differ by a factor of ten -- a plausible-looking number, not an
error.  Referring the harmonics to whatever the caller means by "output" is
the caller's job.
"""

from . import circuit as circuit_module
from .circuit import gnd, defaultepar


class Harmonic(tuple):
    """A frequency component, as integer multiples of the input tones.

    ``Harmonic((2,))`` is the second harmonic of a single tone;
    ``Harmonic((2, -1))`` is the ``2*w1 - w2`` intermodulation product of two.
    Kept as a tuple rather than an int from the outset because the whole
    two-tone extension is "replace the scalar index by a pair" (2013 paper
    eq. 21 onward) -- hardcoding a single harmonic number here would turn
    that extension into a rewrite.
    """

    __slots__ = ()

    def frequency(self, tones):
        """Angular frequency of this component, given the tone frequencies."""
        if len(tones) != len(self):
            raise ValueError('%d tone(s) given for a %d-tone index %r'
                             % (len(tones), len(self), tuple(self)))
        return sum(m * w for m, w in zip(self, tones))

    def order(self):
        """Total order -- ``|m| + |n| + ...``, the power of input it scales as."""
        return sum(abs(m) for m in self)

    def __repr__(self):
        return 'Harmonic(%r)' % (tuple(self),)


def taylor_coefficients(element, x0, epar=defaultepar, toolkit=None):
    """Quadratic and cubic Taylor coefficients of ``element`` about ``x0``.

    Returns ``(b, c)`` such that the element's *strictly nonlinear* current is
    ``b*v**2 + c*v**3 + O(v**4)`` in the deviation ``v`` from the operating
    point, i.e. ``b = i''(x0)/2`` and ``c = i'''(x0)/6``.

    These are **derived from the element's own pure model function**, never
    declared separately.  An element that stated its Taylor coefficients
    directly would be stating its model twice, and the two statements could
    silently disagree -- exactly the defect ``VCVS_limited`` carried for years
    (``doc/architecture.md`` P9) and that ``test_element_jacobians`` exists to
    catch.  Deriving them means they cannot drift, because they *are* the
    model.

    Args:
        element: A circuit element exposing ``eval_i_pure``.
        x0: Operating-point node voltages for this element's terminals.
        epar: Environment parameters.
        toolkit: Override the element's toolkit (for testing).

    Returns:
        ``(b, c)`` scalars, in the element's port variable.
    """
    toolkit = toolkit or element.toolkit

    if not hasattr(element, 'eval_i_pure'):
        raise TypeError(
            '%s has no eval_i_pure, so its Taylor coefficients cannot be '
            'derived from its model. Elements whose nonlinearity is not '
            'expressed as a pure function are not supported.'
            % type(element).__name__)

    ## Build the parameter dict from the element's declared instparams.
    ## Passing {} instead would *silently* fall back to the defaults baked
    ## into eval_i_pure -- the coefficients would come out plausible and
    ## wrong, which is the failure mode this whole module is built to avoid.
    if hasattr(element, '_model_params'):
        params = element._model_params()
    else:
        params = {p.name: getattr(element.iparv, p.name)
                  for p in element.instparams}

    ## The port variable is the differential voltage across the element's
    ## first two terminals; scan along it and read the current out of the
    ## first terminal.  Stage 1 is a single two-terminal nonlinearity; the
    ## multi-port generalisation (2013 eq. 6's matrix form) is stage 2.
    def port_current(v):
        x = list(x0)
        x[0] = x0[0] + v
        return element.eval_i_pure(toolkit.array(x), params, epar, toolkit)[0]

    d2 = toolkit.nth_derivative(port_current, 0, 2)
    d3 = toolkit.nth_derivative(port_current, 0, 3)
    return d2 / 2, d3 / 6


def scalar_nonlinearity(port, b, c, n, toolkit):
    """Coefficient matrices for one nonlinearity at ``port`` in an ``n``-node circuit.

    Convenience for the common single-device case; the matrices it builds are
    the same objects :func:`harmonic_response` takes in general.
    """
    B = toolkit.zeros((n, n))
    C = toolkit.zeros((n, n))
    B[port, port] = b
    C[port, port] = c
    return B, C


class CubicNonlinearity:
    """A polynomial nonlinearity ``f_s = sum_k b_sk x_k**2 + c_sk x_k**3``.

    The matrices follow the 2013 paper's eq. (6): entry ``[s, k]`` is the
    coefficient of the current injected at node ``s`` and controlled by node
    ``k``.  Note the restriction that comes with that form -- **self-terms
    only**.  A device whose current depends on the product of two *different*
    controlling voltages is outside this formulation.
    """

    def __init__(self, b, c):
        self.b = b
        self.c = c

    def harmonic_sources(self, X1, toolkit):
        """Fourier coefficients of the nonlinearity driven by ``X1``.

        For ``x = Re[V e^{j0}]``, ``x**2`` has second-harmonic coefficient
        ``V**2/2`` and ``x**3`` has third-harmonic coefficient ``V**3/4``;
        the derivative ``df/dx = 2b x + 3c x**2`` has fundamental coefficient
        ``2 b V``.  Those three are all the recurrence needs.
        """
        dot = toolkit.dot
        F2 = _scale(dot(self.b, _elementwise_pow(X1, 2)), _reciprocal(toolkit, 2))
        F3 = _scale(dot(self.c, _elementwise_pow(X1, 3)), _reciprocal(toolkit, 4))
        ## B1[s,k] = 2 b[s,k] X1[k]; folded into the mixing term below as a
        ## plain elementwise product, which is what it reduces to.
        B1_apply = lambda vec: _scale(
            dot(self.b, _elementwise_mul(X1, vec)), 2)
        return F2, F3, B1_apply

    def intermodulation_sources(self, X10, X01, toolkit):
        """Fourier coefficients for the ``2*w1 - w2`` product (2013 eqs. 24, 29).

        Note which coefficient reaches which term.  Only the **cubic**
        coefficient produces ``2*w1 - w2`` directly: squaring two tones gives
        components at ``2w1``, ``2w2``, ``w1 +/- w2`` and DC, never at
        ``2w1 - w2``.  The quadratic coefficient reaches the intermodulation
        product only at second order, by mixing the ``2w1`` and ``w1 - w2``
        products it *does* create back against a tone.  That is why an
        analysis truncated at first order can badly misjudge IM3 in a circuit
        whose dominant nonlinearity is quadratic -- and, in the worked example,
        why the two contributions very nearly cancel.
        """
        dot = toolkit.dot
        conj = toolkit.conj
        X01c = [conj(v) for v in X01]

        ## eq. 24 -- cubic only.
        F_2m1 = _scale(dot(self.c, _elementwise_mul(X01c,
                                                    _elementwise_pow(X10, 2))),
                       3 * _reciprocal(toolkit, 4))
        ## eq. 29 -- the quadratic mixing products, both from b.
        F_20 = _scale(dot(self.b, _elementwise_pow(X10, 2)),
                      _reciprocal(toolkit, 2))
        F_1m1 = dot(self.b, _elementwise_mul(X10, X01c))

        def mix(u, v):
            """``2 b diag(u) v`` -- the derivative coefficient applied to v."""
            return _scale(dot(self.b, _elementwise_mul(u, v)), 2)

        return F_2m1, F_20, F_1m1, mix, X01c


class PolynomialNonlinearity:
    """A single-port polynomial nonlinearity of arbitrary order.

    ``f(v) = c[0]*v**2 + c[1]*v**3 + c[2]*v**4 + ...`` -- coefficients start at
    the quadratic term, since a strictly nonlinear ``f`` has no linear part.

    Orders above three matter even though the recurrence keeps only harmonics
    up to the third: ``v**4`` contributes to the *second* harmonic and
    ``v**5`` to the third.  Raising the polynomial order is therefore a
    different axis from raising the harmonic order -- it sharpens the
    harmonics already being computed rather than adding new ones.

    For a drive ``x = Re[X e^{jt}]`` with ``A = |X|`` the standard cosine
    power expansions give

    .. code-block:: text

        F2 = (X**2/2)(c2 + c4 A**2 + ...)
        F3 = (X**3/4)(c3 + (5/4) c5 A**2 + ...)
        B1 =  X      (2 c2 + 3 c4 A**2 + ...)

    which this implements up to fifth order.  Beyond that it raises rather
    than silently dropping terms.
    """

    def __init__(self, coefficients):
        self.c = list(coefficients)
        if len(self.c) > 4:
            raise NotImplementedError(
                'PolynomialNonlinearity is derived to fifth order; a %d-order '
                'term would also feed the harmonics kept here and cannot be '
                'ignored silently.' % (len(self.c) + 1))

    def _padded(self):
        c = self.c + [0] * (4 - len(self.c))
        return c[0], c[1], c[2], c[3]        # c2, c3, c4, c5

    def harmonic_sources(self, X1, toolkit):
        c2, c3, c4, c5 = self._padded()
        X = X1[0]
        A2 = abs(X) ** 2
        F2 = [(X ** 2 / 2) * (c2 + c4 * A2)]
        F3 = [(X ** 3 / 4) * (c3 + 1.25 * c5 * A2)]
        b1 = X * (2 * c2 + 3 * c4 * A2)
        return F2, F3, (lambda vec: [b1 * vec[0]])

    def intermodulation_sources(self, X10, X01, toolkit):
        raise NotImplementedError(
            'Two-tone for a general polynomial is not derived here; use '
            'CubicNonlinearity, whose two-tone case the sources cover.')

    @classmethod
    def taylor_of_exponential(cls, I_S, V_T, order=5):
        """Coefficients of ``I_S(e^{v/V_T} - 1 - v/V_T)`` up to ``order``.

        Handy for asking how much of a distortion error is the *device model*
        being truncated rather than the perturbation series: compare this
        against :class:`ExponentialNonlinearity`, which has no such
        truncation.
        """
        import math
        return cls([I_S / (math.factorial(n) * V_T ** n)
                    for n in range(2, order + 1)])


class ExponentialNonlinearity:
    r"""A junction nonlinearity ``f(v) = I_S (e^{v/V_T} - 1 - v/V_T)``.

    This is the diode and bipolar shape, and it is where the cubic
    approximation stops being adequate: the 2005 paper measures a
    cubic-truncated exponential converging to a second harmonic roughly 20%
    away from the true value.  Handled here **exactly at every harmonic
    order**, because the Fourier coefficients of ``exp(a cos t)`` are modified
    Bessel functions:

    .. math::

        F_{m} = 2 I_S\, I_m(|X_1| / V_T), \quad m > 1

    (2005 paper eqs. 46-47).  No truncation, no polynomial fit.

    The linear and constant terms are subtracted so that ``f(0) = 0`` and
    ``f'(0) = 0`` -- the formulation requires the nonlinearity to be strictly
    nonlinear, with any linear part already absorbed into the circuit's
    admittance as the junction conductance ``I_S / V_T``.

    .. note::

       The Bessel form assumes the drive is a *real* cosine, i.e. that the
       phase of ``X1`` has been absorbed into the time origin.  That is the
       convention the source uses (it chooses the input phase to make the
       controlling voltage real), and it is exact for a single tone.
    """

    def __init__(self, I_S, V_T, port, n):
        self.I_S = I_S
        self.V_T = V_T
        self.port = port
        self.n = n

    def _bessel(self, order, arg, toolkit):
        try:
            return toolkit.besseli(order, arg)
        except AttributeError:
            from scipy.special import iv
            return iv(order, arg)

    def harmonic_sources(self, X1, toolkit):
        drive = X1[self.port]
        amplitude = abs(drive) / self.V_T

        ## The Jacobi-Anger expansion is stated for a *real* cosine drive,
        ## exp(a cos t) = sum_m I_m(a) exp(j m t).  A drive X1 = A exp(j phi)
        ## is A cos(t + phi), so harmonic m picks up exp(j m phi).
        ##
        ## Omitting that factor leaves every |F_m| right and every argument
        ## wrong.  For *one tone and one device* that turns out to be
        ## unobservable: the same factor multiplies F_m and the second-order
        ## mixing term, so it comes out as a common multiplier on the whole
        ## harmonic -- exactly a choice of time origin, and it cancels in
        ## every magnitude ratio (pinned by
        ## test_single_tone_magnitudes_are_invariant_to_the_drive_phase).
        ##
        ## It is carried anyway because it stops being free the moment the
        ## phases cannot all be absorbed into one time origin: a second
        ## nonlinear device seeing a different phase, or a second tone.  That
        ## is the real reason two-tone support needs complex amplitudes here,
        ## not any difficulty in the Bessel expansion itself.
        phase = 1 if amplitude == 0 else drive / abs(drive)

        I1 = self._bessel(1, amplitude, toolkit)
        I2 = self._bessel(2, amplitude, toolkit)
        I3 = self._bessel(3, amplitude, toolkit)

        ## Built as lists, not toolkit.zeros(): that returns a *real* array,
        ## and assigning a complex harmonic into it silently discards the
        ## imaginary part -- the same class of defect as the missing phase.
        F2 = [0] * self.n
        F3 = [0] * self.n
        F2[self.port] = 2 * self.I_S * I2 * phase ** 2
        F3[self.port] = 2 * self.I_S * I3 * phase ** 3

        ## B1 is the fundamental Fourier coefficient of df/dv, which for the
        ## exponential is (I_S/V_T)*exp(v/V_T) -- again a Bessel coefficient.
        b1 = 2 * (self.I_S / self.V_T) * I1 * phase

        def B1_apply(vec):
            out = [0] * self.n
            out[self.port] = b1 * vec[self.port]
            return out

        return F2, F3, B1_apply

    def intermodulation_sources(self, X10, X01, toolkit):
        """Not implemented -- see ``doc/distortion_plan.md`` section 7.

        Not for the reason it might appear.  The mathematics is easy: the
        exponential *factorises* over a sum of tones, so
        ``exp(a1 cos t1 + a2 cos t2)`` expands as
        ``sum_{m,n} I_m(a1) I_n(a2) exp(j(m t1 + n t2))`` -- an ordinary
        product of Bessel functions, verified to machine precision against a
        2-D numerical Fourier transform.  There is no two-argument special
        function involved.

        What is actually missing is phase.  For a single tone the drive can
        be taken real by absorbing its phase into the time origin, which is
        what this class does (it uses ``abs(X1[port])``).  With two tones only
        one phase can be absorbed; the *relative* phase is physical, and the
        first- and second-order contributions to IM3 nearly cancel in
        practice, so the total depends on it.  Supporting this means carrying
        complex amplitudes through the Bessel path rather than magnitudes.

        Refused rather than guessed because no reference in the source set
        derives it and nothing would currently check the result.  A numerical
        two-tone Fourier extraction is a perfectly good oracle -- that is the
        route, not more reading.
        """
        raise NotImplementedError(
            'Two-tone intermodulation for an exponential device is not '
            'implemented: the Bessel path here carries magnitudes only, and '
            'two tones need relative phase. See doc/distortion_plan.md '
            'section 7.')


class CompositeNonlinearity:
    """Several nonlinearities in one circuit, summed.

    A real circuit mixes device types -- a junction here, a transconductor
    there -- and their harmonic contributions simply add, each injected at its
    own node.  This composes any number of nonlinearity objects into one.

    It is also what makes a *phase* error in any single member observable.
    For one tone and one device the factor ``exp(j m phi)`` is common to every
    term and cancels out of all magnitudes (it is just a choice of time
    origin).  With two devices whose controlling nodes sit at different
    phases, no single time origin makes both real, so their relative phase is
    physical and a member that drops it produces a wrong sum.  See
    ``test_distortion``.
    """

    def __init__(self, *parts):
        if not parts:
            raise ValueError('CompositeNonlinearity needs at least one part')
        self.parts = parts

    def harmonic_sources(self, X1, toolkit):
        F2_total = F3_total = None
        appliers = []
        for part in self.parts:
            F2, F3, B1_apply = part.harmonic_sources(X1, toolkit)
            F2_total = F2 if F2_total is None else _add(F2_total, F2)
            F3_total = F3 if F3_total is None else _add(F3_total, F3)
            appliers.append(B1_apply)

        def B1_apply_all(vec):
            total = None
            for apply in appliers:
                contribution = apply(vec)
                total = (contribution if total is None
                         else _add(total, contribution))
            return total

        return F2_total, F3_total, B1_apply_all

    def intermodulation_sources(self, X10, X01, toolkit):
        raise NotImplementedError(
            'Two-tone for a composite nonlinearity is not implemented; the '
            'members would have to agree on the phase convention first. See '
            'doc/distortion_plan.md section 7.')


def harmonic_response(apply_G, U, nonlinearity, tones, toolkit):
    """Harmonics of a weakly nonlinear circuit, to second perturbation order.

    Implements the 2013 paper's eqs. (18)-(20), written in the form that holds
    for *any* analytic nonlinearity rather than only a polynomial one:

    .. code-block:: text

        X1 =  G(w)  U
        X2 = -G(2w) F2
        X3 = -G(3w) [ F3 + (1/2) B1 . X2 ]

    where ``F2``, ``F3`` are the second- and third-harmonic Fourier
    coefficients of the nonlinearity evaluated on the *known* linear solution,
    and ``B1`` is the fundamental coefficient of its derivative.  A
    :class:`CubicNonlinearity` supplies those from polynomial coefficients; an
    :class:`ExponentialNonlinearity` supplies them in closed form from
    modified Bessel functions.  The recurrence does not care which.

    Every line is a solve against the *same* linear operator at a different
    frequency -- ``apply_G`` is called at the fundamental, at twice it and at
    three times it, and nothing else about the circuit changes.  That is the
    property that makes this an analysis rather than a solver.

    Args:
        apply_G: ``apply_G(harmonic, rhs) -> vector``, solving ``Y(s) x = rhs``
            with ``s`` the frequency of ``harmonic``.  The caller's single
            point of contact with the circuit.
        U: Source vector driving the fundamental.
        nonlinearity: An object with ``harmonic_sources(X1, toolkit)``
            returning ``(F2, F3, B1_apply)``.
        tones: Angular frequencies of the input tones.
        toolkit: The active toolkit.

    Returns:
        :class:`DistortionSolution` holding the full harmonic *vectors*.
    """
    pad = (0,) * (len(tones) - 1)
    one, two, three = (Harmonic((m,) + pad) for m in (1, 2, 3))

    ## Order 0 (eq. 18): the plain linear response.  This is also the final
    ## answer for the fundamental -- the perturbation never corrects it at
    ## this truncation, so no gain compression is predicted.
    X1 = apply_G(one, U)

    F2, F3, B1_apply = nonlinearity.harmonic_sources(X1, toolkit)

    ## Order 1 (eq. 19).  The nonlinearity is evaluated on a solution that is
    ## already known, so its harmonics follow in closed form.  That Fourier
    ## step is where harmonic balance enters -- and each harmonic is then
    ## pushed back through G at *its own* frequency, not at the fundamental.
    X2 = apply_G(two, _negate(F2))

    ## Order 2 (eq. 20).  Two contributions at the third harmonic: the
    ## nonlinearity's own third-harmonic content, and -- through the second
    ## perturbation order -- its *derivative* mixing the fundamental with the
    ## second harmonic it just produced.  For a cubic that second term carries
    ## no cubic coefficient at all: it is the quadratic nonlinearity acting
    ## twice.  Dropping it is the classic way to under-predict HD3.
    mixing = _scale(B1_apply(X2), _reciprocal(toolkit, 2))
    X3 = apply_G(three, _negate(_add(F3, mixing)))

    return DistortionSolution({one: X1, two: X2, three: X3}, tones, toolkit)


def intermodulation_response(apply_G, U1, U2, nonlinearity, tones, toolkit):
    """The ``2*w1 - w2`` intermodulation product, to second perturbation order.

    Implements the 2013 paper's eqs. (21)-(30).  This is the *same* recurrence
    as :func:`harmonic_response` -- evaluate the nonlinearity on a known
    solution, push each component back through ``G`` at its own frequency --
    with the scalar harmonic index replaced by a pair.  Carrying
    :class:`Harmonic` as a tuple from the first commit is what makes this an
    extension rather than a rewrite.

    Second order is written here as an explicit two-term sum rather than as a
    general convolution: ``2*w1 - w2`` decomposes exactly two ways, as
    ``(2*w1) + (-w2)`` and as ``(w1) + (w1 - w2)``, and those are the two
    products that can mix back up to it.  The 2005 paper writes the same thing
    as a double convolution over all index pairs; enumerating it is equivalent
    once the sum is truncated to terms of third order in the input.

    Args:
        apply_G: ``apply_G(harmonic, rhs) -> vector``.
        U1, U2: Source vectors for the two tones.
        nonlinearity: Must provide ``intermodulation_sources``.
        tones: ``(w1, w2)``.
        toolkit: The active toolkit.

    Returns:
        :class:`DistortionSolution` carrying the fundamentals, the two mixing
        products, and the intermodulation product itself under
        ``Harmonic((2, -1))``.  The first- and second-order contributions are
        kept separately under ``('first',)`` and ``('second',)`` -- in the
        reference circuit they very nearly cancel, so the split is what makes
        a disagreement diagnosable rather than merely visible.
    """
    if len(tones) != 2:
        raise ValueError('intermodulation_response needs exactly two tones')

    h10, h01 = Harmonic((1, 0)), Harmonic((0, 1))
    h2m1, h20, h1m1 = Harmonic((2, -1)), Harmonic((2, 0)), Harmonic((1, -1))

    ## eq. 21 -- each tone's linear response, at its own frequency.
    X10 = apply_G(h10, U1)
    X01 = apply_G(h01, U2)

    F_2m1, F_20, F_1m1, mix, X01c = nonlinearity.intermodulation_sources(
        X10, X01, toolkit)

    ## eq. 25 -- first order: the cubic coefficient reaching 2w1-w2 directly.
    first = apply_G(h2m1, _negate(F_2m1))

    ## eq. 29 -- the quadratic products that exist at 2w1 and w1-w2 ...
    X20 = apply_G(h20, _negate(F_20))
    X1m1 = apply_G(h1m1, _negate(F_1m1))

    ## eq. 28 -- ... mixed back up to 2w1-w2 against a tone.
    second = apply_G(h2m1, _negate(_scale(
        _add(mix(X01c, X20), mix(X10, X1m1)), _reciprocal(toolkit, 2))))

    return DistortionSolution({
        h10: X10, h01: X01, h20: X20, h1m1: X1m1,
        h2m1: _add(first, second),
        Harmonic(('first',)): first,
        Harmonic(('second',)): second,
    }, tones, toolkit)


def harmonic_response_spectral(apply_G, U, f, fprime, tones, toolkit,
                               n_harmonics=3, n_samples=1024, order=2):
    """Perturbation with explicit harmonic cutoff and perturbation order.

    :func:`harmonic_response` implements the published recurrence: second
    perturbation order, truncated at the third power of the drive, keeping
    only the one second-order convolution term that reaches ``U**3``.  This
    variant makes both truncations adjustable so their effects can be
    measured separately.

    ``n_harmonics`` is the Fourier cutoff.  ``order`` is the number of
    perturbation steps.

    **How ``order`` is realised.**  ``order=2`` reproduces the two-step
    structure explicitly (nonlinearity on the linear solution, then derivative
    times that correction).  Higher orders are taken by *Picard iteration* on
    ``x = G(u - f(x))``, which the 2005 reference proves agrees with the
    perturbation series term by term at each order (its Theorems 1 and 2).
    Each additional iteration therefore adds a perturbation order.

    **What raising it does and does not buy.**  The iteration converges only
    while ``|G f'(x)| < 1``; above that threshold it diverges and additional
    orders make matters worse rather than better.  And note that neither this
    nor ``n_harmonics`` produces a *consistent* truncation in the drive
    amplitude the way the published second-order form does -- see
    ``doc/src/circuit/distortion.rst``, "What the truncation is actually
    truncating".  This function exists to measure those effects, not to
    supersede :func:`harmonic_response`.

    Returns:
        :class:`DistortionSolution` carrying harmonics 1..``n_harmonics``.
    """
    import numpy as np

    if len(tones) != 1:
        raise ValueError('harmonic_response_spectral is single-tone')
    if order < 1:
        raise ValueError('order must be at least 1')

    X1 = np.asarray(apply_G(Harmonic((1,)), U), dtype=complex)
    n_nodes = len(X1)
    theta = 2 * np.pi * np.arange(n_samples) / n_samples

    def waveform(spectrum):
        out = np.zeros((n_nodes, n_samples))
        for m, vec in spectrum.items():
            for i in range(n_nodes):
                out[i] += np.real(vec[i] * np.exp(1j * m * theta))
        return out

    def coefficients(signal, m):
        c = 2 * (signal * np.exp(-1j * m * theta)).mean(axis=1)
        return c / 2 if m == 0 else c

    def solve_all(F_of_m):
        """Push each harmonic of a drive back through G at its own frequency."""
        return {m: np.asarray(apply_G(Harmonic((m,)), -F_of_m(m)), dtype=complex)
                for m in range(0, n_harmonics + 1)}

    ## Order 1: the nonlinearity on the linear solution.
    x0 = waveform({1: X1})
    f0 = f(x0)
    first = solve_all(lambda m: coefficients(f0, m))

    if order == 1:
        total = {Harmonic((m,)): (X1 if m == 1 else first[m])
                 for m in range(1, n_harmonics + 1)}
        return DistortionSolution(total, tones, toolkit)

    ## Order 2: derivative times the whole first-order correction.
    x1 = waveform({m: v for m, v in first.items() if m > 0})
    x1 += np.real(first[0])[:, None]
    f1 = fprime(x0) * x1
    second = solve_all(lambda m: coefficients(f1, m))

    total = {}
    for m in range(1, n_harmonics + 1):
        base = X1 if m == 1 else first[m]
        total[Harmonic((m,))] = base + (0 if m == 1 else second[m])

    if order == 2:
        return DistortionSolution(total, tones, toolkit)

    ## Orders beyond the second: continue by Picard iteration on the full
    ## equation, each pass adding one perturbation order.
    spectrum = {0: np.zeros(n_nodes, dtype=complex)}
    spectrum.update({m: total[Harmonic((m,))] for m in
                     range(1, n_harmonics + 1)})
    Uvec = np.asarray(U, dtype=complex)

    for _ in range(order - 2):
        xt = waveform({m: v for m, v in spectrum.items() if m > 0})
        xt += np.real(spectrum[0])[:, None]
        ft = f(xt)
        drive = {m: coefficients(ft, m) for m in range(0, n_harmonics + 1)}
        drive[1] = drive[1] - Uvec
        spectrum = {m: np.asarray(apply_G(Harmonic((m,)), -drive[m]),
                                  dtype=complex)
                    for m in range(0, n_harmonics + 1)}

    return DistortionSolution(
        {Harmonic((m,)): spectrum[m] for m in range(1, n_harmonics + 1)},
        tones, toolkit)


def _negate(vec):
    return [-v for v in vec]


def _add(a, b):
    return [x + y for x, y in zip(a, b)]


def _scale(vec, factor):
    return [v * factor for v in vec]


def _reciprocal(toolkit, n):
    """``1/n`` in the toolkit's own number type.

    Exact (a sympy Rational) under a symbolic toolkit, so the stage-1 gate can
    compare expressions rather than floats; a plain float otherwise.
    """
    try:
        return toolkit.integer(1) / n
    except AttributeError:
        return 1.0 / n


def _elementwise_pow(vec, power):
    return [v ** power for v in vec]


def _elementwise_mul(a, b):
    return [x * y for x, y in zip(a, b)]


class DistortionSolution:
    """Harmonics of one node, and the distortion figures derived from them.

    Holds the complex amplitude of each :class:`Harmonic`.  ``HD2``/``HD3`` are
    ratios to the fundamental, which -- note -- is the *linear* response: the
    method gives the fundamental no nonlinear correction at this truncation
    (2013 eq. 18 carries no perturbation term), so gain compression is not
    predicted.  All four papers in this line share that property.
    """

    def __init__(self, harmonics, tones, toolkit):
        self.harmonics = dict(harmonics)
        self.tones = tuple(tones)
        self.toolkit = toolkit

    def __getitem__(self, index):
        """The whole node vector for one frequency component."""
        return self.harmonics[Harmonic(index)]

    def amplitude(self, index, node):
        """Complex amplitude of component ``index`` at ``node``."""
        return self.harmonics[Harmonic(index)][node]

    def fundamental(self, node):
        pad = (0,) * (len(self.tones) - 1)
        return self.amplitude((1,) + pad, node)

    def ratio(self, index, node):
        """Amplitude of component ``index`` at ``node``, relative to the fundamental.

        The denominator is the *linear* fundamental -- see the class docstring.
        """
        return abs(self.amplitude(index, node)) / abs(self.fundamental(node))

    def HD2(self, node):
        """Second-harmonic distortion at ``node``."""
        return self.ratio((2,), node)

    def HD3(self, node):
        """Third-harmonic distortion at ``node``."""
        return self.ratio((3,), node)

    def perturbation_ratio(self, node):
        """How hard the perturbation is working -- a validity diagnostic.

        The method assumes the nonlinear term is small beside the linear
        drive.  **No paper in this line quantifies "small"**, and only the
        2005 one acknowledges that the convergence radius is hard to predict
        at all.  Without something reported alongside the answer, the failure
        mode is silent: a circuit outside the valid range returns a confident,
        plausible, wrong number.

        This returns the largest correction relative to the fundamental,
        ``max_m>=2 |X_m| / |X_1|``.  It is not a proof of validity -- there
        isn't one available -- but it is measured on the same solve and it
        tracks the true error closely.  Against a transient simulation of a
        biased diode:

        ==========  ==============  ==============
        ratio       HD2 error       HD3 error
        ==========  ==============  ==============
        0.01        0.04%           0.07%
        0.05        1.0%            1.9%
        0.21        21%             44%
        1.04        51%             93%
        ==========  ==============  ==============

        So the relative error is of the same order as this ratio.  Treat
        anything above a few percent as a warning and anything approaching 1
        as meaningless -- at that point the "correction" is the size of the
        thing it corrects, and the series is not converging.

        **Above the convergence bound, more terms do not help.**  The
        perturbation series is the fixed-point iteration
        ``x <- G(u - f(x))``, whose terms shrink only while ``|G f'(x)| < 1``.
        For an exponential device on a junction-dominated node that factor is
        ``exp(a) - 1`` with ``a = |X_1|/V_T``, so it crosses 1 at
        ``a = ln 2 ~ 0.693`` -- a swing of about 17 mV at room temperature.
        Below it the terms fall geometrically and a higher-order
        implementation would pay off; above it they *grow*, and no number of
        additional orders recovers anything.  Two further consequences, both
        measured: raising the *harmonic* order is a separate axis that always
        converges but is subordinate to this bound, and past the bound a more
        exact device model can give a *worse* answer than a crude one, so
        apparent agreement there carries no information.  See
        ``doc/src/circuit/distortion.rst``, "When more terms stop helping".
        """
        fundamental = abs(self.fundamental(node))
        corrections = [abs(vec[node]) for index, vec in self.harmonics.items()
                       if all(isinstance(m, int) for m in index)
                       and sum(abs(m) for m in index) >= 2]
        if not corrections or fundamental == 0:
            return 0.0
        return max(corrections) / fundamental

    def __repr__(self):
        return 'DistortionSolution(%d harmonics, %d tone(s))' % (
            len(self.harmonics), len(self.tones))
