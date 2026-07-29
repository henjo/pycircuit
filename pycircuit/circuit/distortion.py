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

Scope of this module today (stages 1-3 of the plan): single tone, any number
of nonlinear devices, with either cubic-polynomial or exponential
(diode/bipolar) nonlinearities.  Two-tone intermodulation is a later stage;
the frequency index is deliberately kept general (a tuple, not an integer) so
that it is an extension rather than a rewrite -- see :class:`Harmonic`.

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
        amplitude = abs(X1[self.port]) / self.V_T
        I2 = self._bessel(2, amplitude, toolkit)
        I3 = self._bessel(3, amplitude, toolkit)
        I1 = self._bessel(1, amplitude, toolkit)

        ## Built as lists, not toolkit.zeros(): that returns a *real* array,
        ## and assigning a complex harmonic into it silently discards the
        ## imaginary part -- a phase error that only shows up off DC.
        F2 = [0] * self.n
        F3 = [0] * self.n
        F2[self.port] = 2 * self.I_S * I2
        F3[self.port] = 2 * self.I_S * I3

        ## B1 is the fundamental Fourier coefficient of df/dv, which for the
        ## exponential is (I_S/V_T)*exp(v/V_T) -- again a Bessel coefficient.
        b1 = 2 * (self.I_S / self.V_T) * I1

        def B1_apply(vec):
            out = [0] * self.n
            out[self.port] = b1 * vec[self.port]
            return out

        return F2, F3, B1_apply


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

    def __repr__(self):
        return 'DistortionSolution(%d harmonics, %d tone(s))' % (
            len(self.harmonics), len(self.tones))
