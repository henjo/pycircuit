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

Scope of this module today (stages 1-2 of the plan): cubic polynomial
nonlinearities, single tone, any number of nonlinear devices, plus an
optional nonlinearity acting on the input.  Exponential devices and two-tone
intermodulation are later stages.  The frequency index is deliberately kept
general (a tuple, not an integer) so that two-tone is an extension rather
than a rewrite -- see :class:`Harmonic`.
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


def harmonic_response(apply_G, U, b, c, tones, toolkit,
                      b_h=None, c_h=None, X_in=None):
    """Harmonics of a weakly nonlinear circuit, to second perturbation order.

    Implements the 2013 paper's eqs. (18)-(20) in their general matrix form:
    any number of nonlinear devices, plus an optional nonlinearity acting on
    the *input* itself.

    Every line below is a solve against the *same* linear operator at a
    different frequency -- ``apply_G`` is called at the fundamental, at twice
    it and at three times it, and nothing else about the circuit changes.
    That is the property that makes this an analysis rather than a solver.

    The coefficient matrices follow the paper's eq. (6): entry ``b[s, k]`` is
    the quadratic coefficient of the nonlinear current injected at node ``s``
    and controlled by node ``k``.  Note the model restriction that comes with
    it -- **self-terms only**.  A device whose current depends on the product
    of two *different* controlling voltages is outside this formulation.

    Args:
        apply_G: ``apply_G(harmonic, rhs) -> vector``, solving ``Y(s) x = rhs``
            with ``s`` the frequency of ``harmonic``.  The caller's single
            point of contact with the circuit.
        U: Source vector driving the fundamental (the paper's ``G_h X_in``).
        b, c: ``n x n`` quadratic and cubic coefficient matrices.
        tones: Angular frequencies of the input tones.
        toolkit: The active toolkit.
        b_h, c_h: Optional coefficient matrices for a nonlinearity acting on
            the input vector rather than on the solution.
        X_in: The raw input vector, required if ``b_h``/``c_h`` are given.

    Returns:
        :class:`DistortionSolution` holding the full harmonic *vectors*.
    """
    pad = (0,) * (len(tones) - 1)
    one, two, three = (Harmonic((m,) + pad) for m in (1, 2, 3))

    dot = toolkit.dot

    ## Order 0 (eq. 18): the plain linear response.  This is also the final
    ## answer for the fundamental -- the perturbation never corrects it at
    ## this truncation, so no gain compression is predicted.
    X1 = apply_G(one, U)

    ## The input nonlinearity contributes only here: it does not depend on the
    ## solution, so it never re-enters at a later order.
    def input_drive(coeffs, power):
        if coeffs is None:
            return None
        if X_in is None:
            raise ValueError('b_h/c_h given without X_in')
        return dot(coeffs, _elementwise_pow(X_in, power))

    ## Order 1 (eq. 19).  The nonlinearity is evaluated on the *known* linear
    ## solution, so its harmonics follow in closed form: for x = Re[V e^jt],
    ## the second-harmonic coefficient of x**2 is V**2/2.  That Fourier step
    ## is where harmonic balance enters -- and each harmonic is then pushed
    ## back through G at *its own* frequency, not at the fundamental.
    drive2 = -dot(b, _elementwise_pow(X1, 2))
    h2 = input_drive(b_h, 2)
    if h2 is not None:
        drive2 = drive2 + h2
    ## G is linear, so halving the drive is the same as halving the answer --
    ## and keeps the arithmetic elementwise rather than on the solve result.
    X2 = apply_G(two, _scale(drive2, _reciprocal(toolkit, 2)))

    ## Order 2 (eq. 20).  Two contributions at the third harmonic: the cubic
    ## coefficient acting on the fundamental directly, and -- through the
    ## second perturbation order -- the *quadratic* coefficient mixing the
    ## fundamental with the second harmonic it produced.  Dropping the second
    ## term is the classic way to under-predict HD3 in a circuit whose
    ## dominant nonlinearity is quadratic.
    ##
    ## The paper writes that mixing term as
    ##   -2 b diag(X1) G(j2w) (b_h Xin^2 - b X1^2),
    ## and the inner factor is exactly 2*X2 by eq. (19), so it collapses to
    ## -4 b (X1 elementwise* X2).
    drive3 = -dot(c, _elementwise_pow(X1, 3)) - 4 * dot(b, _elementwise_mul(X1, X2))
    h3 = input_drive(c_h, 3)
    if h3 is not None:
        drive3 = drive3 + h3
    X3 = apply_G(three, _scale(drive3, _reciprocal(toolkit, 4)))

    return DistortionSolution({one: X1, two: X2, three: X3}, tones, toolkit)


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
