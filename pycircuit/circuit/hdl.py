# -*- coding: utf-8 -*-
# Copyright (c) 2008-2026 Pycircuit Development Team
# See LICENSE for details.
"""Behavioural: a Verilog-A-shaped element DSL compiled to MNA stamps.

An element's behaviour is written as *contribution statements* --
``Contribution(b.I, expr)`` and ``Contribution(b.V, expr)``, the DSL forms
of Verilog-A's ``I(a,b) <+ expr`` and ``V(a,b) <+ expr`` -- inside an
``analog()`` method whose arguments name the terminals.  The metaclass
compiles those equations, symbolically, into everything a pycircuit
element supplies: ``i``/``q``/``u`` vectors, exact ``G``/``C`` Jacobians,
``dudt``, the DC variants for pinned states, ``eval_i_pure``/
``eval_q_pure`` for the JAX batched path, ``state_ic`` seeding and
``periodic_states`` declarations.  See ``hdl.md`` (repository root) for
theory, the capability map against Verilog-A, and usage.

The compilation model (the same set of ideas every equation-defined
device compiler uses -- Qucs EDD, Xyce ABM, ngspice B-sources -- done
here with sympy so the Jacobians are exact and closed-form):

* the element's unknown vector is ``[terminal nodes, internal nodes,
  state nodes, branch currents]``, matching ``Circuit.update_node_map``'s
  layout exactly;
* an ``I``-contribution adds ``+expr``/``-expr`` to the two KCL rows;
* a ``V``-contribution creates a branch-current unknown ``i_b`` (added to
  the KCL rows) and the branch row ``-(v_plus - v_minus) + expr = 0``;
* each additive term routes by kind: ``ddt(arg)`` (optionally scaled by
  an x-independent factor) into the charge vector ``q``, terms free of
  circuit quantities into the source vector ``u(t)`` (they may depend on
  the time symbol ``TIME``), everything else into ``i``;
* each ``idt``/``idtmod`` application becomes a *state*: an internal node
  ``s`` with the extra equation ``ds/dt = arg`` (``q``-row ``s``,
  ``i``-row ``-arg``), the application replaced by ``s`` (``idt``) or by
  the floored wrap of ``s`` (``idtmod``, idtmod.md sec. 1.1);
* ``G = di/dx`` and ``C = dq/dx`` are sympy Jacobians of those vectors,
  compiled with ``lambdify(..., cse=True)``.

Residual convention (transient.py): ``i(x) + dq/dt + u(t) = 0``.

Two public instruments come with it: `explain(element)` prints what the
compiler did -- terminals, the ``x``-vector layout, the parameter values
the generated code will actually read, the symbolic vectors and the
generated source -- and `check_jacobians(element, x)` differentiates
``i`` and ``q`` numerically against ``G`` and ``C`` and scans everything
for non-finite entries.  It has THREE verdicts, not two: an entry the
finite difference cannot resolve at this point -- a kink, a signal below
the value's own quantum, a stiff card -- comes back UNRESOLVED rather
than FAILED.
"""

from pycircuit.circuit import circuit
from pycircuit.circuit.circuit import defaultepar
import pycircuit.utilities.param as param

import sympy
from sympy.codegen.cfunctions import log1p as _log1p
from sympy.core.symbol import Str
import numpy as np

import os
import collections
import inspect
import keyword
import itertools
import re
import types
import warnings


#: Stack of node registries, pushed while an `analog()` body runs.  It
#: exists for ONE check: an internal `Node('out')` in an element that
#: already has an `out` terminal.  Nodes are identified by name
#: everywhere downstream, so the two used to merge in silence -- and the
#: collision cannot be found after the fact, because sympy's expression
#: cache freely hands back an equal atom built during an EARLIER
#: compilation, carrying that compilation's node objects.  Recording the
#: constructions themselves is the only reading that stays true.
_NODE_TRACE = []


class Node(circuit.Node):
    def __init__(self, name=None, isglobal=False):
        super().__init__(name, isglobal)
        if _NODE_TRACE:
            _NODE_TRACE[-1].append(self)

    @property
    def V(self):
        return Quantity('V', self)


class Branch(circuit.Branch):
    """A branch between two terminals, optionally named.

    ``Branch(a, b)`` is the usual anonymous branch.  ``Branch(a, b, 'ch')``
    is Verilog-A's ``branch (a,b) ch;`` -- a DISTINCT branch across the
    same node pair, with its own current unknown when it is
    V-contributed or current-probed.  Two branches differing only in name
    contribute to the same KCL rows but never share an unknown.
    """

    @property
    def V(self):
        return Quantity('V', self)

    @property
    def I(self):
        return Quantity('I', self)

    def __repr__(self):
        return 'Branch(%s, %s%s)' % (
            self.plus.name, self.minus.name,
            '' if self.name is None else ', %r' % self.name)


## NOT a sympy.Symbol subclass, deliberately.  It used to be
## `class Parameter(param.Parameter, sympy.Symbol)`, but sympy interns
## Symbols by name: two Parameters named 'c0' in two different elements
## were THE SAME OBJECT, and whichever was declared last set `default` for
## both -- a silent cross-element parameter leak.  The mixin bought
## nothing, since `generate_code` builds its own `Symbol(name)` for each
## declared parameter and binds values from `iparv` at call time.
Parameter = param.Parameter


class ddt(sympy.Function):
    """Time derivative: routes its argument into the charge vector."""
    nargs = (1,)


class idt(sympy.Function):
    """Time integral -- Verilog-A ``idt(expr[, ic])``.

    Compiles to an extra state unknown ``s`` with equation ``ds/dt =
    expr``; the application evaluates to ``s``.  With ``ic`` given, the
    DC solve pins ``s = ic`` (the LRM: "the initial condition shall force
    the DC solution") and ``uic=True`` seeds it; without, the element has
    a DC solution only if feedback zeroes ``expr`` (idtmod.md sec. 5.1).
    """
    nargs = (1, 2)


class idtmod(sympy.Function):
    """Circular integral -- Verilog-A ``idtmod(expr, ic, modulus, offset)``.

    The state integrates unbounded within a step and evaluates through the
    floored wrap into ``[offset, offset+modulus)``; the compiled element
    declares the state through ``periodic_states()`` so the Phase-2 gauge
    shift keeps it bounded between steps (idtmod.md sec. 5.2), and seeds/
    pins from the wrapped ``ic``.
    """
    nargs = (1, 2, 3, 4)


#: The time symbol usable in source terms: ``sin(2*pi*f*TIME)`` etc.
TIME = sympy.Symbol('_hdl_time', real=True)

#: Circuit temperature in kelvin -- Verilog-A's ``$temperature``; bound from
#: ``epar.T`` at evaluation time.
TEMP = sympy.Symbol('_hdl_temp', positive=True)

#: Noise frequency in Hz inside ``flicker_noise`` PSDs; bound from the noise
#: analysis's ``w`` at evaluation time.
FREQ = sympy.Symbol('_hdl_freq', positive=True)

#: Boltzmann constant and elementary charge, for ``vt()`` and noise PSDs.
#: Taken from the numeric toolkit so a generated noise stamp equals the
#: hand-written one exactly (``elements.R.CY``) rather than nearly.
from pycircuit.circuit.constants import kboltzmann as KBOLTZMANN
from pycircuit.circuit.constants import qelectron as QELECTRON


def vt(T=TEMP):
    """Thermal voltage ``k*T/q`` -- Verilog-A's ``$vt``."""
    return KBOLTZMANN * T / QELECTRON


class _wrapfloor(sympy.Function):
    """``floor`` whose symbolic derivative is 0 -- exact almost everywhere.

    The idtmod wrap ``s - m*floor((s-o)/m)`` has slope 1 a.e. and a jump
    of exactly ``m`` at the boundary that no Jacobian can state (idtmod.md
    sec. 5.1); an unevaluated ``Derivative(floor)`` would kill lambdify
    instead.  Printed as the toolkit's ``floor``.
    """
    nargs = (1,)

    def fdiff(self, argindex=1):
        return sympy.Integer(0)


def _conj(expr):
    """Complex conjugate of a model expression, by flipping ``I``.

    Every symbol a model can write -- a parameter, a node voltage, a
    branch current, the temperature, the frequency -- is REAL, so the
    only imaginary thing in an expression got there as an explicit
    `sympy.I`.  Substituting ``I -> -I`` is therefore exact, and unlike
    `sympy.conjugate` it leaves a real expression textually unchanged
    instead of burying every symbol in an unevaluated ``conjugate(...)``
    that the compiled matrix would then carry.
    """
    e = sympy.sympify(expr)
    return e.subs(sympy.I, -sympy.I) if e.has(sympy.I) else e


def _noise_name(args):
    """Sympify a trailing noise NAME as `Str`, not as a `Symbol`.

    Plain `sympify('igid')` yields ``Symbol('igid')``, which then shows up
    in the statement's ``free_symbols`` and is hunted for as a parameter --
    the name would have to be declared to be usable.  `Str` is atomic with
    an empty ``free_symbols``, so a name is inert everywhere except where
    this module deliberately reads it.
    """
    if args and isinstance(args[-1], str):
        return args[:-1] + (Str(args[-1]),)
    return args


class _NamedNoise(sympy.Function):
    """Shared plumbing for the noise functions: an optional trailing name.

    Two sources contributed with the SAME name are PERFECTLY CORRELATED
    (LRM 2.4 sec. 4.5.16) -- one physical fluctuation reaching the circuit
    through more than one branch.  Contributed without a name, a source is
    independent of every other.
    """

    def __new__(cls, *args, **kwargs):
        return super().__new__(cls, *_noise_name(args), **kwargs)

    @property
    def noise_name(self):
        """The correlation group, or None when this source stands alone."""
        last = self.args[-1]
        return last.name if isinstance(last, Str) else None


class white_noise(_NamedNoise):
    """``white_noise(pwr[, name])`` -- a flat current PSD on the branch.

    Contributes ONLY to the noise correlation matrix ``CY`` (LRM: noise
    functions return zero in DC/transient); ``pwr`` may use ``TEMP``.

    A scale factor multiplies the AMPLITUDE, so ``k * white_noise(p)``
    has PSD ``k**2 * p``.  That is what makes a named group composable:
    the factors carry the (possibly complex, possibly signed) transfer
    from the one fluctuation to each branch it reaches.
    """
    nargs = (1, 2)


class flicker_noise(_NamedNoise):
    """``flicker_noise(pwr, exp[, name])`` -- a ``pwr/f^exp`` current PSD."""
    nargs = (2, 3)


def limexp(x, x0=80.0):
    """Verilog-A ``limexp``: ``exp`` continued linearly above ``x0``.

    C1-continuous at the seam, so Newton keeps a bounded derivative
    instead of overflowing -- the standard convergence aid.

    Deliberately NOT both-arms-safe, unlike `expl`.  Both arms of a
    compiled conditional are evaluated, so the discarded ``exp(x)`` does
    overflow to ``inf`` above ``x0`` and raises a floating-point warning;
    the returned value is unaffected.  Clamping the argument would fix
    that, but it also hides the exponential: PCNR's shape detector reads
    ``exp(arg)`` straight out of the expression to recover the device's
    ``IS`` and ``VT`` (sec. "PCNR participation"), and an
    ``exp(min(arg, 80))`` is not a junction it can recognise -- the
    device silently loses its limiting.

    That trade is worth taking only where PCNR is what bounds the
    argument: a limited junction voltage does not reach 80 thermal
    voltages.  So the condition is **PCNR eligibility, not "is it a
    junction"** -- and the two are much further apart than they look.
    `explain()` prints which one a model is.

    A device qualifies only if it introduces no state and no branch
    unknown and EVERY current it contributes is an exponential function
    of its own branch voltage alone.  (Charge does not disqualify a
    device: `pcnr.augmented_system` used to refuse one and no longer
    does -- `test_pcnr_charge.py`.  Neither does `var()`, since
    roadmap 10.2: the detector walks the let-chain instead of reading
    the assembled expression, so a chained junction qualifies exactly
    as its flat twin does.)

    Measured on `elements_hdl.DiodeSpiceHdl`, a `var()` chain with an
    optional series node: `hasattr(el, 'pcnr_i')` is True for the
    default card and for `cjo = tt = rs = 0`, and False for
    `rs = 2`.  The refusal is worth stating exactly, because the
    obvious reason is the wrong one: with `rs` present, EACH of the two
    contributions is still a function of its own branch voltage alone
    -- what refuses the device is that the series resistor's current is
    not exponential, and the rule is every current, not some current.

    `compact.PspMosLongChannel` is refused for the same reason and at
    its very first contribution: `I(g,gi) <+ V(g,gi)/rg`, a plain
    resistor.  Its drain current would fail anyway -- it is
    irreducibly `f(vgs, vds, vbs)` and no scalar branch unknown
    expresses it (roadmap 10.3) -- but that is not the refusal that
    fires.

    A model that cannot qualify pays `limexp`'s deliberate unsafety for
    nothing: 120 `overflow encountered in exp` warnings on a single 5 V
    self-heating DC solve, none under `expl`, and the same solution to
    the last digit.  **Reach for `expl` unless the model is a bare
    exponential junction that PCNR actually takes** -- FET models
    included, and any junction that carries a series resistance as a
    second contribution.
    """
    x = sympy.sympify(x)
    return sympy.Piecewise((sympy.exp(x), x < x0),
                           (sympy.exp(x0) * (1 + x - x0), True))


## ----------------------------------------------------------------------
## The compact-model math kernel.
##
## These are the functions a surface-potential model is written in.  PSP103
## reaches for `expl` 23 times and `hypsmooth` throughout; every one of
## them exists to keep an expression finite where the physics would
## otherwise divide by zero or overflow.  They are provided here rather
## than left to each model because getting them both-arms-safe is fiddly
## and the failure is silent (a warning, or a NaN under stricter error
## settings), not a wrong answer.
## ----------------------------------------------------------------------

#: Argument magnitude beyond which `exp` is continued rather than
#: evaluated: PSP's `se05`, which is ``ln(1e100)``.  Both PSP103
#: (`Common103_macrodefs.include:69`) and the PSP-based `mosvar` use this
#: value, NOT the 80 that `limexp` conventionally uses.
EXPL_THRESHOLD = 2.3025850929940458e+02

#: ``exp(EXPL_THRESHOLD)`` and its reciprocal, as PSP spells them.
_KE05 = 1.0e-100
_KE05INV = 1.0e100


## ----------------------------------------------------------------------
## Selection, sign and magnitude.
##
## `sympy.Max`, `sympy.Min`, `sympy.sign` and `sympy.Abs` all exist and
## all four of them are traps in a compiled model:
##
##   * `Max`/`Min` differentiate to a `Heaviside`, which sympy's numpy
##     printer lowers to `numpy.select` -- 4.4 us per call on this
##     machine against 0.6 us for `numpy.maximum`, scalar, and a
##     compiled model calls it in the Newton inner loop.  An eagerly
##     compiled `expl` derivative emits FIVE of them; the chained one
##     emits 42 `numpy.where` calls, and 6 once the clamps are the
##     functions below;
##   * their SECOND derivative is a `DiracDelta`, which no numeric
##     backend prints -- so a distortion analysis or a sensitivity over
##     a model containing a `Max` fails to compile;
##   * applying `Max` to a COMPOUND expression is separately hazardous.
##     The rule of thumb was "`Max`/`Min` go on ATOMS", justified by a
##     Jacobian that divides by zero.  That symptom did NOT reproduce on
##     sympy 1.14 (`test_hdl_kernel.py` records the attempt); a worse one
##     did -- `sympy.Max` of a kernel compound does not finish `lambdify`
##     at all, where `maxc` of the same compound takes 6 ms;
##   * `sign` differentiates to `DiracDelta`, so a model using it fails
##     to COMPILE;
##   * `Abs` differentiates through `re`/`im` and is `0/0` at zero --
##     not because sympy cannot know a model symbol is real, but because
##     only SOME of them are.  `Quantity` declares `is_real`; a PARAMETER
##     is a bare `sympy.Symbol` and does not, so `Abs` of an expression
##     mixing an unknown with a parameter takes the complex path.
##
## The functions below are the same values with derivatives chosen to be
## the almost-everywhere-correct ones: defined at the kink,
## `DiracDelta`-free at every order, and lowered to a cheap primitive.
## ----------------------------------------------------------------------


class _step(sympy.Function):
    """``1.0`` where ``a >= b``, ``0.0`` below -- the a.e. slope selector.

    This is what `maxc`, `minc` and `sign` differentiate to, and it
    exists so that they do not differentiate to `Heaviside` instead.
    Two reasons, and the second is the one that matters:

    * `Heaviside`'s own derivative is a `DiracDelta`, so a model that
      merely took a SECOND derivative of a `Max` -- a distortion
      analysis, a sensitivity -- would stop compiling.  `_step` has
      ``fdiff -> 0``, which is exact everywhere except the measure-zero
      tie, and there no finite Jacobian exists to be right about;
    * sympy's numpy printer lowers `Heaviside` to `numpy.select`.

    The implementation is one comparison and one multiplication, written
    in plain Python operators so that the SAME function serves numpy
    arrays, numpy scalars and jax tracers.  Comparisons are total on the
    extended reals, so it is exact at ``+-inf`` as well, and on two
    sympy expressions the multiplication raises `TypeError` promptly --
    which is what `evalf` wants to hear, and is why this one needs no
    explicit symbolic guard where `_maxc_numpy` does.
    """
    nargs = (2,)

    def fdiff(self, argindex=1):
        ## Zero on both sides.  The true derivative is zero everywhere
        ## except the tie, where it is not a function at all.
        return sympy.Integer(0)

    @classmethod
    def eval(cls, a, b):
        if a.is_number and b.is_number and a.is_real and b.is_real:
            return sympy.Float(1.0) if bool(a >= b) else sympy.Float(0.0)

    ## `_imp_` is sympy's hook for "this has a numeric implementation",
    ## so a PLAIN `lambdify` works with no modules map.  Safe to give
    ## here -- and only here in this block -- because the implementation
    ## is operator-only and therefore backend-neutral; `_imp_` takes
    ## PRECEDENCE over a modules map, so a numpy-flavoured one would
    ## quietly break the jax path.
    @staticmethod
    def _imp_(a, b):
        return 1.0 * (a >= b)


def _step_numpy(a, b):
    """Runtime form of `_step`.  Kept beside it so the two cannot drift."""
    return 1.0 * (a >= b)


def _maxc_numpy(a, b):
    """Runtime form of `maxc` -- ``numpy.maximum``, or jax's on a tracer.

    Three requirements pull in different directions and this is the only
    shape that meets all three.

    * it must be the backend's own ``maximum``, not arithmetic.  The
      operator form ``a*s + b*(1 - s)`` would be backend-neutral for
      free, and it returns ``nan`` whenever the LOSING argument is
      infinite (``0 * inf``) -- exactly the case a clamp is written to
      survive;
    * it must serve numpy AND jax from ONE function, because sympy's
      ``_imp_`` hook takes precedence over any modules map, so a
      numpy-only implementation would not merely be slow inside a traced
      model, it would stop it compiling
      (``TracerArrayConversionError``);
    * it must cost nothing on the numpy path, which is the scalar Newton
      inner loop.

    Hence the try/except rather than a type test: an untaken ``except``
    is free against ``numpy.maximum``'s own ~0.7 us, where inspecting
    both arguments' types on every call is not.  A jax TRACER raises
    (``TracerArrayConversionError`` is a ``TypeError``) and is caught
    once per trace, not once per evaluation -- after that the
    ``jnp.maximum`` is baked into the jaxpr.

    The SYMBOLIC guard is not optional and is not cheap to omit.  Sympy's
    ``evalf`` calls ``_imp_`` speculatively with the function's own
    sympy arguments; ``numpy.maximum`` on two expressions then evaluates
    ``a <= b``, which builds a `Relational`, which stringifies the whole
    expression to compose its error message.  On a `safe_sqrt` derivative
    that does not finish.  Measured: an unguarded version turned one
    `lambdify` into a hang.  Raising `TypeError` is exactly what evalf
    expects for "no numeric value here".
    """
    if isinstance(a, sympy.Basic) or isinstance(b, sympy.Basic):
        raise TypeError('maxc has no symbolic evaluation')
    try:
        return np.maximum(a, b)
    except TypeError:
        import jax.numpy as jnp
        return jnp.maximum(a, b)


def _minc_numpy(a, b):
    """Runtime form of `minc`; see `_maxc_numpy` for the shape."""
    if isinstance(a, sympy.Basic) or isinstance(b, sympy.Basic):
        raise TypeError('minc has no symbolic evaluation')
    try:
        return np.minimum(a, b)
    except TypeError:
        import jax.numpy as jnp
        return jnp.minimum(a, b)


class maxc(sympy.Function):
    """``max(a, b)`` that may be applied to a COMPOUND expression.

    The two-argument `sympy.Max` with a derivative that does not expand
    its arguments: ``d/da = _step(a, b)``, ``d/db = 1 - _step(a, b)``.
    That is the whole point.  The rule it lifts -- "`Max`/`Min` go on
    ATOMS, bind the argument with `var()` first", undocumented and
    obeyed by construction at about fifty sites -- is real, though not
    for the reason usually given.  The stated symptom, a Jacobian that
    divides by zero, did not reproduce on sympy 1.14.  What does happen,
    measured, is that ``sympy.Max(<a kernel expression>, c)`` no longer
    finishes `lambdify`: the `Heaviside` printer is left comparing
    `Float`s inside the compound and does not come back within 200
    seconds, where `maxc` of the same compound compiles in 6 ms.

    The two partials sum to 1 at every point INCLUDING the tie, where
    `Heaviside` would give 1/2 to each and `maxc(x, x)` would come out
    with slope 2.

    RANGE: total.  Exact for every pair of doubles, ``+-inf`` included,
    because it lowers to the backend's own ``maximum`` rather than to
    arithmetic.  That distinction is not cosmetic: the obvious
    operator-only form ``a*s + b*(1 - s)`` returns ``nan`` whenever the
    LOSING argument is infinite (``0 * inf``), which is exactly the case
    a clamp is most often written to survive.
    """
    nargs = (2,)

    def fdiff(self, argindex=1):
        a, b = self.args
        s = _step(a, b)
        return s if argindex == 1 else 1 - s

    ## `_imp_` so that a PLAIN `sympy.lambdify` works with no modules
    ## map.  That is not a nicety: `hypsmooth` and the `expl` family are
    ## PUBLIC and are clamped with these, so without it every outside
    ## caller who lambdified `hdl.safe_sqrt(x)` with `modules='numpy'`
    ## -- which worked while the clamp was `sympy.Max` -- would get a
    ## `NameError` from inside `<lambdifygenerated-N>`.
    ##
    ## `_imp_` takes PRECEDENCE over a modules map, so it has to serve
    ## every backend by itself; see `_maxc_numpy` for how.
    @staticmethod
    def _imp_(a, b):
        return _maxc_numpy(a, b)

    @classmethod
    def eval(cls, a, b):
        ## Fold two literals, so `maxc(1e-15, 0)` does not survive into
        ## the generated source as a call.
        if a.is_number and b.is_number and a.is_real and b.is_real:
            return sympy.Max(a, b)


class minc(sympy.Function):
    """``min(a, b)`` on a compound expression -- the twin of `maxc`.

    ``d/da = _step(b, a)`` (one where ``a <= b``), ``d/db`` its
    complement, so the partials again sum to exactly 1 at the tie.

    RANGE: total, for the same reason as `maxc`.
    """
    nargs = (2,)

    def fdiff(self, argindex=1):
        a, b = self.args
        s = _step(b, a)
        return s if argindex == 1 else 1 - s

    @staticmethod
    def _imp_(a, b):
        return _minc_numpy(a, b)

    @classmethod
    def eval(cls, a, b):
        if a.is_number and b.is_number and a.is_real and b.is_real:
            return sympy.Min(a, b)


class sign(sympy.Function):
    """``sign(x)`` that COMPILES: ``-1``, ``0``, ``+1``, with ``fdiff -> 0``.

    `sympy.sign` differentiates to a `DiracDelta`, which no numeric
    backend prints -- so a model that used it failed at the lambdify
    step, not at the Newton step, and the failure mode was a `NameError`
    out of ``<lambdifygenerated-N>``.  The measured consequence is that
    conditions get rewritten algebraically to dodge it (``x*x < 1e-6``
    standing in for ``Abs(x) < 1e-3``), which is correct and unreadable.

    The derivative here is 0, which is exact everywhere except the origin
    -- and at the origin no finite derivative exists, so 0 is the only
    answer a Jacobian can carry.  Being `DiracDelta`-free at EVERY order
    also means a second derivative still compiles.

    The class is deliberately named ``sign``: sympy's numpy printer maps
    that name to the backend's own ``sign``, so it lowers to
    ``numpy.sign`` on numpy and ``jnp.sign`` on jax with no registration
    at all.

    RANGE: total.  Exact for every double including ``+-inf``; ``nan``
    in gives ``nan`` out, as it should.
    """
    nargs = (1,)

    def fdiff(self, argindex=1):
        return sympy.Integer(0)

    @classmethod
    def eval(cls, x):
        if x.is_number and x.is_real:
            return sympy.sign(x)


class Abs(sympy.Function):
    """``|x|`` on the REALS, with the real derivative ``sign(x)``.

    `sympy.Abs` is a quiet trap: the value is always right, and it
    differentiates to ``(re(u)*re'(u) + im(u)*im'(u))/Abs(u)``, which is
    ``0/0`` at the origin.  The reason is narrower than "sympy does not
    know a model symbol is real", and the narrowing matters: `Quantity`
    DOES declare ``is_real = True``, so ``Abs`` of a bare branch voltage
    is fine.  A PARAMETER is a plain ``sympy.Symbol`` with no
    assumptions (see `generate_code`), so ``Abs`` of anything MIXING an
    unknown with a parameter -- which is most model expressions -- takes
    the complex path.  A model is then finite everywhere and has a NaN
    Jacobian at the one bias where its argument vanishes -- the worse of
    the two failures, because nothing looks wrong.  Raising an `Abs` to a power
    compounds it, the printed derivative carrying ``sign(u)/(u*Abs(u))``.

    This one differentiates to `sign`, above: defined at the origin,
    where it takes the value 0, and `DiracDelta`-free at every order.

    Choose between this and `safe_abs` by what the consumer does with
    the result.  `safe_abs` is ``sqrt(x*x + eps*eps)``, SMOOTH at the
    origin and therefore the right thing when the result is divided by,
    logged, or raised to a power; it costs half the exponent range (it
    squares its argument, so ``|x| < ~1e150``) and it is not exactly
    ``|x|``.  `Abs` is exact, total over the whole double range, costs
    nothing, and has a kink -- the right thing when the result is added,
    compared, or multiplied.

    RANGE: total.  Exact for every double including ``+-inf``.

    The class is named ``Abs`` so that sympy's numpy printer lowers it to
    the backend's own ``abs``, on numpy and on jax alike, with no
    registration.  It is a drop-in replacement for ``sympy.Abs`` in a
    model.
    """
    nargs = (1,)

    def fdiff(self, argindex=1):
        if argindex != 1:
            return sympy.Integer(0)
        return sign(self.args[0])

    @classmethod
    def eval(cls, x):
        if x.is_number and x.is_real:
            return sympy.Abs(x)



def p3(u):
    """PSP's ``P3`` -- the third-order Taylor polynomial of ``exp``.

    ``1 + u + u**2/2 + u**3/6``.  Continuing ``exp`` with its own Taylor
    series about the seam is what makes the `expl` family **C-3**
    continuous rather than merely C1: three derivatives match, not one.

    Public because it is the generic way to continue ANY bounded
    exponential: a model that needs its own seam gets the same C-3 join
    without re-deriving the coefficients, and `psp_kernel` had a second,
    identical copy of it for exactly that reason.

    RANGE: a plain cubic, so it is finite for ``|u| < ~1e102`` and
    overflows above that.  Its derivative ``1 + u + u**2/2`` overflows a
    little later, at ``~1e154`` -- the one member of the kernel whose
    derivative outranges its value.  It is meant to be evaluated on the
    OVERSHOOT past a seam, not on the raw argument.
    """
    u = sympy.sympify(u)
    return 1.0 + u * (1.0 + 0.5 * (u * (1.0 + u / 3.0)))


#: Historic private spelling of `p3`, kept because the `expl` family
#: below is written in terms of it.
_p3 = p3


def expl_high(x, x0=EXPL_THRESHOLD):
    """``exp`` continued above ``+x0`` -- PSP's ``expl_high``.

    Above the seam the value is ``1e100 * P3(x - x0)``.  Note this is a
    genuinely different function from `limexp`, which continues linearly
    from 80: at ``x = 100`` this returns ``exp(100)`` exactly, where
    `limexp` returns ``exp(80)*21``, seven orders of magnitude away.  Use
    this one when translating a PSP-family model, and `limexp` when you
    want the LRM's function.
    """
    x = sympy.sympify(x)
    ## Each arm sees an argument its own branch makes valid, so the
    ## discarded one cannot overflow -- see the module's kernel notes.
    ## The clamps are `maxc`/`minc` rather than `sympy.Max`/`Min`: the
    ## VALUE is identical and the derivative is the same function, but
    ## `Max`'s derivative is a `Heaviside` per occurrence, and this
    ## family carries four of them.  Measured on a chained element's `G`,
    ## `expl` went from 42 dispatched `numpy.where` calls at 140 us to 6
    ## at 58 us, and `hypsmooth` from 19 at 56 us to 3 at 32 us.
    return sympy.Piecewise(
        (sympy.exp(minc(x, x0)), x < x0),
        (_KE05INV * _p3(maxc(x, x0) - x0), True))


def expl_low(x, x0=EXPL_THRESHOLD):
    """``exp`` continued below ``-x0`` -- PSP's ``expl_low``.

    Below the seam the value is ``1e-100 / P3(-x0 - x)``, which stays
    strictly POSITIVE as ``x -> -inf`` where plain ``exp`` underflows to
    exactly zero around -745.  That is the property the surface-potential
    solver relies on: it divides by the result.
    """
    x = sympy.sympify(x)
    return sympy.Piecewise(
        (_KE05 / _p3(-x0 - minc(x, -x0)), x < -x0),
        (sympy.exp(maxc(x, -x0)), True))


def expl(x, x0=EXPL_THRESHOLD):
    """Two-sided clamped exponential -- PSP's ``expl``.

    ``exp`` between the seams, its own third-order Taylor continuation
    outside them, C-3 at both.  Strictly positive everywhere, and finite
    over the whole range a circuit is ever solved in, which is what
    makes a surface-potential model evaluable at the absurd biases
    Newton visits on its way to the solution.

    This docstring used to say "finite for every double", and that is
    FALSE: the cubic continuation overflows at an argument of 4.76e69,
    where it returns 1.798e308.  Measured 2026-08-25 while writing a
    level-1 MOSFET, which stays finite through 1e68 V on its bulk node,
    so the gap is academic -- but a guarantee that is nearly true is
    the kind that gets relied on at the one bias where it is not.

    PSP103 and JUNCAP200 call the family 31 times between them.
    """
    x = sympy.sympify(x)
    return sympy.Piecewise(
        (_KE05 / _p3(-x0 - minc(x, -x0)), x < -x0),
        (_KE05INV * _p3(maxc(x, x0) - x0), x > x0),
        (sympy.exp(maxc(minc(x, x0), -x0)), True))


def hypsmooth(x, eps):
    """Smooth ``max(x, 0)`` -- PSP's ``hypsmooth``.

    Mathematically ``0.5*(x + sqrt(x*x + 4*eps^2))``, which is strictly
    positive for every real ``x`` and tends to ``eps^2/|x|`` as
    ``x -> -inf``.  That positivity is the whole value of the function:
    downstream code divides by it and takes its logarithm.

    Written literally, floating point destroys it.  At ``x = -100`` with
    ``eps = 1e-12``, ``sqrt(x*x + 4e-24)`` rounds to exactly ``100.0``,
    the sum cancels to exactly ``0.0``, and every "guaranteed positive"
    consumer downstream gets a zero.  (Measured: this is what turned
    `safe_ln` into 2012 infinities over a 4001-point sweep.)

    So the negative side is evaluated through the conjugate form
    ``2*eps^2 / (sqrt(x*x + 4*eps^2) - x)``, algebraically identical but
    a sum of positives in the denominator, with nothing to cancel.  Both
    arms are valid for any input, so it does not matter that both are
    evaluated.

    RANGE: the radicand squares ``x``, so the value is finite for
    ``|x| < ~1e154`` and overflows above it.  Circuit quantities DO
    reach there -- see `safe_ln`, where an intermediate hit 4.4e222 --
    so the bound is a real one and not a theoretical one; clamping the
    argument to buy it back would cost a comparison on every evaluation
    of the most frequently called function in the kernel.  The
    derivative is bounded by 1 and costs no further range, which makes
    `hypsmooth` the one member of the family whose Jacobian does NOT
    die before its value.

    ``eps`` is squared, so it cannot go below about 1.5e-154 without
    underflowing to exactly zero and removing the regularisation
    entirely.
    """
    x, eps = sympy.sympify(x), sympy.sympify(eps)
    ## Each arm sees an argument clamped to ITS OWN side.  Without that
    ## the conjugate arm, evaluated at large POSITIVE x, has `r - x`
    ## cancel to exactly zero and divides by it -- the mirror image of
    ## the cancellation it was written to avoid.
    xp, xn = maxc(x, 0.0), minc(x, 0.0)
    return sympy.Piecewise(
        (0.5 * (xp + sympy.sqrt(xp * xp + 4 * eps * eps)), x >= 0),
        (2 * eps * eps / (sympy.sqrt(xn * xn + 4 * eps * eps) - xn), True))


def safe_sqrt(x, eps=1e-12):
    """``sqrt(x)`` made finite and differentiable for every real ``x``.

    Clamping the argument -- ``sqrt(max(x, 0))`` -- fixes the VALUE but
    not the derivative: sympy differentiates it to
    ``d(max)/dx / (2*sqrt(max(x, 0)))``, which at ``x <= 0`` is ``0/0``,
    and a NaN in the Jacobian kills the whole Newton step rather than one
    entry.  (Measured: 2013 non-finite derivatives over a 4001-point
    sweep before this.)

    So the argument is smoothed rather than clipped:
    ``sqrt(hypsmooth(x, eps))``.  ``hypsmooth`` is strictly positive
    everywhere -- it decays like ``eps^2/|x|`` as ``x -> -inf`` instead of
    reaching zero -- so the square root and its derivative are finite for
    every input, with no branch at all.  For ``x >> eps`` it agrees with
    ``sqrt(x)`` to relative order ``(eps/x)^2``.

    RANGE: inherits `hypsmooth`'s, which SQUARES its argument -- the
    value is finite for ``|x| < ~1e154`` and overflows above that, where
    plain ``sqrt`` would have gone to 1e77 without blinking.  The
    derivative is ``1/(2*sqrt(hypsmooth))``, bounded by ``1/(2*eps)`` --
    with the default ``eps = 1e-12``, 5e11 -- and it reaches that bound
    at ``x = 0``, not at the extremes.  Below zero the value decays like
    ``eps/sqrt(|x|)`` and the derivative with it, so the negative tail
    costs nothing.

    ``eps`` is squared inside `hypsmooth`, so it cannot go below about
    1.5e-154 without underflowing to zero and taking the regularisation
    with it.
    """
    return sympy.sqrt(hypsmooth(x, eps))


def safe_abs(x, eps=1e-30):
    """``abs(x)`` made DIFFERENTIABLE at zero, as ``sqrt(x*x + eps^2)``.

    `sympy.Abs` is the trap here, and it is a quiet one because the
    VALUE is always right: it differentiates to
    ``(re(u)*re'(u) + im(u)*im'(u)) / Abs(u)``, which at ``u = 0`` is
    ``0/0``, whenever ``u`` mixes an unknown with a PARAMETER -- see
    `Abs`, above, for exactly when that is and why.  A model can then be
    finite everywhere and have a NaN Jacobian at the one bias where its
    argument vanishes, which is the worse of the two failures: the value
    looks right and Newton is poisoned.  Raising an ``Abs`` to a power
    compounds it, the printed derivative carrying ``sign(u)/(u*Abs(u))``.

    Choose between this and `Abs` by what the consumer does with the
    result: this one is SMOOTH and costs half the exponent range, `Abs`
    is exact and total and has a kink.

    ``sqrt(x*x + eps^2)`` has derivative ``x/sqrt(x*x + eps^2)``, which
    is bounded by 1 and equals 0 at the origin -- finite for every real
    input, with no branch.  For ``|x| >> eps`` it agrees with ``abs(x)``
    to relative order ``(eps/x)^2 / 2``.

    Squares its argument, so it halves the exponent range: keep ``|x|``
    under about 1e150, which for a voltage-derived quantity means
    clamping the voltage first (`juncap.VJUN_CLAMP` is what does that
    for the junction).
    """
    return sympy.sqrt(x * x + eps * eps)


def safe_ln(x, eps=1e-30):
    """``ln(x)`` made finite for every real ``x``, by the same smoothing.

    ``ln(hypsmooth(x, eps))``: strictly positive argument, so no
    ``-inf`` at zero and no NaN below it, and no branch.  Note the result
    is genuinely large and negative for small ``x`` -- that is ``ln``
    doing its job, not an overflow.

    RANGE: this is the helper whose cost was MEASURED and is worth
    stating twice.  It goes through `hypsmooth`, so it SQUARES its
    argument and is finite only for ``|x| < ~1e154`` -- where plain
    ``log`` is finite for every positive double up to 1.8e308.  PSP's
    gate tunnelling reaches an argument of 4.4e222 at ``Vds = 1e12`` and
    the drain current came back ``-inf`` from an expression every
    ingredient of which was finite.  The derivative is
    ``1/hypsmooth(x, eps)``, bounded by ``1/eps`` -- 1e30 with the
    default -- and reached at ``x = 0``.

    **Ask what the guard is for before buying it.**  The arguments in
    that measurement were ``1 + expl(...)``, at least 1 by construction,
    so the guard could never fire and the 154 decades it cost were a
    pure loss.  If the argument is provably positive, plain
    ``sympy.log`` is both correct and 154 decades wider.
    """
    return sympy.log(hypsmooth(x, eps))


class _recip2(sympy.Function):
    """``1/(b^2 + eps^2)``, with a derivative written in terms of ITSELF.

    This exists for one reason: exponent range.  Any rational
    regularisation has its denominator SQUARED in its own derivative, so
    the naive ``a*b/(b^2 + eps^2)`` is differentiated to something
    carrying ``(b^2 + eps^2)^2``.  That overflows at ``|b| ~ 1e77``,
    while the VALUE survives to ``1e154`` -- so the Jacobian dies two
    decades of exponent before the residual does, which is the worse way
    round: the value looks right and Newton is poisoned.

    Writing the derivative as ``-2*b*self**2`` instead forms the square
    of the RECIPROCAL, which for a large denominator underflows to zero
    gracefully rather than overflowing to ``inf`` and then to ``NaN``.
    The two are identical mathematically; only the intermediate differs.
    Measured: the Jacobian's range goes from ``1e77`` to ``1e154``, which
    is the value's own limit and cannot be bettered without changing the
    value too.

    Printed as ``1/(b*b + e*e)`` through the modules map, the same way
    `_wrapfloor` is printed as the toolkit's ``floor``.
    """
    nargs = (2,)

    def fdiff(self, argindex=1):
        if argindex != 1:
            ## `eps` is a constant of the model, never differentiated.
            return sympy.Integer(0)
        b = self.args[0]
        return -2 * b * self ** 2

    ## `_imp_` is sympy's own hook for "this function has a numeric
    ## implementation".  With it, a PLAIN `sympy.lambdify` works with no
    ## modules map at all -- which matters because `safe_div` is public
    ## and callers outside this module lambdify its result directly.
    ## Requiring a namespace entry would have made the primitive a
    ## breaking change to a published helper rather than an internal
    ## improvement to it.
    @staticmethod
    def _imp_(b, e):
        return 1.0 / (b * b + e * e)


def _recip2_numpy(b, e):
    """Runtime form of `_recip2`.  Kept beside it so the two cannot drift."""
    return 1.0 / (b * b + e * e)


class _rdiv(sympy.Function):
    """``b/(b^2 + eps^2)`` -- the whole of `safe_div` bar its numerator.

    A primitive rather than an expression so its DERIVATIVE can be
    GROUPED, which is the entire point.  Written out, the derivative is
    ``(eps^2 - b^2)/(b^2 + eps^2)^2``; whichever way that is arranged as
    a sum it contains ``inv^2``, and for a large denominator ``inv^2``
    underflows to zero while a surviving ``+inv`` term does not -- so the
    derivative comes out ``+1/b^2`` where it should be ``-1/b^2``.
    **Finite, plausible, and the wrong sign**, which is worse than the
    overflow it replaced, because nothing looks wrong.

    Grouped as ``inv * (1 - 2*b*b*inv)`` there is no ``inv^2`` at all:
    ``b*b*inv`` is ``b^2/(b^2 + eps^2)``, which is bounded in ``[0, 1]``
    for every real ``b``, so the product is bounded by ``inv`` itself.
    Correct sign and correct magnitude for every ``b`` the VALUE can be
    computed at, which is the most that can be asked.
    """
    nargs = (2,)

    def fdiff(self, argindex=1):
        if argindex != 1:
            return sympy.Integer(0)
        b, e = self.args
        inv = _recip2(b, e)
        return inv * (1 - 2 * b * b * inv)

    @staticmethod
    def _imp_(b, e):
        return b / (b * b + e * e)


def _rdiv_numpy(b, e):
    """Runtime form of `_rdiv`."""
    return b / (b * b + e * e)


def safe_div(a, b, eps=1e-30):
    """``a/b`` regularised to ``a*b/(b^2 + eps^2)``.

    The obvious guard -- shift the denominator by ``eps*sign(b)`` -- is
    unusable: sympy differentiates ``sign`` to a ``DiracDelta`` that no
    numeric backend can print, so the model fails to compile rather than
    to converge.

    This form has no ``sign`` and no branch.  It equals ``a/b`` to
    relative order ``(eps/b)^2`` wherever ``|b| >> eps``, tends smoothly
    to zero as ``b -> 0`` instead of blowing up, and its derivative is
    finite everywhere.  The price is that it is NOT ``a/b`` when ``b`` is
    genuinely of order ``eps`` -- pick ``eps`` below any denominator the
    model can legitimately produce.

    The VALUE needs ``|b| < ~1e154``, where ``b*b`` leaves double range,
    and the DERIVATIVE now needs the same rather than ``|b| < ~1e77``:
    see `_recip2`, which exists exactly to buy back those two decades.

    **``eps`` is SQUARED**, so it cannot be made arbitrarily small: below
    about 1.5e-154 the square underflows to exactly zero and the
    regularisation silently disappears, leaving the ``0/0`` it exists to
    prevent.  Measured -- `psp_kernel` first used ``eps = 1e-300`` and got
    a NaN Jacobian at flat band, where the denominator is exactly zero.
    Refused here rather than left to surface as a NaN.
    """
    a, b = sympy.sympify(a), sympy.sympify(b)
    if float(eps) * float(eps) == 0.0:
        raise ValueError(
            'safe_div eps=%r is too small: it is squared, and %r**2 '
            'underflows to exactly 0.0, so the denominator is left '
            'unregularised and b == 0 gives 0/0. Use eps >= 1e-150.'
            % (eps, eps))
    ## Built on `_rdiv` rather than written out, so the DERIVATIVE stays
    ## in double range as far as the value does -- see that class.  The
    ## value printed is identical.
    return a * _rdiv(b, sympy.Float(eps))


def safe_pow(b, e, lo=1e-30, hi=None):
    """``b**e`` with the base CONFINED to ``[lo, hi]`` -- a documented floor.

    ``b**e`` is not a total function.  For ``b < 0`` and a non-integer
    ``e`` it is undefined and the backend returns ``nan``; at ``b = 0``
    with ``e < 1`` the VALUE may be fine while the derivative
    ``e*b**e/b`` diverges.  Both cases are reachable at ordinary bias --
    PSP's Coulomb-scattering ratio is exactly zero at flat band and its
    exponent is 0.59 -- so every power in a real model is hand-floored,
    and this is that floor with its reasoning attached.

    The clamp is a HARD one (`maxc`, and `minc` when ``hi`` is given),
    not a smoothing, and that is deliberate.  ``sqrt(max(x, 0))`` is the
    classic broken clamp because the floor is at zero and the outer
    function is singular there; with a floor strictly INSIDE the domain
    the derivative below it is ``e*lo**(e-1) * 0``, an ordinary finite
    number times zero.  Below the floor the term is constant, which is
    usually the physics as well: no charge, no scattering off it.

    ``lo`` MUST be given a value relative to the quantity's own scale.
    An absolute constant encodes an assumption about a parameter you do
    not control; the default 1e-30 is a backstop for a normalised ratio,
    not a recommendation.

    RANGE:

    * the value is bounded by ``max(lo**e, hi**e)`` and by ``b**e`` for
      an in-range base, so it is finite exactly while those are.  With
      the default floor ``lo = 1e-30`` that means ``e > -10.2``;
    * **the derivative dies first**, as it does for every regulariser
      here.  It carries one further factor ``1/base``, so it is bounded
      by ``|e| * lo**(e - 1)`` -- thirty decades above the value with the
      default floor, and finite only for ``e > -9.2``.  That is the
      number to budget against.  Both bounds are CHECKED at build time
      below for a numeric exponent, and an exponent the floor cannot
      carry is refused rather than left to surface as an ``inf``.  A
      SYMBOLIC bound -- ``lo = 1e-3*vp``, which is the recommendation --
      cannot be checked at compile time and is not;
    * outside ``[lo, hi]`` the result is CONSTANT and its derivative is
      exactly zero.  This is an approximation, not an identity: check
      that the model does not need the tail it truncates.
    """
    b, e = sympy.sympify(b), sympy.sympify(e)
    lo = sympy.sympify(lo)
    hi = None if hi is None else sympy.sympify(hi)
    ## `lo` MAY be an expression -- `1e-3*vp` rather than `1e-3` is the
    ## recommendation above, and a smoothing constant written relative to
    ## a parameter is a symbol at compile time.  Everything below that
    ## needs a number therefore checks first, and a symbolic bound simply
    ## goes unchecked rather than raising.
    if lo.is_number and float(lo) <= 0.0:
        raise ValueError(
            'safe_pow lo=%r must be strictly positive: the floor exists '
            'to keep the base off zero, where b**e has a divergent '
            'derivative for e < 1 and is undefined for b < 0.' % (lo,))
    if (hi is not None and hi.is_number and lo.is_number
            and float(hi) <= float(lo)):
        raise ValueError('safe_pow hi=%r is not above lo=%r' % (hi, lo))
    ## A build-time check that the guard's OWN value is representable.
    ## `lo` is chosen for the base and the exponent is chosen for the
    ## physics, and it is easy to pick a pair whose combination leaves
    ## double range -- at which point the floor that was supposed to keep
    ## the term finite is what makes it infinite.
    if e.is_number and e.is_real:
        ev = float(e)
        bounds = (lo,) if hi is None else (lo, hi)
        for lim in (q for q in bounds if q.is_number and q.is_real):
            ## Budgeted TWICE, value and derivative: the derivative
            ## carries one more `1/base` and so leaves double range
            ## thirty decades of exponent before the value does.
            for what, power, scale in (('value', ev, 1.0),
                                       ('derivative', ev - 1.0, abs(ev))):
                try:
                    val = scale * float(lim) ** power
                except OverflowError:
                    val = float('inf')
                if not np.isfinite(val):
                    raise ValueError(
                        'safe_pow: the clamp itself overflows -- the %s at '
                        'base %r with exponent %r is not representable, so '
                        'the floor that is supposed to keep this term '
                        'finite is what makes it infinite. Raise lo (or '
                        'lower hi) until it is in range.'
                        % (what, float(lim), ev))
    base = maxc(b, lo)
    if hi is not None:
        base = minc(base, hi)
    return base ** e


def softplus(z):
    """``log(1 + exp(z))``, branch-free and finite for every double.

    Written ``max(z, 0) + log1p(exp(min(z, 0) - max(z, 0)))``, so the
    exponential only ever sees a NON-POSITIVE argument and is therefore
    in ``(0, 1]``.  The literal ``log(1 + exp(z))`` overflows at
    ``z = 710`` and the value it was going to return was ``710``.

    ``log1p`` rather than ``log(1 + ...)`` for the other end of the same
    range.  At ``z = -50`` the exponential is 1.9e-22, ``1 + 1.9e-22``
    rounds to exactly 1.0, and ``log`` of that is exactly 0 -- so the
    branch-free form, having been built to keep the large end finite,
    silently loses the small end entirely.  ``log1p`` carries it to
    ``z = -700`` (measured: 9.86e-305, correct to every digit).  This is
    the one respect in which the kernel version differs from
    `psp_kernel._softplus`, which it otherwise reproduces.

    The reason this belongs in the kernel rather than in one model is
    the rule it embodies: **a PRODUCT of two bounded exponentials is not
    bounded.**  `expl` is finite for every double, but past its seam it
    continues polynomially rather than saturating, so ``expl(a)*expl(b)``
    reaches 1e186 apiece and overflows.  Where only the LOGARITHM of
    such a product is needed, carry the exponents and never form the
    product: ``log(1 + e^a * e^b)`` is ``softplus(a + b)``, bounded by
    ``|a + b| + 0.7``.

    ``z`` appears three times in the result, so bind it with `var()`
    first -- see that function on why a compact model needs a let-chain.

    RANGE: total.  Finite and non-negative for every double, ``+-inf``
    included; the derivative is the logistic function and is bounded by
    1 everywhere.  It costs no exponent range at all, which makes it the
    cheapest member of the kernel to reach for.
    """
    z = sympy.sympify(z)
    zp, zm = maxc(z, 0.0), minc(z, 0.0)
    return zp + _log1p(sympy.exp(zm - zp))


def mne(x, y, a=0.0):
    """PSP's smooth MINIMUM of ``x`` and ``y`` (``macrodefs:40-43``).

    ``2/(4-a) * (s - sqrt(s^2 - (4-a)*x*y))`` with ``s = x + y``: exact
    at ``a = 0``, and with smoothing RELATIVE to the product ``x*y``
    rather than absolute, which is what suits a pair whose scale varies
    over decades.  Not the same function as `minc`, which is the exact
    minimum with a kink, nor as PSP's ``MINA``, whose smoothing is
    absolute.

    Evaluated in CONJUGATE form, ``2*x*y / (s + sqrt(...))``.  Written
    literally it subtracts two nearly equal numbers: at ``s = 1e9`` the
    difference is 1e-9 of the terms and rounds to zero, so the smooth
    minimum of two large positives comes out exactly 0.  The conjugate
    is a quotient of positives with nothing to cancel -- the same fix,
    for the same reason, as `hypsmooth`'s negative arm.

    RANGE: the radicand squares ``s``, so ``|x| + |y| < ~1e154``; beyond
    that it overflows.  For ``x, y > 0`` the result is in
    ``[0, min(x, y)]`` and the ``safe_div`` guard (``eps = 1e-30``) never
    fires, since the denominator is a sum of positives.  For a MIXED
    sign pair the radicand can exceed ``s^2`` and `safe_sqrt` is what
    keeps it real.
    """
    x, y, a = sympy.sympify(x), sympy.sympify(y), sympy.sympify(a)
    t1 = 4.0 - a
    s = x + y
    rad = safe_sqrt(s * s - t1 * (x * y))
    return 2.0 * (x * y) * safe_div(1.0, s + rad, eps=1e-30)


def mxe(x, y, a=0.0):
    """PSP's smooth MAXIMUM (``macrodefs:45-48``), the ``+`` twin of `mne`.

    ``2/(4-a) * (s + sqrt(s^2 - (4-a)*x*y))``.  No conjugate form is
    needed or possible: the two terms are ADDED, so there is nothing to
    cancel.

    RANGE: as `mne` -- the radicand squares ``s``, so ``|x| + |y| <
    ~1e154``.  For ``x, y > 0`` the result is in ``[max(x, y), 4/(4-a) *
    max(x, y)]``.
    """
    x, y, a = sympy.sympify(x), sympy.sympify(y), sympy.sympify(a)
    t1 = 4.0 - a
    s = x + y
    return (2.0 / t1) * (s + safe_sqrt(s * s - t1 * (x * y)))


## ----------------------------------------------------------------------
## Runtime registration of the kernel's primitives.
##
## A primitive is only half-built when the sympy class exists: it prints
## as a CALL, so every namespace a compiled model is executed in has to
## bind the name.  There are SIX such sites: the eager `NUMPY_MODULES`,
## the let-chain's `_mods`, the crossing-event map, the chain's own
## hand-built namespace in `_chain_compile`, and the two jax `lambdify`
## maps (PCNR's junction and the batched pure forms).  A name missing
## from one of them fails only for the models that take that path -- the
## chain namespace is reached only by a model that uses `var()`, which is
## every production model and no toy one, and the crossing map only by a
## model that asks the transient solver for a timepoint.
##
## So the map is written ONCE, here, and every registration site is
## derived from it.  `_wrapfloor` stays per-site because it is the one
## entry whose implementation genuinely differs between backends.
##
## Every value here serves EVERY backend: either pure Python arithmetic,
## valid on numpy scalars, numpy arrays and jax tracers alike, or the
## dispatching pair `_maxc_numpy`/`_minc_numpy`.  One implementation per
## primitive, so a dict entry cannot disagree with the `_imp_` that
## overrides it -- `_wrapfloor` is the single exception, and it is
## per-site because `numpy.floor` and `jnp.floor` are genuinely
## different functions with no `_imp_` to reconcile them.
## ----------------------------------------------------------------------

_KERNEL_NUMPY = {
    '_recip2': _recip2_numpy,
    '_rdiv': _rdiv_numpy,
    '_step': _step_numpy,
    'maxc': _maxc_numpy,
    'minc': _minc_numpy,
}


def _kernel_jax(jnp):
    """`_KERNEL_NUMPY` plus the one entry that differs on jax.

    `_wrapfloor` has no `_imp_`, so its binding comes from this map and
    has to name the right module's `floor`.  Everything else in the
    kernel dispatches on its argument and needs no swap.
    """
    return dict(_KERNEL_NUMPY, _wrapfloor=jnp.floor)



#: What ``simparam`` can genuinely answer.  A name absent here is not a
#: gap to be filled with a plausible number -- the LRM's contract is that
#: the caller's default is used, and a model that supplied one has already
#: said what it wants when the simulator has no opinion.
_SIMPARAMS = {
    ## pycircuit has no standing gmin.  Its gmin is a CONTINUATION
    ## schedule inside the DC solver (`dcanalysis.GminSteppingNewton`,
    ## `JunctionGminSteppingNewton`), ramped to zero before the answer is
    ## returned -- not a conductance models should shunt themselves with.
    ## Answering 0.0 is therefore the truth, and it is also what every
    ## caller in PSP103 asks for as its default.
    'gmin': 0.0,
    ## No die shrink or layout scaling layer exists; 1.0 is not a
    ## placeholder, it is the identity these would multiply by.
    'scale': 1.0,
    'shrink': 1.0,
    'sourceScaleFactor': 1.0,
}


def simparam(name, default=None):
    """Verilog-A's ``$simparam(name, default)`` -- ask the simulator.

    Resolved at COMPILE time, because every parameter pycircuit can answer
    is fixed for the run; nothing here varies per iteration. A name the
    simulator has no opinion on returns ``default``, per the LRM. Asking
    for an unknown name with no default is an error, not a silent zero --
    the model would be reading something that does not exist.

    ``$temperature`` is not routed through here: it is a separate LRM
    construct and varies per analysis, so it stays the `TEMP` symbol that
    `vt()` uses.

    PSP103 calls this exactly once, ``$simparam("gmin", 0.0)``
    (`PSP103_module.include:642`), which is the case this exists for.
    """
    if name in _SIMPARAMS:
        return sympy.sympify(_SIMPARAMS[name])
    if default is None:
        raise ValueError(
            'simparam(%r) is not a simulator parameter pycircuit supplies '
            '(it knows %s) and no default was given. Verilog-A leaves an '
            'unknown name with no default undefined; supply the value the '
            'model should use when the simulator has no opinion, as '
            'simparam(%r, <value>).'
            % (name, ', '.join(sorted(_SIMPARAMS)), name))
    return sympy.sympify(default)


def discontinuity(degree=0):
    """Verilog-A's ``$discontinuity(degree)``: tell the solver the model
    has a discontinuity in its ``degree``-th derivative here.

    **Parsed and ignored, deliberately.**  It is advisory: it lets a
    simulator drop integration order rather than discover the corner by
    rejection, and pycircuit does that from breakpoints instead (see
    ``Cross``).  gnucap accepts it and generates an empty body for the
    same reason.  Accepting it costs nothing and lets a model written for
    another simulator compile unchanged -- 41 call sites in the vacask
    device library -- whereas rejecting it would fail those models for a
    hint they can do without.  If it ever becomes load-bearing, the
    honest implementation is to schedule a breakpoint, which is what
    ``Cross`` already does.
    """
    return sympy.Integer(0)


class _Laplace(sympy.Function):
    """Marker for a laplace_* application; see :func:`laplace_nd`."""
    nargs = (3,)


def laplace_nd(expr, num, den):
    """Verilog-A's ``laplace_nd(expr, num[], den[])``: apply the rational
    transfer function ``H(s) = sum(num[k] s^k) / sum(den[k] s^k)``.

    Coefficients ascend in powers of ``s``, as the LRM specifies.  The
    filter is realised as **state equations**, not as a convolution: an
    order-``N`` denominator introduces ``N`` unknowns in controllable
    canonical form, so the simulator integrates them with its own method
    and order and the Jacobian is exact through them.  That is the same
    construction ``idt`` uses, generalised.

    ``H`` must be proper (``len(num) <= len(den)``): an improper transfer
    function differentiates its input, which is ``ddt``'s job and needs a
    charge, not a state.
    """
    num = [sympy.sympify(c) for c in num]
    den = [sympy.sympify(c) for c in den]
    if len(den) < 2:
        raise ValueError('laplace_nd needs a denominator of order >= 1; '
                         'a constant gain is just multiplication')
    if len(num) > len(den):
        raise NotImplementedError(
            'laplace_nd requires a proper transfer function '
            '(len(num) <= len(den)); an improper one differentiates its '
            'input, which is what ddt is for')
    return _Laplace(sympy.sympify(expr), sympy.Tuple(*num),
                    sympy.Tuple(*den))


def laplace_zp(expr, zeros, poles):
    """``laplace_zp(expr, zeros[], poles[])`` -- zeros and poles instead of
    coefficients.  Each is a flat list of (real, imaginary) pairs, as the
    LRM specifies, and the pair ``(r, i)`` contributes a factor
    ``1 - s/(r + j i)`` (or ``1 - s/r`` when ``i`` is zero).

    Converted to :func:`laplace_nd` by expanding the products, so the
    realisation and everything downstream is identical.
    """
    sv = sympy.Symbol('_lap_s')

    def poly_of(pairs):
        expr_ = sympy.Integer(1)
        it = list(pairs)
        for k in range(0, len(it), 2):
            re_, im_ = sympy.sympify(it[k]), sympy.sympify(it[k + 1])
            ## `.is_zero`, NOT `== 0`: sympy's `==` is STRUCTURAL, so
            ## `Float(0.0) == 0` is False and a real pole written with a
            ## float imaginary part took the conjugate-pair branch --
            ## silently producing a filter of twice the intended order.
            if im_.is_zero:
                if re_.is_zero:
                    expr_ *= sv
                else:
                    expr_ *= (1 - sv / re_)
            else:
                ## A conjugate pair contributes a real quadratic.
                mag2 = re_ ** 2 + im_ ** 2
                expr_ *= (1 - 2 * re_ * sv / mag2 + sv ** 2 / mag2)
        return sympy.Poly(sympy.expand(expr_), sv).all_coeffs()[::-1]

    return laplace_nd(expr, poly_of(zeros), poly_of(poles))


class _Limit(sympy.Function):
    """Base marker for a limited probe; see :func:`limit_pnj`.

    Never instantiated itself -- one subclass per limiter KIND, so that
    `atoms(_Limit)` still collects every limited probe in one pass while
    each kind carries its own argument list.  `kind` is the string the
    generated `limit()` dispatches on.

    **Each kind has TWO arities**: the probe plus that kind's own
    parameters, and the same again with a three-integer GROUP TAG
    appended -- `(gid, slot, sequential)` -- which is what
    :func:`limit_together` writes onto a marker to say "this probe is
    limited with the others in group `gid`, and it is the `slot`-th of
    them".  The tag rides in `args` rather than in a side table because
    the markers are sympy expressions embedded in the contribution: a
    registry keyed on the marker object would be keyed on sympy's cache,
    which merges structurally equal instances.
    """
    kind = None


class _LimitPnj(_Limit):
    """`$limit(probe, "pnjlim", IS, VT)`; see :func:`limit_pnj`."""
    nargs = (3, 6)
    kind = 'pnj'


class _LimitFet(_Limit):
    """`$limit(probe, "fetlim", vto)`; see :func:`limit_fet`."""
    nargs = (2, 5)
    kind = 'fet'


class _LimitVds(_Limit):
    """`$limit(probe, "limvds")`; see :func:`limit_vds`."""
    nargs = (1, 4)
    kind = 'vds'


class _LimitId(_Limit):
    """An identity probe -- a limited unknown with NO law; see
    :func:`limit_identity`."""
    nargs = (1, 4)
    kind = 'id'


## How many arguments after the probe belong to the LAW, per kind.  What
## follows them, if anything, is the group tag.
_LIMIT_NPAR = {'pnj': 2, 'fet': 1, 'vds': 0, 'id': 0}

#: The SPICE name of each kind's law, as `explain()` prints it, and the
#: DSL function that declares it.  One table, because the three-entry
#: dict literal used to be spelled at four sites and a fourth kind would
#: have had to be added to each.
_LIMIT_LAW = {'pnj': 'pnjlim', 'fet': 'fetlim', 'vds': 'limvds',
              'id': 'identity (no law)'}
_LIMIT_FN = {'pnj': 'limit_pnj', 'fet': 'limit_fet', 'vds': 'limit_vds',
             'id': 'limit_identity'}

## Group ids only ever have to be unique WITHIN one `analog()` body, but a
## process-wide counter is what makes them so without threading state
## through the DSL -- and it keeps two classes defined in one module from
## sharing a tag through sympy's expression cache.
_LIMIT_GROUP_IDS = itertools.count(1)


def _limit_parts(app):
    """Split a `_Limit` marker into `(probe, params, group)`.

    `group` is `(gid, slot, sequential)` or None.  One place that knows
    the tag's layout, because `generate_code` and `limit_together` both
    have to agree about it.
    """
    npar = _LIMIT_NPAR[app.kind]
    tail = app.args[1 + npar:]
    grp = (int(tail[0]), int(tail[1]), bool(int(tail[2]))) if tail else None
    return app.args[0], tuple(app.args[1:1 + npar]), grp


def limit_together(*probes, sequential=False):
    """Limit several of a device's probes AS ONE, roadmap 10.3(b).

    Takes the ordinary per-probe declarations and returns them tagged as
    one group, in the order given::

        vgs, vds = limit_together(limit_fet(bgs.V, VTO), limit_vds(bds.V))

    The per-probe forms are unchanged and keep working on their own; this
    is an envelope around them, not a replacement, so a model adopts it by
    wrapping one line.

    **What the grouping buys is the WRITE-BACK.**  A per-probe limiter
    moves one endpoint per probe, so two probes sharing a terminal have to
    negotiate for it and one of them is pushed onto a node it would not
    have chosen (roadmap 12.1).  A group is written back as a whole: the
    probes are a graph over the device's terminals, the node that drifted
    LEAST from the last accepted point is held, and every other node is
    derived from it, so **every probe carries exactly its own limited
    value and none is undone** -- see
    :func:`pycircuit.circuit._limiting.device_writeback`.  A probe that
    did not bite is still a constraint: the device holds `vds` where
    Newton put it while it compresses `vgs`, instead of letting `vds`
    follow the source node the other probe moved.

    ``sequential=True`` asks for **SPICE's coupling** instead of
    independent limiting.  `mos1load.c` limits `vgs`, then recomputes
    `vds = vgs_lim - vgd` from the UNLIMITED `vgd` before calling
    `limvds`, so its `vds` is shifted by exactly the compression the gate
    took.  Here that is each probe in turn moving its MINUS terminal by
    the correction it just applied, and the next probe reading the shifted
    vector -- which is the same thing, since holding `vgd` while `vgs`
    moves IS moving the source.  It is deliberately ORDER DEPENDENT; that
    is what "in this order" means.

    **Measured, and it is not an improvement here.**  On the cascode of
    `test_device_limiter.py`, over 48 operating points, independent
    grouping and SPICE's sequencing cost the same to within a couple of
    iterations, and at the reference point (20 V, 2 V, 0.8 V) sequencing
    is worse.  It is offered because it is SPICE's answer and a model
    that wants bit-comparability with SPICE's ITERATES needs it, not
    because it converges better; the default is independent.
    """
    if len(probes) < 2:
        raise ValueError('limit_together needs at least two probes; one '
                         'probe on its own is already a per-probe limiter')
    gid = next(_LIMIT_GROUP_IDS)
    out = []
    for slot, app in enumerate(probes):
        if not isinstance(app, _Limit):
            raise ValueError('limit_together takes limit_pnj/limit_fet/'
                             'limit_vds declarations; got %r' % (app,))
        if app.kind == 'id':
            ## An identity probe never bites, and the grouped write-back
            ## treats a probe that did not bite as a CONSTRAINT it holds
            ## while the others move (`device_writeback`).  Holding a
            ## branch that asked for nothing would be a law it does not
            ## have.  Declare it beside the group, not in it.
            raise ValueError('limit_identity has no law and cannot be '
                             'grouped; declare it beside the group')
        if _limit_parts(app)[2] is not None:
            raise ValueError('this probe is already in a $limit group; a '
                             'probe carries one limiter and one group')
        out.append(type(app)(*app.args, sympy.Integer(gid),
                             sympy.Integer(slot),
                             sympy.Integer(1 if sequential else 0)))
    return tuple(out)


def _limit_probe(probe, who):
    """Validate a limiter's probe argument, shared by the four kinds."""
    if not (isinstance(probe, Quantity) and probe.isbranch
            and probe.quantity == 'V'):
        raise ValueError('%s limits a branch potential (b.V); '
                         'got %r' % (who, probe))
    return probe


def limit_pnj(probe, IS, VT):
    """Verilog-A's ``$limit(probe, "pnjlim", ...)``: ask the simulator to
    bound this branch voltage's per-iteration excursion.

    ``probe`` must be a branch potential (``b.V``); ``IS`` and ``VT`` are
    the junction's saturation current and thermal voltage, which is what
    SPICE's ``pnjlim`` needs to place its critical voltage.  The
    expression evaluates to the probe itself, so it is written inline::

        Contribution(b.I, IS * (limexp(limit_pnj(b.V, IS, vt()) / vt()) - 1))

    The generated element implements the STATE-FREE limiter convention
    (``SubCircuit.limit``): it returns a limited copy of the solution
    sub-vector rather than keeping a private linearisation point.  That
    is what the codebase prefers and what a traced backend could ever
    support -- and it means the limiter moves the *solution*, so the
    residual and Jacobian are never taken at different points, which is
    the inconsistency classical limiting is criticised for.

    Prefer PCNR (``pcnr=True``) where the element qualifies: it makes the
    limited quantity an explicit unknown owned by this device, so two
    devices sharing a node cannot fight over it.  ``limit_pnj`` is the
    fallback for models PCNR cannot take, and for ordinary runs, where
    ``pcnr=False`` is the default.

    ``IS`` and ``VT`` MAY be `var()` symbols (2026-08-25, roadmap 12.4)::

        isT = var(area * IS * expl(factlog), 'isT')
        vbe = var(limit_pnj(b.V, isT, nf * vt()), 'vbe')

    They are resolved against the let-chain, not inlined, so a limiter
    parameter written in terms of a deep chain costs what the chain
    costs.  Before this they were lambdified over the parameter namespace
    alone, where an intermediate symbol is not in scope: `lambdify`
    then returns a function that hands BACK THE SYMBOL, and the first
    Newton iteration raises ``TypeError: Cannot convert expression to
    float`` from inside `limit()`.  That is why the Gummel-Poon BJT
    spells its temperature-scaled saturation current twice.

    What they may NOT do is depend on the solution: a limiter's
    parameters are evaluated from the card, before the device is.  A
    chain that reaches the iterate is refused at compile time, by name.
    """
    _limit_probe(probe, 'limit_pnj')
    return _LimitPnj(probe, sympy.sympify(IS), sympy.sympify(VT))


def limit_fet(probe, vto):
    """Verilog-A's ``$limit(probe, "fetlim", vto)``: bound a FET gate
    voltage's per-iteration excursion.

    ``probe`` must be a branch potential (``b.V``) -- SPICE's ``vgs``, or
    ``vgd`` for a device run in reverse.  ``vto`` is the voltage the law
    measures the excursion against; see :func:`pycircuit.circuit._limiting._fetlim`
    for the four regimes.  The expression evaluates to the probe itself,
    so it is written inline exactly as :func:`limit_pnj` is::

        vgs = limit_fet(bgs.V, VTO)

    **``vto`` may depend on the bias, and is then evaluated at the LAST
    ACCEPTED iterate** -- which is SPICE's semantics: ``mos1load.c``
    passes ``fetlim`` a ``von`` recomputed from the previous iterate's
    bulk bias, not a card constant.  So a model with a body effect
    writes its real turn-on::

        von = vtoT + gamma * (sqrt(phiT + vsb) - sqrt(phiT))
        vgs = limit_fet(bgs.V, von)

    and the limiter is placed where the device actually turns on.  A
    parameter that reads only the card is evaluated once per parameter
    set, exactly as before.

    **History, kept because both halves of it were wrong in turn
    (2026-08-25).**  This used to say a bias-dependent ``vto`` "is not
    expressible: a limiter runs BEFORE the device is evaluated, so its
    parameters cannot read the iterate".  The order was right and the
    conclusion did not follow -- a limiter is handed ``x0`` precisely so
    it can measure against the last accepted point, and a parameter
    evaluated THERE is as well-defined as ``vold`` is.  Before that, the
    reason was given as a lambdify namespace, which was the mechanism
    and not the rule.  The cost of the refusal was measured on the EKV
    card at 2 V of body bias: true turn-on 1.06 V against ``vto`` =
    0.50, every clamp 565 mV low.  Looser, not wrong -- and now gone.

    Ordering: see :func:`limit_vds`.
    """
    _limit_probe(probe, 'limit_fet')
    return _LimitFet(probe, sympy.sympify(vto))


def limit_vds(probe):
    """Verilog-A's ``$limit(probe, "limvds")``: bound a FET's ``vds``.

    ``probe`` must be a branch potential (``b.V``).  Takes no parameters
    -- SPICE's ``limvds`` is a bare piecewise clamp with hard-coded
    breakpoints at 3.5 V and 4 V; see
    :func:`pycircuit.circuit._limiting._limvds`.

    **Ordering, and where this differs from SPICE.**  The generated
    ``limit()`` walks the declared probes in a canonical order -- largest
    correction first -- each reading the partially-limited vector, and
    writes each probe's limited value by moving ONE terminal: whichever
    end drifted further from the last accepted point, or the other end if
    an earlier probe already moved that one (the rule generalises the BJT
    case ``limit_junctions`` handles with its ``move`` field -- two
    junctions sharing a base would otherwise have the second write undo
    the first).  The order is derived from the data, not from the
    declaration, so **writing the two calls the other way round changes
    nothing**.

    **What is NOT true, and used to be claimed here (corrected
    2026-08-25):** that each probe therefore ends up carrying exactly its
    own limited value.  The old wording reasoned from the compile-time
    write-back, where ``limit_fet(V(g,s))`` moved ``g`` and
    ``limit_vds(V(d,s))`` moved ``d``.  Now that the terminal is chosen
    at run time both probes can want the shared source, and the one that
    gets it moves a node the other probe's branch hangs off -- so the
    earlier probe's branch follows and lands somewhere its own law never
    chose.  Measured over 813 random steps in which both probes bite: 27
    of them (3.3%) end with a branch its own law would still move, the
    worst by 36 V.  A per-probe write-back applies each correction as a
    DISPLACEMENT of one node, which is not the same thing as satisfying a
    branch constraint.  :func:`limit_together` is what satisfies them all
    at once; the numbers are in `test_device_limiter.py`.

    SPICE's order DOES matter, and this is the difference.
    ``mos1load.c`` limits ``vgs`` first and then recomputes
    ``vds = vgs - vgd`` from the *unlimited* ``vgd`` before calling
    ``limvds`` -- so its ``vds`` is shifted by exactly the amount the
    gate was compressed, ``delta = vgs_lim - vgs``, and its ``vgd`` is
    preserved instead.  Here ``limvds`` sees the unshifted ``vds`` and it
    is ``vgd`` that moves.  Both are limited points and both leave the
    solution untouched -- a limiter only chooses where the next Jacobian
    is taken -- but they are not the same point, and no per-probe
    limiter can produce SPICE's: the coupling runs through a THIRD
    branch (``vgd``) that neither probe names.

    **2026-08-25:** :func:`limit_together` now expresses it --
    ``limit_together(limit_fet(bgs.V, VTO), limit_vds(bds.V),
    sequential=True)`` is exactly `mos1load.c`'s sequence.  What has NOT
    changed is which one to prefer: measured over 48 operating points of
    the `test_device_limiter.py` cascode, SPICE's ordering is worse than
    the independent one (1029 Jacobian evaluations against 927), and at
    the reference point it is the same.  The paragraph above described
    the gap as a limitation; it is better read as a difference, and the
    default is still not SPICE's.
    """
    _limit_probe(probe, 'limit_vds')
    return _LimitVds(probe)


def limit_identity(probe):
    """A ``$limit`` probe with NO law: the branch potential is declared
    as a limited unknown and left exactly where Newton put it.

    Written inline like the other kinds and evaluating to the probe::

        vsb = limit_identity(bsb.V)

    **On the ordinary path it does nothing, and the nothing is exact.**
    `apply_limit` returns ``vnew`` unchanged, the generated ``limit()``
    sees ``vlim == vnew`` and writes nothing -- the same "a limiter that
    did not bite must touch nothing" rule the real kinds obey, so the
    convergence signal "did limiting fire?" is untouched by declaring
    one.  It may not be put in a `limit_together` group, where a probe
    that did not bite is held as a constraint; an identity probe would
    then be holding a branch it has no opinion about.

    **Under vector PCNR it is what admits a device whose current reads a
    branch no other probe covers.**  The declared-probe route
    (`_pcnr_declared_route`) qualifies a device only when every
    resistive current reaches the solution through its probes, and it
    is a structural rule, not a convergence one: a MOSFET reads
    ``vgs``, ``vds`` AND ``vsb``, and `limit_fet`/`limit_vds` cover two
    of the three.  ``limit_identity(bsb.V)`` gives the third its own
    PCNR unknown, with the identity as its law.  The device's block is
    then ``m = 3`` over four rows, each unknown owned by the device, and
    the tail-node clash of roadmap sec. 15 cannot form on it any more
    than on the other two.  The probe is NOT a junction: it never
    appears in the gmin pair view (`pcnr.pcnr_junction_pairs` lists
    ``'pnj'`` probes only), so no ladder puts a conductance across it.

    `EkvNmosHdl` is the motivating case (``var(vsb) reads b, g, which no
    $limit probe limits``), and its measurement is in
    ``test_limit_identity.py``.
    """
    _limit_probe(probe, 'limit_identity')
    return _LimitId(probe)


class _IntermediateSymbol(sympy.Symbol):
    """A named intermediate; see :func:`var`."""

    def __new__(cls, name):
        return super().__new__(cls, name, real=True)

    ## `Symbol` pickles as `cls(name, **assumptions)`, which this
    ## one-argument `__new__` refuses; the assumption is fixed, so the
    ## name alone rebuilds it.  Needed by the compile cache, where a
    ## chained model's PCNR and limiter records carry its intermediates.
    def __getnewargs_ex__(self):
        return ((self.name,), {})


#: Intermediates declared by the model body currently being compiled.
#: Compilation happens once per class, at class-creation time, on one
#: thread, so a module-level list is the whole mechanism.
_VAR_STACK = []


def var(expr, name=None):
    """Name an intermediate quantity, and return a symbol standing for it.

    **This is what makes a real compact model compilable.**  Without it,
    every reference to an intermediate substitutes its whole definition,
    so an expression that mentions a previous result twice doubles the
    tree at each step: a model nested `n` deep has `2**n` tree
    occurrences, and every sympy traversal -- `subs`, `diff`, printing --
    walks all of them.  Measured on a MOSFET-shaped chain: depth 3
    compiled in a second, depth 8 did not finish.

    With `var`, the definition is recorded once and referenced by symbol,
    so the compiler emits a **let-chain** -- assignments in dependency
    order -- and differentiates it by forward accumulation, which is
    linear in the number of definitions rather than exponential in their
    depth.  The same chain at depth 256 compiles in seconds.

    Use it for anything a Verilog-A model would have written to a local
    variable::

        sp = var(surface_potential(vgs, vsb), 'sp')
        qi = var(inversion_charge(sp), 'qi')
        return Contribution(b.I, W / L * mu * qi * ...)

    Naming is optional and only affects generated-code readability.
    """
    if not _VAR_STACK:
        raise RuntimeError(
            'var() may only be called while an element is being compiled '
            '(inside analog())')
    reg = _VAR_STACK[-1]
    sym = _IntermediateSymbol('_v%d_%s' % (len(reg), name or 'i'))
    reg.append((sym, sympy.sympify(expr)))
    return sym


def _param_symbol(name):
    """The symbol standing for the parameter called `name`.

    `Symbol(name)` for every name that can be spelled in generated code.
    A keyword-named parameter (SPICE's ``lambda``, ``as``; reachable
    only through `ParamNamespace`) gets a mangled name, because the
    let-chain printer emits ``def _f(x, lambda, ...)`` VERBATIM -- a
    SyntaxError -- and every consumer that rebuilds the symbol from the
    parameter name (four JAX sites) must agree on the spelling.
    """
    if name.isidentifier() and not keyword.iskeyword(name):
        return sympy.Symbol(name)
    return sympy.Symbol('_hdl_kw_' + re.sub(r'\W', '_', name))


def param_given(name):
    """Verilog-A's ``$param_given(name)``: 1.0 when the instance was given
    that parameter explicitly, 0.0 when it fell back to its default.

    A *runtime* value, not a compile-time one, because the element is
    compiled once per class while givenness is a property of each
    instance.  Use it in a ``Piecewise`` exactly as Verilog-A uses it in
    an ``if``::

        rs_eff = sympy.Piecewise((rs, param_given('rs') > 0.5), (0.0, True))

    By call count this is the most-used system function in real compact
    models, which is why it is here before more glamorous operators.
    """
    return sympy.Symbol('_hdl_given_%s' % name, real=True)


class ParamNamespace(object):
    """The parameters of one class, as an object handed to ``analog()``.

    Opted into per class with ``params_as = 'p'`` (any identifier), which
    makes the FIRST argument of ``analog()`` this object rather than a
    terminal::

        class Res(Behavioural):
            params_as = 'p'
            instparams = [Parameter(name='r', ...),
                          Parameter(name='lambda', ...)]

            @staticmethod
            def analog(p, plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, b.V / (p.r * (1 + p['lambda'])))

    ``p.r`` and ``p['r']`` are the SAME `sympy.Symbol('r')` the bare-name
    style binds -- this is a compile-time object whose attributes are
    the parameter symbols, so the generated code, the compilation record
    and `explain` are identical between the two spellings.  What it
    changes is scope and spelling:

    * a helper shared by several classes receives ``p`` as an ordinary
      argument and reads ``p.bf``; the bare-name style cannot reach a
      helper at all, since a helper has its own ``__globals__``
      (`elements_hdl._gp_core` used to be rebound by hand for this);
    * a parameter may carry ANY name, including a Python keyword --
      SPICE's ``lambda`` and ``as`` -- reached as ``p['lambda']``;
    * a linter sees an attribute access, not an undefined global, so
      the ``# noqa: F821`` the bare style needs disappears.

    ``p.given('rs')`` is `param_given` with the name checked against
    the declaration; ``p.names`` is the declared order.  A name that
    is not declared raises at compile time, naming the class and the
    declared names, exactly as an undeclared bare name does.

    In a ``params_as`` class the bare names are NOT injected: the two
    styles coexist across classes, not inside one.
    """
    __slots__ = ('_cls', '_names', '_syms')

    def __init__(self, clsname, names, syms):
        object.__setattr__(self, '_cls', clsname)
        object.__setattr__(self, '_names', tuple(names))
        object.__setattr__(self, '_syms', dict(zip(names, syms)))

    @property
    def names(self):
        """The declared parameter names, in declaration order."""
        return self._names

    def _missing(self, name):
        return ('%s.analog() reads parameter %r, which %s does not '
                'declare (its instparams are %s). Declare it: '
                'instparams = [..., Parameter(name=%r, desc=..., unit=..., '
                'default=...)].'
                % (self._cls, name, self._cls,
                   ', '.join(repr(a) for a in self._names)
                   or '(none declared)', name))

    def __getattr__(self, name):
        ## `__slots__` fields and `names` resolve before this is reached,
        ## so only parameter lookups arrive here.  Dunder probes
        ## (`__array__`, `_sympy_`, `__iter__`...) from numpy and sympy
        ## must see a plain AttributeError, and do.
        try:
            return self._syms[name]
        except KeyError:
            if name.startswith('__') or name.startswith('_'):
                raise AttributeError(name)
            raise AttributeError(self._missing(name)) from None

    def __setattr__(self, name, value):
        raise AttributeError(
            '%s.analog(): the parameter namespace is read-only; %r cannot '
            'be assigned. Parameters are declared in instparams and their '
            'values come from the instance.' % (self._cls, name))

    def __getitem__(self, name):
        try:
            return self._syms[name]
        except KeyError:
            raise KeyError(self._missing(name)) from None

    def __contains__(self, name):
        return name in self._syms

    def __iter__(self):
        return iter(self._names)

    def __len__(self):
        return len(self._names)

    def given(self, name):
        """`param_given`, with the name checked against the declaration:
        a typo in the bare ``param_given('rz')`` compiles silently and
        reads as never-given at every instance."""
        if name not in self._syms:
            raise KeyError(self._missing(name))
        return param_given(name)

    def __repr__(self):
        return '<parameters of %s: %s>' % (self._cls, ', '.join(self._names))


def ac_stim(mag=1.0, phase=0.0):
    """Verilog-A's ``ac_stim``: a small-signal excitation, live in AC
    analysis and identically zero everywhere else.

    Contributed like any other term -- ``Contribution(b.I, ac_stim(1.0))``
    -- and routed to an AC-only source vector, so a behavioural element
    can BE a small-signal stimulus rather than merely transparent to one.
    """
    return _AcStim(sympy.sympify(mag), sympy.sympify(phase))


class _AcStim(sympy.Function):
    nargs = (2,)


def ddx(expr, probe):
    """Verilog-A ``ddx(expr, probe)``: the symbolic partial derivative of
    ``expr`` with respect to a ``V``/``I`` probe, other probes held fixed.

    Pure sympy -- ``Quantity`` is an atom, so this is a substitution and a
    ``diff``.  The probe must appear in ``expr`` as the SAME quantity
    object (``b.V`` matched against ``b.V``, not against an expanded
    node-voltage difference).
    """
    d = sympy.Dummy()
    return sympy.sympify(expr).subs(probe, d).diff(d).subs(d, probe)


class Quantity(sympy.AtomicExpr):
    """Reference to a voltage or current of a branch or node.

    A sympy atom, so it takes part in symbolic expressions (arithmetic,
    ``.subs()``, ``Matrix.jacobian()``, ...).
    """
    ## A voltage or a current is REAL.  Without this assumption sympy
    ## refuses to build `Quantity < 80` and every Piecewise over an
    ## operating point -- limexp included -- dies with "Invalid comparison
    ## of non-real".
    is_real = True
    is_commutative = True

    def __new__(cls, quantity, branch_or_node):
        if quantity not in ('V', 'I'):
            raise ValueError("quantity must be either 'V' or 'I'")
        if not isinstance(branch_or_node, (Node, Branch)):
            raise ValueError('branch_or_node must be a Branch or Node object')
        if quantity == 'I' and isinstance(branch_or_node, Node):
            raise ValueError('Current can only be taken on branches')

        obj = super().__new__(cls)
        obj.quantity = quantity
        obj.branch_or_node = branch_or_node
        return obj

    def _hashable_content(self):
        if isinstance(self.branch_or_node, Branch):
            ## The branch NAME is part of the identity.  Without it,
            ## `I(Branch(p,m,'a'))` and `I(Branch(p,m,'b'))` are the same
            ## sympy atom, so both resolve to whichever branch's current
            ## symbol was substituted last -- a wrong stamp, not a
            ## missing one.
            key = ('branch', self.branch_or_node.plus.name,
                   self.branch_or_node.minus.name,
                   getattr(self.branch_or_node, 'name', None))
        else:
            key = ('node', self.branch_or_node.name)
        return (self.quantity,) + key

    @property
    def name(self):
        ## sympy's StrPrinter dispatches on the CLASS NAME and hits
        ## _print_Quantity (written for sympy.physics.units.Quantity), which
        ## reads .name -- without this property, printing/sorting any Mul
        ## holding a Quantity next to a Function dies at class creation.
        return repr(self)

    @property
    def isnode(self):
        return isinstance(self.branch_or_node, Node)

    @property
    def isbranch(self):
        return isinstance(self.branch_or_node, Branch)

    def __repr__(self):
        if self.isbranch:
            nm = getattr(self.branch_or_node, 'name', None)
            return '%s(%s,%s%s)' % (
                self.quantity, self.branch_or_node.plus.name,
                self.branch_or_node.minus.name,
                '' if nm is None else ':%s' % nm)
        return '%s(%s)' % (self.quantity, self.branch_or_node.name)

    __str__ = __repr__


## ----------------------------------------------------------------------
## select() -- the compiler's own both-arms domain propagation
## ----------------------------------------------------------------------


class SelectRefused(ValueError):
    """`select` was handed a condition whose shape it cannot clamp.

    Raised at BUILD time, from `select`, naming the condition and the
    sympy class that made it unsupported.  It is a refusal, not a
    failure: the model can always write the `Piecewise` and clamp by
    hand, and `hdl.unclamped` says so in one word.
    """


class _Unclamped(object):
    """The value `unclamped` returns; see there."""

    __slots__ = ('expr',)

    def __init__(self, expr):
        self.expr = sympy.sympify(expr)


def unclamped(expr):
    """Mark one arm of a `select` as deliberately NOT clamped.

    The escape hatch.  An arm wrapped in this is passed to `Piecewise`
    exactly as written, and its condition still constrains the LATER
    arms.  Use it when the arm is finite over the whole real line
    anyway, when the clamp would cost more than it buys, or when
    `select` refuses a shape it turns out you do not need clamped::

        hdl.select((hdl.unclamped(1.0 + t * v), t < 0), (f(v), True))
    """
    return _Unclamped(expr)


#: The atoms `select` will substitute for a clamped copy: a parameter or
#: intermediate symbol, or a branch/node quantity.  Anything else on one
#: side of a relational is a compound, and a compound cannot be clamped
#: by substitution -- clamping `a*b` would mean choosing which factor to
#: move.
_CLAMPABLE = (sympy.Symbol, Quantity)


def _select_atoms(expr):
    """Every atom in `expr` that `select` is willing to clamp."""
    out = set()
    for a in expr.atoms(*_CLAMPABLE):
        out.add(a)
    return out


def _abs_arg(e):
    """`x` if `e` is ``|x|`` or ``x*x`` or ``x**2`` for a clampable `x`."""
    if isinstance(e, (Abs, sympy.Abs)) and isinstance(e.args[0], _CLAMPABLE):
        return e.args[0], False
    if isinstance(e, sympy.Pow) and e.exp == 2 and \
            isinstance(e.base, _CLAMPABLE):
        return e.base, True
    if isinstance(e, sympy.Mul) and len(e.args) == 2 and \
            e.args[0] == e.args[1] and isinstance(e.args[0], _CLAMPABLE):
        return e.args[0], True
    return None, False


def _refuse(cond, why):
    raise SelectRefused(
        'select cannot clamp the condition `%s`: %s.  Write the arm as '
        'hdl.unclamped(...) if it needs no clamp, clamp it by hand in a '
        'sympy.Piecewise, or pass strict=False to accept an unclamped '
        'arm silently.' % (cond, why))


def _select_constraints(cond, negate, out, bad):
    """Append `(atom, side, bound, is_strict)` constraints for `cond`.

    `side` is ``'lo'`` (the atom is bounded BELOW by `bound`),
    ``'hi'``, or ``'out'`` (``|atom| >= bound``).  With `negate` the
    constraints describe the complement, which is what the arms AFTER
    this one are selected in.

    Nothing is raised here.  Every fragment that could not be turned
    into a box is appended to `bad` as ``(fragment, why)``, and
    `select` decides whether the silence matters -- which it does
    exactly when the ARM uses one of the fragment's symbols.  A refusal
    over a condition an arm does not read would be noise, and noise is
    how a strict mode gets turned off.
    """
    if cond is True or cond is sympy.true:
        ## `True` constrains nothing, and its complement is empty --
        ## which also constrains nothing, since no arm is selected there.
        return
    if cond is False or cond is sympy.false:
        return
    if isinstance(cond, sympy.Not):
        return _select_constraints(cond.args[0], not negate, out, bad)
    if isinstance(cond, (sympy.And, sympy.Or)):
        ## An And is an intersection of boxes and composes; an Or is a
        ## UNION, and the union of two boxes is not a box, so no
        ## substitution can be the identity on all of it.  De Morgan
        ## turns each into the other under negation.
        conj = isinstance(cond, sympy.And) != bool(negate)
        if not conj:
            bad.append((cond, 'it is a %s in a position where the region '
                        'is a UNION of boxes, and no min/max substitution '
                        'is the identity on a union'
                        % ('disjunction' if isinstance(cond, sympy.Or)
                           else 'negated conjunction'),
                        _select_atoms(cond)))
            return
        for a in cond.args:
            _select_constraints(a, negate, out, bad)
        return
    if isinstance(cond, (sympy.Equality, sympy.Unequality)):
        bad.append((cond, 'an equality constrains a measure-zero set (or '
                    'its complement), which no clamp can reach',
                    _select_atoms(cond)))
        return
    if not isinstance(cond, sympy.core.relational.Relational):
        bad.append((cond, 'it is a %s, not a relational'
                    % type(cond).__name__, _select_atoms(cond)))
        return

    ## Normalise to `lhs < rhs` or `lhs <= rhs`, folding the negation in.
    ## `not (a < b)` is `a >= b`, i.e. `b <= a` -- the strictness flips
    ## with the direction, which is the whole reason the two are tracked
    ## separately.
    lhs, rhs = cond.lhs, cond.rhs
    isstrict = isinstance(cond, (sympy.StrictLessThan,
                                 sympy.StrictGreaterThan))
    if isinstance(cond, (sympy.StrictGreaterThan, sympy.GreaterThan)):
        lhs, rhs = rhs, lhs
    if negate:
        lhs, rhs = rhs, lhs
        isstrict = not isstrict
    ## Now the region is `lhs < rhs` (or `<=`).
    hit = False
    if isinstance(lhs, _CLAMPABLE):
        out.append((lhs, 'hi', rhs, isstrict))
        hit = True
    if isinstance(rhs, _CLAMPABLE):
        out.append((rhs, 'lo', lhs, isstrict))
        hit = True
    ## `|x| < c` and its algebraic spelling `x*x < c`: an INTERVAL
    ## around zero, which is a box and composes like one.  `c < |x|` is
    ## the complement of one -- a union of two half-lines, not a box,
    ## but still reachable by a substitution because the hole is
    ## bounded; see `_punch`.
    ##
    ## Tried on BOTH sides, and after the atom rule rather than instead
    ## of it, so that `|x| < c` with `c` a symbol clamps `c` from below
    ## AND `x` into the interval.  A square whose bound is not a literal
    ## is reported even when the other side DID yield a clamp: `x*x < c`
    ## clamps `c`, which is not what an arm that divides by `x` needed.
    for e, other, side in ((lhs, rhs, 'in'), (rhs, lhs, 'out')):
        a, squared = _abs_arg(e)
        if a is None:
            continue
        bound = other
        if squared:
            ## `x*x < c` bounds `|x|` by `sqrt(c)`, and a symbolic `c`
            ## would need `sqrt` of a quantity whose sign nothing here
            ## knows.  A non-negative literal is the case that occurs.
            if not (bound.is_number and bound.is_real and bound >= 0):
                bad.append((cond,
                            'the square of `%s` is bounded by `%s`, which '
                            'is not a non-negative literal, so the bound '
                            'on `%s` itself is not available without '
                            'assuming a sign' % (a, bound, a), {a}))
                continue
            bound = sympy.sqrt(bound)
        if side == 'in':
            out.append((a, 'lo', -bound, isstrict))
            out.append((a, 'hi', bound, isstrict))
        else:
            out.append((a, 'out', bound, isstrict))
        hit = True
    if hit:
        return
    bad.append((cond, 'neither side is a symbol or a branch quantity '
                '(`%s` is a %s, `%s` is a %s), so there is no atom to '
                'substitute for a clamped copy.  Clamping the COMPOUND '
                'is not offered, because a substitution cannot reliably '
                'find one: sympy stores `1.0/(a*b)` as '
                '`1.0*a**-1*b**-1`, in which `a*b` is not a node, so '
                'the clamp would silently do nothing.  Bind the '
                'compound with var() first and the condition has an '
                'atom'
                % (lhs, type(lhs).__name__, rhs, type(rhs).__name__),
                _select_atoms(cond)))


def _punch(v, c):
    """``v`` pushed OUT of ``(-c, c)``, the identity where ``|v| >= c``.

    The complement of an interval is not a box, so `minc`/`maxc` cannot
    express it; this can, in four operations::

        v + pm(v) * maxc(0, c - |v|)

    with ``pm(v) = 2*_step(v, 0) - 1``, which is ``+-1`` and never 0.
    Where ``|v| >= c`` the `maxc` is exactly ``0.0`` and the sum returns
    `v` bit for bit; inside the hole the result is ``+-c``.

    The argument ORDER of the `maxc` is load-bearing.  `maxc(0, u)`
    differentiates to ``(1 - _step(0, u)) * du``, which is exactly zero
    on ``u <= 0`` -- the TIE included, i.e. at ``|v| = c``, which is a
    selected point whenever the arm's condition is non-strict.  Written
    `maxc(u, 0)` the tie would go the other way and the derivative at
    ``|v| = c`` would come out 0 instead of 1.
    """
    pm = 2.0 * _step(v, sympy.Integer(0)) - 1.0
    return v + pm * maxc(sympy.Integer(0), c - Abs(v))


def _clamped_symbol(atom, lo, hi, out, margin):
    """The expression `atom` is replaced by inside one arm."""
    rep = atom
    if lo is not None:
        b, isstrict = lo
        rep = maxc(rep, b + margin if isstrict else b)
    if hi is not None:
        b, isstrict = hi
        rep = minc(rep, b - margin if isstrict else b)
    if out is not None:
        b, isstrict = out
        rep = _punch(rep, b + margin if isstrict else b)
    return rep


def _var_dependencies():
    """symbol -> the symbols its `var()` definition transitively uses.

    Empty outside a compile.  Used only by the shadow check below, so it
    is built on demand and thrown away.
    """
    if not _VAR_STACK:
        return {}
    deps = {}
    for sym, expr in _VAR_STACK[-1]:
        direct = _select_atoms(expr)
        acc = set(direct)
        for d in direct:
            acc |= deps.get(d, set())
        deps[sym] = acc
    return deps


def select(*arms, **kw):
    """A `Piecewise` that clamps each arm to its OWN condition's domain.

    Both arms of a compiled conditional are evaluated -- the selection
    happens afterwards -- so an arm that is never chosen at this bias
    still runs, still raises floating-point flags, and still has its
    derivative taken.  The rule the model has to obey is therefore
    "every arm must be finite, and have a finite derivative, EVERYWHERE",
    not merely where it is selected, and the way that rule is obeyed by
    hand is a clamped copy of the input fed to the arm that does not
    want it::

        sympy.Piecewise((0.0, p <= 0.0), (1.0 / maxc(p, 1e-30), True))

    That `maxc` is domain propagation the compiler could have done: the
    second arm is selected exactly where ``p > 0``, so inside it `p` may
    be replaced by any expression that equals `p` there.  `select` does
    the substitution::

        hdl.select((0.0, p <= 0.0), (1.0 / p, True), margin=1e-30)

    Each arm is `Piecewise`'s ``(expr, cond)`` pair, last condition
    ``True``.  For every arm, the region it is selected in -- its own
    condition AND the complement of every earlier one -- is turned into
    a box of bounds on individual symbols, and each bounded symbol is
    replaced INSIDE THAT ARM by a `minc`/`maxc` clamped copy.

    **The clamp is the identity where the arm is selected, in value and
    in derivative, exactly.**  `minc(v, c)` returns `v` bit for bit
    wherever ``v <= c``, and differentiates to ``_step(c, v) = 1``
    there, with exactly zero weight on `c` -- so a clamped arm and a
    hand-written one agree to the last bit of the residual AND of the
    Jacobian.  The tie ``v == c`` is included on purpose: it is a
    selected point when the condition is non-strict.

    Supported conditions, and nothing else:

    ==========================  ==================================
    ``v < c``, ``v <= c``       ``v -> minc(v, c)``
    ``v > c``, ``v >= c``       ``v -> maxc(v, c)``
    ``a < v`` etc.              the mirror of the above
    ``u < v`` (both atoms)      both, simultaneously
    ``|v| < c``, ``v*v < c``    ``v -> minc(maxc(v, -c), c)``
    ``|v| > c``                 pushed out of the hole; see `_punch`
    ``And(...)``                the intersection of the boxes
    ``Or(...)`` in a NEGATED
    position (an earlier arm)   ditto, by De Morgan
    ``True``                    the complement of the earlier arms
    ==========================  ==================================

    `c` may be any expression, including one that mentions `v` itself
    (``vds < vdsat``): the substitution is simultaneous, so the bound is
    evaluated at the UNCLAMPED `v` and no cycle is created.  It costs
    nothing in finiteness either, because the bound was already being
    evaluated -- by the condition.

    Everything else is REFUSED, loudly, with `SelectRefused` naming the
    condition and the reason: an `Or` in a positive position (a union of
    boxes is not a box), an equality, a relational whose sides are both
    compound (``a*b < c`` constrains neither factor).  A `select` that
    quietly did not clamp would be worse than the `Piecewise` it
    replaced, because the author would trust it.  `strict=False` turns
    the refusals into silence, and `hdl.unclamped(expr)` marks ONE arm
    as deliberately raw.

    ``margin`` (default ``0.0``) is subtracted from every bound that
    came from a STRICT inequality, on the side that keeps it inside the
    selected region.  It exists because clamping to the boundary is not
    always enough: ``p > 0`` clamps to ``maxc(p, 0)``, and an arm that
    divides by `p` still divides by zero.  A margin makes the clamp
    ``maxc(p, 1e-30)``, which is finite -- at the price that the clamp
    is no longer the identity on the sliver ``0 < p < margin``.  That
    is a real trade and it is why the default is zero: with
    ``margin=0`` a selected value can never change.  A margin is only
    ever applied to a strict bound, because a non-strict one has its
    boundary INSIDE the selected region, where moving it would change
    an answer.

    The margin is an absolute number in the units of whatever it bounds,
    so write it relative to that scale (``margin=1e-9*leff``) rather
    than picking a constant that happens to work on one card.

    **What it cannot do**, and refuses rather than pretending: a clamp
    can only reach a symbol the arm MENTIONS.  If the arm uses an
    intermediate that `var()` bound earlier from the same symbol --
    ``vdsc = leff*vmax/us`` used in an arm conditioned on ``vmax > 0``
    -- substituting `vmax` inside the arm does not change `vdsc`, and
    the clamp is a no-op.  `select` detects exactly that case (the arm
    depends on a constrained symbol ONLY through a `var()` intermediate)
    and refuses.
    """
    margin = kw.pop('margin', 0.0)
    strict = kw.pop('strict', True)
    if kw:
        raise TypeError('select() got unexpected keyword arguments %s'
                        % sorted(kw))
    if not arms:
        raise TypeError('select() needs at least one arm')
    margin = sympy.sympify(margin)

    pairs = []
    for arm in arms:
        if not (isinstance(arm, (tuple, list)) and len(arm) == 2):
            raise TypeError(
                'select() takes (expr, cond) pairs, as sympy.Piecewise '
                'does; got %r' % (arm,))
        pairs.append(arm)

    deps = None
    built = []
    for k, (expr, cond) in enumerate(pairs):
        raw = isinstance(expr, _Unclamped)
        if raw:
            expr = expr.expr
        expr = sympy.sympify(expr)
        cond = sympy.sympify(cond)
        if raw:
            built.append((expr, cond))
            continue
        cons, bad = [], []
        ## The arm's own condition, then the EARLIER ones negated.  An
        ## earlier condition is a refinement: dropping it only widens
        ## the region the clamp has to be the identity on, which is
        ## safe.  For a `True` arm it is the OTHER way round -- the
        ## complement of the earlier conditions is the only thing that
        ## can clamp it -- so a fragment lost there is reported too.
        _select_constraints(cond, False, cons, bad)
        isdefault = cond is sympy.true
        for j in range(k):
            _select_constraints(sympy.sympify(pairs[j][1]), True, cons,
                                bad if isdefault else [])
        bysym = {}
        for atom, side, bound, isstrict in cons:
            slot = bysym.setdefault(atom, {'lo': None, 'hi': None,
                                           'out': None})
            have = slot[side]
            if have is None:
                slot[side] = (bound, isstrict)
            else:
                ## Two bounds on one symbol from a conjunction: keep the
                ## TIGHTER, which is the one that makes the box smaller.
                ## `maxc`/`minc` fold two literals at build time, so the
                ## common case costs nothing at run time.
                b, s = have
                comb = (maxc(b, bound) if side == 'lo' else
                        minc(b, bound) if side == 'hi' else
                        maxc(b, bound))
                slot[side] = (comb, s or isstrict)
        rep = {}
        atoms = _select_atoms(expr)
        ## A condition fragment that is not a box leaves the arm
        ## unclamped in the symbols it mentions.  That is a silent
        ## no-op only if the arm USES one of them; when it does, this
        ## is the refusal the whole design exists for.
        if strict:
            for frag, why, risk in bad:
                miss = sorted(str(a) for a in risk & atoms)
                if miss:
                    _refuse(frag, '%s -- and this arm uses %s'
                            % (why, ', '.join(miss)))
        for atom, slot in bysym.items():
            if slot['out'] is not None and (slot['lo'] is not None or
                                            slot['hi'] is not None):
                if strict:
                    _refuse(cond, 'symbol `%s` is constrained both to a '
                            'box and out of a hole, and the composition '
                            'of the two is not the identity on the '
                            'region' % atom)
                continue
            if atom not in atoms:
                ## The arm does not mention it.  Either it genuinely
                ## does not use it -- fine, nothing to clamp -- or it
                ## reaches it through a `var()` intermediate, where a
                ## substitution cannot follow.  The second is a silent
                ## no-op and is refused.
                if deps is None:
                    deps = _var_dependencies()
                shadow = [s for s in atoms
                          if isinstance(s, _IntermediateSymbol)
                          and atom in deps.get(s, ())]
                if shadow and strict:
                    raise SelectRefused(
                        'select cannot clamp `%s` in the arm selected by '
                        '`%s`: the arm does not mention it, it reaches it '
                        'through the var() intermediate %s, and a '
                        'substitution cannot follow a symbol into a '
                        'definition that was already bound.  Clamp that '
                        'intermediate where it is DEFINED, inline it into '
                        'the arm, or mark the arm hdl.unclamped().'
                        % (atom, cond,
                           ', '.join(sorted(str(s) for s in shadow))))
                continue
            rep[atom] = _clamped_symbol(atom, slot['lo'], slot['hi'],
                                        slot['out'], margin)
        built.append((expr.xreplace(rep) if rep else expr, cond))
    return sympy.Piecewise(*built)

class Statement:
    pass


class Cross(Statement):
    """Verilog-A's ``@(cross(expr, direction))``: ask the transient solver
    to place a timepoint where ``expr`` crosses zero.

    ``direction`` is ``+1`` (rising), ``-1`` (falling) or ``0`` (either),
    as in the LRM.  Returned alongside contributions::

        return (Contribution(b_out.V, ...),
                Cross(b_in.V - vref, +1))

    The generated element predicts the crossing by linear extrapolation
    from the last two ACCEPTED points and publishes it through
    ``next_event``, which is the contract ``SubCircuit.next_event`` polls
    and the same first-order prediction gnucap uses.  A prediction only
    has to BRACKET the corner: the solver's break-step machinery lands on
    it and drops the integration order across it.

    Without this the controller finds a comparator's edge by rejecting
    steps until the local error is small enough -- correct, but it smears
    the corner and costs a rejection storm at every transition.
    """

    def __init__(self, expr, direction=0):
        if direction not in (-1, 0, 1):
            raise ValueError('cross direction must be -1, 0 or +1')
        self.expr = sympy.sympify(expr)
        self.direction = direction

    def nodes(self):
        nodes = set()
        for atom in self.expr.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch:
                    nodes.add(atom.branch_or_node.plus)
                    nodes.add(atom.branch_or_node.minus)
                else:
                    nodes.add(atom.branch_or_node)
        return nodes


class Collapse(Statement):
    """Collapse a branch's nodes when a PARAMETER condition holds.

    Verilog-A's optional-parasitic idiom, which PSP103 uses seven times::

        if ((R) > 0.0) I(N1,N2) <+ (G) * V(N1,N2);
        else           V(N1,N2) <+ 0.0;

    is written here as the contribution plus a `Collapse` saying when the
    other arm applies::

        br = Branch(n1, n2, 'rd')
        return (Contribution(br.I, g * br.V),
                Collapse(br, rd <= 0))

    When the condition holds for an instance's parameters, that branch's
    nodes become one and every contribution to it is dropped -- so a
    ``g = 1/rd`` that would be infinite is never compiled, let alone
    evaluated.  When it does not, the `Collapse` is simply absent and the
    contributions stand.

    The condition may mention **parameters only**, never a node voltage
    or branch current.  That restriction is what makes this tractable: a
    collapse decided by the operating point would change the matrix
    sparsity every Newton iteration, whereas a parameter-driven one is
    fixed the moment the instance is built.  Each distinct combination of
    conditions is compiled once and cached on the class, and instances
    are retargeted to the variant matching their parameters.

    The cost of NOT doing this is measured in
    ``benchmarks/collapse_cost.py``: carrying PSP's seven parasitics as
    branches instead is 14x slower on a 100-device DC solve.
    """

    def __init__(self, branch, when):
        if not isinstance(branch, Branch):
            raise TypeError('Collapse takes a Branch, got %r' % (branch,))
        self.branch = branch
        self.when = sympy.sympify(when)
        bad = [a for a in self.when.atoms() if isinstance(a, Quantity)]
        if bad:
            raise ValueError(
                'Collapse condition may use parameters only, but mentions '
                '%s. A collapse that depends on the operating point would '
                'change the matrix sparsity every Newton iteration; keep '
                'the branch and use a Piecewise in the contribution '
                'instead.' % ', '.join(sorted(repr(b) for b in bad)))

    def nodes(self):
        return {self.branch.plus, self.branch.minus}


class Contribution(Statement):
    """``Contribution(b.I, expr)`` or ``Contribution(b.V, expr)`` --
    the DSL form of Verilog-A's ``<+``."""

    def __init__(self, lhs, rhs):
        if not isinstance(lhs, Quantity) or not lhs.isbranch:
            raise ValueError(
                'the contribution target must be the I or V of a Branch '
                '(got %r); node contributions are not part of the language'
                % (lhs,))
        self.lhs = lhs
        self.rhs = sympy.sympify(rhs)

    def nodes(self):
        """The set of node objects referred to in lhs and rhs."""
        nodes = set()
        for atom in self.lhs.atoms() | self.rhs.atoms():
            if isinstance(atom, Quantity):
                if atom.isbranch:
                    nodes.add(atom.branch_or_node.plus)
                    nodes.add(atom.branch_or_node.minus)
                else:
                    nodes.add(atom.branch_or_node)
        return nodes


def _quantity_free(expr):
    """True when nothing in ``expr`` refers to a circuit quantity or state."""
    return not any(isinstance(a, Quantity) for a in expr.atoms()) and \
        not expr.atoms(_StateSymbol)


class _StateSymbol(sympy.Symbol):
    """Symbol standing for an idt/idtmod state unknown."""

    def __new__(cls, name):
        ## real=True for the same reason as Quantity: a Piecewise over a
        ## state must be constructible.
        return super().__new__(cls, name, real=True)

    def __getnewargs_ex__(self):        # see _IntermediateSymbol
        return ((self.name,), {})


def _split_terms(rhs, xset, tset=frozenset()):
    """Route an expression's additive terms into (i, q, u) parts.

    Called AFTER ``resolve()``, so circuit quantities and states appear as
    the x-symbols in ``xset``; a term's kind is decided by which of those
    it touches.  ``ddt(arg)`` terms -- bare, or scaled by an x-independent
    factor -- become charge; an x-DEPENDENT scaling of ``ddt`` is refused
    loudly (``g(v)*ddt(v)`` is not the derivative of any charge; move the
    factor inside the ``ddt``).  x-free terms become sources (they may
    contain ``TIME``); the rest is static current.
    """
    def xfree(e):
        return not (e.free_symbols & xset)

    ## NOT `sympy.expand(rhs)`.  Expansion distributes every product of
    ## sums, and a compact model is nothing but nested products of sums:
    ## measured on a MOSFET-shaped expression, expansion multiplied the
    ## operation count by ~6 PER NESTING LEVEL (12 -> 22 -> 109 -> 294 ->
    ## 828 ...), so a real device model never finished compiling.
    ##
    ## Expansion was only ever needed to separate additive terms of
    ## different KINDS (charge from current from source from noise), and
    ## the top-level `Add` already does that: a model writes
    ## `I <+ f(v) + ddt(q(v))`, not a product that has to be multiplied
    ## out before the `ddt` becomes visible.  A `ddt` buried inside a
    ## product is refused anyway (it is not the derivative of a charge),
    ## so nothing is lost -- and the one thing expansion did buy, putting
    ## `q` in a canonical form, was never required.
    terms = rhs.args if rhs.is_Add else (rhs,)

    iterms, qterms, uterms, nterms, acterms = [], [], [], [], []
    for term in terms:
        acs = term.atoms(_AcStim)
        if acs:
            if isinstance(term, _AcStim):
                app, coeff = term, sympy.Integer(1)
            elif term.is_Mul:
                apps = [f for f in term.args if isinstance(f, _AcStim)]
                rest = [f for f in term.args if not isinstance(f, _AcStim)]
                coeff = sympy.Mul(*rest)
                if len(apps) != 1 or not xfree(coeff):
                    raise NotImplementedError(
                        'ac_stim may appear as a term or scaled by an '
                        'x-independent factor; %r is neither' % (term,))
                app = apps[0]
            else:
                raise NotImplementedError(
                    'ac_stim may only appear additively: %r' % (term,))
            mag, phase = app.args
            acterms.append(coeff * mag * sympy.exp(sympy.I * phase
                                                   * sympy.pi / 180))
            continue
        noises = term.atoms(white_noise) | term.atoms(flicker_noise)
        if noises:
            if isinstance(term, (white_noise, flicker_noise)):
                app, coeff = term, sympy.Integer(1)
            elif term.is_Mul:
                napps = [f for f in term.args
                         if isinstance(f, (white_noise, flicker_noise))]
                rest = [f for f in term.args
                        if not isinstance(f, (white_noise, flicker_noise))]
                coeff = sympy.Mul(*rest)
                ## The scale factor MAY depend on x, unlike `ddt`'s.
                ## `ddt`'s restriction is a conservation argument --
                ## `g(v)*ddt(v)` is not the derivative of any charge.  No
                ## such argument applies here: `CY` is a function of the
                ## operating point by construction (`CY(x, w)`), and the
                ## power INSIDE the noise call has always been allowed to
                ## depend on x, so forbidding it outside was a rule with
                ## nothing behind it.
                if len(napps) != 1:
                    raise NotImplementedError(
                        'a term may scale at most ONE noise function; %r '
                        'has %d.  Two correlated sources are written as two '
                        'contributions sharing a name.'
                        % (term, len(napps)))
                app = napps[0]
            else:
                raise NotImplementedError(
                    'noise functions may only appear additively: %r' % (term,))
            if isinstance(app, white_noise):
                pwr = app.args[0]
            else:
                pwr = app.args[0] / FREQ ** app.args[1]
            ## `coeff` scales the AMPLITUDE, so it enters the PSD squared.
            ## Kept factored out rather than pre-squared: a named group
            ## needs the amplitude itself, to build the CROSS terms.
            nterms.append((app.noise_name, coeff, pwr))
            continue
        ddts = term.atoms(ddt)
        if ddts:
            if isinstance(term, ddt):
                qterms.append(term.args[0])
                continue
            if term.is_Mul:
                dd = [f for f in term.args if isinstance(f, ddt)]
                rest = [f for f in term.args if not isinstance(f, ddt)]
                coeff = sympy.Mul(*rest)
                if len(dd) == 1 and xfree(coeff) and not coeff.atoms(ddt):
                    qterms.append(coeff * dd[0].args[0])
                    continue
            raise NotImplementedError(
                'ddt may appear as a term or scaled by an x-independent '
                'factor; %r is neither. A state-dependent factor is not '
                'the derivative of any charge -- move it inside the ddt.'
                % (term,))
        if xfree(term):
            ## `u` carries only genuinely TIME-VARYING excitation; a
            ## constant belongs to the device's static characteristic
            ## (a diode's -IS), which keeps generated stamps identical to
            ## the hand-written ones and keeps `u` meaning one thing.
            ## `tset` names the intermediates that carry time.  Without
            ## it a source written through `var()` would be an opaque
            ## symbol, `has(TIME)` would say no, and the drive would be
            ## filed as a static characteristic -- a dead source.
            timed = term.has(TIME) or bool(term.free_symbols & tset)
            (uterms if timed else iterms).append(term)
        else:
            iterms.append(term)
    return (sympy.Add(*iterms), sympy.Add(*qterms), sympy.Add(*uterms),
            nterms, sympy.Add(*acterms))


## ----------------------------------------------------------------------
## The syntactic surface.
##
## The semantic checks below (switch branches, state-scaled `ddt`,
## collapse conditions, instparams drift) were always here; what a
## newcomer actually reaches was not.  Three of the mistakes in this
## section used to be SILENT -- an argument named `self` became a
## terminal, a class-body `terminals` was discarded, an internal
## `Node('out')` merged into the `out` pin -- and a silent wrong answer
## is the expensive failure, not a bad message.  Each check names the
## class, says what is wrong, and says what to write instead.
## ----------------------------------------------------------------------

def _analog_function(func, paramnames, paramsyms):
    """`func` rebound to a PRIVATE copy of its globals, with the
    parameter symbols added.

    The bare-name convention -- `analog()` writing `r` for the parameter
    named `r` -- is load-bearing (63 `# noqa: F821` in `compact.py`
    alone), so the names still have to resolve as globals.  What they
    must NOT do is stay there: this used to be
    `analogfunc.__globals__.update(...)`, and since `copy()` of a
    function returns the same object, `__globals__` IS the defining
    module's dict.  Two measured consequences: a parameter named `vt`
    replaced `hdl.vt`, the FUNCTION, in `elements_hdl` for good; and a
    class could reference a parameter that a DIFFERENT class had
    declared and compile silently, failing later as a bare `NameError`
    from inside `<lambdifygenerated-12>`.

    A snapshot copy is exact here because the function is called once,
    immediately, by `generate_code`: every module global it can legally
    see already exists.  Writes from inside `analog()` land in the copy,
    and `globals()` read from inside it still shows the parameters.
    """
    ns = dict(func.__globals__)
    ns.update(zip(paramnames, paramsyms))
    out = types.FunctionType(func.__code__, ns, func.__name__,
                             func.__defaults__, func.__closure__)
    out.__kwdefaults__ = func.__kwdefaults__
    out.__dict__.update(func.__dict__)
    return out


def _check_analog_declaration(cls):
    """Signature-level mistakes, checked before `analog()` is run."""
    name = cls.__name__
    func = cls.__dict__.get('analog', cls.analog)
    if isinstance(func, classmethod):
        raise TypeError(
            '%s.analog() is a classmethod. It must be a staticmethod: its '
            'argument names ARE the element terminals, so a `cls` first '
            'argument would become a terminal called "cls".' % name)
    if isinstance(func, staticmethod):
        func = func.__func__
    spec = inspect.getfullargspec(func)
    args = list(spec.args)
    pas = _params_as(cls)
    if pas is not None:
        if not args or args[0] != pas:
            raise TypeError(
                '%s sets params_as = %r, so analog()\'s FIRST argument must '
                'be named %r and receives the parameter namespace; it is '
                'analog(%s). Write "def analog(%s, %s)", or remove '
                'params_as to use the bare-name style.'
                % (name, pas, pas, ', '.join(args), pas,
                   ', '.join(args) or 'plus, minus'))
        args = args[1:]
    if args[:1] in (['self'], ['cls']):
        rest = args[1:] or ['plus', 'minus']
        raise TypeError(
            '%s.analog() takes %r as its first argument, so %r would become '
            'the element\'s FIRST TERMINAL and every connection after it '
            'would shift by one pin. analog() is a staticmethod whose '
            'argument names are the terminals: write "@staticmethod" above '
            'it and "def analog(%s)".'
            % (name, args[0], args[0], ', '.join(rest)))
    bad = []
    if spec.varargs:
        bad.append('*%s' % spec.varargs)
    if spec.varkw:
        bad.append('**%s' % spec.varkw)
    if spec.kwonlyargs:
        bad.append(', '.join('%s=' % a for a in spec.kwonlyargs))
    if bad:
        raise TypeError(
            '%s.analog() declares %s, which cannot become terminals: a '
            'terminal is one named positional argument, and the pin order '
            'is the argument order. Write the terminals out.'
            % (name, ' and '.join(bad)))
    if not args:
        raise TypeError(
            '%s.analog() declares no %sarguments, so the element has no '
            'terminals and nothing can be connected to it. Name the '
            'terminals as arguments: "def analog(%splus, minus)".'
            % (name, '' if pas is None else 'terminal ',
               '' if pas is None else pas + ', '))
    ## A class-body `terminals` is DISCARDED -- the metaclass overwrites
    ## it with analog()'s argument names.  Refused only when the two
    ## DISAGREE, which is the case where the discarding changes the
    ## element: a declaration that merely restates the signature is
    ## redundant, not wrong, and several existing models carry one.
    declared = cls.__dict__.get('terminals')
    if declared is not None and list(declared) != args:
        raise TypeError(
            '%s declares terminals = %r in the class body, which disagrees '
            'with analog(%s). The declaration is DISCARDED -- an HDL '
            'element\'s pin order is analog()\'s argument order -- so the '
            'element would have been built with %r and every connection '
            'silently placed by that list instead. Delete the declaration '
            'and rename or reorder analog()\'s arguments.'
            % (name, tuple(declared), ', '.join(spec.args),
               tuple(args)))
    return args


def _params_as(cls):
    """The class's ``params_as`` declaration, validated: `None` for the
    bare-name style, else the identifier naming analog()'s first
    argument."""
    pas = getattr(cls, 'params_as', None)
    if pas is None:
        return None
    if not isinstance(pas, str) or not pas.isidentifier() \
            or keyword.iskeyword(pas):
        raise TypeError(
            '%s sets params_as = %r, which is not a Python identifier. It '
            'names analog()\'s first argument, the one that receives the '
            'parameter namespace: params_as = \'p\'.'
            % (cls.__name__, pas))
    return pas


def _run_analog(name, func, args, paramnames, params_as=None):
    """Call the model's `analog()`, translating the errors an
    ordinary-looking Python expression raises out of it: an undefined
    bare name, a comparison used as a condition, and a float-domain
    math function applied to a symbol."""
    try:
        return func(*args)
    except NameError as e:
        unknown = getattr(e, 'name', None)
        if unknown is None:
            m = re.search(r"name '([^']+)' is not defined", str(e))
            unknown = m.group(1) if m else None
        if unknown is None:
            raise
        if params_as is not None:
            raise NameError(
                '%s.analog() uses the name %r, which is not defined: it is '
                'not a global of the module it is written in, and %s sets '
                'params_as = %r, so its parameters (%s) are NOT bound as '
                'bare names -- they are read as %s.%s or %s[%r]. If %r is '
                'meant to be a parameter, declare it and read it through '
                '%s; if the name comes from a helper analog() calls, pass '
                '%s to the helper.'
                % (name, unknown, name, params_as,
                   ', '.join(repr(a) for a in paramnames)
                   or '(none declared)',
                   params_as, unknown if unknown.isidentifier() else 'name',
                   params_as, unknown, unknown, params_as,
                   params_as)) from None
        raise NameError(
            '%s.analog() uses the name %r, which is not defined: it is not '
            'a declared instance parameter of %s (those are %s), and not a '
            'global of the module it is written in. Parameter symbols are '
            'bound from THIS class\'s instparams only -- a parameter '
            'declared by another class used to be in scope here, so a model '
            'like this compiled silently and failed later as a bare '
            'NameError from inside generated code. Declare it: '
            'instparams = [..., Parameter(name=%r, desc=..., unit=..., '
            'default=...)]; if the name comes from a helper analog() calls, '
            'fix it there.'
            % (name, unknown, name,
               ', '.join(repr(a) for a in paramnames) or '(none declared)',
               unknown)) from None
    except TypeError as e:
        text = str(e)
        if 'truth value of Relational' in text:
            raise TypeError(
                '%s.analog(): a circuit quantity was used where Python '
                'wanted a bool (%s). Node voltages, branch currents and '
                'parameters are SYMBOLS -- an "if", "and", "or", "not" or '
                'a chained comparison over them has no answer at compile '
                'time, because the element is compiled once and the '
                'quantity is different at every operating point. Select '
                'between arms at the VALUE level instead: '
                'sympy.Piecewise((a, cond), (b, True)) for an exact switch, '
                'or the smooth kernel forms (maxc, minc, mne, mxe) where '
                'the Jacobian has to stay continuous. A condition on '
                'parameters alone belongs in Collapse(branch, when) or in '
                'param_given().' % (name, text)) from None
        if 'convert expression to float' in text or \
                ('loop of ufunc' in text and 'callable' in text):
            raise TypeError(
                '%s.analog(): a math or numpy function was applied to a '
                'symbolic quantity (%s). Those evaluate floats; the model '
                'builds an expression that is differentiated and compiled '
                'later. Use the sympy or hdl kernel equivalents -- '
                'sympy.exp (or hdl.expl / limexp), sympy.sqrt (or '
                'safe_sqrt), sympy.log (or safe_ln), sympy.tanh, hdl.Abs, '
                'hdl.sign -- which build expressions rather than numbers.'
                % (name, text)) from None
        raise


def _check_node_collisions(name, terminalnames, terminalnodes, built):
    """An internal `Node('out')` in an element with an `out` terminal.

    Nodes are identified by NAME everywhere downstream (`Node.__hash__`
    is `hash(self.name)`), so the internal node used to be absorbed into
    the terminal without a word and its equation added to the pin's.

    `built` is every `Node` CONSTRUCTED while `analog()` ran, from
    `_NODE_TRACE`; identity against the terminal objects (which this
    compiler made, before the trace opened) is then exact.  Reading the
    finished statements instead does not work: sympy's cache returns an
    equal `Quantity` atom built during a previous compilation, so the
    node object in the expression is not necessarily the one the model
    just wrote.
    """
    tnames = set(terminalnames)
    for nd in built:
        if nd.name in tnames and not any(nd is t for t in terminalnodes):
            raise ValueError(
                '%s.analog() builds its own Node(%r), but %r is already a '
                'terminal of this element (it is an argument of analog()). '
                'Nodes are identified by name, so the two would silently '
                'become ONE node and the internal equation would be added '
                'to the pin\'s. Rename the internal node, or use the '
                'terminal argument itself if that is what you meant.'
                % (name, nd.name, nd.name))


def generate_code(cls):
    """Compile the class's ``analog()`` into the element interface.

    Returns a dict consumed by :class:`BehaviouralMeta`; see the module
    docstring for the compilation model.
    """
    terminalnames = _check_analog_declaration(cls)
    terminalnodes = [Node(name) for name in terminalnames]
    params_as = _params_as(cls)

    ## Bind the parameter names as sympy symbols in a PRIVATE copy of the
    ## analog() globals, so the body refers to instance parameters by bare
    ## name; the compiled functions take them as trailing arguments,
    ## supplied from the RESOLVED self.iparv at call time.  The copy is
    ## what keeps the injection out of the defining module -- see
    ## `_analog_function`.
    ##
    ## A `params_as` class gets the SAME symbols through a
    ## `ParamNamespace` as analog()'s first argument, and no bare names
    ## at all -- which is what lets it declare a keyword-named parameter
    ## (the bare style cannot bind `lambda` as a name, so it is refused
    ## there rather than left to fail as a SyntaxError-shaped NameError).
    paramnames = [p.name for p in cls.instparams]
    paramsyms = [_param_symbol(name) for name in paramnames]
    if params_as is None:
        unreachable = [nm for nm in paramnames
                       if not nm.isidentifier() or keyword.iskeyword(nm)]
        if unreachable:
            raise TypeError(
                '%s declares the parameter%s %s, which cannot be read as '
                'bare name%s inside analog(): %s. Set params_as = \'p\' '
                'on the class and read %s as %s -- the namespace reaches '
                'any declared name.'
                % (cls.__name__, '' if len(unreachable) == 1 else 's',
                   ', '.join(repr(nm) for nm in unreachable),
                   '' if len(unreachable) == 1 else 's',
                   'a Python keyword is not a valid name'
                   if all(keyword.iskeyword(nm) for nm in unreachable)
                   else 'not a valid Python identifier',
                   'it' if len(unreachable) == 1 else 'them',
                   ', '.join('p[%r]' % nm for nm in unreachable)))
        analogfunc = _analog_function(cls.analog, paramnames, paramsyms)
        callargs = terminalnodes
    else:
        analogfunc = _analog_function(cls.analog, [], [])
        callargs = [ParamNamespace(cls.__name__, paramnames, paramsyms)] \
            + terminalnodes

    _VAR_STACK.append([])
    _NODE_TRACE.append([])
    try:
        statements = _run_analog(cls.__name__, analogfunc,
                                 callargs, paramnames, params_as)
    finally:
        intermediates = _VAR_STACK.pop()
        builtnodes = _NODE_TRACE.pop()
    if statements is None:
        raise TypeError(
            '%s.analog() returned None. Verilog-A has no "return", but this '
            'DSL is Python: analog() must RETURN its contribution '
            'statements -- "return Contribution(b.I, g * b.V)", or a tuple '
            'of statements. Written as bare expressions they are built and '
            'thrown away, and the element would have no equations at all.'
            % cls.__name__)
    returned = statements
    if isinstance(statements, Statement):
        statements = (statements,)
    else:
        try:
            statements = tuple(statements)
        except TypeError:
            statements = None
    bad = returned if statements is None else \
        next((st for st in statements if not isinstance(st, Statement)), None)
    if bad is not None:
        raise TypeError(
            '%s.analog() returned %r, which is not a contribution '
            'statement. An analog body returns Contribution(b.I, expr) or '
            'Contribution(b.V, expr) -- optionally with Cross(...) and '
            'Collapse(...) -- or a tuple of them. A bare expression says '
            'nothing about which branch it belongs to.'
            % (cls.__name__, bad))
    _check_node_collisions(cls.__name__, terminalnames, terminalnodes,
                           builtnodes)
    crossings = [st for st in statements if isinstance(st, Cross)]
    collapses = [st for st in statements if isinstance(st, Collapse)]
    statements = [st for st in statements
                  if not isinstance(st, (Cross, Collapse))]

    ## PARAMETER-DRIVEN COLLAPSE.  `collapse_mask` says, per `Collapse`
    ## declaration and in declaration order, whether THIS variant takes
    ## the collapsed arm.  The class compiles with an all-False mask (so
    ## a model that declares collapses still has a well-defined default
    ## topology and parameter list); each instance is retargeted to the
    ## variant its own parameters select.
    collapse_mask = tuple(getattr(cls, '_hdl_collapse_mask', None)
                          or (False,) * len(collapses))
    if len(collapse_mask) != len(collapses):
        raise ValueError('collapse mask %r does not match %d Collapse '
                         'declarations' % (collapse_mask, len(collapses)))
    _collapsed_keys = set()
    for taken, st in zip(collapse_mask, collapses):
        if not taken:
            continue
        br = st.branch
        _collapsed_keys.add((br.plus.name, br.minus.name,
                             getattr(br, 'name', None)))
        ## Feed the existing unconditional machinery: a taken collapse IS
        ## `V(a,b) <+ 0`, and every contribution to that branch goes away
        ## -- including one whose coefficient is 1/R and would be
        ## infinite, which is the whole point of the model's conditional.
        statements.append(Contribution(Quantity('V', br), sympy.Integer(0)))
    if _collapsed_keys:
        statements = [
            st for st in statements
            if not (isinstance(st.lhs.branch_or_node, Branch)
                    and (st.lhs.branch_or_node.plus.name,
                         st.lhs.branch_or_node.minus.name,
                         getattr(st.lhs.branch_or_node, 'name', None))
                    in _collapsed_keys
                    and not sympy.sympify(st.rhs).is_zero)]

    ## ------------------------------------------------------------------
    ## Pass 1 -- inventory: nodes, branches with current unknowns (every
    ## V-contributed or I-probed branch), and idt/idtmod states.
    nodes = set()
    vbranches = []      # branches that get a current unknown, in order
    ibranch_keys = set()

    def branch_key(branch):
        ## The NAME is part of the identity.  Verilog-A's
        ## `branch (a,b) br1, br2;` declares two DISTINCT branches across
        ## one node pair, each with its own current, and PSP relies on
        ## that.  Keying on the node pair alone silently merged them into
        ## one unknown -- the second declaration's constitutive relation
        ## simply vanished.  An unnamed branch keeps its old identity, so
        ## nothing that worked before changes.
        return (branch.plus.name, branch.minus.name,
                getattr(branch, 'name', None))

    for st in crossings:
        nodes.update(st.nodes())
    for st in statements:
        nodes.update(st.nodes())
        if st.lhs.quantity == 'V':
            key = branch_key(st.lhs.branch_or_node)
            if key not in ibranch_keys:
                ibranch_keys.add(key)
                vbranches.append(st.lhs.branch_or_node)
        for atom in st.rhs.atoms(Quantity):
            if atom.isbranch and atom.quantity == 'I':
                key = branch_key(atom.branch_or_node)
                if key not in ibranch_keys:
                    ibranch_keys.add(key)
                    vbranches.append(atom.branch_or_node)

    ## ------------------------------------------------------------------
    ## NODE COLLAPSE.  `V(a,b) <+ 0` is Verilog-A's way of shorting a
    ## branch, and the standard idiom for an optional series resistance:
    ##
    ##     if (rs) I(a, ia) <+ V(a, ia)/rs;  else V(a, ia) <+ 0;
    ##
    ## Rather than carry a zero-volt source (a real branch unknown and a
    ## row for it), the two nodes become one.  Restricted, as gnucap
    ## restricts it, to UNCONDITIONAL contributions -- a collapse that
    ## depends on the operating point would change the sparsity pattern
    ## per iteration, which is a different and much larger feature.
    ##
    ## Only an INTERNAL node can be absorbed: terminals belong to the
    ## parent circuit's node map and this element cannot merge them.  A
    ## `V <+ 0` between two terminals therefore stays a zero-volt source,
    ## which is exactly what it means and costs one unknown.
    collapse = {}
    kept_statements = []
    for st in statements:
        b = st.lhs.branch_or_node
        ## `.is_zero`, not `== 0`.  sympy's `==` is STRUCTURAL equality,
        ## and `Float(0.0) == 0` is False -- so `V(a,b) <+ 0.0`, which is
        ## how PSP's own CollapsableR macro spells it, silently did not
        ## collapse while `V(a,b) <+ 0` did.  Same trap that once doubled
        ## a filter's order by taking the complex-pole branch for a real
        ## pole.
        if st.lhs.quantity == 'V' and \
                sympy.sympify(st.rhs).is_zero:
            p_, m_ = b.plus, b.minus
            p_int = p_ not in terminalnodes
            m_int = m_ not in terminalnodes
            if m_int:
                collapse[m_.name] = p_.name
                continue
            if p_int:
                collapse[p_.name] = m_.name
                continue
        kept_statements.append(st)
    if collapse:
        ## Chase chains (a -> b -> terminal) to a fixed point.
        for _ in range(len(collapse) + 1):
            collapse = {k: collapse.get(v, v) for k, v in collapse.items()}
        statements = kept_statements
        alias = {}
        for nd in list(nodes):
            tgt = collapse.get(nd.name)
            if tgt is not None:
                alias[nd] = next(o for o in nodes if o.name == tgt)
        nodes = {alias.get(nd, nd) for nd in nodes}

        def _realias(expr):
            """Point every node/branch reference in `expr` at its
            surviving node.  Returns None when nothing changed."""
            subs_q = {}
            for atom in expr.atoms(Quantity):
                bn = atom.branch_or_node
                if isinstance(bn, Branch):
                    np_, nm_ = alias.get(bn.plus, bn.plus), \
                        alias.get(bn.minus, bn.minus)
                    if (np_, nm_) != (bn.plus, bn.minus):
                        ## Carry the NAME through the rewrite.  A collapse
                        ## renames nodes; it does not merge two named
                        ## branches into one, and dropping the name here
                        ## would do exactly that.
                        nb = Branch(np_, nm_, name=bn.name)
                        subs_q[atom] = Quantity(atom.quantity, nb)
                elif bn in alias:
                    subs_q[atom] = Quantity(atom.quantity, alias[bn])
            return expr.subs(subs_q) if subs_q else None

        ## THE INTERMEDIATES HAVE TO BE REWRITTEN TOO, and for a while
        ## they were not.  On the let-chain path a statement's RHS is
        ## mostly `var()` SYMBOLS, and the branch voltages that mention
        ## the collapsed node live in those symbols' DEFINITIONS rather
        ## than in the statement.  Rewriting only the statements left the
        ## definitions pointing at a node that no longer existed, and the
        ## printer emitted a bare `V` -- a `NameError` at call time, from
        ## a model that compiled without complaint.
        ##
        ## It went unnoticed because `Collapse` was built and tested
        ## against the flatten path, where there are no intermediates to
        ## miss.  Any model that uses both features hits it immediately.
        for k, (sym_, expr_) in enumerate(intermediates):
            new_ = _realias(expr_)
            if new_ is not None:
                intermediates[k] = (sym_, new_)

        ## Rewrite the surviving statements' node references.
        for st in statements:
            new_ = _realias(st.rhs)
            if new_ is not None:
                st.rhs = new_
            lb = st.lhs.branch_or_node
            np_, nm_ = alias.get(lb.plus, lb.plus), alias.get(lb.minus,
                                                              lb.minus)
            if (np_, nm_) != (lb.plus, lb.minus):
                st.lhs = Quantity(st.lhs.quantity,
                                  Branch(np_, nm_, name=lb.name))
        ## Re-derive the branch inventory: a collapsed branch is gone.
        ## Re-derive the branch inventory.  A branch is GONE only if the
        ## collapse brought its two ends onto the same node -- that is the
        ## `V <+ 0` branch itself.  Every other branch survives with its
        ## endpoints renamed; the earlier rule dropped any branch merely
        ## TOUCHING a collapsed node, which silently deleted the
        ## constitutive relation of whatever hung off the absorbed node.
        _renamed = [Branch(alias.get(br.plus, br.plus),
                           alias.get(br.minus, br.minus), name=br.name)
                    for br in vbranches]
        vbranches = [br for br in _renamed
                     if br.plus.name != br.minus.name]
        ibranch_keys = {branch_key(br) for br in vbranches}

    ## SWITCH BRANCHES are refused, per hdl.md sec. 9 phase D.  In
    ## Verilog-A a branch may be a potential source in one operating
    ## region and a flow source in another, the two contributions being
    ## mutually exclusive under a condition.  Implementing that means
    ## re-deciding the stamp -- and the sparsity pattern -- every Newton
    ## iteration, which is gnucap's whole `_pot` re-stamping machine.
    ##
    ## Accepting BOTH unconditionally, as this compiler silently did,
    ## produces a defined but different element: a voltage source with a
    ## conductance in parallel, which is not what the model says.  A
    ## plausible wrong answer is the expensive failure, so it is named.
    _v_branches = {branch_key(st.lhs.branch_or_node)
                   for st in statements if st.lhs.quantity == 'V'}
    _i_branches = {branch_key(st.lhs.branch_or_node)
                   for st in statements if st.lhs.quantity == 'I'}
    _both = _v_branches & _i_branches
    if _both:
        raise NotImplementedError(
            'branch %s has both V and I contributions (a Verilog-A "switch '
            'branch"). Selecting between them per operating point changes '
            'the matrix structure each iteration and is not implemented; '
            'accepting both unconditionally would give a voltage source '
            'with a conductance in parallel, which is not what the model '
            'says. Split the regions into separate elements, or use a '
            'Piecewise inside ONE contribution.'
            % ', '.join(_fmt_branch(k) for k in
                        sorted(_both, key=lambda k: (k[0], k[1], k[2] or ''))))

    internalnodes = sorted(nodes - set(terminalnodes), key=lambda n: n.name)

    ## States: each distinct idt/idtmod APPLICATION is one state.  Walk in
    ## a deterministic order.
    states = []          # (state_symbol, kind, args); kind idt/idtmod/lap
    state_subst = {}

    ## Intermediates FIRST, and this is not tidiness -- it is the same
    ## reasoning as `_strip_limits` below.  A `var()` definition is never
    ## visited by the statement loop, so before this a `var(idt(...))`
    ## carried its `idt` all the way into `_chain_compile` and died there
    ## with "Unsupported ... idt".  It was in one way WORSE than the
    ## $limit case: `resolve()` does apply `state_subst` to intermediates,
    ## so mentioning the same application in a statement as well made it
    ## work -- whether a model compiled depended on where ELSE the author
    ## had happened to write it.  Any real model is chained -- that is
    ## what `var()` is for -- so state operators were unusable by exactly
    ## the models they exist to serve.
    _state_srcs = ([_e for _s, _e in intermediates]
                   + [st.rhs for st in statements])

    ## LAPLACE -> STATE EQUATIONS.  For H(s) = N(s)/D(s) with
    ## D = d0 + d1 s + ... + dN s^N, controllable canonical form gives
    ##
    ##     z1' = z2,  z2' = z3,  ...,
    ##     zN' = (u - d0 z1 - d1 z2 - ... - d_{N-1} zN) / dN
    ##     y   = n0 z1 + n1 z2 + ...
    ##
    ## i.e. N first-order states, which is exactly what this compiler
    ## already emits for `idt` -- one per integrator in the chain.  The
    ## simulator then integrates them with its own method and order, and
    ## the Jacobian is exact through them, which a convolution-based
    ## implementation cannot offer.
    for _src in _state_srcs:
        for app in sorted(_src.atoms(_Laplace), key=sympy.default_sort_key):
            if app in state_subst:
                continue
            u_expr, num, den = app.args
            num = list(num)
            den = list(den)
            order = len(den) - 1
            zs = []
            for _k in range(order):
                sym = _StateSymbol('_state%d' % len(states))
                states.append((sym, 'lap', None))    # rhs filled in below
                zs.append(sym)
            ## z_k' = z_{k+1} for k < N; the last row carries the input.
            rhs = []
            for k in range(order - 1):
                rhs.append(zs[k + 1])
            last = u_expr
            for k in range(order):
                last = last - den[k] * zs[k]
            rhs.append(last / den[order])
            for k in range(order):
                idx = len(states) - order + k
                sym, _kind, _ = states[idx]
                states[idx] = (sym, 'lap', (rhs[k],))
            out = sympy.Integer(0)
            for k, nk in enumerate(num):
                if k < order:
                    out += nk * zs[k]
                elif nk != 0:
                    ## num of the same order as den: the direct term
                    ## n_N/d_N, plus its correction through the chain.
                    out += nk * rhs[order - 1]
            state_subst[app] = out
    for _src in _state_srcs:
        for func_cls, kind in ((idt, 'idt'), (idtmod, 'idtmod')):
            for app in sorted(_src.atoms(func_cls),
                              key=sympy.default_sort_key):
                if app in state_subst:
                    continue
                sym = _StateSymbol('_state%d' % len(states))
                states.append((sym, kind, app.args))
                if kind == 'idt':
                    state_subst[app] = sym
                else:
                    args = app.args
                    m = args[2] if len(args) > 2 else sympy.Integer(1)
                    o = args[3] if len(args) > 3 else sympy.Integer(0)
                    ## the floored wrap, idtmod.md sec. 1.1
                    state_subst[app] = \
                        sym - m * _wrapfloor((sym - o) / m)

    for _sym, _kind, args in states:
        if args is None:
            continue
        if any(a.atoms(idt) or a.atoms(idtmod) for a in args):
            raise NotImplementedError(
                'nested idt/idtmod applications are not supported yet')

    ## ------------------------------------------------------------------
    ## The unknown vector, in Circuit.update_node_map order:
    ## terminals, internal nodes, state nodes, branch currents.
    statenodenames = ['_state%d' % k for k in range(len(states))]
    xlabels = ([n.name for n in terminalnodes]
               + [n.name for n in internalnodes]
               + statenodenames
               + ['_i_br%d' % k for k in range(len(vbranches))])
    xsyms = [sympy.Symbol('_x%d' % k, real=True)
             for k in range(len(xlabels))]

    index_of = {}
    for k, node in enumerate(terminalnodes + internalnodes):
        index_of[('node', node.name)] = k
    off = len(terminalnodes) + len(internalnodes)
    for k, (sym, _kind, _args) in enumerate(states):
        index_of[('state', sym.name)] = off + k
    off += len(states)
    for k, br in enumerate(vbranches):
        index_of[('branch', branch_key(br))] = off + k

    n = len(xlabels)

    ## Substitutions: node voltages, branch-current probes, state symbols.
    subst = {}
    for node in terminalnodes + internalnodes:
        subst[Quantity('V', node)] = xsyms[index_of[('node', node.name)]]
    for br in vbranches:
        subst[Quantity('I', br)] = xsyms[index_of[('branch', branch_key(br))]]
    for sym, _kind, _args in states:
        subst[sym] = xsyms[index_of[('state', sym.name)]]

    def resolve(expr):
        ## STATES FIRST: the recorded idt/idtmod applications are in their
        ## original quantity form, and substituting branch voltages before
        ## them would rewrite the applications out from under the mapping
        ## (idt(V(a,b)) becoming idt(V(a)-V(b)), matching nothing).
        expr = expr.subs(state_subst)
        bsubs = {}
        for atom in expr.atoms(Quantity):
            if atom.isbranch and atom.quantity == 'V':
                b = atom.branch_or_node
                bsubs[atom] = Quantity('V', b.plus) - Quantity('V', b.minus)
            elif atom.isbranch and atom.quantity == 'I' and \
                    branch_key(atom.branch_or_node) not in ibranch_keys:
                raise ValueError(
                    'I%s is probed but that branch has no current '
                    'unknown; probe only V-contributed branches.'
                    % _fmt_branch(branch_key(atom.branch_or_node)))
        expr = expr.subs(bsubs)
        return expr.subs(subst)

    ## $limit: record (branch, KIND, that kind's parameter expressions) for
    ## each limited probe, then replace the marker by the probe itself -- it
    ## is a request to the SOLVER, not a change to the equations.  The kind
    ## is carried rather than assumed: `pnjlim` takes (IS, VT), `fetlim`
    ## takes (vto) and `limvds` takes nothing, and the generated `limit()`
    ## dispatches on it.
    limits = []

    def _strip_limits(expr):
        """Record every `$limit` marker in `expr` and return it unmarked."""
        subs_l = {}
        for app in sorted(expr.atoms(_Limit), key=sympy.default_sort_key):
            probe, pars, grp = _limit_parts(app)
            b_l = probe.branch_or_node
            key = (index_of[('node', b_l.plus.name)],
                   index_of[('node', b_l.minus.name)])
            prev = [k for k in limits if k[0] == key]
            if prev and prev[0][1] != app.kind:
                raise ValueError(
                    'branch %s is limited twice, as %r and as %r; one probe '
                    'carries one limiter'
                    % (_fmt_branch(branch_key(b_l)), prev[0][1], app.kind))
            if not prev:
                limits.append((key, app.kind, pars, grp))
            subs_l[app] = probe
        return expr.subs(subs_l) if subs_l else expr

    ## Intermediates FIRST, and this is not tidiness: a `var()` definition
    ## is never visited by the statement loop, so before this a
    ## `var(limit_pnj(...))` carried its marker all the way into
    ## `_chain_compile` and died there with "Unsupported ... _Limit".  Any
    ## real MOSFET is chained -- that is what `var()` is for -- so $limit
    ## was unusable by exactly the models 10.3(a) exists to serve.
    for _k, (_sym, _expr) in enumerate(intermediates):
        _new = _strip_limits(_expr)
        if _new is not _expr:
            intermediates[_k] = (_sym, _new)
    for st in statements:
        st.rhs = _strip_limits(st.rhs)

    ## ------------------------------------------------------------------
    ## Intermediates (`var`): resolve each definition in order, and work
    ## out which of them depend on the solution.  That second part is not
    ## cosmetic -- `_split_terms` decides "is this a source?" by asking
    ## whether a term touches an x-symbol, and an intermediate is an
    ## opaque symbol, so without this every term mentioning one would be
    ## routed into `u(t)` and the element would have no conductance at
    ## all.
    chain_defs = []
    xdep_syms, tdep_syms = set(), set()
    for sym, expr in intermediates:
        rexpr = resolve(expr)
        chain_defs.append((sym, rexpr))
        fs = rexpr.free_symbols
        if (fs & set(xsyms)) or (fs & xdep_syms):
            xdep_syms.add(sym)
        if rexpr.has(TIME) or (fs & tdep_syms):
            tdep_syms.add(sym)
    xset_split = set(xsyms) | xdep_syms

    ## ------------------------------------------------------------------
    ## Pass 2 -- assemble the residual vectors.
    ivec = [sympy.Integer(0)] * n
    qvec = [sympy.Integer(0)] * n
    uvec = [sympy.Integer(0)] * n
    acvec = [sympy.Integer(0)] * n
    CYacc = {}
    ## name -> [shared power, {row: summed amplitude}].  See `_NamedNoise`:
    ## every member of a group is the SAME fluctuation, so the group has one
    ## power and one injection vector, not one of each per contribution.
    named_noise = {}

    for st in statements:
        b = st.lhs.branch_or_node
        p = index_of[('node', b.plus.name)]
        m_ = index_of[('node', b.minus.name)]
        if st.lhs.quantity == 'I':
            ipart, qpart, upart, npart, acpart = _split_terms(
                resolve(st.rhs), xset_split, tdep_syms)
            ivec[p] += ipart; ivec[m_] -= ipart
            qvec[p] += qpart; qvec[m_] -= qpart
            uvec[p] += upart; uvec[m_] -= upart
            acvec[p] += acpart; acvec[m_] -= acpart
            for name, amp, pwr in npart:
                if name is None:
                    ## An I-contributed noise PSD stamps like a noisy
                    ## conductance: the R.CY pattern (elements.py).
                    psd = amp ** 2 * pwr
                    CYacc[(p, p)] = CYacc.get((p, p), 0) + psd
                    CYacc[(m_, m_)] = CYacc.get((m_, m_), 0) + psd
                    CYacc[(p, m_)] = CYacc.get((p, m_), 0) - psd
                    CYacc[(m_, p)] = CYacc.get((m_, p), 0) - psd
                    continue
                group = named_noise.setdefault(name, [pwr, {}])
                if group[0] != pwr:
                    raise NotImplementedError(
                        'noise sources named %r contribute different powers '
                        '(%r and %r).  Members of a correlation group are ONE '
                        'fluctuation, so they share one power; put what '
                        'differs between branches in the scale factor, which '
                        'multiplies the amplitude.' % (name, group[0], pwr))
                group[1][p] = group[1].get(p, 0) + amp
                group[1][m_] = group[1].get(m_, 0) - amp
        else:
            ## V-contribution: KCL rows carry the branch current; the
            ## branch row is -(v_p - v_m) + rhs = 0 (the elements.py
            ## convention -- compare Idtmod's output row).
            bi = index_of[('branch', branch_key(b))]
            ivec[p] += xsyms[bi]
            ivec[m_] -= xsyms[bi]
            ipart, qpart, upart, npart, acpart = _split_terms(
                resolve(st.rhs), xset_split, tdep_syms)
            if acpart != 0:
                raise NotImplementedError(
                    'ac_stim on a V-contributed branch is not supported '
                    'yet; contribute it with I(...)')
            if npart:
                raise NotImplementedError(
                    'noise contributions on a V-contributed branch are not '
                    'supported yet; contribute noise with I(...)')
            vp, vm = xsyms[p], xsyms[m_]
            ivec[bi] += -(vp - vm) + ipart
            qvec[bi] += qpart
            uvec[bi] += upart

    ## ------------------------------------------------------------------
    ## Correlated groups: one fluctuation reaching several branches.
    ##
    ## A group injects the current vector ``w * n``, with ``n`` a unit
    ## process of power ``pwr`` and ``w`` the summed stamp of its members,
    ## so its correlation matrix is the RANK-ONE ``pwr * w w^dagger``.  The
    ## lone-source stamp above is the same expression with a single member:
    ## ``w = amp * (e_p - e_m)`` gives back the 2x2 conductance pattern.
    ##
    ## Off-diagonal blocks are what a group is FOR.  They are the only way
    ## to say that two branches carry the same noise -- summing two
    ## independent sources of the same power says something different, and
    ## says it wrong wherever the two paths interfere.
    for name, (pwr, wvec) in named_noise.items():
        rows = sorted(wvec)
        for r_ in rows:
            for c_ in rows:
                CYacc[(r_, c_)] = (CYacc.get((r_, c_), 0)
                                   + pwr * wvec[r_] * _conj(wvec[c_]))

    ## State equations: ds/dt = arg  ->  q-row s, i-row -arg.
    dc_pins = {}     # x-index -> (pin expression for i, seed expression)
    periodic = []    # (x-index, modulus expr, offset expr)  [idtmod only]
    for sym, kind, args in states:
        k = index_of[('state', sym.name)]
        arg_i, arg_q, arg_u, arg_n, arg_ac = _split_terms(
            resolve(args[0]), xset_split, tdep_syms)
        if arg_ac != 0:
            raise NotImplementedError('ac_stim inside idt/idtmod is not '
                                      'supported')
        if arg_n:
            raise NotImplementedError('noise inside idt/idtmod is not '
                                      'supported')
        if arg_q != 0:
            raise NotImplementedError('ddt inside idt/idtmod is not supported')
        qvec[k] += subst[sym]
        ivec[k] += -arg_i
        uvec[k] += -arg_u
        if kind == 'lap':
            continue
        if len(args) > 1:
            ic = resolve(args[1])
            if not _quantity_free(ic):
                raise NotImplementedError(
                    'an idt/idtmod ic must not depend on circuit quantities')
            if kind == 'idtmod':
                m = resolve(args[2]) if len(args) > 2 else sympy.Integer(1)
                o = resolve(args[3]) if len(args) > 3 else sympy.Integer(0)
                ic = ic - m * _wrapfloor((ic - o) / m)
            dc_pins[k] = ic
        if kind == 'idtmod':
            m = resolve(args[2]) if len(args) > 2 else sympy.Integer(1)
            o = resolve(args[3]) if len(args) > 3 else sympy.Integer(0)
            if not (_quantity_free(m) and _quantity_free(o)):
                raise NotImplementedError(
                    'idtmod modulus/offset must not depend on circuit '
                    'quantities')
            periodic.append((k, m, o))

    ## The DC variant: pinned state rows become `s - ic` (fold-into-i, the
    ## elements.py convention), everything else identical.
    ivec_dc = list(ivec)
    uvec_dc = list(uvec)
    for k, ic in dc_pins.items():
        ivec_dc[k] = subst[[s for s, _k2, _a in states
                            if index_of[('state', s.name)] == k][0]] - ic
        uvec_dc[k] = sympy.Integer(0)

    ## `$param_given` symbols actually referenced, in a stable order --
    ## they become trailing arguments of every compiled function and are
    ## bound per instance from `ipar`.
    ## The CHAIN is scanned too, and for the same reason the state
    ## operators and `$limit` had to be: on the let-chain path the
    ## assembled vectors are mostly `var()` symbols, so a
    ## `$param_given` used inside a definition appeared in none of them,
    ## became no argument, and the generated function raised
    ## `NameError: _hdl_given_x` at call time -- from a model that
    ## compiled without a word.
    given_syms = sorted(
        {a for vec in (ivec, qvec, uvec, acvec, ivec_dc,
                       [e_ for _s2, e_ in chain_defs])
         for e in vec for a in e.free_symbols
         if str(a).startswith('_hdl_given_')},
        key=lambda a: str(a))
    given_names = [str(a)[len('_hdl_given_'):] for a in given_syms]

    ## ------------------------------------------------------------------
    ## Jacobians and compilation.
    ## `Matrix.jacobian` differentiates the assembled expression TREE.
    ## For a model that reuses intermediates -- every surface-potential
    ## model does -- that tree is exponential in nesting depth even though
    ## the DAG is linear, so this line alone is the wall.  Models that
    ## declare intermediates with `var()` take the let-chain path instead
    ## and never build these.
    if chain_defs:
        G = C = G_dc = sympy.zeros(n)
    else:
        G = sympy.Matrix([ivec]).jacobian(xsyms)
        C = sympy.Matrix([qvec]).jacobian(xsyms)
        G_dc = sympy.Matrix([ivec_dc]).jacobian(xsyms) if dc_pins else G
    CY = sympy.zeros(n)
    for (r_, c_), psd in CYacc.items():
        CY[r_, c_] = psd
    dudt = [sympy.diff(u_, TIME) for u_ in uvec]

    x = sympy.DeferredVector('x')
    xsubs = dict(zip(xsyms, [x[i] for i in range(n)]))

    NUMPY_MODULES = [dict(_KERNEL_NUMPY, _wrapfloor=np.floor), 'numpy']

    def compile_x(expr_or_list, extra=()):
        e = expr_or_list.subs(xsubs) if hasattr(expr_or_list, 'subs') else \
            [e_.subs(xsubs) for e_ in expr_or_list]
        return sympy.lambdify(
            [x] + paramsyms + [TEMP] + given_syms + list(extra), e,
            modules=NUMPY_MODULES, cse=True)

    def compile_t(exprs):
        return sympy.lambdify([TIME] + paramsyms + [TEMP] + given_syms,
                              exprs, modules=NUMPY_MODULES, cse=True)

    if chain_defs:
        ## THE LET-CHAIN PATH.  Four chain compilations replace the eager
        ## stamps; each walks the DAG once and accumulates derivatives
        ## forward, so cost is linear in the number of intermediates
        ## rather than exponential in their nesting.
        chain_args = [x] + paramsyms + [TEMP] + given_syms
        _mods = dict(_KERNEL_NUMPY, _wrapfloor=np.floor)
        ## The compiled signature takes the solution VECTOR, while the
        ## chain is written in terms of the scalar unknowns, so the body
        ## opens by unpacking one into the other.
        _unpack = [(xs.name, 'x[%d]' % k) for k, xs in enumerate(xsyms)]
        ## The C text is printed for every chain that takes `x`, whatever
        ## backend is selected: printing adds ~19% to a cold compile
        ## (see `EMIT_C_SOURCE`), and carrying it means a class served
        ## from the compile cache can switch backend with no symbolic
        ## work at all -- the alternative, keying the cache on the
        ## backend, would make `PYCIRCUIT_HDL_BACKEND=c` a 41 s cold
        ## compile of the MOSFET again.
        cc = lambda out, jac: _chain_compile(
            chain_defs, out, chain_args, want_jacobian_of=jac,
            xsyms=xsyms, modules_map=_mods, unpack=_unpack,
            emit_c=EMIT_C_SOURCE)

        ## The source vectors are x-free BY CONSTRUCTION -- that is what
        ## made them sources -- so they compile against the x-free part
        ## of the chain and take no `x` argument.  They still have to go
        ## through the chain: an intermediate is an opaque symbol, and
        ## lambdify handed one it cannot resolve produces a function that
        ## quietly evaluates to nothing useful rather than complaining.
        free_defs = [(sym, e_) for sym, e_ in chain_defs
                     if sym not in xdep_syms]
        cu = lambda out, extra=(): _chain_compile(
            free_defs, out, list(extra) + paramsyms + [TEMP] + given_syms,
            modules_map=_mods)
        funcs = dict(i=cc(ivec, None), q=cc(qvec, None),
                     G=cc(ivec, True), C=cc(qvec, True),
                     CY=_chain_compile(
                         chain_defs,
                         [[CY[r_, c_] for c_ in range(n)]
                          for r_ in range(n)],
                         chain_args + [FREQ],
                         modules_map=_mods, unpack=_unpack),
                     u=cu(uvec, (TIME,)),
                     dudt=cu(dudt, (TIME,)), uac=cu(acvec))
        ## The DC variants differ from the transient ones ONLY where a
        ## state is pinned (`ivec_dc` is built from `ivec` by replacing
        ## the pinned rows).  With nothing pinned they are the same
        ## vectors, and compiling them again produced the same source
        ## text at the same cost -- for `G_dc` a second forward-mode
        ## Jacobian of the whole chain, measured at about a fifth of the
        ## library's import time.  Share the function instead.
        if dc_pins:
            funcs['i_dc'] = cc(ivec_dc, None)
            funcs['G_dc'] = cc(ivec_dc, True)
            funcs['u_dc'] = cu(uvec_dc, (TIME,))
        else:
            funcs['i_dc'], funcs['G_dc'], funcs['u_dc'] = \
                funcs['i'], funcs['G'], funcs['u']
    else:
      funcs = dict(
        i=compile_x(ivec),
        q=compile_x(qvec),
        G=compile_x(G),
        C=compile_x(C), CY=compile_x(CY, extra=(FREQ,)),
        u=compile_t(uvec), dudt=compile_t(dudt),
        uac=sympy.lambdify(paramsyms + [TEMP] + given_syms, acvec,
                           modules=NUMPY_MODULES, cse=True),
      )
      if dc_pins:
          funcs['i_dc'] = compile_x(ivec_dc)
          funcs['G_dc'] = compile_x(G_dc)
          funcs['u_dc'] = compile_t(uvec_dc)
      else:
          funcs['i_dc'], funcs['G_dc'], funcs['u_dc'] = \
              funcs['i'], funcs['G'], funcs['u']

    ## Pure single-device forms for the JAX batched path (eval_i_pure /
    ## eval_q_pure) reuse the SAME symbolic vectors, compiled lazily with
    ## sympy's jax printer on first use (jax may not be installed).
    pure_spec = None if chain_defs else \
        dict(ivec=ivec, qvec=qvec, xsyms=xsyms, n=n)

    ## The SYMBOLIC forms, kept alongside the compiled ones.  A symbolic
    ## toolkit asked for `i(x)` with symbols in `x` was being served by a
    ## numpy lambda, which works by duck typing for plain arithmetic and
    ## fails the moment the expression contains a `floor`, a `Piecewise`
    ## or anything else numpy evaluates eagerly.  Substituting into these
    ## is exact and cannot fail that way.
    sym_spec = None if chain_defs else \
        dict(i=ivec, q=qvec, G=list(G), C=list(C),
             xsyms=xsyms, paramsyms=paramsyms,
             given_syms=given_syms, shape=(n, n))

    ## ------------------------------------------------------------------
    ## PCNR participation (Aadithya/Keiter/Mei).  A device joins PCNR by
    ## making its limited quantity an explicit unknown and evaluating at
    ## ITS OWN copy, so no two devices can fight over a shared
    ## linearisation point.  The layer needs three things from the device
    ## (`pcnr.augmented_system`): which terminal pair the quantity spans,
    ## the terminal currents as a function of that quantity ALONE, and
    ## their derivative -- because it REMOVES the device's whole i/G
    ## contribution and re-stamps it at the limited voltage.
    ##
    ## So an element qualifies only if its entire current is a function of
    ## one branch voltage, it stores no charge (the layer refuses charge
    ## on a participant in transient: the charge term would have to move
    ## to the limited unknown too), and it introduces no states.  The
    ## limiter itself is generated when the current is recognisably
    ## exponential -- the case limiting exists for -- and the exponential
    ## scale is read off the expression rather than assumed.
    ## PCNR participation, detected PER CONTRIBUTION rather than by
    ## reverse-engineering the assembled `ivec`.  A device joins by making
    ## each limited quantity an explicit unknown and evaluating at ITS OWN
    ## copy, so the layer needs, per junction: the terminal pair it spans,
    ## the terminal currents as a function of that quantity ALONE, and
    ## their derivative -- it REMOVES the device's whole i/G contribution
    ## and re-stamps it at `v_lim`.
    ##
    ## So the rule is: every I-contribution must be a function of its own
    ## branch voltage alone, and at least one must be exponential (the
    ## case limiting exists for).  A device may own SEVERAL such branches
    ## -- a two-junction device is the point of this pass -- and each gets
    ## its own limited unknown, which is exactly what removes the clash
    ## PCNR is about.  Anything else (a V-contribution, a generated state,
    ## a current spanning two branch voltages) disqualifies the element:
    ## better to declare nothing than to claim a capability falsely.
    ## THE CHAINED PATH IS ADMITTED TOO, and the detector never flattens
    ## the chain to do it.  Reading `f.atoms(sympy.exp)` off an assembled
    ## expression is what a `var()` defeats -- behind an intermediate the
    ## exponential is an opaque symbol, so the detector saw none and the
    ## gate refused rather than declare the device un-limitable silently.
    ## Restoring the old reading by INLINING the chain is the one thing
    ## `var()` exists to prevent (measured on PSP's `sp_s`: the flattened
    ## form did not finish in ten minutes).  It is also unnecessary --
    ## only two facts about the flattened expression are wanted, and both
    ## are available by walking the chain:
    ##
    ##   (a) which branch voltages the contribution reaches.  Prune the
    ##       chain to what this contribution reads (`_chain_prune`), put
    ##       V(kp) = v + V(km) into each definition SEPARATELY, and
    ##       require that no x-symbol survives in any of them.  That is
    ##       the same statement as the flat `f.free_symbols & xsyms`
    ##       test, made one definition at a time.
    ##   (b) the exponential scale.  Forward-accumulate d/dv over the
    ##       pruned chain (`_chain_forward_d`), which yields d(arg)/dv
    ##       for every `exp(arg)` in it as a SMALL expression over the
    ##       derivative chain -- the machinery `_chain_compile`'s
    ##       Jacobian already uses.
    ##
    ## Both are linear in the number of definitions.  With no chain the
    ## pruned chain is empty, `_chain_d` degenerates to `sympy.diff` and
    ## this is the flat detector unchanged.
    ##
    ## The refusals that are still correct are kept: a generated state, a
    ## branch-current unknown, a V-contribution, a current spanning two
    ## branch voltages, no single exponential scale -- and, new here, a
    ## contribution reading anything `pcnr_i(v, params, T)` is not given
    ## (a `$param_given` flag, TIME).  Better to decline than to claim a
    ## capability falsely.
    ## TWO ROUTES, tried in order (roadmap sec. 15, Stage 1):
    ##
    ##   1. the DECLARED-PROBE route -- the model carries `$limit` probes,
    ##      and they ARE its limited unknowns: one block per device,
    ##      `pcnr_i(v) -> R^t`, `pcnr_didv(v) -> R^{t x m}`, the declared
    ##      limiter as the law.  This is what admits a transistor, whose
    ##      transport current reads two junction voltages at once;
    ##   2. the LEGACY scalar route below, unchanged, when there is no
    ##      `limit_spec`: one junction per exponential branch, `IS`/`VT`
    ##      read off the expression.  Every diode in the tree takes this
    ##      route and produces the same numbers it always did.
    pcnr_spec = None
    pcnr_refusal = None
    pcnr_vector = None
    if limits:
        pcnr_vector, pcnr_refusal = _pcnr_declared_route(
            limits, ivec, chain_defs, xsyms, xlabels, paramsyms, given_syms,
            states, vbranches)
    elif not states and not vbranches and len(terminalnodes) >= 2:
        cands, ok = [], True
        _xset = set(xsyms)
        _allowed = set(paramsyms) | {TEMP}
        for st in statements:
            if st.lhs.quantity != 'I':
                ok = False
                pcnr_refusal = ('a V-contribution (%s) -- PCNR limits '
                                'currents' % st.lhs)
                break
            b = st.lhs.branch_or_node
            kp = index_of[('node', b.plus.name)]
            km = index_of[('node', b.minus.name)]
            ipart, _q, _u, _nz, _ac = _split_terms(resolve(st.rhs),
                                                   xset_split, tdep_syms)
            if ipart == 0:
                continue                     # pure charge/noise: harmless
            vsym = sympy.Symbol('_pcnr_v', real=True)
            _vsub = {xsyms[kp]: vsym + xsyms[km]}

            def _sub_x(e_, _s=_vsub, _x=_xset):
                ## `expand` only where there is something left to cancel.
                ## On a chain each definition is small, but expanding one
                ## that has nothing to gain is still the operation
                ## `_split_terms` refuses to perform on a whole model.
                e2 = e_.subs(_s)
                return sympy.expand(e2) if (e2.free_symbols & _x) else e2

            defs_j = [(sym, _sub_x(e_))
                      for sym, e_ in _chain_prune(chain_defs, [ipart])]
            f = _sub_x(ipart)
            _fs = set(f.free_symbols)
            for _s2, _e2 in defs_j:
                _fs |= _e2.free_symbols
            if _fs & _xset:
                ok = False                   # not a function of V(b) alone
                pcnr_refusal = ('%s reads other node voltages (%s) -- its '
                                'current is not a function of its own branch '
                                'voltage alone' % (st.lhs, ', '.join(
                                    sorted(xlabels[xsyms.index(q)]
                                           for q in _fs & _xset))))
                break
            if _fs - _allowed - {vsym} - {sym for sym, _e2 in defs_j}:
                ok = False        # reads what pcnr_i is not handed
                pcnr_refusal = ('%s reads %s, which pcnr_i is not handed'
                                % (st.lhs, ', '.join(sorted(str(q) for q in
                                   _fs - _allowed - {vsym}
                                   - {sym for sym, _e2 in defs_j}))))
                break

            ## (b) -- the derivative chain, built once and reused for
            ## every exponential in the contribution AND for `didv`.
            ddefs, dmap = _chain_forward_d(defs_j, vsym)
            _defset = {sym for sym, _e2 in defs_j}
            _comb = defs_j + ddefs
            dfdv = _chain_d(f, vsym, _defset, dmap)

            def _vfree(e_, _c=_comb, _v=vsym):
                """Is `e_` constant in v, reading through the chain?"""
                if _v in e_.free_symbols:
                    return False
                return not any(_v in ee.free_symbols
                               for _s2, ee in _chain_prune(_c, [e_]))

            scales = set()
            for _e2 in [e2 for _s2, e2 in defs_j] + [f]:
                for ex in sorted(_e2.atoms(sympy.exp),
                                 key=sympy.default_sort_key):
                    a = _chain_d(ex.args[0], vsym, _defset, dmap)
                    if a != 0 and _vfree(a):
                        scales.add(sympy.simplify(1 / a))
            if len(scales) != 1:
                ok = False                   # no single exponential scale
                pcnr_refusal = ('%s has %s exponential scale%s -- a linear '
                                'contribution refuses the whole device, and '
                                'so does one with two different scales'
                                % (st.lhs, len(scales) or 'no',
                                   '' if len(scales) == 1 else 's'))
                break
            VT_eff = scales.pop()
            if defs_j:
                ## `IS = VT * di/dv|_(v=0)`, but the derivative lives in
                ## the chain, so v=0 is substituted into the CHAIN rather
                ## than into a flattened expression; the compile below
                ## uses the zeroed copy.
                cands.append(dict(terminals=(kp, km), vsym=vsym,
                                  VT=VT_eff, IS=VT_eff * dfdv,
                                  chain=dict(defs=defs_j, ddefs=ddefs,
                                             out=f, dout=dfdv)))
            else:
                cands.append(dict(terminals=(kp, km), vsym=vsym, f=f,
                                  dfdv=dfdv, VT=VT_eff, chain=None,
                                  IS=sympy.simplify(
                                      VT_eff * dfdv.subs(vsym, 0))))
        if ok and cands:
            pcnr_spec = cands
        elif ok and not cands:
            pcnr_refusal = 'no resistive current contribution at all'
    else:
        pcnr_refusal = ('%s -- the device carries %s'
                        % ('a generated state' if states else
                           'a branch-current unknown' if vbranches else
                           'fewer than two terminals',
                           '%d state(s)' % len(states) if states else
                           '%d V-contributed branch(es)' % len(vbranches)
                           if vbranches else 'one terminal'))

    pcnr_funcs = None
    if pcnr_spec is not None:
        pcnr_funcs = []

        def _chain_one(defs_, out_, args_):
            """One scalar-valued function compiled through the chain."""
            fn_ = _chain_compile(
                defs_, [out_], args_,
                modules_map=dict(_KERNEL_NUMPY, _wrapfloor=np.floor))
            return _first_of(fn_)

        for spec_j in pcnr_spec:
            vs_ = spec_j['vsym']
            ch_ = spec_j['chain']
            if ch_ is None:
                ## The traced backend calls these under `jit`, so they
                ## need a jax-printed twin; built lazily so importing the
                ## module does not require jax.
                _sym_j = dict(f=spec_j['f'], dfdv=spec_j['dfdv'], vsym=vs_)
                pcnr_funcs.append(dict(
                    sym=_sym_j, chain=None,
                    terminals=spec_j['terminals'],
                    i=sympy.lambdify([vs_] + paramsyms + [TEMP], spec_j['f'],
                                     modules=NUMPY_MODULES, cse=True),
                    didv=sympy.lambdify([vs_] + paramsyms + [TEMP],
                                        spec_j['dfdv'], modules=NUMPY_MODULES,
                                        cse=True),
                    VT=sympy.lambdify(paramsyms + [TEMP], spec_j['VT'],
                                      modules=NUMPY_MODULES),
                    IS=sympy.lambdify(paramsyms + [TEMP], spec_j['IS'],
                                      modules=NUMPY_MODULES),
                ))
                continue
            ## The chained twin of the four, compiled through
            ## `_chain_compile` in place of `sympy.lambdify`: the chain
            ## is never flattened, here or in the detector above.  `VT`
            ## is v-free by construction (that is what admitted the
            ## contribution), so pruning drops every v-dependent
            ## definition from its body and it needs no `v` argument;
            ## `IS` reads the SAME chain with v = 0 substituted into each
            ## definition.
            _args_v = [vs_] + paramsyms + [TEMP]
            _comb = ch_['defs'] + ch_['ddefs']
            _zeroed = [(sym, e_.subs(vs_, 0)) for sym, e_ in _comb]
            _ch_j = dict(defs=ch_['defs'], ddefs=ch_['ddefs'],
                         out=ch_['out'], dout=ch_['dout'], vsym=vs_)
            pcnr_funcs.append(dict(
                sym=None, chain=_ch_j,
                terminals=spec_j['terminals'],
                i=_chain_one(ch_['defs'], ch_['out'], _args_v),
                didv=_chain_one(_comb, ch_['dout'], _args_v),
                VT=_chain_one(_comb, spec_j['VT'], paramsyms + [TEMP]),
                IS=_chain_one(_zeroed, spec_j['IS'], paramsyms + [TEMP]),
            ))

    state_meta = dict(
        statenames=statenodenames,
        dc_pins={k: sympy.lambdify(paramsyms, ic, modules=NUMPY_MODULES)
                 for k, ic in dc_pins.items()},
        periodic=[(k, sympy.lambdify(paramsyms, m, modules=NUMPY_MODULES),
                   sympy.lambdify(paramsyms, o, modules=NUMPY_MODULES))
                  for k, m, o in periodic],
        has_time_dependence=any(u_.has(TIME) for u_ in uvec),
    )

    ## An x-independent G/C is a CONSTANT stamp: a hand-written linear
    ## element returns its stored matrix by reference, and recomputing one
    ## per Newton iteration is the whole measured gap (benchmarks/
    ## hdl_overhead.py: 28x on a resistor's G before this).  Recorded here,
    ## cached per parameter set in the metaclass.
    xset_ = set(xsyms)
    const_G = (not chain_defs
               and not any((e.free_symbols & xset_) for e in G))
    const_C = (not chain_defs
               and not any((e.free_symbols & xset_) for e in C))

    cross_spec = None
    if crossings:
        if chain_defs:
            cross_f = _chain_compile(
                chain_defs, [resolve(st.expr) for st in crossings],
                [x] + paramsyms + [TEMP] + given_syms,
                modules_map=dict(_KERNEL_NUMPY, _wrapfloor=np.floor),
                unpack=[(xs.name, 'x[%d]' % k)
                        for k, xs in enumerate(xsyms)])
        else:
            cross_exprs = [resolve(st.expr).subs(xsubs) for st in crossings]
            cross_f = sympy.lambdify(
                [x] + paramsyms + [TEMP] + given_syms, cross_exprs,
                modules=NUMPY_MODULES, cse=True)
        cross_spec = dict(f=cross_f,
                          directions=[st.direction for st in crossings])

    ## Which terminal each probe's limited value is written back to.
    ## Default: the PLUS terminal, which is what the single-probe case has
    ## always done.  If an earlier probe already moved that terminal, move
    ## the MINUS one instead -- otherwise the second write undoes the first,
    ## which is the same hazard `limit_junctions` handles with its `move`
    ## field for a BJT's two junctions sharing a base.  With both terminals
    ## already pinned there is no choice left that keeps every probe's value,
    ## and saying so beats writing a silently wrong one.
    ## GROUPED probes (`limit_together`) are exempt from all of it: their
    ## write-back is solved over the whole device at runtime, where a
    ## constraint that cannot be honoured is dropped by size rather than by
    ## declaration position, so there is nothing for a compile-time rule to
    ## decide and no configuration it has to refuse.  `move` is still filled
    ## in for them -- the plus terminal, the historical default -- so that
    ## every `limit_spec` entry has the same shape; the grouped path never
    ## reads it.
    limit_spec, _moved = [], set()
    _groups = collections.OrderedDict()
    for k, kind_l, pars_l, grp_l in limits:
        ra, rb = int(k[0]), int(k[1])
        if grp_l is not None:
            move = ra
            _groups.setdefault(grp_l[0], (grp_l[2], []))[1].append(
                (grp_l[1], len(limit_spec)))
        elif ra not in _moved:
            move = ra
            _moved.add(move)
        elif rb not in _moved:
            move = rb
            _moved.add(move)
        else:
            raise ValueError(
                'the $limit probes over-determine this device: both '
                'terminals of branch (%s,%s) have already been moved by '
                'earlier probes, so this one cannot be written back without '
                'undoing them.  Declare them with limit_together(), whose '
                'write-back is solved over the whole device at once.'
                % (xlabels[ra], xlabels[rb]))
        ## `given_syms` BELONGS IN THIS SIGNATURE.  Every other compiled
        ## function in this file is called with `_args_of`, which is
        ## `params + [T] + givenness flags`; this one alone was lambdified
        ## over `paramsyms + [TEMP]` and then called with the wider list.
        ## A model declaring BOTH a `$limit` probe and `$param_given`
        ## therefore compiled clean and raised `TypeError: takes N
        ## positional arguments but N+1 were given` from inside `limit()`
        ## on the first Newton iteration.  Nothing had ever used the two
        ## together until the Gummel-Poon BJT wanted `param_given('rbm')`
        ## for SPICE's `RBM defaults to RB`, which a default VALUE cannot
        ## express.  Widened rather than sliced at the call site, so a
        ## limiter's parameters may legitimately read a givenness flag.
        limit_spec.append(((ra, rb), kind_l, move,
                           tuple(_limit_par_fn(resolve(e), kind_l, chain_defs,
                                               xsyms, xlabels,
                                               paramsyms + [TEMP]
                                               + given_syms, NUMPY_MODULES,
                                               x=x)
                                 for e in pars_l)))

    ## `limit_groups`: `(sequential, [spec index, in DECLARATION order])`.
    ## Kept beside `limit_spec` rather than inside it because every reader
    ## of a spec entry -- `explain`, `check_jacobians`, four test files --
    ## unpacks it as a 4-tuple, and widening it would have been a change to
    ## an interface for the sake of a field only `limit()` reads.
    limit_groups = [(seq, [i for _slot, i in sorted(items)])
                    for _gid, (seq, items) in _groups.items()]
    for _seq, _idx in limit_groups:
        if len(_idx) < 2:                                # pragma: no cover
            raise ValueError('a $limit group needs at least two probes')

    ## The declared-probe route's compiled twin, built HERE because its
    ## limiter reads `limit_spec`'s parameter callables (compiled just
    ## above; compiling them a second time would be the same work twice).
    ## `_chain_compile`'s forward-accumulation Jacobian over the probe
    ## symbols gives the `t x m` block directly; a flat model takes
    ## `lambdify` + `Matrix.jacobian`, as the legacy route does.
    pcnr_vec = None
    if pcnr_vector is not None:
        _vs = pcnr_vector['vsyms']
        _args_v = _vs + paramsyms + [TEMP] + given_syms
        _mods = dict(_KERNEL_NUMPY, _wrapfloor=np.floor)
        if chain_defs:
            _vec_i = _chain_compile(pcnr_vector['defs'], pcnr_vector['i_exprs'],
                                    _args_v, modules_map=_mods)
            _vec_d = _chain_compile(pcnr_vector['defs'], pcnr_vector['i_exprs'],
                                    _args_v, want_jacobian_of=True, xsyms=_vs,
                                    modules_map=_mods)
            _vec_sym = None
        else:
            _jac = sympy.Matrix(pcnr_vector['i_exprs']).jacobian(_vs).tolist()
            _vec_i = sympy.lambdify(_args_v, pcnr_vector['i_exprs'],
                                    modules=NUMPY_MODULES, cse=True)
            _vec_d = sympy.lambdify(_args_v, _jac, modules=NUMPY_MODULES,
                                    cse=True)
            _vec_sym = dict(i=pcnr_vector['i_exprs'], didv=_jac)
        pcnr_vec = dict(probes=pcnr_vector['probes'], kinds=pcnr_vector['kinds'],
                        spec_idx=pcnr_vector['spec_idx'],
                        redundant=pcnr_vector['redundant'],
                        m=len(_vs), t=pcnr_vector['t'], vsyms=_vs,
                        i=_vec_i, didv=_vec_d, sym=_vec_sym,
                        chain=(dict(defs=pcnr_vector['defs'],
                                    out=pcnr_vector['i_exprs'])
                               if chain_defs else None),
                        limit_pars=[tuple(pf) for _k, _kind, _m, pf
                                    in limit_spec])

    branchpairs = [branch_key(br) for br in vbranches]
    internalnames = [nd.name for nd in internalnodes]

    return dict(terminalnames=terminalnames, paramnames=paramnames,
                funcs=funcs, pure_spec=pure_spec, state_meta=state_meta,
                branchpairs=branchpairs, internalnames=internalnames,
                const_G=const_G, const_C=const_C, pcnr_funcs=pcnr_funcs,
                pcnr_refusal=pcnr_refusal, pcnr_vector=pcnr_vec,
                given_names=given_names, limit_spec=limit_spec,
                limit_groups=limit_groups,
                cross_spec=cross_spec, sym_spec=sym_spec,
                has_ac=any(e != 0 for e in acvec),
                chained=bool(chain_defs),
                collapse_mask=collapse_mask,
                collapse_conds=[
                    sympy.lambdify(paramsyms, st.when,
                                   modules=NUMPY_MODULES)
                    for st in collapses])


def _fmt_branch(key):
    """Render a `(plus, minus, name)` branch key for a message."""
    plus, minus, name = key
    return '(%s,%s)' % (plus, minus) if name is None else \
        '(%s,%s) named %r' % (plus, minus, name)


def compile_chain(builder, args, wrt=None, modules_map=None):
    """Compile a `var`-using expression builder OUTSIDE a model.

    `builder` is a zero-argument callable that returns a sympy
    expression, using `var()` freely for its intermediates.  `args` is
    the resulting function's parameter list.  With `wrt` given (a list of
    symbols), the compiled function returns the gradient with respect to
    them instead of the value.

    This exists because a kernel routine is worth testing on its own,
    and the fallback of "just build one expression when no model is
    compiling" is not available: the reuse that makes `var()` necessary
    inside a model makes the standalone expression equally unbuildable.
    (Measured on PSP's `sp_s`: the flattened form did not finish in ten
    minutes; through the chain it compiles in under a second.)
    """
    _VAR_STACK.append([])
    try:
        expr = builder()
    finally:
        defs = _VAR_STACK.pop()
    if wrt is None:
        return _chain_compile(defs, [expr], args, modules_map=modules_map)
    fn = _chain_compile(defs, [expr], args, want_jacobian_of=True,
                        xsyms=list(wrt), modules_map=modules_map)
    return fn


class _ChainPrinter(object):
    """`NumPyPrinter`, but `Piecewise` becomes nested `numpy.where`.

    sympy prints a `Piecewise` as `numpy.select([conds], [values])`, and
    `select` is built for arrays: for the SCALAR arguments a compiled
    device model is called with, it spends its time broadcasting and
    allocating.  Profiled on the surface-potential MOSFET's Jacobian --
    950 `select` calls per evaluation, 60% of the total runtime.

    `numpy.where` has exactly the same both-arms-evaluated semantics,
    which the kernel's safety work depends on, and none of the
    machinery.  Nested right-to-left, it reproduces `Piecewise`'s
    first-match-wins ordering.
    """

    def __init__(self):
        self._p = self._printer_class()()

    @classmethod
    def _printer_class(cls):
        """The `NumPyPrinter` subclass this printer drives.

        A classmethod rather than a body of `__init__` so that the C
        printer below can derive from the SAME class and inherit every
        decision made here -- `Piecewise` nesting, `Min`/`Max` as binary
        calls, the DSL primitives as calls -- overriding only the leaf
        syntax.  Two printers that agree on structure and differ only in
        spelling is what makes bitwise agreement between the backends a
        property of the emission rather than a hope.
        """
        from sympy.printing.numpy import NumPyPrinter
        outer = cls

        class _P(NumPyPrinter):
            def _print_Piecewise(self, expr):
                return outer._piecewise(self, expr)

            def _print_Min(self, expr):
                return outer._minmax(self, expr, 'minimum')

            def _print_Max(self, expr):
                return outer._minmax(self, expr, 'maximum')

            ## `NumPyPrinter` is strict about functions it does not know,
            ## so the two `safe_div` primitives print themselves.  Kept
            ## as CALLS rather than inlined: the argument may be a large
            ## sub-expression and inlining would emit it three times.
            def _print__recip2(self, expr):
                return '_recip2(%s, %s)' % tuple(
                    self._print(a) for a in expr.args)

            def _print__rdiv(self, expr):
                return '_rdiv(%s, %s)' % tuple(
                    self._print(a) for a in expr.args)

            ## The selection primitives, same treatment.  A missing entry
            ## here fails ONLY for a model that uses `var()` -- which is
            ## every production model -- because the eager path prints
            ## through sympy's own NumPyPrinter and this one prints its
            ## own source.
            def _print_maxc(self, expr):
                return 'maxc(%s, %s)' % tuple(
                    self._print(a) for a in expr.args)

            def _print_minc(self, expr):
                return 'minc(%s, %s)' % tuple(
                    self._print(a) for a in expr.args)

            def _print__step(self, expr):
                return '_step(%s, %s)' % tuple(
                    self._print(a) for a in expr.args)

            ## The idtmod wrap.  It reaches this printer only through a
            ## chained model -- and until `var()` definitions were
            ## searched for state operators, no chained model could HAVE
            ## an idtmod, so the gap could not show.  Fixing that
            ## uncovered this one immediately: `var(idtmod(...))`
            ## compiled the state correctly and then died here.
            def _print__wrapfloor(self, expr):
                return '_wrapfloor(%s)' % self._print(expr.args[0])

        return _P

    @staticmethod
    def _minmax(printer, expr, op):
        """`Min`/`Max` as nested binary calls.

        sympy prints these as `functools.reduce(numpy.minimum, [...])`,
        which allocates a list and enters `reduce` for what is almost
        always a two-argument call.  Measured at 18% of the Jacobian's
        runtime on the surface-potential MOSFET.
        """
        parts = [printer.doprint(a) for a in expr.args]
        out = parts[0]
        for nxt in parts[1:]:
            out = 'numpy.%s(%s, %s)' % (op, out, nxt)
        return out

    @staticmethod
    def _piecewise(printer, expr):
        import sympy as _s
        args = list(expr.args)
        ## A trailing `(value, True)` is the fall-through; without one
        ## sympy leaves the result undefined outside the conditions, and
        ## `nan` is the honest rendering of that.
        if args and args[-1].cond == _s.true:
            out = printer.doprint(args[-1].expr)
            args = args[:-1]
        else:
            out = 'numpy.nan'
        for e, c in reversed(args):
            out = 'numpy.where(%s, %s, %s)' % (printer.doprint(c),
                                               printer.doprint(e), out)
        return out

    def doprint(self, expr):
        return self._p.doprint(expr)


class CUnsupported(Exception):
    """The C printer met an expression it has no faithful rendering for.

    Raised at compile time, never at call time: a class whose chain
    cannot be printed to C compiles on the numpy path and says so in
    `cls._hdl_backend_status`.
    """


## What the numpy printer's fully-qualified names become in C.  Every
## `numpy.<name>` the emitted Python calls must appear here; one that
## does not raises `CUnsupported` at print time rather than emitting a
## call the C compiler would resolve to something else.  `numpy.abs` is
## `fabs`, `numpy.sign` the prelude's NaN-propagating `_sign`, and the
## constants are C99's.
_C_NAMES = {
    'numpy.exp': 'exp', 'numpy.log': 'log', 'numpy.sqrt': 'sqrt',
    'numpy.log1p': 'log1p', 'numpy.expm1': 'expm1',
    'numpy.sin': 'sin', 'numpy.cos': 'cos', 'numpy.tan': 'tan',
    'numpy.arcsin': 'asin', 'numpy.arccos': 'acos', 'numpy.arctan': 'atan',
    'numpy.arctan2': 'atan2', 'numpy.hypot': 'hypot',
    'numpy.sinh': 'sinh', 'numpy.cosh': 'cosh', 'numpy.tanh': 'tanh',
    'numpy.floor': 'floor', 'numpy.ceil': 'ceil',
    'numpy.abs': 'fabs', 'numpy.sign': '_sign',
    'numpy.minimum': '_npmin', 'numpy.maximum': '_npmax',
    'numpy.pi': 'M_PI', 'numpy.e': 'M_E',
    'numpy.inf': 'INFINITY', 'numpy.nan': 'NAN',
    ## The pycode printer emits the BUILTIN for these (numpy accepts
    ## the builtin protocol), so they arrive without the numpy prefix.
    'abs': 'fabs', 'sign': '_sign', 'math.floor': 'floor',
}

_C_RELOPS = {'==': '==', '!=': '!=', '<': '<', '<=': '<=', '>': '>',
             '>=': '>='}


class _CChainPrinter(_ChainPrinter):
    """`_ChainPrinter`'s structure, C99's spelling.

    Derives from the SAME `NumPyPrinter` subclass `_ChainPrinter` builds,
    so `Add`, `Mul`, the sign handling, the numerator/denominator split,
    parenthesisation and `Piecewise` nesting are the code the numpy path
    runs.  What is overridden is the leaf syntax: numbers carry a decimal
    point (C's `1/2` is integer division), symbols are looked up in
    `symmap` (unknowns as `x[k]`, parameters as `p[k]`, intermediates as
    `L_<name>`), comparisons and the Boolean connectives become the
    prelude's `_land`/`_lor`, `Piecewise` becomes `_sel(c, a, b)` --
    a CALL, so both arms are evaluated before the selection exactly as
    `numpy.where(c, a, b)` evaluates them -- and every `numpy.<f>` is
    mapped through `_C_NAMES`.

    `Pow` is where fidelity is decided.  numpy's scalar `x ** n` is libm
    `pow(x, n)` for every `n` (`benchmarks/backend_spike.py` measured
    it; glibc's `pow(x, 2)` differs from `x*x` in 0.07% of arguments),
    so the printer emits a real `pow` and the compiler is told not to
    fold it (`-fno-builtin-pow`).  ONE exception, and it is a fact about
    numpy rather than about C: `numpy.where` returns a 0-d ARRAY, and
    `ndarray.__pow__` with exponent 2 is `numpy.square` -- `x*x`, not
    `pow`.  So a square whose base is a `Piecewise`, or a chain symbol
    bound to one (`array_syms`), is printed as `_sq`.  Exponents -1 and
    1/2 take the same fast path in numpy, and land on `1.0/x` and
    `sqrt` here, which agree with `pow` bitwise anyway.
    """

    def __init__(self, symmap, array_syms=()):
        base = _ChainPrinter._printer_class()
        array_syms = frozenset(array_syms)
        import sympy as _s

        class _CP(base):
            def _module_format(self, fqn, register=True):
                try:
                    return _C_NAMES[fqn]
                except KeyError:
                    raise CUnsupported('no C rendering of %s' % fqn)

            def emptyPrinter(self, expr):
                raise CUnsupported('no C rendering of %s (%s)'
                                   % (type(expr).__name__, expr))

            _print_not_supported = emptyPrinter

            def _print_Symbol(self, expr):
                try:
                    return symmap[expr.name]
                except KeyError:
                    raise CUnsupported('symbol %s is not an argument or a '
                                       'chain definition' % expr.name)

            def _print_Indexed(self, expr):
                return '%s[%d]' % (expr.args[0], int(expr.args[1]))

            def _print_Integer(self, expr):
                return '%d.0' % expr.p

            def _print_Rational(self, expr):
                return '%d.0/%d.0' % (expr.p, expr.q)

            _print_Half = _print_Rational

            def _print_BooleanTrue(self, expr):
                return '1.0'

            def _print_BooleanFalse(self, expr):
                return '0.0'

            def _print_NaN(self, expr):
                return 'NAN'

            _print_ComplexInfinity = _print_NaN

            def _print_Infinity(self, expr):
                return 'INFINITY'

            def _print_NegativeInfinity(self, expr):
                return '(-INFINITY)'

            def _hprint_Pow(self, expr, rational=False, sqrt=None):
                from sympy.printing.precedence import precedence
                PREC = precedence(expr)
                base_ = expr.base
                if expr.exp == _s.S.Half:
                    return 'sqrt(%s)' % self._print(base_)
                if -expr.exp == _s.S.Half:
                    return '1.0/sqrt(%s)' % self._print(base_)
                if expr.exp is _s.S.NegativeOne:
                    return '1.0/%s' % self.parenthesize(base_, PREC,
                                                        strict=False)
                is_array = isinstance(base_, _s.Piecewise) or \
                    (isinstance(base_, _s.Symbol) and base_ in array_syms)
                if is_array and expr.exp == 2:
                    return '_sq(%s)' % self._print(base_)
                if is_array and expr.exp == -1:
                    ## `ndarray.__pow__(-1.0)` is `reciprocal`, i.e. a
                    ## correctly rounded division, where glibc's `pow`
                    ## need not be.
                    return '1.0/(%s)' % self._print(base_)
                return 'pow(%s, %s)' % (self._print(base_),
                                        self._print(expr.exp))

            def _print_Relational(self, expr):
                op = _C_RELOPS.get(expr.rel_op)
                if op is None:
                    raise CUnsupported('relation %s' % expr.rel_op)
                return '(%s %s %s)' % (self._print(expr.lhs), op,
                                       self._print(expr.rhs))

            def _print_And(self, expr):
                return self._nest('_land', expr.args)

            def _print_Or(self, expr):
                return self._nest('_lor', expr.args)

            def _print_Not(self, expr):
                return '_lnot(%s)' % self._print(expr.args[0])

            def _nest(self, fn, args):
                out = self._print(args[0])
                for a in args[1:]:
                    out = '%s(%s, %s)' % (fn, out, self._print(a))
                return out

            def _print_Min(self, expr):
                return self._nest('_npmin', expr.args)

            def _print_Max(self, expr):
                return self._nest('_npmax', expr.args)

            def _print_Piecewise(self, expr):
                args = list(expr.args)
                if args and args[-1].cond == _s.true:
                    out = self._print(args[-1].expr)
                    args = args[:-1]
                else:
                    out = 'NAN'
                for e, c in reversed(args):
                    out = '_sel(%s, %s, %s)' % (self._print(c),
                                                self._print(e), out)
                return out

            def _print_maxc(self, expr):
                return self._nest('_npmax', expr.args)

            def _print_minc(self, expr):
                return self._nest('_npmin', expr.args)

            def _print__wrapfloor(self, expr):
                return 'floor(%s)' % self._print(expr.args[0])

        self._p = _CP()


## The kernel's runtime contract in C, one definition per primitive,
## kept beside `_KERNEL_NUMPY` for the same reason the numpy forms are
## kept beside their sympy classes: so the two cannot drift apart.
## Every helper is a FUNCTION.  C evaluates all arguments before the
## call, so `_sel(c, a, b)` computes both arms exactly as Python does for
## `numpy.where(c, a, b)`, and only then picks; a `?:` with the arms
## inline would skip one, and `a*c + b*(1-c)` would turn a losing
## infinite arm into NaN (`_maxc_numpy` records why that matters).
## `_npmax`/`_npmin` are numpy's `maximum`/`minimum` verbatim: NaN in
## either argument wins.  `c != 0.0` is numpy's truth test -- nonzero,
## and NaN, are true.
_KERNEL_C = r"""
#include <math.h>
static inline double _sel(double c, double a, double b) { return (c != 0.0) ? a : b; }
static inline double _npmax(double a, double b) { return (a >= b || a != a) ? a : b; }
static inline double _npmin(double a, double b) { return (a <= b || a != a) ? a : b; }
static inline double _step(double a, double b) { return 1.0 * (a >= b); }
static inline double _rdiv(double b, double e) { return b / (b * b + e * e); }
static inline double _recip2(double b, double e) { return 1.0 / (b * b + e * e); }
static inline double _sq(double a) { return a * a; }
static inline double _sign(double a) { return (a != a) ? a : (a > 0.0) ? 1.0 : (a < 0.0) ? -1.0 : 0.0; }
static inline double _land(double a, double b) { return (a != 0.0) && (b != 0.0); }
static inline double _lor(double a, double b) { return (a != 0.0) || (b != 0.0); }
static inline double _lnot(double a) { return (a == 0.0); }
"""

#: The exported symbol of every compiled chain; each lives in its own
#: shared object, loaded RTLD_LOCAL, so the name cannot clash.
_C_ENTRY = 'hdl_fn'

#: Whether `generate_code` prints the C text alongside the numpy source
#: for every x-taking chain.  On, whatever backend is selected --
#: measured at +19% on a COLD library compile (8.4 -> 10.0 s for the
#: 37-class library), paid once and then amortised by the compile
#: cache, and in exchange `PYCIRCUIT_HDL_BACKEND=c` on a warm cache
#: re-runs no sympy at all.  Nothing is COMPILED unless the C backend
#: is chosen.
EMIT_C_SOURCE = True

#: Process-wide backend choice: `None` reads `$PYCIRCUIT_HDL_BACKEND`
#: (default `'numpy'`); `'numpy'` or `'c'` overrides it.  Read at class
#: compile time; `set_backend` changes it for classes built afterwards
#: and re-attaches any class it is handed.  Per-class override: a class
#: attribute `hdl_backend = 'c'` (inherited like any attribute).
BACKEND = None

_BACKENDS = ('numpy', 'c')


def _backend_requested(cls):
    """Which backend `cls` asked for: the class attribute, else the
    module flag, else the environment, else numpy.  An unknown name
    raises -- a typo in `PYCIRCUIT_HDL_BACKEND` must not quietly run
    numpy while the user believes C is on."""
    which = getattr(cls, 'hdl_backend', None)
    if which is None:
        which = BACKEND
    if which is None:
        which = os.environ.get('PYCIRCUIT_HDL_BACKEND', '').strip() or 'numpy'
    if which not in _BACKENDS:
        raise ValueError('unknown HDL backend %r (one of %s; check '
                         'PYCIRCUIT_HDL_BACKEND)' % (which,
                                                     '/'.join(_BACKENDS)))
    return which


def set_backend(which, cls=None):
    """Select the evaluation backend: `'numpy'` (the default) or `'c'`.

    With `cls` given, re-attaches that one class immediately and pins it
    (sets `cls.hdl_backend`); without, sets the process-wide default for
    classes compiled from then on.  `which=None` removes the pin (or the
    process default) so the environment decides again.
    `cls._hdl_backend_status` afterwards says what actually happened --
    `'c'`, or `'numpy (<why not>)'`.
    """
    if which is not None and which not in _BACKENDS:
        raise ValueError('unknown HDL backend %r' % (which,))
    global BACKEND
    if cls is None:
        BACKEND = which
        return
    if which is None:
        try:
            del cls.hdl_backend
        except AttributeError:
            pass
    else:
        cls.hdl_backend = which
    from pycircuit.circuit import _hdl_cbackend
    _hdl_cbackend.attach(cls, cls._hdl_info)
    ## A collapsing model RUNS as a compiled variant subclass with its
    ## own `_hdl_info` (`_collapse_variant`); a pin on the base must
    ## follow, or `set_backend('c', cls)` would silently leave every
    ## existing instance on numpy.  Variants built later inherit
    ## `hdl_backend` and attach themselves at creation.
    base = getattr(cls, '_hdl_collapse_base', cls)
    for variant in (base.__dict__.get('_hdl_mask_classes') or {}).values():
        if variant is not cls:
            _hdl_cbackend.attach(variant, variant._hdl_info)


def _render_c(stmts, cells, args, xsyms):
    """The C function for a chain: `(stmts, cells)` as `_chain_compile`
    collected them, `args` the numpy signature `[x, *trailing]`.

    Returns `(text, shape, layout)`: the function body without the
    prelude, the output array's shape, and `(n_p, t_index)` -- the
    length of the packed trailing-argument vector `p` and the index in
    it of the temperature, which the element writes per call.
    """
    trailing = list(args[1:])
    symmap = {}
    for k, xs in enumerate(xsyms or ()):
        symmap[xs.name] = 'x[%d]' % k
    for k, a in enumerate(trailing):
        symmap[a.name] = 'p[%d]' % k
    array_syms = set()
    for sym, expr in stmts:
        symmap[sym.name] = 'L_' + sym.name
        if isinstance(expr, sympy.Piecewise):
            array_syms.add(sym)
    printer = _CChainPrinter(symmap, array_syms)
    lines = ['void %s(const double *x, const double *p, double *out) {'
             % _C_ENTRY]
    for sym, expr in stmts:
        lines.append('  const double L_%s = %s;'
                     % (sym.name, printer.doprint(expr)))
    if cells and isinstance(cells[0], (list, tuple)):
        ncol = len(cells[0])
        if any(len(r) != ncol for r in cells):
            raise CUnsupported('ragged output')
        shape = (len(cells), ncol)
        flat = [e_ for r in cells for e_ in r]
    else:
        shape = (len(cells),)
        flat = list(cells)
    for k, e_ in enumerate(flat):
        lines.append('  out[%d] = %s;' % (k, printer.doprint(e_)))
    lines.append('}')
    t_index = [k for k, a in enumerate(trailing) if a == TEMP]
    layout = (len(trailing), t_index[0] if t_index else None)
    return '\n'.join(lines) + '\n', shape, layout


def _leaves(o):
    """Free symbols of an expression, or of a nested list of them."""
    if isinstance(o, (list, tuple)):
        out = set()
        for e_ in o:
            out |= _leaves(e_)
        return out
    return set(getattr(o, 'free_symbols', ()))


def _var_name(sym):
    """The name the model author gave a `var()`, from its symbol."""
    return re.sub(r'^_v\d+_', '', sym.name)


def _limit_par_fn(expr, kind, chain_defs, xsyms, xlabels, args, modules,
                  x=None):
    """Compile one `$limit` parameter expression.

    `args` is the ordinary `params + [T] + givenness` signature, and an
    expression that mentions nothing else is lambdified over it exactly
    as before.  An `_IntermediateSymbol` from `var()` is resolved against
    the chain -- READ, not inlined, through `_chain_compile`, the same
    path the PCNR detector's `VT` and `IS` take -- so a limiter parameter
    written in terms of a hundred-deep MOSFET chain costs what the chain
    costs rather than exponentially more.

    **A PARAMETER MAY READ THE SOLUTION, AND IS THEN EVALUATED AT THE
    LAST ACCEPTED ITERATE** (2026-08-25).  This used to be refused by
    name -- "a limiter runs BEFORE the device is, so its parameters
    cannot depend on the solution" -- and the refusal was correct about
    the ORDER and wrong about the CONCLUSION.  A limiter is handed `x0`,
    the last accepted point, precisely so that it can measure the new
    step against it; a parameter evaluated there is as well-defined as
    `vold` is.  That is also exactly SPICE's semantics: `mos1load.c`
    passes `fetlim` a `von` recomputed from the PREVIOUS iterate's bulk
    bias, not a card constant.  Two production models had already paid
    for the refusal -- EKV limited a body-biased device against a
    threshold 565 mV low, and the self-heating BJT carried a second,
    ambient-temperature copy of its saturation current for the limiter
    alone.

    The returned callable carries `_wants_x`: when True its signature is
    `(x0, *args)` and the generated `limit()` supplies the element's own
    last-accepted sub-vector.  Time is still refused -- a limiter has no
    `t` to offer.
    """
    reach = _chain_prune(chain_defs, [expr])
    xset = set(xsyms)
    what = _LIMIT_FN.get(kind, kind)
    late = sorted(_var_name(s_) for s_, e_ in reach if e_.has(TIME))
    if late or expr.has(TIME):
        raise ValueError(
            "%s's parameters cannot depend on time%s."
            % (what, (', and var(%s) does' % ') / var('.join(
                repr(v) for v in late)) if late else ''))
    touches_x = bool(expr.free_symbols & xset) or any(
        e_.free_symbols & xset for _s, e_ in reach)
    mods = dict(_KERNEL_NUMPY, _wrapfloor=np.floor)
    if touches_x:
        inner = _chain_compile(
            reach, [expr], [x] + args, modules_map=mods,
            unpack=[(xs.name, 'x[%d]' % k) for k, xs in enumerate(xsyms)])
        return _first_of(inner, wants_x=True)
    if not reach:
        fn = sympy.lambdify(args, expr, modules=modules)
        fn._wants_x = False
        return fn
    inner = _chain_compile(reach, [expr], args, modules_map=mods)
    return _first_of(inner, wants_x=False)


def _pcnr_declared_route(limits, ivec, chain_defs, xsyms, xlabels, paramsyms,
                         given_syms, states, vbranches):
    """Vector PCNR's detector: the DECLARED-PROBE route (roadmap sec. 15,
    Stages 1 and 2).  Returns ``(spec, refusal)``, exactly one of them set.

    A device with ``$limit`` probes qualifies when every resistive
    current reaches the solution ONLY through those probes: the probe
    vector is the device's block of limited unknowns, and the limiter it
    declared is the law -- no exponential-scale test is run, because a
    pnjlim on a probe already says what the current does there, and a
    fetlim/limvds probe says the device asked for that law on that
    branch.

    The check is the legacy detector's substitution generalised from one
    branch to a probe SET.  The probes are edges over the device's rows;
    each connected component gets a root potential and every other row is
    written as ``root +/- v_j`` along a spanning TREE, so that ``x_a - x_b``
    becomes ``v_j`` wherever a probe's branch voltage is read, and any
    x-symbol that survives (a root, a row no probe reaches) names a
    voltage the current reads that no probe limits.

    **A probe that closes a cycle is REDUNDANT, not refused** (Stage 2).
    SPICE's own MOSFET set -- `(g,s)`, `(d,s)`, `(b,s)`, `(b,d)` -- has
    `(b,d) = (b,s) - (d,s)`, and a fourth unknown for it would be a
    linear combination of the other three, not a quantity of its own.
    The device's PCNR unknowns are the tree's probes (``m = t - 1`` per
    connected component); the redundant probe's branch voltage is read
    off the tree (``coeffs`` is its signed combination of the unknowns),
    and its LIMITER is still applied: `pcnr.limit_block` resolves all of
    a device's laws over the tree by the same rule `limit_together`'s
    write-back uses.  Which probe is the redundant one is decided by
    declaration order (the first to close a cycle) and only affects
    which unknowns are carried, never which laws are honoured.

    Charge, noise and sources are not looked at: charge stays in the MNA
    block (`pcnr.augmented_system`).
    """
    if states or vbranches:
        return None, ('%s -- the device carries %s'
                      % ('a generated state' if states
                         else 'a branch-current unknown',
                         '%d state(s)' % len(states) if states
                         else '%d V-contributed branch(es)' % len(vbranches)))
    all_probes = [(int(k[0]), int(k[1])) for k, _kind, _p, _g in limits]
    all_kinds = [kind for _k, kind, _p, _g in limits]
    vsyms, probes, kinds, spec_idx = [], [], [], []
    redundant = []           # (spec index, (a, b), kind, coeffs over vsyms)

    ## The spanning tree, by union of root-relative potentials.  A
    ## probe's own symbol is allocated only when it becomes a tree edge.
    pot, root_of = {}, {}
    for j, (a, b) in enumerate(all_probes):
        if a in pot and b in pot and root_of[a] == root_of[b]:
            expr = sympy.expand(pot[a] - pot[b])
            coeffs = [int(expr.coeff(vs)) for vs in vsyms]
            if sympy.expand(expr - sum(c_ * vs for c_, vs
                                       in zip(coeffs, vsyms))) != 0:
                raise AssertionError(       # pragma: no cover
                    'redundant probe (%s,%s) is not a +/-1 combination of '
                    'the tree: %s' % (xlabels[a], xlabels[b], expr))
            redundant.append((j, (a, b), all_kinds[j], coeffs))
            continue
        vj = sympy.Symbol('_pcnr_v%d' % len(vsyms), real=True)
        vsyms.append(vj)
        probes.append((a, b))
        kinds.append(all_kinds[j])
        spec_idx.append(j)
        if a in pot and b in pot:
            ## Merge b's tree into a's: root_b = pot[a] - vj - (pot[b] - root_b)
            rb_ = root_of[b]
            rb_sym = xsyms[rb_]
            new_root = pot[a] - vj - (pot[b] - rb_sym)
            for node in [n_ for n_, r_ in root_of.items() if r_ == rb_]:
                pot[node] = pot[node].subs(rb_sym, new_root)
                root_of[node] = root_of[a]
        elif a in pot:
            pot[b] = pot[a] - vj
            root_of[b] = root_of[a]
        elif b in pot:
            pot[a] = pot[b] + vj
            root_of[a] = root_of[b]
        else:
            pot[a] = xsyms[a]
            pot[b] = xsyms[a] - vj
            root_of[a] = root_of[b] = a
    sd = {xsyms[k]: e_ for k, e_ in pot.items() if e_ != xsyms[k]}
    xset = set(xsyms)

    def _sub_x(e_):
        e2 = e_.subs(sd)
        return sympy.expand(e2) if (e2.free_symbols & xset) else e2

    defs_v = [(sym, _sub_x(e_)) for sym, e_ in _chain_prune(chain_defs, ivec)]
    i_exprs = [_sub_x(e_) for e_ in ivec]
    defsyms = {sym for sym, _e in defs_v}
    ## `given_syms` are handed to `pcnr_i` beside the parameters (the
    ## `$given:<name>` flags in the device's param dict), so a chain that
    ## reads a `$param_given` -- Level 1's `gamma` -- qualifies.
    allowed = set(paramsyms) | {TEMP} | set(vsyms) | defsyms | set(given_syms)
    for label, e_ in [('var(%s)' % _var_name(sym), e2)
                      for sym, e2 in defs_v] + \
            [('I(%s)' % xlabels[k], e2) for k, e2 in enumerate(i_exprs)]:
        fs = set(e_.free_symbols)
        if fs & xset:
            return None, ('%s reads %s, which no $limit probe limits -- under '
                          'vector PCNR every resistive current must reach the '
                          'solution only through the declared probes'
                          % (label, ', '.join(sorted(
                              xlabels[xsyms.index(q)] for q in fs & xset))))
        extra = fs - allowed
        if extra:
            return None, ('%s reads %s, which pcnr_i is not handed'
                          % (label, ', '.join(sorted(str(q) for q in extra))))
    return dict(vsyms=vsyms, probes=probes, kinds=kinds, spec_idx=spec_idx,
                redundant=redundant, defs=defs_v, i_exprs=i_exprs,
                t=len(xsyms)), None


def _chain_prune(defs, outputs):
    """The definitions `outputs` actually reach, in `defs` order.

    A model's chain is shared by every compiled vector -- `i`, `q`, `u`,
    the AC stamp -- but each of those reads only part of it, and the PCNR
    detector reads a still smaller part (one contribution's).  Walking
    the DAG backwards from the outputs is linear in the number of
    definitions; it is the alternative to INLINING the chain, which is
    the one thing `var()` exists to prevent.
    """
    defmap = {sym: e_ for sym, e_ in defs}
    wanted, stack = set(), []
    for o in outputs:
        stack.extend(_leaves(o))
    while stack:
        sym = stack.pop()
        if sym in wanted or sym not in defmap:
            continue
        wanted.add(sym)
        stack.extend(defmap[sym].free_symbols)
    return [(sym, e_) for sym, e_ in defs if sym in wanted]


def _chain_d(expr, v, defset, dmap):
    """One forward-accumulation step: d(expr)/dv through the chain.

    `dmap` maps each already-processed definition symbol to the SYMBOL
    holding its own derivative, so the result stays a small expression
    referring to those symbols rather than the flattened derivative --
    which is the same trick `_chain_compile`'s Jacobian uses, and the
    reason the work is linear in the number of definitions instead of
    exponential in their nesting depth.
    """
    g = sympy.diff(expr, v)
    for d in expr.free_symbols:
        if d in defset:
            g = g + sympy.diff(expr, d) * dmap[d]
    return g


def _chain_forward_d(defs, v):
    """Forward-accumulate d/dv over a whole chain.

    Returns `(ddefs, dmap)`: `ddefs` is the derivative chain in
    dependency order, ready to be concatenated after `defs` and handed to
    `_chain_compile`, and `dmap[sym]` names the symbol carrying
    d(sym)/dv.
    """
    defset = {sym for sym, _ in defs}
    ddefs, dmap = [], {}
    for sym, expr in defs:
        gs = sympy.Symbol('_dpcnr_%s' % sym.name, real=True)
        ddefs.append((gs, _chain_d(expr, v, defset, dmap)))
        dmap[sym] = gs
    return ddefs, dmap


def _chain_compile(defs, outputs, args, want_jacobian_of=None, xsyms=None,
                   modules_map=None, unpack=(), emit_c=False):
    """Compile a let-chain into one Python function.

    `defs` is `[(symbol, expr)]` in dependency order; `outputs` the list
    of expressions to return; `args` the function's parameters.  With
    `want_jacobian_of` set, the returned function yields the Jacobian of
    those outputs with respect to `xsyms`, obtained by **forward
    accumulation over the chain**: each definition's gradient is
    expressed in terms of the gradients of the definitions it actually
    mentions, so the work is linear in the number of definitions rather
    than exponential in their nesting depth.

    The numpy source is always produced and is the reference: it is what
    `explain()` shows, what the compile cache stores, and what the C
    backend is tested against.  With `emit_c` the SAME statement list is
    also printed to C (`_render_c`), and the function carries it as
    `_csrc` with `_cshape` and `_clayout`; `_hdl_cbackend` compiles and
    binds it when the backend is selected.  C is only emitted for the
    `(x, *trailing)` signature -- `unpack` non-empty -- because that is
    the one whose per-call cost matters and the one the packed
    parameter vector fits.  A chain the C printer cannot render
    (`CUnsupported`) leaves the attributes off and the reason in
    `_creason`; the numpy function is unaffected.
    """
    printer = _ChainPrinter()

    ## PRUNE to what the outputs actually reach.  Emitting the whole
    ## chain would put symbols in the body that this signature cannot
    ## supply (the source sub-chain reads TIME, which `i` is not given),
    ## and would differentiate definitions nothing downstream uses.
    defs = _chain_prune(defs, outputs)

    ## The statements `(sym, expr)` in emission order and the output
    ## cells (a list, or a list of rows), collected ONCE so that every
    ## printer renders the same list.
    stmts = []
    if want_jacobian_of is None:
        stmts.extend(defs)
        cells = list(outputs)
    else:
        ## Values first -- the gradients reference them.
        stmts.extend(defs)
        ## Forward accumulation.  The chain rule is applied only along
        ## the edges that exist: `expr.free_symbols` names exactly the
        ## upstream definitions this one reads, so the total work is
        ## proportional to the number of DAG EDGES, not to the number of
        ## (definition, definition) pairs.
        defset = {sym for sym, _ in defs}
        dsym, dexpr = {}, {}
        for sym, expr in defs:
            parents = [d for d in expr.free_symbols if d in defset]
            partials = {d: sympy.diff(expr, d) for d in parents}
            for j, xj in enumerate(xsyms):
                g = sympy.diff(expr, xj)
                for d in parents:
                    dj = dsym[(d, j)]
                    if dexpr[(d, j)].is_zero:
                        continue      # that path carries no sensitivity
                    g = g + partials[d] * dj
                gs = sympy.Symbol('_d_%s_%d' % (sym.name, j), real=True)
                dsym[(sym, j)] = gs
                dexpr[(sym, j)] = g
                if not g.is_zero:
                    stmts.append((gs, g))
        cells = []
        for out in outputs:
            parents = [d for d in out.free_symbols if d in defset]
            partials = {d: sympy.diff(out, d) for d in parents}
            row = []
            for j, xj in enumerate(xsyms):
                g = sympy.diff(out, xj)
                for d in parents:
                    if dexpr[(d, j)].is_zero:
                        continue
                    g = g + partials[d] * dsym[(d, j)]
                row.append(g)
            cells.append(row)

    lines, body = [], []
    for name, code in unpack:
        body.append('    %s = %s' % (name, code))
    for sym, expr in stmts:
        body.append('    %s = %s' % (sym.name, printer.doprint(expr)))

    def render(o):
        if isinstance(o, (list, tuple)):
            return '[%s]' % ', '.join(render(e_) for e_ in o)
        return printer.doprint(o)
    ret = render(cells)

    lines.append('def _f(%s):' % ', '.join(a.name if hasattr(a, 'name')
                                           else str(a) for a in args))
    lines.extend(body)
    lines.append('    return %s' % ret)
    src = '\n'.join(lines)
    ns = _chain_namespace(modules_map)
    exec(compile(src, '<hdl-chain>', 'exec'), ns)
    fn = ns['_f']
    fn._src = src
    if emit_c and unpack and xsyms:
        try:
            csrc, shape, layout = _render_c(stmts, cells, args, xsyms)
        except CUnsupported as e:
            fn._creason = str(e)
        else:
            fn._csrc, fn._cshape, fn._clayout = csrc, shape, layout
    return fn


def _chain_namespace(modules_map):
    """The namespace a chain-compiled function runs in.

    Built in ONE place because two callers must agree on it: the
    compiler above, and the on-disk cache (`_hdl_cache`), which
    re-executes the stored source in a namespace rebuilt by this same
    function and checks, before recording anything, that every global
    the function loads is bound to the same object here as it was at
    compile time.
    """
    ns = dict(modules_map or {})
    import functools
    import numpy
    ## The DSL's own primitives, merged unconditionally.  They print as
    ## calls, so a namespace missing one is a `NameError` at call time
    ## from a model that compiled clean -- and `modules_map` comes from
    ## the CALLER, so relying on every call site to pass them is relying
    ## on every future call site too.  Set here, once, where the
    ## namespace is actually built.
    ns.setdefault('_wrapfloor', numpy.floor)
    for _k, _v2 in _KERNEL_NUMPY.items():
        ns.setdefault(_k, _v2)
    ## NumPyPrinter prints Min/Max as `functools.reduce(numpy.minimum, ...)`,
    ## so both names have to be in the namespace, not just numpy's.
    ns.setdefault('numpy', numpy)
    ns.setdefault('functools', functools)
    for k in ('sqrt', 'exp', 'log', 'sin', 'cos', 'tan', 'sinh', 'cosh',
              'tanh', 'floor', 'ceil', 'abs', 'sign', 'pi', 'select',
              'less', 'greater', 'logical_and', 'logical_or', 'nan',
              'inf', 'amin', 'amax', 'minimum', 'maximum', 'where'):
        ns.setdefault(k, getattr(numpy, k, None))
    return ns


def _first_of(fn, wants_x=None):
    """`fn`'s first output as a scalar-valued callable.

    A named wrapper rather than an inline lambda so the compile cache
    can see through it: the wrapper carries the wrapped function as
    `_hdl_inner` and, when given, the `_wants_x` flag `limit()` reads.
    """
    out = lambda *a_, _f=fn: _f(*a_)[0]      # noqa: E731
    out._hdl_inner = fn
    if wants_x is not None:
        out._wants_x = wants_x
    return out


def _params_of(self):
    ## Values from the RESOLVED iparv, plus the givenness flags from
    ## `ipar` -- givenness is a property of what the user wrote, and only
    ## `ipar` records that (see ParameterDict.update_values).
    vals = [getattr(self.iparv, name) for name in self._hdl_paramnames]
    ipar = self.ipar
    vals += [1.0 if ipar.is_given(nm) else 0.0
             for nm in self._hdl_given_names]
    return vals


def _epar_T(epar):
    return getattr(epar, 'T', 300.0)


def _symbolic_eval(self, which, x, epar):
    """Exact substitution for a symbolic toolkit -- see `sym_spec`."""
    spec = self._hdl_info['sym_spec']
    subs = dict(zip(spec['xsyms'], list(x)))
    subs.update(zip(spec['paramsyms'],
                    [getattr(self.iparv, nm)
                     for nm in self._hdl_paramnames]))
    subs[TEMP] = _epar_T(epar)
    ipar = self.ipar
    subs.update({sym: (1.0 if ipar.is_given(nm) else 0.0)
                 for sym, nm in zip(spec['given_syms'],
                                    self._hdl_given_names)})
    vec = [e.subs(subs) for e in spec[which]]
    if which in ('G', 'C'):
        n = spec['shape'][0]
        return sympy.Matrix(n, n, vec)
    return vec


def _collapse_mask_of(cls, paramvals):
    """Which collapses this parameter set takes, in declaration order."""
    conds = cls._hdl_info['collapse_conds']
    n_par = len(cls._hdl_paramnames)
    return tuple(bool(f(*paramvals[:n_par])) for f in conds)


def _collapse_variant(cls, mask):
    """The class compiled for `mask`, built once and cached.

    A SUBCLASS rather than a per-instance table of functions: the
    metaclass already compiles everything a class needs from its
    `analog`, so re-entering it with the mask baked in gets correct
    stamps, branches, PCNR participation, state metadata and pure forms
    with no second code path to keep in step.  Instances are then
    retargeted by assigning `__class__`.
    """
    base = getattr(cls, '_hdl_collapse_base', cls)
    cache = base.__dict__.get('_hdl_mask_classes')
    if cache is None:
        cache = {}
        base._hdl_mask_classes = cache
    if mask in cache:
        return cache[mask]
    ## `analog` must appear in the class body or the metaclass will not
    ## recompile, and `instparams` must be carried across or the
    ## instparams consistency check would reject the variant.
    var = BehaviouralMeta(
        '%s_collapse%s' % (base.__name__,
                           ''.join('1' if b else '0' for b in mask)),
        (base,),
        dict(analog=staticmethod(base.analog),
             instparams=list(base.instparams),
             _hdl_collapse_mask=mask,
             _hdl_collapse_base=base,
             __module__=base.__module__,
             __doc__='%s with collapses %r taken.' % (base.__name__, mask)))
    cache[mask] = var
    return var


def _args_of(self, epar):
    """The trailing argument list every compiled function expects:
    parameter values, then T, then the givenness flags."""
    n_par = len(self._hdl_paramnames)
    vals = _params_of(self)
    return vals[:n_par] + [_epar_T(epar)] + vals[n_par:]


def _dc(epar):
    return getattr(epar, 'analysis_kind', None) == 'dc'


class BehaviouralMeta(type):
    def __init__(cls, name, bases, dct):
        if 'analog' not in dct:
            return
        ## The signature checks run on every start; the compile itself
        ## is served from the on-disk cache when its key matches (see
        ## `_hdl_cache`), and `cls._hdl_cache_status` says which.
        _check_analog_declaration(cls)
        from pycircuit.circuit import _hdl_cache
        info = _hdl_cache.compiled_info(cls, generate_code)
        funcs = info['funcs']
        ## The evaluation backend (numpy by default; C when selected).
        ## Attached AFTER the compile so a cache hit can bind too, and
        ## before the methods below close over `funcs` -- they consult
        ## each function's `_hdl_c` at call time, so `set_backend` can
        ## re-attach without touching the class.
        from pycircuit.circuit import _hdl_cbackend
        _hdl_cbackend.attach(cls, info)
        cls._hdl_paramnames = info['paramnames']
        cls._hdl_given_names = info['given_names']
        cls._hdl_info = info
        cls.terminals = info['terminalnames']
        ## Branch objects at class level, like every element: the base
        ## Circuit.__init__ counts them into n; plus/minus resolve through
        ## the terminal mapping when the branch spans terminals, and to
        ## internal nodes otherwise.
        cls.branches = tuple(
            circuit.Branch(circuit.Node(p), circuit.Node(m), name=nm)
            for p, m, nm in info['branchpairs'])

        state_meta = info['state_meta']
        if state_meta['dc_pins']:
            cls.IC_KIND = 'state'

        def i(self, x, epar=defaultepar, params_tree=None):
            if info['chained']:
                ## The DC pin applies on the chained path too.  Until
                ## 2026-08-26 this returned `funcs['i']` unconditionally,
                ## so a `var()` model with `idt(expr, ic)` or `idtmod(...,
                ## ic, ...)` kept its state equation at DC and the
                ## operating point was singular; `IdtmodHdl` is flat and
                ## never showed it.  Found by `VcoHdl` (fifth batch).
                f = funcs['i_dc'] if _dc(epar) and state_meta['dc_pins'] \
                    else funcs['i']
                ck = f.__dict__.get('_hdl_c')
                if ck is not None:
                    return ck(self, x, epar)
                return np.asarray(f(x, *_args_of(self, epar)), dtype=float)
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'i', x, epar)
            f = funcs['i_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['i']
            return f(x, *_args_of(self, epar))

        const_G, const_C = info['const_G'], info['const_C']

        def update(self, subject):
            ## The iparv observer every element implements (R.update,
            ## elements.py): parameters changed, so any cached constant
            ## stamp is stale.  This is what makes the cache below a bare
            ## attribute test rather than a per-call dict key -- and it
            ## also means a late-resolved parameter expression reaches the
            ## generated code, since values are read from iparv at call
            ## time and the cache is dropped whenever they move.
            self.__dict__.pop('_hdl_Gc', None)
            self.__dict__.pop('_hdl_Cc', None)
            ## The C backend's packed parameter vector is a cache of
            ## iparv too.
            self.__dict__.pop('_hdl_cp', None)
            ## A parameter that gates a collapse decides the element's
            ## SIZE, and the size was baked into the circuit's node map
            ## when this instance was built.  Changing it afterwards
            ## would leave the map describing a different element, so say
            ## so rather than solve a system nobody wrote.
            seen = getattr(self, '_hdl_collapse_seen', None)
            if seen is not None:
                now = _collapse_mask_of(type(self), _params_of(self))
                if now != seen:
                    raise ValueError(
                        '%s: a parameter that gates a node collapse changed '
                        'after the element was built (%r -> %r). That '
                        'changes the number of unknowns, which the circuit '
                        'node map already fixed. Build a new instance '
                        'instead of re-parameterising this one.'
                        % (type(self).__name__, seen, now))

        def G(self, x, epar=defaultepar, params_tree=None):
            if info['chained']:
                f = funcs['G_dc'] if _dc(epar) and state_meta['dc_pins'] \
                    else funcs['G']
                ck = f.__dict__.get('_hdl_c')
                if ck is not None:
                    return ck(self, x, epar)
                return np.asarray(f(x, *_args_of(self, epar)), dtype=float)
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'G', x, epar)
            dc = _dc(epar) and state_meta['dc_pins']
            if const_G and not dc:
                ## Constant stamp: computed once and handed back by
                ## reference, as a hand-written linear element does
                ## (benchmarks/hdl_overhead.py measured the gap this
                ## closes).  Dropped by update() on any parameter change.
                cached = self.__dict__.get('_hdl_Gc')
                if cached is None:
                    cached = np.asarray(funcs['G'](x, *_args_of(self, epar)))
                    self.__dict__['_hdl_Gc'] = cached
                return cached
            f = funcs['G_dc'] if dc else funcs['G']
            return np.asarray(f(x, *_args_of(self, epar)))

        def q(self, x, epar=defaultepar, params_tree=None):
            if info['chained']:
                ck = funcs['q'].__dict__.get('_hdl_c')
                if ck is not None:
                    return ck(self, x, epar)
                return np.asarray(funcs['q'](x, *_args_of(self, epar)),
                                  dtype=float)
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'q', x, epar)
            return funcs['q'](x, *_args_of(self, epar))

        def C(self, x, epar=defaultepar, params_tree=None):
            if info['chained']:
                ck = funcs['C'].__dict__.get('_hdl_c')
                if ck is not None:
                    return ck(self, x, epar)
                return np.asarray(funcs['C'](x, *_args_of(self, epar)),
                                  dtype=float)
            if getattr(self.toolkit, 'symbolic', False):
                return _symbolic_eval(self, 'C', x, epar)
            if const_C:
                cached = self.__dict__.get('_hdl_Cc')
                if cached is None:
                    cached = np.asarray(funcs['C'](x, *_args_of(self, epar)))
                    self.__dict__['_hdl_Cc'] = cached
                return cached
            return np.asarray(funcs['C'](x, *_args_of(self, epar)))

        def CY(self, x, w=0, epar=defaultepar):
            f_hz = np.abs(w) / (2 * np.pi)
            return np.asarray(funcs['CY'](x, *_args_of(self, epar), f_hz))

        def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
            if analysis == 'ac':
                ## ONLY `ac_stim` terms drive a small-signal analysis.
                ## The DC/transient source terms must NOT come through --
                ## injecting a device's bias constants as an AC drive is
                ## the classic ABM bias-leak (sec. 2.4).
                if not info['has_ac']:
                    return np.zeros(self.n)
                return np.asarray(funcs['uac'](*_args_of(self, epar)),
                                  dtype=complex)
            f = funcs['u_dc'] if _dc(epar) and state_meta['dc_pins'] \
                else funcs['u']
            return np.asarray(f(t, *_args_of(self, epar)), dtype=float)

        def dudt(self, t=0.0, epar=defaultepar, analysis=None,
                 params_tree=None):
            return np.asarray(funcs['dudt'](t, *_args_of(self, epar)),
                              dtype=float)

        def state_ic(self):
            out = []
            for k, icf in state_meta['dc_pins'].items():
                out.append((k, float(icf(*_params_of(self)))))
            return out

        def periodic_states(self):
            out = []
            for k, mf, of in state_meta['periodic']:
                m = float(mf(*_params_of(self)))
                if m > 0:
                    out.append((k, m, float(of(*_params_of(self)))))
            return out

        ## An `if/else`, NOT an early `return` -- and that is a FIX, not a
        ## refactor.  `pure_spec is None` is exactly the let-chain path
        ## (`var()`), and the `return` that used to stand here jumped over
        ## everything below it: a chained model got no `@cross` events and
        ## no `$limit` limiter, silently, however loudly it declared them.
        ## `explain()` still listed both, because it reads `info` and not
        ## the class.  Any production compact model is chained -- that is
        ## what `var()` is for -- so the two features were unavailable to
        ## precisely the models they exist for.
        if info['pure_spec'] is None:
            cls.linear = False
            cls.i, cls.G, cls.q, cls.C, cls.CY = i, G, q, C, CY
            cls.u, cls.dudt, cls.update = u, dudt, update
            if state_meta['dc_pins']:
                cls.state_ic = state_ic
            if state_meta['periodic']:
                cls.periodic_states = periodic_states
        else:
            xset2 = set(info['pure_spec']['xsyms'])
            ## `sym_spec['G']` IS this Jacobian, computed once in
            ## generate_code; differentiating `ivec` a second time here
            ## bought nothing and is not available on a cache hit.
            Gmat = info['sym_spec']['G']
            cls.linear = (not any((e.free_symbols & xset2) for e in Gmat)
                          ## a wrap or Piecewise is discontinuous even where
                          ## its slope is constant almost everywhere
                          and not any(e.atoms(_wrapfloor) or
                                      e.atoms(sympy.Piecewise)
                                      for e in info['pure_spec']['ivec']))

            cls.update = update
            cls.i, cls.G, cls.q, cls.C, cls.CY = i, G, q, C, CY
            cls.u, cls.dudt = u, dudt
            if state_meta['dc_pins']:
                cls.state_ic = state_ic
            if state_meta['periodic']:
                cls.periodic_states = periodic_states

        ## @cross -- predict where each declared expression crosses zero
        ## and publish it as a breakpoint.  Two accepted points give the
        ## rate; that is all a first-order prediction needs, and a
        ## prediction only has to BRACKET the corner.
        cspec = info['cross_spec']
        if cspec is not None:
            def accept_step(self, t, x, epar=defaultepar, _cs=cspec):
                vals = [float(v) for v in
                        np.atleast_1d(_cs['f'](np.asarray(x, dtype=float),
                                               *_args_of(self, epar)))]
                prev = self.__dict__.get('_hdl_cross')
                self.__dict__['_hdl_cross'] = (
                    (prev[1] if prev else None), (float(t), vals))

            def reset_state(self, epar=None):
                self.__dict__.pop('_hdl_cross', None)

            def next_event(self, t, _cs=cspec):
                hist = self.__dict__.get('_hdl_cross')
                if not hist or hist[0] is None:
                    return self.toolkit.inf
                (t0, v0), (t1, v1) = hist
                if t1 <= t0:
                    return self.toolkit.inf
                best = self.toolkit.inf
                for k, direction in enumerate(_cs['directions']):
                    rate = (v1[k] - v0[k]) / (t1 - t0)
                    if rate == 0.0:
                        continue
                    if direction and np.sign(rate) != direction:
                        continue
                    ## Time from the LAST accepted point to value zero.
                    dt = -v1[k] / rate
                    if dt <= 0.0:
                        continue
                    tc = t1 + dt
                    if tc > t and tc < best:
                        best = tc
                return best

            cls.accept_step = accept_step
            cls.reset_state = reset_state
            cls.next_event = next_event

        ## $limit -- the state-free convention: return a LIMITED COPY of
        ## the sub-vector.  The limiter moves the solution, so nothing is
        ## ever evaluated at a point other than the one it linearises
        ## about, and no private iterate memory is needed (x0 is handed
        ## in).  Only generated when the model asked for it.
        lspec = info['limit_spec']
        if lspec:
            from pycircuit.circuit._limiting import (apply_limit as _lim,
                                                     device_writeback as _dwb)
            lgroups = info.get('limit_groups') or []
            _grouped = set()
            for _s, _ix in lgroups:
                _grouped.update(_ix)
            lsingle = [i for i in range(len(lspec)) if i not in _grouped]

            def limit(self, x, x0, epar=defaultepar, _ls=lspec,
                      _lg=lgroups, _l1=lsingle):
                out = np.array(x, dtype=float, copy=True)
                x0a = np.asarray(x0, dtype=float)
                args = _args_of(self, epar)

                def _pv(pfs_):
                    ## A parameter that reads the solution is evaluated
                    ## at the LAST ACCEPTED iterate -- SPICE's `von`
                    ## semantics -- never at the unlimited `x`, which is
                    ## the point the limiter exists to distrust.
                    return [float(f(x0a, *args)
                                  if getattr(f, '_wants_x', False)
                                  else f(*args)) for f in pfs_]
                ## WHICH TERMINAL ABSORBS THE CORRECTION IS A RUNTIME
                ## QUESTION, and deciding it at compile time made the
                ## limiter a divergence GENERATOR.
                ##
                ## The write-back is `out[ra] = out[rb] + vlim`, which
                ## copies rb's ABSOLUTE value across the branch.  `vlim` is
                ## bounded; `out[rb]` is not.  Measured on a stacked pair,
                ## with the drain held near its rail by a 100k load:
                ##
                ##     it0  nm  +5.17e+07 -> +3.20e+01   (the lower device
                ##                                        limits its own
                ##                                        wild node)
                ##     it1  nd  +4.00e+01 -> +6.14e+08   (the upper device
                ##     it2  nd  +4.00e+01 -> +6.14e+09    then DESTROYS a
                ##     it3  nd  +4.00e+01 -> +5.78e+10    perfectly good
                ##                                        one)
                ##
                ## Newton had the drain at a sane 40 V and the limiter
                ## dragged it out by a decade per iteration, because the
                ## source node was wild and elements are limited in
                ## instance order -- the upper device reads a node the
                ## lower one has not fixed yet.  225 Jacobian evaluations
                ## against 25 unlimited: limiting made it NINE TIMES worse.
                ## The same mechanism moves a gate that a source pins.
                ##
                ## So move the terminal that has drifted FURTHER from the
                ## last accepted point.  That is the node Newton is being
                ## wild about; the other one -- a pinned gate, a
                ## rail-clamped drain -- is information worth keeping
                ## rather than overwriting.  The branch still ends at
                ## exactly `vlim` either way; only which node pays changes.
                ##
                ## ORDER INDEPENDENCE IS PRESERVED, and it is not free:
                ## two probes may not move the same terminal (that is what
                ## the compile-time `move` was for), so the choice has to
                ## be resolved against the OTHER probes, and doing that in
                ## spec-list order would make the result depend on it.
                ## The probes are therefore visited in a CANONICAL order
                ## derived from the data -- widest drift first, ties broken
                ## by row index -- so reversing the spec list changes
                ## nothing.  `test_the_vgs_vds_star_is_order_independent`
                ## asserts exactly that, by reversing it.
                drift = np.abs(np.asarray(out, dtype=float) - x0a)
                _moved = set()

                def _vlim_of(k, i0, i1, pfs_):
                    vn = float(out[i0] - out[i1])
                    vo = float(x0a[i0] - x0a[i1])
                    pv_ = _pv(pfs_)
                    return vn, _lim(k, vn, vo, pv_, self.toolkit)

                ## DEVICE-LEVEL groups first (`limit_together`, roadmap
                ## 10.3(b)).  They see a pristine vector -- a group is a
                ## statement about the whole device and the per-probe
                ## probes below already know how to keep off a row another
                ## probe wrote, which is not true the other way round.
                for _seq, _idx in _lg:
                    ## EACH PROBE READS WHAT THE EARLIER ONES DID, and it
                    ## is the whole difference between a group that helps
                    ## and one that hurts.  Limiting every probe from the
                    ## SAME unlimited vector and then enforcing all of the
                    ## results looks like the natural device-level thing
                    ## and OVER-CORRECTS: probes share terminals, so one
                    ## write often satisfies the next probe outright.
                    ## Measured on the cascode of test_device_limiter, at
                    ## (5 V, 3 V, 0.8 V): `vgs = 57.6` and `vds = 59.6`
                    ## off a source-pinned gate and drain, with only the
                    ## middle node wild.  Read together and enforced, the
                    ## two clamps disagree about the shared node by 2.5 V,
                    ## and the only way to honour both is to move the
                    ## GATE -- the very failure 12.1 removed.  Read in
                    ## sequence, the `vds` clamp alone brings the middle
                    ## node back to 1.0 V and `vgs` no longer bites at
                    ## all, so nothing but the wild node moves.
                    ##
                    ## The order is CANONICAL (largest correction first,
                    ## ties by row index) exactly as the per-probe path's
                    ## is, so reversing the declaration changes nothing.
                    ## `sequential=True` replaces it with SPICE's, which
                    ## is the declaration order and is meant to be.
                    targets, shift, taken = [], {}, set()
                    if _seq:
                        seq_order = list(_idx)
                    else:
                        _rk = {j: _vlim_of(_ls[j][1], _ls[j][0][0],
                                           _ls[j][0][1], _ls[j][3])
                               for j in _idx}
                        seq_order = sorted(
                            _idx, key=lambda j: (-abs(_rk[j][1] - _rk[j][0]),
                                                 _ls[j][0][0], _ls[j][0][1]))
                    for j in seq_order:
                        (ra, rb), kind, _mv, pfs = _ls[j]
                        pv = _pv(pfs)
                        vorig = float(out[ra] - out[rb])
                        vold = float(x0a[ra] - x0a[rb])
                        vin = vorig + shift.get(ra, 0.0) - shift.get(rb, 0.0)
                        vlim = _lim(kind, vin, vold, pv, self.toolkit)
                        if vlim != vin:
                            if _seq:
                                ## SPICE's coupling: `mos1load.c` limits
                                ## `vgs`, then recomputes `vds` from the
                                ## UNLIMITED `vgd` -- which is holding the
                                ## gate and the drain and moving the
                                ## SOURCE, i.e. always the minus terminal.
                                n = rb
                            else:
                                n = ra if drift[ra] >= drift[rb] else rb
                                if n in taken:
                                    n = rb if n == ra else ra
                            taken.add(n)
                            shift[n] = (shift.get(n, 0.0)
                                        + (vlim - vin) * (1 if n == ra
                                                          else -1))
                        targets.append((ra, rb, vorig, vlim))
                    ## The shifts above are a READING order, not the
                    ## write-back: which node finally pays is decided over
                    ## the whole device, from the targets, below.
                    _moved |= _dwb(out, targets, drift, _moved)

                ## WHO GETS THE SHARED TERMINAL when two probes want it:
                ## the one applying the LARGER correction.  A BJT's two
                ## junctions both hang off the base and both see the same
                ## drift when only the base has moved, so drift alone ties
                ## and the tie-break decides real behaviour -- ordering by
                ## `(ra, rb)` handed the base to the base-COLLECTOR probe
                ## and left the base-emitter one writing the emitter,
                ## which is backwards for a forward-biased device.
                ##
                ## The correction size is data-derived, so the order still
                ## does not depend on the spec list's own order.  These
                ## `vlim`s are for RANKING only; each probe recomputes its
                ## own below against whatever earlier probes wrote.
                _rank = {i: _vlim_of(_ls[i][1], _ls[i][0][0], _ls[i][0][1],
                                     _ls[i][3])
                         for i in _l1}
                order = sorted(
                    _l1,
                    key=lambda i: (-abs(_rank[i][1] - _rank[i][0]),
                                   _ls[i][0][0], _ls[i][0][1]))
                for i in order:
                    (ra, rb), kind, move, pfs = _ls[i]
                    vnew = float(out[ra] - out[rb])
                    vold = float(x0a[ra] - x0a[rb])
                    pv = _pv(pfs)
                    vlim = _lim(kind, vnew, vold, pv, self.toolkit)
                    ## A LIMITER THAT DID NOT BITE MUST TOUCH NOTHING.
                    ## Writing `out[rb] = out[ra] - vlim` when `vlim` is
                    ## already `vnew` does not round-trip -- `a - (a - b)`
                    ## is not `b` in floating point -- so the write alone
                    ## broke "a step of at most 2*VT passes through
                    ## exactly", which is what lets "did limiting fire?"
                    ## be a convergence signal.  Skipping also leaves the
                    ## terminal free for the other probe, which is what
                    ## it would have wanted anyway.
                    if vlim == vnew:
                        continue
                    cand = ra if drift[ra] >= drift[rb] else rb
                    if cand in _moved:
                        cand = rb if cand == ra else ra
                    if cand in _moved:
                        cand = move          # both spoken for
                    _moved.add(cand)
                    if cand == ra:
                        out[ra] = out[rb] + vlim
                    else:
                        out[rb] = out[ra] - vlim
                return out

            cls.limit = limit

        ## PCNR: hand the layer the three things it needs, plus the
        ## limiter for this device's own quantity (pcnr.refine's hook).
        ## Generated only for a recognisably exponential single-branch
        ## element -- see generate_code.  With PCNR enabled the device's
        ## limiting is then PCNR's job, which is the point of the method.
        pf = info['pcnr_funcs']
        if pf is not None:
            from pycircuit.circuit._limiting import _pnjlim
            _pn_pcnr = info['paramnames']

            cls.pcnr_junctions = tuple(tuple(f['terminals']) for f in pf)

            def _pcnr_compiled(jn, which, toolkit, _pf=pf,
                               _pn=_pn_pcnr):
                """The numpy form, or its jax twin for a traced call."""
                if not getattr(toolkit, 'jax', False):
                    return _pf[jn]['i' if which == 'i' else 'didv']
                cache = _pf[jn].setdefault('_jax', {})
                if which not in cache:
                    import jax.numpy as _jnp
                    sym = _pf[jn]['sym']
                    if sym is None:
                        ## The chained twin: the same let-chain printed
                        ## into a jax namespace.  `_ChainPrinter` emits
                        ## `numpy.<f>` calls, so binding the NAME `numpy`
                        ## to `jax.numpy` is the whole swap.
                        ch = _pf[jn]['chain']
                        defs_ = ch['defs'] if which == 'i' else \
                            ch['defs'] + ch['ddefs']
                        out_ = ch['out'] if which == 'i' else ch['dout']
                        f_ = _chain_compile(
                            defs_, [out_],
                            [ch['vsym']]
                            + [_param_symbol(q) for q in _pn] + [TEMP],
                            modules_map=dict(_kernel_jax(_jnp),
                                             numpy=_jnp))
                        cache[which] = lambda *a_, _f=f_: _f(*a_)[0]
                    else:
                        expr = sym['f'] if which == 'i' else sym['dfdv']
                        cache[which] = sympy.lambdify(
                            [sym['vsym']]
                            + [_param_symbol(q) for q in _pn] + [TEMP],
                            expr, modules=[_kernel_jax(_jnp), 'jax'],
                            cse=True)
                return cache[which]

            def pcnr_i(v, params, epar, toolkit, jn=0, _pn=_pn_pcnr):
                f = _pcnr_compiled(jn, 'i', toolkit)
                cur = f(v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([cur, -cur])

            def pcnr_didv(v, params, epar, toolkit, jn=0, _pn=_pn_pcnr):
                f = _pcnr_compiled(jn, 'didv', toolkit)
                g = f(v, *[params[q] for q in _pn], _epar_T(epar))
                return toolkit.array([g, -g])

            def pcnr_scales(params, epar, jn=0, _pf=pf, _pn=_pn_pcnr):
                """(VT, IS) for this junction -- what pnjlim needs, and
                what the traced backend cannot read off an expression."""
                args = [params[q] for q in _pn] + [_epar_T(epar)]
                return (float(_pf[jn]['VT'](*args)),
                        float(_pf[jn]['IS'](*args)))

            cls.pcnr_scales = staticmethod(pcnr_scales)

            def pcnr_limit(vnew, vold, params, epar, toolkit, jn=0,
                           _pf=pf, _pn=_pn_pcnr):
                args = [params[q] for q in _pn] + [_epar_T(epar)]
                return _pnjlim(vnew, vold, float(_pf[jn]['VT'](*args)),
                               float(_pf[jn]['IS'](*args)), toolkit)

            cls.pcnr_i = staticmethod(pcnr_i)
            cls.pcnr_didv = staticmethod(pcnr_didv)
            cls.pcnr_limit = staticmethod(pcnr_limit)

        ## VECTOR PCNR, the declared-probe route: the device's `$limit`
        ## probes are its block of limited unknowns.  `pcnr_probes` (not
        ## `pcnr_junctions`) is the declaration `pcnr.PcnrDevice` reads for
        ## this protocol; the three callables take and return whole blocks.
        pv = info.get('pcnr_vector')
        if pv is not None:
            from pycircuit.circuit._limiting import apply_limit as _apply_lim
            from pycircuit.circuit import pcnr as _pcnr_mod
            _pn_vec = info['paramnames']
            _gn_vec = info['given_names']
            _lpars = pv['limit_pars']
            _lkinds = pv['kinds']
            cls.pcnr_probes = tuple((tuple(p), k)
                                    for p, k in zip(pv['probes'], pv['kinds']))
            ## The redundant probes (a cycle in the declaration): local
            ## rows, kind, and the signed combination of the tree unknowns
            ## that IS their branch voltage.  `pcnr.PcnrDevice` lists them
            ## as gmin targets when they are junctions; `pcnr_limit` below
            ## honours their laws over the tree.
            cls.pcnr_redundant = tuple((tuple(p), k, tuple(c))
                                       for _j, p, k, c in pv['redundant'])
            ## `(kind, parameter callables)` for every probe whose law is
            ## applied, TREE probes first (in unknown order), then the
            ## redundant ones -- the order `limit_block` takes them in.
            _laws = [(k, _lpars[j]) for j, k in zip(pv['spec_idx'], _lkinds)] \
                + [(k, _lpars[j]) for j, _p, k, _c in pv['redundant']]
            _all_probes = list(pv['probes']) + [p for _j, p, _k, _c
                                                in pv['redundant']]
            _m_tree = len(pv['probes'])
            _coeffs = [[1 if q == j else 0 for q in range(_m_tree)]
                       for j in range(_m_tree)] \
                + [list(c) for _j, _p, _k, c in pv['redundant']]

            def _pcnr_vec_compiled(which, toolkit, _pv=pv, _pn=_pn_vec):
                """The numpy form, or its jax twin for a traced call."""
                if not getattr(toolkit, 'jax', False):
                    return _pv[which]
                cache = _pv.setdefault('_jax', {})
                if which not in cache:
                    import jax.numpy as _jnp
                    args_ = list(_pv['vsyms']) \
                        + [_param_symbol(q) for q in _pn] + [TEMP] \
                        + [sympy.Symbol('_hdl_given_%s' % q, real=True)
                           for q in _gn_vec]
                    if _pv['sym'] is None:
                        ch = _pv['chain']
                        cache[which] = _chain_compile(
                            ch['defs'], ch['out'], args_,
                            want_jacobian_of=(True if which == 'didv'
                                              else None),
                            xsyms=(_pv['vsyms'] if which == 'didv'
                                   else None),
                            modules_map=dict(_kernel_jax(_jnp), numpy=_jnp))
                    else:
                        cache[which] = sympy.lambdify(
                            args_, _pv['sym'][which],
                            modules=[_kernel_jax(_jnp), 'jax'], cse=True)
                return cache[which]

            def _vec_args(params, epar, _pn=_pn_vec, _gn=_gn_vec):
                return [params[q] for q in _pn] + [_epar_T(epar)] \
                    + [float(params.get('$given:' + nm, 0.0)) for nm in _gn]

            def pcnr_i(v, params, epar, toolkit):
                f = _pcnr_vec_compiled('i', toolkit)
                return toolkit.array(f(*list(v), *_vec_args(params, epar)))

            def pcnr_didv(v, params, epar, toolkit):
                f = _pcnr_vec_compiled('didv', toolkit)
                return toolkit.array(f(*list(v), *_vec_args(params, epar)))

            def pcnr_limit(v_new, v_old, params, epar, toolkit, x_old_sub,
                           _laws=_laws, _probes=_all_probes,
                           _coeffs=_coeffs):
                """Each probe's own declared law, over the block.  A
                parameter that reads the solution (`_wants_x`) is
                evaluated at the device's LAST ACCEPTED sub-vector --
                the generated `limit()`'s semantics, SPICE's `von`.
                `pcnr.limit_block` applies the laws: directly when the
                probes are a tree, by `limit_together`'s rule (largest
                corrections win, the smallest on a cycle yields) when a
                redundant probe's law is in play too."""
                args = _vec_args(params, epar)
                fns = []
                for kind, pfs in _laws:
                    pars = [float(f(x_old_sub, *args)
                                  if getattr(f, '_wants_x', False)
                                  else f(*args)) for f in pfs]
                    fns.append(lambda vn, vo, _k=kind, _p=pars:
                               _apply_lim(_k, vn, vo, _p, toolkit))
                ## Looked up at call time so that a measurement can swap
                ## the rule (`test_device_limiter.py` measures the redundant
                ## law against its absence).
                return _pcnr_mod.limit_block(_probes, _coeffs, v_new, v_old,
                                             fns)

            cls.pcnr_i = staticmethod(pcnr_i)
            cls.pcnr_didv = staticmethod(pcnr_didv)
            cls.pcnr_limit = staticmethod(pcnr_limit)

        ## eval_i_pure / eval_q_pure: compiled on first use with sympy's
        ## jax printer; the staticmethods exist only if that succeeds, so
        ## the class is admitted to the vmap groups exactly when the pure
        ## forms are real.  Params arrive as a dict on the batched path.
        cls._hdl_pure_cache = {}
        spec = info['pure_spec']
        if spec is None:
            return
        pnames = info['paramnames']

        def _compiled_pure(which):
            cache = cls._hdl_pure_cache
            if which not in cache:
                import jax.numpy as jnp
                xv = sympy.DeferredVector('x')
                xs = dict(zip(spec['xsyms'],
                              [xv[i2] for i2 in range(spec['n'])]))
                vec = [e.subs(xs) for e in spec[which]]
                cache[which] = sympy.lambdify(
                    [xv] + [_param_symbol(p2) for p2 in pnames] + [TEMP],
                    vec, modules=[_kernel_jax(jnp), 'jax'],
                    cse=True)
            return cache[which]

        def eval_i_pure(x, params, epar, toolkit):
            f = _compiled_pure('ivec')
            return toolkit.array(f(x, *[params[p2] for p2 in pnames],
                                   _epar_T(epar)))

        def eval_q_pure(x, params, epar, toolkit):
            f = _compiled_pure('qvec')
            return toolkit.array(f(x, *[params[p2] for p2 in pnames],
                                   _epar_T(epar)))

        try:
            import sympy.printing.numpy as _sp_np
            have_jax = hasattr(_sp_np, 'JaxPrinter')
        except Exception:
            have_jax = False
        if have_jax and not state_meta['dc_pins']:
            ## DC-pinned elements are kept OFF the batched path for now:
            ## the pure form has no epar-driven variant selection, and a
            ## silently-unpinned batched DC is the failure mode item 2 of
            ## the roadmap just removed.
            cls.eval_i_pure = staticmethod(eval_i_pure)
            cls.eval_q_pure = staticmethod(eval_q_pure)


class Behavioural(circuit.Circuit, metaclass=BehaviouralMeta):
    """Behavioural circuit element -- equations in, element out.

    Write the element's terminal behaviour as contribution statements in
    an ``analog()`` staticmethod whose argument names declare the
    terminals; the metaclass compiles them into the full pycircuit element
    interface (``i``/``q``/``G``/``C``/``u``/``dudt``, DC pinning,
    ``uic`` seeding, the JAX pure forms, gauge-shift declarations).

    Example::

        class MyConductor(Behavioural):
            instparams = [Parameter(name='g', desc='Conductance', unit='S',
                                    default=1e-3)]
            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, g * b.V)

    See ``hdl.md`` for the full language: V-contributions, ``ddt``,
    ``idt``/``idtmod``, internal nodes, current probes, and time-dependent
    sources via the ``TIME`` symbol.
    """

    #: Verilog-A's ``aliasparam``: alternative accepted spellings,
    #: ``{alias: canonical}``.  Compact models carry these by the dozen
    #: (217 declarations in the vacask device library) because one device
    #: is known by different parameter names in different foundries'
    #: decks.  Declared per class::
    #:
    #:     aliasparams = {'isat': 'IS', 'js': 'IS'}
    aliasparams = {}

    #: The parameter namespace: ``params_as = 'p'`` makes analog()'s
    #: FIRST argument a `ParamNamespace` (``p.bf``, ``p['lambda']``,
    #: ``p.given('rb')``) instead of a terminal, and stops binding the
    #: parameters as bare names.  `None` is the bare-name style.  The
    #: compiled element is identical either way; see `ParamNamespace`.
    params_as = None

    def __init__(self, *args, **kvargs):
        aliases = type(self).aliasparams
        if aliases:
            for alias, canonical in aliases.items():
                if alias in kvargs:
                    if canonical in kvargs:
                        raise ValueError(
                            '%r was given both as %r and as its alias %r'
                            % (canonical, canonical, alias))
                    kvargs[canonical] = kvargs.pop(alias)
        ## CONNECTION COUNT.  `Circuit.__init__` builds its terminal hook
        ## with `zip(self.terminals, args)`, and zip is silent about a
        ## length mismatch: a three-terminal element given two nodes was
        ## accepted and simply left its third pin unconnected, and a
        ## fourth node was dropped.  Nothing downstream can tell that
        ## from a deliberate name-based connection, which is why zero
        ## arguments -- the `add_instance(..., plus=n1)` form -- is still
        ## allowed.
        terminals = list(type(self).terminals)
        if args and len(args) != len(terminals):
            raise TypeError(
                '%s has %d terminals (%s) but was given %d connection '
                'node%s. An HDL element\'s pins are analog()\'s argument '
                'names, in that order; a mismatched list used to be zipped '
                'silently, leaving pins unconnected or dropping nodes. '
                'Connect all of them positionally, or none of them here and '
                'all by name via add_instance(name, instance, %s=...).'
                % (type(self).__name__, len(terminals),
                   ', '.join(terminals), len(args),
                   '' if len(args) == 1 else 's', terminals[0]))
        try:
            super().__init__(*args, **kvargs)
        except KeyError as e:
            ## `ParameterDict` raises a bare `KeyError: 'parameter R not in
            ## parameter dictionary'` -- no class, no list of what would
            ## have been accepted, and it is reached by every typo and
            ## every case-wrong name.
            detail = e.args[0] if e.args else ''
            if not (isinstance(detail, str)
                    and 'not in parameter dictionary' in detail):
                raise
            unknown = detail.split()[1]
            names = [pp.name for pp in self.instparams]
            hint = ''
            low = {nm.lower(): nm for nm in names}
            if unknown.lower() in low:
                hint = (' Parameter names are case sensitive: did you mean '
                        '%r?' % low[unknown.lower()])
            elif unknown in terminals:
                hint = (' %r is a TERMINAL of this element, not a parameter; '
                        'terminals are connected positionally.' % unknown)
            elif aliases and unknown in aliases:
                hint = ' (%r is an alias for %r.)' % (unknown,
                                                      aliases[unknown])
            raise TypeError(
                '%s has no parameter %r.%s It accepts: %s.%s'
                % (type(self).__name__, unknown, hint,
                   ', '.join(names) or '(none)',
                   '' if not aliases else ' Aliases: %s.'
                   % ', '.join('%s -> %s' % kv
                               for kv in sorted(aliases.items())))) from None
        info = getattr(self, '_hdl_info', None)
        if info is None:
            return
        ## Compilation is keyed on `analog` appearing in the class body, so
        ## a SUBCLASS that changes only `instparams` would inherit the
        ## parent's compiled functions -- which read a different parameter
        ## list -- and answer confidently with the wrong numbers.  Caught
        ## here rather than left to produce plausible waveforms.
        declared = [p.name for p in self.instparams]
        if declared != info['paramnames']:
            raise TypeError(
                '%s changed instparams (%r) without redefining analog(), so '
                'it would run its base class\'s compiled code for %r. '
                'Redeclare analog() in this class (a one-line '
                '`analog = staticmethod(Base.analog)` is enough to trigger '
                'recompilation).'
                % (type(self).__name__, declared, info['paramnames']))
        ## PARAMETER-DRIVEN COLLAPSE: pick the variant this instance's
        ## parameters select, before any node or branch is registered.
        ## `super().__init__` has already copied the class's branch list,
        ## so retargeting here means replacing it too.
        if info['collapse_conds']:
            mask = _collapse_mask_of(type(self), _params_of(self))
            if mask != info['collapse_mask']:
                var = _collapse_variant(type(self), mask)
                self.__class__ = var
                info = var._hdl_info
                self.branches = list(var.branches)
            self._hdl_collapse_seen = mask
        for name in info['internalnames'] + info['state_meta']['statenames']:
            self.add_node(name)

    def next_event(self, t):
        return self.toolkit.inf

    def explain(self, **kvargs):
        """This instance's compilation record as text; see
        :func:`explain`."""
        return explain(self, **kvargs)

    def check_jacobians(self, x, epar=None, **kvargs):
        """Finite-difference this instance's G and C; see
        :func:`check_jacobians`."""
        return check_jacobians(self, x, epar=epar, **kvargs)


## ----------------------------------------------------------------------
## Instruments.
##
## Two public entry points for things that were private knowledge: what
## the compiler DID with a model (`explain`), and whether the Jacobians
## it derived still agree with the vectors they were derived from
## (`check_jacobians`).  Both existed as ad-hoc code in tests and
## notebooks -- an `x`-vector layout reverse-engineered from `el.nodes`
## and `el.branches` (and got wrong: `IndexError: index 5 is out of
## bounds for axis 0 with size 5`), and a 25-line hand-rolled finite
## difference per model.
## ----------------------------------------------------------------------

def _hdl_info_of(target):
    """`(class, _hdl_info)` for a compiled element class or instance."""
    cls = target if isinstance(target, type) else type(target)
    info = getattr(cls, '_hdl_info', None)
    if info is None:
        raise TypeError(
            '%s is not a compiled HDL element: it has no _hdl_info, so it '
            'was not built from an analog() body. explain() and '
            'check_jacobians() read the compiler\'s own record; a '
            'hand-written element has none.'
            % getattr(cls, '__name__', repr(target)))
    return cls, info


def x_layout(target):
    """The element's unknown vector, entry by entry.

    Returns ``[(index, name, kind)]`` in ``Circuit.update_node_map``
    order -- terminals, internal nodes, ``idt``/``idtmod``/laplace
    states, branch currents -- which is the order of the ``x`` that
    ``i``, ``q``, ``G`` and ``C`` take and of the rows and columns they
    return.

    Worth having as a function because the layout is *derivable* from
    ``el.nodes`` and ``el.branches`` and deriving it is where people get
    it wrong: the state nodes sit between the internal nodes and the
    branch currents, and they are not in the model's source at all.
    """
    _cls, info = _hdl_info_of(target)
    sm = info['state_meta']
    out = []
    for nm in info['terminalnames']:
        out.append((len(out), nm, 'terminal'))
    for nm in info['internalnames']:
        out.append((len(out), nm, 'internal node'))
    periodic = {k for k, _m, _o in sm['periodic']}
    for nm in sm['statenames']:
        k = len(out)
        kind = 'state'
        if k in periodic:
            kind += ', bounded (idtmod)'
        if k in sm['dc_pins']:
            kind += ', pinned at DC by its ic'
        out.append((k, nm, kind))
    for k, key in enumerate(info['branchpairs']):
        out.append((len(out), '_i_br%d' % k,
                    'branch current through %s' % _fmt_branch(key)))
    return out


def _generated_source(f):
    """The text `lambdify` or the let-chain compiler produced, or None."""
    src = getattr(f, '_src', None)
    if src is not None:
        return src
    try:
        return inspect.getsource(f)
    except (OSError, TypeError):
        return None


def _clip(text, maxlines, what):
    lines = text.splitlines()
    if maxlines is None or len(lines) <= maxlines:
        return lines
    return lines[:maxlines] + ['... %d more lines of %s (pass maxlines=None '
                               'for all of it)' % (len(lines) - maxlines,
                                                   what)]


def explain(target, source=True, symbolic=True, maxlines=40):
    """Everything the compiler knows about an element, as text.

    Accepts a `Behavioural` subclass or an instance of one; with an
    instance the parameter values the generated code will actually read
    are shown too, so a parameter expression that nothing has resolved is
    visible as a value rather than as a wrong answer later.

    Reports, in order: the terminals, the full ``x``-vector layout (see
    :func:`x_layout`), the parameters, which compilation path the model
    took (eager `lambdify` or the ``var()`` let-chain), what the compiler
    found in it (states, collapses, crossings, ``$limit``, PCNR
    junctions, AC stimulus, constant stamps, JAX pure forms), the
    symbolic ``i``/``q``/``G``/``C``, and the generated source.

    ``source``/``symbolic`` switch the last two sections off; ``maxlines``
    truncates each of them (``None`` for the whole text) -- a compact
    model's generated ``G`` is tens of thousands of lines.
    """
    cls, info = _hdl_info_of(target)
    inst = None if isinstance(target, type) else target
    sm = info['state_meta']
    lines = []
    layout = x_layout(target)
    ## An instance with collapses runs as a compiled VARIANT subclass,
    ## so `type(el).__name__` is not the name anyone wrote.  Report both.
    base = getattr(cls, '_hdl_collapse_base', None)
    title = cls.__name__ if base is None else \
        '%s (compiled variant %s)' % (base.__name__, cls.__name__)
    lines.append('%s: %d terminal%s, %d unknown%s'
                 % (title, len(info['terminalnames']),
                    '' if len(info['terminalnames']) == 1 else 's',
                    len(layout), '' if len(layout) == 1 else 's'))
    lines.append('')
    lines.append('x vector (the order i/q/G/C use):')
    for k, nm, kind in layout:
        lines.append('  x[%d]  %-12s %s' % (k, nm, kind))

    lines.append('')
    if info['paramnames']:
        lines.append('parameters (%d):' % len(info['paramnames']))
        for nm in info['paramnames']:
            if inst is None:
                lines.append('  %s' % nm)
                continue
            resolved = getattr(inst.iparv, nm, None)
            written = inst.ipar.get(nm) if hasattr(inst.ipar, 'get') else None
            note = ''
            if written is not None and written != resolved:
                ## The value the generated code reads is `iparv`.  A
                ## difference means the parameter is an unresolved
                ## expression, and the element is quietly running on the
                ## DEFAULT -- the most confusing wrong answer this layer
                ## produces, and invisible everywhere else.
                note = ('   <- ipar says %r; update_iparv() has not resolved '
                        'it, so the value on the left is what the generated '
                        'code reads' % (written,))
            lines.append('  %-12s = %r%s%s'
                         % (nm, resolved,
                            '' if nm not in info['given_names'] else
                            '   [param_given]', note))
    else:
        lines.append('parameters: none')

    lines.append('')
    lines.append('compiled as: %s'
                 % ('a let-chain (the model uses var(); i/q/G/C are '
                    'generated Python, differentiated by forward '
                    'accumulation)' if info['chained'] else
                    'flat lambdify expressions'))
    ## Which backend the class actually RUNS -- not which was asked
    ## for.  The status carries the reason whenever they differ.
    lines.append('backend: %s'
                 % getattr(cls, '_hdl_backend_status', 'numpy'))
    feats = []
    if sm['statenames']:
        feats.append('%d state%s' % (len(sm['statenames']),
                                     '' if len(sm['statenames']) == 1 else 's'))
    if sm['has_time_dependence']:
        feats.append('time-dependent source u(t)')
    if info['has_ac']:
        feats.append('AC stimulus')
    if info['collapse_conds']:
        feats.append('%d collapse condition%s (this variant takes %r)'
                     % (len(info['collapse_conds']),
                        '' if len(info['collapse_conds']) == 1 else 's',
                        info['collapse_mask']))
    if info['cross_spec']:
        feats.append('%d @cross event%s'
                     % (len(info['cross_spec']['directions']),
                        '' if len(info['cross_spec']['directions']) == 1
                        else 's'))
    if info['limit_spec']:
        _xl = {k: nm for k, nm, _kind in layout}
        ## Name the KIND: with three of them a bare count no longer says
        ## what the model asked for, and this line is the only place a
        ## reader can see it without recompiling.
        _grp = {}
        for _g, (_seq, _idx) in enumerate(info.get('limit_groups') or []):
            for _i in _idx:
                _grp[_i] = ' [group %d%s]' % (_g + 1,
                                              ', SPICE order' if _seq else '')
        feats.append('%d $limit probe%s (%s)'
                     % (len(info['limit_spec']),
                        '' if len(info['limit_spec']) == 1 else 's',
                        ', '.join('%s on (%s,%s)%s%s'
                                  % (_LIMIT_LAW[k],
                                     _xl[ra], _xl[rb], _grp.get(_i, ''),
                                     ' [params at last iterate]'
                                     if any(getattr(f, '_wants_x', False)
                                            for f in _p) else '')
                                  for _i, ((ra, rb), k, _m, _p)
                                  in enumerate(info['limit_spec']))))
    ## Name WHICH rule refused.  The reason used to exist only inside
    ## generate_code, so "does not qualify" was all a reader ever saw and
    ## the first refusal had to be replayed by hand to find out why
    ## (2026-08-25: PSP is refused by its GATE RESISTOR, not its drain
    ## current, and that took a re-measurement to learn).
    _pv = info.get('pcnr_vector')
    if _pv:
        _xl = {k: nm for k, nm, _kind in layout}
        feats.append('PCNR: vector route, %d probe%s over %d row%s (%s)'
                     % (_pv['m'], '' if _pv['m'] == 1 else 's',
                        _pv['t'], '' if _pv['t'] == 1 else 's',
                        ', '.join('%s on (%s,%s)'
                                  % (_LIMIT_LAW[k], _xl[ra], _xl[rb])
                                  for (ra, rb), k in zip(_pv['probes'],
                                                         _pv['kinds']))))
        for _p, _k, _c in [(p, k, c) for _j, p, k, c in _pv['redundant']]:
            feats.append('PCNR redundant probe: %s on (%s,%s) = %s -- its '
                         'branch voltage is read off the tree unknowns and '
                         'its law is applied over them'
                         % (_LIMIT_LAW[_k], _xl[_p[0]], _xl[_p[1]],
                            ' '.join(('%s(%s,%s)' % ('+' if c_ > 0 else '-',
                                                     _xl[pa], _xl[pb]))
                                     for c_, (pa, pb) in zip(_c, _pv['probes'])
                                     if c_).lstrip('+')))
    else:
        feats.append('PCNR: %s'
                     % ('%d junction(s)' % len(info['pcnr_funcs'])
                        if info['pcnr_funcs'] else
                        'does not qualify -- %s. (Rule: with $limit probes, '
                        'every current reaches the solution only through '
                        'pnjlim probes; without, every current an exponential '
                        'function of its own branch voltage alone; no states, '
                        'no branch unknowns; charge and var() are allowed.)'
                        % (info.get('pcnr_refusal') or 'unknown')))
    feats.append('JAX pure forms: %s'
                 % ('yes' if info['pure_spec'] else
                    'no (the let-chain path has none)'))
    feats.append('constant G: %s, constant C: %s'
                 % (info['const_G'], info['const_C']))
    for f in feats:
        lines.append('  ' + f)

    spec = info['sym_spec']
    if symbolic:
        lines.append('')
        if spec is None:
            lines.append('symbolic i/q/G/C: not kept -- this model compiles '
                         'through the let-chain, where the vectors are '
                         'assignments rather than one expression. Read the '
                         'generated source below instead.')
        else:
            n = spec['shape'][0]
            text = []
            for which in ('i', 'q'):
                for k, e in enumerate(spec[which]):
                    if e != 0:
                        text.append('%s[%d] = %s' % (which, k, e))
            for which in ('G', 'C'):
                for k, e in enumerate(spec[which]):
                    if e != 0:
                        text.append('%s[%d,%d] = %s'
                                    % (which, k // n, k % n, e))
            lines.append('symbolic i/q/G/C (nonzero entries):')
            lines += ['  ' + t for t in
                      _clip('\n'.join(text), maxlines, 'symbolic entries')]

    if source:
        lines.append('')
        src = _generated_source(info['funcs']['i'])
        if src is None:
            lines.append('generated source for i(x): unavailable')
        else:
            lines.append('generated source for i(x):')
            lines += ['  ' + t for t in _clip(src, maxlines,
                                              'generated source')]
    return '\n'.join(lines)


## Per-entry verdicts.  UNRESOLVED sits BETWEEN the two, and it is not a
## softer failure: it says the finite difference at this point carries
## more noise than the discrepancy being reported, so the comparison has
## no power here.  See `check_jacobians` for what puts an entry there.
JAC_PASS, JAC_UNRESOLVED, JAC_FAIL = 0, 1, 2

## Double precision's unit roundoff, used to size the roundoff floor.
_JAC_EPS = float(np.finfo(float).eps)

#: One entry of a :class:`JacobianCheck`, as listed by its
#: :attr:`~JacobianCheck.unresolved` and :attr:`~JacobianCheck.failures`.
#: ``reason`` is ``'roundoff'``, ``'truncation'`` or ``'kink'`` for an
#: unresolved entry and ``''`` for a failing one.
JacEntry = collections.namedtuple(
    'JacEntry', 'which row col ana fd err floor reason')


def _jac_widen(col, hk, frozen, xk, fd, rnd, trunc, rungs=8, clear=32.0):
    """Re-probe a column at a wider step, for entries the step ``hk``
    could not move at all, and write the result back into ``fd``/``rnd``.

    ``col(s)`` evaluates the differentiated vector at ``x + s*hk*e_k``.
    ``frozen`` selects the rows to work on.  Modifies ``fd``, ``rnd`` and
    ``trunc`` in place, for those rows only.

    The ladder climbs by a decade at a time and stops for a row once its
    two-sided difference is ``clear`` times the SMALLEST non-zero
    difference that row has shown -- which is what "the signal is above
    the quantisation" means when the quantisation cannot be computed.
    The roundoff floor is then that quantum divided by the widened step,
    and the truncation floor comes from the two rungs either side, whose
    steps differ by ten (so the ``h^2`` term differs by 99).

    A row that never moves is left exactly as it was: ``fd = 0`` and no
    floor, so a non-zero analytic entry against it still FAILS.
    """
    ## Widening past a fraction of the variable's own size stops being a
    ## measurement of the derivative here, so the ladder is capped.
    cap = 1e-2 * max(1.0, abs(xk))
    todo = np.array(frozen, dtype=bool)
    quantum = np.zeros_like(fd)
    prev = np.zeros_like(fd)
    have_prev = np.zeros_like(fd, dtype=bool)
    for m in range(1, rungs + 1):
        H = hk * 10.0 ** m
        if H > cap:
            break
        d = col(10.0 ** m) - col(-10.0 ** m)
        g = d / (2 * H)
        moved = todo & (d != 0.0) & np.isfinite(d)
        fresh = moved & (quantum == 0.0)
        quantum[fresh] = np.abs(d[fresh])
        done = moved & (np.abs(d) >= clear * quantum)
        if np.any(done):
            fd[done] = g[done]
            rnd[done] = quantum[done] / (2 * H)
            trunc[done] = np.where(have_prev[done],
                                   2.0 * np.abs(prev[done] - g[done]) / 99.0,
                                   0.0)
            todo &= ~done
        ## Whatever the last rung reached is the fallback for a row that
        ## runs out of ladder: it is still better than a difference of
        ## exactly zero, and its floor says how little it is worth.
        keep = todo & moved
        fd[keep] = g[keep]
        rnd[keep] = quantum[keep] / (2 * H)
        prev[moved] = g[moved]
        have_prev[moved] = True
        if not np.any(todo):
            break


class JacobianCheck(object):
    """The result of :func:`check_jacobians`; print it."""

    def __init__(self, name, x, layout, results, nonfinite, rtol, atol):
        #: Element class name.
        self.name = name
        #: The point the check was made at.
        self.x = x
        #: `x_layout` of the element, used to name rows and columns.
        self.layout = layout
        #: ``{'G': dict(...), 'C': dict(...)}`` -- per pair, the analytic
        #: matrix, the finite difference, the worst entry and its error.
        self.results = results
        #: ``[(what, index)]`` for every non-finite entry found.
        self.nonfinite = nonfinite
        self.rtol, self.atol = rtol, atol

    @property
    def ok(self):
        """No FAILING entry and nothing non-finite.

        An UNRESOLVED entry does not make this false -- see
        :attr:`unresolved`, which is where a model author has to look to
        find out what the check did *not* answer.
        """
        return not self.nonfinite and all(r['ok'] for r in
                                          self.results.values())

    @property
    def resolved(self):
        """True when every entry got a real verdict.

        ``ok and resolved`` is the strong statement; ``ok and not
        resolved`` means "nothing was caught, and here is what could not
        have been".
        """
        return not any(r['nunresolved'] for r in self.results.values())

    def _entries(self, want):
        out = []
        for which, r in self.results.items():
            st = r['status']
            for i, k in zip(*np.nonzero(st == want)):
                i, k = int(i), int(k)
                out.append(JacEntry(which, i, k, float(r['ana'][i, k]),
                                    float(r['fd'][i, k]),
                                    float(r['err_mat'][i, k]),
                                    float(r['floor'][i, k]),
                                    r['reason'][i][k] if want ==
                                    JAC_UNRESOLVED else ''))
        return sorted(out, key=lambda e: -e.err)

    @property
    def unresolved(self):
        """``[JacEntry]``, worst first, for entries the finite difference
        could not resolve."""
        return self._entries(JAC_UNRESOLVED)

    @property
    def failures(self):
        """``[JacEntry]``, worst first, for entries that genuinely
        disagree."""
        return self._entries(JAC_FAIL)

    @property
    def verdict(self):
        """``'ok'``, ``'UNRESOLVED'``, ``'FAILED'`` or ``'NOT
        COMPARABLE'`` -- the one-word summary the ``__repr__`` prints."""
        if self.nonfinite:
            return 'NOT COMPARABLE'
        if not self.ok:
            return 'FAILED'
        return 'ok' if self.resolved else 'UNRESOLVED'

    def _label(self, k):
        return self.layout[k][1] if k < len(self.layout) else '?'

    def __repr__(self):
        head = self.verdict
        nu = sum(r['nunresolved'] for r in self.results.values())
        if nu and head in ('ok', 'UNRESOLVED', 'FAILED'):
            head += ' (%d entr%s not resolvable here)' % (
                nu, 'y' if nu == 1 else 'ies')
        out = ['JacobianCheck(%s): %s' % (self.name, head)]
        out.append('  x = %s' % np.array2string(self.x, precision=6))
        for which, r in self.results.items():
            against = 'i' if which == 'G' else 'q'
            if not r['finite']:
                out.append('  %s vs %s: NOT COMPARABLE -- non-finite entries '
                           '(listed below); a finite-difference of a '
                           'non-finite function says nothing'
                           % (which, against))
                continue
            if r['scale'] == 0.0:
                out.append('  %s vs %s: both identically zero' % (which,
                                                                  against))
                continue
            k = r['worst_resolved']
            out.append(
                '  %s vs %s: worst |%s - d%s/dx| = %.3e at [%d,%d] '
                '(%s / %s), %s = %.6e, fd = %.6e  %s'
                % (which, against, which, against, r['err_mat'][k], k[0], k[1],
                   self._label(k[0]), self._label(k[1]), which,
                   r['ana'][k], r['fd'][k], 'ok' if r['ok'] else
                   'FAILS rtol=%g atol=%.3e' % (self.rtol, r['atol'])))
            if r['nunresolved']:
                e = [u for u in self.unresolved if u.which == which][0]
                out.append(
                    '  %s vs %s: %d entr%s UNRESOLVED -- the finite '
                    'difference is noisier than the discrepancy, so the '
                    'comparison says nothing there.' % (
                        which, against, r['nunresolved'],
                        'y is' if r['nunresolved'] == 1 else 'ies are'))
                out.append(
                    '      worst at [%d,%d] (%s / %s): %s = %.6e, fd = %.6e, '
                    '|diff| = %.3e, noise floor %s (%s)'
                    % (e.row, e.col, self._label(e.row), self._label(e.col),
                       which, e.ana, e.fd, e.err,
                       'UNBOUNDED -- the difference is not in its '
                       'asymptotic regime at this step' if
                       not np.isfinite(e.floor) else '%.3e' % e.floor,
                       e.reason))
        if self.nonfinite:
            for what, idx in self.nonfinite[:12]:
                out.append('  NON-FINITE: %s[%s]  (%s)' % (
                    what, ','.join(str(i) for i in idx),
                    ', '.join(self._label(i) for i in idx)))
            if len(self.nonfinite) > 12:
                out.append('  ... and %d more non-finite entries'
                           % (len(self.nonfinite) - 12))
        else:
            out.append('  i, q, G, C: all finite')
        return '\n'.join(out)


def check_jacobians(element, x, epar=None, h=None, rtol=1e-5, atol=None):
    """Differentiate ``i`` and ``q`` numerically and compare against the
    ``G`` and ``C`` the compiler derived, at the point ``x``.

    This is the instrument the NaN-Jacobian warnings in the
    documentation never handed anyone.  It answers two different
    questions at once, and they fail in different ways:

    * **is the derivative right?**  Central differences of ``i`` against
      ``G``, of ``q`` against ``C``.  For a generated element this
      cannot fail by transcription -- both come from one expression --
      but it does fail when a kernel primitive's ``fdiff`` is wrong, and
      that is invisible in the value.
    * **is anything non-finite?**  ``i``, ``q``, ``G`` and ``C`` are
      scanned entry by entry, because a Jacobian usually goes NaN a long
      way before the value does, and Newton then walks off silently.

    ``x`` is the element's own unknown vector; :func:`x_layout` says what
    each entry is.  ``epar`` defaults to ``defaultepar`` (300 K,
    transient); pass a DC ``epar`` to check the DC variant of a model
    with pinned states, since ``i``/``G`` switch together on it.

    The step is ``max(1e-7, 1e-7*|x[k]|)`` per column unless ``h`` says
    otherwise.  The tolerance is ONE band for the whole matrix --
    ``atol`` defaults to ``1e-7`` of its largest entry -- deliberately: a
    per-entry relative tolerance passes an entry that is small *because
    it is wrong*.  Tighten `rtol` and watch a known-bad model fail before
    trusting a pass.

    **A third verdict: UNRESOLVED.**  A finite difference is not an
    oracle, and three separate mechanisms make it report a large
    discrepancy on a model that is right.  Each has been hit by a real
    model in this tree, so each is measured per entry rather than
    accommodated by a hand-written ``atol``:

    * **roundoff** -- the entry's own signal is below the representable
      step of the value it differentiates, and the difference returns
      literal zero.  The ordinary floor is ``eps * max|f| / h``, but that
      is derived from the OUTPUTS and a cancellation inside the model is
      not visible there: EKV's ``q`` at deep cutoff is 4.2e-27 with a
      real quantum of 4.1e-31, `eps` times an internal magnitude of
      1.9e-15, so ``eps*|q|/h`` misses by ten decades.  For an entry that
      does not move AT ALL the floor is therefore MEASURED -- the step is
      widened by a decade at a time until the value clears its own
      quantisation.  (Deep-cutoff EKV: ``C`` entries of 1.6e-25 F against
      a difference of exactly 0.0, which widening turns into -1.26e-25
      against an analytic -1.28e-25.)
    * **truncation** -- a stiff card where the second difference is not
      small.  Estimated per entry by RICHARDSON, not bounded a priori:
      the same column is differenced at ``h`` and at ``2h``, and their
      disagreement is three times the ``h^2`` term.  (A memristor with
      ``dR/dx = 1e9``, where a 1e-7 step is a 1% change in ``R``.)
    * **kink** -- the value is C0 but not C1 here, so a central
      difference returns the AVERAGE of the two one-sided slopes while
      the analytic Jacobian returns one of them.  No ``h`` helps: a jump
      has no scale.  Detected from ``f`` ALONE, by whether the one-sided
      disagreement ``|f(x+h) - 2f(x) + f(x-h)| / h`` shrinks with ``h``
      (smooth: it halves; kink: it does not).

    Note what that last one can and cannot say.  It detects that the
    VALUE is kinked, which is a fact about the model independent of any
    Jacobian, and it therefore does not excuse an error larger than the
    jump.  It does hide an error SMALLER than the jump, and there is no
    way around that: at a kink the derivative is genuinely
    two-valued and a difference cannot pick.  So the honest report is
    "not resolvable here", which is what this returns.

    An UNRESOLVED entry does not fail: :attr:`JacobianCheck.ok` stays
    true.  :attr:`JacobianCheck.resolved` is the stronger statement, and
    :attr:`JacobianCheck.unresolved` lists what was skipped and why.
    **The floors are per entry and never widen a resolvable one** -- on a
    diode at 0.42 V the band is 4.4e-12 and the floors are 2.5e-15
    (roundoff) and 7.1e-15 (truncation), three decades under -- so a
    deliberately wrong Jacobian still FAILS.  There is a fourth
    condition, with no floor at all: when the ``h^2`` term is a quarter
    of the difference itself the expansion is not converging at this
    step, and the entry is unresolved outright.

    Returns a :class:`JacobianCheck`; ``print()`` it.
    """
    cls, _info = _hdl_info_of(element)
    if isinstance(element, type):
        raise TypeError('check_jacobians needs an element INSTANCE (it '
                        'evaluates it); got the class %s.' % cls.__name__)
    if epar is None:
        epar = defaultepar
    layout = x_layout(element)
    x = np.asarray(x, dtype=float)
    n = len(layout)
    if x.shape != (n,):
        raise ValueError(
            '%s takes an x of length %d, got %s. The layout is:\n%s'
            % (cls.__name__, n, list(np.shape(x)),
               '\n'.join('  x[%d]  %-12s %s' % e for e in layout)))

    def ev(which, xv):
        return np.asarray(getattr(element, which)(xv, epar), dtype=float)

    results, nonfinite = {}, []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        with np.errstate(all='ignore'):
            for which, against in (('G', 'i'), ('C', 'q')):
                ana = ev(which, x).reshape(n, n)
                f0 = ev(against, x)
                fd = np.zeros((n, n))          # central difference at h
                fd2 = np.zeros((n, n))         # ... and at 2h, for Richardson
                rnd = np.zeros((n, n))         # the roundoff floor
                trunc = np.zeros((n, n))       # the truncation floor
                asym1 = np.zeros((n, n))       # |f(x+h) - 2f(x) + f(x-h)| / h
                asym2 = np.zeros((n, n))       # the same at 2h
                for k in range(n):
                    hk = h if h is not None else max(1e-7, 1e-7 * abs(x[k]))

                    def col(s, _k=k, _h=hk):
                        xs = x.copy()
                        xs[_k] += s * _h
                        return ev(against, xs)

                    fp1, fm1, fp2, fm2 = col(1.), col(-1.), col(2.), col(-2.)
                    fd[:, k] = (fp1 - fm1) / (2 * hk)
                    fd2[:, k] = (fp2 - fm2) / (4 * hk)
                    ## The rounding of a DIFFERENCE is set by the size of
                    ## the operands, not of the result -- that is exactly
                    ## the case where the difference comes back 0.0.
                    mag = np.max(np.abs(np.vstack([f0, fp1, fm1, fp2, fm2])),
                                 axis=0)
                    rnd[:, k] = _JAC_EPS * mag / hk
                    trunc[:, k] = 2.0 * np.abs(fd2[:, k] - fd[:, k]) / 3.0
                    ## One-sided disagreement.  Smooth: |f''|*h, so it
                    ## HALVES with h.  Kinked: |f'(+) - f'(-)|, constant.
                    asym1[:, k] = np.abs(fp1 - 2 * f0 + fm1) / hk
                    asym2[:, k] = np.abs(fp2 - 2 * f0 + fm2) / (2 * hk)
                    ## ---- FROZEN entries: WIDEN THE STEP AND MEASURE ---
                    ## `f` came back bitwise identical at all four probe
                    ## points while the analytic entry is not zero.  The
                    ## roundoff formula above cannot size this, and
                    ## MEASUREMENT SAYS SO: for EKV's normalised charges
                    ## in deep cutoff the true quantum is 4.1e-31, which
                    ## is `eps` times an INTERNAL magnitude of 1.9e-15 --
                    ## ten decades above `eps*|q|` and 46x above
                    ## `eps*max|q|`.  A cancellation inside the model is
                    ## not visible from its outputs, so the floor has to
                    ## be probed rather than derived.
                    ##
                    ## And widening is a real discriminator, not a
                    ## whitewash: a value that is genuinely INDEPENDENT
                    ## of `x[k]` stays frozen at every step up to the
                    ## cap, and then a non-zero analytic entry still
                    ## FAILS.  Only a value that does move, once its
                    ## signal is clear of its own quantum, is excused.
                    frozen = ((fp1 == f0) & (fm1 == f0) & (fp2 == f0)
                              & (fm2 == f0) & (ana[:, k] != 0.0)
                              & np.isfinite(f0))
                    if np.any(frozen):
                        _jac_widen(col, hk, frozen, x[k],
                                   fd[:, k], rnd[:, k], trunc[:, k])
                scale = float(max(np.max(np.abs(ana[np.isfinite(ana)]))
                                  if np.any(np.isfinite(ana)) else 0.0,
                                  np.max(np.abs(fd[np.isfinite(fd)]))
                                  if np.any(np.isfinite(fd)) else 0.0))
                at = atol if atol is not None else 1e-7 * scale
                err = np.abs(ana - fd)
                err[~np.isfinite(err)] = np.inf
                worst = np.unravel_index(int(np.argmax(err)), err.shape)
                finite = bool(np.all(np.isfinite(ana))
                              and np.all(np.isfinite(fd)))
                ## --- the three noise floors, per entry -----------------
                ## A kink is a fact about `f` alone: the one-sided
                ## disagreement stays PUT as `h` shrinks, where a smooth
                ## function's halves.  The band is around 1 rather than
                ## merely below 2, and that matters: a violently curved
                ## value gives 0.5, and calling that a kink excused a
                ## real 200x discrepancy on a stiff memristor card.
                kink = (asym1 > 16.0 * rnd) & (asym1 > 0.0) \
                    & (asym2 <= 1.4 * asym1) & (asym2 >= 0.7 * asym1)
                ## The error a kink puts into a central difference is HALF the
                ## jump -- the difference returns the average of the two arms
                ## and the Jacobian returns one of them -- so the floor is
                ## that with 50% of margin, not the whole jump.  Measured on
                ## the memristor's clamp corners: |diff| = 2.17e-1 against an
                ## `asym1` of 4.35e-1, exactly the half.
                kinkf = np.where(kink, 0.75 * asym1, 0.0)
                ## OUTSIDE THE ASYMPTOTIC REGIME there is no floor to
                ## quote.  Richardson measures the `h^2` term; when that
                ## term is a sizeable fraction of the difference itself,
                ## the Taylor expansion the whole method rests on is not
                ## converging at this step, and the difference is not an
                ## estimate of anything.  Measured on a memristor with
                ## `ron = 1, roff = 1e9` at the clamp corner: `R` moves
                ## by 100x across one 1e-7 step, the difference reports
                ## 3.5e6 against a correct 7.0e8, and no finite bound
                ## covers that honestly.
                nonasym = (trunc > 0.25 * np.abs(fd)) & (trunc > at) \
                    & (trunc > 4.0 * rnd)
                floors = np.dstack([rnd, np.where(nonasym, np.inf, trunc),
                                    kinkf])
                floors[np.isnan(floors)] = 0.0         # a NaN excuses nothing
                extra = np.max(floors, axis=2)
                which_f = np.argmax(floors, axis=2)
                names = ('roundoff', 'truncation', 'kink')
                reason = [[names[int(which_f[i, j])] for j in range(n)]
                          for i in range(n)]
                tol_pass = at + rtol * np.abs(fd)
                tol_unres = np.maximum(at, extra) + rtol * np.abs(fd)
                status = np.where(err <= tol_pass, JAC_PASS,
                                  np.where(err <= tol_unres, JAC_UNRESOLVED,
                                           JAC_FAIL))
                if scale == 0.0:
                    status[:] = JAC_PASS
                nunres = int(np.count_nonzero(status == JAC_UNRESOLVED))
                ## The headline entry is the worst RESOLVED one: reporting
                ## an unresolved entry as "the worst" would put a number
                ## nobody should read at the top of the report.
                errr = np.where(status == JAC_UNRESOLVED, -1.0, err)
                wres = (worst if np.all(errr < 0) else
                        np.unravel_index(int(np.argmax(errr)), err.shape))
                results[which] = dict(
                    ana=ana, fd=fd, scale=scale, atol=at, finite=finite,
                    err=float(err[worst]), worst=tuple(int(k) for k in worst),
                    err_mat=err, status=status, floor=tol_unres,
                    noise=extra, roundoff=rnd, truncation=trunc, kink=kink,
                    reason=reason, nunresolved=nunres,
                    worst_resolved=tuple(int(k) for k in wres),
                    ok=bool(finite and not np.any(status == JAC_FAIL)))
                for idx in zip(*np.nonzero(~np.isfinite(ana))):
                    nonfinite.append((which, tuple(int(k) for k in idx)))
                ## The finite difference too, or a pair can be reported
                ## NOT COMPARABLE with nothing listed under it: `i` can
                ## be finite at `x` and infinite one step away, which is
                ## itself worth seeing.
                for idx in zip(*np.nonzero(~np.isfinite(fd))):
                    nonfinite.append(('d%s/dx (finite difference)' % against,
                                      tuple(int(k) for k in idx)))
            for which in ('i', 'q'):
                v = ev(which, x)
                for idx in zip(*np.nonzero(~np.isfinite(v))):
                    nonfinite.append((which, tuple(int(k) for k in idx)))
    return JacobianCheck(cls.__name__, x, layout, results, nonfinite,
                         rtol, atol if atol is not None else
                         max(r['atol'] for r in results.values()))


def isconstant(expr):
    for atom in expr.atoms():
        if isinstance(atom, Quantity):
            return False
    return True


if __name__ == "__main__":
    import doctest
    doctest.testmod()
