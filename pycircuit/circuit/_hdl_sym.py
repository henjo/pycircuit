"""The HDL language's sympy classes and symbols.

Defined here, not in `hdl`, so that importing the language does not import
sympy: a numeric run of cached models builds no expression, and sympy's import
was ~0.15 s of a short deck's whole process.  `hdl` binds each of these names
as a `_lazy` stand-in, replaced by the object defined here the moment sympy is
imported; `hdl` is the place to import them from.
"""

import sympy
from sympy.core.symbol import Str

from pycircuit.circuit.hdl import (
    Branch,
    Node,
    _maxc_numpy,
    _minc_numpy,
    _noise_name,
    _quartic_root_numpy,
    safe_div,
)


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


class quartic_root(sympy.Function):
    """SPICE's MOS level-2 velocity-saturation root (`_quartic_root_scalar`):
    the smallest positive real root ``x`` of ``x^4 + a1 x^3 + b1 x^2 + c1 x
    + d1`` that SPICE's closed form finds within ``|P| <= 1e-6``, or 0.0
    where there is none (the caller keeps its other ``vdsat`` then).

    Its derivative is IMPLICIT, ``P(x(a), a) = 0``: ``dx/da_k =
    -(dP/da_k) / P'(x)``, ``dP/d(a1, b1, c1, d1) = (x^3, x^2, x, 1)``,
    ``P'(x) = 4x^3 + 3 a1 x^2 + 2 b1 x + c1`` -- exact where the selected
    root moves smoothly, through `safe_div` so the arm a caller does not
    select (``x = 0``) stays finite.  (2026-10-09, the SPICE benchmark
    plan's stage 9.)"""
    nargs = (4,)

    def fdiff(self, argindex=1):
        a1, b1, c1, _d1 = self.args
        x = self
        dp = 4 * x ** 3 + 3 * a1 * x ** 2 + 2 * b1 * x + c1
        num = (x ** 3, x ** 2, x, sympy.Integer(1))[argindex - 1]
        return -safe_div(num, dp, eps=1e-100)

    @staticmethod
    def _imp_(a1, b1, c1, d1):
        return _quartic_root_numpy(a1, b1, c1, d1)


class _Laplace(sympy.Function):
    """Marker for a laplace_* application; see :func:`laplace_nd`."""
    nargs = (3,)


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


class _LimitDelta(_Limit):
    """`$limit(probe, "deltalim", vmax)`; see :func:`limit_delta`."""
    nargs = (2, 5)
    kind = 'delta'


class _LimitId(_Limit):
    """An identity probe -- a limited unknown with NO law; see
    :func:`limit_identity`."""
    nargs = (1, 4)
    kind = 'id'


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


class _AcStim(sympy.Function):
    nargs = (2,)


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


#: The atoms `select` will substitute for a clamped copy: a parameter or
#: intermediate symbol, or a branch/node quantity.  Anything else on one
#: side of a relational is a compound, and a compound cannot be clamped
#: by substitution -- clamping `a*b` would mean choosing which factor to
#: move.
_CLAMPABLE = (sympy.Symbol, Quantity)


class _StateSymbol(sympy.Symbol):
    """Symbol standing for an idt/idtmod state unknown."""

    def __new__(cls, name):
        ## real=True for the same reason as Quantity: a Piecewise over a
        ## state must be constructible.
        return super().__new__(cls, name, real=True)

    def __getnewargs_ex__(self):        # see _IntermediateSymbol
        return ((self.name,), {})
