"""`hdl.select` -- the compiler's own both-arms domain propagation.

Both arms of a compiled conditional are evaluated: `Piecewise` becomes
`where`, which picks AFTERWARDS.  So an arm that is not selected at this
bias still runs, still raises floating-point flags, and still has its
derivative taken -- and the derivative is the half that matters, because
the chain rule multiplies a zero partial by that arm's derivative and
``0 * nan`` is ``nan``.  The model obeys that today by hand, with a
clamped copy of the input fed to the arm that does not want it.

`select` derives the clamp from the condition instead.  The claims under
test, in the order they are made:

* **the clamp is the identity where the arm is selected** -- in value
  AND in derivative, bit for bit, for every supported shape.  That is
  what makes a conversion a no-op on the answer, and it is asserted by
  sweeping the two spellings against each other and comparing bytes on
  the selected mask;
* **the discarded arm is finite** in value and in derivative, under
  ``numpy.errstate(all='raise')``, where the unclamped spelling is not.
  The unclamped spelling is run in the same test, and must fail -- a
  test that cannot fail is not evidence;
* **the refusals are loud and specific.**  A `select` that quietly did
  not clamp would be worse than the `Piecewise` it replaced, because
  the author would trust it, so every unsupported shape raises
  `SelectRefused` naming the shape.  One refusal per unsupported class,
  message pinned;
* **both code generators and both backends** produce the clamped form:
  the eager `lambdify` path, the `var()` let-chain, and the C printer;
* **the six adopted sites in `elements_hdl` did not move a bit.**  The
  record below was taken from 8db5398 -- the commit before `select`
  existed -- over twelve class/card combinations and 39 points each
  including ``+-1e30``.
"""
import hashlib
import warnings

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   SelectRefused, check_jacobians, maxc,
                                   minc, select, unclamped, var)
from pycircuit.utilities.param import Parameter

try:
    from pycircuit.circuit import _hdl_cbackend as cb
    _CC = cb.find_compiler()[0]
except Exception:                                   # pragma: no cover
    _CC = None
needs_cc = pytest.mark.skipif(_CC is None, reason='no C compiler')


## ----------------------------------------------------------------------
## Helpers.

def _strict():
    """``errstate`` that raises on divide, overflow and invalid.

    Underflow is deliberately NOT raised.  `exp` of a large negative
    argument underflows to zero, which is the right answer and happens
    at every bias a clamped exponential is discarded at; raising on it
    would make this suite fail on correct code.  The three that are
    raised are the ones that mean "this arm produced a number the
    Jacobian cannot use"."""
    return np.errstate(divide='raise', over='raise', invalid='raise')


def _f(expr, sym):
    """`expr` and ``d expr/d sym`` as plain numpy callables.

    A PLAIN `lambdify` with no modules map, on purpose: the kernel's
    primitives carry `_imp_`, and a clamp that only worked inside the
    compiler's namespace would be a clamp no outside caller could use.
    """
    return (sympy.lambdify(sym, expr, 'numpy'),
            sympy.lambdify(sym, sympy.diff(expr, sym), 'numpy'))


def _vd(expr, sym, xs):
    """Value and derivative arrays, warnings off."""
    fv, fd = _f(expr, sym)
    with np.errstate(all='ignore'):
        return (np.asarray(fv(xs), float) * np.ones_like(xs),
                np.asarray(fd(xs), float) * np.ones_like(xs))


V = sympy.Symbol('v', real=True)
C = sympy.Symbol('c', real=True)

#: A sweep that lands ON every boundary the shapes below use, on both
#: sides of it, and at the magnitudes the kernel's guards exist for.
XS = np.array([-1e30, -1e10, -3.0, -1.0, -0.5, -1e-3, -1e-30, 0.0,
               1e-30, 1e-3, 0.5, 0.9999, 1.0, 1.0001, 2.0, 3.0,
               1e10, 1e30])

_M = 1e-6
_IN = minc(maxc(V, -1.0 + _M), 1.0 - _M)

#: (id, condition, margin, the arm AS WRITTEN, the hand-clamped arm the
#: model would have written instead, the region the identity is claimed
#: on).  Every arm is singular outside its own condition -- "the clamp
#: did nothing" is never a passing state.
SHAPES = [
    ('lt', V < 1.0, _M,
     sympy.log(1.0 - V), sympy.log(1.0 - minc(V, 1.0 - _M)),
     lambda x: x <= 1.0 - _M),
    ('le', V <= 1.0, 0.0,
     1.0 / (2.0 - V), 1.0 / (2.0 - minc(V, 1.0)),
     lambda x: x <= 1.0),
    ('gt', V > 1.0, _M,
     sympy.log(V - 1.0), sympy.log(maxc(V, 1.0 + _M) - 1.0),
     lambda x: x >= 1.0 + _M),
    ('ge', V >= 1.0, 0.0,
     1.0 / V, 1.0 / maxc(V, 1.0),
     lambda x: x >= 1.0),
    ('reversed', 1.0 < V, _M,
     sympy.log(V - 1.0), sympy.log(maxc(V, 1.0 + _M) - 1.0),
     lambda x: x >= 1.0 + _M),
    ('abs', sympy.Abs(V) < 1.0, _M,
     sympy.log(1.0 - V * V), sympy.log(1.0 - _IN ** 2),
     lambda x: np.abs(x) <= 1.0 - _M),
    ('square', V * V < 1.0, _M,
     sympy.log(1.0 - V * V), sympy.log(1.0 - _IN ** 2),
     lambda x: np.abs(x) <= 1.0 - _M),
    ('and', sympy.And(V > -1.0, V < 1.0), _M,
     sympy.log(1.0 - V * V), sympy.log(1.0 - _IN ** 2),
     lambda x: np.abs(x) <= 1.0 - _M),
]


@pytest.mark.parametrize('name,cond,margin,arm,hand,sel', SHAPES,
                         ids=[s[0] for s in SHAPES])
def test_select_emits_the_clamp_the_model_would_have_written(name, cond,
                                                             margin, arm,
                                                             hand, sel):
    """The emitted tree, not merely the emitted number: a conversion is
    only a no-op on the answer if the expression is the same one."""
    got = select((arm, cond), (sympy.Float(0.0), True), margin=margin)
    assert got.args[0].expr == hand, (name, got.args[0].expr)


@pytest.mark.parametrize('name,cond,margin,arm,hand,sel', SHAPES,
                         ids=[s[0] for s in SHAPES])
def test_clamp_is_the_identity_where_the_arm_is_selected(name, cond, margin,
                                                         arm, hand, sel):
    """Value and derivative, BYTE for byte, against the arm AS WRITTEN,
    on the region the arm is selected in.

    The tie is deliberately in the sweep: `minc(v, c)` differentiates to
    ``_step(c, v)``, which is 1 AT ``v == c``, so a non-strict condition
    keeps its boundary exactly -- and so does the bound's own partial,
    which is exactly zero there."""
    got = select((arm, cond), (sympy.Float(0.0), True), margin=margin)
    mask = sel(XS)
    assert mask.any()
    cv, cd = _vd(got.args[0].expr, V, XS)
    rv, rd = _vd(arm, V, XS)
    assert cv[mask].tobytes() == rv[mask].tobytes(), name
    assert cd[mask].tobytes() == rd[mask].tobytes(), name
    ## ...and the identity is not vacuous: the clamp DOES change the
    ## arm somewhere outside the region.
    assert cv.tobytes() != rv.tobytes(), name


@pytest.mark.parametrize('name,cond,margin,arm,hand,sel', SHAPES,
                         ids=[s[0] for s in SHAPES])
def test_discarded_arm_is_finite_in_value_and_derivative(name, cond, margin,
                                                         arm, hand, sel):
    """...everywhere, under an errstate that raises -- and the arm AS
    WRITTEN is not, which is what makes this a test."""
    got = select((arm, cond), (sympy.Float(0.0), True), margin=margin)
    fv, fd = _f(got.args[0].expr, V)
    with _strict():
        v = np.asarray(fv(XS), float)
        d = np.asarray(fd(XS), float)
    assert np.isfinite(v).all(), (name, v)
    assert np.isfinite(d).all(), (name, d)

    ## The control.
    rv, rd = _f(arm, V)
    bad = False
    try:
        with _strict():
            v2 = np.asarray(rv(XS), float) * np.ones_like(XS)
            d2 = np.asarray(rd(XS), float) * np.ones_like(XS)
        bad = not (np.isfinite(v2).all() and np.isfinite(d2).all())
    except FloatingPointError:
        bad = True
    assert bad, '%s: the unclamped arm is finite anyway' % name


def test_the_bound_may_be_an_expression_and_may_mention_the_symbol():
    """``vds < vdsat`` with `vdsat` a compound, and the pathological
    ``v < 2*v``: the substitution is simultaneous, so each bound is
    evaluated at the UNCLAMPED symbol and no cycle is created."""
    got = select((1.0 / (C - V), V < C), (sympy.Float(0.0), True))
    ## Both sides are atoms, so BOTH are clamped -- each by the other,
    ## as they stood before the substitution.
    assert got.args[0].expr == 1.0 / (maxc(C, V) - minc(V, C))
    ## Self-referential bound: one level, no recursion.
    got = select((V, V < 2 * V), (sympy.Float(0.0), True))
    assert got.args[0].expr == minc(V, 2 * V)


def test_both_sides_atoms_is_the_identity_on_the_region():
    e = select((1.0 / (C - V), V < C), (sympy.Float(0.0), True),
               margin=1e-6)
    f = sympy.lambdify((V, C), e.args[0].expr, 'numpy')
    v = np.array([-1e30, -1.0, 0.0, 0.5, 1.0, 2.0, 1e30])
    c = np.full_like(v, 1.0)
    with _strict():
        got = f(v, c)
    sel = v < c
    with np.errstate(all='ignore'):
        want = 1.0 / (c - v)
    assert got[sel].tobytes() == want[sel].tobytes()
    assert np.isfinite(got).all()


def test_the_default_arm_is_clamped_by_the_complement():
    """A `True` arm has no condition of its own; its region is the
    complement of every earlier one, and that is what clamps it."""
    got = select((sympy.Float(0.0), V <= 0.0), (1.0 / V, True),
                 margin=1e-30)
    assert got.args[1].expr == 1.0 / maxc(V, 1e-30)
    ## Three arms: the middle one is clamped on BOTH sides, by its own
    ## condition and by the negation of the first.
    got = select((sympy.Float(0.0), V < 1.0), (1.0 / V, V < 2.0),
                 (sympy.Float(1.0), True))
    assert got.args[1].expr == 1.0 / minc(maxc(V, 1.0), 2.0)


def test_the_punched_hole_is_exact_at_its_own_boundary():
    """``|v| < c`` in an earlier arm leaves the default arm selected on
    a UNION of two half-lines -- not a box.  `_punch` reaches it, and
    is the identity in value and derivative on the CLOSED region, tie
    included; that is the argument order of its `maxc`."""
    got = select((sympy.Float(0.0), sympy.Abs(V) < 1e-3),
                 (1.0 / V, True))
    fv, fd = _f(got.args[1].expr, V)
    with _strict():
        v = np.asarray(fv(XS), float)
        d = np.asarray(fd(XS), float)
    assert np.isfinite(v).all()
    assert np.isfinite(d).all()
    sel = np.abs(XS) >= 1e-3
    assert v[sel].tobytes() == (1.0 / XS[sel]).tobytes()
    assert d[sel].tobytes() == (-1.0 / XS[sel] ** 2).tobytes()


def test_punch_is_exactly_the_identity_at_the_tie():
    """The tie ``|v| == c`` is a SELECTED point (the earlier arm's
    condition was strict), and `_punch` written with the arguments the
    other way round would give it derivative 0 instead of 1."""
    p = hdl._punch(V, sympy.Float(1e-3))
    fv, fd = _f(p, V)
    xs = np.array([-1e-3, 1e-3])
    assert fv(xs).tobytes() == xs.tobytes()
    assert list(fd(xs)) == [1.0, 1.0]
    wrong = V + (2.0 * hdl._step(V, 0) - 1.0) * maxc(1e-3 - hdl.Abs(V), 0)
    assert list(sympy.lambdify(V, sympy.diff(wrong, V), 'numpy')(xs)) \
        == [0.0, 0.0]


## ----------------------------------------------------------------------
## The margin, and the trade it makes.

def test_margin_is_what_makes_a_divided_by_arm_finite():
    """Clamping to the boundary is not always enough: ``p > 0`` clamps
    to ``maxc(p, 0)``, and an arm that divides by `p` still divides by
    zero.  Both halves are asserted, because the default is zero and an
    author who does not think about it gets the first."""
    at_bound = select((sympy.Float(0.0), V <= 0.0), (1.0 / V, True))
    f = sympy.lambdify(V, at_bound.args[1].expr, 'numpy')
    with np.errstate(all='ignore'):
        assert not np.isfinite(f(np.array([-1.0, 0.0]))).all()
    with_margin = select((sympy.Float(0.0), V <= 0.0), (1.0 / V, True),
                         margin=1e-30)
    f = sympy.lambdify(V, with_margin.args[1].expr, 'numpy')
    with _strict():
        assert np.isfinite(f(np.array([-1.0, 0.0, 1e30]))).all()


def test_margin_costs_the_identity_on_its_own_sliver_and_nothing_else():
    """The recorded price of a margin: the clamp is no longer the
    identity on ``0 < p < margin``.  Named so that it is a decision and
    not a surprise."""
    e = select((sympy.Float(0.0), V <= 0.0), (V, True), margin=1e-3)
    f = sympy.lambdify(V, e.args[1].expr, 'numpy')
    xs = np.array([1e-6, 1e-4, 1e-3, 1e-2, 1.0])
    got = f(xs)
    assert (got[xs >= 1e-3] == xs[xs >= 1e-3]).all()
    assert (got[xs < 1e-3] == 1e-3).all()


def test_a_margin_is_never_applied_to_a_non_strict_bound():
    """A non-strict condition has its boundary INSIDE the selected
    region, where moving it would change an answer."""
    e = select((V, V <= 1.0), (sympy.Float(0.0), True), margin=0.25)
    assert e.args[0].expr == minc(V, 1.0)
    e = select((V, V < 1.0), (sympy.Float(0.0), True), margin=0.25)
    assert e.args[0].expr == minc(V, 0.75)


def test_a_non_strict_condition_cannot_be_given_a_margin_and_sqrt_needs_one():
    """The limit of the whole idea, pinned rather than glossed.

    ``select((sqrt(v), v >= 0), ...)`` clamps to ``sqrt(maxc(v, 0))``,
    whose derivative is ``_step(v, 0)/(2*sqrt(maxc(v, 0)))`` -- ``0/0``
    for every ``v < 0`` and ``1/0`` at zero.  A margin would fix it and
    `select` refuses to apply one, because the condition is NON-STRICT:
    ``v = 0`` is a point the arm is SELECTED at, and moving the clamp
    there would change an answer rather than a discarded one.

    So this shape is not a clamping problem at all.  The arm is
    singular at a bias it is chosen at, and the remedy is the one
    `differentiable-numerics` gives -- smooth instead of clip, i.e.
    `safe_sqrt` -- which is asserted here so that the next author does
    not go looking for a margin that cannot exist."""
    at_bound = select((sympy.sqrt(V), V >= 0.0), (sympy.Float(0.0), True),
                      margin=1e-30)
    assert at_bound.args[0].expr == sympy.sqrt(maxc(V, 0.0))
    _fv, fd = _f(at_bound.args[0].expr, V)
    with np.errstate(all='ignore'):
        assert not np.isfinite(fd(np.array([-1.0, 0.0]))).all()
    ## The remedy, and it needs no clamp at all.
    _fv, fd = _f(hdl.safe_sqrt(V), V)
    with _strict():
        assert np.isfinite(fd(np.array([-1e30, -1.0, 0.0, 1e30]))).all()
    ## A STRICT condition does take the margin, and then the clamp is
    ## enough on its own.
    ok = select((sympy.sqrt(V), V > 0.0), (sympy.Float(0.0), True),
                margin=1e-30)
    assert ok.args[0].expr == sympy.sqrt(maxc(V, 1e-30))
    _fv, fd = _f(ok.args[0].expr, V)
    with _strict():
        assert np.isfinite(fd(np.array([-1.0, 0.0, 1e30]))).all()


## ----------------------------------------------------------------------
## The refusals.  One per unsupported class, message pinned.

REFUSALS = [
    ('or', sympy.Or(V < 0.0, V > 3.0), 'UNION of boxes'),
    ('eq', sympy.Eq(V, 0.0), 'measure-zero set'),
    ('ne', sympy.Ne(V, 0.0), 'measure-zero set'),
    ('product', V * C < 3.0, 'neither side is a symbol'),
    ('sum', V + C < 3.0, 'neither side is a symbol'),
    ('symbolic-square', V * V < 3 * C, 'not a non-negative literal'),
    ('not-a-relational', V, 'not a relational'),
]


@pytest.mark.parametrize('name,cond,msg', REFUSALS,
                         ids=[r[0] for r in REFUSALS])
def test_unsupported_shapes_are_refused_by_name(name, cond, msg):
    with pytest.raises(SelectRefused) as exc:
        select((1.0 / V, cond), (sympy.Float(0.0), True))
    assert msg in str(exc.value), str(exc.value)
    ## ...and the escape hatch takes the same arm without complaint.
    got = select((unclamped(1.0 / V), cond), (sympy.Float(0.0), True))
    assert got.args[0].expr == 1.0 / V


@pytest.mark.parametrize('name,cond,msg', REFUSALS,
                         ids=[r[0] for r in REFUSALS])
def test_strict_false_accepts_every_refused_shape_unclamped(name, cond, msg):
    got = select((1.0 / V, cond), (sympy.Float(0.0), True), strict=False)
    assert got.args[0].expr == 1.0 / V


def test_a_default_arm_over_a_hole_it_uses_is_refused():
    """``|v| < c`` on an earlier arm is fine -- the punch reaches it --
    but ``And(...)`` negates to an `Or`, and no substitution is the
    identity on a union.  Refused only when the default arm actually
    USES one of the symbols, because otherwise nothing was lost."""
    with pytest.raises(SelectRefused) as exc:
        select((sympy.Float(0.0), sympy.And(V > 0.0, C > 0.0)),
               (1.0 / V, True))
    assert 'UNION of boxes' in str(exc.value)
    assert 'this arm uses v' in str(exc.value)
    ## The same shape, an arm that does not touch `v` or `c`: allowed.
    other = sympy.Symbol('w', real=True)
    got = select((sympy.Float(0.0), sympy.And(V > 0.0, C > 0.0)),
                 (1.0 / other, True))
    assert got.args[1].expr == 1.0 / other


def test_a_box_and_a_hole_on_one_symbol_are_refused():
    with pytest.raises(SelectRefused) as exc:
        select((sympy.Float(0.0), sympy.Abs(V) < 1.0),
               (1.0 / V, V < 5.0), (sympy.Float(0.0), True))
    assert 'box and out of a hole' in str(exc.value)


def test_or_is_supported_in_a_negated_position():
    """De Morgan: the complement of a union IS an intersection of
    boxes, so an earlier `Or` clamps the arms after it."""
    got = select((sympy.Float(0.0), sympy.Or(V < 1.0, C < 1.0)),
                 (1.0 / (V * C), True))
    assert got.args[1].expr == 1.0 / (maxc(V, 1.0) * maxc(C, 1.0))


def test_nested_select_clamps_at_each_level():
    inner = select((1.0 / V, V > 1.0), (sympy.Float(0.0), True))
    got = select((inner, V < 5.0), (sympy.Float(0.0), True))
    arm = got.args[0].expr
    ## The outer clamp reached INSIDE the inner Piecewise, arms and
    ## conditions alike -- which is what makes a nest safe.
    assert arm.args[0].expr == 1.0 / maxc(minc(V, 5.0), 1.0)
    assert arm.args[0].cond == (minc(V, 5.0) > 1.0)
    f = sympy.lambdify(V, got, 'numpy')
    with _strict():
        assert np.isfinite(f(XS)).all()


def test_bad_arguments_are_type_errors_not_silence():
    with pytest.raises(TypeError):
        select(V, (sympy.Float(0.0), True))
    with pytest.raises(TypeError):
        select((V, V < 1.0), nosuchkw=1)
    with pytest.raises(TypeError):
        select()


## ----------------------------------------------------------------------
## The one refusal that is about the COMPILER, not about the condition.

def _shadow_class(name, body):
    """A `Behavioural` built INSIDE a test.

    The metaclass compiles at class-creation time, so a class whose
    `select` is meant to be refused cannot sit at module level: it would
    raise during collection and take the whole file with it."""
    return type(name, (Behavioural,),
                dict(instparams=[Parameter(name='k', desc='k', unit='',
                                           default=1.0)],
                     analog=staticmethod(body), __module__=__name__))


def test_a_clamp_shadowed_by_an_intermediate_is_refused():
    """An arm that reaches its constrained symbol only through a
    `var()` intermediate.  `select` cannot follow a symbol into a
    definition that was already bound, so the clamp would be a silent
    no-op -- and that is the failure this whole exercise exists to
    prevent."""
    def analog(plus, minus):
        b = Branch(plus, minus)
        d = var(1.0 / k, 'd')                                  # noqa: F821
        return Contribution(b.I, select((d * b.V, k > 0.0),    # noqa: F821
                                        (sympy.Float(0.0), True)))
    with pytest.raises(SelectRefused) as exc:
        _shadow_class('_Shadowed', analog)(Node('a'), Node('b'))
    m = str(exc.value)
    assert 'through the var() intermediate' in m
    assert '_v0_d' in m
    ## The escape hatch takes it, and so does strict=False.
    def analog2(plus, minus):
        b = Branch(plus, minus)
        d = var(1.0 / k, 'd')                                  # noqa: F821
        return Contribution(b.I, select((unclamped(d * b.V),
                                         k > 0.0),             # noqa: F821
                                        (sympy.Float(0.0), True)))
    el = _shadow_class('_ShadowedRaw', analog2)(Node('a'), Node('b'))
    el.update_iparv()
    assert np.isfinite(np.asarray(el.i([0.5, 0.0]), float)).all()


def test_the_shadow_check_does_not_fire_when_the_clamp_bites():
    """The same shape with the symbol ALSO used directly: the clamp
    does something real, so it is allowed (and the intermediate stays
    unclamped, which the arm's author can see)."""
    def analog(plus, minus):
        b = Branch(plus, minus)
        d = var(1.0 / k, 'd')                                  # noqa: F821
        return Contribution(b.I, select((d * b.V * k,          # noqa: F821
                                         k > 0.0),             # noqa: F821
                                        (sympy.Float(0.0), True)))
    el = _shadow_class('_ShadowedButMentioned', analog)(Node('a'), Node('b'))
    el.update_iparv()
    assert np.isfinite(np.asarray(el.i([0.5, 0.0]), float)).all()


## ----------------------------------------------------------------------
## Both code generators, both backends.

_LIM = 700.0


class _OverflowEager(Behavioural):
    """The both-arms overflow in one line: ``exp(v)`` is finite only
    below ~709, and it is evaluated at EVERY bias whether the arm is
    chosen or not.  Written eagerly -- no `var()` -- so it compiles
    through `lambdify` rather than through the let-chain."""
    instparams = [Parameter(name='k', desc='k', unit='', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, select((sympy.exp(b.V), b.V < _LIM),
                                        (b.V, True)))


class _OverflowChained(_OverflowEager):
    """The same model through `var()`, which is a different printer with
    a different table."""

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        u = var(b.V, 'u')
        return Contribution(b.I, select((sympy.exp(u), u < _LIM),
                                        (u, True)))


class _OverflowUnclamped(_OverflowEager):
    """The mutation: the identical model with the clamp taken out.  It
    must FAIL the finiteness test below, or that test proves nothing."""

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        u = var(b.V, 'u')
        return Contribution(b.I, sympy.Piecewise((sympy.exp(u), u < _LIM),
                                                 (u, True)))


def _mk(cls, **kw):
    el = cls(Node('a'), Node('b'), **kw)
    el.update_iparv()
    return el


BIAS = [0.0, 1.0, 100.0, 699.0, 700.0, 701.0, 1e4, 1e30, -1e30]


def _numpy_i(cls, el, x):
    """`i` through the GENERATED NUMPY SOURCE, whatever backend is bound.

    `el.i()` runs whichever kernel `$PYCIRCUIT_HDL_BACKEND` selected, and
    the C one does not honour `numpy.errstate`: an `exp` that overflows
    inside a C kernel sets the host FPU flag and returns `inf` without
    ever reaching numpy, so the discarded arm's overflow is SILENT
    there.  Pinned by `test_the_c_backend_does_not_honour_errstate`
    below.  Tests about raising therefore go through the numpy function,
    which exists in both backend states -- the C kernel is printed FROM
    it."""
    args = [float(v) for v in hdl._args_of(el, defaultepar)]
    return np.asarray(cls._hdl_info['funcs']['i'](np.asarray(x, float),
                                                  *args), float)


@pytest.mark.parametrize('cls', [_OverflowEager, _OverflowChained],
                         ids=['eager', 'chained'])
def test_both_generators_emit_the_clamp(cls):
    """Through `el.i`, i.e. through whichever backend is bound: the
    finiteness claim holds in both states, and on the numpy path the
    errstate makes it a raise rather than an assertion."""
    el = _mk(cls)
    for v in BIAS:
        with _strict():
            i = np.asarray(el.i([v, 0.0]), float)
            G = np.asarray(el.G([v, 0.0]), float)
        assert np.isfinite(i).all(), v
        assert np.isfinite(G).all(), v
    ## ...and the answer is still the model's: below the seam it is
    ## exp, above it is v.
    assert np.isclose(el.i([1.0, 0.0])[0], np.exp(1.0))
    assert np.isclose(el.i([1e4, 0.0])[0], 1e4)


def test_the_finiteness_test_can_fail():
    """Mutation check.  `_OverflowUnclamped` is `_OverflowChained` with
    the clamp removed; it overflows in the discarded arm at every bias
    above the seam, so `test_both_generators_emit_the_clamp` would fail
    on it.  Named tests that would fail: that one, and
    `test_both_backends_agree_bitwise` below."""
    el = _mk(_OverflowUnclamped)
    raised = False
    try:
        with _strict():
            i = _numpy_i(_OverflowUnclamped, el, [1e4, 0.0])
        raised = not np.isfinite(i).all()
    except FloatingPointError:
        raised = True
    assert raised, 'the unclamped model did not overflow -- the ' \
                   'finiteness tests above prove nothing'


@pytest.mark.parametrize('cls', [_OverflowEager, _OverflowChained],
                         ids=['eager', 'chained'])
def test_the_jacobian_of_a_clamped_model_is_right(cls):
    el = _mk(cls)
    for v in (0.5, 1.0, 100.0, 690.0, 1e4):
        res = check_jacobians(el, [v, 0.0])
        assert res.ok, (v, str(res))


@needs_cc
def test_the_c_backend_does_not_honour_errstate():
    """Measured, and worth knowing before trusting a green suite.

    `numpy.errstate` is numpy's, and a C kernel's arithmetic never goes
    through numpy: the same unclamped model that RAISES on the numpy
    path returns `inf` from the discarded arm in C, silently, and the
    `where` then throws it away.  So on the C backend the overflow of a
    discarded arm is invisible -- what remains dangerous is the same
    arm's DERIVATIVE, which is `nan` on both paths.

    The consequence for this file: every test that asserts something
    RAISES goes through `_numpy_i`, not through `el.i`."""
    cls = _OverflowUnclamped
    el = _mk(cls)
    with pytest.raises(FloatingPointError):
        with _strict():
            _numpy_i(cls, el, [1e4, 0.0])
    hdl.set_backend('c', cls)
    try:
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        with _strict():
            i = np.asarray(el.i([1e4, 0.0]), float)
        assert np.isfinite(i).all()
    finally:
        hdl.set_backend(None, cls)


@needs_cc
def test_both_backends_agree_bitwise():
    """The C printer has its own table; a primitive that prints on one
    path and not the other is a compile error at best and a different
    number at worst."""
    cls = _OverflowChained
    el = _mk(cls)
    args = [float(v) for v in hdl._args_of(el, defaultepar)]
    funcs = cls._hdl_info['funcs']
    hdl.set_backend('c', cls)
    try:
        assert cls._hdl_backend_status == 'c', cls._hdl_backend_status
        f = funcs['i']
        kern = f.__dict__['_hdl_c']
        for v in BIAS:
            x = np.ascontiguousarray([v, 0.0])
            with np.errstate(all='ignore'):
                ref = np.asarray(f(x, *args), float)
            got = np.asarray(kern(el, x, defaultepar), float)
            assert ref.tobytes() == got.tobytes(), v
    finally:
        hdl.set_backend(None, cls)


def test_the_clamp_appears_in_the_generated_source():
    """Assert on what it EMITS, not only on what it evaluates to."""
    el = _mk(_OverflowChained)
    src = _OverflowChained._hdl_info['funcs']['i']._src
    assert 'minc' in src
    assert 'exp' in src
    del el


## ----------------------------------------------------------------------
## The adopted sites: byte identity against the commit before `select`.

_T0 = 300.15
_BJT = dict(IS=2e-15, bf=150.0, br=4.0, nf=1.02, nr=1.05, vaf=80.0, var=12.0,
            ikf=0.05, ikr=0.01, ise=3e-14, ne=1.6, isc=1e-13, nc=1.8,
            rb=120.0, rbm=20.0, re=0.8, rc=2.5, cje=1.2e-12, vje=0.75,
            mje=0.35, cjc=0.5e-12, vjc=0.6, mjc=0.3, xcjc=0.7, tf=3e-10,
            xtf=2.0, vtf=5.0, itf=0.02, tr=5e-8, xtb=1.5, eg=1.12, xti=3.2,
            tnom=_T0, area=1.3, rth=150.0, cth=1e-4)
#: `vaf`/`ikf`/`vtf` = 0 is SPICE's "infinite"/"absent" -- the arm the
#: clamp discards, and therefore the card the record most needs.
_BJT_OFF = dict(_BJT, vaf=0.0, var=0.0, ikf=0.0, ikr=0.0, vtf=0.0)
_NMOS = dict(vto=0.75, kp=8e-5, gamma=0.55, phi=0.70, lambd=0.03,
             tox=2e-8, w=20e-6, l=2e-6, ld=1.5e-7,
             cgso=3.5e-10, cgdo=2.5e-10, cgbo=2.0e-10,
             IS=1e-14, pb=0.85, cj=3e-4, cjsw=3e-10, mj=0.5, mjsw=0.33,
             fc=0.5, js=0.0, ad=1.6e-10, asrc=1.6e-10, pd=4.6e-5, ps=4.6e-5,
             rd=0.0, rs=0.0, rsh=0.0, nrd=2.0, nrs=2.0,
             kf=1e-27, af=1.0, tnom=_T0)
_MOS3 = dict(vto=0.7, u0=500.0, tox=2e-8, nsub=1e16, xj=3e-7, nfs=1e11,
             eta=0.05, delta=0.5, theta=0.1, vmax=1e5, kappa=0.3,
             w=10e-6, l=1e-6, ld=0.05e-6, tnom=_T0,
             IS=1e-14, cj=3e-4, cjsw=3e-10, ad=1.6e-11, pd=1.2e-5,
             ps=1.2e-5, cgso=3.5e-10, cgdo=2.5e-10, cgbo=2.0e-10)
_MOS3['as'] = 1.6e-11

#: (label, class, card, pins) -- the classes the six adopted sites are
#: compiled into, each on a card that reaches the adopted arm and on one
#: that takes the OTHER arm.
_ADOPTERS = [
    ('MosLevel1Hdl', 'MosLevel1Hdl', _NMOS, ('d', 'g', 's', 'b')),
    ('MosLevel1Hdl-r', 'MosLevel1Hdl', dict(_NMOS, rsh=25.0),
     ('d', 'g', 's', 'b')),
    ('MosLevel1Hdl-g0', 'MosLevel1Hdl', dict(_NMOS, gamma=0.0),
     ('d', 'g', 's', 'b')),
    ('MosLevel1PmosHdl', 'MosLevel1PmosHdl', _NMOS, ('d', 'g', 's', 'b')),
    ('GummelPoonNpnHdl', 'GummelPoonNpnHdl', _BJT, ('c', 'b', 'e')),
    ('GummelPoonNpnHdl-off', 'GummelPoonNpnHdl', _BJT_OFF, ('c', 'b', 'e')),
    ('GummelPoonPnpHdl', 'GummelPoonPnpHdl', _BJT, ('c', 'b', 'e')),
    ('GummelPoonNpnThermalHdl', 'GummelPoonNpnThermalHdl', _BJT,
     ('c', 'b', 'e', 'th', 'tha')),
    ('GummelPoonNpnThermalHdl-off', 'GummelPoonNpnThermalHdl', _BJT_OFF,
     ('c', 'b', 'e', 'th', 'tha')),
    ('MosLevel3Hdl', 'MosLevel3Hdl', _MOS3, ('d', 'g', 's', 'b')),
    ('MosLevel3Hdl-off', 'MosLevel3Hdl', dict(_MOS3, vmax=0.0, kappa=0.0),
     ('d', 'g', 's', 'b')),
    ('MosLevel3PmosHdl', 'MosLevel3PmosHdl', _MOS3, ('d', 'g', 's', 'b')),
]

#: One digest over the concatenated `i, G, q, C` bytes at each of 39
#: points, RECORDED AT 8db5398 -- before `select` existed -- and
#: unchanged by the conversion of `sarg`, `sargl` (x2), `sqphbs`, `tfx`
#: and `delxl`.  Not a tolerance: the bytes.
_RECORDED = {
    'MosLevel1Hdl': 'e2c6241eb16203a4',
    'MosLevel1Hdl-r': '03330d0a2a9e611b',
    'MosLevel1Hdl-g0': 'a4bdc5d3bf9ff3f7',
    'MosLevel1PmosHdl': 'e62f01daab99cb38',
    'GummelPoonNpnHdl': 'd25a34d69d5efbf9',
    'GummelPoonNpnHdl-off': 'a5b1739425c7d7d8',
    'GummelPoonPnpHdl': 'dcf3f2d4f8a587a4',
    'GummelPoonNpnThermalHdl': 'f33c131ee2e925cd',
    'GummelPoonNpnThermalHdl-off': '7bd7455df9b16475',
    ## ⚠ RE-RECORDED 2026-08-27 for the three MOS level 3 rows, and only
    ## those.  `_autohold` (roadmap sec. 36) makes the regularisers hold
    ## their own arguments, which stops sympy flattening across the
    ## boundary -- so the emitted arithmetic reassociates and the last bit
    ## moves.  MEASURED before re-recording, on a card that reaches the
    ## changed arms (nsub = 1e16, since `nsub` defaults to 0 and a default
    ## card never gets there -- sec. 28's lesson):
    ##
    ##     3200 samples, 32 of 1120 comparable differ
    ##     max relative 2.535e-16, median 1.669e-16  -- ONE ULP
    ##
    ## Old values, kept so the size of the shift stays checkable:
    ##     MosLevel3Hdl       c45f4d5d4a95a6ab
    ##     MosLevel3Hdl-off   3ea9cdcb0bbad2aa
    ##     MosLevel3PmosHdl   7555edb9a06fcc40
    ##
    ## THIS TEST CAUGHT WHAT THE 37-CLASS DIGEST DID NOT.  That digest
    ## uses default cards and random biases; these rows carry cards that
    ## reach the arms.  It is the better instrument and it earned its
    ## keep here.
    'MosLevel3Hdl': 'da716b9e8c723dcd',
    'MosLevel3Hdl-off': 'ca00b7ee7fae7500',
    'MosLevel3PmosHdl': '7b8575aaf29e3367',
}


def _digest(b):
    return hashlib.sha256(b).hexdigest()[:16]


def _sweep_points(n):
    rng = np.random.RandomState(11)
    pts = [0.7 * rng.uniform(-1, 1, n) for _ in range(30)]
    for v in (0.0, 1.0, -1.0, 5.0, -5.0, 100.0, -100.0, 1e30, -1e30):
        pts.append(np.full(n, v))
    return pts


def _adopter_digest(clsname, card, pins):
    cls = getattr(eh, clsname)
    declared = [p.name for p in cls.instparams]
    el = cls(*[Node(p) for p in pins],
             **{k: v for k, v in card.items() if k in declared})
    el.update_iparv()
    n = len(hdl.x_layout(el))
    out = []
    with np.errstate(all='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for x in _sweep_points(n):
            out.append(_digest(b''.join(
                np.ascontiguousarray(np.asarray(getattr(el, m)(x), float)
                                     ).tobytes()
                for m in ('i', 'G', 'q', 'C'))))
    return _digest(''.join(out).encode())


@pytest.mark.parametrize('label,clsname,card,pins', _ADOPTERS,
                         ids=[a[0] for a in _ADOPTERS])
def test_the_adopted_sites_did_not_move_a_bit(label, clsname, card, pins):
    assert _adopter_digest(clsname, card, pins) == _RECORDED[label]


def test_the_adoption_record_can_fail():
    """Mutation check on the record itself: one parameter moved, and
    every digest moves with it."""
    assert _adopter_digest('GummelPoonNpnHdl', dict(_BJT, vtf=0.5),
                           ('c', 'b', 'e')) != _RECORDED['GummelPoonNpnHdl']


def test_select_with_a_margin_reproduces_the_hand_clamped_piecewise():
    """The `alpha` shape from MOS level 3 (roadmap sec. 21 adoption).

    Hand-written::

        Piecewise((k / maxc(n, 1.0), n > 0.0), (0.0, True))

    and the `select` form that replaces it::

        select((k / n, n > 0.0), (0.0, True), margin=1.0)

    `margin` exists for exactly this: `n > 0` alone clamps to
    ``maxc(n, 0)``, and an arm that divides by `n` still divides by zero.
    Both forms are compiled and compared over the region where the arm is
    SELECTED, over the boundary, and well outside it -- value and
    derivative.
    """
    import numpy as np
    import sympy
    from pycircuit.circuit import hdl

    n = sympy.Symbol('n', real=True)
    k = 2.0

    def hand():
        return sympy.Piecewise((k / hdl.maxc(n, 1.0), n > 0.0),
                               (sympy.Float(0.0), True))

    def viaselect():
        return hdl.select((k / n, n > 0.0), (sympy.Float(0.0), True),
                          margin=1.0)

    f_hand = hdl.compile_chain(hand, [n])
    f_sel = hdl.compile_chain(viaselect, [n])
    d_hand = hdl.compile_chain(hand, [n], wrt=[n])
    d_sel = hdl.compile_chain(viaselect, [n], wrt=[n])

    ## selected region (well inside), the boundary, and outside it
    pts = [1e21, 1e16, 1e6, 100.0, 2.0, 1.0, 0.5, 0.0, -1.0, -1e16]
    for v in pts:
        a = float(np.asarray(f_hand(v)).ravel()[0])
        b = float(np.asarray(f_sel(v)).ravel()[0])
        assert a == b, (v, a, b)
        da = float(np.asarray(d_hand(v)).ravel()[0])
        db = float(np.asarray(d_sel(v)).ravel()[0])
        assert da == db, ('derivative', v, da, db)

    ## and the point of the whole thing: neither form is ever non-finite,
    ## including at n = 0 where the unclamped arm would divide by zero.
    for v in pts:
        assert np.isfinite(float(np.asarray(f_sel(v)).ravel()[0])), v
        assert np.isfinite(float(np.asarray(d_sel(v)).ravel()[0])), v


def test_the_mos3_alpha_arm_is_reached_only_with_nsub_set():
    """⚠ A digest over DEFAULT cards cannot see this conversion.

    `nsm3 = maxc(nsub, 0) * 1e6`, and `nsub` defaults to 0 -- so at
    default parameters the `nsm3 > 0` arm is never selected and any
    margin whatsoever gives identical numbers.  Verified while making
    this change: a deliberately wrong `margin=1e-30` digested BIT-
    IDENTICALLY across all 37 library classes over 30 random biases.

    So the card below is not decoration; it is the only reason this
    change is tested at all.  See `validation-design`: a bias no sweep
    visits tests nothing.
    """
    import numpy as np
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.circuit import defaultepar

    default = eh.MosLevel3Hdl('d', 'g', 's', 'b')
    default.update_iparv()
    assert float(default.iparv.nsub) == 0.0, 'default nsub moved'

    card = dict(nsub=1e16, vto=0.7, kp=2e-5, tox=2e-8, ld=1e-7, u0=600.0,
                xj=2e-7, delta=0.5, eta=0.1, vmax=1e5, nfs=1e11)
    el = eh.MosLevel3Hdl('d', 'g', 's', 'b', **card)
    el.update_iparv()
    rng = np.random.default_rng(11)
    for _ in range(20):
        x = rng.uniform(-2.0, 3.0, el.n)
        for meth in ('i', 'G', 'q', 'C'):
            v = np.asarray(getattr(el, meth)(x, defaultepar), dtype=float)
            assert np.all(np.isfinite(v)), (meth, x)
