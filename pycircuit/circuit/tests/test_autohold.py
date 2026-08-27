# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`_autohold`: the regularisers hold their own arguments (roadmap sec. 36).

Roadmap S4 proposed that the compiler hold **every** non-atomic
sub-expression, so a model author never has to remember `var()`.
Measured, that is the wrong shape:

* **52%** of the 342 holds `psp_kernel.py` writes by hand are ≤3 ops --
  naming, not structure -- so "hold everything" mostly holds nothing
  that needed holding;
* but **434 of 1184** calls into the regularisers receive a *bare*
  non-atomic argument, `safe_div` one of **44 ops**.

That is where an unheld expression gets substituted into something else
and evaluated in an order the author never wrote -- what `_sigma_body`'s
docstring describes, and what S4 exists to prevent.  So the regularisers
hold their own arguments and the author does not have to.

Measured on PSP with a real IHP card, interleaved:

===================  ==========  ==========  ==========
..                   build       G source    G evaluate
===================  ==========  ==========  ==========
without `_autohold`    104-117 s   2675 lines   22.6/24.8 ms
with                   108-110 s   2860 lines   18.6/20.2 ms
===================  ==========  ==========  ==========

**Evaluation is ~1.22x FASTER** -- a held sub-expression is computed once
instead of being substituted and recomputed -- and the build does not get
slower despite 185 more definitions.
"""

import warnings

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import hdl, gnd

cm.default_toolkit = numeric

X = sympy.Symbol('x', real=True)


## ---------------------------------------------------------------------
## 1. When it declines to hold -- each condition is load-bearing
## ---------------------------------------------------------------------

def test_outside_a_compile_it_is_the_identity():
    """`var()` RAISES outside `analog()`, and every regulariser is also
    called directly -- from tests, from `apply_limit`, from the kernel's
    own checks.  Without this guard, importing would break."""
    assert not hdl._VAR_STACK, 'this test assumes no compile in progress'
    e = X * X + X + sympy.sin(X)
    assert hdl._autohold(e, 'q') is e
    ## ...and the regularisers themselves still work out here
    assert hdl.safe_div(1.0, X * X + X) is not None
    assert hdl.hypsmooth(X * X + X, 1e-9) is not None


def test_an_atom_is_never_held():
    assert hdl._autohold(X, 'q') is X
    assert hdl._autohold(sympy.Float(2.5), 'q') == 2.5


def test_a_plain_number_is_never_held():
    """A float or int is already a leaf; `count_ops` on one is not even
    meaningful."""
    for v in (1.0, 0, -3, 1e-30):
        assert hdl._autohold(v, 'q') is v


def test_a_tiny_expression_is_not_held():
    """Naming a two-op expression adds a chain line and buys nothing.
    52% of the holds a real model writes are this small, which is why
    "hold everything" is the wrong rule."""
    small = X + 1.0
    assert sympy.count_ops(small) < hdl.AUTOHOLD_MIN_OPS
    assert hdl._autohold(small, 'q') is small


## ---------------------------------------------------------------------
## 2. When it does hold -- inside a compile
## ---------------------------------------------------------------------

def _chain_of(builder):
    """Compile `builder` through the DSL and return its definitions."""
    seen = {}

    def wrapped():
        e = builder()
        seen['defs'] = list(hdl._VAR_STACK[-1])
        return e

    hdl.compile_chain(wrapped, [X])
    return seen['defs']


BIG = X ** 4 + 3 * X ** 3 + 2 * X ** 2 + X


def test_it_declines_until_the_body_has_held_something_itself():
    """⚠⚠ **The condition that stops this being a regression.**

    `var()` is what makes a model CHAINED, and a chained model has no
    `eval_i_pure` -- so it loses `solve_batched` (20.7x at 512 lanes) and
    the JAX path with it.  Sec. 22 measured the two fast paths as
    DISJOINT.  Auto-holding unconditionally would move any eager model
    that divides or exponentiates out of the batchable set.

    Found by `KernelAllEager` in `test_hdl_kernel.py`, which lost
    `eval_i_pure` outright -- and NOT by the four gates this change was
    planned against, nor by a flip-check across the 37 library classes,
    because the 15 batchable ones are passives that never call a
    regulariser.

    So a body that has held nothing yet is left alone.
    """
    held = [str(s) for s, _ in _chain_of(lambda: hdl.expl(BIG))]
    assert held == [], held


def test_it_holds_once_the_body_is_already_chained():
    """After one `var()` of the body's own, the model is chained however
    this behaves, so holding more is free."""
    def body():
        hdl.var(BIG * 2.0, 'mine')          # the body commits first
        return hdl.expl(BIG)
    held = [str(s) for s, _ in _chain_of(body)]
    assert any('expl_x' in h for h in held), held


def test_each_regulariser_holds_its_argument():
    for fn, tag in ((lambda: hdl.expl(BIG), 'expl_x'),
                    (lambda: hdl.hypsmooth(BIG, 1e-9), 'hyp_x'),
                    (lambda: hdl.safe_ln(BIG), 'sln_x')):
        def body(fn=fn):
            hdl.var(BIG * 2.0, 'mine')
            return fn()
        held = [str(s) for s, _ in _chain_of(body)]
        assert any(tag in h for h in held), (tag, held)


def test_safe_div_holds_too_and_the_exclusion_was_a_misreading():
    """It is the site where the value is largest -- the biggest bare
    arguments in the tree, 44 ops, 100 of the 434.

    It was EXCLUDED for a day (roadmap sec. 36) because holding here made
    `TestPspBitIdentity` report `{'equal': 989, 'zero-sign': 16,
    'value': 0}` on `G`, read as "16 signed zeros".  A signed zero is
    observable -- `1/-0.0` is `-inf` -- so refusing to erode the C
    backend's bitwise guarantee was the right instinct on that reading.

    The reading was wrong (sec. 37).  Traced, they are **240 NaN sign
    bits** at 16 bias points, all at extreme biases where `G` is already
    NaN in both paths.  IEEE-754 does not define that bit.  `_compare`
    had folded "signed zero" and "NaN sign" into one label; it now
    separates them, and the strict assertion still bans the first.

    PSP `G`: 17.680 ms unheld, 15.784 without `safe_div`, **14.121 with
    it**.
    """
    def body():
        hdl.var(BIG * 2.0, 'mine')
        return hdl.safe_div(1.0, BIG)
    held = [str(s) for s, _ in _chain_of(body)]
    assert any('sd_' in h for h in held), held


def test_holding_does_not_change_the_value():
    """The point is evaluation ORDER, not arithmetic.  A held chain and a
    flattened expression must agree wherever both are finite."""
    def body():
        hdl.var(BIG * 2.0, 'mine')
        return hdl.expl(BIG)
    f_held = hdl.compile_chain(body, [X])
    was = hdl._autohold
    hdl._autohold = lambda x, name: x
    try:
        f_flat = hdl.compile_chain(body, [X])
    finally:
        hdl._autohold = was
    for v in (-3.0, -0.5, 0.0, 0.25, 1.0, 7.0, 1e3):
        a = float(np.asarray(f_held(v)).ravel()[0])
        b = float(np.asarray(f_flat(v)).ravel()[0])
        assert a == pytest.approx(b, rel=1e-12, abs=1e-300), v


## ---------------------------------------------------------------------
## 3. The library
## ---------------------------------------------------------------------

def test_the_library_still_computes_what_it_did():
    """36 of 37 classes are bit-identical with `_autohold` on and off;
    `DiodeSpiceThermalHdl` differs by **2.08e-15** relative (median
    1.76e-16, one ulp) because holding stops sympy flattening across the
    boundary.  That is reassociation at the last bit, and this pins the
    magnitude so a real change cannot hide behind it.
    """
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.circuit import defaultepar
    el = eh.DiodeSpiceThermalHdl(
        *['n%d' % i for i in range(len(eh.DiodeSpiceThermalHdl.terminals))])
    el.update_iparv()
    rng = np.random.default_rng(7)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for _ in range(20):
            x = rng.uniform(-2.0, 2.0, el.n)
            for meth in ('i', 'G', 'q', 'C'):
                v = np.asarray(getattr(el, meth)(x, defaultepar), float)
                assert np.all(np.isfinite(v)), (meth, x)
