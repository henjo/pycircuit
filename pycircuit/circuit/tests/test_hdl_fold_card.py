# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`hdl.fold_card` -- card parameters as compile-time constants.

A SPICE `.model` card is fixed for every device that references it,
while geometry is set per instance.  Telling the compiler which is which
lets sympy fold every card-only subexpression as the tree is built.

Measured on PSP with a real IHP sg13g2 card (roadmap sec. 30,
`benchmarks/card_constant_folding.py`): generated ``G`` 2675 -> 1971
lines, its ``minc``/``maxc``/``_step``/``_rdiv`` calls 11 565 -> 7 251,
evaluation **17.73 -> 11.55 ms (1.54x)** on numpy and **83.5 -> 56.5 us
(1.48x)** on the C backend, which also compiles in half the time.

⚠ **It is not bit-identical, and cannot be**: folding reassociates
arithmetic.  On PSP's drain current, banded by magnitude, the difference
is 8.5e-16 relative in strong inversion and 2.7e-15 in the 1e-9..1e-6 A
window the model is validated in -- machine precision where it counts --
while at random biases on all six nodes it reaches 1.6e-6.  A model
validated against a reference must be re-validated with its card folded.

The tests below use cheap models; the PSP numbers live in the benchmark.
"""

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd, hdl, elements_hdl as eh
from pycircuit.circuit.circuit import defaultepar

cm.default_toolkit = numeric


def _src(cls, which='G'):
    return cls._hdl_info['funcs'][which].__dict__.get('_src', '')


def _prims(src):
    return sum(src.count(k) for k in
               ('minc(', 'maxc(', '_step(', '_rdiv(', '_recip2('))


CARD = dict(IS=2.5e-14, bf=180.0, br=3.0, vaf=60.0, ne=1.6, nc=2.0,
            ise=4e-15, isc=8e-15, rb=25.0, re=1.0, rc=8.0)


def test_a_folded_class_gives_the_same_answers():
    base = eh.GummelPoonNpnHdl
    fold = hdl.fold_card(base, instance=(), **CARD)
    a = base('c', 'b', 'e', **CARD); a.update_iparv()
    b = fold('c', 'b', 'e');         b.update_iparv()
    rng = np.random.default_rng(5)
    worst = 0.0
    for _ in range(30):
        x = rng.uniform(-0.2, 0.9, a.n)
        for meth in ('i', 'G', 'q', 'C'):
            va = np.asarray(getattr(a, meth)(x, defaultepar), float)
            vb = np.asarray(getattr(b, meth)(x, defaultepar), float)
            m = np.abs(va) > 1e-25
            if m.any():
                worst = max(worst, float(np.max(np.abs(va[m] - vb[m])
                                                / np.abs(va[m]))))
    ## reassociation only -- not bit-identical, but far inside any
    ## tolerance a circuit solve uses
    assert worst < 1e-9, worst


def test_everything_outside_instance_is_folded():
    """The design error this test exists to prevent.

    Folding only the parameters a card happens to MENTION left 87 of
    PSP's 153 symbolic, almost nothing folded, and measured 0.96x -- a
    slowdown.  A parameter sitting at its default is just as constant as
    one the card sets, so `instance` names what stays variable and
    everything else folds.
    """
    base = eh.GummelPoonNpnHdl
    fold = hdl.fold_card(base, instance=(), **CARD)
    names = [p.name for p in base.instparams]
    assert len(fold._hdl_fold_values) == len(names), (
        'only %d of %d parameters folded' % (len(fold._hdl_fold_values),
                                             len(names)))
    ## ...and folding must actually reach the emitted code.  The
    ## PRIMITIVE count is what the speed depends on (roadmap sec. 29: the
    ## cost is the regularisers), and it is what a real card moves --
    ## the raw line count can sit still while the work drops.
    assert _prims(_src(fold)) < _prims(_src(base))
    ## The BODY must no longer mention the folded names.  (The signature
    ## line still lists them: the compiled function keeps its full
    ## argument list on purpose, so `_args_of` and every caller are
    ## unchanged and the folded arguments are simply unused.)
    import re
    body = '\n'.join(_src(fold).splitlines()[1:])
    base_body = '\n'.join(_src(base).splitlines()[1:])
    for nm in ('bf', 'vaf'):
        assert re.search(r'\b%s\b' % nm, base_body), nm
        assert not re.search(r'\b%s\b' % nm, body), (
            '%s survived folding as a symbol' % nm)


def test_an_instance_parameter_stays_variable():
    base = eh.GummelPoonNpnHdl
    card = dict(CARD)
    card.pop('bf')                      # it is the instance parameter here
    fold = hdl.fold_card(base, instance=('bf',), **card)
    assert 'bf' not in fold._hdl_fold_values
    ## Structural, deliberately.  The obvious test -- evaluate at two
    ## values of `bf` and require the currents to differ -- tests the
    ## chosen BIAS as much as the feature: at a saturated operating point
    ## the terminal currents barely depend on the forward beta, and the
    ## test fails while the code is correct.  What the design guarantees
    ## is that `bf` is still a SYMBOL in the emitted body.
    import re
    body = '\n'.join(_src(fold).splitlines()[1:])
    assert re.search(r'\bbf\b', body), 'bf was folded despite instance='
    assert not re.search(r'\bvaf\b', body), 'vaf should have folded'
    ## and it is still settable, which a folded parameter is not
    a = fold('c', 'b', 'e', bf=50.0)
    a.update_iparv()
    assert float(a.iparv.bf) == 50.0


def test_changing_a_folded_parameter_raises():
    """It cannot be honoured -- the value is in the compiled code, not in
    `iparv` -- so it must not be silently ignored."""
    fold = hdl.fold_card(eh.DiodeHdl, instance=(), IS=2.5e-14)
    ok = fold('p', gnd)
    ok.update_iparv()
    with pytest.raises(ValueError) as e:
        bad = fold('p', gnd, IS=9e-14)
        bad.update_iparv()
    assert 'IS' in str(e.value) and 'fold_card' in str(e.value)


def test_a_folded_instance_reports_the_value_it_evaluates():
    fold = hdl.fold_card(eh.DiodeHdl, instance=(), IS=2.5e-14)
    el = fold('p', gnd)
    el.update_iparv()
    assert float(el.iparv.IS) == 2.5e-14


def test_variants_are_cached_per_card():
    D = eh.DiodeHdl
    a = hdl.fold_card(D, instance=(), IS=2.5e-14)
    assert hdl.fold_card(D, instance=(), IS=2.5e-14) is a
    assert hdl.fold_card(D, instance=(), IS=1e-15) is not a
    ## and a different instance split is a different variant
    assert hdl.fold_card(D, instance=('IS',)) is not a


def test_two_cards_do_not_share_a_compile_cache_entry():
    """Two cards must not collide in the on-disk compile cache.

    The variants have identical `analog` source and identical parameter
    NAMES, so a key that saw only those would make the second card
    silently run the first card's compiled code.

    ⚠ **What separates them is the DEFAULTS, not the `fold=` entry.**
    `fold_card` rewrites each folded parameter's default to its card
    value (so an instance reports what it evaluates), and `_param_decl`
    hashes `p.default` -- so `params=` already discriminates. Verified by
    neutralising the `fold=` entry: this test still passed. It is kept as
    defence in depth, and the assertion below pins the mechanism that is
    actually load-bearing, so that removing the default rewrite cannot
    quietly reintroduce the collision.
    """
    from pycircuit.circuit import _hdl_cache
    D = eh.DiodeHdl
    f1 = hdl.fold_card(D, instance=(), IS=2.5e-14)
    f2 = hdl.fold_card(D, instance=(), IS=1e-15)
    assert _hdl_cache.key_for(f1) != _hdl_cache.key_for(f2)
    ## the load-bearing mechanism: the card reaches `instparams`
    d1 = {p.name: p.default for p in f1.instparams}
    d2 = {p.name: p.default for p in f2.instparams}
    assert d1['IS'] == 2.5e-14 and d2['IS'] == 1e-15, (d1['IS'], d2['IS'])
    ## ...and the numbers really do differ, so the key must
    x = np.array([0.55, 0.0])
    i1 = float(np.asarray(hdl.fold_card(D, instance=(), IS=2.5e-14)
                          ('p', gnd).i(x, defaultepar), float)[0])
    i2 = float(np.asarray(hdl.fold_card(D, instance=(), IS=1e-15)
                          ('p', gnd).i(x, defaultepar), float)[0])
    assert abs(i1 - i2) > 0.1 * abs(i1)


def test_fold_card_refuses_what_it_cannot_fold():
    D = eh.DiodeHdl
    with pytest.raises(TypeError, match='no parameter'):
        hdl.fold_card(D, instance=(), nosuch=1.0)
    with pytest.raises(TypeError, match='real number'):
        hdl.fold_card(D, instance=(), IS='big')
    with pytest.raises(TypeError, match='one or the other'):
        hdl.fold_card(D, instance=('IS',), IS=1e-14)


def test_no_instance_split_and_nothing_to_fold_returns_the_class():
    D = eh.DiodeHdl
    allnames = tuple(p.name for p in D.instparams)
    assert hdl.fold_card(D, instance=allnames) is D


def test_folding_at_bare_defaults_is_refused_with_a_reason():
    """A value that is fine as a SYMBOL can be singular as a NUMBER.

    A model guards `1/p` for the p it expects to be zero; a folded card
    evaluates every such quotient at build time, and a parameter left at
    a 0.0 default makes one ComplexInfinity.  Sympy reports that as
    `KeyError: 'ComplexInfinity'` from inside its printer, which says
    nothing useful -- so `fold_card` names it.  Defaults are not a
    physical card, and this is the test that says so.
    """
    with pytest.raises(ValueError) as e:
        hdl.fold_card(eh.GummelPoonNpnHdl, instance=())
    msg = str(e.value)
    assert 'singular' in msg and '0.0 default' in msg and 'instance' in msg
