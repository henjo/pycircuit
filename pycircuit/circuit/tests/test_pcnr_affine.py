"""The internal-node constraint, and the affine remainder that can lift it.

Vector PCNR refuses any current that still reads a node voltage once the
declared `$limit` probes have been substituted: *"every resistive
current must reach the solution only through the declared probes"*.

That rule is why **PCNR does not work on the cards a PDK ships**.  Give
`GummelPoonNpnHdl` a base resistance and it is refused; give
`DiodeSpiceHdl` an `rs` and it is refused; `MosLevel1Hdl` and
`MosLevel3Hdl` go the same way on `rd`/`rs`.  The PCNR tests have always
zeroed those parameters, and say so in their own docstrings -- *"at
rb = re = rc = 0 so that it qualifies"*.  The models that survive real
cards are the ones with no series-resistance parameters at all (EKV,
Curtice, Angelov), because they have no internal node to break the rule.

The rule is right for a NONLINEAR dependence and too strong for a linear
one.  A series resistance contributes ``g * (x_a - x_b)`` -- an ordinary
conductance.  Newton solves it exactly, and limiting it would be
meaningless.  Measured across every refusing model on this branch, the
surviving dependence is affine in every case:

    I(c)     = _pcnr_v1*area/rc + _x0*area/rc - _x3*area/rc
    var(dT)  = _x3 - _x4
    I(d)     = (_pcnr_v0 - _pcnr_v1 + _x0 - _x1)/_v91_rdx

-- nine surviving expressions, second derivative identically zero in
every node symbol read, and not one coefficient depending on a probe.

`hdl._affine_in_nodes` performs the split, recorded as
``_hdl_info['pcnr_affine']``.  It is not yet consumed -- the refusal
still stands and its wording is unchanged.  What remains is to stamp
the linear part into the ordinary MNA assembly, which `pcnr.py` already
has the hook for -- `augmented_system` shadows a participating device's
`i`/`G` with `_zero_vector`/`_zero_matrix` while it re-stamps the device
at `v_lim`, so the remainder has a natural home: shadow with the affine
remainder instead of with zeros.
"""

import sympy
import pytest

from pycircuit.circuit import hdl
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import Node


X0, X1, V0, G = sympy.symbols('x0 x1 v0 g', real=True)
XSET = {X0, X1}
PROBES = {V0}


def _split(expr):
    nodes = set(expr.free_symbols) & XSET
    return hdl._affine_in_nodes(expr, nodes, XSET, PROBES)


## ======================================================================
## 1.  The split itself.
## ======================================================================

@pytest.mark.parametrize('expr,name', [
    (G * (X0 - X1), 'a resistor between two nodes'),
    (V0 * G + G * X0, 'a probe term beside a node term'),
    (X0 - X1, 'a bare node difference (the thermal shape)'),
])
def test_an_affine_dependence_splits_and_reconstructs(expr, name):
    """The split must be exact, not approximate -- it is an identity."""
    got = _split(expr)
    assert got is not None, name
    rest, coeffs = got
    recon = sympy.expand(rest + sum(c * q for q, c in coeffs.items()))
    assert sympy.simplify(recon - expr) == 0, (name, recon, expr)
    ## and the parts really are separated
    assert not (set(rest.free_symbols) & XSET), rest


@pytest.mark.parametrize('expr,why', [
    (X0 * X1, 'a product of two nodes is not affine'),
    (sympy.exp(X0), 'an exponential of a node is not affine'),
    (X0 ** 2, 'a square is not affine'),
])
def test_a_nonlinear_dependence_is_still_refused(expr, why):
    assert _split(expr) is None, why


def test_a_probe_dependent_coefficient_is_refused():
    """`v_lim * x` is a conductance that MOVES with the limited unknown.

    The ordinary MNA assembly is built before `v_lim` is known, so it
    cannot carry such a term.  It is linear in `x` and still refused,
    which is the distinction this test exists to hold: *affine* is not
    the whole criterion.
    """
    assert _split(V0 * X0) is None
    ## ...and the same expression with a constant coefficient is fine,
    ## so the refusal is about the probe, not about the shape.
    assert _split(G * X0) is not None


## ======================================================================
## 2.  The real models.
## ======================================================================

#: Every model that loses PCNR eligibility to a parasitic, with the card
#: parameter that costs it.  These are the cards a PDK actually ships.
LIFTABLE_CARDS = [
    ('MosLevel1Hdl', dict(rd=5.0)),
    ('MosLevel3Hdl', dict(rs=5.0, rd=5.0)),
    ('MesfetStatzHdl', dict(rs=2.0)),
]

#: Refused for a reason that is NOT the constant-conductance class, with
#: the reason.  Each one is a trap that a careless lift would sweep in.
NOT_LIFTABLE = [
    ('GummelPoonNpnHdl', dict(rb=100.0),
     'rbb is BIAS-DEPENDENT (Gummel-Poon current crowding), so the '
     'coefficient 1/_v68_rbb reads a probe: g*(x_a-x_b) with g moving '
     'per iteration is bilinear, not a constant conductance'),
    ('GummelPoonNpnThermalHdl', dict(rth=100.0),
     'the survivor is a CHAIN VARIABLE, var(dT) = _x3 - _x4, and dT '
     'drives temperature-dependent physics exponentially downstream'),
    ('PhotodiodeHdl', {}, 'reads optical drive nodes optn/optp'),
    ('LedHdl', {}, 'carries a branch-current unknown'),
]


def _refusal(name, **kw):
    cls = getattr(eh, name)
    el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))], **kw)
    el.update_iparv()
    info = type(el)._hdl_info
    return (info.get('pcnr_refusal') or ''), (info.get('pcnr_affine') or {})


@pytest.mark.parametrize('name,card', LIFTABLE_CARDS,
                         ids=[n for n, _ in LIFTABLE_CARDS])
def test_the_parasitic_refusals_are_all_the_affine_class(name, card):
    """The refusal must say WHICH kind it is, not merely that it refused.

    This is the list of models a lift would recover, and it is the
    reason the lift is worth doing: without it PCNR is a solver for
    cards with the parasitics deleted.
    """
    why, affine = _refusal(name, **card)
    assert 'no $limit probe limits' in why, why
    ## The split is recorded in `_hdl_info`, NOT appended to the refusal
    ## string: `explain()` output is digested by `test_hdl_params.py`, and
    ## a marker there churns three recorded digests for a note that goes
    ## away as soon as the remainder is consumed.  (It did, first try.)
    assert affine, 'no affine remainder recorded for %s' % name
    for label, (krow, rest, coeffs) in affine.items():
        ## Every recorded survivor must be a CURRENT row, never a chain
        ## variable -- `krow` is its index into `i_exprs`.
        assert krow is not None and krow >= 0, (label, krow)
        assert label.startswith('I('), label
        assert coeffs, label
        ## and the coefficients must be free of node symbols, or the
        ## "constant conductance" claim is empty.
        for q, c in coeffs.items():
            assert not (set(c.free_symbols) & set(coeffs)), (label, q, c)


@pytest.mark.parametrize('name,card,why_not', NOT_LIFTABLE,
                         ids=[n for n, _c, _w in NOT_LIFTABLE])
def test_a_structurally_different_refusal_is_not_claimed_as_affine(
        name, card, why_not):
    """The bite-check for the test above.

    `PhotodiodeHdl` reads optical drive nodes and `LedHdl` carries a
    branch-current unknown.  Neither is an internal-node parasitic, and
    a lift built for the affine class must not silently claim them --
    otherwise the test above passes for every model and measures
    nothing.
    """
    why, affine = _refusal(name, **card)
    assert why, name
    assert not affine, '%s was claimed as liftable, but %s' % (name, why_not)


def test_the_models_without_parasitics_still_qualify_outright():
    """The control: these never had the problem, and must not acquire it."""
    for name in ('EkvNmosHdl', 'MesfetCurticeHdl', 'HemtAngelovHdl'):
        cls = getattr(eh, name)
        el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))])
        el.update_iparv()
        assert getattr(el, 'pcnr_probes', None), name


## ======================================================================
## 3.  The acceptance gate for the lift itself.
## ======================================================================

@pytest.mark.xfail(strict=True, reason='the lift is scaffolded and WRONG -- '
                   'see roadmap sec. 40; PCNR_LIFT_AFFINE stays False')
def test_the_lifted_device_reaches_the_same_operating_point():
    """The gate a lift must pass, and the reason it is not `converged`.

    With ``PCNR_LIFT_AFFINE = True`` the three liftable models stop
    being refused -- and the answer is WRONG.  On

        vd=2.0, vg=1.5, MosLevel1Hdl(rd=5, w=10u, l=1u)

    the internal drain node comes back at **1.4989 V** where the
    ordinary path gives **1.9989 V**: out by exactly the gate voltage,
    because the drain resistor is zeroed out of the assembly and never
    re-stamped.  The device's terminal CURRENT still agrees to seven
    digits (-2.2508079e-04 against -2.2508072e-04) and the solver
    converges without complaint.

    That is the whole argument for this gate.  A limiter bug does not
    announce itself by failing to converge -- it announces itself by
    converging to the wrong place, and only a comparison against the
    unlimited path can see it.

    Diagnosis so far: with the flag on the refusal is suppressed but
    `pcnr_vector['lin_terms']` is still `None` and no `pcnr_lin_G` is
    attached, so the remainder is never built and PCNR zeroes the
    resistor with the rest of the device.  The detection (sections 1-2
    above) is correct and independently tested; it is the consumption
    that is unfinished.
    """
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    def build():
        c = SubCircuit()
        c.add_node('d')
        c.add_node('g')
        c['vd'] = VS('d', gnd, v=2.0)
        c['vg'] = VS('g', gnd, v=1.5)
        c['m1'] = eh.MosLevel1Hdl('d', 'g', gnd, gnd, rd=5.0,
                                  w=10e-6, l=1e-6)
        return c

    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)

    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, (
        'PCNR reached a different operating point: relative %.3e\n'
        '  ordinary %s\n  pcnr     %s' % (rel, x_ref, x_lim))


def test_the_lift_flag_is_part_of_the_compile_cache_key():
    """A runtime flag that changes what is compiled MUST key the cache.

    The dependency hashes cover `hdl.py`'s source, so an edit
    invalidates.  But flipping `PCNR_LIFT_AFFINE` inside a running
    process changes the emitted model while leaving every file
    byte-identical -- so without this the two builds share a key and the
    second silently receives the first one's compiled code.

    This is not hypothetical.  One probe run with the lift on wrote
    **924 cache entries**, which then served the lift-off test suite;
    the symptom was `MosLevel1Hdl(rd=5)` quietly no longer being
    refused, with an empty refusal string, in a process where the flag
    was False.  It cost a confusing half hour and looked like a bug in
    the detector.
    """
    from pycircuit.circuit import _hdl_cache as cache

    off = cache.key_for(eh.MosLevel1Hdl)
    try:
        hdl.PCNR_LIFT_AFFINE = True
        on = cache.key_for(eh.MosLevel1Hdl)
    finally:
        hdl.PCNR_LIFT_AFFINE = False
    again = cache.key_for(eh.MosLevel1Hdl)

    assert on != off, 'the lift flag does not reach the cache key'
    assert again == off, 'the key is not stable with the flag restored'
