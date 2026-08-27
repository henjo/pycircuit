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
def test_a_parasitic_card_now_qualifies(name, card):
    """The lift, from the outside: these were refused, and are not now.

    Each of these is a card a PDK would ship -- a MOSFET with a drain
    resistance, a MESFET with a source resistance.  Before the lift,
    declaring `rd` cost the device PCNR entirely.
    """
    cls = getattr(eh, name)
    el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))], **card)
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), \
        '%s still refused: %s' % (name, _refusal(name, **card)[0])
    assert hasattr(type(el), 'pcnr_lin_G'), \
        '%s qualifies but carries no conductance block -- the remainder '\
        'would be dropped from the assembly' % name


@pytest.mark.parametrize('name,card', LIFTABLE_CARDS,
                         ids=[n for n, _ in LIFTABLE_CARDS])
def test_the_remainder_is_still_recorded_as_affine(name, card):
    """The detection must survive the lift, or a refusal cannot explain
    itself and the next reader has no data to work from."""
    _why, affine = _refusal(name, **card)
    assert affine, name
    for label, (krow, _rest, coeffs) in affine.items():
        assert label.startswith('I('), label
        assert krow is not None and krow >= 0, (label, krow)
        assert coeffs, label


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

@pytest.mark.parametrize('name,card', LIFTABLE_CARDS,
                         ids=[n for n, _ in LIFTABLE_CARDS])
def test_the_lifted_device_reaches_the_same_operating_point(name, card):
    """THE GATE, and the reason it is not `did it converge`.

    A lifted device must land on the same operating point as the
    unlimited path, because a wrong lift changes answers silently.  The
    first version of this lift did exactly that: on

        vd=2.0, vg=1.5, MosLevel1Hdl(rd=5, w=10u, l=1u)

    it put the internal drain node at **1.4989 V** against **1.9989 V**
    -- out by exactly the gate voltage, the drain resistor zeroed out of
    the assembly and never re-stamped -- while the device's terminal
    CURRENT still agreed to seven digits and the solver converged
    without a murmur.  Convergence saw nothing.  This comparison saw it
    immediately.

    Measured now: 1.5e-18, 1.1e-16, 2.8e-17 relative.
    """
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    cls = getattr(eh, name)
    npins = len(cls.terminals)

    def build():
        c = SubCircuit()
        c.add_node('d')
        c.add_node('g')
        c['vd'] = VS('d', gnd, v=2.0)
        c['vg'] = VS('g', gnd, v=1.5)
        c['m1'] = cls(*(['d', 'g', gnd, gnd][:npins]),
                      **dict(card, **({'w': 10e-6, 'l': 1e-6}
                                      if 'w' in [p.name for p in cls.instparams]
                                      else {})))
        return c

    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)

    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, (
        '%s: PCNR reached a different operating point, relative %.3e\n'
        '  ordinary %s\n  pcnr     %s' % (name, rel, x_ref, x_lim))


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

    was = hdl.PCNR_LIFT_AFFINE
    here = cache.key_for(eh.MosLevel1Hdl)
    try:
        hdl.PCNR_LIFT_AFFINE = not was
        flipped = cache.key_for(eh.MosLevel1Hdl)
    finally:
        hdl.PCNR_LIFT_AFFINE = was
    again = cache.key_for(eh.MosLevel1Hdl)
    off, on = here, flipped

    assert on != off, 'the lift flag does not reach the cache key'
    assert again == off, 'the key is not stable with the flag restored'


def test_a_collapsed_variant_does_not_inherit_the_wrong_block():
    """A conductance block belongs to ONE geometry.

    Card-collapsed variants are subclasses, so an attribute set only on
    the classes that have a remainder is INHERITED by those that do not
    -- and the inherited block is sized for the parent's node count.
    Enabling the lift with that hole failed three tests at once with
    `shapes (6,6) and (4,) not aligned`, from a `(6,6)` block dotted
    into a four-node variant's local voltages.

    So `pcnr_lin_G` is set on every class that reaches the vector route,
    `None` included.  This checks the two states really are independent
    rather than one shadowing the other.
    """
    lifted = eh.MosLevel1Hdl(*[Node('n%d' % i) for i in range(4)], rd=5.0)
    lifted.update_iparv()
    plain = eh.MosLevel1Hdl(*[Node('m%d' % i) for i in range(4)])
    plain.update_iparv()

    assert type(lifted) is not type(plain), 'the cards did not fork a variant'
    assert getattr(type(lifted), 'pcnr_lin_G', None) is not None

    ## The un-lifted variant must report NO block -- not the other one's.
    other = getattr(type(plain), 'pcnr_lin_G', None)
    if other is not None:
        import numpy as np
        from pycircuit.circuit.circuit import defaultepar
        from pycircuit.circuit import pcnr as P
        dev = P._device_of('m1', plain, list(range(plain.n)))
        gl = np.asarray(other(dev.params(), defaultepar, plain.toolkit))
        assert gl.shape == (plain.n, plain.n), (
            'inherited a block shaped %s for a %d-node device'
            % (gl.shape, plain.n))


## ======================================================================
## 4.  Gummel-Poon: the cancellation `var()` hides.
## ======================================================================

#: A real BJT card.  `fold_card` refuses bare defaults, and rightly --
#: "defaults are not a physical card".
GP_CARD = dict(IS=2.5e-14, bf=180.0, br=3.0, vaf=60.0, ne=1.6, nc=2.0,
               ise=4e-15, isc=8e-15, rb=25.0, re=1.0, rc=8.0)


def _folded_gp():
    return hdl.fold_card(eh.GummelPoonNpnHdl, instance=(), **GP_CARD)


def test_gummel_poon_is_refused_while_its_base_resistance_is_free():
    """Unfolded, `rbb` really is bias-dependent and must be refused.

    SPICE's base resistance is

        rbmx = var(Piecewise((rb, rbm < 0), (rbm, True)))
        rbb  = var((rbmx + (rb - rbmx) / qb) / area)

    and `qb` is the base-charge factor, which moves with bias.  For a
    general `rbm` that is a conductance that changes every iteration --
    bilinear in `v_lim` and `x` -- and the ordinary assembly, built
    before `v_lim` is known, cannot carry it.
    """
    el = eh.GummelPoonNpnHdl(Node('c'), Node('b'), Node('e'), **GP_CARD)
    el.update_iparv()
    assert not getattr(el, 'pcnr_probes', None)
    assert 'no $limit probe limits' in (
        type(el)._hdl_info.get('pcnr_refusal') or '')


def test_a_folded_gummel_poon_card_qualifies():
    """Folded, `rbm` is known -- and SPICE's "unset" makes `rbb` CONSTANT.

    `rbm < 0` is SPICE's sentinel for "not given", and it means
    `rbmx = rb`, so `(rb - rbmx)` is an exact zero and the whole `qb`
    term vanishes: `rbb = rb / area`.

    ⚠ sympy cannot see that while `rbmx` is a held `var()` -- it reads
    `25.0 - _v67_rbmx`, which is not zero.  A probe test that follows
    chain REFERENCES therefore over-refuses; one that RESOLVES the chain
    and asks what survives gets it right.  That is the difference
    between this test and the one above, and both must hold: folding
    changes the physics of the coefficient, not merely its spelling.
    """
    el = _folded_gp()(Node('c'), Node('b'), Node('e'))
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), (
        type(el)._hdl_info.get('pcnr_refusal') or '')
    assert getattr(type(el), 'pcnr_lin_G', None) is not None


def test_a_folded_gummel_poon_reaches_the_same_operating_point():
    """The gate, for the model this whole exercise started from."""
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    fold = _folded_gp()

    def build():
        c = SubCircuit()
        c.add_node('c')
        c.add_node('b')
        c['vcc'] = VS('c', gnd, v=5.0)
        c['vb'] = VS('b', gnd, v=0.75)
        c['q1'] = fold('c', 'b', gnd)
        return c

    assert len(P.pcnr_devices(build())) == 1, 'the BJT is not participating'
    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)
    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, (
        'folded Gummel-Poon under PCNR reached a different operating '
        'point, relative %.3e\n  ordinary %s\n  pcnr     %s'
        % (rel, x_ref, x_lim))
