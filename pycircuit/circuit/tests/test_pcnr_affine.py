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

## ⚠ `NOT_LIFTABLE` USED TO LIVE HERE and is gone because it emptied.
##
## It listed the models refused for a reason that was not the
## constant-conductance class -- Gummel-Poon's bias-dependent `rbb`, the
## self-heating thermal node, the LED's branch unknown, the photodiode's
## optical port -- and its test asserted none of them was wrongly
## claimed as affine.  Every one has since been declared as a probe or
## had its gate relaxed (roadmap sec. 42-45), so the list is empty and
## an empty `parametrize` SKIPS rather than fails: a test that looks
## green and measures nothing.
##
## The refusal paths it guarded are still covered, deliberately:
##   * `test_a_nonlinear_dependence_is_still_refused` -- the affine
##     split itself, on synthetic expressions;
##   * `test_a_probe_dependent_coefficient_is_refused` -- bilinear;
##   * `test_a_generated_state_is_still_refused` -- device memory, the
##     half of the old branch/state gate that survived.


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
## 4.  Gummel-Poon: solved by a probe, not by the affine remainder.
## ======================================================================

#: A real BJT card.  `fold_card` refuses bare defaults, and rightly --
#: "defaults are not a physical card".
GP_CARD = dict(IS=2.5e-14, bf=180.0, br=3.0, vaf=60.0, ne=1.6, nc=2.0,
               ise=4e-15, isc=8e-15, rb=25.0, re=1.0, rc=8.0)


def test_gummel_poons_base_resistance_is_a_probe_not_a_remainder():
    """The base resistance was the case the affine remainder CANNOT take.

    SPICE's base resistance is

        rbmx = var(Piecewise((rb, rbm < 0), (rbm, True)))
        rbb  = var((rbmx + (rb - rbmx) / qb) / area)

    and `qb` is the base-charge factor, which moves with bias.  With
    `rbm` unset or equal to `rb` the `(rb - rbmx)` factor is an exact
    zero and `rbb` is constant -- but on a real high-injection card
    (`rbm < rb`) it is not, and `g(v_lim) * (x_a - x_b)` is BILINEAR:
    a conductance that changes every iteration, which the ordinary MNA
    assembly cannot carry because it is built before `v_lim` is known.

    The fix is not more solver machinery.  Declaring the branch itself
    as a `limit_identity` probe puts the whole current behind declared
    unknowns, so it qualifies on ANY card -- folded or not, whatever
    `rbm` says.  One extra PCNR unknown per device buys it.
    """
    el = eh.GummelPoonNpnHdl(Node('c'), Node('b'), Node('e'), **GP_CARD)
    el.update_iparv()
    probes = getattr(el, 'pcnr_probes', None)
    assert probes, type(el)._hdl_info.get('pcnr_refusal')
    kinds = [k for _p, k in probes]
    assert kinds.count('id') == 1, kinds
    assert kinds.count('pnj') == 2, kinds
    ## The two mechanisms COMPOSE, and this is the interesting part.
    ## `rbb` is bias-dependent and needs the probe; `rc` and `re` are
    ## plain constants and are taken by the affine remainder.  So this
    ## device carries an identity probe AND a conductance block, each
    ## handling the resistance the other cannot.
    assert getattr(type(el), 'pcnr_lin_G', None) is not None, \
        'the constant rc/re remainder should still be lifted'


#: `rbm` decides whether `rbb` is constant, and once the branch is a
#: probe that no longer decides ELIGIBILITY -- every card qualifies.
#: Kept parametrized because the arithmetic is still the interesting
#: part and a regression would most likely show up on one card only.
RBM_CASES = [
    ({}, 'unset: SPICE reads rbm < 0 as "not given", so rbmx = rb'),
    (dict(rbm=25.0), 'rbm == rb: the same exact cancellation'),
    (dict(rbm=10.0), 'a real high-injection card: rb - rbm = 15, rbb '
                     'is genuinely bias-dependent'),
    (dict(rbm=0.0), 'rb - rbm = 25'),
]


@pytest.mark.parametrize('extra,why', RBM_CASES,
                         ids=['unset', 'rbm_eq_rb', 'rbm_10', 'rbm_0'])
def test_every_rbm_card_qualifies_now(extra, why):
    """Including the ones the affine remainder had to refuse."""
    el = eh.GummelPoonNpnHdl(Node('c'), Node('b'), Node('e'),
                             **dict(GP_CARD, **extra))
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), '%s -- %s' % (
        why, type(el)._hdl_info.get('pcnr_refusal'))


@pytest.mark.parametrize('extra', [{}, dict(rbm=10.0)],
                         ids=['unset', 'high_injection'])
def test_a_gummel_poon_stage_reaches_the_same_operating_point(extra):
    """The gate, on the model this whole exercise started from.

    `rbm=10` is the case that motivated it: a real high-injection card,
    which the affine remainder could not take at all.
    """
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    card = dict(GP_CARD, **extra)

    def build():
        c = SubCircuit()
        c.add_node('c')
        c.add_node('b')
        c['vcc'] = VS('c', gnd, v=5.0)
        c['vb'] = VS('b', gnd, v=0.78)
        c['q1'] = eh.GummelPoonNpnHdl('c', 'b', gnd, **card)
        return c

    assert len(P.pcnr_devices(build())) == 1, 'the BJT is not participating'
    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)
    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, (
        'Gummel-Poon under PCNR reached a different operating point, '
        'relative %.3e\n  ordinary %s\n  pcnr     %s' % (rel, x_ref, x_lim))


def test_the_identity_probe_did_not_move_a_single_number():
    """`limit_identity` must be numerically INERT, or this is not a
    free win but a silent model change.

    Measured when the probe was added: bit-identical `i`/`G`/`q`/`C`
    over 25 random biases for `GummelPoonNpnHdl`, `GummelPoonPnpHdl`
    and `GummelPoonNpnThermalHdl` -- max relative difference exactly
    0.000e+00, not merely small.

    This holds the property going forward: a device with an identity
    probe must evaluate exactly as the same device evaluated at the
    node voltages, because an identity probe applies no law.
    """
    import numpy as np
    from pycircuit.circuit.circuit import defaultepar

    el = eh.GummelPoonNpnHdl(*[Node('n%d' % i) for i in range(3)], **GP_CARD)
    el.update_iparv()
    probes = [p for p, k in el.pcnr_probes if k == 'id']
    assert len(probes) == 1
    (la, lb), = probes

    rng = np.random.default_rng(11)
    for _ in range(20):
        x = rng.uniform(-0.3, 0.9, el.n)
        ## the identity probe's value IS the branch potential: nothing
        ## between the unknown and the node voltages.
        v_branch = float(x[la]) - float(x[lb])
        got = np.asarray(el.pcnr_i(np.array([0.0, 0.0, v_branch]),
                                   {p.name: getattr(el.iparv, p.name)
                                    for p in el.instparams},
                                   defaultepar, el.toolkit), dtype=float)
        assert np.all(np.isfinite(got)), got


## ======================================================================
## 5.  Self-heating and the series-resistance branch.
## ======================================================================

def test_a_self_heating_device_takes_pcnr():
    """`dT` is the deepest dependence there is, and it is now a probe.

    `SelfHeating` sets `dT = var(b.V)` -- a raw branch voltage -- and
    `dT` then sets `Tj`, which enters every current through
    `exp(v/(n*k*Tj/q))`.  That made self-heating disqualifying: not for
    a parasitic on the side, but for the temperature the whole model is
    evaluated at.

    Declared as a `limit_identity` probe, the chain sits behind declared
    unknowns and the device qualifies.
    """
    el = eh.GummelPoonNpnThermalHdl(
        *[Node('n%d' % i) for i in range(5)], rth=100.0, rb=25.0)
    el.update_iparv()
    probes = getattr(el, 'pcnr_probes', None)
    assert probes, type(el)._hdl_info.get('pcnr_refusal')
    assert [k for _p, k in probes].count('id') >= 1


def test_the_isothermal_limit_takes_pcnr_too():
    """⚠ I asserted the opposite here an hour ago, from reasoning.

    This test was written as `..._is_still_refused_and_should_be`, and
    it said the refusal at `rth = 0` was "correct rather than
    incidental: it is the isothermal limit, where a thermal unknown
    means nothing".  That is a fine argument and it was not a
    measurement.  Relaxing the branch-unknown gate (roadmap sec. 44) and
    actually solving the circuit gives agreement with the unlimited path
    to **9.025e-17**.

    The thermal branch does collapse to a zero-volt source, and the
    resulting probe is pinned to zero -- so it is one unknown that
    cannot move.  That makes it USELESS, which is what I reasoned, and
    not HARMFUL, which is what I asserted.  The two are different and
    only one of them justifies a refusal.
    """
    el = eh.GummelPoonNpnThermalHdl(*[Node('m%d' % i) for i in range(5)],
                                    rb=25.0)
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), \
        type(el)._hdl_info.get('pcnr_refusal')


def test_a_generated_state_is_still_refused():
    """What the relaxation did NOT touch, and why.

    `vbranches` was relaxed; `states` was not.  An `idt`/`idtmod`/
    `laplace` state is device MEMORY -- the device's output depends on
    its own history, not only on the present unknowns -- and nothing in
    PCNR has been measured against one.  This pins the distinction so
    the next relaxation has to argue for it separately.
    """
    el = eh.IdtHdl(*[Node('k%d' % i) for i in range(len(eh.IdtHdl.terminals))])
    el.update_iparv()
    assert not getattr(el, 'pcnr_probes', None)
    assert 'a generated state' in (
        type(el)._hdl_info.get('pcnr_refusal') or '')


def test_a_device_carrying_a_branch_unknown_reaches_the_same_point():
    """The gate for the relaxation itself.

    `LedHdl` carries a V-contributed branch outright -- not a collapsed
    one -- and is the model the relaxation actually buys.
    """
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    def build():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('op')
        c['v1'] = VS('a', gnd, v=2.0)
        c['d1'] = eh.LedHdl('a', gnd, 'op', gnd)
        return c

    assert len(P.pcnr_devices(build())) == 1
    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)
    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, 'LED under PCNR: relative %.3e' % rel


@pytest.mark.parametrize('card', [
    dict(rs=2.0), dict(rs=2.0, cjo=1e-12, tt=1e-9), dict(),
], ids=['rs', 'rs_with_charge', 'plain'])
def test_the_spice_diode_takes_pcnr_with_a_series_resistance(card):
    """`rs` used to disqualify the SPICE diode outright.

    Not because of the resistor -- because `pdiss` is
    `Branch(a, c).V * idio`, the TERMINAL voltage, which spans both the
    series resistance and the junction.  With the series branch
    declared, the spanning tree resolves that terminal voltage into
    probe symbols.
    """
    el = eh.DiodeSpiceHdl(Node('a'), Node('c'), **card)
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), \
        type(el)._hdl_info.get('pcnr_refusal')


@pytest.mark.parametrize('name,pins,card', [
    ('DiodeSpiceHdl', 2, dict(rs=2.0)),
    ('DiodeSpiceThermalHdl', 4, dict(rth=200.0, rs=1.0)),
    ('GummelPoonNpnThermalHdl', 5, dict(rth=100.0, rb=25.0)),
], ids=['diode_rs', 'diode_thermal', 'bjt_thermal'])
def test_the_new_probes_do_not_move_the_operating_point(name, pins, card):
    """The gate.  An identity probe carries no law, so the solution it
    reaches must be the one the unlimited path reaches."""
    import numpy as np
    from pycircuit.circuit.elements import SubCircuit, VS
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC

    cls = getattr(eh, name)

    def build():
        c = SubCircuit()
        c.add_node('a')
        if pins >= 4:
            c.add_node('tj')     # only when the device has a thermal port
        c['v1'] = VS('a', gnd, v=0.8 if pins < 5 else 5.0)
        wiring = {2: ['a', gnd],
                  4: ['a', gnd, 'tj', gnd],
                  5: ['a', 'b', gnd, 'tj', gnd]}[pins]
        if pins == 5:
            c.add_node('b')
            c['vb'] = VS('b', gnd, v=0.78)
        c['d1'] = cls(*wiring, **card)
        return c

    assert len(P.pcnr_devices(build())) == 1, 'device is not participating'
    x_ref = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)
    den = max(1e-30, float(np.max(np.abs(x_ref))))
    rel = float(np.max(np.abs(x_lim - x_ref))) / den
    assert rel < 1e-9, '%s: relative %.3e' % (name, rel)


## ======================================================================
## 6.  The photodiode, and the gate that nearly hid the point.
## ======================================================================

SOLAR = dict(IS=1e-9, n=1.3, rs=0.1, rsh=1000.0, resp=0.035)


def _cell(nsun, rload):
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit import gnd
    c = SubCircuit()
    c.add_node('o')
    c.add_node('sun')
    c['vs'] = VS('sun', gnd, v=nsun)
    c['d1'] = eh.PhotodiodeHdl('o', gnd, 'sun', gnd, **SOLAR)
    c['rl'] = R('o', gnd, r=rload)
    return c


def _string(n, rload):
    from pycircuit.circuit.elements import SubCircuit, VS, R
    from pycircuit.circuit import gnd
    c = SubCircuit()
    c.add_node('sun')
    c['vs'] = VS('sun', gnd, v=1.0)
    names = [gnd] + ['n%d' % k for k in range(1, n + 1)]
    for nm in names[1:]:
        c.add_node(nm)
    c['rl'] = R(names[-1], gnd, r=rload)
    for k in range(n):
        c['d%d' % k] = eh.PhotodiodeHdl(names[k + 1], names[k], 'sun', gnd,
                                        **SOLAR)
    return c


def test_the_photodiode_takes_pcnr():
    """The optical port is a declared probe, so the device qualifies.

    ⚠ It is also the one declaration that is a boundary CONDITION
    rather than a device quantity.  Gummel-Poon's base resistance,
    `SelfHeating`'s junction temperature and the diode's series branch
    are all things the device HAS; `V(optp,optn)` is an input someone
    drives.  Recorded because "declare it" stops discriminating if it is
    the answer to every refusal.
    """
    el = eh.PhotodiodeHdl(Node('a'), Node('c'), Node('p'), Node('n'))
    el.update_iparv()
    assert getattr(el, 'pcnr_probes', None), \
        type(el)._hdl_info.get('pcnr_refusal')
    assert [k for _p, k in el.pcnr_probes].count('id') >= 1


@pytest.mark.parametrize('build,label', [
    (lambda: _cell(1.0, 1e5), '1 sun into 100k, near Voc'),
    (lambda: _cell(50.0, 1e7), '50 suns into 10M'),
    (lambda: _string(4, 1e5), 'a 4-cell series string'),
], ids=['near_voc', 'extreme', 'string'])
def test_pcnr_solves_a_solar_cell_more_accurately(build, label):
    """⚠ THE GATE THAT NEARLY HID THE POINT.

    Every other model on this branch is gated on AGREEING with the
    ordinary path to 1e-9.  On these circuits that gate FAILS -- the two
    differ by about 3e-9 -- and the first reading was "PCNR is wrong
    here".

    It is the other way round.  Measured residuals:

        1 sun, 100k    ordinary 2.933e-09    PCNR 5.424e-16
        50 sun, 10M    ordinary 3.113e-09    PCNR 6.545e-15
        4-cell string  ordinary 4.231e-09    PCNR 7.932e-16
        12-cell string ordinary 3.007e-09    PCNR 5.940e-15

    -- five to six orders of magnitude, on a Jacobian conditioned at
    41, so the residual maps almost one-to-one onto voltage error.  The
    ordinary path stops at its tolerance; PCNR drives the residual to
    machine precision.  The `disagreement` was the REFERENCE being
    loose.

    So this device is gated on its RESIDUAL rather than on agreement,
    which is the honest comparison when neither side is known-correct
    in advance.  It is also the convergence benefit I reported as
    absent when investigating whether to declare this probe at all: I
    checked that both converged and that they agreed, and never asked
    which one was right.
    """
    from pycircuit.circuit import gnd, pcnr as P
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.circuit import defaultepar
    import numpy as np

    probe = build()
    ref = probe.get_node_index(gnd)

    def residual(x):
        f = (np.asarray(probe.i(x, defaultepar), dtype=float)
             + np.asarray(probe.u(0, analysis='dc'), dtype=float))
        return float(np.max(np.abs(np.delete(f, ref))))

    x_ord = np.asarray(DC(build(), refnode=gnd).solve().x, dtype=float)
    out = P.solve_dc(build(), gnd)
    x_lim = np.asarray(out[0] if isinstance(out, tuple) else out, dtype=float)

    r_ord, r_lim = residual(x_ord), residual(x_lim)
    assert r_lim <= max(r_ord, 1e-14), (
        '%s: PCNR residual %.3e is worse than the ordinary path\'s %.3e'
        % (label, r_lim, r_ord))


## ======================================================================
## 7.  Nothing that needs limiting may be refused.
## ======================================================================

#: Branch volts to probe.  The top is deliberately outside any operating
#: point: the question is what Newton sees when it OVERSHOOTS.
_VOLTS = (1.0, 10.0)

#: `max|G|(10 V) / max|G|(1 V)` above which a device is taken to need
#: limiting.  Measured, the two populations sit at **1.95e+87** and
#: **7.92e+03** -- eighty-four orders of magnitude apart -- so this
#: threshold is not a tuned parameter.  1e12 is far above a plausible
#: high-order polynomial (a 10th-order term gives 1e10 across that
#: decade) and far below any of the exponentials.
_NEEDS_LIMITING = 1e12


def _max_g_growth(cls):
    """``max|G|(10 V) / max|G|(1 V)``, sweeping EVERY row.

    Sweeping only row 0 makes `GummelPoonNpnHdl` and `MosLevel1Hdl` look
    flat, because their exponential lives on an internal junction.  That
    nearly went into a report as data.

    Mirrors `benchmarks/limiting_need.py`, which is the human-readable
    version with the full table and the two criteria that do NOT work.
    Deliberately reimplemented rather than imported: no test in this
    tree imports from `benchmarks/`, and a guard should not acquire a
    dependency on a report.
    """
    import numpy as np
    from pycircuit.circuit.circuit import defaultepar

    el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))])
    el.update_iparv()
    out = []
    for v in _VOLTS:
        best = 0.0
        for row in range(el.n):
            x = np.zeros(el.n)
            x[row] = v
            try:
                g = np.asarray(el.G(x, defaultepar), dtype=float)
            except Exception:
                continue
            m = float(np.max(np.abs(g)))
            if np.isfinite(m):
                best = max(best, m)
        out.append(best)
    if out[0] <= 0.0:
        return None
    return out[1] / out[0]


def _library():
    from pycircuit.circuit.hdl import Behavioural
    return sorted(n for n, c in vars(eh).items()
                  if isinstance(c, type) and issubclass(c, Behavioural)
                  and c is not Behavioural and '_hdl_info' in c.__dict__)


def test_no_model_that_needs_limiting_is_refused_by_pcnr():
    """The property the whole PCNR arc (sec. 40-45) was aiming at.

    Limiting exists because Newton steps into a region where the
    linearisation is worthless, and that happens when the CONDUCTANCE
    grows super-polynomially with a branch voltage the solver controls.
    Not because a device is nonlinear, and not because its current
    overflows -- `DiodeSpiceHdl` never overflows at any voltage, because
    `expl` bounds it, and it needs limiting all the same.

    Measured, the two populations are eighty-four orders of magnitude
    apart, so this asserts a fact rather than a tuned bound.

    ⚠ One direction only.  A LOW growth ratio is not proof a device is
    safe: `EkvNmosHdl` reads bounded here because a probe from zero bias
    never reaches its exponential region, and it does need limiting.
    This test therefore says "everything that visibly needs limiting is
    eligible", which is the half that can fail usefully.
    """
    offenders = []
    for name in _library():
        cls = getattr(eh, name)
        try:
            ratio = _max_g_growth(cls)
        except Exception:
            continue
        if ratio is None or ratio <= _NEEDS_LIMITING:
            continue
        el = cls(*[Node('m%d' % i) for i in range(len(cls.terminals))])
        el.update_iparv()
        if not (getattr(el, 'pcnr_probes', None)
                or getattr(el, 'pcnr_junctions', ())):
            offenders.append((name, ratio,
                              cls._hdl_info.get('pcnr_refusal') or ''))
    assert not offenders, (
        'model(s) whose conductance grows super-polynomially and which '
        'PCNR refuses:\n' + '\n'.join(
            '  %s  G(10)/G(1) = %.2e  -- %s' % o for o in offenders))


def test_the_growth_probe_can_tell_the_two_populations_apart():
    """The bite-check for the test above.

    A guard that has never been shown to separate anything is not
    evidence.  This asserts the gap is real and enormous -- an
    exponential device against a bounded one -- so a future change that
    flattened the probe (sweeping one row again, say) would fail here
    rather than silently making the guard vacuous.
    """
    diode = _max_g_growth(eh.DiodeSpiceHdl)
    tanh_block = _max_g_growth(eh.ComparatorHdl)
    assert diode > 1e60, 'the exponential device no longer reads as one: %r' % diode
    assert tanh_block is not None and tanh_block < 1e3, tanh_block
    assert diode / tanh_block > 1e50, (diode, tanh_block)
