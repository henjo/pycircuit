"""FET Newton limiting on the ordinary path -- roadmap 10.3(a).

`_fetlim` / `_limvds` (SPICE's `DEVfetlim` / `DEVlimvds`), the `limit_fet`
and `limit_vds` declarations that put them on an HDL element, and the
`limit_spec` restructure that made `limit()` dispatch on a limiter KIND
instead of assuming `pnjlim`.

The reference followed is **ngspice-47**,
`src/spicelib/devices/devsup.c`, with the device-level call sequence read
from `src/spicelib/devices/mos1/mos1load.c`.  Every expected value below
was worked out BY HAND from that C, branch by branch, and the derivation
is written next to the number -- not captured from this implementation's
own output, which would only prove it is self-consistent.
"""
import numpy as np
from numpy.testing import assert_allclose
import pytest

import pycircuit.circuit.circuit
from pycircuit.circuit import numeric
from pycircuit.circuit.circuit import SubCircuit, gnd, defaultepar
from pycircuit.circuit.elements import VS
from pycircuit.circuit._limiting import _fetlim, _limvds, _pnjlim
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every limiter here is a numeric Newton aid; fix the toolkit once
    centrally rather than at each construction site."""
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


## ------------------------------------------------------------------------
## 1.  The laws themselves, against hand-worked values.

## (vnew, vold, vto, expected, why) -- `why` is the branch of DEVfetlim
## taken and the arithmetic that gives `expected`.
FETLIM_CASES = [
    ## vold < vto: OFF.  Turning on past `vto + 0.5` is clamped there.
    (5.0, 0.0, 1.0, 1.5,
     'off, delv=5>0, vnew=5 > vtemp=vto+0.5=1.5 -> vnew=vtemp'),
    ## Same branch, but landing under vtemp: the `delv > vtstlo` cap is
    ## the one that never fires (see test_fetlim_vtstlo_is_dead_code_in_devfetlim).
    (1.2, 0.0, 1.0, 1.2,
     'off, delv=1.2>0, vnew=1.2 <= vtemp=1.5, delv < vtstlo=|0-1|+1=2'),
    ## OFF, going further off: capped at vtsthi.
    (-10.0, 0.0, 1.0, -4.0,
     'off, delv=-10<=0, vtsthi=|2*(0-1)|+2=4, -delv=10>4 -> vold-vtsthi'),
    (-3.0, 0.0, 1.0, -3.0,
     'off, -delv=3 < vtsthi=4 -> untouched'),
    ## vto <= vold < vto+3.5: the MIDDLE region, a plain box.
    (-6.0, 2.0, 1.0, 0.5,
     'middle (1<=2<4.5), delv<=0 -> max(vnew, vto-0.5) = 0.5'),
    (1.0, 2.0, 1.0, 1.0,
     'middle, delv<=0, vnew=1.0 > vto-0.5=0.5 -> untouched'),
    (9.0, 2.0, 1.0, 5.0,
     'middle, delv>0 -> min(vnew, vto+4) = 5.0'),
    (4.0, 2.0, 1.0, 4.0,
     'middle, delv>0, vnew=4 < vto+4=5 -> untouched'),
    ## vold >= vto+3.5: strongly ON.
    (30.0, 6.0, 1.0, 18.0,
     'on (6>=4.5), delv=24>0, vtsthi=|2*(6-1)|+2=12, '
     'delv>=vtsthi -> vold+vtsthi=18'),
    (15.0, 6.0, 1.0, 15.0,
     'on, delv=9 < vtsthi=12 -> untouched'),
    (-100.0, 20.0, 1.0, 3.0,
     'on (20>=4.5), delv<0, vnew=-100 < vtox=4.5 -> max(vnew, vto+2)=3.0'),
    (5.0, 20.0, 1.0, 5.0,
     'on, delv<0, vnew=5 >= vtox=4.5, -delv=15 < vtstlo=|20-1|+1=20'),
    ## A negative vto (a depletion device / a p-channel folded to n) uses
    ## exactly the same law -- nothing here assumes vto > 0.
    (-6.0, -1.0, -2.0, -2.5,
     'middle (-2<=-1<1.5), delv<=0 -> max(vnew, vto-0.5) = -2.5'),
]

## (vnew, vold, expected, why) for DEVlimvds.
LIMVDS_CASES = [
    (100.0, 1.0, 4.0, 'vold<3.5, rising -> min(vnew, 4)'),
    (2.0, 1.0, 2.0, 'vold<3.5, rising, vnew=2 < 4 -> untouched'),
    (-50.0, 1.0, -0.5, 'vold<3.5, falling -> max(vnew, -0.5)'),
    (0.25, 1.0, 0.25, 'vold<3.5, falling, vnew > -0.5 -> untouched'),
    (100.0, 5.0, 17.0, 'vold>=3.5, rising -> min(vnew, 3*vold+2) = 17'),
    (1.0, 5.0, 2.0, 'vold>=3.5, falling, vnew=1 < 3.5 -> max(vnew, 2)'),
    (4.0, 5.0, 4.0, 'vold>=3.5, falling, vnew=4 not < 3.5 -> untouched'),
    (3.0, 3.5, 3.0, 'vold>=3.5 exactly, falling, vnew<3.5, max(3,2)=3'),
    ## The vold < 0 fold: mos1load.c's `-DEVlimvds(-vds, -vds_old)`.
    (-100.0, -1.0, -4.0, 'fold -> -min(100, 4) = -4'),
    (-100.0, -5.0, -17.0, 'fold -> -min(100, 3*5+2) = -17'),
    (50.0, -1.0, 0.5, 'fold -> -max(-50, -0.5) = 0.5'),
]


@pytest.mark.parametrize('vnew,vold,vto,expected,why', FETLIM_CASES)
def test_fetlim_matches_hand_worked_devfetlim(vnew, vold, vto, expected, why):
    got = _fetlim(vnew, vold, vto, numeric)
    assert got == expected, '%s: expected %r, got %r' % (why, expected, got)


@pytest.mark.parametrize('vnew,vold,expected,why', LIMVDS_CASES)
def test_limvds_matches_hand_worked_devlimvds(vnew, vold, expected, why):
    got = _limvds(vnew, vold, numeric)
    assert got == expected, '%s: expected %r, got %r' % (why, expected, got)


def test_fetlim_vtstlo_is_dead_code_in_devfetlim():
    """Both `vtstlo` caps in `DEVfetlim` are UNREACHABLE, so ngspice's
    documented "new definition for vtstlo" changes nothing.

    Proof, for the OFF branch: it needs `vold < vto`, `vnew <= vto + 0.5`
    and `delv > vtstlo`.  But `delv <= (vto - vold) + 0.5` while
    `vtstlo = (vto - vold) + 1`, so `delv < vtstlo` always.  For the ON
    branch: it needs `vold >= vto + 3.5` and `vnew >= vto + 3.5`, so
    `-delv <= (vold - vto) - 3.5`, while `vtstlo = (vold - vto) + 1`.

    That makes the SPICE3f5 value (`vtsthi/2 + 2 = |vold-vto| + 3`) and
    ngspice-47's (`|vold-vto| + 1`) indistinguishable.  The test asserts
    the OBSERVABLE consequence -- the two definitions agree everywhere --
    which is what a reader of the two sources actually needs to know.
    """
    def f5(vnew, vold, vto):
        """SPICE3f5's DEVfetlim: the pre-Gillespie vtstlo, nothing else."""
        vtsthi = abs(2 * (vold - vto)) + 2.0
        vtstlo = vtsthi / 2 + 2.0
        vtox = vto + 3.5
        delv = vnew - vold
        if vold >= vto:
            if vold >= vtox:
                if delv <= 0.0:
                    if vnew >= vtox:
                        if -delv > vtstlo:
                            vnew = vold - vtstlo
                    else:
                        vnew = max(vnew, vto + 2.0)
                elif delv >= vtsthi:
                    vnew = vold + vtsthi
            elif delv <= 0.0:
                vnew = max(vnew, vto - 0.5)
            else:
                vnew = min(vnew, vto + 4.0)
        else:
            if delv <= 0.0:
                if -delv > vtsthi:
                    vnew = vold - vtsthi
            else:
                vtemp = vto + 0.5
                if vnew <= vtemp:
                    if delv > vtstlo:
                        vnew = vold + vtstlo
                else:
                    vnew = vtemp
        return vnew

    n = 0
    for vto in np.linspace(-3.0, 3.0, 13):
        for vold in np.linspace(-40.0, 40.0, 161):
            for vnew in np.linspace(-200.0, 200.0, 201):
                assert _fetlim(vnew, vold, vto, numeric) == \
                    f5(vnew, vold, vto), (vnew, vold, vto)
                n += 1
    assert n > 400000


def test_fetlim_is_a_strict_no_op_near_the_solution():
    """A converged element must not be moved -- the property that lets
    "was limiting applied?" be a convergence signal.

    The threshold is not asymptotic: `_fetlim`'s tightest clamps sit
    0.5 V from `vto` and its step caps are never below 1 V, so **every**
    step of 0.45 V or less is passed through EXACTLY, whatever `vold` and
    `vto` are.  That band is asserted rather than "small", and it is
    where the band actually is: at 0.6 V the middle region's `vto - 0.5`
    clamp starts to bite, which the second half of this test pins.
    """
    moved = 0
    for vto in np.linspace(-3.0, 3.0, 25):
        for vold in np.linspace(-20.0, 20.0, 161):
            for d in (-0.45, -0.2, -1e-9, 0.0, 1e-9, 0.2, 0.45):
                vnew = vold + d
                assert _fetlim(vnew, vold, vto, numeric) == vnew, \
                    (vold, vto, d)
            ## and the band is not vacuous: a step of 0.6 IS limited
            ## somewhere in this sweep.
            for d in (-0.6, 0.6):
                if _fetlim(vold + d, vold, vto, numeric) != vold + d:
                    moved += 1
    assert moved > 0, 'the 0.45 band would be untested if nothing bigger moved'


def test_limvds_is_a_strict_no_op_near_the_solution():
    """Same property for `limvds`, including at NEGATIVE vds -- which is
    what the `vold < 0` fold is for.  Without it every reverse-mode
    solution below -0.5 V would be dragged to -0.5 on every iteration.
    """
    for vold in np.linspace(-20.0, 20.0, 401):
        for d in (-0.4, -0.1, 0.0, 0.1, 0.4):
            vnew = vold + d
            assert _limvds(vnew, vold, numeric) == vnew, (vold, d)
    ## The unfolded routine does NOT have the property: this is the
    ## measured size of what the fold buys.
    def unfolded(vnew, vold):
        if vold >= 3.5:
            if vnew > vold:
                return min(vnew, 3.0 * vold + 2.0)
            return max(vnew, 2.0) if vnew < 3.5 else vnew
        return min(vnew, 4.0) if vnew > vold else max(vnew, -0.5)
    assert unfolded(-3.1, -3.0) == -0.5      # dragged 2.6 V, every iteration
    assert _limvds(-3.1, -3.0, numeric) == -3.1


def test_limvds_compresses_a_wild_step():
    """The other half: a step that IS wild is bounded, and bounded to
    SPICE's own numbers rather than merely to something finite."""
    ## From a 1 V bias, one iteration can reach 4 V and no further.
    assert _limvds(1e6, 1.0, numeric) == 4.0
    ## and from there the ceiling triples each iteration: 4 -> 14 -> 44.
    assert _limvds(1e6, 4.0, numeric) == 14.0
    assert _limvds(1e6, 14.0, numeric) == 44.0


def test_fetlim_compresses_a_wild_step():
    ## An off device asked to jump to 100 V lands 0.5 V above threshold.
    assert _fetlim(100.0, 0.0, 1.0, numeric) == 1.5
    ## An on device asked to collapse to -100 V lands 2 V above threshold.
    assert _fetlim(-100.0, 20.0, 1.0, numeric) == 3.0


## ------------------------------------------------------------------------
## 2.  The declarations, on real elements.

def _fet(which='both', chained=False):
    """A saturating exponential FET: subthreshold in `vgs`, `Vdsat`-like
    saturation in `vds`.  Both are what a limiter exists for -- the first
    overflows, the second collapses `dIds/dVds` to nothing.

    `which` selects the declaration under test: `'none'`, `'fet'`, `'vds'`,
    `'both'` (two independent per-probe limiters), `'group'`
    (`limit_together`, the device-level form of roadmap 10.3(b)) or
    `'seq'` (the same group with SPICE's `mos1load.c` sequencing).  The
    last two are used by `test_device_limiter.py`, which measures them
    against `'both'` on this same circuit -- the comparison only means
    anything if the MODEL is identical, which is why they live here beside
    it rather than in a second copy.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, limexp, limit_fet, limit_vds,
                                       limit_together, vt)

    class _M(Behavioural):
        instparams = [Parameter(name='VTO', desc='threshold', unit='V',
                                default=0.5),
                      Parameter(name='IS0', desc='current scale', unit='A',
                                default=1e-6),
                      Parameter(name='N', desc='subthreshold slope', unit='',
                                default=1.3),
                      Parameter(name='VE', desc='saturation voltage',
                                unit='V', default=0.5)]

        @staticmethod
        def analog(d, g, s):
            bgs, bds = Branch(g, s), Branch(d, s)
            if which in ('group', 'seq'):
                vgs, vds = limit_together(               # noqa
                    limit_fet(bgs.V, VTO), limit_vds(bds.V),
                    sequential=(which == 'seq'))
            else:
                vgs = limit_fet(bgs.V, VTO) if which in ('both', 'fet') \
                    else bgs.V                                      # noqa
                vds = limit_vds(bds.V) if which in ('both', 'vds') \
                    else bds.V                                      # noqa
            ids = IS0 * limexp((vgs - VTO) / (N * vt())) \
                * (1.0 - limexp(-vds / VE))                         # noqa
            if chained:
                ids = var(ids, 'ids')
            return Contribution(bds.I, ids)
    return _M


def test_limit_fet_is_generated_only_when_asked():
    from pycircuit.circuit.hdl import explain
    assert hasattr(_fet('both'), 'limit')
    assert hasattr(_fet('fet'), 'limit')
    assert hasattr(_fet('vds'), 'limit')
    assert not hasattr(_fet('none'), 'limit')
    ## `explain()` names the KIND, not just a count -- with three kinds a
    ## bare "2 $limit probes" no longer says what the model asked for.
    assert '2 $limit probes (fetlim on (g,s), limvds on (d,s))' \
        in explain(_fet('both'))
    assert '1 $limit probe (limvds on (d,s))' in explain(_fet('vds'))


def test_limit_fet_and_limit_vds_reject_a_non_branch_probe():
    from pycircuit.circuit.hdl import limit_fet, limit_vds, Node, Quantity
    v = Quantity('V', Node('a'))
    with pytest.raises(ValueError, match='branch potential'):
        limit_fet(v, 0.5)
    with pytest.raises(ValueError, match='branch potential'):
        limit_vds(v)
    with pytest.raises(ValueError, match='branch potential'):
        limit_fet(1.0, 0.5)


def test_limit_fet_compiles_through_both_code_generators():
    """Eager (flattened) and chained (`var()`) are two separate code
    generators; a feature that works in one is not thereby in the other.

    Before this change the chained path silently had NO limiter at all:
    the metaclass returned early for `pure_spec is None`, jumping over
    the `$limit` install, and `var(limit_pnj(...))` additionally carried
    its marker into `_chain_compile` and died there.  Both are fixed, so
    the two generators must now agree numerically -- limiting is
    arithmetic on the solution vector and has no business differing.
    """
    eager = _fet('both', chained=False)('d', 'g', 's')
    chained = _fet('both', chained=True)('d', 'g', 's')
    assert chained._hdl_info['chained'] and not eager._hdl_info['chained']
    for el in (eager, chained):
        el.update_iparv()
    rng = np.random.default_rng(20260825)
    for _ in range(300):
        x = rng.uniform(-30.0, 30.0, 3)
        x0 = rng.uniform(-30.0, 30.0, 3)
        a = eager.limit(x, x0)
        b = chained.limit(x, x0)
        assert (a == b).all(), (x, x0, a, b)


def test_limit_fet_on_an_element_is_the_bare_law_on_each_branch():
    """`limit()` must carry each probe's own limited value exactly, and
    those values must be `_fetlim`/`_limvds` of that probe -- nothing in
    the element wrapper is allowed to add its own arithmetic."""
    el = _fet('both')('d', 'g', 's', VTO=0.8)
    el.update_iparv()
    ## terminal order is (d, g, s) -- the analog() argument names -- and
    ## the two probes move `g` and `d`, leaving the source alone.
    assert el.terminals == ['d', 'g', 's']
    assert [(s[0], s[1], s[2]) for s in el._hdl_info['limit_spec']] == \
        [((1, 2), 'fet', 1), ((0, 2), 'vds', 0)]
    ## THE EXACT LAW is asserted on SINGLE-PROBE elements.  With two
    ## probes sharing the source, the loop is sequential -- the second
    ## probe reads the vector the first already wrote -- so a two-probe
    ## element cannot be checked against a law computed from the input
    ## alone without re-implementing the algorithm in the test.  The
    ## shared case is covered by
    ## `test_limit_probes_sharing_a_plus_terminal_do_not_undo_each_other`,
    ## which asserts the invariant that survives any order.
    rng = np.random.default_rng(7)
    for which, law, probe in (('fet', lambda a, b: _fetlim(a, b, 0.8, numeric),
                               (1, 2)),
                              ('vds', lambda a, b: _limvds(a, b, numeric),
                               (0, 2))):
        one = _fet(which)('d', 'g', 's', VTO=0.8)
        one.update_iparv()
        ra, rb = probe
        moved_minus = 0
        for _ in range(500):
            x = rng.uniform(-40.0, 40.0, 3)
            x0 = rng.uniform(-40.0, 40.0, 3)
            out = one.limit(x, x0)
            v = law(x[ra] - x[rb], x0[ra] - x0[rb])
            ## Asserted on the WRITE, not on the recomputed difference:
            ## `(s + v) - s` is not `v` in floating point, and a rounding
            ## error is not what this test is trying to catch.  WHICH end
            ## moves is a runtime choice now, so both spellings are legal
            ## -- but exactly one of them must hold, and the other
            ## terminal must be untouched.
            if out[ra] == x[ra] and out[rb] == x[rb]:
                ## Nothing moved, which is legal ONLY when the limiter
                ## did not bite.  `(a - b)` recovered from a write is not
                ## bit-equal to the original difference, so the no-op
                ## case has to be recognised rather than folded into one
                ## of the two write spellings below.
                assert v == x[ra] - x[rb], (x, x0, out, v)
            elif out[ra] != x[ra]:
                assert out[rb] == x[rb] and out[ra] == x[rb] + v, (x, x0, out)
            else:
                moved_minus += 1
                assert out[rb] == x[ra] - v, (x, x0, out)
            ## the third node is not this probe's business at all
            third = ({0, 1, 2} - {ra, rb}).pop()
            assert out[third] == x[third]
            ## State-free: the input vector is not mutated.
            assert not np.shares_memory(out, x)
        ## And the minus terminal is no longer sacrosanct: under the old
        ## fixed rule it never moved, which is exactly what let a wild
        ## source drag a sane drain out with it.
        assert moved_minus > 20, (which, moved_minus)


def test_limit_probes_sharing_a_plus_terminal_do_not_undo_each_other():
    """The ordering rule, and the case it exists for.

    Two probes sharing a terminal would have the second write undo the
    first if both moved the same node (the hazard `limit_junctions`
    handles with its `move` field for a BJT's two junctions on one
    base).  Which node each probe moves is decided at RUN TIME -- the
    terminal that drifted further from the last accepted point, the
    shared one going to the larger correction (roadmap 12.1) -- and the
    invariant is that the two corrections land on two DIFFERENT rows.
    (This docstring used to describe the compile-time rule, "the second
    probe takes its MINUS terminal"; that rule was replaced on
    2026-08-25 and the assertions below were already about the
    invariant rather than the node identities.)
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, limit_fet, vt)

    class _Pair(Behavioural):
        instparams = [Parameter(name='VTO', desc='vt', unit='V',
                                default=0.5)]

        @staticmethod
        def analog(d, g, s):
            bgs, bgd = Branch(g, s), Branch(g, d)
            vgs = limit_fet(bgs.V, VTO)                             # noqa
            vgd = limit_fet(bgd.V, VTO)                             # noqa
            return (Contribution(bgs.I, 1e-6 * limexp(vgs / vt())),
                    Contribution(bgd.I, 1e-6 * limexp(vgd / vt())))

    ## `move` for the two probes: the plus (g) for the first, the minus
    ## (d) for the second, because g is already claimed.
    spec = _Pair._hdl_info['limit_spec']
    assert [s[1] for s in spec] == ['fet', 'fet']
    assert spec[0][2] == spec[0][0][0]          # first moves its plus, g
    assert spec[1][2] == spec[1][0][1]          # second moves its minus, d

    el = _Pair('d', 'g', 's', VTO=0.5)
    el.update_iparv()
    rng = np.random.default_rng(11)
    distinct = 0
    for _ in range(400):
        x = rng.uniform(-30.0, 30.0, 3)
        x0 = rng.uniform(-30.0, 30.0, 3)
        out = el.limit(x, x0)
        ## THE PROPERTY THAT MATTERS, and the only one that survives the
        ## write-back becoming a runtime choice: no probe is undone.
        ## Each of the two branches must end at a value its own limiter
        ## produced -- not at whatever the other probe's write left
        ## behind.  Asserted structurally, because WHICH node each probe
        ## moved is no longer fixed at compile time.
        ##
        ## The loop is still SEQUENTIAL, as `limit_junctions` is, so the
        ## second probe legitimately sees the first probe's write; what
        ## it may not do is have its own correction overwritten.
        movers = [i for i in range(3) if out[i] != x[i]]
        ## Two probes, two corrections, and they may not land on the
        ## same node -- that IS the no-undo property.
        assert len(set(movers)) == len(movers)
        assert len(movers) <= 2, (movers, x, x0)
        if len(movers) == 2:
            distinct += 1
        ## Every node that did not move is untouched, exactly.
        for i in range(3):
            if i not in movers:
                assert out[i] == x[i]
        ## And the vector is finite: the failure this whole change
        ## exists to remove was a write-back copying an unbounded
        ## absolute value across a branch.
        assert np.all(np.isfinite(out))
        assert np.max(np.abs(out)) <= np.max(np.abs(x)) + 60.0, (x, out)
    ## Both probes really do bite on a decent fraction of the draws --
    ## otherwise the no-undo property above would be vacuous.
    assert distinct > 50, distinct


def test_the_vgs_vds_star_is_order_independent_and_is_not_spice_s_order():
    """The ordering question, recorded as a measurement rather than prose.

    `limit_fet(V(g,s))` and `limit_vds(V(d,s))` write to different
    terminals (`g` and `d`), so each probe ends up with exactly its own
    limited value and DECLARATION ORDER DOES NOT MATTER.  SPICE's order
    does matter: `mos1load.c` limits `vgs`, then recomputes
    `vds = vgs - vgd` from the UNLIMITED `vgd` before calling `limvds`,
    so its `limvds` argument is shifted by `delta = vgs_lim - vgs`.

    This test pins the size of that difference on a concrete step rather
    than describing it.  The device-level limiter (roadmap 10.3(b)) can
    now produce SPICE's number --
    `limit_together(..., sequential=True)`, measured against the same
    step in `test_device_limiter.py` -- but the PER-PROBE form measured
    here is unchanged and is meant to stay so: this is still the
    difference between the two, not a gap waiting to be closed.
    """
    el = _fet('both')('d', 'g', 's', VTO=0.5)
    el.update_iparv()

    ## Order independence, asserted rather than argued.  Source order is
    ## not the knob -- `generate_code` sorts the markers into sympy's
    ## canonical order, so writing the two calls the other way round
    ## produces the same spec.  The spec list itself is the knob: it is
    ## what the generated `limit()` closes over.  (Each `_fet(...)` call
    ## builds a fresh class, so this touches nothing else, and it is put
    ## back either way.)
    spec = el._hdl_info['limit_spec']
    rng = np.random.default_rng(3)
    try:
        for _ in range(200):
            xa = rng.uniform(-30.0, 30.0, 3)
            x0a = rng.uniform(-30.0, 30.0, 3)
            spec.reverse()
            flipped = el.limit(xa, x0a)
            spec.reverse()
            assert (el.limit(xa, x0a) == flipped).all(), (xa, x0a)
    finally:
        if spec is not el._hdl_info['limit_spec']:      # pragma: no cover
            raise AssertionError('spec object was replaced, not reordered')

    ## x = (d, g, s).  A wild gate step from an off device.
    x0 = np.array([0.2, 0.0, 0.0])
    x = np.array([9.0, 40.0, 0.0])
    out = el.limit(x, x0)

    vgs, vgs0 = x[1] - x[2], x0[1] - x0[2]
    vds, vds0 = x[0] - x[2], x0[0] - x0[2]
    vgs_lim = _fetlim(vgs, vgs0, 0.5, numeric)
    assert vgs_lim == 1.0                       # off -> clamped to vto+0.5
    delta = vgs_lim - vgs                       # -39.0

    ours = _limvds(vds, vds0, numeric)
    ## mos1load.c: `vds = vgs - vgd` with the UNLIMITED vgd, which is the
    ## same number as `vds + delta`.
    vgd = vgs - vds
    assert vgs_lim - vgd == vds + delta == -30.0
    spice = _limvds(vgs_lim - vgd, vds0, numeric)

    assert out[1] == x[2] + vgs_lim             # both agree on vgs
    assert out[0] == x[2] + ours
    ## The measured gap, on this step: 4.0 V here against SPICE's -0.5 V,
    ## because SPICE's limvds is handed -30 and clamps it at -0.5.  A
    ## 4.5 V difference in where the next Jacobian is taken; the solution
    ## the iteration is converging to is untouched by either.
    assert ours == 4.0 and spice == -0.5


## ------------------------------------------------------------------------
## 3.  The `limit_pnj` regression -- the restructure had to be invisible.

def _limited_diode(chained=False):
    """`chained` is three-valued, because the chained path had TWO
    separate pre-existing failures and they are not the same bug:

    * ``'wrapped'`` -- `var(limit_pnj(...))`.  The `$limit` marker was
      stripped from contribution statements only, never from `var()`
      definitions, so it reached `_chain_compile` and raised
      "Unsupported by ... _Limit".  A crash, at class creation.
    * ``'elsewhere'`` -- `limit_pnj` in the contribution, `var()` used
      for something else.  That compiled cleanly and produced a class
      with NO `limit()` method, because the metaclass returned early for
      `pure_spec is None`.  Silent.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, limexp, limit_pnj, vt)

    class _D(Behavioural):
        instparams = [Parameter(name='IS', desc='Is', unit='A',
                                default=1e-13)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            v = limit_pnj(b.V, IS, vt())                            # noqa
            if chained == 'wrapped':
                v = var(v, 'vl')
            scale = var(IS, 'scale') if chained == 'elsewhere' else IS  # noqa
            return Contribution(b.I, scale * (limexp(v / vt()) - 1))
    return _D


@pytest.mark.parametrize('ISv', [1e-18, 1e-13, 1e-6, 1.0, 0.0])
def test_limit_pnj_is_unchanged_by_the_limit_spec_restructure(ISv):
    """`limit_spec` now carries a limiter KIND and its own parameters,
    and `limit()` dispatches on it.  For `pnjlim` the result must be
    BIT-IDENTICAL to what the single hardcoded `(IS, VT)` pair produced:
    move the plus terminal, carry `_pnjlim` exactly.

    `_pnjlim` itself is untouched by this change, so asserting against it
    is asserting against the old code.  (Checked once more directly: the
    generated `limit()` was run against the one compiled from `hdl.py` at
    8e7e538 over 525310 (vold, vnew, IS, offset) points on both code
    paths, with zero mismatches.  That comparison cannot be committed --
    it needs the old file -- so this is its standing form.)
    """
    from pycircuit.circuit.constants import kboltzmann, qelectron
    el = _limited_diode()('p', 'n', IS=ISv)
    el.update_iparv()
    ## VT is read back from the spec's OWN compiled `vt()` rather than
    ## recomputed here.  `kboltzmann * T / qelectron` written out in
    ## Python differs from what `lambdify` constant-folds by one ulp, and
    ## `_pnjlim` puts VT inside a log -- so a hand-written VT turns a
    ## bit-identity test into a 1e-15 one and stops being the claim.
    ## What the restructure touched is the DISPATCH and the write-back,
    ## and those are what is asserted exactly.
    (ra, rb), kind, move, (isf, vtf) = el._hdl_info['limit_spec'][0]
    assert (ra, rb, kind, move) == (0, 1, 'pnj', 0)
    VT = float(vtf(ISv, defaultepar.T))
    assert float(isf(ISv, defaultepar.T)) == ISv
    assert VT == pytest.approx(float(kboltzmann * defaultepar.T / qelectron),
                               rel=1e-15)
    limited = 0
    for vold in np.linspace(-5.0, 1.5, 66):
        for vnew in np.linspace(-300.0, 300.0, 121):
            for shift in (0.0, 2.7):
                out = el.limit(np.array([vnew + shift, shift]),
                               np.array([vold + shift, shift]))
                want = _pnjlim(vnew, vold, VT, ISv, numeric)
                assert out[1] == shift          # the plus terminal moves
                assert out[0] == shift + want, (vold, vnew, ISv)
                limited += (want != vnew)
    ## Not vacuous: `IS = 0` is the "no junction" pass-through and is
    ## expected to move nothing; every other card must actually limit.
    assert (limited == 0) == (ISv == 0.0), limited


@pytest.mark.parametrize('chained', ['wrapped', 'elsewhere'])
def test_limit_pnj_now_reaches_the_chained_code_generator(chained):
    """Regression on the two chained-path failures `_limited_diode`
    documents -- the crash and the silent one.

    Both were verified against `hdl.py` at 8e7e538: ``'wrapped'`` raised
    `PrintMethodNotImplementedError` at class creation, and
    ``'elsewhere'`` produced a class with no `limit` in its `__dict__`
    while `explain()` cheerfully reported "1 $limit probe" -- because
    `explain` reads the compiled info, not the class.
    """
    from pycircuit.circuit.hdl import explain
    cls = _limited_diode(chained=chained)
    assert cls._hdl_info['chained']
    assert 'limit' in cls.__dict__, 'chained model lost its $limit'
    assert '1 $limit probe (pnjlim on (plus,minus))' in explain(cls)

    a = _limited_diode(False)('p', 'n', IS=1e-13)
    b = cls('p', 'n', IS=1e-13)
    a.update_iparv(); b.update_iparv()
    moved = 0
    for vnew in np.linspace(-50.0, 50.0, 201):
        x = np.array([vnew, 0.0])
        x0 = np.array([0.1, 0.0])
        out = a.limit(x, x0)
        assert (out == b.limit(x, x0)).all()
        moved += (out[0] != vnew)
    assert moved > 0


def test_cross_now_reaches_the_chained_code_generator_too():
    """The same early `return` also swallowed `@cross`, and this pins the
    other half of that fix.

    Verified against `hdl.py` at 8e7e538: a `var()`-using model with a
    `Cross` statement had `cross_spec` compiled and `explain()` reporting
    "1 @cross event", but no `accept_step`/`next_event` of its own -- it
    silently inherited `Circuit`'s.  The prediction is not new code; it
    was built and then never installed.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Cross, var, explain)

    def comparator(chained):
        class _C(Behavioural):
            instparams = [Parameter(name='G', desc='g', unit='S',
                                    default=1e-3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                v = var(b.V, 'vv') if chained else b.V
                return (Contribution(b.I, G * v),                   # noqa
                        Cross(v - 0.5, +1))
        return _C

    for chained in (False, True):
        cls = comparator(chained)
        assert cls._hdl_info['chained'] is chained
        assert '1 @cross event' in explain(cls)
        assert 'accept_step' in cls.__dict__, chained
        assert 'next_event' in cls.__dict__, chained
        el = cls('a', 'b')
        el.update_iparv()
        el.accept_step(0.0, np.array([0.0, 0.0]))
        el.accept_step(1.0, np.array([0.25, 0.0]))
        ## rising at 0.25 V/s from 0.25 V at t=1 reaches 0.5 V at t=2.
        assert el.next_event(1.0) == 2.0


## ------------------------------------------------------------------------
## 4.  The point of the whole thing: a DC solve that needs it.

def _cascode(cls, vdd, vg2, vg1):
    """Two saturating FETs in series, the middle node held by nothing but
    their two drain conductances.  This is the case `compact.py`'s
    `PspMosLongChannel.limit` was written for: once both devices
    saturate, `dIds/dVds` collapses and the middle row goes numerically
    empty, so a Newton step that walks out there stays."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vg2'] = VS('g2', gnd, v=vg2)
    c['vg1'] = VS('g1', gnd, v=vg1)
    c['M2'] = cls('vdd', 'g2', 'mid')
    c['M1'] = cls('mid', 'g1', gnd)
    return c


def _plain_newton(c, maxiter=100, reltol=1e-4, vabstol=1e-6, iabstol=1e-12):
    """`StandardNewton` on the circuit, WITHOUT `DCAnalysis`'s rescue
    chain (junction-gmin, gshunt, source stepping, pseudo-transient).
    That chain would rescue this circuit either way, which is exactly why
    it is not used here: it would hide what the limiter is worth."""
    from pycircuit.circuit.nrsolver import StandardNewton
    iref = c.get_node_index(gnd)
    epar = defaultepar

    def rfunc(xr):
        x = np.insert(xr, iref, 0.0)
        F = c.i(x, epar) + c.u(0, analysis='dc', epar=epar)
        J = c.G(x, epar)
        return (np.delete(F, iref),
                np.delete(np.delete(J, iref, 0), iref, 1))

    def limiter(xr, x0r):
        x = c.limit(np.insert(xr, iref, 0.0), np.insert(x0r, iref, 0.0),
                    epar)
        return np.delete(x, iref)

    nn, nb = len(c.nodes), len(c.branches)
    abstol = np.delete(np.concatenate((iabstol * np.ones(nn),
                                       vabstol * np.ones(nb))), iref)
    xtol = np.delete(np.concatenate((vabstol * np.ones(nn),
                                     iabstol * np.ones(nb))), iref)
    x, its = StandardNewton().solve_system(
        np.zeros(c.n - 1), rfunc, numeric, reltol, abstol, xtol, maxiter,
        limiter=limiter)
    return np.insert(x, iref, 0.0), its


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_dc_solve_needs_the_limiter_and_gets_the_same_answer():
    """The only test here that shows the feature doing its job.

    Same circuit, same solver, same starting point; the ONLY difference
    is whether the element declares `limit_fet`/`limit_vds`.  Without
    them plain Newton does not converge; with them it does, and it lands
    on the operating point the full `DC` rescue chain finds for the
    unlimited element -- a limiter moves the PATH, never the answer.
    """
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix

    cond = dict(vdd=20.0, vg2=2.0, vg1=0.8)

    ## Unlimited: plain Newton runs out of iterations.
    with pytest.raises((NoConvergenceError, SingularMatrix)):
        _plain_newton(_cascode(_fet('none'), **cond))

    ## Limited: it converges, and well inside the iteration budget.
    c = _cascode(_fet('both'), **cond)
    x, its = _plain_newton(c)
    mid = x[c.get_node_index('mid')]
    assert its <= 40, its

    ## The answer is a real cascode bias, not a rail.
    assert 0.5 < mid < cond['vdd'] - 0.5

    ## and it is the SAME point the rescue chain reaches with no limiter
    ## at all.  `vabstol` is 1e-6 and the node sits at ~1.2 V, so 1e-6
    ## relative is the solver's own floor, not a slack tolerance.
    ref = DC(_cascode(_fet('none'), **cond), toolkit=numeric).solve()
    assert_allclose(mid, float(ref.v('mid')), rtol=1e-6)


@pytest.mark.filterwarnings('ignore:overflow encountered')
@pytest.mark.parametrize('chained', [False, True])
def test_dc_solve_attributes_the_rescue_to_each_limiter(chained):
    """Which limiter does the work -- measured, not assumed.

    At (vdd, vg2, vg1) = (20, 2, 0.8), on this machine:

        none   NoConvergenceError
        fet    12 iterations
        vds    57
        both   13

    **`both` was recorded as 30 here and in roadmap 12.1, and it never
    was 30 after the write-back fix** (corrected 2026-08-25, by re-running
    it at the commit that introduced the table).  30 is the number from
    the row above it -- the BEFORE table -- copied down without being
    re-measured, and it made a one-iteration difference read as a 2.5x
    penalty.  The conclusion it was used for ("adding the second probe
    costs iterations") is true at this point and by exactly one; across
    the 48-point grid of `test_device_limiter.py` it is false, and `both`
    is the cheapest of the four.  A ninth instance in this campaign of a
    claim right in its conclusion and wrong in its stated reason -- here
    the reason was a number nobody re-ran.

    **The attribution changed completely when the write-back stopped
    being fixed at compile time**, and the old numbers are worth keeping
    because of what they said. They were:

        none   NoConvergenceError
        fet    NoConvergenceError    <- and the docstring explained why
        vds    78
        both   30

    The explanation given was that SPICE's `fetlim` bounds a step to
    about a volt, and a volt is forty thermal voltages to a subthreshold
    exponential -- plausible, self-consistent, and NOT the reason.
    `limit_fet` was failing because its write-back always moved the
    GATE, so a wild source dragged the gate out with it; once the
    correction goes to whichever terminal actually drifted, `fetlim`
    alone is the BEST of the four.

    That `both` (13) is worse than `fet` alone (12) is real, and it is
    one iteration: two probes may not move the same terminal, so the
    second is pushed onto a node it would not have chosen.  The
    device-level limiter (roadmap 10.3(b), `limit_together`) removes that
    constraint and measures 13 as well -- see
    `test_device_limiter.py::test_the_grouped_form_costs_about_what_the_
    per_probe_one_costs`.  What it buys is not this number.

    Run on BOTH code generators, since the limiter must not depend on
    which one compiled the element.
    """
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    cond = dict(vdd=20.0, vg2=2.0, vg1=0.8)

    ## Unlimited still fails -- otherwise none of the rest is evidence.
    with pytest.raises((NoConvergenceError, SingularMatrix)):
        _plain_newton(_cascode(_fet('none', chained), **cond))

    its = {}
    for which in ('fet', 'vds', 'both'):
        c = _cascode(_fet(which, chained), **cond)
        x, its[which] = _plain_newton(c)
        ## Every variant that converges must converge to the SAME point:
        ## a limiter moves the path, never the solution.
        assert_allclose(x[c.get_node_index('mid')], 1.2031737, rtol=1e-6)

    ## `fetlim` alone is now the cheapest rescue, by a wide margin.
    assert its['fet'] < its['vds'] / 3.0, its
    ## ... and adding the second probe COSTS iterations, which is the
    ## shared-terminal constraint showing up as a number.
    assert its['both'] > its['fet'], its


## ---------------------------------------------------------------------
## Limiter parameters that read the solution (2026-08-25)
## ---------------------------------------------------------------------

def _von_fet(chained, grouped=False):
    """A FET whose `limit_fet` threshold moves with the bulk bias --
    `von = VTO + K * V(b,s)` -- so the parameter reads the solution."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       var, limexp, limit_fet, limit_vds,
                                       limit_together, vt)

    class _V(Behavioural):
        instparams = [Parameter(name='VTO', desc='threshold', unit='V',
                                default=0.5),
                      Parameter(name='K', desc='body factor', unit='',
                                default=0.5)]

        @staticmethod
        def analog(d, g, s, b):
            bgs, bds, bbs = Branch(g, s), Branch(d, s), Branch(b, s)
            von = VTO + K * bbs.V                                  # noqa
            if chained:
                von = var(von, 'von')
            if grouped:
                vgs, vds = limit_together(limit_fet(bgs.V, von),   # noqa
                                          limit_vds(bds.V))
            else:
                vgs, vds = limit_fet(bgs.V, von), limit_vds(bds.V)
            ids = 1e-6 * limexp((vgs - von) / (1.3 * vt())) \
                * (1.0 - limexp(-vds / 0.5))
            return Contribution(bds.I, var(ids, 'ids') if chained else ids)
    return _V


@pytest.mark.parametrize('chained', [False, True])
@pytest.mark.parametrize('grouped', [False, True])
def test_a_solution_dependent_parameter_is_taken_at_x0_not_x(chained,
                                                             grouped):
    """The discriminating check for roadmap 12.6(c).

    `fetlim` from an OFF device compresses a wild gate step to
    `von + 0.5`.  Put the bulk at one bias in the last accepted iterate
    and at another in the proposed one: which `von` the clamp used is
    then visible in the landing point.  SPICE's semantics -- and the
    only well-defined choice -- is the last accepted iterate.
    """
    el = _von_fet(chained, grouped)('d', 'g', 's', 'b')
    el.update_iparv()
    fet = [sp for sp in el._hdl_info['limit_spec'] if sp[1] == 'fet'][0]
    assert fet[3][0]._wants_x is True
    ## x = (d, g, s, b).  vbs = 0.6 at x0, 2.4 at x.
    x0 = np.array([1.0, 0.0, 0.0, 0.6])
    x = np.array([1.0, 100.0, 0.0, 2.4])
    out = el.limit(x, x0)
    von_x0 = 0.5 + 0.5 * 0.6
    von_x = 0.5 + 0.5 * 2.4
    assert_allclose(out[1] - out[2], von_x0 + 0.5, rtol=1e-12)
    assert abs((out[1] - out[2]) - (von_x + 0.5)) > 0.5    # and NOT x
    ## the bulk itself is not a probe terminal and is never written
    assert out[3] == x[3]


def test_a_solution_dependent_parameter_on_both_generators_agrees():
    """Flat and chained spellings of the same `von` must give the same
    limiter -- the parameter goes through `_chain_compile` on one path
    and through an unpack on the other."""
    a = _von_fet(False)('d', 'g', 's', 'b'); a.update_iparv()
    b = _von_fet(True)('d', 'g', 's', 'b'); b.update_iparv()
    rng = np.random.default_rng(5)
    for _ in range(300):
        x0 = rng.uniform(-3.0, 3.0, 4)
        x = rng.uniform(-40.0, 40.0, 4)
        assert (a.limit(x, x0) == b.limit(x, x0)).all(), (x, x0)


def test_a_time_dependent_limiter_parameter_is_still_refused():
    """`x0` is well-defined for a limiter; `t` is not offered to it."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, limit_pnj, vt, TIME)
    with pytest.raises(ValueError, match='cannot depend on time'):
        class _T(Behavioural):
            instparams = [Parameter(name='IS', desc='IS', unit='A',
                                    default=1e-14)]

            @staticmethod
            def analog(plus, minus):
                bb = Branch(plus, minus)
                v = limit_pnj(bb.V, IS * (1.0 + TIME), vt())     # noqa
                return Contribution(bb.I, IS * (limexp(v / vt()) - 1))


def test_explain_names_a_solution_dependent_limiter_parameter():
    from pycircuit.circuit.hdl import explain
    el = _von_fet(True)('d', 'g', 's', 'b'); el.update_iparv()
    text = explain(el)
    assert 'fetlim on (g,s) [params at last iterate]' in text, text
    assert 'limvds on (d,s)' in text and \
        'limvds on (d,s) [params' not in text, text
