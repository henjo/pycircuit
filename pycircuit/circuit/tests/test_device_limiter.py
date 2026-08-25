"""Device-level Newton limiting: `limit_together`, roadmap 10.3(b).

A per-probe limiter writes each probe back by moving ONE endpoint, so two
probes that share a terminal have to negotiate for it and the loser is
pushed onto a node it would not have chosen -- and four probes on a
four-terminal MOSFET cannot be written back at all, which `generate_code`
refuses at compile time.  `limit_together` hands the whole device to one
write-back: the probes are a graph over the terminals, the least-drifted
node is held, and every other node is derived from it.

**What it is worth, measured, because the roadmap's reason for wanting it
was a stale number.**  Roadmap 12.1 recorded `both` = 30 iterations
against `fet` = 12 and called the gap the thing a device-level limiter
would close.  Re-measured here on the same circuit at the same operating
point, `both` is **13**, not 30 -- the 30 is the number from BEFORE the
write-back fix, copied into the after-table without being re-run.  Over
48 operating points the two-probe `both` is in fact the best per-probe
variant there is, and the grouped form does not beat it: 927 Jacobian
evaluations against 909, never better at a point, worse by exactly one
iteration at 14 of the 48.

The justification that survives measurement is a CAPABILITY, not a count.
SPICE's own MOSFET declaration -- `fetlim` on `vgs`, `limvds` on `vds`,
`pnjlim` on each bulk junction -- is four probes on four terminals with a
cycle in it, and the per-probe form refuses to compile it.  Grouped, it
compiles, and on a body-biased cascode it takes 896 Jacobian evaluations
against 3478 unlimited, with six operating points that do not converge at
all without it.

Sections: 1 the declaration, 2 the write-back's properties, 3 SPICE's
ordering as a measured difference, 4 the DC solves.
"""
import numpy as np
from numpy.testing import assert_allclose
import pytest

import pycircuit.circuit.circuit
from pycircuit.circuit import numeric
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.elements import VS
from pycircuit.circuit._limiting import _fetlim, _limvds, _pnjlim
from pycircuit.circuit.tests.test_limit_fet import (_fet, _cascode,
                                                    _plain_newton)
from pycircuit.utilities.param import Parameter


@pytest.fixture(autouse=True)
def _numeric_toolkit():
    """Every limiter here is a numeric Newton aid; fix the toolkit once
    centrally rather than at each construction site."""
    old = pycircuit.circuit.circuit.default_toolkit
    pycircuit.circuit.circuit.default_toolkit = numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = old


## The reference operating point, shared by everything below.
COND = dict(vdd=20.0, vg2=2.0, vg1=0.8)
## and the grid the aggregate measurements run on.  Four supplies, three
## cascode gate biases, four input biases -- chosen to span cutoff,
## subthreshold and full saturation, not chosen to flatter a variant.
GRID = [(vdd, vg2, vg1)
        for vdd in (2.0, 5.0, 20.0, 40.0)
        for vg2 in (1.0, 2.0, 3.0)
        for vg1 in (0.6, 0.8, 1.2, 2.0)]


def _mos4(mode='group'):
    """SPICE's OWN MOSFET limiting declaration: `fetlim` on `vgs`,
    `limvds` on `vds`, and `pnjlim` on each of the two bulk junctions.

    Four probes over four terminals, and `(b,s)`, `(b,d)` and `(d,s)`
    close a triangle -- so the constraints need not be consistent and one
    of them may have to go.  `mode` is `'none'`, `'probe'` (per-probe,
    which does not compile) or `'group'`.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, limit_fet, limit_vds,
                                       limit_pnj, limit_together, vt)

    class _M4(Behavioural):
        instparams = [Parameter(name='VTO', desc='threshold', unit='V',
                                default=0.5),
                      Parameter(name='IS0', desc='current scale', unit='A',
                                default=1e-6),
                      Parameter(name='ISB', desc='bulk saturation current',
                                unit='A', default=1e-15),
                      Parameter(name='N', desc='subthreshold slope', unit='',
                                default=1.3),
                      Parameter(name='VE', desc='saturation voltage',
                                unit='V', default=0.5)]

        @staticmethod
        def analog(d, g, s, b):
            bgs, bds = Branch(g, s), Branch(d, s)
            bbs, bbd = Branch(b, s), Branch(b, d)
            if mode == 'group':
                vgs, vds, vbs, vbd = limit_together(   # noqa
                    limit_fet(bgs.V, VTO), limit_vds(bds.V),
                    limit_pnj(bbs.V, ISB, vt()), limit_pnj(bbd.V, ISB, vt()))
            elif mode == 'probe':
                vgs = limit_fet(bgs.V, VTO)                         # noqa
                vds = limit_vds(bds.V)                              # noqa
                vbs = limit_pnj(bbs.V, ISB, vt())                   # noqa
                vbd = limit_pnj(bbd.V, ISB, vt())                   # noqa
            else:
                vgs, vds, vbs, vbd = bgs.V, bds.V, bbs.V, bbd.V     # noqa
            ids = IS0 * limexp((vgs - VTO) / (N * vt())) \
                * (1.0 - limexp(-vds / VE))                         # noqa
            return (Contribution(bds.I, ids),
                    Contribution(bbs.I, ISB * (limexp(vbs / vt()) - 1.0)),
                    Contribution(bbd.I, ISB * (limexp(vbd / vt()) - 1.0)))
    return _M4


def _cascode4(cls, vdd, vg2, vg1, vb=-0.5):
    """`_cascode` with a fourth terminal, the bulk, held below ground so
    both junctions are reverse biased and the body effect is real."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vg2'] = VS('g2', gnd, v=vg2)
    c['vg1'] = VS('g1', gnd, v=vg1)
    c['vb'] = VS('bulk', gnd, v=vb)
    c['M2'] = cls('vdd', 'g2', 'mid', 'bulk')
    c['M1'] = cls('mid', 'g1', gnd, 'bulk')
    return c


## ------------------------------------------------------------------------
## 1.  The declaration.

def test_limit_together_wraps_the_per_probe_declarations_unchanged():
    """The grouping is an envelope, not a replacement: a model adopts it
    by wrapping the line it already had, and the per-probe forms keep
    producing exactly the `limit_spec` they always did."""
    plain = _fet('both')._hdl_info
    grouped = _fet('group')._hdl_info
    assert [(s[0], s[1]) for s in plain['limit_spec']] == \
           [(s[0], s[1]) for s in grouped['limit_spec']] == \
           [((1, 2), 'fet'), ((0, 2), 'vds')]
    ## The group is carried BESIDE the spec, in declaration order.  It is
    ## not inside the spec entries because every reader of one unpacks it
    ## as a 4-tuple.
    assert plain['limit_groups'] == []
    assert grouped['limit_groups'] == [(False, [0, 1])]
    assert _fet('seq')._hdl_info['limit_groups'] == [(True, [0, 1])]


def test_explain_names_the_group():
    from pycircuit.circuit.hdl import explain
    assert 'fetlim on (g,s) [group 1], limvds on (d,s) [group 1]' \
        in explain(_fet('group'))
    assert 'fetlim on (g,s) [group 1, SPICE order]' in explain(_fet('seq'))
    ## and says nothing at all when nothing is grouped.
    assert '[group' not in explain(_fet('both'))


def test_limit_together_rejects_what_is_not_a_limited_probe():
    from pycircuit.circuit.hdl import (limit_together, limit_fet, limit_vds,
                                       Branch, Node)
    b = Branch(Node('g'), Node('s'))
    with pytest.raises(ValueError, match='at least two probes'):
        limit_together(limit_fet(b.V, 0.5))
    with pytest.raises(ValueError, match='limit_pnj/limit_fet/limit_vds'):
        limit_together(b.V, limit_vds(b.V))
    ## A probe belongs to ONE group -- nesting would give it two, and the
    ## second tag would silently win.
    one, two = limit_together(limit_fet(b.V, 0.5), limit_vds(b.V))
    with pytest.raises(ValueError, match='already in a \\$limit group'):
        limit_together(one, two)


def test_the_four_terminal_mosfet_is_refused_per_probe_and_accepted_grouped():
    """THE CAPABILITY, and the reason this feature is worth having.

    SPICE limits a MOSFET with four probes -- `fetlim(vgs)`,
    `limvds(vds)`, `pnjlim(vbs)`, `pnjlim(vbd)`.  Per-probe, each one
    claims a terminal, and the fourth finds both of its own already
    claimed: there is no single node it can move that does not undo an
    earlier probe, so `generate_code` refuses rather than writing a
    silently wrong one.  **The declaration every real MOSFET wants is not
    expressible at all in the per-probe form.**

    Grouped it compiles, and the triangle it contains -- `(b,s)`, `(b,d)`
    and `(d,s)` -- is resolved at run time by size instead of at compile
    time by position.
    """
    with pytest.raises(ValueError, match='over-determine this device'):
        _mos4('probe')
    ## and the message points at the fix rather than at the roadmap.
    with pytest.raises(ValueError, match='limit_together'):
        _mos4('probe')

    info = _mos4('group')._hdl_info
    assert [(s[0], s[1]) for s in info['limit_spec']] == \
        [((1, 2), 'fet'), ((0, 2), 'vds'), ((3, 2), 'pnj'), ((3, 0), 'pnj')]
    assert info['limit_groups'] == [(False, [0, 1, 2, 3])]


def test_a_group_compiles_through_both_code_generators():
    """Eager and chained are two code generators; a feature that works in
    one is not thereby in the other (that is how `$limit` reached the
    chained path dead, roadmap 12.1)."""
    eager = _fet('group', chained=False)('d', 'g', 's')
    chained = _fet('group', chained=True)('d', 'g', 's')
    assert chained._hdl_info['chained'] and not eager._hdl_info['chained']
    for el in (eager, chained):
        el.update_iparv()
    rng = np.random.default_rng(20260825)
    for _ in range(300):
        x = rng.uniform(-30.0, 30.0, 3)
        x0 = rng.uniform(-30.0, 30.0, 3)
        a, b = eager.limit(x, x0), chained.limit(x, x0)
        assert (a == b).all(), (x, x0, a, b)


## ------------------------------------------------------------------------
## 2.  The write-back, and the properties it has to keep.

def _both_probes_bite(x, x0, vto=0.5):
    """Would each probe of the `(vgs, vds)` star bite, read independently?"""
    vgs, vgs0 = x[1] - x[2], x0[1] - x0[2]
    vds, vds0 = x[0] - x[2], x0[0] - x0[2]
    return (_fetlim(vgs, vgs0, vto, numeric) != vgs
            and _limvds(vds, vds0, numeric) != vds)


def test_every_probe_ends_where_its_own_law_leaves_it():
    """THE INVARIANT THE PER-PROBE FORM CANNOT GIVE, asserted as a number.

    "Each probe carries its own limited value" is checkable without
    re-implementing the write-back: `fetlim` and `limvds` are IDEMPOTENT
    -- applying either to its own output, against the same `vold`, gives
    that output back -- so a branch that its own law would still move is
    a branch the limiter did not actually leave at a limited point.

    The per-probe form fails that.  It applies each correction as a
    DISPLACEMENT of one node, so when a later probe moves a shared
    terminal the earlier probe's branch follows it and lands somewhere its
    law never chose.  Measured over 813 draws in which both probes bite:
    27 of them (3.3%) escape, the worst by **36 V** -- rare because it
    needs both probes to want the same terminal, and large when it
    happens.  The device-level write-back solves the branches as
    CONSTRAINTS, so every one of them ends at a fixed point of its own
    law, on every draw.

    (Idempotence is a property of these two laws, not of limiters in
    general -- `pnjlim` compresses logarithmically and re-compresses its
    own output.  This test is written on the pair that has it, and says
    so rather than generalising.)
    """
    grouped = _fet('group')('d', 'g', 's', VTO=0.5)
    perprobe = _fet('both')('d', 'g', 's', VTO=0.5)
    for el in (grouped, perprobe):
        el.update_iparv()

    rng = np.random.default_rng(20260826)
    checked = {'group': 0, 'both': 0}
    escaped = {'group': 0, 'both': 0}
    worst = 0.0
    for _ in range(2000):
        x = rng.uniform(-60.0, 60.0, 3)
        x0 = rng.uniform(-60.0, 60.0, 3)
        if not _both_probes_bite(x, x0):
            continue
        for name, el in (('group', grouped), ('both', perprobe)):
            out = el.limit(x, x0)
            checked[name] += 1
            vgs, vds = out[1] - out[2], out[0] - out[2]
            ## A TOLERANCE, and it is the one the skill warns about:
            ## the write is `out[v] = out[u] + s`, so the difference read
            ## back is `(u + s) - u`, which is not `s` in floating point.
            ## An ulp of slack is not what this test is trying to catch --
            ## the per-probe gap below is measured in VOLTS.
            slack = 1e-9 * max(1.0, abs(vgs), abs(vds))
            fixed = (abs(_fetlim(vgs, x0[1] - x0[2], 0.5, numeric) - vgs)
                     <= slack
                     and abs(_limvds(vds, x0[0] - x0[2], numeric) - vds)
                     <= slack)
            if not fixed:
                escaped[name] += 1
                if name == 'both':
                    worst = max(
                        worst,
                        abs(_limvds(vds, x0[0] - x0[2], numeric) - vds))
    ## The draws are not vacuous.
    assert checked['group'] > 300, checked
    ## The grouped form: every probe, every time.
    assert escaped['group'] == 0, escaped
    ## The per-probe form: a few per cent of the draws, by VOLTS.
    ## Asserted as a range rather than a count, so it measures the effect
    ## and not the random seed -- but the lower end is what matters, and
    ## it is what fails if the grouped write-back is quietly replaced by
    ## the per-probe one.
    assert 0.01 < escaped['both'] / checked['both'] < 0.10, escaped
    assert worst > 5.0, worst


def test_a_group_is_independent_of_the_order_its_probes_were_declared_in():
    """Order independence, asserted by literally reversing the group.

    The reading order inside a group is CANONICAL -- largest correction
    first, ties by row index -- and the write-back's spanning forest is
    built the same way, so nothing in the result may depend on the order
    the model author wrote the probes in.  `limit_groups` is what the
    generated `limit()` closes over, so reversing its index list is the
    knob; source order is not, since `generate_code` sorts the markers
    into sympy's canonical order before it ever sees them.
    """
    for cls in (_fet('group'), _mos4('group')):
        el = cls(*['d', 'g', 's', 'b'][:len(cls.terminals)])
        el.update_iparv()
        n = len(el.terminals)
        groups = el._hdl_info['limit_groups']
        idx = groups[0][1]
        rng = np.random.default_rng(5)
        for _ in range(200):
            x = rng.uniform(-40.0, 40.0, n)
            x0 = rng.uniform(-40.0, 40.0, n)
            straight = el.limit(x, x0)
            idx.reverse()
            flipped = el.limit(x, x0)
            idx.reverse()
            assert (straight == flipped).all(), (cls, x, x0)
        assert idx == el._hdl_info['limit_groups'][0][1]


def test_a_group_that_did_not_bite_writes_nothing_at_all():
    """`a - (a - b)` is not `b`, so even a no-op write costs the property
    that lets "did limiting fire?" be a convergence signal.  A group is
    written back as a whole, so the rule has to be stated over the whole
    group: if NO probe bit, not one row is touched.

    Both laws here pass a small step through exactly -- `fetlim` any step
    of 0.45 V or less, `limvds` any step that stays inside its box -- so
    a small perturbation must come back bit-identical.
    """
    el = _fet('group')('d', 'g', 's', VTO=0.5)
    el.update_iparv()
    rng = np.random.default_rng(17)
    for _ in range(400):
        x0 = rng.uniform(0.0, 3.0, 3)
        x = x0 + rng.uniform(-0.2, 0.2, 3)
        out = el.limit(x, x0)
        assert (out == x).all(), (x, x0, out)
        assert not np.shares_memory(out, x)


def test_a_group_never_undoes_a_probe_and_stays_finite():
    """The invariant that survived the write-back becoming a run-time
    choice (roadmap 12.1), restated for a group: no row is written twice,
    every row that did not move is untouched EXACTLY, and the result is
    bounded -- the failure this whole area exists to remove was a
    write-back copying an unbounded absolute value across a branch.
    """
    el = _mos4('group')('d', 'g', 's', 'b', VTO=0.5)
    el.update_iparv()
    rng = np.random.default_rng(23)
    moved_hist = set()
    for _ in range(500):
        x = rng.uniform(-50.0, 50.0, 4)
        x0 = rng.uniform(-50.0, 50.0, 4)
        out = el.limit(x, x0)
        movers = [i for i in range(4) if out[i] != x[i]]
        moved_hist.add(len(movers))
        ## At least one row is held: a group anchors on the least drifted
        ## node and derives the rest, so a device can never have every
        ## terminal rewritten.
        assert len(movers) <= 3, (movers, x, x0)
        assert np.all(np.isfinite(out))
        assert np.max(np.abs(out)) <= np.max(np.abs(x)) + 120.0, (x, out)
    ## and the cases are not all trivial.
    assert {0, 2, 3} <= moved_hist, moved_hist


def test_the_write_back_drops_the_smallest_constraint_of_a_cycle():
    """`(b,s)`, `(b,d)` and `(d,s)` close a triangle, and three limited
    values around a loop need not sum to zero -- so one of the three
    cannot be honoured.  This is the case the per-probe form refuses at
    compile time; here it is a run-time trade, and the trade is by SIZE.

    Asserted on `device_writeback` directly, with targets written out by
    hand.  Reading them off an element would mean re-implementing the
    reading pass in the test, which is how a test comes to assert that
    the code does what the code does.
    """
    from pycircuit.circuit._limiting import device_writeback

    ## rows 0..2, all at 10 V, none drifted.  Three constraints round the
    ## loop: 0-1 wants 1, 1-2 wants 2, 0-2 wants 9 -- and 1 + 2 != 9, so
    ## they cannot all hold.
    out = np.array([10.0, 10.0, 10.0])
    drift = np.array([0.0, 5.0, 7.0])
    targets = [(0, 1, 0.0, 1.0),        # correction 1.0  <- the smallest
               (1, 2, 0.0, 2.0),        # correction 2.0
               (0, 2, 0.0, 9.0)]        # correction 9.0
    written = device_writeback(out, targets, drift)
    ## The anchor is row 0, the least drifted, and it is not written.
    assert written == {1, 2} and out[0] == 10.0
    ## The two largest constraints hold exactly ...
    assert out[0] - out[2] == 9.0
    assert out[1] - out[2] == 2.0
    ## ... and the smallest is the one that gave way.
    assert out[0] - out[1] != 1.0

    ## The choice is made on SIZE, not on position: reversing the list
    ## drops the same constraint and produces the same vector.
    out2 = np.array([10.0, 10.0, 10.0])
    device_writeback(out2, targets[::-1], drift)
    assert (out2 == out).all()


def test_the_write_back_holds_the_least_drifted_row():
    """The per-probe rule -- "move the terminal that has drifted further
    from the last accepted point" -- stated for a whole tree: the anchor
    is the least drifted row in the component, and every other row is
    derived from it.  That is what keeps a bounded branch voltage from
    being added to a wild node (roadmap 12.1's 5e48 V gate).
    """
    from pycircuit.circuit._limiting import device_writeback

    targets = [(0, 1, 40.0, 4.0), (1, 2, -17.0, -3.0)]
    for anchor, drift in ((0, np.array([0.0, 9.0, 9.0])),
                          (1, np.array([9.0, 0.0, 9.0])),
                          (2, np.array([9.0, 9.0, 0.0]))):
        out = np.array([50.0, 10.0, 27.0])
        written = device_writeback(out, targets, drift)
        assert anchor not in written
        assert out[anchor] == [50.0, 10.0, 27.0][anchor]
        ## Whichever row is held, BOTH constraints come out exact -- that
        ## is the difference from a per-probe write-back, where the second
        ## probe's displacement drags the first probe's branch with it.
        assert out[0] - out[1] == 4.0 and out[1] - out[2] == -3.0

    ## A row another probe already wrote outranks drift: rewriting it
    ## would undo that probe.
    out = np.array([50.0, 10.0, 27.0])
    device_writeback(out, targets, np.array([0.0, 9.0, 9.0]), pinned={2})
    assert out[2] == 27.0
    assert out[0] - out[1] == 4.0 and out[1] - out[2] == -3.0

    ## TWO pinned rows in one component is the case that has to give
    ## something up: one of them anchors, and the constraint that would
    ## have written the other is dropped rather than undoing it.  (Only
    ## reachable with several groups on one device, which is why it is
    ## asserted here and not through an element.)
    out = np.array([50.0, 10.0, 27.0])
    written = device_writeback(out, targets, np.array([0.0, 9.0, 9.0]),
                               pinned={1, 2})
    assert out[1] == 10.0 and out[2] == 27.0, out
    assert written <= {0}


def test_the_write_back_writes_nothing_when_no_target_bit():
    """Stated at the level it is implemented at, as well as through an
    element: not one row, not even a round trip."""
    from pycircuit.circuit._limiting import device_writeback
    out = np.array([1.0, 2.0, 4.5])
    before = out.copy()
    assert device_writeback(out, [(0, 1, -1.0, -1.0), (1, 2, -2.5, -2.5)],
                            np.array([9.0, 0.0, 3.0])) == set()
    assert (out == before).all()


## ------------------------------------------------------------------------
## 3.  SPICE's ordering -- offered, and measured rather than assumed.

def test_sequential_is_mos1load_s_order_and_the_difference_is_a_number():
    """`sequential=True` reproduces `mos1load.c` exactly: limit `vgs`,
    then recompute `vds` from the UNLIMITED `vgd` before `limvds`.

    `test_limit_fet.py::test_the_vgs_vds_star_is_order_independent_and_is_
    not_spice_s_order` pinned the gap this closes -- 4.0 V against SPICE's
    -0.5 V on this step, because SPICE's `limvds` is handed -30 V and
    clamps it at -0.5.  Both are limited points and both leave the
    solution untouched; they are different places to take the next
    Jacobian.  That test still measures the per-probe form; this one
    measures that the grouped sequential form lands on SPICE's number.
    """
    x0 = np.array([0.2, 0.0, 0.0])
    x = np.array([9.0, 40.0, 0.0])

    vgs, vds, vds0 = x[1] - x[2], x[0] - x[2], x0[0] - x0[2]
    vgs_lim = _fetlim(vgs, x0[1] - x0[2], 0.5, numeric)
    assert vgs_lim == 1.0
    vgd = vgs - vds                          # UNLIMITED, the third branch
    spice = _limvds(vgs_lim - vgd, vds0, numeric)
    ours = _limvds(vds, vds0, numeric)
    assert (spice, ours) == (-0.5, 4.0)

    for which, expect in (('group', ours), ('seq', spice)):
        el = _fet(which)('d', 'g', 's', VTO=0.5)
        el.update_iparv()
        out = el.limit(x, x0)
        assert out[1] - out[2] == vgs_lim, which
        assert out[0] - out[2] == expect, (which, out)


def test_sequential_is_order_dependent_and_that_is_the_point():
    """"In this order" has to mean something.  The independent group is
    order-independent by construction; the sequential one is not, and
    reversing it must actually change the answer -- otherwise the flag
    is decorative."""
    x0 = np.array([0.2, 0.0, 0.0])
    x = np.array([9.0, 40.0, 0.0])
    el = _fet('seq')('d', 'g', 's', VTO=0.5)
    el.update_iparv()
    idx = el._hdl_info['limit_groups'][0][1]
    straight = el.limit(x, x0)
    idx.reverse()
    try:
        flipped = el.limit(x, x0)
    finally:
        idx.reverse()
    assert not (straight == flipped).all(), (straight, flipped)
    ## vds first, then vgs: `limvds` sees the unshifted -- and therefore
    ## unclamped -- 8.8 V, so the drain lands where the independent group
    ## puts it, and it is `vgs` that reads a shifted vector.
    assert straight[0] - straight[2] == -0.5
    assert flipped[0] - flipped[2] == 4.0


## ------------------------------------------------------------------------
## 4.  The DC solves.  A limiter has no correct output, only a better or
##     worse trajectory, so this is the only section that is evidence.

def _count(cls, cond, maxiter=100):
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    c = _cascode(cls, **cond)
    try:
        x, its = _plain_newton(c, maxiter=maxiter)
        return its, float(x[c.get_node_index('mid')])
    except (NoConvergenceError, SingularMatrix):
        return None, None


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_the_grouped_form_costs_about_what_the_per_probe_one_costs():
    """THE ACCEPTANCE MEASUREMENT, and it does not say what roadmap 12.1
    expected it to say.

    At (20 V, 2 V, 0.8 V), on this machine:

        none   NoConvergenceError
        fet    12 iterations
        vds    57
        both   13          <- roadmap 12.1 recorded 30; see the module
        group  13             docstring.  30 was the pre-fix number.
        seq    13

    Over the 48-point grid, summed with a non-convergence counted as 100:

        none  3517  (7 failures)      both   909  (0)
        fet   1222  (1)               group  927  (0)
        vds   1404  (0)               seq   1029  (0)

    So the device-level write-back is not a convergence improvement on
    this circuit: it ties `both` at 34 of the 48 points and costs exactly
    one iteration at the other 14.  Recorded as a number rather than
    argued away, because "a helper that measurably makes things worse is
    a defect report, not a characteristic" -- and here it is neither: it
    is a 2% price for a write-back that is correct where the cheap one is
    approximate, and for a declaration the cheap one cannot express at
    all (see the four-terminal test above).
    """
    classes = {w: _fet(w) for w in ('fet', 'vds', 'both', 'group', 'seq')}
    at_point = {w: _count(c, COND)[0] for w, c in classes.items()}
    assert _count(_fet('none'), COND)[0] is None
    assert at_point['fet'] == 12 and at_point['vds'] == 57
    assert at_point['both'] == at_point['group'] == 13

    total, fails = {}, {}
    for w, cls in classes.items():
        total[w] = fails[w] = 0
        for vdd, vg2, vg1 in GRID:
            its, mid = _count(cls, dict(vdd=vdd, vg2=vg2, vg1=vg1))
            if its is None:
                fails[w] += 1
                total[w] += 100
            else:
                total[w] += its
    ## Nothing grouped may FAIL where the per-probe form succeeds ...
    assert fails['group'] == fails['both'] == 0, fails
    ## ... and the cost is small and bounded, in both directions.  A wide
    ## window on purpose: this is a measurement being held in place, not a
    ## threshold anything was tuned to.
    assert 0.9 < total['group'] / total['both'] < 1.1, total
    ## The one relation that IS structural: two probes beat one that only
    ## watches vds, whichever way they are declared.
    assert total['group'] < total['vds'] and total['both'] < total['vds']


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_the_four_terminal_mosfet_solves_only_because_of_the_group():
    """The intervention against its ABSENCE, on a real solve, on the
    declaration that only exists because of it.

    Same circuit, same solver, same start; the only difference is whether
    the model declares its four probes.  Plain `StandardNewton` -- no
    rescue chain, so nothing else can be doing the work -- and the
    reference comes from `DC` at **gmin = 0**, because the gmin anchor of
    roadmap 12.3 would otherwise be what is being measured.
    """
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix

    bad = dict(vdd=20.0, vg2=1.0, vg1=0.8)
    with pytest.raises((NoConvergenceError, SingularMatrix)):
        _plain_newton(_cascode4(_mos4('none'), **bad))

    c = _cascode4(_mos4('group'), **bad)
    x, its = _plain_newton(c)
    assert its <= 25, its
    mid = x[c.get_node_index('mid')]
    assert 0.0 < mid < bad['vdd'] - 0.5

    ## The SAME point the rescue chain reaches with no limiter at all: a
    ## limiter moves the path, never the answer.
    ref = DC(_cascode4(_mos4('none'), **bad), toolkit=numeric,
             gmin=0.0).solve()
    assert_allclose(mid, float(ref.v('mid')), rtol=1e-6)


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_the_four_terminal_mosfet_is_cheaper_with_the_group_across_the_grid():
    """The single point above could be a lucky one.  Across the same
    48-point grid, with a non-convergence counted as 100:

        none   3478 Jacobian evaluations, 6 of the 48 not converging
        group   896                     , none

    and every point that both reach lands on the same operating point.
    """
    none_cls, group_cls = _mos4('none'), _mos4('group')
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix

    total = {'none': 0, 'group': 0}
    fails = {'none': 0, 'group': 0}
    for vdd, vg2, vg1 in GRID:
        mids = {}
        for name, cls in (('none', none_cls), ('group', group_cls)):
            c = _cascode4(cls, vdd, vg2, vg1)
            try:
                x, its = _plain_newton(c)
                total[name] += its
                mids[name] = float(x[c.get_node_index('mid')])
            except (NoConvergenceError, SingularMatrix):
                fails[name] += 1
                total[name] += 100
        if len(mids) == 2:
            ## Same answer, limited or not -- the assertion that catches a
            ## limiter that moved the solution.
            assert_allclose(mids['group'], mids['none'], rtol=1e-5,
                            atol=1e-6)
    assert fails['group'] == 0 and fails['none'] >= 4, fails
    assert total['group'] < total['none'] / 3.0, total


## ------------------------------------------------------------------------
## 5.  BETWEEN devices: the inter-device clash (roadmap 15, Stage 0) and
##     the circuit-level forest that was priced against it (the spike).

def _diffpair(cls, vin, m2_first=False, vdd=5.0):
    """Differential pair: ideal 200 uA tail sink, 5 kOhm loads, inputs
    `2.5 +- vin`.  `m2_first` flips the INSERTION order, which is the
    order `SubCircuit.limit` visits the two devices in."""
    from pycircuit.circuit.elements import IS, R
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=vdd)
    c['vp'] = VS('inp', gnd, v=2.5 + vin)
    c['vn'] = VS('inn', gnd, v=2.5 - vin)
    c['rl1'] = R('vdd', 'd1', r=5e3)
    c['rl2'] = R('vdd', 'd2', r=5e3)
    c['itail'] = IS('tail', gnd, i=200e-6)
    devs = [('M1', ('d1', 'inp', 'tail')), ('M2', ('d2', 'inn', 'tail'))]
    for name, nodes in (devs[::-1] if m2_first else devs):
        c[name] = cls(*nodes)
    return c


def _diffpair_row(cls, vin):
    """`[M1 first, M2 first]` Jacobian counts and tail voltages."""
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    row = []
    for m2_first in (False, True):
        c = _diffpair(cls, vin, m2_first)
        try:
            x, its = _plain_newton(c)
            row.append((its, float(x[c.get_node_index('tail')])))
        except (NoConvergenceError, SingularMatrix):
            row.append((None, None))
    return row


DIFF_VIN = (-1.0, -0.3, 0.0, 0.3, 1.0)
DIFF_TAIL = {1.0: 2.818678, 0.3: 2.121207, 0.0: 1.845121}


@pytest.fixture
def _circuit_level(request):
    """Switch `SubCircuit.limit` into the circuit-level forest for one
    test and back out again, whatever happens."""
    from pycircuit.circuit import _limiting
    old = _limiting.CIRCUIT_LEVEL
    _limiting.CIRCUIT_LEVEL = request.param
    yield request.param
    _limiting.CIRCUIT_LEVEL = old


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_two_devices_on_one_tail_are_resolved_by_instance_order():
    """THE DOCUMENTED CLASH (roadmap 15, Stage 0), pinned as it stands.

    Same equations, same `x` layout, both devices limited; only the order
    `SubCircuit.limit` visits them differs.  Measured on this machine,
    `[M1 first, M2 first]`:

        vin        -1.0        -0.3        0        0.3       1.0
        both     [45, 14]   [26, 14]  [15, 15]  [14, 26]  [14, FAIL]

    Antisymmetric in `vin`, as a clash must be.  The mechanism, traced at
    iteration 0 for `vin = 1.0` from a zero start (Newton proposes tail =
    -1.45e8): M1's `fetlim` alone wants `tail = 2.5` (M1 is the device
    that is ON) and M2's alone wants `tail = 0.5`.  Visited FIRST, M2
    writes `tail = 0.5`; M1 then reads a tail that has "drifted" only
    0.5 V against a gate that drifted 3.5 V from zero, and its own
    move-the-wilder-node rule moves the pinned GATE to 1.5 V.  Newton
    restores the gate, M1 sits at `vgs = 3 V` on a subthreshold
    exponential, and the solve falls over.  Visited first, M1 writes
    `tail = 2.5`, M2 then sees `vgs = -1 V` and does not bite at all.

    `limit_together` cannot help: it is per DEVICE and this is between
    devices.  `DC()`'s ladder rescues both orders to the same point,
    which is why this runs on plain Newton.
    """
    from pycircuit.circuit.dcanalysis import DC
    cls = _fet('both')
    rows = {vin: _diffpair_row(cls, vin) for vin in DIFF_VIN}

    ## Every order that converges lands on the ladder's answer.
    ref = DC(_diffpair(_fet('none'), 1.0), toolkit=numeric).solve()
    assert_allclose(float(ref.v('tail')), DIFF_TAIL[1.0], rtol=1e-6)
    for vin, row in rows.items():
        for its, tail in row:
            if its is not None:
                assert_allclose(tail, DIFF_TAIL[abs(vin)], rtol=1e-5)

    ## The clash, as measured.  At vin = +1 the losing order fails
    ## outright; at vin = -1 it converges, at 45 -- not a mirror image,
    ## because the two exponentials are not evaluated at mirrored points
    ## (the tail is 2.5 or 0.5 whoever wrote it).  Pinned as measured.
    assert rows[1.0][1][0] is None, rows
    assert rows[-1.0][0][0] is not None and rows[-1.0][0][0] > 30, rows
    ## The winning order converges cheaply.
    assert rows[1.0][0][0] <= 20 and rows[-1.0][1][0] <= 20, rows
    ## ... and at |vin| = 0.3 the losing order costs a wide margin more.
    assert rows[0.3][1][0] > 1.5 * rows[0.3][0][0], rows
    assert rows[-0.3][0][0] > 1.5 * rows[-0.3][1][0], rows
    ## Balanced, there is nothing to fight over.
    assert rows[0.0][0][0] == rows[0.0][1][0] <= 20, rows


@pytest.mark.filterwarnings('ignore:overflow encountered')
@pytest.mark.parametrize('_circuit_level', [True, 'pinned'], indirect=True)
def test_a_circuit_level_forest_is_order_independent_and_worse(
        _circuit_level):
    """THE SPIKE'S RESULT, and it is a refusal.

    Lifting `device_writeback` to `SubCircuit.limit` -- every element's
    probe targets, ONE spanning forest over the whole circuit -- makes
    the diff pair order-INDEPENDENT and does not make it converge:

        vin          -1.0       -0.3        0        0.3       1.0
        forest     [F, F]    [43, 43]  [15, 15]  [43, 43]   [F, F]
        +pinned    [F, F]    [32, 32]  [15, 15]  [32, 32]   [F, F]

    and over the 48-point cascode grid (`both`, 909 with no failure as
    the baseline) it costs 1878 with 8 failures (forest) and 1052 with
    none (pinned anchors).  The only thing it cures is the Level-1
    cascode's bulk tug-of-war from a uniform 20 V start, and the pinned
    variant loses that again at 40 V.

    WHY, traced at iteration 0, `vin = 1.0`: two devices asking one node
    for 2.5 V and 0.5 V is a CONFLICT, not a cycle.  M1's `(inp, tail)`
    and M2's `(inn, tail)` are two edges of a star, so the forest keeps
    BOTH -- and the only way to honour both is to anchor on the least-
    drifted node (`inn`) and derive `tail = 0.5` and then `inp = 1.5`: a
    source-pinned gate moved by 2 V, which is precisely the failure the
    runtime write-back was built to remove, now reached identically in
    both orders.  Holding source-pinned nodes as anchors (the `'pinned'`
    variant) turns the conflict into a dropped edge, but which edge is
    dropped is decided by anchor choice -- least drift, i.e. `inn` at
    1.5 V against `inp` at 3.5 V from a zero start -- and that hands the
    tail to M2, the OFF device, again.  Nothing in `|x - x0|` knows that
    M1 is the device that is on; the correction sizes differ by 2 V in
    1.45e8 and the tie is decided by an artefact of the start vector.

    So a circuit-level forest cannot resolve a genuine conflict; it can
    only choose a loser without the information to choose the right one.
    That is the number that decides whether vector PCNR gets built.  The
    switch stays OFF; this test keeps the measurement honest.
    """
    cls = _fet('both')
    rows = {vin: _diffpair_row(cls, vin) for vin in DIFF_VIN}
    ## Order-independent: both orders do the same thing, at every point.
    for vin, row in rows.items():
        assert row[0][0] == row[1][0], (vin, row)
    ## ... and the thing they do is fail where the better order succeeded.
    assert rows[1.0][0][0] is None and rows[-1.0][0][0] is None, rows
    assert rows[0.3][0][0] > 25 and rows[-0.3][0][0] > 25, rows
    assert rows[0.0][0][0] <= 20, rows
    for vin, row in rows.items():
        for its, tail in row:
            if its is not None:
                assert_allclose(tail, DIFF_TAIL[abs(vin)], rtol=1e-5)


## --- Stage 2: the same clash under vector PCNR ---------------------------

def _pcnr_solve(c, x0=None, maxiter=100):
    """`pcnr.solve_dc` -- no ladder, no anchor, no node limiting.  NOT
    `DC(pcnr=True)`: that entry gates on `pcnr_junctions()`, the
    pnj-only gmin pair view, so a circuit of `fetlim`/`limvds` devices
    and no junction falls through to the ordinary solver there."""
    from pycircuit.circuit import pcnr as P
    try:
        x, _v, its = P.solve_dc(c, gnd, x0=x0, reltol=1e-4, abstol=1e-6,
                                maxiter=maxiter)
    except (RuntimeError, np.linalg.LinAlgError, FloatingPointError,
            OverflowError, ValueError):
        return None, None
    return its, x


def _diffpair_row_pcnr(cls, vin):
    row = []
    for m2_first in (False, True):
        c = _diffpair(cls, vin, m2_first)
        its, x = _pcnr_solve(c)
        row.append((its, None if x is None
                    else float(x[c.get_node_index('tail')])))
    return row


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_the_diff_pair_under_pcnr_is_order_independent_and_converges():
    """THE NUMBER THE VECTOR-PCNR EFFORT EXISTS FOR (roadmap 15, Stage 2).
    Same circuit as `test_two_devices_on_one_tail_are_resolved_by_
    instance_order`, `[M1 first, M2 first]`, Jacobian evaluations::

        vin          -1.0       -0.3        0        0.3       1.0
        plain      [45, 14]   [26, 14]  [15, 15]  [14, 26]  [14, FAIL]
        PCNR       [22, 22]   [23, 23]  [15, 15]  [23, 23]  [22, 22]

    Two claims, asserted separately:

    1. **Order-independent**: equal counts in both orders at every
       `vin`, and a symmetric table in `vin` -- there is no shared node
       for the two `fetlim`s to fight over, by construction.
    2. **Converges** everywhere plain Newton did, AND at `vin = +1.0`
       where the losing order failed.  Every tail equals `DC()`'s.

    The price: 22 against the winning order's 14 -- the winning order
    got lucky (M1 wrote the tail first), PCNR pays for not depending on
    luck.  Windows are wide; the shape is what is pinned.
    """
    cls = _fet('both')
    rows = {vin: _diffpair_row_pcnr(cls, vin) for vin in DIFF_VIN}
    for vin, row in rows.items():
        assert row[0][0] == row[1][0], (vin, row)                  # claim 1
        assert row[0][0] is not None and row[0][0] <= 40, (vin, row)  # 2
        assert_allclose(row[0][1], DIFF_TAIL[abs(vin)], rtol=1e-5)
        assert_allclose(row[1][1], DIFF_TAIL[abs(vin)], rtol=1e-5)
    assert rows[1.0][0][0] == rows[-1.0][0][0], rows
    assert rows[0.3][0][0] == rows[-0.3][0][0], rows
    assert rows[0.0][0][0] <= 20, rows


@pytest.mark.filterwarnings('ignore:overflow encountered')
def test_the_level_1_cascode_from_20_v_under_pcnr():
    """Stage 0's bulk tug-of-war (`[17,10]`, `[F,10]`, `[F,12]`, `[F,21]`
    at the time; re-measured on plain Newton below) under PCNR, both
    orders, uniform 20 V start::

        (vdd, vg2, vg1[, vb])   plain            PCNR
        (20, 2, 0.8)            [13, 8]          [12, 12]
        (40, 2, 0.8)            [17, 8]          [12, 12]
        (5, 2, 0.8)             [31, 9]          [19, 19]
        (20, 2, 0.8, vb = -2)   [13, 8]          [FAIL, FAIL]

    Order-independent at every point; converges at three; fails in
    BOTH orders at the fourth.  Traced: from the all-off start the
    middle node is held by nothing, the first Newton step proposes
    `mid = 4.9e9`, PCNR limits no node, and by the fourth iteration
    `mid` is 3e117 and the matrix singular.  The ordinary path's
    limiter writes the bounded branch voltage INTO the node; PCNR has
    no such write, which is Stage 1's second obstacle on a FET.  The
    plain Newton column differs from Stage 0's record because `von`
    now follows the bulk bias (roadmap 12.6(c)), which is what Stage 0
    measured without.
    """
    from pycircuit.circuit.elements_hdl import MosLevel1Hdl

    def mk(vdd, vg2, vg1, vb=-0.5, m2_first=False):
        c = SubCircuit()
        c['vdd'] = VS('vdd', gnd, v=vdd)
        c['vg2'] = VS('g2', gnd, v=vg2)
        c['vg1'] = VS('g1', gnd, v=vg1)
        c['vb'] = VS('bulk', gnd, v=vb)
        devs = [('M2', ('vdd', 'g2', 'mid', 'bulk')),
                ('M1', ('mid', 'g1', gnd, 'bulk'))]
        for nm, nodes in (devs[::-1] if m2_first else devs):
            c[nm] = MosLevel1Hdl(*nodes)
        return c

    table = {}
    for pt in ((20.0, 2.0, 0.8), (40.0, 2.0, 0.8), (5.0, 2.0, 0.8),
               (20.0, 2.0, 0.8, -2.0)):
        row = []
        for m2 in (False, True):
            c = mk(*pt, m2_first=m2)
            x0 = np.full(c.n, 20.0)
            x0[c.get_node_index(gnd)] = 0.0
            its, x = _pcnr_solve(c, x0)
            row.append((its, None if x is None
                        else float(x[c.get_node_index('mid')])))
        table[pt] = row
    for pt, row in table.items():
        assert row[0][0] == row[1][0], (pt, row)          # order-independent
    for pt in ((20.0, 2.0, 0.8), (40.0, 2.0, 0.8), (5.0, 2.0, 0.8)):
        assert table[pt][0][0] is not None and table[pt][0][0] <= 30, table
        assert_allclose(table[pt][0][1], table[pt][1][1], rtol=1e-6)
    assert table[(20.0, 2.0, 0.8, -2.0)][0][0] is None, table


@pytest.mark.filterwarnings('ignore:overflow encountered')
@pytest.mark.parametrize('cls_name', ['fet-both', 'mos4-group'])
def test_the_grid_under_pcnr_converges_everywhere_plain_newton_did(cls_name):
    """The 48-point cascode grid, PCNR against plain Newton, measured::

        variant       plain          PCNR (maxiter 100)   PCNR (200)
        fet both      909 / 0 fail   1281 / 0 fail        1281 / 0
        mos4 group    896 / 0 fail   1278 / 1 fail        1430 / 0
                                     (20, 2, 1.2) at 152

    PCNR costs ~40% more over the grid and every answer agrees on the
    nodes.  The one point over 100 is a `limvds` limit cycle (`vds`
    21 -> clamp 2 -> 4 -> 14 -> 21 -> ...) riding Newton's `N*VT`-per-
    step creep down a 1e25 A subthreshold exponential; it is slow, not
    divergent.  The window is wide; what is pinned is "nothing that
    converged stops converging" at a 200 budget, and the answers.
    """
    from pycircuit.circuit.tests.test_limit_fet import _cascode
    if cls_name == 'fet-both':
        cls, mk = _fet('both'), _cascode
    else:
        cls, mk = _mos4('group'), _cascode4
    tot_p, tot_q = 0, 0
    for cond in GRID:
        ip, xp = _count_x(cls, cond, mk)
        c = mk(cls, *cond)
        iq, xq = _pcnr_solve(c, maxiter=200)
        assert ip is not None and iq is not None, (cond, ip, iq)
        nn = len(c.nodes)
        assert_allclose(xq[:nn], xp[:nn], rtol=0, atol=1e-3)
        tot_p += ip
        tot_q += iq
    assert 800 < tot_p < 1000, tot_p
    assert 1000 < tot_q < 1700, tot_q


def _count_x(cls, cond, mk):
    from pycircuit.circuit.nrsolver import NoConvergenceError
    from pycircuit.circuit.analysis import SingularMatrix
    c = mk(cls, *cond)
    try:
        x, its = _plain_newton(c, maxiter=100)
        return its, x
    except (NoConvergenceError, SingularMatrix):
        return None, None
