"""Optional series resistances: the branch formulation, and its limit.

PSP103 writes each of its seven parasitic resistances as

    if ((R) > 0.0) I(N1,N2) <+ (G) * V(N1,N2);
    else           V(N1,N2) <+ 0.0;

Read as topology that is a node collapse decided by a parameter VALUE,
and the element's unknown count depends on its parameters.  Read as
algebra, both arms are one constitutive relation,

    V(N1,N2) <+ I(N1,N2) * R

-- the resistor for R > 0, an exact short for R = 0, because the MNA row
is ``-(v1 - v2) + i_br*R = 0``.  These tests pin that the second reading
holds numerically, including at exactly zero, and pin the reason it is
not the whole answer: it costs an internal node and a branch current per
resistance whether or not the resistance exists.
"""
import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import R, SubCircuit, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.hdl import (Behavioural, Branch, Collapse,
                                   Contribution, Node)
from pycircuit.utilities.param import Parameter


class SeriesR(Behavioural):
    """A resistance written as a branch: ``V = I*R``."""
    instparams = [Parameter(name='r', desc='R', unit='ohm', default=0.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        return Contribution(b.V, b.I * r)                      # noqa: F821


class SeriesG(Behavioural):
    """The same written as a conductance: ``I = V/R``."""
    instparams = [Parameter(name='r', desc='R', unit='ohm', default=1.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        return Contribution(b.I, b.V / r)                      # noqa: F821


def _divider(cls, rser, rload=1e3, vin=1.0):
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['vs'] = VS(na, gnd, v=vin)
    c['rs'] = cls(na, nb, r=rser)
    c['rl'] = R(nb, gnd, r=rload)
    c.update_iparv()
    return float(DC(c, toolkit=numeric).solve().v(nb, gnd))


class TestTheBranchFormCoversBothArms(object):

    @pytest.mark.parametrize('rser', [0.0, 1e-12, 1e-6, 1e-3, 1.0, 100.0,
                                      1e3, 1e6, 1e9])
    def test_exact_at_every_resistance_including_zero(self, rser):
        got = _divider(SeriesR, rser)
        want = 1e3 / (1e3 + rser)
        assert got == pytest.approx(want, rel=1e-12)

    def test_zero_really_is_a_short_not_a_small_resistance(self):
        """Bit-exact, because the row degenerates to `v1 = v2`."""
        assert _divider(SeriesR, 0.0) == 1.0

    def test_the_conductance_form_cannot_express_the_zero_arm(self):
        """Why the model needs a conditional in the first place."""
        with pytest.raises(Exception):
            _divider(SeriesG, 0.0)

    def test_it_matches_a_hand_written_resistor(self):
        for rser in (1.0, 47.0, 1e4):
            c = SubCircuit()
            na, nb = c.add_node('a'), c.add_node('b')
            c['vs'] = VS(na, gnd, v=1.0)
            c['rs'] = R(na, nb, r=rser)
            c['rl'] = R(nb, gnd, r=1e3)
            c.update_iparv()
            hand = float(DC(c, toolkit=numeric).solve().v(nb, gnd))
            assert _divider(SeriesR, rser) == pytest.approx(hand, rel=1e-12)

    def test_it_works_in_a_transient_too(self):
        """A short must stay a short once a capacitor is in the loop."""
        from pycircuit.circuit.elements import C as Cap

        got = []
        for rser in (0.0, 1e-9):
            c = SubCircuit()
            na, nb = c.add_node('a'), c.add_node('b')
            c['vs'] = VS(na, gnd, v=1.0)
            c['rs'] = SeriesR(na, nb, r=rser)
            c['rl'] = R(nb, gnd, r=1e3)
            c['cl'] = Cap(nb, gnd, c=1e-9)
            c.update_iparv()
            res = Transient(c, toolkit=numeric).solve(tend=1e-6,
                                                      timestep=1e-8)
            got.append(np.asarray(res.v(nb, gnd).y, float)[-1])
        assert got[0] == pytest.approx(1.0, abs=1e-9)
        assert got[1] == pytest.approx(got[0], abs=1e-9)


class TestTheCostOfNotCollapsing(object):
    """The reason the branch form is a workaround and not the answer.

    Quantified in `benchmarks/collapse_cost.py`: a chain of 100 devices
    each carrying PSP's seven parasitics goes from 103 unknowns to 1503,
    and from 7.0 ms to 100.6 ms per DC solve -- 14x, for the same
    answer.
    """

    def test_a_zero_resistance_still_costs_two_unknowns(self):
        e = SeriesR(Node('p'), Node('m'), r=0.0)
        e.update_iparv()
        ## one branch current; the node pair is NOT merged
        assert len(e.branches) == 1
        assert e.n == 3

    def test_whereas_an_unconditional_collapse_costs_none(self):
        """`V(a,b) <+ 0` written literally IS collapsed, today.

        The gap is only that the decision cannot depend on a parameter:
        the compiler collapses what it can see is zero, and `r` is a
        symbol at compile time even when the instance passes 0.
        """
        class Shorted(Behavioural):
            instparams = []

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                return (Contribution(Branch(p, mid).V, 0.0),
                        Contribution(Branch(mid, m).I, 1e-3
                                     * Branch(mid, m).V))

        e = Shorted(Node('p'), Node('m'))
        e.update_iparv()
        assert len(e.branches) == 0
        assert e.n == 2, 'the internal node should have been absorbed'

    def test_the_collapse_still_produces_the_right_conductance(self):
        class Shorted(Behavioural):
            instparams = []

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                return (Contribution(Branch(p, mid).V, 0.0),
                        Contribution(Branch(mid, m).I, 1e-3
                                     * Branch(mid, m).V))

        e = Shorted(Node('p'), Node('m'))
        e.update_iparv()
        G = np.asarray(e.G(np.array([1.0, 0.0])), float)
        assert np.allclose(G, 1e-3 * np.array([[1., -1.], [-1., 1.]]))


class TestCollapseKeepsBranchNames(object):
    """A collapse renames nodes; it must not merge named branches.

    The rewrite rebuilt every touched `Branch` from its node pair alone,
    which silently dropped the name -- so two named branches that a
    collapse brought onto the same node pair became one.
    """

    def test_two_named_branches_survive_a_collapse(self):
        class Coll(Behavioural):
            instparams = [Parameter(name='r1', desc='R1', unit='ohm',
                                    default=1e3),
                          Parameter(name='r2', desc='R2', unit='ohm',
                                    default=2e3)]

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                ## `mid` is absorbed into `p`, bringing both named
                ## branches onto the (p, m) pair.
                return (Contribution(Branch(p, mid).V, 0.0),
                        Contribution(Branch(mid, m, 'a').V,
                                     Branch(mid, m, 'a').I * r1),   # noqa
                        Contribution(Branch(mid, m, 'b').V,
                                     Branch(mid, m, 'b').I * r2))   # noqa

        e = Coll(Node('p'), Node('m'), r1=1e3, r2=2e3)
        e.update_iparv()
        assert len(e.branches) == 2, e.branches
        assert {b.name for b in e.branches} == {'a', 'b'}

        ## and they behave as two resistors in parallel
        c = SubCircuit()
        na = c.add_node('a')
        c['vs'] = VS(na, gnd, v=1.0)
        c['e'] = Coll(na, gnd, r1=1e3, r2=2e3)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        i_vs = float(res.i('vs.plus'))
        assert -i_vs == pytest.approx(1.0 / 1e3 + 1.0 / 2e3, rel=1e-9)


class TestParameterDrivenCollapse(object):
    """`Collapse(branch, when)` -- PSP's CollapsableR, expressed.

    The condition may mention parameters only, so the topology is fixed
    the moment an instance is built.  Each distinct combination is
    compiled once and cached; instances are retargeted to theirs.
    """

    class Parasitic(Behavioural):
        """A conductance behind an optional series resistance.

        Written exactly as the model does: the resistor arm divides by
        `rd`, which is only safe because the collapsed variant never
        compiles that contribution at all.
        """
        instparams = [Parameter(name='rd', desc='series R', unit='ohm',
                                default=0.0),
                      Parameter(name='gcore', desc='core G', unit='S',
                                default=1e-3)]

        @staticmethod
        def analog(p, m):
            di = Node('di')
            br = Branch(p, di, 'rd')
            core = Branch(di, m)
            return (Contribution(core.I, gcore * core.V),        # noqa: F821
                    Contribution(br.I, br.V / rd),               # noqa: F821
                    Collapse(br, rd <= 0))                       # noqa: F821

    def test_the_collapsed_instance_costs_nothing(self):
        e = self.Parasitic(Node('p'), Node('m'), rd=0.0, gcore=1e-3)
        e.update_iparv()
        assert e.n == 2, 'internal node should be absorbed'
        assert len(e.branches) == 0

    def test_the_uncollapsed_instance_keeps_the_node(self):
        e = self.Parasitic(Node('p'), Node('m'), rd=50.0, gcore=1e-3)
        e.update_iparv()
        assert e.n == 3, 'internal node should survive'

    def test_a_zero_resistance_never_divides_by_zero(self):
        """The collapsed variant does not compile `V/rd` at all."""
        import warnings
        e = self.Parasitic(Node('p'), Node('m'), rd=0.0, gcore=1e-3)
        e.update_iparv()
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            G = np.asarray(e.G(np.array([1.0, 0.0])), float)
        assert np.allclose(G, 1e-3 * np.array([[1., -1.], [-1., 1.]]))

    @pytest.mark.parametrize('rd', [0.0, 1e-9, 1.0, 50.0, 1e4])
    def test_both_variants_give_the_series_combination(self, rd):
        c = SubCircuit()
        na = c.add_node('a')
        c['vs'] = VS(na, gnd, v=1.0)
        c['e'] = self.Parasitic(na, gnd, rd=rd, gcore=1e-3)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        want = -1.0 / (rd + 1e3)
        assert float(res.i('vs.plus')) == pytest.approx(want, rel=1e-9)

    def test_variants_are_compiled_once_and_shared(self):
        """One class per mask, not one per instance."""
        base = self.Parasitic
        a = base(Node('p'), Node('m'), rd=0.0)
        b = base(Node('p'), Node('m'), rd=0.0)
        c = base(Node('p'), Node('m'), rd=5.0)
        assert type(a) is type(b), 'same mask should share a class'
        assert type(a) is not type(c)
        ## rd > 0 is the class's own all-False mask, so it is not a
        ## variant at all -- only the collapsed arm needed building.
        assert type(c) is base
        assert type(a) is not base
        assert base.__dict__['_hdl_mask_classes'][(True,)] is type(a)

    def test_both_variants_coexist_in_one_circuit(self):
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=1.0)
        c['e0'] = self.Parasitic(na, nb, rd=0.0, gcore=1e-3)
        c['e1'] = self.Parasitic(nb, gnd, rd=250.0, gcore=1e-3)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        want = -1.0 / (1e3 + 250.0 + 1e3)
        assert float(res.i('vs.plus')) == pytest.approx(want, rel=1e-9)

    def test_moving_a_gating_parameter_afterwards_is_refused(self):
        """It would change `n` behind the circuit's node map."""
        e = self.Parasitic(Node('p'), Node('m'), rd=0.0, gcore=1e-3)
        e.update_iparv()
        with pytest.raises(ValueError, match='gates a node collapse'):
            e.ipar.rd = 50.0
            e.update_iparv()

    def test_a_non_gating_parameter_may_still_move(self):
        e = self.Parasitic(Node('p'), Node('m'), rd=0.0, gcore=1e-3)
        e.update_iparv()
        e.ipar.gcore = 4e-3
        e.update_iparv()
        G = np.asarray(e.G(np.array([1.0, 0.0])), float)
        assert G[0, 0] == pytest.approx(4e-3, rel=1e-12)

    def test_an_operating_point_condition_is_refused(self):
        """A collapse that moves per iteration is a different feature."""
        with pytest.raises(ValueError, match='parameters only'):
            class Bad(Behavioural):
                instparams = []

                @staticmethod
                def analog(p, m):
                    b = Branch(p, m, 'x')
                    return (Contribution(b.I, 1e-3 * b.V),
                            Collapse(b, b.V > 0))


class TestManyCollapsesTogether(object):
    """PSP has seven of them; the mask is a combination, not a flag."""

    class Seven(Behavioural):
        instparams = ([Parameter(name='r%d' % k, desc='R%d' % k,
                                 unit='ohm', default=0.0)
                       for k in range(3)]
                      + [Parameter(name='gcore', desc='G', unit='S',
                                   default=1e-3)])

        @staticmethod
        def analog(p, m):
            rs = (r0, r1, r2)                                    # noqa: F821
            node, out = p, []
            for k in range(3):
                nxt = Node('i%d' % k)
                br = Branch(node, nxt, 'r%d' % k)
                out.append(Contribution(br.I, br.V / rs[k]))
                out.append(Collapse(br, rs[k] <= 0))
                node = nxt
            core = Branch(node, m)
            out.append(Contribution(core.I, gcore * core.V))      # noqa: F821
            return tuple(out)

    @pytest.mark.parametrize('rvals,nodes_kept', [
        ((0.0, 0.0, 0.0), 0),
        ((10.0, 0.0, 0.0), 1),
        ((0.0, 20.0, 0.0), 1),
        ((10.0, 20.0, 0.0), 2),
        ((10.0, 20.0, 30.0), 3),
    ])
    def test_the_mask_selects_the_right_topology(self, rvals, nodes_kept):
        e = self.Seven(Node('p'), Node('m'), gcore=1e-3,
                       **{'r%d' % k: v for k, v in enumerate(rvals)})
        e.update_iparv()
        assert e.n == 2 + nodes_kept

    @pytest.mark.parametrize('rvals', [
        (0.0, 0.0, 0.0), (10.0, 0.0, 0.0), (0.0, 20.0, 0.0),
        (10.0, 20.0, 0.0), (10.0, 20.0, 30.0), (0.0, 0.0, 30.0)])
    def test_every_topology_gives_the_series_sum(self, rvals):
        c = SubCircuit()
        na = c.add_node('a')
        c['vs'] = VS(na, gnd, v=1.0)
        c['e'] = self.Seven(na, gnd, gcore=1e-3,
                            **{'r%d' % k: v for k, v in enumerate(rvals)})
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        want = -1.0 / (sum(rvals) + 1e3)
        assert float(res.i('vs.plus')) == pytest.approx(want, rel=1e-9)
