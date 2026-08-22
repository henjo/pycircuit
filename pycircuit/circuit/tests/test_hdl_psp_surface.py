"""Language surface PSP103 needs beyond what the DSL already had.

Two items, both found by reading the shipped PSP 103.8.2 source in the
IHP Open PDK rather than by reading the LRM:

* **named branches.** ``branch (a,b) br1, br2;`` declares two DISTINCT
  branches across one node pair, each with its own current. The compiler
  keyed branches on the node pair alone, so the second declaration was
  silently merged into the first and its constitutive relation vanished.

* **`$simparam`.** PSP calls it exactly once,
  ``$simparam("gmin", 0.0)`` at `PSP103_module.include:642`.

`$mfactor` is deliberately absent -- see `TestMfactorIsNotNeeded`.
"""
import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import R, SubCircuit, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Node, simparam, var)
from pycircuit.utilities.param import Parameter


def _mk(cls, *nodes, **kw):
    e = cls(*[Node(n) for n in nodes], **kw)
    e.update_iparv()
    return e


class TestNamedBranches(object):

    def test_two_named_branches_get_two_currents(self):
        """Both constitutive relations must survive.

        Two voltage sources in series between the same pair of terminals
        is the smallest thing that needs this: keyed on the node pair
        alone, one of them disappears.
        """
        class TwoV(Behavioural):
            instparams = [Parameter(name='v1', desc='V1', unit='V',
                                    default=1.0),
                          Parameter(name='v2', desc='V2', unit='V',
                                    default=2.0)]

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                return (Contribution(Branch(p, mid, 'a').V, v1),   # noqa
                        Contribution(Branch(mid, m, 'b').V, v2))   # noqa

        e = _mk(TwoV, 'p', 'm', v1=1.0, v2=2.0)
        assert len(e.branches) == 2
        assert {b.name for b in e.branches} == {'a', 'b'}

        c = SubCircuit()
        np_ = c.add_node('p')
        c['E'] = TwoV(np_, gnd, v1=1.0, v2=2.0)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        assert float(res.v(np_, gnd)) == pytest.approx(3.0, abs=1e-9)

    def test_an_unnamed_branch_keeps_its_old_identity(self):
        """Two anonymous `Branch(a,b)` objects are still ONE branch.

        This is what every existing element relies on -- a model writes
        `Branch(p, m)` in two places meaning the same branch -- so naming
        had to be additive.  Written as a V-contribution plus a probe of
        the same branch: contributing V and I to one branch would be a
        switch branch, which is refused for unrelated reasons.
        """
        class OneV(Behavioural):
            instparams = [Parameter(name='v', desc='V', unit='V',
                                    default=1.0)]

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                return (Contribution(Branch(p, mid).V, v),         # noqa
                        Contribution(Branch(mid, m).I,
                                     Branch(p, mid).I))

        e = _mk(OneV, 'p', 'm', v=1.5)
        assert len(e.branches) == 1
        assert e.branches[0].name is None

    def test_named_and_unnamed_across_the_same_pair_are_distinct(self):
        class Mixed(Behavioural):
            instparams = []

            @staticmethod
            def analog(p, m):
                return (Contribution(Branch(p, m).V, 1.0),
                        Contribution(Branch(p, m, 'x').V, 2.0))

        e = _mk(Mixed, 'p', 'm')
        assert len(e.branches) == 2
        assert sorted(b.name or '' for b in e.branches) == ['', 'x']

    def test_the_switch_branch_refusal_still_names_the_branch(self):
        """The refusal message had to learn the new key shape."""
        with pytest.raises(NotImplementedError, match='switch branch'):
            class Sw(Behavioural):
                instparams = []

                @staticmethod
                def analog(p, m):
                    b = Branch(p, m, 'ch')
                    return (Contribution(b.V, 1.0),
                            Contribution(b.I, 1e-3))

    def test_only_v_contributed_or_probed_named_branches_get_unknowns(self):
        """Naming does not by itself create an unknown.

        A branch earns a current unknown by being V-contributed or
        current-probed; an I-contributed branch is just a stamp on two
        KCL rows and needs none.  Naming had to leave that rule alone --
        otherwise every named branch would cost a row.
        """
        class Probe(Behavioural):
            instparams = []

            @staticmethod
            def analog(p, m):
                mid = Node('mid')
                a = Branch(p, mid, 'a')       # V-contributed -> unknown
                b = Branch(mid, m, 'b')       # I-contributed -> none
                return (Contribution(a.V, 1.0),
                        Contribution(b.I, a.I * 2.0))

        e = _mk(Probe, 'p', 'm')
        assert {br.name for br in e.branches} == {'a'}

    def test_a_probe_reads_the_named_branch_it_names(self):
        """Two named branches across ONE pair, only one of them probed.

        Merged under a node-pair key, `I(a)` would read whichever branch
        won the merge and the current-controlled contribution would be
        wrong rather than absent -- the expensive kind of failure.
        """
        class Pick(Behavioural):
            instparams = [Parameter(name='r1', desc='R1', unit='ohm',
                                    default=1e3),
                          Parameter(name='r2', desc='R2', unit='ohm',
                                    default=4e3),
                          Parameter(name='k', desc='gain', unit='',
                                    default=3.0)]

            @staticmethod
            def analog(p, m, out):
                a = Branch(p, m, 'a')
                b = Branch(p, m, 'b')
                return (Contribution(a.V, a.I * r1),               # noqa
                        Contribution(b.V, b.I * r2),               # noqa
                        Contribution(Branch(out, m).I, k * a.I))   # noqa

        ## Named differently from the model's parameters: a local `r1`
        ## here would shadow the parameter symbol the DSL injects into
        ## `analog`'s scope, and the class would not compile.
        R1, R2, K, VIN, RL = 1e3, 4e3, 3.0, 2.0, 1e3
        c = SubCircuit()
        n_in, n_out = c.add_node('in'), c.add_node('out')
        c['vs'] = VS(n_in, gnd, v=VIN)
        c['E'] = Pick(n_in, gnd, n_out, r1=R1, r2=R2, k=K)
        c['rl'] = R(n_out, gnd, r=RL)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()

        ## The probe must see branch a's current (VIN/R1), not b's
        ## (VIN/R2) and not their sum.
        want = -K * (VIN / R1) * RL
        assert float(res.v(n_out, gnd)) == pytest.approx(want, rel=1e-9)


class TestSimparam(object):

    def test_known_names_resolve_to_what_pycircuit_actually_does(self):
        ## gmin is a continuation schedule in the DC solver, not a
        ## standing shunt, so the honest answer is zero.
        ## float(), not ==: `sympy.Float(0.0) == 0` is False, because
        ## sympy's == is STRUCTURAL equality, not numeric.
        assert float(simparam('gmin', 1e-12)) == 0.0
        assert float(simparam('scale', 0.5)) == 1.0
        assert float(simparam('shrink')) == 1.0

    def test_an_unknown_name_takes_the_caller_default(self):
        assert float(simparam('imelt', 42.0)) == 42.0
        assert float(simparam('iteration', 7)) == 7.0

    def test_an_unknown_name_with_no_default_is_an_error(self):
        """Not a silent zero: the model would be reading nothing."""
        with pytest.raises(ValueError, match='no default'):
            simparam('completely_made_up')

    def test_it_reaches_the_stamps(self):
        """Compile-time resolution, so it lands in the conductance."""
        class Shunt(Behavioural):
            instparams = [Parameter(name='g0', desc='G', unit='S',
                                    default=1e-3)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                gm = simparam('gmin', 1e-9)
                return Contribution(b.I, (g0 + gm) * b.V)          # noqa

        e = _mk(Shunt, 'p', 'm', g0=1e-3)
        G = np.asarray(e.G(np.array([1.0, 0.0])), float)
        ## gmin resolved to 0, so the stamp is exactly g0
        assert G[0, 0] == pytest.approx(1e-3, rel=1e-15)

    def test_it_composes_with_the_let_chain(self):
        class ShuntChained(Behavioural):
            instparams = [Parameter(name='g0', desc='G', unit='S',
                                    default=1e-3)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                g = var((g0 + simparam('gmin', 1e-9)) * b.V, 'g')  # noqa
                return Contribution(b.I, g)

        e = _mk(ShuntChained, 'p', 'm', g0=1e-3)
        assert ShuntChained._hdl_info['chained'] is True
        G = np.asarray(e.G(np.array([1.0, 0.0])), float)
        assert G[0, 0] == pytest.approx(1e-3, rel=1e-15)


class TestMfactorIsNotNeeded(object):
    """Why `$mfactor` is absent, recorded so it is a decision.

    Reading the shipped PSP 103.8.2 source rather than the LRM: PSP
    references `$mfactor` at exactly two lines
    (`PSP103_module.include:2183-2184`), both inside the noise
    OPERATING-POINT output block, and its own Changelog line 3 reads
    "Remove $mfactor from OP output variables calculation".  Device
    multiplicity in PSP flows through its own `MULT_i` parameter, which
    is an ordinary model parameter the DSL already handles.

    So there is nothing to build here for the PSP target.  Adding an
    auto-scaling multiplicity would change the meaning of every existing
    Behavioural element to serve a construct the model does not use.
    """

    def test_multiplicity_is_expressible_as_an_ordinary_parameter(self):
        class Mult(Behavioural):
            instparams = [Parameter(name='g0', desc='G', unit='S',
                                    default=1e-3),
                          Parameter(name='mult', desc='devices in parallel',
                                    unit='', default=1.0)]

            @staticmethod
            def analog(p, m):
                b = Branch(p, m)
                return Contribution(b.I, mult * g0 * b.V)          # noqa

        for n in (1.0, 4.0, 17.0):
            e = _mk(Mult, 'p', 'm', g0=1e-3, mult=n)
            G = np.asarray(e.G(np.array([1.0, 0.0])), float)
            assert G[0, 0] == pytest.approx(n * 1e-3, rel=1e-12)
