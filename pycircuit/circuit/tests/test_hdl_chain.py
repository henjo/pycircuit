"""The let-chain compilation path (`hdl.var`).

A compact model of any realism names its intermediate quantities and then
reuses them.  Written as one expression, the reuse turns a linear DAG into
an exponential TREE -- differentiate that with ``Matrix.jacobian`` and the
compiler never returns.  ``var()`` keeps each intermediate as a symbol
with its own defining equation, and the compiler emits a topological
let-chain plus a forward-accumulated Jacobian.

These tests pin three things: the chain agrees with the eager path where
both can run, its Jacobian agrees with finite differences where only it
can run, and the classification of terms into i/q/u still holds when the
term is hidden behind an intermediate.
"""
import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import SubCircuit, VS, R
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, TIME,
                                   ddt, var, vt)
from pycircuit.utilities.param import Parameter


def _fd_jacobian(elem, x, h=1e-6):
    """Central-difference Jacobian of the compiled residual."""
    J = np.zeros((len(x), len(x)))
    for j in range(len(x)):
        step = h * max(1.0, abs(float(x[j])))
        xp, xm = np.array(x, float), np.array(x, float)
        xp[j] += step
        xm[j] -= step
        J[:, j] = (np.asarray(elem.i(xp), float)
                   - np.asarray(elem.i(xm), float)) / (2 * step)
    return J


class _DiodeEager(Behavioural):
    instparams = [Parameter(name='IS', desc='Is', unit='A', default=1e-14),
                  Parameter(name='N', desc='ideality', unit='', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        u = b.V / (N * vt())                                   # noqa: F821
        return Contribution(b.I, IS * (sympy.exp(u) - 1))      # noqa: F821


class _DiodeChained(Behavioural):
    """The same device, written through an intermediate."""
    instparams = [Parameter(name='IS', desc='Is', unit='A', default=1e-14),
                  Parameter(name='N', desc='ideality', unit='', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        u = var(b.V / (N * vt()), 'u')                         # noqa: F821
        return Contribution(b.I, IS * (sympy.exp(u) - 1))      # noqa: F821


def _mk(cls, **kw):
    e = cls(cm.Node('a'), cm.Node('b'), **kw)
    e.update_iparv()
    return e


class TestAgreesWithEagerPath(object):
    """Where both paths can run, they must produce the same numbers."""

    @pytest.mark.parametrize('v', [-0.5, -0.02, 0.0, 0.3, 0.6, 0.75])
    def test_i_and_G_match(self, v):
        eager = _mk(_DiodeEager, IS=2e-14, N=1.3)
        chained = _mk(_DiodeChained, IS=2e-14, N=1.3)
        x = np.array([v, 0.0])
        assert np.allclose(np.asarray(eager.i(x), float),
                           np.asarray(chained.i(x), float),
                           rtol=1e-12, atol=1e-300)
        assert np.allclose(np.asarray(eager.G(x), float),
                           np.asarray(chained.G(x), float),
                           rtol=1e-12, atol=1e-300)

    def test_the_two_paths_are_actually_different_paths(self):
        """Guard against the test passing because nothing was chained."""
        assert _DiodeChained._hdl_info['chained'] is True
        assert _DiodeEager._hdl_info['chained'] is False

    def test_dc_solution_matches(self):
        got = []
        for cls in (_DiodeEager, _DiodeChained):
            c = SubCircuit()
            na, nb = c.add_node('a'), c.add_node('b')
            c['vs'] = VS(na, gnd, v=1.0)
            c['R'] = R(na, nb, r=1e3)
            c['D'] = cls(nb, gnd, IS=2e-14, N=1.3)
            c.update_iparv()
            got.append(float(DC(c, toolkit=numeric).solve().v(nb, gnd)))
        assert abs(got[0] - got[1]) < 1e-12


class TestForwardAccumulation(object):
    """The Jacobian, checked against finite differences."""

    def test_deeply_nested_reuse(self):
        """Each level consumes the previous one twice.

        The eager path builds 2**depth copies of the innermost expression
        and does not finish; the chain builds `depth` definitions.
        """
        def analog(d, s):
            b = Branch(d, s)
            u = b.V / vt()
            for k in range(24):
                u = var(sympy.sqrt(1 + u * u) + u * sympy.tanh(u) / (2 + k),
                        'u%d' % k)
            return Contribution(b.I, IS * u)                    # noqa: F821

        cls = type('Deep', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-6)],
            'analog': staticmethod(analog)})
        e = _mk(cls, IS=1e-6)
        x = np.array([0.7, 0.0])
        J = np.asarray(e.G(x), float)
        fd = _fd_jacobian(e, x)
        assert np.allclose(J, fd, rtol=1e-6, atol=1e-9 * np.abs(J).max())

    def test_multi_terminal_jacobian(self):
        """Four terminals, three contributions, shared intermediates."""
        def analog(d, g, s, b):
            bds, bgs, bbs = Branch(d, s), Branch(g, s), Branch(b, s)
            vg = var(bgs.V / vt(), 'vg')
            vd = var(bds.V / vt(), 'vd')
            vb = var(bbs.V / vt(), 'vb')
            surf = var(sympy.sqrt(1e-9 + vg * vg + vb * vb), 'surf')
            chan = var(surf * sympy.tanh(vd / (1 + surf)), 'chan')
            return (Contribution(bds.I, IS * chan),             # noqa: F821
                    Contribution(bgs.I, IS * 1e-3 * surf),      # noqa: F821
                    Contribution(bbs.I, IS * 1e-4 * vb * surf))  # noqa: F821

        cls = type('Quad', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-6)],
            'analog': staticmethod(analog)})
        e = cls(cm.Node('d'), cm.Node('g'), cm.Node('s'), cm.Node('b'),
                IS=1e-6)
        e.update_iparv()
        for x in (np.array([0.9, 1.1, 0.0, -0.3]),
                  np.array([0.05, 0.4, 0.0, 0.0]),
                  np.array([1.8, 1.8, 0.2, -1.0])):
            J = np.asarray(e.G(x), float)
            fd = _fd_jacobian(e, x)
            assert np.allclose(J, fd, rtol=1e-5,
                               atol=1e-8 * max(1.0, np.abs(J).max())), x

    def test_an_intermediate_no_x_depends_on_gives_zero_columns(self):
        """A constant intermediate must not manufacture conductance."""
        def analog(p, m):
            b = Branch(p, m)
            k = var(IS * 3.0, 'k')                              # noqa: F821
            return Contribution(b.I, k + b.V * 1e-3)
        cls = type('ConstInter', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-6)],
            'analog': staticmethod(analog)})
        e = _mk(cls, IS=1e-6)
        G = np.asarray(e.G(np.array([0.5, 0.0])), float)
        assert np.allclose(G, 1e-3 * np.array([[1., -1.], [-1., 1.]]))


class TestTermClassification(object):
    """`i` / `q` / `u` routing must survive an opaque intermediate.

    `_split_terms` decides "is this a source?" by asking whether a term
    mentions a solution symbol.  An intermediate is just a symbol, so
    unless the compiler tracks which intermediates depend on the
    solution, a term written through one is silently routed into `u(t)`
    and the device loses its conductance entirely.
    """

    def test_conductance_is_not_misfiled_as_a_source(self):
        def analog(p, m):
            b = Branch(p, m)
            g = var(b.V * 1e-3, 'g')
            return Contribution(b.I, g)
        cls = type('ChainedR', (Behavioural,), {'instparams': [],
                                                'analog': staticmethod(analog)})
        e = _mk(cls)
        x = np.array([1.0, 0.0])
        assert np.allclose(np.asarray(e.i(x), float), [1e-3, -1e-3])
        assert np.allclose(np.asarray(e.u(0.0), float), [0.0, 0.0])
        assert np.allclose(np.asarray(e.G(x), float),
                           1e-3 * np.array([[1., -1.], [-1., 1.]]))

    def test_charge_behind_an_intermediate_reaches_C(self):
        def analog(p, m):
            b = Branch(p, m)
            q = var(1e-12 * b.V * (1 + 0.1 * b.V), 'q')
            return Contribution(b.I, ddt(q))
        cls = type('ChainedC', (Behavioural,), {'instparams': [],
                                                'analog': staticmethod(analog)})
        e = _mk(cls)
        x = np.array([2.0, 0.0])
        ## dq/dv = 1e-12 * (1 + 0.2*v)
        want = 1e-12 * 1.4
        C = np.asarray(e.C(x), float)
        assert np.allclose(C, want * np.array([[1., -1.], [-1., 1.]]),
                           rtol=1e-9)
        assert np.allclose(np.asarray(e.i(x), float), 0.0)

    def test_a_bias_constant_behind_an_intermediate_lands_where_it_does_eagerly(self):
        """`u` is reserved for time-dependent drive.

        A contribution that is merely constant in the solution is still
        part of the static residual, so both paths put it in `i` with a
        zero conductance.  Pinned on both so the chain cannot drift.
        """
        def eager(p, m):
            return Contribution(Branch(p, m).I, IS * 2.0)        # noqa: F821

        def chained(p, m):
            return Contribution(Branch(p, m).I, var(IS * 2.0, 's'))  # noqa

        out = []
        for name, fn in (('EagerI', eager), ('ChainedI', chained)):
            cls = type(name, (Behavioural,), {
                'instparams': [Parameter(name='IS', desc='I', unit='A',
                                         default=1e-6)],
                'analog': staticmethod(fn)})
            e = _mk(cls, IS=1e-6)
            out.append((np.asarray(e.i(np.array([0.5, 0.0])), float),
                        np.asarray(e.u(0.0), float),
                        np.asarray(e.G(np.array([0.5, 0.0])), float)))
        for i_, u_, G_ in out:
            assert np.allclose(i_, [2e-6, -2e-6])
            assert np.allclose(u_, 0.0)
            assert np.allclose(G_, 0.0)

    def test_a_time_dependent_source_behind_an_intermediate_reaches_u(self):
        """This is the case the chain could silently drop.

        `u` is compiled without the solution vector, so it sees only the
        x-free part of the chain.  If that sub-chain were not emitted the
        source would evaluate to nothing and the device would be dead.
        """
        def analog(p, m):
            amp = var(IS * 3.0, 'amp')                          # noqa: F821
            drive = var(amp * sympy.sin(2e6 * TIME), 'drive')
            return Contribution(Branch(p, m).I, drive)
        cls = type('ChainedSrc', (Behavioural,), {
            'instparams': [Parameter(name='IS', desc='I', unit='A',
                                     default=1e-6)],
            'analog': staticmethod(analog)})
        e = _mk(cls, IS=1e-6)
        for t in (0.0, 3e-7, 1.1e-6):
            want = 3e-6 * np.sin(2e6 * t)
            assert np.allclose(np.asarray(e.u(t), float), [want, -want],
                               atol=1e-18)
        assert np.allclose(np.asarray(e.G(np.array([0.5, 0.0])), float), 0.0)


class TestCapabilityRefusals(object):
    """What the chain path gives up, it must give up loudly."""

    def test_pcnr_is_declined_not_silently_missed(self):
        """A diode behind a `var` has no visible exponential.

        The PCNR shape detector reads the exponential straight out of the
        expression.  Behind an intermediate it would find none and would
        register the device as "nothing to limit" -- which looks
        identical to a linear device.  The compiler must decline instead.
        """
        assert _DiodeChained._hdl_info['pcnr_funcs'] is None
        assert _DiodeEager._hdl_info['pcnr_funcs'] is not None

    def test_no_constant_stamp_claim(self):
        assert _DiodeChained._hdl_info['const_G'] is False
        assert _DiodeChained._hdl_info['const_C'] is False

    def test_pure_and_symbolic_specs_are_absent(self):
        assert _DiodeChained._hdl_info['pure_spec'] is None
        assert _DiodeChained._hdl_info['sym_spec'] is None


class TestUsableInAnalyses(object):

    def test_transient_runs_and_matches_the_eager_device(self):
        got = []
        for cls in (_DiodeEager, _DiodeChained):
            c = SubCircuit()
            na, nb = c.add_node('a'), c.add_node('b')
            c['vs'] = VS(na, gnd, v=0.9)
            c['R'] = R(na, nb, r=1e3)
            c['D'] = cls(nb, gnd, IS=2e-14, N=1.3)
            c.update_iparv()
            res = Transient(c, toolkit=numeric).solve(tend=1e-7,
                                                      timestep=1e-9)
            got.append(np.asarray(res.v(nb, gnd).y, float)[-1])
        assert abs(got[0] - got[1]) < 1e-9


class TestVarMisuse(object):

    def test_var_outside_a_model_is_an_error(self):
        with pytest.raises(RuntimeError):
            var(sympy.Symbol('z'), 'z')
