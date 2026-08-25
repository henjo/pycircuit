"""The chained code path as a first-class citizen.

Every production model in this tree calls ``var()`` and therefore
compiles through the let-chain.  Downstream machinery, though, keeps
being written against the EAGER path: it reads a fact off an assembled
expression, finds an opaque intermediate symbol there instead, and
either refuses or silently does nothing.  Three instances of that were
found and fixed in one week (``@cross`` and ``$limit`` install skipped by
an early return; ``$limit`` markers not stripped from ``var()``
definitions).  This file pins two more, plus the two they uncovered:

1. **State operators inside ``var()``.**  ``generate_code`` collected
   ``idt`` / ``idtmod`` / ``laplace_nd`` / ``laplace_zp`` applications
   from the STATEMENTS only, so ``var(idt(b.V))`` never had a state
   allocated and died in the printer -- unless the same application also
   appeared in a statement, in which case ``resolve()``'s
   ``state_subst`` rescued it and the model compiled.  Whether a model
   worked depended on where ELSE its author had written the operator,
   which is why `test_the_case_that_used_to_work_by_accident` is here.

2. **PCNR admission (roadmap 10.2).**  The shape detector recovered a
   junction's ``IS`` and ``VT`` by reading ``exp(arg)`` out of the
   assembled expression, and a let-chain hides that behind a symbol, so
   the gate refused every chained model outright.  It now walks the
   chain instead: prune to what one contribution reaches, substitute the
   branch voltage into each definition separately, and
   forward-accumulate ``d/dv`` for the exponentials.  Both steps are
   linear in the number of definitions -- **the chain is never
   flattened**, which is the whole point of ``var()``.

3. ``_ChainPrinter`` had no ``_wrapfloor`` entry.  It could not show
   before (1), because no chained model could carry an ``idtmod``.

4. ``$param_given`` inside a ``var()`` definition became no argument at
   all -- ``given_syms`` was read off the assembled vectors too -- so the
   generated function raised ``NameError`` at call time from a model that
   compiled clean.  Found while writing (2), which needs to know what a
   contribution reads.

The discriminating table for (2), from roadmap sections 9 and 10.1:

    flat, no charge         chained=False   pcnr_i=True
    flat, WITH charge       chained=False   pcnr_i=True
    chained, with charge    chained=True    pcnr_i=True   <- was False

`TestTheThreeElementTable` is that table, and it is the regression test
for the whole change.
"""
import numpy as np
import pytest
import sympy

from pycircuit.circuit import gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import IS as ISrc, R, VS
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import pcnr as pcnr_mod
from pycircuit.circuit.hdl import (Behavioural, BehaviouralMeta, Branch,
                                   Contribution, Node, ddt, explain, idt,
                                   idtmod, laplace_nd, laplace_zp, limexp,
                                   param_given, var, vt)
from pycircuit.utilities.param import Parameter


def _make(name, analog, params, terminals=None):
    """Compile a class from a plain function inside a test."""
    body = dict(instparams=[Parameter(name=n, desc=n, unit='', default=d)
                            for n, d in params],
                analog=staticmethod(analog))
    if terminals is not None:
        body['terminals'] = terminals
    return BehaviouralMeta(name, (Behavioural,), body)


def _nx(cls):
    """Length of the unknown vector -- `explain`'s layout, as a number."""
    info = cls._hdl_info
    return (len(info['terminalnames']) + len(info['internalnames'])
            + len(info['state_meta']['statenames'])
            + len(info['branchpairs']))


# ----------------------------------------------------------------------
# 1.  State operators inside var()
# ----------------------------------------------------------------------

_TWOPORT = ['ip', 'inn', 'op', 'on']
_SPAR = (('g', 2.0), ('ic', 0.3))


def _flat_idt(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    return Contribution(bo.V, g * idt(bi.V, ic))                 # noqa: F821


def _chained_idt(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    z = var(idt(bi.V, ic), 'z')                                  # noqa: F821
    return Contribution(bo.V, g * z)                             # noqa: F821


def _flat_idtmod(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    return Contribution(bo.V, g * idtmod(bi.V, ic, 2.0, -1.0))   # noqa: F821


def _chained_idtmod(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    z = var(idtmod(bi.V, ic, 2.0, -1.0), 'z')                    # noqa: F821
    return Contribution(bo.V, g * z)                             # noqa: F821


_NUM, _DEN = [1.0], [1.0, 1e-3, 1e-7]


def _flat_lap_nd(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    return Contribution(bo.V, g * laplace_nd(bi.V, _NUM, _DEN))  # noqa: F821


def _chained_lap_nd(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    z = var(laplace_nd(bi.V, _NUM, _DEN), 'z')
    return Contribution(bo.V, g * z)                             # noqa: F821


_POLES = [-1e3, 0.0, -2e3, 0.0]


def _flat_lap_zp(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    return Contribution(bo.V, g * laplace_zp(bi.V, [], _POLES))  # noqa: F821


def _chained_lap_zp(ip, inn, op, on):
    bi, bo = Branch(ip, inn), Branch(op, on)
    z = var(laplace_zp(bi.V, [], _POLES), 'z')
    return Contribution(bo.V, g * z)                             # noqa: F821


_STATE_CASES = [
    ('idt', _flat_idt, _chained_idt, 1),
    ('idtmod', _flat_idtmod, _chained_idtmod, 1),
    ('laplace_nd', _flat_lap_nd, _chained_lap_nd, 2),
    ('laplace_zp', _flat_lap_zp, _chained_lap_zp, 2),
]


class TestStateOperatorsInsideVar(object):
    """`var(idt(...))` must compile, and to the same numbers as flat.

    Compiling is not enough on its own: a compiler that dropped the state
    silently would also "compile".  Every case therefore asserts the
    STATE COUNT as well, which is the thing a silent drop changes.
    """

    @pytest.mark.parametrize('name, flat, chained, nstates', _STATE_CASES)
    def test_compiles_with_the_same_states_as_the_flat_model(
            self, name, flat, chained, nstates):
        A = _make('Flat_' + name, flat, _SPAR, _TWOPORT)
        B = _make('Chained_' + name, chained, _SPAR, _TWOPORT)
        assert A._hdl_info['chained'] is False
        assert B._hdl_info['chained'] is True
        ## The count, not merely "some states": before the fix a
        ## `var()`-only application allocated NONE.
        assert len(A._hdl_info['state_meta']['statenames']) == nstates
        assert len(B._hdl_info['state_meta']['statenames']) == nstates
        ## And the same through the public diagnostic, so a reader of
        ## explain() is told the truth about a chained model too.
        el = B('a', 'b', 'c', 'd')
        el.update_iparv()
        assert '\n  %d state' % nstates in \
            explain(el, source=False, symbolic=False)

    @pytest.mark.parametrize('name, flat, chained, nstates', _STATE_CASES)
    def test_the_numbers_match_the_flat_model(self, name, flat, chained,
                                              nstates):
        """Same i/q/G/C at random points.  Both models have the same
        unknown layout precisely BECAUSE the state count matches, so a
        dropped state would fail here as a shape error rather than a
        wrong number -- which is why the count is asserted above too."""
        A = _make('FlatN_' + name, flat, _SPAR, _TWOPORT)
        B = _make('ChainedN_' + name, chained, _SPAR, _TWOPORT)
        a, b = A('a', 'b', 'c', 'd'), B('a', 'b', 'c', 'd')
        a.update_iparv()
        b.update_iparv()
        assert _nx(A) == _nx(B)
        rng = np.random.RandomState(20260825)
        for _ in range(5):
            x = rng.randn(_nx(A)) * 2.0
            for which in ('i', 'q', 'G', 'C'):
                fa = np.asarray(getattr(a, which)(x), float)
                fb = np.asarray(getattr(b, which)(x), float)
                assert np.allclose(fa, fb, rtol=1e-12, atol=1e-15), \
                    '%s.%s disagrees' % (name, which)

    def test_the_case_that_used_to_work_by_accident(self):
        """The same `idt` application inside a `var` AND in a statement.

        `resolve()` applies `state_subst` to intermediates, so this
        combination compiled before the fix -- one state, correct
        numbers -- while the `var`-only spelling raised
        `PrintMethodNotImplementedError`.  The failure therefore depended
        on where else the author had written the operator.  Both
        spellings now work, and this one must still allocate exactly ONE
        state: the application is one application however many places
        mention it.
        """
        def analog(ip, inn, op, on):
            bi, bo = Branch(ip, inn), Branch(op, on)
            z = var(idt(bi.V, ic), 'z')                          # noqa: F821
            ## `idt(bi.V, ic)` again, in the statement this time.
            return Contribution(bo.V, g * z + idt(bi.V, ic))     # noqa: F821

        B = _make('Both', analog, _SPAR, _TWOPORT)
        assert B._hdl_info['chained'] is True
        assert len(B._hdl_info['state_meta']['statenames']) == 1

        def flat(ip, inn, op, on):
            bi, bo = Branch(ip, inn), Branch(op, on)
            return Contribution(bo.V, (g + 1) * idt(bi.V, ic))   # noqa: F821

        A = _make('BothFlat', flat, _SPAR, _TWOPORT)
        a, b = A('a', 'b', 'c', 'd'), B('a', 'b', 'c', 'd')
        a.update_iparv()
        b.update_iparv()
        rng = np.random.RandomState(11)
        for _ in range(4):
            x = rng.randn(_nx(A))
            assert np.allclose(np.asarray(a.i(x), float),
                               np.asarray(b.i(x), float))
            assert np.allclose(np.asarray(a.G(x), float),
                               np.asarray(b.G(x), float))

    def test_a_chained_integrator_integrates_in_a_real_transient(self):
        """End to end: the state has to be INTEGRATED by the simulator,
        not merely allocated.  A constant 1 V into `g * idt(V)` is a
        ramp of slope `g`, so a state that was allocated and then never
        advanced shows up as a flat line -- which the closed-form
        assertion below catches and an equality against the flat model
        would not."""
        import pycircuit.circuit.circuit as _cm
        from pycircuit.circuit.integrator import EulerIntegrator
        A = _make('FlatTran', _flat_idt, _SPAR, _TWOPORT)
        B = _make('ChainedTran', _chained_idt, _SPAR, _TWOPORT)
        out = []
        for cls in (A, B):
            _cm.default_toolkit = numeric
            c = SubCircuit()
            n1, n2 = c.add_node('n1'), c.add_node('n2')
            c['vs'] = VS(n1, gnd, v=1.0)
            c['rl'] = R(n2, gnd, r=1e3)
            c['x'] = cls(n1, gnd, n2, gnd, g=2.0, ic=0.3)
            c.update_iparv()
            wf = Transient(c, toolkit=numeric, uic=True,
                           integrator=EulerIntegrator()).solve(
                               tend=0.5, timestep=1e-2).v(n2, gnd)
            out.append((np.asarray(wf.x, float).ravel(),
                        np.asarray(wf.y, float).ravel()))
        t, y = out[1]
        assert np.allclose(out[0][1], y, rtol=1e-12, atol=1e-14)
        ## `g * (ic + t)`, the closed form -- `uic` seeds the state from
        ## the instance parameter literally named `ic`.  A state that was
        ## allocated but never advanced is a flat line at `g * ic`, which
        ## this catches and an equality against the flat model does not.
        assert t[-1] == pytest.approx(0.5)
        assert y[-1] == pytest.approx(2.0 * (0.3 + 0.5), rel=1e-9)
        assert np.all(np.diff(y[1:]) > 0)


# ----------------------------------------------------------------------
# 2.  PCNR admission for chained models  (roadmap 10.2)
# ----------------------------------------------------------------------

_DPAR = (('isat', 1e-14), ('n', 1.0), ('cj0', 1e-12))


def _diode_flat_noq(plus, minus):
    b = Branch(plus, minus)
    return Contribution(b.I, isat * (limexp(b.V / (n * vt())) - 1))


def _diode_flat_q(plus, minus):
    b = Branch(plus, minus)
    return Contribution(b.I, isat * (limexp(b.V / (n * vt())) - 1)
                        + ddt(cj0 * b.V))                        # noqa: F821


def _diode_chained_q(plus, minus):
    """The same device through a let-chain, with the exponential TWO
    definitions deep so the detector has to walk, not peek.

    `limexp` rather than `exp`, which is what its docstring tells a
    PCNR-eligible junction to use -- and now that a chained model can be
    eligible, a chained model can take that advice.  It also keeps the
    Newton iterates that PCNR does NOT limit (the MNA vector) out of
    overflow, which is the whole reason the primitive exists.
    """
    b = Branch(plus, minus)
    vd = var(b.V, 'vd')
    arg = var(vd / (n * vt()), 'arg')                            # noqa: F821
    i = var(isat * (limexp(arg) - 1), 'i')                       # noqa: F821
    return Contribution(b.I, i + ddt(cj0 * vd))                  # noqa: F821


_FlatNoQ = _make('DiodeFlatNoQ', _diode_flat_noq, _DPAR)
_FlatQ = _make('DiodeFlatQ', _diode_flat_q, _DPAR)
_ChainedQ = _make('DiodeChainedQ', _diode_chained_q, _DPAR)


def _params(el):
    return {k: getattr(el.iparv, k) for k in el._hdl_paramnames}


def _inst(cls, **kw):
    el = cls('a', 'b', **kw)
    el.update_iparv()
    return el


class TestTheThreeElementTable(object):
    """Roadmap 10.1's table, which is what named `var()` as the blocker.

    Three elements differing by one thing each; only the last line
    changes, and it is the line the whole change is about.
    """

    def test_flat_no_charge(self):
        el = _inst(_FlatNoQ)
        assert _FlatNoQ._hdl_info['chained'] is False
        assert hasattr(el, 'pcnr_i')

    def test_flat_with_charge(self):
        el = _inst(_FlatQ)
        assert _FlatQ._hdl_info['chained'] is False
        assert hasattr(el, 'pcnr_i')

    def test_chained_with_charge(self):
        el = _inst(_ChainedQ)
        assert _ChainedQ._hdl_info['chained'] is True
        assert hasattr(el, 'pcnr_i')
        assert el.pcnr_junctions == ((0, 1),)
        assert 'PCNR: 1 junction' in explain(el, source=False,
                                             symbolic=False)


class TestTheRecoveredScales(object):

    def test_vt_and_is_match_the_flat_model(self):
        """`pcnr_scales` is what `pnjlim` is driven by, and it is
        recovered from the shape rather than from a parameter, so the
        chained and flat spellings of one device must agree.

        Not bit-identical, and the reason is worth recording: the flat
        detector simplifies ``1/(1/(n*vt))`` back to ``n*vt``
        symbolically, while the chained one evaluates the reciprocal of
        a forward-accumulated derivative at run time.  That is one extra
        round trip through a float division -- a few ulp, not a
        different number.
        """
        a, b = _inst(_FlatQ), _inst(_ChainedQ)
        for kw in ({}, dict(isat=3e-15, n=1.7)):
            a, b = _inst(_FlatQ, **kw), _inst(_ChainedQ, **kw)
            va, ia = a.pcnr_scales(_params(a), _epar())
            vb, ib = b.pcnr_scales(_params(b), _epar())
            assert vb == pytest.approx(va, rel=1e-14)
            assert ib == pytest.approx(ia, rel=1e-14)
            ## And they are the physical numbers, not merely equal.
            assert ia == pytest.approx(kw.get('isat', 1e-14), rel=1e-9)


def _epar():
    from pycircuit.circuit import defaultepar
    return defaultepar


class TestTheCompiledJunctionFunctions(object):

    _SWEEP = np.concatenate([np.linspace(-1.0, 0.0, 11),
                             np.linspace(0.0, 0.85, 35)])

    def test_i_and_didv_agree_with_the_flat_equivalent(self):
        a, b = _inst(_FlatQ), _inst(_ChainedQ)
        pa, pb = _params(a), _params(b)
        for v in self._SWEEP:
            ia = np.asarray(a.pcnr_i(v, pa, _epar(), a.toolkit), float)
            ib = np.asarray(b.pcnr_i(v, pb, _epar(), b.toolkit), float)
            ga = np.asarray(a.pcnr_didv(v, pa, _epar(), a.toolkit), float)
            gb = np.asarray(b.pcnr_didv(v, pb, _epar(), b.toolkit), float)
            assert np.allclose(ia, ib, rtol=1e-12, atol=1e-20), v
            assert np.allclose(ga, gb, rtol=1e-12, atol=1e-20), v
            ## The stamp is a pair (+i, -i) -- the shape `pcnr.py`
            ## re-stamps with.
            assert ib.shape == (2,) and ib[0] == -ib[1]

    def test_didv_is_the_derivative_of_i(self):
        """Finite differences against the CHAINED functions themselves,
        so this passes only if the forward accumulation is right -- not
        merely if it matches whatever the flat path did."""
        b = _inst(_ChainedQ)
        pb = _params(b)
        for v in np.linspace(-0.4, 0.8, 25):
            h = 1e-6 * max(1.0, abs(v))
            f = lambda vv: float(np.asarray(          # noqa: E731
                b.pcnr_i(vv, pb, _epar(), b.toolkit), float)[0])
            fd = (f(v + h) - f(v - h)) / (2 * h)
            an = float(np.asarray(b.pcnr_didv(v, pb, _epar(), b.toolkit),
                                  float)[0])
            assert an == pytest.approx(fd, rel=2e-5, abs=1e-18), v


class TestTheTracedTwin(object):
    """`pcnr_i` is also called under `jit` by the traced backend.

    The flat path serves that with a jax-printed `lambdify` of the same
    expression; a chained model has no such expression, so its twin is
    the SAME let-chain compiled with `numpy` bound to `jax.numpy` --
    `_ChainPrinter` emits `numpy.<f>` calls and nothing else.  Without
    this the capability would be declared and then fail the moment
    anything traced it, which is the failure mode this file is about.
    """

    def test_the_jax_twin_agrees_and_traces(self):
        jax = pytest.importorskip('jax')
        import jax.numpy as jnp

        class _TracedToolkit(object):
            jax = True

            @staticmethod
            def array(a):
                return jnp.array(a)

        el = _inst(_ChainedQ)
        p = _params(el)
        for v in (-0.3, 0.0, 0.5, 0.75):
            for which in ('pcnr_i', 'pcnr_didv'):
                ref = np.asarray(getattr(el, which)(v, p, _epar(),
                                                    el.toolkit), float)
                got = np.asarray(getattr(el, which)(v, p, _epar(),
                                                    _TracedToolkit), float)
                assert np.allclose(ref, got, rtol=1e-12, atol=1e-20), \
                    (which, v)
        ## And it must survive tracing, not merely evaluate eagerly.
        f = jax.jit(lambda vv: el.pcnr_i(vv, p, _epar(), _TracedToolkit))
        assert np.allclose(np.asarray(f(0.7), float),
                           np.asarray(el.pcnr_i(0.7, p, _epar(),
                                                el.toolkit), float),
                           rtol=1e-6)


class TestTheRefusalsThatSurvive(object):
    """Admitting chained models must not admit everything.

    Each of these is a chained model that is refused, and each names the
    refusal it is testing.  Declining beats claiming a capability
    falsely, and that principle is older than this change.
    """

    def test_a_state_still_refuses(self):
        def analog(plus, minus):
            b = Branch(plus, minus)
            i = var(isat * (sympy.exp(b.V / vt()) - 1), 'i')     # noqa: F821
            return Contribution(b.I, i + idt(b.V, 0.0))

        cls = _make('ChainedWithState', analog, (('isat', 1e-14),))
        assert cls._hdl_info['chained'] is True
        assert cls._hdl_info['state_meta']['statenames']
        assert not hasattr(_inst(cls), 'pcnr_i')

    def test_a_branch_current_unknown_still_refuses(self):
        def analog(plus, minus):
            b = Branch(plus, minus)
            e = var(isat * (sympy.exp(b.V / vt()) - 1), 'e')     # noqa: F821
            ## A V-contribution: PCNR's unknown is a branch VOLTAGE, and
            ## this device's constitutive relation is the other way up.
            return Contribution(b.V, e * 1.0)

        cls = _make('ChainedVContrib', analog, (('isat', 1e-14),))
        assert cls._hdl_info['chained'] is True
        assert cls._hdl_info['branchpairs']
        assert not hasattr(_inst(cls), 'pcnr_i')

    def test_a_current_spanning_two_branch_voltages_still_refuses(self):
        """The chained form of `f.free_symbols & xsyms`.  The chain is
        walked, so the second branch voltage is found inside a `var()`
        definition rather than in the statement.

        The exponential deliberately depends on BOTH branch voltages.  A
        controlled source whose exponent is a function of the OTHER
        branch alone is refused by a different check (its `exp` has zero
        derivative in `v`, so no scale is found), which would leave this
        test passing for the wrong reason -- checked by mutation.
        """
        def analog(a_, b_, c_, d_):
            bout, bctl = Branch(a_, b_), Branch(c_, d_)
            e = var(isat * (sympy.exp((bout.V + bctl.V) / vt()) - 1),
                    'e')                                         # noqa: F821
            return Contribution(bout.I, e)

        cls = _make('ChainedTwoBranch', analog, (('isat', 1e-14),),
                    terminals=['a_', 'b_', 'c_', 'd_'])
        assert cls._hdl_info['chained'] is True
        el = cls('a', 'b', 'c', 'd')
        el.update_iparv()
        assert not hasattr(el, 'pcnr_i')

        ## The control: the SAME model with the second branch dropped
        ## qualifies, so what is being measured is the second branch
        ## voltage and not the four terminals or the chain.
        def one_branch(a_, b_, c_, d_):
            bout = Branch(a_, b_)
            e = var(isat * (sympy.exp(bout.V / vt()) - 1), 'e')  # noqa: F821
            return Contribution(bout.I, e)

        ok = _make('ChainedOneBranch', one_branch, (('isat', 1e-14),),
                   terminals=['a_', 'b_', 'c_', 'd_'])
        el2 = ok('a', 'b', 'c', 'd')
        el2.update_iparv()
        assert hasattr(el2, 'pcnr_i')

    def test_a_series_resistance_refuses_and_removing_it_qualifies(self):
        """The `DiodeSpiceHdl` shape, and the roadmap's stated reason for
        it is wrong.

        10.2 says the series-`rs` diode is refused because "its current
        spans two branch voltages and the free-symbol test still refuses
        it".  Measured, each of the two contributions IS a function of
        its own branch voltage alone; what refuses the device is that
        the resistor's current is not exponential, and the rule is EVERY
        current, not some current.  The paired assertion localises it:
        delete the resistor contribution, keep everything else, and the
        same chained diode qualifies.
        """
        def with_rs(plus, minus):
            mid = Node('mid')
            bd, br = Branch(mid, minus), Branch(plus, mid)
            i = var(isat * (sympy.exp(bd.V / vt()) - 1), 'i')    # noqa: F821
            return (Contribution(bd.I, i),
                    Contribution(br.I, br.V / rs))               # noqa: F821

        def without_rs(plus, minus):
            bd = Branch(plus, minus)
            i = var(isat * (limexp(bd.V / vt()) - 1), 'i')       # noqa: F821
            return Contribution(bd.I, i)

        pars = (('isat', 1e-14), ('rs', 2.0))
        a = _inst(_make('ChainedSeriesRs', with_rs, pars))
        b = _inst(_make('ChainedNoRs', without_rs, pars))
        assert not hasattr(a, 'pcnr_i')
        assert hasattr(b, 'pcnr_i')

    def test_two_different_exponential_scales_refuse(self):
        """`len(scales) != 1`, not `scales` -- the limiter is driven by
        ONE `VT`, and a contribution carrying two of them has no single
        answer to give `pnjlim`.  A recombination current at `2*n*VT`
        alongside the ideal one is the realistic shape; both exponents
        are found by walking the chain, and disagreeing is what refuses
        the device."""
        def analog(plus, minus):
            b = Branch(plus, minus)
            vd = var(b.V, 'vd')
            i1 = var(isat * (sympy.exp(vd / vt()) - 1), 'i1')    # noqa: F821
            i2 = var(isat * 1e3 * (sympy.exp(vd / (2 * vt())) - 1), 'i2')
            return Contribution(b.I, i1 + i2)                    # noqa: F821

        cls = _make('ChainedTwoScales', analog, (('isat', 1e-14),))
        assert cls._hdl_info['chained'] is True
        assert not hasattr(_inst(cls), 'pcnr_i')

        ## The control: one scale, same shape, qualifies.
        def one_scale(plus, minus):
            b = Branch(plus, minus)
            vd = var(b.V, 'vd')
            i1 = var(isat * (sympy.exp(vd / vt()) - 1), 'i1')    # noqa: F821
            i2 = var(isat * 1e3 * (sympy.exp(vd / vt()) - 1), 'i2')
            return Contribution(b.I, i1 + i2)                    # noqa: F821

        assert hasattr(_inst(_make('ChainedOneScale', one_scale,
                                   (('isat', 1e-14),))), 'pcnr_i')

    def test_a_non_exponential_chained_current_refuses(self):
        def analog(plus, minus):
            b = Branch(plus, minus)
            u = var(b.V * b.V, 'u')
            return Contribution(b.I, gm * u)                     # noqa: F821

        cls = _make('ChainedSquare', analog, (('gm', 1e-3),))
        assert cls._hdl_info['chained'] is True
        assert not hasattr(_inst(cls), 'pcnr_i')

    def test_reading_a_param_given_flag_refuses(self):
        """`pcnr_i(v, params, epar, toolkit)` is handed the parameters
        and the temperature and nothing else.  A contribution that reads
        a `$param_given` flag cannot be evaluated from that, so it is
        declined rather than compiled into a NameError at solve time --
        which is what the flat path did before this change, quietly.

        The `given_names` assertion is a regression test for a THIRD
        instance of this file's root cause, found while writing the
        second: `given_syms` was collected from the assembled vectors
        only, so a `$param_given` used inside a `var()` definition became
        no argument at all and the compiled function raised
        `NameError: _hdl_given_isat` at call time -- from a model that
        compiled without a word.  The chain is scanned now, which is why
        the two spellings below agree.
        """
        ## Two spellings of one device.  Written out rather than shared
        ## through a helper: parameter names are bound into `analog`
        ## itself, and a nested function does not see them.
        def chained(plus, minus):
            b = Branch(plus, minus)
            sc = sympy.Piecewise((2.0, param_given('isat') > 0.5),
                                 (1.0, True))
            return Contribution(b.I, var(
                sc * isat * (limexp(b.V / vt()) - 1), 'i'))      # noqa: F821

        def flat_(plus, minus):
            b = Branch(plus, minus)
            sc = sympy.Piecewise((2.0, param_given('isat') > 0.5),
                                 (1.0, True))
            return Contribution(
                b.I, sc * isat * (limexp(b.V / vt()) - 1))       # noqa: F821

        cls = _make('ChainedParamGiven', chained, (('isat', 1e-14),))
        flat = _make('FlatParamGiven', flat_, (('isat', 1e-14),))
        assert cls._hdl_info['chained'] is True
        assert cls._hdl_info['given_names'] == ['isat']
        assert not hasattr(_inst(cls), 'pcnr_i')
        ## ... and neither does the flat spelling, for the same reason.
        assert not hasattr(_inst(flat), 'pcnr_i')
        ## The residual still evaluates, and to the flat model's number,
        ## with the flag both set and unset.
        x = np.array([0.6, 0.0])
        for kw in (dict(isat=2e-14), {}):
            a, b_ = _inst(flat, **kw), _inst(cls, **kw)
            assert np.allclose(np.asarray(a.i(x), float),
                               np.asarray(b_.i(x), float), rtol=1e-12)
        ## and `given` really does change the answer, so the comparison
        ## above is not comparing two copies of the same branch.
        assert not np.allclose(
            np.asarray(_inst(cls, isat=1e-14).i(x), float),
            np.asarray(_inst(cls).i(x), float))


class TestItIsNotInlined(object):
    """The chain must be walked, never flattened.

    Inlining is what `var()` exists to prevent: on a PSP-sized chain the
    flattened expression does not finish compiling.  A detector that
    inlined "just for the detection" would be quadratic-or-worse in the
    nesting depth, so a deep chain is the discriminating case -- it costs
    nothing to compile if the chain is walked and does not terminate in
    useful time if it is not.
    """

    def test_a_deep_chain_compiles_in_ordinary_time(self):
        import time

        def analog(plus, minus):
            b = Branch(plus, minus)
            u = var(b.V / vt(), 'u0')
            ## 60 definitions, each reading the previous one TWICE.
            ## Flattened, that is 2**60 leaves; walked, it is 60 lines.
            for k in range(60):
                u = var(u * 0.5 + u * 0.5 + 1e-9, 'u%d' % (k + 1))
            i = var(isat * (limexp(u) - 1), 'i')                 # noqa: F821
            return Contribution(b.I, i)

        t0 = time.time()
        cls = _make('DeepChain', analog, (('isat', 1e-14),))
        el = _inst(cls)
        dt = time.time() - t0
        assert cls._hdl_info['chained'] is True
        assert hasattr(el, 'pcnr_i'), 'the deep chain must still qualify'
        ## Generous by two orders of magnitude against the measured
        ## ~0.1 s; an inlining detector never gets here at all.
        assert dt < 30.0, 'compiling the deep chain took %.1f s' % dt
        vtv = el.pcnr_scales(_params(el), _epar())[0]
        assert vtv == pytest.approx(0.025852, rel=1e-3)


class TestARealPcnrSolve(object):
    """A capability that is declared but never exercised end to end is
    exactly the failure mode this whole package is about."""

    def test_solve_dc_converges_with_a_chained_participant(self):
        c = SubCircuit()
        n1 = c.add_node('n1')
        c['is'] = ISrc(gnd, n1, i=1e-3)
        c['d'] = _ChainedQ(n1, gnd, isat=1e-14, n=1.0)
        c.update_iparv()

        junctions = pcnr_mod.pcnr_junctions(c)
        assert len(junctions) == 1, 'the chained diode must opt in'

        x, v_lim, iters = pcnr_mod.solve_dc(c, gnd)
        vd = float(x[c.get_node_index(n1)])
        ## Against the closed form: 1 mA through IS = 1e-14, at the
        ## recovered VT (which is the number under test -- a wrong one
        ## would put this several tens of mV out).
        vtherm = el_vt = c['d'].pcnr_scales(_params(c['d']), _epar())[0]
        want = float(np.log(1e-3 / 1e-14 + 1.0)) * vtherm
        assert el_vt == pytest.approx(0.02585, rel=1e-2)
        assert vd == pytest.approx(want, rel=2e-3)
        ## PCNR's own unknown must have converged ONTO the branch
        ## voltage; that agreement is what the augmented system solves
        ## for, and it is the thing a wrong `pcnr_i` would break.
        assert float(v_lim[0]) == pytest.approx(vd, rel=1e-6)
        assert iters < 60

    def test_it_matches_the_flat_participant_and_plain_dc(self):
        """Three routes to one operating point: PCNR with the chained
        diode, PCNR with its flat twin, and the ordinary Newton solve."""
        got = []
        for cls in (_FlatQ, _ChainedQ):
            c = SubCircuit()
            n1 = c.add_node('n1')
            c['is'] = ISrc(gnd, n1, i=2e-3)
            c['d'] = cls(n1, gnd)
            c.update_iparv()
            x, _v, _it = pcnr_mod.solve_dc(c, gnd)
            got.append(float(x[c.get_node_index(n1)]))
        assert got[0] == pytest.approx(got[1], rel=1e-9)

        c = SubCircuit()
        n1 = c.add_node('n1')
        c['is'] = ISrc(gnd, n1, i=2e-3)
        c['d'] = _ChainedQ(n1, gnd)
        c.update_iparv()
        res = DC(c, toolkit=numeric).solve()
        assert float(res.v(n1, gnd)) == pytest.approx(got[1], rel=1e-6)

    def test_a_cold_start_that_needs_the_limiting(self):
        """Two chained junctions in series off a 5 V rail through a
        large resistor -- the clash over a shared linearisation point
        PCNR exists to remove.  It is also the case where a wrong
        recovered `IS`/`VT` shows up: `pnjlim` is driven by them."""
        c = SubCircuit()
        n1, n2 = c.add_node('n1'), c.add_node('n2')
        c['vs'] = VS(n1, gnd, v=5.0)
        c['r'] = R(n1, n2, r=1e3)
        c['d1'] = _ChainedQ(n2, gnd)
        c.update_iparv()
        x, v_lim, iters = pcnr_mod.solve_dc(c, gnd, maxiter=200)
        vd = float(x[c.get_node_index(n2)])
        assert 0.5 < vd < 1.0, vd
        ## The current through the resistor and through the diode agree.
        el = c['d1']
        i_d = float(np.asarray(el.pcnr_i(float(v_lim[0]), _params(el),
                                         _epar(), el.toolkit), float)[0])
        assert i_d == pytest.approx((5.0 - vd) / 1e3, rel=1e-4)
