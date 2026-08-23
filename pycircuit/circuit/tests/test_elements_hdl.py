# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""elements_hdl.py: every claim proven against the hand-written elements.

Two proof styles, per hdl.md:
* STAMP equivalence -- where the HDL formulation matches the hand-written
  one row for row, i/G/q/C agree at arbitrary x;
* WAVEFORM equivalence -- where only the physics matches (an HDL branch
  row may be the hand-written row times -1), a simulated circuit built
  from HDL elements reproduces the elements.py circuit's waveforms.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import (SubCircuit, VS, VSin, R, C, L, G,
                                        VCCS, Diode, Idt, Idtmod)
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl


def _n(el):
    return el.n


def _rand_x(n, seed):
    rng = np.random.RandomState(seed)
    return rng.uniform(-0.4, 0.4, n)


## ------------------------------------------------------------------------
## Stamp equivalence


@pytest.mark.parametrize('hdl_el,ref_el', [
    (lambda: eh.RHdl('p', 'n', r=2.2e3), lambda: R('p', 'n', r=2.2e3)),
    (lambda: eh.GHdl('p', 'n', g=3.3e-3), lambda: G('p', 'n', g=3.3e-3)),
    (lambda: eh.CHdl('p', 'n', c=4.7e-9), lambda: C('p', 'n', c=4.7e-9)),
    (lambda: eh.VCCSHdl('a', 'b', 'c', 'd', gm=2e-3),
     lambda: VCCS('a', 'b', 'c', 'd', gm=2e-3)),
])
def test_stamp_equivalence(hdl_el, ref_el):
    pycircuit.circuit.circuit.default_toolkit = numeric
    e1, e2 = hdl_el(), ref_el()
    e1.update_iparv(); e2.update_iparv()
    assert e1.n == e2.n
    x = _rand_x(e1.n, seed=7)
    for meth in ('i', 'q'):
        assert_allclose(np.asarray(getattr(e1, meth)(x), float),
                        np.asarray(getattr(e2, meth)(x), float), atol=1e-15)
    for meth in ('G', 'C'):
        assert_allclose(np.asarray(getattr(e1, meth)(x), float),
                        np.asarray(getattr(e2, meth)(x), float), atol=1e-15)


def test_diode_stamp_equivalence():
    """DiodeHdl equals the hand-written diode's UNLIMITED characteristic.

    A fresh reference element per operating point, deliberately:
    ``elements.Diode`` carries Newton limiting state (``_vlim``, seeded on
    first evaluation and thereafter advanced only by ``limit()``), so
    repeated ``i()`` calls at different x return *limited* values.  The
    generated element is pure -- no ``$limit`` equivalent in the DSL yet
    (hdl.md, capability map) -- which is exactly the difference this test
    pins: same equation, no limiting.
    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    from pycircuit.circuit.circuit import defaultepar
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    e1 = eh.DiodeHdl('p', 'n', IS=1e-13, vt=VT)
    e1.update_iparv()
    for v in (0.1, 0.3, 0.55):
        x = np.array([v, 0.0])
        e2 = Diode('p', 'n', IS=1e-13)      # fresh: _vlim seeds at THIS x
        e2.update_iparv()
        assert_allclose(np.asarray(e1.i(x), float),
                        np.asarray(e2.i(x), float), rtol=1e-12)
        assert_allclose(np.asarray(e1.G(x), float),
                        np.asarray(e2.G(x), float), rtol=1e-12)
        ## and the analytic value, so the test does not merely agree with
        ## another implementation of the same possible mistake
        assert_allclose(np.asarray(e1.i(x), float)[0],
                        1e-13 * (np.exp(v / VT) - 1), rtol=1e-12)

    ## The stateful-vs-pure difference itself, asserted rather than implied.
    e3 = Diode('p', 'n', IS=1e-13); e3.update_iparv()
    e3.i(np.array([0.1, 0.0]))                       # seeds _vlim = 0.1
    limited = float(np.asarray(e3.i(np.array([0.55, 0.0])), float)[0])
    pure = float(np.asarray(e1.i(np.array([0.55, 0.0])), float)[0])
    assert limited < pure / 100, (limited, pure)


def test_generated_jacobian_matches_finite_difference():
    """The symbolic G must be the derivative of the generated i -- checked
    by central differences, the test_element_jacobians discipline."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    el = eh.DiodeHdl('p', 'n', IS=1e-13)
    el.update_iparv()
    x = np.array([0.42, 0.1])
    eps = 1e-7
    J = np.zeros((2, 2))
    for k in range(2):
        dx = np.zeros(2); dx[k] = eps
        J[:, k] = (np.asarray(el.i(x + dx), float)
                   - np.asarray(el.i(x - dx), float)) / (2 * eps)
    assert_allclose(np.asarray(el.G(x), float), J, rtol=1e-5)


## ------------------------------------------------------------------------
## In-circuit equivalence


def _tran(c, node, tend=1e-3, dt=1e-6, **kw):
    import warnings
    tran = Transient(c, toolkit=numeric, **kw)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=dt, fixed_timestep=True)
    return (np.asarray(res.v(node).x[0], float),
            np.asarray(res.v(node).y, float))


def test_rc_lowpass_transient_equivalence():
    """An RC lowpass built from HDL elements reproduces the elements.py
    circuit's waveform -- HDL elements assemble, bias and integrate inside
    a real SubCircuit (the capability the old DSL lacked entirely)."""
    def build(use_hdl):
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VSin(na, gnd, va=1.0, freq=2e3)
        if use_hdl:
            c['R1'] = eh.RHdl(na, nb, r=1e3)
            c['C1'] = eh.CHdl(nb, gnd, c=1e-7)
        else:
            c['R1'] = R(na, nb, r=1e3)
            c['C1'] = C(nb, gnd, c=1e-7)
        return c

    t1, y1 = _tran(build(True), 'b', uic=True)
    t2, y2 = _tran(build(False), 'b', uic=True)
    assert_allclose(y1, y2, atol=1e-12)


def test_inductor_transient_equivalence():
    """RL circuit: LHdl's branch row is the elements.L row times -1; the
    waveforms must nevertheless agree exactly."""
    def build(use_hdl):
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VSin(na, gnd, va=1.0, freq=2e3)
        c['R1'] = R(na, nb, r=100.0)
        c['L1'] = (eh.LHdl if use_hdl else L)(nb, gnd, L=1e-3)
        return c

    t1, y1 = _tran(build(True), 'b', uic=True)
    t2, y2 = _tran(build(False), 'b', uic=True)
    assert_allclose(y1, y2, atol=1e-9)


def test_vcvs_dc_equivalence():
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    ni, no = c.add_node('i'), c.add_node('o')
    c['vs'] = VS(ni, gnd, v=0.25)
    c['E'] = eh.VCVSHdl(ni, gnd, no, gnd, g=4.0)
    c['Rl'] = R(no, gnd, r=1e3)
    res = DC(c, toolkit=numeric).solve()
    assert abs(res.v('o') - 1.0) < 1e-9


def test_hdl_source_transient():
    """VSinHdl -- a time-dependent V-contribution -- drives a divider to
    the analytic waveform, and its generated dudt matches the numeric
    derivative of its u."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['vs'] = eh.VSinHdl(na, gnd, va=2.0, freq=1e3)
    c['R1'] = R(na, nb, r=1e3)
    c['R2'] = R(nb, gnd, r=1e3)
    t, y = _tran(c, 'b', tend=2e-3, dt=2e-6, uic=True)
    assert_allclose(y[1:], np.sin(2 * np.pi * 1e3 * t[1:]), atol=1e-9)

    el = eh.VSinHdl(0, 1, va=2.0, freq=1e3)
    el.update_iparv()
    h = 1e-9
    for t0 in (0.0, 1.3e-4):
        num = (el.u(t0 + h) - el.u(t0 - h)) / (2 * h)
        assert_allclose(el.dudt(t0), num, rtol=1e-5, atol=1e-4)


def test_idt_transient_equivalence():
    def build(use_hdl):
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        nin, nout = c.add_node('in'), c.add_node('out')
        c['vin'] = VS(nin, gnd, v=1.0)
        c['R1'] = R(nout, gnd, r=1e3)
        c['X'] = (eh.IdtHdl if use_hdl else Idt)(nin, gnd, nout, gnd)
        return c

    from pycircuit.circuit.integrator import EulerIntegrator
    t1, y1 = _tran(build(True), 'out', tend=0.5, dt=1e-2, uic=True,
                   integrator=EulerIntegrator())
    t2, y2 = _tran(build(False), 'out', tend=0.5, dt=1e-2, uic=True,
                   integrator=EulerIntegrator())
    assert_allclose(y1, y2, atol=1e-12)


def test_idtmod_full_treatment():
    """The one-line IdtmodHdl gets everything the hand-built Idtmod has:
    congruence-correct wrapped output, DC pinned from ic WITHOUT uic, and
    a bounded state through the generated periodic_states declaration."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nin, nout = c.add_node('in'), c.add_node('out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['X'] = eh.IdtmodHdl(nin, gnd, nout, gnd, modulus=1.0, offset=-0.5,
                          ic=1.7)

    ## DC pin: wrap(1.7) into [-0.5, 0.5) is -0.3.
    res = DC(c, toolkit=numeric).solve()
    assert abs(res.v('out') - (-0.3)) < 1e-9

    ## Transient from the pin, several wraps, congruence + bounded state.
    from pycircuit.circuit.integrator import EulerIntegrator
    t, y = _tran(c, 'out', tend=3.0, dt=1e-2, integrator=EulerIntegrator())
    ref = ((t + (-0.3) + 0.5) % 1.0) - 0.5
    d = np.abs(y[1:] - ref[1:])
    d = np.minimum(d, 1.0 - d)
    assert d.max() < 1e-9
    names = [str(nd) for nd in c.nodes]
    srow = [i for i, nm in enumerate(names) if '_state0' in nm][0]
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        full = Transient(c, toolkit=numeric,
                         integrator=EulerIntegrator()).solve(
            tend=3.0, timestep=1e-2, fixed_timestep=True)
    s = np.asarray(full.x[srow], float)
    assert s.max() - s.min() < 1.5    # 3 wraps happened; unbounded would span 3


## ------------------------------------------------------------------------
## Language-level features


def test_noise_contribution_generates_CY():
    """white_noise in an I-contribution reproduces elements.R's thermal
    CY stamp; flicker_noise carries the 1/f shape."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, white_noise,
                                       flicker_noise, KBOLTZMANN, TEMP)

    class RNoisy(Behavioural):
        instparams = [Parameter(name='r', desc='R', unit='ohm', default=1e3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, b.V / r),
                    Contribution(b.I, white_noise(4 * KBOLTZMANN * TEMP / r)))

    pycircuit.circuit.circuit.default_toolkit = numeric
    e1 = RNoisy('p', 'n', r=2e3)
    e2 = R('p', 'n', r=2e3)
    e1.update_iparv(); e2.update_iparv()
    x = np.zeros(2)
    from pycircuit.circuit.circuit import defaultepar
    assert_allclose(np.asarray(e1.CY(x, 1.0), float),
                    np.asarray(e2.CY(x, 1.0), float), rtol=1e-12)

    class INoise1f(Behavioural):
        instparams = [Parameter(name='kf', desc='k', unit='', default=1e-12)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, flicker_noise(kf, 1)),)

    e3 = INoise1f('p', 'n', kf=1e-12)
    e3.update_iparv()
    w1, w2 = 2 * np.pi * 10.0, 2 * np.pi * 1000.0
    cy1 = np.asarray(e3.CY(x, w1), float)[0, 0]
    cy2 = np.asarray(e3.CY(x, w2), float)[0, 0]
    assert_allclose(cy1 / cy2, 100.0, rtol=1e-12)   # 1/f: x100 at f/100


def test_named_noise_sources_are_correlated():
    """Two contributions sharing a NAME are ONE fluctuation reaching two
    branches (LRM 2.4 sec. 4.5.16), so their ``CY`` block is rank one.

    The distinction this pins is not decorative.  Two INDEPENDENT
    sources of the same powers give the same diagonal and a zero
    off-diagonal, and every circuit in which the two paths interfere
    then gets the wrong answer -- which is precisely the regime the
    feature exists for.
    """
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, white_noise)

    class Pair(Behavioural):
        terminals = ['a', 'b', 'c']
        instparams = [Parameter(name='s', desc='p', unit='', default=1e-20),
                      Parameter(name='k', desc='k', unit='', default=3.0),
                      Parameter(name='u', desc='u', unit='', default=0.0)]

        @staticmethod
        def analog(a, b, c):
            return (Contribution(Branch(a, c).I, white_noise(s, 'g')),  # noqa
                    Contribution(Branch(b, c).I,
                                 k * white_noise(s, 'g')),              # noqa
                    Contribution(Branch(a, c).I, white_noise(u)))       # noqa

    pycircuit.circuit.circuit.default_toolkit = numeric
    e = Pair('a', 'b', 'c', s=1e-20, k=3.0, u=4e-20)
    e.update_iparv()
    CY = np.asarray(e.CY(np.zeros(e.n), 2 * np.pi), float) / 1e-20

    ## The group's stamp vector is `1*(a-c) + 3*(b-c) = (1, 3, -4)`, so
    ## its block is `w w^T`; the lone source adds its own 2x2 pattern.
    w = np.array([1.0, 3.0, -4.0])
    lone = np.zeros((3, 3))
    lone[0, 0] = lone[2, 2] = 4.0
    lone[0, 2] = lone[2, 0] = -4.0
    assert_allclose(CY, np.outer(w, w) + lone, rtol=1e-12, atol=1e-12)

    ## A scale factor multiplies the AMPLITUDE, so it appears SQUARED in
    ## the density -- the b-b entry is `k^2`, not `k`.
    assert_allclose(CY[1, 1], 9.0, rtol=1e-12)

    ## Rank one, and positive semi-definite: a correlation matrix that
    ## is neither is not a noise description at all.
    grp = CY - lone
    assert np.linalg.matrix_rank(grp, tol=1e-9) == 1
    assert np.all(np.linalg.eigvalsh(CY) > -1e-9)


def test_unnamed_noise_sources_are_independent():
    """The same two contributions WITHOUT a name share nothing, so the
    off-diagonal between the two branches is zero."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, white_noise)

    class Loose(Behavioural):
        terminals = ['a', 'b', 'c']
        instparams = [Parameter(name='s', desc='p', unit='', default=1e-20)]

        @staticmethod
        def analog(a, b, c):
            return (Contribution(Branch(a, c).I, white_noise(s)),   # noqa
                    Contribution(Branch(b, c).I, white_noise(s)))   # noqa

    pycircuit.circuit.circuit.default_toolkit = numeric
    e = Loose('a', 'b', 'c', s=1e-20)
    e.update_iparv()
    CY = np.asarray(e.CY(np.zeros(e.n), 2 * np.pi), float) / 1e-20
    assert_allclose(CY[0, 1], 0.0, atol=1e-12)
    assert_allclose(np.diag(CY), [1.0, 1.0, 2.0], rtol=1e-12)


def test_a_correlation_group_must_share_one_power():
    """Members of a group are ONE fluctuation, so two different powers
    under one name is a contradiction rather than a shorthand -- with
    two powers the cross term is `sqrt(S1 S2)` up to a SIGN that nothing
    in the source text determines.  Refused, with the fix named."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, white_noise)

    with pytest.raises(NotImplementedError, match='scale factor'):
        class Mismatch(Behavioural):
            terminals = ['a', 'b', 'c']
            instparams = [Parameter(name='s', desc='p', unit='',
                                    default=1e-20)]

            @staticmethod
            def analog(a, b, c):
                return (Contribution(Branch(a, c).I,
                                     white_noise(s, 'g')),          # noqa
                        Contribution(Branch(b, c).I,
                                     white_noise(2 * s, 'g')))      # noqa

        Mismatch('a', 'b', 'c').update_iparv()


def test_a_noise_name_is_not_a_parameter():
    """The name is inert: it must not become a `Symbol` that the
    parameter machinery then goes looking for."""
    from pycircuit.circuit.hdl import white_noise
    import sympy
    app = white_noise(sympy.Symbol('p'), 'igid')
    assert app.free_symbols == {sympy.Symbol('p')}
    assert app.noise_name == 'igid'
    assert white_noise(sympy.Symbol('p')).noise_name is None


def test_ddx_and_limexp():
    from pycircuit.circuit.hdl import ddx, limexp, Branch, Node
    import sympy
    a, b = Node('a'), Node('b')
    br = Branch(a, b)
    expr = 3 * br.V ** 2 + sympy.exp(br.V)
    d = ddx(expr, br.V)
    assert sympy.simplify(d - (6 * br.V + sympy.exp(br.V))) == 0

    f = sympy.lambdify(sympy.Symbol('z'),
                       limexp(sympy.Symbol('z')), modules='numpy')
    assert_allclose(f(1.0), np.e, rtol=1e-12)
    assert np.isfinite(f(1000.0))                     # linear continuation
    h = 1e-6
    slope = (f(80.0 + h) - f(80.0 - h)) / (2 * h)     # C1 at the seam
    assert_allclose(slope, np.exp(80.0), rtol=1e-3)


def test_rejects_state_dependent_ddt_scaling():
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, ddt)
    with pytest.raises(NotImplementedError, match='ddt'):
        class Bad(Behavioural):
            instparams = []

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return (Contribution(b.I, b.V * ddt(b.V)),)


def test_scaled_ddt_and_nonlinear_charge():
    """The two historic crashers: c*ddt(v) (misclassified) and a scaled
    nonlinear charge (printer name collision) now compile and give the
    right C."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Parameter, ddt)
    import sympy

    class C2(Behavioural):
        instparams = [Parameter(name='c0', desc='c', unit='F', default=2.0)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, c0 * ddt(b.V)),)      # noqa: F821

    class CTanh(Behavioural):
        instparams = [Parameter(name='c0', desc='c', unit='F', default=3.0)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, ddt(c0 * sympy.tanh(b.V))),)  # noqa: F821

    pycircuit.circuit.circuit.default_toolkit = numeric
    e = C2('p', 'n'); e.update_iparv()
    assert_allclose(np.asarray(e.C(np.zeros(2)), float),
                    [[2.0, -2.0], [-2.0, 2.0]], atol=1e-15)
    e = CTanh('p', 'n'); e.update_iparv()
    assert_allclose(np.asarray(e.C(np.zeros(2)), float),
                    [[3.0, -3.0], [-3.0, 3.0]], atol=1e-12)


def test_eval_pures_jax_parity():
    """Generated eval_i_pure/eval_q_pure (sympy jax printer) agree with the
    numeric path -- the admission ticket to the vmap groups."""
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit
    el = eh.DiodeHdl('p', 'n', IS=1e-13)
    el.update_iparv()
    assert hasattr(eh.DiodeHdl, 'eval_i_pure')
    x = np.array([0.31, 0.05])
    params = {'IS': 1e-13, 'vt': el.iparv.vt}
    from pycircuit.circuit.circuit import defaultepar
    got = np.asarray(eh.DiodeHdl.eval_i_pure(
        jaxtoolkit.array(x), params, defaultepar, jaxtoolkit), float)
    want = np.asarray(el.i(x), float)
    assert_allclose(got, want, rtol=1e-12)


def test_constant_stamp_cache_follows_parameters():
    """The x-independent G/C fast path (benchmarks/hdl_overhead.py) must
    not serve a stale matrix: the update() observer drops it whenever
    parameters move, including late resolution inside a SubCircuit."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    el = eh.RHdl('p', 'n', r=1e3)
    el.update_iparv()
    x = np.zeros(2)
    assert_allclose(np.asarray(el.G(x), float)[0, 0], 1e-3, rtol=1e-12)
    el.ipar.set(r=250.0)          # notify -> update() -> cache dropped
    el.update_iparv()
    assert_allclose(np.asarray(el.G(x), float)[0, 0], 4e-3, rtol=1e-12)

    ## Hierarchical resolution: the value is a STRING expression that only
    ## becomes a number at update_iparv, and the generated stamp must
    ## follow it (the defect the old DSL had by reading self.ipar).
    from pycircuit.utilities.param import Parameter, ParameterDict
    glob = ParameterDict(Parameter(name='rbase', desc='base R', unit='ohm'))
    glob.set(rbase=500.0)
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['R1'] = eh.RHdl(na, nb, r='2*rbase')
    c.update_iparv(globalparams=glob, ignore_errors=True)
    assert_allclose(np.asarray(c['R1'].G(x), float)[0, 0], 1e-3, rtol=1e-12)


def test_hdl_element_ac_has_no_bias_leak():
    """A quantity-free constant belongs to the device characteristic (i),
    not to the independent-source vector -- so a nonlinear HDL element
    injects NO spurious AC excitation (the classic ABM bias leak)."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    el = eh.DiodeHdl('p', 'n', IS=1e-13)
    el.update_iparv()
    assert_allclose(np.asarray(el.u(0.0, analysis='ac'), float),
                    np.zeros(2), atol=0)
    assert_allclose(np.asarray(el.u(0.0), float), np.zeros(2), atol=0)
    ## and the -IS constant really is in i:
    assert_allclose(np.asarray(el.i(np.zeros(2)), float), np.zeros(2),
                    atol=1e-30)
    assert float(np.asarray(el.i(np.array([-1.0, 0.0])), float)[0]) < 0


def test_internal_node_and_internal_branch():
    """A V-contribution on a branch that spans an INTERNAL node: the
    generated element must register the internal node, place the branch
    current after it (the update_node_map order), and solve in-circuit."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Node,
                                       Contribution, ddt)

    class SeriesRL(Behavioural):
        instparams = [Parameter(name='rr', desc='R', unit='ohm',
                                default=100.0),
                      Parameter(name='LL', desc='L', unit='H', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            mid = Node('mid')
            b1, b2 = Branch(plus, mid), Branch(mid, minus)
            return (Contribution(b1.I, b1.V / rr),          # noqa: F821
                    Contribution(b2.V, ddt(LL * b2.I)))     # noqa: F821

    pycircuit.circuit.circuit.default_toolkit = numeric
    el = SeriesRL('p', 'n')
    el.update_iparv()
    assert el.n == 4                       # 2 terminals + mid + branch i
    assert 'mid' in [str(nd) for nd in el.nodes]

    c = SubCircuit()
    na = c.add_node('a')
    c['vs'] = VS(na, gnd, v=2.0)
    c['X'] = SeriesRL(na, gnd, rr=100.0, LL=1e-3)
    res = DC(c, toolkit=numeric).solve()
    assert abs(float(res.v('a')) - 2.0) < 1e-9
    ## At DC the inductor is a short, so the series current is v/r.
    x = np.asarray(res.x, float)
    assert np.any(np.isclose(np.abs(x), 0.02, rtol=1e-6)), x


## ------------------------------------------------------------------------
## limexp in circuit, and its relationship with PCNR limiting.


def _limexp_diode_class(fn=None):
    """A diode element whose exponential is limexp (default) or plain exp."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, vt)
    f = limexp if fn is None else fn

    class _D(Behavioural):
        instparams = [Parameter(name='IS', desc='Sat current', unit='A',
                                default=1e-13)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * (f(b.V / vt()) - 1))  # noqa: F821
    return _D


def test_limexp_inside_an_element_and_in_circuit():
    """limexp must be usable in a real element, not just as a standalone
    sympy expression.

    Regression: the DSL's atoms carried no `real` assumption, so sympy
    refused to build the `Piecewise` condition at all ("Invalid comparison
    of non-real") and EVERY limexp element -- indeed every conditional
    model -- died at class creation.  The standalone-function test did not
    catch it because no Quantity was involved.
    """
    pycircuit.circuit.circuit.default_toolkit = numeric
    D = _limexp_diode_class()
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['vs'] = VS(na, gnd, v=5.0)
    c['R'] = R(na, nb, r=1e3)
    c['D'] = D(nb, gnd, IS=1e-13)
    res = DC(c, toolkit=numeric).solve()
    v = float(res.v('b'))
    ## Physical forward drop, and it must satisfy the diode equation.
    assert 0.5 < v < 0.8, v
    VT = numeric.kboltzmann * 300.0 / numeric.qelectron
    assert_allclose(1e-13 * (np.exp(v / VT) - 1), (5.0 - v) / 1e3, rtol=1e-6)


def test_limexp_keeps_residual_finite_where_exp_overflows():
    """What limexp buys: identical values in the physical region, and a
    FINITE residual and Jacobian where plain exp returns inf -- which is
    what keeps Newton alive when an iterate wanders."""
    import sympy
    pycircuit.circuit.circuit.default_toolkit = numeric
    lim = _limexp_diode_class()('p', 'n', IS=1e-13)
    raw = _limexp_diode_class(sympy.exp)('p', 'n', IS=1e-13)
    lim.update_iparv(); raw.update_iparv()

    x = np.array([0.7, 0.0])                       # physical: identical
    assert_allclose(np.asarray(lim.i(x), float),
                    np.asarray(raw.i(x), float), rtol=1e-12)

    x = np.array([50.0, 0.0])                      # wandered iterate
    with np.errstate(over='ignore'):
        assert not np.all(np.isfinite(np.asarray(raw.i(x), float)))
        assert not np.all(np.isfinite(np.asarray(raw.G(x), float)))
    assert np.all(np.isfinite(np.asarray(lim.i(x), float)))
    assert np.all(np.isfinite(np.asarray(lim.G(x), float)))


def test_generated_pcnr_protocol_matches_hand_written_diode():
    """The DSL generates PCNR participation: an element whose whole current
    is an exponential of one branch voltage gets `pcnr_junctions`,
    `pcnr_i`, `pcnr_didv` and `pcnr_limit` emitted from its own equation --
    and they must equal the hand-written Diode's to the last digit."""
    import sympy
    from pycircuit.circuit.circuit import defaultepar
    pycircuit.circuit.circuit.default_toolkit = numeric

    for fn in (sympy.exp, None):            # plain exp and limexp alike
        D = _limexp_diode_class(fn)
        assert D.pcnr_junctions == ((0, 1),), D.pcnr_junctions
        for v in (0.2, 0.6, 0.75):
            got_i = np.asarray(D.pcnr_i(v, {'IS': 1e-13}, defaultepar,
                                        numeric), float)
            ref_i = np.asarray(Diode.pcnr_i(v, {'IS': 1e-13}, defaultepar,
                                            numeric), float)
            assert_allclose(got_i, ref_i, rtol=1e-12)
            got_g = np.asarray(D.pcnr_didv(v, {'IS': 1e-13}, defaultepar,
                                           numeric), float)
            ref_g = np.asarray(Diode.pcnr_didv(v, {'IS': 1e-13}, defaultepar,
                                               numeric), float)
            assert_allclose(got_g, ref_g, rtol=1e-12)

    ## The generated limiter is pnjlim with the scale read off the
    ## expression, so it reproduces the reference limiter exactly.
    from pycircuit.circuit._limiting import _pnjlim
    D = _limexp_diode_class()
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    for vnew, vold in ((0.9, 0.6), (5.0, 0.7), (0.3, 0.2)):
        got = float(D.pcnr_limit(vnew, vold, {'IS': 1e-13}, defaultepar,
                                 numeric))
        ref = _pnjlim(vnew, vold, VT, 1e-13, numeric)
        assert_allclose(got, ref, rtol=1e-12)


def test_pcnr_handles_limiting_for_generated_elements():
    """With PCNR enabled it does the limiting, for generated elements as
    for hand-written ones -- the operating point is the pcnr=False one,
    and it no longer depends on whether the model used limexp or a raw
    exp (PCNR keeps the iterate in range either way)."""
    import sympy
    from pycircuit.circuit.pcnr import pcnr_junctions
    pycircuit.circuit.circuit.default_toolkit = numeric

    def solve(cls, pcnr):
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=5.0)
        c['R'] = R(na, nb, r=1e3)
        c['D'] = cls(nb, gnd, IS=1e-13)
        c.update_iparv()
        return (float(DC(c, toolkit=numeric, pcnr=pcnr).solve().v('b')),
                len(pcnr_junctions(c)))

    v_ref, _ = solve(Diode, False)
    for fn in (None, sympy.exp):
        D = _limexp_diode_class(fn)
        v_off, n_off = solve(D, False)
        v_on, n_on = solve(D, True)
        assert n_on == 1 and n_off == 1        # the DSL element participates
        assert_allclose(v_off, v_ref, rtol=1e-6)
        assert_allclose(v_on, v_ref, rtol=1e-5)

    ## Under PCNR the limexp and raw-exp models agree exactly: the
    ## limiting is PCNR's, not the model's.
    v_lim, _ = solve(_limexp_diode_class(None), True)
    v_raw, _ = solve(_limexp_diode_class(sympy.exp), True)
    assert_allclose(v_lim, v_raw, rtol=1e-12)


def test_pcnr_not_generated_for_unqualified_elements():
    """Participation is claimed only where it is true: an element with
    charge, with a state, or whose current is not one branch voltage's
    exponential declares no junction and is solved the ordinary way."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ddt, idt)

    class WithCharge(Behavioural):
        instparams = [Parameter(name='cc', desc='c', unit='F', default=1e-9)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, ddt(cc * b.V))          # noqa: F821

    class Polynomial(Behavioural):
        instparams = [Parameter(name='k', desc='k', unit='', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, k * b.V ** 3)           # noqa: F821

    for cls in (WithCharge, Polynomial, eh.RHdl, eh.LHdl, eh.IdtmodHdl):
        assert not hasattr(cls, 'pcnr_junctions'), cls.__name__


def test_piecewise_conditional_model():
    """A Piecewise (Verilog-A's if/else in an analog block) compiles and
    gives the right branch and slope on each side -- the same assumption
    fix limexp needed."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution)

    class SoftClip(Behavioural):
        instparams = [Parameter(name='g0', desc='g', unit='S', default=1e-3),
                      Parameter(name='vk', desc='knee', unit='V',
                                default=0.5)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            i = sympy.Piecewise((g0 * b.V, b.V < vk),          # noqa: F821
                                (g0 * vk, True))               # noqa: F821
            return Contribution(b.I, i)

    pycircuit.circuit.circuit.default_toolkit = numeric
    el = SoftClip('p', 'n', g0=1e-3, vk=0.5)
    el.update_iparv()
    lo = np.array([0.2, 0.0]); hi = np.array([0.9, 0.0])
    assert_allclose(float(np.asarray(el.i(lo), float)[0]), 2e-4, rtol=1e-12)
    assert_allclose(float(np.asarray(el.i(hi), float)[0]), 5e-4, rtol=1e-12)
    assert_allclose(float(np.asarray(el.G(lo), float)[0, 0]), 1e-3, rtol=1e-12)
    assert abs(float(np.asarray(el.G(hi), float)[0, 0])) < 1e-15
    assert el.linear is False        # a Piecewise is not a linear element


def test_subclass_changing_params_without_analog_is_refused():
    """Compilation is keyed on `analog` in the class body, so a subclass
    that changes only instparams would silently run its parent's compiled
    code against a different parameter list.  Refused loudly; redeclaring
    analog() recompiles and works."""
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution

    class Base(Behavioural):
        instparams = [Parameter(name='ga', desc='g', unit='S', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, ga * b.V)            # noqa: F821

    class BadChild(Base):
        instparams = [Parameter(name='gb', desc='g', unit='S', default=2e-3)]

    class GoodChild(Base):
        instparams = [Parameter(name='ga', desc='g', unit='S', default=2e-3)]
        analog = staticmethod(Base.analog.__func__
                              if hasattr(Base.analog, '__func__')
                              else Base.analog)

    pycircuit.circuit.circuit.default_toolkit = numeric
    with pytest.raises(TypeError, match='instparams'):
        BadChild('p', 'n')
    el = GoodChild('p', 'n')
    el.update_iparv()
    assert_allclose(np.asarray(el.G(np.zeros(2)), float)[0, 0], 2e-3,
                    rtol=1e-12)


## ------------------------------------------------------------------------
## Phase B (hdl.md sec. 9): the model surface -- $param_given, parameter
## ranges, aliasparam, node collapse, AC excitation.


def test_param_given_selects_a_formulation():
    """`$param_given` is a runtime value, so one compiled class serves
    instances that supplied the parameter and instances that did not --
    which is the whole point, since compilation happens once per class."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       param_given)

    class Rsel(Behavioural):
        instparams = [Parameter(name='g1', desc='g', unit='S', default=1e-3),
                      Parameter(name='g2', desc='g', unit='S', default=5e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            g = sympy.Piecewise((g2, param_given('g2') > 0.5),   # noqa: F821
                                (g1, True))                      # noqa: F821
            return Contribution(b.I, g * b.V)

    pycircuit.circuit.circuit.default_toolkit = numeric
    x = np.array([1.0, 0.0])

    default = Rsel('p', 'n'); default.update_iparv()
    assert_allclose(float(np.asarray(default.i(x), float)[0]), 1e-3,
                    rtol=1e-12)

    ## Same value as the default, but GIVEN -- givenness is not "differs
    ## from the default", which is exactly why the operator exists.
    same = Rsel('p', 'n', g2=5e-3); same.update_iparv()
    assert_allclose(float(np.asarray(same.i(x), float)[0]), 5e-3, rtol=1e-12)

    other = Rsel('p', 'n', g2=2e-3); other.update_iparv()
    assert_allclose(float(np.asarray(other.i(x), float)[0]), 2e-3, rtol=1e-12)


def test_parameter_range_is_enforced():
    from pycircuit.utilities.param import Parameter

    p = Parameter(name='w', desc='width', unit='m', default=1e-6,
                  minval=0.0, maxval=1e-3)
    from pycircuit.utilities.param import ParameterDict
    d = ParameterDict(p)
    d.set(w=5e-4)
    assert d.get('w') == 5e-4
    with pytest.raises(ValueError, match='outside its declared range'):
        d.set(w=2e-3)
    with pytest.raises(ValueError, match='outside its declared range'):
        d.set(w=-1.0)


def test_aliasparam():
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution

    class D(Behavioural):
        instparams = [Parameter(name='IS', desc='Is', unit='A',
                                default=1e-13)]
        aliasparams = {'isat': 'IS', 'js': 'IS'}

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, IS * b.V)          # noqa: F821

    pycircuit.circuit.circuit.default_toolkit = numeric
    x = np.array([1.0, 0.0])
    for kw in ({'IS': 3e-13}, {'isat': 3e-13}, {'js': 3e-13}):
        el = D('p', 'n', **kw); el.update_iparv()
        assert_allclose(float(np.asarray(el.i(x), float)[0]), 3e-13,
                        rtol=1e-12)
    with pytest.raises(ValueError, match='alias'):
        D('p', 'n', IS=1e-13, isat=2e-13)


def test_node_collapse_removes_the_internal_node():
    """`V(a,b) <+ 0` on an internal node merges it away, so the optional
    series resistance idiom costs nothing when the resistance is absent."""
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Node,
                                       Contribution)

    def build(collapse):
        class _D(Behavioural):
            instparams = [Parameter(name='IS', desc='Is', unit='A',
                                    default=1e-13)]

            @staticmethod
            def analog(plus, minus):
                ia = Node('ia')
                b_rs, b_j = Branch(plus, ia), Branch(ia, minus)
                junction = Contribution(
                    b_j.I, IS * (sympy.exp(b_j.V / 0.02585) - 1))  # noqa
                if collapse:
                    return (junction, Contribution(b_rs.V, 0))
                return (junction,
                        Contribution(b_rs.I, b_rs.V / 1e-9))   # tiny R
        return _D

    pycircuit.circuit.circuit.default_toolkit = numeric
    col = build(True)('p', 'n'); col.update_iparv()
    unc = build(False)('p', 'n'); unc.update_iparv()

    ## Collapsed: the internal node is gone entirely.
    assert col.n == 2, col.n
    assert 'ia' not in [str(nd) for nd in col.nodes]
    ## Uncollapsed: the node (and its row) survives.
    assert unc.n == 3, unc.n

    ## And the collapsed element is the plain junction.
    v = 0.3
    got = float(np.asarray(col.i(np.array([v, 0.0])), float)[0])
    assert_allclose(got, 1e-13 * (np.exp(v / 0.02585) - 1), rtol=1e-9)


def test_ac_stim_drives_a_small_signal_analysis():
    """`ac_stim` is live in AC and identically zero elsewhere -- so a
    behavioural element can BE a small-signal source."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       ac_stim)
    from pycircuit.circuit.analysis_ss import AC

    class IAC(Behavioural):
        instparams = []

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, ac_stim(1.0))

    pycircuit.circuit.circuit.default_toolkit = numeric
    el = IAC('p', 'n'); el.update_iparv()
    ## Zero outside AC ...
    assert_allclose(np.asarray(el.u(0.0), float), np.zeros(2), atol=0)
    ## ... and a unit current in AC.
    uac = np.asarray(el.u(0.0, analysis='ac'), complex)
    assert_allclose(uac, np.array([1.0, -1.0], dtype=complex), atol=1e-15)

    ## In a circuit: 1 A into 1 kohm is 1 kV.
    c = SubCircuit()
    na = c.add_node('a')
    c['I1'] = IAC(gnd, na)
    c['R1'] = R(na, gnd, r=1e3)
    res = AC(c, toolkit=numeric).solve(freqs=1e3)
    assert_allclose(abs(complex(res.v('a'))), 1e3, rtol=1e-9)


## ------------------------------------------------------------------------
## Phase A2 (hdl.md sec. 9): $limit, the fallback where PCNR cannot go.


def _limited_diode(use_limit):
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, limit_pnj, vt)

    class _D(Behavioural):
        instparams = [Parameter(name='IS', desc='Is', unit='A',
                                default=1e-13)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            v = limit_pnj(b.V, IS, vt()) if use_limit else b.V  # noqa: F821
            return Contribution(b.I, IS * (limexp(v / vt()) - 1))  # noqa
    return _D


def test_limit_pnj_is_generated_only_when_asked():
    assert hasattr(_limited_diode(True), 'limit')
    assert not hasattr(_limited_diode(False), 'limit')


def test_limit_pnj_compresses_a_wild_step():
    """The point of a limiter: an undamped Newton step from a lightly
    biased point overshoots enormously, and the model is then evaluated
    somewhere it has no business being."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    el = _limited_diode(True)('p', 'n', IS=1e-13)
    el.update_iparv()
    x0 = np.array([0.1, 0.0])
    x = np.array([5.0, 0.0])                    # the wild step
    out = el.limit(x, x0)
    v = float(out[0] - out[1])
    assert 0.1 < v < 0.5, v                     # compressed, not clamped
    ## State-free: the input is not mutated, a limited COPY comes back.
    assert x[0] == 5.0
    ## A small step passes through untouched (the paper's 2*VT escape).
    small = el.limit(np.array([0.11, 0.0]), x0)
    assert_allclose(float(small[0] - small[1]), 0.11, rtol=1e-12)


def test_limit_pnj_changes_the_path_not_the_solution():
    """The defining property of limiting: it moves where the next
    Jacobian is taken, never where the solution is."""
    pycircuit.circuit.circuit.default_toolkit = numeric

    def solve(cls, vsrc):
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=vsrc)
        c['R'] = R(na, nb, r=1e3)
        c['D'] = cls(nb, gnd, IS=1e-13)
        return float(DC(c, toolkit=numeric).solve().v('b'))

    for vsrc in (5.0, 100.0):
        plain = solve(_limited_diode(False), vsrc)
        limited = solve(_limited_diode(True), vsrc)
        assert_allclose(limited, plain, rtol=1e-7)
        ## and it is a real forward-biased junction, not a trivial point
        assert 0.5 < plain < 0.9


## ------------------------------------------------------------------------
## Phase C1 (hdl.md sec. 9): @cross events.


def _comparator(use_cross):
    import sympy
    from pycircuit.utilities.param import Parameter
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Cross)

    class _Comp(Behavioural):
        instparams = [Parameter(name='vref', desc='ref', unit='V',
                                default=0.0),
                      Parameter(name='gain', desc='g', unit='', default=50.0)]

        @staticmethod
        def analog(inp, inn, outp, outn):
            bi, bo = Branch(inp, inn), Branch(outp, outn)
            out = sympy.tanh(gain * (bi.V - vref))          # noqa: F821
            sts = (Contribution(bo.V, out),)
            if use_cross:
                sts = sts + (Cross(bi.V - vref, 0),)        # noqa: F821
            return sts
    return _Comp


def _run_comparator(use_cross):
    import warnings
    from pycircuit.circuit.elements import VSin
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    ni, no = c.add_node('i'), c.add_node('o')
    c['vs'] = VSin(ni, gnd, va=1.0, freq=1e3)
    c['X'] = _comparator(use_cross)(ni, gnd, no, gnd, vref=0.0)
    c['Rl'] = R(no, gnd, r=1e3)
    tran = Transient(c, toolkit=numeric, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=2e-3, timestep=1e-6)
    t = np.asarray(res.v('o').x[0], float)
    edges = (0.5e-3, 1.0e-3, 1.5e-3)          # zero crossings of the sine
    nearest = [float(np.min(np.abs(t - te))) for te in edges]
    return res.statistics, nearest


def test_cross_lands_timepoints_on_the_crossing():
    """A declared crossing becomes a breakpoint, so the solver puts a
    timepoint on the corner instead of discovering it by rejection.

    Measured at introduction: without ``Cross`` the nearest sample to an
    edge is 5.5e-6 to 1.4e-5 s away; with it, 5.4e-7 to 1.9e-6 s -- about
    an order of magnitude closer, for ONE extra accepted step and no
    rejections.  It is a first-order prediction, so it brackets the
    corner rather than landing exactly on it, which is all the
    break-step machinery needs.
    """
    s_off, near_off = _run_comparator(False)
    s_on, near_on = _run_comparator(True)

    assert s_off.breakpoints_hit == 0
    assert s_on.breakpoints_hit >= 3            # one per edge crossed
    assert max(near_on) < max(near_off) / 3.0, (near_on, near_off)
    ## and it does not cost a step explosion
    assert s_on.accepted_steps < s_off.accepted_steps * 1.5
    assert s_on.rejected_steps <= s_off.rejected_steps + 2


def test_cross_predictor_contract():
    """`next_event` must be inf before two accepted points exist, must
    respect the requested direction, and must return strictly future
    times (the contract SubCircuit.next_event relies on)."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       Cross)
    from pycircuit.utilities.param import Parameter

    class Rising(Behavioural):
        instparams = [Parameter(name='g', desc='g', unit='S', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return (Contribution(b.I, g * b.V),        # noqa: F821
                    Cross(b.V, +1))

    pycircuit.circuit.circuit.default_toolkit = numeric
    el = Rising('p', 'n')
    el.update_iparv()
    assert el.next_event(0.0) == np.inf            # nothing accepted yet
    el.accept_step(0.0, np.array([-2.0, 0.0]))
    assert el.next_event(0.0) == np.inf            # only one point
    el.accept_step(1.0, np.array([-1.0, 0.0]))     # rising 1 V/s
    assert_allclose(el.next_event(1.0), 2.0, rtol=1e-9)   # reaches 0 at t=2

    ## Falling motion must NOT trigger a rising-only crossing.
    el.reset_state()
    el.accept_step(0.0, np.array([2.0, 0.0]))
    el.accept_step(1.0, np.array([1.0, 0.0]))
    assert el.next_event(1.0) == np.inf


## ------------------------------------------------------------------------
## Phase C2 (hdl.md sec. 9): the symbolic toolkit.


def test_symbolic_toolkit_gets_exact_expressions():
    """Under a symbolic toolkit the element must return sympy, exactly --
    not the output of a numpy lambda that happens to survive symbols by
    duck typing.  The difference is invisible for plain arithmetic and
    fatal the moment the expression contains a `floor` or a `Piecewise`
    (the wrap of an idtmod, the branch of a limexp)."""
    import sympy
    from pycircuit.circuit.toolkit import symbolic
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       limexp, ddt)
    from pycircuit.utilities.param import Parameter

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = symbolic

        class RC(Behavioural):
            instparams = [Parameter(name='g', desc='g', unit='S',
                                    default=1e-3),
                          Parameter(name='c', desc='c', unit='F',
                                    default=1e-9)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, g * b.V + ddt(c * b.V))  # noqa

        el = RC('p', 'n', toolkit=symbolic)
        el.update_iparv()
        v1, v2 = sympy.symbols('v1 v2')
        i = el.i([v1, v2])
        assert sympy.simplify(i[0] - 1e-3 * (v1 - v2)) == 0
        q = el.q([v1, v2])
        assert sympy.simplify(q[0] - 1e-9 * (v1 - v2)) == 0
        G = el.G([v1, v2])
        assert sympy.simplify(G[0, 0] - 1e-3) == 0
        C = el.C([v1, v2])
        assert sympy.simplify(C[0, 0] - 1e-9) == 0

        ## A Piecewise survives symbolically -- this is what the numpy
        ## path cannot do.
        class Lim(Behavioural):
            instparams = [Parameter(name='IS', desc='Is', unit='A',
                                    default=1e-13)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I, IS * limexp(b.V))    # noqa: F821

        el2 = Lim('p', 'n', toolkit=symbolic)
        el2.update_iparv()
        expr = el2.i([v1, v2])[0]
        assert expr.has(sympy.Piecewise)
        ## and it evaluates to the right number when the symbols are fixed
        got = float(expr.subs({v1: 0.5, v2: 0.0}))
        assert_allclose(got, 1e-13 * np.exp(0.5), rtol=1e-12)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved


## ------------------------------------------------------------------------
## Phase D (hdl.md sec. 9): the deferred table -- what was implemented from
## it, and what is refused on purpose.


def test_discontinuity_is_parsed_and_ignored():
    """`$discontinuity` is advisory. Accepting and ignoring it costs
    nothing and lets a model written for another simulator compile
    unchanged (41 call sites in the vacask library); rejecting it would
    fail those models over a hint they can do without."""
    import sympy
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       discontinuity)
    from pycircuit.utilities.param import Parameter

    class D(Behavioural):
        instparams = [Parameter(name='g', desc='g', unit='S', default=1e-3)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            return Contribution(b.I, g * b.V + discontinuity(1))  # noqa

    pycircuit.circuit.circuit.default_toolkit = numeric
    el = D('p', 'n', g=2e-3)
    el.update_iparv()
    assert_allclose(float(np.asarray(el.i(np.array([1.0, 0.0])), float)[0]),
                    2e-3, rtol=1e-12)


def _laplace_amp(cls, f, tend_cycles=30, pts=400):
    import warnings
    from pycircuit.circuit.elements import VSin
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    c['vs'] = VSin(na, gnd, va=1.0, freq=f)
    c['F'] = cls(na, gnd, nb, gnd)
    c['Rl'] = R(nb, gnd, r=1e6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(c, toolkit=numeric, uic=True).solve(
            tend=tend_cycles / f, timestep=1.0 / (f * pts))
    y = np.asarray(res.v('b').y, float)
    t = np.asarray(res.v('b').x[0], float)
    return float(np.max(np.abs(y[t > (tend_cycles - 10) / f])))


TAU = 1e-4


def _laplace_filter(num, den):
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       laplace_nd)

    class _F(Behavioural):
        instparams = []

        @staticmethod
        def analog(inp, inn, outp, outn):
            bi, bo = Branch(inp, inn), Branch(outp, outn)
            return Contribution(bo.V, laplace_nd(bi.V, num, den))
    return _F


def test_laplace_nd_matches_the_analytic_response():
    """A rational transfer function realised as STATE EQUATIONS -- N
    unknowns in controllable canonical form for an order-N denominator,
    so the simulator integrates them with its own method and the Jacobian
    is exact through them.  Checked against |H(jw)| across three decades.
    """
    lp1 = _laplace_filter([1], [1, TAU])
    fc = 1.0 / (2 * np.pi * TAU)
    for f in (0.1 * fc, fc, 10 * fc):
        got = _laplace_amp(lp1, f)
        want = 1.0 / np.sqrt(1.0 + (2 * np.pi * f * TAU) ** 2)
        assert_allclose(got, want, rtol=0.03), (f, got, want)

    ## Second order: one state per order, so this one carries two.
    lp2 = _laplace_filter([1], [1, 2 * 0.7 * TAU, TAU ** 2])
    el = lp2('a', 'b', 'c', 'd')
    el.update_iparv()
    assert sum('_state' in str(nd) for nd in el.nodes) == 2
    w = 2 * np.pi * fc
    want = abs(1.0 / (1 + 2 * 0.7 * TAU * 1j * w + (TAU * 1j * w) ** 2))
    assert_allclose(_laplace_amp(lp2, fc), want, rtol=0.03)


def test_laplace_zp_expands_to_the_same_filter():
    """Zero/pole form is converted to coefficients, so the realisation is
    identical: a single real pole at -1/tau is 1/(1 + s*tau)."""
    from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                       laplace_zp)

    class ZP(Behavioural):
        instparams = []

        @staticmethod
        def analog(inp, inn, outp, outn):
            bi, bo = Branch(inp, inn), Branch(outp, outn)
            ## one real pole at -1/TAU, no zeros
            return Contribution(bo.V, laplace_zp(bi.V, [], [-1 / TAU, 0.0]))

    fc = 1.0 / (2 * np.pi * TAU)
    got = _laplace_amp(ZP, fc)
    assert_allclose(got, 1.0 / np.sqrt(2.0), rtol=0.03)


def test_improper_laplace_is_refused():
    from pycircuit.circuit.hdl import laplace_nd
    with pytest.raises(NotImplementedError, match='proper'):
        laplace_nd(1.0, [1, 1, 1], [1, 1])


def test_switch_branch_is_refused_clearly():
    """V and I contributed to ONE branch is Verilog-A's switch branch:
    the two are mutually exclusive under a condition, and selecting
    between them per operating point changes the matrix structure every
    iteration.  This compiler used to accept both unconditionally and
    silently build a voltage source with a conductance in parallel --
    defined, but not what the model says."""
    from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
    from pycircuit.utilities.param import Parameter

    pycircuit.circuit.circuit.default_toolkit = numeric
    with pytest.raises(NotImplementedError, match='switch branch'):
        class Sw(Behavioural):
            instparams = [Parameter(name='g', desc='g', unit='S',
                                    default=1e-3)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return (Contribution(b.I, g * b.V),      # noqa: F821
                        Contribution(b.V, 0.5))
