from pycircuit.circuit import *
from pycircuit.circuit.shooting import *
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest

class myC(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(np.zeros(2))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    >>> c.C(np.zeros(2))
    array([[  1.0000e-12,  -1.0000e-12],
           [ -1.0000e-12,   1.0000e-12]])

    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c0', desc='Capacitance', 
                            unit='F', default=1e-12),
                  Parameter(name='c1', desc='Nonlinear capacitance', 
                            unit='F', default=0.5e-12),
                  Parameter(name='v0', desc='Voltage for nominal capacitance', 
                            unit='V', default=1),
                  Parameter(name='v1', desc='Slope voltage ...?', 
                            unit='V', default=1)]

    def C(self, x, epar=defaultepar): 
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        c = c0+c1*self.toolkit.tanh((v-v0)/v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        q = c0*v+c1*v1*self.toolkit.ln(self.toolkit.cosh((v-v0)/v1))
        return self.toolkit.array([q, -q])

def test_shooting():
    """PSS of a linear RC network must match the AC steady state.

    The shooting solver uses backward Euler, which is first order, so the
    error in the periodic steady state scales like O(1/N) in the number of
    timesteps per period.  N is chosen so the discretisation error is a few
    per-mille; the assertions below allow a comfortable margin on top of that
    rather than demanding near-exact equality (which no first-order method
    reaches at a practical N).
    """
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()

    N = 500
    period = 1e-3

    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=1/period, phase=20)
    cir['R'] = R(1,2, r=1e4)
    cir['C'] = C(2,gnd, c=1e-8)

    resac = AC(cir).solve(1/period)

    pss = PSS(cir)

    res = pss.solve(period=period, timestep = period/N)

    v2ac = resac.v(2,gnd)
    v2pss = res['tpss'].v(2,gnd)

    t,dt = numeric.linspace(0,period,num=N,endpoint=True,
                            retstep=True)

    v2ref = numeric.imag(v2ac * numeric.exp(2j*numeric.pi*1/period*t))

    w2ref = Waveform(t,v2ref,ylabel='reference', yunit='V',
                     xunits=('s',), xlabels=('vref(2,gnd!)',))

    ## Check amplitude of the fundamental against the AC result
    v2rms_ac = np.abs(v2ac) / np.sqrt(2)
    v2rms_pss = np.abs(res['fpss'].v(2,gnd)).value(1/period)
    relerr = np.abs(v2rms_pss - v2rms_ac) / v2rms_ac
    assert relerr < 1e-2, 'amplitude rel. error=%g too high' % relerr

    ## Check error of the full waveform against the AC-reconstructed reference
    rmserror = np.sqrt(average((v2pss-w2ref)**2))
    assert rmserror < 1e-2, 'rmserror=%f too high' % rmserror

 
def test_PSS_nonlinear_C():
    """PSS with a tanh-nonlinear capacitor: the result must be PERIODIC.

    STAGE 10.2.  This test asserted NOTHING -- it called `pss.solve(...)` and
    checked no property of the result, so it passed whatever PSS returned,
    including nothing sensible.  It is the same class of defect as
    `test_sparse_toolkit` (passed while never exercising the sparse path) and
    `test_PAC` (skipped): a test that reads as coverage and is not.

    Periodicity is the right assertion because it is exactly what PSS claims to
    deliver -- x(0) = x(T) -- and it needs no external reference to check
    against, so it cannot be fitted. Measured at 2.2e-05 relative; asserted at
    1e-3, which is ~45x margin.

    The second assertion matters as much: without it the test would also pass on
    a degenerate solution that never leaves the capacitor's linear region, where
    `myC` is just a capacitor and the nonlinearity under test is not exercised.
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    pss = PSS(c)
    res = pss.solve(period=1/50e3,timestep=1/50e3/20)

    X = np.asarray(res['tpss'].x, dtype=float)
    assert np.isfinite(X).all(), 'PSS returned non-finite values'

    scale = np.abs(X).max()
    periodicity = np.abs(X[:, 0] - X[:, -1]).max() / scale
    assert periodicity < 1e-3, \
        'the steady state is not periodic: |x(0)-x(T)|/|x| = %.3e' % periodicity

    ## The nonlinearity must actually be traversed: myC's capacitance is
    ## c0 + c1*tanh((v - v0)/v1) with v0 = 1 V, so a solution confined near v0
    ## would exercise a plain capacitor and prove nothing about this circuit.
    v = X[c.get_node_index(2)]
    assert v.max() - v.min() > 2.0, \
        'v(C) spans only %.3f V -- the tanh knee is not being crossed' % (v.max() - v.min())


def test_PAC_is_withdrawn():
    """Stage 11: PAC says it is unimplemented instead of allocating 420 GiB.

    It was `@unittest.skip("Skip failing test")` -- advertised in the analysis
    inventory, never validated, and forming the whole (N*M)x(N*M) operator densely.
    An analysis that announces its absence is strictly better than one that fails
    somewhere deep in an allocation.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    pss = PSS(cir)
    res = pss.solve(period=1e-6, timestep=1e-6/10)
    with pytest.raises(NotImplementedError, match='withdrawn as unimplemented'):
        PAC(cir).solve(pss, np.array([1e6]))


@unittest.skip("Superseded by test_PAC_is_withdrawn; kept for the rewrite")
def test_PAC():
    circuit.default_toolkit = circuit.numeric
    N = 10
    fc = 1e6

    cir = SubCircuit()
    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=fc, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2,gnd)
    cir['C'] = C(2,gnd, c=1e-12)
    
    pss = PSS(cir)

    res = pss.solve(period=1/fc, timestep = 1/(fc*N))

    pac = PAC(cir)
    res = pac.solve(pss, freqs = fc + np.array([1e3, 2e3, 4e3]))
    
    assert False, "Test should compare with spectre simulation"


## ---------------------------------------------------------------------------
## STAGE 11 -- PSS: `method` now selects something, and the inverse is a solve.
## ---------------------------------------------------------------------------

def _series_rlc(Lv=1e-3, Cv=1e-9, Rv=50.0, va=1.0):
    """Series RLC driven AT resonance, where |v(C)| = Q * va analytically."""
    import numpy as _np
    circuit.default_toolkit = circuit.numeric
    f0 = 1.0 / (2 * _np.pi * _np.sqrt(Lv * Cv))
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, va=va, freq=f0)
    c['R'] = R(1, 2, r=Rv)
    c['L'] = L(2, 3, L=Lv)
    c['C'] = C(3, gnd, c=Cv)
    return c, f0, (1.0 / Rv) * _np.sqrt(Lv / Cv)


def _pss_peak(method, steps=20):
    import warnings
    cir, f0, Q = _series_rlc()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = PSS(cir, method=method).solve(period=1.0 / f0,
                                            timestep=1.0 / (f0 * steps))
    v = np.asarray(res['tpss'].v(3, gnd), dtype=float)
    return 0.5 * (v.max() - v.min()), Q


def test_pss_method_parameter_is_actually_read():
    """It was declared with default='euler' and never read anywhere in the file.

    A knob that advertises a choice it does not make is the "thin advertised
    feature" 0.1c warns about -- the same defect `lte_formula` was removed for in
    9(f).  Here the choice is worth having, so it is wired rather than deleted.
    """
    euler_peak, _Q = _pss_peak('euler')
    trap_peak, _Q = _pss_peak('trap')
    assert abs(trap_peak - euler_peak) > 0.1 * euler_peak, \
        'method= changes nothing: euler %.4f, trap %.4f' % (euler_peak, trap_peak)

    with pytest.raises(ValueError, match='method must be'):
        _pss_peak('bogus')


def test_backward_euler_damps_the_limit_cycle_and_trapezoidal_does_not():
    """The defect 0.1c names, measured against an ANALYTIC answer.

    Backward Euler's numerical damping attenuates exactly the limit cycle PSS
    exists to find.  On a Q=20 resonator driven at resonance, where the analytic
    peak is Q*va = 20 V:

        steps/period   euler          trapezoidal
                  20   2.63 V (13%)   19.32 V (97%)
                 200  12.20 V (61%)   19.23 V (96%)

    Euler is not merely less accurate, it is WRONG BY A FACTOR and gets worse as
    the step coarsens -- silently.
    """
    trap_peak, Q = _pss_peak('trap', steps=20)
    euler_peak, _ = _pss_peak('euler', steps=20)
    assert trap_peak > 0.9 * Q, \
        'trapezoidal recovers only %.1f%% of the analytic amplitude' % (100 * trap_peak / Q)
    assert euler_peak < 0.5 * Q, \
        'the euler damping this test documents has gone: %.1f%%' % (100 * euler_peak / Q)


def test_pss_uses_a_solve_not_an_explicit_inverse():
    """`inv(Jf) @ C @ Jshoot` formed a dense inverse per timestep per iteration.

    Asserted structurally: the quantity wanted is the solution of
    `Jf X = C @ Jshoot`, and at N=137/M=1000 with 20 shooting iterations the old
    form was 20,000 dense inversions.  A timing test would be a flake; the source
    check says exactly what changed.
    """
    import inspect
    from pycircuit.circuit import shooting
    ## Both, because the accumulation moved out of `solve` into `_traverse`
    ## when the autonomous system began sharing the period map -- a source
    ## check has to follow the code it is about.
    src = (inspect.getsource(shooting.PSS.solve)
           + inspect.getsource(shooting.PSS._traverse))
    assert 'linalg.inv' not in src, 'the explicit inverse is back'
    assert 'linearsolver' in src


def test_pss_still_matches_the_ac_reference_with_a_fine_step():
    """Both methods must converge to the same, correct answer.

    At a coarse step neither is reliable and Euler's closeness is coincidence --
    measured, it is 0.9886 of the AC answer at dt = RC but 1.3283 at dt = RC/4.
    With a fine enough step both land on 1.0000, which is what makes the
    resonator comparison above a statement about damping rather than about luck.
    """
    import warnings
    from pycircuit.circuit.analysis_ss import AC
    circuit.default_toolkit = circuit.numeric

    def build():
        c = SubCircuit()
        c['VSin'] = VSin(gnd, 1, va=10, freq=50e3, vac=10)
        c['R1'] = R(1, 2, r=1e6)
        c['C'] = C(2, gnd, c=1e-12)
        c['L'] = L(2, gnd, L=1e-3)
        return c

    f = 50e3
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ac = AC(build()).solve(freqs=np.array([f]))
    ref = abs(complex(np.asarray(ac.v(2, gnd)).ravel()[0]))

    for method in ('euler', 'trap'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = PSS(build(), method=method).solve(period=1 / f,
                                                    timestep=1 / f / 1280)
        v = np.asarray(res['tpss'].v(2, gnd), dtype=float)
        amp = 0.5 * (v.max() - v.min())
        assert amp == pytest.approx(ref, rel=0.02), \
            '%s gives %.6f against the AC reference %.6f' % (method, amp, ref)


# ---------------------------------------------------------------------------
# Phase 1: the shooting Newton was not a Newton
# ---------------------------------------------------------------------------

def _q20_rlc(f0=1e3, Q=20.0):
    """A resonator whose per-period decay is exp(-pi/Q) = 0.8546.

    That number is the whole diagnostic: successive substitution converges at
    exactly the circuit's own decay rate, so observing 0.855 per iteration is
    how the missing Jacobian was found.
    """
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def _shooting_trace(method, reltol=1e-4, maxiterations=30):
    """Run PSS, returning (residual per outer iteration, non-converged?)."""
    import warnings as _w
    import pycircuit.circuit.analysis as _an
    trace, orig = [], _an.fsolve

    def spy(f, x0, *a, **kw):
        if f.__qualname__ != 'PSS.solve.<locals>.func':
            return orig(f, x0, *a, **kw)

        def logged(x, *aa):
            F, J = f(x, *aa)
            trace.append((float(np.max(np.abs(F))),
                          float(np.max(np.abs(np.eye(len(x)) - np.asarray(J))))))
            return F, J
        logged.__qualname__ = f.__qualname__
        return orig(logged, x0, *a, **kw)

    circuit.default_toolkit = circuit.numeric
    _an.fsolve = spy
    try:
        with _w.catch_warnings(record=True) as caught:
            _w.simplefilter('always')
            res = PSS(_q20_rlc(), method=method, reltol=reltol).solve(
                period=1e-3, timestep=1e-5, maxiterations=maxiterations)
        nonconv = any('did not converge' in str(c.message) for c in caught)
    finally:
        _an.fsolve = orig
    return trace, nonconv, res


def test_the_shooting_jacobian_is_not_identically_zero():
    """The regression on the defect itself.

    `Jshoot = solve(Jf, C @ Jshoot)` used the RAW capacitance matrix where
    backward Euler's per-step sensitivity is `Jf^-1 C(x_{n-1})/h`.  C is
    singular, so the accumulated product collapsed to EXACTLY zero and the
    Jacobian handed to fsolve was the identity -- making the "shooting
    Newton" plain successive substitution, on every circuit, silently.
    """
    trace, _nc, _res = _shooting_trace('euler')
    jmax = max(j for _f, j in trace)
    assert jmax > 1e-3, \
        'the monodromy is ~zero (max %.3e): the shooting Newton has ' \
        'degenerated to successive substitution' % jmax


def test_euler_shooting_converges_like_a_newton():
    """Few iterations, and tightening the tolerance costs few more.

    Successive substitution on this circuit contracts by 0.8546 per
    iteration, so reaching 1e-9 would need ~130.  A Newton reaches it in a
    handful; measured at landing, 5 iterations at reltol 1e-4 and 10 at
    1e-9.
    """
    loose, nc_loose, _r = _shooting_trace('euler', reltol=1e-4)
    tight, nc_tight, _r2 = _shooting_trace('euler', reltol=1e-9)

    assert not nc_loose and not nc_tight
    assert len(loose) <= 8, 'euler took %d shooting iterations' % len(loose)
    assert len(tight) <= 15, 'euler took %d at reltol 1e-9' % len(tight)
    ## five extra decades for a handful of iterations is the Newton signature
    assert tight[-1][0] < 1e-9, 'final residual %.3e' % tight[-1][0]


def test_non_convergence_is_reported():
    """It used to be silent, which is why the missing Jacobian survived.

    `fsolve` builds the "No convergence" message and discards it whenever
    `full_output=False` -- how PSS called it -- so a shooting solve that
    never converged returned a plausible waveform with no diagnostic.

    This test used to pin the report on TRAPEZOIDAL, which did not converge
    because its monodromy was structurally incomplete.  Phase 3 gave it the
    `(x, iq)` monodromy and it converges in 6 iterations, so that case
    is gone -- as it should be.  The property being protected was never
    "trap fails"; it is "a capped solve says so".  An iteration cap is the
    honest way to produce one, because it cannot expire when a method is
    repaired.
    """
    _trace, nonconv, _res = _shooting_trace('trap', maxiterations=2)
    assert nonconv, 'a non-converged shooting solve returned silently'

    ## and the same solve, uncapped, must NOT warn -- otherwise the
    ## assertion above would pass on a warning that fires unconditionally
    _t2, still, _r2 = _shooting_trace('trap', maxiterations=30)
    assert not still, 'the warning fires even on a converged solve'


def test_pss_tolerance_parameters_reach_the_shooting_solve():
    """`reltol` was a dead knob for the outer Newton.

    `solve()` passed neither tolerance to `fsolve`, so the shooting solve ran
    on library defaults while the inner solves used `par.reltol` -- the two
    were unrelated, which is exactly what the inner/outer ordering rule
    exists to prevent.
    """
    loose, _a, _b = _shooting_trace('euler', reltol=1e-4)
    tight, _c, _d = _shooting_trace('euler', reltol=1e-9)
    assert tight[-1][0] < loose[-1][0] / 100.0, \
        'tightening reltol did not tighten the shooting residual ' \
        '(%.3e vs %.3e)' % (tight[-1][0], loose[-1][0])


def test_pss_tolerances_mean_what_they_mean_in_transient():
    """`reltol`/`iabstol`/`vabstol` are advertised with the same words on
    `Transient`, `JAXTransient` and `PSS`, so they must mean the same thing.

    They did not.  PSS's per-timestep Newton passed neither absolute
    tolerance to `fsolve`, so it ran on library scalar defaults and the two
    Parameters this class documents did nothing to it at all -- while on
    `Transient` they set the per-unknown floors of exactly the same solve.

    Asymmetric values here on purpose: `iabstol` and `vabstol` share the
    default 1e-12, so a swapped FLAVOUR is invisible with the defaults.
    """
    from pycircuit.circuit.analysis import newton_tolerance_vectors
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric

    cir = _q20_rlc()
    n_nodes, n_branches = len(cir.nodes), len(cir.branches)
    IAB, VAB = 3e-9, 7e-6          # distinct, and distinct from the defaults

    abstol, xtol = newton_tolerance_vectors(n_nodes, n_branches, IAB, VAB,
                                            circuit.numeric)
    ## residual flavour: currents on node rows, volts on branch rows
    assert np.allclose(np.asarray(abstol)[:n_nodes], IAB)
    assert np.allclose(np.asarray(abstol)[n_nodes:], VAB)
    ## increment flavour: transposed
    assert np.allclose(np.asarray(xtol)[:n_nodes], VAB)
    assert np.allclose(np.asarray(xtol)[n_nodes:], IAB)

    ## and the transient reaches that same definition
    tr = Transient(cir, toolkit=circuit.numeric, iabstol=IAB, vabstol=VAB)
    assert np.allclose(np.asarray(tr._newton_abstol_vector()),
                       np.asarray(abstol))
    assert np.allclose(np.asarray(tr._newton_xtol_vector()),
                       np.asarray(xtol))


def test_pss_absolute_tolerances_reach_the_inner_solve():
    """Anti-dead-knob: loosening the absolute floors must change the work.

    With the floors set absurdly wide the per-timestep Newton should accept
    almost immediately, so the run differs from one at the defaults.  A knob
    that is accepted and ignored is this codebase's most-paid-for defect
    class, and these two were exactly that here.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric

    def run(**kw):
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            res = PSS(_q20_rlc(), method='euler', **kw).solve(
                period=1e-3, timestep=1e-5, maxiterations=8)
        return np.asarray(res['tpss'].v('c'), dtype=float).ravel()

    tight = run()
    loose = run(iabstol=1e-2, vabstol=1e-2)
    assert not np.allclose(tight, loose), \
        'iabstol/vabstol made no difference to the inner solve'


def test_steadyratio_relates_the_shooting_criterion_to_reltol():
    """`reltol` is the TRANSIENT tolerance, in every analysis; `steadyratio`
    expresses the shooting criterion against it.

    Default 1 means the shooting solve is held to the same relative
    tolerance as the transient -- not tighter, which it could not achieve
    (the period map is only known to the accuracy of the per-timestep
    solves), and not looser by some hidden constant either.  Raising it buys
    fewer shooting iterations for a looser periodic steady state.

    Measured at landing on this resonator at reltol 1e-9: steadyratio 1 took
    10 iterations to |F| = 3.1e-11, and 100 took 8 to 6.4e-9.
    """
    tight, _nc1, _r1 = _shooting_trace('euler', reltol=1e-9,
                                       maxiterations=40)
    circuit.default_toolkit = circuit.numeric
    import warnings as _w
    import pycircuit.circuit.analysis as _an
    trace, orig = [], _an.fsolve

    def spy(f, x0, *a, **kw):
        if f.__qualname__ != 'PSS.solve.<locals>.func':
            return orig(f, x0, *a, **kw)

        def logged(x, *aa):
            F, J = f(x, *aa)
            trace.append(float(np.max(np.abs(F))))
            return F, J
        logged.__qualname__ = f.__qualname__
        return orig(logged, x0, *a, **kw)

    _an.fsolve = spy
    try:
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            PSS(_q20_rlc(), method='euler', reltol=1e-9,
                steadyratio=100.0).solve(period=1e-3, timestep=1e-5,
                                         maxiterations=40)
    finally:
        _an.fsolve = orig

    ## looser criterion => stops earlier, at a larger residual
    assert len(trace) < len(tight), \
        'steadyratio did not relax the shooting solve (%d vs %d iterations)' \
        % (len(trace), len(tight))
    assert trace[-1] > tight[-1][0]


def test_steadyratio_below_one_is_refused():
    """A shooting tolerance tighter than the transient's asks the outer
    residual to resolve the inner solves' own noise.  Refused, not silently
    accepted -- the period map is simply not known that well."""
    circuit.default_toolkit = circuit.numeric
    with pytest.raises(ValueError, match='steadyratio must be >= 1'):
        PSS(_q20_rlc(), method='euler', steadyratio=0.01).solve(
            period=1e-3, timestep=1e-5, maxiterations=4)


# ---------------------------------------------------------------------------
# Phase 2: PSS drives Transient
# ---------------------------------------------------------------------------

def test_pss_finds_the_conducting_solution_of_a_rectifier():
    """The payoff for driving `Transient` instead of a private step.

    PSS carried its own transcription of one integrator step, with no
    limiting.  On a rectifier that Newton never gets the diode to turn on,
    and the non-conducting solution IS periodic -- so the shooting solve
    converged to it and reported success.  Measured before the change: a
    40 V drive returned v(c) spanning +-2.4e-07 V, i.e. reverse leakage,
    with no diagnostic of any kind.  A silently wrong answer, not a
    failure, which is the worse of the two.

    Validated against the circuit integrated to steady state rather than
    against arithmetic: 40 periods of transient on the same grid, last
    period compared point for point.  Measured at landing, 5.8e-04 -- 0.01%
    of the 3.94 V ripple.
    """
    import warnings
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import EulerIntegrator
    circuit.default_toolkit = circuit.numeric

    per, n = 1e-3, 200

    def rect():
        c = SubCircuit()
        c['vs'] = VSin('a', gnd, va=10.0, freq=1 / per)
        c['R'] = R('a', 'b', r=1e3)
        c['D'] = Diode('b', 'c')
        c['RL'] = R('c', gnd, r=1e4)
        c['CL'] = C('c', gnd, c=1e-7)
        return c

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = PSS(rect(), method='euler', reltol=1e-6).solve(
            period=per, timestep=per / n, maxiterations=20)
    t_p = np.asarray(res['tpss'].sweep_values, dtype=float)
    v_p = np.asarray(res['tpss'].v('c'), dtype=float).ravel()

    ## the diode must actually conduct -- the defect this replaces returned
    ## a waveform six orders smaller than this bound
    assert v_p.max() > 1.0, \
        'the rectifier never conducted: v(c) peaks at %.3e' % v_p.max()

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rt = Transient(rect(), toolkit=circuit.numeric,
                       integrator=EulerIntegrator(), reltol=1e-6).solve(
            tend=40 * per, timestep=per / n, fixed_timestep=True)
    t_t = np.asarray(rt.v('c').x, dtype=float).ravel()
    v_t = np.asarray(rt.v('c').y, dtype=float).ravel()
    last = t_t >= 39 * per
    t_l, v_l = t_t[last] - 39 * per, v_t[last]

    dev = float(np.max(np.abs(v_p - np.interp(t_p, t_l, v_l))))
    ripple = float(v_l.max() - v_l.min())
    assert dev < 0.01 * ripple, \
        'PSS differs from the settled transient by %.3e (%.2f%% of ripple)' \
        % (dev, 100 * dev / ripple)


def test_pss_uses_the_transient_integrator_not_a_private_copy():
    """One integrator definition, reached through the real class.

    The private step is gone; `method` now selects an `Integrator` object
    that `Transient.get_diff` drives, which is why `_effective_method`
    reports what actually ran -- including the order drop the integrator
    applies on the first step of each period.
    """
    circuit.default_toolkit = circuit.numeric
    import warnings
    pss = PSS(_q20_rlc(), method='trap')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-5, maxiterations=3)
    tr = pss._transient()
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    assert isinstance(tr.base_integrator, TrapezoidalIntegrator)
    assert tr._effective_method in ('TrapezoidalIntegrator',
                                    'EulerIntegrator')


# ---------------------------------------------------------------------------
# Phase 3: the (x, iq) monodromy
# ---------------------------------------------------------------------------

def test_trapezoidal_shooting_converges_with_the_x_iq_monodromy():
    """Trapezoidal's period map carries `iq`, so the monodromy must too.

    With an x-only monodromy trap did not converge at all -- worse, applying
    the Euler form to it converged SLOWER than no Jacobian (0.90 against
    0.855 per iteration), because the Jacobian was wrong rather than absent.
    Differentiating the two recursions together costs one extra matrix
    product and makes it a real Newton: measured 6 iterations at reltol 1e-4
    and 13 at 1e-9, residual 2.9e-11.
    """
    loose, nc_loose, _r = _shooting_trace('trap', reltol=1e-4,
                                          maxiterations=40)
    tight, nc_tight, _r2 = _shooting_trace('trap', reltol=1e-9,
                                           maxiterations=40)
    assert not nc_loose and not nc_tight
    assert len(loose) <= 10, 'trap took %d shooting iterations' % len(loose)
    assert tight[-1][0] < 1e-9, 'final residual %.3e' % tight[-1][0]


def test_trapezoidal_is_now_right_on_both_axes():
    """The point of the whole repair.

    The two failure modes were orthogonal and neither method escaped both:
    Euler CONVERGED the shooting equation and landed at 8.815 V against a
    20 V analytic peak (a discretisation error), while trapezoidal did not
    converge the shooting equation and landed at 19.848 V.  With the
    `(x, iq)` monodromy, trap converges AND lands at 19.990 V -- 0.05% of
    analytic.

    Euler is asserted unchanged in the same breath, because the unified
    propagation must reduce to the old one when the `iq` row is identically
    zero; if this drifts, the "one formula, two methods" claim is false.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    peaks = {}
    for method in ('euler', 'trap'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = PSS(_q20_rlc(), method=method, reltol=1e-6).solve(
                period=1e-3, timestep=1e-5, maxiterations=40)
        peaks[method] = float(np.max(np.abs(
            np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## trapezoidal does not damp the limit cycle: within 0.5% of 20 V
    assert abs(peaks['trap'] - 20.0) < 0.1, \
        'trap peak %.4f V, analytic 20 V' % peaks['trap']
    ## and Euler still damps it, by its own documented amount -- this is
    ## the level-2 error the shooting solve cannot see
    assert 8.5 < peaks['euler'] < 9.2, \
        'euler peak moved to %.4f V' % peaks['euler']


def test_gear2_shooting_converges_and_damps_between_euler_and_trap():
    """BDF-2 reaches back TWO steps, and the propagation follows it.

    Every method here writes its companion as
    `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`, and the integrator now states
    those coefficients rather than `shooting.py` transcribing them -- this
    file has recorded three times what transcribing an integration constant
    costs.  Euler is `b=0` with one past term, trapezoidal `b=-1` with one,
    Gear-2 `b=0` with two; one recursion serves all three.

    The physics is the check: on a Q=20 resonator against a 20 V analytic
    peak, numerical damping should order euler >> gear2 > trap.  Measured
    at landing -- euler 8.815 V, gear 19.766 V, trap 19.990 V.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    peaks, iters = {}, {}
    for method in ('euler', 'gear', 'trap'):
        trace, nonconv, res = _shooting_trace(method, reltol=1e-9,
                                              maxiterations=40)
        assert not nonconv, '%s did not converge' % method
        iters[method] = len(trace)
        peaks[method] = float(np.max(np.abs(
            np.asarray(res['tpss'].v('c'), dtype=float).ravel())))

    ## all three are real Newtons: a handful of iterations to 1e-9
    for m, k in iters.items():
        assert k <= 20, '%s took %d shooting iterations' % (m, k)

    ## and the damping orders as the methods do
    assert peaks['euler'] < peaks['gear'] < peaks['trap'], \
        'damping does not order euler < gear < trap: %r' % peaks
    assert abs(peaks['gear'] - 20.0) < 1.0, \
        'gear2 peak %.4f V, analytic 20 V' % peaks['gear']


def test_pss_method_selection_cannot_fall_through_silently():
    """A name that is accepted must select the integrator it names.

    ⚠ Found while adding 'gear': the selection was
    `EulerIntegrator() if method == 'euler' else TrapezoidalIntegrator()`,
    so a newly accepted name ran TRAPEZOIDAL -- and it looked like it
    worked, producing numbers identical to trap's to the last digit.  This
    class has already paid once for a `method` Parameter that selected
    nothing at all.
    """
    from pycircuit.circuit.integrator import (EulerIntegrator,
                                              TrapezoidalIntegrator,
                                              Gear2Integrator)
    circuit.default_toolkit = circuit.numeric
    want = {'euler': EulerIntegrator, 'trap': TrapezoidalIntegrator,
            'trapezoidal': TrapezoidalIntegrator,
            'gear': Gear2Integrator, 'gear2': Gear2Integrator}
    for name, cls in want.items():
        tr = PSS(_q20_rlc(), method=name)._transient()
        assert isinstance(tr.par.integrator, cls), \
            'method=%r selected %s' % (name, type(tr.par.integrator).__name__)

    with pytest.raises(ValueError, match="'euler', 'trap' or 'gear'"):
        PSS(_q20_rlc(), method='bdf3').solve(period=1e-3, timestep=1e-5,
                                             maxiterations=2)


# ---------------------------------------------------------------------------
# Phase 4 (arc 5): a phase circuit is an AUTONOMOUS oscillator
# ---------------------------------------------------------------------------

def _phase_circuit():
    """A quadrature phase accumulator driven by a DC source.

    `IdtmodQuadrature` was built so a phase circuit could be handed to a
    shooting analysis: over one output period its state vector returns
    exactly to itself, which the scalar `Idtmod` phase cannot do.  But the
    only excitation is DC -- the oscillation is self-sustaining -- so the
    circuit is autonomous, and that is the property that decides whether
    fixed-period shooting applies.
    """
    from pycircuit.circuit.elements import VS, IdtmodQuadrature
    c = SubCircuit()
    c.add_node('in'); c.add_node('o'); c.add_node('s')
    c['vin'] = VS('in', gnd, v=1e3)
    c['X'] = IdtmodQuadrature('in', gnd, 'o', gnd, 's', gnd, modulus=1.0,
                              amplitude=1.0, ic=0.0)
    c['Ro'] = R('o', gnd, r=1e6)
    c['Rs'] = R('s', gnd, r=1e6)
    return c


def test_an_autonomous_circuit_is_solved_for_its_own_period():
    """Autonomous shooting: the period is an unknown, not an argument.

    ⚠ This test used to assert the opposite -- that a self-oscillating
    circuit could only be DIAGNOSED, returning the trivial orbit with a
    warning.  That was honest while only the fixed-period system existed,
    and expired when the free-period one landed.  The claim it protects is
    the same underneath: such a circuit must not come back silently wrong.

    Why a fixed period cannot work, which is what the free-period system is
    for.  At the nominal period the discretisation precesses -- measured
    2.1e-3 rad per cycle at 100 steps/period, falling as h^2 -- so the
    period map is a rotation by slightly less than 2*pi whose only fixed
    point is the ORIGIN.  Push the period to where the orbit closes and the
    starting phase goes free instead: |eig(M)| 0.9615 -> 1.000226,
    sigma_min(I-M) 2.3e-02 -> 1.6e-04.  Neither has a solution a fixed
    period could find, so `(x0, T)` are solved together with a phase
    condition removing the rotational freedom.

    Validated against an INDEPENDENT measurement: integrating one nominal
    period and reading the angle actually turned gives a precession
    implying T = +83.37 ppm at 200 steps/period; the solver finds
    +83.08 ppm.  The error falls x4 per halving of the step, matching the
    h^2 precession that causes it.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def run(n):
        pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            res = pss.solve(period=1e-3, timestep=1e-3 / n,
                            maxiterations=30)
        nonconv = any('did not converge' in str(c.message) for c in caught)
        return pss, res, nonconv

    pss, res, nonconv = run(200)
    assert pss.autonomous is True
    assert not nonconv, 'the free-period system did not converge'

    ## the orbit closes -- radius is the amplitude, not zero and not drifting
    vo = np.asarray(res['tpss'].v('o'), dtype=float).ravel()
    vs = np.asarray(res['tpss'].v('s'), dtype=float).ravel()
    rad = np.hypot(vo, vs)
    assert abs(rad.max() - 1.0) < 1e-4 and abs(rad.min() - 1.0) < 1e-4, \
        'orbit radius %.7f..%.7f, expected 1' % (rad.min(), rad.max())

    ## and the period it solved for is the precession-corrected one
    ppm = 1e6 * (pss.period - 1e-3) / 1e-3
    assert 70.0 < ppm < 95.0, 'period came out %+.3f ppm, expected ~+83' % ppm

    ## second order: halving the step quarters the correction
    pss2, _r2, _nc2 = run(400)
    ppm2 = 1e6 * (pss2.period - 1e-3) / 1e-3
    assert 3.0 < ppm / ppm2 < 5.0, \
        'period error is not second order in h: %+.3f -> %+.3f' % (ppm, ppm2)


def test_the_free_phase_eigenvalue_appears_at_the_solved_period():
    """The degeneracy is real, and the phase condition is what handles it.

    Solved AT its own period the monodromy has an eigenvalue on the unit
    circle -- that is the rotational freedom, and it is why `I - M` alone is
    singular and the bordered row is not optional.  At the nominal period
    the same circuit reads 0.9615, which is why no spectral threshold can
    detect autonomy (a Q=1000 DRIVEN resonator sits at 0.99686).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
    assert pss.spectral_radius > 0.99, \
        'no unit-circle eigenvalue at the solved period: %.6f' \
        % pss.spectral_radius


@pytest.mark.parametrize('kind', ['rc', 'rectifier'])
def test_driven_circuits_are_not_called_autonomous(kind):
    """Anti-false-positive, and the reason the test is structural.

    ⚠ The spectral radius cannot make this call.  An autonomous orbit gives
    an eigenvalue at 1 only AT its own period -- 0.9615 at the nominal one
    here -- while a merely lightly damped DRIVEN circuit sits near 1 as
    well: a Q=1000 resonator has exp(-pi/Q) = 0.99686.  No threshold
    separates them.  Whether anything depends on `t` does, exactly.
    """
    import warnings
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric

    if kind == 'rc':
        c = SubCircuit()
        c.add_node('a'); c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
        c['R'] = R('a', 'b', r=1e3)
        c['C'] = C('b', gnd, c=1e-7)
    else:
        c = SubCircuit()
        c['vs'] = VSin('a', gnd, va=10.0, freq=1e3)
        c['R'] = R('a', 'b', r=1e3)
        c['D'] = Diode('b', 'c')
        c['RL'] = R('c', gnd, r=1e4)
        c['CL'] = C('c', gnd, c=1e-7)

    pss = PSS(c, method='trap', reltol=1e-6)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=20)
    assert pss.autonomous is False
    assert not any('AUTONOMOUS' in str(x.message) for x in caught)


def _pss_lte(method, timestep=1e-5, reltol=1e-3, **kw):
    """Run the Q=20 resonator and return (peak amplitude, the pss object)."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method=method, reltol=reltol, **kw)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res = pss.solve(period=1e-3, timestep=timestep, maxiterations=40)
    pss._caught = [str(x.message) for x in caught
                   if issubclass(x.category, RuntimeWarning)]
    assert pss.converged, '%s did not converge' % method
    peak = float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel())))
    return peak, pss


def _pss_plain(method, timestep=1e-5, reltol=1e-3, **kw):
    """Force the pre-augmentation formulation, where a seam can exist.

    Gear-2 now solves for its entering history, so its seam is gone by
    construction -- which is the fix, and which leaves the seam machinery
    with nothing to observe unless the old formulation can still be run.
    This is deliberately a test-level override rather than a Parameter: a
    user has no reason to ask for the formulation that measured 1.266e-01 V
    of avoidable error.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method=method, reltol=reltol, **kw)
    pss._solves_history = lambda: False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        res = pss.solve(period=1e-3, timestep=timestep, maxiterations=40)
    pss._caught = [str(x.message) for x in caught
                   if issubclass(x.category, RuntimeWarning)]
    assert pss.converged and not pss.solved_history
    peak = float(np.max(np.abs(
        np.asarray(res['tpss'].v('c'), dtype=float).ravel())))
    return peak, pss


def test_pss_reports_the_truncation_error_neither_newton_can_see():
    """THE THIRD LEVEL, and the only one that ranks the answers.

    Three convergence criteria stand between a PSS run and its answer: the
    inner Newton, the shooting Newton, and the discretisation.  The first
    two are checked while it runs and BOTH SAY YES for all three
    integrators here -- while their amplitudes are 8.815 V, 19.766 V and
    19.990 V against 20 V analytic.  A 56% disagreement between two
    "converged" answers is invisible to every check that asks whether an
    equation was solved, because each integrator solved its own to
    tolerance.

    The report is what sees it, and this test pins the ordering rather than
    the digits: whichever integrator is furthest from the analytic answer
    must carry the largest truncation error.
    """
    peaks, total, peak_lte = {}, {}, {}
    for m in ('euler', 'gear', 'trap'):
        peaks[m], pss = _pss_lte(m)
        total[m], peak_lte[m] = pss.total_lte, pss.max_lte

    ## the physics, unchanged: damping orders euler >> gear2 > trap
    assert peaks['euler'] < peaks['gear'] < peaks['trap'], peaks
    ## and the report orders the same way, both per step and per period
    assert total['euler'] > total['gear'] > total['trap'], total
    assert peak_lte['euler'] > peak_lte['gear'] > peak_lte['trap'], peak_lte


def test_pss_lte_per_step_peak_is_not_enough_for_a_limit_cycle():
    """⚠ THE PER-STEP CRITERION PASSES THE 56%-LOW ANSWER.

    A transient controls its grid on the peak per-step error, and by that
    criterion euler at 100 points/period is IN TOLERANCE at reltol=1e-3
    (0.288).  Its amplitude is 8.815 V against 20 V.  Nothing is wrong with
    the estimate -- it bounds one step, and the 56% is what 99 of them do
    together -- which is exactly why a periodic analysis cannot report only
    that number.  `total_lte` sums the period and reads ~26.

    Deliberately expires: if the estimator ever bounds accumulated error
    directly, this test's premise is gone and it should be rewritten, not
    deleted -- the property to keep is that the report flags the 8.815 V.
    """
    peak, pss = _pss_lte('euler')
    assert abs(peak - 20.0) / 20.0 > 0.4, \
        'euler no longer damps this hard: %.4f V' % peak
    assert pss.max_lte < 1.0, \
        'the per-step peak now flags it; rewrite this test, do not delete it'
    assert pss.total_lte > 1.0, \
        'the period total must flag an answer this far off: %r' % pss.total_lte
    ## and the accurate one must NOT be flagged, or the report is just noise
    peak_t, trap = _pss_lte('trap')
    assert abs(peak_t - 20.0) / 20.0 < 1e-3, peak_t
    assert not trap._caught, \
        'trapezoidal lands at %.5f V and must not be warned about: %r' % (
            peak_t, trap._caught)
    assert pss._caught, 'nothing warned about a 56%-low "converged" answer'
    assert 'accumulated over the period' in pss._caught[0], pss._caught


def test_pss_lte_seam_is_separated_because_it_does_not_follow_the_grid():
    """The cold-start seam is a property of the map, not of the timestep.

    Every shooting iteration re-integrates the period from its own `x0`
    with a fabricated flat history -- that is what keeps phi a function of
    `x0` alone -- so the discrete period map really does open with an
    order-dropped step off a past that never happened.  For the multistep
    methods that seam dwarfs the interior and, unlike the interior, it does
    not fall when the grid is refined.  Reporting one number would let it
    hide the one a smaller timestep can fix.
    """
    _p1, coarse = _pss_plain('gear', timestep=1e-5)
    _p2, fine = _pss_plain('gear', timestep=5e-6)

    ## the interior behaves like a transient's error: refine, it falls
    assert fine.max_lte < 0.5 * coarse.max_lte, \
        'interior LTE did not fall with the grid: %r -> %r' % (
            coarse.max_lte, fine.max_lte)
    assert fine.total_lte < coarse.total_lte

    ## the seam does not -- and it is much larger, so a single max would
    ## have reported only this and called a finer grid the remedy
    assert coarse.max_lte_seam > 50 * coarse.max_lte
    assert fine.max_lte_seam > 0.5 * coarse.max_lte_seam, \
        'the seam now improves with the grid (%r -> %r); if that is a real ' \
        'fix, this test should assert the fix, not the old behaviour' % (
            coarse.max_lte_seam, fine.max_lte_seam)


def test_pss_lte_seam_is_reported_only_where_a_method_can_have_one():
    """⚠ THE SEAM WAS FLAGGED FOR TWO METHODS THAT CANNOT HAVE ONE.

    The first condition was `h_last2 is None` -- the reach of the LTE
    estimator's third divided difference, not of the integrator.  Euler's
    companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and `iq_{n-1}`,
    which the order-dropped opening step supplies consistently.  Neither
    touches the fabricated charge, so neither can pay for it -- measured in
    `benchmarks/pss_seam_cost.py` by comparing PSS's fixed point against the
    limit cycle the same grid and method reach with a real history: euler
    5.1e-12 V, trapezoidal 1.3e-11 V, both zero, while the report was
    calling them 0.286 and 15.1 times tolerance.

    Gear-2 reads `q_{n-2}`, which at that step is the entering unknown, and
    the shooting condition does not constrain that to be the orbit's own
    `x(-dt)`.  Its cost there is 1.266e-01 V against an interior
    contribution of 1.070e-01.

    So the test is mechanistic, not numerical: a seam is reported exactly
    for methods whose companion reaches two charges back.
    """
    for method in ('euler', 'trap'):
        _peak, pss = _pss_lte(method)
        assert pss.max_lte_seam is None, \
            '%s cannot have a seam -- its companion never reads the ' \
            'entering unknown -- but one was reported: %r' % (
                method, pss.max_lte_seam)
        assert pss.max_lte is not None, \
            '%s reported no interior LTE at all' % method

    ## Gear-2 DOES read the entering point, so on the plain formulation it
    ## must still be flagged...
    _peak, plain = _pss_plain('gear')
    assert plain.max_lte_seam is not None and plain.max_lte_seam > 1.0, \
        'gear2 reads the entering stand-in and must report it: %r' \
        % plain.max_lte_seam
    ## ...and must NOT be, once that point is an unknown the solve closed.
    ## A flag that survived its own fix would be the worst of both.
    _peak, aug = _pss_lte('gear')
    assert aug.solved_history and aug.max_lte_seam is None, \
        'the solved history is not a seam: %r' % aug.max_lte_seam


def test_pss_lte_floors_are_the_lte_ones_not_the_newton_ones():
    """`lte_vabstol`/`lte_iabstol` move the report and nothing else.

    They exist as separate parameters for the reason `Transient` records:
    one knob must not move the truncation criterion and the Newton
    criterion together.  Raising them relaxes what the report calls
    resolved; `reltol`/`iabstol`/`vabstol` are untouched, so the same
    solution comes back.
    """
    _peak, tight = _pss_lte('euler')
    peak, loose = _pss_lte('euler', lte_vabstol=1.0, lte_iabstol=1.0)

    assert tight._caught and not loose._caught, \
        'the LTE floors did not silence the report: %r' % (loose._caught,)
    assert loose.total_lte < tight.total_lte
    ## same answer, only the accounting moved
    assert abs(peak - _peak) < 1e-9 * max(1.0, abs(peak))
    assert loose.par.reltol == tight.par.reltol


def test_pss_lte_measurement_does_not_touch_the_solution():
    """A measurement that changes what it measures is not one.

    `step_lte` runs inside the timestep loop, before the history push, and
    reads `_qlast`/`_iqlast`/`_q_at` -- all of which the next step depends
    on.  With the estimator neutralised the waveform must come back bit for
    bit, or the report is participating in the answer.
    """
    import warnings
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric

    def run():
        pss = PSS(_q20_rlc(), method='gear', reltol=1e-3)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=1e-3, timestep=1e-5, maxiterations=40)
        return np.asarray(res['tpss'].x, dtype=float)

    measured = run()
    orig = Transient.step_lte
    try:
        Transient.step_lte = lambda self, *a, **kw: None
        unmeasured = run()
    finally:
        Transient.step_lte = orig
    assert_array_equal(measured, unmeasured)


def test_pss_solves_history_exactly_where_the_companion_needs_it():
    """The formulation follows the integrator's reach, not its name.

    MEASURED, in `benchmarks/pss_seam_cost.py`: a companion reading one
    charge back cannot see the fabricated opening history at all -- euler's
    seam costs 5.1e-12 V and trapezoidal's 1.3e-11 V, both zero -- so
    enlarging their system would double the unknowns to fix nothing.
    Gear-2 reads `q_{n-2}`, which in the plain formulation is the entering
    stand-in, and pays 1.266e-01 V at 100 points per period.
    """
    circuit.default_toolkit = circuit.numeric
    want = {'euler': 1, 'trap': 1, 'trapezoidal': 1, 'gear': 2, 'gear2': 2}
    for name, reach in want.items():
        got = PSS(_q20_rlc(), method=name)._companion_reach()
        assert got == reach, '%s reaches %d charges back, not %d' % (
            name, got, reach)

    for name in ('euler', 'trap'):
        _p, pss = _pss_lte(name)
        assert pss.solved_history is False, \
            '%s has no seam to fix and must not pay for a solved ' \
            'history' % name
    _p, gear = _pss_lte('gear')
    assert gear.solved_history is True


def test_pss_solved_history_removes_the_gear2_seam():
    """The fix, against the number that predicted it.

    The seam was measured by continuing 80 periods past the converged solve
    with no re-seed, which reaches the limit cycle the same grid and method
    produce from a real history: 19.89297 V at 100 points per period,
    19.98524 at 200, 19.99735 at 400, against a plain-formulation PSS that
    lands at 19.76639 / 19.95451 / 19.99008 and an analytic 20 V.

    Making the entering history an unknown must reach the FIRST of those
    numbers -- that is what "the seam is the only difference" means -- so
    this asserts the prediction, not just an improvement.
    """
    circuit.default_toolkit = circuit.numeric
    predicted = {100: 19.89297, 200: 19.98524, 400: 19.99735}
    plain_err = {100: 2.336e-01, 200: 4.549e-02, 400: 9.924e-03}
    for npts, target in predicted.items():
        peak, pss = _pss_lte('gear', timestep=1e-3 / npts, reltol=1e-9)
        assert pss.solved_history and pss.converged
        assert abs(peak - target) < 5e-5, \
            'npts=%d landed at %.5f, not the primed limit cycle %.5f' % (
                npts, peak, target)
        ## and it is an improvement, by the factor the measurement predicted
        assert abs(peak - 20.0) < 0.5 * plain_err[npts], \
            'npts=%d error %.3e is not a clear gain on the plain %.3e' % (
                npts, abs(peak - 20.0), plain_err[npts])
        ## the seam is gone, and the report agrees it is gone
        assert pss.max_lte_seam is None, \
            'the solved history is not a seam: %r' % pss.max_lte_seam


def test_pss_solved_history_jacobian_is_the_exact_one():
    """⚠ THE PLAIN PATH'S JACOBIAN CARRIES THE ASSUMPTION IT IS FIXING.

    It seeds both sensitivity rings with `I`, which says `d x_{-1}/d x_0 =
    I` -- the flat history written into the derivative.  Newton tolerates
    that and converges anyway, which is why it was never visible.  With the
    history solved for, the Jacobian is exact and the same circuit needs a
    handful of residual evaluations instead of a dozen.

    So solving for the history is not a cost: it doubles the unknowns of a
    solve that is not the expensive part, and removes iterations from the
    part that is.
    """
    import pycircuit.circuit.analysis as _an
    circuit.default_toolkit = circuit.numeric

    def evals(force_plain):
        calls = [0]
        orig = _an.fsolve

        def spy(f, x0, *a, **kw):
            def wrapped(*aa):
                calls[0] += 1
                return f(*aa)
            wrapped.__qualname__ = getattr(f, '__qualname__', '')
            return orig(wrapped, x0, *a, **kw)
        pss = PSS(_q20_rlc(), method='gear', reltol=1e-9)
        if force_plain:
            pss._solves_history = lambda: False
        _an.fsolve = spy
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=1e-3, timestep=1e-5, maxiterations=40)
        finally:
            _an.fsolve = orig
        assert pss.converged
        return calls[0]

    plain, aug = evals(True), evals(False)
    assert aug < plain, \
        'the exact Jacobian took %d residual evaluations against the ' \
        'approximate one\'s %d' % (aug, plain)


def test_pss_solved_history_refuses_a_companion_it_cannot_seed():
    """⚠ THE REASON THIS TEST PROTECTED WAS WRONG, AND THE REFUSAL SURVIVES.

    It used to say a `b != 0` companion reads an `iq_{n-1}` that "no charge
    determines".  The DAE determines it exactly -- `iq_{-1} =
    -(i(x_{-1}) + u)`, item 4d -- and seeding it that way was built and
    tried.

    It fails for the derivative running the OTHER way.  A one-step companion
    depends on `x_{-1}` ONLY through `iq_{-1}`, and
    `d(iq_{-1})/d x_{-1} = -G` is singular at every purely reactive node --
    most of a resonator.  Admitting `x_{-1}` as m unknowns leaves the
    2m x 2m system rank-deficient: measured, `LinAlgError: Singular matrix`
    on 25 tests at once.

    So the refusal stands on a sharper reason, and the property under test
    is now that reason: such a method must stay OFF this formulation, and
    the second unknown it would need is `iq_{-1}` itself.
    """
    from pycircuit.circuit.integrator import Gear2Integrator
    circuit.default_toolkit = circuit.numeric

    ## trapezoidal has `b != 0`, so it must NOT be routed here
    assert PSS(_q20_rlc(), method='trap')._solves_history() is False, \
        'a b != 0 companion was admitted to the solved-history formulation; ' \
        'its enlarged system is rank-deficient (dq/dx = -G is singular at ' \
        'reactive nodes), so this must stay on the plain path until the ' \
        'iq_{-1} formulation exists'
    assert PSS(_q20_rlc(), method='gear')._solves_history() is True

    ## and reaching the installer with such a companion is refused, loudly
    pss = PSS(_q20_rlc(), method='gear')
    orig = Gear2Integrator.companion_coefficients
    try:
        Gear2Integrator.companion_coefficients = \
            lambda self, h, hl: (orig(self, h, hl)[0], -1.0)
        with pytest.raises(NotImplementedError, match='rank-deficient'):
            pss.solve(period=1e-3, timestep=1e-5, maxiterations=4)
    finally:
        Gear2Integrator.companion_coefficients = orig


def test_the_composed_autonomous_system_removes_the_seam_too():
    """⚠ THIS TEST EXPIRED AS WRITTEN, AND WAS REWRITTEN AS IT ASKED.

    It used to assert that autonomous Gear-2 KEPT its seam, with the note
    that if the two enlargements were ever composed "the wobble assertion
    below is what must flip -- rewrite it, do not delete it".  They were,
    the same day, and this is that rewrite: the property under protection
    is unchanged -- the seam is visible as an orbit that does not close in
    radius -- only its expected value moved.

    An autonomous circuit under a two-step method needs BOTH enlargements:
    the period because it is not given, the history because the companion
    reads it.  The composed system solves `(x_0, x_{-1}, T)` against both
    states closing plus a phase condition.

    The gate is the free-running measurement, as it was for the driven
    case: continuing past the solve on a real history gives +332.184 ppm at
    200 steps and +82.652 at 400, against a plain formulation's +329.682
    and +82.342.  The composed solve must LAND there, not merely improve.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def solved(method, n):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=1e-3, timestep=1e-3 / n, maxiterations=30)
        assert pss.converged and pss.autonomous
        rad = np.hypot(
            np.asarray(res['tpss'].v('o'), dtype=float).ravel(),
            np.asarray(res['tpss'].v('s'), dtype=float).ravel())
        return pss, float(np.max(rad) - np.min(rad)), \
            1e6 * (pss.period - 1e-3) / 1e-3

    ## a two-step companion needs the history even when the period is free
    gear200, wobble200, ppm200 = solved('gear', 200)
    assert gear200.solved_history is True, \
        'a free period does not remove the need for a readable history'

    ## the orbit closes now -- this is the assertion that flipped
    assert wobble200 < 1e-9, \
        'the seam is back: orbit radius spread %.3e (was 2.095e-04 on the ' \
        'plain formulation, and must now be gone)' % wobble200

    ## and the period is the free-running one, not merely nearer to it
    for n, target in ((200, 332.184), (400, 82.652)):
        _p, _w, ppm = solved('gear', n)
        assert abs(ppm - target) < 0.05, \
            'n=%d solved %+.3f ppm, not the free-running %+.3f' % (
                n, ppm, target)

    ## the free-phase eigenvalue survives the enlargement, so the autonomous
    ## diagnostic still says what it said -- on a 2m x 2m spectrum that now
    ## also carries the two-step method's parasitic roots
    assert gear200.spectral_radius > 0.99, \
        'no unit-circle eigenvalue in the composed monodromy: %.6f' \
        % gear200.spectral_radius

    ## ⚠ AND THE VALUE CAVEAT, PINNED.  On this circuit Gear-2 is still the
    ## worse choice -- its own phase error is ~4x trapezoidal's -- so the
    ## default stays `trap`.  What composing buys is that a two-step method
    ## is CORRECT when it is the right tool, i.e. a stiff oscillator.
    _t, _tw, trap_ppm = solved('trap', 200)
    assert abs(ppm200) > 2.0 * abs(trap_ppm), \
        "gear2's phase error (%+.1f ppm) is no longer dominant against " \
        "trapezoidal's (%+.1f ppm); the recommendation to default to trap " \
        "rested on that, so re-measure it" % (ppm200, trap_ppm)


def test_an_autonomous_period_that_is_a_multiple_is_reported_as_one():
    """⚠ `k*T` SOLVES THE PERIODICITY CONDITION WHENEVER `T` DOES.

    So the free-period system has a solution at every integer multiple and
    converges to whichever the seed is nearest -- measured on this element
    (true period 1.000e-03): seeds of 1e-3, 2e-3 and 3e-3 return
    1.000083e-03, 2.000665e-03 and 3.002245e-03, and ALL report
    `converged`.  Each waveform is a correct periodic solution; the
    FUNDAMENTAL FREQUENCY is wrong by the factor, which is usually the
    thing a PSS user wanted.  Nothing said so until this landed.

    The detector needs no extra solve: a k-fold orbit comes back near
    `x_0` partway through.  It must find the EARLIEST such return -- a
    three-fold orbit passes close at `T/3` and `2T/3`, and reporting the
    nearer of the two named a fundamental that was itself a multiple.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def run(seed, method='trap'):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=seed, timestep=seed / 200, maxiterations=30)
        assert pss.converged
        return pss, [str(c.message) for c in caught
                     if 'MULTIPLE of the fundamental' in str(c.message)]

    ## the fundamental itself must NOT warn -- a report that fires on the
    ## right answer is noise
    one, hits = run(1.0e-3)
    assert not hits, 'warned on the fundamental itself: %r' % hits
    assert one.fundamental_period is None

    ## two and three times it must, and must name the right factor
    for seed, factor in ((2.0e-3, 2.0), (3.0e-3, 3.0)):
        pss, hits = run(seed)
        assert hits, 'no warning at %gx the fundamental' % factor
        assert pss.fundamental_period is not None
        ratio = pss.period / pss.fundamental_period
        assert abs(ratio - factor) < 0.25, \
            'seed %g: called it %.2f times the fundamental, expected %g -- ' \
            'if this drifted to a smaller factor the earliest-recurrence ' \
            'rule has regressed to a nearest-recurrence one' % (
                seed, ratio, factor)

    ## a DRIVEN run is exempt: its period is the caller's, and asking for
    ## two source periods is a legitimate request rather than a mistake
    pss = PSS(_q20_rlc(), method='gear', reltol=1e-6)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        pss.solve(period=2e-3, timestep=2e-3 / 200, maxiterations=40)
    assert pss.autonomous is False
    assert not [c for c in caught if 'MULTIPLE' in str(c.message)], \
        'a driven run was told its own requested period is a multiple'


def test_the_trivial_period_root_is_named_not_returned_bare():
    """⚠ `T = 0` IS A REGULAR ROOT OF EVERY AUTONOMOUS SHOOTING SYSTEM.

    `x0 - phi_T(x0)` vanishes identically at `T = 0`, and the phase
    condition does not exclude it -- it constrains `x0`, not the period.
    So a seed below the fundamental is drawn there, and the run either
    returns a period of ~1e-18 or dies with a singular Jacobian on the way
    down.  Measured on BOTH autonomous elements, so it belongs to the
    formulation and not to a circuit: from a 1e-4 seed against a 1e-3
    fundamental, Gear-2 returned -1.5e-20 and 3.9e-19, and trapezoidal
    raised a bare `LinAlgError` from three seeds of five.

    Neither was a SILENT wrong answer -- the collapse reports
    `converged = False` and the exception is loud -- but neither said
    anything about the cause, and the generic advice attached to
    non-convergence ("raise maxiterations") is actively wrong here: no
    number of iterations reaches a fundamental from below.

    Which of the two failure modes appears is method-dependent, so this
    accepts either and insists only that it be NAMED.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    for method in ('trap', 'gear'):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                pss.solve(period=1e-4, timestep=1e-4 / 200, maxiterations=30)
            except np.linalg.LinAlgError as exc:
                assert 'seed BELOW the fundamental' in str(exc), \
                    '%s raised a bare singular-matrix error with no cause: ' \
                    '%s' % (method, exc)
                continue
        named = [c for c in caught if 'TRIVIAL root' in str(c.message)]
        assert named, \
            '%s returned period %r from a seed below the fundamental ' \
            'without naming the trivial root' % (method, pss.period)

    ## and a sound seed is untouched -- the guard must not fire on the
    ## answer it exists to distinguish from
    for method in ('trap', 'gear'):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
        assert pss.converged
        assert not [c for c in caught if 'TRIVIAL root' in str(c.message)], \
            '%s: the guard fired on a good solve' % method
        assert abs(pss.period - 1e-3) / 1e-3 < 1e-3


def test_the_monodromy_is_the_pair_map_not_a_corner_of_it():
    """⚠ A SUB-BLOCK OF A SENSITIVITY IS NOT A MONODROMY.

    For a two-step method the one-period map acts on the PAIR
    `(x_n, x_{n-1})`, so the monodromy is `2m x 2m`.  The solved-history
    path used to hand back the `d x_{N-1}/d x_0` corner of it instead, and
    the eigenvalues of a corner mean nothing: it reported a spectral radius
    of 1.2796 for this resonator -- ABOVE ONE, i.e. an unstable orbit --
    where the circuit decays by `exp(-pi/Q)` every period.

    The analytic decay is the check, so this cannot drift with a
    reimplementation.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    analytic = float(np.exp(-np.pi / 20.0))

    def rho(method, force_plain):
        pss = PSS(_q20_rlc(), method=method, reltol=1e-9)
        if force_plain:
            pss._solves_history = lambda: False
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=1e-3, timestep=1e-5, maxiterations=40)
        assert pss.converged
        return pss

    for method in ('trap', 'gear'):
        for plain in (True, False):
            pss = rho(method, plain)
            assert abs(pss.spectral_radius - analytic) < 0.01, \
                'method=%r plain=%r reports rho=%.6f against the analytic ' \
                'per-period decay %.6f -- a value above 1 here means a ' \
                'corner of the sensitivity is being read as a monodromy' % (
                    method, plain, pss.spectral_radius, analytic)

    ## and the shape follows the method's reach, mechanically
    m = _q20_rlc().n - 1
    assert np.asarray(rho('gear', False)._monodromy).shape == (2 * m, 2 * m)
    assert np.asarray(rho('trap', False)._monodromy).shape == (m, m)


def test_the_parasitic_roots_stay_far_from_the_physical_ones():
    """The k-step method's spurious roots, and why they are not a problem.

    A k-step method turns an m-dimensional system into a k*m-dimensional
    discrete one, so the monodromy's spectrum carries (k-1)*m PARASITIC
    roots beside the physical multipliers -- controlling them is the whole
    subject of the boundary-value-methods literature cited in the class
    docstring.  Reading `max |eig|` off such a spectrum is only safe while
    the parasitic roots stay small.

    For Gear-2 they do, and the reason is quantitative rather than hopeful:
    its parasitic root is 1/3 per STEP (the roots of `1.5z^2 - 2z + 0.5`
    are 1 and 1/3), so over a period it is `(1/3)^N` -- about 1e-95 at 200
    points.  Measured, the autonomous 16x16 spectrum is the physical unit
    eigenvalue and nothing else above 1e-5.

    ⚠ THIS NO LONGER GUARDS `spectral_radius`, and it is kept for what it
    still says.  Since 2026-09-02 the multipliers ARE separated rather than
    maximised over (`_spectral_report`), so a near-unit parasitic root no
    longer reads as stability -- that case is
    `test_a_parasitic_root_near_the_unit_circle_is_not_read_as_stability`.
    What this test still pins is the QUANTITATIVE claim the class docstring
    makes about Gear-2 specifically: its spurious roots are ~1e-95 over a
    period, so on every method in this tree today the separation changes no
    number and only documents why the old maximum was safe.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    ## BDF-2's own roots, so the claim above is checked and not asserted
    assert np.allclose(sorted(np.roots([1.5, -2.0, 0.5])), [1.0 / 3.0, 1.0])

    pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
    ev = np.sort(np.abs(np.linalg.eigvals(np.asarray(pss._monodromy))))[::-1]

    assert abs(ev[0] - 1.0) < 1e-3, \
        'the physical free-phase eigenvalue is not on the unit circle: %r' \
        % ev[0]
    assert ev[1] < 1e-3, \
        'a second eigenvalue at %.3e is no longer negligible -- with the ' \
        'parasitic roots this close to the physical one, `spectral_radius` ' \
        'is a maximum over a mixed spectrum and needs separating' % ev[1]


def _grid_fracs(kind, n):
    if kind == '2:1':
        f = np.where(np.arange(n) % 2 == 0, 2.0, 1.0)
    elif kind == 'smooth':
        f = 1.0 + 0.8 * np.sin(2 * np.pi * np.arange(n) / n)
    else:
        f = np.ones(n)
    return f / f.sum()


def test_a_callers_grid_is_validated_before_it_is_used():
    """The contract is step FRACTIONS of the period, summing to one.

    Fractions and not absolute times, because an autonomous period is an
    unknown: every step has to scale with `T` or `dh/dT = h/T` -- the
    identity the period column rests on -- stops holding.  A grid that
    silently did not sum to one would shorten or lengthen the period the
    solve believes it integrated, which no residual could detect.
    """
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method='trap')
    for bad, why in ((np.array([0.5, 0.4]), 'sum to 1'),
                     (np.array([0.5, 0.6]), 'sum to 1'),
                     (np.array([-0.5, 1.5]), 'positive'),
                     (np.array([1.0]), 'at least two')):
        with pytest.raises(ValueError, match=why):
            pss._period_grid(1e-3, 200, bad)

    ## and a good one lands where it says
    fr = _grid_fracs('2:1', 100)
    times, hs = pss._period_grid(1e-3, 100, fr)
    assert len(hs) == 100 and len(times) == 101
    assert abs(times[-1] - 1e-3) < 1e-15
    assert np.allclose(np.diff(times), hs)


def test_a_non_uniform_grid_solves_a_driven_circuit():
    """RECORDED SCOPE ITEM 5, the driven half.

    ⚠ The recorded blocker -- "blocked on `Transient` accepting a
    non-uniform grid; `fixed_timestep` is uniform-only" -- was stale.
    `Transient.solve`'s loop is uniform-only and always was, but PSS never
    uses that loop: it drives `solve_timestep` one step at a time, and
    non-uniform steps went through unchanged.

    The analytic per-period decay is the check, because it does not care
    what grid produced it.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    analytic = float(np.exp(-np.pi / 20.0))
    for method in ('trap', 'gear'):
        for kind in ('2:1', 'smooth'):
            pss = PSS(_q20_rlc(), method=method, reltol=1e-9)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = pss.solve(period=1e-3, timestep=1e-5, maxiterations=40,
                                grid=_grid_fracs(kind, 200))
            assert pss.converged, '%s/%s did not converge' % (method, kind)
            peak = float(np.max(np.abs(np.asarray(
                res['tpss'].v('c'), dtype=float).ravel())))
            assert abs(peak - 20.0) < 0.5, \
                '%s/%s peak %.5f against 20 V analytic' % (method, kind, peak)
            assert abs(pss.spectral_radius - analytic) < 0.01, \
                '%s/%s rho %.6f against exp(-pi/Q) %.6f' % (
                    method, kind, pss.spectral_radius, analytic)


def test_a_non_uniform_grid_works_when_the_period_is_an_unknown():
    """The autonomous half, which is where the FRACTIONS matter.

    The grid is rebuilt at the current `T` on every residual evaluation, so
    each step scales with the unknown and `dh/dT = h/T` still holds.  If the
    grid were frozen in absolute time instead, the period column would be
    wrong and the solve would converge to the wrong `T` -- which is what
    this asserts against.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'gear'):
        for kind in ('2:1', 'smooth'):
            pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30,
                          grid=_grid_fracs(kind, 200))
            assert pss.converged, '%s/%s did not converge' % (method, kind)
            ## the true period is 1.000e-03; a deliberately awkward grid
            ## costs accuracy but must not move the answer by a factor
            err = abs(pss.period - 1e-3) / 1e-3
            assert err < 5e-3, \
                '%s/%s solved T=%.9f, %.1f ppm from the true 1e-3 -- a grid ' \
                'frozen in absolute time rather than in fractions would ' \
                'fail exactly here' % (method, kind, pss.period, 1e6 * err)


def test_a_grid_that_opens_coarse_is_subdivided_but_a_benign_one_is_not():
    """RECORDED SCOPE ITEM 5: the opening step is MANUFACTURED.

    `_traverse` builds `x(0)` from the unknown with one order-dropped Euler
    step of `hs[0]`, so a grid taken from an adaptive transient opens
    wherever that transient's window happened to start -- which has nothing
    to do with what a good opening step is.  On van der Pol that step is
    3200x the grid's median.

    ⚠ The guard is what this pins as much as the subdivision.  A grid that
    already works must come back EXACTLY as the caller wrote it, or every
    recorded non-uniform result silently moves onto a different grid.
    """
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_q20_rlc(), method='trap')

    ## benign grids are returned untouched -- '2:1' opens at 2x its finest,
    ## 'smooth' at 5x, and a uniform grid at 1x
    for kind, n in (('2:1', 100), ('smooth', 100), ('uniform', 100)):
        fr = _grid_fracs(kind, n)
        times, hs = pss._period_grid(1e-3, n, fr)
        assert len(hs) == n, '%s grid was resized' % kind
        assert np.allclose(hs, fr * 1e-3), '%s grid was rewritten' % kind

    ## a grid opening far coarser than its finest step gains ONE step, and
    ## opens on that finest step
    fr = np.concatenate(([0.5], np.full(500, 0.001)))
    fr = fr / fr.sum()
    times, hs = pss._period_grid(1e-3, len(fr), fr)
    assert len(hs) == len(fr) + 1
    assert hs[0] == pytest.approx(fr.min() * 1e-3, rel=1e-12)
    assert hs[0] + hs[1] == pytest.approx(fr[0] * 1e-3, rel=1e-12)
    ## the period is preserved -- a subdivision that moved it would change
    ## the interval the solve believes it integrated
    assert times[-1] == pytest.approx(1e-3, rel=1e-12)
    assert np.allclose(np.diff(times), hs)

    ## and it is idempotent: the result opens at its own finest step, so
    ## feeding it back changes nothing
    fr2 = hs / hs.sum()
    _t2, hs2 = pss._period_grid(1e-3, len(fr2), fr2)
    assert len(hs2) == len(hs) and np.allclose(hs2, hs)


def _van_der_pol(mu=100.0):
    """The canonical stiff relaxation oscillator; no sources, so autonomous.

        C dv/dt = -i_L + mu (v - v^3/3),    L di_L/dt = v
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    return c


def test_the_lte_chosen_grid_solves_van_der_pol_through_the_analysis():
    """RECORDED SCOPE ITEM 5's PAYOFF CASE, and it is the whole point of it.

    A transient adapts because it cannot see the future; PSS re-solves the
    SAME interval repeatedly, so it can be handed a grid chosen well once.
    The prize on van der Pol at mu=100 is ~18x fewer points than the
    uniform grid that converges (1106 against 20000).

    ⚠ THIS CASE FAILED THROUGH THE ANALYSIS FOR A REASON THE RECORD NAMED
    WRONG.  It was attributed to the plain path's ~30% Jacobian, on the
    evidence that a prototype with finite differences converged.  Measured
    afterwards: an exact finite-difference Jacobian does NOT fix it, and
    with the opening step subdivided the analytic and finite-difference
    Jacobians agree to six digits.  The blocker was the manufactured
    opening step, which the prototype did not have because its unknown was
    `x_0` itself.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 162.842412                      # measured free-running period

    ## one period of ACCEPTED steps from a settled adaptive transient
    cir = _van_der_pol()
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-7).solve(refnode=gnd, tend=1200.0 + T,
                                                timestep=0.05, x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    xs = np.asarray(res.x, dtype=float)
    j0 = int(np.searchsorted(t, t[-1] - T))
    win_t, win_x = t[j0:], xs[:, j0:]
    fr = np.diff(win_t)
    fr = fr / fr.sum()
    iref = cir.get_node_index(gnd)
    seed = np.concatenate((win_x[:iref, 0], win_x[iref + 1:, 0]))

    ## the pathology this case exists for: the window opens on a coarse step
    assert fr[0] > 100 * np.median(fr), \
        'the LTE grid no longer opens coarse (%.3e against a median of ' \
        '%.3e), so this case no longer tests what it was written for' \
        % (fr[0], np.median(fr))

    pss = PSS(_van_der_pol(), method='trap', reltol=1e-7)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=T, x0=seed, grid=fr, maxiterations=25)

    assert pss.converged, \
        'van der Pol did not solve through PSS.solve(grid=...) on its own ' \
        'LTE-chosen grid -- item 5 has no payoff case without this'
    err = 1e6 * (pss.period - T) / T
    assert abs(err) < 200, \
        'solved T=%.6f, %.1f ppm from the measured %.6f' % (pss.period, err, T)
    ## fewer than 1200 steps, against the 20000 a uniform grid needs
    assert len(pss.times) < 1200


def _rc_ladder(sections):
    """A driven RC ladder whose `m` grows linearly with `sections`."""
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, freq=1e3)
    for k in range(sections):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-9)
    return c


def _varying_c_ladder(sections=6):
    """The ladder with a STATE-DEPENDENT capacitance on its last node.

    ⚠ WITHOUT THIS THE MATVEC TESTS CANNOT FAIL, and that was measured, not
    supposed.  `_step_sensitivity` carries a RING of capacitances -- `Cs[0]`
    and `Cs[1]`, the two steps a Gear-2 companion reaches back -- and on a
    linear RC ladder every `C` along the period is the SAME matrix, so
    corrupting the ring (`Cs = [C_new, C_new]` instead of
    `[C_new, Cs[0]]`) changes nothing and the tests stayed green through
    the mutation.  A `q_func` makes `C` a function of the solution, the
    ring entries genuinely differ (measured: 2.99e-09 against a 1e-09
    linear part), and the same mutation is caught.
    """
    c = _rc_ladder(sections)
    last = 'n%d' % sections
    c['Q'] = BSource(last, gnd, last, gnd,
                     q_func=lambda v: 2e-9 * np.tanh(v / 0.5))
    return c


def test_the_matrix_free_matvec_is_the_dense_monodromy():
    """RECORDED SCOPE ITEM 6: `M v` without ever forming `M`.

    The dense path propagates `2m` columns per step; the matrix-free path
    replays the SAME recursion at width 1 on a stored factorisation.  If
    those two disagree, everything built on the matvec is measuring its own
    bug, so this pins them against each other directly rather than against
    a converged answer that could absorb the difference.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_varying_c_ladder(), method='gear', reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 25, maxiterations=20)
    m = pss.cir.n - 1
    times, hs = pss._period_grid(1e-3, 25, None)
    ## ⚠ THE TWO HISTORY STATES MUST DIFFER, and that was also measured.
    ## Seeded with `x_0 = x_{-1} = 0`, the two opening capacitances
    ## `C(x_0)` and `C(x_{-1})` are the same matrix, so seeding the ring
    ## with the WRONG one is invisible -- a mutation replacing
    ## `_C_at(xm1_in)` with `_C_at(x0_in)` stayed green.  Two genuinely
    ## different entering states are what make the opening testable.
    xa = np.full(m, 0.05)
    xb = np.full(m, -0.03)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _xl, _xp, P_last, P_prev = pss._traverse_solved_history(xa, xb, times, hs)
        M = np.vstack((P_last, P_prev))
        C0, steps, _a, _b = pss._traverse_factored(xa, xb, times, hs)

    rng = np.random.default_rng(0)
    for _ in range(3):
        v = rng.standard_normal(2 * m)
        got = pss._monodromy_matvec(C0, steps, v)
        want = M @ v
        err = np.linalg.norm(got - want) / np.linalg.norm(want)
        assert err < 1e-11, 'matvec differs from the dense monodromy by %.3e' % err


def test_a_matrix_free_solve_agrees_with_the_dense_one():
    """The two paths must answer the same, or the fast one is not an option.

    ⚠ The convergence TEST is not identical between them -- `fsolve` scales
    its residual by `|J| . |x|`, which matrix-free has no way to form -- so
    this asserts on the converged ANSWER and the converged/not verdict,
    which are what a caller sees, and not on the iteration count, which
    measurably differs (matrix-free takes one more at m >= 502).
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    out = []
    for mf in (False, True):
        pss = PSS(_varying_c_ladder(), method='gear', reltol=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=1e-3, timestep=1e-3 / 25, maxiterations=20,
                            matrix_free=mf)
        out.append((pss, np.asarray(res['tpss'].x, dtype=float).ravel()))
    (pd, xd), (pm, xm) = out
    assert pd.converged and pm.converged
    rel = np.max(np.abs(xd - xm)) / np.max(np.abs(xd))
    assert rel < 1e-9, 'matrix-free and dense waveforms differ by %.3e' % rel

    ## ⚠ NOT FORMING THE MONODROMY MEANS NOT REPORTING ONE.  `_monodromy`
    ## survives on the object between solves, so a matrix-free run that left
    ## it alone would hand back the PREVIOUS run's radius as its own.
    assert pm.spectral_radius is None
    assert pd.spectral_radius is not None


def test_matrix_free_covers_every_shooting_system():
    """All four systems take `matrix_free=True`; none falls back silently.

    ⚠ THIS TEST USED TO ASSERT A REFUSAL and is kept, inverted, on purpose.
    It guarded the composed autonomous system raising `NotImplementedError`
    rather than quietly taking the dense route -- a performance surprise
    with no symptom, since the answer would be right and only the reason for
    asking would be gone. Now that the system is converted, the same risk
    runs the other way: a later refactor could reintroduce a silent dense
    fallback, and the shape of that failure is `spectral_radius` coming back
    NOT None, because only the dense path forms a monodromy.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cases = ((_q20_rlc(), 'trap', 1e-3, 100),        # driven plain
             (_rc_ladder(6), 'gear', 1e-3, 50),      # driven solved-history
             (_phase_circuit(), 'trap', 1e-3, 200),  # autonomous plain
             (_phase_circuit(), 'gear', 1e-3, 100))  # autonomous composed
    for cir, method, period, npts in cases:
        pss = PSS(cir, method=method, reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=period, timestep=period / npts,
                      maxiterations=25, matrix_free=True)
        assert pss.spectral_radius is None, \
            '%s/%s reported a spectral radius after a matrix-free solve, ' \
            'which means a monodromy was formed -- the dense path ran' \
            % (method, 'autonomous' if pss.autonomous else 'driven')


def test_a_parasitic_root_near_the_unit_circle_is_not_read_as_stability():
    """RECORDED SCOPE ITEM 3, and the case the whole item exists for.

    Gear-2's spurious root is `(1/3)^N` over a period, so for every method
    in this tree today `max |eig|` over the composed spectrum happens to
    pick the physical multiplier.  A method whose parasitic root sat near
    the unit circle would make that maximum report the DISCRETISATION'S OWN
    ARTEFACT as the orbit's stability -- a decaying orbit read as marginal,
    with nothing in the output to say so.

    No such method is in the tree, so the monodromy is SYNTHESISED with the
    structure the theory gives it: physical modes whose two halves agree
    (`v_{-1} = v_0`, one timestep apart on a smooth trajectory) and
    parasitic modes whose halves are opposed (`v_{-1} = -v_0`, the
    alternating root).  Waiting for such a method to exist before testing
    the separation would mean shipping it untested on the day it arrives.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_rc_ladder(2), method='gear', reltol=1e-6)
    m = pss.cir.n - 1

    rng = np.random.default_rng(7)
    Vp = rng.standard_normal((m, m))     # physical directions
    Vq = rng.standard_normal((m, m))     # parasitic directions
    ## halves EQUAL for physical, OPPOSED for parasitic
    V = np.block([[Vp, Vq], [Vp, -Vq]])
    phys_vals = np.linspace(0.80, 0.30, m)          # a decaying orbit
    para_vals = np.linspace(0.99, 0.95, m)          # spurious, near the circle
    lam = np.concatenate((phys_vals, para_vals))
    M = V @ np.diag(lam) @ np.linalg.inv(V)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        rho, phys, para = pss._spectral_report(M)

    ## the naive maximum -- what this used to report -- is the artefact
    naive = float(np.max(np.abs(np.linalg.eigvals(M))))
    assert abs(naive - 0.99) < 1e-8, \
        'the synthesised spectrum is not the case under test (max %r)' % naive

    assert abs(rho - 0.80) < 1e-8, \
        'spectral_radius is %r; the physical multiplier is 0.80 and 0.99 ' \
        'is the parasitic root -- the separation did not happen' % rho
    assert len(phys) == m and len(para) == m
    assert abs(para[0] - 0.99) < 1e-8

    ## and it says so, rather than separating silently
    assert any('PARASITIC' in str(w.message) for w in caught), \
        'a parasitic root within a decade of the physical spectrum was ' \
        'separated without warning'


def test_the_physical_block_split_shrinks_with_the_step():
    """Why the eigenvector's block structure is the right discriminator.

    The claim is not that physical and parasitic modes happen to differ; it
    is that they differ FOR A REASON THAT SCALES.  A physical mode's two
    halves are one timestep apart on a smooth trajectory, so they converge
    as `v_{-1} = e^{-lambda h} v_0 -> v_0`; a parasitic mode's are related
    by the method's spurious root and stay O(1) apart however fine the grid
    gets.  So refining `h` must shrink the physical split and leave the
    parasitic one alone -- which is a prediction, and this checks it.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def split_at(npts):
        pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=1e-3, timestep=1e-3 / npts, maxiterations=30)
        M = np.asarray(pss._monodromy)
        m = pss.cir.n - 1
        _ev, V = np.linalg.eig(M)
        s = np.sort(np.linalg.norm(V[m:, :] - V[:m, :], axis=0))
        return s[0], s[-1]          # smallest physical, largest parasitic

    p50, q50 = split_at(50)
    p200, q200 = split_at(200)

    ## the physical split falls with h, and at first order
    assert p200 < p50 / 2.5, \
        'the physical block split did not shrink with h (%.4f at N=50, ' \
        '%.4f at N=200) -- the discriminator rests on it doing so' % (p50, p200)
    ratio = p50 / p200
    assert 2.5 < ratio < 6.0, \
        'the physical split scaled by %.2f for a 4x refinement; first order ' \
        'in h predicts ~4 and this is the evidence the split means what it ' \
        'is taken to mean' % ratio
    ## while the parasitic one does not
    assert q200 > 0.5 and q50 > 0.5, \
        'parasitic splits collapsed too (%.3f, %.3f); nothing separates' \
        % (q50, q200)


def _scaled_vdp(mu=1.0, gain=1e4):
    """van der Pol with a VCVS copy of `v` scaled by `gain`.

    The copy is perfectly slaved, so the ORBIT is unchanged and only the
    arithmetic inside the phase-pin's `argmax` moves.  That makes it the
    instrument for asking whether comparing coordinates in different units
    is a defect.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    c.add_node('big')
    c['E'] = VCVS('v', gnd, 'big', gnd, g=gain)
    c['Rb'] = R('big', gnd, r=1e9)
    return c


def test_the_phase_pin_compares_units_on_purpose():
    """Why `argmax |x2 - x1|` is right to be dimensionally inconsistent.

    The phase row `e_k` exists to remove the orbit's own tangent direction
    from the singular `I - M`, and it can only do that in proportion to
    `|e_k . fhat|` -- its alignment with that direction.  `argmax |dx_k|`
    over ONE STEP is `argmax |f_k|` up to `h`, and `argmax |f_k|` is exactly
    `argmax |e_k . fhat|`.  So the rule maximises the very quantity the row
    needs, and it does so BECAUSE it compares raw magnitudes rather than
    despite it.

    ⚠ WRITTEN BECAUSE THE OBVIOUS FIX IS WRONG AND WAS ALMOST MADE.  Reading
    `argmax` over a vector mixing volts, amperes and a scaled copy looks
    like a units bug, and the natural repair -- normalise each coordinate by
    its own swing -- was measured here and picks a row **704x worse
    aligned** on this circuit.  The scaling that lets a large coordinate win
    the argmax is the same scaling that makes it dominate the flow vector;
    the two cancel exactly, which is why the raw comparison is the correct
    one.

    Seeded at `v`'s TURNING POINT, the hardest case for the rule: the
    pinned value then sits ~1e-3 of the way into the coordinate's range,
    i.e. nearly tangent to the orbit in that coordinate's own terms, and
    the row is STILL fully aligned because the scaled copy dominates `f`.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    T = 6.663293                      # measured free-running period, mu=1

    cir = _scaled_vdp()
    iref = cir.get_node_index(gnd)
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = Transient(cir, reltol=1e-8).solve(refnode=gnd, tend=35.0,
                                                timestep=0.02, x0=x0)
    t = np.asarray(res.sweep_values, dtype=float).ravel()
    X = np.asarray(res.x, dtype=float)
    W = X[:, t > 21.0]                 # a settled window, ~2 periods
    red = lambda f: np.concatenate((f[:iref], f[iref + 1:]))
    seed = red(W[:, int(np.argmax(W[iv]))])          # v at its turning point
    swing = np.array([np.ptp(W[i]) for i in range(W.shape[0]) if i != iref])

    ## replay the production rule for choosing the row
    pss = PSS(_scaled_vdp(), method='trap', reltol=1e-9)
    times, hs = pss._period_grid(T, 200, None)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss._begin_period(seed)
        x1 = np.asarray(pss.solve_timestep(seed, times[0], hs[0]), float)
        x2 = np.asarray(pss.solve_timestep(x1, times[1], hs[0],
                                           iq_last=pss._iq), float)
    step = np.abs(x2 - x1)
    k = int(np.argmax(step))
    k_norm = int(np.argmax(step / np.maximum(swing, 1e-300)))

    ## the flow at the seed, and each candidate row's alignment with it
    f = x1 - seed
    f = f / np.linalg.norm(f)
    align, align_norm = abs(f[k]), abs(f[k_norm])

    assert align > 0.9, \
        'the row `argmax |dx|` chose is only %.3e aligned with the flow; ' \
        'the whole defence of the rule is that this stays near 1' % align
    assert k_norm != k, \
        'the unit-normalised argmax picked the same row, so this circuit ' \
        'no longer separates the two rules and the test proves nothing'
    assert align_norm < 0.01, \
        'the unit-normalised row is %.3e aligned -- it was 1.4e-03 when ' \
        'this was written, and the point is that it is far worse' % align_norm
    assert align / align_norm > 100, \
        'the two rules are now within %.0fx; the measured gap was 704x' \
        % (align / align_norm)


def test_the_plain_matrix_free_matvec_carries_each_step_s_coefficients():
    """The plain path's `M v` and `dphi/dT`, against the dense traversal.

    ⚠ WRITTEN AROUND TWO BUGS THIS CAUGHT, and `trap` and `gear` are in the
    list because of them.  `_coeffs` is LIVE STATE and the MANUFACTURING
    step is order-dropped to Euler (`b = 0`) while the loop steps run the
    method's own pair -- measured, trapezoidal opens at
    `((49000, -49000), 0.0)` and runs at `((98000, -98000), -1.0)`.

    Applying the OPENING's pair to every step put the period column 40-50%
    out for `trap` and `gear`; applying the LOOP's to the opening seeded
    `Pq` non-zero where `_traverse` seeds it at zero, a 100% error for
    `trap`.  Both times the TRAJECTORY matched to zero, so nothing but a
    direct comparison against the dense sensitivity would have found them,
    and `euler` was exact under both -- a one-method test would have passed.

    ⚠ What this actually guards is the OPENING pair.  Inside the loop the
    coefficients are constant for every method in the tree, so a mutation
    swapping them for a post-run snapshot does NOT fail this -- checked.
    They are still stored per step, against a future variable-order method.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'euler', 'gear'):
        pss = PSS(_rc_ladder(12), method=method, reltol=1e-6)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            pss.solve(period=1e-3, timestep=1e-3 / 50, maxiterations=2)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(1e-3, 50, None)
        rng = np.random.default_rng(3)
        x_in = 0.01 * rng.standard_normal(m)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _x0, _xe, Mx, Mt = pss._traverse(x_in, 1e-3, times, hs,
                                             want_dT=True)
            opening, steps, x0f, xef, Mtf = pss._traverse_factored_plain(
                x_in, 1e-3, times, hs, want_dT=True)

        ## the trajectory must be the same one, or the rest is meaningless
        assert np.allclose(np.asarray(_x0), np.asarray(x0f), rtol=0, atol=0)
        assert np.allclose(np.asarray(_xe), np.asarray(xef), rtol=0, atol=0)

        ## the period column: one vector, computed once per Newton iteration
        want_t = np.asarray(Mt).ravel()
        err_t = (np.linalg.norm(np.asarray(Mtf).ravel() - want_t)
                 / max(np.linalg.norm(want_t), 1e-300))
        assert err_t < 1e-11, \
            '%s: dphi/dT differs from the dense column by %.3e' % (method, err_t)

        for _ in range(3):
            v = rng.standard_normal(m)
            got = pss._monodromy_matvec_plain(opening, steps, v)
            want = np.asarray(Mx) @ v
            err = np.linalg.norm(got - want) / np.linalg.norm(want)
            assert err < 1e-11, \
                '%s: matvec differs from the dense monodromy by %.3e' \
                % (method, err)


def test_a_matrix_free_plain_solve_agrees_with_the_dense_one():
    """RECORDED SCOPE ITEM 6 on the driven plain path (m columns, not 2m).

    Measured share 42.5% at m=502 and 60.7% at m=1002 -- about half the
    solved-history path's, as the column count says it must be -- and it
    delivers 1.34x and 1.63x end to end there.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'euler'):
        out = []
        for mf in (False, True):
            pss = PSS(_q20_rlc(), method=method, reltol=1e-6)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = pss.solve(period=1e-3, timestep=1e-3 / 100,
                                maxiterations=40, matrix_free=mf)
            out.append((pss, np.asarray(res['tpss'].x, dtype=float).ravel()))
        (pd, xd), (pm, xm) = out
        assert pd.converged and pm.converged, method
        rel = np.max(np.abs(xd - xm)) / np.max(np.abs(xd))
        assert rel < 1e-9, '%s: waveforms differ by %.3e' % (method, rel)


def test_a_matrix_free_autonomous_solve_agrees_with_the_dense_one():
    """The BORDERED system, matrix-free: unknowns `(x_0, T)`.

    `dphi/dT` is one column and does not depend on the Krylov direction, so
    the trajectory pass computes it once per Newton iteration and the matvec
    reads it:  `J [v; s] = [ (I - M) v - s dphi/dT ; v_k ]`.

    The PERIOD is asserted as well as the waveform, because it is the
    unknown the border exists for -- a matvec that dropped the `dphi/dT`
    term entirely would still close the state equations and return a wrong
    period.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    out = []
    for mf in (False, True):
        pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=1e-3, timestep=1e-3 / 200,
                            maxiterations=30, matrix_free=mf)
        out.append((pss, np.asarray(res['tpss'].x, dtype=float).ravel()))
    (pd, xd), (pm, xm) = out
    assert pd.converged and pm.converged
    assert np.max(np.abs(xd - xm)) / np.max(np.abs(xd)) < 1e-9
    assert abs(pm.period - pd.period) < 1e-12 * pd.period, \
        'matrix-free solved T=%.12g against the dense %.12g' \
        % (pm.period, pd.period)


def test_a_matrix_free_composed_autonomous_solve_agrees_with_the_dense_one():
    """BOTH enlargements at once, matrix-free: `(x_0, x_{-1}, T)`.

    The last of the four systems.  With `w = (v_0, v_{-1}, s)` and `M v` the
    `2m` pair map,

        J w = [ v - M v - s (dx_{N-1}/dT, dx_{N-2}/dT) ;  v[k] ]

    -- ONE phase row, pinning the `x_0` block only, because time translation
    slides both states along the orbit together and the freedom stays
    one-dimensional.

    ⚠ BOTH period columns are asserted against the dense traversal before
    the solve is.  They are the part with no analogue in the other three
    systems, and a solve can absorb a wrong one by moving `x_0` instead --
    returning a converged answer at a period that is quietly off.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=25)
    m = pss.cir.n - 1
    T = pss.period
    times, hs = pss._period_grid(T, 100, None)
    rng = np.random.default_rng(11)
    a, b = 0.01 * rng.standard_normal(m), 0.01 * rng.standard_normal(m)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        xl, xp, Pl, Pp, Ptl, Ptp = pss._traverse_solved_history(
            a, b, times, hs, T=T, want_dT=True)
        C0, st, xlf, xpf, Ptlf, Ptpf = pss._traverse_factored(
            a, b, times, hs, T=T, want_dT=True)

    rel = lambda g, w: (np.linalg.norm(np.asarray(g).ravel()
                                       - np.asarray(w).ravel())
                        / max(np.linalg.norm(np.asarray(w).ravel()), 1e-300))
    assert rel(Ptlf, Ptl) < 1e-11, 'dx_{N-1}/dT differs by %.3e' % rel(Ptlf, Ptl)
    assert rel(Ptpf, Ptp) < 1e-11, 'dx_{N-2}/dT differs by %.3e' % rel(Ptpf, Ptp)
    M = np.vstack((Pl, Pp))
    for _ in range(3):
        v = rng.standard_normal(2 * m)
        assert rel(pss._monodromy_matvec(C0, st, v), M @ v) < 1e-9

    ## and the solve itself, period included
    out = []
    for mf in (False, True):
        p = PSS(_phase_circuit(), method='gear', reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = p.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=25,
                        matrix_free=mf)
        out.append((p, np.asarray(r['tpss'].x, dtype=float).ravel()))
    (pd, xd), (pm, xm) = out
    assert pd.converged and pm.converged
    assert np.max(np.abs(xd - xm)) / np.max(np.abs(xd)) < 1e-9
    assert abs(pm.period - pd.period) < 1e-12 * pd.period, \
        'matrix-free solved T=%.12g against the dense %.12g' \
        % (pm.period, pd.period)


def test_pss_forwards_its_solver_strategies_to_the_inner_transient():
    """`nrsolver`, `linearsolver` and `scaler` must reach the thing that solves.

    ⚠ THEY DID NOT, AND THE CONSTRUCTOR ACCEPTED THEM ANYWAY.  All three are
    declared on the base `Analysis`, so `PSS(cir, linearsolver=SuperLUSolver())`
    has always been valid input -- and `_transient()` forwarded the tolerances
    and not these, so the inner `Transient` resolved to `DenseSolver` and
    `StandardNewton` whatever the caller asked for.

    That is the THIRD parameter this class has accepted and never read
    (`method` declared and never read; `analysis='PSS'` matching nothing),
    and each one was invisible for the same reason: the constructor validates
    the name, nothing checks the boundary.  So this asserts on the RESOLVED
    OBJECT the inner analysis will actually use -- a test that the string
    'linearsolver' appears in the source passes with the parameter dropped.
    """
    from pycircuit.circuit.linearsolver import SuperLUSolver, DenseSolver
    from pycircuit.circuit.nrsolver import DampedNewton, StandardNewton
    circuit.default_toolkit = circuit.numeric

    pss = PSS(_q20_rlc(), method='trap', reltol=1e-6,
              linearsolver=SuperLUSolver(), nrsolver=DampedNewton())
    tran = pss._transient()
    assert isinstance(tran._get_linearsolver(), SuperLUSolver), \
        'the inner Transient resolved to %r, not the SuperLUSolver asked for' \
        % tran._get_linearsolver()
    assert isinstance(tran._get_nrsolver(), DampedNewton), \
        'the inner Transient resolved to %r, not the DampedNewton asked for' \
        % tran._get_nrsolver()

    ## and the default is still the historical dense path, unmoved
    plain = PSS(_q20_rlc(), method='trap', reltol=1e-6)._transient()
    assert isinstance(plain._get_linearsolver(), DenseSolver)
    assert isinstance(plain._get_nrsolver(), StandardNewton)


def test_a_failed_inner_krylov_solve_says_so():
    """A Krylov failure must name itself, not become 'No convergence'.

    ⚠ THE INNER SOLVE'S VERDICT USED TO BE DISCARDED.  `gmres` returns an
    `info` saying whether it converged, broke down, or ran out of cycles,
    and dropping it meant a matrix-free run could spend its whole budget --
    scipy's `maxiter` counts RESTART CYCLES, so the pair multiplies, and the
    earlier `restart=200, maxiter=400` was a worst case near 80 000 matvecs,
    each a full replay of the period -- and then report the same generic
    outer message as a circuit that simply needed another iteration.

    This file's standard is that a failure says what happened, which is why
    `T = 0` and the singular free-period Jacobian are named.  The inner
    solve is now held to it too.  Starved of cycles on purpose, because the
    honest way to test a failure path is to cause the failure.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_rc_ladder(20), method='gear', reltol=1e-6)
    old_cycles, old_restart = pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART
    try:
        ## one cycle of a one-dimensional Krylov space cannot solve this
        pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART = 1, 1
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=1e-3, timestep=1e-3 / 25, maxiterations=5,
                      matrix_free=True)
    finally:
        pss.KRYLOV_MAX_CYCLES, pss.KRYLOV_RESTART = old_cycles, old_restart

    assert not pss.converged, \
        'a starved inner solve still reported convergence'
    msgs = [str(w.message) for w in caught]
    assert any('inner solve did not converge' in m for m in msgs), \
        'the Krylov failure was not named; warnings were %r' % msgs
    assert any('matrix_free=False' in m for m in msgs), \
        'the message does not say what to do instead'


def test_pss_refuses_a_circuit_carrying_hidden_state():
    """A `TLine` under PSS was SILENTLY A SHORT, and reported success.

    ⚠ THIS IS KUNDERT'S HIDDEN STATE, and the tree already knew the defect
    by name -- `transient.py` carries a comment about the identical failure
    being found and fixed THERE ("TLine.G saw an empty buffer and stamped
    the line as a DC SHORT").  PSS reproduced it, because `TLine.history` is
    filled by `cir.accept_step`, which a forward transient calls at every
    accepted step and which PSS never calls: it drives `solve_timestep`
    directly.

    Measured on a quarter-wave open stub before the refusal went in:

        PSS  converged       True
        PSS  spectral_radius 0.0        <- a circuit reporting no state
        PSS  warnings        NONE
        PSS  amp v(b)        0.999969   <- line absent, source sees an open
        TRAN amp v(b)        0.244201   <- line active

    ⚠ AND FILLING THE BUFFER IS NOT THE FIX, which is why this asserts a
    REFUSAL and not a number.  Calling `accept_step` per step would populate
    the history and make `phi` genuinely history-dependent, so the monodromy
    would be the derivative of a neighbouring problem -- the exact thing
    `_begin_period` exists to prevent.  `_begin_period` resets what is in
    `x`; no reset of the integrator rings can make a period map that is not
    a function of `x_0` into one that is.
    """
    from pycircuit.circuit import TLine, Node
    circuit.default_toolkit = circuit.numeric
    per = 1e-9
    c = SubCircuit()
    c['vs'] = VSin(Node('a'), gnd, va=1.0, freq=1.0 / per)
    c['rs'] = R(Node('a'), Node('b'), r=50.0)
    c['t1'] = TLine(Node('b'), gnd, Node('c'), gnd, Z0=50.0, TD=per / 4.0)

    ## the element declares it, and the circuit finds it by instance name
    assert c.hidden_state_elements() == ['t1']

    with pytest.raises(NotImplementedError, match='HIDDEN STATE'):
        PSS(c).solve(period=per, timestep=per / 200)

    ## ⚠ and the refusal must not sweep in ordinary circuits.  Overriding
    ## `accept_step` is NOT the criterion -- `_WrapEvents` and the HDL
    ## `@cross` wrapper do, and feed only `next_event`, which PSS ignores
    ## because it imposes its own grid.  A diode's `_vlim` is hidden state
    ## too and is self-erasing (the inner Newton drives it to `v` at
    ## convergence), so it is not declared either.
    assert _q20_rlc().hidden_state_elements() == []
    pss = PSS(_q20_rlc(), method='trap', reltol=1e-6)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=20)
    assert pss.converged


def test_the_period_column_is_a_total_derivative_not_a_partial():
    """RECORDED SCOPE: Gear-2's `d(phi)/dT` was exactly 3/2 too large.

    ⚠ THE CAUSE IS A RESULT CARRIED ACROSS CONTEXTS.  `residual_dh` is
    Fang's `p`, `d/dh_m` with the past steps HELD FIXED -- correct for the
    coupled time-stepping method it was written for, where they are fixed.
    `bdf2_alphas_dh`'s own docstring says so: "h2 is a past step and is held
    fixed".  A shooting analysis solving for the period rebuilds its grid at
    the current `T`, so every step is `c_k T` and `h_{n-1}` moves too.

    Euler and trapezoidal never noticed -- their coefficients depend on
    `h_n` alone, so the partial IS the total -- which is exactly why only
    Gear-2 was hit, and why a one-method test would have found nothing.
    Measured before the fix, the ratio of the code's column to the true one
    was 1.4859 / 1.4939 / 1.4972 at 100 / 200 / 400 points, converging on
    3/2, with the residue after dividing by 3/2 halving per refinement.

    The fix is one term and no new derivative: the `alphas` are homogeneous
    of degree -1 in the step sizes, so `T d(iq)/dT = -sum_k a_k q_{n-k}`
    exactly (`Integrator.companion_dT`).

    ⚠ TRAPEZOIDAL IS STILL O(h) HERE AND THAT IS A DIFFERENT DEFECT.  Its
    autonomous route is the PLAIN path, whose `Pt` opens at zero although
    `x_0` is manufactured by a step of size `c_0 T` and does depend on T.
    That is the plain path's seeding, not this; asserted loosely on purpose
    so it does not silently tighten.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cap = {}
    orig = PSS._free_period_solve

    def grab(self, func, z0, abstol, xtol, reltol, maxiter, seed_period,
             solver=None):
        cap['func'] = func
        cap['z0'] = np.asarray(z0, dtype=float).copy()
        return orig(self, func, z0, abstol, xtol, reltol, maxiter,
                    seed_period, solver=solver)

    try:
        PSS._free_period_solve = grab
        got = {}
        for method in ('gear', 'trap'):
            pss = PSS(_phase_circuit(), method=method, reltol=1e-9)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
            func, z = cap['func'], cap['z0'].copy()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _F, J = func(z)
                d = 1e-7 * abs(z[-1])
                zp, zm = z.copy(), z.copy()
                zp[-1] += d
                zm[-1] -= d
                fp, _ = func(zp)
                fm, _ = func(zm)
            fd = (np.asarray(fp, float) - np.asarray(fm, float)) / (2 * d)
            col = np.asarray(J, float)[:, -1]
            got[method] = (np.linalg.norm(col - fd)
                           / max(np.linalg.norm(fd), 1e-300))
    finally:
        PSS._free_period_solve = orig

    ## Gear-2's column is now the total derivative, to the FD floor
    assert got['gear'] < 1e-7, \
        "Gear-2's period column is %.3e from finite differences; it was " \
        '4.9e-01 when it used the partial, and 3/2 of the truth' % got['gear']
    ## trapezoidal unchanged, and still carrying the plain path's seed error
    assert got['trap'] < 5e-2, \
        "trapezoidal's period column is %.3e, worse than the O(h) seeding " \
        'error it is expected to carry' % got['trap']
