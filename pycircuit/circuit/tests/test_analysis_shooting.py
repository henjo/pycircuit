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
    src = inspect.getsource(shooting.PSS.solve)
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
    augmented (x, iq) state and it converges in 6 iterations, so that case
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
# Phase 3: the augmented (x, iq) monodromy
# ---------------------------------------------------------------------------

def test_trapezoidal_shooting_converges_with_the_augmented_state():
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
    augmented state, trap converges AND lands at 19.990 V -- 0.05% of
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
