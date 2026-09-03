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


def test_PAC_runs_on_a_nonlinear_circuit_without_forming_the_operator():
    """⚠ PAC IS NO LONGER WITHDRAWN — this test used to assert that it was.

    Stage 11 withdrew it for forming the whole `(N m) x (N m)` operator
    densely: 419.5 GiB at `N = 137`, `m = 1000`. The withdrawal note said
    the body was "the starting point for a matrix-free rewrite", and that
    turned out to be exactly right — the operator in it was CORRECT, just
    un-preconditioned. See `PAC`'s docstring and
    `test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`.

    This keeps the withdrawal's own circuit — a diode, so genuinely
    nonlinear and genuinely time-varying, which is the case PAC exists for
    and the one the AC gate cannot cover — and asserts the two things the
    rewrite promises: it runs, and it never allocates the operator.

    ⚠ THE MEMORY ASSERTION IS THE POINT, not the numbers. The dense route
    would allocate `(N m)^2` complex entries; this circuit is small enough
    that doing so would succeed and pass a value check silently. So the
    test counts what is stored: `N` factorisations of an `m x m` matrix,
    which is `O(N m^2)` and not `O(N^2 m^2)`.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    pss = PSS(cir, method='gear')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / 40)
    assert pss.converged

    pac = PAC(cir)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pac.solve(pss, np.array([1e6, 3e6]))
    X = np.asarray(res.x, dtype=complex)
    assert np.isfinite(X).all(), 'PAC returned non-finite entries'
    assert np.abs(X).max() > 0, \
        'PAC returned all zeros on a driven circuit -- the source vector is ' \
        'probably being read positionally again (see the AC gate)'

    ## the operator is never formed: what is stored is per-STEP, m x m
    fp = pss.factored_period()
    m = cir.n - 1
    N = len(fp.steps)
    assert N > 1 and fp.width <= 2 * m, \
        'the factored period should hold N per-step m x m factorisations, ' \
        'not one (N m) x (N m) matrix'
    assert pac.matvecs < 4 * fp.width, \
        'PAC used %d matvecs for 2 frequencies on an m=%d circuit; forming ' \
        'the monodromy would take %d, and the point is not to' \
        % (pac.matvecs, m, fp.width)


@unittest.skip("Superseded by the PAC gate tests at the end of this file")
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


def test_autonomy_is_decided_on_every_grid_point_not_a_stride():
    """A narrow pulse must not read as a DC circuit.

    ⚠ `_is_autonomous` sampled `times[1::len(times)//8]` -- about nine
    points -- while its docstring called the test "exact where a spectral
    test is not".  Nine points cannot see a narrow pulse.  Measured on an RC
    driven by a `VPulse` placed BETWEEN two of those samples: 40% and 20%
    duty were read correctly, and 5%, 1% and 0.5% all came back AUTONOMOUS.

    What that costs is not a warning: an autonomous verdict routes the solve
    to the FREE-PERIOD system, which solves for `T` and DISCARDS the period
    the caller asked for.  `DEGENERATE_PERIOD_FACTOR` cannot catch it either
    -- it tests the magnitude of `T`, not whether the circuit was driven.
    PWM, sampling clocks, S/H and mixer LOs are core PSS workload and are
    exactly the shapes a stride misses.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def pulsed(duty):
        ## td places the pulse between the samples the old stride took
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=per * 0.19,
                         tr=per * 1e-4, tf=per * 1e-4,
                         pw=per * duty, per=per)
        c['r'] = R('a', 'b', r=1e3)
        c['c'] = C('b', gnd, c=1e-9)
        return c

    for duty in (0.40, 0.05, 0.01, 0.005):
        pss = PSS(pulsed(duty), method='trap', reltol=1e-6)
        times, _hs = pss._period_grid(per, 200, None)
        assert not pss._is_autonomous(times), \
            'a %.1f%% duty pulse read as autonomous; it would be solved at ' \
            'a period of the analysis\'s own choosing' % (100 * duty)

    ## and the verdict still holds for a genuinely autonomous circuit
    pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
    times, _hs = pss._period_grid(1e-3, 200, None)
    assert pss._is_autonomous(times)

    ## end to end: the requested period is the one that comes back
    pss = PSS(pulsed(0.01), method='trap', reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 200, maxiterations=30)
    assert pss.converged
    assert abs(pss.period - per) < 1e-15 * per, \
        'a driven solve returned T=%.9g for a requested %.9g' \
        % (pss.period, per)


def test_one_reference_node_per_analysis():
    """The traversal and the reported waveform must agree on which node is 0.

    `self.irefnode` is fixed in `__init__` and is what the TRAVERSAL uses --
    `_transient.irefnode`, every `remove_row_col`, the monodromy's shape.
    `solve()` computed its OWN from `refnode=` and used that to reinsert the
    zero row into the result.  The two were never compared, so
    `PSS(cir).solve(refnode='b')` solved against ground and reported against
    `b`: each row sensible alone, the set incoherent, with ground itself
    coming back non-zero.

    Refused rather than silently rotated, because there is no answer to
    give -- the two choices disagree about which variable was eliminated
    before the solve began.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['r'] = R('a', 'b', r=1e3)
    c['c'] = C('b', gnd, c=1e-7)

    with pytest.raises(ValueError, match='different reference node'):
        PSS(c).solve(refnode='b', period=1e-3, timestep=1e-5)

    ## the default still works, and so does agreeing explicitly
    pss = PSS(c)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-5, maxiterations=20)
    assert pss.converged


def test_the_alternating_mode_is_l_stability_not_the_iq_recursion():
    """Why trapezoidal shooting dies at even step counts -- the real cause.

    The class docstring blamed trapezoidal's `iq` RECURSION
    (`iq_n = ... - iq_{n-1}`, homogeneous mode `(-1)^n`).  That is the wrong
    cause: trapezoidal is A-stable but NOT L-stable, so it maps the null
    space of the singular MNA `C` by exactly `-1` per step, and an X-ONLY
    formulation carrying no `iq` variable at all is singular in the same
    way.  For `C x' + G x + u = 0`,

        A_trap  = (C/h + G/2)^-1 (C/h - G/2)   ->  -I  on null(C)
        A_euler = (C/h + G)^-1 (C/h)           ->   0  on null(C)

    so the count of `-1` modes is exactly `m - rank(C)`, one per ALGEBRAIC
    variable, on every MNA circuit.  That is a prediction with a number in
    it, which is why it is worth a test: it says the opening step must be
    L-STABLE and nothing more specific, so any L-stable opener rescues
    exactly those modes and Euler is not privileged.
    """
    circuit.default_toolkit = circuit.numeric
    h = 1e-3 / 200
    for name, cir in (('Q=20 RLC', _q20_rlc()), ('RC', _rc_ladder(1))):
        n = cir.n
        Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)
        Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)
        iref = cir.get_node_index(gnd)
        keep = [i for i in range(n) if i != iref]
        Cm = Cm[np.ix_(keep, keep)]
        Gm = Gm[np.ix_(keep, keep)]
        m = len(keep)
        rank_c = np.linalg.matrix_rank(Cm)
        algebraic = m - rank_c
        assert algebraic > 0, '%s has no algebraic variables to test' % name

        e_trap = np.linalg.eigvals(np.linalg.solve(Cm / h + Gm / 2.0,
                                                   Cm / h - Gm / 2.0))
        e_eul = np.linalg.eigvals(np.linalg.solve(Cm / h + Gm, Cm / h))
        n_minus1 = int(np.sum(np.abs(e_trap + 1.0) < 1e-9))
        n_zero = int(np.sum(np.abs(e_eul) < 1e-9))
        assert n_minus1 == algebraic, \
            '%s: trapezoidal has %d modes at -1, m - rank(C) = %d' \
            % (name, n_minus1, algebraic)
        assert n_zero == algebraic, \
            '%s: Euler has %d modes at 0, m - rank(C) = %d' \
            % (name, n_zero, algebraic)


def _resonator_at_resonance():
    """A Q=20 series resonator DRIVEN AT ITS RESONANCE, analytic peak 20 V."""
    Lv, Cv, Q = 1e-3, 1e-9, 20.0
    f0 = 1.0 / (2 * np.pi * np.sqrt(Lv * Cv))
    c = SubCircuit()
    c.add_node('n1')
    c.add_node('n2')
    c['vs'] = VSin(gnd, 'n1', va=1.0, freq=f0)
    c['L'] = L('n1', 'n2', L=Lv)
    c['C'] = C('n2', gnd, c=Cv)
    c['R'] = R('n1', 'n2', r=Q * np.sqrt(Lv / Cv))
    return c, 1.0 / f0


def test_a_grid_that_outruns_zero_stability_says_so():
    """A caller's grid can demote Gear-2 to first order, silently.

    `_period_grid` validated positivity and sum-to-1 and NOTHING about the
    interior ratios.  A two-step method is zero-stable only to
    `h_n/h_{n-1} = 1 + sqrt(2)`; past it the integrator drops the step to
    Euler -- correct, and invisible.  Measured on this resonator with an
    alternating 3:1 grid, half the steps demoted:

          npts   uniform    3:1 grid     (analytic peak 20 V)
           100   19.91489    7.99821
           800   20.02443   16.85280

    60% low at 100 points, crawling up at FIRST order, `converged = True`
    every time.  ⚠ And refining does not fix it -- a refined 3:1 grid is
    still 3:1 -- which is why the warning names the RATIO rather than
    advising a smaller step.

    ⚠ This is the premise item 5 removed.  The class docstring's literature
    note argues Wambacq's objections to non-uniform BDF "do not bite inside
    a run" because the grid is uniform and frozen; a caller's grid is frozen
    but not uniform.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir, per = _resonator_at_resonance()
    n = 100
    bad = np.where(np.arange(n) % 2 == 0, 3.0, 1.0)
    bad = bad / bad.sum()
    smooth = _grid_fracs('smooth', n)

    def warns(method, fr):
        pss = PSS(cir, method=method, reltol=1e-8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss._period_grid(per, n, fr)
        return [w for w in caught if 'zero-stable' in str(w.message)]

    assert warns('gear', bad), \
        'a 3:1 grid on a two-step method was accepted without a word'
    ## and it does not cry wolf
    assert not warns('gear', smooth), 'a smooth grid should not warn'
    assert not warns('trap', bad), \
        'a ONE-step method is not subject to the two-step bound'

    ## the damage the warning is about, so the number is pinned too
    peaks = {}
    for label, fr in (('uniform', None), ('3:1', bad)):
        pss = PSS(cir, method='gear', reltol=1e-8)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=per / n, grid=fr,
                            maxiterations=30)
        v = np.asarray(res['tpss'].v('n2', gnd), dtype=float).ravel()
        peaks[label] = 0.5 * (v.max() - v.min())
        assert pss.converged, '%s did not converge' % label
    assert abs(peaks['uniform'] - 20.0) < 0.2, \
        'the uniform reference is %.4f, not the analytic 20 V' % peaks['uniform']
    assert peaks['3:1'] < 0.6 * peaks['uniform'], \
        'the 3:1 grid gave %.4f against the uniform %.4f -- if these now ' \
        'agree, the demotion is gone and this test has lost its subject' \
        % (peaks['3:1'], peaks['uniform'])


def test_an_index_2_solve_that_converges_is_correct():
    """Index-2 circuits: what actually happens, measured on four topologies.

    An external review reported that index-2 circuits break the DEFAULT
    method, with "89% error on the v-source branch current", and proposed
    detecting index > 1 and refusing.  Reproduced and then measured
    further, and the conclusion does not survive:

        circuit      index>1   trap   euler   gear
        CV-loop        yes       F      F      T
        CV-loop2       yes       X      F      F
        LI-cutset      yes       T      T      T
        LI-cutset2     yes       X      X      F

    ⚠ INDEX > 1 IS NOT PREDICTIVE OF FAILURE.  Every one of those is index
    2 by the projector test (`N^T G N` singular, `N` a null basis of `C`),
    and on `LI-cutset` all three methods converge.  A refusal keyed on index
    would reject circuits that work.

    ⚠ AND GEAR IS NOT THE WORKAROUND.  It succeeds on `CV-loop` -- the one
    circuit the review tested -- and FAILS on two of the other three.  Had
    the refusal named `method='gear'` as the remedy, as proposed, it would
    have sent users to a method that does not generalise.

    ⚠ WHAT THIS PINS is the part that would matter if it broke: when an
    index-2 solve reports CONVERGED, its answer is right.  Measured against
    a settled transient over every state, algebraic ones included, the worst
    relative error was 6.6e-05 (CV-loop/gear) and 1.3e-05 / 1.7e-04 /
    7.4e-05 (LI-cutset trap/euler/gear).  The failures are LOUD --
    `converged = False` or an exception -- not silent wrong answers, which
    is why this is documented rather than refused.  A regression to a
    silently wrong algebraic variable is exactly what this test would catch.

    ⚠ WHAT IS PINNED IS PINNED FOR THE LINEAR CASE.  All four circuits are
    LINEAR, so `rank(C)` and the constraint structure are constant along the
    orbit.  On a NONLINEAR index-2 circuit both can vary with the operating
    point, so the index itself can change around the period -- and a solve
    could then be consistent at the points it checks and inconsistent
    between them, which is the one shape a converged-and-wrong answer could
    still take.  Offered as the gap rather than as a finding: no such
    circuit has been built here, and neither reviewer nor this test has
    evidence one exists in this tree.
    """
    import warnings
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def cv_loop():
        c = SubCircuit()
        c.add_node('a')
        c.add_node('b')
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per)
        c['c1'] = C('a', 'b', c=1e-9)
        c['c2'] = C('b', gnd, c=1e-9)
        c['r'] = R('b', gnd, r=1e9)
        return c

    ## it really is index 2: N^T G N is singular for N a null basis of C
    cir = cv_loop()
    n = cir.n
    Cm = np.asarray(cir.C(np.zeros(n)), dtype=float)
    Gm = np.asarray(cir.G(np.zeros(n)), dtype=float)
    iref = cir.get_node_index(gnd)
    keep = [i for i in range(n) if i != iref]
    Cm, Gm = Cm[np.ix_(keep, keep)], Gm[np.ix_(keep, keep)]
    _U, sv, Vt = np.linalg.svd(Cm)
    d = int(np.sum(sv > len(keep) * sv[0] * np.finfo(float).eps))
    N = Vt[d:].T
    s2 = np.linalg.svd(N.T @ Gm @ N, compute_uv=False)
    ## ⚠ the projection can be identically ZERO, which is singular in the
    ## strongest sense -- so the ratio needs a guarded denominator rather
    ## than `s2[-1]/s2[0]`, which is 0/0 there.  Measured: this circuit
    ## gives exactly that.
    ratio = float(s2[-1] / max(s2[0], 1e-300))
    assert ratio < 1e-10, \
        'this circuit is no longer index 2 (sigma_min/sigma_max = %.3e), ' \
        'so the test has lost its subject' % ratio

    ## reference: a settled transient
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ## ⚠ a SHORT reference on purpose.  R*C here is ~0.5 s against a
        ## 1 ms period, so no affordable transient settles the slow DC
        ## drift -- but the quantity compared is the per-period AMPLITUDE,
        ## which the capacitive divider sets within a couple of periods.
        ## Sixty periods cost 52 s and bought no accuracy over eight.
        rt = Transient(cv_loop(), reltol=1e-10).solve(
            refnode=gnd, tend=per * 8, timestep=per / 200)
    t = np.asarray(rt.sweep_values, dtype=float).ravel()
    X = np.asarray(rt.x, dtype=float)
    W = X[:, t > t[-1] - per]
    ref = np.array([0.5 * (W[i].max() - W[i].min()) for i in range(W.shape[0])])

    pss = PSS(cv_loop(), method='gear', reltol=1e-8)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = pss.solve(period=per, timestep=per / 200, maxiterations=40)
    assert pss.converged, 'gear no longer converges on the CV-loop'
    Xp = np.asarray(res['tpss'].x, dtype=float)
    got = np.array([0.5 * (Xp[i].max() - Xp[i].min())
                    for i in range(Xp.shape[0])])
    k = min(len(got), len(ref))
    worst = np.max(np.abs(got[:k] - ref[:k])
                   / np.maximum(np.abs(ref[:k]), 1e-12))
    assert worst < 1e-3, \
        'a CONVERGED index-2 solve disagrees with a settled transient by ' \
        '%.3e over its states -- that is the silent wrong answer this ' \
        'exists to catch' % worst


def test_the_returned_waveform_closes_on_a_non_uniform_grid():
    """The replay must walk the same `(t, h)` pairs the traversal did.

    The two traversals pair them differently: `_traverse_solved_history`
    walks `times[1:]` with `hs[_j]`, while the plain `_traverse` takes the
    MANUFACTURING step at `(times[0], hs[0])` first and only then walks
    `times[1:]` with `hs[_j]` -- so in the plain case the step after the
    opening one uses `hs[0]` again, not `hs[1]`.  The replay set
    `walk = times` and indexed `hs[min(_j, ...)]`, pairing `times[k]` with
    `hs[k]` from k=1 on: off by one against the traversal.

    ⚠ A UNIFORM GRID HIDES IT COMPLETELY, since every `hs` is the same
    number -- which is why it survived every uniform test here.  Measured on
    this resonator at 200 points, closure `|x(T) - x(0)|` of the RETURNED
    waveform, with `converged = True` in every row:

          grid      trap (before)   trap (after)   gear (control)
          uniform     5.61e-13        5.61e-13       4.62e-14
          4:1         1.70e-02        1.44e-11       5.33e-15
          16:1        4.88e-03        3.55e-14       1.78e-14

    Gear closes on every grid either way because it takes the
    solved-history branch, whose pairing was already right -- which is what
    made it a clean control rather than a second unknown.

    ⚠ AND THE HEADLINE ppm FIGURES WERE NEVER AFFECTED, which is worth
    stating because the opposite was assumed when this was found.  The
    period comes from the shooting SOLVE, not the replay: van der Pol still
    reads -73.8 ppm (trap) and -100.6 ppm (gear) after the fix.  What the
    shift corrupted is the returned waveform, `times`, and the LTE reports
    -- and every recorded LTE figure was taken on a uniform grid, so no
    recorded number moved.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir, per = _resonator_at_resonance()
    npts = 200

    def closure(method, ratio):
        if ratio == 1.0:
            g = None
        else:
            fr = np.where(np.arange(npts) % 2 == 0, ratio, 1.0)
            g = fr / fr.sum()
        pss = PSS(cir, method=method, reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=per / npts, grid=g,
                            maxiterations=40)
        assert pss.converged, '%s/%s did not converge' % (method, ratio)
        X = np.asarray(res['tpss'].x, dtype=float)
        return float(np.max(np.abs(X[:, -1] - X[:, 0])))

    ## the plain path, where the pairing was wrong
    for ratio in (1.0, 4.0, 16.0):
        c = closure('trap', ratio)
        assert c < 1e-8, \
            'the returned waveform does not close on a %.0f:1 grid ' \
            '(|x(T)-x(0)| = %.4e) -- it is not the solution the residual ' \
            'was driven to zero on' % (ratio, c)

    ## and the solved-history control, which was always right
    for ratio in (1.0, 4.0, 16.0):
        assert closure('gear', ratio) < 1e-8


def test_x0_unknown_solves_for_the_period_s_own_start():
    """`x0_unknown=True`: solve for `x_0`, manufacture nothing.

    The default plain path hands `fsolve` the PRE-IMAGE of a manufactured
    opening step while returning a Jacobian taken with respect to `x_0` --
    a frame error, and the reason the true `dF/dx_in` is singular and the
    iteration is a contraction rather than a Newton.

    ⚠ TRAPEZOIDAL STILL NEEDS AN L-STABLE OPENER, which is why this moves
    the Euler step INSIDE the period rather than removing it.  Without one
    the period map is `A_trap^K`, singular at EVEN K on every MNA circuit --
    the `(-1)^n` obstruction, and the fourth design to meet it.  Checked on
    the model problem before any of this was built: `sigma_min(I - A^K)` is
    exactly 0.0 for bare trapezoidal at K=100 and 200 and 6.0e-03 with the
    Euler opener.  Hence the even/odd point counts below: they are the
    recorded falsifier for that whole family of reformulations.

    What it costs and buys, both measured:

        Q=20 at resonance, analytic 20 V     default    x0_unknown
          100 points                         20.01273    19.76939
          200 points                         20.02208    19.96123

    -- the in-period Euler step degrades the ORBIT, worse on a benign
    uniform grid.  And on van der Pol's own 1105-step LTE grid the sign
    flips: -47.3 ppm against the default's -73.8.

    ⚠ THAT GAIN IS NOT THE FORMULATION ALONE and the first attribution here
    was wrong.  It comes from the formulation making the opening-step
    SUBDIVISION unnecessary: that subdivision exists only to protect a
    manufactured step taken from an iterate that may be far from the orbit,
    and measured on the same grid it costs -47.3 -> -73.8.  With `x_0` the
    unknown the first step starts ON the orbit and the raw grid is solvable.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir, per = _resonator_at_resonance()

    got = {}
    for method in ('trap', 'euler'):
        for npts in (100, 101, 200):
            for flag in (False, True):
                pss = PSS(cir, method=method, reltol=1e-10)
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    res = pss.solve(period=per, timestep=per / npts,
                                    maxiterations=40, x0_unknown=flag)
                assert pss.converged, \
                    '%s/%d/x0_unknown=%s did not converge -- an EVEN point ' \
                    'count failing here is the (-1)^n mode returning' \
                    % (method, npts, flag)
                X = np.asarray(res['tpss'].x, dtype=float)
                closure = float(np.max(np.abs(X[:, -1] - X[:, 0])))
                assert closure < 1e-8, \
                    '%s/%d/x0_unknown=%s: the returned waveform does not ' \
                    'close (%.3e)' % (method, npts, flag, closure)
                v = np.asarray(res['tpss'].v('n2', gnd), dtype=float).ravel()
                got[(method, npts, flag)] = 0.5 * (v.max() - v.min())

    ## EULER IS THE CONTROL: its manufacturing step IS an Euler step, so the
    ## two formulations describe the same map and must agree to the digit.
    for npts in (100, 101, 200):
        a, b = got[('euler', npts, False)], got[('euler', npts, True)]
        assert abs(a - b) < 1e-9 * max(abs(a), 1.0), \
            'euler at %d points gives %.6f without and %.6f with the flag; ' \
            'these describe the same map and must agree' % (npts, a, b)

    ## trapezoidal genuinely differs, and in the direction measured
    assert got[('trap', 100, True)] < got[('trap', 100, False)], \
        'the in-period Euler step no longer costs accuracy on a uniform ' \
        'grid -- that cost is why this is an option and not the default'
    assert abs(got[('trap', 200, True)] - 20.0) \
        < abs(got[('trap', 100, True)] - 20.0), \
        'the x0_unknown error does not shrink with refinement'

    ## and the combination with nothing to change is refused, not ignored
    with pytest.raises(NotImplementedError, match='nothing to change'):
        PSS(cir, method='gear').solve(period=per, timestep=per / 100,
                                      x0_unknown=True)

    ## ⚠ AND IT COMPOSES WITH `matrix_free`, which was written and NOT
    ## exercised: the matrix-free traversal read `self._coeffs` for its
    ## opening seed, and on this path no step has run to set it, so the
    ## pair raised `AttributeError`.  Two flags that are each tested alone
    ## and never together is how that survives.
    for mf in (False, True):
        pss = PSS(cir, method='trap', reltol=1e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=per / 100, maxiterations=40,
                            x0_unknown=True, matrix_free=mf)
        assert pss.converged, 'x0_unknown + matrix_free=%s did not converge' % mf
        v = np.asarray(res['tpss'].v('n2', gnd), dtype=float).ravel()
        peak = 0.5 * (v.max() - v.min())
        assert abs(peak - got[('trap', 100, True)]) < 1e-6, \
            'matrix_free=%s changed the x0_unknown answer: %.6f against ' \
            '%.6f' % (mf, peak, got[('trap', 100, True)])


def test_the_period_column_agrees_with_the_vector_field_at_T():
    """Aprille & Trick's closed form, used as an independent check.

    [AT-O] (IEEE TCT 19(4) 354-360, 1972) eq.(19) gives the period column in
    closed form: `dH/dT = -f(x(T))`, the vector field at the END of the
    period -- so `dx_end/dT = xdot(T)`.

    ⚠ THAT IS THE CONTINUOUS DERIVATIVE AND IT IS NOT WHAT THIS CODE
    COMPUTES, deliberately.  PSS rebuilds the grid at the current `T`, so
    every step scales and the exact derivative of the DISCRETE map carries
    the `dh/dT = h/T` terms `companion_dT` accumulates.  The two agree only
    to O(h), and the accumulated one is the more exact object -- replacing
    it with `-f(x(T))` would be a step backwards.

    ⚠ WHAT MAKES IT VALUABLE IS THE INDEPENDENCE.  The check does not
    re-derive the accumulation; it compares against a quantity taken from
    the converged waveform itself, so an error IN the accumulation cannot
    hide in it.  The assertion is therefore about the RATE, not the value:
    the gap must be O(h) and must FALL under refinement.

    Measured (autonomous phase circuit, relative gap, ratio between
    successive doublings):

        npts    trap                gear
         100    9.51e-02            3.25e-02
         200    4.73e-02  (2.01x)   1.61e-02  (2.02x)
         400    2.36e-02  (2.00x)   7.99e-03  (2.01x)
         800    1.18e-02  (2.00x)   3.99e-03  (2.01x)

    ⚠ AND IT REJECTS THE DEFECT IT WAS ADDED FOR.  With Gear-2's period
    column back on the PARTIAL `d/dh_n` -- the 3/2 error fixed earlier --
    the gap stops falling and plateaus at exactly 0.5, which is
    `|1.5x - x| / |x|`: 0.486 / 0.495 / 0.498 / 0.499, ratio 1.00.  A
    constant-factor error is invisible to any test that only asks whether a
    number is small; it is unmissable to one that asks whether it falls.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'gear'):
        gaps = []
        for npts in (100, 200, 400):
            pss = PSS(_phase_circuit(), method=method, reltol=1e-10)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = pss.solve(period=1e-3, timestep=1e-3 / npts,
                                maxiterations=30)
            assert pss.converged
            T = pss.period
            Xw = np.asarray(res['tpss'].x, dtype=float)
            ir = pss.irefnode
            red = lambda col: np.concatenate((Xw[:ir, col], Xw[ir + 1:, col]))
            x0 = red(0)
            times, hs = pss._period_grid(T, npts, None)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                if pss._solves_history():
                    _a, _b, _c, _d, Mt, _e = pss._traverse_solved_history(
                        x0, red(-2), times, hs, T=T, want_dT=True)
                else:
                    _x0, _xe, _M, Mt = pss._traverse(x0, T, times, hs,
                                                     want_dT=True)
            ## xdot(T) from the converged waveform -- nothing borrowed from
            ## the accumulation being checked
            xdot = (red(-1) - red(-2)) / hs[-1]
            Mt = np.asarray(Mt, dtype=float).ravel()
            gaps.append(float(np.linalg.norm(Mt - xdot)
                              / max(np.linalg.norm(xdot), 1e-300)))

        assert gaps[0] < 0.25, \
            '%s: the period column is %.3e from the vector field at T on ' \
            'the coarsest grid -- an O(1) gap is a wrong column, not a ' \
            'discretisation difference' % (method, gaps[0])
        for a, b in zip(gaps, gaps[1:]):
            assert b < a / 1.6, \
                '%s: the gap to -f(x(T)) went %.3e -> %.3e, a factor of ' \
                '%.2f where O(h) needs ~2. A gap that does not FALL is a ' \
                'constant-factor error in the column -- exactly what a 3/2 ' \
                'partial-derivative bug looks like here' % (method, a, b, a / b)


def test_the_monodromy_transpose_is_a_reverse_replay():
    """`M^T v` from the stored factors, backwards -- no reverse integrator.

    A shooting monodromy is a product of per-step solves, so its transpose
    is that product replayed in REVERSE ORDER with each solve transposed.
    `_traverse_factored` already stores every step's factorisation and every
    factorisation here can solve transposed, so the reverse pass needs no
    new integrator, no refactorisation and no second traversal.

    ⚠ WHY THAT IS WORTH A TEST FOR MACHINERY NOTHING YET USES.  Demir &
    Roychowdhury (TCAD 22(2) 188-196) treat reverse integration as the
    barrier to computing a PPV -- "often unavailable even in existing
    time-domain simulators and may require significant changes to core
    simulation routines" -- and adjoint noise needs the same object.  That
    barrier is real for a forward-only DENSE implementation.  It is not real
    here, and this pins the reason rather than leaving it in a scratch file:
    both remaining capabilities bottleneck on this one piece.

    Measured when it was written: agreement 1.8e-15 with the dense `M^T`,
    and the reverse pass costs 0.75x the forward one -- CHEAPER, because it
    does two `C^T` products against one shared transposed solve where the
    forward does two `C` products and a solve.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_rc_ladder(6), method='gear', reltol=1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 50, maxiterations=2)
    m = pss.cir.n - 1
    times, hs = pss._period_grid(1e-3, 50, None)
    z = np.zeros(m)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        C0, steps, _xl, _xp = pss._traverse_factored(z, z, times, hs)

    ## dense M from the FORWARD matvec, itself pinned against the traversal
    ## by test_the_matrix_free_matvec_is_the_dense_monodromy
    M = np.zeros((2 * m, 2 * m))
    for i in range(2 * m):
        e = np.zeros(2 * m)
        e[i] = 1.0
        M[:, i] = pss._monodromy_matvec(C0, steps, e)

    rng = np.random.default_rng(5)
    for _ in range(4):
        v = rng.standard_normal(2 * m)
        got = pss._monodromy_matvec_transposed(C0, steps, v)
        want = M.T @ v
        err = np.linalg.norm(got - want) / np.linalg.norm(want)
        assert err < 1e-11, \
            'the reverse replay differs from the dense M^T by %.3e' % err

    ## ⚠ and the identity that makes it useful: <M^T u, w> == <u, M w>.
    ## A transpose that is only checked against a matrix built from the
    ## forward pass could share an error with it; this cannot.
    u, w = rng.standard_normal(2 * m), rng.standard_normal(2 * m)
    lhs = float(np.dot(pss._monodromy_matvec_transposed(C0, steps, u), w))
    rhs = float(np.dot(u, pss._monodromy_matvec(C0, steps, w)))
    assert abs(lhs - rhs) < 1e-9 * max(abs(rhs), 1.0), \
        'adjoint identity fails: <M^T u, w> = %.12g against <u, M w> = ' \
        '%.12g' % (lhs, rhs)


def test_every_flag_combination_is_walked():
    """method x `x0_unknown` x `matrix_free` x driven/autonomous, all 24.

    ⚠ WRITTEN BECAUSE A PAIR THAT EACH WORKED ALONE CRASHED TOGETHER.
    `solve(x0_unknown=True, matrix_free=True)` raised `AttributeError` on a
    line written in the same commit as the feature: the matrix-free
    traversal read `self._coeffs` for its opening seed, and on the
    `open_at_x0` path no step has run to set it.  Each flag had a test; the
    PAIR had none, and threading a flag through a second path "so they
    cannot disagree" is exactly when the untested combination looks safest.

    Three things are asserted, and the second is the one with teeth:

      * nothing CRASHES -- an unexpected exception type is a failure, while
        a documented `NotImplementedError` is a result;
      * `matrix_free` does not change the ANSWER.  It is an implementation
        of the same solve, so any difference beyond the convergence
        tolerance is a defect in the matvec, not a numerical detail;
      * the REFUSALS are exactly the expected set -- `gear` with
        `x0_unknown`, which has nothing to change.  A refusal that quietly
        spreads to another combination would look like a passing test.
    """
    import warnings
    import itertools
    circuit.default_toolkit = circuit.numeric

    driven, per_d = _resonator_at_resonance(), None
    per_d = 1.0 / (1.0 / (2 * np.pi * np.sqrt(1e-3 * 1e-9)))
    systems = (('driven', lambda: _resonator_at_resonance()[0],
                _resonator_at_resonance()[1]),
               ('autonomous', _phase_circuit, 1e-3))

    refused, answers = set(), {}
    for sysname, cf, per in systems:
        for method, x0u, mf in itertools.product(
                ('euler', 'trap', 'gear'), (False, True), (False, True)):
            key = (sysname, method, x0u, mf)
            pss = PSS(cf(), method=method, reltol=1e-8)
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    pss.solve(period=per, timestep=per / 100,
                              maxiterations=30, x0_unknown=x0u,
                              matrix_free=mf)
            except NotImplementedError:
                refused.add(key)
                continue
            except Exception as exc:                  # noqa: BLE001
                raise AssertionError(
                    '%r crashed with %s: %s -- a combination that raises '
                    'anything but NotImplementedError is a defect, not a '
                    'documented limit' % (key, type(exc).__name__, exc))
            assert pss.converged, '%r did not converge' % (key,)
            answers[key] = float(pss.period)

    ## `matrix_free` is an implementation of the same solve
    for (sysname, method, x0u, mf), per in list(answers.items()):
        if mf:
            continue
        other = (sysname, method, x0u, True)
        if other in answers:
            assert abs(answers[other] - per) < 1e-9 * max(abs(per), 1e-30), \
                'matrix_free changed the answer for %r: %.12g against ' \
                '%.12g' % ((sysname, method, x0u), answers[other], per)

    expected = {(s, 'gear', True, mf)
                for s in ('driven', 'autonomous') for mf in (False, True)}
    assert refused == expected, \
        'the refused set moved: %r were refused, expected %r' \
        % (sorted(refused), sorted(expected))


def test_the_outer_newton_is_damped_and_the_damping_is_nearly_free():
    """A line search on the shooting Newton, and what it does and does not buy.

    All three `fsolve` call sites took the FULL step with `limiter=None` and
    no line search.  Brachtendorf et al. (TCAD 33(6) 867-878) describe
    "shooting, finite difference, or harmonic balance techniques IN
    CONJUNCTION WITH A DAMPED NEWTON METHOD" as what is "widely employed"
    for limit cycles, so the absence was a departure from standard practice
    rather than a neutral choice.

    ⚠ WHAT IT DOES NOT BUY, measured before it was kept: it does NOT rescue
    a far seed.  Van der Pol at mu=1 from seeds at 4x/10x/30x the orbit
    amplitude still fails -- the failure MODE changes (10x and 30x now raise
    loudly where they used to return a wrong period) but seed dependence is
    untouched.  That matches the literature: the higher the Q, "the tighter
    are the constraints for the initial estimate", and even the probe
    technique "is still not always obtained".  Damping is the baseline, not
    the fix for the trivial-root basin.

    ⚠ AND THE COST HAD TO BE ENGINEERED AWAY.  The obvious implementation
    evaluates `F` at the trial point and lets the NEXT iteration evaluate it
    again at the same point, which DOUBLES every converging solve -- measured
    10 -> 20 residual evaluations here, each a full pass over the period.
    Carrying the trial evaluation forward makes the common case cost what
    the undamped loop cost: 10 -> 11 and 2 -> 3, the +1 being the last
    iteration's trial that the loop exits before reusing.
    """
    import warnings
    import numpy as _np
    from pycircuit.circuit import analysis as _an
    from pycircuit.circuit import numeric as _tk
    circuit.default_toolkit = circuit.numeric

    ## the damping itself, on the textbook overshoot problem
    calls = [0]

    def f(x):
        calls[0] += 1
        v = float(_np.clip(x[0], -1e8, 1e8))
        return (_np.array([_np.arctan(v)]),
                _np.array([[1.0 / (1.0 + v * v)]]))

    _x, _i, ier_off, _m = _an.fsolve(f, _np.array([1.5]), maxiter=40,
                                     toolkit=_tk, full_output=True,
                                     line_search=False)
    assert ier_off == 2, \
        'undamped Newton converged on arctan from 1.5, so this problem no ' \
        'longer separates the two and the test proves nothing'
    calls[0] = 0
    x_on, _i, ier_on, _m = _an.fsolve(f, _np.array([1.5]), maxiter=40,
                                      toolkit=_tk, full_output=True,
                                      line_search=True)
    assert ier_on == 1 and abs(float(x_on[0])) < 1e-8, \
        'the damped solve did not reach the root (ier=%r, x=%r)' \
        % (ier_on, x_on)

    ## and the cost, on a real shooting solve: within one evaluation of
    ## undamped, not double it
    cir, per = _resonator_at_resonance()
    counts = {}
    for tag, ls in (('undamped', False), ('damped', True)):
        n = [0]
        orig = _an.fsolve

        def counting(fn, *a, **k):
            def g(z, *aa):
                n[0] += 1
                return fn(z, *aa)
            k['line_search'] = ls
            return orig(g, *a, **k)

        import pycircuit.circuit.shooting as _sh
        _sh.analysis.fsolve = counting
        try:
            pss = PSS(cir, method='trap', reltol=1e-10)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=per, timestep=per / 200, maxiterations=40)
            assert pss.converged
            counts[tag] = n[0]
        finally:
            _sh.analysis.fsolve = orig

    assert counts['damped'] <= counts['undamped'] + 2, \
        'damping cost %d residual evaluations against %d undamped -- the ' \
        'trial evaluation is not being carried forward, which doubles ' \
        'every converging solve' % (counts['damped'], counts['undamped'])

    ## ⚠ AND THAT PSS ACTUALLY ASKS FOR IT.  The two checks above force the
    ## flag themselves, so both pass with every shipped call site setting
    ## `line_search=False` -- verified by mutation, which is how this gap
    ## was found.  This one RECORDS what PSS passes instead of overriding
    ## it.
    seen = []
    orig = _an.fsolve

    def recording(fn, *a, **k):
        seen.append(bool(k.get('line_search', False)))
        return orig(fn, *a, **k)

    import pycircuit.circuit.shooting as _sh
    _sh.analysis.fsolve = recording
    try:
        for cf, per_, method in ((lambda: _resonator_at_resonance()[0],
                                  _resonator_at_resonance()[1], 'trap'),
                                 (_phase_circuit, 1e-3, 'trap')):
            pss = PSS(cf(), method=method, reltol=1e-8)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=per_, timestep=per_ / 100, maxiterations=30)
    finally:
        _sh.analysis.fsolve = orig

    assert seen and all(seen),         'PSS called fsolve with line_search=%r -- the damping is '         'implemented and not asked for' % (seen,)


def test_the_monodromy_is_correct_across_a_switching_boundary():
    """The saltation concern, run as a falsifier and NOT confirmed.

    Bizzarri & Wei (ECCTD 2011) state that for hybrid systems "the monodromy
    matrix is not defined at impact events", so the transition matrix needs
    a SALTATION correction `S` -- and with no state jump `S` still differs
    from the identity whenever the VECTOR FIELD jumps, which an ideal switch
    does.  If that reached this code, every switching circuit's monodromy
    would be accumulated through a boundary where the linearisation does not
    exist.

    ⚠ IT DOES NOT REACH IT, and the reason is structural rather than lucky.
    PSS's monodromy is the exact derivative of the DISCRETE period map: each
    step uses its own converged `Jf` and `C`, which already describe
    whichever side of the switch that step is on.  The saltation matrix is a
    CONTINUOUS-time construct for correcting `Phi(t2,t0)` when the flow is
    discontinuous; a discrete map has no instant at which the field is
    undefined.

    ⚠ THE ASSERTION IS ON THE RATE, not the size.  A missing saltation term
    is an O(1) structural error that does not vanish under refinement, so
    "the gap is small" would not distinguish it from discretisation error.
    Measured on a switched RC whose control crosses twice per period, with
    Ron/Roff spanning six decades:

        npts    rel err (analytic monodromy vs finite differences)
         400    4.858e-04
         800    2.426e-04   (2.00x)
        1600    1.212e-04   (2.00x)
        3200    6.058e-05   (2.00x)

    Exactly first order, so the residue is discretisation and there is no
    O(1) term hiding under it.

    ⚠ AND THE CIRCUIT HAD TO BE BUILT TWICE.  The first two attempts put the
    switch in series with the source, which CLAMPS the node each period and
    erases the state: |M| came out 1.4e-22 and then 7.4e-51, so the test was
    comparing numerical zero against numerical zero and reporting a
    meaningless "rel err 1.000".  The switched element must change the decay
    RATE without clamping, or there is no monodromy to check.
    """
    import warnings
    from pycircuit.circuit import VSwitch
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def switched(ron, roff):
        c = SubCircuit()
        c.add_node('ctl')
        c.add_node('a')
        c.add_node('b')
        c['vc'] = VSin('ctl', gnd, va=2.0, freq=1.0 / per)
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per, phase=90.0)
        c['rs'] = R('a', 'b', r=1e5)
        c['c'] = C('b', gnd, c=1e-6)
        ## a switched LOAD, not a switched path to the source
        c['sw'] = VSwitch('b', gnd, 'ctl', gnd, Ron=ron, Roff=roff,
                          Von=1.0, Voff=0.0)
        return c

    gaps = []
    for npts in (200, 400, 800):
        pss = PSS(switched(1e3, 1e9), method='trap', reltol=1e-11)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(per, npts, None)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=per / npts, maxiterations=40)
        assert pss.converged
        ir = pss.irefnode
        Xw = np.asarray(res['tpss'].x, dtype=float)
        ctl = Xw[pss.cir.get_node_index('ctl')]
        toggles = int(np.sum(np.diff((ctl > 0.5).astype(int)) != 0))
        assert toggles >= 2, \
            'the switch no longer toggles (%d crossings), so this test has ' \
            'lost its subject' % toggles
        x0 = np.concatenate((Xw[:ir, 0], Xw[ir + 1:, 0]))

        def phi(v):
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                _a, xe, _b, _c = pss._traverse(np.asarray(v, dtype=float),
                                               per, times, hs, want_dT=False)
            return np.asarray(xe, dtype=float)

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _a, _e, M, _c = pss._traverse(x0, per, times, hs, want_dT=False)
        M = np.asarray(M, dtype=float)
        assert np.linalg.norm(M) > 0.1, \
            'the monodromy is %.3e -- the circuit erases its state each ' \
            'period, so there is nothing to check' % np.linalg.norm(M)

        base = phi(x0)
        Mfd = np.zeros((m, m))
        for j in range(m):
            d = 1e-7 * max(abs(x0[j]), 1.0)
            xp = x0.copy()
            xp[j] += d
            Mfd[:, j] = (phi(xp) - base) / d
        gaps.append(float(np.linalg.norm(M - Mfd)
                          / max(np.linalg.norm(Mfd), 1e-300)))

    for a, b in zip(gaps, gaps[1:]):
        assert b < a / 1.6, \
            'the monodromy gap across the switch went %.3e -> %.3e, a ' \
            'factor of %.2f where O(h) needs ~2. A gap that does NOT fall ' \
            'is the O(1) signature of a missing saltation correction' \
            % (a, b, a / b)


def test_tstab_rescues_the_trivial_root_basin():
    """`tstab`: pre-integrate a transient, then shoot from where it lands.

    The stabilisation time every commercial PSS offers, and the standard
    answer to a seed that is not close enough. It is the remedy for the
    failure this analysis fails most often — an unseeded autonomous run
    starts at the operating point, which sits at the bottom of the
    trivial-root basin.

    Measured on van der Pol seeded near the unstable DC point:

        circuit                        without tstab      periods needed
        mu = 1  (strongly attracting)  LinAlgError              1
        mu = 0.05 (high-Q)             not converged           ~24

    The count is the `1/mu` amplitude-envelope constant — a property of how
    strongly the cycle attracts, not of how bad the seed is. From 4x and
    even 20x the orbit amplitude, one period suffices at mu = 1.

    ⚠ THE STOPPING POINT IS THE CALLER'S, DELIBERATELY. De Luca et al. give
    a criterion for detecting the handoff; the probe it rests on was
    measured here and does NOT identify it — near the DC point the monodromy
    is nearly constant, so the probe settles into its own eigenvector and
    reports convergence while the state is stuck at the trivial root. Every
    obvious substitute shares that defect, because the trivial root IS a
    fixed point of the period map and passes every periodicity test.

    ⚠ AND `tstab` MUST OUTRANK THE OPERATING-POINT SEED. The autonomous path
    seeds from DC when `x0 is None`; a pre-integration run precisely to
    leave that basin must not then be replaced by the point at the bottom of
    it. This asserts that ordering, because getting it backwards would
    disable the option exactly where it earns its place.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    def vdp(mu):
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: mu * (u - u ** 3 / 3.0))
        return c

    per = 6.663293
    seed = np.zeros(vdp(1.0).cir.n - 1 if hasattr(vdp(1.0), 'cir')
                    else PSS(vdp(1.0)).cir.n - 1)
    seed[0] = 0.04                      # near the unstable DC point

    def run(mu, period, tstab, x0=None):
        pss = PSS(vdp(mu), method='trap', reltol=1e-9)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                pss.solve(period=period, timestep=period / 200, x0=x0,
                          maxiterations=40, tstab=tstab)
            return pss, ('converged' if pss.converged else 'not-converged')
        except Exception as exc:                          # noqa: BLE001
            return pss, type(exc).__name__

    ## without it, the strongly-attracting case cannot be solved at all
    _p, cold = run(1.0, per, None, seed)
    assert cold != 'converged', \
        'the cold solve now converges from a DC-adjacent seed, so this test ' \
        'has lost the failure it exists to rescue (got %r)' % cold

    ## one period of stabilisation is enough there
    warm, out = run(1.0, per, 1.0 * per, seed)
    assert out == 'converged', \
        'tstab of one period did not rescue the DC-adjacent seed: %r' % out
    assert abs(warm.period - per) < 5e-3 * per, \
        'tstab converged to T=%.6f against a measured %.6f' \
        % (warm.period, per)

    ## and the state actually moved -- tstab ran, rather than being ignored
    assert hasattr(warm, 'tstab_state'), 'tstab did not record its state'
    assert np.linalg.norm(np.asarray(warm.tstab_state) - seed) > 0.1, \
        'the pre-integration returned essentially the seed it was given'

    ## ⚠ AND THE LIMIT IT CANNOT PASS, asserted rather than left to be
    ## rediscovered.  With `x0=None` the seed is the OPERATING POINT, and on
    ## an autonomous circuit that is an equilibrium -- a transient started
    ## exactly there never leaves, so no amount of `tstab` escapes the basin.
    ## The pre-integration needs somewhere to go: an `x0` off the
    ## equilibrium, or a device `ic`.  This is the same reason the option is
    ## not a substitute for the probe technique.
    _auto, out = run(1.0, per, 3.0 * per, None)
    assert out != 'converged', \
        'tstab from the operating point now escapes an equilibrium (%r) -- ' \
        'if that is real the docstring limit is wrong and should be fixed, ' \
        'not the test' % out


def test_tstab_also_runs_on_the_driven_path():
    """`tstab` is NOT autonomous-only, and this test exists to pin that.

    The pre-integration sits outside the `if self.autonomous:` branch on
    purpose: a driven circuit gets a warm start too, and that is the class
    every commercial tool applies it to most. Nesting it one level deeper
    would be invisible -- the autonomous test would still pass, driven runs
    would silently ignore `tstab=`, and the option would be half a feature.

    ⚠ It is also the class the AUTOMATIC criterion is for. De Luca, Bolcato
    & Schilders (2019) is titled for NON-autonomous circuits and offers the
    autonomous case only as conditional future work, so when that criterion
    is built it is gated here, not on van der Pol -- a driven circuit has no
    trivial root for the test to be attracted to. See the A4 entry in
    `doc/pss_roadmap_260902.md`; the earlier gate ran on the wrong class.

    The check is that the state actually MOVED and the answer did not: a
    warm start may not change where a converging solve lands.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def run(tstab):
        pss = PSS(_q20_rlc(), method='trap', reltol=1e-9)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pss.solve(period=per, timestep=1e-5, maxiterations=40,
                            tstab=tstab)
        assert pss.converged, 'tstab=%r did not converge' % (tstab,)
        return pss, float(np.max(np.abs(np.asarray(
            res['tpss'].v('c'), dtype=float).ravel())))

    cold_pss, cold_peak = run(None)
    warm_pss, warm_peak = run(3.0 * per)

    assert getattr(cold_pss, 'tstab_state', None) is None, \
        'tstab_state was recorded on a run that asked for no tstab'
    warm_seed = getattr(warm_pss, 'tstab_state', None)
    assert warm_seed is not None, \
        'tstab= was accepted on the driven path but no pre-integration ran ' \
        '-- the block is probably nested inside `if self.autonomous:`'
    assert float(np.max(np.abs(np.asarray(warm_seed, dtype=float)))) > 1e-9, \
        'the driven pre-integration returned the zero seed unchanged'

    ## the warm start moves the SEED, never the SOLUTION
    assert abs(warm_peak - cold_peak) < 1e-3, \
        'tstab changed a converged answer: %.6f warm against %.6f cold' \
        % (warm_peak, cold_peak)


def _pac_L_and_B(pss, N):
    """`L` and `B` exactly as the withdrawn `PAC.solve` body builds them.

    Kept as a helper rather than inlined because the point of the test is
    that THIS construction -- the one in the tree -- has the properties
    DAC'96 claims for it, and only for the method it was written for.
    """
    times = pss.times[:-1]
    hs = np.diff(pss.times)
    M = len(times)
    L = np.zeros((N * M, N * M))
    B = np.zeros_like(L)
    for i, (_t, h, J, Cm) in enumerate(zip(times, hs, pss.Jtvec, pss.Cvec)):
        L[i * N:(i + 1) * N, i * N:(i + 1) * N] = J
        if i > 0:
            L[i * N:(i + 1) * N, (i - 1) * N:i * N] = -np.asarray(Cm) / h
    B[0:N, (M - 1) * N:M * N] = -np.asarray(pss.Cvec[-1]) / hs[0]
    return L, B, M


def _pac_operator_pieces(method, npts=40):
    """One trajectory, two descriptions of it: `L`/`B`, and our own matvec.

    Both sides are driven from the SAME seed and the SAME grid on purpose.
    Convergence is irrelevant here — these are linear-algebra identities on
    whatever `Jtvec`/`Cvec` hold — but the two sides must describe one
    trajectory or the comparison means nothing.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _q20_rlc()
    per, N = 1e-3, cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-9)
    pss._open_at_x0 = False
    pss.autonomous = False
    times, hs = pss._period_grid(per, npts, None)
    x0 = np.zeros(N)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss._traverse(x0, per, times, hs, False)
        Jt = [np.asarray(j).copy() for j in pss.Jtvec]
        Cv = [np.asarray(c).copy() for c in pss.Cvec]
        tms = np.asarray(pss.times).copy()
        opening, steps, _x0, _x, _ = pss._traverse_factored_plain(
            x0, per, times, hs)
    pss.Jtvec, pss.Cvec, pss.times = Jt, Cv, tms
    L, B, M = _pac_L_and_B(pss, N)
    Mmv = np.column_stack([pss._monodromy_matvec_plain(opening, steps, e)
                           for e in np.eye(N)])
    ## the last block of `L^-1 B`, as a map on the wrapped state
    H = np.zeros((N, N))
    for j in range(N):
        e = np.zeros(N)
        e[j] = 1.0
        b = np.zeros(N * M)
        b[0:N] = B[0:N, (M - 1) * N:M * N] @ e
        H[:, j] = np.linalg.solve(L, b)[(M - 1) * N:M * N]
    return pss, L, B, N, M, H, Mmv


def test_the_pac_operator_is_the_monodromy_and_L_is_never_formed():
    """⚠ THE WITHDRAWN PAC BODY HOLDS THE RIGHT OPERATOR, and this pins it.

    DAC'96 reaches the iterative form by "reinterpreting the use of `L^-1`
    … as a preconditioner". The algebra is one line:

        (L + αB)v = -u   ⟺   (I + α L^-1 B) v = -L^-1 u

    `L` is block lower bidiagonal, so applying `L^-1` is forward
    substitution through the timesteps — which is the recursion
    `_monodromy_matvec_plain` already runs against stored factors. So
    `H := L^-1 B` IS THE MONODROMY, the 419.5 GiB that withdrew PAC was
    entirely the cost of FORMING `L`, and the rewrite keeps the operator
    and never builds the matrix.

    This test exists because that identification is the load-bearing claim
    under the whole PAC item, and it arrived by relay from a reading of the
    paper. A claim that is not in the suite is a claim that drifts.
    """
    _pss, L, B, N, M, H, Mmv = _pac_operator_pieces('euler')

    ## (1) the two structural claims DAC'96 makes -- about OUR matrices
    upper = 0.0
    for i in range(M):
        for j in range(M):
            if j > i or j < i - 1:
                upper = max(upper, np.abs(
                    L[i * N:(i + 1) * N, j * N:(j + 1) * N]).max())
    assert upper == 0.0, \
        'L is not block lower bidiagonal (max %g outside) -- the forward ' \
        'substitution that makes L^-1 cheap does not apply' % upper
    Bout = B.copy()
    Bout[0:N, (M - 1) * N:M * N] = 0.0
    assert np.abs(Bout).max() == 0.0, \
        'B is not confined to the first N rows and last N columns; the ' \
        'periodic wrap is the only thing it is supposed to carry'

    ## (2) the algebraic identity the whole reformulation rests on
    rng = np.random.default_rng(0)
    v = rng.standard_normal(N * M)
    alpha = np.exp(-2j * np.pi * 3.0)
    lhs = (L + alpha * B) @ v
    rhs = L @ (v + alpha * np.linalg.solve(L, B @ v))
    assert np.linalg.norm(lhs - rhs) / np.linalg.norm(lhs) < 1e-12, \
        '(L + aB)v != L(I + a L^-1 B)v -- the preconditioned form is not ' \
        'the same operator'

    ## (3) and the identification itself: L^-1 B IS the monodromy, up to
    ##     the sign B carries (`B = -C/h`)
    rel = np.linalg.norm(H + Mmv) / np.linalg.norm(Mmv)
    assert rel < 1e-11, \
        'the last block of L^-1 B is not -M (rel %.3g) -- PAC cannot be ' \
        'built on the traversal if the two disagree' % rel


@pytest.mark.parametrize('method', ['trap', 'gear'])
def test_the_pac_L_is_backward_euler_only(method):
    """⚠ AND THE OPERATOR IS EULER-SHAPED, WHICH IS NOT A DETAIL.

    The withdrawn body says so in its own comment — "create LHS matrix
    using backward Euler discretization" — and builds `L` with exactly two
    terms per row: `L[i,i] = J`, `L[i,i-1] = -C/h`. A TWO-STEP method's
    variational system is block TRI-diagonal, so for `trap` or `gear` that
    `L` describes a different recursion than the trajectory it was built
    from, and `L^-1 B` is not the monodromy at all.

    ⚠ THIS IS THE TRAP A REWRITE WOULD FALL INTO. The structural checks in
    the sibling test PASS for every method — `L` really is block lower
    bidiagonal whatever `Jtvec` holds, because that is a property of how
    the loop writes it, not of the physics. Verifying the structure and
    then assuming the identification is exactly the mistake: measured on
    the Q=20 resonator at 40 points, our monodromy has ρ = 0.8545 (trap) /
    0.8412 (gear) against the analytic 0.854636, while `-L^-1 B` has ρ = 0.

    So the rewrite must take `H` from the traversal, which carries each
    step's own `(alphas, b)`, and must NOT rebuild `L`. That is the same
    conclusion the memory cost forces, reached independently.
    """
    _pss, _L, _B, _N, _M, H, Mmv = _pac_operator_pieces(method)

    rho_ours = max(abs(np.linalg.eigvals(Mmv)))
    rho_H = max(abs(np.linalg.eigvals(-H)))
    assert abs(rho_ours - 0.854636) < 0.02, \
        '%s monodromy rho %.6f is not the analytic exp(-pi/Q) -- this test ' \
        'compares against it, so it has to be right first' % (method, rho_ours)
    assert rho_H < 1e-3, \
        '%s: -L^-1 B now has rho %.6g. If the Euler-shaped L has started ' \
        'agreeing with a two-step monodromy, either L gained the third term ' \
        'or the monodromy lost it -- find out which before trusting either' \
        % (method, rho_H)


def _pac_circuit(f0=1e3, Q=20.0):
    """The Q=20 resonator with an AC amplitude on its source.

    `vac` and `va` are different knobs: `va` drives the LARGE signal that
    the periodic operating point is a response to, `vac` is the small
    signal PAC linearises for. A circuit with only `va` set has nothing
    for PAC to analyse, which `PAC.solve` refuses rather than returning
    zeros.
    """
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * f0) ** 2 * 1e-3)
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=f0)
    c['R'] = R('a', 'b', r=(1.0 / Q) * np.sqrt(L_ / C_))
    c['L'] = L('b', 'c', L=L_)
    c['C'] = C('c', gnd, c=C_)
    return c


def _pac_y0(pss, freq, per):
    """`y_0` by a DENSE solve of the same `m x m` system PAC solves.

    Dense on purpose: this is the reference the matrix-free path is
    measured against, so it must not share its solver.
    """
    fp = pss.factored_period()
    irn = pss.irefnode
    (u_ac,) = remove_row_col((pss.cir.u(0, analysis='ac'),), irn,
                             pss.toolkit)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()
    w, _ = pss._forced_replay(fp, freq, u_ac)
    alpha = np.exp(-2j * np.pi * freq * per)
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    return np.linalg.solve(np.eye(n) - alpha * M, alpha * np.asarray(w)), M


def _pac_vs_ac(method, npts, freq=700.0, per=1e-3, x0_unknown=False):
    import warnings
    from pycircuit.circuit.analysis_ss import AC
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    m = cir.n - 1
    pss = PSS(cir, method=method, reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40,
                  x0_unknown=x0_unknown)
    assert pss.converged, '%s at %d points did not converge' % (method, npts)
    y0, _M = _pac_y0(pss, freq, per)
    irn = pss.irefnode
    xac = np.asarray(AC(cir, toolkit=circuit.numeric).solve(freqs=freq).x,
                     dtype=complex).ravel()
    xac = np.concatenate((xac[:irn], xac[irn + 1:]))
    return (np.linalg.norm(np.asarray(y0)[:m] - xac) / np.linalg.norm(xac),
            pss)


def test_pac_agrees_with_the_ac_analysis_on_a_linear_circuit():
    """⚠ THE GATE FOR PAC, against a reference it cannot influence.

    On a LINEAR circuit the periodic operating point is irrelevant to the
    small signal: the linearisation is constant, every sideband but `l = 0`
    vanishes, and the LPTV response collapses to the LTI transfer function.
    So `AC` — a different analysis, a different code path, an `(sC + G)`
    solve that never touches a monodromy — is the answer PAC must produce.

    This is the check the withdrawn implementation never had: its only test
    was `@unittest.skip('Skip failing test')`.

    ⚠ AND IT IS WHAT CAUGHT THE DEAD BODY'S REAL DEFECT. That body read its
    source vector as `self.cir.u(0, analysis_name)` — POSITIONALLY, into a
    signature whose second parameter is `epar`. It took the transient source
    at `t = 0`, which is zero for every sinusoid, so the whole analysis
    would have returned zeros. Reproduced here before it was fixed: `|PAC|`
    exactly 0 against `|AC|` of 1.01, at every frequency and every method.
    """
    rel, _pss = _pac_vs_ac('gear', 1000)
    assert rel < 1e-4, \
        'PAC disagrees with the AC analysis by %.3e on a LINEAR circuit, ' \
        'where the two are solving the same problem by different routes' % rel


@pytest.mark.parametrize('method,x0_unknown,expect', [
    ('trap', False, 'first'),
    ('trap', True, 'second'),
    ('euler', False, 'first'),
    ('euler', True, 'first'),
])
def test_pac_order_is_lost_to_the_manufacturing_step(method, x0_unknown,
                                                     expect):
    """⚠ PAC ON THE PLAIN PATH IS FIRST ORDER WHATEVER THE METHOD.

    `_traverse_factored_plain` takes one step OUTSIDE its loop to
    manufacture a history, and that step is not in `steps`. For the
    homogeneous map its effect is folded into the `opening` triple — the
    documented flat-history approximation. For the DRIVEN map it also means
    THE SOURCE IS NEVER APPLIED THERE: one step of `u` out of `N`, a
    relative O(h), which drags a second-order method down to first.

    ⚠ THE RATE IS THE ASSERTION, NOT THE SIZE. A constant-factor error is
    invisible to "is it small" and unmissable to "does it fall" — the same
    reason the saltation falsifier and the A&T period column are rate
    checks. Measured against the AC analysis at 700 Hz, per doubling:

        trap, plain            2.00x   4.13e-03 at 250 points
        trap, x0_unknown=True  4.00x   1.09e-04 at 250 points
        euler, either          2.00x   1.40e-02, identical to five digits

    ⚠ EULER IS THE CONTROL AND IT IS WHY THIS IS THE MANUFACTURING STEP.
    `x0_unknown` does not move euler's rate at all — it is a first-order
    method either way — so the trapezoidal gain cannot be some other thing
    the formulation does. The trajectory is not implicated either: trap's
    WAVEFORM converges at ~4.2x per doubling on this circuit with or
    without the manufacturing step.
    """
    r1, _ = _pac_vs_ac(method, 250, x0_unknown=x0_unknown)
    r2, _ = _pac_vs_ac(method, 500, x0_unknown=x0_unknown)
    ratio = r1 / r2
    if expect == 'second':
        assert 3.4 < ratio < 4.6, \
            '%s/x0_unknown=%s: %.4e -> %.4e is %.2fx per doubling, not the ' \
            "method's own second order. If the manufacturing step now " \
            'carries its source this test should be REWRITTEN, not relaxed' \
            % (method, x0_unknown, r1, r2, ratio)
    else:
        assert 1.7 < ratio < 2.4, \
            '%s/x0_unknown=%s: %.4e -> %.4e is %.2fx per doubling, not the ' \
            'first order the dropped manufacturing source predicts' \
            % (method, x0_unknown, r1, r2, ratio)


def test_the_forced_replay_superposes():
    """`y_end = M y0 + w` — the property the whole `m x m` reduction rests on.

    PAC solves an `m x m` system instead of an `(N m) x (N m)` one because
    the driven replay is LINEAR in its initial state and its source
    separately. If that ever stopped holding, `(I - alpha M) y_0 = alpha w`
    would be solving the wrong equation, and the answer would still look
    entirely plausible.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 200, maxiterations=40)
    fp = pss.factored_period()
    irn = pss.irefnode
    (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), irn, pss.toolkit)
    u_ac = np.asarray(u_ac, dtype=complex).ravel()

    rng = np.random.default_rng(3)
    y0 = (rng.standard_normal(fp.width)
          + 1j * rng.standard_normal(fp.width))
    w, _ = pss._forced_replay(fp, 700.0, u_ac)
    both, _ = pss._forced_replay(fp, 700.0, u_ac, y0=y0)
    pred = np.asarray(fp.matvec(y0)) + np.asarray(w)
    rel = np.linalg.norm(np.asarray(both) - pred) / np.linalg.norm(pred)
    assert rel < 1e-12, \
        'the driven replay does not superpose (rel %.3e): y_end != M y0 + w' \
        % rel


def test_pac_recycling_matches_the_per_frequency_solve():
    """ONE Krylov subspace for the whole sweep — Telichevesky et al. Thm 1.

    `A(alpha) = I - alpha M`, so `span{r, Ar, A^2 r, …}` is the Krylov
    space of `M` and does not depend on `alpha` at all. A basis built once
    therefore serves every frequency in the sweep, and each frequency costs
    a small dense least-squares over it rather than its own run of
    full-period replays.

    ⚠ THE RIGHT-HAND SIDE IS NOT SHARED, and that is why the implementation
    checks rather than assumes: `w(f)` genuinely differs per frequency, so
    the shared span is not the space GMRES would have picked for any but
    the first. It minimises the TRUE residual over the span and extends the
    basis — one matvec, kept for every later frequency — until every
    frequency is inside tolerance. The answer is therefore never worse than
    the per-frequency solve; only the count moves. Measured on RC ladders,
    24 frequencies:

        m=4    72 matvecs -> 6    12.0x
        m=14  168 matvecs -> 15   11.2x
        m=42  302 matvecs -> 26   11.6x

    with the two routes agreeing to 1.2e-13.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, vac=1.0, freq=1e3)
    for k in range(6):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-7)
    per = 1e-3
    pss = PSS(c, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 150, maxiterations=40)
    assert pss.converged

    freqs = np.logspace(1, 4, 12)
    out = {}
    for rec in (False, True):
        pac = PAC(c, toolkit=circuit.numeric)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = pac.solve(pss, freqs, recycle=rec)
        out[rec] = (np.asarray(res.x, dtype=complex), pac.matvecs)

    a, na = out[False]
    b, nb = out[True]
    rel = np.linalg.norm(a - b) / np.linalg.norm(a)
    assert rel < 1e-9, \
        'recycled sweep disagrees with the per-frequency solve by %.3e -- ' \
        'the shared subspace is supposed to change the COST and not the ' \
        'answer' % rel
    assert nb < na, \
        'recycling used %d matvecs against %d for the per-frequency solve, ' \
        'so it is not recycling anything' % (nb, na)


def test_pac_refuses_what_it_cannot_answer():
    """The two ways to ask PAC a question that has no answer.

    Both were reachable in the withdrawn body, and neither announced
    itself: a source with no `vac` gave a zero right-hand side and returned
    zeros, and nothing checked that the operating point had converged.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    ## (a) no operating point yet
    cir = _pac_circuit()
    pss = PSS(cir, method='gear', reltol=1e-10)
    with pytest.raises(RuntimeError, match='call solve'):
        pss.factored_period()

    ## (b) an operating point, but no small-signal source
    L_, C_ = 1e-3, 1.0 / ((2 * np.pi * 1e3) ** 2 * 1e-3)
    quiet = SubCircuit()
    quiet.add_node('a'); quiet.add_node('b')
    quiet['vs'] = VSin('a', gnd, va=1.0, vac=0.0, freq=1e3)
    quiet['R'] = R('a', 'b', r=(1.0 / 20.0) * np.sqrt(L_ / C_))
    quiet['L'] = L('b', 'c', L=L_)
    quiet['C'] = C('c', gnd, c=C_)
    pss2 = PSS(quiet, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss2.solve(period=per, timestep=per / 150, maxiterations=40)
    assert pss2.converged
    with pytest.raises(ValueError, match='identically zero'):
        PAC(quiet, toolkit=circuit.numeric).solve(pss2, [700.0])


def test_the_matvecs_take_a_complex_vector():
    """PAC needs `I + alpha(f) H` with `alpha` complex; `M` is real.

    So a complex product is TWO REAL REPLAYS against the same stored
    factors, exactly — not a complex refactorisation, which would double
    the stored factors for a map with no imaginary part. The three matvecs
    used to cast with `dtype=float`, which does not refuse a complex vector,
    it DISCARDS its imaginary half.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 150, maxiterations=40)
    fp = pss.factored_period()
    rng = np.random.default_rng(7)
    vr = rng.standard_normal(fp.width)
    vi = rng.standard_normal(fp.width)

    for name, mv in (('forward', fp.matvec),
                     ('transposed', fp.matvec_transposed)):
        got = mv(vr + 1j * vi)
        want = np.asarray(mv(vr)) + 1j * np.asarray(mv(vi))
        assert np.iscomplexobj(got), \
            '%s matvec returned a real result for a complex vector -- the ' \
            'imaginary half was discarded' % name
        assert np.array_equal(got, want), \
            '%s matvec is not exactly two real replays' % name
    assert not np.iscomplexobj(np.asarray(fp.matvec(vr))), \
        'the real path became complex; a real replay must stay real'


def _adjoint_ladder(sections=8):
    c = SubCircuit()
    c.add_node('n0')
    c['vs'] = VSin('n0', gnd, va=1.0, vac=1.0, freq=1e3)
    for k in range(sections):
        a, b = 'n%d' % k, 'n%d' % (k + 1)
        c.add_node(b)
        c['R%d' % k] = R(a, b, r=1e3)
        c['C%d' % k] = C(b, gnd, c=1e-7)
    return c


def test_the_adjoint_row_is_m_forward_solves_in_one():
    """⚠ THE MANY-TO-ONE IDENTITY, which is what makes pnoise affordable.

    pnoise is hundreds of sources into one output. Forward, that is one
    solve PER SOURCE, and recycling does not help because the right-hand
    side is what changes. Okumura et al. (1993) reach for the adjoint for
    exactly this reason -- "it is efficient to use the adjoint method …
    because circuits have many noise sources" -- and

        d^T y_0 = alpha ((I - alpha M)^-T d)^T W u

    turns the whole row into ONE transposed solve plus a reverse replay.

    ⚠ THE ASSERTION IS THE COUNT, NOT THE CLOCK. Measured speedup grows
    linearly with `m` (7.0x / 16.9x / 40.0x at m = 6 / 14 / 32, agreement
    ~1e-15), which is the right shape -- but this machine runs more than
    one agent and a wall-clock ratio here would be measuring the neighbours.
    So the test compares matvec counts, which no concurrent load can move.

    ⚠ AND IT COSTS NO NEW MACHINERY. `_traverse_factored` already stores
    every step's factorisation and every factorisation already solves
    transposed, so the reverse pass needs no reverse integrator -- the
    thing Demir & Roychowdhury call "often unavailable even in existing
    time-domain simulators".
    """
    import warnings
    import scipy.sparse.linalg as spla
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder()
    per, f = 1e-3, 700.0
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 120, maxiterations=40)
    assert pss.converged
    fp = pss.factored_period()
    n = fp.width
    alpha = np.exp(-2j * np.pi * f * per)

    pac = PAC(cir, toolkit=circuit.numeric)
    row = pac.adjoint_transfer_row(pss, f, 1)
    adjoint_matvecs = pac.matvecs

    ## the forward route: one solve per source
    d = np.zeros(n)
    d[1] = 1.0
    fwd = np.zeros(m, dtype=complex)
    fwd_matvecs = [0]

    def _mv(v):
        fwd_matvecs[0] += 1
        return np.asarray(v) - alpha * fp.matvec(v)

    for i in range(m):
        u = np.zeros(m)
        u[i] = 1.0
        w, _ = pss._forced_replay(fp, f, u)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        y0, info = spla.gmres(A, alpha * np.asarray(w), rtol=1e-13,
                              restart=min(n, 50), maxiter=min(n, 200))
        assert info == 0
        fwd[i] = d @ y0

    rel = np.linalg.norm(row - fwd) / np.linalg.norm(fwd)
    assert rel < 1e-11, \
        'the adjoint row disagrees with %d forward solves by %.3e -- one ' \
        'of the two transposes is wrong, and the forward route is the one ' \
        'with an independent reference behind it' % (m, rel)
    assert adjoint_matvecs < fwd_matvecs[0], \
        'the adjoint used %d matvecs against %d for the forward route at ' \
        'm=%d; the whole point is that it does not scale with the number ' \
        'of sources' % (adjoint_matvecs, fwd_matvecs[0], m)


def test_the_adjoint_row_refuses_the_plain_path():
    """It needs the reverse replay, which is Gear-2 only.

    Refused in words rather than reaching for `alphas[2]` on a one-step
    companion, which is what it used to do one frame deeper.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    per = 1e-3
    pss = PSS(cir, method='trap', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 120, maxiterations=40)
    assert pss.converged
    with pytest.raises(NotImplementedError, match='solved-history'):
        PAC(cir, toolkit=circuit.numeric).adjoint_transfer_row(pss, 700.0, 1)


def _sideband_forward(pss, fp, freq, k, l, N, alpha, A):
    """`H_l` for every source, the expensive way: one driven solve each."""
    m = pss.cir.n - 1
    tms = np.asarray(fp.times, dtype=float)
    w0 = 2.0 * np.pi / float(fp.T)
    out = np.zeros(m, dtype=complex)
    for i in range(m):
        u = np.zeros(m)
        u[i] = 1.0
        w, _ = pss._forced_replay(fp, freq, u)
        y0 = np.linalg.solve(A, alpha * np.asarray(w))
        _e, ys = pss._forced_replay(fp, freq, u, y0=y0, collect=True)
        y = [np.asarray(y0)[:m]] + [np.asarray(v)[:m] for v in ys]
        ## the DFT of `v = y exp(-j w t)`, which is what is T-periodic
        out[i] = sum(np.exp(-1j * (l * w0 + 2.0 * np.pi * freq) * tms[j])
                     * y[j][k] for j in range(N)) / N
    return out


def _sideband_pair(cir, per, freq, k, ls, npts=60, reltol=1e-11):
    """The adjoint rows and the `m` forward driven solves they replace."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method='gear', reltol=reltol)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    fp = pss.factored_period()
    n = fp.width
    N = len(fp.steps)
    alpha = np.exp(-2j * np.pi * freq * per)
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    A = np.eye(n) - alpha * M
    rows = PAC(cir, toolkit=circuit.numeric).adjoint_sideband_row(
        pss, freq, k, ls)
    fwds = np.array([_sideband_forward(pss, fp, freq, k, l, N, alpha, A)
                     for l in np.atleast_1d(ls)])
    return pss, rows, fwds


@pytest.mark.parametrize('l', [0, 1, -2])
def test_the_sideband_adjoint_needs_both_terms(l):
    """⚠ A SIDEBAND IS A FUNCTIONAL OVER THE WHOLE PERIOD, not a value at
    one instant, and that changes the adjoint in two ways.

    `H_l` is the `l`-th Fourier coefficient of `v(t) = y(t) exp(-j w t)`,
    the part that is actually T-periodic. Every state on the trajectory
    contributes, so the reverse pass takes an injection at EVERY step
    rather than a seed at the end — and the initial state `y_0` is itself a
    function of the source through the periodic boundary condition, which
    is a SECOND term:

        dH/du = [forced part, from the injected reverse pass]
              + [alpha * W^T z,  z = (I - alpha M)^-T g]

    ⚠ DROPPING THE SECOND TERM WOULD NOT LOOK WRONG. Measured, the two are
    comparable — 140 against 749 at `l = 0` — so an implementation with
    only the first returns a plausible number rather than an obviously
    broken one. Hence a comparison against `m` independent forward driven
    solves, not a reasonableness check.

    ⚠ AND THE CIRCUIT HAS TO BE NONLINEAR FOR THIS TO MEAN ANYTHING. On a
    linear circuit the linearisation is time-INvariant and every `l != 0`
    is zero, so an `l = 1` comparison would be two noise-level numbers
    agreeing. The diode makes the sidebands real: measured
    `|H_1|/|H_0| = 0.50` and `|H_-2|/|H_0| = 0.13`.
    """
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e4)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    _pss, rows, fwds = _sideband_pair(cir, 1e-6, 1e6, 1, [l])

    assert np.abs(fwds[0]).max() > 1e-3, \
        'sideband l=%d is numerically absent (%.3g), so this comparison ' \
        'would be two zeros agreeing -- the circuit is not mixing' \
        % (l, np.abs(fwds[0]).max())
    rel = np.linalg.norm(rows[0] - fwds[0]) / np.linalg.norm(fwds[0])
    assert rel < 1e-11, \
        'sideband l=%d: the adjoint row disagrees with %d forward driven ' \
        'solves by %.3e' % (l, cir.n - 1, rel)


def test_sidebands_vanish_on_a_linear_circuit():
    """⚠ THE CHECK THAT CAUGHT A CONVENTION ERROR, and could only have.

    A linear circuit's linearisation is time-INvariant however hard it is
    driven, so it converts nothing: `H_l = 0` for every `l != 0`. Measured
    here at ~10 orders below `H_0`.

    The first implementation decomposed `y(t)` rather than
    `v(t) = y(t) exp(-j w t)`. That is self-consistent — it agreed with a
    forward reference written the same way to 1e-15 — but it smears a
    Dirichlet kernel across every `l` whenever `f` is not a multiple of
    `1/T`, and reported `|H_1|` comparable to `|H_0|` on this very circuit.
    Only a case whose answer is known independently separates the two.
    """
    cir = _adjoint_ladder(4)
    _pss, rows, _f = _sideband_pair(cir, 1e-3, 700.0, 1, [0, 1, -2, 3])
    h0 = np.abs(rows[0]).max()
    assert h0 > 1.0, 'H_0 is empty (%.3g); nothing is being measured' % h0
    for li, l in enumerate([0, 1, -2, 3][1:], start=1):
        rel = np.abs(rows[li]).max() / h0
        assert rel < 1e-8, \
            'a LINEAR circuit converted %.3g of H_0 into sideband l=%d. ' \
            'It has no time-varying linearisation to convert with, so this ' \
            'is the sideband convention, not the circuit' % (rel, l)


def test_H0_is_the_ac_transfer_function_and_converges_like_gear():
    """`H_0` against a reference the adjoint path cannot influence.

    On a linear circuit the `l = 0` row IS the LTI transfer function from
    each source to the output — an `(sC + G)` solve, no monodromy, no
    replay, no adjoint. That makes it the strongest available check of the
    whole `H_l` path, and it is the one the pnoise build should be able to
    lean on.

    ⚠ THE RATE IS THE ASSERTION. The residual at 60 points is 1.2e-03,
    which in isolation says nothing — it could be a defect of any size.
    Per doubling of the grid it falls 4.06x / 4.03x / 4.02x, which is
    Gear-2's own O(h^2) and therefore pure discretisation.
    """
    from pycircuit.circuit.analysis import remove_row_col
    circuit.default_toolkit = circuit.numeric
    per, freq, k = 1e-3, 700.0, 1
    rels = []
    for npts in (60, 120, 240):
        cir = _adjoint_ladder(4)
        m = cir.n - 1
        _pss, rows, _f = _sideband_pair(cir, per, freq, k, [0], npts=npts,
                                        reltol=1e-12)
        irn = _pss.irefnode
        G = np.asarray(cir.G(np.zeros(cir.n)))
        Cm = np.asarray(cir.C(np.zeros(cir.n)))
        G, Cm = remove_row_col((G, Cm), irn, _pss.toolkit)
        sM = 2j * np.pi * freq * np.asarray(Cm) + np.asarray(G)
        ac = np.array([-np.linalg.solve(sM, np.eye(m)[:, i])[k]
                       for i in range(m)])
        rels.append(np.linalg.norm(rows[0] - ac) / np.linalg.norm(ac))

    for a, b in zip(rels, rels[1:]):
        ratio = a / b
        assert 3.4 < ratio < 4.6, \
            'H_0 against the AC transfer converges at %.2fx per doubling ' \
            '(%.4e -> %.4e), not the O(h^2) Gear-2 gives everywhere else. ' \
            'A wrong constant is invisible to a size check and obvious ' \
            'here' % (ratio, a, b)


def test_the_sideband_bound_is_hard():
    """⚠ THE TRUNCATION BOUND IS A REFUSAL, NOT A TOLERANCE.

    Okumura et al. eq. (32): the maximum frequency the analysis can speak
    about is the grid's own, so `|l| <= (w_max - w0)/ws`. Nothing aliases
    down from above what the grid represents. The abstract's "accumulated
    until their contributions become negligible" is a ratio test operating
    INSIDE that ceiling — an implementation carrying only the ratio test
    stops for the wrong reason, and on a coarse grid it would stop having
    summed harmonics the grid cannot represent.

    So this is refused rather than clamped: the remedy is a finer period
    grid, and silently returning the nearest representable sideband would
    answer a question the caller did not ask.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 40, maxiterations=40)
    fp = pss.factored_period()
    lmax = len(fp.steps) // 2
    pac = PAC(cir, toolkit=circuit.numeric)

    pac.adjoint_sideband_row(pss, 700.0, 1, lmax)      # at the edge: fine
    with pytest.raises(ValueError, match='Nyquist'):
        pac.adjoint_sideband_row(pss, 700.0, 1, lmax + 1)


def test_the_sideband_row_costs_one_solve_per_sideband():
    """The many-to-one property has to survive the distributed output.

    The injected reverse pass and the transposed solve both depend on `l`,
    so the cost is per SIDEBAND — but still not per SOURCE, which is the
    asymmetry pnoise is shaped by. This asserts the count, not the clock:
    this machine runs more than one agent.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(10)
    per, freq = 1e-3, 700.0
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / 60, maxiterations=40)
    pac = PAC(cir, toolkit=circuit.numeric)
    rows = pac.adjoint_sideband_row(pss, freq, 1, [0, 1, 2])
    assert rows.shape == (3, m)
    per_sideband = pac.matvecs / 3.0
    assert per_sideband < m, \
        'the sideband row took %.1f matvecs per sideband at m=%d; it is ' \
        'supposed to be independent of the number of sources' \
        % (per_sideband, m)


def _divider():
    """A linear divider with noisy resistors, and an independent answer."""
    c = SubCircuit()
    n1, n2 = c.add_nodes('net1', 'net2')
    c['vs'] = VSin(n1, gnd, va=1.0, vac=1.0, freq=1e3)
    c['R1'] = R(n1, n2, r=9e3)
    c['R2'] = R(n2, gnd, r=1e3)
    c['C'] = C(n2, gnd, c=1e-9)
    return c


def _diode_mixer():
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, vac=1.0, va=2.0, freq=1e6, phase=20)
    c['R'] = R(1, 2, r=1e4)
    c['D'] = Diode(2, gnd)
    c['C'] = C(2, gnd, c=1e-12)
    return c


def _pnoise_at(cir, per, fout, node, npts, **kw):
    import warnings
    circuit.default_toolkit = circuit.numeric
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index(node)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        S, used = pac.pnoise(pss, fout, k, **kw)
    return S, used, pac


def test_pnoise_reduces_to_the_stationary_analysis_on_a_linear_circuit():
    """⚠ OKUMURA'S `p = 1` CASE, against a reference pnoise cannot influence.

    A linear circuit converts nothing, so every sideband but `l = 0`
    vanishes and the fold collapses to the ordinary stationary formula —
    "exactly the same as that derived for a stationary noise". That makes
    `analysis_ss.Noise` the answer: a different analysis, an `(sC + G)`
    solve, no monodromy and no period anywhere in it.

    Measured: the `l != 0` terms come back ~1e-32 of the total, and the
    ratio against `Svnout` is 1.000000.

    ⚠ THE RATE, NOT THE SIZE. Per doubling of the period grid the residual
    falls 7.80x / 7.33x / 6.75x, reaching 3.9e-09 at 400 points. That is at
    least second order — the exact exponent is not asserted here because it
    has not been established, only that the disagreement is discretisation
    and not a defect.
    """
    from pycircuit.circuit.analysis_ss import Noise
    circuit.default_toolkit = circuit.numeric
    per, fout = 1e-3, 700.0
    rels = []
    for npts in (50, 100, 200):
        cir = _divider()
        ref = complex(Noise(cir, inputsrc='vs',
                            outputnodes=(cir.get_node('net2'), gnd)
                            ).solve(fout)['Svnout']).real
        S, used, _pac = _pnoise_at(cir, per, fout, cir.get_node('net2'), npts)
        assert S > 0, 'pnoise returned %r' % S
        rels.append(abs(S - ref) / ref)
        assert 0 in used and len(used) > 1, \
            'the accumulation never looked past l=0, so the fold is untested'

    assert rels[-1] < 1e-6, \
        'pnoise disagrees with the AC noise analysis by %.3e on a LINEAR ' \
        'circuit, where the two compute the same quantity' % rels[-1]
    for a, b in zip(rels, rels[1:]):
        assert a / b > 4.0, \
            'the residual against the stationary analysis falls only %.2fx ' \
            'per doubling (%.3e -> %.3e). A constant offset would look ' \
            'small and never move' % (a / b, a, b)


def test_pnoise_folds_and_the_fold_is_not_a_rounding_term():
    """On a mixer the sidebands carry a large share of the output noise.

    62% here, so an implementation that quietly summed only `l = 0` would
    be wrong by a factor, not by a tolerance — and would still return a
    plausible-looking PSD. The linear test above cannot catch that, because
    there the fold is *supposed* to contribute nothing.
    """
    cir = _diode_mixer()
    S_all, used, _p = _pnoise_at(cir, 1e-6, 3e5, 2, 80)
    S_l0, _u, _p2 = _pnoise_at(_diode_mixer(), 1e-6, 3e5, 2, 80,
                               maxsidebands=0)
    assert S_l0 > 0 and S_all > S_l0
    share = (S_all - S_l0) / S_all
    assert share > 0.25, \
        'the sidebands contribute only %.1f%% of the output noise on a ' \
        'driven diode; either the fold is not working or this circuit ' \
        'stopped mixing' % (100 * share)
    assert max(abs(np.asarray(used))) > 5, \
        'the accumulation stopped after |l|=%d on a switching circuit' \
        % max(abs(np.asarray(used)))


def test_pnoise_says_when_the_grid_stopped_it_rather_than_the_series():
    """⚠ WHICH STOPPING RULE FIRED IS PART OF THE ANSWER.

    Ending on the ratio test means the series converged. Ending on the
    grid's Nyquist means the grid ran out first, and every sideband above
    it is MISSING rather than small — the number is a lower bound. A
    strongly switching circuit reaches that readily: measured on this
    diode, 80 and 160 points per period both end on the bound, 320 and 640
    end on the ratio test, and the totals differ by only 0.04% — so the
    warning is not proof of a bad answer, it is a statement that the
    accumulation cannot vouch for itself.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _diode_mixer()
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / 80, maxiterations=40)
    irn = pss.irefnode
    k = cir.get_node_index(2)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    with pytest.warns(RuntimeWarning, match='Nyquist'):
        pac.pnoise(pss, 3e5, k)
    assert pac.alias_stop == 'bound'

    ## the linear circuit's series does converge, and says so
    _S, _u, pac2 = _pnoise_at(_divider(), 1e-3, 700.0,
                              _divider().get_node('net2'), 100)
    assert pac2.alias_stop == 'ratio', \
        'a linear circuit folds nothing and must stop on the ratio test, ' \
        'not on the grid'


def test_pnoise_refuses_cyclostationary_sources():
    """⚠ A BIAS-DEPENDENT `CY` IS A DIFFERENT MODEL, NOT A HARDER SUM.

    Okumura's cyclostationary treatment windows each source to a single
    timestep, and the windows' Fourier coefficients CORRELATE the
    sidebands: they stop adding in power and pick up cross terms needing
    the `R_{m,n}` construction. Summing powers anyway would answer a
    different question and look entirely normal doing it.

    Every noise source in the discrete element library is bias-independent
    — a resistor's `4kT/R` does not read `x` at all — so this refusal is
    unreachable there. A compact device's `CY` does read `x`, which is why
    the check samples the orbit rather than reasoning from element types.
    """
    class _BiasNoisyR(R):
        """A resistor whose noise follows its own terminal voltage."""
        def CY(self, x, w, epar=circuit.defaultepar):
            base = 4 * self.toolkit.kboltzmann * epar.T / self.iparv.r
            iPSD = base * (1.0 + 10.0 * abs(float(np.asarray(x).ravel()[0])))
            return self.toolkit.array([[iPSD, -iPSD], [-iPSD, iPSD]])

    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    n1, n2 = c.add_nodes('net1', 'net2')
    c['vs'] = VSin(n1, gnd, va=1.0, vac=1.0, freq=1e3)
    c['R1'] = _BiasNoisyR(n1, n2, r=9e3)
    c['R2'] = R(n2, gnd, r=1e3)
    c['C'] = C(n2, gnd, c=1e-9)
    pss = PSS(c, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
    assert pss.converged
    with pytest.raises(NotImplementedError, match='cyclostationary'):
        PAC(c, toolkit=circuit.numeric).pnoise(pss, 700.0, 1)


def _vdp_ppv(npts, mu=1.0):
    """A converged van der Pol and its PPV."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / npts,
                  x0=np.array([2.0, 0.0]), maxiterations=50)
    assert pss.converged, 'van der Pol did not converge at %d points' % npts
    v, info = pss.ppv()
    return c, pss, v, info


def test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has():
    """⚠ THE GATE FOR A2, and it is a physical experiment, not an identity.

    A PPV is only worth having if `v . delta` is the phase shift a real
    perturbation causes. So: displace van der Pol's state by `eps delta`,
    integrate the FULL NONLINEAR system for several periods until the
    transverse components have died (the second multiplier is 8.6e-4 per
    period), and read the surviving displacement along the orbit tangent.
    Nothing in that measurement touches the monodromy, the adjoint, or the
    bordered solve.

    ⚠ THE SCALE IS THE ASSERTION, NOT JUST THE DIRECTION. A direction check
    passes for any normalisation, and the normalisation is exactly what was
    in doubt. Measured against the true shift, per doubling of the period
    grid:

        npts   worst |1-ratio|   rel resid   fitted scale
         200      5.52e-02        1.90e-02     1.006266
         400      2.45e-02        8.94e-03     1.003290
         800      1.16e-02        4.35e-03     1.001677

    The fitted scale converges to ONE -- not to some other constant that a
    direction-only test would have accepted -- and the residual falls at
    O(h), which is the discretisation of the map the PPV is computed from.
    """
    import warnings
    from pycircuit.circuit.transient import Transient
    rng = np.random.default_rng(0)
    dirs = [d / np.linalg.norm(d) for d in
            (rng.standard_normal(2) for _ in range(4))]

    out = []
    for npts in (200, 400):
        cir, pss, v, info = _vdp_ppv(npts)
        m = cir.n - 1
        irn = pss.irefnode
        T = pss.period
        xdot = info['xdot']
        x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))

        def integrate(xi, nper=3, ppp=2000):
            ## ⚠ `reltol` HERE IS THE COST, not `ppp`. `Transient` adapts,
            ## so `timestep` is a first step and the tolerance sets the
            ## step count. At 1e-12 this test took 447 s; the signal being
            ## measured is ~9e-3, so 1e-9 is still three orders finer than
            ## anything it has to resolve.
            tran = Transient(cir, toolkit=circuit.numeric, reltol=1e-9,
                             iabstol=1e-13, vabstol=1e-11)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(refnode=gnd, tend=nper * T,
                                 timestep=T / ppp, x0=xi)
            return np.asarray(res.x, dtype=float)[:, -1]

        ref = integrate(x0f)
        meas, pred = [], []
        for d in dirs:
            eps = 1e-5
            dr = np.concatenate((d[:irn], np.zeros(1), d[irn:]))
            dx = np.delete(integrate(x0f + eps * dr) - ref, irn)
            meas.append(float(dx @ xdot) / float(xdot @ xdot) / eps)
            pred.append(float(v[:m] @ d))
        meas, pred = np.array(meas), np.array(pred)
        out.append((np.linalg.norm(meas - pred) / np.linalg.norm(meas),
                    float((pred @ meas) / (pred @ pred))))

    (r_coarse, _s_coarse), (r_fine, s_fine) = out
    assert abs(s_fine - 1.0) < 0.02, \
        'the PPV predicts the phase shift up to a factor of %.4f, not 1. ' \
        'The direction is right and the NORMALISATION is not -- which is ' \
        'the whole question this test exists to settle' % s_fine
    ratio = r_coarse / r_fine
    assert 1.6 < ratio < 2.8, \
        'the disagreement with the measured phase shift falls %.2fx per ' \
        'doubling (%.3e -> %.3e), not the O(h) that says it is ' \
        'discretisation of the period map' % (ratio, r_coarse, r_fine)


def test_the_ppv_normalisation_is_not_the_one_transcribed():
    """⚠ `v . q = 1` IS THE WRONG SCALE HERE, and it looks right.

    Demir's Remark 3.1 reads `v_1^T C(0) u_1(0) = 1`, and bordering the
    augmented system with `q = C xdot` makes `v . q = 1` fall out for free.
    It is a 7% error. The vector this bordered solve returns behaves as
    `C^T v_1` -- it contracts with a STATE perturbation directly -- so the
    two statements are both true of different objects, and using one where
    the other belongs is a silent scale error in every phase-noise number
    downstream.

    ⚠ THE OBVIOUS REPAIR IS ALSO WRONG, and was measured so before this
    normalisation was chosen: predicting the shift as `v^T C delta`
    gives residuals of 0.36 / 0.40 / 0.42 that GROW under refinement, with
    per-direction ratios scattering from -0.44 to 28.7. `v . delta` with
    `v . xdot = 1` converges at O(h). Do not "fix" this back without
    re-running that experiment.
    """
    cir, _pss, v, info = _vdp_ppv(400)
    m = cir.n - 1
    assert abs(float(v[:m] @ info['xdot']) - 1.0) < 1e-9, \
        'the returned PPV is not normalised by v . xdot = 1'
    vq = float(v[:m] @ info['q'])
    assert abs(vq - 1.0) > 0.5, \
        'v . q is %.6f, i.e. indistinguishable from the normalisation this ' \
        'test exists to rule out. Either C is now the identity on this ' \
        'circuit (in which case pick another) or the scale has been ' \
        'changed back' % vq


def test_the_ppv_border_residual_is_the_free_check():
    """`y` coming back zero is D&R's own correctness check, and it is real.

    With a zero first block on the right-hand side, `(I - M^T)v + y q = 0`
    forces `y q = 0`, so a nonzero `y` means the border absorbed a residual
    that belongs to the null space -- the computed `v` is not in it. Both
    solves report it, and both come back at ~1e-11.
    """
    _cir, _pss, _v, info = _vdp_ppv(200)
    for key in ('border_residual', 'tangent_border_residual'):
        assert abs(info[key]) < 1e-7, \
            '%s is %.3e; the bordered system absorbed a null-space ' \
            'residual, so the vector it returned is not the PPV' \
            % (key, info[key])
    assert info['null_residual'] < 1e-7, \
        '||v - M^T v|| / ||v|| is %.3e' % info['null_residual']


def test_the_ppv_refuses_what_has_no_phase():
    """A driven circuit's phase is its source's, not its own."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    assert pss.converged and not pss.autonomous
    with pytest.raises(ValueError, match='FREE-RUNNING'):
        pss.ppv()


def test_the_ppv_samples_are_the_propagated_adjoint():
    """`v(t)` over the whole period, checked by FORWARD machinery.

    Oscillator phase noise is an integral over the orbit — Demir's
    diffusion constant is `c = (1/T) ∫ v₁ᵀ(t) B(t)Bᵀ(t) v₁(t) dt` — so the
    PPV at `t = 0` is not enough. `Phi(T,s)ᵀ v(T) = v(s)`, and the reverse
    replay computes that sequence on its way to the answer; it was being
    discarded.

    ⚠ CHECKED AGAINST THE FORWARD STEP MAPS, not against itself. The
    propagator `Phi(T,s_j)` is rebuilt here by applying the forward
    recursion from step `j` to the end, and its transpose is compared with
    what the reverse pass recorded. Same shape as the existing `M` vs `Mᵀ`
    test, but per step, and it is what a future oscillator pnoise would be
    resting on. Measured worst case over all steps: 1.8e-15.
    """
    _cir, pss, v, info = _vdp_ppv(60)
    fp = pss.factored_period()
    n = fp.width
    m = pss.cir.n - 1
    samples = np.asarray(info['samples'])
    assert samples.shape == (len(fp.steps), n), \
        'expected one PPV sample per step, got %r' % (samples.shape,)

    cs0, cs1, ring = [], [], list(fp.opening)
    for _lu, C_new, _a, _b in fp.steps:
        cs0.append(ring[0])
        cs1.append(ring[1])
        ring = [C_new, ring[0]]

    def fwd_from(j, p):
        p0, p1 = p[:m].copy(), p[m:].copy()
        for k in range(j, len(fp.steps)):
            lu, _Cn, alphas, b = fp.steps[k]
            assert not b, 'this rebuild assumes the Gear-2 companion'
            pn = -lu.solve(alphas[1] * (cs0[k] @ p0)
                           + alphas[2] * (cs1[k] @ p1))
            p0, p1 = pn, p0
        return np.concatenate((p0, p1))

    ## the rebuild must be the monodromy, or it is checking nothing
    Mf = np.column_stack([fwd_from(0, e) for e in np.eye(n)])
    Mm = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    assert np.linalg.norm(Mf - Mm) / np.linalg.norm(Mm) < 1e-12, \
        'the forward rebuild is not the monodromy, so it cannot check the ' \
        'reverse pass'

    worst = 0.0
    for j in range(len(fp.steps)):
        Psi = np.column_stack([fwd_from(j, e) for e in np.eye(n)])
        pred = Psi.T @ v
        worst = max(worst, np.linalg.norm(pred - samples[j])
                    / max(np.linalg.norm(pred), 1e-300))
    assert worst < 1e-11, \
        'the reverse pass states are not the propagated adjoint (worst ' \
        '%.3e), so they are not the PPV over the period' % worst


def test_no_periodic_covariance_exists_for_an_oscillator():
    """⚠ WHY OSCILLATOR NOISE NEEDS THE PPV AND NOT A COVARIANCE SHOOT.

    Time-varying noise statistics are a Lyapunov ODE alongside the
    transient, and for a periodic large signal its periodic solution is a
    shooting problem whose monodromy is the KRONECKER SQUARE of the
    circuit's: `Phi_lyap = M ⊗ M`, so its multipliers are the pairwise
    products `lambda_i lambda_j`.

    For a DRIVEN circuit that is a single linear solve — the Lyapunov
    equation is linear in `K`, so no Newton iteration. For an AUTONOMOUS
    one it does not exist: `lambda_1 = 1` gives `lambda_1^2 = 1`, and
    `I - M ⊗ M` is exactly as singular as `I - M`. The covariance does not
    settle, it GROWS — variance linear in `t` is a random walk, which is
    phase diffusion, which is the linewidth.

    ⚠ SO THE UNIT MULTIPLIER IS NOT AN INCONVENIENCE HERE, IT IS THE
    ANSWER. The same near-unit-eigenvalue obstruction this codebase keeps
    meeting — eigen-selection, the bordered phase row, the PPV — appears
    once more as `lambda_1^2 = 1`, and this time what it obstructs is the
    wrong method for the question. Relayed measurements on three LTP
    systems put the kron identity at 2.2e-14 and the phase-mode growth at
    dead linear (trace K 1.825 -> 1793.4 over 1000 periods); this checks
    the structural half on our own monodromy.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    _cir, pss, _v, _info = _vdp_ppv(60)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    K = np.kron(M, M)

    ev = np.linalg.eigvals(M)
    outer = np.sort_complex(np.array([a * b for a in ev for b in ev]))
    got = np.sort_complex(np.linalg.eigvals(K))
    scale = max(float(np.max(np.abs(outer))), 1e-30)
    assert np.linalg.norm(got - outer) / scale < 1e-9, \
        'the covariance monodromy is not the Kronecker square of the ' \
        'circuit monodromy, so its multipliers are not the pairwise products'

    s_lyap = np.linalg.svd(np.eye(n * n) - K, compute_uv=False)[-1]
    s_circ = np.linalg.svd(np.eye(n) - M, compute_uv=False)[-1]
    assert s_lyap < 1e-7, \
        'I - M kron M has sigma_min %.3e on an AUTONOMOUS circuit. If the ' \
        'covariance shooting problem has become solvable, the unit ' \
        'multiplier has gone, and so has the oscillator' % s_lyap
    assert s_lyap < 100 * max(s_circ, 1e-16), \
        'the covariance system is singular to a different degree than the ' \
        'circuit one (%.3e against %.3e); the obstruction should be the ' \
        'same unit multiplier squared' % (s_lyap, s_circ)

    ## the contrast: a DRIVEN circuit has no such obstruction
    cir2 = _adjoint_ladder(4)
    pss2 = PSS(cir2, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    fp2 = pss2.factored_period()
    n2 = fp2.width
    M2 = np.column_stack([fp2.matvec(e) for e in np.eye(n2)])
    s2 = np.linalg.svd(np.eye(n2 * n2) - np.kron(M2, M2),
                       compute_uv=False)[-1]
    assert s2 > 1e-3, \
        'I - M kron M is near-singular (%.3e) on a DRIVEN circuit too, so ' \
        'the obstruction is not the unit multiplier after all' % s2


@pytest.mark.parametrize('mu', [0.5, 1.0])
def test_the_coloured_noise_functional_is_exactly_zero_here(mu):
    """⚠ A REGRESSION TEST WITH AN EXACT ANSWER OF ZERO, and it separates
    two functionals that share a PPV.

    Demir 2002 defines two different scalars from the same `v_1(t)`:

        c_w  = (1/T) ∫ v_1ᵀ B_w B_wᵀ v_1 dt     WHITE   — QUADRATIC
        V_0m = (1/T) ∫ v_1ᵀ B_cm dt             COLOURED — LINEAR

    The white one is the time-average of a SQUARE; the coloured one is the
    plain time-average — the zeroth Fourier coefficient of a periodic
    scalar. Using the quadratic form for a coloured source returns a
    plausible non-zero number from the same PPV, and nothing that did not
    know to look would catch it.

    §VIII gives the discriminating case: on a parallel-RLC oscillator with
    a nonlinear current source, "the time-average of [the Floquet vector
    entry] for the capacitor voltage is 0! Thus … any stationary …
    colored-noise source … connected across the capacitor has NO
    contribution to the oscillator spectrum due to phase noise, because
    V_0m = 0." Van der Pol is that circuit.

    ⚠ THE SECOND ASSERTION IS WHAT MAKES THIS A TEST RATHER THAN A
    TAUTOLOGY. A vector of zeros would pass the first one. The RMS of the
    same entries is ~0.40, so the WHITE functional is emphatically not
    zero in the same position — measured `|mean|/rms ~ 1e-11`. Zero and
    non-zero from one PPV, which is exactly the discrimination the
    coloured functional needs and the quadratic one destroys.

    (Both entries come back zero here, not just the capacitor's, which van
    der Pol's symmetry under `(v,i) -> (-v,-i)` predicts: the PPV is odd
    over the period.)
    """
    period = 6.6634 if mu >= 1.0 else 6.35
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=period, timestep=period / 200,
                  x0=np.array([2.0, 0.0]), maxiterations=60)
    assert pss.converged
    m = c.n - 1
    _v, info = pss.ppv()

    S = np.asarray(info['samples'])[:, :m]
    tms = np.asarray(info['times'], dtype=float)
    h = np.diff(tms)
    T = pss.period
    mean = (S * h[:, None]).sum(axis=0) / T           # the COLOURED scalar
    rms = np.sqrt(((S ** 2) * h[:, None]).sum(axis=0) / T)   # ~ the WHITE one

    assert (rms > 0.1).all(), \
        'the PPV samples are ~zero (rms %s), so the zero below would be ' \
        'vacuous' % np.array2string(rms, precision=4)
    ratio = np.abs(mean) / rms
    assert (ratio < 1e-8).all(), \
        'the time-average of the PPV is %s (relative %s), not zero. Demir ' \
        'section VIII says a coloured source on this oscillator ' \
        'contributes NO phase noise; a non-zero mean here means either ' \
        'the PPV or the averaging is wrong' \
        % (np.array2string(mean, precision=4),
           np.array2string(ratio, precision=4))


def _vdp_with_slow_node(tau_over_T=None, T=6.6634, mu=1.0):
    """van der Pol, optionally with one weakly coupled slow RC node.

    The coupling resistor is large, so the node barely loads the orbit —
    what it adds is a Floquet multiplier at `exp(-T/tau)`, which is the
    only thing under test.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    if tau_over_T is not None:
        c.add_node('w')
        rbig = 1e6
        c['Rs'] = R('v', 'w', r=rbig)
        c['Cs'] = C('w', gnd, c=tau_over_T * T / rbig)
    return c


def _solve_slow(tau_over_T):
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _vdp_with_slow_node(tau_over_T)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-11)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / 300, x0=x0,
                  maxiterations=60)
    assert pss.converged, 'tau/T=%r did not converge' % (tau_over_T,)
    return cir, pss


def test_a_slow_node_degrades_the_ppv_border_silently():
    """⚠ THE BORDER FIXES THE PHASE MODE AND NOTHING ELSE.

    `ppv()`'s bordered solve removes the singularity the UNIT multiplier
    causes. A second multiplier approaching 1 — which one slow node puts
    there — is untouched, and the conditioning goes with it. Measured on
    van der Pol with one weakly coupled RC node:

        tau/T   |lambda_2|   sigma_min(bordered)   null residual
        none     0.000856         8.62e-01            4.1e-11
        1e2      0.990049         4.49e-03            4.6e-11
        1e4      0.999900         4.47e-05            4.6e-11
        1e6      0.999999         4.47e-07            4.4e-11

    ⚠ `sigma_min` tracks `T/tau` over six decades WHILE THE RESIDUAL DOES
    NOT MOVE. GMRES converges, every diagnostic reads clean, and six
    digits of conditioning are gone. That is why `ppv()` estimates
    `|lambda_2|` explicitly instead of trusting a small residual — a
    residual cannot see this, and neither could this test if it only
    checked residuals.
    """
    smins, lam2s = [], []
    for tt in (1e2, 1e4):
        cir, pss = _solve_slow(tt)
        m = cir.n - 1
        fp = pss.factored_period()
        n = fp.width
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            _v, info = pss.ppv()
        q = info['q']
        qp = np.concatenate((q, np.zeros(n - m)))
        A = np.zeros((n + 1, n + 1))
        A[:n, :n] = np.eye(n) - M.T
        A[:n, n] = qp
        A[n, :n] = qp
        smins.append(np.linalg.svd(A, compute_uv=False)[-1])
        lam2s.append(info['second_multiplier'])
        assert info['null_residual'] < 1e-7, \
            'the residual should stay clean -- that is the whole point'

    assert lam2s[0] > 0.98 and lam2s[1] > 0.999, \
        'the slow node did not produce a near-unit multiplier (%s), so ' \
        'nothing is being tested' % lam2s
    decay = smins[0] / smins[1]
    assert 30 < decay < 300, \
        'sigma_min of the bordered system fell %.1fx for 100x the time ' \
        'constant (%.3e -> %.3e); the degradation is supposed to track ' \
        'T/tau' % (decay, smins[0], smins[1])


def test_the_ppv_says_when_a_second_multiplier_is_near_one():
    """The detector, because no residual can report this.

    Deflated power iteration against the eigenvectors `ppv()` already has:
    `u` spans the unit mode and `v` is its left partner, so the projection
    removes it exactly and what is left converges to `|lambda_2|`.
    Recovered to six digits against a dense eigendecomposition.

    ⚠ THE WARNING IS ABOUT MORE THAN CONDITIONING. The phase equation
    treats the oscillator's frequency response as INSTANTANEOUS, so a slow
    node that filters a nearby device's noise is invisible to it and phase
    noise comes out OVER-ESTIMATED. Better extraction does not fix that —
    Lai (Cadence) is explicit that the PPV "can be extracted correctly"
    and the analysis is "still inaccurate". So the result is an upper
    bound, and the warning says so rather than implying a tolerance would
    help.
    """
    import warnings
    for tt, expect in ((None, False), (1e0, False), (1e2, True)):
        cir, pss = _solve_slow(tt)
        fp = pss.factored_period()
        n = fp.width
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        true2 = np.sort(np.abs(np.linalg.eigvals(M)))[::-1][1]
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            _v, info = pss.ppv()
        got = info['second_multiplier']
        assert abs(got - true2) < 1e-4 * max(true2, 1e-3), \
            'the deflated power iteration reports |lambda_2| = %.6f ' \
            'against %.6f from a dense eigendecomposition' % (got, true2)
        fired = any('SECOND Floquet multiplier' in str(w.message)
                    for w in rec)
        assert fired is expect, \
            'tau/T=%r: |lambda_2| = %.6f and warned=%s' % (tt, got, fired)


def test_the_pac_operator_is_singular_at_every_harmonic_of_an_oscillator():
    """⚠ AT DC *AND* AT EVERY HARMONIC, AND ONLY FOR AN OSCILLATOR.

    `I - exp(-j w T) M` is what PAC solves. At `w = k w0` the factor is 1
    and the operator is `I - M`, which an autonomous circuit's unit
    multiplier makes singular. Measured on van der Pol:

        offset/f0    0      0.25   0.5    0.75    1       2       3
        sigma_min  2.8e-11  0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

    and LINEAR in the distance to the nearest harmonic — 2.5e-1, 2.6e-2,
    2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way. The same sweep
    on a DRIVEN ladder never drops below 0.78.

    ⚠ IT IS PHYSICS, NOT CONDITIONING, which is why PAC refuses rather
    than tightening a tolerance: a perturbation at a harmonic is a
    perturbation ALONG the orbit, and an oscillator answers that with
    unbounded phase drift. There is no bounded periodic response to
    return, so a number there would be a wrong answer rather than an
    imprecise one.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    _cir, pss = _solve_slow(None)                 # van der Pol, autonomous
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])

    def smin(ratio):
        a = np.exp(-2j * np.pi * ratio)
        return np.linalg.svd(np.eye(n) - a * M, compute_uv=False)[-1]

    for k in (0.0, 1.0, 2.0, 3.0):
        assert smin(k) < 1e-7, \
            'harmonic %g is not singular (%.3e); the unit multiplier has ' \
            'gone and so has the oscillator' % (k, smin(k))
    for r in (0.25, 0.5, 0.75):
        assert smin(r) > 0.1, \
            'offset %g is singular too (%.3e), so the singularity is not ' \
            'specific to the harmonics' % (r, smin(r))
    ## linear approach, not quadratic and not a cliff
    near = [smin(1.0 - d) for d in (1e-2, 1e-3, 1e-4)]
    for a, b in zip(near, near[1:]):
        assert 5 < a / b < 20, \
            'sigma_min falls %.1fx per decade of distance to the harmonic ' \
            '(%.3e -> %.3e), not the linear 10x' % (a / b, a, b)

    ## and PAC refuses there rather than returning something
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    f0 = 1.0 / pss.period
    with pytest.raises(ValueError, match='harmonic'):
        pac.adjoint_sideband_row(pss, f0, 0, 0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pac.adjoint_sideband_row(pss, 0.37 * f0, 0, 0)   # off-harmonic: fine

    ## a DRIVEN circuit has no such structure and is not refused
    cir2 = _adjoint_ladder(3)
    pss2 = PSS(cir2, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    fp2 = pss2.factored_period()
    n2 = fp2.width
    M2 = np.column_stack([fp2.matvec(e) for e in np.eye(n2)])
    for k in (0.0, 1.0, 2.0):
        s2 = np.linalg.svd(np.eye(n2) - np.exp(-2j * np.pi * k) * M2,
                           compute_uv=False)[-1]
        assert s2 > 1e-3, \
            'a DRIVEN circuit is singular at harmonic %g (%.3e)' % (k, s2)
    PAC(cir2, toolkit=circuit.numeric).adjoint_sideband_row(
        pss2, 1.0 / pss2.period, 1, 0)             # not refused


def test_pac_refuses_an_operating_point_from_another_circuit():
    """⚠ THE DRIVEN-OSCILLATOR TRAP, WHICH IS NOT A TYPO GUARD.

    The natural way to model a driven oscillator is to solve the PSS of
    the bare oscillator and treat the injection as a perturbation. It is
    wrong: the injection DEVICE is present even when its SIGNAL is zero.
    Buonomo & Lo Schiavo — "in absence of the injection signal, the
    injection circuit affects the basic LC oscillator by CHANGING THE
    NONLINEARITY OF THE FEEDBACK LOOP … [it] can affect the start-up
    condition … OR ITS OSCILLATION AMPLITUDE, or both."

    So the free-running orbit of the circuit-with-the-device is not the
    orbit of the circuit-without-it, and every Floquet quantity built on
    the wrong one inherits the error. The analysis would converge and
    report a plausible number.

    ⚠ THE REFERENCE-NODE CHECK CANNOT CATCH THIS. Two circuits differing
    by one device have the same reference node and, as here, the same node
    count — so the existing guard passes and the answer is quietly about
    the wrong orbit.
    """
    import warnings
    circuit.default_toolkit = circuit.numeric

    ## same topology, same node count, same refnode — different objects
    bare = _adjoint_ladder(3)
    other = _adjoint_ladder(3)
    pss = PSS(bare, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    assert pss.converged
    assert bare.n == other.n, 'the two circuits must look alike for this ' \
        'test to be about identity rather than shape'

    pac_other = PAC(other, toolkit=circuit.numeric)
    for call in (lambda: pac_other.solve(pss, [700.0]),
                 lambda: pac_other.adjoint_transfer_row(pss, 700.0, 1),
                 lambda: pac_other.adjoint_sideband_row(pss, 700.0, 1, 0),
                 lambda: pac_other.pnoise(pss, 700.0, 1)):
        with pytest.raises(ValueError, match='different circuit'):
            call()

    ## the matching pair is accepted
    pac_same = PAC(bare, toolkit=circuit.numeric)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pac_same.adjoint_sideband_row(pss, 700.0, 1, 0)


def test_the_am_pm_split_is_exact_and_the_conjugate_is_load_bearing():
    """⚠ ONE CONJUGATE SEPARATES AM FROM PM, and `a ± b` looks just as right.

    The two sidebands COUNTER-ROTATE about the carrier phasor, so the sum
    traces an ellipse: the component along the carrier is amplitude
    modulation, perpendicular is phase modulation. Pure AM keeps the
    envelope on the carrier's axis, forcing `a = conj(b)`; pure PM keeps it
    perpendicular, `a = −conj(b)`. Hence `m_am = a + conj(b)` and
    `m_pm = a − conj(b)`, each vanishing exactly when the other case holds.

    ⚠ WITHOUT THE CONJUGATE the split still produces two numbers and they
    are wrong for any modulation whose sidebands are not real relative to
    the carrier — it reports a rotating ellipse as pure AM. This test pins
    that by checking the naive form does NOT vanish where the correct one
    does, so a "simplification" back to `a ± b` fails here rather than in
    somebody's phase-noise number.
    """
    rng = np.random.default_rng(0)
    for _ in range(5):
        a = complex(rng.standard_normal(), rng.standard_normal())

        m_am, m_pm = PAC.am_pm_indices(a, np.conj(a))        # pure AM
        assert abs(m_pm) < 1e-14 * max(abs(m_am), 1.0), \
            'pure AM leaked %.3e into the PM index' % abs(m_pm)
        m_am2, m_pm2 = PAC.am_pm_indices(a, -np.conj(a))     # pure PM
        assert abs(m_am2) < 1e-14 * max(abs(m_pm2), 1.0), \
            'pure PM leaked %.3e into the AM index' % abs(m_am2)

    ## a case where the naive split is unambiguously wrong: pure AM with a
    ## complex modulation phase. `a + b` is then NOT zero-PM.
    a = complex(0.3, 0.7)
    b = np.conj(a)                       # pure AM by construction
    _m_am, m_pm = PAC.am_pm_indices(a, b)
    assert abs(m_pm) < 1e-14
    assert abs(a - b) > 0.5, \
        'the naive `a - b` should be far from zero on this pure-AM case, ' \
        'which is exactly why the conjugate cannot be dropped'


def test_am_pm_on_a_driven_mixer_is_predominantly_am():
    """A diode detector converts amplitude to amplitude; it has no free phase.

    So the modulation a small signal imposes on the carrier should come out
    overwhelmingly AM — measured `|m_pm|/|m_am| = 0.005` — and should be
    stable across modulation offsets, since neither index has a `1/f`
    mechanism behind it here. (An oscillator is the opposite case: its
    phase response goes as `1/ω_m`, so PM dominates near the carrier.)
    """
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _diode_mixer()
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / 120, maxiterations=40)
    assert pss.converged
    irn = pss.irefnode
    k = cir.get_node_index(2)
    k = k - 1 if k > irn else k
    pac = PAC(cir, toolkit=circuit.numeric)
    f0 = 1.0 / pss.period

    ## the carrier has to exist before it can be modulated
    C1 = pac.carrier_phasor(pss, k, 1)
    assert abs(C1) > 1e-3, 'no fundamental at the output (%s)' % C1

    ratios = []
    for r in (0.3, 0.1, 0.03):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            am, pm = pac.am_pm(pss, r * f0, k, 1)
        i = int(np.argmax(np.abs(am)))
        assert abs(am[i]) > 1.0, \
            'the AM index is %.3e -- nothing is being modulated' % abs(am[i])
        ratios.append(abs(pm[i]) / abs(am[i]))

    assert max(ratios) < 0.05, \
        'a diode detector came out %.4f PM/AM; it has no free phase to ' \
        'modulate' % max(ratios)
    assert max(ratios) / min(ratios) < 1.2, \
        'the PM/AM ratio moved %.2fx across offsets (%s); neither index has ' \
        'a 1/f mechanism here, so it should be flat' \
        % (max(ratios) / min(ratios), ratios)


def test_am_pm_refuses_a_harmonic_with_no_carrier():
    """AM/PM of nothing is not a small number, it is undefined."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _adjoint_ladder(3)
    pss = PSS(cir, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    pac = PAC(cir, toolkit=circuit.numeric)
    with pytest.raises(ValueError, match='no component at harmonic'):
        pac.am_pm(pss, 137.0, 1, carrier=57)     # far above anything present


@pytest.mark.parametrize('tau_over_T', [10.0, 100.0])
def test_the_oscillator_Q_is_recovered_from_the_second_multiplier(tau_over_T):
    """⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS, and it is already computed.

    An amplitude perturbation decays to `|lambda_2|` of its size each
    cycle, so the cycles needed to fall below a threshold IS the
    oscillator's Q: `Q = log(threshold)/log|lambda_2|`. The usual
    definitions do not apply to an autonomous circuit — `f_r/df` presumes a
    Bode plot of a BIBO-stable linear system, stored/dissipated presumes
    damping a self-sustaining oscillator does not have, and it is NOT the
    Q of the resonator inside it.

    ⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS MODULE WARNS
    ABOUT. "High Q", "a second multiplier near 1", "slow amplitude
    restoration" and "a long time constant" are four vocabularies for one
    thing — which is why the same circuits defeat the phase row, the
    eigen-split, the probe's continuation and the PPV's
    instantaneous-response assumption.

    ⚠ THE TEST IS A CONSTRUCTION CHECK, WHICH IS WHY IT IS WORTH HAVING.
    The slow node is built with a KNOWN `tau/T`, and its multiplier is
    `exp(-T/tau)`, so `Q` must come back as `tau/T` itself. Measured 10.00
    and 99.99 for 10 and 100 — an independent quantity reproducing an input
    the computation never sees.
    """
    _cir, pss = _solve_slow(tau_over_T)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _v, info = pss.ppv()
    Q = info['Q']
    assert abs(Q - tau_over_T) < 0.02 * tau_over_T, \
        'Q came back %.3f for a node built with tau/T = %g; Q is defined ' \
        'as cycles-to-1/e and the multiplier is exp(-T/tau), so the two ' \
        'are the same number' % (Q, tau_over_T)
    assert info['second_multiplier'] > 0.9, \
        'the slow node did not produce a near-unit multiplier'


def test_a_fast_oscillator_has_a_small_Q():
    """The control: plain van der Pol restores amplitude within a cycle.

    `|lambda_2| = 8.5e-04`, so Q ~ 0.14 — the perturbation is gone before
    the cycle ends, which is what "low Q" means for an oscillator and is
    the opposite corner from the case that breaks every method here.
    """
    _cir, pss = _solve_slow(None)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        _v, info = pss.ppv()
    assert info['Q'] < 1.0, \
        'plain van der Pol reports Q = %.3f; it restores amplitude in ' \
        'well under a cycle' % info['Q']


def _vdp_with_noise(psd=1e-6, mu=1.0):
    """van der Pol with a stationary white current source at its core.

    `i=0` so the orbit is untouched; only `CY` changes. This is the exact
    configuration the Monte Carlo gate measured, which is what lets the
    closed form be checked against a physical number.
    """
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=1.0)
    c['L'] = L('v', gnd, L=1.0)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: mu * (u - u ** 3 / 3.0))
    c['n'] = IS('v', gnd, i=0.0, noisePSD=psd)
    return c


def _solve_vdp_noise(npts=240, psd=1e-6):
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _vdp_with_noise(psd)
    m = cir.n - 1
    pss = PSS(cir, method='gear', reltol=1e-12)
    x0 = np.zeros(m)
    x0[0] = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=6.6634, timestep=6.6634 / npts, x0=x0,
                  maxiterations=60)
    assert pss.converged
    return cir, pss, PAC(cir, toolkit=circuit.numeric)


def test_the_diffusion_constant_matches_the_monte_carlo_measurement():
    """⚠ THE SCALAR THE WHOLE SPECTRUM IS BUILT ON, against a physical number.

    `c = (1/T) ∫ v₁ᵀ (CY/2) v₁ dt`, measured a completely different way —
    Monte Carlo on the full nonlinear circuit, phase read from
    zero-crossing timing so no PPV appears in the measurement at all —
    which gives `7.7083e-08` against `7.9516e-08` predicted, ratio 1.0316
    inside the measurement's 4.1% uncertainty at N=1200.

    ⚠ AN EARLIER VERSION OF THIS TEST PASSED AT 1.5903e-07, EXACTLY TWICE,
    against a Monte Carlo injecting `Var(i) = CY/h`. That measurement
    carried the same one-sided-as-two-sided convention the code did, so it
    confirmed the bug rather than catching it, at a ratio of 0.9965.
    `kT/C` — external to both — showed that injection reproducing 1.92× the
    right answer over ten runs. A measurement built on the assumption under
    test cannot test it.
    """
    _cir, pss, pac = _solve_vdp_noise()
    c = pac.diffusion_constant(pss)
    assert abs(c - 7.9516e-08) < 0.05 * 7.9516e-08, \
        'c = %.6e against 7.9516e-08, which a Monte Carlo with the ' \
        'CORRECT injection (Var(i) = CY/2h) measured at 7.7083e-08 -- ' \
        'ratio 1.0316, inside its 4.1%% uncertainty' % c
    ## ⚠ THE OLD VALUE HERE WAS 1.5903e-07, EXACTLY TWICE THIS, and it
    ## passed against a Monte Carlo that injected `Var(i) = CY/h`. That
    ## measurement carried the same one-sided-as-two-sided error as the
    ## code, so it confirmed the bug instead of catching it. `kT/C` --
    ## external to both -- showed that injection reproducing 1.92x the
    ## right answer over ten runs. A measurement built on the assumption
    ## under test cannot test it.
    assert abs(c * 2 - 1.5903e-07) < 0.05 * 1.5903e-07, \
        'the pre-fix value is no longer exactly twice this one; if the ' \
        'convention changed again, re-derive rather than re-fit'


def test_the_lorentzian_conserves_the_carrier_power_exactly():
    """⚠ THE INVARIANT THAT SEPARATES THIS FROM AN LTV TREATMENT.

    Noise spreads the carrier's power into a line of finite width; it does
    not create any. `∫ S_i df = 1` exactly, for every harmonic and every
    `c`. Analyses "based on linear time-invariant or linear time-varying
    concepts erroneously predict infinite noise power [at the carrier] as
    well as infinite total integrated power" — so this is the property
    that says the closed form is doing the nonlinear thing.

    Integrated numerically over the IMPLEMENTED function, not re-derived:
    a re-derivation would only be checking the algebra against itself.
    """
    from scipy.integrate import quad
    f0 = 150.0
    for c in (1e-9, 1e-7, 1e-5):
        for i in (1, 2, 5):
            tot, err = quad(lambda f: float(PAC.lorentzian(f, c, f0, i)),
                            -np.inf, np.inf, limit=400)
            assert abs(tot - 1.0) < 1e-6, \
                'harmonic %d at c=%g integrates to %.9f, not 1 — the ' \
                'carrier power is not conserved' % (i, c, tot)


def test_the_lineshape_is_lorentzian_where_it_should_be():
    """Finite at the carrier, `1/f²` far out, and the corner where predicted.

    Three properties, each of which a wrong constant would break
    differently: the peak is `1/(π² i² f₀² c)`, the far skirt falls as
    `1/f²` (a RATE, asserted as one), and the half-width is `π i² f₀² c`.
    """
    f0, c, i = 150.0, 1e-7, 1
    peak = float(PAC.lorentzian(0.0, c, f0, i))
    want = 1.0 / (np.pi ** 2 * i * i * f0 * f0 * c)
    assert abs(peak / want - 1.0) < 1e-12, \
        'peak %.6e against the analytic %.6e' % (peak, want)

    ## far out: doubling the offset must quarter the density
    far = [float(PAC.lorentzian(f, c, f0, i)) for f in (1e5, 2e5, 4e5)]
    for a, b in zip(far, far[1:]):
        assert abs(a / b - 4.0) < 0.02, \
            'the skirt falls %.3fx per doubling, not the 4x of 1/f²' % (a / b)

    ## half-width
    hw = np.pi * i * i * f0 * f0 * c
    assert abs(float(PAC.lorentzian(hw, c, f0, i)) / peak - 0.5) < 1e-12


def test_higher_harmonics_are_noisier_by_20log10i():
    """The skirt scales as `i²` and the corner as `i⁴`.

    So harmonic `i` sits `20 log₁₀(i)` dB above the fundamental far from
    the carrier — 6.02 dB for the second, 9.54 for the third. A designer
    reads that as "the divider makes it better, the multiplier worse", and
    it is a free consequence of the closed form rather than a separate
    calculation.
    """
    f0, c = 150.0, 1e-7
    far = 1e6
    base = float(PAC.lorentzian(far, c, f0, 1))
    for i in (2, 3, 5):
        db = 10.0 * np.log10(float(PAC.lorentzian(far, c, f0, i)) / base)
        want = 20.0 * np.log10(i)
        assert abs(db - want) < 0.05, \
            'harmonic %d is %.3f dB above the fundamental, not %.3f' \
            % (i, db, want)
        ## and its corner is i^4 wider
        hw_i = np.pi * i * i * f0 * f0 * c
        assert abs(float(PAC.lorentzian(hw_i, c, f0, i))
                   / float(PAC.lorentzian(0.0, c, f0, i)) - 0.5) < 1e-12


def test_the_oscillator_spectrum_is_built_and_refuses_a_driven_circuit():
    """End to end, and the one circuit class it does not describe."""
    import warnings
    _cir, pss, pac = _solve_vdp_noise()
    offs = np.array([1e-4, 1e-3, 1e-2, 1e-1])
    Sv, L = pac.oscillator_spectrum(pss, offs, 0, harmonic=1)
    assert np.all(np.isfinite(Sv)) and np.all(np.diff(Sv) < 0), \
        'the spectrum should be finite and falling with offset: %s' % Sv
    assert L[0] > L[-1], 'L(f) should fall with offset'
    ## far out it is 1/f^2, which in dB is -20 dB/decade
    slope = (L[-1] - L[-2]) / np.log10(offs[-1] / offs[-2])
    assert abs(slope + 20.0) < 1.0, \
        'the far skirt is %.2f dB/decade, not the -20 of 1/f²' % slope

    circuit.default_toolkit = circuit.numeric
    driven = _adjoint_ladder(3)
    p2 = PSS(driven, method='gear', reltol=1e-11)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        p2.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    with pytest.raises(ValueError, match='FREE-RUNNING'):
        PAC(driven, toolkit=circuit.numeric).diffusion_constant(p2)


def test_the_deflated_solve_removes_the_harmonic_singularity():
    """⚠ A7 — THE POLE IS THE ANSWER'S, AND IT IS NOW CARRIED ANALYTICALLY.

    At `α = 1` the operator `I − αM` is singular by the unit multiplier,
    and the solution really does diverge: an oscillator's phase response
    goes as `1/Δf`. What is wrong is *computing* a `1/ε` quantity through a
    system whose conditioning is also `1/ε` — the answer is genuinely large
    and the digits are genuinely gone.

    Bordering with the tangent and the PPV — both of which `ppv()` already
    returns — makes the border variable bounded, `s = (vᵀb)/(vᵀu)`, and the
    pole comes back through `u/(1 − α)` in closed form.

    Measured here: the plain operator's `σ_min` tracks the offset over nine
    decades while the bordered one is FLAT, and the solutions agree to
    ~1e-11 where the plain solve is still trustworthy.

    ⚠ THE TEST ASSERTS BOTH HALVES. Flat conditioning alone would be
    satisfied by an operator that had stopped solving the right problem, so
    agreement with the plain solve — in the regime where the plain solve
    can be believed — is what says it is still the same equation.
    """
    import warnings
    _cir, pss = _solve_slow(None)                # van der Pol, autonomous
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    rng = np.random.default_rng(0)
    b = rng.standard_normal(n) + 1j * rng.standard_normal(n)

    smins, rels = [], []
    for r in (1e-2, 1e-4, 1e-6):
        alpha = np.exp(-2j * np.pi * r)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            y = pac._deflated_solve(pss, alpha, b)
        smins.append(np.linalg.svd(np.eye(n) - alpha * M,
                                   compute_uv=False)[-1])
        yp = np.linalg.solve(np.eye(n) - alpha * M, b)
        rels.append(np.linalg.norm(y - yp) / np.linalg.norm(y))

    ## the plain operator really is degrading, or this proves nothing
    assert smins[0] / smins[-1] > 1e3, \
        'the plain operator only degraded %.1fx over four decades of ' \
        'offset (%s); the singularity being removed is not being exercised' \
        % (smins[0] / smins[-1], smins)

    ## and where the plain solve can still be believed, the two agree
    assert rels[0] < 1e-7, \
        'the deflated solve disagrees with the plain one by %.3e at an ' \
        'offset where the plain one is well conditioned — it is solving a ' \
        'different equation' % rels[0]

    ## the physical pole survives: |y| grows as 1/df
    mags = []
    for r in (1e-3, 1e-4, 1e-5):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mags.append(np.linalg.norm(
                pac._deflated_solve(pss, np.exp(-2j * np.pi * r), b)))
    for a, b_ in zip(mags, mags[1:]):
        assert 8.0 < b_ / a < 12.0, \
            'the response grows %.2fx per decade closer to the harmonic, ' \
            'not the 10x of a 1/df pole — the pole has been removed from ' \
            'the ANSWER rather than from the conditioning' % (b_ / a)


def test_the_deflated_solve_still_refuses_an_exact_harmonic():
    """Because the physical response there is unbounded.

    `1/(1 − α)` is a division by zero at a harmonic. The pole is removed
    from the CONDITIONING, not from the answer, so an exact harmonic is
    still a question with no finite answer.
    """
    _cir, pss = _solve_slow(None)
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    n = pss.factored_period().width
    with pytest.raises(ValueError, match='EXACT harmonic'):
        pac._deflated_solve(pss, 1.0 + 0.0j, np.ones(n))


def test_the_deflated_solve_works_transposed_too():
    """The adjoint rows need `(I − αMᵀ)`, whose borders swap.

    Its null space is spanned by the PPV and its left null space by the
    tangent, so `u` and `v` exchange roles. Checked against a dense
    transposed solve where that is still trustworthy.
    """
    import warnings
    _cir, pss = _solve_slow(None)
    pac = PAC(pss.cir, toolkit=circuit.numeric)
    fp = pss.factored_period()
    n = fp.width
    M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
    rng = np.random.default_rng(3)
    d = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    alpha = np.exp(-2j * np.pi * 1e-2)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        x = pac._deflated_solve(pss, alpha, d, transposed=True)
    xp = np.linalg.solve(np.eye(n) - alpha * M.T, d)
    rel = np.linalg.norm(x - xp) / np.linalg.norm(x)
    assert rel < 1e-7, \
        'the transposed deflated solve disagrees with a dense transposed ' \
        'solve by %.3e; the borders are probably not swapped' % rel


def _rc_noisy(Cval=1e-7, Rval=1e3, per=1e-3):
    """An RC lowpass with a noisy resistor, driven so a PSS exists.

    Linear, so its linearisation is time-invariant and the covariance is
    constant — and its capacitor-voltage variance is the exact `kT/C`,
    famously independent of `R`.
    """
    c = SubCircuit()
    c.add_node('a')
    c.add_node('b')
    c['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=1.0 / per)
    c['R'] = R('a', 'b', r=Rval)
    c['C'] = C('b', gnd, c=Cval)
    return c


def _cov_ratio(npts, Cval=1e-7, per=1e-3):
    import warnings
    from pycircuit.circuit.constants import kboltzmann
    circuit.default_toolkit = circuit.numeric
    cir = _rc_noisy(Cval=Cval, per=per)
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=per, timestep=per / npts, maxiterations=40)
    assert pss.converged
    pac = PAC(cir, toolkit=circuit.numeric)
    K0 = pac.covariance(pss)
    irn = pss.irefnode
    k = cir.get_node_index(cir.get_node('b'))
    k = k - 1 if k > irn else k
    T = float(circuit.defaultepar.T)
    return K0[k, k] / (kboltzmann * T / Cval), K0


def test_the_periodic_covariance_converges_to_kTC():
    """⚠ AN EXACT, FAMOUS, INDEPENDENT ANSWER: `Var(v_C) = kT/C`.

    Independent of `R`, which makes it a strong check — nothing about the
    resistor, the drive or the grid should appear in it. The covariance
    here is built from per-step Lyapunov accumulation and one Kronecker
    solve, and shares no machinery with the closed form it is checked
    against.

    ⚠ THE ASSERTION IS THE RATE. Measured 0.931 / 0.964 / 0.982 / 0.991 at
    100 / 200 / 400 / 800 points — the error halving each doubling, which
    is the O(h) a piecewise-constant approximation to white noise gives.
    A wrong constant would sit at a fixed offset and never move; this is
    what distinguishes the two.

    ⚠ AND IT SETTLED A FACTOR OF TWO BY MEASUREMENT. With the full `CY` the
    same sequence converges to 2.0, not 1.0 — so `CY` is a ONE-SIDED
    density and the per-step injection carries `CY/2`. Argument alone
    would not have decided it; the two candidate conventions differ by
    exactly the factor the test resolves.
    """
    rs = [_cov_ratio(n)[0] for n in (100, 200, 400)]
    assert rs[-1] > 0.95, \
        'the covariance is %.4f of kT/C at 400 points; it should be ' \
        'approaching 1' % rs[-1]
    errs = [abs(1.0 - r) for r in rs]
    for a, b in zip(errs, errs[1:]):
        assert 1.6 < a / b < 2.6, \
            'the error falls %.2fx per doubling (%s), not the ~2x of O(h). ' \
            'A wrong CONSTANT would not move at all' % (a / b, errs)


def test_the_covariance_grid_must_resolve_the_noise_bandwidth():
    """⚠ A PRECONDITION, NOT AN ACCURACY NOTE — and it cost a wrong reading.

    The first attempt at the gate above came back at 0.517 and looked like
    a factor-of-two bug. It was not: the RC pole sat at 159 kHz while the
    grid's Nyquist was 100 kHz, so the DISCRETE system genuinely does not
    carry the noise the continuous one does. `kT/C` requires integrating
    past the pole.

    Reproduced here deliberately: the same circuit with the pole above
    Nyquist reads far low, and moving the pole below it recovers the
    answer. A `kT/C` that comes back low is the grid, not the code.
    """
    ## pole at 159 kHz, Nyquist at 100 kHz — under-resolved
    under, _K = _cov_ratio(200, Cval=1e-9)
    ## pole at 1.6 kHz, Nyquist at 100 kHz — resolved
    ok, _K2 = _cov_ratio(200, Cval=1e-7)
    assert under < 0.7, \
        'the under-resolved case reads %.4f; it is supposed to be visibly ' \
        'low, or this test is not demonstrating the precondition' % under
    assert ok > 0.95, 'the resolved case reads %.4f' % ok


def test_the_covariance_refuses_an_oscillator():
    """`I − M⊗M` is singular there, and that is the physics.

    The unit multiplier squares to one, the covariance grows without bound
    rather than settling, and that growth IS the phase diffusion. An
    oscillator's output noise is stationary, not cyclostationary — there is
    no object here to compute, so it refuses and names the right route.
    """
    _cir, pss = _solve_slow(None)
    with pytest.raises(ValueError, match='no periodic covariance'):
        PAC(pss.cir, toolkit=circuit.numeric).covariance(pss)


def test_the_covariance_samples_are_periodic_and_positive():
    """The time-varying statistic itself: `K(t)` over one period.

    Two structural properties that a wrong recursion breaks differently —
    it must return to itself after a period (that is what the Kronecker
    solve enforces), and every `K_j` must be a positive-semidefinite
    covariance, which nothing in the solve guarantees a priori.
    """
    _r, _K = _cov_ratio(100)
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = _rc_noisy()
    pss = PSS(cir, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
    K0, seq = PAC(cir, toolkit=circuit.numeric).covariance(pss, samples=True)

    rel = np.linalg.norm(seq[-1] - K0) / max(np.linalg.norm(K0), 1e-300)
    assert rel < 1e-8, \
        'K(T) differs from K(0) by %.3e -- the solve is supposed to ' \
        'enforce exactly that periodicity' % rel
    for j, K in enumerate(seq):
        w = np.linalg.eigvalsh(K)
        assert w.min() > -1e-12 * max(abs(w).max(), 1e-300), \
            'K at step %d has eigenvalue %.3e; a covariance cannot be ' \
            'negative definite' % (j, w.min())


def _osc_cov(npts=240):
    _cir, pss, pac = _solve_vdp_noise(npts=npts)
    K_orb, d, info = pac.oscillator_covariance(pss, samples=True)
    return pss, pac, K_orb, d, info


def test_the_oscillator_covariance_predicts_the_walk_forty_periods_out():
    """⚠ THE ASSERTION IS A PREDICTION, NOT A PROPERTY OF THE SOLVE.

    `covariance` refuses an oscillator because `I − M⊗M` is singular there
    and no periodic covariance exists. The claim this makes instead is a
    split — a bounded part plus a random walk along the orbit:

        K(t₀ + nT) = K_orb + n·d·u·uᵀ

    which is falsifiable in a way "the residual is small" is not: run the
    real Lyapunov recursion forward for forty periods (9,600 steps, from
    `K = 0`, touching nothing the bordered solve produced) and compare.
    A bordered system can always be solved; only this says the answer
    means anything.

    ⚠ THE FIRST PERIOD IS THE LOOSEST, AND THAT IS PHYSICS. Starting from
    `K = 0` rather than `K_orb` leaves a transient in the bounded part; it
    decays with `|λ₂| = 8.5e-4` and is gone by period two. Measured
    2.6e-07, 1.1e-10, 2.2e-10, 5.4e-10, 1.2e-09, 2.5e-09 at periods
    1/2/5/10/20/40 — the walk dominates 99.9% of the trace by then, so the
    late numbers test the growth term and the early one tests the bound.
    """
    pss, pac, K_orb, d, info = _osc_cov()
    u = info['tangent_pair']
    As, Qs, _K1, M, _m, n = pac._lyapunov_pieces(pss, 'test')
    ## the split's own periodicity statement: P is periodic UP TO the
    ## growth, which appears exactly once per period -- not periodic
    P = info['orbital_samples']
    grow = d * np.outer(u, u)
    assert np.max(np.abs(P[-1] - K_orb - grow)) < 1e-10 * np.max(np.abs(grow))

    K = np.zeros((n, n))
    seen = {}
    for p in range(1, 41):
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        if p in (1, 2, 40):
            pred = K_orb + p * grow
            seen[p] = np.max(np.abs(K - pred)) / np.max(np.abs(pred))
    assert seen[1] < 1e-5, \
        'one period off by %.3e; the bounded part is wrong' % seen[1]
    assert seen[40] < 1e-7, \
        'forty periods off by %.3e; the growth RATE is wrong -- that ' \
        'error accumulates where the first-period one does not' % seen[40]
    assert seen[2] < seen[1], \
        'the initial transient does not decay, so the split is not a ' \
        'decomposition into a settling part plus a walk'


def test_the_growth_rate_is_the_diffusion_constant_by_another_route():
    """⚠ TWO SEPARATELY ANCHORED QUANTITIES, CLOSED INTO A LOOP.

    A phase deviation `α` displaces the state by `α·u`, so the growing
    covariance is `Var(α)·u uᵀ = c·t·u uᵀ` and therefore `d = c·T`. The two
    sides share the `CY/2` convention and nothing else: `c` is a quadratic
    form in the ADJOINT-replayed PPV, `d` comes from a FORWARD Lyapunov
    recursion closed by a bordered Kronecker solve.

    ⚠ AND THEIR ANCHORS ARE INDEPENDENT TOO, which is the property that
    was missing when a 2× error survived a 0.9965 agreement. The
    injection behind `d` is pinned by `kT/C`; `c` is pinned by a nonlinear
    Monte Carlo reading phase from zero crossings. Neither measurement can
    influence the other's answer.

    ⚠ ASSERTED AS A CONVERGENCE, NOT A TOLERANCE. Both carry an O(h)
    piecewise-constant approximation to white noise, so at any single grid
    they differ by a real amount; what must hold is that the difference
    HALVES per doubling. Measured 1.87% → 1.03% → 0.54%. A shared error
    would cancel and give a flat ratio of 1; a wrong one would not shrink.
    """
    errs = []
    for npts in (120, 240, 480):
        _cir, pss, pac = _solve_vdp_noise(npts=npts)
        _K, _d, info = pac.oscillator_covariance(pss)
        c = pac.diffusion_constant(pss)
        errs.append(abs(info['c_from_growth'] / c - 1.0))
    assert errs[0] < 0.03, 'coarsest grid off by %.3f' % errs[0]
    for a, b in zip(errs, errs[1:]):
        assert b < 0.7 * a, \
            'the disagreement %.4f -> %.4f is not first-order; a ' \
            'difference that does not converge away is a defect, not ' \
            'discretisation' % (a, b)


def test_the_bordered_kronecker_solve_matches_its_closed_form():
    """`d = (vᵀK₁v)/(v·u)²` — the `O(n²)` contraction behind the `O(n⁴)` solve.

    Left-multiplying the bordered system by `(v⊗v)ᵀ` annihilates the
    singular block, so the growth rate never needed the Kronecker at all.
    They are the same quantity by construction, which makes a disagreement
    diagnostic rather than a precision question: it would mean the border
    pair is not the null pair.

    ⚠ AND THE DEFLATION IS WHAT MAKES EITHER COMPUTABLE. `I − M⊗M` has
    `σ_min = 2.3e-11` against a next singular value of 0.997 — a cleanly
    one-dimensional null space, which is why bordering with a single pair
    is the right repair. Bordered, `σ_min` comes back to 5.4e-02: nine
    orders recovered, the same shape as A7's deflated PAC solve.
    """
    _pss, _pac, _K, d, info = _osc_cov()
    assert info['d_residual'] < 1e-9, \
        'the solve and the closed form differ by %.3e; the border pair ' \
        'is not the null pair' % info['d_residual']
    assert info['sigma_min'] < 1e-8, \
        'sigma_min = %.3e; I - M kron M is supposed to be SINGULAR here ' \
        '-- if it is not, this circuit is not autonomous' % info['sigma_min']
    assert info['sigma_min_bordered'] > 1e-3, \
        'bordered sigma_min = %.3e; the deflation did not recover the '\
        'conditioning' % info['sigma_min_bordered']
    assert info['null_residual'] < 1e-8


def test_the_pair_inner_product_is_not_one_and_d_scales_with_the_tangent():
    """⚠ THE 2.31× THAT WAS CHASED AS A CODE DEFECT, PINNED AS A NUMBER.

    `ppv()` normalises on the FIRST BLOCK, `v[:m]·ẋ = 1` — right for a
    perturbation entering the first block, which is where an injected
    current lands and what every shipped path does. The FULL PAIR
    contraction is a different number, ≈0.663, so `(v·u)⁻² ≈ 2.27`.
    Assuming the pair product is 1 is exactly the error that produced a
    2.31× discrepancy between two Monte Carlo routes.

    ⚠ AND `d` ALONE IS NOT AN INVARIANT. Rescaling `u → s·u` sends
    `d → d/s²`; only the product `d·u uᵀ` is a property of the circuit.
    Asserted directly, because a future change to how the tangent is
    scaled would silently move `d` while leaving every structural check
    passing — `info['growth']` is what downstream code should read.
    """
    _pss, _pac, _K, d, info = _osc_cov()
    vu = info['pair_inner']
    assert 0.5 < abs(vu) < 0.9, \
        'v.u = %.4f; if this has become 1.0 the pair normalisation ' \
        'changed and every d is off by %.3f' % (vu, 1.0 / vu ** 2)
    u = info['tangent_pair']
    assert np.allclose(info['growth'], d * np.outer(u, u))
    ## the invariant, stated as a rescaling that must not move it
    s = 3.0
    assert abs((d / s ** 2) * np.outer(s * u, s * u)
               - d * np.outer(u, u)).max() < 1e-18


def test_the_oscillator_covariance_refuses_a_driven_circuit():
    """The mirror of `covariance`'s refusal, and it is not symmetry for its
    own sake. A driven circuit's `I − M⊗M` is nonsingular, so there is no
    null direction to border with and no walk to split off; bordering it
    anyway would return a `d` near zero and an arbitrary `K_orb`, which is
    a plausible wrong answer rather than an error.
    """
    import warnings
    cir = _rc_noisy()
    pss = PSS(cir, method='gear')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-3, timestep=1e-3 / 400, refnode=gnd)
    with pytest.raises(ValueError, match='GROWS'):
        PAC(cir, toolkit=circuit.numeric).oscillator_covariance(pss)


def _vdp_scaled(cval, lval, period, npts=400):
    """van der Pol with reactances that are NOT unity — see the test below."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c.add_node('v')
    c['C'] = C('v', gnd, c=cval)
    c['L'] = L('v', gnd, L=lval)
    c['B'] = BSource('v', gnd, gnd, 'v',
                     i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    pss = PSS(c, method='gear', reltol=1e-12)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=period, timestep=period / npts,
                  x0=np.array([2.0, 0.0]), maxiterations=80)
    assert pss.converged
    v, info = pss.ppv()
    irn = pss.irefnode
    x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
    x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
    Cf = np.asarray(c.C(x0f))
    Cm = np.delete(np.delete(Cf, irn, 0), irn, 1).astype(float)
    return c, pss, v, info, Cm


def test_the_ppv_fixture_is_blind_to_a_missing_C_from_one_direction():
    """⚠ THIS TESTS THE TEST, and it found a real hole in the fixture.

    `v·δ` (right) and `vᵀCδ` (the transcription of Demir's Remark 3.1 that
    cost 7% before it was measured out) are distinguished by
    `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has` —
    but only because it perturbs in four RANDOM directions. Van der Pol as
    shipped has `C = diag(1, −1)`, so along the capacitor node the two
    formulations are **numerically identical**: `v·e₀ = vᵀCe₀` exactly,
    ratio 1.0000. Any single-direction probe at the capacitor — which is
    what two earlier gates in this campaign actually did — cannot see the
    difference at all.

    ⚠ "VERIFIED TO 1e-15" SAYS NOTHING ABOUT WHICH ERRORS A FIXTURE CAN
    SEE. A unit reactance makes `C` the identity up to a sign, and an
    implementation that drops the `C` weighting then passes every check
    exactly. The repair is to run the same circuit at a reactance that is
    not 1: with `c = 2` the same blind direction separates the two
    formulations by exactly the capacitance.

    Asserted on both fixtures on purpose — the blindness is recorded as a
    measured property of the shipped one, not as a hypothesis about it, so
    that a future test written against van der Pol knows what it is
    choosing when it perturbs at the capacitor.
    """
    _c1, _p1, v1, _i1, Cm1 = _vdp_scaled(1.0, 1.0, 6.6634)
    m = 2
    e0 = np.array([1.0, 0.0])
    plain1 = float(v1[:m] @ e0)
    weighted1 = float(v1[:m] @ (Cm1 @ e0))
    assert abs(weighted1 / plain1 - 1.0) < 1e-12, \
        'the shipped fixture was expected to be BLIND here (ratio 1); it ' \
        'now reads %.6f, so C is no longer diag(1,-1) and this test has ' \
        'stopped describing the fixture' % (weighted1 / plain1)

    _c2, _p2, v2, _i2, Cm2 = _vdp_scaled(2.0, 3.0, 6.6634 * np.sqrt(6.0))
    plain2 = float(v2[:m] @ e0)
    weighted2 = float(v2[:m] @ (Cm2 @ e0))
    assert abs(weighted2 / plain2 - 2.0) < 1e-9, \
        'ratio %.6f; at c = 2 the missing-C error must show as exactly ' \
        'the capacitance from the SAME direction the unit fixture cannot ' \
        'see it from' % (weighted2 / plain2)
    ## and the normalisation itself must survive the rescaling -- the
    ## property under test is the fixture's discriminating power, not a
    ## claim that a scaled circuit is solved differently
    assert abs(float(v2[:m] @ _i2['xdot']) - 1.0) < 1e-9


def test_the_ppv_gate_probes_the_direction_of_maximum_sensitivity():
    """⚠ THE DESIGNED PROBE, replacing four random directions and hope.

    `v·δ` and `vᵀCδ` agree exactly when `vᵀ(C − I)δ = 0`, so the set of
    directions blind to a dropped `C` is a HYPERPLANE with a known normal:

        blind  ⟺  δ ⊥ (C − I)ᵀv

    which makes `(C − I)ᵀv` itself the direction of maximum sensitivity,
    available for one matrix-vector product from quantities already in
    hand. Van der Pol's `C = diag(1, −1)` gives `(C − I)ᵀv = (0, −2v₁)`, so
    `e₀` — the capacitor node — is orthogonal to it EXACTLY. That is why
    the companion test measures a ratio of 1.0000 there: not a near miss,
    an exact one.

    ⚠ ALONG THE DESIGNED PROBE THE TWO HYPOTHESES PREDICT OPPOSITE SIGNS,
    so no tolerance is needed to separate them — the circuit picks one:

        npts   measured      v·δ (right)   vᵀCδ (wrong)   ratio
         400   −5.4215e-01   −5.4131e-01   +5.4131e-01    0.998443
         800   −5.4153e-01   −5.4110e-01   +5.4110e-01    0.999213

    The magnitude converges at O(h) as a bonus; the SIGN alone already
    decides it. Four random directions scored 0.43/0.78/0.81/0.96 on this
    fixture — they work, and the spread is exactly the luck this removes.

    ⚠ THIS IS A TARGETED PROBE AND ITS OPTIMALITY IS ABOUT ONE ERROR. It
    maximises sensitivity to a dropped `C` weighting and says nothing
    about any other defect; the random-direction gate stays because it is
    not aimed at a hypothesis. Recording which errors a check can see is
    the whole point of §D shape 0c, and that applies to this one too.

    ⚠ AND THE STRUCTURAL RULE IS WORTH MORE THAN THIS CIRCUIT. `C` is ZERO
    on algebraic rows, so `C − I = −I` there and the discrepancy is
    maximal: in an MNA-shaped system the rows a DAE solver already treats
    specially are the BEST probes, and the capacitor nodes everyone
    reaches for first are the blind ones.
    """
    import warnings
    from pycircuit.circuit.transient import Transient
    out = []
    for npts in (400, 800):
        cir, pss, v, info = _vdp_ppv(npts)
        m = cir.n - 1
        irn = pss.irefnode
        T = pss.period
        xdot = info['xdot']
        x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        Cm = np.delete(np.delete(np.asarray(cir.C(x0f)), irn, 0),
                       irn, 1).astype(float)

        normal = (Cm - np.eye(m)).T @ v[:m]
        assert abs(float(normal @ np.eye(m)[0])) < 1e-12 * max(
            np.linalg.norm(normal), 1e-300), \
            'the capacitor axis is no longer exactly in the blind ' \
            'hyperplane, so this fixture has changed shape'
        d = normal / np.linalg.norm(normal)

        def integrate(xi, nper=3, ppp=2000):
            tran = Transient(cir, toolkit=circuit.numeric, reltol=1e-9,
                             iabstol=1e-13, vabstol=1e-11)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tran.solve(refnode=gnd, tend=nper * T,
                                 timestep=T / ppp, x0=xi)
            return np.asarray(res.x, dtype=float)[:, -1]

        eps = 1e-5
        ref = integrate(x0f)
        dr = np.concatenate((d[:irn], np.zeros(1), d[irn:]))
        dx = np.delete(integrate(x0f + eps * dr) - ref, irn)
        meas = float(dx @ xdot) / float(xdot @ xdot) / eps
        p_ok = float(v[:m] @ d)
        p_bug = float(v[:m] @ (Cm @ d))

        assert np.sign(p_ok) != np.sign(p_bug), \
            'the designed probe no longer separates the two formulations ' \
            'by sign; it is not the maximum-sensitivity direction'
        assert np.sign(meas) == np.sign(p_ok), \
            'the CIRCUIT chose the vTCd formulation (measured %+.4e, ' \
            'v.d %+.4e); the PPV normalisation is wrong' % (meas, p_ok)
        out.append(abs(p_ok / meas - 1.0))

    assert out[0] < 5e-3, 'coarse grid off by %.3e' % out[0]
    assert out[1] < 0.75 * out[0], \
        'the designed probe does not converge (%.3e -> %.3e); a residue ' \
        'that does not shrink is a defect, not discretisation'\
        % (out[0], out[1])


class _TransCap(Circuit):
    """A two-node element with a deliberately NON-SYMMETRIC `C`.

    `q = C v` with `C[0,1] != C[1,0]`. Physically this is a
    transcapacitance: `∂q_a/∂v_b` need not equal `∂q_b/∂v_a` once the
    charge is partitioned between terminals, which is what Ward-Dutton
    does in a MOS channel. It exists here because every other fixture in
    this file has a symmetric `C` and therefore cannot express the
    difference — see the test below.
    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c11', desc='', unit='F', default=1e-12),
                  Parameter(name='c12', desc='', unit='F', default=0.0),
                  Parameter(name='c21', desc='', unit='F', default=0.0),
                  Parameter(name='c22', desc='', unit='F', default=1e-12)]

    def update(self, subject):
        p = self.iparv
        self._C = self.toolkit.array([[p.c11, p.c12], [p.c21, p.c22]])

    def C(self, x, epar=None):
        return self._C

    def G(self, x, epar=None):
        return self.toolkit.zeros((2, 2))

    def i(self, x, epar=None):
        return self.toolkit.zeros(2)

    def q(self, x, epar=None):
        return self._C @ x


def _transcap_pss(c12, c21, npts=60):
    """A driven pair whose `C` is symmetric or not, as asked."""
    import warnings
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('a')
    cir.add_node('b')
    cir['vs'] = VSin('a', gnd, va=1.0, vac=1.0, freq=1e6)
    cir['R1'] = R('a', 'b', r=1e3)
    cir['R2'] = R('b', gnd, r=1e3)
    cir['Cb'] = C('b', gnd, c=1e-12)
    cir['T'] = _TransCap('a', 'b', c11=2e-12, c12=c12, c21=c21, c22=2e-12)
    pss = PSS(cir, method='gear', reltol=1e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        pss.solve(period=1e-6, timestep=1e-6 / npts, refnode=gnd)
    assert pss.converged
    return cir, pss


def _replay_transposed_by_hand(pss, use_T):
    """The reverse recursion, with the `C` transpose optionally omitted."""
    fp = pss.factored_period()
    m = pss.cir.n - 1
    n = fp.width
    cs0, cs1, ring = [], [], list(fp.opening)
    for _lu, C_new, _a, _b in fp.steps:
        cs0.append(ring[0])
        cs1.append(ring[1])
        ring = [C_new, ring[0]]

    def one(v):
        w1, w2 = v[:m].copy(), v[m:].copy()
        for j in range(len(fp.steps) - 1, -1, -1):
            lu, _Cn, alphas, _b = fp.steps[j]
            t = lu.solve_transposed(w1)
            A0 = cs0[j].T if use_T else cs0[j]
            A1 = cs1[j].T if use_T else cs1[j]
            w1, w2 = (-alphas[1] * (A0 @ t) + w2,
                      -alphas[2] * (A1 @ t))
        return np.concatenate((w1, w2))

    return np.column_stack([one(e) for e in np.eye(n)])


def test_the_transposed_replay_gate_cannot_see_a_dropped_transpose():
    """⚠⚠ THE GATE UNDER THE WHOLE ADJOINT STACK IS BLIND ON EVERY FIXTURE
    IN THIS FILE, and the code it certifies is correct only by inspection.

    `_monodromy_matvec_transposed` does `cs0[j].T @ t`, and its recorded
    gate is "MEASURED against the dense `Mᵀ` … agreement 1.8e-15". That
    gate is passed IDENTICALLY by an implementation with no transpose in
    it, because van der Pol and the RC pair both have a symmetric `C`:

        fixture            shipped code    with the .T DROPPED
        symmetric C          3.994e-16       3.994e-16   ← blind
        NON-symmetric C      3.994e-16       4.667e-01   ← caught

    ⚠ AND THE ERROR CLASS IS REAL IN THIS TREE, not hypothetical.
    `compact.PspMosLongChannel`'s `C` is non-symmetric — `Cgd = −4.31 fF`
    against `Cdg = −0.13 fF`, a factor of 33, which is Ward-Dutton charge
    partition and the model being RIGHT. So the moment a PSS is run on a
    MOS circuit the transpose matters, and until then nothing in the suite
    would change state to say so.

    ⚠ QUOTED AGAINST `Cox`, BECAUSE THAT IS THE COMPARABLE NORMALISATION
    AND THE FIRST ATTEMPT USED A DIFFERENT ONE. McAndrew's figure is a
    nonreciprocity over `Cox` (`|C_ij − C_ji| ≲ 0.01·Cox`); the 0.44 this
    once quoted was a ratio to `max|C|`, and 33 is a spread between two
    entries — three different denominators. `Cox` is anchored OUTSIDE the
    `C()` code by the model's own geometry: `eps_ox W L / tox` with
    `W = L = 1e-6`, `tox = 2.2e-9` gives 15.70 fF, 2.2% from the measured
    max `Cgg`.

    ⚠⚠ AND FIXING THE DENOMINATOR DID NOT FIX THE COMPARISON, BECAUSE THE
    CONDITIONS WERE ALSO MISMATCHED. McAndrew's bound is stated at
    `VDS = 0`, under his eq. (2), which is derived there. The 4.97 fF
    above is a WORST CASE over `Vg, Vd in [0, 1.2]` — a box containing
    saturation, where Ward-Dutton partition makes `Cgd`/`Cdg` asymmetric
    BY DESIGN. Measured at his condition instead:

        Vds = 0, worst over Vg      0.844% of Cox   =  0.84x his 0.01
        Vds = 1.2, Vg = 1.2        31.66% of Cox    = 31.66x

    with the growth monotone in `Vds` — 0.47, 1.63, 3.77, 9.47, 20.9,
    30.6, 31.7 (× 0.01) at `Vds` = 0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.2.

    ⚠ SO HIS 1% IS VERIFIED, NOT CONTRADICTED, and the earlier note here
    claiming it was "a floor for the ideal case, not a typical value" was
    wrong. A real compact model, gated to 1.3e-6 against a compiled
    PSP103, satisfies his bound at his stated condition to within a factor
    of 1.2. The 32x is a statement about SATURATION — which is what this
    gate cares about, and not a claim about the paper.

    ⚠ THE NONRECIPROCITY IS ESSENTIALLY CREATED BY `Vds`: 68x growth from
    `Vds = 0` to `Vds = 1.2`. That is why a DC-biased fixture would be a
    weak transpose gate and a swinging one is a strong one.

    ⚠ A SYMMETRIC `C` MAKES THIS A TOTAL BLIND SPOT RATHER THAN A WEAK
    PROBE. `C = Cᵀ` means no perturbation direction, random or designed,
    separates the two implementations — the discriminating quantity is
    `C − Cᵀ` and it is identically zero. Contrast the dropped-`C`
    hyperplane, where a bad direction exists alongside good ones.
    """
    for label, (c12, c21), blind in (
            ('symmetric', (-1e-12, -1e-12), True),
            ('transcapacitive', (-1.7e-12, -0.3e-12), False)):
        _cir, pss = _transcap_pss(c12, c21)
        fp = pss.factored_period()
        n = fp.width
        Mf = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        scale = max(float(np.max(np.abs(Mf))), 1e-300)
        shipped = float(np.max(np.abs(
            _replay_transposed_by_hand(pss, True) - Mf.T))) / scale
        dropped = float(np.max(np.abs(
            _replay_transposed_by_hand(pss, False) - Mf.T))) / scale
        assert shipped < 1e-12, \
            '%s: the shipped recursion is wrong (%.3e)' % (label, shipped)
        if blind:
            assert dropped < 1e-12, \
                'the symmetric fixture was expected to be BLIND; it now ' \
                'separates by %.3e, so its C is no longer symmetric and ' \
                'this test has stopped describing it' % dropped
        else:
            assert dropped > 1e-3, \
                'the transcapacitive fixture no longer sees a dropped ' \
                'transpose (%.3e); it is the only thing in this file ' \
                'that can' % dropped
        ## and the shipped path itself, through the public entry point
        Mt = np.column_stack([fp.matvec_transposed(e) for e in np.eye(n)])
        assert float(np.max(np.abs(Mt - Mf.T))) / scale < 1e-12


def test_the_compact_mos_models_really_are_transcapacitive():
    """The premise of the test above, pinned so it cannot go stale quietly.

    If a future change symmetrised these models' `C`, the blindness note
    would become harmless and nobody would know to remove it — and, worse,
    a real physical effect would have been lost. Ward-Dutton partition in
    saturation decouples the drain from the channel charge, so `∂q_g/∂v_d`
    is small while `∂q_d/∂v_g` is not. That asymmetry is the model being
    right, not a defect.
    """
    from pycircuit.circuit import compact
    circuit.default_toolkit = circuit.numeric
    worst = 0.0
    for name in ('PspMosLongChannel', 'PspPmosLongChannel'):
        cls = getattr(compact, name)
        inst = cls(*cls.terminals)
        for vg in np.linspace(0.0, 1.2, 7):
            for vd in np.linspace(0.0, 1.2, 7):
                Cm = np.asarray(inst.C(np.array([vd, vg, 0.0, 0.0])),
                                dtype=float)
                s = float(np.max(np.abs(Cm)))
                if s == 0.0:
                    continue
                worst = max(worst, float(np.max(np.abs(Cm - Cm.T))) / s)
    assert worst > 0.1, \
        'max|C-C^T|/max|C| = %.3e; the compact MOS models are no longer ' \
        'transcapacitive, so the transpose blind spot recorded against ' \
        'them needs re-deriving rather than deleting' % worst

    ## ⚠ AND `Cox` ITSELF, ANCHORED BY GEOMETRY RATHER THAN BY `C()`.
    ## The comparable form of McAndrew's bound is `|C_ij - C_ji| / Cox`,
    ## so `Cox` has to come from somewhere the code under test cannot
    ## move it: `eps_ox W L / tox` from the model's own parameters.
    inst = compact.PspMosLongChannel(*compact.PspMosLongChannel.terminals)
    p = inst.iparv
    cox_geom = 3.9 * 8.8541878128e-12 * p.w * p.l / p.tox
    cgg = max(abs(float(np.asarray(inst.C(np.array([0.0, vg, 0.0, 0.0])),
                                   dtype=float)[1, 1]))
              for vg in np.linspace(0.0, 2.5, 26))
    assert abs(cgg / cox_geom - 1.0) < 0.05, \
        'max Cgg = %.4e against eps_ox W L / tox = %.4e; if these have ' \
        'parted company the Cox normalisation below is no longer ' \
        'anchored by geometry' % (cgg, cox_geom)
    nonrecip = 0.0
    for vg in np.linspace(0.0, 1.2, 13):
        for vd in np.linspace(0.0, 1.2, 13):
            Cm = np.asarray(inst.C(np.array([vd, vg, 0.0, 0.0])),
                            dtype=float)
            nonrecip = max(nonrecip, float(np.max(np.abs(Cm - Cm.T))))
    assert 25.0 < nonrecip / cgg / 0.01 < 40.0, \
        'nonreciprocity is %.1fx McAndrew 1%% of Cox; the recorded 32x ' \
        'no longer describes this model' % (nonrecip / cgg / 0.01)

    ## ⚠ AND THE SAME QUANTITY AT McANDREW'S OWN CONDITION, `Vds = 0`,
    ## which is where his bound is stated and where it HOLDS. Pinned
    ## separately from the saturation figure on purpose: they are two
    ## different claims and conflating them is what produced a wrong
    ## reading of the paper.
    at_vds0 = 0.0
    for vg in np.linspace(0.0, 1.2, 13):
        Cm = np.asarray(inst.C(np.array([0.0, vg, 0.0, 0.0])), dtype=float)
        at_vds0 = max(at_vds0, float(np.max(np.abs(Cm - Cm.T))))
    assert at_vds0 / cgg < 0.01, \
        'at Vds = 0 the nonreciprocity is %.4f of Cox, above McAndrew 1%% ' \
        'where he states it. That would be the STRONGER claim -- his ' \
        'bound failing for a real compact model -- and it needs saying ' \
        'so, not silently absorbing into the saturation figure'\
        % (at_vds0 / cgg)
    assert nonrecip / at_vds0 > 20.0, \
        'the nonreciprocity is no longer created by Vds (ratio %.1f); the ' \
        'transpose gate depends on the fixture SWINGING, not on its DC ' \
        'bias' % (nonrecip / at_vds0)


def test_the_transcap_fixtures_sensitivity_is_linear_not_thresholded():
    """⚠ HOW SMALL A TRANSCAPACITANCE `_TransCap` CAN STILL CATCH.

    A designed probe's discriminating power degrades with weak asymmetry,
    which raises the fair question of whether this fixture is a gross-error
    gate. It is not, and the reason is structural rather than lucky: the
    dropped-transpose error is proportional to `C − Cᵀ` itself, so it
    scales LINEARLY with the asymmetry and has no threshold.

        asym    dropped-.T error      ratio to asym
        0.700      4.6667e-01            0.6667
        0.350      2.3333e-01            0.6667
        0.100      6.6667e-02            0.6667
        0.020      1.3333e-02            0.6667
        0.005      3.3333e-03            0.6667
        0.000      3.9942e-16            (the floor)

    ⚠ AND THAT IS A DIFFERENT STRUCTURE FROM A PROBE-DIRECTION FAILURE,
    which is why the "weak asymmetry defeats it" result does NOT transfer
    here. A probe degrades because a DIRECTION becomes misaligned — a
    geometric effect with a distribution over draws. This degrades because
    the QUANTITY shrinks, deterministically, with no distribution at all.
    At McAndrew's own 1% the error is still 1e-2 against a 4e-16 floor:
    fourteen orders of margin. The discriminator dies only at exactly
    zero, which is what makes it a §D shape 0d rather than a weak test.
    """
    base = -1.0e-12
    out = []
    for asym in (0.35, 0.1, 0.02):
        d = asym * 2e-12 / 2.0
        _cir, pss = _transcap_pss(base - d, base + d)
        fp = pss.factored_period()
        n = fp.width
        Mf = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        scale = max(float(np.max(np.abs(Mf))), 1e-300)
        err = float(np.max(np.abs(
            _replay_transposed_by_hand(pss, False) - Mf.T))) / scale
        out.append(err / asym)
        assert err > 1e-4, \
            'asym %.3f gives only %.3e; the fixture has become a ' \
            'gross-error gate' % (asym, err)
    assert max(out) / min(out) - 1.0 < 1e-6, \
        'the sensitivity is not linear in the asymmetry (%s); a ' \
        'threshold would mean weakly transcapacitive models slip past'\
        % np.round(out, 6)
