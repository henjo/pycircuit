"""Shooting tests: shooting pss.  Split out of test_analysis_shooting.py on
2026-09-27 (in its original order); shared helpers are in
`_shooting_fixtures.py`, the HDL elements in `_shooting_elements.py`.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    ConvergenceWarning,
    UsageWarning,
)
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools
import warnings
from pycircuit.circuit.tests._shooting_elements import (_PllMultPd,
    _PllPhaseDiv)
from pycircuit.circuit.tests._shooting_fixtures import (_a10_vdp,
    _diode_fixture,
    _diode_mixer,
    _force_plain_map,
    _pac_circuit,
    _pss_lte,
    _pss_peak,
    _pwm_loop,
    _q20_rlc,
    _rc_ladder,
    _scaled_vdp,
    _shooting_trace,
    _varying_c_ladder,
    _vdp_with_noise)


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

    ## (N steps, N + 1 points: `timestep = period / N` is N steps)
    t,dt = numeric.linspace(0,period,num=N + 1,endpoint=True,
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


def test_pss_uses_a_solve_not_an_explicit_inverse(monkeypatch):
    """`inv(Jf) @ C @ Jshoot` formed a dense inverse per timestep per iteration.

    The quantity wanted is the solution of `Jf X = C @ Jshoot`, and at
    N=137/M=1000 with 20 shooting iterations the old form was 20,000 dense
    inversions.  A timing test would be a flake, so this asserts what the
    walk does: it runs with `np.linalg.inv` refusing, and each step's
    sensitivity is ONE matrix right-hand side solved through the caller's
    `linearsolver` (a sparse matrix-free path must not be measured against a
    hard-wired dense baseline).  The plain map and gear's pair both.
    """
    from pycircuit.circuit.linearsolver import DenseSolver

    class _Counting(DenseSolver):
        def __init__(self):
            self.rhs = []

        def solve(self, A, b, toolkit):
            self.rhs.append(np.shape(b))
            return DenseSolver.solve(self, A, b, toolkit)

    def _refuse(*a, **k):
        raise AssertionError('the explicit inverse is back')

    circuit.default_toolkit = circuit.numeric
    for method, kind, width in (('trap', 'plain', 1), ('gear2', 'pair', 2)):
        cir = _q20_rlc()
        per, m = 1e-3, cir.n - 1
        solver = _Counting()
        pss = PSS(cir, method=method, reltol=1e-9, linearsolver=solver)
        pss._open_at_x0 = False
        pss.autonomous = False
        assert pss._map_kind() == kind
        times, hs = pss._period_grid(per, 40, None)
        with monkeypatch.context() as mp, quiet():
            mp.setattr(np.linalg, 'inv', _refuse)
            M = pss._walk(kind, np.zeros(width * m), times, hs,
                          T=per).monodromy()
        assert np.all(np.isfinite(M)) and M.shape == (width * m, width * m)
        ## one matrix solve per step of the loop (the plain map's opening
        ## step is the manufacturing step, not a step of the map)
        matrix = [r for r in solver.rhs if len(r) == 2]
        assert matrix == [(m, width * m)] * (len(times) - 1), (method, matrix)


def test_pss_still_matches_the_ac_reference_with_a_fine_step():
    """Both methods must converge to the same, correct answer.

    At a coarse step neither is reliable and Euler's closeness is coincidence --
    measured, it is 0.9886 of the AC answer at dt = RC but 1.3283 at dt = RC/4.
    With a fine enough step both land on 1.0000, which is what makes the
    resonator comparison above a statement about damping rather than about luck.
    """
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
    with quiet():
        ac = AC(build()).solve(freqs=np.array([f]))
    ref = abs(complex(np.asarray(ac.v(2, gnd)).ravel()[0]))

    for method in ('euler', 'trap'):
        with quiet(AccuracyWarning, ConvergenceWarning):
            res = PSS(build(), method=method).solve(period=1 / f,
                                                    timestep=1 / f / 1280)
        v = np.asarray(res['tpss'].v(2, gnd), dtype=float)
        amp = 0.5 * (v.max() - v.min())
        assert amp == pytest.approx(ref, rel=0.02), \
            '%s gives %.6f against the AC reference %.6f' % (method, amp, ref)


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
    circuit.default_toolkit = circuit.numeric

    def run(**kw):
        with quiet(AccuracyWarning):
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
    pss = PSS(_q20_rlc(), method='euler', reltol=1e-9, steadyratio=100.0)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-5, maxiterations=40, trace=True)
    trace = [float(np.max(np.abs(F))) for _z, F, _J in pss.shooting_trace]

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

    with quiet(AccuracyWarning):
        res = PSS(rect(), method='euler', reltol=1e-6).solve(
            period=per, timestep=per / n, maxiterations=20)
    t_p = np.asarray(res['tpss'].sweep_values, dtype=float)
    v_p = np.asarray(res['tpss'].v('c'), dtype=float).ravel()

    ## the diode must actually conduct -- the defect this replaces returned
    ## a waveform six orders smaller than this bound
    assert v_p.max() > 1.0, \
        'the rectifier never conducted: v(c) peaks at %.3e' % v_p.max()

    with quiet():
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
    pss = PSS(_q20_rlc(), method='trap')
    with quiet(AccuracyWarning, ConvergenceWarning):
        pss.solve(period=1e-3, timestep=1e-5, maxiterations=3)
    tr = pss._transient()
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    assert isinstance(tr.base_integrator, TrapezoidalIntegrator)
    assert tr._effective_method in ('TrapezoidalIntegrator',
                                    'EulerIntegrator')


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
                                              ThetaIntegrator,
                                              Gear2Integrator, TRBDF2Integrator,
                                              RadauIIA3Integrator,
                                              ESDIRK43Integrator)
    circuit.default_toolkit = circuit.numeric
    want = {'euler': EulerIntegrator, 'trap': TrapezoidalIntegrator,
            'trapezoidal': TrapezoidalIntegrator,
            'theta': ThetaIntegrator,
            'gear': Gear2Integrator, 'gear2': Gear2Integrator,
            'trbdf2': TRBDF2Integrator, 'radau': RadauIIA3Integrator,
            'esdirk43': ESDIRK43Integrator}
    for name, cls in want.items():
        tr = PSS(_q20_rlc(), method=name)._transient()
        assert isinstance(tr.par.integrator, cls), \
            'method=%r selected %s' % (name, type(tr.par.integrator).__name__)

    with pytest.raises(ValueError,
                       match="'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                             "'radau'"):
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
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
    with quiet(AccuracyWarning):
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


def _pss_plain(method, timestep=1e-5, reltol=1e-3, **kw):
    """Force the pre-augmentation formulation, where a seam can exist.

    Gear-2 now solves for its entering history, so its seam is gone by
    construction -- which is the fix, and which leaves the seam machinery
    with nothing to observe unless the old formulation can still be run.
    This is deliberately a test-level override rather than a Parameter: a
    user has no reason to ask for the formulation that measured 1.266e-01 V
    of avoidable error.
    """
    circuit.default_toolkit = circuit.numeric
    pss = _force_plain_map(PSS(_q20_rlc(), method=method, reltol=reltol, **kw))
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
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric

    def run():
        pss = PSS(_q20_rlc(), method='gear', reltol=1e-3)
        with quiet():
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
    circuit.default_toolkit = circuit.numeric

    def evals(force_plain):
        pss = PSS(_q20_rlc(), method='gear', reltol=1e-9)
        if force_plain:
            _force_plain_map(pss)
        with quiet(AccuracyWarning):
            pss.solve(period=1e-3, timestep=1e-5, maxiterations=40,
                      trace=True)
        assert pss.converged
        return len(pss.shooting_trace)

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
    circuit.default_toolkit = circuit.numeric

    def solved(method, n):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
        with quiet(AccuracyWarning):
            ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
            res = pss.solve(period=1e-3, timestep=1e-3 / (n - 1), maxiterations=30)
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

    ## ⚠ THE ORBIT IS ITS DIFFERENTIAL STATE, FOLDED; THE SCALE IS ITS SPEED
    ## AT x_0 (2026-09-23).  The detector used to compare every row against
    ## three times the LARGEST step on the orbit.  An idtmod VCO's wrapped
    ## output jumps a whole modulus in one step, and its unfolded phase
    ## never recurs, so every one-fold solve warned "19.8 times" and true
    ## multiples were named 19.9 / 8.65 / 11.96.  A van der Pol seeded at
    ## its turning point, moving slowly there, warned "24.0 times".
    from pycircuit.circuit.elements_hdl import VcoHdl

    def vco(k):
        c = SubCircuit()
        for nd in ('vco', 'ph', 'ctl'):
            c.add_node(nd)
        c['X1'] = VcoHdl('ctl', gnd, 'vco', gnd, 'ph', f0=1e6, kvco=0.0,
                         va=1.0, modulus=1.0)
        c['Rc'] = R('ctl', gnd, r=1e3)
        c['Rl'] = R('vco', gnd, r=1e3)
        names = [str(nd) for nd in c.nodes if str(nd) != 'gnd!']
        x0 = np.zeros(c.n - 1)
        x0[names.index('X1._state0')] = x0[names.index('ph')] = 0.1
        x0[names.index('vco')] = np.sin(0.2 * np.pi)
        pss = PSS(c, method='radau', reltol=1e-9)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            pss.solve(period=k * 1.03e-6, timestep=1.03e-6 / 100, x0=x0,
                      maxiterations=40)
        assert pss.converged and abs(pss.period - k * 1e-6) < 1e-12
        return pss, [c_ for c_ in caught if 'MULTIPLE' in str(c_.message)]

    one, hits = vco(1)
    assert not hits and one.fundamental_period is None, \
        'the one-fold VCO was called a multiple: %r' % (hits,)
    two, hits = vco(2)
    assert hits and abs(two.period / two.fundamental_period - 2.0) < 0.05, \
        (hits, two.fundamental_period)

    cir, _mu = _a10_vdp()
    vdp = PSS(cir, method='radau', reltol=1e-12)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        vdp.solve(period=2 * np.pi, timestep=2 * np.pi / 50,
                  x0=np.array([2.0, 0.0]), maxiterations=100)
    assert vdp.converged
    assert not [c_ for c_ in caught if 'MULTIPLE' in str(c_.message)], \
        'a van der Pol seeded at its slow turning point was called a multiple'

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
    circuit.default_toolkit = circuit.numeric
    ## BDF-2's own roots, so the claim above is checked and not asserted
    assert np.allclose(sorted(np.roots([1.5, -2.0, 0.5])), [1.0 / 3.0, 1.0])

    pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
    with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    analytic = float(np.exp(-np.pi / 20.0))
    for method in ('trap', 'gear'):
        for kind in ('2:1', 'smooth'):
            pss = PSS(_q20_rlc(), method=method, reltol=1e-9)
            with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'gear'):
        for kind in ('2:1', 'smooth'):
            pss = PSS(_phase_circuit(), method=method, reltol=1e-8)
            with quiet(AccuracyWarning):
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

    The plain walk builds `x(0)` from the unknown with one order-dropped Euler
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

    ## a grid opening far coarser than its finest step gains a RAMP of
    ## doubling steps from its finest step up to the coarse one (2026-09-21:
    ## ONE step of `fr.min()` before the coarse remainder was a growth of
    ## `fr[0]/fr.min()` -- 138x on a relaxation van der Pol's own `lte_grid`,
    ## beyond a two-step method's zero-stability bound, so the integrator
    ## dropped two steps to Euler and gear's transposed replay refused the
    ## grid gear had produced).  Every ratio in the ramp is at most 2.
    fr = np.concatenate(([0.5], np.full(500, 0.001)))
    fr = fr / fr.sum()
    times, hs = pss._period_grid(1e-3, len(fr), fr)
    nramp = len(hs) - len(fr) + 1
    assert 5 <= nramp <= 12, nramp                       # ~log2(500) halvings
    ## top-down (2026-09-21): the first piece is the coarse step halved down
    ## to its finest scale -- within a factor 2 of `fr.min()` -- doubled up,
    ## so the hand-off to the next step is a ratio of 2 too
    assert fr.min() * 1e-3 / 2.0 < hs[0] <= fr.min() * 1e-3 * (1.0 + 1e-12)
    assert np.sum(hs[:nramp]) == pytest.approx(fr[0] * 1e-3, rel=1e-12)
    assert np.all(hs[1:nramp + 1] / hs[:nramp] <= 2.0 + 1e-12)
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
    circuit.default_toolkit = circuit.numeric
    T = 162.842412                      # measured free-running period

    ## one period of ACCEPTED steps from a settled adaptive transient
    cir = _van_der_pol()
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0
    with quiet():
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
    with quiet(AccuracyWarning):
        pss.solve(period=T, x0=seed, grid=fr, maxiterations=25)

    assert pss.converged, \
        'van der Pol did not solve through PSS.solve(grid=...) on its own ' \
        'LTE-chosen grid -- item 5 has no payoff case without this'
    err = 1e6 * (pss.period - T) / T
    assert abs(err) < 200, \
        'solved T=%.6f, %.1f ppm from the measured %.6f' % (pss.period, err, T)
    ## fewer than 1200 steps, against the 20000 a uniform grid needs
    assert len(pss.times) < 1200


def test_a_matrix_free_solve_agrees_with_the_dense_one():
    """The two paths must answer the same, or the fast one is not an option.

    ⚠ The convergence TEST is not identical between them -- `fsolve` scales
    its residual by `|J| . |x|`, which matrix-free has no way to form -- so
    this asserts on the converged ANSWER and the converged/not verdict,
    which are what a caller sees, and not on the iteration count, which
    measurably differs (matrix-free takes one more at m >= 502).
    """
    circuit.default_toolkit = circuit.numeric
    out = []
    for mf in (False, True):
        pss = PSS(_varying_c_ladder(), method='gear', reltol=1e-6)
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    cases = ((_q20_rlc(), 'trap', 1e-3, 100),        # driven plain
             (_rc_ladder(6), 'gear', 1e-3, 50),      # driven solved-history
             (_phase_circuit(), 'trap', 1e-3, 200),  # autonomous plain
             (_phase_circuit(), 'gear', 1e-3, 100))  # autonomous composed
    for cir, method, period, npts in cases:
        pss = PSS(cir, method=method, reltol=1e-8)
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric

    def split_at(npts):
        pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    T = 6.663293                      # measured free-running period, mu=1

    cir = _scaled_vdp()
    iref = cir.get_node_index(gnd)
    iv = cir.get_node_index('v')
    x0 = np.zeros(cir.n)
    x0[iv] = 2.0
    with quiet():
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
    with quiet():
        pss._begin_period(seed)
        x1 = np.asarray(pss.solve_timestep(seed, times[0], hs[0]), float)
        x2 = np.asarray(pss.solve_timestep(x1, times[1], hs[0]), float)
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
    `Pq` non-zero where the plain walk seeds it at zero, a 100% error for
    `trap`.  Both times the TRAJECTORY matched to zero, so nothing but a
    direct comparison against the dense sensitivity would have found them,
    and `euler` was exact under both -- a one-method test would have passed.

    ⚠ What this actually guards is the OPENING pair.  Inside the loop the
    coefficients are constant for every method in the tree, so a mutation
    swapping them for a post-run snapshot does NOT fail this -- checked.
    They are still stored per step, against a future variable-order method.
    """
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'euler', 'gear'):
        pss = PSS(_rc_ladder(12), method=method, reltol=1e-6)
        with quiet(AccuracyWarning, ConvergenceWarning):
            pss.solve(period=1e-3, timestep=1e-3 / 50, maxiterations=2)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(1e-3, 50, None)
        rng = np.random.default_rng(3)
        x_in = 0.01 * rng.standard_normal(m)
        with quiet():
            wk = pss._walk('plain', x_in, times, hs, T=1e-3, want_dT=True)
            _x0, _xe, Mx, Mt = (wk.x0, wk.x_end, wk.monodromy(),
                                wk.period_column())
            _wf = pss._walk('plain', x_in, times, hs, T=1e-3, dense=False,
                            keep=True, want_dT=True)
            opening, steps, x0f, xef, Mtf = (_wf.opening, _wf.steps, _wf.x0,
                                             _wf.x_end, _wf.period_column())

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
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'euler'):
        out = []
        for mf in (False, True):
            pss = PSS(_q20_rlc(), method=method, reltol=1e-6)
            with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    out = []
    for mf in (False, True):
        pss = PSS(_phase_circuit(), method='trap', reltol=1e-8)
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    pss = PSS(_phase_circuit(), method='gear', reltol=1e-8)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=25)
    m = pss.cir.n - 1
    T = pss.period
    times, hs = pss._period_grid(T, 100, None)
    rng = np.random.default_rng(11)
    a, b = 0.01 * rng.standard_normal(m), 0.01 * rng.standard_normal(m)
    with quiet():
        wk = pss._walk('pair', np.concatenate((a, b)), times, hs, T=T,
                       want_dT=True)
        Pl, Pp, Ptl, Ptp = wk.P[0], wk.P[1], wk.Pt[0], wk.Pt[1]
        _wf = pss._walk('pair', np.concatenate((a, b)), times, hs, T=T,
                        dense=False, keep=True, want_dT=True)
        C0, st, Ptlf, Ptpf = _wf.opening, _wf.steps, _wf.Pt[0], _wf.Pt[1]

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
        with quiet(AccuracyWarning):
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
    from pycircuit.circuit.linearsolver import AutoSolver, SuperLUSolver
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

    ## and the default is the analyses' own: `AutoSolver` since 2026-10-02
    ## (dense below 250 unknowns, the historical path; `DenseSolver` until then)
    plain = PSS(_q20_rlc(), method='trap', reltol=1e-6)._transient()
    assert isinstance(plain._get_linearsolver(), AutoSolver)
    assert isinstance(plain._get_nrsolver(), StandardNewton)


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
    with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    import types
    cap = {}

    def grab(self, func, z0, *args, **kwargs):
        cap['func'] = func
        cap['z0'] = np.asarray(z0, dtype=float).copy()
        return PSS._free_period_solve(self, func, z0, *args, **kwargs)

    got = {}
    for method in ('gear', 'trap'):
        pss = PSS(_phase_circuit(), method=method, reltol=1e-9)
        pss._free_period_solve = types.MethodType(grab, pss)
        with quiet(AccuracyWarning):
            pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
        func, z = cap['func'], cap['z0'].copy()
        with quiet():
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
    with quiet(AccuracyWarning):
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
    with quiet():
        pss.solve(period=1e-3, timestep=1e-5, maxiterations=20)
    assert pss.converged


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
        with quiet(AccuracyWarning):
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


def test_the_returned_waveform_closes_on_a_non_uniform_grid():
    """The replay must walk the same `(t, h)` pairs the traversal did.

    The two walks pair them differently: gear's pair walk (`_walk('pair',
    ...)`) walks `times[1:]` with `hs[_j]`, while the plain walk takes the
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
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    cir, per = _resonator_at_resonance()

    got = {}
    for method in ('trap', 'euler'):
        for npts in (100, 101, 200):
            for flag in (False, True):
                pss = PSS(cir, method=method, reltol=1e-10)
                with quiet(AccuracyWarning):
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
        with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    for method in ('trap', 'gear'):
        gaps = []
        for npts in (100, 200, 400):
            pss = PSS(_phase_circuit(), method=method, reltol=1e-10)
            with quiet(AccuracyWarning):
                res = pss.solve(period=1e-3, timestep=1e-3 / npts,
                                maxiterations=30)
            assert pss.converged
            T = pss.period
            Xw = np.asarray(res['tpss'].x, dtype=float)
            ir = pss.irefnode
            red = lambda col: np.concatenate((Xw[:ir, col], Xw[ir + 1:, col]))
            x0 = red(0)
            times, hs = pss._period_grid(T, npts, None)
            with quiet():
                if pss._solves_history():
                    Mt = pss._walk('pair', np.concatenate((x0, red(-2))),
                                   times, hs, T=T, want_dT=True).Pt[0]
                else:
                    Mt = pss._walk('plain', x0, times, hs, T=T,
                                   want_dT=True).Pt[0]
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
                with quiet(AccuracyWarning):
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

        _an.fsolve = counting
        try:
            pss = PSS(cir, method='trap', reltol=1e-10)
            with quiet(AccuracyWarning):
                pss.solve(period=per, timestep=per / 200, maxiterations=40)
            assert pss.converged
            counts[tag] = n[0]
        finally:
            _an.fsolve = orig

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

    _an.fsolve = recording
    try:
        for cf, per_, method in ((lambda: _resonator_at_resonance()[0],
                                  _resonator_at_resonance()[1], 'trap'),
                                 (_phase_circuit, 1e-3, 'trap')):
            pss = PSS(cf(), method=method, reltol=1e-8)
            with quiet(AccuracyWarning):
                pss.solve(period=per_, timestep=per_ / 100, maxiterations=30)
    finally:
        _an.fsolve = orig

    assert seen and all(seen),         'PSS called fsolve with line_search=%r -- the damping is '         'implemented and not asked for' % (seen,)


def test_the_phase_pin_is_reselected_every_iterate_and_rescues_a_far_seed():
    """B3, built as `phase_rule='reselect'` (2026-09-16): Aprille & Trick's
    Step 3 -- at every iterate pin `k = argmax |dphi/dT|` at the iterate's
    OWN value -- against the frozen seed pin it replaces.

    Van der Pol at mu = 1 from 4x the orbit amplitude: a frozen pin names a
    value the orbit never attains and the solve fails (0/6 seeds measured,
    trap/radau/gear alike); re-selected, it reaches the on-orbit period
    (6/6).  ⚠ The substitution A&T write it as is NOT the gain -- with `k`
    frozen it is the bordered solve to <= 1e-9 -- so the frozen control is
    the old rule, not a strawman.  Matrix-free keeps less of the gain (4/6);
    the seed below is one it solves, and a DENSE frozen solve re-seeded at
    its answer must reproduce the period, so the matrix-free row is the same
    orbit's, not a lucky collapse.

    ⚠⚠ AND THE DEFAULT IS ASSERTED HERE BECAUSE RE-SELECTION IS NOT FREE.
    With it as the default the full suite failed six tests (2026-09-16): the
    grid-aligned `Idtmod` wrap stops converging, and re-selection lands on a
    different PHASE of the same orbit, which moves every phase-sensitive
    surface (`frequency_aware_ppv` mode content at 0.1 f0, 1.64e-6 against
    2.45e-6) and breaks the dense-vs-matrix-free `lambda_2` bit-equality.
    So it ships opt-in, and a future flip of the default has to fail this.
    """
    circuit.default_toolkit = circuit.numeric

    def vdp():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: u - u ** 3 / 3.0)
        return c

    def run(seed, rule, **kw):
        pss = PSS(vdp(), method='trap', reltol=1e-9)
        try:
            with quiet(AccuracyWarning, ConvergenceWarning):
                pss.solve(period=6.3, timestep=6.3 / 200, x0=seed,
                          maxiterations=40, phase_rule=rule, **kw)
            return pss, pss.converged
        except Exception:                                 # noqa: BLE001
            return pss, False

    ## on the orbit the two rules agree on the answer
    on = 2.0 * np.array([np.cos(0.3), np.sin(0.3)])
    p_f, ok_f = run(on, 'frozen')
    p_r, ok_r = run(on, 'reselect')
    assert ok_f and ok_r
    T_on = float(p_r.period)
    assert abs(T_on - 6.6633) < 5e-3
    assert abs(float(p_f.period) - T_on) < 1e-9 * T_on, (p_f.period, T_on)
    assert p_r.phase_rule == 'reselect' and p_f.phase_rule == 'frozen'

    ## ⚠ THE DEFAULT IS THE FROZEN RULE, and this is where a flip of it has
    ## to fail: six suite tests broke under a re-selecting default.
    p_d = PSS(vdp(), method='trap', reltol=1e-9)
    with quiet(AccuracyWarning):
        p_d.solve(period=6.3, timestep=6.3 / 200, x0=on, maxiterations=40)
    assert p_d.phase_rule == 'frozen', \
        'the autonomous phase rule now defaults to %r; with re-selection as ' \
        'the default the full suite failed six tests (the grid-aligned ' \
        'Idtmod wrap stops converging, and phase-sensitive PPV surfaces ' \
        'move), so this is opt-in on purpose' % (p_d.phase_rule,)

    ## from 4x, the frozen pin fails and re-selection reaches the orbit
    for ang in (0.3, 1.35):
        far = 8.0 * np.array([np.cos(ang), np.sin(ang)])
        p_f, ok_f = run(far, 'frozen')
        assert not (ok_f and abs(float(p_f.period) - T_on) < 1e-6 * T_on), \
            'the frozen pin now solves the 4x seed at angle %g, so this ' \
            'fixture no longer separates the rules' % ang
        p_r, ok_r = run(far, 'reselect')
        assert ok_r and abs(float(p_r.period) - T_on) < 1e-9 * T_on, \
            're-selection did not reach the orbit from 4x at angle %g: ' \
            'converged=%r, T=%r against %r' % (ang, ok_r, p_r.period, T_on)

    ## matrix-free: the builders' row re-selects too
    far = 8.0 * np.array([np.cos(0.3), np.sin(0.3)])
    p_m, ok_m = run(far, 'frozen', x0_unknown=True, matrix_free=True)
    assert not (ok_m and abs(float(p_m.period) - 6.6635) < 1e-3)
    p_m, ok_m = run(far, 'reselect', x0_unknown=True, matrix_free=True)
    assert ok_m and abs(float(p_m.period) - 6.6635) < 1e-3, (ok_m, p_m.period)
    x_start = np.asarray(p_m._period_state[1], dtype=float).ravel()
    p_d, ok_d = run(x_start, 'frozen', x0_unknown=True)
    assert ok_d and abs(float(p_d.period) - float(p_m.period)) < 1e-12, \
        (p_d.period, p_m.period)

    with pytest.raises(ValueError, match='phase_rule'):
        PSS(vdp(), method='trap').solve(period=6.3, timestep=6.3 / 50,
                                        x0=on, phase_rule='moving')


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
            with quiet(AccuracyWarning):
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
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def run(tstab):
        pss = PSS(_q20_rlc(), method='trap', reltol=1e-9)
        with quiet(AccuracyWarning):
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


@pytest.mark.parametrize('method', ['trbdf2', 'radau', 'glm2', 'glm3'])
def test_forward_replay_is_the_exact_transpose_of_the_adjoint(method):
    """``<xa, W u> == <W^T xa, u>`` to machine precision for the stage methods
    and the Nordsieck GLMs' maps on the state (`PSS._state_map`: the source
    enters every stage, the output rows and the opening startup's
    substages -- `_GLMStateStep`; 2026-09-25).

    The forward driven replay ``_forced_replay`` and the adjoint
    ``_forced_replay_transposed`` are built from the SAME per-step source
    coupling, so they must be exact transposes -- the step-level identity that
    pins the sign convention (a flipped source term negates the whole driven
    response, which this catches while an end-to-end magnitude check might
    not).  Checked on a converting diode mixer, where every abscissa carries a
    non-trivial coupling.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _diode_mixer()
    pss = PSS(cir, method=method, reltol=1e-11)
    with quiet():
        pss.solve(period=1e-6, timestep=1e-6 / 160, maxiterations=40)
    fp = pss._state_map()
    m = cir.n - 1
    rng = np.random.default_rng(1)
    u = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    xa = rng.standard_normal(m) + 1j * rng.standard_normal(m)
    Wu, _ = pss._forced_replay(fp, 3e5, u, y0=np.zeros(m))
    WTxa = pss._forced_replay_transposed(fp, 3e5, xa)
    lhs = xa @ Wu
    rhs = WTxa @ u
    assert abs(lhs - rhs) / abs(lhs) < 1e-12, \
        '%s forward replay is not the transpose of the adjoint: %.2e' \
        % (method, abs(lhs - rhs) / abs(lhs))


def test_the_forced_replay_superposes():
    """`y_end = M y0 + w` — the property the whole `m x m` reduction rests on.

    PAC solves an `m x m` system instead of an `(N m) x (N m)` one because
    the driven replay is LINEAR in its initial state and its source
    separately. If that ever stopped holding, `(I - alpha M) y_0 = alpha w`
    would be solving the wrong equation, and the answer would still look
    entirely plausible.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
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


def test_the_matvecs_take_a_complex_vector():
    """PAC needs `I + alpha(f) H` with `alpha` complex; `M` is real.

    So a complex product is TWO REAL REPLAYS against the same stored
    factors, exactly — not a complex refactorisation, which would double
    the stored factors for a map with no imaginary part. The three matvecs
    used to cast with `dtype=float`, which does not refuse a complex vector,
    it DISCARDS its imaginary half.
    """
    circuit.default_toolkit = circuit.numeric
    cir = _pac_circuit()
    per = 1e-3
    pss = PSS(cir, method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
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
    with quiet(AccuracyWarning):
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


def test_saltation_is_unneeded_for_a_discontinuous_INJECTION_too():
    """⚠⚠ C5 GENERALISES, AND THAT REMOVES A6'S REASON FOR SALTATION.

    C5 measured a switched *conductance* (`VSwitch`) and found the
    monodromy-vs-finite-difference gap falling at exactly **2.00× per
    doubling** — O(h) discretisation, not the O(1) a missing saltation
    term leaves. The reason it gave is general: *"each step uses its own
    converged `Jf` and `C`, which already describe whichever side of the
    switch that step is on. The saltation matrix is a CONTINUOUS-time
    construct … a discrete map has no instant at which the field is
    undefined."*

    A PFD is not a switched conductance — it is a discontinuous **current
    injection**, `K·sign(v_ctl − V_th)`, whose Jacobian is *zero* either
    side of the crossing and undefined at it. So C5's result had to be
    re-run with the element changed and everything else held, and it
    reproduces:

        element                    200        400          800
        VSwitch (conductance)   9.74e-04   4.86e-04   2.43e-04   2.01x 2.00x
        sign source (injection) 7.54e-05   3.76e-05   1.88e-05   2.01x 2.00x

    Both toggle twice per period with `|M|` = 0.653 and 0.990, so neither
    is the numerical-zero trap C5's own docstring records falling into.

    ⚠⚠ SO A6'S SALTATION CLAIM IS NOW WRONG TWICE OVER. Its stated reason
    — a locked orbit sitting at zero phase error, leaving the loop open —
    describes an unmitigated dead zone, which real designs avoid by
    offsetting the charge pump. And the replacement reason — that the
    field is discontinuous at the switching instants — is falsified here,
    measured, at exactly first order.

    ⚠ WHAT WOULD STILL NEED IT, recorded as the open question rather than
    a finding: C5's argument turns on the discrete map having no undefined
    instant, which holds when the discontinuity is in the ALGEBRAIC part
    — the current or the conductance — because the per-step Newton
    resolves it. It would not obviously hold for a discontinuity in the
    STATE: a genuine reset of `x`, which is what a divider or counter
    rollover is. **So the hard part of a PLL may be the divider rather
    than the PFD**, and that is the next thing to measure rather than
    assume.
    """
    from pycircuit.circuit import VSwitch
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def build(kind):
        c = SubCircuit()
        for nn in ('ctl', 'a', 'b'):
            c.add_node(nn)
        c['vc'] = VSin('ctl', gnd, va=2.0, freq=1.0 / per)
        c['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / per, phase=90.0)
        c['rs'] = R('a', 'b', r=1e5)
        c['c'] = C('b', gnd, c=1e-6)
        if kind == 'vswitch':
            c['sw'] = VSwitch('b', gnd, 'ctl', gnd, Ron=1e3, Roff=1e9,
                              Von=1.0, Voff=0.0)
        else:
            c['sw'] = BSource('b', gnd, 'ctl', gnd,
                              i_func=lambda u: 2e-5 * np.sign(u - 0.5))
        return c

    def gap(kind, npts):
        pss = PSS(build(kind), method='trap', reltol=1e-11)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(per, npts, None)
        with quiet(AccuracyWarning):
            res = pss.solve(period=per, timestep=per / npts, maxiterations=60)
        assert pss.converged
        ir = pss.irefnode
        Xw = np.asarray(res['tpss'].x, dtype=float)
        ctl = Xw[pss.cir.get_node_index('ctl')]
        toggles = int(np.sum(np.diff((ctl > 0.5).astype(int)) != 0))
        assert toggles >= 2, \
            '%s no longer toggles (%d), so the test has lost its subject' \
            % (kind, toggles)
        x0 = np.concatenate((Xw[:ir, 0], Xw[ir + 1:, 0]))

        def phi(v):
            with quiet():
                xe = pss._walk('plain', np.asarray(v, dtype=float), times,
                               hs, T=per).x_end
            return np.asarray(xe, dtype=float)

        with quiet():
            M = pss._walk('plain', x0, times, hs, T=per).monodromy()
        M = np.asarray(M, dtype=float)
        assert np.linalg.norm(M) > 0.1, \
            '%s erases its state (|M| = %.3e); comparing zero against ' \
            'zero reports a meaningless agreement' \
            % (kind, np.linalg.norm(M))
        base = phi(x0)
        Mfd = np.zeros((m, m))
        eps = 1e-7
        for j in range(m):
            d = np.zeros(m)
            d[j] = eps
            Mfd[:, j] = (phi(x0 + d) - base) / eps
        return np.max(np.abs(M - Mfd)) / max(np.max(np.abs(Mfd)), 1e-300)

    for kind in ('vswitch', 'signsrc'):
        gaps = [gap(kind, n) for n in (200, 400, 800)]
        for a, b in zip(gaps, gaps[1:]):
            assert 1.7 < a / b < 2.3, \
                '%s: the gap falls %.2fx per doubling, not ~2. A rate near ' \
                '1 is an O(1) structural term -- i.e. saltation IS needed ' \
                'here -- and that would overturn C5. Gaps: %s' \
                % (kind, a / b, gaps)


def test_a_state_reset_needs_no_saltation_but_grid_alignment_is_a_cliff():
    """⚠⚠ THE DIVIDER QUESTION, ANSWERED — and the answer is not saltation.

    A PFD's discontinuity is in the ALGEBRAIC part, which the per-step
    Newton resolves; C5's argument covers it and the companion test above
    measures it. A divider is different: a counter **resets its state**.
    `Idtmod` is that object — and contrary to a note in this record, the
    wrap folds the STATE, not only the output map (`I.idt_node` stays
    inside one modulus and is periodic to 2.6e-15).

    ⚠ OFF A GRID POINT, THE FLOW MAP IS DIFFERENTIABLE AND SALTATION IS
    UNNECESSARY. Variational monodromy against finite differences, wrap
    placed off-grid by `ic = 0.31`:

        npts    rel err     rate    FD noise floor
         250   2.754e-06            7.7e-08
         500   1.349e-06   2.04x    1.3e-07
        1000   6.738e-07   2.00x    2.7e-07
        2000   5.551e-07   1.21x    5.5e-07   <- floor reached

    Exactly first order, same as the switched conductance and the
    discontinuous injection. The last row's 1.21x is the FD instrument's
    own noise meeting the signal, not an O(1) term — which is why the
    eps-stability is measured alongside rather than assumed away.

    ⚠⚠ BUT ON A GRID POINT THE MAP IS GENUINELY DISCONTINUOUS, AND THE
    PSS CONVERGES ANYWAY. With the wrap landing exactly on a grid point,
    `|Δφ|/ε` does not settle to a derivative — it scales as `1/ε`:

        eps       npts=600 (off)   npts=1200 (ON a grid point)
        1e-10     1.732085         8.202463e+07
        1e-08     1.732026         8.202453e+05
        1e-06     1.732025         8.201463e+03
        1e-05     1.732025         8.192475e+02

    Constant across six decades on the left; `|Δφ|` a CONSTANT ≈8.2e-3
    independent of `ε` on the right. **A perturbation of any size produces
    the same finite jump**, because an infinitesimal change flips which
    step the reset lands in and a whole modulus propagates.

    ⚠ SO A6'S REAL PROBLEM IS EVENT LOCALISATION, NOT SALTATION. The
    machinery exists — `_WrapEvents.next_event` predicts the crossing —
    and the question is whether a fixed-grid traversal uses it. Recorded
    as the finding rather than fixed here, because "the monodromy is
    wrong when the reset is grid-aligned" and "the PSS reports
    convergence there" are two separate defects and the second is the
    dangerous one.
    """
    circuit.default_toolkit = circuit.numeric
    per = 1e-3

    def build(ic):
        c = SubCircuit()
        for nn in ('in', 'out', 'f'):
            c.add_node(nn)
        c['vin'] = VS('in', gnd, v=2000.0)
        c['I'] = Idtmod('in', gnd, 'out', gnd, modulus=1.0, ic=ic)
        c['Rf'] = R('out', 'f', r=1e3)
        c['Cf'] = C('f', gnd, c=1e-7)
        c['Rl'] = R('out', gnd, r=1e5)
        return c

    def pieces(npts, ic):
        pss = PSS(build(ic), method='trap', reltol=1e-11)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(per, npts, None)
        with quiet(AccuracyWarning, ConvergenceWarning):
            ## (`npts` POINTS, the walk's `_period_grid` above: npts - 1 steps)
            res = pss.solve(period=per, timestep=per / (npts - 1), maxiterations=60)
        assert pss.converged
        ir = pss.irefnode
        Xw = np.asarray(res['tpss'].x, dtype=float)
        x0 = np.concatenate((Xw[:ir, 0], Xw[ir + 1:, 0]))

        def phi(v):
            with quiet():
                xe = pss._walk('plain', np.asarray(v, dtype=float), times,
                               hs, T=per).x_end
            return np.asarray(xe, dtype=float)

        with quiet():
            M = pss._walk('plain', x0, times, hs, T=per).monodromy()
        return np.asarray(M, dtype=float), phi, phi(x0), m

    ## the state really does reset -- otherwise this tests nothing new
    M0, _phi0, _base0, m0 = pieces(500, 0.31)
    assert m0 >= 4 and np.linalg.norm(M0) > 0.1, \
        'the monodromy is %.3e; a circuit that erases its state has ' \
        'nothing to check' % np.linalg.norm(M0)

    def gap(npts, ic, eps=1e-7):
        pss = PSS(build(ic), method='trap', reltol=1e-11)
        m = pss.cir.n - 1
        times, hs = pss._period_grid(per, npts, None)
        with quiet(AccuracyWarning, ConvergenceWarning):
            ## (`npts` POINTS, the walk's `_period_grid` above: npts - 1 steps)
            res = pss.solve(period=per, timestep=per / (npts - 1), maxiterations=60)
        assert pss.converged
        ir = pss.irefnode
        Xw = np.asarray(res['tpss'].x, dtype=float)
        x0 = np.concatenate((Xw[:ir, 0], Xw[ir + 1:, 0]))

        def phi(v):
            with quiet():
                xe = pss._walk('plain', np.asarray(v, dtype=float), times,
                               hs, T=per).x_end
            return np.asarray(xe, dtype=float)

        with quiet():
            M = pss._walk('plain', x0, times, hs, T=per).monodromy()
        M = np.asarray(M, dtype=float)
        base = phi(x0)
        Mfd = np.zeros((m, m))
        for j in range(m):
            d = np.zeros(m)
            d[j] = eps
            Mfd[:, j] = (phi(x0 + d) - base) / eps
        return (np.max(np.abs(M - Mfd)) / max(np.max(np.abs(Mfd)), 1e-300),
                x0, phi, base, m)

    ## OFF-GRID: a genuine derivative, and the gap falls at O(h)
    g = [gap(n, 0.31)[0] for n in (250, 500)]
    assert 1.7 < g[0] / g[1] < 2.4, \
        'off-grid the gap falls %.2fx per doubling, not ~2. A rate near 1 ' \
        'would be an O(1) term -- saltation genuinely needed for a state ' \
        'reset -- which would overturn this result. Gaps: %s' \
        % (g[0] / g[1], g)

    ## ON-GRID: not a derivative at all -- |dphi| is constant in eps
    _rel, x0g, phig, baseg, mg = gap(1200, 0.0)
    norms = []
    for eps in (1e-9, 1e-7, 1e-5):
        d = np.zeros(mg)
        d[3] = eps
        norms.append(float(np.linalg.norm(phig(x0g + d) - baseg)))
    spread = max(norms) / min(norms)
    assert spread < 10.0, \
        'on a grid-aligned reset |dphi| should be ~constant in eps (a ' \
        'DISCONTINUITY, not a derivative); it spread %.1fx over four ' \
        'decades of eps: %s' % (spread, norms)


def test_the_closing_step_period_column_matches_its_own_derivative():
    """⚠ The `'closing'` period column, pinned to what it CLAIMS and no more.

    `_period_column = 'closing'` implements the convention commercial PSS
    engines use: the inner transient owns the steps and the LAST one is
    placed on the period boundary, so `dh_i/dT = 0` inside and
    `dh_N/dT = 1` at the close. It uses `residual_dh`, THE PARTIAL, rather
    than `residual_dT` -- the total is 3/2 of the partial for Gear-2
    precisely because it assumes every step scales, and here only one step's
    `h` moves.

    ⚠⚠ THIS TEST DOES NOT CLAIM THE CONVENTION IS BETTER. An earlier
    measurement said the shipped proportional column was `O(h)` wrong; that
    was RETRACTED -- the `O(h)` was in the finite-difference instrument, and
    the shipped column converges at `O(h^2)` (roadmap section D shape 0j).
    What is asserted here is only that the analytic `'closing'` column
    agrees with a finite difference of the SAME convention, and that
    selecting it does not disturb the default.
    """
    circuit.default_toolkit = circuit.numeric
    Q, a = 8.0, 0.25
    mu = 1.0 / (2 * np.pi * Q)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0)
                       + a * mu * (u ** 2 - 2.0))
    cir.add_node('x')
    cir['L'] = L('v', 'x', L=1.0)
    cir['Rs'] = R('x', gnd, r=0.2 * mu)
    pss = PSS(cir, method='trap', reltol=1e-12)
    x0 = np.zeros(cir.n - 1)
    x0[0] = 2.0
    with quiet(AccuracyWarning):
        pss.solve(period=2 * np.pi, timestep=2 * np.pi / 240, x0=x0,
                  maxiterations=250)
    assert pss.converged
    _s, x0s, _xm1, times, hs, T, _xu = pss._period_state
    xin = np.asarray(x0s, dtype=float).ravel()
    assert pss._period_column == 'proportional', \
        'the default convention must be unchanged'

    def endpoint(Tv, t_, h_):
        with quiet():
            o = pss._walk('plain', xin, t_, h_, T=Tv)
        return np.asarray(o.x_end, dtype=float)

    def analytic(mode):
        pss._period_column = mode
        try:
            with quiet():
                o = pss._walk('plain', xin, times, hs, T=T, want_dT=True)
        finally:
            pss._period_column = 'proportional'
        return np.asarray(o.period_column(), dtype=float).ravel()

    ## the finite difference of the SAME convention: only the closing step
    ## lengthens, and it is hs[len(times) - 2] reaching times[-1]
    base = endpoint(T, times, hs)
    dT = T * 1e-7
    tl = np.array(times, dtype=float).copy()
    tl[-1] = T + dT
    hl = np.array(hs, dtype=float).copy()
    hl[len(times) - 2] = hl[len(times) - 2] + dT
    fd_close = (endpoint(T + dT, tl, hl) - base) / dT

    ana_close = analytic('closing')
    scale = max(float(np.linalg.norm(fd_close)), 1e-30)
    assert np.linalg.norm(ana_close - fd_close) / scale < 1e-5, \
        'the analytic closing column must match a finite difference of the ' \
        'same convention; got %.3e relative' \
        % (np.linalg.norm(ana_close - fd_close) / scale)

    ## and the two conventions must actually DIFFER, or the flag does nothing
    ana_prop = analytic('proportional')
    assert np.linalg.norm(ana_prop - ana_close) / scale > 1e-9, \
        'the two conventions produced the same column, so the flag is inert'


def test_the_manufactured_opening_step_is_INCONSISTENT_on_an_L_I_cutset():
    """⚠⚠ THE DEFAULT METHOD RETURNS EXACTLY 2x, AND ON HALF THE GRIDS IT
    REPORTS SUCCESS WHILE DOING SO.

    An `ISin` forcing an inductor to ground is an L-I cutset with a CLOSED
    FORM: the source fixes `i = I sin(wt)`, so `|v1| = w L I` needs no
    integrator. Measured (`timestep = per/N` gives N points and N-1 STEPS):

        method  pts  steps  parity  converged   |v1|/exact   |v(T)-v(0)|/V
        trap    100    99   odd      False        2.000672      2.0e+00
        trap    101   100   even     TRUE         2.000658      2.7e-13
        trap    400   399   odd      False        2.000041      2.0e+00
        trap    401   400   even     TRUE         2.000041      1.3e-12

    ⚠⚠⚠ ON AN EVEN NUMBER OF STEPS THE 2x IS CONVERGED **AND** PERIODIC TO
    1e-13. The flag does not merely fail to discriminate -- it AFFIRMS the
    wrong answer. Whether a caller is warned depends on the parity of the
    point count, which nobody would think to vary. An earlier version of this
    record said "converged reports False, so this is not silent"; that was
    measured on odd-step grids only and is WRONG.

    ⚠ THE MECHANISM IS AN INCONSISTENT INITIAL VALUE, NOT A BAD MODE. The
    samples are exactly `v_n = V (cos(w t_n) - (-1)^n)`: the smooth part is
    RIGHT and a unit ripple rides on it, so peak-to-peak doubles. The index-2
    constraint FIXES `v(0) = wLI`, but the plain path MANUFACTURES `x(0)`
    from the entering state with one order-dropped step and starts the
    algebraic variable at ZERO -- an error of exactly `-V`. What each method
    then does with that seed is its stability function at the algebraic limit
    `|sh| -> inf`: trap `(2+sh)/(2-sh) -> -1` carries it forever at constant
    amplitude, Euler `1/(1-sh) -> 0` kills it in one step, Gear-2 `-> 1/3` in
    a few. That predicts the whole table, and predicts the ripple amplitude
    is EXACTLY `V` rather than an arbitrary null-space coefficient.

    ⚠ SO EULER'S `False` IS HONEST, not a false alarm: it damps the seed
    within the period, so its endpoints genuinely differ by that one seed and
    `|v(T)-v(0)|/V = 1.0` exactly.

    ✅ AND `x0_unknown=True` FIXES IT AT EVERY PARITY, because it makes
    `x(0)` a genuine unknown instead of manufacturing it -- which is why
    Gear-2, whose solved-history path already solves for `x(0)`, was never
    affected. This is a second and stronger reason for B1.
    """
    from pycircuit.circuit.shooting import topological_index
    from pycircuit.circuit.elements import ISin
    circuit.default_toolkit = circuit.numeric
    freq, ia, ll = 1.0, 1.0, 1e-3
    per = 1.0 / freq
    exact = 2 * np.pi * freq * ll * ia

    def build():
        c = SubCircuit()
        c.add_node('1')
        c['is'] = ISin(gnd, '1', ia=ia, freq=freq)
        c['l'] = L('1', gnd, L=ll)
        return c

    ## the netlist alone says L-I cutset, before any solve -- which is what
    ## makes the risk stateable in advance
    idx, info = topological_index(build())
    assert idx == 2 and set(info['cutset']) == {'is', 'l'}
    assert not info['ill_posed']

    def run(method, npts, x0_unknown=False):
        cir = build()
        pss = PSS(cir, method=method, reltol=1e-10)
        with quiet(AccuracyWarning, ConvergenceWarning):
            ## (the grid this was measured on: `T / N` gave N - 1 steps until 2026-09-30)
            pss.solve(period=per, timestep=per / (npts - 1),
                      x0=np.zeros(cir.n - 1), maxiterations=60,
                      x0_unknown=x0_unknown)
        X = np.asarray(pss.waveform[1], dtype=float)
        row = X[[str(nd) for nd in cir.nodes].index('1')]
        return (float(np.max(np.abs(row))) / exact, pss.converged,
                abs(float(row[-1]) - float(row[0])) / exact,
                float(row[0]) / exact)

    ## ⚠ THE DANGEROUS CASE: an EVEN number of steps (401 points).
    r, conv, pres, v0 = run('trap', 401)
    assert abs(r - 2.0) < 1e-3, 'expected exactly 2x; got %.6f' % r
    assert conv, \
        'the point of this test is that the 2x is REPORTED AS CONVERGED on ' \
        'an even number of steps'
    assert pres < 1e-9, \
        'and periodic to machine precision, so a residual check passes too; ' \
        'got %.2e' % pres
    assert abs(v0) < 1e-3, \
        'the algebraic variable starts at ZERO, which is inconsistent -- the ' \
        'constraint fixes it at wLI; got v(0)/exact = %.5f' % v0

    ## ⚠ and the odd-step grid DOES report failure, which is why the defect
    ## hid: whether you are warned depends on the parity of the point count
    _r2, conv2, pres2, _v2 = run('trap', 400)
    assert not conv2 and pres2 > 0.5, \
        'the odd-step grid should fail loudly (%r, %.2e)' % (conv2, pres2)

    ## ✅ x0_unknown=True fixes it at BOTH parities
    for npts in (400, 401):
        rf, convf, presf, v0f = run('trap', npts, x0_unknown=True)
        assert abs(rf - 1.0) < 1e-2, \
            'x0_unknown must give the closed form at %d points; got %.6f' \
            % (npts, rf)
        assert convf and presf < 1e-9, \
            'and must converge periodically (%r, %.2e)' % (convf, presf)
        assert abs(abs(v0f) - 1.0) < 1e-2, \
            'with a CONSISTENT v(0) = wLI; got %.5f' % v0f

    ## gear was never affected -- its solved-history path solves for x(0)
    rg, convg, presg, v0g = run('gear', 401)
    assert abs(rg - 1.0) < 1e-2 and convg and presg < 1e-9
    assert abs(abs(v0g) - 1.0) < 1e-2, \
        'gear starts consistent already; got %.5f' % v0g


def test_x0_unknown_defaults_from_the_topology_and_only_where_it_is_proved():
    """⚠ `x0_unknown=None` means "decide from the topology", and the
    CONDITIONALITY is the design, not a hedge.

    ⚠⚠ IT IS NOT A NEW GLOBAL DEFAULT, because `x0_unknown` is NOT FREE.
    Trapezoidal still needs an L-stable opener, so switching it on moves the
    Euler step INSIDE the period, where it degrades the ORBIT rather than
    just the opening -- measured on a `Q = 20` resonator against its analytic
    20 V peak at 20.01273 (default) against 19.76939 (`x0_unknown`) at 100
    points. Turning it on everywhere would trade a real defect on a few
    circuits for a real regression on most.

    On an index-2 netlist the trade reverses: the manufactured opening step
    is INCONSISTENT there, and trapezoidal carries the seed forever.
    """
    from pycircuit.circuit.elements import ISin, VSin
    circuit.default_toolkit = circuit.numeric
    freq, ia, ll = 1.0, 1.0, 1e-3
    per = 1.0 / freq
    exact = 2 * np.pi * freq * ll * ia

    def cutset():
        c = SubCircuit()
        c.add_node('1')
        c['is'] = ISin(gnd, '1', ia=ia, freq=freq)
        c['l'] = L('1', gnd, L=ll)
        return c

    def solve_it(cir, npts, **kw):
        pss = PSS(cir, method='trap', reltol=1e-10)
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            pss.solve(period=per, timestep=per / npts,
                      x0=np.zeros(cir.n - 1), maxiterations=60, **kw)
        fired = any('index 2' in str(x.message) for x in w)
        return pss, fired

    ## the DEFAULT now gives the closed form and says why
    cir = cutset()
    pss, fired = solve_it(cir, 401)
    row = np.asarray(pss.waveform[1], dtype=float)[
        [str(nd) for nd in cir.nodes].index('1')]
    assert abs(float(np.max(np.abs(row))) / exact - 1.0) < 1e-2, \
        'the default must now give the closed form on an index-2 netlist'
    assert abs(abs(float(row[0]) / exact) - 1.0) < 1e-2, \
        'and a CONSISTENT v(0)'
    assert fired, 'and must say so -- a silently different formulation is ' \
                  'what makes a later measurement inexplicable'
    assert pss._open_at_x0 is True

    ## ⚠ an explicit False is HONOURED, defect and all -- the escape hatch
    ## has to actually work or the default is a trap of its own
    cir2 = cutset()
    pss2, fired2 = solve_it(cir2, 401, x0_unknown=False)
    row2 = np.asarray(pss2.waveform[1], dtype=float)[
        [str(nd) for nd in cir2.nodes].index('1')]
    assert abs(float(np.max(np.abs(row2))) / exact - 2.0) < 1e-2, \
        'x0_unknown=False must be honoured untouched'
    assert not fired2, 'and must not warn about a choice the caller made'

    ## ⚠⚠ and an INDEX-1 circuit must be left completely alone, or the
    ## measured regression above is inflicted on every ordinary netlist
    c3 = SubCircuit()
    c3.add_node('a')
    c3.add_node('b')
    c3['vs'] = VSin('a', gnd, va=1.0, freq=freq)
    c3['r'] = R('a', 'b', r=1e3)
    c3['c1'] = C('b', gnd, c=1e-4)
    pss3, fired3 = solve_it(c3, 200)
    assert not fired3 and pss3._open_at_x0 is False, \
        'an index-1 netlist must keep the shipped formulation'

    ## ⚠ and a bad `method` must still raise the ValueError the caller is
    ## owed -- the defaulting helper runs BEFORE argument validation and
    ## must never change which exception an invalid call produces
    c4 = cutset()
    p4 = PSS(c4, method='bogus')
    with pytest.raises(ValueError):
        p4.solve(period=per, timestep=per / 50, x0=np.zeros(c4.n - 1))


def _resonant_driven(npts, method, Lv=1e-3, Cv=1e-9, Rs=10.0):
    """A DRIVEN RLC, lightly damped and driven ON resonance.

    ⚠ THE DAMPING IS THE POINT OF THE FIXTURE, and the first version of it
    got this wrong in a way that PASSED. With `Rs = 1k` against this `L`
    and `C` the quality factor is `~1e-3`, so the monodromy decays to
    numerical ZERO over one period — and a transpose check then compares
    the zero matrix against the zero matrix and agrees to `1.9e-46`. The
    tell was that the number was too good: a real comparison of a real
    quantity does not come back at `1e-46`. At `Rs = 10` the dominant
    multiplier is `≈ 0.94`, so `M` is `O(1)` and the check has something
    to fail on.
    """
    from pycircuit.circuit.elements import VSin
    circuit.default_toolkit = circuit.numeric
    T = 2.0 * np.pi * np.sqrt(Lv * Cv)
    cir = SubCircuit()
    cir.add_node('a')
    cir.add_node('b')
    cir['vs'] = VSin('a', gnd, va=1.0, freq=1.0 / T)
    cir['r'] = R('a', 'b', r=Rs)
    cir['l'] = L('b', gnd, L=Lv)
    cir['c1'] = C('b', gnd, c=Cv)
    pss = PSS(cir, method=method, reltol=1e-11)
    with quiet(AccuracyWarning):
        ## the LTE advisory fires here -- this fixture is a transpose
        ## check, not an accuracy one, and the grid is deliberately coarse
        pss.solve(period=T, timestep=T / npts, x0=np.zeros(cir.n - 1),
                  maxiterations=250, x0_unknown=False)
    assert pss.converged
    return cir, pss


def test_the_plain_transposed_replay_matches_a_dense_transpose():
    """B8: `M^T v` on the PLAIN path, against `M` built column by column.

    ⚠ THIS IS A FROM-SCRATCH ADJOINT DERIVATION, and the last one in this
    file came out SIGN-INVERTED (roadmap §0h) because the derivative in
    Demir's eq (24) acts on `y` alone and not on the product `C^T y`. So
    the reference is not another derivation: it is the FORWARD replay,
    transposed densely. If the two disagree the new recursion is wrong,
    full stop.

    ⚠ GEAR IS THE CONTROL, not a subject. It goes through the
    solved-history path that shipped long before this, so a disagreement
    on gear means the HARNESS is broken rather than the new code — which
    is the distinction that makes a green result on euler and trap worth
    anything.

    ⚠ BOTH BRANCHES ARE EXERCISED. Euler has `b = 0`, so `Pq` never
    re-enters and the recursion has ONE term; trapezoidal has `b = -1`,
    the pair `(P, Pq)` stays coupled, and the two solves collapse into one
    through a shared bracket. Testing only trap would leave the simpler
    branch unrun, and it is the branch every adjoint surface reaches
    first.
    """
    worst = {}
    for method in ('euler', 'trap', 'gear'):
        for npts in (120, 240):
            _cir, pss = _resonant_driven(npts, method)
            fp = pss.factored_period()
            n = fp.width
            Mf = np.column_stack([np.asarray(fp.matvec(e), float)
                                  for e in np.eye(n)])
            Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), float)
                                  for e in np.eye(n)])
            scale = float(np.max(np.abs(Mf)))
            ## a reference that decayed to nothing agrees with anything
            assert scale > 1e-3, \
                '%s/%d: ||M|| = %.3e, so this fixture has become ' \
                'degenerate and the comparison below is VACUOUS -- the ' \
                'damping must be light enough to leave a monodromy' \
                % (method, npts, scale)
            rel = float(np.max(np.abs(Mt - Mf.T))) / scale
            worst[(method, npts)] = rel
            assert rel < 1e-9, \
                '%s/%d: the transposed replay disagrees with the dense ' \
                'transpose of the forward one by %.3e (relative). A sign ' \
                'or an index in the reverse recursion is wrong.' \
                % (method, npts, rel)
    ## and the multipliers must actually be O(1), so the agreement above is
    ## a statement about a real map
    assert max(worst.values()) < 1e-9, worst


def test_the_plain_transposed_replay_carries_the_autonomous_multiplier():
    """The same check on an AUTONOMOUS oscillator, where `M` has `λ₁ = 1`.

    The driven fixture above is a contraction; an oscillator is not, and
    the unit multiplier along the trajectory is what every PPV and phase-
    noise surface actually reads out of the transpose. A recursion can be
    right on a decaying map and wrong on this one.

    ⚠ EULER IS ABSENT ON PURPOSE and is not a gap in coverage: at these
    grids it does not converge on a `Q = 8` van der Pol at all, which is a
    property of a first-order method on a limit cycle and not of the
    adjoint. It is covered on the driven fixture above, where it converges.
    """
    ## built inline rather than through a shared fixture: this needs the
    ## SAME circuit under two methods, and `_vdp_at_Q` pins gear
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    lam = {}
    for method in ('trap', 'gear'):
        for npts in (120, 240):
            cir = SubCircuit()
            cir.add_node('v')
            cir['C'] = C('v', gnd, c=1.0)
            cir['B'] = BSource('v', gnd, gnd, 'v',
                               i_func=lambda u: mu * (u - u ** 3 / 3.0)
                               + 0.25 * mu * (u ** 2 - 2.0))
            cir.add_node('x')
            cir['L'] = L('v', 'x', L=1.0)
            cir['Rs'] = R('x', gnd, r=0.2 * mu)
            pss = PSS(cir, method=method, reltol=1e-12)
            z = np.zeros(cir.n - 1)
            z[0] = 2.0
            with quiet(AccuracyWarning):
                pss.solve(period=2 * np.pi, timestep=2 * np.pi / npts, x0=z,
                          maxiterations=250, x0_unknown=False)
            assert pss.converged
            pss.monodromy = 'native'   # the PLAIN replay is the object under test
            fp = pss.factored_period()
            n = fp.width
            Mf = np.column_stack([np.asarray(fp.matvec(e), float)
                                  for e in np.eye(n)])
            Mt = np.column_stack([np.asarray(fp.matvec_transposed(e), float)
                                  for e in np.eye(n)])
            scale = float(np.max(np.abs(Mf)))
            rel = float(np.max(np.abs(Mt - Mf.T))) / scale
            assert rel < 1e-9, \
                '%s/%d: autonomous transpose disagrees by %.3e' \
                % (method, npts, rel)
            ## `M^T` must carry the unit multiplier too -- it is the
            ## same spectrum, and it is the one the phase surfaces read.
            ##
            ## ⚠ AGAINST THE FORWARD MAP, NOT AGAINST 1.0. The first
            ## version of this asserted `|λ| ≈ 1` to 5e-3 and FAILED at
            ## trap/120 with 0.992007 -- and the failure was the test's,
            ## not the code's: `M` itself carries that deficit, because
            ## 120 points on a `Q = 8` limit cycle is a coarse grid. What
            ## the transpose owes us is the FORWARD map's spectrum,
            ## whatever the grid made of it.
            lf = float(np.sort(np.abs(np.linalg.eigvals(Mf)))[-1])
            lt = float(np.sort(np.abs(np.linalg.eigvals(Mt)))[-1])
            lam[(method, npts)] = lf
            assert abs(lf - lt) < 1e-10, \
                '%s/%d: the transposed map has a different dominant ' \
                'multiplier from the forward one (%.12f vs %.12f) -- ' \
                'a transpose cannot change the spectrum' \
                % (method, npts, lt, lf)

    ## and the deficit is the DISCRETISATION converging, which is what
    ## licenses reading it as the grid rather than as a defect: trap is
    ## second-order, so halving `h` must quarter the distance to 1, and
    ## gear reaches it outright.
    e120 = abs(lam[('trap', 120)] - 1.0)
    e240 = abs(lam[('trap', 240)] - 1.0)
    assert 3.0 < e120 / e240 < 6.0, \
        'trap\'s unit-multiplier deficit is not converging at O(h^2) ' \
        '(%.3e then %.3e, ratio %.2f); if it has stopped converging the ' \
        'deficit is no longer the grid and this reading is wrong' \
        % (e120, e240, e120 / e240)
    assert abs(lam[('gear', 240)] - 1.0) < 1e-9, \
        'gear used to reach the unit multiplier exactly (%.12f)' \
        % lam[('gear', 240)]


def test_the_plain_transposed_replay_is_exact_for_a_multistep_companion():
    """The plain map's reverse recursion used to be derived for a ONE-STEP
    companion and REFUSED a step with a second history term (it would have
    silently dropped it).  Since 2026-09-23 every linear multistep step
    transposes through one derivation (`_LMMStep.adjoint`: trap's shared
    bracket and gear's pair are its two special cases), which carries the
    `a_2 C_{n-2}` term -- so a forged step with a NONZERO third coefficient
    is transposed exactly: ``<M^T u, w> = <u, M w>`` to rounding, with the
    forward replay taking the same term through `_step_sensitivity`."""
    _cir, pss = _resonant_driven(120, 'trap')
    fp = pss.factored_period()
    lu, C_new, alphas, b = fp.steps[3]
    forged = list(fp.steps)
    forged[3] = (lu, C_new, (alphas[0], alphas[1], 0.3 * alphas[1]), b)
    rng = np.random.default_rng(7)
    u = rng.standard_normal(fp.width)
    w = rng.standard_normal(fp.width)
    Mw = pss._monodromy_matvec_plain(fp.opening, forged, w)
    MTu = pss._monodromy_matvec_transposed_plain(fp.opening, forged, u)
    ## the forged term is live: the map moved
    assert np.linalg.norm(Mw - fp.matvec(w)) > 1e-6 * np.linalg.norm(Mw)
    lhs, rhs = float(np.dot(MTu, w)), float(np.dot(u, Mw))
    assert abs(lhs - rhs) < 1e-12 * max(abs(lhs), 1.0), (lhs, rhs)


def test_lte_grid_derives_a_frozen_nonuniform_grid_that_solve_accepts():
    """B7a: the DERIVATION side, promoted from `benchmarks/pss_lte_grid.py`.

    `_period_grid` has consumed caller-supplied step fractions since item
    5; what was missing was deriving them. `PSS.lte_grid` runs an adaptive
    transient, takes the accepted steps of one settled period, and returns
    them as fractions plus the state that starts the window — so the pair
    feeds straight back in:

        fracs, seed = pss.lte_grid(period=T)
        pss.solve(period=T, grid=fracs, x0=seed)

    ⚠ FRACTIONS, NOT TIMES. An autonomous period is an unknown, so every
    step must scale with `T` or `dh/dT = h/T` — the identity the period
    column rests on — stops holding.

    ⚠⚠ AND THIS IS FOR STIFF SMOOTH PROBLEMS, NOT EVENTS. On a wrapping
    `Idtmod` the derived grid is measurably WORSE than a uniform grid of
    the same count (max LTE 2.64e+05 against 1.67e+05 times tolerance),
    because the LTE peak sits at the RESET on every grid and no step size
    makes a discontinuity's truncation error small. That half of B7 needs
    the event time to be a Newton unknown and is filed against A6.
    """
    circuit.default_toolkit = circuit.numeric
    mu, T = 4.0, 11.0

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0))
        return cir

    cir = vdp()
    pss = PSS(cir, method='gear', reltol=1e-6)
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('v')] = 2.0
    with quiet(UsageWarning):
        fr, seed = pss.lte_grid(period=T, x0=x0, tstab=12 * T, reltol=1e-5)

    ## it is a grid: positive fractions of exactly one period
    assert fr.ndim == 1 and len(fr) > 10
    assert np.all(fr > 0.0)
    assert abs(float(fr.sum()) - 1.0) < 1e-12, \
        'the fractions must sum to exactly one period; got %.15f' \
        % float(fr.sum())

    ## ⚠ AND IT IS GENUINELY NON-UNIFORM, which is the only reason to
    ## derive one. A uniform result would mean the transient never
    ## adapted and the whole exercise bought nothing -- measured spreads
    ## of 4x here and 170x at mu=10.
    spread = float(fr.max() / fr.min())
    assert spread > 2.0, \
        'the derived grid is essentially uniform (max/min = %.2f), so ' \
        'the adaptive run contributed nothing' % spread

    ## and `solve` takes it, on a FRESH circuit, reaching the same period
    ## a uniform grid of the same count reaches
    p2 = PSS(vdp(), method='gear', reltol=1e-6)
    with quiet(AccuracyWarning):
        ## at the period `lte_grid` OBSERVED (2026-09-21): the grid is
        ## fractions of it, and this test's hint T = 11 is 8 % above the true
        ## 10.2 -- at the hint the folded grid's seam (mid-edge) meets the
        ## period error at the fastest dynamics
        p2.solve(period=pss.lte_period, grid=fr, x0=seed, maxiterations=60)
    assert p2.converged, 'the derived grid did not converge'

    p3 = PSS(vdp(), method='gear', reltol=1e-6)
    with quiet(AccuracyWarning):
        p3.solve(period=T, timestep=T / len(fr), x0=seed, maxiterations=60)
    assert p3.converged

    ## ⚠ AGAINST A FINE REFERENCE, NOT AGAINST EACH OTHER. The first
    ## version of this demanded the two periods agree to 1e-4 and failed
    ## at 2.5e-3 -- and the demand was WRONG, not merely tight: the two
    ## grids are supposed to differ, that is the entire point of deriving
    ## one. What they must both do is bracket the true period.
    ##
    ## Measured at mu = 4 with 174 steps, against a 3000-point run
    ## (gear is second order, so that reference is ~300x finer than the
    ## grids under test -- ample, and it keeps this test under 30s):
    ##     derived 1.214e-03      uniform 1.282e-03
    ## -- the derived grid is BETTER, but only by 5%, because mu = 4 is
    ## not stiff enough to separate them. The headline win (converging
    ## where the same count of uniform steps does NOT, and beating a
    ## 20000-point grid at 18x fewer points) is at mu = 100 in
    ## `benchmarks/pss_lte_grid.py`, which is far too slow for this suite.
    ## So this asserts NOT-WORSE, which is what is cheaply checkable here.
    pf = PSS(vdp(), method='gear', reltol=1e-11)
    with quiet(AccuracyWarning):
        pf.solve(period=T, timestep=T / 3000, x0=seed, maxiterations=60)
    assert pf.converged
    e_der = abs(p2.period - pf.period) / pf.period
    e_uni = abs(p3.period - pf.period) / pf.period
    assert e_der < 5e-3 and e_uni < 5e-3, \
        'neither grid is resolving the period (derived %.3e, uniform ' \
        '%.3e); the reference or the fixture has moved' % (e_der, e_uni)
    assert e_der < 1.5 * e_uni, \
        'the DERIVED grid is now materially worse than a uniform grid of ' \
        'the same step count (%.3e against %.3e). Adapting the steps is ' \
        'supposed to be free or better on a smooth stiff problem; if it ' \
        'has become a cost the derivation is picking the wrong window' \
        % (e_der, e_uni)


def test_lte_grid_refuses_what_it_cannot_derive():
    """The two ways to hand it something that is not a period.

    A bad period is refused in words rather than returning a grid for the
    wrong interval, which `solve` would accept without complaint — the
    fractions sum to 1 whatever window they came from.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: 1.0 * (u - u ** 3 / 3.0))
    pss = PSS(cir, method='gear', reltol=1e-6)
    with pytest.raises(ValueError, match='period must be positive'):
        pss.lte_grid(period=0.0)
    with pytest.raises(ValueError, match='tstab must not be negative'):
        pss.lte_grid(period=1.0, tstab=-1.0)


def test_warm_start_finds_the_linear_region_and_the_handoff_works():
    """`PSS.find_initial_solution` -- De Luca, Bolcato & Schilders (TCAS-I 2019)
    Algorithm 2: pre-integrate until the fixed-point iteration has entered its
    LINEAR region, then hand the iterate to shooting.

    Three things are pinned, and the first is a number PREDICTED BEFORE IT WAS
    MEASURED rather than read off the code:

    1. ON A LINEAR CIRCUIT THE ANSWER IS FORCED.  `phi` is affine, so
       `J_phi(x_khat) u` is EXACT for every `khat` and the linear prediction
       equals the actual shooting error to roundoff.  So the check must pass on
       the very first iteration and stop after exactly `n_iter` of them:
       `khat == 0` and `periods == n_iter`.  The paper says the same of its own
       RLC ("we expect the linear region to be found at the first iteration");
       its reported `khat = 4` is four preliminary integrations it performs for
       implementation reasons, not detection.  The circuit here IS the paper's:
       R = 1 ohm, L = 20 mH, C = 2 uF, T = 1.256 ms, 9 sin -- Q = 100, i.e. the
       slow-settling case that motivates a warm start at all.

    2. THE CRITERION MUST DISCRIMINATE, or it proves nothing.  A criterion that
       always answered "linear" would give `khat == 0` on every circuit, linear
       or not, and still pass check 1.  So on a NONLINEAR circuit the first
       check is required to FAIL -- the run resets and `khat` moves off 0 --
       which is only meaningful because the same code returns `khat == 0` on the
       linear one.

    3. THE HANDOFF MUST EARN ITS PLACE: shooting that does NOT converge from a
       cold start must converge from the returned iterate, at the same
       iteration budget.  Otherwise the detector is correct and useless.

    ⚠ NON-AUTONOMOUS ONLY -- see `find_initial_solution`.  For an autonomous
    oscillator the equilibrium is a fixed point of the period map and the map is
    linear around it, so this criterion certifies the TRIVIAL ROOT; the van der
    Pol case in `benchmarks/pss_warm_start.py` is outside the paper's scope and
    is deliberately not tested here as if it were solved.
    """
    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric
    T = 1.256e-3
    N_ITER = 7

    def rlc(diode):
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, va=9.0, freq=1.0 / T)
        c['R'] = R(1, 2, r=1.0)
        c['L'] = L(2, 3, L=2e-2)
        c['C'] = C(3, gnd, c=2e-6)
        if diode:
            c['D'] = Diode(3, gnd)
        return c

    ## (1) linear: the answer is forced by phi being affine
    ## vabstol=1e-12: the PREMISE below (a cold start that fails in 8 iterations) was
    ## established at the pre-2026-09-19 default; at 1e-6 the cold start gets there
    ## in 8 and the fixture demonstrates nothing.  Asked for by name.
    pss = PSS(rlc(False), method='euler', reltol=1e-9, vabstol=1e-12)
    with quiet():
        _x, info = pss.find_initial_solution(period=T, npts=200,
                                             n_iter=N_ITER, max_periods=60)
    assert info['found'], 'linear circuit: no linear region found at all'
    assert info['khat'] == 0 and info['periods'] == N_ITER, \
        'a linear phi makes the linear generator EXACT, so the region must be ' \
        'found immediately: expected khat=0, periods=%d, got khat=%d, ' \
        'periods=%d' % (N_ITER, info['khat'], info['periods'])

    ## (2) nonlinear: the criterion must REFUSE the first iterate, or (1) is
    ## satisfied by a criterion that never says no
    pssd = PSS(rlc(True), method='euler', reltol=1e-9, vabstol=1e-12)
    with quiet():
        xw, infod = pssd.find_initial_solution(period=T, npts=200,
                                               n_iter=N_ITER, max_periods=80)
    assert infod['found'], 'nonlinear circuit: no linear region found'
    assert not infod['history'][0]['ok'] and infod['khat'] > 0, \
        'the criterion did not discriminate: it accepted the FIRST iterate of ' \
        'a nonlinear circuit (khat=%d), so khat==0 on the linear one means ' \
        'nothing' % infod['khat']
    ## and the gap must actually collapse, not merely dip under a loose bound
    assert infod['history'][0]['gap'] > 100.0 * infod['history'][-1]['gap'], \
        'the linear-region gap did not collapse: %.3e -> %.3e' \
        % (infod['history'][0]['gap'], infod['history'][-1]['gap'])

    ## (3) the handoff has to be worth taking
    def solves_from(x0):
        p = PSS(rlc(True), method='euler', reltol=1e-9, vabstol=1e-12)
        try:
            with quiet(AccuracyWarning, ConvergenceWarning):
                p.solve(period=T, timestep=T / 200, x0=x0, maxiterations=8)
            return bool(p.converged)
        except Exception:
            return False
    assert not solves_from(None), \
        'the cold start now converges, so this circuit no longer demonstrates ' \
        'anything about the warm start -- pick a harder one'
    assert solves_from(xw), \
        'shooting did not converge from the warm-start iterate, which is the ' \
        'whole point of finding it'


def test_an_autonomous_collapse_onto_the_trivial_root_reports_not_converged():
    """⚠ THE MODULE ASSERTED THIS IN PROSE FOR TWO TURNS OF THE RECORD AND
    NOTHING ENFORCED IT.

    `_free_period_solve`'s docstring and `solve`'s both said "the collapse
    reports `converged = False`".  It did not.  `self.converged` is
    `(_ier == 1)` and nothing else, and `T = 0` is a REGULAR root of every
    autonomous shooting system -- `x0 - phi_T(x0)` vanishes identically there
    and the phase condition constrains `x0`, not the period -- so `fsolve`
    reaches it cleanly and reports SUCCESS.  Measured: Gear-2 returned
    `T = 5.42e-18` with `converged = True` on a circuit with no orbit in it.

    ⚠ The trivial-root warning fired correctly the whole time, and that is
    what let this survive: a reader who checks the documented flag rather than
    catching warnings got `True`.  A correct diagnostic beside a wrong status
    flag is worse than no diagnostic, because the flag is the machine-readable
    one.

    The fix demotes `ier` inside `_free_period_solve`, so all three autonomous
    call sites -- plain, solved-history and matrix-free -- inherit it, and so
    does any path added later.

    ⚠ WHAT IS NOT FIXED, AND CANNOT BE HERE: the collapse itself.  It is a
    property of the formulation, not of a method or a circuit, and the seed
    sweep below is the evidence -- at and above the fundamental both methods
    find the orbit, below it they do not, and no iteration count reaches a
    fundamental from below.  The remedy is the seed (or the PROBE technique,
    which widens the basin rather than removing the dependence).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    T0 = 1e-3

    def folding():
        c = SubCircuit()
        ## 0.5 modulus per T0, so the true fundamental is 2*T0.
        c['vin'] = VS('in', gnd, v=0.5 / T0)
        c['X'] = Idtmod('in', gnd, 'o', gnd, modulus=1.0, ic=0.31)
        c['Ro'] = R('o', gnd, r=1e6)
        return c

    ## (1) A seed below the fundamental collapses -- and must SAY so, in the
    ## flag as well as the warning.
    p = PSS(folding(), method='gear')
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter('always')
        p.solve(period=T0, timestep=T0 / 300)
    Tp = float(p._period_state[5])
    assert abs(Tp) < p.DEGENERATE_PERIOD_FACTOR * T0, \
        'the fixture stopped collapsing (T = %.6g) -- it no longer tests ' \
        'what it is named for' % Tp
    assert any('TRIVIAL root' in str(c.message) for c in caught), \
        'the collapse warning stopped firing'
    assert p.converged is False, \
        'a solve that collapsed onto T = %.6g reported converged = %r -- ' \
        'this is the defect: the documented flag says the non-orbit is an ' \
        'answer' % (Tp, p.converged)

    ## (2) ⚠ AND A GENUINE SOLVE MUST STILL REPORT True, or the "fix" is just
    ## a flag wired to False.
    q = PSS(folding(), method='gear')
    with quiet():
        q.solve(period=2 * T0, timestep=2 * T0 / 300)
    assert q.converged is True, \
        'the demotion leaked into a healthy autonomous solve'
    assert abs(float(q._period_state[5]) / (2 * T0) - 1.0) < 1e-6, \
        'expected the fundamental 2*T0, got %.6g' % q._period_state[5]

    ## (3) The seed sweep that shows the collapse is the FORMULATION, so the
    ## docstring's "not fixable here" is measured rather than asserted.
    found = {}
    for seed in (T0, 1.5 * T0, 2.0 * T0):
        r = PSS(folding(), method='gear')
        try:
            with quiet(ConvergenceWarning):
                r.solve(period=seed, timestep=seed / 300)
            found[seed] = (float(r._period_state[5]), bool(r.converged))
        except np.linalg.LinAlgError:
            found[seed] = (float('nan'), False)
    assert found[T0][1] is False, 'a seed at half the fundamental should not converge'
    assert found[1.5 * T0][1] is True and found[2.0 * T0][1] is True, \
        'seeds at and above the fundamental should find it: %r' % (found,)


def test_warping_estimate_reproduces_the_period_error_at_one_grid_with_no_reference():
    """B7, BUILT 2026-09-08.  `PSS.warping_estimate` -- defect correction on
    the solve's own grid, no refinement, no analytic reference -- against the
    true period error `T_h - T_ref`, radau at 800 points as the reference
    (its own error is at the 1e-10 ppm floor, four orders below the coarsest
    figure pinned here).

    Pinned, from the stack gate (2026-09-08): trap 1.0002; radau with the
    septic `IDEC_DEGREE` 1.0003; esdirk43 with its quintic 0.9998; and the
    EXACTNESS-CLASS ZERO -- radau with a CUBIC estimates 0.0001 of the true
    error, because a cubic spline lies inside a 3-stage collocation method's
    exactness class and the neighbouring problem is solved exactly.  That
    zero is a property, pinned so that a future "helpful" degree change
    announces itself: the failure it represents is a clean small number.
    """
    circuit.default_toolkit = circuit.numeric
    with quiet(AccuracyWarning, ConvergenceWarning):
        cir, _ = _a10_vdp()
        ref = PSS(cir, method='radau', reltol=1e-14)
        ref.solve(period=2 * np.pi, timestep=2 * np.pi / 800,
                  x0=np.array([2.0, 0.0]), maxiterations=400)
        T_ref = float(ref.period)

        def ratio(method, npts, degree=None):
            cir, _ = _a10_vdp()
            p = PSS(cir, method=method, reltol=1e-14)
            p.solve(period=T_ref, timestep=T_ref / npts,
                    x0=np.array([2.0, 0.0]), maxiterations=400)
            delta = float(p.period) - T_ref
            est = p.warping_estimate(degree=degree)
            assert est['autonomous'] is True
            return est['period_error'] / delta, est['degree']

        r, k = ratio('trap', 100)
        assert k == 3 and abs(r - 1.0) < 5e-3, 'trap/cubic %.4f (measured 1.0002 at 200-400 pts, 0.9996 at 100 in the prototype)' % r
        r, k = ratio('radau', 50)
        assert k == 7 and abs(r - 1.0) < 5e-3, 'radau/septic %.4f (measured 1.0003)' % r
        r, k = ratio('esdirk43', 100)
        assert k == 5 and abs(r - 1.0) < 5e-3, 'esdirk43/quintic %.4f (measured 0.9998)' % r
        r, k = ratio('radau', 50, degree=3)
        assert k == 3 and abs(r) < 1e-2, \
            'radau/CUBIC %.4f -- must be ~0 (inside the exactness class; measured 0.0001)' % r


def test_warping_estimate_refuses_a_period_reading_on_a_driven_circuit():
    """A forcing at `T` pins the period, so warping cannot present as a period
    change and the lag against the interpolant is bounded (entrained).  The
    guard is the `analysis='tran'` flag on `Circuit.u`: without it every
    source reads as DC and the first gate's driven control came back
    `autonomous=True` with a 'period error' read off an entrained lag.
    """
    from pycircuit.circuit.elements import ISin
    circuit.default_toolkit = circuit.numeric
    with quiet(AccuracyWarning, ConvergenceWarning):
        cir, mu = _a10_vdp()
        T = 2 * np.pi
        cir['inj'] = ISin('v', gnd, ia=0.2 * mu * 2.0, freq=1.0 / T)
        p = PSS(cir, method='trap', reltol=1e-12)
        p.solve(period=T, timestep=T / 100, x0=np.array([2.0, 0.0]),
                maxiterations=200, x0_unknown=False)
        est = p.warping_estimate(periods=8)
    assert est['autonomous'] is False
    assert est['period_error'] is None and est['ppm'] is None
    assert np.all(np.isfinite(est['lag'])) and len(est['lag']) == 8


def test_an_autonomous_solve_that_returns_an_equilibrium_is_not_converged():
    """2026-09-08, found while gating `warping_estimate` on an index-2
    oscillator.  `_free_period_solve` demoted the `T -> 0` root; the
    EQUILIBRIUM is a second trivial root, periodic at every T, and radau
    seeded 10 % below the fundamental landed on it -- amplitude 0.0000,
    state 1e-27, period near the seed, `converged = True`.  Pinned: the low
    seed reports False with the equilibrium warning; the right seed still
    reports True with the 2 V orbit; the same low seed under trap fails
    honestly (it did before, and must keep doing so).
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    cir, _ = _a10_vdp()
    T0 = 2 * np.pi
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        p = PSS(cir, method='radau', reltol=1e-14)
        p.solve(period=0.9 * T0, timestep=0.9 * T0 / 200, x0=np.array([2.0, 0.0]),
                maxiterations=400)
    W = np.asarray(p.waveform[1], dtype=float)
    amp = (W[cir.get_node_index('v')].max() - W[cir.get_node_index('v')].min()) / 2
    assert amp < 1e-6, 'the low seed no longer collapses; the test needs a new basin'
    assert p.converged is False, 'an equilibrium reported as a converged orbit'
    assert any('EQUILIBRIUM' in str(r.message) for r in rec), \
        'the demotion must say why: %r' % [str(r.message)[:60] for r in rec]
    with quiet():
        cir, _ = _a10_vdp()
        p = PSS(cir, method='radau', reltol=1e-14)
        p.solve(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]), maxiterations=400)
    W = np.asarray(p.waveform[1], dtype=float)
    amp = (W[cir.get_node_index('v')].max() - W[cir.get_node_index('v')].min()) / 2
    assert p.converged is True and abs(amp - 2.0) < 1e-2, (p.converged, amp)


def _current_sense_relaxation_oscillator(k=20.0):
    """The circuit class that cannot be integrated undamped at ANY grid: a
    comparator whose input is a branch current through a capacitor.  The
    sensed current is the companion difference quotient of the ramp voltage,
    stage sensitivity `tau/h`, undamped basin `0.94 h/(k tau)` (READING-LOG
    2.156).  Op-amp-style relaxation oscillator, comparator gain `k`, output
    time constant tau/20, the current sensed through a CCVS whose 0 V input
    branch pins the far end of a small capacitor from the ramp node.
    """
    R_, C_ = 1e3, 1e-6
    c1 = 0.01 * C_
    Ctot = C_ + c1
    tau = R_ * Ctot
    Ro = 1.0
    Co = (tau / 20) / Ro
    Rs = R_ * Ctot / c1
    cir = SubCircuit()
    for n in ('vo', 'p', 'c', 'd', 'b', 's'):
        cir.add_node(n)
    cir['Ro'] = R('vo', gnd, r=Ro)
    cir['Co'] = C('vo', gnd, c=Co)
    cir['cmp'] = BSource('d', gnd, gnd, 'vo',
                         i_func=lambda u: (1.0 / Ro) * np.tanh(k * u))
    cir['R1'] = R('vo', 'p', r=1e3)
    cir['R2'] = R('p', gnd, r=1e3)
    cir['Rc'] = R('vo', 'c', r=R_)
    cir['Cc'] = C('c', gnd, c=C_)
    cir['c1'] = C('c', 'b', c=c1)
    cir['pin'] = CCVS('b', gnd, 's', gnd, r=Rs)
    cir['sense'] = VCVS('s', 'p', 'd', gnd, g=1.0)
    T_ideal = 2 * tau * np.log(3.0)
    return cir, T_ideal


def test_an_unimprovable_step_is_counted_and_named_instead_of_committed_in_silence():
    """⚠⚠ PEER REPORT, 2026-09-16, REPRODUCED -- and it is a DIAGNOSIS defect,
    not a Jacobian one.

    The report was "the free-period Jacobian is not the derivative of its
    residual".  It is not -- `dF/dx_in` is SINGULAR (shooting.py's own
    2026-09-02 measurement: sigma_min exactly 0, rank 1/3, 2/4, 1/3), so an
    exact Newton for that formulation does not exist and `I - dx_end/dx_0` is
    an approximation to a DIFFERENT, well-posed derivative.  The reporter
    withdrew that headline.  What survives is narrower and real: on a weakly
    limited tank (second multiplier 0.99) the step is UPHILL against the true
    derivative, the line search's four halvings cannot improve it, `fsolve`
    commits it anyway -- there is no other candidate -- and NOTHING SAID SO.
    Measured on the reporter's fixture, rebuilt here: the residual RISES with
    budget, 6.765e-07 at 25 iterations to 9.011e-07 at 200.

    Pinned here: the count reaches `infodict`, a failing solve NAMES it, a
    converging solve is untouched, and the message no longer claims that
    iterations never help -- which was measured for GEAR's solved-history
    stall and written as though it held for every multistep solve (a
    trapezoidal solve at multiplier 0.9 converges with a bigger budget on the
    reporter's circuit and at the default 25 on the van der Pol here).

    ⚠ NOT FIXED BY, measured on that fixture: `x0_unknown=True` (Jacobian
    becomes the derivative, 6e-07 against 1.000, step becomes descent -- and
    `||F||` still plateaus at 3.069e-05, identical at 25 and 200 iterations)
    or `tstab` (50 periods: still non-converged).  `radau` DOES converge there
    in 25 iterations, which is what the message steers to.
    """
    from pycircuit.circuit import analysis as _an
    circuit.default_toolkit = circuit.numeric

    ## the counter exists and is silent on a healthy solve
    def arctan(x):
        return (np.array([float(np.arctan(x[0]))]),
                np.array([[1.0 / (1.0 + x[0] ** 2)]]))
    _x, info, ier, _m = _an.fsolve(arctan, np.array([0.2]), maxiter=40,
                                   toolkit=circuit.numeric, full_output=True,
                                   line_search=True)
    assert ier == 1 and info['ls_unimproved'] == 0, (ier, info)

    ## the reporter's lambda_2 = 0.99 tank, rebuilt: series-resonant LC with a
    ## cubic negative resistance; their switch sits at a fixed control voltage,
    ## so it is the gl = 1e-4 S tank loss and is written as a resistor here.
    T_TANK = 6.28318530718e-07

    def tank():
        c = SubCircuit()
        c.add_node('p')
        c['L0'] = L('p', gnd, L=1e-4)
        c['C0'] = C('p', gnd, c=1e-10)
        c['Rl'] = R('p', gnd, r=1e4)
        c['N0'] = BSource('p', gnd, gnd, 'p',
                          i_func=lambda u: 1.01599560631e-04 * u
                          - 2.13274750776e-06 * u ** 3)
        return c

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        cir = tank()
        x0 = np.zeros(cir.n - 1)
        x0[0] = 1.0                      # the kick is ON the orbit (A = 1)
        pss = PSS(cir, method='trap', reltol=1e-10)
        pss.solve(period=T_TANK, timestep=T_TANK / 800, x0=x0,
                  maxiterations=25)
    assert not pss.converged, \
        'the weakly damped tank now converges under trap, so this fixture no ' \
        'longer exercises the uphill step it exists for'
    named = [str(w.message) for w in caught
             if 'line search could not improve' in str(w.message)]
    assert named, \
        'a solve whose step the line search could not improve said nothing: ' \
        '%r' % ([str(w.message)[:80] for w in caught],)
    msg = named[0]
    ## the count is real (24 of 25 measured), and the claim is narrowed
    import re
    n_named = int(re.search(r'could not improve (\d+) of', msg).group(1))
    assert n_named >= 1, msg
    assert 'UPHILL' in msg and 'radau' in msg, msg
    assert 'does not help' not in msg, \
        'the blanket "iterations do not help" claim is back; it was measured ' \
        'for gear\'s solved-history stall only'

    ## and a converging autonomous solve is untouched by the bookkeeping
    def vdp():
        c = SubCircuit()
        c.add_node('v')
        c['C'] = C('v', gnd, c=1.0)
        c['L'] = L('v', gnd, L=1.0)
        c['B'] = BSource('v', gnd, gnd, 'v',
                         i_func=lambda u: u - u ** 3 / 3.0)
        return c
    with quiet(AccuracyWarning):
        p2 = PSS(vdp(), method='trap', reltol=1e-9)
        p2.solve(period=6.6634, timestep=6.6634 / 200,
                 x0=np.array([2.0, 0.0]), maxiterations=40)
    assert p2.converged and abs(p2.period / 6.663571642 - 1.0) < 1e-6, \
        (p2.converged, p2.period)


def test_the_line_search_is_the_last_resort_and_reaches_the_shooting_path():
    """Owner decision 2026-09-08 ("Do 2"): the inner step Newton on the PSS
    path had no damping, and a nonlinearity fed by a branch current through
    a capacitor cannot be integrated undamped at any grid.  Pinned here:

    * the current-sense relaxation oscillator's INNER steps converge under
      radau and trap -- the solve no longer fails with the stage Newton's
      `NoConvergenceError`; what it reaches instead is the OUTER free-period
      Jacobian going singular (a separate, recorded obstacle), so the
      assertion is on the exception TYPE;
    * the search is the LAST resort: a tanh charge circuit at h = 100 tau,
      where the undamped stage Newton fails and the rescue ladder recovers
      the physical root, still gives v in [0, 1] through the plain
      transient -- a first-retry line search beat the ladder to v = -0.010.
    """
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    from pycircuit.circuit.nrsolver import NoConvergenceError
    circuit.default_toolkit = circuit.numeric
    with quiet(AccuracyWarning, ConvergenceWarning):
        for method in ('radau', 'trap'):
            cir, T = _current_sense_relaxation_oscillator()
            p = PSS(cir, method=method, reltol=1e-12)
            x0 = np.zeros(cir.n)
            x0[cir.get_node_index('vo')] = 1.0
            x0[cir.get_node_index('p')] = 0.5
            x0[cir.get_node_index('d')] = 0.5
            x0 = np.delete(x0, cir.get_node_index(gnd))
            try:
                p.solve(period=T, timestep=T / 400, x0=x0, maxiterations=60)
                reached = 'converged'
            except NoConvergenceError as e:
                reached = 'stage Newton: %s' % str(e)[:60]
            except np.linalg.LinAlgError:
                reached = 'outer Jacobian'
            ## `startswith`, not `!=`: the first version compared against the
            ## bare label and PASSED ON THE OLD CODE (the message carries a
            ## suffix) -- a regression test that cannot fail before the fix.
            assert not reached.startswith('stage Newton'), \
                '%s: the inner step Newton still fails -- %s' % (method, reached)
        ## the last-resort ordering, on the circuit that exposed a first-retry search
        cir = SubCircuit()
        cir.add_node('s'); cir.add_node('v')
        cir['vs'] = VPulse('s', gnd, v1=0.0, v2=1.0, td=0.0, tr=0.0, tf=0.0, pw=0.5, per=1.0)
        cir['g'] = BSource('s', 'v', gnd, 'v', i_func=lambda u: 4.0 * np.tanh(u / 0.01))
        cir['C'] = C('v', gnd, c=1.0)
        tr = Transient(cir, integrator=RadauIIA3Integrator())
        res = tr.solve(tend=2.0, timestep=0.25, x0=np.zeros(cir.n), fixed_timestep=True)
        v = np.asarray(res.x, dtype=float)[cir.get_node_index('v')]
    assert v.min() >= -1e-9 and v.max() <= 1.0 + 1e-9, \
        'the search pre-empted the ladder: v in [%.4f, %.4f]' % (v.min(), v.max())


@pytest.mark.fast_tier
@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.AccuracyWarning')
def test_warping_estimate_refuses_a_reading_the_interpolant_sets():
    """2026-09-08 (item 3): Part I's "only if" as a self-check.  The
    defect-correction estimate is the METHOD's error only while the
    interpolant's own defect is asymptotically smaller, and then it does
    not depend on the interpolant: the same pass through every second
    sample of the same solution must read the same slope.  Measured (radau
    and trap, 50 and 100 points): van der Pol half/full 1.0000; the
    relaxation orbit with a comparator edge (the case that read 0.09 and
    0.65 of the truth at 200 and 400 points): 0.0056 at 200 and 0.72 at
    400 under radau, 0.41 at 200 under trap -- all refused at the 5 %
    tolerance, the smooth case accepted with four orders of margin.
    Pinned on the cheap end of each: van der Pol radau at 50 points
    (trusted, ratio within 1e-3 of 1, the estimate unchanged by the check)
    and the relaxation orbit under trap at 200 points (trusted False, a
    warning naming the interpolant, the number still returned).
    """
    import numpy as np
    from pycircuit.circuit import circuit
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import R, C, VS, VCVS, BSource
    from pycircuit.circuit.shooting import PSS
    circuit.default_toolkit = circuit.numeric
    cir, _mu = _a10_vdp()
    p = PSS(cir, method='radau', reltol=1e-12)
    p.solve(period=2 * np.pi, timestep=2 * np.pi / 50, x0=np.array([2.0, 0.0]), maxiterations=100)
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        est = p.warping_estimate(periods=10)
    assert est['trusted'] is True and abs(est['check_ratio'] - 1.0) < 1e-3, est['check_ratio']
    plain = p.warping_estimate(periods=10, check=False)
    assert plain['trusted'] is None and plain['period_error'] == est['period_error']
    ## the relaxation orbit: a tanh comparator around an RC ramp, edge a few points wide
    R_ = 1e3; C_ = 1e-6; beta = 0.5; Vsat = 1.0; k = 20.0; c1 = 0.01 * C_
    tau = R_ * (C_ + c1); T = 2 * tau * np.log((1 + beta) / (1 - beta)); Ro = 1.0
    cir = SubCircuit()
    for n in ('vo', 'p', 'c', 'd', 'b'):
        cir.add_node(n)
    cir['Ro'] = R('vo', gnd, r=Ro); cir['Co'] = C('vo', gnd, c=(tau / 20) / Ro)
    cir['cmp'] = BSource('d', gnd, gnd, 'vo', i_func=lambda u: (Vsat / Ro) * np.tanh(k * u))
    cir['R1'] = R('vo', 'p', r=1e3); cir['R2'] = R('p', gnd, r=1e3)
    cir['Rc'] = R('vo', 'c', r=R_); cir['Cc'] = C('c', gnd, c=C_); cir['c1'] = C('c', 'b', c=c1)
    cir['pin'] = VS('b', gnd, v=0.0); cir['sense'] = VCVS('p', 'c', 'd', gnd, g=1.0)
    x0 = np.zeros(cir.n)
    x0[cir.get_node_index('vo')] = Vsat; x0[cir.get_node_index('p')] = beta * Vsat; x0[cir.get_node_index('d')] = beta * Vsat
    x0 = np.delete(x0, cir.get_node_index(gnd))
    p = PSS(cir, method='trap', reltol=1e-12)
    p.solve(period=T, timestep=T / 200, x0=x0, maxiterations=100)
    assert p.converged
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        est = p.warping_estimate(periods=10)
    assert est['trusted'] is False and abs(est['check_ratio'] - 1.0) > PSS.WARPING_CHECK_TOL, est['check_ratio']
    assert est['period_error'] is not None
    assert any('interpolant' in str(x.message) for x in w), [str(x.message) for x in w]


def test_no_analysis_answer_depends_on_a_stale_limiting_state():
    """⚠⚠ `cir.G(x)` IS NOT A PURE FUNCTION OF `x`.  `Diode` linearises around
    a stored `_vlim`, so any site evaluating `G` outside a converged solve
    reads whatever the last solve left behind.  The defence in this tree is
    PER-CALL-SITE (shooting routes `_G_at` through PCNR), which means a NEW
    site inherits the bug silently -- and one did, the same day the hazard was
    written down.

    This gate is the general form: poison the limiting state and assert every
    analysis answer is unmoved.  Measured before the fix:

        DC operating point      0.0e+00
        AC small signal         0.0e+00
        Transient final state   0.0e+00
        algebraic_conditioning  1.0e+00   <-- inherited it
    """
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.analysis_ss import AC
    poison = 5.0

    def poisoned():
        c = _diode_fixture()
        c['d'].__dict__['_vlim'] = poison
        return c

    ## ⚠⚠ THE INSTRUMENT ARM, AND IT IS NOT OPTIONAL.  A first version of this
    ## measurement poisoned via `limit()` with a FULL-CIRCUIT x vector; that
    ## takes the ELEMENT's local x, so it set `_vlim = 0.0` -- the unpoisoned
    ## value -- and every analysis read 0.0e+00 including the raw `G` that was
    ## known to move.  A clean sweep of nulls from a DEAD POISON.  Prove the
    ## poison moves something before believing it moves nothing.
    x = np.zeros(_diode_fixture().n)
    x[_diode_fixture().get_node_index('b')] = 0.7
    g_clean = np.asarray(_diode_fixture().G(x, defaultepar), dtype=float)
    g_poisoned = np.asarray(poisoned().G(x, defaultepar), dtype=float)
    assert np.max(np.abs(g_poisoned - g_clean)) > 1.0, (
        'the poison does not move a raw G evaluation, so every null below '
        'would be measuring nothing')

    def unmoved(name, fn):
        a = np.atleast_1d(np.asarray(fn(_diode_fixture()), dtype=complex)).ravel()
        b = np.atleast_1d(np.asarray(fn(poisoned()), dtype=complex)).ravel()
        scale = max(float(np.max(np.abs(a))), 1e-300)
        moved = float(np.max(np.abs(a - b))) / scale
        assert moved < 1e-12, (name, moved)

    unmoved('DC', lambda c: DC(c, refnode=gnd).solve().x)
    unmoved('AC', lambda c: np.asarray(
        AC(c).solve(np.array([1e3]), refnode=gnd).v('b', 'gnd')))
    unmoved('transient', lambda c: Transient(c).solve(
        refnode=gnd, tend=1e-3, timestep=1e-3 / 50, fixed_timestep=True).x[-1])
    unmoved('algebraic_conditioning', lambda c: algebraic_conditioning(c)[0])

    ## ⚠ AND AT AN OPERATING POINT HANDED IN (`dcx`, 2026-09-28): the small-
    ## signal analyses read `G` at the diode's stored state, not at the point
    ## -- after a DC sweep, its last point.  Measured with this poison before
    ## the fix: AC 1.0 (the gain collapsed), noise 1.0.
    from pycircuit.circuit.analysis_ss import Noise
    x_op = np.asarray(DC(_diode_fixture(), refnode=gnd).solve().x, dtype=float)
    unmoved('AC at a given operating point', lambda c: np.asarray(
        AC(c, dcx=x_op).solve(np.array([1e3]), refnode=gnd).v('b', 'gnd')))
    unmoved('noise at a given operating point', lambda c: Noise(
        c, inputsrc='vs', outputnodes=('b', gnd), dcx=x_op).solve(1e3)['Svnout'])


def test_the_accuracy_estimate_reads_the_devices_at_its_spline():
    """`PSS.warping_estimate`'s defect source reads `C` and `i` at points of
    the periodic spline, inside a transient of its own, and a stateful
    limiter (`Diode`) there read the tangent at that transient's state
    (`_limiting.devices_at` now puts it at the spline point).  Measured
    against the state-free twin of the diode (`_state_free_diode`), radau on
    this fixture: the lag 2.3e-5 off before the fix (1e-4 .. 2.4e-3 on a
    half-wave rectifier; 2026-09-28)."""
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.tests.test_analysis_transient import _state_free_diode
    circuit.default_toolkit = circuit.numeric

    def lag(cls):
        c = _diode_fixture()
        c['d'] = cls('b', gnd)
        pss = PSS(c, method='radau', reltol=1e-9)
        with quiet(AccuracyWarning):
            pss.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=60)
            return np.asarray(pss.warping_estimate(periods=5)['lag'], dtype=float)
    a, b = lag(_state_free_diode()), lag(Diode)
    assert np.max(np.abs(a - b)) <= 1e-9 * np.max(np.abs(a)), \
        np.max(np.abs(a - b)) / np.max(np.abs(a))


def test_a_shooting_solve_sitting_on_its_answer_says_why_it_failed_its_step_test():
    """⚠⚠ "PLAIN PSS ENDS AT N ~ 32 ON A PLL" WAS AN ARTEFACT OF THE STEP TEST.

    A /N PLL (VCO at N f_ref, the loop held at f_c = 100 Hz, 32 points per VCO
    cycle).  The A5 spike read its cost as ~N^1.5 ending in non-convergence at
    N = 32, and a seeding experiment then read "a fragile Newton basin" off
    which seeds happened to pass -- BOTH WRONG, and both mine: I never looked
    at a residual history.  It says: Newton reaches the solution in TWO
    iterations (|F| 7e-13) and the iterate never moves again, while the STEP of
    the VCO's output node -- a unit sine sampled at its zero crossing
    (-1.5e-05), so `reltol |x|` vanishes and only `vabstol` is left -- stays
    RANDOM at 4e-12 .. 1e-09: the rounding of a 1024-step traversal carrying an
    unknown of 3.2e7 (the VCO's frequency node; one ulp is 3.7e-09).  At
    `vabstol = 1e-12` the solve burned its budget and reported failure, or
    passed when a draw landed under 1e-12.  Fixed seed, |lambda| = 0.999371:

        N      vabstol = 1e-12            vabstol = 1e-9    1e-6
        32     142 s, NOT CONVERGED        5.7 s             5.7 s
        64     297 s                      12.9 s
        128    (never reached)            16.3 s

    ⚠ AND THE FIRST REPAIR WAS WRONG.  `fsolve` declared success on the
    signature (residual test held three iterations, step not contracting); the
    suite refused it twice, on fixtures whose `I - M` is SINGULAR -- there the
    residual is met on a whole manifold, the step wanders along the null
    direction, and "not converged" is the designed answer.  Both are rounding
    amplified by conditioning and the iteration cannot tell them apart.  So the
    signature is COUNTED, never acted on: the iteration is bit for bit what it
    was, and a failed solve now says why, naming both causes.
    """
    import warnings as _w
    from pycircuit.circuit.elements_hdl import VcoHdl
    from pycircuit.circuit import analysis as _an
    circuit.default_toolkit = circuit.numeric
    fref, K, N = 1e6, 1e-2, 32
    T = 1.0 / fref

    def solve(maxiterations, **kw):
        c = SubCircuit()
        for nd in ('ref', 'ph', 'dv', 'pd', 'ctl'):
            c.add_node(nd)
        c['Vref'] = VSin('ref', gnd, va=1.0, freq=fref, phase=0.0)
        c['X1'] = VcoHdl('ctl', gnd, 'vco', gnd, 'ph', f0=N * fref,
                         kvco=N * 2e4, va=1.0, modulus=float(N))
        c['DV'] = _PllPhaseDiv('ph', gnd, 'dv', gnd, n=float(N))
        c['PD'] = _PllMultPd('dv', gnd, 'ref', gnd, 'pd', gnd, k=K)
        c['Rf'] = R('pd', 'ctl', r=1e3)
        c['Cf'] = C('ctl', gnd, c=1e-9)
        names = [str(n) for n in c.nodes if str(n) != 'gnd!']
        x0 = np.zeros(c.n - 1)
        x0[names.index('X1._state0')] = 0.25 * N
        x0[names.index('ph')] = 0.25 * N
        pss = PSS(c, method='gear', reltol=1e-9, **kw)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            pss.solve(period=T, timestep=T / (32 * N), x0=x0,
                      maxiterations=maxiterations)
        said = [r for r in rec if 'STOPPED MOVING' in str(r.message)]
        lam = float(np.max(np.abs(np.linalg.eigvals(np.asarray(pss._monodromy)))))
        return pss, said, lam, names

    expected = np.exp(-np.pi * 2e4 * K * T)
    ## a tolerance above the floor: converged at once, nothing to say
    pss, said, lam, names = solve(60, vabstol=1e-6)
    assert pss.converged and not said and not pss.step_floor
    assert abs(lam / expected - 1.0) < 1e-5, (lam, expected)
    ## below it: the SAME orbit, reported as not converged -- and now with why
    pss2, said2, lam2, _n = solve(10, vabstol=1e-12)
    assert not pss2.converged
    assert len(said2) == 1, said2
    sf = pss2.step_floor
    assert names[sf['index'] % len(names)] == 'vco', (sf, names)
    assert sf['since'] <= 4 and sf['ratio'] > 3.0 and sf['step'] < 1e-7, sf
    assert abs(lam2 / lam - 1.0) < 1e-8, (lam2, lam)     # it WAS the solution
    msg = str(said2[0].message)
    assert 'ARITHMETIC FLOOR' in msg and 'SINGULAR' in msg

    ## ⚠ COUNTED, NOT ACTED ON: with the signature present the iteration is
    ## what it would have been without the flag.
    calls = []

    def stuck(x):
        ## residual met at once; a step that never contracts (a noisy J^-1 F)
        calls.append(1)
        k = len(calls)
        return (np.array([1e-4 * (1.0 if k % 2 else -1.2)]),
                np.array([[1.0]]))
    kw = dict(maxiter=9, reltol=1e-12, abstol=1e-3, xtol=1e-12,
              toolkit=circuit.numeric, full_output=True)
    ref = _an.fsolve(stuck, np.array([0.0]), **kw)
    n_ref = len(calls)
    del calls[:]
    got = _an.fsolve(stuck, np.array([0.0]), floor_detect=True, **kw)
    assert got[2] == ref[2] == 2 and len(calls) == n_ref
    assert float(got[0][0]) == float(ref[0][0])
    assert got[1]['step_floor'] and got[1]['step_floor']['index'] == 0
    assert ref[1]['step_floor'] is None
    ## and a LINEARLY converging Newton never shows the signature: its step
    ## contracts every iteration (a frozen Jacobian of 2 for 1 halves the error)
    def lin(x):
        return np.array([x[0] - 1.0]), np.array([[2.0]])
    got = _an.fsolve(lin, np.array([0.0]), maxiter=200, reltol=1e-12,
                     abstol=1e-3, xtol=1e-12, toolkit=circuit.numeric,
                     full_output=True, floor_detect=True)
    assert got[2] == 1 and got[1]['step_floor'] is None


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.ModelWarning')
def test_the_free_period_and_matrix_free_solves_record_the_stalled_step_signature_too():
    """`57600a8` gave the DENSE driven shooting solves a diagnosis for "sitting on
    the answer and still failing the step test" and said the free-period and
    matrix-free solves were not covered.  Now they are, the same way: COUNTED,
    never acted on.  Unit-level, because no cheap fixture stalls those paths on
    demand: the free-period solve hands `fsolve` the flag, and the matrix-free
    Newton records the signature with an iteration that is bit for bit what it
    was (same iterate, same verdict, same number of `build` calls).
    """
    from pycircuit.circuit import analysis as _an
    circuit.default_toolkit = circuit.numeric

    ## (1) the free-period solve asks fsolve for the signature
    seen = {}
    real = _an.fsolve

    def spy(*a, **kw):
        seen.update(kw)
        return real(*a, **kw)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    pss = PSS(cir, method='gear', reltol=1e-9)

    def f(z):
        return (np.array([z[0] - 1.0, z[1] - 2.0]),
                np.array([[1.0, 0.0], [0.0, 1.0]]))
    _an.fsolve = spy
    try:
        z, info, ier, _m = pss._free_period_solve(
            f, np.array([0.0, 3.0]), np.array([1e-9, 1e-9]),
            np.array([1e-9, 1e-15]), 1e-9, 20, 3.0)
    finally:
        _an.fsolve = real
    assert seen.get('floor_detect') is True
    assert 'step_floor' in info

    ## (2) the matrix-free Newton: a stuck system shows the signature, a
    ## converging one does not, and the iteration is unchanged either way
    calls = []

    def stuck(z):
        calls.append(1)
        k = len(calls)
        F = np.array([1e-4 * (1.0 if k % 2 else -1.2), 0.0])
        ## a COPY: GMRES works in place, and a matvec that returns its own
        ## input aliases its work vectors (a real matvec never does)
        return F, (lambda v: np.array(v, dtype=float))
    ## reltol 1e-6, not 1e-12: GMRES's rtol is a factor below reltol, and at
    ## 1e-12 it hit its iteration limit before the Newton loop ran once.
    ## UNDAMPED (`line_search=False`, the default since 2026-09-25): the
    ## signature is recorded on the iteration as it is, and this system's
    ## alternating |F| would otherwise spend its `build`s on halvings
    z1, info1, ier1, _m = pss._matrix_free_newton(
        stuck, np.zeros(2), np.array([1e-3, 1e-3]), np.array([1e-12, 1e-12]),
        1e-6, 9, line_search=False)
    assert ier1 == 2 and len(calls) == 9
    assert info1['step_floor'] and info1['step_floor']['index'] == 0 \
        and info1['step_floor']['since'] <= 2, info1

    def fine(z):
        return np.array([z[0] - 1.0, z[1] + 1.0]), (lambda v: np.array(v, dtype=float))
    z2, info2, ier2, _m = pss._matrix_free_newton(
        fine, np.zeros(2), np.array([1e-12, 1e-12]), np.array([1e-12, 1e-12]),
        1e-6, 20)
    assert ier2 == 1 and info2['step_floor'] is None
    np.testing.assert_allclose(z2, [1.0, -1.0], rtol=0, atol=1e-12)


def test_the_period_column_carries_the_constant_source_vector_of_an_autonomous_circuit():
    """`dphi/dT` for the DIRK and coupled-Radau kinds is built from the stage
    derivative `-(i(Y) + u)`, and `u` is the CONSTANT source vector of an
    autonomous circuit -- not zero (2026-09-20).

    ⚠ It used to be `-i(Y)`, "no source term", and both FD checks in the
    build sat on a source-free van der Pol where that is the same thing.
    On a row a source pins `i = -u` at convergence, so the column read
    `-u/T` there: on this van der Pol with a decoupled 1 kV node the
    oscillator rows agreed with the finite difference to 4 digits while the
    supply node read -157.9 = -1e3/T and its branch current +i_R/T.  On the
    tree's own phase fixture (a 1 kV DC supply) radau's column was
    `[-1e6, ~0, -1.8e-4, ...]` against the FD's `[0, -3.18, 6.28e3, ...]`,
    the first Newton step threw `x0` to 2.8e6 and the stage matrix went
    singular there -- the DEFAULT method failing the free-period solve on a
    shipped autonomous fixture, mislabelled "seed below the fundamental".
    Fixed: max |Mt - FD| / max |FD| = 4.6e-11 (radau) / 7.9e-11 (trbdf2) on
    the phase fixture, 1.2e-10 here; radau then finds the phase fixture's
    period at order 5: +1.52e-6 / +2.36e-8 / +6.5e-10 ppm at N = 100 / 200 /
    400.
    """
    circuit.default_toolkit = circuit.numeric
    mu = 1.0

    def vdp_vs():
        cir = SubCircuit()
        cir.add_node('v')
        cir.add_node('p')
        cir['C'] = C('v', gnd, c=4.0)
        cir['L'] = L('v', gnd, L=0.25)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: mu * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['vp'] = VS('p', gnd, v=1e3)
        cir['Rp'] = R('p', gnd, r=1e6)
        return cir
    T0 = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    x0 = np.zeros(4)
    x0[0], x0[1] = 2.0, 1e3
    g = PSS(vdp_vs(), method='gear', reltol=1e-8)
    with quiet(AccuracyWarning):
        g.solve(period=T0, timestep=T0 / 200, maxiterations=40,
                break_events=False, x0=x0)
    assert g.converged
    xs = np.asarray(g._period_state[1], float).ravel()
    T = float(g.period)

    for method in ('radau', 'trbdf2'):
        def trav(TT, want):
            q = PSS(vdp_vs(), method=method, reltol=1e-10)
            q._grid_fracs = None
            tms, hs = q._period_grid(TT, 200, None)
            return q._walk('stage', xs, tms, hs, T=TT, want_dT=want)
        Mt = np.asarray(trav(T, True).Pt, float).ravel()
        xe = lambda TT: np.asarray(trav(TT, False).x_end, float).ravel()
        fd = (xe(T * (1 + 1e-6)) - xe(T * (1 - 1e-6))) / (2e-6 * T)
        rel = np.max(np.abs(Mt - fd)) / np.max(np.abs(fd))
        assert rel < 1e-7, (method, rel, Mt, fd)
        ## the source rows are exactly the ones that used to be -u/T
        assert abs(Mt[1]) < 1e-6 and abs(Mt[3]) < 1e-9, (method, Mt)

    ## and the default method solves the phase fixture it used to fail on
    prev = None
    for n in (100, 200):
        pss = PSS(_phase_circuit(), method='radau', reltol=1e-10)
        with quiet():
            pss.solve(period=1e-3, timestep=1e-3 / n, maxiterations=60,
                      break_events=False)
        assert pss.converged and pss.autonomous
        e = abs(1e6 * (pss.period - 1e-3) / 1e-3)
        assert e < 1e-4, (n, e)
        if prev is not None:
            assert prev / e > 16.0, (prev, e)
        prev = e


def test_the_closing_period_column_is_the_default_again_and_its_polished_answer_does_not_depend_on_the_seed():
    """The 'closing' period column (B7c completed on every kind, 2026-09-21)
    and what its evidence was.  On a RAW single window of an adaptive run
    (`lte_grid(fold=False)`), whose seam sits on the slow branch, a free-period
    solve seeded 16 % low fails its per-step Newton under 'proportional'
    (every step rescaled with the trial period: the fine regions slide off
    the edges) and converges under 'closing' (inner steps frozen, only the
    last follows T) -- the closing step then holds the whole 3 s correction,
    and a proportional polish on the caller's fractions at the solved period
    lands on the reference-seeded answer.

    ⚠ 'auto' WAS REVERSED TO 'proportional' THE SAME DAY, AND RESTORED
    LATER THAT DAY ON THE DIAGNOSIS.  The reversal's evidence ("+690 ppm
    off on the fold, fails from 2 %") was the polish firing only beyond
    the zero-stability bound: a stretch of the last step INSIDE the bound
    is a change of discretisation and stayed in the answer (gear on its
    fold from seeds 0 / 0.5 / 2 % low: +1615 / +2470 / +3646 ppm at
    stretches 1.04 / 1.19 / 1.61, proportional +1431 at every seed; this
    raw window at 5 %: +3129 ppm at 2.28x, and its -33.5 ppm at 2 % was a
    cancellation), and the 2 % failure was the old mid-edge seam.  The
    polish is unconditional now, so the closing answer is proportional's
    at the solved period, seed-independent (+1408.5 / +1408.6 / +1408.9
    on the fold; -1319.1 / -1319.0 / -1318.7 / -1318.1 here at 0 .. 5 %).
    Pinned: 'auto' is bit-identical to 'closing' on the caller's grid;
    closing from 0, 2 and 5 % low agrees with itself to 5e-6 and with
    proportional to 1e-4;
    proportional fails from 16 % low where closing takes it.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    MU = 10.0
    T_REF = 19.098600502

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: MU * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    cir = vdp()
    p = PSS(cir, method='gear')
    xfull = np.zeros(cir.n)
    xfull[[str(n_) for n_ in cir.nodes].index('v')] = 2.0
    with quiet():
        fr, seed = p.lte_grid(19.1, x0=xfull, reltol=1e-5, fold=False)   # the raw window
    fr = np.asarray(fr, float)
    N = len(fr)
    T_LOW = 0.84 * p.lte_period

    def run(pc, Tseed):
        c_ = vdp()
        q = PSS(c_, method='gear', reltol=1e-9, period_column=pc)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=Tseed, timestep=Tseed / N, x0=seed, maxiterations=80,
                    break_events=False, grid=fr)
        if not q.converged:
            raise RuntimeError('did not converge')       # the expected failure, catchable
        return float(q.period), any('closing step ended' in str(w_.message) for w_ in rec)
    ## 'auto' is 'closing', bit for bit; the polished answer does not move
    ## with the seed, and it is proportional's
    Ta, _ = run('auto', p.lte_period)
    Tc0, _ = run('closing', p.lte_period)
    Tp, _ = run('proportional', p.lte_period)
    assert Ta == Tc0, (Ta, Tc0)
    Tc2, _ = run('closing', 0.98 * p.lte_period)
    Tc5, _ = run('closing', 0.95 * p.lte_period)
    assert abs(Tc2 / Tc0 - 1.0) < 5e-6 and abs(Tc5 / Tc0 - 1.0) < 5e-6, (Tc0, Tc2, Tc5)   # 0.8 / 2.1 ppm measured
    assert abs(Tc0 / Tp - 1.0) < 1e-4, (Tc0, Tp)
    ## proportional fails from 16 % low on the raw window ...
    try:
        run('proportional', T_LOW)
        raise AssertionError('proportional converged from the low seed on the raw window')
    except AssertionError:
        raise
    except Exception:
        pass
    ## ... and closing takes it there, with the polish, to the same answer
    Tc, second = run('closing', T_LOW)
    assert second
    assert abs(Tc / Tp - 1.0) < 1e-4, (Tc, Tp)


def test_lte_grid_measures_its_own_period_and_the_solve_seeds_from_it():
    """`lte_grid` cuts its window at the period the adaptive run SHOWS, not
    at the hint, and leaves it in `pss.lte_period` (2026-09-21, Andreas:
    "Do 1").  The grid is fractions of the period, so a hint 16 % off (my
    relaxation estimate (3 - 2 ln 2) mu at mu = 10 against the true 19.10)
    put the fine regions 16 % away from the edges: every method's per-step
    Newton failed, and the free-period solve needed the closing convention
    plus a second pass just to recover -- landing at gear -1464 ppm.  The
    run already contains the period: rising crossings of the fastest
    state, interpolated, read 19.1003 (+9e-5, the transient's own
    discretisation), consistent to 5e-5.  Cut there and seeded from it::

        method   period_column   T err (ppm)   c rel     second pass
        gear     proportional    -17           -4.5 %    no
        gear     auto            -22           -4.6 %    no
        trbdf2   proportional    +10           +0.36 %   no
        trbdf2   auto            +10           +0.36 %   no

    against -1436 / -195 ppm with the window cut at a right period but the
    grid built from the hint's length.  A hint more than 1 % off is warned
    with the observed period; an inconsistent detection keeps the hint and
    says so.
    """
    import warnings as _w
    circuit.default_toolkit = circuit.numeric
    MU = 10.0
    T_REF = 19.098600502

    def vdp():
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=1.0)
        cir['L'] = L('v', gnd, L=1.0)
        cir['B'] = BSource('v', gnd, gnd, 'v',
                           i_func=lambda u: MU * (u - u ** 3 / 3.0) + 0.3 * u * u)
        cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
        return cir
    T_HINT = (3.0 - 2.0 * np.log(2.0)) * MU                # 16 % low
    cir = vdp()
    p = PSS(cir, method='gear')
    xfull = np.zeros(cir.n)
    xfull[[str(n_) for n_ in cir.nodes].index('v')] = 2.0
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter('always')
        fr, seed = p.lte_grid(T_HINT, x0=xfull, reltol=1e-5)
    assert p.lte_period is not None
    assert abs(p.lte_period / T_REF - 1.0) < 3e-4, p.lte_period
    assert any('recurs every' in str(w_.message) for w_ in rec), 'the off hint must be named'
    fr = np.asarray(fr, float)
    N = len(fr)
    for method, bound in (('gear', 500.0), ('trbdf2', 100.0)):
        c_ = vdp()
        q = PSS(c_, method=method, reltol=1e-9)
        with _w.catch_warnings(record=True) as rec:
            _w.simplefilter('always')
            q.solve(period=p.lte_period, timestep=p.lte_period / N, x0=seed,
                    maxiterations=80, break_events=False, grid=fr)
        assert q.converged
        assert not any('closing step ended' in str(w_.message) for w_ in rec), method
        e = 1e6 * (float(q.period) - T_REF) / T_REF
        assert abs(e) < bound, (method, e)     # folded grid: gear O(300), trbdf2 +7..+10


def test_a_switching_window_is_re_cut_whatever_the_base_step():
    """`_land_fractions` cuts the window between a switch's two landed edges
    into `EVENT_WINDOW_STEPS` equal steps -- and until 2026-09-25 it did so
    only when NO base node fell inside (a window narrower than one base
    step).  A window spanning 2-7 base steps kept them, so the switch's
    whole transition was resolved by a handful of steps that depended on
    where the edges fell: on the PWM loop (window 0.0047 T) the staged
    radau error against radau at 3200 points read 6.5e-5 / 1.0e-5 / 7.0e-5
    / 4.2e-6 / 7.5e-5 / 3.9e-6 of the swing at N = 780..805, non-monotone
    for every method, and a staged 800-point radau "reference" was itself
    7.5e-5 off.  Re-cut, 5e-6 to 1.1e-5 at the same N, and radau / glm3 /
    glm2 / trbdf2 converge monotonically (a floor at the window's own 8
    steps).  Pinned on the fractions: a window over three base steps comes
    back as `EVENT_WINDOW_STEPS` equal steps; one narrower than a base step
    keeps the old rule's nodes; one wider than `EVENT_WINDOW_STEPS` base
    steps is left to the base grid.

    The count is the solve's `event_window_steps` Parameter (Andreas,
    2026-09-25: "Cut by 16. can this be controlled?"), 16 by default -- it
    was a class constant (8) that one instance could not change.  Pinned:
    the default, an explicit count reaching the cut, a run's own count
    reaching its stage, and a count below 2 refused."""
    circuit.default_toolkit = circuit.numeric
    K = PSS.EVENT_WINDOW_STEPS
    assert K == 16
    base = np.full(100, 0.01)

    def window_steps(a, b, ws=None):
        fr, landed = PSS._land_fractions(base, [a, b], window_steps=ws)
        pts = np.concatenate(([0.0], np.cumsum(fr)))
        i, j = int(np.argmin(np.abs(pts - a))), int(np.argmin(np.abs(pts - b)))
        return np.diff(pts[i:j + 1])
    ## over three base steps (0.4033 .. 0.4331): re-cut
    w = window_steps(0.4033, 0.4331)
    assert len(w) == K and np.allclose(w, (0.4331 - 0.4033) / K, rtol=1e-12), w
    ## narrower than a base step: the old rule's nodes, unchanged
    w = window_steps(0.5021, 0.5063)
    assert len(w) == K and np.allclose(w, (0.5063 - 0.5021) / K, rtol=1e-12), w
    ## wider than K base steps: the base grid's own nodes
    w = window_steps(0.2033, 0.2033 + 0.01 * (K + 2))
    assert len(w) > K and np.max(w) > 0.9 * 0.01, w
    ## an explicit count
    assert len(window_steps(0.4033, 0.4331, ws=8)) == 8

    ## a run's own count reaches its stage: the PWM loop's two windows are
    ## each cut into `event_window_steps`
    T = 1e-5
    pts = {}
    for ws in (8, 16):
        cir = _pwm_loop(T)
        p = PSS(cir, method='radau', reltol=1e-8, event_window_steps=ws)
        assert p.par.event_window_steps == ws
        with quiet():
            p.solve(period=T, timestep=T / 100, x0=np.zeros(cir.n - 1),
                    maxiterations=100)
        assert p.converged and len(p._state_event_fracs) == 4
        pts[ws] = len(p.waveform[0])
    assert pts[16] - pts[8] == 2 * (16 - 8), pts
    with pytest.raises(ValueError, match='event_window_steps'):
        PSS(_pwm_loop(T), method='radau', event_window_steps=1).solve(
            period=T, timestep=T / 20, x0=np.zeros(_pwm_loop(T).n - 1))

def test_the_frozen_phase_pin_is_the_raw_rule_kept_on_measurement():
    """The autonomous solve pins the coordinate moving fastest over the seed's
    first step, as a RAW `argmax |x_2 - x_1|` -- volts against amps.  A peer
    session showed the cost and it reproduces: on a van der Pol in LC form
    with real units (1 nF, 1 uH) seeded at v's turning point the rule pins v,
    nearly along the flow, and the bordered Jacobian (scaled by swing) has
    cond 120-160 against 3.75-9.3 for a pin chosen relative to each
    coordinate's swing.  That swing-scaled rule was built and REVERTED on
    measurement (2026-09-23, Andreas): on the comparator relaxation oscillator
    seeded at ten phases of its exact orbit it converged from 4/10 against the
    raw rule's 7/10.  Pinned: the raw rule's choice on the LC (v), that every
    method still converges there, and `phase_k` exposed for diagnosis."""
    circuit.default_toolkit = circuit.numeric
    Cv, Lv, mu = 1e-9, 1e-6, 0.3
    Z = np.sqrt(Lv / Cv)
    T0 = 2.0 * np.pi * np.sqrt(Lv * Cv)
    for method in ('gear', 'radau', 'trbdf2'):
        cir = SubCircuit()
        cir.add_node('v')
        cir['C'] = C('v', gnd, c=Cv)
        cir['L'] = L('v', gnd, L=Lv)
        cir['B'] = BSource('v', gnd, gnd, 'v', i_func=lambda u: (mu / Z) * (u - u ** 3 / 3.0))
        p = PSS(cir, method=method, reltol=1e-9)
        with quiet(AccuracyWarning):
            p.solve(period=T0, timestep=T0 / 200, x0=np.array([2.0, 0.0]), maxiterations=60)
        assert p.converged, method
        assert p.phase_k == 0, (method, p.phase_k)            # v: the raw rule's choice


def test_a_pss_shoots_at_the_temperature_it_was_given(monkeypatch):
    """The PSS's `epar` reaches every analysis it builds (the transients
    of `_new_transient`, the monodromy twin, the seed's operating point).
    ⚠ Measured before (2026-09-28):
    every analysis gets its own copy of `defaultepar`, so a `PSS(epar=...)`
    at 400 K shot at 300.15 K -- its waveform BIT-EQUAL to the default's,
    where forward transients at the two temperatures differ by 4.2e-4 V."""

    from pycircuit.circuit import dcanalysis
    from pycircuit.circuit.circuit import defaultepar
    from pycircuit.circuit.elements import Diode
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    seen = []
    _dc_solve = dcanalysis.DC.solve

    def recording_solve(self, *a, **k):
        seen.append(float(self.epar.T))
        return _dc_solve(self, *a, **k)
    monkeypatch.setattr(dcanalysis.DC, 'solve', recording_solve)

    def build():
        c = SubCircuit()
        c['vs'] = VSin(1, gnd, va=1.0, freq=1e3)
        c['R'] = R(1, 2, r=1e3)
        c['D'] = Diode(2, gnd)
        c['C'] = C(2, gnd, c=1e-7)
        return c
    hot = defaultepar.copy()
    hot.T = 400.0
    W = {}
    with quiet(AccuracyWarning):
        for label, kw in (('default', {}), ('hot', {'epar': hot})):
            p = PSS(build(), method='radau', reltol=1e-9, **kw)
            p.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
            assert p.converged
            W[label] = (np.asarray(p.waveform[1], float), p)
        c = build()
        res = Transient(c, reltol=1e-9, epar=hot).solve(
            tend=5e-3, timestep=1e-5, x0=np.zeros(c.n))
        ## the seed's operating point, which only an autonomous run solves
        del seen[:]
        PSS(_phase_circuit(), method='trap', reltol=1e-8, epar=hot).solve(
            period=1e-3, timestep=1e-3 / 100, maxiterations=30)
    p_hot = W['hot'][1]
    assert p_hot._transient().epar.T == 400.0
    assert seen and set(seen) == {400.0}, seen         # the seed's DC
    assert np.max(np.abs(W['hot'][0] - W['default'][0])) > 1e-5
    ## and it is the 400 K orbit: the forward run, settled after five
    ## periods (tau = 0.1 ms), ends where the hot PSS starts
    io = [str(n) for n in c.nodes].index('2')
    assert abs(float(np.asarray(res.x, float)[io, -1])
               - float(W['hot'][0][io, 0])) < 1e-5


def test_event_grid_does_not_read_an_earlier_runs_element_state():
    """`event_grid` walks `cir.next_event`, and an Idtmod's is a prediction
    from its last ACCEPTED point.  ⚠ Measured before (2026-09-28): after a
    forward `Transient` to 0.3 T on the same circuit object the grid gained
    7 spurious events (0.3529, 0.4529, ... -- that run's predicted wraps);
    after `lte_grid` it did not.  The element state is now reset first."""

    from pycircuit.circuit.elements import Idtmod, VPulse
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    T = 1e-6

    def build():
        c = SubCircuit()
        for nn in ('in', 'out'):
            c.add_node(nn)
        c['vs'] = VPulse('in', gnd, v1=0.0, v2=1.0, td=0.05 * T, tr=0.01 * T,
                         tf=0.01 * T, pw=0.4 * T, per=T)
        c['I'] = Idtmod('in', gnd, 'out', gnd, modulus=0.1 * T)
        c['Rl'] = R('out', gnd, r=1e3)
        return c
    with quiet():
        fresh = PSS(build(), method='radau')
        fresh.event_grid(T, npts=100)
        c = build()
        Transient(c).solve(tend=0.3 * T, timestep=T / 200, x0=np.zeros(c.n))
        after = PSS(c, method='radau')
        after.event_grid(T, npts=100)
    assert list(after.event_times) == list(fresh.event_times), \
        (after.event_times, fresh.event_times)
    assert len(fresh.event_times) == 4, fresh.event_times


def test_the_noise_and_the_ppv_are_read_at_the_pss_temperature(monkeypatch):
    """The noise densities (`CY`, `noise_amplitudes`) and the PPV's
    linearisation (`C`, `G`, `i`) are read at the PSS's `epar`, the
    temperature its orbit was solved at.  ⚠ Measured before (2026-09-28):
    every one used `defaultepar`, so a 400 K PSS's resistor noise was
    300.15 K's (the pnoise ratio 1.000000 where 400/300.15 is right) and its
    PPV linearised the devices at 300.15 K."""

    from pycircuit.circuit.circuit import defaultepar
    circuit.default_toolkit = circuit.numeric
    hot = defaultepar.copy()
    hot.T = 400.0

    def rc():
        c = SubCircuit()
        for nn in ('in', 'out'):
            c.add_node(nn)
        c['vs'] = VSin('in', gnd, va=0.1, freq=1e3)
        c['R'] = R('in', 'out', r=1e3)
        c['C'] = C('out', gnd, c=1e-7)
        return c
    S = {}
    with quiet():
        for label, kw in (('default', {}), ('hot', {'epar': hot})):
            cir = rc()
            pss = PSS(cir, method='radau', reltol=1e-9, **kw)
            pss.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=40)
            assert pss.converged
            k = cir.get_node_index('out')
            d = np.zeros(cir.n - 1)
            d[k - 1 if k > pss.irefnode else k] = 1.0
            S[label], _used = PAC(cir, toolkit=circuit.numeric).pnoise(
                pss, 300.0, d, maxsidebands=2)
    ## a linear circuit: the orbit is the same at both temperatures, so the
    ## ratio is 4kT's alone
    assert abs(S['hot'] / S['default'] - 400.0 / defaultepar.T) < 1e-9, \
        S['hot'] / S['default']

    ## the PPV's device reads, on an oscillator
    seen = []
    for name, at in (('C', 0), ('G', 0), ('i', 0), ('CY', 1)):
        orig = getattr(SubCircuit, name)

        def rec(self, x, *a, _o=orig, _at=at, **k):
            ep = k.get('epar', a[_at] if len(a) > _at else defaultepar)
            seen.append(float(ep.T))
            return _o(self, x, *a, **k)
        monkeypatch.setattr(SubCircuit, name, rec)
    with quiet():
        p = PSS(_phase_circuit(), method='radau', reltol=1e-8, epar=hot)
        p.solve(period=1e-3, timestep=1e-3 / 100, maxiterations=30)
        assert p.converged
        del seen[:]
        p.ppv()
    assert seen and set(seen) == {400.0}, sorted(set(seen))


def test_the_monodromy_twin_takes_every_setting_of_its_pss():
    """The twin a one-step LMM oscillator reads its monodromy from is the
    same analysis under another method: every Parameter but `method`.
    ⚠ Before (2026-09-28) it took four tolerances and `epar`, so `maxiter`,
    `pcnr`, the solvers, `relref` and the LTE floors never reached the map
    its PPV was read on."""
    circuit.default_toolkit = circuit.numeric
    with quiet(AccuracyWarning):
        p = PSS(_phase_circuit(), method='trap', reltol=1e-8, maxiter=77,
                relref='pointlocal', lte_vabstol=1e-9, TRTOL=5.0)
        p.solve(period=1e-3, timestep=1e-3 / 200, maxiterations=30)
        assert p.converged
        tw = p.monodromy_twin()
    assert tw is not p and tw.par.method == 'radau'
    for name in ('reltol', 'maxiter', 'relref', 'lte_vabstol', 'TRTOL',
                 'steadyratio', 'epar'):
        assert getattr(tw.par, name) == getattr(p.par, name), \
            (name, getattr(tw.par, name), getattr(p.par, name))


def test_a_pss_keeps_the_toolkit_it_was_given():
    """`PSS(cir, toolkit=...)` passes its toolkit on.  ⚠ Before (2026-09-29)
    it was accepted and dropped: every PSS ran on the default toolkit, and
    so did every twin built with `toolkit=self.toolkit` -- found by the
    dead-argument guard over the shooting package
    (`test_no_dead_knobs.py`)."""
    from pycircuit.circuit import toolkit as _toolkit
    circuit.default_toolkit = circuit.numeric
    tk = _toolkit.sparse_numeric
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, va=1.0, freq=1e3)
    c['R'] = R(1, gnd, r=1e3)
    assert PSS(c, toolkit=tk).toolkit is tk
    assert PSS(c).toolkit is circuit.default_toolkit


def _review_vdp(method='trap', npts=100, ref=None, ac=False):
    """van der Pol (Q = 8) solved autonomously, on ground or on its own node
    `ref` as the reference; `ac` adds a small-signal current into the tank."""
    circuit.default_toolkit = circuit.numeric
    mu = 1.0 / (2.0 * np.pi * 8.0)
    cir = SubCircuit()
    cir.add_node('v')
    cir['C'] = C('v', gnd, c=1.0)
    cir['L'] = L('v', gnd, L=1.0)
    cir['B'] = BSource('v', gnd, gnd, 'v',
                       i_func=lambda u: mu * (u - u ** 3 / 3.0))
    cir['n'] = IS('v', gnd, i=0.0, noisePSD=1e-6)
    if ac:
        cir['ac'] = IS('v', gnd, i=0.0, iac=1.0)
    T = 2.0 * np.pi / np.sqrt(1.0 - mu ** 2 / 4.0)
    rn = gnd if ref is None else cir.get_node(ref)
    pss = PSS(cir, method=method, reltol=1e-12, irefnode=rn)
    with quiet(AccuracyWarning):
        pss.solve(period=T, timestep=T / npts, refnode=rn, maxiterations=200,
                  x0=np.array([-2.0 if ref is not None else 2.0, 0.0]),
                  state_events=False, break_events=False)
    assert pss.converged
    return cir, pss


def test_every_re_solve_keeps_the_callers_event_settings(monkeypatch):
    """The review of 2026-09-30 (F3): a re-solve of this analysis -- the
    monodromy twin a trap/euler oscillator's PPV is read on, `grid_error`'s
    refined runs -- is the same analysis on another grid or integrator, so
    it takes the caller's `state_events` and `break_events`.  They took the
    DEFAULTS (both True) until then: the twin of an unstaged solve was
    staged, and `grid_error` compared a staged refinement with an unstaged
    base."""
    _cir, pss = _review_vdp()
    seen = []
    solve = PSS.solve

    def recording(self, *a, **kw):
        seen.append((kw.get('state_events', True), kw.get('break_events')))
        return solve(self, *a, **kw)

    monkeypatch.setattr(PSS, 'solve', recording)
    with quiet(AccuracyWarning):
        pss.ppv()
        pss.grid_error(lambda q: float(q.period), levels=2)
    assert len(seen) >= 2, seen
    assert all(se is False and be is False for se, be in seen), seen


def test_a_period_of_fewer_than_two_steps_is_refused():
    """The review of 2026-09-30 (F15): `timestep > period / 2` leaves fewer
    than two steps; one step did not converge and none raised an IndexError
    from inside the grid (measured, every method).  Refused up front."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['V'] = VSin('in', gnd, va=1.0, freq=1e5)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-9)
    for frac in (2.0, 1.0, 0.51):
        with pytest.raises(ValueError, match='at least two'):
            PSS(c, method='trap').solve(period=1e-5, timestep=frac * 1e-5)


def test_a_re_solve_after_a_settings_change_is_the_solve_those_settings_give():
    """The review of 2026-09-30 (F9): the inner Transient is cached for the
    analysis's life, and was built from the settings of the FIRST solve --
    `pss.par.method = 'radau'` on a solved gear PSS crashed the re-solve (a
    Gear-2 transient under a stage walk), trap -> gear raised the solved-
    history refusal, and a new `reltol` never reached the steps.  The cache
    is dropped when the settings it was built from moved (`_settings_key`):
    the re-solve is bit for bit a fresh analysis's."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['V'] = VSin('in', gnd, va=1.0, freq=1e5)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-9)

    def solved(p):
        with quiet(AccuracyWarning):
            p.solve(period=1e-5, timestep=1e-5 / 40)
        return np.asarray(p.waveform[1], float).copy()

    for attr, v0, v1 in (('method', 'gear', 'radau'), ('method', 'trap', 'gear'),
                         ('reltol', 1e-4, 1e-12)):
        p = PSS(c, **{attr: v0})
        solved(p)
        setattr(p.par, attr, v1)
        assert np.array_equal(solved(p), solved(PSS(c, **{attr: v1}))), (attr, v1)
    assert p._transient().par.reltol == 1e-12


def test_a_timestep_of_t_over_n_gives_n_steps():
    """F21 of the review of 2026-09-30 (Andreas: "Fix it so that one gets the
    actual requested steps"): `T / (T / N)` is `N - 1e-14` in floating point
    for 4-6 % of N, and `int(period / timestep)` made those N - 1 -- 154 of
    the suite's 1397 solves ran a step short (400 points were 399).  The
    count is floored with a relative slack (`steps_in`): `T / N` gives N,
    `T / (N + 0.5)` still gives N, a timestep that does not divide the
    period still gives the floor."""
    from pycircuit.circuit.shooting._numerics import steps_in
    for T in (1e-6, 1e-5, 1.391837294386364e-06, 2.0 * np.pi, 6.2835, 1e-3):
        for N in range(2, 5001):
            assert steps_in(T, T / N) == N, (T, N)
            assert steps_in(T, T / (N + 0.5)) == N, (T, N)
    assert steps_in(1.0, 0.3) == 3
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['V'] = VSin('in', gnd, va=1.0, freq=1e6)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-12)
    for N in (80, 100, 400):
        p = PSS(c, method='gear')
        with quiet(AccuracyWarning):
            p.solve(period=1e-6, timestep=1e-6 / N)
        assert len(p.factored_period().times) == N + 1, (N, len(p.factored_period().times))


@pytest.mark.parametrize('method', ['gear', 'radau', 'trap', 'glm3'])
def test_a_solved_pss_is_freed_when_its_last_reference_goes(method):
    """The review of 2026-09-30 (M1): a solved PSS was a REFERENCE CYCLE --
    `shooting_residual` is a closure over the analysis, and the factored
    period pointed back at it -- so dropping it freed nothing until a full
    garbage collection, and a sweep of 30 solves held all 30 (+7 MB each).
    The closures and the period hold the analysis weakly now: with the
    collector OFF, `del` frees a solved PSS whose factored period and PPV
    were built."""
    import gc
    import weakref
    _cir, pss = _review_vdp(method=method, npts=60)
    with quiet():
        pss.factored_period()
        pss.ppv()
    res = pss.shooting_residual
    z0 = np.asarray(pss._period_state[1], dtype=float)
    gc.collect()
    gc.disable()
    try:
        r = weakref.ref(pss)
        del pss
        assert r() is None, 'the solved PSS is still alive after del'
    finally:
        gc.enable()
    ## (the residual it kept says what happened, rather than a bare
    ## ReferenceError from deep inside)
    with pytest.raises(RuntimeError, match='keep the PSS alive'):
        res(np.append(z0, 1.0) if method != 'glm3' else z0)


def test_pss_shares_the_transient_s_parameters():
    """The tolerances, `maxiter`, `pcnr`, the LTE floors, `TRTOL`, `relref`
    and `analysis` -- the names forwarded to the inner transient -- are
    `Transient`'s own Parameter objects (`pcnr` a copy worded for the
    shooting), not re-declarations whose defaults matched by hand (the
    review's O18, 2026-10-01)."""
    from pycircuit.circuit.transient import Transient
    T = {p.name: p for p in Transient.parameters}
    P = {p.name: p for p in PSS.parameters}
    for name in ('analysis', 'reltol', 'iabstol', 'vabstol', 'maxiter',
                 'lte_vabstol', 'lte_iabstol', 'TRTOL', 'relref'):
        assert P[name] is T[name], name
    assert P['pcnr'] is not T['pcnr'] and 'inner transient' in P['pcnr'].desc
    assert (P['pcnr'].default, P['pcnr'].unit) == (T['pcnr'].default,
                                                   T['pcnr'].unit)
    ## (and radau's cost transform and the chord Jacobian, worded for the
    ## shooting)
    for name in ('radau_transform', 'chord_jacobian'):
        assert P[name] is not T[name] and P[name].default is T[name].default
    ## ... and the inner steps take them (`_new_transient`)
    p = PSS(_q20_rlc(), method='gear', chord_jacobian=True)
    assert p._transient().par.chord_jacobian is True


def test_the_converged_replay_is_the_factored_period():
    """The review's S6 (2026-10-01): a converged multistep solve's replay
    walks the period FACTORED and keeps it, so `factored_period()` does not
    walk the converged period a second time -- 15 % of a van der Pol gear
    solve-and-factor, 11 % of a compact MOSFET's.  The map is the one a
    second walk builds, bit for bit, on gear's pair and on a driven plain
    map opened either way; a trap oscillator's stays its twin's (none is
    kept at the solve), and past `REPLAY_FACTOR_BUDGET` none is kept."""
    circuit.default_toolkit = circuit.numeric

    def dense(fp):
        return np.column_stack([np.asarray(fp.matvec(e), float)
                                for e in np.eye(fp.width)])

    def vdp():
        cir = _vdp_with_noise(1e-6)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        return cir, {'period': 6.6634, 'timestep': 6.6634 / 100, 'x0': x0,
                     'maxiterations': 60}

    def q20(**kw):
        return _q20_rlc(), dict(period=1e-3, timestep=1e-3 / 60,
                                maxiterations=30, **kw)

    maps = {}
    for label, (cir, skw), method in (
            ('vdp gear', vdp(), 'gear'),
            ('driven trap', q20(), 'trap'),
            ('driven trap at x0', q20(x0_unknown=True), 'trap'),
            ('driven theta', q20(), 'theta')):
        p = PSS(cir, method=method, reltol=1e-10)
        with quiet(AccuracyWarning):
            p.solve(**skw)
        assert p.converged and p._factored_period_cache is not None, label
        steps = []
        orig = p.solve_timestep
        p.solve_timestep = (lambda *a, _s=steps, _o=orig, **k:
                            _s.append(1) or _o(*a, **k))
        fp = p.factored_period()
        assert not steps, f'{label}: factored_period walked again'
        M1 = dense(fp)
        p._factored_period_cache = None
        M2 = dense(p.factored_period())
        assert steps and np.array_equal(M1, M2), label
        maps[label] = M1

    ## a trap oscillator's period is its twin's: nothing kept at the solve
    cir, skw = vdp()
    p = PSS(cir, method='trap', reltol=1e-10)
    with quiet(AccuracyWarning):
        p.solve(**skw)
    assert p.converged and p._factored_period_cache is None
    assert p.factored_period().T != p.period      # the twin's own period

    ## past the budget the replay keeps nothing, and the map is unchanged
    cir, skw = vdp()
    p = PSS(cir, method='gear', reltol=1e-10)
    p.REPLAY_FACTOR_BUDGET = 0
    with quiet(AccuracyWarning):
        p.solve(**skw)
    assert p.converged and p._factored_period_cache is None
    assert np.array_equal(dense(p.factored_period()), maps['vdp gear'])


def test_a_stage_methods_converged_replay_is_its_factored_period():
    """The speed plan's P3 (2026-10-02): a stage method's converged replay
    IS the walk `factored_period()` takes -- a transient of its own on the
    solved grid's fractions -- with the waveform read off it and each step's
    stage states recorded, so the first small-signal call factors that
    record instead of walking the period again (radau is the default).
    Nothing is factored at the solve.  The map is the on-demand walk's bit
    for bit: radau on van der Pol, on a driven stateful diode, and TR-BDF2
    (diagonally implicit) on a pulse-driven non-uniform grid.  The waveform
    is a rounding apart from a replay on the PSS's own transient at the
    solve's `(t, h)`; a GLM records nothing."""
    from copy import copy

    from pycircuit.circuit.elements import VPulse
    circuit.default_toolkit = circuit.numeric

    def dense(fp):
        return np.column_stack([np.asarray(fp.matvec(e), float)
                                for e in np.eye(fp.width)])

    def vdp():
        cir = _vdp_with_noise(1e-6)
        x0 = np.zeros(cir.n - 1)
        x0[0] = 2.0
        return cir, {'period': 6.6634, 'timestep': 6.6634 / 100, 'x0': x0,
                     'maxiterations': 60}

    def diode():
        return _diode_fixture(), {'period': 1e-3, 'timestep': 1e-3 / 60,
                                  'maxiterations': 30}

    def pulsed():
        T = 1e-3
        c = SubCircuit()
        c['vp'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=0.0, tr=T / 50,
                         tf=T / 50, pw=T / 2 - T / 50, per=T)
        c['r'] = R('a', 'b', r=1e3)
        c['c'] = C('b', gnd, c=1e-7)
        c['d'] = Diode('b', gnd)
        return c, {'period': T, 'timestep': T / 60, 'maxiterations': 30}

    for label, fixture, method, uniform in (
            ('vdp radau', vdp, 'radau', True),
            ('diode radau', diode, 'radau', True),
            ('pulsed trbdf2', pulsed, 'trbdf2', False)):
        cir, skw = fixture()
        p = PSS(cir, method=method, reltol=1e-10)
        with quiet(AccuracyWarning):
            p.solve(**skw)
        assert p.converged, label
        assert p._factored_period_cache is None, label
        _, x0, _, times, hs, _, _ = p._period_state
        assert len(p._stage_replay[2]) == len(times) - 1, label
        _hs = np.asarray(hs, dtype=float)
        assert (float(np.ptp(_hs)) / float(_hs.max()) < 1e-9) == uniform, label
        steps = []
        orig = p.solve_timestep
        p.solve_timestep = (lambda *a, _s=steps, _o=orig, **k:
                            _s.append(1) or _o(*a, **k))
        M1 = dense(p.factored_period())
        assert not steps, f'{label}: factored_period walked again'
        assert p._stage_replay is None, label
        p._factored_period_cache = None
        M2 = dense(p.factored_period())
        assert steps and np.array_equal(M1, M2), label
        p.solve_timestep = orig

        ## the waveform: the solve's time axis, and a rounding from a replay
        ## on the PSS's own transient at the solve's (t, h)
        t1, X1 = p.waveform
        assert np.array_equal(t1, times), label
        p._begin_period(x0)
        X = [x0]
        for t, h in zip(times[1:], hs[:len(times) - 1]):
            X.append(copy(p.solve_timestep(X[-1], t, h)))
        X0 = np.array([np.asarray(p._insert_refnode(x), float) for x in X]).T
        assert X0.shape == X1.shape, label
        scale = np.max(np.abs(X0), axis=1)[:, None] + 1e-30
        assert np.max(np.abs(X1 - X0) / scale) < 1e-12, label

    ## a GLM's factored walk is not its replay: nothing recorded
    cir, skw = vdp()
    p = PSS(cir, method='glm2', reltol=1e-10)
    with quiet(AccuracyWarning):
        p.solve(**skw)
    assert p.converged and p._stage_replay is None


def test_the_factored_period_walks_the_solved_grid_itself():
    """A self-starting method's factored period walks the SOLVE'S OWN
    grid, not its step fractions re-summed (2026-10-02, found by the speed
    plan's P3 gate).  Rebuilt from fractions, a node can come out an ulp off
    the solved one, and where a landed edge with no ramp sits on it the step
    reads the source on the other side of the jump: on this pulsed RC the
    factored walk then did not close -- ``x_N - x_0`` 5.7e-11 at 100 points
    and 6.9e-9 at 200 (1.4e-20 at 400, where the node happened to land) --
    so the map linearised was not the one solved.  Its times are the
    waveform's exactly, and it closes to the solve's own residual."""
    from pycircuit.circuit.elements import VPulse
    circuit.default_toolkit = circuit.numeric
    T = 1e-6
    for N in (100, 200, 400):
        c = SubCircuit()
        c['vs'] = VPulse(1, gnd, v1=0.0, v2=1.0, td=0.0125 * T, tr=0.0,
                         tf=0.0, pw=0.4 * T, per=T)
        c['R'] = R(1, 2, r=1e3)
        c['C'] = C(2, gnd, c=3e-10)
        p = PSS(c, method='radau', reltol=1e-10)
        with quiet():
            p.solve(period=T, timestep=T / N, maxiterations=40)
        assert p.converged and p.break_events, N
        fp = p.factored_period()
        assert np.array_equal(np.asarray(fp.times, dtype=float),
                              np.asarray(p.waveform[0], dtype=float)), N
        x0 = np.asarray(p._period_state[1], dtype=float)
        gap = float(np.max(np.abs(np.asarray(fp.x_last, dtype=float) - x0)))
        assert gap < 1e-13 * float(np.max(np.abs(x0))), (N, gap)


def test_pss_forwards_radau_s_cost_transform():
    """The review's S15 (2026-10-01): a chord Jacobian for compact models was
    measured first -- `G` and `C` are 92 % of a compact MOSFET's PSS solve --
    and radau's already exists: `Transient`'s `radau_transform` (simplified
    Newton through eig(A^-1)), 107 -> 34 s on that solve.  PSS forwards it
    to every inner transient ('auto' by default since 2026-10-01: on where the
    compiled device Jacobians are expensive, so off here).  The inner steps take the transform; the orbit agrees with
    the full Newton's to its tolerance; a changed setting rebuilds the
    inner transient."""
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric
    calls = []
    orig = Transient._rk_step_transformed
    had = '_rk_step_transformed' in Transient.__dict__

    def counted(self, *a, **k):
        calls.append(1)
        return orig(self, *a, **k)

    def solve(p):
        x0 = np.zeros(p.cir.n - 1)
        x0[0] = 2.0
        with quiet(AccuracyWarning):
            p.solve(period=6.6634, timestep=6.6634 / 100, x0=x0,
                    maxiterations=60)
        assert p.converged
        return np.asarray(p.waveform[1], dtype=float)

    Transient._rk_step_transformed = counted
    try:
        p = PSS(_vdp_with_noise(1e-6), method='radau', reltol=1e-10)
        X0 = solve(p)
        ## (the default is 'auto', off on van der Pol's hand-written
        ## elements)
        assert not calls and p._transient().par.radau_transform == 'auto'
        X1 = solve(PSS(_vdp_with_noise(1e-6), method='radau', reltol=1e-10,
                       radau_transform=True))
        assert calls
        assert np.max(np.abs(X1 - X0)) / np.max(np.abs(X0)) < 1e-8
        ## a changed setting reaches the steps (`_settings_key`)
        del calls[:]
        p.par.radau_transform = True
        solve(p)
        assert calls and p._transient().par.radau_transform is True
    finally:
        if had:
            Transient._rk_step_transformed = orig
        else:
            del Transient._rk_step_transformed
