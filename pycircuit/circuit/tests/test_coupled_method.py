"""STAGE 12B -- `coupled_method`, how the coupled path corrects `h`.

    'approx'   Fang sec. 3.4 -- the step comes from the error RATIO (eq 17) and
               the solution is corrected by eq (18). The default, and since
               2026-09-27 the only one.
    'bordered' Fang eq (12)/(14) -- a linearised Newton step on the LTE
               equation.  RETIRED 2026-09-27: once its double-counted `q^T dv0`
               term was removed (2026-09-19) it took the same steps as 'approx'
               to every printed digit; the record below is of why it was built.

**Eq (12) is usable only because its denominator is computed analytically.** The
paper forms it as `q^T dxh + d`, and those two terms are how the solution moves
with the step size and how the extrapolation moves; both are approximately
`dv/dt`, so the difference is the derivative of the truncation error -- tiny by
construction and computed as a difference of two large numbers.

Measured at h = 3.48e-7 on a driven RC, against a ground truth obtained by
RE-SOLVING the circuit at perturbed `h` (`benchmarks/transient_review/`):

    ground truth   +4.678e6
    analytic       +4.392e6   ratio  0.939
    subtraction    -9.680e5   ratio -0.207   <- wrong SIGN

The subtracted terms are -1.310e8 and +1.300e8; the result is 0.74% of the
larger. That wrong sign is why an earlier attempt at eq (12) drove the step size
down four decades while the error sat far below its band.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.transient import Transient

TAU = 1e-4
W = 2 * np.pi * 1e3


def _rc():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def _analytic(t):
    A = 1.0 / np.sqrt(1.0 + (W * TAU) ** 2)
    phi = np.arctan(W * TAU)
    return A * (np.sin(W * t - phi) + np.sin(phi) * np.exp(-t / TAU))


def _run(method=None):
    tran = Transient(_rc(), toolkit=numeric, reltol=1e-5)
    if method is not None:
        tran.par.coupled_method = method
    with quiet():
        res = tran.solve(tend=5e-4, timestep=1e-5, coupled_lte=True)
    t = np.asarray(res.v('b').x, dtype=float).ravel()
    v = np.asarray(res.v('b').y, dtype=float).ravel()
    return tran.statistics, float(np.max(np.abs(v - _analytic(t))[2:]))


def test_the_default_is_the_method_with_the_measured_record():
    """`approx` is sec. 3.4, and it is what every stage-12 number was taken on."""
    tran = Transient(_rc(), toolkit=numeric)
    assert tran.par.coupled_method == 'approx'


def test_the_coupled_step_solves_the_circuit_and_takes_no_rejections():
    """Figure 3 has no rejection branch."""
    st, err = _run('approx')
    assert st.rejected_steps == 0, 'took %d rejections' % st.rejected_steps
    assert st.accepted_steps > 10
    assert err < 5e-3, 'max error %g' % err


def test_bordered_is_retired_and_says_why():
    """`coupled_method='bordered'` (Fang eq 12/14) raises, naming why: once its
    double-counted `q^T dv0` term was removed (2026-09-19) it took the same
    steps as 'approx' -- measured before retiring it, on this smooth RC and
    on the pulsed one: the same accepted steps (100 / 264), Newton
    iterations (222 / 1064), rejections and error, to every printed digit."""
    tran = Transient(_rc(), toolkit=numeric)
    tran.par.coupled_method = 'bordered'
    with pytest.raises(ValueError, match='retired'):
        tran.solve(tend=5e-5, timestep=1e-5, coupled_lte=True)


def test_an_unknown_method_is_refused():
    """A typo must not silently select the default."""
    tran = Transient(_rc(), toolkit=numeric)
    tran.par.coupled_method = 'schur'
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=5e-5, timestep=1e-5, coupled_lte=True)
    assert 'schur' in str(exc.value)


def _pulsed_rc():
    from pycircuit.circuit.elements import VPulse
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=1e-6, tf=1e-6,
                     pw=2e-5, per=5e-5)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def _pulse_run(method):
    tran = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5)
    tran.par.coupled_method = method
    with quiet():
        tran.solve(tend=6e-5, timestep=1e-6, coupled_lte=True)
    return tran.statistics


def test_the_coupled_step_rejects_few_steps_on_a_pulsed_circuit():
    """Figure 3's promise has to survive real discontinuities, not just smooth
    drives -- which is the whole reason a breakpoint circuit is in this file.

    RE-DERIVED at F13 (doc/transient_review_260820.md): the old `== 0`
    passed vacuously -- the coupled path never counted rejections at all.
    With the counter live, the honest statement of Figure 3 is narrower:
    the steps the method SOLVES for take no rejections (the smooth test
    above asserts exactly 0), while HELD steps -- breakpoint- or
    tend-truncated, whose size was never the method's to choose -- may
    retry when their imposed size fails the error test.  Measured at
    re-derivation: 16/985 and 14/1035 rejections, all at edges.  Bound
    them to a small fraction rather than pretending they are zero."""
    for method in ('approx',):
        st = _pulse_run(method)
        assert st.rejected_steps <= 0.05 * st.accepted_steps, \
            '%s: %d rejections against %d accepted -- held-step retries ' \
            'should be a few percent at worst' \
            % (method, st.rejected_steps, st.accepted_steps)
        assert st.breakpoints_hit > 0, '%s never hit an edge' % method

def test_coupled_tline_matches_standard_path():
    """The CPU coupled path runs delay lines now -- three fixes, each traced:

    - `TLine.dudt` written (derivative of the history interpolation), so the
      coupled residual's `p` vector carries the source term it used to refuse.
    - Kink discipline ported from the JAX fix: the step ring is emptied on a
      breakpoint landing, source corners are echoed as wavefront arrivals
      (corner + k*TD), and the solve's growth is capped at the breakpoint --
      without the cap the entry-h truncation test cleared a corner that the
      solved h then straddled (measured: entry 6.78e-11 under the 1.2e-9
      corner, solved 8.97e-11, landing 8.9e-12 past it).
    - The history interpolation is monotone-limited: a quadratic stencil
      spanning a recorded kink overshot the reflected EMF to 1.009 against
      samples bounded by 1.000, and a band-blind step accepted a solution
      against that phantom, which no later step could reconcile (an
      h-independent LTE floor of exactly the pollution).

    Gate: the pulsed matched line livelocked at t=2.01e-9 before the fixes
    (NoConvergenceError, h collapsed to 1e-16); it now lands on tend and
    matches the standard Gear2 path to 5.6e-16.  The mismatched RC load
    (its far-end reflection exercises the limiter) completes to the correct
    steady level: Gamma = 1/3, so vb settles at 2/3 of the 1 V swing.
    """
    from pycircuit.circuit.elements import R as _R, C as _C, VPulse, TLine
    from pycircuit.circuit.integrator import Gear2Integrator

    def line(rc_load):
        c = SubCircuit()
        c.add_node('a'); c.add_node('b')
        c['V1'] = VPulse('s', gnd, v1=0.0, v2=1.0, td=1e-9, tr=2e-10,
                         tf=2e-10, pw=1e-8, per=1e-7)
        c['Rs'] = _R('s', 'a', r=50.0)
        c['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
        c['Rl'] = _R('b', gnd, r=100.0 if rc_load else 50.0)
        if rc_load:
            c['Cl'] = _C('b', gnd, c=2e-12)
        return c

    ## Matched line: bit-close to the standard path.
    ref = Transient(line(False), toolkit=numeric, reltol=1e-4,
                    integrator=Gear2Integrator(), uic=True,
                    timestep_max=2e-10)
    with quiet():
        rr = ref.solve(gnd, tend=8e-9, timestep=2e-10)
    tr = np.asarray(rr.sweep_values, float)
    vr = np.asarray(rr.v('b'), float).reshape(-1)

    tran = Transient(line(False), toolkit=numeric, reltol=1e-4,
                     uic=True, timestep_max=2e-10)
    with quiet():
        res = tran.solve(gnd, tend=8e-9, timestep=2e-10, coupled_lte=True)
    t = np.asarray(res.sweep_values, float)
    vb = np.asarray(res.v('b'), float).reshape(-1)
    assert t[-1] >= 8e-9 * (1.0 - 1e-9)
    dev = float(np.max(np.abs(np.interp(tr, t, vb) - vr)))
    ## Measured 5.551e-16 at landing; 1e-12 leaves margin without letting a
    ## controller regression hide.
    assert dev < 1e-12, 'coupled+TLine drifted from standard: %.3e' % dev

    ## Mismatched RC load: must complete and settle at (1 + Gamma)/2 = 2/3.
    tran2 = Transient(line(True), toolkit=numeric, reltol=1e-4,
                      uic=True, timestep_max=2e-10)
    with quiet():
        res2 = tran2.solve(gnd, tend=8e-9, timestep=2e-10, coupled_lte=True)
    t2 = np.asarray(res2.sweep_values, float)
    vb2 = np.asarray(res2.v('b'), float).reshape(-1)
    assert t2[-1] >= 8e-9 * (1.0 - 1e-9)
    assert abs(vb2[-1] - 2.0 / 3.0) < 5e-3



def test_the_coupled_step_count_does_not_depend_on_how_tightly_newton_is_converged():
    """⚠⚠ THE `q^T dv0` TERM WAS COUNTED TWICE, and `vabstol = 1e-12` hid it.

    Eq (12) is `dh = -(f + q^T dv0)/denom` with the LTE residual `f` taken at the
    iterate BEFORE the Newton update; this code takes `err` at `x_stage1 = x +
    dx0`, where the update is already in it, and kept the term.  `denom = err
    w'/w` is tiny wherever the error is, so the spurious term decided the SIGN of
    `dh`: per time point the step grew 15 % on the first iteration and shrank
    15 % on the second (0.9775, on 8821 of 8828 points).  Measured, driven RC:

        vabstol     approx    bordered (before)    bordered (after)
        1e-12         93           93                  93
        1e-6         100         8828                 100      same error, 2.5e-4

    At 1e-12 the loop kept iterating until `dx0` had decayed to ~1e-13 and the
    term with it, so nothing showed until the default became 1e-6.  The property
    pinned is the one that was violated: a step controller's step COUNT must not
    hang on the Newton tolerance.
    """
    counts = {}
    for method in ('approx',):
        for va in (1e-12, 1e-9, 1e-6):
            tran = Transient(_rc(), toolkit=numeric, reltol=1e-5, vabstol=va)
            tran.par.coupled_method = method
            with quiet():
                res = tran.solve(tend=5e-4, timestep=1e-5, coupled_lte=True)
            t = np.asarray(res.v('b').x, dtype=float).ravel()
            v = np.asarray(res.v('b').y, dtype=float).ravel()
            err = float(np.max(np.abs(v - _analytic(t))[2:]))
            counts[method, va] = tran.statistics.accepted_steps
            assert err < 3.5e-4, (method, va, err)
    n = [counts['approx', va] for va in (1e-12, 1e-9, 1e-6)]
    assert max(n) <= 1.15 * min(n), n


def test_fang_coupled_stepping_refuses_a_stage_method_or_glm_by_name():
    """`coupled_lte=True` is built on a linear multistep companion.  Under
    radau, trbdf2, esdirk43 or a Nordsieck GLM it died inside
    `compute_derivatives` with a NotImplementedError that never named
    `coupled_lte` (measured 2026-09-24).  Fang's path solves the step from
    eq (6), a solution-space LTE over the LMM's step history; a stage
    method or GLM judges its step by its own embedded estimate.  So it is
    refused at entry, by name, and the LMMs run -- without the
    `ComplexWarning` every coupled run raised (`_state_row_mask` cast a
    complex `toMatrix(C)` with a zero imaginary part to float).  Kept
    refused on measurement (2026-09-25): a prototype around the stage step
    re-solved more often than the standard run rejects, and on the stiff
    RLC its error was flat across two decades of reltol (see the comment at
    the refusal)."""
    from pycircuit.circuit.integrator import (
        Gear2Integrator, TrapezoidalIntegrator, EulerIntegrator,
        RadauIIA3Integrator, TRBDF2Integrator, ESDIRK43Integrator,
        GLM2Integrator, GLM3Integrator)
    for integ in (RadauIIA3Integrator(), TRBDF2Integrator(),
                  ESDIRK43Integrator(), GLM2Integrator(), GLM3Integrator()):
        tran = Transient(_rc(), toolkit=numeric, integrator=integ)
        with pytest.raises(NotImplementedError, match='coupled_lte'):
            tran.solve(tend=5e-4, timestep=1e-5, coupled_lte=True)
    for integ in (Gear2Integrator(), TrapezoidalIntegrator(),
                  EulerIntegrator()):
        tran = Transient(_rc(), toolkit=numeric, reltol=1e-5, integrator=integ)
        with warnings.catch_warnings():
            warnings.simplefilter('error', np.exceptions.ComplexWarning)
            res = tran.solve(tend=5e-4, timestep=1e-5, coupled_lte=True)
        t = np.asarray(res.v('b').x, dtype=float).ravel()
        v = np.asarray(res.v('b').y, dtype=float).ravel()
        assert np.max(np.abs(v - _analytic(t))[2:]) < 5e-3, type(integ).__name__
