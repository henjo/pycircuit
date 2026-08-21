"""STAGE 12B, GATE 12B-1 -- ``p = df_ckt/dh``, against the residual it differentiates.

This is the gate that decides whether the coupled system is solving the right
problem. `p` is assembled from two independently-tested pieces -- the integrator's
`companion_dh` and the circuit's `dudt` -- so the thing that can still be wrong is
the ASSEMBLY: a missing term, a double-counted one, or a sign.

Checked by perturbing `h` and re-forming the WHOLE residual through
`solve_timestep`'s own code path, so the comparison is against what the solver
actually evaluates rather than against a re-derivation that could share a mistake
with the code under test.
"""
import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VSin, VPulse, IS
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import (EulerIntegrator, TrapezoidalIntegrator,
                                          Gear2Integrator)


def _driven_rc():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.5, freq=2e3, phase=25.0)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def _residual(tran, x, t, h):
    """The residual `solve_timestep` forms, at an explicit step size."""
    tran._dt = h
    q = tran.cir.q(x, tran.epar)
    Cm = tran.cir.C(x, tran.epar)
    iq, _geq = tran.get_diff(q, Cm)
    return (tran.cir.i(x, tran.epar) + iq
            + tran.cir.u(t, tran.epar, analysis=tran.par.analysis))


def _prepared(integrator, h=1.1e-6, h_last=7.0e-7):
    """A Transient with a plausible history, short of running a full solve."""
    cir = _driven_rc()
    tran = Transient(cir, toolkit=numeric, integrator=integrator)
    tran.irefnode = cir.get_node_index(gnd)
    tran.base_integrator = tran._get_integrator()
    n = cir.n

    rng = np.random.default_rng(20260801)
    x = 0.3 * rng.standard_normal(n)
    hist = max(2, tran.base_integrator.get_required_history())
    tran._qlast = np.array([cir.q(x + 0.05 * (k + 1), tran.epar)
                            for k in range(hist)])
    tran._iqlast = np.zeros((hist, n))
    tran._dt_last = h_last
    tran._dt_last2 = 9.0e-7
    tran._is_first_step = False
    tran._dt = h
    ## `get_diff` selects the active integrator, and `residual_dh` delegates to
    ## it -- so it must be primed exactly as a real step would.
    tran.get_diff(cir.q(x, tran.epar), cir.C(x, tran.epar))
    return tran, x, h


CASES = [('euler', EulerIntegrator), ('trapezoidal', TrapezoidalIntegrator),
         ('gear2', Gear2Integrator)]


@pytest.mark.parametrize('name,make', CASES, ids=[c[0] for c in CASES])
def test_p_matches_a_central_difference_of_the_whole_residual(name, make):
    tran, x, h = _prepared(make())
    t = 4.3e-4

    analytic = tran.residual_dh(x, t, h)

    eps = 1e-6 * h
    fp = _residual(tran, x, t + eps, h + eps)
    fm = _residual(tran, x, t - eps, h - eps)
    fd = (fp - fm) / (2 * eps)
    tran._dt = h

    scale = max(np.max(np.abs(fd)), 1e-12)
    assert np.max(np.abs(analytic - fd)) / scale < 1e-5, \
        '%s: analytic %r, finite difference %r' % (name, analytic, fd)
    ## Not vacuously zero -- both terms must actually be present.
    assert np.max(np.abs(analytic)) > 0.0


def test_the_source_term_is_actually_in_p():
    """Drop `du/dt` and `p` must change -- otherwise the test above is vacuous.

    On a driven circuit the source term is usually the LARGER contribution, so a
    `p` that silently omitted it would still look plausible.
    """
    tran, x, h = _prepared(Gear2Integrator())
    t = 4.3e-4
    full = tran.residual_dh(x, t, h)
    q = tran._q_at(x)
    integrator_only = tran.active_integrator.companion_dh(
        q, tran._qlast, h, tran._dt_last)
    assert np.max(np.abs(full - integrator_only)) > 0.0
    ## And it is not a rounding-level difference.
    assert np.max(np.abs(full - integrator_only)) > \
        1e-3 * np.max(np.abs(integrator_only))


def test_a_time_invariant_circuit_has_no_source_term():
    """With a constant source, `p` is the integrator term alone."""
    c = SubCircuit()
    c['is'] = IS(1, gnd, i=1e-3)
    c['R'] = R(1, gnd, r=1e3)
    c['C'] = C(1, gnd, c=1e-9)
    tran = Transient(c, toolkit=numeric)
    assert np.all(tran.cir.dudt(1e-4, tran.epar, analysis='tran') == 0.0)


def test_a_delay_line_dudt_is_the_derivative_of_its_u():
    """RE-DERIVED from a refusal test: `TLine.dudt` used to raise because the
    term was unwritten, and the raise was the honest alternative to inheriting
    the base class's silent zero.  It is written now -- the derivative of the
    SAME interpolation polynomial `u` evaluates (limiter included), so the two
    must agree with a central finite difference of `u` to rounding.  Populated
    history comes from a real standard-path run, not a synthetic ring.
    """
    from pycircuit.circuit.elements import TLine, VPulse
    import warnings
    c = SubCircuit()
    c.add_node('a'); c.add_node('b')
    c['vs'] = VPulse('s', gnd, v1=0.0, v2=1.0, td=1e-9, tr=2e-10,
                     tf=2e-10, pw=1e-8, per=1e-7)
    c['Rs'] = R('s', 'a', r=50.0)
    c['T'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
    c['R'] = R('b', gnd, r=50.0)
    tran = Transient(c, toolkit=numeric, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(gnd, tend=4e-9, timestep=2e-10)
    tl = tran.cir['T']
    assert len(tl.history) > 3
    worst = 0.0
    for t in (2.3e-9, 2.75e-9, 3.1e-9, 3.6e-9):
        eps = 1e-14
        fd = (np.asarray(tl.u(t + eps, analysis='tran'))
              - np.asarray(tl.u(t - eps, analysis='tran'))) / (2 * eps)
        an = np.asarray(tl.dudt(t, analysis='tran'))
        scale = max(float(np.max(np.abs(fd))), 1.0)
        worst = max(worst, float(np.max(np.abs(an - fd))) / scale)
    ## Measured 0.0 at landing; 1e-6 absorbs the FD's own rounding.
    assert worst < 1e-6, 'dudt drifted from d/dt of u: %.3e' % worst
