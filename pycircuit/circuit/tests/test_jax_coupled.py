"""Fang's coupled (x, h) stepping on the JAX backend (parity item P19).

The 'approx' branch of the CPU's fang_timestep, traced.  Measured at landing:
rc step counts 127/127 CPU/JAX, cross-backend deviation 6.8e-3 (rc) and
4.3e-3 (vsin) at reltol 1e-4, JAX-coupled MORE accurate than CPU-coupled at
every tolerance probed (1.0e-3/1.0e-3/2.9e-4 vs 1.45e-2/5.8e-3/1.9e-3 at
reltol 1e-3/1e-4/1e-5 on the uic RC), and per-lane batched coupled errors of
1.3e-3 / 1.1e-5 against analytic references.

The vsin STEP-COUNT comparison is deliberately absent: CPU takes ~3x the
steps there because Sin.next_event fires quarter-period breakpoints on the
CPU and returns inf on this backend -- a breakpoint-feature divergence
(parity doc P14 note), not a property of the coupled method.
"""

import warnings

import numpy as np
import pytest

jax = pytest.importorskip('jax')
import jax.numpy as jnp

from pycircuit.circuit import circuit as circuit_mod, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin


def _build(kind='rc'):
    c = SubCircuit()
    c.add_node('in'); c.add_node('out')
    c['V1'] = VS('in', gnd, v=1.0) if kind == 'rc' else \
        VSin('in', gnd, va=1.0, freq=1e3)
    c['R1'] = R('in', 'out', r=1e3)
    c['C1'] = C('out', gnd, c=1e-6)
    return c


def _jax_coupled(kind, uic, reltol):
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(_build(kind), reltol=reltol, coupled_lte=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=5e-3, timestep=1e-4, uic=uic)
        return (np.asarray(res.sweep_values, float),
                np.asarray(res.v('out'), float).reshape(-1))
    finally:
        circuit_mod.default_toolkit = saved


def _cpu_coupled(kind, uic, reltol):
    from pycircuit.circuit.transient import Transient
    tran = Transient(_build(kind), toolkit=numeric, reltol=reltol, uic=uic)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=5e-3, timestep=1e-4, coupled_lte=True)
    return (np.asarray(res.sweep_values, float),
            np.asarray(res.v('out'), float).reshape(-1))


@pytest.mark.parametrize('kind,uic,bound', [('rc', True, 2e-2),
                                            ('vsin', False, 2e-2)])
def test_jax_coupled_matches_cpu_coupled(kind, uic, bound):
    tc, vc = _cpu_coupled(kind, uic, 1e-4)
    tj, vj = _jax_coupled(kind, uic, 1e-4)
    for t in (tc, tj):
        assert t[0] == 0.0
        assert t[-1] == pytest.approx(5e-3, rel=1e-12)   # F3's guarantee holds
        assert np.all(np.diff(t) > 0)
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < bound, 'coupled backends diverge: %.3e' % dev
    if kind == 'rc':
        ## no breakpoint divergence on this drive: step parity is tight
        ## (127/127 measured at landing)
        assert 0.8 < len(tc) / len(tj) < 1.25


def test_jax_coupled_reltol_governs():
    errs = []
    for rt in (1e-3, 1e-5):
        tj, vj = _jax_coupled('rc', True, rt)
        errs.append(float(np.max(np.abs(vj - (1 - np.exp(-tj / 1e-3))))))
    assert errs[1] < errs[0], \
        'two decades of reltol did not move the error: %r' % errs


def test_jax_coupled_batched_lanes():
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(_build('rc'), reltol=1e-4, coupled_lte=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve_batched(
                gnd, override_params_tree={'R': {'r': jnp.array([[1e2], [1e4]])}},
                tend=5e-3, timestep=1e-4, uic=True)
        for b, tau, bound in ((0, 1e-4, 5e-3), (1, 1e-2, 1e-3)):
            t = np.asarray(res[b].sweep_values, float)
            v = np.asarray(res[b].v('out'), float).reshape(-1)
            assert t[0] == 0.0 and np.all(np.diff(t) > 0)
            assert t[-1] == pytest.approx(5e-3, rel=1e-12)
            err = float(np.max(np.abs(v - (1 - np.exp(-t / tau)))))
            assert err < bound, 'lane %d err %.2e' % (b, err)
    finally:
        circuit_mod.default_toolkit = saved


def test_jax_coupled_band_sentinel():
    """The 'auto' band resolves to Fang's values; explicit values survive --
    the F5 contract, on this backend."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        assert JAXTransient(_build())._coupled_band() == (0.7, 3.0, 0.15)
        assert JAXTransient(_build(), lte_gamma_min=0.0, lte_gamma_max=1.0,
                            lte_eta=None)._coupled_band() == (0.0, 1.0, None)
    finally:
        circuit_mod.default_toolkit = saved


def test_jax_coupled_tline_matches_cpu_standard():
    """RE-DERIVED from a refusal test (this used to expect
    NotImplementedError).  The refusal's cause was a livelock, and the fix
    took three composed pieces, each earned by a traced measurement:

    - a from-zero kink makes every signal scale together, so the coupled
      LTE's dev is proportional to its ref and h CANCELS -- err =
      1/(TRTOL*reltol) = 1428.6 measured, constant from h=4.4e-11 down to
      the excursion floor.  No step size can satisfy the band; the outer
      retry collapses dt to dt_min.
    - one graced step is not enough: the NEXT step extrapolates degree-2
      through a PRE-kink history point, and a quadratic through two flat
      points and one ramp point misses a line by err = 139.2e11*h --
      in-band only below the within-point floor.  Hence the coupled accept
      path ZEROES h_history on a breakpoint landing (band skipped one step,
      degree held at 1 for two, cold-start semantics).
    - the far-end wavefront onset is the same kink NOT at any element-
      reported breakpoint, so collect_breakpoints registers source corners
      + k*TD arrivals.  (Arrivals alone were tried first and falsified --
      a registered kink without the history reset still livelocks.)

    The gate: the pulsed matched line that livelocked at every tend past
    the first arrival now lands on tend and matches the CPU standard path
    to 5e-16 (measured 4.996e-16 at tend=8e-9, 112 points)."""
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.integrator import Gear2Integrator
    from pycircuit.circuit.elements import TLine, VPulse

    def line():
        c = SubCircuit()
        c.add_node('a'); c.add_node('b')
        c['V1'] = VPulse('s', gnd, v1=0.0, v2=1.0, td=1e-9, tr=2e-10,
                         tf=2e-10, pw=1e-8, per=1e-7)
        c['Rs'] = R('s', 'a', r=50.0)
        c['T1'] = TLine('a', gnd, 'b', gnd, Z0=50.0, TD=1e-9)
        c['Rl'] = R('b', gnd, r=50.0)
        return c

    cpu = Transient(line(), toolkit=numeric, reltol=1e-4,
                    integrator=Gear2Integrator(), uic=True,
                    timestep_max=2e-10)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res_c = cpu.solve(gnd, tend=8e-9, timestep=2e-10)
    tc = np.asarray(res_c.sweep_values, float)
    vc = np.asarray(res_c.v('b'), float).reshape(-1)

    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(line(), reltol=1e-4, coupled_lte=True,
                            timestep_max=2e-10)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=8e-9, timestep=2e-10, uic=True)
        tj = np.asarray(res.sweep_values, float)
        vj = np.asarray(res.v('b'), float).reshape(-1)
    finally:
        circuit_mod.default_toolkit = saved

    ## The run must reach tend (the livelock died mid-ramp at ~1.05e-9)...
    assert tj[-1] > 8e-9 * (1.0 - 1e-6)
    ## ...must carry the delayed wavefront (the DC-stamp defect gave 24.5 V
    ## garbage; a dropped history gives vb == 0 on a matched line)...
    vb_delay = tj[np.argmax(vj > 0.25)]
    assert 1.0e-9 < vb_delay < 2.5e-9
    ## ...and must match the CPU standard path.  5e-13 leaves three orders
    ## over the measured 5e-16 without letting a controller regression hide.
    dev = float(np.max(np.abs(vc - np.interp(tc, tj, vj))))
    assert dev < 5e-13, 'coupled+TLine drifted from CPU: %.3e' % dev


# ---------------------------------------------------------------------------
# The coupled path's scope, and one mismatch inside it (2026-09-01)
# ---------------------------------------------------------------------------

def _euler_coupled_pair(reltol=1e-6, tend=2e-5, ts=2e-7):
    """The same circuit, coupled, integrated with EULER on both backends."""
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit.integrator import EulerIntegrator
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.jaxtransient import JAXTransient

    def build(tk):
        _cm.default_toolkit = tk
        c = SubCircuit(toolkit=tk)
        c['vs'] = VSin('a', gnd, va=1.0, freq=1e5, toolkit=tk)
        c['R'] = R('a', 'b', r=1e3, toolkit=tk)
        c['C'] = C('b', gnd, c=1e-9, toolkit=tk)
        return c

    saved = _cm.default_toolkit
    try:
        cpu = Transient(build(numeric), toolkit=numeric,
                        integrator=EulerIntegrator(), reltol=reltol, uic=True)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rc = cpu.solve(tend=tend, timestep=ts, coupled_lte=True)
        tc = np.asarray(rc.v('b').x[0], float)
        vc = np.asarray(rc.v('b').y, float)

        j = JAXTransient(build(jaxtoolkit), coupled_lte=True,
                         integrator='euler', reltol=reltol)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rj = j.solve(gnd, tend=tend, timestep=ts, uic=True)
        tj = np.asarray(rj.sweep_values, float).reshape(-1)
        vj = np.asarray(rj.v('b'), float).reshape(-1)
    finally:
        _cm.default_toolkit = saved
    return tc, vc, tj, vj


def test_euler_coupled_uses_an_euler_dh_derivative():
    """The Fang `p = d(iq)/dh` must describe the method the residual used.

    It was hardcoded euler-or-gear2 while the residual is assembled with
    `method=eval_method`, so `integrator='euler'` with `coupled_lte=True`
    integrated with Euler and then differentiated a GEAR-2 companion with
    respect to h on every non-first-order step.  Fang's step-size Newton was
    solving a slightly wrong equation -- and it still converged on x, which
    is exactly why nothing failed and the mismatch survived.

    Measured across the fix on this circuit (2026-09-01):

        before   max|jax - cpu| = 9.1162e-03 V (1.058% of span), 16522 steps
        after    max|jax - cpu| = 3.8516e-03 V (0.447% of span),  1011 steps

    The step count is the louder half: a wrong dh-derivative made the
    step-size Newton flail, taking 16x the steps to agree half as well.  The
    bounds below sit above the measured values with room, and would both
    have failed before the fix.
    """
    tc, vc, tj, vj = _euler_coupled_pair()
    span = float(np.max(np.abs(vc)))
    dev = float(np.max(np.abs(vj - np.interp(tj, tc, vc))))
    assert dev < 6e-3, 'euler+coupled drifted from the CPU: %.4e V' % dev
    assert dev < 0.01 * span
    ## the flailing signature, pinned independently of the waveform
    assert len(tj) < 4000, 'step count regressed to the flailing regime: %d' \
        % len(tj)


def test_vector_pcnr_runs_inside_the_coupled_path():
    """The DEVICE view of PCNR, inside Fang -- wired 2026-09-01.

    PCNR-inside-Fang understood only the PAIR view `(ra, rb, IS, VT, fns)`.
    Handed the device-view 3-tuple `('vector', devices, epar)` it used the
    literal string 'vector' as an array index and died with "JAX does not
    support string indexing" several frames deep.  Since every real compact
    model takes the device view, that was reachable by anyone combining a
    MOSFET with `coupled_lte`; it was refused at setup as a stopgap, and is
    now actually implemented.

    What made the port small: the Schur-reduced `(F_eff, J_eff)` is an
    n-sized system whose Newton step equals predict's `dx_mna`, so fang's
    machinery -- eq (18)'s second solve against the same factors included --
    works on it unchanged.  That argument never depended on WHICH view
    produced the blocks, which is why the device view drops in.

    Measured at landing against the CPU (EKV NMOS, uic, reltol 1e-8/1e-9):
    standard-path vector PCNR 2.000e-05, coupled 2.030e-05 -- the coupled
    path is as accurate as the path already trusted, which is the claim.
    """
    import pycircuit.circuit.circuit as _cm
    from pycircuit.circuit import elements_hdl as _eh
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    from pycircuit.circuit.transient import Transient

    def build(tk):
        _cm.default_toolkit = tk
        c = SubCircuit(toolkit=tk)
        nd, ng, nv = c.add_node('d'), c.add_node('g'), c.add_node('vdd')
        c['vg'] = VSin(ng, gnd, va=0.4, vo=1.1, freq=1e6, toolkit=tk)
        c['vd'] = VS(nv, gnd, v=1.8, toolkit=tk)
        c['Rd'] = R(nv, nd, r=1e4, toolkit=tk)
        c['Cd'] = C(nd, gnd, c=1e-13, toolkit=tk)
        c['M'] = _eh.EkvNmosHdl(nd, ng, gnd, gnd, toolkit=tk)
        return c

    tend, ts = 1e-6, 2e-8
    saved = _cm.default_toolkit
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            rc = Transient(build(numeric), toolkit=numeric, uic=True,
                           reltol=1e-9).solve(refnode=gnd, tend=tend,
                                              timestep=ts)
        tc = np.asarray(rc.v('d').x[0], float)
        vc = np.asarray(rc.v('d').y, float)

        _cm.default_toolkit = jaxtoolkit
        tran = JAXTransient(build(jaxtoolkit), pcnr=True, reltol=1e-8,
                            coupled_lte=True)
        ## the routing this test exists for
        meta, _vt = tran._pcnr_setup()
        assert meta[0] == 'vector', 'test no longer covers the device view'
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(refnode=gnd, tend=tend, timestep=ts, uic=True)
        t = np.asarray(res.sweep_values, float).reshape(-1)
        v = np.asarray(res.v('d'), float).reshape(-1)
    finally:
        _cm.default_toolkit = saved

    assert len(t) > 10
    assert np.all(np.isfinite(v))
    span = float(np.max(np.abs(vc)))
    dev = float(np.max(np.abs(v - np.interp(t, tc, vc))))
    assert dev < 1e-3 * span, \
        'coupled vector PCNR drifted from the CPU: %.3e (span %.3f)' % (dev,
                                                                        span)


# ---------------------------------------------------------------------------
# grid_locked: the caller's grid wins over the method's step-size solve
# ---------------------------------------------------------------------------

_TD, _TR, _PW, _TF, _PER = 1e-5, 1e-7, 2e-5, 1e-7, 5e-5
_STEP, _GTEND = 1e-6, 2e-5


def _pulsed_rc(tk=None):
    from pycircuit.circuit.elements import VPulse
    kw = {'toolkit': tk} if tk is not None else {}
    c = SubCircuit(**kw)
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=_TD, tr=_TR, tf=_TF,
                     pw=_PW, per=_PER, **kw)
    c['R'] = R('a', 'b', r=1e3, **kw)
    c['C'] = C('b', gnd, c=1e-9, **kw)
    return c


def _jax_pulsed(fixed):
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient
    saved = circuit_mod.default_toolkit
    circuit_mod.default_toolkit = jaxtoolkit
    try:
        tran = JAXTransient(_pulsed_rc(jaxtoolkit), coupled_lte=True,
                            reltol=1e-5)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=_GTEND, timestep=_STEP,
                             fixed_timestep=fixed)
        return (np.asarray(res.sweep_values, float).reshape(-1),
                np.asarray(res.v('b'), float).reshape(-1))
    finally:
        circuit_mod.default_toolkit = saved


def test_fixed_timestep_keeps_the_grid_on_the_coupled_path():
    """The JAX mirror of the CPU's gate 12-4 test.

    `fixed_timestep` and Fang's method are not in conflict: one exists to
    CHOOSE the step size, the other says the caller already has.  So the grid
    is kept and the LTE equation is dropped on every step, exactly as for a
    breakpoint-truncated one -- the circuit is still solved coupled, it just
    has nothing left to solve for.

    The whole semantic content is one flag: an over-band HELD step is normally
    reported unconverged so the caller shrinks and retries, and under a locked
    grid shrinking is not an option we have.  On the CPU, conflating the two
    made the retry shrink h and the uniform grid disappear.
    """
    t, _v = _jax_pulsed(True)
    dt = np.diff(t)
    assert np.allclose(dt, _STEP, rtol=1e-9), \
        'grid was not uniform: steps ranged %g .. %g' % (dt.min(), dt.max())


def test_fixed_timestep_and_adaptive_differ_on_the_coupled_path():
    """Guard against the test above passing because nothing adapts anyway."""
    t, _v = _jax_pulsed(False)
    dt = np.diff(t)
    assert not np.allclose(dt, _STEP, rtol=1e-9), \
        'the adaptive coupled path produced a uniform grid, so the ' \
        'fixed-step test above proves nothing'


def test_fixed_grid_coupled_matches_the_cpu():
    """And the waveform on that grid is the CPU's.

    Measured at landing: 21 steps on both backends, max deviation 4.4e-16 --
    the two solvers are doing the same arithmetic on the same points, which
    is the strongest form this comparison can take (no interpolation).
    """
    from pycircuit.circuit.transient import Transient

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        rc = Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-5).solve(
            tend=_GTEND, timestep=_STEP, coupled_lte=True,
            fixed_timestep=True)
    vc = np.asarray(rc.v('b').y, float).ravel()

    tj, vj = _jax_pulsed(True)
    assert len(tj) == len(vc), 'step counts differ: %d vs %d' % (len(tj),
                                                                 len(vc))
    dev = float(np.max(np.abs(vj - vc)))
    assert dev < 1e-12, 'fixed-grid coupled backends disagree: %.3e' % dev
