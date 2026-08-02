"""STAGE 13(13-3) -- PCNR: each limited quantity as its own unknown.

Aadithya, Keiter & Mei, *"Predictor/Corrector Newton-Raphson (PCNR)"* (Sandia).
The circuit used throughout is the paper's own Fig. 1: two parallel diodes with
very different saturation currents, fed through a resistor from a voltage source.
It is the example the paper introduces to show what goes wrong when two devices
limit the same branch voltage.
"""
import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit, defaultepar
from pycircuit.circuit.elements import R, VS, Diode
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit import pcnr as P


def _fig1(is1=1e-15, is2=1e-9, v=1.0, r=1e3):
    c = SubCircuit()
    c['vs'] = VS('e1', gnd, v=v)
    c['R'] = R('e1', 'e2', r=r)
    c['D1'] = Diode('e2', gnd, IS=is1)
    c['D2'] = Diode('e2', gnd, IS=is2)
    return c


def test_both_diodes_get_their_own_unknown_on_the_same_branch():
    """The clash the paper describes: two devices, one branch voltage."""
    c = _fig1()
    js = P.pcnr_junctions(c)
    assert len(js) == 2
    assert js[0][2:] == js[1][2:], 'the two junctions should share both rows'


def test_the_limited_unknowns_agree_with_the_branch_voltage_and_each_other():
    """GATE 13-3's declared success.

    `g_lim = v_k - (e_a - e_b)` must hold at convergence, for every device, which
    is exactly the consistency classic limiting cannot offer: there, each device
    sees whatever the last limiter left in the shared node.
    """
    c = _fig1()
    x, v_lim, _its = P.solve_dc(c, gnd)
    branch = float(x[c.get_node_index('e2')] - x[c.get_node_index(gnd)])
    assert v_lim[0] == pytest.approx(branch, abs=1e-9)
    assert v_lim[1] == pytest.approx(branch, abs=1e-9)
    assert v_lim[0] == pytest.approx(v_lim[1], abs=1e-12)


def test_pcnr_finds_the_same_operating_point_as_the_classic_solver():
    """Consistency is worth nothing if it changes the answer.

    The limiter only moves the point the next Jacobian is taken at; the converged
    solution must be a solution of the unmodified equations, so both paths must
    land in the same place.
    """
    c = _fig1()
    x, _v, _its = P.solve_dc(c, gnd)
    ref = DC(_fig1(), toolkit=numeric).solve()
    assert float(x[c.get_node_index('e2')]) == pytest.approx(float(ref.v('e2')),
                                                             rel=1e-7)


def test_the_currents_split_by_the_saturation_current_ratio():
    """A cross-check that is not just the residual: at one shared voltage the
    two diode currents must be in the ratio of their IS."""
    c = _fig1(is1=1e-15, is2=1e-12)
    _x, v_lim, _its = P.solve_dc(c, gnd)
    d1, d2 = c.elements['D1'], c.elements['D2']
    i1 = float(d1.pcnr_i(v_lim[0], {'IS': 1e-15}, defaultepar, numeric)[0])
    i2 = float(d2.pcnr_i(v_lim[1], {'IS': 1e-12}, defaultepar, numeric)[0])
    assert i2 / i1 == pytest.approx(1e3, rel=1e-6)


def test_the_solved_matrix_stays_the_size_of_the_original_mna_system():
    """The paper's third idea, and the whole reason the extra unknowns are
    affordable: `J_lim/lim` is the identity, so the border collapses into a
    rank-k update and the `Ax=b` actually solved is unchanged in size.

    Asserted on the matrix `predict` factorises, not on a timing.
    """
    c = _fig1()
    js = P.pcnr_junctions(c)
    x = np.zeros(c.n)
    v_lim = np.zeros(len(js))
    g_mna, g_lim, J_mm, J_ml, J_lm, _didv = P.augmented_system(c, x, v_lim, js)

    n = c.n
    assert J_mm.shape == (n, n)
    assert J_ml.shape == (n, len(js))
    assert J_lm.shape == (len(js), n)

    schur = J_mm - J_ml @ J_lm
    assert schur.shape == (n, n), \
        'the Schur complement grew the system to %r' % (schur.shape,)


def test_a_device_contributes_nothing_to_the_mna_jacobian_block():
    """Under PCNR the diode's conductance moves ENTIRELY into `J_MNA/lim`.

    Its current depends on its own unknown and not on the node voltages, so its
    contribution to `J_MNA/MNA` is zero -- which is what makes the two blocks add
    up to the same physics as before rather than double-counting it.
    """
    c = _fig1()
    js = P.pcnr_junctions(c)
    x = np.zeros(c.n)
    x[c.get_node_index('e2')] = 0.4
    v_lim = np.array([0.4, 0.4])
    _g, _gl, J_mm, J_ml, _J_lm, _dd = P.augmented_system(c, x, v_lim, js)

    ## The resistor and source still stamp; the diodes must not.
    e2 = c.get_node_index('e2')
    assert abs(J_mm[e2, e2] - 1e-3) < 1e-9, \
        'J_MNA/MNA at e2 is %g, so a diode is still stamping there' % J_mm[e2, e2]
    ## and their conductance is present in the other block instead
    assert np.max(np.abs(J_ml)) > 0.0


def test_refine_limits_each_device_independently():
    """The CORRECT phase: each device limits only what it owns.

    Choosing the operating point for this test needs care, and the first attempt
    got it wrong. `pnjlim` only depends on `IS` through the critical voltage
    `vc = VT log(VT / (1.414 IS))`, and only when `vold > 0` -- from `vold = 0` it
    takes the `VT log(vnew/VT)` branch, which is IS-independent, so two devices
    four decades apart limit to exactly the same value and the test asserted
    something the algorithm does not do.

    Here `vc` is 0.788 V for IS=1e-15 and 0.432 V for IS=1e-9, so a step to
    0.6 V is above one and below the other: one device limits, the other passes
    it through. That is the case where independence is observable.
    """
    c = _fig1(is1=1e-15, is2=1e-9)
    js = P.pcnr_junctions(c)
    v_old = np.array([0.5, 0.5])
    v_new = np.array([0.6, 0.6])
    out = P.refine(js, v_old, v_new)

    assert out[0] == pytest.approx(0.6), \
        'IS=1e-15 has vc=0.788 V, so 0.6 V should pass through unlimited'
    assert out[1] < 0.6, \
        'IS=1e-9 has vc=0.432 V, so 0.6 V should be limited'
    assert out[0] != out[1], \
        'the two devices produced the same value, so they are not independent'


## ---------------------------------------------------------------------------
## STAGE 13 -- PCNR on the transient path.
##
## The transient residual is `f = i(x) + iq + u(t)` with `J = G + Geq`, so the
## companion terms enter the coupled system as extra blocks. Everything else is
## the DC flow unchanged.

import warnings

from pycircuit.circuit.elements import C, VSin
from pycircuit.circuit.transient import Transient


def _rectifier():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b', IS=1e-15)
    c['R'] = R('b', gnd, r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def _run_tran(pcnr, reltol=1e-5, tend=2e-3):
    tran = Transient(_rectifier(), toolkit=numeric, reltol=reltol, pcnr=pcnr)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=1e-5)
    t = np.asarray(res.v('b').x, dtype=float).ravel()
    v = np.asarray(res.v('b').y, dtype=float).ravel()
    return t, v, tran.statistics.accepted_steps


def test_pcnr_and_limiting_agree_on_a_rectifier_waveform():
    """Consistency must not move the answer, on the transient path either."""
    ts, vs, _ = _run_tran(False)
    tp, vp, _ = _run_tran(True)
    d = np.abs(np.interp(ts, tp, vp) - vs)
    assert float(np.median(d)) < 1e-3, \
        'waveforms differ, median %g V' % float(np.median(d))


## RETIRED 2026-08-02 -- `test_pcnr_reaches_the_same_accuracy_in_far_fewer_steps`
## asserted `steps_pcnr < steps_limiting / 2`, which was TRUE ONLY BECAUSE OF A
## DEFECT.  `_solve_timestep_pcnr` handed the step controller `cir.G(x) + Geq`,
## and on the PCNR path `Diode.G` linearises around a `_vlim` that nothing
## updates, so that matrix had no diode conductance in it.  The controller formed
## `J^-1 Eg` from it and took steps 6.6x too large.
##
## The test therefore encoded the bug as a requirement, and would have blocked the
## fix.  Its replacement asserts the opposite and correct property -- that the two
## paths agree step for step, since they differ only in the iteration path and a
## converged Newton solution does not depend on the route taken to it.  See
## `test_gate_13_6_pcnr_and_limiting_take_the_same_steps` below.
##
## Kept as a comment rather than deleted because "PCNR is several times faster in
## transient" was a documented headline result, and a silent deletion leaves no
## trace of why it stopped being one.

def test_a_circuit_with_no_pcnr_device_falls_through():
    """`pcnr=True` on a linear circuit must solve it, not refuse it."""
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    tran = Transient(c, toolkit=numeric, reltol=1e-5, pcnr=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=2e-4, timestep=1e-5)
    assert tran.statistics.accepted_steps > 5
    assert np.max(np.abs(np.asarray(res.v('b').y, dtype=float))) > 0.1


def test_a_charge_storing_pcnr_device_is_refused():
    """A participating device with charge storage would contribute to `iq`/`Geq`
    too, and that term would have to move to the limited unknown as well.

    Left in the MNA block it would be evaluated at the NODE voltage while the
    current is evaluated at `v_lim` -- exactly the inconsistency PCNR exists to
    remove, reintroduced through the charge. Refused rather than silently wrong.
    """
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['D'] = Diode('a', 'b', IS=1e-15)
    c['R'] = R('b', gnd, r=1e3)

    js = P.pcnr_junctions(c)
    d = c.elements['D']
    ## Give the diode a non-zero capacitance stamp, which it does not have by
    ## default, and check the guard fires rather than quietly proceeding.
    d.C = lambda x, epar=None: np.array([[1e-12, -1e-12], [-1e-12, 1e-12]])
    x = np.zeros(c.n)
    with pytest.raises(NotImplementedError) as exc:
        P.augmented_system(c, x, np.zeros(len(js)), js,
                           J_extra=np.zeros((c.n, c.n)))
    assert 'charge storage' in str(exc.value)


## ---------------------------------------------------------------------------
## GATE 13-6.  The Jacobian `_solve_timestep_pcnr` hands to the step controller
## must be the one it actually solved.
##
## `Diode.G` linearises around `self._vlim`, written only by `Diode.limit` --
## which PCNR never calls, because limiting is what PCNR replaces.  `_vlim`
## therefore froze at its initial value (measured: 0.0 V across 2283 `G`
## evaluations while the junction swung -18.47 to 0.75 V), so `cir.G(x)` carried
## no diode conductance at all.
##
## It cancelled inside `augmented_system` -- added by `cir.G` and subtracted
## again by the per-device loop -- so the solved system and every converged
## voltage were right, and every DC test passed.  It escaped through exactly one
## door: the matrix returned to the step controller, which forms `J^-1 Eg`.
## These two tests watch that door from both sides.
## ---------------------------------------------------------------------------

def _mains_rectifier():
    """NOT `_rectifier` above: a 50 Hz 10 V mains rectifier with a 100 uF
    reservoir, whose diode swings -18.5 to +0.75 V. The wide reverse excursion is
    what makes a frozen `_vlim` show up as a 18.5 V discrepancy rather than a
    rounding difference."""
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit import gnd
    from pycircuit.circuit.elements import R, C, VSin, Diode
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=10.0, freq=50.0)
    c['D1'] = Diode(1, 2)
    c['C1'] = C(2, gnd, c=100e-6)
    c['R1'] = R(2, gnd, r=1e3)
    return c


def test_gate_13_6_controller_jacobian_carries_live_diode_conductance():
    """The returned J must hold IS/VT exp(v/VT) at the CONVERGED junction voltage.

    Fails against the defect by 12 orders of magnitude: frozen at v=0 the entry
    is IS/VT ~ 4e-12 S, while forward-biased it is ~1e-2 S.  Asserting a ratio
    rather than an absolute value keeps the test independent of `IS` and `VT`.
    """
    import warnings
    import numpy as np
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.transient import Transient

    cir = _mains_rectifier()
    tran = Transient(cir, toolkit=numeric, reltol=1e-5, pcnr=True)

    ia = cir.get_node_index('1')
    ib = cir.get_node_index('2')
    seen = []

    real = Transient._solve_timestep_pcnr

    def spy(self, *a, **k):
        x, feval, J, f = real(self, *a, **k)
        v = float(x[ia] - x[ib])
        ## The diode's own contribution to the node-2 diagonal.  Node 2 also
        ## carries R1 and C1's companion, so compare the OFF-DIAGONAL entry,
        ## which only the diode writes.
        seen.append((v, -float(J[ib, ia])))
        return x, feval, J, f

    Transient._solve_timestep_pcnr = spy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tran.solve(tend=0.01, timestep=1e-4)
    finally:
        Transient._solve_timestep_pcnr = real

    assert seen, 'the PCNR timestep path never ran'

    ## VT is k*T/q at the *analysis* temperature, so recompute it exactly as
    ## `pcnr_didv` does rather than hardcoding a value -- a constant here was
    ## wrong by 2.2% in `g`, which at this v is a 0.5 mV error in VT.
    from pycircuit.circuit.circuit import defaultepar
    epar = tran.epar if hasattr(tran, 'epar') else defaultepar
    VT = numeric.kboltzmann * epar.T / numeric.qelectron
    IS = cir['D1'].iparv.IS

    ## Check where the conductance is large enough to be unambiguous; near
    ## reverse bias every model agrees on "approximately zero".
    checked = 0
    for v, g in seen:
        g_exact = IS * np.exp(v / VT) / VT
        if g_exact < 1e-6:
            continue
        assert g == pytest.approx(g_exact, rel=1e-6), (
            'J[2,1] = %.6e but the diode at v=%.6f V has g=%.6e; the Jacobian '
            'given to the step controller is not the one that was solved' %
            (g, v, g_exact))
        checked += 1
    assert checked > 20, 'only %d forward-biased points -- test is not exercising it' % checked


def test_gate_13_6_pcnr_and_limiting_take_the_same_steps():
    """PCNR and classic limiting must agree -- step for step, not merely closely.

    They differ only in the ITERATION PATH; a converged Newton solution is
    whatever the residual says it is, independent of the route.  Equal solutions
    give equal LTE and hence equal step sequences, so any divergence in step
    count is a defect signature.  Against the gate 13-6 defect this reported
    4.1-6.8x fewer steps for PCNR.
    """
    import warnings
    import numpy as np
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.transient import Transient

    out = {}
    for pcnr in (False, True):
        tran = Transient(_mains_rectifier(), toolkit=numeric, reltol=1e-5, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=0.01, timestep=1e-4)
        w = res.v(2, gnd)
        out[pcnr] = (tran.statistics.accepted_steps,
                     np.asarray(w.x, dtype=float).ravel(),
                     np.asarray(w.y, dtype=float).ravel())

    n_lim, t_lim, v_lim = out[False]
    n_pc, t_pc, v_pc = out[True]

    assert n_pc == n_lim, ('PCNR took %d steps, limiting %d -- the two paths '
                           'solve the same equations and must agree' % (n_pc, n_lim))
    assert len(t_pc) == len(t_lim)
    assert np.allclose(t_pc, t_lim, rtol=1e-9, atol=0.0), 'step sequences differ'
    assert np.max(np.abs(v_pc - v_lim)) < 1e-9 * max(1.0, np.max(np.abs(v_lim)))


def test_gate_13_7_pcnr_is_honoured_on_the_coupled_path():
    """`pcnr=True` must not be silently ignored when `coupled_lte=True`.

    `_solve` dispatches on `coupled_lte` before it looks at `pcnr`, so the
    coupled path used to run the classic limiter whatever `pcnr` said -- 0 PCNR
    steps against 4869 `Diode.limit` calls, and results bit-identical to
    `pcnr=False`. A parameter that silently does nothing is the same class of
    defect as an orphaned `h_next`: nothing fails, the run is simply not the one
    that was asked for.

    Asserted by counting CALLS rather than by comparing waveforms, because the
    two paths agree to within their own truncation error -- which is the correct
    outcome, and therefore cannot distinguish "PCNR ran" from "PCNR was ignored".
    """
    import warnings
    from pycircuit.circuit import numeric
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.elements import Diode
    import pycircuit.circuit.pcnr as _pcnr

    counts = {}
    for pcnr in (False, True):
        hits = {'limit': 0, 'aug': 0}
        real_limit, real_aug = Diode.limit, _pcnr.augmented_system

        def spy_limit(self, *a, **k):
            hits['limit'] += 1
            return real_limit(self, *a, **k)

        def spy_aug(*a, **k):
            hits['aug'] += 1
            return real_aug(*a, **k)

        Diode.limit, _pcnr.augmented_system = spy_limit, spy_aug
        try:
            tran = Transient(_mains_rectifier(), toolkit=numeric,
                             reltol=1e-5, pcnr=pcnr)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tran.solve(tend=2e-3, timestep=1e-4, coupled_lte=True)
        finally:
            Diode.limit, _pcnr.augmented_system = real_limit, real_aug
        counts[pcnr] = hits

    assert counts[False]['aug'] == 0, \
        'pcnr=False assembled the PCNR system %d times' % counts[False]['aug']
    assert counts[False]['limit'] > 10, \
        'pcnr=False did not use the classic limiter at all'

    assert counts[True]['aug'] > 10, \
        'pcnr=True on the coupled path assembled the PCNR system only %d times ' \
        '-- the parameter is being ignored' % counts[True]['aug']
    ## The DC operating point still runs the classic path, so this is "no
    ## limiting in the transient loop", not "never called".
    assert counts[True]['limit'] < counts[False]['limit'] / 100, \
        'pcnr=True still called the classic limiter %d times against %d ' \
        '-- device limiting is active under PCNR' \
        % (counts[True]['limit'], counts[False]['limit'])
