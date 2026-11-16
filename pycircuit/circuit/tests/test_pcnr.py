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
    i1 = float(d1.pcnr_i(v_lim[[0]], {'IS': 1e-15}, defaultepar, numeric)[0])
    i2 = float(d2.pcnr_i(v_lim[[1]], {'IS': 1e-12}, defaultepar, numeric)[0])
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


## ---------------------------------------------------------------------------
## GATE 13-8 -- the PCNR devices leave the ordinary assembly.
##
## `augmented_system` used to assemble every device and then subtract the PCNR
## ones out again.  That survived floating point only because `Diode` read a
## stale `_vlim` and so presented a tiny companion; an honest diode at an
## unlimited node voltage overflows `exp` and the subtraction becomes nan
## (benchmarks/transient_review/pcnr_narrow_fix.py).  These tests pin the
## replacement: a `skip` set that assembly honours, and refuses where it cannot.
## ---------------------------------------------------------------------------

def test_gate_13_8_pcnr_devices_are_not_assembled():
    """Counted, not inferred: the two forms agree to rounding, so no answer can
    tell "skipped" from "assembled and subtracted".  Entry: 2 calls of each per
    device per iteration (fig1: 32 and 32)."""
    hits = {'G': 0, 'i': 0}
    real_G, real_i = Diode.G, Diode.i

    def spy_G(self, *a, **k):
        hits['G'] += 1
        return real_G(self, *a, **k)

    def spy_i(self, *a, **k):
        hits['i'] += 1
        return real_i(self, *a, **k)

    Diode.G, Diode.i = spy_G, spy_i
    try:
        x, _v, its = P.solve_dc(_fig1(), gnd)
    finally:
        Diode.G, Diode.i = real_G, real_i
    assert its > 1
    assert hits == {'G': 0, 'i': 0}, \
        'PCNR assembled its own devices through Diode.G/i: %r' % hits


def test_gate_13_8_skip_leaves_out_exactly_the_named_instances():
    """The reference is built by stamping the REMAINING elements directly, not
    by subtracting D2 from the full assembly -- this test's first version did
    that, and failed at 8.9e-11 relative: D1's 1.6e-20 S does not survive
    `(D1 + D2) - D2` when D2 is larger.  Which is the 13-8 point in miniature."""
    c = _fig1()
    x = np.array([0.0, 1.0, 0.62, -1e-3])[:c.n]
    G_ref = np.zeros((c.n, c.n))
    i_ref = np.zeros(c.n)
    for name, element in c.elements.items():
        if name == 'D2':
            continue
        rows = c.elementnodemap[name]
        G_ref[np.ix_(rows, rows)] += np.asarray(element.G(x[rows], defaultepar))
        i_ref[rows] += np.asarray(element.i(x[rows], defaultepar))
    G_skip = np.asarray(c.G(x, defaultepar, skip={'D2'}), dtype=float)
    i_skip = np.asarray(c.i(x, defaultepar, skip={'D2'}), dtype=float)
    np.testing.assert_allclose(G_skip, G_ref, rtol=1e-15, atol=1e-300)
    np.testing.assert_allclose(i_skip, i_ref, rtol=1e-15, atol=1e-300)
    ## And an empty skip is the ordinary assembly, exactly.
    assert np.array_equal(np.asarray(c.G(x, defaultepar, skip=frozenset())),
                          np.asarray(c.G(x, defaultepar)))


def test_gate_13_8_skip_of_an_unknown_instance_is_an_error():
    """A typo would otherwise skip nothing and look like it worked."""
    with pytest.raises(KeyError):
        _fig1().G(np.zeros(4), defaultepar, skip={'D3'})


## (`test_gate_13_8_skip_refuses_a_batched_class` lived here.  13-8 refused a
## skip of a batched class; 13-11 implements it, so the refusal test became the
## three gate 13-11 tests at the end of this module.)


## ---------------------------------------------------------------------------
## GATE 13-9 -- `Diode` reads `_vlim` only inside a classic-limiting scope.
##
## Outside the three Newton loops that call `cir.limit`, nothing moves `_vlim`,
## so every other caller of `cir.G`/`cir.i` used to get a frozen companion: a
## valid tangent of a linearisation nobody asked for.  13-6 met it as a 6.6x
## step error; 13-9 met it as PSS returning 0.000 V for a rectifier.
## ---------------------------------------------------------------------------

def test_gate_13_9_cir_G_is_honest_after_a_classic_dc_op():
    """Entry: after the DC op `_vlim` = 0.0 and `cir.G` carried g(0 V) =
    3.87e-12 S at EVERY voltage -- 15.5 S was the right answer at 0.75 V."""
    import warnings
    c = _mains_rectifier()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        DC(c, toolkit=numeric).solve()
    d = c['D1']
    rows = c.elementnodemap['D1']
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    for v in (-1.0, 0.3, 0.6, 0.75):
        x = np.zeros(c.n)
        x[rows[0]] = v
        G = np.asarray(c.G(x, defaultepar), dtype=float)
        i = np.asarray(c.i(x, defaultepar), dtype=float)
        g_true = d.iparv.IS * np.exp(v / VT) / VT
        i_true = d.iparv.IS * (np.exp(v / VT) - 1.0)
        assert G[rows[0], rows[0]] == pytest.approx(g_true, rel=1e-12), v
        assert i[rows[0]] == pytest.approx(i_true, rel=1e-12, abs=1e-30), v


def test_gate_13_9_diode_writes_no_state_outside_limiting():
    """A plain evaluation must not plant the `_vlim` the next one would read."""
    d = Diode(1, 2, toolkit=numeric)
    d.G(np.array([0.3, 0.0]), defaultepar)
    d.i(np.array([0.6, 0.0]), defaultepar)
    assert '_vlim' not in d.__dict__
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    assert float(d.G(np.array([0.6, 0.0]), defaultepar)[0, 0]) == \
        pytest.approx(d.iparv.IS * np.exp(0.6 / VT) / VT, rel=1e-12)


def test_gate_13_9_inside_the_scope_the_limiter_state_is_used():
    """The classic path is unchanged: inside the scope `G` is taken at `_vlim`."""
    from pycircuit.circuit._limiting import classic_limiting
    d = Diode(1, 2, toolkit=numeric)
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    with classic_limiting(d):
        d.G(np.array([0.0, 0.0]), defaultepar)          # seeds _vlim = 0
        d.limit(np.array([5.0, 0.0]), np.array([0.0, 0.0]), defaultepar)
        g = float(d.G(np.array([5.0, 0.0]), defaultepar)[0, 0])
    assert 0.0 < d._vlim < 5.0, 'the limiter did not hold the junction back'
    assert g == pytest.approx(d.iparv.IS * np.exp(d._vlim / VT) / VT, rel=1e-12)


def test_gate_13_9_scope_is_depth_counted():
    from pycircuit.circuit._limiting import classic_limiting
    c = _fig1()
    with classic_limiting(c):
        with classic_limiting(c):
            pass
        assert c['D1']._classic_limiting == 1, 'an inner scope ended the outer one'
    assert c['D1']._classic_limiting == 0
    with pytest.raises(RuntimeError):
        c['D1'].set_classic_limiting(False)


## ---------------------------------------------------------------------------
## GATE 13-11 -- `skip` on the batched (JAX) path.
##
## 13-8 refused it, which made PCNR on a JAX-toolkit circuit fail outright --
## and before 13-8 PCNR on JAX went nan on any overshoot, because autodiff makes
## the diode honest and the old add-then-subtract then cancels inf.  These need
## JAX and are skipped without it.
## ---------------------------------------------------------------------------

def _jax_circuit(builder):
    from pycircuit.circuit import circuit as cm
    from pycircuit.circuit.toolkit import jaxtoolkit
    saved = cm.default_toolkit
    cm.default_toolkit = jaxtoolkit
    try:
        c = builder()
        c.toolkit = jaxtoolkit
        c.update_iparv()
        return c
    finally:
        cm.default_toolkit = saved


def _one_amp():
    from pycircuit.circuit.elements import IS as ISrc
    c = SubCircuit()
    c['is'] = ISrc(gnd, 1, i=1.0)
    c['D1'] = Diode(1, gnd)
    return c


def test_gate_13_11_pcnr_runs_on_a_jax_circuit():
    """Including the 1 A case, which went nan here before 13-8 and raised after."""
    pytest.importorskip('jax')
    for builder in (_fig1, _one_amp):
        cj = _jax_circuit(builder)
        assert cj._eval_groups, 'nothing is batched; the test proves nothing'
        xj, _vj, itj = P.solve_dc(cj, gnd)
        xn, _vn, itn = P.solve_dc(builder(), gnd)
        assert itj == itn
        np.testing.assert_allclose(np.asarray(xj, dtype=float), xn,
                                   rtol=1e-12, atol=1e-15)


def test_gate_13_11_skipped_instances_are_not_evaluated_in_a_batch():
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit
    sizes = []
    real = jaxtoolkit.generate_batched_eval

    def spy(cls, method='i'):
        f = real(cls, method)

        def g(X, p, epar):
            sizes.append((cls.__name__, int(np.shape(X)[0])))
            return f(X, p, epar)
        return g

    jaxtoolkit.generate_batched_eval = spy
    try:
        P.solve_dc(_jax_circuit(_fig1), gnd)
    finally:
        jaxtoolkit.generate_batched_eval = real
    assert not [s for s in sizes if s[0] == 'Diode'], sizes


def test_gate_13_11_batched_skip_matches_per_element_skip():
    """A PARTIAL skip -- some instances of a class kept -- is the case where node
    maps, parameters and scatter indices must be sliced alike.  Checked against
    the numeric toolkit's per-element loop, on every stamp."""
    pytest.importorskip('jax')
    from pycircuit.circuit.elements import C as Cap

    def build():
        cir = SubCircuit()
        cir['R1'] = R('a', 'b', r=1e3)
        cir['D1'] = Diode('a', gnd, IS=1e-14)
        cir['D2'] = Diode('b', gnd, IS=3e-13)
        cir['D3'] = Diode('a', 'b', IS=2e-12)
        cir['C1'] = Cap('a', gnd, c=1e-9)
        return cir

    batched, per_element = _jax_circuit(build), build()
    assert Diode in batched._eval_groups
    x = np.array([0.61, 0.33, 0.0])[:batched.n]
    for skip in ({'D2'}, {'D1', 'D3'}, {'D1', 'D2', 'D3'}, {'C1'}, {'D3', 'C1'}):
        for method in ('G', 'i'):
            a = np.asarray(getattr(batched, method)(x, defaultepar, skip=skip),
                           dtype=float)
            b = np.asarray(getattr(per_element, method)(x, defaultepar, skip=skip),
                           dtype=float)
            np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-18,
                                       err_msg='%s skip=%r' % (method, skip))
    ## No skip is exactly the ordinary batched assembly.
    for method in ('G', 'i'):
        assert np.array_equal(
            np.asarray(getattr(batched, method)(x, defaultepar, skip=frozenset())),
            np.asarray(getattr(batched, method)(x, defaultepar)))


## ---------------------------------------------------------------------------
## STAGE 14-2 -- a BJT under PCNR, in DC.
##
## `Semiconductor` derives `pcnr_i`/`pcnr_didv` from its own `eval_i_pure`
## through `x = M v`, so the model still lives once.  Asserted on the MATRIX as
## well as on answers (design study section 9): answers within Newton tolerance
## cannot see a Jacobian that is wrong but convergent.
## ---------------------------------------------------------------------------

_BJT_P = {'IS': 1e-14, 'BF': 100.0, 'BR': 1.0, 'VA': 100.0}


def _ebers_moll_didv(vbe, vbc):
    """Analytic d(ic, ib, ie)/d(vbe, vbc) of BJT.eval_i_pure, Early factor included."""
    VT = numeric.kboltzmann * defaultepar.T / numeric.qelectron
    IS, BF, BR, VA = (_BJT_P[k] for k in ('IS', 'BF', 'BR', 'VA'))
    ebe, ebc = np.exp(vbe / VT), np.exp(vbc / VT)
    i_be, i_bc = IS * (ebe - 1.0), IS * (ebc - 1.0)
    g_be, g_bc = IS * ebe / VT, IS * ebc / VT
    early = 1.0 + (vbe - vbc) / VA
    dbe = g_be * early + (i_be - i_bc) / VA
    dbc = -g_bc * early - (i_be - i_bc) / VA
    return np.array([[dbe, dbc - g_bc / BR],
                     [g_be / BF, g_bc / BR],
                     [-dbe - g_be / BF, -dbc]])


_BIASES = [(0.70, -4.30), (0.75, 0.60), (-1.0, -5.0), (-3.0, 0.70)]


@pytest.mark.parametrize('vbe,vbc', _BIASES,
                         ids=['forward', 'saturated', 'cutoff', 'reverse'])
def test_gate_14_2_bjt_pcnr_didv_is_the_device_jacobian_through_the_map(vbe, vbc):
    """Exactly the toolkit's own Jacobian chained through M -- the adapter adds
    nothing -- and that Jacobian is Ebers-Moll's to central-difference accuracy.
    In cut-off every conductance is ~6e-30 S and the difference quotient returns
    0 (`IS*(exp(-38.7) - 1)` cannot resolve the exponential), so there the check
    is absolute."""
    from pycircuit.circuit.semiconductors import BJT
    v = np.array([vbe, vbc])
    M = BJT._pcnr_map()
    D = np.asarray(BJT.pcnr_didv(v, _BJT_P, defaultepar, numeric), dtype=float)
    GM = np.asarray(numeric.jacobian(BJT.eval_i_pure, M @ v, _BJT_P, defaultepar),
                    dtype=float) @ M
    assert np.array_equal(D, GM)
    A = _ebers_moll_didv(vbe, vbc)
    if np.max(np.abs(A)) > 1e-20:
        assert np.max(np.abs(D - A)) / np.max(np.abs(A)) < 1e-6
    else:
        assert np.max(np.abs(D - A)) < 1e-20


@pytest.mark.parametrize('vbe,vbc', _BIASES,
                         ids=['forward', 'saturated', 'cutoff', 'reverse'])
def test_gate_14_2_bjt_pcnr_didv_is_ebers_moll_under_autodiff(vbe, vbc):
    """With exact differentiation the chain rule itself is checked, with no
    difference quotient in the way: measured <= 1.1e-14 at all four biases."""
    pytest.importorskip('jax')
    from pycircuit.circuit.semiconductors import BJT
    from pycircuit.circuit.toolkit import jaxtoolkit
    D = np.asarray(BJT.pcnr_didv(np.array([vbe, vbc]), _BJT_P, defaultepar,
                                 jaxtoolkit), dtype=float)
    A = _ebers_moll_didv(vbe, vbc)
    assert np.max(np.abs(D - A)) / np.max(np.abs(A)) < 1e-12


def test_gate_14_2_pcnr_i_reproduces_the_model_at_any_common_offset():
    from pycircuit.circuit.semiconductors import BJT
    for vbe, vbc in _BIASES:
        xb = 1.234
        x = np.array([xb - vbc, xb, xb - vbe])
        ref = np.asarray(BJT.eval_i_pure(x, _BJT_P, defaultepar, numeric), dtype=float)
        got = np.asarray(BJT.pcnr_i(np.array([vbe, vbc]), _BJT_P, defaultepar,
                                    numeric), dtype=float)
        np.testing.assert_allclose(got, ref, rtol=1e-13, atol=1e-30)


def test_gate_14_2_a_non_spanning_declaration_is_refused_at_discovery():
    from pycircuit.circuit.semiconductors import BJT

    class TwoOnOnePair(BJT):
        pcnr_junctions = ((1, 2), (2, 1))     # the same pair twice: c is unreached

    c = SubCircuit()
    c['Q'] = TwoOnOnePair('c', 'b', gnd)
    with pytest.raises(ValueError, match='spanning tree'):
        P.pcnr_junctions(c)


def _ce(rb):
    from pycircuit.circuit.semiconductors import BJT
    cir = SubCircuit(toolkit=numeric)
    cir['VCC'] = VS('vcc', gnd, v=5.0)
    cir['VIN'] = VS('vin', gnd, v=0.8)
    cir['RB'] = R('vin', 'b', r=rb)
    cir['RC'] = R('vcc', 'c', r=1e3)
    cir['Q1'] = BJT('c', 'b', gnd)
    return cir


@pytest.mark.parametrize('rb', [10.0, 100.0, 1000.0])
def test_gate_14_2_pcnr_solves_the_common_emitter_stage(rb):
    """Entry: `P.solve_dc` refused this circuit outright.  It must now land on
    classic's operating point (measured <= 1.5e-12, in the same iteration
    counts), satisfy the unlimited equations, and never touch the device's
    classic limiter or its node-voltage stamps."""
    import warnings
    from pycircuit.circuit.semiconductors import BJT
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ref = np.asarray(DC(_ce(rb), toolkit=numeric).solve().x, dtype=float).ravel()
    c = _ce(rb)
    calls = {'limit': 0, 'G': 0, 'i': 0}
    real = {k: getattr(BJT, k) for k in calls}

    def counted(name):
        def spy(*a, **k):
            calls[name] += 1
            return real[name](*a, **k)
        return spy
    for k in calls:
        setattr(BJT, k, counted(k))
    try:
        x, _v, _its = P.solve_dc(c, gnd)
    finally:
        for k in calls:
            setattr(BJT, k, real[k])
    assert calls == {'limit': 0, 'G': 0, 'i': 0}, calls
    x = np.asarray(x, dtype=float)
    np.testing.assert_allclose(x, ref, rtol=1e-9, atol=1e-12)
    f = (np.asarray(c.i(x, defaultepar), dtype=float)
         + np.asarray(c.u(0, analysis='dc', epar=defaultepar), dtype=float))
    assert np.max(np.abs(np.delete(f, c.get_node_index(gnd)))) < 1e-9
    assert 0.4 < float(x[c.get_node_index('b')]) < 0.9


def test_gate_14_2_dc_with_pcnr_runs_pcnr_on_a_bjt_circuit():
    """Entry: `DC(pcnr=True)` silently ran classic limiting here -- 0 PCNR
    assemblies against 18 limiter calls."""
    import warnings
    real = P.augmented_system
    hits = {'n': 0}

    def spy(*a, **k):
        hits['n'] += 1
        return real(*a, **k)
    P.augmented_system = spy
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            DC(_ce(100.0), toolkit=numeric, pcnr=True).solve()
    finally:
        P.augmented_system = real
    assert hits['n'] > 5


## ---------------------------------------------------------------------------
## STAGE 14-3 -- charge on the limited unknowns: a BJT transient under PCNR.
##
## The charge is evaluated where the current is, at the device's own unknowns:
## `q`/`C` are assembled without it, its `q(v_lim)` is added back before
## `get_diff`, and `alpha0 * dq/dv` joins its J_ml columns.  Entry: this
## transient silently ran classic limiting (then, after 14-2, refused).
## ---------------------------------------------------------------------------

def _bjt_switch():
    from pycircuit.circuit.elements import VPulse
    from pycircuit.circuit.semiconductors import BJT
    per = 1e-6
    cir = SubCircuit(toolkit=numeric)
    cir['VCC'] = VS('vcc', gnd, v=5.0)
    cir['VB'] = VPulse('vb', gnd, v1=0.0, v2=0.9, td=per / 10, tr=per / 1000,
                       tf=per / 1000, pw=per / 2, per=per)
    cir['RB'] = R('vb', 'b', r=1e3)
    cir['RC'] = R('vcc', 'c', r=1e3)
    cir['Q1'] = BJT('c', 'b', gnd)
    return cir


def _bjt_run(pcnr, coupled=False):
    import warnings
    from pycircuit.circuit.transient import Transient
    tran = Transient(_bjt_switch(), toolkit=numeric, reltol=1e-5, pcnr=pcnr)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=1.6e-6, timestep=2.5e-9, coupled_lte=coupled)
    w = res.v('c')
    return (tran, np.asarray(w.x, dtype=float).ravel(),
            np.asarray(w.y, dtype=float).ravel())


def test_gate_14_3_pcnr_and_limiting_take_the_same_steps_on_a_bjt_switch():
    """13-6's standard: the route to a converged solution must not change the
    solution, so equal solutions give equal LTE and equal steps.  Measured: the
    same 1136 accepted and 36 rejected steps, time points within 4.7e-10
    relative, waveform within 1.1e-10."""
    tc, t_c, y_c = _bjt_run(False)
    tp, t_p, y_p = _bjt_run(True)
    assert tp.statistics.accepted_steps == tc.statistics.accepted_steps
    assert tp.statistics.rejected_steps == tc.statistics.rejected_steps
    np.testing.assert_allclose(t_p, t_c, rtol=1e-8, atol=0)
    assert np.max(np.abs(y_p - y_c)) / np.max(np.abs(y_c)) < 1e-6


def test_gate_14_3_the_bjt_is_never_evaluated_at_node_voltages_under_pcnr():
    """From the first PCNR timestep on (the DC operating point still limits
    classically, as recorded at 13-7): no limiter, no node-voltage stamp."""
    from pycircuit.circuit.semiconductors import BJT
    from pycircuit.circuit.transient import Transient
    names = ('limit', 'G', 'i', 'C', 'q')
    calls = dict.fromkeys(names, 0)
    real = {k: getattr(BJT, k) for k in names}
    real_step = Transient._solve_timestep_pcnr
    armed = {'on': False}

    def step(self, *a, **k):
        armed['on'] = True
        return real_step(self, *a, **k)

    def counted(name):
        def spy(*a, **k):
            if armed['on']:
                calls[name] += 1
            return real[name](*a, **k)
        return spy
    Transient._solve_timestep_pcnr = step
    for k in names:
        setattr(BJT, k, counted(k))
    try:
        _bjt_run(True)
    finally:
        Transient._solve_timestep_pcnr = real_step
        for k in names:
            setattr(BJT, k, real[k])
    assert armed['on'], 'the PCNR timestep path never ran'
    assert calls == dict.fromkeys(names, 0), calls


def test_gate_14_3_controller_jacobian_is_the_residual_jacobian():
    """Assert on the MATRIX (design study section 9): at converged time points
    the Jacobian handed to the step controller must equal a central difference
    of the full residual i(x) + iq(x) + u.  Measured worst 7.5e-10 over 40
    points; a dropped `alpha0 * dq/dv` column is invisible to every waveform
    test that agrees within truncation error, but not to this."""
    from pycircuit.circuit.transient import Transient
    real = Transient._solve_timestep_pcnr
    worst = {'rel': 0.0, 'n': 0, 'seen': 0}

    def spy(self, x0, t, *a, **k):
        x, feval, J, f = real(self, x0, t, *a, **k)
        worst['seen'] += 1
        if worst['seen'] % 97 or worst['n'] >= 6:
            return x, feval, J, f
        cir, ep = self.cir, self.epar
        h_last = getattr(self, '_dt_last', self._dt)

        def F(xx):
            q = np.asarray(cir.q(xx, ep), dtype=float)
            iq, _g = self.active_integrator.compute_derivatives(
                q_curr=q, C_curr=1.0, h_curr=self._dt, q_last=self._qlast,
                iq_last=self._iqlast, h_last=h_last,
                is_first_step=getattr(self, '_is_first_step', False),
                toolkit=self.toolkit)
            return (np.asarray(cir.i(xx, ep), dtype=float) + np.asarray(iq, dtype=float)
                    + np.asarray(cir.u(t, ep, analysis=self.par.analysis), dtype=float))
        xx = np.asarray(x, dtype=float)
        n = len(xx)
        FD = np.zeros((n, n))
        for j in range(n):
            dl = 1e-7 * max(1.0, abs(xx[j]))
            xp, xm = xx.copy(), xx.copy()
            xp[j] += dl
            xm[j] -= dl
            FD[:, j] = (F(xp) - F(xm)) / (2 * dl)
        keep = [i for i in range(n) if i != self.irefnode]
        Jr = np.asarray(J, dtype=float)[np.ix_(keep, keep)]
        FDr = FD[np.ix_(keep, keep)]
        worst['rel'] = max(worst['rel'], float(np.max(np.abs(Jr - FDr)) / np.max(np.abs(FDr))))
        worst['n'] += 1
        return x, feval, J, f
    Transient._solve_timestep_pcnr = spy
    try:
        _bjt_run(True)
    finally:
        Transient._solve_timestep_pcnr = real
    assert worst['n'] >= 6
    assert worst['rel'] < 1e-6, worst


def test_gate_14_3_pcnr_runs_on_the_coupled_path_with_a_bjt():
    from pycircuit.circuit.semiconductors import BJT
    hits = {'aug': 0, 'limit': 0}
    ra, rl = P.augmented_system, BJT.limit

    def sa(*a, **k):
        hits['aug'] += 1
        return ra(*a, **k)

    def sl(*a, **k):
        hits['limit'] += 1
        return rl(*a, **k)
    P.augmented_system, BJT.limit = sa, sl
    try:
        _tc, t_p, y_p = _bjt_run(True, coupled=True)
    finally:
        P.augmented_system, BJT.limit = ra, rl
    assert hits['aug'] > 100
    assert hits['limit'] < 10, 'the coupled loop is limiting classically'
    _tc2, t_c, y_c = _bjt_run(False, coupled=True)
    assert t_p[-1] == t_c[-1] == 1.6e-6
    assert abs(y_p[-1] - y_c[-1]) < 1e-6


def test_gate_14_3_a_caller_that_does_not_move_the_charge_is_refused():
    """`augmented_system` with `J_extra` but no `charge_scale` on a charge-
    storing PCNR device would leave its charge at the node voltages."""
    c = _bjt_switch()
    js = P.pcnr_junctions(c)
    x = np.zeros(c.n)
    v = np.zeros(len(js))
    with pytest.raises(NotImplementedError, match='caller did not'):
        P.augmented_system(c, x, v, js, defaultepar, J_extra=np.zeros((c.n, c.n)))


def test_the_pcnr_timestep_path_counts_its_newton_iterations():
    """It counted none, so a PCNR run reported 0 Newton iterations.  Measured
    after the fix: 9647 on the 5 V rectifier against classic's 9654."""
    tran = Transient(_rectifier(), toolkit=numeric, reltol=1e-5, pcnr=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tran.solve(tend=2e-3, timestep=1e-5)
    s = tran.statistics
    assert s.newton_iterations >= s.accepted_steps > 0
