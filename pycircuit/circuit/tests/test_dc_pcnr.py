import pytest

from pycircuit.circuit.tests._warnpolicy import quiet
def test_dc_pcnr_diode():
    """
    Test PCNR limiting in DC simulation with a diode driven by a large source.
    Without PCNR pnjlim algorithm, solving a diode often causes Python OverflowErrors 
    during early Newton-Raphson iterations when predictions overshoot massively.
    With PCNR, the solver clamps the exponential argument gracefully.
    """
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import Diode, IS
    from pycircuit.circuit.dcanalysis import DC

    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1.0) # 1A current
    c['D1'] = Diode(1, gnd)
    
    dc = DC(c)
    # The solver should converge without OverflowError due to PCNR pnjlim
    # The limiter limits x_next during StandardNewton iteration
    res = dc.solve()
    
    # Check that diode forward voltage settles properly (approx 0.7 - 0.9V)
    v_diode = res.v(1)
    assert 0.5 < v_diode < 1.5, f"Diode voltage {v_diode} outside expected forward bias range"


def test_dc_pcnr_true_takes_a_circuit_with_no_junction_at_all():
    """`DC(pcnr=True)` used to gate PARTICIPATION on `pcnr_junctions()`,
    the pnj-only pair view built for the gmin ladders.  A circuit whose
    only participants carry `fetlim`/`limvds` probes has an EMPTY pair
    view, so a MOSFET differential pair asked for PCNR silently got the
    ordinary solver instead -- and the vector-PCNR acceptance tables
    were taken through `pcnr.solve_dc` directly for exactly that reason
    (Stage 2, 2026-08-26).  Gated on `pcnr_devices()` now, and this test
    would have caught the fall-through: it counts the calls.
    """
    import numpy as np
    from pycircuit.circuit import gnd, pcnr
    from pycircuit.circuit.elements import SubCircuit, VS, R, IS
    from pycircuit.circuit.dcanalysis import DC
    import pycircuit.circuit.tests.test_limit_fet as T
    c = SubCircuit()
    for n in ('vdd', 'inp', 'inn', 'o1', 'o2', 'tail'):
        c.add_node(n)
    c['vdd'] = VS('vdd', gnd, v=5.0)
    c['vp'] = VS('inp', gnd, v=3.5)
    c['vn'] = VS('inn', gnd, v=1.5)
    c['r1'] = R('vdd', 'o1', r=5e3)
    c['r2'] = R('vdd', 'o2', r=5e3)
    c['it'] = IS('tail', gnd, i=200e-6)
    F = T._fet('both')
    c['M1'] = F('o1', 'inp', 'tail')
    c['M2'] = F('o2', 'inn', 'tail')
    ## The precondition that made the bug invisible: no pnj pair at all.
    assert len(pcnr.pcnr_junctions(c)) == 0
    assert len(pcnr.pcnr_devices(c)) == 2
    calls = [0]
    orig = pcnr.solve_dc

    def spy(*a, **k):
        calls[0] += 1
        return orig(*a, **k)

    pcnr.solve_dc = spy
    try:
        with quiet():
            r = DC(c, pcnr=True).solve()
    finally:
        pcnr.solve_dc = orig
    assert calls[0] >= 1, 'pcnr=True fell through to the ordinary solver'
    ## 2.818593 with the exact Boltzmann constant (2.818678 with 1.38e-23,
    ## re-pinned 2026-09-05 -- the thermal voltage moved by 4.7e-4)
    assert abs(float(r.v('tail', gnd)) - 2.818523) < 1e-5


def _bjt_mirror_20v(rev=False):
    """The Stage 0/1 BJT mirror: from a uniform 20 V start PCNR's undamped
    first step proposes ~9e45 V and pnjlim's law cannot recover it, in
    either instance order; the ordinary chain solves it in ~146."""
    import numpy as np
    from pycircuit.circuit import gnd, elements_hdl as eh
    from pycircuit.circuit.elements import SubCircuit, VS, R
    import pycircuit.circuit.tests.test_elements_hdl_library3 as T
    c = SubCircuit()
    vcc = c.add_node('vcc'); nb = c.add_node('nb'); no = c.add_node('no')
    c['vcc'] = VS(vcc, gnd, v=5.0)
    c['rref'] = R(vcc, nb, r=4.3e3)
    c['rl'] = R(vcc, no, r=1e3)
    devs = [('q1', lambda: eh.GummelPoonNpnHdl(nb, nb, gnd, **T.NPN)),
            ('q2', lambda: eh.GummelPoonNpnHdl(no, nb, gnd, **T.NPN))]
    for n, f in (devs[::-1] if rev else devs):
        c[n] = f()
    x0 = np.full(c.n, 20.0)
    x0[c.get_node_index(gnd)] = 0.0
    return c, x0


def test_dc_pcnr_true_falls_back_to_the_ordinary_chain_when_pcnr_fails():
    """`DC(pcnr=True)` used to return straight from `pcnr.solve_dc`, so a
    PCNR failure was a HARD failure on circuits the default path solves.
    So: fall back, flag it, keep the answer.

    ⚠⚠ **THE FAILURE IS NOW FORCED, because nothing fails any more.**

    This test has run on three different circuits.  It began on the BJT
    mirror, which `v_lim_init`'s limited seed fixed (sec. 27); it moved
    to the EKV differential pair, which `limit_delta` on the bulk fixed
    (sec. 34).  Probing 72 configurations afterwards -- four circuit
    families, starts from -40 V to +100 V, both instance orders --
    **`pcnr.solve_dc` converged on every one.**

    A precondition that stops holding does not make a test pass; it
    makes it VACUOUS, and both times this one failed loudly with `DID
    NOT RAISE` rather than passing quietly.  Having run out of circuits,
    it now raises from `solve_dc` deliberately.  That is a weaker test of
    *when* the fallback fires and a better test of *what it does* -- and
    it cannot rot again.

    The fallback stays because "nothing I can find fails" is not "nothing
    fails".
    """
    import logging
    import numpy as np
    from pycircuit.circuit import gnd, pcnr
    from pycircuit.circuit.dcanalysis import DC

    for rev in (False, True):
        c, x0 = _bjt_mirror_20v(rev)
        ## the precondition, forced: PCNR raises the way it used to
        real = pcnr.solve_dc

        def boom(*a, **k):
            raise RuntimeError('forced PCNR failure')

        pcnr.solve_dc = boom
        logging.disable(logging.CRITICAL)
        try:
            a = DC(c, pcnr=True)
            r = a.solve(x0=x0)
        finally:
            pcnr.solve_dc = real
            logging.disable(logging.NOTSET)
        ref = DC(_bjt_mirror_20v(rev)[0]).solve(x0=x0)
        assert a.pcnr_fell_back is True
        ## and the richer form agrees -- 'fell-back' is distinguishable
        ## from 'no-participants', which the boolean cannot do.
        assert a.pcnr_status == 'fell-back'
        assert abs(float(r.v('no', gnd)) - float(ref.v('no', gnd))) < 1e-6


def test_pcnr_now_converges_where_it_used_to_need_the_fallback():
    """The other half, and the one that would catch a regression.

    Both circuits this module's fallback test used to run on now go
    through PCNR itself, from the wild starts that used to break them.
    """
    import numpy as np
    from pycircuit.circuit import gnd, pcnr
    for rev in (False, True):
        c, x0 = _bjt_mirror_20v(rev)
        x, v, its = pcnr.solve_dc(c, gnd, x0=x0)
        ## 20 with the exact Boltzmann constant (19 with 1.38e-23, 2026-09-05);
        ## the budget is 200 and the claim is convergence, not the count
        assert its < 25
        assert abs(float(x[c.get_node_index('no')]) - 3.9956) < 5e-3

    from pycircuit.circuit.tests.test_limit_identity import _diffpair
    for rev in (False, True):
        c = _diffpair(0.0, rev)
        x0 = np.full(c.n, 20.0)
        x0[c.get_node_index(gnd)] = 0.0
        x, v, its = pcnr.solve_dc(c, gnd, x0=x0)
        ## 20 with the exact Boltzmann constant (19 with 1.38e-23, 2026-09-05);
        ## the budget is 200 and the claim is convergence, not the count
        assert its < 25


def test_dc_pcnr_true_does_not_fall_back_when_pcnr_converges():
    """The flag is False on a solve PCNR handled itself -- otherwise the
    test above could pass on a fallback that always fires."""
    import numpy as np
    from pycircuit.circuit import gnd
    from pycircuit.circuit.dcanalysis import DC
    c, x0 = _bjt_mirror_20v(False)
    a = DC(c, pcnr=True)
    r = a.solve()                     # zero start: PCNR converges in ~9
    assert a.pcnr_fell_back is False
    assert abs(float(r.v('no', gnd)) - 3.996) < 5e-3


## ======================================================================
## `pcnr=True` is a REQUEST -- what actually happened is `pcnr_status`.
## ======================================================================

def _linear_only():
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit import elements_hdl as eh
    c = SubCircuit()
    c.add_node('a')
    c['v1'] = VS('a', gnd, v=1.0)
    c['r1'] = eh.RHdl('a', gnd, r=1e3)
    return c


def _with_diode():
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit import elements_hdl as eh
    c = SubCircuit()
    c.add_node('a')
    c['v1'] = VS('a', gnd, v=1.0)
    c['d1'] = eh.DiodeSpiceHdl('a', gnd)
    return c


def test_pcnr_status_is_always_defined_even_when_pcnr_is_off():
    """It was assigned only inside `if self.par.pcnr`.

    So a caller that checked `dc.pcnr_fell_back` after an ordinary solve
    got an `AttributeError` rather than an answer -- the attribute
    existed only when it was least interesting.
    """
    from pycircuit.circuit.circuit import gnd
    from pycircuit.circuit.dcanalysis import DC

    dc = DC(_with_diode(), refnode=gnd)          # pcnr defaults to False
    dc.solve()
    assert dc.pcnr_status == 'off'
    assert dc.pcnr_fell_back is False


def test_pcnr_status_separates_ran_from_did_nothing():
    """⚠ THE DEFECT THIS EXISTS FOR.

    `pcnr_fell_back` reads **False** both when PCNR solved the circuit
    and when `pcnr=True` was asked for and no device declared a probe.
    A caller could not tell "PCNR solved this" from "PCNR did nothing",
    which is the silent fallthrough this module's own comment records
    fixing for the MOSFET diff pair -- the participation GATE was fixed
    and the REPORTING was not.
    """
    from pycircuit.circuit.circuit import gnd
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit import pcnr as P

    ran = DC(_with_diode(), refnode=gnd, pcnr=True)
    ran.solve()
    assert P.pcnr_devices(_with_diode()), 'the fixture stopped participating'
    assert ran.pcnr_status == 'used'

    nothing = DC(_linear_only(), refnode=gnd, pcnr=True)
    nothing.solve()
    assert not P.pcnr_devices(_linear_only())
    assert nothing.pcnr_status == 'no-participants'

    ## The old boolean cannot tell them apart, which is why the string
    ## exists.  Asserted so a future change that tried to fold them back
    ## into one flag has to notice.
    assert ran.pcnr_fell_back is nothing.pcnr_fell_back is False


def test_a_circuit_with_no_limited_device_still_solves_with_pcnr_on():
    """`pcnr=True` must never make a solvable circuit unsolvable.

    `pcnr.solve_dc` raises "no device in this circuit declares a PCNR
    junction" when called directly; through `DC(pcnr=True)` that case
    falls through to the ordinary chain and returns the answer.
    """
    import numpy as np
    from pycircuit.circuit.circuit import gnd
    from pycircuit.circuit.dcanalysis import DC

    off = np.asarray(DC(_linear_only(), refnode=gnd).solve().x, dtype=float)
    on = DC(_linear_only(), refnode=gnd, pcnr=True)
    got = np.asarray(on.solve().x, dtype=float)
    assert on.pcnr_status == 'no-participants'
    assert np.allclose(got, off, rtol=0, atol=0), (got, off)


## ======================================================================
## The transient counterpart: a transient cannot have ONE status.
## ======================================================================

def _tran_diode():
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS, R
    from pycircuit.circuit import elements_hdl as eh
    c = SubCircuit()
    c.add_node('a')
    c['v1'] = VS('a', gnd, v=1.0)
    c['d1'] = eh.DiodeSpiceHdl('a', gnd)
    c['r1'] = R('a', gnd, r=1e3)
    return c


def test_a_transient_counts_the_steps_pcnr_carried():
    """DC reports one `pcnr_status`; a transient cannot.

    PCNR is attempted per timestep and can fall through on any of them,
    so the honest report is a COUNT.  `transient.py` had the same
    fallback DC did and recorded nothing at all -- a run could lose PCNR
    on some steps and leave no trace but log lines.

    The counts are SOLVER INVOCATIONS, not accepted steps: a rejected
    step is solved and thrown away, and it is still a step PCNR did or
    did not carry.
    """
    from pycircuit.circuit.circuit import gnd
    from pycircuit.circuit.transient import Transient

    off = Transient(_tran_diode(), pcnr=False)
    off.solve(refnode=gnd, tend=2e-9, timestep=2e-10)
    assert off.pcnr_status == 'off'
    assert off.pcnr_solves == 0 and off.pcnr_fallbacks == 0

    on = Transient(_tran_diode(), pcnr=True)
    on.solve(refnode=gnd, tend=2e-9, timestep=2e-10)
    assert on.pcnr_status == 'used'
    assert on.pcnr_solves > 0, 'pcnr=True carried no step at all'
    assert on.pcnr_fallbacks == 0


def test_a_transient_says_when_nothing_participates():
    """`pcnr=True` on an all-linear circuit is indistinguishable from
    `pcnr=False` unless it says so.  Same defect DC carried."""
    from pycircuit.circuit.circuit import SubCircuit, gnd
    from pycircuit.circuit.elements import VS
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.transient import Transient

    c = SubCircuit()
    c.add_node('a')
    c['v1'] = VS('a', gnd, v=1.0)
    c['r1'] = eh.RHdl('a', gnd, r=1e3)
    t = Transient(c, pcnr=True)
    t.solve(refnode=gnd, tend=2e-9, timestep=2e-10)
    assert t.pcnr_status == 'no-participants'
    assert t.pcnr_solves == 0


@pytest.mark.filterwarnings('ignore::pycircuit.circuit.simwarnings.ConvergenceWarning')
def test_a_partly_failing_transient_reports_partial():
    """⚠ The case the counter exists for.

    A transient that loses PCNR on SOME steps is the interesting one --
    it is neither 'used' nor 'fell-back', and before this it looked
    exactly like a clean PCNR run from outside.  Injected by making the
    PCNR step solver fail on every third call.
    """
    import logging
    from pycircuit.circuit.circuit import gnd
    from pycircuit.circuit.transient import Transient

    t = Transient(_tran_diode(), pcnr=True)
    real = t._solve_timestep_pcnr
    calls = {'n': 0}

    def flaky(*a, **kw):
        calls['n'] += 1
        if calls['n'] % 3 == 0:
            raise RuntimeError('injected PCNR failure')
        return real(*a, **kw)

    t._solve_timestep_pcnr = flaky
    logging.disable(logging.CRITICAL)
    try:
        t.solve(refnode=gnd, tend=2e-9, timestep=2e-10)
    finally:
        logging.disable(logging.NOTSET)

    assert t.pcnr_fallbacks > 0, 'the injection did not fire'
    assert t.pcnr_solves > 0, 'nothing was carried by PCNR either'
    assert t.pcnr_status == 'partial', (
        t.pcnr_status, t.pcnr_solves, t.pcnr_fallbacks)


def _failing(which):
    """`pcnr.solve_dc` raising for the calls `which(kw)` names -- the plain
    solve, a gshunt rung -- and running every other; returns the patched
    function and the calls made, as `(gshunt, source)` pairs."""
    from pycircuit.circuit import pcnr
    real, calls = pcnr.solve_dc, []

    def solve_dc(*a, **k):
        calls.append((k.get('gshunt', 0.0), k.get('source', 1.0)))
        if which(k):
            raise RuntimeError('forced PCNR failure')
        return real(*a, **k)
    return solve_dc, calls


def _plain_call(k):
    return (k.get('gshunt', 0.0) == 0.0 and k.get('source', 1.0) == 1.0
            and k.get('v_lim0') is None)


@pytest.mark.parametrize('ladder', ['gshunt', 'source'])
def test_pcnr_continues_on_its_own_ladders_where_the_plain_solve_fails(ladder, monkeypatch):
    """`pcnr.solve_dc_continued`: the plain PCNR solve forced to fail, the
    gshunt ladder carries it -- and with every gshunt rung failing too,
    source stepping does; PCNR solves it either way (`pcnr_status`
    'used', `pcnr_path` the ladder), the answer the ordinary chain's."""
    from pycircuit.circuit import gnd, pcnr
    from pycircuit.circuit.dcanalysis import DC
    if ladder == 'gshunt':
        fails = _plain_call
    else:
        def fails(k):
            return _plain_call(k) or k.get('gshunt', 0.0) > 0.0
    patched, calls = _failing(fails)
    monkeypatch.setattr(pcnr, 'solve_dc', patched)
    c, _x0 = _bjt_mirror_20v()
    a = DC(c, pcnr=True)
    r = a.solve()
    ref = DC(_bjt_mirror_20v()[0]).solve()
    assert a.pcnr_status == 'used' and a.pcnr_path == ladder, (a.pcnr_status, a.pcnr_path)
    assert any(g > 0.0 for g, _s in calls) and (ladder == 'gshunt'
                                                 or any(s < 1.0 for _g, s in calls))
    assert abs(float(r.v('no', gnd)) - float(ref.v('no', gnd))) < 1e-6


def _failing_kinds(which):
    """`_failing`, each call recorded as `(kind, value, source)` -- kind 'p'
    a plain or pure solve, 'g' a gshunt rung, 'j' a junction-gmin rung."""
    from pycircuit.circuit import pcnr
    real, calls = pcnr.solve_dc, []

    def solve_dc(*a, **k):
        g, j = k.get('gshunt', 0.0), k.get('jgmin', 0.0)
        calls.append(('g' if g else ('j' if j else 'p'), g or j, k.get('source', 1.0)))
        if which(k):
            raise RuntimeError('forced PCNR failure')
        return real(*a, **k)
    return solve_dc, calls


def test_a_ladder_that_came_down_is_followed_by_the_junction_ladder(monkeypatch):
    """Every solve at full sources without a deformation forced to fail --
    the plain solve, each ladder's pure solve -- every rung converging: the
    gshunt ladder runs its eleven decades once (1e-2 down to the last at
    ~1e-12, not one more rung at exactly GSHUNT_END: ten divisions by 10
    leave 1.0000000000000002e-12), its pure solve fails, then the junction-
    gmin ladder runs its ten (1e-2 to 1e-11) and its own pure solve fails,
    then source stepping -- which never reaches full sources here -- and the
    plain solve's error is raised.  (Until the junction-gmin ladder, a gshunt
    ladder that came down ended the route: only the unshunted system had
    failed, which source stepping at no shunt could not change; the
    junction ladder deforms the junctions instead, and with it the route
    is the ordinary chain's.)"""
    from pycircuit.circuit import gnd, pcnr
    patched, calls = _failing_kinds(lambda k: k.get('gshunt', 0.0) == 0.0 and k.get('jgmin', 0.0) == 0.0
                                    and k.get('source', 1.0) == 1.0)
    monkeypatch.setattr(pcnr, 'solve_dc', patched)
    c, x0 = _bjt_mirror_20v()
    with pytest.raises(RuntimeError, match='forced PCNR failure'):
        pcnr.solve_dc_continued(c, gnd, x0=x0)
    head = calls[:24]
    assert [k for k, _v, _s in head] == ['p'] + ['g'] * 11 + ['p'] + ['j'] * 10 + ['p'], head
    assert all(s == 1.0 for _k, _v, s in head), head
    g = [v for k, v, _s in head if k == 'g']
    assert g[0] == pcnr.GSHUNT_START and abs(g[-1] / pcnr.GSHUNT_END - 1.0) < 1e-9, g
    jg = [v for k, v, _s in head if k == 'j']
    assert jg[0] == pcnr.JGMIN_START and abs(jg[-1] / pcnr.JGMIN_END - 1.0) < 1e-9, jg
    assert calls[24][2] < 1.0, calls[24]


def test_source_stepping_gives_up_after_four_failing_rungs(monkeypatch):
    """Every PCNR solve forced to fail: the plain solve, the gshunt ladder's
    first rung, the junction-gmin ladder's, then source stepping's four steps
    from the start (0.1 halved to SOURCE_STEP_MIN), each tried plain and by
    both ladders' first rungs at its level, and no more."""
    from pycircuit.circuit import gnd, pcnr
    patched, calls = _failing_kinds(lambda k: True)
    monkeypatch.setattr(pcnr, 'solve_dc', patched)
    c, x0 = _bjt_mirror_20v()
    with pytest.raises(RuntimeError, match='forced PCNR failure'):
        pcnr.solve_dc_continued(c, gnd, x0=x0)
    step = [('p', 0.0), ('g', pcnr.GSHUNT_START), ('j', pcnr.JGMIN_START)]
    assert calls == [(k, v, s) for s in (1.0, 0.1, 0.05, 0.025, 0.0125) for k, v in step], calls


def test_a_gold_deck_pcnr_solves_itself():
    """slowlatch, which used to fall back to the ordinary chain, is solved
    by PCNR itself, and the point satisfies KCL.  Which path is a knife
    edge at rounding level: its gshunt ladder at this commit, its plain
    solve when the view's stubs sat after every other element (the
    ladders themselves are pinned by the forced-failure tests above).
    Skips where the benchmark data is not fetched."""
    import numpy as np

    from pycircuit._testing import benchdata
    from pycircuit.circuit import circuit, spice_import
    from pycircuit.circuit.dcanalysis import DC
    text = benchdata.spice_data('Netlists/CircuitSim90/MOS2/slowlatch.cir')
    if text is None:
        pytest.skip('benchmark data not fetched (benchmarks/fetch_spice_suite.py)')
    tr, _kw = spice_import.import_netlist(text).transient()
    cir = tr.cir
    a = DC(cir, toolkit=circuit.numeric, epar=tr.epar, pcnr=True)
    with quiet():
        r = a.solve()
    assert a.pcnr_status == 'used' and a.pcnr_path in ('plain', 'gshunt', 'source'), (
        a.pcnr_status, a.pcnr_path)
    x = np.asarray(r.x, dtype=float)
    f = np.asarray(cir.i(x, tr.epar), dtype=float) + np.asarray(
        cir.u(0, analysis='dc', epar=tr.epar), dtype=float)
    keep = np.arange(cir.n) != a.irefnode
    assert np.abs(f[keep]).max() < 1e-9

