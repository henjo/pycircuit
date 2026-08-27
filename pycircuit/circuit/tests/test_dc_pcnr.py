import pytest
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
    import warnings
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
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            r = DC(c, pcnr=True).solve()
    finally:
        pcnr.solve_dc = orig
    assert calls[0] >= 1, 'pcnr=True fell through to the ordinary solver'
    assert abs(float(r.v('tail', gnd)) - 2.818678) < 1e-5


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
        assert its < 20
        assert abs(float(x[c.get_node_index('no')]) - 3.9956) < 5e-3

    from pycircuit.circuit.tests.test_limit_identity import _diffpair
    for rev in (False, True):
        c = _diffpair(0.0, rev)
        x0 = np.full(c.n, 20.0)
        x0[c.get_node_index(gnd)] = 0.0
        x, v, its = pcnr.solve_dc(c, gnd, x0=x0)
        assert its < 20


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
