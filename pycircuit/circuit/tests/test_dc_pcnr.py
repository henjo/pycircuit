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
