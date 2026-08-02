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
    g_mna, g_lim, J_mm, J_ml, J_lm = P.augmented_system(c, x, v_lim, js)

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
    _g, _gl, J_mm, J_ml, _J_lm = P.augmented_system(c, x, v_lim, js)

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
