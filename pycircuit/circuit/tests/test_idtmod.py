# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Phase-1 idtmod semantics (idtmod.md sec. 5.1).

The LRM contract under test: ``idtmod`` returns ``k`` with
``offset <= k < offset + modulus`` and ``integral + ic = n*modulus + k`` --
a *floored* wrap, which the old ``(-s % m) + offset`` violated for every
``offset`` that is not a multiple of ``modulus``.  Plus the two element
mechanics that hang off it: ``ic`` pinning the DC solution (LRM: "the
initial condition shall force the DC solution"), and ``i()`` being pure
array arithmetic so the same code runs under the JAX toolkit.
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

import pycircuit.circuit.circuit
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.analysis import NoConvergenceError, SingularMatrix
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import (SubCircuit, VS, R, Idt, Idtmod,
                                        floored_wrap)
from pycircuit.circuit.transient import Transient


def _ramp_circuit(**idtmod_params):
    """1 V into an Idtmod: the output is a unit-slope ramp, wrapped."""
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nin, nout = c.add_node('in'), c.add_node('out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, **idtmod_params)
    return c


def test_floored_wrap_congruence():
    """Range and congruence for negative arguments and non-multiple offsets."""
    for m, o in [(1.0, -0.5), (1.0, 0.0), (2.5, 0.3)]:
        for y in [-3.7, -0.5, -1e-3, 0.0, 0.2, 0.5, 5.3, 17.0]:
            k = floored_wrap(y, m, o, numeric)
            ## The ulp slack: at an exact multiple the float result may land
            ## one rounding step outside the half-open range (documented on
            ## floored_wrap); periodic consumers cannot see the difference.
            ulp = 1e-12 * max(abs(y), 1.0)
            assert o - ulp <= k < o + m + ulp, (y, m, o, k)
            n = (y - k) / m
            assert abs(n - round(n)) < 1e-12, (y, m, o, k)


def test_offset_congruence_transient():
    """Kundert's canonical arguments: modulus=1, offset=-0.5.

    The output must be the ramp folded into [-0.5, 0.5) -- i.e. congruent to
    the integral mod 1 -- not ``(ramp % 1) - 0.5``, which is what the old
    truncating-style wrap produced (a constant -0.5 error half the time).
    """
    c = _ramp_circuit(modulus=1.0, offset=-0.5)
    from pycircuit.circuit.integrator import EulerIntegrator
    ## Euler for the same reason as test_Idtmod_modulo: samples landing
    ## exactly ON the wrap take the right-limit convention under Euler.
    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=EulerIntegrator())
    result = tran.solve(tend=2.0, timestep=1e-2, fixed_timestep=True)
    y = result.v('out').y
    t = result.v('out').x[0]
    expected = ((t[1:] + 0.5) % 1.0) - 0.5
    ## Two-sided: at a sample landing exactly ON the wrap the sawtooth is
    ## double-valued and sub-ulp rounding picks the limit (see
    ## test_Idtmod_modulo).  Congruence is the well-defined comparison.
    d = np.abs(y[1:] - expected)
    d = np.minimum(d, 1.0 - d)
    assert_array_almost_equal(d, np.zeros_like(d))
    assert np.all(y >= -0.5 - 1e-12) and np.all(y < 0.5 + 1e-12)


def test_dc_pins_output_to_ic():
    """LRM: "the initial condition shall force the DC solution"."""
    ## Idt: output exactly ic.
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nin, nout = c.add_node('in'), c.add_node('out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['Idt'] = Idt(nin, gnd, nout, gnd, ic=2.75)
    res = DC(c, toolkit=numeric).solve()
    assert abs(res.v('out') - 2.75) < 1e-9

    ## Idtmod: output wrap(ic); ic=1.7 with modulus=1, offset=-0.5 -> -0.3.
    c = _ramp_circuit(modulus=1.0, offset=-0.5, ic=1.7)
    res = DC(c, toolkit=numeric).solve()
    assert abs(res.v('out') - (-0.3)) < 1e-9

    ## The flag is scoped: a leaked 'dc' would pin every integrator for the
    ## whole transient that follows (idtmod.md sec. 5.5).
    from pycircuit.circuit.circuit import defaultepar
    assert getattr(defaultepar, 'analysis_kind', None) is None


def test_dc_without_ic_stays_singular():
    """No ic: the LRM's no-ic branch -- a driven integrator has no DC
    operating point, and that must still fail loudly, not silently pin."""
    c = _ramp_circuit(modulus=1.0)
    with pytest.raises((NoConvergenceError, SingularMatrix)):
        DC(c, toolkit=numeric).solve()


def test_transient_from_dc_ic():
    """With ic given, `uic=True` is no longer mandatory: the transient's
    operating point comes from the pin and integration continues from it."""
    c = _ramp_circuit(modulus=1.0, offset=-0.5, ic=0.25)
    from pycircuit.circuit.integrator import EulerIntegrator
    tran = Transient(c, toolkit=numeric, integrator=EulerIntegrator())
    result = tran.solve(tend=0.2, timestep=1e-2, fixed_timestep=True)
    y = result.v('out').y
    t = result.v('out').x[0]
    ## Ramp from the pinned value; stays below the +0.5 boundary throughout.
    assert_array_almost_equal(y[1:], 0.25 + t[1:])


def test_uic_seeds_state_from_ic():
    """`uic=True` + element ic: the IC_KIND='state' route writes the private
    state row, wrapped into range, so the output starts at wrap(ic)."""
    c = _ramp_circuit(modulus=1.0, offset=-0.5, ic=1.25)  # wrap(1.25) = 0.25
    from pycircuit.circuit.integrator import EulerIntegrator
    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=EulerIntegrator())
    result = tran.solve(tend=0.2, timestep=1e-2, fixed_timestep=True)
    y = result.v('out').y
    t = result.v('out').x[0]
    assert_array_almost_equal(y[1:], 0.25 + t[1:])


def test_i_is_pure_under_jax_toolkit():
    """The old ``i()`` mutated its result vector in place, which raises on
    JAX's immutable arrays; the rewrite is pure array arithmetic and must
    return identical values under both toolkits."""
    pytest.importorskip('jax')
    from pycircuit.circuit.toolkit import jaxtoolkit

    x_vals = [0.1, 0.0, 0.3, 0.0, -1.7, 0.05]
    results = {}
    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        for key, tk in (('numeric', numeric), ('jax', jaxtoolkit)):
            pycircuit.circuit.circuit.default_toolkit = tk
            cir = SubCircuit(toolkit=tk)
            cir['E'] = Idtmod('iplus', 'iminus', 'oplus', 'ominus',
                              modulus=1.0, offset=-0.5, toolkit=tk)
            cir.update_iparv()
            el = cir['E']
            results[key] = np.asarray(el.i(tk.array(x_vals)), dtype=float)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved

    vals = list(results.values())
    assert_array_almost_equal(vals[0], vals[1])
    ## And the wrap itself is right: -s = 1.7 -> wrap = -0.3, so the branch
    ## row is -(v_op - v_on) + wrap = -0.3 - 0.3 = -0.6.
    assert abs(vals[0][-1] - (-0.6)) < 1e-12


def test_solve_batched_sweeps_modulus():
    """`eval_i_pure`/`eval_q_pure` admit Idtmod to the vmap groups, so a
    per-lane `modulus` override is accepted -- it used to be refused with
    NotImplementedError -- and each lane wraps at its own modulus."""
    pytest.importorskip('jax')
    import warnings
    import jax.numpy as jnp
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = jaxtoolkit
        c = SubCircuit(toolkit=jaxtoolkit)
        nin, nout = c.add_node('in'), c.add_node('out')
        c['vin'] = VS(nin, gnd, v=1.0, toolkit=jaxtoolkit)
        c['R1'] = R(nout, gnd, r=1e3, toolkit=jaxtoolkit)
        c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0,
                             toolkit=jaxtoolkit)

        tran = JAXTransient(c, reltol=1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve_batched(
                refnode=gnd,
                override_params_tree={
                    'Idtmod': {'modulus': jnp.array([[0.5], [1.0]])}},
                tend=1.4, timestep=1e-2, uic=True)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved

    assert len(res) == 2
    for lane, m in zip(res, (0.5, 1.0)):
        y = np.asarray(lane.v('out'), float).reshape(-1)
        assert y.max() < m + 1e-3
        assert y.min() >= -1e-3


## ------------------------------------------------------------------------
## Phase 2 (idtmod.md sec. 5.2): the post-acceptance gauge shift.


def _idt_state_row(c):
    """Row index of the element's private state in the assembled system."""
    rows = [i for i, node in enumerate(c.nodes) if 'idt_node' in str(node)]
    assert len(rows) == 1
    return rows[0]


def test_gauge_shift_keeps_state_bounded():
    """The declared row is rewrapped after every accepted step: over many
    wraps the recorded state spans about one modulus instead of growing
    with the integral, and the output is untouched."""
    c = _ramp_circuit(modulus=1.0, offset=-0.5)
    assert c.periodic_states() == [(_idt_state_row(c), 1.0, -0.5)]

    from pycircuit.circuit.integrator import EulerIntegrator
    tran = Transient(c, toolkit=numeric, uic=True,
                     integrator=EulerIntegrator())
    result = tran.solve(tend=5.0, timestep=1e-2, fixed_timestep=True)
    y = result.v('out').y
    t = result.v('out').x[0]
    d = np.abs(y[1:] - (((t[1:] + 0.5) % 1.0) - 0.5))
    d = np.minimum(d, 1.0 - d)
    assert d.max() < 1e-9
    s = np.asarray(result.x[_idt_state_row(c)], float)
    ## 5 wraps happened; unshifted the state would span ~5 moduli.
    assert s.max() - s.min() < 1.5


def test_long_run_precision_payoff():
    """The measured claim the Phase-2 recommendation stands on (idtmod.md
    5.4 item 4): the bounded state holds the phase to ~1 ulp regardless of
    run length, while the unbounded state's error grows with the integral.

    dt = 0.5 keeps the time grid exact in binary; v = 0.2 makes the
    per-step phase increment (0.1) inexact, so state arithmetic must round
    -- at the state's own magnitude.  The increment is 1/10 exactly, so
    the analytic phase at step k is (k mod 10)/10: an exact, cyclic,
    non-growing reference.  Measured at introduction: shift 2.2e-16 (flat
    from tend=500 to 2000), no-shift 1.4e-12 at tend=500 growing to
    2.2e-11 at tend=2000.
    """
    from pycircuit.circuit.integrator import EulerIntegrator

    def run(shift):
        pycircuit.circuit.circuit.default_toolkit = numeric
        c = SubCircuit()
        nin, nout = c.add_node('in'), c.add_node('out')
        c['vin'] = VS(nin, gnd, v=0.2)
        c['R1'] = R(nout, gnd, r=1e3)
        c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)
        tran = Transient(c, toolkit=numeric, uic=True,
                         integrator=EulerIntegrator())
        if not shift:
            ## Phase-1 behaviour: output wrap over an unbounded state.
            tran._apply_periodic_shifts = lambda x, X: None
        res = tran.solve(tend=500.0, timestep=0.5, fixed_timestep=True)
        y = res.v('out').y
        k = np.arange(len(y))
        expected = (k % 10) / 10.0
        d = np.abs(y[1:] - expected[1:])
        d = np.minimum(d, 1.0 - d)
        return d.max()

    err_shift = run(True)
    err_noshift = run(False)
    assert err_shift < 1e-14, err_shift
    assert err_noshift > 10 * err_shift, (err_noshift, err_shift)


def test_jax_gauge_shift():
    """The branchless shift inside the traced accept: bounded state and
    congruence-correct output on the JAX backend."""
    pytest.importorskip('jax')
    import warnings
    from pycircuit.circuit.toolkit import jaxtoolkit
    from pycircuit.circuit.jaxtransient import JAXTransient

    saved = pycircuit.circuit.circuit.default_toolkit
    try:
        pycircuit.circuit.circuit.default_toolkit = jaxtoolkit
        c = SubCircuit(toolkit=jaxtoolkit)
        nin, nout = c.add_node('in'), c.add_node('out')
        c['vin'] = VS(nin, gnd, v=1.0, toolkit=jaxtoolkit)
        c['R'] = R(nout, gnd, r=1e3, toolkit=jaxtoolkit)
        c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0, offset=-0.5,
                             toolkit=jaxtoolkit)
        tran = JAXTransient(c, reltol=1e-4)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(gnd, tend=3.0, timestep=1e-2, uic=True,
                             fixed_timestep=True)
    finally:
        pycircuit.circuit.circuit.default_toolkit = saved

    y = np.asarray(res.v('out'), float).reshape(-1)
    t = np.asarray(res.v('out').x[0], float).reshape(-1)
    d = np.abs(y[1:] - (((t[1:] + 0.5) % 1.0) - 0.5))
    d = np.minimum(d, 1.0 - d)
    assert d.max() < 1e-6
    s = np.asarray(res.x, float)[_idt_state_row(c)]
    assert s.max() - s.min() < 1.5


def test_periodic_contract_violation_raises():
    """A declaration whose row lacks the unit C diagonal must fail loudly at
    collection -- shifting it would corrupt the LTE estimate silently."""
    class BadIdtmod(Idtmod):
        def periodic_states(self):
            ## The output BRANCH row: its C row is all zeros, not a unit
            ## diagonal, so q[row] != x[row] and the shift contract fails.
            return [(self._branch_index, 1.0, -0.5)]

    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    nin, nout = c.add_node('in'), c.add_node('out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R1'] = R(nout, gnd, r=1e3)
    c['bad'] = BadIdtmod(nin, gnd, nout, gnd, modulus=1.0)
    tran = Transient(c, toolkit=numeric, uic=True)
    with pytest.raises(ValueError, match='periodic_states'):
        tran.solve(tend=0.1, timestep=1e-2, fixed_timestep=True)
