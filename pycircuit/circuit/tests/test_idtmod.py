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
    assert_array_almost_equal(y[1:], expected)
    assert np.all(y >= -0.5) and np.all(y < 0.5)


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
