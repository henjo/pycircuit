# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Pins the P8 fix (architecture.md): the transient/Newton strategy layer
(``Integrator``, ``NonLinearSolver``, ``Scaler``) now follows the same
convention the toolkit layer already used -- pass an instance, not a string
through an if/elif factory keyed by name.

``StepController`` already worked this way before this fix (``transient.py``
constructs ``IntegralController()`` directly and only overrides it if the
caller injected one); this brings the other three ABCs in line with it,
a precedent already established in this codebase.

The string form is not kept as compatibility sugar: this project has never
had a version bump past ``0.0`` and no external distribution, so there is no
downstream API contract to protect. Passing the old string form now raises a
clear TypeError immediately, rather than failing confusingly deep inside a
Newton loop or a transient step (e.g. ``AttributeError: 'str' object has no
attribute 'solve_system'``).
"""
import pytest

from pycircuit.circuit import SubCircuit, R, C, VS, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.integrator import EulerIntegrator, Gear2Integrator
from pycircuit.circuit.nrsolver import StandardNewton, DampedNewton
from pycircuit.circuit.scaler import NoneScaler, RowMaxScaler


def _rc_circuit():
    c = SubCircuit()
    c['VS'] = VS(1, gnd, v=1.0)
    c['R'] = R(1, 2, r=1e3)
    c['C'] = C(2, gnd, c=1e-9)
    return c


def test_transient_defaults_to_euler_integrator():
    tran = Transient(_rc_circuit())
    assert isinstance(tran._get_integrator(), EulerIntegrator)


def test_transient_accepts_an_integrator_instance():
    tran = Transient(_rc_circuit(), integrator=Gear2Integrator())
    integrator = tran._get_integrator()
    assert isinstance(integrator, Gear2Integrator)


def test_transient_rejects_the_old_string_form_clearly():
    tran = Transient(_rc_circuit(), integrator='euler')
    with pytest.raises(TypeError, match='Integrator instance'):
        tran._get_integrator()


def test_analysis_defaults_to_standard_newton_and_no_scaling():
    dc = DC(_rc_circuit())
    assert isinstance(dc._get_nrsolver(), StandardNewton)
    assert isinstance(dc._get_scaler(), NoneScaler)


def test_analysis_accepts_strategy_instances():
    dc = DC(_rc_circuit(), nrsolver=DampedNewton(), scaler=RowMaxScaler())
    assert isinstance(dc._get_nrsolver(), DampedNewton)
    assert isinstance(dc._get_scaler(), RowMaxScaler)


def test_analysis_rejects_the_old_string_form_clearly():
    dc = DC(_rc_circuit(), nrsolver='standard')
    with pytest.raises(TypeError, match='NonLinearSolver instance'):
        dc._get_nrsolver()

    dc2 = DC(_rc_circuit(), scaler='none')
    with pytest.raises(TypeError, match='Scaler instance'):
        dc2._get_scaler()
