# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Stage 5 gate: does the algebra predict an actual simulated circuit?

Every other gate in ``doc/distortion_plan.md`` checks that Buonomo's algebra
was implemented faithfully -- against the Volterra series, against his own
closed forms, against his published curves.  **This one is the only check that
the algebra predicts a circuit at all.**

It runs pycircuit's own numeric transient engine on a biased diode, FFTs the
steady state, and compares the harmonics against the perturbation prediction.
Nothing on the measurement path shares code with the analysis under test: a
different solver, a different representation, a different failure mode.

The gate has two halves, and the second matters as much as the first.  The
prediction must agree *inside* the method's validity range, and it must
visibly *diverge* outside it -- otherwise the diagnostic that reports validity
would be unfalsifiable.

Getting here took two bug fixes, both found because this comparison refused to
agree: a phantom 1 mA DC offset in ``ISin`` (see ``test_source_dc_offset``) and
the need to solve the operating point for the whole circuit rather than the
diode alone.  Neither raised anything; both produced plausible wrong numbers.
"""

import numpy as np
import pytest
from scipy.optimize import brentq

from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit import circuit as circuit_module

## This whole module compares the analysis against transient simulation, so
## every test here pays for an integration.  Deselected by default
## via pytest.ini; run them with `-m slow` or `-m ""`.  See architecture P15.
pytestmark = pytest.mark.slow
from pycircuit.circuit.elements import Diode, ISin, R, C
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.distortion import harmonic_response, CubicNonlinearity


_IS_D = 1e-13
_IBIAS, _RLOAD, _CLOAD, _F0 = 1e-3, 1e3, 1e-9, 1e4


def _thermal_voltage():
    """kT/q from the toolkit's own constants -- what the Diode model uses.

    Not the round 25 mV: the difference is 3%, which exp() turns into a
    2.2x error in the Taylor coefficients and looks exactly like a defect.
    """
    return numeric.kboltzmann * 300 / numeric.qelectron


def _operating_point():
    """Junction voltage with the bias split between diode and resistor.

    Solving for the diode alone puts this ~20 mV too high, which the
    exponential magnifies into a 2.25x error in every coefficient.
    """
    VT = _thermal_voltage()
    return brentq(lambda v: v / _RLOAD + _IS_D * np.expm1(v / VT) - _IBIAS,
                  0.0, 1.0)


def _predict(signal_amplitude):
    VT = _thermal_voltage()
    v0 = _operating_point()
    gain = np.exp(v0 / VT)
    alpha = (_IS_D / VT) * gain                     # junction conductance
    b = (_IS_D / (2 * VT ** 2)) * gain
    c = (_IS_D / (6 * VT ** 3)) * gain

    w0 = 2 * np.pi * _F0
    Y = lambda s: 1 / _RLOAD + s * _CLOAD + alpha

    def apply_G(harmonic, rhs):
        return [np.asarray(rhs, dtype=complex)[0]
                / Y(1j * harmonic.frequency((w0,)))]

    return harmonic_response(apply_G, [signal_amplitude],
                             CubicNonlinearity(np.array([[b]]), np.array([[c]])),
                             tones=(w0,), toolkit=numeric)


_MEASUREMENT_CACHE = {}


def _measure(signal_amplitude, periods=40, points=512, keep=16):
    """HD2/HD3 and the DC level from a transient simulation and an FFT.

    Memoised: several tests interrogate the same operating point, and a
    transient run is by far the most expensive thing in this file.
    """
    key = (signal_amplitude, periods, points, keep)
    if key in _MEASUREMENT_CACHE:
        return _MEASUREMENT_CACHE[key]

    circuit_module.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)
    n1 = c.add_node('n1')
    c['I'] = ISin(gnd, n1, io=_IBIAS, ia=signal_amplitude, freq=_F0)
    c['R'] = R(n1, gnd, r=_RLOAD)
    c['C'] = C(n1, gnd, c=_CLOAD)
    c['D'] = Diode(n1, gnd, IS=_IS_D)

    ## timestep_max pins the old timestep-as-cap sampling density -- the
    ## FFT-based harmonic floors here need it (cap decoupled 2026-08-21).
    _ts = 1.0 / (_F0 * points)
    res = Transient(c, toolkit=numeric, timestep_max=_ts).solve(
        refnode=gnd, tend=periods / _F0, timestep=_ts)
    wave = res.v(n1)
    t = np.asarray(wave.x[0])
    v = np.asarray(wave.y)

    ## Resample the last whole periods onto an exact grid so the FFT bins
    ## land on the harmonics; the solver's own steps are not uniform.
    grid = np.linspace(t[-1] - keep / _F0, t[-1], keep * points, endpoint=False)
    spectrum = np.fft.rfft(np.interp(grid, t, v)) / (keep * points)
    k = keep
    result = (abs(spectrum[2 * k]) / abs(spectrum[k]),
              abs(spectrum[3 * k]) / abs(spectrum[k]),
              spectrum[0].real,
              2 * abs(spectrum[k]))
    _MEASUREMENT_CACHE[key] = result
    return result


def test_operating_point_and_fundamental_agree():
    """Before the harmonics: the linear picture must be right."""
    signal = 2e-5
    sol = _predict(signal)
    _, _, dc, fundamental = _measure(signal)

    assert abs(dc - _operating_point()) < 1e-4, (
        'simulated DC %.4f V vs predicted %.4f V' % (dc, _operating_point()))
    predicted_fundamental = abs(sol.amplitude((1,), 0))
    assert abs(fundamental / predicted_fundamental - 1) < 0.01


def test_harmonics_agree_inside_the_validity_range():
    """The gate. Independent solver, independent representation, 1% tolerance."""
    signal = 2e-5
    sol = _predict(signal)
    hd2, hd3, _, _ = _measure(signal)

    assert sol.perturbation_ratio(0) < 0.05, 'test point is not weakly nonlinear'
    assert abs(hd2 / sol.HD2(0) - 1) < 0.01, (
        'HD2: perturbation %.6e, transient %.6e' % (sol.HD2(0), hd2))
    assert abs(hd3 / sol.HD3(0) - 1) < 0.01, (
        'HD3: perturbation %.6e, transient %.6e' % (sol.HD3(0), hd3))


def test_prediction_visibly_fails_outside_the_validity_range():
    """The other half of the gate.

    A method that agreed everywhere would mean the diagnostic measures
    nothing.  Driven well past weak nonlinearity the prediction must be
    clearly wrong -- and it is, by tens of percent.
    """
    signal = 2e-3                      # ~100x the valid drive
    sol = _predict(signal)
    hd2, _, _, _ = _measure(signal)

    assert sol.perturbation_ratio(0) > 1.0, (
        'this drive should be far outside the valid range')
    assert abs(hd2 / sol.HD2(0) - 1) > 0.2, (
        'the prediction did not visibly fail where it should have')


@pytest.mark.parametrize('signal,max_ratio,max_error', [
    (2e-5, 0.02, 0.01),
    (1e-4, 0.10, 0.05),
])
def test_diagnostic_bounds_the_error(signal, max_ratio, max_error):
    """``perturbation_ratio`` must actually track the true error.

    A diagnostic nobody has calibrated is decoration.  This pins the
    relationship the docstring claims: the relative error stays at or below
    the order of the reported ratio.
    """
    sol = _predict(signal)
    hd2, _, _, _ = _measure(signal)

    ratio = sol.perturbation_ratio(0)
    error = abs(hd2 / sol.HD2(0) - 1)

    assert ratio < max_ratio, 'ratio %.4f above the band for this point' % ratio
    assert error < max_error, (
        'error %.3f%% exceeds the band the diagnostic implies (ratio %.4f)'
        % (100 * error, ratio))
