# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Independent sources must inject only the current they are asked to.

``IS``'s DC parameter ``i`` defaulted to **1 mA**, not zero, while ``VS``'s
``v`` defaulted to zero.  Because ``IS.u`` returns ``iparv.i +
function.f(t)``, and ``ISin``'s waveform *already* carries its own offset
``io``, every sinusoidal current source silently added a spurious 1 mA of DC
on top of whatever it was configured with.

That is the kind of defect this suite exists to catch: nothing raised, every
waveform still looked like a sinusoid on the right frequency, and the only
symptom was that a biased nonlinear circuit sat at the wrong operating point.
It was found while cross-checking a distortion prediction against a transient
simulation, when the simulated DC level implied a diode current three times
larger than the source could supply.
"""

import numpy as np
import pytest

from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit.elements import IS, ISin, VS, VSin, R


def _source_drive(element):
    """The source vector a one-node circuit sees at t = 0."""
    circuit_module.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)
    n1 = c.add_node('n1')
    c['S'] = element(n1)
    return np.asarray(c.u(0.0, analysis='tran'), dtype=float)


def test_dc_current_source_defaults_to_no_current():
    """An unconfigured source must be inert, as ``VS`` already was."""
    assert _source_drive(lambda n: IS(gnd, n))[0] == 0.0
    assert _source_drive(lambda n: VS(n, gnd))[0] == 0.0


def test_dc_current_source_injects_what_it_is_given():
    drive = _source_drive(lambda n: IS(gnd, n, i=1e-3))
    assert abs(drive[0] + 1e-3) < 1e-15


@pytest.mark.parametrize('offset', [0.0, 1e-3, -2e-3])
def test_sinusoidal_current_source_offset_is_not_doubled(offset):
    """``ISin``'s DC must be exactly ``io``, not ``io`` plus a hidden default.

    The regression: with ``io = 1 mA`` the source delivered 2 mA of DC.  At
    ``io = 0`` it delivered 1 mA out of nothing, which is the more insidious
    case -- a source configured as pure AC quietly biased the circuit.
    """
    drive = _source_drive(lambda n: ISin(gnd, n, io=offset, ia=1e-4, freq=1e4))
    assert abs(drive[0] + offset) < 1e-15, (
        'ISin(io=%g) drove %g A of DC' % (offset, -drive[0]))


def test_sinusoidal_sources_agree_on_their_offset_convention():
    """``VSin`` and ``ISin`` must treat their offsets the same way.

    ``VSin`` was always correct because ``VS.v`` defaulted to zero; the
    asymmetry between the two defaults is what made only one of them wrong.
    """
    circuit_module.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)
    n1 = c.add_node('n1')
    c['V'] = VSin(n1, gnd, vo=1.5, va=0.1, freq=1e4)
    assert abs(c['V'].iparv.v) < 1e-15, 'VSin picked up a stray DC v'

    c2 = SubCircuit(toolkit=numeric)
    n2 = c2.add_node('n1')
    c2['I'] = ISin(gnd, n2, io=1.5e-3, ia=1e-4, freq=1e4)
    assert abs(c2['I'].iparv.i) < 1e-15, 'ISin picked up a stray DC i'


def test_biased_resistor_sits_where_ohms_law_says():
    """End-to-end: the defect showed up as a wrong operating point.

    A 1 mA source into 1 kOhm must sit at 1 V.  With the doubled offset it sat
    at 2 V -- a perfectly plausible number, which is why nothing caught it.
    """
    circuit_module.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)
    n1 = c.add_node('n1')
    c['I'] = ISin(gnd, n1, io=1e-3, ia=0.0, freq=1e4)
    c['R'] = R(n1, gnd, r=1e3)

    from pycircuit.circuit.dcanalysis import DC
    v = DC(c, refnode=gnd, toolkit=numeric).solve().x[c.get_node_index(n1)]
    assert abs(float(v) - 1.0) < 1e-9, 'node sits at %g V, expected 1 V' % v
