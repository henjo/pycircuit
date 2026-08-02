"""STAGE 10.3 -- `ic`, initial node voltages for `uic=True`.

`uic=True` used to mean "start from a vector of zeros", which is not what SPICE
means by it and leaves a class of circuit **unsimulable rather than merely
inconvenient**: an LC tank at zero is *at* an equilibrium and stays there
forever, and a latch at zero sits on its metastable point. Neither could be
started at all.

Both are tested here, because they are the two shapes the feature exists for.
"""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, L, VCVS
from pycircuit.circuit.transient import Transient


def _lc_tank(l=1e-6, c=1e-9):
    ck = SubCircuit()
    ck['L'] = L(1, gnd, L=l)
    ck['C'] = C(1, gnd, c=c)
    return ck


def test_an_lc_tank_cannot_start_from_zero():
    """The defect, stated as a test: zero is an equilibrium.

    Without this the feature looks optional. It is not -- there is no argument to
    `solve()` that makes this circuit oscillate.
    """
    tran = Transient(_lc_tank(), toolkit=numeric, reltol=1e-6, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=2e-7, timestep=1e-9)
    v = np.asarray(res.v(1, gnd).y, dtype=float).ravel()
    assert np.max(np.abs(v)) < 1e-12, \
        'a tank started at zero moved to %g V' % np.max(np.abs(v))


def test_an_lc_tank_started_by_ic_oscillates_at_its_own_frequency():
    """And the frequency is checked against `1/(2*pi*sqrt(LC))`, not a fitted
    number -- so the test cannot pass on a circuit that merely wobbles."""
    L_, C_ = 1e-6, 1e-9
    f0 = 1.0 / (2 * np.pi * np.sqrt(L_ * C_))
    tran = Transient(_lc_tank(L_, C_), toolkit=numeric, reltol=1e-7,
                     uic=True, ic={'1': 1.0})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=4.0 / f0, timestep=1e-10)
    t = np.asarray(res.v(1, gnd).x, dtype=float).ravel()
    v = np.asarray(res.v(1, gnd).y, dtype=float).ravel()

    assert v[0] == pytest.approx(1.0, abs=1e-9), 'the ic was not applied'
    assert np.max(v) > 0.5 and np.min(v) < -0.5, 'it did not oscillate'

    ## Zero crossings give the period without fitting anything.
    sign = np.sign(v)
    crossings = t[:-1][sign[:-1] * sign[1:] < 0]
    assert len(crossings) >= 4, 'only %d zero crossings' % len(crossings)
    half_periods = np.diff(crossings)
    measured = 1.0 / (2.0 * np.mean(half_periods))
    assert measured == pytest.approx(f0, rel=0.02), \
        'measured %.4g Hz against %.4g Hz' % (measured, f0)


def test_a_latch_started_by_ic_leaves_its_metastable_point():
    """Cross-coupled inverting stages, each driving an RC. Zero is metastable.

    THE TEST CIRCUIT WAS WRONG TWICE BEFORE IT WAS RIGHT, and both mistakes
    produce the same symptom -- v(1) sits at zero -- which is also what a broken
    `ic` would produce. Recorded because that ambiguity is the whole difficulty:

    * cross-coupling two ideal `VCVS` directly makes the node voltages
      ALGEBRAICALLY determined, so `v1 = 9*v1` forces `v1 = 0` and there is no
      state for an initial condition to set. The capacitors must sit behind a
      resistance to be states at all.
    * `VCVS` is `(inp, inn, outp, outn)`. Written the other way round, each
      source DRIVES the node it was meant to sense: the first solve then pushed
      1 A into node 1 and discharged its capacitor completely in a single
      1e-15 s step.

    In both cases `_initial_state` was returning the right vector the whole time.
    """
    ck = SubCircuit()
    ## Each stage: a gain of -3 into an RC, so v(1) and v(2) are real states.
    ## VCVS is (inp, inn, outp, outn) -- INPUT first. Writing it the other way
    ## round makes each source DRIVE the node it was meant to sense, which pins
    ## the state node instead of feeding it and looks exactly like "the latch
    ## stayed at zero".
    ck['E1'] = VCVS(2, gnd, 'a', gnd, g=-3.0)     # sense v(2), drive 'a'
    ck['R1'] = R('a', 1, r=1e3)
    ck['C1'] = C(1, gnd, c=1e-12)
    ck['E2'] = VCVS(1, gnd, 'b', gnd, g=-3.0)     # sense v(1), drive 'b'
    ck['R2'] = R('b', 2, r=1e3)
    ck['C2'] = C(2, gnd, c=1e-12)

    tran = Transient(ck, toolkit=numeric, reltol=1e-6, uic=True,
                     ic={'1': 1e-3})
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=5e-9, timestep=1e-12)
    v1 = np.asarray(res.v(1, gnd).y, dtype=float).ravel()

    assert abs(v1[0]) > 0, 'the ic was not applied at all'
    assert abs(v1[-1]) > 10 * abs(v1[0]), \
        'the latch did not diverge from metastability: %g -> %g' % (v1[0], v1[-1])


def test_naming_a_node_that_does_not_exist_raises():
    """A typo in a node name must not silently apply nothing."""
    tran = Transient(_lc_tank(), toolkit=numeric, uic=True,
                     ic={'nosuch': 1.0})
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=1e-9, timestep=1e-10)
    assert 'nosuch' in str(exc.value)


def test_setting_the_reference_node_raises():
    """It is held at 0 V by construction, so this is a request that cannot be
    honoured -- and silently dropping it would leave the caller believing an
    initial condition was applied."""
    tran = Transient(_lc_tank(), toolkit=numeric, uic=True, ic={gnd: 1.0})
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=1e-9, timestep=1e-10)
    assert 'reference node' in str(exc.value)


def test_ic_without_uic_raises():
    """SPICE's `.ic` without UIC constrains the operating point and then
    releases it -- a different feature, and not this one. Silently doing neither
    is the defect class this project keeps finding."""
    tran = Transient(_lc_tank(), toolkit=numeric, ic={'1': 1.0})
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=1e-9, timestep=1e-10)
    assert 'uic=True' in str(exc.value)
