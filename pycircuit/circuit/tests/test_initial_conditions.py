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


## ---------------------------------------------------------------------------
## STAGE 10.3, element initial conditions: `L(..., ic=...)`.
##
## Unlike a capacitor's initial voltage -- which constrains a DIFFERENCE of two
## node unknowns and needs a spanning-tree solve -- an inductor's is a branch
## CURRENT whose unknown already exists in the MNA vector. It is an assignment
## once the row is known, and finding the row is the whole of the work.

def test_an_inductor_ic_lands_on_its_own_branch_row():
    ck = SubCircuit()
    ck['L'] = L(1, gnd, L=1e-6, ic=0.5)
    ck['C'] = C(1, gnd, c=1e-9)
    tran = Transient(ck, toolkit=numeric, uic=True)
    x0 = tran._initial_state(gnd)
    row = ck.instance_branch_indices('L')[0]
    assert x0[row] == pytest.approx(0.5)
    ## and nothing else was touched
    assert np.count_nonzero(x0) == 1


def test_parallel_inductors_get_their_own_rows():
    """THE REASON THE SPAN MAP IS RECORDED RATHER THAN RECONSTRUCTED.

    `Branch.__eq__` compares node pairs, so two inductors between the same two
    nodes produce EQUAL branches and `branches.index()` returns the first for
    both -- measured, before the map existed, as `get_branch_index` giving row 2
    for each. An initial current given to the second would have landed on the
    first one's unknown with nothing to indicate it.
    """
    ck = SubCircuit()
    ck['L1'] = L(1, gnd, L=1e-6, ic=0.25)
    ck['L2'] = L(1, gnd, L=2e-6, ic=-0.75)
    ck['C'] = C(1, gnd, c=1e-9)

    r1 = ck.instance_branch_indices('L1')[0]
    r2 = ck.instance_branch_indices('L2')[0]
    assert r1 != r2, 'parallel inductors share a row'

    tran = Transient(ck, toolkit=numeric, uic=True)
    x0 = tran._initial_state(gnd)
    assert x0[r1] == pytest.approx(0.25)
    assert x0[r2] == pytest.approx(-0.75)


def test_an_element_ic_without_uic_raises():
    """Same rule as the analysis-level `ic`: a starting value that the operating
    point would overwrite must not be silently accepted."""
    ck = SubCircuit()
    ck['L'] = L(1, gnd, L=1e-6, ic=0.5)
    ck['C'] = C(1, gnd, c=1e-9)
    tran = Transient(ck, toolkit=numeric)
    with pytest.raises(ValueError) as exc:
        tran.solve(tend=1e-9, timestep=1e-10)
    assert 'uic=True' in str(exc.value)


def test_an_ic_inside_a_subcircuit_is_refused_not_ignored():
    """A nested instance owns a span covering all its children's branches, so a
    flat walk cannot place an ic inside one. Refused with the reason, because
    accepting the parameter and ignoring it is the defect this stage exists to
    stop repeating."""
    ## A SubCircuit used as an element must be a class with `terminals`,
    ## instantiated with the nodes it connects to -- assigning a bare
    ## `SubCircuit()` raises KeyError from `_instance_nodes`, because nothing
    ## maps its terminal names onto parent nodes.
    class Inner(SubCircuit):
        terminals = ('p', 'n')

        def __init__(self, *args, **kvargs):
            super().__init__(*args, **kvargs)
            self['L'] = L('p', 'n', L=1e-6, ic=0.5)

    ck = SubCircuit()
    ck['X1'] = Inner(1, gnd)
    ck['C'] = C(1, gnd, c=1e-9)

    tran = Transient(ck, toolkit=numeric, uic=True)
    with pytest.raises(NotImplementedError) as exc:
        tran._initial_state(gnd)
    assert 'subcircuit' in str(exc.value).lower()


def test_a_tank_started_by_inductor_current_has_the_analytic_amplitude():
    """All the energy in L, so v starts at ZERO and swings to I*sqrt(L/C).

    Distinguishable from the voltage-started tank by both initial value and
    amplitude, and the amplitude is analytic -- a misassigned branch row or a
    dropped ic shows up immediately and cannot be fitted away.
    """
    L_, C_, I0 = 1e-6, 1e-9, 0.5
    f0 = 1.0 / (2 * np.pi * np.sqrt(L_ * C_))
    expected = I0 * np.sqrt(L_ / C_)

    ck = SubCircuit()
    ck['L'] = L(1, gnd, L=L_, ic=I0)
    ck['C'] = C(1, gnd, c=C_)

    tran = Transient(ck, toolkit=numeric, reltol=1e-7, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=3.0 / f0, timestep=1e-10)
    v = np.asarray(res.v(1, gnd).y, dtype=float).ravel()

    ## The result excludes x0 -- its first sample is the first STEPPED point, at
    ## t = timestep*1e-3 -- so this is `small relative to the swing`, not zero.
    ## At 1e-13 s with 0.5 A into 1 nF the exact value is 5e-5 V, and asserting
    ## `< 1e-9` failed on the correct answer.
    assert abs(v[0]) < 1e-3 * expected, \
        'v should start near zero, got %g against a %g swing' % (v[0], expected)
    assert np.max(np.abs(v)) == pytest.approx(expected, rel=0.03), \
        'amplitude %.4g against an analytic %.4g' % (np.max(np.abs(v)), expected)


## ---------------------------------------------------------------------------
## STAGE 10.3, capacitor initial voltages: `C(..., ic=...)`.
##
## Not an assignment. A capacitor has NO state variable of its own -- `q` is
## derived from the node voltages -- so `ic` constrains a DIFFERENCE of two
## unknowns and a set of them is solved as a spanning tree. See
## doc/initial_conditions.md sec. 4a.

def test_a_grounded_capacitor_chain_propagates_from_ground():
    ck = SubCircuit()
    ck['C1'] = C('a', 'b', c=1e-9, ic=2.0)
    ck['C2'] = C('b', gnd, c=1e-9, ic=3.0)
    ck['R'] = R('a', gnd, r=1e9)
    tran = Transient(ck, toolkit=numeric, uic=True)
    x0 = tran._initial_state(gnd)
    assert x0[ck.get_node_index('b')] == pytest.approx(3.0)
    assert x0[ck.get_node_index('a')] == pytest.approx(5.0)


def test_a_floating_group_raises_rather_than_choosing_a_reference():
    """DECISION (a). The differences asked for are satisfied by infinitely many
    assignments, and the absolute values reach the output waveform -- so picking
    one silently is the quiet-wrong-answer shape this stage exists to avoid."""
    ck = SubCircuit()
    ck['Ca'] = C('p', 'q', c=1e-9, ic=1.5)   # touches neither ground nor an ic
    ck['C0'] = C(1, gnd, c=1e-9)
    tran = Transient(ck, toolkit=numeric, uic=True)
    with pytest.raises(ValueError) as exc:
        tran._initial_state(gnd)
    assert 'up to a constant' in str(exc.value)
    assert 'Ca' in str(exc.value)


def test_a_node_ic_can_anchor_an_otherwise_floating_group():
    """And it is the documented way out of the error above."""
    ck = SubCircuit()
    ck['Ca'] = C('p', 'q', c=1e-9, ic=1.5)
    ck['C0'] = C(1, gnd, c=1e-9)
    tran = Transient(ck, toolkit=numeric, uic=True, ic={'q': 0.5})
    x0 = tran._initial_state(gnd)
    assert x0[ck.get_node_index('q')] == pytest.approx(0.5)
    assert x0[ck.get_node_index('p')] == pytest.approx(2.0)


def test_a_consistent_loop_is_accepted():
    """v(a)-v(b) = 1, v(b)-v(c) = 2, v(a)-v(c) = 3 -- the loop sums to zero."""
    ck = SubCircuit()
    ck['Cab'] = C('a', 'b', c=1e-9, ic=1.0)
    ck['Cbc'] = C('b', 'c', c=1e-9, ic=2.0)
    ck['Cac'] = C('a', 'c', c=1e-9, ic=3.0)
    ck['Cg'] = C('c', gnd, c=1e-9, ic=0.0)
    tran = Transient(ck, toolkit=numeric, uic=True)
    x0 = tran._initial_state(gnd)
    assert x0[ck.get_node_index('c')] == pytest.approx(0.0)
    assert x0[ck.get_node_index('b')] == pytest.approx(2.0)
    assert x0[ck.get_node_index('a')] == pytest.approx(3.0)


def test_an_inconsistent_loop_raises():
    """The same loop with 3.5 instead of 3.0 is a contradiction the caller
    cannot have meant, and resolving it by insertion order would make the answer
    depend on the order elements were added."""
    ck = SubCircuit()
    ck['Cab'] = C('a', 'b', c=1e-9, ic=1.0)
    ck['Cbc'] = C('b', 'c', c=1e-9, ic=2.0)
    ck['Cac'] = C('a', 'c', c=1e-9, ic=3.5)
    ck['Cg'] = C('c', gnd, c=1e-9, ic=0.0)
    tran = Transient(ck, toolkit=numeric, uic=True)
    with pytest.raises(ValueError) as exc:
        tran._initial_state(gnd)
    assert 'contradictory' in str(exc.value)


def test_parallel_capacitors_that_disagree_raise():
    ck = SubCircuit()
    ck['C1'] = C('a', gnd, c=1e-9, ic=1.0)
    ck['C2'] = C('a', gnd, c=2e-9, ic=1.5)
    tran = Transient(ck, toolkit=numeric, uic=True)
    with pytest.raises(ValueError) as exc:
        tran._initial_state(gnd)
    assert 'contradictory' in str(exc.value)


def test_a_capacitor_ic_disagreeing_with_a_node_ic_raises():
    ck = SubCircuit()
    ck['C1'] = C('a', gnd, c=1e-9, ic=1.0)
    tran = Transient(ck, toolkit=numeric, uic=True, ic={'a': 2.0})
    with pytest.raises(ValueError) as exc:
        tran._initial_state(gnd)
    assert 'contradictory' in str(exc.value)


def test_both_terminals_on_one_node_must_have_a_zero_ic():
    ck = SubCircuit()
    ck['C1'] = C('a', 'a', c=1e-9, ic=1.0)
    ck['C2'] = C('a', gnd, c=1e-9)
    tran = Transient(ck, toolkit=numeric, uic=True)
    with pytest.raises(ValueError) as exc:
        tran._initial_state(gnd)
    assert '0 ==' in str(exc.value)


def test_inductor_and_capacitor_ics_coexist():
    """The two are solved by different mechanisms -- an assignment and a tree
    walk -- and must not disturb each other."""
    ck = SubCircuit()
    ck['L'] = L(1, gnd, L=1e-6, ic=0.25)
    ck['C'] = C(1, gnd, c=1e-9, ic=3.0)
    tran = Transient(ck, toolkit=numeric, uic=True)
    x0 = tran._initial_state(gnd)
    assert x0[ck.get_node_index(1)] == pytest.approx(3.0)
    assert x0[ck.instance_branch_indices('L')[0]] == pytest.approx(0.25)


def test_a_precharged_capacitor_discharges_with_the_analytic_time_constant():
    """End to end, against `v0*exp(-t/RC)` -- so a wrong sign or a dropped ic
    shows up in the waveform rather than only in `x0`."""
    R_, C_, V0 = 1e3, 1e-9, 2.0
    tau = R_ * C_
    ck = SubCircuit()
    ck['C'] = C(1, gnd, c=C_, ic=V0)
    ck['R'] = R(1, gnd, r=R_)
    tran = Transient(ck, toolkit=numeric, reltol=1e-7, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=4 * tau, timestep=tau / 50)
    t = np.asarray(res.v(1, gnd).x, dtype=float).ravel()
    v = np.asarray(res.v(1, gnd).y, dtype=float).ravel()
    assert np.max(np.abs(v - V0 * np.exp(-t / tau))) < 2e-3 * V0
