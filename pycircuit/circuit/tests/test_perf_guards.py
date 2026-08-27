"""A bound around the one published figure nothing else guards.

Roadmap sec. 38 recorded that no performance figure on this branch is
asserted by any test.  That was half right, and the half that was wrong
matters for reading this file:

* `fold_card` and `_autohold` already have DIRECTION guards.
  `test_hdl_fold_card.py::test_everything_outside_instance_is_folded`
  asserts the folded class emits FEWER primitives than the symbolic one,
  and `test_autohold.py::test_each_regulariser_holds_its_argument`
  asserts the regularisers really do hold.  Neither optimisation can
  silently stop working.
* Their MAGNITUDE is unguarded, and deliberately stays so.  The obvious
  strengthening -- turn `<` into `< 0.75 *` -- was measured and refused:
  on the Gummel-Poon card that test uses, folding moves the primitive
  count 275 -> 263, a 4.4% margin.  A bound there would pin noise.  The
  38% reduction that produces the published 1.54x lives on PSP, whose
  folded compile costs 19.5 s -- too slow to spend on every run.
* The END-TO-END DSL overhead had no guard of either kind, because it is
  call overhead and has no emitted-work proxy.  That is what this file
  adds.

Why this figure and not the others: `hdl.rst` advertised the overhead as
1.14x for weeks; re-measured on the same machine it was 1.25x.  The
benchmark that produced it was in the tree the whole time and passing --
it asserted that the two waveforms AGREED, and never that the ratio
held.  Agreement alone is not a performance test, so this file asserts
both, and neither assertion is load-bearing without the other: equality
without the ratio is what drifted, and the ratio without equality could
be "won" by quietly doing less work.

MEASUREMENT DESIGN.  Absolute timings are worthless here -- the machine
this was written on ranged from load 7 to load 17 within the hour.
Three choices make the number mean something, and a fourth
(`_require_a_still_machine`) decides when it means nothing at all:

* Compare a RATIO measured inside one process, not a wall time.  Even
  that is not fully load-immune: an early run of the benchmark at load 7
  reported 0.86x where a quiet machine reports 1.01x, because the
  benchmark times all of A and then all of B and slow drift lands
  entirely on B.  So the two sides are INTERLEAVED here, A B A B, and
  the minimum of each is taken -- noise only ever adds time, so the
  minimum is the cleanest estimator of the underlying cost.
* Keep the work large enough that the ratio is stable.  Measured at
  three sizes, interleaved, five pairs each:

      1000 steps   min-ratio 0.966   pairs 0.96 0.96 1.11 1.02 1.01
       300 steps   min-ratio 1.080   pairs 0.99 1.17 1.05 1.08 1.00
       150 steps   min-ratio 1.131   pairs 0.85 1.14 1.28 1.17 1.06

  Shrinking the run is not a free saving: it is noisier AND biased
  upward, because per-solve setup is amortised over less work.  1000
  steps is the smallest size whose ratio is worth asserting on.
* Confirm a breach before believing it -- a single reading flaked once
  in 27 runs, and the bound cannot be widened to absorb that without
  passing 1.22x, the state being guarded against.  See
  `_confirmed_ratio`.
"""

import time
import warnings

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import SubCircuit, VSin, R, C
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.transient import Transient

#: Steps in the ladder solve.  See the module docstring: below this the
#: ratio is both noisier and biased, so this is a floor, not a knob.
TEND = 1e-3
TIMESTEP = 1e-6

#: Interleaved A B pairs.  The minimum of each side is taken.
PAIRS = 3

#: The bound.  Measured min-ratio is 0.97-1.01 on a quiet machine and
#: individual pairs range 0.96-1.11; the state BEFORE the sec. 23-26
#: performance work was **1.22x**.  0.15 of headroom sits clear of the
#: observed spread and still fails on a return to that state.  If this
#: fires spuriously the fix is more PAIRS, not a bigger bound -- raising
#: it past 1.22 would make the test unable to fail for its own reason.
MAX_OVERHEAD = 1.15


def _ladder(mk_r, mk_c, stages=8):
    ## Setting a module global from a test is how the deleted root
    ## scratch files poisoned each other (see the 2026-08-27 cleanup), so
    ## it is worth saying why this one is safe: `numeric` is already the
    ## baseline every test starts from, and `conftest.py`'s autouse
    ## `reset_global_toolkit` restores it after each test regardless.
    ## This line only removes a dependence on execution order -- it is
    ## setting the global to the value it is required to already hold.
    cm.default_toolkit = numeric
    c = SubCircuit()
    prev = c.add_node('n0')
    c['vs'] = VSin(prev, gnd, va=1.0, freq=2e3)
    for k in range(stages):
        nxt = c.add_node('n%d' % (k + 1))
        c['R%d' % k] = mk_r(prev, nxt, r=1e3)
        c['C%d' % k] = mk_c(nxt, gnd, c=1e-8)
        prev = nxt
    return c, 'n%d' % stages


def _solve(mk_r, mk_c):
    """One timed solve.  Returns (wall, waveform, accepted steps)."""
    c, node = _ladder(mk_r, mk_c)
    tran = Transient(c, toolkit=numeric, uic=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        res = tran.solve(tend=TEND, timestep=TIMESTEP, fixed_timestep=True)
        wall = time.perf_counter() - t0
    return (wall, np.asarray(res.v(node).y, float),
            res.statistics.accepted_steps)


def _interleaved_ratio(a, b, pairs=PAIRS):
    """min(b)/min(a), measured A B A B so drift cannot land on one side.

    Returns the ratio and both sides' results, so the caller can check
    that the two did the SAME work before believing the number.
    """
    _solve(*a), _solve(*b)                       # warm caches, discard
    ta, tb, last_a, last_b = [], [], None, None
    for _ in range(pairs):
        wa, ya, na = _solve(*a)
        wb, yb, nb = _solve(*b)
        ta.append(wa)
        tb.append(wb)
        last_a, last_b = (ya, na), (yb, nb)
    return min(tb) / min(ta), last_a, last_b, (ta, tb)


#: How far the reference side's own timings may spread before the
#: measurement is abandoned as unreadable.  See `_require_a_still_machine`.
MAX_DRIFT = 1.25


def _require_a_still_machine(ta, tb):
    """Skip rather than assert when the machine is not holding still.

    Measured while the test suite itself was running under `-n 8`, the
    reference side's eight timings came out

        2.384  1.206  1.208  0.807  0.795  0.866  0.756  0.777

    -- a **3.15x** ramp inside a single measurement, as the parallel run
    finished and handed back cores.  A ratio taken across that is not a
    slow measurement, it is a meaningless one: the two sides sample
    different parts of the ramp, and the reading swung 0.86x to 1.35x on
    consecutive trials with nothing changed.  1.35x would have failed the
    1.15x bound and reported a regression that does not exist.

    Note what does NOT work as a gate: comparing the reference against
    ITSELF.  That was tried first and reads 0.994-0.996 even in the
    middle of the ramp, because both sides are the same computation and
    the bias cancels exactly where it needs to show.  The instrument
    looks perfect precisely when it is least trustworthy.

    So the gate is on the SPREAD of the reference side's own timings,
    which is a direct measurement of whether the machine held still.
    """
    drift = max(max(ta) / min(ta), max(tb) / min(tb))
    if drift > MAX_DRIFT:
        pytest.skip(
            'machine not stable enough to measure a %.0f%% effect: the '
            'reference timings spread %.2fx within one measurement '
            '(limit %.2fx). This is a skip, not a pass -- re-run on an '
            'idle machine, or serially rather than under -n 8.'
            % ((MAX_OVERHEAD - 1) * 100, drift, MAX_DRIFT))


def _confirmed_ratio(a, b, bound, want='below'):
    """Measure, and re-measure once before believing a breach.

    Measured flake rate of a single reading against the bound below:
    **one failure in 27 runs**, roughly 4%.  That is low enough to look
    fine in development and high enough to fire every dozen or so full
    suite runs -- which is how a performance test earns a reputation for
    crying wolf and gets deleted.  Widening the bound instead is not
    available: 1.22x is the state being guarded against, so there is
    nowhere above 1.15x to go without making the test unable to fail for
    its own reason.

    So a breach is confirmed rather than believed.  Two independent
    excursions in the same direction take the false-alarm rate to the
    square of one, while a real regression -- which is not an excursion
    but a level shift -- still breaches both times.  The second
    measurement is only paid for when the first says something
    interesting, so the common path stays at one measurement.

    `want` is where the caller needs the ratio to sit: 'below' for the
    guard, 'above' for the bite-check.  On a re-measure the LESS damning
    of the two readings is kept, which is the conservative choice for
    both -- it can only make this function slower to complain.
    """
    assert want in ('below', 'above'), want
    ratio, ra, rb, (ta, tb) = _interleaved_ratio(a, b)
    breached = ratio >= bound if want == 'below' else ratio <= bound
    if breached:
        ## Only gate when about to complain: a reading that is already
        ## comfortable needs no defending, and skipping a happy test
        ## would only hide the guard.
        _require_a_still_machine(ta, tb)
        again, ra, rb, (ta2, tb2) = _interleaved_ratio(a, b)
        _require_a_still_machine(ta2, tb2)
        ratio = min(ratio, again) if want == 'below' else max(ratio, again)
    return ratio, ra, rb


@pytest.mark.slow
def test_the_dsl_stays_at_parity_with_hand_written_elements():
    """The published claim: DSL overhead is at parity (1.01x).

    Asserted together with waveform equality, because either alone is
    the bug that let `hdl.rst` drift from 1.14x to 1.25x unnoticed.
    """
    ratio, (y_ref, n_ref), (y_hdl, n_hdl) = _confirmed_ratio(
        (R, C), (eh.RHdl, eh.CHdl), MAX_OVERHEAD, want='below')

    ## Equal work, or the ratio compares two different computations.
    assert n_ref == n_hdl, (n_ref, n_hdl)
    assert float(np.max(np.abs(y_hdl - y_ref))) < 1e-12

    assert ratio < MAX_OVERHEAD, (
        'DSL overhead %.3fx exceeds the %.2fx bound; published parity is '
        '1.01x and the pre-optimisation state was 1.22x. Re-run '
        'benchmarks/hdl_overhead.py before adjusting this bound.'
        % (ratio, MAX_OVERHEAD))


class _SlowR(R):
    """A resistor with a deliberate, calibrated overhead per evaluation.

    Exists only so the guard above can be shown to bite.  The work is
    real numpy on a real array -- not `sleep` -- so it costs the same
    kind of time the thing being guarded costs.

    Calibrated, because the first attempt was too small to breach the
    bound and the bite-check said so (1.096x against a 1.15x bound).
    What it took to fix is worth knowing, since it is the same effect
    this whole file is guarding:

        1 dot of 2000 elements   ->  1.089x
        3 dots of 400 elements   ->  1.167x
        5 dots of 400 elements   ->  1.398x

    Twenty-five times less arithmetic, and it costs MORE.  The price is
    the per-call numpy entry, not the flops -- which is exactly why the
    DSL's overhead is a call-count question and has no emitted-work
    proxy to assert on instead.
    """

    _WASTE = np.ones(400)
    _CALLS = 5

    def _burn(self):
        for _ in range(self._CALLS):
            float(np.dot(self._WASTE, self._WASTE))

    def G(self, x, epar=None):
        self._burn()
        return super(_SlowR, self).G(x, epar) if epar is not None \
            else super(_SlowR, self).G(x)

    def i(self, x, epar=None):
        self._burn()
        return super(_SlowR, self).i(x, epar) if epar is not None \
            else super(_SlowR, self).i(x)


@pytest.mark.slow
def test_the_parity_guard_can_actually_detect_a_regression():
    """Bite-check: the measurement must fail when there IS a slowdown.

    Five tests on this branch passed with the feature they were meant to
    protect deleted, so a guard that has never been shown to fire is not
    yet evidence.  This drives the SAME measurement with a deliberately
    slowed element and requires it to breach the same bound -- so a
    green run of the test above means the harness looked and found
    nothing, not that it cannot see.
    """
    ratio, (_, n_ref), (_, n_slow) = _confirmed_ratio(
        (R, C), (_SlowR, C), MAX_OVERHEAD, want='above')
    assert n_ref == n_slow, (n_ref, n_slow)
    assert ratio > MAX_OVERHEAD, (
        'the injected slowdown measured only %.3fx, which is under the '
        '%.2fx bound -- this harness cannot currently detect a regression '
        'of the size it claims to guard against' % (ratio, MAX_OVERHEAD))
