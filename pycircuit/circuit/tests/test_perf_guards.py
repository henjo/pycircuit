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
(reproducibility, read by `_verdict`) decides when it means nothing at
all:

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

⚠ REVISED 2026-08-31, and the revision contradicts a claim above.  The
skip message this file used to print advised re-running "on an idle
machine", which assumed one would be available; on a shared machine it
is not, and a guard that only works on a quiet box does not guard.
Measured under four concurrent openEMS runs (load 17 of 24 cores), the
minimum-of-N estimator DOES converge there -- it just needs samples:
three independent trials read 1.008/0.916/1.054 at 3 samples a side and
1.018/1.007/1.003 at 12.  So `PAIRS` rose to 12, and the readability
gate moved off "did the machine hold still" (which reads 1.86x-3.35x in
exactly those converged trials) and onto "did the estimate reproduce".
The load is a reason to take more samples, not a reason to abstain.
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

#: Interleaved A B pairs per measurement.  The minimum of each side is
#: taken; noise only ever adds time, so the minimum converges from above
#: and the only question is how many samples it costs to reach the floor.
#:
#: MEASURED 2026-08-31 on a machine at load 17 of 24 cores (four openEMS
#: runs), which is the adverse condition worth calibrating against
#: rather than an idle machine that was not available.  Three
#: independent trials, min-of-N ratio against sample count:
#:
#:     samples/side    trial 0   trial 1   trial 2    spread
#:        3             1.008     0.916     1.054      14%
#:       12             1.018     1.007     1.003     1.5%
#:
#: At 3 the spread is the ENTIRE headroom of the 1.15 bound, and trial 1
#: sat at 0.916 until its twelfth sample -- a 1.20x regression would
#: have read 1.09x there and passed.  At 12 the three agree to 1.5% and
#: centre on the published 1.01x.  Under load the floor is found late,
#: and finding it is the whole measurement.
PAIRS = 12

#: ⚠ TRIED AND REFUSED, same day.  A `PAIRS_COARSE = 3` was given to the
#: bite-check on the reasoning that it resolves a ~1.40x injected
#: slowdown rather than a 15% effect, so it should need fewer samples.
#: It skipped on its first run: two measurements disagreed by **1.147x**
#: against the 1.05x limit.
#:
#: The reasoning conflated two different things.  The readability gate
#: constrains the PRECISION OF THE ESTIMATE, and at 3 samples a side
#: that precision is +/-14% (measured, see `PAIRS`) no matter how large
#: the effect being measured is.  A coarse effect needs fewer samples to
#: SEPARATE; it needs exactly as many to REPRODUCE.  Both callers
#: therefore use the same PAIRS.

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


#: How far two INDEPENDENT measurements of the same ratio may disagree
#: before the pair is abandoned as unreadable.
#:
#: This replaced a gate on `_drift` -- the spread of the raw timings
#: within one measurement -- on 2026-08-31, because that statistic
#: measures the MACHINE and the question is whether the ESTIMATE
#: reproduces.  Measured, they disagree completely: in the three trials
#: recorded at `PAIRS`, raw-sample spread read 1.86x to 3.35x while the
#: min-of-12 estimate those same samples produced agreed across trials
#: to 1.5%.  Gating on drift would have abstained on every one of them,
#: and gating on it at PAIRS=3 abstained on a real regression (below).
#:
#: 1.05 is not a taste.  The assertion is `ratio < 1.15` against a true
#: ~1.01, so the estimate must be good to well inside 14%; two estimates
#: agreeing to 5% make it so.  Loosening to 1.15 would admit both
#: failure modes at once -- a true 1.01 reading 1.16 (false alarm) and a
#: true 1.20 reading 1.05 (false pass) -- which is the check paying for
#: itself with the thing it was bought for.  So this stays, and SAMPLES
#: are the knob.
MAX_DISAGREEMENT = 1.05

#: Measurements attempted before abstaining.  Measured 2026-08-31 under
#: load 17, disagreement between consecutive N-sample estimates:
#:
#:     N= 3   1.098 1.020 1.163 1.045 1.104 1.049 1.118 1.020 ...  8/16 within 1.05
#:     N= 6   1.131 1.084 1.025 1.054 1.009 1.058 1.013 1.022     4/8
#:     N=12   1.034 1.049 1.029 1.004                             4/4
#:     N=24   1.081 1.032                                         1/2
#:
#: Agreement improves with N but does not become RELIABLE: reaching 1.05
#: dependably would take ~48 samples a side, minutes per test, on a
#: suite that already runs 9-10 min.  Retrying is the cheaper lever --
#: a busy machine has quiet moments, and a second or third attempt
#: usually lands in one.  Cost is paid only when the first pair
#: disagrees, so a quiet machine still pays for exactly two.
MAX_MEASUREMENTS = 4


def _drift(ta, tb):
    """Spread of the raw timings within one measurement.

    Kept because it is worth REPORTING when a measurement is thrown
    away -- it says whether the machine or the method was at fault --
    but it is no longer what any decision is taken on.  See
    `MAX_DISAGREEMENT` for why, and note what does NOT work as a gate:
    comparing the reference against ITSELF reads 0.994-0.996 even in the
    middle of a ramp, because both sides are the same computation and
    the bias cancels exactly where it needs to show.
    """
    return max(max(ta) / min(ta), max(tb) / min(tb))


def _disagreement(r1, r2):
    """How far two independent estimates of the same ratio sit apart."""
    return max(r1, r2) / min(r1, r2)


def _verdict(ratio, disagreement, bound, want):
    """'unreadable' | 'breach' | 'ok', decided without touching a clock.

    READABILITY IS DECIDED FIRST, AND IN BOTH DIRECTIONS.  That ordering
    is the fix of 2026-08-31.  Until then the stillness gate sat inside
    the `if breached:` arm of `_confirmed_ratio`, which left two holes
    pulling in opposite directions:

    * A REAL REGRESSION WAS CONVERTED INTO A SKIP.  Breach, then gate,
      then abstain -- so on any machine too busy to measure, the one
      outcome this file exists to report was the one it could not
      report.  Demonstrated on a machine at load 17: driving the guard's
      own path with the calibrated `_SlowR` raised `Skipped` where a
      failure was due.
    * A COMFORTABLE READING WAS NEVER CHECKED AT ALL, on the argument
      that "a reading that is already comfortable needs no defending".
      That argument is wrong, because the bias runs BOTH ways: measured
      at PAIRS=3 under load, the same parity read 0.916x to 1.054x
      around a true ~1.01x, so a 1.20x regression reads 1.09x and
      passes -- silently, exactly when the machine is least able to see
      it.

    ⚠ And the bite-check covers neither hole.  It drives `want='above'`,
    where a successful bite means breached=False and returns before any
    gate is consulted -- a different path from the one a real regression
    takes.  `test_the_verdict_*` below cover this decision directly with
    synthetic values, which is also the only way to exercise it without
    owning the machine's load.
    """
    assert want in ('below', 'above'), want
    if disagreement > MAX_DISAGREEMENT:
        return 'unreadable'
    if want == 'below':
        return 'breach' if ratio >= bound else 'ok'
    return 'breach' if ratio <= bound else 'ok'


def _abstain(disagreement, drift):
    """Decline to report, and make the declining visible.

    A silent skip on a guard is indistinguishable from a skip for a
    boring reason -- an uninstalled optional dependency, say -- and this
    one means something quite different: the published figure went
    UNCHECKED on that run.  So it warns before it skips, which puts it
    in the warnings summary of a plain run rather than only behind
    `-rs`.

    Both numbers are reported because they say different things: the
    disagreement is why the measurement was thrown away, the raw drift
    says whether the machine was the reason.
    """
    reason = (
        'the %.0f%% overhead bound was NOT checked on this run: two '
        'independent measurements disagreed by %.3fx (limit %.2fx), so '
        'neither is assertable. Raw timing spread within a measurement '
        'was %.2fx. This is a skip, not a pass -- re-run when the '
        'machine is quieter, or raise PAIRS.'
        % ((MAX_OVERHEAD - 1) * 100, disagreement, MAX_DISAGREEMENT, drift))
    warnings.warn(reason, RuntimeWarning, stacklevel=2)
    pytest.skip(reason)


def _confirmed_ratio(a, b, bound, want='below', pairs=PAIRS):
    """Measure twice, and assert only on a number that reproduced.

    The original measured once and re-measured only to confirm a breach,
    which made the second measurement a defence against false ALARMS
    alone.  Measured flake rate of a single reading against the bound
    was one in 27 runs, and that reasoning still holds -- but it is
    half the problem.  The other half is a reading that is wrong in the
    comfortable direction, which no amount of re-measuring-on-breach can
    see, because it never breaches.

    So both measurements are taken always, and their agreement is the
    readability gate.  Two independent excursions in the same direction
    are what a false alarm would need; a real regression is not an
    excursion but a level shift, and reproduces.

    `want` is where the caller needs the ratio to sit: 'below' for the
    guard, 'above' for the bite-check.  The LESS damning of the two
    readings is kept, which is conservative for both -- it can only make
    this function slower to complain.
    """
    assert want in ('below', 'above'), want
    seen, drifts, disagreement = [], [], None
    for _attempt in range(MAX_MEASUREMENTS):
        r, ra, rb, (ta, tb) = _interleaved_ratio(a, b, pairs)
        seen.append((r, ra, rb))
        drifts.append(_drift(ta, tb))
        if len(seen) < 2:
            continue
        ## The two MOST RECENT, not the best pair out of all of them:
        ## picking the closest two from a growing set would let the
        ## harness shop for agreement, and two readings agree most
        ## easily when both are wrong in the same direction.
        r_prev, r_now = seen[-2][0], seen[-1][0]
        disagreement = _disagreement(r_prev, r_now)
        if _verdict(r_now, disagreement, bound, want) != 'unreadable':
            return ((min(r_prev, r_now) if want == 'below'
                     else max(r_prev, r_now)), ra, rb)
    _abstain(disagreement, max(drifts))


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


## ---------------------------------------------------------------------
## The guard's own decision, tested without a clock.
##
## Everything above needs a machine quiet enough to measure a 15% effect.
## These do not: `_verdict` is a pure function of (ratio, drift), so the
## paths that a real regression takes can be exercised on any machine, at
## any load, in microseconds.  That matters beyond convenience -- the
## bite-check above drives `want='above'`, which returns before the gate
## is consulted, so until these existed NOTHING covered the path a real
## regression actually takes.  See `_verdict.__doc__`.
## ---------------------------------------------------------------------

@pytest.mark.parametrize('ratio,disagreement,want,expected', [
    ## A regression that reproduced: report it.
    (1.30, 1.01, 'below', 'breach'),
    ## THE HOLE THIS CLOSES.  Same regression, but the two measurements
    ## disagree, so neither is assertable: abstain.  Before the fix the
    ## gate sat downstream of the breach and turned this into a SKIP --
    ## the one outcome the file exists to report was the one it could
    ## not report.
    (1.30, 1.40, 'below', 'unreadable'),
    ## Parity that reproduced: pass.
    (0.95, 1.01, 'below', 'ok'),
    ## THE OTHER HOLE.  A comfortable reading that did NOT reproduce is
    ## not a pass, it is an unread instrument -- the bias runs both ways,
    ## and at PAIRS=3 under load a true 1.01x read anywhere in
    ## 0.916-1.054.  This used to be an unconditional pass.
    (0.95, 1.40, 'below', 'unreadable'),
    ## The bite-check's direction: a successful bite is 'ok'...
    (1.40, 1.01, 'above', 'ok'),
    ## ...a bite that failed to bite is the breach...
    (1.00, 1.01, 'above', 'breach'),
    ## ...and it is gated on reproducibility too, in both outcomes.
    (1.40, 1.40, 'above', 'unreadable'),
    (1.00, 1.40, 'above', 'unreadable'),
])
def test_the_verdict_gates_on_readability_before_comparing(
        ratio, disagreement, want, expected):
    """Readability is decided first, and for every outcome.

    The parametrisation is the specification: two of these eight rows
    are the 2026-08-31 fix, and both were silent passes or silent skips
    before it.
    """
    assert _verdict(ratio, disagreement, MAX_OVERHEAD, want) == expected


def test_the_verdict_treats_the_bound_as_a_breach_not_a_pass():
    """Exactly at the bound is a breach, not a pass.

    The assertion the guard finally makes is `ratio < MAX_OVERHEAD`, so
    a verdict of 'ok' at exactly the bound would hand the test a value
    its own assert then rejects -- a disagreement that would surface as
    a confusing failure rather than a confirmed one.
    """
    assert _verdict(MAX_OVERHEAD, 1.0, MAX_OVERHEAD, 'below') == 'breach'
    assert _verdict(MAX_OVERHEAD, 1.0, MAX_OVERHEAD, 'above') == 'breach'


def test_disagreement_is_symmetric_and_bottoms_out_at_one():
    """Order of the two measurements must not change the verdict."""
    assert _disagreement(1.0, 1.0) == 1.0
    assert _disagreement(1.10, 1.00) == pytest.approx(1.10)
    assert _disagreement(1.00, 1.10) == pytest.approx(1.10)


def test_the_drift_statistic_reads_the_worse_of_the_two_sides():
    """Either side drifting makes the measurement unreadable.

    The ratio is a quotient of the two minima, so instability on the DSL
    side corrupts it exactly as much as instability on the reference
    side -- and the reference side is the one the original docstring
    described watching.
    """
    steady = [1.0, 1.0, 1.0]
    assert _drift(steady, steady) == 1.0
    assert _drift([1.0, 2.0, 1.0], steady) == pytest.approx(2.0)
    assert _drift(steady, [1.0, 2.0, 1.0]) == pytest.approx(2.0)


def test_an_abstention_is_loud_enough_to_find_afterwards():
    """A silent skip on a guard is indistinguishable from a boring one.

    A skipped optional-dependency test and a skipped performance guard
    mean very different things -- the second says the published figure
    went unchecked on that run.  So abstaining warns as well as skips,
    which puts it in the warnings summary of a plain run rather than
    only behind `-rs`.

    ⚠ `pytest.raises(Exception)` does NOT catch this.  `Skipped` derives
    from `BaseException` only, so the first draft of this test let the
    skip escape and BECAME a skip -- passing its own review by silently
    not running, which is the exact failure this file is about.  Catch
    `pytest.skip.Exception` by name.
    """
    with pytest.warns(RuntimeWarning, match='NOT checked on this run'):
        with pytest.raises(pytest.skip.Exception) as exc:
            _abstain(1.4, 3.0)
    ## Both numbers survive into the reason: the disagreement is why the
    ## reading was dropped, the drift says whether the machine was why.
    assert '1.400x' in str(exc.value), str(exc.value)
    assert '3.00x' in str(exc.value), str(exc.value)


def _canned(monkeypatch, ratios):
    """Drive `_confirmed_ratio` off a scripted list of estimates.

    The retry loop is live-path logic and the clock cannot be asked to
    produce a disagreement on demand, so the measurement is replaced and
    the DECISIONS are what get tested.  `calls` records how many
    measurements were actually taken, which is the cost half of the
    design.
    """
    calls = []

    def fake(a, b, pairs=None):
        r = ratios[len(calls)]
        calls.append(r)
        ## Timings are irrelevant here; only their drift is reported.
        return r, ('wave_a', 7), ('wave_b', 7), ([1.0, 1.0], [1.0, 1.0])

    monkeypatch.setattr('pycircuit.circuit.tests.test_perf_guards'
                        '._interleaved_ratio', fake)
    return calls


def test_a_reproducible_pair_costs_exactly_two_measurements(monkeypatch):
    """The quiet-machine path must not pay for the busy-machine one."""
    calls = _canned(monkeypatch, [1.01, 1.02, 9.99, 9.99])
    ratio, _, _ = _confirmed_ratio((R, C), (R, C), MAX_OVERHEAD, want='below')
    assert len(calls) == 2, calls
    ## 'below' keeps the less damning of the pair.
    assert ratio == pytest.approx(1.01)


def test_a_disagreeing_pair_is_retried_rather_than_abstained_on(monkeypatch):
    """One bad measurement must not cost the run its only check.

    Measured under load 17, consecutive 12-sample estimates disagreed
    about half the time -- so abstaining on the first disagreement would
    abstain on half of all runs, while a machine that is busy on average
    is still quiet in patches.
    """
    calls = _canned(monkeypatch, [1.01, 1.40, 1.03, 1.02])
    ratio, _, _ = _confirmed_ratio((R, C), (R, C), MAX_OVERHEAD, want='below')
    assert len(calls) == 4, calls
    assert ratio == pytest.approx(1.02)


def test_agreement_is_judged_on_the_two_most_recent_only(monkeypatch):
    """The harness must not shop the set for the friendliest pair.

    1.01 and 1.02 agree, but they are not consecutive: 1.40 sits between
    them.  Taking the closest two from a growing set would let a run pick
    whichever pair suited it, and two readings agree most easily when
    both are wrong the same way.
    """
    calls = _canned(monkeypatch, [1.01, 1.40, 1.02, 1.99])
    with pytest.warns(RuntimeWarning):
        with pytest.raises(pytest.skip.Exception):
            _confirmed_ratio((R, C), (R, C), MAX_OVERHEAD, want='below')
    assert len(calls) == MAX_MEASUREMENTS, calls


def test_a_machine_that_never_settles_abstains_and_says_so(monkeypatch):
    """Bounded, loud, and not a pass."""
    calls = _canned(monkeypatch, [1.00, 1.50, 1.00, 1.50])
    with pytest.warns(RuntimeWarning, match='NOT checked on this run'):
        with pytest.raises(pytest.skip.Exception):
            _confirmed_ratio((R, C), (R, C), MAX_OVERHEAD, want='below')
    assert len(calls) == MAX_MEASUREMENTS, calls
