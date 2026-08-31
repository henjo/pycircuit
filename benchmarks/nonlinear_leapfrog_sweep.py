"""Fully symbolic nonlinear analysis of the leapfrog, 127 unknowns:
amplitude swept, orders to U^13, accuracy against **transient simulation** -- via
**two-tone IM3** rather than HD3.

Stage 14 attached a nonlinearity *inside* the perturbation solve, as a callable.
Nothing in the circuit was nonlinear, so a transient run would have been perfectly
linear and produced no harmonics -- there would have been nothing to compare against.

Here the nonlinearity is a real element, a `BSource` cubic ``i = g3*v**3`` across one
amplifier's input node, so **the transient engine and the perturbation analysis see
the same circuit**.  Two independent paths, no shared code on the measurement side --
a different solver, a different representation, a different failure mode -- which is
the property the existing distortion-vs-transient gate was built around.

One simplification makes this exact rather than approximate: ``d/dv (g3 v^3) = 0`` at
``v = 0``, so the cubic contributes nothing to the linearised matrix.  The matrix the
perturbation solves against is therefore *exactly* the linear leapfrog's, and the
nonlinearity enters only through ``f(x)``.  No separate linearisation is needed and
none can drift.

WHY IM3 AND NOT HD3
-------------------
HD3 cannot be measured on this circuit by any transient, at any cost.  The frequency
that makes the distortion appear at all is the frequency that kills the measurement:
the drive has to sit where the amplifier's loop gain has fallen (100 kHz), and that
is 2.6 octaves into this 5th-order lowpass's stopband.  The fundamental is attenuated
106x from the nonlinear node to the output; the **third harmonic, at 300 kHz, is
attenuated 160 000x**.  Measured, at a 1 V drive:

    3rd harmonic at s0_e1   6.21e-06 V        HD3 = 9.84e-04
    3rd harmonic at output  3.89e-11 V        HD3 = 5.83e-07

An integrator with ``reltol`` r leaves roughly ``r*|signal|`` of error behind, and the
output signal is 6.7e-05 V, so even ``reltol = 1e-8`` leaves ~7e-13 V against a 3.89e-11 V
harmonic.  The quantity is not small because the analysis is imprecise; it is small
because five lowpass stages removed it.

Two tones at 100 and 110 kHz put the IM3 product ``2*f1 - f2`` at **90 kHz -- beside
the fundamentals**, where it takes essentially the same filter loss (a 5th-order
rolloff across 90->100 kHz is a factor of 1.7, against 160 000x for the harmonic).
Same circuit, same nonlinearity, same drive level; only the observable changes:

    IM3 at output, 1 V drive  1.08e-08 V      IM3/fund = 1.61e-04

277x larger than HD3 at the same node and drive, which lifts it clear of the floor.

WHY 100 AND 110 kHz SPECIFICALLY
--------------------------------
Tone spacing sets the run length, because the FFT has to resolve the tones from the
product.  100/101 kHz would need a >= 1 ms window for 1 kHz resolution.  With 10 kHz
spacing everything is a multiple of ``FB = 10 kHz``, so a 100 us window puts the three
signals of interest exactly on bins 9, 10 and 11 -- no leakage, no windowing, and a
tenth of the integration.

RUN IT SINGLE-THREADED.  Measured 2026-08-31, and it is the largest single cost
factor in this file -- larger than the integrator choice, larger than the settle:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/nonlinear_leapfrog_sweep.py

numpy threads the dense solve across every core it can find, and at n = 139 that
is pure overhead.  Three interleaved pairs on the linear circuit, identical 292
output points each time:

    single-threaded   11.1 s   10.3 s    8.2 s
    default          159.5 s  145.1 s  164.9 s

Minima 8.2 s against 145.1 s -- **17.7x**, per-pair 14.4x/14.1x/20.1x.  With the
cubic live the ratio holds: **0.73 s per simulated microsecond single-threaded
against 10.35 s, 14.2x**.  So the "hours per amplitude" this measurement has
been parked behind since July was substantially an environment variable: at
0.73 s/us a full 5-tau run is ~14 min per amplitude, and a seeded 2-tau run is
~4 min.

Run:  PYTHONPATH=<repo>:<repo>/benchmarks python3 benchmarks/nonlinear_leapfrog_sweep.py
      (add --settle-convergence for the settle/tolerance calibration)
"""

import sys
import time

import numpy as np
import sympy

from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.benchmark_circuits import (_UA741_NODE_NAMES,
                                                 build_leapfrog_network)
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.distortion import (GradedSpectrum, GradedVector,
                                          graded_response_mimo)
from pycircuit.circuit.integrator import TrapezoidalIntegrator
from pycircuit.circuit.transient import Transient, resample_uniform


## The bin base.  Every frequency here is an integer multiple of it, so one 1/FB
## window resolves all of them exactly.
FB = 10.0e3
F1, F2 = 10 * FB, 11 * FB       # 100 kHz, 110 kHz
BIN_IM3, BIN_F1, BIN_F2 = 9, 10, 11        # 2*f1 - f2 = 90 kHz, then the tones
## Matched to `order_convergence.py` so the two paths describe the same circuit:
## same node, same coefficient, same drive levels.
G3 = 5.0e-2
STAGES = 5
NODE_NAME = 's0_e1'
## Largest drive FIRST.  IM3 grows faster than the fundamental, so 3 V has the
## widest margin over the integrator floor (778x at reltol 1e-6, against 162x at
## 1 V) -- and each transient costs hours, so if only one row completes it should
## be the one with the most headroom.
AMPLITUDES = (3.0, 1.0)
ORDERS = (3, 5, 7, 9, 11, 13)

## SETTLING, NOT STEPPING, IS THE COST DRIVER -- and the figure it is set from is a
## property of the repaired circuit, not a round number.
##
## The slowest pole is Re = -4.815e+03, so tau = 208 us.  That is NOT the 10 us
## integrator time constant an earlier version of this file assumed: the slowest mode
## is the Q ~ 16.8 resonance at 25.7 kHz that the topology repair disclosed and could
## not remove (omega0/(2Q) = 2*pi*25700/(2*16.8) = 4806, which is that pole).  Settling
## to 5 tau is ~1.04 ms, an order of magnitude more simulated time than the 200 us
## previously used -- and at ~9 ns steps that is hours per amplitude.
SETTLE = 5 * 208e-6
MEASURE = 1.0 / FB
## vabstol comes back DOWN here, well below the library default of 1e-6 (Spectre's
## value, right for general use).  This measurement needs the integrator's residual
## error under the ~1e-8 V IM3 amplitude at the output, and the default leaves
## etol = TRTOL*(reltol*|x| + 1e-6) ~ 7e-6 -- larger than the 6.7e-05 V signal's own
## precision requirement.  Measured on a fixed 2.5 us window, |vout| error against a
## tight reference: defaults +2.03%, vabstol 1e-9 +0.45%, vabstol 1e-9 + reltol 1e-6
## -0.02%.  The last is what IM3 needs and what it costs.
## THE INTEGRATOR IS AN EXPLICIT CHOICE, and it is worth 10x.  This harness never
## passed one, so it inherited `EulerIntegrator()` -- 1st order, LTE ~ h^2.  On the
## repaired circuit at these tolerances, over a fixed 2.5 us window:
##
##     euler (the inherited default)   2896 steps   482.5 s
##     trapezoidal                      288 steps    53.6 s   <- 10x fewer
##     gear2 (ywr)                      347 steps    60.7 s
##     gear2 (classic)                  376 steps    67.3 s
##
## The three 2nd-order methods agree within 30% of each other; Trapezoidal wins by
## roughly the ratio of its error constant to Gear-2's (1/12 against 2/9).  The gap to
## Euler is this large *because* the cost here is resolving a lightly-damped Q ~ 16.8
## resonance -- a smooth oscillation, which is precisely what a 2nd-order method
## integrates well and backward Euler does not.
##
## This only became a safe choice once the LTE work landed: before it,
## `Gear2Integrator`'s estimate returned ~1e-15 of the true error, so it ran at
## `max_step` with zero rejections and produced numbers that looked fine and controlled
## nothing.  See `doc/transient_repair_plan.md`.
TRAN_OPTS = {'vabstol': 1e-9, 'reltol': 1e-6,
             'integrator': TrapezoidalIntegrator()}
## Samples per 1/F1 period.  This is NOT just an output grid -- `timestep` also sets
## `max_step`, which here is a CORRECTNESS knob, not only a cost knob.
##
## `IntegralController.evaluate_step` accepts unconditionally when `is_first_step`
## (stepcontroller.py:26), and `transient.py:314` re-arms that flag after every
## breakpoint.  `Sin.next_event` (func.py:30) fires every quarter period, so each
## tone contributes a periodic, drive-synchronous step that gets NO error check and
## runs at the full `max_step`.  Capping `max_step` is what bounds the damage.
##
## 512 costs nothing at the tolerances in TRAN_OPTS: the controller settles at
## h ~ 4.4 ns there, well under the 19.5 ns cap, so the cap never binds on an
## accuracy-controlled step and only limits the unchecked breakpoint ones.  It is
## not free at loose tolerance -- with the library defaults the controller wants
## ~96 ns and the cap then forces >6000 steps -- so a cheap structural check should
## lower it rather than inherit this value.
POINTS = 512


def build_transient_leapfrog(amplitude):
    """The same topology as `leapfrog_5th_order`, numeric, driven by TWO tones.

    The *build* is separate from the fixture's because the fixture uses the symbolic
    toolkit and forces a symbolic `VS`, while a transient needs numeric elements and
    a time-domain source.  But nothing about the circuit is replicated: the devices
    come from `_stamp_ua741_devices` and the coupling network from
    `wire_leapfrog_couplings`, both shared with the fixture.

    An earlier version of this docstring claimed the circuit "cannot drift from the
    fixture's" on the strength of sharing `_stamp_ua741_devices` alone.  That was
    true of the transistors and false of everything connecting them, and the
    difference mattered: when the sign error that put two poles in the right half
    plane was fixed in the fixture, the copy here kept the broken wiring.

    The two tones are two sources in series, so the input node carries their sum.
    """
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        ## The ENTIRE circuit, from the same builder the fixture uses -- including
        ## the source.  The builder dispatches on the toolkit: numeric here, so it
        ## injects one `VSin` per tone in series rather than the symbolic `VS` the
        ## fixture gets.  There is nothing left for this function to get wrong.
        node, amp, _names = build_leapfrog_network(
            cir, stages=STAGES, tones=((amplitude, F1), (amplitude, F2)))

        ## THE nonlinearity, and the only nonlinear element in the circuit: a
        ## cubic conductance to ground on stage 0's input-pair emitter.
        nl = amp[0][NODE_NAME.split('_')[-1]]
        cir['nl'] = BSource(nl, gnd, nl, gnd, i_func=lambda v: G3 * v ** 3)
        return cir, amp
    finally:
        circuit_module.default_toolkit = saved


def linear_steady_state_x0(cir, source_names, refnode=gnd):
    """The LINEAR circuit's state at t=0, by AC superposition -- the settling seed.

    SETTLING, NOT STEPPING, IS THIS MEASUREMENT'S COST.  The slowest pole is the
    Q ~ 16.8 resonance the topology repair disclosed, tau = 208 us, so a
    from-rest run must integrate ~5 tau = 1.04 ms of physics before the first
    microsecond of it is worth measuring -- against a 100 us measurement window.
    Eleven twelfths of the run exists to throw away a transient whose value is
    known in closed form.

    It is known in closed form because **this circuit is linear apart from the
    cubic**: the amplifiers are small-signal stamps (`add_small_signal_bjt`) and
    the only nonlinear element is `BSource(i = g3 v^3)`, whose derivative is
    ZERO at v = 0.  So the matrix an AC analysis assembles at the origin is
    *exactly* the linear leapfrog's -- the same property that lets the
    perturbation series linearise once and never drift.  The cubic contributes
    ~1e-4 of the response, and that part still has to settle; everything else
    can simply be evaluated.

    THE CONVENTION, which is the part worth checking rather than assuming.
    `func.Sin` with the harness's parameters (`vo=td=theta=0`) is

        v(t) = va * sin(omega t + phase)
             = Im{ (va e^{j phase}) e^{j omega t} }

    and `VS.u(analysis='ac')` builds its stimulus as `vac * e^{j phase}`.  Those
    are the SAME complex amplitude, so the AC solution X at that frequency
    satisfies x_ss(t) = Im{X e^{j omega t}}, and therefore

        x_ss(0) = Im{X}.

    ⚠ Note `phase` is ONE parameter serving both roles here -- `VS` declares it
    as the AC phase and `VSin` re-declares it as the sine phase.  That is a
    latent trap in general and is exactly right here, because the two
    conventions agree; it is why only `vac` needs setting below.

    Superposition over the tones is legitimate for the same reason: the system
    the seed describes is linear, so the two tones' responses add.

    ⚠ Every `VS` in this tree defaults to **`vac = 1`**, not 0.  A source left
    alone would therefore inject a spurious unit AC stimulus into every solve,
    so this zeroes the AC amplitude of every source in the circuit and restores
    it afterwards -- the function must not leave the circuit changed.

    Validated on an RC whose steady state is analytic: the seed reproduced
    `v_out(0)` to all printed digits, and over three periods the seeded run
    deviated from the analytic steady state by **7.0e-07** against **4.5e-01**
    from a zero start.

    Args:
        cir: the circuit, already built.
        source_names: element names of the tone sources, e.g. ``('vs0', 'vs1')``.
        refnode: reference node, as passed to the transient.

    Returns:
        A full-length ``x`` vector in the circuit's own ordering, with the
        reference row present -- what ``Transient.solve(x0=...)`` expects.
    """
    from pycircuit.circuit.analysis_ss import AC

    ## Zero every AC stimulus first, so superposition drives exactly one tone
    ## per solve and no defaulted `vac=1` leaks in.
    saved = {}
    for name, el in cir.elements.items():
        if 'vac' in el.iparv:
            saved[name] = el.iparv.vac
            el.iparv.vac = 0.0

    ## Linearise at the ORIGIN rather than at whatever a DC solve returns.  At
    ## v = 0 the cubic's derivative is exactly zero, which is the property this
    ## whole construction rests on; leaving `dcx` unset would run a DC analysis
    ## and linearise whereever it happened to land.
    dcx = np.zeros(cir.n)

    x0 = np.zeros(cir.n)
    try:
        for name in source_names:
            el = cir.elements[name]
            va, freq = float(el.iparv.va), float(el.iparv.freq)
            el.iparv.vac = va
            try:
                res = AC(cir, toolkit=numeric, dcx=dcx).solve(
                    freqs=np.array([freq]), refnode=refnode)
                X = np.asarray(res.x)
                X = X[:, 0] if X.ndim == 2 else X
                x0 = x0 + np.imag(X)
            finally:
                el.iparv.vac = 0.0
    finally:
        for name, value in saved.items():
            cir.elements[name].iparv.vac = value

    return x0


def transient_im3(amplitude, opts=None, settle=SETTLE, measure=MEASURE,
                  seed=True):
    """IM3 at the output and at the nonlinear node, from a transient plus an FFT.

    ``seed`` starts the run from the linear circuit's steady state
    (`linear_steady_state_x0`) instead of from rest, which is what makes a
    short ``settle`` legitimate: the transient that ``settle`` exists to
    discard is, to within the cubic's ~1e-4 contribution, already gone at t=0.

    ⚠ ``seed`` does NOT license ``settle=0`` on its own.  The seed removes the
    LINEAR settling; what remains is the nonlinear part reaching its own
    periodic state, and only a measurement can say how long that takes.  Use
    `settle_convergence` to take it rather than assuming, and keep the settle
    the measurement supports.

    Returns a dict per node of ``(im3_ratio, fund, im3)`` plus cost.
    """
    cir, amp = build_transient_leapfrog(amplitude)
    probes = {'out': amp[STAGES - 1]['out'], NODE_NAME: amp[0]['e1']}
    x0 = linear_steady_state_x0(cir, ('vs0', 'vs1')) if seed else None
    t0 = time.time()
    res = Transient(cir, toolkit=numeric,
                    **(TRAN_OPTS if opts is None else opts)).solve(
        refnode=gnd, tend=settle + measure, timestep=1.0 / (F1 * POINTS),
        x0=x0)
    secs = time.time() - t0

    out = {}
    npt = int(round(measure * F1 * POINTS))
    for label, nd in probes.items():
        wave = res.v(nd)
        t = np.asarray(wave.x[0])
        v = np.asarray(wave.y)
        ## Resample the measured window onto an exact grid; the solver's own steps
        ## are not uniform, and the bins only land if the grid is.
        grid = np.linspace(t[-1] - measure, t[-1], npt, endpoint=False)
        ## QUADRATIC, matching the integrator's order -- `resample_uniform`
        ## (transient.py, stage 10.2) rather than `np.interp`.  That function
        ## was built to fix exactly this line and names this file in its
        ## docstring; it went unadopted because its grid contract could not
        ## express a TRAILING window with `endpoint=False`, which is what keeps
        ## the tones on exact bins.  It takes an explicit `grid=` as of
        ## 2026-08-31.
        ##
        ## ⚠ Measured, on a synthetic two-tone at this file's own IM3 ratio:
        ## how much this matters depends ENTIRELY on the step-size spread of
        ## the adaptive grid, not on the mean step.
        ##
        ##     spread    linear IM3 error    quadratic
        ##      10x           +0.16%           +0.00%
        ##     100x          +39.78%          -10.51%
        ##    1000x        +97055%          +33767%
        ##
        ## At 1000x neither interpolant helps: the grid is not resolving the
        ## signal and the fix is a smaller `max_step`, not a cleverer resample
        ## (`resample_uniform`'s own caveat).  `POINTS` caps `max_step` here,
        ## which is what keeps this run in the benign regime -- so that cap is
        ## load-bearing for ACCURACY, not only for the unchecked breakpoint
        ## steps its own comment describes.
        _g, v_uniform = resample_uniform(t, v, grid=grid)
        spec = np.fft.rfft(v_uniform) / npt
        ## The window is exactly 1/FB, so bin k IS k*FB and the constants can be
        ## used directly.  Asserted rather than assumed: a window that is not a
        ## whole number of 1/FB puts the tones between bins and the IM3 reading
        ## becomes leakage from the fundamental, which would look like distortion.
        assert abs(measure * FB - round(measure * FB)) < 1e-12, \
            'measurement window must be a whole number of 1/FB'
        fund, im3 = abs(spec[BIN_F1]), abs(spec[BIN_IM3])
        out[label] = (im3 / fund if fund else float('nan'), fund, im3)
    return out, secs, len(np.asarray(res.v(probes['out']).x[0]))


def settle_convergence(amplitude=1.0, settles=(0.0, 208e-6, 2 * 208e-6,
                                                5 * 208e-6), seed=True):
    """How much settle the run actually needs -- measured, not assumed.

    THE SEED DOES NOT BY ITSELF LICENSE A SHORTER SETTLE, and this is the
    measurement that says what it does license.  `linear_steady_state_x0`
    removes the LINEAR settling exactly (verified: the seeded run starts on the
    AC-predicted trajectory with e(t=0) = 0 and thereafter deviates only by the
    integrator's own h^2 error -- 2.455e-05 at the harness tolerance against
    2.670e-07 at 1000x tighter, a 92x fall for a 9.6x step reduction, which is
    the trapezoidal rule's order and nothing else).

    What the seed does NOT remove is the CUBIC's own approach to periodic
    steady state.  That part is ~1e-4 of the response and its settling is a
    property of the nonlinearity, not of the linear poles, so the only honest
    way to choose a settle is to watch IM3 stop moving as the settle grows.

    Prints IM3 at both probes against settle, so the shortest settle whose IM3
    has converged can be read off rather than guessed.  Run it before quoting
    any T4 number from a shortened run.

    MEASURED 2026-08-31, amplitude 1.0, seeded, single-threaded:

        settle    tau    IM3/fund out    IM3/fund nl     secs
        0         0.00   7.630947e-04    2.755349e-03     59.9
        208 us    1.00   5.488567e-04    2.808578e-03    106.6
        416 us    2.00   7.157515e-04    2.811437e-03    222.3
        1040 us   5.00   7.047518e-04    2.811553e-03    671.4

    **The nl probe settles; the out probe does not converge at all.**  The nl
    node moves monotonically to within 0.004% of the 5 tau reference by 2 tau
    (0.1% by 1 tau), which says the seeded circuit reaches periodic steady
    state long before the settle this harness pays for.  The out node reads
    +8.28%, -22.12%, +1.56% against 5 tau -- not monotone, so not a decaying
    transient, and both numbers come from the SAME runs, the same window and
    the same FFT, so the circuit cannot be settled for one probe and not the
    other.

    ⚠ **THE OUT PROBE IS INTEGRATOR-LIMITED AT THIS HARNESS'S TOLERANCES, and
    that is a defect in the tolerance choice, not in the settle.**  Holding the
    settle fixed at 2 tau and tightening `reltol` 1e-6 -> 3e-8 and `vabstol`
    1e-9 -> 1e-11:

        harness tol   IM3/fund out 7.157515e-04   IM3 out 3.378931e-08   nl 2.811437e-03
        30x tighter   IM3/fund out 1.760965e-04   IM3 out 8.334364e-09   nl 2.840985e-03

    The output number moves by **4.07x** while the nl number moves 1.05%.  It
    is still moving at 30x, so the true output IM3 is NOT bracketed -- it is
    somewhere at or below 8.3e-09 V and this harness has not yet measured it.

    ⚠ **The mechanism is a metric chosen on the wrong signal**, and it is
    recorded in this file's own `TRAN_OPTS` comment: the tolerance was
    calibrated by driving `|vout|` error to -0.02% against a tight reference.
    But `IM3/fund` at the output is **1.6e-04**, so a 2e-4 relative error on
    the fundamental is the same order as the entire quantity being measured.
    Accuracy of a large signal does not bound accuracy of a component four
    decades beneath it.  The 162x margin this module's header claims over the
    integrator floor is not supported at these tolerances.

    ⚠ So: **the seeding and the settle are solved; the output tolerance is
    not.**  Before T4 quotes an output IM3, run this at a fixed settle over a
    tolerance ladder and find where the number stops moving.  The nl-node
    figure is already trustworthy and can be quoted now.

    THE LADDER WAS RUN, 2026-08-31, against the perturbation series (U^13:
    out 1.816614e-04, nl 2.844758e-03), amplitude 1.0, settle 2 tau:

        tolerance        IM3/f @out     vs pert    IM3/f @nl      vs pert
        harness          7.157515e-04   +294.00%   2.811437e-03    -1.17%
        30x              1.760965e-04     -3.06%   2.840985e-03    -0.13%
        300x             1.696682e-04     -6.60%   2.843814e-03    -0.03%

    **The nl node converges monotonically onto the series** -- T4 passes there,
    two independent methods to 3 parts in 10^4.

    ⚠ **The output CROSSES the series value and keeps going** (+294 -> -3.06 ->
    -6.60), so it is not simply integrator-limited converging onto the oracle,
    and the tolerance ladder alone does not explain it.  All three rungs held
    settle at 2 tau, leaving two variables unseparated: the cubic's own
    settling as seen through five filter stages (the seed removes only the
    LINEAR settling), and an error `reltol` does not control -- the FFT
    resamples the solver's non-uniform steps onto a uniform grid by LINEAR
    interpolation, on a quantity reached through a five-stage cancellation.
    The next measurement is 300x at 2 tau against 5 tau.

    ⚠ Run this single-threaded.  numpy threads a 139-unknown dense solve across
    every core it can find, which is overhead at this size and antisocial on a
    shared machine:

        OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
            python benchmarks/nonlinear_leapfrog_sweep.py --settle-convergence
    """
    print('settle sweep at amplitude %.3g, seed=%s' % (amplitude, seed))
    print('%-12s %-9s %-13s %-13s %-13s %s'
          % ('settle', 'tau', 'IM3/fund out', 'IM3 out', 'IM3/fund nl', 'secs'))
    rows = []
    for settle in settles:
        out, secs, npts = transient_im3(amplitude, settle=settle, seed=seed)
        r_out, _f_out, im3_out = out['out']
        r_nl, _f_nl, _im3_nl = out[NODE_NAME]
        rows.append((settle, r_out, im3_out, r_nl, secs))
        print('%-12.4g %-9.2f %-13.6e %-13.6e %-13.6e %.1f'
              % (settle, settle / 208e-6, r_out, im3_out, r_nl, secs),
              flush=True)
    ## The useful reading is the SPREAD against the longest settle, which is
    ## what says whether a short run is quoting a settled number.
    if len(rows) > 1:
        ref = rows[-1][1]
        print('\nrelative to the longest settle (%.4g s):' % rows[-1][0])
        for settle, r_out, _im3, _r_nl, _secs in rows[:-1]:
            print('  settle %-10.4g IM3/fund differs by %+.2f%%'
                  % (settle, 100.0 * (r_out - ref) / ref))
    return rows


def perturbation_im3(system, node_row, out_row, drive_row, amplitude, order):
    """IM3 from the two-tone graded perturbation series truncated at ``U^order``.

    The solver here is numeric per amplitude, deliberately: this script's job is to
    check the series against an independent time-domain solver, not to demonstrate
    the symbolic representation -- `order_convergence.py` is the one that keeps
    ``Adrv`` and ``kk`` symbolic and substitutes last.

    ``A`` is affine in ``s``, so ``A0 + s*C`` is built once numerically instead of
    substituting into 16 129 sympy entries per solve.  Two tones need far more solves
    than one tone does, and that shortcut is what makes this affordable; it is
    checked against the recorded single-tone HD3 numbers by `main`.
    """
    n, s_sym = system.dim, system.s
    A = sympy.Matrix(system.A)
    env = {s: complex(system.params[s]) for s in system.A.free_symbols
           if s is not s_sym}
    A0 = np.array([[complex(A[i, j].subs(s_sym, 0).subs(env)) for j in range(n)]
                   for i in range(n)], dtype=complex)
    Cm = np.array([[complex(sympy.diff(A[i, j], s_sym).subs(env))
                    for j in range(n)] for i in range(n)], dtype=complex)

    def solve(s, rhs):
        return list(np.linalg.solve(A0 + s * Cm, np.asarray(rhs, dtype=complex)))

    def f(x):
        outv = [GradedSpectrum() for _ in range(n)]
        cube = ((x[node_row] * x[node_row]).truncated(order)
                * x[node_row]).truncated(order)
        outv[node_row] = cube.scaled(G3)
        return GradedVector(outv)

    vec = GradedVector([GradedSpectrum() for _ in range(n)])
    for idx in ((1, 0), (0, 1)):
        for k, v in GradedSpectrum.from_phasor(idx, 1,
                                               complex(amplitude)).terms.items():
            vec.components[drive_row].terms[k] = \
                vec.components[drive_row].terms.get(k, 0) + v

    t0 = time.time()
    g = graded_response_mimo(solve, vec, f, (2 * np.pi * F1, 2 * np.pi * F2),
                             max_power=order)
    secs = time.time() - t0
    out = {}
    for label, row in (('out', out_row), (NODE_NAME, node_row)):
        fund = abs(g[row].phasor((1, 0)))
        im3 = abs(g[row].phasor((2, -1)))
        out[label] = (im3 / fund if fund else float('nan'), fund, im3)
    return out, secs


def main():
    t0 = time.time()
    system = bc.leapfrog_5th_order()
    names = ['in'] + ['s%d_%s' % (k, b) for k in range(STAGES)
                      for b in _UA741_NODE_NAMES]
    node_row = names.index(NODE_NAME)
    out_row = system.out_index
    drive_row = [i for i in range(system.dim) if system.b[i] != 0][0]
    print('leapfrog: %d unknowns, built in %.1f s' % (system.dim,
                                                      time.time() - t0))
    print('nonlinearity: BSource i = %g*v^3 at %s (row %d); drive row %d,'
          ' output row %d' % (G3, NODE_NAME, node_row, drive_row, out_row))
    print('two tones %g/%g kHz; IM3 at 2f1-f2 = %g kHz (bins %d/%d/%d of a'
          ' %g us window)' % (F1 / 1e3, F2 / 1e3, (2 * F1 - F2) / 1e3,
                              BIN_IM3, BIN_F1, BIN_F2, MEASURE * 1e6))
    print('transient options: %s (library defaults are vabstol 1e-6/reltol 1e-4;'
          ' see TRAN_OPTS)' % TRAN_OPTS)
    print('the cubic has zero slope at v=0, so the linearised matrix IS the'
          ' linear leapfrog\'s')
    print()

    ## The perturbation side, first: it is cheap, and it says what the transient is
    ## being asked to resolve before any integration is paid for.
    print('== perturbation IM3 (independent of the transient) ==')
    print('%-7s %-6s %-11s %-11s %-11s %-11s %-11s'
          % ('amp', 'order', 'fund@out', 'IM3@out', 'IM3/f@out', 'IM3@%s' % 's0_e1',
             'IM3/f@e1'))
    print('-' * 76)
    pert = {}
    for amp in AMPLITUDES:
        for order in ORDERS:
            try:
                r, secs = perturbation_im3(system, node_row, out_row, drive_row,
                                           amp, order)
            except Exception as exc:
                print('%-7.3g U^%-4d FAILED: %s: %s'
                      % (amp, order, type(exc).__name__, exc))
                break
            pert[(amp, order)] = r
            ro, re = r['out'], r[NODE_NAME]
            print('%-7.3g U^%-4d %.4e %.4e %.4e %.4e %.4e  (%.0f s)'
                  % (amp, order, ro[1], ro[2], ro[0], re[2], re[0], secs))
            sys.stdout.flush()
        print()

    ## Self-convergence of the series, which needs no external oracle.
    print('== does the series settle?  d(U^k) = |IM3(U^k) - IM3(U^k-2)| /'
          ' |IM3(U^k)| ==')
    head = ('%-7s %-11s' % ('amp', 'IM3/f@out')
            + ''.join('%-11s' % ('d(U^%d)' % o) for o in ORDERS[1:])
            + ' agrees from')
    print(head)
    print('-' * len(head))
    for amp in AMPLITUDES:
        have = [o for o in ORDERS if (amp, o) in pert]
        if len(have) < 2:
            continue
        row = '%-7.3g %-11.4e' % (amp, pert[(amp, have[-1])]['out'][0])
        agrees = None
        for prev, order in zip(have, have[1:]):
            a = pert[(amp, prev)]['out'][2]
            c = pert[(amp, order)]['out'][2]
            rel = abs(c - a) / abs(c) if abs(c) else float('nan')
            row += '%-11s' % ('%.1e' % rel)
            if agrees is None and rel < 1e-6:
                agrees = prev
        row += ' %s' % ('U^%d' % agrees if agrees else
                        'not by U^%d' % have[-1])
        print(row)
    print()

    ## The transient, last, and only where the perturbation says IM3 clears the
    ## integrator's floor.  Running it where it cannot would produce a number, and
    ## the number would be integration error.
    print('== transient IM3, and the comparison ==')
    print('estimated cost at these tolerances: ~4-5 h per amplitude for a'
          ' %g us run.' % ((SETTLE + MEASURE) * 1e6))
    print('Each amplitude prints as it completes, so a kill leaves partial results.')
    print()
    header = ('%-7s %-11s %-11s %-9s %-8s %-7s'
              % ('amp', 'IM3/f tran', 'IM3/f pert', 'rel diff', 'steps', 'hours'))
    print(header)
    print('-' * len(header))
    for amp in AMPLITUDES:
        best = max(o for o in ORDERS if (amp, o) in pert) if any(
            (amp, o) in pert for o in ORDERS) else None
        if best is None:
            continue
        try:
            tr, secs, nsteps = transient_im3(amp)
        except Exception as exc:
            print('%-7.3g TRANSIENT FAILED: %s: %s'
                  % (amp, type(exc).__name__, exc))
            break
        for label in ('out', NODE_NAME):
            p = pert[(amp, best)][label][0]
            t = tr[label][0]
            rel = abs(t - p) / abs(p) if p else float('nan')
            print('%-7.3g %-11.4e %-11.4e %-9.2e %-8d %-7.2f  @%s (U^%d)'
                  % (amp, t, p, rel, nsteps, secs / 3600.0, label, best))
        sys.stdout.flush()

    print()
    print('The transient is an independent solver; nothing on that path shares')
    print('code with the analysis under test.  A transient IM3 that matches to a')
    print('few percent is a real cross-check; one that matches to 1e-6 would mean')
    print('the two paths are not actually independent and should be suspected.')
    return 0


if __name__ == '__main__':
    ## `--settle-convergence` takes the measurement that says how short a
    ## seeded run's settle may be.  Kept as a separate entry point because it
    ## is a calibration of the harness, not a row of its results.
    if '--settle-convergence' in sys.argv:
        settle_convergence()
        sys.exit(0)
    sys.exit(main())
