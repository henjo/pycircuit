# JAX transient — limitations & refinement roadmap

Status (2026-07): the JAX transient (`pycircuit/circuit/jaxtransient.py`) is now
**algorithmically aligned** with the CPU transient (`transient.py`) on the
LTE / step-control math — integration-method dispatch, method order, LTE-band
step rejection, and the YWR DAE LTE (default) all match, verified in
`tests/test_jaxtransient.py`.

What remains is the CPU transient's **convergence-robustness machinery**. The
JAX engine is currently a "happy-path" solver: on par for linear / mildly
non-linear circuits, but the CPU path is still the reliable one for stiff,
strongly non-linear circuits. Gaps below, worst first.

## Robustness / convergence (highest value)

1. **Device limiting.** CPU applies `cir.limit(...)` (pn-junction / `pnjlim`)
   inside every Newton iteration; JAX only has a uniform `max_dv=0.5` voltage
   clamp in `newton_inner_loop`. This is the biggest convergence gap for
   diodes/BJTs. *Port target #1.*
2. **Breakpoint order reset.** JAX clamps `dt` onto breakpoints but does not
   drop to 1st order + reset integrator history after crossing a discontinuity
   (CPU: `was_break_step -> _is_first_step`). Risk: ringing/overshoot on sharp
   edges (VPulse/PWL). *Port target #2.*
3. **Nonlinear continuation.** ~~No gmin-stepping / source-stepping (CPU
   `nrsolver.py` strategies). Hard operating points that need homotopy
   fail.~~ **DONE** — DC got the full adaptive ladder
   (`dc_with_continuation`: junction-gmin -> gshunt -> Psi-tc, P18/P25),
   and the transient got the failed-time-point rescue on 2026-09-01
   (`transient_point_rescue`, gated at the dt floor, on the plain Newton
   and both PCNR views; `continuation` Parameter, default off because a
   traced graph carries the branch whether or not it runs).  SOURCE
   stepping stays absent on purpose: scaling u(t) mid-transient would
   scale the integrator's companion history with it.
4. **Jacobian scaling / equilibration.** No scaler strategies (CPU
   `_get_scaler`: RowMax / L2 / Sinkhorn). Worse conditioning on wide
   dynamic-range matrices.
5. **Order-drop on aggressive shrink.** JAX falls back to Euler only for the
   first 2 steps / `dt_prev==0`; lacks Gear2's `h_curr/h_last < 0.1` protection
   (CPU `check_order_drop`).
6. **Per-domain tolerances.** CPU uses separate `iabstol`/`vabstol` for node vs
   branch rows and a `dx`+residual convergence test; JAX uses a single `abstol`
   and an L1 residual-only test.

## Features / API

7. **Integration-method selection.** JAX hardcodes `eval_method='gear'`;
   expose `method=` (euler/trap/gear2) like the CPU path.
8. **Coupled solver.** No `coupled_lte=True` (Fang DAC'13 Option A) equivalent.
9. **Pluggable step controllers.** CPU has the `StepController` strategy
   (`IntegralController`, `PIController`); JAX inlines one predictor.
10. **Device model bypassing** (`bypass` / `bypasstol`) — absent.
11. **Hooks:** no `provided_function` (feedback / Volterra) and no
    `cir.accept_step(...)` callback for stateful elements.
12. **DC init in `solve_batched`.** Starts at `x=0` ("For now, use 0"); only
    scalar `solve()` computes a DC operating point.

## Minor

`fixed_timestep`, `minbreak` handling, and a flexible `analysis` name (JAX
hardcodes `'tran'`).

## Suggested order of work

Port **#1 (device limiting)** and **#2 (breakpoint order reset)** first — those
two close most of the real-world robustness gap. Then #3/#4 for hard DC/steps.
The limiter must be JAX-friendly (branch-free / `jnp.where`-based) to stay
inside the compiled `jax.lax.while_loop`.
