"""E2's spike: is a PER-STAGE Newton worth putting on the JAX backend?

A MEASUREMENT, not a build.  Nothing here ships as a backend; the point is to
let a number decide whether any stage method (and then GLM) belongs on
`jaxtransient.py`, before anyone commits weeks to the port.  Asked for by the
tree's owner, 2026-09-16, with GLM explicitly NOT started until this reports.

WHY IT IS WRITTEN AGAINST THE CIRCUIT CALLABLES RATHER THAN THE BACKEND: the
backend reaches its integrator only through `compute_integration`, which knows
'euler' / 'gear' / 'trap' and nothing else, so a stage method cannot be asked
for there at all.  Both marches below therefore use the same hand-written
Newton on `cir.i/G/q/C/u(..., params_tree=...)`, which is what makes the
comparison fair: identical machinery, identical convergence tests, identical
`vmap`, identical fixed grid -- the ONLY difference is the companion.

    TR-BDF2   two implicit stages sharing one diagonal (gamma = 2 - sqrt 2):
              i(Y) + u(t_i) + (q(Y) - target)/(h a) = 0,  a = gamma/2
    Gear-2    one implicit solve per step:
              i(x) + u(t) + (3 q(x) - 4 q_n + q_{n-1})/(2h) = 0
              (first step backward Euler, as a BDF2 run starts)

⚠⚠ TWO TRAPS THIS HARNESS EXISTS TO NOT REPEAT.

1. `params_tree` IS CONSUMED ONLY THROUGH THE EVALUATION GROUPS, and
   `_eval_groups` is built LAZILY inside `Circuit.G` (circuit.py).  A march that
   never calls `G` leaves it empty, `batched_contributions` takes its
   `if not groups: return None`, and every lane silently runs the BUILD value.
   Measured before the fix: 8 bit-identical lanes, 3.48 V from the CPU sweep,
   `nonconv = 0`, and a meaningless "20x speedup" over the CPU loop -- exactly
   the failure `solve_batched`'s own comment records ("a parameter sweep that
   sweeps nothing, with no symptom").  `jax_setup()` builds the groups.

2. COMPARING AGAINST `solve_batched` IS NOT A COMPARISON.  It has no
   `fixed_timestep`, so gear runs ADAPTIVE: measured 7.3e-4 against the CPU
   reference while the fixed-step TR-BDF2 march reached 5e-14 -- different work
   at different accuracy, and the resulting "20-47x" says nothing about the
   integrators.  Hence the work-precision curves below.

FINDINGS (2026-09-16, 128 lanes, one CUDA device, the half-wave rectifier;
error is the worst lane against a validated TR-BDF2 reference at 12800 steps):

      steps   TR-BDF2 warm   err (V)    gear-2 warm   err (V)
        100        0.060s   1.571e-03       0.038s   1.585e-02
        200        0.100s   3.823e-04       0.053s   3.613e-03
        400        0.171s   9.432e-05       0.093s   8.442e-04
        800        0.314s   2.336e-05       0.173s   2.020e-04

  * both are cleanly second order (error ratios 4.04-4.4 per doubling);
  * gear-2 is 1.7-1.8x CHEAPER PER STEP (216-377 us vs 393-597 us), which is
    what one solve against two stages should cost;
  * TR-BDF2's error constant is ~10x smaller, so gear-2 needs ~sqrt(10) = 3.2x
    the steps for equal accuracy;
  * net, AT EQUAL ACCURACY, TR-BDF2 is ~1.4-1.6x faster (checked at two points:
    matching TR-BDF2@200 needs ~590 gear steps ~ 0.136s vs 0.100s; matching
    TR-BDF2@800 needs ~2350 gear steps ~ 0.51s vs 0.314s).

⚠ WHAT THIS DOES NOT SAY.  One fixture with FOUR unknowns, where per-step cost
is dominated by fixed overheads (a 3x3 dense solve, the traced while_loop), so
the per-step ratio need not transfer to a circuit with hundreds of unknowns --
and that is exactly where a batched backend earns its keep.  Fixed step, no LTE
control, no rejection/retry, no breakpoints, no rescue ladder: all of those are
what a production port must add, and they are where this backend's complexity
already lives.  A ~1.5x per-step-efficiency win is not by itself a reason to
port; measure at a larger `m` before deciding.

Run:  XLA_PYTHON_CLIENT_PREALLOCATE=false python benchmarks/stage_method_batched.py [lanes] [ref_steps]

`ref_steps` defaults to 12800 -- the figures above -- and is the expensive
compile; pass a smaller one (e.g. `8 1600`) for a quick check of the harness
itself.  ⚠ A reference only 2x finer than the finest competitor grid carries its
own error into the error column, so quote the defaults, not a quick run: the
`8 1600` check reads 1.761e-05 at 800 steps where the 12800-step reference reads
2.336e-05, and its validation gate passes at 1.79e-13.
"""
import sys
import time
import warnings
import numpy as np

TEND = 2e-3
LANES = 128
GRIDS = (100, 200, 400, 800)
NSTEP_REF = 12800
RELTOL = 1e-6
ABSTOL = 1e-12
XTOL = 1e-12
MAXITER = 30
R_LO, R_HI = 3e2, 3e4
CHECK_LANES = 4

GAMMA = 2.0 - np.sqrt(2.0)
A_DIAG = GAMMA / 2.0


def lane_values(n):
    return np.logspace(np.log10(R_LO), np.log10(R_HI), n)


def build(R, C, VSin, Diode, SubCircuit, gnd, r):
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b')
    c['R'] = R('b', gnd, r=float(r))
    c['C'] = C('b', gnd, c=1e-7)
    return c


def jax_setup():
    from pycircuit.circuit import circuit as circuit_mod, gnd
    from pycircuit.circuit.toolkit import jaxtoolkit
    circuit_mod.default_toolkit = jaxtoolkit
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VSin, Diode
    cir = build(R, C, VSin, Diode, SubCircuit, gnd, 1e3)
    ## trap 1 above: without this every lane runs the build value.
    cir._eval_groups = cir.toolkit.evaluation_groups(cir)
    return cir, cir.get_node_index(gnd), cir.get_node_index('b')


def make_marches(cir, iref, nstep):
    """`(trbdf2, gear)` marches over `nstep` fixed steps, sharing one Newton."""
    import jax
    import jax.numpy as jnp
    h = TEND / nstep
    m = cir.n

    def src(t, p):
        return cir.u(t, analysis='tran', params_tree=p)

    def newton(x_seed, make_FJ):
        def cond(st):
            _x, it, done = st
            return jnp.logical_and(jnp.logical_not(done), it < MAXITER)

        def body(st):
            x, it, _d = st
            F, J = make_FJ(x)
            J_sub = jnp.delete(jnp.delete(J, iref, axis=0), iref, axis=1)
            d = jnp.insert(jnp.linalg.solve(J_sub, -jnp.delete(F, iref)),
                           iref, 0.0)
            x_new = x + d
            I_scale = jnp.abs(J) @ jnp.abs(x_new) + jnp.abs(F)
            conv_f = jnp.all(jnp.delete(
                jnp.abs(F) < RELTOL * I_scale + ABSTOL, iref))
            conv_x = jnp.all(jnp.abs(d) < RELTOL * jnp.maximum(
                jnp.abs(x_new), jnp.abs(x)) + XTOL)
            return (x_new, it + 1, jnp.logical_and(conv_f, conv_x))

        x_out, _it, done = jax.lax.while_loop(
            cond, body, (x_seed, jnp.asarray(0), jnp.asarray(False)))
        return x_out, done

    def implicit(x_seed, target, t_i, a, p):
        scale = h * a

        def make_FJ(x):
            F = (cir.i(x, params_tree=p) + src(t_i, p)
                 + (cir.q(x, params_tree=p) - target) / scale)
            J = cir.G(x, params_tree=p) + cir.C(x, params_tree=p) / scale
            return F, J
        return newton(x_seed, make_FJ)

    def march_trbdf2(p):
        def step(carry, _):
            x0, t0, bad = carry
            q0 = cir.q(x0, params_tree=p)
            K0 = -(cir.i(x0, params_tree=p) + src(t0, p))
            x1, ok1 = implicit(x0, q0 + (GAMMA * h / 2.0) * K0,
                               t0 + GAMMA * h, A_DIAG, p)
            q1 = cir.q(x1, params_tree=p)
            c2 = 1.0 / (GAMMA * (2.0 - GAMMA))
            c0 = (1.0 - GAMMA) ** 2 / (GAMMA * (2.0 - GAMMA))
            x2, ok2 = implicit(x1, c2 * q1 - c0 * q0, t0 + h, A_DIAG, p)
            bad = bad + jnp.where(jnp.logical_and(ok1, ok2), 0, 1)
            return (x2, t0 + h, bad), None
        (xf, _t, bad), _ = jax.lax.scan(
            step, (jnp.zeros(m), jnp.asarray(0.0), jnp.asarray(0)), None,
            length=nstep)
        return xf, bad

    def march_gear(p):
        x0 = jnp.zeros(m)
        q0 = cir.q(x0, params_tree=p)
        x1, ok1 = implicit(x0, q0, h, 1.0, p)      # backward Euler opener

        def step(carry, _):
            xn, q_n, q_nm1, t_n, bad = carry
            ## the BDF2 companion as the same DC-flow form: a = 2/3,
            ## target = (4 q_n - q_{n-1}) / 3.
            target = (4.0 * q_n - q_nm1) / 3.0
            xnp1, ok = implicit(xn, target, t_n + h, 2.0 / 3.0, p)
            q_np1 = cir.q(xnp1, params_tree=p)
            return (xnp1, q_np1, q_n, t_n + h,
                    bad + jnp.where(ok, 0, 1)), None
        (xf, _q, _qm, _t, bad), _ = jax.lax.scan(
            step, (x1, cir.q(x1, params_tree=p), q0, h,
                   jnp.where(ok1, 0, 1)), None, length=nstep - 1)
        return xf, bad

    return march_trbdf2, march_gear


def timed(fn, tree, ib):
    """Cold (with compile) and warm wall time, and the per-lane finals."""
    import jax
    batched = jax.jit(jax.vmap(fn))
    t0 = time.perf_counter()
    xf, bad = batched(tree)
    xf.block_until_ready()
    cold = time.perf_counter() - t0
    t0 = time.perf_counter()
    xf, bad = batched(tree)
    xf.block_until_ready()
    warm = time.perf_counter() - t0
    return cold, warm, np.asarray(xf)[:, ib], int(np.sum(np.asarray(bad)))


def cpu_same_grid(lanes_subset, nstep):
    """The CPU's OWN TR-BDF2 on the SAME grid -- the yardstick's validation.

    ⚠ IT MUST BE THE SAME GRID.  Comparing the 12800-step march against the CPU
    at 1600 steps measured the COARSER grid's discretisation error and read as a
    disagreement: 5.75e-6, against 5.8e-6 predicted from the 800-step row.  On
    an identical grid the two backends agree to 1.9e-12 (measured, 4 lanes).
    """
    from pycircuit.circuit import gnd, numeric
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VSin, Diode
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit.integrator import TRBDF2Integrator
    out = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for r in lanes_subset:
            c = build(R, C, VSin, Diode, SubCircuit, gnd, r)
            tr = Transient(c, toolkit=numeric, integrator=TRBDF2Integrator(),
                           reltol=RELTOL)
            res = tr.solve(tend=TEND, timestep=TEND / nstep,
                           x0=np.zeros(c.n), fixed_timestep=True)
            out.append(float(np.asarray(res.v('b'))[-1]))
    return np.asarray(out)


def main():
    import jax
    import jax.numpy as jnp
    lanes = int(sys.argv[1]) if len(sys.argv) > 1 else LANES
    ## the reference grid is the expensive compile; a second argument makes the
    ## harness usable for a quick check without waiting for 12800 steps.
    ref_steps = int(sys.argv[2]) if len(sys.argv) > 2 else NSTEP_REF
    print('devices:', jax.devices())
    print('lanes=%d  tend=%g  reltol=%g   (fixed grids, identical Newton)'
          % (lanes, TEND, RELTOL))

    cir, iref, ib = jax_setup()
    rs = lane_values(lanes)
    tree = {'R': {'r': jnp.asarray(rs).reshape(-1, 1)}}

    ## the reference, and its validation on a grid the CPU also runs
    ref_trbdf2, _gear = make_marches(cir, iref, ref_steps)
    _c, _w, ref, ref_bad = timed(ref_trbdf2, tree, ib)
    idx = np.linspace(0, lanes - 1, CHECK_LANES).astype(int)
    val_march, _g = make_marches(cir, iref, 1600)
    _c2, _w2, val, _b2 = timed(val_march, {'R': {'r': jnp.asarray(rs[idx]).reshape(-1, 1)}}, ib)
    err = float(np.max(np.abs(val - cpu_same_grid(rs[idx], 1600))))
    print('reference: TR-BDF2 at %d steps, nonconv %d; the same march at 1600 '
          'steps against the CPU at 1600 steps, %d lanes: %.2e'
          % (ref_steps, ref_bad, CHECK_LANES, err))
    if err > 1e-9:
        print('  ⚠ THE MARCH DOES NOT REPRODUCE THE CPU ON AN IDENTICAL GRID -- '
              'the numbers below measure nothing until that is explained.')

    print('%8s %10s %12s %12s %12s %12s %10s' % (
        'steps', 'method', 'cold', 'warm', 'max err (V)', 'per-step us',
        'nonconv'))
    for nstep in GRIDS:
        tr_march, gear_march = make_marches(cir, iref, nstep)
        for name, fn in (('trbdf2', tr_march), ('gear2', gear_march)):
            cold, warm, finals, bad = timed(fn, tree, ib)
            e = float(np.max(np.abs(finals - ref)))
            print('%8d %10s %11.3fs %11.3fs %12.3e %12.2f %10d' % (
                nstep, name, cold, warm, e, warm / nstep * 1e6, bad),
                flush=True)


if __name__ == '__main__':
    main()
