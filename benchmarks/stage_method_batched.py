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

LARGER-m FOLLOW-UP (2026-09-17), which the entry above asked for: the same
harness at nsec = 1 / 10 / 25, i.e. m = 4 / 13 / 28, on grids 400/800/1600
against a 12800-step reference.  Every row converged (nonconv 0):

        m    per-step   err ratio   gear steps for   TR-BDF2 net
             tr/gear     gear/tr    equal accuracy   at equal acc
        4      1.97        8.55          2.92           1.48
       13      1.88        8.63          2.94           1.57
       28      1.77        8.67          2.94           1.66

  * BOTH HALVES SURVIVE.  Gear-2 is still cheaper per step at every size
    (1.44-1.97 over all nine cells, never approaching 1), and TR-BDF2 still
    wins at equal accuracy at every size (1.48-2.13).  The m=4 result was not
    an artefact of fixed overheads, which is exactly what was in doubt.
  * THE ERROR-CONSTANT RATIO IS THE STABLE PART: 8.55-9.32 across every cell,
    and 8.55/8.63/8.67 at the finest grid.  "Gear needs ~3x the steps" is a
    property of the two methods, not of the fixture's size.
  * ⚠ NO TREND IN m IS CLAIMED.  The per-step column looks monotone at the
    finest grid (1.97 -> 1.88 -> 1.77) but is NOT monotone at 800 steps
    (1.81 -> 1.87 -> 1.69), and every cell is a SINGLE timing: the
    grid-to-grid scatter at fixed m (m=28 reads 1.44 / 1.69 / 1.77) is as
    large as the variation across m.  Repeats would be needed to say more.
  * Cost grows superlinearly in m, as a dense solve should: at 1600 steps
    TR-BDF2 runs 390 / 784 / 2388 us per step for m = 4 / 13 / 28.

⚠ TWO FIXTURE DEFECTS WERE FOUND GETTING THERE, both documented at their sites:
the chain's far sections were DEAD at the original drive (see `build`), so a
first version of this sweep padded m with unknowns that did no nonlinear work
and printed identical error columns at nsec=10 and nsec=25; and `params_tree`
is keyed by CLASS, needing one column per element (see `lane_tree`).

⚠ WHAT THIS STILL DOES NOT SAY.  Twenty-eight unknowns is not the "hundreds"
where a batched backend earns its keep, and every size here is the same
rectifier chain -- one topology, one device model.  Fixed step, no LTE control,
no rejection/retry, no breakpoints, no rescue ladder: all of those are what a
production port must add, and they are where this backend's complexity already
lives.  A ~1.5-1.7x net win, now shown to hold from m=4 to m=28, is a real
argument for the stage method but still not by itself a schedule.

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


def lane_tree(jnp, rs, nsec=1):
    """The lane sweep as a `params_tree`: ONE COLUMN PER RESISTOR IN THE GROUP.

    ⚠ `params_tree` IS KEYED BY CLASS NAME, NOT BY INSTANCE NAME, and a hit
    REPLACES THE WHOLE GROUP'S params (`toolkit.batched_contributions`:
    `if cls.__name__ in params_tree: params = params_tree[cls.__name__]`).  So
    `{'R': {'r': ...}}` addresses EVERY resistor at once and must carry one
    column per resistor -- `(lanes, nsec)`, not `(lanes, 1)`.

    At `nsec=1` the two shapes coincide, which is why the original fixture ran
    for a day without showing this.  At `nsec=10` it is a LOUD failure, not a
    silent wrong answer: `vmap got inconsistent sizes ... one axis had size 10
    ... one axis had size 1`.  That is the opposite of the trap the sweep guard
    exists for, and worth knowing -- a shape that is merely too small fails
    closed here, while a KEY that stops matching fails open.

    Every section gets the same lane value, which also makes column order
    irrelevant: the group's order is its instance order, not necessarily the
    build order one would assume.
    """
    col = jnp.asarray(rs).reshape(-1, 1)
    return {'R': {'r': jnp.repeat(col, int(nsec), axis=1)}}


def build(R, C, VSin, Diode, SubCircuit, gnd, r, nsec=1):
    """`nsec` diode-RC sections in a chain.  `nsec=1` is the original fixture.

    ⚠ THE SECTIONS KEEP THE DIODE, and an RC ladder was rejected for this.
    A linear ladder converges in ONE Newton iteration, so it would measure
    linear-solve cost only and hide exactly the per-stage NONLINEAR work that
    separates TR-BDF2 (two implicit stages) from gear-2 (one).  Scaling `m`
    with a linear circuit would answer a different question than the one the
    spike asked.

    ⚠ THE SWEPT PARAMETER IS ADDRESSED BY CLASS, NOT BY INSTANCE NAME.  The
    section names (`R1..R{nsec-1}`) are irrelevant to `params_tree`, which
    overrides the whole `R` group at once; what matters is that its array
    carries one column per resistor.  See `lane_tree`, which is the only place
    that shape is built.  If the override ever stops matching BY KEY, every
    lane silently runs the build value -- the eight-identical-lanes trap this
    harness already caught once -- which is why `timed()` returns the per-lane
    finals and both `main()` and the sweep assert they differ.
    """
    ## ⚠ THE DRIVE SCALES WITH THE CHAIN, and this is not cosmetic.  At the
    ## original va = 5 the series drops (MEASURED: a steady 0.42 V per section)
    ## exhaust the source by section 8 -- nodes b9..b24 sat at 5e-11 down to
    ## 1e-161 V.  So a nsec=25 chain carried NINE live nonlinear sections and
    ## sixteen dead unknowns, and nsec=10 and nsec=25 printed error columns
    ## identical to every digit because they were solving the SAME live
    ## subcircuit.  Dead unknowns add linear-algebra and element-eval work but
    ## no per-stage NONLINEAR work -- exactly the defect the RC ladder was
    ## rejected for, rebuilt in another shape.  `nsec=1` keeps va = 5.0 to the
    ## bit, so the one-section findings quoted above are unchanged.
    ##
    ## THE SLOPE IS MEASURED WHERE IT IS USED.  The first 0.42 V/section was
    ## read off a chain that DIED at section 8 -- the shallow end -- and
    ## extrapolating it to 25 sections left the deepest node at 1.5e-6 V (the
    ## liveness guard caught that, which is what it is for).  Swept over
    ## va = 15.1..55 V the drop is FLAT: 0.431 V near the source, 0.474 V
    ## mid-chain, unchanged at every drive, with v(b) ~= 0.775*va - 1.06.  A
    ## chain therefore needs v(b) > ~0.45*(nsec-1) of headroom.
    ##
    ## ⚠ AND THE DRIVE IS BOUNDED ON BOTH SIDES: at va = 55 the fixed-step
    ## Newton starts failing (nonconv 44 at 400 steps) where 15.1..40 all give
    ## 0, so "more drive" is not a free fix and the smallest sufficient one is
    ## wanted.  0.85 V/section puts nsec=25 at va = 25.4, deepest node ~39 % of
    ## the first, nonconv 0.
    va = 5.0 + 0.85 * (int(nsec) - 1)
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=va, freq=1e3)
    prev = 'a'
    for k in range(int(nsec)):
        nd = 'b' if k == 0 else 'b%d' % k
        c['D%s' % ('' if k == 0 else k)] = Diode(prev, nd)
        name = 'R' if k == 0 else 'R%d' % k
        c[name] = R(nd, gnd, r=float(r))
        c['C%s' % ('' if k == 0 else k)] = C(nd, gnd, c=1e-7)
        prev = nd
    return c


def jax_setup(nsec=1):
    from pycircuit.circuit import circuit as circuit_mod, gnd
    from pycircuit.circuit.toolkit import jaxtoolkit
    circuit_mod.default_toolkit = jaxtoolkit
    from pycircuit.circuit.circuit import SubCircuit
    from pycircuit.circuit.elements import R, C, VSin, Diode
    cir = build(R, C, VSin, Diode, SubCircuit, gnd, 1e3, nsec)
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


def cpu_same_grid(lanes_subset, nstep, nsec=1):
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
            c = build(R, C, VSin, Diode, SubCircuit, gnd, r, nsec)
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
    ## E2 FOLLOW-UP: the entry's own next measurement is "the same harness at a
    ## larger m".  At four unknowns the per-step cost is dominated by fixed
    ## overheads (a 3x3 dense solve, the traced while_loop), so the 1.7-1.8x
    ## per-step gear advantage need not survive where a batched backend earns
    ## its keep.  `nsec` diode-RC sections give m = nsec + 3 (MEASURED, and not
    ## the 2*nsec + 2 first written here: R and C hang on the SAME section node,
    ## so a section adds ONE unknown, not two; nsec 1/10/25 -> m 4/13/28).
    nsec = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    print('devices:', jax.devices())
    print('lanes=%d  tend=%g  reltol=%g  nsec=%d   (fixed grids, identical Newton)'
          % (lanes, TEND, RELTOL, nsec))

    cir, iref, ib = jax_setup(nsec)
    print('m = cir.n = %d unknowns' % cir.n)
    rs = lane_values(lanes)
    tree = lane_tree(jnp, rs, nsec)

    ## the reference, and its validation on a grid the CPU also runs
    ref_trbdf2, _gear = make_marches(cir, iref, ref_steps)
    _c, _w, ref, ref_bad = timed(ref_trbdf2, tree, ib)
    ## ⚠ the eight-identical-lanes guard: a sweep that does not reach the
    ## elements gives every lane the build value and looks fast.
    spread = float(np.max(ref) - np.min(ref))
    assert spread > 1e-9, (
        'the %d lane finals span only %.3e V -- params_tree is not reaching '
        'the circuit and every lane is running the build value' % (lanes, spread))
    print('lane finals span %.4e V (the identical-lanes guard)' % spread)
    idx = np.linspace(0, lanes - 1, CHECK_LANES).astype(int)
    val_march, _g = make_marches(cir, iref, 1600)
    _c2, _w2, val, _b2 = timed(val_march, lane_tree(jnp, rs[idx], nsec), ib)
    err = float(np.max(np.abs(val - cpu_same_grid(rs[idx], 1600, nsec))))
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
