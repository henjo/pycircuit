"""`PSS`, the periodic steady state by shooting: the class, its parameters and
`solve` in its phases.  The class inherits its themes from the `_pss_*`
modules.
"""
from copy import copy
import numpy as np
import warnings
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.analysis import Parameter
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
from pycircuit.post import InternalResultDict
import pycircuit.circuit.analysis as analysis
from ._numerics import freq_analysis
from ._numerics import periodic_spline_weights
from ._pss_accuracy import _AccuracyChecks
from ._pss_events import _StateEvents
from ._pss_grids import _PeriodGrids
from ._pss_inner import _InnerTransient
from ._pss_newton import _ShootingNewton
from ._pss_periodic import _PeriodicStates
from ._pss_ppv import _PPVFloquet
from ._pss_replays import _FactoredReplays
from ._pss_walks import _PeriodWalks
from .diagnostics import algebraic_conditioning
from .diagnostics import topological_index


## Autonomy is decided STRUCTURALLY, not spectrally: a circuit is autonomous
## when nothing in it depends on `t` (`_is_autonomous` compares `u(t)` across
## the period against this relative tolerance), so a phase accumulator driven
## by a DC source is autonomous however energetically it oscillates.  No
## spectral-radius threshold separates the cases: an autonomous orbit reads a
## unit multiplier only AT its own period, and a lightly damped DRIVEN circuit
## reads near 1 too (a Q=1000 resonator: `exp(-pi/Q) = 0.99686`).
## History: `doc/shooting_history.md`, `AUTONOMOUS_U_TOL`.
AUTONOMOUS_U_TOL = 1e-12


class _SolveRun(object):
    """The state one `PSS.solve` carries between its phases (see `solve`):
    set by `_solve_prepare`, extended by `_shoot`, read by the rest."""

    def __init__(self, **fields):
        self.__dict__.update(fields)


class PSS(_ShootingNewton, _PeriodGrids, _StateEvents,
          _InnerTransient, _PeriodWalks, _FactoredReplays, _PPVFloquet, _AccuracyChecks, _PeriodicStates,
          Analysis):
    """Periodic Steady-State using shooting Newton iterations

    The algorithm is described in [1] p65.

     1. Kenneth S. Kundert, Jacob K. White, Alberto Sangiovanni-Vincentelli
        (1990)
        Steady-State Methods for Simulating Analog and Microwave Circuits
        Kluwer Academic Publishers
        ISBN 0792390695

    WHAT IS SOLVED.  `phi_T(x_0)` is one period of the inner transient
    (`Transient.solve_timestep`, driven one step at a time on a frozen
    grid), and shooting finds the fixed point `x_0 = phi_T(x_0)` by Newton
    on the monodromy `M = dphi/dx_0`.  A DRIVEN circuit is solved at the
    caller's period.  An AUTONOMOUS one -- nothing in it depends on `t`
    (`_is_autonomous`) -- has a one-parameter family of periodic solutions,
    so `I - M` is singular along the orbit and it is solved on the
    FREE-PERIOD system: unknowns `(x_0, T)`, bordered by the period column
    and one phase row (`phase_rule`).  `solve` runs in phases
    (`_solve_prepare`, `_shoot`, `_report_convergence`, `_replay_orbit`,
    `_report_lte`, `_check_fundamental`, `_orbit_results`,
    `_closing_polish`); the themes live in the `_pss_*` mixins.

    **THREE CONVERGENCE CHECKS NEST HERE, AND THEY ARE NOT INTERCHANGEABLE.**

    1. the per-timestep Newton, which solves the discretised circuit
       equations at one time point;
    2. the local truncation error, which decides how far the DISCRETE
       trajectory is from the true one;
    3. the shooting Newton, which finds the periodic point of the discrete
       map.

    Two rules order them.  **(3) cannot be tighter than (1)**: the period
    map is only KNOWN to the accuracy of the per-timestep solves, so the
    shooting residual has a floor there whatever its Jacobian is.  And
    **LTE must not run per shooting iteration**: an
    adaptive grid makes the step sequence a function of `x0`, so the period
    map stops being smooth and (3) loses its quadratic rate.  Choose the
    grid once, freeze it, shoot on it.

    So (2) is not a controller inside a solve.  It is a MEASUREMENT of
    ACCURACY on the converged period (`max_lte`, `total_lte`,
    `max_lte_seam`; see `_report_lte`), and grid adaptivity enters BETWEEN
    solves (`grid=`, `lte_grid`).  ⚠ A converged shooting solve is not by
    itself evidence of a correct answer: each method's own numerical
    damping is invisible to (1) and (3) -- on a Q=20 resonator at 100
    points/period euler, gear2 and trap all converge, at 56 %, 1.2 % and
    0.05 % below the analytic peak.  Driving `Transient` buys one
    integrator definition, the limiting/PCNR machinery, breakpoints and the
    order drop -- not (2) as a controller.

    THE SHOOTING JACOBIAN.  Every linear multistep method here writes its
    companion as

        iq_n = sum_k a_k q_{n-k}  +  b iq_{n-1}

    so ONE recursion differentiates all of them:

        S    = sum_{k>=1} a_k C_{n-k} Px_{n-k}  +  b Pq
        Px_n = -Jf_n^-1 S
        Pq_n = a_0 C_n Px_n + S

    Euler is `b = 0` reaching back one step, trapezoidal `b = -1` reaching
    back one, Gear-2 `b = 0` reaching back two.  The coefficients come from
    `Integrator.companion_coefficients` -- from the integrator that ACTUALLY
    ran, so an order-dropped opening step contributes its own -- never
    transcribed here.  A wrong Jacobian is worse than none: an x-only
    monodromy for trapezoidal converges slower than no Jacobian at all.

    FOUR THINGS ENLARGE THE SYSTEM, and they are named apart so that a
    statement about one is not applied to another:

      the `(x, iq)` MONODROMY   what trapezoidal's period map differentiates.
                                Its recursion carries a companion current, so
                                an x-only monodromy is structurally
                                incomplete.  About the DERIVATIVE, not the
                                unknowns.
      the FREE-PERIOD system    what an autonomous circuit solves: unknowns
                                `(x0, T)` with a phase condition, because the
                                period is not given.  `func_autonomous` in
                                `_shoot`, on every kind's `_pmap`.
      a SOLVED ENTERING HISTORY unknowns `(x0, x_{-1})`, 4b.  About where
                                the period map STARTS, not how long it runs.
                                DRIVEN circuits.  `self.solved_history`.
      the COMPOSED system       unknowns `(x0, x_{-1}, T)`, 4c -- the second
                                and third TOGETHER, for an autonomous
                                circuit under a two-step method.  NOT a
                                synonym for either half.

    THE PERIOD MAP PER KIND (`_map_kind`):

      * 'plain' (euler, trap, theta): one entering unknown.  By default it
        is `x_in`, and `x(0)` is manufactured from it by one order-dropped
        Euler step.  That L-STABLE OPENER is what makes trapezoidal's
        shooting problem well posed: trapezoidal is A-stable but not
        L-stable, so it maps the null space of the singular MNA `C` by
        exactly `-1` per step (`m - rank(C)` such modes, one per algebraic
        variable), and ANY period map `A_trap^K` without an L-stable
        opening step is singular at even `K` on every MNA circuit.  (Three
        reformulations without the opener -- solving for `x_{-1}`, for
        `iq_0`, or making `iq_0` dependent -- were built and fail; see
        history.)  The cost is the frame: the Jacobian is taken with
        respect to `x_0` while the unknown is `x_in`, whose true `dF/dx_in`
        is singular, so the iteration is a CONTRACTION with a linear rate,
        not a Newton.  `x0_unknown=True` makes `x_0` the unknown -- Aprille
        & Trick's canonical formulation (Proc IEEE 60(1) 108-114, 1972) --
        with an exact Jacobian; see `solve` for its trade-off.  `theta`
        (trapezoidal biased by `C h`, `ThetaIntegrator`) takes the opener
        out.
      * 'pair' (Gear-2), 4b: the companion reads `q_{n-2}`, which on the
        plain path is the entering stand-in -- the SEAM, 54 % of Gear-2's
        error at 100 points/period and a growing share under refinement.
        (Euler's and trapezoidal's seams measure zero: they read only
        `q_{n-1}` and `iq_{n-1}`.)  So `(x_0, x_{-1})` are unknowns
        together and BOTH close (the pair walk, `_walk_lmm`), with an exact
        Jacobian.  A k-step method turns a first-order problem into a k-th
        order discrete one that needs k conditions, not one (boundary value
        methods: Brugnano & Trigiante; Ascher, Mattheij & Russell, SIAM
        1995).
      * 4c, the pair on an AUTONOMOUS circuit, is the COMPOSED system:
        without the second closure the missing condition lands in the
        PERIOD.  Its 2m x 2m map carries Gear-2's parasitic root, `(1/3)^N`
        over a period, far below the physical multipliers, so `max |eig|`
        picks the physical one (a test pins the gap).
      * 4d: the DAE gives `qdot = -(i + u)` exactly with no solve (a
        converged step satisfies `i(x) + iq + u = 0`), while `qddot` needs
        the singular `C`: which is why trapezoidal's `iq_{-1}` is exactly
        initialised and Gear-2's second charge is not.  A pseudo-history
        built from `x_0` ("method H") fixes the orbit's shape but not an
        autonomous circuit's frequency, and is not shipped (see history).
      * 'stage' (radau, trbdf2, esdirk43): self-starting, the unknown IS
        `x_0`, `M` the dense stage product; no opener, no seam.
      * 'glm' (Nordsieck GLM2-4): the unknown is `x_0`, the map ``x_0 ->
        x_N`` exact through the linearised startup (`_walk_glm`); the
        factored period acts on the multivalue state, and `ppv` /
        `floquet_modes` read it on the state (`_GLMPeriod.state_map`).

    ⚠ "SEAM-FREE" IS NOT "ERROR-FREE": the solved history removes the seam,
    not ordinary interior discretisation error (Gear-2 at 100 points/period
    still returns 19.89297 V against 20 V).

    The circuit literature (Kundert, ICCAD 1997; Gourary et al., MES 2019)
    writes the shooting map for one-step methods, which the plain path is;
    the gap opens only when a two-step companion is handed to it.  The
    standard case against Gear-2 (Wambacq et al., "CAD for RF circuits")
    is about ADAPTIVE step-ratio changes; inside one frozen solve it does
    not bite, and `_period_grid` warns when a caller's ratios leave the
    zero-stability bound.

    ⚠ HIGH-Q OSCILLATORS, whose physical multipliers cluster near 1, are
    where this degrades, from one cause surfacing in three places: `max
    |eig|` and `_spectral_report`'s split assume the unit root is
    identifiable; the phase row removes the singularity from the unit
    multiplier only, so the bordered system stays nonsingular but badly
    conditioned (a slow, loud Newton rather than a `LinAlgError`); and PPV
    eigen-selection would pick the same root.  Independent sources say so
    (Demir, IJCTA 2000; Bizzarri et al.; Demir & Roychowdhury).  Demir's
    2003 remedy is why a PPV built here goes to the augmented solve --
    sampled `C(t) u_1(t)` as the bordering row `q` -- and NOT to the
    eigenvectors this method returns.  Not measured here on a high-Q case.

    ⚠ FREE-PERIOD PITFALLS, properties of the free-period system and not of
    any element (4e): `k*T` solves the periodicity condition whenever `T`
    does and the solve follows its seed (`_check_fundamental` warns and
    sets `fundamental_period`; driven runs are exempt); and `T = 0` is a
    REGULAR root of every autonomous shooting system, which any seed below
    the fundamental is drawn to (`_free_period_solve` demotes and names it
    -- seed at or above the expected period).  A scalar `Idtmod` closes
    through its wrap (the gauge shift and `_fold_periodic`).  Still open: a
    pure phase accumulator has no amplitude to pin, so its orbit is a ramp
    and only the phase condition distinguishes one starting phase from
    another.

    GRIDS (recorded scope item 5).  `solve(grid=...)` takes step FRACTIONS
    of the period and freezes them (`_period_grid`); an autonomous run
    rebuilds the grid at the current `T` on every residual evaluation, so
    `dh/dT = h/T` holds.  The uniform grid is what costs on a stiff
    relaxation oscillator: van der Pol at mu = 100 needs 20000 uniform
    points and solves on its own 1105-step LTE-chosen grid (`lte_grid`),
    pinned by
    `test_the_lte_chosen_grid_solves_van_der_pol_through_the_analysis`.
    `_period_grid` opens on the grid's own finest step (a coarse
    manufactured opening step fails the inner Newton), and
    `_install_history` takes the entering step `h_prev = hs[-1]`, not
    `hs[0]`.

    MATRIX-FREE (recorded scope item 6; Telichevesky, Kundert & White, DAC
    1995): `solve(matrix_free=True)` never forms `J_phi`.  What it removes
    is the sensitivity propagation (`N` steps of `_step_sensitivity`, each
    O(m^3) twice over), which passes ~30 % of a traversal near m = 220;
    GMRES iterations track the number of SLOW modes, not m.  Worth it above
    m ~ 250; see `solve` for its costs.

    TOLERANCES.  `reltol`, `iabstol` and `vabstol` mean exactly what they
    mean on `Transient` -- the tolerances of the TRANSIENT solution,
    applied to the per-timestep Newton, with the two absolute floors
    applied PER UNKNOWN in both flavours by
    `analysis.newton_tolerance_vectors`, the single definition all three
    analyses read.  Nothing here rescales them.  The shooting criterion is
    expressed against that one: `steadyratio` (>= 1, default 1) multiplies
    it, so by default the shooting solve is held to the SAME relative
    tolerance as the transient, and raising it buys fewer shooting
    iterations for a looser periodic steady state.

    CHOOSING `method`.  The default is chosen for ACCURACY, not cost (owner
    decision: "we do not want to fool the user; instead they should change
    integrator and know the impact").  Every entry is MEASURED in this tree
    (doc/pss_roadmap_260902.md A10, doc/pss_log_260902.md), none quoted
    from a textbook order, and every method's error is estimable by
    refinement (`grid_error`):

      method    order   period error   above its monotonicity limit
                        (ppm @ 400pt,  (h_FE = 1/steepest slope)
                        Q=1e4 vdP)
      radau     5       5.7e-10        OUTPUT stays in range; only the
                                       STAGES leave it, <= 0.3 % of the
                                       swing above ~4 h_FE (seen only on a
                                       scalar square-wave fixture)
      esdirk43  4       5.3e-05        output in range; stages leave it
                                       <= 3 % above ~2.5 h_FE
      trbdf2    2       10.1           OUTPUT rings above 2.4 h_FE
      trap      2       20.8           OUTPUT rings above 2 h_FE
      gear      2       83.1           a small ring at h >= 10 tau
      euler     1       --             never rings; damps the orbit it is
                                       asked to find (13 % of the amplitude
                                       at 20 pts/period)

      * radau: self-starting, L-stable, carries its own monodromy.  Index
        2: the period keeps classical order, the algebraic unknowns
        converge at the stage order (3), and a period timed by such a
        variable still converges at order 5.  A comparator's inner step
        Newton can fail undamped at every grid (hence the line search as
        the last resort), and `warping_estimate` is unreliable across a
        comparator edge (0.09-0.65 of the true error: the interpolant does
        not resolve it).
      * esdirk43: the MOST expensive method at every size (four sequential
        stage solves), never the cheap alternative.
      * trbdf2: order 2, contractive; the alternative to price against
        `grid_error` above a few hundred unknowns.
      * trap: its PPV surfaces are its TR-BDF2 twin's (`monodromy_twin`).
      * gear: the PPV, the diffusion constant, the period (free or driven),
        the adjoint modes and every noise fold are second order on any grid
        whose step ratios stay inside the zero-stability bound 1+sqrt(2);
        beyond it (3:1) the period is FIRST order and `_period_grid` warns.
        With `break_events` it is second order across a landed edge, at a
        CONSTANT cost at the step after a hard corner (30x trap's there).
        Index 2: order 2 in both subspaces.  Adjoints: the exact transpose
        for the PPV and the noise folds, a second-order continuous adjoint
        for the mode shapes on non-uniform grids (dense to 8 unknowns,
        Ritz-certified Arnoldi above).  One real factorisation per step,
        3.0-11.8x cheaper than radau as n grows; TR-BDF2 matched it
        1.5-1.7x better at equal accuracy on E2's fixture, so price both
        with `grid_error`.  On a given grid gear's error is a steady ~8x
        trbdf2's, both second order -- choose by that factor, and refine
        before believing a good number (a window can sit on a SIGN CHANGE
        of the error).  For NOISE on a coarse or folded grid choose trbdf2
        or radau (exact per-stage injection); gear's covariance is first
        order in h/tau.  It does not CERTIFY a free-period solve at 1e-14
        below ~400 pts at Q=1e4.
      * STATE EVENTS (a comparator, a threshold switch; `state_events`):
        use radau; trbdf2's and gear's own second-order error on the
        switch-off decay dominates there.
      * ADAPTIVE GRIDS (`lte_grid`): the adaptive run's reltol must be
        tighter than the transient's habit (1e-5 is too coarse above mu ~
        10 for every method); pass the period it observed
        (`pss.lte_period`) to `solve`, since the grid is fractions of it.

    Cost: radau is one real plus one complex factorisation per step (a
    coupled 3n system) against one per step for the others.  At small n it
    is cheaper on wall-clock at equal accuracy anyway (60 points beat
    trap's 480 on both axes); per-step assembly dominates to n ~ 100, and
    the 3n stage factorisation dominates from a few hundred unknowns
    (radau/gear ~3-4x up to n ~ 300, 11.8x at n = 1002).  Above ~300
    unknowns the accuracy is bought at an order of magnitude in
    wall-clock: price trbdf2 against `grid_error` on your circuit.
    `grid_error` and `warping_estimate` price the trade on YOUR circuit.

    History: `doc/shooting_history.md`, `PSS` (class docstring).
    """

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name',
                   ## Sources supply their time-domain waveform only for an
                   ## analysis name in timedomain_analyses (('dc','tran')); any
                   ## other name makes cir.u(t) zero and the solve unexcited.
                   default='tran'),
         Parameter(name='reltol', 
                   desc='Relative tolerance', unit='', 
                   default=1e-4),
         Parameter(name='iabstol', 
                   desc='Absolute current error tolerance', unit='A', 
                   default=1e-12),
         ## DC, Transient, JAXTransient and PSS share one meaning and one
         ## default; the reason is at `Transient.vabstol`
         Parameter(name='vabstol', 
                   desc='Absolute voltage error tolerance', unit='V', 
                   default=1e-6),
         Parameter(name='maxiter',
                   desc='Maximum number of iterations', unit='',
                   default=100),
         ## Forwarded to the inner Transient so PCNR (the junction-continuation
         ## limiting) reaches the shooting per-step solve too -- it lives in
         ## `Transient.solve_timestep`, which PSS DOES call, so no per-accepted-
         ## step machinery is needed (unlike breakpoints / continuation rescue,
         ## which are armed in `Transient.solve` and stay out of reach).
         Parameter(name='pcnr',
                   desc='Use Predictor/Corrector Newton-Raphson instead of '
                        'limiting in the inner transient; off by default',
                   unit='', default=False),
         ## Which step lengths move with an unknown period.  'proportional'
         ## rescales every step with T; 'closing' keeps a caller's inner steps
         ## at their absolute lengths and lets the last step close the period.
         ## 'auto' is 'closing' on a caller's grid for an autonomous run and
         ## 'proportional' otherwise -- see `_period_grid`, `_solve_prepare`
         ## and `_closing_polish`.
         Parameter(name='period_column',
                   desc="'auto' (= 'closing' + proportional polish on a caller's grid, 'proportional' on a uniform one), 'proportional' or 'closing': "
                        "which step lengths depend on an unknown period; see "
                        "the note at the policy in solve()",
                   unit='', default='auto'),
         ## `method='theta'`'s one knob, as the DIMENSIONLESS `C T`: the
         ## transferable quantity (see `ThetaIntegrator.DEFAULT_CT` for why
         ## `h` cancels), which `_theta_biased` turns into the rate THIS period
         ## needs.  `None` takes the measured knee.
         Parameter(name='theta_ct',
                   desc="method='theta' only: the null(C) damping one PERIOD "
                        'applies, as the dimensionless product C*T. None '
                        'takes ThetaIntegrator.DEFAULT_CT (the measured '
                        'knee). Ignored by every other method.',
                   unit='', default=None),
         ## `reltol` MEANS THE SAME THING IN EVERY ANALYSIS: the relative
         ## tolerance of the transient solution, applied to the per-timestep
         ## Newton here exactly as `Transient` applies it, never rescaled.
         ## `steadyratio` expresses the SHOOTING criterion against it
         ## (shooting reltol = reltol * steadyratio).  It is >= 1 because the
         ## period map is only KNOWN to the accuracy of the inner solves;
         ## raise it for a looser periodic steady state in fewer iterations.
         ## The LTE floors are separate from the Newton ones (one knob must not
         ## move both criteria -- see `Transient`); same names, defaults and
         ## meaning, so this reports the number a transient would control on.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage floor for the truncation-error '
                        'estimate', unit='V', default=1e-12),
         Parameter(name='lte_iabstol',
                   desc='Absolute current floor for the truncation-error '
                        'estimate', unit='A', default=1e-12),
         Parameter(name='TRTOL',
                   desc='Truncation error over-estimation factor (SPICE '
                        'TRTOL / lteratio in a commercial simulator)', unit='', default=7.0),
         Parameter(name='relref',
                   desc="What the relative LTE tolerance is measured "
                        "against: 'pointlocal', 'alllocal' or 'sigglobal'",
                   unit='', default='sigglobal'),
         Parameter(name='steadyratio',
                   desc='Shooting tolerance as a multiple of reltol (>= 1); '
                        '1 holds the shooting solve to the same relative '
                        'tolerance as the transient, larger relaxes it',
                   unit='', default=1.0),
         ## The default is `radau` (owner decision), for ACCURACY: the floor of
         ## this stack is DISCRETISATION, growing linearly in Q, and radau --
         ## order 5, self-starting (no manufactured opener, so no seam in the
         ## period map), L-stable -- carries its own monodromy, so an
         ## autonomous run takes NO TR-BDF2 twin and reads its own spectrum
         ## (`monodromy_twin`, `carries_own_monodromy`).  It costs more per
         ## step on a fine grid, but for any oscillator surface the twin
         ## dominates `trap`'s cost, and at equal accuracy radau wins on both
         ## axes (60 points beat trap's 480).  `trap` remains one argument away
         ## for a cheap coarse answer.
         ## History: `doc/shooting_history.md`, `PSS.parameters`.
         Parameter(name='method',
                   desc="Integration method for the inner transient: 'radau' "
                        "(default, order 5), 'esdirk43', 'trbdf2', 'theta', "
                        "'gear' (BDF-2), 'trap' or 'euler'. The default is "
                        "chosen for ACCURACY; the class docstring's 'CHOOSING "
                        "method' table states, from measurement, what each "
                        "alternative gives up in order, error estimability, "
                        "monotonicity and cost",
                   unit='',
                   default="radau")]        

    
    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(PSS, self).parameters + self.parameters
        super(PSS, self).__init__(cir, **kvargs)
        ## The reference row is fixed for the analysis, and both the shooting
        ## loop and the Transient this drives need it.
        self.irefnode = self.cir.get_node_index(
            gnd if irefnode is None else irefnode)
        self._tran = None
        ## Only assembled when the period is an unknown; an extra assembly
        ## per timestep is not worth paying on the fixed-period path.
        self._want_dfdh = False
        self._dfdT = None
        ## WHICH PERIOD-COLUMN CONVENTION `want_dT` USES; `_solve_prepare`
        ## sets it per solve from `period_column`.  'proportional': every
        ## step scales with `T`, `dh_i/dT = h_i/T`.  'closing' (the commercial
        ## convention): the inner transient owns the steps and the LAST one
        ## is placed on the period boundary, so `dh_i/dT = 0` inside and
        ## `dh_N/dT = 1` at the close.  The proportional column is `O(h)`
        ## wrong on a smooth uniform grid and worse at large step ratios.
        self._period_column = 'proportional'
        self._dfdh = None
        self._want_lte = False
        ## the state-event fractions the last solve landed on (None: none),
        ## and `dtheta/dx_0` from the bordered solve (K x m)
        self._state_event_fracs = None
        self._event_sensitivity = None
        self._event_columns = None
        self.stage_one_converged = None
        self._captured = {}
        ## The caller's step fractions, or None for the uniform grid.  Read
        ## by the autonomous closures, which rebuild the grid at the current
        ## `T` on every residual evaluation.
        self._grid_fracs = None
        ## the period `lte_grid` observed in its adaptive run (see there)
        self.lte_period = None
        self._lte = None
        self._lte_seam = False
        self._lte_valid = True
        self._history_is_solved = False
        ## Set by `solve`: whether the entering history joined the unknowns.
        self.solved_history = False
        ## Reported by `solve`: the peak normalised truncation error over the
        ## converged period, and where in the period it fell.  None until a
        ## solve has run, or when the grid was too short to difference.
        self.max_lte = None
        self.max_lte_time = None
        self.max_lte_seam = None
        self.total_lte = None
        ## Set by `solve` when an autonomous run lands on a multiple of the
        ## fundamental; None when it did not (or on a driven run).
        self.fundamental_period = None
        ## What `factored_period()` needs to replay the CONVERGED period:
        ## which seed, which grid, which opening.  Written at the end of
        ## `solve`, and cleared at its start so a failed or interrupted run
        ## cannot leave a previous solution's state readable as this one's.
        self._period_state = None
        self._factored_period_cache = None
        ## ⚠ WHICH MONODROMY THE OSCILLATOR SURFACES READ -- see
        ## `monodromy_twin`.  Selects the method that supplies the PPV,
        ## Floquet modes and factored period when a one-step LMM (trap/euler)
        ## solved an autonomous circuit, whose OWN monodromy is first-order
        ## on a limit cycle (the opener seam).  'trbdf2' (DEFAULT): a TR-BDF2
        ## twin on the same grid -- self-starting, no opener, 12-32x more
        ## accurate on lambda2 than the Gear-2 twin at practical step counts.
        ## 'gear': a Gear-2 twin.  'radau' / 'esdirk43': those twins (radau
        ## the most accurate by far, at about trbdf2's cost -- see
        ## `monodromy_twin`).  'native': the run's OWN plain factorisation --
        ## ⚠ the WORST under a one-step method (its `Q` DIVERGES under
        ## refinement); it is for the gates that measure that defect, not for
        ## results.  gear and stage-method runs are self-sufficient
        ## (second-order native monodromy) and ignore this.
        ## History: `doc/shooting_history.md`, `PSS.__init__`.
        self.monodromy = 'trbdf2'
        self._monodromy_twin = None
        self._twins = {}
        self._solve_kwargs = {}
        ## Set to None by `solve` on the tstab path only -- see there.
        self.tstab_state = None

    def _is_autonomous(self, times):
        """True when nothing in the circuit depends on `t`.

        Exact where a spectral test is not: see `AUTONOMOUS_U_TOL`.  The
        source vector is evaluated at EVERY point of the grid and compared
        with the first; a circuit driven only by DC -- a VCO macromodel, a
        phase accumulator, an LC or ring oscillator -- has a constant `u`
        and a one-parameter family of periodic solutions, which is what
        makes fixed-period shooting ill-posed for it.
        """
        u0 = np.asarray(self.cir.u(times[0], analysis=self.par.analysis),
                        dtype=float)
        scale = max(float(np.max(np.abs(u0))), 1.0)
        ## ⚠ EVERY POINT ON THE GRID, NOT A STRIDE THROUGH IT: a sparse sample
        ## misses a narrow pulse (PWM, sampling clocks, S/H, mixer LOs), and a
        ## clock misread as autonomous is routed to the free-period system,
        ## which solves for `T` and DISCARDS the period the caller asked for
        ## -- `DEGENERATE_PERIOD_FACTOR` tests the magnitude of `T`, not
        ## whether the circuit was driven, so it cannot catch it.  Cost: `N`
        ## evaluations of `u` once per solve, exiting at the first sample
        ## that differs (the common case for every driven circuit).
        ## History: `doc/shooting_history.md`, `PSS._is_autonomous`.
        for t in times[1:]:
            u = np.asarray(self.cir.u(t, analysis=self.par.analysis),
                           dtype=float)
            if float(np.max(np.abs(u - u0))) > AUTONOMOUS_U_TOL * scale:
                return False
        return True

    ## The one mapping from `method` to a class.  Read by `_transient` to
    ## build the integrator and by `_companion_reach` to ask how far it
    ## reaches, so the two can never disagree about which method is running.
    @classmethod
    def _integrator_for(cls, method):
        from pycircuit.circuit.integrator import (EulerIntegrator,
                                                  TrapezoidalIntegrator,
                                                  ThetaIntegrator,
                                                  Gear2Integrator)
        from pycircuit.circuit.integrator import (TRBDF2Integrator,
                                                  RadauIIA3Integrator,
                                                  ESDIRK43Integrator,
                                                  GLM2Integrator,
                                                  GLM3Integrator,
                                                  GLM4Integrator)
        ## THE single method -> integrator map, and the one place method names
        ## are validated: an unknown name raises the ValueError here rather than
        ## a KeyError three frames down.  The polymorphic predicates
        ## (`_companion_reach`, `needs_x0_unknown`, ...) call this before the
        ## solve()-level whitelist runs, so the validation must live where the
        ## lookup does.
        table = {'euler': EulerIntegrator,
                 'trap': TrapezoidalIntegrator,
                 'trapezoidal': TrapezoidalIntegrator,
                 ## `theta` is trapezoidal biased by `C h` -- see
                 ## `ThetaIntegrator`. It takes the L-stable opener OUT, which
                 ## is the one thing the other one-step LMMs cannot do.
                 'theta': ThetaIntegrator,
                 'gear': Gear2Integrator,
                 'gear2': Gear2Integrator,
                 'trbdf2': TRBDF2Integrator,
                 'radau': RadauIIA3Integrator,
                 'esdirk43': ESDIRK43Integrator,
                 ## Nordsieck GLMs: stage order = order, one factorisation per
                 ## step, no index-2 order split.  Driven PSS only (no free
                 ## period), and the period map is on the MULTIVALUE state --
                 ## see `factored_period_glm`.
                 'glm2': GLM2Integrator,
                 'glm3': GLM3Integrator,
                 'glm4': GLM4Integrator}
        try:
            return table[method]()
        except KeyError:
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau', 'esdirk43', 'glm2', 'glm3' or 'glm4', not %r"
                % (method,))

    def _resolve_break_events(self, requested):
        """`break_events`, ON for every method when not given.

        Landing source discontinuities on grid points helps every method,
        gear included: on a pulsed RC gear + events wins at most N and stays
        second order.  Its one cost is a CONSTANT, not an order: a two-step
        formula takes an O(h^2 [x'']) hit at the ONE step after a corner,
        where its history straddles the jump in x'' (about 30x trap's
        there).  At a TRUE jump (tr = 0) the gain needs `event_grid` to keep
        both ends of the clamped ramp (see there).

        ⚠ An explicit `True`/`False` is honoured untouched; this only fills in
        `None`.

        History: `doc/shooting_history.md`, `PSS._resolve_break_events`.
        """
        if requested is not None:
            return bool(requested)
        return True

    def _resolve_x0_unknown(self, requested):
        """`x0_unknown`, defaulted from the circuit's TOPOLOGY when not given.

        Conditional, not a global default, because `x0_unknown` is NOT free:
        trapezoidal still needs an L-stable opener, so switching it on moves
        the Euler step INSIDE the period, where it degrades the ORBIT rather
        than just the opening (a `Q = 20` resonator against its analytic
        20 V peak: 19.76939 with it, 20.01273 without, at 100 points).

        ⚠ ON AN INDEX-2 NETLIST THE TRADE REVERSES. The manufactured opening
        step is INCONSISTENT there: the constraint fixes the algebraic
        variable at a value the step cannot produce, so trapezoidal returns
        EXACTLY 2x on an L-I cutset -- and on an even number of steps reports
        `converged` while doing it (roadmap section 0k). `x0_unknown`
        removes it at every parity because `x(0)` becomes a genuine unknown.

        Three refusals, each deliberate:

          * **an explicit `True`/`False` is honoured untouched** -- this only
            fills in `None`;
          * **a two-step method is left alone**, because its solved-history
            formulation already solves for `x(0)` and was never affected;
          * ⚠⚠ **a PROVISIONAL verdict does not trigger it, and that is a
            REFUSAL ON THE THEORY RATHER THAN CAUTION.** With a controlled
            source in the loop or cutset the index is not bounded by 2 and
            need not be a function of the topology at all (see
            `topological_index`), and the remedy is justified for an
            INDEX-2 algebraic row only.  **If the true index is 3 the
            remedy is not known to apply, and switching it on would mask a
            worse problem while reporting a fix.**  Leaving a known defect
            visible is the better failure mode.  The same goes for a
            structurally singular netlist, which has no index at all.

            Index 3 IS reachable on this element set and is VALUE-dependent:
            a VCVS of gain `g` inside a C-V loop (C1 v->b, C2 v->gnd, source
            b->gnd = g*v) is index 3 exactly on `g* = 1 + C2/C1`, and a gain
            near `g*` carries a nilpotent tail -- the NEIGHBOURHOOD is what
            bites.  `topological_index` reads 1 there (the `idx != 2` guard
            fires before `provisional` is consulted), and **the index-3
            netlist DC-solves cleanly and silently**: the case to guard
            against is the one that solves and looks healthy (see
            `_warn_if_the_block_disagrees`).

        Warns when it fires, because a silently different formulation is the
        kind of thing that makes a later measurement inexplicable.

        History: `doc/shooting_history.md`, `PSS._resolve_x0_unknown`.
        """
        if requested is not None:
            return bool(requested)
        ## ⚠ EVERYTHING BELOW IS BEST-EFFORT AND MUST NEVER RAISE.  This runs
        ## BEFORE `solve` validates its own arguments, and a defaulting helper
        ## has no business changing which exception an invalid call raises (a
        ## bad `method` must reach `solve`'s `ValueError`, not a `KeyError`).
        try:
            if self._solves_history():
                return False
            idx, info = topological_index(self.cir)
        except Exception:
            return False
        if idx != 2 or info['provisional'] or info['ill_posed']:
            self._warn_if_the_block_disagrees(idx, info)
            return False
        where = (('C-V loop: ' + ', '.join(info['loop'])) if info['loop']
                 else ('L-I cutset: ' + ', '.join(info['cutset'])))
        warnings.warn(
            'PSS: this netlist is index 2 (%s), where the manufactured '
            'opening step is INCONSISTENT -- it starts an algebraic variable '
            'at a value the constraint forbids, and trapezoidal carries that '
            'seed forever (exactly 2x on an L-I cutset, reported as CONVERGED '
            'on an even number of steps). Solving for x_0 directly instead; '
            'pass x0_unknown=False to override.' % where,
            RuntimeWarning, stacklevel=3)
        return True

    def _warn_if_the_block_disagrees(self, idx, info):
        """Say so when the TOPOLOGICAL index reads below 2 and the NUMERIC
        algebraic block says otherwise.  Never raises, never changes behaviour.

        ⚠⚠ THE CASE THIS EXISTS FOR IS THE ONE THAT SOLVES AND LOOKS HEALTHY:
        the index-3 VCVS netlist of `_resolve_x0_unknown` DC-solves cleanly
        and silently, while a structurally singular one is loud.  A netlist
        the classifier cannot read is declined HERE in silence.

        On that fixture (a VCVS of gain `g` inside a C-V loop, index 2 off
        `g* = 1 + C2/C1` and index 3 on it) `topological_index` reads 1
        (provisional -- the VCVS is outside its covered class), and
        :func:`algebraic_conditioning` reads the algebraic block as SINGULAR
        at every gain, which is index >= 2.  The numeric one is right; the
        disagreement is the whole signal, and it is free here, since this
        path has already decided to decline.

        ⚠ NOT A BEHAVIOUR CHANGE, DELIBERATELY.  `x0_unknown` stays off: the
        remedy is justified for index 2 and NOT known to apply at index 3
        (`_resolve_x0_unknown`).  This says what was seen and names the
        explicit override; it does not take it.

        History: `doc/shooting_history.md`, `PSS._warn_if_the_block_disagrees`.
        """
        try:
            if idx is None or idx >= 2:
                return
            _sigma, ac = algebraic_conditioning(self.cir)
            if ac.get('verdict') != 'singular':
                return
            extra = ''
            if info.get('unclassified'):
                extra = (' %d element(s) are outside the classifier\'s covered '
                         'class (%s), so the topological reading is partial.'
                         % (len(info['unclassified']),
                            ', '.join(info['unclassified'][:4])))
            warnings.warn(
                'PSS: the topological index reads %s for this netlist, but its '
                'ALGEBRAIC BLOCK is numerically SINGULAR '
                '(sigma_min(d g_2/d y) = 0), which means index >= 2.%s A '
                'circuit like this solves cleanly and reports nothing, so the '
                'disagreement is the only signal you get. If it has a C-V loop '
                'or an L-I cutset through an element the classifier does not '
                'recognise, the manufactured opening step may be INCONSISTENT; '
                'pass x0_unknown=True explicitly to apply the index-2 remedy, '
                'having checked that the index really is 2 and not 3 -- the '
                'remedy is not known to apply at index 3.'
                % (idx, extra), RuntimeWarning, stacklevel=4)
        except Exception:                                      # noqa: BLE001
            ## ⚠ Best-effort, like everything else on this path: a DIAGNOSTIC
            ## that raises inside a defaulting helper would change which
            ## exception an invalid call reports, which is the defect the
            ## comment above records.  The gates call the helper DIRECTLY so a
            ## bug in it cannot hide behind this.
            return

    def _companion_reach(self):
        """How many charges back the chosen integrator's companion reads.

        The mechanistic property that decides whether this analysis needs
        the entering history as an unknown: a method reaching one charge
        back can be started from a single unknown, one reaching two
        cannot.  Asked of the
        integrator rather than inferred from `method`, so a fourth method
        arrives with the right answer instead of the default one.
        """
        ## Polymorphic: the method answers.  A stage method returns 1 (self
        ## starting, reads only x_n); an LMM computes it from its companion
        ## coefficients.  No isinstance/name branch to extend per method.
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        return integ.companion_reach()

    def _solves_history(self):
        """Whether the period map needs the entering history as an unknown.

        MEASURED, NOT ASSUMED (`benchmarks/pss_seam_cost.py`).  A method
        whose companion reads only `q_{n-1}` cannot see the fabricated
        opening history at all (euler's and trapezoidal's seams measure
        zero), so enlarging their system would fix nothing.  Gear-2 reads
        `q_{n-2}`, which in the plain formulation is the entering stand-in,
        and pays 54% of its total error there at 100 points/period, rising
        as the seam falls one order slower than the interior.

        Autonomous runs take it too, through the COMPOSED system (4c in the
        class docstring): a free period does not remove the need for a
        history the companion can read.

        History: `doc/shooting_history.md`, `PSS._solves_history`.
        """
        ## ⚠ TRAPEZOIDAL CANNOT JOIN THIS FORMULATION BY SOLVING FOR
        ## `x_{-1}`: a one-step companion reads ONLY `iq_{-1}`, and
        ## `d(iq_{-1})/d x_{-1} = -G` is SINGULAR at every purely reactive
        ## node, so the 2m x 2m system is rank-deficient.  (Solving for `iq`
        ## itself fails too; see the class docstring's 'plain' kind.)
        return self._map_kind() == 'pair'

    def _map_kind(self):
        """Which period map the PSS's method shoots on -- decided HERE, once;
        `solve`, `factored_period` and the state-event check all read it:

        * 'plain': a linear-multistep companion reaching one charge back --
          one entering unknown (`_walk_lmm`);
        * 'pair':  one reaching two (Gear-2) -- the solved history ``(x_0,
          x_{-1})`` (`_walk_lmm`; see `_solves_history` for why);
        * 'stage': a self-starting Runge-Kutta method (`_walk_stage`);
        * 'glm':   a Nordsieck multivalue method (`_glm_period_blocks`).

        Asked of the integrator, never inferred from the method's name.
        `_solves_history` is this decision's answer for 'pair', and this is
        the one place to override it: the tests that run Gear-2 on the old
        plain formulation do it here (`_force_plain_map`).

        History: `doc/shooting_history.md`, `PSS._map_kind`."""
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if integ.is_stage_method():
            return ('glm' if getattr(integ, 'is_multivalue', lambda: False)()
                    else 'stage')
        return 'pair' if integ.companion_reach() >= 2 else 'plain'


    def find_initial_solution(self, period, x0=None, npts=60, method=None,
                              eps_rel_lin=1e-2, eps_abs_lin=1e-3, n_iter=7,
                              max_periods=200, zeta=1e-7):
        """A PROPER initial solution to start shooting from, by pre-integrating
        until the fixed-point iteration has entered its LINEAR region.

        De Luca, Bolcato & Schilders, *Proper Initial Solution to Start Periodic
        Steady-State-Based Methods*, IEEE TCAS-I 2019 -- their Algorithm 2.  The
        paper is in `~/docs/07-shooting-methods/`.

        Shooting-Newton needs a start inside its contraction region, and the
        usual remedy is to GUESS a number of pre-integration periods; if the
        guess is wrong the solve diverges and the guess is repeated with no clue
        for the next one.  This detects the handoff point instead, from
        quantities the integration already produces.

        ⚠ WHAT THE CRITERION COMPARES -- AND WHAT IT IS NOT.  It is NOT "the
        iterate stopped moving" and NOT "a carried probe settled" (a settled
        probe only says the Jacobian stopped changing, which is equally true
        at an equilibrium; `benchmarks/pss_warm_start.py`).  The paper
        compares TWO SEQUENCES: the LINEAR prediction of the next shooting error
        against the one the ACTUAL nonlinear integration produces.  They agree
        only where the fixed-point map really has become linear, which is
        exactly the region a Newton-type method needs::

            u_k         = x_k - phi(x_k)                        (eq. 4)
            u_{k+1}     = J_phi(x_khat) u_k                     (eq. 12)
            utilde_{k+1} = x_{k+1} - phi(x_{k+1})               (eq. 13)

        accepted, componentwise, when (eq. 16)::

            |u_{k+1,j} - utilde_{k+1,j}| <= eps_rel_lin |u_{khat,j}| + eps_abs_lin

        holds for ``n_iter`` CONSECUTIVE iterations (the paper's defaults, used
        for every experiment in it: ``eps_rel_lin=1e-2``, ``eps_abs_lin=1e-3``,
        ``n_iter=7``).  The scale on the right is the shooting error at the
        DETECTION index ``khat``, not at the current one.  A failed check resets
        the run AND re-freezes ``J_phi`` at the new index, which is why ``khat``
        can move.

        ⚠ NON-AUTONOMOUS ONLY, AND THAT IS THE PAPER'S SCOPE, NOT AN OVERSIGHT
        HERE.  Its title, abstract and index terms all say non-autonomous, and
        the reason bites: for a FORCED circuit the DC point is not a fixed point
        of the period map, so "the map became linear" can only mean the orbit.
        For an AUTONOMOUS oscillator the equilibrium IS a fixed point of the
        period map and the map is linear in a neighbourhood of it, so this
        criterion will happily certify the trivial root.  **That case is NOT
        solved by this method and must not be handed to it.**

        ``J_phi(x_khat) u`` is taken by the paper's own alternative, the
        directional derivative of eq. (15),
        ``[phi(x + zeta u) - phi(x)] / zeta``, rather than by its Alg. 1
        left-product.  Both are in the paper; this one costs one extra period
        integration per iteration and buys freedom from the opening-frame
        question (which state a stored factorisation is the derivative *about*).

        Returns ``(x, info)``: the reduced state to start shooting from -- the
        paper's line 21, "the last computed x_k" -- and a dict with ``khat``,
        ``periods``, ``found`` and the per-period ``history``.  ``found=False``
        means ``max_periods`` ran out with no linear region; the returned ``x``
        is then simply the last iterate and carries no promise.

        History: `doc/shooting_history.md`, `PSS.find_initial_solution`.
        """
        T = float(period)
        n = self.cir.n
        iref = self.irefnode
        m = n - 1
        npts = int(npts)
        if npts < 1:
            raise ValueError('find_initial_solution: npts must be >= 1, got %r'
                             % (npts,))
        if n_iter < 1:
            raise ValueError('find_initial_solution: n_iter must be >= 1')
        times = np.linspace(0.0, T, npts + 1)
        hs = np.diff(times)
        integ = method if method is not None else getattr(self.par, 'method',
                                                          'euler')

        def phi(xr):
            """One period of the inner transient from `xr` -- the map the
            shooting residual is built on, run on a fixed uniform grid."""
            tr_saved = getattr(self, '_tran', None)
            self._tran = self._new_transient(self._integrator_for(integ))
            try:
                self._want_dfdh = False
                self._want_lte = False
                self._begin_period(np.asarray(xr, dtype=float))
                x = copy(np.asarray(xr, dtype=float))
                for j, t in enumerate(times[1:]):
                    x = copy(self.solve_timestep(x, t, hs[j]))
                return np.asarray(x, dtype=float).ravel()
            finally:
                self._tran = tr_saved

        if x0 is None:
            x = np.zeros(m)
        else:
            x = np.asarray(x0, dtype=float).ravel()
            if x.shape[0] == n:
                x = np.concatenate((x[:iref], x[iref + 1:]))

        phi_x = phi(x)
        k = 0
        i_iter = 0
        khat = 0
        x_khat = x.copy()
        phi_khat = phi_x.copy()
        u = x - phi_x
        u_khat = u.copy()
        history = []
        found = False
        while k < int(max_periods):
            if i_iter == 0:
                ## re-freeze the linear generator at the current index (Alg. 2
                ## lines 4-8): khat moves whenever the run of successes breaks.
                khat = k
                x_khat = x.copy()
                phi_khat = phi_x.copy()
                u = x - phi_x
                u_khat = u.copy()
            x_next = phi_x
            phi_next = phi(x_next)
            u_tilde_next = x_next - phi_next            ## eq. (13)
            ## eq. (12) via the paper's eq. (15) directional derivative, with
            ## the step scaled to the iterate so `zeta` is a RELATIVE size.
            un = float(np.linalg.norm(u))
            if un <= 0.0:
                u_next = np.zeros(m)
            else:
                z = zeta * max(float(np.linalg.norm(x_khat)), 1.0) / un
                u_next = (phi(x_khat + z * u) - phi_khat) / z
            gap = np.abs(u_next - u_tilde_next)
            ok = bool(np.all(gap <= eps_rel_lin * np.abs(u_khat) + eps_abs_lin))
            k += 1
            i_iter = i_iter + 1 if ok else 0
            history.append({'k': k, 'khat': khat, 'ok': ok, 'run': i_iter,
                            'gap': float(np.max(gap)),
                            'shooting_error': float(np.max(np.abs(u_tilde_next)))})
            x, phi_x, u = x_next, phi_next, u_next
            if i_iter >= int(n_iter):
                found = True
                break
        return x, {'khat': khat, 'periods': k, 'found': found,
                   'history': history}

    ## Set by `solve`; declared here so a caller may read them on a PSS that
    ## has not solved yet, and so `event_times` is never a stale leftover.
    break_events = False
    event_times = []

    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20, grid=None, matrix_free=False,
              x0_unknown=None, tstab=None, break_events=None,
              phase_rule='frozen', state_events=True, trace=False):
        """Solve for the periodic steady state.

        Returns an `InternalResultDict` with `tpss`, the orbit in time (the
        reference node reinserted), and `fpss`, its spectrum -- ⚠ RMS,
        energy-folded, positive frequency: multiply by `sqrt(2)` to compare
        against a PEAK-convention spectrum.  Among what it sets on `self`:
        `converged` (⚠ the waveform is returned either way -- it is the last
        iterate when this is False), `period` (the solved period on an
        autonomous circuit), `spectral_radius` and `floquet_multipliers`
        (None after a matrix-free solve), `max_lte` / `total_lte` /
        `max_lte_seam` (`_report_lte`), `fundamental_period` (set when an
        autonomous solve returned a multiple), `autonomous`,
        `solved_history`, `waveform`, `tstab_state`, `shooting_residual`
        and `shooting_trace`.

        `period` is the period of a DRIVEN circuit and the SEED of an
        AUTONOMOUS one (nothing in it depends on `t`), whose period joins
        the unknowns with a phase condition.  ⚠ Seed an oscillator at or
        above its expected period: `k*T` solves the periodicity condition
        whenever `T` does and the solve follows its seed
        (`_check_fundamental` warns), and `T = 0` is a regular root that a
        seed below the fundamental is drawn to (`_free_period_solve`).

        `x0` is the seed state (reference node excluded); `None` shoots from
        zeros on a driven circuit and from the operating point on an
        autonomous one.  `timestep` sets the uniform grid
        (`int(period / timestep)` points) unless `grid` is given.
        `maxiterations` bounds the shooting Newton.

        `grid` is RECORDED SCOPE ITEM 5: a sequence of step FRACTIONS of the
        period, summing to 1, used in place of the uniform `timestep` grid.
        Fractions rather than absolute times because an autonomous period is
        an unknown and every step has to scale with it.  See
        `_period_grid`, and `lte_grid` for a grid chosen by an adaptive
        transient; the grid is still frozen for the whole solve, so the
        shooting Newton stays exact.

        `break_events` lands the circuit's source discontinuities on grid
        points (`event_grid`).  `None` -- the default -- turns it on for
        every method (`_resolve_break_events`).  A circuit whose sources
        declare no discontinuity is untouched bit-for-bit either way.

        `state_events` (default True): STATE EVENTS AS NEWTON UNKNOWNS.  A
        comparator's, a threshold switch's or a latch's edge happens where
        the SOLUTION crosses a condition -- `Circuit.state_events()`,
        `row . x = threshold` -- at a time no frozen grid can hold, and with
        the crossing inside a step every method is FIRST order.  With
        `state_events` a converged solve runs a second stage
        (`_state_event_stage`): the crossings of the first stage's orbit
        become unknowns `theta_k` (fractions of the period) alongside `x_0`,
        the grid between consecutive events scales with its segment (a
        proportional column per event), each event contributes the row
        `row . x(theta_k T) = threshold`, and the bordered Newton lands the
        grid on the crossing exactly.  The solved fractions become the grid
        every consumer replays on and the event nodes are breaks for the
        period quadrature.  It runs on every kind -- the stage methods,
        gear's pair, a Nordsieck GLM (its map on the state) and the plain map
        (euler, trap, theta) OPENED AT `x(0)` (`x0_unknown`, its default when
        the circuit declares state events) -- driven or free period: on an
        AUTONOMOUS circuit `z = [x_0, theta, T]` (gear's `x_0` a pair), the
        period one more column of the event algebra (`d h_j / d T =
        fraction_j`), the phase row closing the system, the polish
        convention proportional.
        A first stage that did not converge is not the end: unless it
        collapsed onto a trivial root, the stage runs from its last iterate
        and its own Newton decides `converged` --
        `stage_one_converged` records whether the first stage did (None
        when no stage ran).  Radau is the method for these circuits: gear's
        and trbdf2's own second-order error on the switch-off decay
        dominates.  The stage is for crossings sharper than the grid.
        `state_events=False` keeps the one-stage solve.

        `tstab` runs a TRANSIENT for that many seconds before shooting and
        uses its final state as the seed -- the stabilisation time every
        commercial PSS offers, and the remedy for the failure this analysis
        fails most often: a seed in the trivial-root basin near the DC
        point.  `None` (the default) shoots from `x0`, or from the operating
        point.  The periods needed are a property of how strongly the limit
        cycle attracts, not of how bad the seed is (the `1/mu`
        amplitude-envelope constant: van der Pol needs one at mu = 1 and ~24
        at mu = 0.05) -- time constants much longer than the period need
        proportionally more (Kundert, *Introduction to RF Simulation*; cited,
        not verified here).
        ⚠ THE STOPPING POINT IS THE CALLER'S, DELIBERATELY: near the DC
        point the trivial root IS a fixed point of the period map and passes
        every periodicity test, so no automatic handoff criterion measured
        here identifies it on an oscillator (`benchmarks/pss_warm_start.py`;
        for a driven circuit see `find_initial_solution`).
        ⚠ IT CANNOT ESCAPE AN EQUILIBRIUM IT IS STARTED ON: with `x0=None`
        an autonomous seed is the operating point, which a transient never
        leaves -- pass an `x0` off the equilibrium, or an `ic` on a device.

        `x0_unknown` solves for `x_0` itself instead of for `x_in`, the
        pre-image of the plain path's manufactured opening Euler step.  The
        plain path's Jacobian is taken with respect to `x_0` while its
        unknown is `x_in` -- a frame error; the true `dF/dx_in` is SINGULAR
        -- so its iteration is a CONTRACTION with a linear rate rather than
        a Newton.  With this set there is no manufacturing step, the unknown
        is the period's own start, and the Jacobian is exact.  `None` (the
        default) DECIDES FROM THE TOPOLOGY: on for a netlist the index
        criterion proves is index 2, where the manufactured opening step is
        inconsistent, off otherwise; an explicit `True` or `False` is
        honoured untouched (`_resolve_x0_unknown`).
        ⚠ IT IS NOT FREE, AND THE TRADE-OFF IS GRID-DEPENDENT.  Trapezoidal
        still needs an L-stable opener (without one its period map is
        singular at EVEN step counts on every MNA circuit), so the Euler
        step moves INSIDE the period, where it degrades the ORBIT (Q=20
        resonator, 20 V analytic peak, 100 points: 20.01273 without, 19.76939
        with; the gap closes as the one first-order step is diluted).  On a
        stiff circuit on a non-uniform grid it is the other way round (van
        der Pol on its LTE grid: -47.3 ppm against -73.8, quadratic at even
        and odd point counts alike).  So: uniform grid on a benign circuit,
        leave it off; non-uniform grid on a stiff one, turn it on
        (`benchmarks/pss_x0_unknown.py`).  Euler is unchanged either way --
        its manufacturing step IS an Euler step.  A two-step (solved-history)
        method refuses it: it already solves for `x_0`.

        `matrix_free` is RECORDED SCOPE ITEM 6: solve the outer system
        without ever forming the monodromy, propagating ONE vector per
        Krylov iteration instead of `2m` columns per step.  Worth asking for
        on LARGE circuits only -- against a dense-solver dense path it loses
        at m=40 and wins above roughly m=250 (1.36x at m=242, 2.13x at
        m=1002 on the driven solved-history system; systems with `m` rather
        than `2m` columns gain less).  ⚠ It buys that with memory, `2 N m^2`
        doubles of stored factorisations (~800 MB at m=1002, 50 points), and
        it does NOT produce a monodromy, so `spectral_radius` is `None`
        after a matrix-free solve.  Built for every kind, driven or free
        period (a Nordsieck GLM through its map on the state, seeded by the
        linearised startup).  ⚠ The
        state-event stage needs the dense map and does not run under
        `matrix_free`: a circuit that declares state events warns.  ⚠ Those
        figures were taken with a DENSE linear solver on both sides; with a
        sparse one the m~250 gate may move --
        `benchmarks/pss_matrix_free_sparse.py` is the harness, to be run
        quiet before quoting any of this as the sparse answer.

        `phase_rule` chooses how an AUTONOMOUS solve removes the orbit's
        time-translation freedom.  `'frozen'` (the DEFAULT) picks `k` and
        its pinned value once at the seed.  `'reselect'` is Aprille &
        Trick's own Step 3: at every Newton iterate the pinned coordinate
        is `k = argmax |dphi/dT|` -- the vector field at the period's end --
        and it is held at the iterate's OWN value, so the Newton step moves
        every other coordinate and the period (`dz[k] = 0`).  It is built
        inside the residual, as the bordered row with a zero residual and a
        pure function of the iterate, so the damped Newton's carried trial
        evaluation stays consistent.  A driven circuit ignores it.

        ⚠ RESELECT IS OPT-IN, AND THE REASON IS MEASURED, NOT CAUTION.  Ask
        for it when a seed is far from the orbit and the default reports
        non-convergence (`tstab` is the other remedy and they compose): a
        frozen pin is a VALUE the orbit must attain, and a far seed names
        one it never reaches (van der Pol at 4x the orbit amplitude: 0/6
        seeds frozen, 6/6 reselected).  Matrix-free keeps less of the gain,
        and far seeds can still reach the `T = 0` root.  What it costs:

          * **A discontinuous period map stops converging.** `Idtmod` with
            the wrap landing exactly ON a grid point solves frozen and fails
            reselected: the moving pin wanders on a map with no derivative.
          * **It lands on a different PHASE of the same orbit**, so any
            surface whose value depends on where `t = 0` sits on the orbit
            (e.g. `frequency_aware_ppv`'s mode content) moves with it --
            and such a surface is phase-specific under the frozen rule too.

        Dense and matrix-free Newtons pick `k` from `dphi/dT` values that
        differ in the last ulp, so their `lambda_2` agree to ~1e-13 rather
        than to the bit.  Where the converged period differs between the
        rules it is the GRID, not the rule (trap with `x0_unknown=True` has
        a start-phase-dependent discrete orbit).

        `trace` records the shooting Newton: every evaluation of its
        residual, line-search trials included, as `(z, F, J)` in
        `shooting_trace` (`J` is None on the matrix-free path; a closing
        second pass appends to the same list).  `shooting_residual` is always
        kept -- the function the solve drove to zero, returning `(F, J)` --
        so a caller can re-evaluate it, e.g. for a finite-difference check of
        `J`.

        History: `doc/shooting_history.md`, `PSS.solve`.
        """
        ## ONE SHOOTING SOLVE, IN ITS PHASES.  Each phase is a method with its
        ## own record; they share the run's state through `run`.
        run = self._solve_prepare(refnode, period, x0, timestep,
                                  maxiterations, grid, matrix_free,
                                  x0_unknown, tstab, break_events,
                                  phase_rule, state_events, trace)
        self._shoot(run)
        self._report_convergence(run)
        X, walk, lte_seen = self._replay_orbit(run)
        self._report_lte(run, lte_seen)
        self._check_fundamental(X, walk, run.period)
        tpss, fpss = self._orbit_results(run, X)
        polished = self._closing_polish(run)
        if polished is not None:
            return polished
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

    def _solve_prepare(self, refnode, period, x0, timestep, maxiterations,
                       grid, matrix_free, x0_unknown, tstab, break_events,
                       phase_rule, state_events, trace=False):
        """`solve`, phase 1: validate, build the grid, decide autonomy and
        the period column, seed (operating point, `tstab`), pin the phase,
        choose the formulation.  Returns the run's state."""
        self._solve_kwargs = dict(refnode=refnode, maxiterations=maxiterations,
                                  matrix_free=matrix_free, tstab=tstab,
                                  period_seed=float(period))
        self._monodromy_twin = None
        self._twins = {}
        ## ⚠ HIDDEN STATE IS REFUSED, NOT INTEGRATED AND HOPED OVER.
        ## `TLine.history` is filled by `cir.accept_step`, which the
        ## TRANSIENT calls at every accepted step and which this analysis
        ## never calls -- PSS drives `solve_timestep` directly.  With the
        ## buffer empty `TLine.G`/`u` stamp a DC SHORT, so the line is
        ## silently absent and the solve reports `converged` with a wrong
        ## answer.  Filling the buffer is not the fix either: `phi` would
        ## become history-dependent and the monodromy the derivative of a
        ## neighbouring problem, which `_begin_period` exists to prevent (it
        ## resets what is IN `x`; no reset can fix a period map that is not
        ## a function of `x_0`).  The refusal is per ELEMENT
        ## (`Circuit.hidden_state` defaults False), not per class: a
        ## distributed component with a known frequency-domain description is
        ## tractable in other analyses (Yang & Phillips, DAC 2002; cited, not
        ## verified here), and an
        ## element that declared its state properly would pass.
        ## History: `doc/shooting_history.md`, `PSS._solve_prepare`.
        _hidden = self.cir.hidden_state_elements()
        if _hidden:
            raise NotImplementedError(
                'PSS: these elements carry HIDDEN STATE -- %s -- so THIS '
                'formulation cannot solve this circuit. The period map must '
                'be a function of x_0 alone, and they stamp from state that '
                'lives outside x and is filled by accept_step, which only a '
                'forward transient calls. Left alone the answer would be '
                'silently wrong rather than slow: an empty TLine history '
                'stamps the line as a DC short and the solve reports '
                'converged. Use Transient for this circuit. ⚠ This is a '
                'limit of the ELEMENT as implemented here, not of shooting '
                'or of distributed components as a class: a component with '
                'a KNOWN frequency-domain description (a transmission line, '
                'an S-parameter block) is tractable in a time-domain '
                'steady-state solve by either admitting the delay state '
                'into the unknowns, or applying the component spectrally -- '
                'the Fourier transform diagonalises the convolution, so its '
                'action becomes a multiply by Y/Z/S while the state stays '
                'finite. Both are different analyses than this one.'
                % ', '.join(sorted(_hidden)))

        self.period = period
        toolkit = self.toolkit

        ## ⚠ ONE REFERENCE NODE PER ANALYSIS, CHECKED.  `self.irefnode` is
        ## fixed in `__init__` and is what the TRAVERSAL eliminates; this
        ## local one comes from `solve`'s own `refnode=` and is what
        ## reinserts the zero row into the RESULT.  If they differ the rows
        ## are incoherent (ground itself comes back non-zero) and there is no
        ## answer to give, so it is refused rather than silently rotated.
        irefnode = self.cir.get_node_index(refnode)
        if irefnode != self.irefnode:
            raise ValueError(
                'PSS: solve(refnode=...) names a different reference node '
                '(index %d) than the analysis was constructed with (index '
                '%d). The traversal eliminated one and the result would '
                'reinsert the other, so the waveform would be reported '
                'against a node the solve never used -- ground itself comes '
                'back non-zero. Pass the same node to both, or construct '
                'the analysis with PSS(cir, irefnode=...) and leave '
                "solve()'s refnode at its default."
                % (irefnode, self.irefnode))
        ## ⚠ CLEARED BEFORE THE RUN, not after it.  These describe the
        ## period this call is about to solve for; leaving the previous
        ## call's behind would let `factored_period()` hand back an operator
        ## for the LAST solve after this one failed, and `converged` alone
        ## would not catch it.
        self._period_state = None
        self._factored_period_cache = None
        self.waveform = None

        ## ⚠ `theta`'s bias is PER-PERIOD (`_theta_biased`), so it must be
        ## known before anything builds or reuses the inner transient.
        ## `_new_transient` only runs when there is no cache, so a SECOND
        ## `solve()` at a different period would otherwise silently keep the
        ## first one's bias -- re-bias the cached integrator in place rather
        ## than dropping the cache, which would rebuild a `Transient` per
        ## solve for every method that does not care.  A no-op for all of
        ## them (`_theta_biased` type-checks, and tolerates `None`).
        self._theta_period = float(period)
        _tr_cached = getattr(self, '_tran', None)
        if _tr_cached is not None:
            self._theta_biased(getattr(_tr_cached.par, 'integrator', None))

        ## Everything `grid_error` needs to repeat THIS solve on a finer grid.
        ## Recorded rather than re-derived so the refinement differs from the
        ## original in the timestep and in nothing else.
        self._solve_args = dict(refnode=refnode, period=period, x0=x0,
                                timestep=timestep,
                                maxiterations=maxiterations, grid=grid,
                                matrix_free=matrix_free,
                                x0_unknown=x0_unknown, tstab=tstab,
                                break_events=break_events,
                                phase_rule=phase_rule)

        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = toolkit.zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !


        #create vector with timepoints and a more fitting dt
        ## ⚠ the flag must be set BEFORE the grid is built, because
        ## `_period_grid` consults it to decide whether to subdivide a
        ## coarse opening step -- see the note there.
        ## ⚠ `None` means "decide from the topology" -- see
        ## `_resolve_x0_unknown`.  Resolved to a concrete bool HERE, before
        ## anything reads it, so every downstream use sees one value.
        ## ⚠ STATE EVENTS ON A PLAIN MAP OPEN IT AT `x(0)`.  The event stage
        ## needs the map opened there: the manufactured opener's own step
        ## moves with the first crossing, and that dependence is not carried.
        ## Measured, trap on the comparator oscillator: with the opener the
        ## stage failed to converge at 200 points and read -1.2e-3 at 800
        ## (unstaged -7.2e-4); opened at `x(0)`, +1.4e-4 / +9.1e-6 against
        ## +6.4e-3 / -7.2e-4 unstaged.  So a plain-map run with declared
        ## state events defaults `x0_unknown` to True (an explicit value is
        ## honoured, and False warns that the stage does not run).
        _se_rows = (self.cir.state_events()
                    if state_events and hasattr(self.cir, 'state_events') else [])
        if (x0_unknown is None and _se_rows and not matrix_free
                and self._map_kind() == 'plain'):
            x0_unknown = True
        x0_unknown = self._resolve_x0_unknown(x0_unknown)
        if self._integrator_for(getattr(self.par, 'method', 'euler')
                                ).needs_x0_unknown():
            ## Self-starting stage methods: `x_in` IS `x_0`, there is no
            ## manufacturing step to differentiate `x_0` back through, so the
            ## unknown is always `x_0` itself.  The method says so
            ## (`needs_x0_unknown`), which keeps the phase pin and every
            ## open-at-x0 branch consistent without a name check here.
            x0_unknown = True
        self._open_at_x0 = bool(x0_unknown)
        ## Break the traversal's steps at the source discontinuities -- see
        ## `_resolve_break_events`, and `event_grid` for the snap that keeps
        ## it from manufacturing slivers.
        self.break_events = self._resolve_break_events(break_events)
        if self.break_events:
            _ev = (self.event_grid(period, grid=grid) if grid is not None
                   else self.event_grid(period, npts=int(period / dt)))
            ## ⚠ ONLY replace the grid when there ARE events.  `event_grid`
            ## rebuilds a uniform grid from `linspace` even when it finds none,
            ## and that differs from `_period_grid`'s own in the last bit --
            ## enough to move every event-free solve in the suite for nothing.
            ## Touching the grid only when an event exists keeps every circuit
            ## without one BIT-IDENTICAL, the same guarantee `_fold_periodic`
            ## gives a circuit with no periodic state.
            if self.event_times:
                grid = _ev
        times, hs = self._period_grid(period, int(period / dt), grid)
        npts = len(times)
        self._grid_fracs = (None if grid is None
                            else np.asarray(grid, dtype=float))
        ## read by `_period_grid`, which is called from the residual
        ## closures and so cannot take it as an argument
        self._open_at_x0 = bool(x0_unknown)
        ## The fold gauge, collected once per solve (late-bound moduli are
        ## resolved by now).  See `_fold_periodic` for why the residual needs
        ## it and the Jacobian does not.
        self._periodic_fold = self._collect_periodic_fold()
        self._wrap_jump_warned = False
        alpha = 1

        ## AUTONOMY IS DECIDED BEFORE THE SOLVE, because it decides which
        ## system is solved.  Structural and exact -- see `_is_autonomous`.
        self.autonomous = self._is_autonomous(times)
        ## the state-event stage runs on every kind but a matrix-free solve
        ## and a plain map not opened at `x(0)` -- see the docstring; say so
        ## once when the circuit declares events
        self._state_event_fracs = None
        self._event_columns = None
        if state_events:
            _rows = self.cir.state_events() if hasattr(self.cir, 'state_events') else []
            _method_se = getattr(self.par, 'method', 'euler')
            if _rows and matrix_free:
                warnings.warn(
                    'PSS: this circuit declares %d state event(s) (a threshold '
                    'switch or comparator), but the state-event stage needs '
                    'the dense period map and does not run under '
                    'matrix_free=True: the crossings stay inside their steps '
                    'and the solve is first order there. Drop matrix_free to '
                    'land them, or pass state_events=False to silence this.'
                    % len(_rows), RuntimeWarning, stacklevel=3)
            elif (_rows and self._map_kind() == 'plain'
                  and not self._open_at_x0):
                warnings.warn(
                    'PSS: this circuit declares %d state event(s), and on the '
                    'plain map (method %r) the state-event stage needs the '
                    'map opened at x(0): with x0_unknown=False the crossings '
                    'stay inside their steps and the solve is first order '
                    'there. Pass x0_unknown=True (the default with state '
                    "events), or use method='radau'." % (len(_rows), _method_se),
                    RuntimeWarning, stacklevel=3)
        ## the period-column convention for this solve (see the Parameter)
        _pc = str(getattr(self, '_force_period_column', None)
                  or getattr(self.par, 'period_column', 'auto'))
        if _pc not in ('auto', 'proportional', 'closing'):
            raise ValueError("period_column must be 'auto', 'proportional' "
                             "or 'closing', got %r" % (_pc,))
        ## ⚠ 'auto' IS 'closing' WITH THE PROPORTIONAL POLISH, ON A CALLER'S
        ## GRID FOR AN AUTONOMOUS RUN.  Closing keeps its basin (it converges
        ## from seeds several percent off where proportional's per-step
        ## Newton fails), and the unconditional polish (`_closing_polish`)
        ## makes the answer proportional's at the solved period,
        ## seed-independent.  'proportional' stays selectable by name -- one
        ## solve, the seed's own grid.  On a UNIFORM grid 'auto' is
        ## proportional: the closing column there makes one step of N absorb
        ## the period correction for no reason, and uniform-grid gear solves
        ## that converge proportionally did not converge closing.
        self._period_column = ('closing' if (_pc == 'closing' or
                               (_pc == 'auto' and grid is not None
                                and getattr(self, 'autonomous', False)))
                               else 'proportional')
        self._closing_inner = None
        self._closing_warned = False
        phase_k, phase_pin = 0, 0.0
        if phase_rule not in ('reselect', 'frozen'):
            raise ValueError("phase_rule must be 'reselect' or 'frozen', not %r"
                             % (phase_rule,))
        self.phase_rule = phase_rule

        if self.autonomous:
            ## An unseeded autonomous run starts at the origin, which IS a
            ## periodic solution -- the trivial one -- and the free-period
            ## system would sit there just as contentedly as the fixed one
            ## did.  The operating point is the honest default: for a phase
            ## accumulator `ic` pins it on the orbit.
            ##
            ## ⚠ `tstab` RUNS AFTER THIS, NOT BEFORE, and the order is the
            ## whole of what makes it work.  A pre-integration seeded from
            ## `zeros` starts AT the equilibrium of an autonomous circuit and
            ## a transient from an exact equilibrium never leaves it, so the
            ## warm start would return the basin it was asked to escape.
            ## Starting it from the operating point is the honest version of
            ## the same statement.
            if x0 is None:
                from pycircuit.circuit.dcanalysis import DC
                xdc = np.asarray(DC(self.cir, toolkit=self.toolkit).solve().x,
                                 dtype=float).reshape(-1)
                x = np.concatenate((xdc[:irefnode], xdc[irefnode + 1:]))

        if tstab:
            ## ⚠ THE PRE-INTEGRATION IS A PLAIN TRANSIENT, and it has to be:
            ## its whole value is that it is NOT a shooting solve, so it
            ## cannot be captured by the basin that the shooting Newton is
            ## stuck in.  It runs on its own adaptive grid -- `timestep` is
            ## a first step, not an imposed one -- because nothing here
            ## needs `phi` to be a function of `x_0`; that requirement
            ## starts when the shooting does.
            from pycircuit.circuit.transient import Transient
            _xred = np.asarray(x, dtype=float).reshape(-1)
            _xfull = np.concatenate((_xred[:irefnode],
                                     np.zeros(1), _xred[irefnode:]))
            _pre = Transient(
                self.cir, toolkit=self.toolkit, reltol=self.par.reltol,
                iabstol=self.par.iabstol, vabstol=self.par.vabstol,
                nrsolver=self.par.nrsolver,
                linearsolver=self.par.linearsolver, scaler=self.par.scaler)
            _res = _pre.solve(refnode=refnode, tend=float(tstab),
                              timestep=dt, x0=_xfull)
            _last = np.asarray(_res.x, dtype=float)[:, -1]
            x = np.concatenate((_last[:irefnode], _last[irefnode + 1:]))
            self.tstab_state = x

        if self.autonomous:
            ## THE PHASE CONDITION pins the coordinate moving FASTEST at the
            ## seed, so the orbit crosses the pinning hyperplane
            ## transversally.  Pin a slow one and the last row of the
            ## bordered Jacobian is nearly parallel to the null direction it
            ## exists to remove, which is a singular system wearing an extra
            ## equation.
            ##
            ## ⚠ The row removes the singularity from the UNIT multiplier
            ## only, so an oscillator whose other multipliers cluster near 1
            ## gives a bordered system that is nonsingular and ill
            ## conditioned (the high-Q note in the class docstring).
            ##
            ## ⚠ THE `argmax` COMPARES VOLTS WITH AMPERES ON PURPOSE: the row
            ## can only remove the orbit's tangent in proportion to
            ## `|e_k . fhat|`, and one step of `|dx_k|` IS `|f_k|` up to `h`,
            ## so this argmax maximises the quantity the row needs.  The
            ## scaling that lets a large coordinate win the argmax is the
            ## same scaling that makes it dominate `f`; the two cancel.
            ## Pinned by `test_the_phase_pin_compares_units_on_purpose`.
            ## (A swing-normalised pin and an orthogonality (Poincare) row
            ## were both measured and are not better here; see history.)
            ##
            ## `phase_pin` is a VALUE the orbit must attain, so a seed far off
            ## the orbit can pin one outside its range and the system is then
            ## INCONSISTENT rather than merely hard, reporting ordinary
            ## non-convergence.  That is about SEEDS (`tstab`,
            ## `phase_rule='reselect'`), not the formulation.
            self._begin_period(x)
            _x1 = self.solve_timestep(x, times[0], hs[0])
            _x2 = self.solve_timestep(_x1, times[1], hs[0], iq_last=self._iq)
            phase_k = int(np.argmax(np.abs(np.asarray(_x2) - np.asarray(_x1))))
            self.phase_k = phase_k                  # which coordinate is pinned, for diagnosis
            ## The rule is Aprille & Trick's (oscillator paper, Step 3:
            ## select k by max |f_k|), kept raw on measurement.  Their Step 3
            ## sits INSIDE the iteration, re-choosing `k` and pinning the
            ## iterate's own value: that is `phase_rule='reselect'`
            ## (`_phase_row`), OPT-IN -- gains and costs in `solve`'s
            ## docstring.
            ##
            ## ⚠ THE PHASE ROW SITS OUTSIDE THE INTEGRATOR ON PURPOSE: it
            ## augments the OUTER shooting system, so it cannot raise the
            ## index of the DAE actually being integrated (Brachtendorf et
            ## al., TCAD 33(6) 867-878, warn that adding an algebraic
            ## equation to the integrated system does).
            ##
            ## ⚠ THE PIN MUST BE IN THE UNKNOWN'S OWN FRAME.  `_x1` is the
            ## state one step AFTER the seed, which is the right thing to
            ## pin when the unknown is `x_in` and `x_0` is manufactured from
            ## it -- and the wrong thing when the unknown IS `x_0`.  With a
            ## fine opening step the two are nearly equal and the mismatch
            ## hides; on a coarse one it pins a value the orbit need never
            ## attain and the solve dies with a bare non-convergence.
            ## History: `doc/shooting_history.md`, `PSS._solve_prepare`.
            phase_pin = float(np.asarray(x if x0_unknown else _x1)[phase_k])

        ## Resolved here as well as in `solve_timestep`, because the SHOOTING
        ## Jacobian depends on which integrator the inner steps used.
        ##
        ## ⚠ THE NAME IS VALIDATED BEFORE ANYTHING ASKS THE INTEGRATOR A
        ## QUESTION, so an unknown name raises this `ValueError` rather than
        ## a `KeyError` from several frames down (`_solves_history`).
        method = getattr(self.par, 'method', 'euler')
        if method not in ('euler', 'trap', 'trapezoidal', 'theta', 'gear',
                          'gear2', 'trbdf2', 'radau', 'esdirk43',
                          'glm2', 'glm3', 'glm4'):
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau', 'esdirk43', 'glm2', 'glm3' or 'glm4', not %r"
                % (method,))

        ## Whether the entering history joins the unknowns.  Decided once,
        ## here, because it chooses which system is solved -- like autonomy.
        solved_history = self._solves_history()
        self.solved_history = solved_history
        xm1_ss = None

        return _SolveRun(
            refnode=refnode, period=period, x=x, dt=dt,
            maxiterations=maxiterations, matrix_free=matrix_free,
            x0_unknown=x0_unknown, phase_rule=phase_rule,
            state_events=state_events, irefnode=irefnode, n=n, times=times,
            hs=hs, npts=npts, alpha=alpha, phase_k=phase_k,
            phase_pin=phase_pin, method=method,
            solved_history=solved_history, xm1_ss=xm1_ss, trace=trace)

    def _shoot(self, run):
        """`solve`, phase 2: THE SHOOTING NEWTON -- the fixed-period or
        free-period system on the method's period map (dense or matrix-free),
        then the state-event stage.  Leaves the solution in `run`."""
        toolkit = self.toolkit
        (n, x, period, times, hs, npts, alpha, irefnode) = (
            run.n, run.x, run.period, run.times, run.hs, run.npts, run.alpha,
            run.irefnode)
        (x0_unknown, solved_history, method, phase_rule, phase_k,
         phase_pin) = (run.x0_unknown, run.solved_history, run.method,
                       run.phase_rule, run.phase_k, run.phase_pin)
        (maxiterations, matrix_free, state_events, xm1_ss) = (
            run.maxiterations, run.matrix_free, run.state_events, run.xm1_ss)

        def _phase_row(x0_vec, tcol):
            """The autonomous phase row at THIS iterate: `(k, residual)`.

            `'reselect'` pins `k = argmax |dphi/dT|` over the `x_0` block at
            the iterate's own value, so the residual is zero and the row
            only fixes the step (`dz[k] = 0`); `'frozen'` compares the seed's
            `k` against the seed's value.  See `solve`'s docstring."""
            if phase_rule == 'reselect':
                tc = np.abs(np.asarray(tcol, dtype=float).ravel()[:n - 1])
                return int(np.argmax(tc)), 0.0
            return (phase_k,
                    float(np.asarray(x0_vec, dtype=float)[phase_k]) - phase_pin)

        def _closing(x0, x_end, M, tms_):
            """The fixed-period system at one iterate: ``F = x_0 - phi(x_0)``
            (folded on the idtmod rows, see `_close_periodic`), ``J = I -
            alpha M``."""
            F = self._close_periodic(x0, x_end, tms_)
            D = np.asarray(toolkit.eye(F.shape[0]))
            return F, D - alpha * M

        def _bordered(F, J, tcol, x0_vec, phase_col):
            """The FREE-PERIOD system: the fixed-period one `(F, J)` bordered
            by the period column and the phase row --

                F = [ F ,  x0[k] - pinned ]
                J = [[ J , -dphi/dT ],
                     [ e_k^T ,  0   ]]

            -- because without a phase condition the system is singular by
            construction (every point on the orbit is a solution, so `I - M`
            has a null direction along it).  `k` is chosen by `_phase_row`
            from `phase_col`; the row pins only the `x_0` block, whatever the
            width of `F` (gear's pair: one phase row still suffices)."""
            w = len(F)
            Jb = np.zeros((w + 1, w + 1))
            Jb[:w, :w] = J
            Jb[:w, w] = -np.asarray(tcol).ravel()
            _k, _r = _phase_row(x0_vec, phase_col)
            Jb[w, _k] = 1.0
            Fb = np.zeros(w + 1)
            Fb[:w] = F
            Fb[w] = _r
            return Fb, Jb
        ## THE SHOOTING JACOBIAN FOLLOWS THE INTEGRATOR'S OWN COEFFICIENTS
        ## (the general recursion is in the class docstring).
        ##
        ## Backward Euler's per-step sensitivity is
        ##     dx_n/dx_{n-1} = Jf_n^-1 * C(x_{n-1})/h
        ## -- the COMPANION CONDUCTANCE at the previous point, not the raw
        ## capacitance matrix.  Without the `/h`, a singular C collapses the
        ## accumulated product to EXACTLY ZERO, the Jacobian is `I`, and the
        ## "shooting Newton" is silently successive substitution.
        ##
        ## TRAPEZOIDAL'S MONODROMY MUST CARRY `iq` AS WELL AS `x`: its
        ## period map is a function of (x, iq), and an x-only monodromy
        ## converges SLOWER than no Jacobian at all.  Differentiating
        ##
        ##     iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}
        ##     0    = i(x_n) + iq_n + u(t_n)
        ##
        ## gives a propagation of `d(x,iq)/dx0` that costs one extra matrix
        ## product over the Euler form:
        ##
        ##     rhs   = Geq_{n-1} Px + Pq
        ##     Px_n  = Jf_n^-1 rhs
        ##     Pq_n  = Geq_n Px_n - rhs
        ##
        ## Euler is the SAME recursion with `Pq == 0`: one formula, two
        ## methods, not a second code path.
        ## History: `doc/shooting_history.md`, `PSS._shoot`.

        ## THE PERIOD MAP, per kind -- the only thing about the Newton that
        ## depends on the method.
        _kind = self._map_kind()

        def _pmap(z, T, tms_, hs_, want_dT):
            """One period from the unknown `z`: ``(z_0, z_end, M, Mt)`` with
            `M = dz_end/dz_0` and, on request, the period column `Mt`.

            * PLAIN (one-step LMM): `z` is the entering state; with a
              manufacturing step the period opens one step in, so `z_0` is
              the state it opened at.
            * gear's PAIR: `z = (x_0, x_{-1})`, and BOTH close --
              ``F = [x_0 - x_{N-1}, x_{-1} - x_{N-2}]``, the rows of
              ``M = [[A(N-1,0), A(N-1,-1)], [A(N-2,0), A(N-2,-1)]]`` from
              the pair walk (`_walk_lmm`).  A two-step companion needs two
              states to be continued, so periodicity of ONE is an
              under-determined statement about the orbit.  ⚠ THE HISTORY
              POINT MOVES WITH T: `x_{-1}` sits at `-T/(N-1)`, and `x_{-1}`,
              `x_{N-2}` are the same phase of the orbit at every `T`, so the
              residual is still right and its `T` column is the propagation
              to step N-2.
            * STAGE (Radau, TR-BDF2, ESDIRK): self-starting, `z` IS `x_0`,
              `M` the dense stage product (`_walk_stage`).
            * GLM: the method's own map (the startup at the top of the
              period, then N multivalue steps); ⚠ `M` is APPROXIMATE -- the
              residual is exact, the Jacobian drops the startup's derivative
              (see `_walk_glm`); its period column carries the two
              explicit `T` dependences a multivalue method has.

            ⚠ THE PERIOD COLUMN IS TRACTABLE ONLY FOR AN AUTONOMOUS CIRCUIT:
            the grid is rebuilt at the current `T` (``dh/dT = h/T`` for every
            step, uniform or not) and the stage derivatives carry no time of
            their own."""
            w = self._walk(_kind, z, tms_, hs_, T=T, want_dT=want_dT,
                           open_at_x0=x0_unknown)
            M = w.monodromy()
            ## kept for the checks after the solve: the spectrum is the only
            ## place a free period announces itself.  (A GLM's too since
            ## 2026-09-24: its map on `x` is exact once the startup is
            ## linearised -- see `_walk_glm`.)
            self._monodromy = M
            return (w.z0(), w.end(), M,
                    w.period_column() if want_dT else None)

        def func(z):
            """The fixed-period system, ``x_0 - phi(x_0) = 0``."""
            z0_, z_end, M, _Mt = _pmap(z, period, times, hs, False)
            return _closing(z0_, z_end, M, times)

        def func_autonomous(zT):
            """The FREE-PERIOD system: unknowns `(z, T)`, the fixed-period
            equations bordered by the period column and the phase row (see
            `_bordered`) -- rebuilt at the CURRENT `T`, which is what keeps
            `dh/dT = h/T` true of every step."""
            z, T = zT[:-1], float(zT[-1])
            tms_, hs_T = self._period_grid(T, npts, self._grid_fracs)
            z0_, z_end, M, Mt = _pmap(z, T, tms_, hs_T, True)
            return _bordered(*_closing(z0_, z_end, M, tms_), Mt, z0_, Mt)

        ## the residual this solve drives to zero, kept for the caller; with
        ## `trace` every evaluation is recorded (a closing second pass
        ## appends to the first pass's list)
        self.shooting_residual = func_autonomous if self.autonomous else func
        if run.trace:
            if (not getattr(self, '_closing_second_pass', False)
                    or getattr(self, 'shooting_trace', None) is None):
                self.shooting_trace = []

            def _traced(f):
                def g(z, *a):
                    F, J = f(z, *a)
                    self.shooting_trace.append(
                        (np.array(z, dtype=float), np.array(F, dtype=float),
                         None if callable(J) else np.array(J, dtype=float)))
                    return F, J
                return g
            func, func_autonomous = _traced(func), _traced(func_autonomous)
        elif not getattr(self, '_closing_second_pass', False):
            self.shooting_trace = None

        ## THE SHOOTING RESIDUAL IS IN SOLUTION UNITS, NOT KCL UNITS.
        ## `x0 - phi(x0)` is a difference of SOLUTIONS -- volts on node rows,
        ## amps on branch rows -- so its absolute floor is the `xtol` flavour
        ## (vabstol on nodes, iabstol on branches), not the residual flavour
        ## the transient's Newton uses for `i(x)`.  Getting that backwards is
        ## F6(a)'s defect, and it is easy to walk into here because the
        ## quantity is called a residual.
        _tol = analysis.newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[1]
        (_tol,) = remove_row_col((_tol,), irefnode, self.toolkit)

        ## The shooting criterion, expressed against the transient one.
        _ratio = float(self.par.steadyratio)
        if _ratio < 1.0:
            raise ValueError(
                'steadyratio must be >= 1 (got %g): the period map is only '
                'known to the accuracy of the per-timestep solves, so a '
                'shooting tolerance tighter than reltol asks the outer '
                'residual to resolve its own noise.' % _ratio)
        _shoot_reltol = self.par.reltol * _ratio
        _tol = _tol * _ratio

        ## ⚠ REFUSED RATHER THAN SILENTLY IGNORED (as `matrix_free` is below
        ## on a kind without it).  A solved-history method already solves for
        ## `x_0` and `x_{-1}` directly and manufactures nothing, so the flag
        ## would be a no-op -- and a no-op flag that the caller believes
        ## changed something is worse than an error.
        if x0_unknown and solved_history:
            raise NotImplementedError(
                'PSS: x0_unknown=True has nothing to change for a two-step '
                "method (method=%r). The solved-history formulation already "
                'solves for x_0 and x_{-1} as real trajectory states and '
                'manufactures no opening step, so its Jacobian is exact '
                'without this. Drop the flag, or use a one-step method '
                "(method='trap' or 'euler') where the manufactured opening "
                'is what this replaces.' % method)

        ## ⚠ THE OUTER NEWTON IS DAMPED, as is standard for limit cycles
        ## (Brachtendorf et al., TCAD 33(6) 867-878: "shooting ... IN
        ## CONJUNCTION WITH A DAMPED NEWTON METHOD").  The full step is tried
        ## first and kept whenever it improves the residual, so a converging
        ## solve is unchanged; the halving only runs where the undamped
        ## iteration would have moved uphill.

        ## Find the periodic steady state: ONE Newton for every kind.  The
        ## unknown is the entering state -- gear's PAIR `(x_0, x_{-1})`,
        ## seeded `x_{-1} = x_0` (the plain formulation's assumption), so a
        ## pair run starts where a plain one does -- and on an autonomous
        ## circuit the period joins it.  Its row is the phase condition, in
        ## the units of the coordinate it pins (`_tol[phase_k]`), while the
        ## UNKNOWN it adds is a time whose own floor must be a time: mixing
        ## the two is flavour error F6(a) one row further out.
        m = n - 1
        _width = 2 if _kind == 'pair' else 1
        ## ⚠ EVERY KIND RUNS MATRIX-FREE.  The stage kinds' factored map has
        ## the mat-vec and the walk carries the period column without the
        ## dense map (radau, trbdf2 and esdirk43 matched their dense solves
        ## to 1e-15; refused until 2026-09-24 as "a dense stage product").  A
        ## Nordsieck GLM's factored map acts on its Nordsieck state; the
        ## Newton's operator on `x_0` is `x_matvec`, the recursion seeded by
        ## the linearised startup (`_GLMStartup`) -- without it, matrix-free
        ## glm2 and glm3 diverged on a driven RLC.
        z0 = np.concatenate([np.asarray(x, dtype=float)] * _width)
        tol_z = np.concatenate([_tol] * _width)
        _mf = None
        if matrix_free:
            ## ⚠ THE MONODROMY IS NOT FORMED, so it must not be REPORTED
            ## either: `_monodromy` survives from any earlier traversal and
            ## `spectral_radius` reads it without knowing which run wrote it.
            self._monodromy = None

            def _mf_build(zz):
                """RECORDED SCOPE ITEM 6: the Newton's residual and its
                Jacobian as a MAT-VEC, from the factored period (`m` columns
                on the plain map, `2m` on the pair, never formed).  With the
                period an unknown, `dphi/dT` is ONE column independent of the
                Krylov direction, computed once per Newton iteration:

                    J [v; s] = [ (I - M) v - s dphi/dT ; v_k ]

                -- one phase row, pinning the `x_0` block only, as the dense
                system does."""
                if self.autonomous:
                    z, T_ = zz[:-1], float(zz[-1])
                    tms_, hs_ = self._period_grid(T_, npts, self._grid_fracs)
                else:
                    z, T_, tms_, hs_ = zz, period, times, hs
                w_ = self._walk(_kind, z, tms_, hs_, T=T_, dense=False,
                                keep=True, want_dT=self.autonomous,
                                open_at_x0=x0_unknown)
                fp_ = w_.factored(self)
                _mv = fp_.x_matvec if _kind == 'glm' else fp_.matvec
                z0_ = w_.z0()
                Mt_ = w_.period_column() if self.autonomous else None
                F_ = self._close_periodic(z0_, w_.end(), tms_)
                if not self.autonomous:
                    return F_, (lambda v: v - alpha * _mv(v))
                Mt_ = np.asarray(Mt_, dtype=float).ravel()
                k_, r_ = _phase_row(z0_, Mt_)

                def mv_(w):
                    v_, s_ = w[:-1], float(w[-1])
                    top = (v_ - alpha * _mv(v_)) - s_ * Mt_
                    return np.concatenate((top, [v_[k_]]))
                return np.concatenate((F_, [r_])), mv_

            if run.trace:
                _mf_build = _traced(_mf_build)

            def _mf(z0_, ab_, xt_, rt_, mi_):
                return self._matrix_free_newton(_mf_build, z0_, ab_, xt_,
                                                rt_, mi_)
        ## a state-event stage may follow (below): stage 1's own stall
        ## diagnosis then waits for its outcome
        _staged = (state_events and not matrix_free
                   and (_kind in ('stage', 'pair', 'glm')
                        or (_kind == 'plain' and self._open_at_x0)))
        if self.autonomous:
            zT0 = np.concatenate((z0, [period]))
            abstol_z = np.concatenate((tol_z, [_tol[phase_k]]))
            xtol_z = np.concatenate((tol_z, [1e-15 * period]))
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous, zT0, abstol_z, xtol_z, _shoot_reltol,
                maxiterations, period, solver=_mf, defer_diagnosis=_staged)
            self.period = period = float(z_ss[-1])
            z_ss = z_ss[:-1]
            ## the grid follows the solved period; everything downstream --
            ## the replay, the waveform, the DFT -- must use it, or the
            ## answer is reported on a period the solver rejected
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif matrix_free:
            z_ss, _info, _ier, _mesg = _mf(z0, tol_z, tol_z, _shoot_reltol,
                                           maxiterations)
        else:
            z_ss, _info, _ier, _mesg = analysis.fsolve(
                func, z0, maxiter=maxiterations, reltol=_shoot_reltol,
                abstol=tol_z, xtol=tol_z, toolkit=self.toolkit,
                full_output=True, line_search=True, floor_detect=True)
        ## THE STATE EVENTS AS NEWTON UNKNOWNS: a second, bordered stage from
        ## stage 1's orbit, for every kind that has one -- the stage methods
        ## and gear's pair, driven or free period (see `_state_event_stage`).
        ## ⚠ ALSO FROM A STAGE 1 THAT DID NOT CONVERGE, unless it collapsed
        ## onto a trivial root.  Across a sharp switch the UNSTAGED map is
        ## nearly discontinuous in the state -- where the crossing falls in
        ## its step moves with every iterate -- and its Newton can fail where
        ## the staged system, the crossings landed, converges.  The stage's
        ## own Newton decides `converged`; `stage_one_converged` records
        ## whether stage 1 did.  History: `doc/shooting_history.md`,
        ## `PSS._shoot`.
        _info_one, _ier_one = _info, _ier
        if _staged and (_ier == 1 or not (isinstance(_info, dict)
                                          and _info.get('collapsed'))):
            (z_ss, _info, _ier, _mesg, period, times,
             hs) = self._state_event_stage(
                _kind, z_ss, _info, _ier, _mesg, period, times, hs,
                maxiterations, _tol, _shoot_reltol, alpha,
                *((_phase_row, phase_k) if self.autonomous else ()))
            if self.autonomous:
                self.period = period
        self.stage_one_converged = (
            (_ier_one == 1) if self._state_event_fracs is not None else None)
        _stall = (_info_one.get('stall_diagnosis')
                  if isinstance(_info_one, dict) else None)
        if _stall is not None and _ier != 1:
            _stall()
        x0_ss = z_ss[:m]
        if _kind == 'pair':
            xm1_ss = z_ss[m:]
        run.period, run.times, run.hs = period, times, hs
        run.x0_ss, run.xm1_ss = x0_ss, xm1_ss
        run.info, run.ier = _info, _ier

    def _report_convergence(self, run):
        """`solve`, phase 3: the convergence flag, the stalled-step and
        non-convergence diagnostics, the Floquet report."""
        _ier, _info = run.ier, run.info
        maxiterations, method = run.maxiterations, run.method
        self.converged = (_ier == 1)
        ## ⚠ WHY A SOLVE THAT HAD STOPPED MOVING STILL FAILED (see `fsolve`'s
        ## `floor_detect`: counted there, never acted on).  The generic
        ## non-convergence warning cannot tell a solve that is lost from one
        ## that is sitting ON its answer; this one can, and names both causes.
        self.step_floor = (_info.get('step_floor')
                           if isinstance(_info, dict) else None)
        if not self.converged and self.step_floor:
            _sf = self.step_floor
            warnings.warn(
                'PSS: the shooting solve STOPPED MOVING AND STILL FAILED ITS STEP '
                'TEST: the periodicity residual has met its tolerance since '
                'iteration %d, but the Newton step of unknown %d stays at %.1e '
                'against a tolerance of %.1e (%.0fx) and does not contract.  Two '
                'causes look like this. (1) THE ARITHMETIC FLOOR: the step is '
                'rounding -- typically a long traversal carrying an unknown of '
                'very different magnitude, tested on a node near a zero crossing '
                'where only the absolute tolerance is left -- and the waveform IS '
                'the solution; an absolute tolerance (vabstol / iabstol) at or '
                'above ~%.0e ends the solve as soon as it is there. (2) A '
                'SINGULAR I - M: a Floquet multiplier at 1 (an autonomous or '
                'marginally stable circuit solved at a fixed period), so the '
                'periodic solution is NOT UNIQUE and the step wanders along the '
                'null direction; check `spectral_radius`.'
                % (_sf['since'], _sf['index'], _sf['step'], _sf['tol'],
                   _sf['ratio'], _sf['step']),
                RuntimeWarning, stacklevel=3)
        self.shooting_iterations = maxiterations if not self.converged else None
        ## ⚠ AN AUTONOMOUS OSCILLATOR CANNOT BE SOLVED AT A FIXED PERIOD.  A
        ## self-sustaining oscillation (a VCO macromodel, an LC or ring
        ## oscillator, a DC-driven phase accumulator) has a one-parameter
        ## family of periodic solutions, so its monodromy has an eigenvalue at
        ## exactly 1 and `I - M` is singular AT the true period; away from it
        ## the discretised orbit does not close and the period map's only
        ## fixed point is the origin.  No fixed period gives both a solution
        ## and an invertible Jacobian, which is why an autonomous circuit is
        ## solved on the free-period system.
        rho, self.floquet_multipliers, self.parasitic_roots = \
            self._spectral_report(getattr(self, '_monodromy', None))
        self.spectral_radius = rho
        ## `self.autonomous` was decided before the solve and chose which
        ## system ran; the period was an unknown and `self.period` holds what
        ## it came to.
        if not self.converged:
            ## ⚠ NEVER SILENT: an unconverged solve still returns a
            ## plausible-looking waveform (the last iterate).  The advice is
            ## what is backed: the plain path is a contraction with a linear
            ## rate (the true `dF/dx_in` is singular, so no method solves a
            ## true Newton in its frame), while the solved-history route has
            ## an exact Jacobian and converges quadratically.  ⚠ THE GEAR
            ## ADVICE IS FOR A DRIVEN SOLVE ONLY: on an autonomous oscillator
            ## near a unit second multiplier gear is the method that STALLS
            ## (see `_diagnose_lmm_free_period_stall`).
            ## History: `doc/shooting_history.md`, `PSS._report_convergence`.
            if getattr(self, 'autonomous', False):
                _advice = ("On an oscillator, method='radau' (the default) is "
                           'the robust choice; see any preceding diagnosis.')
            else:
                _advice = ("Raise maxiterations, or use method='gear', whose "
                           'solved-history formulation has an exact Jacobian '
                           'and converges quadratically where the plain path '
                           'is a contraction with a linear rate.')
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r). ⚠ The returned waveform IS STILL A FULL '
                'RESULT -- it is the last iterate, not a periodic steady '
                'state -- so a reader who does not check `converged` gets '
                'an array that looks like an answer and is not. %s'
                % (maxiterations, method, _advice),
                RuntimeWarning, stacklevel=3)
        

    def _replay_orbit(self, run):
        """`solve`, phase 4: replay the CONVERGED period the way the solve
        opened it, collecting the per-step LTE.  Returns ``(X, walk,
        lte_seen)`` -- the states, the `(t, h)` pairs, the LTE readings."""
        (solved_history, x0_ss, xm1_ss, times, hs, period, x0_unknown,
         method) = (run.solved_history, run.x0_ss, run.xm1_ss, run.times,
                    run.hs, run.period, run.x0_unknown, run.method)
        ## THE THIRD LEVEL, MEASURED ON THE WAY OUT.  The inner Newton and the
        ## shooting Newton both ask whether an EQUATION was solved; neither
        ## asks whether it was the right one -- the discrete period map is
        ## not the continuous one, and the gap is truncation error.  PSS
        ## imposes its grid, so this cannot be a CONTROL signal (shrinking a
        ## step would change the period map between shooting iterations and
        ## destroy the monodromy); it is a MEASUREMENT on the converged
        ## solution over the final replay, `|J^-1 Eg| / (TRTOL (reltol ref +
        ## lte_abstol))`, the quantity a transient would reject a step on.
        ## The nesting works because shooting is a multilevel Newton: the
        ## inner transient absorbs the nonlinearity, so the outer Newton sees
        ## a nearly linear `phi_T` (Kundert, *Introduction to RF
        ## Simulation*, 2003; cited, not verified here).
        ##
        ## ⚠ THE REPLAY MUST OPEN THE WAY THE SOLVE DID, or the waveform is
        ## not the solution: a plain replay of a solved-history answer would
        ## reintroduce the seam that formulation removes.
        ## ⚠ AND IT MUST WALK THE SAME (t, h) PAIRS.  The plain walk
        ## takes the MANUFACTURING step at `(times[0], hs[0])` FIRST and
        ## only then walks `times[1:]` with `hs[_j]`, so the step after the
        ## opening one uses `hs[0]` again; indexing `times` and `hs` in
        ## parallel is off by one against it, which a uniform grid hides
        ## completely.  Built as explicit `(t, h)` PAIRS, which cannot be
        ## misaligned by one.
        ##
        ## `_period_state` keeps what a later factored replay cannot
        ## re-derive (which seed, which grid, which opening), and
        ## `factored_period()` runs the traversal on demand: the matrix-free
        ## Newton's last factorisations belong to its last TRIAL iterate, and
        ## retaining `N` of them from every solve would cost `2 N m^2`
        ## doubles on every run.
        ## History: `doc/shooting_history.md`, `PSS._replay_orbit`.
        self._period_state = (bool(solved_history), copy(x0_ss),
                              None if xm1_ss is None else copy(xm1_ss),
                              times, hs, float(period), bool(x0_unknown))
        self._factored_period_cache = None

        if solved_history:
            self._install_history(x0_ss, xm1_ss, hs[0], h_prev=hs[-1])
            tr = self._transient()
            X = [np.asarray(x0_ss, dtype=float)]
            walk = list(zip(times[1:], hs[:len(times) - 1]))
        else:
            X = [x0_ss]
            tr = self._begin_period(x0_ss)
            ## the manufacturing step, then the loop -- exactly the plain walk
            ## ... unless there was no manufacturing step, in which case the
            ## replay opens AT `x_0` and walks the period alone.  Getting
            ## this wrong is the same class of defect as the grid shift
            ## above: a replay that does not reproduce its own traversal.
            walk = list(zip(times[1:], hs[:len(times) - 1]))
            if not x0_unknown:
                walk = [(times[0], hs[0])] + walk
        ## Fresh probe, so `relref='sigglobal'`'s running signal maximum is
        ## the period's, not something an earlier shooting iteration saw.
        tr._lte_probe = None
        ## A stage method has no LMM divided-difference LTE (compute_lte
        ## refuses), and the seam/interior split is a property of a manufactured
        ## opener it does not have -- so the replay collects no per-step LTE for
        ## it, and the three LTE figures below report None (honestly: the
        ## diagnostic does not apply to a self-starting stage method).  The
        ## method says which it is.
        self._want_lte = not self._integrator_for(method).is_stage_method()
        lte_seen = []
        for t, dt in walk:
            x = self.solve_timestep(X[-1], t, dt)
            if self._lte is not None:
                lte_seen.append((float(self._lte), float(t), self._lte_seam,
                                 self._lte_valid))
            X.append(copy(x))
        self._want_lte = False

        return X, walk, lte_seen

    def _report_lte(self, run, lte_seen):
        """`solve`, phase 5: the three truncation-error figures and their
        warning."""
        method, npts = run.method, run.npts
        ## THREE NUMBERS, BECAUSE THEY HAVE DIFFERENT REMEDIES.
        ##
        ## `max_lte` is the INTERIOR per-step peak -- steps whose estimator
        ## saw only real past charges -- exactly the quantity a transient
        ## controls its grid on.
        ##
        ## `total_lte`, the SUM over the interior, is the one that sees what
        ## a WHOLE PERIOD does: a per-step peak can be in tolerance while the
        ## orbit is badly damped (euler on a Q=20 resonator: peak 0.288, sum
        ## 26.27, amplitude 56% low).  It is an upper bound -- it adds
        ## magnitudes, so it cannot see cancellation -- which is the right
        ## direction for a diagnostic to be wrong in.
        ##
        ## `max_lte_seam` is the peak over the opening steps of a method
        ## whose COMPANION reads the entering unknown -- Gear-2 on the plain
        ## formulation, and `None` for euler and trapezoidal, which cannot
        ## have a seam, and on the solved-history path, which removes it.
        ## See `solve_timestep`.  ⚠ IT IS A FLAG, NOT A MAGNITUDE: the
        ## estimator differences a fabricated charge while the solution merely
        ## reads one, so the number overstates the seam by orders of
        ## magnitude; `benchmarks/pss_seam_cost.py` measures what it costs.
        ## It does not improve with a smaller timestep.
        ##
        ## An unsound estimate is reported as neither: for trapezoidal that
        ## step is the ONLY one whose number was ever wrong, and dropping it
        ## is what keeps it out of both figures.
        ## History: `doc/shooting_history.md`, `PSS._report_lte`.
        interior = [p for p in lte_seen if not p[2] and p[3]]
        seam = [p for p in lte_seen if p[2]]
        if interior:
            self.max_lte, self.max_lte_time = max(interior)[:2]
            self.total_lte = float(sum(p[0] for p in interior))
        else:                                            # pragma: no cover
            self.max_lte = self.max_lte_time = self.total_lte = None
        self.max_lte_seam = max(seam)[0] if seam else None

        ## Named so the warning can lead with whichever is actually
        ## limiting: the three have three different answers.
        _limits = [
            (self.total_lte, 'accumulated over the period',
             'use a smaller timestep or a less damping method -- this is '
             'the figure that sets a limit cycle, and a per-step criterion '
             'can be in tolerance while it is not'),
            (self.max_lte, 'in one interior step',
             'use a smaller timestep or a higher-order method'),
            (self.max_lte_seam, 'over the opening steps',
             "this is the period map's own seam, where each shooting "
             'iteration cold-starts from a fabricated history; it does '
             'NOT improve with a smaller timestep'),
        ]
        _over = [(v, where, fix) for v, where, fix in _limits
                 if v is not None and v > 1.0]
        if _over:
            v, where, fix = max(_over, key=lambda r: r[0])
            ## ⚠ DO NOT ASSERT CONVERGENCE HERE: the LTE report is produced
            ## whether or not the solve converged, and must not contradict
            ## the non-convergence warning beside it.
            warnings.warn(
                'PSS: the shooting solve %s, and the periodic '
                'solution is not resolved at this accuracy (method=%r, %d '
                'points per period). Local truncation error reaches %.3g '
                'times tolerance %s: %s. Neither Newton criterion can see '
                'this -- they ask whether the discrete equations were '
                'solved, not whether the discretisation is the right one. '
                '(peak interior %s at t=%.6g s, period total %s, opening '
                'steps %s; relax lte_vabstol/lte_iabstol/TRTOL if this '
                'accuracy is intended.)'
                % ('converged' if self.converged else 'did NOT converge',
                   method, npts, v, where, fix,
                   'n/a' if self.max_lte is None else '%.3g' % self.max_lte,
                   -1.0 if self.max_lte_time is None else self.max_lte_time,
                   'n/a' if self.total_lte is None
                   else '%.3g' % self.total_lte,
                   'n/a' if self.max_lte_seam is None
                   else '%.3g' % self.max_lte_seam),
                RuntimeWarning, stacklevel=3)

    def _check_fundamental(self, X, walk, period):
        """`solve`, phase 6: warn when an autonomous solve returned a
        MULTIPLE of the fundamental (sets `fundamental_period`)."""
        ## ⚠ AN AUTONOMOUS PERIOD IS ONLY DETERMINED UP TO AN INTEGER
        ## MULTIPLE, AND THE SOLVE FOLLOWS THE SEED.  `k*T` is a period
        ## whenever `T` is, so the free-period system converges to whichever
        ## multiple the seed is nearest and reports `converged` with a
        ## correct periodic waveform whose FUNDAMENTAL is wrong by the factor.
        ##
        ## The detector is cheap and needs no extra solve: an orbit
        ## traversed k times comes back near `x_0` partway through.  Grid
        ## points do not land on `T/k` in general, so this is a
        ## NEAREST-APPROACH test against the orbit's own diameter rather than
        ## an equality, and the endpoints are excluded because every orbit is
        ## near `x_0` there.  Driven runs are exempt: their period is the
        ## caller's, and asking for two source periods is legitimate.
        ## History: `doc/shooting_history.md`, `PSS._check_fundamental`.
        if self.autonomous and len(X) > 8:
            ## ⚠ THE ORBIT IS ITS DIFFERENTIAL STATE, WITH THE PERIODIC ROWS
            ## FOLDED.  Two points of an autonomous orbit coincide exactly
            ## when their differential states do -- the algebraic rows are
            ## functions of them -- and a phase coincides modulo its modulus.
            ## Over the whole vector an idtmod VCO could never recur (its
            ## unfolded phase state advances one modulus per fundamental), and
            ## its wrapped output jumps by a whole modulus mid-period.
            _pts = np.array([np.asarray(v, dtype=float).ravel() for v in X])
            _dyn = np.any(np.asarray(self._C_at(_pts[0]), dtype=float) != 0.0,
                          axis=0)
            if not np.any(_dyn):
                _dyn[:] = True
            _dev = _pts - _pts[0]
            for _r, _m, _o in self._periodic_fold:
                _dev[:, _r] -= _m * np.round(_dev[:, _r] / _m)
            _d = np.max(np.abs(_dev[:, _dyn]), axis=1)
            _diam = float(np.max(_d))
            ## ⚠ THE THRESHOLD IS THE ORBIT'S OWN SPEED AT `x_0`.  A k-fold
            ## orbit passes `x_0` again at the SAME speed, between grid
            ## points, so its nearest point is within about half a step's
            ## displacement there -- and an orbit moving away from `x_0` is
            ## two steps' displacement off by the excluded edge.  Not the
            ## LARGEST step on the orbit: that is an output wrap's jump or a
            ## fast edge.  `_h` is each step of the replay, so non-uniform
            ## grids scale locally.
            _h = np.array([h for _t, h in walk], dtype=float)
            _v0 = _d[1] / _h[0]
            ## ⚠ THE EXCLUDED EDGES ARE TIME, NOT POINTS: 5 % of the period
            ## (at least two steps' worth), which on a uniform grid is the
            ## same `max(2, N // 20)` points.  Counted in points, a grid that
            ## OPENS WITH A RAMP (gear's, after a landed state event: ten
            ## doubling steps from 1e-5 T) excluded 3e-4 of the period, and
            ## the orbit still leaving `x_0` -- displacement about the speed
            ## times the step, since the elapsed time IS about the step --
            ## read as a recurrence ("traverses it about 3182 times").
            _tj = np.concatenate(([0.0], np.cumsum(_h)))
            _span = float(_tj[-1])
            _edge_t = max(2, len(_d) // 20) * _span / len(_h) * (1.0 - 1e-9)
            ## ⚠ THE EARLIEST RECURRENCE, NOT THE NEAREST.  A three-fold
            ## orbit passes close to `x_0` at both `T/3` and `2T/3`, and the
            ## later one is itself a multiple.
            _near = [j for j in range(1, len(_h))
                     if _edge_t <= _tj[j] <= _span - _edge_t
                     and _d[j] < _v0 * max(_h[j - 1], _h[j])
                     and _d[j] < 0.25 * _diam]
            if _diam > 0.0 and _near:
                ## The closest approach WITHIN THE FIRST cluster: the first
                ## point over the threshold can be up to a step early.
                _run = [_near[0]]
                for _c in _near[1:]:
                    if _c != _run[-1] + 1:
                        break
                    _run.append(_c)
                _j = min(_run, key=lambda i: _d[i])
                ## the time from the reference point, summed over the
                ## replay's own steps: with a caller's grid the points are
                ## not evenly spaced, and on the plain path `X[0]` sits one
                ## step before `t = 0`
                self.fundamental_period = float(np.sum(_h[:_j]))
                warnings.warn(
                    'PSS: this autonomous solve returned a period that '
                    'is a MULTIPLE of the fundamental. The orbit comes '
                    'back within %.2g of its own diameter at t=%.6g s, '
                    'so the fundamental is about %.6g s and the '
                    'returned %.6g s traverses it about %.1f times. '
                    'k*T solves the periodicity condition whenever T '
                    'does, so the solve follows its seed -- re-run with '
                    'period=%.6g to get the fundamental. The waveform '
                    'is a correct periodic solution either way; its '
                    'FUNDAMENTAL FREQUENCY is what is off.'
                    % (_d[_j] / _diam, self.fundamental_period,
                       self.fundamental_period, period,
                       period / self.fundamental_period,
                       self.fundamental_period),
                    RuntimeWarning, stacklevel=3)

    def _orbit_results(self, run, X):
        """`solve`, phase 7: the reported waveform (`self.waveform`) and the
        two results, `tpss` in time and `fpss` in frequency."""
        toolkit = self.toolkit
        (times, irefnode, solved_history, x0_unknown) = (
            run.times, run.irefnode, run.solved_history, run.x0_unknown)
        ## ⚠ THE FIRST ENTRY IS DROPPED ONLY WHEN IT IS NOT PART OF THE
        ## PERIOD.  On the default plain path `X[0]` is `x_in`, the
        ## pre-image of the manufactured step, which sits one step BEFORE
        ## t=0 and is not a point of the orbit.  With `x0_unknown` -- and on
        ## the solved-history path -- `X[0]` IS `x(0)`, so dropping it would
        ## discard a real sample, shift the waveform by a step and leave it
        ## one column short of `times`.
        X = toolkit.array(X if (solved_history or x0_unknown) else X[1:]).T

        # Insert reference node voltage
        X = toolkit.concatenate((X[:irefnode], 
                                 toolkit.zeros((1,len(times))), 
                                 X[irefnode:]))

        tpss = analysis.CircuitResult(self.cir, x=X, xdot=None,
                                      sweep_values=times, sweep_label='time', 
                                      sweep_unit='s')

        ## ⚠ KEPT FOR THE CARRIER PHASOR, which AM/PM needs and which
        ## `fpss` below cannot supply: `freq_analysis` returns an RMS,
        ## energy-folded, positive-frequency spectrum -- right for
        ## reporting and wrong for a phasor, because folding destroys the
        ## phase relationship between a carrier and its sidebands, which is
        ## the entire content of an AM/PM decomposition.
        self.waveform = (np.asarray(times, dtype=float),
                         np.asarray(X, dtype=float))

        freqs, FX = freq_analysis(X[:,:-1], times[:-1])
        ## ⚠ ON A NON-UNIFORM GRID (`grid=`) THE INDEX DFT ABOVE IS NOT A
        ## FOURIER COEFFICIENT and does not converge.  Same layout and RMS
        ## fold, taken as a weighted sum at the true times (uniform grids
        ## never reach this branch).
        ## History: `doc/shooting_history.md`, `PSS._orbit_results`.
        _h = np.diff(np.asarray(times, dtype=float))
        if len(_h) >= 2 and float(np.max(_h)) / float(np.min(_h)) - 1.0 \
                > self.UNIFORM_GRID_TOL:
            _tt = np.asarray(times, dtype=float)
            _Tp = float(_tt[-1] - _tt[0])
            ## the SAME rule `_period_quadrature` gives every consumer: a
            ## periodic cubic spline on an event-free grid, a piecewise one
            ## breaking at the landed event nodes -- `fpss` and
            ## `carrier_phasor` are pinned equal to 1e-12
            _wq = periodic_spline_weights(_tt[:-1], _Tp, self._event_nodes(_tt[:-1], _Tp)) / _Tp
            _ks = np.arange(len(freqs))
            freqs = _ks / _Tp
            _E = np.exp(-2j * np.pi * np.outer(_ks, (_tt[:-1] - _tt[0]) / _Tp)) \
                * _wq[None, :]
            FX = np.asarray(X[:, :-1], dtype=float) @ _E.T
            FX[:, 1:] *= np.sqrt(2)

        ## ⚠ `fpss` IS RMS, AND THE USUAL THING TO COMPARE IT AGAINST IS NOT.
        ## `freq_analysis` returns an RMS, energy-folded, positive-frequency
        ## spectrum; a commercial simulator's frequency-domain PSS output is
        ## conventionally PEAK, so a harmonic compared across the two differs
        ## by `sqrt(2)` -- 3.01 dB -- with nothing announcing it.  Multiply
        ## `fpss` by `sqrt(2)` for a peak-convention comparison, or divide
        ## theirs.
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return tpss, fpss

    def _closing_polish(self, run):
        """`solve`, phase 8: after a 'closing' free-period solve, solve once
        more proportionally on the caller's fractions from the converged
        state, and return THAT result; None when no second pass is due."""
        (refnode, period, hs, x0_ss, maxiterations, matrix_free, x0_unknown,
         phase_rule) = (run.refnode, run.period, run.hs, run.x0_ss,
                        run.maxiterations, run.matrix_free, run.x0_unknown,
                        run.phase_rule)
        ## ⚠ THE SECOND PASS, ALWAYS.  'closing' keeps a caller's inner steps
        ## where the transient validated them, which widens the free-period
        ## Newton's basin (it converges from a seed period 16 % off where
        ## 'proportional' fails its per-step Newton).  But the closing step
        ## then absorbs the whole period correction, and ANY stretch of the
        ## last step is a change of discretisation that stays in the answer
        ## -- even from an exact seed, where the last step absorbs the
        ## discretisation's own period error.  So once converged the grid is
        ## re-fractioned at the solved period and solved once more
        ## proportionally from the converged state: the sane grid the answer
        ## is reported on.  Cost: one Newton from a converged state, one or
        ## two iterations.
        ## History: `doc/shooting_history.md`, `PSS._closing_polish`.
        if (getattr(self, '_period_column', 'proportional') == 'closing'
                and self.converged and getattr(self, 'autonomous', False)
                and not getattr(self, '_closing_second_pass', False)
                and len(hs) > 2):
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            _hs = np.asarray(hs, dtype=float)
            _r = float(_hs[-1] / _hs[-2])
            if _r > ZERO_STABILITY_RATIO or _r < 1.0 / ZERO_STABILITY_RATIO:
                warnings.warn(
                    'PSS: the closing step ended %.2fx its neighbour after the '
                    'free-period solve moved the period from %.6g to %.6g s; '
                    'solving once more on that grid re-fractioned at the '
                    'solved period (proportional), from the converged state.'
                    % (_r, float(self._solve_kwargs.get('period_seed', period)),
                       float(period)), RuntimeWarning, stacklevel=3)
            ## ⚠ ON THE CALLER'S FRACTIONS, not the closing-distorted grid:
            ## re-fractioning THAT grid keeps the giant last step.
            _fr_caller = (np.asarray(self._grid_fracs, dtype=float)
                          if self._grid_fracs is not None else None)
            self._closing_second_pass = True
            self._force_period_column = 'proportional'
            try:
                return self.solve(refnode=refnode, period=float(period),
                                  x0=copy(x0_ss),
                                  timestep=float(period) / len(_hs),
                                  maxiterations=maxiterations,
                                  grid=_fr_caller,
                                  matrix_free=matrix_free,
                                  x0_unknown=x0_unknown, tstab=None,
                                  break_events=self.break_events,
                                  phase_rule=phase_rule, trace=run.trace)
            finally:
                self._closing_second_pass = False
                self._force_period_column = None
