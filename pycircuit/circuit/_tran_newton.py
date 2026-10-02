"""The per-step Newton, its continuation rescue and both tolerance flavours
(the Newton residual and the LTE floor).  A theme of `Transient` (see
`transient.py`).
"""

import numpy as np
from numpy.linalg import LinAlgError

from pycircuit.circuit._limiting import stateful_limiters
from pycircuit.circuit.analysis import (
    NoConvergenceError,
    SingularMatrix,
    reduced_row_names,
    remove_row_col,
)
from pycircuit.circuit.dcanalysis import refnode_removed

#: A compiled function running as C (`hdl.set_backend('c')`) counts this
#: fraction of its bytecode in `compiled_jacobian_size`: measured
#: 2026-10-01, the C kernel evaluates `G` 53x (MosLevel1, 206 -> 3.9 us) to
#: 290x (PSP, 14 ms -> 49 us) faster than the numpy path.  2026-10-02,
#: against the CSE twin that runs and the C kernel with its libm calls
#: declared const: 11x (SPICE diode, 37 -> 3.3 us) to 142x (PSP, 3.2 ms ->
#: 22 us), and the placement it sets measured again (`AUTO_JACOBIAN_CODE`).
C_KERNEL_SHARE = 1.0 / 100.0


def compiled_jacobian_size(cir):
    """The bytecode, in bytes, of the compiled `G` and `C` of every hdl
    element of `cir` at any depth (`Behavioural` classes: `_hdl_info`) --
    a deterministic measure of what one Jacobian evaluation costs, which
    orders the models as their time does: `DiodeHdl` 0.2 KB, a HEMT 2.2 KB,
    the SPICE diode 20 KB (37 us a `G`), `MosLevel1Hdl` 26 KB (68 us),
    `MosLevel3Hdl` 92 KB (0.37 ms), the PSP MOSFET 1.8 MB (3.2 ms).  A
    function bound to its C kernel counts `C_KERNEL_SHARE` of that (PSP 18
    KB, MosLevel3 0.9 KB).  Hand-written elements and `BSource` count 0.

    The size is the REFERENCE function's -- the chain as printed, before
    the bit-identical CSE (`_hdl_cse`, `fn._hdl_codelen`) -- not that of
    the twin that runs (PSP's is a third of it).  Re-calibrated 2026-10-02
    with the twins and the fused passes running (`AUTO_JACOBIAN_CODE`), the
    reference size still parts the models where the options gain from
    those where they do not, so it was kept; and it does not move with
    `PYCIRCUIT_HDL_CSE`, so turning the CSE off stays bit-identical end to
    end -- a size read off the twin could flip a model near the threshold,
    and the run's Newton with it.  Read by the 'auto' Newton options
    (`Transient._newton_option`)."""
    total = 0
    stack = [cir]
    while stack:
        c = stack.pop()
        for e in (getattr(c, 'elements', None) or {}).values():
            if getattr(e, 'elements', None):
                stack.append(e)
            info = getattr(type(e), '_hdl_info', None)
            funcs = info.get('funcs') if isinstance(info, dict) else None
            if not funcs:
                continue
            for k in ('G', 'C'):
                fn = funcs.get(k)
                code = getattr(fn, '__code__', None)
                if code is None:
                    continue
                ## (the REFERENCE function's size when the bit-identical CSE
                ## replaced it: `_hdl_cse`; why, above)
                size = getattr(fn, '__dict__', {}).get('_hdl_codelen') \
                    or len(code.co_code)
                if getattr(fn, '__dict__', {}).get('_hdl_c') is not None:
                    size = int(size * C_KERNEL_SHARE)
                total += size
    return total


class _StepNewton:
    """The per-step Newton, its continuation rescue and both tolerance flavours
    (the Newton residual and the LTE floor).  A theme of `Transient`
    (see `transient.py`)."""

    def _rescue_solver(self, base):
        """P18 phase 3 + P25: the failed-time-point continuation chain --
        junction-gmin -> gshunt -> pseudo-transient, for exactly one
        point.  No source stepping here: scaling u(t) mid-transient would
        scale the integrator's companion history too -- ill-posed.
        Psi-tc IS mid-transient-safe: its anchor is the last accepted
        state (physical continuity) and it scales nothing.  Every ladder
        ends with a PURE solve, so an accepted rescued point carries no
        residue (the P22 rule); Psi-tc's rungs are solved by the PLAIN
        base, never the chain (the same reasoning as the DC wiring).
        Junction rows are reduced-system indices, as in DC.  Extracted
        so the wiring is testable -- the P18 scope finding stands: no
        legitimate triggering circuit could be fabricated, so the chain's
        behavior is gated at the nrsolver level and its topology here."""
        from pycircuit.circuit.nrsolver import (
            GminSteppingNewton,
            JunctionGminSteppingNewton,
            PseudoTransientNewton,
        )
        from pycircuit.circuit.pcnr import pcnr_junctions
        _jrows = []
        for _i, _e, _ra, _rb in pcnr_junctions(self.cir):
            if self.irefnode in (_ra, _rb):
                continue
            _jrows.append((_ra - (_ra > self.irefnode),
                           _rb - (_rb > self.irefnode)))
        chain = GminSteppingNewton(
            JunctionGminSteppingNewton(base, _jrows))
        return PseudoTransientNewton(chain, rung_solver=base)

    ## import it from there instead.
    ## But it's an object method requiring a DC as self
    ## so using DC._newton doesn't work
    def _honours_continuation_rescue(self):
        """Whether THIS step path can actually reach the continuation ladder.

        The stepping loop (`_rescue_step`) arms `_continuation_rescue` once the
        step has shrunk to `minstep`, as the last resort before giving up.  Every
        path now reaches a ladder, by one of three routes:

        * the LMM companions and the DIRK/ESDIRK stages solve through
          :meth:`_newton`, which reads the flag and wraps the rescue chain;
        * the FULL coupled stage solve carries its OWN gshunt ladder, because
          `_newton` is MNA-SIZED (it reduces an `n`-vector at `irefnode`, limits a
          full `n`-vector, and carries per-MNA-row tolerances) and cannot take a
          `3m` block system;
        * the PCNR variants of both stage paths have no ladder of their own -- a
          gshunt one and a junction-gmin one were built and MEASURED not to
          rescue them, because PCNR's bottleneck is the junction limiter's slew
          and no deformation of the circuit accelerates that -- so they FALL BACK
          to the device-limiting solve that does carry one, exactly as the LMM
          step and DC have always done on a PCNR failure.

        So this is `True` for everything.  `_rescue_step` consults it, so the
        failure message never says a continuation "could not rescue the point"
        on a path where none was attempted.  A new step path that reaches
        neither a ladder nor a fallback should return `False` here rather than
        inherit a message that claims a rescue it never tried.

        History: `doc/transient_history.md`, `Transient._honours_continuation_rescue`.
        """
        return True

    ## WHERE THE 'auto' NEWTON OPTIONS TURN ON (`chord_jacobian`,
    ## `radau_transform`): a circuit whose compiled Jacobian is at least this
    ## many bytes (`compiled_jacobian_size`).  A driven stage's PSS at 40
    ## points, gear chord / radau transform against the full Newton, measured
    ## 2026-10-02 with the bit-identical CSE twins and the fused passes
    ## running -- they cut the full Newton's `G` and `C` 3-5x, so the options
    ## gain less than at the 2026-10-01 calibration (17-68 %): the PSP MOSFET
    ## (1.8 MB) -25 / -68 %, switching -32 / -25 %; MosLevel3 (92 KB) -14 /
    ## -27 % (switching -19 / -16 %); Gummel-Poon (57 KB) -10 / 0 %; EKV and
    ## MosLevel1 (26 KB) -5 / -12 % and -6 / -15 % (-9 / +3 %); the SPICE
    ## diode (20 KB) -6 / -10 % -- the answers 0 to 4e-8 apart at reltol
    ## 1e-8, no fallback.  Below: the HEMT (2.2 KB) 0 / -3 %, `DiodeHdl` (0.2
    ## KB) -1 / +2 %; and where the options LOST (2026-10-01), every circuit
    ## was of hand-written elements (0 bytes): van der Pol (`BSource`) -6 /
    ## +18 %, a switching PWM loop -- built-in switches -- +15 / +79 %.  No
    ## crossover moved, so the threshold did not.
    ## ⚠ ON THE hdl C BACKEND the Jacobian is cheap, and the two part ways
    ## (2026-10-02, the libm calls declared const): the chord never lost
    ## (MosLevel1 / MosLevel3 / EKV / Gummel-Poon / SPICE diode -2 to -4 %,
    ## PSP -12 / -18 %), the transform did on the mid-sized models
    ## (switching MosLevel1 +16 %, MosLevel3 +9 %, Gummel-Poon +15 %) while
    ## PSP kept -63 / -9 %.  A C-bound function counts `C_KERNEL_SHARE` of
    ## its bytecode, so PSP stays on (18 KB) and the mid-sized models go off
    ## (under 1 KB).
    ## History: `doc/transient_history.md`, `ChordNewton`.
    AUTO_JACOBIAN_CODE = 10000

    def _newton_option(self, value, name):
        """The Newton option `name` (`chord_jacobian`, `radau_transform`),
        given as `value`, as a bool: True / False as given, 'auto' where the
        circuit's compiled Jacobian is expensive (`AUTO_JACOBIAN_CODE`;
        decided once per circuit, again after `_memo_clear`)."""
        if isinstance(value, str):
            if value != 'auto':
                raise ValueError(
                    f"{name} must be True, False or 'auto', not {value!r}")
            on = getattr(self, '_jacobian_expensive', None)
            if on is None:
                on = self._jacobian_expensive = (
                    compiled_jacobian_size(self.cir) >= self.AUTO_JACOBIAN_CODE)
            return on
        if value in (True, False):
            return bool(value)
        raise ValueError(
            f"{name} must be True, False or 'auto', not {value!r}")

    def _newton_abstol_vector_reduced(self):
        a = self._newton_abstol_vector()
        (a,) = remove_row_col((a,), self.irefnode, self.toolkit)
        return a

    def _newton_xtol_vector_reduced(self):
        t = self._newton_xtol_vector()
        (t,) = remove_row_col((t,), self.irefnode, self.toolkit)
        return t

    def _newton_limiter(self):
        """The device limiting the step's Newton applies, in reduced
        coordinates: `cir.limit` on the full vector.  One definition for the
        step's own Newton and the branch check's speculative one, which must
        solve the same equation.  `None` when nothing in the circuit
        limits: `cir.limit` would hand back its argument, and the pass
        (two inserts and a concatenate per iteration) was 6-10 % of a
        solve (2026-09-30)."""
        _has = getattr(self, '_has_limiters', None)
        if _has is None:
            from pycircuit.circuit._limiting import has_limiters
            _has = self._has_limiters = has_limiters(self.cir)
        if not _has:
            return None

        def limiter_func(xr, x0r):
            x = self.toolkit.insert(xr, self.irefnode, 0.0)
            x0_full = self.toolkit.insert(x0r, self.irefnode, 0.0)

            x = self.cir.limit(x, x0_full, self.epar)
            return self.toolkit.concatenate((x[:self.irefnode], x[self.irefnode+1:]))
        return limiter_func

    def _newton(self, func, x0, residual=None):
        ## `residual(x)`: the step residual alone, for the chord iterations
        ## (`chord_jacobian`; the multistep step passes it, see
        ## `ChordNewton`).  The continuation rescue and a caller's own
        ## strategy keep the full Newton.
        ## ⚠ THE FIRST EVALUATION IS AT THE SEED, not at the tangent of the
        ## point before it.  A stateful limiter (`Diode`) reads `i` / `G` as
        ## the tangent at its stored `_vlim`, which sits at the last solved
        ## point (`x_n`, or the previous stage), so its Newton's first
        ## iteration evaluated the predictor as that tangent, spent an
        ## iteration, and stopped at the tolerance where a state-free device
        ## converges to rounding.  Measured against the state-free twin:
        ## 22-88 % more Newton iterations (gear2 681 / 557, TR-BDF2 on a
        ## rectifier 6517 / 3462), and TR-BDF2's embedded estimate turned
        ## the tolerance-level stage error, through the diode's 700 S, into
        ## 301 rejections against 36 (1023 steps against 863).  So the
        ## limiters are moved toward the seed first -- a LIMITED move from
        ## the stored state, which a seed past the knee cannot defeat.
        lims = getattr(self, '_stateful_lims', None)
        if lims is None:
            lims = self._stateful_lims = stateful_limiters(self.cir)
        if lims:
            x0s = np.array(x0, dtype=float)
            self.cir.limit(x0s, x0s, self.epar)
        abstol = self._newton_abstol_vector()
        xtol = self._newton_xtol_vector()
        
        (x0, abstol, xtol) = remove_row_col((x0, abstol, xtol), self.irefnode, self.toolkit)
        
        limiter_func = self._newton_limiter()

        solver = self._get_nrsolver()
        chord = None
        if getattr(self, '_continuation_rescue', False):
            solver = self._rescue_solver(solver)
        elif (residual is not None and self.par.nrsolver is None
              and self._newton_option(self.par.chord_jacobian,
                                      'chord_jacobian')):
            from pycircuit.circuit.nrsolver import ChordNewton
            iref, tk = self.irefnode, self.toolkit

            def residual_reduced(xr):
                f = residual(tk.concatenate((xr[:iref], tk.array([0.0]),
                                             xr[iref:])))
                (f,) = remove_row_col((f,), iref, tk)
                return f
            solver = chord = ChordNewton(residual_reduced, solver)
        scaler = self._get_scaler()
        linsolver = self._get_linearsolver()
        try:
            x_res, _iters = solver.solve_system(
                x0,
                refnode_removed(func, self.irefnode, self.toolkit),
                self.toolkit,
                self.par.reltol,
                abstol,
                xtol,
                self.par.maxiter,
                limiter=limiter_func,
                scaler=scaler,
                linsolver=linsolver,
                ## Stage 6: lets the solver name a node instead of a row index.
                row_names=reduced_row_names(self.cir, self.irefnode),
            )
        ## NARROW, deliberately.  A broad `except Exception` would report every
        ## failure inside a device model (a `ZeroDivisionError`, a `TypeError`, an
        ## `AttributeError`) as a convergence failure, which sends the reader to
        ## look at the bias point of a circuit whose real problem is a bug three
        ## frames down.  The solvers already classify what they mean, so only their
        ## own exceptions and genuine linear-algebra failures are translated here;
        ## everything else propagates with its original type and traceback.
        ## History: `doc/transient_history.md`, `Transient._newton`.
        except SingularMatrix:
            raise
        except NoConvergenceError as e:
            ## The solvers wrap a singular factorisation as NoConvergenceError
            ## ("Singular Jacobian: ..."); promote that to SingularMatrix so callers
            ## can tell "no solution here" from "could not get there".  Matching on
            ## the message is weak.
            ## History: `doc/transient_history.md`, `Transient._newton`.
            if 'Singular' in str(e) or 'linalgerror' in str(e).lower():
                raise SingularMatrix(str(e)) from e
            ## ⚠ THE LINE SEARCH, AS A RETRY (owner decision 2026-09-08, "Do 2";
            ## the coupled Radau path carries the same retry in
            ## `_rk_step_coupled`).  The default `StandardNewton` has no
            ## damping, and a nonlinearity fed by a branch current through a
            ## capacitor cannot be integrated undamped at ANY grid (stage
            ## sensitivity tau/h, undamped basin 0.94 h/(k tau), READING-LOG
            ## 2.156).  `DampedNewton` is tried once, from the same seed, only
            ## after the plain Newton has failed and only when the caller left
            ## the strategy at its default -- so every step that converged
            ## before converges to the same numbers, and a caller-chosen
            ## strategy keeps its own failure.
            from pycircuit.circuit.nrsolver import DampedNewton
            ## Only once the rescue ladder is ARMED (it wraps `solver`, so the
            ## plain Newton and the ladder have both had their turn) and only
            ## for the default strategy -- the search is the last resort, for
            ## the same reason as on the coupled path above.
            if (getattr(self, '_continuation_rescue', False)
                    or getattr(self, '_damped_last_resort', False)) \
                    and self.par.nrsolver is None:
                try:
                    x_res, _iters = DampedNewton().solve_system(
                        x0,
                        refnode_removed(func, self.irefnode, self.toolkit),
                        self.toolkit,
                        self.par.reltol,
                        abstol,
                        xtol,
                        self.par.maxiter,
                        limiter=limiter_func,
                        scaler=scaler,
                        linsolver=linsolver,
                        row_names=reduced_row_names(self.cir, self.irefnode),
                    )
                except NoConvergenceError:
                    raise e
            else:
                raise
        except LinAlgError as e:
            raise SingularMatrix(str(e)) from e
        
        ## Stage 6(c): the Newton iterations, counted.
        ## History: `doc/transient_history.md`, `Transient._newton`.
        stats = getattr(self, 'statistics', None)
        if stats is not None:
            stats.newton_iterations += int(_iters)
            if chord is not None:
                stats.chord_fallbacks += chord.fallbacks

        ## BRANCH DETECTION -- the screen is O(m) unless something collapsed
        if self._branch_on():
            ## ⚠ the REDUCED function, the one the solver actually solved --
            ## handing it the full-width `func` with a reduced seed is a
            ## dimension mismatch that the diagnostic's own except would eat
            self._branch_after_solve(
                refnode_removed(func, self.irefnode, self.toolkit), x_res)

        x = x_res
        
        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

    ## STAGE 12B -- small helpers the coupled path needs, factored out of `_solve`
    ## rather than re-derived, so the two paths cannot drift apart on tolerances.

    ## The LTE tolerance multiplier.  `TRTOL` in this module, `lteratio` in a
    ## commercial simulator: the LTE estimate is deliberately conservative, so
    ## the allowed truncation error is this many times the Newton-solve
    ## tolerance.
    ## A property reading the `TRTOL` Parameter (the JAX backend's too, P2),
    ## so every `self.LTERATIO` read follows the Parameter and the two cannot
    ## drift.
    ## History: `doc/transient_history.md`, `Transient.LTERATIO`.
    @property
    def LTERATIO(self):
        return float(self.par.TRTOL)

    def _lte_abstol_vector(self):
        """The absolute floor of the LTE tolerance, per unknown.

        Note this is the `lte_*` pair, NOT the Newton `abstol` -- the two were
        split by stage 0.3d precisely because they are different quantities, and
        the coupled path must take the step-control flavour like every other
        controller here.
        """
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        return self.toolkit.concatenate(
            (self.par.lte_vabstol * ones_nodes,
             self.par.lte_iabstol * ones_branches))

    def _newton_abstol_vector(self):
        """The Newton RESIDUAL tolerance, per unknown.

        A node row of the residual is a current and a branch row is a voltage,
        hence `iabstol` on nodes and `vabstol` on branches.  This is the flavour
        for testing ``f``, not for testing an increment ``dx`` -- see
        :meth:`_newton_xtol_vector`, and `_newton` which builds both.
        """
        from pycircuit.circuit.analysis import newton_tolerance_vectors
        return newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[0]

    def _newton_xtol_vector(self):
        """The Newton SOLUTION tolerance, per unknown -- the other flavour.

        An increment on a node is a voltage and on a branch is a current, so the
        two vectors are transposed with respect to each other.  Getting this
        backwards is the same class of error stage 0.3d separated for the LTE
        tolerances: the numbers are dimensionally different quantities, with
        different defaults (`iabstol` 1e-12, `vabstol` 1e-6 -- which is what
        makes a swap visible).

        History: `doc/transient_history.md`, `Transient._newton_xtol_vector`.
        """
        from pycircuit.circuit.analysis import newton_tolerance_vectors
        return newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[1]
