"""The LMM companion model and the per-step device memo.  A theme of
`Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit import _evalhint, _paths, _tran_core

_PC = _paths.COUNTS


def _first_equal(a, b):
    """False where two exact float64 vectors of one shape differ at their
    first entry -- `a == b` cannot then be all true (NaN never equal, signed
    zeros equal, as `==` has them) -- so a lookup that misses is decided at
    one entry: 83 % of the vdP PSS's C lookups miss, after a comparison of
    ~20 k instructions each -- the C lookup's and the charge's (`_q_at`)
    (speed round 12, stage 3).  True otherwise: the
    full comparison follows."""
    if (type(a) is np.ndarray and type(b) is np.ndarray and a.dtype is _F64
            and b.dtype is _F64 and a.ndim == 1 and a.shape[0]):
        return a[0] == b[0]
    return True

#: `u(t)` ONCE PER STEP (speed round 4, stage B; 2026-10-02): the source
#: vector assembled at a time serves every later request at that exact
#: time within the step (`_source_at`).  False is the old behaviour -- the
#: sources re-assembled at every Newton iteration -- for the byte-identity
#: test and as an escape.
U_MEMO = True

_F64 = np.dtype('float64')


def _as_float(a, toolkit):
    """`toolkit.array(a, dtype=float)` without the copy when `a` is already
    a float64 ndarray: the sums a step builds (`i + iq + u`, `G + Geq`) are
    fresh and nobody else holds them (speed round 4, stage C)."""
    if type(a) is np.ndarray and a.dtype is _F64:
        return a
    return toolkit.array(a, dtype=float)


class _CompanionModel:
    """The LMM companion model and the per-step device memo.  A theme of
    `Transient` (see `transient.py`)."""

    ## The integrator is selected by the `integrator` Parameter, not per call.
    ## History: `doc/transient_history.md`, `Transient.get_diff`.
    def get_diff(self, q, C):
        """Method used to calculate time derivative for charge storing elements (i_eq and g_eq)."""
        # Determine the active integrator based on step size variations
        h_last = self._dt_last if self._dt_last is not None else self._dt
        self.active_integrator = self.base_integrator.check_order_drop(
            self._dt, h_last, self._is_first_step
        )
        
        iq, geq = self.active_integrator.compute_derivatives(
            q_curr=q, C_curr=C, h_curr=self._dt, 
            q_last=self._qlast, iq_last=self._iqlast, h_last=h_last,
            is_first_step=self._is_first_step,
            toolkit=self.toolkit
        )
        
        self._iq = iq
        ## The companion CONDUCTANCE, stored beside the companion current for
        ## the same reason: a caller that needs the per-step sensitivity
        ## `Jf^-1 Geq` -- shooting's monodromy -- cannot recompute it without
        ## repeating the whole assembly.  `_iq` has been kept here since
        ## stage 11; this is its other half.
        self._Geq = geq
        ## And the two things that turn one step into a DERIVATIVE: the
        ## capacitance matrix this step used, and the companion coefficients
        ## of the integrator that actually ran -- `active_integrator`, not
        ## `base_integrator`, so an order drop is reflected rather than
        ## assumed away.  A caller differentiating the step (shooting) needs
        ## the past `C` matrices as well as the current one, which `_Geq`
        ## alone cannot supply once a method reaches back more than one step.
        self._Cmat = C
        self._companion_coeffs = self.active_integrator.companion_coefficients(
            self._dt, h_last)
        ## The class actually running this step -- check_order_drop() may have
        ## dropped to a lower-order integrator than the one requested, so this
        ## is derived from the live object rather than compared against a
        ## removed user-supplied method name.
        self._effective_method = type(self.active_integrator).__name__
        return iq, geq

    def _memo_clear(self):
        """Forget the device-evaluation memo (`_memo_get`) and the 'auto'
        Newton options' verdict on the circuit (`_newton_option`): at the
        start of every solve, the one place a caller can change the circuit
        between two steps."""
        self._dev_memo = ({}, {})
        self._memo_rolling = False
        self._jacobian_expensive = None
        self._u_memo = None
        self._pred_wmemo = None

    def _memo_step(self):
        """A new step: the current generation becomes the previous one (the
        step's START is the previous step's last stage)."""
        memo = getattr(self, '_dev_memo', None) or ({}, {})
        self._dev_memo = ({}, memo[0])
        self._memo_rolling = True

    def _memo_put(self, x, rec):
        """Record the device evaluations `rec` (some of `q`, `i`, `C`, `G`)
        at `x`, merged into what is recorded there already (a record the
        previous step holds is carried into this one)."""
        memo = getattr(self, '_dev_memo', None)
        if memo is None:
            memo = self._dev_memo = ({}, {})
        ## (an exact float64 array is its own `asarray`: 0.6 k for the key
        ## where `asarray` and `tobytes` cost 1.7 k -- speed round 12)
        key = (x if type(x) is np.ndarray and x.dtype is _F64
               else np.asarray(x, dtype=float)).tobytes()
        cur = memo[0].get(key)
        if cur is None:
            cur = memo[0][key] = dict(memo[1].get(key, ()))
        cur.update(rec)

    def _memo_ok(self):
        """Whether an evaluation may be recorded (`_memo_put`): on a path
        whose memo rolls per step (the coupled stage steps, `_memo_step`;
        elsewhere it would only grow), with no stateful limiter and no
        bypass -- the guards the dense stage Newton records under."""
        return (bool(getattr(self, '_memo_rolling', False))
                and not getattr(self, '_stateful_lims', None)
                and float(getattr(self.epar, 'bypasstol', -1.0) or -1.0) < 0.0)

    def _memo_get(self, x):
        """The device evaluations (`q`, `i`, `C`, `G`, full width) the coupled
        stage Newton made AT this exact state in this step or the previous
        one, or None.  A MEMOISATION, as `_q_at`: keyed on the state's bytes,
        recorded only without a shunt and without a stateful limiter, so the
        value is the one recomputing gives, bit for bit.  Radau re-read `C`
        and `G` at states its Newton had just evaluated -- 28 % of a compact
        MOSFET's PSS solve (2026-09-30)."""
        memo = getattr(self, '_dev_memo', None)
        if memo is None or x is None:
            return None
        if not memo[0] and not memo[1]:
            ## (the multistep path never records: no key to build)
            return None
        key = (x if type(x) is np.ndarray and x.dtype is _F64
               else np.asarray(x, dtype=float)).tobytes()
        rec = memo[0].get(key)
        return memo[1].get(key) if rec is None else rec

    def _C_at_state(self, x):
        """``cir.C(x)``, reusing the last assembly's (`_companion_at`) when it
        was at this state -- the branch screen reads `C` at the converged
        point, where `jacobian_only` has just assembled it: 360 of a compact
        MOSFET's 1500 `C` evaluations on a 40-point gear PSS (2026-09-30).
        A MEMOISATION, as `_q_at`: identity, then full equality; and never
        across a stateful limiter (`Diode`), whose device reads its stored
        state as well as `x`."""
        C = self._C_lookup(x)
        if C is None:
            C = self.cir.C(x, self.epar)
            ## (recorded where the memo rolls: the stage paths' branch screen
            ## reads `C` at each stage, and the shooting's stage step then
            ## reads it there again -- `_C_at`)
            if self._memo_ok():
                self._memo_put(x, {'C': C})
        return C

    def _C_lookup(self, x):
        """The `C` `_C_at_state` would serve at `x` without evaluating, or
        None -- what an evaluation session asks before naming `C`
        (`_evalhint`: a session names exactly what will be evaluated)."""
        rec = self._memo_get(x)
        if rec is not None and 'C' in rec:
            _PC['memo.C:memo'] += 1
            return rec['C']
        cached = getattr(self, '_C_cache', None)
        if (cached is not None and not getattr(self, '_stateful_lims', None)
                and float(getattr(self.epar, 'bypasstol', -1.0) or -1.0) < 0.0):
            x_cached, C_cached = cached
            if x_cached is x:
                _PC['memo.C:same'] += 1
                return C_cached
            if (x_cached is not None and x is not None
                    and getattr(x_cached, 'shape', None) == getattr(x, 'shape', None)
                    and _first_equal(x_cached, x)
                    and bool(self.toolkit.alltrue(x_cached == x))):
                _PC['memo.C:equal'] += 1
                return C_cached
        return None                             # paths: not a decline (a lookup's miss)

    def _q_at(self, x):
        """``cir.q(x)``, reusing the value computed during the last assembly.

        STAGE 2c.  This is a *memoisation*, not an approximation: the cached value
        was produced by the same function at the same state, so it is bit-identical
        to recomputing, and the whole of stage 2 is defined as behaviour-preserving.
        The guard is deliberately strict -- identity first, then full equality --
        because serving a charge vector from the wrong state would corrupt the LTE
        estimate silently, which is precisely the failure class stage 1 removed.
        """
        cached = getattr(self, '_q_cache', None)
        if cached is not None:
            x_cached, q_cached = cached
            if x_cached is x:
                _PC['memo.q:same'] += 1
                return q_cached
            if (x_cached is not None and x is not None
                    and getattr(x_cached, 'shape', None) == getattr(x, 'shape', None)
                    and _first_equal(x_cached, x)
                    and bool(self.toolkit.alltrue(x_cached == x))):
                _PC['memo.q:equal'] += 1
                return q_cached
        _PC['memo.q:miss'] += 1
        return self.cir.q(x, self.epar)

    def _companion_at(self, x, C=None):
        """``(iq, Geq)``: the step's companion current and conductance at `x`
        (the current ``self._dt``), with the charge cached against the state
        it belongs to.  One assembly for every step's residual and Jacobian.
        `C` is the capacitance a caller has already looked up at `x`
        (`_C_lookup`), so the lookup is not made twice.

        `self.epar`, not the module-level `defaultepar`: without it every
        device is evaluated at defaultepar's T = 300 K whatever the caller
        asked for, and -- because `Analysis.__init__` attaches `bypasstol` to
        the analysis's own epar and nowhere else -- the `bypass` parameter
        does nothing at all.

        STAGE 2c.  The charge vector is stashed alongside the state it belongs
        to.  `solve()` needs `q` at the converged point twice more -- once for
        the step controller and once for the history roll -- and `_q_at`
        serves both from here instead of repeating the assembly.
        Keyed by the state so a stale value can never be served: the check is
        identity-then-equality on x, not a bare "did we cache".

        History: `doc/transient_history.md`, `Transient._companion_at`."""
        ## (`C` the branch screen may just have read at this state:
        ## `_C_at_state`, a memoisation)
        if C is None:
            C = self._C_at_state(x)
        q = self.cir.q(x, self.epar)
        self._q_cache = (x, q)
        self._C_cache = (x, C)
        return self.get_diff(q, C)

    def _source_at(self, t, provided_function=None):
        """`u(t)`: the circuit's sources, plus `provided_function(t)`.

        ONE CONTRACT: `provided_function(t)` is an extra source term, on every
        path (F4).  A caller written for a post-solve callback
        `provided_function(f, J, C)` breaks loudly on arity.

        ONCE PER STEP AND TIME (speed round 4, stage B; 2026-10-02).  The
        sources are a function of `t`, `epar` and the analysis name, and
        none of them moves inside a step: nothing on the numeric path
        writes `epar.t` (it is only read), `analysis_kind` is scoped around
        a whole solve, and an element's state moves only at `accept_step`
        / `reset_state`.  So the vector assembled at the first request at
        a time serves every later request at that exact `t` -- the
        Newton's iterations, the chord's residual-only ones, the branch
        confirmation's re-solve, a coupled method's repeated stage times --
        and the memo lives exactly as long as the step (`solve_timestep`
        opens and closes it; outside one there is none).  The same
        function on the same inputs: the same bits.  `provided_function(t)`
        is called every time, as before (a caller may count it), and the
        sum below is a new vector.  Measured: a gear step assembled `u`
        2.03x on the PSP stage (31 us a call in-run), 1.89x on a 20-PSP
        chain (126 us).  `U_MEMO` switches it off.

        History: `doc/transient_history.md`, `Transient._source_at`."""
        analysis = self.par.analysis
        memo = self.__dict__.get('_u_memo') if U_MEMO else None
        key = (t, analysis)
        u = None
        if memo is not None:
            try:
                u = memo.get(key)
            except TypeError:
                ## a time that cannot be a key -- a JAX array under that
                ## toolkit, a 0-d numpy array: assembled every time, as before
                memo = None
                _PC['umemo:unhashable'] += 1
        if u is None:
            u = self.cir.u(t, self.epar, analysis=analysis)
            if memo is not None:
                _PC['umemo:miss'] += 1
                memo[key] = u
        elif memo is not None:
            _PC['umemo:hit'] += 1
        if provided_function is not None:
            u = u + provided_function(t)
        return u

    def _residual_and_jacobian(self, x, t, provided_function=None):
        """``(f, J)`` at ``(x, t)`` using the current ``self._dt`` -- the step
        residual `solve_timestep`'s Newton drives to zero, and reachable
        without one: the coupled method needs one residual per iteration of
        its OWN loop.

        One evaluation session (`_evalhint`): a compiled model computes the
        `q`, `i`, `G` (and `C`, unless a cache serves it) it is asked for
        here in one fused pass -- the same bits, a compact MOSFET's four
        evaluations at a state 6.3 -> 3.4 ms."""
        ## THE EVALUATE CORE (`_tran_core`, speed round 7): the passes, the
        ## companion and the residual in one C call where it serves, the
        ## same state left behind; None, and the path below, where not
        r = _tran_core.evaluate(self, x, t, provided_function, 'fj')
        if r is not None:
            return r
        C = self._C_lookup(x)
        need = ('q', 'i', 'G') if C is not None else ('C', 'q', 'i', 'G')
        with _evalhint.evaluating(*need):
            iq, Geq = self._companion_at(x, C)
            u = self._source_at(t, provided_function)
            f = self.cir.i(x, self.epar) + iq + u
            J = self.cir.G(x, self.epar) + Geq
        return _as_float(f, self.toolkit), _as_float(J, self.toolkit)
