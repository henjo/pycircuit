"""The stage predictor every implicit family starts its Newton from.  A theme
of `Transient` (see `transient.py`).
"""

import math
import os
import struct

import numpy as np
from numpy._core import umath as _um

from pycircuit.circuit import _paths

#: THE WEIGHTS MEMO (speed round 4, stage D; 2026-10-02): the Vandermonde
#: solve of `_fit` is a pure function of the normalised node times, and
#: the shooting re-walks the same grid every iteration.  False recomputes
#: it every time, for the byte-identity test and as an escape.
PRED_WEIGHT_MEMO = True
PRED_WEIGHT_MEMO_SIZE = 512

#: THE MULTISTEP FAST PATH (speed round 8, stage 2b; 2026-10-04): a
#: prediction with no extra nodes from a history wholly behind the target
#: -- every multistep step, Radau's stages -- in 90 k instructions where the
#: general path takes 168 k.  Bit for bit the general path: behind the
#: target, nearest-first IS newest-first, so the order and the 1e-13 dedupe
#: depend on the history's times alone and are kept per history list
#: (checked per call: the list, its length, each entry -- tests edit it in
#: place); the fit's times in Python floats (the same IEEE subtractions and
#: divisions numpy makes elementwise); the weights memo under the same key
#: bytes; the clamp the same `clip` ufunc, without `np.clip`'s Python
#: wrapper.  A history not strictly ascending (one a test edited) takes the
#: general path.  `PRED_FAST` (env `PYCIRCUIT_PRED_FAST=0`) is read per call.
PRED_FAST = os.environ.get('PYCIRCUIT_PRED_FAST', '1') != '0'

_PC = _paths.COUNTS
_CLIP = _um.clip
#: what the fast path returns to send a call to the general path
_GENERAL = object()


class _StagePredictor:
    """The stage predictor every implicit family starts its Newton from.  A
    theme of `Transient` (see `transient.py`)."""

    ## ------------------------------------------------------------------
    ## THE STAGE PREDICTOR (all families)
    ##
    ## Every implicit integrator here starts its Newton from a value at the
    ## WRONG TIME.  MEASURED on a state-free exponential, as a fraction of one
    ## step's own state motion: coupled Radau IIA(3) 1.00 (`x_n` for all three
    ## stages at once), ESDIRK43 0.50 and TR-BDF2 0.59 (the previous stage),
    ## Gear-2 and trapezoidal 1.00 (`x_n`, one solve per step).  The first
    ## iterations of every solve are then spent travelling `O(h)` rather than
    ## converging, and the stages ARE the per-step work.
    ##
    ## One predictor serves all of them because all of them ask the same
    ## question: what is the trajectory at time `t`?  Nodes are ABSOLUTE TIMES
    ## -- accepted states and, on the sequential paths, this step's already
    ## converged stages -- so a VARIABLE STEP needs no special case, which the
    ## per-tableau formulation this replaces could not do.
    ## `'on'` (the default) or `'off'`, the pre-predictor seed, which is the
    ## control every measurement of this feature is made against
    stage_predictor = 'on'
    ## fit order; None means the method's own order, capped at 4
    stage_predictor_degree = None
    PRED_HIST_MAX = 8
    ## how many linear-extrapolation displacements a prediction may take;
    ## MEASURED across 1.0 / 1.5 / 2.0 / 3.0 / unbounded (see
    ## `benchmarks/stage_predictor.py`) -- 1.5 is where the gain has
    ## arrived on every family and the worst case has not yet started to
    ## grow (Gear-2 on a coarse grid: +3.0% at 1.5, +10.4% at 3.0)
    PRED_CLAMP = 1.5

    def _pred_reset(self):
        """Forget the predictor's node history.  The transient start and a
        shooting period seam, where the trajectory is DISCONTINUOUS and a node
        from before it is a silently wrong guess, not an error."""
        self._pred_hist = []
        self._pred_pending = None

    def _pred_promote(self, x):
        """Turn the last step's pending record into predictor nodes.

        ⚠ Called from the ACCEPT site (through :meth:`_push_history`) and
        nowhere else, so a REJECTED step's stages never become nodes -- they
        are samples of a trajectory the run then threw away.  Every family
        accepts through it -- a stage method reads no charge rings, but its
        idtmod rows need the periodic shifts.

        History: `doc/transient_history.md`, `Transient._pred_promote`.
        """
        pend = getattr(self, '_pred_pending', None)
        if pend is not None:
            self._pred_note(pend[0], x, pend[1])
            self._pred_pending = None

    def _pred_note(self, t, x, stages=()):
        """Record one accepted point, and the stages that produced it, as
        predictor nodes.  Called from :meth:`_push_history` -- the single
        ACCEPT choke point, so a rejected step's stages never become nodes."""
        h = getattr(self, '_pred_hist', None)
        if h is None:
            h = self._pred_hist = []
        ## ⚠⚠ COPY, DO NOT VIEW.  `np.asarray` on an array that is already
        ## float64 returns THE SAME OBJECT, so these entries would alias the
        ## live state and stage vectors -- and the periodic gauge shift
        ## subtracts `n*modulus` from every live history it knows about, so an
        ## aliased entry takes the shift TWICE and the accepted state is
        ## corrupted.  The bookkeeping runs with the predictor switched OFF
        ## too, so 'off' is a control for the SEED, not for this.
        ## History: `doc/transient_history.md`, `Transient._pred_note`.
        for ts, ys in stages:
            h.append((float(ts), np.array(ys, dtype=float)))
        h.append((float(t), np.array(x, dtype=float)))
        ## keep the newest nodes; ties (a stiffly accurate method's last stage
        ## IS the step) are harmless, the nearest-node pick drops duplicates
        h.sort(key=lambda e: e[0])
        uniq = []
        for e in h:
            if uniq and abs(e[0] - uniq[-1][0]) <= 1e-14 * max(abs(e[0]), 1.0):
                uniq[-1] = e
            else:
                uniq.append(e)
        self._pred_hist = uniq[-self.PRED_HIST_MAX:]

    def _pred_degree(self):
        """How many nodes the prediction fits.  The method's own order, so the
        predictor is as accurate as what it is predicting for; capped at 4,
        which is where the measured gain stops improving."""
        d = self.stage_predictor_degree
        if d:
            return int(d)
        return min(4, max(1, int(getattr(self.base_integrator, 'order', 2))))

    def _pred_or(self, fallback, ttarget, extra=()):
        """:meth:`_predict_state`, or ``fallback`` when it declines."""
        p = self._predict_state(ttarget, extra=extra)
        return fallback if p is None else p

    def _predict_state(self, ttarget, extra=(), deg=None):
        """The trajectory at ``ttarget``, from the ``deg + 1`` recorded nodes
        NEAREST it, or ``None`` when there are too few to fit a line.

        ``extra`` is ``(t, x)`` pairs known only within the current step -- the
        already converged stages of a sequential method.  They are nodes like
        any other, and being the nearest ones they are what turns a whole-step
        extrapolation into a one-stage-gap one.

        ⚠⚠ THE CLAMP IS NOT OPTIONAL.  A polynomial continued past its last
        node can leave the region the circuit actually visits, and on an
        exponential device a 3x overshoot is ``exp(3 dV / VT)``: the worst
        stage then costs several times the old seed's Newton iterations while
        the MEAN still improves -- which is how such a heuristic passes its
        own gate and fails in use.  The prediction is therefore confined
        COMPONENTWISE to the range its own nodes span, widened by the motion
        between the two newest.  The scale comes from the data.

        History: `doc/transient_history.md`, `Transient._predict_state`.
        """
        if self.stage_predictor == 'off':
            return None
        if PRED_FAST and not extra:
            out = self._predict_fast(ttarget, deg)
            if out is not _GENERAL:
                return out
        nodes = [(float(tt), np.asarray(xx, dtype=float)) for tt, xx in extra]
        nodes.extend(getattr(self, '_pred_hist', ()) or ())
        if len(nodes) < 2:
            return None
        if deg is None:
            deg = self._pred_degree()
        ## (speed round 4, stage D: the dedupe as a plain loop with the same
        ## rule and the same first-wins order, one sort of the newest, and
        ## the weights served from `_fit`'s memo -- 26.5 -> 14 us a step on
        ## a hit, 22.5 on a miss, bit-identical)
        ## ⚠ DUPLICATE TIMES MAKE THE VANDERMONDE SINGULAR, and they are the
        ## normal case here, not a corner one: a stiffly accurate method's
        ## last stage IS the step it ends, and an ESDIRK's explicit first
        ## stage IS the state it starts from, so `extra` routinely repeats a
        ## recorded node.  Left in, the solve raises and the predictor
        ## declines -- measured, 37% of ESDIRK43's stages silently kept the
        ## old seed.  `extra` is listed first so it wins a tie, being the
        ## value from the step actually in progress.
        nodes.sort(key=lambda e: abs(e[0] - ttarget))
        uniq = []
        for e in nodes:
            t0 = e[0]
            tol = 1e-13 * max(abs(t0), 1.0)
            for u in uniq:
                if abs(t0 - u[0]) <= tol:
                    break
            else:
                uniq.append(e)
        nodes = uniq
        if len(nodes) < 2:
            return None
        nfit = min(int(deg) + 1, len(nodes))
        take = nodes[:nfit]

        newest_first = sorted(take, key=lambda e: -e[0])
        if len(newest_first) < 2:
            return None
        motion = np.abs(newest_first[0][1] - newest_first[1][1])
        ## ⚠ NO SELF-VALIDATION GATE (retrodict the newest node from the older
        ## ones, decline on a miss): it changes nothing measurable, because it
        ## looks BACKWARD, and the case it would be for is a knee that has not
        ## happened yet.  What actually bounds the damage is the clamp below.
        ## History: `doc/transient_history.md`, `Transient._predict_state`.
        pred = self._fit(take, ttarget)
        if pred is None:
            return None
        ## ⚠⚠ THE CLAMP, AND IT IS THE WHOLE DIFFERENCE BETWEEN A SPEED-UP AND
        ## A REGRESSION.  A polynomial continued past its last node can leave
        ## the region the circuit visits, and on an exponential device a 3x
        ## overshoot is `exp(3 dV / VT)`.
        ##
        ## The bound is the LINEAR prediction: from the newest node, moving at
        ## the rate the last step moved, the target is `ratio * motion` away.
        ## A higher-order fit is allowed a multiple of that -- real curvature
        ## needs the room -- but not an unbounded one.  Both scales come from
        ## the data: `motion` is the last step's own displacement and `ratio`
        ## is how far ahead this target is in units of the last step.
        xref = newest_first[0][1]
        dt_last = abs(newest_first[0][0] - newest_first[1][0])
        ratio = (abs(ttarget - newest_first[0][0]) / dt_last) if dt_last > 0 \
            else 1.0
        w = self.PRED_CLAMP * motion * max(ratio, 1.0)
        out = np.clip(pred, xref - w, xref + w)
        ## ⚠⚠ A WRAPPING STATE IS NOT A TRAJECTORY THIS CAN FIT.  A periodic
        ## row folds by its modulus, and that is a DISCONTINUITY in exactly the
        ## curve a polynomial is being put through -- the gauge shift keeps the
        ## recorded nodes in one gauge, but the fold can also fall between the
        ## newest node and the target, and then the fit runs straight across
        ## it (with the wrap exactly ON a grid point, predicting these rows
        ## stops the shooting Newton converging at all).  Those rows keep the
        ## old seed -- the newest node's value -- and every other row still
        ## gets the prediction.
        ## History: `doc/transient_history.md`, `Transient._predict_state`.
        for row, _m, _o in (getattr(self, '_periodic_rows', None) or ()):
            out[row] = xref[row]
        return out

    def _predict_fast(self, ttarget, deg):
        """`_predict_state` with no extra nodes, from a history wholly
        behind `ttarget` (`PRED_FAST`), or `_GENERAL`: the general path."""
        hist = getattr(self, '_pred_hist', None) or ()
        if len(hist) < 2:
            return None
        if type(hist) is not list or not hist[-1][0] < ttarget:
            _PC['pred:general'] += 1
            return _GENERAL
        ## (a `_fit` patched on the class or shadowed on the instance is a
        ## caller that expects its calls -- the core's rule for its methods)
        if type(self)._fit is not _FIT or '_fit' in self.__dict__:
            _PC['pred:patched'] += 1
            return _GENERAL
        found = self._pred_fast_nodes(hist)
        if found is None:
            _PC['pred:general'] += 1
            return _GENERAL
        nodes, gap = found
        ## (nearest-first is newest-first unless two distances to the target
        ## ROUND EQUAL: the general path's stable sort then lists the older
        ## of the two first, its dedupe keeps that one, and its fit takes the
        ## nodes in that order.  A tie needs two times closer than the
        ## largest distance's rounding (2.2e-16 of it; this test has a margin
        ## of two) -- never the transient's own history, whose times are
        ## 1e-14 apart, but a history a caller wrote: the deep hypothesis
        ## profile drew 0 and 1.6e-45 behind a target 1.4e-3 away, 2026-10-04)
        if not gap > 4.5e-16 * (ttarget - hist[0][0]):
            _PC['pred:tie'] += 1
            return _GENERAL
        _PC['pred:fast'] += 1
        if len(nodes) < 2:
            return None
        if deg is None:
            deg = self._pred_degree()
        take = nodes[:min(int(deg) + 1, len(nodes))]
        if len(take) < 2:
            return None
        ## (`take` is newest-first: the general path's `newest_first`)
        motion = np.abs(take[0][1] - take[1][1])
        pred = self._fit_fast(take, ttarget)
        if pred is None:
            return None
        xref = take[0][1]
        dt_last = abs(take[0][0] - take[1][0])
        ratio = (abs(ttarget - take[0][0]) / dt_last) if dt_last > 0 else 1.0
        w = self.PRED_CLAMP * motion * max(ratio, 1.0)
        out = _CLIP(pred, xref - w, xref + w)
        for row, _m, _o in (getattr(self, '_periodic_rows', None) or ()):
            out[row] = xref[row]
        return out

    def _pred_fast_nodes(self, hist):
        """The history nearest-first -- newest-first, every node being
        behind the target and no two distances tied -- with
        `_predict_state`'s 1e-13 dedupe, and the smallest gap between two of
        its times (the tie test), kept per history list: `(nodes, gap)`, or
        None where the times are not strictly ascending."""
        c = self.__dict__.get('_pred_fast')
        if c is not None and c[0] is hist and len(c[1]) == len(hist):
            for a, b in zip(hist, c[1]):
                if a is not b:
                    break
            else:
                return c[2]
        ents = tuple(hist)
        uniq = []
        prev = prev_kept = None
        gap = float('inf')
        ## newest first: every kept node is later than the one at hand, the
        ## nearest of them the last kept -- the general path's dedupe against
        ## all kept nodes is this one comparison
        for e in reversed(ents):
            t0 = e[0]
            if prev is not None:
                if not t0 < prev[0]:
                    uniq = None
                    break
                gap = min(gap, prev[0] - t0)
                if abs(t0 - prev_kept) <= 1e-13 * max(abs(t0), 1.0):
                    prev = e
                    continue
            uniq.append(e)
            prev = e
            prev_kept = t0
        out = None if uniq is None else (uniq, gap)
        self.__dict__['_pred_fast'] = (hist, ents, out)
        return out

    def _fit_fast(self, sub, tat):
        """`_fit` with the node times in Python floats: the same IEEE
        operations, the same memo key bytes, the same solve on a miss."""
        tat = float(tat)
        d = [e[0] - tat for e in sub]
        for v in d:
            if math.isnan(v):
                return None
        scale = max(abs(v) for v in d)
        if not scale < float('inf') or scale <= 0.0:
            return None
        n = len(d)
        tau = [v / scale for v in d]
        memo = self.__dict__.get('_pred_wmemo') if PRED_WEIGHT_MEMO else None
        key = (n, struct.pack(f'{n}d', *tau))
        w = memo.get(key) if memo is not None else None
        if w is None:
            rhs = np.zeros(n)
            rhs[0] = 1.0
            try:
                w = np.linalg.solve(np.vander(np.array(tau), n, increasing=True).T, rhs)
            except np.linalg.LinAlgError:
                return None
            if not np.all(np.isfinite(w)):
                return None
            if PRED_WEIGHT_MEMO:
                if memo is None:
                    memo = self.__dict__['_pred_wmemo'] = {}
                elif len(memo) >= PRED_WEIGHT_MEMO_SIZE:
                    memo.clear()
                memo[key] = w
        return w @ np.array([e[1] for e in sub], dtype=float)

    def _fit(self, sub, tat):
        """The polynomial through ``sub`` evaluated at ``tat``, or None.

        Fitted in a coordinate centred on the TARGET and scaled by the node
        spread: an absolute-time Vandermonde at t ~ 1e-3 with 1e-6 spacing
        is hopeless, and centring makes the right-hand side exactly e_0.

        The weights are a pure function of the normalised times `tau`, so
        they are kept on the instance keyed by their bits (`PRED_WEIGHT_MEMO`;
        bounded, forgotten at `_memo_clear` -- a `solve` starts afresh, the
        shooting's period walks share it): the shooting re-walks the same
        grid every iteration, and an exact uniform grid repeats them too.
        `np.vander` stays -- its `multiply.accumulate` powers are the bits.
        """
        tv = np.array([e[0] for e in sub], dtype=float)
        scale = float(np.max(np.abs(tv - tat)))
        if not np.isfinite(scale) or scale <= 0.0:
            return None
        tau = (tv - tat) / scale
        n = len(tau)
        memo = self.__dict__.get('_pred_wmemo') if PRED_WEIGHT_MEMO else None
        key = (n, tau.tobytes())
        w = memo.get(key) if memo is not None else None
        if w is None:
            rhs = np.zeros(n)
            rhs[0] = 1.0
            try:
                w = np.linalg.solve(np.vander(tau, n, increasing=True).T, rhs)
            except np.linalg.LinAlgError:
                return None
            if not np.all(np.isfinite(w)):
                return None
            if PRED_WEIGHT_MEMO:
                if memo is None:
                    memo = self.__dict__['_pred_wmemo'] = {}
                elif len(memo) >= PRED_WEIGHT_MEMO_SIZE:
                    memo.clear()
                memo[key] = w
        return w @ np.array([e[1] for e in sub], dtype=float)


#: the general fit, as defined: `_predict_fast` declines where it is not
#: the one a call would reach
_FIT = _StagePredictor._fit
