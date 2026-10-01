"""Breakpoints, state events and the excursion veto of the stepping loop.  A
theme of `Transient` (see `transient.py`).
"""

import numpy as np


class _StepEvents:
    """Breakpoints, state events and the excursion veto of the stepping loop.
    A theme of `Transient` (see `transient.py`)."""

    ## E7: a state event is landed when the crossing sits within this
    ## fraction of the step's end, in at most this many secant re-solves.
    ## Measured on the comparator oscillator (gear, reltol 1e-6): 1e-6 / 6
    ## took 2.6 cuts per landing, 1e-3 / 4 takes 2.1 for the same
    ## period-to-period spread (1.6e-5 against 1.3e-4 unlanded); one cut
    ## alone leaves the spread where it was.
    EVENT_LAND_RTOL = 1e-3
    EVENT_LAND_MAXITER = 4
    ## whether a landed crossing restarts the multistep history (measured
    ## both ways on the comparator oscillator -- see the E7 entry)
    EVENT_RESTART_HISTORY = True

    def _source_signal_scales(self):
        """(v_scale, i_scale): the largest swing any source declares, via
        the signal_scale element hook, walked recursively."""
        vs, is_ = 0.0, 0.0
        def walk(circuit):
            nonlocal vs, is_
            for element in getattr(circuit, 'elements', {}).values():
                hook = getattr(element, 'signal_scale', None)
                if hook is not None:
                    v, i = hook()
                    vs, is_ = max(vs, float(v)), max(is_, float(i))
                walk(element)
        walk(self.cir)
        return vs, is_

    def _dv_step_bounds(self):
        """((bv_static, cv_rel), (bi_static, ci_rel)) for the excursion
        check: effective bound = max(static, rel * running unit-group max).

        Manual factor f: static = max(f, 1) * lte_abstol (the LTE family's
        abstols -- the owner's rationale is that the LTE scales with these,
        so this check must too; clamp-at-1 = the solver-noise floor),
        rel = 0.  'auto' (sampling theory, owner request): rel = 2*pi/N for
        N = points_per_period, static = max(rel * source swing,
        lte_abstol) -- source-anchored at signal birth, where a running
        reference h-cancels.  None: (inf, 0), disabled."""
        import math
        c_rel = 2.0 * math.pi / float(self.par.points_per_period)
        v_src, i_src = (self._source_signal_scales()
                        if 'auto' in (self.par.max_dv_step,
                                      self.par.max_di_step) else (0.0, 0.0))

        def resolve(knob, abstol, src):
            if knob is None:
                return float('inf'), 0.0
            if knob == 'auto':
                return max(c_rel * src, abstol), c_rel
            return max(float(knob), 1.0) * abstol, 0.0

        bv = resolve(self.par.max_dv_step,
                     float(self.par.lte_vabstol), v_src)
        bi = resolve(self.par.max_di_step,
                     float(self.par.lte_iabstol), i_src)
        return bv, bi

    def _excursion_ratio(self, x_new, x_prev):
        """The excursion check (`max_dv_step` / `max_di_step`): the largest
        per-step change over the node rows against the voltage bound, and
        over the branch rows against the current bound -- a ratio above 1
        vetoes the step, and a proportional retry (``0.9 / ratio``) follows
        (the industry semantics).  None when both knobs are off.

        The running unit-group maxima behind 'auto' are updated from
        `x_prev`, the ACCEPTED state, before this candidate is looked at --
        so the relative term cannot h-cancel at a signal birth.

        The one stepping loop applies it to every family.

        History: `doc/transient_history.md`, `Transient._excursion_ratio`."""
        if self.par.max_dv_step is None and self.par.max_di_step is None:
            return None
        (_bvs, _cvr), (_bis, _cir_) = self._dv_step_bounds()
        _nn = len(self.cir.nodes)
        _xa = np.abs(np.asarray(x_prev, dtype=float))
        self._dv_run_v = max(getattr(self, '_dv_run_v', 0.0),
                             float(np.max(_xa[:_nn])))
        self._dv_run_i = max(getattr(self, '_dv_run_i', 0.0),
                             float(np.max(_xa[_nn:]))
                             if _nn < len(_xa) else 0.0)
        _bv = max(_bvs, _cvr * self._dv_run_v)
        _bi = max(_bis, _cir_ * self._dv_run_i)
        _d = np.abs(np.asarray(x_new, dtype=float)
                    - np.asarray(x_prev, dtype=float))
        ratio = float(np.max(_d[:_nn])) / _bv
        if _nn < len(_d):
            ratio = max(ratio, float(np.max(_d[_nn:])) / _bi)
        return ratio

    def _next_breakpoint(self, t):
        """The next source breakpoint after `t`, strictly advancing: a
        corner closer than `minbreak` (relative) is skipped so a step of
        `dt = 0` cannot loop."""
        nb = self.cir.next_event(t)
        if nb <= t + self.par.minbreak * max(abs(t), 1.0):
            nb = self.cir.next_event(t + (self.par.minbreak * 1e3) * max(abs(t), 1.0))
        return nb

    def _state_event_step(self, x_prev, x_new, h, t_end, cuts, minstep):
        """E7: a declared crossing inside the accepted step `x_prev -> x_new`
        (length `h`, ending at `t_end`)?  Returns `('cut', f)` when the step
        should be re-solved at `f * h` -- the secant on the fraction, `cuts`
        re-solves already taken -- `'landed'` once the crossing sits within
        `EVENT_LAND_RTOL` of the step's end (the time recorded in
        `event_times`), or `None` with no crossing.  One hook for every
        stepping loop; what a landing does to the history is the loop's
        business (a multistep method restarts it, a one-step one has none,
        the coupled solve had held its step while cutting)."""
        if self._ev_rows is None:
            return None
        w = self._ev_rows.shape[1]
        s0 = self._ev_rows @ np.asarray(x_prev, dtype=float)[:w] - self._ev_thr
        s1 = self._ev_rows @ np.asarray(x_new, dtype=float)[:w] - self._ev_thr
        cross = np.flatnonzero(s0 * s1 < 0.0)
        if not len(cross):
            return None
        f = float(np.min(s0[cross] / (s0[cross] - s1[cross])))
        if (f < 1.0 - self.EVENT_LAND_RTOL and cuts < self.EVENT_LAND_MAXITER
                and f * h >= minstep):
            self.statistics.state_event_cuts += 1
            return 'cut', f
        self.event_times.append(float(t_end))
        self.statistics.state_events_hit += 1
        return 'landed'

    def _init_state_events(self, n):
        """E7 (2026-09-22): the circuit's declared state events as rows on
        the FULL state (node voltages first; the branch currents that follow
        get zero weight), and an empty `event_times` for the run."""
        self.event_times = []
        self._ev_rows, self._ev_thr = None, None
        if self.par.state_events and hasattr(self.cir, 'state_events'):
            _rows = self.cir.state_events()
            if _rows:
                self._ev_rows = np.array([np.pad(np.asarray(r, dtype=float).ravel(),
                                                 (0, max(0, n - len(np.asarray(r).ravel()))))
                                          for r, _t in _rows])
                self._ev_thr = np.array([float(t_) for _r, t_ in _rows])
