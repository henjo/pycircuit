from .toolkit import numeric
from scipy import interpolate

class TimeFunction():
    """Time dependent function"""
    
    def __init__(self, toolkit=numeric):
        self.toolkit = toolkit
    
    def f(self, t):
        return 0        
        
    def next_event(self, t):
        """Returns the time of the next event given the current time t"""
        return self.toolkit.inf

    def dfdt(self, t):
        """Time derivative of the source, for Fang's ``p = df_ckt/dh``.

        STAGE 12B.  The circuit residual is evaluated at ``t_{m-1} + h``, so
        differentiating it with respect to the step size needs ``du/dt``.  With
        the LTE made an unknown (DAC 2013 eq 11) this is no longer optional --
        the source term is often the LARGEST part of ``p`` on a driven circuit,
        and dropping it does not make the coupled system slightly wrong, it makes
        it solve a different problem.

        The default is a central finite difference, so a source that does not
        override this still works.  Subclasses with a closed form override it:
        a finite difference costs two extra evaluations and, more importantly,
        is meaningless across the kinks of a piecewise source -- which is exactly
        where the step size is under the most pressure.

        **The right-hand limit at a corner is not a compromise, it is the only
        value ever asked for.**  `Pulse` and `PWL` are piecewise linear, so
        `du/dt` genuinely does not exist at their corners -- but every one of
        those corners is reported by :meth:`next_event`, so the step is truncated
        to land exactly on it and the breakpoint re-arms `_is_first_step`, which
        drops the following step to backward Euler.  A step therefore always
        STARTS at a corner and moves forward into the segment ahead.  That is the
        segment these derivatives report.

        The invariant this depends on -- every discontinuity in `du/dt` is a
        reported breakpoint -- is asserted in
        `test_source_derivatives.py::test_every_kink_is_reported_as_a_breakpoint`
        rather than assumed, because a source that grew a new kink without
        reporting it would hand `p` a derivative from the wrong segment and
        nothing else would notice.
        """
        eps = 1e-9 * max(abs(t), 1.0)
        return (self.f(t + eps) - self.f(t - eps)) / (2.0 * eps)


class Sin(TimeFunction):
    def __init__(self, offset=0, amplitude=0, 
                 freq=0, td=0, 
                 theta=0, phase=0,
                 toolkit=numeric):
        self.offset = offset
        self.amplitude = amplitude
        self.omega = 2 * toolkit.pi * freq
        self.phase = phase * toolkit.pi / 180
        self.theta = theta
        self.td = td
        self.toolkit = toolkit
    
    def next_event(self, t):
        """The start of the source, and nothing after it.

        STAGE 4g(a).  This used to return an event at every peak and every zero
        crossing -- one per quarter period, forever.  **A breakpoint means a
        discontinuity**, and a sine has none: it is C-infinity everywhere after
        ``td``.  Peaks and zero crossings are not discontinuities, they are just
        places where a derivative happens to vanish.

        The cost of pretending otherwise was not cosmetic.  The transient truncates
        the step to land exactly on each breakpoint and then re-arms an order drop
        across it, so a plain sine drive produced a periodic, drive-synchronous
        order drop -- and for the trapezoidal rule each restart re-seeds the
        alternating `(-1)^n` mode in the companion current.  It also forced a
        step-size truncation four times per period regardless of what the error
        estimate wanted.

        ``td`` IS a real event: the waveform is only C-infinity *after* the source
        starts, and at ``td`` the derivative is generally discontinuous.  So return
        it once, then nothing.
        """
        ## Under a symbolic toolkit `t` is a sympy symbol and `t < td` is a
        ## relational, not a bool -- asking for its truth value raises.  There is
        ## no meaningful breakpoint schedule for a symbolic time anyway, so fall
        ## through to "no further events".
        try:
            if bool(t < self.td):
                return self.td
        except TypeError:
            pass
        return self.toolkit.inf
        
    def f(self, t):
        # --- DAMPED SINE WAVE EQUATION ---
        # V(t) = VO + VA * exp(-theta * (t - td)) * sin(omega * (t - td) + phase)
        #
        # STAGE 8(b) -- CLAMPED AT td, AND THAT IS A CORRECTNESS FIX, NOT A TIDY-UP.
        #
        # This used to evaluate `t - td` unclamped, which is wrong twice over:
        #
        #   * the sine ran from t = 0 instead of from `td` -- measured at 0.9511 of
        #     full amplitude at t = 0.2*td, where SPICE holds the source at its
        #     offset;
        #   * and for t < td the exponent `-theta*(t - td)` is POSITIVE, so the
        #     "damping" term GREW backwards in time.  Measured with theta = 5000 and
        #     a 1 V amplitude: **f(2e-4) = 51.93 V**, and it is unbounded in
        #     theta*td -- the review recorded 2835 V from a 1 V source.
        #
        # Clamping `t - td` at zero gives SPICE's rule for free: at dt = 0 the
        # expression collapses to `VO + VA*sin(phase)`, which is exactly what SPICE
        # specifies before the source starts.  No second branch is needed, so there
        # is no eager-evaluation hazard of the kind 8(a) has.
        # THE CLAMP IS TIME-DOMAIN ONLY, and that distinction is not cosmetic.
        # Under a symbolic toolkit `toolkit.where` produces
        # `Piecewise((t - td, t > td), (0.0, True))`, and the symbolic consumers of
        # this expression -- AC, transfer functions, the DDD machinery -- need a
        # SMOOTH function of `t`; a start delay has no meaning there in the first
        # place.  Two existing tests (`test_func.test_sin`, `test_elements.test_vsin`)
        # assert the smooth form, and they caught this when the clamp was applied
        # unconditionally.  They were protecting the symbolic path, not pinning the
        # defect above.
        toolkit = self.toolkit
        if getattr(toolkit, 'symbolic', False):
            dt = t - self.td
        else:
            dt = toolkit.where(t > self.td, t - self.td, 0.0)
        return self.offset + \
            self.amplitude * toolkit.exp(-self.theta * dt) * \
            toolkit.sin(self.omega * dt + self.phase)

    def dfdt(self, t):
        """Derivative of the damped sine, zero before ``td``.

        Follows `f`'s own clamping exactly: for ``t <= td`` the source is held at
        its offset, so the derivative is zero, and applying the chain rule to the
        unclamped expression there would report a slope the source does not have.
        """
        toolkit = self.toolkit
        if getattr(toolkit, 'symbolic', False):
            dt = t - self.td
            gate = 1.0
        else:
            dt = toolkit.where(t > self.td, t - self.td, 0.0)
            gate = toolkit.where(t > self.td, 1.0, 0.0)
        env = self.amplitude * toolkit.exp(-self.theta * dt)
        ph = self.omega * dt + self.phase
        return gate * env * (self.omega * toolkit.cos(ph)
                             - self.theta * toolkit.sin(ph))


class Pulse(TimeFunction):
    ## STAGE 12B -- the floor on the rise and fall times.
    ##
    ## SPICE's own default for `tr`/`tf` is zero, which makes the edge a true
    ## discontinuity: `f` jumps, and `du/dt` does not exist there at all.  That is
    ## awkward for Fang's `p = df_ckt/dh`, which needs a source derivative.
    ## Clamping the edge to the smallest step the analysis would ever take keeps
    ## the waveform a function rather than a jump, so the derivative is finite --
    ## very large, but finite and correct for the ramp actually being integrated.
    ##
    ## The value matches `Transient`'s `minstep` default.  It is a class
    ## attribute rather than a constructor argument because the sensible floor is
    ## a property of the analysis, not of the waveform, and a caller who wants a
    ## different one is usually setting it for every source at once.
    MIN_EDGE = 1e-18

    def __init__(self, v1, v2, td, tr, tf, pw, per, toolkit=numeric):
        self.v1, self.v2, self.td, self.pw, self.per = v1, v2, td, pw, per
        ## Clamped here rather than at each use, so `f`, `dfdt` and `next_event`
        ## cannot disagree about where the edges are -- which is exactly the class
        ## of defect that made `f` and `next_event` diverge before stage 8(a).
        self.tr = max(tr, self.MIN_EDGE)
        self.tf = max(tf, self.MIN_EDGE)
        self.toolkit = toolkit

    def next_event(self, t):
        """The next discontinuity STRICTLY AFTER ``t``.

        Strictly: ``next_event(t) > t`` for every ``t``.  That is the invariant
        every other ``TimeFunction`` here already keeps -- ``Sin`` returns ``td``
        only while ``t < td``, ``PWL`` does ``if pt > t``, ``Exp`` does
        ``if t < td1`` -- and ``Pulse`` was the one violator: it ended with
        ``if tmod == 0: return t``, handing back ``t`` itself at ``t = 0`` and at
        every exact period boundary.

        A caller that enumerates events by feeding the result back in then never
        advances.  ``jaxtransient.solve_batched`` does exactly that and **hung**;
        ``JAXTransient.solve`` escaped only because a second bug meant it never
        reached this method.  ``transient.py:762`` avoids both by re-calling at a
        nudged ``t`` -- a workaround for this defect, not a fix for it.

        EVENT TIMES ARE BUILT ABSOLUTELY, from the period index, rather than as
        ``t + (edge - t % per)``.  The offset form looks equivalent and is not: the
        enumeration feeds each result back in, so ``t % per`` drifts off the edge
        values, and once an edge sits less than one ULP of ``t`` above ``t % per``
        the sum rounds back to ``t`` and the loop stalls.  Measured: at
        ``t = 1.60999999999999989e-3`` the increment came to **1.08e-19** against a
        ULP of ~2.2e-19.  Anchoring to ``k*per + edge`` keeps every returned value
        on the same grid however many times it is fed back.
        """
        edges = (self.td,
                 self.td + self.tr,
                 self.td + self.tr + self.pw,
                 self.td + self.tr + self.pw + self.tf)

        ## Under a symbolic toolkit `t` is a sympy symbol and these comparisons are
        ## relationals, not bools.  There is no breakpoint schedule for symbolic
        ## time; fall through to "no further events", as `Sin` does.
        try:
            if self.per != 0:
                k = int(t // self.per)
                ## One period back as well as forward: `t` may sit a hair below a
                ## boundary that floor() has already rounded past.
                cands = [j * self.per + e
                         for j in (k - 1, k, k + 1)
                         for e in edges + (self.per,)]
            else:
                cands = list(edges)
            later = [c for c in cands if bool(c > t)]
            return min(later) if later else self.toolkit.inf
        except TypeError:
            return self.toolkit.inf

    def f(self, t):
        # --- PERIODIC PULSE WAVEFORM EVALUATION ---
        toolkit = self.toolkit
        
        # Modulo arithmetic folds continuous time t into the [0, PER) base period
        if self.per != 0:
            t = t % self.per
        
        ## STAGE 8(a) -- tr/tf OF ZERO ARE SPICE'S OWN DEFAULTS, and `where()` builds
        ## every branch EAGERLY, so `(v2 - v1)/tr` was evaluated even where the branch
        ## is not selected: a plain `VPulse()` raised ZeroDivisionError.
        ##
        ## The substituted denominator cannot leak, and that is why this is a fix
        ## rather than a mask: with `tr == 0` the rise branch's condition is
        ## `t < td + 0`, and the branch AFTER it re-selects `v1` on exactly that
        ## interval, so the ramp's value can never survive.  Same for `tf == 0`
        ## against the `v2` branch.  The edge becomes a step, which is the limit of
        ## the ramp; SPICE instead substitutes TSTEP, which is not reachable from
        ## here -- and since 9(d) the breakpoint machinery lands a step exactly on
        ## the edge either way.
        ##
        ## `per == 0` needs no such care: it is already guarded above, and the
        ## review's claim that it divides by zero does not reproduce.
        ## No zero guard needed any more: `__init__` floors both at MIN_EDGE, so
        ## these denominators are always positive.  The guard that used to live
        ## here substituted 1.0 for a zero `tr` while the branch CONDITION still
        ## used the real zero, so the rise branch was never selected -- correct,
        ## but only because two places agreed to disagree.
        tr = self.tr
        tf = self.tf

        # Phase 1: Initial Delay (td)
        v_out = self.v1
        v_out = self.toolkit.where(t < self.td + self.tr + self.pw + self.tf,
                                   self.v2 + (self.v1 - self.v2) / tf * (t - (self.td+self.tr+self.pw)),
                                   v_out)
        v_out = self.toolkit.where(t < self.td + self.tr + self.pw,
                                   self.v2,
                                   v_out)
        v_out = self.toolkit.where(t < self.td + self.tr,
                                   self.v1 + ((self.v2 - self.v1) / tr) * (t - self.td),
                                   v_out)
        v_out = self.toolkit.where(t < self.td,
                                   self.v1,
                                   v_out)
        return v_out

    def dfdt(self, t):
        """Slope of the segment being entered at ``t``: 0, rise, 0, fall, 0.

        Mirrors `f`'s branch structure exactly, including the period fold and the
        `tr`/`tf` zero guards -- a derivative that disagrees with its own function
        about which segment `t` is in is worse than no derivative at all.

        The branches take the right-hand limit at each corner, which is the
        segment a step starting there moves into -- and since all four edges are
        reported by `next_event`, a step always does start there.

        With `tr` or `tf` zero the edge is an instantaneous jump, so the rise (or
        fall) segment has zero duration and the segment actually ahead is the
        flat one.  Reporting zero is therefore the correct right-hand limit, not
        a guard against a large number.
        """
        toolkit = self.toolkit
        if self.per != 0:
            t = t % self.per

        rise = (self.v2 - self.v1) / self.tr
        fall = (self.v1 - self.v2) / self.tf

        d = 0.0 * t
        d = toolkit.where(t < self.td + self.tr + self.pw + self.tf, fall, d)
        d = toolkit.where(t < self.td + self.tr + self.pw, 0.0 * t, d)
        d = toolkit.where(t < self.td + self.tr, rise, d)
        d = toolkit.where(t < self.td, 0.0 * t, d)
        return d


class PWL(TimeFunction):
    def __init__(self, t_v_pairs, toolkit=numeric):
        """t_v_pairs is a flat list/array of alternating time and voltage/current: [t0, v0, t1, v1, ...]"""
        self.toolkit = toolkit
        if len(t_v_pairs) % 2 != 0:
            raise ValueError("PWL requires an even number of values [t0, v0, t1, v1, ...]")
        self.times = t_v_pairs[0::2]
        self.values = t_v_pairs[1::2]
        
    def next_event(self, t):
        for pt in self.times:
            if pt > t:
                return pt
        return self.toolkit.inf
        
    def f(self, t):
        if t <= self.times[0]:
            return self.values[0]
        if t >= self.times[-1]:
            return self.values[-1]
            
        # Interpolate
        for i in range(len(self.times) - 1):
            if self.times[i] <= t <= self.times[i+1]:
                t0, v0 = self.times[i], self.values[i]
                t1, v1 = self.times[i+1], self.values[i+1]
                if t1 == t0:
                    return v1
                return v0 + (v1 - v0) * (t - t0) / (t1 - t0)

    def dfdt(self, t):
        """Slope of the segment being ENTERED at ``t``.

        The right-hand limit, which is the only one the solver needs: every
        table point is reported by `next_event`, so a step lands on the corner
        and moves forward into the segment ahead.  Outside the table the source
        is held flat and the derivative is zero.
        """
        if t < self.times[0] or t >= self.times[-1]:
            return 0.0
        for i in range(len(self.times) - 1):
            t0, t1 = self.times[i], self.times[i + 1]
            if t0 <= t < t1:
                if t1 == t0:
                    return 0.0
                return (self.values[i + 1] - self.values[i]) / (t1 - t0)
        return 0.0
        return self.values[-1]

class Exp(TimeFunction):
    def __init__(self, v1, v2, td1, tau1, td2, tau2, toolkit=numeric):
        self.v1, self.v2 = v1, v2
        self.td1, self.tau1 = td1, tau1
        self.td2, self.tau2 = td2, tau2
        self.toolkit = toolkit
        
    def next_event(self, t):
        if t < self.td1: return self.td1
        if t < self.td2: return self.td2
        return self.toolkit.inf
        
    def f(self, t):
        if t <= self.td1:
            return self.v1
        elif t <= self.td2:
            return self.v1 + (self.v2 - self.v1) * (1 - self.toolkit.exp(-(t - self.td1) / self.tau1))
        else:
            v_at_td2 = self.v1 + (self.v2 - self.v1) * (1 - self.toolkit.exp(-(self.td2 - self.td1) / self.tau1))
            return v_at_td2 + (self.v1 - v_at_td2) * (1 - self.toolkit.exp(-(t - self.td2) / self.tau2))

    def dfdt(self, t):
        exp = self.toolkit.exp
        if t <= self.td1:
            return 0.0 * t
        if t <= self.td2:
            return (self.v2 - self.v1) / self.tau1 * exp(-(t - self.td1) / self.tau1)
        v_at_td2 = self.v1 + (self.v2 - self.v1) * (1 - exp(-(self.td2 - self.td1) / self.tau1))
        return (self.v1 - v_at_td2) / self.tau2 * exp(-(t - self.td2) / self.tau2)

class AM(TimeFunction):
    def __init__(self, vo, va, fc, fm, m, toolkit=numeric):
        self.vo, self.va = vo, va
        self.fc, self.fm = fc, fm
        self.m = m
        self.toolkit = toolkit
        
    def next_event(self, t):
        return self.toolkit.inf
        
    def f(self, t):
        # vo + va * (1 + m * sin(2*pi*fm*t)) * sin(2*pi*fc*t)
        pi = self.toolkit.pi
        sin = self.toolkit.sin
        mod = 1.0 + self.m * sin(2.0 * pi * self.fm * t)
        return self.vo + self.va * mod * sin(2.0 * pi * self.fc * t)

    def dfdt(self, t):
        pi = self.toolkit.pi
        sin, cos = self.toolkit.sin, self.toolkit.cos
        wm, wc = 2.0 * pi * self.fm, 2.0 * pi * self.fc
        mod = 1.0 + self.m * sin(wm * t)
        dmod = self.m * wm * cos(wm * t)
        return self.va * (dmod * sin(wc * t) + mod * wc * cos(wc * t))

class SFFM(TimeFunction):
    def __init__(self, vo, va, fc, mdi, fm, toolkit=numeric):
        self.vo, self.va = vo, va
        self.fc, self.fm = fc, fm
        self.mdi = mdi
        self.toolkit = toolkit
        
    def next_event(self, t):
        return self.toolkit.inf
        
    def f(self, t):
        # vo + va * sin(2*pi*fc*t + mdi * sin(2*pi*fm*t))
        pi = self.toolkit.pi
        sin = self.toolkit.sin
        return self.vo + self.va * sin(2.0 * pi * self.fc * t + self.mdi * sin(2.0 * pi * self.fm * t))

    def dfdt(self, t):
        pi = self.toolkit.pi
        sin, cos = self.toolkit.sin, self.toolkit.cos
        wm, wc = 2.0 * pi * self.fm, 2.0 * pi * self.fc
        phase = wc * t + self.mdi * sin(wm * t)
        return self.va * cos(phase) * (wc + self.mdi * wm * cos(wm * t))

class ScalarFunction():
    """Scalar function"""
    
    def __init__(self, toolkit=numeric):
        self.toolkit = toolkit
    
    def f(self,x):
        return 0        
    
    def fprime(self,x):
        return 0

    def F(self,x):
        return 0

class Tanh(ScalarFunction):
    """Scalar function"""
    
    def __init__(self, offset = 0, level = 0, 
                 toolkit = numeric):
        self.offset  = offset # Offset from zero
        self.level   = level  # Level of limiting
        self.toolkit = toolkit 
    
    ## A smooth symmetric limiter:  f(x) = level * tanh((x - offset) / level).
    ##
    ## The `level` factor is what makes this a limiter with a *unit* slope
    ## rather than a bare `tanh`: f'(offset) == 1, so a driver multiplying by a
    ## gain g gets small-signal gain exactly g, and saturation at g*level.
    ## Without it the slope at the origin is 1/level and every gain is off by
    ## that factor.
    ##
    ## The three methods must stay mutually consistent -- fprime is the
    ## derivative of f and F its integral -- because callers mix them:
    ## VCVS_limited stamps its Jacobian from fprime while computing its
    ## residual from f.  test_element_jacobians checks exactly that.

    # Function tanh
    def f(self,x):
        return self.level * self.toolkit.tanh((x-self.offset)/self.level)

    # Derivate
    def fprime(self,x):
        ## d/dx [level * tanh(u)] with u = (x-offset)/level  ->  1 - tanh(u)**2.
        ##
        ## This body was unreachable behind a `return 0` from 2009 until 2026,
        ## so VCVS_limited's Jacobian carried no input-to-output coupling at all
        ## and Newton had to converge without it.  The shadowing was not
        ## careless: the expression called `toolkit.power`, which no backend
        ## has ever provided, so unshadowing it alone raises AttributeError.
        ##
        ## `**` is used instead of a `power` primitive because numpy arrays,
        ## sympy expressions and JAX arrays all implement it -- adding the
        ## primitive would have meant adding it to every backend to get the
        ## same result.  Reach for an operator before a primitive.
        return 1 - self.toolkit.tanh((x-self.offset)/self.level)**2

    # Integral
    def F(self,x):
        return self.level**2 * \
            self.toolkit.log(self.toolkit.cosh((x-self.offset)/self.level))

                    
class TabFunc(ScalarFunction):
    """Return interpolated values from a lookup table

    >>> xvec=numeric.linspace(-2,2,100)
    >>> yvec=numeric.tanh(xvec)
    >>> myFunc=TabFunc(xvec,yvec)
    >>> myFunc.f(numeric.pi)
    1.02465815883
    >>> myFunc.fprime(0)
    1.00000009622
    >>> 
    
    """

    def __init__(self, xvec, yvec, s=None):
        self.xyspline=interpolate.splrep(xvec, yvec, s=s)

    def f(self,x):
        return interpolate.splev(x,self.xyspline,der=0)

    def fprime(self,x,order=1):
        return interpolate.splev(x,self.xyspline,der=order)

