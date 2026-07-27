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
        """Return events at peaks and zero-crossings"""
        
#        phase = self.toolkit.simplify(self.omega * (t - self.td) + self.phase)
        phase = self.omega * (t - self.td) + self.phase
        nextevent_phase = (self.toolkit.floor(phase / (self.toolkit.pi / 2)) + 1) * self.toolkit.pi / 2
        return t + (nextevent_phase - phase) / self.omega
        
    def f(self, t):
        # --- DAMPED SINE WAVE EQUATION ---
        # The SPICE sine wave is defined as a decaying/growing sinusoid:
        # V(t) = V_offset + V_amplitude * exp(-theta * (t - td)) * sin(2*pi*freq*(t - td) + phase)
        # 
        # Parameters:
        # - theta: The damping factor (1/seconds). If theta > 0, the sine wave exponentially decays.
        # - td: Time delay before the sine wave begins.
        # - omega: Angular frequency (2 * pi * freq).
        toolkit = self.toolkit
        return self.offset + \
            self.amplitude * toolkit.exp(-self.theta*(t-self.td)) * \
            toolkit.sin(self.omega * (t - self.td) + self.phase)

class Pulse(TimeFunction):
    def __init__(self, v1, v2, td, tr, tf, pw, per, toolkit=numeric):
        self.v1, self.v2, self.td, self.tr, self.tf, self.pw, self.per = \
            v1, v2, td, tr, tf, pw, per
        self.toolkit = toolkit

    def next_event(self, t):
        """
        Dynamically calculates the very next breakpoint after time t.
        Because Pulse is periodic, this avoids pre-allocating an infinite 
        number of breakpoints. The modulo arithmetic ensures that only the
        immediate next edge is returned.
        """
        if self.per != 0:
            tmod = t % self.per
        else:
            tmod = t
        
        if tmod == 0:
            return t
        elif tmod < self.td:
            return t + self.td - tmod
        elif tmod < self.td + self.tr:
            return t + self.td + self.tr - tmod
        elif tmod < self.td + self.tr + self.pw:
            return t + self.td + self.tr + self.pw - tmod
        elif tmod < self.td + self.tr + self.pw + self.tf:
            return t + self.td + self.tr + self.pw + self.tf - tmod
        else:
            return self.toolkit.ceil(t / self.per) * self.per

    def f(self, t):
        # --- PERIODIC PULSE WAVEFORM EVALUATION ---
        toolkit = self.toolkit
        
        # Modulo arithmetic folds continuous time t into the [0, PER) base period
        if self.per != 0:
            t = t % self.per
        
        # Phase 1: Initial Delay (td)
        v_out = self.v1
        v_out = self.toolkit.where(t < self.td + self.tr + self.pw + self.tf,
                                   self.v2 + (self.v1 - self.v2) / self.tf * (t - (self.td+self.tr+self.pw)),
                                   v_out)
        v_out = self.toolkit.where(t < self.td + self.tr + self.pw,
                                   self.v2,
                                   v_out)
        v_out = self.toolkit.where(t < self.td + self.tr,
                                   self.v1 + ((self.v2 - self.v1) / self.tr) * (t - self.td),
                                   v_out)
        v_out = self.toolkit.where(t < self.td,
                                   self.v1,
                                   v_out)
        return v_out

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
    
    # Function tanh
    def f(self,x):
        return self.toolkit.tanh((x-self.offset)/self.level)        
    
    # Derivate
    def fprime(self,x):
        ## d/dx tanh(u) = 1 - tanh(u)**2, with u = (x-offset)/level.
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
        return (1 - self.toolkit.tanh((x-self.offset)/self.level)**2)/self.level

    # Integral
    def F(self,x):
        return self.toolkit.log(self.toolkit.cosh((x-self.offset)/self.level))*self.level

                    
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

