import numeric
import sympy

class TimeFunction(object):
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
        toolkit = self.toolkit
        
        if self.per != 0:
            t = t % self.per
        
        if t < self.td:
            return self.v1
        elif t < self.td + self.tr:
            return self.v1 + ((self.v2 - self.v1) / self.tr) * (t - self.td)
        elif t < self.td + self.tr + self.pw:
            return self.v2
        elif t < self.td + self.tr + self.pw + self.tf:
            return self.v2 + \
                (self.v1 - self.v2) / self.tf * (t - (self.td+self.tr+self.pw))
        else:
            return self.v1

        
    
        
