import numeric

class TimeFunction(object):
    """Time dependent function"""
    
    def f(self, t):
        return 0        
        
    def next_event(self, t):
        """Returns the time of the next event given the current time t"""
        return inf
    
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
    
    def f(self, t):
        toolkit = self.toolkit
        return self.offset + \
            self.amplitude * toolkit.exp(-self.theta*(t-self.td)) * \
            toolkit.sin(self.omega * (t - self.td) + self.phase)
