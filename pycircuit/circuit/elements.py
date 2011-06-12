# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from __future__ import division
from circuit import *
import func

class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['R']
    R('plus','minus',r=1000.0,noisy=True)
    >>> c.G(self.toolkit.zeros(2))
    array([[ 0.001, -0.001],
           [-0.001,  0.001]])
    >>> c = SubCircuit()
    >>> n2=c.add_node('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(self.toolkit.zeros(2))
    array([[ 0.001, -0.001],
           [-0.001,  0.001]])

    """
    terminals = ('plus', 'minus')
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm', 
                            default=1e3),
                  Parameter(name='noisy', desc='No noise', unit='', 
                            default=True),
                  ]

    def update(self, subject):
        g = 1/self.iparv.r
        self._G = self.toolkit.array([[g, -g],
                                      [-g, g]])

    def G(self, x, epar=defaultepar): return self._G

    def CY(self, x, w, epar=defaultepar):
        if self.iparv.noisy:
            iPSD = 4 * self.toolkit.kboltzmann * epar.T / self.iparv.r
        else:
            iPSD = 0

        return  self.toolkit.array([[iPSD, -iPSD],
                                    [-iPSD, iPSD]])

class G(Circuit):
    """Conductor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['G'] = G(n1, gnd, g=1e-3)
    >>> c['G']
    G('plus','minus',g=0.001,nonoise=False)
    >>> c.G(self.toolkit.zeros(2))
    array([[ 0.001, -0.001],
           [-0.001,  0.001]])
    """
    terminals = ('plus', 'minus')
    instparams = [Parameter(name='g', desc='Conductance', unit='S', 
                            default=1e-3),
                  Parameter(name='noisy', desc='No noise', unit='', 
                            default=False)
                  ]

    def update(self, subject):
        g = self.iparv.g
        self._G = self.toolkit.array([[g, -g],
                                      [-g, g]])

    def G(self, x, epar=defaultepar): return self._G

    def CY(self, x, w, epar=defaultepar):
        if self.iparv.noisy:
            iPSD = 4*self.toolkit.kboltzmann * epar.T*self.iparv.g
            return  self.toolkit.array([[iPSD, -iPSD],
                                        [-iPSD, iPSD]])
        else:
            return super(G, self).CY(x, w, epar=epar)
        

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(self.toolkit.zeros(2))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    >>> c.C(self.toolkit.zeros(2))
    array([[  1.0000e-12,  -1.0000e-12],
           [ -1.0000e-12,   1.0000e-12]])

    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c', desc='Capacitance', 
                            unit='F', default=1e-12)]

    def update(self, subject):
        c = self.iparv.c
        self._C =  self.toolkit.array([[c, -c],
                                       [-c, c]])

    def C(self, x, epar=defaultepar): return self._C

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
    >>> c.G(self.toolkit.zeros(3))
    array([[ 0.,  0.,  1.],
           [ 0.,  0., -1.],
           [ 1., -1.,  0.]])
    >>> c.C(self.toolkit.zeros(3))
    array([[  0.0000e+00,   0.0000e+00,   0.0000e+00],
           [  0.0000e+00,   0.0000e+00,   0.0000e+00],
           [  0.0000e+00,   0.0000e+00,  -1.0000e-09]])
   """
    terminals = ('plus', 'minus')
    branches = (Branch(Node('plus'), Node('minus')),)

    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9)]

    def __init__(self, *args, **kvargs):
        super(L, self).__init__(*args, **kvargs)
        self._G = self.toolkit.array([[0 , 0, 1],
                                      [0 , 0, -1],
                                      [1 , -1, 0]])

    def update(self, subject):
        n = self.n
        C = self.toolkit.zeros((n,n))
        C[-1,-1] = -self.iparv.L
        self._C = C

    def G(self, x, epar=defaultepar): return self._G
    def C(self, x, epar=defaultepar): return self._C

class VS(Circuit):
    """Independent DC voltage source

    >>> from dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c,refnode=gnd).solve().x
    array([ 1.5   ,  0.    , -0.0015])
    
    """
    terminals = ('plus', 'minus')
    branches = (Branch(Node('plus'), Node('minus')),)
    instparams = [Parameter(name='v', desc='Source DC voltage', 
                            unit='V', default=0),
                  Parameter(name='vac', desc='AC analysis amplitude', 
                            unit='V', default=1),
                  Parameter(name='phase', desc='AC analysis phase', 
                            unit='deg', default=0),
                  Parameter(name='noisePSD', 
                            desc='Voltage noise power spectral density', 
                            unit='V^2/Hz', default=0)]
    function = func.TimeFunction()

    def update(self, subject):
        self._G = self.toolkit.array([[0 ,  0,  1],
                                      [0 ,  0, -1],
                                      [1 , -1,  0]])

    def G(self, x, epar=defaultepar): return self._G 

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'ac':
            phase = self.iparv.phase * self.toolkit.pi / 180
            vac = self.iparv.vac * self.toolkit.exp(1j*phase)
            return self.toolkit.array([0, 0, -vac])
        elif analysis in timedomain_analyses:
            v = self.iparv.v + self.function.f(t)
            return self.toolkit.array([0, 0, -v])
        else:
            return self.toolkit.array([0, 0, 0])

    def CY(self, x, w, epar=defaultepar):
        CY = super(VS, self).CY(x, w)
        CY[2, 2] = self.iparv.noisePSD
        return CY

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

class VSin(VS):
    """ Independent sinus volatge source

    """

    instparams = VS.instparams + [
        Parameter(name='vo', desc='Offset voltage', 
                  unit='V', default=0),
        Parameter(name='va', desc='Voltage voltage', 
                  unit='V', default=0),
        Parameter(name='freq', desc='Frequency', 
                  unit='Hz', default=0),
        Parameter(name='td', desc='Time delay', 
                  unit='s', default=0),
        Parameter(name='theta', desc='Damping factor', 
                  unit='1/s', default=0),
        Parameter(name='phase', desc='Phase in degrees', 
                  unit='deg', default=0)]
    def __init__(self, *args, **kvargs):
        super(VSin, self).__init__(*args, **kvargs)
        self.function = func.Sin(self.iparv.vo,
                                 self.iparv.va,
                                 self.iparv.freq,
                                 self.iparv.td,
                                 self.iparv.theta,
                                 self.iparv.phase,
                                 toolkit = self.toolkit
                                 )

class IS(Circuit):
    """Independent DC current source

    >>> from dcanalysis import DC, gnd as gnd2
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c,refnode=gnd).solve().x
    array([ 1.,  0.])

    """
    instparams = [Parameter(name='i', desc='DC Current', 
                            unit='A', default=1e-3),
                  Parameter(name='iac', desc='AC analysis current amplitude', 
                            unit='A', default=0),
                  Parameter(name='phase', desc='AC analysis phase', 
                            unit='deg', default=0),
                  Parameter(name='noisePSD', 
                            desc='Current noise power spectral density', 
                            unit='A^2/Hz', default=0.0)]
    terminals = ('plus', 'minus')
    function = func.TimeFunction()

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'ac':
            phase = self.iparv.phase * self.toolkit.pi / 180.            
            iac = self.iparv.iac * self.toolkit.exp(1j*phase)
            return self.toolkit.array([iac, -iac])
        elif analysis in timedomain_analyses:
            i = self.iparv.i + self.function.f(t)
            return self.toolkit.array([i, -i])
        else:
            return self.toolkit.array([0, 0])

    def CY(self, x, w, epar=defaultepar):
        return  self.toolkit.array([[self.iparv.noisePSD, -self.iparv.noisePSD],
                                    [-self.iparv.noisePSD, self.iparv.noisePSD]])

class ISin(IS):
    """ Independent sinus current source

    """

    instparams = IS.instparams + [
        Parameter(name='io', desc='Offset current', 
                  unit='A', default=0),
        Parameter(name='ia', desc='Current amplitude', 
                  unit='V', default=0),
        Parameter(name='freq', desc='Frequency', 
                  unit='Hz', default=0),
        Parameter(name='td', desc='Time delay', 
                  unit='s', default=0),
        Parameter(name='theta', desc='Damping factor', 
                  unit='1/s', default=0),
        Parameter(name='phase', desc='Phase in degrees', 
                  unit='deg', default=0)]
    def __init__(self, *args, **kvargs):
        super(ISin, self).__init__(*args, **kvargs)
        self.function = func.Sin(self.iparv.io,
                                 self.iparv.ia,
                                 self.iparv.freq,
                                 self.iparv.td,
                                 self.iparv.theta,
                                 self.iparv.phase,
                                 toolkit = self.toolkit
                                 )

class VPulse(VS):
    """Independent pulse voltage source

    """
    instparams = VS.instparams + [
        Parameter(name='v1', desc='Initial voltage', 
                  unit='V', default=0),
        Parameter(name='v2', desc='Pulsed value', 
                  unit='V', default=0),
        Parameter(name='td', desc='Delat time', 
                  unit='s', default=0),
        Parameter(name='tr', desc='Rise time', 
                  unit='s', default=0),
        Parameter(name='tf', desc='Fall time', 
                  unit='s', default=0),
        Parameter(name='pw', desc='Pulse width', 
                  unit='s', default=0),
        Parameter(name='per', desc='Period', 
                  unit='s', default=0)]

    def __init__(self, *args, **kvargs):
        super(VPulse, self).__init__(*args, **kvargs)
        self.function = func.Pulse(self.iparv.v1,
                                   self.iparv.v2,
                                   self.iparv.td,
                                   self.iparv.tr,
                                   self.iparv.tf,
                                   self.iparv.pw,
                                   self.iparv.per,
                                   toolkit = self.toolkit
                                 )


class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> from dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'), Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(self.toolkit.zeros(4))
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  1.],
           [ 0.,  0.,  0.,  0., -1.],
           [ 2., -2., -1.,  1.,  0.]])


    """
    instparams = [Parameter(name='g', desc='Voltage gain',unit='V/V', 
                            default=1)]

    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('outp'), Node('outn')),)
               
    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name])
             for name in self.terminals)

        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] += 1
        G[branchindex, inpindex] += self.iparv.g
        G[branchindex, innindex] += -self.iparv.g                       
        self._G = G

    def G(self, x, epar=defaultepar): return self._G


class SVCVS(Circuit):
    """Voltage controlled voltage source with frequency dependent transfer

    
    >>> from dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'), Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(self.toolkit.zeros(4))
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  0.],
           [ 0.,  0.,  0.,  0.,  1.],
           [ 0.,  0.,  0.,  0., -1.],
           [ 2., -2., -1.,  1.,  0.]])


    """
    instparams = [Parameter(name='numerator', 
                            desc='Numerator coefficients of laplace defined '
                            'transfer function',unit=None, default=(1,)),
                  Parameter(name='denominator', 
                            desc='Denominator coefficients of laplace defined '
                            'transfer function',unit=None, default=(1,0)),
                  Parameter(name='realisation', desc='State space realisation' 
                            'form for transfer function, '
                            'values \"observable\" and \"controlable\""',
                            unit=None, default='observable')]

    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('outp'), Node('outn')),)

    def __init__(self, *args, **kvargs):
        super(SVCVS, self).__init__(*args, **kvargs)

        if not(self.iparv.numerator) and not(self.iparv.denominator):
            raise Exception("Numerator and denominator not defined")
        elif not(self.iparv.numerator):
            raise Exception("Numerator not defined")
        elif not(self.iparv.denominator):
            raise Exception("Denominator not defined")

        # Ckeck that he order of the denominator is at least one less than
        # the orderof the denoiminator
        if not(len(self.iparv.numerator) < len(self.iparv.denominator)):
            raise Exception("Numerator order not less than denominator order")
        elif len(self.iparv.numerator) == 0:
            raise Exception("Numerator not defined")

        if self.iparv.denominator[0] == 0:
            raise Exception("The first coefficient in the denominator must" + 
                            " not be equal to 0")

        self.den = self.toolkit.array(self.iparv.denominator)
        self.num = self.toolkit.array(self.iparv.numerator)

        # Normalize
        if self.den[0] != 1:
            self.num = self.num/self.den[0]
            self.den = self.den/self.den[0]

        self.denlen = len(self.den)
        self.numlen = len(self.num)

        # store number of nodes/states in inital G and C matrix
        self.first_state_node = len(self.nodes)

        # Add nodes for new states, one for each pole in denominator
        newnodes = [Node("_a%d"%state) for state in range(self.denlen-1)]
        self.nodes.extend(newnodes)

        n = self.n
        G = self.toolkit.zeros((n,n), dtype = int)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name])
             for name in self.terminals)

        G[outpindex, branchindex] +=  1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] +=  1

        first = self.first_state_node

        if self.iparv.realisation == 'observable':
            # Observable canonical state space form
            # Add denominator coefficiencts
            G[first:first + self.denlen-1, first] = -self.den[1:]
            # Add states
            if self.denlen-1==1:
                G[first+1,first+1] = 1
            else:
                G[first:first+self.denlen-2, first+1:first+1+self.denlen-2] = \
                    self.toolkit.eye(self.denlen-2, dtype=int)
            # Input numerator coefficients
            index = first + self.denlen-1 - self.numlen
            G[index:index+self.numlen, inpindex] =  self.num
            G[index:index+self.numlen, innindex] = -self.num
            # Output
            G[branchindex, first] = 1
        else:
            # Controllable canonical state space form
            # Add denominator coefficiencts and states
            if self.denlen-1==1:
                G[first,first] = -self.den[1]
                G[first,first+1] = 1
            else:
                G[first, first:first + self.denlen-1 ] = -self.den[1:]
                G[first+1:first+1+self.denlen-2, first:first+self.denlen-2] = \
                    self.toolkit.eye(self.denlen-2, dtype=int)
            # Input
            G[first, inpindex] =  1
            G[first, innindex] = -1
            # Output, all numerator coefficients
            index = first + self.denlen-1 - self.numlen
            G[branchindex, index:index+self.numlen] = self.num

        self._G = G

        C = self.toolkit.zeros((n,n), dtype=int)
        C[first:first+self.denlen-1, first:first+self.denlen-1] = \
            -1*self.toolkit.eye(self.denlen-1, dtype=int)
        self._C = C

    def G(self, x, epar=defaultepar): return self._G

    def C(self, x, epar=defaultepar): return self._C

class VCCS(Circuit):
    """Voltage controlled current source

    >>> from dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1,n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vccs'] = VCCS(n1, gnd, n2, gnd, gm=1e-3)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> DC(c,refnode=gnd).solve().x
    array([ 1.5, -1.5,  0. ,  0. ])

    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    instparams = [Parameter(name='gm', desc='Transconductance', 
                            unit='A/V', default=1e-3)]
    
    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        gm=self.iparv.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, inpindex] += gm
        G[outpindex, innindex] += -gm
        G[outnindex, inpindex] += -gm
        G[outnindex, innindex] += gm
        self._G = G

    def G(self, x, epar=defaultepar): return self._G

class Nullor(Circuit):
    """Nullor

    From Wikipedia, the free encyclopedia

     A nullor is a theoretical two-port network comprised of a nullator at 
     its input and a norator at its output.[1]
     Nullors represent an ideal amplifier, having infinite current, 
     voltage, transconductance and transimpedance gain.[2] 
     Its transmission parameters are all zero.

     1. The name "nullor" was introduced by H.J. Carlin
      Singular network elements, 
      IEEE Trans. Circuit Theory, March 1965, vol. CT-11, pp. 67-72.
 
     2. Verhoeven C J M van Staveren A Monna G L E Kouwenhoven
       M H L & Yildiz E (2003). 
       Structured electronic design: negative feedback amplifiers.
       Boston/Dordrecht/London: Kluwer Academic, §2.2.2 pp. 32-34. 
       ISBN 1402075901.

    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('outp'), Node('outn')),)

    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))

        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, inpindex] += 1
        G[branchindex, innindex] += -1
        self._G = G

    def G(self, x, epar=defaultepar): return self._G

class Transformer(Circuit):
    """Ideal transformer

    >>> from dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = Transformer(n1, gnd, n2, gnd, n=2)
    >>> c['vcvs'].nodes
    [Node('inp'), Node('inn'), Node('outp'), Node('outn')]
    >>> c['vcvs'].branches
    (Branch(Node('outp'),Node('outn')),)
    >>> c['vcvs'].G(self.toolkit.zeros(4))
    array([[ 0.,  0.,  0.,  0.,  2.],
           [ 0.,  0.,  0.,  0., -2.],
           [ 0.,  0.,  0.,  0.,  1.],
           [ 0.,  0.,  0.,  0., -1.],
           [-1.,  1.,  2., -2.,  0.]])

    """
    instparams = [Parameter(name='n', desc='Winding ratio', unit='', default=1)]
    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('outp'), Node('outn')),)

    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[inpindex, branchindex] += self.iparv.n
        G[innindex, branchindex] += -self.iparv.n
        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += self.iparv.n
        G[branchindex, outnindex] += -self.iparv.n
        G[branchindex, inpindex] += -1
        G[branchindex, innindex] += 1
        self._G = G

    def G(self, x, epar=defaultepar): return self._G

class Gyrator(Circuit):
    """Gyrator

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['Gyrator'] = Gyrator(n1, gnd, gm=1)
    >>> c.G(self.toolkit.zeros(4))
    array([[  0.,   0.,  1., -1.],
           [  0.,   0., -1.,  1.],
           [ -1.,   1.,  0.,  0.],
           [  1.,  -1.,  0.,  0.]])
   """

    terminals = ('inp', 'inn', 'outp', 'outn')
    instparams = [Parameter(name='gm', desc='Transconductance', 
                            unit='A/V', default=1e-3)]
    
    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        gm=self.iparv.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        # 
        G[outpindex, inpindex] +=  gm
        G[outpindex, innindex] += -gm
        G[outnindex, inpindex] += -gm
        G[outnindex, innindex] +=  gm
        #        
        G[inpindex,  outpindex] += -gm
        G[innindex,  outpindex] +=  gm
        G[inpindex,  outnindex] +=  gm
        G[innindex,  outnindex] += -gm
        self._G = G
        
    def G(self, x, epar=defaultepar): return self._G

class Diode(Circuit):
    """ Nonlinear diode
    """
    terminals = ('plus', 'minus')
    instparams = [Parameter(name='IS', desc='Saturation current', 
                  unit='A', default=1e-13)]
    linear = False
    def G(self, x, epar=defaultepar):
        VD = x[0]-x[1]

        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        g = self.iparv.IS * self.toolkit.exp(VD/VT) / VT
        return self.toolkit.array([[g, -g],
                                   [-g, g]])

    def i(self, x, epar=defaultepar):
        """
        
        """
        VD = x[0]-x[1]
        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        I = self.iparv.IS * (self.toolkit.exp(VD/VT)-1)
        return self.toolkit.array([I, -I])

class VCVS_limited(Circuit):
    """Voltage controlled voltage source with limited output voltage.

    The output voltage is limited by a $Than$ function
    
    """
    instparams = [Parameter(name='g', desc='Voltage gain',unit='V/V',
                            default=1),
                  Parameter(name='level', desc='Limit voltage',unit='V',
                            default=0.5),
                  Parameter(name='offset', desc='offset voltage',unit='V',
                            default=0)]

    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('outp'), Node('outn')),)
    linear = False

    def __init__(self, *args, **kvargs):
        super( VCVS_limited, self).__init__(*args, **kvargs)
        self.function = func.Tanh(self.iparv.offset,
                                       self.iparv.level,
                                       toolkit = self.toolkit)                                       
    
    def G(self, x, epar=defaultepar):
        n = self.n
        G = self.toolkit.zeros((n,n))
        g_limit = self.function.fprime(x[1]-x[0])
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
        (self.nodes.index(self.nodenames[name])
        for name in self.terminals)
        G[outpindex,   branchindex] +=  1
        G[outnindex,   branchindex] += -1
        G[branchindex, outpindex]   += -1
        G[branchindex, outnindex]   +=  1
        G[branchindex, inpindex]    +=  g_limit*self.iparv.g
        G[branchindex, innindex]    += -g_limit*self.iparv.g
        return G

    def i(self, x, epar=defaultepar):
        vout = x[3] - x[2] - self.function.fprime(x[1]-x[0])*self.function.f(x[1]-x[0])
        return self.toolkit.array([0,0,x[4],-x[4],vout])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
