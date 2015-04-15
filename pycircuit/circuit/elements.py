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
    branches = (BranchI(Node('plus'), Node('minus'), input=True, output='i', noisy=True),)
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm', 
                            default=1e3),
                  Parameter(name='noisy', desc='No noise', unit='', 
                            default=True),
                  ]

    def update(self, subject):
        self.g = 1/self.iparv.r

    def eval_iqu(self, x, epar):
        branch_v = x[0]
        return self.g * branch_v,

    def eval_noise(self, epar):
        if self.iparv.noisy:
            iPSD = 4 * self.toolkit.kboltzmann * epar.T / self.iparv.r
        else:
            iPSD = 0
        return iPSD,

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
    branches = (BranchI(Node('plus'), Node('minus'), input=True, output='q'),)
    instparams = [Parameter(name='c', desc='Capacitance', 
                            unit='F', default=1e-12)]

    def eval_iqu(self, x, epar):
        branch_v = x[0]
        q = self.iparv.c * branch_v
        return q,

class nC(Circuit):
    """Nonlinear capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(np.zeros(2))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    >>> c.C(np.zeros(2))
    array([[  1.0000e-12,  -1.0000e-12],
           [ -1.0000e-12,   1.0000e-12]])

    """

    terminals = ('plus', 'minus')
    branches = (BranchI(Node('plus'), Node('minus'), input=True, output='q'),)
    instparams = [Parameter(name='c0', desc='Capacitance', 
                            unit='F', default=1e-12),
                  Parameter(name='c1', desc='Nonlinear capacitance', 
                            unit='F', default=0.5e-12),
                  Parameter(name='v0', desc='Voltage for nominal capacitance', 
                            unit='V', default=1),
                  Parameter(name='v1', desc='Slope voltage ...', 
                            unit='V', default=1)
                  ]

    linear = False

    def eval_iqu(self, x, epar):
        branch_v = x[0]

        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1

        q = c0*branch_v+c1*v1*self.toolkit.log(self.toolkit.cosh((branch_v-v0)/v1))
        return q,

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
   """
    terminals = ('plus', 'minus')
    branches = (BranchV(Node('plus'), Node('minus'), input=True, output='q'),)

    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9)]

    def eval_iqu(self, x, epar):
        branch_i = x[0]

        q = self.iparv.L * branch_i

        return q,

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
    branches = (BranchV(Node('plus'), Node('minus'), output='u'),)
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

    def eval_iqu(self, x, epar):
        if epar.analysis == 'ac':
            phase = self.iparv.phase * self.toolkit.pi / 180
            return self.iparv.vac * self.toolkit.exp(1j*phase),
        elif epar.analysis in timedomain_analyses:
            return self.iparv.v + self.function.f(epar.t),
        else:
            return 0,

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
        Parameter(name='va', desc='Voltage amplitude', 
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
    branches = (BranchI(Node('plus'), Node('minus'), output='u'),)
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

    def eval_iqu(self, x, epar):
        if epar.analysis == 'ac':
            phase = self.iparv.phase * self.toolkit.pi / 180.            
            return self.iparv.iac * self.toolkit.exp(1j*phase),
        elif epar.analysis in timedomain_analyses:
            return self.iparv.i + self.function.f(epar.t),
        else:
            return 0,

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
    branches = ( Branch( Node('inp'),  Node('inn'),  input=True),
                 BranchV(Node('outp'), Node('outn'), output='i') )
               
    def eval_iqu(self, x, epar):
        inbranch_v = x[0]
        return self.iparv.g * inbranch_v,

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

        tk = self.toolkit
        
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

        self.den = tk.array(self.iparv.denominator)
        self.num = tk.array(self.iparv.numerator)

        # Normalize
        if self.den[0] != 1:
            self.num = self.num/self.den[0]
            self.den = self.den/self.den[0]

        self.denlen = len(self.den)
        self.numlen = len(self.num)

        # store number of nodes/states in inital G and C matrix
        self.first_state_node = len(self.nodes)

        # Add nodes for new states, one for each pole in denominator
        newnodes = ["_a%d"%state for state in range(self.denlen-1)]
        self.add_nodes(*newnodes)

        n = self.n
        G = tk.zeros((n,n), dtype = int)
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
                    tk.eye(self.denlen-2, dtype=int)
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
                    tk.eye(self.denlen-2, dtype=int)
            # Input
            G[first, inpindex] =  1
            G[first, innindex] = -1
            # Output, all numerator coefficients
            index = first + self.denlen-1 - self.numlen
            G[branchindex, index:index+self.numlen] = self.num

        self._G = G

        C = tk.zeros((n,n), dtype=int)
        C[first:first+self.denlen-1, first:first+self.denlen-1] = \
            -1*tk.eye(self.denlen-1, dtype=int)
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
    branches = ( Branch( Node('inp'),  Node('inn'),  input=True),
                 BranchI(Node('outp'), Node('outn'), output='ib') )
    instparams = [Parameter(name='gm', desc='Transconductance', 
                            unit='A/V', default=1e-3)]
    
    def eval_iqu(self, x, epar):
        inbranch_v = x[0]
        return self.iparv.gm * inbranch_v,

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
       Boston/Dordrecht/London: Kluwer Academic, �2.2.2 pp. 32-34. 
       ISBN 1402075901.

    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch('inp', 'inn', input=True), BranchV('outp', 'outn', indirect=True, output='i'))

    def eval_iqu(self, x, epar):
        input_v = x[0]
        ## The ouput branch voltage should adjust it value so that the equation
        ## invput_v == 0 is fulfilled
        eq = input_v # input_v == 0
        return eq,

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
    branches = (BranchI(Node('plus'), Node('minus'), input=True, output='i', 
                        linear=False),)
    instparams = [Parameter(name='IS', desc='Saturation current', 
                  unit='A', default=1e-13)]
    linear = False
    
    def eval_iqu(self, x, epar):
        VD = x[0]
        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        return self.iparv.IS * (self.toolkit.exp(VD/VT)-1),

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

class Idt(Circuit):
    """Integrator
    
    Output voltage is the time integral of input voltage.
    
    """
    
    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)
        
    def __init__(self, *args, **kvargs):
        super(Idt, self).__init__(*args, **kvargs)
        branchindex = -1 ## add last in matrix
        idt_index = self.nodes.index(self.add_node('idt_node')) #note side effect
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in self.terminals)
        G = self.toolkit.zeros((self.n,self.n))
        G[idt_index, inpindex] +=  1
        G[idt_index, innindex] += -1
        G[outpindex, branchindex] +=  1
        G[outnindex, branchindex] += -1
        G[branchindex, idt_index] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] +=  1
        self._G = G
        
        C = self.toolkit.zeros((self.n,self.n))
        C[idt_index, idt_index] +=  1
        self._C = C

    def C(self, x, epar=defaultepar):
        return self._C
    
    def G(self, x, epar=defaultepar):
        return self._G

class Idtmod(Circuit):
    """Modulus integrator
    
    Output voltage is the time integral of input voltage,
    modulus "modulus", and and offset.
    
    """
    instparams = [Parameter(name='modulus', desc='Output modulus',unit='V/V',
                            default=1.),
                  Parameter(name='offset', desc='offset voltage',unit='V',
                            default=0.)]
    
    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)
        
    def __init__(self, *args, **kvargs):
        super(Idtmod, self).__init__(*args, **kvargs)
        branchindex = -1 ## add last in matrix
        idt_index = self.nodes.index(self.add_node('idt_node')) #note side effect
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in self.terminals)
        G = self.toolkit.zeros((self.n,self.n))
        G[idt_index, inpindex] +=  1
        G[idt_index, innindex] += -1
        G[outpindex, branchindex] +=  1
        G[outnindex, branchindex] += -1
        G[branchindex, idt_index] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] +=  1
        self._G = G
        
        C = self.toolkit.zeros((self.n,self.n))
        C[idt_index, idt_index] +=  1
        self._C = C
        self.modulus = self.iparv.modulus
        self.offset = self.iparv.offset
        
    def C(self, x, epar=defaultepar):
        return self._C
    
    def G(self, x, epar=defaultepar):
        return self._G

    def q(self, x, epar=defaultepar): # q == -v_out in current implementation
        #self._C is constant and _q is non-zero at one index only
        _q = (self.toolkit.dot(self._C, x)  % -self.modulus)
        _qmask = np.sign(np.abs(_q))
        _q += _qmask*self.offset
        return _q

if __name__ == "__main__":
    import doctest
    doctest.testmod()
