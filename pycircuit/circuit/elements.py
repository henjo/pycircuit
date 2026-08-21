# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from __future__ import division
from .circuit import *
from . import func as func
from ._limiting import _pnjlim

class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['R']
    R('plus','minus',r=1000.0,noisy=True)
    >>> c.G(numeric.zeros(2))
    array([[ 0.001, -0.001],
           [-0.001,  0.001]])
    >>> c = SubCircuit()
    >>> n2=c.add_node('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(numeric.zeros(2))
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

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        ## F2(b) (doc/transient_review_260820.md): the pure current equation,
        ## which is what admits R to the vmap evaluation groups -- without it,
        ## a batched per-lane override {'R': {'r': ...}} was SILENTLY IGNORED
        ## and every lane returned the same waveform (F2(a) turned that into
        ## a loud refusal; this makes the feature real).  i = g*(V+ - V-),
        ## with the stamp recovered by autodiff exactly as C's is.
        g = 1.0 / params.get('r', 1e3)
        i = g * (x[0] - x[1])
        return toolkit.array([i, -i])

    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        ## (On the JAX batched path the same stamp comes from autodiff of
        ## eval_i_pure above; this stored form serves every other toolkit.)
        return self._G

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
    G('plus','minus',g=0.001,noisy=False)
    >>> c.G(numeric.zeros(2))
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
        # --- MODIFIED NODAL ANALYSIS (MNA) STAMP ---
        # For a conductor G connected between nodes 'plus' and 'minus',
        # the current leaving the 'plus' node is I = G * (V_plus - V_minus).
        # The current entering the 'minus' node is -I.
        # This yields the classic 2x2 admittance stamp matrix in the G matrix:
        #        plus    minus
        # plus  [  G,     -G  ]
        # minus [ -G,      G  ]
        g = self.iparv.g
        self._G = self.toolkit.array([[g, -g],
                                      [-g, g]])


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

    def CY(self, x, w, epar=defaultepar):
        if self.iparv.noisy:
            iPSD = 4*self.toolkit.kboltzmann * epar.T*self.iparv.g
            return  self.toolkit.array([[iPSD, -iPSD],
                                        [-iPSD, iPSD]])
        else:
            return super().CY(x, w, epar=epar)
        

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(numeric.zeros(2))
    array([[0., 0.],
           [0., 0.]])
    >>> c.C(numeric.zeros(2))
    array([[ 1.e-12, -1.e-12],
           [-1.e-12,  1.e-12]])

    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c', desc='Capacitance', 
                            unit='F', default=1e-12),
                  ## STAGE 10.3 -- SPICE's `C ... IC=`, honoured only under
                  ## `uic=True`.  Unlike `L`'s, this is NOT an assignment: a
                  ## capacitor has no state variable of its own -- `q` is derived
                  ## from the node voltages -- so an initial voltage constrains a
                  ## DIFFERENCE of two unknowns and a set of them is solved as a
                  ## spanning tree.  See `doc/initial_conditions.md`.
                  Parameter(name='ic',
                            desc='Initial voltage for uic=True (V); None means '
                                 'unconstrained',
                            unit='V', default=None)]

    IC_KIND = 'voltage'

    def update(self, subject):
        # --- CAPACITOR MNA STAMP ---
        # A capacitor's current equation is I = C * d(V_plus - V_minus)/dt.
        # Instead of going in the conductance matrix G, it goes into the 
        # capacitance matrix C, which the Integrator will multiply by the 
        # inverse time-step (e.g. 1/dt) to convert it into an algorithmic 
        # equivalent conductance (Geq) during transient simulation.
        c = self.iparv.c
        self._C =  self.toolkit.array([[c, -c],
                                       [-c, c]])

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        # The pure charge equation for JAX
        # q = C * (V_plus - V_minus)
        c = params.get('c', 1e-12)
        v = x[0] - x[1]
        q = c * v
        return toolkit.array([q, -q])

    def C(self, x, epar=defaultepar): 
        return self.toolkit.jacobian(self.eval_q_pure, x, {'c': self.iparv.c},
                                     epar, self._C)

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
    >>> c.G(numeric.zeros(3))
    array([[ 0.,  0.,  1.],
           [ 0.,  0., -1.],
           [ 1., -1.,  0.]])
    >>> c.C(numeric.zeros(3))
    array([[ 0.e+00,  0.e+00,  0.e+00],
           [ 0.e+00,  0.e+00,  0.e+00],
           [ 0.e+00,  0.e+00, -1.e-09]])
   """
    terminals = ('plus', 'minus')
    branches = (Branch(Node('plus'), Node('minus')),)

    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9),
                  ## STAGE 10.3 -- SPICE's `L ... IC=`, honoured only under
                  ## `uic=True`, like the node `ic` on the analysis.
                  ##
                  ## An inductor's initial condition is a branch CURRENT, and its
                  ## unknown already exists in the MNA vector -- so unlike a
                  ## capacitor's initial voltage, which constrains a difference of
                  ## two node unknowns and needs a spanning-tree solve, this is a
                  ## direct assignment once the row is known.
                  Parameter(name='ic',
                            desc='Initial current for uic=True (A); None means 0',
                            unit='A', default=None)]

    ## STAGE 10.3 -- what an `ic` on this element MEANS.  Declared rather than
    ## inferred from "does it own a branch": that correlation holds for every
    ## element today and would silently mis-handle the first one that breaks it.
    IC_KIND = 'current'

    def update(self, subject):
        # --- INDUCTOR MNA STAMP (EXTRA BRANCH EQUATION) ---
        # Unlike R and C, the voltage across an inductor depends on the *derivative* 
        # of its current: V_plus - V_minus = L * di/dt.
        # Therefore, we cannot express it purely in terms of nodal voltages.
        # We must introduce an extra "branch" variable representing the current `I_L`.
        # This makes the element matrix 3x3 (V_plus, V_minus, I_L).
        n = self.n
        
        # The G matrix handles the KCL/KVL structural links:
        # Row 1 (plus node KCL): I_L flows out -> +1
        # Row 2 (minus node KCL): I_L flows in -> -1
        # Row 3 (Inductor KVL branch eqn): V_plus - V_minus ... 
        self._G = self.toolkit.array([[0 , 0, 1],
                                      [0 , 0, -1],
                                      [1 , -1, 0]])
                                      
        # The C matrix handles the dynamic part:
        # ... = L * d(I_L)/dt -> Row 3, Col 3 gets -L (moved to RHS conceptually)
        L = self.iparv.L
        self._C = self.toolkit.array([[0, 0, 0],
                                      [0, 0, 0],
                                      [0, 0, -L]])


    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        L = params.get('L', 1e-9)
        i_L = x[2]
        return toolkit.array([0.0, 0.0, -L * i_L])

    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

    def C(self, x, epar=defaultepar): 
        return self.toolkit.jacobian(self.eval_q_pure, x, {'L': self.iparv.L},
                                     epar, self._C)

class VS(Circuit):
    """Independent DC voltage source

    >>> from pycircuit.circuit.dcanalysis import DC
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


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

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

    def dudt(self, t=0.0, epar=defaultepar, analysis=None):
        """STAGE 12B -- ``du/dt`` for Fang's ``p``. The DC offset drops out."""
        if analysis in timedomain_analyses:
            return self.toolkit.array([0, 0, -self.function.dfdt(t)])
        return self.toolkit.array([0, 0, 0])

    def next_event(self, t):
        if hasattr(self, 'function') and hasattr(self.function, 'next_event'):
            return self.function.next_event(t)
        return self.toolkit.inf

    def CY(self, x, w, epar=defaultepar):
        CY = super().CY(x, w)
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
        super().__init__(*args, **kvargs)
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

    >>> from pycircuit.circuit.dcanalysis import DC, gnd as gnd2
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c,refnode=gnd).solve().x
    array([1., 0.])

    """
    instparams = [Parameter(name='i', desc='DC Current',
                            unit='A', default=0),
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

    def dudt(self, t=0.0, epar=defaultepar, analysis=None):
        """STAGE 12B -- ``du/dt`` for Fang's ``p``. The DC offset drops out."""
        if analysis in timedomain_analyses:
            di = self.function.dfdt(t)
            return self.toolkit.array([di, -di])
        return self.toolkit.array([0, 0])

    def next_event(self, t):
        if hasattr(self, 'function') and hasattr(self.function, 'next_event'):
            return self.function.next_event(t)
        return self.toolkit.inf

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
                  unit='A', default=0),
        Parameter(name='freq', desc='Frequency', 
                  unit='Hz', default=0),
        Parameter(name='td', desc='Time delay', 
                  unit='s', default=0),
        Parameter(name='theta', desc='Damping factor', 
                  unit='1/s', default=0),
        Parameter(name='phase', desc='Phase in degrees', 
                  unit='deg', default=0)]
    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.Sin(self.iparv.io,
                                 self.iparv.ia,
                                 self.iparv.freq,
                                 self.iparv.td,
                                 self.iparv.theta,
                                 self.iparv.phase,
                                 toolkit = self.toolkit
                                 )

class IPulse(IS):
    """Independent pulse current source

    """
    instparams = IS.instparams + [
        Parameter(name='i1', desc='Initial current', 
                  unit='A', default=0),
        Parameter(name='i2', desc='Pulsed value', 
                  unit='A', default=0),
        Parameter(name='td', desc='Delay time', 
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
        super().__init__(*args, **kvargs)
        self.function = func.Pulse(self.iparv.i1,
                                   self.iparv.i2,
                                   self.iparv.td,
                                   self.iparv.tr,
                                   self.iparv.tf,
                                   self.iparv.pw,
                                   self.iparv.per,
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
        Parameter(name='td', desc='Delay time', 
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
        super().__init__(*args, **kvargs)
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

    >>> from pycircuit.circuit.dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'),Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(numeric.zeros(4))
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
        branchindex = -1 ## add last in matrix
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


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G


class SVCVS(Circuit):
    """Voltage controlled voltage source with frequency dependent transfer

    
    >>> from pycircuit.circuit.dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'),Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(numeric.zeros(4))
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
        super().__init__(*args, **kvargs)

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

class CCVS(Circuit):
    """Current Controlled Voltage Source

    >>> from pycircuit.circuit.dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['is'] = IS(n1, gnd, i=1)
    >>> c['ccvs'] = CCVS(n1, gnd, n2, gnd, r=1)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'),Node('gnd', isglobal=True))]
    >>> c['ccvs'].G(numeric.zeros(6))
    array([[ 0.,  0.,  0.,  0.,  1.,  0.],
           [ 0.,  0.,  0.,  0., -1.,  0.],
           [ 0.,  0.,  0.,  0.,  0.,  1.],
           [ 0.,  0.,  0.,  0.,  0., -1.],
           [ 1., -1.,  0.,  0.,  0., -1.],
           [ 0.,  0.,  1., -1.,  0.,  0.]])


    """
    instparams = [Parameter(name='r', desc='Transresistance',unit='V/I', 
                            default=1)]

    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('inp'), Node('inn')),Branch(Node('outp'), Node('outn')))
               
    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n,n))
        branchindexK = 5 ## add last in matrix
        branchindexJ = 4 ## add last in matrix
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name])
             for name in self.terminals)

        G[outpindex, branchindexK] += 1
        G[outnindex, branchindexK] += -1
        G[branchindexK, outpindex] += 1
        G[branchindexK, outnindex] += -1

        G[inpindex, branchindexJ] += 1
        G[innindex, branchindexJ] += -1
        G[branchindexJ, inpindex] += 1
        G[branchindexJ, innindex] += -1


        G[branchindexJ, branchindexK] += -self.iparv.r
        self._G = G


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G



class VCCS(Circuit):
    """Voltage controlled current source

    >>> from pycircuit.circuit.dcanalysis import DC
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


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

class Nullor(Circuit):
    """Nullor

    From Wikipedia, the free encyclopedia

     A nullor is a theoretical two-port network comprised of a nullator at 
     its input and a norator at its output.[1]
     Nullors represent an ideal amplifier, having infinite current, 
     voltage, transconductance and transimpedance gain.[2] 
     Its transmission parameters are all zero.

     [1] The name "nullor" was introduced by H.J. Carlin
      Singular network elements,
      IEEE Trans. Circuit Theory, March 1965, vol. CT-11, pp. 67-72.

     [2] Verhoeven C J M van Staveren A Monna G L E Kouwenhoven
       M H L & Yildiz E (2003). 
       Structured electronic design: negative feedback amplifiers.
       Boston/Dordrecht/London: Kluwer Academic, 2.2.2 pp. 32-34. 
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


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

class Transformer(Circuit):
    """Ideal transformer

    >>> from pycircuit.circuit.dcanalysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = Transformer(n1, gnd, n2, gnd, n=2)
    >>> c['vcvs'].nodes
    [Node('inp'), Node('inn'), Node('outp'), Node('outn')]
    >>> c['vcvs'].branches
    [Branch(Node('outp'),Node('outn'))]
    >>> c['vcvs'].G(numeric.zeros(4))
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


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

class Gyrator(Circuit):
    """Gyrator

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> n2=c.add_node('2')
    >>> n3=c.add_node('3')
    >>> n4=c.add_node('4')
    >>> c['Gyrator'] = Gyrator(n1, n2, n3, n4, gm=1)
    >>> c['Gyrator'].G(numeric.zeros(4))
    array([[ 0.,  0.,  1., -1.],
           [ 0.,  0., -1.,  1.],
           [-1.,  1.,  0.,  0.],
           [ 1., -1.,  0.,  0.]])
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
        G[outpindex, inpindex] += -gm
        G[outpindex, innindex] +=  gm
        G[outnindex, inpindex] +=  gm
        G[outnindex, innindex] += -gm
        #        
        G[inpindex,  outpindex] +=  gm
        G[inpindex,  outnindex] += -gm
        G[innindex,  outpindex] += -gm
        G[innindex,  outnindex] +=  gm
        self._G = G
        

    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

class Diode(Circuit):
    """ Nonlinear diode
    """
    terminals = ('plus', 'minus')
    instparams = [Parameter(name='IS', desc='Saturation current', 
                  unit='A', default=1e-13)]
    linear = False

    ## STAGE 13(13-2) -- ONE `_pnjlim`, TWO CONVENTIONS, AND THAT IS THE RESULT.
    ##
    ## The declared plan was to move `Diode` onto `Semiconductor`'s state-free
    ## form: return a limited vector, let `SubCircuit.limit` write it back, and
    ## take `G` at the same point `i` is evaluated at.  That removes a real
    ## inconsistency -- on a toolkit without automatic differentiation the old
    ## `G` linearised around `_vlim` while `i` used the node voltage, so the
    ## Jacobian was not the derivative of the residual being solved.
    ##
    ## IT WAS IMPLEMENTED AND REVERTED, because it fails gate 13-5 on the first
    ## test it meets: `test_dc_pcnr_diode`, a 1 A source into a diode, stops
    ## converging.  Not diverging -- CRAWLING.  After 100 iterations the update
    ## is 3.5x over tolerance and the residual 1e4x over.  The limiter itself is
    ## fine, advancing the junction 0.136 V per iteration exactly as `pnjlim`
    ## should; what changed is that a consistent `G` at a limited point is a much
    ## smaller conductance than the residual implies, so each Newton step buys
    ## less than the limiter gives back.
    ##
    ## THE INTERESTING PART IS THAT THIS IS AN ARGUMENT *FOR* PCNR, not against
    ## the unification.  Limiting cannot be made consistent by moving a shared
    ## node voltage, because the node is not the device's to move -- everything
    ## else attached to it sees the change too.  PCNR's first key idea is exactly
    ## that each device gets its OWN unknown for the quantity it limits.  Until
    ## that exists, the stateful form stays.
    ##
    ## What is kept: both devices now call one `_pnjlim`, in `_limiting.py`.
    def reset_state(self, epar=None):
        """Forget the limiter's remembered junction voltage -- stage 10.1.

        `limit` stores `_vlim` on the instance and `G` linearises around it, so
        the state outlives the analysis that produced it: running DC twice on the
        same circuit gave **15 residual evaluations the first time and 2 the
        second**, and that difference is device state, not the circuit.
        """
        if hasattr(self, '_vlim'):
            del self._vlim

    ## STAGE 13(13-3) -- this device's PCNR unknown: one junction, plus to minus.
    pcnr_junctions = ((0, 1),)

    @staticmethod
    def pcnr_i(v, params, epar, toolkit):
        """Terminal currents as a function of the device's OWN unknown."""
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        i = params.get('IS', 1e-13) * (toolkit.exp(v / VT) - 1.0)
        return toolkit.array([i, -i])

    @staticmethod
    def pcnr_didv(v, params, epar, toolkit):
        """d(terminal currents)/dv -- the whole of this device's Jacobian
        contribution, which under PCNR lands in J_MNA/lim rather than in
        J_MNA/MNA."""
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        g = params.get('IS', 1e-13) * toolkit.exp(v / VT) / VT
        return toolkit.array([g, -g])

    def limit(self, x, x0, epar=defaultepar):
        """Stateful: remembers `_vlim` and returns None. See the note above."""
        if not hasattr(self, '_vlim'):
            self._vlim = x0[0] - x0[1]

        vnew = x[0] - x[1]
        vold = self._vlim
        try:
            bypasstol = epar.bypasstol
        except (AttributeError, KeyError):
            bypasstol = 1e-12
        try:
            if abs(vnew - vold) < bypasstol:
                self._vlim = vnew
                return
        except TypeError:
            pass

        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        ## The one shared implementation -- this is what 13-2 actually achieved.
        self._vlim = _pnjlim(vnew, vold, VT, self.iparv.IS, self.toolkit)

    def G(self, x, epar=defaultepar):
        if self.toolkit.supports('autodiff'):
            return self.toolkit.jacobian(self.eval_i_pure, x,
                                         {'IS': self.iparv.IS}, epar)
        if not hasattr(self, '_vlim'):
            self._vlim = x[0] - x[1]

        VD = self._vlim
        try:
            bypasstol = epar.bypasstol
        except (AttributeError, KeyError):
            bypasstol = 1e-12
        try:
            if hasattr(self, '_vlim_cached_G') and abs(VD - self._vlim_cached_G) < bypasstol:
                return self._G_cached
        except TypeError:
            pass

        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        g = self.iparv.IS * self.toolkit.exp(VD / VT) / VT

        self._G_cached = self.toolkit.array([[g, -g], [-g, g]])
        self._vlim_cached_G = VD
        return self._G_cached

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v = x[0] - x[1]
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        IS = params.get('IS', 1e-13)
        i = IS * (toolkit.exp(v / VT) - 1.0)
        return toolkit.array([i, -i])

    def i(self, x, epar=defaultepar):
        if self.toolkit.supports('autodiff'):
            return self.eval_i_pure(x, {'IS': self.iparv.IS}, epar, self.toolkit)
        if not hasattr(self, '_vlim'):
            self._vlim = x[0]-x[1]

        VD = self._vlim
        try:
            bypasstol = epar.bypasstol
        except (AttributeError, KeyError):
            bypasstol = 1e-12
        bypass = False
        try:
            bypass = hasattr(self, '_vlim_cached_i') and abs(VD - self._vlim_cached_i) < bypasstol
        except TypeError:
            pass
            
        if bypass:
            I, g = self._I_cached, self._g_cached
        else:
            VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
            I = self.iparv.IS * (self.toolkit.exp(VD/VT)-1)
            g = self.iparv.IS * self.toolkit.exp(VD/VT) / VT
            self._I_cached = I
            self._g_cached = g
            self._vlim_cached_i = VD
            
        v_nodes = x[0] - x[1]
        I_eff = I - g * (VD - v_nodes)
        
        return self.toolkit.array([I_eff, -I_eff])

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
        super().__init__(*args, **kvargs)
        ## A pure limiter, saturating at +/-level.  The offset is applied by
        ## this element to its *input*, not by the limiter, because `offset` is
        ## an input-referred offset voltage while `level` clamps the output.
        self.function = func.Tanh(0.0, self.iparv.level,
                                  toolkit = self.toolkit)
    
    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v_inp = x[0]
        v_inn = x[1]
        v_outp = x[2]
        v_outn = x[3]
        i_branch = x[4]
        g = params.get('g', 1.0)
        level = params.get('level', 0.5)
        offset = params.get('offset', 0.0)
        
        ## func.Tanh.f inlined, so the expression stays traceable by an autodiff
        ## toolkit.  It must track func.Tanh *exactly*: G()'s hand-written path
        ## differentiates that object, and the two forms are only
        ## interchangeable while they agree.  An earlier inlining here used
        ## ``offset + level*tanh(...)``, which func.Tanh does not, so the
        ## differentiated and stamped paths were limiting differently.
        ## vout = level * tanh(g * vin / level):  the gain sits *inside* the
        ## limiter, so the output saturates at +/-level whatever the gain is,
        ## while the small-signal gain is still exactly g.  With the gain
        ## outside it instead, the clamp would scale with g.
        din = v_inp - v_inn - offset
        f = level * toolkit.tanh(g * din / level)

        ## Branch equation  v_outn - v_outp - g*f(v_inn - v_inp) = 0.
        ##
        ## This is fixed by the hand-written G(), which stamps -1/+1 on the
        ## output nodes and ±g*f'(u) on the inputs; the residual with those
        ## derivatives is the one above, and it agrees with the stamp to
        ## finite-difference precision at every operating point tested.
        ##
        ## It previously read ``- fprime(u)*f(u)`` with no g at all, which
        ## matched neither the stamp nor the documented behaviour -- for
        ## v_inp - v_inn beyond the knee it even carried the opposite sign.
        vout = v_outn - v_outp + f
        return toolkit.array([0.0, 0.0, i_branch, -i_branch, vout])

    def i(self, x, epar=defaultepar):
        params = {'g': self.iparv.g, 'level': self.iparv.level, 'offset': self.iparv.offset}
        return self.eval_i_pure(x, params, epar, self.toolkit)

    def G(self, x, epar=defaultepar):
        if self.toolkit.supports('autodiff'):
            return self.toolkit.jacobian(self.eval_i_pure, x,
                                         {'g': self.iparv.g, 'level': self.iparv.level, 'offset': self.iparv.offset}, epar)
        n = self.n
        G = self.toolkit.zeros((n,n))
        ## f' is evaluated at the *amplified, offset-shifted* input, since that
        ## is what the limiter sees -- d/dv_inp [level*tanh(g*din/level)] is
        ## g*f'(g*din).  Evaluating it at the raw input instead would be right
        ## only in the small-signal limit.
        din = x[0] - x[1] - self.iparv.offset
        g_limit = self.function.fprime(self.iparv.g * din)
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

class Idt(Circuit):
    """Integrator
    
    Output voltage is the time integral of input voltage.
    
    """
    
    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)
        
    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
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
    modulo "modulus", and an offset.

    >>> import pycircuit.circuit._numeric as numeric
    >>> from pycircuit.circuit.transient import Transient
    >>> c = SubCircuit()
    >>> nin, nout = c.add_nodes('in', 'out')
    >>> c['vin'] = VS(nin, gnd, v=1.0)
    >>> c['R'] = R(nout, gnd, r=1e3)
    >>> c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)

    ``uic=True`` is required, not incidental: an ideal integrator has no DC
    operating point -- its output is the unbounded integral of a constant input --
    so the bias solve is singular and the run must be told to start from zeros.

    ``fixed_timestep=True`` is also deliberate: this example is documenting the
    modulo wrap, so it wants output on the uniform grid it asks for. An adaptive
    run opens at ``timestep*1e-3`` and grows from there (see ``firststep``), which
    is right for accuracy and wrong for a doctest about output values.

    >>> tran = Transient(c, toolkit=numeric, uic=True)
    >>> result = tran.solve(tend=1.5, timestep=0.5, fixed_timestep=True)
    >>> result.v(nout).y
    array([0. , 0.5, 0. , 0.5])
    """
    instparams = [Parameter(name='modulus', desc='Output modulus',unit='V/V',
                            default=1.),
                  Parameter(name='offset', desc='offset voltage',unit='V',
                            default=0.)]
    
    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)
    linear = False
        
    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        branchindex = -1 ## add last in matrix
        self._idt_index = self.nodes.index(self.add_node('idt_node')) #note side effect
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in self.terminals)
        G = self.toolkit.zeros((self.n,self.n))
        G[self._idt_index, inpindex] +=  1
        G[self._idt_index, innindex] += -1
        G[outpindex, branchindex] +=  1
        G[outnindex, branchindex] += -1
        G[branchindex, self._idt_index] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] +=  1
        self._G = G
        
        C = self.toolkit.zeros((self.n,self.n))
        C[self._idt_index, self._idt_index] +=  1
        self._C = C
        self.modulus = self.iparv.modulus
        self.offset = self.iparv.offset
        
    def C(self, x, epar=defaultepar):
        return self._C
    
    def G(self, x, epar=defaultepar):
        return self._G

    def i(self, x, epar=defaultepar):
        _i = self.toolkit.dot(self._G, x)
        branchindex = -1
        
        # Remove the linear term for v_idt
        _i[branchindex] -= self._G[branchindex, self._idt_index] * x[self._idt_index]
        
        # Add the nonlinear modulo term
        v_mod = (-x[self._idt_index] % self.modulus) + self.offset
        _i[branchindex] += v_mod
        
        return _i

class VPWL(VS):
    """Independent piecewise linear voltage source"""
    instparams = VS.instparams + [
        Parameter(name='tvpairs', desc='List of time/voltage pairs', default=[])]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.PWL(self.iparv.tvpairs, toolkit=self.toolkit)

class IPWL(IS):
    """Independent piecewise linear current source"""
    instparams = IS.instparams + [
        Parameter(name='tvpairs', desc='List of time/current pairs', default=[])]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.PWL(self.iparv.tvpairs, toolkit=self.toolkit)

class VExp(VS):
    """Independent exponential voltage source"""
    instparams = VS.instparams + [
        Parameter(name='v1', desc='Initial voltage', unit='V', default=0),
        Parameter(name='v2', desc='Pulsed voltage', unit='V', default=0),
        Parameter(name='td1', desc='Rise delay time', unit='s', default=0),
        Parameter(name='tau1', desc='Rise time constant', unit='s', default=0),
        Parameter(name='td2', desc='Fall delay time', unit='s', default=0),
        Parameter(name='tau2', desc='Fall time constant', unit='s', default=0)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.Exp(self.iparv.v1, self.iparv.v2, self.iparv.td1, 
                                 self.iparv.tau1, self.iparv.td2, self.iparv.tau2, 
                                 toolkit=self.toolkit)

class IExp(IS):
    """Independent exponential current source"""
    instparams = IS.instparams + [
        Parameter(name='i1', desc='Initial current', unit='A', default=0),
        Parameter(name='i2', desc='Pulsed current', unit='A', default=0),
        Parameter(name='td1', desc='Rise delay time', unit='s', default=0),
        Parameter(name='tau1', desc='Rise time constant', unit='s', default=0),
        Parameter(name='td2', desc='Fall delay time', unit='s', default=0),
        Parameter(name='tau2', desc='Fall time constant', unit='s', default=0)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.Exp(self.iparv.i1, self.iparv.i2, self.iparv.td1, 
                                 self.iparv.tau1, self.iparv.td2, self.iparv.tau2, 
                                 toolkit=self.toolkit)

class VAM(VS):
    """Independent Amplitude Modulated voltage source"""
    instparams = VS.instparams + [
        Parameter(name='vo', desc='Offset voltage', unit='V', default=0),
        Parameter(name='va', desc='Amplitude', unit='V', default=1),
        Parameter(name='fc', desc='Carrier frequency', unit='Hz', default=1e3),
        Parameter(name='fm', desc='Modulation frequency', unit='Hz', default=1e2),
        Parameter(name='m', desc='Modulation index', unit='', default=1.0)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.AM(self.iparv.vo, self.iparv.va, self.iparv.fc, 
                                self.iparv.fm, self.iparv.m, toolkit=self.toolkit)

class IAM(IS):
    """Independent Amplitude Modulated current source"""
    instparams = IS.instparams + [
        Parameter(name='io', desc='Offset current', unit='A', default=0),
        Parameter(name='ia', desc='Amplitude', unit='A', default=1),
        Parameter(name='fc', desc='Carrier frequency', unit='Hz', default=1e3),
        Parameter(name='fm', desc='Modulation frequency', unit='Hz', default=1e2),
        Parameter(name='m', desc='Modulation index', unit='', default=1.0)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.AM(self.iparv.io, self.iparv.ia, self.iparv.fc, 
                                self.iparv.fm, self.iparv.m, toolkit=self.toolkit)

class VSFFM(VS):
    """Independent Single Frequency FM voltage source"""
    instparams = VS.instparams + [
        Parameter(name='vo', desc='Offset voltage', unit='V', default=0),
        Parameter(name='va', desc='Amplitude', unit='V', default=1),
        Parameter(name='fc', desc='Carrier frequency', unit='Hz', default=1e3),
        Parameter(name='mdi', desc='Modulation index', unit='', default=1.0),
        Parameter(name='fm', desc='Modulation frequency', unit='Hz', default=1e2)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.SFFM(self.iparv.vo, self.iparv.va, self.iparv.fc, 
                                  self.iparv.mdi, self.iparv.fm, toolkit=self.toolkit)

class ISFFM(IS):
    """Independent Single Frequency FM current source"""
    instparams = IS.instparams + [
        Parameter(name='io', desc='Offset current', unit='A', default=0),
        Parameter(name='ia', desc='Amplitude', unit='A', default=1),
        Parameter(name='fc', desc='Carrier frequency', unit='Hz', default=1e3),
        Parameter(name='mdi', desc='Modulation index', unit='', default=1.0),
        Parameter(name='fm', desc='Modulation frequency', unit='Hz', default=1e2)]

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self.function = func.SFFM(self.iparv.io, self.iparv.ia, self.iparv.fc, 
                                  self.iparv.mdi, self.iparv.fm, toolkit=self.toolkit)

class CoupledInductors(Circuit):
    """Coupled Inductors (Mutual Inductance)"""
    terminals = ('p1', 'm1', 'p2', 'm2')
    branches = (Branch(Node('p1'), Node('m1')), Branch(Node('p2'), Node('m2')))

    instparams = [Parameter(name='L1', desc='Inductance 1', unit='H', default=1e-9),
                  Parameter(name='L2', desc='Inductance 2', unit='H', default=1e-9),
                  Parameter(name='K', desc='Coupling Coefficient', unit='', default=0.99)]

    def update(self, subject):
        n = self.n
        # Nodes: p1, m1, p2, m2, br1, br2
        # G matrix inserts 1, -1 for branch currents into node KCLs, and 1, -1 for node voltages into branch equations
        G = self.toolkit.zeros((n, n))
        # Branch 1 KCL contributions
        G[0, 4] = 1; G[1, 4] = -1
        # Branch 1 voltage equation
        G[4, 0] = 1; G[4, 1] = -1
        
        # Branch 2 KCL contributions
        G[2, 5] = 1; G[3, 5] = -1
        # Branch 2 voltage equation
        G[5, 2] = 1; G[5, 3] = -1
        
        self._G = G
        
        C = self.toolkit.zeros((n, n))
        L1, L2, K = self.iparv.L1, self.iparv.L2, self.iparv.K
        M = K * self.toolkit.sqrt(L1 * L2)
        
        C[4, 4] = -L1
        C[5, 5] = -L2
        C[4, 5] = -M
        C[5, 4] = -M
        self._C = C


    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        v_p1, v_m1, v_p2, v_m2, i1, i2 = x[0], x[1], x[2], x[3], x[4], x[5]
        L1 = params.get('L1', 1e-9)
        L2 = params.get('L2', 1e-9)
        K = params.get('K', 0.99)
        M = K * toolkit.sqrt(L1 * L2)
        return toolkit.array([
            0.0, 0.0, 0.0, 0.0,
            -L1 * i1 - M * i2,
            -M * i1 - L2 * i2
        ])

    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

    def C(self, x, epar=defaultepar): 
        return self.toolkit.jacobian(
            self.eval_q_pure, x, {'L1': self.iparv.L1, 'L2': self.iparv.L2, 'K': self.iparv.K}, epar, self._C)

class VSwitch(Circuit):
    """Voltage Controlled Switch"""
    terminals = ('plus', 'minus', 'cp', 'cm')
    instparams = [Parameter('Ron', 'On resistance', default=1.0),
                  Parameter('Roff', 'Off resistance', default=1e6),
                  Parameter('Von', 'Control voltage for on state', default=1.0),
                  Parameter('Voff', 'Control voltage for off state', default=0.0)]

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v = x[0] - x[1]
        vc = x[2] - x[3]
        Ron = params.get('Ron', 1.0)
        Roff = params.get('Roff', 1e6)
        Von = params.get('Von', 1.0)
        Voff = params.get('Voff', 0.0)
        
        Gon = 1.0 / Ron
        Goff = 1.0 / Roff
        
        Vmid = (Von + Voff) / 2.0
        Vscale = Von - Voff
        if Vscale == 0:
            Vscale = 1e-6
            
        x_norm = (vc - Vmid) / Vscale
        factor = (toolkit.tanh(x_norm * 2.0) + 1.0) / 2.0
        g = Goff + (Gon - Goff) * factor
        
        i_val = v * g
        return toolkit.array([i_val, -i_val, 0.0, 0.0])

    def i(self, x, epar=defaultepar):
        params = {'Ron': self.iparv.Ron, 'Roff': self.iparv.Roff, 'Von': self.iparv.Von, 'Voff': self.iparv.Voff}
        return self.eval_i_pure(x, params, epar, self.toolkit)

    def G(self, x, epar=defaultepar):
        if self.toolkit.supports('autodiff'):
            return self.toolkit.jacobian(self.eval_i_pure, x,
                                         {'Ron': self.iparv.Ron, 'Roff': self.iparv.Roff, 'Von': self.iparv.Von, 'Voff': self.iparv.Voff}, epar)
        v = x[0] - x[1]
        vc = x[2] - x[3]
        Ron = self.iparv.Ron
        Roff = self.iparv.Roff
        Von = self.iparv.Von
        Voff = self.iparv.Voff
        
        Gon = 1.0 / Ron
        Goff = 1.0 / Roff
        
        Vmid = (Von + Voff) / 2.0
        Vscale = Von - Voff
        if Vscale == 0:
            Vscale = 1e-6
            
        x_norm = (vc - Vmid) / Vscale
        factor = (self.toolkit.tanh(x_norm * 2.0) + 1.0) / 2.0
        g = Goff + (Gon - Goff) * factor
        
        dfactor = 1.0 / self.toolkit.cosh(x_norm * 2.0)**2 * 2.0 / Vscale / 2.0
        dg = (Gon - Goff) * dfactor
        
        g_vc = v * dg
        
        return self.toolkit.array([[g, -g, g_vc, -g_vc],
                                   [-g, g, -g_vc, g_vc],
                                   [0.0, 0.0, 0.0, 0.0],
                                   [0.0, 0.0, 0.0, 0.0]])

class CCCS(Circuit):
    """Current Controlled Current Source"""
    instparams = [Parameter(name='F', desc='Current gain', unit='A/A', default=1.0)]
    terminals = ('inp', 'inn', 'outp', 'outn')
    branches = (Branch(Node('inp'), Node('inn')),)
    
    def update(self, subject):
        n = self.n
        G = self.toolkit.zeros((n, n))
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in self.terminals)
        branchindex = 4
        
        # Branch KCL (Input acts as short, measuring current)
        G[inpindex, branchindex] += 1
        G[innindex, branchindex] += -1
        
        # Branch equation: V_inp - V_inn = 0
        G[branchindex, inpindex] += 1
        G[branchindex, innindex] += -1
        
        # Output current injections
        G[outpindex, branchindex] += self.iparv.F
        G[outnindex, branchindex] += -self.iparv.F
        
        self._G = G


    def G(self, x, epar=defaultepar):
        ## Linear element: the Jacobian *is* this stored matrix, so there
        ## is nothing to differentiate and no second form to keep in step
        ## with it.  The base class derives i(x) as dot(G(x), x).
        return self._G

class ISwitch(Circuit):
    """Current Controlled Switch"""
    terminals = ('plus', 'minus', 'cp', 'cm')
    branches = (Branch(Node('cp'), Node('cm')),)
    instparams = [Parameter('Ron', 'On resistance', default=1.0),
                  Parameter('Roff', 'Off resistance', default=1e6),
                  Parameter('Ion', 'Control current for on state', default=1e-3),
                  Parameter('Ioff', 'Control current for off state', default=0.0)]
                  
    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v = x[0] - x[1]
        i_ctrl = x[4]
        Ron = params.get('Ron', 1.0)
        Roff = params.get('Roff', 1e6)
        Ion = params.get('Ion', 1e-3)
        Ioff = params.get('Ioff', 0.0)
        
        Gon, Goff = 1.0/Ron, 1.0/Roff
        Imid = (Ion + Ioff) / 2.0
        Iscale = Ion - Ioff
        if Iscale == 0: Iscale = 1e-9
        
        x_norm = (i_ctrl - Imid) / Iscale
        factor = (toolkit.tanh(x_norm * 2.0) + 1.0) / 2.0
        g = Goff + (Gon - Goff) * factor
        i_val = v * g
        return toolkit.array([i_val, -i_val, i_ctrl, -i_ctrl, 0.0])

    def i(self, x, epar=defaultepar):
        params = {'Ron': self.iparv.Ron, 'Roff': self.iparv.Roff, 'Ion': self.iparv.Ion, 'Ioff': self.iparv.Ioff}
        return self.eval_i_pure(x, params, epar, self.toolkit)

    def G(self, x, epar=defaultepar):
        if self.toolkit.supports('autodiff'):
            return self.toolkit.jacobian(self.eval_i_pure, x,
                                         {'Ron': self.iparv.Ron, 'Roff': self.iparv.Roff, 'Ion': self.iparv.Ion, 'Ioff': self.iparv.Ioff}, epar)
        v = x[0] - x[1]
        i_ctrl = x[4]
        Ron, Roff, Ion, Ioff = self.iparv.Ron, self.iparv.Roff, self.iparv.Ion, self.iparv.Ioff
        Gon, Goff = 1.0/Ron, 1.0/Roff
        Imid = (Ion + Ioff) / 2.0
        Iscale = Ion - Ioff
        if Iscale == 0: Iscale = 1e-9
        x_norm = (i_ctrl - Imid) / Iscale
        
        tanh_val = self.toolkit.tanh(x_norm * 2.0)
        factor = (tanh_val + 1.0) / 2.0
        g = Goff + (Gon - Goff) * factor
        d_factor = (1.0 - tanh_val**2) / Iscale
        dg_di = (Gon - Goff) * d_factor
        
        G_mat = self.toolkit.zeros((5, 5))
        
        # d(i_val)/dv = g
        G_mat[0, 0] = g; G_mat[0, 1] = -g
        G_mat[1, 0] = -g; G_mat[1, 1] = g
        
        # d(i_val)/di_ctrl = v * dg_di
        g_i = v * dg_di
        G_mat[0, 4] = g_i
        G_mat[1, 4] = -g_i
        
        # cp/cm KCL contributions
        G_mat[2, 4] = 1.0
        G_mat[3, 4] = -1.0
        
        # Branch equation V_cp - V_cm = 0
        G_mat[4, 2] = 1.0
        G_mat[4, 3] = -1.0
        
        return G_mat

class BSource(Circuit):
    """Behavioral Source for non-linear current and charge (capacitance)
    
    Can evaluate both a static current function i_out = i_func(v_ctrl)
    and a charge function q_out = q_func(v_ctrl).
    """
    terminals = ('inp', 'inn', 'outp', 'outn')
    instparams = [Parameter('i_func', 'Function i_out = f(v_ctrl)', default=None),
                  Parameter('q_func', 'Function q_out = f(v_ctrl)', default=None)]
    
    def i(self, x, epar=defaultepar):
        if self.iparv.i_func is None:
            return self.toolkit.zeros(4)
        v_ctrl = x[0] - x[1]
        i_val = self.iparv.i_func(v_ctrl)
        return self.toolkit.array([0.0, 0.0, i_val, -i_val])

    def G(self, x, epar=defaultepar):
        if self.iparv.i_func is None:
            return self.toolkit.zeros((4, 4))
        v_ctrl = x[0] - x[1]
        
        di_dv = self.toolkit.derivative(self.iparv.i_func, v_ctrl)

        G_mat = self.toolkit.zeros((4, 4))
        G_mat[2, 0] = di_dv
        G_mat[2, 1] = -di_dv
        G_mat[3, 0] = -di_dv
        G_mat[3, 1] = di_dv
        return G_mat

    def q(self, x, epar=defaultepar):
        if self.iparv.q_func is None:
            return self.toolkit.zeros(4)
        v_ctrl = x[0] - x[1]
        q_val = self.iparv.q_func(v_ctrl)
        return self.toolkit.array([0.0, 0.0, q_val, -q_val])

    def C(self, x, epar=defaultepar):
        if self.iparv.q_func is None:
            return self.toolkit.zeros((4, 4))
        v_ctrl = x[0] - x[1]
        
        dq_dv = self.toolkit.derivative(self.iparv.q_func, v_ctrl)

        C_mat = self.toolkit.zeros((4, 4))
        C_mat[2, 0] = dq_dv
        C_mat[2, 1] = -dq_dv
        C_mat[3, 0] = -dq_dv
        C_mat[3, 1] = dq_dv
        return C_mat

NonLinearVCCS = BSource

if __name__ == "__main__":
    import doctest
    doctest.testmod()

class TLine(Circuit):
    """
    Ideal Lossless Transmission Line using Method of Characteristics.
    
    JAX VECTORIZATION NOTE:
    TLine intentionally does NOT implement `eval_i_pure` and is completely 
    excluded from JAX `vmap` vectorization mapping. 
    
    Why?
    Unlike all other elements, TLine depends on *historical* states (t - delay).
    To vectorize this in JAX, we would need to pre-allocate massive fixed-shape
    2D history buffers and pass them into the `vmap` on every single iteration,
    then perform JAX-native temporal interpolation. 
    Doing so would dramatically slow down the evaluation of all other simple 
    elements (resistors, diodes) by exhausting GPU memory bandwidth just to carry 
    around unused historical time-series data. It would also break Symbolic evaluations.
    
    By leaving TLine on the legacy Python fallback path, the rest of the circuit 
    remains blazingly fast in pure JAX, while Python handles the messy dynamic 
    memory lookups for the delayed states!
    """
    terminals = ('p1', 'm1', 'p2', 'm2')
    branches = (Branch(Node('p1'), Node('m1')), Branch(Node('p2'), Node('m2')))
    
    instparams = [
        Parameter('Z0', 'Characteristic Impedance', default=50.0, unit='Ohm'),
        Parameter('TD', 'Time Delay', default=1e-9, unit='s')
    ]
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.history = []  # List of (t, V1, V2, I1, I2) tuples
        
    def reset_state(self, epar=None):
        """Drop the history from any previous analysis -- stage 8(d).

        `G()` selects the DC stamp on `len(self.history) == 0`, so leftover history
        made a DC solve return a transient-flavoured answer: **v(b) = 0.500000
        before** a transient and **0.000000 after**, on a matched line where 0.5 is
        correct.  Resetting per analysis makes the stamp depend on the analysis
        being run rather than on what happened to run before it.
        """
        self.history = []

    def max_timestep(self):
        """``TD/2`` -- the line cannot be resolved by steps as long as its delay.

        Measured under `fixed_timestep` with TD = 1e-9, the observed delay was
        2.00x TD at dt = 1e-9, 4.00x at 2e-9 and 8.00x at 5e-9, with no warning.
        Half the delay is the coarsest sampling that leaves the interpolation two
        points to work with; the plan asks for exactly this bound.
        """
        return float(self.iparv.TD) / 2.0

    def accept_step(self, t, x, epar):
        # Extract voltages and currents from x
        # x is [Vp1, Vm1, Vp2, Vm2, I1, I2]
        v1 = x[0] - x[1]
        v2 = x[2] - x[3]
        i1 = x[4]
        i2 = x[5]
        self.history.append((float(t), float(v1), float(v2), float(i1), float(i2)))
        
        # Prune history older than t - TD (keep at least one point before t-TD for interpolation)
        td = self.iparv.TD
        while len(self.history) > 2 and self.history[1][0] < t - td:
            self.history.pop(0)

    def _interpolate_history(self, t_target):
        # --- TRANSMISSION LINE HISTORY INTERPOLATION ---
        # A lossless transmission line uses the Method of Characteristics.
        # This requires looking up the exact voltage/current states from exactly 
        # `t - TD` (Time Delay) ago. Because the adaptive timestep `dt` rarely 
        # lands exactly on a multiple of TD, we must mathematically interpolate 
        # the historical points.
        #
        # We use a 3-point Quadratic Lagrange Interpolation (rather than linear) 
        # because the Integrator (Gear2/Trapezoidal) fits 2nd order polynomials.
        # If the transmission line only used linear interpolation, it would inject 
        # 1st-order truncation errors back into a 2nd-order solver, destroying 
        # the convergence and creating numerical ringing.
        if len(self.history) == 0:
            return 0.0, 0.0, 0.0, 0.0
            
        if len(self.history) == 1 or t_target <= self.history[0][0]:
            h = self.history[0]
            return h[1], h[2], h[3], h[4]
            
        if t_target >= self.history[-1][0]:
            h = self.history[-1]
            return h[1], h[2], h[3], h[4]
            
        # Find the interval
        for i in range(len(self.history) - 1):
            if self.history[i][0] <= t_target <= self.history[i+1][0]:
                # Exact match avoids interpolation (perfect for coupled solver)
                if abs(t_target - self.history[i][0]) < 1e-15:
                    h = self.history[i]
                    return h[1], h[2], h[3], h[4]
                if abs(t_target - self.history[i+1][0]) < 1e-15:
                    h = self.history[i+1]
                    return h[1], h[2], h[3], h[4]
                
                # Use quadratic interpolation if we have a point before the interval
                if i > 0:
                    t0, v1_0, v2_0, i1_0, i2_0 = self.history[i-1]
                    t1, v1_1, v2_1, i1_1, i2_1 = self.history[i]
                    t2, v1_2, v2_2, i1_2, i2_2 = self.history[i+1]
                    
                    # Lagrange basis polynomials
                    L0 = (t_target - t1) * (t_target - t2) / ((t0 - t1) * (t0 - t2))
                    L1 = (t_target - t0) * (t_target - t2) / ((t1 - t0) * (t1 - t2))
                    L2 = (t_target - t0) * (t_target - t1) / ((t2 - t0) * (t2 - t1))
                    
                    v1 = v1_0 * L0 + v1_1 * L1 + v1_2 * L2
                    v2 = v2_0 * L0 + v2_1 * L1 + v2_2 * L2
                    i1 = i1_0 * L0 + i1_1 * L1 + i1_2 * L2
                    i2 = i2_0 * L0 + i2_1 * L1 + i2_2 * L2
                    return v1, v2, i1, i2
                else:
                    # Linear interpolation (fallback)
                    t0, v1_0, v2_0, i1_0, i2_0 = self.history[i]
                    t1, v1_1, v2_1, i1_1, i2_1 = self.history[i+1]
                    alpha = (t_target - t0) / (t1 - t0)
                    v1 = v1_0 + alpha * (v1_1 - v1_0)
                    v2 = v2_0 + alpha * (v2_1 - v2_0)
                    i1 = i1_0 + alpha * (i1_1 - i1_0)
                    i2 = i2_0 + alpha * (i2_1 - i2_0)
                    return v1, v2, i1, i2
                
        h = self.history[-1]
        return h[1], h[2], h[3], h[4]

    def G(self, x, epar=None):
        Z0 = self.iparv.Z0
        
        # Build G_mat in standard python list to avoid JAX immutability issues
        G = [[0.0 for _ in range(6)] for _ in range(6)]
        
        # Check if we are in DC (no history)
        if len(self.history) == 0:
            G[4][0] = 1.0
            G[4][1] = -1.0
            G[4][2] = -1.0
            G[4][3] = 1.0
            G[5][4] = 1.0
            G[5][5] = 1.0
        else:
            G[4][0] = 1.0
            G[4][1] = -1.0
            G[4][4] = -Z0
            G[5][2] = 1.0
            G[5][3] = -1.0
            G[5][5] = -Z0
            
        G[0][4] = 1.0
        G[1][4] = -1.0
        G[2][5] = 1.0
        G[3][5] = -1.0
        
        return self.toolkit.array(G)

        
    def dudt(self, t=0.0, epar=None, analysis=None):
        """STAGE 12B -- deliberately NOT implemented.

        `TLine`'s source vector is built from the stored history interpolated at
        ``t - TD``, so ``du/dt`` is the derivative of that interpolation and is a
        real term, not a zero.  Inheriting `Circuit.dudt`'s zero would silently
        drop it from ``p`` and leave the coupled solve working on a subtly wrong
        problem -- the failure mode this whole stage has been chasing.

        Raising instead means a coupled run on a circuit containing a delay line
        stops and says why.  The standard adaptive path is unaffected.
        """
        if analysis not in ('tran', 'transient', 'time'):
            return self.toolkit.zeros(6)
        raise NotImplementedError(
            'TLine does not provide du/dt: its source vector comes from the '
            'delayed history, so the derivative is that of the history '
            'interpolation. The coupled (Fang) time-stepping method cannot be '
            'used on a circuit containing a delay line until this is written.')

    def u(self, t=0.0, epar=None, analysis=None):
        if len(self.history) == 0 or analysis not in ('tran', 'transient', 'time'):
            return self.toolkit.zeros(6)
            
        td = self.iparv.TD
        Z0 = self.iparv.Z0
        
        v1_past, v2_past, i1_past, i2_past = self._interpolate_history(t - td)
        
        e1 = v2_past + Z0 * i2_past
        e2 = v1_past + Z0 * i1_past
        
        out = [0.0]*6
        out[4] = -e1
        out[5] = -e2
        return self.toolkit.array(out)

            
        if analysis not in ('tran', 'transient', 'time'):
            return out
            
        td = self.iparv.TD
        Z0 = self.iparv.Z0
        
        # E1(t) = V2(t-TD) + Z0 * I2(t-TD)
        # E2(t) = V1(t-TD) + Z0 * I1(t-TD)
        
        v1_past, v2_past, i1_past, i2_past = self._interpolate_history(t - td)
        
        e1 = v2_past + Z0 * i2_past
        e2 = v1_past + Z0 * i1_past
        
        # u is added to i(x), so F = Gx + u = 0.
        # We want Gx = E, so u = -E
        out[4] = -e1
        out[5] = -e2
        
        return out
