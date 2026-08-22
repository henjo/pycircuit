# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from __future__ import division
import math
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
    def signal_scale(self):
        ## A uic run steps this node from 0 to v at the start.
        return (abs(self.iparv.v), 0.0)

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

    def signal_scale(self):
        ## The 'auto' excursion bound's static anchor: amplitude plus
        ## offset excursion from a uic start.
        return (abs(self.iparv.va) + abs(self.iparv.vo), 0.0)

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
    def signal_scale(self):
        return (0.0, abs(self.iparv.i))

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

    def signal_scale(self):
        return (0.0, abs(self.iparv.ia) + abs(self.iparv.io))

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

    def signal_scale(self):
        i1, i2 = self.iparv.i1, self.iparv.i2
        return (0.0, max(abs(i1), abs(i2), abs(i2 - i1)))
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

    def signal_scale(self):
        v1, v2 = self.iparv.v1, self.iparv.v2
        return (max(abs(v1), abs(v2), abs(v2 - v1)), 0.0)
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
        branchindex = -1 ## add last in matrix
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name])
             for name in self.terminals)

        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (outpindex, branchindex, 1),
             (outnindex, branchindex, -1),
             (branchindex, outpindex, -1),
             (branchindex, outnindex, 1),
             (branchindex, inpindex, self.iparv.g),
             (branchindex, innindex, -self.iparv.g),
            ])
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
        branchindexK = 5 ## add last in matrix
        branchindexJ = 4 ## add last in matrix
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name])
             for name in self.terminals)




        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (outpindex, branchindexK, 1),
             (outnindex, branchindexK, -1),
             (branchindexK, outpindex, 1),
             (branchindexK, outnindex, -1),
             (inpindex, branchindexJ, 1),
             (innindex, branchindexJ, -1),
             (branchindexJ, inpindex, 1),
             (branchindexJ, innindex, -1),
             (branchindexJ, branchindexK, -self.iparv.r),
            ])
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
        ## Stamp construction is the TOOLKIT'S job (matrix_from_entries):
        ## this element's stamp broke twice in one day when it assumed the
        ## backend -- in-place assignment on the (immutable) JAX arrays,
        ## then a float numpy intermediate that rejected symbolic gm.  Each
        ## toolkit now builds its own matrix from the (row, col, value)
        ## entries, which is the whole decision the element should make.
        gm=self.iparv.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        self._G = self.toolkit.matrix_from_entries(
            (n, n),
            [(outpindex, inpindex, gm), (outpindex, innindex, -gm),
             (outnindex, inpindex, -gm), (outnindex, innindex, gm)])


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
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))

        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (outpindex, branchindex, 1),
             (outnindex, branchindex, -1),
             (branchindex, inpindex, 1),
             (branchindex, innindex, -1),
            ])
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
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (inpindex, branchindex, self.iparv.n),
             (innindex, branchindex, -self.iparv.n),
             (outpindex, branchindex, 1),
             (outnindex, branchindex, -1),
             (branchindex, outpindex, self.iparv.n),
             (branchindex, outnindex, -self.iparv.n),
             (branchindex, inpindex, -1),
             (branchindex, innindex, 1),
            ])
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
        gm=self.iparv.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        # 
        #        
        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (outpindex, inpindex, -gm),
             (outpindex, innindex, gm),
             (outnindex, inpindex, gm),
             (outnindex, innindex, -gm),
             (inpindex, outpindex, gm),
             (inpindex, outnindex, -gm),
             (innindex, outpindex, -gm),
             (innindex, outnindex, gm),
            ])
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
        G = self.toolkit.matrix_from_entries(
            (n,n),
            [
             (outpindex, branchindex, 1),
             (outnindex, branchindex, -1),
             (branchindex, outpindex, -1),
             (branchindex, outnindex, 1),
             (branchindex, inpindex, g_limit*self.iparv.g),
             (branchindex, innindex, -g_limit*self.iparv.g),
            ])
        return G

def floored_wrap(y, modulus, offset, toolkit):
    """The Verilog-A ``idtmod`` wrap: fold ``y`` into ``[offset, offset+modulus)``.

    ``y - modulus*floor((y - offset)/modulus)`` -- a *floored* modulo, so the
    result lands in the half-open range for negative ``y`` too, and the LRM
    congruence ``y = n*modulus + k`` holds with integer ``n``.  A C-style
    truncating ``%`` breaks both properties for negative integrals (see
    idtmod.md section 1.1).

    When ``y - offset`` sits within rounding of an exact multiple of the
    modulus, the float result can land one ulp outside the half-open range;
    consumers of a wrapped value are periodic in it, so this is harmless
    and deliberately not patched over with clamping branches.
    """
    return y - modulus * toolkit.floor((y - offset) / modulus)


class _WrapEvents:
    """Phase-3 wrap breakpoints (idtmod.md 5.3), CPU solvers only.

    The wrapped output is a sawtooth whose corner the adaptive controller
    would otherwise discover by rejection; predicting the crossing lets the
    solver LAND a step boundary on it and order-drop across the corner.
    Shared by the scalar (``Idtmod``) and phasor-pair (``IdtmodCircular``)
    circular integrators: both expose the wrapped value via
    ``_event_output(x)`` and both have output rate ``v_ip - v_in`` exactly.
    """

    def _event_output(self, x):
        """The wrapped output value at the element-local state ``x``."""
        raise NotImplementedError

    def accept_step(self, t, x, epar):
        """Cache what the crossing prediction needs, in gauge-invariant
        form: the WRAPPED output (unchanged by the Phase-2 shift, whenever
        it runs relative to this call) and the rate read off the element's
        own ODE right-hand side -- ``d(out)/dt = v_ip - v_in`` exactly, no
        finite differencing, no two-point history."""
        m = self.iparv.modulus
        if m is None or m <= 0:
            self._bp_cache = None
            return
        self._bp_cache = (float(t),
                          float(self._event_output(x)),
                          float(x[0] - x[1]))

    def reset_state(self, epar=None):
        self._bp_cache = None

    def next_event(self, t):
        """Predicted time the output next crosses a wrap boundary.

        Linear prediction from the last ACCEPTED point; it only needs to
        BRACKET the corner -- the solver's break-step machinery order-drops
        across whatever lands there.  ``inf`` before the first accepted
        step, in idt-degradation mode, and for a (near-)stalled phase.
        Strictly greater than ``t``, per the ``next_event`` contract; a
        prediction already at or behind ``t`` (the solver is sitting on the
        corner) advances by whole periods.
        """
        cache = getattr(self, '_bp_cache', None)
        if cache is None:
            return self.toolkit.inf
        t0, phase, rate = cache
        m = self.iparv.modulus
        offset = self.iparv.offset
        if rate == 0.0:
            return self.toolkit.inf
        if rate > 0.0:
            gap = (offset + m) - phase   # in (0, m]: the range is half-open
        else:
            gap = phase - offset         # in [0, m): can sit ON the corner
            if gap <= 0.0:
                gap = m
        period = m / abs(rate)
        t_cross = t0 + gap / abs(rate)
        if t_cross <= t:
            t_cross += (math.floor((t - t_cross) / period) + 1.0) * period
            if t_cross <= t:     # float residue at extreme t/period ratios
                return self.toolkit.inf
        return t_cross


class _IdtBase(Circuit):
    """Shared topology of the ``Idt``/``Idtmod`` integrators.

    One private node carries the state ``s = -integral(v_iplus - v_iminus)dt``
    (row ``idt_node``: ``(v_ip - v_in) + ds/dt = 0`` via the unit C-diagonal),
    and the output branch row defines ``v_out = wrap(-s)``, where ``wrap`` is
    the identity for ``Idt`` and the LRM floored modulo for ``Idtmod``.

    DC semantics follow the LRM (see idtmod.md section 5.1): with an ``ic``
    given, the DC solve pins the state so the output is ``wrap(ic)`` -- the
    row is switched when the analysis marks ``epar.analysis_kind == 'dc'``.
    Without an ``ic`` the integrator keeps its integral row at DC, which is
    the LRM's no-ic branch: the operating point exists only if feedback
    forces the integrand to zero (otherwise use ``uic=True``).
    """

    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)

    instparams = [Parameter(name='ic',
                            desc='Initial output value; forces the DC solution '
                                 'and seeds uic=True. None keeps the LRM no-ic '
                                 'behaviour (integrand forced to zero at DC)',
                            unit='V', default=None)]

    ## The initial condition lives on the element's own private state row --
    ## neither a branch current (`L`) nor a node-voltage difference (`C`) --
    ## so it is applied by direct assignment through `state_ic()`.
    IC_KIND = 'state'

    ## `update()` is first called from `Circuit.__init__` before the private
    ## node exists; this sentinel makes that call a no-op (see `__init__`).
    _idt_index = None

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self._idt_index = self.nodes.index(self.add_node('idt_node'))
        self._branch_index = self.n - 1
        self.update(self.ipar)

    def _wrap(self, y):
        """Output map applied to the integral ``-s``; identity here."""
        return y

    def update(self, subject):
        if self._idt_index is None:
            return
        tk = self.toolkit
        n = self.n
        idt, br = self._idt_index, self._branch_index
        ip, in_, op, on = (self.nodes.index(self.nodenames[name])
                           for name in self.terminals)

        ## Built via matrix_from_entries so each toolkit gets its native
        ## entry type -- the symbolic toolkit needs exact integers (1/s, not
        ## 1.0/s).  The entries: integral row `(v_ip - v_in) + ds/dt = 0`,
        ## branch current into the output nodes, and the output row
        ## `-v_op + v_on + wrap(-s) = 0`, whose (br, idt) slope -1 is the
        ## exact derivative of the wrap almost everywhere.
        common = [(op, br, 1), (on, br, -1),
                  (br, op, -1), (br, on, 1)]
        tran_rows = [(idt, ip, 1), (idt, in_, -1)]
        dc_rows = [(idt, idt, 1)]   # DC pin: `s + wrap(ic) = 0` with u
        slope = [(br, idt, -1)]

        self._G = tk.matrix_from_entries((n, n), tran_rows + common + slope)
        self._G_dc = tk.matrix_from_entries((n, n), dc_rows + common + slope)

        ## `i()` computes the branch row's wrap term itself; these carry the
        ## linear remainder (the (br, idt) entry removed) plus a one-hot to
        ## place the wrap -- pure array arithmetic, no in-place mutation, so
        ## the same code runs under the numeric and JAX toolkits.
        self._G0 = tk.matrix_from_entries((n, n), tran_rows + common)
        self._G0_dc = tk.matrix_from_entries((n, n), dc_rows + common)
        e = [0.0] * n
        e[br] = 1.0
        self._e_branch = tk.array(e)

        self._C = tk.matrix_from_entries((n, n), [(idt, idt, 1)])

    def _dc_pinned(self, epar):
        return (self.iparv.ic is not None
                and getattr(epar, 'analysis_kind', None) == 'dc')

    def state_ic(self):
        """``[(local_row, value)]`` seeds for ``uic=True`` (IC_KIND='state')."""
        ic = self.iparv.ic
        if ic is None:
            return []
        return [(self._idt_index, -self._wrap(ic))]

    def G(self, x, epar=defaultepar, params_tree=None):
        if self._dc_pinned(epar):
            return self._G_dc
        return self._G

    def C(self, x, epar=defaultepar, params_tree=None):
        return self._C

    def i(self, x, epar=defaultepar, params_tree=None):
        G0 = self._G0_dc if self._dc_pinned(epar) else self._G0
        return (self.toolkit.dot(G0, x)
                + self._e_branch * self._wrap(-x[self._idt_index]))

    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        if self._dc_pinned(epar):
            ## Residual row idt is `s + wrap(ic) = 0` -> s = -wrap(ic), so the
            ## branch row's `wrap(-s)` lands the output on wrap(ic) exactly.
            u = [0.0] * self.n
            u[self._idt_index] = self._wrap(self.iparv.ic)
            return self.toolkit.array(u)
        return self.toolkit.zeros(self.n)


class Idt(_IdtBase):
    """Integrator

    Output voltage is the time integral of input voltage.
    """

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        ## q on the private state row only: q[idt_node] = s.
        return toolkit.array([0.0, 0.0, 0.0, 0.0, x[4], 0.0])

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        i_br = x[5]
        return toolkit.array([0.0, 0.0, i_br, -i_br,
                              x[0] - x[1],
                              -(x[2] - x[3]) - x[4]])


class Idtmod(_WrapEvents, _IdtBase):
    """Circular (modulo) integrator -- Verilog-A ``idtmod``

    Output voltage is the time integral of the input voltage folded into
    ``[offset, offset + modulus)`` with the LRM floored wrap (see idtmod.md).
    ``modulus <= 0`` degrades to a plain ``Idt``, per the LRM's "not
    specified" clause.

    >>> import pycircuit.circuit._numeric as numeric
    >>> from pycircuit.circuit.transient import Transient
    >>> c = SubCircuit()
    >>> nin, nout = c.add_nodes('in', 'out')
    >>> c['vin'] = VS(nin, gnd, v=1.0)
    >>> c['R'] = R(nout, gnd, r=1e3)
    >>> c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)

    Without an ``ic``, ``uic=True`` is required: the LRM's no-ic branch
    demands the integrand be forced to zero at DC, so a driven integrator's
    bias solve is singular. Giving ``ic=`` pins the DC solution instead.

    ``fixed_timestep=True`` is deliberate: this example is documenting the
    modulo wrap, so it wants output on the uniform grid it asks for. An
    adaptive run opens at ``timestep*1e-3`` and grows from there (see
    ``firststep``), which is right for accuracy and wrong for a doctest
    about output values.

    >>> tran = Transient(c, toolkit=numeric, uic=True)
    >>> result = tran.solve(tend=1.5, timestep=0.5, fixed_timestep=True)
    >>> result.v(nout).y
    array([0. , 0.5, 0. , 0.5])
    """
    instparams = _IdtBase.instparams + [
                  Parameter(name='modulus', desc='Output modulus', unit='V/V',
                            default=1.),
                  Parameter(name='offset', desc='offset voltage', unit='V',
                            default=0.)]

    linear = False

    def _wrap(self, y):
        m = self.iparv.modulus
        if m is None or m <= 0:
            return y
        return floored_wrap(y, m, self.iparv.offset, self.toolkit)

    def periodic_states(self):
        """The state ``s = -phase`` is defined only up to ``n*modulus``: the
        output is ``wrap(-s)``, periodic in ``s``, so the solvers may keep
        ``s`` bounded by exact ``n*modulus`` translations (idtmod.md 5.2).
        The window ``[-(offset+modulus), -offset)`` mirrors the output range
        through the sign convention, so the wrap in ``i()`` never sees more
        than one modulus of excursion between accepted steps.

        NOT declared when ``modulus <= 0``: that is idt-degradation mode, a
        plain integral whose absolute value is meaningful and must never be
        wrapped.
        """
        m = self.iparv.modulus
        if m is None or m <= 0:
            return []
        return [(self._idt_index, m, -(self.iparv.offset + m))]

    def _event_output(self, x):
        return self._wrap(-x[self._idt_index])

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        return toolkit.array([0.0, 0.0, 0.0, 0.0, x[4], 0.0])

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        ## Batched path: `modulus > 0` is required here -- the degrade-to-idt
        ## branch is data-dependent and cannot be traced.
        m = params.get('modulus', 1.0)
        o = params.get('offset', 0.0)
        i_br = x[5]
        return toolkit.array([0.0, 0.0, i_br, -i_br,
                              x[0] - x[1],
                              -(x[2] - x[3]) + floored_wrap(-x[4], m, o,
                                                            toolkit)])


class IdtmodCircular(_WrapEvents, Circuit):
    """Two-state circular integrator -- ``idtmod`` on the unit phasor.

    Where ``Idtmod`` integrates a scalar phase and folds it (idtmod.md sec.
    5.2), this element integrates the phase's unit phasor
    ``(c, s) = (cos(2*pi*k/modulus), sin(2*pi*k/modulus))``::

        dc/dt = -w*s - gamma*|w|*c*(r^2 - 1)
        ds/dt = +w*c - gamma*|w|*s*(r^2 - 1),   w = 2*pi*(v_ip - v_in)/modulus

    and recovers the wrapped output ``k = floored_wrap(modulus *
    atan2(s, c) / (2*pi), modulus, offset)``.  The STATE is smooth for all
    time -- no wrap, no gauge shift, no discontinuity anywhere in the
    trajectory; the sawtooth exists only in the output map, where the
    ``atan2`` branch cut is exactly cancelled by the wrap and the remaining
    corner is the one ``idtmod`` semantics require (handled by the shared
    ``_WrapEvents`` breakpoints).

    The ``gamma`` terms are the Baumgarte-style orbit correction (idtmod.md
    sec. 7): numerical integration does not conserve the circle invariant
    ``r^2 = c^2 + s^2 = 1`` (Gear-type methods damp the radius, forward
    rules grow it), and the correction makes the unit circle exponentially
    attracting with rate ``2*gamma*|w|``: ``d(r^2)/dt = -2*gamma*|w|*r^2*
    (r^2 - 1)``.  Scaling by ``|w|`` makes the gain dimensionless --
    ``gamma`` is a per-radian-travelled recovery rate, so behaviour does not
    depend on the oscillation frequency.  ``gamma=0`` disables it (useful to
    demonstrate the drift; trapezoidal integration conserves the circle for
    constant ``w`` even then).

    Semantics notes vs ``Idtmod``: ``modulus`` must be positive (a phasor
    cannot represent an unbounded plain integral, so there is no
    idt-degradation mode); ``ic`` defaults to 0 and ALWAYS defines the DC
    solution and the ``uic`` seed -- the phasor must start on the circle
    (``(c, s) = (0, 0)`` is a degenerate equilibrium the dynamics cannot
    leave).  Not supported on the symbolic toolkit (``atan2``/``floor``).
    Under ``JAXTransient`` use ``uic=True`` (or ``x0``): its internal DC
    solve does not carry the ``epar.analysis_kind`` pin and would settle on
    the degenerate centre.
    """

    terminals = ('iplus', 'iminus', 'oplus', 'ominus')
    branches = (Branch(Node('oplus'), Node('ominus')),)
    linear = False

    instparams = [Parameter(name='ic',
                            desc='Initial output value; forces the DC '
                                 'solution and seeds uic=True',
                            unit='V', default=0.),
                  Parameter(name='modulus', desc='Output modulus (> 0)',
                            unit='V/V', default=1.),
                  Parameter(name='offset', desc='offset voltage', unit='V',
                            default=0.),
                  Parameter(name='gamma',
                            desc='Baumgarte orbit-correction gain, per '
                                 'radian of phase travel (0 disables)',
                            unit='', default=1.)]

    IC_KIND = 'state'

    _cos_index = None

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self._cos_index = self.nodes.index(self.add_node('cos_node'))
        self._sin_index = self.nodes.index(self.add_node('sin_node'))
        self._branch_index = self.n - 1
        self.update(self.ipar)

    def update(self, subject):
        if self._cos_index is None:
            return
        m = self.iparv.modulus
        if m is None or m <= 0:
            raise ValueError(
                'IdtmodCircular requires modulus > 0: a phasor pair cannot '
                'represent an unbounded plain integral (use Idt/Idtmod for '
                'the degradation mode).')
        n = self.n
        self._C = self.toolkit.matrix_from_entries(
            (n, n), [(self._cos_index, self._cos_index, 1),
                     (self._sin_index, self._sin_index, 1)])

    def _params(self):
        return {'modulus': self.iparv.modulus, 'offset': self.iparv.offset,
                'gamma': self.iparv.gamma, 'ic': self.iparv.ic}

    def _theta0(self):
        """Phasor angle whose recovered output is ``wrap(ic)``."""
        m = self.iparv.modulus
        kw = floored_wrap(self.iparv.ic, m, self.iparv.offset, self.toolkit)
        return 2.0 * self.toolkit.pi * kw / m

    def state_ic(self):
        th = self._theta0()
        return [(self._cos_index, math.cos(th)),
                (self._sin_index, math.sin(th))]

    def _event_output(self, x):
        m = self.iparv.modulus
        return floored_wrap(
            m * math.atan2(float(x[self._sin_index]),
                           float(x[self._cos_index])) / (2.0 * math.pi),
            m, self.iparv.offset, self.toolkit)

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        return toolkit.array([0.0, 0.0, 0.0, 0.0, x[4], x[5], 0.0])

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        m = params.get('modulus', 1.0)
        o = params.get('offset', 0.0)
        c, s, i_br = x[4], x[5], x[6]
        two_pi = 2.0 * toolkit.pi
        k = floored_wrap(m * toolkit.arctan2(s, c) / two_pi, m, o, toolkit)
        out_row = -(x[2] - x[3]) + k
        if getattr(epar, 'analysis_kind', None) == 'dc':
            ## The ic pin, LRM-style: the phasor lands on the circle at the
            ## angle whose output is wrap(ic).  Folded into i() (no u part),
            ## so the finite-difference/autodiff Jacobian of this same
            ## function is the DC stamp too.
            ic = params.get('ic', 0.0)
            kw = floored_wrap(ic, m, o, toolkit)
            th = two_pi * kw / m
            return toolkit.array([0.0, 0.0, i_br, -i_br,
                                  c - toolkit.cos(th),
                                  s - toolkit.sin(th),
                                  out_row])
        gamma = params.get('gamma', 1.0)
        w = two_pi * (x[0] - x[1]) / m
        corr = gamma * toolkit.abs(w) * (c * c + s * s - 1.0)
        return toolkit.array([0.0, 0.0, i_br, -i_br,
                              w * s + corr * c,
                              -w * c + corr * s,
                              out_row])

    def i(self, x, epar=defaultepar, params_tree=None):
        return self.eval_i_pure(x, self._params(), epar, self.toolkit)

    def G(self, x, epar=defaultepar, params_tree=None):
        ## Single source of truth: autodiff under JAX, central differences
        ## under the numeric toolkit (the Semiconductor pattern).  The
        ## |w| kink and the sawtooth corner make G exact a.e. only, same
        ## story as Idtmod's wrap slope.
        return self.toolkit.jacobian(self.eval_i_pure, x, self._params(),
                                     epar)

    def C(self, x, epar=defaultepar, params_tree=None):
        return self._C

    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        return self.toolkit.zeros(self.n)


class IdtmodQuadrature(Circuit):
    """Quadrature phasor integrator -- the two states ARE the outputs.

    The same Baumgarte-corrected phasor dynamics as ``IdtmodCircular``
    (idtmod.md sec. 7), but with NO wrapped output and NO ``atan2``: the
    two output branches carry the states directly,

        v(cplus,cminus) = amplitude * cos(2*pi*phase/modulus)
        v(splus,sminus) = amplitude * sin(2*pi*phase/modulus)

    where ``d(phase)/dt = v(iplus,iminus)``.  Every node voltage and every
    state is smooth for all time -- there is no sawtooth anywhere, hence no
    wrap breakpoints, no events, and nothing for a step controller to
    negotiate.  This is the natural macromodel for a sinusoidal VCO, an
    I/Q source, or a mixer LO, and -- the reason it exists as a separate
    element -- the representation suited to periodic-steady-state /
    shooting analysis when one is implemented: over one output period the
    state vector returns EXACTLY to itself, which the scalar ``Idtmod``
    phase (advancing by one modulus per period, or jumping at the wrap)
    structurally cannot do.

    Same accuracy regime as ``IdtmodCircular``: the method's per-cycle
    phase lag at the carrier rate (idtmod.md 7.6) -- use ``Idtmod`` when
    accumulated phase is the measurement.  Same semantics: ``modulus`` must
    be positive; ``ic`` (an initial phase, in output-value units: the
    starting angle is ``2*pi*ic/modulus``) always pins the DC solution and
    seeds ``uic=True``; ``gamma`` is the per-radian orbit-correction gain;
    under ``JAXTransient`` use ``uic=True`` or ``x0``.  Not supported on
    the symbolic toolkit.
    """

    terminals = ('iplus', 'iminus', 'cplus', 'cminus', 'splus', 'sminus')
    branches = (Branch(Node('cplus'), Node('cminus')),
                Branch(Node('splus'), Node('sminus')))
    linear = False

    instparams = [Parameter(name='ic',
                            desc='Initial phase (output-value units); '
                                 'forces the DC solution and seeds uic',
                            unit='V', default=0.),
                  Parameter(name='modulus',
                            desc='Phase units per full rotation (> 0)',
                            unit='V/V', default=1.),
                  Parameter(name='gamma',
                            desc='Baumgarte orbit-correction gain, per '
                                 'radian of phase travel (0 disables)',
                            unit='', default=1.),
                  Parameter(name='amplitude', desc='Output amplitude',
                            unit='V', default=1.)]

    IC_KIND = 'state'

    _cos_index = None

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self._cos_index = self.nodes.index(self.add_node('cos_node'))
        self._sin_index = self.nodes.index(self.add_node('sin_node'))
        self.update(self.ipar)

    def update(self, subject):
        if self._cos_index is None:
            return
        m = self.iparv.modulus
        if m is None or m <= 0:
            raise ValueError(
                'IdtmodQuadrature requires modulus > 0: a phasor pair '
                'cannot represent an unbounded plain integral.')
        n = self.n
        self._C = self.toolkit.matrix_from_entries(
            (n, n), [(self._cos_index, self._cos_index, 1),
                     (self._sin_index, self._sin_index, 1)])

    def _params(self):
        return {'modulus': self.iparv.modulus, 'gamma': self.iparv.gamma,
                'ic': self.iparv.ic, 'amplitude': self.iparv.amplitude}

    def state_ic(self):
        th = 2.0 * math.pi * self.iparv.ic / self.iparv.modulus
        return [(self._cos_index, math.cos(th)),
                (self._sin_index, math.sin(th))]

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        return toolkit.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              x[6], x[7], 0.0, 0.0])

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        ## Local x: [ip, in, cp, cn, sp, sn, c, s, i_cbr, i_sbr].
        m = params.get('modulus', 1.0)
        A = params.get('amplitude', 1.0)
        c, s = x[6], x[7]
        i_cbr, i_sbr = x[8], x[9]
        out_c = -(x[2] - x[3]) + A * c
        out_s = -(x[4] - x[5]) + A * s
        if getattr(epar, 'analysis_kind', None) == 'dc':
            th = 2.0 * toolkit.pi * params.get('ic', 0.0) / m
            return toolkit.array([0.0, 0.0, i_cbr, -i_cbr, i_sbr, -i_sbr,
                                  c - toolkit.cos(th),
                                  s - toolkit.sin(th),
                                  out_c, out_s])
        gamma = params.get('gamma', 1.0)
        w = 2.0 * toolkit.pi * (x[0] - x[1]) / m
        corr = gamma * toolkit.abs(w) * (c * c + s * s - 1.0)
        return toolkit.array([0.0, 0.0, i_cbr, -i_cbr, i_sbr, -i_sbr,
                              w * s + corr * c,
                              -w * c + corr * s,
                              out_c, out_s])

    def i(self, x, epar=defaultepar, params_tree=None):
        return self.eval_i_pure(x, self._params(), epar, self.toolkit)

    def G(self, x, epar=defaultepar, params_tree=None):
        return self.toolkit.jacobian(self.eval_i_pure, x, self._params(),
                                     epar)

    def C(self, x, epar=defaultepar, params_tree=None):
        return self._C

    def u(self, t=0.0, epar=defaultepar, analysis=None, params_tree=None):
        return self.toolkit.zeros(self.n)

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
        # Branch 1 KCL contributions
        # Branch 1 voltage equation
        
        # Branch 2 KCL contributions
        # Branch 2 voltage equation
        
        G = self.toolkit.matrix_from_entries(
            (n, n),
            [
             (0, 4, 1),
             (1, 4, -1),
             (4, 0, 1),
             (4, 1, -1),
             (2, 5, 1),
             (3, 5, -1),
             (5, 2, 1),
             (5, 3, -1),
            ])
        self._G = G
        
        L1, L2, K = self.iparv.L1, self.iparv.L2, self.iparv.K
        M = K * self.toolkit.sqrt(L1 * L2)
        
        C = self.toolkit.matrix_from_entries(
            (n, n),
            [
             (4, 4, -L1),
             (5, 5, -L2),
             (4, 5, -M),
             (5, 4, -M),
            ])
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
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) for name in self.terminals)
        branchindex = 4
        
        # Branch KCL (Input acts as short, measuring current)
        
        # Branch equation: V_inp - V_inn = 0
        
        # Output current injections
        
        G = self.toolkit.matrix_from_entries(
            (n, n),
            [
             (inpindex, branchindex, 1),
             (innindex, branchindex, -1),
             (branchindex, inpindex, 1),
             (branchindex, innindex, -1),
             (outpindex, branchindex, self.iparv.F),
             (outnindex, branchindex, -self.iparv.F),
            ])
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
        ## Row 4 is the control branch's DEFINING equation, v_cp - v_cm = 0
        ## -- it was returned as a constant 0.0, so on any autodiff toolkit
        ## the KVL row vanished from the jacobian (measured: G[4,2]/G[4,3]
        ## = 0 under JAX against the manual stamp's 1/-1), and the residual
        ## never enforced the equation the manual G claims.  Found by the
        ## P24 cross-toolkit stamp gate.
        return toolkit.array([i_val, -i_val, i_ctrl, -i_ctrl, x[2] - x[3]])

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
        
        
        # d(i_val)/dv = g
        
        # d(i_val)/di_ctrl = v * dg_di
        g_i = v * dg_di
        
        # cp/cm KCL contributions
        
        # Branch equation V_cp - V_cm = 0
        
        G_mat = self.toolkit.matrix_from_entries(
            (5, 5),
            [
             (0, 0, g),
             (0, 1, -g),
             (1, 0, -g),
             (1, 1, g),
             (0, 4, g_i),
             (1, 4, -g_i),
             (2, 4, 1.0),
             (3, 4, -1.0),
             (4, 2, 1.0),
             (4, 3, -1.0),
            ])
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

        G_mat = self.toolkit.matrix_from_entries(
            (4, 4),
            [
             (2, 0, di_dv),
             (2, 1, -di_dv),
             (3, 0, -di_dv),
             (3, 1, di_dv),
            ])
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

        C_mat = self.toolkit.matrix_from_entries(
            (4, 4),
            [
             (2, 0, dq_dv),
             (2, 1, -dq_dv),
             (3, 0, -dq_dv),
             (3, 1, dq_dv),
            ])
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
                    
                    ## MONOTONE LIMITER, per channel.  A quadratic whose third
                    ## point sits across a kink in the recorded waveform
                    ## OVERSHOOTS inside the bracket: measured on a pulsed line
                    ## (tr=2e-10, TD=1e-9), the stencil (ramp point, corner,
                    ## flat point) put the reflected EMF at 1.009 where the
                    ## recorded samples bound it by 1.000 -- and the coupled
                    ## path accepted a solution against that phantom, which
                    ## the NEXT solve (its stencil changed by accept_step's
                    ## pruning) could never reconcile: the solution-LTE hit an
                    ## h-independent floor of exactly the pollution, and the
                    ## run livelocked.  A value the quadratic puts OUTSIDE the
                    ## bracket's own [min, max] is evidence the stencil spans
                    ## a kink, so that channel falls back to the linear value
                    ## -- which also makes the result indifferent to whether
                    ## the pre-bracket point has been pruned yet.  Smooth
                    ## regions keep the quadratic and its accuracy argument
                    ## (see the note at the top of this method).
                    alpha = (t_target - t1) / (t2 - t1)
                    out = []
                    for y0, y1, y2 in ((v1_0, v1_1, v1_2), (v2_0, v2_1, v2_2),
                                       (i1_0, i1_1, i1_2), (i2_0, i2_1, i2_2)):
                        yq = y0 * L0 + y1 * L1 + y2 * L2
                        lo, hi = (y1, y2) if y1 <= y2 else (y2, y1)
                        if yq < lo or yq > hi:
                            yq = y1 + alpha * (y2 - y1)
                        out.append(yq)
                    return tuple(out)
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

        
    def _interpolate_history_slope(self, t_target):
        ## The DERIVATIVE of the polynomial `_interpolate_history` evaluates --
        ## the same branch structure, differentiated, so `dudt` is exactly
        ## d/dt of `u` and the coupled residual's `p` vector stays consistent
        ## with the residual itself.  Outside the recorded span the value
        ## interpolation clamps to a constant, so the slope there is zero.
        if len(self.history) < 2:
            return 0.0, 0.0, 0.0, 0.0
        if t_target <= self.history[0][0] or t_target >= self.history[-1][0]:
            return 0.0, 0.0, 0.0, 0.0
        for i in range(len(self.history) - 1):
            if self.history[i][0] <= t_target <= self.history[i+1][0]:
                if i > 0:
                    t0, v1_0, v2_0, i1_0, i2_0 = self.history[i-1]
                    t1, v1_1, v2_1, i1_1, i2_1 = self.history[i]
                    t2, v1_2, v2_2, i1_2, i2_2 = self.history[i+1]
                    L0 = (t_target - t1) * (t_target - t2) / ((t0 - t1) * (t0 - t2))
                    L1 = (t_target - t0) * (t_target - t2) / ((t1 - t0) * (t1 - t2))
                    L2 = (t_target - t0) * (t_target - t1) / ((t2 - t0) * (t2 - t1))
                    dL0 = (2.0*t_target - t1 - t2) / ((t0 - t1) * (t0 - t2))
                    dL1 = (2.0*t_target - t0 - t2) / ((t1 - t0) * (t1 - t2))
                    dL2 = (2.0*t_target - t0 - t1) / ((t2 - t0) * (t2 - t1))
                    ## Same per-channel limiter as `_interpolate_history`:
                    ## where the VALUE falls back to linear, the slope must
                    ## too, or `dudt` stops being d/dt of `u`.
                    inv = 1.0 / (t2 - t1)
                    out = []
                    for y0, y1, y2 in ((v1_0, v1_1, v1_2), (v2_0, v2_1, v2_2),
                                       (i1_0, i1_1, i1_2), (i2_0, i2_1, i2_2)):
                        yq = y0 * L0 + y1 * L1 + y2 * L2
                        lo, hi = (y1, y2) if y1 <= y2 else (y2, y1)
                        if yq < lo or yq > hi:
                            out.append((y2 - y1) * inv)
                        else:
                            out.append(y0*dL0 + y1*dL1 + y2*dL2)
                    return tuple(out)
                else:
                    t0, v1_0, v2_0, i1_0, i2_0 = self.history[i]
                    t1, v1_1, v2_1, i1_1, i2_1 = self.history[i+1]
                    inv = 1.0 / (t1 - t0)
                    return ((v1_1 - v1_0) * inv, (v2_1 - v2_0) * inv,
                            (i1_1 - i1_0) * inv, (i2_1 - i2_0) * inv)
        return 0.0, 0.0, 0.0, 0.0

    def dudt(self, t=0.0, epar=None, analysis=None):
        """d/dt of the delayed-history source vector -- STAGE 12B, at last.

        ``u`` builds the reflected EMFs from the history interpolated at
        ``t - TD``; ``du/dt`` is the time derivative of that interpolation,
        evaluated by differentiating the SAME polynomial the value uses
        (`_interpolate_history_slope`), so the coupled path's ``p`` vector is
        consistent with the residual to machine precision rather than to a
        finite-difference step.  This used to raise: inheriting the base
        class's zero would have dropped a real term from ``p`` silently, and
        the raise stood until the term was written.  The refusal test in
        test_residual_dh.py was re-derived into a slope-correctness test when
        this landed.
        """
        if len(self.history) == 0 or analysis not in ('tran', 'transient',
                                                      'time'):
            return self.toolkit.zeros(6)
        td = self.iparv.TD
        Z0 = self.iparv.Z0
        dv1, dv2, di1, di2 = self._interpolate_history_slope(t - td)
        de1 = dv2 + Z0 * di2
        de2 = dv1 + Z0 * di1
        out = [0.0]*6
        out[4] = -de1
        out[5] = -de2
        return self.toolkit.array(out)

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
