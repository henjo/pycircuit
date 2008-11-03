from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, alltrue, maximum, conj, dot, imag, eye
from scipy import optimize
from circuit import Circuit, SubCircuit, VS,IS,R,C,L,Diode, gnd, defaultepar
from pycircuit.waveform import Waveform
from pycircuit.internalresult import InternalResultDict
import numpy
from copy import copy

class NoConvergenceError(Exception):
    pass

def removeRowCol(matrices, n):
    result = []
    for A in matrices:
        for axis in range(len(A.shape)):
            A=delete(A, [n], axis=axis)
        result.append(A)
    return tuple(result)

class Analysis(object):
    @staticmethod
    def linearsolver(*args):

        args = [A.astype('complex') for A in args]
        
        return linalg.solve(*args)

    @staticmethod
    def toMatrix(array): return array.astype('complex')
        
    def __init__(self, circuit):
        self.c = circuit

        self.result = None

    def run(self, *args, **kvargs):
        """Start the analysis and return the result as a Result object"""

        x = self.solve(*args, **kvargs)
        
        result = InternalResultDict()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes),0], nodes):
            result[self.c.getNodeName(node)] =  xvalue
        for i, data in enumerate(zip(x[len(nodes):, 0], self.c.branches)):
            xvalue, branch = data
            result['i' + str(i)] = xvalue

        self.result = result
        return result
    
    def getResult(self):
        return self.result


def fsolve(f, x0, fprime=None, args=(), full_output=False, maxiter=200,
           xtol=1e-6, reltol=1e-4, abstol=1e-12):
    """Solve a multidimensional non-linear system of equations with Newton-Raphson's method

    In each iteration the linear system

    M{J(x_n)(x_{n+1}-x_n) + F(xn) = 0

    is solved and a new value for x is obtained x_{n+1}
    
    """
    
    converged = False
    ier = 2
    for i in xrange(maxiter):
        J = fprime(x0, *args) # TODO: Make sure J is never 0, e.g. by gmin (stepping)
        F = f(x0, *args)
        xdiff = linalg.solve(J, -F)# TODO: Limit xdiff to improve convergence
        x = x0 + xdiff

        if alltrue(abs(xdiff) < reltol * maximum(x, x0) + xtol):
            ier = 1
            mesg = "Success"
            break
        if alltrue(abs(F) < reltol * max(F) + abstol):
            ier = 1
            mesg = "Success"
            break
            
        x0 = x

    if ier == 2:
        mesg = "No convergence. xerror = "+str(xdiff)
    
    infodict = {}
    if full_output:
        return x, infodict, ier, mesg
    else:
        return x
           

class DC(Analysis):
    """DC analyis class
    
    Linear circuit example:
    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> dc = DC(c)
    >>> res = dc.run()
    >>> dc.result.keys()
    ['i0', 'gnd', 'net1']
    >>> dc.result['net1']
    1.5

    Non-linear example:

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['is'] = IS(gnd, n1, i=57e-3)
    >>> c['D'] = Diode(n1, gnd)
    >>> dc = DC(c)
    >>> res = dc.run()
    >>> dc.result.keys()
    ['i0', 'gnd', 'net1']
    >>> dc.result['net1']
    0.7

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['is'] = IS(gnd, n1, i=57e-3)
    >>> c['R'] = R(n1, n2, r=1e1)
    >>> c['D'] = Diode(n2, gnd)
    >>> dc = DC(c)
    >>> res = dc.run()
    >>> dc.result.getSignalNames()
    ['gnd', 'net2', 'net1']
    >>> dc.result['net2']
    0.7

    """
    
    def solve(self, refnode=gnd):
        n=self.c.n
        x = zeros(n)
        G=self.c.G(x)
        # G += eye(G.shape[0])*1e-2 # Could use this to add a Gmin to every node
        U=self.c.U(x)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,U=removeRowCol((G,U), irefnode)

        # Solve i(x) + u = 0
        def func(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            f =  self.c.i(x) + self.c.U(0)
            (f,) = removeRowCol((f,), irefnode)
#            print x,f.T[0]
            return array(f, dtype=float)
        def fprime(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
#            J = self.c.G(array([x]).T)
            J = self.c.G(x)
            (J,) = removeRowCol((J,), irefnode)
#            print "J:",array(J, dtype=float)
            return array(J, dtype=float)
        x0 = zeros(n-1) # Would be good with a better initial guess
        #x0 = zeros(n-1)+0.75 # works good with this guess, for example

        rtol = 1e-4

        self.__SI__ = numpy.zeros(x0.shape,'d')
        iwk = numpy.zeros((100*len(self.__SI__)),'i')
        rwk = numpy.zeros((100*len(self.__SI__)),'d')
        iopt = numpy.zeros((50),'i')
        s_scale = copy(x0)
        
        iopt[2] =  1 # self.nleq2_jacgen #2
        iopt[8] =  0 # self.nleq2_iscal  #0
        iopt[10] = 1 # self.nleq2_mprerr #1
        iopt[30] = 4 # self.nleq2_nonlin #4
        iopt[31] = 0 #self.nleq2_qrank1 #0
        iopt[34] = 0 # self.nleq2_qnscal #0
        iopt[37] = 0 #self.nleq2_ibdamp #0
        iopt[38] = 0 # self.nleq2_iormon #0

        iwk[30] = 50
        
#        import nleq2.nleq2 as nleq2
#        res,s_scale,rtol,iopt,ierr = nleq2.nleq2(func,fprime,x0,s_scale,rtol,iopt,iwk,rwk)
#        print res, ierr

        x, infodict, ier, mesg = fsolve(func, x0, fprime=fprime, full_output=True)

        if ier != 1:
            raise NoConvergenceError(mesg)

        x = x.reshape((n-1,1))
        # Insert reference node voltage
        x = concatenate((x[:irefnode, :], array([[0.0]]), x[irefnode:,:]))
        return x


class Tran_spec(Analysis):
    """Simple transient analyis class.

    Starting a transient analysis class with support for linear
    dynamic elements.

    Using an explicit method (forward euler) because of it's 
    simplicity. Without any other tricks, backward euler (being
    an implicit method) is otherwise more suited for stiff systems 
    of ODEs.
    
    The usual companion models are used.

    i(t) = c*dv/dt
    v(t) = L*di/dt

    forward euler:
    i(n+1) = c/dt*(v(n) - v(n-1)) = geq*v(n) + Ieq
    v(n+1) = L/dt*(i(n) - i(n-1)) = req*i(n) + Veq

    For each time-step:
    (G+Geq)*x(n) + U + Ueq

    Linear circuit example:
    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['Is'] = IS(gnd, n1, i=10)
    >>> c['R1'] = R(n1, gnd, r=1)
    >>> c['R2'] = R(n1, n2, r=1e3)
    >>> c['R3'] = R(n2, gnd, r=100e3)
    >>> c['C'] = C(n2, gnd, c=1e-5)
    >>> tran = Tran_spec(c)
    >>> res = tran.run(tend=20e-3)
    >>> tran.result
    ' roughly 10,6,0?'

    Linear circuit example:
    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['Is'] = IS(gnd, n1, i=10)
    >>> c['R'] = R(n1, gnd, r=1)
    >>> c['C'] = C(n1, gnd, c=1e-6)
    >>> c['L'] = L(n1, gnd, L=1e-3)
    >>> tran = Tran_spec(c)
    >>> nodeset = array([1,0,0])
    >>> res = tran.run(tend=2e-6,timestep=1e-8)
    >>> tran.result
    ' haven't checked yet'

    """
    # Very incomplete TODO-list:
    # solve with implicit method (with fsolve)
    # generalise to using different coefficients for other methods
    # try using telescopic projective solver with explicit method
    # much more work with presenting and storing results
    # try differential evolution for determining steady-state solution
    def run(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-5):

        irefnode=self.c.getNodeIndex(refnode)
        n = self.c.n
        dt = timestep

        if x0 is None:
            x = zeros(n)
        else:
            x = x0

        xold = x #use x(0) as x(-1) in forward euler
        G=self.c.G(x)
        Geq=self.c.C(x)/dt
        U=self.c.U(x)
        Ueq=-dot(Geq,xold) #use x0 as x(-1) in forward euler
        
        #create vector with timepoints and a more fitting dt
        times,dt=numpy.linspace(0,tend,num=int(tend/dt),endpoint=True,retstep=True)

        for t in times:
            G=self.c.G(x)
            Geq=self.c.C(x)/dt
            U=self.c.U(x)
            Ueq=-dot(Geq,xold)
            G+=Geq
            U+=Ueq
            # Refer the voltages to the reference node by removing
            # the rows and columns that corresponds to this node
            G,U=removeRowCol((G,U), irefnode)
            xold=x
            x=linalg.solve(G,-U)
            # Insert reference node voltage
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            #print(x,t)
        self.result=x
        return x

class AC(Analysis):
    """
    AC analysis class
    
    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, vac=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> ac = AC(c)
    >>> res = ac.run(freqs=array([1e6, 2e6]))
    >>> ac.result.keys()
    ['i0', 'gnd', 'net1']
    >>> ac.result['net1']
    Waveform([ 1000000.  2000000.],[ 1.5+0.j  1.5+0.j])
    >>> ac.v(n1, gnd)
    Waveform([ 1000000.  2000000.],[ 1.5+0.j  1.5+0.j])
    >>> ac.i('vs.minus')
    Waveform([ 1000000.  2000000.],[ 0.0015 +9.42477796e-06j  0.0015 +1.88495559e-05j])

    """

    numeric = True

    def run(self, freqs, **kvargs):
        x = self.solve(freqs, **kvargs)
        
        result = InternalResultDict()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes)], nodes):
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue
            result[self.c.getNodeName(node)] = wave
        for i, data in enumerate(zip(x[len(nodes):], self.c.branches)):
            xvalue, branch = data
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue
            result['i' + str(i)] = wave

        self.result = result

        return result

    def v(self, node1, node2):
        """Return the voltage v(node1, node2) from last run"""

        if self.result != None:
            return self.result[self.c.getNodeName(node1)] - \
                   self.result[self.c.getNodeName(node2)]

    def i(self, term):
        """Return terminal current i(term) of a circuit element from last run"""
        if self.result != None:
            branch, sign = self.c.getTerminalBranch(term)
            ibranch = self.c.branches.index(branch)
            return sign * self.result['i%d'%ibranch]
        
    def solve(self, freqs, refnode=gnd, complexfreq = False, U = None):
        n=self.c.n
        
        x = zeros(n) ## FIXME, this should be calculated from the dc analysis
        
        G=self.c.G(x)
        C=self.c.C(x)

        ## Allow for custom stimuli, mainly used by other analyses
        if U == None:
            U=self.c.U(x, analysis='ac')
        else:
            U = U

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,C,U = removeRowCol((G,C,U), irefnode)

        if complexfreq:
            ss = freqs
        else:
            ss = 2j*pi*freqs

        def solvecircuit(s):
            x = self.linearsolver(s*C + G, -U) 

            # Insert reference node voltage
            return concatenate((x[:irefnode], array([0.0]), x[irefnode:]))

        if isiterable(freqs):
            out = [solvecircuit(s) for s in ss]
            # Swap frequency and x-vector dimensions
            return array(out).swapaxes(0,1)
        else:
            return solvecircuit(ss)

class Noise(Analysis):
    """Noise analysis that calculates input and output referred noise.
    
    The analysis is using the adjoint admittance matrix method to calculate the transfers from
    each noise source to the output.
    
    Example, calculate input referred noise of a voltage divider:

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['vs'] = VS(n1, gnd, vac=1.0)
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = Noise(c, inputsrc=c['vs'], outputnodes=(n2, gnd)).run(0)
    >>> res['Svnout']
    (1.4904e-17+0j)
    >>> res['Svninp']
    (1.4904e-15+0j)
    >>> res['gain']
    (0.1+0j)
    
    """
    def __init__(self, circuit, inputsrc=None, outputnodes=None, outputsrc=None):
        """
        Initiate a noise analysis.

        Parameters
        ----------
        circuit : Circuit instance
            The circuit to be analyzed
        inputsrc : VS or IS instance
            A voltage or current source in the circuit where the input noise should be referred to
        outputnodes : tuple
            A tuple with the output nodes (outputpos outputneg)
        outputsrc: VS instance
            The voltage source where the output current noise is measured
        """

        Analysis.__init__(self, circuit)
    
        if not (outputnodes != None or outputsrc != None):
            raise ValueError('Output is not specified')
        elif outputnodes != None and outputsrc != None:
            raise ValueError('Cannot measure both output current and voltage noise')
        
        self.inputsrc = inputsrc
        self.outputnodes = outputnodes
        self.outputsrc = outputsrc
        self.epar = defaultepar

    def run(self, freqs, refnode=gnd, complexfreq=False):
        n = self.c.n
        x = zeros(n) # This should be the x-vector at the DC operating point

        ## Complex frequency variable
        if complexfreq:
            s = freqs
        else:
            s = 2j*pi*freqs

        epar = self.epar
        G = self.c.G(x, epar)
        C = self.c.C(x, epar)
        CY = self.c.CY(x, imag(s),epar)

        # Calculate output voltage noise
        if self.outputnodes != None:
            U = zeros(n, dtype=int)
            ioutp, ioutn = (self.c.getNodeIndex(node) for node in self.outputnodes)
            U[ioutp] = -1
            U[ioutn] = 1
        # Calculate output current noise
        else:
            U = zeros(n, dtype=int)
            ibranch = self.c.getBranchIndex(self.outputsrc.branch)
            U[ibranch] = -1

        ## Refer the voltages to the gnd node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.nodes.index(refnode)
        G,C,U,CY = removeRowCol((G,C,U,CY), irefnode)

        # Calculate the reciprocal G and C matrices
        Yreciprocal = G.T + s*C.T

        Yreciprocal, U = (self.toMatrix(A) for A in (Yreciprocal, U))

        ## Calculate transimpedances from currents in each nodes to output
        zm = self.linearsolver(Yreciprocal, -U)

        xn2out = dot(dot(zm.reshape(1,size(zm)), CY), conj(zm))

        xn2out = xn2out[0]

        # Store results
        result = InternalResultDict()

        if self.outputnodes != None:
            result['Svnout'] = xn2out
        elif self.outputsrc != None:
            result['Sinout'] = xn2out

        # Calculate the gain from the input voltage source by using the transimpedance vector
        # to find the transfer from the branch voltage of the input source to the output
        gain = None
        if isinstance(self.inputsrc, VS):
            gain = self.c.extractI(zm, self.inputsrc.branch, refnode=refnode, refnode_removed=True)
            result['gain'] = gain
            result['Svninp'] = xn2out / gain**2
        elif isinstance(self.inputsrc, IS):
            gain = self.c.extractV(zm, self.inputsrc.getNode('plus'), refnode=refnode, refnode_removed=True)
            result['gain'] = gain
            result['Sininp'] = xn2out / gain**2

        return result


def isiterable(object):
    return hasattr(object,'__iter__')

if __name__ == "__main__":
    import doctest
    doctest.testmod()

#    c = SubCircuit()
#    n1 = c.addNode('net1')
#    c['is'] = IS(gnd, n1, i=100e-6)
#    c['D'] = Diode(n1, gnd)
#    c['D'].G(array([[0,0]]).T)
#    dc = DC(c)
#    res = dc.run()



