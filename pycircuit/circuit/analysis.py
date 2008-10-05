from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, alltrue, maximum
from scipy import optimize
from circuit import Circuit, SubCircuit, VS,IS,R,C,Diode, gnd
from pycircuit.result import Waveform
from pycircuit.internalresult import InternalResultSet, InternalResult
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
    def __init__(self, circuit):
        self.c = circuit

        self.result = None

    def run(self, *args, **kvargs):
        """Start the analysis and return the result as a Result object"""

        x = self.solve(*args, **kvargs)
        
        result = InternalResult()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes),0], nodes):
            result.storeSignal(self.c.getNodeName(node), xvalue)
        for i, data in enumerate(zip(x[len(nodes):, 0], self.c.branches)):
            xvalue, branch = data
            result.storeSignal('i' + str(i), xvalue)

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
        J = fprime(x0, *args)
        F = f(x0, *args)
        xdiff = linalg.solve(J, -F)
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
    >>> dc.result.getSignalNames()
    ['i0', 'gnd', 'net1']
    >>> dc.result.getSignal('net1')
    1.5

    Non-linear example:

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['is'] = IS(gnd, n1, i=100e-6)
    >>> c['D'] = Diode(n1, gnd)
    >>> c['D'].G(array([[0,0]]).T)
    >>> dc = DC(c)
    >>> res = dc.run()
    >>> dc.result.getSignalNames()
    ['i0', 'gnd', 'net1']
    >>> dc.result.getSignal('net1')
    0.7
    """
    
    def solve(self, refnode=gnd):
        n=self.c.n
        x = zeros((n,1))
        G=self.c.G(x)
        U=self.c.U(x)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,U=removeRowCol((G,U), irefnode)

        # Solve i(x) + u = 0
        def func(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            f =  self.c.i(array([x]).T) + self.c.U(0)
            (f,) = removeRowCol((f,), irefnode)
#            print x,f.T[0]
            return array(f.T[0], dtype=float)
        def fprime(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            J = self.c.G(array([x]).T)
            (J,) = removeRowCol((J,), irefnode)
#            print "J:",array(J, dtype=float)
            return array(J, dtype=float)
        x0 = zeros(n-1)

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


class AC(Analysis):
    """
    AC analysis class
    
    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> ac = AC(c)
    >>> res = ac.run(freqs=array([1e6, 2e6]))
    >>> ac.result.getSignalNames()
    ['i0', 'gnd', 'net1']
    >>> ac.result.getSignal('net1')
    Waveform([ 1000000.  2000000.],[ 1.5+0.j  1.5+0.j])
    >>> ac.v(n1, gnd)
    Waveform([ 1000000.  2000000.],[ 1.5+0.j  1.5+0.j])
    >>> ac.i('vs.minus')
    Waveform([ 1000000.  2000000.],[ 0.0015 +9.42477796e-06j  0.0015 +1.88495559e-05j])

    """

    @staticmethod
    def linearsolver(*args):
        return linalg.solve(*args)
        
    numeric = True

    def run(self, freqs, **kvargs):
        x = self.solve(freqs, **kvargs)
        
        result = InternalResult()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes)], nodes):
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue[0]
            result.storeSignal(self.c.getNodeName(node), wave)
        for i, data in enumerate(zip(x[len(nodes):], self.c.branches)):
            xvalue, branch = data
            if isiterable(freqs):
                wave = Waveform(freqs, xvalue)
            else:
                wave = xvalue[0]
            result.storeSignal('i' + str(i),wave)

        self.result = result

        return result

    def v(self, node1, node2):
        """Return the voltage v(node1, node2) from last run"""

        if self.result != None:
            return self.result.getSignal(self.c.getNodeName(node1)) - \
                   self.result.getSignal(self.c.getNodeName(node2))

    def i(self, term):
        """Return terminal current i(term) of a circuit element from last run"""
        if self.result != None:
            branch, sign = self.c.getTerminalBranch(term)
            ibranch = self.c.branches.index(branch)
            return sign * self.result.getSignal('i%d'%ibranch)
        
    def solve(self, freqs, refnode=gnd, complexfreq = False):
        n=self.c.n
        
        x = zeros((n,1)) ## FIXME, this should be calculated from the dc analysis
        
        G=self.c.G(x)
        C=self.c.C(x)
        U=self.c.U(x)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,C,U = removeRowCol((G,C,U), irefnode)
        if self.numeric:
            G,C,U = (A.astype(complex) for A in (G,C,U))

        out = []

        

        if complexfreq:
            ss = freqs
        else:
            ss = 2j*pi*freqs

        def solvecircuit(s):
            solver = self.linearsolver

            x = solver(s*C + G, -U) 

            # Insert reference node voltage
            return concatenate((x[:irefnode, :], array([[0.0]]), x[irefnode:,:]))

        if isiterable(freqs):
            out = [solvecircuit(s) for s in ss]
            # Swap frequency and x-vector dimensions
            return array(out)[:,:,0].swapaxes(0,1)
        else:
            return solvecircuit(ss)

class TwoPort(Analysis):
    """Analysis to find the 2-ports parameters of a circuit

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = TwoPort(c, n1, gnd, n2, gnd).run(freqs = array([0]))
    >>> res.getSignal('mu').y[0]
    (0.1+0j)
    >>> res.getSignal('gamma').y[0]
    (0.000111111111111+0j)
    >>> res.getSignal('zeta').y[0]
    (1000+0j)
    >>> res.getSignal('beta').y[0]
    (1+0j)
    
    """
    
    ACAnalysis = AC
    
    def __init__(self, circuit, inp, inn, outp, outn):
        self.c = circuit

        self.ports = inp, inn, outp, outn
        
    def run(self, freqs, **kvargs):
        result = InternalResult()

        abcd = self.solve(freqs, **kvargs)
        result.storeSignal('ABCD', abcd)

        result.storeSignal('mu', 1/abcd[0,0])
        result.storeSignal('gamma', 1/abcd[0,1])
        result.storeSignal('zeta', 1/abcd[1,0])
        result.storeSignal('beta', 1/abcd[1,1])

        self.result = result

        return result

    def solve(self, freqs, refnode = gnd, complexfreq = False):
        inp, inn, outp, outn = self.ports
                
        ## Add voltage source at input port and create
        ## copies with output open and shorted respectively
        circuit_vs_open = copy(self.c)

        circuit_vs_open['VS_TwoPort'] = VS(inp, inn, v=1.0)

        circuit_vs_shorted = copy(circuit_vs_open)

        circuit_vs_shorted['VL_TwoPort'] = VS(outp, outn, v=0.0)

        ## Run AC-analysis on the two circuits
        ac_open = self.ACAnalysis(circuit_vs_open)
        ac_shorted = self.ACAnalysis(circuit_vs_shorted)

        ac_open.run(freqs, refnode = refnode, complexfreq=complexfreq)

        ac_shorted.run(freqs, refnode = refnode, complexfreq=complexfreq)
        
        A = ac_open.v(inp, inn) / ac_open.v(outp, outn)
        B = ac_shorted.v(inp, inn) / ac_shorted.i('VL_TwoPort.plus')
        C = ac_open.i('VS_TwoPort.minus') / ac_open.v(outp, outn)
        D = ac_shorted.i('VS_TwoPort.minus') / ac_shorted.i('VL_TwoPort.plus')

        return array([[A,B],[C,D]])

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



