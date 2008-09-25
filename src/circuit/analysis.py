from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros, alltrue, maximum
from scipy import optimize
from circuit import Circuit, SubCircuit, VS,IS,R,C,Diode, gnd
from pycircuit.result import Waveform
from pycircuit.internalresult import InternalResultSet, InternalResult
import numpy
import copy

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
        n=self.c.n()
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
        s_scale = copy.copy(x0)
        
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
    """
    
    def run(self, freqs, **kvargs):
        x = self.solve(freqs, **kvargs)
        
        result = InternalResult()

        nodes = self.c.nodes
        for xvalue, node in zip(x[:len(nodes)], nodes):
            wave = Waveform(freqs, xvalue)
            result.storeSignal(self.c.getNodeName(node), wave)
        for i, data in enumerate(zip(x[len(nodes):], self.c.branches)):
            xvalue, branch = data
            wave = Waveform(freqs, xvalue)
            result.storeSignal('i' + str(i),wave)

        self.result = result

        return result

    def solve(self, freqs, refnode=gnd):
        n=self.c.n()
        x = zeros((n,1))
        G=self.c.G(x)
        C=self.c.C(x)
        U=self.c.U(x)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,C,U = removeRowCol((G,C,U), irefnode)
        G,C,U = (A.astype(complex) for A in (G,C,U))

        out = []
        for f in freqs:
            x = linalg.solve(2j*pi*f*C + G, -U) 
            
            # Insert reference node voltage
            x = concatenate((x[:irefnode, :], array([[0.0]]), x[irefnode:,:]))

            out.append(x)

        # Swap frequency and x-vector dimensions

        return array(out)[:,:,0].swapaxes(0,1)

if __name__ == "__main__":
    c = SubCircuit()
    n1 = c.addNode('net1')
    c['is'] = IS(gnd, n1, i=100e-6)
    c['D'] = Diode(n1, gnd)
    c['D'].G(array([[0,0]]).T)
    dc = DC(c)
    res = dc.run()


    import doctest
    doctest.testmod()
