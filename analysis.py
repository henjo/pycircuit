from numpy import array, delete, linalg, size, zeros, concatenate, pi, zeros
from circuit import Circuit, SubCircuit, VS,R,C, gnd
from result import Waveform
from internalresult import InternalResultSet, InternalResult

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



class DC(Analysis):
    """DC analyis class
    
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

    """
    
    def solve(self, refnode=gnd):
        n=self.c.n()
        G=self.c.G(zeros((n,1)))
        U=self.c.U()

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.c.getNodeIndex(refnode)
        G,U=removeRowCol((G,U), irefnode)
        # Calculate state vector
        x = linalg.solve(G,-U)

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

        G=self.c.G(zeros((n,1)))
        C=self.c.C(zeros((n,1)))
        U=self.c.U()

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
    import doctest
    doctest.testmod()
