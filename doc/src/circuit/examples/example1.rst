Example 1
---------

Simple example - RC low-pass filter
```````````````````````````````````
.. image:: vdiv.*

.. sympy::

    import numpy, pylab
    from pycircuit.circuit import *

    circuit.default_toolkit = symbolic
    
    ## Create circuit
    cir = SubCircuit()
    cir['VS'] = VS(1, gnd, vac=1)
    cir['R1'] = R(1, 2, r=Symbol('R1'))
    cir['C1'] = C(2, gnd, c=Symbol('C1'))

    ## Run symbolic AC analysis
    ac = AC(cir, toolkit=symbolic)
    result = ac.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print transfer function from the voltage source to net 2
    simplify(result.v(2, gnd) / result.v(1, gnd))

Calculate ABCD parameters:

.. sympy::

    import numpy, pylab
    from pycircuit.circuit import *
    
    ## Create circuit
    cir = SubCircuit()
    ## n1,n2 = nodes('1','2')
    cir['R1'] = R(1, 2, r=Symbol('R1'))
    cir['C1'] = C(2, gnd, c=Symbol('C1'))

    ## Run symbolic 2-port analysis
    twoport_ana = TwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd, toolkit=symbolic)
    result = twoport_ana.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print ABCD parameter matrix
    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()
    ABCD

