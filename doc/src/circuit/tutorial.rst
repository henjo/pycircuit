Tutorial
--------

Simple example - RC low-pass filter
```````````````````````````````````
.. image:: examples/vdiv.*

Find transfer function from V1 to V(2,0):

.. plot::
    :include-source: True
    :width: 10cm

    import numpy, pylab
    from pycircuit.circuit import *
    from pycircuit.post.functions import *
    
    ## Create circuit
    cir = SubCircuit()
    cir['VS'] = VS(1, gnd, vac=1.0)
    cir['R1'] = R(1, 2, r=1e3)
    cir['C1'] = C(2, gnd, c=1e-12)

    ## Run AC analysis
    ac = AC(cir)
    result = ac.solve(freqs=numpy.logspace(6,9))

    ## Plot voltage between net 2 and ground    
    v2 = db20(result.v(2, gnd))
    v2.semilogx()
    pylab.grid(True)

And now symbolically using a symbolic ac analysis:

.. sympy::

    import numpy, pylab
    from pycircuit.circuit import *
    
    ## Create circuit
    cir = SubCircuit()
    cir['VS'] = VS(1, gnd, vac=1)
    cir['R1'] = R(1, 2, r=Symbol('R1'))
    cir['C1'] = C(2, gnd, c=Symbol('C1'))

    ## Run symbolic AC analysis
    ac = SymbolicAC(cir)
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
    twoport_ana = SymbolicTwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd)
    result = twoport_ana.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print ABCD parameter matrix
    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()
    ABCD
