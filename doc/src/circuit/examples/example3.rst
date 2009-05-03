Example 3
----------

Simple example - VCCS with resistor load 
```````````````````````````````````````````````
.. image:: vccs_res.*

Find symbolic expression of transfer function from input voltage to output voltage:

.. sympy::
     :persistent:

    import numpy, pylab
    from pycircuit.circuit import *
    from pycircuit.post.functions import *
    from sympy import Symbol, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect    

    pycircuit.circuit.circuit.default_toolkit = symbolic

    ## Create circuit
    cir = SubCircuit()

    n1, n2 = cir.add_nodes('1', '2')
    R1, gm = [Symbol(symname, real=True) for symname in 'R1,gm'.split(',')]

    cir['VS'] = VS(n1, gnd, v=1)
    cir['R1'] = R(n2, gnd, r = R1)
    cir['VCCS'] = VCCS(n1, gnd, n2, gnd,gm=gm)

    ## Run symbolic AC analysis     
    ac = AC(cir, toolkit=symbolic)
    result = ac.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print transfer function from the voltage source to output
    simplify(result.v(2, gnd))
