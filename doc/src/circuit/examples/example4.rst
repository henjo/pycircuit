Example 4
----------

Simple example - Input impedance of gyrator loaded with capacitance
```````````````````````````````````````````````````````````````````
###.. image:: gyrator.*

Find symbolic expression of input impedance of gyrator loaded with gounded capacitor:

.. sympy::
    :include-source: False
    :persistent:

    import numpy, pylab
    from pycircuit.circuit import *
    from pycircuit.post.functions import *
    from sympy import Symbol, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect

    circuit.default_toolkit = symbolic

    C1, gm1 = [Symbol(symname, real=True , positive=True ) for symname in 'C1,gm1'.split(',')]
    s = Symbol('s', complex = True)

    ## Create circuit object
    cir = SubCircuit()

    ## Add nodes to circuit
    n1, n2 = cir.add_nodes('1', '2')

    ## Add circuit elements
    cir['cap']  = C(n2, gnd, c = C1)
    # Gyrator
    cir['Gyrator'] = Gyrator(n1, gnd, n2, gnd, gm = 1)
    # Current source for AC stimuli
    cir['ISource'] = IS(gnd,n1, iac=1)

    ## Run symbolic AC analysis
    ac = AC(cir, toolkit=symbolic)    
    result = ac.solve(freqs=s, complexfreq=True)

    # Input impedance
    impedance = simplify(result.v(n1, gnd))
    impedance

ABCD matrix:

.. sympy::
    :persistent:

    # Delet current source used to calculate input impedance 
    
    del cir['ISource']

    ## Run symbolic 2-port analysis
    twoport_ana = TwoPortAnalysis(cir, Node('1'), gnd, Node('2'), gnd, toolkit=symbolic)
    result = twoport_ana.solve(freqs=s, complexfreq=True)

    ## Print ABCD parameter matrix
    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()
    ABCD
