Example 2
----------

Simple example - Multi feedBack (MFB) filter 
```````````````````````````````````````````````
.. image:: mfb.*

Find symbolic expression of transfer function from Is to V(3,0):

.. sympy::
     :persistent:


    import numpy, pylab
    from pycircuit.circuit import *
    from pycircuit.post.functions import *
    from sympy import Symbol, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect    

    ## Create circuit
    cir = SubCircuit()

    n1, n2, n3 = cir.add_nodes('1', '2', '3')
    R1, R2, R3, C1, C2, i_s = [Symbol(symname, real=True) for symname in 'R1,R2,R3,C1,C2,is'.split(',')]

    cir['R3'] = R(n1, gnd, r = R3)
    cir['R2'] = R(n1, n2, r = R2)
    cir['R1'] = R(n1, n3, r = R1)
    cir['C1'] = C(n1, gnd, c = C1)
    cir['C2'] = C(n2, n3, c = C2)
    cir['Nullor'] = Nullor(n2, gnd, n3, gnd)

    # Current source for AC stimuli
    cir['ISource'] = IS(n1, gnd, iac=i_s)

    ## Run symbolic AC analysis     
    ac = SymbolicAC(cir)
    result = ac.solve(freqs=Symbol('s'), complexfreq=True)

    ## Print transfer function from the voltage source to net 2
    simplify(result.v(3, gnd) / i_s)

Get denominator of transfer function:

.. sympy::
   :persistent:
   
   denominator = fraction(simplify(result.v(3, gnd) / i_s))[1]
   denominator

Find poles of transfer function:

.. sympy::
   :persistent:

   s = Symbol('s', complex = True)   
   solve(denominator.subs(s,s),s)[0]

.. sympy::
   :persistent:

   s = Symbol('s', complex = True)   
   solve(denominator.subs(s,s),s)[1]
   
