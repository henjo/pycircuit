Example 2
----------

Simple example - Multi feedBack (MFB) filter 
```````````````````````````````````````````````
.. image:: mfb.*

Find symbolic expression of transfer function from Is to V(3,0):

.. sympy::
     :persistent:


    from pycircuit.circuit import *
    from sympy import symbols, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect    

    ## Create circuit
    cir = SubCircuit()

    R1,R2,R3,C1,C2,i_s = symbols('R1 R2 R3 C1 C2 is', real=True)
    s = Symbol('s', complex = True)   

    cir['R3'] = R(1, gnd, r = R3)
    cir['R2'] = R(1, 2, r = R2)
    cir['R1'] = R(1, 3, r = R1)
    cir['C1'] = C(1, gnd, c = C1)
    cir['C2'] = C(2, 3, c = C2)
    cir['Nullor'] = Nullor(2, gnd, 3, gnd)

    # Current source for AC stimuli
    cir['ISource'] = IS(1, gnd, iac=i_s)

    ## Run symbolic AC analysis     
    ac = SymbolicAC(cir)
    result = ac.solve(freqs=s, complexfreq=True)

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

   solve(denominator, s)[0]

.. sympy::
   :persistent:

   solve(denominator, s)[1]
   
