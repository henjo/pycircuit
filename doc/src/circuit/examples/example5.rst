Example 5
---------

Noise analysis of cascaded resistor V->I and I->V amplifier
```````````````````````````````````````````````````````````
.. image:: viiv.*

In this example we are analyzing a balanced resistor V->I 
converter and a I->V amplifier in cascade. The active element of
the I->V amplifier is a noise nullor with imput referred noise
sources :math:`v_{n,1}` and :math:`i_{n,1}`.

The goal is to use the symbolic 2-port analysis to calculate gain and
input referred noise.

.. sympy::

    from sympy import *
    import numpy, pylab
    from pycircuit.circuit import *

    Ri,Rfb,Sv1,Si1 = symbols('Ri Rfb Sv1 Si1', real=True, positive=True)

    ## Create circuit
    cir = SubCircuit(toolkit=symbolic)
    cir['Ri1'] = R('vinp', 1, r=Ri, toolkit=symbolic)
    cir['Ri2'] = R('vinn', 2, r=Ri, toolkit=symbolic)
    cir['Rfb1'] = R(1, 'voutp', r=Rfb, toolkit=symbolic)
    cir['Rfb2'] = R(2, 'voutn', r=Rfb, toolkit=symbolic)
    cir['Vn'] = VS(1, 3, vac=0, noisePSD = Sv1, toolkit=symbolic)
    cir['In'] = IS(3, 2, iac=0, noisePSD = Si1, toolkit=symbolic)
    cir['nullor'] = Nullor(3,2,'voutp','voutn', toolkit=symbolic)

    ## Run symbolic 2-port analysis
    twoport_ana = TwoPortAnalysis(cir, 'vinp', 'vinn', 'voutp', 'voutn', method='sparam', noise=True)
    result = twoport_ana.solve(freqs=Symbol('s'), complexfreq=True, refnode=Node('vinn'))

    ## Print ABCD parameter matrix
    ABCD = Matrix(result['twoport'].A)
    ABCD.simplify()
    ABCD

    ## Calculate Voltage gain (mu)
    mu = 1 / ABCD[0,0]
    mu

    ## Input referred voltage noise power spectral density
    ## using |Rfb| = |mu| * Ri
    expand(result['Svn']).subs(Rfb, abs(mu)*Ri)

    ## Input referred current noise power spectral density
    ## using |Rfb| = |mu| * Ri
    expand(result['Sin']).subs(Rfb, abs(mu)*Ri)

    
