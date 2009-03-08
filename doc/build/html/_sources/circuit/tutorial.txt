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
