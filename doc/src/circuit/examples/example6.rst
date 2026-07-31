Example 6
---------

Modulus integrator (Idtmod)
```````````````````````````

This example shows the modulus integrator in a transient simulation.
The output voltage is the time integral of the input voltage, modulo the modulus.
With a constant input of 1V and a modulus of 1.0, the output will ramp up to 1.0, wrap around to 0.0, and ramp up again.

.. plot::
    :include-source: True
    :width: 10cm

    import pylab
    import pycircuit.circuit._numeric as numeric
    from pycircuit.circuit.transient import Transient
    from pycircuit.circuit import SubCircuit, VS, R, Idtmod, gnd

    c = SubCircuit()
    nin, nout = c.add_nodes('in', 'out')
    c['vin'] = VS(nin, gnd, v=1.0)
    c['R'] = R(nout, gnd, r=1e3)
    c['Idtmod'] = Idtmod(nin, gnd, nout, gnd, modulus=1.0)

    ## An ideal integrator has no DC operating point -- there is no resistive path
    ## setting the output -- so the bias solve is skipped deliberately.  `uic` is a
    ## Transient() argument, not a solve() one.
    tran = Transient(c, toolkit=numeric, uic=True)
    result = tran.solve(tend=3.0, timestep=0.01)

    t = result.v(nout).x[0]
    vin = result.v(nin).y
    vout = result.v(nout).y

    pylab.plot(t, vin, label='Input Voltage')
    pylab.plot(t, vout, label='Output Voltage (mod 1.0)')
    pylab.xlabel('Time (s)')
    pylab.ylabel('Voltage (V)')
    pylab.title('Idtmod Transient Response')
    pylab.legend()
    pylab.grid(True)
