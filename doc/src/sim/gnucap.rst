Gnucap interface
================

Simple example
--------------

This example shows a transient analysis of a CMOS inverter where
the circuit is given as a native gnucap netlist.

.. plot::
     :include-source: True
      
     import numpy as np
     import pylab
     import pycircuit.sim.gnucap as gnucap
     
     netlist = """
               .MODEL CMOSP PMOS (level=2 KP=28U VTO=0.4 LAMBDA=0.01 GAMMA=0.9 PHI=0.5)
	       .MODEL CMOSN NMOS (level=2 KP=28U VTO=0.4 LAMBDA=0.01 GAMMA=0.9 PHI=0.5)
               .subckt inv psup nsup out in 
               Mp out in psup psup CMOSP l=0.25u w=1.0u
	       Mn out in nsup nsup CMOSN l=0.25u w=1.0u
	       .ends inv
	  
	       vsup psup 0 1.8     
	       x1 psup 0 out in inv
	       
	       vin in 0 pwl (0 0, 10n 0, 10.1n 1.8, 20n 1.8, 20.1n 0, 30n, 0, 30.1n 1.8)
	       """


     sim = gnucap.Simulation(netlist, direct=True)
     res = sim.run_transient(t=np.linspace(0, 40e-9, 100))

     pylab.hold()
     res.v('in').plot()
     res.v('out').plot()
     pylab.legend()



