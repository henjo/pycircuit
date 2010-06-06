
#Circuit from netlist
circuit=Circuit('/here/input.scs')
circuit.elements
circuit.par.VAC=1 #sim.par would also exist and would override circuit.par
#Add element
circuit['VS']=VS(1,gnd,v=0,vac='VAC')
circuit['IDUT.Rout'].par.r=100e3
del circuit['IDUT.Rbias']
del circuit['IDUT']['Rbias'] #change subcircuit
circuit['IDUT.Rbias']=R(IDUT.1,IDUT.2,r=1e3) #change in in instance (not in subcircuit)
sim=Simulation(circuit)
sim['AC1']=AC(LogSweep(1e4, 1e5, decade=10)+LinSweep(1e4,3e4,10)) #+ indicates union of all sweep points
sim['AC1'].options.reltol=1e-3
result=sim.run() #run all analyses
result['AC1'].v("IDUT.vOut")
result['AC1'].v("IDUT.bOut[5:1]") #bus slice
for v in result['AC1'].v("*Out"): #voltage of all nodes ending with Out (consider regexp)
    v.plot()


#Build circuit in pycircuit
circuit=Circuit()
circuit['VS']=VS(1,gnd,v=0,vac=1)
circuit['R']=R(2,gnd,r=100)
circuit['VCVS']=VCVS(1,gnd,2,gnd,gm=1)
sim=Simulation(circuit)
sim.options.reltol=1e-3
sim['AC1']=AC(LogSweep(1e4, 1e5, decade=10)+LinSweep(1e4,3e4,10))
sim['AC2'] = Alter(SetVariables(v2=2), SetOptions(reltol=1e-4), AC(LogSweep(1e4, 1e5, decade=10)+LinSweep(1e4,3e4,10))) #?

#Problems in batch
result=sim.run()
result['AC1'].v(2)


#
circuit=Circuit()
circuit['VS']=VS(1,gnd,v=0,vac=1)
circuit['R']=R(2,gnd,r=100)
circuit['VCVS']=VCVS(1,gnd,2,gnd,gm=1)
sim1=Simulation(circuit)
sim1.options.reltol=1e-3
sim1['AC1']=AC(LogSweep(1e4, 1e5, decade=10))
sim1['DC1']=DC()
sim2=sim1.copy() #deep copy? ...
sim2.options.reltol=1e-4
sim2.par.v2=2
sim=SimGroup(sim1,sim2) #possible to use different circuits in sim1, sim2
sim.run()

#sweep
