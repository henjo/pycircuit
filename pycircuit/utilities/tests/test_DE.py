# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

# Placed into the public domain by:
# James R. Phillips
# 2548 Vera Cruz Drive
# Birmingham, AL 35235 USA
# email: zunzun@zunzun.com

from pycircuit.utilities import DESolver
import numpy, time

class TestSolver(DESolver.DESolver):

    #some test data
    xData = numpy.array([5.357, 9.861, 5.457, 5.936, 6.161, 6.731])
    yData = numpy.array([0.376, 7.104, 0.489, 1.049, 1.327, 2.077])


    def externalEnergyFunction(self, trial):
        # inverse exponential with offset, y = a * exp(b/x) + c
        predicted = trial[0] * numpy.exp(trial[1] / self.xData) + trial[2]

        # sum of squared error
        error = predicted - self.yData
        return numpy.sum(error*error)


# for profiling and debugging
def Test500():
    for i in range(500):
        # parameterCount, populationSize, maxGenerations, minInitialValue, maxInitialValue, deStrategy, diffScale, crossoverProb, cutoffEnergy, useClassRandomNumberMethods, polishTheBestTrials
        solver = TestSolver(3, 600, 600, -10, 10, "Rand2Exp", 0.7, 0.6, 0.01, True, True)
        solver.Solve()



if __name__ == '__main__':

    # parameterCount, populationSize, maxGenerations, minInitialValue, maxInitialValue, deStrategy, diffScale, crossoverProb, cutoffEnergy, useClassRandomNumberMethods, polishTheBestTrials
    solver = TestSolver(3, 600, 600, -10, 10, "Rand2Exp", 0.7, 0.6, 0.01, True, True)
    tStart = time.time()
    solver.Solve()
    tElapsed = time.time() - tStart
    print "Best energy:", solver.bestEnergy, ": Elapsed time", tElapsed, 'seconds for', solver.generation, 'generation(s)'
    print tElapsed / solver.generation, 'seconds per generation with polisher on'
    print
    
    # parameterCount, populationSize, maxGenerations, minInitialValue, maxInitialValue, deStrategy, diffScale, crossoverProb, cutoffEnergy, useClassRandomNumberMethods, polishTheBestTrials
    solver = TestSolver(3, 600, 600, -10, 10, "Rand2Exp", 0.7, 0.6, 0.01, True, False)
    tStart = time.time()
    solver.Solve()
    tElapsed = time.time() - tStart
    print "Best energy:", solver.bestEnergy, ": Elapsed time", tElapsed, 'seconds for', solver.generation, 'generation(s)'
    print tElapsed / solver.generation, 'seconds per generation with polisher off'
    print

    #Test500()
    #print "Elapsed time for 500 tests", time.time()-tStart
