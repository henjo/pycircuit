import pylab

def plotall(*waveforms):
    """Plot waveforms in a single plot"""
    
    pylab.hold(True)
    for wave in waveforms:
        wave.plot()
    pylab.legend([wave.ylabel for wave in waveforms])
    
    pylab.hold(False)
