def plotall(*waveforms, **args):
    """Plot waveforms in a single plot"""
    import pylab

    plotargs = []
    plotkvargs = {}

    if 'plotargs' in args:
        plotargs = args['plotargs']
    if 'plotkvargs' in args:
        plotkvargs = args['plotkvargs']

    pylab.hold(True)
    for wave in waveforms:
        wave.plot(*plotargs, **plotkvargs)

    pylab.legend([wave.ylabel for wave in waveforms])
    
    pylab.hold(False)
