def plotall(*waveforms, **args):
    """Plot waveforms in a single plot"""
    import matplotlib.pyplot as plt

    plotkvargs = dict(args.get('plotkvargs', {}))
    ax = plotkvargs.pop('ax', None) or plt.gca()

    for wave in waveforms:
        wave = wave.numeric()
        wave.plot(ax=ax, **({'label': wave.yname} | plotkvargs))

    ax.legend()
    return ax
