********************
Simulator interfaces
********************

External simulators (gnucap, Spectre) are run by `simdeck <https://github.com/henjo/simdeck>`_,
whose results are polars-waveform waveforms like pycircuit's own numeric results::

    from simdeck.gnucap import Gnucap

    ac = Gnucap(netlist).ac(1e3, 1e9, decade=20)
    ac.v("out").bandwidth()

Spectre results on disk are read by :mod:`pycircuit.post.cds` (polars-psf).
