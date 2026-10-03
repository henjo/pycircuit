Waveforms
=========

Numeric results (``res.v('out')`` of an AC, DC sweep or transient analysis) are
`polars-waveform <https://github.com/henjo/polars-waveform>`_ waveforms: lazy measurements,
families over parameter sweeps and plotting. The :class:`Waveform` below holds symbolic results
(sympy expressions); its measurements evaluate on ``numeric()``, the polars-waveform equivalent.

Classes
-------

.. module:: pycircuit.post
.. autoclass:: Waveform
   :members: 

Functions
---------

.. autofunction:: astable
.. autofunction:: compatible
.. autofunction:: compose
