Waveforms
=========

Swept results (``res.v('out')`` of an AC, DC sweep or transient analysis) are
`polars-waveform <https://github.com/henjo/polars-waveform>`_ waveforms:

* numeric results: ``polars_waveform.Waveform``, with lazy measurements, families over
  parameter sweeps and plotting;
* symbolic results (sympy values) and complex-frequency sweeps:
  ``polars_waveform.PandasWaveform``. Elementwise math stays symbolic; ``subs()`` substitutes
  values and the measurements run on ``numeric()``.

.. code-block:: python

    import sympy
    from pycircuit.circuit import AC, symbolic

    R, C = sympy.symbols('R C', positive=True)
    v = AC(cir, toolkit=symbolic).solve(freqs).v('out')   # PandasWaveform
    v.db20()                                             # symbolic
    v.subs({R: 1e3, C: 1e-9}).bandwidth()                # number

``pycircuit.post.from_arrays(x, y, ...)`` builds either kind from arrays.
