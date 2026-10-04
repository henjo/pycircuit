"""Results and post-processing.

Waveforms come from polars-waveform: numeric results are ``polars_waveform.Waveform`` (lazy
measurements, families, plotting), symbolic ones ``polars_waveform.PandasWaveform`` (sympy
values; ``subs()`` then ``numeric()`` for measurements). ``pycircuit.post.functions`` re-exports
the calculator functions.
"""

from polars_waveform import PandasWaveform, Waveform, WaveformBase, from_arrays

from .result import *
from .internalresult import *
from .functions import *
from .plot import *
