# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""Functions that operate on waveforms or scalars.

These are the calculator functions of :mod:`polars_waveform.functions`. They work on numeric
results (polars_waveform.Waveform), symbolic ones (polars_waveform.PandasWaveform: elementwise
functions stay symbolic, measurements evaluate numerically) and plain numbers and numpy arrays.

Changes from earlier pycircuit versions: ``cross(w, threshold, edge=1, type="either")`` counts
crossings from 1 (negative from the end) and the edge type is ``"rising"``, ``"falling"`` or
``"either"``.
"""

from polars_waveform.functions import *  # noqa: F401,F403
from polars_waveform.functions import __all__  # noqa: F401
