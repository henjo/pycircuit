import numpy as np
import polars_waveform as pw
import pytest
import sympy

from pycircuit.circuit import AC, C, R, SubCircuit, VS, gnd, symbolic
from pycircuit.post import PandasWaveform, from_arrays
from pycircuit.post.functions import bandwidth, db20

Rs, Cs = sympy.symbols('R C', positive=True)
FREQS = np.logspace(3, 7, 81)


def rc(r, c):
    cir = SubCircuit()
    cir['VS'] = VS('1', gnd, vac=1.0)
    cir['R1'] = R('1', '2', r=r)
    cir['C1'] = C('2', gnd, c=c)
    return cir


def test_symbolic_sweep_is_a_pandas_waveform():
    v2 = AC(rc(Rs, Cs), toolkit=symbolic).solve(FREQS).v('2')
    assert isinstance(v2, PandasWaveform) and v2.xname == 'frequency' and v2.yunit == 'V'
    assert db20(v2).to_pandas().iloc[0].has(Rs)                      # elementwise: still symbolic
    with pytest.raises(ValueError, match='C'):
        v2.numeric()
    num = v2.subs({Rs: 1e3, Cs: 1e-9})
    assert bandwidth(num) == pytest.approx(1 / (2 * np.pi * 1e-6), rel=1e-2)
    ref = AC(rc(1e3, 1e-9)).solve(FREQS).v('2')                     # the numeric analysis agrees
    assert isinstance(ref, pw.Waveform)
    np.testing.assert_allclose(num.numeric().to_numpy(), ref.to_numpy(), rtol=1e-9)


def test_from_arrays_kinds():
    assert isinstance(from_arrays(FREQS, np.ones(81)), pw.Waveform)
    assert isinstance(from_arrays(FREQS, [Rs] * 81), PandasWaveform)
