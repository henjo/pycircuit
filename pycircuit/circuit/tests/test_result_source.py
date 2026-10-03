import numpy as np
import polars_waveform as pw
import pytest

from pycircuit.circuit import AC, C, R, SubCircuit, VS, gnd, symbolic


def rc():
    c = SubCircuit()
    c['VS'] = VS('1', gnd, vac=1.0)
    c['R1'] = R('1', '2', r=1e3)
    c['C1'] = C('2', gnd, c=1e-9)
    return c


def test_numeric_result_is_a_result_source():
    res = AC(rc()).solve(np.logspace(3, 7, 201))
    assert isinstance(res, pw.ResultSource) and res.names == ['1', '2']
    df = res.scan().collect()
    assert df.columns == ['frequency', 'v(1)', 'v(2)'] and df.height == 201
    w = res.v('2')
    assert isinstance(w, pw.Waveform) and w.yunit == 'V' and w.xunit == 'Hz'
    assert w.bandwidth() == pytest.approx(1 / (2 * np.pi * 1e-6), rel=1e-2)


def test_symbolic_result_stays_symbolic():
    import sympy
    s = sympy.Symbol('s')
    res = AC(rc(), toolkit=symbolic).solve(s, complexfreq=True)
    assert not isinstance(res.v('2'), pw.Waveform)
