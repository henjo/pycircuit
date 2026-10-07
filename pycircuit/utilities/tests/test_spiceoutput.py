"""`spiceoutput`: the simulators' output readers and Xyce's regression
metric, transcribed (the SPICE benchmark plan's stage 4) -- on synthetic
files and hand-computed values; the 4049 oscillator against Xyce's gold
where the benchmark data is fetched."""
import textwrap

import numpy as np
import pytest

from pycircuit.utilities import spiceoutput as so


def test_a_prn_its_header_rows_and_end(tmp_path):
    p = tmp_path / 'a.prn'
    p.write_text(textwrap.dedent("""\
        Index       TIME        {V(8)+4}        {I(VPULSE) + 3.0E-3}
        0   0.0e+00   4.0e+00   1.5e+00
        1   1.0e-09   4.5e+00  -2.0e-01

        End of Xyce(TM) Simulation
        """))
    names, a = so.read_prn(str(p))
    assert names == ['Index', 'TIME', '{V(8)+4}', '{I(VPULSE) + 3.0E-3}']
    assert a.tolist() == [[0.0, 0.0, 4.0, 1.5], [1.0, 1e-9, 4.5, -0.2]]


def test_ngspice_print_tables_split_by_width_and_repeated_per_page():
    text = textwrap.dedent("""\
        Circuit: test
        Index   time            v(a)            v(b)
        ---------------------------------------------
        0       0.000000e+00    1.000000e+00    2.000000e+00
        1       1.000000e-09    1.500000e+00    2.500000e+00

        Index   time            v(a)            v(b)
        2       2.000000e-09    1.750000e+00    2.750000e+00

        Index   time            v(c)
        0       0.000000e+00    3.000000e+00
        1       1.000000e-09    3.500000e+00
        2       2.000000e-09    3.750000e+00
        """)
    cols, a = so.read_ngspice_print(text)
    assert cols == ['time', 'v(a)', 'v(b)', 'v(c)']
    assert a.tolist() == [[0.0, 1.0, 2.0, 3.0], [1e-9, 1.5, 2.5, 3.5], [2e-9, 1.75, 2.75, 3.75]]


def test_the_ibm_reference_output():
    got = so.read_ibm_output(['', 'Node: n0_1_2', '', ' 0.000e+00 3.5e-04', ' 1.000e-11 3.6e-04',
                              'Node: n9', ' 0 1.8'])
    assert list(got) == ['n0_1_2', 'n9']
    assert got['n0_1_2'][0].tolist() == [0.0, 1e-11] and got['n9'][1].tolist() == [1.8]


def test_comp_lines_set_a_columns_tolerances():
    assert so.comp_tolerances(['{v(5)+1.0} reltol=0.025', 'V(8) RELTOL=0.02 offset=1',
                               'I(V1) numfail=3']) == {
        '{v(5)+1.0}': {'reltol': 0.025}, 'V(8)': {'reltol': 0.02, 'offset': 1.0}, 'I(V1)': {}}


def test_the_gold_is_interpolated_at_the_test_times_as_the_script_does():
    """The first gold time at or after each test time: its value where the
    times are equal (the first of a repeated time), else the line from the
    one before; a test outside the gold's span is refused."""
    gt = np.array([0.0, 1.0, 1.0, 3.0])
    g = np.array([0.0, 10.0, 20.0, 40.0])
    assert so.interpolate_at([0.0, 0.5, 1.0, 2.0, 3.0], gt, g).tolist() == [
        0.0, 5.0, 10.0, 30.0, 40.0]
    with pytest.raises(ValueError, match='does not span'):
        so.interpolate_at([-1.0, 1.0], gt, g)
    with pytest.raises(ValueError, match='does not span'):
        so.interpolate_at([0.0, 4.0], gt, g)


def test_the_metric_is_the_rms_relative_error_in_units_of_reltol():
    """Equal series: 0.  A test 1 % below the gold everywhere: 1 in units of
    reltol 0.01 (the edge of passing), 0.5 at reltol 0.02.  A hand-computed
    trapezoid: integrand 0, 1, 0 at t = 0, 1, 2 -> sqrt(1/2).  absdifftol
    zeroes a small difference; zerotol and offset act on both series."""
    t = np.linspace(0.0, 1.0, 11)
    g = 1.0 + t
    assert so.xyce_verify(t, g, t, g) == 0.0
    assert so.xyce_verify(t, g * 0.99, t, g, abstol=0.0) == pytest.approx(1.0, rel=1e-12)
    assert so.xyce_verify(t, g * 0.99, t, g, abstol=0.0, reltol=0.02) == pytest.approx(0.5,
                                                                                     rel=1e-12)
    tt = np.array([0.0, 1.0, 2.0])
    assert so.xyce_verify(tt, [1.0, 0.99, 1.0], tt, [1.0, 1.0, 1.0], abstol=0.0) == \
        pytest.approx(np.sqrt(0.5), rel=1e-12)
    assert so.xyce_verify(tt, [1.0, 1.0 - 1e-13, 1.0], tt, [1.0] * 3) == 0.0
    assert so.xyce_verify(tt, [1e-13, 0.0, 0.0], tt, [0.0, 0.0, 0.0]) == 0.0
    assert so.xyce_verify(tt, [0.5, 0.5, 0.5], tt, [0.5, 0.5, 0.5], offset=-0.5) == 0.0
    assert so.xyce_verify([0.0], [0.98], [0.0], [1.0], abstol=0.0) == pytest.approx(2.0)
    with pytest.raises(TypeError, match='unknown tolerances'):
        so.xyce_verify(tt, g[:3], tt, g[:3], numfail=1)


def test_crossing_times():
    t = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    v = np.array([0.0, 2.0, 0.0, 2.0, 2.0])
    assert so.crossing_times(t, v, 1.0, rising=True).tolist() == [0.5, 2.5]
    assert so.crossing_times(t, v, 1.0, rising=False).tolist() == [1.5]
    assert so.crossing_times(t, v, 1.0).tolist() == [0.5, 1.5, 2.5]
