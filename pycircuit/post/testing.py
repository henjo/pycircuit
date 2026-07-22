import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

def assert_waveform_equal(a, b):
    for ax, bx in zip(a.x, b.x):
        assert_array_equal(ax, bx)

    assert_array_equal(np.array(a.y), np.array(b.y))
#    assert a.xlabels == b.xlabels
#    assert a.ylabel == b.ylabel
#    assert a.xunits == b.xunits
#    assert a.yunit == b.yunit

def assert_waveform_almost_equal(a, b, decimal=6):
    for ax, bx in zip(a.x, b.x):
        assert_array_almost_equal(ax, bx, decimal=decimal)
    assert_array_almost_equal(np.array(a.y), np.array(b.y), decimal=decimal)
#    assert a.xlabels == b.xlabels
#    assert a.ylabel == b.ylabel
#    assert a.xunits == b.xunits
#    assert a.yunit == b.yunit

