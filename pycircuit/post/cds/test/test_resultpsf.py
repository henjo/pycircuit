import os

import numpy as np
import polars_waveform as pw
import pytest

from pycircuit.post.cds import PSFResultSet

here = os.path.dirname(__file__)


def test_operating_point_and_info():
    rs = PSFResultSet(os.path.join(here, 'dcop.raw'))
    assert 'dcOp-dc' in rs.keys() and 'designParamVals-info' in rs
    assert rs['dcOp-dc']['vout'] == pytest.approx(2.5)
    assert rs['dcOp']['vin'] == pytest.approx(5.0)
    assert rs['designParamVals-info']['top-level']['k'] == pytest.approx(0.5)
    with pytest.raises(KeyError):
        rs['nope']


def test_sweep_and_family():
    out = PSFResultSet(os.path.join(here, 'dcsweep.raw'))['dc1-dc']['out']
    assert isinstance(out, pw.Waveform)
    np.testing.assert_allclose(out.x.to_numpy(), [1.0, 3.66666667, 6.33333333, 9.0])
    np.testing.assert_allclose(out.y.to_numpy(), [2.0, 4.66666667, 7.33333333, 10.0])

    fam = PSFResultSet(os.path.join(here, 'pardcsweep.raw'))['dc1'].v('out')
    assert fam.groups == ['vdc3', 'vdc2'] and fam.ymax().height == 16
