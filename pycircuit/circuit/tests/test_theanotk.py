from nose.tools import *
from ..circuit import *
from ..elements import *
from ..na import MNA
from .. import theanotk

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
from numpy.testing.decorators import slow

def test_eval_iqu_and_der():
    r = R(0,1,r=10e3)

    theanotk.generate_eval_iqu_and_der(r)
    
    iqu, GC = r.eval_iqu_and_der([0.1], defaultepar)

    assert_equal(iqu, (0.1/10e3,))

    assert_equal(GC, (1/10e3,))




