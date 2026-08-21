"""Regression tests for the Phase-0 fixes of doc/transient_review_260820.md.

One test (or small group) per finding, appended in the order the fixes landed
so each commit carries the test that fails before it and passes after it.
"""

import numpy as np
import pytest

from pycircuit.circuit import gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, IS
from pycircuit.circuit.transient import Transient


def _rc():
    c = SubCircuit()
    c['is'] = IS(gnd, 'a', i=1e-3)
    c['R'] = R('a', gnd, r=1e3)
    c['C'] = C('a', gnd, c=1e-9)
    return c


## F8 / F18 -- dead arguments deleted; passing them must now fail loudly
## rather than being silently discarded.

def test_f8_analytical_eh_is_gone():
    with pytest.raises(TypeError):
        Transient(_rc()).solve(tend=1e-6, timestep=1e-8, analytical_eh=True)


def test_f18_init_irefnode_is_gone():
    ## Not TypeError: the deleted argument now falls into **kvargs, where the
    ## parameter dictionary rejects it by name -- louder and more specific.
    with pytest.raises(KeyError, match='irefnode'):
        Transient(_rc(), irefnode=0)


def test_f18_get_diff_method_is_gone():
    tran = Transient(_rc())
    with pytest.raises(TypeError):
        tran.get_diff(np.zeros(2), np.zeros((2, 2)), method='euler')
