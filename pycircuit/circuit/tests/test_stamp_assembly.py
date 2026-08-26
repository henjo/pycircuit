# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The stamping loop stores VIEWS of element stamps, not copies.

`SubCircuit._add_element_submatrices` and `_add_element_subvectors`
collect each element's contribution into a pending list and scatter the
lot in one pass.  They used to do that with ``np.asarray(rhs).flatten()``
-- and ``flatten`` *always* copies, even when the stamp is already a
contiguous array, which a generated 2-D stamp always is.  The copy was
pure waste (412 ns against 100 ns per element) because the only consumer
is ``np.concatenate``, which allocates its own buffer regardless.

Changed to ``ravel``, which returns a view when it can.  That makes the
pending list hold references into element-owned memory for the duration
of one assembly call, so the property this file pins is that **nothing
observable aliases**:

* the assembled matrix must not share memory with any element's stamp,
  so writing to the result cannot corrupt an element (`test_..._is_not_a
  _view`) -- the failure this would cause is silent and arbitrarily far
  from its cause;
* an element serving a CACHED constant stamp (`_hdl_Gc`, the same array
  object returned every call) must survive being stamped repeatedly;
* the assembled answer must be unchanged from the copying version, which
  is asserted directly against a hand-computed reference rather than
  against a golden number.
"""

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd, elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit, defaultepar

cm.default_toolkit = numeric


def _ladder(n=4):
    c = SubCircuit()
    c['vs'] = eh.VSinHdl('n0', gnd, va=5.0, freq=1e3)
    for k in range(n):
        c['D%d' % k] = eh.DiodeHdl('n%d' % k, 'n%d' % (k + 1))
        c['R%d' % k] = eh.RHdl('n%d' % (k + 1), gnd, r=1e3)
        c['C%d' % k] = eh.CHdl('n%d' % (k + 1), gnd, c=1e-7)
    c.update_iparv()
    return c


def _element_stamp_buffers(c, x, meth):
    """Every element's raw stamp, as the assembly loop sees it."""
    out = []
    for instance, element in c.elements.items():
        sub = x[c.elementnodemap[instance]]
        out.append(np.asarray(getattr(element, meth)(sub, defaultepar)))
    return out


@pytest.mark.parametrize('meth', ['G', 'C'])
def test_the_assembled_matrix_is_not_a_view_of_any_element(meth):
    c = _ladder()
    x = np.zeros(c.n)
    lhs = np.asarray(getattr(c, meth)(x, defaultepar))
    for buf in _element_stamp_buffers(c, x, meth):
        assert not np.shares_memory(lhs, buf), (
            'the assembled %s aliases an element stamp' % meth)


@pytest.mark.parametrize('meth', ['i', 'q'])
def test_the_assembled_vector_is_not_a_view_of_any_element(meth):
    c = _ladder()
    x = np.zeros(c.n)
    rhs = np.asarray(getattr(c, meth)(x, defaultepar))
    for buf in _element_stamp_buffers(c, x, meth):
        assert not np.shares_memory(rhs, buf), (
            'the assembled %s aliases an element stamp' % meth)


def test_writing_to_the_result_cannot_corrupt_an_element():
    ## The concrete damage an aliasing bug would do, and the reason the
    ## test above is worth having: a caller that scales or clears the
    ## assembled matrix would silently rewrite the element's own stamp,
    ## and the next iteration would use the corrupted one.
    c = _ladder()
    x = np.zeros(c.n)
    before = np.asarray(c.elements['R0'].G(x[c.elementnodemap['R0']],
                                           defaultepar)).copy()
    lhs = np.asarray(c.G(x, defaultepar))
    lhs *= -1234.0
    after = np.asarray(c.elements['R0'].G(x[c.elementnodemap['R0']],
                                          defaultepar))
    assert np.array_equal(before, after), 'element stamp was corrupted'


def test_a_cached_constant_stamp_survives_repeated_assembly():
    ## A linear element serves the SAME array object every call
    ## (`_hdl_Gc`).  Under `ravel` that object is what lands in the
    ## pending list, so stamping repeatedly is the case most likely to
    ## expose an aliasing mistake.
    c = _ladder()
    x = np.zeros(c.n)
    first = np.asarray(c.G(x, defaultepar)).copy()
    for _ in range(5):
        again = np.asarray(c.G(x, defaultepar))
        assert np.array_equal(first, again), 'repeated assembly drifted'


def test_the_assembled_conductance_matches_a_hand_sum():
    ## Not a golden number: the assembly is checked against the sum it is
    ## supposed to be computing, built independently from each element's
    ## own stamp and index map.
    c = _ladder()
    x = np.zeros(c.n)
    lhs = np.asarray(c.G(x, defaultepar))

    ref = np.zeros((c.n, c.n))
    for instance, element in c.elements.items():
        rc = c._map_indices_2d.get(instance)
        if rc is None:
            continue
        rows, cols = rc
        sub = x[c.elementnodemap[instance]]
        val = np.asarray(element.G(sub, defaultepar)).ravel()
        np.add.at(ref, (rows, cols), val)

    assert np.array_equal(lhs, ref), 'assembled G differs from the hand sum'
