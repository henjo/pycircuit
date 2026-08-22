# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""``_pnjlim`` against Listing 1 of the PCNR paper.

Aadithya, Keiter & Mei, *"Predictor/Corrector Newton-Raphson (PCNR)"*
(doi 10.1007/978-3-030-44101-2_19), Listing 1::

    def pnjlim(vold, vnew, Is, VT):
        vc = VT*math.log(VT/(Is*math.sqrt(2.0)))
        if vnew <= vc or abs(vnew - vold) <= 2.0*VT:
            ans, limiting_applied = vnew, False
        else:
            limiting_applied = True
            if vold > 0.0:
                arg = 1.0 + (vnew - vold)/VT
                if arg > 0.0:
                    ans = (vold + VT*math.log(arg))
                else:
                    ans = vc
            else:
                ans = VT * math.log(vnew/VT)
        return ans, limiting_applied

This module pins our implementation against that reference, states the
one deliberate divergence, and holds the two implementations (scalar and
branchless/traced) to each other.
"""

import math

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pycircuit.circuit._limiting import _pnjlim, _pnjlim_branchless
from pycircuit.circuit.toolkit import numeric

VT = 0.025852
IS = 1e-13


def _paper_pnjlim(vold, vnew, Is, VT):
    """Listing 1, transcribed verbatim (argument order included)."""
    vc = VT * math.log(VT / (Is * math.sqrt(2.0)))
    if vnew <= vc or abs(vnew - vold) <= 2.0 * VT:
        ans, limiting_applied = vnew, False
    else:
        limiting_applied = True
        if vold > 0.0:
            arg = 1.0 + (vnew - vold) / VT
            if arg > 0.0:
                ans = (vold + VT * math.log(arg))
            else:
                ans = vc
        else:
            ans = VT * math.log(vnew / VT)
    return ans, limiting_applied


CASES = [
    ## (vold, vnew) -- spanning both sides of vc and of the 2*VT escape
    (0.6, 0.62), (0.6, 0.65), (0.6, 0.9), (0.6, 5.0), (0.6, 50.0),
    (0.0, 0.9), (0.0, 5.0), (-0.5, 1.0), (-0.5, 0.05),
    (0.7, 0.70001), (0.2, 0.25), (0.75, 0.2), (0.75, -3.0),
    (0.9, 0.6), (1.0, -1.0),
]


@pytest.mark.parametrize('vold,vnew', CASES)
def test_matches_paper_listing(vold, vnew):
    """Our scalar limiter equals Listing 1 wherever the listing is
    well-defined (see test_positive_guard for the one exception)."""
    ref, _flag = _paper_pnjlim(vold, vnew, IS, VT)
    got = _pnjlim(vnew, vold, VT, IS, numeric)
    assert_allclose(got, ref, rtol=1e-12, atol=0)


@pytest.mark.parametrize('vold,vnew', CASES)
def test_scalar_and_branchless_agree(vold, vnew):
    """The traced form must be the same law: PCNR on the JAX backend uses
    the branchless one, the CPU the scalar one, and a divergence would be
    a backend-dependent operating point."""
    a = _pnjlim(vnew, vold, VT, IS, numeric)
    b = float(_pnjlim_branchless(vnew, vold, VT, IS, np.log))
    assert_allclose(b, a, rtol=1e-12, atol=1e-15)


def test_small_step_escape_is_exact():
    """The paper's escape: a step of at most 2*VT passes through
    UNCHANGED even far above vc.  Before this was implemented such a step
    was compressed by log(1+e) ~ e -- an O(e^2) perturbation that made
    limiting never quite a no-op near the solution."""
    vold = 0.75                      # well above vc
    for dv in (0.0, 1e-9, 0.5 * VT, 1.999 * VT, -1.999 * VT):
        vnew = vold + dv
        assert _pnjlim(vnew, vold, VT, IS, numeric) == vnew
        assert float(_pnjlim_branchless(vnew, vold, VT, IS, np.log)) == vnew
    ## Just outside the escape the limiter engages again.
    vnew = vold + 2.001 * VT
    assert _pnjlim(vnew, vold, VT, IS, numeric) < vnew


def test_large_step_is_compressed():
    """Above vc and beyond the escape, the step is compressed
    logarithmically -- the property the whole device relies on."""
    vold, vnew = 0.6, 50.0
    got = _pnjlim(vnew, vold, VT, IS, numeric)
    assert vold < got < vold + 10 * VT, got
    ref, applied = _paper_pnjlim(vold, vnew, IS, VT)
    assert applied
    assert_allclose(got, ref, rtol=1e-12)


def test_positive_guard_is_a_deliberate_divergence():
    """Our one departure from Listing 1, and why it is kept.

    With a large IS, vc goes NEGATIVE, so a negative vnew can sit above
    vc; the listing's `vold <= 0` branch then evaluates log(vnew/VT) at a
    negative argument, which raises.  We guard with `vnew > 0`.
    """
    big_IS = 1.0
    vc = VT * math.log(VT / (big_IS * math.sqrt(2.0)))
    assert vc < 0.0
    vold, vnew = -0.5, -0.05         # vnew above vc, still negative
    assert vnew > vc and abs(vnew - vold) > 2 * VT

    with pytest.raises(ValueError):  # the listing's failure mode
        _paper_pnjlim(vold, vnew, big_IS, VT)

    got = _pnjlim(vnew, vold, VT, big_IS, numeric)
    assert got == vnew and np.isfinite(got)
    assert np.isfinite(float(_pnjlim_branchless(vnew, vold, VT, big_IS,
                                                np.log)))


def test_no_junction_passes_through():
    """IS <= 0 means there is no junction to limit."""
    assert _pnjlim(5.0, 0.1, VT, 0.0, numeric) == 5.0
    assert float(_pnjlim_branchless(5.0, 0.1, VT, 0.0, np.log)) == 5.0
