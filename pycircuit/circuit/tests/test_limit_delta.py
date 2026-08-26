# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""`limit_delta`: the limiter law with no device physics in it.

`pnjlim` compresses a junction logarithmically, `fetlim` works about a
threshold, `limvds` has hard-coded breakpoints -- each needs a model that
HAS such a quantity.  `deltalim` only says *do not move this branch more
than `vmax` volts in one Newton step*, so it applies to any branch at
all.  Roadmap sec. 3 asked for it ("35 hand-written lines of FET Newton
limiting per model"); sec. 33 is why it was finally built.

**It is not invented here.**  `compact.PspMosLongChannel.limit` already
bounds d, g and b by a symmetric ``|delta| <= 1 V`` about the source, by
hand.  This is that law, named and shared.

The measurement that motivated it is in sec. 33 and repeated at the
bottom of this file: `EkvNmosHdl` declares `limit_identity` on its bulk
-- a probe with no law -- and is the one model still falling back to the
ordinary solver from a wild start, because sec. 27's fix limits the SEED
through each device's own law and an identity law clamps nothing.
"""

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd, hdl
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit._limiting import apply_limit
from pycircuit.utilities.param import Parameter

cm.default_toolkit = numeric


## ---------------------------------------------------------------------
## 1. The law
## ---------------------------------------------------------------------

@pytest.mark.parametrize('vnew,vold,vmax,want', [
    (5.0, 0.0, 1.0, 1.0),          # rising, bites
    (-5.0, 0.0, 1.0, -1.0),        # falling, bites symmetrically
    (0.3, 0.0, 1.0, 0.3),          # inside the bound, untouched
    (0.0, 0.0, 1.0, 0.0),          # no motion at all
    (1.0, 0.0, 1.0, 1.0),          # exactly at the bound is NOT a bite
    (12.0, 10.0, 0.5, 10.5),       # about a non-zero previous value
    (-12.0, -10.0, 0.5, -10.5),
])
def test_the_law_bounds_the_excursion(vnew, vold, vmax, want):
    got = apply_limit('delta', vnew, vold, (vmax,), numeric)
    assert got == pytest.approx(want, rel=0, abs=0)


def test_not_biting_returns_vnew_itself():
    """The write-back decides "did limiting fire?" by comparing
    ``vlim == vnew``, and `pcnr.limit_block` by ``lim != raw``.  A law
    that reconstructed an equal value would report a bite that did not
    happen -- so this is identity, not equality."""
    v = 0.3
    assert apply_limit('delta', v, 0.0, (1.0,), numeric) is v


def test_the_law_is_symmetric():
    """A blunt limiter must not smuggle a direction preference into a
    model that has none."""
    for d in (0.2, 1.0, 3.0, 1e6):
        up = apply_limit('delta', d, 0.0, (1.0,), numeric)
        dn = apply_limit('delta', -d, 0.0, (1.0,), numeric)
        assert up == -dn


## ---------------------------------------------------------------------
## 2. Declaring it
## ---------------------------------------------------------------------

class DeltaDiode(Behavioural):
    """A diode limited by the blunt law instead of `pnjlim`."""
    terminals = ['plus', 'minus']
    instparams = [Parameter('IS', default=1e-14),
                  Parameter('vt', default=0.025865),
                  Parameter('vmax', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus, 'd')
        v = hdl.var(hdl.limit_delta(b.V, vmax), 'vd')            # noqa: F821
        return (Contribution(b.I, IS * (hdl.limexp(v / vt) - 1.0)),)  # noqa


def test_explain_names_the_law():
    el = DeltaDiode('p', gnd)
    el.update_iparv()
    text = el.explain()
    assert 'deltalim' in text
    assert '1 $limit probe' in text


def test_a_delta_probe_qualifies_the_device_for_vector_pcnr():
    """The reason it was built: an identity probe admits a device to
    PCNR but gives its seed nothing to clamp (sec. 27, 33).  A delta
    probe admits it AND carries a law."""
    el = DeltaDiode('p', gnd)
    el.update_iparv()
    assert 'vector route' in el.explain()
    assert 'deltalim' in el.explain().split('PCNR:')[1]


def test_the_generated_limiter_bites_through_the_element():
    el = DeltaDiode('p', gnd)
    el.update_iparv()
    x0 = np.array([0.0, 0.0])
    for target, want in ((0.5, 0.5), (5.0, 1.0), (-5.0, -1.0)):
        lim = np.asarray(el.limit(np.array([target, 0.0]), x0, defaultepar),
                         float)
        assert lim[0] - lim[1] == pytest.approx(want, rel=1e-12)


def test_vmax_is_carried_as_a_parameter_not_baked_in():
    a = DeltaDiode('p', gnd, vmax=0.25)
    a.update_iparv()
    lim = np.asarray(a.limit(np.array([5.0, 0.0]), np.array([0.0, 0.0]),
                             defaultepar), float)
    assert lim[0] - lim[1] == pytest.approx(0.25, rel=1e-12)


## ---------------------------------------------------------------------
## 3. What it refuses
## ---------------------------------------------------------------------

@pytest.mark.parametrize('bad', [0.0, -1.0, -1e-30])
def test_a_non_positive_bound_is_refused(bad):
    """It would pin the branch at its previous value for ever -- a stall,
    not a limit -- and a solver that never moves is far harder to
    diagnose than a build-time error."""
    b = Branch(cm.Node('p'), cm.Node('m'), 'd')
    with pytest.raises(ValueError, match='POSITIVE'):
        hdl.limit_delta(b.V, bad)


def test_a_parameter_bound_is_allowed():
    """`vmax` may be a parameter expression -- only a NUMBER can be
    checked at build time, and a symbol is not one."""
    b = Branch(cm.Node('p'), cm.Node('m'), 'd')
    assert hdl.limit_delta(b.V, sympy.Symbol('vmax', real=True)) is not None


def test_it_may_be_grouped_unlike_an_identity_probe():
    """`limit_identity` cannot join a `limit_together` group -- the
    grouped write-back holds a non-biting probe as a CONSTRAINT, and an
    identity probe has no opinion to hold.  A delta probe has one."""
    bd = Branch(cm.Node('d'), cm.Node('s'), 'ds')
    bg = Branch(cm.Node('g'), cm.Node('s'), 'gs')
    grouped = hdl.limit_together(hdl.limit_delta(bd.V, 1.0),
                                 hdl.limit_fet(bg.V, 0.7))
    assert len(grouped) == 2
    with pytest.raises(ValueError, match='cannot be grouped'):
        hdl.limit_together(hdl.limit_identity(bd.V), hdl.limit_fet(bg.V, 0.7))
