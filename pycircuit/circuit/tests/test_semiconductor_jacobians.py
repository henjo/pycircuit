# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""A semiconductor's Jacobian must come from its toolkit, and must be exact.

``test_element_jacobians.py`` pins the invariant for elements that write a
``G()`` stamp by hand.  The devices in :mod:`pycircuit.circuit.semiconductors`
take the opposite approach -- they write *no* stamp, and the Jacobian is
obtained by differentiating ``eval_i_pure``/``eval_q_pure``.  That makes a
different thing worth pinning: *which* differentiation runs, and whether it is
exact when it could be.

The defect this guards against: ``Semiconductor.G`` used to choose its method
with an inline ``try: import jax`` plus ``isinstance(x, jax.numpy.ndarray)``,
i.e. by whether JAX happened to be *installed* rather than by which toolkit was
running -- the same defect ``Toolkit.jacobian`` and ``Toolkit.derivative`` were
introduced to remove from ``elements.py``, which this module was never migrated
onto.  Under the symbolic toolkit it fell through to a hand-rolled central
difference that perturbed a *sympy symbol* by ``1e-6`` and divided by ``2e-6``,
producing a numerical approximation -- simultaneously enormous and inexact -- of
an expression sympy differentiates exactly.

Two further gaps sat behind it, both found only by running the code: neither
backend defined ``minimum``, so ``Varactor``'s junction clamp reached for
``jax.numpy.minimum`` under *every* toolkit (crashing on sympy, and silently
pulling JAX into a plain numpy circuit), and ``_symbolic`` lacked ``maximum``.
"""

import numpy as np
import pytest
import sympy

from pycircuit.circuit import _numeric, _symbolic
from pycircuit.circuit.semiconductors import BJT, JFET, ZenerDiode, Varactor
from pycircuit.circuit.toolkit import numeric, symbolic


## The symbolic toolkit carries k and q as symbols; substitute the numeric
## toolkit's own values so the two sides are actually comparable.  (Using
## full CODATA constants instead gives a ~6% disagreement that looks like a
## real bug and is not: exp(v/VT) amplifies a 0.04% VT error at v/VT ~ 25.)
CONSTANTS = {_symbolic.kboltzmann: _numeric.kboltzmann,
             _symbolic.qelectron: _numeric.qelectron}


def _symbolic_jacobian_at(device_cls, x, method='G'):
    """Evaluate the symbolic Jacobian of ``device_cls`` at operating point ``x``."""
    syms = sympy.symbols('n0:%d' % len(x))
    device = device_cls(*range(1, len(x) + 1), toolkit=symbolic)
    J = getattr(device, method)(list(syms))
    at = dict(zip(syms, x))
    return np.array([[float(sympy.N(J[k, j].subs(at).subs(CONSTANTS)))
                      for j in range(len(x))]
                     for k in range(J.shape[0])])


def _finite_difference(device_cls, x, method='i', eps=1e-6):
    device = device_cls(*range(1, len(x) + 1), toolkit=numeric)
    f = getattr(device, method)
    x = np.asarray(x, dtype=float)
    cols = []
    for j in range(len(x)):
        x_plus, x_minus = x.copy(), x.copy()
        x_plus[j] += eps
        x_minus[j] -= eps
        cols.append((f(x_plus) - f(x_minus)) / (2 * eps))
    return np.column_stack(cols)


@pytest.mark.parametrize('device_cls,x,label', [
    (BJT, [0.7, 0.65, 0.0], 'forward active'),
    (BJT, [0.3, 0.65, 0.0], 'saturation'),
    (JFET, [1.0, -1.0, 0.0], 'saturation'),
    (JFET, [0.1, -0.5, 0.0], 'triode'),
    (ZenerDiode, [0.7, 0.0], 'forward'),
    (ZenerDiode, [-5.2, 0.0], 'breakdown'),
])
def test_symbolic_jacobian_matches_finite_difference(device_cls, x, label):
    """The exact symbolic Jacobian must agree with a numeric finite difference.

    Two independent routes to the same matrix: sympy differentiating the model
    symbolically, and numpy finite-differencing ``i()``.  Operating points are
    chosen to straddle each device's distinct regions, since a Jacobian can be
    right in one region and wrong in another.
    """
    exact = _symbolic_jacobian_at(device_cls, x)
    approx = _finite_difference(device_cls, x)

    scale = np.abs(approx).max()
    assert np.abs(exact - approx).max() / scale < 1e-4, (
        '%s %s: symbolic Jacobian disagrees with finite difference'
        % (device_cls.__name__, label))


@pytest.mark.parametrize('device_cls,x', [
    (BJT, [0.7, 0.65, 0.0]),
    (JFET, [1.0, -1.0, 0.0]),
    (ZenerDiode, [0.7, 0.0]),
])
def test_symbolic_jacobian_is_exact_not_differenced(device_cls, x):
    """No finite-difference artifacts may appear in a *symbolic* Jacobian.

    The regression: a central difference over symbols leaves the perturbation
    ``1e-6`` and its reciprocal ``500000.0`` littered through the expression.
    An exact derivative contains neither.
    """
    syms = sympy.symbols('n0:%d' % len(x))
    device = device_cls(*range(1, len(x) + 1), toolkit=symbolic)
    text = str(device.G(list(syms)))

    for artifact in ('1.0e-6', '500000.0'):
        assert artifact not in text, (
            '%s: symbolic G looks central-differenced (found %r)'
            % (device_cls.__name__, artifact))


def test_varactor_capacitance_is_symbolic():
    """``Varactor.C`` must work under the symbolic toolkit at all.

    It did not: ``eval_q_pure`` reached for ``jax.numpy.minimum``, which
    rejects sympy objects.  Pins both the backend-parity fix and the
    charge-model differentiation path (``eval_q_pure``, not ``eval_i_pure``).
    """
    va, vb = sympy.symbols('va vb')
    C = Varactor(1, 2, toolkit=symbolic).C([va, vb])

    assert C[0, 0].free_symbols, 'Varactor.C should depend on its bias'
    ## Charge is conserved: the two terminals carry equal and opposite charge.
    assert sympy.simplify(C[0, 0] + C[0, 1]) == 0


def test_varactor_capacitance_matches_finite_difference():
    """The exact ``dq/dv`` must agree with differencing ``q()`` numerically."""
    x = [0.5, 0.0]
    exact = _symbolic_jacobian_at(Varactor, x, method='C')
    approx = _finite_difference(Varactor, x, method='q')

    scale = np.abs(approx).max()
    assert np.abs(exact - approx).max() / scale < 1e-4


@pytest.mark.parametrize('name', ['minimum', 'maximum', 'where'])
def test_backend_parity(name):
    """Both backends must define the primitives the device models call.

    ``minimum`` was missing from both and ``maximum`` from ``_symbolic``, so
    ``Varactor`` silently fell through to ``jax.numpy`` -- under the *numeric*
    toolkit too, which is how an optional accelerator ends up stamping a plain
    numpy circuit.  Same shape as the ``tanh`` and ``real`` gaps before it.
    """
    assert hasattr(_numeric, name), '_numeric lacks %s' % name
    assert hasattr(_symbolic, name), '_symbolic lacks %s' % name


def test_semiconductors_do_not_need_jax(monkeypatch):
    """A numeric semiconductor circuit must not import JAX.

    The old code chose its differentiation path on whether JAX was importable,
    so the numeric answer changed depending on an unrelated install.  Blocking
    the import is the cheapest way to state that it must not matter -- the same
    mechanism ``test_jax_optional.py`` uses.
    """
    import sys
    monkeypatch.setitem(sys.modules, 'jax', None)
    monkeypatch.setitem(sys.modules, 'jax.numpy', None)

    x = np.array([0.7, 0.65, 0.0])
    G = BJT(1, 2, 3, toolkit=numeric).G(x)
    assert np.isfinite(G).all()

    ## Varactor is the one that reached for jax.numpy.minimum by name.
    C = Varactor(1, 2, toolkit=numeric).C(np.array([0.5, 0.0]))
    assert np.isfinite(C).all()
