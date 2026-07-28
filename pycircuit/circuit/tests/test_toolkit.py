# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""``Toolkit.__getattr__`` delegates unknown names to the backend module --
that's the main extension mechanism (architecture.md section 2) and earns
its keep. But a name missing from *both* the toolkit class and the backend
used to raise naming the backend module, e.g. "module
'pycircuit.circuit._numeric' has no attribute 'sparse'" -- which sends a
reader to the wrong file, since the toolkit is what a caller actually named
(architecture.md P2). This pins the corrected message, plus the related fix
that made ``sparse`` itself safe to read on any toolkit rather than only on
``SparseNumericToolkit``.
"""
import numpy as np
import pytest
import sympy

from pycircuit.circuit.toolkit import (numeric, sparse_numeric,
                                        SymbolicPolyToolkit, NumDenMixin,
                                        symbolic_poly, ddd_toolkit)


def test_missing_attribute_names_the_toolkit_not_just_the_backend():
    with pytest.raises(AttributeError) as excinfo:
        numeric.this_primitive_does_not_exist

    message = str(excinfo.value)
    assert 'NumericToolkit' in message
    assert 'this_primitive_does_not_exist' in message


def test_existing_backend_primitive_still_delegates():
    """The fix must not disturb the delegation itself, only its errors."""
    assert numeric.cos(0.0) == 1.0


def test_sparse_flag_is_safe_on_any_toolkit():
    """``sparse`` used to be declared only on ``SparseNumericToolkit``, unlike
    ``symbolic``/``poly``/``jax``, which are always safe to read (architecture.md
    P2's secondary note). Now on the base like the others.
    """
    assert numeric.sparse is False
    assert sparse_numeric.sparse is True


def test_ddd_toolkit_shares_the_contract_not_the_class():
    """``DDDToolkit`` used to subclass ``SymbolicPolyToolkit``, overriding
    nearly everything it inherited except a two-line ``linearsolver`` wrapper
    -- claiming a kinship the code didn't use (architecture.md P3). Both now
    implement ``NumDenMixin``'s contract (``linearsolver_num_den`` ->
    ``linearsolver`` for free, plus ``supports('num_den')``) as siblings,
    each free to answer ``ac_solution``/``noise_psd`` its own way.
    """
    assert isinstance(symbolic_poly, NumDenMixin)
    assert isinstance(ddd_toolkit, NumDenMixin)
    assert not isinstance(ddd_toolkit, SymbolicPolyToolkit)


def test_num_den_mixin_linearsolver_matches_across_siblings():
    """The mixin's ``linearsolver`` gives the same, correct answer for both
    a plain algebraic solve, regardless of which toolkit's own
    ``linearsolver_num_den`` produced the numerator/denominator.
    """
    A = np.array(sympy.Matrix([[2, 0], [0, 2]]))
    b = np.array([sympy.Integer(4), sympy.Integer(6)], dtype=object)

    expected = [sympy.Integer(2), sympy.Integer(3)]
    assert list(symbolic_poly.linearsolver(A, b)) == expected
    assert list(ddd_toolkit.linearsolver(A, b)) == expected
