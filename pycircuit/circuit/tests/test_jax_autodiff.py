import pytest
import numpy as np
from pycircuit.circuit.toolkit import numeric as num_tk, jaxtoolkit as jax_tk
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.elements import Diode


@pytest.mark.skipif(jax_tk is None, reason='jax not installed')
def test_jax_autodiff_diode_G_matches_numeric():
    """A JAX-toolkit ``Diode``'s ``G`` and ``i`` -- produced by autodiff over
    ``eval_i_pure`` through ``toolkit.jacobian`` -- must agree with the
    numeric-toolkit stamp.

    This is the live autodiff path: ``Diode.G``/``Diode.i`` guard on
    ``self.toolkit.supports('autodiff')`` and call ``eval_i_pure`` /
    ``toolkit.jacobian(self.eval_i_pure, ...)`` directly, with no per-instance
    machinery involved. It replaces a prior version of this test that
    exercised a different, per-instance mechanism
    (``toolkit.generate_eval_i_and_G`` installing ``eval_i_and_G`` on the
    element) which no shipped element ever triggers automatically -- see
    architecture.md P7, which this test's rewrite closes out.
    """
    x = np.array([0.7, 0.0])  # 0.7V forward bias

    d_num = Diode(0, 1)
    d_num.toolkit = num_tk

    d_jax = Diode(0, 1)
    d_jax.toolkit = jax_tk

    manual_i = np.asarray(d_num.i(x, defaultepar), dtype=float)
    manual_G = np.asarray(d_num.G(x, defaultepar), dtype=float)

    jax_i = np.asarray(d_jax.i(x, defaultepar), dtype=float)
    jax_G = np.asarray(d_jax.G(x, defaultepar), dtype=float)

    np.testing.assert_array_almost_equal(manual_i, jax_i)
    np.testing.assert_array_almost_equal(manual_G, jax_G)
