"""Tests for the experimental GiNaC-backed symbolic toolkit.

Skipped unless the compiled ``_ginac_ext`` extension is built (see
``pycircuit/circuit/build_ginac_ext.sh``).
"""
import numpy as np
import pytest
import sympy

pytest.importorskip("pycircuit.circuit._ginac_ext")

from pycircuit.circuit import symbolic_poly, SubCircuit, R, C, VS, gnd
from pycircuit.circuit.toolkit import ginac_toolkit
from pycircuit.circuit.analysis_ss import AC


def _rc_poles(tk):
    s = sympy.Symbol('s', imaginary=True)
    cir = SubCircuit(toolkit=tk)
    cir['VS'] = VS('in', gnd, vac=1)
    cir['R'] = R('in', 'out', r=sympy.Symbol('R', positive=True))
    cir['C'] = C('out', gnd, c=sympy.Symbol('C', positive=True))
    res = AC(cir, toolkit=tk).solve(s, complexfreq=True)
    return res.tf('out', gnd).poles()


def test_ginac_toolkit_available():
    assert ginac_toolkit is not None


def test_ginac_rc_poles_match_symbolic_poly():
    R_, C_ = sympy.symbols('R C', positive=True)
    poles = _rc_poles(ginac_toolkit)
    assert poles == {-1 / (C_ * R_): 1}
    assert poles == _rc_poles(symbolic_poly)


def test_ginac_linearsolver_num_den_matches():
    s = sympy.Symbol('s')
    Rs, Cs = sympy.symbols('R C', positive=True)
    A = np.array([[1 / Rs + s * Cs, -1 / Rs], [-1 / Rs, 1 / Rs]], dtype=object)
    b = np.array([0, 1], dtype=object)

    n_poly, d_poly = symbolic_poly.linearsolver_num_den(A, b)
    n_g, d_g = ginac_toolkit.linearsolver_num_den(A, b)

    for i in range(2):
        assert sympy.simplify(n_poly[i] / d_poly - n_g[i] / d_g) == 0


def test_ginac_large_system_falls_back():
    # Above ginac_max_dim the toolkit must fall back to sympy (guards the CLN
    # coefficient-explosion hang) and still return a correct result.
    n = ginac_toolkit.ginac_max_dim + 3
    s = sympy.Symbol('s')
    A = np.array(sympy.eye(n).tolist(), dtype=object)
    A[0, 0] = 1 + s
    b = np.zeros(n, dtype=object)
    b[0] = 1
    num, den = ginac_toolkit.linearsolver_num_den(A, b)
    assert sympy.simplify(num[0] / den - 1 / (1 + s)) == 0


def test_eval_sweep_rc_lowpass_matches_numpy():
    # H(jw) = 1 / (1 + jwRC): compile once (native), evaluate over a sweep.
    from pycircuit.circuit import _ginac
    w, Rs, Cs = sympy.symbols('w R C', real=True)
    H = 1 / (1 + sympy.I * w * Rs * Cs)
    ws = np.linspace(1e3, 1e7, 16)
    got = _ginac.eval_sweep(H, {w: ws, Rs: 1000.0, Cs: 1e-9})
    ref = 1.0 / (1.0 + 1j * ws * 1000.0 * 1e-9)
    assert got.shape == ws.shape
    assert np.max(np.abs(got - ref)) < 1e-9


def test_eval_sweep_toolkit_method():
    w = sympy.Symbol('w', real=True)
    got = ginac_toolkit.eval_sweep(2 * w + 1, {w: np.array([0.0, 1.0, 2.0])})
    assert np.allclose(got, [1.0, 3.0, 5.0])


def test_eval_sweep_preserves_symbol_assumptions_in_fallback():
    # A symbol carrying assumptions (imaginary s) must still substitute on the
    # sympy fallback path -- force the fallback with a tiny char budget.
    from pycircuit.circuit import _ginac
    s = sympy.Symbol('s', imaginary=True)
    w = sympy.Symbol('w', real=True)
    H = 1 / (1 + s)
    ws = np.array([1.0, 2.0, 3.0])
    orig = _ginac.MAX_COMPILE_CHARS
    _ginac.MAX_COMPILE_CHARS = 1  # force the lambdify fallback (no compiler)
    try:
        got = _ginac.eval_sweep(H.subs(s, sympy.I * w), {w: ws})
    finally:
        _ginac.MAX_COMPILE_CHARS = orig
    ref = 1.0 / (1.0 + 1j * ws)
    assert np.max(np.abs(got - ref)) < 1e-12


def test_eval_sweep_missing_symbol_raises():
    from pycircuit.circuit import _ginac
    a, b_ = sympy.symbols('a b', real=True)
    with pytest.raises(ValueError):
        _ginac.eval_sweep(a + b_, {a: np.array([1.0, 2.0])})


def _rc_ladder(nstage):
    """Symbolic RC ladder A x = b (node voltages, unit source at node 0)."""
    Rs = sympy.symbols('R1:%d' % (nstage + 1), positive=True)
    Cs = sympy.symbols('C1:%d' % (nstage + 1), positive=True)
    s = sympy.Symbol('s', imaginary=True)
    n = nstage
    A = sympy.zeros(n, n)
    for k in range(n):
        g = 1 / Rs[k]
        gnext = 1 / Rs[k + 1] if k + 1 < n else 0
        A[k, k] = g + gnext + s * Cs[k]
        if k + 1 < n:
            A[k, k + 1] = -gnext
            A[k + 1, k] = -gnext
    b = np.zeros(n, dtype=object)
    b[0] = 1 / Rs[0]
    return np.array(A.tolist(), dtype=object), b, s


def test_solve_native_len_and_component_match():
    A, b, s = _rc_ladder(3)
    res = ginac_toolkit.solve_native(A, b)
    assert len(res) == 3
    num, den = symbolic_poly.linearsolver_num_den(A, b)
    for i in range(3):
        assert sympy.simplify(res.component(i) - num[i] / den) == 0


def test_solve_native_tf_matches_symbolic_poly():
    # Transfer function v2/v0 -- cancelled natively in GiNaC, must equal the
    # sympy reference (shared denominator cancels).
    A, b, s = _rc_ladder(3)
    res = ginac_toolkit.solve_native(A, b)
    num, den = symbolic_poly.linearsolver_num_den(A, b)
    tf_ref = sympy.cancel((num[2] / den) / (num[0] / den))
    assert sympy.simplify(res.tf(2, 0) - tf_ref) == 0


def test_solve_native_denominator_matches_up_to_scalar():
    A, b, s = _rc_ladder(3)
    res = ginac_toolkit.solve_native(A, b)
    _, den = symbolic_poly.linearsolver_num_den(A, b)
    ratio = sympy.simplify(sympy.cancel(res.denominator() / den))
    assert ratio.free_symbols == set()  # differ only by a constant scale


def test_solve_native_eval_tf_numeric():
    from pycircuit.circuit import _ginac
    A, b, s = _rc_ladder(3)
    res = _ginac.solve_native(A, b)
    tf = res.tf(2, 0)
    subs = {sy: 1000.0 for sy in tf.free_symbols if sy.name.startswith('R')}
    subs.update({sy: 1e-9 for sy in tf.free_symbols if sy.name.startswith('C')})
    w = sympy.Symbol('w', real=True)
    Hjw = tf.subs(subs).subs(s, sympy.I * w)
    ws = np.logspace(3, 8, 6)
    Hnum = _ginac.eval_sweep(Hjw, {w: ws})
    ref = sympy.lambdify(w, Hjw, 'numpy')(ws)
    assert np.max(np.abs(Hnum - ref)) < 1e-9
    # low-pass ladder: |H| ~ 1 at DC, rolls off at high frequency
    assert abs(Hnum[0]) > 0.9 and abs(Hnum[-1]) < 1e-3


def test_solve_native_preserves_symbol_assumptions():
    # The imaginary-s assumption must be re-applied on the returned sympy expr.
    A, b, s = _rc_ladder(2)
    res = ginac_toolkit.solve_native(A, b)
    returned_s = [sy for sy in res.denominator().free_symbols if sy.name == 's']
    assert returned_s and returned_s[0].is_imaginary
