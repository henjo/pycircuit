"""Tests for the Sequence-of-Expressions (SoE) solver and result object."""
import numpy as np
import sympy

from pycircuit.circuit.soe import solve_soe, SoESolution


def _ladder(N):
    s = sympy.Symbol('s')
    R = [sympy.Symbol('R%d' % i, positive=True) for i in range(N)]
    C = [sympy.Symbol('C%d' % i, positive=True) for i in range(N)]
    A = sympy.zeros(N, N)
    b = sympy.zeros(N, 1)
    for i in range(N):
        A[i, i] = 1 / R[i] + s * C[i] + (1 / R[i - 1] if i > 0 else 0)
        if i + 1 < N:
            A[i, i + 1] = -1 / R[i]
            A[i + 1, i] = -1 / R[i]
    b[0] = 1 / R[0]
    return A, b, s, R, C


def _numeric_env(R, C, s, freq):
    env = {R[i]: 100.0 * (i + 1) for i in range(len(R))}
    env.update({C[i]: 1e-9 * (i + 1) for i in range(len(C))})
    env[s] = 1j * 2 * np.pi * freq
    return env


def test_solve_soe_returns_solution_object():
    A, b, s, R, C = _ladder(4)
    sol = solve_soe(A, b)
    assert isinstance(sol, SoESolution)
    assert len(sol) == 4
    assert len(sol.assignments) > 0


def test_soe_assignment_count_is_linear_for_ladder():
    # The SoE stays compact: assignment count grows ~linearly in N (the whole
    # point of the representation vs the exponential exploded N/D).
    counts = [len(solve_soe(*_ladder(N)[:2]).assignments) for N in (4, 8, 12)]
    # near-linear: doubling N roughly doubles the count (not squaring/exploding)
    assert counts[1] < 2.5 * counts[0]
    assert counts[2] < 2.5 * counts[1]


def test_soe_eval_matches_direct_solve():
    N = 5
    A, b, s, R, C = _ladder(N)
    sol = solve_soe(A, b)
    env = _numeric_env(R, C, s, 1e6)
    ws = env[s] * np.logspace(0, 3, 4)  # a small frequency sweep
    params = {k: v for k, v in env.items() if k is not s}
    params[s] = ws
    vals = sol.eval(params)                       # (N, 4)
    assert vals.shape == (N, 4)
    for p, wk in enumerate(ws):
        An = np.array([[complex(A[i, j].subs({**params, s: wk}))
                        for j in range(N)] for i in range(N)])
        bn = np.array([complex(b[i].subs({**params, s: wk})) for i in range(N)])
        xn = np.linalg.solve(An, bn)
        assert np.max(np.abs(vals[:, p] - xn)) < 1e-6 * max(1, np.max(np.abs(xn)))


def test_soe_eval_scalar_params():
    A, b, s, R, C = _ladder(4)
    sol = solve_soe(A, b)
    env = _numeric_env(R, C, s, 1e6)
    vals = sol.eval(env)
    An = np.array([[complex(A[i, j].subs(env)) for j in range(4)]
                   for i in range(4)])
    bn = np.array([complex(b[i].subs(env)) for i in range(4)])
    xn = np.linalg.solve(An, bn)
    assert np.max(np.abs(vals.ravel() - xn)) < 1e-6 * max(1, np.max(np.abs(xn)))


def test_soe_eval_tf_matches():
    N = 5
    A, b, s, R, C = _ladder(N)
    sol = solve_soe(A, b)
    env = _numeric_env(R, C, s, 1e6)
    ws = env[s] * np.logspace(0, 2, 3)
    params = {k: v for k, v in env.items() if k is not s}
    params[s] = ws
    vals = sol.eval(params)
    tf = sol.eval_tf(N - 1, 0, params)
    assert np.max(np.abs(tf - vals[N - 1] / vals[0])) < 1e-9


def test_soe_eval_accepts_symbol_names():
    A, b, s, R, C = _ladder(3)
    sol = solve_soe(A, b)
    env = _numeric_env(R, C, s, 1e6)
    by_name = {str(k): v for k, v in env.items()}
    assert np.allclose(sol.eval(env), sol.eval(by_name))


def test_soe_to_ratio_reconstructs_small_circuit():
    # to_ratio inlines the assignments; on a small circuit it reconstructs the
    # transfer function (and bridges to the GiNaC N/D tools).
    N = 3
    A, b, s, R, C = _ladder(N)
    sol = solve_soe(A, b)
    H = sol.to_ratio(N - 1, 0)
    env = _numeric_env(R, C, s, 1e6)
    ws = env[s] * np.logspace(0, 2, 3)
    params = {k: v for k, v in env.items() if k is not s}
    params[s] = ws
    f = sympy.lambdify([s] + R + C, H, 'numpy')
    args = [ws] + [np.full_like(ws, 100.0 * (i + 1)) for i in range(N)] \
                + [np.full_like(ws, 1e-9 * (i + 1)) for i in range(N)]
    assert np.max(np.abs(f(*args) - sol.eval_tf(N - 1, 0, params))) < 1e-6
