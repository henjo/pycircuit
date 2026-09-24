"""Small numerical helpers of the shooting analyses: the DFT, the complex and
transposed solves, GMRES by Arnoldi, the periodic spline weights.
"""
import numpy as np


def freq_analysis(x, t, rms = True, axis=-1, freqoffset = 0):
    """Return dft of equidistant sampled signal x"""
    
    npoints = np.size(x, axis)

    dt = t[1] - t[0]

    if x.dtype in (np.cdouble, np.cdouble):
        X = np.fft.fftshift(np.fft.fft(x, axis=axis),axes=(axis,)) / npoints
        freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
    else:
        freqs = np.fft.fftfreq(npoints, d=dt)[:int(np.ceil(npoints / 2.))]
        slices = [slice(None)] * x.ndim
        slices[axis] = slice(0, len(freqs))
        X = np.fft.fft(x, axis=axis)[tuple(slices)] / npoints
        ## Fold energy from negative frequencies
        X[:,1:] *= np.sqrt(2)

    if not rms:
        X *= np.sqrt(2)

    return freqs, X


def _complex_solve(lu, b):
    """`lu.solve(b)` for a complex `b` against a REAL factorisation.

    Two back-substitutions, not a complex refactorisation: the step
    Jacobians are real, so the solve is real and splits exactly.
    """
    b = np.asarray(b)
    if b.dtype.kind == 'c':
        return lu.solve(b.real) + 1j * lu.solve(b.imag)
    return lu.solve(b)


def _complex_solve_transposed(lu, b):
    """`lu.solve_transposed(b)` for complex `b` against a REAL factorisation --
    two transposed back-substitutions.  `None` if the solver cannot
    transpose (mirrors `solve_transposed`)."""
    b = np.asarray(b)
    if b.dtype.kind == 'c':
        re = lu.solve_transposed(b.real)
        im = lu.solve_transposed(b.imag)
        if re is None or im is None:
            return None
        return re + 1j * im
    return lu.solve_transposed(b)


def _arnoldi_gmres(matvec, b, rtol=1e-12, maxiter=None, reortho=True):
    """GMRES that keeps its Hessenberg matrix and judges its own residual.

    Returns `(x, relres, H, k)`: the solution, the RELATIVE residual, the
    `k x k` Hessenberg matrix of the Krylov basis actually built, and `k`.

    ⚠ WRITTEN RATHER THAN IMPORTED FOR TWO REASONS, AND SPEED IS NEITHER.

    FIRST, `H` IS THE POINT.  Garcia, Romero & Acha (IEEE Trans. Power
    Systems 37(1), 2022) read the Floquet multipliers off exactly this
    matrix -- Ritz values `theta` of `I - M` map back as `lam = 1 - theta`
    -- so a GMRES that discards `H` throws away the spectrum it just
    computed.  `scipy.sparse.linalg.gmres` discards it.

    SECOND, AND THE REASON THIS IS A CORRECTNESS CHANGE: SciPy REPORTS
    BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  When the Krylov space is
    exhausted the next basis vector is numerically zero -- a HAPPY
    breakdown, where the answer is EXACT -- and it comes back as
    `info = 4`.  Trusting that flag turns an exact answer into a
    `RuntimeError`, which is what it did for AM/PM at small offsets and
    why `PAC._gmres_checked` exists to overrule it.  Here the breakdown is
    detected where it happens and returned as the converged answer it is.

    ⚠ REORTHOGONALISED ONCE BY DEFAULT.  Modified Gram-Schmidt loses
    orthogonality as the basis grows, and the Ritz values are read off `H`
    -- so a basis that has drifted gives multipliers that are wrong in a
    way the residual cannot see.  One extra pass is `O(k n)` against the
    matvec's cost, which here is a full replay of the period.

    ⚠ NO RESTARTS.  Restarting discards the basis, which is the object
    this exists to keep.  For the systems here -- `2m` unknowns, `k`
    bounded by `n` -- the full basis is affordable; a caller that needs
    restarts needs a different function and should not silently get one.
    """
    b = np.asarray(b)
    n = b.shape[0]
    kmax = int(min(n, maxiter if maxiter else n))
    beta = float(np.linalg.norm(b))
    if beta == 0.0 or kmax < 1:
        return np.zeros_like(b), 0.0, np.zeros((0, 0)), 0
    Q = [b / beta]
    H = np.zeros((kmax + 1, kmax), dtype=b.dtype)
    for j in range(kmax):
        w = np.asarray(matvec(Q[j]))
        for i in range(j + 1):
            H[i, j] = np.vdot(Q[i], w)
            w = w - H[i, j] * Q[i]
        if reortho:
            for i in range(j + 1):
                c = np.vdot(Q[i], w)
                H[i, j] += c
                w = w - c * Q[i]
        H[j + 1, j] = float(np.linalg.norm(w))
        rhs = np.zeros(j + 2, dtype=b.dtype)
        rhs[0] = beta
        y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], rhs, rcond=None)
        relres = float(np.linalg.norm(H[:j + 2, :j + 1] @ y - rhs)) / beta
        happy = H[j + 1, j] <= 1e-14 * max(beta, 1.0)
        if relres <= rtol or happy or j + 1 == kmax:
            x = np.zeros_like(b)
            for i in range(j + 1):
                x = x + y[i] * Q[i]
            return x, relres, np.array(H[:j + 1, :j + 1]), j + 1
        Q.append(w / H[j + 1, j])
    raise AssertionError('unreachable')


def periodic_spline_weights(t, T, breaks=None):
    """Weights `w` with `sum_j w_j y_j` = the integral over one period of the
    PERIODIC CUBIC SPLINE through `(t_j, y_j)`, `j = 0..n-1`, `t_n = t_0 + T`
    -- or, with `breaks` (node indices where the integrand KINKS), of the
    PIECEWISE not-a-knot cubic spline whose pieces meet at those nodes and
    at node 0.

    ⚠ UNDER LANDED EVENTS THE PIECES BREAK AT THE EVENT NODES (2026-09-21,
    item 2 of the non-uniform-grid list).  A landed edge puts a kink in the
    integrand at its node; a spline that is C^2 across it rings, which is
    why the callers used to keep the trapezoid whenever `event_times` was
    non-empty -- and the trapezoid caps every method's period integrals at
    second order on exactly the grids a clocked circuit gets.  A cubic
    spline fitted PER SEGMENT (not-a-knot ends; a two-node segment is the
    trapezoid, a three-node one the parabola `CubicSpline` builds) never
    crosses a kink.  Node 0 is always a break when any event exists: an
    edge at the drive's t = 0 is dropped by `event_grid` (the period
    boundary is not its to move) and would otherwise sit inside the
    periodic seam.  Measured on `exp(sin) + triangle` (kinks landed, a
    1 + 0.15 sin grid): the periodic spline +4.0e-5 .. -2.2e-8 with no order
    (ringing), the trapezoid +1.8e-5 .. -7.8e-9 at second order, the
    piecewise rule -8.6e-7 .. -1.8e-15 over N = 52 .. 1602; on a smooth
    integrand with only node 0 broken it is fourth order too (-8.9e-7 ..
    -6.1e-13 against the periodic rule's +1.7e-8 .. +2.5e-13), so the seam
    break costs a constant, not an order.

    The higher-order period quadrature for a NON-UNIFORM, EVENT-FREE grid
    (2026-09-21).  The periodic trapezoid rule is spectrally accurate on a
    uniform grid and on an alternating one (two interleaved uniform sums) and
    genuinely O(h^2) on a smoothly varying grid -- measured: radau's diffusion
    constant on a 1 + 0.5 sin grid sat at +4.9e-5 / 1.3e-5 / 3.3e-6 / 8.4e-7
    (N = 100 .. 800) while its period, multipliers and mode invariant were at
    order 5 / 1e-11 / 1e-15, so every method's noise was capped at second
    order by the quadrature on exactly the grids `lte_grid`
    produces.  With these weights: -5.5e-9 / -2.4e-10 / -7.3e-12 / +2.6e-12,
    the reference's floor from N = 400.

    ⚠ THE PERIODIC RULE IS FOR EVENT-FREE GRIDS.  A spline is C^2 across
    every node; a landed source edge puts a KINK in the integrand at that
    node, and a spline through a kink rings where the trapezoid is
    exact-ish -- under events the callers pass `breaks` and get the
    piecewise rule above (they used to keep the trapezoid).  ⚠ ON A UNIFORM
    GRID these equal the trapezoid weights to 1.7e-18 (a periodic spline's
    `sum M_j = 0`), and the callers keep their uniform path, which is
    bit-identical to before.  The spline system is cyclic tridiagonal and
    is solved sparse (O(n)); the weights are `wt - B^T A^{-T} d`, with
    `wt` the trapezoid weights and `d_j = (h_j^3 + h_{j-1}^3) / 24` the
    coefficient of the node's second derivative in the spline integral."""
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spla
    t = np.asarray(t, dtype=float).ravel()
    n = len(t)
    if breaks is not None and len(breaks) > 0:
        from scipy.interpolate import CubicSpline as _CS
        w = np.zeros(n)
        b = sorted(set([0] + [int(i) % n for i in breaks]))
        for k in range(len(b)):
            i0 = b[k]
            i1 = b[k + 1] if k + 1 < len(b) else n
            idx = list(range(i0, i1)) + [i1 % n]
            tt = t[idx].copy()
            if i1 == n:
                tt[-1] = t[0] + float(T)
            if len(tt) == 2:
                np.add.at(w, idx, [0.5 * (tt[1] - tt[0])] * 2)
                continue
            ## the integral of every cardinal spline of the segment at once
            cs = _CS(tt, np.eye(len(tt)), axis=0)
            np.add.at(w, idx, cs.integrate(tt[0], tt[-1]))
        return w
    h = np.empty(n)
    h[:-1] = np.diff(t)
    h[-1] = float(T) + t[0] - t[-1]
    if n < 4 or np.any(h <= 0.0):
        w = 0.5 * h
        w[1:] += 0.5 * h[:-1]
        w[0] += 0.5 * h[-1]
        return w
    j = np.arange(n)
    jm, jp = (j - 1) % n, (j + 1) % n
    hm = h[jm]
    A = _sp.csc_matrix((np.concatenate((hm / 6.0, (hm + h) / 3.0, h / 6.0)),
                        (np.concatenate((j, j, j)), np.concatenate((jm, j, jp)))),
                       shape=(n, n))
    B = _sp.csc_matrix((np.concatenate((1.0 / h, -1.0 / h - 1.0 / hm, 1.0 / hm)),
                        (np.concatenate((j, j, j)), np.concatenate((jp, j, jm)))),
                       shape=(n, n))
    wt = 0.5 * h
    wt[jp] += 0.5 * h
    d = h ** 3 / 24.0
    d = d + d[jm]
    x = _spla.spsolve(A.T.tocsc(), d)
    return wt - np.asarray(B.T @ x).ravel()


## Bound once: `_lu_solve_split` runs once per step of every coupled replay
## (~20k calls per PAC + adjoint row on the PWM fixture), and a function-local
## import plus `np.iscomplexobj` there was 6 % of that row's profile.
## scipy.linalg is loaded by the package's import already.
from scipy.linalg import lu_solve as _sla_lu_solve


def _lu_solve_split(lu, b, trans=0):
    """`scipy.linalg.lu_solve` against a REAL `lu_factor` pair -- a complex
    `b` (an ndarray) is two real back-substitutions (`_complex_solve` for a
    factor that is not a solver object)."""
    if b.dtype.kind == 'c':
        return (_sla_lu_solve(lu, b.real, trans=trans)
                + 1j * _sla_lu_solve(lu, b.imag, trans=trans))
    return _sla_lu_solve(lu, b, trans=trans)


def _cx_collect(a, b):
    """`a + 1j*b` over possibly NESTED lists of arrays (collected adjoint
    samples; a DIRK nests per-stage solves inside per-step entries)."""
    if a is None:
        ## an explicit first stage stores no solve (`D_0 = I`)
        return None
    if isinstance(a, (list, tuple)):
        return type(a)(_cx_collect(x, y) for x, y in zip(a, b))
    return np.asarray(a) + 1j * np.asarray(b)
