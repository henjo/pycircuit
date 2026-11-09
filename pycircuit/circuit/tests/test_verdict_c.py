"""The sparse verdicts (`csrc/pcx_verdict.h`): a test the Python makes on
numpy's dense product -- `StructuralKLU`'s refactor residual, nrsolver's
residual test -- decided from the sparse sum and its rounding bound where
the bound decides it, and never against numpy's verdict there; where it does
not, the dense test (the units count it: `vdense`).  Rows are drawn AT the
threshold, so the bound is seen undecided as well as deciding both ways.
The premise -- numpy's dgemv within the bound of the sum in order -- is
test_pcnr_c's `test_numpys_dgemv_rounds_a_row_once_a_nonzero_term`."""
import re

import numpy as np
import pytest

from pycircuit.circuit import _cbuild

#: a probe unit: the two verdicts behind one entry
_SRC = '''#include <math.h>
#include "pcx_minmax.h"
#include "pcx_verdict.h"
long pcx_verdict_probe(long which, const int *Ap, const int *Ai, const double *Ax, double tol,
                       const double *b, const double *x, const double *abstol, double *w, long m)
{
    if (which == 0) return pcx_resid_sparse(Ap, Ai, Ax, tol, b, x, w, w + m, w + 2 * m, m);
    if (which == 2) return pcx_resid_gradual(Ap, Ai, Ax, tol, b, x, w, w + m, w + 2 * m, m);
    return pcx_scale_sparse(Ap, Ai, Ax, x, b, tol, abstol, w, w + m, m);
}
'''
_CDEF = ('long pcx_verdict_probe(long, const int *, const int *, const double *, double, '
         'const double *, const double *, const double *, double *, long);')
_HOLDS, _FAILS, _UNDECIDED = 0, 1, 2


@pytest.fixture(scope='module')
def probe():
    src = re.sub(r'^#include "([\w.]+)"\n', lambda m: _cbuild.csrc(m.group(1)), _SRC,
                 flags=re.MULTILINE)
    try:
        ffi, cfn = _cbuild.load_unit(src, _CDEF, 'pcx_verdict_probe')[:2]
    except _cbuild.CompileError as e:
        pytest.skip(f'no C compiler: {e}')

    def call(which, A, x, b, tol, abstol=None, diag=True):
        m = A.shape[0]
        Ap, Ai, Ax = _csc(A, diag)
        ab = np.zeros(m) if abstol is None else abstol
        w = np.empty(3 * m)
        fb = ffi.from_buffer
        return cfn(which, fb('int *', Ap), fb('int *', Ai), fb('double *', Ax), float(tol),
                   fb('double *', b), fb('double *', x), fb('double *', ab), fb('double *', w), m)
    return call


def _csc(A, diag=True):
    """A's CSC as KLU holds it: `(Ap, Ai, Ax)` on its nonzeros and its
    diagonal (`diag`)."""
    m = A.shape[0]
    rows, cols = np.nonzero(A != 0.0)
    ## (the pattern holds some exact zeros too, as KLU's may: a zero
    ## product joins exactly)
    keep = np.zeros_like(A, dtype=bool)
    keep[rows, cols] = True
    if diag:
        keep[np.arange(m), np.arange(m)] = True
    order = np.lexsort(np.nonzero(keep.T)[::-1])
    cc, rr = np.nonzero(keep.T)
    cc, rr = cc[order], rr[order]
    Ap = np.zeros(m + 1, dtype=np.int32)
    np.add.at(Ap, cc + 1, 1)
    Ap = np.cumsum(Ap).astype(np.int32)
    return Ap, rr.astype(np.int32), np.ascontiguousarray(A[rr, cc])


def _python(A, x, b, tol, diag=True):
    """The Python twin's verdict (`linearsolver._resid_sparse`) as the C's
    code."""
    from pycircuit.circuit import linearsolver
    got = linearsolver._resid_sparse(*_csc(A, diag), tol, b, x)
    return _UNDECIDED if got is None else (_FAILS if got else _HOLDS)


def _sparse(rng, m, signed=True):
    A = np.zeros((m, m))
    for i in range(m):
        k = int(rng.integers(1, min(m, 12) + 1))
        c = rng.choice(m, size=k, replace=False)
        v = rng.normal(size=k) * 10.0 ** rng.integers(-4, 5, size=k)
        A[i, c] = v if signed else np.abs(v)
    return A


def _tally(seen, v, dense):
    if v != _UNDECIDED:
        assert (v == _HOLDS) == dense, 'the sparse verdict contradicts the dense one'
    seen[v] += 1


@pytest.mark.parametrize('m', [3, 17, 64, 65, 250])
def test_the_residual_verdict_never_contradicts_numpys(probe, m):
    """`max |A x - b| / max(max |b|, 1e-300) <= tol`, numpy's way (`arr @ x`),
    against the bound's verdict, b drawn so the residual sits at the
    tolerance (a hair either side) and away from it."""
    rng = np.random.default_rng(20261022 + m)
    tol = 1e-10
    seen = {_HOLDS: 0, _FAILS: 0, _UNDECIDED: 0}
    for _ in range(60):
        A = _sparse(rng, m)
        x = rng.normal(size=m) * 10.0 ** rng.integers(-3, 4, size=m)
        y = A @ x
        u = rng.choice([-1.0, 1.0], size=m) * rng.random(m)
        u[int(rng.integers(m))] = 1.0
        for s in (0.25, 1.0 - 1e-9, 1.0, 1.0 + 1e-9, 4.0):
            b = y + 0.0
            scale = max(float(np.max(np.abs(b))), 1e-300)
            b = y - s * tol * scale * u
            dense = bool(np.max(np.abs(A @ x - b)) / max(float(np.max(np.abs(b))), 1e-300)
                         <= tol)
            _tally(seen, probe(0, A, x, b, tol), dense)
    assert seen[_HOLDS] and seen[_FAILS] and seen[_UNDECIDED], seen


@pytest.mark.parametrize('m', [3, 17, 64, 65, 250])
def test_the_convergence_verdict_never_contradicts_numpys(probe, m):
    """nrsolver's `aF < reltol * (dot(abs(J), ax_next) + aF) + abstol`, all
    rows, against the bound's verdict: residuals drawn at each row's
    threshold (a hair either side), below it and above it."""
    rng = np.random.default_rng(20261023 + m)
    reltol = 1e-3
    seen = {_HOLDS: 0, _FAILS: 0, _UNDECIDED: 0}
    for _ in range(60):
        J = _sparse(rng, m)
        ax = np.abs(rng.normal(size=m) * 10.0 ** rng.integers(-3, 4, size=m))
        abstol = np.full(m, 1e-12) * 10.0 ** rng.integers(0, 4, size=m)
        Is = np.dot(np.abs(J), ax)
        edge = (reltol * Is + abstol) / (1.0 - reltol)
        for mode in ('below', 'at', 'one above', 'edge'):
            f = edge * rng.uniform(0.0, 0.5, size=m)
            if mode == 'at':
                f = edge * (1.0 + 1e-15 * rng.integers(-4, 5, size=m))
            elif mode == 'one above':
                f[int(rng.integers(m))] = edge[0] * 3.0 + 1.0
            elif mode == 'edge':
                i = int(rng.integers(m))
                f[i] = edge[i] * (1.0 + 1e-15 * int(rng.integers(-4, 5)))
            F = f * rng.choice([-1.0, 1.0], size=m)
            aF = np.abs(F)
            dense = bool(np.all(aF < reltol * (np.dot(np.abs(J), ax) + aF) + abstol))
            _tally(seen, probe(1, J, ax, F, reltol, abstol), dense)
    assert seen[_HOLDS] and seen[_FAILS] and seen[_UNDECIDED], seen


#: the smallest subnormal, and offsets from a threshold in its multiples
_TINY = float(np.nextafter(0.0, 1.0))
_AROUND = (-40000, -4000, -400, -24, -2, -1, 0, 1, 2, 24, 400, 4000, 40000)


@pytest.mark.parametrize('m', [3, 17, 64, 250])
def test_the_verdicts_hold_where_products_underflow(probe, m):
    """Rows whose products are subnormal -- each rounded by up to half the
    smallest subnormal, an ABSOLUTE error a relative bound does not cover --
    drawn at each test's threshold, from a few smallest subnormals either
    side (undecided: the bound's absolute term, 8 (k + 2) of them) to
    thousands (decided): neither verdict contradicts numpy's.  (The term is
    for a dgemv that fuses a product into its sum; with this BLAS a row whose
    T_i is below ~7e-309 -- where the relative term's slack no longer covers
    an underflow -- is the sum in order exactly, so the term's removal is not
    seen here.)"""
    rng = np.random.default_rng(20261105 + m)
    tol, reltol = 1e-10, 1e-3
    seen = {_HOLDS: 0, _FAILS: 0, _UNDECIDED: 0}
    seen2 = {_HOLDS: 0, _FAILS: 0, _UNDECIDED: 0}
    for _ in range(40):
        A = _sparse(rng, m)
        x = rng.normal(size=m) * 1e-312 * 10.0 ** rng.integers(-2, 3, size=m)
        y = A @ x
        u = rng.choice([-1.0, 1.0], size=m) * rng.random(m)
        u[int(rng.integers(m))] = 1.0
        ## (max |b| below 1e-300: the scale is 1e-300, the threshold
        ## tol * 1e-300 subnormal)
        for j in _AROUND:
            b = y - (tol * 1e-300 + j * _TINY) * u
            dense = bool(np.max(np.abs(A @ x - b)) / max(float(np.max(np.abs(b))), 1e-300)
                         <= tol)
            _tally(seen, probe(0, A, x, b, tol), dense)
        ax = np.abs(x)
        Is = np.dot(np.abs(A), ax)
        i = int(rng.integers(m))
        ## (an absolute tolerance of ~2 000 smallest subnormals: the rows
        ## at zero hold outright)
        abstol = np.full(m, 1e-320)
        for j in _AROUND:
            f = np.zeros(m)
            f[i] = max((reltol * Is[i] + abstol[i]) / (1.0 - reltol) + j * _TINY, 0.0)
            F = f * rng.choice([-1.0, 1.0], size=m)
            aF = np.abs(F)
            dense = bool(np.all(aF < reltol * (np.dot(np.abs(A), ax) + aF) + abstol))
            _tally(seen2, probe(1, A, ax, F, reltol, abstol), dense)
    assert seen[_HOLDS] and seen[_FAILS], seen
    assert seen2[_HOLDS] and seen2[_FAILS], seen2


def _cancelling(rng, m, delta):
    """A sparse A and an x whose rows' sums cancel to about `delta` of their
    terms (|A| |x| a 1/`delta` times |A x|): a solve's residual test near
    convergence, where the double sums' band reaches the tolerance."""
    A = _sparse(rng, m)
    x = rng.normal(size=m) * 10.0 ** rng.integers(-2, 3, size=m)
    for i in range(m):
        cols = np.flatnonzero(A[i])
        if len(cols) < 2:
            cols = np.array([i, (i + 1) % m]) if m > 1 else cols
            A[i, cols] = rng.normal(size=len(cols))
        last = cols[-1]
        rest = float(np.dot(A[i, cols[:-1]], x[cols[:-1]]))
        A[i, last] = -rest * (1.0 + delta * rng.uniform(0.5, 2.0)) / x[last]
    return A, x


@pytest.mark.parametrize('m', [3, 17, 64, 250])
def test_the_wide_pass_decides_where_the_double_sums_band_does_not(probe, m):
    """(W) (pcx_verdict.h's WHERE NEITHER DECIDES): rows whose |A| |x| is
    1e5-1e8 times |A x|, residuals from a hundredth of the tolerance to ten
    times it and a hair either side of it: the verdict -- the wide pass
    where the double sums' band does not decide -- never contradicts numpy's,
    decides where that band alone did not, and is the Python twin's."""
    rng = np.random.default_rng(20261109 + m)
    tol = 1e-8
    seen = {_HOLDS: 0, _FAILS: 0, _UNDECIDED: 0}
    only_wide = 0
    for _ in range(30):
        A, x = _cancelling(rng, m, 10.0 ** rng.uniform(-8, -5))
        y = A @ x
        u = rng.choice([-1.0, 1.0], size=m) * rng.random(m)
        u[int(rng.integers(m))] = 1.0
        scale = max(float(np.max(np.abs(y))), 1e-300)
        for s in (0.01, 0.1, 0.5, 1.0 - 1e-9, 1.0 + 1e-9, 2.0, 10.0):
            b = y - s * tol * scale * u
            dense = bool(np.max(np.abs(A @ x - b)) / max(float(np.max(np.abs(b))), 1e-300)
                         <= tol)
            got = probe(0, A, x, b, tol)
            _tally(seen, got, dense)
            only_wide += got != _UNDECIDED and probe(2, A, x, b, tol) == _UNDECIDED
            assert _python(A, x, b, tol) == got
    assert only_wide > 0 and seen[_HOLDS] and seen[_UNDECIDED], (only_wide, seen)


def _snan():
    return np.frombuffer(np.uint64(0x7FF0000000000001).tobytes(), dtype=np.float64)[0]


def _signalling(v):
    """A NaN whose quiet bit is clear."""
    return bool(np.isnan(v)) and not int(np.float64(v).view(np.uint64)) >> 51 & 1


@pytest.mark.parametrize('m', [3, 17, 64])
def test_a_value_not_finite_is_decided_where_numpy_raises_nothing(probe, m):
    """(N): a quiet or signalling NaN or an infinity placed in x, b or A,
    rows ordinary or one past 1e300: where numpy's dense test raises a flag it
    reports (its matmul's 0 x inf, a signalling NaN, an overflow) the verdict
    is undecided -- the dense test hands back for numpy's warning; elsewhere
    it is numpy's (it fails), and it is made for every quiet NaN and every
    infinite b on ordinary rows; the Python twin's the same."""
    rng = np.random.default_rng(20261110 + m)
    tol = 1e-8
    made = want = 0
    for _ in range(12):
        A0 = _sparse(rng, m)
        x0 = rng.normal(size=m)
        b0 = A0 @ x0
        for where in ('x', 'b', 'A'):
            for what in (np.nan, np.inf, -np.inf, _snan()):
                for huge in (False, True):
                    A, x, b = A0.copy(), x0.copy(), b0.copy()
                    i = int(rng.integers(m))
                    if huge:
                        A[i] *= 1e301
                    if where == 'x':
                        x[int(rng.integers(m))] = what
                    elif where == 'b':
                        b[i] = what
                    else:
                        r, c = np.nonzero(A)
                        k = int(rng.integers(len(r)))
                        A[r[k], c[k]] = what
                    try:
                        with np.errstate(over='raise', invalid='raise', divide='raise',
                                         under='ignore'):
                            res = np.abs(A @ x - b)
                        raised = False
                    except FloatingPointError:
                        raised = True
                    got = probe(0, A, x, b, tol)
                    assert _python(A, x, b, tol) == got
                    if raised:
                        assert got == _UNDECIDED, (where, what, huge)
                        continue
                    dense = float(res.max()) / max(float(np.abs(b).max()), 1e-300) <= tol
                    if got != _UNDECIDED:
                        assert (got == _HOLDS) == dense, (where, what, huge)
                    if not huge and (where == 'b' or np.isnan(what)) and not _signalling(what):
                        want += 1
                        made += got == _FAILS
    assert want and made == want, (made, want)
    ## (an infinity in x where its column holds no entry: no row's sum sees
    ## it, but the dense product's 0 x inf raises invalid in every row)
    A = _sparse(rng, m)
    A[:, 0] = 0.0
    x = rng.normal(size=m)
    b = A @ x
    x[0] = np.inf
    with pytest.raises(FloatingPointError), np.errstate(invalid='raise'):
        A @ x
    assert probe(0, A, x, b, tol, diag=False) == _UNDECIDED
    assert _python(A, x, b, tol, diag=False) == _UNDECIDED
