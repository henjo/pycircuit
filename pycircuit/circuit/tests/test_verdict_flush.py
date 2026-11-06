"""THE RESIDUAL'S SUMS WITH SUBNORMALS FLUSHED (pcx_verdict.h's
`pcx_resid_sparse`; plan 261's L11): the refactor's residual verdict is made
with the SSE unit's flush-to-zero and denormals-are-zero set, its rounding
bound widened by what a flush can drop (DBL_MIN a term's operand, product or
sum) -- where its sums took a sixth of the C transform loop in microcode
assists on the subnormal entries a converging update holds; where it does not
decide, the gradual sums do (`pcx_resid_gradual`).  A verdict it makes is
still the dense test's: numpy's `max|A @ x - b| / max(max|b|, 1e-300) <=
tol`, made with gradual underflow.  (test_verdict_c's products-underflow
test is the gradual sums' at their thresholds.)"""
import numpy as np
import pytest

from pycircuit.circuit import _cbuild

_C = '''
#include <stdint.h>
#include <stddef.h>
''' + _cbuild._read('pcx_minmax.h') + _cbuild._read('pcx_verdict.h') + '''
long pcx_verdict_probe(const int *Ap, const int *Ai, const double *Ax, double tol,
                       const double *b, const double *x, double *y, double *t, double *cnt,
                       long m)
{
    return (long) pcx_resid_sparse(Ap, Ai, Ax, tol, b, x, y, t, cnt, m);
}
'''
_CDEF = '''long pcx_verdict_probe(const int *Ap, const int *Ai, const double *Ax, double tol,
                       const double *b, const double *x, double *y, double *t, double *cnt,
                       long m);'''


@pytest.fixture(scope='module')
def verdict():
    try:
        ffi, cfn, _key, _cold, _secs = _cbuild.load_kernel(_C, _CDEF, prelude=False,
                                                           entry='pcx_verdict_probe')
    except (_cbuild.CompileError, OSError) as e:
        pytest.skip(f'no C build here: {e}')
    fb = ffi.from_buffer

    def call(Ap, Ai, Ax, tol, b, x):
        m = len(b)
        y, t, cnt = np.empty(m), np.empty(m), np.empty(m)
        return int(cfn(fb('int *', Ap), fb('int *', Ai), fb('double *', Ax), float(tol),
                       fb('double *', b), fb('double *', x), fb('double *', y), fb('double *', t),
                       fb('double *', cnt), m))
    return call


def _case(rng, m, sub_x, sub_a, scale):
    """A sparse m x m A (CSC), a solution x with `sub_x` of its entries
    subnormal and some just above DBL_MIN, `sub_a` of A's subnormal; all of
    x and A times 2^`scale`."""
    pat = rng.random((m, m)) < 0.03
    pat[np.arange(m), np.arange(m)] = True
    cols, rows = np.nonzero(pat.T)
    Ap = np.zeros(m + 1, dtype=np.int32)
    Ap[1:] = np.cumsum(np.bincount(cols, minlength=m))
    Ai = rows.astype(np.int32)
    Ax = rng.standard_normal(len(rows)) * 10.0 ** rng.integers(-3, 3, len(rows))
    x = rng.standard_normal(m) * 10.0 ** rng.integers(-6, 1, m)
    tiny = np.finfo(float).tiny
    x[rng.choice(m, sub_x, replace=False)] = rng.uniform(-1, 1, sub_x) * tiny * 2.0 ** -20
    x[rng.choice(m, 4, replace=False)] = rng.uniform(1, 4, 4) * tiny
    if sub_a:
        Ax[rng.choice(len(Ax), sub_a, replace=False)] = tiny * 2.0 ** -30
    x = np.ldexp(x, scale)
    Ax = np.ldexp(Ax, scale)
    return Ap, Ai, Ax, rows, cols, x


def _dense(Ap, Ai, Ax, rows, cols, m, x, b, tol):
    A = np.zeros((m, m))
    A[rows, cols] = Ax
    return 1 if not (float(np.abs(A @ x - b).max()) / max(float(np.abs(b).max()), 1e-300)
                     <= tol) else 0


@pytest.mark.parametrize('scale', [0, -400, -900], ids=['normal', 'small', 'tiny'])
@pytest.mark.parametrize('seed', range(4))
def test_a_flushed_verdict_is_the_dense_tests(verdict, seed, scale):
    """Wherever the flushed verdict decides, it is numpy's dense test:
    solutions with subnormal entries and entries just past DBL_MIN, A with
    subnormal entries, the whole scaled toward the subnormal range, the
    residual far under, at and far over the tolerance; at the normal scale it
    decides the clear cases."""
    rng = np.random.default_rng(1000 * seed - scale)
    m = 300
    tol = 1e-10
    Ap, Ai, Ax, rows, cols, x = _case(rng, m, sub_x=12, sub_a=5, scale=scale)
    A = np.zeros((m, m))
    A[rows, cols] = Ax
    b0 = A @ x
    decided = []
    for ratio in (0.0, 1e-3, 0.5, 1.0, 2.0, 1e3):
        r = rng.standard_normal(m)
        b = b0 + ratio * tol * max(np.abs(b0).max(), 1e-300) * r / np.abs(r).max()
        want = _dense(Ap, Ai, Ax, rows, cols, m, x, b, tol)
        got = verdict(Ap, Ai, Ax, tol, b, x)
        assert got in (0, 1, 2)
        if got != 2:
            assert got == want, (seed, scale, ratio)
        decided.append((ratio, got))
    if scale == 0:
        assert dict(decided)[0.0] == 0 and dict(decided)[1e3] == 1


def test_a_subnormal_residual_is_never_decided_wrong(verdict):
    """A right-hand side whose largest entry is far below the flushed
    bound's absolute part (~1e-305) -- a residual the flush could hide --
    is decided by the gradual sums, or left to the dense test: never
    decided wrong."""
    rng = np.random.default_rng(7)
    m = 100
    Ap, Ai, Ax, rows, cols, x = _case(rng, m, sub_x=5, sub_a=0, scale=-1000)
    A = np.zeros((m, m))
    A[rows, cols] = Ax
    b = A @ x
    for tol in (1e-10, 1e-3):
        got = verdict(Ap, Ai, Ax, tol, b, x)
        want = _dense(Ap, Ai, Ax, rows, cols, m, x, b, tol)
        assert got == 2 or got == want
