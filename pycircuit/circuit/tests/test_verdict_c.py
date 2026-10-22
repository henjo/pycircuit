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

    def call(which, A, x, b, tol, abstol=None):
        m = A.shape[0]
        rows, cols = np.nonzero(A != 0.0)
        ## (the pattern holds some exact zeros too, as KLU's may: a zero
        ## product joins exactly)
        keep = np.zeros_like(A, dtype=bool)
        keep[rows, cols] = True
        keep[np.arange(m), np.arange(m)] = True
        order = np.lexsort(np.nonzero(keep.T)[::-1])
        cc, rr = np.nonzero(keep.T)
        cc, rr = cc[order], rr[order]
        Ap = np.zeros(m + 1, dtype=np.int32)
        np.add.at(Ap, cc + 1, 1)
        Ap = np.cumsum(Ap).astype(np.int32)
        Ai = rr.astype(np.int32)
        Ax = np.ascontiguousarray(A[rr, cc])
        ab = np.zeros(m) if abstol is None else abstol
        w = np.empty(3 * m)
        fb = ffi.from_buffer
        return cfn(which, fb('int *', Ap), fb('int *', Ai), fb('double *', Ax), float(tol),
                   fb('double *', b), fb('double *', x), fb('double *', ab), fb('double *', w), m)
    return call


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
