
#include <stdint.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <fenv.h>
#include "pcx_abi.h"
#include "radau.h"

#define RC_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)

static long rc_bail(radau_t *s, long why)
{
    feclearexcept(FE_ALL_EXCEPT);
    s->status = why;
    return why;
}

/* one assembly at the stages Ya[a]: each stage's passes through the core
   (formula 4: G in binG, i in bini, C and q in place), then the residual and
   the blocks in `_coupled_stage_system`'s operations over the full width
   (each result stored, so each operation precedes the flag test) */
static long assemble(radau_t *s, long a)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, nn = n * n;
    long i, j, k, rr, cc, ri, ci;
    core_fn_t core = (core_fn_t) s->core_fn;
    double *Y = s->Ya + a * 3 * n, *q = s->qa + a * 3 * n, *iv = s->ia + a * 3 * n;
    double *C = s->Ca + a * 3 * nn, *G = s->Ga + a * 3 * nn;
    double h = s->h;
    for (j = 0; j < 3; j++) {
        /* (a non-finite stage: `passes` declines and the Python evaluates) */
        for (k = 0; k < n; k++) if (!isfinite(Y[j * n + k])) return 2;
        core(s->core, Y + j * n, s->T, 7, 4, 0.0, 0.0, 0.0, 1.0, 0.0,
             s->dummy, s->dummy, s->dummy, s->dummy, s->dummy, C + j * nn, q + j * n,
             s->dummy, s->dummy, s->dummy, s->dummy);
        memcpy(G + j * nn, s->binG, (size_t) nn * sizeof(double));
        memcpy(iv + j * n, s->bini, (size_t) n * sizeof(double));
    }
    feclearexcept(FE_ALL_EXCEPT);
    for (j = 0; j < 3; j++)
        for (k = 0; k < n; k++) s->K[j * n + k] = -(iv[j * n + k] + s->u[j * n + k]);
    for (i = 0; i < 3; i++) {
        const double a0 = s->A[3 * i], a1 = s->A[3 * i + 1], a2 = s->A[3 * i + 2];
        double *f = s->F + i * n;
        for (k = 0; k < n; k++) {
            double sm = 0.0 + a0 * s->K[k];
            sm = sm + a1 * s->K[n + k];
            sm = sm + a2 * s->K[2 * n + k];
            f[k] = (q[i * n + k] - s->qn[k]) - h * sm;
        }
        for (k = 0, rr = 0; k < n; k++) {
            if (k == r) continue;
            s->R[i * m + rr] = f[k];
            rr++;
        }
    }
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            const double hA = h * s->A[3 * i + j];
            const double *Gj = G + j * nn, *Ci = C + i * nn;
            double *blk = s->blk;
            if (i == j) for (k = 0; k < nn; k++) blk[k] = Ci[k] + hA * Gj[k];
            else for (k = 0; k < nn; k++) blk[k] = hA * Gj[k];
            for (ri = 0, rr = 0; ri < n; ri++) {
                if (ri == r) continue;
                for (ci = 0, cc = 0; ci < n; ci++) {
                    if (ci == r) continue;
                    s->Jf[(j * m + cc) * M + (i * m + rr)] = blk[ri * n + ci];
                    cc++;
                }
                rr++;
            }
        }
    }
    if (fetestexcept(RC_FLAGS)) return 3;
    return 0;
}

long hdl_fn(radau_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, M = 3 * m, a = 0, it, i, k, e;
    int64_t M64 = M, one64 = 1, info64 = 0;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    s->assemblies = 0;
    s->iters = 0;
    s->status = 0;
    /* every input finite: anything else is the Python Newton's */
    if (!(isfinite(s->h) && isfinite(s->reltol) && isfinite(s->abstol))) return rc_bail(s, 2);
    for (k = 0; k < 9; k++) if (!isfinite(s->A[k])) return rc_bail(s, 2);
    for (k = 0; k < n; k++) if (!isfinite(s->qn[k])) return rc_bail(s, 2);
    for (k = 0; k < 3 * n; k++) if (!isfinite(s->u[k])) return rc_bail(s, 2);
    memcpy(s->Ya, s->seed, (size_t) (3 * n) * sizeof(double));
    e = assemble(s, 0);
    if (e) return rc_bail(s, e);
    for (it = 0; it < s->maxiter; it++) {
        double scale = 0.0, ynorm = 0.0;
        double *Yp = s->Ya + a * 3 * n, *Yn;
        if (a + 1 >= s->maxa) return rc_bail(s, 6);
        Yn = Yp + 3 * n;
        /* np.linalg.solve(Jbig, -R), then sum|R| (the undamped Newton reads
           it for its overflow alone: a bound, numpy's pairwise sum not
           mirrored) */
        for (k = 0; k < M; k++) {
            if (!(fabs(s->R[k]) < 0x1p1000)) return rc_bail(s, 4);
            s->b[k] = -s->R[k];
        }
        for (k = 0; k < M * M; k++) if (!isfinite(s->Jf[k])) return rc_bail(s, 4);
        ((dgesv64_t) s->dgesv)(&M64, &one64, s->Jf, &M64, s->ipiv, s->b, &M64, &info64);
        /* (numpy's solve clears what LAPACK raised; an info is its LinAlgError) */
        feclearexcept(FE_ALL_EXCEPT);
        if (info64 != 0) return rc_bail(s, 5);
        for (k = 0; k < M; k++) if (!isfinite(s->b[k])) return rc_bail(s, 9);
        for (i = 0; i < 3; i++) {
            const double *yp = Yp + i * n, *d = s->b + i * m;
            double *yn = Yn + i * n, mx = 0.0;
            for (k = 0; k < r; k++) yn[k] = yp[k] + d[k];
            yn[r] = yp[r] + 0.0;
            for (k = r + 1; k < n; k++) yn[k] = yp[k] + d[k - 1];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
            if (walk) {
                e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, yn, yp, 0, s->wn);
                if (e < s->wn) return rc_bail(s, 8);
                feclearexcept(FE_ALL_EXCEPT);
            }
            for (k = 0; k < n; k++) s->stp[k] = yn[k] - yp[k];
            if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(s->stp[k]);
                if (!isfinite(v)) return rc_bail(s, 9);
                if (v > mx) mx = v;
            }
            if (mx > scale) scale = mx;
        }
        a++;
        e = assemble(s, a);
        if (e) return rc_bail(s, e);
        for (i = 0; i < 3; i++)
            for (k = 0; k < n; k++) {
                double v;
                if (k == r) continue;
                v = fabs(Yn[i * n + k]);
                if (v > ynorm) ynorm = v;
            }
        feclearexcept(FE_ALL_EXCEPT);
        s->tolv[0] = s->reltol * ynorm;
        s->tolv[1] = s->tolv[0] + s->abstol;
        if (fetestexcept(RC_FLAGS)) return rc_bail(s, 3);
        s->iters = it + 1;
        if (scale <= s->tolv[1]) {
            s->assemblies = a + 1;
            s->status = 1;
            return 1;
        }
    }
    return rc_bail(s, 11);
}
