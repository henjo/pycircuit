
#include <stdint.h>
#include <stddef.h>
#include <math.h>
typedef long (*core_fn_t)(const void *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
typedef long (*walk_fn_t)(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n);
typedef void (*dgesv_t)(const int64_t *, const int64_t *, double *, const int64_t *,
                        int64_t *, double *, const int64_t *, int64_t *);
typedef void (*dgetrf_t)(const int *, const int *, double *, const int *, int *, int *);
typedef void (*dgetrs_t)(const char *, const int *, const int *, const double *, const int *,
                         const int *, double *, const int *, int *, size_t);
typedef void (*dgemv_t)(int, int, int64_t, int64_t, double, const double *, int64_t,
                        const double *, int64_t, double, double *, int64_t);
typedef struct {
    const void *core; void *core_fn;
    double T; long formula; double a0, a1, a2, h, theta;
    const double *q1, *q2, *iq1, *u;
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    void *dgesv, *dgetrf, *dgetrs, *dgemv;
    long n, iref, maxiter, chord, first_bits;
    double reltol; const double *abstol, *xtol;
    const double *Cin0; const double *xseed; const double *cx, *cC; long have_cache, hit;
    double *C0, *q0, *iq0, *Geq0, *F0, *J0;
    double *C, *q, *iq, *Geq, *F, *J;
    double *xf, *xnf, *x0f;
    double *A, *absJ;
    double *b, *Fr, *x, *xn, *xd, *ax, *Is;
    int64_t *ipiv64; int *ipiv32;
    long iters, last_eval, status;
    long fold_j; double *xc, *Cj, *qj, *iqj, *Geqj, *Jj;
} newton_t;

long hdl_fn(newton_t *s)
{
    long n = s->n, m = n - 1, r = s->iref, it, i, j;
    core_fn_t core = (core_fn_t) s->core_fn;
    walk_fn_t walk = (walk_fn_t) s->walk_fn;
    dgemv_t gemv = (dgemv_t) s->dgemv;
    double prev = 0.0;
    int have_prev = 0;
    int64_t m64 = m, one64 = 1, info64 = 0;
    int m32 = (int) m, one32 = 1, info32 = 0;
    char trans = 'N';
    s->iters = 0;
    s->status = 0;
    s->last_eval = -1;
    s->hit = 0;
    /* the seed, reduced (remove_row_col), and _C_lookup's equality test of
       the full seed it evaluates against the cached state */
    for (i = 0, j = 0; i < n; i++) { if (i == r) continue; s->x[j++] = s->xseed[i]; }
    s->Cin0 = 0;
    if (s->have_cache) {
        int eq = 1;
        for (i = 0; i < r; i++) if (!(s->cx[i] == s->x[i])) { eq = 0; break; }
        if (eq && !(s->cx[r] == 0.0)) eq = 0;
        if (eq) for (i = r; i < m; i++) if (!(s->cx[i + 1] == s->x[i])) { eq = 0; break; }
        if (eq) { s->Cin0 = s->cC; s->hit = 1; }
    }
    s->first_bits = s->Cin0 ? 5 : 7;
    for (it = 0; it < s->maxiter; it++) {
        const double *Cin;
        long bits;
        double *F, *J, *C, *q, *iq, *Geq;
        int full = (it == 0) || !s->chord;
        if (it == 0) { F = s->F0; J = s->J0; C = s->C0; q = s->q0; iq = s->iq0; Geq = s->Geq0; }
        else { F = s->F; J = s->J; C = s->C; q = s->q; iq = s->iq; Geq = s->Geq; }
        /* insert_row: the full x of this evaluation */
        for (i = 0; i < r; i++) s->xf[i] = s->x[i];
        s->xf[r] = 0.0;
        for (i = r; i < m; i++) s->xf[i + 1] = s->x[i];
        for (i = 0; i < n; i++) if (!isfinite(s->xf[i])) { s->status = 2; return 2; }
        if (full) {
            bits = (it == 0) ? s->first_bits : 7;
            Cin = (it == 0 && s->Cin0) ? s->Cin0 : C;
        } else {
            bits = 4;
            Cin = s->Cin0 ? s->Cin0 : s->C0;
        }
        core(s->core, s->xf, s->T, bits, s->formula, s->a0, s->a1, s->a2, s->h, s->theta,
             s->q1, s->q2, s->iq1, s->u, Cin, C, q, iq, Geq, F, J);
        s->last_eval = it;
        for (i = 0, j = 0; i < n; i++) {
            if (i == r) continue;
            s->Fr[j] = F[i];
            s->b[j] = -F[i];
            j++;
        }
        if (full) {
            long ri, ci, rr, cc;
            for (ri = 0, rr = 0; ri < n; ri++) {
                if (ri == r) continue;
                for (ci = 0, cc = 0; ci < n; ci++) {
                    double v;
                    if (ci == r) continue;
                    v = J[ri * n + ci];
                    s->A[cc * m + rr] = v;
                    s->absJ[rr * m + cc] = fabs(v);
                    cc++;
                }
                rr++;
            }
        }
        if (s->chord) {
            if (it == 0) {
                for (i = 0; i < m * m; i++) if (!isfinite(s->A[i])) { s->status = 3; return 3; }
                ((dgetrf_t) s->dgetrf)(&m32, &m32, s->A, &m32, s->ipiv32, &info32);
                if (info32 != 0) { s->status = 4; return 4; }
            }
            for (i = 0; i < m; i++) if (!isfinite(s->b[i])) { s->status = 5; return 5; }
            ((dgetrs_t) s->dgetrs)(&trans, &m32, &one32, s->A, &m32, s->ipiv32, s->b, &m32,
                                   &info32, 1);
            if (info32 != 0) { s->status = 6; return 6; }
        } else {
            ((dgesv_t) s->dgesv)(&m64, &one64, s->A, &m64, s->ipiv64, s->b, &m64, &info64);
            if (info64 != 0) { s->status = 7; return 7; }
        }
        for (i = 0; i < m; i++) s->xn[i] = s->x[i] + s->b[i];
        if (walk) {
            long e;
            for (i = 0; i < r; i++) { s->xnf[i] = s->xn[i]; s->x0f[i] = s->x[i]; }
            s->xnf[r] = 0.0;
            s->x0f[r] = 0.0;
            for (i = r; i < m; i++) { s->xnf[i + 1] = s->xn[i]; s->x0f[i + 1] = s->x[i]; }
            e = walk(s->wF, s->wPR, s->wTI, s->T, s->wNM, s->wOFF, s->wK, s->xnf, s->x0f, 0, s->wn);
            if (e < s->wn) { s->status = 8; return 8; }
            for (i = 0, j = 0; i < n; i++) { if (i == r) continue; s->xn[j++] = s->xnf[i]; }
            for (i = 0; i < m; i++) s->xd[i] = s->xn[i] - s->x[i];
        } else {
            for (i = 0; i < m; i++) s->xd[i] = s->b[i];
        }
        for (i = 0; i < m; i++) if (!isfinite(s->xn[i])) { s->status = 9; return 9; }
        for (i = 0; i < m; i++) s->ax[i] = fabs(s->xn[i]);
        gemv(101, 111, m, m, 1.0, s->absJ, m, s->ax, 1, 0.0, s->Is, 1);
        {
            int conv_x = 1, conv_f = 1;
            for (i = 0; i < m; i++) {
                double aF = fabs(s->Fr[i]);
                double Is = s->Is[i] + aF;
                double axi = s->ax[i], axo = fabs(s->x[i]);
                double mx = (axi >= axo) ? axi : axo;
                if (!(fabs(s->xd[i]) < s->reltol * mx + s->xtol[i])) conv_x = 0;
                if (!(aF < s->reltol * Is + s->abstol[i])) conv_f = 0;
            }
            s->iters = it + 1;
            if (conv_x && conv_f) {
                for (i = 0; i < m; i++) s->x[i] = s->xn[i];
                if (s->fold_j) {
                    /* jacobian_only at the converged point: G, C, q -- C and
                       q alone where J is unread (fold_j holds the passes) */
                    for (i = 0; i < r; i++) s->xc[i] = s->x[i];
                    s->xc[r] = 0.0;
                    for (i = r; i < m; i++) s->xc[i + 1] = s->x[i];
                    core(s->core, s->xc, s->T, s->fold_j, s->formula, s->a0, s->a1, s->a2, s->h,
                         s->theta, s->q1, s->q2, s->iq1, s->u, s->Cj, s->Cj, s->qj, s->iqj,
                         s->Geqj, s->F, s->Jj);
                }
                s->status = 1;
                return 1;
            }
        }
        if (s->chord) {
            double size = -INFINITY;
            for (i = 0; i < m; i++) {
                double v = fabs(s->xd[i]);
                if (isnan(v)) { size = v; break; }
                if (v > size) size = v;
            }
            if (!isfinite(size) || (have_prev && !(size < prev))) { s->status = 10; return 10; }
            prev = size;
            have_prev = 1;
        }
        for (i = 0; i < m; i++) s->x[i] = s->xn[i];
    }
    s->status = 11;
    return 11;
}
