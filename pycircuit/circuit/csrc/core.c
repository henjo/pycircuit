
#include <stdint.h>
#include <string.h>
#include "pcx_abi.h"
#include "core.h"
/* passes 0 G, 1 C, 2 i, 3 q; their want bits pcx_pass (pcx_abi.h).  A class
   with a fused kernel (z_fn) is evaluated once a call, before the passes,
   for every pass the call wants, into its staging (zst); its batches then
   scatter from there -- the passes' writes and bins as before.  ztab, a
   group's row: ne, k, ti, mask, nmoff, proff, then per pass G C i q its
   output length and its staging offset; bftab, a batch's row: its group
   (-1: none) and its offset in fidxall (each element's place in it; -1:
   the group's own order, the usual case) */
#define ZW 14
long hdl_fn(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J)
{
    double X[64], o[4096], y[64], xg[64];
    double *bins[4] = {c->binG, C, c->bini, q};
    const long bits[4] = {PCX_PASS_G, PCX_PASS_C, PCX_PASS_I, 0};
    dgemv64_t gemv = (dgemv64_t) c->dgemv;
    long n = c->n, m, b, g, e, j, s, z;
    const double *Cs;
    for (z = 0; z < c->nz; z++) {
        const int64_t *zt = c->ztab + ZW * z;
        fused_fn_t ff = (fused_fn_t) c->z_fn[z];
        long zw = ((want & PCX_WANT_ALL) | PCX_PASS_Q) & zt[3];
        if (!ff || !zw) continue;
        {
            long ne = zt[0], k = zt[1], ti = zt[2];
            const int64_t *NM = c->NMall + zt[4];
            double **PR = c->PRall + zt[5];
            double *sG = c->zst + zt[7], *sC = c->zst + zt[9], *si = c->zst + zt[11];
            double *sq = c->zst + zt[13];
            for (e = 0; e < ne; e++) {
                for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
                if (ti >= 0) PR[e][ti] = T;
                ff(X, PR[e], zw, sG + e * zt[6], sC + e * zt[8], si + e * zt[10], sq + e * zt[12]);
            }
        }
    }
    for (m = 0; m < 4; m++) {
        double *buf = c->work;
        long L = c->p_L[m];
        if (bits[m] && !(want & bits[m])) continue;
        if (m < 2) memcpy(buf, c->tmplall + c->p_toff[m], (size_t) L * sizeof(double));
        else memset(buf, 0, (size_t) L * sizeof(double));
        for (b = 0; b < c->nb; b++) {
            if (c->b_pass[b] != m) continue;
            pass_fn_t f = (pass_fn_t) c->b_fn[b];
            long ne = c->b_ne[b], k = c->b_k[b], so = c->b_so[b], ti = c->b_ti[b];
            long fz = c->nz ? c->bftab[2 * b] : -1;
            const int64_t *NM = c->NMall + c->b_nmoff[b];
            const int64_t *dst = c->bdstall + c->b_dstoff[b];
            double **PR = c->PRall + c->b_proff[b];
            if (fz >= 0 && c->z_fn[fz]) {
                const double *sz = c->zst + c->ztab[ZW * fz + 7 + 2 * m];
                long fo = c->bftab[2 * b + 1];
                const int64_t *fi = fo >= 0 ? c->fidxall + fo : 0;
                for (e = 0; e < ne; e++) {
                    const double *oz = sz + (fi ? fi[e] : e) * so;
                    for (j = 0; j < so; j++) buf[dst[e * so + j]] = oz[j];
                }
                continue;
            }
            for (e = 0; e < ne; e++) {
                for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
                if (ti >= 0) PR[e][ti] = T;
                f(X, PR[e], o);
                for (j = 0; j < so; j++) buf[dst[e * so + j]] = o[j];
            }
        }
        for (g = 0; g < c->ng; g++) {
            if (c->g_pass[g] != m) continue;
            long E = c->g_E[g], k = c->g_k[g], nsel = c->g_nsel[g];
            const double *S = c->Sall + c->g_soff[g];
            const int64_t *I = c->Iall + c->g_ioff[g];
            const int64_t *sel = c->selall + c->g_seloff[g];
            const int64_t *gdst = c->gdstall + c->g_seloff[g];
            for (e = 0; e < E; e++) {
                for (j = 0; j < k; j++) xg[j] = x[I[e * k + j]];
                gemv(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, k, k, 1.0, S + e * k * k, k, xg, 1, 0.0, y, 1);
                for (j = 0; j < k; j++) c->ywork[e * k + j] = y[j];
            }
            for (s = 0; s < nsel; s++) buf[gdst[s]] = c->ywork[sel[s]];
        }
        {
            const int64_t *flat = c->flatall + c->p_toff[m];
            double *bin = bins[m];
            long nbins = c->p_nbins[m];
            memset(bin, 0, (size_t) nbins * sizeof(double));
            for (s = 0; s < L; s++) bin[flat[s]] += buf[s];
        }
    }
    /* formula 4: the passes alone (`passes`: G in binG, i in bini, C and q
       in the caller's buffers) -- no companion, no residual, no Jacobian */
    if (formula == PCX_FORMULA_PASSES) return 0;
    Cs = (want & PCX_PASS_C) ? C : Cin;
    if (formula == PCX_FORMULA_GEAR2) {
        for (j = 0; j < n; j++) iq[j] = (a0 * q[j] + a1 * q1[j]) + a2 * q2[j];
        for (s = 0; s < n * n; s++) Geq[s] = a0 * Cs[s];
    } else if (formula == PCX_FORMULA_EULER) {
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / h;
        for (s = 0; s < n * n; s++) Geq[s] = Cs[s] / h;
    } else if (formula == PCX_FORMULA_TRAPEZOIDAL) {
        for (j = 0; j < n; j++) iq[j] = (2.0 * (q[j] - q1[j])) / h - iq1[j];
        for (s = 0; s < n * n; s++) Geq[s] = (2.0 * Cs[s]) / h;
    } else {  /* PCX_FORMULA_THETA */
        double th = theta * h, w = (1.0 - theta) / theta;
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / th - w * iq1[j];
        for (s = 0; s < n * n; s++) Geq[s] = Cs[s] / th;
    }
    if (want & PCX_PASS_I) for (j = 0; j < n; j++) F[j] = (c->bini[j] + iq[j]) + u[j];
    if (want & PCX_PASS_G) for (s = 0; s < n * n; s++) J[s] = c->binG[s] + Geq[s];
    return 0;
}
