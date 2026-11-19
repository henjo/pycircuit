"""THE SPD REDUCTION ON CHOLMOD: KLU's entries over CHOLMOD's supernodal
Cholesky, for a linear grid's transient matrix.

A circuit of resistors, capacitors, inductors, voltage sources to ground,
current sources and banks (`_bulk`) has a transient matrix J = G +
(alpha / h) C that is exactly symmetric, with a node block K = G_N +
(alpha / h) C_N, each voltage source's branch row coupling its node alone
(J[vb, a] = J[a, vb] = +-1, J[vb, vb] = 0) and each inductor's its two
nodes (J[l, l] = -(alpha / h) L).  Solving J x = r then reduces exactly:

- a voltage source's node is fixed, x_a = r_vb / J[vb, a];
- an inductor's current is i_l = (r_l - sum_n J[l, n] x_n) / J[l, l];
- the free nodes solve K' x_F = r', K'[f, m] = J[f, m] + sum_l g_l
  J[f, l] J[l, m] with g_l = -1 / J[l, l] > 0, r'_f = r_f + sum_l g_l
  J[f, l] r_l - sum_(m fixed) K'[f, m] x_m -- K' symmetric positive
  definite wherever every free node reaches ground or a fixed node;
- each voltage source's current from its node's row,
  i_vb = (r_a - sum_m J[a, m] x_m - sum_l J[a, l] i_l) / J[a, vb].

The C below does it on CHOLMOD's supernodal LL' (the ordering CHOLMOD
chooses, once a pattern) behind KLU's own signatures -- `klu_factor`,
`klu_refactor`, `klu_solve`, `klu_tsolve` (J symmetric: the solve),
`klu_free_numeric`, `klu_free_symbolic`, `klu_rgrowth` (1.0 at the common's
`rgrowth`: a Cholesky of a positive definite matrix has no pivot order to
go stale) -- so a `StructuralKLU` holding it (`linearsolver`, its
`cholmod` plan) is the same solver to the Python and to every C unit that
takes KLU's entries from it.  A factor CHOLMOD refuses (not positive
definite, an inductor's diagonal zero, a non-finite value) is KLU's
singular one: NULL, or 0 from a refactor.

`plan(cir, iref)`: the rows' kinds -- 0 a node, 1 a voltage source's
branch, 2 an inductor's -- of the reduced unknowns where every leaf
element is one of the kinds above, a voltage source with one end on the
reference and no node fixed twice; else the reason (a str).  Switch:
`ENABLED` (env `PYCIRCUIT_CHOLMOD=0`); taken by `AutoSolver` for a
transient's matrix from `CHOLMOD_FROM` unknowns.
"""
import ctypes
import os

import numpy as np

from pycircuit.circuit import _blas, _cbuild

ENABLED = os.environ.get('PYCIRCUIT_CHOLMOD', '1') != '0'
#: the unknowns from which a linear grid's transient matrix is CHOLMOD's
CHOLMOD_FROM = 500_000
STATUS = 'not loaded'

_SRC = r'''#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dlfcn.h>
#include <suitesparse/cholmod.h>

typedef int (*ck_start_t)(cholmod_common *);
typedef cholmod_sparse *(*ck_asp_t)(size_t, size_t, size_t, int, int, int, int, cholmod_common *);
typedef int (*ck_fsp_t)(cholmod_sparse **, cholmod_common *);
typedef cholmod_factor *(*ck_an_t)(cholmod_sparse *, cholmod_common *);
typedef int (*ck_fz_t)(cholmod_sparse *, cholmod_factor *, cholmod_common *);
typedef cholmod_factor *(*ck_cpf_t)(cholmod_factor *, cholmod_common *);
typedef int (*ck_ff_t)(cholmod_factor **, cholmod_common *);
typedef cholmod_dense *(*ck_ad_t)(size_t, size_t, size_t, int, cholmod_common *);
typedef int (*ck_fd_t)(cholmod_dense **, cholmod_common *);
typedef cholmod_dense *(*ck_sv_t)(int, cholmod_factor *, cholmod_dense *, cholmod_common *);

static struct {
    int tried, ok;
    ck_start_t start, finish;
    ck_asp_t asp; ck_fsp_t fsp; ck_an_t an; ck_fz_t fz; ck_cpf_t cpf; ck_ff_t ff;
    ck_ad_t ad; ck_fd_t fd; ck_sv_t sv;
} C_;

static int ck_load(void)
{
    void *h;
    if (C_.tried) return C_.ok;
    C_.tried = 1;
    h = dlopen("libcholmod.so.5", RTLD_NOW | RTLD_LOCAL);
    if (!h) h = dlopen("libcholmod.so", RTLD_NOW | RTLD_LOCAL);
    if (!h) return 0;
    C_.start = (ck_start_t) dlsym(h, "cholmod_start");
    C_.finish = (ck_start_t) dlsym(h, "cholmod_finish");
    C_.asp = (ck_asp_t) dlsym(h, "cholmod_allocate_sparse");
    C_.fsp = (ck_fsp_t) dlsym(h, "cholmod_free_sparse");
    C_.an = (ck_an_t) dlsym(h, "cholmod_analyze");
    C_.fz = (ck_fz_t) dlsym(h, "cholmod_factorize");
    C_.cpf = (ck_cpf_t) dlsym(h, "cholmod_copy_factor");
    C_.ff = (ck_ff_t) dlsym(h, "cholmod_free_factor");
    C_.ad = (ck_ad_t) dlsym(h, "cholmod_allocate_dense");
    C_.fd = (ck_fd_t) dlsym(h, "cholmod_free_dense");
    C_.sv = (ck_sv_t) dlsym(h, "cholmod_solve");
    C_.ok = C_.start && C_.finish && C_.asp && C_.fsp && C_.an && C_.fz && C_.cpf && C_.ff
            && C_.ad && C_.fd && C_.sv;
    return C_.ok;
}

/* the symbolic: the rows' kinds, the reduction's maps, K' (the free nodes' upper
   triangle, CSC) and CHOLMOD's analysis of it */
typedef struct {
    long m, nnz, nf, nv, nl, nd, nk, rg_off;
    int32_t *Ap, *Ai, *kind, *fidx;
    int32_t *v_row, *v_node, *l_row, *l_n1, *l_n2;
    int64_t *v_pav, *v_pva, *l_pll, *l_p1, *l_p2, *l_k11, *l_k22, *l_k12;
    int64_t *d_src, *d_dst;
    int32_t *Kp, *Ki;
    cholmod_common cm;
    cholmod_factor *Lsym;
} ck_sym;

typedef struct {
    double *Ax, *Kx, *x, *rhs;
    cholmod_common cm;
    cholmod_sparse *K;
    cholmod_factor *L;
    cholmod_dense *B;
} ck_num;

static void *ck_alloc(size_t n, size_t s) { return calloc(n ? n : 1, s); }

/* the position of row r in column c of (Ap, Ai) -- sorted or not -- or -1 */
static int64_t ck_at(const int32_t *Ap, const int32_t *Ai, long c, long r)
{
    int64_t p;
    for (p = Ap[c]; p < Ap[c + 1]; p++) if (Ai[p] == r) return p;
    return -1;
}

static int ck_cmp64(const void *a, const void *b)
{
    int64_t x = *(const int64_t *) a, y = *(const int64_t *) b;
    return (x > y) - (x < y);
}

/* K' position of the free pair (i, j), i <= j, or -1 */
static int64_t ck_kpos(const ck_sym *s, long i, long j)
{
    long lo = s->Kp[j], hi = s->Kp[j + 1] - 1;
    while (lo <= hi) {
        long mid = (lo + hi) / 2;
        if (s->Ki[mid] == i) return mid;
        if (s->Ki[mid] < i) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

static int ck_free_symbolic(void **ps, void *common)
{
    ck_sym *s = (ck_sym *) *ps;
    (void) common;
    if (!s) return 1;
    if (s->Lsym) C_.ff(&s->Lsym, &s->cm);
    C_.finish(&s->cm);
    free(s->Ap); free(s->Ai); free(s->kind); free(s->fidx);
    free(s->v_row); free(s->v_node); free(s->l_row); free(s->l_n1); free(s->l_n2);
    free(s->v_pav); free(s->v_pva); free(s->l_pll); free(s->l_p1); free(s->l_p2);
    free(s->l_k11); free(s->l_k22); free(s->l_k12); free(s->d_src); free(s->d_dst);
    free(s->Kp); free(s->Ki);
    free(s);
    *ps = NULL;
    return 1;
}

/* the reduction's classification of the rows `kind` (0 node, 1 V branch, 2 L
   branch) checked against the pattern, its maps, K' and CHOLMOD's analysis;
   NULL where the pattern is not a linear grid's (a V coupling more than its
   one node or another branch, a node fixed twice, an L touching another
   branch or more than two nodes) or CHOLMOD is missing */
static void *ck_analyze(int m_, const int32_t *Ap, const int32_t *Ai, const int32_t *kind,
                        long rg_off)
{
    long m = m_, c, r, k, nnz = Ap[m_], nf = 0, nv = 0, nl = 0, npair = 0, j;
    int64_t p, *pairs = NULL, *key = NULL;
    ck_sym *s;
    if (!ck_load()) return NULL;
    s = (ck_sym *) calloc(1, sizeof(ck_sym));
    if (!s) return NULL;
    C_.start(&s->cm);
    s->m = m; s->nnz = nnz; s->rg_off = rg_off;
    s->Ap = (int32_t *) ck_alloc(m + 1, sizeof(int32_t));
    s->Ai = (int32_t *) ck_alloc(nnz, sizeof(int32_t));
    s->kind = (int32_t *) ck_alloc(m, sizeof(int32_t));
    s->fidx = (int32_t *) ck_alloc(m, sizeof(int32_t));
    if (!s->Ap || !s->Ai || !s->kind || !s->fidx) goto fail;
    memcpy(s->Ap, Ap, (m + 1) * sizeof(int32_t));
    memcpy(s->Ai, Ai, nnz * sizeof(int32_t));
    for (r = 0; r < m; r++) {
        s->kind[r] = kind[r];
        s->fidx[r] = -1;
        if (kind[r] == 1) nv++;
        else if (kind[r] == 2) nl++;
        else if (kind[r] != 0) goto fail;
    }
    s->nv = nv; s->nl = nl;
    s->v_row = (int32_t *) ck_alloc(nv, sizeof(int32_t));
    s->v_node = (int32_t *) ck_alloc(nv, sizeof(int32_t));
    s->v_pav = (int64_t *) ck_alloc(nv, sizeof(int64_t));
    s->v_pva = (int64_t *) ck_alloc(nv, sizeof(int64_t));
    s->l_row = (int32_t *) ck_alloc(nl, sizeof(int32_t));
    s->l_n1 = (int32_t *) ck_alloc(nl, sizeof(int32_t));
    s->l_n2 = (int32_t *) ck_alloc(nl, sizeof(int32_t));
    s->l_pll = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    s->l_p1 = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    s->l_p2 = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    s->l_k11 = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    s->l_k22 = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    s->l_k12 = (int64_t *) ck_alloc(nl, sizeof(int64_t));
    if (!s->v_row || !s->v_node || !s->v_pav || !s->v_pva || !s->l_row || !s->l_n1 || !s->l_n2
        || !s->l_pll || !s->l_p1 || !s->l_p2 || !s->l_k11 || !s->l_k22 || !s->l_k12) goto fail;
    /* each V: its one node (fixed), its two incidences */
    nv = 0;
    for (c = 0; c < m; c++) {
        long node = -1;
        if (s->kind[c] != 1) continue;
        for (p = Ap[c]; p < Ap[c + 1]; p++) {
            r = Ai[p];
            if (r == c) continue;
            if (s->kind[r] != 0 || node >= 0) goto fail;
            node = r;
        }
        if (node < 0 || s->fidx[node] == -2) goto fail;
        s->fidx[node] = -2;                       /* fixed */
        s->v_row[nv] = (int32_t) c; s->v_node[nv] = (int32_t) node;
        s->v_pav[nv] = ck_at(Ap, Ai, c, node);
        s->v_pva[nv] = ck_at(Ap, Ai, node, c);
        if (s->v_pav[nv] < 0 || s->v_pva[nv] < 0) goto fail;
        nv++;
    }
    for (r = 0; r < m; r++) if (s->kind[r] == 0 && s->fidx[r] == -1) s->fidx[r] = (int32_t) nf++;
    s->nf = nf;
    /* each L: its diagonal and its (at most two) nodes */
    nl = 0;
    for (c = 0; c < m; c++) {
        long n1 = -1, n2 = -1;
        if (s->kind[c] != 2) continue;
        s->l_pll[nl] = -1;
        for (p = Ap[c]; p < Ap[c + 1]; p++) {
            r = Ai[p];
            if (r == c) { s->l_pll[nl] = p; continue; }
            if (s->kind[r] != 0) goto fail;
            if (n1 < 0) n1 = r; else if (n2 < 0) n2 = r; else goto fail;
        }
        if (s->l_pll[nl] < 0) goto fail;
        s->l_row[nl] = (int32_t) c; s->l_n1[nl] = (int32_t) n1; s->l_n2[nl] = (int32_t) n2;
        s->l_p1[nl] = n1 >= 0 ? ck_at(Ap, Ai, c, n1) : -1;
        s->l_p2[nl] = n2 >= 0 ? ck_at(Ap, Ai, c, n2) : -1;
        nl++;
    }
    /* a node column coupling a V row that is not its own: not a grid's */
    for (c = 0; c < m; c++) {
        if (s->kind[c] != 0) continue;
        for (p = Ap[c]; p < Ap[c + 1]; p++) {
            r = Ai[p];
            if (s->kind[r] == 1 && s->fidx[c] != -2) goto fail;
        }
    }
    /* K' pattern: the free node block's upper triangle and each L's free pairs */
    pairs = (int64_t *) ck_alloc(nnz + 3 * nl, sizeof(int64_t));
    if (!pairs) goto fail;
    for (c = 0; c < m; c++) {
        long fc = s->fidx[c];
        if (s->kind[c] != 0 || fc < 0) continue;
        for (p = Ap[c]; p < Ap[c + 1]; p++) {
            long fr;
            r = Ai[p];
            if (s->kind[r] != 0) continue;
            fr = s->fidx[r];
            if (fr < 0 || fr > fc) continue;
            pairs[npair++] = (int64_t) fc * nf + fr;
        }
    }
    for (k = 0; k < nl; k++) {
        long a = s->l_n1[k] >= 0 ? s->fidx[s->l_n1[k]] : -1;
        long b = s->l_n2[k] >= 0 ? s->fidx[s->l_n2[k]] : -1;
        if (a >= 0) pairs[npair++] = (int64_t) a * nf + a;
        if (b >= 0) pairs[npair++] = (int64_t) b * nf + b;
        if (a >= 0 && b >= 0) pairs[npair++] = (int64_t) (a > b ? a : b) * nf + (a < b ? a : b);
    }
    qsort(pairs, (size_t) npair, sizeof(int64_t), ck_cmp64);
    s->Kp = (int32_t *) ck_alloc(nf + 1, sizeof(int32_t));
    s->Ki = (int32_t *) ck_alloc(npair, sizeof(int32_t));
    if (!s->Kp || !s->Ki) goto fail;
    s->nk = 0;
    for (p = 0; p < npair; p++) {
        if (p && pairs[p] == pairs[p - 1]) continue;
        s->Ki[s->nk] = (int32_t) (pairs[p] % nf);
        s->Kp[pairs[p] / nf + 1]++;
        s->nk++;
    }
    for (j = 0; j < nf; j++) s->Kp[j + 1] += s->Kp[j];
    free(pairs); pairs = NULL;
    /* the direct map: each free upper node entry of J to its K' place */
    key = (int64_t *) ck_alloc(2 * (size_t) nnz, sizeof(int64_t));
    if (!key) goto fail;
    s->nd = 0;
    for (c = 0; c < m; c++) {
        long fc = s->fidx[c];
        if (s->kind[c] != 0 || fc < 0) continue;
        for (p = Ap[c]; p < Ap[c + 1]; p++) {
            long fr;
            r = Ai[p];
            if (s->kind[r] != 0) continue;
            fr = s->fidx[r];
            if (fr < 0 || fr > fc) continue;
            key[2 * s->nd] = p;
            key[2 * s->nd + 1] = ck_kpos(s, fr, fc);
            s->nd++;
        }
    }
    s->d_src = (int64_t *) ck_alloc(s->nd, sizeof(int64_t));
    s->d_dst = (int64_t *) ck_alloc(s->nd, sizeof(int64_t));
    if (!s->d_src || !s->d_dst) goto fail;
    for (k = 0; k < s->nd; k++) { s->d_src[k] = key[2 * k]; s->d_dst[k] = key[2 * k + 1]; }
    free(key); key = NULL;
    for (k = 0; k < nl; k++) {
        long a = s->l_n1[k] >= 0 ? s->fidx[s->l_n1[k]] : -1;
        long b = s->l_n2[k] >= 0 ? s->fidx[s->l_n2[k]] : -1;
        s->l_k11[k] = a >= 0 ? ck_kpos(s, a, a) : -1;
        s->l_k22[k] = b >= 0 ? ck_kpos(s, b, b) : -1;
        s->l_k12[k] = (a >= 0 && b >= 0) ? ck_kpos(s, a < b ? a : b, a > b ? a : b) : -1;
    }
    /* CHOLMOD's analysis of K' (its pattern; the values come with a factor) */
    {
        cholmod_sparse *K;
        s->cm.supernodal = CHOLMOD_SUPERNODAL;
        K = C_.asp((size_t) nf, (size_t) nf, (size_t) s->nk, 1, 1, 1, CHOLMOD_REAL, &s->cm);
        if (!K) goto fail;
        memcpy(K->p, s->Kp, (size_t) (nf + 1) * sizeof(int32_t));
        memcpy(K->i, s->Ki, (size_t) s->nk * sizeof(int32_t));
        for (p = 0; p < s->nk; p++) ((double *) K->x)[p] = 1.0;
        s->Lsym = C_.an(K, &s->cm);
        C_.fsp(&K, &s->cm);
        if (!s->Lsym) goto fail;
    }
    return s;
fail:
    free(pairs); free(key);
    {
        void *v = s;
        ck_free_symbolic(&v, NULL);
    }
    return NULL;
}

static int ck_free_numeric(void **pn, void *common)
{
    ck_num *n = (ck_num *) *pn;
    (void) common;
    if (!n) return 1;
    if (n->B) C_.fd(&n->B, &n->cm);
    if (n->L) C_.ff(&n->L, &n->cm);
    if (n->K) C_.fsp(&n->K, &n->cm);
    C_.finish(&n->cm);
    free(n->Ax); free(n->x); free(n->rhs);
    free(n);
    *pn = NULL;
    return 1;
}

/* K' values from J's (`Ax`, kept on the numeric) and CHOLMOD's factor of them:
   1, or 0 -- a non-finite value, an inductor's diagonal not negative, K' not
   positive definite */
static int ck_values(const ck_sym *s, ck_num *n, const double *Ax)
{
    long k;
    double *Kx = (double *) n->K->x;
    for (k = 0; k < s->nnz; k++) {
        if (!isfinite(Ax[k])) return 0;
        n->Ax[k] = Ax[k];
    }
    memset(Kx, 0, (size_t) s->nk * sizeof(double));
    for (k = 0; k < s->nd; k++) Kx[s->d_dst[k]] += Ax[s->d_src[k]];
    for (k = 0; k < s->nl; k++) {
        double dll = Ax[s->l_pll[k]], g, a1, a2;
        if (!(dll < 0.0)) return 0;
        g = -1.0 / dll;
        a1 = s->l_p1[k] >= 0 ? Ax[s->l_p1[k]] : 0.0;
        a2 = s->l_p2[k] >= 0 ? Ax[s->l_p2[k]] : 0.0;
        if (s->l_k11[k] >= 0) Kx[s->l_k11[k]] += g * a1 * a1;
        if (s->l_k22[k] >= 0) Kx[s->l_k22[k]] += g * a2 * a2;
        if (s->l_k12[k] >= 0) Kx[s->l_k12[k]] += g * a1 * a2;
    }
    C_.fz(n->K, n->L, &n->cm);
    return n->cm.status == CHOLMOD_OK && (long) n->L->minor >= s->nf;
}

static void *ck_factor(int *Ap, int *Ai, double *Ax, void *sym, void *common)
{
    ck_sym *s = (ck_sym *) sym;
    ck_num *n;
    void *v;
    (void) Ap; (void) Ai;
    if (!s) return NULL;
    n = (ck_num *) calloc(1, sizeof(ck_num));
    if (!n) return NULL;
    C_.start(&n->cm);
    n->cm.supernodal = CHOLMOD_SUPERNODAL;
    n->Ax = (double *) ck_alloc(s->nnz, sizeof(double));
    n->x = (double *) ck_alloc(s->m, sizeof(double));
    n->rhs = (double *) ck_alloc(s->nf, sizeof(double));
    n->K = C_.asp((size_t) s->nf, (size_t) s->nf, (size_t) s->nk, 1, 1, 1, CHOLMOD_REAL, &n->cm);
    n->L = C_.cpf(s->Lsym, &n->cm);
    n->B = C_.ad((size_t) s->nf, 1, (size_t) s->nf, CHOLMOD_REAL, &n->cm);
    if (n->Ax && n->x && n->rhs && n->K && n->L && n->B) {
        memcpy(n->K->p, s->Kp, (size_t) (s->nf + 1) * sizeof(int32_t));
        memcpy(n->K->i, s->Ki, (size_t) s->nk * sizeof(int32_t));
        if (ck_values(s, n, Ax)) return n;
    }
    v = n;
    ck_free_numeric(&v, common);
    return NULL;
}

static int ck_refactor(int *Ap, int *Ai, double *Ax, void *sym, void *num, void *common)
{
    (void) Ap; (void) Ai; (void) common;
    if (!sym || !num) return 0;
    return ck_values((const ck_sym *) sym, (ck_num *) num, Ax);
}

/* J x = b for each of the `nrhs` columns of `B` (n x nrhs, column-major), in place */
static int ck_solve(void *sym, void *num, int n_, int nrhs, double *B, void *common)
{
    const ck_sym *s = (const ck_sym *) sym;
    ck_num *nm = (ck_num *) num;
    const double *Ax;
    double *x, *rhs;
    long j, k, r, c;
    int64_t p;
    (void) common;
    if (!s || !nm || n_ != s->m) return 0;
    Ax = nm->Ax; x = nm->x; rhs = nm->rhs;
    for (j = 0; j < nrhs; j++) {
        double *b = B + (size_t) j * s->m;
        cholmod_dense *X;
        memset(x, 0, (size_t) s->m * sizeof(double));
        for (k = 0; k < s->nv; k++) x[s->v_node[k]] = b[s->v_row[k]] / Ax[s->v_pva[k]];
        for (r = 0; r < s->m; r++) if (s->kind[r] == 0 && s->fidx[r] >= 0) rhs[s->fidx[r]] = b[r];
        for (k = 0; k < s->nl; k++) {
            double g = -1.0 / Ax[s->l_pll[k]], bl = b[s->l_row[k]];
            long n1 = s->l_n1[k], n2 = s->l_n2[k];
            if (n1 >= 0 && s->fidx[n1] >= 0) rhs[s->fidx[n1]] += g * Ax[s->l_p1[k]] * bl;
            if (n2 >= 0 && s->fidx[n2] >= 0) rhs[s->fidx[n2]] += g * Ax[s->l_p2[k]] * bl;
            /* (an L between a free and a fixed node: its term of K'[f, fixed] x_fixed) */
            if (n1 >= 0 && n2 >= 0) {
                if (s->fidx[n1] >= 0 && s->fidx[n2] == -2)
                    rhs[s->fidx[n1]] -= g * Ax[s->l_p1[k]] * Ax[s->l_p2[k]] * x[n2];
                if (s->fidx[n2] >= 0 && s->fidx[n1] == -2)
                    rhs[s->fidx[n2]] -= g * Ax[s->l_p2[k]] * Ax[s->l_p1[k]] * x[n1];
            }
        }
        for (k = 0; k < s->nv; k++) {
            c = s->v_node[k];
            for (p = s->Ap[c]; p < s->Ap[c + 1]; p++) {
                r = s->Ai[p];
                if (s->kind[r] == 0 && s->fidx[r] >= 0) rhs[s->fidx[r]] -= Ax[p] * x[c];
            }
        }
        memcpy(nm->B->x, rhs, (size_t) s->nf * sizeof(double));
        X = C_.sv(CHOLMOD_A, nm->L, nm->B, &nm->cm);
        if (!X) return 0;
        for (r = 0; r < s->m; r++)
            if (s->kind[r] == 0 && s->fidx[r] >= 0) x[r] = ((double *) X->x)[s->fidx[r]];
        C_.fd(&X, &nm->cm);
        for (k = 0; k < s->nl; k++) {
            double v = b[s->l_row[k]];
            if (s->l_n1[k] >= 0) v -= Ax[s->l_p1[k]] * x[s->l_n1[k]];
            if (s->l_n2[k] >= 0) v -= Ax[s->l_p2[k]] * x[s->l_n2[k]];
            x[s->l_row[k]] = v / Ax[s->l_pll[k]];
        }
        for (k = 0; k < s->nv; k++) {
            double v;
            c = s->v_node[k];
            v = b[c];
            for (p = s->Ap[c]; p < s->Ap[c + 1]; p++) {
                r = s->Ai[p];
                if (r != s->v_row[k] && s->kind[r] != 1) v -= Ax[p] * x[r];
            }
            x[s->v_row[k]] = v / Ax[s->v_pav[k]];
        }
        memcpy(b, x, (size_t) s->m * sizeof(double));
    }
    return 1;
}

static int ck_rgrowth(int *Ap, int *Ai, double *Ax, void *sym, void *num, void *common)
{
    const ck_sym *s = (const ck_sym *) sym;
    double one = 1.0;
    (void) Ap; (void) Ai; (void) Ax;
    if (!s || !num || s->rg_off < 0) return 0;
    memcpy((char *) common + s->rg_off, &one, sizeof one);
    return 1;
}

/* the factor's statistics: [the free nodes, K's entries, L's entries, the ordering CHOLMOD chose] */
static void ck_info(void *sym, double *out)
{
    const ck_sym *s = (const ck_sym *) sym;
    out[0] = (double) s->nf; out[1] = (double) s->nk;
    out[2] = s->Lsym ? (double) s->cm.lnz : -1.0;
    out[3] = (double) s->cm.selected;
}

void *pcx_cholklu(int which)
{
    switch (which) {
    case 0: return (void *) ck_analyze;
    case 1: return (void *) ck_factor;
    case 2: return (void *) ck_refactor;
    case 3: return (void *) ck_solve;
    case 4: return (void *) ck_free_numeric;
    case 5: return (void *) ck_free_symbolic;
    case 6: return (void *) ck_rgrowth;
    case 7: return (void *) ck_info;
    case 8: return (void *) (long) ck_load();
    }
    return NULL;
}
'''

_CDEF = 'void *pcx_cholklu(int which);'
_LIB = []


class _Lib:
    """KLU's entries over the reduction (`klu_factor`, `klu_refactor`,
    `klu_solve`, `klu_tsolve`, `klu_free_numeric`, `klu_free_symbolic`,
    `klu_rgrowth`, KLU's prototypes), `klu_defaults` the KLU library's own,
    `analyze(m, Ap, Ai, kind, rg_off)` the symbolic and `info(sym)` its
    statistics."""

    def __init__(self, ffi, entry, klu):
        ct = ctypes
        vp, ip, dp = ct.c_void_p, ct.POINTER(ct.c_int), ct.POINTER(ct.c_double)

        def fn(k, res, *args):
            return ct.CFUNCTYPE(res, *args)(int(ffi.cast('uintptr_t', entry(k))))
        self._keep = (ffi, entry)
        self.analyze = fn(0, vp, ct.c_int, ip, ip, ip, ct.c_long)
        self.klu_factor = fn(1, vp, ip, ip, dp, vp, vp)
        self.klu_refactor = fn(2, ct.c_int, ip, ip, dp, vp, vp, vp)
        self.klu_solve = fn(3, ct.c_int, vp, vp, ct.c_int, ct.c_int, dp, vp)
        self.klu_tsolve = self.klu_solve
        self.klu_free_numeric = fn(4, ct.c_int, vp, vp)
        self.klu_free_symbolic = fn(5, ct.c_int, vp, vp)
        self.klu_rgrowth = fn(6, ct.c_int, ip, ip, dp, vp, vp, vp)
        self.info = fn(7, None, vp, dp)
        self.klu_defaults = klu.klu_defaults


def library(klu):
    """The `_Lib` (built once a process), or None (`STATUS`): no CHOLMOD, no
    compiler, no header."""
    global STATUS
    if _LIB:
        return _LIB[0]
    try:
        ffi, entry, _key, _cold, _secs = _cbuild.load_unit(_SRC, _CDEF, 'pcx_cholklu')
    except Exception as e:  # noqa: BLE001 -- no CHOLMOD header or compiler: KLU stays
        STATUS = f'off ({type(e).__name__})'
        return None
    ## CHOLMOD brings a BLAS of its own (the system OpenBLAS, which reads
    ## OPENBLAS_NUM_THREADS as it loads), started at one thread unless the
    ## environment says: started on every core and held at one, CHOLMOD's
    ## solve on ibmpg3t ran 6.5 % slower, the same bits
    unset = 'OPENBLAS_NUM_THREADS' not in os.environ
    if unset:
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
    try:
        loaded = int(ffi.cast('long', entry(8)))
    finally:
        if unset:
            os.environ.pop('OPENBLAS_NUM_THREADS', None)
    ## (the BLAS hold must see that library, or a held solve leaves it
    ## threaded)
    _blas.library_loaded()
    if not loaded:
        STATUS = 'off (no libcholmod)'
        return None
    _LIB.append(_Lib(ffi, entry, klu))
    STATUS = 'on'
    return _LIB[0]


#: the element classes whose stamps keep a grid's matrix the reduction's
_GRID = {}


def _grid_kinds():
    if not _GRID:
        from pycircuit.circuit import elements
        _GRID['plain'] = (elements.R, elements.C)
        _GRID['L'] = elements.L
        _GRID['VS'], _GRID['IS'] = elements.VS, elements.IS
    return _GRID


def plan(cir, iref):
    """The reduced unknowns' kinds (int32: 0 a node, 1 a voltage source's
    branch, 2 an inductor's) for `cir`'s transient matrix (its reference row
    `iref` removed), or why not (a str): every leaf element a resistor,
    capacitor, inductor, current source or bank, or a voltage source with one
    end on the reference -- its classes' own G and C, none overridden on
    the instance."""
    g = _grid_kinds()
    view = cir._flat_view() if hasattr(cir, '_flat_view') else cir
    if view is not cir:
        return 'hierarchy'
    view._ensure_node_map()
    kind = np.zeros(view.n, dtype=np.int32)
    fixed = set()
    for inst, el in view.elements.items():
        cls = type(el)
        if 'G' in el.__dict__ or 'C' in el.__dict__:
            return 'instance'
        nm = [int(i) for i in view.elementnodemap[inst]]
        if cls in g['plain'] and len(nm) == 2:
            continue
        if issubclass(cls, g['IS']) and cls.G is g['IS'].G and cls.C is g['IS'].C:
            continue
        if issubclass(cls, g['L']) and cls.G is g['L'].G and cls.C is g['L'].C and len(nm) == 3:
            kind[nm[2]] = 2
            continue
        if issubclass(cls, g['VS']) and cls.G is g['VS'].G and cls.C is g['VS'].C and len(nm) == 3:
            ends = [i for i in nm[:2] if i != iref]
            if len(ends) != 1 or ends[0] in fixed:
                return 'vsource'
            fixed.add(ends[0])
            kind[nm[2]] = 1
            continue
        return 'element'
    return np.ascontiguousarray(np.delete(kind, iref) if iref is not None else kind)
