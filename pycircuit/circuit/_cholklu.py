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

THE KEPT FACTORS.  A numeric may keep up to `keep - 1` earlier factors
beside its own (`ck_keep` on the symbolic: `StructuralKLU.keep_factors`;
`KEEP`, `KEEP_BYTES` the defaults a solver starts with -- 1, none kept).
A refactor of values one of them was made on swaps it in -- L's values and
J's, by pointer -- and factors nothing: CHOLMOD's factor of the same values
on the same symbolic is the same bytes, so every caller (the Python's
`solve`, the C units through KLU's entries, a hand-back's restore) has the
answer a refactor gives.  Any other refactor keeps the numeric's factor in
a slot (a new one while `keep - 1`, `KEEP_BYTES` and THE BUDGET allow, else
the least recently used) and factors into the buffers that come back.  A slot is L's
values (`xsize` doubles, CHOLMOD's allocator) and J's -- hundreds of MB on
a grid of a million unknowns -- for the factors a step-size ladder
revisits around every corner (plan 261 L24).  Counted on the numeric
(`num_info`: `StructuralKLU.kept_counts`).

THE BUDGET.  Where a limit is set (`keep_factors`' `limit`; the step
ladder's is `memory_limit()`: the least of the `memory.max` up this
process's cgroup tree and its resident memory plus the machine's
MemAvailable, read as the run starts), a new slot is made only while the
process's resident memory (/proc/self/statm, read at that slot), the slot
and `KEEP_MARGIN` fit under it -- the bytes of `KEEP_BYTES` then
unbounded.  A ladder run makes its slots in its first steps, one matrix a
step as the sizes double up from the first; on the IBM grids (10 ns,
eleven slots) the run grows by 0.53 GiB (ibmpg3t) and 0.58 GiB (ibmpg6t)
after its last slot, so the margin leaves about 0.9 GiB.  The answers are
a refactor's however many slots the budget allows.

THE OWN SOLVE.  K' x = b on CHOLMOD's supernodal factor L L' = P K' P' is
solved here, not by `cholmod_solve`, which calls the BLAS for every
supernode -- twice each way -- where a grid's factor has hundreds of
thousands of them, most a few columns wide: a supernode below
`CK_BLAS_FROM` columns in plain loops (its rows below the triangle four
columns at a time), one of `CK_BLAS_FROM` columns or more by `dtrsv` and
`dgemv` of the BLAS CHOLMOD loads, as `cholmod_solve` calls them.  The
answer is `cholmod_solve`'s to rounding (another order of the sums), not
its bytes.  Switch: `OWN_SOLVE` (env `PYCIRCUIT_CHOLMOD_OWN_SOLVE=0`:
`cholmod_solve`), set on the symbolic at its analysis (`own`); a factor
the loops do not read (not supernodal, not int32, no BLAS entries found)
is `cholmod_solve`'s.  Counted on the numeric (`num_info`:
`StructuralKLU.solve_counts`); the factor's supernodes, and those of
`CK_BLAS_FROM` columns or more, in `info`.  A one-column supernode is taken
inline, and the backward pass's rows-below sums eight columns at a time:
each column's sum the same operations in the same order.

THE THREADS.  The own solve on a factor of `THREADS_FROM` values or more
runs on the C units' pool (`_tran_core.pool_address`, the count
`set_threads`'s) where it has two threads or more -- the bytes of the solve
on one.  The supernodal tree is cut into subtrees (32 a thread, none above
twice its share, or a leaf) below a top of their ancestors (a plan per
symbolic and count, made at its first solve on threads).  The backward pass
makes the top, then the subtrees on the threads: a supernode reads only its
ancestors' rows and writes only its own, so any order after its parent is
the serial one's bytes.  The forward pass makes the subtrees on the threads
first, each keeping its updates to the rows above it in its own order, then
in the serial order of positions each subtree's kept updates and each top
supernode -- every row taking its updates as the serial pass gives them.
Switch: `THREADS` (env `PYCIRCUIT_CHOLMOD_THREADS=0`: the solve on the
caller alone).  Counted on the numeric (`StructuralKLU.solve_counts`'
`threaded`: the own solves whose passes the pool took).
"""
import ctypes
import os

import numpy as np

from pycircuit.circuit import _blas, _cbuild

ENABLED = os.environ.get('PYCIRCUIT_CHOLMOD', '1') != '0'
#: THE OWN SOLVE (the module's note) in place of `cholmod_solve`
OWN_SOLVE = os.environ.get('PYCIRCUIT_CHOLMOD_OWN_SOLVE', '1') != '0'
#: the unknowns from which a linear grid's transient matrix is CHOLMOD's
CHOLMOD_FROM = 500_000
#: THE THREADS (the module's note): the own solve on the C units' pool, from
#: a factor of `THREADS_FROM` values
THREADS = os.environ.get('PYCIRCUIT_CHOLMOD_THREADS', '1') != '0'
THREADS_FROM = 1 << 20
#: THE KEPT FACTORS a solver starts with: the factors a numeric holds, its
#: own among them (1: none kept), and the most bytes its kept ones may take
#: where no limit is set (THE BUDGET)
KEEP = 1
KEEP_BYTES = 3.5 * 2 ** 30
#: THE BUDGET's margin: the bytes a new slot leaves under the limit
KEEP_MARGIN = 1.5 * 2 ** 30
STATUS = 'not loaded'

_SRC = r'''#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <unistd.h>
#include <dlfcn.h>
#include <suitesparse/cholmod.h>

/* the most factors a numeric keeps, and the entries of J's values a
   comparison reads before the whole */
#define CK_KEEP_MAX 32
#define CK_NSIG 64

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
typedef void *(*ck_mal_t)(size_t, size_t, cholmod_common *);
typedef void *(*ck_fre_t)(size_t, size_t, void *, cholmod_common *);
/* the BLAS CHOLMOD loads (Fortran: the characters' lengths last) */
typedef void (*ck_trsv_t)(const char *, const char *, const char *, const int *, const double *,
                          const int *, double *, const int *, size_t, size_t, size_t);
typedef void (*ck_gemv_t)(const char *, const int *, const int *, const double *, const double *,
                          const int *, const double *, const int *, const double *, double *,
                          const int *, size_t);

static struct {
    int tried, ok;
    ck_start_t start, finish;
    ck_asp_t asp; ck_fsp_t fsp; ck_an_t an; ck_fz_t fz; ck_cpf_t cpf; ck_ff_t ff;
    ck_ad_t ad; ck_fd_t fd; ck_sv_t sv;
    /* CHOLMOD's own allocator, for the kept factors' values (NULL: none kept) */
    ck_mal_t mal; ck_fre_t fre;
    /* THE OWN SOLVE's BLAS entries, from CHOLMOD's (NULL: cholmod_solve) */
    ck_trsv_t trsv; ck_gemv_t gemv;
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
    C_.mal = (ck_mal_t) dlsym(h, "cholmod_malloc");
    C_.fre = (ck_fre_t) dlsym(h, "cholmod_free");
    C_.trsv = (ck_trsv_t) dlsym(h, "dtrsv_");
    C_.gemv = (ck_gemv_t) dlsym(h, "dgemv_");
    C_.ok =C_.start && C_.finish && C_.asp && C_.fsp && C_.an && C_.fz && C_.cpf && C_.ff
            && C_.ad && C_.fd && C_.sv;
    return C_.ok;
}

/* THE THREADS: a row above a subtree and its update, kept for the replay; the plan */
typedef struct { int32_t r; double v; } ck_upd;
struct ck_par;

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
    /* THE KEPT FACTORS (`ck_keep`): the factors a numeric keeps, its own
       among them (1: none kept), the bytes its kept ones may take, and the
       entries of J's values a comparison reads first */
    long keep;
    double keep_bytes;
    long nsig;
    long sig[CK_NSIG];
    /* THE OWN SOLVE for this symbolic's numerics (`ck_own`; 0: cholmod_solve) */
    int own;
    /* THE THREADS' plan for its factors (`ck_par_for`, made at the first solve
       on threads, for one count of them), and the count a plan could not be
       made for (0: none) */
    struct ck_par *par;
    long par_failed;
    /* THE BUDGET (`ck_room`): the bytes under which a new slot leaves the
       process's resident memory, and the margin it keeps (0: no limit) */
    double keep_limit, keep_margin;
} ck_sym;

/* a kept factor: L's values and J's (the solve reads J's), and whether they
   are a factor of those values, and when it was last the numeric's */
typedef struct {
    double *x, *Ax;
    int valid;
    int64_t used;
} ck_slot;

typedef struct {
    double *Ax, *Kx, *x, *rhs;
    cholmod_common cm;
    cholmod_sparse *K;
    cholmod_factor *L;
    cholmod_dense *B;
    /* THE KEPT FACTORS: the slots (`nslot` made, each `slot_xsize` of L's
       values), whether the numeric's own values are a factor of its `Ax`,
       the clock, the refactors served by a slot and the factorizations made */
    ck_slot slot[CK_KEEP_MAX];
    long nslot;
    size_t slot_xsize;
    int valid;
    int64_t clock, hits, facts;
    /* THE OWN SOLVE: its work vectors (the free nodes; the most rows below a
       supernode's triangle), made at its first solve, and the solves made by
       it and by cholmod_solve */
    double *w, *e;
    int64_t own_solves, ck_solves;
    /* THE THREADS: the updates each subtree keeps for rows above it (at the
       plan's `soff`) and their counts, made for the plan `ub_for` (its
       `gen`); the own solves whose passes the pool took */
    ck_upd *ub;
    int64_t *ucnt;
    long ub_for;
    int64_t par_solves;
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

/* THE THREADS' plan for a symbolic's factor and a count of threads (`T`):
   the subtrees of its supernodal tree (each its first descendant .. its
   root in CHOLMOD's postorder), heaviest first, and each one's offset in a
   numeric's buffer of the updates to rows above it (`soff`, `nsub` + 1); the
   top (the supernodes in no subtree: their ancestors), descending; the
   forward pass's replay by position -- a subtree (>= 0) or a top supernode
   (-1 - it); the most rows below a triangle; the factor it was made on
   (`ns` supernodes, `xsize` values) and its generation */
typedef struct ck_par {
    long T, nsub, ntop, nitem, most, ns, gen;
    size_t xsize;
    int32_t *sfirst, *sroot, *top, *item;
    int64_t *soff;
} ck_par;

static void ck_par_free(ck_par *P)
{
    if (!P) return;
    free(P->sfirst); free(P->sroot); free(P->top); free(P->item); free(P->soff);
    free(P);
}

static int ck_free_symbolic(void **ps, void *common)
{
    ck_sym *s = (ck_sym *) *ps;
    (void) common;
    if (!s) return 1;
    if (s->Lsym) C_.ff(&s->Lsym, &s->cm);
    ck_par_free(s->par);
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

/* THE KEPT FACTORS: the slots freed (L's values by CHOLMOD's allocator, which
   made them -- or made the numeric's own, swapped into a slot) */
static void ck_drop(ck_num *n)
{
    long i;
    for (i = 0; i < n->nslot; i++) {
        if (n->slot[i].x) C_.fre(n->slot_xsize, sizeof(double), n->slot[i].x, &n->cm);
        free(n->slot[i].Ax);
        memset(&n->slot[i], 0, sizeof(ck_slot));
    }
    n->nslot = 0;
}

static int ck_free_numeric(void **pn, void *common)
{
    ck_num *n = (ck_num *) *pn;
    (void) common;
    if (!n) return 1;
    ck_drop(n);
    if (n->B) C_.fd(&n->B, &n->cm);
    if (n->L) C_.ff(&n->L, &n->cm);
    if (n->K) C_.fsp(&n->K, &n->cm);
    C_.finish(&n->cm);
    free(n->Ax); free(n->x); free(n->rhs);
    free(n->w); free(n->e);
    free(n->ub); free(n->ucnt);
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
        if (ck_values(s, n, Ax)) {
            n->valid = 1;
            n->facts = 1;
            return n;
        }
    }
    v = n;
    ck_free_numeric(&v, common);
    return NULL;
}

/* THE KEPT FACTORS.  Whether J's values `a` and `b` are the same, bit for bit
   (-0.0 is not 0.0) -- the signature's entries first */
static int ck_same(const ck_sym *s, const double *a, const double *b)
{
    long k;
    for (k = 0; k < s->nsig; k++)
        if (memcmp(a + s->sig[k], b + s->sig[k], sizeof(double)) != 0) return 0;
    return memcmp(a, b, (size_t) s->nnz * sizeof(double)) == 0;
}

/* the numeric's factor and a slot's exchanged: L's values, J's, and whether
   they are a factor of them -- pointers, no copy */
static void ck_swap(ck_num *n, ck_slot *t)
{
    double *p = (double *) n->L->x;
    int v = n->valid;
    n->L->x = t->x; t->x = p;
    p = n->Ax; n->Ax = t->Ax; t->Ax = p;
    n->valid = t->valid; t->valid = v;
}

/* the slots a numeric may make: `keep` - 1, CK_KEEP_MAX at most, within
   `keep_bytes` (a slot is L's values and J's); none where CHOLMOD's
   allocator was not found or L is not a supernodal numeric */
static long ck_cap(const ck_sym *s, const ck_num *n)
{
    long cap = s->keep - 1;
    double one;
    if (cap <= 0 || !C_.mal || !C_.fre || !n->L || !n->L->x || !n->L->is_super) return 0;
    if (cap > CK_KEEP_MAX) cap = CK_KEEP_MAX;
    one = 8.0 * ((double) n->L->xsize + (double) s->nnz);
    if (s->keep_bytes < one * (double) cap) cap = (long) (s->keep_bytes / one);
    return cap > 0 ? cap : 0;
}

/* THE BUDGET: whether a new slot (L's values and J's) leaves the process's
   resident memory and `keep_margin` under `keep_limit` -- 1 where no limit
   is set, 0 where the resident memory cannot be read */
static int ck_room(const ck_sym *s, const ck_num *n)
{
    long size = 0, pages = -1;
    FILE *f;
    if (s->keep_limit <= 0.0) return 1;
    f = fopen("/proc/self/statm", "r");
    if (!f) return 0;
    if (fscanf(f, "%ld %ld", &size, &pages) != 2) pages = -1;
    fclose(f);
    if (pages < 0) return 0;
    return (double) pages * (double) sysconf(_SC_PAGESIZE)
           + 8.0 * ((double) n->L->xsize + (double) s->nnz) + s->keep_margin <= s->keep_limit;
}

/* a refactor where factors are kept: a slot holding a factor of exactly
   `Ax` swapped in, no factorization (CHOLMOD's factor of the same values on
   the same symbolic: the same bytes); else the numeric's factor swapped into
   a new slot (room) or the least recently used one -- one holding none
   first -- and CHOLMOD's factor of `Ax` made in place on the buffers that
   came back.  1, or 0 as `ck_values` */
static int ck_kept(const ck_sym *s, ck_num *n, const double *Ax)
{
    long i, pick = -1, cap;
    int ok;
    /* (L's values another size: the slots do not fit it) */
    if (n->nslot && n->L->xsize != n->slot_xsize) ck_drop(n);
    cap = ck_cap(s, n);
    n->clock++;
    for (i = 0; i < n->nslot; i++)
        if (n->slot[i].valid && ck_same(s, n->slot[i].Ax, Ax)) {
            ck_swap(n, &n->slot[i]);
            n->slot[i].used = n->clock;
            n->hits++;
            return 1;
        }
    if (n->nslot < cap && ck_room(s, n)) {
        ck_slot *t = &n->slot[n->nslot];
        t->x = (double *) C_.mal(n->L->xsize, sizeof(double), &n->cm);
        t->Ax = (double *) ck_alloc((size_t) s->nnz, sizeof(double));
        if (t->x && t->Ax) {
            t->valid = 0;
            t->used = 0;
            n->slot_xsize = n->L->xsize;
            pick = n->nslot++;
        } else {
            if (t->x) C_.fre(n->L->xsize, sizeof(double), t->x, &n->cm);
            free(t->Ax);
            t->x = NULL;
            t->Ax = NULL;
        }
    }
    if (pick < 0)
        for (i = 0; i < n->nslot; i++)
            if (pick < 0 || (n->slot[pick].valid
                             && (!n->slot[i].valid || n->slot[i].used < n->slot[pick].used)))
                pick = i;
    if (pick >= 0) {
        ck_swap(n, &n->slot[pick]);
        n->slot[pick].used = n->clock;
    }
    ok = ck_values(s, n, Ax);
    n->facts++;
    n->valid = ok;
    return ok;
}

static int ck_refactor(int *Ap, int *Ai, double *Ax, void *sym, void *num, void *common)
{
    const ck_sym *s = (const ck_sym *) sym;
    ck_num *n = (ck_num *) num;
    int ok;
    (void) Ap; (void) Ai; (void) common;
    if (!sym || !num) return 0;
    if (s->keep > 1) return ck_kept(s, n, Ax);
    ok = ck_values(s, n, Ax);
    n->facts++;
    n->valid = ok;
    return ok;
}

/* THE KEPT FACTORS set on the symbolic: `keep` factors a numeric holds, its
   own among them (1: none kept), `bytes` the most its kept ones take, THE
   BUDGET's `limit` (0: none) and `margin`; the signature's entries spread
   over J's values */
static void ck_keep(void *sym, long keep, double bytes, double limit, double margin)
{
    ck_sym *s = (ck_sym *) sym;
    long k;
    if (!s) return;
    s->keep = keep < 1 ? 1 : keep;
    s->keep_bytes = bytes > 0.0 ? bytes : 0.0;
    s->keep_limit = limit > 0.0 ? limit : 0.0;
    s->keep_margin = margin > 0.0 ? margin : 0.0;
    s->nsig = s->nnz < CK_NSIG ? s->nnz : CK_NSIG;
    for (k = 0; k < s->nsig; k++)
        s->sig[k] = (long) (((double) k + 0.5) * (double) s->nnz / (double) s->nsig);
}

/* a numeric's kept factors: [the refactors a slot served, the factorizations
   made, the slots made, those holding a factor, L's values a slot], then its
   solves: [by THE OWN SOLVE, by cholmod_solve, the own ones the pool took
   (THE THREADS)] */
static void ck_num_info(void *num, double *out)
{
    const ck_num *n = (const ck_num *) num;
    long i, v = 0;
    if (!n) {
        out[0] = out[1] = out[2] = out[3] = out[4] = out[5] = out[6] = out[7] = -1.0;
        return;
    }
    for (i = 0; i < n->nslot; i++) v += n->slot[i].valid != 0;
    out[0] = (double) n->hits; out[1] = (double) n->facts; out[2] = (double) n->nslot;
    out[3] = (double) v; out[4] = (double) n->slot_xsize;
    out[5] = (double) n->own_solves; out[6] = (double) n->ck_solves;
    out[7] = (double) n->par_solves;
}

/* THE OWN SOLVE (the module's note): a supernode of CK_BLAS_FROM columns or
   more by the BLAS, a narrower one in plain loops */
#define CK_BLAS_FROM 32

/* whether the loops read L: supernodal, real doubles, int32, and the BLAS
   entries found */
static int ck_own_reads(const cholmod_factor *L)
{
    return L && L->is_super && L->is_ll && L->xtype == CHOLMOD_REAL && L->dtype == CHOLMOD_DOUBLE
           && L->itype == CHOLMOD_INT && L->x && L->Perm && C_.trsv && C_.gemv;
}

/* the work vectors, made at the numeric's first own solve: 1, or 0 (no memory) */
static int ck_own_ready(ck_num *n)
{
    const cholmod_factor *L = n->L;
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    long k, most = 0;
    if (n->w && n->e) return 1;
    for (k = 0; k < (long) L->nsuper; k++) {
        long below = (long) (pi[k + 1] - pi[k]) - (long) (super[k + 1] - super[k]);
        if (below > most) most = below;
    }
    n->w = (double *) ck_alloc(L->n, sizeof(double));
    n->e = (double *) ck_alloc((size_t) most, sizeof(double));
    if (n->w && n->e) return 1;
    free(n->w); free(n->e);
    n->w = n->e = NULL;
    return 0;
}

/* L y = w on one supernode of two columns or more (`Lx` its values,
   column-major with `nsrow` rows; its columns k1 ..): the triangle, then its
   rows below it summed into `e` -- how many, returned: the caller takes them
   out of w (`ck_out`) */
static long ck_own_down(const double *Lx, long k1, long nscol, long nsrow, double *w, double *e)
{
    long below = nsrow - nscol, i, j;
    if (nscol >= CK_BLAS_FROM) {
        const int nc = (int) nscol, nr = (int) nsrow, nb = (int) below, one = 1;
        const double d1 = 1.0, d0 = 0.0;
        C_.trsv("L", "N", "N", &nc, Lx, &nr, w + k1, &one, 1, 1, 1);
        if (nb <= 0) return 0;
        C_.gemv("N", &nb, &nc, &d1, Lx + nscol, &nr, w + k1, &one, &d0, e, &one, 1);
        return below;
    }
    for (j = 0; j < nscol; j++) {
        const double *c = Lx + j * nsrow;
        const double xj = w[k1 + j] / c[j];
        w[k1 + j] = xj;
        for (i = j + 1; i < nscol; i++) w[k1 + i] -= c[i] * xj;
    }
    if (below <= 0) return 0;
    for (i = 0; i < below; i++) e[i] = 0.0;
    for (j = 0; j + 4 <= nscol; j += 4) {
        const double *c0 = Lx + j * nsrow + nscol, *c1 = c0 + nsrow, *c2 = c1 + nsrow, *c3 = c2 + nsrow;
        const double x0 = w[k1 + j], x1 = w[k1 + j + 1], x2 = w[k1 + j + 2], x3 = w[k1 + j + 3];
        for (i = 0; i < below; i++) e[i] += (c0[i] * x0 + c1[i] * x1) + (c2[i] * x2 + c3[i] * x3);
    }
    for (; j < nscol; j++) {
        const double *c = Lx + j * nsrow + nscol;
        const double xj = w[k1 + j];
        for (i = 0; i < below; i++) e[i] += c[i] * xj;
    }
    return below;
}

/* the terms `e` of a supernode's rows below (`Lr`, `below` of them) out of
   w -- or, a row at or past `lim` (above a subtree: THE THREADS), the row
   and its term kept in `ub` in this order for the replay */
static inline void ck_out(const int32_t *Lr, long below, const double *e, double *w, long lim,
                          ck_upd *ub, int64_t *nu)
{
    long i;
    for (i = 0; i < below; i++) {
        const int32_t r = Lr[i];
        if (r >= lim) {
            ub[*nu].r = r;
            ub[*nu].v = e[i];
            (*nu)++;
        } else {
            w[r] -= e[i];
        }
    }
}

/* L y = w over the supernodes a .. b in order: a one-column one inline --
   its unknown over the diagonal, each row below's term the one a wider
   supernode's sum makes (0.0 + l x) -- the others by `ck_own_down`; the rows
   below out of w as `ck_out` takes them (`lim`, `ub`, `nu`) */
static void ck_down_range(const cholmod_factor *L, long a, long b, double *w, double *e, long lim,
                          ck_upd *ub, int64_t *nu)
{
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    const int32_t *px = (const int32_t *) L->px, *ls = (const int32_t *) L->s;
    const double *lx = (const double *) L->x;
    long sn, i;
    for (sn = a; sn <= b; sn++) {
        const long k1 = super[sn], nc = super[sn + 1] - k1, nr = pi[sn + 1] - pi[sn];
        const double *c = lx + px[sn];
        const int32_t *Ls = ls + pi[sn];
        if (nc == 1) {
            const double xj = w[k1] / c[0];
            w[k1] = xj;
            for (i = 1; i < nr; i++) {
                const double ei = 0.0 + c[i] * xj;
                const int32_t r = Ls[i];
                if (r >= lim) {
                    ub[*nu].r = r;
                    ub[*nu].v = ei;
                    (*nu)++;
                } else {
                    w[r] -= ei;
                }
            }
            continue;
        }
        ck_out(Ls + nc, ck_own_down(c, k1, nc, nr, w, e), e, w, lim, ub, nu);
    }
}

/* L' z = w on one supernode of two columns or more: its rows below the
   triangle gathered into `e` and their sums out of its columns -- each
   column's sum 0.0 then its rows in order, eight columns at a time, then
   four, then one -- then the triangle transposed */
static void ck_own_up(const double *Lx, const int32_t *Ls, long k1, long nscol, long nsrow,
                      double *w, double *e)
{
    long below = nsrow - nscol, i, j;
    for (i = 0; i < below; i++) e[i] = w[Ls[nscol + i]];
    if (nscol >= CK_BLAS_FROM) {
        const int nc = (int) nscol, nr = (int) nsrow, nb = (int) below, one = 1;
        const double d1 = 1.0, dm1 = -1.0;
        if (nb > 0) C_.gemv("T", &nb, &nc, &dm1, Lx + nscol, &nr, e, &one, &d1, w + k1, &one, 1);
        C_.trsv("L", "T", "N", &nc, Lx, &nr, w + k1, &one, 1, 1, 1);
        return;
    }
    for (j = 0; j + 8 <= nscol; j += 8) {
        const double *c0 = Lx + j * nsrow + nscol, *c1 = c0 + nsrow, *c2 = c1 + nsrow, *c3 = c2 + nsrow;
        const double *c4 = c3 + nsrow, *c5 = c4 + nsrow, *c6 = c5 + nsrow, *c7 = c6 + nsrow;
        double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0, t4 = 0.0, t5 = 0.0, t6 = 0.0, t7 = 0.0;
        for (i = 0; i < below; i++) {
            const double ei = e[i];
            t0 += c0[i] * ei; t1 += c1[i] * ei; t2 += c2[i] * ei; t3 += c3[i] * ei;
            t4 += c4[i] * ei; t5 += c5[i] * ei; t6 += c6[i] * ei; t7 += c7[i] * ei;
        }
        w[k1 + j] -= t0; w[k1 + j + 1] -= t1; w[k1 + j + 2] -= t2; w[k1 + j + 3] -= t3;
        w[k1 + j + 4] -= t4; w[k1 + j + 5] -= t5; w[k1 + j + 6] -= t6; w[k1 + j + 7] -= t7;
    }
    for (; j + 4 <= nscol; j += 4) {
        const double *c0 = Lx + j * nsrow + nscol, *c1 = c0 + nsrow, *c2 = c1 + nsrow, *c3 = c2 + nsrow;
        double t0 = 0.0, t1 = 0.0, t2 = 0.0, t3 = 0.0;
        for (i = 0; i < below; i++) {
            const double ei = e[i];
            t0 += c0[i] * ei; t1 += c1[i] * ei; t2 += c2[i] * ei; t3 += c3[i] * ei;
        }
        w[k1 + j] -= t0; w[k1 + j + 1] -= t1; w[k1 + j + 2] -= t2; w[k1 + j + 3] -= t3;
    }
    for (; j < nscol; j++) {
        const double *c = Lx + j * nsrow + nscol;
        double t = 0.0;
        for (i = 0; i < below; i++) t += c[i] * e[i];
        w[k1 + j] -= t;
    }
    for (j = nscol - 1; j >= 0; j--) {
        const double *c = Lx + j * nsrow;
        double t = w[k1 + j];
        for (i = j + 1; i < nscol; i++) t -= c[i] * w[k1 + i];
        w[k1 + j] = t / c[j];
    }
}

/* L' z = w over the supernodes b down to a: a one-column one inline (its
   rows' sum 0.0 then in order, out of its unknown, over the diagonal) */
static void ck_up_range(const cholmod_factor *L, long a, long b, double *w, double *e)
{
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    const int32_t *px = (const int32_t *) L->px, *ls = (const int32_t *) L->s;
    const double *lx = (const double *) L->x;
    long sn, i;
    for (sn = b; sn >= a; sn--) {
        const long k1 = super[sn], nc = super[sn + 1] - k1, nr = pi[sn + 1] - pi[sn];
        const double *c = lx + px[sn];
        const int32_t *Ls = ls + pi[sn];
        if (nc == 1) {
            double t = 0.0, u;
            for (i = 1; i < nr; i++) t += c[i] * w[Ls[i]];
            u = w[k1] - t;
            w[k1] = u / c[0];
            continue;
        }
        ck_own_up(c, Ls, k1, nc, nr, w, e);
    }
}

/* x = P' L'^-1 L^-1 P b on L (`ck_own_reads`), `w` and `e` the work vectors */
static void ck_own_solve(const cholmod_factor *L, const double *b, double *x, double *w, double *e)
{
    const int32_t *perm = (const int32_t *) L->Perm;
    const long n = (long) L->n, ns = (long) L->nsuper;
    int64_t nu = 0;
    long k;
    for (k = 0; k < n; k++) w[k] = b[perm[k]];
    ck_down_range(L, 0, ns - 1, w, e, LONG_MAX, NULL, &nu);
    ck_up_range(L, 0, ns - 1, w, e);
    for (k = 0; k < n; k++) x[perm[k]] = w[k];
}

/* THE THREADS (the module's note): pool.c's entry, the threads a solve may
   take, and the factor's values from which it takes them (`ck_threads`) */
typedef long (*ck_pool_t)(void *range, const void *ctx, long n, long nthreads);
static ck_pool_t ck_pool;
static long ck_nthreads = 1, ck_threads_from = 1L << 20, ck_par_gen;
/* the most rows below a triangle a plan takes (its passes' `e` on the stack) */
#define CK_PAR_ESIZE 8192

static void ck_threads(void *pool, long n, long from)
{
    ck_pool = (ck_pool_t) pool;
    ck_nthreads = n < 1 ? 1 : n;
    ck_threads_from = from < 0 ? 0 : from;
}

/* the supernodal tree: each supernode's parent (the one holding its first
   row below the triangle; -1 a root), its subtree's first supernode and its
   work (its values and a supernode's cost, summed over the subtree); 0
   where a parent does not follow its child (not CHOLMOD's postorder) or a
   subtree is not the range first .. it */
static int ck_tree(const cholmod_factor *L, int32_t *parent, int32_t *first, double *work)
{
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    const int32_t *ls = (const int32_t *) L->s;
    const long ns = (long) L->nsuper, n = (long) L->n;
    int32_t *col = (int32_t *) malloc((size_t) (n ? n : 1) * sizeof(int32_t)), *cnt;
    long sn, k;
    int ok = 1;
    if (!col) return 0;
    for (sn = 0; sn < ns; sn++) for (k = super[sn]; k < super[sn + 1]; k++) col[k] = (int32_t) sn;
    for (sn = 0; sn < ns; sn++) {
        const long nc = super[sn + 1] - super[sn], nr = pi[sn + 1] - pi[sn];
        parent[sn] = nr > nc ? col[ls[pi[sn] + nc]] : -1;
        if (parent[sn] >= 0 && parent[sn] <= sn) ok = 0;
        first[sn] = (int32_t) sn;
        work[sn] = (double) nc * (double) nr + 40.0;
    }
    free(col);
    cnt = (int32_t *) calloc((size_t) (ns ? ns : 1), sizeof(int32_t));
    if (!ok || !cnt) { free(cnt); return 0; }
    for (sn = 0; sn < ns; sn++) {
        cnt[sn]++;
        if (parent[sn] < 0) continue;
        work[parent[sn]] += work[sn];
        cnt[parent[sn]] += cnt[sn];
        if (first[sn] < first[parent[sn]]) first[parent[sn]] = first[sn];
    }
    for (sn = 0; sn < ns && ok; sn++) ok = cnt[sn] == sn - first[sn] + 1;
    free(cnt);
    return ok;
}

/* a max-heap of subtree roots by their work */
static void ck_heap_push(int32_t *h, long *nh, const double *work, int32_t r)
{
    long k = (*nh)++;
    while (k > 0 && work[h[(k - 1) / 2]] < work[r]) { h[k] = h[(k - 1) / 2]; k = (k - 1) / 2; }
    h[k] = r;
}

static int32_t ck_heap_pop(int32_t *h, long *nh, const double *work)
{
    const int32_t top = h[0], last = h[--(*nh)];
    long k = 0, c;
    while ((c = 2 * k + 1) < *nh) {
        if (c + 1 < *nh && work[h[c + 1]] > work[h[c]]) c++;
        if (work[h[c]] <= work[last]) break;
        h[k] = h[c];
        k = c;
    }
    if (*nh > 0) h[k] = last;
    return top;
}

/* the subtrees: from the roots, the heaviest split into its children (its
   root to the top) while fewer than 32 a thread, or one weighs more than
   twice its share of that many, and it has children; their roots in `h`
   (`nh`), the top's in `top` (`ntop`) */
static void ck_split(long ns, const int32_t *parent, const double *work, long T, int32_t *h, long *nh,
                     int32_t *top, long *ntop)
{
    int32_t *kid0 = (int32_t *) calloc((size_t) ns + 2, sizeof(int32_t));
    int32_t *kids = (int32_t *) malloc((size_t) (ns ? ns : 1) * sizeof(int32_t));
    const long want = 32 * T;
    double total = 0.0;
    long sn, c;
    *nh = *ntop = 0;
    if (!kid0 || !kids) { free(kid0); free(kids); return; }
    for (sn = 0; sn < ns; sn++) if (parent[sn] >= 0) kid0[parent[sn] + 2]++;
    for (sn = 0; sn < ns; sn++) kid0[sn + 2] += kid0[sn + 1];
    for (sn = 0; sn < ns; sn++) if (parent[sn] >= 0) kids[kid0[parent[sn] + 1]++] = (int32_t) sn;
    for (sn = 0; sn < ns; sn++)
        if (parent[sn] < 0) { ck_heap_push(h, nh, work, (int32_t) sn); total += work[sn]; }
    while (*nh > 0) {
        const int32_t r = h[0];
        if ((*nh >= want && work[r] <= 2.0 * total / (double) want) || kid0[r] == kid0[r + 1]) break;
        ck_heap_pop(h, nh, work);
        top[(*ntop)++] = r;
        for (c = kid0[r]; c < kid0[r + 1]; c++) ck_heap_push(h, nh, work, kids[c]);
    }
    free(kid0); free(kids);
}

/* sort keys: a subtree by its work (heaviest first, then its root), a
   position (ascending) */
typedef struct { double key; int32_t a, b; } ck_key;

static int ck_by_work(const void *x, const void *y)
{
    const ck_key *p = (const ck_key *) x, *q = (const ck_key *) y;
    if (p->key != q->key) return p->key > q->key ? -1 : 1;
    return (p->a > q->a) - (p->a < q->a);
}

static int ck_by_key_up(const void *x, const void *y)
{
    const ck_key *p = (const ck_key *) x, *q = (const ck_key *) y;
    return (p->key > q->key) - (p->key < q->key);
}

/* the plan's subtrees heaviest first, each one's buffer offset (the rows
   above it its supernodes update), the top descending and the replay by
   position, from the split's roots `h` and top; 1, or 0 (no memory) */
static int ck_plan_fill(const cholmod_factor *L, ck_par *P, const int32_t *first, const double *work,
                        const int32_t *h, const int32_t *top)
{
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    const int32_t *ls = (const int32_t *) L->s;
    ck_key *k = (ck_key *) malloc((size_t) (P->nsub + P->ntop + 1) * sizeof(ck_key));
    long q, sn, i;
    if (!k) return 0;
    for (q = 0; q < P->nsub; q++) { k[q].key = work[h[q]]; k[q].a = h[q]; k[q].b = 0; }
    qsort(k, (size_t) P->nsub, sizeof(ck_key), ck_by_work);
    P->soff[0] = 0;
    for (q = 0; q < P->nsub; q++) {
        const long lim = super[k[q].a + 1];
        int64_t c = 0;
        P->sroot[q] = k[q].a;
        P->sfirst[q] = first[k[q].a];
        for (sn = P->sfirst[q]; sn <= P->sroot[q]; sn++)
            for (i = pi[sn] + (super[sn + 1] - super[sn]); i < pi[sn + 1]; i++) c += ls[i] >= lim;
        P->soff[q + 1] = P->soff[q] + c;
    }
    for (q = 0; q < P->ntop; q++) { k[q].key = -(double) top[q]; k[q].a = top[q]; }
    qsort(k, (size_t) P->ntop, sizeof(ck_key), ck_by_key_up);
    for (q = 0; q < P->ntop; q++) P->top[q] = k[q].a;
    for (q = 0; q < P->nsub; q++) { k[q].key = (double) P->sfirst[q]; k[q].a = (int32_t) q; }
    for (q = 0; q < P->ntop; q++) { k[P->nsub + q].key = (double) P->top[q]; k[P->nsub + q].a = -1 - P->top[q]; }
    qsort(k, (size_t) P->nitem, sizeof(ck_key), ck_by_key_up);
    for (q = 0; q < P->nitem; q++) P->item[q] = k[q].a;
    free(k);
    return 1;
}

/* THE THREADS' plan for L on `T` threads (`ck_par`), or NULL: no memory, not
   CHOLMOD's postorder, a triangle with more rows below it than
   CK_PAR_ESIZE */
static ck_par *ck_plan(const cholmod_factor *L, long T)
{
    const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
    const long ns = (long) L->nsuper;
    int32_t *parent = (int32_t *) malloc((size_t) (ns ? ns : 1) * sizeof(int32_t));
    int32_t *first = (int32_t *) malloc((size_t) (ns ? ns : 1) * sizeof(int32_t));
    int32_t *h = (int32_t *) malloc((size_t) (ns ? ns : 1) * sizeof(int32_t));
    int32_t *top = (int32_t *) malloc((size_t) (ns ? ns : 1) * sizeof(int32_t));
    double *work = (double *) malloc((size_t) (ns ? ns : 1) * sizeof(double));
    ck_par *P = (ck_par *) calloc(1, sizeof(ck_par));
    long k, nh = 0, ntop = 0;
    int ok = parent && first && h && top && work && P && ck_tree(L, parent, first, work);
    if (ok) {
        for (k = 0; k < ns; k++) {
            const long below = (pi[k + 1] - pi[k]) - (super[k + 1] - super[k]);
            if (below > P->most) P->most = below;
        }
        ck_split(ns, parent, work, T, h, &nh, top, &ntop);
        P->T = T; P->nsub = nh; P->ntop = ntop; P->nitem = nh + ntop; P->ns = ns; P->xsize = L->xsize;
        P->gen = ++ck_par_gen;
        P->sfirst = (int32_t *) malloc((size_t) (nh + 1) * sizeof(int32_t));
        P->sroot = (int32_t *) malloc((size_t) (nh + 1) * sizeof(int32_t));
        P->soff = (int64_t *) malloc((size_t) (nh + 1) * sizeof(int64_t));
        P->top = (int32_t *) malloc((size_t) (ntop + 1) * sizeof(int32_t));
        P->item = (int32_t *) malloc((size_t) (nh + ntop + 1) * sizeof(int32_t));
        ok = P->most <= CK_PAR_ESIZE && nh > 0 && P->sfirst && P->sroot && P->soff && P->top && P->item
             && ck_plan_fill(L, P, first, work, h, top);
    }
    free(parent); free(first); free(h); free(top); free(work);
    if (ok) return P;
    ck_par_free(P);
    return NULL;
}

/* the plan a numeric's solve takes, its buffers made -- or NULL, the solve on
   the caller alone: fewer than two threads, a factor of fewer than
   `ck_threads_from` values, no plan for this count (tried once), fewer than
   two subtrees, no memory */
static const ck_par *ck_par_for(ck_sym *s, ck_num *nm)
{
    const cholmod_factor *L = nm->L;
    const ck_par *P;
    if (ck_nthreads < 2 || (double) L->xsize < (double) ck_threads_from) return NULL;
    if (s->par && (s->par->T != ck_nthreads || s->par->ns != (long) L->nsuper || s->par->xsize != L->xsize)) {
        ck_par_free(s->par);
        s->par = NULL;
    }
    if (!s->par && s->par_failed != ck_nthreads) {
        s->par = ck_plan(L, ck_nthreads);
        if (!s->par) s->par_failed = ck_nthreads;
    }
    P = s->par;
    if (!P || P->nsub < 2) return NULL;
    if (nm->ub_for != P->gen) {
        free(nm->ub); free(nm->ucnt);
        nm->ub = (ck_upd *) malloc((size_t) (P->soff[P->nsub] + 1) * sizeof(ck_upd));
        nm->ucnt = (int64_t *) malloc((size_t) P->nsub * sizeof(int64_t));
        nm->ub_for = nm->ub && nm->ucnt ? P->gen : 0;
        if (!nm->ub_for) return NULL;
    }
    return P;
}

/* a solve on the pool: the factor, the plan, the work vector, the buffers */
typedef struct {
    const cholmod_factor *L;
    const ck_par *P;
    double *w;
    ck_upd *ub;
    int64_t *ucnt;
} ck_job;

/* the forward pass on a job's subtrees lo .. hi - 1 (pool.c's range): each in
   order, its updates to the rows above it kept in its buffer */
static void ck_fwd_items(const void *ctx, long lo, long hi)
{
    const ck_job *j = (const ck_job *) ctx;
    const int32_t *super = (const int32_t *) j->L->super;
    double e[j->P->most + 1];
    long k;
    for (k = lo; k < hi; k++) {
        int64_t nu = 0;
        ck_down_range(j->L, j->P->sfirst[k], j->P->sroot[k], j->w, e, super[j->P->sroot[k] + 1],
                      j->ub + j->P->soff[k], &nu);
        j->ucnt[k] = nu;
    }
}

/* the backward pass on a job's subtrees lo .. hi - 1: each from its root down */
static void ck_bwd_items(const void *ctx, long lo, long hi)
{
    const ck_job *j = (const ck_job *) ctx;
    double e[j->P->most + 1];
    long k;
    for (k = lo; k < hi; k++) ck_up_range(j->L, j->P->sfirst[k], j->P->sroot[k], j->w, e);
}

/* x = P' L'^-1 L^-1 P b with the subtrees on the pool, the bytes of
   `ck_own_solve`: the forward pass's subtrees first, then in position order
   each subtree's kept updates and each top supernode (every row takes its
   updates in the serial order); the backward pass's top descending, then
   its subtrees (each reads only its ancestors, writes only its columns).  1
   where the pool took a pass, 0 where the caller made both alone */
static int ck_par_solve(const cholmod_factor *L, const ck_par *P, const double *b, double *x,
                        double *w, double *e, ck_upd *ub, int64_t *ucnt)
{
    const int32_t *perm = (const int32_t *) L->Perm;
    const long n = (long) L->n;
    ck_job j;
    long k, q;
    int took = 0;
    int64_t nu = 0;
    j.L = L; j.P = P; j.w = w; j.ub = ub; j.ucnt = ucnt;
    for (k = 0; k < n; k++) w[k] = b[perm[k]];
    if (ck_pool && ck_pool((void *) ck_fwd_items, &j, P->nsub, P->T)) took = 1;
    else ck_fwd_items(&j, 0, P->nsub);
    for (q = 0; q < P->nitem; q++) {
        const long it = P->item[q];
        if (it < 0) {
            ck_down_range(L, -1 - it, -1 - it, w, e, LONG_MAX, NULL, &nu);
        } else {
            const ck_upd *u = ub + P->soff[it];
            for (k = 0; k < ucnt[it]; k++) w[u[k].r] -= u[k].v;
        }
    }
    for (q = 0; q < P->ntop; q++) ck_up_range(L, P->top[q], P->top[q], w, e);
    if (ck_pool && ck_pool((void *) ck_bwd_items, &j, P->nsub, P->T)) took = 1;
    else ck_bwd_items(&j, 0, P->nsub);
    for (k = 0; k < n; k++) x[perm[k]] = w[k];
    return took;
}

/* J x = b for each of the `nrhs` columns of `B` (n x nrhs, column-major), in place */
static int ck_solve(void *sym, void *num, int n_, int nrhs, double *B, void *common)
{
    ck_sym *s = (ck_sym *) sym;
    ck_num *nm = (ck_num *) num;
    const double *Ax;
    double *x, *rhs, *sol;
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
        X = NULL;
        if (s->own && ck_own_reads(nm->L)) {
            if (!ck_own_ready(nm)) return 0;
            const ck_par *P = ck_par_for(s, nm);
            sol = (double *) nm->B->x;
            if (P) nm->par_solves += ck_par_solve(nm->L, P, rhs, sol, nm->w, nm->e, nm->ub, nm->ucnt);
            else ck_own_solve(nm->L, rhs, sol, nm->w, nm->e);
            nm->own_solves++;
        } else {
            memcpy(nm->B->x, rhs, (size_t) s->nf * sizeof(double));
            X = C_.sv(CHOLMOD_A, nm->L, nm->B, &nm->cm);
            if (!X) return 0;
            sol = (double *) X->x;
            nm->ck_solves++;
        }
        for (r = 0; r < s->m; r++)
            if (s->kind[r] == 0 && s->fidx[r] >= 0) x[r] = sol[s->fidx[r]];
        if (X) C_.fd(&X, &nm->cm);
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

/* the factor's statistics: [the free nodes, K's entries, L's entries, the
   ordering CHOLMOD chose, its supernodes, those of CK_BLAS_FROM columns or
   more, of those the ones with rows below their triangle] (-1: none known) */
static void ck_info(void *sym, double *out)
{
    const ck_sym *s = (const ck_sym *) sym;
    const cholmod_factor *L = s->Lsym;
    out[0] = (double) s->nf; out[1] = (double) s->nk;
    out[2] = L ? (double) s->cm.lnz : -1.0;
    out[3] = (double) s->cm.selected;
    out[4] = out[5] = out[6] = -1.0;
    if (L && L->is_super && L->itype == CHOLMOD_INT && L->super && L->pi) {
        const int32_t *super = (const int32_t *) L->super, *pi = (const int32_t *) L->pi;
        long k, wide = 0, below = 0;
        for (k = 0; k < (long) L->nsuper; k++) {
            const long nscol = super[k + 1] - super[k];
            if (nscol < CK_BLAS_FROM) continue;
            wide++;
            below += pi[k + 1] - pi[k] > nscol;
        }
        out[4] = (double) L->nsuper; out[5] = (double) wide; out[6] = (double) below;
    }
}

/* THE OWN SOLVE (1) or cholmod_solve (0) for the symbolic's numerics */
static void ck_own(void *sym, int on)
{
    ck_sym *s = (ck_sym *) sym;
    if (s) s->own = on != 0;
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
    case 9: return (void *) ck_keep;
    case 10: return (void *) ck_num_info;
    case 11: return (void *) ck_own;
    case 12: return (void *) ck_threads;
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
    `analyze(m, Ap, Ai, kind, rg_off)` the symbolic and `info(sym, out)` its
    statistics (seven doubles: `ck_info`); `num_info(num, out)` a numeric's
    counts (seven: `ck_num_info`)."""

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
        ## THE KEPT FACTORS: `keep(sym, keep, bytes, limit, margin)` (THE
        ## BUDGET's limit and margin) and a numeric's counts `num_info(num,
        ## out)` (the module's note)
        self.keep = fn(9, None, vp, ct.c_long, ct.c_double, ct.c_double, ct.c_double)
        self.num_info = fn(10, None, vp, dp)
        ## THE OWN SOLVE: `own(sym, on)` (the module's note)
        self.own = fn(11, None, vp, ct.c_int)
        ## THE THREADS: `threads(pool, n, from)` -- pool.c's entry, the count,
        ## the factor's values from which a solve takes them
        self.threads = fn(12, None, vp, ct.c_long, ct.c_long)
        self.klu_defaults = klu.klu_defaults


def memory_limit(proc='/proc', cgroup='/sys/fs/cgroup'):
    """THE BUDGET's limit (the module's note), in bytes: the least of the
    `memory.max` up this process's cgroup (v2) tree and its resident memory
    plus the machine's MemAvailable; None where neither can be read."""
    got = []
    try:
        with open(os.path.join(proc, 'self', 'statm')) as f:
            rss = int(f.read().split()[1]) * os.sysconf('SC_PAGESIZE')
        with open(os.path.join(proc, 'meminfo')) as f:
            avail = next(int(ln.split()[1]) * 1024 for ln in f if ln.startswith('MemAvailable:'))
        got.append(rss + avail)
    except (OSError, ValueError, IndexError, StopIteration, AttributeError):
        pass
    try:
        with open(os.path.join(proc, 'self', 'cgroup')) as f:
            parts = next(ln[3:].strip().strip('/') for ln in f if ln.startswith('0::')).split('/')
    except (OSError, StopIteration):
        parts = None
    if parts is not None:
        for k in range(len(parts), -1, -1):
            try:
                with open(os.path.join(cgroup, *[q for q in parts[:k] if q], 'memory.max')) as f:
                    v = f.read().strip()
                if v != 'max':
                    got.append(int(v))
            except (OSError, ValueError):
                pass
    return float(min(got)) if got else None


class _Threads:
    """THE THREADS' count, as `_tran_core.set_threads` writes every live
    struct's (`nthreads`): handed to the C with pool.c's entry and
    `THREADS_FROM` -- one thread where `THREADS` is off or the pool cannot
    be had."""

    def __init__(self, lib):
        self._lib = lib
        self._n = 1

    @property
    def nthreads(self):
        return self._n

    @nthreads.setter
    def nthreads(self, n):
        from pycircuit.circuit import _tran_core
        self._n = int(n)
        pool = _tran_core.pool_address() if THREADS and self._n > 1 else 0
        self._lib.threads(pool or None, self._n if pool else 1, int(THREADS_FROM))


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
    lib = _Lib(ffi, entry, klu)
    ## THE THREADS: the count `set_threads` keeps (held on the library)
    from pycircuit.circuit import _tran_core
    lib.threads_count = _tran_core._threaded(_Threads(lib))
    _LIB.append(lib)
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
