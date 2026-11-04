
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "pcx_abi.h"
#include "core.h"
/* passes 0 G, 1 C, 2 i, 3 q; their want bits pcx_pass (pcx_abi.h).  A class
   with a fused group (z_fn) is evaluated once a call, before the passes,
   for every pass the call wants, into its staging (zst, in vals); each
   pass then pulls its bins from the values its slots read (`one_pass`) --
   the bins as the slots' bincount makes them.  b_fz, a batch's group (-1:
   none) */
/* a ztab row: core.h's `core_ztab` names its columns */
/* `body` for each matrix entry the call writes (`patt`: the pattern's; else
   all n x n), `s` the entry.  It captures from the function: `patt`, `pk`
   (the pattern's cursor), `s` (assigned: the entry's row-major index -- an
   entry index here, not a struct), `c` (the core: npat, pat) and `n`.  The
   entries off the pattern are left as the caller's buffers hold them
   (zero) */
#define EACH_ENTRY(body) do { \
        if (patt) { for (pk = 0; pk < c->npat; pk++) { s = c->pat[pk]; body; } } \
        else { for (s = 0; s < n * n; s++) { body; } } \
    } while (0)
/* J packed (PCX_WANT_PACKED): for each pattern entry, `J[pk] = binG[s] +
   (geq)`, `geq` the companion's entry at `s` as the formula makes it -- the
   same operations as `Geq[s] = geq` then `J[s] = binG[s] + Geq[s]` -- and
   no Geq.  It captures `pk`, `s`, `c` and `J` */
#define PACKED_J(geq) do { \
        for (pk = 0; pk < c->npat; pk++) { s = c->pat[pk]; J[pk] = c->binG[s] + (geq); } \
    } while (0)
/* one element of fused group z into its staging rows, for the passes `zw`:
   its inputs gathered into X, the temperature into its parameter row (as
   the batches' driver, pass.c, writes it), its kernel -- a function of
   those alone, writing its own rows (any thread may make it) */
static void fused_element(const hdl_core_t *c, long z, long zw, long e, const double *x,
                          double T, double *X)
{
    const int64_t *zt = c->ztab + ZT_WIDTH * z;
    long k = zt[ZT_K], ti = zt[ZT_TI], j;
    const int64_t *NM = c->NMall + zt[ZT_NM];
    double **PR = c->PRall + zt[ZT_PR], *st = c->zst;
    for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
    if (ti >= 0) PR[e][ti] = T;
    ((fused_fn_t) c->z_fn[z])(X, PR[e], zw, st + zt[ZT_G_AT] + e * zt[ZT_G_LEN],
                              st + zt[ZT_C_AT] + e * zt[ZT_C_LEN],
                              st + zt[ZT_I_AT] + e * zt[ZT_I_LEN],
                              st + zt[ZT_Q_AT] + e * zt[ZT_Q_LEN]);
}

/* the passes this call makes that group z covers (q always); 0 where its
   fused kernel does not serve or the call makes none of them */
static long group_want(const hdl_core_t *c, long z, long want)
{
    if (!c->z_fn[z]) return 0;
    return ((want & PCX_WANT_ALL) | PCX_PASS_Q) & c->ztab[ZT_WIDTH * z + ZT_PASSES];
}

/* THE THREADS: a call's fused-group elements on `nthreads` threads (the
   pool, pool.c's `pcx_pool`, through `pool_fn`) where they reach POOL_MIN.
   Each element is `fused_element` on the same inputs whatever the thread,
   into its own rows; the bins are summed after by the caller, in their
   order: every byte as one thread makes it, the flags the threads raise
   raised in the caller's (pool.c). */
/* (POOL_MIN measured on the gold decks, four threads against one: ring's
   34 elements a call 1.65x, nand's 25 1.31x, gm3's 30 1.14x; at 8, the
   decks of 11-14 elements 1.00-1.15x -- below 16 no thread is started) */
enum { POOL_MIN = 16 };
_Static_assert(sizeof(((core_pool_job_t *) 0)->zl) == CORE_POOL_ZMAX * sizeof(long),
               "the pool job's rows are CORE_POOL_ZMAX wide");

/* the pool job's elements [lo, hi) (pool_range_t): each its group's kernel
   (a chunk's elements come in rising order: the group cursor only moves on) */
static void fused_range(const void *jv, long lo, long hi)
{
    const core_pool_job_t *j = (const core_pool_job_t *) jv;
    double X[CORE_MAX_INPUTS];
    long r = 0, i;
    for (i = lo; i < hi; i++) {
        while (i >= j->off[r + 1]) r++;
        fused_element(j->c, j->zl[r], j->zw[r], i - j->off[r], j->x, j->T, X);
    }
}

/* the fused groups on the pool: 1 where it made them, 0 where the caller
   makes them alone (too few elements, too many groups, no pool, the pool
   another caller's or not running) */
static int pool_groups(const hdl_core_t *c, const double *x, double T, long want)
{
    core_pool_job_t j;
    long z, zw, total = 0;
    if (!c->pool_fn || c->nz > CORE_POOL_ZMAX) return 0;
    for (z = 0; z < c->nz; z++)
        if (group_want(c, z, want)) total += c->ztab[ZT_WIDTH * z + ZT_NE];
    if (total < POOL_MIN) return 0;
    j.c = c;
    j.x = x;
    j.T = T;
    for (z = 0, j.nr = 0, total = 0; z < c->nz; z++) {
        if (!(zw = group_want(c, z, want))) continue;
        j.zl[j.nr] = z;
        j.zw[j.nr] = zw;
        j.off[j.nr++] = total;
        total += c->ztab[ZT_WIDTH * z + ZT_NE];
    }
    j.off[j.nr] = j.total = total;
    return (int) ((pool_fn_t) c->pool_fn)((void *) fused_range, &j, total, c->nthreads);
}

/* the fused groups this call wants into their staging (`zst`): each
   group's kernel once an element (`fused_element`), for every pass it
   covers that the call makes (q always) -- on the pool where it serves */
static void fused_groups(const hdl_core_t *c, const double *x, double T, long want, double *X)
{
    long z, e, zw;
    if (c->nthreads > 1 && pool_groups(c, x, T, want)) return;
    for (z = 0; z < c->nz; z++) {
        if (!(zw = group_want(c, z, want))) continue;
        for (e = 0; e < c->ztab[ZT_WIDTH * z + ZT_NE]; e++) fused_element(c, z, zw, e, x, T, X);
    }
}

/* batch b of pass m where it makes its outputs itself -- outside a group,
   or its group's kernel not serving the call (`fused_groups` made the rest):
   its class's kernel per element for pass m alone (its inputs gathered, the
   temperature into its parameter row), the outputs into the element's rows
   of vals (`orowall`), where the pass's slots read them */
static void batch_into(const hdl_core_t *c, long b, long m, const double *x, double T,
                       double *X)
{
    long e, j;
    long ne = c->b_ne[b], k = c->b_k[b], ti = c->b_ti[b], fz = c->b_fz[b];
    fused_fn_t f = (fused_fn_t) c->b_fn[b];
    const int64_t *NM = c->NMall + c->b_nmoff[b], *orow = c->orowall + c->b_ooff[b];
    double **PR = c->PRall + c->b_proff[b], *o;
    if (fz >= 0 && c->z_fn[fz]) return;
    for (e = 0; e < ne; e++) {
        for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
        if (ti >= 0) PR[e][ti] = T;
        o = c->vals + orow[e];
        f(X, PR[e], 1L << m, o, o, o, o);
    }
}

/* THE SMALL PRODUCTS: a k x k row-major product, k 2 or 3, as OpenBLAS's
   dgemv makes it (the row's middle term, the first fused onto it, then the
   third; added into the zeroed y) -- taken only where `small_gemv` says this
   process's dgemv was checked to make exactly that (`_tran_core`); the
   multiplication alone is never contracted (-ffp-contract=off) and `0.0 +`
   keeps a -0.0's sign as the zeroed y's addition does */
static void small_gemv(const double *S, const double *x, long k, double *y)
{
    long r;
    for (r = 0; r < k; r++) {
        const double *a = S + r * k;
        double t = a[1] * x[1];
        t = fma(a[0], x[0], t);
        if (k == 3) t = fma(a[2], x[2], t);
        y[r] = 0.0 + t;
    }
}

/* pass m's constant groups: each element's product into the group's
   products in vals (`g_yoff`), where the slots they fill read them */
static void groups_into(const hdl_core_t *c, long m, const double *x, double *xg)
{
    dgemv64_t gemv = (dgemv64_t) c->dgemv;
    long g, e, j;
    for (g = 0; g < c->ng; g++) {
        /* a constant group: `np.dot(S, x[I])` per element through
           numpy's own dgemv, the routine np.dot calls (the plan's
           batched np.matmul is taken only where `_matmul_ok` found it
           equal to that; else the core is Unservable) -- or its
           operation order where checked (`small_gemv`) */
        if (c->g_pass[g] != m) continue;
        long E = c->g_E[g], k = c->g_k[g];
        const double *S = c->Sall + c->g_soff[g];
        const int64_t *I = c->Iall + c->g_ioff[g];
        double *Y = c->vals + c->g_yoff[g];
        for (e = 0; e < E; e++) {
            for (j = 0; j < k; j++) xg[j] = x[I[e * k + j]];
            if ((k == 2 || k == 3) && (c->small_gemv & (1L << k)))
                small_gemv(S + e * k * k, xg, k, Y + e * k);
            else
                gemv(PCX_CBLAS_ROW_MAJOR, PCX_CBLAS_NO_TRANS, k, k, 1.0, S + e * k * k, k, xg, 1, 0.0,
                     Y + e * k, 1);
        }
    }
}

/* a run of `nb` bins of `cnt` sources each (core.h's THE PULL): each bin
   from +0.0, its sources' values added in order, one store.  The counts up
   to four written out, a longer one a loop of `cnt` turns -- the same for
   every bin of the run: no branch that a bin's count decides (a loop
   whose length changed bin by bin mispredicted its exit at most bins, and
   left the predictor that much worse for the code after it) */
static void pull_run(const double *V, const int64_t *ix, const int64_t *wb, double *bin, long cnt,
                     long nb)
{
    long q, k;
    switch (cnt) {
    case 0:
        for (q = 0; q < nb; q++) bin[wb[q]] = 0.0;
        break;
    case 1:
        for (q = 0; q < nb; q++) bin[wb[q]] = 0.0 + V[ix[q]];
        break;
    case 2:
        for (q = 0; q < nb; q++, ix += 2) bin[wb[q]] = (0.0 + V[ix[0]]) + V[ix[1]];
        break;
    case 3:
        for (q = 0; q < nb; q++, ix += 3)
            bin[wb[q]] = ((0.0 + V[ix[0]]) + V[ix[1]]) + V[ix[2]];
        break;
    case 4:
        for (q = 0; q < nb; q++, ix += 4)
            bin[wb[q]] = (((0.0 + V[ix[0]]) + V[ix[1]]) + V[ix[2]]) + V[ix[3]];
        break;
    default:
        for (q = 0; q < nb; q++, ix += cnt) {
            double t = 0.0;
            for (k = 0; k < cnt; k++) t += V[ix[k]];
            bin[wb[q]] = t;
        }
    }
}

/* pass m (0 G, 1 C, 2 i, 3 q) into its bins `bin`: the values its slots
   read made -- the batches that make their own (`batch_into`), the constant
   groups (`groups_into`); the templates are set once, the groups' staging by
   `fused_groups` -- then THE PULL, run by run (`pull_run`): each written bin
   from +0.0, its slots' values added in slot order, one store --
   `np.bincount(flat, weights=slots, minlength=nbins)`'s sums in its order:
   a bin started at its first value would keep a lone -0.0 (bincount's is
   +0.0), and another order rounds differently.  A matrix pass whose
   caller's buffer is not zero off the pattern (`patt` clear where there is
   one) has every entry zeroed first */
static void one_pass(const hdl_core_t *c, long m, const double *x, double T, double *bin,
                     int patt, double *X, double *xg)
{
    const int64_t *run = c->runall + 2 * c->pr_off[m];
    const int64_t *wb = c->wball + c->pl_off[m], *ix = c->ixall + c->px_off[m];
    long n = c->n, b, r, nr = c->pr_off[m + 1] - c->pr_off[m];
    for (b = 0; b < c->nb; b++)
        if (c->b_pass[b] == m) batch_into(c, b, m, x, T, X);
    groups_into(c, m, x, xg);
    if (m < 2 && c->npat && !patt) memset(bin, 0, (size_t) (n * n) * sizeof(double));
    for (r = 0; r < nr; r++) {
        long cnt = run[2 * r], nb = run[2 * r + 1];
        pull_run(c->vals, ix, wb, bin, cnt, nb);
        ix += cnt * nb;
        wb += nb;
    }
}

/* the entry (contract: core.h) */
long pcx_core(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J)
{
    /* sizes `_Core.__init__` holds every element to (it refuses, as
       Unservable, core.h's `core_cap`): an element's inputs, a group's
       x[I] */
    double X[CORE_MAX_INPUTS], xg[CORE_MAX_INPUTS];
    /* each pass's bins: G's and i's the core's own, C's and q's the
       caller's */
    double *bins[4] = {c->binG, C, c->bini, q};
    /* each pass's want bit; q's 0: always made */
    const long bits[4] = {PCX_PASS_G, PCX_PASS_C, PCX_PASS_I, 0};
    /* (`s`: an entry of a matrix -- an index) */
    long n = c->n, m, j, s, pk;
    /* the matrices' entries this call writes: the pattern's where the caller
       asks (its matrices zero off the pattern, which no stamp reaches), or
       every one */
    const int patt = (want & PCX_WANT_PATTERN) && c->npat;
    /* J packed on the pattern, no Geq (PCX_WANT_PACKED): only with it */
    const int packed = patt && (want & PCX_WANT_PACKED);
    const double *Cs;
    fused_groups(c, x, T, want, X);
    for (m = 0; m < 4; m++)
        if (!bits[m] || (want & bits[m])) one_pass(c, m, x, T, bins[m], patt, X, xg);
    /* formula 4: the passes alone (`passes`: G in binG, i in bini, C and q
       in the caller's buffers) -- no companion, no residual, no Jacobian */
    if (formula == PCX_FORMULA_PASSES) return 0;
    /* the companion's C: the one just made, else the caller's held one
       (`_C_lookup`'s hit, the chord's held C, `residual_only`'s `_Cmat`) */
    Cs = (want & PCX_PASS_C) ? C : Cin;
    /* the companions as `_lte_kernels` writes them, in Python's
       left-to-right order: the parentheses below are the ones Python's
       evaluation puts there, and another grouping rounds differently */
    if (formula == PCX_FORMULA_GEAR2) {
        /* `a0 * q + a1 * q1 + a2 * q2`, `a0 * C` */
        for (j = 0; j < n; j++) iq[j] = (a0 * q[j] + a1 * q1[j]) + a2 * q2[j];
        if (!packed) EACH_ENTRY(Geq[s] = a0 * Cs[s]);
        else if (want & PCX_PASS_G) PACKED_J(a0 * Cs[s]);
    } else if (formula == PCX_FORMULA_EULER) {
        /* `(q - q1) / h`, `C / h` */
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / h;
        if (!packed) EACH_ENTRY(Geq[s] = Cs[s] / h);
        else if (want & PCX_PASS_G) PACKED_J(Cs[s] / h);
    } else if (formula == PCX_FORMULA_TRAPEZOIDAL) {
        /* `2.0 * (q - q1) / h - iq1`, `2.0 * C / h` */
        for (j = 0; j < n; j++) iq[j] = (2.0 * (q[j] - q1[j])) / h - iq1[j];
        if (!packed) EACH_ENTRY(Geq[s] = (2.0 * Cs[s]) / h);
        else if (want & PCX_PASS_G) PACKED_J((2.0 * Cs[s]) / h);
    } else {  /* PCX_FORMULA_THETA */
        /* `th = theta * h`, `(q - q1) / th - ((1.0 - theta) / theta) * iq1`,
           `C / th`: the two scalars once, as Python computes them before
           the array operations */
        double th = theta * h, w = (1.0 - theta) / theta;
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / th - w * iq1[j];
        if (!packed) EACH_ENTRY(Geq[s] = Cs[s] / th);
        else if (want & PCX_PASS_G) PACKED_J(Cs[s] / th);
    }
    /* `f = (i + iq) + u`, `J = G + Geq` (packed: made above) */
    if (want & PCX_PASS_I) for (j = 0; j < n; j++) F[j] = (c->bini[j] + iq[j]) + u[j];
    if ((want & PCX_PASS_G) && !packed) EACH_ENTRY(J[s] = c->binG[s] + Geq[s]);
    return 0;
}
