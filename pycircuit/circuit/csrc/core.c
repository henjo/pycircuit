
#define _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#include <dlfcn.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include "pcx_abi.h"
#include "core.h"
/* passes 0 G, 1 C, 2 i, 3 q; their want bits pcx_pass (pcx_abi.h).  A class
   with a fused kernel (z_fn) is evaluated once a call, before the passes,
   for every pass the call wants, into its staging (zst); its batches then
   scatter from there -- the passes' writes and bins as before.  bftab, a
   batch's row: its group (-1: none) and its offset in fidxall (each
   element's place in it; -1: the group's own order, the usual case) */
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

/* THE POOL: a call's fused-group elements on `nthreads` threads -- the
   caller and workers kept from call to call (spinning a while for the
   next call, then asleep), each taking chunks of the elements from one
   counter -- where they reach POOL_MIN.  Each element is `fused_element`
   on the same inputs whatever the thread, into its own rows; the bins are
   summed after by the caller, in their order: every byte as one thread
   makes it.  The workers take the caller's FP environment each call and
   the flags they raise are raised in the caller's after (the union, as
   one thread leaves them).  One call at a time on the pool (a concurrent
   caller makes its own alone); a forked child starts workers of its own;
   the library pins itself (RTLD_NODELETE) before a worker runs in it. */
/* (POOL_MIN measured on the gold decks, four threads against one: ring's
   34 elements a call 1.65x, nand's 25 1.31x, gm3's 30 1.14x; at 8, the
   decks of 11-14 elements 1.00-1.15x -- below 16 no thread is started) */
enum { POOL_MAX = 31, POOL_CHUNK = 4, POOL_MIN = 16, POOL_SPIN = 1 << 16 };
_Static_assert(sizeof(((core_pool_job_t *) 0)->zl) == CORE_POOL_ZMAX * sizeof(long),
               "the pool job's rows are CORE_POOL_ZMAX wide");

static pthread_mutex_t pool_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t pool_wake = PTHREAD_COND_INITIALIZER;
static long pool_started;
static pid_t pool_pid;
static atomic_long pool_gen, pool_next, pool_left;
static atomic_int pool_flags;
static atomic_flag pool_busy = ATOMIC_FLAG_INIT;
static core_pool_job_t pool_job;
static fenv_t pool_env;

static void pool_relax(void)
{
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#endif
}

/* the job's elements, a chunk at a time from the shared counter, until
   none is left (a thread's chunks come in rising order: its group cursor
   only moves on) */
static void pool_take(const core_pool_job_t *j)
{
    double X[CORE_MAX_INPUTS];
    long r = 0, i, end;
    for (;;) {
        i = atomic_fetch_add_explicit(&pool_next, POOL_CHUNK, memory_order_relaxed);
        if (i >= j->total) return;
        end = i + POOL_CHUNK < j->total ? i + POOL_CHUNK : j->total;
        for (; i < end; i++) {
            while (i >= j->off[r + 1]) r++;
            fused_element(j->c, j->zl[r], j->zw[r], i - j->off[r], j->x, j->T, X);
        }
    }
}

/* a worker's wait for the call after `seen`: spinning, then asleep */
static long pool_wait(long seen)
{
    long spin;
    for (spin = 0; spin < POOL_SPIN; spin++) {
        if (atomic_load_explicit(&pool_gen, memory_order_acquire) != seen) goto out;
        pool_relax();
    }
    pthread_mutex_lock(&pool_mu);
    while (atomic_load_explicit(&pool_gen, memory_order_acquire) == seen)
        pthread_cond_wait(&pool_wake, &pool_mu);
    pthread_mutex_unlock(&pool_mu);
out:
    return atomic_load_explicit(&pool_gen, memory_order_acquire);
}

/* a worker: each call (the caller waits for every worker before the
   next, so none is missed), in the caller's FP environment, its flags
   handed back */
static void *pool_worker(void *arg)
{
    long seen = (long) (intptr_t) arg;
    int fl;
    for (;;) {
        seen = pool_wait(seen);
        fesetenv(&pool_env);
        feclearexcept(FE_ALL_EXCEPT);
        pool_take(&pool_job);
        fl = fetestexcept(FE_ALL_EXCEPT);
        if (fl) atomic_fetch_or_explicit(&pool_flags, fl, memory_order_relaxed);
        atomic_fetch_sub_explicit(&pool_left, 1, memory_order_release);
    }
    return NULL;
}

/* this library kept loaded while workers may run in it */
static int pool_pin(void)
{
    Dl_info info;
    if (!dladdr((void *) pool_worker, &info) || !info.dli_fname) return 0;
    return dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE) != NULL;
}

/* `n - 1` workers (at most POOL_MAX, named "pcx-core") running in this
   process: 1 where any is.  In a forked child the parent's are gone: its
   own, from none */
static int pool_ready(long n)
{
    long want = n - 1 < POOL_MAX ? n - 1 : POOL_MAX;
    pthread_attr_t at;
    pthread_t t;
    if (pool_pid != getpid()) {
        pool_pid = getpid();
        pool_started = 0;
        pthread_mutex_init(&pool_mu, NULL);
        pthread_cond_init(&pool_wake, NULL);
        if (!pool_pin()) pool_pid = 0;
    }
    if (!pool_pid) return 0;
    while (pool_started < want) {
        pthread_attr_init(&at);
        pthread_attr_setdetachstate(&at, PTHREAD_CREATE_DETACHED);
        if (pthread_create(&t, &at, pool_worker,
                           (void *) (intptr_t) atomic_load(&pool_gen)) != 0) {
            pthread_attr_destroy(&at);
            break;
        }
        pthread_attr_destroy(&at);
        pthread_setname_np(t, "pcx-core");
        pool_started++;
    }
    return pool_started > 0;
}

/* one call on the pool: the job published, the caller's share taken, every
   worker's waited for, their flags raised here */
static void pool_run(const core_pool_job_t *j)
{
    int fl;
    fegetenv(&pool_env);
    atomic_store_explicit(&pool_next, 0, memory_order_relaxed);
    atomic_store_explicit(&pool_flags, 0, memory_order_relaxed);
    atomic_store_explicit(&pool_left, pool_started, memory_order_relaxed);
    pthread_mutex_lock(&pool_mu);
    atomic_fetch_add_explicit(&pool_gen, 1, memory_order_release);
    pthread_cond_broadcast(&pool_wake);
    pthread_mutex_unlock(&pool_mu);
    pool_take(j);
    while (atomic_load_explicit(&pool_left, memory_order_acquire) > 0) pool_relax();
    fl = atomic_load_explicit(&pool_flags, memory_order_relaxed);
    if (fl) feraiseexcept(fl);
}

/* the fused groups on the pool: 1 where it made them, 0 where the caller
   makes them alone (too few elements, too many groups, the pool another
   caller's or not running) */
static int pool_groups(const hdl_core_t *c, const double *x, double T, long want)
{
    core_pool_job_t *j = &pool_job;
    long z, zw, total = 0;
    if (c->nz > CORE_POOL_ZMAX) return 0;
    for (z = 0; z < c->nz; z++)
        if (group_want(c, z, want)) total += c->ztab[ZT_WIDTH * z + ZT_NE];
    if (total < POOL_MIN || atomic_flag_test_and_set_explicit(&pool_busy, memory_order_acquire))
        return 0;
    if (!pool_ready(c->nthreads)) {
        atomic_flag_clear_explicit(&pool_busy, memory_order_release);
        return 0;
    }
    j->c = c;
    j->x = x;
    j->T = T;
    for (z = 0, j->nr = 0, total = 0; z < c->nz; z++) {
        if (!(zw = group_want(c, z, want))) continue;
        j->zl[j->nr] = z;
        j->zw[j->nr] = zw;
        j->off[j->nr++] = total;
        total += c->ztab[ZT_WIDTH * z + ZT_NE];
    }
    j->off[j->nr] = j->total = total;
    pool_run(j);
    atomic_flag_clear_explicit(&pool_busy, memory_order_release);
    return 1;
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

/* batch b of pass m into its slots `buf`: from the staging where its
   class's fused kernel served (`fused_groups`), else its kernel per
   element -- its inputs gathered, the temperature into its parameter row
   -- each output scattered to its destination */
static void batch_into(const hdl_core_t *c, long b, long m, double *buf, const double *x,
                       double T, double *X, double *o)
{
    long e, j;
    pass_fn_t f = (pass_fn_t) c->b_fn[b];
    long ne = c->b_ne[b], k = c->b_k[b], so = c->b_so[b], ti = c->b_ti[b];
    long fz = c->nz ? c->bftab[2 * b] : -1;
    const int64_t *NM = c->NMall + c->b_nmoff[b];
    const int64_t *dst = c->bdstall + c->b_dstoff[b];
    double **PR = c->PRall + c->b_proff[b];
    if (fz >= 0 && c->z_fn[fz]) {
        /* from the staging: the fused kernel wrote this pass's
           rows (so wide: `_fuse_groups` checked the widths) */
        const double *sz = c->zst + c->ztab[ZT_WIDTH * fz + ZT_G_AT + 2 * m];
        long fo = c->bftab[2 * b + 1];
        const int64_t *fi = fo >= 0 ? c->fidxall + fo : 0;
        for (e = 0; e < ne; e++) {
            const double *oz = sz + (fi ? fi[e] : e) * so;
            for (j = 0; j < so; j++) buf[dst[e * so + j]] = oz[j];
        }
        return;
    }
    for (e = 0; e < ne; e++) {
        for (j = 0; j < k; j++) X[j] = x[NM[e * k + j]];
        if (ti >= 0) PR[e][ti] = T;
        f(X, PR[e], o);
        for (j = 0; j < so; j++) buf[dst[e * so + j]] = o[j];
    }
}

/* pass m's constant groups into its slots `buf` */
static void groups_into(const hdl_core_t *c, long m, double *buf, const double *x,
                        double *y, double *xg)
{
    dgemv64_t gemv = (dgemv64_t) c->dgemv;
    long g, e, j, s;
    for (g = 0; g < c->ng; g++) {
        /* a constant group: `np.dot(S, x[I])` per element through
           numpy's own dgemv, the routine np.dot calls (the plan's
           batched np.matmul is taken only where `_matmul_ok` found it
           equal to that; else the core is Unservable), then
           `buf[dst] = Y.reshape(-1)[sel]` */
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
}

/* pass m (0 G, 1 C, 2 i, 3 q) into its bins `bin`: the slots -- G's and
   C's from the template (the constant stamps' entries), i's and q's
   from zeros -- each batch's outputs (`batch_into`), the constant
   groups' (`groups_into`), then the bins */
static void one_pass(const hdl_core_t *c, long m, const double *x, double T, double *bin,
                     int patt, double *X, double *o, double *y, double *xg)
{
    double *buf = c->work;
    long L = c->p_L[m], n = c->n, b, s, pk;
    if (m < 2) memcpy(buf, c->tmplall + c->p_toff[m], (size_t) L * sizeof(double));
    else memset(buf, 0, (size_t) L * sizeof(double));
    for (b = 0; b < c->nb; b++)
        if (c->b_pass[b] == m) batch_into(c, b, m, buf, x, T, X, o);
    groups_into(c, m, buf, x, y, xg);
    /* `np.bincount(flat, weights=buf, minlength=nbins)`: every bin
       from +0.0, the slots added in their order -- a bin started
       at its first value would keep a lone -0.0 (bincount's is
       +0.0), and another order rounds differently */
    const int64_t *flat = c->flatall + c->p_toff[m];
    long nbins = c->p_nbins[m];
    if (patt && nbins == n * n)
        for (pk = 0; pk < c->npat; pk++) bin[c->pat[pk]] = 0.0;
    else
        memset(bin, 0, (size_t) nbins * sizeof(double));
    for (s = 0; s < L; s++) bin[flat[s]] += buf[s];
}

/* the entry (contract: core.h) */
long pcx_core(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J)
{
    /* sizes `_Core.__init__` holds every element to (it refuses, as
       Unservable, core.h's `core_cap`): an element's inputs, a kernel's
       outputs, a group's x[I] and its product */
    double X[CORE_MAX_INPUTS], o[CORE_MAX_OUTPUTS], y[CORE_MAX_INPUTS], xg[CORE_MAX_INPUTS];
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
    const double *Cs;
    fused_groups(c, x, T, want, X);
    for (m = 0; m < 4; m++)
        if (!bits[m] || (want & bits[m])) one_pass(c, m, x, T, bins[m], patt, X, o, y, xg);
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
        EACH_ENTRY(Geq[s] = a0 * Cs[s]);
    } else if (formula == PCX_FORMULA_EULER) {
        /* `(q - q1) / h`, `C / h` */
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / h;
        EACH_ENTRY(Geq[s] = Cs[s] / h);
    } else if (formula == PCX_FORMULA_TRAPEZOIDAL) {
        /* `2.0 * (q - q1) / h - iq1`, `2.0 * C / h` */
        for (j = 0; j < n; j++) iq[j] = (2.0 * (q[j] - q1[j])) / h - iq1[j];
        EACH_ENTRY(Geq[s] = (2.0 * Cs[s]) / h);
    } else {  /* PCX_FORMULA_THETA */
        /* `th = theta * h`, `(q - q1) / th - ((1.0 - theta) / theta) * iq1`,
           `C / th`: the two scalars once, as Python computes them before
           the array operations */
        double th = theta * h, w = (1.0 - theta) / theta;
        for (j = 0; j < n; j++) iq[j] = (q[j] - q1[j]) / th - w * iq1[j];
        EACH_ENTRY(Geq[s] = Cs[s] / th);
    }
    /* `f = (i + iq) + u`, `J = G + Geq` */
    if (want & PCX_PASS_I) for (j = 0; j < n; j++) F[j] = (c->bini[j] + iq[j]) + u[j];
    if (want & PCX_PASS_G) EACH_ENTRY(J[s] = c->binG[s] + Geq[s]);
    return 0;
}
