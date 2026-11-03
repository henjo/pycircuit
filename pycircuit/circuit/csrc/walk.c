
#include <stdint.h>
#include <string.h>
#include "pcx_abi.h"
#include "walk.h"
/* THE WALK: each element's limiter (lim_fn_t) on its rows of the state,
   in the elements' order, its limited values written back before the next
   element reads -- an element sharing a node with an earlier one reads
   what that one wrote.

   AHEAD ON THE POOL, where the walk's context wires one and the run is
   long enough: every element's limiter made first on the state as the walk
   found it (`walk_range`, any thread: a limiter is a function of its
   inputs, its parameter row and its node rows alone), then the walk in its
   order (`walk_ahead`): an element none of whose nodes an earlier element
   changed -- its written value's bits not the value it found -- takes its
   result as made, which IS what it makes here, from the same inputs; any
   other element is made here, as the walk alone makes it.  Every byte as
   the walk alone leaves it.  The flags the limiters raise are the
   caller's: every caller clears them untested after the walk (newton.c
   tests none), and so the union the pool raises -- more than the walk alone
   raises where an element is made twice or past a stop -- reaches no test. */

/* the walk alone, from `start` */
static long walk_alone(void **F, double **PR, const int64_t *TI, double T,
                       const int64_t *NM, const int64_t *OFF, const int64_t *K,
                       double *xg, const double *x0g, long start, long n)
{
    double o[WALK_MAX_NODES], b[WALK_MAX_NODES];
    long e, j, k;
    for (e = start; e < n; e++) {
        const int64_t *nm;
        if (F[e] == 0) return e;
        k = K[e];
        nm = NM + OFF[e];
        for (j = 0; j < k; j++) { o[j] = xg[nm[j]]; b[j] = x0g[nm[j]]; }
        if (TI[e] >= 0) PR[e][TI[e]] = T;
        if (((lim_fn_t) F[e])(b, PR[e], o)) return e;
        for (j = 0; j < k; j++) xg[nm[j]] = o[j];
    }
    return n;
}

/* the job's elements [lo, hi) from its first (pool_range_t): each limiter
   on the state as the walk found it, into its buffer and its return */
static void walk_range(const void *jv, long lo, long hi)
{
    const walk_job_t *j = (const walk_job_t *) jv;
    double b[WALK_MAX_NODES];
    long e, q, k;
    for (e = j->start + lo; e < j->start + hi; e++) {
        const int64_t *nm;
        double *o;
        if (j->F[e] == 0) continue;
        k = j->K[e];
        nm = j->NM + j->OFF[e];
        o = j->so + j->OFF[e];
        for (q = 0; q < k; q++) { o[q] = j->xg[nm[q]]; b[q] = j->x0g[nm[q]]; }
        if (j->TI[e] >= 0) j->PR[e][j->TI[e]] = j->T;
        j->rc[e] = ((lim_fn_t) j->F[e])(b, j->PR[e], o);
    }
}

static int same_bits(double a, double b)
{
    uint64_t u, v;
    memcpy(&u, &a, sizeof u);
    memcpy(&v, &b, sizeof v);
    return u == v;
}

/* the walk in its order over the results made ahead: an element on no
   node an earlier one changed takes its own; the rest are made here (and
   past WALK_DIRTY_MAX changed nodes, every later element) */
static long walk_ahead(const walk_job_t *j, double *xg, long n)
{
    int64_t dirty[WALK_DIRTY_MAX];
    double o[WALK_MAX_NODES], b[WALK_MAX_NODES];
    long nd = 0, e, q, d, k;
    int all = 0;
    for (e = j->start; e < n; e++) {
        const int64_t *nm;
        const double *r;
        int64_t rc;
        int clean = !all;
        if (j->F[e] == 0) return e;
        k = j->K[e];
        nm = j->NM + j->OFF[e];
        for (q = 0; clean && q < k; q++)
            for (d = 0; d < nd; d++)
                if (dirty[d] == nm[q]) { clean = 0; break; }
        if (clean) {
            r = j->so + j->OFF[e];
            rc = j->rc[e];
        } else {
            for (q = 0; q < k; q++) { o[q] = xg[nm[q]]; b[q] = j->x0g[nm[q]]; }
            if (j->TI[e] >= 0) j->PR[e][j->TI[e]] = j->T;
            rc = ((lim_fn_t) j->F[e])(b, j->PR[e], o);
            r = o;
        }
        if (rc) return e;
        for (q = 0; q < k; q++) {
            if (!same_bits(xg[nm[q]], r[q])) {
                if (nd < WALK_DIRTY_MAX) dirty[nd++] = nm[q];
                else all = 1;
            }
            xg[nm[q]] = r[q];
        }
    }
    return n;
}

long pcx_walk(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n, void *ctx)
{
    const walk_ctx_t *cx = (const walk_ctx_t *) ctx;
    if (cx && cx->pool_fn && cx->nthreads > 1 && n - start >= WALK_POOL_MIN && n <= cx->n) {
        walk_job_t j;
        j.F = F;
        j.PR = PR;
        j.TI = TI;
        j.T = T;
        j.NM = NM;
        j.OFF = OFF;
        j.K = K;
        j.xg = xg;
        j.x0g = x0g;
        j.so = cx->so;
        j.rc = cx->rc;
        j.start = start;
        if (((pool_fn_t) cx->pool_fn)((void *) walk_range, &j, n - start, cx->nthreads))
            return walk_ahead(&j, xg, n);
    }
    return walk_alone(F, PR, TI, T, NM, OFF, K, xg, x0g, start, n);
}
