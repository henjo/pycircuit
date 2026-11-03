/* pool.h -- the thread pool's entry: the interface `_tran_core.pool_address`
 * loads pool.c through; pool.c includes it, so the compiler checks the entry
 * against it.  Its callers reach it through a `void *` of their structs (the
 * core's `pool_fn`, the walk's context), cast to pool_fn_t (pcx_abi.h). */

/* THE POOL'S JOB (pool.c): the items [0, n) of a range function's context */
typedef struct {
    void *range;          /* pool_range_t (pcx_abi.h): the items [lo, hi) */
    const void *ctx;      /* its context */
    long n;               /* the items */
    long chunk;           /* the items a thread takes at a time (pcx_pool) */
} pool_job_t;

long pcx_pool(void *range, const void *ctx, long n, long nthreads);
