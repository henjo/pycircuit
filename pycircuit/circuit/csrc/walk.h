/* walk.h -- the limiter walk's entry: the interface `_hdl_climit` loads
 * walk.c through; walk.c includes it, so the compiler checks the
 * entry against it. */

/* THE WALK'S CONTEXT (`_hdl_climit._Walk`): the pool its elements' limiters
   are made on ahead of the walk, and their results' buffers */
typedef struct {
    void *pool_fn;        /* pool.c's `pcx_pool` (pool_fn_t, pcx_abi.h), or
                             NULL: the walk alone */
    long nthreads;        /* the threads, the caller one of them (1: alone) */
    long n;               /* the elements the buffers hold */
    double *so;           /* each element's result, at its node map's offset
                             (the node maps' total) */
    int64_t *rc;          /* (n) each element's limiter's return */
} walk_ctx_t;

/* ONE WALK'S JOB ON THE POOL (walk.c): its tables, the state as the walk
   found it, the buffers, its first element */
typedef struct {
    void **F;             /* (n) each element's limiter (lim_fn_t; 0: none) */
    double **PR;          /* (n) its parameter row */
    const int64_t *TI;    /* (n) the temperature's slot in it (-1: none) */
    double T;             /* the temperature */
    const int64_t *NM,    /* the node maps, concatenated */
                  *OFF,   /* (n) each one's offset in NM */
                  *K;     /* (n) its node count */
    const double *xg,     /* the state as the walk found it */
                 *x0g;    /* the point the limiters limit against */
    double *so;           /* the context's results' buffer */
    int64_t *rc;          /* the context's returns' buffer */
    long start;           /* the walk's first element */
} walk_job_t;

long pcx_walk(void **F, double **PR, const int64_t *TI, double T, const int64_t *NM, const int64_t *OFF, const int64_t *K, double *xg, const double *x0g, long start, long n, void *ctx);

/* an element's node map, for the walk's stack buffers (walk.c) --
   `_hdl_climit` walks no element with more nodes; the walk's elements
   below which it runs alone; the nodes a walk tracks as written before it
   takes every later element alone */
enum walk_cap {
    WALK_MAX_NODES = 64,
    WALK_POOL_MIN = 16,
    WALK_DIRTY_MAX = 64
};
