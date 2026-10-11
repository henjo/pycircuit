/* walk.h -- the limiter walk's entry: the interface `_hdl_climit` loads
 * walk.c through; walk.c includes it, so the compiler checks the
 * entry against it. */
long pcx_walk(void **F, double **PR, const int64_t *TI, double T, const int64_t *NM, const int64_t *OFF, const int64_t *K, double *xg, const double *x0g, long start, long n);

/* an element's node map, for the walk's stack buffers (walk.c) --
   `_hdl_climit` walks no element with more nodes */
enum walk_cap {
    WALK_MAX_NODES = 64
};
