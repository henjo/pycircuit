/* pass.h -- the batched pass driver's entry: the interface `_hdl_batch` loads
 * pass.c through; pass.c includes it, so the compiler checks the
 * entry against it. */
void hdl_fn(const void *fn, const double *X, double **PR, double *OUT, long n, long sx, long so, long ti, double T);
