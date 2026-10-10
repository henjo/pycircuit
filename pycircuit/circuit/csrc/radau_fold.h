/* radau_fold.h -- the transformed Radau fold's state and entry: the interface `_tran_radau_tc` loads
 * radau_fold.c through.  cffi lays the struct out from this text and
 * radau_fold.c includes it, so the compiler checks the entry against it. */
typedef struct {
    long m, fused, nnz, nz_at;
    double s0, c1r, c1i;
    const double *Cr, *Gr;
    const uint8_t *nz;
    double *rf, *cf, *ax;
    long status, same;
} radau_fold_t;
long hdl_fn(radau_fold_t *f);
