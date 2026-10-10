/* lte.h -- the C error test's state and entry: the interface `_tran_lte_c` loads
 * lte.c through.  cffi lays the struct out from this text and
 * lte.c includes it, so the compiler checks the entry against it. */
typedef struct {
    long n, iref, nn, relref, have_run, abs_vec, has_mask;
    double h1, h2, h3, s12, s23, s123, coef, reltol, abstol, trtol;
    const double *xc, *xl, *q, *q1, *q2, *q3, *J, *run, *abst, *mask;
    double *A, *b, *lte, *ref, *out, *tq;
    int64_t *ipiv;
    void *dgesv;
    double err;
    long status;
} lte_t;
long hdl_fn(lte_t *s);
