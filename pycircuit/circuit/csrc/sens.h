/* sens.h -- the shooting sensitivities' state and entry: the interface `_sens_c` loads
 * sens.c through.  cffi lays the struct out from this text and
 * sens.c includes it, so the compiler checks the entry against it. */
typedef struct {
    void *dgetrf, *dgetrs;
    long s, m, k, nf, iref;
    double h;
    const double *A, *Cf, *Gf, *base;
    double *J; int *ipiv; double *Z;
    long status, info;
} sens_t;
long hdl_fn(sens_t *t);
