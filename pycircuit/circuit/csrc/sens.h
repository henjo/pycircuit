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

/* hdl_fn's result, left in `status` too (`_sens_c.BAIL` names the bails
   from these): */
enum sens_status {
    SENS_DONE = 1,
    SENS_BAIL_NONFINITE = 2,     /* a reduced entry not finite */
    SENS_BAIL_FLAGS = 3,         /* a floating-point flag numpy would raise on */
    SENS_BAIL_RHS = 4,           /* C_n P not finite */
    SENS_BAIL_SINGULAR = 5,      /* the block's factorisation: singular */
    SENS_BAIL_LAPACK = 6         /* the block's solve failed */
};
