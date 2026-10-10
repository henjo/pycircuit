/* lte.h -- the C error test's state and entry: the interface `_tran_lte_c` loads
 * lte.c through.  cffi lays the struct out from this text and
 * lte.c includes it, so the compiler checks the entry against it. */
typedef struct {
    long n, iref, nn, relref, have_run, abs_vec, has_mask, formula;
    double h1, h2, h3, s12, s23, s123, coef, reltol, abstol, trtol;
    const double *xc, *xl, *q, *q1, *q2, *q3, *iq1, *J, *run, *abst, *mask;
    double *A, *b, *lte, *ref, *out, *tq;
    int64_t *ipiv;
    void *dgesv;
    double err;
    long status;
} lte_t;
long hdl_fn(lte_t *s);

/* hdl_fn's result, left in `status` too (`_tran_lte_c.BAIL` names the
   bails from these): */
enum lte_status {
    LTE_DONE = 1,
    LTE_BAIL_NONFINITE = 2,      /* an input not finite */
    LTE_BAIL_FLAGS = 3,          /* a floating-point flag numpy would raise on */
    LTE_BAIL_DGESV = 4,          /* the divided difference's solve: singular */
    LTE_BAIL_LTE = 5,            /* the estimate not finite */
    LTE_BAIL_TOLERANCE = 6       /* a tolerance not positive and finite */
};
/* `formula`, the estimate (the integrator's `compute_lte`): `coef` the
   glue's, in the integrator's own scalar arithmetic */
enum lte_formula {
    LTE_FORMULA_DD3 = 0,         /* gear, trapezoid: coef * the third divided
                                    difference of the charge (q, q1..q3) */
    LTE_FORMULA_EULER = 1        /* Euler: coef * ((q - q1) / h1 - iq1) */
};
/* `relref`, the reference the tolerance scales (Transient's `relref`) */
enum lte_relref {
    LTE_RELREF_POINTLOCAL = 0,   /* this point's own |x| */
    LTE_RELREF_ALLLOCAL = 1,     /* each unknown's running maximum */
    LTE_RELREF_SIGGLOBAL = 2     /* each unit group's running maximum (any other) */
};
