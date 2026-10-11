/* lte.h -- the C error test's state and entry: the interface `_tran_lte_c` loads
 * lte.c through.  cffi lays the struct out from this text and
 * lte.c includes it, so the compiler checks the entry against it.
 *
 * WHO SETS WHAT.  `_tran_lte_c._Ctx` makes one struct per controller, size
 * and reference row: it sets `n` and `iref`, points every pointer below at a
 * buffer of its own (the sizes are given per field; m = n - 1 is the reduced
 * size) and sets `dgesv`; `_tran_lte_c._klu_bind` sets the sparse block.
 * Then one of two callers fills the rest:
 *   - `_tran_lte_c.max_error` (one Python error test): every [run] and
 *     [step] field at each call, copying the inputs into the context's
 *     buffers (`J` it may instead point at the caller's array where it
 *     stands);
 *   - the C loop: `_tran_loop_c._bind` sets the [run] fields once a run, on
 *     a context of the run's own, and loop.c sets the [step] fields at each
 *     attempt, re-pointing the [step] pointers at its own and the Newton's
 *     buffers.
 * hdl_fn reads the inputs and writes the work buffers, `out`, `err` and
 * `status` -- lte.c says what holds after a bail. */
typedef struct {
    long n,                     /* [ctx] unknowns: J is n x n, each vector n */
        iref,                   /* [ctx] the reference row, 0 <= iref < n: cut
                                   from J and Eg for the solve, restored as +0.0 */
        nn,                     /* [run] `n_nodes`: sigglobal's first unit group
                                   is the first nn unknowns (the node voltages),
                                   one group where <= 0 or >= n (the glue's -1:
                                   None) */
        relref,                 /* [run] the reference: enum lte_relref */
        have_run,               /* [step] 1: `run` holds the controller's running
                                   maximum, 0: none yet (the reference starts
                                   from this point's own) */
        abs_vec,                /* [run] 1: the absolute tolerance per unknown,
                                   `abst`, 0: the scalar `abstol` */
        has_mask,               /* [run] 1: only the unknowns `mask` marks are
                                   judged (`Transient.lte_algebraic`) */
        formula;                /* [step] the estimate: enum lte_formula */
    double h1, h2, h3,          /* [step] h_curr, h_last, h_last2 (Euler reads h1
                                   alone, the callers pass 0.0 for the others) */
        s12, s23, s123,         /* [step] the divided difference's step sums
                                   h1 + h2, h2 + h3, s12 + h3, made by the caller
                                   in that order (0.0 for Euler) */
        coef,                   /* [step] the estimate's coefficient, made by the
                                   caller in the integrator's scalar arithmetic */
        reltol, abstol, trtol;  /* [run] the tolerance TRTOL (reltol ref + abstol):
                                   trtol is TRTOL; `abstol` is 0.0 where abs_vec */
    const double *xc, *xl,      /* [step] x_curr, x_last (n each) */
        *q, *q1, *q2, *q3,      /* [step] q_curr and the three past charges, most
                                   recent first (n each) */
        *iq1,                   /* [step] the last companion current iq_last[0]
                                   (n; read for Euler only) */
        *J,                     /* [step] the Jacobian, n x n row-major; only read */
        *run,                   /* [step] the running maximum (n), read where
                                   have_run: finite and not negative, or handed
                                   back */
        *abst,                  /* [run] the absolute tolerance per unknown (n),
                                   read where abs_vec */
        *mask;                  /* [run] n: an unknown is judged where its entry
                                   is nonzero (a NaN counts); read where has_mask */
    /* work buffers and the one output vector: the context's own (`_Ctx`),
       pointed at once; their contents are this call's */
    double *A,                  /* m x m: the reduced J column-major, LAPACK's
                                   layout (the dense path) */
        *b,                     /* m: the reduced Eg, solved in place into the
                                   reduced J^-1 Eg */
        *lte,                   /* n: Eg, then the full J^-1 Eg (+0.0 at iref) */
        *ref,                   /* n: the tolerance's reference, by relref */
        *out,                   /* OUT, n: the new running maximum (relref not
                                   POINTLOCAL), which the caller keeps as
                                   `_ref_running`; read on LTE_DONE only */
        *tq;                    /* n: the tolerance, then |lte| / it */
    int64_t *ipiv;              /* m: dgesv's pivots */
    void *dgesv;                /* [ctx] numpy's own `scipy_dgesv_64_` (dgesv64_t) */
    double err;                 /* OUT: the largest normalised error; written on
                                   LTE_DONE only */
    long status;                /* OUT: hdl_fn's result, enum lte_status */
    /* THE SPARSE SOLVE: where `klu_sym` is set, J^-1 Eg is
       `StructuralKLU.solve`'s -- the reduced J row-major
       (`R`, remove_row_col's copy), its pattern's values (`kAx`, by `ktake`:
       a matrix with a nonzero off the pattern is handed back, its fallback
       solver is the Python's), `klu_refactor` on the run's numeric
       (`klu_num`), `klu_solve`, and the refactor's validation -- max|R x - b|
       / max(max|b|, 1e-300) <= `resid_tol`, R x numpy's `dgemv`; a failure is
       handed back too (the Python alone factors again).  Set by
       `_tran_lte_c._klu_bind` once a binding (the buffers are the binding's
       own), but `klu_num` and `resid_tol`: `max_error` sets both at each
       call, `_tran_loop_c` both once a run and `klu_num` again at each loop
       entry.  `klu_sym` NULL (`_Ctx`; `max_error` for a run with no sparse
       solver) selects the dense path */
    long knnz;                  /* the pattern's entries */
    int *kAp, *kAi;             /* the pattern in CSC (StructuralKLU's Ap, Ai) */
    const int64_t *ktake;       /* knnz: each pattern entry's place in R */
    double *kAx,                /* work, knnz (at least 1): the pattern's values,
                                   KLU's Ax */
        *R,                     /* work, m x m: the reduced J row-major */
        *x,                     /* work, m: the solution, KLU's in place */
        *y,                     /* work, m: R x, for the validation */
        resid_tol;              /* the validation's bound,
                                   StructuralKLU.REFACTOR_RESIDUAL_TOL */
    void *klu_sym,              /* KLU's symbolic analysis; NULL: the dense path */
        *klu_num,               /* the run's KLU numeric, refactored in place */
        *klu_common,            /* the solver's KLU common block */
        *klu_refactor, *klu_solve,  /* KLU's entries (klu_refactor_t, klu_solve_t) */
        *dgemv;                 /* numpy's own cblas dgemv (dgemv64_t), for R x */
} lte_t;
long hdl_fn(lte_t *s);

/* hdl_fn's result, left in `status` too (`_tran_lte_c.BAIL` names the
   bails from these): */
enum lte_status {
    LTE_DONE = 1,                /* `err` and `out` are the answer */
    LTE_BAIL_NONFINITE = 2,      /* an input not finite */
    LTE_BAIL_FLAGS = 3,          /* a floating-point flag numpy would raise on */
    LTE_BAIL_DGESV = 4,          /* the divided difference's solve: singular */
    LTE_BAIL_LTE = 5,            /* the estimate not finite */
    LTE_BAIL_TOLERANCE = 6,      /* a tolerance not positive and finite */
    LTE_BAIL_FOREIGN = 7,        /* the sparse solve: a nonzero off the pattern */
    LTE_BAIL_KLU = 8,            /* the sparse solve: KLU's refactor or solve failed */
    LTE_BAIL_RESIDUAL = 9        /* the sparse solve: the refactor's residual too large */
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
