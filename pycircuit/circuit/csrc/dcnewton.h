/* dcnewton.h -- the DC operating point's base Newton solve in C: its struct,
 * entry and codes, the interface `_dc_newton_c` loads dcnewton.c through
 * (cffi lays the struct out from this text, after pcx_klu.h's; dcnewton.c
 * includes both).
 *
 * One `StandardNewton.solve_system` of the DC's continuation chain, for an
 * evaluation the ladders compose from the DC's own residual -- `cut(i + u)`
 * or source stepping's `cut(i + lam * u)` and `cut(G)` (`dcanalysis`'s
 * `func`, `source_callback`, `refnode_removed`) -- under gshunt (`F + x g`,
 * `J + 0.0` and g on the diagonal: `nrsolver._plus_diagonal`) and junction
 * gmin over it (`nrsolver._junction_gmin`'s per-junction terms, in order).
 * The solve: numpy's dgesv (dense), or StructuralKLU's own sequence (a
 * refactor on the run's numeric, a fresh factor on a zero pivot, klu_solve,
 * the residual against REFACTOR_RESIDUAL_TOL, a fresh factor where it fails);
 * x_next, the limiter walk, `xdiff = x_next - x`, nrsolver's convergence test.
 *
 * THE NUMERIC (pcx_klu.h): a hand-back (any DCNEWTON_BAIL_*) leaves
 * StructuralKLU's numeric as it was at entry; on DCNEWTON_CONVERGED and
 * DCNEWTON_MAXITER, `klu.klu_out` is the numeric the Python's sequence ends
 * with and the counters are this call's. */

typedef struct {
    /* THE CORE (CTX: the record's; T per call) */
    const void *core;     /* the core's hdl_core_t (`_tran_core._Core.cs`) */
    void *core_fn;        /* core.c's `pcx_core` (core_fn_t) */
    long pattern;         /* PCX_WANT_PATTERN where the core writes only its
                             pattern (its binG zero off it), else 0 */
    double T;             /* CALL: the temperature (`_Core.probe`'s) */
    /* THE RESIDUAL (CALL) */
    const double *u;      /* (n) `cir.u(0, analysis='dc')`, full */
    long has_lam;         /* 1: source stepping's `i + lam * u`; 0: `i + u` */
    double lam;           /* its factor */
    /* THE LADDERS' LAYERS (CALL): junction gmin over gshunt over the base */
    long has_gshunt;      /* 1: gshunt's `F + x g`, `J + 0.0`, g on the diagonal */
    double gshunt;        /* its conductance */
    long has_jgmin;       /* 1: junction gmin's terms */
    double gjgmin;        /* its conductance */
    long nj;              /* the junctions */
    const int64_t *ja, *jb;   /* (nj) their reduced rows, in the loop's order */
    /* THE LIMITER WALK (CALL: `_tran_newton_c._walk_ready`'s tables, as
       newton_t's); walk_fn NULL and wn 0 where nothing limits */
    void *walk_fn;
    void **wF;
    double **wPR;
    const int64_t *wTI, *wNM, *wOFF, *wK;
    long wn;
    /* LAPACK AND BLAS (CTX): numpy's own `scipy_dgesv_64_` and
       `scipy_cblas_dgemv64_` */
    void *dgesv, *dgemv;
    long n,               /* CTX: the full size (n >= 2) */
         iref,            /* CTX: the reference row */
         maxiter;         /* CALL */
    double reltol;        /* CALL */
    const double *abstol, *xtol;   /* CALL (m): the reduced tolerances */
    /* THE ITERATE (CALL, OWN): x0 at entry; the solution on CONVERGED, the
       last iteration's x on MAXITER */
    double *x;            /* (m) */
    /* THE LAST ITERATION (OWN): what `_worst_row_report` reads on MAXITER */
    double *Fr,           /* (m) F, reduced, deformed */
           *xn,           /* (m) x_next, limited */
           *xd;           /* (m) xdiff = x_next - x */
    /* SCRATCH (CTX buffers) */
    double *xf,           /* (n) insert_row(x) */
           *Ff,           /* (n) F full: `i + u` (or `i + lam * u`) */
           *xnf, *x0f,    /* (n) the walk's x_next and x, full */
           *Cs, *q,       /* (n x n), (n): the core's C and q bins (C unread) */
           *R,            /* (m x m, row-major) J: zero off the active entries */
           *absJ,         /* (m x m, row-major) |J|: zero off them */
           *A,            /* (m x m, column-major) the dense solve's copy */
           *b,            /* (m) -F: the solve's right-hand side; the dense
                             solve's step in place */
           *kx,           /* (m) KLU's solution (the step) */
           *ax, *Is,      /* (m) |x_next|, |J| |x_next| */
           *ky;           /* (m) the residual's J x */
    int64_t *ipiv64;      /* (m) dgesv's pivots */
    /* THE ACTIVE ENTRIES (CTX): where J can be nonzero -- KLU's pattern in
       its CSC order (sparse), every entry (dense) */
    long nact;
    const int64_t *act_r, /* (nact) the entry in the reduced, row-major J */
                  *act_f; /* (nact) its place in the core's full n x n G */
    /* KLU (pcx_klu.h; klu.klu_sym NULL: the dense solve; kAx (nact) J's
       values in its order; klu_out, the counters: on CONVERGED and MAXITER) */
    pcx_klu_t klu;
    /* OUT */
    long iters,           /* the iterations made: the Python's `i + 1` */
         status;          /* enum dcnewton_status, as returned */
} dcnewton_t;

long pcx_dcnewton(dcnewton_t *s);

enum dcnewton_status {
    DCNEWTON_CONVERGED = 1,       /* `x` the solution, `iters` the count */
    DCNEWTON_MAXITER = 2,         /* no convergence in `maxiter`: the last
                                     iteration's F, x_next, xdiff and x */
    DCNEWTON_BAIL_FLAGS = 3,      /* a floating-point flag numpy would warn
                                     of (or a non-finite value): the
                                     Python's solve, its warnings */
    DCNEWTON_BAIL_SINGULAR = 4,   /* a factor or a solve failed: the
                                     Python's solve raises */
    DCNEWTON_BAIL_WALKSTOP = 5    /* a limiter declined (the walk stopped) */
};
