/* newton.h -- the C Newton's record and entry: the interface `_tran_newton_c`
 * loads newton.c through.  cffi lays the struct out from this text and
 * newton.c includes it, so the compiler checks the entry against it.
 *
 * One `newton_t` per `_tran_newton_c._Ctx` (one evaluation core, one size
 * n, one reference row iref).  Sizes: n is the full system, m = n - 1 the
 * reduced one (the reference row and column removed, as
 * `analysis.remove_row_col` removes them); an n x n matrix is row-major, as
 * numpy's C-ordered arrays.  A buffer is a numpy array the `_Ctx` allocates
 * once and keeps (`_Ctx.buf`; its n x n ones start at zero) unless its
 * field says otherwise.
 *
 * WHO SETS WHAT -- the tag that opens a field's comment:
 *   CTX   `_tran_newton_c._Ctx.__init__`, once a record;
 *   WALK  `_tran_newton_c._walk_ready` (`solve` clears walk_fn and wn where
 *         nothing limits);
 *   KLU   `_tran_newton_c._klu_bind` (`solve` clears klu_sym on a dense
 *         choice);
 *   CALL  `_tran_newton_c.solve`, before every entry;
 *   OWN   newton.c, during the call: its scratch, its outputs, the status.
 * LOOP marks the 16 fields csrc/loop.c sets before each of its attempts:
 * formula, a0, a1, a2, h, theta (0.0), q1, q2, iq1, u, xseed, have_cache,
 * cx, cC, Cj, Geqj.  It points those pointers at its own buffers (its
 * rings, the source table's output, the predictor's output or its last
 * point, its two C and Geq buffers); `_tran_loop_c` saves the ten pointer
 * and flag ones (q1, q2, iq1, u, xseed, cx, cC, Cj, Geqj, have_cache)
 * before its call and restores them after, and `solve` sets the six
 * scalars anew every call.  Every other field the loop takes as the last
 * served `solve` left it (T, maxiter, reltol, chord, fold_j, abstol, xtol,
 * the walk, the sparse chord). */
typedef struct {
    /* THE CORE (CTX): `_tran_core._Core.cs` and core.c's entry */
    const void *core;     /* the core's hdl_core_t */
    void *core_fn;        /* core.c's `pcx_core` (core_fn_t, pcx_abi.h) */
    /* THE COMPANION (CALL; LOOP sets all but T): the integrator
       `check_order_drop` chose and its scalars, as `_tran_core.evaluate`
       sets them up */
    double T;             /* the temperature (`_Core.probe`'s), written into
                             each element's parameter row at its slot */
    long formula;         /* PCX_FORMULA_* (pcx_abi.h) */
    double a0, a1, a2,    /* gear2's `bdf2_alphas` (0.0 for the others) */
           h, theta;      /* the step; theta's `theta_at(h)` (else 0.0) */
    /* THE HISTORY AND THE SOURCES (CALL: CTX buffers `solve` copies into;
       LOOP points them at its rings and its source table's output) */
    const double *q1,     /* (n) the charge one step back */
                 *q2,     /* (n) two steps back: read by gear2 alone */
                 *iq1,    /* (n) the companion current one step back: read
                             by trapezoidal and theta alone */
                 *u;      /* (n) the sources at the step's end */
    /* THE LIMITER WALK (WALK): walk.c's entry and the tables of the
       circuit's walk object (`_hdl_climit`: these point into its arrays,
       not the record's); walk_fn NULL and wn 0 where nothing limits */
    void *walk_fn;        /* walk.c's `pcx_walk` (walk_fn_t, pcx_abi.h) */
    void **wF;            /* (wn) each element's C limiter (lim_fn_t; 0
                             stops the walk there) */
    double **wPR;         /* (wn) its parameter row */
    const int64_t *wTI,   /* (wn) the temperature's slot in it (-1: none) */
                  *wNM,   /* the node maps, concatenated */
                  *wOFF,  /* (wn) its node map's offset in wNM */
                  *wK;    /* (wn) its node count */
    long wn;              /* the walk's elements */
    /* LAPACK AND BLAS (CTX, from `_tran_newton_c.driver`): numpy's own
       `scipy_dgesv_64_` (numpy.linalg.solve's: the full Newton) and
       `scipy_cblas_dgemv64_` (np.dot's), SciPy's `scipy_dgetrf_` and
       `scipy_dgetrs_` (lu_factor's and lu_solve's: the dense chord) */
    void *dgesv, *dgetrf, *dgetrs, *dgemv;
    long n,               /* CTX: the full size (n >= 2) */
         iref,            /* CTX: the reference row, 0 <= iref < n */
         maxiter,         /* CALL: `par.maxiter` */
         chord;           /* CALL: 1 the chord (`ChordNewton`), 0 the full
                             Newton (`StandardNewton`) */
    double reltol;        /* CALL: `par.reltol` */
    const double *abstol, *xtol;  /* (m) CTX buffers: `_newton_tolerances()`
                                     [2] and [3], copied in by `solve` when
                                     it reads another tuple than the last */
    const double *xseed;  /* CALL, LOOP: (n) the seed, full */
    /* `_C_lookup`'s cached state (CALL: the transient's `_C_cache`, the
       last evaluated point; LOOP: its last attempt's) */
    const double *cx,     /* (n) its x */
                 *cC;     /* (n x n) its C, read only: `solve` points it at
                             the transient's `_C_cache` array itself (n x n
                             >= INPLACE_MIN_ENTRIES and C-ordered) or at the
                             record's copy */
    long have_cache;      /* CALL, LOOP: 1 where cx and cC hold a state */
    /* THE EVALUATIONS' OUTPUTS (OWN, CTX buffers; C, Geq, J n x n and
       q, iq, F n): the first iteration's (C0 .. J0) apart from the later
       ones' (C .. J) -- a chord's later residuals read C0 as the held C
       where the seed missed the cache */
    double *C0, *q0, *iq0, *Geq0, *F0, *J0;
    double *C, *q, *iq, *Geq, *F, *J;
    /* NEWTON.C'S SCRATCH (OWN, CTX buffers) */
    double *xf,           /* (n) the full x of an evaluation (insert_row) */
           *xnf, *x0f;    /* (n) the walk's full x_next (limited in place)
                             and x */
    double *A,            /* (m x m, column-major) the reduced J as numpy's
                             solve and SciPy's lu_factor copy it for LAPACK,
                             then its LU in place (the chord's: held) */
           *absJ;         /* (m x m, row-major) |J| reduced, the convergence
                             test's (the chord's: held from its first
                             iteration); zero where nothing writes it */
    double *b,            /* (m) -F reduced, then the solve's dx in place */
           *Fr,           /* (m) F reduced */
           *x,            /* (m) the iterate: the reduced seed at entry, the
                             solution on NEWTON_CONVERGED (`solve` reads it) */
           *xn,           /* (m) x + dx, limited where the walk runs */
           *xd,           /* (m) the step taken: xn - x (limited), else dx */
           *ax,           /* (m) |xn| */
           *Is;           /* (m) |J| |xn| (dgemv) */
    int64_t *ipiv64;      /* (m) dgesv's pivots */
    int *ipiv32;          /* (m) dgetrf's, held for the chord's dgetrs */
    /* OUTPUTS AND STATUS (OWN) */
    long iters,           /* the iterations made: the Python's count */
         status;          /* the result (enum newton_status), as returned */
    /* THE CONVERGED POINT: `jacobian_only` there, in the same call.  xc ..
       Jj are written only on NEWTON_CONVERGED with fold_j nonzero, and hold
       until the next Newton call (`solve` copies them out; loop.c copies
       xc, qj, iqj and alternates its own two Cj and Geqj buffers) */
    long fold_j;          /* CALL: its passes, a PCX_PASS_* mask, not a
                             flag: G | C (3) where the step's J is read, C
                             (2) where nothing reads it; 0: no evaluation */
    double *xc,           /* OWN (n): the converged point, full: insert_row
                             of x, +0.0 at iref */
           *Cj,           /* OWN (n x n): C there (LOOP points it at one
                             of its own two) */
           *qj,           /* OWN (n): q there */
           *iqj,          /* OWN (n): the companion's iq there */
           *Geqj,         /* OWN (n x n): its Geq there (LOOP likewise) */
           *Jj;           /* OWN (n x n): J = G + Geq there, written only
                             where fold_j has PCX_PASS_G */
    long pattern;   /* CTX: PCX_WANT_PATTERN where the core writes only the
                       pattern's entries (these buffers zero off it), else 0 */
    /* THE SPARSE CHORD (KLU): where `klu_sym` is set, the chord factors and
       solves with KLU on the circuit's structural pattern --
       `StructuralKLU`'s own analysis and common, so the factors are the
       ones its `factor` makes -- in place of LAPACK.  (The chord only:
       `_tran_newton_c` declines a full Newton on KLU.)  The arrays are the
       StructuralKLU's own or the bind's (`_Ctx.klu_keep`). */
    long knnz;            /* the pattern's entries */
    int *kAp, *kAi;       /* the pattern in CSC (StructuralKLU's Ap, Ai) */
    const int64_t *kat,   /* (knnz) each entry's place in the core's full
                             n x n J */
                  *kabs;  /* (knnz) and in the reduced, row-major |J| (absJ:
                             zero off the pattern) */
    double *kAx;          /* (knnz) OWN: the entries' values, J's at the
                             first iteration */
    void *klu_sym,        /* its symbolic analysis; NULL: LAPACK's path */
         *klu_common,     /* its klu_common (KLU's calls write it) */
         *klu_factor,     /* KLU's klu_factor (klu_factor_t, pcx_abi.h) */
         *klu_solve,      /* klu_solve (klu_solve_t) */
         *klu_free;       /* klu_free_numeric (klu_free_numeric_t) */
} newton_t;
/* THE SOLVE: `_newton`'s Newton on the reduced system -- `StandardNewton`
   (chord 0), or `ChordNewton` up to its fallback (chord 1) -- and
   `jacobian_only` at the converged point (fold_j), bit for bit.
   Assumes every CTX and CALL field set (WALK and KLU where used), the
   core ready for this call (`_Core.probe`, the caller's), and the n x n
   buffers zero off the circuit's pattern where `pattern` is set.
   Returns the status, also left in `status`.  NEWTON_CONVERGED: x, iters
   and (fold_j) the converged point's fields written.  A bail writes none
   of xc .. Jj: only the OWN scratch, x, iters, status and the other OWN
   fields.  Either way the call writes, outside the record, what the
   Python's evaluation writes -- each element's parameter row at its
   temperature slot (the core's and the walk's kernels) -- the core's own
   scratch, and on the sparse chord KLU's common; the numeric it factors
   it frees before it returns. */
long pcx_newton(newton_t *s);

/* pcx_newton's result, left in `status` too (`_tran_newton_c.BAIL` names the
   bails from these).  Every code from 2 is a hand-back: `_tran_newton_c`
   rolls the attempt's trace back and the Python Newton repeats the solve
   from the same seed (loop.c stops, LOOP_BAIL_NEWTON: the attempt is the
   Python's). */
enum newton_status {
    NEWTON_RUNNING = 0,          /* in `status` while it runs; never returned */
    NEWTON_CONVERGED = 1,
    NEWTON_BAIL_X = 2,           /* an iterate not finite */
    NEWTON_BAIL_JNONFINITE = 3,  /* the chord's Jacobian not finite */
    NEWTON_BAIL_GETRF = 4,       /* the chord's factorisation: singular */
    NEWTON_BAIL_BNONFINITE = 5,  /* a chord step's right-hand side not finite */
    NEWTON_BAIL_GETRS = 6,       /* the chord's solve failed */
    NEWTON_BAIL_DGESV = 7,       /* a full step's solve: singular */
    NEWTON_BAIL_WALKSTOP = 8,    /* a limiter declined: the walk stopped */
    NEWTON_BAIL_XNEXT = 9,       /* the next iterate not finite */
    NEWTON_BAIL_CHORDBREAK = 10, /* the chord's steps stopped shrinking */
    NEWTON_BAIL_MAXITER = 11,    /* no convergence in maxiter */
    NEWTON_BAIL_KLU_FACTOR = 12, /* the sparse chord's factorisation: singular */
    NEWTON_BAIL_KLU_SOLVE = 13   /* the sparse chord's solve failed */
};
