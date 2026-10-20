/* radau_tc.h -- the transformed Radau stage solve's state and entry: the interface `_tran_radau_tc` loads
 * radau_tc.c through.  cffi lays the struct out from this text and
 * radau_tc.c includes it, so the compiler checks the entry against it.
 *
 * One struct per transient, core, size and reference row: `_tran_radau_tc._Ctx`
 * makes it and its buffers (kept on the transient as `_radau_tc`), `solve` fills
 * the per-call fields, writes the factorisations' state back and reads the
 * outputs.  Sizes: n the full unknowns, m = n - 1 the reduced ones (row
 * `iref` dropped).  Arrays are C order (row-major) unless said otherwise;
 * a complex array is (re, im) pairs of doubles. */
typedef struct {
    /* wiring, set once by `_Ctx.__init__`: the evaluate core's struct
       (`core.cs`) and its entry (core.c's `pcx_core`, a core_fn_t, called with
       formula 4: the passes alone) */
    const void *core;
    void *core_fn;
    /* per call (`solve`): the temperature the core and the walk evaluate
       at (`core.probe`) */
    double T;
    /* wiring, set once: the core's own i bin (formula 4 leaves i there,
       n); max(n*n, 1) zeros for the core's arguments a call neither reads
       nor writes; `bits`, the iterations' want bits (`_tran_core._PASS_BITS`
       'q' | 'i': i -- q is always made) */
    const double *bini;
    double *dummy;
    long bits;
    /* wiring: the limiter walk (walk.c's `pcx_walk`, a walk_fn_t, and its
       tables; `wn` its elements), set by `_tran_newton_c._walk_ready` where
       the circuit limits; walk_fn NULL and wn 0 where it does not */
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    /* wiring: numpy's own OpenBLAS `dgesv` and `dgetrs` (64-bit integers:
       dgesv64_t, dgetrs64_t), set once (`linearsolver._numpy_lapack`); the
       complex KLU's `klu_z_refactor` and `klu_z_solve` (klu_refactor_t,
       klu_solve_t), set when the solver's library or handles change */
    void *dgesv; void *dgetrs; void *klu_refactor; void *klu_solve;
    /* the complex factor's `ComplexKLUSolver` handles (`_symbolic`,
       `_numeric`, its `_common`), set when they change (`rec.klu_kk`).
       (newton.h and lte.h name their KLU bindings klu_sym, kAp, ...) */
    void *Symbolic; void *Numeric; void *Common;
    /* the step's record (`_FrozenTransform.prep`): the complex factor in
       CSC, Ap (m + 1) and Ai (nnz) set when those arrays change, Ax (nnz
       complex) every call -- what KLU refactors and the residual check
       multiplies */
    int *Ap; int *Ai; double *Ax;
    /* n, iref (0 <= iref < n, n >= 2) set once by `_Ctx.__init__`; maxiter
       per call; fused set once: numpy's complex product's form (`cmul_mode`:
       1 fused, 0 two products and a sum) */
    long n, iref, maxiter, fused;
    /* per call: the step; the convergence test's tolerances
       (`_stages_converged`); rtol, the residual check's tolerance
       (`ComplexKLUSolver.REFACTOR_RESIDUAL_TOL`) */
    double h, reltol, abstol, rtol;
    /* the tableau's A, A[3 i + j] = A_ij, set when its bytes change */
    double A[9];
    /* `_FrozenTransform.p`: P's entries P00 P01 P02 P10 P11 P12, complex,
       set when the step's tuple changes */
    double p[12];
    /* `_FrozenTransform.v`, per stage i: V[i,0].real, then V[i,1]'s re and
       im, set with p */
    double v[9];
    /* per-call inputs, the glue's buffers (`_Ctx.buf`) filled by `solve`:
       qn (n) the charge at x_n; u (3 x n) the sources at the three stage
       times; seed (3 x n) the stages the loop starts from */
    const double *qn, *u, *seed;
    /* the real factor's kept LU (`linearsolver.NumpyLU`): lu_a (m x m,
       column-major) and lu_ipiv (m) by their addresses (`lu.addr`), set when
       they change.  lu_info, in and out each call: -1 no factors yet (lu_a
       holds the factor; the first solve's dgesv factors it in place), 0
       factors kept (dgetrs); after a dgesv, its info -- `solve` writes it
       back to `NumpyLU._info`, on a bail too */
    double *lu_a; int64_t *lu_ipiv; int64_t lu_info;
    /* klu_fresh, in: 1 where the numeric already holds this record's
       refactor (`ComplexKLUSolver._fresh is prep`), so none is made; set 1
       by the refactor.  refactored, out: the refactors made (0 or 1) --
       `solve` adds them to `refactors` and sets `_fresh`, on a bail too */
    long klu_fresh, refactored;
    /* work and outputs (`_Ctx.buf`): Y (3 x n) the stages -- the seed
       copied in, the converged stages out; Yp (3 x n) the stages before the
       update; q, iv, K, F (3 x n) the stages' q and i, -(i + u), the full
       residual; R (3 m) the reduced residual; b0 (m) the real right-hand
       side; x0 (m) rhs0's imaginary part (made, then dropped), then the real
       solve w0; r1 (m complex) the complex right-hand side; x1 (m complex)
       its solve w1; y1 (m complex) A w1 and ab (m) sum_j |A_kj| |w1_j|, the
       residual check's; dY (3 m) the update; stp (n) one stage's limited
       step; tolv (2) reltol * ynorm, then plus abstol (stored: see
       radau_tc.c) */
    double *Y, *Yp, *q, *iv, *K, *F, *R, *b0, *x0, *r1, *x1, *y1, *ab, *dY, *stp, *tolv;
    /* outputs, set by pcx_radau_tc on every call: iters (the iterations
       completed), status (a radau_tc_status, also pcx_radau_tc's result) */
    long iters, status;
    /* the step's end passes (`_stage_end_passes`'s calls, made after
       convergence).  binG, wiring set once: the core's own G bin.  end_bits,
       per call: the last stage's want bits (7: G, C, i; 0: none made);
       end_cg, per call: 3 (C, G) at the first two stages where the shooting
       reads them, else 0; end_done, out: 1 where they were made */
    const double *binG; long end_bits, end_cg, end_done;
    /* the end passes' outputs (`_Ctx.buf`), read where end_done: the last
       stage's e_q, e_i (n), e_C, e_G (n x n); the first two stages' e_C01,
       e_G01 (2 x n x n); e_q01 (n) the q those two calls make (q is always
       made), not read */
    double *e_q, *e_i, *e_C, *e_G, *e_C01, *e_G01, *e_q01;
} radau_tf_t;
long pcx_radau_tc(radau_tf_t *s);

/* pcx_radau_tc's result, left in `status` too (`_tran_radau_tc.BAIL` names the
   bails from these).  CONVERGED and MAXITER are the loop's own outcomes,
   served.  A bail hands the loop back: the glue reads the status and the
   factorisations' state (lu_info, refactored), and the Python loop runs
   from the same seed: */
enum radau_tc_status {
    RADAU_TC_RUNNING = 0,
    RADAU_TC_CONVERGED = 1,
    RADAU_TC_BAIL_NONFINITE = 2, /* an input or a stage not finite */
    RADAU_TC_BAIL_FLAGS = 3,     /* a floating-point flag in the C's own arithmetic (TC_FLAGS) */
    RADAU_TC_BAIL_RHS = 4,       /* a real right-hand side not finite (the analysis solver's) */
    RADAU_TC_BAIL_LU = 5,        /* the real block's factorisation: singular (dgesv's info) */
    RADAU_TC_BAIL_REFACTOR = 6,  /* KLU's refactorisation failed */
    RADAU_TC_BAIL_KLUSOLVE = 7,  /* KLU's solve failed */
    RADAU_TC_BAIL_RESIDUAL = 8,  /* the complex block's residual check failed or not decided */
    RADAU_TC_BAIL_STEP = 9,      /* a step not finite */
    RADAU_TC_BAIL_WALKSTOP = 10, /* a limiter declined: the walk stopped */
    RADAU_TC_MAXITER = 11        /* no convergence in maxiter: the Python loop's False */
};
