/* radau.h -- the Radau IIA stage solve's state and entry: the interface `_tran_radau_c` loads
 * radau.c through.  cffi lays the struct out from this text and
 * radau.c includes it, so the compiler checks the entry against it.
 *
 * One struct per transient, core, size and reference row: `_tran_radau_c._Ctx`
 * makes it and its buffers (kept on the transient as `_radau_c`), `solve` fills
 * the per-call fields and reads the outputs.  Sizes: n the full unknowns, m = n - 1
 * the reduced ones (row and column `iref` dropped), M = 3 m the coupled system.
 * Arrays are C order (row-major) unless said otherwise. */
typedef struct {
    /* wiring, set once by `_Ctx.__init__`: the evaluate core's struct
       (`core.cs`) and its entry (core.c's `pcx_core`, a core_fn_t, called with
       formula 4: the passes alone) */
    const void *core;
    void *core_fn;
    /* per call (`solve`): the temperature the core and the walk evaluate
       at (`core.probe`) */
    double T;
    /* wiring, set once: the core's own bins, where formula 4 leaves G
       (n x n) and i (n) of the call just made -- copied out straight after
       each stage's call */
    const double *binG, *bini;
    /* wiring, set once: max(n*n, 1) zeros for the core's arguments formula
       4 neither reads nor writes */
    double *dummy;
    /* wiring: the limiter walk (walk.c's `pcx_walk`, a walk_fn_t, and its
       tables; `wn` its elements), set by `_tran_newton_c._walk_ready` where
       the circuit limits; walk_fn NULL and wn 0 where it does not */
    void *walk_fn; void **wF; double **wPR; const int64_t *wTI, *wNM, *wOFF, *wK; long wn;
    /* wiring, set once: numpy's own OpenBLAS `scipy_dgesv_64_` (64-bit
       integers, a dgesv64_t), by address (`_tran_newton_c.driver`) */
    void *dgesv;
    /* n, iref (0 <= iref < n, n >= 2) and maxa (the assembly slots of Ya..Ga:
       `MAXA`) set once by `_Ctx.__init__`; maxiter (the Newton's iterations)
       per call */
    long n, iref, maxiter, maxa;
    /* per call: the step and the convergence test's tolerances
       (`_stages_converged`: the step within reltol of the largest reduced
       stage plus abstol) */
    double h, reltol, abstol;
    /* per call: the tableau's A, A[3 i + j] = A_ij */
    double A[9];
    /* per-call inputs, the glue's buffers (`_Ctx.buf`) filled by `solve`:
       qn (n) the charge at x_n; u (3 x n) the sources at the three stage
       times; seed (3 x n) the stages the Newton starts from */
    const double *qn, *u, *seed;
    /* outputs, maxa slots each (`_Ctx.buf`): slot a holds assembly a's
       stages Ya (3 x n) and their passes qa, ia (3 x n), Ca, Ga (3 x n x n).
       On RADAU_CONVERGED slots 0 .. assemblies-1 are every assembly the
       Python Newton makes, in its order (the device memo's records), and
       slot assemblies-1's Ya the converged stages; on a bail not read */
    double *Ya, *qa, *ia, *Ca, *Ga;
    /* work (`_Ctx.buf`): K (3 x n) -(i + u); F (3 x n) the full residual;
       R (M) the reduced residual; Jf (M x M, column-major: numpy's solve's
       LAPACK copy) the reduced block Jacobian, dgesv's LU after the solve;
       b (M) -R, then the solution dY; blk (n x n) one block before its
       reduction; stp (n) one stage's limited step; tolv (2) reltol * ynorm,
       then plus abstol (stored: see radau.c) */
    double *K, *F, *R, *Jf, *b, *blk, *stp, *tolv;
    /* work: dgesv's pivots (max(M, 1), `_Ctx.buf['ipiv']`) */
    int64_t *ipiv;
    /* outputs, set by pcx_radau on every call: assemblies (0 until
       RADAU_CONVERGED, then the assemblies made), iters (the iterations
       completed), status (a radau_status, also pcx_radau's result) */
    long assemblies, iters, status;
} radau_t;
long pcx_radau(radau_t *s);

/* pcx_radau's result, left in `status` too (`_tran_radau_c.BAIL` names the
   bails from these; 7 and 10 unused).  A bail hands the solve back: the
   glue reads nothing but the status, and the Python Newton runs from the
   same seed: */
enum radau_status {
    RADAU_RUNNING = 0,
    RADAU_CONVERGED = 1,
    RADAU_BAIL_NONFINITE = 2,    /* an input or a stage not finite */
    RADAU_BAIL_FLAGS = 3,        /* a floating-point flag in the C's own arithmetic (RC_FLAGS) */
    RADAU_BAIL_RBOUND = 4,       /* a residual entry not below 2^1000, or a Jacobian entry not finite */
    RADAU_BAIL_DGESV = 5,        /* the stage system's solve: singular (dgesv's info) */
    RADAU_BAIL_ASSEMBLIES = 6,   /* more assemblies than the buffers hold (maxa) */
    RADAU_BAIL_WALKSTOP = 8,     /* a limiter declined: the walk stopped */
    RADAU_BAIL_STEP = 9,         /* a step not finite */
    RADAU_BAIL_MAXITER = 11      /* no convergence in maxiter */
};
