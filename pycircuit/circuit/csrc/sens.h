/* sens.h -- the shooting sensitivities' state and entry: the interface `_sens_c` loads
 * sens.c through.  cffi lays the struct out from this text and
 * sens.c includes it, so the compiler checks the entry against it.
 *
 * One struct per inner transient and size: `_sens_c._Ctx` makes it and its
 * buffers (kept on the transient as `_sens_c`, keyed (s, nf, k, iref)),
 * `_sens_c._map` fills the per-call fields.  Sizes: s stages, nf the full
 * unknowns, m = nf - 1 the reduced ones, n = s m the stage system, k the
 * monodromy's columns. */
typedef struct {
    /* wiring, set once by `_Ctx.__init__`: SciPy's own OpenBLAS
       `scipy_dgetrf_` and `scipy_dgetrs_` (32-bit integers: dgetrf32_t,
       dgetrs32_t), by address (`_sens_c.driver`) */
    void *dgetrf, *dgetrs;
    /* the sizes, set once by `_Ctx.__init__`: s, m, k, nf as above; iref the
       reference row and column the readers drop from the full C and G */
    long s, m, k, nf, iref;
    /* per call (`_map`): the step */
    double h;
    /* inputs, the glue's buffers (`_Ctx`), copied in by `_map`: A (s x s)
       the tableau, A[i s + j] = A_ij, copied when its bytes change; Cf, Gf
       (s x nf x nf) each stage's FULL C and G (the step's device memo);
       base (m x k) C_n P, numpy's product, made in Python */
    const double *A, *Cf, *Gf, *base;
    /* J (n x n, column-major), the stage block, then dgetrf's LU of it;
       ipiv (n) its pivots -- `_Ctx`'s buffers.  Z (n x k, column-major:
       Z[c n + row]) the solution, a fresh array `_map` sets each call and
       clears after; its last stage block is the new P */
    double *J; int *ipiv; double *Z;
    /* outputs: status (a sens_status, also hdl_fn's result), set on every
       return; info, LAPACK's info of the last call made (dgetrf, then
       dgetrs) -- not written where the call bails before dgetrf, not read
       by the glue */
    long status, info;
} sens_t;
long hdl_fn(sens_t *t);

/* hdl_fn's result, left in `status` too (`_sens_c.BAIL` names the bails
   from these).  A bail hands the step back: the glue reads the status
   alone and the Python step runs: */
enum sens_status {
    SENS_DONE = 1,
    SENS_BAIL_NONFINITE = 2,     /* a reduced entry not finite */
    SENS_BAIL_FLAGS = 3,         /* a flag `_stage_block`'s errstate(all='raise') raises on */
    SENS_BAIL_RHS = 4,           /* C_n P not finite */
    SENS_BAIL_SINGULAR = 5,      /* the block's factorisation: singular */
    SENS_BAIL_LAPACK = 6         /* the block's solve failed */
};
