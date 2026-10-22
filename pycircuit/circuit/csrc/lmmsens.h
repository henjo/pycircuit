/* lmmsens.h -- the dense multistep walk's monodromy step: its state and
 * entry, the interface `shooting/_lmmsens_c` loads lmmsens.c through.  cffi
 * lays the struct out from this text and lmmsens.c includes it, so the
 * compiler checks the entry against it.
 *
 * One struct per walk object and size: `_lmmsens_c._Ctx` makes it and its
 * work buffers (kept on the walk as `_lmmsens_c`, keyed (m, its widest w)),
 * `_lmmsens_c.step` fills the per-call fields.  Sizes: m the reduced
 * unknowns, w the monodromy's live columns this step (w <= wcap), L the
 * companion's coefficients (2 or 3). */
/* numpy's own `scipy_cblas_dgemm64_` (CBLAS, 64-bit integers): the one
   unit that multiplies matrices, so its type is declared here */
typedef void (*dgemm64_t)(int, int, int, int64_t, int64_t, int64_t, double, const double *,
                          int64_t, const double *, int64_t, double, double *, int64_t);
typedef struct {
    /* wiring, set once by `_Ctx.__init__`: numpy's own `scipy_dgesv_64_`
       (dgesv64_t) and `scipy_cblas_dgemm64_` (dgemm64_t), by address */
    void *dgesv, *dgemm;
    /* the sizes: m, and the work buffers' widest w */
    long m, wcap;
    /* per call: the live columns, the coefficients alpha[0 .. L-1] and b */
    long w, L;
    double alpha[3], b;
    /* inputs, the caller's arrays, row-major: P0 (m x w) the newest
       monodromy, P1 (m x w) the one before, Pq (m x w) the companion
       derivative's (read where b != 0); C1, C2 (m x m) the capacitances the
       recursion reads them through, Cn (m x m) the step's own, Jf (m x m)
       the step's Jacobian */
    const double *P0, *P1, *Pq, *C1, *C2, *Cn, *Jf;
    /* outputs, fresh arrays the glue sets each call, row-major (m x w):
       Pn the new monodromy, Pqn the new companion derivative */
    double *Pn, *Pqn;
    /* work, `_Ctx`'s: Ck (m x m) a capacitance's nonzero columns, Pk
       (m x wcap) the monodromy's rows they read, T (m x wcap) a product, S
       (m x wcap) the right-hand side, A (m x m, column-major) dgesv's
       matrix then its LU, X (m x wcap, column-major) its right-hand side
       then the solution, ipiv (m) its pivots */
    double *Ck, *Pk, *T, *S, *A, *X;
    int64_t *ipiv;
    /* the nonzero columns of the step's capacitances, kept from call to
       call in `_Ctx`'s buffer: kcs (3 m) C1's, C2's and Cn's lists, nks
       their counts.  reuse, per call: bit 1 -- C1 is the last call's Cn,
       bit 2 -- C2 is the last call's C1 (the same array: the walk never
       writes one), their lists taken from there instead of read again */
    int64_t *kcs;
    long nks[3];
    long reuse;
    /* outputs: status (an lmmsens_status, also pcx_lmmsens' result); info
       dgesv's info -- not written where the call bails before it */
    long status, info;
} lmmsens_t;
long pcx_lmmsens(lmmsens_t *t);

/* pcx_lmmsens' result, left in `status` too (`_lmmsens_c.BAIL` names the
   bails).  A bail hands the step back: the glue reads the status alone and
   the Python step runs on the same state: */
enum lmmsens_status {
    LMMSENS_DONE = 1,
    LMMSENS_BAIL_SINGULAR = 2,      /* dgesv's info: numpy's solve raises */
    LMMSENS_BAIL_NONFINITE = 3      /* a result not finite: the Python step's full-width fallback */
};
