/* sensklu.h -- the shooting's coupled stage step map on the walk's KLU: its
 * state and entry, the interface `_sensklu_c` loads sensklu.c through.  cffi
 * lays the struct out from this text and sensklu.c includes it, so the
 * compiler checks the entry against it.
 *
 * One struct per StructuralKLU a walk factors its blocks on (`_sensklu_c._Ctx`,
 * kept on the solver), its buffers the glue's; `_sensklu_c.block_step` fills
 * the per-call fields.  Sizes: M = s m the stage system, m the reduced
 * unknowns, s the stages, k the monodromy's columns, nnz the pattern's
 * entries. */
typedef struct {
    /* THE KLU (`_Ctx.__init__`, once a solver): StructuralKLU's own pattern,
       analysis and common, KLU's entries, `REFACTOR_RESIDUAL_TOL` */
    int *kAp, *kAi;       /* the pattern in CSC (StructuralKLU's Ap, Ai) */
    const int64_t *cols;  /* (nnz) each entry's column (`_KLUWalk._cols`) */
    void *klu_sym,        /* the analysis (`_symbolic`) */
         *klu_common,     /* the common KLU's calls write */
         *klu_factor,     /* klu_factor (klu_factor_t) */
         *klu_refactor,   /* klu_refactor (klu_refactor_t) */
         *klu_solve,      /* klu_solve (klu_solve_t) */
         *klu_free;       /* klu_free_numeric (klu_free_numeric_t) */
    double resid_tol;     /* the refactor check's tolerance */
    long M, nnz;          /* the block's size and its pattern's entries */
    /* PER CALL (`block_step`) */
    double *Ax;           /* (nnz) the block's values in the pattern's order
                             (KLU's entries take them unconst; none writes) */
    void *num;            /* in: the walk's numeric (`_KLUWalk._num`, NULL:
                             none yet); out: the one it holds after */
    const double *base;   /* (m x k, row-major) C_n P, numpy's product */
    long m, k, s;         /* as above, M = s m */
    double *X;            /* (M x k, column-major) out: every stage's block
                             `base`, solved -- its last stage block the new P */
    /* SCRATCH (`_Ctx`'s buffers, M each): the refactor check's probe */
    double *pb,           /* the block times ones */
           *px,           /* its solution */
           *pr;           /* the block times that */
    /* OUT, this call's: `_KLUWalk`'s counts and the result */
    long refactors, fallbacks, factors;
    long status;          /* a sensklu_status, also pcx_sensklu's result */
} sensklu_t;
long pcx_sensklu(sensklu_t *t);

/* pcx_sensklu's result, left in `status` too (`_sensklu_c.BAIL` names the
   bails from these).  SENSKLU_SINGULAR and SENSKLU_SOLVE_FAILED end the step
   as the Python's does -- `num` and the counts made, the glue raising its
   LinAlgError.  A bail hands the step back before any count: the glue reads
   `num` alone (a refactor in place leaves its pivot order, so the Python's
   repeat makes the same numeric) and the Python step runs: */
enum sensklu_status {
    SENSKLU_DONE = 1,
    SENSKLU_SINGULAR = 2,        /* a fresh factor failed: `num` NULL */
    SENSKLU_SOLVE_FAILED = 3,    /* klu_solve failed on the stage system */
    SENSKLU_BAIL_FLAGS = 4       /* a flag numpy's arithmetic of the check would warn of */
};
