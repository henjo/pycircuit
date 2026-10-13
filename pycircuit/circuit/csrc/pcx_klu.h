/* pcx_klu.h -- the run's StructuralKLU as a unit holds it: its pattern,
 * its analysis and common, KLU's entries, its numeric and its counters.
 * Included BEFORE a unit's own header (which embeds the struct; cffi reads
 * this text before that one), with the sequence of `pcx_kluseq.h`.
 *
 * THE NUMERIC.  `klu_num` is StructuralKLU's `_numeric` at entry, kept
 * unfreed until the solve ends (`pcx_klu_finish`): a hand-back leaves it as
 * it was -- a refactor in place changes its values, never its pivot order,
 * and every use of it starts with a refactor -- so the Python's own solve
 * repeats the same refactors on the same pivot order.  A served outcome
 * leaves `klu_out`, the numeric the Python's sequence ends with (the
 * entry's freed where another replaced it; NULL after a failed factor, as
 * `_numeric` is None after its LinAlgError). */

typedef struct {
    int *kAp, *kAi;       /* CTX: StructuralKLU's CSC pattern */
    double *kAx;          /* OWN: the matrix's values in its order */
    void *klu_sym, *klu_common;   /* CTX: its analysis (NULL: the dense solve), its common */
    void *klu_factor, *klu_refactor, *klu_solve, *klu_free;   /* CTX: KLU's entries */
    double resid_tol;     /* CTX: StructuralKLU.REFACTOR_RESIDUAL_TOL */
    void *klu_num;        /* CALL: its `_numeric` at entry (NULL: none yet) */
    void *klu_out;        /* OUT: its `_numeric` after; through the solve the one
                             the Python's sequence holds now */
    long refactors, factors, pivot_fallbacks, residual_fallbacks;   /* OUT */
} pcx_klu_t;
