/* radau_fold.h -- the transformed Radau fold's state and entry: the interface `_tran_radau_tc` loads
 * radau_fold.c through.  cffi lays the struct out from this text and
 * radau_fold.c includes it, so the compiler checks the entry against it.
 *
 * One struct a process (`_tran_radau_tc._fold_driver`), its buffers kept
 * for the last size (`_fold_bufs`); `fold` fills the per-call fields and
 * copies the outputs out.  m x m arrays are C order (row-major); a complex
 * array is (re, im) pairs of doubles. */
typedef struct {
    /* m, the factors' size: set with the buffers (`_fold_bufs`); fused,
       numpy's complex product's form (`cmul_mode`: 1 fused, 0 two products
       and a sum): set once (`_fold_driver`); nnz, the last record's
       nonzeros, and nz_at, nz's layout (1: laid out as `A.T`, 0: as `A`):
       set by `fold` where it passes a pattern, else stale and not read */
    long m, fused, nnz, nz_at;
    /* per call (`fold`): the real factor's scalar s0 = lam0.real / h, the
       complex one's c1 = lam1 / h as c1r + i c1i */
    double s0, c1r, c1i;
    /* inputs (m x m): the glue's buffers, set with them (`_fold_bufs`),
       `Cr` and `Gr` copied in each call */
    const double *Cr, *Gr;
    /* per call: the last record's pattern, `A.T != 0` (m x m bytes,
       `ComplexKLUSolver._csc_last`), its handle made once a pattern
       (`_fold_nz`); NULL where there is none (no `ax` made) */
    const uint8_t *nz;
    /* outputs, the glue's buffers (`_fold_bufs`): rf (m x m) the real
       factor; cf (m x m complex) the complex one; ax (room for m * m
       complex) the complex one's nonzeros in the record's order, the first
       nnz read where `same` */
    double *rf, *cf, *ax;
    /* outputs, set every call: status (a radau_fold_status); same, 1 where
       the nonzeros fall exactly on the pattern (ax then complete) */
    long status, same;
} radau_fold_t;
long pcx_radau_fold(radau_fold_t *f);

/* `status` after pcx_radau_fold (which returns 0): */
enum radau_fold_status {
    FOLD_STATUS_MADE = 0,
    FOLD_STATUS_FLAGGED = 1,     /* a flag `errstate(all='raise')` raises on */
    FOLD_STATUS_NONFINITE = 2    /* C or G not finite: nothing made */
};
