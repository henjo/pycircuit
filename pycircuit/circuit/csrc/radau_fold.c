/* radau_fold.c -- `_RadauStages._radau_frozen`'s two factors in one call:
   the real one (lam0.real / h) Cr + Gr and the complex one (lam1 / h) Cr +
   Gr, numpy's bytes and the flags its `errstate(all='raise')` raises on,
   and the complex one's nonzeros where they fall in the last record's
   pattern.  The glue, `_tran_radau_tc.fold`, fills radau_fold_t
   (radau_fold.h) and returns what it read; every decision stays in
   `_radau_frozen`. */
#include <stdint.h>
#include <math.h>
#include <fenv.h>
#include "pcx_fenv.h"
#include "radau_fold.h"

/* the four flags `errstate(all='raise')` raises on (inexact is not one).
   MUST KEEP: every result is stored to rf or cf before the one
   fetestexcept.  The compiler does not model the flags (no FENV_ACCESS),
   so an operation whose result stays in a register may move past the
   test; one whose result is stored to memory the call could read may
   not. */

/* numpy's complex product, in the form its loops take here (`cmul_mode`):
   fused as `cmul_mode` reads it from numpy (re rounds ai bi first and
   fuses ar br into it, im rounds ai br first and fuses ar bi), or two
   products and a sum */
static void cmul(long fused, double ar, double ai, double br, double bi, double *re, double *im)
{
    if (fused) {
        *re = fma(ar, br, -(ai * bi));
        *im = fma(ar, bi, ai * br);
    } else {
        *re = ar * br - ai * bi;
        *im = ar * bi + ai * br;
    }
}

/* `_radau_frozen`'s factors as numpy makes them: the real one (s0 * Cr) +
   Gr, two roundings as its two loops; the complex one c1 * Cr -- Cr cast to
   complex, numpy's product -- + Gr cast to complex (the imaginary part
   + 0.0).  `status`: 0 made; 1 a flag `errstate(all='raise')` raises on,
   in either; 2 a non-finite input, none made (a NaN's payload follows
   numpy's loops, not the operands' order: numpy's to make).  Finite inputs
   and no flag: both factors finite.  Then, given the last record's pattern
   (`nz`: `A.T != 0`, laid out as `A.T` -- `nz_at` -- or as `A`), the complex
   one's nonzeros column by column, rows ascending, into `ax`; `same`: every
   nonzero where the pattern has one and nowhere else.
   Assumes the buffers of radau_fold.h's sizes for `m`.  Writes `same` (0
   unless the nonzeros fall exactly on the pattern) and `status` on every
   call; rf and cf unless status 2 (read by the glue only on status 0); ax
   only past status 0 with a pattern, partly where `same` stays 0 (not
   read then).  Returns 0 always: the outcome is `status`. */
long hdl_fn(radau_fold_t *f)
{
    const long m = f->m, mm = m * m;
    long i, j, k;
    f->same = 0;
    for (k = 0; k < mm; k++)
        if (!isfinite(f->Cr[k]) || !isfinite(f->Gr[k])) {
            f->status = FOLD_STATUS_NONFINITE;
            return 0;
        }
    feclearexcept(FE_ALL_EXCEPT);
    /* (s0 * Cr) + Gr: the product rounded, then the sum -- numpy's two
       loops, never one fused operation */
    for (k = 0; k < mm; k++) {
        const double t = f->s0 * f->Cr[k];
        f->rf[k] = t + f->Gr[k];
    }
    for (k = 0; k < mm; k++) {
        double re, im;
        /* c1 * (Cr + 0j), then + (Gr + 0j): the product on the left, as the
           Python's `c1 * Cr + Gr` */
        cmul(f->fused, f->c1r, f->c1i, f->Cr[k], 0.0, &re, &im);
        f->cf[2 * k] = re + f->Gr[k];
        /* Gr's imaginary part is +0.0: added, not skipped -- a -0.0
           product becomes +0.0, as numpy's complex sum makes it */
        f->cf[2 * k + 1] = im + 0.0;
    }
    f->status = fetestexcept(PCX_FP_FLAGS) != 0;
    feclearexcept(FE_ALL_EXCEPT);
    if (f->status || !f->nz) return 0;
    /* the CSC gather: a nonzero is `!= 0` in either part (a -0.0 is
       zero), as the pattern's mask tests it */
    k = 0;
    for (j = 0; j < m; j++)
        for (i = 0; i < m; i++) {
            const double re = f->cf[2 * (i * m + j)], im = f->cf[2 * (i * m + j) + 1];
            const long nz = re != 0.0 || im != 0.0;
            const long was = f->nz[f->nz_at ? j * m + i : i * m + j] != 0;
            if (nz != was) return 0;
            if (nz) {
                if (k == f->nnz) return 0;
                f->ax[2 * k] = re;
                f->ax[2 * k + 1] = im;
                k++;
            }
        }
    f->same = k == f->nnz;
    return 0;
}
