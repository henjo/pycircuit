
#include <stdint.h>
#include <math.h>
#include <fenv.h>
#include "radau_fold.h"

#define FOLD_FLAGS (FE_OVERFLOW | FE_UNDERFLOW | FE_INVALID | FE_DIVBYZERO)

/* numpy's complex product, in the form its loops take here (`cmul_mode`) */
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
   nonzero where the pattern has one and nowhere else */
long hdl_fn(radau_fold_t *f)
{
    const long m = f->m, mm = m * m;
    long i, j, k;
    f->same = 0;
    for (k = 0; k < mm; k++)
        if (!isfinite(f->Cr[k]) || !isfinite(f->Gr[k])) {
            f->status = 2;
            return 0;
        }
    feclearexcept(FE_ALL_EXCEPT);
    for (k = 0; k < mm; k++) {
        const double t = f->s0 * f->Cr[k];
        f->rf[k] = t + f->Gr[k];
    }
    for (k = 0; k < mm; k++) {
        double re, im;
        cmul(f->fused, f->c1r, f->c1i, f->Cr[k], 0.0, &re, &im);
        f->cf[2 * k] = re + f->Gr[k];
        f->cf[2 * k + 1] = im + 0.0;
    }
    f->status = fetestexcept(FOLD_FLAGS) != 0;
    feclearexcept(FE_ALL_EXCEPT);
    if (f->status || !f->nz) return 0;
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
