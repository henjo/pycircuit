/* pcx_abi.h -- the function types pycircuit's runtime C units share, declared
 * once.
 *
 * A unit includes it as `#include "pcx_abi.h"`; `_cbuild.csrc` expands that
 * line in place, so the text a unit is keyed by, compiled from and digested
 * by its pinned pair is the text with this one in it: change a type here and
 * every unit using it rebuilds and its pair names it.  Every function below
 * reaches a unit as a `void *` field of its struct and is cast at the call.
 */
#include <stddef.h>
#include <stdint.h>

/* ---- the other units' entries ---------------------------------------- */

/* the evaluation core, core.c's `pcx_core` (its first argument a
   `const hdl_core_t *` there): the elements' C, q, i, G at x, stamped */
typedef long (*core_fn_t)(const void *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
/* the limiter walk, walk.c's `pcx_walk`: each limited element's limiter on xg
   from `start`; returns n, or the element whose limiter declined */
typedef long (*walk_fn_t)(void **F, double **PR, const int64_t *TI, double T,
            const int64_t *NM, const int64_t *OFF, const int64_t *K,
            double *xg, const double *x0g, long start, long n);

/* ---- the kernels the HDL printer emits (hdl.py) ---------------------- */

/* one element's printed pass: its inputs, its parameter row, its outputs */
typedef void (*pass_fn_t)(const double *, const double *, double *);
/* a class's fused kernel, every pass a call wants at once (bind_fused) */
typedef void (*fused_fn_t)(const double *, const double *, long, double *, double *, double *,
                           double *);
/* one element's limiter (_hdl_climit): nonzero when it declines */
typedef int (*lim_fn_t)(const double *, const double *, double *);

/* ---- LAPACK and BLAS, named by the width of their integers: 64 from
        numpy's own OpenBLAS (libscipy_openblas64_: scipy_dgesv_64_,
        scipy_dgetrs_64_, scipy_cblas_dgemv64_, scipy_cblas_ddot64_), 32
        from SciPy's (libscipy_openblas: scipy_dgetrf_, scipy_dgetrs_) ----- */
typedef void (*dgesv64_t)(const int64_t *, const int64_t *, double *, const int64_t *,
                          int64_t *, double *, const int64_t *, int64_t *);
typedef void (*dgetrs64_t)(const char *, const int64_t *, const int64_t *, const double *,
                           const int64_t *, const int64_t *, double *, const int64_t *,
                           int64_t *, size_t);
typedef void (*dgemv64_t)(int, int, int64_t, int64_t, double, const double *, int64_t,
                          const double *, int64_t, double, double *, int64_t);
/* (np.dot of two float64 vectors: numpy's DOUBLE_dot calls it) */
typedef double (*ddot64_t)(int64_t, const double *, int64_t, const double *, int64_t);
/* CBLAS's own values for dgemv64_t's first two arguments */
enum pcx_cblas { PCX_CBLAS_ROW_MAJOR = 101, PCX_CBLAS_NO_TRANS = 111, PCX_CBLAS_TRANS = 112 };
typedef void (*dgetrf32_t)(const int *, const int *, double *, const int *, int *, int *);
typedef void (*dgetrs32_t)(const char *, const int *, const int *, const double *, const int *,
                           const int *, double *, const int *, int *, size_t);

/* ---- KLU's int interface: klu_factor (a numeric, NULL where singular);
        klu_refactor, klu_solve (klu_tsolve, A^T x = b, its type too),
        klu_free_numeric (nonzero on success) -------------------------- */
typedef void *(*klu_factor_t)(int *, int *, double *, void *, void *);
typedef int (*klu_refactor_t)(int *, int *, double *, void *, void *, void *);
typedef int (*klu_solve_t)(void *, void *, int, int, double *, void *);
typedef int (*klu_free_numeric_t)(void **, void *);

/* ---- the evaluation core's codes (core.c), said by its callers too and
        read by `_tran_core` ------------------------------------------- */

/* the passes: a call's `want` (q is always made, so it has no want bit)
   and a fused kernel's pass mask (where q has its own) */
enum pcx_pass {
    PCX_PASS_G = 1,
    PCX_PASS_C = 2,
    PCX_PASS_I = 4,
    PCX_PASS_Q = 8,
    PCX_WANT_ALL = PCX_PASS_G | PCX_PASS_C | PCX_PASS_I,
    /* (not a pass: the caller's matrices are zero off the circuit's
       pattern, so the core writes only the pattern's entries -- the
       core's `npat`, `pat`) */
    PCX_WANT_PATTERN = 16
};
/* the companion the core forms: an integrator's, or none at all (the
   passes alone: G in binG, i in bini, C and q in the caller's buffers) */
enum pcx_formula {
    PCX_FORMULA_GEAR2 = 0,
    PCX_FORMULA_EULER = 1,
    PCX_FORMULA_TRAPEZOIDAL = 2,
    PCX_FORMULA_THETA = 3,
    PCX_FORMULA_PASSES = 4
};
