/* libm.h -- the C library's functions over an array: the interface the libm
 * probe (`_hdl_cbackend._libm_probe`) loads libm.c through; libm.c includes
 * it, so the compiler checks the entry against it. */
long pcx_libm(long which, const double *a, const double *b, long n, double *out);

/* the functions, as `_hdl_cbackend.LIBM_PAIRS` names them (its C names) */
enum libm_fn {
    LIBM_EXP = 0,
    LIBM_LOG = 1,
    LIBM_LOG1P = 2,
    LIBM_EXPM1 = 3,
    LIBM_SIN = 4,
    LIBM_COS = 5,
    LIBM_TAN = 6,
    LIBM_ASIN = 7,
    LIBM_ACOS = 8,
    LIBM_ATAN = 9,
    LIBM_SINH = 10,
    LIBM_COSH = 11,
    LIBM_TANH = 12,
    LIBM_ATAN2 = 13,
    LIBM_HYPOT = 14,
    LIBM_POW = 15
};
