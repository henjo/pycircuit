
typedef void (*hdl_batch_fn_t)(const double *, const double *, double *);
void hdl_fn(const void *fn, const double *X, double **PR, double *OUT,
            long n, long sx, long so, long ti, double T)
{
    hdl_batch_fn_t f = (hdl_batch_fn_t) fn;
    long e;
    if (ti >= 0) {
        for (e = 0; e < n; e++) {
            PR[e][ti] = T;
            f(X + e * sx, PR[e], OUT + e * so);
        }
    } else {
        for (e = 0; e < n; e++)
            f(X + e * sx, PR[e], OUT + e * so);
    }
}
