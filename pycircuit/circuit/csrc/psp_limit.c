
#include <math.h>
/* compact.PspMosLongChannel.limit: d, g, b about the source (rows 0, 1, 3
   about row 2), in Python's order; any non-finite value declines */
int hdl_fn(const double *x, const double *p, double *out)
{
    const double lim = %(lim)s;
    const double vs = out[2], vs0 = x[2];
    static const int ks[3] = {0, 1, 3};
    int j;
    (void) p;
    if (!(isfinite(vs) && isfinite(vs0))) return 1;
    for (j = 0; j < 3; j++) {
        const int k = ks[j];
        double vold, delta;
        if (!(isfinite(out[k]) && isfinite(x[k]))) return 1;
        vold = x[k] - vs0;
        delta = (out[k] - vs) - vold;
        if (!(isfinite(vold) && isfinite(delta))) return 1;
        if (fabs(delta) > lim) {
            const double nv = (vs + vold) + (delta > 0.0 ? lim : -lim);
            if (!isfinite(nv)) return 1;
            out[k] = nv;
        }
    }
    return 0;
}
