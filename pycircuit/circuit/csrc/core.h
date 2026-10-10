/* core.h -- the evaluation core's plan and entry: the interface `_tran_core` loads
 * core.c through.  cffi lays the struct out from this text and
 * core.c includes it, so the compiler checks the entry against it. */
typedef struct {
    long n;
    const int64_t *p_L; const int64_t *p_toff; const double *tmplall;
    const int64_t *flatall; const int64_t *p_nbins;
    long nb; const int64_t *b_pass; void **b_fn; const int64_t *b_ne;
    const int64_t *b_k; const int64_t *b_so; const int64_t *b_ti;
    const int64_t *NMall; const int64_t *b_nmoff; double **PRall;
    const int64_t *b_proff; const int64_t *bdstall; const int64_t *b_dstoff;
    long ng; const int64_t *g_pass; const int64_t *g_E; const int64_t *g_k;
    const double *Sall; const int64_t *g_soff; const int64_t *Iall;
    const int64_t *g_ioff; const int64_t *g_nsel; const int64_t *selall;
    const int64_t *gdstall; const int64_t *g_seloff;
    void *dgemv;
    double *work; double *ywork; double *binG; double *bini;
    long nz; void **z_fn; const int64_t *ztab; double *zst;
    const int64_t *bftab; const int64_t *fidxall;
    long npat; const int64_t *pat;
} hdl_core_t;
long hdl_fn(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
