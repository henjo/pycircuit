/* sources.h -- the independent sources' waveforms in C (the C stepping loop
 * plan's stage 1): the interface `_tran_sources` loads sources.c through.
 * cffi lays the struct out from this text and sources.c includes it, so the
 * compiler checks the entry against it. */
typedef struct {
    long n, ns, mode;
    const int64_t *kind, *isv, *row0, *row1, *row2, *poff, *npts;
    const double *a, *par;
    double t, t_next;
    double *out;
    long status, bad;
} sources_t;
long hdl_fn(sources_t *s);

/* a source's waveform (its `function`'s class), and its parameters' order
   in `par` from `poff` */
enum sources_kind {
    SOURCES_DC = 0,     /* the base `TimeFunction`: f = 0 */
    SOURCES_SIN = 1,    /* td, theta, omega, phase, offset, amplitude */
    SOURCES_PULSE = 2,  /* v1, v2, td, tr, tf, pw, per */
    SOURCES_PWL = 3,    /* npts times, then npts values */
    SOURCES_EXP = 4,    /* v1, v2, td1, tau1, td2, tau2 */
    SOURCES_AM = 5,     /* vo, va, fc, fm, m, pi */
    SOURCES_SFFM = 6    /* vo, va, fc, mdi, fm, pi */
};
/* what hdl_fn does: `mode` */
enum sources_mode {
    SOURCES_U = 0,      /* u(t) into out[n], each bin from +0.0 */
    SOURCES_NEXT = 1    /* the earliest event strictly after t, into t_next */
};
/* hdl_fn's result, left in `status` too (`_tran_sources.BAIL` names the
   bails from these; `bad` the source): */
enum sources_status {
    SOURCES_DONE = 1,
    SOURCES_BAIL_VALUE = 2,     /* a value the Python would make otherwise:
                                   non-finite, or past the fast add's 1e300 */
    SOURCES_BAIL_PWL = 3        /* PWL found no segment (a NaN time) */
};
