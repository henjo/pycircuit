/* sources.h -- the independent sources' waveforms in C: the interface
 * `_tran_sources` loads sources.c through.  cffi lays the struct out from
 * this text and sources.c includes it, so the compiler checks the entry
 * against it.
 *
 * WHO SETS WHAT.  `_tran_sources.Table` builds one struct per direct source
 * pass: `n`, `ns`, the per-source arrays and `out` (numpy arrays the table
 * keeps), once -- `Table.holds` re-checks what they were read from.  At
 * each call `Table._call`, or in the C loop loop.c (the breakpoint search
 * and each attempt's u), sets `mode` and `t`; hdl_fn writes `out` or
 * `t_next`, `status` and `bad`. */
typedef struct {
    long n,                 /* the unknowns: `out` holds n */
        ns,                 /* the sources: every per-source array holds ns */
        mode;               /* what this call makes: enum sources_mode */
    /* per source, in the pass's element order (int64) */
    const int64_t *kind,    /* the waveform: enum sources_kind */
        *isv,               /* 1: a voltage source (-v into row2); 0: a current
                               source (+v into row0, -v into row1) */
        *row0, *row1,       /* a current source's two node rows (0 otherwise) */
        *row2,              /* a voltage source's branch row (0 otherwise) */
        *poff,              /* its first parameter's place in `par` */
        *npts;              /* a PWL's points (0 otherwise) */
    const double *a,        /* ns: each source's DC value */
        *par;               /* the parameters, source after source, each kind's
                               in enum sources_kind's order */
    double t,               /* the time, each call */
        t_next;             /* OUT (SOURCES_NEXT): the earliest event strictly
                               after t, inf where none */
    double *out;            /* OUT (SOURCES_U), n: u(t) */
    long status,            /* OUT: hdl_fn's result, enum sources_status */
        bad;                /* OUT: the source a bail stopped at, -1 otherwise */
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
                                   non-finite, past the fast add's 1e300, past
                                   `_sin_scalar`'s guards or exp's range; in
                                   SOURCES_NEXT a pulse's period index past
                                   2^52 */
    SOURCES_BAIL_PWL = 3        /* PWL found no segment (a NaN time) */
};
