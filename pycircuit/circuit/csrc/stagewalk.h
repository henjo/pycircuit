/* stagewalk.h -- the shooting's plain stage walk in C: its state and entry,
 * the interface `_stagewalk_c` loads stagewalk.c through.  cffi lays the
 * struct out from this text and stagewalk.c includes it, so the compiler
 * checks the entry against it.
 *
 * The walk drives the units the Python step calls, through their entries
 * and their structs, as `_PeriodWalks._walk_stage` and the inner transient's
 * step drive the Python pieces they stand in for: the stage predictor (the
 * seeds), the C stage Newton (radau.c), the branch screen (branch.c, each
 * stage) and the walk's step map on its KLU (sensklu.c), the monodromy's
 * product `C_n P` numpy's own BLAS.  It holds what the Python walk carries
 * between steps -- the point and its charge and capacitance, the history
 * rings, the predictor's nodes, the screen's running maximum, the live
 * columns -- and stops before any step it cannot make as the Python would,
 * the state left as the Python walk leaves it after the last step it made
 * (`_stagewalk_c` writes it back).
 *
 * THE UNITS' STRUCTS ARE BORROWED: the radau record's (`_tran_radau_c`'s,
 * bound to the stage KLU), a predictor state of the glue's own (its `ht`,
 * `hx` this struct's), the screen's (`_tran_branch_c`'s binding), the stage
 * KLU's sensklu context (`_sensklu_c`'s).  stagewalk.c sets their per-step
 * fields; the glue sets the rest once a walk and restores what it must. */
typedef struct {
    /* THE UNITS (the glue, once a walk): their entries and structs */
    void *radau_fn, *radau;     /* radau.c's pcx_radau, its radau_t */
    void *pred_fn, *pred;       /* pred.c's pcx_pred, its pred_t */
    void *branch_fn, *branch;   /* branch.c's pcx_branch, its branch_t; NULL:
                                   nothing screened */
    void *sens_fn, *sens;       /* sensklu.c's pcx_sensklu, its sensklu_t */
    void *dgemm, *dgemv;        /* numpy's own scipy_cblas_dgemm64_ (dgemm64_t)
                                   and scipy_cblas_dgemv64_ (dgemv64_t) */
    /* SIZES (once a walk) */
    long n, m, iref;            /* the full unknowns, the reduced ones, the
                                   reference row */
    long steps;                 /* the walk's steps */
    long hist_len;              /* the history rings' depth */
    long live_on;               /* `_steps.LIVE_COLUMNS` and the walk's live
                                   columns not switched off */
    /* THE GRID AND THE SOURCES (once a walk, the glue's arrays) */
    const double *hs;           /* (steps) each step's length */
    const double *tend;         /* (steps) each step's end */
    const double *tst;          /* (steps x 3) each step's stage times */
    const double *U;            /* (steps x 4 x n) the sources at the three
                                   stage times and at the step's end */
    /* THE STATE, in and out: the glue writes it from the Python walk's
       before the call and reads it back after */
    long j;                     /* in: the first step; out: the next step to
                                   make (`steps`: the walk made) */
    double *x;                  /* (n) the point, full */
    double *xr;                 /* (m) the point, reduced */
    double *qn;                 /* (n) its charge */
    double *Cn;                 /* (m x m) its capacitance, reduced */
    long hk;                    /* the predictor's nodes, oldest first, <= 8 */
    double *ht, *hx;            /* their times (8), their states (8 x n) */
    double *qring, *iqring;     /* (hist_len x n each) the rings, newest first */
    double cmax;                /* the screen's running maximum */
    double *P;                  /* the monodromy's live columns: (m x k),
                                   column-major with leading dimension `ld`
                                   (`trans` 1), or row-major (`trans` 0) */
    long k, ld, trans;          /* its columns, leading dimension, layout */
    int64_t *cols;              /* (m) the full width's columns kept, where
                                   `have_cols` */
    long have_cols;             /* 0: every column (`_LiveColumns.cols` None) */
    long nadv;                  /* the steps advanced (`_LiveColumns.n`) */
    void *walk_num;             /* the walk's numeric (`_KLUWalk._num`) */
    /* SCRATCH (the glue's, once a walk) */
    double *seed;               /* (3 x n) the stage seeds */
    double *base;               /* (m x m) `C_n P`, row-major (k columns used) */
    double *Xs, *Xs2;           /* (3 m x m each) the step map's solution,
                                   column-major (k columns used): the one
                                   the carried P is not in */
    double *Pc;                 /* (m x m) the live columns compacted,
                                   column-major */
    /* THE RECORDS (out): each step made's point, reduced, and its stages */
    double *Xout;               /* (steps x m) */
    double *Yout;               /* (steps x 3 x n) */
    /* COUNTS (out, this call's): the stage KLU's and the walk numeric's, as
       the Python's calls add them */
    long klu_refactors, klu_factors, klu_pivot_fallbacks, klu_residual_fallbacks;
    long walk_refactors, walk_fallbacks, walk_factors;
    long assemblies;            /* the C stage Newton's, summed */
    long last_assemblies;       /* the last step made's (its final block) */
    void *klu_out;              /* the stage KLU's numeric after the last step
                                   made */
    long status;                /* a stagewalk_status, also pcx_stagewalk's
                                   result */
} stagewalk_t;
long pcx_stagewalk(stagewalk_t *t);

/* pcx_stagewalk's result, left in `status` too (`_stagewalk_c.BAIL` names the
   bails from these).  A bail stops before step `j`, nothing of it
   committed: the units' own hand-backs leave their state as they found it
   (a refactor in place keeps its pivot order, so the Python's repeat of the
   step makes the same numerics); the Python walk goes on from `j`: */
enum stagewalk_status {
    STAGEWALK_DONE = 1,         /* every step made: `j` is `steps` */
    /* before the step's Newton */
    STAGEWALK_BAIL_LIVE = 2,    /* no live column */
    STAGEWALK_BAIL_NODES = 3,   /* the step's nodes not after the history's */
    STAGEWALK_BAIL_PRED = 4,    /* a seed the predictor's general path makes */
    STAGEWALK_BAIL_RADAU = 5,   /* the C stage Newton handed the step back */
    /* after it converged: its result left in the radau record (the stage
       KLU's numeric moved with it), which `_stagewalk_c` hands to
       `_tran_radau_c.solve` for the Python's repeat of the step */
    STAGEWALK_BAIL_FINAL = 6,   /* no final block (not finite) */
    STAGEWALK_BAIL_REF = 7,     /* the point's reference entry not +0.0 */
    STAGEWALK_BAIL_SCREEN = 8,  /* a stage the screen did not find QUIET */
    STAGEWALK_BAIL_SENS = 9,    /* the step map handed back (a check flag) */
    /* step `j` made, its step map's numeric and counts with it, `j` past
       it: what the Python does next is the glue's */
    STAGEWALK_NONFINITE = 10,   /* the step map not finite on dropped columns:
                                   the live columns not advanced
                                   (`_LiveColumns.advance`'s full width) */
    STAGEWALK_SINGULAR = 11,    /* the step map's factor singular: the
                                   Python's LinAlgError */
    STAGEWALK_SOLVE_FAILED = 12 /* the step map's solve failed: the same */
};
