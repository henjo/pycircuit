/* fixed.h -- the fixed-step multistep loop in C (the C stepping loop plan's
 * stage 4): its state and entry, the interface `_tran_fixed_c` loads fixed.c
 * through.  cffi lays the struct out from this text and fixed.c includes
 * it, so the compiler checks the entry against it.
 *
 * The loop drives four units through their entries and their structs, as
 * `_SteppingLoop` drives the Python pieces they stand in for: the sources'
 * table (where: the next breakpoint; take: u(t)), the predictor (the seed),
 * the C Newton (the step and its converged point) and the branch screen
 * (after it).  It holds what the Python loop carries between steps -- the
 * charge and companion-current rings, the predictor's nodes, the step
 * record, the cached C -- and stops before any step it cannot take as the
 * Python would, the state left as the Python loop leaves it after the last
 * step it took (`_tran_fixed_c` writes it back). */
typedef struct {
    /* the units: their entries, their structs */
    void *newton_fn, *newton;   /* newton.c, its newton_t (the C Newton's record) */
    void *pred_fn, *pred;       /* pred.c, its pred_t */
    void *branch_fn, *branch;   /* branch.c, its branch_t; NULL: nothing screened */
    void *src_fn, *src;         /* sources.c, its sources_t */
    long n, hist_len;
    /* the run */
    double tend, timestep, minbreak;
    long base;                  /* the base integrator, PCX_FORMULA_*: GEAR2, EULER
                                   or TRAPEZOIDAL */
    double zs_ratio;            /* Gear-2's ZERO_STABILITY_RATIO */
    long shrink_guard;          /* Gear-2's shrink branch on */
    long cap;                   /* the outputs' rows */
    /* the state, in and out */
    double t, h, dt_last, dt_last2;
    long have_dt_last, have_dt_last2, landing, is_first;
    double *qring, *iqring;     /* hist_len x n each, newest first */
    double *xlast;              /* the last accepted point, full */
    long hk;                    /* the predictor's nodes, oldest first, <= 8 */
    double *ht, *hx;            /* their times (8), their states (8 x n) */
    double *cC[2], *Geq[2];     /* the converged point's C and Geq, alternating:
                                   a step handed back after its Newton leaves
                                   the last accepted step's in place */
    long cur;                   /* which holds the last accepted step's; -1:
                                   no cached C yet */
    double cmax;                /* the branch screen's running maximum */
    /* out: the steps taken */
    long k;
    double *X, *T;              /* k x n, k */
    double t_break;             /* where's breakpoint for the step it stopped at */
    long last_formula;          /* the last step's companion (PCX_FORMULA_*) */
    double last_h_last;         /* and the h_last it was made with */
    long iters, euler_steps, breaks, screens, pred_fast;
    long have_hm;               /* the statistics' min/max step are set */
    double hmin, hmax;
    long status;
} fixed_t;
long hdl_fn(fixed_t *s);

/* hdl_fn's result, left in `status` too: why the loop stopped */
enum fixed_status {
    FIXED_DONE = 1,             /* the run's end: `tend` reached */
    FIXED_FULL = 2,             /* the outputs are full: take them, come back */
    FIXED_BAIL_SOURCES = 3,     /* the table handed a time back */
    FIXED_BAIL_PRED = 4,        /* the predictor handed back (general path, a flag) */
    FIXED_BAIL_NEWTON = 5,      /* the Newton bailed: the Python takes the step */
    FIXED_BAIL_SCREEN = 6,      /* the screen fired or handed back: the SVD is Python's */
    FIXED_BAIL_FORMULA = 7      /* a companion this loop does not make */
};
