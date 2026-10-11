/* loop.h -- the stepping loop in C: its state and entry, the interface
 * `_tran_loop_c` loads loop.c through.  cffi lays the struct out from this text and loop.c includes it,
 * so the compiler checks the entry against it.
 *
 * The loop drives the units the Python steps call, through their entries
 * and their structs, as `_SteppingLoop` drives the Python pieces they stand
 * in for: the sources' table (where: the next breakpoint; take: u(t)), the
 * predictor (the seed), the C Newton (the step and its converged point),
 * the branch screen (after it) and, adaptive, the error test and the
 * controller's law (judge).  It holds what the Python loop carries between
 * attempts -- the charge and companion-current rings, the predictor's
 * nodes, the step record, the cached C, the error test's running
 * reference, the last attempt's point -- and stops before any attempt it
 * cannot make as the Python would, the state left as the Python loop leaves
 * it after the last attempt it made (`_tran_loop_c` writes it back).
 *
 * THE UNITS' STRUCTS ARE BORROWED.  For each attempt loop.c sets sixteen of
 * the Newton record's inputs: the companion's scalars (`formula`, `a0`-`a2`,
 * `h`, `theta`: `_tran_newton_c.solve` sets them anew on every call) and ten
 * that point into this struct's buffers (`q1`, `q2`, `iq1`, `u`, `xseed`,
 * `cx`, `cC`, `Cj`, `Geqj`, `have_cache`), which `_tran_loop_c.run` saves
 * before the call and restores after it, whatever the status.  It sets the
 * error test's, the predictor's, the screen's and the law's per-call inputs
 * the same way; their run constants and buffers are `_tran_loop_c._bind`'s.
 * The units, the sizes and the buffers below are `_bind`'s, once a run; the
 * run's constants and the state `_tran_loop_c.run` writes on every entry. */
typedef struct {
    /* the units: their entries, their structs */
    void *newton_fn, *newton;   /* newton.c, its newton_t (the C Newton's record) */
    void *pred_fn, *pred;       /* pred.c, its pred_t */
    void *branch_fn, *branch;   /* branch.c, its branch_t; NULL: nothing screened */
    void *src_fn, *src;         /* sources.c, its sources_t */
    void *lte_fn, *lte;         /* lte.c, its lte_t: the error test (adaptive) */
    void *ctrl_fn, *ctrl;       /* ctrl.c, its ctrl_t: the controller's law (adaptive) */
    long n, hist_len;           /* the unknowns (the full width); the rings' depth */
    /* the run (`_tran_loop_c.run`, every entry) */
    long adaptive;              /* 0: `fixed_timestep` */
    double tend, timestep, minbreak, max_step, minstep;
    long max_retries, max_reject;   /* the family's budget at one point */
    long base;                  /* the base integrator, PCX_FORMULA_*: GEAR2, EULER
                                   or TRAPEZOIDAL */
    double zs_ratio;            /* Gear-2's ZERO_STABILITY_RATIO */
    long shrink_guard;          /* Gear-2's shrink branch on */
    long keeps_run;             /* the error test keeps a running reference (relref) */
    long cap;                   /* the outputs' rows */
    /* the state, in and out: `_tran_loop_c.run` writes it from the Python
       loop's before each call and reads it back after */
    double t, h, dt_last, dt_last2;
    long have_dt_last, have_dt_last2, landing, is_first, order_drop;
    long rejects, point_retries;    /* this point's, as `judge` counts them */
    long npoint;                /* the accepted points (`X`'s length) */
    double *qring, *iqring;     /* hist_len x n each, newest first */
    double *xlast;              /* the last accepted point, full */
    long hk;                    /* the predictor's nodes, oldest first, <= 8 */
    double *ht, *hx;            /* their times (8), their states (8 x n) */
    double *cC[2], *Geq[2];     /* the last attempt's converged C and Geq,
                                   alternating: an attempt handed back after
                                   its Newton leaves the last one's in place */
    long cur;                   /* which holds the last attempt's; -1: no
                                   cached C yet */
    double cmax;                /* the branch screen's running maximum */
    double *run;                /* the error test's running reference (n) */
    long have_run;
    double *xatt, *qatt, *iqatt;    /* the last attempt's point, charge, current */
    double att_h, att_h_last, att_t, last_err;
    long att_formula, att_rejected;
    /* out: the steps taken, and the counts the Python's tallies are made
       from (zeroed at entry) */
    long k;                     /* accepted */
    double *X, *T;              /* k x n, k */
    double t_break;             /* the last attempt's breakpoint */
    long attempts, rejected, iters, euler_steps, breaks, screens, pred_fast, ltes, lte_euler;
    long have_hm;               /* the statistics' min/max step are set */
    double hmin, hmax;
    long status;
} loop_t;
long hdl_fn(loop_t *s);

/* hdl_fn's result, left in `status` too: why the loop stopped */
enum loop_status {
    LOOP_DONE = 1,              /* the run's end: `tend` reached */
    LOOP_FULL = 2,              /* the outputs are full: take them, come back */
    LOOP_BAIL_SOURCES = 3,      /* the table handed a time back */
    LOOP_BAIL_PRED = 4,         /* the predictor handed back (general path, a flag) */
    LOOP_BAIL_NEWTON = 5,       /* the Newton bailed: the Python takes the attempt */
    LOOP_BAIL_SCREEN = 6,       /* the screen fired or handed back: the SVD is Python's */
    LOOP_BAIL_FORMULA = 7,      /* a companion this loop does not make */
    LOOP_BAIL_LTE = 8,          /* the error test handed back (or has no history) */
    LOOP_BAIL_FORCE = 9         /* the budget spent: the force-accept is the Python's */
};
