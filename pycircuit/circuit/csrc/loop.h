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
    double early_th;            /* `Transient.early_abort` where the Python arms it
                                   (`_SteppingLoop.early_arm`), else 0.0: off */
    long base;                  /* the base integrator, PCX_FORMULA_*: GEAR2, EULER
                                   or TRAPEZOIDAL */
    double zs_ratio;            /* Gear-2's ZERO_STABILITY_RATIO */
    long shrink_guard;          /* Gear-2's shrink branch on */
    long keeps_run;             /* the error test keeps a running reference (relref) */
    long newton_reject;         /* the full Newton's iteration limit rejected here, as
                                   `_SteppingLoop.take` rejects its exception (adaptive,
                                   no chord, no damped retry or rescue armed); 0: the
                                   Newton's failure handed back */
    long cap;                   /* the outputs' rows */
    /* the state, in and out: `_tran_loop_c.run` writes it from the Python
       loop's before each call and reads it back after */
    double t, h, dt_last, dt_last2;
    long have_dt_last, have_dt_last2, landing, is_first, order_drop;
    long imposed;               /* `_SteppingLoop.imposed`: a failed Newton's smaller
                                   retry of a held step (the multistep family reads
                                   no hold; the glue writes it back) */
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
    long early_aborts;          /* attempts abandoned at their first iterate (each
                                   one of `rejected`, not of `attempts`) */
    long newton_rejects;        /* attempts rejected at their Newton's iteration limit
                                   (each one of `rejected`, not of `attempts`) */
    long early_checks;          /* the early abort's estimates the attempts made
                                   or abandoned made (a handed-back attempt's are
                                   the Python's: it makes the attempt again) */
    long have_hm;               /* the statistics' min/max step are set */
    double hmin, hmax;
    long status;
    /* THE WALK (`shooting/_lmmwalk_c`; `walk` 0: the stepping loop above):
       a shooting walk's multistep steps on a grid of its own -- each made by
       `take` with its end the grid's time and its sources the table's (the
       Python's own vectors), accepted as `accept` accepts it (no breakpoint,
       no judge), then the walk's step map on the monodromy's live columns
       (lmmsens.c, `_LiveMonodromy`).  The state above is the walk's inner
       transient's; X and T take each step's point */
    long walk;                  /* 1: the walk's steps (`wsteps` of the table) */
    long wsteps, wj;            /* the table's steps; in: the first to make, out:
                                   the next (`wsteps`: the table made) */
    const double *wh, *wt;      /* (wsteps) each step's length and end */
    const double *wu;           /* (wsteps x n) each step's sources at its end */
    void *map_fn, *map;         /* lmmsens.c, its lmmsens_t (the walk's) */
    long m, iref;               /* the reduced unknowns, the reference row */
    double *Jf;                 /* (m x m) the step's J, reduced */
    double *Cr[3];              /* (m x m) the reduced C ring: Cr[cn] C_n, Cr[cp]
                                   C_{n-1}, the third the step's own */
    long cn, cp;
    double *Pb[3];              /* (m x wcap) P_n (Pb[p0]), P_{n-1} (Pb[p1]) and
                                   the map's output, row-major on the live columns */
    double *Pqb[2];             /* (m x wcap) Pq_n (Pqb[pq]) and the map's */
    long p0, p1, pq;
    long w;                     /* the live columns (the call keeps them all) */
    long mono_n;                /* the steps the monodromy took (`_LiveMonodromy.n`) */
    long maps;                  /* out: the step maps made this call */
    long map_status;            /* out, LOOP_WALK_MAP: the map's status */
    /* THE STAGE FAMILY (`family` 1, `_tran_loop_c`'s stage binding):
       `_StageSteps` on a Radau IIA(3) step -- the dense coupled Newton
       `_rk_step_coupled` makes, radau.c on the transient's record.  An
       attempt is the stage walk's step on the loop's grid (pcx_stagestep.h:
       the seeds, each stage's screen; the nodes' merge at accept), the
       sources at the stage times, the end and the start (the table, each
       time once an attempt, in the order the Python asks), the embedded
       estimate (`_radau_error_estimate`), lte.c's STAGE judgement and the
       family's law (`_StageSteps.judge`).  The device memo the Python
       steps keep, two generations rolled at each attempt: an attempt writes
       the output set the last made attempt's is not in (`gen0`), so a stop
       after its Newton leaves both; the made attempt before that keeps only
       its stages (`ygen1`), from which the glue rebuilds the older
       generation by the core's passes at them (the very calls its Newton
       made).  The step's start is looked up by its bytes, as `_memo_get`
       keys it: in this attempt's set, then the last made attempt's (the
       first attempt: the entry's memo, `first_rec`) */
    long family;                /* enum loop_family */
    void *radau_fn, *radau;     /* radau.c, its radau_t (`_tran_radau_c`'s record) */
    double *oset[2][5];         /* the two output sets: Ya, qa, ia, Ca, Ga, each
                                   `maxa` slots (radau.h) */
    long gen0;                  /* the set the last made attempt wrote; -1: none
                                   (the entry's memo is the Python's) */
    long A0;                    /* its assemblies */
    double *ygen1;              /* (maxa x 3 x n) the made attempt before it: its
                                   stages, A1 assemblies (0: none since entry) */
    long A1;
    long made;                  /* the attempts made since entry */
    long first_rec;             /* 1: the start's record is in the glue's `xr_*`
                                   (the entry's: the memo the Python left holds it) */
    double ct[3], a33;          /* the tableau's c and A[2][2] */
    double gamma_r, dd1, dd2, dd3;  /* the estimate's constants (`GAMMA_REAL`, dd1..3) */
    double *ust;                /* (5 x n) work: the sources at the stage times, the
                                   end and the start */
    double *useed;              /* (3 x n) work: the stage seeds */
    double *xr_q, *xr_C, *xr_G, *xr_i;  /* (n, n x n, n x n, n) the start's passes
                                   where no record holds them (the core's) */
    long pass_bits, pass_zero;  /* those passes as `_tran_core.passes` asks them:
                                   the want bits (the pattern's where the core has
                                   one), C zeroed first where it writes the pattern */
    double *eF;                 /* (n) work: the estimate's F1, then its rhs */
    double *eA;                 /* (m x m) work: the real factor, column-major */
    double *eb;                 /* (m) work: the reduced rhs, then the solution */
    int64_t *eipiv;             /* (m) work: dgesv's pivots */
    void *dgesv, *dgemv;        /* numpy's own (dgesv64_t, dgemv64_t) */
    long under_ok;              /* numpy ignores underflow (`func._underflow_ignored`):
                                   the estimate's flag test does not count it */
    double safety, kclamp, pst; /* the law: SAFETY, K, p (EMBEDDED_ORDER + 1) */
    double min_shrink, rq_min;  /* MIN_SHRINK_RATIO, RETRY_ORDER_MIN */
    long retry_order;           /* enum loop_retry_order (`retry_order`) */
    long have_sretry, sretry_point;  /* the family's `_retry`: its point (an index
                                   into X, as the loop's ctrl keeps it), */
    double sretry_h, sretry_worst;   /* its step and worst */
    long have_sq;               /* `_q_carried` set */
    double sq;                  /* `_q_carried` */
    void *klu_cur;              /* the stage KLU's numeric (NULL: dense, or none yet) */
    double *hq;                 /* (n) the attempt's start charge, as its Newton reads it */
    double *liq;                /* (n) the last made attempt's `iq = -(i + u(t_end))` */
    double *lest;               /* (n) the last made attempt's estimate (adaptive) */
    long held;                  /* out: the attempt the loop stopped in converged its
                                   Newton -- its result for the glue's hold */
    double cur_h;               /* out: that attempt's step */
    /* out: the last attempt made's -- its set, stage times, end, step and
       estimate (in the error test's `lte`) -- and the family's counts */
    long att_set;
    double att_ts[3], att_tend;
    long st_fast, st_umiss, st_uhit, st_passes, st_estmemo, st_assemblies;
    long st_klu_refactors, st_klu_factors, st_klu_pivot, st_klu_resid, st_klu_on;
    /* THE TRANSFORM (`tmode` 1): the stage family's attempts by the cost
       transform, `_rk_step_transformed` -- x_n's C and G frozen for the
       step; the two factors (radau_fold.c; P from diag(lam) Tinv by numpy's
       complex division by h); the real one's LU in the loop's buffer
       (`NumpyLU`: dgesv at its first solve, then dgetrs); the complex one's
       values in the last record's pattern for the transient's
       ComplexKLUSolver, which radau_tc.c refactors at its first solve;
       radau_tc.c's Newton and its end passes on the loop's own radau_tf_t;
       each stage's screen on C by a core pass at the first two stages
       (`_C_at_state`'s `cir.C`, which the core makes bit for bit) and the
       end's at the last; the estimate on the frozen real solve.  The device
       memo as the Python's transform step keeps it, by state and passes in
       two generations (`tg`): the glue makes the records at write-back */
    long tmode;                 /* 1: the transform's attempts; 0: the coupled Newton's */
    void *tc_fn, *tc;           /* radau_tc.c, the loop's radau_tf_t */
    void *fold_fn, *fold;       /* radau_fold.c, its radau_fold_t */
    double tcm[12];             /* rows 0 and 1 of diag(lam) @ Tinv, numpy's (complex) */
    double lam0r, lam1r, lam1i; /* lam[0].real and lam[1]: the factors' scalars over h */
    double *lua;                /* (m x m, column-major) the real factor, then its LU */
    int64_t *lupiv;             /* (m) its pivots */
    double *crr, *grr;          /* (m x m) x_n's C and G reduced, the fold's inputs */
    double *tcC;                /* (2 x n x n) C at the first two stages (the screen's) */
    double *tcq;                /* (n) work: a pass's q not kept */
    double *tcY;                /* (3 x n) the last made attempt's stages */
    double *te[2][4];           /* the end passes' two sets: q, i (n), C, G (n x n); an
                                   attempt writes the one x_n's record is not in */
    long te_x;                  /* the set holding x_n's record; -1: `xr_*` holds what
                                   the C made or the entry's memo had (`xfields`) */
    long xfields;               /* the passes `xr_*` holds (LOOP_TG_*) */
    double *tg;                 /* (LOOP_TG_SLOTS x LOOP_TG_MAX x n) the memo's
                                   generations by state: an attempt's own in a slot
                                   the made attempts' two do not hold, a fallback's
                                   second in the fourth (loop.c's tg_*) */
    long tg_n[4], tg_f[4][32];  /* each slot's records and their passes (4:
                                   LOOP_TG_SLOTS, 32: LOOP_TG_MAX) */
    long tg_gen, tg_prev;       /* the slots of the last made attempt's generation
                                   and the one before (the attempt's own while made) */
    double tc_h;                /* the step `p` (P's rows) was made for; `tc_hp` 0: none */
    long tc_hp;
    double ph_made;             /* the last made attempt's step: `_radau_Pc`'s (`have_ph`) */
    long have_ph;
    long tc_maxit;              /* the Newton's `maxiter` */
    /* the early abort (`radau_early_abort`, the stepping loop's `radau_arm`):
       its threshold (0: never armed); f0 = -(i_n + u(t_n)) once a point
       (`_StageSteps.early_f0`): made for the point `f0_point` (`have_f0`) */
    double tc_th;
    double *tcf0;               /* (n) */
    long have_f0, f0_point;
    /* a new pattern of the complex factor (the fold's nonzeros off the last
       record's: `ComplexKLUSolver.prepare` gathers a new record and
       `solve_prepared` analyses and factors it afresh): the loop's two
       pattern buffers (Ap m + 1, Ai m x m, the mask m x m laid out as A.T),
       the one holding the current pattern (`pcur`, -1: the solver's own,
       the entry's) and KLU's entries for the analysis, the complex factor
       and the frees (`klu_analyze`, `klu_z_factor`, `klu_free_symbolic`,
       `klu_z_free_numeric`) */
    int *pAp[2], *pAi[2];
    uint8_t *pnz[2];
    long pcur;
    void *klu_analyze, *klu_zfactor, *klu_free_sym, *klu_free_znum;
    /* the counts the Python transform attempts' calls make, this call's --
       of every attempt its frozen factors reached (`tc_frozen`: made, or
       abandoned at the early abort): the Newton's iterations, its source
       memo hits and its complex refactors; its serves, its probes; the
       early abort's checks and aborts; the core passes at x_n (the
       context's q, the frozen C and G, the estimate's CGi, the check's i)
       and at the first two stages (`cir.C`); the end's records (made,
       handed, with a pass counted, unhanded); the screens and their memo
       hits; the estimates; the fold taken (`radau.frozen:fold`) */
    long tc_frozen, tc_iters, tc_uhit, tc_refactors, tc_served, tc_probes, tc_checks,
         tc_aborts, tc_xpasses, tc_cpasses, tc_ends, tc_handed, tc_endpass, tc_unhanded,
         tc_screened, tc_mcmemo, tc_ests, tc_folds;
    /* ... and of the attempts on a new pattern: the analyses and fresh
       factors, the folds whose values fell off the pattern, the C Newton's
       declines (`radau_tc:klu`: the Python loop then iterates), that loop's
       iterations (each stage's 'qi' pass and limiter walk), the end records
       its passes made (no hand-off) */
    long tc_analyses, tc_factors, tc_foldoff, tc_klu, tc_pyiters, tc_pyends;
    /* the analyses and factors a made attempt's new pattern freed (the
       solver's `analyses_freed`, `factors_freed`) */
    long tc_symfreed, tc_numfreed;
    /* THE FALLBACK (`_rk_step_coupled` after the transform's
       NoConvergenceError: the dense coupled Newton, radau.c on the
       transient's record and its stage KLU, the coupled family's sets
       above): `fb_ok` 1 where the glue bound it (else the loop stops there,
       the Python's); x_n's record in the coupled set the last made fallback
       wrote (`fbx`: its q, i, C, G at `fbv`); the last fallback's stages
       made this call (`fbY`, `have_fbY`); the transform Newton's
       non-convergences the Python counts (`radau_tc:noconv`), the fallbacks
       made, the certificates their screens made */
    long fb_ok, fbx;
    const double *fbv[4];
    double *fbY;
    long have_fbY;
    long tc_noconv, tc_fb, tc_fbcert;
    /* THE REAL FACTOR ON KLU (`rk_klu` 1: the analysis solver's
       StructuralKLU, radau_tc.h's `rk`, held by the loop -- `rk_hold`): its
       pattern's flat indices into the reduced m x m factor (`rk_take`,
       `rk_nnz`), where `tc_frozen` takes the fold's real factor's values
       once a step (a nonzero off them is the Python's: `_ours` hands it to
       the fallback solver); `rkx` (m) the solves' work.  An attempt starts
       the numeric's sequence from the one it finds; made, its sequence
       stands (`pcx_klu_finish`: the one it started from freed where another
       replaced it); not made, taken back.  The counts its made attempts'
       solves made (`tc_rkdense`: the residual's dense tests) */
    long rk_klu, rk_nnz;
    const int64_t *rk_take;
    double *rkx;
    long tc_rkref, tc_rkfac, tc_rkpiv, tc_rkres, tc_rkdense;
} loop_t;
long pcx_loop(loop_t *s);

/* one attempt's own values while loop.c makes it -- loop.c's working state,
   not the glue's: `where` sets its start and step, `take` its companion and
   buffer, `judge` its error and the next step */
typedef struct {
    double t, h, tn;            /* where: its start, its step, its end */
    int have_tn;                /* the walk's: `tn` the grid's, not t + h */
    double t_break;             /* where: the next breakpoint (`_next_breakpoint`) */
    int is_first, landing;      /* where: the order drop asked; an adaptive step
                                   truncated onto t_break, or a grid's landing */
    int clamped;                /* where: its end decided its size (a landing, the
                                   run's end): `_SteppingLoop.clamped` */
    double h_last;              /* take: the last accepted step (h on a run's first) */
    long formula;               /* take: the companion, PCX_FORMULA_* */
    double a0, a1, a2;          /* take: Gear-2's alphas (`bdf2_alphas`), else 0 */
    long nxt;                   /* take: the C and Geq buffer its Newton writes */
    long early_checks;          /* take: the early abort's estimates its Newton
                                   made -- counted where it is made or abandoned,
                                   not where it is handed back */
    double err, h_next;         /* judge: the error test's error (adaptive);
                                   the next step */
    /* the stage family's (stage_take): the set its Newton writes, its
       assemblies, the stage times, the start time (`t_end - h`), the step's
       nodes in merge order; the start's record (q, C, G, i: a set's slot or
       the glue's `xr_*`); the sources it asked (each time once: `ut`, the
       row `us` of `ust` holding it) and its counts */
    long set, A;
    double ts[3], t0, nodes[4];
    long ord[4];
    const double *rq, *rC, *rG, *ri;
    double ut[5];
    long us[5], nu, umiss, uhit, passes, estmemo, fast;
    double cmax;                /* the screen's running maximum after its stages */
    /* its outcome, either take's: the three converged stages (3 x n), and
       the last one's record -- q, i (n), C, G (n x n) */
    const double *Ys, *q3, *i3, *C3, *G3;
    /* the transform's (tc_take): the memo's slots before its roll; its
       Newton's iterations and refactors; the core passes the Python's calls
       would make at x_n and at the first two stages; the end's record
       (`fuse`: 1 handed, 2 handed with a pass counted, 3 unhanded); the
       screen's memo hits; the fold taken */
    long tg_gen0, tg_prev0, had_f0, f0_point0;
    long iters, refact, xpass, cpass, fuse, mcmemo, fold;
    /* ... and its outcome: armed (the early abort), its Newton's source memo
       hits, served, probed (not converged at the first iterate), checked,
       screened, estimated; the check's worst unknown.  `tg_gen0`, `tg_prev0`,
       `had_f0`, `f0_point0`: the memo's slots and the point's f0 before it,
       taken back where it is not made (`tc_undo`) */
    long armed, tuhit, served, probed, checked, screened, estimated;
    double worst;
    /* a new pattern (`tc_repattern`): set; the Python loop's iterations, its
       end record made by its passes; what it replaced, put back where the
       attempt is not made (`tc_undo`) */
    long repat, pyit, pyend, old_pcur, old_nnz, old_at;
    void *old_sym, *old_num, *old_Ap, *old_Ai, *old_nz;
    /* the fallback's (`tc_fallback`): made; the transform Newton's
       non-convergence counted as the C Newton's (`noconv`; the Python
       loop's own is not, `pyloop`); its seeds' count; the coupled set its
       Newton wrote; its screen's certificates */
    long fb, noconv, pyloop, fast2, fbset, fbcert;
} loop_attempt_t;

/* KLU's analysis (`klu_analyze(n, Ap, Ai, Common)`: a symbolic, NULL where
   it fails), a new pattern's in the transform's take; its factor, refactor,
   solve and free are pcx_abi.h's */
typedef void *(*klu_analyze_t)(int, int *, int *, void *);

/* the device memo's records a generation holds (`tg`): a transform
   attempt puts four at most -- x_n's, the first two stages' C, the end's;
   a fallback's Newton three a stage's assembly (radau.h's maxa: 24) */
enum loop_tg {
    LOOP_TG_MAX = 32,
    /* a record's passes */
    LOOP_TG_Q = 1, LOOP_TG_I = 2, LOOP_TG_C = 4, LOOP_TG_G = 8
};

/* the memo's slots (`tg`): the made attempts' two generations, an
   attempt's own, a fallback's second roll */
enum loop_tg_slots { LOOP_TG_SLOTS = 4 };

/* `family`: whose attempts the loop makes */
enum loop_family {
    LOOP_FAMILY_LMM = 0,        /* the multistep family (`_LMMSteps`) */
    LOOP_FAMILY_STAGE = 1       /* the stage family (`_StageSteps`, Radau IIA(3)) */
};

/* `retry_order`: the order a rejection's retry assumes (`_retry_order`) */
enum loop_retry_order {
    LOOP_RETRY_ASSUMED = 0,     /* the method's */
    LOOP_RETRY_MEASURED = 1,    /* from a point's second rejection, measured */
    LOOP_RETRY_CARRIED = 2      /* measured, and carried to later first rejections */
};

/* `judge`'s outcome inside loop.c, never a status: the attempt accepted, or
   rejected and retried at the point (the loop takes another attempt) */
enum loop_judged {
    LOOP_JUDGE_ACCEPT = 0,
    LOOP_JUDGE_RETRY = -1,
    LOOP_JUDGE_EARLY = -2,      /* take's: the attempt abandoned at its first iterate */
    LOOP_JUDGE_NEWTON = -3      /* take's: the Newton's iteration limit, rejected here */
};

/* pcx_loop's result, left in `status` too: why the loop stopped */
enum loop_status {
    LOOP_DONE = 1,              /* the run's end: `tend` reached */
    LOOP_FULL = 2,              /* the outputs are full: take them, come back */
    LOOP_BAIL_SOURCES = 3,      /* the table handed a time back */
    LOOP_BAIL_PRED = 4,         /* the predictor handed back (general path, a flag) */
    LOOP_BAIL_NEWTON = 5,       /* the Newton bailed: the Python takes the attempt (the
                                   full Newton's limit on an adaptive step above
                                   minstep is rejected here: LOOP_JUDGE_NEWTON) */
    LOOP_BAIL_SCREEN = 6,       /* the screen fired or handed back: the SVD is Python's */
    LOOP_BAIL_FORMULA = 7,      /* a companion this loop does not make */
    LOOP_BAIL_LTE = 8,          /* the error test handed back (or has no history) */
    LOOP_BAIL_FORCE = 9,        /* the budget spent: the force-accept is the Python's */
    LOOP_WALK_MAP = 10,         /* the walk: the step made, its map handed back
                                   (`map_status`: singular, not finite) -- the
                                   Python's `_LiveMonodromy.step` makes it */
    LOOP_WALK_PRUNE = 11,       /* the walk: the step and its map made, a column
                                   of the monodromy dead -- the prune is
                                   `_LiveMonodromy._prune`'s (its views) */
    /* the stage family's, before the attempt's Newton (nothing of it made) */
    LOOP_BAIL_NODES = 12,       /* the step's nodes not after the history's */
    LOOP_BAIL_RADAU = 13,       /* the C stage Newton handed the attempt back */
    /* after its converged Newton: the result held for the Python's repeat */
    LOOP_BAIL_ESTIMATE = 14,    /* the estimate: a flag numpy would warn of, or
                                   its real factor singular */
    /* the transform's (`tmode`): before its Newton, nothing of the attempt
       made; after it, nothing held -- the Python's repeat refactors the same
       values (deterministic) and the LU is the loop's own */
    LOOP_BAIL_FOLD = 15,        /* the factors: a flag, a non-finite input, or the
                                   complex one's nonzeros off the last pattern */
    LOOP_BAIL_TRANSFORM = 16,   /* radau_tc.c handed back, or did not converge where
                                   no fallback is bound (the Python's dense one) */
    LOOP_BAIL_END = 17          /* the end passes not made (a non-finite stage) */
};
