/* pcnrdc.h -- PCNR's DC solve in C: its struct, entry and codes, the
 * interface `_pcnr_c` loads pcnrdc.c through (cffi lays the struct out from
 * this text, after pcx_klu.h's; pcnrdc.c includes both).
 *
 * One `pcnr.solve_dc`: the ordinary part on the circuit's VIEW (the view's
 * C core: `pcnr.pcnr_view`), each participant at its own limited unknowns
 * (its compiled `i` and `G` kernels at the local point its probe tree
 * builds, the affine remainder off, `block = G_nl T`), the Schur step, the
 * refine (each class's kernel: `_pcnr_c.refine_kernel`) and the three
 * convergence tests -- statement by statement the Python's.  The solve is
 * the DC's: numpy's dgesv, or StructuralKLU's own sequence (pcx_kluseq.h;
 * THE NUMERIC: pcx_klu.h); `choose` where the AutoSolver has not chosen
 * yet stops at the first Schur matrix (PCNRDC_CHOOSE).  Rows are the
 * circuit's (full n, the reference row included); `m = n - 1` the reduced
 * size the solve sees.  Per-device tables are concatenated, each device's
 * slice at its offset. */

typedef struct {
    /* THE VIEW'S CORE (CTX): `_tran_core._Core` on the view */
    const void *core;     /* its hdl_core_t */
    void *core_fn;        /* core.c's `pcx_core` (core_fn_t) */
    long pattern;         /* PCX_WANT_PATTERN where the core writes only its pattern */
    double T;             /* CALL: the temperature */
    /* THE SOLVE */
    long n,               /* CTX: unknowns (full) */
         iref,            /* CTX: the reference row */
         maxiter;         /* CALL */
    long m;               /* CALL: the reduced system's rows (n - 1 - held) */
    const int64_t *red;   /* CALL (n): each row's reduced index, -1 the
                             reference and each held row (`held`) */
    double reltol, abstol, iabstol;   /* CALL: `solve_dc`'s */
    double gshunt;        /* CALL: the deformation's conductance, 0 none */
    long nn;              /* CTX: the node rows (gshunt's: each but iref) */
    double gjgmin;        /* CALL: the junctions' added conductance, 0 none */
    long nj;              /* CTX: the junction pairs */
    const int64_t *ja, *jb;   /* CTX (nj): each pair's rows */
    const double *u;      /* CALL (n): the sources, scaled (`u * source`) */
    void *dgesv, *dgemv;  /* CTX: numpy's own LAPACK and BLAS */
    long choose;          /* CALL: 1 -- stop at the first Schur matrix */
    /* KLU (pcx_klu.h; klu.klu_sym NULL: the dense solve); the reduced Schur
       matrix's entries in its CSC order (`take`) */
    pcx_klu_t klu;
    long nact;
    const int64_t *act_r; /* (nact) the entry in the reduced, row-major matrix */
    /* THE PARTICIPANTS (CTX) */
    long nd,              /* devices */
         k;               /* limited unknowns, all devices */
    const int64_t *d_t, *d_m, *d_off;     /* (nd) rows, unknowns, offset in v */
    const int64_t *d_rows;                /* (nd) offset in `rows` */
    const int64_t *rows;                  /* each device's t circuit rows */
    const int64_t *prb;                   /* (2 k) each unknown's probe rows (ra, rb) */
    const int64_t *d_walk, *n_walk;       /* (nd) offset in `walk`, steps */
    const int64_t *walk;                  /* 4 a step: row, parent, unknown, sign (+-1) */
    const int64_t *d_mir, *n_mir;         /* (nd) offset in `mir`, pairs */
    const int64_t *mir;                   /* 2 a pair: row, the row it mirrors */
    const int64_t *d_tm;                  /* (nd) offset in `Tm` and `blk` (t m) */
    const double *Tm;                     /* each device's T (t x m, row-major) */
    const int64_t *d_lin;                 /* (nd) offset in `lin` (t t), -1 none */
    const double *lin;                    /* each device's affine remainder */
    void **ifn, **Gfn, **rfn;             /* (nd) its i, G and refine kernels */
    double **prow;                        /* (nd) its packed parameters */
    const int64_t *ti;                    /* (nd) the temperature's slot, -1 none */
    /* THE ITERATE (CALL in, OWN out): x0 and v_lim0 in; the solution (CONVERGED)
       or the last iterate (MAXITER) out */
    double *x;            /* (n) */
    double *v;            /* (k) */
    /* THE LAST ITERATION (OWN): what the glue's message reads on MAXITER */
    double *fe,           /* (n) f_eff */
           *dx,           /* (n) dx_MNA */
           *gl;           /* (k) g_lim */
    /* SCRATCH (CTX buffers) */
    double *F,            /* (n) g_MNA */
           *J,            /* (n x n) J_MNA, then the Schur matrix in place */
           *Jabs,         /* (n x n) |the Schur matrix| */
           *rc,           /* (n) the rhs corrections */
           *A,            /* (m x m, column-major) the dense solve's copy */
           *R,            /* (m x m, row-major) KLU's: the reduced Schur matrix */
           *b,            /* (m) the solve's rhs, then the dense solve's solution */
           *kx, *ky,      /* (m) KLU's solution, the residual's A x */
           *xn, *ax, *Is, /* (n) x_new, |x_new|, |J| |x_new| */
           *vn,           /* (k) v_new */
           *blk,          /* each device's block (t x m, at d_tm) */
           *pot, *cur,    /* (tmax) the local point, the device's current */
           *Gl,           /* (tmax x tmax) the device's G */
           *xs,           /* (tmax) the device's last accepted sub-vector */
           *rio,          /* (2 mmax) the refine's in/out */
           *Cs, *q;       /* (n x n, n) the core's C and q bins (q always made) */
    int64_t *ipiv64;      /* (m) dgesv's pivots */
    /* OUT */
    long iters,           /* the iterations made: the Python's `it + 1` */
         status;          /* enum pcnrdc_status, as returned */
} pcnrdc_t;

long pcx_pcnrdc(pcnrdc_t *s);

enum pcnrdc_status {
    PCNRDC_CONVERGED = 1,       /* x, v the solution, iters the count */
    PCNRDC_MAXITER = 2,         /* the iteration limit: fe, dx, gl the last
                                   iteration's (the glue raises) */
    PCNRDC_SINGULAR = 3,        /* dgesv's info, or KLU's failed factor: the
                                   Python's LinAlgError at iteration `iters + 1`
                                   (the glue raises) */
    PCNRDC_NONFINITE = 4,       /* x_new not finite: the Python's RuntimeError
                                   (`pcnr.NONFINITE`) at iteration `iters + 1`
                                   (the glue raises) */
    PCNRDC_CHOOSE = 5,          /* `choose`: the first Schur matrix in J (this
                                   and each code after it: not an outcome) */
    PCNRDC_BAIL_FLAGS = 6,      /* a floating-point flag numpy would report,
                                   or a value that is not finite where the
                                   C does not carry one: the Python's solve,
                                   its warnings */
    PCNRDC_BAIL_REFINE = 7,     /* a refine kernel declined (a NaN key) */
    PCNRDC_BAIL_FOREIGN = 8,    /* a nonzero off KLU's pattern: the Python's
                                   fallback solve */
    PCNRDC_BAIL_KLU = 9         /* klu_solve failed */
};
