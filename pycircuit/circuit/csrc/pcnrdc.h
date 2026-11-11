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
    /* THE PATTERN (CTX; np 0: none): every entry of J that can be nonzero --
       the core's pattern (binG is zero off it), each Schur position and the
       junction pairs' entries.  With one J is PACKED: its entries on P in
       P's order, then one zero (`J[np]`: written by none), the n x n matrix
       zero off P -- so a pass over it is the dense pass.  J'S POSITIONS
       (`smap`, `dmap`, `jmap`, `kfull`, `ffull`) index that layout: packed
       with a pattern, row-major n x n without one */
    long np;
    const int64_t *P;     /* (np) row-major indices of n x n, ascending */
    /* KLU ON THE PATTERN (CALL; kfull NULL: the reduced R made dense -- with
       a pattern and KLU always set: `_pcnr_c`'s `on_pattern`) */
    const int64_t *kfull; /* (nact) each KLU entry's position in J (an entry
                             off P: J[np], the zero) */
    long nfor;            /* P's entries in the reduced matrix off KLU's pattern */
    const int64_t *ffull; /* (nfor) their positions in J: `_ours`' foreign test */
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
    void **kfn, **rfn;                    /* (nd) its class's kernel (i and G in one
                                             call) and its refine kernel */
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
           *J,            /* J_MNA, then the Schur matrix in place: packed
                             (np + 1) with a pattern, n x n without (THE
                             PATTERN) */
           *Jabs,         /* (n x n) |the Schur matrix|: with a pattern zero off it
                             (allocated zero, written on it alone) */
           *rc,           /* (n) the rhs corrections */
           *A,            /* (m x m, column-major) the dense solve's copy */
           *R,            /* (m x m, row-major) KLU's: the reduced Schur matrix;
                             on the pattern (kfull) zero between its dense
                             tests (allocated zero: PCX_A_ZERO) */
           *b,            /* (m) the solve's rhs, then the dense solve's solution */
           *kx, *ky, *kt, /* (m, m, 2 m) KLU's solution, the residual's A x, its
                             |A| |x| and its rows' entry counts */
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
    long dense_conv,      /* the convergence tests the dense dgemv decided */
         dense_resid;     /* the residual tests the dense one decided */
    long iters,           /* the iterations made: the Python's `it + 1` */
         status;          /* enum pcnrdc_status, as returned */
    /* J'S POSITIONS (CTX; THE PATTERN's layout) */
    const int64_t *smap;  /* each device's Schur positions in `schur`'s order:
                             per unknown its t (row, ra) entries, then its t
                             (row, rb) -- 2 t m a device, device by device */
    const int64_t *dmap;  /* (n) each row's diagonal */
    const int64_t *jmap;  /* (4 nj) each junction pair's (a, a), (a, b), (b, a),
                             (b, b) */
    /* THE SEED (CALL) */
    long seed;            /* 1: v_lim0 is not given -- the solve makes it from
                             x (`v_lim_init(limit=True)`) before its first
                             iteration; 0: v holds it */
    /* THE PACKED G (CTX): with the pattern, the core makes G in its packed
       bins (PCX_WANT_SETS in `pattern`: binGp, the same values at the same
       entries) and J is gathered from them */
    const int64_t *gmap;  /* (np) each entry of P's index in the core's
                             pattern, -1 off it (G zero there); set where P
                             is -- NULL (no pattern): J gathered from binG */
    /* THE HELD VALUES (pcx_held.h): StructuralKLU's `_held_v` and
       `_held_st` where the record is bound to its run numeric (NULL: the
       rule off); OWN `held_own`, whether this solve opened the checkpoint */
    double *held_v;
    int64_t *held_st;
    long held_own;
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
