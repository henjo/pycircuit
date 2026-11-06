/* core.h -- the evaluation core's plan and entry: the interface `_tran_core` loads
 * core.c through.  cffi lays the struct out from this text and
 * core.c includes it, so the compiler checks the entry against it.
 *
 * One `hdl_core_t` per circuit plan (`_tran_core._Core`): the circuit's
 * four passes -- m = 0 G, 1 C, 2 i, 3 q -- flattened into tables.
 * `_Core.__init__` sets every field once, from numpy arrays it keeps
 * (`_Core.arrays`; index arrays int64); `_Core.probe` fills the address
 * tables b_fn, PRall and z_fn before the calls it readies.  An n x n
 * matrix is row-major.  A call reads the tables and writes only its own
 * scratch (vals beyond its templates, binG, bini), the caller's output
 * buffers, each element's parameter row (its temperature slot, its held
 * statements) and, the bypass armed, its references and counts. */
typedef struct {
    long n;                       /* the circuit's size (unknowns) */
    /* THE VALUES every slot of every pass is made from (`_Flat.values`):
       G's and C's templates (`_MatrixPlan`: the constant stamps' entries,
       set once), the fused groups' staging (zst), the outputs of the
       batches outside a group, and the constant groups' products */
    double *vals;
    /* THE PULL (`np.bincount` of each pass's slots): pass m's bins written
       a call -- G's and C's the pattern's entries (every entry where npat
       is 0), i's and q's every row -- each from +0.0 with its slots'
       values added in slot order (a slot whose value is a constant zero
       left out: no addition of a zero changes a sum made from +0.0).  The
       bins are laid out in RUNS of one source count -- the count the same
       for every bin of a run, so no branch depends on a bin's -- run r of
       pass m being runall[2 r] sources a bin over runall[2 r + 1] bins;
       pass m's runs are runall[2 pr_off[m] ..), its bins wball[pl_off[m]
       ..) and its sources, bin after bin, ixall[px_off[m] ..) */
    const int64_t *pl_off;        /* (5) pass m's first written bin; the
                                     last the end */
    const int64_t *wball;         /* each written bin: its row-major entry
                                     of n x n (G, C) or its row (i, q) */
    const int64_t *pr_off;        /* (5) pass m's first run */
    const int64_t *runall;        /* (2 a run) its sources a bin, its bins */
    const int64_t *px_off;        /* (5) pass m's first source in ixall */
    const int64_t *ixall;         /* each source: its index in vals */
    /* THE BATCHES: one generated C kernel over a class's elements in one
       pass (`_hdl_batch.Batch`; a lone C-bound element a batch of one) */
    long nb;                      /* the batch count */
    const int64_t *b_pass;        /* (nb) its pass m */
    void **b_fn;                  /* (nb) its class's kernel (fused_fn_t),
                                     asked for the batch's pass alone;
                                     filled by `_Core.probe` */
    const int64_t *b_ne;          /* (nb) its element count ne */
    const int64_t *b_k;           /* (nb) inputs per element k (<= 64) */
    const int64_t *b_so;          /* (nb) outputs per element so (<= 4096) */
    const int64_t *b_ti;          /* (nb) the temperature's slot in a
                                     parameter row (-1: none) */
    const int64_t *NMall;         /* the node maps, ne x k each, of every
                                     batch and then every fused group:
                                     element e's input j is x[NM[e*k + j]] */
    const int64_t *b_nmoff;       /* (nb) its node map's offset in NMall */
    double **PRall;               /* each element's parameter row (its
                                     pack's array), per batch and then per
                                     fused group: filled by `_Core.probe` */
    const int64_t *b_proff;       /* (nb) its rows' offset in PRall */
    const int64_t *b_fz;          /* (nb) its fused group (-1: none) */
    const int64_t *orowall;       /* element e's outputs where the batch
                                     makes them itself: their offset in
                                     vals (its group's staging rows, or
                                     the batch's own) */
    const int64_t *b_ooff;        /* (nb) its elements' offset in orowall */
    /* THE CONSTANT GROUPS of the vector passes (i, q): the elements whose
       stamp is a constant k x k matrix S, one group per pass and k
       (`_stamp_plan._VectorPlan.groups`): y = S x[I] per element into the
       group's products in vals, read by the slots they fill */
    long ng;                      /* the group count */
    const int64_t *g_pass;        /* (ng) its pass m */
    const int64_t *g_E;           /* (ng) its element count E */
    const int64_t *g_k;           /* (ng) its width k (<= 64) */
    const double *Sall;           /* the stamps: E k x k row-major matrices
                                     a group */
    const int64_t *g_soff;        /* (ng) its stamps' offset in Sall */
    const int64_t *Iall;          /* the node maps: E x k a group */
    const int64_t *g_ioff;        /* (ng) its node map's offset in Iall */
    const int64_t *g_yoff;        /* (ng) its E x k products' offset in
                                     vals */
    void *dgemv;                  /* numpy's own `scipy_cblas_dgemv64_`
                                     (`_tran_core._numpy_dgemv`): np.dot's */
    /* THE CORE'S BINS, written by the calls that make them */
    double *binG;                 /* (n x n, zero-initialised) G's bins: G
                                     stays here (`_tran_core.passes` reads
                                     it) */
    double *bini;                 /* (n) i's bins */
    /* THE FUSED GROUPS: a class whose kernel makes every pass the call
       wants at once, before the passes, into the staging zst (in vals),
       which the passes' slots read.  None (nz 0): ztab and zst NULL */
    long nz;                      /* the group count */
    void **z_fn;                  /* (nz) its fused kernel (fused_fn_t); 0
                                     where it does not serve (`_Core.probe`):
                                     its batches call their own */
    const int64_t *ztab;          /* nz rows of ZW = 14 (core.c names the
                                     columns) */
    double *zst;                  /* the staging (in vals): per group and
                                     pass, ne rows of that pass's output
                                     length */
    /* THE PATTERN: every entry a matrix pass can write, and the diagonal */
    long npat;                    /* its entry count; 0: unused (above
                                     `PATTERN_MAX_FILL` of n x n, or switched
                                     off) -- every entry written */
    const int64_t *pat;           /* its entries, row-major indices of n x n,
                                     ascending */
    /* THE THREADS: the fused groups' kernels of a call made on this many
       threads -- the caller one of them -- where its elements reach the
       pool's minimum (core.c); 1 (or 0): the caller alone.  Every byte as
       one thread makes them */
    long nthreads;
    /* the pool those threads are (pool.c's `pcx_pool`, a pool_fn_t:
       `_tran_core.pool_address`), or NULL: the caller alone */
    void *pool_fn;
    /* THE SMALL PRODUCTS: bit k (1 << k) set, a constant group's k x k
       product (k 2 or 3) by OpenBLAS's dgemv's own operation order (core.c's
       `small_gemv`) -- set where `_tran_core._small_gemv_ok` found this
       process's dgemv to make exactly that; clear, the dgemv itself */
    long small_gemv;
    /* THE BYPASS (`_Core.bypass`; bp_at NULL: off): a fused group's element
       whose first input, and each other input's difference to it, lie
       within bp_rt * max(|now|, |then|) + bp_vn of those of its REFERENCE
       -- its last full evaluation (every pass its group has) -- at the
       reference's temperature is not evaluated: it is linearized there,
       i = Ir + G dX and q = Qr + C dX, G and C the reference's (SPICE's
       bypass, loaded as SPICE loads a bypassed device: its companion at the
       reference).  Per group z: its references' offset in bp_ref (-1: the
       group never bypasses), element e's at bp_at[z] + e (3 k + 2 k k): its
       k inputs, k i, k q, k x k G, k x k C.  Per element ordinal (its row
       in PRall): the reference's temperature, its state (0 none; 1 held,
       the staging's G and C its own; 2 held, the staging's G or C another
       call's) and the evaluations bypassed */
    const int64_t *bp_at;
    double *bp_ref;
    double *bp_T;
    int64_t *bp_st;
    int64_t *bp_n;
    double bp_rt, bp_vn;
    /* THE PACKED SETS (PCX_WANT_SETS): the matrix passes' bins in the
       pattern -- G's in `binGp`, C's in the caller's packed buffer -- each
       written bin's place in `pat`; the vector passes' rows as in wball.
       npat 0: unused (one entry each) */
    const int64_t *wpall;         /* (as wball) each written bin's index in
                                     pat (G, C) or its row (i, q) */
    double *binGp;                /* (npat) G's bins, packed: G stays here
                                     after a call with the sets */
} hdl_core_t;
/* THE EVALUATION: the passes `want` asks for (q always) at x, then the
   companion `formula` forms, f = (i + iq) + u and J = G + Geq -- bit for
   bit `_residual_and_jacobian` (want G | C | I), `jacobian_only` (G | C,
   or C alone where J is unread) and `residual_only` (I); with
   PCX_FORMULA_PASSES the passes alone (`_tran_core.passes`, the Radau
   units), returning after the bins.
   Assumes: x (n) finite (its callers check it; the core does not); the
   core ready (`_Core.probe` this call, or as stamped); q1, q2, iq1 and u
   (n) set where the formula and want read them; Cin (n x n) the held C
   where want lacks PCX_PASS_C (it may be C itself); with
   PCX_WANT_PATTERN in want and npat nonzero, C, Geq and J zero off the
   pattern -- only its entries are written; with PCX_WANT_PACKED as well, J
   (want G) is packed -- entry k is the pattern's k-th, `pat[k]`, npat values
   -- and Geq is not written (nor read: it may be any buffer); with
   PCX_WANT_SETS as well (THE PACKED SETS), C, Cin and Geq are packed too
   and G is made in binGp, binG untouched -- the same values at the same
   entries, each entry's operations the same.
   Writes: C (want has C), q (always), iq and Geq (unless the formula is
   PCX_FORMULA_PASSES; Geq not packed), F (want has I), J (want has G; npat
   values packed); binG (G) and bini
   (i), the core's own; each element's parameter row at its temperature
   slot; the bypass armed (bp_at), the references, states and counts of
   the elements evaluated.  F and J are touched only where written, so a caller may pass
   any buffer of the size otherwise (`_tran_core.evaluate` passes q and
   Geq).  Returns 0. */
long pcx_core(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);

/* what an element may be, for the core's stack buffers (core.c) --
   `_Core.__init__` refuses a wider one as Unservable: its inputs (a batch's
   `k`, a fused group's x), and with them a kernel's outputs (a batch's `so`:
   k x k at most) */
enum core_cap {
    CORE_MAX_INPUTS = 64,
    CORE_MAX_OUTPUTS = 4096,
    /* the fused groups a call makes on the pool (core.c), at most: more and
       the caller makes them alone -- `core_pool_job_t`'s rows */
    CORE_POOL_ZMAX = 64
};

/* THE POOL'S JOB (core.c): one call's fused groups made on the threads --
   the context of its range function (`fused_range`) on pool.c's pool */
typedef struct {
    const hdl_core_t *c;          /* the call's core */
    const double *x;              /* (n) its state */
    double T;                     /* its temperature */
    long nr;                      /* the groups it makes, at most
                                     CORE_POOL_ZMAX (the rows' width) */
    long total;                   /* their elements together */
    long zl[64];                  /* (nr) each group's index */
    long zw[64];                  /* (nr) the passes it makes */
    long off[65];                 /* (nr + 1) its first element in the
                                     job's numbering; the last the total */
} core_pool_job_t;

/* THE PULL ON THE POOL's job (core.c's `pull_pool`): a pass's bins split in
   blocks across the pool's threads -- the values the slots read (`V`), the
   bins, the pass's runs (CORE_PULL_RMAX at most) with their sources and
   written bins, and each run's first bin and first source */
enum core_pull { CORE_PULL_RMAX = 64 };
typedef struct {
    const double *V;                 /* the values the pass's slots read (`vals`) */
    double *bin;                     /* the pass's bins */
    const int64_t *run, *ix, *wb;    /* its runs (count, bins), sources, written bins */
    long nr, total;                  /* its runs and bins */
    long boff[CORE_PULL_RMAX + 1];   /* each run's first bin, in the pass's order */
    long xoff[CORE_PULL_RMAX + 1];   /* each run's first source index */
} core_pull_job_t;

/* a fused group's row in `ztab` (`_Core.__init__` writes it, core.c reads
   it): the group's elements, an element's inputs, the temperature's slot
   in a parameter row (-1: none), its passes (a PCX_PASS_* mask, q's bit
   included), its node map's offset in NMall, its rows' offset in PRall;
   then per pass -- G, C, i, q -- its output length per element (`_LEN`)
   and its staging's offset in zst (`_AT`); a pass m's pair is
   `ZT_G_LEN + 2 * m` and `ZT_G_AT + 2 * m` */
enum core_ztab {
    ZT_NE = 0, ZT_K = 1, ZT_TI = 2, ZT_PASSES = 3, ZT_NM = 4, ZT_PR = 5,
    ZT_G_LEN = 6, ZT_G_AT = 7, ZT_C_LEN = 8, ZT_C_AT = 9,
    ZT_I_LEN = 10, ZT_I_AT = 11, ZT_Q_LEN = 12, ZT_Q_AT = 13,
    ZT_WIDTH = 14
};
