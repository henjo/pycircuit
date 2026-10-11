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
 * scratch (work, ywork, binG, bini, zst), the caller's output buffers, and
 * each element's parameter row at its temperature slot. */
typedef struct {
    long n;                       /* the circuit's size (unknowns) */
    /* THE PASSES, each indexed by m: its slots (L of them) are a work
       buffer the constant stamps, the batches and the groups fill, then
       summed into its bins (`np.bincount`) */
    const int64_t *p_L;           /* (4) its slot count L */
    const int64_t *p_toff;        /* (4) its offset in tmplall and flatall */
    const double *tmplall;        /* G's and C's templates (`_MatrixPlan`:
                                     the constant stamps' entries
                                     pre-filled); i's and q's zeros, unread */
    const int64_t *flatall;       /* each slot's bin: its row-major entry
                                     of n x n (G, C) or its row (i, q) */
    const int64_t *p_nbins;       /* (4) its bin count: n x n or n */
    /* THE BATCHES: one generated C kernel over a class's elements in one
       pass (`_hdl_batch.Batch`; a lone C-bound element a batch of one) */
    long nb;                      /* the batch count */
    const int64_t *b_pass;        /* (nb) its pass m */
    void **b_fn;                  /* (nb) its kernel (pass_fn_t), filled by
                                     `_Core.probe` */
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
    const int64_t *bdstall;       /* output j of element e: its slot,
                                     dst[e*so + j] */
    const int64_t *b_dstoff;      /* (nb) its slots' offset in bdstall */
    /* THE CONSTANT GROUPS of the vector passes (i, q): the elements whose
       stamp is a constant k x k matrix S, one group per pass and k
       (`_stamp_plan._VectorPlan.groups`): y = S x[I] per element, then
       the selected entries of the y's into their slots */
    long ng;                      /* the group count */
    const int64_t *g_pass;        /* (ng) its pass m */
    const int64_t *g_E;           /* (ng) its element count E */
    const int64_t *g_k;           /* (ng) its width k (<= 64) */
    const double *Sall;           /* the stamps: E k x k row-major matrices
                                     a group */
    const int64_t *g_soff;        /* (ng) its stamps' offset in Sall */
    const int64_t *Iall;          /* the node maps: E x k a group */
    const int64_t *g_ioff;        /* (ng) its node map's offset in Iall */
    const int64_t *g_nsel;        /* (ng) its selected entries' count */
    const int64_t *selall;        /* each selected entry's index in the
                                     group's E x k products */
    const int64_t *gdstall;       /* its slot (paired with selall) */
    const int64_t *g_seloff;      /* (ng) its offset in selall and gdstall */
    void *dgemv;                  /* numpy's own `scipy_cblas_dgemv64_`
                                     (`_tran_core._numpy_dgemv`): np.dot's */
    /* THE CORE'S SCRATCH, written by the calls that make them */
    double *work;                 /* (max L) a pass's slots */
    double *ywork;                /* (max E*k) a group's products */
    double *binG;                 /* (n x n, zero-initialised) G's bins: G
                                     stays here (`_tran_core.passes` reads
                                     it) */
    double *bini;                 /* (n) i's bins */
    /* THE FUSED GROUPS: a class whose fused kernel makes every pass the
       call wants at once, before the passes, into the staging zst; its
       batches scatter from there.  None (nz 0): ztab, zst, bftab and
       fidxall NULL */
    long nz;                      /* the group count */
    void **z_fn;                  /* (nz) its fused kernel (fused_fn_t); 0
                                     where it does not serve (`_Core.probe`):
                                     its batches call their own */
    const int64_t *ztab;          /* nz rows of ZW = 14 (core.c names the
                                     columns) */
    double *zst;                  /* the staging: per group and pass, ne
                                     rows of that pass's output length */
    const int64_t *bftab;         /* (2 nb) per batch: its group (-1:
                                     none), its offset in fidxall (-1: the
                                     group's own element order) */
    const int64_t *fidxall;       /* per batch in another order: each of its
                                     elements' row in its group */
    /* THE PATTERN: every entry a matrix pass can write, and the diagonal */
    long npat;                    /* its entry count; 0: unused (above
                                     `PATTERN_MAX_FILL` of n x n, or switched
                                     off) -- every entry written */
    const int64_t *pat;           /* its entries, row-major indices of n x n,
                                     ascending */
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
   pattern -- only its entries are written.
   Writes: C (want has C), q (always), iq and Geq (unless the formula is
   PCX_FORMULA_PASSES), F (want has I), J (want has G); binG (G) and bini
   (i), the core's own; each element's parameter row at its temperature
   slot.  F and J are touched only where written, so a caller may pass
   any buffer of the size otherwise (`_tran_core.evaluate` passes q and
   Geq).  Returns 0. */
long hdl_fn(const hdl_core_t *c, const double *x, double T, long want, long formula,
            double a0, double a1, double a2, double h, double theta,
            const double *q1, const double *q2, const double *iq1, const double *u,
            const double *Cin, double *C, double *q, double *iq, double *Geq,
            double *F, double *J);
