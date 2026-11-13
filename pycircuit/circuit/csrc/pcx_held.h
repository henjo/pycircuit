/* pcx_held.h -- THE HELD VALUES of StructuralKLU (`linearsolver`): the
 * values its run numeric (`_numeric`) was last REFACTORED on, so that a
 * refactor of the very same values is not made again -- KLU's refactor of
 * the same values on the same symbolic and pivot order is the same
 * numeric.  The C units that repeat `StructuralKLU.solve`'s sequence on
 * that numeric (newton.c's sparse full Newton, lte.c's sparse solve,
 * pcx_kluseq.h) keep StructuralKLU's own buffers, bound for the run:
 *
 *   hv  (float64, 2 nnz): the held values, then a checkpoint's backup;
 *   st  (int64, HELD_N): the state -- `linearsolver.HELD_*`, the same
 *       indices: OK (the held values are the numeric's last refactor),
 *       REUSED (refactors not made, StructuralKLU.reused), and a C call's
 *       checkpoint: OPEN, CK_OK, MOVED, CK_REUSED, CK_NUM.
 *
 * A NULL `st` is the rule off: every refactor made, nothing held.
 *
 * THE CHECKPOINT.  A C call the Python may repeat (a hand-back) opens one
 * where it starts (`pcx_held_open`; an enclosing call's open one is kept:
 * the C loop opens one an attempt) and closes it served
 * (`pcx_held_commit`) or handed back (`pcx_held_back`).  A hand-back leaves
 * the state, the count and the numeric as the checkpoint found them: where
 * the held values moved since and were valid then, the numeric recorded at
 * the open is refactored on the backup again (the same values on the same
 * pivots: the same numeric, whatever the call did to it in place) -- so the
 * Python's repeat decides, refactor by refactor, as a run without the C.
 * A numeric REPLACED since (a fallback's fresh factor: the one it found is
 * freed) cannot be restored: the hand-back then leaves no values held, and
 * the repeat refactors as without the rule -- its answer the same, one
 * refactor counted where a run without the C may reuse.
 * Included by a unit's .c before pcx_kluseq.h (a header includes none). */
#include <stdint.h>
#include <string.h>
#include <fenv.h>

enum pcx_held_index {
    PCX_HELD_OK = 0, PCX_HELD_REUSED, PCX_HELD_OPEN, PCX_HELD_CK_OK, PCX_HELD_MOVED,
    PCX_HELD_CK_REUSED, PCX_HELD_CK_NUM, PCX_HELD_KEPT, PCX_HELD_CK_KEPT, PCX_HELD_N
};

/* whether the numeric holds a refactor of exactly `Ax` (bit for bit: -0.0
   is not 0.0) -- counted where it does, no refactor to make */
static inline int pcx_held_same(int64_t *st, const double *hv, const double *Ax, long nnz)
{
    if (!st || !st[PCX_HELD_OK] || memcmp(hv, Ax, (size_t) nnz * sizeof(double)) != 0) return 0;
    st[PCX_HELD_REUSED]++;
    return 1;
}

/* the held state about to change: an open checkpoint keeps what it found */
static inline void pcx_held_move(int64_t *st, double *hv, long nnz)
{
    if (st[PCX_HELD_OPEN] && !st[PCX_HELD_MOVED]) {
        if (st[PCX_HELD_CK_OK]) memcpy(hv + nnz, hv, (size_t) nnz * sizeof(double));
        st[PCX_HELD_MOVED] = 1;
    }
}

/* a refactor of `Ax` made on the numeric (`ok`: KLU's verdict -- a failed
   refactor leaves no values held) */
static inline void pcx_held_refactored(int64_t *st, double *hv, const double *Ax, long nnz,
                                       int ok)
{
    if (!st) return;
    pcx_held_move(st, hv, nnz);
    if (ok) memcpy(hv, Ax, (size_t) nnz * sizeof(double));
    st[PCX_HELD_OK] = ok != 0;
}

/* the numeric made afresh, or replaced: no values held -- and an open
   checkpoint's numeric, which its replacement frees, no longer restored */
static inline void pcx_held_lost(int64_t *st, double *hv, long nnz)
{
    if (!st) return;
    pcx_held_move(st, hv, nnz);
    st[PCX_HELD_OK] = 0;
    if (st[PCX_HELD_OPEN]) st[PCX_HELD_CK_OK] = 0;
}

/* a checkpoint opened on the numeric `num` as it stands; 1 where this call
   owns it (none was open) */
static inline int pcx_held_open(int64_t *st, void *num)
{
    if (!st || st[PCX_HELD_OPEN]) return 0;
    st[PCX_HELD_OPEN] = 1;
    st[PCX_HELD_CK_OK] = st[PCX_HELD_OK] && num;
    st[PCX_HELD_MOVED] = 0;
    st[PCX_HELD_CK_REUSED] = st[PCX_HELD_REUSED];
    st[PCX_HELD_CK_KEPT] = st[PCX_HELD_KEPT];
    st[PCX_HELD_CK_NUM] = (int64_t) (intptr_t) num;
    return 1;
}

/* the owner's checkpoint closed on a served outcome: what the call did stands */
static inline void pcx_held_commit(int64_t *st, int own)
{
    if (st && own) st[PCX_HELD_OPEN] = 0;
}

/* the owner's checkpoint closed on a hand-back: the state, the count and
   the numeric as it found them (THE CHECKPOINT); the floating-point flags
   as they stood */
static inline void pcx_held_back(int64_t *st, double *hv, long nnz, int own, void *refactor,
                                 int *Ap, int *Ai, void *sym, void *common)
{
    if (!st || !own) return;
    st[PCX_HELD_REUSED] = st[PCX_HELD_CK_REUSED];
    st[PCX_HELD_KEPT] = st[PCX_HELD_CK_KEPT];
    if (st[PCX_HELD_MOVED]) {
        void *num = (void *) (intptr_t) st[PCX_HELD_CK_NUM];
        int ok = 0;
        if (st[PCX_HELD_CK_OK] && num) {
            fexcept_t flags;
            fegetexceptflag(&flags, FE_ALL_EXCEPT);
            ok = ((klu_refactor_t) refactor)(Ap, Ai, hv + nnz, sym, num, common);
            fesetexceptflag(&flags, FE_ALL_EXCEPT);
        }
        if (ok) memcpy(hv, hv + nnz, (size_t) nnz * sizeof(double));
        st[PCX_HELD_OK] = ok != 0;
    }
    st[PCX_HELD_OPEN] = 0;
    st[PCX_HELD_MOVED] = 0;
}
