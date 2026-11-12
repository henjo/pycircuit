"""KLU'S ORDER WHERE THE STRUCTURAL PATTERN MISLEADS IT.

`linearsolver.StructuralKLU` factors on the circuit's structural pattern
(`analysis.jacobian_structure`): every place an element can stamp and the
whole diagonal.  Some of those places hold zero in every matrix the circuit
makes -- the diagonal of a voltage source's branch row, which only the DC's
ladders add to, and a constant stamp's zeros.  KLU's analysis (BTF's maximum
transversal, then AMD) takes them for entries, so it can choose an order
whose pivots are those zeros; the factor then pivots off the diagonal and
fills far past KLU's own estimate.  Measured on the IBM power grid ibmpg1t
(54 266 unknowns, 14 308 voltage sources, nearly all 0 V shorts): a
refactor 86 ms, against 6.6 ms on the same matrix's nonzeros.

So the analysis is TRIED.  Where the solver knows which of the pattern's
places an element can stamp (`analysis.stamped_places`, the stamp plan's G
and C), KLU's own symbolic is factored once on values made from that alone:
seeded off-diagonals in [-1, 1] at the stamped places, each stamped
diagonal its row's sum of magnitudes plus one (a circuit matrix's
dominance), zero at every other place.  Where KLU had to reallocate in that
factor, or filled more than `FILL` times its symbolic estimate, the order
is made again: a maximum matching on the stamped places (each column's
pivot a place an element stamps), KLU's own BTF and AMD on the pattern with
those columns, the composed order given to `klu_analyze_given` on the
pattern itself, and kept where its own trial fills less than KLU's order
did.  Else KLU's own analysis stands, untouched -- as where the circuit's
branch rows are under `BRANCH_SHARE` of its unknowns: no trial is made -- the same
symbolic, the same bytes.  Sized on 20 benchmark decks (41-519 unknowns),
chains and a mesh of 4 000 and both IBM grids: the trial's fill 0.89-1.22
of the estimate and no reallocation on all but the 54 266-unknown grid
(4.69, 20 reallocations), where the matched order's refactor is 9.8 ms.

KLU's structs are read at offsets the compiler gives on the installed
header (`suitesparse/klu.h`), as `pycircuit._testing.klumem` reads its one
field, and only once KLU's own behaviour confirms them; without the header,
a compiler or that confirmation, no trial is made and KLU's own order stands.

FROM A SIZE, NESTED DISSECTION.  From `ND_FROM` unknowns a third order is
tried: a nested dissection of A + A' (CHOLMOD's: METIS's node bisector
applied recursively, then CCOLAMD on the constraint sets -- `pcx_ndorder`,
CHOLMOD opened at its first use), given to KLU as P = Q, kept where its own
trial fills less than the order otherwise kept.  On a large 2-D mesh AMD
fills past a dissection: the IBM grid ibmpg2t (164 897 unknowns, few
always-zero places -- no reorder) 9.04 M entries and a 441 ms refactor on
KLU's order, 6.83 M and 231 ms on the dissection; ibmpg1t's matched order
(1.02 M) stays.  At that size the mask of stamped places is made whatever the
branch share (`analysis.stamped_places`).

`_paths` prefix `klu_order:`: `trial`, `reordered`, `kept_own` (the matched
order filled no less), `unread` (the offsets not had), `no_matching`,
`failed` (the trial's factor); `dissection` (tried), `dissected` (kept),
`no_cholmod`."""
import ctypes

import numpy as np

from pycircuit.circuit import _paths

#: the fill, over KLU's symbolic estimate of it, past which the trial's
#: factor reads as an order built on zeros
FILL = 2.0
#: the trial values' seed
SEED = 7919
#: the unknowns from which a nested dissection is tried (`dissected`): above
#: every circuit of the gate's size set -- the IBM grids
ND_FROM = 32768
#: the share of a circuit's unknowns its branch rows (a voltage source's,
#: an inductor's: a diagonal no element need stamp) must reach for the
#: trial to be made -- below it so few zero diagonals cannot mislead KLU
#: far, and no stamp plan is read (`analysis.stamped_places`)
BRANCH_SHARE = 0.01

_PC = _paths.COUNTS

_SRC = '''#include <stddef.h>
#include <stdint.h>
#include <suitesparse/klu.h>
long pcx_kluord(int k)
{
    switch (k) {
    case 0: return (long) offsetof(klu_common, nrealloc);
    case 1: return (long) offsetof(klu_symbolic, lnz);
    case 2: return (long) offsetof(klu_symbolic, unz);
    case 3: return (long) offsetof(klu_symbolic, P);
    case 4: return (long) offsetof(klu_symbolic, Q);
    case 5: return (long) offsetof(klu_numeric, lnz);
    case 6: return (long) offsetof(klu_numeric, unz);
    case 7: return (long) sizeof(((klu_numeric *) 0)->lnz);
    case 8: return (long) sizeof(((klu_symbolic *) 0)->P[0]);
    case 9: return (long) offsetof(klu_symbolic, n);
    case 10: return (long) sizeof(((klu_symbolic *) 0)->lnz);
    case 11: return (long) sizeof(((klu_common *) 0)->nrealloc);
    }
    return -1;
}
'''

_OFF = []

_ND_SRC = '''#include <stddef.h>
#include <stdint.h>
#include <dlfcn.h>
#include <suitesparse/cholmod.h>
typedef int (*c_start_t)(cholmod_common *);
typedef cholmod_sparse *(*c_alloc_t)(size_t, size_t, size_t, int, int, int, int, cholmod_common *);
typedef cholmod_sparse *(*c_trans_t)(cholmod_sparse *, int, cholmod_common *);
typedef cholmod_sparse *(*c_add_t)(cholmod_sparse *, cholmod_sparse *, double *, double *, int, int,
                                   cholmod_common *);
typedef int64_t (*c_nd_t)(cholmod_sparse *, int32_t *, size_t, int32_t *, int32_t *, int32_t *,
                          cholmod_common *);
typedef int (*c_free_t)(cholmod_sparse **, cholmod_common *);
/* a nested-dissection order of A + A' (CHOLMOD's) for the m x m CSC pattern
   (Ap, Ai) into Perm; CParent, Cmember (m) scratch.  1, or 0: CHOLMOD not
   had or failed */
long pcx_ndorder(long m, const int32_t *Ap, const int32_t *Ai, int32_t *Perm, int32_t *CParent,
                 int32_t *Cmember)
{
    cholmod_common cm;
    cholmod_sparse *A, *At = NULL, *B = NULL;
    long j, p, ok = 0;
    void *h = dlopen("libcholmod.so.5", RTLD_NOW | RTLD_LOCAL);
    c_start_t start, finish; c_alloc_t alloc; c_trans_t trans; c_add_t add; c_nd_t nd; c_free_t fr;
    if (!h) h = dlopen("libcholmod.so", RTLD_NOW | RTLD_LOCAL);
    if (!h) return 0;
    start = (c_start_t) dlsym(h, "cholmod_start"); finish = (c_start_t) dlsym(h, "cholmod_finish");
    alloc = (c_alloc_t) dlsym(h, "cholmod_allocate_sparse");
    trans = (c_trans_t) dlsym(h, "cholmod_transpose"); add = (c_add_t) dlsym(h, "cholmod_add");
    nd = (c_nd_t) dlsym(h, "cholmod_nested_dissection"); fr = (c_free_t) dlsym(h, "cholmod_free_sparse");
    if (!(start && finish && alloc && trans && add && nd && fr) || !start(&cm)) return 0;
    A = alloc((size_t) m, (size_t) m, (size_t) Ap[m], 1, 1, 0, CHOLMOD_PATTERN, &cm);
    if (A) {
        for (j = 0; j <= m; j++) ((int32_t *) A->p)[j] = Ap[j];
        for (p = 0; p < Ap[m]; p++) ((int32_t *) A->i)[p] = Ai[p];
        At = trans(A, 0, &cm);
        if (At) B = add(A, At, NULL, NULL, 0, 1, &cm);
        if (B) {
            B->stype = 1;
            ok = nd(B, NULL, 0, Perm, CParent, Cmember, &cm) >= 0;
        }
    }
    fr(&A, &cm); fr(&At, &cm); fr(&B, &cm);
    finish(&cm);
    return ok;
}
'''
_ND = []


def _nd_fn():
    """`(ffi, pcx_ndorder)` once built, or None (no header or compiler)."""
    if not _ND:
        _ND.append(None)
        try:
            from pycircuit.circuit import _cbuild
            ffi, fn, _key, _cold, _s = _cbuild.load_unit(
                _ND_SRC, 'long pcx_ndorder(long m, const int32_t *Ap, const int32_t *Ai, '
                         'int32_t *Perm, int32_t *CParent, int32_t *Cmember);', 'pcx_ndorder')
            _ND[0] = (ffi, fn)
        except Exception:       # noqa: BLE001 -- no header or compiler: no dissection tried
            return None
    return _ND[0]


def dissected(sk):
    """The pattern's symbolic on a nested dissection of A + A' (`pcx_ndorder`)
    given to KLU as P = Q -- or None (CHOLMOD not had, an analysis failed)."""
    got = _nd_fn()
    m = sk.m
    perm = np.empty(m, dtype=np.int32)
    if got is not None:
        ffi, fn = got
        cp, cmem = np.empty(m, dtype=np.int32), np.empty(m, dtype=np.int32)
        fb = ffi.from_buffer
        Ap = np.ascontiguousarray(sk.Ap, dtype=np.int32)
        Ai = np.ascontiguousarray(sk.Ai, dtype=np.int32)
        ok = fn(m, fb('int32_t *', Ap), fb('int32_t *', Ai), fb('int32_t *', perm),
                fb('int32_t *', cp), fb('int32_t *', cmem))
    if got is None or not ok or not np.array_equal(np.sort(perm), np.arange(m, dtype=np.int32)):
        _PC['klu_order:no_cholmod'] += 1
        return None
    lib, ct = sk._lib, sk._ct
    given = lib.klu_analyze_given
    given.argtypes = [ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                      ct.POINTER(ct.c_int), ct.c_void_p]
    given.restype = ct.c_void_p
    sym = given(m, sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(perm, ct.c_int),
                sk._ptr(perm, ct.c_int), ct.byref(sk._common))
    return sym or None


def _offsets():
    """`{name: offset or size}` for the fields read here, from the compiler
    on the installed header and confirmed (`_confirmed`) -- or None."""
    if _OFF:
        return _OFF[0]
    _OFF.append(None)
    try:
        from pycircuit.circuit import _cbuild
        _ffi, fn, _key, _cold, _s = _cbuild.load_unit(_SRC, 'long pcx_kluord(int k);', 'pcx_kluord')
        names = ('c_nrealloc', 's_lnz', 's_unz', 's_P', 's_Q', 'n_lnz', 'n_unz', 'n_lnz_size',
                 's_P_size', 's_n', 's_lnz_size', 'c_nrealloc_size')
        off = {nm: int(fn(k)) for k, nm in enumerate(names)}
    except Exception:       # noqa: BLE001 -- no header, compiler or KLU: no trial, not a failure
        return None
    if (off['n_lnz_size'] in (4, 8) and off['s_P_size'] == 4 and off['s_lnz_size'] == 8
            and off['c_nrealloc_size'] == 4 and _confirmed(off)):
        _OFF[0] = off
    return _OFF[0]


def _int(addr, size):
    return (ctypes.c_int32 if size == 4 else ctypes.c_int64).from_address(addr).value


def _perm(sym, off, key, n):
    p = ctypes.c_void_p.from_address(sym + off[key]).value
    return np.ctypeslib.as_array(ctypes.cast(p, ctypes.POINTER(ctypes.c_int32)), (n,)).copy()


def _confirmed(off):
    """KLU's own behaviour at the offsets: a 3 x 3 tridiagonal's symbolic
    has n 3, permutations for P and Q and a positive estimate; its factor
    no reallocation and at least n entries in L and U."""
    from pycircuit.circuit.linearsolver import KLUSolver
    try:
        sk = KLUSolver()
    except ImportError:
        return False
    lib, ct = sk._lib, sk._ct
    Ap = np.array([0, 2, 5, 7], dtype=np.int32)
    Ai = np.array([0, 1, 0, 1, 2, 1, 2], dtype=np.int32)
    Ax = np.array([4.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0])
    com = ct.create_string_buffer(4096)
    lib.klu_defaults(ct.byref(com))
    sym = lib.klu_analyze(3, sk._ptr(Ap, ct.c_int), sk._ptr(Ai, ct.c_int), ct.byref(com))
    if not sym:
        return False
    num = lib.klu_factor(sk._ptr(Ap, ct.c_int), sk._ptr(Ai, ct.c_int), sk._ptr(Ax, ct.c_double),
                         ct.c_void_p(sym), ct.byref(com))
    try:
        ok = (ctypes.c_int32.from_address(sym + off['s_n']).value == 3
              and sorted(_perm(sym, off, 's_P', 3)) == [0, 1, 2]
              and sorted(_perm(sym, off, 's_Q', 3)) == [0, 1, 2]
              and ctypes.c_double.from_address(sym + off['s_lnz']).value > 0.0
              and bool(num) and ctypes.c_int32.from_buffer(com, off['c_nrealloc']).value == 0
              and _int(num + off['n_lnz'], off['n_lnz_size']) >= 3
              and _int(num + off['n_unz'], off['n_lnz_size']) >= 3)
    finally:
        if num:
            lib.klu_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(com))
        lib.klu_free_symbolic(ct.byref(ct.c_void_p(sym)), ct.byref(com))
    return ok


def trial_values(m, rows, cols, stamped):
    """The trial's values on the pattern `(rows, cols)` (CSC order):
    seeded off-diagonals in [-1, 1] at the stamped places, each stamped
    diagonal its row's sum of magnitudes plus one, zero elsewhere."""
    v = np.where(stamped, np.random.default_rng(SEED).uniform(-1.0, 1.0, len(rows)), 0.0)
    dg = rows == cols
    off = ~dg
    rowsum = np.bincount(rows[off], weights=np.abs(v[off]), minlength=m)
    v[dg] = np.where(stamped[dg], rowsum[rows[dg]] + 1.0, 0.0)
    return v


def _trial(sk, sym, off, stamped):
    """The trial's factor of the symbolic `sym`: `(reallocations, entries
    in L and U, KLU's estimate of them)`, or None where the factor failed."""
    lib, ct = sk._lib, sk._ct
    rows = sk.Ai.astype(np.intp)
    cols = np.repeat(np.arange(sk.m), np.diff(sk.Ap))
    Ax = np.ascontiguousarray(trial_values(sk.m, rows, cols, stamped))
    com = ct.create_string_buffer(4096)
    lib.klu_defaults(ct.byref(com))
    num = lib.klu_factor(sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(Ax, ct.c_double),
                         ct.c_void_p(sym), ct.byref(com))
    if not num:
        return None
    try:
        realloc = ctypes.c_int32.from_buffer(com, off['c_nrealloc']).value
        fill = (_int(num + off['n_lnz'], off['n_lnz_size'])
                + _int(num + off['n_unz'], off['n_lnz_size']))
        est = (ctypes.c_double.from_address(sym + off['s_lnz']).value
               + ctypes.c_double.from_address(sym + off['s_unz']).value)
    finally:
        lib.klu_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(com))
    return realloc, fill, est


def _fails(sk, sym, off, stamped):
    """Whether the symbolic `sym` misleads KLU on the pattern: the trial's
    factor reallocated, or filled past `FILL` times the estimate; None
    where the factor failed."""
    got = _trial(sk, sym, off, stamped)
    if got is None:
        return None
    realloc, fill, est = got
    return realloc > 0 or fill > FILL * max(est, 1.0)


def _matched(sk, off, stamped):
    """The pattern's symbolic on the matched order: each column's pivot a
    stamped place (a maximum matching), KLU's BTF and AMD on the pattern
    with those columns, the composed order given on the pattern -- or None
    (no full matching, an analysis failed)."""
    import scipy.sparse as sp
    from scipy.sparse.csgraph import maximum_bipartite_matching
    lib, ct, m = sk._lib, sk._ct, sk.m
    rows = sk.Ai.astype(np.intp)
    cols = np.repeat(np.arange(m), np.diff(sk.Ap))
    PL = sp.csr_matrix((np.ones(int(stamped.sum())), (rows[stamped], cols[stamped])), shape=(m, m))
    match = maximum_bipartite_matching(PL, perm_type='column')
    if (match < 0).any():
        _PC['klu_order:no_matching'] += 1
        return None
    Q0 = np.ascontiguousarray(match, dtype=np.int32)
    S = sp.csc_matrix((np.ones(sk.nnz), sk.Ai, sk.Ap), shape=(m, m))[:, Q0].tocsc()
    S.sort_indices()
    Ap2 = np.ascontiguousarray(S.indptr, dtype=np.int32)
    Ai2 = np.ascontiguousarray(S.indices, dtype=np.int32)
    com = ct.create_string_buffer(4096)
    lib.klu_defaults(ct.byref(com))
    s1 = lib.klu_analyze(m, sk._ptr(Ap2, ct.c_int), sk._ptr(Ai2, ct.c_int), ct.byref(com))
    if not s1:
        return None
    try:
        P = np.ascontiguousarray(_perm(s1, off, 's_P', m))
        Q = np.ascontiguousarray(Q0[_perm(s1, off, 's_Q', m)], dtype=np.int32)
    finally:
        lib.klu_free_symbolic(ct.byref(ct.c_void_p(s1)), ct.byref(com))
    given = lib.klu_analyze_given
    given.argtypes = [ct.c_int, ct.POINTER(ct.c_int), ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                      ct.POINTER(ct.c_int), ct.c_void_p]
    given.restype = ct.c_void_p
    sym = given(m, sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int), sk._ptr(P, ct.c_int),
                sk._ptr(Q, ct.c_int), ct.byref(sk._common))
    return sym or None


def better(sk, sym, stamped):
    """The symbolic `sk` should keep in place of KLU's own `sym` -- the
    matched order's, where the trial says `sym` misleads KLU; from `ND_FROM`
    unknowns a nested dissection's, where it fills less than the order kept
    -- or None: keep `sym`.  `stamped`: per entry of `sk`'s pattern, whether
    an element stamps it."""
    if stamped is None or len(stamped) != sk.nnz:
        return None
    off = _offsets()
    if off is None:
        _PC['klu_order:unread'] += 1
        return None
    _PC['klu_order:trial'] += 1
    stamped = np.asarray(stamped, dtype=bool)
    own = _trial(sk, sym, off, stamped)
    if own is None:
        _PC['klu_order:failed'] += 1
        return None
    realloc, fill, est = own
    got = None
    if realloc > 0 or fill > FILL * max(est, 1.0):
        got = _matched(sk, off, stamped)
        ## (kept only where its own trial fills less than KLU's order did)
        mine = _trial(sk, got, off, stamped) if got is not None else None
        if got is not None and (mine is None or not mine[1] < fill):
            _free(sk, got)
            got = None
            _PC['klu_order:kept_own'] += 1
        elif got is not None:
            fill = mine[1]
            _PC['klu_order:reordered'] += 1
    if sk.m >= ND_FROM:
        _PC['klu_order:dissection'] += 1
        nd = dissected(sk)
        mine = _trial(sk, nd, off, stamped) if nd is not None else None
        if nd is not None and mine is not None and mine[1] < fill:
            if got is not None:
                _free(sk, got)
            got = nd
            _PC['klu_order:dissected'] += 1
        elif nd is not None:
            _free(sk, nd)
    return got


def _free(sk, sym):
    sk._lib.klu_free_symbolic(sk._ct.byref(sk._ct.c_void_p(sym)), sk._ct.byref(sk._common))
