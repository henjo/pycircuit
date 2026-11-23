"""Branch detection: did the step equation have more than one root.  A theme of
`Transient` (see `transient.py`).
"""

import contextlib
import math

import numpy as np

from pycircuit.circuit import (
    _evalhint,
    _packed,
    _paths,
    _stamp_plan,
    _tran_branch_c,
    _tran_core,
)
from pycircuit.circuit._limiting import state_restore, state_snapshot
from pycircuit.circuit.analysis import insert_row
from pycircuit.circuit.simwarnings import (
    ModelWarning,
    warn,
)
from pycircuit.utilities import param as _param

#: the screen's verdicts where no direction comes with them: quiet, and C
#: collapsed entirely (`_branch_screen`)
_QUIET, _COLLAPSED = (False, None), (True, None)
#: the doubles' unit roundoff, as the C's DBL_EPSILON (`_certify`)
_EPS = float(np.finfo(float).eps)


def _components(n, r, c):
    """`scipy.sparse.csgraph.connected_components`' answer, ``directed=True,
    connection='weak'``, for the graph on `n` nodes with the edges ``(r[k],
    c[k])``: ``(count, labels)``, each component labelled by its rank among
    the components ordered by their smallest node, `int32` -- without scipy,
    whose import was 0.12 s of every run's start-up.  Each node's parent is
    hooked to the smaller of an edge's two roots and pointers jumped to the
    roots until no edge joins two: a component's root is its smallest node,
    which no hook can move."""
    parent = np.arange(n, dtype=np.intp)
    r = np.asarray(r, dtype=np.intp)
    c = np.asarray(c, dtype=np.intp)
    while r.size:
        pr, pc = parent[r], parent[c]
        apart = pr != pc
        if not apart.any():
            break
        np.minimum.at(parent, np.maximum(pr, pc)[apart], np.minimum(pr, pc)[apart])
        while True:
            up = parent[parent]
            if np.array_equal(up, parent):
                break
            parent = up
    roots, labels = np.unique(parent, return_inverse=True)
    return roots.size, labels.astype(np.int32)


#: the block size from which `_rank_by_blocks` tries THE FULL-RANK PROOF
#: before a block's SVD (`_full_rank_bound`)
PROOF_MIN = 64

#: libklu through ctypes (`_klu`): `KLUSolver`'s library, bound without its
#: scipy -- a run imports none
_KLU = {}


def _klu():
    """`(ctypes, lib, common)` for libklu's analysis, factor, solves and
    frees, or None where it does not load: bound once a process."""
    got = _KLU.get('lib')
    if got is None:
        import ctypes
        got = False
        try:
            lib = ctypes.CDLL('libklu.so.2')
        except OSError:
            lib = None
        if lib is not None:
            ip, dp, vp, ci = (ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
                              ctypes.c_void_p, ctypes.c_int)
            for nm, argt, res in (('klu_defaults', [vp], ci),
                                  ('klu_analyze', [ci, ip, ip, vp], vp),
                                  ('klu_factor', [ip, ip, dp, vp, vp], vp),
                                  ('klu_solve', [vp, vp, ci, ci, dp, vp], ci),
                                  ('klu_tsolve', [vp, vp, ci, ci, dp, vp], ci),
                                  ('klu_free_numeric', [vp, vp], ci),
                                  ('klu_free_symbolic', [vp, vp], ci)):
                f = getattr(lib, nm)
                f.argtypes, f.restype = argt, res
            common = ctypes.create_string_buffer(4096)
            if lib.klu_defaults(ctypes.byref(common)):
                got = (ctypes, lib, common)
        _KLU['lib'] = got
    return got or None


def _m_solves(k, r, c, a):
    """`(u, v)`: M u = e and M^T v = e by KLU for the k x k matrix with the
    entries `a` at `(r, c)` (no repeats), or None where KLU does not load or
    the matrix does not factor."""
    klu = _klu()
    if klu is None:
        return None
    ct, lib, common = klu
    o = np.lexsort((r, c))
    Ai = np.ascontiguousarray(r[o], dtype=np.int32)
    Ax = np.ascontiguousarray(a[o], dtype=float)
    Ap = np.zeros(k + 1, dtype=np.int32)
    np.cumsum(np.bincount(c, minlength=k), out=Ap[1:])
    ip, dp = ct.POINTER(ct.c_int), ct.POINTER(ct.c_double)
    sym = lib.klu_analyze(k, Ap.ctypes.data_as(ip), Ai.ctypes.data_as(ip), ct.byref(common))
    if not sym:
        return None
    num = None
    try:
        num = lib.klu_factor(Ap.ctypes.data_as(ip), Ai.ctypes.data_as(ip), Ax.ctypes.data_as(dp),
                             sym, ct.byref(common))
        if not num:
            return None
        u, v = np.ones(k), np.ones(k)
        if not (lib.klu_solve(sym, num, k, 1, u.ctypes.data_as(dp), ct.byref(common))
                and lib.klu_tsolve(sym, num, k, 1, v.ctypes.data_as(dp), ct.byref(common))):
            return None
        return u, v
    finally:
        if num:
            lib.klu_free_numeric(ct.byref(ct.c_void_p(num)), ct.byref(common))
        lib.klu_free_symbolic(ct.byref(ct.c_void_p(sym)), ct.byref(common))


def _inverse_norm_bound(k, r, c, a, x):
    """An upper bound on ||M^-1||_inf for the Z-matrix M (entries `a` at
    `(r, c)`, off the diagonal <= 0) from `x` >= 0 with M x > 0 -- which
    makes M a nonsingular M-matrix, M^-1 >= 0, so M^-1 e <= x / min(M x):
    `max x / min M x`, M x bounded below by its sums' rounding (each row's
    `m` products and sums within (m + 1) eps |M| x).  inf where it fails."""
    if not (x.size == k and np.all(np.isfinite(x)) and np.all(x >= 0.0)):
        return math.inf
    m = int(np.bincount(r, minlength=k).max())
    mx = np.bincount(r, weights=a * x[c], minlength=k)
    ax = np.bincount(r, weights=np.abs(a) * x[c], minlength=k)
    low = float(np.min(mx - (m + 2) * 1.01 * _EPS * ax))
    if not low > 0.0:
        return math.inf
    return float(np.max(x)) / low * (1.0 + 4.0 * _EPS)


def _full_rank_bound(k, r, c, b):
    """THE FULL-RANK PROOF of a block (`_rank_by_blocks`): a lower bound on
    the smallest singular value of the k x k matrix B (entries `b` at
    `(r, c)`, row-major, no repeats), or 0.0 where none is shown.

    B an H-matrix: its comparison matrix M (|b_ii| on the diagonal, -|b_ij|
    off it) a nonsingular M-matrix -- shown by u, v >= 0 with M u > 0 and
    M^T v > 0, u and v KLU's solves of M u = e and M^T v = e (their accuracy
    only sharpens the bound) -- gives |B^-1| <= M^-1 (Ostrowski), so
    ||B^-1||_inf <= ||M^-1||_inf, ||B^-1||_1 <= ||M^-1||_1, and
    sigma_min(B) = 1 / ||B^-1||_2 >= 1 / sqrt(||B^-1||_1 ||B^-1||_inf).
    A capacitance matrix whose every piece reaches ground through a
    capacitor is one; a block missing its diagonal, or coupled more
    strongly than it is grounded, is not."""
    diag = r == c
    if int(np.count_nonzero(diag)) != k:
        return 0.0
    a = np.where(diag, np.abs(b), -np.abs(b))
    got = _m_solves(k, r, c, a)
    if got is None:
        return 0.0
    n_inf = _inverse_norm_bound(k, r, c, a, got[0])
    n_one = _inverse_norm_bound(k, c, r, a, got[1])
    if not math.isfinite(n_inf * n_one):
        return 0.0
    return (1.0 - 4.0 * _EPS) / math.sqrt(n_inf * n_one)


def _unproved(C, rows, cols, label, members, keep, tol, band, count=True):
    """The rows of `members` (blocks of one size k, `_rank_by_blocks`) whose
    block THE FULL-RANK PROOF does not settle: proved, a block's bound
    exceeds the band's top (`tol` without a band) by 8 k eps ||B||_F, its
    SVD's rounding as the band reckons it -- every singular value its SVD
    finds is then above the top, none in the band, all k counted."""
    k = members.shape[1]
    top = float(band[1]) if band is not None else float(tol)
    lab_r = label[rows]
    ent = np.argsort(lab_r, kind='stable')
    lab_s = lab_r[ent]
    loc = np.empty(label.shape[0], dtype=np.intp)
    left = []
    for i, mem in enumerate(members):
        b0 = label[mem[0]]
        sel = ent[np.searchsorted(lab_s, b0, 'left'):np.searchsorted(lab_s, b0, 'right')]
        loc[mem] = np.arange(k)
        gr, gc = rows[sel], cols[sel]
        b = (C[gr, gc] if keep is None else C[keep[gr], keep[gc]]).astype(float)
        bound = _full_rank_bound(k, loc[gr], loc[gc], b)
        need = (top + 8.0 * k * _EPS * math.sqrt(float(np.dot(b, b)))) * (1.0 + 4.0 * _EPS)
        if bound > need:
            if count is True:
                _paths.COUNTS['branch.rank:proved'] += 1
            elif count:
                count[0] += 1
        else:
            left.append(i)
    return np.asarray(left, dtype=np.intp)


def _rank_by_blocks(C, tol, band=None, nz=None, keep=None, count=True):
    """The number of singular values of the square `C` above `tol`, one
    block at a time; ``None`` where one lies in the open-closed interval
    `band` ``(low, high]`` (the caller's ambiguous range).  `keep`: the
    matrix is `C`'s rows and columns `keep` (``C[ix_(keep, keep)]``, read
    without the copy); `nz`: its nonzeros' positions ``(rows, cols)`` in
    `np.nonzero`'s row-major order where the caller holds them (else
    `np.nonzero`'s scan).  `count`: the full-rank proofs counted in `_paths`,
    or in `count[0]` (a list: a caller's own tally).

    `C` permuted by the connected components of its nonzero pattern is
    block diagonal, and a block-diagonal matrix's singular values are its
    blocks' together: the count is the one a dense SVD of `C` gives, made
    without it.  A 1 x 1 block's singular value is its entry's magnitude;
    blocks of one size go to `np.linalg.svd` stacked.  ⚠ MEASURED
    (2026-10-07): the dense SVD of a power grid's `C` was 72 s of a 77.6 s
    run at n = 4681, three calls at O(n^3) -- the time the check took, not
    the step's; by blocks (capacitors to ground are blocks of one once the
    reference is out, `_reduced_count`) 1.0 s.  A `C` whose
    pattern connects everything is one block -- an SVD of n x n, 2.5 s a
    state at 2 001 unknowns, ~20 s at 4 001 -- so a block of `PROOF_MIN` or
    more is first offered THE FULL-RANK PROOF (`_unproved`,
    `_full_rank_bound`: milliseconds, linear in its nonzeros), its SVD made
    only where that settles nothing.
    """
    if keep is None:
        n = C.shape[0]
        rows, cols = np.nonzero(C) if nz is None else nz
    else:
        n = keep.shape[0]
        rows, cols = np.nonzero(C[np.ix_(keep, keep)]) if nz is None else nz
    _ncomp, label = _components(n, rows, cols)
    order = np.argsort(label, kind='stable')
    starts = np.flatnonzero(np.r_[True, np.diff(label[order]) != 0])
    sizes = np.diff(np.r_[starts, n])
    rank = 0
    for k in np.unique(sizes):
        members = order[(starts[sizes == k])[:, None] + np.arange(k)]   # (blocks, k)
        if k >= PROOF_MIN:
            ## (THE FULL-RANK PROOF: a block shown to have every singular
            ## value past the band's top by more than its SVD's rounding is
            ## counted whole, as its SVD would count it)
            left = _unproved(C, rows, cols, label, members, keep, tol, band, count)
            rank += int(k) * (members.shape[0] - left.size)
            members = members[left]
            if not members.shape[0]:
                continue
        if keep is not None:
            members = keep[members]
        if k == 1:
            sv = np.abs(C[members[:, 0], members[:, 0]])
        else:
            sv = np.linalg.svd(C[members[:, :, None], members[:, None, :]],
                               compute_uv=False)                      # (blocks, k)
        if band is not None and np.any((sv > band[0]) & (sv <= band[1])):
            return None
        rank += int((sv > tol).sum())
    return rank


#: THE SCREEN BY BLOCKS (`_svd_by_blocks`): a component of this many unknowns
#: or more is not taken -- its SVD k x k (~1 s at 2 000); the check then off
SCREEN_BLOCK_MAX = 2000


def _svd_by_blocks(C, tol, kmax=None, skip=-1):
    """The screen's SVD of the `_packed.PatternMatrix` `C` where no dense
    copy is made (THE PACKED LAYOUT past `_packed.DENSE_MAX`), by the
    connected components of its nonzeros: `C` permuted by them is block
    diagonal, its singular values its blocks' together (`_rank_by_blocks`)
    -- `(count above tol, u)`, `u` a left singular vector of the smallest
    singular value (its block's, zero elsewhere; among equal smallest the
    last block in component order: the dense SVD's `U[:, -1]` is one of a
    degenerate set there too); a 1 x 1 block's value its entry's magnitude,
    blocks of one size one stacked SVD.  None where a component has `kmax`
    unknowns or more (`SCREEN_BLOCK_MAX`).  `skip` (or -1): C without that
    node's row and column -- the reference's, past DENSE_MAX
    (`_reduced_count`) -- its node out of the count and of the smallest, the
    vector 0 there."""
    kmax = SCREEN_BLOCK_MAX if kmax is None else kmax
    n = C.pat.n
    on = C.v != 0.0
    if skip >= 0:
        on &= (C.pat.rows != skip) & (C.pat.cols != skip)
    rows, cols, vals = C.pat.rows[on], C.pat.cols[on], C.v[on]
    _ncomp, label = _components(n, rows, cols)
    order = np.argsort(label, kind='stable')
    starts = np.flatnonzero(np.r_[True, np.diff(label[order]) != 0])
    sizes = np.diff(np.r_[starts, n])
    if int(sizes.max()) >= kmax:
        return None
    nb = starts.size
    block = np.empty(n, dtype=np.intp)
    block[order] = np.repeat(np.arange(nb), sizes)
    local = np.empty(n, dtype=np.intp)
    local[order] = np.arange(n) - np.repeat(starts, sizes)
    count, low, last, u = 0, np.inf, -1, None
    for k in np.unique(sizes):
        k = int(k)
        sel = np.flatnonzero(sizes == k)
        slot = np.full(nb, -1, dtype=np.intp)
        slot[sel] = np.arange(sel.size)
        b = slot[block[rows]]
        m = b >= 0
        A = np.zeros((sel.size, k, k))
        A[b[m], local[rows[m]], local[cols[m]]] = vals[m]
        if k == 1:
            sv, U = np.abs(A[:, :, 0]), np.ones((sel.size, 1, 1))
        else:
            U, sv, _Vt = np.linalg.svd(A)
        mins = sv[:, -1]
        if skip >= 0 and k == 1 and slot[block[skip]] >= 0:
            ## (the skipped node: its entries out, a 1 x 1 block of 0 -- no
            ## singular value of C's, so neither counted nor the smallest)
            sv = sv.copy()
            sv[slot[block[skip]]] = 0.0
            mins = sv[:, -1].copy()
            mins[slot[block[skip]]] = np.inf
            if not np.isfinite(mins).any():
                continue
        count += int((sv > tol).sum())
        m0 = float(mins.min())
        j = int(np.flatnonzero(mins == m0)[-1])
        if m0 < low or (m0 == low and sel[j] > last):
            low, last = m0, int(sel[j])
            u = np.zeros(n)
            u[order[starts[last] + np.arange(k)]] = U[j, :, -1]
    return count, u


#: `_plan_nonzeros`' last plan's positions: `(its flat indices, np.unique of them)`
_UPOS = [None]



def _into_held(C, held):
    """`(C to read, the array to hold)`: at size the states' C each on the
    same pattern pages (`_tran_core.pattern_of`), each written into the
    first one's -- zero off that pattern in both, the same bytes -- so the
    dense reads after fault the untouched pages in once, not once a state;
    else `C` itself (and a `_packed.PatternMatrix`: no pages)."""
    if type(C) is _packed.PatternMatrix:
        return C, None
    pat = _tran_core.pattern_of(C)
    if pat is None:
        return C, None
    if held is not None and _tran_core.pattern_of(held) is pat:
        held.reshape(-1)[pat] = C.reshape(-1)[pat]
        return held, held
    return C, C

def _plan_nonzeros(cir, C):
    """`C`'s nonzeros from the stamp plan's positions -- ``(rows, cols,
    values)``, the positions in `np.nonzero`'s row-major order -- or None.
    The plan makes C as `bincount` over its positions (`_stamp_plan.
    assemble_matrix`), so its nonzeros lie among them: gathered there, and
    proven all of them by `count_nonzero` (a NaN a nonzero in both); a
    nonzero off them -- a matrix the plan did not make -- is None."""
    if type(C) is _packed.PatternMatrix:
        ## (THE PACKED LAYOUT: the nonzeros are its pattern's, in the order
        ## `np.nonzero` gives)
        return C.values_rowmajor()
    plan = cir.__dict__.get('_stamp_plan')
    mp = None if plan is None else plan.methods.get('C')
    n = C.shape[0]
    if mp is None or C.ndim != 2 or C.shape[1] != n or not C.flags.c_contiguous:
        return None
    ## (its positions, sorted and unique: kept for the plan's flat indices
    ## they were made of)
    kept = _UPOS[0]
    if kept is not None and kept[0] is mp.flat:
        upos = kept[1]
    else:
        upos = np.unique(mp.flat)
        _UPOS[0] = (mp.flat, upos)
    if upos.size and not 0 <= int(upos[0]) <= int(upos[-1]) < n * n:
        return None
    vals = C.ravel()[upos]
    on = vals != 0.0
    vals = vals[on]
    ## (counted at its pattern where `C` is on pattern pages)
    if vals.size != _tran_core.pattern_count(C):
        return None
    rows, cols = np.divmod(upos[on], n)
    return rows, cols, vals


def _no_dense_svd(C):
    """`C` a `_packed.PatternMatrix` of `_packed.DENSE_MAX` unknowns or more:
    no dense copy of it is made (`to_dense` refuses), so no dense SVD."""
    return type(C) is _packed.PatternMatrix and C.pat.n >= _packed.DENSE_MAX


def _off_reference(Ca, g):
    """`Ca` -- a fresh |C|, dense or a `_packed.PatternMatrix` -- with the
    reference node's diagonal (`g`, or -1: none) zeroed where it is finite:
    out of the screen's scale (|C|'s largest entry) and of its cheap level's
    diagonal.  That one entry sums every capacitance to ground; each other
    entry of the reference's row and column is one node's capacitance to
    ground, under that node's own diagonal (passive capacitors), so off the
    diagonal alone the largest entry is the reduced C's -- on all 27 gold
    decks and the IBM grids.  A non-finite one stays, as branch.c hands it
    back."""
    if g >= 0:
        if type(Ca) is _packed.PatternMatrix:
            k = Ca.pat.diagonal()[g]
            if math.isfinite(Ca.v[k]):
                Ca.v[k] = 0.0
        elif math.isfinite(Ca[g, g]):
            Ca[g, g] = 0.0
    return Ca


def _reduced_count(C, nz, tol_rel, g):
    """THE CHECK READS C WITHOUT THE REFERENCE NODE (`g`, or -1: none):
    `(the scale, the count)` of one state's `C` (dense or a PatternMatrix,
    its nonzeros `nz`) -- the scale |C|'s largest entry off the reference's
    diagonal (`_off_reference`), the count the reduced C's singular values
    above `tol_rel` of it, by blocks (`_rank_by_blocks`, the reference's
    row and column out), with no band: the reduced C is the one counted,
    not a stand-in for the full one, so no dense SVD is there to agree
    with, and the reduction's residual -- the reference's row sums every
    capacitance to ground, its rounding growing with the circuit -- does
    not enter.  The screen's SVD reads the same reduced C (`_branch_svd`).
    Counted `branch.rank:reduced`; a non-finite entry raises."""
    rows, cols, vals = nz
    a = np.abs(vals)
    if not np.isfinite(a).all():
        raise FloatingPointError('the structural rank: C not finite at a generic state')
    if g >= 0:
        a[(rows == g) & (cols == g)] = 0.0
    nrm = float(np.max(a)) if a.size else 0.0
    _paths.COUNTS['branch.rank:reduced'] += 1
    if not nrm > 0.0:
        return nrm, 0
    tol = tol_rel * nrm
    if g < 0:
        return nrm, _rank_by_blocks(C, tol, nz=(rows, cols))
    on = (rows != g) & (cols != g)
    rr, rc = rows[on], cols[on]
    keep = np.r_[0:g, g + 1:C.shape[0]]
    return nrm, _rank_by_blocks(C, tol, nz=(rr - (rr > g), rc - (rc > g)), keep=keep)


class _Cert:
    """THE CERTIFICATE's data for a run (`_cert_data`): B = C[ci, cj], r0 x r0,
    block diagonal by `bstart` (the blocks `bid`; `off` marks B's entries
    between two blocks), and the C's work for its largest block (`kmax`)."""

    __slots__ = ('bid', 'bstart', 'ci', 'cj', 'ipiv', 'kmax', 'n', 'nb', 'off', 'r0', 'work')

    def __init__(self, n, ci, cj, bstart, bid):
        self.n, self.r0, self.nb = n, len(ci), len(bstart) - 1
        self.ci = np.ascontiguousarray(ci, dtype=np.int64)
        self.cj = np.ascontiguousarray(cj, dtype=np.int64)
        self.bstart = np.ascontiguousarray(bstart, dtype=np.int64)
        self.bid = np.ascontiguousarray(bid, dtype=np.int64)
        self.off = self.bid[:, None] != self.bid[None, :]
        self.kmax = int(np.diff(self.bstart).max())
        self.work = np.empty(2 * self.kmax * self.kmax)
        self.ipiv = np.empty(self.kmax, dtype=np.int64)


def _cert_data(generic, r0, iref):
    """THE CERTIFICATE's data from the structural rank's generic `C` and its
    pattern (`generic`, `_BranchCheck._branch_structural_rank`): `_Cert`, or
    None where fewer than r0 rows and columns are structurally nonzero.

    B's rows and columns are C's structurally nonzero ones without the
    reference row; where they are more than r0, an r0 x r0 submatrix of the
    generic C by pivoted QR (r0 independent columns, then r0 independent
    rows of those) -- interlacing makes ANY submatrix a lower bound, the QR
    picks one likely to be well conditioned.  Its blocks are the connected
    pieces of the generic pattern (the three random states' union):
    singular values of a block-diagonal matrix are its blocks' together."""
    import scipy.linalg as sl
    C0, pat = generic
    nz = pat.any(axis=1) & pat.any(axis=0)
    if 0 <= iref < len(nz):
        nz[iref] = False
    core = np.flatnonzero(nz)
    if r0 <= 0 or core.size < r0:
        return None
    rows = cols = core
    if core.size > r0:
        piv = sl.qr(C0[np.ix_(core, core)], mode='r', pivoting=True)[1]
        cols = np.sort(core[piv[:r0]])
        piv = sl.qr(C0[np.ix_(core, cols)].T, mode='r', pivoting=True)[1]
        rows = np.sort(core[piv[:r0]])
    r, c = np.nonzero(pat[np.ix_(rows, cols)])
    nb, lab = _components(r0, r, c)
    order = np.argsort(lab, kind='stable')
    sizes = np.bincount(lab, minlength=nb)
    return _Cert(C0.shape[0], rows[order], cols[order], np.r_[0, np.cumsum(sizes)], lab[order])


def _certify(C, cert, tol, growth):
    """THE CERTIFICATE on `C` (`_BranchCheck._branch_certify`): `(proved,
    the smallest block's bound or None)`.  Proved: every entry of
    B = C[ci, cj] between two blocks is zero and every block's bound is
    above tol (1 + delta_k) -- |b| for a block of one, else 1 / ||B_k^-1||_F,
    delta_k = growth k eps ||B_k||_F ||B_k^-1||_F (the inverse's rounding at
    the block's own conditioning).  Then sigma_r0(C) >= sigma_min(B) =
    min_k sigma_min(B_k) >= min_k 1 / ||B_k^-1||_F > tol.
    The C's (branch.c's `certify`) operation for operation: the inverse is
    numpy's inv = dgesv on the identity, each sum of squares sequential in
    column-major order (cumsum), errors ignored as the C tests no flag."""
    B = C[np.ix_(cert.ci, cert.cj)]
    if B[cert.off].any():
        return False, None
    worst, bs = math.inf, cert.bstart
    with np.errstate(all='ignore'):
        for b in range(cert.nb):
            o, e = int(bs[b]), int(bs[b + 1])
            k, delta = e - o, 0.0
            if k == 1:
                bound = abs(float(B[o, o]))
            else:
                a = B[o:e, o:e].ravel(order='F')
                a2 = float(np.cumsum(a * a)[-1])
                try:
                    x = np.linalg.inv(B[o:e, o:e]).ravel(order='F')
                except np.linalg.LinAlgError:
                    return False, None
                f = float(np.cumsum(x * x)[-1])
                bound = 1.0 / math.sqrt(f) if f > 0.0 else math.inf
                delta = growth * k * _EPS * math.sqrt(a2) * math.sqrt(f)
            if not bound > tol * (1.0 + delta):
                return False, None
            worst = min(worst, bound)
    return True, worst


class _BranchCheck:
    """Branch detection: did the step equation have more than one root.  A
    theme of `Transient` (see `transient.py`)."""

    ## ------------------------------------------------------------------
    ## BRANCH DETECTION: did the step equation have more than one root?
    ##
    ## At a point where `rank C(x)` DROPS and the surviving dynamics REPEL,
    ## the DAE's solution is genuinely non-unique (Lamour, März & Tischendorf
    ## Thm 3.53) -- and numerically that becomes MULTIPLICITY OF ROOTS OF THE
    ## STEP EQUATION, which the Newton resolves by taking whichever one its
    ## seed is nearest, silently.  MEASURED (`benchmarks/branch_selection.py`,
    ## 2026-09-10): one netlist, one grid, one tolerance returns -1.414744, 0
    ## or +1.414744 depending on the seed alone, with a spread that does NOT
    ## shrink under refinement.  Both ingredients are needed -- with a passive
    ## conductance the equilibrium attracts and uniqueness is safe however
    ## badly `C` degenerates.  A negative conductance is not exotic: it is what
    ## an oscillator's active device is.
    ##
    ## ⚠⚠ THE OBVIOUS SCREEN DOES NOT WORK HERE: firing on the STEP MATRIX
    ## ACQUIRING A NEGATIVE EIGENVALUE fires on EVERYTHING on real MNA.  Two
    ## independent reasons, both structural: MNA WITH A VOLTAGE SOURCE IS A
    ## SADDLE-POINT SYSTEM and is indefinite by construction, and an
    ## OSCILLATOR'S `G` HAS A NEGATIVE EIGENVALUE BY DESIGN with no rank drop
    ## anywhere.  4 false fires out of 4 ordinary circuits, so it cannot gate
    ## anything.
    ##
    ## What is screened instead is the condition itself: `rank C(x)` below the
    ## STRUCTURAL rank -- what `C` has at a generic operating point.  A
    ## structurally zero row (a resistive node, a source branch) is not a drop;
    ## a capacitance that VANISHES is.  That is exactly "im D(t) is not
    ## time-invariant", and it is quiet on all four circuits above.
    ##
    ## ⚠ THE ASYMMETRY IS THE POINT: the confirmation perturbs the seed by a
    ## HEURISTIC amount, so it can MISS a second root -- but when it fires it
    ## has an actual second solution in hand and converged to it.  A warning is
    ## evidence; silence is not.
    ## History: `doc/transient_history.md`, `Transient.BRANCH_SCREEN_TOL`.
    BRANCH_SCREEN_TOL = 1e-9
    ## THE CERTIFICATE (2026-10-14, `_branch_certify`): where the cheap level
    ## fires, a proof of `rank C >= r0` before the SVD.  MEASURED on the 27
    ## gold decks: four (mosrect, ab_integ, e1480, g1310) fired the cheap
    ## level on EVERY step -- a constant small capacitance under its
    ## `1e-6 * ref` -- and the SVD never fired (0 of 5 834), yet it was most
    ## of their time (ab_integ 3.3 of 5.5 s) and it kept them out of the C
    ## loop.  The certificate proved all 5 834.  `BRANCH_CERT` False: every
    ## fired cheap level goes to the SVD, as before; `BRANCH_CERT_GROWTH`
    ## the rounding factor of its margin; a `C` above `BRANCH_CERT_MAX_N`
    ## unknowns is not kept for it (the SVD decides there).
    BRANCH_CERT = True
    BRANCH_CERT_GROWTH = 10.0
    BRANCH_CERT_MAX_N = 1000

    def _branch_structural_rank(self):
        """`rank C` at a GENERIC operating point, computed once per run.

        The comparison has to be against this and not against a running
        maximum: on the degenerate branch `C` can be identically zero for the
        whole run, so its rank never "drops" -- it was never up.  Measured:
        screening against a running maximum reads QUIET on the very fixture
        that motivated this.

        ⚠ THE RANK OF C WITHOUT THE REFERENCE NODE, at every size
        (`_reduced_count`): the step solves C without the reference's row
        and column, so its rank is the one a vanishing capacitance drops.
        Where it and the full matrix's count part -- a capacitor to ground
        within a factor of the node count of the tolerance (1 pF to ground
        against 1 mH: its singular value 2 pF in the full matrix, 1 pF in
        the reduced one, the tolerance 1 pF) -- the reduced count is taken.
        History: `doc/transient_history.md`, `_tran_branch._reduced_count`.
        """
        cached = getattr(self, '_branch_rank0', None)
        if cached is not None:
            return cached
        rng = np.random.RandomState(20260910)
        best = 0
        scale = 0.0
        ## (the certificate's generic C -- the first state's -- and the
        ## three states' pattern, `_branch_cert_data`)
        generic = None
        held = None
        for _ in range(3):
            xr = rng.uniform(-1.0, 1.0, self.cir.n)
            try:
                ## (`_C_fresh`: at size the core's pass, on pattern pages; a
                ## `_packed.PatternMatrix` as it comes -- THE PACKED LAYOUT)
                C = self._C_fresh(xr)
                if type(C) is not _packed.PatternMatrix:
                    C = np.asarray(C, dtype=float)
            except Exception as exc:                           # noqa: BLE001
                ## (the check switches itself off for this analysis: said
                ## once, not silently -- the review's F11, 2026-10-01)
                warn(
                    'transient: the branch check could not read the '
                    f'structural rank of C ({type(exc).__name__}: '
                    f'{str(exc)[:80]}) and is OFF for this analysis.', ModelWarning)
                self._branch_rank0 = (0, 0.0)
                return self._branch_rank0
            if C.size == 0:
                continue
            if self.BRANCH_CERT and C.shape[0] <= self.BRANCH_CERT_MAX_N:
                ## (the certificate's data dense: BRANCH_CERT_MAX_N is below
                ## `_packed.DENSE_MAX`)
                Cd = C.to_dense() if type(C) is _packed.PatternMatrix else C
                generic = (Cd, Cd != 0) if generic is None else (generic[0], generic[1] | (Cd != 0))
            else:
                C, held = _into_held(C, held)
            ## (C without the reference node, by blocks: `_reduced_count` --
            ## its nonzeros from the stamp plan's positions where they are all
            ## there, no n x n scan; else the dense C's own)
            nz = _plan_nonzeros(self.cir, C)
            if nz is None:
                rows, cols = np.nonzero(C)
                nz = (rows, cols, C[rows, cols])
            nrm, r = _reduced_count(C, nz, self.BRANCH_SCREEN_TOL,
                                    _tran_branch_c.ref_index(self, C.shape[0]))
            scale, best = max(scale, nrm), max(best, r)
        self._branch_generic = generic
        self._branch_rank0 = (best, scale)
        return self._branch_rank0

    def _branch_cert_kept(self):
        """`(known, data)`: THE CERTIFICATE's data where this run made them
        (keyed by the structural rank's tuple: `_branch_cert_data`), else
        `(False, None)`."""
        got = self.__dict__.get('_branch_cert')
        if got is not None and got[0] is not None and got[0] is self.__dict__.get('_branch_rank0'):
            return True, got[1]
        return False, None

    def _branch_cert_data(self):
        """THE CERTIFICATE's data for this run (`_cert_data`, `_Cert`), or
        None (off, no generic C, or no r0 x r0 submatrix) -- made once a
        run, at the first screen whose cheap level fires."""
        known, cert = self._branch_cert_kept()
        if known:
            return cert
        rk = self.__dict__.get('_branch_rank0')
        generic = self.__dict__.get('_branch_generic')
        cert = None
        if self.BRANCH_CERT and rk is not None and generic is not None:
            cert = _cert_data(generic, rk[0], self.irefnode)
        self._branch_cert = (rk, cert)
        return cert

    def _branch_certify(self, C, tol):
        """THE CERTIFICATE: True where it proves `rank C >= r0` -- the SVD's
        quiet verdict -- without the SVD (`_certify`; History:
        `doc/transient_history.md`, `Transient._branch_certify`).

        It proves the rank in exact arithmetic up to its inverses' rounding,
        so it can differ from the SVD only where the SVD's own rounding puts
        sigma_r0 across `tol` -- a rank within rounding of the threshold,
        which the screen does not resolve either.  It never fires: what it
        cannot prove, the SVD decides."""
        cert = self._branch_cert_data()
        if cert is not None and type(C) is _packed.PatternMatrix:
            C = C.to_dense()
        return (cert is not None and cert.n == C.shape[0]
                and _certify(C, cert, tol, self.BRANCH_CERT_GROWTH)[0])

    def _branch_screen(self, x):
        """`(fired, null_direction)` -- has `rank C(x)` fallen below the
        structural rank?

        ⚠ The cheap proxy runs first and is what keeps this affordable: an SVD
        every step would be the same order as the factorisation.  `C`'s
        sparsity pattern is fixed by topology, so a vanishing reactance shows
        up as a structurally nonzero DIAGONAL entry collapsing, which is
        `O(m)` to check.  The SVD only runs when that fires.
        """
        r0, scale0 = self._branch_structural_rank()
        if r0 <= 0 or scale0 <= 0.0:
            return False, None
        ## (inside the converged point's evaluation session, where the step
        ## sets one: its `jacobian_only` reads `q` and `G` here next --
        ## `Transient.solve_timestep`, `_evalhint`)
        conv = getattr(self, '_conv_session', None)
        cached = self.__dict__.get('_C_cache')
        with (_evalhint.evaluating(session=conv) if conv is not None
              else contextlib.nullcontext()):
            _C = self._C_at_state(x)
        ## (kept for the step's own assembly at this state, which follows --
        ## unless it IS the kept one: the converged point was evaluated
        ## first, `_tran_newton_c`)
        if cached is None or _C is not cached[1]:
            self._C_cache = (x, _C)
        C = _C if type(_C) is _packed.PatternMatrix else np.asarray(_C, dtype=float)
        if C.size == 0:
            return False, None
        ## ⚠ A CONSTANT `C` (every stamp in the plan constant) has one largest
        ## entry, and ONE verdict for given `r0` and scale `ref` (the running
        ## maximum with it): kept, keyed on the plan's object (a parameter
        ## write or a topology change makes a new one), the running maximum
        ## moved as the reductions would have.  MEASURED (2026-10-08): the
        ## reductions were 64 ms of an 829 ms step on a 4681-unknown power
        ## grid, its `C` assembly 14 ms.  ⚠ The lookup is not free: ~11.5 k
        ## instructions a screen in situ, +1.6 % on the 4049 oscillator
        ## (whose MOS `C` is never constant), so a NO is kept for the
        ## parameter epoch -- a stale no only forgoes the shortcut.
        const = None
        epoch = _param.ParameterDict._epoch
        if self.__dict__.get('_branch_varying') is not epoch:
            const = _stamp_plan.constant_matrix(self.cir, 'C', x)
            if const is None:
                self._branch_varying = epoch
        kept = self.__dict__.get('_branch_const')
        if const is not None and kept is not None and kept[0] is const and kept[1] == r0:
            ref = max(getattr(self, '_branch_cmax', 0.0), kept[2])
            if ref == kept[3]:
                self._branch_cmax = ref
                return kept[4]
        ## (the cheap level in C where it serves, the certificate with it:
        ## the same verdict, the SVD the Python's -- `_tran_branch_c`, the C
        ## loop plan's stage 3)
        got = _tran_branch_c.screen(self, C)
        if got is None:
            verdict, nrm = self._branch_screen_C(C, r0)
        else:
            code, nrm, ref = got
            self._branch_cmax = ref
            verdict = (_QUIET if code == _tran_branch_c.QUIET
                       else _COLLAPSED if code == _tran_branch_c.COLLAPSED
                       else self._branch_svd(C, r0, ref))
        if const is not None:
            self._branch_const = (const, r0, nrm, self._branch_cmax, verdict)
        return verdict

    def _branch_screen_C(self, C, r0):
        """`_branch_screen`'s verdict on the dense `C` and `C`'s largest
        entry, `_branch_cmax` moved: the cheap level, where it fires the
        certificate (`_branch_certify`), where neither decides the SVD
        (`_branch_svd`)."""
        ## ⚠ THE SCALE IS THE LARGEST `C` THIS RUN HAS SEEN, not the random
        ## states' (`scale0`).  Those sit anywhere in +-1 V, where a forward
        ## junction's diffusion capacitance is e^38 above any state the
        ## circuit visits, so a real junction reads as ZERO against it and
        ## the screen fires (and the confirmation re-solves) on most steps.
        ## The RANK stays structural -- a running maximum of the rank reads
        ## quiet where `C` is zero for the whole run -- and for a linear `C`
        ## the two scales are one number.
        ## History: `doc/transient_history.md`, `Transient._branch_screen`.
        ## (each reduction made once -- `|C|`, its largest entry, its
        ## diagonal, the positive part -- where they were made up to three
        ## times: the same reductions of the same arrays, so the same
        ## values; speed round 10, B3.2)
        ## ⚠ C WITHOUT THE REFERENCE NODE (`_off_reference`, `_reduced_count`):
        ## its diagonal sums every capacitance to ground, so a scale holding
        ## it grows with the circuit -- on ibmpg6t 2.28e-7 against 1e-9 for
        ## any other entry, the cheap level's threshold, 1e-6 of it, 3 % under
        ## the grid's smallest capacitor; an 8 000-node ladder of 1 pF with a
        ## 1 fF fired the cheap level and took a 174 s SVD.
        ## History: `doc/transient_history.md`, `_tran_branch._reduced_count`.
        g = _tran_branch_c.ref_index(self, C.shape[0])
        Ca = _off_reference(np.abs(C), g)
        ## (on THE PACKED LAYOUT its largest entry and diagonal are the
        ## pattern's: zero off it, the diagonal on it)
        nrm = float(Ca.max())
        ref = max(getattr(self, '_branch_cmax', 0.0), nrm)
        self._branch_cmax = ref
        d = Ca.diagonal()
        tol = self.BRANCH_SCREEN_TOL * ref
        if nrm > tol and float(d.max()) > tol:
            ## nothing has collapsed at the cheap level
            pos = d[d > 0.0]
            if (float(pos.min()) if pos.size else 0.0) > 1e-6 * ref:
                return (False, None), nrm
        if nrm <= tol:
            return (True, None), nrm     # C has collapsed entirely
        ## (the rank proved without the SVD where it can be: quiet)
        if self._branch_certify(C, tol):
            return (False, None), nrm
        return self._branch_svd(C, r0, ref), nrm

    def _branch_svd(self, C, r0, ref):
        """The screen's verdict where its cheap level fired: the rank of C
        without the reference node (`_reduced_count`'s) by the SVD against
        the running maximum `ref` -- fired with the null direction (0 at the
        reference) where it fell below `r0`: on THE PACKED LAYOUT a small
        circuit's dense copy (`branch.screen:dense`); past
        `_packed.DENSE_MAX`, where none is made, by the components of C's
        nonzeros, `_svd_by_blocks` (`branch.screen:blocks`: a component of
        `SCREEN_BLOCK_MAX` or more raises, the check off as the caller
        says)."""
        tol = self.BRANCH_SCREEN_TOL * ref
        n = C.shape[0]
        g = _tran_branch_c.ref_index(self, n)
        if _no_dense_svd(C):
            got = _svd_by_blocks(C, tol, skip=g)
            if got is None:
                raise MemoryError(f'the screen: a component of {SCREEN_BLOCK_MAX} unknowns or more')
            _paths.COUNTS['branch.screen:blocks'] += 1
            r, u = got
        else:
            Cd = C.to_dense() if type(C) is _packed.PatternMatrix else np.asarray(C)
            keep = np.r_[0:g, g + 1:n] if g >= 0 else np.arange(n)
            U, sv, _Vt = np.linalg.svd(Cd[np.ix_(keep, keep)])
            _paths.COUNTS['branch.screen:dense'] += 1
            r = int((sv > tol).sum())
            u = np.zeros(n)
            u[keep] = U[:, -1]
        if r >= r0:
            return (False, None)
        return (True, u)

    def _branch_confirm(self, func, x_res, direction):
        """Re-solve the SAME step from a perturbed seed; return the other root
        if the Newton lands somewhere materially different.

        ⚠ The perturbation magnitude is a HEURISTIC.  The other roots of the
        step equation sit an `O(sqrt(h))` distance away -- measured on the
        reference fixture, they are at 0.061 at 800 points per period and a
        seed of 0.1 reached them while 0.01 did not -- and that distance is not
        knowable in general, so several magnitudes are tried.  A miss is
        therefore possible and silence proves nothing; a fire has an actual
        second solution in hand.
        """
        xr = np.asarray(x_res, dtype=float)
        n = len(xr)
        if direction is None or len(direction) != n:
            direction = np.ones(n) / np.sqrt(n)
        scale = max(float(np.max(np.abs(xr))), 1.0)
        tol = self.par.reltol * scale
        for mag in (1.0, 0.1):
            for sign in (1.0, -1.0):
                seed = xr + sign * mag * scale * np.asarray(direction,
                                                            dtype=float)
                ## ⚠ WITH THE STEP'S OWN LIMITER.  Without it no
                ## `cir.limit` runs, a stateful device (`Diode`) stays
                ## linearised at the ACCEPTED point, and the solve and the
                ## root test below both see that linearisation: measured, a
                ## diode on the rank-dropping node at 0.85 V gave an
                ## "alternative" with residual 0 there and 489 A in the
                ## step equation itself.
                try:
                    alt, _it = self._get_nrsolver().solve_system(
                        seed, func, self.toolkit, self.par.reltol,
                        self._newton_abstol_vector_reduced(),
                        self._newton_xtol_vector_reduced(), self.par.maxiter,
                        limiter=self._newton_limiter())
                except Exception:                              # noqa: BLE001
                    continue
                alt = np.asarray(alt, dtype=float)
                if float(np.max(np.abs(alt - xr))) <= 1e3 * tol:
                    continue
                ## ⚠⚠ AND IT MUST ACTUALLY BE A ROOT.  A solver that hands the
                ## SEED BACK -- converged at iteration zero, or bailed -- looks
                ## exactly like a second solution to a pure distance test, and
                ## that is a FALSE ALARM on a default-on diagnostic.
                ## History: `doc/transient_history.md`, `Transient._branch_confirm`.
                if not self._branch_is_root(func, alt):
                    continue
                return alt
        return None

    def _branch_is_root(self, func, x):
        """Does `x` satisfy the step residual?  The alternative has to be a
        SOLUTION, not merely a different vector."""
        try:
            F, J = func(x)
        except Exception:                                      # noqa: BLE001
            return False
        F = np.asarray(F, dtype=float)
        J = np.asarray(J, dtype=float)
        scale = float(np.max(np.abs(J) @ np.abs(np.asarray(x, dtype=float)))) \
            if J.size else 0.0
        return bool(np.all(np.abs(F)
                           <= self.par.reltol * max(scale, 1.0)
                           + float(np.max(self._newton_abstol_vector()))))

    def _branch_after_solve(self, func, x_res):
        """Screen the converged point, and confirm only if the screen fires.

        ⚠ `x_res` is the REDUCED vector the solver returned; the circuit's
        `C` wants the full one, so the reference row goes back in before the
        screen and comes out again for the re-solve.
        """
        ## ⚠ THE FINDING IS WARNED AFTER THE `try`, NOT IN IT: raised inside,
        ## a warning made an error (`-W error`, the suite's own policy) was
        ## swallowed and reported as the check failing (the review's W3)
        found = None
        try:
            xf = insert_row(x_res, self.irefnode, self.toolkit)
            fired, direction = self._branch_screen(xf)
            if not fired:
                return
            self._branch_count('branch_screens')
            if direction is not None:
                direction = self.toolkit.concatenate(
                    (direction[:self.irefnode], direction[self.irefnode + 1:]))
            ## ⚠⚠ A DIAGNOSTIC THAT CHANGES THE SIMULATION IS A DEFECT.  Every
            ## speculative solve writes the devices' limiting state (`Diode`'s
            ## `_vlim`, which its `i` and `G` read), so the NEXT step would
            ## start from the alternative's linearisation.  It goes back
            ## EXACTLY, from a snapshot: a `limit(x, x)` re-sync clamps against
            ## the stored state and lands short above a junction's critical
            ## voltage (`state_snapshot`).
            ## History: `doc/transient_history.md`, `Transient._branch_after_solve`.
            _snap = self._branch_hold()
            try:
                alt = self._branch_confirm(func, x_res, direction)
            finally:
                self._branch_put_back(_snap)
            if alt is None:
                return
            self._branch_count('branch_points')
            t = float(getattr(self.epar, 't', 0.0) or 0.0)
            if not getattr(self, '_branch_warned', False):
                self._branch_warned = True
                found = (t, float(np.max(np.abs(
                    alt - np.asarray(x_res, dtype=float)))))
        except Exception as exc:                               # noqa: BLE001
            ## ⚠⚠ A DIAGNOSTIC MUST NOT BE ABLE TO FAIL A SOLVE -- it runs
            ## after the answer is in hand and only reports.  But a bare
            ## `pass` here would hide the check's own bugs (it would silently
            ## do nothing while looking healthy), so the failure is recorded
            ## and announced ONCE.
            ## History: `doc/transient_history.md`, `Transient._branch_after_solve`.
            self._warn_branch_error(exc)
        if found is not None:
            warn(
                f'transient: THE STEP EQUATION HAD MORE THAN ONE ROOT at '
                f't={found[0]:.6g} s. rank C fell below its structural value '
                'there, and re-solving the same step from a different seed '
                'converged to a DIFFERENT solution (largest component '
                f'differs by {found[1]:.3g}). The answer returned is one of '
                'several valid ones and the choice was made by the Newton '
                'seed. See branch_points in the statistics; set '
                'branch_check="off" to skip this test.', ModelWarning)

    def _branch_hold(self):
        """What a confirmation's speculative solves write, held to be put
        back EXACTLY (`_branch_put_back`): the elements' state
        (`state_snapshot`: a stateful limiter's `_vlim`) and the evaluation
        core's bypass (`_tran_core`'s `bypass_state`: each device's last full
        evaluation, which a bypassed one is linearized at).  ⚠ Under the
        bypass the confirmation's evaluations moved those references, and the
        run after the check took other values: rich3 under gear2 with the
        bypass, 3 524 steps with the check on against 3 519 with it off.
        History: `doc/transient_history.md`, `Transient._branch_hold`."""
        rec = self.__dict__.get('_tran_core')
        core = rec[1] if rec is not None else None
        hold = getattr(core, 'bypass_state', None)
        return state_snapshot(self.cir), core, (hold() if hold is not None else None)

    def _branch_put_back(self, held):
        """Put back what `_branch_hold` held."""
        snap, core, bypass = held
        state_restore(snap)
        if bypass is not None:
            core.bypass_restore(bypass)

    def _warn_branch_error(self, exc):
        """The branch check's own failure, warned ONCE a run."""
        if not getattr(self, '_branch_error', None):
            self._branch_error = repr(exc)
            warn(f'transient: the branch check itself failed '
                 f'({self._branch_error}); it is disabled for this run and the '
                 'solve is unaffected', ModelWarning)

    def _branch_after_coupled(self, stage_newton, seed0, Y, residual,
                              build=None):
        """The coupled-path branch check: screen every stage, and if any is at
        a rank drop, re-solve the WHOLE BLOCK from a perturbed seed.

        ⚠ Perturbing one stage is not enough and would be a different question:
        the three stages are unknowns of ONE Newton, so a second root of the
        block is what "the step equation had several roots" means here.
        """
        if not self._branch_on():
            return
        ## (the finding warned after the `try`, as `_branch_after_solve`'s)
        found = None
        try:
            ## ⚠ "FIRED" AND "GAVE A DIRECTION" ARE DIFFERENT ANSWERS: when
            ## `C` collapses ENTIRELY the screen returns `(True, None)` --
            ## there is no null direction to name because every direction is
            ## one -- and a `fired = direction` idiom reads it as "did not
            ## fire".
            ## History: `doc/transient_history.md`, `Transient._branch_after_coupled`.
            fired = False
            direction = None
            for Yi in Y:
                hit, dvec = self._branch_screen(Yi)
                if hit:
                    fired, direction = True, dvec
                    break
            if not fired:
                return
            self._branch_count('branch_screens')
            if stage_newton is None:
                ## built only now the screen has fired: a path that solved
                ## with its own Newton (PCNR, the transform) confirms with
                ## the dense coupled one, on the same step equation
                _ctx, stage_newton, residual, seed0 = build()
            ## ⚠⚠ THE PERTURBATION MUST NOT BE A GAUGE SHIFT.  These are
            ## FULL-WIDTH stage vectors, so a direction of `ones` moves the
            ## REFERENCE NODE too -- a common-mode shift the circuit cannot
            ## see.  The solve leaves the pinned row alone and hands it back
            ## unchanged, the "alternative" differs from the base only in that
            ## row, and its residual is EXACTLY ZERO because it is the same
            ## physical solution -- which the residual check cannot catch.
            ## History: `doc/transient_history.md`, `Transient._branch_after_coupled`.
            n = len(np.asarray(Y[0], dtype=float))
            if direction is not None:
                d = np.asarray(direction, dtype=float).copy()
            else:
                d = np.zeros(n)
                d[(self.irefnode + 1) % n] = 1.0
            if 0 <= self.irefnode < n:
                d[self.irefnode] = 0.0
            dn = float(np.linalg.norm(d))
            if dn <= 0.0:
                return
            d = d / dn
            scale = max(float(np.max([np.max(np.abs(np.asarray(Yi,
                                                               dtype=float)))
                                      for Yi in Y])), 1.0)
            base = np.array([np.asarray(Yi, dtype=float) for Yi in Y])
            _snap = self._branch_hold()
            try:
                gap = self._branch_coupled_scan(stage_newton, seed0, base, d,
                                                scale, residual)
            finally:
                ## ⚠ ALWAYS, however the scan leaves: every speculative solve
                ## has written the devices' limiting state and the bypass's
                self._branch_put_back(_snap)
            if gap is None:
                return
            self._branch_count('branch_points')
            if not getattr(self, '_branch_warned', False):
                self._branch_warned = True
                found = (float(getattr(self.epar, 't', 0.0) or 0.0), gap)
        except Exception as exc:                               # noqa: BLE001
            self._warn_branch_error(exc)
        if found is not None:
            warn(
                f'transient: THE COUPLED STAGE SYSTEM HAD MORE THAN ONE ROOT '
                f'at t={found[0]:.6g} s. rank C fell below its structural '
                'value at one of the stages, and re-solving the block from a '
                'different seed converged to a DIFFERENT solution (largest '
                f'component differs by {found[1]:.3g}). The answer returned is '
                'one of several valid ones. Set branch_check="off" to skip '
                'this test.', ModelWarning)

    def _branch_coupled_scan(self, stage_newton, seed0, base, d, scale,
                             residual):
        """Try the perturbations; return the gap to a genuine second root, or
        `None`.  Split out so the caller can restore limiting state in a
        `finally` however this leaves."""
        for mag in (1.0, 0.1):
            for sign in (1.0, -1.0):
                seed = [np.asarray(seed0[i], dtype=float)
                        + sign * mag * scale * d for i in range(len(seed0))]
                try:
                    alt = stage_newton(seed)
                except Exception:                              # noqa: BLE001
                    continue
                altm = np.array([np.asarray(a, dtype=float) for a in alt])
                gap = float(np.max(np.abs(altm - base)))
                if gap <= 1e3 * self.par.reltol * scale:
                    continue
                ## ⚠⚠ AND IT MUST ACTUALLY BE A ROOT.  `_stage_newton` can
                ## hand the SEED BACK, which a pure distance test reads as a
                ## second solution.  ⚠ A FIXED-POINT test does not catch that:
                ## re-solving from a seed the solve hands back returns it again
                ## and the test passes vacuously.  So the BLOCK RESIDUAL is
                ## assembled and measured -- the alternative has to solve the
                ## equations, not merely be a different vector.
                try:
                    r_alt = residual([a.copy() for a in altm])
                    r_base = residual([b.copy() for b in base])
                except Exception:                              # noqa: BLE001
                    continue
                if r_alt > max(1e3 * r_base, 1e-9 * scale):
                    continue
                return gap
        return None

    def _branch_on(self):
        """`branch_check` is 'on' and has not failed in this run -- a
        failure is announced once and disables the check until the next
        `solve` (`_solve` clears `_branch_error`; PSS, which never calls it,
        keeps a failure for the analysis' life)."""
        return (getattr(self, 'branch_check', 'on') == 'on'
                and not getattr(self, '_branch_error', None))

    def _branch_count(self, name):
        """Count on `statistics` when there is one, on the instance otherwise
        -- a hand-driven march has no `statistics` object."""
        stats = getattr(self, 'statistics', None)
        target = stats if stats is not None else self
        setattr(target, name, getattr(target, name, 0) + 1)
    ## `'on'` (the default) or `'off'`: after each Newton, test whether the
    ## step equation had more than one root -- see `_branch_after_solve`
    branch_check = 'on'
