"""Structural diagnostics a shooting solve consults: the topological index,
the conditioning of the algebraic block, noise on the constraints.
"""
import numpy as np
from pycircuit.circuit.analysis import defaultepar
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
import pycircuit.circuit.analysis as analysis
from pycircuit.circuit._limiting import state_restore, state_snapshot


## Element classes the topological index criterion recognises.  Anything not
## listed makes the answer PROVISIONAL rather than wrong -- see
## `topological_index`, which reports what it could not classify.
_TI_CAPACITIVE = ('C',)


_TI_VOLTAGE = ('VS', 'VSin', 'VPulse', 'VSquare')


_TI_INDUCTIVE = ('L',)


_TI_CURRENT = ('IS', 'ISin', 'IPulse', 'ISquare')


_TI_RESISTIVE = ('R', 'G', 'Gyrator', 'Transformer')


def noise_enters_constraints(C_reduced, CY_reduced, tol=1e-12):
    """`(bad, residual)` — Winkler's index-1 SDAE precondition, `im B ⊆ im C`.

    Winkler (JCAM 163:435–463, 2004) Definition 2: an SDAE is **index 1**
    when "the noise sources do not appear in the constraints" -- in circuit
    terms, the noise input must lie in the image of the capacitance matrix.
    Otherwise it is an SDAE **WITH DIRECT NOISE**, outside the class that
    theory covers.  `CY = B B^T`, so `im B = im CY`: `CY`'s columns are
    projected onto `im C` and the relative residual is compared with `tol`.

    ⚠ When it fails, white noise drives a variable fixed by a CONSTRAINT and
    filtered by nothing, so **that node has no finite variance** (e.g. a
    series tank-loss resistor with no capacitance at its node; a parasitic
    capacitor there restores the class).  It bears on the COVARIANCE, not on
    the phase (`diffusion_constant`).

    History: `doc/shooting_history.md`, `noise_enters_constraints`.
    """
    Cr = np.asarray(C_reduced, dtype=float)
    Cy = np.real(np.asarray(CY_reduced))
    scale = float(np.max(np.abs(Cy))) if Cy.size else 0.0
    if scale == 0.0:
        return False, 0.0
    sol, *_ = np.linalg.lstsq(Cr, Cy, rcond=None)
    resid = float(np.max(np.abs(Cy - Cr @ sol))) / scale
    return resid > tol, resid


def topological_index(cir):
    """`(index, info)` — the DAE index from the netlist, WITHIN A STATED CLASS.

    Estevez Schwarz & Tischendorf (IJCTA 28(2):131–162, 2000), Thm 4.1, and
    Thm 4.2 for the charge-oriented MNA this tree uses: for a nonlinear
    time-independent network **without controlled sources**, and assuming
    positive-definite element Jacobians,

        the index is 2 IF AND ONLY IF the network contains a C-V loop or an
        L-I cutset; otherwise it is 1.

    Lamour, März & Tischendorf Theorem 3.47 adds index 0: a capacitive path
    from every node to datum AND no voltage sources.  `info['cap_path_to_datum']`
    says whether that test was reached and what it found.  Cross-checks: index
    0 IFF the reduced `C` is nonsingular, and Lemma 3.45's incidence-rank form
    of both criteria
    (`test_the_topological_index_agrees_with_an_incidence_RANK_criterion`).

    ⚠ IT READS AN INDEX; IT DOES NOT CERTIFY THAT ONE EXISTS.  No frozen-`t`
    test can (Brenan, Campbell & Petzold's counterexample has every local
    pencil regular and an infinite-dimensional solution family; relayed
    from a source reading, reproduced by that session, not here).  Use it to
    READ an index on a circuit already known to be solvable, never to
    ESTABLISH solvability.

    ⚠ C-ONLY LOOPS ARE NOT COUNTED: they make `C` singular without making the
    index 2 -- index 2 needs a VOLTAGE SOURCE fixing the loop.  (The paper's
    C-only remark, a relayed quote, is about Chua & Lin's variable set, not
    MNA.)

    ⚠ A DIAGNOSTIC, NOT A REFUSAL: `index > 1` does not predict which
    integrator fails (roadmap C4).  What it buys is a better message, so it
    LOCALISES: `info['loop']` and `info['cutset']` name the ELEMENTS, not just
    the verdict.

    ⚠ **V-ONLY LOOPS AND I-ONLY CUTSETS ARE NOT INDEX 2.**  They over-determine
    KVL / KCL, so MNA is STRUCTURALLY SINGULAR -- a netlist error, not a
    solver problem.  They come back as `info['v_loop']`, `info['i_cutset']`
    and `info['ill_posed']`, and **`index` is then `None`**.

    `info['unclassified']` lists elements outside the covered class —
    controlled sources above all, which the theorem excludes.  When it is
    non-empty `info['provisional']` is True, and ⚠ `provisional` is NOT an
    index-2 verdict at all: with controlled sources the index may exceed 2
    and may depend on element VALUES.  The criterion is reported rather than
    withheld, because a named assumption beats a silent refusal.

    History: `doc/shooting_history.md`, `topological_index`.
    """
    nodes = list(cir.nodes)
    nn = len(nodes)
    nmap = cir.elementnodemap
    kinds, terms, unclassified = {}, {}, []
    for name in cir.elements:
        el = cir[name]
        cls = type(el).__name__
        idx = [int(i) for i in np.asarray(nmap[name]).ravel() if int(i) < nn]
        terms[name] = sorted(set(idx))
        if cls in _TI_CAPACITIVE:
            kinds[name] = 'C'
        elif cls in _TI_VOLTAGE:
            kinds[name] = 'V'
        elif cls in _TI_INDUCTIVE:
            kinds[name] = 'L'
        elif cls in _TI_CURRENT:
            kinds[name] = 'I'
        elif cls in _TI_RESISTIVE:
            kinds[name] = 'R'
        else:
            kinds[name] = '?'
            unclassified.append('%s (%s)' % (name, cls))

    def forest(names, only_close_on=None):
        """Union-find over `names`; returns (parent, closing edge, tree, find).

        `only_close_on` restricts which KIND may be reported as the closing
        element: a cycle closed by anything else is skipped rather than
        reported, so the caller gets a loop guaranteed to contain that kind.
        """
        par = list(range(nn))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        closing = None
        tree = []
        for nm in names:
            t = terms[nm]
            if len(t) < 2:
                continue
            ra, rb = find(t[0]), find(t[1])
            if ra == rb:
                if closing is None and (only_close_on is None
                                        or kinds[nm] == only_close_on):
                    closing = nm
            else:
                par[ra] = rb
                tree.append(nm)
        return par, closing, tree, find

    ## ---- C-V loop (C-only loops included by construction) --------------
    ## ⚠ CAPACITORS FIRST, THEN SOURCES ONE AT A TIME, so the closing
    ## element is ALWAYS a voltage source: union-find returns only the FIRST
    ## closing edge, and a C-only loop closing first would hide a C-V loop.
    cv = ([nm for nm in cir.elements if kinds[nm] == 'C']
          + [nm for nm in cir.elements if kinds[nm] == 'V'])
    _p, closing, tree, _f = forest(cv, only_close_on='V')
    def close_the_loop(closing, tree):
        """The closing element plus the forest path joining its endpoints.

        That path IS the loop, and naming it is the localisation the whole
        criterion exists for -- "the opportunity to LOCALIZE critical element
        modellings", in the authors' words.
        """
        if closing is None:
            return []
        adj = {}
        for nm in tree:
            a, b = terms[nm][0], terms[nm][1]
            adj.setdefault(a, []).append((b, nm))
            adj.setdefault(b, []).append((a, nm))
        src, dst = terms[closing][0], terms[closing][1]
        seen, stack = {src: None}, [src]
        while stack:
            u = stack.pop()
            if u == dst:
                break
            for v, nm in adj.get(u, ()):
                if v not in seen:
                    seen[v] = (u, nm)
                    stack.append(v)
        node, path = dst, []
        while seen.get(node):
            u, nm = seen[node]
            path.append(nm)
            node = u
        return [closing] + list(reversed(path))

    loop = close_the_loop(closing, tree)

    ## ---- L-I cutset ----------------------------------------------------
    ## An L-I cutset exists exactly when deleting every L and I branch
    ## disconnects something the full network joins.
    def components(names):
        par = list(range(nn))

        def find(a):
            while par[a] != a:
                par[a] = par[par[a]]
                a = par[a]
            return a
        for nm in names:
            t = terms[nm]
            if len(t) >= 2:
                par[find(t[0])] = find(t[1])
        return {find(i) for i in range(nn)}, find

    allnm = list(cir.elements)
    full, _ = components(allnm)
    kept = [nm for nm in allnm if kinds[nm] not in ('L', 'I')]
    reduced, rfind = components(kept)
    cutset = []
    if len(reduced) > len(full):
        ## the L/I branches that bridge two different reduced components are
        ## the offending ones
        cutset = [nm for nm in allnm if kinds[nm] in ('L', 'I')
                  and len(terms[nm]) >= 2
                  and rfind(terms[nm][0]) != rfind(terms[nm][1])]

    ## ---- V-only loops and I-only cutsets: NOT index 2, ILL-POSED --------
    ## A loop of voltage sources over-determines KVL and a cutset of current
    ## sources over-determines KCL: MNA is STRUCTURALLY SINGULAR (a netlist
    ## error, not a higher index), so these are reported separately.  The V
    ## loop is localised the same way as the C-V loop, naming every source
    ## in it.
    _pv, v_close, v_tree, _vf = forest(
        [nm for nm in cir.elements if kinds[nm] == 'V'], only_close_on='V')
    v_loop = close_the_loop(v_close, v_tree)
    i_only = [nm for nm in cir.elements if kinds[nm] == 'I']
    i_cutset = []
    if i_only:
        kept_i = [nm for nm in allnm if kinds[nm] != 'I']
        red_i, rfind_i = components(kept_i)
        if len(red_i) > len(full):
            i_cutset = [nm for nm in i_only if len(terms[nm]) >= 2
                        and rfind_i(terms[nm][0]) != rfind_i(terms[nm][1])]

    ## ⚠ AN ILL-POSED NETLIST HAS NO INDEX: `index` is `None`, never 2, and
    ## the offending set is reported as `v_loop` / `i_cutset`, never as a
    ## C-V loop (a pure-V loop closes on a source in the C-V search too).
    if v_loop or i_cutset:
        return None, {'loop': [], 'cutset': [],
                      'v_loop': v_loop, 'i_cutset': i_cutset,
                      'ill_posed': True,
                      'kinds': kinds, 'unclassified': unclassified,
                      'provisional': bool(unclassified)}
    index = 2 if (loop or cutset) else 1
    ## ⚠ THE INDEX-0 RUNG, Theorem 3.47: a capacitive path from every node
    ## to datum AND no voltage sources.  With no voltage sources the only
    ## branch-current unknowns are inductive, and each carries `L di/dt` in
    ## `q`, so every row of `C` has a reactive entry and `C` is nonsingular.
    ## An INDUCTOR does not spoil it: the flux term makes that row
    ## differential, not algebraic.
    cap_to_datum = None
    ## ⚠ A NEGATIVE CLAIM CANNOT BE MADE FROM A PARTIAL READING.  "No
    ## voltage sources" passes VACUOUSLY for an unclassified element (`'?'`),
    ## and a `VCVS` is exactly an unrecognised voltage source.  Blocking on
    ## ANY unclassified element is too strict, though: the van der Pol's
    ## nonlinear conductance is unclassified and it is genuinely index 0.
    ## So the absence is established from the MNA dimension, which is
    ## complete: a voltage source of any kind adds a branch-current unknown
    ## (a VCCS, a current source or a conductance does not), so
    ## `cir.n - len(cir.nodes)` in excess of the `V` and `L` count is an
    ## unclassified element that could be a voltage source.
    ## ⚠ Conservative: an unclassified branch unknown that is NOT a voltage
    ## source (a transformer, an ammeter) blocks the rung too and reports 1
    ## where 0 is true.  The presence criteria (C-V loop, L-I cutset) need
    ## no such guard; `provisional` then says only that there may be more.
    branch_unknowns = cir.n - len(cir.nodes)
    accounted = sum(1 for nm in cir.elements if kinds[nm] in ('V', 'L'))
    if (index == 1 and branch_unknowns <= accounted
            and not any(kinds[nm] == 'V' for nm in cir.elements)):
        _pc, _cc, _ct, cfind = forest(
            [nm for nm in cir.elements if kinds[nm] == 'C'])
        try:
            datum = cir.get_node_index(gnd)
        except Exception:                                      # noqa: BLE001
            datum = None
        if datum is not None:
            cap_to_datum = all(cfind(i) == cfind(datum) for i in range(nn))
            if cap_to_datum:
                index = 0
    return index, {'loop': loop, 'cutset': cutset,
                   'cap_path_to_datum': cap_to_datum,
                   'v_loop': [], 'i_cutset': [],
                   'ill_posed': False,
                   'kinds': kinds, 'unclassified': unclassified,
                   'provisional': bool(unclassified)}


def algebraic_conditioning(cir, x=None, epar=None, refnode=gnd,
                           decades=8, flat_tol=1e-2, floor_k=1e3):
    """`(sigma, info)` — how well conditioned the circuit's ALGEBRAIC block is.

    `sigma` is `sigma_min(d g_2 / d y)`, the smallest singular value of the
    block that Bächle 2007 Thm 2.26 requires to have a BOUNDED INVERSE before
    a stiffly accurate method with `R(inf) = 0` is entitled to its classical
    order on an index-1 DAE.  `None` when there is no such block, or when the
    question cannot be answered on this circuit — the info dict says which.

    ⚠ It linearises at `x`, and `G(x)` is not always a pure function of `x`:
    an element carrying Newton LIMITING state stamps from solver history as
    well.  The body resets that state, starts it at `x`, and restores it
    afterwards, so the reading does not depend on solver history.

    HOW, with no index-1 splitting.  With `N = ker C` (the algebraic
    unknowns) and `Z = ker C^T` (the algebraic equations),
    `d g_2/d y = Z^T G N`, and

        sigma_min(C + h G) / h  ->  sigma_min(Z^T G N)    as h -> 0

    so SVDs down a ladder of `h` answer the question from `C` and `G` alone,
    with no basis extraction and no `(x, y)` form.  ⚠ THE CONVENTION IS OURS:
    `J = C + a h G`, not the literature's `C/h + G`, and it flips every
    exponent here.

    ⚠ THE VERDICT IS FLATNESS, NOT MAGNITUDE.  `sigma_min(C + hG) ~ h` holds
    whenever an algebraic block EXISTS (a voltage source is enough), and
    existence is not a defect.  The hypothesis fails only when
    `sigma_min(Z^T G N) -> 0`, which shows as the ratio FAILING TO SETTLE.
    Never gate on the ratio turning up again on a singular block either:
    that turn-up is pure roundoff, and on MNA's exact structural zeros it is
    not there to see.

    ⚠ THERE IS A WINDOW AND IT CAN BE EMPTY.  Above `sigma_r(C)/sigma`, where
    `sigma_r(C)` is the smallest NONZERO singular value of `C`, the ratio is
    reading the differential directions instead; below `eps*||C||/sigma` it is
    reading roundoff.  A badly conditioned `C` leaves no window at all, and
    this returns `verdict='no-window'` rather than a number.

    `info` carries `verdict`, one of:

      ``well-conditioned``   the ratio settled; `sigma` is the plateau, and
                             Thm 2.26's hypothesis holds.
      ``singular``           an algebraic block exists and its smallest
                             singular value is zero to working precision.  This
                             is `theta_0 > 0`, i.e. index >= 2 — cross-check it
                             against :func:`topological_index`.
      ``no-algebraic-block`` `C` is nonsingular, so there is nothing to condition
                             (the ratio grows like `1/h`; the ODE case).  An
                             absence claim, safe because it is read from the
                             ASSEMBLED `C` by SVD against a relative tolerance:
                             no classifier, complete data.
      ``no-window``          `C` is too ill-conditioned for a plateau to exist
                             between the turn and the roundoff floor.

    plus `h`, `ratio` and `usable` (the probe ladder and which points were not
    roundoff), `spread` (max/min over the usable points), and `window`.

    ⚠ `spread - 1` BOUNDS THE RELATIVE ERROR OF `sigma` and is the only
    accuracy statement on offer, so `flat_tol` is the worst-case accuracy the
    caller agrees to accept.

    ⚠ RESOLUTION LIMIT: at the defaults the smallest block still resolved
    (within 10 %) is about `1e-4 * ||G||`.  Both `flat_tol` (what fires) and
    `floor_k` (how many points the flatness test sees: probes below
    `floor_k * eps * ||C||` are dropped) bind, about a decade each; `decades`
    does not, because the floor guard caps the usable points.  `floor_k = 1e3`
    stays high deliberately -- the other side is a CONFIDENT WRONG NUMBER
    built from roundoff, worse than `singular`.  The limit is set by `||G||`
    alone, invariant to the capacitance unit, because the null space of `C`
    is taken by a RELATIVE tolerance and the verdict is the scale-free
    FLATNESS -- both are load-bearing.

    History: `doc/shooting_history.md`, `algebraic_conditioning`.
    """
    n = cir.n
    if epar is None:
        epar = defaultepar
    xv = np.zeros(n) if x is None else np.asarray(x, dtype=float)
    ## ⚠ RE-SYNC THE LIMITING STATE TO `xv`, AND PUT IT BACK AFTERWARDS.
    ## `G(x)` is NOT a pure function of `x` for a device with a Newton
    ## limiter -- `Diode` linearises around a stored `_vlim` -- so without
    ## this the answer depends on whatever solve ran last.  The restore is
    ## not optional: a diagnostic that changes the simulation is a defect.
    ## ⚠ RESET FIRST, THEN `limit(x, x)`.  `limit` clamps against the
    ## STORED state, so from a stale one it lands short of `x` wherever the
    ## junction is above its critical voltage (measured, a diode on an
    ## algebraic node at 0.8 V: `_vlim` 0.089 from a stored 0 V, and
    ## `sigma` 0.0210 against 0.990).  After `reset_state` the first
    ## `limit` starts the state AT `x`.
    _snap = state_snapshot(cir)
    try:
        cir.reset_state(epar)
        try:
            cir.limit(xv, xv, epar)
        except Exception:                                      # noqa: BLE001
            pass
        Cm = np.asarray(cir.C(xv, epar), dtype=float)
        Gm = np.asarray(cir.G(xv, epar), dtype=float)
    finally:
        state_restore(_snap)
    irn = cir.get_node_index(refnode)
    if irn is not None:
        Cm, Gm = remove_row_col((Cm, Gm), irn, analysis.numeric)
        Cm = np.asarray(Cm, dtype=float)
        Gm = np.asarray(Gm, dtype=float)

    sv = np.linalg.svd(Cm, compute_uv=False)
    cnorm = float(sv[0]) if len(sv) else 0.0
    gnorm = float(np.linalg.norm(Gm, 2)) if Gm.size else 0.0
    info = {'h': [], 'ratio': [], 'usable': [], 'spread': None,
            'window': None, 'c_norm': cnorm, 'g_norm': gnorm}
    if cnorm == 0.0 or gnorm == 0.0:
        info['verdict'] = 'no-window'
        return None, info

    ## The null space of `C` by a RELATIVE tolerance -- an absolute one makes
    ## the rank a function of the capacitance unit.
    tol = max(Cm.shape) * np.finfo(float).eps * cnorm
    nz = sv[sv > tol]
    if len(nz) == len(sv):
        info['verdict'] = 'no-algebraic-block'
        return None, info
    sigma_r = float(nz[-1]) if len(nz) else cnorm

    ## Probe a decade ladder below the turn.  `h_hi` is conservative because
    ## the turn sits at `sigma_r / sigma` and `sigma` is what we are after;
    ## bounding `sigma <= ||G||` puts the turn at or above `sigma_r/||G||`.
    h_hi = 1e-2 * sigma_r / gnorm
    floor = floor_k * np.finfo(float).eps * cnorm
    hs, ratios, usable = [], [], []
    for k in range(decades):
        h = h_hi * (0.1 ** k)
        smin = float(np.linalg.svd(Cm + h * Gm, compute_uv=False).min())
        hs.append(h)
        ratios.append(smin / h)
        usable.append(smin > floor)
    info['h'], info['ratio'], info['usable'] = hs, ratios, usable
    info['window'] = (floor, h_hi)

    good = [r for r, u in zip(ratios, usable) if u]
    if len(good) < 3:
        ## Either the block is zero (every probe is roundoff) or `C` is too
        ## ill-conditioned to leave a window.  Those are different answers and
        ## the rank of the block separates them without another probe.
        _u, sC, vtC = np.linalg.svd(Cm)
        Nb = vtC[sC <= tol].T
        Zb = _u[:, sC <= tol]
        blk = Zb.T @ Gm @ Nb
        sb = np.linalg.svd(blk, compute_uv=False)
        if len(sb) and sb.max() <= tol * max(1.0, gnorm / cnorm):
            info['verdict'] = 'singular'
            return 0.0, info
        info['verdict'] = 'no-window'
        return None, info

    ## ⚠ THE ERROR IS U-SHAPED, SO NEITHER END OF THE LADDER IS THE ANSWER:
    ## the limit's truncation falls like `h`, while `C`'s near-null singular
    ## values and roundoff contaminate more as `h` falls.  The FLATTEST
    ## 3-POINT WINDOW finds the turn wherever it is, and its own spread is
    ## the error bound.  (Not the median, not the last point, and not
    ## flatness over the whole ladder, whose coarse end is unconverged.)
    win = 3 if len(good) >= 3 else len(good)
    best = min((max(good[i:i + win]) / min(good[i:i + win]), i)
               for i in range(len(good) - win + 1))
    spread, at = best
    info['spread'] = spread
    info['window_at'] = at
    if spread <= 1.0 + flat_tol:
        info['verdict'] = 'well-conditioned'
        return float(np.median(good[at:at + win])), info
    info['verdict'] = 'singular'
    return 0.0, info
