"""Structural diagnostics a shooting solve consults: the topological index,
the conditioning of the algebraic block, noise on the constraints.
"""
import numpy as np
from pycircuit.circuit.analysis import defaultepar
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
import pycircuit.circuit.analysis as analysis


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
    when "the noise sources do not appear in the constraints", `im G ⊆ im A`
    — in circuit terms, the noise input must lie in the image of the
    capacitance matrix.  Otherwise it is an SDAE **WITH DIRECT NOISE**, which
    is outside the class that theory covers.

    ⚠ `CY = B B^T`, so `im B = im CY` and no factorisation is needed: project
    `CY`'s columns onto `im C` and look at the residual.

    ⚠⚠ WHAT IT MEANS WHEN IT FAILS, MEASURED (roadmap §0j).  White noise
    applied to a variable fixed by a CONSTRAINT rather than by an integrator
    is filtered by nothing, so **that node has no finite variance** — not a
    missing term and not a discretisation artefact.  A series tank-loss
    resistor with no capacitance at its node fails this; adding a parasitic
    capacitor moves the circuit back into the class.

    ⚠ The PHASE is a different question and is NOT covered by this: `c` on
    the failing fixture still agrees with an equivalent parallel-loss circuit
    to 0.999973 and with the Lyapunov route.  So this bears on the
    COVARIANCE, not on `diffusion_constant`.
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

    ⚠ A second source with the same class boundary (docs session, 2026-09-09,
    from the rendered pages of Lamour, März & Tischendorf): Lemma 3.45 gives
    the same two criteria as RANK conditions on incidence matrices --
    `[A_C A_R A_V]` full row rank iff no L-I cutset, `Q_C^T A_V` full column
    rank iff no C-V loop -- under the same hypothesis, "let all current and
    voltage sources be independent"; Theorem 3.47 adds the index-0 case
    (a capacitive path from every node to datum AND no voltage sources) and
    ⚠⚠⚠ AND IT READS AN INDEX; IT DOES NOT CERTIFY THAT ONE EXISTS.  No
    frozen-`t` / structural test can, and that is a theorem rather than a
    caution.  Brenan, Campbell & Petzold's counterexample (quoted in Estevez
    Schwarz, Lamour & März, "The common ground of DAE approaches"):

        E(t) = [[-t, t^2], [-1, t]],  F = I,  t in [-1, 1]

    has `det(lam E(t) + F(t)) = 1` exactly -- EVERY LOCAL PENCIL IS REGULAR --
    and yet `x(t) = gamma(t) [t, 1]^T` solves the homogeneous DAE for ARBITRARY
    smooth `gamma`.  An infinite-dimensional solution family, so no index is
    meaningful at all.  The failure is NON-LOCAL: the pair is pre-regular and
    the REDUCED pair has `im[E_1 F_1] = {0}`, which lives in the reduction
    sequence and not in any quantity computable at fixed `t`.  BCP's own text
    says this regularity notion "does not imply solvability".

    So: use this to READ an index on a circuit already known to be solvable,
    never to ESTABLISH solvability.  A simulator that reports "index 1"
    pointwise has not shown the problem is well posed.  (Relayed from a source
    reading; the counterexample was reproduced by that session to 4.4e-16 for
    three unrelated `gamma`, not by this one.)

    ⚠⚠ IT USED TO BE FLOORED AT 1 BY CONSTRUCTION, and that is fixed as of
    2026-09-11.  Estevez Schwarz & Tischendorf's criterion is "index 2 IF AND
    ONLY IF the network contains a C-V loop or an L-I cutset, OTHERWISE 1", so
    it could not return 0 and answered 1 for an implicit ODE -- SILENTLY, while
    the line below already recorded Theorem 3.47's index-0 case as something
    the theory "adds".  Found because a numerical probe read index 0 for a van
    der Pol and was CLAMPED to agree with this function, which was agreeing for
    the wrong reason.  Theorem 3.47's case is now implemented: `index` can be
    0, and `info['cap_path_to_datum']` says whether that test was reached and
    what it found.  ⚠ The rank condition is the independent cross-check --
    index 0 IFF the reduced `C` is NONSINGULAR, an implicit ODE -- and it is
    gated as such.

    closed-form projectors (3.61)/(3.62).  ⚠ THIS LINE SAID "not implemented"
    UNTIL 2026-09-10 AND WAS STALE BY A DAY: the rank form IS implemented, as
    `test_the_topological_index_agrees_with_an_incidence_RANK_criterion`,
    which cross-checks this graph traversal on four topologies and makes each
    condition FAIL on the one it names (an all-pass comparison would prove
    nothing).

    ⚠ A THIRD SOURCE, corroborating the class boundary below rather than
    extending it (docs session, 2026-09-10, Lamour, März & Tischendorf
    Remark 3.49): Theorem 3.48's constant-projector structure "remains valid
    also for CONTROLLED current and voltage sources IF they do not belong to
    C-V loops or L-I cutsets and their controlling voltages and currents do
    not belong to C-V loops or L-I cutsets".  That is the same exclusion the
    Estevez Schwarz & Tischendorf note below already states, reached from a
    different book -- it does not widen or narrow what this function claims.

    ⚠⚠ NOT "FROM THE NETLIST ALONE", WHICH AN EARLIER VERSION OF THIS LINE
    CLAIMED AND WHICH IS FALSE FOR CONTROLLED SOURCES.  Estevez Schwarz &
    Tischendorf close the paper by giving up BOTH halves of the criterion for
    that case: "if arbitrary controlling elements for the controlled sources
    are considered then THE INDEX OF THE NETWORK EQUATIONS MAY DEPEND ON THE
    PARAMETERS", and "if controlled sources are allowed to form a part of L-I
    cutsets or C-V loops then IT IS POSSIBLE TO BE CONFRONTED WITH HIGHER
    INDEX (> 2) PROBLEMS".

    So `provisional` is NOT a lower-confidence index-2 verdict — **it is not
    an index-2 verdict at all.**  Two independent failures at once: the index
    is no longer bounded by 2, and it is no longer a function of the topology
    at all, because it can turn on element VALUES.

    Estevez Schwarz & Tischendorf (IJCTA 28(2):131–162, 2000): for a
    nonlinear time-independent network **without controlled sources**, and
    assuming positive-definite element Jacobians,

        the index is 2 IF AND ONLY IF the network contains a C-V loop or an
        L-I cutset; otherwise it is 1.

    ⚠⚠ C-ONLY LOOPS ARE *NOT* COUNTED HERE, AGAINST THE RELAYED QUOTE, AND
    THE MEASUREMENT IS WHY.  The paper is quoted as saying "C-only loops have
    to be added to the class of C-V loops since the currents through C-only
    loops belong to the network variables whereas these currents are excluded
    in MNA formulations", and a first version of this function counted them.
    **Measured against a direct computation on our own MNA matrices — `C`'s
    null basis `N`, then the rank of `N^T G N` — all three C-only topologies
    come out INDEX 1:**

        capacitor ring touching ground     topological 2, measured 1
        floating capacitor triangle        topological 2, measured 1
        the same with every node grounded  topological 2, measured 1

    and the arithmetic is checkable by hand.  For the grounded ring
    `C = [[c1+c3, -c1], [-c1, c1+c2]]` has determinant
    `c1 c2 + c1 c3 + c2 c3 != 0` — NOT SINGULAR, so the system is not even a
    DAE.  For the floating triangle `C` IS singular (it is the triangle's
    Laplacian, null vector all-ones) but `N^T G N = (1/R)/3 != 0`, so the
    constraint is uniquely solvable and the index is 1.

    **A C-only loop makes `C` singular WITHOUT making the index 2; index 2
    needs a VOLTAGE SOURCE fixing the loop.**  ✅ RESOLVED AT THE SOURCE
    (docs session, 2026-09-08; Estévez Schwarz & Tischendorf, IJCTA 28(2)
    2000, on disk): the quote is FAITHFUL and describes Chua & Lin's
    variable set, not MNA.  Their MNA theorem, Thm 4.1 p.141, has no C-only
    clause -- "the conventional MNA leads to an index-1 DAE if and only if
    the network contains neither L-I cutsets nor C-V loops.  Otherwise ...
    index-2" -- which is exactly what this function implements.  The
    C-only sentence is Remark 4 p.143, comparing with Table 10-3-1 of Chua
    & Lin (Reference [10], the normal-tree / state-variable formulation),
    and states its own reason: "in this case, C-only loops have to be
    added to the class of C-V loops SINCE THE CURRENTS THROUGH C-ONLY LOOPS
    BELONG TO THE NETWORK VARIABLES WHEREAS THESE CURRENTS ARE EXCLUDED IN
    MNA FORMULATIONS."  So the measurement above AGREES with the theorem;
    the loop test requires a voltage source, and nothing is split.  Thm 4.2,
    immediately below, extends the same conclusions to the CHARGE-ORIENTED
    MNA -- this tree's formulation, `d/dt q(x) + i(x) + u(t) = 0` -- so the
    tree is covered by name, under 4.1's hypotheses (positive-definite C, L
    and conductance matrices, the controlled-source conditions of their
    Tables I-VI), which this docstring already carries.

    ⚠⚠ THIS IS A DIAGNOSTIC, NOT A REFUSAL, AND THE DIFFERENCE IS MEASURED.
    Roadmap C4 closed index-2 detect-and-refuse because `index > 1` is NOT
    PREDICTIVE: all three integrators converge on an L-I cutset and Gear-2
    fails on 2 of 4 index-2 topologies.  Deciding the index was never the
    obstacle — knowing it exactly and still not knowing which method to use
    is the actual state of things.  What this buys is a BETTER MESSAGE when
    something does fail.

    ⚠ WHICH IS WHY IT LOCALISES.  The authors' stated design goal is
    "topological criteria that can be checked very fast ... based on LOCAL
    assumptions, i.e. we want to provide the opportunity to LOCALIZE
    critical element modellings", for networks of ~1e7 elements where "it is
    often difficult to find the circuit configurations that lead to
    numerical difficulties".  So `info['loop']` and `info['cutset']` name the
    ELEMENTS, not just the verdict.

    ⚠ **V-ONLY LOOPS AND I-ONLY CUTSETS ARE REPORTED SEPARATELY, AND THEY ARE
    NOT INDEX 2.**  A loop of voltage sources over-determines KVL and a cutset
    of current sources over-determines KCL, so the MNA system is
    STRUCTURALLY SINGULAR — it has no solution at all, barring an exact
    cancellation.  The index criterion presumes a well-posed network.  Calling
    those "index 2" would send a reader hunting a solver problem instead of a
    netlist error, so they come back as `info['v_loop']`, `info['i_cutset']`
    and `info['ill_posed']` — and **`index` is then `None`**, because the DAE
    index presumes a solvable system and there is no honest value to give.

    `info['unclassified']` lists elements outside the covered class —
    controlled sources above all, which the theorem excludes.  When it is
    non-empty the verdict is PROVISIONAL and `info['provisional']` is True;
    the criterion is reported rather than withheld, because a named
    assumption beats a silent refusal.
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
    ## element is ALWAYS a voltage source.  A first version searched C and V
    ## together and discarded any loop that turned out to have no `V` in it
    ## -- which is wrong on a netlist carrying BOTH a C-only loop and a C-V
    ## loop, because union-find returns only the FIRST closing edge and the
    ## C-only one can close first, hiding the real one.
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
    ## ⚠ ANDREAS ASKED FOR THIS AND IT IS A DIFFERENT CATEGORY.  A loop of
    ## voltage sources over-determines KVL and a cutset of current sources
    ## over-determines KCL: the MNA system is STRUCTURALLY SINGULAR and has
    ## no solution at all (barring an exact cancellation), rather than having
    ## a higher index.  The index criterion presumes a well-posed network, so
    ## these are reported separately -- calling them "index 2" would send a
    ## reader looking for a solver problem instead of a netlist error.
    ## ⚠ the V loop is localised the SAME way as the C-V loop -- naming only
    ## the closing source would point at one of three parallel sources and
    ## leave the reader to find the rest, which is the opposite of the point.
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

    ## ⚠⚠ AN ILL-POSED NETLIST HAS NO INDEX, AND REPORTING ONE IS WORSE THAN
    ## REPORTING NOTHING.  The DAE index presumes a solvable system; a V loop
    ## or an I cutset makes MNA structurally singular, so `index` comes back
    ## `None`.  A first version returned 2 here and ALSO mislabelled the
    ## offending set -- a loop of three voltage sources was reported as a
    ## "C-V loop" containing no capacitor, because the C-V search unions
    ## capacitors then sources and a pure-V loop closes on a source.  Both
    ## symptoms point a reader at the solver when the netlist is the error,
    ## which is precisely what this split exists to prevent.
    if v_loop or i_cutset:
        return None, {'loop': [], 'cutset': [],
                      'v_loop': v_loop, 'i_cutset': i_cutset,
                      'ill_posed': True,
                      'kinds': kinds, 'unclassified': unclassified,
                      'provisional': bool(unclassified)}
    index = 2 if (loop or cutset) else 1
    ## ⚠⚠ THE INDEX-0 RUNG, Theorem 3.47 (added 2026-09-11).  Estevez Schwarz &
    ## Tischendorf's criterion is "2 iff a C-V loop or an L-I cutset, otherwise
    ## 1" and is FLOORED AT 1 BY CONSTRUCTION, so this function used to answer
    ## 1 for an implicit ODE -- SILENTLY, and its own docstring had recorded
    ## the omission all along.  MEASURED on a van der Pol (C, L and a BSource,
    ## no voltage source): its reduced `C` is NONSINGULAR, rank 2 of 2 with
    ## sigma_min 1.0, which is index 0.
    ##
    ## The condition is "a capacitive path from every node to datum AND no
    ## voltage sources": with no voltage sources the only branch-current
    ## unknowns are inductive, and each carries `L di/dt` in `q`, so every row
    ## of `C` has a reactive entry and `C` is nonsingular.
    ##
    ## ⚠ An INDUCTOR does not spoil it, which is the case a reading of the
    ## theorem's wording alone might get wrong -- the flux term makes that row
    ## differential, not algebraic.  The van der Pol is exactly that shape and
    ## the rank cross-check agrees.
    cap_to_datum = None
    ## ⚠⚠ A NEGATIVE CLAIM CANNOT BE MADE PROVISIONALLY, and that is why
    ## `not unclassified` is a condition of this rung and not of the others.
    ## Theorem 3.47's index-0 test asserts an ABSENCE -- a capacitive path from
    ## every node to datum AND NO VOLTAGE SOURCES.  An element outside the
    ## covered class is classified `'?'`, so `kinds[nm] == 'V'` is FALSE for it
    ## and the absence test passes VACUOUSLY.  A `VCVS` is exactly that: a
    ## voltage source this classifier does not recognise.
    ## ⚠ MEASURED 2026-09-11, and this rung shipped WITH the defect that
    ## morning: on the documented P1 fixture (a VCVS of gain `g` inside a C-V
    ## loop, index 2 off `g* = 1 + C2/C1` and index 3 on it) this returned
    ## INDEX 0 at every gain, where the docstring above records that every
    ## provisional fixture returned 1.  `algebraic_conditioning` says the
    ## algebraic block is SINGULAR at every gain, i.e. index >= 2, and it is
    ## right.
    ## The other criteria survive `unclassified` because they assert
    ## PRESENCE -- "I found a C-V loop" stands whatever else is in the netlist,
    ## and `provisional` then says only that there may be MORE.  An absence
    ## cannot be established from a partial reading at all.
    ##
    ## ⚠⚠ BUT "BLOCK ON ANY UNCLASSIFIED ELEMENT" IS TOO STRICT, and the van
    ## der Pol fixture is the proof: its nonlinear conductance is unclassified
    ## and it is GENUINELY INDEX 0 (recorded 2026-09-10, after a clamp-at-zero
    ## was reverted for exactly that reason).  Blocking there trades a vacuous
    ## TRUE for an avoidable FALSE.
    ##
    ## The absence can be established from COMPLETE data instead of a partial
    ## reading, which is the actual requirement.  A voltage source -- of any
    ## kind, recognised or not -- contributes a BRANCH-CURRENT unknown to MNA;
    ## a VCCS, a current source or a nonlinear conductance does not.  Measured:
    ## VS 1, VCVS 1, L 1, VCCS 0, IS 0, R 0.  So `cir.n - len(cir.nodes)` is
    ## the total branch-unknown count, `V` and `L` are the classified elements
    ## that carry one, and any EXCESS is an unclassified element that could be
    ## a voltage source.  No excess means no unrecognised voltage source can
    ## exist -- established from the MNA dimension, which is complete.
    ## ⚠ An unclassified element carrying a branch unknown that is NOT a
    ## voltage source (a transformer, an ammeter) blocks the rung too, so this
    ## reports 1 where 0 is true.  Conservative, and `provisional` already
    ## says the reading is partial; an absence asserted wrongly is the failure
    ## that has no floor.
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


def _limit_state_snapshot(cir):
    """Every element's instance dict, shallowly, so a diagnostic can put the
    limiting state back exactly as it found it."""
    out = []

    def walk(c):
        elems = getattr(c, 'elements', None)
        if not elems:
            return
        for e in elems.values():
            out.append((e, dict(e.__dict__)))
            walk(e)
    walk(cir)
    return out


def _limit_state_restore(snap):
    for elem, saved in snap:
        elem.__dict__.clear()
        elem.__dict__.update(saved)


def algebraic_conditioning(cir, x=None, epar=None, refnode=gnd,
                           decades=8, flat_tol=1e-2, floor_k=1e3):
    """`(sigma, info)` — how well conditioned the circuit's ALGEBRAIC block is.

    ⚠ IT LINEARISES AT `x`, AND `G(x)` IS NOT ALWAYS A PURE FUNCTION OF `x`.
    An element carrying Newton LIMITING state stamps from solver history as
    well: measured on `Diode`, `|G(x)|` moves by 3.6e+02 after a `limit()`
    call.  On a diode fixture the verdict and `sigma` here were UNMOVED to
    2.2e-16, because the diode sat in the differential part -- but that is one
    circuit, not a guarantee.  For a reading that cannot depend on solver
    history, pass a converged `x` or reset the circuit's state first.

    `sigma` is `sigma_min(d g_2 / d y)`, the smallest singular value of the
    block that Bächle 2007 Thm 2.26 requires to have a BOUNDED INVERSE before
    a stiffly accurate method with `R(inf) = 0` is entitled to its classical
    order on an index-1 DAE.  `None` when there is no such block, or when the
    question cannot be answered on this circuit — the info dict says which.

    HOW, and why no index-1 splitting is needed.  With `N = ker C` (the
    algebraic unknowns) and `Z = ker C^T` (the algebraic equations),
    `d g_2/d y = Z^T G N`, and

        sigma_min(C + h G) / h  ->  sigma_min(Z^T G N)    as h -> 0

    so two SVDs at different `h` answer the question from `C` and `G` alone,
    with no basis extraction and no `(x, y)` form.  ⚠ THE CONVENTION IS OURS:
    `J = C + a h G`, not the literature's `C/h + G`, and it flips every
    exponent here.  (Identification relayed from the docs session, 2026-09-11;
    verified here on six random systems with NON-symmetric `G` and mixed-rank
    `C`, recovered to six significant figures, plus both negative controls.)

    ⚠⚠ THE VERDICT IS FLATNESS, NOT MAGNITUDE, AND THAT IS THE WHOLE POINT.
    `sigma_min(C + hG) ~ h` fires whenever an algebraic block EXISTS, which is
    essentially every circuit — a voltage source alone is enough.  So the
    screen that fires on "`sigma_min ~ h`" is testing EXISTENCE, and existence
    is not a defect.  The hypothesis fails only when `sigma_min(Z^T G N) -> 0`,
    and THAT shows up as the ratio FAILING TO SETTLE.  Reading the magnitude
    instead of the flatness is what made an earlier screen point at the
    healthiest fixture we had.

    ⚠ AND THE RATIO DOES NOT SIMPLY FALL WHEN THE BLOCK IS SINGULAR.  On a
    singular block it falls while the leading term is resolved and then GROWS,
    because what is left is roundoff divided by `h`: measured 1.3e-07,
    5.2e-09, 7.7e-07, 3.0e-05 over four decades.  Anything keyed on "is it
    decreasing" reports a singular block as healthy at small enough `h`.
    ⚠⚠ BUT THOSE NUMBERS ARE FROM RANDOM DENSE SYSTEMS AND YOU WILL NOT SEE
    THEM ON A NETLIST.  `sigma_min(C + hG)` is ORTHOGONALLY INVARIANT, so in
    exact arithmetic no rotation of the same problem can change it -- the
    turn-up is therefore PURE ROUNDOFF and depends on how `ker C` happens to
    be represented.  A capacitor-free node gives MNA an EXACT structural zero,
    and there the ratio falls cleanly with no turn-up at all (docs-46,
    2026-09-11, measured on these same assembled matrices: axis-aligned zero
    no turn-up, the same problem rotated turn-up, MNA as assembled no
    turn-up).  So gate on FLATNESS, which is right in either basis, and never
    on SEEING the turn-up -- on real circuits it is not there to see.

    ⚠ THERE IS A WINDOW AND IT CAN BE EMPTY.  Above `sigma_r(C)/sigma`, where
    `sigma_r(C)` is the smallest NONZERO singular value of `C`, the ratio is
    reading the differential directions instead; below `eps*||C||/sigma` it is
    reading roundoff.  A badly conditioned `C` leaves no window at all, and
    this returns `verdict='no-window'` rather than a number — the one answer
    that must never be silently replaced by a plausible-looking value.

    `info` carries `verdict`, one of:

      ``well-conditioned``   the ratio settled; `sigma` is the plateau, and
                             Thm 2.26's hypothesis holds.
      ``singular``           an algebraic block exists and its smallest
                             singular value is zero to working precision.  This
                             is `theta_0 > 0`, i.e. index >= 2 — cross-check it
                             against :func:`topological_index`.
      ``no-algebraic-block`` `C` is nonsingular, so there is nothing to condition
                             (the ratio grows like `1/h`; the ODE case).
                             ⚠ This is an ABSENCE claim, and absence claims
                             cannot be made from a PARTIAL reading -- that is
                             what broke `topological_index`'s index-0 rung,
                             whose "no voltage sources" test passed vacuously
                             on an element it could not classify.  This one is
                             safe for a reason worth stating rather than
                             assuming: it is established from the ASSEMBLED
                             `C` by SVD against a relative tolerance, so there
                             is no classifier and nothing can be outside its
                             covered class.  Complete data, not a partial
                             reading.
      ``no-window``          `C` is too ill-conditioned for a plateau to exist
                             between the turn and the roundoff floor.

    plus `h`, `ratio` and `usable` (the probe ladder and which points were not
    roundoff), `spread` (max/min over the usable points), and `window`.

    ⚠ `spread - 1` BOUNDS THE RELATIVE ERROR OF `sigma`, and it is the only
    accuracy statement on offer -- so `flat_tol` is not a cosmetic threshold,
    it is the worst-case accuracy the caller is agreeing to accept.  Measured
    against an explicitly formed `Z^T G N`:

        fixture        sigma          explicit    rel err     spread-1
        RC             0.99005        0.99005     1.98e-12    1.96e-11
        ladder 1e16    0.001          0.001       2.04e-10    1.98e-07
        ladder 1e18    0.001          0.001       2.04e-11    1.98e-08
        e5 G22=2e-04   0.0002         0.0002      2.98e-08    2.95e-07
        e5 G22=2e-06   1.99999e-06    2e-06       3.09e-06    3.06e-05

    ⚠⚠ THERE IS A RESOLUTION LIMIT, AND A WIDTH WITHOUT ITS PARAMETERS IS
    NOT A RESULT.  On the `e5` sweep (a VCCS cancelling a node's
    self-conductance, so `G_22` passes through zero with `rank C` FIXED) the
    last value still reported is `G_22 = 2e-07` against `||G|| = 2e-03` --
    but that `1e-4 * ||G||` holds only AT `flat_tol = 1e-2` AND a 10%
    acceptance criterion, and both move it:

        flat_tol   1e-1    1e-2    1e-3    1e-4       accepted within 10%
        last G_22  2e-07   2e-07   2e-06   2e-06
        accepted   50%     10%     1e-3    1e-6       at flat_tol = 1e-2
        last G_22  2e-07   2e-07   2e-07   2e-05

    ⚠ WHAT DOES *NOT* SET IT IS THE LADDER DEPTH.  Extending `decades` from 8
    to 12, 16, 20, 24 and 30 leaves the limit at 2e-07, unmoved to every
    digit, because THE FLOOR GUARD CAPS THE USABLE POINT COUNT: 7 points
    whether `decades` is 8 or 30, every extra probe falling below
    `floor_k * eps * ||C||` and never reaching the flatness test.  That is the
    mechanism; the depth-insensitivity is its consequence.  (docs-46 predicted
    one decade of limit per decade of ladder, offered the falsification
    explicitly, and it FAILED here -- then read this code and identified the
    guard, which is the arm that did hold.)

    ⚠⚠ BOTH KNOBS BIND, ABOUT A DECADE EACH, AND A SINGLE-CAUSE STORY IS
    WRONG.  An earlier version of this docstring named `flat_tol` as the
    cause; docs-46 then named the floor guard instead, on a transcription
    where loosening `flat_tol` bought nothing.  Measured here on an
    independent transcription with both exposed, at `decades = 30`:

        shipped (floor_k=1e3, flat_tol=1e-2)    2e-07
        flat_tol loosened 100x to 1e0           2e-08
        floor guard floor_k 1e3 -> 1            2e-08
        both                                    2e-09

    so each is worth a decade and they compose.  `flat_tol` is what FIRES --
    every run stops with the spread over threshold -- and `floor_k` decides
    HOW MANY points the flatness test ever sees.  `floor_k = 1e3` is a
    three-decade safety margin over the natural roundoff floor `eps*||C||`,
    and it costs about a decade of resolution: a judgement call, exposed as a
    parameter so it can be measured rather than argued.  The default stays
    high deliberately -- the failure mode on the other side is a CONFIDENT
    WRONG NUMBER built from roundoff, which is worse than `singular`.

    ⚠ One might expect that price to vary by circuit, since the guard is keyed
    on `||C||` while the quantity it protects (`smin ~ h * sigma` near the
    plateau) carries no `||C||` at all -- so the margin, expressed in units of
    what is actually being guarded, moves with `||C||/||G||`.  It does not
    bite: the capacitance sweep above varies `||C||/||G||` over twelve decades
    at fixed `||G||` and the limit is 2e-07 throughout.  (Raised by docs-46 as
    untested; it was already covered by that sweep.)

    ⚠⚠ AND THAT LIMIT IS INVARIANT TO THE CAPACITANCE UNIT, which is the whole
    question.  An absolute rank test on these blocks smears the index-2
    crossing into a false window of width `~ tol * ||G||/||C||`, so it widens
    as `1/||C||` and at picofarads it is enormous.  Measured here across
    twelve decades of `C` -- 1e6, 1e3, 1, 1e-3, 1e-6 times nominal -- the
    limit sits at `G_22 = 2e-07` in EVERY case, exactly `1e-4` of `||G||`.
    The window is set by `||G||` alone.  Two things make that so and both are
    load-bearing: the null space of `C` is taken by a RELATIVE tolerance, and
    the verdict is the FLATNESS of the ratio, which is scale-free, rather than
    a magnitude compared against a fixed number.
    """
    n = cir.n
    if epar is None:
        epar = defaultepar
    xv = np.zeros(n) if x is None else np.asarray(x, dtype=float)
    ## ⚠⚠ RE-SYNC THE LIMITING STATE TO `xv`, AND PUT IT BACK AFTERWARDS.
    ## `G(x)` is NOT a pure function of `x` for a device with a Newton
    ## limiter -- `Diode` linearises around a stored `_vlim` -- so without
    ## this the answer depends on whatever solve ran last.  MEASURED: with a
    ## poisoned `_vlim` this routine's `sigma` moved by a RELATIVE 1.0, while
    ## DC, AC and transient were all unaffected to 0.0e+00 exactly, because a
    ## converged solve leaves the state consistent by construction.  This was
    ## the ONLY site in the tree that inherited it, and it was shipped the
    ## same day the hazard was written down -- which is the argument for the
    ## gate rather than for vigilance.
    ## `limit(x, x)` at ZERO DELTA is the documented re-sync (the same one
    ## `Transient._branch_restore_limits` uses, and the PCNR coupled step at
    ## its own convergence).  ⚠ And the restore is not optional: "A DIAGNOSTIC
    ## THAT CHANGES THE SIMULATION IS A DEFECT, AND THIS ONE DID" is recorded
    ## on that method about `branch_check`, which left `_vlim` at a
    ## speculative solve's value and moved the NEXT step's Jacobian.
    _snap = _limit_state_snapshot(cir)
    try:
        try:
            cir.limit(xv, xv, epar)
        except Exception:                                      # noqa: BLE001
            pass
        Cm = np.asarray(cir.C(xv, epar), dtype=float)
        Gm = np.asarray(cir.G(xv, epar), dtype=float)
    finally:
        _limit_state_restore(_snap)
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

    ## ⚠⚠ THE ERROR IS U-SHAPED, SO NEITHER END OF THE LADDER IS THE ANSWER.
    ## Two error sources move in opposite directions: the identification is a
    ## LIMIT, so its truncation falls like `h`; and `C`'s near-null singular
    ## values (plus roundoff) contaminate more as `h` falls.  On a clean
    ## fixture the error decreases all the way down and the LAST point is
    ## best; on a circuit carrying a singular value just under the rank
    ## tolerance the ladder passes THROUGH the true value and climbs again,
    ## and the last point is the WORST.  Two earlier versions read the median
    ## (1000x worse than attainable) and then the last point (wrong by 1e-4 on
    ## the second shape).  Take the FLATTEST 3-POINT WINDOW instead: it finds
    ## the turn wherever it is, and its own spread is the error bound.
    ## ⚠ Flatness over the WHOLE ladder is not the test either -- the coarse
    ## end is simply unconverged, and judging on it reported a perfectly
    ## healthy block (ratio converging to 2.0e-04 to ten digits) as SINGULAR
    ## because the top of its ladder was 3% off.
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
