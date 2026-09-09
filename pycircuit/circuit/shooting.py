from pycircuit.post import InternalResultDict
from .circuit import gnd
from pycircuit.circuit.analysis import *
from copy import copy
import warnings
import pycircuit.circuit.analysis as analysis
import numpy as np

def freq_analysis(x, t, rms = True, axis=-1, freqoffset = 0):
    """Return dft of equidistant sampled signal x"""
    
    npoints = np.size(x, axis)

    dt = t[1] - t[0]

    if x.dtype in (np.cdouble, np.cdouble):
        X = np.fft.fftshift(np.fft.fft(x, axis=axis),axes=(axis,)) / npoints
        freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
    else:
        freqs = np.fft.fftfreq(npoints, d=dt)[:int(np.ceil(npoints / 2.))]
        slices = [slice(None)] * x.ndim
        slices[axis] = slice(0, len(freqs))
        X = np.fft.fft(x, axis=axis)[tuple(slices)] / npoints
        ## Fold energy from negative frequencies
        X[:,1:] *= np.sqrt(2)

    if not rms:
        X *= np.sqrt(2)

    return freqs, X

## ⚠ THE SPECTRAL RADIUS CANNOT DECIDE THIS, and trying it first is the
## instructive part.  An autonomous orbit gives an eigenvalue at exactly 1
## -- but only AT its own period, and a run at any other period reads well
## below (measured 0.9615 on the quadrature phase element at the nominal
## period against 1.000226 at the corrected one).  Worse, a merely
## lightly-damped DRIVEN circuit sits near 1 too: a Q=1000 resonator has
## `exp(-pi/Q) = 0.99686`.  So no threshold separates the two -- one
## setting misses the autonomous case where users will actually run it, the
## other fires on every high-Q filter.
##
## The distinction is structural, not spectral, and it is exact: a circuit
## is autonomous when nothing in it depends on `t`.  `u(t)` is sampled
## across the period and compared; a phase accumulator driven by a DC
## source is autonomous however energetically it oscillates, which is
## precisely the case arc 5 asks about.
AUTONOMOUS_U_TOL = 1e-12


def _complex_solve(lu, b):
    """`lu.solve(b)` for a complex `b` against a REAL factorisation.

    Two back-substitutions, not a complex refactorisation: the step
    Jacobians are real, so the solve is real and splits exactly.
    """
    b = np.asarray(b)
    if np.iscomplexobj(b):
        return lu.solve(b.real) + 1j * lu.solve(b.imag)
    return lu.solve(b)


def _complex_solve_transposed(lu, b):
    """`lu.solve_transposed(b)` for complex `b` against a REAL factorisation --
    two transposed back-substitutions.  `None` if the solver cannot
    transpose (mirrors `solve_transposed`)."""
    b = np.asarray(b)
    if np.iscomplexobj(b):
        re = lu.solve_transposed(b.real)
        im = lu.solve_transposed(b.imag)
        if re is None or im is None:
            return None
        return re + 1j * im
    return lu.solve_transposed(b)



def _arnoldi_gmres(matvec, b, rtol=1e-12, maxiter=None, reortho=True):
    """GMRES that keeps its Hessenberg matrix and judges its own residual.

    Returns `(x, relres, H, k)`: the solution, the RELATIVE residual, the
    `k x k` Hessenberg matrix of the Krylov basis actually built, and `k`.

    ⚠ WRITTEN RATHER THAN IMPORTED FOR TWO REASONS, AND SPEED IS NEITHER.

    FIRST, `H` IS THE POINT.  Garcia, Romero & Acha (IEEE Trans. Power
    Systems 37(1), 2022) read the Floquet multipliers off exactly this
    matrix -- Ritz values `theta` of `I - M` map back as `lam = 1 - theta`
    -- so a GMRES that discards `H` throws away the spectrum it just
    computed.  `scipy.sparse.linalg.gmres` discards it.

    SECOND, AND THE REASON THIS IS A CORRECTNESS CHANGE: SciPy REPORTS
    BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  When the Krylov space is
    exhausted the next basis vector is numerically zero -- a HAPPY
    breakdown, where the answer is EXACT -- and it comes back as
    `info = 4`.  Trusting that flag turns an exact answer into a
    `RuntimeError`, which is what it did for AM/PM at small offsets and
    why `PAC._gmres_checked` exists to overrule it.  Here the breakdown is
    detected where it happens and returned as the converged answer it is.

    ⚠ REORTHOGONALISED ONCE BY DEFAULT.  Modified Gram-Schmidt loses
    orthogonality as the basis grows, and the Ritz values are read off `H`
    -- so a basis that has drifted gives multipliers that are wrong in a
    way the residual cannot see.  One extra pass is `O(k n)` against the
    matvec's cost, which here is a full replay of the period.

    ⚠ NO RESTARTS.  Restarting discards the basis, which is the object
    this exists to keep.  For the systems here -- `2m` unknowns, `k`
    bounded by `n` -- the full basis is affordable; a caller that needs
    restarts needs a different function and should not silently get one.
    """
    b = np.asarray(b)
    n = b.shape[0]
    kmax = int(min(n, maxiter if maxiter else n))
    beta = float(np.linalg.norm(b))
    if beta == 0.0 or kmax < 1:
        return np.zeros_like(b), 0.0, np.zeros((0, 0)), 0
    Q = [b / beta]
    H = np.zeros((kmax + 1, kmax), dtype=b.dtype)
    for j in range(kmax):
        w = np.asarray(matvec(Q[j]))
        for i in range(j + 1):
            H[i, j] = np.vdot(Q[i], w)
            w = w - H[i, j] * Q[i]
        if reortho:
            for i in range(j + 1):
                c = np.vdot(Q[i], w)
                H[i, j] += c
                w = w - c * Q[i]
        H[j + 1, j] = float(np.linalg.norm(w))
        rhs = np.zeros(j + 2, dtype=b.dtype)
        rhs[0] = beta
        y, *_ = np.linalg.lstsq(H[:j + 2, :j + 1], rhs, rcond=None)
        relres = float(np.linalg.norm(H[:j + 2, :j + 1] @ y - rhs)) / beta
        happy = H[j + 1, j] <= 1e-14 * max(beta, 1.0)
        if relres <= rtol or happy or j + 1 == kmax:
            x = np.zeros_like(b)
            for i in range(j + 1):
                x = x + y[i] * Q[i]
            return x, relres, np.array(H[:j + 1, :j + 1]), j + 1
        Q.append(w / H[j + 1, j])
    raise AssertionError('unreachable')


class FactoredPeriod(object):
    """One converged period, kept FACTORED -- the hook PAC/PPV/pnoise share.

    `PSS.solve` throws its factorisations away: they live inside the Newton's
    `build` closure and go out of scope with it.  Every small-signal analysis
    over a periodic operating point needs the same thing back -- products
    with the monodromy, never the monodromy itself -- so this is what
    `PSS.factored_period()` hands out.

    ⚠ `matvec` IS THE WHOLE INTERFACE, and deliberately so.  The withdrawn
    `PAC` read `pss.Jtvec` / `pss.Cvec` and rebuilt the `(NM)x(NM)` operator
    from them, which is where its 419.5 GiB went.  Two things are wrong with
    that route and only one of them is the memory: those two lists are
    written by `_traverse` and `_traverse_solved_history` and by NEITHER
    factored traversal, so after `solve(matrix_free=True)` they are stale or
    absent -- an analysis reading them would rebuild an operator for a
    DIFFERENT trajectory than the one that converged, silently.

    ⚠ AND THE REBUILT OPERATOR WAS EULER-SHAPED.  See
    `test_the_pac_L_is_backward_euler_only`: the old `L` has two terms per
    row, so for `trap` or `gear` it is not the discretisation the trajectory
    was produced by and `L^-1 B` is not the monodromy at all -- measured,
    spectral radius 0 against the analytic 0.8546.  A `matvec` cannot make
    that mistake, because each step carries its own `(alphas, b)`.
    """

    __slots__ = ('kind', 'opening', 'steps', 'x_last', 'x_prev', 'width',
                 'times', 'T', 'open_at_x0', '_pss')

    def __init__(self, kind, opening, steps, x_last, x_prev, pss,
                 times=None, T=None, open_at_x0=False):
        self.kind, self.opening, self.steps = kind, opening, steps
        self.x_last, self.x_prev, self._pss = x_last, x_prev, pss
        ## the grid the steps were taken on -- a forced replay needs the
        ## time of each step to evaluate `exp(j w t)` there, and reading it
        ## off `pss.times` later is exactly the parallel-indexing trap the
        ## final replay's `(t, h)` pairing was rewritten to remove
        self.times, self.T = times, T
        ## whether the period opened AT `x_0` -- no manufacturing step.  A
        ## driven replay needs to know: the manufacturing step is not in
        ## `steps`, so its share of the source is not applied, and that
        ## costs an order.  See `PAC.solve`.
        self.open_at_x0 = bool(open_at_x0)
        m = pss.cir.n - 1
        ## the solved-history map acts on the PAIR, the plain map on one
        ## state -- the same distinction `_monodromy` carries
        self.width = 2 * m if kind == 'solved_history' else m

    def matvec(self, v):
        """`M v`, real or complex, replaying the stored factors."""
        if self.kind == 'solved_history':
            return self._pss._monodromy_matvec(self.opening, self.steps, v)
        if self.kind == 'full':
            return self._pss._monodromy_matvec_full(self.steps, v)
        if self.kind == 'dirk':
            return self._pss._monodromy_matvec_dirk(self.steps, v)
        return self._pss._monodromy_matvec_plain(self.opening, self.steps, v)

    def matvec_transposed(self, v, collect=False, inject=None):
        """`M^T v` -- see `_monodromy_matvec_transposed{,_plain}`.

        ⚠ B8: the PLAIN path now has one too, for the ONE-STEP companions.
        It used to refuse outright, which made every adjoint surface --
        `ppv`, PAC, `pnoise`, `covariance` -- Gear-2 only.

        ⚠ `collect` AND `inject` FORWARD TO BOTH, which is what lets the
        callers stop naming a map. They used to reach past this method to
        `_monodromy_matvec_transposed` DIRECTLY, so B8 shipped without
        reaching any of them -- the machinery existed and every surface
        still refused. A caller that goes through here gets whichever
        recursion its `kind` calls for and cannot acquire a Gear-2
        assumption by accident.

        ⚠ THE COLLECTED STATE HAS THE MAP'S OWN WIDTH: `2m` for the pair,
        `m` for the plain one. `st[:m]` is the differential block under
        both, which is the slice every consumer wants.
        """
        if self.kind == 'solved_history':
            return self._pss._monodromy_matvec_transposed(
                self.opening, self.steps, v, collect=collect, inject=inject)
        if self.kind == 'full':
            return self._pss._monodromy_matvec_transposed_full(
                self.steps, v, collect=collect, inject=inject)
        if self.kind == 'dirk':
            return self._pss._monodromy_matvec_transposed_dirk(
                self.steps, v, collect=collect, inject=inject)
        return self._pss._monodromy_matvec_transposed_plain(
            self.opening, self.steps, v, collect=collect, inject=inject)



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
    return index, {'loop': loop, 'cutset': cutset,
                   'v_loop': [], 'i_cutset': [],
                   'ill_posed': False,
                   'kinds': kinds, 'unclassified': unclassified,
                   'provisional': bool(unclassified)}


def _cx_collect(a, b):
    """`a + 1j*b` over possibly NESTED lists of arrays (collected adjoint
    samples; a DIRK nests per-stage solves inside per-step entries)."""
    if a is None:
        ## an explicit first stage stores no solve (`D_0 = I`)
        return None
    if isinstance(a, (list, tuple)):
        return type(a)(_cx_collect(x, y) for x, y in zip(a, b))
    return np.asarray(a) + 1j * np.asarray(b)


class PSS(Analysis):
    """Periodic Steady-State using shooting Newton iterations

    The algorithm is described in [1] p65.

     1. Kenneth S. Kundert, Jacob K. White, Alberto Sangiovanni-Vincentelli
        (1990)
        Steady-State Methods for Simulating Analog and Microwave Circuits
        Kluwer Academic Publishers
        ISBN 0792390695

    **THREE CONVERGENCE CHECKS NEST HERE, AND THEY ARE NOT INTERCHANGEABLE.**

    1. the per-timestep Newton, which solves the discretised circuit
       equations at one time point;
    2. the local truncation error, which decides how far the DISCRETE
       trajectory is from the true one;
    3. the shooting Newton, which finds the periodic point of the discrete
       map.

    Two rules order them.  **(3) cannot be tighter than (1)**: the period
    map is only KNOWN to the accuracy of the per-timestep solves, so the
    shooting residual has a floor there whatever its Jacobian is.  And
    **LTE must not run per shooting iteration**: an
    adaptive grid makes the step sequence a function of `x0`, so the period
    map stops being smooth and (3) loses its quadratic rate.  Choose the
    grid once, freeze it, shoot on it.

    **(3) is a true Newton for both methods.**  It was neither, for years:
    the monodromy accumulation was missing a factor of `1/h`, and since `C`
    is singular the product collapsed to exactly zero -- so the Jacobian was
    the identity and this was successive substitution, `x0 <- phi(x0)`, on
    every circuit and without saying so.  It was found by its RATE: the
    residual fell 0.855 per iteration on a Q=20 resonator, and `exp(-pi/Q)
    = 0.8546` is that circuit's own per-period decay.  Non-convergence is
    reported now; it used to be discarded with `full_output=False`.

    Trapezoidal needed more than the factor, and Gear-2 more again.  Every
    method here writes its companion as

        iq_n = sum_k a_k q_{n-k}  +  b iq_{n-1}

    so ONE recursion differentiates all of them:

        S    = sum_{k>=1} a_k C_{n-k} Px_{n-k}  +  b Pq
        Px_n = -Jf_n^-1 S
        Pq_n = a_0 C_n Px_n + S

    Euler is `b = 0` reaching back one step, trapezoidal `b = -1` reaching
    back one, Gear-2 `b = 0` reaching back two.  The coefficients come from
    `Integrator.companion_coefficients` -- from the integrator that ACTUALLY
    ran, so an order-dropped opening step contributes its own -- rather than
    being transcribed here, which this tree has paid for three times.

    An x-only monodromy is not merely less accurate for the methods with
    memory: measured, applying the Euler form to trapezoidal converged
    SLOWER than no Jacobian at all (0.90 against 0.855).  A wrong Jacobian
    is worse than none.

    ⚠ **THE TWO FAILURE MODES ARE ORTHOGONAL, and level 3 cannot see level
    2.**  On the Q=20 resonator against a 20 V analytic peak, all three now
    converge -- and they do not agree:

        euler   5 iterations,  8.815 V   (56% low)
        gear2   7 iterations, 19.766 V   (1.2% low)
        trap    6 iterations, 19.990 V   (0.05% low)

    That is each method's own numerical damping, invisible to levels 1 and
    3, and it is why the LTE report in the recorded scope below is worth
    more than it looks: a converged shooting solve is not by itself
    evidence of a correct answer.

    (2) DOES NOT EXIST HERE, AND UNDER A SHOOTING METHOD IT CANNOT BE A
    CONTROLLER.  The grid is a fixed `linspace`; under `fixed_timestep` both
    transient backends skip the LTE verdict outright, keeping only the
    order drop that protects the integrator across a breakpoint.  Nor is
    that a wiring gap to be closed: if the step sequence adapts to `x0`
    then phi is a DIFFERENT discrete map for each `x0`, so it is not smooth
    in `x0` and the accumulated monodromy is the derivative of a
    neighbouring problem.  Freezing the grid is what makes (3) a Newton.

    So (2) changes kind here.  The estimator is still computable on a
    frozen grid, and what it measures is not convergence but ACCURACY:
    which of the three levels is limiting the answer.  That question is
    live rather than academic -- on a Q=20 resonator with `method='euler'`
    the shooting solve converges completely (5 iterations, residual
    3.9e-05) and lands at 8.815 V against a 20 V analytic peak.  The answer
    is 56% low for a reason (1) and (3) cannot see, and nothing currently
    says so.

      4. LTE AS A REPORT -- DONE (2026-09-01).  `max_lte`, `total_lte` and
         `max_lte_seam`, measured on the final replay through
         `Transient.step_lte`; see the block in `solve` for what each one
         means and why one number was not enough.  Two things it taught,
         neither of them anticipated by the paragraph above:

           - THE PER-STEP PEAK PASSES THE 56%-LOW ANSWER.  At reltol=1e-3
             euler's peak LTE on that resonator is 0.288, in tolerance,
             because the estimator bounds ONE step and the 56% is what 99
             of them do together.  A transient is right to control on that
             number; a periodic analysis cannot report only it.  The sum
             over the period reads 26.27, against gear2's 0.941 and trap's
             0.340 -- tracking amplitude errors of 55.9%, 1.17% and 0.05%.
           - THE COLD-START SEAM IS PART OF THE MAP.  `_begin_period`
             re-seeds a flat history every shooting iteration -- which is
             exactly what keeps phi a function of `x0` alone -- so the
             discrete period map opens off a past that never happened, and
             that defect is inside the map the solve converged on.

         ⚠ AND A THIRD THING, WHICH CORRECTED THE SECOND (same day, by
         measurement -- `benchmarks/pss_seam_cost.py`).  The report first
         flagged a seam for ALL THREE methods and called it the dominant
         term for the multistep pair.  Measured against the limit cycle the
         same grid and method reach with a real history, the seam costs
         **5.1e-12 V for euler and 1.3e-11 V for trapezoidal -- zero** --
         while the report was calling them 0.286 and 15.1 times tolerance.
         Only Gear-2 pays: 1.266e-01 V at 100 points/period, 54% of its
         total error, rising to 73% at 400 points because the seam falls as
         h^2 while the interior falls faster.

         The discriminator is HOW FAR THE COMPANION REACHES, not how far the
         estimator does.  Euler reads `q_{n-1}`; trapezoidal reads `q_{n-1}`
         and `iq_{n-1}`, which the order-dropped opening step supplies
         consistently.  Gear-2 reads `q_{n-2}`, which at that step is the
         entering unknown -- and shooting constrains `x(0) = x(P)`, NOT
         `x_in` to be the orbit's `x(-dt)`, so it is an O(h^2) stand-in
         being read as a history point.  Removing the order drop makes
         Gear-2 WORSE (2.34e-1 -> 4.39e-1 at 100 points): the drop is
         protective and 1.266e-01 is the residue it leaves.

         ⚠ AND ON EVERY SHIPPED METHOD IT IS `None`, so the treatment below
         describes a number the code as configured cannot produce.  Verified
         2026-09-02 across all three: the seam is collected only for a
         method whose companion reaches two charges back, which is Gear-2
         alone -- and Gear-2 always takes the solved-history path, which has
         no manufactured opening and clears the flag.  Euler and trapezoidal
         reach one.  So the figures below are reachable only through the
         superseded plain-Gear route (a test gets at them by monkeypatching
         it), and the `_limits` warning entry that quotes them is dead code
         on the shipped paths.  Kept because the reasoning is what the
         figure is FOR, and because a fourth method reaching two charges
         back on the plain path would revive both.

         So `max_lte_seam` is a FLAG, not a magnitude -- at 100 points the
         estimator's seam/interior ratio is 505x and the answer's is 1.18x.

    ⚠ THREE THINGS IN THIS CLASS ENLARGE SOMETHING, AND ALL THREE WERE
    ONCE CALLED "AUGMENTED".  One word, three referents, in one file --
    which is how a reader ends up applying a statement about one of them to
    another.  They are now named apart, and the names are worth learning:

      the `(x, iq)` MONODROMY   what trapezoidal's period map differentiates.
                                Its recursion carries a companion current, so
                                an x-only monodromy is structurally
                                incomplete.  About the DERIVATIVE, not the
                                unknowns.
      the FREE-PERIOD system    what an autonomous circuit solves: unknowns
                                `(x0, T)` with a phase condition, because the
                                period is not given.  `func_solved_history`'s
                                sibling `func_autonomous`.
      a SOLVED ENTERING HISTORY unknowns `(x0, x_{-1})`, 4b.  About where
                                the period map STARTS, not how long it runs.
                                DRIVEN circuits.  `self.solved_history`.
      the COMPOSED system       unknowns `(x0, x_{-1}, T)`, 4c -- the second
                                and third TOGETHER, for an autonomous
                                circuit under a two-step method.  ⚠ NOT a
                                synonym for either half: quoting 4b's
                                numbers as "the composed system's" is an
                                error this docstring has already made.

    ⚠ AND "THE SEAM IS REMOVED" IS NOT "THE ANSWER IS RIGHT".  Below,
    `exact` means the SEAM is gone -- the solve lands on the limit cycle the
    same grid and method reach from a real history.  It does NOT mean the
    analytic answer: 4b's Gear-2 at 100 points/period returns 19.89297 V
    against 20 V, still 1.070e-01 V out, and that residue is ordinary
    interior discretisation error which no history fix touches.  Read
    `exact` as "seam-free", never as "error-free".

      4b. A SOLVED ENTERING HISTORY FOR A TWO-STEP COMPANION (DRIVEN
         circuits; the autonomous composition is 4c) -- DONE (2026-09-01),
         and it is the remedy that measurement pointed at.  `(x_0, x_{-1})`
         are unknowns together and both must close; see
         `_traverse_solved_history`.  Applied where the COMPANION reaches two
         charges back (`_companion_reach`), which is Gear-2 alone -- euler
         and trapezoidal keep the plain path because their seam measured
         zero, and enlarging their system would double the unknowns to fix
         nothing.

         THE GATE WAS THE PREDICTION, NOT AN IMPROVEMENT.  If the seam is
         the only difference between PSS's answer and the cycle a real
         history produces, removing it must LAND on that cycle.  It does,
         to 2.5e-07 V:

              points   plain      primed     solved-hist error     gain
                 100   19.76639   19.89297   19.89297    2.34e-1 -> 1.07e-1  2.18x
                 200   19.95451   19.98524   19.98524    4.55e-2 -> 1.48e-2  3.08x
                 400   19.99008   19.99735   19.99735    9.92e-3 -> 2.65e-3  3.74x

         ⚠ Read the `error` column: 19.89297 is NOT 20 V.  The seam is gone;
         1.070e-01 V of interior discretisation error remains, untouched.

         ⚠ AND IT IS CHEAPER, WHICH WAS NOT THE EXPECTATION.  Two residual
         evaluations against twelve, 4.3x faster wall-clock on that circuit.
         The plain path seeds BOTH sensitivity rings with `I` -- which is
         the flat-history assumption written into the Jacobian -- so its
         Newton was inexact and nobody could see it, because Newton
         converges anyway from an approximate Jacobian.  The
         solved-history one is exact.

      4c. THE COMPOSED SYSTEM -- `(x_0, x_{-1}, T)`, AUTONOMOUS ONLY -- DONE
         (2026-09-01), `func_autonomous_solved_history`.  It is 4b's
         enlargement AND the free-period one at once, because an autonomous
         circuit under a two-step method needs both: the period because it
         is not given, the history because the companion reads it.  4b alone
         does not cover such a circuit and 4b's numbers are not this one's.

         ⚠ REFUSED FIRST, ON THE MEASUREMENT BELOW, THEN BUILT BECAUSE IT
         WAS ASKED FOR.  The evidence did not change; the decision did, and
         it was the owner's to make.  The measurement still stands and is
         still the reason `trap` is the default -- see the value caveat at
         the end of this item.

         WHAT THE SEAM DOES TO AN OSCILLATOR IS NOT WHAT IT DOES TO A DRIVEN
         CIRCUIT, and the guess was backwards.  The expectation was that it
         would matter MORE, the period being an unknown a per-period kick
         could land in.  Measured on the quadrature phase element it is the
         opposite: the seam moved the solved period by 2.5 ppm of a 332 ppm
         error (0.75%, and 0.38% at 400 steps) against 54% for the driven
         resonator.  The reason is one item up -- a shooting fixed point
         ABSORBS a once-per-period perturbation, and a free period is one
         more degree of freedom to absorb it into -- so the kick landed in
         the orbit's SHAPE, as a radius wobble of 2.095e-04 on an orbit of
         radius 1, instead of its frequency.

         Composing fixes both, and the wobble is the visible half:

              points   plain       free-running   composed
                 200   +329.682    +332.184       +332.185   ppm
                 400    +82.342     +82.652        +82.651   ppm
              radius wobble at 200: 2.095e-04  ->  6.6e-12

         The free-phase eigenvalue survives the enlargement -- the composed
         run reads `spectral_radius` 1.000000 -- so the autonomous
         diagnostic still says what it said.  It now reads the FULL 2m x 2m
         map, whose spectrum carries the parasitic roots of the two-step
         discretisation alongside the physical multipliers.

         ⚠ THAT WAS RECORDED AS AN OPEN WORRY AND HAS BEEN MEASURED: the
         parasitic roots are not a problem, and the reason is quantitative.
         Gear-2's parasitic root is 1/3 per STEP -- the roots of
         `1.5z^2 - 2z + 0.5` are 1 and 1/3 -- so over a period it is
         `(1/3)^N`, about 1e-95 at 200 points.  The autonomous 16x16
         spectrum measures as [1.000, 3.5e-06, 5.2e-16, 4.7e-17, 0, ...]:
         one physical unit eigenvalue and nothing else above rounding, so
         `max |eig|` picks the physical one.  A method whose parasitic root
         sat nearer the unit circle would need them separated; Gear-2's does
         not, and a test pins the gap so a future method cannot inherit the
         assumption silently.

         ⚠ THAT CLEAN SEPARATION IS A PROPERTY OF THE CIRCUITS TESTED, NOT
         OF THE METHOD, and it is known to degrade in one specific place:
         HIGH-Q OSCILLATORS, whose PHYSICAL multipliers cluster near 1.
         FOUR INDEPENDENT WITNESSES, and the earliest is Demir (IJCTA 2000)
         on a Colpitts oscillator: "the four largest eigenvalues of the
         monodromy matrix ... all four eigenvalues are +1.  THIS IS USUALLY
         THE CASE FOR HIGH-Q OSCILLATORS.  In fact, several eigenvalues can
         become very close to 1 such that they are NOT NUMERICALLY
         DISTINGUISHABLE from the one that is supposed to be equal to 1
         theoretically."  Then Bizzarri et al. -- "[shooting] is not suited
         to simulate oscillators based on very high quality resonators since
         these lead to fundamental matrices with eigenvalues very close to
         1" -- and Demir & Roychowdhury on the PPV -- "the oscillatory-mode
         eigenvalue of 1 ... cannot be distinguished from other eigenvalues
         of the system that are close to 1.  This is particularly true for
         many LC oscillators." -- and this codebase's own measurement below.

         ⚠ AND THE HISTORY SAYS WHAT TO DO ABOUT IT, which is why it is
         recorded rather than just cited.  Demir's 2000 remedy was to SELECT
         the eigenvector with the largest inner product against
         `C(0) xdot(0)` -- measured 0.2 against 1e-5, 1e-7, 2e-5.  His own
         2003 paper rejects that heuristic: "no guarantee that any of the
         candidate eigenvectors will be appreciably more orthonormal than
         the others, leading to a potential breakdown."  The SAME VECTOR
         then changes role: sampled `C(t) u_1(t)` becomes the augmented row
         `q`, so no selection happens at all.  The quantity used to CHOOSE
         among candidates becomes the CONSTRAINT that makes the candidate
         unique.  That is the whole of the 2003 improvement, and it is why
         a PPV built here must go to the augmented solve and NOT to the
         eigenvectors this method returns.

         ⚠ ONE CAUSE, AND IT SURFACES IN THREE PLACES HERE, which is why it
         is written down once instead of three times as coincidences:
           - `max |eig|` above, and `_spectral_report`'s split, both assume
             the unit root is identifiable.
           - THE PHASE ROW below removes the singularity from the unit
             eigenvalue and only that one.  Bizzarri et al. again: "this is
             not enough if any other eigenvalue is close to 1 and in this
             case possibly ill conditioned matrices must be managed."  The
             bordered system stays formally nonsingular and gets badly
             conditioned, so the failure is a slow, loud Newton rather than
             a `LinAlgError` -- the shape that invites the wrong diagnosis.
           - PPV eigen-selection, if it is ever built, picks the same root
             and inherits the same limit.
         Measured here only that the unbordered composed null space is
         EXACTLY 1-D on the circuits tested (sigma_min/sigma_next = 1.2e-11)
         -- true, and NOT a claim about the high-Q case, which has not been
         measured.

         ⚠ WHAT THE SAME MEASUREMENT DID FIND, in the DRIVEN path rather
         than the autonomous one: `_traverse_solved_history` was handing
         back `Px[0][:, :m]`, the `d x_{N-1}/d x_0` CORNER of the
         sensitivity, as the monodromy.  A corner of a sensitivity is not a
         monodromy and its eigenvalues mean nothing -- it reported
         `spectral_radius` 1.279605 for the Q=20 resonator, ABOVE ONE and so
         reading as an unstable orbit, where the analytic per-period decay
         is exp(-pi/Q) = 0.854636 and every other path reports 0.855.  Fixed
         to the pair map, which now reads 0.854833 -- closer to analytic
         than the plain path's 0.853369 -- and pinned against the analytic
         decay rather than against itself.

         ⚠ THE VALUE CAVEAT, UNCHANGED BY BUILDING IT.  Gear-2's own phase
         error is +332 ppm against trapezoidal's +83 ppm at the same grid,
         both second order, so on THIS circuit composed Gear-2 is still the
         worse choice and `method='trap'` -- the default -- remains the
         right answer.

         ⚠⚠ AND THE STIFF-CIRCUIT JUSTIFICATION DOES NOT REPRODUCE.  This
         item used to say the composition earns its place because a stiff
         autonomous circuit needs a two-step method, citing
         `doc/transient_review.md` sec. 4.6 -- trapezoidal ringing at
         `|e_n/e_{n-1}| ~ 0.9960` at `h*lambda = -1000` where Gear-2 sits at
         0.0972.  ⚠ THOSE ARE RINGDOWN NUMBERS.  They measure a TRANSIENT,
         and a periodic steady state has no transient to ring; the citation
         was carried across contexts without checking that it transfers.

         Measured on two stiff autonomous circuits: the phase element plus a
         fast RC at exactly `h*lambda = -1000`, and a diode peak detector
         whose orbit has a fast edge every period.  Trapezoidal shows NO
         ringing in either -- the alternating signature is identical between
         the methods and falls at ~h^3 under refinement, so it is the sharp
         edge being resolved, not an undamped mode -- and trapezoidal is 4x
         BETTER on frequency at both grids (+83.084 against +332.180 ppm at
         200 points; +1.287 against +5.126 at 1600).

         So no circuit measured so far prefers Gear-2 for autonomous PSS.
         The composition still earns its place on its OWN evidence and does
         not need that story: without it an autonomous Gear-2 run is
         silently biased in the period by 2.5 ppm and its orbit does not
         close in radius, so `method='gear'` was quietly WRONG there rather
         than merely inferior.  Making an offered method correct is the
         justification; "and it is the better method for stiff oscillators"
         was mine and is unsupported.

         ⚠ THE UNTRIED CASE IS NOW TRIED AND THE CLAIM IS REFUTED, not
         merely unsupported.  Van der Pol at `mu = 100` -- the canonical
         stiff relaxation oscillator, fast mode IN the orbit, measured
         stiffness ratio 5443 (edge timescale 0.0299 against a period of
         162.842412):

              method   npts        h   outcome         period      err ppm
              trap     2000   0.0815   NoConvergence        -            -
              gear     2000   0.0815   NoConvergence        -            -
              trap     8000   0.0204   NOT converged   162.813755        -
              gear     8000   0.0204   NoConvergence        -            -
              trap    20000   0.0081   converged       162.832543    -60.6
              gear    20000   0.0081   converged       162.823215   -117.9

         Trapezoidal wins on both counts: it is the only method that
         produced a finite answer at 8000 points, and at 20000 it is TWICE
         as accurate.  Across three circuits no case has been found where
         Gear-2 is the better choice for autonomous PSS, so the default
         stands and this half of the justification is closed as refuted.

         ⚠ AND THE BINDING CONSTRAINT TURNED OUT TO BE THE GRID, NOT THE
         METHOD.  Neither method runs that circuit below 20000 points,
         because this analysis freezes a UNIFORM grid (which is what makes
         (3) a Newton) and the edge needs `h < 0.01` against a period of
         162.8.  The adaptive transient that produced the reference used
         ~1160 points per period, so the uniform grid costs about 17x the
         points on this circuit class.  That is a measured argument for
         RECORDED SCOPE ITEM 5, the LTE-chosen grid -- and the first one it
         has had.

    4d. THE CHEAP APPROXIMATE ALTERNATIVE, measured and NOT shipped.  The
        4b system is seam-free but doubles the unknowns, so it is worth
        knowing what an approximation buys.  (Method H was measured against
        4b, the DRIVEN system -- not against 4c.)  `q_{-1}` can be BUILT from
        `x_0` instead of solved for, keeping the system at m unknowns.

        ⚠ AND THE PLAIN PATH IS ALREADY THE FIRST-ORDER MEMBER OF THAT
        FAMILY, which is the fact that reframes the whole question.  Its
        entering charge is EXACTLY `q_0 - h qdot_0` -- checked against the
        converged iterate to 1.6e-38 relative -- because backward Euler on
        the opening step says exactly that.  So "add a pseudo-history" is
        not an alternative to what shipped before 4b; it IS what shipped
        before 4b, at first order, and it is what measured 1.266e-01 V.

        THE STRUCTURAL FACT UNDERNEATH: a converged step satisfies
        `i(x) + iq + u = 0`, so `qdot = -(i + u)` is EXACTLY available with
        no solve, while `qddot` needs `xdot = C^-1(...)` and C is singular
        in MNA.  The DAE gives away the first derivative of the charge and
        refuses the second.  ⚠ That also explains 4b's measurement rather
        than merely restating it: TRAPEZOIDAL needs only `iq_{-1}`, which is
        that free derivative, so it is exactly initialised and has no seam
        (1.3e-11 V).  Gear-2 needs a second CHARGE, which no residual
        equation supplies.  That asymmetry is the whole story.

        Carrying it one term further, using derivatives rather than fitting
        charges (better conditioned: a quadratic fit `q_{-1} = 3q_0 - 3q_1 +
        q_2` has coefficients summing to 7 in magnitude and amplifies
        inner-solve noise):

            q_{-1} = q_0 - (3h/2) qdot_0 + (h/2) qdot_1,  error (5/12) h^3

        `x_1` is one throwaway backward-Euler predictor; its O(h^2) error
        enters with coefficient h/2 and lands at O(h^3), so a first-order
        predictor suffices.  Measured on the Q=20 resonator
        (`benchmarks/pss_seam_cost.py`, `solve_back_extrapolated`):

              points   seam plain   seam H     share plain   share H
                 100   1.266e-01    1.985e-03      1.183      0.0185
                 200   3.074e-02    1.175e-04      2.083      0.0080
                 400   7.272e-03    8.769e-06      2.742      0.0033

        64x smaller at 100 points, and falling 16.9x then 13.4x per halving
        against the plain path's 4.1x -- at least the h^3 predicted, better
        than that on these grids.  The share is what matters: the plain
        seam GROWS as a fraction of the error and method H's VANISHES.

        ⚠ AND IT DOES NOT SUBSTITUTE FOR 4c.  Measured on the autonomous
        phase element at 200 steps/period, where the frequency is the
        unknown and so the thing to watch:

              formulation        its   period ppm   err vs ref   wobble
              plain (x_in, T)     18     329.682      2.50      2.095e-04
              method H (x0, T)     3     329.720      2.46      5.8e-06
              composed 4c         14     332.185      0.001     6.6e-12

        H removes the orbit's WOBBLE -- 5.8e-06 against 5.0e-04 for a flat
        seed, so the construction itself works -- and leaves 98% of the
        FREQUENCY error.  It fixes the shape of the orbit and not its
        period.

        ⚠ THAT SPLIT CORRECTS THIS DOCSTRING'S OWN EARLIER ATTRIBUTION.  4c
        above called the 2.5 ppm "what the seam does to the period".  It is
        not: the two are separate defects that one subtraction
        (`|plain - free-running|`) had lumped together, because that
        difference contains both and can be attributed to whichever one the
        reader has in mind.

              the SEAM (entering-history error)  ->  the radius WOBBLE.
                  H fixes it, 86x.
              the MISSING CLOSURE (one equation
                  where a two-step map needs two) ->  the FREQUENCY, 2.5 ppm.
                  H does not touch it; only 4c does.

        It also sharpens the literature note: "k conditions for a k-step
        method" is NOT a statement about initialising history.  History can
        be initialised perfectly -- H does -- and the answer is still wrong
        if k-1 conditions are missing.  The two are independent
        requirements.

        THE REASON IS A MISSING EQUATION, NOT AN INACCURATE HISTORY, and
        method H is what PROVES it rather than merely suggesting it, because
        it changes exactly one variable:

              formulation   |q_-1 - q_N-2|/|q|   order   freq err   order
              plain              4.984e-04        h^2     2.502      h^3
              method H           3.005e-05        h^3     2.464      h^3
              4c (solved)        0                 --     0.001       --

        H cuts the history mismatch 16.6x and lifts its order from h^2 to
        h^3 -- and the frequency moves 1.5%.  A 16.6x better history buys
        nothing, so the history is not what is wrong.

        What is wrong is the equation count.  For a two-step method the
        discrete state is the PAIR, so periodicity is a condition on the
        pair: 4c imposes BOTH closures (`x_0 = x_{N-1}` and `x_{-1} =
        x_{N-2}`), while H imposes one and CONSTRUCTS the other state.  The
        period that closes the first component is not the period that closes
        the pair, and that gap is O(h^3) however well the history is built.  This is the literature note's
        "k conditions for a k-step method" arriving as a number: H supplies
        one (plus the phase row), 4c supplies two.  On a DRIVEN circuit the
        period is given, so the missing condition has nowhere to go and H
        does fine; on an autonomous one it goes straight into the period.

        ⚠ AND H CONVERGES FASTEST, TO THE WRONG ANSWER -- three Newton
        iterations against 4c's fourteen and the plain path's eighteen, with
        `ier == 1` and the residual satisfied.  This class already records
        that a converged shooting solve is not evidence of a correct answer
        (see the successive-substitution defect above); here it is again, in
        a formulation built the same day.  Fast convergence on a smaller
        system is not a merit when the system is missing an equation.

        SO WHY IS IT NOT SHIPPED.  ⚠ NOT ON SPEED -- that argument was made
        here and was WRONG.  Method H's 15 residual evaluations are an
        artefact of the flat Jacobian it was handed for the ACCURACY
        measurement; given a good one (finite differences) it converges in
        3, against 4b's 2.  That is not a difference worth a decision, and
        comparing an exact-Jacobian formulation against a deliberately
        crippled one was not a fair comparison.

        What is actually left against it: it would be a SECOND formulation
        for a job 4b already does seam-free, in a tree that has paid more
        than once for duplicate paths -- and its exact Jacobian is
        unwritten.  `d q_{-1}/d x_0 = C_0 + (3h/2) G_0 - (h/2) G_1
        dx_1/dx_0` needs deriving, including the predictor's own
        sensitivity; until it exists H is either slow (flat Jacobian, 15) or
        expensive (finite differences, m+1 traversals per iteration).

        For the DRIVEN case that leaves a maintenance judgement -- a second
        formulation for a job 4b already does seam-free -- which belongs to
        whoever owns the trade.  For the AUTONOMOUS case it is no longer a
        judgement at all: H is measurably wrong there, by 98% of the
        frequency error, and 4c is not optional.

        WHAT WOULD REOPEN IT: a circuit large enough that the 2m x 2m dense
        `J_phi` factorisation dominates -- Kundert puts that above a few
        hundred unknowns.  ⚠ But note where that argument really points:
        with the matrix-free Krylov solve of recorded scope item 6, the
        enlargement costs 2x (vector length), not 8x (factorisation).  So
        the scaling case is an argument for item 6 first, and only then for
        approximating the history.

    WHAT THE LITERATURE SAYS ABOUT 4b (checked 2026-09-01, because the fix
    above looked like something that ought to be standard):

      THE GENERAL RESULT IS CLASSICAL, and it is exactly what 4b hit.  A
      k-step linear multistep method turns a first-order continuous problem
      into a k-th ORDER DISCRETE one, which introduces parasitic (spurious)
      solutions and needs **k conditions** to determine the discrete
      solution -- not one.  That is why a single periodicity condition on
      `x0` is under-determined for Gear-2.  It is the founding observation
      of BOUNDARY VALUE METHODS (Brugnano & Trigiante), which supply the k
      conditions as "one initial and k-1 final", chosen at both ends
      deliberately because it improves stability.  A periodic BVP hands
      them over for free: requiring the whole k-tuple to close is what
      `func_solved_history` does.  The parasitic roots are also the spurious
      eigenvalues that appear in the 2m x 2m monodromy, which is why the
      autonomous eigenvalue-at-1 diagnostic would need redefining if the two
      systems were ever composed.  Standard reference for the BVP side:
      Ascher, Mattheij & Russell, "Numerical Solution of Boundary Value
      Problems for ODEs" (SIAM, 1995).

      THE CIRCUIT LITERATURE ASSUMES THE PROBLEM AWAY, consistently and
      reasonably.  Kundert ("Simulation Methods for RF Integrated Circuits",
      ICCAD 1997) writes the shooting map as `phi_T(v0, 0)` with the state
      `v` alone and gives the sensitivity's component pieces as
      `Jf(v(ts)) = G(v(ts)) + C(v(ts))/hs` -- a ONE-STEP, backward-Euler
      shaped Jacobian, no `q_{n-2}` term, no history in the state.  Gourary,
      Rusakov, Ulyanov & Zharov (MES 2019,
      doi:10.31114/2078-7707-2019-1-25-30) likewise solve "with respect to
      the state vector at the beginning of one period".  Both are correct
      for the one-step methods they are written for; the plain path IS that
      formulation.  The gap opens only when a two-step companion is handed
      to it.  ⚠ Note the 2019 paper's title promises more than it delivers
      here: its case for single-step is A-STABILITY ("the common drawback of
      BDF methods is the lack of A-stability for order higher than 2"), not
      history at the period boundary.

      THE CLOSEST ANYONE COMES is Wambacq, Vandersteen, Phillips,
      Roychowdhury, Eberle, Yang, Long & Demir, "CAD for RF circuits", which
      argues for one-step Chebyshev-IRK discretisation and says of it:
      "Each step is independent of the ones before and after".  That
      independence is precisely the property whose ABSENCE is the seam --
      but they argue it from stability and step adaptivity, not from
      initialisation.

      ⚠ AND THE STANDARD CASE AGAINST GEAR-2 DOES NOT APPLY HERE, which is
      the part worth keeping.  The same paper's objections are all about an
      ADAPTIVE grid: BDF "not actually as numerically stable as is popularly
      believed"; "the second order Gear method is not A-stable for
      nonuniform steps, and in fact it is not stable for any timestep if the
      ratios between consecutive steps exceed about 2.4"; "a rapid change of
      timestep in a multistep code also necessarily comes with a loss of
      order".  PSS freezes a UNIFORM grid for the whole solve -- that is
      what makes (3) a Newton, see above -- so none of those bite inside a
      run.  Shooting is the one place Gear-2 is on its best behaviour, which
      is an argument for having repaired it rather than refused it.

      WHAT IS NOT IN ANY OF IT: a number.  No source found quantifies what a
      mis-initialised history costs on a circuit.  The figures in 4b are
      this tree's own; `benchmarks/pss_seam_cost.py` is the measurement.

    4e. WHAT THE SCALAR `Idtmod` FORMS DO -- idtmod arc 5, investigated
        2026-09-01, and it turned up a defect that is NOT about `Idtmod`.

        The scalar element runs: with `ic` set, DC solves and a transient
        integrates it correctly, the state advancing one modulus per output
        period and the gauge shift (`_apply_periodic_shifts`) wrapping it
        back.  ⚠ That wrap is why shooting does not simply refuse it -- the
        state DOES return to itself once wrapped, so `x0 - phi(x0)` is
        satisfiable, and `IdtmodQuadrature`'s docstring claim that the
        scalar form "structurally cannot" close is true of the RAW state
        and not of the wrapped one this tree actually integrates.

        But asking whether the answer was right exposed the general defect:

        ⚠ AN AUTONOMOUS PERIOD IS DETERMINED ONLY UP TO AN INTEGER MULTIPLE,
        AND THE SOLVE FOLLOWS ITS SEED.  `k*T` satisfies the periodicity
        condition whenever `T` does.  Measured on the quadrature element,
        true period 1.000e-03:

              seed      solved         converged?
              1e-3      1.000083e-03   yes
              2e-3      2.000665e-03   yes
              3e-3      3.002245e-03   yes

        Every one is a correct periodic solution and every one reports
        success; the FUNDAMENTAL FREQUENCY -- usually the thing a PSS user
        wanted -- is wrong by the factor.  This is a property of the
        free-period system (4c), not of `Idtmod`, and it was silent.
        Detected now: see the recurrence check in `solve`, which needs no
        extra solve and sets `fundamental_period`.  Driven runs are exempt,
        their period being the caller's.

        ⚠ AND THE SECOND DEFECT IS ALSO NOT ABOUT `Idtmod`.  What looked
        like the scalar form's seed-fragility -- Gear-2 collapsing to
        `T ~ 1e-18`, trapezoidal raising a bare `LinAlgError` -- reproduces
        on the QUADRATURE element too: from a 1e-4 seed against a 1e-3
        fundamental, Gear-2 returns -1.5e-20 there and trapezoidal dies from
        three seeds of five.  `T = 0` is a REGULAR ROOT of every autonomous
        shooting system, because `x0 - phi_T(x0)` vanishes identically
        there and the phase condition constrains `x0`, not the period.  So
        any seed below the fundamental is drawn to it.

        Neither outcome is silent -- the collapse reports
        `converged = False` (⚠ since 2026-09-06: it previously reported
        `converged = True`, because `T = 0` is a REGULAR root that `fsolve`
        reaches cleanly and reports success on; `_free_period_solve` now
        demotes `ier` when it detects the collapse) and the exception is
        loud -- but neither named
        its cause, and the generic non-convergence advice ("raise
        maxiterations") is actively wrong for it: no number of iterations
        reaches a fundamental from below.  `_free_period_solve` now names
        both, and says what to do instead (seed at or above the expected
        period; a short transient and the interval between two output
        recurrences gives one).

        ⚠ THREE TIMES IN THIS ITEM the defect looked like the element and
        was the formulation.  The scalar form's raw state not closing is
        answered by the gauge shift; the multiple-period ambiguity and the
        trivial root are properties of the free-period system that any
        autonomous circuit has.  Arc 5's deliverable turned out to be two
        diagnostics on 4c, not a change to `Idtmod` at all.

        STILL OPEN, and genuinely about the scalar form: seeded correctly
        it converges, but a pure phase accumulator has no amplitude to pin,
        so its orbit is a ramp and nothing distinguishes one starting phase
        from another beyond the phase condition itself.  Whether that is a
        limitation worth removing has not been established.

    RECORDED SCOPE, in order, neither of these planned work yet:

      5. LTE-CHOSEN grid -- THE MECHANISM IS BUILT (2026-09-01), the
         payoff case is not yet reached.  `PSS.solve(grid=...)` takes step
         FRACTIONS of the period and freezes them; see `_period_grid`.

         ⚠ THE RECORDED BLOCKER WAS STALE.  This item said "blocked on
         `Transient` accepting a non-uniform grid; `fixed_timestep` is
         uniform-only".  That loop IS uniform-only -- and PSS never uses
         it, driving `solve_timestep` one step at a time instead, where
         non-uniform steps worked unchanged.  Verified before anything was
         written.

         Verified on benign circuits: driven and autonomous, both methods,
         2:1 and smoothly-varying grids all converge, with the driven
         spectral radius still matching exp(-pi/Q).  The autonomous case is
         what makes FRACTIONS the contract rather than times -- the grid is
         rebuilt at the current `T` on every residual evaluation, so
         `dh/dT = h/T` keeps holding.

         VAN DER POL, THE CIRCUIT THAT MOTIVATED IT, NOW SOLVES THROUGH IT
         (2026-09-02).  On its own 1105-step LTE-chosen grid, one added
         opening step: trapezoidal converges at -73.8 ppm and Gear-2 at
         -100.6, against the 20000 UNIFORM points a uniform grid needs for
         -60.6.  18x fewer points, and Gear-2 solves it for the first time
         on any grid.  Pinned by
         `test_the_lte_chosen_grid_solves_van_der_pol_through_the_analysis`.

         ⚠ THE BLOCKER WAS THE MANUFACTURED OPENING STEP, AND THIS
         DOCSTRING PREVIOUSLY NAMED IT WRONG.  It read: "what is left is
         that the throwaway used FINITE DIFFERENCES ... a stiff relaxation
         oscillator is the first circuit that cannot tolerate the plain
         path's 30%".  MEASURED: an exact finite-difference Jacobian on the
         real analysis does NOT fix it -- it fails identically -- and once
         the opening step is subdivided the analytic and finite-difference
         Jacobians agree to SIX DIGITS (-73.8 ppm both).  The ~30% was
         never what stopped this case.

         What stopped it: `_traverse` manufactures `x(0)` with ONE
         order-dropped Euler step of `hs[0]`, and a grid taken from an
         adaptive transient opens wherever that transient's window happened
         to start.  Here `h[0] = 1.4845` against a MEDIAN step of 4.62e-04
         -- 3200x coarser -- and the INNER Newton fails on that one step
         (100 iterations, residual 1e7x over tolerance at `v`).  The
         throwaway never met this because its unknown was `x_0` itself, so
         it had no manufacturing step to take.  `_period_grid` now opens on
         the grid's own finest step; see the guard there.

         ⚠ THE OLD DIAGNOSIS HAD ALREADY SEEN THE CAUSE AND FILED IT AS AN
         EXONERATION -- "the opening step is large not small" was written
         down as one of the obvious causes ELIMINATED.  It was checked
         against the wrong worry (that the opening step might be too small
         to open the trajectory) and never against its own size.

         ⚠ THE FIX PROPOSED HERE WAS TRIED AND IS WRONG AS STATED.  It
         read: `iq_{-1} = -(i(x_0) + u(t_0))` is exactly available from the
         DAE (item 4d), so trapezoidal can take a solved-history
         formulation seeded that way -- an exact Jacobian, clear of the
         fold.  The seeding half is true and was built.

         The formulation half is not.  A one-step companion depends on
         `x_{-1}` ONLY through `iq_{-1}`, and `d(iq_{-1})/d x_{-1} = -G` is
         SINGULAR at every purely reactive node -- most of a resonator.  So
         admitting `x_{-1}` as m unknowns leaves the 2m x 2m system
         rank-deficient, and it fails exactly as it should:
         `LinAlgError: Singular matrix`, on 25 tests at once.

         ⚠ THE CORRECTED DESIGN WAS ALSO WRONG, and was also built and
         measured before that was known.  It read: for a `b != 0` method
         the second unknown is `iq_0` itself -- the `(x, iq)` state its
         monodromy already uses -- closed by `iq_0 = iq_{N-1}`.

         It works, and then it does not.  On the Q=20 resonator it cut the
         outer solve to TWO residual evaluations (against euler's ten) and
         returned a peak identical to the plain path's 19.98968, confirming
         both that the prize is real and that trapezoidal's fixed point was
         never in doubt.  But `spectral_radius` came back 1.000000 where
         the circuit decays by exp(-pi/Q) = 0.854636 -- and that is the
         tell:

              npts   steps           outcome
               200     199 (odd)     converged
               201     200 (even)    LinAlgError
               202     201 (odd)     converged
               203     202 (even)    LinAlgError

         Trapezoidal's `iq` recursion is `iq_n = ... - iq_{n-1}`, whose
         homogeneous mode is `(-1)^n` -- UNDAMPED, exactly on the unit
         circle.  Over an even number of steps it returns `+1`, so `I - M`
         is singular and the solve dies; over an odd number it returns `-1`
         and nothing shows.  A user choosing 201 points instead of 200 gets
         a crash.

         ⚠ BOTH FAILED ATTEMPTS SHARE ONE CAUSE: trapezoidal's `iq` is not
         a coordinate a periodicity condition can pin.  Solving for a
         previous STATE is rank-deficient (`d(iq)/dx = -G` is singular at
         reactive nodes); solving for `iq` itself is degenerate (its mode is
         marginally stable, so closing it is vacuous or singular).

         SO THE PLAIN PATH IS CORRECT FOR TRAPEZOIDAL, not a legacy wart.
         Its Euler manufacturing step SUPPLIES `iq_0` from the DAE
         deterministically -- `-(i(x_0) + u(t_0))`, item 4d -- instead of
         asking a marginal mode to close.  That is why trapezoidal has no
         seam (1.3e-11 V) and converges where the enlargements do not.  The
         inexact Jacobian (item 4b, ~30%) is the price, and well-posedness
         is what it buys.

         A THIRD DESIGN WAS DERIVED AND ALSO FAILS, IDENTICALLY.  Make
         `iq_0` DEPENDENT rather than free -- `iq_0 := -(i(x_0) + u(t_0))`,
         unknown `x_0` alone, seeds `P_x(0) = I` and `P_q(0) = -G(x_0)`.
         It looked strictly better: m x m, no manufacturing step, exact
         Jacobian, and it needs only the FORWARD derivative `-G`, never the
         inverse that killed the first attempt.

         Measured, it dies on the same line:

              npts   steps           outcome        evals      rho
               200     199 (odd)     converged          2   1.000000
               201     200 (even)    LinAlgError        1        --
               202     201 (odd)     converged          2   1.000000
               203     202 (even)    LinAlgError        1        --

         Making `iq_0` dependent does not avoid the alternating mode -- it
         EXCITES it.  A perturbation `dx_0` gives `diq_0 = -G dx_0`, which
         drives the `(-1)^n` mode, which returns undamped with multiplier
         `(-1)^N`.  So the unit eigenvalue lands in `d x_end/d x_0` itself:
         even the m x m monodromy reads 1.000000, against the circuit's
         true 0.854636.

         ⚠ THREE DESIGNS, ONE OBSTRUCTION, AND IT NAMES WHAT THE PLAIN PATH
         IS FOR.  Every reformulation couples `iq_0` to `x_0` through `-G`,
         and that coupling excites a mode that never decays.  The plain
         path couples them through the EULER COMPANION instead -- its seed
         is `a_0 C`, not `-G` -- and that coupling is not degenerate: it
         reports 0.855 where all three of these report 1.000000.  The
         manufacturing step is not scaffolding to be removed; it is what
         makes trapezoidal's shooting problem well-posed, and the ~30%
         Jacobian error is what that costs.

         ⚠ THE "~30% JACOBIAN ERROR" IS WRONG IN BOTH DIRECTIONS, and the
         real statement is stronger.  Measured 2026-09-02 against a
         finite-difference `dF/dx_in` on three circuits at 100 points:

               circuit   method   relerr(J_code, J_true)   rank(J_true)/m
               RC        euler/trap        1.56                 1/3
               Q=20 RLC  euler             6.6e-03              2/4
               Q=20 RLC  trap              3.7e-03              2/4
               nonlinear euler/trap        4.73                 1/3

         So it is 0.4% on one circuit and 470% on another -- "~30%" was a
         single-circuit number carried as if it were a property of the
         method.

         ⚠ AND THE FRAMING WAS WRONG TOO: `J_code` IS NOT AN APPROXIMATION
         TO `dF/dx_in`.  That derivative is SINGULAR on every circuit above
         -- `sigma_min` is exactly 0.0 and the rank is 1/3, 2/4, 1/3 -- so
         AN EXACT NEWTON FOR THIS FORMULATION DOES NOT EXIST, and a solver
         handed the true Jacobian fails on the first step rather than
         converging faster.  What `_traverse` returns is `I - dx_end/dx_0`
         with a flat-history seed: an approximation to a DIFFERENT and
         WELL-POSED derivative, taken with respect to `x_0` while the
         unknown is `x_in`.  That is why the plain path works at all, and
         why its convergence is LINEAR rather than quadratic -- measured on
         the autonomous phase circuit, trapezoidal's residual falls
         3.91e-03 -> 3.14e-04 -> 2.66e-05 at a constant ratio ~0.076 where
         Gear-2's solved-history route is quadratic.  It is a contraction,
         not a Newton.

         ⚠ AND THE LITERATURE SETTLES THE DESIGN QUESTION: `x_0` IS THE
        CANONICAL UNKNOWN, so this is a RETURN and not an invention.
        Aprille & Trick (Proc IEEE 60(1) 108-114, 1972), "The Problem":
        "determine the periodic state w(0) such that integrating (1) from
        the initial state w(0) over the interval [0, T] we obtain the
        periodic solution", with Step 1 "for the given initial state x_0^i
        compute the solution x^i(t; x_0^i), 0 <= t <= T".  The trajectory
        begins AT the unknown; there is no pre-image and nothing
        manufactured before the first step.  Same in the oscillator paper
        (IEEE TCT 19(4) 354-360) and in Trick, Colon & Fan (TCAS 22(5)
        391-396) eq.(2).

        ⚠ AND THE IDENTITY SEED IS THEIRS TOO, which is what makes the frame
        error precise.  [AT-P] gives the discrete sensitivity for backward
        Euler as `Phi(T,0;x_0) = PROD [I - h F(x)]^-1` -- k factors seeded by
        `z(0)` at `t = 0`, i.e. the identity AT `x_0`, which in their
        formulation IS the unknown.  This file kept A&T's JACOBIAN and
        changed A&T's UNKNOWN.  `Px = [I, I]` is not wrong in isolation; it
        is correct for the formulation it came from and wrong for the one it
        now sits in.

        ⚠ WHAT PRECEDENT DOES NOT SETTLE is the COST, because A&T never paid
        it: they used BACKWARD EULER, which is L-stable, so the `(-1)^n`
        obstruction that forces an L-stable opener here never arose for
        them.  Their Phi is derived at fixed step for one method; the k-step
        and non-uniform-grid questions are later work.

        ⚠ WHICH MAKES THE FIX NAMEABLE: make `x_0` the unknown.  The
         throwaway driver that solved van der Pol did exactly that (its
         unknown was `x_0` itself), and item 5's note already half-records
         it.  Not built here; written down so the next attempt starts from
         the right statement rather than from "the Jacobian is 30% off".

         NOT PROPOSING A FOURTH.  Bordering the `(x, iq)` system to remove
         the alternating mode remains formally available -- the analogy to
         the phase condition is exact -- but the null direction there is a
         property of the discretisation with no closed form to pin, unlike
         `xdot(0)`, and three derivations in this item have now looked sound
         and failed on contact.  Anyone taking it up should measure before
         building: the falsifier is cheap and is even/odd step counts.

         ⚠ AND THE CONCLUSION IS NOW A THEOREM, NOT A TALLY OF THREE FAILED
         ATTEMPTS.  This item blamed trapezoidal's `iq` RECURSION for the
         `(-1)^n` mode.  That is the wrong cause, and an external review
         (2026-09-02) found the right one: trapezoidal is A-stable but NOT
         L-STABLE, so it maps the null space of the singular MNA `C` by
         exactly `-1` per step.  For `C x' + G x + u = 0` the one-step
         amplification is

             A_trap  = (C/h + G/2)^-1 (C/h - G/2)   ->  -I  on null(C)
             A_euler = (C/h + G)^-1 (C/h)           ->   0  on null(C)

         so the count of `-1` modes is exactly `m - rank(C)` -- one per
         ALGEBRAIC variable -- on every MNA circuit, and it has nothing to
         do with `iq`.  Verified here to 1e-09: the Q=20 resonator has
         `m = 4`, `rank(C) = 2` and exactly two eigenvalues at `-1` under
         trapezoidal and two at `0` under Euler; a plain RC has `m = 3`,
         `rank(C) = 1` and two of each.  A review reproduced the same
         singularity in an X-ONLY formulation containing no `iq` variable at
         all (cond 6.9e+18 at even K, 8.2e+04 at odd), which is the
         falsifier for the old attribution.

         What that buys, beyond a correct cause: the conclusion generalises.
         ANY formulation whose period map is `A_trap^K` without an L-stable
         opening step is singular at even `K` on every MNA circuit -- so
         "the plain path is correct for trapezoidal" is provable rather than
         observed.  And the cure is not specific to Euler: any L-STABLE
         opening step works, and it rescues exactly `m - rank(C)` modes.
         Corroborated in the literature the review checked: Houben (2003,
         App. A) biases theta off 1/2 so that "the numerical oscillations
         due to the DAE character of the equations are damped out during
         the 'insensitive time'" (p. 22, verbatim; an earlier rendering here
         elided three words without an ellipsis, caught by the docs
         session's quotation audit 2026-09-08) -- the same mechanism,
         attributed there to the DAE.

         Kept from the attempt: `_install_history` now takes the entering
         step size `h_prev` separately.  `x_{-1}` sits one step BEFORE
         `x_0`, and on a periodic grid that step is `hs[-1]`, not `hs[0]`.
         Uniform grids hide this; a caller's grid with a 16438:1 spread
         does not.

      5b. LTE-CHOSEN grid, the original wording -- pick the step sequence from an adaptive run and
         freeze it, refining BETWEEN shooting solves.  The grid still never
         moves inside one, so (3) stays exact.  Blocked on `Transient`
         accepting a non-uniform grid; `fixed_timestep` is uniform-only.
         Note this is the better structure anyway: a transient adapts
         because it cannot see the future, while PSS re-solves the same
         interval repeatedly and can therefore choose its grid once, well.
      6. MATRIX-FREE variational shooting (Telichevesky, Kundert & White,
         DAC 1995) -- the only structure that permits per-iteration
         adaptivity, because `M v` is obtained by integrating the
         variational system along the stored trajectory on that iteration's
         own grid.  The outer solve becomes an INEXACT Newton, phi shifting
         slightly between iterations.  A rewrite, not an increment on the
         above, and it should not leak into one.

         ITS CASE GOT STRONGER THREE TIMES ON 2026-09-01, from work that was
         not about it:

           - 4b/4c DOUBLE THE UNKNOWNS for a two-step method, so whatever
             the dense `J_phi` costs, this analysis now pays it on a vector
             twice as long.
           - THE LITERATURE NOTE puts a number on when that bites: Kundert
             (ICCAD'97) has forming and factoring `J_phi` at O(N^2 S) and
             O(N^3), "intractable when N exceeds several hundred", and names
             matrix-implicit Krylov as the answer -- this item.
           - 4d's SCALING ARGUMENT for method H pointed here instead.  H's
             only advantage over 4b was keeping the system at m unknowns for
             a dense solve; with a matrix-free solve the enlargement costs
             2x (vector length) rather than 8x (factorisation), so item 6
             removes H's reason to exist rather than competing with it.

         ⚠ AND THE COST TO ATTACK IS NOT THE FACTORISATION.  The final
         `J_phi` factorisation is one O(m^3).  The SENSITIVITY PROPAGATION
         that builds it is `N` steps of `_step_sensitivity`, and each of
         those is O(m^3) TWICE OVER: `linearsolver(Jf, S)` with a 2m-column
         right-hand side, and the `C @ P` products that form `S` and `Pq`
         -- m x m against m x 2m.  Matrix-free replaces the whole step, not
         half of it: it never forms `P`, so both go, leaving one matvec per
         Krylov iteration.

         MEASURED (2026-09-02, `benchmarks/pss_matrix_free_ceiling.py`), on
         a quiet box, single-threaded BLAS, every reading reproducible to
         better than 0.5%:

               m   propagation   (solve alone)   ceiling k=20
              40          4.8%            2.3%          1.02x
             110         15.0%            6.2%          1.14x
             242         38.1%           14.1%          1.54x
             502         63.8%           21.3%          2.58x
            1002         79.1%           24.9%          4.44x

         THE GATE IS PASSED.  It was "does the propagation share pass
         ~30%": it crosses near m=220 and reaches 79% at m=1002.  Item 6 is
         JUSTIFIED above m~250 and stays POINTLESS at m=40.

         AND IT IS NOW BUILT (2026-09-02) for the DRIVEN solved-history
         path: `solve(matrix_free=True)`, on `_traverse_factored` +
         `_monodromy_matvec` + `_matrix_free_solve`.  The recursion is not
         duplicated -- `_step_sensitivity` took a `solve` argument and is
         still the one recursion, run at width 1 against a stored
         factorisation instead of width 2m.  End to end: 1.36x at m=242,
         1.51x at m=502, 2.13x at m=1002, and the answers agree with the
         dense path to 1.1e-16.

         ⚠ k WAS ALSO A GUESS AND IS NOW MEASURED: GMRES takes 2/4/7/12
         iterations at m=40/110/242/502, because `I - M` has almost every
         eigenvalue within 1% of 1.0 -- the fast modes decay to nothing
         over a period, so k tracks the number of SLOW MODES, not m.  That
         is the property the item actually rests on, and it was assumed
         rather than checked until now.

         ⚠ WHAT IT COSTS: `2 N m^2` doubles of stored factorisations
         (~800 MB at m=1002, 50 points) where the dense path holds O(m^2);
         one MORE Newton iteration at m>=502 (3 against 2) from a
         convergence test that cannot be identical; and no monodromy, so
         `spectral_radius` is None on this path.  The autonomous and plain
         systems are NOT converted and raise rather than silently going
         dense.

         ⚠ THIS OVERTURNS THIS DOCSTRING'S OWN 2026-09-01 RECORD, and the
         box being quiet is NOT why.  That record read "the propagation is
         2.2% of a traversal at m=40, ceiling 1.01x-1.03x", and it is what
         happens when a measurement is named for the thing it was meant to
         settle rather than the thing it timed: it accumulated
         `toolkit.linearsolver` alone, which is under a third of the step
         at every size measured.  The 2.2% was reproducible, stable, and
         answering a different question.  ⚠ It also produced the confident
         wrong sentence "most of the time is ASSEMBLY -- which matrix-free
         does not touch"; at m=502 the propagation is 64% and assembly is
         not what dominates.  The m=40 VERDICT survives on the corrected
         number (4.8%, 1.02x) -- being right about m=40 is what made the
         error cheap to keep.

         ⚠ QUOTE THE THREADING CONDITION WITH THE NUMBER.  With BLAS
         threads free the same sizes read 4.8 / 7.2 / 12.4 / 18.2%: the
         `C @ P` products thread and the Python-level assembly does not.
         Single-threaded is the trustworthy column here -- the threaded
         traversal moved 2.64 s -> 4.02 s at m=502 between runs on this
         box while every single-threaded reading held to 0.1%.  A machine
         giving BLAS all its cores sees a smaller prize.

         ⚠ AND `k=20` IS AN ASSUMPTION, NOT A MEASUREMENT.  The ceiling
         charges matrix-free `k/m` of the propagation and nothing else, so
         it is an upper bound on an upper bound.  What the Krylov solve
         actually costs on these systems is unmeasured.

    Driving `Transient` -- done -- buys one integrator definition, the
    limiting/PCNR machinery, breakpoints and the order drop.  It does NOT
    buy (2) as a CONTROLLER; saying otherwise was this docstring's own
    earlier error, and item 4 above is what (2) turned out to be able to be
    instead.

    `reltol`, `iabstol` and `vabstol` mean exactly what they mean on
    `Transient` -- the tolerances of the TRANSIENT solution, applied to the
    per-timestep Newton, with the two absolute floors applied PER UNKNOWN in
    both flavours by `analysis.newton_tolerance_vectors`, the single
    definition all three analyses read.  Nothing here rescales them.

    The shooting criterion is expressed against that one: `steadyratio`
    (>= 1, default 1) multiplies it, so by default the shooting solve is
    held to the SAME relative tolerance as the transient, and raising it
    buys fewer shooting iterations for a looser periodic steady state.

    CHOOSING `method`: WHAT EACH ALTERNATIVE GIVES UP (owner decision,
    2026-09-08: the default is chosen for ACCURACY, not cost -- "we do not
    want to fool the user; instead they should change integrator and know
    the impact").  Every entry below is MEASURED in this tree
    (doc/pss_roadmap_260902.md, A10, the radau-default section and its
    monotonicity fixtures); nothing is quoted from a textbook order alone.

      method    order   period error   error        above its monotonicity limit
                        (ppm @ 400pt,  ESTIMABLE    (h_FE = 1/steepest slope)
                        Q=1e4 vdP)     by refinement
      radau     5 (6.1  5.7e-10        yes          OUTPUT stays in range at every
                on a                                step measured; only the STAGES
                smooth                              leave it, <= 0.3 % of the swing,
                orbit)                              above ~4 h_FE -- seen ONLY on a
                                                    scalar square-wave fixture (the
                                                    widest practical margin of the
                                                    four); on a driven detector
                                                    circuit the stages never leave
                                                    the hull at 20-320 pts/period,
                                                    with a diode OR a tanh device
      esdirk43  4       5.3e-05        yes          output in range; stages leave it
                                                    <= 3 % above ~2.5 h_FE
      trbdf2    2       10.1           yes          OUTPUT rings above 2.4 h_FE
      trap      2       20.8           NO -- its    OUTPUT rings above 2 h_FE
                                       error
                                       CHANGES SIGN
                                       near Q~100,
                                       so a two-grid
                                       estimate can
                                       under-state
                                       it 300x
      gear      2       83.1           yes          not an RK; a small ring at
                                                    h >= 10 tau measured; does not
                                                    CERTIFY a free-period solve at
                                                    1e-14 below ~400 pts at Q=1e4
      euler     1       --             yes          never rings; damps the orbit it
                                                    is asked to find (13 % of the
                                                    amplitude at 20 pts/period)

    Cost: radau is one real plus one complex factorisation per step (a
    coupled 3n system) against one per step for the others; at the point
    counts above it is cheaper on wall-clock anyway (60 points beat trap's
    480 on both axes, roadmap radau-default section), and `grid_error` /
    `warping_estimate` price the trade on YOUR circuit.  AS `n` GROWS
    (measured 2026-09-08, ladder oscillator, 200 points, same grid for
    every method): radau/gear on one transient period 3.0 / 3.5 / 3.9 /
    5.7 at n = 12 / 32 / 102 / 302, and on the whole PSS solve 3.0 / 3.9 /
    4.2 / 4.2 -- per-step assembly dominates to n ~ 100 and the 3n stage
    system only starts to show at 300; radau/trbdf2 1.5 -> 2.6; esdirk43
    is the MOST expensive method at every size (four sequential stage
    solves), never the cheap alternative.  Absolute: n = 302 at 200 points
    is 20 s per period and 271 s per PSS solve under radau, 64 s under
    gear.  At n = 1002 (transient period only, same grid): gear 26 s,
    trbdf2 55 s, radau 311 s -- radau/gear 11.8x, radau/trbdf2 5.7x -- so
    the 3n stage factorisation dominates from a few hundred unknowns and
    the ratio roughly doubles per 3.3x in n there; above ~300 unknowns
    the accuracy is bought at an order of magnitude in wall-clock, and
    trbdf2 (order 2, contractive, R = 1 + sqrt 2) is the alternative to
    price against `grid_error` on your circuit.  Index-2: the
    period keeps classical order under radau (6.1 measured on a smooth
    orbit); the algebraic unknowns converge at the stage order (3) -- and a
    relaxation oscillator whose PERIOD is timed by such a variable (a
    comparator on a current sense through a C-V loop) STILL keeps order 5
    (5.66 / 5.06 / 5.05 measured, identical errors to its voltage-sensed
    twin from 400 pts on), so no caveat applies.  ⚠ What that fixture DID
    show: its inner step Newton fails undamped at every grid (stage
    sensitivity tau/h, basin 0.94 h/(k tau)) -- hence the line search as
    the last resort -- and `warping_estimate` reads 0.09-0.65 of the true
    error on it: the septic interpolant does not resolve a comparator
    edge, and the estimate is only as good as the interpolant (Part I's
    "moderately smooth" boundary).  Roadmap, radau-default section.
    """

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name',
                   ## Sources supply their time-domain waveform only for an
                   ## analysis name in timedomain_analyses (('dc','tran')); the
                   ## old default 'PSS' matched nothing, so cir.u(t) returned 0
                   ## and the whole shooting solve had no excitation.
                   default='tran'),
         Parameter(name='reltol', 
                   desc='Relative tolerance', unit='', 
                   default=1e-4),
         Parameter(name='iabstol', 
                   desc='Absolute current error tolerance', unit='A', 
                   default=1e-12),
         Parameter(name='vabstol', 
                   desc='Absolute voltage error tolerance', unit='V', 
                   default=1e-12),
         Parameter(name='maxiter',
                   desc='Maximum number of iterations', unit='',
                   default=100),
         ## Forwarded to the inner Transient so PCNR (the junction-continuation
         ## limiting) reaches the shooting per-step solve too -- it lives in
         ## `Transient.solve_timestep`, which PSS DOES call, so no per-accepted-
         ## step machinery is needed (unlike breakpoints / continuation rescue,
         ## which are armed in `Transient.solve` and stay out of reach).
         Parameter(name='pcnr',
                   desc='Use Predictor/Corrector Newton-Raphson instead of '
                        'limiting in the inner transient; off by default',
                   unit='', default=False),
         ## ⚠⚠ THE ONE KNOB `method='theta'` HAS, AND IT WAS UNREACHABLE.
         ## `_integrator_for` builds `table[method]()`, so every shooting run
         ## took `ThetaIntegrator.DEFAULT_C` -- a RATE, calibrated on ONE
         ## fixture's period.  The transferable quantity is the DIMENSIONLESS
         ## `C T` (see `ThetaIntegrator.DEFAULT_CT` for why `h` cancels), so
         ## that is what this parameter is, and `_theta_biased` turns it into
         ## the rate THIS period needs.  `None` takes the measured knee.
         Parameter(name='theta_ct',
                   desc="method='theta' only: the null(C) damping one PERIOD "
                        'applies, as the dimensionless product C*T. None '
                        'takes ThetaIntegrator.DEFAULT_CT (the measured '
                        'knee). Ignored by every other method.',
                   unit='', default=None),
         ## `reltol` MEANS THE SAME THING IN EVERY ANALYSIS: the relative
         ## tolerance of the transient solution.  It is applied to the
         ## per-timestep Newton here exactly as `Transient` applies it, and
         ## nothing rescales it.
         ##
         ## `steadyratio` is how the SHOOTING criterion is expressed relative
         ## to it: shooting reltol = reltol * steadyratio, with 1 meaning the
         ## two are equal.  It is >= 1 because the period map is only KNOWN
         ## to the accuracy of the inner solves, so asking the shooting
         ## residual to beat that is asking it to resolve its own noise --
         ## refused rather than silently accepted.  Raise it to accept a
         ## looser periodic steady state for fewer shooting iterations.
         ## The LTE floors, separate from the Newton ones for the reason
         ## `Transient` records: one knob must not move both criteria.  Same
         ## names, same defaults, same meaning -- this analysis reports the
         ## number a transient would have controlled on.
         Parameter(name='lte_vabstol',
                   desc='Absolute voltage floor for the truncation-error '
                        'estimate', unit='V', default=1e-12),
         Parameter(name='lte_iabstol',
                   desc='Absolute current floor for the truncation-error '
                        'estimate', unit='A', default=1e-12),
         Parameter(name='TRTOL',
                   desc='Truncation error over-estimation factor (SPICE '
                        'TRTOL / lteratio in a commercial simulator)', unit='', default=7.0),
         Parameter(name='relref',
                   desc="What the relative LTE tolerance is measured "
                        "against: 'pointlocal', 'alllocal' or 'sigglobal'",
                   unit='', default='sigglobal'),
         Parameter(name='steadyratio',
                   desc='Shooting tolerance as a multiple of reltol (>= 1); '
                        '1 holds the shooting solve to the same relative '
                        'tolerance as the transient, larger relaxes it',
                   unit='', default=1.0),
         ## ⚠⚠ THE DEFAULT IS `radau` (owner decision, 2026-09-07), CHANGED
         ## FROM `trap`.  The floor of this stack is DISCRETISATION and it
         ## grows LINEARLY IN Q; at 240 points per period the relative error
         ## in the diffusion constant against the analytic high-Q reference is
         ##
         ##     Q      gear        trap        radau
         ##      100   1.79e-03    1.49e-06    6.97e-10
         ##      500   9.02e-03    1.04e-04    3.48e-09
         ##     1000   1.82e-02    2.36e-04    6.97e-09
         ##
         ## `trap` is not merely less accurate than `radau` here -- its error
         ## CHANGES SIGN near Q = 100, which is why it fits no clean law in
         ## that table and why `grid_error` has to refuse it (its two-grid
         ## difference under-states the true error by up to 300x there).  A
         ## default whose error estimate cannot be trusted is a poor default.
         ## `radau` is order 5, self-starting (no manufactured opener, so no
         ## seam in the period map), L-stable, and carries its own monodromy,
         ## so an autonomous run takes NO TR-BDF2 twin and reads its own
         ## spectrum -- see `monodromy_twin` and `carries_own_monodromy`.
         ##
         ## ⚠ It IS more expensive per step at a fine grid (3-stage fully
         ## implicit, through the 1-real/1-complex transform): PURE SOLVE time
         ## on van der Pol at Q=100, 480 points, is 3.538 s against trap's
         ## 2.586 s.  At a coarse grid it is cheaper (1.088 vs 1.519 at 120),
         ## the shooting Newton needing fewer iterations without an opener
         ## seam in the period map.
         ##
         ## ⚠⚠ BUT FOR ANY OSCILLATOR SURFACE THE TWIN DOMINATES, AND THAT IS
         ## WHAT SETTLES THE COST QUESTION.  `trap` is not self-sufficient: an
         ## autonomous run must solve a SECOND, TR-BDF2 PSS for its monodromy
         ## (`monodromy_twin`), and `radau` carries its own.  Measured
         ## `diffusion_constant` cost, which pays for that twin:
         ##
         ##     npts   trap solve / c-eval   radau solve / c-eval
         ##      120     1.519 / 1.255         1.088 / 0.187   (no twin)
         ##      480     2.586 / 4.161         3.538 / 0.737   (no twin)
         ##
         ## ⚠ AND AT EQUAL ACCURACY IT IS NOT CLOSE.  Relative error in `c`
         ## against the analytic high-Q reference, with total wall-clock:
         ##
         ##     Q~100   trap  480 pts  4.061e-06   6.747 s
         ##             radau  60 pts  7.599e-07   0.908 s
         ##     Q~500   trap  480 pts  9.230e-06  12.434 s
         ##             radau  60 pts  3.802e-06   0.633 s
         ##
         ## Radau at SIXTY points beats trap at four hundred and eighty, on
         ## both axes at once.  ⚠ Note also `trap` at Q~100 going 2.913e-06 at
         ## 240 to 4.061e-06 at 480 -- it does not even improve monotonically
         ## here, which is the sign change again.
         ##
         ## `trap` remains one argument away for a cheap coarse answer.
         Parameter(name='method',
                   desc="Integration method for the inner transient: 'radau' "
                        "(default, order 5), 'esdirk43', 'trbdf2', 'theta', "
                        "'gear' (BDF-2), 'trap' or 'euler'. The default is "
                        "chosen for ACCURACY; the class docstring's 'CHOOSING "
                        "method' table states, from measurement, what each "
                        "alternative gives up in order, error estimability, "
                        "monotonicity and cost",
                   unit='',
                   default="radau")]        

    
    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(PSS, self).parameters + self.parameters
        super(PSS, self).__init__(cir, **kvargs)
        ## The reference row is fixed for the analysis, and both the shooting
        ## loop and the Transient this drives need it.  It was recomputed in
        ## every method from a `refnode` argument that no caller ever varied.
        self.irefnode = self.cir.get_node_index(
            gnd if irefnode is None else irefnode)
        self._tran = None
        ## Only assembled when the period is an unknown; an extra assembly
        ## per timestep is not worth paying on the fixed-period path.
        self._want_dfdh = False
        self._dfdT = None
        ## ⚠ WHICH PERIOD-COLUMN CONVENTION `want_dT` USES.  'proportional'
        ## is the shipped one: every step scales with `T`, `dh_i/dT = h_i/T`.
        ## 'closing' is the commercial one Andreas described -- the inner
        ## transient owns the steps and the LAST one is placed on the period
        ## boundary, so `dh_i/dT = 0` inside and `dh_N/dT = 1` at the close.
        ## MEASURED (roadmap B7c gates 1 and 4): the proportional column is
        ## `O(h)` wrong on a smooth uniform grid and 4.2% RELATIVE wrong on
        ## van der Pol at `mu = 100` (step ratio 16438x), where 'closing' is
        ## 46x closer.  Default unchanged pending the rest of B7c.
        self._period_column = 'proportional'
        self._dfdh = None
        self._want_lte = False
        ## The caller's step fractions, or None for the uniform grid.  Read
        ## by the autonomous closures, which rebuild the grid at the current
        ## `T` on every residual evaluation.
        self._grid_fracs = None
        self._lte = None
        self._lte_seam = False
        self._lte_valid = True
        self._history_is_solved = False
        ## Set by `solve`: whether the entering history joined the unknowns.
        self.solved_history = False
        ## Reported by `solve`: the peak normalised truncation error over the
        ## converged period, and where in the period it fell.  None until a
        ## solve has run, or when the grid was too short to difference.
        self.max_lte = None
        self.max_lte_time = None
        self.max_lte_seam = None
        self.total_lte = None
        ## Set by `solve` when an autonomous run lands on a multiple of the
        ## fundamental; None when it did not (or on a driven run).
        self.fundamental_period = None
        ## What `factored_period()` needs to replay the CONVERGED period:
        ## which seed, which grid, which opening.  Written at the end of
        ## `solve`, and cleared at its start so a failed or interrupted run
        ## cannot leave a previous solution's state readable as this one's.
        self._period_state = None
        self._factored_period_cache = None
        ## ⚠ WHICH MONODROMY THE OSCILLATOR SURFACES READ -- see
        ## `monodromy_twin`.  Selects the method that supplies the PPV,
        ## Floquet modes and factored period when a one-step LMM (trap/euler)
        ## solved an autonomous circuit, whose OWN monodromy is first-order
        ## on a limit cycle (the opener seam).  'trbdf2' (DEFAULT): a TR-BDF2
        ## twin on the same grid -- self-starting, no opener, measured 12-32x
        ## more accurate on lambda2 than the Gear-2 twin at practical step
        ## counts (see `monodromy_twin`).  'gear': the former default, a
        ## Gear-2 twin, kept selectable.  'native': the run's OWN plain
        ## factorisation.  ⚠ 'native' UNDER A ONE-STEP METHOD IS THE WORST OF
        ## THE THREE, not an "exact, no twin" escape: trapezoidal's own
        ## monodromy is FIRST order on a limit cycle and its `Q` DIVERGES
        ## under refinement (11.1 / 28.4 / 63.9 vs an exact 5.9083), euler's
        ## likewise -- it is for the gates that measure that defect, not for
        ## results.  gear and trbdf2 runs are self-sufficient (second-order
        ## native monodromy) and ignore this -- they never twin.
        self.monodromy = 'trbdf2'
        self._monodromy_twin = None
        self._twins = {}
        self._solve_kwargs = {}
        ## Set to None by `solve` on the tstab path only -- see there.
        self.tstab_state = None

    def _is_autonomous(self, times):
        """True when nothing in the circuit depends on `t`.

        Exact where a spectral test is not: see `AUTONOMOUS_U_TOL`.  The
        source vector is evaluated at EVERY point of the grid and compared
        with the first; a circuit driven only by DC -- a VCO macromodel, a
        phase accumulator, an LC or ring oscillator -- has a constant `u`
        and a one-parameter family of periodic solutions, which is what
        makes fixed-period shooting ill-posed for it.
        """
        u0 = np.asarray(self.cir.u(times[0], analysis=self.par.analysis),
                        dtype=float)
        scale = max(float(np.max(np.abs(u0))), 1.0)
        ## ⚠ EVERY POINT ON THE GRID, NOT A STRIDE THROUGH IT.  This used to
        ## sample `times[1::len(times)//8]` -- about nine points -- while
        ## calling itself "exact where a spectral test is not".  Nine points
        ## cannot see a narrow pulse: measured on an RC driven by a `VPulse`
        ## positioned BETWEEN two samples, 40% and 20% duty were read
        ## correctly and 5%, 1% and 0.5% all came back AUTONOMOUS.  A clock
        ## misread that way is routed to the free-period system, which
        ## solves for `T` and DISCARDS the period the caller asked for --
        ## and `DEGENERATE_PERIOD_FACTOR` cannot catch it, because it tests
        ## the magnitude of `T`, not whether the circuit was driven.  PWM,
        ## sampling clocks, S/H and mixer LOs are core PSS workload and are
        ## exactly the shapes a stride misses.
        ##
        ## The cost is `N` evaluations of `u` ONCE per solve, against `N`
        ## Newton solves in the traversal it decides -- and the loop exits
        ## at the first sample that differs, which is the common case for
        ## every driven circuit.
        for t in times[1:]:
            u = np.asarray(self.cir.u(t, analysis=self.par.analysis),
                           dtype=float)
            if float(np.max(np.abs(u - u0))) > AUTONOMOUS_U_TOL * scale:
                return False
        return True

    ## The one mapping from `method` to a class.  Read by `_transient` to
    ## build the integrator and by `_companion_reach` to ask how far it
    ## reaches, so the two can never disagree about which method is running.
    @classmethod
    def _integrator_for(cls, method):
        from pycircuit.circuit.integrator import (EulerIntegrator,
                                                  TrapezoidalIntegrator,
                                                  ThetaIntegrator,
                                                  Gear2Integrator)
        from pycircuit.circuit.integrator import (TRBDF2Integrator,
                                                  RadauIIA3Integrator,
                                                  ESDIRK43Integrator)
        ## THE single method -> integrator map, and the one place method names
        ## are validated: an unknown name raises the ValueError here rather than
        ## a KeyError three frames down.  The polymorphic predicates
        ## (`_companion_reach`, `needs_x0_unknown`, ...) call this before the
        ## solve()-level whitelist runs, so the validation must live where the
        ## lookup does.
        table = {'euler': EulerIntegrator,
                 'trap': TrapezoidalIntegrator,
                 'trapezoidal': TrapezoidalIntegrator,
                 ## `theta` is trapezoidal biased by `C h` -- see
                 ## `ThetaIntegrator`. It takes the L-stable opener OUT, which
                 ## is the one thing the other one-step LMMs cannot do.
                 'theta': ThetaIntegrator,
                 'gear': Gear2Integrator,
                 'gear2': Gear2Integrator,
                 'trbdf2': TRBDF2Integrator,
                 'radau': RadauIIA3Integrator,
                 'esdirk43': ESDIRK43Integrator}
        try:
            return table[method]()
        except KeyError:
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau' or 'esdirk43', not %r" % (method,))

    ## Below this fraction of the seed, a solved period is the trivial
    ## root rather than an orbit.  Deliberately loose: a real fundamental
    ## reached from a seed one decade high is ~0.1 of it, and the trivial
    ## root lands 15 orders down, so nothing sits near this line.
    DEGENERATE_PERIOD_FACTOR = 1e-6
    ## The equilibrium test above: a returned state whose DC residual is
    ## within this factor of `iabstol` is an equilibrium, not an orbit.  1e3
    ## because the shooting Newton's own residual sits at `abstol`, and an
    ## orbit's DC residual at t = 0 is a CIRCUIT-scale current (measured:
    ## the van der Pol at 2 V has |i + u| ~ 1 A there against 1e-27 for the
    ## collapsed state -- 27 orders apart, so the factor is not delicate).
    TRIVIAL_ORBIT_FACTOR = 1e3

    def _free_period_solve(self, func, z0, abstol, xtol, reltol, maxiter,
                           seed_period, solver=None):
        """Solve a free-period system, with its degenerate root named.

        ⚠ `T = 0` IS A REGULAR ROOT OF EVERY AUTONOMOUS SHOOTING SYSTEM.
        `x0 - phi_T(x0)` vanishes identically at `T = 0`, and the phase
        condition does not exclude it -- it constrains `x0`, not the
        period -- so Newton reaches it from any seed below the fundamental
        and the run returns a period of ~1e-18 with no orbit in it.

        Measured on BOTH autonomous elements in the tree, so it is a
        property of the formulation and not of any circuit: from a 1e-4
        seed against a 1e-3 fundamental, Gear-2 returned -1.5e-20 on the
        quadrature element and 3.9e-19 on the scalar `Idtmod`, and
        trapezoidal raised a bare `LinAlgError` from three seeds of five as
        its Jacobian went singular on the way down.

        Neither outcome is a silent wrong answer -- the collapse reports
        `converged = False` (⚠ ENFORCED HERE, by demoting `ier`; asserting it
        in prose was not enough -- see the note at the demotion) and the
        exception is loud -- but both told the
        user nothing about the cause, and the generic non-convergence
        advice ("raise maxiterations") is wrong for it: no number of
        iterations reaches a fundamental from below.
        """
        try:
            if solver is None:
                z, info, ier, mesg = analysis.fsolve(
                    func, z0, maxiter=maxiter, reltol=reltol, abstol=abstol,
                    xtol=xtol, toolkit=self.toolkit, full_output=True,
                    line_search=True)
            else:
                ## ⚠ THE MATRIX-FREE ROUTE COMES THROUGH HERE TOO, so the
                ## trivial-root diagnosis below covers it.  Routing it around
                ## this wrapper would have lost the one message that makes an
                ## autonomous collapse readable.
                z, info, ier, mesg = solver(z0, abstol, xtol, reltol, maxiter)
        except np.linalg.LinAlgError as exc:
            raise np.linalg.LinAlgError(
                'PSS: the free-period Jacobian went singular while solving '
                'for an autonomous period seeded at %.6g s (%s). The usual '
                'cause is a seed BELOW the fundamental: `T = 0` solves the '
                'periodicity condition identically, so the iteration is '
                'drawn to it and the Jacobian degenerates on the way. Seed '
                'at or above the expected period -- a short transient and '
                'the interval between two output recurrences is the usual '
                'way to get one.' % (seed_period, exc)) from exc

        T = float(z[-1])
        ## ⚠⚠ THE SECOND TRIVIAL ROOT, found 2026-09-08 while gating
        ## `warping_estimate` on an index-2 oscillator: the guard below
        ## catches `T -> 0`, and an autonomous solve has ANOTHER root that it
        ## cannot see -- the EQUILIBRIUM, `x(t) = x_dc`, which is periodic at
        ## EVERY `T`.  Radau seeded 10 % below the fundamental on the
        ## index-1 van der Pol and on the index-2 fixture returned
        ## amplitude 0.0000 (state 1e-27) at a period near the seed with
        ## `converged = True`; trapezoidal on the same seed failed honestly.
        ## `T` stays finite, so the period test passes, and the periodicity
        ## residual is exactly zero because the equilibrium IS periodic --
        ## at every `T`, so it is a whole LINE of roots in `(x0, T)`, not a
        ## point (peer's sharpening): that is why the residual is exactly
        ## zero rather than small, why no residual-based guard could have
        ## caught it, and why the DC residual does in one evaluation.
        ## The test that sees it is the DC residual of the returned state:
        ## an orbit has `C x' != 0` somewhere at t = 0, so `i(x) + u` is far
        ## from zero there; an equilibrium has it at solver tolerance.
        trivial_orbit = False
        try:
            xr = np.asarray(z[:-1], dtype=float)
            irn = self.irefnode
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            r_dc = (np.asarray(self.cir.i(xf, self.epar), dtype=float).ravel()
                    + np.asarray(self.cir.u(0.0, self.epar, analysis='dc'),
                                 dtype=float).ravel())
            r_dc = np.delete(r_dc, irn)
            tol = float(getattr(self.par, 'iabstol', 1e-12))
            trivial_orbit = bool(np.abs(r_dc).max() <= self.TRIVIAL_ORBIT_FACTOR * tol)
        except Exception:
            trivial_orbit = False
        if trivial_orbit:
            ier = 5
            mesg = ('collapsed onto the EQUILIBRIUM (a trivial orbit, periodic '
                    'at every T) at T = %.6g s from a seed of %.6g s' % (T, seed_period))
            warnings.warn(
                'PSS: this autonomous solve returned an EQUILIBRIUM, not an '
                'orbit: the state at t = 0 satisfies the DC equations to '
                '%.1e (max |i(x) + u|), so the periodicity residual is zero '
                'at ANY period and the solver reported success at T = %.6g s '
                'from a seed of %.6g s. The basin of the equilibrium is '
                'entered from a seed below the fundamental (measured: 10 %% '
                'low under radau); seed at or above the expected period, or '
                'from a transient that is already on the orbit. '
                '`converged` is False.' % (np.abs(r_dc).max(), T, seed_period),
                RuntimeWarning, stacklevel=3)
        elif not np.isfinite(T) or abs(T) < self.DEGENERATE_PERIOD_FACTOR * abs(
                seed_period):
            ## ⚠⚠ THE COLLAPSE MUST BE DEMOTED HERE, and for two turns of this
            ## record it was not.  The docstrings above and on `solve` both
            ## asserted "the collapse reports `converged = False`" -- and
            ## NOTHING ENFORCED IT.  `self.converged` is `(_ier == 1)` and
            ## nothing else, while `T = 0` is a REGULAR root: the solver
            ## reaches it cleanly and reports success, so Gear-2 returned
            ## `T = 5.42e-18` with `converged = True` on a circuit with no
            ## orbit in it.  The warning fired correctly the whole time, which
            ## is exactly what made this survive -- a reader who checks the
            ## documented flag instead of catching warnings got `True`.
            ##
            ## Demoting `ier` rather than assigning `self.converged` is
            ## deliberate: all three autonomous call sites already feed this
            ## return value into `self.converged`, so one demotion covers the
            ## plain, solved-history and matrix-free paths, and any future
            ## path inherits it by construction.  `ier = 5` is `fsolve`'s
            ## "not making good progress" code -- the closest existing
            ## meaning, and already handled everywhere `ier` is read.
            ier = 5
            mesg = ('collapsed onto the trivial root T = %.6g s from a seed '
                    'of %.6g s' % (T, seed_period))
            warnings.warn(
                'PSS: this autonomous solve collapsed onto the TRIVIAL root, '
                'returning a period of %.6g s from a seed of %.6g s. `T = 0` '
                'satisfies `x0 - phi_T(x0) = 0` identically and the phase '
                'condition does not exclude it (it constrains x0, not the '
                'period), so a seed below the fundamental is drawn there. '
                'The returned waveform is not a periodic steady state. '
                'Raising maxiterations will not help -- seed at or above the '
                'expected period instead; a short transient and the interval '
                'between two output recurrences gives one. The literature '
                'remedy for this basin is the PROBE technique -- a periodic '
                'voltage source that feeds the oscillator until its own '
                'current reaches zero, widening the basin so a Newton is '
                'less likely to fall into the DC solution (Bizzarri et al., '
                '"Probe Based Shooting Method ..."); it is not implemented '
                'here, and it widens the basin rather than removing the '
                'seed dependence.'
                % (T, seed_period), RuntimeWarning, stacklevel=3)
        return z, info, ier, mesg

    def _resolve_break_events(self, requested):
        """`break_events`, defaulted from the METHOD when not given.

        ⚠⚠ THIS IS ON FOR ONE-STEP METHODS AND OFF FOR A MULTISTEP ONE, AND
        THE SPLIT IS MEASURED RATHER THAN ASSUMED.  Landing a source's
        discontinuities on grid points helps a one-step method and HURTS
        Gear-2, on the same circuit, at the same step count:

            method               uniform     + events   jittered, no events
            gear   (multistep)   8.23e-03    1.29e-02   1.24e-02   lost 7 of 9
            trap   (one-step)    4.98e-03    3.15e-03   6.82e-03   lost 0 of 9
            radau  (one-step)    --          1.02-1.89x gain       lost 0 of 9

        ⚠ **The jittered column is the control that makes this a cause.**  A
        grid of the same step COUNT and comparable non-uniformity, with the
        events deliberately NOT landed, hurts gear just as much as the event
        grid does (1.24e-2 against 1.29e-2).  So gear's loss is NON-UNIFORMITY
        ITSELF, not a defect in `event_grid`: a multistep method's companion
        coefficients depend on the step-size RATIO, so a uniform grid is its
        best case and any insertion is a real cost.  `trap` pays that cost too
        (jittered 6.82e-3 against uniform 4.98e-3) and the event alignment is
        worth MORE than the cost, so it nets out ahead.

        The predicate is `companion_reach() == 1` -- the method's own statement
        of how many charges back its companion reads, which is exactly the
        property that makes step ratios matter.  `RungeKuttaIntegrator` says
        the mechanism in its own words: *"a one-step method carries no
        zero-stability step-ratio limit"*.  Asked of the method, never inferred
        from a name.

        ⚠ An explicit `True`/`False` is honoured untouched; this only fills in
        `None`.
        """
        if requested is not None:
            return bool(requested)
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        return int(integ.companion_reach()) == 1

    def _resolve_x0_unknown(self, requested):
        """`x0_unknown`, defaulted from the circuit's TOPOLOGY when not given.

        ⚠⚠ WHY THIS IS CONDITIONAL AND NOT A NEW GLOBAL DEFAULT. `x0_unknown`
        is NOT free: trapezoidal still needs an L-stable opener, so switching
        it on moves the Euler step INSIDE the period, where it degrades the
        ORBIT rather than just the opening. Measured on a `Q = 20` resonator
        against its analytic 20 V peak, `x0_unknown` is WORSE --
        20.01273 against 19.76939 at 100 points, 20.02208 against 19.96123 at
        200. Turning it on everywhere would trade a real defect on a few
        circuits for a real regression on most.

        ⚠ ON AN INDEX-2 NETLIST THE TRADE REVERSES, and not marginally. The
        manufactured opening step is INCONSISTENT there: the constraint fixes
        the algebraic variable at a value the step cannot produce, so
        trapezoidal returns EXACTLY 2x on an L-I cutset -- and on an even
        number of steps reports `converged` and a periodicity residual of
        1e-13 while doing it. See the roadmap's section 0k. `x0_unknown`
        removes it at every parity because `x(0)` becomes a genuine unknown.

        Three refusals, each deliberate:

          * **an explicit `True`/`False` is honoured untouched** -- this only
            fills in `None`;
          * **a two-step method is left alone**, because its solved-history
            formulation already solves for `x(0)` and was never affected;
          * ⚠⚠ **a PROVISIONAL verdict does not trigger it, and that is a
            REFUSAL ON THE THEORY RATHER THAN CAUTION.** With a controlled
            source in the loop or cutset the index is not bounded by 2 and
            need not be a function of the topology at all (see
            `topological_index`).  So the premise "the criterion PROVES index
            2, therefore switch" is unavailable — and so is the REMEDY's
            justification, because `x0_unknown` fixes an inconsistent opening
            step on an INDEX-2 algebraic row.  **If the true index is 3 the
            remedy is not known to apply, and switching it on would mask a
            worse problem while reporting a fix.**  Leaving a known defect
            visible is the better failure mode.  The same goes for a
            structurally singular netlist, which has no index at all.

            ⚠ MEASURED 2026-09-08 (docs session, two instruments sharing
            only C and G -- InitDAE eq. 8's 1-fullness of the derivative
            array, and the Kronecker index of the pencil -- 5/5 against
            `topological_index` where that is valid, exact on nilpotent
            pencils of degree 1..5): **index 3 IS reachable on this element
            set, it is VALUE-dependent, and on every provisional fixture
            `topological_index` returned 1 -- not a low-confidence 2.**  A
            VCVS of gain `g` inside a C-V loop (C1 v->b, C2 v->gnd, source
            b->gnd = g*v) is index 3 exactly on `g* = 1 + C2/C1` (six
            (C1, C2) pairs verified) -- the controlled source cancels the
            node's total capacitance, `[C2 + C1(1-g)] dv/dt`, and one more
            differentiation is needed -- and index 2 off it, with
            `|M^3| = 224 |g - 2|`, so a gain within 1e-3 of `g*` still
            carries a ~1 % nilpotent tail: the surface is measure-zero, the
            NEIGHBOURHOOD is what bites.  ⚠ A CCCS half of the same report
            was WITHDRAWN by its author the same day: its fixture put the
            CCCS's ammeter input ACROSS the inductor instead of in series
            (DC-singular at every gain), and rebuilt correctly the L-I
            cutset is genuine, `topological_index` says 2 (provisional), the
            index is 2 at every gain but `F = 1`, and at `F = 1` the row at
            the output node is all-zero so `_structural_singularity` fires
            and `GminAnchorNewton` refuses the rescue in as many words --
            the DC layer catches it loudly; nothing for `ill_posed` to do.
            So: on the VCVS netlists (P1, its self-controlled and
            CCVS-controlled variants) the topological number is 1, wrong by
            two, and `idx != 2` short-circuits before `provisional` is
            consulted -- two independent guards, both firing; and the
            emphasis is INVERTED from "the case that fails to solve": **the
            index-3 netlist DC-solves cleanly and silently at `g*`** (v =
            1e-3 A x 10 Ohm, b = g v, no warning), while the singular one
            is loud.  The case to guard against is the one that solves and
            looks healthy.  What is NOT established: that any of this is
            reachable through the analog blocks people write; the claim is
            about what the element set PERMITS, which is what a refusal has
            to be justified against.  The numerical test (one SVD of an n(k+1) matrix) is
            an OFFER, not built: exact for constant C and G only, so on a
            nonlinear netlist it is the index of the linearisation at one
            operating point, and controlled sources move the index with the
            operating point too.  Code and fixtures:
            `~/docs/.corpus/checks/numindex.py`, `dae_index_probe.py`
            (`--quick` asserts the surface and both pencil cases).

        Warns when it fires, because a silently different formulation is the
        kind of thing that makes a later measurement inexplicable.
        """
        if requested is not None:
            return bool(requested)
        ## ⚠ EVERYTHING BELOW IS BEST-EFFORT AND MUST NEVER RAISE.  This runs
        ## BEFORE `solve` validates its own arguments, so a bad `method` was
        ## reaching `_solves_history` and coming back as `KeyError: 'bogus'`
        ## instead of the `ValueError('method must be ...')` the caller is
        ## owed -- two tests caught exactly that.  A defaulting helper has no
        ## business changing which exception an invalid call raises.
        try:
            if self._solves_history():
                return False
            idx, info = topological_index(self.cir)
        except Exception:
            return False
        if idx != 2 or info['provisional'] or info['ill_posed']:
            return False
        where = (('C-V loop: ' + ', '.join(info['loop'])) if info['loop']
                 else ('L-I cutset: ' + ', '.join(info['cutset'])))
        warnings.warn(
            'PSS: this netlist is index 2 (%s), where the manufactured '
            'opening step is INCONSISTENT -- it starts an algebraic variable '
            'at a value the constraint forbids, and trapezoidal carries that '
            'seed forever (exactly 2x on an L-I cutset, reported as CONVERGED '
            'on an even number of steps). Solving for x_0 directly instead; '
            'pass x0_unknown=False to override.' % where,
            RuntimeWarning, stacklevel=3)
        return True

    def event_grid(self, period, npts=None, grid=None, min_sep=0.05):
        """A step grid with the circuit's EVENT TIMES landed on exactly.

        A6/B7b.  `Transient` breaks its steps at `cir.next_event`; the PSS
        traversal does not, so a pulse edge inside a step is integrated straight
        through.  This returns step FRACTIONS for `solve(grid=...)` with each
        event in the period placed ON a grid point -- by SNAPPING the nearest
        point onto it when one is close, and INSERTING otherwise, so no
        arbitrarily small step is ever created (the B7c lesson).

        Measured on an RC driven by a `VPulse`, against a 4000-point reference,
        a 40-step uniform grid versus the same grid with its 3-4 event times
        landed on::

            edge offset   uniform      + events        gain
            td = 0        6.787e-03    8.032e-04       8.5x
            td = 0.0125T  5.720e-03    2.099e-04        27x
            td = 0.0092T  2.483e-03    7.910e-04       3.1x

        ⚠ TIME-DRIVEN EVENTS ONLY, AND THE LIMIT IS STRUCTURAL.  `next_event(t)`
        is parameterised by time, so a source's edges can be walked out once and
        placed.  A STATE-DEPENDENT reset -- `Idtmod`'s wrap -- cannot: its
        `next_event` is a linear prediction from the last accepted point and
        returns `inf` before a traversal has started, so there is nothing to walk.
        ⚠ THE REASON RECORDED HERE WAS WRONG, AND IS CORRECTED (2026-09-06).
        This used to say the wrap time "has to become an unknown the Newton
        solves for", citing a map "discontinuous by |dphi| ~ 8.2e-3 on a grid
        point".  Measured: that jump is **grid-INDEPENDENT** -- 1.414214e+09
        at seven different grids, with the wrap ON a node and OFF it alike,
        and unchanged when the exact wrap times are added to the grid.  A
        quantity that does not move when the grid moves is not a grid
        artefact.  The real defect was the fold at the period ENDPOINT, in the
        OUTPUT map, and it is fixed in the RESIDUAL by `_fold_periodic`.
        The conclusion for THIS method is unchanged and now for the right
        reason: a state-dependent reset is not a grid feature, so `event_grid`
        does not help it -- but nor does it need the traversal surgery that
        sentence implied.

        ⚠ AND THIS IS NOT SALTATION.  Saltation was measured and falsified twice
        for this codebase (a switched conductance, a discontinuous injection and
        an `Idtmod` wrap all give a monodromy-vs-FD gap falling at 2.00x per
        doubling, i.e. O(h)); each step already uses its own `Jf`/`C`, which
        describe whichever side of the switch that step is on.  The problem was
        only ever that the grid could not BREAK at the event.
        """
        T = float(period)
        if not T > 0.0:
            raise ValueError('event_grid: period must be positive, got %g' % T)
        if grid is not None:
            fr = np.asarray(grid, dtype=float).ravel()
            tot = float(np.sum(fr))
            if not np.isclose(tot, 1.0, rtol=0, atol=1e-9):
                raise ValueError('event_grid: `grid` fractions must sum to 1, '
                                 'they sum to %.12g' % tot)
            pts = np.concatenate(([0.0], np.cumsum(fr)))
            pts[-1] = 1.0
        else:
            if npts is None:
                raise ValueError('event_grid: give either `npts` or `grid`')
            pts = np.linspace(0.0, 1.0, int(npts) + 1)

        ## walk the events across one period
        ev = []
        t = 0.0
        for _ in range(10 * len(pts) + 100):
            e = float(self.cir.next_event(t))
            if not np.isfinite(e) or e >= T * (1.0 - 1e-15):
                break
            if e > T * 1e-15:
                ev.append(e / T)
            if e <= t:
                break
            t = e
        self.event_times = list(ev)
        if not ev:
            return list(np.diff(pts))

        for f in ev:
            j = int(np.argmin(np.abs(pts - f)))
            if j == 0 or j == len(pts) - 1:
                ## never move an endpoint: the period boundary is not ours
                k = 1 if j == 0 else len(pts) - 2
                h = abs(pts[k] - pts[j])
                if abs(pts[k] - f) < min_sep * h:
                    pts[k] = f
                    continue
            else:
                h = min(pts[j] - pts[j - 1], pts[j + 1] - pts[j])
                if abs(pts[j] - f) < min_sep * h:
                    pts[j] = f          ## SNAP -- no tiny step created
                    continue
            pts = np.append(pts, f)
            pts = np.sort(pts)
        pts = np.unique(pts)
        return list(np.diff(pts))

    def refine_grid(self, grid, x0, period=None, refnode=gnd, gamma=2.0,
                    delta=0.25, reltol=None, timestep=None):
        """REPAIR an under-resolved step grid, given a solution solved on it.

        B7c.  Returns a new fraction list: `grid` with points ADDED wherever an
        adaptive run from `x0` asks for a step more than `gamma` times finer
        than what is there.  Feed it back to `solve(grid=..., x0=...)`.

        ⚠⚠ THIS IS A REPAIR PATH, NOT A REPLACEMENT FOR `lte_grid`.  Measured on
        van der Pol at mu=100: refining a deliberately decimated grid recovers
        **+424.8 ppm -> +18.2 ppm in ONE pass** for 2.5x the points, and then
        REACHES A FIXED POINT (1407 -> 1411 -> 1415, +4 points a stage, error
        unchanged).  But handed a grid that is already good it has nothing to do:
        a `lte_grid`-quality grid went 1137 pts / -3.81 ppm -> 1153 pts /
        -3.98 ppm, adding ~16 points and drifting marginally WORSE.  Use it when
        you suspect a grid is too coarse, not as a matter of course.

        ⚠ **SOLVE FIRST, THEN REFINE -- THE ORDER IS THE WHOLE DESIGN.**  Doing
        this from scratch, refining at every shooting iteration, costs **2.7x -
        5.2x** the points for no accuracy gain, and a warmup does NOT fix it:
        with iterates that are shrinking perturbations (3e-2 -> 0) of the settled
        point the grids match in SIZE (1238, 1209, 1185, 1169, 1137) but their
        points barely coincide, because a tiny perturbation of `x_0` shifts every
        step boundary.  Adaptive step POSITIONS are not stable under small state
        perturbations.  Once the solve has converged the iterates stop moving,
        the criterion stops firing, and the grid settles -- which is why this
        takes an `x0` that is already a solution.

        ⚠ `delta` IS THE SEPARATION RADIUS, SCALED TO THE STEP THE CONTROLLER
        ASKED FOR -- not to the gap that happens to be there.  A new point is
        refused if it would sit within `delta` times the WANTED local step of an
        existing one; without that, merged grids acquire arbitrarily small steps.
        ⚠⚠ AND ITS SETTING DEPENDS ON WHICH RULE IS USING IT -- a constant
        transplanted between the two is a silent no-op.  For the UNION scheme
        (merge whole per-iterate grids) the sweep put the knee at `delta = 1`:
        0 -> 1 takes the cost 5.21x -> 2.74x at no accuracy cost, and at 1.5 the
        count falls further (1.50x) while the error collapses forty-fold
        (+339.6 ppm).  **For the SUBDIVISION rule used here `delta = 1` rejects
        almost everything** -- the points it inserts are already spaced about one
        WANTED step apart, so demanding a full step of clearance from the
        interval ends refuses them.  Measured: `delta = 1` here recovered
        +424.8 -> +423.5 ppm, i.e. nothing, adding 4 points a stage.  `0.25` is
        the measured value for THIS rule and is the default.

        ⚠ Scaling the radius to the EXISTING grid's local gap instead of the
        wanted step was also tried and admitted 33 of 1158 points on a coarse
        grid, producing grids that did not converge at all.

        `gamma` is the weaker knob (1.5 -> 5.0 moves the cost only 3.24x ->
        2.79x and starts costing accuracy at 5); 2 is the measured default.

        All of the above is ONE fixture (van der Pol mu=100, `gear`), one seed.
        """
        import warnings as _warnings
        from pycircuit.circuit.transient import Transient
        fr = np.asarray(grid, dtype=float).ravel()
        if fr.ndim != 1 or fr.size < 1:
            raise ValueError('refine_grid: `grid` must be a list of step '
                             'FRACTIONS, got %r' % (grid,))
        tot = float(np.sum(fr))
        if not np.isclose(tot, 1.0, rtol=0, atol=1e-9):
            raise ValueError(
                'refine_grid: the step fractions must sum to 1 (they are '
                'fractions of the period, as `solve(grid=...)` takes); they '
                'sum to %.12g. Pass `np.diff(points)`, not the points.' % tot)
        T = float(self.par.period if period is None else period)
        cur = np.concatenate(([0.0], np.cumsum(fr)))
        cur[-1] = 1.0

        ## the grid an adaptive run WANTS from this solution, as fractions
        xf = np.asarray(x0, dtype=float).ravel()
        if xf.shape[0] == self.cir.n - 1:
            xf = self._insert_refnode(xf)
        rt = self.par.reltol if reltol is None else float(reltol)
        h0 = (T / 200.0) if timestep is None else float(timestep)
        tr = Transient(self.cir, toolkit=self.toolkit, reltol=rt)
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            res = tr.solve(refnode=refnode, tend=T, timestep=h0, x0=xf)
        t = np.asarray(res.sweep_values, dtype=float).ravel()
        if len(t) < 3:
            raise RuntimeError(
                'refine_grid: the adaptive run put only %d points in the '
                'period, which is not a grid to refine against.' % len(t))
        cand = np.clip((t - t[0]) / (t[-1] - t[0]), 0.0, 1.0)
        hc = np.diff(cand)
        hloc = np.minimum(np.r_[hc[0], hc], np.r_[hc, hc[-1]])

        ## subdivide only the intervals that are coarser than asked
        want = np.interp(0.5 * (cur[:-1] + cur[1:]), cand, hloc)
        have = np.diff(cur)
        out = [cur[0]]
        for i, (lo, hi) in enumerate(zip(cur[:-1], cur[1:])):
            if have[i] > gamma * want[i] and want[i] > 0.0:
                k = int(np.ceil(have[i] / want[i]))
                for q in lo + (hi - lo) * np.arange(1, k) / k:
                    ## the separation rule, on the WANTED step
                    if min(q - lo, hi - q) >= delta * want[i]:
                        out.append(q)
            out.append(hi)
        pts = np.unique(np.asarray(out, dtype=float))
        return list(np.diff(pts))

    def lte_grid(self, period, x0=None, refnode=gnd, tstab=None,
                 reltol=None, timestep=None):
        """Step FRACTIONS for `solve(grid=...)`, derived from an adaptive run.

        B7a.  A transient adapts because it cannot see the future; PSS
        re-solves the SAME interval over and over, so it can be handed a
        grid that was chosen well ONCE and then frozen.  This is the
        derivation side; `_period_grid` is the consumption side and has
        been shipped since item 5.

        Returns `(fracs, seed)` -- the accepted steps of one settled
        period as fractions summing to 1, and the state at the start of
        that window, ready to pass straight back in:

            fracs, seed = pss.lte_grid(period=T)
            pss.solve(period=T, grid=fracs, x0=seed)

        MEASURED on van der Pol at `mu = 100`
        (`benchmarks/pss_lte_grid.py`, the gate this was promoted from):
        1105 derived steps converge where 1105 UNIFORM steps do not, and
        beat a 20000-point uniform grid -- 18x fewer points and -47.3 ppm
        against -60.6.

        ⚠⚠ IT IS FOR STIFF SMOOTH PROBLEMS AND NOT FOR EVENTS, and that
        boundary is measured rather than cautionary.  On a wrapping
        `Idtmod` the derived grid is WORSE than a uniform grid of the same
        count -- max LTE 2.64e+05 against 1.67e+05 times tolerance at
        ~1429 steps -- because the LTE peak sits at the RESET on every
        grid, and no step size makes a discontinuity's local truncation
        error small.  The event half of B7 is a different item: it needs
        the event time to be an unknown the Newton solves for, because a
        grid frozen from a PAST traversal cannot represent an event whose
        time MOVES as the Newton iterates.

        ⚠ FRACTIONS, NOT TIMES, and that is load-bearing rather than a
        convenience.  An autonomous period is an unknown, so every step
        must scale with `T` or `dh/dT = h/T` -- the identity the period
        column rests on -- stops holding.  See `_period_grid`.

        `tstab` is how long to run before the window is taken; it defaults
        to 200 periods, which is a settling heuristic and not a
        convergence criterion.  ⚠ A RUN THAT HAS NOT SETTLED YIELDS A GRID
        FOR THE WRONG TRAJECTORY, silently -- the fractions will still sum
        to 1 and `solve` will still accept them.  Pass a longer `tstab`,
        or seed `x0` on the orbit, when the answer matters.
        """
        import warnings as _warnings
        from pycircuit.circuit.transient import Transient
        T = float(period)
        if not T > 0.0:
            raise ValueError('lte_grid: period must be positive, got %g' % T)
        tstab = 200.0 * T if tstab is None else float(tstab)
        if tstab < 0.0:
            raise ValueError('lte_grid: tstab must not be negative, got %g'
                             % tstab)
        rt = self.par.reltol if reltol is None else float(reltol)
        h0 = (T / 200.0) if timestep is None else float(timestep)

        tr = Transient(self.cir, toolkit=self.toolkit, reltol=rt)
        with _warnings.catch_warnings():
            _warnings.simplefilter('ignore')
            res = tr.solve(refnode=refnode, tend=tstab + T, timestep=h0,
                           x0=x0)
        t = np.asarray(res.sweep_values, dtype=float).ravel()
        xs = np.asarray(res.x, dtype=float)
        ## a settled window of exactly one period, taken from the END
        j0 = int(np.searchsorted(t, t[-1] - T))
        win_t, win_x = t[j0:], xs[:, j0:]
        if len(win_t) < 3:
            raise RuntimeError(
                'lte_grid: the adaptive run put only %d points in the last '
                'period, which is not a grid. Either the transient took '
                'steps larger than the period (raise tstab or lower '
                'timestep) or the period given is wrong.' % len(win_t))
        hs = np.diff(win_t)
        total = float(hs.sum())
        if total <= 0.0:
            raise RuntimeError('lte_grid: the derived window has zero span.')
        fr = hs / total
        ## the seed is the state at the window's start, with the reference
        ## row removed -- the shape `solve(x0=...)` takes
        iref = self.cir.get_node_index(refnode)
        seed = np.concatenate((win_x[:iref, 0], win_x[iref + 1:, 0]))
        return fr, seed

    def _period_grid(self, period, npts, grid):
        """`(times, hs)` for one period -- uniform, or a caller's own grid.

        RECORDED SCOPE ITEM 5.  A transient adapts because it cannot see
        the future; PSS re-solves the SAME interval over and over, so it can
        be handed a grid that was chosen well ONCE and then frozen.  The
        grid still never moves inside a solve, so the shooting Newton stays
        exact -- freezing is what makes it a Newton, and this changes only
        WHICH frozen grid.

        ⚠ THE RECORDED BLOCKER WAS STALE.  Item 5 said this was "blocked on
        `Transient` accepting a non-uniform grid; `fixed_timestep` is
        uniform-only".  `Transient.solve`'s loop is uniform-only and always
        was -- but PSS never uses that loop.  It drives `solve_timestep`
        directly, one step at a time, and non-uniform steps worked through
        that path unchanged.  Verified before any of this was written.

        `grid` is a sequence of step FRACTIONS of the period, summing to 1.
        Fractions rather than times because an autonomous period is an
        unknown: every step must scale with `T`, or `dh/dT = h/T` -- the
        identity the period column rests on -- stops holding.

        Measured on van der Pol at mu=100 (`benchmarks/pss_lte_grid.py`):
        1105 steps taken from an adaptive transient converge where 1105
        UNIFORM steps do not, and beat a 20000-point uniform grid on
        accuracy -- 18x fewer points and -47.3 ppm against -60.6.
        """
        if grid is None:
            times, dt = self.toolkit.linspace(0.0, period, num=npts,
                                              endpoint=True, retstep=True)
            return times, np.full(len(times), float(dt))
        fr = np.asarray(grid, dtype=float)
        if fr.ndim != 1 or len(fr) < 2:
            raise ValueError('grid must be a 1-D sequence of at least two '
                             'step fractions, got shape %r' % (fr.shape,))
        if not np.all(fr > 0.0):
            raise ValueError('every grid step fraction must be positive; '
                             'the smallest given is %g' % float(fr.min()))
        total = float(fr.sum())
        if abs(total - 1.0) > 1e-9:
            raise ValueError(
                'grid step fractions must sum to 1 (they are fractions of '
                'the period, so that every step scales with T when the '
                'period is an unknown); they sum to %.12g' % total)
        ## ⚠ THE OPENING STEP IS MANUFACTURED, SO IT MUST NOT BE THE GRID'S
        ## COARSE END.  `_traverse` builds `x(0)` from the unknown with ONE
        ## order-dropped step of `hs[0]`, and a grid taken from an adaptive
        ## transient opens wherever that transient's window happened to
        ## start -- on van der Pol at mu=100, `h[0] = 1.4845` against a
        ## MEDIAN of 4.62e-04, 3200x coarser.  That single Euler step moves
        ## the state 7.4%, and the shooting Newton then has to invert a map
        ## whose first act is that step.  Opening at the grid's own finest
        ## step instead costs ONE extra step in 1105 and is the difference
        ## between not converging and converging.
        ##
        ## ⚠ THE 8x IS MEASURED, NOT DERIVED, and it is a guard rather than
        ## a threshold: it exists so grids that already work are left
        ## exactly as the caller wrote them ('2:1' opens at 2x its finest,
        ## 'smooth' at 5x, and neither needs this).  The falsifier is in
        ## the record: on van der Pol's grid, an opening step of 1e-1 still
        ## fails for gear and 1e-2 converges, against a ratio here of
        ## 13939.  Anything between those bounds separates the two cases.
        ## ⚠ AND IT IS ONLY NEEDED WHEN THERE IS A MANUFACTURED STEP TO
        ## PROTECT.  The subdivision exists because `_traverse` builds `x(0)`
        ## with one order-dropped Euler step of `hs[0]` FROM `x_in` -- an
        ## iterate that may be far from the orbit -- and a coarse `hs[0]`
        ## there defeats the inner Newton.  With `x0_unknown` the first step
        ## starts ON the orbit and the same coarse step is solvable:
        ## measured on van der Pol's own LTE grid, the raw 1105-step grid
        ## converges and reaches -47.3 ppm where the subdivided 1106-step
        ## one reaches -73.8.  So the subdivision COSTS accuracy, and it is
        ## skipped where it buys nothing.
        if fr[0] > 8.0 * fr.min() and not getattr(self, '_open_at_x0', False):
            d = float(fr.min())
            fr = np.concatenate(([d, fr[0] - d], fr[1:]))

        ## ⚠ A CALLER'S GRID CAN SILENTLY DEMOTE GEAR-2 TO FIRST ORDER.
        ## `_period_grid` validated positivity and sum-to-1 and nothing
        ## about the INTERIOR ratios.  A two-step method is zero-stable only
        ## up to `h_n / h_{n-1} = 1 + sqrt(2)`, and past it the integrator's
        ## own guard drops the step to Euler -- correct, and invisible.
        ## Measured on a Q=20 resonator driven at resonance with an
        ## alternating 3:1 grid, where half the steps are demoted:
        ##
        ##       npts   uniform    3:1 grid
        ##        100   19.91489    7.99821
        ##        200   20.00960   11.42923
        ##        400   20.02218   14.54985
        ##        800   20.02443   16.85280
        ##
        ## against an analytic peak of 20 V -- 60% low at 100 points,
        ## crawling up at FIRST order, and `converged = True` every time.
        ## Refining does not fix it, because refining a 3:1 grid keeps it
        ## 3:1.  So the warning names the ratio rather than suggesting a
        ## smaller step, which is the advice that does not work here.
        ##
        ## ⚠ This is what item 5 removed the premise for: the literature
        ## note in the class docstring argues Wambacq's objections to
        ## non-uniform BDF "do not bite inside a run" because the grid is
        ## UNIFORM and frozen.  A caller's grid is frozen but not uniform.
        if len(fr) > 1 and self._companion_reach() >= 2:
            from pycircuit.circuit.integrator import ZERO_STABILITY_RATIO
            ratios = fr[1:] / fr[:-1]
            worst = float(np.max(ratios))
            if worst > ZERO_STABILITY_RATIO:
                n_bad = int(np.sum(ratios > ZERO_STABILITY_RATIO))
                warnings.warn(
                    'PSS: this grid steps up by %.3fx where a two-step '
                    'method is zero-stable only to %.3fx, at %d of %d '
                    'interior ratios. Those steps are dropped to Euler by '
                    'the integrator, so the run is first-order there and '
                    'the answer can be far low while reporting converged -- '
                    'measured 60%% low on a Q=20 resonator with an '
                    'alternating 3:1 grid. Refining will NOT fix it: a '
                    'refined 3:1 grid is still 3:1. Smooth the grid so '
                    'adjacent steps stay within %.3fx, or use a one-step '
                    "method (method='trap')."
                    % (worst, ZERO_STABILITY_RATIO, n_bad, len(ratios),
                       ZERO_STABILITY_RATIO), RuntimeWarning, stacklevel=3)

        hs = fr * period
        times = np.concatenate(([0.0], np.cumsum(hs)))
        return times, hs

    def _companion_reach(self):
        """How many charges back the chosen integrator's companion reads.

        The mechanistic property that decides whether this analysis needs
        the entering history as an unknown: a method reaching one charge
        back can be started from a single unknown, one reaching two
        cannot.  Asked of the
        integrator rather than inferred from `method`, so a fourth method
        arrives with the right answer instead of the default one.
        """
        ## Polymorphic: the method answers.  A stage method returns 1 (self
        ## starting, reads only x_n); an LMM computes it from its companion
        ## coefficients.  No isinstance/name branch to extend per method.
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        return integ.companion_reach()

    def _solves_history(self):
        """Whether the period map needs the entering history as an unknown.

        MEASURED, NOT ASSUMED (`benchmarks/pss_seam_cost.py`).  A method
        whose companion reads only `q_{n-1}` cannot see the fabricated
        opening history at all -- euler's seam costs 5.1e-12 V and
        trapezoidal's 1.3e-11 V, both zero -- so enlarging their system
        would quadruple the shooting solve to fix nothing.  Gear-2 reads
        `q_{n-2}`, which in the plain formulation is the entering stand-in,
        and pays 1.266e-01 V at 100 points/period: 54% of its total error,
        rising to 73% at 400 as the seam falls one order slower than the
        interior.

        Autonomous runs take it too, through the composed system
        (`func_autonomous_solved_history`): a free period does not remove
        the need for a history the companion can read.  See 4c in the class
        docstring for what the seam does to an oscillator, which is NOT what
        it does to a driven circuit -- 0.75% of the error against 54%,
        landing in the orbit's shape rather than its frequency.
        """
        ## ⚠ TRAPEZOIDAL CANNOT JOIN THIS FORMULATION BY SOLVING FOR
        ## `x_{-1}`, and trying it is how that was learned.  `iq_{-1}` is
        ## exactly derivable from `x_{-1}` -- `-(i(x_{-1}) + u)`, item 4d --
        ## but the derivative that matters runs the other way: a one-step
        ## companion reads ONLY `iq_{-1}`, so the trajectory depends on
        ## `x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
        ## SINGULAR wherever a node carries no conductance -- every purely
        ## reactive node, which is most of a resonator.  Adding `x_{-1}` as
        ## m unknowns then makes the 2m x 2m system rank-deficient, and it
        ## fails exactly as it should: `LinAlgError: Singular matrix`, on 25
        ## tests at once.
        ##
        ## The right second unknown for a `b != 0` method is `iq_{-1}`
        ## ITSELF -- the `(x, iq)` state its monodromy already uses -- with
        ## the closure `iq_{-1} = iq_{N-1}`.  See item 5's note; not built.
        return self._companion_reach() >= 2

    def _step_sensitivity(self, Px, Cs, Pq, Jf, C_new, solve=None,
                          coeffs=None, source=None):
        """One step of the sensitivity recursion, for ANY seed width.

        ONE RECURSION FOR EVERY METHOD (and now for either formulation).
        Each writes its companion as `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`,
        so differentiating the step gives

            S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            P_n  = -Jf_n^-1 S
            Pq_n = a_0 C_n P_n + S

        `P` is `d x_j / d(unknowns)`: one block wide in the plain
        formulation, two when the history is solved for.  Nothing here
        depends on that
        width, which is why the two systems share this and not a copy.

        ⚠ A SOLVE, NOT AN INVERSE (stage 11).  `inv(Jf) @ ...` formed a dense
        inverse per timestep per iteration and squared the condition number
        it then multiplied through.

        `solve` overrides how that solve is taken, and exists so the
        MATRIX-FREE path (item 6) can hand in a PRE-FACTORED `Jf` without
        this recursion being copied.  It is the same recursion either way,
        which is the point: a second copy would be a second thing to get
        wrong, and this one is already shared by every method and both
        formulations.  Without it, matrix-free would refactor `Jf` once per
        step PER KRYLOV ITERATION -- `k` times the factorisations the dense
        path takes, which is worse than the problem it set out to fix.
        """
        ## `coeffs` overrides the LIVE `_coeffs` for the same reason `solve`
        ## overrides the solve: a matrix-free replay happens after the run
        ## that produced the steps, when `_coeffs` no longer describes the
        ## step being replayed.  See `_traverse_factored_plain`.
        ## ⚠ `source` ENTERS THE SOLVE AND NOT `Pq`, and the asymmetry is
        ## the physics rather than a convenience.  A small-signal source
        ## appears in the step's residual -- `Jf dx + S + du = 0` -- but NOT
        ## in the companion, because `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`
        ## is built from CHARGES, and an injected current is not one.
        ## Adding it to `Pq` as well would feed a fictitious charge forward
        ## into every later step, and the error would grow along the period
        ## rather than announce itself.
        ##
        ## This is what makes PAC share the recursion instead of copying it:
        ## the homogeneous propagation (`source=None`) is the monodromy and
        ## the driven one is the forced response, and they differ by this
        ## one term.
        alphas, b = self._coeffs if coeffs is None else coeffs
        S = b * Pq if b else np.zeros_like(Px[0])
        for k in range(1, len(alphas)):
            S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
        S_solve = S if source is None else S + source
        if solve is None:
            ## ⚠ THROUGH THE CALLER'S SOLVER, not `toolkit.linearsolver`.
            ## This is the DENSE propagation -- the thing matrix-free is
            ## measured against -- and it used to be hardcoded to the
            ## toolkit, so `linearsolver=SuperLUSolver()` reached the inner
            ## Newton (once forwarded) and never the propagation.  Comparing
            ## a sparse matrix-free path against a dense baseline would have
            ## flattered it; both sides go through the same strategy now.
            ## `DenseSolver` IS `toolkit.linearsolver`, so the default is
            ## unchanged.
            Px_new = -self._get_linearsolver().solve(Jf, S_solve, self.toolkit)
        else:
            Px_new = -solve(S_solve)
        Pq_new = alphas[0] * (C_new @ Px_new) + S
        return Px_new, Pq_new

    def _traverse_solved_history(self, x0_in, xm1_in, times, hs,
                                 T=None, want_dT=False):
        """One period from a SOLVED history, for a two-step companion.

        The plain `_traverse` solves for a single entering state `x_in` and
        manufactures `x(0)` from it with one order-dropped Euler step.  That
        is sound for a companion reaching one charge back and measurably
        wrong for one reaching two: the shooting condition constrains
        `x(0) = x(P)`, it does NOT constrain `x_in` to be the orbit's own
        `x(-dt)`, so `x_in` is an O(h^2) stand-in -- and Gear-2 reads it as
        a history point.

        Here BOTH `x(0)` and `x(-dt)` are unknowns and both are required to
        close:

            F = [ x_{N-1} - x_0 ,  x_{N-2} - x_{-1} ]

        which is periodicity of the whole state a two-step method needs,
        rather than of one slice of it.  The trajectory then opens at full
        order off a history the solve is responsible for, so there is no
        opening step to drop and no fabricated charge to read.

        Returns `(x_last, x_prev, P_last, P_prev)`, the `P` being
        `d x / d(x_0, x_{-1})` as one `n x 2n` block -- and with `want_dT`,
        two more entries: `d x_{N-1}/dT` and `d x_{N-2}/dT`, for the
        autonomous system, where the period is an unknown too.  BOTH rows
        need a period column, which is the one thing the composed system
        needs that neither enlargement carried alone.
        """
        toolkit = self.toolkit
        self._want_dfdh = want_dT
        toolkit = self.toolkit
        m = self.cir.n - 1

        ## THE HISTORY IS INSTALLED, NOT SEEDED.  `_begin_run(x_{-1})` opens
        ## the rings on the earlier point and the push puts `x_0` in front of
        ## it, so the first real step reads `q(x_0)` and `q(x_{-1})` -- two
        ## genuine solved points.  The flags then say what is true of them:
        ## a step of `dt` has been taken (`_dt_last`), the run is no longer
        ## opening (`_is_first_step`, `_no_history`), and `_dt_last2` stays
        ## None because the THIRD charge is still `q(x_{-1})` repeated -- the
        ## LTE estimator differences three, so its opening reading remains
        ## unsound and the report goes on discarding it.
        tr = self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])

        ## `P_0 = [I 0]`, `P_{-1} = [0 I]` -- the two unknowns, exactly.  The
        ## plain path seeds BOTH rings with `I`, which is the flat-history
        ## assumption written into the Jacobian; here there is nothing to
        ## assume.
        eye = np.asarray(toolkit.eye(m))
        zero = np.zeros((m, m))
        Px = [np.hstack((eye, zero)), np.hstack((zero, eye))]
        Cs = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        ## `Pq` is `d(iq_{-1})/d(x_0, x_{-1})`.  For `b = 0` the recursion
        ## never reads it and zero is right.  For `b != 0` it is NOT zero,
        ## and it is exactly differentiable because `_install_history` seeds
        ## `iq_{-1} = -(i(x_{-1}) + u)` from the DAE:
        ##
        ##     d(iq_{-1})/d x_{-1} = -G(x_{-1}),   d(iq_{-1})/d x_0 = 0
        ##
        ## Leaving it at zero would make the Jacobian wrong for trapezoidal
        ## in the same way the plain path's flat seed is wrong -- which is
        ## the whole defect this formulation removes.
        Pq = np.zeros((m, 2 * m))

        ## Kept for PAC as the plain path does -- but opening EMPTY, because
        ## on this path `x_0` is an unknown rather than the result of a step,
        ## so there is no solved `(C, Jf)` pair at it to record.  The lists
        ## therefore align with `times[1:]`, one shorter than `X`.  PAC is
        ## withdrawn (`test_PAC_is_withdrawn`) and nothing else reads them;
        ## stated here so its rewrite does not read a stale alignment out of
        ## the plain path's shape.
        self.Cvec = []
        self.Jtvec = []
        self.times = times

        ## The period column, propagated by the SAME recursion with one
        ## extra source term.  Zero at the start: neither unknown depends on
        ## T -- they are states, and the solve owns them.  (In the plain
        ## autonomous system the single entering unknown is likewise
        ## T-independent, so this opens the same way.)
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)

        x, x_prev = copy(x0_in), copy(xm1_in)
        P_prev = Px[1]
        for _j, t in enumerate(times[1:]):
            dt = hs[_j]
            x_prev = x
            x = copy(self.solve_timestep(x, t, dt))
            self.Cvec.append(copy(self._C))
            self.Jtvec.append(copy(self._Jf))

            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C)
            Px_new, Pq = self._step_sensitivity(Px, Cs, Pq, Jf, C_new)

            if want_dT:
                ## Every step scales together (`h = T/(N-1)`), so `dh/dT =
                ## h/T`, and `df/dh` at fixed solution is Fang's `p` --
                ## `residual_dh`, already shared.  For an AUTONOMOUS circuit
                ## its `du/dt` half vanishes, which is what makes solving
                ## for the period tractable at all.
                alphas, b = self._coeffs
                St = b * Pqt if b else np.zeros_like(Pt[0])
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -toolkit.linearsolver(Jf, St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]

            P_prev = Px[0]
            Px = [Px_new, Px[0]]
            Cs = [copy(C_new), Cs[0]]

        self._want_dfdh = False
        ## ⚠ THE FULL 2m x 2m MAP, not the `d x_{N-1}/d x_0` corner.  For a
        ## two-step method the one-period map acts on the PAIR, and its
        ## spectrum is the physical Floquet multipliers TOGETHER WITH the
        ## parasitic roots the k-step discretisation introduces (see the
        ## literature note in the class docstring -- controlling those roots
        ## is the whole subject).  The autonomous unit-circle eigenvalue is
        ## in there.  ⚠ `spectral_radius` NO LONGER TAKES A MAXIMUM OVER
        ## THIS MIXTURE (2026-09-02): `_spectral_report` separates the
        ## physical multipliers from the parasitic roots by the
        ## eigenvector's block structure and reports the maximum over the
        ## physical ones, so a method whose spurious root sits near the
        ## unit circle cannot have it read back as the orbit's stability.
        ## ⚠ ALWAYS THE FULL 2m x 2m MAP, NEVER THE `d x_{N-1}/d x_0`
        ## CORNER.  This used to hand the corner back on the driven path
        ## (`Px[0][:, :m]`), and a sub-block of a sensitivity is not a
        ## monodromy: its eigenvalues mean nothing.  Measured, it reported
        ## `spectral_radius` 1.279605 for the Q=20 resonator -- ABOVE ONE,
        ## reading as an unstable orbit -- where the analytic per-period
        ## decay is exp(-pi/Q) = 0.854636 and every other path reports
        ## 0.855.  For a two-step method the one-period map acts on the
        ## PAIR, so the monodromy is the pair's.
        ##
        ## Its spectrum carries the parasitic roots of the discretisation
        ## alongside the physical multipliers -- but they are not a problem
        ## here and the reason is quantitative: BDF-2's parasitic root is
        ## 1/3 per STEP (roots of `1.5z^2 - 2z + 0.5` are 1 and 1/3), so
        ## over a period it is (1/3)^N, which at 200 points is ~1e-95.
        ## Measured on the autonomous element the 16x16 spectrum is
        ## [1.000, 3.5e-06, 5.2e-16, 4.7e-17, 0, ...]: one physical unit
        ## eigenvalue and nothing else above rounding, so `max |eig|` picks
        ## the physical one.  A method whose parasitic root sat nearer the
        ## unit circle would need that separated; Gear-2's does not.
        self._monodromy = np.vstack((Px[0], Px[1]))
        if want_dT:
            return x, x_prev, Px[0], Px[1], Pt[0], Pt[1]
        return x, x_prev, Px[0], P_prev

    def _factorise(self, Jf):
        """One step's `Jf`, factored by the CALLER'S linear solver.

        ⚠ THIS USED TO REACH FOR `scipy.linalg.lu_factor` DIRECTLY, which
        made every matrix-free run dense-LAPACK whatever `linearsolver=`
        said -- and a circuit Jacobian at m=1002 is very sparse, so a dense
        LU is ~3e8 flops per step, `N` times over.  The recorded 2.13x at
        m=1002 and the m~250 gate were therefore both measured against a
        DENSE baseline; see `benchmarks/pss_matrix_free_ceiling.py` for what
        they become when both sides get a sparse solver.

        A solver whose `factor` returns `None` (the symbolic toolkits, and
        any solver that has not implemented one) falls back to `solve` per
        replayed step -- correct, and `k` times the factorisations, which is
        the cost matrix-free exists to avoid.  It is a fallback, not a mode
        to run in.
        """
        solver = self._get_linearsolver()
        fac = solver.factor(Jf, self.toolkit)
        if fac is not None:
            return fac
        toolkit, A = self.toolkit, Jf

        class _PerSolve(object):
            def solve(self, b):
                return solver.solve(A, b, toolkit)
        return _PerSolve()

    def _traverse_factored(self, x0_in, xm1_in, times, hs, T=None,
                           want_dT=False):
        """One period WITHOUT the sensitivities, keeping each step factored.

        RECORDED SCOPE ITEM 6, the trajectory half.  `_traverse_solved_history`
        does two things at once: it walks the period, and it propagates a
        `2m`-column sensitivity alongside.  Matrix-free needs the walk
        without the propagation, because the propagation is the thing it
        replaces -- so this returns what a MATVEC needs to replay the same
        steps: the two opening capacitances, and per step a FACTORED `Jf`
        with its `C`.

        ⚠ THE FACTORISATION IS THE POINT, not an optimisation on top.  The
        matvec runs one solve per step per Krylov iteration; refactoring
        `Jf` each time would cost `k` times the factorisations the dense
        path takes and lose before it started.  Factoring once here makes
        every later solve a back-substitution.

        ⚠ AND THE COST OF THAT IS MEMORY: `N` factorisations and `N`
        capacitances, `2 N m^2` doubles, where the dense path holds `O(m^2)`
        at a time.  At `m = 1002` and 50 points that is ~800 MB.  This is
        the trade matrix-free makes here and it is not free; a caller with a
        long period and a large circuit can run out of memory where the
        dense path merely ran slowly.
        """
        m = self.cir.n - 1
        self._want_dfdh = want_dT
        self._install_history(x0_in, xm1_in, hs[0], h_prev=hs[-1])
        C0 = [np.asarray(self._C_at(x0_in)), np.asarray(self._C_at(xm1_in))]
        Cs = list(C0)
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        steps = []
        x, x_prev = copy(x0_in), copy(xm1_in)
        for _j, t in enumerate(times[1:]):
            x_prev = x
            x = copy(self.solve_timestep(x, t, hs[_j]))
            ## ⚠ THE COEFFICIENTS ARE STORED PER STEP, like `Jf` and `C`.
            ## `_coeffs` is live state; on THIS path the history is
            ## installed so no step is order-dropped and the pair is
            ## constant, which is why a single post-run snapshot was right
            ## here.  It was right by luck rather than by construction --
            ## the plain path, whose opening IS dropped to Euler, showed
            ## what that luck is worth.  See `_traverse_factored_plain`.
            alphas, b = self._coeffs
            C_new = np.asarray(self._C).copy()
            lu = self._factorise(np.asarray(self._Jf))
            steps.append((lu, C_new, alphas, b))
            if want_dT:
                ## BOTH ROWS NEED A PERIOD COLUMN -- `dx_{N-1}/dT` and
                ## `dx_{N-2}/dT` -- and the ring carries both, so this is the
                ## same one-column recursion the plain path runs, kept to the
                ## end.  Neither depends on the Krylov direction, so both are
                ## computed once per Newton iteration rather than per GMRES
                ## iteration.
                St = b * Pqt if b else np.zeros(m)
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -lu.solve(St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]
            Cs = [C_new, Cs[0]]
        self._want_dfdh = False
        if want_dT:
            return C0, steps, x, x_prev, Pt[0], Pt[1]
        return C0, steps, x, x_prev

    def _monodromy_matvec(self, C0, steps, v):
        """`M v` for the solved-history map, replaying stored steps at width 1.

        THE SAME RECURSION THE DENSE PATH USES -- `_step_sensitivity`, with
        its solve pointed at the stored factorisation.  The dense path seeds
        it with `P_0 = [I 0]`, `P_{-1} = [0 I]` and carries `2m` columns;
        here it is seeded with the two halves of `v` and carries one, which
        is the whole difference between the two paths.

        `v` is `(v_0, v_{-1})`; the result is `(P_last v, P_prev v)`, so the
        shooting Jacobian's action is `v - M v` -- see `_matrix_free_solve`.
        """
        m = self.cir.n - 1
        ## ⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.
        ## `Jf` and `C` are real, so `M` is a REAL linear map and
        ## `M(a + ib) = Ma + i Mb` exactly.  PAC needs complex products
        ## (`I + alpha(f) H` with `alpha = exp(-2j pi f T)`), and the
        ## alternative -- factoring `Jf` in complex arithmetic -- would
        ## double the stored factors and roughly quadruple the solve cost
        ## for a map that has no imaginary part to begin with.  Splitting
        ## costs two back-substitutions against the SAME factors.
        ##
        ## ⚠ The float cast below is therefore a GUARD, not a convenience:
        ## it used to swallow a complex `v` silently by discarding the
        ## imaginary part, which is a wrong answer rather than an error.
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec(C0, steps, v.real)
                    + 1j * self._monodromy_matvec(C0, steps, v.imag))
        v = v.astype(float)
        Px = [v[:m].copy(), v[m:].copy()]
        Cs = list(C0)
        Pq = np.zeros(m)
        for lu, C_new, alphas, b in steps:
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _l.solve(S), coeffs=(alphas, b))
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        return np.concatenate((Px[0], Px[1]))

    def _traverse_factored_plain(self, x_in, T, times, hs, want_dT=False,
                                 open_at_x0=False):
        """The PLAIN path's trajectory pass, factored -- item 6 for one-step
        methods and for the systems whose unknown is a single entering state.

        Mirrors `_traverse`'s opening exactly, which is the whole
        requirement: the manufacturing step first, then BOTH rings seeded
        with the same `C` and `Pq = a_0 C` for a `b != 0` method.  Getting
        that opening wrong would give a matvec for a DIFFERENT map than the
        dense path's, and the two would disagree only in the third figure --
        the kind of difference a converged answer absorbs.

        With `want_dT` it also propagates the PERIOD COLUMN.  That column
        does not depend on the Newton direction, so it is computed once here
        rather than once per Krylov iteration, and returned as a vector --
        which is exactly what the bordered autonomous matvec needs.
        """
        toolkit = self.toolkit
        m = self.cir.n - 1
        self._want_dfdh = want_dT
        self._begin_period(x_in)
        if open_at_x0:
            ## no manufacturing step -- see `_traverse`, which this mirrors
            x = copy(x_in)
            x0 = copy(x_in)
            C_open = np.asarray(self._C_at(x_in)).copy()
        else:
            x = self.solve_timestep(x_in, times[0], hs[0])
            x0 = copy(x)
            C_open = np.asarray(self._C).copy()

        ## ⚠ THE COEFFICIENTS BELONG TO THE STEP, NOT TO THE RUN.  `_coeffs`
        ## is live state and the MANUFACTURING step is order-dropped, so it
        ## reports Euler's `(alphas, b)` -- `b = 0` -- where the loop steps
        ## report the method's own: measured on this ladder, trapezoidal
        ## opens at `((49000, -49000), 0.0)` and then runs at
        ## `((98000, -98000), -1.0)`.
        ##
        ## Both halves of that were got wrong here first, in opposite
        ## directions, with the trajectory itself matching to ZERO both
        ## times: reading the coefficients ONCE BEFORE THE LOOP applied the
        ## OPENING's to every step and put the period column 40-50% out for
        ## `trap` and `gear`; reading them once in the matvec applied the
        ## LOOP's to the opening and made `Pq` non-zero where `_traverse`
        ## seeds it at zero, which is a 100% error for `trap`.  `euler` was
        ## exact under both and would have passed a one-method test.
        ##
        ## ⚠ WHAT IS LOAD-BEARING IS THE OPENING PAIR.  Inside the loop the
        ## coefficients are CONSTANT for every method in this tree, so
        ## storing them per step is currently belt-and-braces -- a mutation
        ## replacing them with a post-run snapshot does NOT fail the tests,
        ## and that is recorded rather than hidden.  They are stored anyway
        ## because a variable-order method would make the loop vary too, and
        ## that failure would be silent.
        ## ⚠ `b_open = 0` when opening AT `x_0`: no companion current has
        ## been formed, so the matvec must seed `Pq` at zero to match
        ## `_traverse`'s exact seed rather than the manufactured one.
        ## ⚠ `_coeffs` DOES NOT EXIST YET on the `open_at_x0` path -- no step
        ## has run to set it -- so it must not be read at all, not even for
        ## the half that is then discarded.  `b_open = 0` makes the matvec
        ## seed `Pq` at zero, which is what `_traverse` does there, and
        ## `a_open` is unused in consequence.
        a_open, b_open = ((None, 0.0) if open_at_x0 else self._coeffs)
        ## ⚠ AND THE ONE THING `b_open = 0` DOES NOT COVER: a method that
        ## SEEDS a consistent `iq_{-1}` has formed a companion current at
        ## `x_0` after all, and it depends on `x_0`.  Carried in `opening`
        ## rather than read off `self` because a matrix-free replay happens
        ## after the run -- the same reason `steps` stores its coefficients.
        ## `None` for every other method, which keeps them bit-identical.
        pq_open = self._pq_seed_at_x0(x_in) if open_at_x0 else None
        Pt = [np.zeros(m), np.zeros(m)]
        Pqt = np.zeros(m)
        Cs = [C_open, C_open]
        steps = []
        for _j, t in enumerate(times[1:]):
            dt = hs[_j]
            x = copy(self.solve_timestep(x, t, dt))
            alphas, b = self._coeffs
            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C).copy()
            lu = self._factorise(Jf)
            steps.append((lu, C_new, alphas, b))
            if want_dT:
                ## the same recursion with the step size's own source term,
                ## one column wide -- see `_traverse`
                St = b * Pqt if b else np.zeros(m)
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -lu.solve(St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]
            Cs = [C_new, Cs[0]]
        self._want_dfdh = False
        return ((C_open, a_open, b_open, pq_open), steps, x0, x,
                (Pt[0] if want_dT else None))

    def _monodromy_matvec_plain(self, opening, steps, v):
        """`M v` for the plain path, one column through the stored steps."""
        m = self.cir.n - 1
        C_open, a_open, b_open, pq_open = opening
        ## ⚠ COMPLEX `v` IS TWO REAL REPLAYS, NOT A COMPLEX FACTORISATION.
        ## `Jf` and `C` are real, so `M` is a REAL linear map and
        ## `M(a + ib) = Ma + i Mb` exactly.  PAC needs complex products
        ## (`I + alpha(f) H` with `alpha = exp(-2j pi f T)`), and the
        ## alternative -- factoring `Jf` in complex arithmetic -- would
        ## double the stored factors and roughly quadruple the solve cost
        ## for a map that has no imaginary part to begin with.  Splitting
        ## costs two back-substitutions against the SAME factors.
        ##
        ## ⚠ The float cast below is therefore a GUARD, not a convenience:
        ## it used to swallow a complex `v` silently by discarding the
        ## imaginary part, which is a wrong answer rather than an error.
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec_plain(opening, steps, v.real)
                    + 1j * self._monodromy_matvec_plain(opening, steps, v.imag))
        v = v.astype(float)
        Px = [v.copy(), v.copy()]
        Cs = [C_open, C_open]
        ## ⚠ THE OPENING PAIR, NOT THE LOOP'S.  `_traverse` opens `Pq` at
        ## `a_0 C` for a `b != 0` method reading `_coeffs` right after the
        ## MANUFACTURING step -- which is order-dropped to Euler, where
        ## `b = 0`, so in practice `Pq` opens at ZERO for every method here.
        ## Using the loop's pair instead makes it non-zero and the matvec
        ## 100% wrong for `trap`; this is the half that is load-bearing.
        Pq = a_open[0] * (C_open @ v) if b_open else np.zeros(m)
        ## the consistent-`iq_0` seed, applied to this column -- see
        ## `_pq_seed_at_x0`.  `None` leaves the zero above untouched.
        if pq_open is not None:
            Pq = Pq + pq_open @ v
        for lu, C_new, alphas, b in steps:
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _l.solve(S), coeffs=(alphas, b))
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
        return Px[0]

    def _monodromy_matvec_transposed_plain(self, opening, steps, v,
                                           collect=False, inject=None):
        """`M^T v` for the PLAIN path — the one-step companions (B8).

        The forward recursion `_step_sensitivity` implements is

            S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            P_n  = -K S,          K = Jf_n^-1
            Pq_n = a_0 C_n P_n + S

        and for a ONE-STEP companion the sum has a single term, so the two
        shipped one-step methods transpose differently and both are cheap:

        **Euler** (`b = 0`): `Pq` never re-enters, the map is
        `P_n = -a_1 K C_{n-1} P_{n-1}`, and the transpose is one term:

            w <- -a_1 C_{n-1}^T K^T w

        **Trapezoidal** (`b = -1`): `Pq` DOES re-enter, so the state is the
        pair `(P, Pq)` — but a DIFFERENT pair from Gear-2's `(P_n, P_{n-1})`,
        which is why the solved-history replay cannot be reused. Writing
        `S = a_1 C_{n-1} P_{n-1} + b Pq_{n-1}`,

            P_n  = -K S,        Pq_n = (I - a_0 C_n K) S

        so the transpose acting on `(w1, w2)` shares one bracket

            r = w2 - K^T (w1 + a_0 C_n^T w2)
            (w1, w2) <- (a_1 C_{n-1}^T r,  b r)

        ⚠ ONE TRANSPOSED SOLVE PER STEP, the same cost as Gear-2 — the naive
        arrangement takes two (`K^T w1` and `K^T C_n^T w2` separately) and
        the factorisation above avoids it.

        ⚠⚠ `Pq` OPENS AT ZERO, which the forward plain replay documents as
        the load-bearing half: `_traverse` opens `Pq` right after the
        MANUFACTURING step, which is order-dropped to Euler where `b = 0`.
        So the seed enters through `P` alone and the backward pass reads its
        answer out of `w1`.

        ⚠ DERIVED HERE AND GATED AGAINST A DENSE REFERENCE, because a
        from-scratch adjoint derivation in this file has come out
        sign-inverted before (roadmap §0h): the asymmetry between "the
        derivative acts on the product `C x`" and "on `y` alone" is easy to
        carry over wrongly. The test builds `M` column by column from the
        FORWARD replay and compares `M^T`.
        """
        m = self.cir.n - 1
        C_open, _a_open, _b_open, pq_open = opening
        v = np.asarray(v)
        inj = None if inject is None else [np.asarray(z) for z in inject]
        if np.iscomplexobj(v) or (
                inj is not None and any(np.iscomplexobj(z) for z in inj)):
            ii = (None, None) if inj is None else (
                [z.real for z in inj], [z.imag for z in inj])
            re = self._monodromy_matvec_transposed_plain(
                opening, steps, np.asarray(v).real, collect, ii[0])
            im = self._monodromy_matvec_transposed_plain(
                opening, steps, np.asarray(v).imag, collect, ii[1])
            if collect:
                ## ⚠ The collected lists may be NESTED (a DIRK's `ts` holds
                ## per-stage solves per step), so a flat `a + 1j*b`
                ## multiplied a LIST by `1j`: `floquet_modes` under trbdf2
                ## raised "can't multiply sequence by non-int of type
                ## 'complex'" for as long as the path existed.  Recurse.
                return (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                        _cx_collect(re[2], im[2]))
            return re + 1j * im
        v = v.astype(float)
        if not steps:
            return (v.copy(), [], []) if collect else v.copy()
        ## ⚠ `w2` IS THE ADJOINT OF `Pq_0`, AND IT IS ONLY DISCARDABLE WHEN
        ## `Pq_0` DOES NOT DEPEND ON THE SEED.  The forward map opens at
        ## `P_0 = v`, `Pq_0 = pq_open v`, so the transpose closes at
        ## `w1 + pq_open^T w2` -- see `_pq_seed_at_x0`.  With `pq_open` None
        ## (every method but `theta`) the second term is absent and this
        ## returns `w1` exactly as it always did.
        ## `C_{n-1}` for each step: the previous step's `C_new`, or the
        ## opening capacitance for the first
        prevC = [C_open] + [np.asarray(st[1]) for st in steps[:-1]]

        w1 = v.copy()
        w2 = np.zeros(m)
        ts = []
        states = []
        for j in range(len(steps) - 1, -1, -1):
            lu, C_new, alphas, b = steps[j]
            if len(alphas) != 2:
                raise NotImplementedError(
                    'PSS: the plain transposed replay is derived for a '
                    'ONE-STEP companion (two alpha coefficients) and this '
                    'step has %d. A multistep method on the plain path needs '
                    'its own reverse recursion.' % len(alphas))
            Cn = np.asarray(C_new)
            rhs = w1 + (alphas[0] * (Cn.T @ w2) if b else np.zeros(m))
            t = lu.solve_transposed(rhs)
            if t is None:
                raise NotImplementedError(
                    'PSS: this linear solver cannot solve transposed, so the '
                    'monodromy transpose cannot be replayed. Use DenseSolver '
                    'or SuperLUSolver.')
            if collect:
                ts.append(t)
            r = (w2 - t) if b else (-t)
            Cp = np.asarray(prevC[j])
            w1 = alphas[1] * (Cp.T @ r)
            w2 = (b * r) if b else np.zeros(m)
            if inj is not None:
                w1 = w1 + inj[j]
            if collect:
                ## ⚠ THE PLAIN ADJOINT STATE IS WIDTH `m`, NOT `2m`.  The
                ## solved-history replay collects `concat(w1, w2)` because
                ## its state IS the pair; here `w2` is the companion term
                ## `Pq`, not a second state block, and `v(s_j)` is `w1`
                ## alone.  Callers slice `st[:m]`, which is the whole
                ## vector here and the first block there -- the same
                ## quantity under both maps, which is what lets `ppv`
                ## consume either without knowing which it has.
                states.append(w1.copy())
        if pq_open is not None:
            w1 = w1 + pq_open.T @ w2
        if collect:
            ## reversed so `ts[j]`/`states[j]` line up with `steps[j]`,
            ## matching `_monodromy_matvec_transposed`'s contract
            ts.reverse()
            states.reverse()
            return w1, ts, states
        return w1

    def _traverse_factored_full(self, x0_in, times, hs):
        """One period under Radau IIA(3), kept factored -- the `m x m`
        monodromy of a SELF-STARTING fully-implicit collocation method.

        ⚠ COUPLED, NOT A COMPANION SUM AND NOT A DIRK COMPOSITION.  The three
        stages are solved together each step, so the per-step Jacobian is the
        `3m x 3m` coupled operator, not a product of per-stage solves.
        Linearising the stage residuals
        ``F_i = q(Y_i) - q(x_n) - h sum_j A_ij K_j`` w.r.t. the entering
        ``x_n`` (``K_j = -(i(Y_j) + u)``, so ``dK_j/dY_j = -G(Y_j)``):

            J_block (dY/dx_n) = [Cn; Cn; Cn],
            J_block[i][j] = delta_ij C(Y_i) + h A_ij G(Y_j)

        and the step map is ``dx_{n+1}/dx_n = (dY/dx_n)_3`` -- the third
        `m`-block -- by stiff accuracy (``x_{n+1} == Y_3``).  Each step stores
        ``(lu_block, Cn, m)``: the dense `3m x 3m` factor and the entering
        capacitance.  ``C(Y_i)`` and ``G(Y_j)`` are at the THREE distinct
        stage points, which is why `_C_at`/`_G_at` exist rather than a stored
        `Geq`.  Verified against the pencil ``exp(mu T)`` before shipping.

        ⚠ NO OPENER, NO PAIR.  Like TR-BDF2, Radau is self-starting -- every
        step reads only ``x_n`` -- so the map is `m x m` and high-order all
        the way round, with no order-dropped opening seam inside the period.
        The dense `3m x 3m` factor is the coupled real solve; the
        ``A^{-1}``-eigenbasis transform (1 real + 1 complex LU) is the
        efficiency follow-up documented on `RadauIIA3Integrator`.
        """
        import scipy.linalg as sla
        from pycircuit.circuit.integrator import RungeKuttaIntegrator
        tr = self._transient()
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x0_in)
        integ = tr.base_integrator
        if not isinstance(integ, RungeKuttaIntegrator):
            raise ValueError('_traverse_factored_full needs a Runge-Kutta stage '
                             'inner integrator, got %r' % (integ,))
        Amat = np.array(integ.A, dtype=float)
        s = Amat.shape[0]
        iref = self.irefnode
        x = copy(x0_in)
        x_prev = copy(x0_in)
        steps = []
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            x_prev = xn
            Yf = tr._rk_Y
            Ys = [self.toolkit.concatenate((yf[:iref], yf[iref + 1:]))
                  for yf in Yf]
            Cn = np.asarray(self._C_at(xn))
            m = Cn.shape[0]
            Cs = [np.asarray(self._C_at(y)) for y in Ys]
            Gs = [np.asarray(self._G_at(y)) for y in Ys]
            Jb = np.zeros((s * m, s * m))
            for i in range(s):
                for jj in range(s):
                    blk = h * Amat[i, jj] * Gs[jj]
                    if i == jj:
                        blk = Cs[i] + blk
                    Jb[i * m:(i + 1) * m, jj * m:(jj + 1) * m] = blk
            lu = sla.lu_factor(Jb)
            steps.append((lu, Cn, m))
        self._want_dfdh = False
        return steps, x, x_prev

    def _monodromy_matvec_full(self, steps, v):
        """`M v` for the Radau IIA(3) map, replaying the stored coupled
        factors.  Each step: stack the entering direction into the coupled
        RHS ``[Cn v; Cn v; Cn v]``, solve the `3m x 3m` system, and carry the
        THIRD `m`-block forward (``x_{n+1} == Y_3``).  Real map, so a complex
        `v` splits into two real replays exactly (the same guard the LMM and
        TR-BDF2 matvecs use -- a float cast would silently drop the imaginary
        part)."""
        import scipy.linalg as sla
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec_full(steps, v.real)
                    + 1j * self._monodromy_matvec_full(steps, v.imag))
        w = v.astype(float)
        for lu, Cn, m in steps:
            s = lu[0].shape[0] // m
            cw = Cn @ w
            Z = sla.lu_solve(lu, np.concatenate([cw] * s))
            w = Z[(s - 1) * m:s * m]
        return w

    def _monodromy_matvec_transposed_full(self, steps, v, collect=False,
                                           inject=None):
        """`M^T v` for the Radau IIA(3) map -- the adjoint of the coupled
        step, replayed in reverse step order.

        The forward step is ``M_j = E3 J^{-1} (1_3 (x) Cn)`` with ``E3`` the
        third-block extractor and ``1_3 (x) Cn`` the stacking of ``Cn``, so

            p = J^{-T} [0; 0; w]
            M_j^T w = Cn^T (p_1 + p_2 + p_3)

        -- one coupled transposed solve per step.  ``M = M_{N-1} ... M_0`` so
        ``M^T = M_0^T ... M_{N-1}^T`` and the loop runs LAST step to first.

        With `collect`, returns ``(w, ts, states)`` where ``states[j]`` is the
        adjoint state after step `j` (the PPV over the period that `ppv`
        reads) and ``ts[j]`` is the FULL `3m` coupled transposed solve ``p``
        at that step -- what the Radau forced/sideband fold needs, since a
        source injected at step `j` enters through all three stages
        (``A (x) B``), so there is no two-vector shortcut."""
        import scipy.linalg as sla
        v = np.asarray(v)
        if np.iscomplexobj(v):
            ii = (None, None) if inject is None else (
                np.real(inject), np.imag(inject))
            re = self._monodromy_matvec_transposed_full(
                steps, v.real, collect, ii[0])
            im = self._monodromy_matvec_transposed_full(
                steps, v.imag, collect, ii[1])
            if collect:
                ## ⚠ The collected lists may be NESTED (a DIRK's `ts` holds
                ## per-stage solves per step), so a flat `a + 1j*b`
                ## multiplied a LIST by `1j`: `floquet_modes` under trbdf2
                ## raised "can't multiply sequence by non-int of type
                ## 'complex'" for as long as the path existed.  Recurse.
                return (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                        _cx_collect(re[2], im[2]))
            return re + 1j * im
        v = v.astype(float)
        if not steps:
            return (v.copy(), [], []) if collect else v.copy()
        w = v.copy()
        ts = []
        states = []
        for j in range(len(steps) - 1, -1, -1):
            lu, Cn, m = steps[j]
            s = lu[0].shape[0] // m
            b3 = np.concatenate([np.zeros(m)] * (s - 1) + [w])
            p = sla.lu_solve(lu, b3, trans=1)
            w = Cn.T @ sum(p[k * m:(k + 1) * m] for k in range(s))
            if inject is not None:
                w = w + inject[j]
            if collect:
                ts.append(p)
                states.append(w.copy())
        if collect:
            ts.reverse()
            states.reverse()
            return w, ts, states
        return w

    ## ------------------------------------------------------------------
    ## The DIRK / ESDIRK (lower-triangular) shooting family -- SEQUENTIAL,
    ## tableau-generic over any number of stages.  A fully-implicit method
    ## uses the coupled `_*_full` family instead; the DAE forces the split
    ## (an explicit first stage makes the coupled block singular), and the
    ## sequential recursion below never forms that block.  See
    ## doc/integrator_architecture_260906.md.
    ## ------------------------------------------------------------------

    def _traverse_factored_dirk(self, x0_in, times, hs):
        """One period under a lower-triangular (DIRK/ESDIRK) stage method, kept
        factored -- the `m x m` monodromy, solved stage by stage.

        Each step stores, for every stage, the factor
        ``K_i = LU(C(Y_i) + h A_ii G(Y_i))`` (``None`` for an explicit stage),
        the stage conductances ``G(Y_i)``, the entering ``C(x_n)``, and ``h``.
        The monodromy recursion (differentiating the stage residuals w.r.t. the
        entering ``x_n``) is

            D_0 = I                                    (explicit first stage)
            D_i = K_i^{-1} (C_n - h sum_{j<i} A_ij G_j D_j)   (implicit)

        with ``x_{n+1} = Y_s`` (stiff accuracy).  No coupled ``sm`` block is
        ever formed, so an explicit first stage costs nothing and never makes a
        singular block -- the whole reason a DIRK is not routed through the
        `_*_full` coupled family on a DAE.
        """
        from pycircuit.circuit.integrator import RungeKuttaIntegrator
        tr = self._transient()
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x0_in)
        integ = tr.base_integrator
        if not isinstance(integ, RungeKuttaIntegrator):
            raise ValueError('_traverse_factored_dirk needs a Runge-Kutta '
                             'stage inner integrator, got %r' % (integ,))
        Amat = np.array(integ.A, dtype=float)
        s = Amat.shape[0]
        iref = self.irefnode
        x = copy(x0_in)
        x_prev = copy(x0_in)
        steps = []
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            x_prev = xn
            Yf = tr._rk_Y
            Ys = [self.toolkit.concatenate((yf[:iref], yf[iref + 1:]))
                  for yf in Yf]
            Cn = np.asarray(self._C_at(xn))
            Gs = [np.asarray(self._G_at(Ys[i])) for i in range(s)]
            Kfacs = []
            for i in range(s):
                if abs(Amat[i, i]) < 1e-14:
                    Kfacs.append(None)  # explicit stage
                else:
                    Ci = np.asarray(self._C_at(Ys[i]))
                    Kfacs.append(self._factorise(Ci + h * Amat[i, i] * Gs[i]))
            ## carry the tableau IN the step so the replay uses the FP's own
            ## method (a factored_period_dirk built with method= may differ
            ## from self.par.method -- the FULL matvec dodges this by replaying
            ## a stored factor, but the DIRK replay needs A and c).
            steps.append((Kfacs, Gs, Cn, float(h), Amat,
                          np.array(integ.C, dtype=float)))
        self._want_dfdh = False
        return steps, x, x_prev

    def _monodromy_matvec_dirk(self, steps, v):
        """`M v` for a lower-triangular stage map -- the sequential replay of
        the stored per-stage factors.  Real map, so a complex `v` splits into
        two real replays (the guard every stage/LMM matvec uses)."""
        v = np.asarray(v)
        if np.iscomplexobj(v):
            return (self._monodromy_matvec_dirk(steps, v.real)
                    + 1j * self._monodromy_matvec_dirk(steps, v.imag))
        w = v.astype(float)
        for Kfacs, Gs, Cn, h, Amat, _cvec in steps:
            s = Amat.shape[0]
            d = [None] * s
            for i in range(s):
                if Kfacs[i] is None:
                    if i != 0:
                        raise NotImplementedError(
                            'DIRK monodromy: an explicit stage other than the '
                            'first is not supported (its solve would need C^-1, '
                            'singular on a DAE).')
                    d[0] = w  # explicit first stage: Y_0 = x_n, so D_0 = I
                else:
                    rhs = Cn @ w - h * sum(Amat[i, j] * (Gs[j] @ d[j])
                                           for j in range(i))
                    d[i] = Kfacs[i].solve(rhs)
            w = d[s - 1]
        return w

    def _monodromy_matvec_transposed_dirk(self, steps, v, collect=False,
                                          inject=None):
        """`M^T v` for a lower-triangular stage map -- the reverse-mode adjoint
        of the sequential recursion, replayed last step to first.

        Within a step, seed the adjoint on the last stage and process stages in
        REVERSE: for an implicit stage ``i``,
        ``rbar = K_i^{-T} dbar_i``; ``wbar += C_n^T rbar``;
        ``dbar_j += -h A_ij G_j^T rbar`` for ``j<i``.  The explicit first stage
        contributes ``wbar += dbar_0``.  With `collect`, ``states[j]`` is the
        adjoint state after step `j` (the PPV) and ``ts[j]`` is the per-stage
        reverse solves the DIRK forced fold reads."""
        v = np.asarray(v)
        if np.iscomplexobj(v):
            ii = (None, None) if inject is None else (
                np.real(inject), np.imag(inject))
            re = self._monodromy_matvec_transposed_dirk(
                steps, v.real, collect, ii[0])
            im = self._monodromy_matvec_transposed_dirk(
                steps, v.imag, collect, ii[1])
            if collect:
                ## ⚠ The collected lists may be NESTED (a DIRK's `ts` holds
                ## per-stage solves per step), so a flat `a + 1j*b`
                ## multiplied a LIST by `1j`: `floquet_modes` under trbdf2
                ## raised "can't multiply sequence by non-int of type
                ## 'complex'" for as long as the path existed.  Recurse.
                return (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                        _cx_collect(re[2], im[2]))
            return re + 1j * im
        v = v.astype(float)
        if not steps:
            return (v.copy(), [], []) if collect else v.copy()
        w = v.copy()
        ts = []
        states = []
        for j in range(len(steps) - 1, -1, -1):
            Kfacs, Gs, Cn, h, Amat, _cvec = steps[j]
            s = Amat.shape[0]
            dbar = [np.zeros_like(w) for _ in range(s)]
            dbar[s - 1] = w
            wbar = np.zeros_like(w)
            rbars = [None] * s
            for i in range(s - 1, -1, -1):
                if Kfacs[i] is None:
                    if i == 0:
                        wbar = wbar + dbar[0]
                else:
                    rb = Kfacs[i].solve_transposed(dbar[i])
                    if rb is None:
                        raise NotImplementedError(
                            'PSS: this linear solver cannot solve transposed, '
                            'so the DIRK monodromy transpose cannot be '
                            'replayed. Use DenseSolver or SuperLUSolver.')
                    rbars[i] = rb
                    wbar = wbar + Cn.T @ rb
                    for jj in range(i):
                        dbar[jj] = dbar[jj] - h * Amat[i, jj] * (Gs[jj].T @ rb)
            w = wbar
            if inject is not None:
                w = w + inject[j]
            if collect:
                ts.append(rbars)
                states.append(w.copy())
        if collect:
            ts.reverse()
            states.reverse()
            return w, ts, states
        return w

    def _traverse_dirk(self, x_in, T, times, hs, want_dT=False):
        """One period under a lower-triangular (DIRK/ESDIRK) stage method with
        the DENSE sensitivities -- the shooting Newton's monodromy `P` and,
        with `want_dT`, the period column `Pt`, solved stage by stage.

        Same sequential recursion as :meth:`_monodromy_matvec_dirk`, carried on
        the full `m x m` matrix `P`:

            D_0 = I;  D_i = K_i^{-1}(C_n P - h sum_{j<i} A_ij G_j D_j);  P = D_s

        and the period column (autonomous, ``h_j = frac_j T``, ``dh/dT=h/T``,
        ``K_j = -i(Y_j)``):

            Dt_0 = Pt;  Dt_i = K_i^{-1}(C_n Pt + (h/T) S_i - h sum_{j<i} A_ij G_j Dt_j)

        with ``S_i = sum_{j<=i} A_ij K_j``.  FD-checked before use (the dT
        column has been got wrong in this file twice -- roadmap 0j)."""
        toolkit = self.toolkit
        m = self.cir.n - 1
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x_in)
        integ = self._transient().base_integrator
        Amat = np.array(integ.A, dtype=float)
        s = Amat.shape[0]
        iref = self.irefnode
        Tf = float(T)
        x = copy(x_in)
        x0 = copy(x_in)
        P = np.asarray(toolkit.eye(m), dtype=float)
        Pt = np.zeros(m)
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            Yf = self._transient()._rk_Y
            Ys = [toolkit.concatenate((yf[:iref], yf[iref + 1:])) for yf in Yf]
            Cn = np.asarray(self._C_at(xn))
            Gs = [np.asarray(self._G_at(Ys[i])) for i in range(s)]
            Kf = []
            for i in range(s):
                if abs(Amat[i, i]) < 1e-14:
                    Kf.append(None)
                else:
                    Ci = np.asarray(self._C_at(Ys[i]))
                    Kf.append(self._factorise(Ci + h * Amat[i, i] * Gs[i]))
            D = [None] * s
            CnP = Cn @ P
            for i in range(s):
                if Kf[i] is None:
                    D[0] = P  # explicit first: D_0 = I, so D_0 P_in = P_in
                else:
                    rhs = CnP - h * sum(Amat[i, j] * (Gs[j] @ D[j])
                                        for j in range(i))
                    ## solve the whole m x m RHS at once (the factor's solve
                    ## takes a 2D b), not column by column
                    D[i] = np.asarray(Kf[i].solve(rhs))
            P = D[s - 1]
            if want_dT:
                Ks = [-np.asarray(self._i_at(Ys[i])) for i in range(s)]
                Dt = [None] * s
                CnPt = Cn @ Pt
                for i in range(s):
                    if Kf[i] is None:
                        Dt[0] = Pt
                    else:
                        Si = sum(Amat[i, j] * Ks[j] for j in range(i + 1))
                        rhs = CnPt + (h / Tf) * Si \
                            - h * sum(Amat[i, j] * (Gs[j] @ Dt[j])
                                      for j in range(i))
                        Dt[i] = Kf[i].solve(rhs)
                Pt = Dt[s - 1]
        self._want_dfdh = False
        self._monodromy = P
        if want_dT:
            return x0, x, P, Pt
        return x0, x, P, None

    def _forced_replay_dirk(self, fp, freq, u_ac, y0=None, collect=False):
        """One driven period under a lower-triangular stage method -- the
        FORWARD sequential replay.  Source at `freq` enters stage `i`'s residual
        at every abscissa ``k <= i``:

            d_0 = y                                       (explicit first)
            d_i = K_i^{-1}(C_n y - h sum_{j<i} A_ij G_j d_j
                          - h sum_{k<=i} A_ik u e^{jw t_{n,k}})

        so ``y_end = M y0 + w(freq)`` by linearity.  The exact transpose of
        `_forced_replay_transposed_dirk`."""
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        tms = np.asarray(fp.times, dtype=float)
        y = (np.zeros(m, dtype=complex) if y0 is None
             else np.asarray(y0, dtype=complex).ravel().copy())
        ys = []
        for jstep, (Kfacs, Gs, Cn, h, Amat, cvec) in enumerate(fp.steps):
            ts = tms[jstep]
            s = Amat.shape[0]
            d = [None] * s
            Cny = Cn @ y
            for i in range(s):
                if Kfacs[i] is None:
                    d[0] = y
                else:
                    src = -h * sum(Amat[i, k] * u_ac * np.exp(jw * (ts + cvec[k] * h))
                                   for k in range(i + 1))
                    rhs = Cny - h * sum(Amat[i, j] * (Gs[j] @ d[j])
                                        for j in range(i)) + src
                    d[i] = self._csolve_fac(Kfacs[i], rhs)
            y = d[s - 1]
            if collect:
                ys.append(y.copy())
        return y, ys

    def _csolve_fac(self, fac, b):
        """Complex solve against a real factor object (`.solve` real-only):
        two back-substitutions."""
        b = np.asarray(b, dtype=complex)
        return fac.solve(b.real) + 1j * fac.solve(b.imag)

    def _forced_replay_transposed_dirk(self, fp, freq, xa):
        """`W^T xa` for a lower-triangular stage map -- the reverse-mode adjoint
        of `_forced_replay_dirk` (no output injection).

        Per step, the monodromy reverse pass yields the per-stage reverse solves
        ``rbar_i = K_i^{-T} dbar_i``; the source at abscissa `k` couples to every
        stage ``i >= k``, so

            acc += -h sum_k e^{jw t_{n,k}} sum_{i>=k} A_ik rbar_i

        and the costate propagates by ``w <- wbar``.  Exact transpose of the
        forward replay (dual-consistent to machine precision)."""
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        w = np.asarray(xa, dtype=complex).ravel().copy()
        acc = np.zeros(m, dtype=complex)
        tms = np.asarray(fp.times, dtype=float)
        for j in range(len(fp.steps) - 1, -1, -1):
            Kfacs, Gs, Cn, h, Amat, cvec = fp.steps[j]
            s = Amat.shape[0]
            ts = tms[j]
            dbar = [np.zeros(m, dtype=complex) for _ in range(s)]
            dbar[s - 1] = w
            wbar = np.zeros(m, dtype=complex)
            rbars = [None] * s
            for i in range(s - 1, -1, -1):
                if Kfacs[i] is None:
                    if i == 0:
                        wbar = wbar + dbar[0]
                else:
                    rb = self._csolve_fac_T(Kfacs[i], dbar[i])
                    rbars[i] = rb
                    wbar = wbar + Cn.T @ rb
                    for jj in range(i):
                        dbar[jj] = dbar[jj] - h * Amat[i, jj] * (Gs[jj].T @ rb)
            for k in range(s):
                tk = ts + cvec[k] * h
                coup = sum(Amat[i, k] * rbars[i] for i in range(k, s)
                           if rbars[i] is not None)
                if not np.isscalar(coup):
                    acc = acc - h * np.exp(jw * tk) * coup
            w = wbar
        return acc

    def _csolve_fac_T(self, fac, b):
        """Complex transposed solve against a real factor object."""
        b = np.asarray(b, dtype=complex)
        rr = fac.solve_transposed(b.real)
        if rr is None:
            raise NotImplementedError(
                'PSS: this linear solver cannot solve transposed, so the DIRK '
                'forced adjoint cannot be replayed. Use DenseSolver or '
                'SuperLUSolver.')
        return rr + 1j * fac.solve_transposed(b.imag)

    def _sideband_forced_dirk(self, fp, freq, l, d):
        """The forced (source-injected) part of a lower-triangular stage
        method's sideband row `l`, and the final costate `g` -- the injected
        sibling of `_forced_replay_transposed_dirk` (output functional `d`
        added to the costate AFTER each step's update, causality)."""
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        T = float(fp.T)
        w0 = 2.0 * np.pi / T
        N = len(fp.steps)
        tms = np.asarray(fp.times, dtype=float)
        d = np.asarray(d, dtype=complex).ravel()
        lam = np.zeros(m, dtype=complex)
        forced = np.zeros(m, dtype=complex)
        for j in range(N - 1, -1, -1):
            Kfacs, Gs, Cn, h, Amat, cvec = fp.steps[j]
            s = Amat.shape[0]
            ts = tms[j]
            dbar = [np.zeros(m, dtype=complex) for _ in range(s)]
            dbar[s - 1] = lam
            wbar = np.zeros(m, dtype=complex)
            rbars = [None] * s
            for i in range(s - 1, -1, -1):
                if Kfacs[i] is None:
                    if i == 0:
                        wbar = wbar + dbar[0]
                else:
                    rb = self._csolve_fac_T(Kfacs[i], dbar[i])
                    rbars[i] = rb
                    wbar = wbar + Cn.T @ rb
                    for jj in range(i):
                        dbar[jj] = dbar[jj] - h * Amat[i, jj] * (Gs[jj].T @ rb)
            for k in range(s):
                tk = ts + cvec[k] * h
                coup = sum(Amat[i, k] * rbars[i] for i in range(k, s)
                           if rbars[i] is not None)
                if not np.isscalar(coup):
                    forced = forced - h * np.exp(jw * tk) * coup
            lam = wbar + np.exp(-1j * (float(l) * w0
                                       + 2.0 * np.pi * float(freq)) * ts) / N * d
        return forced, lam

    def _i_at(self, x_reduced):
        """The reduced resistive current `i(x)` at a point.  For an
        AUTONOMOUS circuit `dq/dt = -i(x)` (no source term), which is the
        stage derivative the TR-BDF2 period column needs."""
        tr = self._transient()
        i = tr.cir.i(self._insert_refnode(x_reduced), tr.epar)
        iref = self.irefnode
        return self.toolkit.concatenate((i[:iref], i[iref + 1:]))

    def _traverse_full(self, x_in, T, times, hs, want_dT=False):
        """One period under Radau IIA(3) with the DENSE sensitivities -- the
        shooting Newton's monodromy for the fully-implicit collocation method.

        The coupled analogue of `_traverse_trbdf2`.  Each step propagates the
        full `m x m` monodromy `P = dx/dx0` (and, with `want_dT`, the period
        column `Pt = dx/dT`) through the `3m x 3m` coupled stage operator
        ``J_block[i][j] = delta_ij C(Y_i) + h A_ij G(Y_j)``:

            J_block (dY/dx0) = [Cn P; Cn P; Cn P],   P <- (dY/dx0)_3

        -- the monodromy is the third `m`-block by stiff accuracy
        (``x_{n+1} == Y_3``).  Self-starting: `x_in` IS `x_0`, so there is no
        opener seam and `M` is order-5 round the whole period.

        ⚠ THE PERIOD COLUMN IS TRACTABLE ONLY BECAUSE THE CIRCUIT IS
        AUTONOMOUS.  With `T` unknown the grid rebuilds as ``h_j = frac_j T``,
        so ``dh_j/dT = h_j/T`` and the stage derivative ``K_j = -i(Y_j)`` has
        no explicit time dependence.  Differentiating the stage residuals
        ``F_i = q(Y_i) - q(x_n) - h sum_j A_ij K_j`` w.r.t. `T`:

            J_block (dY/dT) = [Cn Pt + (h/T) sum_j A_ij K_j]_i,   Pt <- (dY/dT)_3

        Finite-difference checked before use (the dT column has been got wrong
        in this file twice -- roadmap 0j).
        """
        import scipy.linalg as sla
        toolkit = self.toolkit
        m = self.cir.n - 1
        self._want_dfdh = False
        self._want_lte = False
        self._begin_period(x_in)
        integ = self._transient().base_integrator
        Amat = np.array(integ.A, dtype=float)
        s = Amat.shape[0]
        iref = self.irefnode
        Tf = float(T)
        x = copy(x_in)
        x0 = copy(x_in)
        P = np.asarray(toolkit.eye(m), dtype=float)
        Pt = np.zeros(m)
        for _j, t in enumerate(times[1:]):
            h = hs[min(_j, len(hs) - 1)]
            xn = x
            x = copy(self.solve_timestep(xn, t, h))
            Yf = self._transient()._rk_Y
            Ys = [toolkit.concatenate((yf[:iref], yf[iref + 1:])) for yf in Yf]
            Cn = np.asarray(self._C_at(xn))
            Cs = [np.asarray(self._C_at(y)) for y in Ys]
            Gs = [np.asarray(self._G_at(y)) for y in Ys]
            Jb = np.zeros((s * m, s * m))
            for i in range(s):
                for jj in range(s):
                    blk = h * Amat[i, jj] * Gs[jj]
                    if i == jj:
                        blk = Cs[i] + blk
                    Jb[i * m:(i + 1) * m, jj * m:(jj + 1) * m] = blk
            lu = sla.lu_factor(Jb)
            CnP = Cn @ P
            Z = sla.lu_solve(lu, np.vstack([CnP] * s))
            P = Z[(s - 1) * m:s * m, :]
            if want_dT:
                Ks = [-np.asarray(self._i_at(y)) for y in Ys]
                rhs = np.zeros(s * m)
                CnPt = Cn @ Pt
                for i in range(s):
                    Si = sum(Amat[i, jj] * Ks[jj] for jj in range(s))
                    rhs[i * m:(i + 1) * m] = CnPt + (h / Tf) * Si
                Zt = sla.lu_solve(lu, rhs)
                Pt = Zt[(s - 1) * m:s * m]
        self._want_dfdh = False
        self._monodromy = P
        if want_dT:
            return x0, x, P, Pt
        return x0, x, P, None

    def _monodromy_matvec_transposed(self, C0, steps, v, collect=False,
                                     inject=None):
        """`M^T v` -- the same stored steps, REPLAYED BACKWARDS.

        A shooting monodromy is a product of per-step solves, so its
        transpose is that product in reverse order with each solve
        transposed.  For Gear-2 (`b = 0`, three alphas) the per-step state
        is the PAIR `(Px_n, Px_{n-1})` and the step map is

            B_n = [ -a1 Jf^-1 C_{n-1}   -a2 Jf^-1 C_{n-2} ]
                  [        I                    0         ]

        so `B_n^T (v1; v2)` is
        `(-a1 C_{n-1}^T Jf^-T v1 + v2 ; -a2 C_{n-2}^T Jf^-T v1)`.

        ⚠ THIS IS WHY IT COSTS NOTHING TO HAVE.  Demir & Roychowdhury
        (TCAD 22(2) 188-196) call reverse integration "often unavailable
        even in existing time-domain simulators", requiring "significant
        changes to core simulation routines" -- true of a forward-only
        DENSE implementation.  `_traverse_factored` already stores every
        step's factorisation, and every factorisation here already knows how
        to solve transposed, so the reverse pass needs no new integrator, no
        refactorisation and no second traversal.

        MEASURED against the dense `M^T` built from the forward matvec:
        agreement 1.8e-15, and the reverse pass costs 0.75x the forward one
        -- CHEAPER, because it does two `C^T` products against one shared
        transposed solve where the forward does two `C` products and a
        solve.

        ⚠ NOT A CAPABILITY, A BUILDING BLOCK.  It is the shared dependency
        of a PPV (Demir & Roychowdhury's reverse Jacobian `J_r`) and of
        adjoint noise, neither of which is built.  It exists because the
        spike that established it is worth keeping, and it is pinned by a
        test rather than left in a scratch file.
        """
        m = self.cir.n - 1
        ## complex `v` is two real reverse replays -- see the note in
        ## `_monodromy_matvec`; `M^T` is real for exactly the same reason
        ## `M` is, and adjoint noise needs complex products.
        ## ⚠ `collect` HANDS BACK THE PER-STEP TRANSPOSED SOLVES, and it
        ## exists so ADJOINT NOISE does not get a second copy of this
        ## recursion.  The sensitivity of the final state to a source
        ## injected at step `j` is exactly `-Jf_j^-T` applied to the adjoint
        ## state there -- which is `t` below, already computed.  A driven
        ## reverse replay is therefore this pass plus a weighted sum, not a
        ## new traversal, and the two cannot drift apart because there is
        ## only one of them.  See `_forced_replay_transposed`.
        ## ⚠ `inject` MAKES THE OUTPUT A FUNCTIONAL OVER THE WHOLE PERIOD
        ## RATHER THAN A VALUE AT ITS END, and that is the difference
        ## between "the response at `t = 0`" and a SIDEBAND coefficient.
        ## `H_l` is `(1/N) sum_n exp(-j l w0 t_n) d^T y_n`: every state on
        ## the trajectory contributes, so its adjoint takes `c_n` INTO the
        ## reverse state at every step instead of seeding once at the end.
        ## Same recursion, one added term -- which is why this is a
        ## parameter and not a second function to keep in step.
        ##
        ## `inject[j]` lands on the state AFTER `steps[j]` has been applied
        ## backwards, i.e. on `P_j`, and `v` seeds `P_N`.
        v = np.asarray(v)
        inj = None if inject is None else np.asarray(inject)
        cx = np.iscomplexobj(v) or (inj is not None and np.iscomplexobj(inj))
        if cx:
            ir = None if inj is None else inj.real
            ii = None if inj is None else inj.imag
            re = self._monodromy_matvec_transposed(
                C0, steps, np.asarray(v).real, collect, ir)
            im = self._monodromy_matvec_transposed(
                C0, steps, np.asarray(v).imag, collect, ii)
            if collect:
                ## ⚠ The collected lists may be NESTED (a DIRK's `ts` holds
                ## per-stage solves per step), so a flat `a + 1j*b`
                ## multiplied a LIST by `1j`: `floquet_modes` under trbdf2
                ## raised "can't multiply sequence by non-int of type
                ## 'complex'" for as long as the path existed.  Recurse.
                return (re[0] + 1j * im[0], _cx_collect(re[1], im[1]),
                        _cx_collect(re[2], im[2]))
            return re + 1j * im
        v = v.astype(float)
        w1, w2 = v[:m].copy(), v[m:].copy()
        ts = []
        states = []

        ## The capacitance ring as the FORWARD pass saw it, so the reverse
        ## pass can consume it backwards.  Rebuilt rather than stored: it is
        ## `len(steps)` references to matrices that already exist.
        cs0, cs1, ring = [], [], list(C0)
        for _lu, C_new, _alphas, _b in steps:
            cs0.append(ring[0])
            cs1.append(ring[1])
            ring = [C_new, ring[0]]

        for j in range(len(steps) - 1, -1, -1):
            lu, _C_new, alphas, b = steps[j]
            if len(alphas) < 3:
                ## ⚠ A ONE-STEP COMPANION HAS NO `alphas[2]`, and this used
                ## to read it anyway -- an `IndexError` from inside a
                ## reverse loop, three frames from the cause, where the
                ## `b != 0` case one line below states its refusal
                ## plainly.  Unreachable through `solve`, which takes this
                ## path only for Gear-2, but PPV and adjoint noise both
                ## call the transposed replay directly and a one-step
                ## method is the obvious thing to try first.
                raise NotImplementedError(
                    'PSS: the transposed replay is derived for a two-step '
                    'companion (Gear-2) and this step has %d alpha '
                    'coefficients, so there is no second history term to '
                    'transpose. A one-step method needs its own reverse '
                    'recursion, not this one.' % len(alphas))
            if b:
                ## a `b != 0` companion carries `iq` in the state, so the
                ## step map is not the pair above.  Unreachable today --
                ## `_solves_history` is true only for Gear-2, whose `b` is
                ## zero -- and refused rather than silently transposing a
                ## different operator.
                raise NotImplementedError(
                    'PSS: the transposed replay is derived for a companion '
                    'with b = 0 (Gear-2), and this step reports b = %r. A '
                    'b != 0 method carries iq in the per-step state, so its '
                    'transpose is not the pair map this implements.' % (b,))
            t = lu.solve_transposed(w1)
            if t is None:
                raise NotImplementedError(
                    'PSS: this linear solver cannot solve transposed, so '
                    'the monodromy transpose cannot be replayed. Use '
                    'DenseSolver or SuperLUSolver.')
            if collect:
                ts.append(t)
            w1, w2 = (-alphas[1] * (cs0[j].T @ t) + w2,
                      -alphas[2] * (cs1[j].T @ t))
            if inj is not None:
                w1 = w1 + inj[j]
            if collect:
                ## ⚠ THE ADJOINT STATE AFTER THE STEP, which for a seed of
                ## the PPV IS the PPV at that time -- `Phi(T,s_j)^T v(T) =
                ## v(s_j)`.  Oscillator phase noise needs `v(t)` over the
                ## whole period, not `v(0)`, and it is already being
                ## computed here and thrown away.
                states.append(np.concatenate((w1.copy(), w2.copy())))
        if collect:
            ## reversed, so `ts[j]` lines up with `steps[j]` -- the reverse
            ## loop produced them last-first and a caller weighting them by
            ## a per-step time must not have to remember that
            ts.reverse()
            states.reverse()
            return np.concatenate((w1, w2)), ts, states
        return np.concatenate((w1, w2))

    ## How many deflated power iterations estimate the second multiplier.
    ## Convergence is at |lambda_3|/|lambda_2|, which is fast in the case
    ## that matters -- a lone slow node leaves everything below it tiny.
    ## ⚠ RETIRED 2026-09-03, kept as a name so the history reads.  The
    ## deflated power iteration this sized converged at
    ## `|lambda_3|/|lambda_2|` and lost three digits at a ratio of 1.065;
    ## Arnoldi replaced it at machine precision and fewer matvecs.
    PPV_DEFLATION_ITERS = 30
    ## Arnoldi basis size for the second-multiplier estimate.  Exact at
    ## `k = n`; below that it is a truncation and `lam2` is a lower bound
    ## ONLY for a normal `M` -- see `ppv`, where a circuit monodromy was
    ## measured to break the bound in both directions.  This is now the
    ## STARTING basis: `_ritz_second_multiplier` grows it until the pair's
    ## own Ritz residual certifies it.
    PPV_RITZ_BASIS = 12

    ## ⚠ THE GATE ON A TRUNCATED `lam2`: the SELECTED PAIR's Ritz residual,
    ## `|h_{k+1,k}| |y_i[last]| / ||y_i||`.  It needs no extra matvec -- both
    ## factors are already in the `H` this class forms -- and it is the one
    ## quantity that separates a converged pair from a leaked one.  The
    ## GMRES-style residual cannot: `_arnoldi_gmres`'s own note says a
    ## drifted basis "gives multipliers that are wrong in a way the residual
    ## cannot see", which is true of the SOLVE residual and false of the
    ## EIGENPAIR one.
    ##
    ## MEASURED (`_osc_with_ladder`, k = 12): <= 3.1e-07 at every `nslow`
    ## the truncated path gets right, and 2.8e-04 / 3.5e-04 / 2.1e-03 at the
    ## three it gets wrong -- and 1.5e-16 once `k` is large enough to be
    ## exact.  A peer session's independent sweep puts the two populations
    ## thirteen decades apart at the median, ⚠ TOUCHING at ~1e-5 (right 90th
    ## percentile 1.27e-05 against wrong 10th percentile 1.03e-05).  So the
    ## robust band is BELOW ~1e-6, and that -- not a magic 1e-8 -- is what
    ## this is set to.
    PPV_RITZ_RESIDUAL_TOL = 1e-6

    ## ⚠ A COST CEILING, NOT A CORRECTNESS THRESHOLD.  `k` doubles until the
    ## residual certifies or this is reached; overrunning it produces a
    ## WARNING and an uncertified number, never a silently wrong one, which
    ## is what makes an arbitrary-ish constant safe here.
    ##
    ## The size is set by what `k` has to reach: measured, `k` tracks the
    ## SLOW-MODE COUNT and not `n` (this tree's ladder cannot separate the
    ## two -- it sets `nslow = nladder` -- but a peer's synthetic can, and
    ## reports `k_min` flat under a doubling of `n` at fixed `nslow`).  The
    ## largest published case is Lai's 64-gated-capacitor DCO, so 128 leaves
    ## 2x headroom on it while costing 1/6 of that circuit's `n = 813`.
    PPV_RITZ_MAX_BASIS = 128
    ## Above this, the bordered extraction is losing digits AND the phase
    ## equation's instantaneous-response assumption is in doubt.
    PPV_SECOND_MULTIPLIER_WARN = 0.9

    def _ritz_second_multiplier(self, fp, kk):
        """`(lam2, residual)` from a `kk`-dimensional Arnoldi on `I - M`.

        Garcia, Romero & Acha (IEEE Trans. Power Systems 37(1), 2022): the
        Ritz values of `I - M` map back as `lambda = 1 - theta`.  Returns the
        selected pair's own RITZ RESIDUAL alongside it,
        `|h_{k+1,k}| |y_i[last]| / ||y_i||`, which costs nothing -- both
        factors are already in `H` and its eigenvectors.

        ⚠ `eig`, NOT `eigvals`.  The eigenVECTOR's last component is half the
        residual, so asking only for the values is what made this estimate
        uncheckable for as long as it was.

        ⚠ EVERY LOCAL IS UNDERSCORED ON PURPOSE, inherited from when this was
        inline in `ppv`: the first version used `q` for the Arnoldi start
        vector, silently overwriting `C(0) xdot(0)` -- returned as
        `info['q']` and consumed by two tests -- and the suite caught it as a
        shape mismatch three frames away.

        `residual` is `inf` when no Ritz value survives the deflation, so a
        caller that gates on it cannot read "nothing found" as "certified".
        """
        _n = fp.width
        kk = int(min(_n, kk))
        _rng = np.random.default_rng(12345)
        _q0 = _rng.standard_normal(_n)
        _q0 = _q0 / np.linalg.norm(_q0)
        _Qb = [_q0]
        _H = np.zeros((kk + 1, kk))
        for _j in range(kk):
            _wj = _Qb[_j] - np.asarray(fp.matvec(_Qb[_j]))
            for _i in range(_j + 1):
                _H[_i, _j] = float(_Qb[_i] @ _wj)
                _wj = _wj - _H[_i, _j] * _Qb[_i]
            _H[_j + 1, _j] = float(np.linalg.norm(_wj))
            if _H[_j + 1, _j] < 1e-13:
                ## ⚠ AN INVARIANT SUBSPACE: the basis closed on itself, so
                ## every Ritz pair in it is EXACT and the residual is zero by
                ## construction rather than by convergence.
                kk = _j + 1
                _theta, _Y = np.linalg.eig(_H[:kk, :kk])
                _lams = 1.0 - _theta
                _mask = np.abs(_lams - 1.0) > 1e-6
                if not _mask.any():
                    return 0.0, float('inf')
                return float(max(np.max(np.real(_lams)[_mask]), 0.0)), 0.0
            _Qb.append(_wj / _H[_j + 1, _j])
        _theta, _Y = np.linalg.eig(_H[:kk, :kk])
        _lams = 1.0 - _theta
        ## drop the unit root the border already accounts for
        _mask = np.abs(_lams - 1.0) > 1e-6
        if not _mask.any():
            return 0.0, float('inf')
        _idx = int(np.where(_mask)[0][int(np.argmax(np.real(_lams)[_mask]))])
        _lam2 = float(max(np.real(_lams)[_idx], 0.0))
        _res = (abs(_H[kk, kk - 1]) * abs(_Y[kk - 1, _idx])
                / max(float(np.linalg.norm(_Y[:, _idx])), 1e-300))
        return _lam2, float(_res)

    def _algebraic_adjoint_pattern(self, xf):
        """`(rows, cols)` — the ALGEBRAIC equations and the algebraic states.

        `rows` are the equations with no time-derivative (a zero row of
        `C`); `cols` are the state variables that appear under no
        derivative anywhere (a zero column).  For an index-1 DAE the two
        have the same count and the block between them is invertible,
        which is what makes the fill below well posed.

        Returns `([], [])` for a plain ODE, and the caller then does no
        per-sample work at all -- which is why this costs nothing on every
        fixture that has no algebraic row.
        """
        m = self.cir.n - 1
        Cr, = remove_row_col((np.asarray(self.cir.C(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Cr = np.asarray(Cr, dtype=float)
        rows = [i for i in range(m) if not np.any(Cr[i, :])]
        cols = [j for j in range(m) if not np.any(Cr[:, j])]
        return rows, cols

    def _equation_row_ppv(self, vblock, xf, rows, cols):
        """`v_1` — the adjoint contracted with an EQUATION-ROW input.

        ⚠⚠ THERE ARE TWO ADJOINT VECTORS AND CONFLATING THEM WAS A DEFECT.
        Demir 2000 puts both conventions on one page: a STATE initial
        condition contracts as `v_1^T(0) C(0) x(0)` (eq 41, WITH `C`), while
        an EQUATION-ROW input contracts as `v_1^T(s) b(s)` (eq 42, and the
        phase equation 44, BARE).  `ppv()` returns `C^T v_1`, which is the
        right object for a state perturbation and is documented as such.
        `CY` is an equation-row covariance -- a current injected into a KCL
        row -- so `diffusion_constant` and `colour_projection` need `v_1`,
        and this produces it.

        ⚠ ONE DEFECT, TWO SYMPTOMS, both measured.  `(C^T v_1)_j` is COLUMN
        `j` of `C` dotted with `v_1`, so on DIFFERENTIAL states `C^T`
        multiplies by the capacitance -- `diffusion_constant` was wrong by
        `C^2`, ratios 0.010003 / 1.000334 / 100.033536 over a 100x sweep --
        and on ALGEBRAIC states `C^T` ANNIHILATES, so `c` came back EXACTLY
        0.0 for an oscillator whose only noise was its series tank loss.
        Neither was visible because every fixture here uses `C = 1 F`.

        ⚠⚠ AND `C` IS NEVER INVERTED.  The solve is on `C[D, NZ]^T`, the
        block between the DIFFERENTIAL equations and the NON-algebraic
        states, which is square and invertible by construction; the
        singular `C` as a whole is not touched.  The algebraic entries come
        from the constraint below instead.
        """
        m = self.cir.n - 1
        Cr, = remove_row_col((np.asarray(self.cir.C(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Cr = np.asarray(Cr, dtype=float)
        diff = [i for i in range(m) if i not in rows]
        nz = [j for j in range(m) if j not in cols]
        out = np.zeros(m, dtype=float)
        blkC = Cr[np.ix_(diff, nz)].T
        try:
            out[diff] = np.linalg.solve(blkC, np.asarray(vblock)[nz])
        except np.linalg.LinAlgError:
            warnings.warn(
                'PSS.ppv: the differential block C[D, NZ] is singular, so '
                'the equation-row adjoint cannot be recovered; falling back '
                'to the state-perturbation vector, which is wrong by a '
                'factor of the capacitance in any CY contraction.',
                RuntimeWarning, stacklevel=2)
            return np.array(vblock, dtype=float, copy=True)
        if not rows:
            return out
        return self._algebraic_adjoint_fill(out, xf, rows, cols)

    def _algebraic_adjoint_fill(self, vblock, xf, rows, cols):
        """Fill the adjoint's ALGEBRAIC entries, which are SLAVED, not free.

        ⚠ THE REPLAY LEAVES THEM AT ZERO AND ZERO IS NOT THEIR VALUE.  The
        PPV entry for a row IS the phase sensitivity to a perturbation
        entering that row, and an algebraic row's perturbation reaches the
        dynamics through the CONSTRAINT rather than through its own row.
        On a tank with series loss, eliminating `v_x = r (i_L + b)` puts
        `-r b` into the inductor's equation, so the sensitivity to node
        `x` is `r` times the branch row's -- and a noise current landing
        there was being contracted against a structural zero, which made
        `diffusion_constant` return EXACTLY 0.0 for an oscillator whose
        only noise was its tank loss.  Measured against three independent
        references; see the roadmap's section 0d.

        The adjoint equation's ALGEBRAIC-STATE columns are what determines
        them.  For a column `j` with no `C` entry the equation carries no
        derivative, so it reads `sum_i G_ij v_i = 0`, and splitting `i`
        into algebraic and differential rows gives

            v_A  =  (G[A, Z]^T)^-1 G[D, Z]^T v_D

        ⚠ THE MAGNITUDE IS STRUCTURAL AND THE SIGN IS MEASURED, and saying
        which is which is the point.  On a RESISTIVE DIVIDER between the
        inductor and ground the two algebraic nodes fold into the branch
        row with coefficients `(r1 + r2)` and `r2`, so their entries must
        stand in a ratio the topology fixes -- measured `10.000000` against
        a chosen `10.000000`, on a fixture built so the single-resistor
        degeneracy cannot hide a mistake.  ⚠ THE SINGLE SERIES RESISTOR
        CANNOT SETTLE THIS: there `|integral v_0|` and `|r integral
        v_branch|` agree to 1.5e-4, so BOTH SIGNS FIT and an agreement
        there is no evidence at all.

        The overall sign is then fixed by requiring algebraic and
        differential rows to share ONE convention -- `dT/dA = +integral
        v_j` -- and measured on the divider: with it, the DC-injection
        probe reads +0.9999849 (node v, differential), +0.9999982 and
        +0.9998744 (the two algebraic nodes); with the sign the naive
        derivation gives, the last two come back NEGATIVE.

        Returns `vblock` unchanged, with a warning, when the structure is
        not index-1: `len(rows) != len(cols)` or a singular block.  That
        case is section B4's, and guessing at it would be worse than
        leaving a known zero.
        """
        if not rows:
            return vblock
        m = self.cir.n - 1
        if len(rows) != len(cols):
            warnings.warn(
                'PSS.ppv: %d algebraic equations against %d algebraic '
                'states, so the adjoint\'s algebraic entries are not '
                'determined by a square solve -- this is an index > 1 '
                'structure (roadmap B4). They are left at zero, and noise '
                'entering those rows will be UNDER-COUNTED.'
                % (len(rows), len(cols)), RuntimeWarning, stacklevel=2)
            return vblock
        Gr, = remove_row_col((np.asarray(self.cir.G(xf), dtype=float),),
                             self.irefnode, self.toolkit)
        Gr = np.asarray(Gr, dtype=float)
        diff = [i for i in range(m) if i not in rows]
        blk = Gr[np.ix_(rows, cols)].T
        ## ⚠ THE SIGN IS NOW THE DERIVED ONE, because this acts on `v_1`.
        ## It was FLIPPED while this fill acted on `C^T v_1`: on that
        ## fixture `C = diag(1, 0, -L)`, the INDUCTOR BRANCH ROW CARRIES
        ## `-L`, and the term below is dominated by that branch -- so
        ## reading `C^T v_1` as `v_1` negated it.  The derivation and the
        ## measurement were describing DIFFERENT VECTORS and both were
        ## right.  Pinned by the eq (24) constraint, which returns
        ## 0.0000e+00 exactly here and 1.9870e+00 for either alternative.
        rhs = -(Gr[np.ix_(diff, cols)].T @ np.asarray(vblock)[diff])
        try:
            va = np.linalg.solve(blk, rhs)
        except np.linalg.LinAlgError:
            warnings.warn(
                'PSS.ppv: the algebraic block G[A, Z] is singular, so the '
                'adjoint\'s algebraic entries cannot be recovered; they '
                'are left at zero and noise entering those rows will be '
                'UNDER-COUNTED.', RuntimeWarning, stacklevel=2)
            return vblock
        out = np.array(vblock, dtype=float, copy=True)
        out[rows] = va
        return out

    def _ppv_propagate(self, fp, v, m, xdot, alg_rows, alg_cols):
        """The pair-consistent SECOND-ORDER propagation of an anchor vector
        `v` over the period (the block `ppv()` applies to its null vector,
        lifted 2026-09-08 so `frequency_aware_ppv` can run it on a COMPLEX
        anchor).  Returns `(states, states_pair, ts, Xf)`: `states` the
        per-step state-space samples (`C^T v`, second order, rescaled by
        `v . xdot = 1` at the first sample), `states_pair` the raw
        pair-space replay, `ts` the transposed per-step states, `Xf` the
        orbit.  Only the `solved_history` (LMM) kind carries the
        correction; the others return the replay as it is."""
        _dt = complex if np.iscomplexobj(v) else float
        _alg_rows, _alg_cols = alg_rows, alg_cols
        _end, _ts, states = fp.matvec_transposed(v, collect=True)
        states_pair = [np.array(st, dtype=_dt, copy=True) for st in states]
        _Xf = np.asarray(self.waveform[1], dtype=float)
        ## ⚠⚠ THE PAIR'S FIRST BLOCK IS NOT THE PPV, AND THE ERROR IS FIRST
        ## ORDER AND GROWS WITH Q.  For Gear-2 the adjoint state is the pair
        ## `(w1, w2) = (dphi/dx_k, dphi/dx_{k-1})`, and `w1` alone is the
        ## response to a perturbation of `x_k` WITH `x_{k-1}` HELD -- an
        ## inconsistent history, which the two-step method resolves through
        ## its parasitic root.  A physical state perturbation moves both:
        ## `dx_{k-1} = Phi(t_{k-1}, t_k) dx_k`, so the phase functional is
        ##
        ##     v(t_k) = w1 + Phi(t_{k-1}, t_k)^T w2,   Phi ~ I - h J + O(h^2)
        ##
        ## and with `w2 = C_{k-1}^T z` (`z = -a2 t_k`, exact by the
        ## recursion) that is `w1 + (C_{k-1} + h G)^T z` -- no inverse of
        ## `C`, so it holds for a DAE.  Equivalently `w1` is orthogonal to
        ## the amplitude eigenvector's first block, which is the true
        ## amplitude direction ROTATED by `O(h)`; `v . xdot = 1` then
        ## amplifies that rotation by `|v||xdot|`, the near-cancellation a
        ## non-isochronous oscillator has (its PPV grows with `Q_lambda`).
        ## MEASURED against the exact continuous adjoint (DOP853 at 1e-12,
        ## no shooting code in the reference) on `vdp + 0.3 u^2`, whose
        ## `c` is 100x van der Pol's: the first block gave `c` 16.6 / 8.0 /
        ## 3.9 / 1.9 / 1.0% high at 400..6400 points -- clean first order
        ## -- and violated `v(t) . xdot(t) = 1` along the orbit by 12%
        ## (std 2.7e-2).  This contraction holds the invariant to 8e-5 and
        ## gives `c` to 1.8e-3 at 400 and 8e-5 at 1600, second order.  On
        ## van der Pol both agree to 1e-4: the two rows are in quadrature
        ## there, so the rotation averaged out of `<v^2>` -- the fixture
        ## shared the claim's assumption (failure shape 0b), and `pnoise`,
        ## which contracts in PAIR space, was right all along and 14%
        ## below `c` on the fixture that could see it.
        ## The seed's scale is `w1(0) . xdot = 1`; the consistent object
        ## is renormalised ONCE by its own `v(0) . xdot`, which is why the
        ## per-step invariant is the test and not the definition.
        if (fp.kind == 'solved_history' and len(states) > 0
                and len(states[0]) == 2 * m):
            _cs1, _ring = [], list(fp.opening)
            for _lu, _Cn, _al, _b in fp.steps:
                _cs1.append(_ring[1])
                _ring = [_Cn, _ring[0]]
            _hs = np.diff(np.asarray(fp.times, dtype=float))
            _vphys = []
            for _j, st in enumerate(states):
                _lu, _Cn, _al, _b = fp.steps[_j]
                _z = -_al[2] * np.asarray(_ts[_j], dtype=_dt)
                ## ⚠ DIFFERENTIAL ROWS ONLY.  `w2 = C^T z` does not see the
                ## algebraic rows of `z` (their rows of `C` are zero), so
                ## the decomposition is non-unique there, and those
                ## multipliers are O(1/h): `h G^T z` would carry an O(1)
                ## component along the constraint normal into the
                ## differential entries.  The consistent propagation
                ## `C_D dx_{k-1} = (C_D + h G_D) dx_k` involves only the
                ## differential equations, which is the choice that makes
                ## it unique.  Measured: with the algebraic rows in, a DC
                ## injection probe on a series-loss tank flipped sign.
                if _alg_rows:
                    _z[np.asarray(_alg_rows, dtype=int)] = 0.0
                _xj = _Xf[:, _j if _j < _Xf.shape[1] else -1]
                _Gj, = remove_row_col(
                    (np.asarray(self.cir.G(_xj), dtype=float),),
                    self.irefnode, self.toolkit)
                _Gj = np.asarray(_Gj, dtype=float)
                _Cj = np.asarray(_cs1[_j], dtype=float)
                if _alg_rows:
                    ## ⚠ ON A DAE THE ALGEBRAIC STATE IS SLAVED, AND ITS
                    ## COUPLING INTO THE DIFFERENTIAL PROPAGATION IS O(h).
                    ## A consistent perturbation propagates as
                    ## `C_D dx_{k-1} = (C_D + h G_red) dx_k` on the
                    ## differential states, with `G_red` the Schur
                    ## complement `G[D,NZ] - G[D,Z] G[A,Z]^-1 G[A,NZ]`.
                    ## With the full `G` instead, the series-loss tank had
                    ## `c` 0.6 / 0.3 / 0.15% low at 240/480/960 points
                    ## (first order) and the invariant drifting at 1.1e-3;
                    ## with the complement `c` is 2.7e-4 / 7e-5 / 2e-5 from
                    ## the exact reduced-ODE value and the drift 4e-4 /
                    ## 1.1e-4 / 2.7e-5 -- second order (found through the
                    ## review session's linear-DAE partition, 2026-09-05).
                    _A = np.asarray(_alg_rows, dtype=int)
                    _Zc = np.asarray(_alg_cols, dtype=int)
                    _D = np.array([i for i in range(m) if i not in _alg_rows],
                                  dtype=int)
                    _NZ = np.array([j for j in range(m) if j not in _alg_cols],
                                   dtype=int)
                    ## ⚠ AND `G[A,Z]` NONSINGULAR *IS* THE INDEX-1 CONDITION.
                    ## At index >= 2 (an L-I cutset, a C-V loop) it is
                    ## singular by definition and the algebraic variables
                    ## come from a differentiation, not a solve; the
                    ## complement does not exist.  Same shape as the fill
                    ## above: warn once with the reason and fall back to
                    ## the full-`G` propagation, which is then first order.
                    ## (Boundary named by the review session from the
                    ## pencil: `eig(-G_red, C[D,NZ])` equals the finite
                    ## generalised eigenvalues of `(C, G)` to 1e-12 on the
                    ## series-loss tank, and the reduction is undefined on
                    ## `li_plus_rc` and `cv_plus_rc`.)
                    try:
                        if len(_alg_rows) != len(_alg_cols):
                            raise np.linalg.LinAlgError('not square')
                        _Gred = (_Gj[np.ix_(_D, _NZ)]
                                 - _Gj[np.ix_(_D, _Zc)] @ np.linalg.solve(
                                     _Gj[np.ix_(_A, _Zc)], _Gj[np.ix_(_A, _NZ)]))
                    except np.linalg.LinAlgError:
                        if _j == 0:
                            warnings.warn(
                                'PSS.ppv: the algebraic block G[A,Z] is '
                                'singular (index > 1: an L-I cutset or a '
                                'C-V loop), so the pair-consistent '
                                'propagation cannot eliminate the algebraic '
                                'state and falls back to the full G -- the '
                                'PPV samples can then be FIRST order in the '
                                'step, as they are for the algebraic fill. '
                                'Priced on the fixture that can see the '
                                'dropped term (a non-isochronous core with '
                                'its inductor split, and the same core with '
                                'a C-V loop through a bias rail): within 3e-4 '
                                'of the index-1 object and second order on '
                                'both, so at index 2 this fallback is the '
                                'whole answer.',
                                RuntimeWarning, stacklevel=2)
                        _Gred = _Gj[np.ix_(_D, _NZ)]
                    _corr = np.zeros(m, dtype=_dt)
                    _corr[_NZ] = (_Cj[np.ix_(_D, _NZ)]
                                  + _hs[_j] * _Gred).T @ _z[_D]
                    _vp_j = st[:m] + _corr
                else:
                    _vp_j = st[:m] + (_Cj + _hs[_j] * _Gj).T @ _z
                ## ⚠ AND ZERO ON THE ALGEBRAIC COLUMNS, as `C^T v_1` is:
                ## the state functional contracts a perturbation ON the
                ## constraint manifold, whose algebraic components are
                ## slaved, and `h G^T z` would otherwise leave 4e-3 there
                ## (caught by the full suite's Demir-(24) gate).  The
                ## equation-row conversion never reads these entries.
                if _alg_cols:
                    _vp_j[np.asarray(_alg_cols, dtype=int)] = 0.0
                _vphys.append(_vp_j)
            _scale = (_vphys[0] @ xdot) if _dt is complex else float(_vphys[0] @ xdot)
            if _scale == 0.0:
                raise ValueError(
                    'PSS.ppv: the pair-consistent adjoint is orthogonal to '
                    'the orbit tangent at t = 0.')
            ## ⚠ AND ITS DC CONTENT IS THE CONSISTENT OBJECT'S TOO -- taking
            ## the mean from the raw block was TRIED AND MEASURED WRONG.
            ## The raw block's orbit integral reproduces a same-grid
            ## DC-injection probe to 1e-5 on the divider fixture (node row,
            ## true mean 4e-6 |v|), where the consistent object's O(h^2)
            ## pointwise errors leave an absolute floor of ~1e-5 |v| --
            ## the wrong sign at 480 points.  But on the bias-sensitive
            ## fixture's INDUCTOR row (a DC voltage in series with L, true
            ## dT/dV = 16.20 by a second-order re-solve) the raw block
            ## reads 17.49 / 16.83 / 16.51 at 400/800/1600 -- first order,
            ## 8% off -- while the consistent object holds `v . xdot = 1`
            ## to 8e-5 along the orbit, which pins its mean in EVERY row to
            ## ~1e-5 |v|.  Stitching the raw mean in broke that invariant
            ## by +-0.3.  So the raw block's DC exactness is row- or
            ## fixture-specific (mechanism open, recorded in the roadmap),
            ## and `samples` is one object, second order everywhere, with a
            ## ~1e-5 |v| absolute floor on its mean.  The raw pair is kept
            ## as `samples_pair` for the structural gates that live on its
            ## discrete identities.
            states = [np.concatenate((vp / _scale, st[m:] / _scale))
                      for vp, st in zip(_vphys, states)]
        return states, states_pair, _ts, _Xf

    def ppv(self, tol=None):
        """The perturbation projection vector at `t = 0` (Demir & Roychowdhury).

        Returns `(v, info)`.  `v` is the pair-space left null vector of
        `I - M`, normalised so that `v . xdot(0) = 1` -- see the note on the
        normalisation below, which was MEASURED rather than transcribed.
        The phase shift caused by a state perturbation `delta` at `t = 0` is
        then `v[:m] . delta`.  `info` carries both border residuals, the
        null residual, `q` and the scaled tangent.

        ⚠ AN AUGMENTED SOLVE, NOT AN EIGENVECTOR -- AND IT IS THE FIX FOR A
        NAMED FAILURE OUR OWN FIXTURES SIT INSIDE.  Demir &
        Sangiovanni-Vincentelli, 1998 (the book, read firsthand by the docs
        session 2026-09-04), report the eigenvector route BREAKING on a
        high-Q circuit, with a table of the crowded eigenvalues (their
        Table 6.4): "Since this circuit is a high-Q one, Phi(T,0) has
        eigenvalues with magnitudes close to 1 other than the one which is
        supposed to be equal to 1 ... Because of numerical errors, we can
        not identify the eigenvalue that is supposed to be equal to 1 ...
        so it is not feasible to identify the correct" one.  Not
        ill-conditioned there -- INFEASIBLE.  That reported failure is the
        stated motivation for the single-solve method two years later, so
        the lineage is firsthand end to end: 1998 selection fails at high
        Q; 2000 the single linear solve; 2001 "particularly useful for
        high-Q oscillators"; 2003 the fuller procedure.  ⚠ It also changes
        what the second-multiplier warning below MEANS: not "this result is
        degrading" but "you are in the regime this method was invented to
        escape".  And the same book's eq (6.72), `|exp(eta_i)| << 1`, is
        the closed-form variance's validity condition -- the book says it
        "is satisfied for 'most' oscillator circuits" and defers the rest
        to the high-Q discussion above, so that condition and the crowding
        are ONE condition seen from two sides.

        ⚠ ATTRIBUTION CORRECTED
        2026-09-04: the idea ORIGINATES in Demir, Long & Roychowdhury,
        ICCAD 2000 ("Computing Phase Noise Eigenfunctions Directly from
        Steady-State Jacobian Matrices" -- "a single linear solution of the
        oscillator's ... steady-state Jacobian matrix ... dispenses with the
        need to select the correct one eigenfunction"), with the 2001
        companion carrying it to HB/shooting matrices and noting the
        selection heuristic is worst "for high-Q oscillators".  The 2003
        paper is the fuller procedure and the source of the quote below,
        not the origin.  Demir's IJCTA 2000 method SELECTED the
        right eigenvector by its inner product against `C(0) xdot(0)` --
        measured 0.2 against 1e-5, 1e-7, 2e-5 on a Colpitts.  His 2003
        paper rejects that: "no guarantee that any of the candidate
        eigenvectors will be appreciably more orthonormal than the others,
        leading to a potential breakdown."  The same vector then changes
        role -- it becomes the BORDER, so the candidate is unique and
        nothing is selected:

            [ I - M^T   q ] [v]     [0]
            [   q^T     0 ] [y]  =  [1]

        This matters here specifically.  `_spectral_report`'s eigenvector
        split was measured labelling a parasitic root physical at ~2 points
        per cycle, and multipliers cluster near 1 on exactly the high-Q
        oscillators a PPV is wanted for -- four independent sources say so.
        A bordered solve does not care how close the other multipliers are;
        it never has to tell them apart.

        ⚠ AND THAT IS THE METHOD'S STATED DESIGN DRIVER, not a lucky
        property of it.  Demir, Long & Roychowdhury (ICCAD 2000), who
        introduced the single-solve route: it is "especially advantageous
        for HIGH-Q OSCILLATORS, MONODROMY MATRICES OF WHICH OFTEN HAVE MANY
        EIGENVALUES CLOSE TO 1 THAT ARE NUMERICALLY INDISTINGUISHABLE from
        the oscillatory [unit eigenvalue]", and "a key advantage is that it
        DISPENSES WITH THE NEED TO SELECT THE CORRECT ONE-EIGENFUNCTION
        from amongst a potentially large set of choices".  So the hardest
        case in this codebase is the case the method was aimed at.

        ⚠ WHICH DOES NOT RETRACT THE SLOW-NODE BOUNDARY BELOW, and keeping
        the two apart is the point.  Near-degenerate multipliers make
        EIGENANALYSIS ILL-POSED -- there is no fact of the matter about
        which eigenfunction is the PPV -- while they make this bordered
        solve merely ILL-CONDITIONED, which is measured above and warned
        about.  Ill-conditioned beats ill-posed, and neither is the same as
        the PHASE EQUATION's own limit, which is about the response being
        treated as instantaneous and is not an extraction question at all.
        Three separate things that a "high-Q oscillators are hard" summary
        would blur into one.

        ⚠ THE QUADRATIC RUNG BELOW WOULD NOT FIX THE SLOW-NODE BOUNDARY,
        which is the obvious hope and is wrong.  TWO INDEPENDENT
        approximations are in play: the LINEAR ISOCHRON one is in the
        perturbation's AMPLITUDE -- it treats isochrons as flat
        hyperplanes, and the quadratic rung adds their curvature -- while
        the INSTANTANEOUS-RESPONSE one is in the DYNAMICS, ignoring the
        bandwidth between injection point and core.  Slow nodes are the
        second.  Noise is small by construction, so the linear term
        dominates it by definition; the quadratic rung would earn its cost
        on LARGE perturbations -- injection locking, big supply or
        substrate interferers -- not on phase noise.

        ⚠ AND THE PPV IS ONE RUNG ON A LADDER, worth knowing before it is
        mistaken for exact.  Suvak & Demir (TCAD 2011) place it: an EXACT
        phase equation exists and is "practically unusable"; the PPV
        equation is its LINEAR isochron approximation; a QUADRATIC one is
        more accurate.  Isochrons are the geometric form of asymptotic
        phase, so an oscillator without asymptotic phase is one whose
        isochrons do not exist -- the same fact as the Floquet condition,
        seen in the geometry.  Computing exact isochrons is exponential, so
        the only live question is which local approximation is affordable.

        `c` -- the diffusion constant this vector feeds -- has the
        designer-facing reading "JITTER PER SECOND".

        ⚠ THREE NAMES FOR THIS OBJECT, AND ONE NEAR-MISS THAT IS NOT IT.
        The PPV, Kaertner's adjoint LPTV impulse response, the PRC of
        mathematical biology, and Hajimiri's NUMERICAL ISF are the same
        thing.  His CLOSED-FORM ISF is NOT: it is the normalised tangent,
        and the difference is not a scale factor but a SIGN -- for a noise
        impulse at one point in the cycle "the closed-form ISF predicts a
        POSITIVE phase change, whereas in fact the correct phase change is
        in the NEGATIVE direction and of a different magnitude".  It also
        does not scale with the perturbation, where the PPV does.

        ⚠ AND `xdot` IS NOT A CHEAP SUBSTITUTE FOR IT -- "time-shifts and
        amplitudes are both different ... the two waveforms scale in
        OPPOSITE DIRECTIONS with respect to the RC time constant".  Nothing
        here offers it as one: `xdot` appears only as the NORMALISATION
        (`v . xdot = 1`, the PPV's defining property), as the border `q =
        C(0) xdot(0)`, and in the record above of the SELECTION heuristic
        that was rejected.  Stated because the substitution is a documented
        point of common confusion, and the failure would be a sign error
        rather than a visible one.

        ⚠ `y` COMING BACK ZERO IS A FREE CORRECTNESS CHECK, and it is not
        decoration.  With a zero first block on the right-hand side,
        `(I - M^T) v + y q = 0` forces `y q = 0`, so a nonzero `y` means the
        border absorbed a residual the null space should have taken -- the
        computed `v` is not in the null space.  Measured on van der Pol:
        1.4e-11.

        ⚠ ITS VALIDITY BOUNDARY IS SLOW NODES, AND THE VAN DER POL GATE
        CANNOT SEE IT.  The phase equation `alpha' = v_1^T(t+alpha) b(t)`
        treats the oscillator's frequency response as INSTANTANEOUS; the
        truth is a convolution, and the PPV form is what you get by
        assuming the kernel is `v_1(t) delta(t - tau)`.  Real circuits have
        finite bandwidth, so a slow node FILTERS the noise of devices near
        it, the PPV cannot see the filtering, and phase noise is
        OVER-ESTIMATED.  Lai (Cadence) is explicit that better extraction
        does not help: "although the PPV can be extracted correctly, the
        oscillator noise analysis is still inaccurate: the phase noise is
        always over-estimated."

        ⚠ AND HE NAMES THIS TEST'S OWN REGIME AS THE BLIND SPOT: "the
        phase equation was verified to be correct in many previous works
        ... because it was evaluated on SMALL, SIMPLE OSCILLATORS, and
        perturbations were applied to OSCILLATOR CORES.  Since oscillator
        cores have very wide bandwidth, ignoring the dynamics may not
        compromise the macromodelling accuracy very much."  That is
        `test_the_ppv_predicts_a_phase_shift_the_oscillator_actually_has`
        exactly -- van der Pol, perturbed at its core.  It passes whether
        or not this failure mode is present, so it establishes that the
        extraction and normalisation are right and says NOTHING about the
        model's range.  The fix is a frequency-aware PPV, which is this
        same bordered system at nonzero `w_s` -- the classical PPV is its
        DC point, and `PAC` already solves at nonzero frequency.  Not
        built.  VERIFIED at the source (docs session, 2026-09-08): Lai
        2008 eq. (24) at `w_s = 0` "is the augmented PPV extraction
        equation (6) and (7)", verbatim; two scope limits: (24) is a
        NEAR-DC approximation of (23) (the AC Toeplitz columns dropped,
        "if we are only interested in ... w_s close to DC"), and (23)
        "is very difficult to solve using iterative solvers (such as
        GMRES)" because the border degrades the block preconditioner.
        ⚠ AND `_vdp_with_slow_node` CANNOT SHOW THE EFFECT AT ANY tau/T:
        its slow node couples through `Rs = 1e6` against a tank impedance
        of 1, so its PPV entry is 7.3e-6 of the core's, flat over
        tau/T = 1e2..1e6 while lambda_2 moves four decades -- it tests
        conditioning, not the filtering of noise that REACHES the phase.
        The fixture owed is a slow node IN the phase path (small Rs, large
        Cs), with tau/T and coupling as separate knobs.  BUILT AND MEASURED
        the same day: with an asymmetric core AND tank loss (G_0 != 0 needs
        both) the slow node's PPV entry is DC-dominated (|G0|/|G1| = 30) and
        the Lorentzian over-states a source behind it by 1000x at 0.1 f0,
        predicted to four digits from this PPV's harmonics and the RC
        filter -- so `c` from this PPV is right and the SHAPE above
        T/(2 pi tau) is what the frequency-aware PPV corrects; see
        `oscillator_spectrum`.  A2's "gated at tau/T = 10" was a Monte
        Carlo of `c`, which cannot see it.

        ⚠ `q` IS EXACT, NOT DIFFERENCED.  `q = C(0) xdot(0)` looks like it
        needs the orbit's tangent, and differencing the waveform for it
        would be O(h) at best.  The DAE gives it directly: `dq/dt + i(x) +
        u(t) = 0` and `dq/dt = C xdot`, so `q = -(i(x_0) + u(0))` -- two
        evaluations at the converged solution, no derivative anywhere.
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.ppv(tol)
        import scipy.sparse.linalg as spla
        fp = self.factored_period()
        ## ⚠ NO LONGER GEAR-ONLY (B8). The refusal that stood here said the
        ## transposed replay was "implemented for the solved-history map
        ## only"; since `_monodromy_matvec_transposed_plain` shipped that
        ## sentence is false, and every call below goes through
        ## `fp.matvec_transposed`/`fp.matvec` so the map's kind is the
        ## dispatcher's business rather than this method's.
        if not self.autonomous:
            raise ValueError(
                'PSS.ppv: a perturbation projection vector describes the '
                'phase of a FREE-RUNNING oscillator. This circuit is '
                'driven, so its phase is the source\'s and there is no '
                'unit Floquet multiplier to project onto.')

        m = self.cir.n - 1
        n = fp.width
        irn = self.irefnode
        x0r = np.asarray(self._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        qf = -(np.asarray(self.cir.i(x0f)).ravel()
               + np.asarray(self.cir.u(0.0,
                                       analysis=self.par.analysis)).ravel())
        q = np.delete(np.asarray(qf, dtype=float), irn)
        qp = np.concatenate((q, np.zeros(n - m)))
        nq = float(np.linalg.norm(qp))
        if nq == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) xdot(0) is zero, so the orbit has no tangent '
                'at t=0 and there is nothing to normalise against. That '
                'should not happen on a converged limit cycle.')

        def _mv(z):
            z = np.asarray(z)
            v_, y_ = z[:n], z[n]
            top = v_ - fp.matvec_transposed(v_) + y_ * qp
            return np.concatenate((top, [float(qp @ v_)]))

        rtol = max(self.par.reltol * 1e-2 if tol is None else tol, 1e-14)
        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=float)
        rhs = np.zeros(n + 1)
        rhs[n] = 1.0
        ## ⚠ JUDGED BY ITS RESIDUAL, NOT BY A STATUS CODE.  SciPy returns
        ## `info = 4` on a HAPPY breakdown -- the Krylov space exhausted
        ## because the answer is exact -- and these bordered operators are
        ## small enough to hit that routinely.  `_arnoldi_gmres` detects
        ## the breakdown where it happens and returns the exact answer.
        z, relres, _Ha, _ka = _arnoldi_gmres(
            _mv, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relres > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PSS.ppv: the augmented solve did not converge (relative '
                'residual %.3e). The border is `q` itself; if the orbit '
                'tangent is nearly orthogonal to the null direction the '
                'bordering is poor.' % relres)
        v, y = z[:n], float(z[n])
        resid = float(np.linalg.norm(v - fp.matvec_transposed(v)))

        ## ⚠ THE SCALE NEEDS THE TANGENT, so the RIGHT null vector is solved
        ## for too -- by the same bordering, not by an eigendecomposition,
        ## for the same reason: on a high-Q oscillator the other multipliers
        ## crowd 1 and no selection among candidates is reliable.
        def _mvf(zz):
            zz = np.asarray(zz)
            u_, yy = zz[:n], zz[n]
            top = u_ - fp.matvec(u_) + yy * qp
            return np.concatenate((top, [float(qp @ u_)]))

        Af = spla.LinearOperator((n + 1, n + 1), matvec=_mvf, dtype=float)
        zf, relresf, _Hf, _kf = _arnoldi_gmres(
            _mvf, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relresf > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PSS.ppv: the tangent solve did not converge (relative '
                'residual %.3e).' % relresf)
        u, yf = zf[:n], float(zf[n])

        ## `u` is the tangent's DIRECTION; its scale comes from `C u = q`,
        ## which is the definition of `q` read backwards.  Least squares
        ## because `C` is singular for a DAE and only its range is
        ## determined.
        x0red = np.asarray(self.cir.C(x0f))
        Cm = np.delete(np.delete(x0red, irn, 0), irn, 1)
        Cu = np.asarray(Cm, dtype=float) @ u[:m]
        denom = float(Cu @ Cu)
        if denom == 0.0:
            raise ValueError(
                'PSS.ppv: C(0) annihilates the orbit tangent, so its scale '
                'is not determined by C u = q.')
        xdot = u[:m] * (float(q @ Cu) / denom)

        ## ⚠ NORMALISED BY `v . xdot = 1`, WHICH IS NOT WHAT `v . q = 1`
        ## GIVES, and the difference is not cosmetic: on `_vdp_ppv(400)`
        ## `v . xdot = 1.0` against `v . q = -1.0696` -- a factor 2.07 AND
        ## the opposite sign (an earlier version of this note said "7%",
        ## which was the |v . q| - 1 residual and not the error; corrected
        ## by the review session's audit, 2026-09-04).  The defining
        ## property is that displacing the state
        ## ALONG the orbit by `eps xdot` advances the phase by `eps`, so
        ## `v . xdot = 1` is the normalisation a state perturbation sees.
        ## Demir's Remark 3.1 reads `v_1^T C u_1 = 1`; the vector this
        ## bordered solve returns behaves as `C^T v_1` -- it is contracted
        ## with a state perturbation directly -- so the two statements agree
        ## about different objects.
        ##
        ## ⚠ AND THAT IS NOT A QUIRK OF THIS FORMULATION, WHICH THIS
        ## COMMENT USED TO IMPLY.  The conserved pairing propagates to
        ## `M^T (C(0)^T v_1) = C(0)^T v_1`, so the left eigenvector of the
        ## STATE-SPACE monodromy simply IS `C(0)^T v_1` -- for ANY `C`,
        ## symmetric or not, and whatever the augmentation.  MEASURED on a
        ## limit cycle with a constant NON-SYMMETRIC `C`: alignment with
        ## `C(0)^T v_1` is 1.000000000000 against 0.9657 for `v_1` itself,
        ## and bordering with `xdot(0)` reproduces Demir's normalisation
        ## exactly while bordering with `C(0) xdot(0)` gives 0.805.
        ## (Derived and measured by the docs session, 2026-09-04.)  ⚠ TREATING THEM AS THE SAME OBJECT WAS
        ## MEASURED WRONG: predicting a state jump's phase shift as
        ## `v^T C delta` gives residuals of 0.36/0.40/0.42 that GROW with
        ## refinement and per-direction ratios scattering from -0.44 to
        ## 28.7, while `v . delta` converges at O(h).
        ## ⚠ THE ALGEBRAIC ENTRIES ARE FILLED *AFTER* THIS, AND THAT IS A
        ## DECISION RATHER THAN AN ORDERING ACCIDENT.  Filling first was
        ## tried and MEASURED WORSE: the DC-injection probe went from
        ## 0.9999849 to 0.9992364, because `v` here is also the REPLAY'S
        ## SEED and the algebraic components are SLAVED -- propagating them
        ## through the step map corrupts the differential ones.
        ##
        ## ⚠ AND THE NORMALISATION SHOULD NOT SEE THEM EITHER.  `v . xdot`
        ## is about a STATE perturbation, and a state perturbation of a DAE
        ## lies ON the constraint manifold: its algebraic components are
        ## determined by its differential ones, not free.  The algebraic
        ## entries of `v` answer a different question -- the sensitivity to
        ## a perturbation of an EQUATION ROW, which is what a noise current
        ## injected into an algebraic KCL row is.  So this line is
        ## unchanged, and every PPV number on every circuit is
        ## bit-for-bit what it was.
        _alg_rows, _alg_cols = self._algebraic_adjoint_pattern(x0f)
        vx = float(v[:m] @ xdot)
        if vx == 0.0:
            raise ValueError(
                'PSS.ppv: the null vector is orthogonal to the orbit '
                'tangent, so no normalisation makes it a phase projector.')
        v = v / vx
        ## ⚠ THE PPV OVER THE PERIOD, not just at `t = 0`, because that is
        ## what an oscillator noise calculation needs: Demir's diffusion
        ## constant is `c = (1/T) integral v_1^T(t) B(t) B^T(t) v_1(t) dt`,
        ## an integral over the orbit.  `Phi(T,s)^T v(T) = v(s)`, and the
        ## reverse replay computes exactly that sequence on its way to the
        ## answer -- it was being discarded.
        states, states_pair, _ts, _Xf = self._ppv_propagate(fp, v, m, xdot, _alg_rows, _alg_cols)
        ## ⚠ AND FILL EVERY SAMPLE TOO, at ITS OWN operating point, because
        ## `G` is state-dependent and the algebraic entries are a pointwise
        ## function of the differential ones.  Done here rather than by
        ## seeding the replay: these components are SLAVED, so there is
        ## nothing to propagate, and post-processing leaves the validated
        ## step map untouched.  The pair's SECOND block is the history term
        ## and is deliberately not filled -- `v(t)` is the first block.
        ## ⚠⚠ THE EQUATION-ROW ADJOINT IS A SECOND OBJECT, NOT A CORRECTION
        ## TO THE FIRST.  `states` and `v` stay exactly what they were --
        ## `C^T v_1`, the vector a STATE perturbation contracts with, which
        ## is what this method documents and what every existing gate
        ## measures.  Demir gives both conventions on one page (eq 41 with
        ## `C`, eq 42 and the phase equation 44 bare), so naming both is
        ## the fix; converting one into the other would have silently
        ## changed what `ppv()` returns.
        ## ⚠ ONE WARNING PER CALL, NOT ONE PER SAMPLE: the fill warns when
        ## `G[A,Z]` is singular, and at index 2 it is singular at every
        ## sample -- 240 identical warnings for one call, which trains a
        ## reader to filter this module's warnings and miss a real one.
        with warnings.catch_warnings(record=True) as _caught:
            warnings.simplefilter('always')
            _eq = [self._equation_row_ppv(
                       st[:m], _Xf[:, _sj if _sj < _Xf.shape[1] else -1],
                       _alg_rows, _alg_cols)
                   for _sj, st in enumerate(states)]
        _seen = set()
        for _w in _caught:
            _key = (str(_w.message), _w.category)
            if _key not in _seen:
                _seen.add(_key)
                warnings.warn(str(_w.message), _w.category, stacklevel=2)
        _v_eq = self._equation_row_ppv(v[:m], x0f, _alg_rows, _alg_cols)
        ## ⚠ A SECOND MULTIPLIER NEAR 1 BREAKS THIS SILENTLY, and none of
        ## the residuals above can see it.  The border removes the PHASE
        ## mode's singularity and does nothing about any OTHER root
        ## approaching the unit circle -- which a slow node puts there.
        ## MEASURED on van der Pol with one weakly coupled RC node:
        ##
        ##     tau/T    |lambda_2|    sigma_min(bordered)   null residual
        ##     none      0.000856          8.62e-01            4.1e-11
        ##     1e2       0.990049          4.49e-03            4.6e-11
        ##     1e4       0.999900          4.47e-05            4.6e-11
        ##     1e6       0.999999          4.47e-07            4.4e-11
        ##
        ## `sigma_min` tracks `T/tau` over six decades while the residual
        ## does not move at all: GMRES converges, the answer looks clean,
        ## and the conditioning has lost six digits.  So this estimates
        ## `|lambda_2|` explicitly rather than trusting a small residual.
        ##
        ## ⚠ AND THE ACCURACY COST WAS GATED, WITH A NEGATIVE RESULT worth
        ## recording so nobody re-derives a fix from the warning alone.
        ## Monte Carlo on the FULL NONLINEAR circuit -- 200 realisations,
        ## 150 periods, phase read from zero-crossing timing, so no PPV
        ## appears anywhere in the measurement:
        ##
        ##     core injection (control)   c_ppv/c_meas = 0.9965
        ##     slow node, tau/T = 10      c_ppv/c_meas = 0.8016
        ##
        ## Within 20%, about 2 sigma at this sample count, and in the
        ## UNDER-predicting direction.
        ##
        ## ⚠ BUT THAT IS NOT A FALSIFICATION, AND THIS DOCSTRING SAID IT
        ## WAS.  `tau/T = 10` is OUTSIDE the regime the reported mechanism
        ## needs: it bites through ill-conditioning, and by the table above
        ## `sigma_min` at `tau/T = 10` is ~4.5e-02 -- healthy.  The PPV has
        ## no large entries there and nothing is splitting into two nearly
        ## cancelling components.  Lai's own case is a gated-capacitor
        ## tuning bank (226 MOSFETs, 3.15 GHz) whose off-caps have RC
        ## exceeding ~1 s, i.e. `tau/T ~ 3e9` -- eight orders from what was
        ## tested.  A null result at 10 is what the mechanism PREDICTS, not
        ## evidence against it.
        ##
        ## ⚠⚠ PROVENANCE, AND THE CHAIN IS NOW FULLY TRACED -- A UNIT WAS
        ## MANUFACTURED IN TWO STEPS.  This used to render "larger than 1
        ## second" AS A QUOTATION.  The primary source IS on disk, at
        ## `~/docs/09-phase-macromodels-and-prc/Lai-2008-Frequency-Aware
        ## PPV ... (Cadence).pdf`, and p.4 reads, verbatim:
        ##
        ##     "Since the RC time constants of the "off" gated capacitors is
        ##      very large (LARGER THAN 1), it is safe to assume that these
        ##      gates have very small contribution to the total phase noise
        ##      when offset frequency is reasonably large."
        ##
        ## **NO UNIT.**  Our own reading of the paper
        ## (`~/docs/pycircuit-frequency-aware-ppv.md`) paraphrased it as
        ## "their RC constants exceed 1 s" -- ADDING the unit, and unmarked,
        ## beside that file's properly marked quotations.  This comment then
        ## promoted the paraphrase to a QUOTATION, carrying the added unit
        ## with it.  Two steps, each small, and the result was a quoted unit
        ## the source does not contain.
        ##
        ## ⚠ Seconds remains the natural reading (the `tau/T ~ 3e9` above
        ## follows from it and nothing downstream moves), but it is OURS and
        ## is marked as such.  ⚠⚠ AND THE FIRST VERSION OF THIS CORRECTION
        ## SAID THE PDF WAS "NOT ON DISK AT ALL" -- it is, in a
        ## SUBDIRECTORY, and the search that missed it looked only at the
        ## top level of `~/docs`.  Search a library recursively before
        ## reporting a source missing.
        ##
        ## ⚠ SO THE HONEST RECORD IS: not reproduced at `tau/T = 10`, which
        ## is outside the regime where the mechanism predicts an effect;
        ## UNTESTED at the `tau/T ~ 1e9` where it is reported.  And the
        ## reason the fix is still not built is COST, not falsification:
        ## the measurement needs ~15 time constants of settling, so at
        ## `tau/T = 1e4` that is 150 000 periods per realisation.  That
        ## argument stands on its own; the falsification framing does not,
        ## and this codebase's ledger distinguishes them.
        ##
        ## ⚠ AND THE 0.80 IS IN THE OPPOSITE DIRECTION TO THE REPORTED
        ## EFFECT.  If it survives the ~10% Monte Carlo uncertainty at 200
        ## realisations it is a separate ~20% UNDER-prediction at a `tau/T`
        ## where the conditioning is fine -- not a weak version of Lai's.
        ## At ~2.5 sigma it is not established either way, and it is
        ## recorded rather than resolved.
        ##
        ## ⚠ Larger `tau/T` is untested and the cost is why: the
        ## measurement needs ~15 time constants of settling.
        ## ⚠ AND IT TOOK THREE ATTEMPTS.  A window of 2-4 time constants
        ## read the slow mode's DECAY as diffusion; an impulse test could
        ## not resolve a 1e-11 time shift; and one noise amplitude for both
        ## circuits put a 2.5 V jump per step on an orbit of amplitude 2,
        ## because the slow node's capacitance is 6.7e-5 F against the
        ## core's 1.0.  Each time the number was read before the
        ## MEASUREMENT was shown to be in the regime it assumes.
        ##
        ## ⚠ ARNOLDI RITZ VALUES, NOT A DEFLATED POWER ITERATION -- and the
        ## replacement is BOTH more accurate and cheaper, which is rare
        ## enough to state plainly.  Power iteration converges at
        ## `|lambda_3|/|lambda_2|`, so it fails exactly where a parasitic
        ## multiplier crowds the oscillatory one.  Arnoldi does not care
        ## about that ratio.  MEASURED on van der Pol at `Q = 16` with one
        ## parasitic RC swept through it:
        ##
        ##     tau_p/T   lam2/lam3   POWER err   RITZ err
        ##       1        2.554      3.53e-14    0
        ##       4        1.206      1.52e-06    1.11e-15
        ##       8        1.065      1.41e-03    0
        ##      16        1.000      1.10e-05    2.22e-16
        ##      32        1.032      4.61e-03    2.22e-16
        ##     100        1.054      2.48e-03    2.22e-16
        ##
        ## Machine precision everywhere INCLUDING at exact degeneracy,
        ## against a power iteration losing three digits at a ratio of
        ## 1.065 -- and at `k` matvecs rather than
        ## `PPV_DEFLATION_ITERS = 30`.
        ##
        ## ⚠ THE ROUTE IS GARCIA, ROMERO & ACHA (IEEE Trans. Power
        ## Systems 37(1), 2022): Ritz values of `I - M` map back as
        ## `lambda = 1 - theta`.  They take `H` from the GMRES that
        ## already solved the Newton correction; here it is a small
        ## dedicated Arnoldi, because this call has no GMRES of its own.
        ##
        ## ⚠ EXACT AT `k = n` AND A TRUNCATION OTHERWISE.  The cap keeps a
        ## large circuit from paying `n` matvecs for a diagnostic.
        ##
        ## ⚠ A TRUNCATED `lam2` IS A LOWER BOUND, so the near-unit warning
        ## below can only UNDER-fire -- and that is Cauchy interlacing,
        ## not an accident: `theta_j >= lambda_j(A)` for a Rayleigh-Ritz
        ## projection, hence `1 - theta_2 <= lam2`.  MEASURED on a
        ## synthetic 40x40 with a verified-normal `M`
        ## (`||M^H M - M M^H||/||M||^2 = 8.9e-16`): a lower bound in
        ## 100.0% of 200 draws at `k` = 3, 5, 8 and 12.
        ##
        ## ⚠⚠ AND THIS SURVIVED A ROUND TRIP THROUGH A FALSE REFUTATION,
        ## which is why it is written out.  An intermediate version of
        ## this comment said "over-estimate, not a lower bound", on a
        ## measurement whose selection rule took the largest `|1 - theta|`
        ## after discarding `|lam - 1| < 1e-8`.  On a truncated basis that
        ## picks the UNCONVERGED UNIT-MODE Ritz value -- below 1 but above
        ## `lam2` -- so it measured its own filter and attributed the
        ## result to the phase mode contaminating `lam2`.  Selecting
        ## `theta_2` as the second-smallest Ritz value of `I - M`, which
        ## is what interlacing is about, restores the bound.
        ##
        ## ⚠ THE BOUND IS FOR A NORMAL `M`.  At forced eigenvector
        ## conditioning `cond(V) = 1e6` it fails in 100% of draws at
        ## `k = 20` -- but by exactly `1 - lam2`, i.e. the unit mode being
        ## SELECTED as `lam2`, a selection failure rather than a Ritz one:
        ## at that conditioning `|lam_1 - 1|` is 1.5e-8 to 4.3e-8 and no
        ## value-based rule separates it.  A fixed absolute tolerance is
        ## the weak point; deflating the phase mode explicitly with `q`,
        ## which this method already has, would sidestep it.
        ##
        ## ⚠⚠⚠ AND THE ESCAPE CLAUSE IS LOAD-BEARING: A CIRCUIT MONODROMY IS
        ## NOT NORMAL, AND THE BOUND FAILS ON ONE.  MEASURED 2026-09-07 on
        ## THIS FILE'S OWN `_osc_with_ladder(Q, 14, nslow)` against the dense
        ## spectrum of the SAME operator (`n` matvecs, the route the
        ## `dirk`/`full` branch below already takes), scored in the GAP
        ## because `Q ~ 1/(1 - lam2)`::
        ##
        ##     nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arnoldi
        ##      <=11   0.995706203    0.995706197      1.000       232 / 232
        ##        12   0.996324417    1.000114048     -0.031       271 / inf
        ##        13   0.996818781    0.942674586     18.020       313 / 16.9
        ##        14   0.997220139    0.999318472      0.245       359 / 1467
        ##
        ## **THE ERROR IS NOT ONE-SIGNED**, so "can only UNDER-fire" does not
        ## hold here: 13 under-estimates (19x low in `Q`), 12 and 14
        ## OVER-estimate, and 12 returns `lam2 > 1` -- a spurious UNSTABLE
        ## multiplier, which `Q` reports as `inf`.  The two failures are
        ## different: at 13 the Arnoldi never resolves `0.99682` and selects
        ## the next TRUE eigenvalue down (`0.9427`); at 14 it selects a
        ## SPURIOUS Ritz value at `0.99932` that is no eigenvalue at all.
        ##
        ## ⚠⚠ SO "NOT LIVE ON A CIRCUIT MONODROMY" (below) WAS MEASURED ON THE
        ## WRONG AXIS.  It was checked against EIGENVECTOR CONDITIONING; the
        ## trigger here is the number of distinct near-unit CLUSTERS, which
        ## `_osc_with_ladder` varies BY CONSTRUCTION and which the fixture's
        ## own test already reports reaching ~29 at `nslow = 14`.  It is not
        ## Q-specific either: the same `nslow` fails at Q = 8/16/256.
        ##
        ## ⚠⚠ IT IS A SIZING PROBLEM, AND RAISING THIS CONSTANT IS THE WRONG
        ## FIX -- MEASURED.  `k = 16` is exact on the fixture above and FAILS
        ## on a longer ladder, because the required `k` grows with `n`::
        ##
        ##     ladder/nslow   n    k=12 ratio   k=16 ratio   k=20 ratio
        ##        14 / 14     32      0.245        1.000        1.000
        ##        20 / 20     44      0.277        0.279        1.000
        ##        26 / 26     56      2.313        0.410        0.265
        ##
        ## `k ~ n/2` and rising, against a DENSE route that costs `n` and needs
        ## no threshold at all.
        ##
        ## ⚠⚠ BUT `k ~ n/2` IS AN ARTEFACT OF THIS FIXTURE, AND THE FIXTURE
        ## CANNOT SEE IT.  `_osc_with_ladder` sets `nslow = nladder`, so `n`
        ## and the slow-mode count move together here and no measurement on
        ## it can separate "k tracks n" from "k tracks nslow".  A peer
        ## session's synthetic CAN separate them and reports `k_min` rising
        ## with `nslow` and FLAT under a doubling of `n` at fixed `nslow`
        ## (24->24, 24->16, 48->48, 48->48).  If that transfers, the rule is
        ## **cost tracks the SLOW-MODE COUNT, not the system size** -- a
        ## large fast circuit is cheap and a small one with a big tuning
        ## bank is not, which also says Lai's 813-equation oscillator is
        ## expensive because of the BANK and not the 813.  Recorded with
        ## that provenance: measured on a synthetic, consistent with
        ## everything measured here, and NOT separable on this fixture.
        ##
        ## ⚠⚠ BUT THE DENSE ROUTE IS OUT ON THE CIRCUITS THAT MOTIVATE THIS.
        ## `FLOQUET_DENSE_LIMIT = 400`, and the published cases are LARGER:
        ## Lai's 64-gated-capacitor DCO is "about 200 transistors, and the
        ## system size is more than 500 ... We have trouble to apply direct
        ## harmonic balance in this case due to memory issue" (DAC 2006
        ## p.1021, verified on disk), and [L08]'s tuning oscillator is 813.
        ## So dense is the right default only in the `n <= 400` band this
        ## class already draws, and the RITZ-RESIDUAL gate is what the large
        ## end needs.  ⚠ And a 64-element bank with any realistic fraction
        ## off is an order of magnitude past the >= 3 clusters that break
        ## `k = 12` -- i.e. the extension AT ITS CURRENT BASIS SIZE would
        ## fail on exactly the circuits it exists for.
        ## THE DIAGNOSTIC THAT SEPARATES THEM CLEANLY IS THE PER-PAIR RITZ
        ## RESIDUAL `|h_{k+1,k}| |y_i[last]|`, free from `H`: 1.0e-02 at
        ## k=8, 2.1e-03 at k=12 (both wrong), 1.5e-16 at k=16 (right), and
        ## <=3.1e-07 at every `nslow` the shipped path gets right.  Neither
        ## is built -- see the roadmap; `lam2` and `Q` are REPORTED
        ## DIAGNOSTICS with no non-test consumer, so nothing computes wrong,
        ## but a caller reading `info['Q']` on a bias network with many long
        ## time constants can be off by 4x to 19x, silently.
        ##
        ## ⚠ ONE FAILURE MODE CHECKED AND NOT LIVE HERE -- ⚠⚠ SUPERSEDED BY
        ## THE MEASUREMENT ABOVE, KEPT BECAUSE IT RECORDS WHAT WAS TESTED.  At
        ## `cond(V) >= 1e4` the `|lam - 1|` filter itself fails: the phase
        ## mode stops being resolved to the tolerance, survives the
        ## discard, and is selected as `lam2`, sending `Q` to infinity.
        ## MEASURED on what was then this class's stiffest realistic fixture
        ## -- a Q=60 oscillator with a 10-mode damped bulk, `m = 12` --
        ## `cond(V) = 92` and `|lam_1 - 1| = 3.0e-13`, seven orders inside
        ## the 1e-6 filter.  That axis is still clean; the CLUSTER-COUNT axis
        ## is not.  Sorted by real part, not magnitude, because an
        ## amplitude mode is real and positive while a complex pair of
        ## larger modulus would be an oscillation about the orbit.
        vu = float(v[:m] @ u[:m] + v[m:] @ u[m:])
        lam2 = 0.0
        ## `_resid`/`_certified` describe how `lam2` was obtained; the
        ## degenerate `n < 2` fall-through never enters either branch, and
        ## `lam2 = 0` there is exact rather than estimated.
        _resid, _certified = 0.0, True
        kk = int(min(n, self.PPV_RITZ_BASIS))
        ## ⚠⚠ DENSE WHENEVER IT IS AFFORDABLE, AND THAT IS NOW THE DEFAULT
        ## RATHER THAN A STAGE-METHOD CARVE-OUT.  Forming `M` by `n` matvecs
        ## and taking its exact spectrum has no threshold, no basis size and
        ## no selection ambiguity; the truncated Arnoldi below has all three.
        ##
        ## It used to run only for `dirk`/`full`, on the argument quoted
        ## below -- and that argument was never stage-specific.  MEASURED
        ## 2026-09-07 on `_osc_with_ladder(16, 14, nslow)` (`gear`, so the
        ## Arnoldi path), against this same dense spectrum, scored in the GAP
        ## because `Q ~ 1/(1 - lam2)`::
        ##
        ##     nslow   dense lam2     k=12 Arnoldi   gap ratio   Q dense/Arn
        ##      <=11   0.995706203    0.995706197      1.000       232 / 232
        ##        12   0.996324417    1.000114048     -0.031       271 / inf
        ##        13   0.996818781    0.942674586     18.020       313 / 16.9
        ##        14   0.997220139    0.999318472      0.245       359 / 1467
        ##
        ## Not one-signed, so the Cauchy lower bound recorded above does not
        ## hold on a circuit monodromy (it is stated for a NORMAL `M`, and
        ## this is not one); and at `nslow = 12` it reports `lam2 > 1`, a
        ## spurious UNSTABLE multiplier, which `Q` turns into `inf`.
        ##
        ## ⚠ RAISING `PPV_RITZ_BASIS` IS NOT THE FIX AND WAS MEASURED NOT TO
        ## BE: `k = 16` is exact on that fixture and fails on a longer
        ## ladder (20/20 -> 0.279, 26/26 -> 0.410), because the basis has to
        ## grow with the problem.  A constant cannot.
        ##
        ## ⚠ `FLOQUET_DENSE_LIMIT` is the same cap `floquet_modes` applies to
        ## the same assembly, so the two agree about what "affordable" means.
        ## `dirk`/`full` keep the dense route ABOVE it as well: there it is
        ## expensive, but the alternative is not slower, it is WRONG, and
        ## those paths have never had the truncated one.
        _dense_ok = (fp.kind in ('dirk', 'full')
                     or n <= self.FLOQUET_DENSE_LIMIT)
        if _dense_ok:
            ## ⚠ THE STAGE MAP IS DENSE AND WIDTH `m`, so its exact spectrum
            ## is cheap -- and the Arnoldi below resolves it BADLY here.
            ## `I - M` has `M`'s annihilated modes clustered at eigenvalue 1
            ## and the physical unit root also at 1 after `1 - theta`;
            ## measured, the Arnoldi left the unit root at `1 - 1.5e-6`, past
            ## the `1e-6` deflation, so it reported the ORBIT TANGENT as the
            ## second multiplier and `Q ~ 6e5`.  Forming `M` by `m` matvecs
            ## and taking its eigenvalues directly gives the unit root to
            ## machine precision (it deflates cleanly) and the true second
            ## multiplier -- 8.59e-4 on van der Pol, matching Gear-2's
            ## 8.58e-4.
            _Md = np.column_stack([np.asarray(fp.matvec(_e), dtype=float)
                                   for _e in np.eye(n)])
            _lams = np.linalg.eigvals(_Md)
            _keep = np.real(_lams)[np.abs(_lams - 1.0) > 1e-6]
            if _keep.size:
                lam2 = float(max(np.max(_keep), 0.0))
            ## the spectrum is exact, so there is nothing to certify against
            _resid, _certified = 0.0, True
        elif kk >= 2:
            ## ⚠⚠ THE TRUNCATED PATH, NOW GATED ON THE PAIR'S OWN RITZ
            ## RESIDUAL AND GROWN UNTIL IT CERTIFIES.  This branch runs only
            ## where the dense spectrum is unaffordable -- which is exactly
            ## where the truncation is least trustworthy, since a big circuit
            ## is the one likely to carry the many slow nodes that break the
            ## selection.  A fixed basis cannot work here: `k` has to track
            ## the slow-mode count, so `PPV_RITZ_BASIS = 16` was measured to
            ## be exact on one ladder and wrong on a longer one.  Doubling
            ## until the residual certifies is the same rule at every size.
            ##
            ## ⚠ THE LOOP TERMINATES ON THREE THINGS and only one of them is
            ## a threshold: the residual certifying, the basis reaching `n`
            ## (where the Arnoldi IS the spectrum), or the cost ceiling
            ## `PPV_RITZ_MAX_BASIS` -- which produces a WARNING and an
            ## uncertified number, never a silently wrong one.
            _budget = int(min(n, self.PPV_RITZ_MAX_BASIS))
            _resid = float('inf')
            while True:
                lam2, _resid = self._ritz_second_multiplier(fp, kk)
                if _resid <= self.PPV_RITZ_RESIDUAL_TOL or kk >= _budget:
                    break
                kk = int(min(2 * kk, _budget))
            _certified = _resid <= self.PPV_RITZ_RESIDUAL_TOL
            if not _certified:
                ## ⚠ AND IT WARNS ONLY WHEN IT FAILS.  An unconditional
                ## warning on every truncated call is noise a caller learns
                ## to ignore, which is worse than none: the point of the gate
                ## is that silence now MEANS something.
                warnings.warn(
                    'PSS.ppv: `second_multiplier` (%.6f) is NOT CERTIFIED. '
                    'n = %d exceeds FLOQUET_DENSE_LIMIT = %d, so it comes '
                    'from a truncated Arnoldi, and the selected pair\'s Ritz '
                    'residual is %.2e against a tolerance of %.0e after '
                    'growing the basis to %d (ceiling %d). Measured on a '
                    'ladder oscillator, an uncertified value is wrong by '
                    '4x-19x, in BOTH directions, and can exceed 1. Treat '
                    '`second_multiplier` and `Q` as indicative; '
                    "`info['second_multiplier_certified']` says which. "
                    'Raising PPV_RITZ_BASIS is measured NOT to be the fix -- '
                    'the basis has to grow with the slow-mode count, which '
                    'is what the loop above does; raise PPV_RITZ_MAX_BASIS '
                    'if the cost is acceptable.'
                    % (lam2, n, self.FLOQUET_DENSE_LIMIT, _resid,
                       self.PPV_RITZ_RESIDUAL_TOL, kk, _budget),
                    RuntimeWarning, stacklevel=2)
        if lam2 > self.PPV_SECOND_MULTIPLIER_WARN:
            warnings.warn(
                'PSS.ppv: a SECOND Floquet multiplier sits at %.6f, near '
                'the unit circle. The bordered extraction removes only the '
                'phase mode, so its conditioning degrades as that root '
                'approaches 1 -- measured losing six digits over six '
                'decades of time constant while every residual stayed at '
                '1e-11. ⚠ AND THE PHASE EQUATION ITSELF IS THE DEEPER '
                'ISSUE: it treats the frequency response as instantaneous, '
                'so slow nodes that FILTER a device\'s noise are not seen '
                'and phase noise is OVER-ESTIMATED. Neither a smaller '
                'tolerance nor a better extraction fixes that; it needs a '
                'frequency-aware PPV. Treat this result as an upper bound.'
                % lam2, RuntimeWarning, stacklevel=2)
        ## ⚠ ONE NUMBER THAT SUBSUMES FOUR DIAGNOSTICS.  An amplitude
        ## perturbation decays to `|lambda_2|` of its size each cycle, so
        ## the cycles needed to fall below a threshold IS the oscillator's
        ## Q: `Q = log(threshold)/log|lambda_2|` (Wang & Roychowdhury).
        ## The usual definitions do not apply to an autonomous circuit --
        ## `f_r/df` presumes a Bode plot of a BIBO-stable linear system,
        ## and stored/dissipated presumes damping a self-sustaining
        ## oscillator does not have.  Nor is it the resonator's Q.
        ##
        ## ⚠ AND IT IS THE SAME CONDITION AS EVERY FAILURE THIS CLASS
        ## WARNS ABOUT.  "High Q", "a second multiplier near 1", "slow
        ## amplitude restoration" and "a long time constant" are four
        ## vocabularies for one thing -- which is why the same circuits
        ## defeat the phase row, the eigen-split, the probe's continuation
        ## and the PPV's instantaneous-response assumption.  Not four
        ## coincidences.  It costs nothing here: the Arnoldi above already
        ## produced `|lambda_2|`.
        ##
        ## ⚠⚠ AND ITS NAME IS ONLY RIGHT WHILE THE OSCILLATOR'S AMPLITUDE
        ## MODE IS THE SLOWEST NON-UNIT MODE.  A parasitic with
        ## `tau_p/T > Q_osc` simply IS the second multiplier -- by
        ## definition, not by error -- and then this reports THE
        ## PARASITIC'S DECAY TIME IN PERIODS under the name `Q`.  MEASURED:
        ## at `tau_p/T` = 32 and 100 on a `Q = 16` oscillator, `lam2` is
        ## the parasitic and `Q` returns 32 and 100.
        ##
        ## ⚠ THE NUMBER IS RIGHT AND ITS NAME IS WRONG, which is why
        ## nothing misbehaves: every residual stays clean and the value is
        ## well converged.  A DCO's gated capacitor sits at
        ## `tau_p/T ~ 1e4`, i.e. permanently in that regime, so on exactly
        ## the circuits a hierarchical DCO method exists for, a reported
        ## `Q` would be the gated cap's RC in periods.  Read `Q` as
        ## "cycles for the SLOWEST NON-UNIT MODE to decay by 1/e", which is
        ## what it computes; it is the oscillator's Q only when that mode
        ## is the oscillator's.
        ##
        ## Reported for a `1/e` threshold, so `Q` is cycles-to-1/e.
        ##
        ## ⚠⚠ AND THIS LINE IS WANG & ROYCHOWDHURY'S IDENTITY
        ## `Q = log(threshold)/log|lambda_2|`, which does double duty and
        ## was shipped before either use was noticed.  It is what makes
        ## "bounded by Q" and "bounded by lambda_2" the SAME SENTENCE --
        ## the organising fact of this whole area, since a designer's
        ## objective (raise Q) IS the numerics' failure mode (lambda_2 ->
        ## 1).  It is also an ERROR AMPLIFIER:
        ##
        ##     (dQ/Q) / (dlambda_2/lambda_2)  =  -1/ln(lambda_2)  =  Q
        ##
        ## ⚠ SO THE RELATIVE ERROR IN `Q` IS `Q` TIMES THE RELATIVE ERROR
        ## IN `lambda_2`, and a caller reading `Q` at high Q is reading a
        ## quantity far less accurate than the multiplier behind it.
        ## MEASURED end to end on van der Pol tuned by `mu = 1/(2 pi Q)`,
        ## against the finest grid:
        ##
        ##     Q      npts   rel err lam2   rel err Q   ratio
        ##      3.18   120    1.04e-03      3.31e-03      3.2
        ##     15.92   120    5.75e-04      9.23e-03     16.1
        ##     63.66   120    4.89e-04      3.21e-02     65.7
        ##     63.66   480    7.56e-06      4.81e-04     63.7
        ##
        ## ⚠ THAT IS A RESOLUTION REQUIREMENT SCALING WITH `Q`, NOT A
        ## FIXED ACCURACY: 120 points/period gives `Q` to 0.3% at Q = 3
        ## and only 3.2% at Q = 64.  Payable here because Gear-2's
        ## `lambda_2` converges at better than second order (~8x per
        ## doubling); a method that BIASES `lambda_2` at fixed order has
        ## no such escape, and backward Euler's 5.6e-2 bias would become
        ## 85% in `Q` at Q = 100.
        Q = (-1.0 / np.log(lam2) if 0.0 < lam2 < 1.0 else float('inf'))
        info = {'border_residual': y,
                'tangent_border_residual': yf,
                'Q': Q,
                'null_residual': resid / max(float(np.linalg.norm(v)), 1e-300),
                ## ⚠ MULTIPLY `null_residual` BY THIS TO GET THE RELATIVE
                ## ERROR IN `v` THE RESIDUAL CANNOT EXCLUDE.  `null_residual`
                ## is `||v - M^T v|| / ||v||`, so an error component along the
                ## `lam2` left-eigendirection enters it scaled by `1 - lam2`
                ## and is nearly INVISIBLE exactly when `lam2 -> 1`.
                ## MEASURED on `_vdp_with_slow_node`, injecting a 1% error
                ## into a converged `v` (floor 4.6e-11):
                ##
                ##     lam2        r(random dir)   r(lam2 dir)   0.01*(1-lam2)
                ##     0.000856      1.65e-02       1.003e-02      9.99e-03
                ##     0.990049      1.65e-02       9.950e-05      9.95e-05
                ##     0.999900      1.65e-02       1.000e-06      1.00e-06
                ##     0.999999      1.65e-02       1.000e-08      1.00e-08
                ##
                ## Exact to every digit printed.  A RANDOM error is caught
                ## nine orders above the floor, so `null_residual` is a real
                ## gate and this module's assertions on it can fail -- but it
                ## loses sensitivity in the ONE direction that matters as the
                ## circuit gets better, which is the opposite of the
                ## reassurance a flat residual gives.
                ##
                ## ⚠⚠ THIS IS WHY A FLAT `null_residual` IS NOT EVIDENCE OF
                ## ACCURACY.  A residual that does not move while `lam2`
                ## sweeps toward 1 is not reporting that the answer stayed
                ## good; the bordered system is well-conditioned BY
                ## CONSTRUCTION, and the quantity it fails to see is
                ## precisely the one that grows.  Read the two numbers
                ## together or neither.
                'null_residual_amplification': (
                    1.0 / max(1.0 - lam2, np.finfo(float).eps)),
                'second_multiplier': lam2,
                ## Which route produced it, so a caller can tell an exact
                ## spectrum from a truncated estimate without re-deriving
                ## the rule.  See the branch above.
                'second_multiplier_route': ('dense' if _dense_ok
                                            else 'arnoldi'),
                ## The selected pair's Ritz residual and whether it cleared
                ## `PPV_RITZ_RESIDUAL_TOL`.  Exact on the dense route (0.0,
                ## True).  A caller that reads `Q` should read this too.
                'second_multiplier_residual': float(_resid),
                'second_multiplier_certified': bool(_certified),
                'q': q, 'xdot': xdot, 'tangent_pair': u,
                'samples': np.asarray(states),
                ## ⚠ `samples_eq` IS THE ONE TO CONTRACT `CY` AGAINST.
                ## `samples` is `C^T v_1` (a state perturbation's
                ## sensitivity); this is `v_1` (an equation-row input's).
                ## A noise current injected into a KCL row is the latter.
                'samples_pair': np.asarray(states_pair),
                'monodromy_method': getattr(self.par, 'method', '?'),
                'samples_eq': np.asarray(_eq),
                'v_eq': _v_eq,
                'times': np.asarray(fp.times, dtype=float)}
        return v, info

    def frequency_aware_ppv(self, offset, tol=None):
        """The PPV at a nonzero modulation frequency (Lai 2008, eq. 23).

        ⚠ EQ. 23, NOT 24 (docs session, 2026-09-09): Lai's eq. 24 drops the
        AC columns of the Toeplitz block and is justified only "if we are
        only interested in the transfer functions when w_s is close to
        DC"; this shooting form has no such truncation -- `I - exp(-j w_s T)
        M^T` is the exact sampled LPTV adjoint at ANY offset, the monodromy
        already carrying the full time variation -- so it is eq. 23 for
        what the object IS, exact for the discretised system, and eq. 24
        only for the DC-reduction sentence pinned below.

        The classical PPV is the left null vector of `I - M^T`, bordered by
        `q = C(0) xdot(0)`; it is the phase response to a perturbation that
        is SLOW against every other Floquet mode.  This is the SAME
        bordered system at `alpha = exp(-j w_s T)`,

            [[I - alpha M^T,  q], [q^T, 0]] [v; y] = [0; 1],

        whose solution `v(w_s)` is the phase sensitivity to a perturbation
        modulated at `w_s`: at `w_s = 0` it IS `ppv()` (pinned), and away
        from it the AMPLITUDE mode admixes with weight
        `(1 - alpha)/(1 - alpha mu_2)` -- zero at DC, rising ten-fold per
        decade, cornering where `2 pi f_s T = 1 - mu_2` and flat above
        (docs session, 2026-09-08, on the slow-node fixture; the corner
        tracks the slow multiplier over two decades).  Verified at the
        source: Lai's eq. (24) at `w_s = 0` "is the augmented PPV
        extraction equation", verbatim; his construction is harmonic
        balance, this is the shooting basis, and the two agree on what the
        object is.

        Returns `(v, info)`: `v` the pair-space anchor vector (complex);
        `info['samples_pair']` the T-periodic envelope over the period,
        `lambda_k exp(+j w_s t_k)` with `lambda` the transposed replay of
        `v` -- the sideband rows' own convention, so its Fourier
        coefficient at harmonic `k` is the phase transfer of a source band
        at `k f0 + f_s`; `info['samples']` the same in state space
        (`C^T v`), SECOND order in the step through `ppv()`'s own
        pair-consistent propagation run on the complex anchor
        (`_ppv_propagate`, lifted 2026-09-08; at w_s = 0 it is `ppv()`'s
        `samples`); `info['admixture']` the norm fraction of `v(w_s)`
        orthogonal to the DC PPV; `info['corner']` the predicted corner
        `|1 - mu_2| / (2 pi T)` in Hz; `info['alpha']`; `info['ppv']` the
        DC object's info.

        ⚠ WHAT IT IS FOR.  A source that reaches the phase through a slow
        path (an RC leg, tau >> T) is filtered at its own corner, and the
        DC PPV cannot see that; the harmonic sum built from THIS object's
        coefficients carries the filter inside `c_k(w_s)` with no explicit
        model of the path (A2, roadmap).  MEASURED (2026-09-08): the ratio
        `sum_k |c_k(w_s)|^2 / sum_k |c_k(0)|^2` reproduces `pnoise`'s
        `S_pm/(4 S_v)` for a source behind the slow node within 0.6 % to
        r = 1e-2 and 2 % at 5e-2 (a = 0.4, loss 0.2, tau/T = 100), where the
        DC sum with the filter by hand was 4 % off; at r = 0.1 the two
        differ by -6 %, unchanged at twice the grid -- a gap between PM by
        sideband quadrature and phase-mode projection, both 1e-3 of DC
        there -- and put to a nonlinear Monte Carlo that instantiates
        neither construction (2026-09-09; 4 seeds x 10 000 periods per
        point, phase read two ways).  ⚠ A first reading at one asymmetry
        (a = 0.25) assigned each construction to ONE phase definition
        crosswise (crossing phase to S_pm at 1.027, demodulated phase to
        this sum at 0.977, each +-1.8 %); the asymmetry sweep a = 0.25 /
        0.12 / 0.05 / 0 REPLACED it.  Double ratios (Monte Carlo slow/core
        over predicted slow/core, each estimator calibrated on the core):

            a      crossing vs S_pm / this sum   demod vs S_pm / this sum
            0.25   1.032 / 1.066  (+-1.8 %)      0.954 / 0.985
            0.12   1.093 / 1.121  (+-3.5 %)      0.974 / 0.999
            0.05   1.122 / 1.149  (+-4 %)        1.004 / 1.027
            0.00   1.134 / 1.162  (+-4 %)        1.027 / 1.052

        ⚠⚠ AT 16 SEEDS PER POINT (Andreas, same day; +-1.5-1.7 %):

            a      crossing vs S_pm / this sum   demod vs S_pm / this sum
            0.25   1.038 / 1.072                 0.957 / 0.989
            0.12   1.085 / 1.113                 0.979 / 1.004
            0.05   1.117 / 1.143                 1.011 / 1.035
            0.00   1.129 / 1.156                 1.036 / 1.061

        The statistic is the SLOPE in a, not any one point (the docs
        session's framing): demod vs this sum -0.27 +- 0.09 per unit a
        (3.2 sigma), vs S_pm 3.6 sigma; crossing 4.2 / 4.7 sigma.  So the
        ratio is NOT constant in a at ~3 sigma for the demodulated phase
        and above 4 for the crossing: NEITHER construction describes
        EITHER measured phase across the range.  The two constructions
        track each other to 1 % over the sweep while both estimators --
        a point sample at a crossing and an average over a period -- drift
        together, in the same direction, against both.  The a = 0.25
        agreement of the demodulated phase with this sum is where its
        curve crosses the sum, not a match.  What the crossing carries
        beyond that: two thirds of its excess is waveform content beyond
        0.5 f0 from the carrier (an instantaneous crossing aliases the
        additive noise a one-period demodulation cannot see), the rest
        the demod's own boxcar loss (1/sinc^2 = 1.045 at the band centre)
        plus a common amplitude-to-crossing gain of ~0.6 (band-limited
        coherence).  The common drift of BOTH estimators against BOTH
        linear constructions as a -> 0 is the open object; candidate, a
        second-order amplitude-to-phase conversion the linear theory
        cannot contain -- REFUTED the same evening (PSD/4: the slow
        fixture scales linearly, 0.993 +- 0.013; grid doubling moves the
        constructions < 0.7 %; band conventions identical).  RESOLVED by a
        forward tone-transient route on the MC's own discretised system
        (no adjoint, no sideband assembly): pnoise's S_pm agrees with the
        forward LPTV PM sidebands to 1 % at both asymmetries, and the
        estimators' a-dependent double ratios are REPRODUCED by that
        deterministic linear route (demod 0.938 -> 1.003 against the MC's
        0.957 -> 1.036).  The drift is the ESTIMATORS: a one-period
        fundamental demodulation leaks the other harmonics' sidebands
        through its boxcar (sinc(pi(1-r)) ~ 0.1 for the second harmonic's,
        which is ~a), zero crossings convert every harmonic's; both read a
        given sideband PM with a fixture-dependent gain.  This object and
        S_pm are PM by quadrature of the FUNDAMENTAL'S sidebands; compare
        them with that, not with a demodulated or crossing phase.  ⚠ The
        premise "a -> 0 makes the definitions coincide" was wrong: the
        asymmetry removes even harmonics only, van der Pol's third stays
        at 9.7 % of the fundamental, and the construction gap GROWS as
        a -> 0 (core 1.010 -> 1.055); the sinusoidal limit is mu -> 0.  The
        instrument hypothesis (spectrum analyser <-> this sum, time-
        interval analyser <-> S_pm) is refuted in its crossing half.  The slow multiplier's own coefficient
        (`mode_content[0]`) corners at 1.6e-3 f0 for tau/T = 100 with a
        plateau of 2.45e-6 (the docs session's 2.29e-6), scaling as
        T/tau.  ⚠ DO NOT GATE ON `|v|`: with
        `q^T v = 1` the `1/(1 - alpha)` pole cancels between numerator
        and denominator and the norm is frequency-flat by construction; a
        one-percent orthogonal admixture moves it by 5e-5.  The change is a
        DIRECTION -- read `admixture`, or the per-harmonic coefficients.
        """
        import scipy.sparse.linalg as spla
        v0, info0 = self.ppv(tol)
        fp = self.factored_period()
        m = self.cir.n - 1
        n = fp.width
        irn = self.irefnode
        T = float(fp.T)
        alpha = np.exp(-2j * np.pi * float(offset) * T)
        q = np.asarray(info0['q'], dtype=float)
        qp = np.concatenate((q, np.zeros(n - m))).astype(complex)

        def _mv(z):
            z = np.asarray(z, dtype=complex)
            v_, y_ = z[:n], z[n]
            top = v_ - alpha * np.asarray(fp.matvec_transposed(v_), dtype=complex) + y_ * qp
            return np.concatenate((top, [complex(qp @ v_)]))
        rtol = max(self.par.reltol * 1e-2 if tol is None else tol, 1e-14)
        rhs = np.zeros(n + 1, dtype=complex)
        rhs[n] = 1.0
        z, relres, _H, _k = _arnoldi_gmres(_mv, rhs, rtol=rtol, maxiter=min(n + 1, 200))
        if relres > max(1e3 * rtol, 1e-8):
            ## Lai's own warning about this object: eq. 23 "is very
            ## difficult to solve using iterative solvers (such as GMRES),
            ## because the extra columns and rows from the Toeplitz block
            ## degrade the block diagonal preconditioner" -- a property of
            ## the formulation, not a defect of the fixture
            raise RuntimeError(
                'PSS.frequency_aware_ppv: the bordered solve at offset %g did '
                'not converge (relative residual %.3e) -- a known property of '
                'the bordered LPTV adjoint (Lai 2008), not a fixture defect.'
                % (float(offset), relres))
        v = np.asarray(z[:n], dtype=complex)
        ## `ppv()`'s own normalisation, `v . xdot(0) = 1` -- the border fixes
        ## `q^T v = 1` only, and the two differ by a factor AND a sign on
        ## the van der Pol (measured 2.07 there; -0.94 on the slow-node
        ## fixture): at alpha = 1 this is what makes the object `ppv()`.
        xdot = np.asarray(info0['xdot'], dtype=float)
        vx = complex(v[:m] @ xdot)
        if vx == 0.0:
            raise ValueError('PSS.frequency_aware_ppv: v(w_s) is orthogonal to the orbit tangent.')
        v = v / vx
        ## the envelope over the period, in the sideband rows' convention:
        ## `ppv()`'s own second-order propagation on the complex anchor
        ## (lifted into `_ppv_propagate`), then the per-step phase
        irn = self.irefnode
        x0r = np.asarray(self._period_state[1], dtype=float).ravel()
        x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
        _alg_rows, _alg_cols = self._algebraic_adjoint_pattern(x0f)
        states, states_pair, _ts, _Xf = self._ppv_propagate(fp, v, m, xdot, _alg_rows, _alg_cols)
        st = np.asarray(states_pair, dtype=complex)
        tms = np.asarray(fp.times, dtype=float)[:st.shape[0]]
        phase = np.exp(2j * np.pi * float(offset) * tms)
        samples_pair = st * phase[:, None]
        samples = np.asarray(states, dtype=complex) * phase[:, None]
        ## the admixture: what of v(w_s) is NOT along the DC PPV.  ⚠ This is
        ## dominated by whichever mode has the largest weight, which on a
        ## core with a fast amplitude mode is THAT one (weight ~ 2 pi r,
        ## cornering at r ~ (1 - lambda_2)/2pi ~ 0.16 for the van der Pol),
        ## so a slow node's own admixture -- 1e-6 to 1e-2 of it -- is
        ## invisible in the norm.  `mode_content` reads each mode's own
        ## coefficient (dense route below FLOQUET_DENSE_LIMIT).
        v0c = np.asarray(v0, dtype=complex)
        v0n = v0c / np.linalg.norm(v0c)
        proj = np.vdot(v0n, v) * v0n
        admixture = float(np.linalg.norm(v - proj) / np.linalg.norm(v))
        mode_content = None
        multipliers = None
        if n <= self.FLOQUET_DENSE_LIMIT:
            Md = np.column_stack([np.asarray(fp.matvec(e), dtype=float) for e in np.eye(n)])
            mu, P = np.linalg.eig(Md.T)
            order = np.argsort(-np.abs(mu))
            mu, P = mu[order], P[:, order]
            a = np.linalg.solve(P, v)
            mode_content = np.abs(a[1:]) / abs(a[0])
            multipliers = mu
        mu2 = float(info0.get('second_multiplier', 0.0))
        info = {'alpha': alpha, 'offset': float(offset), 'admixture': admixture,
                'mode_content': mode_content, 'multipliers': multipliers,
                'corner': abs(1.0 - mu2) / (2.0 * np.pi * T), 'second_multiplier': mu2,
                'samples_pair': samples_pair, 'samples': samples, 'times': tms,
                'residual': float(relres), 'ppv': info0}
        return v, info

    def _forced_replay_transposed(self, fp, freq, xa):
        """`W^T xa` -- the transpose of the map `u -> w(freq)`.

        THE MANY-TO-ONE HALF, and the reason adjoint noise is affordable.
        The forward replay answers "what does THIS source do at the
        output"; one run per source.  This answers "what does the output
        owe to EVERY source", in one reverse pass -- which is the shape
        pnoise has, with hundreds of sources and one output.  Okumura et
        al. (1993) choose the adjoint for exactly this: "it is efficient to
        use the adjoint method ... because circuits have many noise
        sources."

        The identity is one line.  The forward replay makes each step
        `Px_j = -Jf_j^-1 (S_j + u e^{jw t_j})`, so the final state's
        sensitivity to `u` through step `j` is `-Jf_j^-T` applied to the
        adjoint state there, weighted by `e^{jw t_j}`.  That solve is
        already taken by the reverse pass -- it is `t` in
        `_monodromy_matvec_transposed` -- so this is that pass plus a
        weighted sum, with no second recursion to keep in step.

        ⚠ NO LONGER SOLVED-HISTORY ONLY (B8).  This said "because the
        reverse pass is", which was true when written and is not now: the
        plain path has its own reverse recursion, so the dispatch belongs
        to `FactoredPeriod` and this reads `ts` from whichever map it was
        handed.  `ts[j]` is the transposed solve at step `j` under BOTH
        recursions -- that is the quantity the identity below needs, and
        it is what makes this method map-agnostic rather than merely
        permitted.
        """
        if fp.kind == 'dirk':
            return self._forced_replay_transposed_dirk(fp, freq, xa)
        if fp.kind == 'full':
            return self._forced_replay_transposed_full(fp, freq, xa)
        _end, ts, _states = fp.matvec_transposed(xa, collect=True)
        jw = 2j * np.pi * float(freq)
        acc = np.zeros(self.cir.n - 1, dtype=complex)
        for tvec, t in zip(ts, fp.times[1:]):
            acc = acc - np.exp(jw * float(t)) * np.asarray(tvec)
        return acc

    def _forced_replay_transposed_full(self, fp, freq, xa):
        """`W^T xa` for Radau IIA(3) -- the COUPLED three-stage adjoint (no
        output injection), the source-coupling sibling of
        `_monodromy_matvec_transposed_full`.

        A source at frequency `freq` enters stage `k` (abscissa
        ``t_{j,k} = t_j + c_k h``) through ``K_k = -(i + u)``, so the forward
        source RHS of the coupled step is ``S_i = -h sum_k A_ik du_k`` with
        ``du_k = exp(jw t_{j,k})``.  The endpoint reads block 3, so the
        adjoint of one step's source-to-endpoint sensitivity is, with
        ``p = J^{-T} [0; 0; w]`` (the same coupled transposed solve the
        monodromy transpose takes),

            acc += -h sum_k exp(jw t_{j,k}) (A_1k p_1 + A_2k p_2 + A_3k p_3)

        and the costate propagates back by ``w <- Cn^T (p_1 + p_2 + p_3)``.
        ⚠ NO TWO-VECTOR SHORTCUT: the source couples through all three stages
        with the full ``A (x) B`` weighting, so all three `m`-blocks of `p`
        are read (a DIRK's ``A1 p`` feed is the two-stage special case of
        this).  Dual-consistent with the forward replay and matched to
        forward driven solves; the injected sibling is
        `_sideband_forced_full`.
        """
        import scipy.linalg as sla
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        Amat, _Bw, cvec = integ.butcher()
        s = Amat.shape[0]
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        w = np.asarray(xa, dtype=complex).ravel().copy()
        acc = np.zeros(m, dtype=complex)
        tms = np.asarray(fp.times, dtype=float)

        def csolveT(lu, b):
            return (sla.lu_solve(lu, b.real, trans=1)
                    + 1j * sla.lu_solve(lu, b.imag, trans=1))

        for j in range(len(fp.steps) - 1, -1, -1):
            lu, Cn, mm = fp.steps[j]
            ts = tms[j]; te = tms[j + 1]; h = te - ts
            b3 = np.concatenate([np.zeros(m)] * (s - 1) + [w])
            p = csolveT(lu, b3)
            pb = [p[i * m:(i + 1) * m] for i in range(s)]
            for k in range(s):
                tk = ts + cvec[k] * h
                coup = sum(Amat[i, k] * pb[i] for i in range(s))
                acc = acc - h * np.exp(jw * tk) * coup
            w = Cn.T @ sum(pb)
        return acc

    def _sideband_forced_full(self, fp, freq, l, d):
        """The forced (source-injected) part of the Radau IIA(3) sideband row
        for sideband `l`, and the final costate `g` for the closure.

        The three-vector analogue of `_sideband_forced_trbdf2`: the reverse
        pass injects the OUTPUT functional `d` (weighted by
        `exp(-j(l w0 + w) t_n)/N`) at each step and reads the SOURCE coupling
        through ALL THREE stages at every step -- the coupled fold
        (`_forced_replay_transposed_full` with the injection added).  ⚠ The
        injection is added AFTER the step's costate update, so the output at
        step `n` couples to the source at steps `< n` (causality); the
        accepted trajectory state at `t_n` IS the entering state `x_n` of step
        `n`, so `d` couples at `t_n` (`ts`) exactly as in the TR-BDF2 fold.
        The source coupling folds the three abscissae with
        `sum_k exp(jw t_{n,k}) sum_i A_ik p_i` -- there is no two-vector
        shortcut (`A (x) B`).  Returns `(forced, g)`.
        """
        import scipy.linalg as sla
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        Amat, _Bw, cvec = integ.butcher()
        s = Amat.shape[0]
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        T = float(fp.T)
        w0 = 2.0 * np.pi / T
        N = len(fp.steps)
        tms = np.asarray(fp.times, dtype=float)
        d = np.asarray(d, dtype=complex).ravel()
        lam = np.zeros(m, dtype=complex)
        forced = np.zeros(m, dtype=complex)

        def csolveT(lu, b):
            return (sla.lu_solve(lu, b.real, trans=1)
                    + 1j * sla.lu_solve(lu, b.imag, trans=1))

        for j in range(N - 1, -1, -1):
            lu, Cn, mm = fp.steps[j]
            ts = tms[j]; te = tms[j + 1]; h = te - ts
            p = csolveT(lu, np.concatenate([np.zeros(m)] * (s - 1) + [lam]))
            pb = [p[i * m:(i + 1) * m] for i in range(s)]
            for k in range(s):
                tk = ts + cvec[k] * h
                coup = sum(Amat[i, k] * pb[i] for i in range(s))
                forced = forced - h * np.exp(jw * tk) * coup
            lam = Cn.T @ sum(pb)
            lam = lam + np.exp(-1j * (float(l) * w0
                                      + 2.0 * np.pi * float(freq)) * ts) / N * d
        return forced, lam

    def monodromy_twin(self):
        """The `PSS` whose monodromy the oscillator surfaces read.

        `self` for a driven circuit, when `self.monodromy == 'native'`, or
        when this run's own method already carries a second-order monodromy
        (`gear`/`trbdf2`).  Otherwise (an autonomous circuit solved with a
        one-step LMM, trap or euler) a twin `PSS` of the same circuit under
        the method `self.monodromy` names, solved once on the SAME grid from
        this orbit's converged state and cached.

        ⚠⚠ THE TWIN DEFAULTS TO TR-BDF2, MEASURED (2026-09-05), and Gear-2
        is one setting away, not retired.  The twin exists because trap's
        and euler's own monodromy is unusable on a limit cycle: trap's
        diverges with refinement (B16, below) and euler's is first order.
        Both TR-BDF2 and Gear-2 give a clean second-order twin, but TR-BDF2
        is more accurate on `lambda2` -- measured against exact references
        (`exp(A T)` on a linear oscillator; Abel's `exp(mu integral(1-v^2))`
        on van der Pol) it is 12-32x better at practical step counts, and
        the advantage GROWS with Q and with coarser grids -- the regime a
        real oscillator PSS sits in.  At Q=100 and 50 points/period the
        Gear-2 twin misreads Q by 22%, the TR-BDF2 twin by 0.27%.  The gap
        is a coarse-grid/high-Q effect, not a fixed factor: refine the grid
        or drop Q and both fall to the ordinary O(h^2) floor where the
        difference is single digits and can even favour Gear-2.  So
        `monodromy = 'trbdf2'` is the default, `'gear'` restores the former
        twin, `'native'` reads the run's own.

        ⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
        bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
        5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
        with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
        at 400/800/1600 points and DIVERGES with refinement, and
        `x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
        multiplier at 1 (the one-step companion's parasitic mode), while
        the state and period are second order either way.  So "the most
        accurate" is not a choice between openers: the STATE keeps the
        method you asked for, and every monodromy-derived quantity -- `Q`,
        the PPV and everything built on it, the Floquet modes, the
        phase-noise surfaces -- comes from the twin on the same orbit.  The
        twin's period differs from this one's by O(h^2); its orbit is
        re-converged, not copied.

        ⚠⚠ THE B16 DECISION, TAKEN ON A MEASUREMENT (2026-09-05).  On the
        bias-sensitive oscillator (`vdp + 0.3 u^2`, exact `Q_lambda =
        5.908`, `c_true = 5.3703e-06`) trapezoidal's monodromy is unusable
        with EITHER opener: the default reads `Q_lambda` 11.1 / 28.4 / 63.9
        at 400/800/1600 points and DIVERGES with refinement, and
        `x0_unknown=True` reads 3086 / 12228 / 48699 -- a spurious
        multiplier at 1 (the one-step companion's parasitic mode), while
        the state and period are second order either way.  Gear-2 reads
        5.9094 / 5.9086 / 5.9084.  So "the most accurate" is not a choice
        between openers: the STATE keeps the method you asked for, and
        every monodromy-derived quantity -- `Q`, the PPV and everything
        built on it, the Floquet modes, the phase-noise surfaces -- comes
        from Gear-2 on the same orbit.  The twin's period differs from
        this one's by O(h^2); its orbit is re-converged, not copied.
        """
        mono = getattr(self, 'monodromy', 'trbdf2')
        if mono not in ('trbdf2', 'gear', 'native'):
            raise ValueError(
                "PSS.monodromy must be 'trbdf2', 'gear' or 'native', not %r"
                % (mono,))
        _integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if (mono == 'native'
                or _integ.carries_own_monodromy()
                or not getattr(self, 'autonomous', False)):
            ## ⚠ SELF-SUFFICIENT METHODS TAKE NO TWIN, and the method says which
            ## it is (`carries_own_monodromy`): Gear-2 (a second-order native
            ## companion monodromy) and every stage method (TR-BDF2, Radau -- no
            ## opener seam, verified against the pencil).  This twin exists only
            ## because a one-step LMM's monodromy is first-order on a limit
            ## cycle (its manufactured opener is dropped to Euler and that seam
            ## sits in the period map); twinning a self-sufficient method would
            ## replace its own map with another on a re-converged orbit -- pure
            ## cost -- and hide the run's own spectrum.  They read `native`
            ## regardless of `monodromy`; the knob governs which twin trap/euler
            ## borrow.
            return self
        if getattr(self, '_period_state', None) is None or not self.converged:
            return self
        twin = self._solve_twin(mono)
        self._monodromy_twin = twin
        return twin

    def _solve_twin(self, method):
        """A converged twin `PSS` of this circuit under `method`, re-solved
        once on the SAME grid from this orbit's converged state, cached per
        method.  The shared machinery behind `monodromy_twin` (which picks
        the method by `self.monodromy`) and `_lyapunov_host` (which forces
        `gear`, because the noise-injection surfaces cannot use a TR-BDF2
        twin yet).  Raises if the re-solve does not converge.
        """
        cache = self._twins
        if method in cache:
            return cache[method]
        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        kw = dict(self._solve_kwargs)
        twin = PSS(self.cir, toolkit=self.toolkit, irefnode=None,
                   method=method, reltol=self.par.reltol,
                   iabstol=self.par.iabstol, vabstol=self.par.vabstol)
        hs = np.asarray(hs, dtype=float)
        ## the same grid: its fractions when it is not uniform, else the
        ## uniform step (a one-step plain state can carry an `hs` whose
        ## sum is a trial period, so the fractions are the safe object)
        nonuniform = float(hs.max() / hs.min()) > 1.0 + 1e-9
        grid = (hs / float(hs.sum())) if nonuniform else None
        x0r = np.asarray(x0, dtype=float)[:self.cir.n - 1]
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                twin.solve(refnode=kw.get('refnode', gnd), period=float(T),
                           x0=x0r, timestep=float(T) / len(hs), grid=grid,
                           maxiterations=max(int(kw.get('maxiterations',
                                                        20)), 20),
                           matrix_free=bool(kw.get('matrix_free', False)))
            _ok = bool(twin.converged)
        except NoConvergenceError:
            _ok = False
        if not _ok:
            raise RuntimeError(
                'PSS.monodromy_twin: the %s re-solve from the converged '
                '%s orbit did not converge, so no second-order monodromy is '
                'available; try pss.monodromy = "gear" (the other twin) or '
                '"native" to read the one-step method\'s own (first-order, '
                'and on an oscillator its second multiplier is not the '
                'physical one).'
                % (method, getattr(self.par, 'method', '?')))

        ## ⚠⚠ THE TWIN MUST HAVE CONVERGED TO THE SAME ORBIT IT WAS SEEDED
        ## ON, and a MORE ROBUST twin makes this check load-bearing rather
        ## than paranoid.  Measured: from a poor seed (euler at 400 pts, its
        ## orbit 55% off) the Gear-2 twin fails to converge -- LOUD -- but
        ## the TR-BDF2 twin, being more robust, CONVERGES to a SPURIOUS limit
        ## cycle and reports `Q = 1.97` against the exact 5.91 with no error.
        ## Improving the method degraded safety: the failure moved from a
        ## refusal to a plausible wrong number.  So the twin's converged
        ## orbit is checked against the seed it was handed: two convergent
        ## methods on the SAME limit cycle agree on period and entering state
        ## to O(h^p) -- measured 6e-6 / 1e-4 (trbdf2) and 3e-5 / 2e-2 (gear)
        ## on a good seed -- while the spurious jump above sits at 2.17 /
        ## 0.91.  The 0.25 gate is ~12x above the worst good case and ~3.6x
        ## below the spurious one; it is set from the (universal, tiny)
        ## good-case agreement, not the (fixture-dependent) failure size, so
        ## it transfers.  The definitive test is refinement (a spurious orbit
        ## does not survive h/2); this cheap consistency check is the
        ## conservative stand-in -- it REFUSES a too-poor seed rather than
        ## risk trusting it, which is the safe direction.
        Th = float(T)
        dT = abs(float(twin.period) - Th) / max(abs(Th), 1e-30)
        x0t = np.asarray(twin._period_state[1],
                         dtype=float)[:self.cir.n - 1]
        dx = (float(np.linalg.norm(x0t - x0r))
              / (float(np.linalg.norm(x0r)) + 1e-30))
        if dT > 0.25 or dx > 0.25:
            raise RuntimeError(
                'PSS.monodromy_twin: the %s twin converged to a DIFFERENT '
                'orbit than the seeding %s orbit (period Delta = %.2e, state '
                'Delta = %.2e, gate 0.25): that orbit is too poor to seed a '
                'monodromy twin -- the free-period Newton reached a spurious '
                'limit cycle.  Refine the grid so the %s state is a good '
                'seed, or set pss.monodromy = "native" to read the run\'s own '
                '(defective) monodromy rather than a wrong number.'
                % (method, getattr(self.par, 'method', '?'), dT, dx,
                   getattr(self.par, 'method', '?')))
        twin.monodromy = 'native'
        cache[method] = twin
        return twin

    def _lyapunov_host(self):
        """The `PSS` the Lyapunov noise surfaces (`covariance`,
        `oscillator_covariance`) read.

        The covariance is propagated on the SAME orbit the Floquet quantities
        use, so this is just the monodromy twin: a trap/euler autonomous run
        hands its covariance to the twin (`TR-BDF2` by default, `gear` if
        `monodromy='gear'`) so both come from one consistent orbit; a
        gear/trbdf2 host is its own host.  TR-BDF2's per-step injection is
        built (`_lyapunov_pieces_trbdf2`, DAE-projected Van Loan), so there
        is no Gear-2 fallback -- the injection follows the chosen twin.
        """
        return self.monodromy_twin()

    def _adjoint_host(self):
        """Host for the ADJOINT SIDEBAND noise surface (`pnoise`, via
        `adjoint_sideband_row`).

        The two-stage sideband fold IS built for TR-BDF2
        (`_sideband_forced_trbdf2`, the two-vector injected reverse pass that
        carries the source coupling through both stages), so this is the
        monodromy twin -- the same orbit the Floquet and Lyapunov surfaces
        use, no Gear-2 fallback.  `monodromy='gear'` still routes to the
        Gear-2 twin if asked.
        """
        return self.monodromy_twin()

    ## Nominal convergence order per `method`, for `grid_error`'s ceiling on
    ## a plausible OBSERVED order.  Sourced from this file's own measured
    ## records rather than from the literature: trap/gear/theta second order,
    ## TR-BDF2 measured at 4.01x/4.01x/4.00x per halving (exact `O(h^2)`),
    ## Radau IIA(3) at 31.50x/31.74x (`O(h^5)`, theoretical 32), euler first.
    ## ⚠ An unlisted method falls back to a generic range and the ceiling is
    ## not applied -- add it here rather than letting it default silently.
    METHOD_ORDER = {'euler': 1, 'trap': 2, 'gear': 2, 'theta': 2,
                    'trbdf2': 2, 'radau': 5}

    def grid_error(self, evaluate, refine=2, levels=3, label=None):
        """How much of a scalar is DISCRETISATION rather than answer.

        Re-solves this circuit on a `refine`x finer grid, with every other
        argument identical to the original `solve()`, and reports how far
        `evaluate` moves.  `evaluate` takes a solved `PSS` and returns a
        float -- `lambda p: PAC(cir).diffusion_constant(p)`, a Floquet
        multiplier, a harmonic amplitude, the period.

        ⚠⚠ WHY THIS EXISTS RATHER THAN A PER-METHOD FORMULA.  The floor of
        this stack is DISCRETISATION, it grows LINEARLY IN Q, and it is a
        METHOD property: measured on the analytic high-Q van der Pol
        reference, the relative error in the diffusion constant at 240
        points per period is

            Q      gear        trap        radau
             100   1.79e-03    1.49e-06    6.97e-10
             500   9.02e-03    1.04e-04    3.48e-09
            1000   1.82e-02    2.36e-04    6.97e-09

        i.e. `~1.8e-05 Q` for gear against `~7.0e-12 Q` for radau -- SIX
        ORDERS at the same cost per step.  Those constants are real but they
        belong to THAT fixture at THAT grid: `gear` converges at `O(h^3)` on
        an autonomous problem for `Q >= 5` and at `O(h^2)` at `mu = 1`, so a
        shipped predictor built from them would extrapolate a fitted constant
        across a regime change (roadmap D.0y).  Refining the actual circuit
        measures the actual number instead, and needs no calibration.

        ⚠⚠ WHY THREE GRIDS AND NOT TWO.  With `f_h = f + C h^p`, two grids
        give `|f_h - f_h/r| = |C| h^p (1 - r^-p)`, which over-states the fine
        grid's own error `|C|(h/r)^p` by `r^p - 1` -- an upper bound, and a
        tempting place to stop.  **IT IS NOT SAFE, AND THIS STACK CONTAINS A
        COUNTEREXAMPLE.**  `trap` on the high-Q van der Pol has an error that
        CHANGES SIGN near `Q = 100`: two terms of opposite sign cancel, the
        two-grid difference collapses, and the estimate UNDER-STATES the true
        error by 3.6x (measured: change 1.15e-06 against a true 4.06e-06 at
        240 points).  A bound that fails silently where the error is
        interesting is worse than none.

        So the third grid is not extra confidence, it is the VALIDITY CHECK.
        From `d1 = |f_h - f_h/r|` and `d2 = |f_h/r - f_h/r^2|`,

            order = log(d1/d2) / log(r)

        is the order the circuit ACTUALLY shows, and it is checked against the
        single-power-law assumption before the error estimate built on it is
        offered.  Measured orders on that fixture: `gear` 2.94 (its `O(h^3)`
        autonomous rate), `radau` ~5, and `trap` failing the check exactly
        where it cancels.  `error` is then `d2 / (r^order - 1)`, and
        `power_law=False` means READ `d2` AS A RAW CHANGE AND NOTHING MORE.

        ⚠ AND IT IS AN ESTIMATE OF THE GRID ERROR ONLY.  It cannot see an
        error both grids share -- a wrong stamp, a wrong tolerance
        convention, a mis-specified circuit.  A small `rel_change` says the
        grid is fine enough; it does NOT say the answer is right.  That is
        the same trap `null_residual_amplification` documents one screen up,
        and it is worth stating twice.

        `levels=3` (the default) costs two extra solves, at `r` and `r^2`
        times the points.  `levels=2` is the cheap two-grid difference with
        no order and no validity check -- use it only where the method's
        order on this problem is already known.

        Returns a dict with `values` (coarse to finest), `deltas`, `order`,
        `error` (of the FINEST value), `rel_error`, `power_law`, `refine`
        and `npts`.
        """
        args = getattr(self, '_solve_args', None)
        if args is None:
            raise RuntimeError(
                'PSS.grid_error: call solve() before grid_error() -- the '
                'refinement repeats THIS solve and there is nothing to '
                'repeat yet.')
        refine = int(refine)
        if refine < 2:
            raise ValueError(
                'PSS.grid_error: refine must be >= 2, got %r; a refinement '
                'that does not refine reports 0.0 and means nothing.'
                % (refine,))

        levels = int(levels)
        if levels not in (2, 3):
            raise ValueError('PSS.grid_error: levels must be 2 or 3, got %r'
                             % (levels,))
        ## ⚠ A CALLER-SUPPLIED GRID CANNOT BE REFINED BY A TIMESTEP.  `grid`
        ## fixes the sample fractions outright, so dividing `timestep` would
        ## change nothing and this would report a confident 0.0.
        if args.get('grid') is not None:
            raise ValueError(
                'PSS.grid_error: this solve used an explicit `grid`, whose '
                'fractions fix the samples regardless of `timestep`. The '
                'refinement would return the same grid and report 0.0 -- '
                'pass a refined `grid` and compare directly instead.')

        ## Same class, same parameters, same solve arguments -- only the
        ## timestep changes.  Fresh instances rather than re-solving `self`,
        ## so the caller's solved state survives the call.
        kv = {}
        for _p in self.parameters:
            try:
                kv[_p.name] = getattr(self.par, _p.name)
            except AttributeError:
                pass

        values = [float(evaluate(self))]
        for _k in range(1, levels):
            ## `self.irefnode` is an INDEX; `__init__` wants the node, and
            ## feeding the index back through `get_node_index` would not
            ## round-trip.
            twin = type(self)(self.cir, toolkit=self.toolkit,
                              irefnode=self.cir.nodes[self.irefnode], **kv)
            sub = dict(args)
            sub['timestep'] = args['timestep'] / float(refine ** _k)
            twin.solve(**sub)
            if not getattr(twin, 'converged', False):
                raise RuntimeError(
                    'PSS.grid_error: the %dx refined solve did not converge, '
                    'so there is no comparison to report.'
                    % (refine ** _k,))
            values.append(float(evaluate(twin)))

        deltas = [abs(values[i + 1] - values[i])
                  for i in range(len(values) - 1)]
        _tiny = np.finfo(float).tiny
        order, power_law = None, None
        if levels == 3:
            d1, d2 = deltas
            ## ⚠ `d2 >= d1` means the sequence is NOT settling: either the
            ## error is not a single power law (`trap`'s sign change) or the
            ## finest grid has reached a roundoff floor.  Either way the
            ## Richardson step below would be arithmetic on noise.
            _sgn = ((values[1] - values[0]) * (values[2] - values[1]) > 0.0)
            if d2 > _tiny and d1 > d2:
                order = float(np.log(d1 / d2) / np.log(float(refine)))
                ## ⚠⚠ THE CEILING IS THE POINT OF THIS CHECK, AND A GENERIC
                ## RANGE IS NOT ENOUGH.  A method cannot converge faster than
                ## its order; an observed order well above it means two error
                ## terms nearly cancelled at this grid, which makes the
                ## deltas shrink faster than the error and the estimate
                ## UNDER-state.  MEASURED: `trap` at 120 points shows an
                ## apparent order of 6.45 -- monotone, same-signed deltas,
                ## nothing else suspicious -- while its estimate under-states
                ## the true error by 300x.  A plain `0.5 <= order <= 8` range
                ## ACCEPTS that case; the ceiling below rejects it.
                ## ⚠ The `+ 1.5` allowance is not slack: on an AUTONOMOUS
                ## problem the period is an unknown that absorbs the leading
                ## frequency error, so `gear` (nominal 2) genuinely converges
                ## at 3.01 here.  Without the allowance this would reject the
                ## shipped default method on its own reference fixture.
                _nom = self.METHOD_ORDER.get(
                    str(getattr(self.par, 'method', '')).lower())
                _hi = 8.0 if _nom is None else (_nom + 1.5)
                power_law = bool(0.5 <= order <= _hi and _sgn)
            else:
                order, power_law = None, False
            if not power_law:
                warnings.warn(
                    'PSS.grid_error: the refinement does not follow a single '
                    'power law (changes %.3e then %.3e over %dx refinements, '
                    'implied order %s). The error estimate is WITHHELD -- '
                    'read the raw change and nothing more. This happens when '
                    'two error terms of opposite sign cancel at some grid, '
                    'which makes a two-grid difference UNDER-state the true '
                    'error, and when the finest grid has hit a roundoff '
                    'floor.'
                    % (d1, d2, refine,
                       ('%.2f' % order) if order is not None else 'none'),
                    RuntimeWarning, stacklevel=2)

        if power_law:
            err = deltas[-1] / (float(refine) ** order - 1.0)
        else:
            ## No validated order: report the raw change, which is what a
            ## two-grid call gets and is honest about being unbounded.
            err = deltas[-1]
        rel = err / max(abs(values[-1]), _tiny)
        _npts_c = int(round(args['period'] / args['timestep']))
        return {'values': values, 'deltas': deltas, 'order': order,
                'power_law': power_law, 'error': err, 'rel_error': rel,
                'refine': refine,
                'npts': [_npts_c * refine ** k for k in range(levels)],
                'label': label}

    ## B7: the interpolant degree the defect-correction estimate needs, per
    ## method.  Two-part rule, MEASURED 2026-09-08 (doc/pss_roadmap_260902.md,
    ## B7's gate): the interpolant's DEGREE must exceed the method's stage
    ## count -- a cubic spline lies INSIDE Radau IIA(3)'s collocation
    ## exactness class (degree s = 3), so the neighbouring problem is solved
    ## EXACTLY and the estimate is 1e-10 ppm against a true 6e-3 -- and its
    ## ORDER must exceed the method's effective order, or a constant bias
    ## remains (a quintic against radau's measured 6.1 left 4.9 % at every
    ## grid; a septic gave 1.0001 / 0.9999 / 0.9998).  Cubic reproduced trap
    ## and TR-BDF2 to 0.9996 -> 1.0000.
    ## ⚠ THE STAGE-COUNT CLAUSE IS A COLLOCATION PROPERTY (peer, measured the
    ## same day): ESDIRK43 has SIX stages and is not a collocation method, so
    ## a quintic fails the clause literally -- and quintic and septic AGREE
    ## through the stack (0.9998 / 1.0000 at 100 pts, 1.0000 / 1.0001 at 200).
    ## For a non-collocation method only the ORDER clause is established;
    ## written as "degree > stage count" the rule would over-constrain every
    ## DIRK ever added.  The failure the clause guards against is SILENT (a
    ## clean small number), which is why it was measured rather than argued.
    WARPING_CHECK_TOL = 0.05     # |half-grid / full-grid - 1| above this: the interpolant sets the reading
    IDEC_DEGREE = {'euler': 3, 'trap': 3, 'gear': 3, 'theta': 3,
                   'trbdf2': 3, 'esdirk43': 5, 'radau': 7}

    def warping_estimate(self, periods=20, degree=None, check=True):
        """Estimate THIS solve's period (warping) error at ITS OWN grid, with
        no reference solution and no refinement -- by defect correction.

        B7's answer, MEASURED (2026-09-08).  A per-step local truncation
        estimate cannot see an accumulating period error because warping is a
        GLOBAL error; defect correction (Sickenberger, Weinmueller & Winkler,
        "Local Error Estimates for Moderately Smooth ODEs and DAEs", Part I,
        Sec. 1) estimates the global error directly:

          1. p(t)  -- a periodic spline through this solve's own grid values
                      (`IDEC_DEGREE[method]`, or `degree`);
          2. r(t)  = C(p) p' + i(p) + u(t)   -- the DEFECT of the interpolant
                      against the circuit's own equations, T-periodic;
          3. the NEIGHBOURING problem  d/dt q(y) + i(y) + u(t) - r(t) = 0,
                      whose exact solution is p by construction, integrated as
                      a TRANSIENT with the SAME method at the SAME step for
                      `periods` periods (`Transient.solve` takes `-r` through
                      `provided_function`, an extra source term on every path);
          4. the phase lag of y against p, per period, by projection onto p':
                      tau_k = <p'.(y - p)> / <p'.p'>;  its slope against the
                      period count is the estimated period error.

        A transient and not a periodic solve, deliberately: a forcing at T
        fixes the period, and warping cannot present as a period change.

        Measured on A10's van der Pol (Q = 1e4), estimate / true period error
        (true = T_h - T_ref, radau at 3200 points), numpy prototype:

            trap   / cubic    0.9996  0.9999  1.0000  1.0000   (100..800 pts)
            trbdf2 / cubic    0.9999  1.0000  1.0000
            radau  / cubic    0.0000  0.0000  0.0000   (25..50 pts) -- INSIDE the
                                                       exactness class: see IDEC_DEGREE
            radau  / quintic  1.0490  1.0494  1.0495   -- a constant bias where the
                                                       orders tie (6 vs 6.1)
            radau  / septic   1.0001  0.9999  0.9998

        Controls: with RADAU solving the neighbouring problem of TRAP's defect
        the estimate is 0.0000 -- the drift is the METHOD's error, not the
        defect's; a LINEAR interpolant gives 0.03 / 2.3 / 3.4 -- the
        interpolant-order wall from below.

        Through THIS method (the stack, same fixture, reference radau at 3200
        points, 2026-09-08): trap 1.0002 at 400 and 200 pts; radau septic
        1.0003 / 1.0002 at 50 / 35 pts and CUBIC 0.0001 (the exactness-class
        zero, reproduced); esdirk43 quintic 0.9998 / 1.0000 and septic
        1.0000 / 1.0001 at 100 / 200 pts -- so for a non-collocation method
        the order clause alone is established, and the stage-count clause is
        a collocation property.  ⚠ The first stack gate's driven control came
        back `autonomous=True`: `Circuit.u(t)` evaluates its time functions
        only when told `analysis='tran'`, and without it every source
        VANISHES (zeros, DC value included) -- so that control ran against a
        circuit with NO source at all; fixed at both call sites; with the
        flag the driven van der
        Pol returns `autonomous=False`, `period_error=None`, and a bounded
        lag series, as it must.

        ⚠⚠ THE INTERPOLANT IS THE LIMIT (measured 2026-09-08): on a
        relaxation oscillator with a comparator edge a few points wide the
        estimate reads 0.09 of the true period error at 200 points per
        period and 0.65 at 400 -- uniformly in every component, so not a
        collapse: the septic spline does not resolve the edge and the
        defect is interpolation error, not the method's.  Part I's own
        scope is "moderately smooth"; a relaxation orbit at PSS grids is
        outside it, and the number returned is then wrong by a factor that
        nothing in it announces.  Trust it on smooth orbits (1.000 to four
        digits on the van der Pol, index 1 and 2); on an orbit with edges,
        refine until the estimate converges in `periods` and grid, or use
        `grid_error`.
        ⚠ THE LITERATURE'S ANSWER IS STRUCTURAL, NOT "REFINE" (docs session,
        Part I p. 9, READING-LOG 2.165).  The gate this instrument should
        test before returning a number is Part I's own "only if": the
        estimate is asymptotically correct ONLY IF the interpolant's defect
        error is o(h^{p+1}) -- asymptotically SMALLER than the truncation
        error it is meant to reveal; on the comparator edge it is not, and
        the number is wrong by a factor nothing announces.  And the fix for
        a non-smooth orbit is to form the defect as a WEIGHTED SUM OF
        f-VALUES with an auxiliary scheme sharing the base scheme's
        left-hand side, so the solution terms cancel identically (their eq.
        2.13, an extra factor h) -- not a higher-degree interpolant of the
        solution, which is exactly the construction this one uses.  Scope:
        their construction is the LOCAL error of an LMM; whether it
        transfers to a period functional is unproven.  THE GATE IS BUILT
        (2026-09-08, `check=True`): the same pass through every second
        sample of the same solution, transient still at the solve's step;
        `check_ratio` = half-grid slope / full-grid slope, `trusted` =
        within `WARPING_CHECK_TOL` (5 %) of 1, else a warning and the number
        still returned.  Measured: van der Pol 1.0000 (radau and trap, 50
        and 100 points); the relaxation orbit 0.0056 / 0.72 at 200 / 400
        points under radau and 0.41 at 200 under trap -- the cases that
        read 0.09 / 0.65 of the truth are refused, the smooth case accepted
        with four orders of margin.  ⚠ The prediction "below 0.5 at 400"
        was wrong (0.72): the ratio approaches 1 as the edge resolves, so
        the tolerance is the gate, not the ratio's distance from 0.  The
        restructured (f-value) defect (Part I eq. 2.13; for trap the Milne
        device) was GATED and REFUTED as the edge fix (2026-09-08): smooth
        0.9990 / 0.9998, but on the edge orbit 0.06 / 0.61 / 1.29 at 200 /
        400 / 800 points against the spline route's 0.04 / 0.29 / 0.63 --
        faster with the grid and NOT monotone, so a reading near 1 is
        indistinguishable from a wrong one; the paper's own remedy is mesh
        adaptation.  Not built.  Cost of the check: it doubles the call (a
        second `periods`-long transient).
        ⚠ Scope and limits.  The period reading needs an AUTONOMOUS solve;
        on a driven circuit the lag is bounded (entrained) and `period_error`
        is returned as None with the per-period lag series still filled.
        The prototype ran on a 2-state ODE; on a DAE the differential and
        algebraic components converge at different orders (H&W VI.7), so the
        interpolant threshold binds per component and a component-wise
        exactness collapse would be invisible in this scalar phase drift --
        `component_rms` is returned so a caller can look.  Index-2 is outside
        Part I's stated scope.  Cost: `periods` periods of transient at the
        working grid -- no refinement sweep, no analytic reference.

        Returns a dict: `period_error` (s, signed: positive = this solve's
        period is LONG), `ppm`, `lag` (per-period phase lag, s), `degree`,
        `periods`, `autonomous`, `component_rms` (RMS of y - p per unknown
        over the last period, the global error estimate in state space),
        `check_ratio` and `trusted` (the half-grid self-check above; both
        None with `check=False` or when a slope is not finite),
        `lag_components` (periods x unknowns: the same lag per component --
        every row's slope should equal `period_error`; a row at ~0 while
        the others agree is the per-component exactness-class collapse the
        DAE caveat names; a row of NaN is a component with no motion, e.g.
        a node pinned by a source).
        """
        import numpy as _np
        from scipy.interpolate import make_interp_spline
        if self.waveform is None:
            raise ValueError('warping_estimate needs a solved PSS -- call solve() first')
        T = float(self.period)
        times = _np.asarray(self.waveform[0], dtype=float)
        X = _np.asarray(self.waveform[1], dtype=float)         # (n, m)
        if abs(times[-1] - T) > 1e-12 * T:
            times = _np.r_[times, T]; X = _np.column_stack([X, X[:, 0]])
        X = X.copy(); X[:, -1] = X[:, 0]                        # close the orbit exactly
        method = str(self.par.method)
        k = int(self.IDEC_DEGREE.get(method, 3) if degree is None else degree)
        if X.shape[1] <= k + 1:
            raise ValueError('warping_estimate: %d points per period cannot carry a degree-%d '
                             'periodic spline' % (X.shape[1] - 1, k))
        cir, epar = self.cir, self.epar
        ## ⚠ `analysis='tran'`, on BOTH calls.  `Circuit.u(t)` evaluates a
        ## time function only when told which analysis is asking (`VS.u`:
        ## `elif analysis in timedomain_analyses`); without it every source
        ## VANISHES -- the else-branch returns zeros, and even the DC value
        ## lives inside the gated branch (`timedomain_analyses = ('dc',
        ## 'tran')`).  The first gate's driven control -- an `ISin` on the
        ## van der Pol -- came back `autonomous=True` for exactly that reason,
        ## and the defect would have omitted the drive on a driven circuit.
        _u = lambda t: _np.asarray(cir.u(t, epar, analysis='tran'), dtype=float)
        u0 = _u(0.0)
        autonomous = all(_np.allclose(u0, _u(f * T)) for f in (0.37, 0.71))
        n_per = X.shape[1] - 1

        def _run(times_i, X_i):
            """One defect-correction pass: the periodic spline through
            (times_i, X_i), its defect, the neighbouring transient at the
            SOLVE's step T/n_per -- the method's error at ITS grid is what is
            read, and only the interpolant's grid differs between the main
            pass and the self-check below -- and the per-period lag."""
            p = make_interp_spline(times_i, X_i.T, k=k, bc_type='periodic')
            dp = p.derivative()
            def _defect_source(t):
                tt = t % T
                x = _np.asarray(p(tt), dtype=float); xd = _np.asarray(dp(tt), dtype=float)
                r = (_np.asarray(cir.C(x, epar), dtype=float) @ xd
                     + _np.asarray(cir.i(x, epar), dtype=float) + _u(t))
                return -r
            tr = self._new_transient(self._integrator_for(self.par.method))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                res = tr.solve(tend=periods * T, x0=X_i[:, 0].copy(), timestep=T / n_per,
                               provided_function=_defect_source, fixed_timestep=True)
            ty = _np.asarray(res.sweep_values, dtype=float)
            Y = _np.asarray(res.x, dtype=float)                     # (n, steps+1)
            if Y.shape[0] != X_i.shape[0]:
                Y = Y.T
            P = _np.asarray(p(ty % T), dtype=float).T; dP = _np.asarray(dp(ty % T), dtype=float).T
            E = Y - P
            lag = []; lag_c = []
            for j in range(periods):
                sl = (ty >= j * T - 1e-12 * T) & (ty < (j + 1) * T - 1e-12 * T)
                num = float(_np.sum(dP[:, sl] * E[:, sl])); den = float(_np.sum(dP[:, sl] ** 2))
                lag.append(num / den if den > 0 else _np.nan)
            ## per component: the same projection restricted to one unknown.
                ## A phase shift moves every component by the same lag, so on a
                ## healthy estimate every row's slope equals the period error;
                ## a row reading ~0 while the others read the period error is
                ## the exactness-class collapse on THAT component (the DAE
                ## caveat), invisible in the scalar `lag` above.
                num_c = _np.sum(dP[:, sl] * E[:, sl], axis=1); den_c = _np.sum(dP[:, sl] ** 2, axis=1)
            ## a RELATIVE threshold: a node pinned by a source has a
                ## derivative of pure roundoff (measured 1e-32 rms), and
                ## `den > 0` let it print a ratio of 96 where NaN was meant.
                with _np.errstate(divide='ignore', invalid='ignore'):
                    lag_c.append(_np.where(den_c > 1e-20 * den_c.max(), num_c / den_c, _np.nan))
            lag = _np.asarray(lag); lag_c = _np.asarray(lag_c)
            ok = _np.isfinite(lag)
            slope = float(_np.polyfit(_np.arange(periods)[ok], lag[ok], 1)[0]) if ok.sum() >= 2 else _np.nan
            return slope, lag, lag_c, E, ty

        slope, lag, lag_c, E, ty = _run(times, X)
        ## SIGN: a positive lag means y is AHEAD of p; a LONG period makes y
        ## fall BEHIND, so the period error is minus the slope.
        period_error = -slope if autonomous else None
        last = ty >= (periods - 1) * T - 1e-12 * T
        component_rms = _np.sqrt(_np.mean(E[:, last] ** 2, axis=1))
        ## THE SELF-DIAGNOSTIC (Part I's "only if", built 2026-09-08): the
        ## estimate is the method's error only while the interpolant's own
        ## defect is asymptotically smaller than it, and then it does NOT
        ## depend on the interpolant: the same pass through EVERY SECOND
        ## sample of the same solution (the transient still at the solve's
        ## step) must read the same slope.  Where the interpolation error
        ## dominates -- an edge a few points wide -- the two passes disagree,
        ## and the reading is refused (`trusted=False`, a warning) instead
        ## of returned as a number wrong by a factor nothing announces.
        check_ratio = None; trusted = None
        if check:
            sub = _np.arange(0, n_per + 1, 2)
            if sub[-1] != n_per:
                sub = _np.r_[sub, n_per]
            if len(sub) > k + 1:
                slope2 = _run(times[sub], X[:, sub])[0]
                if _np.isfinite(slope) and _np.isfinite(slope2) and slope != 0.0:
                    check_ratio = float(slope2 / slope)
                    trusted = bool(abs(check_ratio - 1.0) <= self.WARPING_CHECK_TOL)
                    if not trusted:
                        warnings.warn('warping_estimate: the interpolant, not the method, sets '
                                      'this reading (half-grid pass / full-grid pass = %.3f); '
                                      'refine the grid until the two agree' % check_ratio)
        return dict(period_error=period_error,
                    ppm=(period_error / T * 1e6) if period_error is not None else None,
                    lag=lag, lag_components=lag_c, degree=k, periods=periods,
                    autonomous=autonomous, component_rms=component_rms,
                    check_ratio=check_ratio, trusted=trusted)

    def factored_period(self):
        """The converged period's steps, kept factored -- see `FactoredPeriod`.

        Runs the factored traversal ONCE, at the solution, and caches it.
        Lazy on purpose: `_traverse_factored` stores `N` factorisations and
        `N` capacitances (`2 N m^2` doubles -- ~800 MB at m=1002 and 50
        points), which is a bad trade to impose on every `solve` for the
        callers who never ask.

        ⚠ IT RE-TRAVERSES RATHER THAN REUSING THE NEWTON'S FACTORS, and the
        difference is not efficiency.  The last `build` call inside the
        Newton is at the last TRIAL iterate; the converged answer is the one
        after it.  Reusing those factors would give an operator for a
        trajectory near the solution instead of at it -- a small error, in
        the third figure, of exactly the kind a converged answer absorbs
        without complaint.
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.factored_period()
        if getattr(self, '_period_state', None) is None:
            raise RuntimeError(
                'PSS: no period to factor -- call solve() first. '
                '(factored_period() replays the CONVERGED trajectory, so '
                'there has to be one.)')
        if not self.converged:
            raise RuntimeError(
                'PSS: the shooting solve did not converge, so there is no '
                'periodic operating point to linearise about. A '
                'small-signal analysis over a non-solution is not a '
                'meaningful answer -- fix the PSS run first.')
        if self._factored_period_cache is not None:
            return self._factored_period_cache

        solved, x0, xm1, times, hs, T, x0_unknown = self._period_state
        _integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if _integ.is_stage_method():
            ## A self-starting stage method has its own factored map (no opener,
            ## no pair).  Route BY STRUCTURE: a fully-implicit tableau to the
            ## coupled builder, a lower-triangular one (DIRK/ESDIRK) to the
            ## sequential builder -- the coupled block is singular on a DAE for
            ## an explicit first stage, so the two cannot share a path.  A new
            ## method of either family reaches the right builder with no edit
            ## here.
            if _integ.is_fully_implicit():
                fp = self.factored_period_full(x0, T, len(times) - 1)
            else:
                fp = self.factored_period_dirk(x0, T, len(times) - 1)
            self._factored_period_cache = fp
            return fp
        if solved:
            C0, steps, x_last, x_prev = self._traverse_factored(
                x0, xm1, times, hs, T=T)
            fp = FactoredPeriod('solved_history', C0, steps, x_last, x_prev,
                                self, times=times, T=T)
        else:
            opening, steps, x0_out, x_last, _dT = \
                self._traverse_factored_plain(x0, T, times, hs,
                                              open_at_x0=x0_unknown)
            fp = FactoredPeriod('plain', opening, steps, x_last, x0_out,
                                self, times=times, T=T,
                                open_at_x0=x0_unknown)
        self._factored_period_cache = fp
        return fp

    def factored_period_full(self, x0, T, npts, method=None):
        """The factored period map of ANY FULLY-IMPLICIT (FULL) stage method
        about a periodic point `x0` -- the tableau-generic coupled monodromy.

        Integrates the orbit under the stage method (`method`, or the PSS's
        own) on a uniform `npts`-step grid and returns the `m x m` monodromy as
        a `FactoredPeriod(kind='full')`.  Each step is the COUPLED `sm x sm`
        stage solve; the stored per-step object is one factorisation plus the
        entering capacitance `Cn`, and the matvec reads the LAST `m`-block of the
        coupled solution (`x_{n+1} == Y_s`, stiff accuracy).  Radau IIA(3) is
        the first member; any Gauss/Lobatto/higher-Radau tableau reuses this.

        ⚠ FULLY-IMPLICIT ONLY.  A DIRK/ESDIRK (lower-triangular `A`, explicit
        first stage) must NOT come here: its coupled block `[0][0] = C` is
        singular on a DAE -- it takes the sequential `factored_period_dirk`
        path instead.  `factored_period` routes by structure.

        ⚠ SELF-STARTING, SO NO TWIN.  No order-dropped opener, so this does not
        consult `monodromy_twin`.
        """
        if method is None:
            method = getattr(self.par, 'method', 'euler')
        x0 = np.asarray(x0, dtype=float)
        if x0.shape[0] == self.cir.n:
            x0 = np.concatenate((x0[:self.irefnode], x0[self.irefnode + 1:]))
        times = np.linspace(0.0, float(T), int(npts) + 1)
        hs = np.diff(times)
        tr_saved = getattr(self, '_tran', None)
        self._tran = self._new_transient(self._integrator_for(method))
        try:
            steps, x_last, x_prev = self._traverse_factored_full(
                x0, times, hs)
        finally:
            self._tran = tr_saved
        return FactoredPeriod('full', None, steps, x_last, x_prev,
                              self, times=times, T=float(T))

    def factored_period_dirk(self, x0, T, npts, method=None):
        """The factored period map of ANY lower-triangular (DIRK/ESDIRK) stage
        method about a periodic point `x0` -- the tableau-generic SEQUENTIAL
        monodromy builder (`FactoredPeriod(kind='dirk')`).  TR-BDF2 is the first
        member; a new DIRK/ESDIRK of any stage count reuses this.  Self-starting,
        so no twin."""
        if method is None:
            method = getattr(self.par, 'method', 'euler')
        x0 = np.asarray(x0, dtype=float)
        if x0.shape[0] == self.cir.n:
            x0 = np.concatenate((x0[:self.irefnode], x0[self.irefnode + 1:]))
        times = np.linspace(0.0, float(T), int(npts) + 1)
        hs = np.diff(times)
        tr_saved = getattr(self, '_tran', None)
        self._tran = self._new_transient(self._integrator_for(method))
        try:
            steps, x_last, x_prev = self._traverse_factored_dirk(x0, times, hs)
        finally:
            self._tran = tr_saved
        return FactoredPeriod('dirk', None, steps, x_last, x_prev,
                              self, times=times, T=float(T))

    FLOQUET_DENSE_LIMIT = 400
    ## below this a multiplier is an annihilated algebraic
    ## direction, not a mode -- see `floquet_modes`
    FLOQUET_NULL_TOL = 1e-12

    def floquet_modes(self, pss_unused=None, nmodes=None, fp=None):
        """The Floquet pairs `(λ_l, μ_l, p_l(t), q_l(t))` — A9's prerequisite.

        Returns a list of dicts, one per mode, ordered by `|λ|` descending.
        `nmodes=None` returns EVERY non-null mode, and that default is the
        requirement rather than a convenience:

        ⚠⚠ ALL OF THEM ARE REQUIRED, BY THE SOURCE. Traversa & Bonani, IET
        CDS 2011: "The calculation of orbital fluctuations and of the
        phase-orbital correlation within Floquet-based noise analysis of
        autonomous systems requires the availability of ALL the direct and
        adjoint Floquet eigenvectors associated with the noiseless limit
        cycle."  (Cited, not verified here; relayed from the paper.)  An
        earlier cost estimate for A9 -- "a few more Floquet pairs" -- was
        relayed without checking it against that sentence, and is wrong.

        ⚠ TRUNCATION IS LEGITIMATE ONLY WITH A BOUND. Traversa & Bonani,
        TCAD 2013, compute a CHOSEN number of exponents and both
        eigenvector sets for the linearisation of index-1 DAEs around a
        limit cycle -- this formulation -- with the error "proved to tend
        to zero along with the ratio between the norms of the NEGLECTED
        AND RETAINED ROWS".  So passing `nmodes` is allowed, but a caller
        who does owes that ratio as the gate; this routine does not
        compute it.  The dense route below returns everything anyway, so
        at the sizes it serves the question does not arise.


            lam    the Floquet MULTIPLIER, eigenvalue of the monodromy
            mu     the Floquet EXPONENT, `log(λ)/T` (complex)
            u0,v0  right and left eigenvectors at `t = 0`, biorthonormal
                   (`v_k† u_l = δ_kl`)
            p      `p_l(t_j) = Φ(t_j,0) u_l(0) · exp(−μ_l t_j)` — the
                   T-PERIODIC part, sampled on the PSS grid
            q      the adjoint counterpart from the reverse replay
            times  the grid `p` and `q` are sampled on

        ⚠⚠ **WHY THIS EXISTS: `S_yy` NEEDS THE EIGENVECTORS OVER THE
        PERIOD, NOT JUST THE EXPONENTS.** Traversa & Bonani (TCAS-I 2011)
        Lemma 3.5 makes the orbital spectrum a sum of Lorentzians centred
        at `jω₀ + Im{μ_l}` with half-width `|Re{μ_l}| + ½h²ω₀²c`, weighted
        by `C_lhj` (their eq 22) — and `C_lhj` is built from the FOURIER
        COEFFICIENTS of `u_l(t)` and of `v_l(t)ᵀ B(t)`. Their own text is
        explicit that the exponents alone do not order the result: *"a
        major role in the C and D coefficients is also played by the
        Floquet eigenvectors, which could determine large orbital
        fluctuations contributions even when the Floquet exponents are not
        near to zero."* So `|λ₂|` — the only mode information this class
        used to expose — is not sufficient, by the source's own statement.

        ⚠ THE PERIODIC PART IS THE OUTPUT, NOT `Φ(t,0)u(0)`. Floquet's
        theorem says the solution is `p_l(t)exp(μ_l t)` with `p_l`
        T-periodic; the raw propagated vector is not periodic and its
        Fourier series is not the one eq (22) wants. Dividing out
        `exp(μ_l t)` is what makes `p_l(T) = p_l(0)` — which is also the
        gate below, and the only check here that needs no reference.

        ⚠ DENSE, AND REFUSED ABOVE `FLOQUET_DENSE_LIMIT`. The monodromy is
        assembled column by column (`n` matvecs) and diagonalised. That is
        honest for the sizes this is useful at and wrong to hide at larger
        ones: an Arnoldi route would return Ritz VECTORS rather than only
        the Ritz values `ppv()` currently keeps, and is the extension.

        ⚠⚠ AND THAT EXTENSION CONVERGES TO THE PHYSICAL MODE *LAST*, WORSE
        AS Q RISES -- a structural fact, not a measurement (Garcia, Romero
        & Acha 2022, read firsthand by the docs session).  Arnoldi resolves
        the LARGEST-magnitude eigenvalues first; the Ritz route works on
        `A = I - M` and recovers `lam = 1 - theta`, so the physical
        `lam_2 -> 1` maps to `theta_2 -> 0`, the SMALLEST, while the fast
        parasitic modes (`lam ~ 0`) sit at `theta ~ 1` and are resolved
        first.  The separation to resolve is `1/theta_2 ~ Q_lambda`: 3.7,
        16.4, 64.5, 128.5 at `Q_lambda` = 3.18, 15.9, 64, 128.  The
        difficulty scales with the very quantity being measured.  ⚠ This
        is a DIFFERENT Krylov problem from B13's, which measured GMRES
        iterations for the bordered SOLVE `(I - M) w = b` and found them
        independent of Q -- solving a system and extracting its smallest
        eigenvalue are not the same question, and B13 says nothing about
        the second.  The paper is a sound source for the method and was
        validated on power networks, not RF oscillators, so it reports no
        evidence either way about the high-Q regime; and it states that a
        truncated run "cannot compute ALL the Floquet multipliers" -- which
        is the requirement eq (22) carries (IET CDS 2011, above).
        """
        _tw = self.monodromy_twin()
        if _tw is not self:
            return _tw.floquet_modes(pss_unused, nmodes, fp)
        fp = pss_unused.factored_period() if fp is None else fp
        n = fp.width
        T = float(fp.T)
        if n > self.FLOQUET_DENSE_LIMIT:
            raise NotImplementedError(
                'PSS.floquet_modes: the monodromy is assembled densely and '
                'this one is %d wide, past the %d limit. The extension is '
                'an Arnoldi that keeps its Ritz VECTORS -- ppv() already '
                'builds the basis and discards them.'
                % (n, self.FLOQUET_DENSE_LIMIT))

        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        lam, U = np.linalg.eig(M)
        lam_l, V = np.linalg.eig(M.T)

        order = np.argsort(-np.abs(lam))
        lam, U = lam[order], U[:, order]
        ## pair each right eigenvalue with its left partner by value
        pair = []
        used = set()
        for k in range(len(lam)):
            d = np.abs(lam_l - lam[k])
            for j in np.argsort(d):
                if j not in used:
                    used.add(int(j))
                    pair.append(int(j))
                    break
        V = V[:, pair]

        times = np.asarray(fp.times, dtype=float)
        out = []
        ## ⚠ NULL MODES ARE DROPPED, NOT RETURNED WITH A BAD RESIDUAL. A
        ## DAE's monodromy has exact zeros (the algebraic directions the
        ## step map annihilates); their "eigenvectors" are arbitrary, the
        ## exponent `log(0)` does not exist, and the periodic part comes
        ## back as noise -- measured, residual 0.56 and periodicity 8.9
        ## against 1e-15 for the physical pair. Returning them invites a
        ## caller to average over a mode that means nothing.
        keep = [k for k in range(n) if abs(lam[k]) > self.FLOQUET_NULL_TOL]
        if nmodes is not None and int(nmodes) < len(keep):
            ## ⚠⚠ TRUNCATING BY MULTIPLIER MAGNITUDE IS REFUTED BY THE SOURCE'S
            ## OWN WORKED EXAMPLE.  Traversa & Bonani TCAS-I 2011 Sec. V, on
            ## their Colpitts: "six orders of magnitude separate mu_2 and mu_3,
            ## while the corresponding contribution to orbital noise are not in
            ## the same ratio.  Rather, far from the oscillator harmonics, the
            ## contribution of mu_3 is dominant with respect to mu_2".  The
            ## ordering INVERTS.  Four statements agree: eq (8)'s sum over
            ## k = 2..n (structural), p.4 (asserted), Sec. V (measured on a
            ## real circuit), this repo's concentration sweep (m/n = 0.97).
            ## The caller who truncates owes the dropped weight as a gate;
            ## this says so at the call rather than only in the docstring.
            warnings.warn(
                'PSS.floquet_modes: nmodes=%d keeps %d of %d non-null modes, '
                'selected by |lambda|. Orbital-noise weight does NOT follow '
                'multiplier magnitude -- Traversa & Bonani (TCAS-I 2011, '
                'Sec. V) show the contribution ordering INVERTING across six '
                'orders in mu on their Colpitts, and this repo measured no '
                'concentration (m/n = 0.97). A covariance or spectrum built '
                'from a truncated set is missing weight you have not bounded; '
                'pass nmodes=None for all modes.'
                % (int(nmodes), int(nmodes), len(keep)),
                RuntimeWarning, stacklevel=2)
        ## ⚠ `None` means ALL non-null modes -- the default since the IET CDS
        ## 2011 correction -- and it used to fall into `int(None)` here because
        ## the only test passed a number. A default nobody exercises is not a
        ## default.
        for k in (keep if nmodes is None else keep[:int(nmodes)]):
            lk = complex(lam[k])
            uk, vk = U[:, k].astype(complex), V[:, k].astype(complex)
            nrm = complex(np.vdot(vk, uk))
            if abs(nrm) < 1e-30:
                raise ValueError(
                    'PSS.floquet_modes: mode %d has left and right '
                    'eigenvectors orthogonal to each other (v.u = %.3e), so '
                    'it cannot be biorthonormalised. That happens at a '
                    'defective eigenvalue -- two multipliers have collided.'
                    % (k, abs(nrm)))
            vk = vk / np.conj(nrm)                     ## v_k† u_k = 1
            muk = np.log(lk) / T

            ## ⚠ EVERYTHING BELOW IS THE WIDTH-`m` STATE BLOCK, NOT THE
            ## WIDTH-`n` MAP INPUT. Under a solved-history map `n = 2m`
            ## and the second block is the history term, not a second
            ## state; the replays collect `m`-wide states either way, and
            ## `u_l(t)` in eq (22) is a state-space function. Mixing the
            ## two is a shape error that surfaces three frames away.
            m = self.cir.n - 1

            ## forward: Phi(t_j,0) u_k(0), by an UNFORCED driven replay
            zero = np.zeros(m)
            _end, fwd = self._forced_replay(fp, 0.0, zero, y0=uk, collect=True)
            traj = ([np.asarray(uk, dtype=complex)[:m]]
                    + [np.asarray(z, dtype=complex).ravel()[:m] for z in fwd])
            tt = times[:len(traj)]
            traj = traj[:len(tt)]
            p = np.column_stack([traj[j] * np.exp(-muk * tt[j])
                                 for j in range(len(traj))])

            ## adjoint: Phi(T,s_j)^T v_k(T) -- B8 made this available under
            ## every integrator, not only the solved-history one
            _e2, _ts, st = fp.matvec_transposed(vk, collect=True)
            qtraj = ([np.asarray(z, dtype=complex).ravel()[:m] for z in st]
                     + [np.asarray(vk, dtype=complex)[:m]])
            ts2 = times[:len(qtraj)]
            qtraj = qtraj[:len(ts2)]
            q = np.column_stack([qtraj[j] * np.exp(muk * ts2[j])
                                 for j in range(len(qtraj))])

            ## ⚠⚠⚠ THE REPLAYED VECTOR IS `C^T q`, NOT `q`.  The conserved
            ## bilinear form of the variational DAE is `w^T C delta`, so over
            ## a period `M_a^T C M = C`, which makes the LEFT eigenvector of
            ## the state monodromy `C^T w(0)` -- the adjoint mode in the
            ## "left-eigenvector coordinates", one factor of `C^T` away from
            ## the state-space adjoint `q` that eq (22) and every covariance
            ## here need.  The transposed replay propagates that object, so
            ## every sample of `q` above is `C(t)^T q_true(t)`.
            ##
            ## ⚠ INVISIBLE ON EVERY FIXTURE THIS REPO HAD, for a geometric
            ## reason: van der Pol's reduced `C` is `diag(1, -1)`, and at
            ## `t = 0` the orbit sits at `[2, 0]` where the adjoint is nearly
            ## axis-aligned, so `C^T q` and `q` point the same way up to sign
            ## (`|cos| = 0.9972`).  On an ASYMMETRIC orbit the seed is off-axis
            ## and the two separate -- measured `|cos(v_k, q_true)| = 0.5738`
            ## at `a = 0.30` on van der Pol + `a u^2` -- while the two adjoints
            ## there are nearly PARALLEL (`|cos(q_2, q_1)| = 0.997`), so the
            ## wrong vector is mostly phase adjoint.  Result: the orbital
            ## covariance was 81x LOW against a Monte Carlo (0.0123 of the
            ## truth), and `|cos(C^-T v_k, q_true)| = 1.0000` at both
            ## asymmetries.  Applying `C^-T` here takes it to 1.06 of the
            ## Monte Carlo at `a = 0.30` and 1.0004 at `a = 0`.
            ##
            ## ⚠ THIS ALSO DISSOLVES THE "q IS STORED IN REVERSE TIME"
            ## finding recorded the same day: with the right vector,
            ## `q(t)^T C p(t)` is conserved at the SAME index (4.2e-04) and
            ## NOT the reversed one (2.0).  `diag(1, -1)` flips one
            ## component, which on a half-wave symmetric orbit is exactly the
            ## relation between `q(t)` and `q(T - t)` -- a sign flip read as a
            ## time reversal.
            ##
            ## ⚠ Per sample, because `C` may depend on the state.  `pinv`
            ## rather than `inv` so a singular reduced `C` (an index-2 MNA,
            ## algebraic rows) does not raise; the algebraic components of `q`
            ## are then the minimum-norm choice, which is a SCOPE LIMIT and
            ## not a solution -- recorded, not hidden.
            _Wq = np.delete(np.asarray(self.waveform[1], dtype=float),
                            self.irefnode, axis=0)
            _nw = _Wq.shape[1]
            for _j in range(q.shape[1]):
                _Cj = np.asarray(self._C_at(_Wq[:, min(_j, _nw - 1)]),
                                 dtype=float)
                q[:, _j] = np.linalg.pinv(_Cj.T) @ q[:, _j]

            ## ⚠⚠ RENORMALISE ON THE STATE BLOCK. `v_k` was biorthonormalised
            ## against `u_k` at the map's FULL width `n`; under a
            ## solved-history map that is the pair `[x_n; x_{n-1}]`, and the
            ## width-`m` state block then carries `q(0)^T p(0) = c0 != 1`.
            ## MEASURED on van der Pol under gear: c0 = 1.324143, constant
            ## around the cycle to four digits -- and the orbital covariance
            ## assembled from these parts came out too large by EXACTLY
            ## c0^2 = 1.7535 against two independent routes, because `q`
            ## enters it quadratically. The periodicity gate p(T) = p(0)
            ## cannot see this: periodicity is scale-free. On the plain path
            ## n = m and c0 = 1, so this is a no-op there. The adjoint takes
            ## the scale (the right vector is the physical direction).
            ## ⚠⚠ THE INNER PRODUCT IS `C`-WEIGHTED, AND THE UNWEIGHTED ONE
            ## WAS WRONG BY A FACTOR OF `C` -- INVISIBLE ON EVERY FIXTURE
            ## THIS REPO HAD.  The conserved bilinear form of the variational
            ## DAE is `q(t)^T C(t) p(t)`, not `q(t)^T p(t)`: differentiating
            ## `G p + d(C p)/dt = 0` against the adjoint gives
            ## `d/dt [q^T C p] = 0`, so `q^T C p` is the invariant and the
            ## biorthonormality that eq (22) assumes is `q_k^T C p_l = d_kl`.
            ##
            ## ⚠ ON A UNIT-REACTANCE FIXTURE THE TWO ARE THE SAME NUMBER,
            ## which is exactly why this survived: van der Pol with
            ## `c = L = 1` gives `q^T C p = 0.9992` against `q^T p = 1`.
            ## Sweep the capacitance at fixed `w0` and the two separate --
            ## MEASURED `q^T C p` = 0.2495 / 0.9992 / 3.9982 at
            ## `C` = 0.25 / 1 / 4, i.e. exactly `C`, while `q^T p` stayed
            ## pinned at 1.000000.
            ##
            ## ⚠⚠ AND `q` ENTERS THE COVARIANCE QUADRATICALLY, so the orbital
            ## covariance came out too large by exactly `C^2`.  Measured
            ## against the independent Lyapunov reference before the fix:
            ## ratio 0.0624 / 1.0001 / 16.043 at those same `C` -- right ONLY
            ## at `C = 1`, which is the only place A9's three-way gate ever
            ## ran.  §D 0c, on the very circuit that produced that entry: a
            ## unit reactance makes `C` the identity and the two inner
            ## products indistinguishable.
            ##
            ## The pair-slicing correction this block was written for is
            ## subsumed: normalising on `q^T C p` fixes the slice scale and
            ## the weighting in one step.
            _x0r = np.delete(np.asarray(self.waveform[1], dtype=float)[:, 0],
                             self.irefnode)
            _Cm = np.asarray(self._C_at(_x0r), dtype=float)
            c0 = complex(np.vdot(q[:, 0], _Cm @ p[:, 0]))
            if abs(c0) < 1e-30:
                raise ValueError(
                    'PSS.floquet_modes: mode %d has q(0)^T C p(0) = %.3e on '
                    'the state block, so it cannot be biorthonormalised '
                    'there.' % (k, abs(c0)))
            q = q / np.conj(c0)
            out.append({'lam': lk, 'mu': muk, 'u0': uk, 'v0': vk,
                        'p': p, 'q': q, 'times': tt, 'c0': c0,
                        'residual': float(np.linalg.norm(M @ uk - lk * uk)
                                          / max(abs(lk), 1e-300))})
        return out

    def _forced_replay_full(self, fp, freq, u_ac, y0=None, collect=False):
        """One driven period under Radau IIA(3) -- the FORWARD coupled replay,
        the transpose of :meth:`_forced_replay_transposed_full`.

        Each step maps the entering state and the source to the endpoint by the
        coupled stage solve: with the source at frequency `freq` entering stage
        `k` (abscissa ``t_{n,k} = t_n + c_k h``),

            rhs_i = C_n y_n - h sum_k A_ik u exp(jw t_{n,k})
            y_{n+1} = (J_block^{-1} rhs)_3            (stiff accuracy)

        so ``y_end = M y0 + w(freq)`` by linearity, the superposition PAC
        relies on.  ⚠ THE SOURCE COUPLING IS THE SAME ``-h sum_k A_ik exp(...)``
        the adjoint reads, so the two are exact transposes (dual-consistent to
        machine precision); a real ``J_block`` factor takes a complex rhs as
        two back-substitutions.
        """
        import scipy.linalg as sla
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        Amat, _Bw, cvec = integ.butcher()
        s = Amat.shape[0]
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        tms = np.asarray(fp.times, dtype=float)
        y = (np.zeros(m, dtype=complex) if y0 is None
             else np.asarray(y0, dtype=complex).ravel().copy())

        def csolve(lu, b):
            return (sla.lu_solve(lu, b.real) + 1j * sla.lu_solve(lu, b.imag))

        ys = []
        for j, (lu, Cn, _mm) in enumerate(fp.steps):
            ts = tms[j]; te = tms[j + 1]; h = te - ts
            cw = Cn @ y
            rhs = np.concatenate([cw] * s).astype(complex)
            for i in range(s):
                srci = -h * sum(Amat[i, k] * u_ac * np.exp(jw * (ts + cvec[k] * h))
                                for k in range(s))
                rhs[i * m:(i + 1) * m] += srci
            Z = csolve(lu, rhs)
            y = Z[(s - 1) * m:s * m]
            if collect:
                ys.append(y.copy())
        return y, ys

    def _forced_replay(self, fp, freq, u_ac, y0=None, collect=False):
        """One period of the LINEARISED circuit, driven at `freq`.

        The same recursion as `_monodromy_matvec`, with `source` switched
        on: `_step_sensitivity` is linear in `(y0, source)`, so this returns

            y_end = M y0 + w(freq)

        with `w` the particular response from a zero initial state.  That
        superposition is not an incidental property -- it is what lets PAC
        solve an `m x m` system instead of an `(N m) x (N m)` one, and it is
        asserted in the suite rather than assumed.

        ⚠ THE SOLVE IS REAL, THE REPLAY IS COMPLEX.  `Jf` is real and stays
        factored once; a complex right-hand side costs two back-substitutions
        against those same factors.  See the note in `_monodromy_matvec`.
        """
        if fp.kind == 'dirk':
            return self._forced_replay_dirk(fp, freq, u_ac, y0, collect)
        if fp.kind == 'full':
            return self._forced_replay_full(fp, freq, u_ac, y0, collect)
        m = self.cir.n - 1
        jw = 2j * np.pi * float(freq)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()

        if fp.kind == 'solved_history':
            if y0 is None:
                Px = [np.zeros(m, dtype=complex), np.zeros(m, dtype=complex)]
            else:
                y0 = np.asarray(y0, dtype=complex).ravel()
                Px = [y0[:m].copy(), y0[m:].copy()]
            Cs = list(fp.opening)
            Pq = np.zeros(m, dtype=complex)
        else:
            C_open, a_open, b_open, pq_open = fp.opening
            v = (np.zeros(m, dtype=complex) if y0 is None
                 else np.asarray(y0, dtype=complex).ravel().copy())
            Px = [v.copy(), v.copy()]
            Cs = [C_open, C_open]
            ## the OPENING pair, exactly as `_monodromy_matvec_plain` -- see
            ## the note there; using the loop's makes it wrong for `trap`
            Pq = (a_open[0] * (C_open @ v) if b_open
                  else np.zeros(m, dtype=complex))
            ## ⚠ AND THE CONSISTENT-`iq_0` SEED, for the same reason and by
            ## the same term: `y0` perturbs `x_0`, and a method that seeds
            ## `iq_0 = -(i(x_0) + u(t_0))` carries that perturbation into the
            ## companion current before the first step.  `None` for every
            ## method that does not declare `needs_consistent_iq0`.  See
            ## `_pq_seed_at_x0`; THE FORCED REPLAY MUST MATCH THE MONODROMY
            ## it superposes with, or `y_end = M y0 + w` stops holding.
            if pq_open is not None:
                Pq = Pq + pq_open @ v

        ys = []
        for (lu, C_new, alphas, b), t in zip(fp.steps, fp.times[1:]):
            src = u_ac * np.exp(jw * float(t))
            Px_new, Pq = self._step_sensitivity(
                Px, Cs, Pq, None, C_new,
                solve=lambda S, _l=lu: _complex_solve(_l, S),
                coeffs=(alphas, b), source=src)
            Px = [Px_new, Px[0]]
            Cs = [C_new, Cs[0]]
            if collect:
                ys.append(Px_new.copy())

        end = (np.concatenate((Px[0], Px[1])) if fp.kind == 'solved_history'
               else Px[0])
        return end, ys

    ## How hard GMRES is asked to solve, relative to the shooting tolerance.
    ## An inexact Newton only needs the step accurate enough not to spoil the
    ## outer convergence; measured k is 2-12 on circuits whose `I - M`
    ## clusters at 1 (the fast modes decay over a period, leaving the slow
    ## ones), so k tracks the number of SLOW MODES, not m.
    KRYLOV_TOLERANCE_FACTOR = 1e-2
    ## ⚠ THE BUDGET IS A CHOICE AND SCIPY'S UNITS ARE A TRAP: `maxiter` counts
    ## RESTART CYCLES, not matvecs, so the pair multiplies. 200 x 20 is far
    ## more than a clustered system needs; a circuit that exceeds it does not
    ## cluster, and the answer is the dense path, not a bigger budget.
    KRYLOV_RESTART = 200
    KRYLOV_MAX_CYCLES = 20

    def _matrix_free_solve(self, z0, times, hs, abstol, xtol, reltol,
                           maxiter):
        """The outer Newton with the monodromy never formed (item 6).

        The dense path builds the `2m x 2m` Jacobian and factors it once per
        iteration; here the same iteration runs on a matvec, so the
        `2m`-column propagation never happens.  Measured against the dense
        path on the RC ladder, single-threaded, k=12:

              m     dense traversal   trajectory + 12 matvecs   speedup
             40             0.0843                    0.1025      0.82x
            110             0.2366                    0.2175      1.09x
            242             0.7503                    0.5378      1.40x
            502             3.4709                    1.5457      2.23x
           1002            20.1143                    5.5512      3.62x

        -- 82-87% of the predicted ceiling, and a LOSS at m=40, which the
        ceiling said too.

        ⚠ THOSE ARE TRAVERSAL FIGURES AND THE END-TO-END SOLVE GAINS LESS.
        A `solve` also does its setup, the replay that builds the waveform
        and the DFT, none of which this touches, and matrix-free spends an
        extra Newton iteration (below).  Measured end to end, same circuits:

              m    dense (iters)      matrix-free (iters)     speedup
            242      2.113 s (2)            1.557 s (2)        1.36x
            502      9.255 s (2)            6.131 s (3)        1.51x
           1002     52.402 s (2)           24.636 s (3)        2.13x

        Quote whichever answers the question being asked, and say which it
        is; 2.23x and 1.51x at m=502 are both true and are not the same
        measurement.

        ⚠ THE CONVERGENCE TEST IS NOT BIT-IDENTICAL TO `analysis.fsolve`'s,
        and it cannot be.  `fsolve` scales its residual test by
        `|J| . |x|`, an ELEMENTWISE absolute value of the Jacobian, which no
        matrix-free method has.  The substitute here is `|x| + |M x| + |F|`,
        one extra matvec per iteration.

        ⚠ AND IT IS NOT PROVABLY THE STRICT DIRECTION.  This docstring first
        claimed the substitute was a LOWER bound on `fsolve`'s scale, so
        that the test could only ever be stricter.  That is FALSE: at
        `M = I` the true scale `|I - M| . |x|` is zero while the substitute
        is `2|x|`, so the substitute is the LARGER one there, and at `M = 0`
        they are equal.  Neither dominates the other in general.

        What is measured, on the RC ladder at m=242/502/1002: the two paths
        agree on the converged waveform to 1.1e-16 and on the converged/not
        verdict, and matrix-free takes ONE MORE Newton iteration at m>=502
        (3 against 2) -- so it is stricter in practice here, and still wins
        on wall time while doing 50% more traversals.  One circuit is not a
        proof of direction, and this is the first thing to check if the two
        paths ever disagree on convergence.
        """
        m = self.cir.n - 1

        def build(z):
            C0, steps, x_last, x_prev = self._traverse_factored(
                z[:m], z[m:], times, hs)
            F = np.concatenate((z[:m] - np.asarray(x_last),
                                z[m:] - np.asarray(x_prev)))
            return F, (lambda v: v - self._monodromy_matvec(C0, steps, v))

        return self._matrix_free_newton(build, z0, abstol, xtol, reltol,
                                        maxiter)

    def _matrix_free_newton(self, build, z0, abstol, xtol, reltol, maxiter):
        """The Newton loop every matrix-free system shares.

        `build(z)` returns `(F, matvec)` for the current iterate -- one
        trajectory pass, then a linear operator that never forms its matrix.
        Written once because the four systems differ ONLY in those two
        things: the plain path's `I - M`, the solved-history path's `2m`
        pair, and the bordered autonomous versions of each.
        """
        import scipy.sparse.linalg as spla
        z = np.asarray(z0, dtype=float).copy()
        n = len(z)
        ier, mesg, xdiff = 2, 'No convergence', None
        for _i in range(maxiter):
            F, mv = build(z)

            def _mv(v, _f=mv):
                return _f(v)

            J = spla.LinearOperator((n, n), matvec=_mv, dtype=float)
            xdiff, info = spla.gmres(
                J, -F, rtol=self.KRYLOV_TOLERANCE_FACTOR * reltol,
                restart=min(n, self.KRYLOV_RESTART),
                maxiter=self.KRYLOV_MAX_CYCLES)
            ## ⚠ THE INNER SOLVE'S VERDICT IS NOT DISCARDED.  It used to be,
            ## and a Krylov breakdown then surfaced as the generic outer
            ## 'No convergence' with nothing naming the cause -- in a file
            ## whose whole standard is that a failure says what happened
            ## (`T = 0`, the trivial root, the singular free-period
            ## Jacobian).  An unconverged GMRES makes `xdiff` a direction
            ## the Newton has no reason to trust, so the outer loop is told
            ## to stop rather than iterate on it.
            if info != 0:
                warnings.warn(
                    'PSS: the matrix-free inner solve did not converge at '
                    'outer iteration %d -- GMRES returned info=%d (%s) on a '
                    '%d-unknown system, after at most %d matvecs, each a '
                    'full replay of the period. The Newton step it returned '
                    'is not a direction worth iterating on, so this solve '
                    'stops here and reports not-converged. Measured k on '
                    'well-behaved circuits is 2-12 because `I - M` clusters '
                    'at 1; needing more than %d means this system does not '
                    'cluster, and the dense path (matrix_free=False) is the '
                    'reliable answer for it.'
                    % (_i, info,
                       'breakdown' if info < 0 else 'iteration limit',
                       n, min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES,
                       min(n, self.KRYLOV_RESTART) * self.KRYLOV_MAX_CYCLES),
                    RuntimeWarning, stacklevel=3)
                return z, {}, 2, 'No convergence (inner Krylov solve failed)'
            z_new = z + xdiff

            ## `|J| . |x|` is not available without the matrix; see
            ## `_matrix_free_solve` for what this substitute is and is not.
            I_scale = np.abs(z_new) + np.abs(mv(z_new)) + np.abs(F)
            conv_x = np.all(np.abs(xdiff)
                            < reltol * np.maximum(np.abs(z_new), np.abs(z))
                            + xtol)
            conv_f = np.all(np.abs(F) < reltol * I_scale + abstol)
            z = z_new
            if conv_x and conv_f:
                ier, mesg = 1, 'Success'
                break
        return z, {}, ier, mesg

    ## Below this, a multiplier says the mode decays by six decades in one
    ## period and no stability question turns on it -- so parasitic
    ## contamination at that level is not worth a warning.  Used only to
    ## keep the warning off circuits where the WHOLE spectrum is numerical
    ## noise; it never changes a reported number.
    SPECTRAL_NOISE_FLOOR = 1e-6

    def _spectral_report(self, M):
        """Split a composed spectrum into physical multipliers and parasitics.

        RECORDED SCOPE ITEM 3.  A k-step method turns an m-dimensional
        system into a k*m-dimensional discrete one, so the composed
        monodromy's spectrum carries `(k-1) m` PARASITIC roots beside the
        physical Floquet multipliers.  `max |eig|` over that mixture is only
        a stability verdict while the parasitic roots stay small -- which
        for Gear-2 they emphatically do (`(1/3)^N`, ~1e-95 at 200 points)
        and for a method whose spurious root sits nearer the unit circle
        they would not.  This separates them instead of hoping.

        THE DISCRIMINATOR IS THE EIGENVECTOR'S BLOCK STRUCTURE, not the
        eigenvalue.  The composed map acts on the PAIR `(x_0, x_{-1})`:

          - a PHYSICAL mode follows the linearised ODE, so its two halves
            are one timestep apart on a smooth trajectory and
            `v_{-1} = e^{-lambda h} v_0 -> v_0` as `h -> 0`;
          - a PARASITIC mode is `r^n u` for the method's spurious root `r`,
            so `v_{-1} = u / r` -- three times `v_0` for BDF-2, minus it for
            a trapezoidal-like root -- and the halves differ by O(1)
            whatever `h` is.

        So `||v_{-1} - v_0||` (against a unit-norm eigenvector) is O(h) for
        a physical mode and O(1) for a parasitic one.  MEASURED, and it is
        the h-scaling that makes it a prediction rather than a story: on the
        phase circuit the physical ratio falls 0.1281 -> 0.0316 when the
        grid goes from 50 to 200 points -- a factor of 4.05 for a factor of
        4 in `h` -- while the parasitic ratios sit at 1.0 to 10.  On the
        Q=20 RLC the parasitic ratio is 1.9997 against the 2.0 that BDF-2's
        `v_{-1} = 3 v_0` predicts exactly.

        ⚠ THE MODE COUNT HERE IS AN ODE COUNT AND THE OBJECT IS A DAE, and
        the difference is structural rather than an off-by-`k`.  Demir
        (IJCTA 28:163-185, 2000) gives the DAE monodromy as

            Phi(t,s) = U(t) D(t-s) V(s) C(s)

        with `D = diag[exp(mu_1 (t-s)), ..., exp(mu_d (t-s)), 0, ..., 0]`
        for `d = rank(C)`: "equation (19) has k = n - m Floquet multipliers
        that are 0", and on a real circuit "there are also eigenvalues
        exactly equal to 0 due to the ALGEBRAIC EQUATIONS in the MNA
        formulation".  So the `m - rank(C)` structural zeros are the
        theory's, not an artefact -- which is why `parasitic_roots` comes
        back identically zero on every MNA circuit tried here.

        ⚠ AND THE FACTORISATION CARRIES A TRAILING `C(s)` WITH NO ODE
        ANALOGUE (where `C = I` and it disappears).  A DAE monodromy is not
        simply a product of state-transition blocks, so an ODE-shaped
        count does not merely miscount -- it describes a different object.
        Anyone revisiting this split should start there and not from the
        eigenvector heuristic below.  Relayed from the docs session's read;
        check it against the paper before building on it.

        ⚠ THE SPLIT IS BY RANK, NOT BY A THRESHOLD, and that was measured
        into the design rather than chosen.  A threshold of 0.25 was tried
        first and returned NO physical modes at all on a stiff RC ladder --
        `lambda h ~ 40` there, so every mode's halves differ by O(1) and the
        classifier called the entire spectrum parasitic, handing back a
        `spectral_radius` of `None` where the old code said 6e-15.  A
        `k`-step method on `m` states has EXACTLY `m` physical multipliers
        and `(k-1) m` spurious ones -- that is structural -- so the `m`
        smallest splits are the physical set by construction, and the
        question of where to put a cut never arises.

        ⚠ THE COUNT IS AN ODE COUNT, AND MNA CIRCUITS ARE DAEs.  This
        splits `2m` eigenvalues as `m` physical and `m` parasitic, which is
        right for an ODE.  An index-1 MNA system with `d = rank(C) < m` has
        `d` physical multipliers, `d` parasitic ones and `2(m - d)`
        STRUCTURAL ZEROS from the algebraic variables -- so on a real
        circuit both arrays are mislabelled: measured on the Q=20 resonator
        (`m = 4`, `rank(C) = 2`), `parasitic_roots` comes back identically
        zero and `floquet_multipliers` carries two structural zeros beside
        the two real multipliers.

        ⚠ `spectral_radius` IS UNAFFECTED, which is why this is recorded
        rather than re-engineered.  The physical multipliers have the
        SMALLEST block split by construction, so they are always inside the
        first `m`, and the maximum over that set is the right number --
        0.97531 on that circuit, against the analytic 0.9753.  What is
        unreliable is the LABELLING of the diagnostic arrays.  And it cannot
        be fixed by magnitude either: Gear-2's true parasitic roots are
        `(1/3)^N`, about 1e-95, which is numerically indistinguishable from
        a structural zero -- so on this method the two populations cannot be
        told apart at all, by any test, and saying so is the honest
        position.

        ⚠ ON A STIFF CIRCUIT THE LABELS MAY STILL BE WRONG, and it does not
        matter: when the physical modes are themselves stiff, a parasitic
        root can have the smaller split and swap places with one.  Every
        mode involved then has `|mu|` at the noise floor, so the RADIUS is
        unaffected -- it is the labels, not the number, that degrade.  What
        this buys is the case that motivated the item: a method whose
        spurious root sits NEAR THE UNIT CIRCLE, where the physical modes
        are well resolved, the splits separate cleanly, and taking a
        maximum over the mixture would report the discretisation's own
        artefact as the orbit's stability.

        Returns `(rho, physical, parasitic)`: the spectral radius over the
        PHYSICAL multipliers only, and both sets sorted by magnitude.  An
        `m x m` monodromy (any one-step method, the plain path) has no pairs
        and no parasitic roots, so everything in it is physical.  `None`
        gives `(None, None, None)` -- the matrix-free path forms no
        monodromy at all.
        """
        if M is None:
            return None, None, None
        M = np.asarray(M)
        m = self.cir.n - 1
        try:
            ev, V = np.linalg.eig(M)
        except np.linalg.LinAlgError:                     # pragma: no cover
            return None, None, None

        if M.shape[0] != 2 * m:
            ## one-step method: the monodromy IS the physical map
            phys = np.sort(np.abs(ev))[::-1]
            return float(phys[0]), phys, np.array([])

        ## columns of `V` are unit-norm, so this needs no denominator and
        ## cannot divide by a vanishing block -- a mode living entirely in
        ## one half reads as O(1) here, which is what it is.
        split = np.linalg.norm(V[m:, :] - V[:m, :], axis=0)
        order = np.argsort(split)
        phys = np.sort(np.abs(ev[order[:m]]))[::-1]
        para = np.sort(np.abs(ev[order[m:]]))[::-1]
        rho = float(phys[0])
        ## ⚠ THE POINT AT WHICH THIS STOPS BEING BOOKKEEPING.  While the
        ## parasitic roots are 80 orders down, separating them changes no
        ## number and only documents why the maximum was safe.  Once one
        ## climbs to within a decade of the physical spectrum, the method's
        ## spurious roots are a real part of what the analysis reports and
        ## the user is entitled to know before reading a stability verdict.
        if len(para) and para[0] > 0.1 * rho and rho > self.SPECTRAL_NOISE_FLOOR:
            warnings.warn(
                'PSS: this method\'s PARASITIC roots are no longer '
                'negligible -- the largest is %.4g against a physical '
                'spectral radius of %.4g. A k-step method contributes '
                '(k-1)*m spurious roots to the composed monodromy, and '
                '`spectral_radius` now reports the maximum over the '
                'PHYSICAL multipliers only (separated by eigenvector block '
                'structure). Treat the separation as load-bearing here '
                'rather than cosmetic: check `floquet_multipliers` and '
                '`parasitic_roots` before drawing a stability conclusion.'
                % (para[0], rho), RuntimeWarning, stacklevel=3)
        return rho, phys, para

    def _install_history(self, x0_in, xm1_in, dt, h_prev=None):
        """Open a run ON a solved two-point history rather than a seed.

        `_begin_run(x_{-1})` opens the rings on the earlier point and the
        push puts `x_0` in front of it, so the first real step reads
        `q(x_0)` and `q(x_{-1})` -- two genuine solved points.  The flags
        then say what is true of them: a step of `dt` has been taken, and
        the run is no longer opening, so nothing drops order.

        `_dt_last2` stays None on purpose.  The THIRD charge in the ring is
        still `q(x_{-1})` repeated, and the LTE estimator differences three,
        so its opening reading remains unsound and the report goes on
        discarding it -- a solved history fixes what the SOLUTION reads,
        not what the estimator does.

        Shared by `_traverse_solved_history` and the final replay, because a
        replay that opened differently from the solve would report a
        waveform the residual was never driven to zero on.
        """
        tr = self._transient()
        _alphas, b = tr._get_integrator().companion_coefficients(dt, dt)
        tr._begin_run(self._insert_refnode(xm1_in), self.cir.n)
        tr._dt = dt

        ## ⚠ A `b != 0` COMPANION IS REFUSED HERE, AND THE FIRST REASON
        ## GIVEN FOR IT WAS WRONG.  It said a solved history carries CHARGES
        ## while such a method also reads `iq_{-1}`, "which no charge
        ## determines".  The DAE determines it exactly -- a converged point
        ## satisfies `i(x) + iq + u = 0`, so `iq_{-1} = -(i(x_{-1}) + u)`,
        ## the same identity item 4d rests on -- and seeding it was tried.
        ##
        ## It fails for the derivative running the OTHER way.  A one-step
        ## companion reads only `iq_{-1}`, so the trajectory depends on
        ## `x_{-1}` solely through it, and `d(iq_{-1})/d x_{-1} = -G` is
        ## SINGULAR wherever a node carries no conductance -- every purely
        ## reactive node, which is most of a resonator.  Admitting `x_{-1}`
        ## as m unknowns then leaves the 2m x 2m system rank-deficient:
        ## measured, `LinAlgError: Singular matrix` on 25 tests at once.
        ##
        ## The right second unknown for such a method is `iq_{-1}` ITSELF --
        ## the `(x, iq)` state its monodromy already uses -- closed by
        ## `iq_{-1} = iq_{N-1}`.  That is a different formulation, not a
        ## seeding fix, and it is not built.
        if b:
            raise NotImplementedError(
                'a solved entering history admits `x_{-1}` as the second '
                'unknown, and a companion with a b != 0 term depends on it '
                'only through `iq_{-1} = -(i(x_{-1}) + u)`, whose derivative '
                '-G is singular at every purely reactive node -- so the '
                'enlarged system would be rank-deficient. Such a method '
                'needs `iq_{-1}` itself as the unknown, which is a different '
                'formulation; this refuses rather than solving a singular '
                'one.')
        ## The charge half of `_push_history`, without its `_iq` roll: no
        ## step has been solved yet, so there is no companion current to
        ## push -- and a `b = 0` companion never reads one, which the guard
        ## above is what makes true.
        q0 = tr.cir.q(self._insert_refnode(x0_in), tr.epar)
        tr._qlast = self.toolkit.concatenate(
            (self.toolkit.array([q0]), tr._qlast))[:-1]
        tr._q_cache = None
        tr._is_first_step = False
        tr._no_history = False
        ## ⚠ THE STEP THAT PRODUCED `x_0` IS THE PERIOD'S LAST ONE, NOT ITS
        ## FIRST.  `x_{-1}` sits one step BEFORE `x_0`, and on a periodic
        ## grid that step is `hs[-1]`.  With a uniform grid the two are
        ## equal and this never showed; on a caller's grid (item 5) with a
        ## 16438:1 spread, handing `hs[0]` to a method that reads `h_last`
        ## states a step ratio that never happened.
        tr._dt_last = dt if h_prev is None else h_prev
        tr._dt_last2 = None
        self._history_is_solved = True
        return tr

    def _sync_limit_at(self, x_full):
        """Put every limiting device's internal state AT ``x_full``.

        ⚠ EVALUATING "AT A POINT" REQUIRES THE DEVICE LIMITING STATE TO BE AT
        THAT POINT.  A junction device's `i`/`G` are read at its stored `_vlim`,
        not at the vector handed in, and there is only ONE `_vlim` per device --
        while the monodromy evaluates `C`/`G` at SEVERAL distinct points per step
        (`x_n` and every stage).  Whatever the last step's solve happened to
        leave behind (the LAST stage) was therefore used for all of them, so the
        period map linearised the junction at the wrong voltage.

        Measured against a finite-difference derivative of the discrete period
        map (a reference this code cannot influence) on a diode loaded through a
        series resistor: the analytic monodromy was off by a FIXED 1.65e-3
        (Radau) / 9.1e-4 (TR-BDF2) relative -- flat across four decades of the
        FD step, so a genuine error and not FD noise -- and the error grew with
        how hard the junction was driven, vanishing when it was off.  With this
        sync the same comparison lands at ~3e-9, the FD noise floor.

        `limit(x, x)` sets `_vlim` to `x`'s own branch voltage at zero delta, so
        it moves the state without perturbing the point.  Same defect and same
        remedy as the coupled stage solve in `Transient._rk_step_coupled`.
        """
        tr = self._transient()
        tr.cir.limit(x_full, x_full, tr.epar)

    def _C_at(self, x_reduced):
        """The reduced capacitance at a point, without taking a step.

        ⚠ NO LIMITING SYNC HERE, AND THAT IS MEASURED, NOT ASSUMED.  `_G_at`
        needs the device limiting state to be at the point it is evaluating,
        because a junction's `i`/`G` are read at the stored `_vlim`.  CHARGE IS
        NOT: surveyed across every limiter in the tree,

          * `elements.Diode` is the only STATEFUL one (it keeps `_vlim`), and
            its `C`/`q` do not read it -- with the stored state moved far from
            the evaluation point, `dC = dq = 0` while the control `dG = 15.2`
            and `di = 3.9e-1` confirm the limiting was live;
          * `Semiconductor` (BJT/JFET/ZenerDiode/Varactor) limits STATE-FREE by
            construction -- "Return a limited copy of `x` -- STATE-FREE, and
            that is the point";
          * `compact.PspMosLongChannel` likewise returns a limited copy;
          * the hdl devices keep no `_vlim` at all (it is a codegen local).

        So there is no device whose capacitance a sync could correct.  A sync
        was carried here for a while as "correct in principle" insurance and was
        never exercised by any test -- this tree's own rule is that unexercised
        machinery is a liability.  ⚠ If a stateful limiter whose CHARGE reads its
        state is ever added, this is where the sync goes back; `_sync_limit_at`
        is kept for that, and for `_G_at`'s no-junction path.
        """
        tr = self._transient()
        xf = self._insert_refnode(x_reduced)
        C = tr.cir.C(xf, tr.epar)
        (C,) = remove_row_col((C,), self.irefnode, self.toolkit)
        return C

    def _G_at(self, x_reduced):
        """The reduced conductance `di/dx` at a point, without taking a step.

        The companion to `_C_at`.  The one-step-method monodromy needs the
        PHYSICAL `(C, G)` at each stage point -- not the companion
        `Geq = a h G` the accepted step happens to store -- because a
        two-stage DIRK linearises `q` and `i` at THREE distinct points per
        step (`x_n`, the internal stage, and `x_{n+1}`), each with its own
        stage coefficient.  Recovering a physical `G` from a single stored
        `Geq` would divide out only one of those coefficients and mislabel
        the other two.

        ⚠ THIS GOES THROUGH PCNR WHEN THERE ARE JUNCTIONS, and the reason is
        STATELESSNESS, not accuracy.  `pcnr.augmented_system` + `schur_reduce`
        build `G` from an explicitly-passed `v_lim` instead of from the device's
        stored one, so the answer depends on the POINT ALONE.  The `limit(x, x)`
        route `_C_at` still uses does not: `limit` clamps relative to the STORED
        `_vlim`, so it lands on the true point only when the previous evaluation
        was already nearby.  Measured, varying the prior `_vlim` before
        evaluating at a fixed point: PCNR's `J_eff` moves by 0.0, the limit-sync
        `G` by up to 15.15.  It was right in the traversal only BY LOCALITY
        (steps are small, so the prior state is always close) -- the same
        accident `_begin_period` warns about when it insists the period map be a
        function of `x0` alone, applied to its linearisation.

        Numerically this changes NOTHING today: against a finite difference of
        the discrete period map both routes give the same monodromy to every
        printed digit (3.025e-09 radau / 2.358e-09 trbdf2, identical either
        way).  It removes a latent order-dependence, and it is what lets the
        transient and the monodromy share ONE limiting.

        ⚠ `_C_at` CANNOT JOIN: PCNR re-stamps `i`/`G` at `v_lim` but leaves `q`
        alone (`pcnr.py` treats the algebraic equations; diffusion charge is its
        stated caveat), so the capacitance keeps the limit-sync.
        """
        tr = self._transient()
        xf = self._insert_refnode(x_reduced)
        junctions = self._pcnr_junctions()
        if not junctions:
            ## no junction devices: nothing limits, so the plain read IS the
            ## physical G and PCNR would only add an assembly for no reason.
            self._sync_limit_at(xf)
            G = tr.cir.G(xf, tr.epar)
        else:
            from pycircuit.circuit import pcnr as _pcnr
            xfa = np.asarray(xf, dtype=float)
            v_lim = _pcnr.v_lim_init(junctions, xfa)
            g_mna, g_lim, J_mm, _J_ml, _J_lm, didv = _pcnr.augmented_system(
                tr.cir, xfa, v_lim, junctions, tr.epar,
                u_extra=0.0, dense_blocks=False, J_extra=0.0)
            _f_eff, G = _pcnr.schur_reduce(
                g_mna, g_lim, J_mm, junctions=junctions, didv=didv)
            G = np.asarray(G)
        (G,) = remove_row_col((G,), self.irefnode, self.toolkit)
        return G

    def _pq_seed_at_x0(self, x_reduced):
        """``d(iq_0)/d(x_0)`` when the method SEEDS a consistent companion current.

        ⚠⚠ THE CHAIN RULE THE `open_at_x0` PATH ASSUMED AWAY.  Every branch
        that opens at `x_0` seeds `Pq = 0` and says so in the same words --
        "no companion current has been formed yet".  That was true of every
        method in this tree until `theta`, which refuses the L-stable opener
        and therefore READS `iq_{-1}` on its first step: `_begin_run` seeds it
        at ``iq_0 = -(i(x_0) + u(t_0))``, the DAE's own `dq/dt`, and that is a
        FUNCTION OF `x_0`.  Differentiating it gives `-G(x_0)`, and dropping
        that term is not a small error -- it is the whole `null(C)` mode.

        Measured on the B2 gate resonator at `K = 200`, against a
        finite-difference of the shooting residual (delta-swept over six
        decades, FLAT, so a real error and not FD noise): the analytic
        monodromy mapped `null(C)` to ZERO -- exactly what an L-stable Euler
        opener would do -- where the true map multiplies it by `-0.7778`,
        which is `(-(1-theta)/theta)^K` from `ThetaIntegrator`'s own table.
        Relative Jacobian error 6.344; with this seed, 1.4e-10.

        The cost of that was NOT a wrong answer -- the residual is what it is,
        so the solve still lands on the right orbit -- but the Newton lost its
        quadratic step: on a LINEAR circuit an exact shooting Newton converges
        in ONE iteration (`trap` with `x0_unknown` takes 3 evaluations at every
        `K`), and `theta` was taking 9 / 64 / 99 at `K = 100 / 200 / 400`.

        `None` -- the default for every method that does NOT declare
        `needs_consistent_iq0` -- means the zero seed is exact, and those
        methods stay bit-identical.
        """
        integ = self._integrator_for(getattr(self.par, 'method', 'euler'))
        if not integ.needs_consistent_iq0():
            return None
        ## `u(t_0)` carries no `x`, so only `i` contributes: `d(-i)/dx = -G`.
        return -np.asarray(self._G_at(x_reduced), dtype=float)

    def _pcnr_junctions(self):
        """The circuit's PCNR-participating devices, found once and cached.

        `pcnr_devices` walks every element and rebuilds the node map, which is
        far too much to repeat inside `_G_at` -- the monodromy calls it once per
        stage per step per traversal.
        """
        junc = getattr(self, '_pcnr_junctions_cache', None)
        if junc is None:
            from pycircuit.circuit import pcnr as _pcnr
            junc = _pcnr.pcnr_devices(self.cir)
            self._pcnr_junctions_cache = junc
        return junc

    def _traverse(self, x_in, T, times, hs, want_dT, open_at_x0=False):
        """One pass over the period, with the sensitivities accumulated.

        Returns ``(x0, x_end, dx_end/dx0, dx_end/dT)`` -- the last only when
        asked for.  Shared by the fixed-period and autonomous systems so the
        period map is written once; they differ only in what they build from
        it.

        EVERY SHOOTING ITERATION IS ITS OWN RUN.  phi must be a function of
        its arguments alone; if iteration k+1 inherited the ring buffers
        iteration k ended with, the period map would depend on which
        iteration it was and the monodromy would be the derivative of
        something else.  `_begin_run` also makes the first step
        `is_first_step`, so a multi-step method opens at order 1 -- the same
        restart the old `iq_last=None` produced, now for the integrator's
        own reason.
        """
        toolkit = self.toolkit
        n = self.cir.n
        self._want_dfdh = want_dT
        self._begin_period(x_in)
        if open_at_x0:
            ## ⚠ NO MANUFACTURING STEP: the caller's unknown IS `x_0`.
            ## `_begin_period` has seeded both charge rings from `q(x_in)`
            ## and marked the next step `is_first_step`, so the first step
            ## INSIDE the period is order-dropped to Euler -- which is the
            ## L-stable opener this formulation cannot do without.  See
            ## `x0_unknown` on `solve` for why, and for what it costs.
            x = copy(x_in)
            x0 = copy(x_in)
            C_open = np.asarray(self._C_at(x_in))
        else:
            x = self.solve_timestep(x_in, times[0], hs[0])
            x0 = copy(x)
            C_open = np.asarray(self._C)
        ## ⚠ `None`, not `self._iq`, when nothing has been solved yet: no
        ## companion current exists at `x_0` because no step has formed one.
        ## `solve_timestep` reads `iq_last=None` as "open the run", which is
        ## the same restart `_begin_period` already asked for -- passing a
        ## stale `_iq` from a previous traversal would be the hidden-state
        ## defect this class refuses elsewhere.
        iq_last = None if open_at_x0 else self._iq

        ## `Px[k]` is d(x_{n-k})/d(x0), `Pq` is d(iq_n)/d(x0), `Cs[k]` the
        ## capacitance of step n-k.  Two of each -- as far back as any
        ## method here reaches.  Both rings open seeded with the entering
        ## step, mirroring how the transient seeds `_qlast` with `q0`
        ## repeated: at the start of a period there is no earlier point to
        ## differentiate against.
        eye = np.asarray(toolkit.eye(n - 1))
        Px = [eye, eye]
        Cs = [copy(C_open), copy(C_open)]
        if open_at_x0:
            ## ⚠ AND HERE THE SEED IS EXACT RATHER THAN ASSUMED.  With `x_0`
            ## the unknown and both rings holding `q(x_0)`, `dq_{-1}/dx_0`
            ## really IS `C` -- the same `I` the other branch writes down as
            ## an assumption about a history it did not solve for.  `Pq` is
            ## zero because no companion current has been formed yet: the
            ## opening step has not been taken.  That is the whole reason
            ## this path has an exact Jacobian and the other does not.
            ##
            ## ⚠ UNLESS THE METHOD SEEDS ONE.  `theta` refuses the opener and
            ## so reads `iq_{-1}` on step one, where `_begin_run` puts
            ## `-(i(x_0) + u(t_0))` -- a function of the unknown.  See
            ## `_pq_seed_at_x0`; `None` there restores the zero seed exactly.
            Pq = self._pq_seed_at_x0(x_in)
            if Pq is None:
                Pq = np.zeros((n - 1, n - 1))
        else:
            a_first, b_first = self._coeffs
            Pq = (a_first[0] * Cs[0] if b_first else np.zeros((n - 1, n - 1)))
        ## The period column, propagated the same way with one extra source
        ## term.  Zero at the start: the entering state does not depend on T.
        Pt = [np.zeros(n - 1), np.zeros(n - 1)]
        Pqt = np.zeros(n - 1)

        ## Kept for PAC.
        ## ⚠ `C_open`, not `self._C`: on the `open_at_x0` path no step has
        ## run, so `_C` does not exist yet.  `C_open` is the same matrix the
        ## other branch would have found there, evaluated at `x_0` directly.
        ## (Kept for PAC, which is withdrawn; `Jtvec` opens empty for the
        ## same reason -- there is no solved `Jf` at `x_0`.)
        self.Cvec = [copy(C_open)]
        self.Jtvec = ([] if open_at_x0 else [copy(self._Jf)])
        self.times = times

        for _j, t in enumerate(times[1:]):
            dt = hs[min(_j, len(hs) - 1)]
            x = copy(self.solve_timestep(x, t, dt, iq_last=iq_last))
            iq_last = self._iq
            self.Cvec.append(copy(self._C))
            self.Jtvec.append(copy(self._Jf))

            ## ONE RECURSION FOR EVERY METHOD.  Each writes its companion as
            ## `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`, so differentiating
            ## the step gives
            ##
            ##     S    = sum_{k>=1} a_k C_{n-k} P_{n-k} + b Pq
            ##     P_n  = -Jf_n^-1 S
            ##     Pq_n = a_0 C_n P_n + S
            ##
            ## Euler is `b = 0` reaching back one step, trapezoidal `b = -1`
            ## reaching back one, Gear-2 `b = 0` reaching back two.  The
            ## coefficients come from the integrator that RAN, so an
            ## order-dropped step contributes its own.
            ##
            ## ⚠ A SOLVE, NOT AN INVERSE (stage 11).  `inv(Jf) @ ...` formed
            ## a dense inverse per timestep per iteration and squared the
            ## condition number it then multiplied through.
            alphas, b = self._coeffs
            Jf = np.asarray(self._Jf)
            C_new = np.asarray(self._C)

            ## ⚠ THIS USED TO BE AN INLINE COPY of `_step_sensitivity`,
            ## byte-for-byte, while that method's own docstring said the two
            ## systems "share this and not a copy".  They did not: only the
            ## solved-history path called it.  Found while measuring the
            ## plain path's propagation share -- a timer wrapped around
            ## `_step_sensitivity` reported 0.0%, which is not a small
            ## number but a wrong one.
            Px_new, Pq = self._step_sensitivity(Px, Cs, Pq, Jf, C_new)

            if want_dT:
                ## The SAME recursion, plus the step size's own dependence
                ## on T.  Every step scales together (`h = T/(N-1)`), so
                ## `dh/dT = h/T`, and `df/dh` at fixed solution is Fang's
                ## `p` -- already shared.  For an autonomous circuit its
                ## `du/dt` half vanishes, which is what makes this term the
                ## companion derivative alone.
                St = b * Pqt if b else np.zeros_like(Pt[0])
                for k in range(1, len(alphas)):
                    St = St + alphas[k] * (Cs[k - 1] @ Pt[k - 1])
                if self._period_column == 'closing':
                    ## only the CLOSING step's length depends on `T`, and it
                    ## does so with `dh/dT = 1` -- so the term appears once,
                    ## undivided, on the last step and nowhere else.
                    if _j == len(times) - 2:
                        St = St + np.asarray(self._dfdh).ravel()
                else:
                    St = St + np.asarray(self._dfdT).ravel() / T
                Pt_new = -self.toolkit.linearsolver(Jf, St)
                Pqt = alphas[0] * (C_new @ Pt_new) + St
                Pt = [Pt_new, Pt[0]]

            Px = [Px_new, Px[0]]
            Cs = [copy(C_new), Cs[0]]

        self._want_dfdh = False
        ## Kept for the autonomous check after the solve: its spectrum is
        ## the only place a free period announces itself.
        self._monodromy = Px[0]
        return x0, x, Px[0], (Pt[0] if want_dT else None)

    def _transient(self):
        """The `Transient` this analysis integrates with.

        PSS used to carry its OWN transcription of one integrator step --
        the third in the tree, after `Transient` and `JAXTransient` -- and
        it had already cost the two defects its docstring records: `method`
        declared and never read, and a companion current fed back from the
        iterate before the converged one.  Driving the real thing removes
        the copy and brings what came with it: the limiting machinery, PCNR,
        breakpoint order drops, and the continuation rescue.

        ⚠ THAT LIST WAS 1-FOR-4 AS SHIPPED (external review, 2026-09-02); it is
        now 2-FOR-4.  LIMITING reaches -- `cir.limit` is called on the inner
        Newton and the rectifier measurably conducts.  PCNR now reaches too:
        `PSS(cir, pcnr=True)` is a declared Parameter forwarded to the inner
        `Transient` above, and PCNR lives in `Transient.solve_timestep` (the
        LMM `_solve_timestep_pcnr` and, for stage methods, `_rk_stage_pcnr`),
        which PSS DOES call -- so it needs no per-accepted-step machinery.  The
        remaining two still do not: the continuation rescue (`_rescue_solver`)
        and breakpoints (`cir.next_event`) are armed only inside
        `Transient.solve`, which PSS never calls -- it drives `solve_timestep`
        directly on its own frozen grid, so a breakpoint has nothing to move.
        The same structural fact behind the TLine refusal above: what
        `Transient.solve` does per accepted step, PSS does not do at all.

        The tolerances are handed over unchanged, which is the point of
        `newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
        same thing on both sides, so passing them through is a no-op in
        meaning.
        """
        if getattr(self, '_tran', None) is None:
            ## ⚠ A MAPPING, not an if/else on 'euler'.  Written as
            ## `EulerIntegrator() if method == 'euler' else Trapezoidal...`
            ## it silently ran trapezoidal for every other name -- caught
            ## while adding 'gear', which produced numbers identical to
            ## trap's to the last digit.  This class has already paid once
            ## for a `method` that selected nothing; a dict raises KeyError
            ## on a name nobody wired.
            self._tran = self._new_transient(
                self._integrator_for(self.par.method))
        return self._tran

    def _theta_biased(self, integ):
        """Give a `ThetaIntegrator` the bias THIS period needs, not a fixture's.

        ⚠⚠ `theta - 1/2 = C h` makes `C` a RATE, but the quantity that decides
        anything is the DIMENSIONLESS product `C T`: `null(C)` is damped over a
        period by `((1-theta)/theta)^K ~ exp(-4 C h K) = exp(-4 C T)`, and `h`
        cancels.  So `ThetaIntegrator.DEFAULT_C = 1e4` is not a recipe -- it is
        the measured knee `C T = 0.0628` divided by ONE fixture's period
        (6.283e-6 s).  On a circuit 159x slower it is 159x the calibrated bias,
        and MEASURED on `_q20_rlc` (analytic 20 V) that is a peak of 15.91 at
        K = 100 -- 20% low, `converged=True`, because it did converge: to its
        own over-damped discretisation.  See `ThetaIntegrator.DEFAULT_CT`.

        This is where a shooting run stops inheriting that.  `theta_ct` is the
        dimensionless knob (`None` = the measured knee) and the period is this
        solve's, so `method='theta'` is now correct on any circuit.

        ⚠ THE SEED PERIOD IS ENOUGH, and that is a measurement not a hope: the
        gate's own table has `rcond(I - A^K)` at 4.4e-03 / 4.6e-03 / 4.3e-03
        across `C T` = 0.0063 / 0.0628 / 0.628, i.e. FLAT over two decades.  An
        autonomous solve moving `T` by a few percent moves the bias by the same
        few percent, which the knee does not notice.

        Every other method is returned untouched, so nothing else moves a bit.
        """
        from pycircuit.circuit.integrator import ThetaIntegrator
        T = getattr(self, '_theta_period', None)
        if not isinstance(integ, ThetaIntegrator) or T is None:
            return integ
        ct = getattr(self.par, 'theta_ct', None)
        integ.cbias = float(ThetaIntegrator.DEFAULT_CT if ct is None
                            else ct) / float(T)
        return integ

    def _new_transient(self, integ):
        """A `Transient` on this circuit driven by `integ`, with every
        strategy and tolerance PSS was given handed through.

        Extracted so the TR-BDF2 monodromy can build a transient on a
        DIFFERENT integrator than `self.par.method` without duplicating the
        pass-through -- and so the pass-through cannot drift between the two
        call sites.

        ⚠ THE SOLVER STRATEGIES GO THROUGH TOO, and they used not to.
        `nrsolver`, `linearsolver` and `scaler` are declared on the base
        `Analysis`, so `PSS(cir, linearsolver=...)` has always been ACCEPTED
        -- and then dropped here, with the inner `Transient` resolving to
        `DenseSolver`/`StandardNewton` whatever the caller asked for.  That
        was the third time this class took a parameter it never read
        (`method` declared and never read; `analysis='PSS'` matching
        nothing), and the same shape each time: accepted at the constructor,
        silently discarded at the boundary.
        """
        from pycircuit.circuit.transient import Transient
        ## ⚠ ONE CHOKE POINT for the theta bias, because there are four call
        ## sites building a transient and a per-site fix would drift.  A no-op
        ## for every other integrator -- see `_theta_biased`.
        integ = self._theta_biased(integ)
        tr = Transient(
            self.cir, toolkit=self.toolkit, integrator=integ,
            reltol=self.par.reltol, iabstol=self.par.iabstol,
            vabstol=self.par.vabstol, maxiter=self.par.maxiter,
            analysis=self.par.analysis,
            lte_vabstol=self.par.lte_vabstol,
            lte_iabstol=self.par.lte_iabstol,
            TRTOL=self.par.TRTOL, relref=self.par.relref,
            nrsolver=self.par.nrsolver,
            linearsolver=self.par.linearsolver,
            scaler=self.par.scaler,
            pcnr=self.par.pcnr)
        ## The line search as the last resort on the shooting path, which
        ## never arms the transient's rescue ladder (owner decision
        ## 2026-09-08, "Do 2"; see `_rk_step_coupled` and `solve_timestep`).
        tr._damped_last_resort = True
        tr.irefnode = self.irefnode
        return tr

    def _begin_period(self, x_reduced):
        """Start one traversal of the period from a clean integrator state.

        Every shooting iteration re-integrates the SAME interval from its own
        `x0`, so "begin a run" happens once per iteration here, not once per
        analysis.  Without the reset, iteration k+1 would inherit the ring
        buffers iteration k ended with and the period map would depend on
        which iteration it was -- phi must be a function of `x0` alone or the
        monodromy is the derivative of something else.
        """
        tr = self._transient()
        tr._begin_run(self._insert_refnode(x_reduced), self.cir.n)
        return tr

    def _insert_refnode(self, x):
        return self.toolkit.concatenate(
            (x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

    def solve_timestep(self, x0, t, dt, refnode=gnd, iq_last=None):
        """One timestep of the inner transient, taken by `Transient`.

        This used to be a private transcription of one integrator step --
        the third in the tree -- and it had already cost two defects that
        its own comments recorded: `method` was declared and never read, so
        PSS was backward-Euler only, and the companion current fed back to
        the next step belonged to the iterate BEFORE the converged one.
        Both are structurally impossible now: the integrator is an
        `Integrator` object driven by `Transient.get_diff`, and the
        companion current is the one that class stores at its own converged
        point.

        What came with the change, none of which the copy had: the limiting
        machinery (measured -- a rectifier whose diode never turned on, so
        the non-conducting solution was returned as a converged periodic
        steady state), PCNR when the circuit and Parameters ask for it,
        breakpoint order drops, and the continuation rescue.

        `dt` is imposed by the caller: PSS owns the grid, which is what
        keeps the period map a smooth function of `x0`.  `iq_last` is
        retained in the signature for callers that pass it, but the
        companion history now lives in the `Transient`'s own ring buffers,
        rolled here through `_push_history`.

        The measured cost of backward Euler on a limit cycle is unchanged
        and still the reason `method` matters -- it damps exactly what PSS
        exists to find:

            steps/period    PSS peak    fraction of analytic
                      20      2.63 V       13.2%
                      50      5.61 V       28.1%
                     100      8.81 V       44.1%
                     200     12.20 V       61.0%
        """
        toolkit = self.toolkit
        irefnode = self.irefnode
        tr = self._transient()

        ## ONE INTEGRATOR STEP, taken by the class that owns the definition.
        ## `Transient.solve_timestep` applies the chosen integrator through
        ## `get_diff` (so `method` selects something because the integrator
        ## object does), the limiting machinery, PCNR when asked for, and the
        ## continuation rescue.  None of that existed on the copy this
        ## replaced.
        tr._dt = dt
        x_full, J_full = self._transient_step(tr, x0, t)

        ## d(residual)/dh at the converged point, for the period
        ## derivative -- BEFORE `_push_history`, because it reads `_qlast`
        ## as the PREVIOUS charges, which the push is about to overwrite.
        ## `Transient.residual_dh` is Fang's `p`, already shared: it is
        ## `d(iq)/dh + du/dt`, and for an AUTONOMOUS circuit the second term
        ## is identically zero -- which is exactly why solving for the
        ## period is tractable here and would not be on a driven circuit,
        ## where scaling T also moves every source evaluation.
        if self._want_dfdh:
            ## ⚠ `residual_dT`, NOT `residual_dh`.  The grid is rebuilt at
            ## the current `T` on every residual evaluation, so every step
            ## scales and the partial `d/dh_n` is not the total -- for
            ## Gear-2 it is 3/2 of it on a uniform grid, measured against
            ## finite differences at 1.4859/1.4939/1.4972 for 100/200/400
            ## points, converging on the exact 3/2.  Euler and trapezoidal
            ## were never wrong: their coefficients depend on `h_n` alone,
            ## so the partial IS the total, which is why only Gear-2 was
            ## hit.  See `Integrator.companion_dT`.
            if self._period_column == 'closing':
                ## ⚠ `residual_dh`, THE PARTIAL, NOT `residual_dT`.  The
                ## note above explains why the total is 3/2 of the partial
                ## for Gear-2: `residual_dT` accounts for EVERY step scaling
                ## with `T`.  Under the closing-step convention only ONE
                ## step's `h` moves, so the partial IS the derivative and
                ## the 3/2 would be exactly the error.
                (self._dfdh,) = remove_row_col(
                    (tr.residual_dh(x_full, t, dt),), irefnode, toolkit)
            else:
                (self._dfdT,) = remove_row_col(
                    (tr.residual_dT(x_full, dt),), irefnode, toolkit)
        ## Measured, not controlled: the grid is the caller's, so nothing can
        ## act on this.  Also before the push, for the same reason.
        if self._want_lte:
            self._lte = tr.step_lte(x_full, self._insert_refnode(x0), J_full)
            ## A SEAM STEP IS ONE WHOSE COMPANION READS THE ENTERING
            ## UNKNOWN, not merely one whose ESTIMATOR does.
            ##
            ## ⚠ THIS CONDITION WAS `h_last2 is None` ALONE, AND THAT FLAGGED
            ## A PHANTOM FOR TWO METHODS OF THREE.  `h_last2 is None` is the
            ## transient's statement that the third past charge is not real,
            ## which is the reach of the LTE estimator's third divided
            ## difference -- not the reach of the integrator.  Euler's
            ## companion reads `q_{n-1}`; trapezoidal's reads `q_{n-1}` and
            ## `iq_{n-1}`, and the order-dropped opening step supplies an
            ## `iq` consistent with it, which is what that drop is FOR.
            ## Neither can see the fabricated charge at all.  Measured
            ## (`benchmarks/pss_seam_cost.py`): trapezoidal's seam reading
            ## was 15.1 times tolerance while its cost is 1.3e-11 V, and
            ## euler's 0.286 against 5.1e-12 V.  Both are exactly zero; the
            ## reading was an artefact of the measurement.
            ##
            ## Gear-2 reads `q_{n-2}` -- which at that step IS the entering
            ## unknown -- and the shooting condition constrains `x(0)` to
            ## equal `x(P)`, NOT `x_in` to be the orbit's own `x(-dt)`.  So
            ## `x_in` sits O(h^2) off a real history point and Gear-2 reads
            ## it as one.  That one costs 1.266e-01 V at 100 points/period
            ## against an interior contribution of 1.070e-01 -- the seam is
            ## 54% of its total error -- and it is the term that STOPS
            ## converging: it falls as h^2 while the interior falls faster,
            ## so its share grows to 68% at 200 points and 73% at 400.
            ##
            ## TWO DIFFERENT THINGS ARE TRUE OF AN OPENING STEP, and the
            ## first version of this conflated them -- which showed up as
            ## trapezoidal's phantom simply MOVING from the seam into the
            ## interior total (0.340 -> 15.47) when the seam test was
            ## tightened.  Suppressing a bad number is not the same as
            ## classifying it.
            ##
            ##   the ESTIMATE is invalid when the ESTIMATOR differences a
            ##   charge that was never a real point.  Its divided difference
            ##   reaches `p = ORDER + 1` charges back -- 2 for Euler, 3 for
            ##   the second-order pair -- so at the step with only two real
            ##   past charges, euler's estimate is sound and trap's and
            ##   gear2's are not.  An unsound estimate is DISCARDED: it is
            ##   not an interior reading and, on its own, not a seam either.
            ##
            ##   a SEAM exists when the COMPANION reads the entering
            ##   unknown, i.e. when its charge reach `len(alphas) - 1` is
            ##   deep enough to touch it.  Only Gear-2's is.
            ##
            ## `_dt_last2 is None` says exactly "two real past charges" here
            ## (it is set from `_dt_last` one step later), so it is the step
            ## index in disguise; both tests are written against the count.
            _real_past = 2 if tr._dt_last2 is None else 3
            _p = getattr(tr.active_integrator, 'ORDER', 1) + 1
            _reach = len(tr._companion_coeffs[0]) - 1
            self._lte_valid = _real_past >= _p
            ## `_history_is_solved` is that formulation saying the
            ## deepest charge is an UNKNOWN the solve closed, not a stand-in
            ## -- so there is no seam to report even though the companion
            ## reaches that far.  Without this the fix would go on flagging
            ## the defect it removed.
            self._lte_seam = (_reach >= _real_past
                              and not self._history_is_solved)

        ## The history advance is the accept path's, called rather than
        ## copied -- and `_dt_last` must roll AFTER the step, because
        ## `get_diff` read it as `h_last` while solving.
        tr._push_history(x_full)
        tr._dt_last2 = tr._dt_last
        tr._dt_last = dt
        tr._is_first_step = False
        tr._no_history = False

        ## Reduced-system views for the shooting Jacobian.  `_Geq` is the
        ## companion conductance the step actually used, which is the factor
        ## the monodromy needs; `_iq` is kept for the caller's own bookkeeping
        ## as before.
        (self._Jf, self._Geq, self._C) = remove_row_col(
            (J_full, tr._Geq, tr._Cmat), irefnode, toolkit)
        ## The coefficients of the integrator that ACTUALLY ran this step --
        ## an order drop on the opening step reports Euler's, which is what
        ## the propagation must use for that step.
        self._coeffs = tr._companion_coeffs
        self._iq = tr._iq

        x = toolkit.concatenate((x_full[:irefnode], x_full[irefnode + 1:]))
        return x

    def _transient_step(self, tr, x0_reduced, t):
        """`Transient.solve_timestep` on the FULL vector, returning
        ``(x, J)``.  PSS works on the reduced system throughout; this is the
        one place the two conventions meet."""
        x_full = self._insert_refnode(x0_reduced)
        x, _feval, J, _f = tr.solve_timestep(x_full, t)
        return x, J


    def find_initial_solution(self, period, x0=None, npts=60, method=None,
                              eps_rel_lin=1e-2, eps_abs_lin=1e-3, n_iter=7,
                              max_periods=200, zeta=1e-7):
        """A PROPER initial solution to start shooting from, by pre-integrating
        until the fixed-point iteration has entered its LINEAR region.

        De Luca, Bolcato & Schilders, *Proper Initial Solution to Start Periodic
        Steady-State-Based Methods*, IEEE TCAS-I 2019 -- their Algorithm 2.  The
        paper is in `~/docs/07-shooting-methods/`.

        Shooting-Newton needs a start inside its contraction region, and the
        usual remedy is to GUESS a number of pre-integration periods; if the
        guess is wrong the solve diverges and the guess is repeated with no clue
        for the next one.  This detects the handoff point instead, from
        quantities the integration already produces.

        ⚠ WHAT THE CRITERION COMPARES -- AND WHAT IT IS NOT.  It is NOT "the
        iterate stopped moving" and NOT "a carried probe settled": this project
        measured that guess and it moves the WRONG WAY (drift 1.4e-2 while the
        solve still fails, 1.2e-1 once it succeeds -- `benchmarks/
        pss_warm_start.py`), because a settled probe only says the Jacobian
        stopped changing, which is equally true at an equilibrium.  The paper
        compares TWO SEQUENCES: the LINEAR prediction of the next shooting error
        against the one the ACTUAL nonlinear integration produces.  They agree
        only where the fixed-point map really has become linear, which is
        exactly the region a Newton-type method needs::

            u_k         = x_k - phi(x_k)                        (eq. 4)
            u_{k+1}     = J_phi(x_khat) u_k                     (eq. 12)
            utilde_{k+1} = x_{k+1} - phi(x_{k+1})               (eq. 13)

        accepted, componentwise, when (eq. 16)::

            |u_{k+1,j} - utilde_{k+1,j}| <= eps_rel_lin |u_{khat,j}| + eps_abs_lin

        holds for ``n_iter`` CONSECUTIVE iterations (the paper's defaults, used
        for every experiment in it: ``eps_rel_lin=1e-2``, ``eps_abs_lin=1e-3``,
        ``n_iter=7``).  The scale on the right is the shooting error at the
        DETECTION index ``khat``, not at the current one.  A failed check resets
        the run AND re-freezes ``J_phi`` at the new index, which is why ``khat``
        can move.

        ⚠ NON-AUTONOMOUS ONLY, AND THAT IS THE PAPER'S SCOPE, NOT AN OVERSIGHT
        HERE.  Its title, abstract and index terms all say non-autonomous, and
        the reason bites: for a FORCED circuit the DC point is not a fixed point
        of the period map, so "the map became linear" can only mean the orbit.
        For an AUTONOMOUS oscillator the equilibrium IS a fixed point of the
        period map and the map is linear in a neighbourhood of it, so this
        criterion will happily certify the trivial root -- which is precisely the
        van der Pol failure recorded in `benchmarks/pss_warm_start.py`.  **That
        case is NOT solved by this method and must not be handed to it.**

        ``J_phi(x_khat) u`` is taken by the paper's own alternative, the
        directional derivative of eq. (15),
        ``[phi(x + zeta u) - phi(x)] / zeta``, rather than by its Alg. 1
        left-product.  Both are in the paper; this one costs one extra period
        integration per iteration and buys freedom from the opening-frame
        question (which state a stored factorisation is the derivative *about*),
        a distinction that has already cost this file one wrong answer.

        Returns ``(x, info)``: the reduced state to start shooting from -- the
        paper's line 21, "the last computed x_k" -- and a dict with ``khat``,
        ``periods``, ``found`` and the per-period ``history``.  ``found=False``
        means ``max_periods`` ran out with no linear region; the returned ``x``
        is then simply the last iterate and carries no promise.
        """
        T = float(period)
        n = self.cir.n
        iref = self.irefnode
        m = n - 1
        npts = int(npts)
        if npts < 1:
            raise ValueError('find_initial_solution: npts must be >= 1, got %r'
                             % (npts,))
        if n_iter < 1:
            raise ValueError('find_initial_solution: n_iter must be >= 1')
        times = np.linspace(0.0, T, npts + 1)
        hs = np.diff(times)
        integ = method if method is not None else getattr(self.par, 'method',
                                                          'euler')

        def phi(xr):
            """One period of the inner transient from `xr` -- the map the
            shooting residual is built on, run on a fixed uniform grid."""
            tr_saved = getattr(self, '_tran', None)
            self._tran = self._new_transient(self._integrator_for(integ))
            try:
                self._want_dfdh = False
                self._want_lte = False
                self._begin_period(np.asarray(xr, dtype=float))
                x = copy(np.asarray(xr, dtype=float))
                for j, t in enumerate(times[1:]):
                    x = copy(self.solve_timestep(x, t, hs[j]))
                return np.asarray(x, dtype=float).ravel()
            finally:
                self._tran = tr_saved

        if x0 is None:
            x = np.zeros(m)
        else:
            x = np.asarray(x0, dtype=float).ravel()
            if x.shape[0] == n:
                x = np.concatenate((x[:iref], x[iref + 1:]))

        phi_x = phi(x)
        k = 0
        i_iter = 0
        khat = 0
        x_khat = x.copy()
        phi_khat = phi_x.copy()
        u = x - phi_x
        u_khat = u.copy()
        history = []
        found = False
        while k < int(max_periods):
            if i_iter == 0:
                ## re-freeze the linear generator at the current index (Alg. 2
                ## lines 4-8): khat moves whenever the run of successes breaks.
                khat = k
                x_khat = x.copy()
                phi_khat = phi_x.copy()
                u = x - phi_x
                u_khat = u.copy()
            x_next = phi_x
            phi_next = phi(x_next)
            u_tilde_next = x_next - phi_next            ## eq. (13)
            ## eq. (12) via the paper's eq. (15) directional derivative, with
            ## the step scaled to the iterate so `zeta` is a RELATIVE size.
            un = float(np.linalg.norm(u))
            if un <= 0.0:
                u_next = np.zeros(m)
            else:
                z = zeta * max(float(np.linalg.norm(x_khat)), 1.0) / un
                u_next = (phi(x_khat + z * u) - phi_khat) / z
            gap = np.abs(u_next - u_tilde_next)
            ok = bool(np.all(gap <= eps_rel_lin * np.abs(u_khat) + eps_abs_lin))
            k += 1
            i_iter = i_iter + 1 if ok else 0
            history.append({'k': k, 'khat': khat, 'ok': ok, 'run': i_iter,
                            'gap': float(np.max(gap)),
                            'shooting_error': float(np.max(np.abs(u_tilde_next)))})
            x, phi_x, u = x_next, phi_next, u_next
            if i_iter >= int(n_iter):
                found = True
                break
        return x, {'khat': khat, 'periods': k, 'found': found,
                   'history': history}

    ## Set by `solve`; declared here so a caller may read them on a PSS that
    ## has not solved yet, and so `event_times` is never a stale leftover.
    break_events = False
    event_times = []

    ## No fold until a solve collects one -- `_fold_periodic` is then the
    ## identity, which is exactly right for every circuit without a folding
    ## state and for any caller that reaches a residual outside `solve`.
    _periodic_fold = []

    def _collect_periodic_fold(self):
        """`[(reduced_row, modulus)]` for every state defined up to `n*modulus`.

        `Circuit.periodic_states()` reports GLOBAL rows; the shooting residual
        lives in reduced coordinates, so each is mapped through the same
        `irefnode` deletion the stamps use.  A declared row that IS the
        reference node has no reduced coordinate and is dropped.

        The offset is deliberately discarded: a DIFFERENCE of two states on a
        periodic row is defined up to `n*modulus` whatever window each state
        was folded into, so only the modulus enters.
        """
        cir = self.cir
        if not hasattr(cir, 'periodic_states'):
            return []
        try:
            declared = cir.periodic_states()
        except Exception:
            ## A late-bound modulus that cannot be resolved is not a reason to
            ## refuse the solve -- it is a reason not to fold.
            return []
        iref = self.irefnode
        out = []
        for row, modulus, _offset in declared or []:
            row = int(row)
            m = float(modulus)
            if row == iref or not np.isfinite(m) or m <= 0.0:
                continue
            out.append((row if row < iref else row - 1, m))
        return out

    def _fold_periodic(self, F):
        """Fold a shooting residual's periodic rows into `[-m/2, m/2)`.

        ⚠ THIS IS WHAT MAKES THE PERIOD MAP'S FIXED-POINT PROBLEM WELL POSED
        FOR A FOLDING STATE, and it belongs in the RESIDUAL rather than in the
        traversal.  `Idtmod`'s state is defined only up to `n*modulus` -- it
        says so itself through `periodic_states()`, and the transient engine
        already uses that declaration to keep the state bounded by exact gauge
        translations.  Shooting did not: it asked for `x_0 - phi(x_0) == 0`
        literally, which on a folding row demands the SAME REPRESENTATIVE, not
        the same state.  An orbit that closes after advancing exactly one
        modulus -- the normal case for a phase -- then has NO root at all, and
        near the fold the raw difference jumps by a whole modulus while the
        state moves infinitesimally.

        Measured (see `test_a_state_fold_breaks_the_period_map_at_the_ENDPOINT_
        not_on_the_grid`): that jump is grid-INDEPENDENT -- 1.414214e+09 at
        seven different grids, with the wrap on a node and off it alike -- so
        it is not an event-localisation defect and no refinement of the time
        grid can reach it.  It is the output map, and the output map is what
        this folds.

        ⚠ THE JACOBIAN IS DELIBERATELY NOT TOUCHED.  `d/dx0` of
        `wrap(x0 - phi(x0))` equals `d/dx0 (x0 - phi(x0))` almost everywhere --
        the wrap has unit slope between its jumps -- so `D - alpha*Mx` is
        already the right derivative of the folded residual.  The fold moves
        the residual onto the branch the Jacobian was always describing; that
        is the whole reason this is a local change and not surgery on the six
        traversal loops.
        """
        rows = self._periodic_fold
        if not rows:
            return F
        F = np.asarray(F, dtype=float).copy()
        for r, m in rows:
            if r < F.shape[0]:
                F[r] -= m * np.round(F[r] / m)
        return F

    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20, grid=None, matrix_free=False,
              x0_unknown=None, tstab=None, break_events=None):
        """Solve for the periodic steady state.

        `break_events` lands the circuit's source discontinuities on grid
        points (`event_grid`).  `None` -- the default -- decides from the
        METHOD: on for a one-step method, off for a multistep one, because it
        HELPS the first and HURTS the second.  See `_resolve_break_events` for
        the measurement and the control that pins the cause.  A circuit whose
        sources declare no discontinuity is untouched bit-for-bit either way.

        `grid` is RECORDED SCOPE ITEM 5: a sequence of step FRACTIONS of the
        period, summing to 1, used in place of the uniform `timestep` grid.
        Fractions rather than absolute times because an autonomous period is
        an unknown and every step has to scale with it.  See
        `_period_grid`; the grid is still frozen for the whole solve, so the
        shooting Newton stays exact.

        `tstab` runs a TRANSIENT for that many seconds before shooting and
        uses its final state as the seed -- the stabilisation time every
        commercial PSS offers, and the standard answer to a seed that is not
        close enough.  `None` (the default) shoots from `x0`, or from the
        operating point, exactly as before.

        ⚠ IT IS THE REMEDY FOR THE FAILURE THIS ANALYSIS FAILS MOST OFTEN.
        Seeded near the unstable DC point -- the trivial-root basin, which is
        where an unseeded autonomous run starts -- van der Pol does not
        solve at all, and one period of `tstab` fixes it:

              circuit                    without tstab      periods needed
              mu = 1  (strongly attracting)  LinAlgError            1
              mu = 0.05 (high-Q)             not converged         ~24

        The count is the `1/mu` amplitude-envelope constant, so it is a
        property of how strongly the limit cycle attracts and NOT of how bad
        the seed is: from 4x and even 20x the orbit amplitude, one period
        suffices at `mu = 1`.  A high-Q oscillator needs proportionally
        more, which is the usual guidance stated as a number.

        ⚠ AND THE STOPPING POINT IS THE CALLER'S, DELIBERATELY.  De Luca,
        Bolcato & Schilders (2019) give a criterion for detecting the handoff
        automatically, and the probe it rests on was measured here and does
        NOT identify it -- near the DC point the monodromy is nearly
        constant, so the probe settles into its own eigenvector and reports
        convergence while the state is stuck at the trivial root.  Every
        obvious substitute shares that defect, because the trivial root IS a
        fixed point of the period map and passes every periodicity test.  So
        the number is asked for rather than guessed.

        ⚠ AND THE COMMERCIAL AUTOMATIC CRITERION IS REPORTED NOT TO SURVIVE
        Q EITHER.  A commercial simulator does offer one; this tree's user reports from
        their own practice that it does not work properly on circuits of
        even moderate Q and does not work on high-Q ones.  That is field
        experience and not a measurement, and it is recorded as such --
        but it points the same way as the probe measured here, which
        reported "settled" while the state was still in the trivial-root
        basin, and it is the same axis: the harder the oscillator attracts,
        the longer the approach and the easier it is for a detector to stop
        early.  Measured above, `mu = 0.05` needs ~24 periods and FIVE is
        not enough.

        ⚠ AND IT CANNOT ESCAPE AN EQUILIBRIUM IT IS STARTED ON.  With
        `x0=None` the seed is the operating point, which on an autonomous
        circuit is an equilibrium -- a transient started exactly there never
        leaves, so no amount of `tstab` helps.  The pre-integration needs
        somewhere to go: pass an `x0` off the equilibrium, or an `ic` on a
        device.  That limit is why this is not a substitute for the probe
        technique, which pumps energy in precisely so the solve cannot fall
        to the DC point.

        ⚠ THE SAME STATEMENT EXISTS IN THE LITERATURE, in words, and the
        measurement above is its quantitative form.  Kundert
        (*Introduction to RF Simulation*, v2 2003; relayed, cited not
        verified here) on starting shooting from a plain transient: "this
        is sufficient to get convergence even on troublesome circuits
        EXCEPT WHEN THE TIME CONSTANTS IN THE CIRCUIT ARE MUCH LARGER THAN
        THE PERIOD OF THE SIGNAL."  That is exactly the mu = 1 -> 1 period
        against mu = 0.05 -> ~24 periods above, and the reason the count
        tracks the `1/mu` envelope rather than the seed: the envelope time
        constant IS the time constant he names.

        See `benchmarks/pss_warm_start.py` for the probe that failed and the
        counts that decided the interface.

        ⚠ `x0_unknown` DEFAULTS TO `None`, WHICH MEANS "DECIDE FROM THE
        TOPOLOGY": it is switched ON automatically for a netlist the index
        criterion proves is index 2, where the manufactured opening step is
        inconsistent, and left OFF otherwise -- because it is not free (see
        the measured regression below). An explicit `True` or `False` is
        honoured untouched. See `_resolve_x0_unknown`.

        `x0_unknown` solves for `x_0` itself instead of for `x_in`, the
        pre-image of a manufactured opening step.  The plain path's default
        is to manufacture `x(0)` with one order-dropped Euler step and hand
        `fsolve` that step's PRE-IMAGE, while the Jacobian it returns is
        taken with respect to `x_0` -- a frame error, and the reason the
        true `dF/dx_in` is SINGULAR (rank 1/3, 2/4, 1/3 on three circuits,
        `sigma_min` exactly 0) and the iteration is a CONTRACTION with a
        linear rate rather than a Newton.  With this set there is no
        manufacturing step, the unknown is the period's own start, and the
        Jacobian is exact.

        ⚠ IT IS NOT FREE, AND THE TRADE-OFF IS GRID-DEPENDENT.  Trapezoidal
        still needs an L-stable opener -- without one its period map is
        `A_trap^K`, singular at EVEN K on every MNA circuit (the
        `(-1)^n` obstruction, now a theorem; see item 4d) -- so the Euler
        step moves INSIDE the period, where it degrades the ORBIT and not
        just the opening.  Measured on a Q=20 resonator against its analytic
        20 V peak:

              npts    default (x_in)    x0_unknown
               100       20.01273        19.76939
               200       20.02208        19.96123

        18x worse at 100 points, 1.8x at 200; the order is preserved and the
        gap closes as the single first-order step is diluted, but the
        constant is real.

        ⚠ AND ON A HARD CASE IT IS THE OTHER WAY ROUND, which is why this is
        an option and not a repair.  Van der Pol on its 1105-step LTE grid
        reaches -47.3 ppm in four iterations against the default path's
        -73.8 -- exactly the figure the throwaway driver in
        `benchmarks/pss_lte_grid.py` has been recorded as "the target to
        hit" since item 5 was written, and for the same reason: its unknown
        was `x_0` too.  Convergence is quadratic at even and odd point
        counts alike (2.9e+00 -> 9.1e-08 -> 1.1e-14).

        So: uniform grid on a benign circuit, leave it off; non-uniform grid
        on a stiff one, turn it on, where the quadratic rate is what decides
        whether it converges at all.  `benchmarks/pss_x0_unknown.py` is the
        driver those numbers come from.  Euler is unchanged to the digit
        either way -- its manufacturing step IS an Euler step, so the two
        maps coincide, which is what makes it the control.

        `matrix_free` is RECORDED SCOPE ITEM 6: solve the outer system
        without ever forming the monodromy, propagating ONE vector per
        Krylov iteration instead of `2m` columns per step.  It is worth
        asking for on LARGE circuits only -- measured on the RC ladder it
        LOSES at m=40 (0.82x) and wins above roughly m=250: 1.40x at m=242,
        2.23x at m=502, 3.62x at m=1002.  ⚠ It buys that with memory, `2 N
        m^2` doubles of stored factorisations (~800 MB at m=1002, 50
        points), and it does NOT produce a monodromy, so `spectral_radius`
        is `None` after a matrix-free solve -- not forming that matrix is
        the entire point.

        ⚠ EVERY FIGURE BELOW HAS A DENSE BASELINE ON BOTH SIDES, and that
        is a property of how they were taken rather than a choice.  Until
        2026-09-02 this path reached past `linearsolver=` to
        `scipy.linalg.lu_factor`, and the dense propagation reached past it
        to `toolkit.linearsolver`; both go through the caller's solver now,
        so the comparison can be re-asked with a sparse one -- and it has
        NOT been, because the box was shared with another session
        benchmarking at load 31-43 when it was tried.  A circuit Jacobian at
        m=1002 is very sparse, so the m~250 gate and these ratios may move.
        `benchmarks/pss_matrix_free_sparse.py` is the harness; run it quiet
        before quoting any of this as the sparse answer.

        ALL FOUR SYSTEMS are converted.  Measured end to end,
        single-threaded, against a DENSE-solver dense path:

          driven solved-history (2m columns)   1.36x/1.51x/2.13x  m=242/502/1002
          driven plain          (m columns)    1.10x/1.34x/1.63x  m=242/502/1002
          autonomous plain      (m + border)   1.02x/1.11x/1.34x  m=128/308/608
          autonomous composed   (2m + border)  1.22x/1.65x        m=208/408

        The plain path's share is about HALF the solved-history path's --
        42.5% against 63.8% at m=502 -- which is what `m` columns instead of
        `2m` against the same assembly has to mean, and it is why the two
        `2m` systems win earlier and by more.  Every one of them agrees with
        the dense path to <= 2e-16, and both autonomous ones reproduce the
        period exactly.
        """
        self._solve_kwargs = dict(refnode=refnode, maxiterations=maxiterations,
                                  matrix_free=matrix_free, tstab=tstab)
        self._monodromy_twin = None
        self._twins = {}
        ## ⚠ HIDDEN STATE IS REFUSED, NOT INTEGRATED AND HOPED OVER.
        ## `TLine.history` is filled by `cir.accept_step`, which the
        ## TRANSIENT calls at every accepted step and which this analysis
        ## never calls -- PSS drives `solve_timestep` directly.  With the
        ## buffer empty `TLine.G`/`u` stamp a DC SHORT, so the line is
        ## silently absent: measured on a quarter-wave open stub, PSS
        ## returned `converged = True`, `spectral_radius = 0.0`, NO warning,
        ## and an amplitude of 0.999969 where a transient gives 0.244201.
        ##
        ## ⚠ AND FILLING THE BUFFER IS NOT THE FIX.  Calling `accept_step`
        ## per step would populate it and make `phi` genuinely
        ## history-dependent -- so the monodromy would be the derivative of
        ## a neighbouring problem, which is the exact failure `_begin_period`
        ## exists to prevent.  `_begin_period` resets what is IN `x`; that is
        ## the right scope for the integrator rings and the wrong one for
        ## state living outside the vector, and no reset of the rings can
        ## fix a period map that is not a function of `x_0`.  The honest
        ## answers are to admit the delay state into the unknowns (a
        ## different analysis) or to refuse; this refuses.
        ##
        ## ⚠ AND THE SCOPE OF THAT REFUSAL IS THE ELEMENT, NOT THE CLASS.
        ## Two things look alike and are not.  KUNDERT'S HIDDEN STATE is a
        ## behavioural model carrying internal state the simulator does not
        ## know about -- genuinely broken, and a commercial RF simulator "outlaws [them]
        ## outright".  A DISTRIBUTED COMPONENT has a KNOWN
        ## infinite-dimensional structure described by frequency-dependent
        ## Y/Z/S parameters, and is tractable: "the convolution operation is
        ## diagonalized by the Fourier transform", so the component is
        ## applied spectrally while the STATE stays finite and lumped
        ## (Yang & Phillips, DAC 2002), and an autonomous time-domain
        ## steady-state solve with transmission lines and exact period
        ## derivatives exists independently.
        ##
        ## `TLine` here trips the first test because of HOW IT IS
        ## IMPLEMENTED -- a `history` buffer filled by `accept_step` -- not
        ## because a transmission line is unshootable.  The flag is opt-in
        ## per element (`Circuit.hidden_state` defaults False), so nothing
        ## refuses distributed components as a class, and an element that
        ## declared its state properly would pass.  Cited, not verified
        ## here.
        _hidden = self.cir.hidden_state_elements()
        if _hidden:
            raise NotImplementedError(
                'PSS: these elements carry HIDDEN STATE -- %s -- so THIS '
                'formulation cannot solve this circuit. The period map must '
                'be a function of x_0 alone, and they stamp from state that '
                'lives outside x and is filled by accept_step, which only a '
                'forward transient calls. Left alone the answer would be '
                'silently wrong rather than slow: an empty TLine history '
                'stamps the line as a DC short and the solve reports '
                'converged. Use Transient for this circuit. ⚠ This is a '
                'limit of the ELEMENT as implemented here, not of shooting '
                'or of distributed components as a class: a component with '
                'a KNOWN frequency-domain description (a transmission line, '
                'an S-parameter block) is tractable in a time-domain '
                'steady-state solve by either admitting the delay state '
                'into the unknowns, or applying the component spectrally -- '
                'the Fourier transform diagonalises the convolution, so its '
                'action becomes a multiply by Y/Z/S while the state stays '
                'finite. Both are different analyses than this one.'
                % ', '.join(sorted(_hidden)))

        self.period = period
        toolkit = self.toolkit

        ## ⚠ ONE REFERENCE NODE PER ANALYSIS, CHECKED.  `self.irefnode` is
        ## fixed in `__init__` from `irefnode=` and is what the TRAVERSAL
        ## uses -- `_transient.irefnode`, every `remove_row_col`, the
        ## monodromy's shape.  This local one comes from `solve`'s own
        ## `refnode=` and is what reinserts the zero row into the RESULT.
        ## They were never compared, so `PSS(cir).solve(refnode=b)` solved
        ## against ground and reported against `b`: each row sensible on its
        ## own, the set of them incoherent, with ground itself coming back
        ## non-zero.  Refused rather than silently rotated, because there is
        ## no answer to give -- the two choices disagree about which
        ## variable was eliminated before the solve began.
        irefnode = self.cir.get_node_index(refnode)
        if irefnode != self.irefnode:
            raise ValueError(
                'PSS: solve(refnode=...) names a different reference node '
                '(index %d) than the analysis was constructed with (index '
                '%d). The traversal eliminated one and the result would '
                'reinsert the other, so the waveform would be reported '
                'against a node the solve never used -- ground itself comes '
                'back non-zero. Pass the same node to both, or construct '
                'the analysis with PSS(cir, irefnode=...) and leave '
                "solve()'s refnode at its default."
                % (irefnode, self.irefnode))
        ## ⚠ CLEARED BEFORE THE RUN, not after it.  These describe the
        ## period this call is about to solve for; leaving the previous
        ## call's behind would let `factored_period()` hand back an operator
        ## for the LAST solve after this one failed, and `converged` alone
        ## would not catch it -- a second solve that fails leaves the first
        ## solve's `converged=True` nowhere in sight but its state very much
        ## in reach.
        self._period_state = None
        self._factored_period_cache = None
        self.waveform = None

        ## ⚠ `theta`'s bias is PER-PERIOD (`_theta_biased`), so it must be
        ## known before anything builds or reuses the inner transient.
        ## `_new_transient` only runs when there is no cache, so a SECOND
        ## `solve()` at a different period would otherwise silently keep the
        ## first one's bias -- re-bias the cached integrator in place rather
        ## than dropping the cache, which would rebuild a `Transient` per
        ## solve for every method that does not care.  A no-op for all of
        ## them (`_theta_biased` type-checks, and tolerates `None`).
        self._theta_period = float(period)
        _tr_cached = getattr(self, '_tran', None)
        if _tr_cached is not None:
            self._theta_biased(getattr(_tr_cached.par, 'integrator', None))

        ## Everything `grid_error` needs to repeat THIS solve on a finer grid.
        ## Recorded rather than re-derived so the refinement differs from the
        ## original in the timestep and in nothing else.
        self._solve_args = dict(refnode=refnode, period=period, x0=x0,
                                timestep=timestep,
                                maxiterations=maxiterations, grid=grid,
                                matrix_free=matrix_free,
                                x0_unknown=x0_unknown, tstab=tstab,
                                break_events=break_events)

        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = toolkit.zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !


        #create vector with timepoints and a more fitting dt
        ## ⚠ the flag must be set BEFORE the grid is built, because
        ## `_period_grid` consults it to decide whether to subdivide a
        ## coarse opening step -- see the note there.
        ## ⚠ `None` means "decide from the topology" -- see
        ## `_resolve_x0_unknown`.  Resolved to a concrete bool HERE, before
        ## anything reads it, so every downstream use sees one value.
        x0_unknown = self._resolve_x0_unknown(x0_unknown)
        if self._integrator_for(getattr(self.par, 'method', 'euler')
                                ).needs_x0_unknown():
            ## Self-starting stage methods: `x_in` IS `x_0`, there is no
            ## manufacturing step to differentiate `x_0` back through, so the
            ## unknown is always `x_0` itself.  The method says so
            ## (`needs_x0_unknown`), which keeps the phase pin and every
            ## open-at-x0 branch consistent without a name check here.
            x0_unknown = True
        self._open_at_x0 = bool(x0_unknown)
        ## Break the traversal's steps at the source discontinuities, when the
        ## method is one whose accuracy that helps -- see `_resolve_break_events`
        ## for the measurement, and `event_grid` for the snap that keeps it from
        ## manufacturing slivers.
        self.break_events = self._resolve_break_events(break_events)
        if self.break_events:
            _ev = (self.event_grid(period, grid=grid) if grid is not None
                   else self.event_grid(period, npts=int(period / dt)))
            ## ⚠ ONLY replace the grid when there ARE events.  `event_grid`
            ## rebuilds a uniform grid from `linspace` even when it finds none,
            ## and that differs from `_period_grid`'s own in the last bit --
            ## enough to move every event-free solve in the suite for nothing.
            ## Touching the grid only when an event exists keeps every circuit
            ## without one BIT-IDENTICAL, the same guarantee `_fold_periodic`
            ## gives a circuit with no periodic state.
            if self.event_times:
                grid = _ev
        times, hs = self._period_grid(period, int(period / dt), grid)
        npts = len(times)
        self._grid_fracs = (None if grid is None
                            else np.asarray(grid, dtype=float))
        ## read by `_period_grid`, which is called from the residual
        ## closures and so cannot take it as an argument
        self._open_at_x0 = bool(x0_unknown)
        ## The fold gauge, collected once per solve (late-bound moduli are
        ## resolved by now).  See `_fold_periodic` for why the residual needs
        ## it and the Jacobian does not.
        self._periodic_fold = self._collect_periodic_fold()
        alpha = 1

        ## AUTONOMY IS DECIDED BEFORE THE SOLVE, because it decides which
        ## system is solved.  Structural and exact -- see `_is_autonomous`.
        self.autonomous = self._is_autonomous(times)
        phase_k, phase_pin = 0, 0.0
        if self.autonomous:
            ## An unseeded autonomous run starts at the origin, which IS a
            ## periodic solution -- the trivial one -- and the free-period
            ## system would sit there just as contentedly as the fixed one
            ## did.  The operating point is the honest default: for a phase
            ## accumulator `ic` pins it on the orbit.
            ##
            ## ⚠ `tstab` RUNS AFTER THIS, NOT BEFORE, and the order is the
            ## whole of what makes it work.  A pre-integration seeded from
            ## `zeros` starts AT the equilibrium of an autonomous circuit and
            ## a transient from an exact equilibrium never leaves it, so the
            ## warm start would return the basin it was asked to escape.
            ## Starting it from the operating point is the honest version of
            ## the same statement.
            if x0 is None:
                from pycircuit.circuit.dcanalysis import DC
                xdc = np.asarray(DC(self.cir, toolkit=self.toolkit).solve().x,
                                 dtype=float).reshape(-1)
                x = np.concatenate((xdc[:irefnode], xdc[irefnode + 1:]))

        if tstab:
            ## ⚠ THE PRE-INTEGRATION IS A PLAIN TRANSIENT, and it has to be:
            ## its whole value is that it is NOT a shooting solve, so it
            ## cannot be captured by the basin that the shooting Newton is
            ## stuck in.  It runs on its own adaptive grid -- `timestep` is
            ## a first step, not an imposed one -- because nothing here
            ## needs `phi` to be a function of `x_0`; that requirement
            ## starts when the shooting does.
            from pycircuit.circuit.transient import Transient
            _xred = np.asarray(x, dtype=float).reshape(-1)
            _xfull = np.concatenate((_xred[:irefnode],
                                     np.zeros(1), _xred[irefnode:]))
            _pre = Transient(
                self.cir, toolkit=self.toolkit, reltol=self.par.reltol,
                iabstol=self.par.iabstol, vabstol=self.par.vabstol,
                nrsolver=self.par.nrsolver,
                linearsolver=self.par.linearsolver, scaler=self.par.scaler)
            _res = _pre.solve(refnode=refnode, tend=float(tstab),
                              timestep=dt, x0=_xfull)
            _last = np.asarray(_res.x, dtype=float)[:, -1]
            x = np.concatenate((_last[:irefnode], _last[irefnode + 1:]))
            self.tstab_state = x

        if self.autonomous:
            ## ⚠ AND IT IS SUFFICIENT IN PRINCIPLE, DEGRADED IN PRACTICE ON
            ## HIGH-Q: the row removes the singularity from the UNIT
            ## multiplier and only that one, so an oscillator whose other
            ## multipliers cluster near 1 gives a bordered system that is
            ## nonsingular and ill conditioned.  Same cause as the
            ## eigen-selection limit in the class docstring; read it there.
            ##
            ## THE PHASE CONDITION pins the coordinate moving FASTEST at the
            ## seed, so the orbit crosses the pinning hyperplane
            ## transversally.  Pin a slow one and the last row of the
            ## bordered Jacobian is nearly parallel to the null direction it
            ## exists to remove, which is a singular system wearing an extra
            ## equation.
            ##
            ## ⚠ THE `argmax` COMPARES VOLTS WITH AMPERES ON PURPOSE, and
            ## the obvious repair is WRONG.  The row can only remove the
            ## orbit's tangent in proportion to `|e_k . fhat|`; one step of
            ## `|dx_k|` IS `|f_k|` up to `h`, so this argmax is exactly
            ## `argmax |e_k . fhat|` -- it maximises the quantity the row
            ## needs.  Normalising each coordinate by its own swing, which
            ## is what the vector's mixed units invite, was MEASURED on a
            ## van der Pol carrying a VCVS-scaled copy of `v` and picks a
            ## row 704x WORSE aligned (1.4e-03 against 1.0000).  The
            ## scaling that lets a large coordinate win the argmax is the
            ## same scaling that makes it dominate `f`; the two cancel.
            ## Pinned by `test_the_phase_pin_compares_units_on_purpose`.
            ##
            ## ⚠ A REVIEW DISPUTED THESE CONDITION NUMBERS AND THEN
            ## RETRACTED, and the reason is worth keeping because it is a
            ## trap this file can fall into again.  The reviewer measured
            ## the pin and the orthogonality row as IDENTICAL (edge 1.000x)
            ## and argued a gap was arithmetically impossible, since
            ## `|fhat[k]| = 0.9999` makes `e_k` and `fhat` nearly parallel.
            ## Their harness took the CODE'S analytic `J` and swapped only
            ## the border row, so both readings were the conditioning of
            ## the SAME operator -- and on the plain path that operator is
            ## in the wrong frame (see item 4b's correction below), so its
            ## defect dominated `cond` in both cases.  Identical numbers
            ## were guaranteed by construction.
            ##
            ## ⚠ THE REASONING WAS ALSO WRONG, INDEPENDENTLY: `cond` is a
            ## function of the WHOLE row, not of its projection on `e_k`.
            ## Alignment 0.9999 makes two rows nearly PARALLEL, not equal --
            ## `e_k` is a unit vector and `fhat` is dense -- and it places
            ## no bound on the conditioning gap.  Switching the phase
            ## condition means solving a DIFFERENT system and living with
            ## ITS conditioning, so the comparison has to build both
            ## systems, which is what the numbers above do.
            ##
            ## ⚠ AND THE UPGRADE THIS INVITED WAS TESTED AND REJECTED.  An
            ## orthogonality (Poincare) row `<x0 - x_ref, f(x_ref)> = 0`
            ## looks strictly better -- it is the flow-aligned row by
            ## construction, and it cannot pin an unattainable VALUE.
            ## Measured against this rule on the case built to break it
            ## (seeded at `v`'s turning point, so the pin sits 1e-3 of the
            ## way into its coordinate's range): both converge, to the same
            ## answer, at every grid tried, and the bordered condition
            ## number is 1.2e3/3.0e2/8.2e1 for the pin against
            ## 2.0e2/6.0e1/3.7e1 for orthogonality at 200/800/3200 points --
            ## a 2-6x edge that never decides anything, and it SHRINKS as
            ## the grid refines.  The row's alignment at the solution stayed
            ## 1.0000 throughout: the tangency this was supposed to induce
            ## never materialised.
            ##
            ## What IS real: `phase_pin` is a VALUE the orbit must attain,
            ## so a seed far off the orbit can pin one outside its range and
            ## the system is then INCONSISTENT rather than merely hard --
            ## measured on van der Pol at mu=1, seeds at 4x/10x/30x the
            ## orbit amplitude pin `v` at -5.66/-9.50/-15.46 against an
            ## orbit range of [-2.01, 2.01], and all three report ordinary
            ## non-convergence.  ⚠ That is NOT an argument for the
            ## orthogonality row: its plane through the same far seed misses
            ## the orbit too (checked).  It is an argument about SEEDS, and
            ## the remaining gap is diagnostic, not formulational.
            self._begin_period(x)
            _x1 = self.solve_timestep(x, times[0], hs[0])
            _x2 = self.solve_timestep(_x1, times[1], hs[0], iq_last=self._iq)
            phase_k = int(np.argmax(np.abs(np.asarray(_x2) - np.asarray(_x1))))
            ## ⚠ THE RULE ITSELF IS CANONICAL.  Aprille & Trick's
            ## oscillator paper, Step 3: "select k by
            ## |f_k(x^i(T^i))| = max_k |f_k(x^i(T^i))|" -- argmax of the
            ## vector field, which is what an argmax over one step is up to
            ## `h`, and `h` is common to every coordinate so the mixed units
            ## cancel exactly.  Both the rule and the decision not to
            ## replace it with a Poincare row have precedent as well as
            ## measurement behind them.
            ##
            ## ⚠ WHAT IS NOT CANONICAL IS FREEZING IT.  Their Step 3 sits
            ## INSIDE the iteration (Step 5 returns to Step 1), so `k` and
            ## the pinned value are re-chosen from the CURRENT trajectory
            ## every iterate -- "note that in this method, an initial k and
            ## w_0k are not required".  Pinning `w_0k = x_0k^i`, the
            ## iterate's OWN current value, is attainable by construction,
            ## which removes the failure mode measured below (a far seed
            ## pinning a value the orbit never reaches) structurally rather
            ## than by advice.  Not done here: `analysis.fsolve` exposes no
            ## per-iteration hook, so it needs the autonomous outer solve
            ## restructured, and re-selecting BETWEEN outer iterations does
            ## not threaten `phi` being a function of `x_0` alone.
            ##
            ## ⚠ THE PHASE ROW SITS OUTSIDE THE INTEGRATOR ON PURPOSE, and
            ## that placement is load-bearing rather than incidental.
            ## Brachtendorf et al. (TCAD 33(6) 867-878) warn of the
            ## alternative: "adding an algebraic equation transforms the
            ## system of (implicit) ODEs to a system of DAEs.  Transient
            ## methods may run into severe problems when the index of a
            ## system of DAEs is two or higher."  This augments the OUTER
            ## shooting system; the inner integration is unaugmented, so the
            ## phase condition cannot raise the index of the DAE actually
            ## being integrated -- which matters here in proportion to how
            ## badly index 2 already behaves (trap and euler both fail to
            ## converge on a V-source/C/C/R loop where gear does).
            ##
            ## ⚠ THE RULE IS CANONICAL; RE-SELECTING IT WAS TRIED AND
            ## REJECTED, WITH NUMBERS.  Aprille & Trick's oscillator paper
            ## picks `k` by `argmax |f_k(x^i(T^i))|` -- the same quantity an
            ## argmax over one step is, up to an `h` common to every
            ## coordinate -- and their Step 3 sits INSIDE the loop, pinning
            ## `w_0k = x_0k^i`, the iterate's OWN value, so that "an initial
            ## k and w_0k are not required".  That looks like a free repair
            ## for the far-seed failure recorded below.  It is not.
            ##
            ## Built and measured 2026-09-02: re-selecting `k` and the pin
            ## from the current trajectory between outer iterations REGRESSES
            ## the working case.  Van der Pol at mu=1 from an ON-ORBIT seed
            ## went from converged to NOT converged, and far seeds wandered
            ## to periods of -52, -1088 and +110 against a true 6.6633.
            ##
            ## ⚠ AND THE REASON IS STRUCTURAL, not a bug in the attempt.
            ## Pinning the iterate's own value makes the phase residual
            ## `x_0[k] - pin` IDENTICALLY ZERO at every iterate, so the row
            ## carries no information: it constrains the STEP (`dz[k] = 0`)
            ## and nothing else, and with `k` re-chosen each iteration a
            ## different coordinate is frozen each time, so the orbit slides
            ## along itself.
            ##
            ## A&T do not have this problem because THEY HAVE NO PHASE
            ## EQUATION.  Their unknown vector SUBSTITUTES the period for the
            ## pinned coordinate -- `v = [x_01, ..., x_0(k-1), T,
            ## x_0(k+1), ..., x_0n]` -- an n x n system in which `x_0k` is a
            ## CONSTANT rather than an unknown with a trivially satisfied
            ## equation.  The constraint is structural where ours is
            ## algebraic.  So Step 3 is not portable to a bordered
            ## formulation as a drop-in: taking it means taking the
            ## substitution with it.
            ##
            ## The frozen pin's failure mode below therefore STANDS, and its
            ## fix is that substitution, not a moving pin.
            ##
            ## ⚠ THE PIN MUST BE IN THE UNKNOWN'S OWN FRAME.  `_x1` is the
            ## state one step AFTER the seed, which is the right thing to
            ## pin when the unknown is `x_in` and `x_0` is manufactured from
            ## it -- and the wrong thing when the unknown IS `x_0`.  With a
            ## fine opening step the two are nearly equal and the mismatch
            ## hides; on van der Pol's own LTE grid, where the opening step
            ## is 1.4845 against a median of 4.6e-04, it pins a value the
            ## orbit need never attain and the solve dies with a bare
            ## non-convergence.
            phase_pin = float(np.asarray(x if x0_unknown else _x1)[phase_k])

        ## Resolved here as well as in `solve_timestep`, because the SHOOTING
        ## Jacobian depends on which integrator the inner steps used.
        ##
        ## ⚠ THE NAME IS VALIDATED BEFORE ANYTHING ASKS THE INTEGRATOR A
        ## QUESTION.  `_solves_history` resolves `method` to a class to ask how
        ## far its companion reaches; run first, it turned an unknown name
        ## into a `KeyError` from a dict several frames down, in place of the
        ## `ValueError` this raises.  Two tests caught it, both written for
        ## the class's earlier fall-through defects.
        method = getattr(self.par, 'method', 'euler')
        if method not in ('euler', 'trap', 'trapezoidal', 'theta', 'gear',
                          'gear2', 'trbdf2', 'radau', 'esdirk43'):
            raise ValueError(
                "method must be 'euler', 'trap', 'theta', 'gear', 'trbdf2', "
                "'radau' or 'esdirk43', not %r" % (method,))

        ## Whether the entering history joins the unknowns.  Decided once,
        ## here, because it chooses which system is solved -- like autonomy,
        ## and after it, since the two are not composed yet.
        solved_history = self._solves_history()
        self.solved_history = solved_history
        xm1_ss = None

        ## THE SHOOTING JACOBIAN FOLLOWS THE INTEGRATOR'S OWN COEFFICIENTS.
        ##
        ## Backward Euler's per-step sensitivity is
        ##     dx_n/dx_{n-1} = Jf_n^-1 * C(x_{n-1})/h
        ## -- the COMPANION CONDUCTANCE at the previous point, not the raw
        ## capacitance matrix.  The `/h` was missing, and because C is
        ## singular the accumulated product collapsed to EXACTLY ZERO: the
        ## Jacobian handed to fsolve was `I - 0 = I`, so the "shooting
        ## Newton" was plain successive substitution `x0 <- phi(x0)`.  That
        ## converges at the circuit's own per-period decay -- measured on a
        ## Q=20 resonator as 0.855 per iteration against exp(-pi/Q) = 0.8546,
        ## which is how it was found -- and it never reached fsolve's
        ## tolerance, on any circuit, silently.  With the companion
        ## conductance the same resonator converges in FIVE iterations
        ## (2.64 -> 2.6e-2 -> 1.0e-2 -> 1.9e-4 -> 3.9e-5).
        ##
        ## TRAPEZOIDAL'S MONODROMY MUST CARRY `iq` AS WELL AS `x`.  Its
        ## recursion carries `iq` as well as `x`, so the period map is a
        ## function of (x, iq); an x-only monodromy is structurally
        ## incomplete, and measured, using the Euler form for trap converged
        ## SLOWER than no Jacobian at all (0.90 against 0.855 per iteration).
        ##
        ## Differentiating the two recursions together,
        ##
        ##     iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}
        ##     0    = i(x_n) + iq_n + u(t_n)
        ##
        ## gives a propagation of `d(x,iq)/dx0` that costs one extra matrix
        ## product over the Euler form:
        ##
        ##     rhs   = Geq_{n-1} Px + Pq
        ##     Px_n  = Jf_n^-1 rhs
        ##     Pq_n  = Geq_n Px_n - rhs
        ##
        ## Euler is the SAME recursion with `Pq == 0`: its companion carries
        ## no `iq_{n-1}` term, so dF/diq_{n-1} vanishes and the second row
        ## never enters.  One formula, two methods, which is why this is not
        ## a second code path.
        ## (`newton = True` lived here and fed a message branch that is
        ## gone: it was set unconditionally, so the 'successive
        ## substitution' alternative was dead and the claim it selected was
        ## false anyway -- see the non-convergence warning below.)

        def func(x):
            x0, x_end, Mx, _Mt = self._traverse(x, period, times, hs,
                                                want_dT=False,
                                                open_at_x0=x0_unknown)
            D = np.asarray(toolkit.eye(n - 1))
            return self._fold_periodic(x0 - x_end), D - alpha * Mx

        def func_solved_history(z):
            """Residual and Jacobian when the entering history is an unknown.

            Unknowns are `(x_0, x_{-1})`; the equations are that BOTH close,

                F = [ x_0 - x_{N-1} ,  x_{-1} - x_{N-2} ]
                J = [[ I - A(N-1,0) , -A(N-1,-1) ],
                     [   -A(N-2,0)  , I - A(N-2,-1) ]]

            with `A(j,k) = d x_j / d x_k`, which is what `_traverse_solved_history`
            returns as one `n x 2n` block per row.  A two-step companion needs
            two states to be continued, so periodicity of ONE of them is an
            under-determined statement about the orbit -- which is the defect
            this replaces, measured at 1.266e-01 V for Gear-2 at 100 points
            per period.
            """
            m = n - 1
            x0_in, xm1_in = z[:m], z[m:]
            x_last, x_prev, P_last, P_prev = self._traverse_solved_history(
                x0_in, xm1_in, times, hs)

            D = np.asarray(toolkit.eye(m))
            J = np.zeros((2 * m, 2 * m))
            J[:m, :m] = D - alpha * P_last[:, :m]
            J[:m, m:] = -alpha * P_last[:, m:]
            J[m:, :m] = -alpha * P_prev[:, :m]
            J[m:, m:] = D - alpha * P_prev[:, m:]
            F = np.concatenate((
                self._fold_periodic(np.asarray(x0_in) - np.asarray(x_last)),
                self._fold_periodic(np.asarray(xm1_in) - np.asarray(x_prev))))
            return F, J

        def func_autonomous(z):
            """Residual and Jacobian of the FREE-PERIOD system.

            Unknowns are `(x0, T)`; equations are the period map's fixed
            point plus a phase condition, because without one the system is
            singular by construction -- every point on the orbit is a
            solution, so `I - M` has a null direction along it.

                F = [ x0 - phi_T(x0) ,  x0[k] - pinned ]
                J = [[ I - M , -dphi/dT ],
                     [ e_k^T ,     0    ]]

            `k` is the coordinate moving fastest at the seed, so the orbit
            crosses the pinning hyperplane transversally; pinning a slow
            coordinate makes the last row nearly parallel to the null
            direction it is there to remove.
            """
            x_in, T = z[:-1], float(z[-1])
            ## Rebuilt at the CURRENT T, which is what keeps `dh/dT = h/T`
            ## true of every step -- uniform or not.
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            x0, x_end, Mx, Mt = self._traverse(x_in, T, tms, hs_T,
                                               want_dT=True,
                                               open_at_x0=x0_unknown)

            D = np.asarray(toolkit.eye(n - 1))
            m = n - 1
            J = np.zeros((m + 1, m + 1))
            J[:m, :m] = D - alpha * Mx
            J[:m, m] = -np.asarray(Mt).ravel()
            J[m, phase_k] = 1.0
            F = np.zeros(m + 1)
            F[:m] = self._fold_periodic(x0 - x_end)
            F[m] = x0[phase_k] - phase_pin
            return F, J
        
        def func_autonomous_solved_history(z):
            """Both enlargements at once: `(x_0, x_{-1}, T)`.

            The free-period system and a solved entering history grow the
            SAME unknown vector in different directions, and a two-step
            method on an autonomous circuit needs both -- the period because
            it is not given, the history because the companion reads it.

                F = [ x_0  - x_{N-1} ,  x_{-1} - x_{N-2} ,  x_0[k] - pin ]

                J = [[ I - A(N-1,0) ,   -A(N-1,-1) , -dx_{N-1}/dT ],
                     [   -A(N-2,0)  , I - A(N-2,-1), -dx_{N-2}/dT ],
                     [     e_k^T    ,       0      ,       0      ]]

            ⚠ ONE PHASE ROW STILL SUFFICES, and it pins only the `x_0`
            block.  Time translation slides BOTH states along the orbit
            together -- the null vector is `(xdot(0) ds, xdot(-h) ds, dT)` --
            so the freedom stays one-dimensional and pinning `x_0[k]` kills
            it whenever `xdot(0)[k] != 0`.  That is the same reason `k` is
            chosen as the fastest-moving coordinate at the seed, and the
            same reason a slow one would leave the last row nearly parallel
            to the direction it exists to remove.

            ⚠ THE HISTORY POINT MOVES WITH T.  `x_{-1}` sits at `-T/(N-1)`,
            so its location tracks the unknown.  The residual is still the
            right statement -- `x_{-1}` and `x_{N-2}` are the same phase of
            the orbit at every T -- and the T column is just the propagation
            accumulated to step N-2, which the ring already carries.
            """
            m = n - 1
            x0_in, xm1_in, T = z[:m], z[m:2 * m], float(z[-1])
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            (x_last, x_prev, P_last, P_prev, Pt_last,
             Pt_prev) = self._traverse_solved_history(
                x0_in, xm1_in, tms, hs_T, T=T, want_dT=True)

            D = np.asarray(toolkit.eye(m))
            J = np.zeros((2 * m + 1, 2 * m + 1))
            J[:m, :m] = D - alpha * P_last[:, :m]
            J[:m, m:2 * m] = -alpha * P_last[:, m:]
            J[:m, 2 * m] = -np.asarray(Pt_last).ravel()
            J[m:2 * m, :m] = -alpha * P_prev[:, :m]
            J[m:2 * m, m:2 * m] = D - alpha * P_prev[:, m:]
            J[m:2 * m, 2 * m] = -np.asarray(Pt_prev).ravel()
            J[2 * m, phase_k] = 1.0

            F = np.zeros(2 * m + 1)
            F[:m] = np.asarray(x0_in) - np.asarray(x_last)
            F[m:2 * m] = np.asarray(xm1_in) - np.asarray(x_prev)
            F[2 * m] = np.asarray(x0_in)[phase_k] - phase_pin
            return F, J

        def func_full(x):
            """Driven fixed-period residual and Jacobian for Radau IIA(3).

            `x` IS `x_0` (self-starting), so `F = x_0 - phi(x_0)` and
            `J = I - M` with `M` the dense coupled monodromy from
            `_traverse_full`.  No manufacturing step, no solved history.
            """
            x0, x_end, Mx, _Mt = self._traverse_full(
                x, period, times, hs, want_dT=False)
            D = np.asarray(toolkit.eye(n - 1))
            return (self._fold_periodic(np.asarray(x0) - np.asarray(x_end)),
                    D - alpha * Mx)

        def func_autonomous_full(z):
            """Free-period residual and Jacobian for Radau IIA(3).

            Unknowns `(x0, T)`; `F = [x0 - phi_T(x0), x0[k] - pinned]`,
            `J = [[I - M, -dphi/dT], [e_k^T, 0]]`.  The period column
            `dphi/dT` comes from `_traverse_full(want_dT=True)`, tractable
            because the circuit is autonomous.
            """
            x_in, T = z[:-1], float(z[-1])
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            x0, x_end, Mx, Mt = self._traverse_full(
                x_in, T, tms, hs_T, want_dT=True)
            m = n - 1
            D = np.asarray(toolkit.eye(m))
            J = np.zeros((m + 1, m + 1))
            J[:m, :m] = D - alpha * Mx
            J[:m, m] = -np.asarray(Mt).ravel()
            J[m, phase_k] = 1.0
            F = np.zeros(m + 1)
            F[:m] = self._fold_periodic(np.asarray(x0) - np.asarray(x_end))
            F[m] = np.asarray(x0)[phase_k] - phase_pin
            return F, J

        def func_dirk(x):
            """Driven fixed-period residual/Jacobian for a lower-triangular
            (DIRK/ESDIRK) stage method: `F = x0 - phi(x0)`, `J = I - M`, with
            `M` the sequential dense monodromy from `_traverse_dirk`."""
            x0, x_end, Mx, _Mt = self._traverse_dirk(
                x, period, times, hs, want_dT=False)
            D = np.asarray(toolkit.eye(n - 1))
            return (self._fold_periodic(np.asarray(x0) - np.asarray(x_end)),
                    D - alpha * Mx)

        def func_autonomous_dirk(z):
            """Free-period residual/Jacobian for a DIRK/ESDIRK method; the
            period column comes from `_traverse_dirk(want_dT=True)`."""
            x_in, T = z[:-1], float(z[-1])
            tms, hs_T = self._period_grid(T, npts, self._grid_fracs)
            x0, x_end, Mx, Mt = self._traverse_dirk(
                x_in, T, tms, hs_T, want_dT=True)
            m = n - 1
            D = np.asarray(toolkit.eye(m))
            J = np.zeros((m + 1, m + 1))
            J[:m, :m] = D - alpha * Mx
            J[:m, m] = -np.asarray(Mt).ravel()
            J[m, phase_k] = 1.0
            F = np.zeros(m + 1)
            F[:m] = self._fold_periodic(np.asarray(x0) - np.asarray(x_end))
            F[m] = np.asarray(x0)[phase_k] - phase_pin
            return F, J

        ## THE SHOOTING RESIDUAL IS IN SOLUTION UNITS, NOT KCL UNITS.
        ## `x0 - phi(x0)` is a difference of SOLUTIONS -- volts on node rows,
        ## amps on branch rows -- so its absolute floor is the `xtol` flavour
        ## (vabstol on nodes, iabstol on branches), not the residual flavour
        ## the transient's Newton uses for `i(x)`.  Getting that backwards is
        ## F6(a)'s defect, and it is easy to walk into here because the
        ## quantity is called a residual.
        _tol = analysis.newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[1]
        (_tol,) = remove_row_col((_tol,), irefnode, self.toolkit)

        ## The shooting criterion, expressed against the transient one.
        _ratio = float(self.par.steadyratio)
        if _ratio < 1.0:
            raise ValueError(
                'steadyratio must be >= 1 (got %g): the period map is only '
                'known to the accuracy of the per-timestep solves, so a '
                'shooting tolerance tighter than reltol asks the outer '
                'residual to resolve its own noise.' % _ratio)
        _shoot_reltol = self.par.reltol * _ratio
        _tol = _tol * _ratio

        ## ⚠ REFUSED RATHER THAN SILENTLY IGNORED.  A caller asking for
        ## matrix-free on a path that has not got it wants the cost model it
        ## implies; quietly taking the dense route would be a performance
        ## surprise with no symptom, which is the shape of defect this tree
        ## has paid for before.
        ## ⚠ REFUSED RATHER THAN IGNORED, like `matrix_free` above.  A
        ## solved-history method already solves for `x_0` and `x_{-1}`
        ## directly and manufactures nothing, so there is no frame to
        ## correct and the flag would be a no-op -- and a no-op flag that
        ## the caller believes changed something is worse than an error.
        if x0_unknown and solved_history:
            raise NotImplementedError(
                'PSS: x0_unknown=True has nothing to change for a two-step '
                "method (method=%r). The solved-history formulation already "
                'solves for x_0 and x_{-1} as real trajectory states and '
                'manufactures no opening step, so its Jacobian is exact '
                'without this. Drop the flag, or use a one-step method '
                "(method='trap' or 'euler') where the manufactured opening "
                'is what this replaces.' % method)

        ## ⚠ THE OUTER NEWTON IS DAMPED, which it was not.  All three
        ## `fsolve` calls took the FULL step with `limiter=None` and no line
        ## search -- a departure from standard practice rather than a
        ## neutral choice: Brachtendorf et al. (TCAD 33(6) 867-878) describe
        ## "shooting, finite difference, or harmonic balance techniques IN
        ## CONJUNCTION WITH A DAMPED NEWTON METHOD" as what is "widely
        ## employed" for limit cycles.  The full step is still tried first
        ## and kept whenever it improves the residual, so a solve that was
        ## converging is unchanged; the halving only runs where the
        ## undamped iteration would have moved uphill.

        ## Find periodic steady state x-vector
        if self._integrator_for(method).is_stage_method():
            ## A self-starting stage method: its own dense monodromy, never
            ## solved-history, always self-starting (x0 is the unknown).  Route
            ## BY STRUCTURE -- fully-implicit uses the coupled dense traverse
            ## (`func_full`), a DIRK/ESDIRK the sequential one (`func_dirk`).
            ## Both are tableau-generic; a new method of either family reaches
            ## the right one with no edit here.
            _fully = self._integrator_for(method).is_fully_implicit()
            _fdr = func_full if _fully else func_dirk
            _fda = func_autonomous_full if _fully else func_autonomous_dirk
            _label = method
            if matrix_free:
                raise NotImplementedError(
                    'PSS: matrix-free shooting is not built for %s; its '
                    'monodromy is a dense stage product. Drop '
                    'matrix_free, or use a one-step LMM.' % _label)
            if self.autonomous:
                xa = np.asarray(x, dtype=float)
                z0 = np.concatenate((xa, [period]))
                abstol_z = np.concatenate((_tol, [_tol[phase_k]]))
                xtol_z = np.concatenate((_tol, [1e-15 * period]))
                z_ss, _info, _ier, _mesg = self._free_period_solve(
                    _fda, z0, abstol_z, xtol_z,
                    _shoot_reltol, maxiterations, period)
                x0_ss = z_ss[:-1]
                self.period = period = float(z_ss[-1])
                times, hs = self._period_grid(period, npts, self._grid_fracs)
            else:
                x0_ss, _info, _ier, _mesg = analysis.fsolve(
                    _fdr, x, maxiter=maxiterations,
                    reltol=_shoot_reltol, abstol=_tol, xtol=_tol,
                    toolkit=self.toolkit, full_output=True, line_search=True)
        elif self.autonomous and solved_history:
            ## BOTH unknowns and the period.  The floors follow the same
            ## rule as the plain autonomous system: the two state blocks
            ## take the solution-unit tolerance, and the row that adds a
            ## TIME as an unknown takes a time as its floor -- mixing them
            ## is flavour error F6(a) one row further out.
            m_ = n - 1
            xa = np.asarray(x, dtype=float)
            z0 = np.concatenate((xa, xa, [period]))
            abstol_z = np.concatenate((_tol, _tol, [_tol[phase_k]]))
            xtol_z = np.concatenate((_tol, _tol, [1e-15 * period]))
            _mfc = None
            if matrix_free:
                ## BOTH ENLARGEMENTS AT ONCE, matrix-free.  With
                ## `w = (v_0, v_{-1}, s)` and `M v` the `2m` pair map,
                ##     J w = [ v - M v - s (Pt_last, Pt_prev) ; v[k] ]
                ## -- one phase row, pinning the `x_0` block only, exactly as
                ## the dense system does and for the same reason.
                self._monodromy = None

                def _build_comp(z):
                    x0_, xm1_, T_ = z[:m_], z[m_:2 * m_], float(z[-1])
                    tms_, hsT_ = self._period_grid(T_, npts, self._grid_fracs)
                    (C0_, st_, xl_, xp_, Ptl_,
                     Ptp_) = self._traverse_factored(
                        x0_, xm1_, tms_, hsT_, T=T_, want_dT=True)
                    Ptv_ = np.concatenate((np.asarray(Ptl_).ravel(),
                                           np.asarray(Ptp_).ravel()))
                    F_ = np.concatenate(
                        (np.asarray(x0_, dtype=float) - np.asarray(xl_, dtype=float),
                         np.asarray(xm1_, dtype=float) - np.asarray(xp_, dtype=float),
                         [float(np.asarray(x0_, dtype=float)[phase_k])
                          - phase_pin]))

                    def mv_(w):
                        v_, s_ = w[:2 * m_], float(w[2 * m_])
                        top = (v_ - alpha * self._monodromy_matvec(C0_, st_, v_)
                               - s_ * Ptv_)
                        return np.concatenate((top, [v_[phase_k]]))
                    return F_, mv_

                def _mfc(z0_, ab_, xt_, rt_, mi_):
                    return self._matrix_free_newton(_build_comp, z0_, ab_,
                                                    xt_, rt_, mi_)
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous_solved_history, z0, abstol_z, xtol_z,
                _shoot_reltol, maxiterations, period, solver=_mfc)
            x0_ss, xm1_ss = z_ss[:m_], z_ss[m_:2 * m_]
            self.period = period = float(z_ss[-1])
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif self.autonomous:
            ## The period joins the unknowns.  Its residual row is the phase
            ## condition -- in the units of the coordinate it pins, hence
            ## `_tol[phase_k]` -- while the UNKNOWN it adds is a time, whose
            ## own floor has to be a time; mixing the two is the flavour
            ## error F6(a) names, one row further out.
            z0 = np.concatenate((np.asarray(x, dtype=float), [period]))
            abstol_z = np.concatenate((_tol, [_tol[phase_k]]))
            xtol_z = np.concatenate((_tol, [1e-15 * period]))
            _mf = None
            if matrix_free:
                ## THE BORDERED SYSTEM, matrix-free.  `dphi/dT` is ONE column
                ## and does not depend on the Krylov direction, so the
                ## trajectory pass computes it once per Newton iteration and
                ## the matvec just uses it:
                ##     J [v; s] = [ (I - M) v - s dphi/dT ; v_k ]
                m_ = n - 1

                def _build_auto(z):
                    xin_, T_ = z[:m_], float(z[-1])
                    tms_, hsT_ = self._period_grid(T_, npts, self._grid_fracs)
                    op_, st_, x0_, xe_, Mt_ = self._traverse_factored_plain(
                        xin_, T_, tms_, hsT_, want_dT=True,
                        open_at_x0=x0_unknown)
                    x0_ = np.asarray(x0_, dtype=float)
                    Mt_ = np.asarray(Mt_, dtype=float).ravel()
                    F_ = np.concatenate((x0_ - np.asarray(xe_, dtype=float),
                                         [x0_[phase_k] - phase_pin]))

                    def mv_(w):
                        v_, s_ = w[:m_], float(w[m_])
                        top = (v_ - alpha * self._monodromy_matvec_plain(
                            op_, st_, v_)) - s_ * Mt_
                        return np.concatenate((top, [v_[phase_k]]))
                    return F_, mv_

                def _mf(z0_, ab_, xt_, rt_, mi_):
                    return self._matrix_free_newton(_build_auto, z0_, ab_,
                                                    xt_, rt_, mi_)
            z_ss, _info, _ier, _mesg = self._free_period_solve(
                func_autonomous, z0, abstol_z, xtol_z, _shoot_reltol,
                maxiterations, period, solver=_mf)
            x0_ss = z_ss[:-1]
            self.period = period = float(z_ss[-1])
            ## The grid follows the solved period; everything downstream --
            ## the replay, the waveform, the DFT -- must use it or the
            ## answer is reported on a period the solver rejected.
            times, hs = self._period_grid(period, npts, self._grid_fracs)
        elif solved_history:
            ## THE SEED IS THE OLD FORMULATION'S ASSUMPTION, written down:
            ## `x_{-1} = x_0`.  It is what the plain path silently assumes
            ## (it seeds both charge rings with the entering state), so an
            ## solved-history run starts where a plain one starts and the
            ## comparison between them is about the SOLVE, not the seed.
            xa = np.asarray(x, dtype=float)
            z0 = np.concatenate((xa, xa))
            tol_z = np.concatenate((_tol, _tol))
            if matrix_free:
                ## ⚠ THE MONODROMY IS NOT FORMED, so it must not be
                ## REPORTED either.  `_monodromy` survives on the object
                ## from any earlier traversal, and `spectral_radius` reads
                ## it below without knowing which run wrote it -- so a
                ## matrix-free solve that left it alone would report the
                ## PREVIOUS solve's radius as this one's.  Cleared here,
                ## and `spectral_radius` is documented as None on this path.
                self._monodromy = None
                z_ss, _info, _ier, _mesg = self._matrix_free_solve(
                    z0, times, hs, tol_z, tol_z, _shoot_reltol,
                    maxiterations)
            else:
                z_ss, _info, _ier, _mesg = analysis.fsolve(
                    func_solved_history, z0, maxiter=maxiterations,
                    reltol=_shoot_reltol, abstol=tol_z, xtol=tol_z,
                    toolkit=self.toolkit, full_output=True, line_search=True)
            x0_ss, xm1_ss = z_ss[:n - 1], z_ss[n - 1:]
        elif matrix_free:
            ## RECORDED SCOPE ITEM 6 on the PLAIN path: `m` columns rather
            ## than `2m`, so a lower share and a lower ceiling than the
            ## solved-history route -- 42.5% and 1.71x at m=502, 60.7% and
            ## 2.50x at m=1002.
            def _build(z):
                opening, steps, x0_, xe_, _dT = self._traverse_factored_plain(
                    z, period, times, hs, want_dT=False,
                    open_at_x0=x0_unknown)
                return (np.asarray(x0_) - np.asarray(xe_),
                        lambda v: v - alpha * self._monodromy_matvec_plain(
                            opening, steps, v))
            x0_ss, _info, _ier, _mesg = self._matrix_free_newton(
                _build, np.asarray(x, dtype=float), _tol, _tol,
                _shoot_reltol, maxiterations)
        else:
            x0_ss, _info, _ier, _mesg = analysis.fsolve(
                func, x, maxiter=maxiterations, reltol=_shoot_reltol,
                abstol=_tol, xtol=_tol, toolkit=self.toolkit,
                full_output=True, line_search=True)
        self.converged = (_ier == 1)
        self.shooting_iterations = maxiterations if not self.converged else None
        ## ⚠ AN AUTONOMOUS OSCILLATOR CANNOT BE SOLVED AT A FIXED PERIOD,
        ## and this is the only place it says so.
        ##
        ## A circuit whose oscillation is self-sustaining -- a VCO
        ## macromodel, an LC or ring oscillator, any phase accumulator
        ## driven by a DC source -- has a one-parameter family of periodic
        ## solutions, because rotating the starting point along the orbit
        ## gives another one.  Its monodromy therefore has an eigenvalue at
        ## exactly 1 and `I - M` is singular AT the true period.  Away from
        ## it the orbit does not close at all: measured on the quadrature
        ## phase element, the discretisation precesses by 2.1e-3 rad per
        ## cycle at 100 steps/period (falling as h^2), so the period map is
        ## a rotation by slightly less than 2*pi whose ONLY fixed point is
        ## the origin -- which is what an unseeded run returns, silently.
        ##
        ## Measured either side on that element: |eig(M)| = 0.968 and
        ## sigma_min(I-M) = 2.3e-02 at the nominal period, against
        ## |eig(M)| = 1.000226 and 1.6e-04 at the corrected one.  So there is
        ## no period at which this analysis both has a solution and an
        ## invertible Jacobian, and the answer is not a better seed: it is
        ## the autonomous formulation, which solves for the period jointly
        ## with a phase condition.  Not implemented -- but a run that
        ## returns the origin, or refuses to converge, deserves to be told
        ## why rather than left to look like a tolerance problem.
        rho, self.floquet_multipliers, self.parasitic_roots = \
            self._spectral_report(getattr(self, '_monodromy', None))
        self.spectral_radius = rho
        ## `self.autonomous` was decided before the solve and chose which
        ## system ran; nothing to re-derive here.  It used to WARN at this
        ## point that a self-oscillating circuit could not be solved at all,
        ## which was true of the fixed-period system and is no longer true
        ## of this one -- the period was an unknown and `self.period` holds
        ## what it came to.
        if not self.converged:
            ## ⚠ THIS USED TO BE SILENT.  `fsolve` builds the "No
            ## convergence" message and then discards it whenever
            ## `full_output=False`, which is how this call was written -- so
            ## a shooting solve that never converged returned a
            ## plausible-looking waveform with no diagnostic at all.  It was
            ## non-convergent on EVERY circuit, including a linear RLC whose
            ## answer was visibly close, which is why nobody noticed.
            ## ⚠ THIS MESSAGE USED TO CLAIM A "true Newton", AND THAT IS
            ## FALSE ON THE PLAIN PATH.  `newton` is set True
            ## unconditionally, so the alternative branch was dead and every
            ## non-convergence was reported as a Newton failing.  Measured
            ## 2026-09-02: the true `dF/dx_in` is SINGULAR on every circuit
            ## tried (rank 1/3, 2/4, 1/3; sigma_min exactly 0), so no method
            ## solves a true Newton in the frame the plain path's unknown
            ## lives in -- it is a contraction, and its residual falls
            ## LINEARLY at a constant ratio.  Advising `method='euler'` on
            ## the strength of a distinction that does not exist sent people
            ## sideways.
            ##
            ## The advice that IS backed: the solved-history route has an
            ## exact Jacobian (item 4b) and converges quadratically --
            ## measured 1.69e-01 -> 1.06e-02 -> 3.78e-06 -> 3.48e-12 against
            ## trapezoidal's linear 3.91e-03 -> 3.14e-04 -> 2.66e-05 on the
            ## same circuit.
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r). ⚠ The returned waveform IS STILL A FULL '
                'RESULT -- it is the last iterate, not a periodic steady '
                'state -- so a reader who does not check `converged` gets '
                'an array that looks like an answer and is not. Raise '
                "maxiterations, or use method='gear', whose solved-history "
                'formulation has an exact Jacobian and converges '
                'quadratically where the plain path is a contraction with a '
                'linear rate.'
                % (maxiterations, method),
                RuntimeWarning, stacklevel=2)
        
        ## THE THIRD LEVEL, MEASURED ON THE WAY OUT.
        ##
        ## ⚠ WHY THE NESTING WORKS AT ALL, which this docstring stated the
        ## shape of and never the reason for.  Kundert (*Introduction to RF
        ## Simulation*, v2 2003; relayed from the docs session, cited not
        ## verified here): "the strong convergence properties of shooting
        ## methods result from its nature AS A MULTILEVEL NEWTON METHOD, and
        ## not from the fact it is a time-domain method.  Indeed, it is
        ## possible to formulate harmonic balance as a time-domain method
        ## yet its convergence properties do not fundamentally change."
        ##
        ## The mechanism is that `phi_T` "is a near linear function ... even
        ## when the underlying circuit is behaving in a strongly nonlinear
        ## fashion, because `phi_T` is evaluated over one period of the
        ## large periodic clock signal" -- the nonlinearity is absorbed by
        ## the INNER transient, which is "a natural continuation method,
        ## quite robust".  So the outer Newton sees a nearly linear map and
        ## the arrangement below is not an arbitrary ordering of three
        ## tests: each level exists because the level inside it has already
        ## made the level outside tractable.
        ##
        ## Three convergence criteria stand between a PSS run and its answer,
        ## and the first two are checked while it runs: the inner Newton
        ## (`i(x) + iq + u` under `reltol/iabstol/vabstol`) and the shooting
        ## Newton (`x0 - phi(x0)` under the same, times `steadyratio`).  Both
        ## ask whether an EQUATION was solved.  Neither asks whether the
        ## equation was the right one -- the discrete period map is not the
        ## continuous one, and the gap between them is truncation error.
        ##
        ## PSS imposes its grid (`h = T/(N-1)`, uniform, N from `timestep`),
        ## so this cannot be a CONTROL signal -- nothing here may shrink a
        ## step, and doing so would change the period map between shooting
        ## iterations and destroy the monodromy.  It is a MEASUREMENT, taken
        ## on the converged solution over the final replay and reported.
        ##
        ## ⚠ IT IS THE LEVEL THAT WAS SILENT, and the one that dominates.
        ## On the Q=20 resonator at 100 steps/period all three integrators
        ## report a converged shooting solve, and their amplitudes are
        ## 8.815 V (euler), 19.766 V (gear2) and 19.990 V (trap) against
        ## 20 V analytic -- a 56% disagreement between two "converged"
        ## answers.  Nothing in the two Newton criteria can see that, because
        ## each integrator solved ITS OWN equations to tolerance.  This
        ## number can: it is `|J^-1 Eg| / (TRTOL (reltol ref + lte_abstol))`,
        ## the quantity a transient would have rejected a step on.
        ## ⚠ THE REPLAY MUST OPEN THE WAY THE SOLVE DID, or the waveform is
        ## not the solution: a plain replay of a solved-history answer
        ## would reintroduce the very seam that formulation exists to
        ## remove,
        ## and the reported amplitude would not be the one the residual was
        ## driven to zero on.
        ## ⚠ AND IT MUST WALK THE SAME (t, h) PAIRS, which it did not.
        ## The two traversals pair them differently: `_traverse_solved_history`
        ## walks `times[1:]` with `hs[_j]`, while the plain `_traverse` takes
        ## the MANUFACTURING step at `(times[0], hs[0])` FIRST and only then
        ## walks `times[1:]` with `hs[_j]` -- so in the plain case the step
        ## AFTER the opening one uses `hs[0]` again, not `hs[1]`.  This
        ## replay set `walk = times` and indexed `hs[min(_j, ...)]`, which
        ## pairs `times[k]` with `hs[k]` from `k = 1` on and is off by one
        ## against the traversal.
        ##
        ## A UNIFORM GRID HIDES IT COMPLETELY -- every `hs` is the same
        ## number -- which is why it survived every uniform test in this
        ## file.  Measured on the Q=20 resonator at 200 points, closure
        ## `|x(T) - x(0)|` of the RETURNED waveform:
        ##
        ##       grid      trap            gear (control)
        ##       uniform   5.61e-13        4.62e-14
        ##       4:1       1.70e-02        5.33e-15
        ##       16:1      4.88e-03        1.78e-14
        ##
        ## with `converged = True` in every row.  Gear closes on every grid
        ## because it takes the solved-history branch, whose pairing was
        ## already right; the plain path returned a waveform that is not the
        ## solution its own residual was driven to zero on.
        ##
        ## Now built as explicit `(t, h)` PAIRS rather than two sequences
        ## indexed in parallel, because the parallel indexing is the bug and
        ## a pair cannot be misaligned by one.
        ## ⚠ WHAT A LATER FACTORED REPLAY NEEDS, and the reason it is kept
        ## HERE rather than inside the Newton.  `_traverse_factored*` runs
        ## inside `_matrix_free_solve`'s `build` closure, whose `steps` go
        ## out of scope with the closure -- and the last `build` call is at
        ## the last TRIAL iterate, which is the converged one only by
        ## accident.  `PAC` wants the factors of the SOLUTION.
        ##
        ## So this keeps the four things a replay cannot re-derive from
        ## outside (which seed, which grid, which opening) and
        ## `factored_period()` runs the traversal on demand.  That costs one
        ## traversal for a caller who asks and nothing at all for one who
        ## does not -- where retaining `N` factorisations from every solve
        ## would cost `2 N m^2` doubles on every run, which is the memory
        ## trade `_traverse_factored` documents and most callers never want.
        self._period_state = (bool(solved_history), copy(x0_ss),
                              None if xm1_ss is None else copy(xm1_ss),
                              times, hs, float(period), bool(x0_unknown))
        self._factored_period_cache = None

        if solved_history:
            self._install_history(x0_ss, xm1_ss, hs[0], h_prev=hs[-1])
            tr = self._transient()
            X = [np.asarray(x0_ss, dtype=float)]
            walk = list(zip(times[1:], hs[:len(times) - 1]))
        else:
            X = [x0_ss]
            tr = self._begin_period(x0_ss)
            ## the manufacturing step, then the loop -- exactly `_traverse`
            ## ... unless there was no manufacturing step, in which case the
            ## replay opens AT `x_0` and walks the period alone.  Getting
            ## this wrong is the same class of defect as the grid shift
            ## above: a replay that does not reproduce its own traversal.
            walk = list(zip(times[1:], hs[:len(times) - 1]))
            if not x0_unknown:
                walk = [(times[0], hs[0])] + walk
        ## Fresh probe, so `relref='sigglobal'`'s running signal maximum is
        ## the period's, not something an earlier shooting iteration saw.
        tr._lte_probe = None
        ## A stage method has no LMM divided-difference LTE (compute_lte
        ## refuses), and the seam/interior split is a property of a manufactured
        ## opener it does not have -- so the replay collects no per-step LTE for
        ## it, and the three LTE figures below report None (honestly: the
        ## diagnostic does not apply to a self-starting stage method).  The
        ## method says which it is.
        self._want_lte = not self._integrator_for(method).is_stage_method()
        lte_seen = []
        for t, dt in walk:
            x = self.solve_timestep(X[-1], t, dt)
            if self._lte is not None:
                lte_seen.append((float(self._lte), float(t), self._lte_seam,
                                 self._lte_valid))
            X.append(copy(x))
        self._want_lte = False

        ## THREE NUMBERS, BECAUSE THEY HAVE DIFFERENT REMEDIES.
        ##
        ## `max_lte` is the INTERIOR per-step peak -- steps whose estimator
        ## saw only real past charges.  It is exactly the quantity a
        ## transient controls its grid on, and it ranks the integrators the
        ## way their answers rank: on the Q=20 resonator at 100 points per
        ## period it reads euler 0.2876, gear2 0.0763, trap 0.0239 against
        ## amplitudes of 8.815 / 19.766 / 19.990 V (analytic 20 V).
        ##
        ## ⚠ BUT THE PEAK IS A PER-STEP NUMBER AND A LIMIT CYCLE IS WHAT A
        ## WHOLE PERIOD DOES.  At `reltol=1e-3` euler's peak is 0.288 -- in
        ## tolerance -- while its amplitude is 56% low, because a transient's
        ## criterion bounds each step and says nothing about the 99 of them
        ## that damp the orbit.  `total_lte`, the SUM over the interior, is
        ## the one that sees it: 26.27 for that run against 0.941 for gear2
        ## and 0.340 for trap, tracking the amplitude errors of 55.9%, 1.17%
        ## and 0.05%.  It is an upper bound -- it adds magnitudes, so it
        ## cannot see the cancellation that makes trapezoidal's real error
        ## far smaller than its summed one -- which is the right direction
        ## for a diagnostic to be wrong in.
        ##
        ## `max_lte_seam` is the peak over the opening steps of a method
        ## whose COMPANION reads the entering unknown -- Gear-2 here, and
        ## `None` for euler and trapezoidal, which cannot have a seam.  See
        ## `solve_timestep` for why that is the right condition and what the
        ## looser one reported.
        ##
        ## ⚠ IT IS A FLAG, NOT A MAGNITUDE.  The seam is real -- measured at
        ## 1.266e-01 V on the Q=20 resonator at 100 points/period, against
        ## an interior contribution of 1.070e-01, and it is the term that
        ## stops converging (54% of Gear-2's error there, 73% at 400 points)
        ## -- but the NUMBER printed overstates it by orders of magnitude,
        ## because the estimator differences a fabricated charge while the
        ## solution merely reads one.  At 100 points the estimator's
        ## seam/interior ratio is 505x and the answer's is 1.18x.  So use it
        ## to know the seam is there; use `benchmarks/pss_seam_cost.py` to
        ## know what it costs.
        ##
        ## The fix is not a smaller timestep -- refining makes its SHARE
        ## grow.  It is to make the entering history part of the shooting
        ## unknowns, so the map is a fixed point in the state a two-step
        ## method actually needs, rather than one that opens off a
        ## stand-in.  Not built; the prize
        ## measured on that resonator is Gear-2's error going 2.34e-1 ->
        ## 1.07e-1 at 100 points, at no extra cost per iteration.
        ## An unsound estimate is reported as neither: for trapezoidal that
        ## step is the ONLY one whose number was ever wrong, and dropping it
        ## is what keeps it out of both figures.
        interior = [p for p in lte_seen if not p[2] and p[3]]
        seam = [p for p in lte_seen if p[2]]
        if interior:
            self.max_lte, self.max_lte_time = max(interior)[:2]
            self.total_lte = float(sum(p[0] for p in interior))
        else:                                            # pragma: no cover
            self.max_lte = self.max_lte_time = self.total_lte = None
        self.max_lte_seam = max(seam)[0] if seam else None

        ## Named so the warning can lead with whichever is actually
        ## limiting: the three have three different answers.
        _limits = [
            (self.total_lte, 'accumulated over the period',
             'use a smaller timestep or a less damping method -- this is '
             'the figure that sets a limit cycle, and a per-step criterion '
             'can be in tolerance while it is not'),
            (self.max_lte, 'in one interior step',
             'use a smaller timestep or a higher-order method'),
            (self.max_lte_seam, 'over the opening steps',
             "this is the period map's own seam, where each shooting "
             'iteration cold-starts from a fabricated history; it does '
             'NOT improve with a smaller timestep'),
        ]
        _over = [(v, where, fix) for v, where, fix in _limits
                 if v is not None and v > 1.0]
        if _over:
            v, where, fix = max(_over, key=lambda r: r[0])
            ## ⚠ DO NOT ASSERT CONVERGENCE HERE.  This clause read "the
            ## shooting solve converged, but ..." unconditionally, and the
            ## LTE report is produced whether or not it did -- so a
            ## non-converged run emitted a warning whose first words said it
            ## had converged, directly beside the warning that said it had
            ## not.  Two reviewers read non-converged waveforms as answers
            ## in this file's history; contradictory warnings are not why,
            ## but they are not help either.
            warnings.warn(
                'PSS: the shooting solve %s, and the periodic '
                'solution is not resolved at this accuracy (method=%r, %d '
                'points per period). Local truncation error reaches %.3g '
                'times tolerance %s: %s. Neither Newton criterion can see '
                'this -- they ask whether the discrete equations were '
                'solved, not whether the discretisation is the right one. '
                '(peak interior %s at t=%.6g s, period total %s, opening '
                'steps %s; relax lte_vabstol/lte_iabstol/TRTOL if this '
                'accuracy is intended.)'
                % ('converged' if self.converged else 'did NOT converge',
                   method, npts, v, where, fix,
                   'n/a' if self.max_lte is None else '%.3g' % self.max_lte,
                   -1.0 if self.max_lte_time is None else self.max_lte_time,
                   'n/a' if self.total_lte is None
                   else '%.3g' % self.total_lte,
                   'n/a' if self.max_lte_seam is None
                   else '%.3g' % self.max_lte_seam),
                RuntimeWarning, stacklevel=2)

        ## ⚠ THE PLAIN PATH'S FIRST ENTRY IS A SEED, THE OTHER PATH'S IS A
        ## SOLUTION.  Plain takes N steps from `x0_ss` and reports their
        ## results; the other starts AT `x_0` and takes N-1, so dropping the
        ## first would drop a real point and shift the waveform by a step.
        ## ⚠ AN AUTONOMOUS PERIOD IS ONLY DETERMINED UP TO AN INTEGER
        ## MULTIPLE, AND THE SOLVE FOLLOWS THE SEED.
        ##
        ## `k*T` is a period whenever `T` is, so `x0 - phi_{kT}(x0) = 0` has
        ## solutions at every multiple and the free-period system converges
        ## to whichever one the seed is nearest.  Measured on the quadrature
        ## phase element, whose true period is 1.000e-03: seeds of 1e-3,
        ## 2e-3 and 3e-3 return 1.000083e-03, 2.000665e-03 and 3.002245e-03
        ## and ALL report `converged`.  The reported waveform is a correct
        ## periodic solution in each case -- and its fundamental frequency
        ## is wrong by the factor, which is what a PSS user is usually
        ## after.  Nothing said so.
        ##
        ## The detector is cheap and needs no extra solve: an orbit
        ## traversed k times comes back near `x_0` partway through.  Grid
        ## points do not land on `T/k` in general (`T/2` at 199 steps is
        ## step 99.5), so this is a NEAREST-APPROACH test against the
        ## orbit's own diameter rather than an equality, and the endpoints
        ## are excluded because every orbit is near `x_0` there.
        ##
        ## Driven runs are exempt: their period is the caller's, and asking
        ## for two source periods is a legitimate request, not a mistake.
        if self.autonomous and len(X) > 8:
            _pts = [np.asarray(v, dtype=float) for v in X]
            _d = np.array([float(np.max(np.abs(v - _pts[0]))) for v in _pts])
            _diam = float(np.max(_d))
            ## ⚠ THE THRESHOLD IS THE GRID, NOT A FIXED FRACTION.  A k-fold
            ## orbit returns BETWEEN grid points, so the nearest approach is
            ## bounded below by how far the solution moves in one step: at
            ## 200 points over two periods that is ~1.6% of the diameter,
            ## and a fixed 1% test therefore fired on nothing.  Comparing
            ## against the per-step displacement is scale-free and tightens
            ## automatically as the grid refines.
            _step = float(np.max([np.max(np.abs(_pts[i + 1] - _pts[i]))
                                  for i in range(len(_pts) - 1)]))
            _edge = max(2, len(_d) // 20)
            _inner = _d[_edge:-_edge]
            ## ⚠ THE EARLIEST RECURRENCE, NOT THE NEAREST.  A three-fold
            ## orbit passes close to `x_0` at both `T/3` and `2T/3`, and
            ## `argmin` picked whichever happened to be numerically nearer
            ## -- it reported `2T/3` as "the fundamental", which is wrong by
            ## a factor of two and would have sent the reader to a period
            ## that is itself a multiple.
            _near = [j for j in range(_edge, len(_d) - _edge)
                     if _d[j] < 3.0 * _step and _d[j] < 0.25 * _diam]
            if _diam > 0.0 and _near:
                ## The closest approach WITHIN THE FIRST cluster: the first
                ## point over the threshold is up to a step early, which
                ## read 3% low and made the multiple look like 2.06 rather
                ## than 2.00.
                _run = [_near[0]]
                for _c in _near[1:]:
                    if _c != _run[-1] + 1:
                        break
                    _run.append(_c)
                _j = min(_run, key=lambda i: _d[i])
                if True:
                    ## `times[_j]`, not `period * _j/(N-1)`: with a caller's
                    ## grid the points are not evenly spaced.
                    self.fundamental_period = float(
                        times[min(_j, len(times) - 1)])
                    warnings.warn(
                        'PSS: this autonomous solve returned a period that '
                        'is a MULTIPLE of the fundamental. The orbit comes '
                        'back within %.2g of its own diameter at t=%.6g s, '
                        'so the fundamental is about %.6g s and the '
                        'returned %.6g s traverses it about %.1f times. '
                        'k*T solves the periodicity condition whenever T '
                        'does, so the solve follows its seed -- re-run with '
                        'period=%.6g to get the fundamental. The waveform '
                        'is a correct periodic solution either way; its '
                        'FUNDAMENTAL FREQUENCY is what is off.'
                        % (_d[_j] / _diam, self.fundamental_period,
                           self.fundamental_period, period,
                           period / self.fundamental_period,
                           self.fundamental_period),
                        RuntimeWarning, stacklevel=2)

        ## ⚠ THE FIRST ENTRY IS DROPPED ONLY WHEN IT IS NOT PART OF THE
        ## PERIOD.  On the default plain path `X[0]` is `x_in`, the
        ## pre-image of the manufactured step, which sits one step BEFORE
        ## t=0 and is not a point of the orbit.  With `x0_unknown` -- and on
        ## the solved-history path -- `X[0]` IS `x(0)`, so dropping it both
        ## discards a real sample and leaves the waveform one column short
        ## of `times`.
        X = toolkit.array(X if (solved_history or x0_unknown) else X[1:]).T

        # Insert reference node voltage
        X = toolkit.concatenate((X[:irefnode], 
                                 toolkit.zeros((1,len(times))), 
                                 X[irefnode:]))

        tpss = analysis.CircuitResult(self.cir, x=X, xdot=None,
                                      sweep_values=times, sweep_label='time', 
                                      sweep_unit='s')

        ## ⚠ KEPT FOR THE CARRIER PHASOR, which AM/PM needs and which
        ## `fpss` below cannot supply: `freq_analysis` returns an RMS,
        ## energy-folded, positive-frequency spectrum -- right for
        ## reporting and wrong for a phasor, because folding destroys the
        ## phase relationship between a carrier and its sidebands, which is
        ## the entire content of an AM/PM decomposition.
        self.waveform = (np.asarray(times, dtype=float),
                         np.asarray(X, dtype=float))

        freqs, FX = freq_analysis(X[:,:-1], times[:-1])

        ## ⚠ `fpss` IS RMS, AND THE USUAL THING TO COMPARE IT AGAINST IS NOT.
        ## `freq_analysis` returns an RMS, energy-folded, positive-frequency
        ## spectrum (see the note above).  A commercial simulator's frequency
        ## -domain PSS output is conventionally PEAK, so a harmonic read from
        ## one and compared against the other differs by `sqrt(2)` -- 3.01 dB
        ## -- with nothing in either result announcing it.  Multiply `fpss`
        ## by `sqrt(2)` for a peak-convention comparison, or divide theirs.
        ## Recorded because this campaign has lost time to factor-of-two and
        ## factor-of-pi convention defects more than once, and a 3 dB offset
        ## is small enough to be mistaken for a modelling difference.
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

class SidebandResponse(object):
    """Every input band that lands on ONE output frequency.

    ⚠ THE POINT OF THIS OBJECT IS THAT PAC'S ANSWER IS NOT ONE NUMBER.
    Kundert: *"for a single output frequency there may be many transfer
    functions from a single input"*.  Conversion gain is one entry; image
    rejection, LO feedthrough and supply rejection are others, and a
    caller who takes a single coefficient and calls it "the gain" has
    silently picked one of them.

    ⚠⚠ AND THE BANDS SIT AT DIFFERENT INPUT FREQUENCIES, WHICH IS THE
    THING THAT IS EASY TO GET WRONG.  `adjoint_sideband_row` is indexed by
    the INPUT frequency: a source at `f` reaches the output at
    `f + l f0` through sideband `l`.  Fixing the OUTPUT instead means each
    sideband is fed from its own input band, `f_in = f_out - l f0`.  So
    image rejection is NOT `H_l` against `H_-l` at one input frequency --
    it is two different input bands mapping onto one output, and computing
    it the first way gives a plausible number for a different quantity.

    Attributes: `f_out`, `f0`, `sidebands`, `inputs` (the `f_in` per
    sideband) and `rows` (`(len(sidebands), m)`, source-indexed).
    """

    __slots__ = ('f_out', 'f0', 'sidebands', 'inputs', 'rows')

    def __init__(self, f_out, f0, sidebands, inputs, rows):
        self.f_out = float(f_out)
        self.f0 = float(f0)
        self.sidebands = list(sidebands)
        self.inputs = list(inputs)
        self.rows = np.asarray(rows)

    def _index(self, l):
        try:
            return self.sidebands.index(int(l))
        except ValueError:
            raise KeyError(
                'sideband %r was not computed; this response carries %r'
                % (l, self.sidebands))

    def input_frequency(self, l):
        """The input band feeding sideband `l`: `f_out - l f0`."""
        return self.inputs[self._index(l)]

    def transfer(self, l, source):
        """The complex coefficient from `source` at `f_in(l)` to `f_out`."""
        return complex(self.rows[self._index(l)][int(source)])

    def gain_db(self, l, source):
        """`20 log10 |H_l|` -- one band's conversion gain, named as such."""
        mag = abs(self.transfer(l, source))
        if mag == 0.0:
            return -np.inf
        return 20.0 * np.log10(mag)

    def rejection_db(self, wanted, other, source):
        """How far `other` sits below `wanted`, in dB.

        ⚠ WHICH REJECTION THIS IS DEPENDS ENTIRELY ON WHICH TWO SIDEBANDS
        ARE NAMED, and the method refuses to guess.  Image rejection is
        the wanted band against the one mirrored about the LO; LO
        feedthrough is the wanted band against `l` such that `f_in = 0`.
        Naming them at the call site is the difference between a number
        and a labelled number.
        """
        w = abs(self.transfer(wanted, source))
        o = abs(self.transfer(other, source))
        if o == 0.0:
            return np.inf
        if w == 0.0:
            return -np.inf
        return 20.0 * np.log10(w / o)

    def __repr__(self):
        return ('SidebandResponse(f_out=%g, f0=%g, sidebands=%r)'
                % (self.f_out, self.f0, self.sidebands))


class ProbeShooting:
    """Bizzarri's probe-based shooting -- oscillator amplitude and frequency
    from a DRIVEN solve, plus the 2x2 power-flow instability screen.

    A periodic voltage source of amplitude ``A`` and frequency ``f`` is placed
    across ``node``, and ``(A, f)`` is solved so that the probe's OWN CURRENT
    vanishes.  At that point the probe sources nothing and can be removed
    without changing the steady state, so it is the oscillator's own solution.

    The probe makes the circuit NON-AUTONOMOUS: the period is known (``1/f``),
    so there is no phase condition, no free-period unknown, and no ``T = 0``
    trivial root for a seed below the fundamental to fall into.

    ⚠ IT IS NOT A CONVERGENCE AID, AND THE PAPER SAYS SO ITSELF.  On its own
    flagship high-Q Pierce example the authors report *"it is easy to assign a
    tentative current to the Ls inductor ... and obtain convergence in a few
    iterations (we did this with conventional SH)"*.  What this buys is the
    SWEEP: unstable limit cycles, coexisting solutions, and a stability screen.

    ⚠⚠ ONE TONE GIVES A DESCRIBING-FUNCTION SOLVE, NOT THE ORBIT.  A single
    tone forces a SINUSOID, so a non-sinusoidal orbit can only null the probe's
    FUNDAMENTAL current.  Measured on van der Pol, the frequency error is
    QUADRATIC in harmonic content -- ``df/f = 4.0 THD^2`` to 3% across two
    decades -- and the probe sits at the LC resonance at every ``mu`` because
    ``mu(u - u^3/3)`` is odd and memoryless, so its describing function shifts
    no phase.  ``harmonics=K`` forces K tones and nulls K harmonics, which is
    harmonic balance with a shooting inner solve: at K=3 the van der Pol error
    falls from 5.9e-02 to 3.8e-04.

    ⚠⚠ PROBE PLACEMENT IS CIRCUIT-SPECIFIC, AND ITS FAILURE IS NOT A SOLVER
    FAILURE.  Forcing a node fixes every state the source reaches; any state it
    does NOT reach whose DC level is then unconstrained makes the shooting
    Jacobian SINGULAR, because a whole family satisfies periodicity.  Measured
    across van der Pol's only node with no series resistance: periodicity error
    **2.11e-15** -- already a periodic solution -- reported as
    ``converged = False``.  :meth:`degenerate_placement` names that pairing.
    """

    def __init__(self, factory, node, refnode=gnd, method='gear',
                 reltol=1e-10, npts=300, maxiterations=30, phase=90.0,
                 harmonics=1, tones=None, warm_start=True):
        """`factory()` must return a FRESH circuit WITHOUT the probe.

        `tones` overrides `harmonics` with an explicit list of harmonic
        numbers (``[1, 3, 5]`` to skip the even ones).  ⚠ Do not choose it by
        eye -- see :meth:`even_harmonic_content`.
        """
        self.factory = factory
        self.node = node
        self.refnode = refnode
        self.method = method
        self.reltol = float(reltol)
        self.npts = int(npts)
        self.maxiterations = int(maxiterations)
        self.phase = float(phase)
        self.tones = ([int(k) for k in tones] if tones is not None
                      else list(range(1, int(harmonics) + 1)))
        if self.tones[0] != 1:
            raise ValueError('the fundamental (tone 1) must be present, got %r'
                             % (self.tones,))
        self.warm_start = bool(warm_start)
        self._x0 = None
        self.evaluations = 0

    ## ---- circuit construction -------------------------------------------

    def _build(self, f, amps, phases):
        """The oscillator with K probe sources IN SERIES across `node`.

        In series they share ONE physical current (KCL ties them at the
        intermediate nodes), which is the quantity to null -- so K tones need
        no new element type, only K sources and K-1 nodes.
        """
        from pycircuit.circuit.elements import VSin
        cir = self.factory()
        prev = self.node
        n = len(self.tones)
        for j, k in enumerate(self.tones):
            nxt = self.refnode if j == n - 1 else cir.add_node('__probe_n%d' % j)
            ## ⚠⚠ `vac=0` EXPLICITLY.  `VS.vac` DEFAULTS TO 1, not 0, so every
            ## probe source in the chain would be AC-excited at once -- and in
            ## series they share one current, so the PAC response came back
            ## exactly K times too large (measured ratios 1, 2, 3 at K = 1, 2,
            ## 3).  K=1 validated because there was nothing to contaminate it,
            ## which is why the diagonal alone could not catch this.
            cir['__probe%d' % j] = VSin(prev, nxt, va=float(amps[j]),
                                        freq=float(k) * float(f),
                                        phase=float(phases[j]), vac=0.0)
            prev = nxt
        return cir

    def _probe_row(self, cir):
        """The global row of the probe chain's branch current."""
        rows = cir.elementnodemap['__probe0']
        el = cir['__probe0']
        if len(el.branches) != 1:
            raise ValueError(
                'ProbeShooting expects the probe to carry exactly one branch, '
                'this one carries %d.' % len(el.branches))
        return int(rows[-1])

    ## ---- evaluation ------------------------------------------------------

    def _spectrum(self, f, amps, phases, upto=None):
        """`(I, pss)` -- the probe current's harmonic phasors.

        `upto` defaults to the configured tones; pass a larger count to look at
        harmonics the solve is NOT nulling, which is what
        :meth:`even_harmonic_content` needs.
        """
        import warnings as _w
        cir = self._build(f, amps, phases)
        row = self._probe_row(cir)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        ## ⚠ WARM START: the finite-difference columns perturb a parameter by
        ## ~1e-5, so the trajectory barely moves and solving each from cold is
        ## waste.  Measured 2.09 s cold against 1.31 s warm -- 1.60x -- with
        ## the answers agreeing to 3e-15.  A pure accelerator: it changes which
        ## iterate the Newton starts from and nothing else, and carries no
        ## assumption about the circuit.
        ## ⚠ `PSS.solve` takes the REDUCED state (length n-1), not the full
        ## one -- passing `X[:, 0]` verbatim makes a solved-history run build a
        ## 2(n-1) pair against an n-length vector and raise on the shapes.
        x0 = self._x0 if (self.warm_start and self._x0 is not None
                          and len(self._x0) == cir.n - 1) else None
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts, x0=x0,
                      maxiterations=self.maxiterations)
        self.evaluations += 1
        X = np.asarray(pss.waveform[1], dtype=float)
        if pss.converged:
            self._x0 = np.delete(X[:, 0], pss.irefnode).copy()
        ## ⚠ DROP THE DUPLICATE ENDPOINT: `waveform` carries t=0 and t=T, the
        ## same point on a periodic solution, and keeping both biases every bin.
        i_probe = X[row, :-1]
        N = i_probe.shape[0]
        n = np.arange(N)
        ks = self.tones if upto is None else list(range(1, int(upto) + 1))
        I = np.array([(2.0 / N) * np.sum(i_probe * np.exp(-2j * np.pi * k * n / N))
                      for k in ks])
        return I, pss

    def probe_current(self, A, f):
        """The FUNDAMENTAL phasor of the probe current. `I1 == 0` is the
        oscillation condition.  Returns `(I1, pss)`."""
        amps = [float(A)] + [0.0] * (len(self.tones) - 1)
        phases = [self.phase] * len(self.tones)
        I, pss = self._spectrum(f, amps, phases)
        return complex(I[0]), pss

    def even_harmonic_content(self, f, amps=None, phases=None, upto=6):
        """`|I_even| / |I_odd|` on the probe current -- the ONLY basis on which
        even tones may be dropped.

        ⚠⚠ PRUNING THE EVEN HARMONICS DOES **NOT** HOLD IN GENERAL, and the
        failure is silent: dropping a tone that is really there removes both an
        unknown and the residual row constraining it, so the solve converges to
        the wrong waveform.  Van der Pol is HALF-WAVE SYMMETRIC and its even
        content is 7.1e-16; add an even term ``beta u^2`` to the same
        nonlinearity and it is not::

            beta    H2/H1       H3/H1       H4/H1
            0.00    7.080e-16   1.168e-01   2.710e-16   <- symmetric
            0.05    2.649e-02   1.162e-01   9.239e-03
            0.20    1.058e-01   1.068e-01   3.600e-02   <- H2 EQUALS H3
            0.50    2.613e-01   5.994e-02   7.610e-02   <- H2 is 4x H3

        At ``beta = 0.5`` pruning would discard the LARGEST correction after
        the fundamental, and ``beta = 0.05`` already gives 2.6% -- there is no
        margin to judge by eye.  So this measures rather than assumes, and a
        caller passing `tones=[1, 3, 5]` should check it first.
        """
        n = len(self.tones)
        if amps is None:
            amps = [1.0] + [0.0] * (n - 1)
        if phases is None:
            phases = [self.phase] * n
        I, _pss = self._spectrum(f, amps, phases, upto=upto)
        mag = np.abs(I)
        even = mag[1::2]
        odd = mag[0::2]
        ref = float(np.max(odd)) if odd.size else 0.0
        return (float(np.max(even)) / ref if ref > 0 else float('inf')), mag

    def degenerate_placement(self, A, f, tol=1e-10):
        """`(is_degenerate, periodicity_error, converged)` for this placement.

        ⚠ A placement leaving some state's DC level unconstrained gives a
        SINGULAR shooting Jacobian: every member of a one-parameter family
        satisfies periodicity, so the solve cannot converge even though what it
        returns is already a periodic solution.  The signature is exactly that
        pairing -- a tiny periodicity error with `converged = False`.
        """
        import warnings as _w
        n = len(self.tones)
        cir = self._build(f, [float(A)] + [0.0] * (n - 1),
                          [self.phase] * n)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts,
                      maxiterations=self.maxiterations)
        X = np.asarray(pss.waveform[1], dtype=float)
        perr = float(np.max(np.abs(X[:, -1] - X[:, 0])))
        return (perr < tol and not pss.converged), perr, bool(pss.converged)

    ## ---- solves ----------------------------------------------------------

    def _pack(self, f, amps, phases):
        """Unknowns: `f`, `A_1`, then `(A_k, phi_k)` for the rest.

        `phi_1` is NOT an unknown: it is the time origin, and carrying it would
        make the system singular along the trivial time-shift direction --
        the same marginal mode a phase condition removes in an autonomous
        solve.
        """
        z = [float(f), float(amps[0])]
        for j in range(1, len(self.tones)):
            z += [float(amps[j]), float(phases[j])]
        return np.array(z, dtype=float)

    def _unpack(self, z):
        f = float(z[0])
        amps = [float(z[1])]
        phases = [self.phase]
        for j in range(1, len(self.tones)):
            amps.append(float(z[2 * j]))
            phases.append(float(z[2 * j + 1]))
        return f, amps, phases

    def solve_multitone(self, freq0, amps0, phases0=None, tol=1e-8,
                        maxiter=20, rel_step=1e-5, use_pac=False):
        """Newton driving the probe current to zero at EVERY configured tone.

        `2K` unknowns against `2K` residuals, so the system is square.  Returns
        `(f, amps, phases, info)`.

        Measured on van der Pol (`mu = 1`, autonomous `f = 0.150229`)::

            K   f          df/f         note
            1   0.159134   +5.93e-02
            2   0.159134   +5.93e-02    A_2 = 2.3e-13 -- no change at all
            3   0.150172   -3.77e-04    157x better

        ⚠ K=2 buys NOTHING here because the even harmonics do not exist on a
        half-wave-symmetric circuit -- NOT because two tones cannot help.  See
        :meth:`even_harmonic_content` before concluding the same elsewhere.
        """
        n = len(self.tones)
        if phases0 is None:
            phases0 = [self.phase] * n
        z = self._pack(freq0, amps0, phases0)
        hist = []
        r = None
        for it in range(int(maxiter)):
            f, amps, phases = self._unpack(z)
            I, _p = self._spectrum(f, amps, phases)
            r = np.concatenate([[c.real, c.imag] for c in I])
            hist.append((f, list(amps), float(np.linalg.norm(r))))
            if np.linalg.norm(r) < tol:
                return f, amps, phases, {
                    'iterations': it, 'residual': float(np.linalg.norm(r)),
                    'history': hist, 'converged': True,
                    'evaluations': self.evaluations}
            J = np.zeros((2 * n, z.shape[0]))
            if use_pac:
                ## ⚠ The VOLTAGE block from K linear PAC solves; the FREQUENCY
                ## column stays a finite difference, because `df` moves every
                ## tone at once and is not a small-signal excitation about the
                ## same operating point.  Validated on the FIRST iteration only
                ## -- the conventions do not change as the Newton walks, and
                ## paying the check every iteration would give back the saving.
                Jv, vinfo = self.pac_jacobian(f, amps, phases,
                                              validate=(it == 0))
                ## CHAIN RULE, NOT A POSITIONAL COPY.  `pac_jacobian` returns
                ## `d(Re I, Im I)/d(Re V, Im V)`; the unknowns here are
                ## `(A_k, phi_k)` with `phi` in DEGREES.  Copying the columns
                ## across positionally feeds the Newton a Jacobian for the
                ## WRONG VARIABLES -- it diverged to f = 0.0348 against 0.1502
                ## while `pac_jacobian` itself validated at 1e-04, which is how
                ## a correct derivative and a broken solve coexisted.
                ## With `V_k = A_k exp(j psi_k)`, `psi_k = (phi_k - phase0)*pi/180`:
                ##     dV/dA   = (cos psi, sin psi)
                ##     dV/dphi = A (pi/180) (-sin psi, cos psi)
                dzf = rel_step * max(abs(z[0]), 1e-3)
                zz = z.copy()
                zz[0] += dzf
                f2, a2, p2 = self._unpack(zz)
                I2, _ = self._spectrum(f2, a2, p2)
                r2 = np.concatenate([[c.real, c.imag] for c in I2])
                J[:, 0] = (r2 - r) / dzf
                rad = np.pi / 180.0
                for jx in range(n):
                    psi = (phases[jx] - self.phase) * rad
                    cA, sA = np.cos(psi), np.sin(psi)
                    col_re = Jv[:, 2 * jx]
                    col_im = Jv[:, 2 * jx + 1]
                    dA_col = cA * col_re + sA * col_im
                    dP_col = amps[jx] * rad * (-sA * col_re + cA * col_im)
                    if jx == 0:
                        J[:, 1] = dA_col          # phi_1 is the time origin
                    else:
                        J[:, 2 * jx] = dA_col
                        J[:, 2 * jx + 1] = dP_col
            else:
                for j in range(z.shape[0]):
                    dz = rel_step * max(abs(z[j]), 1e-3)
                    zz = z.copy()
                    zz[j] += dz
                    f2, a2, p2 = self._unpack(zz)
                    I2, _ = self._spectrum(f2, a2, p2)
                    r2 = np.concatenate([[c.real, c.imag] for c in I2])
                    J[:, j] = (r2 - r) / dz
            step, _res, _rank, _sv = np.linalg.lstsq(J, -r, rcond=None)
            z = z + step
        f, amps, phases = self._unpack(z)
        return f, amps, phases, {
            'iterations': maxiter,
            'residual': float(np.linalg.norm(r)) if r is not None else None,
            'history': hist, 'converged': False,
            'evaluations': self.evaluations}

    def solve(self, amp0, freq0, tol=1e-9, maxiter=20, damp=1.0,
              rel_step=1e-4):
        """Single-tone Newton on `(A, f)`.  Returns `(A, f, info)`.

        Kept as its own entry point because the one-tone case is the cheap
        screen -- three solves an iteration -- and because its return shape is
        two scalars rather than the vectors :meth:`solve_multitone` returns.
        """
        A, f = float(amp0), float(freq0)
        hist = []
        r = np.array([np.inf, np.inf])
        for it in range(int(maxiter)):
            I0, _p = self.probe_current(A, f)
            r = np.array([I0.real, I0.imag], dtype=float)
            hist.append((A, f, float(np.linalg.norm(r))))
            if np.linalg.norm(r) < tol:
                return A, f, {'iterations': it,
                              'residual': float(np.linalg.norm(r)),
                              'history': hist, 'converged': True,
                              'evaluations': self.evaluations}
            dA = rel_step * max(abs(A), 1e-12)
            df = rel_step * max(abs(f), 1e-12)
            IA, _ = self.probe_current(A + dA, f)
            IF, _ = self.probe_current(A, f + df)
            J = np.array([[(IA.real - I0.real) / dA, (IF.real - I0.real) / df],
                          [(IA.imag - I0.imag) / dA, (IF.imag - I0.imag) / df]])
            if abs(np.linalg.det(J)) < 1e-300:
                raise np.linalg.LinAlgError(
                    'ProbeShooting: the 2x2 probe Jacobian is singular at '
                    'A=%.6g f=%.6g. Either the probe cannot see the '
                    'oscillation, or the placement leaves a state undetermined '
                    '-- check `degenerate_placement`.' % (A, f))
            step = np.linalg.solve(J, -r)
            A += damp * step[0]
            f += damp * step[1]
        return A, f, {'iterations': maxiter,
                      'residual': float(np.linalg.norm(r)),
                      'history': hist, 'converged': False,
                      'evaluations': self.evaluations}

    #: excitation offset as a fraction of `f`, so the folded sideband pair
    #: lands at distinguishable frequencies -- see `_pac_response`.
    _pac_delta = 1e-4

    def _pac_response(self, pss, cir, row, f, excite_harmonic, want_harmonics):
        """`{m: dI_m}` from ONE linear PAC solve exciting harmonic `j`.

        ⚠⚠ TWO CONVENTION FACTORS AND ONE INDEXING TRAP, all three pinned
        against a circuit whose answer is analytic (a resistor across the
        probe, where `dI/dV = 1/R` exactly) rather than against the finite
        difference this is meant to replace.  Calibrating against FD would
        make the agreement circular and would absorb a genuine sideband-index
        error into the fitted constant.

        * **The frequency list carries DUPLICATES.**  `0.159155` appears twice
          in the returned sweep, one entry near zero and one carrying the
          response; `argmin(|fs - target|)` picks whichever comes first and it
          was the wrong one, reading 5.9e-21 where the answer is 1e-3.  So the
          entry is chosen by LARGEST RESPONSE among those at the target
          frequency, not by proximity alone.
        * **A factor of two**: this file's probe spectrum uses the peak-amplitude
          convention `(2/N) sum(...)`; PAC returns a phasor.
        * **A 90 degree rotation that is OURS, not PAC's**: `VS` builds its AC
          phasor as `vac * exp(j*phase)`, and the operating-point probe sets
          `phase = 90` to make the forcing a cosine -- so the same parameter
          rotates the AC excitation.  Dividing by the excitation phasor removes
          it, which is why the excitation is read from the circuit rather than
          assumed to be 1.
        """
        tk = cir.toolkit
        pac = PAC(cir, toolkit=tk)
        ## ⚠⚠⚠ EXCITE SLIGHTLY OFF THE HARMONIC, WHICH REMOVES THE AMBIGUITY
        ## INSTEAD OF GUESSING IT.  Exciting exactly at `j*f0` sends TWO
        ## sidebands to the same absolute output frequency -- `k = m - j` and
        ## `k = -m - j` -- and `PAC.solve` returns absolute frequencies with the
        ## sideband index folded away, so the two arrive in an order that is not
        ## stable.  Ordering them by magnitude worked AT THE SOLUTION and failed
        ## away from it (validation 1.6e-05 at the solved amplitudes, 1.763 at
        ## the Newton's starting point), which is the kind of heuristic that
        ## passes a gate and then fails in use.
        ##
        ## With the excitation at `j*f0 + delta` every output lands at
        ## `(j+k)*f0 + delta`, all distinct, so the wanted term is simply the
        ## one nearest `m*f0 + delta` and no ordering rule is needed.  `delta`
        ## is small enough to be a negligible perturbation of the response and
        ## large enough to separate the pair.
        delta = self._pac_delta * float(f)
        f_ex = float(excite_harmonic) * float(f) + delta
        res = pac.solve(pss, freqs=np.array([f_ex]))
        fs = np.asarray(res.sweep_values, dtype=float)
        X = np.asarray(res.x)
        ## the excitation phasor actually applied, read from the circuit
        (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), pss.irefnode,
                                 tk)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        scale = u_ac[np.argmax(np.abs(u_ac))]
        out = {}
        for m in want_harmonics:
            ## ⚠⚠ BOTH TERMS ARE PHYSICAL, AND THE OFFSET IS WHAT NAMES THEM.
            ## Exciting at `j*f0 + delta`, the DIRECT sideband `k = m - j`
            ## lands at `m*f0 + delta`, and the IMAGE `k = -m - j` lands at
            ## `-m*f0 + delta`, which PAC folds onto `m*f0 - delta` and
            ## conjugates.  Taking only the direct one halves the answer
            ## (measured: a uniform 0.5 at every K); taking both without being
            ## able to tell them apart is what forced the magnitude-ordering
            ## heuristic that passed at the solution and failed at the Newton's
            ## start.  With the offset they are identified by FREQUENCY, so the
            ## rule is derived rather than guessed.
            base = float(m) * float(f)
            i_dir = np.where(np.abs(fs - (base + delta)) < 0.25 * delta)[0]
            i_img = np.where(np.abs(fs - (base - delta)) < 0.25 * delta)[0]
            if i_dir.size == 0:
                continue
            d_val = complex(X[row, i_dir[int(np.argmin(np.abs(
                fs[i_dir] - (base + delta))))]])
            g_val = (complex(X[row, i_img[int(np.argmin(np.abs(
                fs[i_img] - (base - delta))))]]) if i_img.size else 0.0 + 0.0j)
            out[m] = -(d_val - g_val) / scale
            continue
            near = np.where(np.abs(fs - (base + delta)) < 0.25 * delta)[0]
            if near.size == 0:
                continue
            ## ⚠⚠ THE DUPLICATES ARE A PAIR AND THEY SUBTRACT, NOT ADD.
            ## PAC folds negative output frequencies onto |f| and conjugates
            ## them, so two entries land at each harmonic.  Measured, at the
            ## fundamental:
            ##
            ##     resistor       e1 = -3.4e-21   e2 = +1.0e-03
            ##     van der Pol    e1 = -9.900250e-01   e2 = +9.900507e-01
            ##
            ## They are nearly EQUAL AND OPPOSITE on a circuit with harmonic
            ## content, so SUMMING them cancels (2.6e-05 against a true
            ## 1.98) and taking the LARGEST halves the answer -- which is
            ## exactly the factor 2 that appeared on van der Pol and not on
            ## the resistor, where one member is ~0 and max == sum == diff.
            ## The difference reproduces the finite difference on BOTH
            ## fixtures at once, which is the falsifier a fitted constant
            ## could not have passed.
            ## ⚠⚠⚠ ORDER THE PAIR BY MAGNITUDE, NOT BY ARRAY POSITION.
            ## `PAC.solve` returns ABSOLUTE frequencies -- the sideband index
            ## is folded away -- so the two entries at a harmonic arrive in an
            ## order that is NOT stable across `m`.  Measured, exciting
            ## harmonic 1 at K=3:
            ##
            ##     m=1   [ +0.75j , 3.94e-04 - 1.041251j ]   larger is 2nd
            ##     m=3   [ -0.75j , -0.25j             ]   larger is 1st
            ##
            ## Taking `e[-1] - e[0]` therefore gave m=1 correctly and m=3 with
            ## the RIGHT MAGNITUDE AND THE WRONG SIGN (ratio -0.999985).  The
            ## larger entry is the DIRECT response and the smaller its image,
            ## which is an ordering the array position does not carry.
            ## ⚠ Fragile where the two magnitudes are close; the `validate`
            ## gate is what stands behind it.
            ## exactly one entry now -- the offset separated the pair
            best = near[int(np.argmin(np.abs(fs[near] - target)))]
            out[m] = -complex(X[row, best]) / scale
        return out

    def pac_jacobian(self, f, amps, phases, validate=True, rel_step=1e-5,
                     rtol=5e-2):
        """`dI/dV` from K LINEAR PAC solves instead of 2K nonlinear ones.

        The finite-difference Jacobian recomputes, by `2K` full PSS solves, a
        quantity that is a PERIODIC SMALL-SIGNAL response about the operating
        point the base solve already produced.  One PSS plus K PAC solves gets
        the voltage block, which is where the cost is.

        Returns `(J, info)` with `J` the real `2K x 2K` block
        `d(Re I, Im I)/d(Re V, Im V)`.

        ⚠⚠⚠ **NOT SHIPPED-READY: THE NORMALISATION IS INCOMPLETE, AND THE
        DEFAULT VALIDATION CORRECTLY REFUSES.**  The route is confirmed viable
        -- PAC returns the right quantity, verified against a resistor where
        `dI/dV = 1/R` analytically -- and two of the three discrepancies are
        pinned and removed.  A third is not:

            fixture                 PAC / FD after normalisation
            resistor (linear)       -1.0        (sign only)
            van der Pol (K=1)       -0.499999   (sign AND a factor 2)

        **A constant that differs between two circuits is not a convention,
        it is a missing term**, so the remaining factor is NOT applied by
        fitting it -- that would make the "independent" Jacobian a fit to the
        finite difference it replaces, hide any sideband-index error inside the
        fitted constant, and reproduce exactly the circular-verification
        failure this file already records for the PPV.

        Until it is derived, `validate=True` raises on real circuits and the
        finite-difference Jacobian in :meth:`solve_multitone` remains the
        shipped path.  What IS established and reusable:

        * the response is present and correct (1/R recovered exactly);
        * the 90 degree rotation is OURS -- `VS` builds its AC phasor as
          `vac * exp(j*phase)` and the operating-point probe sets `phase = 90`
          -- and dividing by the excitation phasor removes it, measured: the
          residual ratio is real, not imaginary;
        * PAC's frequency sweep contains DUPLICATE entries at the same
          frequency, one near zero and one carrying the response, so
          `argmin(|fs - target|)` reads 5.9e-21 where the answer is 1e-3.

        ⚠⚠ `validate=True` CHECKS ONE COLUMN AGAINST THE FINITE DIFFERENCE AND
        RAISES ON DISAGREEMENT, and it is on by default deliberately.  A wrong
        Jacobian does not announce itself: the Newton still converges, to the
        wrong orbit -- the same silent failure that even-harmonic pruning
        produces, which this file already has a falsifier for.  The check costs
        ONE extra solve, amortised over the whole solve, and it exercises the
        conventions in `_pac_response` on the circuit actually in hand rather
        than on the resistor they were derived from.

        ⚠ The OFF-DIAGONAL entries `dI_m/dV_j` with `m != j` are the ones that
        exercise the sideband map `k = m - j` and PAC's negative-frequency
        conjugation.  With `harmonics=1` there are none, so a passing K=1
        validation says NOTHING about the index map -- validate at K >= 2
        before trusting a multitone Jacobian.
        """
        import warnings as _w
        n = len(self.tones)
        ## ⚠⚠ ONE PSS SOLVE FOR EVERY COLUMN.  `vac` is read ONLY under
        ## `analysis='ac'` -- it does not enter the transient residual, so it
        ## cannot move the periodic operating point.  Building a fresh circuit
        ## and re-solving the PSS per column therefore recomputed the SAME
        ## orbit K times, which is the whole cost this method exists to avoid:
        ## it made the "cheap" Jacobian K nonlinear solves plus K linear ones,
        ## against the finite difference's 2K.  Solve once, then walk `vac`
        ## across the probes and take K LINEAR PAC solves against that one
        ## operating point.
        cir = self._build(f, amps, phases)
        row = self._probe_row(cir)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts,
                      maxiterations=self.maxiterations)
        if not pss.converged:
            raise RuntimeError(
                'pac_jacobian: the periodic operating point did not converge, '
                'so there is nothing to linearise about.')
        J = np.zeros((2 * n, 2 * n))
        cols = {}
        for jx, k in enumerate(self.tones):
            ## excite exactly one probe; the others must be silent, and `vac`
            ## DEFAULTS TO 1 on a VS, which is why every one is set explicitly
            for jj in range(n):
                cir['__probe%d' % jj].ipar.vac = 1.0 if jj == jx else 0.0
            cir.update_iparv()
            resp = self._pac_response(pss, cir, row, f, k, self.tones)
            cols[jx] = resp
            for ix, m in enumerate(self.tones):
                d = resp.get(m, 0.0 + 0.0j)
                J[2 * ix, 2 * jx] = d.real
                J[2 * ix + 1, 2 * jx] = d.imag
                J[2 * ix, 2 * jx + 1] = -d.imag
                J[2 * ix + 1, 2 * jx + 1] = d.real
        info = {'columns': cols, 'validated': False}
        if validate:
            I0, _ = self._spectrum(f, amps, phases), None
            I0 = I0[0]
            dA = rel_step * max(abs(amps[0]), 1e-12)
            a2 = list(amps)
            a2[0] += dA
            I1, _p = self._spectrum(f, a2, phases)
            fd = (I1 - I0) / dA
            pac_col = np.array([J[2 * ix, 0] + 1j * J[2 * ix + 1, 0]
                                for ix in range(n)])
            num = np.linalg.norm(pac_col - fd)
            den = max(np.linalg.norm(fd), 1e-300)
            info['validation_reldiff'] = float(num / den)
            info['validated'] = bool(num / den < rtol)
            if not info['validated']:
                raise ValueError(
                    'pac_jacobian: the PAC column disagrees with the finite '
                    'difference by %.3e (relative), above rtol=%.2e. PAC gave '
                    '%s against FD %s. Do NOT use this Jacobian -- a wrong one '
                    'converges silently to the wrong orbit.'
                    % (num / den, rtol,
                       np.array2string(pac_col, precision=5),
                       np.array2string(fd, precision=5)))
        return J, info

    def power_flow(self, A, f, rel_step=1e-4):
        """The 2x2 power-flow screen `P = dv_R dy_R + dv_I dy_I`.

        ⚠⚠ ONE-DIRECTIONAL, AND THAT IS THE AUTHORS' OWN STATEMENT.  Only
        `P > 0 => unstable` is proven; they say the converse and
        `P < 0 => stable` "have not been proven", only tested.  So this is a
        cheap INSTABILITY DETECTOR and never a replacement for
        `_spectral_report`: a non-positive `P` means NOT DETECTED, not stable.

        Returns `(P_max, info)`, `P_max` being the largest value over
        perturbation directions -- the largest eigenvalue of the symmetric part
        of `dY/dV`.  It skips the system Jacobian's eigenvalues, which is the
        whole point of the construction.
        """
        I0, _ = self.probe_current(A, f)
        Y0 = I0 / A
        dv = rel_step * max(abs(A), 1e-12)
        Ir, _ = self.probe_current(A + dv, f)
        Yr = Ir / (A + dv)
        ph = dv / max(abs(A), 1e-12)
        Ii, _ = self.probe_current(A, f * (1.0 + ph / (2.0 * np.pi)))
        Yi = Ii / A
        dYdV = np.array([[(Yr.real - Y0.real) / dv, (Yi.real - Y0.real) / dv],
                         [(Yr.imag - Y0.imag) / dv, (Yi.imag - Y0.imag) / dv]])
        sym = 0.5 * (dYdV + dYdV.T)
        eig = np.linalg.eigvalsh(sym)
        return float(np.max(eig)), {'dYdV': dYdV, 'symmetric_part': sym,
                                    'eigenvalues': eig, 'Y0': Y0,
                                    'unstable': bool(np.max(eig) > 0.0)}


class PAC(Analysis):
    """Small-signal analysis over a periodic operating point, matrix-free.

    THE OPERATOR IS THE MONODROMY, and the whole method is one line of
    algebra on the withdrawn implementation's own system.  That system was

        (L + alpha B) v = -u,   alpha = exp(-2j pi f T)

    with `L` the block lower bidiagonal discretisation over the period and
    `B` the periodic wrap.  Telichevesky, Kundert & White (DAC 1996) reach
    the iterative form by "reinterpreting the use of `L^-1` ... as a
    preconditioner":

        (I + alpha L^-1 B) v = -L^-1 u

    `L` is block lower bidiagonal, so applying `L^-1` is forward
    substitution through the timesteps -- which is the recursion PSS already
    runs against stored factors -- and `B` is confined to the first `m` rows
    and last `m` columns, so `L^-1 B` acts only on the LAST block.  Both
    claims are checked against our own matrices in
    `test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`.

    What is left after that is `m x m`:

        (I - alpha M) y_0 = alpha w(f)

    with `M` the monodromy and `w` the forced response over one period from
    a zero initial state.  `y_0` is the small-signal state at `t = 0`; one
    more driven replay gives the rest of the period.

    ⚠ WHY THE 419.5 GiB IS GONE, precisely.  It was never the operator: it
    was the cost of FORMING `L` and `B`, `(N m)^2` complex entries, 279.7 +
    139.8 GiB at `N = 137`, `m = 1000`.  Nothing here forms either.  The
    stored per-step factors PSS already makes are the preconditioner, and
    the only dense object is `m x m` and only if the caller asks for it.

    ⚠ AND THE OLD `L` WAS BACKWARD-EULER-SHAPED, which is the trap a rewrite
    falls into.  It has two terms per row; a two-step method's variational
    system has three.  Rebuilding it for `trap` or `gear` gives an operator
    for a different recursion than the trajectory it came from -- measured,
    spectral radius 0 against the analytic 0.8546
    (`test_the_pac_L_is_backward_euler_only`).  Taking `M` from the
    traversal cannot make that mistake, because every step carries its own
    `(alphas, b)`.
    """

    parameters  = [Parameter(name='analysis', desc='Analysis name',
                             default='ac')]

    ## How hard GMRES is asked to solve the `m x m` system, relative to the
    ## PSS reltol that produced the operating point.  Looser than the
    ## operating point itself would be answering a question the trajectory
    ## cannot support; much tighter buys nothing, because the linearisation
    ## is only as good as the trajectory.
    KRYLOV_FACTOR = 1e-2

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(PAC, self).parameters + self.parameters
        super(PAC, self).__init__(cir, toolkit=toolkit, **kvargs)

    def solve(self, pss, freqs, refnode=gnd, recycle=True):
        """Sideband response at each frequency in `freqs`.

        `pss` must be a CONVERGED `PSS` -- the periodic operating point is
        what this linearises about, and there is no meaningful small-signal
        answer about a non-solution.  `PSS.factored_period()` enforces it.

        `recycle` shares one Krylov subspace across the sweep, which is
        where the sweep's cost goes; see `_solve_subspace`.
        """
        toolkit = self.toolkit
        freqs = np.atleast_1d(np.asarray(freqs, dtype=float))
        fp = pss.factored_period()
        T = float(fp.T)
        m = self.cir.n - 1

        irefnode = self.cir.get_node_index(refnode)
        if irefnode != pss.irefnode:
            raise ValueError(
                'PAC: refnode (index %d) differs from the PSS the operating '
                'point came from (index %d). The monodromy eliminated one '
                'and this would report against the other.'
                % (irefnode, pss.irefnode))
        ## ⚠ `analysis=` BY KEYWORD.  `Circuit.u(t, epar, analysis, ...)`
        ## takes `epar` second, and the withdrawn body wrote
        ## `self.cir.u(0, analysis_name)` -- passing 'ac' as the element
        ## parameter set and taking the TRANSIENT source vector, which is
        ## zero at `t = 0` for every sinusoid.  The whole analysis would
        ## have returned zeros, silently, with no source to speak of.
        (u_ac,) = remove_row_col((self.cir.u(0, analysis=self.par.analysis),),
                                 irefnode, toolkit)
        if not np.any(np.asarray(u_ac)):
            raise ValueError(
                'PAC: the %r source vector is identically zero, so there is '
                'nothing to analyse. Independent sources take their '
                'small-signal amplitude from `vac`/`iac`, not from `va` -- '
                'a source with va= set and vac=0 drives the operating point '
                'and not this.' % self.par.analysis)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()

        ## the forced response at each frequency -- one period replay each,
        ## and unavoidable: the source is what changes across the sweep
        ## ⚠ THE MANUFACTURING STEP IS NOT IN `steps`, AND IT COSTS AN
        ## ORDER.  On the plain path `_traverse_factored_plain` takes one
        ## step OUTSIDE the loop to manufacture a history, and folds its
        ## effect into the `opening` triple as a flat-history assumption.
        ## For the HOMOGENEOUS map that is the documented approximation the
        ## whole plain path is built on.  For the DRIVEN one it also means
        ## the source is never applied at that step -- one step of `u` out
        ## of `N`, i.e. a relative O(h).
        ##
        ## ⚠ MEASURED, on the Q=20 resonator against the AC analysis at
        ## 700 Hz, rel error per doubling of the grid:
        ##
        ##     trap, plain            2.00x  (O(h))   4.13e-03 at 250 pts
        ##     trap, x0_unknown=True  4.00x  (O(h^2)) 1.09e-04 at 250 pts
        ##     euler, either          2.00x  (O(h))   1.40e-02, unchanged
        ##
        ## The euler row is the control: `x0_unknown` does not move it at
        ## all (identical to five digits), so the trapezoidal gain is the
        ## manufacturing step and not something else the formulation does.
        ## The trajectory is NOT the problem -- trap's waveform converges at
        ## 4.2x per doubling either way.
        ##
        ## So this is a silent order loss for a caller who did nothing
        ## wrong, which is the one thing worth a warning.  Gear-2 takes the
        ## solved-history path and has no manufacturing step at all.
        if (fp.kind == 'plain' and not fp.open_at_x0
                and pss.par.method != 'euler'):
            warnings.warn(
                'PAC: this operating point was solved on the PLAIN path '
                'with a manufacturing step (method=%r, x0_unknown=False). '
                'The manufacturing step carries no small-signal source, so '
                'the response is FIRST order in the timestep whatever the '
                "method's own order -- measured 2.00x per doubling against "
                '4.00x for the same run with x0_unknown=True. The answer is '
                'not wrong, it is one order less accurate than the '
                'trajectory it came from. Re-solve with x0_unknown=True, or '
                "with method='gear', to get the method's own order."
                % pss.par.method,
                RuntimeWarning, stacklevel=2)

        self._check_circuit(pss)
        for f in freqs:
            self._check_harmonic(pss, f, 'a sweep point')

        rhs = []
        for f in freqs:
            w, _ = pss._forced_replay(fp, f, u_ac)
            rhs.append(np.exp(-2j * np.pi * f * T) * np.asarray(w))

        alphas = [np.exp(-2j * np.pi * f * T) for f in freqs]
        tol = max(pss.par.reltol * self.KRYLOV_FACTOR, 1e-14)
        ## ⚠ ON AN OSCILLATOR THE OPERATOR HAS THE ANSWER'S OWN POLE at every
        ## harmonic (see `_check_harmonic`), and a plain solve near one
        ## carries relative error `eta / (2 pi df/f0)`, `eta = |lambda_1 - 1|`
        ## the computed unit multiplier's displacement -- measured to four
        ## digits over five decades (docs session, Gourary reading).  The
        ## deflated route (`_deflated_solve`) borders the pole out and is
        ## exact there; it was wired into `adjoint_sideband_row` only, and
        ## this sweep solved plain outside HARMONIC_GUARD (2026-09-08).
        ## Under the radau default eta ~ 1e-12 puts the unguarded band
        ## inside the guard, so this is correctness hygiene, not a fix a
        ## user would see; the subspace recycling across frequencies is
        ## given up on the autonomous path (one bordered solve per point).
        self.deflated = bool(getattr(pss, 'autonomous', False))
        if self.deflated:
            ys = [self._deflated_solve(pss, a, b, transposed=False, tol=tol)
                  for a, b in zip(alphas, rhs)]
            self.matvecs = None
        elif recycle:
            ys, self.matvecs = self._solve_subspace(fp, alphas, rhs, tol)
        else:
            ys, self.matvecs = self._solve_each(fp, alphas, rhs, tol)

        ## one driven replay per frequency turns `y_0` into the period
        outfreq, outV = [], []
        for f, y0 in zip(freqs, ys):
            _end, ysteps = pss._forced_replay(fp, f, u_ac, y0=y0, collect=True)
            y = np.array([np.asarray(y0)[:m]] + [np.asarray(v)[:m]
                                                 for v in ysteps])
            ## `v(t) = y(t) exp(-j w t)` is T-periodic; its DFT is the
            ## sideband set, exactly as the withdrawn body intended
            tms = np.asarray(fp.times, dtype=float)[:len(y)]
            v = y * np.exp(-2j * np.pi * f * tms)[:, None]
            ## ⚠ TWO REPORTING DEFECTS, FOUND BY AN EXTERNAL REFERENCE CROSS-CHECK
            ## (2026-09-05), neither in the solve.  (a) `fp.times` spans
            ## `[0, T]` INCLUSIVE, so the last sample repeats the first on a
            ## T-periodic `v` (|v[0] - v[-1]| / |v[0]| = 7e-18 measured) and
            ## the DFT's `dt = T/(N-1)` put the sidebands at `f0 (N-1)/N`:
            ## 99 500 Hz for 100 000 at N = 200 -- and cost an ORDER, O(h)
            ## for O(h^2), 68x at 800 points.  `PSS.solve` already drops
            ## the endpoint one function away; this did not.  Guarded on
            ## the window rather than sliced blind, since the plain path's
            ## `[:len(y)]` need not be inclusive.  (b) `|sb + f|` folded a
            ## NEGATIVE sideband frequency to positive and left the
            ## coefficient alone; the physical response there is the
            ## CONJUGATE.  Uncorrected, `l = -1` was 166% off and did not
            ## converge under refinement; conjugated it lands on its
            ## positive twin's error to three digits (4.873e-3 / 4.877e-3).
            ## Both defects are invisible on a circuit whose `v(t)` is
            ## constant over the period -- every earlier PAC gate.
            if len(tms) > 1 and np.isclose(tms[-1] - tms[0], T,
                                           rtol=1e-9, atol=0.0):
                v, tms = v[:-1], tms[:-1]
            sb, V = freq_analysis(v, tms, axis=0)
            fs = np.asarray(sb, dtype=float) + f
            V = np.asarray(V)
            neg = fs < 0.0
            if np.any(neg):
                V = V.copy()
                V[neg] = np.conj(V[neg])
            outfreq.extend(np.abs(fs).tolist())
            outV.extend(V.tolist())

        order = np.argsort(np.asarray(outfreq))
        fout = np.asarray(outfreq)[order]
        X = np.asarray(outV)[order]
        X = np.concatenate((X[:, :irefnode],
                            np.zeros((len(fout), 1)),
                            X[:, irefnode:]), axis=1)
        self.result = analysis.CircuitResult(
            self.cir, x=X.T, xdot=None, sweep_values=fout,
            sweep_label='freq', sweep_unit='Hz')
        return self.result

    def adjoint_transfer_row(self, pss, freq, output, recycle_tol=None):
        """Every source to ONE output, in a single transposed solve.

        Returns a row `r` of length `m`: `r[i]` is the small-signal
        response at `output` (at `t = 0`) to a unit source injected at
        reduced coordinate `i` at `freq`. Forward, that is `m` separate
        solves; here it is one.

        ⚠ THIS IS THE ASYMMETRY pnoise IS SHAPED BY, and the reason
        Okumura et al. (1993) reach for the adjoint at all: "it is
        efficient to use the adjoint method ... BECAUSE CIRCUITS HAVE MANY
        NOISE SOURCES." Recycling does not help the forward route, because
        the right-hand side is what changes from source to source.

            output = d^T y_0 = alpha * d^T (I - alpha M)^-1 w(u)
                             = alpha * ((I - alpha M)^-T d)^T W u

        so one transposed solve for `x^a`, then `W^T x^a` from the reverse
        replay, and the whole row falls out. MEASURED against `m` forward
        solves on an RC ladder: agreement 9.6e-16.

        ⚠ WHAT THIS IS NOT, so the next reader does not over-read it. The
        output here is the state at `t = 0`, a single linear functional.
        A SIDEBAND coefficient `H_l` is a functional DISTRIBUTED over the
        period -- `(1/N) sum_n exp(-j l w0 t_n) d^T y_n` -- and its adjoint
        needs the reverse pass to take an injection at every step rather
        than a seed at the end. That extension is the next piece of A3, and
        it is not built.

        ⚠ WAS SOLVED-HISTORY ONLY until B8 gave the one-step companions
        their own reverse recursion; it now runs under every method.
        """
        import scipy.sparse.linalg as spla
        ## ⚠ NO LONGER GEAR-ONLY (B8): the transposed replay exists for the
        ## one-step companions too, and every use below goes through
        ## `fp.matvec_transposed`.
        fp = pss.factored_period()
        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the adjoint row')
        m = pss.cir.n - 1
        n = fp.width
        alpha = np.exp(-2j * np.pi * float(freq) * float(fp.T))

        d = np.zeros(n, dtype=complex)
        if np.isscalar(output):
            d[int(output)] = 1.0
        else:
            out = np.asarray(output, dtype=complex).ravel()
            d[:len(out)] = out

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        tol = (self.KRYLOV_FACTOR * pss.par.reltol if recycle_tol is None
               else recycle_tol)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        ## the same pole as in `solve` and `adjoint_sideband_row`: deflated
        ## on an oscillator, plain (and cheaper) on a driven circuit
        self.deflated = bool(getattr(pss, 'autonomous', False))
        if self.deflated:
            xa = self._deflated_solve(pss, alpha, d, transposed=True,
                                      tol=max(tol, 1e-14))
        else:
            xa = self._gmres_checked(A, d, max(tol, 1e-14), 'the adjoint solve')
        self.matvecs = count[0]
        return alpha * pss._forced_replay_transposed(fp, freq, xa)

    def adjoint_sideband_row(self, pss, freq, output, sidebands=0):
        """`H_l` rows: every source to ONE output's sideband `l`.

        Returns an array of shape `(len(sidebands), m)`. Entry `[li, i]` is
        the coefficient at sideband `l` of the output at `output`, for a
        unit source injected at reduced coordinate `i` at `freq`:

            H_l = (1/N) sum_n exp(-j l w0 t_n) d^T y_n

        ⚠ THIS IS THE ONE `adjoint_transfer_row` IS NOT.  That row is the
        response at a single instant -- one linear functional, adjointed by
        seeding the reverse pass at the end.  A SIDEBAND is a functional
        DISTRIBUTED over the period, so its adjoint takes an injection at
        EVERY step, and the answer comes in two pieces:

            dH/du  =  [the forced part, from the injected reverse pass]
                    + [alpha * W^T z, with z = (I - alpha M)^-T g]

        where `g` is the reverse pass's own final state -- the sensitivity
        of the functional to the initial state `y_0`, which is itself a
        function of the source through the periodic boundary condition.
        Dropping the second term would leave an answer that looks entirely
        reasonable: MEASURED on an RC ladder the two terms are comparable
        in size (303 against 498 at `l = 0`, 79 against 606 at `l = 1`), so
        neither is a correction to the other.

        Still ONE transposed solve per sideband whatever the number of
        sources, which is the property pnoise needs.  Agreement with the
        `m` forward driven solves: 9.2e-16 / 3.3e-16 / 1.3e-15 at
        `l = 0 / 1 / -2`.

        ⚠ WAS SOLVED-HISTORY ONLY, like the reverse pass; B8 lifted both.
        """
        import scipy.sparse.linalg as spla
        ## ⚠ NO LONGER GEAR-ONLY (B8) -- see the adjoint row.
        fp = pss.factored_period()

        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the sideband row')
        m = pss.cir.n - 1
        n = fp.width
        T = float(fp.T)
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        w0 = 2.0 * np.pi / T
        alpha = np.exp(-2j * np.pi * float(freq) * T)
        ls = np.atleast_1d(np.asarray(sidebands, dtype=int))

        ## ⚠ A HARD BOUND, NOT A HEURISTIC.  Okumura et al. eq. (32): the
        ## maximum frequency the analysis can speak about is the grid's own
        ## `w_max`, so `|l| <= (w_max - w0)/ws`.  You cannot alias down from
        ## above what the grid can represent, and a ratio test on the
        ## accumulated power operates INSIDE this ceiling rather than
        ## instead of it -- an implementation carrying only the ratio test
        ## terminates for the wrong reason.  The grid's ceiling here is its
        ## Nyquist, `N/2` harmonics of the period.
        lmax = N // 2
        bad = ls[np.abs(ls) > lmax]
        if len(bad):
            raise ValueError(
                'PAC: sideband %s is above the grid\'s Nyquist (|l| <= %d '
                'at %d points per period). Nothing can alias down from '
                'above the maximum frequency the grid represents, so this '
                'is not a tolerance to relax -- use a finer period grid.'
                % (bad.tolist(), lmax, N))

        d = np.zeros(m, dtype=complex)
        if np.isscalar(output):
            d[int(output)] = 1.0
        else:
            out = np.asarray(output, dtype=complex).ravel()
            d[:len(out)] = out

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        tol = max(self.KRYLOV_FACTOR * pss.par.reltol, 1e-14)
        phase = np.exp(2j * np.pi * float(freq) * tms[1:N + 1])

        rows = np.zeros((len(ls), m), dtype=complex)
        for li, l in enumerate(ls):
            ## ⚠ THE PHASE OF THE INPUT COMES OUT FIRST, and getting this
            ## wrong is self-consistent rather than loud.  What is
            ## T-PERIODIC is `v(t) = y(t) exp(-j w t)`, not `y` -- so the
            ## sideband set is the DFT of `v`, which is what `solve` takes.
            ## Decomposing `y` instead gives a Dirichlet kernel smeared
            ## across every `l` whenever `f` is not a multiple of `1/T`,
            ## and it AGREES with a forward reference written the same way,
            ## so only a check against a circuit whose answer is known
            ## independently catches it.
            if fp.kind == 'dirk':
                ## the source couples through the lower-triangular stages; the
                ## sequential fold carries it (verified vs forward driven solves
                ## and the bespoke trbdf2 fold) -- see `_sideband_forced_dirk`
                forced, g = pss._sideband_forced_dirk(fp, freq, l, d)
            elif fp.kind == 'full':
                ## the source couples through ALL THREE stages (A (x) B), which
                ## needs the coupled three-vector fold -- see
                ## `_sideband_forced_full`
                forced, g = pss._sideband_forced_full(fp, freq, l, d)
            else:
                inject = ((np.exp(-1j * (float(l) * w0 + 2.0 * np.pi
                                         * float(freq)) * tms[:N]) / N)[:, None]
                          * d[None, :])
                g, ts, _st = fp.matvec_transposed(
                    np.zeros(n, dtype=complex), collect=True, inject=inject)
                forced = -np.tensordot(phase, np.asarray(ts), axes=(0, 0))
            ## ⚠ ON AN OSCILLATOR THIS OPERATOR IS SINGULAR AT EVERY
            ## HARMONIC and near-singular around them, which is exactly
            ## where phase noise is measured.  The deflated route borders
            ## the pole out and carries `1/(1 - alpha)` analytically; on a
            ## driven circuit there is no pole and the plain solve is both
            ## correct and cheaper.
            if getattr(pss, 'autonomous', False):
                z = self._deflated_solve(pss, alpha, g, transposed=True,
                                         tol=tol)
            else:
                z = self._gmres_checked(
                    A, g, tol, 'the adjoint solve at sideband %d' % l)
            rows[li] = forced + alpha * pss._forced_replay_transposed(
                fp, freq, z)
        self.matvecs = count[0]
        return rows

    ## How small a sideband's contribution must be, relative to the running
    ## total, before the accumulation stops.  Okumura et al.: powers are
    ## "accumulated until their contributions become negligible".
    ALIAS_RATIO_TOL = 1e-9

    def mixer_response(self, pss, f_out, output, sidebands=(-1, 0, 1)):
        """Every input band landing on `f_out` — a `SidebandResponse`.

        For each `l`, the transfer from a source at `f_in = f_out - l f0`
        to the output at `f_out`, which is one
        `adjoint_sideband_row` per sideband because each has its own input
        frequency.  Cost is `len(sidebands)` rows.

        ⚠ A NEGATIVE INPUT BAND IS REFUSED RATHER THAN FOLDED.  For
        `f_out < l f0` the input frequency comes out negative; physically
        that band is the conjugate of `|f_in|`, and quietly taking the
        absolute value would return the right magnitude attached to the
        wrong label -- exactly the mislabelling this object exists to
        prevent.  Ask for the sidebands whose input bands exist.
        """
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        ls = [int(l) for l in sidebands]
        ins, rows = [], []
        for l in ls:
            f_in = float(f_out) - l * f0
            if f_in < 0.0:
                raise ValueError(
                    'PAC.mixer_response: sideband %d would be fed from '
                    '%.12g Hz, which is negative. That band is the '
                    'conjugate of %.12g Hz, and returning it under the '
                    'label %d would attach a right magnitude to a wrong '
                    'name. Request sidebands whose input bands exist, or '
                    'move f_out.' % (l, f_in, abs(f_in), l))
            row = self.adjoint_sideband_row(pss, f_in, output, sidebands=l)
            ins.append(f_in)
            rows.append(np.asarray(row).reshape(-1))
        return SidebandResponse(f_out, f0, ls, ins, np.array(rows))

    ## ⚠ THE FOLD BELOW IS FOR DRIVEN CIRCUITS.  It is a frequency-conversion
    ## computation and is complete for one; for an AUTONOMOUS oscillator it is
    ## structurally incomplete -- the near-carrier phase-noise skirt is not a
    ## conversion effect (Rizzoli, Mastri & Masotti, MTT 42-807, 1994).  Free-
    ## running phase noise goes through the Floquet/PPV stack instead; see
    ## `oscillator_spectrum` for why the two cannot be unified and why the wrong
    ## one still returns a plausible number.
    def pnoise(self, pss, freq, output, ratio_tol=None, maxsidebands=None,
               modulated=False, cyclostationary=False):
        """TIME-AVERAGED output noise PSD at `freq`, sidebands folded in.

        ⚠ `cyclostationary=True` IS THE CONSTRUCTION FOR A BIAS-DEPENDENT
        `CY` (2026-09-08, on the corrected Okumura reading).  A source whose
        PSD follows the orbit is white noise `xi` MODULATED by
        `B(t) = sqrt(CY(x(t)))`, a T-periodic matrix with Fourier
        coefficients `B_k` -- read off the PSS samples by one DFT, no
        window count `p` at all (Okumura's windows are a piecewise-constant
        approximation of exactly this, and their boxcar coefficients its
        crude version).  `xi`'s band at `g_p = f - p f0` reaches the output
        at `f` through EVERY modulation harmonic `k` and the sideband row
        `a_{p-k}` (source at `g_p + k f0`, output at `f`) -- the SAME rows
        the stationary fold computes -- COHERENTLY over `k` (one white band,
        one realisation) and incoherently over `p`.  Summing the bands
        turns the square root into the PSD's OWN harmonics `P_j` (the DFT
        of `CY(x(t))`, no matrix square root anywhere):

            S(f) = sum_{l,l'} a_l P_{l'-l} a_{l'}^H,

        which is exact on the grid (⚠ the sqrt-modulation form, tried
        first, left a 2.8e-5 residual tied to the modulation's zero
        crossings; this form agrees with the stationary side to 9e-16).
        Constant `CY` gives `P_0 = CY` and nothing else, and the sum
        collapses to the stationary `sum_l a_l CY a_l^H` -- Okumura's
        `p = 1` case, pinned to machine precision.  The cost is the
        stationary fold's (the rows dominate; the double sum is free).
        A coloured source is folded band by band (each white band the
        modulation reaches carries its own `CY`; see `_cyclostationary_fold`).
        Like the stationary fold this is a LOWER bound at a sideband cap.
        ⚠ Coherence is the whole content: `modulated=True` (the cycle-
        averaged `CY`, Hull & Meyer's stationary equivalent) keeps the
        power and drops the correlation between sidebands, and the two
        differ wherever the modulation has harmonics -- measured on a
        driven multiplier, and the identity against the STATIONARY fold of
        the same physics written as a white source through a periodically
        varying gain is the gate (`test_..._cyclostationary_...`).  ON A
        MOS (2026-09-09): an EKV stage switched by a 1 MHz LO whose channel
        noise passes through a second EKV switched by the same LO reads
        0.376 of the cycle average (thermal, every offset) and 0.32-0.44
        with flicker at ten times thermal; the switch's own channel noise
        is largest when its channel shunts it.  ⚠ Noise that reaches the
        output through a time-INVARIANT transfer (a single stage's drain
        into an RC load) gives cyc = cycle average to four digits, since
        only P_0 survives -- the construction shows only where the
        modulated noise crosses a periodically varying transfer.  Cost:
        white = the cycle average's; coloured (any frequency-dependent
        CY, a negligible flicker coefficient included) ~6x.
        ⚠ FLICKER, AND WHAT OKUMURA'S EQ. 23 MEANS HERE (measured
        2026-09-08): a coloured source is folded band by band, and against
        the stationary fold of the same SEPARABLE physics (a stationary
        flicker source through a periodically varying gain) it is exact --
        1.000000 -- as long as the modulation is SIGN-DEFINITE.  When the
        modulation changes sign the two are DIFFERENT physics (0.56 / 1.33,
        grid-independent to six digits): for white noise `m xi` and `|m| xi`
        are one process, for a coloured one whose correlation spans the
        sign change they are not (`R(t,t') = m(t) m(t') R_c(t-t')` keeps
        the sign product), and a PSD cannot carry the sign -- so this fold,
        like the HDL model feeding it, is the `|m|` one.  That is Okumura's
        "cannot be modeled as a cyclostationary process by using this
        method, because it has very long time constants" in concrete form.
        A flicker source with a bias-dependent coefficient gets the `|m|`
        number, correct when its modulation does not change sign.

        Returns `(S, sidebands_used)`.  `S` is the one-sided
        **time-averaged** PSD at the output, in the same units as
        `analysis_ss.Noise`'s `Svnout`.

        ⚠ "TIME-AVERAGED" IS NOT A HEDGE, IT IS THE SPECIFICATION, and
        saying so is the whole of this paragraph's job.  TWO SEPARATE
        MECHANISMS make an output noise cyclostationary, and only one of
        them is about the sources:

          1. bias-dependent sources modulated by the time-varying operating
             point -- this is what `_cy_reduced` refuses, because the
             stationary sum would be the wrong model;
          2. the PERIODIC SOURCE-TO-OUTPUT TRANSFER FUNCTION -- which
             applies even when every source is stationary.  A circuit whose
             only noise is the thermal noise of constant resistors STILL
             has cyclostationary output noise.

        The sideband sum here handles (2) correctly and returns its TIME
        AVERAGE.  That is the right answer for most uses and it is
        incomplete for two ordinary RF topologies, both named by Kundert
        (*Introduction to RF Simulation*, v2 2003 -- relayed from the docs
        session, cited not verified here): a NONLINEAR SUBSEQUENT STAGE
        ("an oscillator drives a limiter ... the same is true when an
        oscillator drives a mixer"), and CASCADED STAGES OFF A SHARED
        REFERENCE, where "the second mixer is synchronous with, and tracks
        the variations in, the cyclostationary noise of the first."  The
        test is whether anything downstream can track the PSD's variation:
        if it cannot, the phase is unknown to it and the time average is
        sufficient.

        ⚠ AND A SCALAR CANNOT CARRY WHAT IS MISSING.  Cyclostationary noise
        is CORRELATED between frequencies separated by `k f0`, where
        stationary noise has no correlation between different frequencies
        at all.  This returns one number per output frequency, so it does
        not represent that correlation -- deliberately, and stated here
        rather than left for a caller to discover by getting a wrong answer
        in one of the two topologies above.

            S(f) = sum_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

        Noise entering at `f - l f0` leaves at `f` through sideband `l`, and
        white sources in disjoint bands are uncorrelated, so the bands add
        in POWER.  Each `h_l` is one adjoint row -- one transposed solve for
        every source in the circuit, which is the whole reason this is
        affordable.

        ⚠ A PRECONDITION FOR THE FIRST COLOURED SOURCE, recorded here
        because it is unreachable today and will be silent when it is not.
        A 1/f source is singular at DC, and folding puts a copy of that
        singularity at EVERY harmonic.  A commercial RF simulator: "place a cluster of
        frequencies near each harmonic ... but AVOID PUTTING FREQUENCY
        POINTS PRECISELY ON THE HARMONICS ... you run the risk of
        generating absurd noise totals because a very narrow noise peak
        artificially has its apparent width greatly magnified by a large
        frequency, and has its amplitude exaggerated by placing a point
        precisely at the singularity."  Plausible nonsense, no error
        raised.  Every source in the discrete library is white, so `freq`
        landing on `k f0` is harmless now; it stops being harmless the day
        one is not.

        ⚠ AND AN OSCILLATOR IS NOT THIS FUNCTION'S PROBLEM AT ALL.  A
        driven circuit's output noise is cyclostationary; an AUTONOMOUS
        one's is STATIONARY, and structurally so -- "cyclostationarity in
        the oscillator's output would, by definition, imply a time
        reference ... noisy autonomous systems cannot provide a perfect
        time reference" (Demir 2002).  That is the physical counterpart of
        `I - M kron M` being exactly singular for an oscillator
        (`test_no_periodic_covariance_exists_for_an_oscillator`): there is
        no cyclostationary object to compute, not a hard one.  Oscillator
        phase noise is a different output shape entirely -- a closed form
        in a few scalars with no frequency sweep -- and is not built.

        ⚠ WHY CYCLOSTATIONARY IS NOT BUILT -- AND THE REASON RECORDED HERE
        FIRST WAS WRONG.  This said the cross terms need "the `R_{m,n}`
        construction from section III-B", unread, as though the window
        Fourier coefficients were an exotic object.  They are not.  The
        model is `c(t) = sum_m n_m(t) w_m(t)` with `w_m` a T-periodic
        RECTANGULAR window over interval `m`, non-overlapping -- so
        `W_{m,k}` is the Fourier series of a BOXCAR, closed form, a `sinc`
        times a phase.  The `n_m` are taken UNCORRELATED across intervals,
        justified because `H(jw,t)` is time-invariant within each one, so
        the sum is INCOHERENT over `m` and coherent only over `k` within a
        single interval.  Nothing there is missing.
        ⚠ THE ACTUAL BARRIER IS COST, which is a different decision -- and
        the cost as first recorded here was OVER-STATED (verified at the
        source by the docs session, 2026-09-08, Okumura et al. 1993).  The
        source count is p x (noisy devices) where p is the number of
        intervals over which "H(jw,t) is time-invariant within each
        interval" -- set by how fast the transfer varies (their Fig. 2 has
        FIVE windows; a switching circuit moves fast only at transitions),
        NOT by the integration grid: the earlier "500-point grid x 50
        devices = 25 000 sources" tied p to the timestep and was high by
        (timepoints)/p, an order or two.  The reported noise analysis ran at
        14.1x the PSS per frequency point (1086 s vs 77 s), "because all
        aliasing components need to be computed" -- and the NEXT sentence,
        elided before: "it is expected that this problem can be greatly
        alleviated using a vectorization technique, because most of the
        computational power is used to solve linear problems" -- which is
        the batched JAX path this tree already carries.  Whether that is
        affordable is unmeasured here.
        ⚠ AND THE METHOD CANNOT MODEL FLICKER (p. 585, verbatim): "Flicker
        noise generated under a periodic large signal excitation cannot be
        modeled as a cyclostationary process by using this method, because
        it has very long time constants and thus equation (23) does not
        hold" -- eq. 23 being the uncorrelated-across-intervals assumption.
        Their fallback is that flicker "may exist as independent noise
        sources which are practically modeled as stationary random
        processes".  So the construction covers cyclostationary thermal
        and shot noise, and NOT one of the three mechanisms `_cy_reduced`
        names it as the precondition for.
        ⚠ AND ITS AUTHORS LEFT THE PHYSICS OPEN: "it is further necessary
        to discuss the correspondence between the actual physical phenomena
        of noises and this modeling".  The windowed-stationary
        decomposition is a numerical construct, and its fidelity to a real
        device is not settled by its numerical validation.
        (Verified at the source 2026-09-08.  The boxcar is the paper's own
        closed form, R_{m,n} = (h_m/T) Sa(n w_s h_m/2) exp(-j n w_s (tau_{m-1}
        + h_m/2)), p. 585 -- an earlier line here claimed the observation as
        ours.)

        ⚠ STATIONARY SOURCES ONLY, AND IT CHECKS -- mechanism (1) above.
        Okumura's cyclostationary model windows each source to a single
        timestep, and the windows'
        Fourier coefficients then CORRELATE the sidebands -- they stop
        adding in power, and the cross terms need the `R_{m,n}` construction
        from his §III-B.  Every noise source in this element library is
        bias-INdependent (a resistor's `4kT/R` does not read `x` at all), so
        the stationary formula is exact for them; a compact device with a
        bias-dependent `CY` is not covered, and this raises rather than
        returning a number that is quietly the wrong model.

        ⚠ `maxsidebands` IS AN ACCURACY KNOB HERE AND A REPORTING KNOB IN
        `PAC.solve`, WHICH IS THE OPPOSITE OF HOW IT READS.  A commercial RF simulator's
        own documentation states the inversion (relayed, cited not verified
        here): reducing sidebands "affects only the amount of information
        generated, not its quality.  HOWEVER, NOISE SOURCES GENERATE
        SIGNALS AT ALL FREQUENCIES, and therefore with PNoise, reducing the
        number of sidebands acts to REDUCE THE NUMBER OF NOISE
        CONTRIBUTIONS in the output and so REDUCES THE ACCURACY of the
        result."  A driven signal lives at the frequencies it is driven at,
        so dropping sidebands drops answers you did not ask for; noise
        lives at all of them, so dropping sidebands drops power that
        belonged in the total.  Capping it here always makes `S` a LOWER
        bound, never a cheaper estimate of the same number.

        ⚠ AN OBSERVABLE SYMPTOM WORTH KNOWING BEFORE IT IS SEEN.  For an
        oscillator `Phi(T) - I` is singular and its null vector IS THE PPV,
        so a near-carrier noise computation is ill-conditioned by
        construction.  Gourary et al. name what that looks like: "the
        standard time domain noise analysis yields FLAT PSD CURVES OR
        CURVES WITH UNEXPECTED SLOPE NEAR THE OSCILLATION FREQUENCY."  If
        oscillator noise ever comes out flat near the carrier, that is the
        singularity -- not the physics, the noise models or the source
        definitions -- which points at the right layer immediately.  The
        published removal (Gourary et al., eq. 27/28: replace the output
        row of J^T by u^T; verified at the source by the docs session,
        2026-09-08) IS built here as `_deflated_solve`, which borders with
        BOTH null vectors and is the better conditioned of the two; it is
        wired into `adjoint_sideband_row` (so into this method) and, since
        2026-09-08, into `PAC.solve` and `adjoint_transfer_row` as well
        (`PAC.deflated` says which route ran).  The plain solve's relative
        error near a harmonic is `eta / (2 pi df/f0)` with `eta =
        |lambda_1 - 1|` the computed unit multiplier's displacement --
        measured 1.1e-12 (Q = 16) and 1.8e-13 (Q = 100) under radau, so
        under the default integrator the unguarded band sat inside
        `HARMONIC_GUARD` and the wiring is hygiene; under gear at
        df/f0 = 1e-10 the plain solve refuses outright (GMRES residual
        1.7e-6) where the deflated one carries the pole to 1 %.

        ⚠ TWO STOPPING RULES, AND THE BOUND IS NOT THE RATIO TEST.  The
        accumulation stops when a sideband pair adds less than `ratio_tol`
        of the running total -- and it can never pass `|l| <= N/2`, the
        grid's own Nyquist, because nothing aliases down from above the
        maximum frequency the grid represents (eq. 32).  An implementation
        with only the ratio test terminates for the wrong reason and, on a
        coarse grid, after summing harmonics its own grid cannot carry.

        GATED against `analysis_ss.Noise` on a linear circuit, where the
        sidebands vanish and this must reduce to the stationary answer --
        Okumura's own `p = 1` case, "exactly the same as that derived for a
        stationary noise".  Measured ratio 1.000000, with every `l != 0`
        contributing ~1e-32 of the total.

        ⚠ AND THE TIME-AVERAGE CHOICE MATCHES THE REFERENCE IMPLEMENTATION,
        which is worth recording because it was documented above as a
        deliberate scope decision and could have been the wrong one.
        a commercial RF simulator's theory notes on PNoise and QPnoise, both: "THE TIME-AVERAGE of
        the noise at the output of the circuit is computed in the form of a
        spectral density versus frequency."  Same quantity, same
        limitation.  (Relayed from the docs session; cited, not verified
        here.)
        """
        self._check_circuit(pss)
        ## pnoise folds sidebands through the ADJOINT (adjoint_sideband_row ->
        ## _forced_replay_transposed), whose two-stage chained transpose is
        ## not built for TR-BDF2, so it falls back to a Gear-2 twin -- see
        ## `_adjoint_host`.  (covariance/oscillator_covariance use the built
        ## TR-BDF2 Lyapunov injection via `_lyapunov_host`.)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        m = pss.cir.n - 1
        N = len(fp.steps)
        T = float(fp.T)
        f0 = 1.0 / T
        tol = self.ALIAS_RATIO_TOL if ratio_tol is None else float(ratio_tol)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)

        w = 2.0 * np.pi * float(freq)
        ## ⚠ `modulated=True` IS HULL & MEYER'S ROUTE, NOT A TOLERANCE
        ## RELAXATION.  Off, a bias-dependent `CY` raises, because the
        ## stationary sum would be the wrong model.  On, the source is
        ## replaced by ONE stationary source at the CYCLE-AVERAGED bias and
        ## the modulation is carried by `H_l` -- which is the standard
        ## treatment of exactly this case, and the only route to MOS
        ## pnoise, since no physically correct MOS noise model has a
        ## state-independent `CY`.
        colour = None
        if cyclostationary:
            ## the stop rule and the harmonic probes below run on the
            ## cycle-averaged power (the modulation's B_0 B_0^H); the fold
            ## itself is the convolution after the rows are gathered.
            ## The colour model (see `_cy_colour_model`) is fitted ONCE
            ## here and serves both: the 34 orbit sweeps of the stop rule
            ## were 1.1 s of a 4.1 s call after the fold itself was cut.
            colour = self._cy_colour_model(pss, float(freq), f0)
            if colour is None:
                cyfn = self._cy_cycle_averaged
            else:
                fp_ = pss.factored_period()
                hs_ = np.diff(np.asarray(fp_.times, dtype=float))
                def cyfn(pss_, w_, _m=colour, _h=hs_):
                    Cs = _m(w_)
                    ns = min(len(_h), Cs.shape[0])
                    return np.einsum('k,kij->ij', _h[:ns], Cs[:ns]) / float(_h[:ns].sum())
        else:
            cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        cy = cyfn(pss, w)

        ## ⚠⚠ ON A HARMONIC, A SIDEBAND FOLDS THE SOURCES TO DC -- AND
        ## SOME DEVICE MODELS ARE NOT DEFINED THERE.  Sideband `l`
        ## evaluates `CY` at `f - l f0`, so `f = k f0` evaluates it at
        ## ZERO.  A `1/f` term is infinite there; and MEASURED on
        ## `PspMosLongChannel`, a flicker term with its coefficient set to
        ## ZERO is `0/0` and returns `nan`:
        ##
        ##     fnt=1, nfa=0        CY(f=0) = nan   <- DISABLED flicker
        ##     fnt=1, nfa=8e22     CY(f=0) = inf   <- the real singularity
        ##
        ## ⚠ THE FIRST IS THE NASTIER ONE: a caller who sets `nfa = 0`
        ## believing flicker is off still gets `nan` out of `pnoise`, with
        ## no exception anywhere.
        ##
        ## ⚠ AND IT IS NOT "HARMONICS ARE BAD".  A driven divider with
        ## white sources returns 1.490351e-17 at exactly `f0`, and at
        ## `2 f0` and `3 f0` -- the fold to DC is harmless when the
        ## sources are defined there.  So this checks the SOURCES at the
        ## frequency that will actually be used, rather than refusing a
        ## harmonic on principle.
        f0_ = 1.0 / float(pss.period)
        lscan = max(1, int(maxsidebands or 8))
        offs = np.abs(float(freq) - np.arange(-lscan, lscan + 1) * f0_)
        near = float(np.min(offs))
        if near <= self.HARMONIC_GUARD * f0_:
            probe = cyfn(pss, 2.0 * np.pi * near)
            if not np.all(np.isfinite(np.asarray(probe))):
                raise ValueError(
                    'PAC.pnoise: %.12g Hz sits on a harmonic of %.12g Hz, '
                    'so a sideband folds the noise sources to DC -- and at '
                    'DC this circuit\'s CY is not finite. A 1/f term is '
                    'infinite there; a flicker term whose COEFFICIENT IS '
                    'ZERO is 0/0 and gives nan, so disabling flicker does '
                    'not avoid this. Offset from the harmonic: a commercial RF simulator\'s '
                    'own advice is to cluster frequencies NEAR each '
                    'harmonic and never place one ON it.'
                    % (float(freq), f0_))

        ## ⚠ AND THE STEEP REGION BESIDE IT IS A SWEEP HAZARD RATHER THAN
        ## A WRONG NUMBER, so it warns instead of raising.  MEASURED with
        ## a real flicker source: the plateau is 1.321483e-16, `f0 + 1` Hz
        ## gives 1.321766e-16 and `f0 + 0.01` Hz gives 1.350047e-16 -- 2%
        ## high, finite, entirely plausible.  The VALUE is right; a grid
        ## that lands there by accident integrates a spike it never
        ## resolved.  A commercial RF simulator: "you run the risk of generating absurd
        ## noise totals because a very narrow noise peak artificially has
        ## its apparent width greatly magnified".
        elif near < 1e-6 * f0_ and float(freq) > 0.0:
            cy_hi = cyfn(pss, 2.0 * np.pi * max(float(freq) * 2.0, f0_))
            if not np.allclose(cy, cy_hi, rtol=1e-9, atol=0.0):
                warnings.warn(
                    'PAC.pnoise: %.12g Hz is %.3g Hz from a harmonic of '
                    '%.12g Hz and a source has a frequency-dependent CY, '
                    'so the folded density varies steeply here. The VALUE '
                    'is correct; a swept grid landing this close will '
                    'misrepresent the integrated total. Cluster near each '
                    'harmonic deliberately rather than by accident.'
                    % (float(freq), near, f0_), RuntimeWarning, stacklevel=2)

        total = 0.0
        used = []
        quiet = 0
        self.alias_stop = 'bound'
        rows = {}
        for l in range(0, lmax + 1):
            step = 0.0
            for sl in ((0,) if l == 0 else (l, -l)):
                fin = float(freq) - sl * f0
                h = self.adjoint_sideband_row(pss, fin, output, sl)[0]
                rows[sl] = np.asarray(h, dtype=complex)
                step += float(np.real(h @ cyfn(
                    pss, 2.0 * np.pi * fin) @ np.conj(h)))
                used.append(sl)
            total += step
            if total > 0 and abs(step) < tol * abs(total):
                ## ⚠ TWO QUIET PAIRS, NOT ONE.  A single sideband can come
                ## back near zero by symmetry while its neighbours do not,
                ## and stopping there would truncate a series that had not
                ## converged.
                quiet += 1
                if quiet >= 2:
                    self.alias_stop = 'ratio'
                    break
            else:
                quiet = 0
        self.sidebands_used = used
        if cyclostationary:
            total = self._cyclostationary_fold(pss, float(freq), rows, model=colour)
        ## ⚠ WHICH RULE STOPPED IT IS PART OF THE ANSWER.  Ending on the
        ## ratio test means the series converged; ending on the Nyquist
        ## bound means the grid ran out before the series did, and the
        ## number is a LOWER bound on the folded noise -- every sideband
        ## above the grid's own maximum frequency is missing, not small.
        ## A strongly switching circuit does this readily: measured on a
        ## driven diode at 80 points per period, the accumulation reached
        ## l = +-39 without the ratio test ever firing, while folding was
        ## already contributing 62% of the total.
        if self.alias_stop == 'bound' and lmax > 0:
            warnings.warn(
                'PAC.pnoise: the sideband accumulation stopped at the '
                "grid's Nyquist (|l| = %d at %d points per period), not "
                'because the contributions became negligible. Sidebands '
                'above the grid\'s maximum frequency are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Re-solve the PSS on a finer period grid and compare.'
                % (lmax, N),
                RuntimeWarning, stacklevel=2)
        return total, used

    def _cy_harmonics(self, pss, w):
        """`P_j`: the Fourier coefficient matrices of `CY(x(t), w)` over the
        orbit, `(N, n, n)` indexed like `numpy.fft.fftfreq`.  `P_0` is the
        cycle average; `P_j` with `j != 0` carry the modulation and vanish
        for a bias-independent source.  No square root: the fold uses the
        PSD's own harmonics (`a P a^H`), which is exact on the grid, where
        a sqrt-modulation route (tried first) left a 2.8e-5 residual: the
        square root of a PSD that crosses zero has a kink, its harmonic
        tail decays slowly, and the convolution's window -- the sidebands
        the ratio stop kept, 7 here -- truncated it (measured -1.6e-4 /
        -2.8e-5 / -3e-7 at 5 / 7 / 17 sidebands).  The PSD's harmonics
        decay fast, so this form is exact at any window."""
        fp = pss.factored_period()
        irn = pss.irefnode
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = len(fp.steps)
        Cs = []
        for k in range(nsamp):
            xr = np.asarray(xs[:, k], dtype=float).ravel()
            xf = xr if xr.shape[0] == pss.cir.n else np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        Cs = np.asarray(Cs, dtype=complex)
        return np.fft.fft(Cs, axis=0) / Cs.shape[0]

    def _cy_samples(self, pss, w):
        """`CY(x(t_k), w)` over the orbit, reduced, `(N, n, n)` complex."""
        fp = pss.factored_period()
        irn = pss.irefnode
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = len(fp.steps)
        Cs = []
        for k in range(nsamp):
            xr = np.asarray(xs[:, k], dtype=float).ravel()
            xf = xr if xr.shape[0] == pss.cir.n else np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    @staticmethod
    def _sqrt_harmonics_of(Cs):
        """The DFT of the symmetric square root of the sampled `CY` (see
        `_cy_sqrt_harmonics`), for a `(N, n, n)` array already in hand."""
        Bs = []
        for cyk in Cs:
            cyk = 0.5 * (cyk + cyk.conj().T)
            lam, U = np.linalg.eigh(cyk)
            lam = np.clip(np.real(lam), 0.0, None)
            Bs.append((U * np.sqrt(lam)[None, :]) @ U.conj().T)
        Bs = np.asarray(Bs, dtype=complex)
        return np.fft.fft(Bs, axis=0) / Bs.shape[0]

    def _cy_colour_model(self, pss, f, f0):
        """Fit `CY(x(t), w) = A(t) + B(t) (w1/w)^ef` entry by entry from three
        frequencies and verify at a fourth; return a callable `w -> (N, n, n)`
        or None when the fit fails anywhere (the caller then evaluates the
        circuit per band, as before).  The exponent is per entry, found by
        a bracketed root find on the ratio of differences, so a mix of
        flicker exponents across sources is fine; a white entry (B = 0)
        needs no exponent."""
        from scipy.optimize import brentq
        f = abs(float(f))
        ## three to fit, two to verify: one BETWEEN the fit points and one
        ## at the FAR end of the band range the fold reaches (up to
        ## ~(N/2 + lmax) f0), so a shape that is not thermal-plus-flicker
        ## is caught where the model would have been extrapolating
        ws = [2.0 * np.pi * x for x in (max(f, 1e-3 * f0), 3.0 * f0 + f, 10.0 * f0 + f,
                                        2.0 * f0 + f, 150.0 * f0 + f)]
        C1, C2, C3, C4, C5 = (self._cy_samples(pss, w) for w in ws)
        N, n, _ = C1.shape
        A = np.zeros_like(C1); B = np.zeros_like(C1); EF = np.zeros((N, n, n))
        scale = max(float(np.max(np.abs(C1))), 1e-300)
        w1, w2, w3, w4, w5 = ws
        for k in range(N):
            for i in range(n):
                for j in range(n):
                    c1, c2, c3 = C1[k, i, j], C2[k, i, j], C3[k, i, j]
                    if abs(c1 - c2) <= 1e-12 * scale and abs(c2 - c3) <= 1e-12 * scale:
                        A[k, i, j] = c1
                        continue
                    ratio = (c1 - c2) / (c2 - c3)
                    def g(ef, ratio=ratio):
                        g1, g2, g3 = 1.0, (w1 / w2) ** ef, (w1 / w3) ** ef
                        return float(np.real((g1 - g2) / (g2 - g3) - ratio))
                    try:
                        ef = brentq(g, 0.05, 4.0, xtol=1e-12)
                    except ValueError:
                        return None
                    g2 = (w1 / w2) ** ef
                    Bv = (c1 - c2) / (1.0 - g2)
                    A[k, i, j] = c1 - Bv; B[k, i, j] = Bv; EF[k, i, j] = ef
        for wv, Cv in ((w4, C4), (w5, C5)):
            if float(np.max(np.abs(A + B * (w1 / wv) ** EF - Cv))) > 1e-8 * scale:
                return None
        def model(w):
            return A + B * (w1 / float(w)) ** EF
        return model

    def _cy_sqrt_harmonics(self, pss, w):
        """`B_k`: the DFT of the symmetric square root of `CY(x(t), w)` over
        the orbit, `(N, n, n)`, for the band-resolved (coloured) fold."""
        fp = pss.factored_period()
        irn = pss.irefnode
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = len(fp.steps)
        Bs = []
        for k in range(nsamp):
            xr = np.asarray(xs[:, k], dtype=float).ravel()
            xf = xr if xr.shape[0] == pss.cir.n else np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            cyk = np.asarray(cyk, dtype=complex)
            cyk = 0.5 * (cyk + cyk.conj().T)
            lam, U = np.linalg.eigh(cyk)
            lam = np.clip(np.real(lam), 0.0, None)
            Bs.append((U * np.sqrt(lam)[None, :]) @ U.conj().T)
        Bs = np.asarray(Bs, dtype=complex)
        return np.fft.fft(Bs, axis=0) / Bs.shape[0]

    def _cyclostationary_fold(self, pss, freq, rows, model=None):
        """`S(f) = sum_{l,l'} a_l Q_{l,l'} a_{l'}^H` over the gathered
        sideband rows (`rows[l]` = the row for a source at `f - l f0`,
        output at `f`).  WHITE source: `Q_{l,l'} = P_{l'-l}`, the DFT of
        `CY(x(t))` itself -- no square root, exact on the grid (9e-16
        against the stationary fold of the same physics).  COLOURED source
        (`CY` depends on `w`; detected by comparing two bands): the white
        band `p = l + k` shared by rows `l` and `l'` carries its OWN `CY`,
        so `Q_{l,l'} = sum_k B_k^{(l+k)} B_{k+l-l'}^{(l+k) H}` with
        `B^{(p)}` the sqrt-DFT at the band's frequency `|f - p f0|`, summed
        over ALL `N` modulation harmonics `k` (which is what makes the
        square root exact here: the 2.8e-5 of the first version came from
        a window on `k`, not from the root).  ⚠ Measured by the docs
        session on a flicker source: `||P_0||` differs 24x across the bands
        the fold sums, so "the band of l" (the first version's shortcut)
        was a 24x approximation on the case the feature exists for; the
        band-resolved form is pinned against the stationary fold of a
        stationary FLICKER source through the same multiplier.

        COST (2026-09-09): the coloured branch was 6x the white one
        because of the circuit's `CY` (231 bands x 230 samples), not the
        algebra.  With the colour model (`_cy_colour_model`, fitted once
        in `pnoise` and shared with the stop rule) and the pair sum
        vectorised it is 2.2x the white call and below the cycle average
        (1.7 s / 0.8 s / 2.0 s on the switched EKV fixture), exact to
        1e-11 against the per-band evaluation, which remains the fallback
        for a colour the model does not fit."""
        f0 = 1.0 / float(pss.period)
        ls = sorted(rows)
        f = float(freq)
        ## coloured or white?  two bands, same test the stationary path uses
        w_a = 2.0 * np.pi * abs(f - ls[0] * f0)
        w_b = 2.0 * np.pi * max(abs(f) * 2.0, f0)
        Pa = self._cy_harmonics(pss, w_a)
        Pb = self._cy_harmonics(pss, w_b)
        coloured = not np.allclose(Pa, Pb, rtol=1e-9, atol=0.0)
        total = 0.0
        if not coloured:
            P = Pa
            N = P.shape[0]
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ P[(lp - l) % N] @ np.conj(rows[lp]))
            return float(np.real(total))
        ## band-resolved: every white band the modulation harmonics reach.
        ## ⚠ THE COST WAS THE CIRCUIT'S CY, NOT THE ALGEBRA (profiled
        ## 2026-09-09 on a switched EKV stage: 231 bands x 230 samples =
        ## 53 000 CY evaluations, 7.8 s of a 10.4 s fold; the eigen-
        ## decompositions 0.85 s).  Every colour in the library is thermal
        ## plus flicker in 1/f^ef, so THREE evaluations per sample fix each
        ## entry's shape (A + B (w1/w)^ef, ef by a root find on the ratio of
        ## differences), a FOURTH frequency verifies the fit to 1e-8, and
        ## all the bands come from the model with no further circuit
        ## calls; a source whose colour is not of that shape fails the
        ## check and gets the full evaluation as before.
        cache = {}
        Nn = Pa.shape[0]
        if model is None:
            model = self._cy_colour_model(pss, f, f0)
        ## ⚠ A SPECIFICATION LIMIT, NOT AN IMPLEMENTATION ONE (docs session,
        ## 2026-09-08): a coloured source under a modulation that CHANGES
        ## SIGN is not representable by any fold built from a PSD -- the
        ## correlation R(t,t') = m(t) m(t') R_c(t-t') keeps the sign product
        ## and CY cannot carry it -- so this fold, like the HDL model
        ## feeding it, computes the |m| process (measured 0.56 / 1.33 of
        ## the signed one on a flicker source through a zero-crossing
        ## gain, grid-independent; 1.000000000 for a sign-definite gain).
        ## The sign is invisible here; its NECESSARY condition is a PSD
        ## that touches zero along the orbit with a KINK in its square
        ## root, so that is warned on.  ⚠ SCOPE (measured on an EKV stage,
        ## 2026-09-09): a DEVICE's own noise has no sign to lose -- its
        ## modulation is a physical intensity, sqrt(PSD(x(t))) >= 0 IS the
        ## process -- so for intrinsic MOS thermal, shot or flicker noise the
        ## PSD-specified model is the physics and this warning does not
        ## apply; the ambiguity belongs to noise passing through a SIGNED
        ## external gain (the multiplier fixtures).  Okumura's eq. 23
        ## objection to flicker is then the separate, physical question of
        ## whether a trap process is "modulated coloured noise" at all.
        ## The proxy's threshold: a zero crossing SAMPLED on an N-point grid
        ## bottoms out near (pi/N)^2 of the maximum (6e-4 at 200 points on
        ## the gate fixture), while a sign-definite PSD with a ten-fold
        ## swing sits at 1e-2 -- so 1e-2 separates them here; a heuristic,
        ## and it is a warning for that reason.
        Cs0 = np.asarray([np.abs(np.diag(np.fft.ifft(Pa, axis=0)[k])) for k in range(Nn)])
        dmax = Cs0.max(axis=0)
        touches = (dmax > 0) & (Cs0.min(axis=0) <= 1e-2 * dmax)
        ## ⚠ THE ORDER OF THE ZERO (peer): a LINEAR sign crossing m ~ a t
        ## gives sqrt(PSD) ~ |a t|, a first-derivative KINK; a sign-definite
        ## quadratic touch m ~ b t^2 gives sqrt(PSD) ~ b t^2, SMOOTH.  The
        ## circular second difference of sqrt(PSD) is 2|a|h at a kink and
        ## 2b h^2 where smooth -- both shrink under refinement, the smooth
        ## one faster -- so a RAW threshold encodes the grid (5e-3 was safe
        ## at 240 points and a false positive below ~100; peer).  Divided
        ## by h/T and by the maximum it is a DERIVATIVE JUMP, grid-
        ## independent at a kink (2|a|T/s_max ~ 4 pi for a sinusoidal
        ## slope, 12.6 here) and falling as h/T where smooth (~2 (2 pi)^2
        ## h/T: 0.33 at 240 points, 1.0 at 80, 2.0 at 40), so 3 separates
        ## them down to ~50 points per period and the separation grows
        ## with refinement.  Clears the squared-gain case (k V_lo^2, exact
        ## to nine digits) that the touch test alone flagged.  ⚠ STILL
        ## NECESSARY, NOT SUFFICIENT, AND THE DETECTOR'S SENSITIVITY RUNS
        ## INVERSE TO THE EFFECT (peer): the indicator is 12.57 for a
        ## sinusoidal crossing, 0.24 for sign|sin|^1.5 and 0 for sign|sin|^2
        ## -- all sign-changing -- while the discrepancy stays O(1):
        ## MEASURED on the flicker identity with the LO shaped to
        ## v |v|^(p-1), B/A = 0.187 / 1.895 (p = 1, warned), 0.204 / 1.779
        ## (p = 1.5, silent), 0.217 / 1.699 (p = 2, silent) at 0.13 / 1.37
        ## f0; steeper crossings (p = 0.5: 0.161 / 2.046, p = 0.8: 0.175 /
        ## 1.970) err MORE and are caught.  So detector and effect are
        ## aligned for p <= 1 and the silent region is exactly p > 1 (the
        ## crossing flatter than linear): the deviation stays O(1) there
        ## while the indicator falls by orders.  A quiet warning is
        ## therefore not evidence of a small discrepancy.
        kinked = np.zeros_like(touches)
        hT = 1.0 / float(Nn)
        for jj in np.where(touches)[0]:
            sq = np.sqrt(Cs0[:, jj])
            d2 = np.abs(sq - 0.5 * (np.roll(sq, 1) + np.roll(sq, -1)))
            kinked[jj] = bool(d2.max() / (sq.max() * hT) > 3.0)
        if bool(np.any(kinked)):
            warnings.warn(
                'PAC.pnoise(cyclostationary=True): a COLOURED source whose PSD '
                'touches zero along the orbit -- if its modulation changes sign '
                '(a switching gain), no PSD-specified model can represent the '
                'coloured process (Okumura eq. 23 in concrete form), and this '
                'fold computes the |m| one (its square root has a first-derivative '
                'kink at the zero, the signature of a LINEAR sign crossing; a '
                'necessary condition -- a shallow crossing shows no kink and '
                'errs MORE): measured 0.56x and 1.33x of the signed '
                'physics at two offsets on a flicker source through a '
                'zero-crossing gain -- EITHER direction, the sign of the '
                'discrepancy is set by the offset, not the mechanism -- and '
                'exact for a sign-definite one. Only the element knows the sign.',
                RuntimeWarning, stacklevel=3)
        def _B(p):
            key = round(abs(f - p * f0) / f0, 12)
            if key not in cache:
                wp = 2.0 * np.pi * abs(f - p * f0)
                if model is not None:
                    cache[key] = self._sqrt_harmonics_of(model(wp))
                else:
                    cache[key] = self._cy_sqrt_harmonics(pss, wp)
            return cache[key]
        ks = np.fft.fftfreq(Nn, d=1.0 / Nn).astype(int)
        ## every band the sum reaches, stacked once: BB[pi, k] = B_k^{(p)}
        ## with pi = p - pmin.  The (l, l') pair sum is then two fancy
        ## indexings and one einsum instead of N small products in Python
        ## (204 000 `_B` calls, 1.9 s of a 2.8 s fold, before).
        pmin = min(ls) + int(ks.min()); pmax = max(ls) + int(ks.max())
        BB = np.asarray([_B(p) for p in range(pmin, pmax + 1)], dtype=complex)
        for l in ls:
            for lp in ls:
                ## (B B^H)_j = sum_k B_k B_{k-j}^H: the partner index is
                ## k + l - l', NOT k + l' - l -- the mirror was invisible to
                ## the constant-modulation reduction (only k = 0 there) and
                ## read 0.49 / 0.17 on the smooth-modulation flicker identity.
                ## ⚠ NO CIRCULAR WRAP HERE: a partner beyond N/2 would be
                ## paired with the wrong BAND (each band carries its own
                ## weight), harmless in the white P-form and wrong here --
                ## it read 0.56 / 1.33 on the kinked (zero-crossing)
                ## modulation whose coefficients reach N/2.
                kp = ks + l - lp
                ok = np.abs(kp) <= Nn // 2
                kk, kk2 = ks[ok], kp[ok]
                pi = (l + kk) - pmin
                X = BB[pi, kk % Nn]
                Y = BB[pi, kk2 % Nn]
                Q = np.einsum('kij,klj->il', X, Y.conj())
                total += complex(rows[l] @ Q @ np.conj(rows[lp]))
        return float(np.real(total))

    def _cy_cycle_averaged(self, pss, w):
        """`CY` time-averaged over the orbit — Hull & Meyer's construction.

        ⚠⚠ VALID FOR GENTLE MODULATION ONLY, AND IT FAILS AS A FACTOR, NOT A
        PERCENTAGE.  Measured by an external reference-simulator cross-check (2026-09-05)
        on a series switch + shunt capacitor, `pnoise` at 10 kHz against
        a reference simulator, swept over the modulation depth `goff/gon`:

            goff/gon   1        1e-1     1e-2     1e-3     1e-6
            ratio      1.000    4.33     13.2     15.7     16.0

        The 1.000 at the top is what makes the 16 readable: with no
        modulation the cycle average IS the value.  The mechanism is not
        subtle -- the averaged source injects `4kT <g>` (about half the
        on-state current noise) for the WHOLE period, including the hold
        phase, where the node it injects into is 1 Gohm in parallel with
        100 pF; Hull & Meyer's own condition ("none of the large-signal
        state variables may change significantly over the decay time of
        the impulse response") fails there by six orders (100 ns closed,
        0.1 s open).  So this route is for a mixer's `gm`, a bias-
        dependent shot noise -- not for a switch.  ⚠ A switch's noise IS
        reachable exactly: `covariance` and `oscillator_covariance`
        evaluate `CY` at every step and need no averaging (the switched
        capacitor's held variance reads `kT/C` to 1e-4 there).

        ⚠ THIS IS WHAT `_cy_reduced` REFUSES, DONE INSTEAD OF REFUSED, and
        the literature's answer rather than ours.  Hull & Meyer (1993):
        *"cyclostationary noise sources, such as shot noise, may be modeled
        as MODULATED STATIONARY NOISE SOURCES.  The impulse response that
        is calculated INCLUDES THE EFFECT OF THIS MODULATION.  In the case
        of shot noise, the hypothetical stationary noise source has
        spectral density `S_i = 2q Ibar_c`"* with `Ibar_c` the
        cycle-averaged current.

        So the modulation is carried by the RESPONSE, which `pnoise`
        already computes as `H_l`, rather than by the SOURCES.  Okumura's
        route puts one independent stationary source per timestep interval
        per device -- `p` per device, ~25,000 sources on a real circuit.
        This is ONE per device.  Same physics, `p` times cheaper.

        ⚠ AND ITS CONDITION IS CHECKABLE RATHER THAN A BLANKET REFUSAL:
        *"valid when the impulse response duration is much less than the
        time it takes for the mixer circuit to significantly change its
        state ... NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY CHANGE
        SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE RESPONSE."*

        ⚠⚠ WHICH IS THE OPPOSITE OF HIGH-Q, AND THEY SAY SO: *"high-Q
        filters should be avoided, since they cause the impulse response to
        ring, and thus require a very large value of M."*  So this
        construction degrades exactly where `lambda_2 -> 1` -- the same
        boundary as everything else in this class, arriving from a fourth
        direction.  That makes the two constructions COMPLEMENTARY rather
        than competing: Hull & Meyer for fast-settling circuits, Okumura's
        expensive one for the high-Q case that needs it.  `info` reports
        `|lambda_2|` so the caller can see which regime they are in.

        ⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT.  The average
        that matters is over the LARGE-SIGNAL waveform, so `CY` is
        evaluated at every stored state and averaged with the step weights
        -- the same quadrature `diffusion_constant` uses, so the two remain
        comparable.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        tms = np.asarray(fp.times, dtype=float)
        hs = np.diff(tms)
        T = float(fp.T)
        xs = np.asarray(pss.waveform[1], dtype=float)
        m = pss.cir.n - 1
        nsamp = min(len(hs), xs.shape[1])
        acc = None
        for k in range(nsamp):
            xr = xs[:m, k]
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            cyk = np.asarray(cyk, dtype=complex) * hs[k]
            acc = cyk if acc is None else acc + cyk
        if acc is None:
            raise NotImplementedError(
                'PAC: the orbit has no stored samples to average CY over.')
        return acc / float(hs[:nsamp].sum())

    def _cy_reduced(self, pss, w):
        """`CY` with the reference node removed, refusing a moving one.

        ⚠ THE CHECK IS THE POINT.  A bias-dependent `CY` makes the sources
        CYCLOSTATIONARY, and then sidebands stop adding in power -- the
        stationary sum this class computes would be the wrong model, not
        merely an inaccurate one, and nothing downstream would say so.
        Sampled at three states on the converged orbit rather than argued
        from the element types, because a compact model's `CY` reads `x`
        and the discrete library's does not.

        ⚠⚠ AND ITS SCOPE IS WIDER THAN "AN UNUSUAL CASE" -- IT IS A
        BLANKET REFUSAL OF MOS pnoise.  ⚠ AN EARLIER VERSION OF THIS NOTE
        SAID IT WAS UNREACHABLE BECAUSE `PspMosLongChannel.CY` IS
        IDENTICALLY ZERO.  That was measured on a DEFAULT-CONSTRUCTED
        element: the model has channel thermal and flicker noise, and
        `fnt = 0` by default because "an element built without a card is
        noiseless".  With `fnt = 1` it is nonzero, white, and
        bias-dependent -- so the refusal is REACHABLE from a real device
        today, and `modulated=True` is the route past it.  There is no physically correct MOS noise model whose
        `CY` is state-independent: thermal channel noise is
        `4kT gamma g_d0` with `g_d0` bias-dependent, flicker goes as
        `I_D^AF`, gate shot noise as `2qI_G`, and Mahmutoglu & Demir
        (2015) are explicit that trap rates "depend on the voltages across
        the MOSFET which can considerably vary with time during
        large-signal operation".  So the answer to "will a real device
        pass this check" is already determined, and it is no.

        ⚠ THE ROUTE OUT IS THE CYCLOSTATIONARY CONSTRUCTION, NOT A
        DIFFERENT DEVICE MODEL, and that reorders the roadmap: the
        cyclostationary path is not an enhancement for MOS pnoise, it is
        the PRECONDITION -- for the thermal and shot mechanisms.  ⚠ NOT
        for flicker: Okumura's own construction excludes it (p. 585,
        "cannot be modeled as a cyclostationary process by using this
        method, because it has very long time constants"; verified at the
        source 2026-09-08), by the same long-time-constant physics as the
        trap-rate caveat beside it, and falls back to a stationary flicker
        source.  Hull & Meyer (1993) make it affordable -- one
        stationary source per device at the cycle-averaged current, with
        the modulation carried by the impulse response `H_l` that A1
        already computes -- and their worked example IS shot noise
        modulated by the collector current, i.e. exactly this case.  Their
        condition is checkable rather than a blanket refusal, and it fails
        in the familiar direction: a ringing impulse response breaks it,
        so it degrades as `lambda_2 -> 1`.

        ⚠ SECOND-ORDER CONSEQUENCE, WORTH KNOWING BEFORE THE MODEL LANDS.
        This same check is what keeps the Ito/Stratonovich choice out of
        reach (`CY = GG^T`, so a state-dependent `CY` is a state-dependent
        `G`).  Relaxing it for MOS makes the two interpretations diverge,
        and Demir's escape -- "the noise signals are small compared with
        the deterministic signals" -- may NOT carry for trap noise: a trap
        occupancy is a two-state Markov chain rather than a small
        perturbation of a large signal, and the same paper says the state
        dependence "in fact makes the equation nonlinear".  The tell would
        be a discrepancy in a MEAN but not in a variance.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        ## ⚠ THREE STATES ON THE ORBIT -- the first used to be the ZERO
        ## VECTOR, which is on the orbit only by accident, and a linear
        ## time-invariant RC held by a DC clock was refused as
        ## cyclostationary because a switch model read `goff` at v(ck) = 0
        ## (found by an external reference-simulator cross-check, 2026-09-05).  The third
        ## probe is now the stored state half a period in.
        _W = np.asarray(pss.waveform[1], dtype=float)
        _mid = np.delete(_W[:, _W.shape[1] // 2], irn, axis=0)
        states = [np.asarray(_mid, dtype=float).ravel()[:pss.cir.n - 1],
                  np.asarray(fp.x_last, dtype=float).ravel(),
                  np.asarray(fp.x_prev, dtype=float).ravel()[:pss.cir.n - 1]]
        mats = []
        for xr in states:
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cy,) = remove_row_col((cy,), irn, pss.toolkit)
            mats.append(np.asarray(cy, dtype=complex))
        scale = max(float(np.max(np.abs(mats[0]))), 1e-300)
        for other in mats[1:]:
            drift = float(np.max(np.abs(other - mats[0]))) / scale
            if drift > 1e-9:
                raise NotImplementedError(
                    'PAC: this circuit has a BIAS-DEPENDENT CY (varies by '
                    '%.3g over the orbit), so its noise sources are '
                    'cyclostationary. The sidebands are then correlated '
                    'through the window Fourier coefficients and no longer '
                    'add in power, so the stationary sum here would be the '
                    'wrong model rather than an imprecise one. '
                    '⚠ THIS IS THE NORMAL CASE FOR ANY COMPACT MOS MODEL, '
                    'not an exotic one: thermal channel noise is '
                    '4kT.gamma.gd0 with gd0 bias-dependent, flicker goes '
                    'as I_D^AF, gate shot noise as 2qI_G, and trap capture '
                    'and emission rates read the terminal voltages. So '
                    'this is in effect a refusal of MOS pnoise, and the '
                    'route out is the CYCLOSTATIONARY construction rather '
                    'than a different device model -- BUILT: pass '
                    'cyclostationary=True (the PSD\'s own harmonics, exact; '
                    'modulated=True is the cycle average, which drops the '
                    'sideband correlation and read 0.53 of the truth on a '
                    'driven multiplier). Hull & Meyer (1993) '
                    'make it affordable -- ONE stationary source per '
                    'device at the CYCLE-AVERAGED current, with the '
                    'modulation carried by the impulse response, valid '
                    'while no large-signal state variable changes much '
                    'over the impulse response decay time. ⚠ THAT ROUTE '
                    'IS BUILT: pass modulated=True to use it. It is a '
                    'MODEL CHOICE with the validity condition above, not '
                    'a tolerance relaxation, which is why it is opt-in '
                    'and why this refusal is the default.'
                    % drift)
        return mats[0]

    ## How close to a harmonic of `f0` counts as "on" it, as a fraction of
    ## `f0`.
    ##
    ## ⚠ TIGHTENED BY A7.  This used to be the conditioning floor being
    ## accepted, because `sigma_min` of the plain operator falls LINEARLY
    ## with the distance and everything nearer was unusable.  The deflated
    ## solve removes that: its conditioning is FLAT (measured 2.04e-01 from
    ## 0.3 down to 1e-9 of `f0`), so the only remaining reason to refuse is
    ## the physical one -- at an EXACT harmonic `1/(1 - alpha)` is a
    ## division by zero and the response is genuinely unbounded.  So the
    ## guard now excludes only what has no finite answer, not what was
    ## merely hard to compute.
    def _cy_at(self, pss, w, xr):
        """`CY` at ONE reduced state `xr`, reference row/column removed, with
        NO cyclostationarity check -- for the routes that evaluate the
        source at every step and therefore model a modulated source
        exactly (`covariance`, `oscillator_covariance`).  The stationary
        sum in `pnoise` cannot, which is why `_cy_reduced` refuses there.
        """
        irn = pss.irefnode
        xr = np.asarray(xr, dtype=float).ravel()[:pss.cir.n - 1]
        xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
        cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
        (cy,) = remove_row_col((cy,), irn, pss.toolkit)
        return np.asarray(cy, dtype=complex)

    HARMONIC_GUARD = 1e-12

    def _check_circuit(self, pss):
        """The operating point must belong to THIS circuit, not a similar one.

        ⚠ A DRIVEN OSCILLATOR MAKES THIS A CORRECTNESS TRAP RATHER THAN A
        TYPO GUARD.  The natural way to model one is to solve the PSS of
        the bare oscillator and then treat the injection as a perturbation
        -- and it is wrong, because the injection DEVICE is present even
        when its SIGNAL is zero.  Buonomo & Lo Schiavo: "in absence of the
        injection signal, the injection circuit affects the basic LC
        oscillator by CHANGING THE NONLINEARITY OF THE FEEDBACK LOOP ...
        [it] can affect the start-up condition of the basic differential LC
        oscillator OR ITS OSCILLATION AMPLITUDE, or both."

        So the free-running orbit of the circuit-with-the-device is not the
        orbit of the circuit-without-it, and every Floquet quantity built
        on the wrong one inherits the error -- monodromy, PPV, phase noise.
        The analysis would converge and report a plausible number.

        The reference-node check below catches a mismatched `refnode`; it
        cannot catch this, because two circuits differing by one device
        have the same reference node and often the same node count.
        """
        if self.cir is not pss.cir:
            raise ValueError(
                'PAC: this analysis was built on a different circuit object '
                'than the PSS it was handed. If that is deliberate -- e.g. '
                'solving the PSS of a bare oscillator and perturbing a '
                'version with an injection device added -- it is a '
                'CORRECTNESS error, not a bookkeeping one: the injection '
                'device changes the free-running orbit even with its signal '
                'at zero, so the base solution is the wrong one to '
                'linearise about. Solve the PSS on the SAME circuit.')

    def _check_harmonic(self, pss, freq, what):
        """Refuse an autonomous small-signal solve sitting on a harmonic.

        ⚠ `I - exp(-j w T) M` IS SINGULAR AT EVERY HARMONIC OF `f0`, NOT
        JUST AT DC, AND ONLY FOR AN OSCILLATOR.  At `w = k w0` the factor
        `exp(-j w T)` is 1 and the operator is `I - M`, which an autonomous
        circuit's unit multiplier makes singular.  MEASURED on van der Pol,
        `sigma_min(I - exp(-jwT) M)`:

            offset/f0    0     0.25    0.5    0.75    1      2      3
            sigma_min   2.8e-11 0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

        and LINEAR in the distance to the nearest one -- 2.5e-1, 2.6e-2,
        2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way there.  The
        same sweep on a DRIVEN ladder never drops below 0.78: no unit
        multiplier, no singularity, harmonics included.

        ⚠ AND IT IS PHYSICS, NOT CONDITIONING.  A perturbation at a
        harmonic is a perturbation along the orbit, and an oscillator's
        response to that is unbounded phase drift -- there is no bounded
        periodic answer to return.  So this refuses rather than tightening
        a tolerance, and says which quantity to ask for instead.
        """
        if not getattr(pss, 'autonomous', False):
            return
        f0 = 1.0 / float(pss.factored_period().T)
        r = abs(float(freq)) / f0
        d = abs(r - round(r))
        if d <= self.HARMONIC_GUARD:
            raise ValueError(
                'PAC: %s at %.6g Hz is on harmonic %d of this OSCILLATOR\'s '
                'own frequency (%.6g Hz), where I - exp(-j w T) M is '
                'singular -- the unit Floquet multiplier makes it exactly '
                'I - M there. That is physics, not conditioning: a '
                'perturbation along the orbit produces unbounded phase '
                'drift, so there is no bounded periodic response to '
                'return. Ask off-harmonic, or ask for the phase quantity '
                'instead (PSS.ppv()).' % (what, float(freq), round(r), f0))

    @staticmethod
    def _gmres_checked(A, b, rtol, what):
        """GMRES, judged by its RESIDUAL rather than by its status flag.

        ⚠ SCIPY REPORTS BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  These
        operators are `2m x 2m` and often tiny, so the Krylov space is
        exhausted in a handful of steps; the next vector is then numerically
        zero, which is a LUCKY breakdown -- the solution is exact -- and it
        comes back as `info = 4` all the same.  Trusting the flag turns an
        exact answer into a `RuntimeError`, which is what it did for AM/PM
        at small offsets.

        So the residual decides.  A genuine failure still fails, and it
        fails with the residual quoted, because the real cause near a
        harmonic is that the operator is nearly singular there and no
        tolerance will fix it.
        """
        ## ⚠ NOW OUR OWN ARNOLDI-GMRES, which returns the residual as its
        ## verdict instead of a status flag that has to be overruled.
        ## The workaround below survives as the TOLERANCE decision; what
        ## has gone is the second opinion about whether the solve failed.
        n = b.shape[0]
        x, relres, _H, _k = _arnoldi_gmres(
            A.matvec, b, rtol=rtol, maxiter=min(n, 200))
        info = 0
        r = float(np.linalg.norm(b - A.matvec(x)))
        scale = max(float(np.linalg.norm(b)), 1e-300)
        if r / scale > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PAC: %s did not converge (info=%r, relative residual '
                '%.3e). Near a harmonic of the oscillator this operator is '
                'genuinely near-singular and a smaller tolerance will not '
                'help -- move the offset, or ask for the phase quantity.'
                % (what, info, r / scale))
        return x

    def _refuse_coloured(self, pss, what):
        """Refuse a coloured source where the machinery assumes WHITE.

        ⚠ THE TRAP IS THAT NOTHING ELSE WOULD OBJECT. The Lyapunov
        recursion, `diffusion_constant` and eq (22)'s collapse all read
        `CY` at ONE frequency and treat it as the noise intensity at every
        frequency; a coloured source folded that way returns a plausible
        number, not an error (A4d names exactly this shape). Detected by
        evaluating the reduced `CY` at two frequencies -- colour is
        frequency dependence, bias dependence is what `_cy_reduced`
        refuses separately.
        """
        w1 = 2.0 * np.pi / float(pss.period)
        ## ⚠ ONE state, two frequencies: the colour question is separable
        ## from the bias question, and asking it through `_cy_reduced`
        ## refused every MODULATED source before the covariance routes
        ## (which evaluate `CY` per step and handle modulation exactly)
        ## could reach it.
        _xl = np.asarray(pss.factored_period().x_last, dtype=float).ravel()
        c1 = self._cy_at(pss, w1, _xl)
        c2 = self._cy_at(pss, 10.0 * w1, _xl)
        sc = max(float(np.max(np.abs(c1))), 1e-300)
        if float(np.max(np.abs(c1 - c2))) > 1e-9 * sc:
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (its CY '
                'differs between w0 and 10 w0), and this routine assumes '
                'white sources -- it would fold CY at one frequency as if '
                'it held at every frequency and return a plausible wrong '
                'number. Use the frequency-resolved surfaces (pnoise, '
                'phase_psd/coloured_diffusion), or the white-through-filter '
                'form of the source.' % what)

    def _lyapunov_pieces(self, pss, what):
        """The per-step maps, injections and one-period accumulation.

        Returns `(As, Qs, K1, M, m, n)`: the step maps `A_j`, the noise
        injections `Q_j`, the covariance `K1` reached after one period
        starting from zero, the monodromy `M`, and the two widths.

        ⚠ SHARED BY THE DRIVEN AND AUTONOMOUS ROUTES ON PURPOSE.  The two
        differ only in what they do with `I - M kron M`: `covariance`
        inverts it, `oscillator_covariance` borders it because it is
        singular there.  Everything upstream -- the `CY/2` convention, the
        `b = 0` restriction, the `C` ring the forward recursion sees -- is
        one implementation, so the pair cannot drift apart in the way that
        `diffusion_constant` and `covariance` once did over exactly this
        factor of two.
        """
        self._refuse_coloured(pss, what)
        fp = pss.factored_period()
        if fp.kind == 'dirk':
            ## the sequential DIRK per-step map + the SAME exact Van Loan
            ## injection -- see `_lyapunov_pieces_dirk`.
            return self._lyapunov_pieces_dirk(pss, fp, what)
        if fp.kind == 'full':
            ## Radau's own injection: the SAME exact Van Loan integral (it is
            ## the continuous per-step covariance, method-independent), with
            ## the coupled Radau per-step map as A_n -- see
            ## `_lyapunov_pieces_full`.
            return self._lyapunov_pieces_full(pss, fp, what)
        if fp.kind != 'solved_history':
            return self._lyapunov_pieces_plain(pss, fp, what)
        m = pss.cir.n - 1
        n = fp.width
        hs = np.diff(np.asarray(fp.times, dtype=float))
        ## ⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
        ## period.  The Lyapunov accumulation was already per step; hoisting
        ## a single `CY` out of it put a MODULATED source (a switch's
        ## `4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
        ## formulation rather than outside the accuracy, and the
        ## cyclostationarity refusal in `_cy_reduced` then closed the door
        ## on exactly the circuits whose noise is the point (found by the
        ## reference-simulator cross-check, 2026-09-05).  Evaluated at the state
        ## the step's companion was factored at (the implicit step's own
        ## solution); the colour refusal still applies -- colour is a
        ## different axis.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._cy_at(pss, w0,
                                   _W[:, min(k + 1, _W.shape[1] - 1)]))
               for k in range(len(fp.steps))]

        ## the C ring as the forward recursion sees it -- see the replays
        cs0, cs1, ring = [], [], list(fp.opening)
        for _lu, C_new, _a, _b in fp.steps:
            cs0.append(ring[0])
            cs1.append(ring[1])
            ring = [C_new, ring[0]]

        def step_map(k):
            lu, _Cn, alphas, b = fp.steps[k]
            if b:
                raise NotImplementedError(
                    'PAC.%s: derived for a b = 0 companion (Gear-2).' % what)
            A = np.zeros((n, n))
            for j in range(n):
                p0 = np.zeros(m)
                p1 = np.zeros(m)
                (p0 if j < m else p1)[j if j < m else j - m] = 1.0
                A[:m, j] = -lu.solve(alphas[1] * (cs0[k] @ p0)
                                     + alphas[2] * (cs1[k] @ p1))
                A[m:, j] = p0
            return A

        As, Qs = [], []
        for k, (lu, _Cn, _a, _b) in enumerate(fp.steps):
            ## Q = Jf^-1 (CY / 2h) Jf^-T, symmetrised against round-off
            half = cys[k] / (2.0 * hs[k])
            left = np.column_stack([lu.solve(half[:, j]) for j in range(m)])
            Q1 = np.column_stack([lu.solve(left[j, :]) for j in range(m)]).T
            Q = np.zeros((n, n))
            Q[:m, :m] = 0.5 * (Q1 + Q1.T)
            Qs.append(Q)
            As.append(step_map(k))

        K = np.zeros((n, n))
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        M = np.column_stack([fp.matvec(e) for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _lyapunov_pieces_plain(self, pss, fp, what):
        """`_lyapunov_pieces` for the PLAIN path — the one-step companions.

        ⚠ THE PER-STEP STATE DEPENDS ON THE METHOD, and that is the whole
        content of this routine.  `_monodromy_matvec_plain` writes every
        one-step companion as

            S    = a1 C_{k-1} x_{k-1} + b iq_{k-1}
            x_k  = -K S,           K = Jf_k^-1
            iq_k = a0 C_k x_k + S

        **Euler** (`b = 0`): `iq` never re-enters, the state is `x` alone,
        `A_k = -a1 K C_{k-1}` is `m x m`, and nothing downstream changes --
        `n = m`, `M = fp.matvec`, and `ppv()`'s width-`m` vectors border
        it directly.

        **Trapezoidal** (`b = -1`): `iq` DOES re-enter, so the per-step
        state is the PAIR `(x, iq)`, `A_k` is `2m x 2m`, and the noise --
        which enters the KCL rows and reaches `x_k` through `K` -- reaches
        `iq_k` through `a0 C_k K` as well:

            G_k = [ K ; a0 C_k K ],     Q_k = G_k (CY/2h_k) G_k^T

        The period map on that pair is the plain product of the `A_k`
        with NO re-seeding of `iq` at the boundary.  ⚠ THAT IS A DIFFERENT
        OBJECT FROM `fp.matvec`, deliberately: the shooting SOLVE re-seeds
        the companion at each period start (the manufactured opener, B16),
        but the discretised noisy system does not, and the covariance is
        a property of the latter.  The tie between the two is exact and is
        the gate: the pair product applied to `(x, 0)` and read out on `x`
        IS `fp.matvec`.

        ⚠ `oscillator_covariance` IS REFUSED FOR TRAP-PLAIN, with the
        reason: it borders `I - M kron M` with `ppv()`'s null vectors,
        which are width `m` on the plain path, and the pair map is
        `2m x 2m`.  The pair's own null vectors would be needed, with a
        normalisation this record has been burned on twice today
        (`floquet_modes`' state-block scale, `ppv()`'s `v . xdot`).
        Named rather than approximated; euler-plain and gear both work.
        """
        m = pss.cir.n - 1
        hs = np.diff(np.asarray(fp.times, dtype=float))
        ## ⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
        ## period.  The Lyapunov accumulation was already per step; hoisting
        ## a single `CY` out of it put a MODULATED source (a switch's
        ## `4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
        ## formulation rather than outside the accuracy, and the
        ## cyclostationarity refusal in `_cy_reduced` then closed the door
        ## on exactly the circuits whose noise is the point (found by the
        ## reference-simulator cross-check, 2026-09-05).  Evaluated at the state
        ## the step's companion was factored at (the implicit step's own
        ## solution); the colour refusal still applies -- colour is a
        ## different axis.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._cy_at(pss, w0,
                                   _W[:, min(k + 1, _W.shape[1] - 1)]))
               for k in range(len(fp.steps))]
        C_open = np.asarray(fp.opening[0], dtype=float)
        prevC = [C_open] + [np.asarray(st[1], dtype=float)
                            for st in fp.steps[:-1]]
        bs = {bool(st[3]) for st in fp.steps}
        if len(bs) != 1:
            raise NotImplementedError(
                'PAC.%s: the plain period mixes b = 0 and b != 0 steps, '
                'which have different per-step states.' % what)
        pair = bs.pop()
        if pair and what == 'oscillator_covariance':
            raise NotImplementedError(
                'PAC.oscillator_covariance: the trapezoidal plain path\'s '
                'per-step state is the pair (x, iq), so its period map is '
                '2m x 2m, and the bordered solve needs THAT map\'s null '
                "vectors -- ppv()'s are width m. Not built (see "
                '_lyapunov_pieces_plain). Use method=\'euler\' for a plain '
                "width-m reference, or method='gear'.")
        n = 2 * m if pair else m
        As, Qs = [], []
        for k, (lu, C_new, alphas, b) in enumerate(fp.steps):
            Ck = np.asarray(C_new, dtype=float)
            Cp = prevC[k]
            A = np.zeros((n, n))
            for j in range(n):
                e = np.zeros(n)
                e[j] = 1.0
                p0, p1 = e[:m], (e[m:] if pair else None)
                S = alphas[1] * (Cp @ p0)
                if pair:
                    S = S + b * p1
                x = -np.asarray(lu.solve(S), dtype=float)
                A[:m, j] = x
                if pair:
                    A[m:, j] = alphas[0] * (Ck @ x) + S
            ## noise: K (CY/2h) K^T on the state block, built the same way
            ## as the solved-history route so the two cannot drift apart
            half = cys[k] / (2.0 * hs[k])
            left = np.column_stack([lu.solve(half[:, j]) for j in range(m)])
            Q1 = np.column_stack([lu.solve(left[j, :]) for j in range(m)]).T
            Q1 = 0.5 * (Q1 + Q1.T)
            Q = np.zeros((n, n))
            Q[:m, :m] = Q1
            if pair:
                Bk = alphas[0] * Ck
                Q[:m, m:] = Q1 @ Bk.T
                Q[m:, :m] = Bk @ Q1
                Q[m:, m:] = Bk @ Q1 @ Bk.T
                Q = 0.5 * (Q + Q.T)
            As.append(A)
            Qs.append(Q)
        K = np.zeros((n, n))
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
        if pair:
            ## ⚠ THE PERIOD MAP RE-SEEDS THE COMPANION, AND THAT IS
            ## LOAD-BEARING.  The plain product of the A_k carries `iq`
            ## across the boundary, and its `I - M kron M` is SINGULAR:
            ## trapezoidal maps an algebraic row's companion by exactly -1
            ## per step, so the un-reset pair carries a marginal mode --
            ## the `(-1)^n` obstruction this file records for every
            ## formulation that keeps `iq` across a period (measured here:
            ## LinAlgError on a driven RLC).  The shooting solve is
            ## well-posed because the manufactured opener re-seeds `iq` at
            ## zero; the covariance's period map must do the same.  With
            ## `iq` zeroed at the start, the x->x block of the product IS
            ## `fp.matvec` (tied to 1e-12 above), and the map on the pair
            ## is the product applied to `(x, 0)`.
            Mp = np.eye(n)
            for A in As:
                Mp = A @ Mp
            M = np.zeros((n, n))
            M[:, :m] = Mp[:, :m]
        else:
            M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                                 for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _vanloan_step_injection(self, Cr, Gr, CYr, h):
        """The per-step process-noise covariance `Q_n` for TR-BDF2, by the
        DAE-projected VAN LOAN integral.

        For ADDITIVE (linearised) noise the injection is the DETERMINISTIC
        integral `Q = integral_0^h Phi(h,s) D Phi(h,s)^T ds` -- the Levy
        areas vanish, so there are no stochastic stage weights to derive
        (Roemisch & Winkler; confirmed by the naive two-stage scheme coming
        out 27% biased on kT/C).  Van Loan evaluates it exactly: the
        upper-right block of `expm([[-A, D],[0, A^T]] h)` premultiplied by
        the flow.

        ⚠ BUT MNA IS A DAE (`C` singular), and the nilpotent block
        DIFFERENTIATES white noise -- discretised white noise has variance
        `S/h`, so a covariance formed on an algebraic row diverges as `1/h`
        (measured).  So Van Loan is applied on the DIFFERENTIAL SUBSPACE
        only (the capacitive nodes -- Demir 1996 propagates exactly there),
        after eliminating the algebraic variables by their Schur complement.
        The algebraic noise is routed to the differential rows through the
        same elimination (`R_proj`), so a source with a capacitive path
        (Winkler's `im A_N subset im A_C`) is handled; a source on a bare
        constraint has no differential image and is dropped rather than
        divergently amplified -- the projection is structurally immune to
        the `1/h` blow-up.

        Verified against kT/C at second order on R||C (ODE) and VS-R-C
        (DAE), matching an independent measurement to the digit; the
        stationary error is the METHOD's O(h^2), NOT machine zero (a
        machine-zero kT/C would mean a method-consistent `Q = P(1-A^2)`
        fudge that corrupts the transient covariance).
        """
        import scipy.linalg as sla
        Cr = np.asarray(Cr, dtype=float)
        Gr = np.asarray(Gr, dtype=float)
        CYr = np.asarray(np.real(CYr), dtype=float)
        m = Cr.shape[0]
        d = [i for i in range(m)
             if np.any(np.abs(Cr[i, :]) > 0) or np.any(np.abs(Cr[:, i]) > 0)]
        a = [i for i in range(m) if i not in d]
        if not d:
            raise NotImplementedError(
                'PAC: this circuit has no capacitive (differential) node, so '
                'there is no covariance to propagate -- every state is '
                'algebraic and a white source on it is differentiated by the '
                'DAE. Add the capacitance that shunts the noise, or ask for a '
                'quantity that does not need a covariance.')
        di = np.ix_(d, d)
        Emb = np.zeros((m, len(d)))
        for k, i in enumerate(d):
            Emb[i, k] = 1.0
        if a:
            Gaa = Gr[np.ix_(a, a)]
            Gai = np.linalg.inv(Gaa)
            Gad = Gr[np.ix_(a, d)]
            Sc = Gr[di] - Gr[np.ix_(d, a)] @ Gai @ Gad
            ## R_proj = [I_d, -G_da G_aa^-1] routes the algebraic-row noise
            ## into the differential rows through the same elimination
            Rproj = np.zeros((len(d), m))
            for k, i in enumerate(d):
                Rproj[k, i] = 1.0
            Rproj[:, a] = -Gr[np.ix_(d, a)] @ Gai
            CYred = Rproj @ CYr @ Rproj.T
            ## the algebraic variables are slaved to the differential ones
            Emb[np.ix_(a, range(len(d)))] = -Gai @ Gad
        else:
            Sc = Gr[di]
            CYred = CYr[di]
        Cdd = Cr[di]
        Cinv = np.linalg.inv(Cdd)
        Ared = -Cinv @ Sc
        ## CY is a ONE-SIDED density; CY/2 is the two-sided intensity, the
        ## same convention `_lyapunov_pieces` and `diffusion_constant` use
        Dred = Cinv @ (0.5 * CYred) @ Cinv.T
        Dred = 0.5 * (Dred + Dred.T)
        md = len(d)
        Z = np.zeros((md, md))
        E = sla.expm(np.block([[-Ared, Dred], [Z, Ared.T]]) * float(h))
        Phi = E[md:, md:].T
        Qd = Phi @ E[:md, md:]
        Qd = 0.5 * (Qd + Qd.T)
        return Emb @ Qd @ Emb.T

    def _lyapunov_pieces_full(self, pss, fp, what):
        """`_lyapunov_pieces` for a Radau IIA(3) Floquet source.

        Identical in structure to `_lyapunov_pieces_trbdf2`: the per-step
        injection `Q_n` is the DAE-projected VAN LOAN integral at that step's
        operating point (`_vanloan_step_injection`), which is the EXACT
        continuous per-step covariance and so is the same object for every
        integrator; only the per-step transition `A_n` differs -- here the
        coupled Radau step map (dense, `m x m`, via
        `_monodromy_matvec_full` one step at a time).  State width `m`, so
        `n = m`.

        ⚠ THE INJECTION IS EXACT, THE METHOD SETS ONLY THE PROPAGATION.  The
        covariance still converges to the stationary target (kT/C on an RC)
        at the injection's O(h^2), not at Radau's O(h^5): Van Loan already
        integrates the step exactly, so refining the grid gains on the
        recursion's discretisation of a continuous Lyapunov flow, which the
        higher-order transition does not change.
        """
        self._refuse_coloured(pss, what)
        m = pss.cir.n - 1
        n = m
        hs = np.diff(np.asarray(fp.times, dtype=float))
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        As, Qs = [], []
        for k, step in enumerate(fp.steps):
            xk = _W[:, min(k + 1, _W.shape[1] - 1)]
            Cn = np.asarray(pss._C_at(xk), dtype=float)
            Gn = np.asarray(pss._G_at(xk), dtype=float)
            CYn = self._cy_at(pss, w0, xk)
            A_k = np.column_stack([
                np.asarray(pss._monodromy_matvec_full([step], e), dtype=float)
                for e in np.eye(m)])
            As.append(A_k)
            Qs.append(self._vanloan_step_injection(Cn, Gn, CYn, hs[k]))
        K = np.zeros((n, n))
        for A_k, Q_k in zip(As, Qs):
            K = A_k @ K @ A_k.T + Q_k
        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _lyapunov_pieces_dirk(self, pss, fp, what):
        """`_lyapunov_pieces` for a lower-triangular (DIRK/ESDIRK) Floquet
        source -- identical to `_lyapunov_pieces_full` except the per-step
        transition `A_n` is the sequential DIRK step map
        (`_monodromy_matvec_dirk`).  The Van Loan injection `Q_n` is the same
        exact continuous per-step covariance (method-independent)."""
        self._refuse_coloured(pss, what)
        m = pss.cir.n - 1
        n = m
        hs = np.diff(np.asarray(fp.times, dtype=float))
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        As, Qs = [], []
        for k, step in enumerate(fp.steps):
            xk = _W[:, min(k + 1, _W.shape[1] - 1)]
            Cn = np.asarray(pss._C_at(xk), dtype=float)
            Gn = np.asarray(pss._G_at(xk), dtype=float)
            CYn = self._cy_at(pss, w0, xk)
            A_k = np.column_stack([
                np.asarray(pss._monodromy_matvec_dirk([step], e), dtype=float)
                for e in np.eye(m)])
            As.append(A_k)
            Qs.append(self._vanloan_step_injection(Cn, Gn, CYn, hs[k]))
        K = np.zeros((n, n))
        for A_k, Q_k in zip(As, Qs):
            K = A_k @ K @ A_k.T + Q_k
        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def covariance(self, pss, samples=False):
        """The periodic (cyclostationary) state covariance — DRIVEN circuits.

        ⚠ A GRID CHOSEN FOR `kT/C` IS NOT A GRID FOR THE PROFILE.  With
        `CY` per step (2026-09-05) a switched capacitor's HELD variance
        reads `kT/C` to 1e-4 at 1600 points and converges at better than
        second order, while the TRACKING phase sits at this routine's
        O(h/tau) floor -- 4% out at 800 points where the held value is
        already 1.6e-4 -- and both agree with a reference simulator's sampled pnoise at
        matched instants to 1e-3 (0.99878 track, 0.99915 edge, 0.99999
        hold).  The tracked variance is 0.957 kT/C, NOT kT/C: a sinusoidal
        clock holds the switch at full `gon` only instantaneously, so the
        capacitor is never in equilibrium with `Ron`; both tools agree on
        that independently.

        Returns `K0`, the covariance at `t = 0`; with `samples=True`,
        `(K0, [K_j])`, the covariance at every step, which is the
        time-varying statistic this exists to produce.

        The noise covariance obeys a Lyapunov recursion alongside the
        trajectory, `K_{j+1} = A_j K_j A_jᵀ + Q_j`, so over one period
        `K_N = M K_0 Mᵀ + K_1`.  Periodicity closes it:

            (I - M ⊗ M) vec(K_0) = vec(K_1)

        ⚠ ONE LINEAR SOLVE, NO NEWTON.  The Lyapunov equation is LINEAR in
        `K`, so shooting on it is exact in a single step — unlike the
        trajectory it rides on.  The monodromy of the covariance system is
        the KRONECKER SQUARE of the circuit's, so its multipliers are the
        pairwise products `lambda_i lambda_j`.

        ⚠ AND THAT IS WHY IT REFUSES AN OSCILLATOR.  There `lambda_1 = 1`
        gives `lambda_1^2 = 1`, so `I - M ⊗ M` is exactly as singular as
        `I - M` — measured 3.1e-11 against 3.8e-11 — and the covariance
        does not settle, it GROWS.  Variance linear in `t` is a random
        walk, which is phase diffusion, which is the linewidth.  Demir 2002
        gives the physical counterpart: an oscillator's output noise is
        STATIONARY, not cyclostationary, because "noisy autonomous systems
        cannot provide a perfect time reference".  There is no
        cyclostationary object there to compute, and `oscillator_spectrum`
        is the right route instead.

        ⚠ `CY/2` IS THE ONE-SIDED-TO-TWO-SIDED CONVERSION AND IT IS NOT
        COSMETIC.  `CY` is a one-sided density (a resistor's `4kT/R`), so
        the per-step injection is `Q_j = Jf_j^-1 (CY_j / 2h_j) Jf_j^-T`.
        MEASURED against `kT/C` — exact, famously independent of `R` — on
        an RC circuit, with and without the half:

            npts        100      200      400      800
            CY          1.861    1.928    1.963    1.981
            CY/2        0.931    0.964    0.982    0.991

        The full-`CY` column converges to 2 and the halved one to 1, so the
        factor is settled by the measurement rather than by argument.  The
        residual halves per grid doubling — O(h), first order, which a
        piecewise-constant approximation to white noise is.

        ⚠ AND THE GRID MUST RESOLVE THE NOISE BANDWIDTH, which is a real
        precondition rather than an accuracy note.  The first attempt at
        that gate read 0.517 because the RC pole at 159 kHz sat ABOVE the
        grid's 100 kHz Nyquist: the discrete system genuinely does not
        carry the noise the continuous one does.  A `kT/C` that comes back
        low is the grid, not the code.

        ⚠ COST: the solve has `(2m)^2` unknowns and is dense here, so it is
        `O(m^4)`.  Small circuits only until that is replaced.
        """
        self._check_circuit(pss)
        if getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.covariance: an OSCILLATOR has no periodic covariance. '
                'Its unit multiplier squares to one, so I - M kron M is '
                'singular and the covariance grows without bound rather '
                'than settling -- that growth IS the phase diffusion, and '
                'its output noise is stationary rather than '
                'cyclostationary. Use oscillator_covariance() for the '
                'split into a bounded orbital part and that growth, or '
                'oscillator_spectrum() for the lineshape it produces.')
        ## the source-injection surfaces use the gear twin when the Floquet
        ## source is TR-BDF2 (its two-stage Q_j is not built) -- see
        ## `_lyapunov_host`
        pss = pss._lyapunov_host()
        As, Qs, K1, M, m, n = self._lyapunov_pieces(pss, 'covariance')
        S = np.eye(n * n) - np.kron(M, M)
        K0 = np.linalg.solve(S, K1.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        if not samples:
            return K0
        seq, K = [K0], K0
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
            seq.append(0.5 * (K + K.T))
        return K0, seq

    ## ⚠⚠⚠ THE NOTE THAT WAS HERE ACCUSED THIS ROUTE AND WAS WRONG.  A
    ## MONTE CARLO SETTLED IT THE OTHER WAY: this route is CORRECT and
    ## `orbital_correlation` is the one that fails on an asymmetric orbit.
    ## Direct SDE simulation of the variational system (trapezoidal, the
    ## calibrated `Var(i) = CY/(2h)` injection, phase projected out every
    ## step), sharing no Lyapunov solve and no modal sum:
    ##
    ##     a      |R| MONTE CARLO  |R| modal    |Lyap| proj   MC/Lyap
    ##     0.00   4.3629e-06       4.4515e-06   4.4437e-06    0.9818
    ##     0.30   2.9992e-04       3.6913e-06   2.9986e-04    1.0002
    ##
    ## `a = 0` is the CONTROL -- both routes agree there, so the MC had a
    ## known answer to hit, and it did (2 %).  At `a = 0.30` it lands on this
    ## route to 0.02 % and is 81x from the modal one.
    ##
    ## ⚠⚠ THE ARGUMENT THAT MISLED ME, RECORDED BECAUSE IT WAS PLAUSIBLE:
    ## `|lam2|` FALLS with asymmetry, so relaxation gets FASTER, so the
    ## transverse variance "should" shrink -- and the modal route did shrink
    ## while this one grew 67x.  That reasoning is WRONG: asymmetry changes
    ## the MODE SHAPES, so the noise projected onto the orbital direction
    ## grows, and the variance rises DESPITE the faster relaxation.  A
    ## physical argument is not a measurement.
    ##
    ## The original (refuted) note follows for the record:
    ## ⚠ SUPERSEDED: THIS ROUTE DEPARTS FROM PHYSICS
    ## ON AN ASYMMETRIC ORBIT, AND `orbital_correlation` DOES NOT.  van der
    ## Pol + `a u^2`, sweeping `a` (orbit asymmetry 0 -> 0.41):
    ##
    ##     a      |R| modal    |Lyap| proj  |K_orb raw|  |lam2|    amp
    ##     0.00   4.4515e-06   4.4437e-06   5.0034e-06   0.882521  2.000
    ##     0.05   4.4702e-06   4.7152e-06   4.7471e-06   0.881719  2.004
    ##     0.15   4.9505e-06   1.7371e-05   2.7846e-05   0.874886  2.034
    ##     0.30   3.6913e-06   2.9986e-04   9.5944e-04   0.844322  2.164
    ##
    ## `|lam2|` FALLS (0.883 -> 0.844), so amplitude relaxation gets FASTER
    ## and the transverse variance should get slightly SMALLER.  The modal
    ## sum does exactly that (4.45e-06 -> 3.69e-06).  This route grows 67x
    ## projected and 192x raw.  ⚠ The RAW covariance grows MORE than the
    ## projected one, so it is not the oblique projection -- it is this
    ## covariance.
    ##
    ## ⚠ FOUR EXPLANATIONS EXCLUDED BY MEASUREMENT, not by argument:
    ##   * harmonic truncation -- the disagreement is FLAT at 9.877e-01 from
    ##     `H = 4` to `H = 128`;
    ##   * the projection's tangent proxy -- `|cos(u, tangent)| = 1.000000`
    ##     at every asymmetry, against the Floquet phase mode;
    ##   * the modal decomposition -- `|lam1| = 1.000000`, `lam2` real and
    ##     well separated, one orbital mode, `p(T)-p(0) ~ 1e-14`,
    ##     `q^T C p = 1.000000`;
    ##   * a defect in `orbital_correlation` -- its two internal routes agree
    ##     to 3.4e-04 independently of `a`.
    ##
    ## ⚠⚠ WHAT THIS DOES **NOT** INVALIDATE.  Every use of this function as a
    ## reference in this file was on a SYMMETRIC orbit, where the two routes
    ## agree to 0.3 % -- including A9 step 3's three-way gate and the C^2
    ## biorthonormalisation defect it caught on 2026-09-07.  Those stand.
    ## ⚠ WHAT IS OPEN: which route is right is NOT settled.  The physical
    ## argument favours the modal one, but that is an argument, and this file
    ## does not close items on arguments.  A transient MONTE CARLO of the
    ## orbital fluctuation is the decisive third route and has not been run.
    ## Until then, treat this on a strongly asymmetric orbit as unvalidated.
    def oscillator_covariance(self, pss, samples=False):
        """The state covariance of a FREE-RUNNING oscillator, split in two.

        Returns `(K_orb, d, info)`.  `K_orb` is the BOUNDED periodic
        (orbital) part of the covariance at `t = 0`; `d` is the growth per
        period along the orbit tangent, so

        ⚠ "BOUNDED" IS NOT "TRANSVERSE".  `K_orb` has the SECULAR growth
        removed and still contains the phase direction's bounded
        within-period variance.  Demir's orbital deviation `y` is the
        OBLIQUE projection `v_1^T y = 0`, so the transverse covariance is
        `Pi K_orb Pi^T` with `Pi = I - u v^T/(v^T u)` -- which is what
        `orbital_correlation`'s eq (23) sum equals (to 1e-4), and what
        `K_orb` itself does NOT equal (2-6 %, falling as 1/Q).  Read
        `K_orb` as the bounded part; project it if you want `R_yy(0)`.
        See `orbital_correlation`.

            K(t_0 + n T) = K_orb + n d u u^T

        exactly, for every integer `n`, with `u` the pair-space tangent
        scaled so its first block is `xdot(0)`.

        ⚠ WITH `samples=True` THE SPLIT MOVES WITH THE ORBIT, AND IT IS
        WORTH SAYING PRECISELY BECAUSE THE OBVIOUS READING IS WRONG.
        `info['orbital_samples'][j]` is `P(t_j)`, the solution started from
        `K_orb` at `t = 0`, and it satisfies

            K(t_j + n T) = P(t_j) + n d u_j u_j^T,   u_j = Phi(t_j, 0) u

        so `P` is periodic UP TO the growth -- `P(T) = P(0) + d u u^T`, not
        `P(T) = P(0)`.  The walk is along the orbit and the orbit turns, so
        the growth DIRECTION is the propagated tangent rather than a fixed
        `u`.  MEASURED: `P(T) - P(0)` matches `d u u^T` to 3.2e-15, and the
        full prediction holds to 2.5e-09 against a brute-force recursion
        run forty periods (9,600 steps) from `K = 0`.

        ⚠ THIS IS THE OBJECT `covariance` REFUSES TO RETURN, and the
        refusal was right: there is no periodic solution, so anything that
        returned one number would be hiding the physics.  `lambda_1 = 1`
        gives `lambda_1^2 = 1`, so `I - M kron M` is exactly singular --
        MEASURED here at `sigma_min` 2.3e-11 with the next singular value
        at 0.997, i.e. a null space that is cleanly ONE-DIMENSIONAL and
        spanned by `u kron u`, with left null `v kron v`.  So it borders
        exactly as the PPV and the deflated PAC solve do, and the border is
        the pair the rest of this class already computes.

        ⚠ THE SPLIT IS NOT A NUMERICAL DEVICE, IT IS THE ANSWER.  Demir
        2002: an oscillator's noise is STATIONARY, not cyclostationary,
        because "noisy autonomous systems cannot provide a perfect time
        reference".  `K_orb` is the part a designer can read as an
        amplitude/orbital noise -- it settles, it is periodic, it is
        finite.  `n d u u^T` is the random walk ALONG the orbit, which
        never settles and which no periodic object can hold.  Reporting
        only their sum at some finite time is what makes an oscillator
        covariance look divergent and useless; reporting the parts makes
        both usable.

            [ I - M kron M    u kron u ] [ vec(K_orb) ]   [ vec(K_1) ]
            [ (v kron v)^T        0    ] [     d      ] = [     0    ]

        ⚠ AND `d` HAS A CLOSED FORM THAT NEEDS NO KRONECKER AT ALL.
        Left-multiplying the first row by `(v kron v)^T` kills the singular
        block, leaving

            d = (v^T K_1 v) / (v . u)^2

        an `O(n^2)` contraction against the `n^4` solve.  Both are computed
        and `info['d_residual']` is their relative difference; they are the
        same quantity by construction, so a disagreement means the border
        pair is wrong rather than that one route is less accurate.

        ⚠ `(v . u)` IS NOT 1 AND ASSUMING IT IS COSTS A FACTOR OF 2.3.
        `ppv()` normalises on the FIRST BLOCK, `v[:m] . xdot = 1`, which is
        the normalisation a state perturbation entering the first block
        sees -- an injected current, and what every other shipped path
        does.  The FULL PAIR contraction is a different number: 0.663 on
        van der Pol, so `(v . u)^-2 = 2.28`.  That exact mistake produced a
        2.31x discrepancy that was chased as a code defect for a while; it
        is why `d` is written with the pair inner product spelled out.

        ⚠ `d` ALONE IS MEANINGLESS WITHOUT PINNING `u`'s SCALE.  Rescaling
        `u -> s u` sends `d -> d / s^2`, so only the PRODUCT `d u u^T` --
        returned as `info['growth']` -- is an invariant of the circuit.
        `u` is pinned here by `C u = q`, the same condition `ppv()` uses to
        scale the tangent, which makes its first block exactly `xdot(0)`
        and gives `d` its physical reading below.

        ⚠ WHICH MAKES `d / T` A COMPLETELY INDEPENDENT ROUTE TO THE
        DIFFUSION CONSTANT, and that is this method's real gate.  A phase
        deviation `alpha` displaces the state by `alpha u`, so the growing
        covariance is `Var(alpha) u u^T = c t u u^T`, giving `d = c T`.
        The two computations share only the `CY/2` convention: `c` is a
        quadratic form in the ADJOINT-replayed PPV, while `d` comes from a
        FORWARD Lyapunov recursion closed by a bordered Kronecker solve.
        Neither touches the other's machinery.

        ⚠ AND THE TWO ANCHORS BEHIND THEM ARE ALSO INDEPENDENT, which is
        the property that was missing when a 2x error survived a 0.9965
        agreement.  `covariance`'s injection is anchored to `kT/C`;
        `diffusion_constant` is anchored to a nonlinear Monte Carlo reading
        phase from zero crossings.  `info['c_from_growth']` against
        `diffusion_constant` therefore closes a loop between two separately
        anchored quantities rather than reproducing one of them.

        ⚠ COST: the bordered solve has `(2m)^2 + 1` unknowns and is dense,
        so it is `O(m^4)` like `covariance`.  Small circuits only.  The
        closed form for `d` is cheap; pass `samples=False` and read
        `info['c_from_growth']` if the orbital part is not wanted.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.oscillator_covariance: this splits a covariance that '
                'GROWS into a bounded part plus a random walk along the '
                'orbit. A driven circuit has neither -- its covariance '
                'settles, and I - M kron M is nonsingular. Use '
                'covariance().')
        ## the source-injection surfaces use the gear twin when the Floquet
        ## source is TR-BDF2 (its two-stage Q_j is not built).  Swapped BEFORE
        ## both the Lyapunov pieces and `ppv` below, so the bordering keeps
        ## them on one host -- see `_lyapunov_host`.
        pss = pss._lyapunov_host()
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'oscillator_covariance')

        v, pinfo = pss.ppv()
        v = np.asarray(v, dtype=float).ravel()
        u = np.asarray(pinfo['tangent_pair'], dtype=float).ravel()
        xdot = np.asarray(pinfo['xdot'], dtype=float).ravel()
        ## rescale the bordered solve's DIRECTION onto the tangent `ppv`
        ## already scaled by `C u = q`; least squares so this is stable
        ## even where `u[:m]` is small, and exact where it is not.
        uu = float(u[:m] @ u[:m])
        if uu == 0.0:
            raise ValueError(
                'PAC.oscillator_covariance: the tangent has no first '
                'block, so its scale cannot be pinned to xdot(0).')
        u = u * (float(u[:m] @ xdot) / uu)

        vu = float(v @ u)
        if vu == 0.0:
            raise ValueError(
                'PAC.oscillator_covariance: the left and right null '
                'directions are orthogonal in the pair space, so the '
                'bordered system is singular. That should not happen on a '
                'converged limit cycle.')
        d_closed = float(v @ K1 @ v) / (vu * vu)

        S = np.eye(n * n) - np.kron(M, M)
        uk = np.kron(u, u)
        vk = np.kron(v, v)
        B = np.zeros((n * n + 1, n * n + 1))
        B[:n * n, :n * n] = S
        B[:n * n, n * n] = uk
        B[n * n, :n * n] = vk
        rhs = np.concatenate((K1.reshape(-1), [0.0]))
        z = np.linalg.solve(B, rhs)
        K_orb = z[:n * n].reshape(n, n)
        K_orb = 0.5 * (K_orb + K_orb.T)
        d = float(z[n * n])

        scale = max(abs(d), abs(d_closed), 1e-300)
        info = {'d_closed_form': d_closed,
                'd_residual': abs(d - d_closed) / scale,
                'growth': d * np.outer(u, u),
                'c_from_growth': d / float(pss.period),
                'tangent_pair': u,
                'ppv_pair': v,
                'pair_inner': vu,
                'sigma_min': float(np.linalg.svd(S, compute_uv=False)[-1]),
                'sigma_min_bordered':
                    float(np.linalg.svd(B, compute_uv=False)[-1]),
                'null_residual': float(np.linalg.norm(S @ uk))
                                 / max(float(np.linalg.norm(uk)), 1e-300),
                'ppv': pinfo}
        if samples:
            ## ⚠ THE GROWTH DIRECTION MOVES WITH THE ORBIT.  The invariant
            ## is `K(t_j + nT) = K_orb(t_j) + n d u_j u_j^T` with `u_j` the
            ## FORWARD-propagated tangent, not `u` held fixed -- the walk
            ## is along the orbit, and the orbit turns.
            orb, grw, K, uj = [K_orb], [d * np.outer(u, u)], K_orb, u
            for A, Q in zip(As, Qs):
                K = A @ K @ A.T + Q
                uj = A @ uj
                orb.append(0.5 * (K + K.T))
                grw.append(d * np.outer(uj, uj))
            info['orbital_samples'] = orb
            info['growth_samples'] = grw
            info['times'] = np.asarray(pss.factored_period().times,
                                       dtype=float)
        return K_orb, d, info

    def orbital_mode_weights(self, pss, nmodes=None):
        """`K_orb` resolved onto the Floquet modes — A9's second step.

        ⚠⚠⚠ READ THIS FIRST: THE BASIS OMITS THE ANNIHILATED MODES, AND WHAT
        THEY CARRY IS A FLOOR NOTHING BELOW CAN GO UNDER.  `floquet_modes`
        returns the NON-NULL directions, so `sum cw[k,k'] u_k u_k'^H` reproduces
        only the part of `K_orb` that lives on them.  How much that is depends
        entirely on WHERE THE NOISE ENTERS -- measured on `_osc_with_ladder`'s
        circuit at `nslow = 4`, moving one current source and changing nothing
        else::

            injected at            ||K_orb||    reconstruction residual
            the oscillator node    2.70e-05     1.80e-03   (0.18%)
            a SLOW ladder node     3.94e-01     3.56e-01   (36%)
            a FAST ladder node     6.87e+02     9.996e-01  (99.96%)
            a faster one           3.33e+03     9.999e-01  (99.99%)

        **When the injection lands in a fast branch the non-null modes capture
        essentially NOTHING of the covariance.**  The annihilated modes are
        killed by the period map, so they enter the stationary covariance only
        through the `j = 0` term -- but that term is not small when the noise
        is injected there, and THAT IS WHERE DEVICE NOISE ACTUALLY IS: every
        resistor in a bias or tuning network.

        ⚠ SO A MODAL ORBITAL SPECTRUM BUILT ON THIS BASIS IS COMPLETE ONLY FOR
        NOISE THAT ENTERS THE SLOW SUBSPACE, and the suite's own gate on this
        (`rel < 1e-2`) holds because its fixture injects at the oscillator
        node.  That is a property of the fixture, not of the method.

        ⚠ AND THE RESIDUAL IS A DETECTOR, NOT A TRUNCATION BOUND.  It catches a
        DROPPED NON-NULL MODE well -- which is what the note below claims for
        it -- but it SATURATES at the floor above, so it cannot certify a
        truncation below whatever the null modes carry, however many modes are
        kept.  Independently reproduced by a peer session on a different
        oscillator with a different `K_orb` route (69% there, mechanism
        identical, magnitude not transferable).

        Returns `(cw, modes, K_orb)` with `cw[k, k'] = v_k† K_orb v_k'`,
        the weight of each pair of Floquet directions in the bounded
        (orbital) part of the state covariance.

        ⚠ **THIS IS THE BRIDGE BETWEEN THE TWO ROUTES WE ALREADY OWN.**
        `oscillator_covariance` gets `K_orb` from a bordered Kronecker
        solve; `floquet_modes` gets the eigen-directions from the
        monodromy. Traversa & Bonani's eq (22) sums over exactly these
        mode pairs, so resolving the covariance we already trust onto the
        modes is the step that connects them — and, unlike the spectrum
        itself, it has an **exact identity** to check against:

            Σ_{k,k'} cw[k,k'] · u_k u_{k'}†  =  K_orb

        because `(u, v)` are biorthonormal. A wrong pairing, a wrong
        normalisation, or a dropped mode all break that reconstruction
        while leaving every individual eigenvector residual clean.

        ⚠ **THE PHASE MODE IS INCLUDED AND ITS WEIGHT SHOULD BE SMALL, NOT
        ZERO.** `K_orb` is the part of the covariance that stays bounded,
        with the along-orbit growth `n·d·uuᵀ` already removed — so the
        `k = k' = 1` entry is what the split left behind rather than a
        quantity that must vanish. Reading it as an error is a
        misinterpretation of `oscillator_covariance`'s own contract.

        ⚠ **NOT THE SPECTRUM.** `S_yy` additionally needs the Fourier
        coefficients of the periodic parts (`floquet_modes` returns them
        as `p`/`q`) and the resolvent `1/(i(j−j')ω₀ − μ_l' − μ_l*)` of
        eq (22), and the OUTPUT spectrum then needs a layer that can carry
        contributions asymmetric about the carrier. Those are not built.
        """
        K_orb, _d, _info = self.oscillator_covariance(pss)
        K = np.asarray(K_orb, dtype=float)
        n = K.shape[0]
        modes = pss.floquet_modes(pss, nmodes=(n if nmodes is None
                                               else int(nmodes)))
        V = np.column_stack([m['v0'] for m in modes])
        cw = V.conj().T @ K @ V
        return cw, modes, K

    ORBITAL_HARMONICS = 32

    ## Half-wave asymmetry above which `orbital_correlation` is known to be
    ## wrong.  MEASURED (below); 0.02 is a decade inside the smallest
    ## asymmetry at which the error was already visible.
    ORBITAL_ASYMMETRY_LIMIT = 0.02

    def _orbit_asymmetry(self, pss):
        """Half-wave asymmetry of the orbit, in [0, ~1].

        `max|x(t) + x(t + T/2)| / max|x|` on the first state row -- zero for a
        half-wave symmetric orbit (van der Pol), growing as the orbit
        distorts.  Cheap: the waveform is already stored.
        """
        W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        if W.size == 0 or W.shape[1] < 4:
            return 0.0
        row = W[0]
        half = len(row) // 2
        den = float(np.max(np.abs(row)))
        if den <= 0.0:
            return 0.0
        return float(np.max(np.abs(row[:half] + row[half:2 * half]))) / den

    def _warn_if_orbit_is_asymmetric(self, pss):
        """⚠ On a strongly asymmetric orbit the modal sum carries an O(h)
        discretisation residual that the symmetric fixtures never show.

        ⚠⚠ THIS WARNING WAS WRITTEN FOR A DEFECT THAT IS NOW FIXED, and kept
        for the residual.  On 2026-09-07 `orbital_correlation` read 81x LOW
        against a Monte Carlo on van der Pol + `a u^2` at asymmetry 0.41.
        The cause was in `floquet_modes`: the replayed adjoint is `C^T q`,
        not `q`, and was used untransformed -- invisible on a unit-reactance
        symmetric orbit, catastrophic off-axis where the two adjoints are
        nearly parallel.  With the `C^-T` transform applied, against the
        same Monte-Carlo-validated Lyapunov reference at `a = 0.30`:

            npts    eq22 / Lyapunov
             400       1.0595
             800       1.0300
            1600       1.0151

        halving per doubling -- an `O(h)` DISCRETISATION residual of the
        adjoint replay, converging to 1, not a defect.  At `a = 0` it is
        1.0001.  So this warns that the residual is grid-limited on such an
        orbit and says how to shrink it; it no longer says the answer is
        wrong, because it is not.
        """
        try:
            asym = self._orbit_asymmetry(pss)
        except Exception:
            return
        if asym <= self.ORBITAL_ASYMMETRY_LIMIT:
            return
        warnings.warn(
            'PAC.orbital_correlation: this orbit has half-wave asymmetry '
            '%.3f. On such an orbit the modal sum carries an O(h) '
            'discretisation residual of the adjoint replay that symmetric '
            'orbits do not show -- measured 6%% high at 400 points per period '
            'and halving per doubling against a Monte-Carlo-validated '
            'reference. Refine the grid to tighten it, or use '
            'PAC.oscillator_covariance (Lyapunov) for the covariance alone.'
            % (asym,),
            RuntimeWarning, stacklevel=3)

    def orbital_correlation(self, pss, H=None):
        """`R_yy(0)` and the `C_lhj` of Traversa & Bonani eq (22) — A9 step 3.

        Returns `(R, C)`.  `R` is the STATIONARY transverse (orbital)
        state covariance, `m x m` real symmetric — eq (23),
        `R = Σ_{l≥2,h,j} C_lhj`.  `C` maps `(l, h, j)` to the `m x m`
        complex coefficient, over every non-null orbital mode `l ≥ 2` and
        harmonics `|h|, |j|, |j'| ≤ H`.  `H` defaults to
        `ORBITAL_HARMONICS`; van der Pol converges by `H = 4`, and a
        strongly non-sinusoidal orbit needs more — check by raising it.

        ⚠ STATIONARY WHITE SOURCES ONLY, and that is what makes it
        computable without `B`.  Eq (22) needs the Fourier coefficients
        of `v_l(t)^T B(t)`; with `CY = B B^T` constant those products
        collapse to `V~_{l'k}^T CY V~*_{lk'}`, so the noise enters only
        through the reduced `CY` that `_cy_reduced` already refuses to
        hand over when it is bias-dependent.

        ⚠ `CY/2`, NOT `CY`.  The library's `CY` is one-sided; eq (22)
        integrates `B B^T` as a two-sided intensity.  Consistent with the
        `kT/C`-calibrated Monte Carlo injection `Var(i) = CY/(2h)` in this
        file's record, and confirmed here by three routes agreeing.

        ⚠⚠ GATED THREE WAYS, because a modal sum transcribed from an image
        of an equation is exactly the object this record distrusts.  On
        van der Pol under gear: (i) this sum against `R_yy(0)` evaluated
        from its DEFINITION as a 1-D Lyapunov integral along the orbital
        mode, no Fourier machinery — agree to 3.5e-4; (ii) both against
        the CYCLE-MEAN transverse part of `oscillator_covariance`'s
        samples (`P(t_j) - (t_j/T)·growth_samples[j]`), which shares no
        machinery with either — magnitude to < 1e-3.  That third route is
        what found the state-block scale defect in `floquet_modes`.

        ⚠ THE REFERENCE IS THE CYCLE MEAN, NOT `K_orb(0)`.  Lemma 3.5's
        `R∞_yy` depends on `τ` only — the stationary part.  At `t = 0`
        van der Pol's amplitude direction is pure-v while this is
        isotropic, which is a rotating radial direction averaged over a
        cycle, not a disagreement.

        ✅ THE 2-3 % SHAPE RESIDUAL WAS THE REFERENCE, NOT THIS SUM -- closed
        2026-09-04.  Subtracting only the SECULAR growth `(t/T) d u u^T`
        from the Lyapunov samples leaves the phase direction's BOUNDED
        within-period variance, which eq (22)'s `l >= 2` sum correctly
        excludes.  Demir's `y` is defined by the OBLIQUE projection
        `v_1^T y = 0`; project the samples with `Pi = I - u v^T/(v^T u)`
        and the three-way agreement is 5.9e-4 / 6.0e-4 / 3.1e-4 / 1.5e-4
        at Q = 4 / 8 / 16 / 32 (euler-plain, n = m), improving with
        refinement -- quadrature.  The old residual fell as 1/Q_lambda
        (5.4 / 2.4 / 1.2 / 0.6 %), which is orbital variance ~ Q against a
        constant phase-bounded part: the same fact, seen from the sweep.
        ⚠ A candidate recorded earlier -- the phase-orbital CORRELATION's
        tau = 0 value -- was DISPROVED from eq (18a) before it was built:
        at tau = 0 its brace is {1 - 1} = 0 identically, and (23) states
        R_yy(0) = sum C_lhj alone.  Named so nobody rebuilds it.
        """
        self._refuse_coloured(pss, 'orbital_correlation')
        H = self.ORBITAL_HARMONICS if H is None else int(H)
        modes = pss.floquet_modes(pss)
        m = self.cir.n - 1
        Tp = float(pss.period)
        w0 = 2.0 * np.pi / Tp
        CY2 = 0.5 * np.real(np.asarray(self._cy_reduced(pss, 0.0)))
        orb = [k for k, md in enumerate(modes)
               if abs(abs(md['lam']) - 1.0) > 1e-6]
        if not orb:
            raise ValueError(
                'PAC.orbital_correlation: no orbital mode -- every non-null '
                'multiplier sits on the unit circle.')

        def fcoef(P):
            X = np.asarray(P)[:, :-1]
            N = X.shape[1]
            return np.fft.fft(X, axis=1) / N, N

        U, V, N = {}, {}, None
        for k in orb:
            U[k], N = fcoef(modes[k]['p'])
            V[k], _ = fcoef(modes[k]['q'])
        H = min(H, N // 2 - 1)
        hs = np.arange(-H, H + 1)
        idx = lambda k: k % N

        C = {}
        R = np.zeros((m, m), dtype=complex)
        for l in orb:
            mul = modes[l]['mu']
            for lp in orb:
                mulp = modes[lp]['mu']
                for j in hs:
                    Ulj = U[l][:, idx(j)]
                    ## the Lambda products for every (h, j') at once
                    Vl_hj = V[l][:, idx(hs - j)]              # m x nh  (h - j)
                    for jp in hs:
                        res = 1.0 / (1j * (j - jp) * w0 - mulp - np.conj(mul))
                        outer = res * np.outer(U[lp][:, idx(jp)], np.conj(Ulj))
                        Vlp_hjp = V[lp][:, idx(hs - jp)]      # m x nh  (h - j')
                        sc = np.einsum('ih,ik,kh->h', Vlp_hjp, CY2, np.conj(Vl_hj))
                        for hi, h in enumerate(hs):
                            term = sc[hi] * outer
                            key = (l, int(h), int(j))
                            C[key] = C.get(key, 0.0) + term
                            R += term
        return np.real(R), C

    def orbital_spectrum(self, pss, offsets, output, harmonic=1, H=None):
        """`S_yy` — the ORBITAL (amplitude) noise spectrum. A9 step 4.

        Returns `S` at `harmonic*f0 + offsets`, in the same V^2/Hz scale as
        `oscillator_spectrum`'s `S_v`, so **the two are summed** — which is
        what Traversa & Bonani (TCAS-I 2011) say to do:

            x(t) = x_s(t + a(t)) + y(t)      a = phase, y = orbital

        with the phase--orbital CROSS term dropped.  ⚠ That is a documented
        approximation with a KNOWN SIGN, not an oversight: the paper reports
        the correlation spectrum negligible on two circuits, and that when
        present it *decreases* the total.  **Dropping it therefore OVER-states
        noise** — conservative for a design margin, wrong in a known
        direction.  It is identically zero with no AM-to-PM coupling.

        **Lemma 3.5**: the orbital spectrum is a sum of Lorentzians centred at
        `j*w0 + Im(mu_l)` with half-width `|Re(mu_l)| + (1/2) h^2 w0^2 c`,
        weighted by the `C_lhj` of eq (22).  Every input already exists:
        `orbital_correlation` returns `C_lhj` (gated three ways), the
        exponents come from `floquet_modes`, and `c` from
        `diffusion_constant`.

        ⚠ UNIT CONVERSION, DONE ONCE HERE.  Lemma 3.5's widths are ANGULAR.
        `(1/2) h^2 w0^2 c` rad/s is `pi h^2 f0^2 c` Hz -- exactly the
        half-width `lorentzian` already uses for the phase line -- and
        `|Re(mu_l)|` rad/s is `|Re(mu_l)|/(2 pi)` Hz.  The two half-widths
        ADD, so an orbital mode's line is the phase line broadened by the
        mode's own relaxation rate.

        ⚠⚠ AND THAT IS WHY IT MATTERS AT LARGE OFFSET, WHICH IS THE WHOLE
        POINT OF THE ITEM.  The phase line's width is `pi h^2 f0^2 c`, which
        for a good oscillator is tiny, so its skirt has fallen as `1/f^2` long
        before the orbital line -- width `|Re(mu_2)|/(2 pi)`, i.e. the
        AMPLITUDE RELAXATION RATE -- has even started to roll off.  The
        crossover therefore sits near

            f_amp = -ln(lam2) f0 / (2 pi) = f0 / (2 pi Q)

        the same pole `oscillator_spectrum` warns about from the other side.
        ⚠ Those two arrived independently -- one from a commercial
        simulator's excess over our phase-only answer, one from this paper's
        modal sum -- and they must land in the same place.  That is the gate
        (`test_the_orbital_spectrum_crosses_the_phase_spectrum_near_f_amp`),
        and it is the check that can actually fail.

        ⚠ `output` follows `oscillator_spectrum`: an integer indexes the
        REDUCED state (the reference row already removed), an array is a
        weight vector over it.

        ⚠ STATIONARY WHITE SOURCES ONLY -- inherited from
        `orbital_correlation`, which needs `CY` constant for eq (22)'s
        products to collapse.
        """
        self._warn_if_orbit_is_asymmetric(pss)
        R, C = self.orbital_correlation(pss, H=H)
        modes = pss.floquet_modes(pss)
        c = float(self.diffusion_constant(pss))
        f0 = 1.0 / float(pss.period)
        m = pss.cir.n - 1

        d = np.asarray(output)
        if d.ndim == 0:
            row = np.zeros(m, dtype=float)
            row[int(d)] = 1.0
        else:
            row = np.asarray(d, dtype=float).ravel()[:m]

        f = float(harmonic) * f0 + np.atleast_1d(
            np.asarray(offsets, dtype=float))
        S = np.zeros_like(f, dtype=float)
        for (l, h, j), Clhj in C.items():
            ## The weight is the output's own share of this term.  It is real
            ## for the total (`R` is real symmetric); an individual `(l,h,j)`
            ## can carry a small imaginary part that cancels against its
            ## conjugate partner, so take the real part per term rather than
            ## asserting each is real.
            w = float(np.real(row @ Clhj @ row))
            if w == 0.0:
                continue
            mul = modes[l]['mu']
            ## Hz, both terms -- see the unit note above.
            gam = abs(float(np.real(mul))) / (2.0 * np.pi) \
                + np.pi * float(h) ** 2 * f0 ** 2 * c
            fc = float(j) * f0 + float(np.imag(mul)) / (2.0 * np.pi)
            if gam <= 0.0:
                continue
            ## Normalised Lorentzian: integrates to 1 over all `f`, so the
            ## total power is `sum(w) = row^T R row` by construction.
            S = S + w * (gam / np.pi) / ((f - fc) ** 2 + gam ** 2)
        return S

    def diffusion_constant(self, pss):
        """`c` — the phase diffusion constant, in seconds.

        `c = (1/T) ∫ v₁ᵀ(t) B(t) Bᵀ(t) v₁(t) dt` with `B Bᵀ = CY`, so this
        is the time-average of a QUADRATIC form in the PPV.  It is the one
        scalar the whole free-running phase-noise spectrum is built from,
        and it reads, for a designer, as JITTER PER SECOND.

        ⚠ QUADRATIC FOR WHITE SOURCES, LINEAR FOR COLOURED ONES, and the
        two are different functionals of the same vector: a coloured
        source contributes `V_0m = (1/T) ∫ v₁ᵀ B_cm dt`, with no square.
        Using this one for a coloured source returns a plausible non-zero
        number from the same PPV.  Only stationary white sources are
        supported here, which `_cy_reduced` enforces.

        ⚠ ITS SCALE WAS WRONG BY 2x AND IS NOW FIXED -- kept because the
        way it survived is the instructive part.  `diffusion_constant` used
        the full `CY` while `covariance` used `CY/2`: two functions in this
        class disagreeing about whether `CY` is one- or two-sided.  It was
        validated against a Monte Carlo injecting `Var(i) = CY/h` per step
        and agreed to 0.9965 -- because that Monte Carlo carried the SAME
        hot convention.  A measurement built on the assumption under test
        cannot test it.

        SETTLED AGAINST `kT/C`, which is external to both: an injection of
        `Var(i) = CY/h` reproduces 1.92x `kT/C` over ten independent runs
        (1.75-2.04).  With `CY/2` throughout, `diffusion_constant` gives
        7.9516e-08 against a correctly scaled Monte Carlo at 7.7083e-08 --
        ratio 1.0316, inside that measurement's 4.1% uncertainty.

        ⚠ AND A SECOND DISCREPANCY WAS NOT A CODE DEFECT AT ALL.  Two Monte
        Carlo routes disagreed by 2.31x, which looked like a third error.
        It was in the DIAGNOSTIC: `ppv()` normalises on the FIRST BLOCK
        (`v[:m] . xdot = 1`), which is right for a perturbation entering
        the first block -- an injected current, and what every shipped path
        does -- but wrong for contracting against a full PAIR deviation,
        where the factor is `1/(v . u_pair) = 1.508`.  Correcting it turned
        a 2.13 variance ratio into 1.07.  The sign difference alongside it
        is a convention, not an error: a later zero crossing means DELAYED,
        while projecting onto the tangent makes positive mean ADVANCED.
        """
        self._check_circuit(pss)
        self._refuse_coloured(pss, 'diffusion_constant')
        self._refuse_driven(pss, 'diffusion_constant')
        return self._white_diffusion_at(pss, 2.0 * np.pi / float(pss.period))

    def _refuse_driven(self, pss, what):
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.%s: phase diffusion is a property of a '
                "FREE-RUNNING oscillator. A driven circuit's phase is its "
                "source's, and its noise is pnoise's problem, not this one."
                % what)

    def _white_diffusion_at(self, pss, w):
        """`(1/T) integral v_1^T (CY(w)/2) v_1 dt` with `CY` FROZEN at `w`.

        The white functional at one frequency, with no refusal: it is `c`
        when the source is white, and for a coloured source it is the
        value `diffusion_constant` used to return silently.  `phase_psd`
        reads it at the carrier for the Lorentzian CORNER, which is a
        white-noise construct whatever the source's colour; the spectrum
        itself comes from `coloured_diffusion_resolved`.
        """
        v, info = pss.ppv()
        ## `lambda_2` is computed here anyway; `oscillator_spectrum` needs it to
        ## report its own validity limit and a second `ppv()` would be a full
        ## extra solve.  Recorded, not returned, so this method's signature is
        ## unchanged -- and read ONLY immediately after a call, which is how
        ## `oscillator_spectrum` uses it.
        self._last_second_multiplier = (
            info.get('second_multiplier'),
            info.get('second_multiplier_certified'))
        m = pss.cir.n - 1
        ## ⚠ `samples_eq`, NOT `samples`.  `CY` is an EQUATION-ROW
        ## covariance and `samples` is `C^T v_1`; contracting it here made
        ## `c` wrong by `C^2` on the differential rows and exactly zero on
        ## the algebraic ones.  See `_equation_row_ppv`.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        h = np.diff(tms)
        T = float(pss.period)
        cy = self._cy_reduced(pss, float(w))
        ## ⚠ `cy/2`, THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION `covariance`
        ## USES.  `CY` is a one-sided density (a resistor's `4kT/R`), and
        ## these two functions disagreed about it until a Monte Carlo was
        ## run against `kT/C`: an injection of `Var(i) = CY/h` per step
        ## reproduces `1.92x kT/C` over ten independent runs (1.75-2.04),
        ## so that convention carries TWICE the physical noise power.
        ## `covariance` was already right; this was not, and its agreement
        ## with a Monte Carlo built on the SAME hot convention is exactly
        ## why the error survived.
        quad = np.einsum('ij,jk,ik->i', S, 0.5 * np.real(cy), S)
        return float((quad * h).sum() / T)

    def colour_projection(self, pss):
        """`<v_1>` — the PPV's TIME AVERAGE, which is a different functional.

        Returns `(vbar, info)`.  `vbar` is `(1/T) integral v_1(t) dt` over
        the orbit; `info` carries the per-entry `rms` and the ratio
        `|mean|/rms`, which is the number that says whether a coloured
        source at that node can upconvert at all.

        ⚠ COLOURED SOURCES CONTRACT THE SQUARE OF THE MEAN; WHITE ONES
        CONTRACT THE MEAN OF THE SQUARE.  `diffusion_constant` computes
        `(1/T) integral v^T (CY/2) v dt`.  A coloured source's low-frequency
        power cannot be modulated away, so what survives is
        `V_0m = (1/T) integral v_1^T B_cm dt` -- LINEAR, no square -- and the
        contraction is `vbar^T (CY/2) vbar`.  Same vector, same matrix,
        the mean and the square exchanged.

        ⚠ AND USING THE QUADRATIC ONE FOR A COLOURED SOURCE RETURNS A
        PLAUSIBLE NUMBER, NOT AN ERROR.  It is never zero where the white
        answer is not, so nothing downstream would look wrong.  The
        measured separation is 22 ORDERS on van der Pol -- `c = 7.95e-08`
        against `Gamma = 1.9e-29` -- so the two functionals are not close
        approximations of each other and cannot be substituted.

        ⚠ TWO INDEPENDENT MECHANISMS FORCE `vbar` TO ZERO, AND ONLY ONE OF
        THEM IS THE ONE DESIGNERS KNOW.  Measured on an LC oscillator,
        sweeping an even term `a (u^2 - 2)` in the nonlinearity and a
        series tank resistance `Rs`:

            a      Rs      Gamma/c
            0.00   0.00    2.4e-22
            0.00   0.20    9.7e-23
            0.25   0.00    4.9e-23
            0.25   0.05    2.1e-04
            0.25   0.20    4.1e-03

        NEITHER ASYMMETRY ALONE NOR LOSS ALONE UPCONVERTS.  `c` is
        7.9e-08 to 1.2e-07 in every row, so the quadratic functional
        cannot produce that pattern.

        * A SYMMETRIC waveform gives `vbar = 0` -- Hajimiri & Lee, and the
          reason symmetry is the first thing a VCO designer reaches for.
        * A LOSSLESS LC TANK gives `vbar[0] = 0` STRUCTURALLY, whatever the
          waveform does.  `v` behaves as `C^T v_1` and `dv/dt = G^T v_1`,
          whose inductor row is exactly `v[0]`; periodicity of `v[1]` then
          forces `integral v[0] dt = 0`.  ⚠ THIS IS A PROPERTY OF THE
          TOPOLOGY, NOT OF THE ORBIT, and it is why van der Pol reports
          zero at every asymmetry -- it makes van der Pol useless as a
          POSITIVE fixture and perfect as a negative one.

        ⚠ `Gamma <= c` ALWAYS, at the same `CY`, by Cauchy-Schwarz on the
        weighted mean -- with equality only if `v` is constant over the
        orbit.  Both use the same quadrature here so the bound holds
        exactly at the discrete level, which makes it an assertion rather
        than an expectation.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.colour_projection: the PPV time-average is the kernel '
                "of a FREE-RUNNING oscillator's coloured-noise upconversion. "
                "A driven circuit's phase is its source's.")
        _v, info = pss.ppv()
        m = pss.cir.n - 1
        ## ⚠ the EQUATION-ROW adjoint, for the same reason
        ## `diffusion_constant` uses it: a coloured source is an
        ## equation-row input too.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        h = np.diff(tms)
        T = float(pss.period)
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES, deliberately:
        ## it is what makes `Gamma <= c` exact rather than approximate.
        vbar = (S * h[:, None]).sum(0) / T
        rms = np.sqrt((S ** 2 * h[:, None]).sum(0) / T)
        with np.errstate(divide='ignore', invalid='ignore'):
            sym = np.where(rms > 0, np.abs(vbar) / rms, 0.0)
        return vbar, {'rms': rms, 'symmetry': sym,
                      'samples': S, 'times': tms}

    def coloured_diffusion(self, pss, freqs):
        """`Gamma(f) = vbar^T (CY(2 pi f)/2) vbar` — the coloured analogue of `c`.

        Returns an array over `freqs`.  `CY` is evaluated at each offset,
        so a source whose density varies with frequency -- which is what
        "coloured" means -- is folded in exactly as a white one is.

        ⚠ NO FILTER, NO EXTRA STATE, NO SDE.  Demir 1996 synthesises 1/f
        from white sources through a Lorentzian network at "one state
        variable per decade", because Ito theory admits only white driving
        noise.  That is an artefact of the SDE formulation.  This path
        never forms an SDE, so a coloured source is just a different
        `S(f)` -- a SLOPE, NOT A STATE.  A commercial RF simulator confirms by omission:
        no filter and no augmentation in its treatment of flicker.

        ⚠ THE `CY/2` IS THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION THE
        REST OF THIS CLASS USES, and it is shared rather than repeated so
        the pair cannot drift the way `diffusion_constant` and
        `covariance` once did over exactly that factor.
        """
        vbar, _ = self.colour_projection(pss)
        out = []
        for f in np.atleast_1d(np.asarray(freqs, dtype=float)):
            cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * float(f)))
            out.append(float(vbar @ (0.5 * cy) @ vbar))
        return np.asarray(out)

    def coloured_diffusion_resolved(self, pss, freqs, harmonics=None):
        """`c(f) = sum_l V_l^H (CY(2 pi |f - l f_0|)/2) V_l` — the fold PER HARMONIC.

        `V_l` are the Fourier coefficients of the equation-row PPV `v_1(t)`
        (the rows `diffusion_constant` contracts), so a source's density is
        read at the SOURCE-SIDE frequency `f - l f_0` for each harmonic it
        folds through -- which is what `pnoise` has done from the start and
        what a coloured source requires.  Returns an array over `freqs`.

        ⚠ THIS IS THE OBJECT `c + Gamma(f)` STOOD IN FOR, and the stand-in
        is wrong in two ways that the fixture could not show: `c` reads
        `CY` at ONE frequency (`2 pi / T`) as if it held at every harmonic,
        and `Gamma` is exactly the `l = 0` term of this sum, so `c + Gamma`
        counts `l = 0` twice.  Neither was visible on van der Pol, whose
        PPV at the tank node averages to zero (`|V_0|/|V_1| = 5e-13`: the
        inductor shorts the node at DC, so no core can bias it) -- the
        fixture shared the claim's assumption, failure shape 0b.

        EXACT FOR WHITE, BY PARSEVAL: with `CY` constant the sum is
        `(1/T) integral v_1^T (CY/2) v_1 dt = c`, and the discrete version
        with the grid's step weights reproduces `diffusion_constant` to
        round-off -- that equality pins the transform's normalisation, and
        it is asserted.  For a DC-centred colour (Lorentzian, flicker) and
        `f << f_0` the `l != 0` terms read `CY(l f_0)` to `O(f/f_0)`, so
        the sum differs from `c + Gamma` only where `V_0` is not small.

        `harmonics` caps `|l|`; by default every harmonic carrying more
        than 1e-14 of the PPV's energy is kept, which is all of them that
        can move the sum at double precision.
        """
        self._check_circuit(pss)
        self._refuse_driven(pss, 'coloured_diffusion_resolved')
        m = pss.cir.n - 1
        v0, info = pss.ppv()
        S = np.asarray(info['samples_eq'], dtype=float)[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        n = S.shape[0]
        T = float(pss.period)
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES: one sample per
        ## step, weighted by that step, so that Parseval closes exactly.
        t = tms[1:1 + n]
        h = np.diff(np.concatenate(([tms[0]], t)))
        w0 = 2.0 * np.pi / T
        L = n // 2 if harmonics is None else int(harmonics)
        ls = np.arange(-L, L + 1) if harmonics is not None else np.arange(-L, L)
        E = np.exp(-1j * np.outer(ls, w0 * t)) * h[None, :]          # (nl, n)
        V = (E @ S) / T                                               # (nl, m)
        energy = np.sum(np.abs(V) ** 2, axis=1)
        keep = energy > 1e-14 * energy.sum()
        ls, V = ls[keep], V[keep]
        out = []
        for f in np.atleast_1d(np.asarray(freqs, dtype=float)):
            tot = 0.0
            for l, vl in zip(ls, V):
                cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * abs(float(f) - l / T)))
                tot += float(np.real(np.conj(vl) @ (0.5 * cy) @ vl))
            out.append(tot)
        return np.asarray(out)

    def phase_psd(self, pss, offsets, harmonic=1):
        """`S_phi(f)` in rad^2/Hz at `offsets` from harmonic `i` — white AND coloured.

            S_phi,i(f) = i^2 f_0^2 c(f) / f^2,   c(f) = sum_l V_l^H (CY(f - l f_0)/2) V_l

        `c(f)` is `coloured_diffusion_resolved`: the phase diffusion with
        each harmonic's colour read at its own source-side frequency.  For
        a white source it is `c` exactly; the earlier `c + Gamma(f)` form
        counted the `l = 0` term twice and is retired.

        ⚠ THE CONVENTION IS PINNED BY `oscillator_spectrum`, NOT ARGUED.
        `lorentzian`'s far skirt is `i^2 f_0^2 c / f^2` exactly, and that
        object was gated by power conservation to 1.000000.  So this
        expression is the same quantity its tail already reports, with the
        coloured term added -- no second convention is introduced, which
        is the only reason a `S_phi` is shipped here at all after a
        one-sided/two-sided error cost this class a factor of two.

        ⚠ THE TWO TERMS ADD BECAUSE THE SOURCES ARE INDEPENDENT, and with
        `CY ~ 1/f` the coloured term gives `S_phi ~ 1/f^3` -- Kundert's
        "S_u(f) is generally pink ... then S_phi(f) would be proportional
        to 1/f^3 at low frequencies".

        ⚠ AND THIS IS THE LINEARISED PHASE MODEL, WHICH IS EXACT ENOUGH
        ONLY BECAUSE THE SOURCES ARE STATIONARY.  Vanassche, Gielen &
        Sansen (ICCAD 2002) locate the split between the exact phase
        equation `theta' = eps Gamma(t + theta) n(t)` and the approximate
        `theta' = eps Gamma(t) n(t)`: for a STATIONARY source the two
        "will, up to 0-th order in eps, predict the same output phase
        noise", and they diverge otherwise.  Their operational form is
        better than "non-stationary" -- "at first, near t = 0, the
        predicted phases are the same. However, when THETA BECOMES TOO
        LARGE [they diverge]".  A stationary source makes `theta` DIFFUSE;
        a driven one makes it grow SECULARLY, which is what carries it out
        of range.  ⚠ So this is sound for free-running noise and must NOT
        be reused for injection locking, a PLL in lock, or coupled
        oscillators -- there the shift has to stay inside the argument.

        ⚠ REFUSED BELOW THE LORENTZIAN CORNER, and this is a validity
        boundary rather than a conditioning one.  There the excess phase is
        a Wiener process whose spectrum is singular at the origin; the
        finite value the real lineshape attains comes from the NONLINEAR
        phase-to-voltage map, which `oscillator_spectrum` carries and this
        does not.  Reporting `S_phi` near the carrier is the mistake this
        object invites, so it raises instead.
        """
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        i = int(harmonic)
        if i < 1:
            raise ValueError('PAC.phase_psd: harmonic must be >= 1.')
        offs = np.atleast_1d(np.asarray(offsets, dtype=float))
        if np.any(offs <= 0.0):
            raise ValueError(
                'PAC.phase_psd: offsets must be positive; S_phi diverges '
                'at zero offset and that divergence is physical.')
        cres = self.coloured_diffusion_resolved(pss, offs)
        ## ⚠ THE CORNER IS THE WHITE LORENTZIAN'S, read at the carrier as it
        ## always was.  For a coloured source `f_h = pi i^2 f0^2 c` is not a
        ## lineshape parameter at all -- there is no Lorentzian -- and
        ## taking the folded value nearest the carrier instead put a 1/f
        ## source's corner ABOVE the offsets, in front of the power bound
        ## below, which is the floor that actually binds for colour.
        c = self._white_diffusion_at(pss, 2.0 * np.pi * f0)
        ## The i-th harmonic's Lorentzian half-width.  `S_i(f) =
        ## i^2 f0^2 c / (pi^2 i^4 f0^4 c^2 + f^2)` is a Lorentzian in `f`
        ## whose denominator is `f_h^2 + f^2`, so `f_h = pi i^2 f0^2 c`.
        corner = np.pi * (i ** 2) * (f0 ** 2) * c
        if offs.min() <= corner:
            raise ValueError(
                'PAC.phase_psd: offset %.6g Hz is at or below the '
                'Lorentzian corner %.6g Hz for harmonic %d, where S_phi is '
                'not the right object -- the excess phase is a Wiener '
                'process and its spectrum is singular at the origin. The '
                'finite value the LINESHAPE attains there comes from the '
                'nonlinear phase-to-voltage map: use oscillator_spectrum().'
                % (float(offs.min()), corner, i))
        sphi = (i ** 2) * (f0 ** 2) * cres / offs ** 2

        ## ⚠ POWER CONSERVATION AS A SECOND, INDEPENDENT FLOOR -- and for a
        ## COLOURED source it is the binding one, by orders.  The
        ## normalised lineshape integrates to 1, and the integral over one
        ## box of width `df` on each side is a lower bound on it, so
        ##
        ##     2 df S_phi(df) <= 1
        ##
        ## is NECESSARY for the linearised skirt to be consistent with
        ## unit power.  Vanassche, Gielen & Sansen (2003) derive the same
        ## statement for a 1/f input and reduce it to
        ## `df_c >= eps f0 sqrt(2 f_1f)`; the form here needs no
        ## assumption about the source's colour, and REPRODUCES their
        ## worked example exactly -- 100.000 Hz against their ">= 100 Hz"
        ## at `eps^2 = 1e-19`, `f0 = 1 GHz`, `f_1f = 50 kHz`.
        ##
        ## ⚠ THE LORENTZIAN CORNER ABOVE DOES NOT CATCH THIS.  It is built
        ## from `c` alone, so it knows nothing about a `Gamma(f)` that
        ## grows as the offset falls.  MEASURED on this class's own
        ## flicker fixture: the power bound bites at 2.5e-06 Hz while the
        ## Lorentzian corner sits at 8.2e-09 Hz -- 306x too permissive,
        ## and the swept spectrum was carrying 3.10x unit power at the
        ## bottom of the range before this check existed.
        ##
        ## ⚠ IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN.
        ## Passing it is not a guarantee: on Vanassche's own example the
        ## observed flattening sits at ~300 Hz, 3x the bound.  So this
        ## refuses what is definitely invalid and admits a band that is
        ## already suspect -- deliberately, because refusing at 3x would
        ## be fitting a threshold to one example.
        ## ⚠ AND THE DERIVATION HAS A PRECONDITION THE BOUND DOES NOT
        ## STATE, so it is checked rather than assumed.  The box argument
        ## is `2 df S(df) <= integral_{-df}^{+df} S <= 1`, and the FIRST
        ## inequality needs `S(f) >= S(df)` for every `|f| <= df` -- the
        ## spectrum must not dip below its edge value anywhere further in.
        ## True of a monotone skirt; TRUE of the flattened near-carrier
        ## shape; true even with a spur, which ADDS power inside rather
        ## than creating a dip.
        ##
        ## ⚠ FALSE FOR A LOCKED PLL, whose phase-noise transfer function
        ## is HIGH-PASS: the spectrum is SUPPRESSED at DC and rises to the
        ## free-running level beyond the loop bandwidth, so it dips below
        ## its edge value everywhere inside.  The bound is not thereby
        ## shown to be violated there -- total power is still 1 -- it is
        ## NO LONGER DERIVED, and a floor that is not derived cannot be
        ## used as one.  Unreachable today because this method refuses a
        ## driven circuit, and squarely in the way of the driven-oscillator
        ## work, which is why it is a check and not a comment.
        probe = np.unique(np.concatenate((
            offs, np.logspace(np.log10(offs.min() / 1e3),
                              np.log10(offs.max()), 32))))
        sprobe = ((i ** 2) * (f0 ** 2)
                  * self.coloured_diffusion_resolved(pss, probe) / probe ** 2)
        if np.any(np.diff(sprobe) > 1e-12 * np.abs(sprobe[:-1])):
            k = int(np.argmax(np.diff(sprobe) > 0)) + 1
            raise ValueError(
                'PAC.phase_psd: the spectrum RISES with offset near '
                '%.6g Hz, so it dips below its edge value further in and '
                'the power bound below is no longer derived -- its box '
                'argument needs S(f) >= S(df) for every |f| <= df. That '
                'happens for a high-pass-shaped spectrum such as a locked '
                'loop, and for a source whose density grows faster than '
                'f^2. The bound may still hold; it is not established '
                'here, so it is refused rather than applied.'
                % float(probe[k]))

        power = 2.0 * offs * sphi
        bad = power >= 1.0
        if np.any(bad):
            k = int(np.argmax(bad))
            raise ValueError(
                'PAC.phase_psd: at offset %.6g Hz the linearised skirt '
                'already carries %.3f times the TOTAL power of the '
                'carrier (2 f S_phi >= 1), so it has broken down there -- '
                'a normalised spectrum integrates to 1. This bound is '
                'independent of the Lorentzian corner (%.6g Hz here) and '
                'for a coloured source it binds far earlier, because '
                'Gamma(f) grows as the offset falls. Sweep above it, or '
                'use oscillator_spectrum() for the lineshape. Note the '
                'TRUE breakdown is higher still: this is a lower bound.'
                % (float(offs[k]), float(power[k]), corner))
        return sphi

    @staticmethod
    def lorentzian(offsets, c, f0, harmonic=1):
        """The `i`-th harmonic's normalised lineshape at `offsets` from it.

            S_i(f) = i² f₀² c / (π² i⁴ f₀⁴ c² + f²)

        ⚠ EXACT FOR WHITE SOURCES, not a limiting form.  With coloured
        sources the transform "does not have a simple closed form" and only
        two-regime approximations exist — which is one more reason this
        module supports white sources only.

        ⚠ AND ITS TOTAL POWER IS EXACTLY 1.  `∫ a/(b²+f²) df = aπ/b`, and
        here `a = i² f₀² c`, `b = π i⁴ f₀⁴ c² ^ ½`… concretely `b = π i²
        f₀² c`, so the integral is exactly one.  **The carrier's power is
        redistributed, never created or destroyed** — which is the
        invariant that separates this from LTV small-signal treatments,
        which "erroneously predict infinite noise power [at the carrier] as
        well as infinite total integrated power".  It is asserted in the
        suite.

        The half-width is `π i² f₀² c` and the peak `1/(π² i² f₀² c)`, so a
        higher harmonic has a skirt scaling as `i²` and a corner as `i⁴` —
        `20 log₁₀(i)` dB noisier far out.
        """
        i = int(harmonic)
        if i == 0:
            return np.zeros_like(np.asarray(offsets, dtype=float))
        f = np.asarray(offsets, dtype=float)
        a = (i * i) * f0 * f0 * c
        b = np.pi * (i * i) * f0 * f0 * c
        return a / (b * b + f * f)

    def oscillator_spectrum(self, pss, offsets, output, harmonic=1):
        """Free-running output spectrum at `offsets` from harmonic `harmonic`.

        ⚠⚠ THIS DOES NOT GO THROUGH `pnoise`'s SIDEBAND FOLD, AND IT CANNOT.
        The fold is a FREQUENCY-CONVERSION computation, and for a driven circuit
        -- a mixer, the diode-mixer fold case -- that is complete.  For an
        AUTONOMOUS oscillator it is structurally incomplete, and what it omits is
        exactly the near-carrier phase-noise skirt this method returns.  Rizzoli,
        Mastri & Masotti (IEEE MTT 42-807, 1994) state it directly: frequency
        conversion alone is insufficient for autonomous circuits, because the
        noise-induced FREQUENCY MODULATION OF THE CARRIER at low offsets is not a
        frequency-conversion effect (verified at the source 2026-09-08: p. 807,
        Introduction, verbatim "frequency-conversion techniques alone are not
        sufficient to solve the noise analysis problem for general autonomous
        circuits (oscillators). An important further aspect that must be taken
        into account is the noise-induced frequency modulation of the carrier
        taking place at low frequency offsets, which is not a
        frequency-conversion effect").  Their Section III (p. 810) NAMES the
        two stacks: CONVERSION noise, power exchanged among the sidebands of
        the unperturbed steady state, "invariably raises as 1/f for f -> 0,
        which is not consistent with the measured behavior"; MODULATION noise,
        "a jitter of the oscillatory steady state", proportional to noise power
        over f^2 so the PSD "raises as 1/f^3 for f -> 0 in agreement with the
        measured performance"; and the two DECOUPLE exactly at the steady state
        (M_BH = M_HB = 0).  They also say the two are "usually nearly equal" in
        an INTERMEDIATE offset band, "so that (20) and (21) are
        interchangeable" -- a cross-stack agreement test this tree does not yet
        have (recorded in the roadmap, not built).  ⚠ Their construction is
        harmonic balance; what transfers is the classification, the two slopes
        and the interchangeability, none of which need HB.  Diagnostic value:
        a FLAT PSD near the carrier is neither slope -- it is the Phi(T) - I
        singularity, not the conversion model being the wrong physics.

        So the two stacks -- the Floquet/PPV one (`ppv`, `diffusion_constant`,
        this method) and the sideband fold (`pnoise`) -- ARE NOT TWO
        IMPLEMENTATIONS OF ONE QUANTITY, and unifying them is not a
        simplification waiting to be made.  ⚠ THE HAZARD IS THAT THE WRONG ONE
        STILL RETURNS A NUMBER: deriving oscillator phase noise from the fold
        alone would produce a spectrum -- the conversion terms are real and
        non-zero -- just one missing the dominant contribution near the carrier.
        A plausible wrong answer, which is the failure shape this whole area
        keeps generating.  That is the completeness argument for the split; the
        efficiency argument (Floquet is cheaper) is the weaker one and was for a
        long time the only one written down.

        Returns `(S_v, L_dBc)`.  ⚠ `S_v` is the Lorentzian lineshape scaled by
        `|X_1|^2 = A^2/4`, the carrier PHASOR's square -- which is HALF the
        carrier power `A^2/2` a one-sided PSD carries, so `S_v` is exactly
        0.5000x a one-sided PSD of the output voltage (measured against
        a reference simulator at every offset over four decades, 2026-09-05).  `L_dBc`
        is unaffected, `|X_1|^2` dividing out of the ratio; the absolute
        V^2/Hz matters to anyone integrating `S_v` to a power, and the
        scale is kept rather than doubled because it is a return value
        that callers may already divide by `|X_1|^2` themselves.  `S_v`
        was documented as the one-sided PSD of the output
        voltage; `L_dBc` is that normalised to the harmonic's own power,
        in dBc/Hz.

        ⚠ NO SWEEP AND NO PER-FREQUENCY SOLVE.  Once the PSS waveform's
        Fourier coefficients and the scalar `c` are known, "we have an
        analytical expression that gives us the spectrum at any frequency.
        The computation of the spectrum is not performed separately for
        every frequency of interest."  Which also means it never meets the
        near-carrier singularity that a swept small-signal computation
        would, and never meets the 1/f sweep-grid trap — there is no sweep
        to place a point on.

        ⚠⚠ SCOPE: A SOURCE BEHIND A SLOW NODE (A2, resolved 2026-09-08).
        The Lorentzian uses the DC PPV, so for a noise source that reaches
        the core through a slow path (RC leg, tau >> T) it holds only
        BELOW the source's corner `T/(2 pi tau)`; above it the true skirt
        is this one scaled by the PPV-harmonic-weighted filter
        `sum_k |G_k|^2 F_k(f) / sum_k |G_k|^2 F_k(0)` (G_k the PPV entry's
        Fourier coefficients at the source node, F_k the path's transfer at
        k f0 + f), which is 1/1000 at 0.1 f0 on a one-RC-leg fixture with
        an asymmetric core AND tank loss (both needed for G_0 != 0; an
        ideal tank inductor shorts DC).  `c` is still right (the filter
        removes only high-frequency content); `pnoise` computes the true
        value at any offset; this method does not, and a Monte Carlo of
        `c` cannot see it.  Measured against `pnoise` to four digits
        through the corner (test ..._behind_a_slow_node_...).

        ⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
        small-signal analysis cannot produce `L(f)` there however well
        conditioned it is: the excess phase is a Wiener process, its
        spectrum has a singularity at the origin and no physical meaning,
        and the finite value `L` attains comes from the NONLINEAR
        phase-to-voltage map — which is what this closed form carries.
        Reporting `S_phi` near the carrier instead is the mistake that
        object invites.
        """
        c = self.diffusion_constant(pss)
        f0 = 1.0 / float(pss.period)
        self._warn_above_amplitude_pole(offsets, f0)
        X = self.carrier_phasor(pss, output, harmonic)
        Sv = abs(X) ** 2 * self.lorentzian(offsets, c, f0, harmonic)
        with np.errstate(divide='ignore'):
            L = 10.0 * np.log10(np.maximum(Sv / max(abs(X) ** 2, 1e-300),
                                           1e-300))
        return Sv, L

    def _warn_above_amplitude_pole(self, offsets, f0):
        """⚠ THE PHASE-ONLY SPECTRUM IS A LOWER BOUND ABOVE `f_amp`.

        `oscillator_spectrum` returns the PHASE contribution only.  A real
        oscillator also carries AMPLITUDE noise, which is suppressed near the
        carrier because the limit cycle restores the amplitude -- but only at
        the amplitude-relaxation rate.  Above the pole where that restoring
        action runs out, amplitude noise stops decaying within a period and
        adds to the total, so this method UNDER-reports.  Relayed measurement
        against a commercial simulator's total noise, as excess over the
        phase-only prediction:

            offset     lam2 = 0.90            lam2 = 0.99
                       (f_amp 26.7 kHz)       (f_amp 2.55 kHz)
            100 Hz     -0.00 dB               -0.01 dB
            1 kHz      -0.00 dB               -0.54 dB
            10 kHz     -0.50 dB               -2.90 dB
            100 kHz    -3.11 dB               -3.27 dB

        ⚠⚠ AND THE VALID REGION SHRINKS AS `1/Q`, which makes this section 0
        again rather than a detail.  With `f_amp = -ln(lam2)/(2 pi T)` and
        `Q = -1/ln(lam2)`,

            f_amp = f0 / (2 pi Q)

        -- verified both ways at 26671.9 / 2544.2 / 253.3 Hz for
        `lam2 = 0.90 / 0.99 / 0.999`.  So the better the oscillator, the
        narrower the band in which its phase-only spectrum is the whole
        answer; at `lam2 = 0.999` it has collapsed below ~253 Hz.

        ⚠ THIS IS THE OPPOSITE SIGN FROM THE ERROR `PSS.ppv` ALREADY WARNS
        ABOUT.  That one says the instantaneous phase equation misses slow
        nodes which FILTER device noise, so phase noise is OVER-estimated.
        This one is a second, independent mechanism in which the phase-only
        answer is UNDER-estimated.  Both are live and they are not the same
        effect.
        """
        lam2, certified = getattr(self, '_last_second_multiplier',
                                  (None, None))
        if lam2 is None:
            return
        lam2 = float(lam2)
        ## `lam2 <= 0` is a real or overdamped mode with no relaxation pole to
        ## speak of, and `lam2 >= 1` is not a decaying mode at all -- in both
        ## cases there is no `f_amp` and inventing one would be worse than
        ## silence.
        if not (0.0 < lam2 < 1.0):
            return
        ## `f_amp = -ln(lam2)/(2 pi T)` and `T = 1/f0`.
        f_amp = -np.log(lam2) * float(f0) / (2.0 * np.pi)
        off = np.atleast_1d(np.asarray(offsets, dtype=float))
        worst = float(np.max(np.abs(off))) if off.size else 0.0
        if worst < f_amp:
            return
        warnings.warn(
            'PAC.oscillator_spectrum: this is a PHASE-ONLY spectrum and %g Hz '
            'is above the amplitude-relaxation pole f_amp = %.4g Hz '
            '(lambda_2 = %.6f, f_amp = f0/(2*pi*Q)). Above f_amp the '
            'amplitude noise no longer decays within a period and adds to the '
            'total, so the value returned here is a LOWER BOUND: measured '
            'excess of a commercial simulator over the phase-only prediction '
            'is -0.54 dB at 1 kHz and -2.90 dB at 10 kHz for lambda_2 = 0.99. '
            '%sThe valid band scales as 1/Q, so it NARROWS as the oscillator '
            'improves.'
            % (worst, f_amp, lam2,
               ('' if certified is not False else
                'lambda_2 itself is NOT certified here (see '
                "info['second_multiplier_certified']), so f_amp is uncertain "
                'too. ')),
            RuntimeWarning, stacklevel=3)

    @staticmethod
    def am_pm_indices(a, b):
        """Split a sideband pair into AM and PM modulation indices.

        `a` and `b` are the upper and lower sideband amplitudes, each
        already divided by the carrier phasor.  Returns `(m_am, m_pm)`.

        THE WHOLE THING IS ONE CONJUGATE.  Write the complex envelope's
        deviation as `a e^{j w_m t} + b e^{-j w_m t}`.  The two sidebands
        COUNTER-ROTATE about the carrier phasor, so the sum traces an
        ellipse; the component ALONG the carrier is amplitude modulation
        and the component PERPENDICULAR to it is phase modulation.

          - pure AM keeps the envelope on the carrier's axis, which forces
            `d = conj(d)` for all `t`, i.e. `a = conj(b)`;
          - pure PM keeps it perpendicular, `d = -conj(d)`, i.e.
            `a = -conj(b)`.

        so `m_am = a + conj(b)` and `m_pm = a - conj(b)` -- each vanishing
        exactly when the other case holds.  No new solve: this is a change
        of basis on transfer functions `adjoint_sideband_row` already
        returns.

        ⚠ `conj(b)`, NOT `b`.  Using `a +- b` looks equally plausible and
        is wrong for any modulation whose sidebands are not real relative
        to the carrier -- it would report a rotating ellipse as pure AM.
        The conjugate is what makes the lower sideband counter-rotate.
        """
        a = np.asarray(a, dtype=complex)
        b = np.asarray(b, dtype=complex)
        return a + np.conj(b), a - np.conj(b)

    def _output_waveform_row(self, pss, output):
        """The steady-state waveform of `output`, as an index OR a direction.

        ⚠ THE REST OF THIS CLASS TAKES A DIRECTION VECTOR AND THIS PAIR
        TOOK AN INTEGER, which is not a style difference -- it meant
        `am_pm` and `carrier_phasor` could not express a DIFFERENTIAL
        output at all.  `pnoise`, `adjoint_transfer_row` and
        `adjoint_sideband_row` all accept `d`; these did `int(output)`.
        For an oscillator the output of interest is very often
        differential, and for the coordinate-invariance property an AM/PM
        split has to have (Kaertner 1990 section 3.2) a
        reference-independent observable is the whole point.

        An integer is still accepted, so callers that name a node keep
        working; an array is contracted against the full waveform with the
        reference row reinserted.
        """
        if getattr(pss, 'waveform', None) is None:
            raise RuntimeError(
                'PAC: the PSS has no stored waveform -- call solve() first.')
        _times, X = pss.waveform
        Xf = np.asarray(X, dtype=float)
        irn = pss.irefnode
        d = np.asarray(output)
        if d.ndim == 0:
            k = int(d)
            return Xf[k if k < irn else k + 1]
        row = np.zeros(Xf.shape[1], dtype=float)
        for i, wgt in enumerate(np.asarray(d, dtype=float)):
            if wgt != 0.0:
                row = row + wgt * Xf[i if i < irn else i + 1]
        return row

    def carrier_phasor(self, pss, output, carrier=1):
        """The `carrier`-th Fourier coefficient of the steady-state output.

        Computed here rather than taken from `fpss`, whose spectrum is RMS
        and energy-folded -- correct for reporting a magnitude and useless
        for a phasor, since folding discards the phase the AM/PM split is
        made of.
        """
        times, _X = pss.waveform
        row = self._output_waveform_row(pss, output)
        t = np.asarray(times, dtype=float)[:-1]
        v = row[:len(t)]
        w0 = 2.0 * np.pi / float(pss.period)
        return complex(np.sum(v * np.exp(-1j * carrier * w0 * t)) / len(t))

    def am_pm(self, pss, freq, output, carrier=1):
        """AM and PM modulation indices at `carrier`, per noise/signal source.

        Returns `(m_am, m_pm)`, each a row of length `m`: the modulation a
        unit source at reduced coordinate `i`, driven at `freq`, imposes on
        the `carrier`-th harmonic of the output.

        ⚠ TWO SOLVES AT ±freq, NOT ONE.  The upper sideband of harmonic `i`
        sits at `i f0 + freq` and the lower at `i f0 - freq`; with the
        convention that an input at `f` produces output at `f + l f0`,
        those are `H_i(freq)` and `H_i(-freq)`.  They are NOT conjugates of
        each other -- that would hold for an LTI circuit, and the whole
        point of an LPTV analysis is that it does not.  Taking one and
        conjugating it would silently force `m_pm = 0` or `m_am = 0`
        depending on which.

        ⚠ AND AN OSCILLATOR IS ALMOST PURE PM NEAR ITS CARRIER, which is
        the physical check: the phase response to a perturbation goes as
        `1/w_m` while the amplitude response stays bounded, so
        `|m_pm|/|m_am|` grows without bound as `freq -> 0`.  A
        decomposition that got the conjugate wrong gives a bounded ratio
        instead.

        ⚠ THE ABSOLUTE MAGNITUDE ON AN OSCILLATOR IS SMALL FOR A REASON
        (established 2026-09-08, the three-leg chain).  These are the p = 0
        band of `am_pm_noise`: a source at BASEBAND `freq` reaching the
        carrier sideband.  A baseband current moves the PHASE through the
        PPV's DC coefficient (Hajimiri-Lee's c_0), and a half-wave-symmetric
        orbit -- odd nonlinearity, `u(t + T/2) = -u(t)` -- has none, so on
        such a fixture the rows measure a symmetry zero (6e-9 .. 1e-13,
        proportional to 1/freq and to mu), the same zero the coloured
        up-conversion gate records for Gamma.  Breaking the symmetry
        (`_lc_osc(a)`) lifts |m_pm| at 1e-3 f0 from 1.2e-8 to 46.6 (a =
        0.05) and 231 (a = 0.25) -- linear in `a`.  The DIRECT rows (source
        at f0 + freq, sideband 0) agree with `pnoise` at every offset, and
        the split lands on the externally certified Lorentzian.  So do not
        read a small `am_pm` on a symmetric oscillator as a defect: it is
        the 1/f^3 up-conversion coefficient, and it is zero there.
        """
        C = self.carrier_phasor(pss, output, carrier)
        ## ⚠ RELATIVE TO THE SIGNAL, NOT AGAINST ZERO.  A harmonic the
        ## circuit does not produce still has a phasor of ~1e-16 rather
        ## than exactly 0, and dividing by it turns "there is no carrier
        ## here" into an enormous, confident modulation index.
        row = self._output_waveform_row(pss, output)
        scale = float(np.max(np.abs(row)))
        if abs(C) <= 1e-9 * max(scale, 1e-300):
            raise ValueError(
                'PAC.am_pm: the output carries no component at harmonic %d '
                '(|C| = %.3e against a signal scale of %.3e), so there is '
                'no carrier to modulate and AM/PM are not defined. Dividing '
                'by it would report a huge modulation of nothing. Pick a '
                'harmonic the circuit actually produces.'
                % (carrier, abs(C), scale))
        upper = self.adjoint_sideband_row(pss, freq, output, carrier)[0]
        lower = self.adjoint_sideband_row(pss, -freq, output, carrier)[0]
        return self.am_pm_indices(upper / C, lower / C)

    def am_pm_noise(self, pss, freq, output, carrier=1, maxsidebands=None,
                    modulated=False):
        """Output NOISE split into its AM and PM parts at `freq` from `carrier`.

        Returns `(S_am, S_pm, bands_used)`.  The two add to the noise in the
        pair of sidebands they decompose -- see the identity below -- and are in
        the same units as :meth:`pnoise`.

        ⚠ THIS NEEDS THE SIDEBAND *CORRELATION*, WHICH IS WHY IT IS NOT
        `|m_am|^2` FROM :meth:`am_pm`.  That method is the TRANSFER pair for a
        deterministic input; noise asks a different question, because whether
        the upper and lower sidebands are CORRELATED is exactly what decides the
        split.  Uncorrelated sidebands carry equal AM and PM -- the classical
        result for narrowband noise through an LTI system -- and it is the
        periodic operating point that correlates them.

        THE BAND BOOKKEEPING, which is the whole of the derivation and the one
        place a sign error would produce a plausible wrong answer.
        `adjoint_sideband_row(pss, g, output, l)` is the coefficient at output
        `g + l f0` for a unit source at `g`.  The two output sidebands sit at
        `carrier*f0 ± freq`, so a REAL noise band whose positive-frequency
        component is at `g = freq + p f0` reaches

            the UPPER output at `+g` through sideband `l = carrier - p`,
            the LOWER output at `-g` through sideband `l = carrier + p`,

        the second because a real process has `N(-g) = conj(N(g))` -- and that
        shared realisation IS the correlation.  Both contributions come from ONE
        band, so they are combined coherently; different `p` are different
        bands and are summed in power.  :meth:`am_pm` is exactly the `p = 0`
        term of this sum.

        The split per band is the same conjugate one :meth:`am_pm_indices`
        makes, `a + conj(b)` and `a - conj(b)` -- ⚠ the CONJUGATE, not `a ± b`:
        the sidebands counter-rotate about the carrier, and dropping it reports
        a rotating ellipse as pure AM.

        ⚠ THE GATE IS AN IDENTITY, NOT A TOLERANCE.  `pnoise` at the upper
        sideband folds precisely the bands `g = freq + p f0`, and at the lower
        precisely their negatives, so with the factor of one half below

            S_am + S_pm  ==  pnoise(carrier*f0 + freq) + pnoise(carrier*f0 - freq)

        exactly, because `|a+c|^2 + |a-c|^2 = 2|a|^2 + 2|c|^2` leaves no cross
        term.  A pairing error breaks it, which is what the test asserts.

        ⚠ THE AUTONOMOUS CAVEAT IS RETIRED (three-leg chain, 2026-09-08).  On a
        free-running oscillator this split sits on the SAME absolute scale as
        `pnoise` (identity to 1e-12 / 1e-16) and as the externally certified
        `oscillator_spectrum` (`S_pm = 4 S_v` at every offset: the PM content
        of the pair IS the Lorentzian, 2 S_v per sideband), with `S_am` rising
        from ~0 below the AM corner `f0/(2 pi Q_lambda)` to `S_pm` above it
        (⚠ this line said `f0/(4 pi Q)` until 2026-09-09 -- a factor of two
        the docs session caught against this very function: the ratio is an
        exact Lorentzian `u^2/(u_c^2 + u^2)` in `u = offset/f0` with
        `u_c = 1/(2 pi Q_lambda)`, `Q_lambda = -1/ln|lambda_2|`, half-power
        0.5007 there and 0.20 at the old corner, at Q = 8 and 16, 240 and
        480 points; the old formula OVERSTATED the AM content at every
        offset, 2.5x at its own corner) -- so the
        pair total is 4 S_v there and 8 S_v far out.  The "~1e-12 rows" were
        `am_pm`'s p = 0 band on a half-wave-symmetric fixture: a symmetry
        zero, see `am_pm`.  Oscillator magnitudes from this are trustworthy.
        """
        self._check_circuit(pss)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        N = len(fp.steps)
        f0 = 1.0 / float(fp.T)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)
        cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        k = int(carrier)
        S_am = 0.0
        S_pm = 0.0
        bands = []
        for p in range(-lmax, lmax + 1):
            g = float(freq) + p * f0
            a = self.adjoint_sideband_row(pss, g, output, k - p)[0]
            b = self.adjoint_sideband_row(pss, -g, output, k + p)[0]
            cy = cyfn(pss, 2.0 * np.pi * g)
            m_am = a + np.conj(b)
            m_pm = a - np.conj(b)
            S_am += 0.5 * float(np.real(m_am @ cy @ np.conj(m_am)))
            S_pm += 0.5 * float(np.real(m_pm @ cy @ np.conj(m_pm)))
            bands.append(p)
        return S_am, S_pm, bands

    def _deflated_solve(self, pss, alpha, b, transposed=False, tol=None):
        """`(I - alpha M) y = b` on an OSCILLATOR, with the pole taken out.

        ⚠ THE SINGULARITY IS THE ANSWER'S OWN POLE, NOT A NUMERICAL DEFECT,
        and that reframing is what makes the fix obvious.  At `alpha = 1`
        the operator is `I - M`, singular by the unit multiplier, and the
        solution really does diverge -- an oscillator's phase response to a
        perturbation goes as `1/df`.  What is wrong is COMPUTING a `1/eps`
        quantity through a system whose conditioning is also `1/eps`: the
        answer is genuinely large and the digits are genuinely gone.

        So the pole is carried ANALYTICALLY.  With `u`, `v` the right and
        left null vectors of `I - M` -- the orbit tangent and the PPV, both
        of which `ppv()` already returns -- border the system:

            [ I - alpha M   u ] [ w ]   [ b ]
            [     v^T       0 ] [ s ] = [ 0 ]

        Because `v^T (I - alpha M) = (1 - alpha) v^T` and `v^T w = 0`, the
        border variable comes out BOUNDED, `s = (v^T b)/(v^T u)`, with no
        `1/eps` in it.  And since `(I - alpha M) u = (1 - alpha) u`, the
        solution is recovered as

            y = w + s u / (1 - alpha)

        where `1 - alpha = 1 - exp(-j w T)` is the factor that vanishes at
        every harmonic, evaluated in closed form rather than inverted
        numerically.

        MEASURED on van der Pol, offsets from 0.3 down to 1e-9 of `f0`:

            offset/f0    sigma_min(plain)   sigma_min(bordered)
            3e-01           5.68e-01            1.17e-01
            1e-03           2.61e-03            2.04e-01
            1e-06           2.61e-06            2.04e-01
            1e-09           2.61e-09            2.04e-01

        The plain operator tracks the offset over nine decades; the
        bordered one is FLAT.  The two solutions agree to 5.7e-12 where the
        plain solve is still trustworthy, and their disagreement grows as
        `1/df` -- that is the PLAIN solve losing digits, not this one.

        ⚠ IT STILL DIVERGES AT AN EXACT HARMONIC, and it should: `1/(1 -
        alpha)` is then a division by zero, and the physical response is
        unbounded.  What changes is that every offset NEAR a harmonic is
        now well conditioned, which is where phase noise is measured.

        `transposed` solves `(I - alpha M^T) x = b`, whose null space is
        spanned by `v` and whose left null space is spanned by `u`, so the
        borders swap.
        """
        import scipy.sparse.linalg as spla
        fp = pss.factored_period()
        n = fp.width
        _v, info = pss.ppv()
        v = np.asarray(_v, dtype=float)
        u = np.asarray(info['tangent_pair'], dtype=float)
        vu = float(v @ u)
        if abs(vu) < 1e-300:
            raise ValueError(
                'PAC: the PPV is orthogonal to the orbit tangent, so the '
                'bordering is singular and the pole cannot be removed.')
        col, row = (v, u) if transposed else (u, v)
        mv = (fp.matvec_transposed if transposed else fp.matvec)
        b = np.asarray(b, dtype=complex).ravel()

        def _mv(z):
            z = np.asarray(z)
            w_, s_ = z[:n], z[n]
            top = w_ - alpha * np.asarray(mv(w_)) + s_ * col
            return np.concatenate((top, [complex(row @ w_)]))

        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=complex)
        rhs = np.concatenate((b, [0.0 + 0.0j]))
        rt = max(self.KRYLOV_FACTOR * pss.par.reltol if tol is None else tol,
                 1e-14)
        z = self._gmres_checked(A, rhs, rt, 'the deflated solve')
        w, s = z[:n], z[n]
        denom = 1.0 - alpha
        if denom == 0:
            raise ValueError(
                'PAC: the deflated solve was asked for an EXACT harmonic, '
                'where 1/(1 - alpha) is a division by zero and the physical '
                'response is unbounded. The pole is removed from the '
                'CONDITIONING, not from the answer.')
        return w + s * col / denom

    def _op(self, fp, alpha):
        """`v -> (I - alpha M) v`, never forming `M`."""
        return lambda v: np.asarray(v) - alpha * fp.matvec(v)

    def _solve_each(self, fp, alphas, rhs, tol):
        """One GMRES per frequency -- the baseline the sweep is measured against."""
        import scipy.sparse.linalg as spla
        n = fp.width
        ys, count = [], [0]

        for alpha, b in zip(alphas, rhs):
            op = self._op(fp, alpha)

            def _mv(v, _op=op):
                count[0] += 1
                return _op(v)

            A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
            ys.append(self._gmres_checked(A, b, tol, 'the m x m solve'))
        return ys, count[0]

    def _solve_subspace(self, fp, alphas, rhs, tol):
        """ONE Krylov subspace for the whole sweep -- the recycling.

        ⚠ THE SUBSPACE IS FREQUENCY-INDEPENDENT AND THAT IS THE WHOLE POINT.
        `A(alpha) = I - alpha M`, so

            span{r, A r, A^2 r, ...} = span{r, M r, M^2 r, ...}

        for every `alpha` -- Telichevesky et al.'s Theorem 1.  A basis of
        `M`'s Krylov space therefore serves every frequency, and each one
        costs a small dense least-squares over it instead of its own run of
        full-period replays.

        ⚠ WHAT IS NOT FREE is the right-hand side: `w(f)` genuinely differs
        per frequency, so a basis grown from one frequency's residual is not
        the space GMRES would have chosen for another.  This does not
        guess -- it minimises the TRUE residual over the shared span, checks
        it, and extends the basis (one matvec, kept for every later
        frequency) until every frequency is inside tolerance.  So the answer
        is never worse than the per-frequency solve; only the matvec count
        varies.
        """
        n = fp.width
        V = np.zeros((n, 0), dtype=complex)
        MV = np.zeros((n, 0), dtype=complex)
        count = [0]

        def extend(seed):
            """One Arnoldi step of `M` from `seed`, orthogonal to `V`."""
            nonlocal V, MV
            v = np.asarray(seed, dtype=complex).ravel().copy()
            if V.shape[1]:
                v = v - V @ (V.conj().T @ v)
                v = v - V @ (V.conj().T @ v)   ## reorthogonalise once
            nv = np.linalg.norm(v)
            if nv < 1e-14:
                return False
            v = v / nv
            count[0] += 1
            Mv = np.asarray(fp.matvec(v), dtype=complex)
            V = np.hstack((V, v[:, None]))
            MV = np.hstack((MV, Mv[:, None]))
            return True

        extend(rhs[0])
        ys = [None] * len(alphas)
        pending = list(range(len(alphas)))
        for _round in range(min(n, 200)):
            still = []
            for i in pending:
                alpha, b = alphas[i], rhs[i]
                AV = V - alpha * MV
                y, *_ = np.linalg.lstsq(AV, b, rcond=None)
                x = V @ y
                r = b - (x - alpha * (MV @ y))
                ys[i] = x
                nb = np.linalg.norm(b)
                if np.linalg.norm(r) > tol * (nb if nb else 1.0):
                    still.append((i, r))
            if not still:
                return ys, count[0]
            if V.shape[1] >= n:
                break
            ## grow the shared basis on the worst residual -- one matvec,
            ## and every frequency gets to use it
            worst = max(still, key=lambda p: np.linalg.norm(p[1]))
            if not extend(worst[1]):
                break
            pending = [i for i, _r in still]
        return ys, count[0]
