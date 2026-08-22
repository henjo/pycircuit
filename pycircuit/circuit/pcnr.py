"""PCNR — Predictor/Corrector Newton-Raphson, Fang's replacement for limiting.

STAGE 13(13-3).  K. V. Aadithya, E. R. Keiter, T. Mei (Sandia), *"Predictor/
Corrector Newton-Raphson (PCNR): A Simple, Flexible, Scalable, Modular, and
Consistent Replacement for Limiting in Circuit Simulation"*.

**The idea.** Classic limiting moves a NODE voltage so a junction sees a bounded
value. That node is not the device's to move -- everything else attached sees the
change -- which is why gate 13-2 could not make it consistent. PCNR instead gives
each limited quantity its own unknown::

    x = [x_MNA ; x_lim]        g = [g_MNA ; g_lim]
    g_lim,k = v_k - (e_a - e_b)

and has the device evaluate its current at **its own** ``v_k`` rather than at
``e_a - e_b``. Nothing is shared, so nothing can clash.

**Why this needs no change to MNA assembly.** Element stamps are additive, so the
"everything except the PCNR devices" part of the residual and Jacobian is just
the ordinary assembly minus each participating device's own stamp. That is what
:func:`augmented_system` does, and it is the reason this module is a layer rather
than a rewrite.

The Jacobian's ``lim/lim`` block is the identity, which is what makes the Schur
elimination in :func:`predict` cost one extra solve against a matrix of the
ORIGINAL size -- the paper's third key idea.
"""
import numpy as np

from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit._limiting import _pnjlim


def pcnr_junctions(circuit):
    """Every PCNR-participating junction in the circuit.

    Returns a list of ``(instance, device, row_anode, row_cathode)`` in circuit
    row coordinates.  A device opts in by declaring ``pcnr_junctions`` as a
    sequence of ``(anode, cathode)`` indices into its own terminals.
    """
    found = []
    elements = getattr(circuit, 'elements', None)
    if not elements:
        return found
    nodemap = circuit.elementnodemap
    for instance, element in elements.items():
        for anode, cathode in getattr(element, 'pcnr_junctions', ()):
            rows = nodemap[instance]
            found.append((instance, element, int(rows[anode]), int(rows[cathode])))
    return found


def augmented_system(circuit, x, v_lim, junctions, epar=defaultepar,
                     u_extra=None, dense_blocks=True, J_extra=None):
    """Residual and Jacobian blocks of the coupled ``[x_MNA ; x_lim]`` system.

    Returns ``(g_mna, g_lim, J_mm, J_ml, J_lm)``.  ``J_ll`` is not returned: it
    is the identity by construction, and materialising it would invite someone to
    use it.
    """
    n = circuit.n
    k = len(junctions)

    ## Ordinary assembly, then remove each PCNR device's own contribution -- the
    ## device is about to be re-stamped at its own `v_lim` instead of at the node
    ## voltage difference.
    g_mna = np.array(circuit.i(x, epar), dtype=float)
    J_mm = np.array(circuit.G(x, epar), dtype=float)
    if u_extra is not None:
        g_mna = g_mna + np.asarray(u_extra, dtype=float)
    if J_extra is not None:
        ## The transient's companion terms: `f = i(x) + iq + u`, `J = G + Geq`.
        ## They are passed in rather than recomputed because `solve_timestep`
        ## already has them and `get_diff` is not free.
        J_mm = J_mm + np.asarray(J_extra, dtype=float)

    nodemap = circuit.elementnodemap
    seen = set()
    for instance, element, _ra, _rb in junctions:
        if instance in seen:
            continue
        seen.add(instance)
        rows = nodemap[instance]
        sub = x[rows]
        ## CHARGE STAYS IN THE MNA BLOCK, at the node voltage; only the
        ## RESISTIVE current moves to the limited unknown.  This used to be
        ## refused outright, on the grounds that leaving it here reintroduces
        ## "the exact inconsistency PCNR exists to remove".  That reasoning
        ## conflated two different things, and the refusal cost every junction
        ## device with capacitance -- which is to say every real one.
        ##
        ## What PCNR removes is a CLASH BETWEEN DEVICES over a shared
        ## linearisation point: two diodes on one branch each limiting
        ## `e_a - e_b`, the second undoing the first, the outcome depending on
        ## visit order.  A charge is not limited by anyone, so it has no
        ## clash to remove and no owner to fight over.
        ##
        ## And the Newton system stays EXACTLY consistent: `g_mna`'s charge
        ## part is a function of `x_MNA` and its derivative `Geq` is in
        ## `J_MNA/MNA` where it belongs, while the device's current is a
        ## function of `v_lim` alone and its derivative is in `J_MNA/lim`.
        ## `J` is `dg/dx` for both blocks; nothing is taken at a different
        ## point from the thing it differentiates.
        ##
        ## What IS given up is that the reactive part is not limited along
        ## the iteration path.  At the answer that costs nothing: convergence
        ## already requires `|g_lim|` below tolerance (`_solve_timestep_pcnr`
        ## checks it), i.e. `v_lim == e_a - e_b`, so the charge was evaluated
        ## at the voltage the current converged to.  The paper does not derive
        ## the reactive case either -- its footnote 1 says PCNR "works for
        ## differential-algebraic equations as well, but for simplicity, we
        ## only consider algebraic equations".
        ##
        ## Proven, not argued: `test_pcnr_charge.py` finite-differences the
        ## augmented residual to confirm `J_eff == df_eff/dx` with a
        ## charge-storing participant, and compares pcnr on/off transients.
        g_mna[rows] -= np.asarray(element.i(sub, epar), dtype=float)
        J_mm[np.ix_(rows, rows)] -= np.asarray(element.G(sub, epar), dtype=float)

    ## `dense_blocks=False` skips the (n,k) and (k,n) allocations entirely.
    ## `predict`'s sparse path never reads them -- each column of J_ml has two
    ## nonzeros and each row of J_lm has two, so the rank-one form carries the
    ## same information in `didv` plus the junction rows.  Allocating and filling
    ## them anyway was 2 x n x k doubles per Newton iteration for nothing.
    J_ml = np.zeros((n, k)) if dense_blocks else None
    J_lm = np.zeros((k, n)) if dense_blocks else None
    g_lim = np.zeros(k)
    didv_list = []

    ## A device may own MORE THAN ONE limited quantity -- a BJT's two
    ## junctions, say -- so the device is told WHICH of its own junctions
    ## is being asked about.  `pcnr_junctions` lists a device's pairs in
    ## declaration order, so counting occurrences as we go gives each its
    ## local index.  Single-junction devices never see anything but 0 and
    ## may ignore the argument, which is why it is keyword-with-default.
    _seen_of = {}

    for idx, (instance, element, ra, rb) in enumerate(junctions):
        jn = _seen_of.get(instance, 0)
        _seen_of[instance] = jn + 1
        v = float(v_lim[idx])
        ## Cached: rebuilding this per junction per Newton iteration is 60 dict
        ## comprehensions over `getattr` on a 60-device circuit, every iteration.
        params = getattr(element, '_pcnr_params', None)
        if params is None:
            params = {p.name: getattr(element.iparv, p.name)
                      for p in element.instparams}
            element._pcnr_params = params
        i_terms = element.pcnr_i(v, params, epar, element.toolkit, jn=jn)
        di_terms = element.pcnr_didv(v, params, epar, element.toolkit,
                                     jn=jn)

        ## The device's current now enters through its OWN unknown, so its
        ## contribution to J_MNA/MNA is zero and all of it lands in J_MNA/lim.
        g_mna[ra] += float(i_terms[0])
        g_mna[rb] += float(i_terms[1])
        if dense_blocks:
            J_ml[ra, idx] += float(di_terms[0])
            J_ml[rb, idx] += float(di_terms[1])
        didv_list.append((float(di_terms[0]), float(di_terms[1])))

        ## g_lim = v - (e_a - e_b)
        g_lim[idx] = v - (float(x[ra]) - float(x[rb]))
        if dense_blocks:
            J_lm[idx, ra] = -1.0
            J_lm[idx, rb] = +1.0

    return g_mna, g_lim, J_mm, J_ml, J_lm, didv_list


def schur_reduce(g_mna, g_lim, J_mm, J_ml=None, J_lm=None, junctions=None,
                 didv=None):
    """The augmented system collapsed onto the original MNA size.

    Returns ``(f_eff, J_eff)`` such that ``J_eff dx = -f_eff`` is the MNA part of
    one Newton step on ``[x_MNA ; x_lim]``::

        J_eff = J_mm - J_ml J_lm        f_eff = g_mna - J_ml g_lim

    **This is the only place that formula is written.** It has two consumers with
    different needs -- :func:`predict`, which solves with it, and the transient
    step controller, which only wants the matrix -- and when the second one
    open-coded it instead, it got `cir.G(x)` and a step count 6.6x too large.
    A shared definition is the mechanism that stops that recurring, which matters
    more here than the handful of lines it saves.

    ``J_ll`` being the identity is what makes the border collapse to a rank-k
    update; with ``junctions``/``didv`` that is exploited directly, at ``O(k)``
    rather than the ``O(n^2 k)`` of a dense ``J_ml @ J_lm``.
    """
    if junctions is not None and didv is not None:
        schur = np.array(J_mm, copy=True)
        rhs_corr = np.zeros(len(g_mna))
        for idx, (_inst, _el, ra, rb) in enumerate(junctions):
            dia, dib = didv[idx]
            ## column k of J_ml is (dia at ra, dib at rb); row k of J_lm is
            ## (-1 at ra, +1 at rb).  Their outer product is these four entries.
            schur[ra, ra] += dia
            schur[ra, rb] -= dia
            schur[rb, ra] += dib
            schur[rb, rb] -= dib
            rhs_corr[ra] += dia * g_lim[idx]
            rhs_corr[rb] += dib * g_lim[idx]
        return g_mna - rhs_corr, schur
    return g_mna - J_ml @ g_lim, J_mm - J_ml @ J_lm


def predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode, junctions=None,
            didv=None):
    """One Newton step on the coupled system, by Schur complement.

        dx_MNA = -(J_mm - J_ml J_lm)^-1 (g_mna - J_ml g_lim)
        dx_lim = -(g_lim + J_lm dx_MNA)

    The matrix factorised is the size of the ORIGINAL MNA system: the ``lim/lim``
    block being the identity is what collapses the border into a rank-k update.

    **That collapse has to be exploited, not merely stated.** Forming
    ``J_ml @ J_lm`` as a dense ``(n,k)·(k,n)`` product costs ``O(n^2 k)`` and
    measured **+62% per iteration** against classic limiting on 60 diodes -- a
    gate-13-4 failure that was entirely the implementation's. Each column of
    ``J_ml`` has two nonzeros and each row of ``J_lm`` has two, so the product is
    a sum of ``k`` rank-one terms touching four entries each: ``O(k)`` work, not
    ``O(n^2 k)``. Pass ``junctions``/``didv`` to take that path.
    """
    f_eff, schur = schur_reduce(g_mna, g_lim, J_mm, J_ml, J_lm, junctions, didv)
    rhs = -f_eff

    keep = [i for i in range(len(g_mna)) if i != irefnode]
    dx_r = np.linalg.solve(schur[np.ix_(keep, keep)], rhs[keep])
    dx_mna = np.zeros(len(g_mna))
    dx_mna[keep] = dx_r

    if junctions is not None:
        dx_lim = np.empty(len(g_lim))
        for idx, (_inst, _el, ra, rb) in enumerate(junctions):
            dx_lim[idx] = -(g_lim[idx] + (-dx_mna[ra] + dx_mna[rb]))
    else:
        dx_lim = -(g_lim + J_lm @ dx_mna)
    return dx_mna, dx_lim


def refine(junctions, v_old, v_new, epar=defaultepar):
    """The CORRECT phase: each device limits only the variables it owns.

    This is where PCNR differs from limiting in the way that matters. Nothing is
    shared, so applying one device's limiter cannot disturb another's -- the
    clash section 2 of the paper describes cannot arise, by construction rather
    than by ordering.
    """
    out = np.array(v_new, dtype=float, copy=True)
    _seen_of = {}
    for idx, (_instance, element, _ra, _rb) in enumerate(junctions):
        jn = _seen_of.get(_instance, 0)
        _seen_of[_instance] = jn + 1
        toolkit = element.toolkit
        ## THE DEVICE OWNS THE LAW FOR ITS OWN QUANTITY.  PCNR supplies the
        ## architecture -- one unknown per limited quantity, nothing shared,
        ## so no device's limiter can disturb another's -- and the device
        ## supplies the limiter, which is the paper's modularity claim.  A
        ## device that declares no `pcnr_limit` gets SPICE's `pnjlim` with
        ## its `IS`, which is what every junction device wants and what the
        ## hand-written `Diode` has always used.
        limiter = getattr(element, 'pcnr_limit', None)
        if limiter is not None:
            params = getattr(element, '_pcnr_params', None)
            if params is None:
                params = {q.name: getattr(element.iparv, q.name)
                          for q in element.instparams}
                element._pcnr_params = params
            out[idx] = float(limiter(float(v_new[idx]), float(v_old[idx]),
                                     params, epar, toolkit, jn=jn))
            continue
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        IS = getattr(element.iparv, 'IS', 0.0)
        out[idx] = _pnjlim(float(v_new[idx]), float(v_old[idx]), VT, IS, toolkit)
    return out


def solve_dc(circuit, refnode, x0=None, epar=defaultepar, maxiter=200,
             reltol=1e-6, abstol=1e-12):
    """DC operating point by PCNR.  Returns ``(x, v_lim, iterations)``.

    The flow is Figure 2's: initialise, then repeat predict-then-correct until
    converged.  There is no limiting of the MNA vector anywhere -- the only thing
    limited is each device's own unknown, in :func:`refine`.
    """
    junctions = pcnr_junctions(circuit)
    if not junctions:
        raise ValueError('no device in this circuit declares a PCNR junction')

    n = circuit.n
    irefnode = circuit.get_node_index(refnode)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float, copy=True)
    x[irefnode] = 0.0

    ## Initialise each device's unknown from the branch voltage it stands for --
    ## the paper's "independent initialization by different devices".
    v_lim = np.array([float(x[ra] - x[rb]) for _i, _e, ra, rb in junctions])

    u = np.asarray(circuit.u(0.0, epar, analysis='dc'), dtype=float)

    for it in range(maxiter):
        g_mna, g_lim, J_mm, J_ml, J_lm, didv = augmented_system(
            circuit, x, v_lim, junctions, epar, u_extra=u, dense_blocks=False)

        dx_mna, dx_lim = predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode,
                                 junctions=junctions, didv=didv)

        x_new = x + dx_mna
        x_new[irefnode] = 0.0
        v_new = v_lim + dx_lim

        ## CORRECT: each device limits only what it owns.
        v_new = refine(junctions, v_lim, v_new, epar)

        converged = (np.max(np.abs(dx_mna)) < reltol * np.max(np.abs(x_new)) + abstol
                     and np.max(np.abs(g_lim)) < reltol * max(np.max(np.abs(v_new)), 1.0) + abstol)
        x, v_lim = x_new, v_new
        if converged:
            return x, v_lim, it + 1

    raise RuntimeError(
        'PCNR did not converge in %d iterations: max|g_lim| = %g'
        % (maxiter, float(np.max(np.abs(g_lim)))))
