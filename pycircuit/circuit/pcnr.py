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
    sequence of ``(anode, cathode)`` indices into its own terminals.  A device
    with several junctions (a BJT) contributes one entry per junction, in its
    declared order and contiguously; :func:`pcnr_devices` regroups them.
    """
    found = Junctions()
    elements = getattr(circuit, 'elements', None)
    if not elements:
        return found
    nodemap = circuit.elementnodemap
    for instance, element in elements.items():
        declared = getattr(element, 'pcnr_junctions', ())
        if not declared:
            continue
        ## STAGE 14-2 -- a multi-junction device is refused HERE, at discovery,
        ## if its junctions cannot reproduce its terminal voltages.
        validate = getattr(element, '_pcnr_map', None)
        if validate is not None:
            validate()
        rows = [int(r) for r in nodemap[instance]]
        k0 = len(found)
        for anode, cathode in declared:
            found.append((instance, element, rows[anode], rows[cathode]))
        found.devices.append((instance, element, rows, k0, len(found)))
    return found


class Junctions(list):
    """The junction list, plus the per-device table built alongside it.

    STAGE 14-1.  ``devices`` holds ``(instance, element, rows, k0, k1)``: the
    device's terminal rows, and the slice ``v_lim[k0:k1]`` of its own unknowns
    (a device's junctions are emitted contiguously).  Built once per solve,
    because regrouping on every Newton iteration made `augmented_system` 40%
    dearer on 60 diodes.  It is still a plain list of 4-tuples to everyone who
    iterates it.
    """

    def __init__(self, *args):
        super().__init__(*args)
        self.devices = []


def pcnr_devices(junctions, circuit=None):
    """STAGE 14-1 -- the junctions regrouped by device.

    Returns ``[(instance, element, rows, k0, k1), ...]``: the device's terminal
    rows and the slice ``[k0:k1]`` of ``junctions`` (and so of ``v_lim``) that
    it owns.  A device is evaluated ONCE per iteration, from the vector of all
    its own unknowns: a BJT's terminal currents depend on both junction
    voltages, so evaluating it junction by junction is not possible, and for a
    diode the vector has one entry and nothing changes.

    Read from :class:`Junctions` when the list came from :func:`pcnr_junctions`;
    rebuilt (needing ``circuit``) for a hand-made list.
    """
    devices = getattr(junctions, 'devices', None)
    if devices is not None:
        return devices
    nodemap = circuit.elementnodemap
    out = []
    for k, (instance, element, _ra, _rb) in enumerate(junctions):
        if out and out[-1][0] == instance:
            out[-1] = out[-1][:4] + (k + 1,)
        else:
            out.append((instance, element,
                        [int(r) for r in nodemap[instance]], k, k + 1))
    return out


def border_update(J, junctions, didv, g_lim=None, rhs_corr=None):
    """``J -= J_ml J_lm`` in place, in the rank-k form; optionally
    ``rhs_corr += J_ml g_lim``.

    STAGE 14-1 -- THE ONE PLACE THIS IS WRITTEN.  It has two consumers, the Schur
    reduction and the Jacobian the transient hands its step controller, and gate
    13-6 was a second consumer that open-coded it.  ``didv[k]`` is column ``k``
    of ``J_ml`` as ``(rows, vals)``: one entry per terminal of the device that
    owns junction ``k`` -- two for a diode, three for a BJT.  Row ``k`` of
    ``J_lm`` is always ``-1`` at the anode and ``+1`` at the cathode, so their
    outer product is two entries per terminal.
    """
    for idx, (_inst, _el, ra, rb) in enumerate(junctions):
        rows, vals = didv[idx]
        for r, val in zip(rows, vals):
            J[r, ra] += val
            J[r, rb] -= val
            if rhs_corr is not None:
                rhs_corr[r] += val * g_lim[idx]


def _params(element):
    ## Cached: rebuilding this per device per Newton iteration is 60 dict
    ## comprehensions over `getattr` on a 60-device circuit, every iteration.
    params = getattr(element, '_pcnr_params', None)
    if params is None:
        params = {p.name: getattr(element.iparv, p.name)
                  for p in element.instparams}
        element._pcnr_params = params
    return params


def has_charge(element):
    """STAGE 14-3 -- does this PCNR device store charge it can evaluate at its
    own unknowns?  (A charge model in `eval_q_pure` and the `pcnr_q` form.)"""
    return (hasattr(element, 'pcnr_q')
            and getattr(type(element), 'eval_q_pure', None) is not None)


def device_charge(q, v_lim, junctions, epar=defaultepar):
    """Add each charge-storing PCNR device's ``q(v_lim)`` into ``q`` in place.

    STAGE 14-3.  The caller assembles ``q`` WITHOUT the PCNR devices (``skip``)
    and calls this before ``get_diff``, so the companion current, the history and
    the LTE all see one charge vector -- with each device's charge where its
    current is, at its own unknowns.  At convergence ``v_lim`` equals the branch
    voltages and this is ``cir.q(x)``.  Returns whether any device had charge.
    """
    v_lim = np.asarray(v_lim, dtype=float)
    charged = False
    for _instance, element, rows, k0, k1 in pcnr_devices(junctions):
        if not has_charge(element):
            continue
        qd = element.pcnr_q(v_lim[k0:k1], _params(element), epar, element.toolkit)
        for r, val in zip(rows, np.asarray(qd, dtype=float).tolist()):
            q[r] += val
        charged = True
    return charged


def augmented_system(circuit, x, v_lim, junctions, epar=defaultepar,
                     u_extra=None, dense_blocks=True, J_extra=None,
                     charge_scale=None):
    """Residual and Jacobian blocks of the coupled ``[x_MNA ; x_lim]`` system.

    Returns ``(g_mna, g_lim, J_mm, J_ml, J_lm)``.  ``J_ll`` is not returned: it
    is the identity by construction, and materialising it would invite someone to
    use it.

    ``charge_scale`` -- STAGE 14-3 -- is the integrator's ``alpha0`` in
    ``geq = alpha0 * C``.  When given, each charge-storing device's
    ``alpha0 * dq/dv`` joins its ``J_ml`` columns; the caller must then have built
    ``J_extra`` from ``C`` assembled WITHOUT the PCNR devices, and their charge
    into ``u_extra`` via :func:`device_charge`.
    """
    n = circuit.n
    k = len(junctions)

    ## STAGE 13-8 -- ordinary assembly WITHOUT the PCNR devices; each is stamped
    ## below at its own `v_lim` instead of at the node voltage difference.
    ##
    ## This used to assemble them and subtract them out again.  That was exact in
    ## exact arithmetic and survived floating point only because `Diode` read a
    ## stale `_vlim` and so presented a tiny affine companion: the cancelled term
    ## was ~1e-13 A on the rectifiers.  An honest diode at an UNLIMITED node
    ## voltage -- 2.58e11 V on the first predict of a 1 A source -- overflows
    ## `exp`, and inf - inf is nan (benchmarks/transient_review/pcnr_narrow_fix.py).
    ## It also assembled every PCNR device twice per iteration, which was most of
    ## gate 13-4's cost.
    skip = frozenset(instance for instance, _e, _ra, _rb in junctions)
    g_mna = np.array(circuit.i(x, epar, skip=skip), dtype=float)
    J_mm = np.array(circuit.G(x, epar, skip=skip), dtype=float)
    if u_extra is not None:
        g_mna = g_mna + np.asarray(u_extra, dtype=float)
    if J_extra is not None:
        ## The transient's companion terms: `f = i(x) + iq + u`, `J = G + Geq`.
        ## They are passed in rather than recomputed because `solve_timestep`
        ## already has them and `get_diff` is not free.
        J_mm = J_mm + np.asarray(J_extra, dtype=float)

    ## A participating device with charge storage contributes to `iq`/`Geq`, and
    ## that contribution has to move to its own unknowns too.  With
    ## `charge_scale` it does (below, STAGE 14-3).  Without it -- or for a device
    ## that has charge but no `pcnr_q` -- it is REFUSED rather than silently left
    ## in the MNA block, where it would be evaluated at the node voltage while the
    ## current is evaluated at `v_lim`: the exact inconsistency PCNR exists to
    ## remove, reintroduced through the charge.
    if J_extra is not None:
        nodemap = circuit.elementnodemap
        for instance in sorted(skip):
            element = circuit.elements[instance]
            if has_charge(element) and charge_scale is not None:
                continue
            Csub = np.asarray(element.C(x[nodemap[instance]], epar), dtype=float)
            if np.any(np.abs(Csub) > 0.0):
                raise NotImplementedError(
                    '%r has charge storage and a PCNR junction; the charge term '
                    'would also have to move to the limited unknown, which %s'
                    % (instance, 'this caller did not do' if has_charge(element)
                       else 'this device cannot do (no pcnr_q)'))

    ## `dense_blocks=False` skips the (n,k) and (k,n) allocations entirely.
    ## `predict`'s sparse path never reads them -- each column of J_ml has one
    ## nonzero per device terminal and each row of J_lm has two, so the rank-k
    ## form carries the same information in `didv` plus the junction rows.
    ## Allocating and filling them anyway was 2 x n x k doubles per Newton
    ## iteration for nothing.
    J_ml = np.zeros((n, k)) if dense_blocks else None
    J_lm = np.zeros((k, n)) if dense_blocks else None
    g_lim = np.zeros(k)
    didv_list = [None] * k

    v_lim = np.asarray(v_lim, dtype=float)
    for instance, element, rows, k0, k1 in pcnr_devices(junctions, circuit):
        ## STAGE 14-1 -- the device is evaluated once, from ALL its own
        ## unknowns, and its current enters through them: its contribution to
        ## J_MNA/MNA is zero and all of it lands in J_MNA/lim.
        v = v_lim[k0:k1]
        params = _params(element)
        i_terms = element.pcnr_i(v, params, epar, element.toolkit)
        di = np.asarray(element.pcnr_didv(v, params, epar, element.toolkit),
                        dtype=float)
        for r, it in zip(rows, np.asarray(i_terms, dtype=float).tolist()):
            g_mna[r] += it
        if k1 - k0 == 1:
            ## The single-junction fast path, which is every diode.  It may
            ## return its one column as a 1-D (terminals,) array; building a
            ## (terminals, 1) array and transposing it per device per
            ## iteration cost ~10% of this function on 60 diodes.
            columns = (di.reshape(-1).tolist(),)
        else:
            columns = di.reshape(len(rows), k1 - k0).T.tolist()
        if charge_scale is not None and J_extra is not None and has_charge(element):
            ## STAGE 14-3 -- d(iq)/dv = alpha0 * dq/dv joins the current's column.
            dq = np.asarray(element.pcnr_dqdv(v, params, epar, element.toolkit),
                            dtype=float).reshape(len(rows), k1 - k0).T.tolist()
            columns = [[a + charge_scale * b for a, b in zip(col, dcol)]
                       for col, dcol in zip(columns, dq)]
        for local, vals in enumerate(columns):
            didv_list[k0 + local] = (rows, vals)
            if dense_blocks:
                for r, val in zip(rows, vals):
                    J_ml[r, k0 + local] += val

    for idx, (_instance, _element, ra, rb) in enumerate(junctions):
        ## g_lim = v - (e_a - e_b)
        g_lim[idx] = float(v_lim[idx]) - (float(x[ra]) - float(x[rb]))
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
        border_update(schur, junctions, didv, g_lim, rhs_corr)
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
    ``O(n^2 k)``. Pass ``junctions``/``didv`` to take that path.  (A BJT's
    column has three nonzeros, not two: still ``O(k)``.)
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
    for idx, (_instance, element, _ra, _rb) in enumerate(junctions):
        toolkit = element.toolkit
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
