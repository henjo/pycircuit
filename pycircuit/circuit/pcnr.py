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
                     u_extra=None):
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

    nodemap = circuit.elementnodemap
    seen = set()
    for instance, element, _ra, _rb in junctions:
        if instance in seen:
            continue
        seen.add(instance)
        rows = nodemap[instance]
        sub = x[rows]
        g_mna[rows] -= np.asarray(element.i(sub, epar), dtype=float)
        J_mm[np.ix_(rows, rows)] -= np.asarray(element.G(sub, epar), dtype=float)

    J_ml = np.zeros((n, k))
    J_lm = np.zeros((k, n))
    g_lim = np.zeros(k)

    for idx, (instance, element, ra, rb) in enumerate(junctions):
        v = float(v_lim[idx])
        params = {p.name: getattr(element.iparv, p.name)
                  for p in element.instparams}
        i_terms = element.pcnr_i(v, params, epar, element.toolkit)
        di_terms = element.pcnr_didv(v, params, epar, element.toolkit)

        ## The device's current now enters through its OWN unknown, so its
        ## contribution to J_MNA/MNA is zero and all of it lands in J_MNA/lim.
        g_mna[ra] += float(i_terms[0])
        g_mna[rb] += float(i_terms[1])
        J_ml[ra, idx] += float(di_terms[0])
        J_ml[rb, idx] += float(di_terms[1])

        ## g_lim = v - (e_a - e_b)
        g_lim[idx] = v - (float(x[ra]) - float(x[rb]))
        J_lm[idx, ra] = -1.0
        J_lm[idx, rb] = +1.0

    return g_mna, g_lim, J_mm, J_ml, J_lm


def predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode):
    """One Newton step on the coupled system, by Schur complement.

        dx_MNA = -(J_mm - J_ml J_lm)^-1 (g_mna - J_ml g_lim)
        dx_lim = -(g_lim + J_lm dx_MNA)

    The matrix inverted is the size of the ORIGINAL MNA system: the ``lim/lim``
    block being the identity is what collapses the border into a rank-k update.
    """
    schur = J_mm - J_ml @ J_lm
    rhs = -(g_mna - J_ml @ g_lim)

    keep = [i for i in range(len(g_mna)) if i != irefnode]
    dx_r = np.linalg.solve(schur[np.ix_(keep, keep)], rhs[keep])
    dx_mna = np.zeros(len(g_mna))
    dx_mna[keep] = dx_r

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
        g_mna, g_lim, J_mm, J_ml, J_lm = augmented_system(
            circuit, x, v_lim, junctions, epar, u_extra=u)

        dx_mna, dx_lim = predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode)

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
