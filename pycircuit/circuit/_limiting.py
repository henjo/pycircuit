"""Newton limiting for exponential devices, shared by every device that has a junction.

STAGE 13(13-2).  This lived in `semiconductors.py` while `elements.Diode` carried
a second, incompatible implementation of the same idea -- one keeping private
`_vlim` state and returning `None`, the other state-free and returning a limited
vector.  Two conventions for one concept is how `SubCircuit.limit` came to accept
both, and the state-carrying one is what made the write-back bug invisible for as
long as it was: the only limiter in the tree was the one that could not expose it.

**This module imports nothing from the package**, so both `elements` and
`semiconductors` can use it without a cycle -- `semiconductors` imports the
package, and the package imports `elements`.
"""
import numpy as np


def _pnjlim(vnew, vold, VT, IS, toolkit):
    """SPICE's `pnjlim`: bound the per-iteration excursion of a junction voltage.

    Newton linearises an exponential, so a step taken from a lightly-biased point
    overshoots enormously -- and the model is then evaluated somewhere it has no
    business being.  `_expl` stops that becoming `nan`, but it does not stop the
    iteration wandering: with the clamp alone and no limiting, the common-emitter
    stage converges to a GENUINE but spurious operating point 200 V below the
    rails, with the base-collector junction forward biased at 0.91 V carrying the
    20 A the base resistor delivers.  That is a solution of the circuit equations;
    it is simply not the one anyone wants, and only limiting keeps Newton out of
    it.

    `vc` is the voltage at which the exponential's curvature starts to dominate.
    Above it the update is compressed logarithmically, which is what bounds the
    step without changing where the solution is -- the limiter only moves the
    point the next Jacobian is taken at, never the equations.
    """
    if IS <= 0.0:
        return vnew
    vc = VT * toolkit.log(VT / (IS * 1.414213562))
    if vnew > vc and vnew > 0.0:
        if vold > 0.0:
            arg = 1.0 + (vnew - vold) / VT
            return vold + VT * toolkit.log(arg) if arg > 0.0 else vc
        return VT * toolkit.log(vnew / VT)
    return vnew


def limit_junctions(x, x0, junctions, VT, IS_for, toolkit):
    """Apply `pnjlim` to each declared junction, returning a limited copy of `x`.

    `junctions` is a sequence of `(anode, cathode, move)` indices into the
    device's terminals.  `move` is not redundant: a BJT's two junctions share the
    base as anode, so limiting both by adjusting the anode would have the second
    undo the first.

    Returns a new vector rather than mutating device state.  That is what lets
    the limiter run on a traced backend, where device-private Python state cannot
    go, and it is the form PCNR's `refine` would need if stage 13 adopts it.
    """
    if not junctions:
        return x

    try:
        out = np.array(x, dtype=float, copy=True)
    except (TypeError, ValueError):
        ## Symbolic x: limiting is a numeric Newton aid and has nothing to
        ## contribute here.  Returning it untouched is correct, not a fallback.
        return x

    x0a = np.asarray(x0, dtype=float)
    for index, (anode, cathode, move) in enumerate(junctions):
        vnew = float(out[anode] - out[cathode])
        vold = float(x0a[anode] - x0a[cathode])
        vlim = _pnjlim(vnew, vold, VT, IS_for(index), toolkit)
        ## Reassign the moving terminal so the junction carries `vlim` exactly.
        if move == anode:
            out[anode] = out[cathode] + vlim
        else:
            out[cathode] = out[anode] - vlim
    return out
