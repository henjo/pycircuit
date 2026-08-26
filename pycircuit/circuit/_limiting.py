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

    **A step of at most `2*VT` passes through unlimited**, per Listing 1 of
    Aadithya, Keiter & Mei (the PCNR paper, doi 10.1007/978-3-030-44101-2_19):
    `if vnew <= vc or abs(vnew - vold) <= 2.0*VT`.  This module previously
    omitted that escape and compressed every step taken above `vc`, however
    small.  The compression of a small step is nearly the identity
    (`log(1+e) ~ e`), so no operating point moved -- but "nearly" is the point:
    the paper's form makes limiting a STRICT no-op near the solution, which is
    what lets a simulator use "was limiting applied?" as a convergence signal
    (the paper's footnote 3) instead of watching an O(e^2) perturbation.

    The `vnew > 0.0` guard is a DELIBERATE addition to the listing, kept: with a
    large `IS`, `vc` goes negative, and the paper's `vold <= 0` branch then
    evaluates `log(vnew/VT)` at a negative `vnew` and returns nan.
    """
    if IS <= 0.0:
        return vnew
    vc = VT * toolkit.log(VT / (IS * 1.414213562))
    if abs(vnew - vold) <= 2.0 * VT:
        return vnew
    if vnew > vc and vnew > 0.0:
        if vold > 0.0:
            arg = 1.0 + (vnew - vold) / VT
            return vold + VT * toolkit.log(arg) if arg > 0.0 else vc
        return VT * toolkit.log(vnew / VT)
    return vnew


def _pnjlim_branchless(vnew, vold, VT, IS, log):
    """`_pnjlim`, as pure operators (P19 stage 2 -- the traced PCNR refine).

    Same law, arithmetic selects instead of Python `if`s, so it runs on numpy
    scalars/arrays and traced jax arrays alike.  `log` is passed in to keep
    this module import-free, as its header requires.  Guards: every `log`
    argument is floored to a tiny positive value BEFORE the call; the floored
    branches are only ever SELECTED when their condition made the argument
    positive, so the floor changes no selected value.  `IS <= 0` (no junction)
    passes vnew through, as the branching original does.
    """
    tiny = 1e-300

    def _floor(a):
        ## Select-based positive floor.  NOT (a+b+|a-b|)/2: that identity
        ## cancels to exactly 0.0 for large-negative `a`, and log(0) = -inf
        ## then rides through the arithmetic selects as 0*inf = nan.
        pos = a > tiny
        return a * pos + tiny * (1 - pos)

    is_pos = IS > 0.0
    vc = VT * log(_floor(VT / (_floor(IS) * 1.414213562)))
    ## The paper's `abs(vnew - vold) <= 2*VT` escape, as a select -- the two
    ## implementations must agree bit for bit (see test_limiting parity).
    small_step = abs(vnew - vold) <= 2.0 * VT
    active = (vnew > vc) * (vnew > 0.0) * (1 - small_step)
    pos_old = vold > 0.0
    arg = 1.0 + (vnew - vold) / VT
    arg_pos = arg > 0.0
    r_old = vold + VT * log(_floor(arg))
    r_old = arg_pos * r_old + (1 - arg_pos) * vc
    r_new = VT * log(_floor(vnew / VT))
    limited = pos_old * r_old + (1 - pos_old) * r_new
    out = active * limited + (1 - active) * vnew
    return is_pos * out + (1 - is_pos) * vnew


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


def _fetlim(vnew, vold, vto, toolkit):
    """SPICE's `fetlim`: bound the per-iteration excursion of a FET gate voltage.

    Transcribed from **ngspice-47**, `src/spicelib/devices/devsup.c`,
    `DEVfetlim` (the same routine as SPICE3f5's `Include/devsup.c`, with the
    one change noted below).  The law is a piecewise *step* bound rather than
    a logarithmic compression: it is not an exponential that is being tamed
    here but a threshold, so what matters is how far the point may move
    RELATIVE TO `vto` in one iteration, and on which side of `vto` it starts.

    Three regimes, on `(vold, delv = vnew - vold)`:

    * `vold >= vto + 3.5` -- strongly on.  Going further on, a step is capped
      at `vtsthi = 2*|vold - vto| + 2`; coming off far enough to cross back
      under `vto + 3.5`, the new point lands no lower than `vto + 2`.
    * `vto <= vold < vto + 3.5` -- the middle region.  The new point is simply
      confined to `[vto - 0.5, vto + 4]`.
    * `vold < vto` -- off.  Going further off, capped at `vtsthi`; turning on,
      never past `vto + 0.5` in one step.

    The bounds are never tighter than 2 V (`vtsthi >= 2`), and the middle
    region's clamps sit 0.5 V and 4 V from `vto`, so **every step of 0.45 V
    or less passes through EXACTLY**, whatever `vold` and `vto` are -- the
    limiter is a strict no-op near the solution, the same property
    `_pnjlim` has from its `2*VT` escape, obtained here by a different
    mechanism.  (0.45 is asserted, not asymptotic: see
    `test_fetlim_is_a_strict_no_op_near_the_solution`.)

    **`vtstlo` does not appear above because both of its branches are
    unreachable.**  The C computes it and tests it twice; neither test can
    fire.  Off, turning on, the branch needs `vnew <= vto + 0.5` with
    `vold < vto`, so `delv <= (vto - vold) + 0.5`, while
    `vtstlo = (vto - vold) + 1`.  On, coming off, it needs `vnew >= vto +
    3.5` with `vold >= vto + 3.5`, so `-delv <= (vold - vto) - 3.5`, again
    below `vtstlo`.  That matters for reading the two sources: ngspice's
    header credits Alan Gillespie with "a new definition for vtstlo"
    (`|vold-vto| + 1` where SPICE3f5 had `vtsthi/2 + 2 = |vold-vto| + 3`),
    and the two are indistinguishable in every input -- checked
    algebraically above and by exhaustive sweep in
    `test_fetlim_vtstlo_is_dead_code_in_devfetlim`.  So the ngspice-47
    lineage claimed here costs nothing and buys nothing; it is stated
    because a reader will find the difference in the C and deserves to
    know it is inert.

    Note that `_pnjlim` above does NOT follow ngspice -- it follows the PCNR
    paper's listing, which is the SPICE3f5 `pnjlim`, and ngspice's `pnjlim`
    IS materially different (Gillespie's negative-voltage limiting).  The two
    functions have different lineages, each stated where it is.

    `vto` is SPICE's `von`, the actual turn-on voltage -- which in a real
    MOSFET model is BIAS DEPENDENT (body effect), recomputed each iteration.
    Nothing in this function assumes otherwise; the restriction that the DSL
    can only pass a parameter-level constant is `hdl.limit_fet`'s, and is
    documented there.

    `toolkit` is unused -- the law is `abs`/`min`/`max` and no transcendental
    -- and is kept in the signature so every limiter in this module is called
    the same way.
    """
    vtsthi = abs(2.0 * (vold - vto)) + 2.0
    ## Dead in both of its uses below -- kept so this stays a transcription
    ## of the C rather than an edit of it.  See the docstring for why.
    vtstlo = abs(vold - vto) + 1.0
    vtox = vto + 3.5
    delv = vnew - vold

    if vold >= vto:
        if vold >= vtox:
            if delv <= 0.0:
                ## Going off.
                if vnew >= vtox:
                    if -delv > vtstlo:        # unreachable
                        vnew = vold - vtstlo
                else:
                    vnew = max(vnew, vto + 2.0)
            else:
                ## Staying on.
                if delv >= vtsthi:
                    vnew = vold + vtsthi
        else:
            ## Middle region.
            if delv <= 0.0:
                vnew = max(vnew, vto - 0.5)
            else:
                vnew = min(vnew, vto + 4.0)
    else:
        ## Off.
        if delv <= 0.0:
            if -delv > vtsthi:
                vnew = vold - vtsthi
        else:
            vtemp = vto + 0.5
            if vnew <= vtemp:
                if delv > vtstlo:             # unreachable
                    vnew = vold + vtstlo
            else:
                vnew = vtemp
    return vnew


def _limvds(vnew, vold, toolkit):
    """SPICE's `limvds`: bound the per-iteration excursion of a FET's `vds`.

    Transcribed from **ngspice-47**, `src/spicelib/devices/devsup.c`,
    `DEVlimvds` -- byte for byte the same routine as SPICE3f5's, unlike
    `fetlim`.  For `vold >= 3.5` a rising `vds` may at most triple (`3*vold +
    2`) and a falling one may not drop below 2 while it is still above 3.5;
    below 3.5, `vds` is confined to `[-0.5, 4]`.

    This is what stops a saturating model being evaluated at hundreds of
    volts.  `compact.py`'s `PspMosLongChannel.limit` records what happens
    without it: `dIds/dVds` falls to 1e-11 by 500 V, so a solver that lands
    out there has a numerically empty row and is reported singular rather
    than slow.  That method is worth reading beside this one, and worth
    reading carefully: it names `DEVfetlim` as the part it plays, but it
    bounds d, g and b alike by a single symmetric `|delta| <= vlimit = 1 V`
    about the source, which is neither `fetlim` nor `limvds` -- it is a
    third, blunter law that happens to cover both jobs.  The failure its
    docstring describes is specifically the one THIS function exists for.

    **The `vold < 0` fold is a DELIBERATE addition**, and it is not invented:
    `DEVlimvds` alone assumes a forward-mode device, and clamping at `-0.5`
    would drag any reverse-biased `vds` up to -0.5 -- so the routine is NOT a
    no-op near a solution with `vds < -0.5`, which is a limiter's one
    indispensable property.  SPICE handles that at the device level rather
    than in `devsup.c`: `mos1load.c` selects on the sign of the previous
    `vds` and calls `-DEVlimvds(-vds, -vds_old)` in the reverse branch.  The
    fold here is exactly that call, moved inside so that a per-probe
    declaration -- which cannot see the device's mode -- still gets it.  With
    the fold the routine is a strict no-op for a small step about any `vold`.

    `toolkit` is unused; see `_fetlim`.
    """
    if vold < 0.0:
        ## Reverse mode: mos1load.c's `vds = -DEVlimvds(-vds, -vds_old)`.
        return -_limvds(-vnew, -vold, toolkit)

    if vold >= 3.5:
        if vnew > vold:
            vnew = min(vnew, 3.0 * vold + 2.0)
        elif vnew < 3.5:
            vnew = max(vnew, 2.0)
    else:
        if vnew > vold:
            vnew = min(vnew, 4.0)
        else:
            vnew = max(vnew, -0.5)
    return vnew


def _deltalim(vnew, vold, vmax, toolkit):
    """Bound the per-iteration EXCURSION of a branch: ``|vnew - vold| <= vmax``.

    The blunt law, and the one with no device physics in it.  Where
    `pnjlim` compresses a junction logarithmically and `fetlim` works
    about a threshold, this only says *do not move this branch more than
    `vmax` volts in one step* -- so it applies to any branch at all,
    including one whose model has no notion of a critical voltage.

    **It is not invented for the DSL.**  `compact.PspMosLongChannel.limit`
    already does exactly this by hand, and its docstring says why: with
    the drain voltage limited to `Vdsat` the current saturates, so
    ``dIds/dVds`` falls to 1e-11 by 500 V and 1e-28 by 1e7 -- a solver
    that lands out there has a numerically empty row and is reported
    singular rather than slow.  That method bounds d, g and b alike by a
    symmetric ``|delta| <= 1 V`` about the source.  This is that law,
    named and shared.

    Two properties the other kinds also have, and which the write-back
    depends on:

    * **it returns `vnew` ITSELF when it does not bite**, not a copy and
      not ``vold + (vnew - vold)`` -- the ordinary path decides "did
      limiting fire?" by comparing `vlim == vnew`, and
      `pcnr.limit_block` by `lim != raw`.  Reconstructing an equal value
      would report a bite that did not happen;
    * it is **symmetric**, so it cannot introduce a direction preference
      into a model that has none.

    `vmax` must be positive; a non-positive bound would pin the branch at
    `vold` for ever, which is a stall rather than a limit, so it is
    rejected where it is declared (:func:`pycircuit.circuit.hdl.limit_delta`)
    rather than here on the hot path.
    """
    d = vnew - vold
    if d > vmax:
        return vold + vmax
    if d < -vmax:
        return vold - vmax
    return vnew


def apply_limit(kind, vnew, vold, pars, toolkit):
    """Dispatch one limiter KIND on its own parameter list.

    One place where a `kind` string becomes a call, shared by the per-probe
    write-back and the device-level one (`limit_together`, roadmap 10.3(b)).
    It existed twice inside the generated `limit()` -- once for ranking and
    once for the real pass -- and a third copy for the device path is how the
    two spellings of `pnjlim`'s argument order (`pars` is `(IS, VT)`; the
    function takes `(VT, IS)`) would eventually have diverged.
    """
    if kind == 'pnj':
        return _pnjlim(vnew, vold, pars[1], pars[0], toolkit)
    if kind == 'fet':
        return _fetlim(vnew, vold, pars[0], toolkit)
    if kind == 'vds':
        return _limvds(vnew, vold, toolkit)
    if kind == 'delta':
        return _deltalim(vnew, vold, pars[0], toolkit)
    if kind == 'id':
        ## `limit_identity`: a limited unknown with no law.  Returning
        ## `vnew` ITSELF (not a copy, not `vold + (vnew - vold)`) is what
        ## makes "did not bite" exact on both paths: the ordinary
        ## write-back compares `vlim == vnew` and `pcnr.limit_block`
        ## compares `lim != raw`.
        return vnew
    raise ValueError('unknown limiter kind %r' % (kind,))    # pragma: no cover


def device_writeback(out, targets, drift, pinned=()):
    """Write a whole DEVICE's limited branch voltages back at once.

    `targets` is `(row_plus, row_minus, v_unlimited, v_limited)` per probe,
    all four already numbers; `drift` is `|x - x0|` per row; `pinned` names
    rows another probe has already written.  `out` is mutated in place and
    the set of rows written is returned.

    **This is the device-level generalisation of the per-probe rule, and
    what it removes is the competition.**  The per-probe write-back moves
    exactly one endpoint per probe, so two probes sharing a terminal have to
    negotiate for it and the loser is pushed onto a node it would not have
    chosen (roadmap 12.1).  Here the probes are a GRAPH -- nodes are
    terminals, edges are branch-voltage constraints -- and a graph with no
    cycle is satisfiable exactly:

      1. take the constraints in order of DECREASING correction, adding each
         to a spanning forest and dropping any that closes a cycle.  A cycle
         is three probes whose limited values need not sum to zero
         (`vgs - vds - vgd != 0`), so one of them cannot be honoured; the
         one dropped is the one asking for the least;
      2. in each component, hold the node that has drifted LEAST from the
         last accepted point and derive every other node from it.  That is
         the same judgement the per-probe rule makes ("move the terminal
         that drifted further") stated for a whole tree instead of one edge,
         and it keeps the property that made the runtime choice necessary:
         a bounded branch voltage is only ever added to the SANEST node in
         the device, never to whichever endpoint the declaration named.

    Every kept probe therefore ends carrying exactly its own limited value,
    with no probe undone and no ordering between them -- the forest and the
    anchors are functions of the data, so reversing the probe list changes
    nothing.

    **If no probe bit, nothing is written at all** -- not even a round-trip
    `out[a] = out[b] + (out[a] - out[b])`, which is not the identity in
    floating point and would cost the "a small step passes through exactly"
    property that lets `did limiting fire?` be a convergence signal.  A
    probe that did not bite is still a CONSTRAINT when some other probe did:
    that is how the device holds `vds` at the value Newton asked for while
    it compresses `vgs`, which a per-probe limiter cannot do -- it would
    move the shared source and let `vds` follow.
    """
    if not any(t[3] != t[2] for t in targets):
        return set()

    parent = {}

    def _find(n):
        parent.setdefault(n, n)
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    ## Kruskal on |correction|, ties by row index -- data-derived, so the
    ## result does not depend on the order the probes were declared in.
    adj, nodes = {}, set()
    for i in sorted(range(len(targets)),
                    key=lambda k: (-abs(targets[k][3] - targets[k][2]),
                                   targets[k][0], targets[k][1])):
        ra, rb, _vn, vl = targets[i]
        ka, kb = _find(ra), _find(rb)
        if ka == kb:
            continue
        parent[ka] = kb
        adj.setdefault(ra, []).append((rb, -vl))
        adj.setdefault(rb, []).append((ra, +vl))
        nodes.update((ra, rb))

    comps = {}
    for n in nodes:
        comps.setdefault(_find(n), []).append(n)

    written = set()
    for members in comps.values():
        ## The anchor: a row another probe already wrote if there is one --
        ## rewriting it would undo that probe -- otherwise the least drifted.
        anchor = min(members, key=lambda n: (0 if n in pinned else 1,
                                             float(drift[n]), n))
        seen, stack = {anchor}, [anchor]
        while stack:
            u = stack.pop()
            for v, s in adj.get(u, ()):
                if v in seen:
                    continue
                seen.add(v)
                if v in pinned:
                    ## A second already-written row in the same component.
                    ## Writing it would undo the probe that put it there, so
                    ## this constraint is dropped and the walk stops here --
                    ## the same trade the cycle rule above makes.
                    continue
                out[v] = out[u] + s
                written.add(v)
                stack.append(v)
    return written


## ----------------------------------------------------------------------
## Circuit-level resolution (roadmap 15, spike 2026-08-25)
## ----------------------------------------------------------------------
##
## `SubCircuit.limit` lets every element write its own nodes back, in
## dictionary order, so two elements that want DIFFERENT values on one
## shared node are resolved by whoever writes last.  This is the adapter
## that lets `SubCircuit.limit` instead COLLECT each element's probe
## targets and run ONE `device_writeback` over the whole circuit.
##
## `CIRCUIT_LEVEL` is the switch; it exists so the intervention can be
## measured against its absence in one process.  See
## `test_device_limiter.py` section 5 for what the measurement said.
CIRCUIT_LEVEL = False


def element_targets(element, subx, subx0, epar, nodemap):
    """Re-derive an element's probe targets from its own `limit()`.

    Returns `(targets, limited)`: `targets` is `(ra, rb, vorig, vlim)` per
    declared probe with rows already mapped into the CIRCUIT's `x`, or
    `None` when the element has no `limit_spec` (a hand-written limiter
    such as `Diode`, which keeps the per-element write-back); `limited` is
    the element's own limited sub-vector, for that fallback.

    `vlim` is read as `out[ra] - out[rb]` on the element's limited output,
    so a probe whose constraint the device itself honoured carries exactly
    its own law's value, and a probe the device let FOLLOW another
    probe's write carries the value it ended with -- the device's whole
    decision, re-stated as constraints.  An element whose output is
    bit-identical to its input contributes nothing, which keeps the
    "no probe bit, nothing written" property at the circuit level.

    This is the duck-typed spike form.  The real integration would have
    the generated `limit()` in `hdl.py` return its `targets` list directly
    (it already builds one for grouped probes) instead of this adapter
    reading them back off the written vector.
    """
    info = getattr(element, '_hdl_info', None)
    spec = info.get('limit_spec') if isinstance(info, dict) else None
    limited = element.limit(subx, subx0, epar)
    if not spec or limited is None:
        return None, limited
    out = np.asarray(limited, dtype=float)
    xin = np.asarray(subx, dtype=float)
    if np.array_equal(out, xin):
        return [], limited
    targets = []
    for (ra, rb), _kind, _move, _pfs in spec:
        vorig = float(xin[ra] - xin[rb])
        vlim = float(out[ra] - out[rb])
        targets.append((int(nodemap[ra]), int(nodemap[rb]), vorig, vlim))
    return targets, limited
