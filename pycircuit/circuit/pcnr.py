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

**The unit is a DEVICE'S PROBE VECTOR** (vector PCNR, roadmap sec. 15, Stage 1,
2026-08-25).  A device with ``m`` limited quantities over ``t`` rows owns a
block ``v`` of ``m`` unknowns and answers with ``pcnr_i(v) -> R^t`` and
``pcnr_didv(v) -> R^{t x m}`` -- a BJT's transport current ``f(vbe, vbc)`` is
the case a scalar-per-branch protocol could not express.  Two device
protocols are accepted:

* **vector** -- the class declares ``pcnr_probes = (((ra, rb), kind), ...)``
  (indices into its own rows, a limiter kind per probe) and
  ``pcnr_i(v, params, epar, toolkit)``, ``pcnr_didv(...)`` and
  ``pcnr_limit(v_new, v_old, params, epar, toolkit, x_old_sub)`` over the
  whole block;
* **scalar** (the original one, unchanged) -- ``pcnr_junctions = ((a, c),
  ...)`` and per-junction ``pcnr_i(v, ..., jn=j) -> [i, -i]``.  This is the
  ``m = 1`` special case, and :class:`PcnrDevice` adapts it: every hand-
  written ``Diode``, every auto-detected DSL diode and the two-junction
  test devices stamp EXACTLY as before, entry for entry, which
  ``test_pcnr_vector.py`` asserts against a copy of the old formulas.
"""
import numpy as np

from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit._limiting import _pnjlim


_MISSING = object()


def _zero_vector(t):
    def _i(x, epar=defaultepar, params_tree=None):
        return np.zeros(t)
    return _i


def _linear_vector(gl):
    """``i(x) = gl @ x`` -- the remainder's own current.

    `x` is already the element's LOCAL sub-vector: `circuit.py` slices
    `subx = x[elementnodemap[instance]]` before calling `element.i`.
    Re-indexing by the device's global `rows` here would read the wrong
    entries, and would do it silently -- the shapes agree.
    """
    def _i(x, epar=defaultepar, params_tree=None):
        return gl.dot(np.asarray(x, dtype=float))
    return _i


def _linear_matrix(gl):
    """``G(x) = gl`` -- constant, so the node voltages are not read."""
    def _G(x, epar=defaultepar, params_tree=None):
        return gl
    return _G


def _zero_matrix(t):
    def _G(x, epar=defaultepar, params_tree=None):
        return np.zeros((t, t))
    return _G


class PcnrDevice(object):
    """One participating device: its rows, its probes, and how to ask it.

    ``rows`` are the device's rows in circuit coordinates (terminals, then
    internal nodes -- ``elementnodemap``'s order, the one ``element.i(sub)``
    is evaluated over); ``probes`` are ``(ra, rb)`` pairs in the same
    coordinates; ``local`` the same pairs as indices into ``rows``;
    ``kinds`` the limiter kind per probe.  ``off`` is the device's offset
    into the flat ``v_lim``/``g_lim`` vectors, assigned by
    :func:`pcnr_devices`.  ``scalar`` says which protocol the element
    speaks; the adapter methods below hide the difference from every
    consumer.
    """
    __slots__ = ('instance', 'element', 'rows', 'probes', 'local', 'kinds',
                 'off', 'scalar', 'redundant')

    def __init__(self, instance, element, rows, local, kinds, scalar,
                 off=0, redundant=()):
        self.instance = instance
        self.element = element
        self.rows = None if rows is None else [int(r) for r in rows]
        self.local = [(int(a), int(b)) for a, b in local]
        self.probes = None if rows is None else \
            [(self.rows[a], self.rows[b]) for a, b in self.local]
        self.kinds = tuple(kinds)
        self.scalar = bool(scalar)
        self.off = int(off)
        ## Probes that closed a cycle in the declaration (Stage 2): no
        ## unknown of their own -- ``((la, lb), kind, coeffs)`` with
        ## ``coeffs`` the signed combination of this device's unknowns that
        ## is their branch voltage.  Their LAW is applied by the device's
        ## `pcnr_limit`; here they matter as gmin targets when they are
        ## junctions.
        self.redundant = tuple(((int(a), int(b)), k, tuple(int(c_) for c_
                                                           in c))
                               for (a, b), k, c in redundant)

    @property
    def m(self):
        return len(self.local)

    @property
    def t(self):
        return len(self.rows)

    def params(self):
        """The device's parameter dict, cached on the element.

        Rebuilding this per junction per Newton iteration was 60 dict
        comprehensions over `getattr` on a 60-device circuit, every
        iteration.  A generated element's `$param_given` flags ride along
        under ``'$given:<name>'`` so a limiter parameter that reads one
        (SPICE's `RBM defaults to RB`) can be evaluated without `self`.
        """
        element = self.element
        params = getattr(element, '_pcnr_params', None)
        if params is None:
            params = {p.name: getattr(element.iparv, p.name)
                     for p in element.instparams}
            ipar = getattr(element, 'ipar', None)
            for nm in getattr(element, '_hdl_given_names', ()) or ():
                params['$given:' + nm] = (1.0 if ipar is not None
                                          and ipar.is_given(nm) else 0.0)
            element._pcnr_params = params
        return params

    ## -- the adapter -----------------------------------------------------
    def stamp(self, v, params, epar, g_mna):
        """Add the device's currents at its OWN ``v`` into ``g_mna`` and
        return its ``t x m`` derivative block.

        The scalar protocol is stamped per junction, in declaration order,
        exactly as the original per-junction loop did -- not accumulated
        into a block first -- so that a device whose junctions share a
        row (the base of `TwoJunction`) adds the same numbers in the same
        order and the result is bit-identical.
        """
        element = self.element
        toolkit = element.toolkit
        block = np.zeros((self.t, self.m))
        if self.scalar:
            for j, (la, lb) in enumerate(self.local):
                vj = float(v[j])
                i_terms = element.pcnr_i(vj, params, epar, toolkit, jn=j)
                di_terms = element.pcnr_didv(vj, params, epar, toolkit,
                                             jn=j)
                g_mna[self.rows[la]] += float(i_terms[0])
                g_mna[self.rows[lb]] += float(i_terms[1])
                block[la, j] += float(di_terms[0])
                block[lb, j] += float(di_terms[1])
            return block
        vv = np.asarray(v, dtype=float)
        cur = np.asarray(element.pcnr_i(vv, params, epar, toolkit),
                         dtype=float).reshape(self.t)
        ## `add.at`, not fancy-index `+=`: a diode-connected transistor
        ## has two terminals on ONE node, and `g[rows] += cur` would keep
        ## only the last of the duplicates.
        np.add.at(g_mna, self.rows, cur)
        block[:, :] = np.asarray(element.pcnr_didv(vv, params, epar,
                                                   toolkit),
                                 dtype=float).reshape(self.t, self.m)
        return block

    def limit(self, v_new, v_old, params, epar, x_old_sub):
        """The CORRECT phase for this device's block."""
        element = self.element
        toolkit = element.toolkit
        limiter = getattr(element, 'pcnr_limit', None)
        out = np.array(v_new, dtype=float, copy=True)
        if not self.scalar:
            if limiter is None:                          # pragma: no cover
                raise TypeError('%r declares pcnr_probes but no pcnr_limit'
                                % (self.instance,))
            return np.asarray(limiter(out, np.asarray(v_old, dtype=float),
                                      params, epar, toolkit, x_old_sub),
                              dtype=float).reshape(self.m)
        for j in range(self.m):
            if limiter is not None:
                out[j] = float(limiter(float(v_new[j]), float(v_old[j]),
                                       params, epar, toolkit, jn=j))
                continue
            ## A device that declares no `pcnr_limit` gets SPICE's `pnjlim`
            ## with its `IS`, which is what every junction device wants and
            ## what the hand-written `Diode` has always used.
            VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
            IS = getattr(element.iparv, 'IS', 0.0)
            out[j] = _pnjlim(float(v_new[j]), float(v_old[j]), VT, IS,
                             toolkit)
        return out


def _device_of(instance, element, rows):
    """The record for one element, or None when it does not participate."""
    probes = getattr(element, 'pcnr_probes', None)
    if probes:
        local = [tuple(p) for p, _k in probes]
        kinds = [k for _p, k in probes]
        return PcnrDevice(instance, element, rows, local, kinds, scalar=False,
                          redundant=getattr(element, 'pcnr_redundant', ()))
    pairs = getattr(element, 'pcnr_junctions', ())
    if pairs:
        return PcnrDevice(instance, element, rows, [tuple(p) for p in pairs],
                          ['pnj'] * len(pairs), scalar=True)
    return None


def pcnr_devices(circuit):
    """Every PCNR-participating device, with its offset into ``v_lim``."""
    found = []
    elements = getattr(circuit, 'elements', None)
    if not elements:
        return found
    nodemap = circuit.elementnodemap
    off = 0
    for instance, element in elements.items():
        dev = _device_of(instance, element, nodemap[instance])
        if dev is None:
            continue
        dev.off = off
        off += dev.m
        found.append(dev)
    return found


def pcnr_junction_pairs(circuit):
    """The p-n junctions, as ``(instance, device, row_anode, row_cathode)``.

    This is the GMIN-TARGET view: `dcanalysis._jrows` and
    `jaxtransient._gmin_junction_rows` put a conductance across each pair
    on the ordinary ``pcnr=False`` path too.  Only ``'pnj'`` probes are
    listed -- a gmin across a FET's ``vgs`` would be a gate leak -- so a
    device's ``fet``/``vds`` probes never reach a ladder.  The order is the
    flat ``v_lim`` order, which is what lets a consumer that still thinks
    in pairs (`fang_timestep`'s ``v_lim`` seed) stay correct while every
    probe is a junction.

    A device's REDUNDANT junctions (a MOSFET's `(b,d)`, which closed the
    cycle and has no unknown of its own) are listed after its own
    unknowns' pairs: a gmin belongs across both bulk junctions, as SPICE
    puts it, and which of the two the declaration happened to name last
    is not a physical distinction.  They carry no ``v_lim`` slot, so a
    pair-shaped ``v_lim`` consumer is correct only while every probe of
    every device is a tree junction -- which is the case it already
    required (a device with any ``fet``/``vds`` probe was never
    pair-shaped).
    """
    found = []
    for dev in pcnr_devices(circuit):
        for (ra, rb), kind in zip(dev.probes, dev.kinds):
            if kind == 'pnj':
                found.append((dev.instance, dev.element, ra, rb))
        for (la, lb), kind, _c in dev.redundant:
            if kind == 'pnj':
                found.append((dev.instance, dev.element, dev.rows[la],
                              dev.rows[lb]))
    return found


## The original name, kept for its consumers: it has always returned the
## pair list, and `jaxtransient` (Stage 3, not touched here) unpacks it as
## one.
pcnr_junctions = pcnr_junction_pairs


def _as_devices(circuit, junctions):
    """Device records from either shape a caller may hold.

    The tests hand `augmented_system` and `refine` the pair list
    `pcnr_junctions` returns; they get the same records `pcnr_devices`
    builds, grouped by instance in order.  A device's OWN declaration
    decides its protocol and its probes, so a pair list is only used to
    say which instances take part.
    """
    if not junctions:
        return []
    if all(isinstance(j, PcnrDevice) for j in junctions):
        return list(junctions)
    nodemap = circuit.elementnodemap if circuit is not None else None
    out, seen, off = [], {}, 0
    for instance, element, _ra, _rb in junctions:
        if instance in seen:
            continue
        rows = nodemap[instance] if nodemap is not None else None
        dev = _device_of(instance, element, rows)
        if dev is None:                                  # pragma: no cover
            raise ValueError('%r is in the junction list but declares no '
                             'PCNR protocol' % (instance,))
        dev.off = off
        off += dev.m
        seen[instance] = dev
        out.append(dev)
    return out


def flat_probes(junctions):
    """``(ra, rb)`` per flat limited unknown, in ``v_lim`` order."""
    if junctions and all(isinstance(j, PcnrDevice) for j in junctions):
        return [p for dev in junctions for p in dev.probes]
    return [(int(ra), int(rb)) for _i, _e, ra, rb in junctions]


def v_lim_init(junctions, x, epar=defaultepar, limit=False):
    """Each device's unknowns seeded from the branch voltages they stand
    for -- the paper's "independent initialization by different devices" --
    and then LIMITED by the device's own law, exactly as :func:`refine`
    limits every later iterate.

    **The seed used to be the raw branch voltage, and that was the whole
    of the "unlimited-node failure" on junction devices.** From a uniform
    20 V start a BJT mirror seeded `vbe = 20 V`, so the very first
    `augmented_system` evaluated the transport current at `exp(20/vt)`:
    `max|g_mna| = 1.6e92`, `cond(schur) = 4.6e94`, and the first step --
    computed from that Jacobian, before any limiting could act -- proposed
    4.5e45 V.  No ladder around the solve could help, because every rung
    began by building the same Jacobian at the same unlimited seed
    (roadmap sec. 15).

    Limiting the seed costs nothing where the start is already sane: at a
    zero start the raw seed is zero and every law returns it unchanged, so
    the vector is bit-identical.  Measured on the mirror: uniform 0 V and
    -5 V starts unchanged at 10 iterations, +5 V **165 -> 9**, +20 V
    **LinAlgError -> 8**, both instance orders, same answer.

    `x_old` is the start itself: a limiter parameter that follows the bias
    (SPICE's `von`) has no earlier iterate to read at initialisation.

    ⚠ **`limit` is off by default, and the distinction is the point: clamp
    an ARBITRARY start, never a converged one.** :func:`solve_dc` passes
    `limit=True` because its `x0` comes from the caller and may be
    anything. The transient's two call sites do NOT, because they seed
    from the previous time point's *accepted solution* -- a real operating
    point, and the best information available. Clamping that would discard
    it to no purpose, and it is measurable: doing so moved the rectifier's
    waveform by 1.7e-8 V and its time points by 6.5e-12 s against the
    classic-limiting path they are required to track.

    The clamp is inert on any sane branch voltage in any case -- seeds up
    to ~0.8 V come back unchanged on the Gummel-Poon card, and only 1 V
    and beyond are pulled in (1 V -> 0.094, 5 V -> 0.136, 20 V -> 0.172).
    """
    raw = np.array([float(x[ra] - x[rb]) for ra, rb in flat_probes(junctions)])
    if not limit:
        return raw
    return refine(junctions, np.zeros_like(raw), raw, epar, x_old=x)


def augmented_system(circuit, x, v_lim, junctions, epar=defaultepar,
                     u_extra=None, dense_blocks=True, J_extra=None):
    """Residual and Jacobian blocks of the coupled ``[x_MNA ; x_lim]`` system.

    Returns ``(g_mna, g_lim, J_mm, J_ml, J_lm, didv)``.  ``J_ll`` is not
    returned: it is the identity by construction, and materialising it would
    invite someone to use it.  ``didv`` is one ``(rows, probes, block)`` per
    device -- the sparse form of ``J_ml`` that :func:`schur_reduce` and
    :func:`predict` consume.
    """
    n = circuit.n
    devices = _as_devices(circuit, junctions)
    k = sum(dev.m for dev in devices)

    ## Ordinary assembly WITHOUT the participating devices: each is about to
    ## be re-stamped at its own `v_lim` instead of at the node voltages, so
    ## its ordinary `i`/`G` are replaced by zeros for the duration of the
    ## two assembly calls (an instance attribute shadowing the class method;
    ## restored in `finally`).  This is Design A of `doc/pcnr_native_design.md`
    ## -- the skip set -- done from inside the layer.
    ##
    ## It used to ASSEMBLE EVERYTHING AND SUBTRACT `element.i(sub)` again,
    ## and that is exact only in arithmetic.  A participant's current at
    ## the NODE voltages is unbounded -- PCNR limits `v_lim`, not the nodes
    ## -- and measured on a two-BJT mirror from a zero start (2026-08-25)
    ## the iterate visits vbe = 5.2 V at the nodes while `v_lim` sits at
    ## 0.8 V: `expl` keeps the current finite at ~1e72, the cancellation
    ## leaves noise of ~1e56 in the base row, and WHICH noise depends on
    ## the assembly's summation order, i.e. on instance order.  One order
    ## converged in 11 iterations and the other in 179, from identical
    ## states -- order dependence manufactured by the layer that exists to
    ## remove it.  `inf - inf = nan` (`test_limexp_is_what_makes_a_pcnr_
    ## participant_robust`) was the same defect at its limit.  Excluding
    ## the device from the sum has no such term to cancel.
    ##
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
    ## only consider algebraic equations".  (A DIFFUSION charge is the
    ## resistive current times a transit time, so on a BJT it is as wild
    ## at the node voltages as the current was; that is the price of this
    ## trade on the transient path, and it is not paid at DC.)
    ##
    ## Proven, not argued: `test_pcnr_charge.py` finite-differences the
    ## augmented residual to confirm `J_eff == df_eff/dx` with a
    ## charge-storing participant, and compares pcnr on/off transients.
    saved = []
    try:
        for dev in devices:
            el, t = dev.element, dev.t
            saved.append((el, el.__dict__.get('i', _MISSING),
                          el.__dict__.get('G', _MISSING)))
            ## THE AFFINE REMAINDER (roadmap sec. 40).  A participating
            ## device is normally zeroed out of the ordinary assembly
            ## and re-stamped at `v_lim`.  When part of its current is a
            ## CONSTANT CONDUCTANCE -- a series resistance to an
            ## internal node -- that part reaches the solution through
            ## node voltages, not through any probe, and PCNR's rule
            ## refuses the whole device for it.  Newton solves a linear
            ## term exactly and limiting it would mean nothing, so the
            ## remainder is shadowed IN rather than zeroed out, and the
            ## ordinary assembly stamps it exactly as on `pcnr=False`.
            lin = getattr(el, 'pcnr_lin_G', None)
            if lin is None:
                el.i = _zero_vector(t)
                el.G = _zero_matrix(t)
            else:
                gl = np.asarray(lin(dev.params, epar, el.toolkit),
                                dtype=float)
                el.i = _linear_vector(gl)
                el.G = _linear_matrix(gl)
        g_mna = np.array(circuit.i(x, epar), dtype=float)
        J_mm = np.array(circuit.G(x, epar), dtype=float)
    finally:
        for el, old_i, old_G in saved:
            for name, old in (('i', old_i), ('G', old_G)):
                if old is _MISSING:
                    el.__dict__.pop(name, None)
                else:
                    el.__dict__[name] = old
    if u_extra is not None:
        g_mna = g_mna + np.asarray(u_extra, dtype=float)
    if J_extra is not None:
        ## The transient's companion terms: `f = i(x) + iq + u`, `J = G + Geq`.
        ## They are passed in rather than recomputed because `solve_timestep`
        ## already has them and `get_diff` is not free.
        J_mm = J_mm + np.asarray(J_extra, dtype=float)

    ## `dense_blocks=False` skips the (n,k) and (k,n) allocations entirely.
    ## `predict`'s sparse path never reads them -- each column of J_ml has t
    ## nonzeros (the device's rows) and each row of J_lm has two, so the
    ## rank-one form carries the same information in `didv`.  Allocating and
    ## filling them anyway was 2 x n x k doubles per Newton iteration for
    ## nothing.
    J_ml = np.zeros((n, k)) if dense_blocks else None
    J_lm = np.zeros((k, n)) if dense_blocks else None
    g_lim = np.zeros(k)
    didv_list = []

    for dev in devices:
        off, m = dev.off, dev.m
        v = np.asarray(v_lim[off:off + m], dtype=float)
        params = dev.params()
        ## The device's current now enters through its OWN unknowns, so its
        ## contribution to J_MNA/MNA is zero and all of it lands in J_MNA/lim.
        block = dev.stamp(v, params, epar, g_mna)
        rows = np.asarray(dev.rows, dtype=int)
        for j, (ra, rb) in enumerate(dev.probes):
            ## g_lim = v - (e_a - e_b)
            g_lim[off + j] = float(v[j]) - (float(x[ra]) - float(x[rb]))
            if dense_blocks:
                np.add.at(J_ml[:, off + j], rows, block[:, j])
                ## ACCUMULATED, not assigned: a diode-connected transistor's
                ## base-collector probe spans ONE row (ra == rb), and its
                ## incidence row is -1 + 1 = 0 -- which is what the sparse
                ## path computes and what an assignment would have left
                ## at +1.
                J_lm[off + j, ra] += -1.0
                J_lm[off + j, rb] += +1.0
        didv_list.append((dev.rows, dev.probes, block))

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
    update; with ``didv`` that is exploited directly.  ``J_ml J_lm`` is a sum
    over UNKNOWNS: column ``off+j`` of ``J_ml`` is the device's ``block[:, j]``
    on its ``t`` rows, row ``off+j`` of ``J_lm`` is ``(-1 at ra_j, +1 at rb_j)``,
    so each probe is one rank-one term touching ``2t`` entries -- ``O(sum m t)``
    rather than the ``O(n^2 k)`` of a dense ``J_ml @ J_lm``.  A scalar diode is
    the ``t = 2, m = 1`` case: the same four entries as before, in the same
    order.  (``junctions`` is accepted and ignored: ``didv`` carries the rows.)
    """
    if didv is not None:
        schur = np.array(J_mm, copy=True)
        rhs_corr = np.zeros(len(g_mna))
        off = 0
        for rows, probes, block in didv:
            rows = np.asarray(rows, dtype=int)
            for j, (ra, rb) in enumerate(probes):
                col = block[:, j]
                np.add.at(schur[:, ra], rows, col)
                np.subtract.at(schur[:, rb], rows, col)
                np.add.at(rhs_corr, rows, col * g_lim[off + j])
            off += len(probes)
        return g_mna - rhs_corr, schur
    return g_mna - J_ml @ g_lim, J_mm - J_ml @ J_lm


def dx_lim_of(junctions, g_lim, dx_mna):
    """``dx_lim = -(g_lim + J_lm dx_MNA)``, per limited unknown.

    Row ``k`` of ``J_lm`` is ``(-1 at ra, +1 at rb)``, so this is the same
    two-entry formula whatever the device's ``m``.  One definition, shared
    by :func:`predict` and the transient's coupled path, which used to copy
    it by hand.
    """
    probes = flat_probes(junctions)
    dx_lim = np.empty(len(probes))
    for k, (ra, rb) in enumerate(probes):
        dx_lim[k] = -(g_lim[k] + (-dx_mna[ra] + dx_mna[rb]))
    return dx_lim


def predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode, junctions=None,
            didv=None, gmin=0.0):
    """One Newton step on the coupled system, by Schur complement.

        dx_MNA = -(J_mm - J_ml J_lm)^-1 (g_mna - J_ml g_lim)
        dx_lim = -(g_lim + J_lm dx_MNA)

    The matrix factorised is the size of the ORIGINAL MNA system: the ``lim/lim``
    block being the identity is what collapses the border into a rank-k update.

    **That collapse has to be exploited, not merely stated.** Forming
    ``J_ml @ J_lm`` as a dense ``(n,k)·(k,n)`` product costs ``O(n^2 k)`` and
    measured **+62% per iteration** against classic limiting on 60 diodes -- a
    gate-13-4 failure that was entirely the implementation's.  Pass
    ``junctions``/``didv`` to take the sparse path.
    """
    f_eff, schur = schur_reduce(g_mna, g_lim, J_mm, J_ml, J_lm, junctions, didv)
    rhs = -f_eff

    ## THE DAMPING, and note WHERE it goes: on the Jacobian only, never on
    ## the residual.  `f_eff` is what convergence is judged on, so the
    ## converged point solves the untouched equations and the anchor cannot
    ## bias the answer -- verified against `DC()` to 1e-11 on all three
    ## circuits that needed it.  What it changes is the DIRECTION of the
    ## step, which is the whole difficulty: on a wild start with every
    ## device cut off, the Schur complement is RANK DEFICIENT (an EKV
    ## differential pair at a uniform 20 V start measured rank 8 of 9, a
    ## tail node with 2e-19 S to anywhere), and a rank-deficient solve does
    ## not produce a large step, it produces a meaningless one -- 3e48 V.
    ##
    ## ⚠ Two cheaper things were measured first and BOTH are wrong:
    ##   * truncating the step (`dx *= cap/max|dx|`) rescues the diff pair
    ##     at every cap from 10 V to 1e6 V and BREAKS the Level-1 cascode
    ##     at every one of them.  A huge step can be better than a
    ##     truncated one, so a magnitude cap is not a safe damper.
    ##   * backtracking on the residual norm improves the cascode (9 -> 5)
    ##     and does NOT fix the diff pair, because the step it accepts is
    ##     still one whose NEXT Jacobian is singular.  A line search
    ##     controls the step; it cannot condition the matrix.
    ## Only regularising the matrix addresses a rank deficiency.
    if gmin:
        schur = np.array(schur, dtype=float, copy=True)
        d = np.arange(schur.shape[0])
        schur[d, d] += gmin

    keep = [i for i in range(len(g_mna)) if i != irefnode]
    dx_r = np.linalg.solve(schur[np.ix_(keep, keep)], rhs[keep])
    dx_mna = np.zeros(len(g_mna))
    dx_mna[keep] = dx_r

    if junctions is not None:
        dx_lim = dx_lim_of(junctions, g_lim, dx_mna)
    else:
        dx_lim = -(g_lim + J_lm @ dx_mna)
    return dx_mna, dx_lim


def limit_block(probes, coeffs, v_new, v_old, laws):
    """A device's laws applied over its tree of unknowns (Stage 2).

    ``probes`` are ``(a, b)`` row pairs for EVERY probe whose law is
    applied -- the ``m`` tree probes first, then any redundant ones --
    and ``coeffs[p]`` is the signed combination of the ``m`` unknowns that
    is probe ``p``'s branch voltage (an identity row for a tree probe).
    ``laws[p](vnew, vold)`` is its limiter.  Returns the ``m`` limited
    unknowns.

    **Without a redundant probe this is the per-probe loop, exactly**:
    each unknown is its own law's output, bit for bit, which is what
    keeps every Stage-1 participant's numbers where they were.

    **With one, the laws over-determine the unknowns** -- SPICE's MOSFET
    set has `(b,d) = (b,s) - (d,s)` and three limited values need not
    sum to zero -- and the rule is `limit_together`'s write-back rule
    (`_limiting.device_writeback`), restated for a tree of branch
    voltages instead of a forest over nodes: take the probes in order of
    DECREASING correction, keep each that does not close a cycle, and
    solve the unknowns from the kept set.  The one dropped is the one
    asking for the least, and the rule is a function of the data, so
    declaration order cannot matter.  There is no anchor to choose here
    -- the unknowns ARE branch voltages -- which is the part of the
    node write-back that PCNR made unnecessary.

    If no law bit, ``v_new`` is returned untouched: the strict no-op
    near the solution that lets "did limiting fire?" stay a signal.  If
    the kept set is the tree itself, the limited values are assigned
    directly (no potentials are formed, so nothing is rounded).
    """
    m = len(v_new)
    vn = np.asarray(v_new, dtype=float)
    vo = np.asarray(v_old, dtype=float)
    C = np.asarray(coeffs, dtype=float).reshape(len(probes), m)
    raw = C @ vn
    old = C @ vo
    lim = np.array([float(law(float(r), float(o)))
                    for law, r, o in zip(laws, raw, old)])
    bit = [p for p in range(len(probes)) if lim[p] != raw[p]]
    if not bit:
        return np.array(vn, copy=True)
    if len(probes) == m:
        return lim
    ## Kruskal on |correction|, ties by row pair -- the data decides.
    parent = {}

    def _find(n):
        parent.setdefault(n, n)
        while parent[n] != n:
            parent[n] = parent[parent[n]]
            n = parent[n]
        return n

    kept = []
    for p in sorted(range(len(probes)),
                    key=lambda q: (-abs(lim[q] - raw[q]),
                                   probes[q][0], probes[q][1])):
        a, b = probes[p]
        ka, kb = _find(a), _find(b)
        if ka == kb:
            continue
        parent[ka] = kb
        kept.append(p)
    if sorted(kept) == list(range(m)):
        return lim[:m]
    ## Potentials over the kept tree, then each unknown as a difference.
    adj = {}
    for p in kept:
        a, b = probes[p]
        adj.setdefault(a, []).append((b, -lim[p]))
        adj.setdefault(b, []).append((a, +lim[p]))
    pot = {}
    for start in adj:
        if start in pot:
            continue
        pot[start] = 0.0
        stack = [start]
        while stack:
            u = stack.pop()
            for w, sgn in adj[u]:
                if w not in pot:
                    pot[w] = pot[u] + sgn
                    stack.append(w)
    return np.array([pot[a] - pot[b] for a, b in probes[:m]])


def refine(junctions, v_old, v_new, epar=defaultepar, x_old=None):
    """The CORRECT phase: each device limits only the variables it owns.

    This is where PCNR differs from limiting in the way that matters. Nothing is
    shared, so applying one device's limiter cannot disturb another's -- the
    clash section 2 of the paper describes cannot arise, by construction rather
    than by ordering.

    A device is handed its WHOLE block at once, plus its own slice of the
    last accepted iterate ``x_old`` (a limiter parameter that follows the
    bias -- SPICE's `von` -- is evaluated there).  There is no sequential
    reading between probes: that was forced on the ordinary path by node
    write-back, and PCNR writes no node.
    """
    devices = _as_devices(None, junctions)
    out = np.array(v_new, dtype=float, copy=True)
    for dev in devices:
        off, m = dev.off, dev.m
        params = dev.params()
        x_sub = None
        if x_old is not None and dev.rows is not None:
            x_sub = np.asarray(x_old, dtype=float)[np.asarray(dev.rows,
                                                             dtype=int)]
        ## THE DEVICE OWNS THE LAW FOR ITS OWN QUANTITY.  PCNR supplies the
        ## architecture -- one unknown per limited quantity, nothing shared,
        ## so no device's limiter can disturb another's -- and the device
        ## supplies the limiter, which is the paper's modularity claim.
        out[off:off + m] = dev.limit(out[off:off + m], v_old[off:off + m],
                                     params, epar, x_sub)
    return out


def lim_converged(g_lim, v_new, reltol, abstol):
    """Per-component: ``|g_lim,j| < reltol*max(|v_j|, 1) + abstol``.

    Per COMPONENT, not against ``max|v_new|`` over the whole vector: with
    vector devices a 40 V ``vds`` in the same vector would loosen a
    ``vbe`` component forty-fold.  Where every device is a scalar diode
    on one branch the two criteria coincide.
    """
    g = np.abs(np.asarray(g_lim, dtype=float))
    scale = np.maximum(np.abs(np.asarray(v_new, dtype=float)), 1.0)
    return bool(np.all(g < reltol * scale + abstol))


def solve_dc(circuit, refnode, x0=None, epar=defaultepar, maxiter=200,
             reltol=1e-6, abstol=1e-12, iabstol=1e-12, gmin=0.0):
    """DC operating point by PCNR.  Returns ``(x, v_lim, iterations)``.

    The flow is Figure 2's: initialise, then repeat predict-then-correct until
    converged.  There is no limiting of the MNA vector anywhere -- the only thing
    limited is each device's own unknown, in :func:`refine`.

    **Convergence is `StandardNewton`'s test, per component** (Stage 2,
    2026-08-25): every ``|dx_i| < reltol*max(|x_i|, |x_i'|) + abstol``,
    every ``|f_i| < reltol*(|J||x| + |f|)_i + iabstol``, and every limited
    unknown consistent (:func:`lim_converged`).  It used to be
    ``max|dx| < reltol*max|x_new|`` over the WHOLE vector -- and the
    vector holds the sources' branch currents.  On a FET cascode from a
    40 V supply the first iterate put the supply's current at 7.5e27 A,
    the tolerance on every node became 7e23 V, and the solve declared
    victory in 5 iterations at a point whose KCL residual was 7.5e27:
    a "converged" non-solution, 688 V from the answer.  Diodes never
    exposed it because a diode's current at a limited ``v_lim`` is
    bounded; a FET's ``vgs`` unknown is limited but the current at it
    still overflows the row.  The residual test is what separates a
    small step from a small residual -- at a 1e29 S point the Newton
    step is small BECAUSE the row is stiff, not because it is solved.

    **`gmin` is a damper on the step, not a term in the equations.** It
    anchors the Schur complement's diagonal inside :func:`predict` and
    never touches the residual, so the point this converges to solves the
    untouched system -- matched against `DC()` to 1e-11. ``gmin=0``
    restores the undamped behaviour, and is the default.

    ⚠⚠ **ITS MOTIVATING CIRCUIT IS GONE (2026-08-26).** This was built
    for the EKV differential pair from a uniform 20 V start, where the
    reduced matrix goes rank deficient (measured: rank 8 of 9, both
    channels cut off, 2e-19 S to the tail) and an unregularised solve
    returns not a large step but a meaningless one. That pair now
    converges undamped: `limit_delta` on EKV's bulk gave the probe a law
    for `v_lim_init`'s limited seed to apply (roadmap sec. 27, 34). A
    72-configuration probe afterwards -- four circuit families, starts
    from -40 V to +100 V, both instance orders -- found **nothing that
    `solve_dc` fails on**.

    So this damper currently rescues no circuit anyone can name, while
    still costing cascode grid point (2, 2, 2) when switched on. It is
    kept because "nothing I can find fails" is not "nothing fails", and
    because it remains the only tool for a genuinely rank-deficient
    start -- but it should not be reached for without a failing case in
    hand.

    ⚠⚠ **IT DEFAULTS TO OFF, because it is not free: it trades one
    circuit for another.** With `gmin = 1e-9` the EKV differential pair
    converges from a uniform 20 V start where it previously raised
    `LinAlgError` -- and cascode grid point `(2, 2, 2)`, which converges
    in 8 with no anchor, stops converging inside 200 iterations. That is
    true at **every value swept, 1e-12 through 1e-6**, so there is no
    window that buys the first without paying the second; the anchor
    perturbs that point into a limit cycle rather than conditioning it.

    The same shape defeated the two cheaper dampers (see :func:`predict`):
    each fixes one circuit and breaks another. `800 < total < 1000` over
    that grid with **nothing that converged stopping** is the property the
    suite pins, and it outranks rescuing a start the fallback already
    handles.

    It is also non-monotone where it does help: a Level-1 cascode from a
    uniform 20 V start fails at ``1e-14``/``1e-12``/``1e-11``, converges
    from ``1e-10`` up, **and converges at exactly zero**. A perturbation
    too small to condition the matrix is still large enough to move the
    trajectory. Treat any iteration count from a wild start as unstable.

    Turn it on deliberately, per solve, when a rank-deficient start is the
    known problem -- and re-measure, do not assume.
    """
    devices = pcnr_devices(circuit)
    if not devices:
        raise ValueError('no device in this circuit declares a PCNR junction')

    n = circuit.n
    irefnode = circuit.get_node_index(refnode)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float, copy=True)
    x[irefnode] = 0.0

    ## `limit=True`: `x0` is the caller's, and may be a wild start.
    v_lim = v_lim_init(devices, x, epar, limit=True)

    u = np.asarray(circuit.u(0.0, epar, analysis='dc'), dtype=float)

    for it in range(maxiter):
        g_mna, g_lim, J_mm, J_ml, J_lm, didv = augmented_system(
            circuit, x, v_lim, devices, epar, u_extra=u, dense_blocks=False)

        dx_mna, dx_lim = predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode,
                                 junctions=devices, didv=didv, gmin=gmin)

        x_new = x + dx_mna
        x_new[irefnode] = 0.0
        v_new = v_lim + dx_lim

        ## CORRECT: each device limits only what it owns.
        v_new = refine(devices, v_lim, v_new, epar, x_old=x)

        keep = np.arange(n) != irefnode
        conv_x = bool(np.all(np.abs(dx_mna)[keep]
                             < reltol * np.maximum(np.abs(x_new), np.abs(x))[keep]
                             + abstol))
        ## The residual the step was computed from, judged the way
        ## `StandardNewton` judges its own: against the currents that
        ## actually flow in the row.
        f_eff, J_eff = schur_reduce(g_mna, g_lim, J_mm, junctions=devices,
                                    didv=didv)
        i_scale = np.abs(J_eff) @ np.abs(x_new) + np.abs(f_eff)
        conv_f = bool(np.all(np.abs(f_eff)[keep]
                             < reltol * i_scale[keep] + iabstol))
        converged = (conv_x and conv_f
                     and lim_converged(g_lim, v_new, reltol, abstol))
        x, v_lim = x_new, v_new
        if converged:
            return x, v_lim, it + 1

    raise RuntimeError(
        'PCNR did not converge in %d iterations: max|g_lim| = %g, '
        'max|dx| = %g, max|f| = %g'
        % (maxiter, float(np.max(np.abs(g_lim))), float(np.max(np.abs(dx_mna))),
           float(np.max(np.abs(f_eff)))))
