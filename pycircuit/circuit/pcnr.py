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
from pycircuit.circuit.linearsolver import AutoSolver
from pycircuit.utilities.param import ParameterDict


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
                 'off', 'scalar', 'redundant', '_kmap')

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
        self._kmap = None

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

        Rebuilt where the parameter epoch moved (every parameter write moves
        it, `ParameterDict.notify`).  ⚠ Until 2026-10-14 it was kept for the
        element's life: a parameter written after the first PCNR solve was
        never seen by anything reading it (`pcnr_view`'s remainder stub, the
        limiter's law parameters).
        """
        element = self.element
        cached = getattr(element, '_pcnr_params', None)
        epoch = ParameterDict._epoch
        if cached is not None and cached[0] is epoch:
            return cached[1]
        params = {p.name: getattr(element.iparv, p.name)
                 for p in element.instparams}
        ipar = getattr(element, 'ipar', None)
        for nm in getattr(element, '_hdl_given_names', ()) or ():
            params['$given:' + nm] = (1.0 if ipar is not None
                                      and ipar.is_given(nm) else 0.0)
        element._pcnr_params = (epoch, params)
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
        ## THE KERNEL ROUTE (2026-10-12, `dc-pcnr-c`): the device's
        ## ORDINARY evaluation -- the compiled kernel where it has one --
        ## at a local point that realises `v`, not the SymPy `pcnr_i` /
        ## `pcnr_didv`: the same numbers to rounding on every gold class
        ## (7e-15), 13x faster, and the evaluation a C unit can call.
        cur, block[:, :] = self._kernel_stamp(np.asarray(v, dtype=float),
                                              params, epar)
        ## `add.at`, not fancy-index `+=`: a diode-connected transistor
        ## has two terminals on ONE node, and `g[rows] += cur` would keep
        ## only the last of the duplicates.
        np.add.at(g_mna, self.rows, cur)
        return block

    def _kernel_map(self, params, epar):
        """`(G_lin, walk, mirror, T)`: the affine remainder (or None), the
        probe tree's walk -- `(row, parent, probe, sign)` in visiting order,
        `pot[row] = pot[parent] + sign * v[probe]`, the root at 0 -- the
        rows it does not reach, each with the reached row it mirrors
        (tied through the remainder, so that branch carries exactly
        zero), and `T = d pot / d v` (entries -1, 0, +1).  Cached per
        temperature."""
        key = float(getattr(epar, 'T', 300.0))
        if self._kmap is not None and self._kmap[0] == key:
            return self._kmap[1]
        el = self.element
        lin = getattr(el, 'pcnr_lin_G', None)
        G_lin = (None if lin is None else
                 np.asarray(lin(params, epar, el.toolkit), dtype=float))
        t, m = self.t, self.m
        T = np.zeros((t, m))
        reached = [False] * t
        walk = []
        ## every component of the probes from its own root at 0 -- a
        ## thermal or optical port is a component apart from the junctions,
        ## and the currents read differences within a component only
        for _j0, (_la0, root) in enumerate(self.local):
            if reached[root] or reached[_la0]:
                continue
            reached[root] = True
            stack = [root]
            while stack:
                a = stack.pop()
                for j, (la, lb) in enumerate(self.local):
                    ## v_j = pot[la] - pot[lb]
                    for b, sign in ((la, +1.0), (lb, -1.0)):
                        other = lb if b == la else la
                        if other == a and not reached[b]:
                            reached[b] = True
                            T[b] = T[a]
                            T[b, j] += sign
                            walk.append((b, a, j, sign))
                            stack.append(b)
        mirror = []
        changed = True
        while changed:
            changed = False
            for r in range(t):
                if reached[r] or G_lin is None:
                    continue
                for c in range(t):
                    if reached[c] and G_lin[r, c] != 0.0:
                        reached[r] = True
                        T[r] = T[c]
                        mirror.append((r, c))
                        changed = True
                        break
        out = (G_lin, tuple(walk), tuple(mirror), T)
        self._kmap = (key, out)
        return out

    def _kernel_stamp(self, v, params, epar):
        """`(current, block)` on the device's rows at its own `v`: its
        ordinary `i` and `G` at the local point `pot` (`_kernel_map`), the
        affine remainder taken off (`G_lin pot`, `G - G_lin`), and
        `block[:, j] = sum_r G_nl[:, r] T[r, j]` -- every sum in a fixed
        order, the C unit's (`csrc/pcnrdc.c`)."""
        el = self.element
        G_lin, walk, mirror, T = self._kernel_map(params, epar)
        t, m = self.t, self.m
        pot = np.zeros(t)
        for row, parent, j, sign in walk:
            pot[row] = pot[parent] + sign * v[j]
        for row, src in mirror:
            pot[row] = pot[src]
        cur = np.array(el.i(pot, epar), dtype=float).reshape(t)
        G = np.array(el.G(pot, epar), dtype=float).reshape(t, t)
        if G_lin is not None:
            for r in range(t):
                acc = 0.0
                for c in range(t):
                    acc += G_lin[r, c] * pot[c]
                cur[r] -= acc
            G -= G_lin
        block = np.zeros((t, m))
        for j in range(m):
            for r in range(t):
                if T[r, j] != 0.0:
                    block[:, j] += G[:, r] * T[r, j]
        return cur, block

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


_REMAINDER_CLASSES = {}


def _remainder_class(nterm, branches):
    """The constant-stamp element a participant stands as in `pcnr_view`:
    `nterm` terminals (`p0`, `p1`, ...), the participant's `branches` as
    `(plus, minus)` terminal indices, `G` its stored matrix and `i` the
    stamp plan's `G x` (`i` is not overridden, so `pair_kind` says
    'cached')."""
    key = (nterm, tuple(branches))
    cls = _REMAINDER_CLASSES.get(key)
    if cls is None:
        from pycircuit.circuit.circuit import Branch, Circuit, Node

        def G(self, x, epar=defaultepar):
            return self._G
        cls = _REMAINDER_CLASSES[key] = type(
            f'_PcnrRemainder{nterm}', (Circuit,),
            {'terminals': tuple(f'p{k}' for k in range(nterm)),
             'branches': [Branch(Node(f'p{a}'), Node(f'p{b}')) for a, b in branches],
             '_constant_stamps': ('G',), 'instparams': [], 'G': G,
             '__doc__': "A participant's affine remainder (`pcnr_view`)."})
    return cls


def _remainder_of(circuit, dev, epar):
    """The stub standing for `dev` in `pcnr_view`, or None where it has
    neither a remainder nor a branch: its node rows as terminals, its
    branch rows re-declared between them, `G` the remainder (zeros where
    none) over all its rows."""
    nn = len(circuit.nodes)
    node_rows = [r for r in dev.rows if r < nn]
    lin = getattr(dev.element, 'pcnr_lin_G', None)
    branch_rows = [r for r in dev.rows if r >= nn]
    if lin is None and not branch_rows:
        return None, None
    index = {}
    for k, r in enumerate(node_rows):
        index.setdefault(circuit.nodes[r], k)
    brs = []
    for r in branch_rows:
        b = circuit.branches[r - nn]
        brs.append((index[b.plus], index[b.minus]))
    stub = _remainder_class(len(node_rows), brs)()
    stub._G = (np.zeros((dev.t, dev.t)) if lin is None else
               np.asarray(lin(dev.params(), epar, dev.element.toolkit), dtype=float))
    return stub, {f'p{k}': circuit.nodes[r] for k, r in enumerate(node_rows)}


def pcnr_view(circuit, devices, epar=defaultepar):
    """`circuit` without its participating `devices`: every other element
    in its order on the SAME nodes in the same order, and in each
    participant's place the constant stamp of its affine remainder
    (`pcnr_lin_G`) carrying its branches (`_remainder_of`) -- so the same
    rows, branches included.  Cached on the circuit per participants and
    temperature, and rebuilt where the stamp plan is (`_stamp_plan._plan_for`):
    the parameter epoch moved or the element node map is a new object.

    ⚠ The participants and temperature alone were the key until 2026-10-14,
    and a circuit edited between two solves kept its first view: an added
    resistor was not in the second solve (the old answer, status 'used'), an
    added unknown left PCNR's C core on a view one size short of its buffers
    (ASan: a heap-buffer-overflow in core.c), and a Gummel-Poon's `re`
    changed moved the device but not its remainder's stub (v(c) -0.682 V,
    a fresh circuit 4.230).  `test_pcnr_c`'s "between two solves" tests."""
    from pycircuit.circuit.circuit import SubCircuit
    key = (tuple(d.instance for d in devices), float(getattr(epar, 'T', 300.0)))
    cached = getattr(circuit, '_pcnr_view', None)
    if (cached is not None and cached[0] == key and cached[2] is circuit.elementnodemap
            and cached[3] is ParameterDict._epoch):
        return cached[1]
    part = {d.instance: d for d in devices}
    view = SubCircuit()
    view.toolkit = circuit.toolkit
    for node in circuit.nodes:
        view.append_node(node)
    for inst, el in circuit.elements.items():
        dev = part.get(inst)
        if dev is None:
            view.add_instance(inst, el, **circuit.term_node_map[inst])
            continue
        stub, nodes = _remainder_of(circuit, dev, epar)
        if stub is not None:
            view.add_instance(inst + '#pcnr_lin', stub, **nodes)
    if view.n != circuit.n:                              # pragma: no cover
        raise ValueError(f'pcnr_view: {view.n} unknowns, the circuit {circuit.n}')
    ## (the epoch read after the build: adding the circuit's elements to the
    ## view writes their parameters, which may move it)
    circuit._pcnr_view = (key, view, circuit.elementnodemap, ParameterDict._epoch)
    return view


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
    ## THE ORDINARY PART ON A VIEW (2026-10-12, `dc-pcnr-c` C1b): the
    ## circuit without its participants, each one's AFFINE REMAINDER (roadmap
    ## sec. 40: a series resistance to an internal node reaches the solution
    ## through node voltages, not through a probe, so it is stamped here
    ## exactly as on `pcnr=False`) in its place as a constant stamp -- the
    ## same nodes in the same order (`pcnr_view`).  It used to SHADOW each
    ## participant's `i` and `G` for the duration of the two calls; the view
    ## sums the same terms in the stamp plan's one order, which a C core
    ## built on it repeats bit for bit.
    view = pcnr_view(circuit, devices, epar)
    g_mna = np.array(view.i(x, epar), dtype=float)
    J_mm = np.array(view.G(x, epar), dtype=float)
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
            didv=None, gmin=0.0, linsolver=None, toolkit=None, held=()):
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

    ``held`` rows leave the solve with the reference's (their step 0): the
    nodes :func:`solve_dc` holds.
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

    out = {irefnode, *held}
    keep = [i for i in range(len(g_mna)) if i not in out]
    ## (`linsolver`: the DC's own -- dense below its floor, StructuralKLU on
    ## the circuit's pattern above -- which the C unit repeats; numpy's
    ## dense solve where none is given)
    if linsolver is None:
        dx_r = np.linalg.solve(schur[np.ix_(keep, keep)], rhs[keep])
    else:
        dx_r = np.asarray(linsolver.solve(schur[np.ix_(keep, keep)], rhs[keep],
                                          toolkit), dtype=float).ravel()
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
    ## (each probe's value as a loop in the unknowns' order -- the C unit's
    ## order; `C @ v` left the summation to BLAS)
    raw = np.zeros(len(probes))
    old = np.zeros(len(probes))
    for p in range(len(probes)):
        for j in range(m):
            if C[p, j] != 0.0:
                raw[p] += C[p, j] * vn[j]
                old[p] += C[p, j] * vo[j]
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


#: :func:`solve_dc`'s RuntimeError where an iterate is not finite (the
#: iteration's number); `_pcnr_c` raises the same
NONFINITE = 'PCNR iterate not finite at iteration %d: it cannot converge from there'


def solve_dc(circuit, refnode, x0=None, epar=defaultepar, maxiter=200,
             reltol=1e-6, abstol=1e-12, iabstol=1e-12, gmin=0.0,
             gshunt=0.0, source=1.0, v_lim0=None, linsolver=None, jgmin=0.0, held=None,
             anchor=None):
    """DC operating point by PCNR.  Returns ``(x, v_lim, iterations)``.

    ``gshunt``, ``jgmin``, ``source`` and ``v_lim0`` are
    :func:`solve_dc_continued`'s: a conductance from every node to ground
    IN the equations (a deformation, unlike ``gmin`` below), a conductance
    across every junction pair (:func:`pcnr_junction_pairs`) -- SPICE's
    gmin stepped, true gmin stepping -- the independent sources scaled,
    and the limited unknowns to start from.  At their defaults nothing is
    added and nothing is scaled: the plain solve is unchanged.

    ``anchor``: the gshunt conductance anchored there, ``gshunt (x - anchor)``
    on the node rows instead of ``gshunt x`` -- a pseudo-transient step
    (:func:`_ptc_ladder`).  Carried in the sources the loop is handed
    (``-gshunt anchor`` on those rows), so the C loop and the Python loop
    take it alike.

    ``held`` ({row: volts}) holds nodes at ``source`` times their volts --
    `DC`'s ``pin`` and ``nodeset``, eliminated as the ordinary chain's
    `_Held` eliminates them: not unknowns of the solve, their step 0 and
    their KCL rows out of the residual test (the hold supplies that current).

    Raises RuntimeError at the iteration limit, and at the first iterate
    that is not finite (:data:`NONFINITE`): `x + dx` keeps a non-finite
    unknown so and a NaN fails every test, so the solve has failed there --
    a rung's outcome is the same, found up to `maxiter` iterations sooner.

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
    ## The view first (cached): building it adds every element to a new
    ## circuit, a parameter notification that moves the epoch -- after the
    ## circuit's own stamp plan was built that plan went stale and a
    ## transient rebuilt it once after its operating point
    pcnr_view(circuit, devices, epar)

    n = circuit.n
    irefnode = circuit.get_node_index(refnode)
    x = np.zeros(n) if x0 is None else np.array(x0, dtype=float, copy=True)
    x[irefnode] = 0.0
    hrows = sorted(held) if held else ()
    for r in hrows:
        x[r] = source * held[r]
    if hrows and type(linsolver) is AutoSolver:
        ## (the held solve's own solver: its reduced matrix is not the
        ## circuit's structure -- `AutoSolver.held`, the rows numbered
        ## without the reference's, as the structure numbers them)
        linsolver = linsolver.held([r - (r > irefnode) for r in hrows])

    ## `limit=True`: `x0` is the caller's, and may be a wild start.  A
    ## rung of `solve_dc_continued` passes `v_lim0`, the previous rung's
    ## converged unknowns -- a real point, not clamped (`v_lim_init`).
    v_lim = (v_lim_init(devices, x, epar, limit=True) if v_lim0 is None
             else np.array(v_lim0, dtype=float, copy=True))

    u = np.asarray(circuit.u(0.0, epar, analysis='dc'), dtype=float)
    if source != 1.0:
        u = u * source
    shunt = (np.array([i for i in range(len(circuit.nodes)) if i != irefnode],
                      dtype=int) if gshunt else None)
    if anchor is not None and shunt is not None:
        ## (the shunt's anchor in the sources: gshunt x + u - gshunt anchor)
        u = u.copy()
        u[shunt] -= gshunt * np.asarray(anchor, dtype=float)[shunt]

    ## THE LOOP IN C (`_pcnr_c`, the PCNR-in-C plan's C2): from this start,
    ## bit for bit the loop below -- or None, declined or handed back, and
    ## the loop below runs
    from pycircuit.circuit import _pcnr_c
    got = _pcnr_c.solve(circuit, irefnode, devices, x, v_lim, u, epar, maxiter, reltol,
                        abstol, iabstol, gmin, gshunt, linsolver, jgmin, hrows)
    if got is not None:
        return got
    jpairs = ([(int(p[2]), int(p[3])) for p in pcnr_junction_pairs(circuit)] if jgmin
              else None)

    for it in range(maxiter):
        g_mna, g_lim, J_mm, J_ml, J_lm, didv = augmented_system(
            circuit, x, v_lim, devices, epar, u_extra=u, dense_blocks=False)
        if shunt is not None:
            ## (a conductance from every node to ground, IN the equations)
            g_mna[shunt] += gshunt * x[shunt]
            J_mm[shunt, shunt] += gshunt
        if jpairs is not None:
            ## (a conductance across every junction pair, IN the equations:
            ## C's `assemble`, in its order)
            for a, b in jpairs:
                vj = x[a] - x[b]
                g_mna[a] = g_mna[a] + jgmin * vj
                g_mna[b] = g_mna[b] - jgmin * vj
                J_mm[a, a] = J_mm[a, a] + jgmin
                J_mm[a, b] = J_mm[a, b] - jgmin
                J_mm[b, a] = J_mm[b, a] - jgmin
                J_mm[b, b] = J_mm[b, b] + jgmin

        dx_mna, dx_lim = predict(g_mna, g_lim, J_mm, J_ml, J_lm, irefnode,
                                 junctions=devices, didv=didv, gmin=gmin,
                                 linsolver=linsolver, toolkit=circuit.toolkit, held=hrows)

        x_new = x + dx_mna
        x_new[irefnode] = 0.0
        ## A non-finite unknown stays so (`x + dx`), and a NaN fails every
        ## test: the rung's failure now, not at the iteration limit
        if not np.all(np.isfinite(x_new)):
            raise RuntimeError(NONFINITE % (it + 1))
        v_new = v_lim + dx_lim

        ## CORRECT: each device limits only what it owns.
        v_new = refine(devices, v_lim, v_new, epar, x_old=x)

        keep = np.arange(n) != irefnode
        keep[list(hrows)] = False
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


#: `solve_dc_continued`'s gshunt ladder: the first rung, the last before the
#: pure solve, and the geometric halvings a failed rung is retried with
GSHUNT_START, GSHUNT_END, GSHUNT_RETRIES = 1e-2, 1e-12, 4

#: its junction-gmin ladder (`jgmin`): the first rung and the last before the
#: pure solve, by decades; a failed rung is retried as gshunt's
JGMIN_START, JGMIN_END = 1e-2, 1e-11

#: its pseudo-transient route (`_ptc_ladder`): the exponents of the anchored
#: conductance -- the first rung, the last before the pure solve, and the
#: largest an unconverged first rung escalates to (the chain's
#: `PseudoTransientNewton`'s)
PTC_E_START, PTC_E_END, PTC_E_MAX = 0.0, -12.0, 6.0


def solve_dc_continued(circuit, refnode, x0=None, epar=defaultepar, maxiter=200,
                       reltol=1e-6, abstol=1e-12, iabstol=1e-12, linsolver=None, held=None,
                       v_lim0=None, raw_seed=False):
    """The DC operating point by PCNR, with ladders of its own.

    Returns ``(x, v_lim, iterations, path)`` -- ``path`` is ``'plain'``,
    ``'gshunt'``, ``'jgmin'``, ``'source'`` or ``'ptc'``.  The plain :func:`solve_dc`
    first; where it raises, a gshunt ladder (:func:`_gshunt_ladder`) and the
    pure solve from its last rung; then a junction-gmin ladder
    (:func:`_jgmin_ladder`, true gmin stepping -- it deforms only the
    junctions, so it can come down where the node shunt's last rungs or its
    pure solve could not); then source stepping (:func:`_source_steps`),
    a step that fails retried by both ladders at its source level (the
    junction ladder first, :func:`_source_step`) -- the
    ordinary chain's shape, which a hard circuit needs: opampal fails the
    first rung of every ladder alone and reaches PCNR's answer only so;
    LAST, pseudo-transient continuation (:func:`_ptc_ladder`), the chain's
    own last route, which PCNR lacked: where every ladder above fails, the
    solve is not handed to the chain yet.
    Every rung is a PCNR solve seeded with the previous rung's ``x`` and
    limited unknowns, with the nodes ``held`` held (:func:`solve_dc`);
    ``v_lim0`` seeds the plain solve's limited unknowns -- a converged point
    is no wild start, and re-seeding from ``x0`` clamped would lose its
    state (`DC`'s nodeset release).  Raises the plain solve's exception
    where the ladders fail.

    THE ACCEPTANCE CHECK: every answer it would return -- the plain
    solve's, a ladder's pure solve -- must solve the DC equations to
    `DC_ACCEPT_RESIDUAL` (`dc_residual`); one that does not counts as that
    solve failing, and the next route is tried.

    ``raw_seed``: `x0` may be a solution -- a caller's operating point, a
    sweep's previous point (`DC`) -- so the plain solve is tried first from
    its UNCLAMPED seed (``v_lim_init(limit=False)``), where the clamp would
    move it; failing, the plain solve from the clamped seed and the ladders
    as above.  The clamp limits each junction's voltage as a step from 0 V
    (pnjlim takes a forward-biased 1.155 V LED to 0.098 V): from a settled
    point it switches the conducting junctions off, and the ladders may
    then converge to another root (an LED load's converter: vout -189 V
    where the chain finds 1.655 V, 2026-10-21).  A wild start's unclamped
    solve fails at once (the 20 V BJT mirror: a singular matrix at its first
    iterate) and the clamped solve follows.  (Ladders were tried once before, 2026-08-26, and failed on a BJT
    mirror from 20 V -- every rung then began from an UNLIMITED seed, since
    fixed in :func:`v_lim_init`.)
    """
    kw = {'epar': epar, 'maxiter': maxiter, 'reltol': reltol, 'abstol': abstol,
          'iabstol': iabstol, 'linsolver': linsolver}
    if held:
        kw['held'] = held
    rejected = []

    def solves(got):
        ## (the acceptance check: `got` where its point solves the DC
        ## equations, else None and the residual kept for the message)
        if got is None:
            return None
        r, row = dc_residual(circuit, refnode, got[0], epar, held)
        if r <= DC_ACCEPT_RESIDUAL:
            return got
        rejected.append((r, row))
        return None
    if raw_seed and v_lim0 is None and x0 is not None:
        ## (the caller's point as `solve_dc` starts from it: the reference
        ## row 0, the held rows at their volts)
        xs = np.array(x0, dtype=float, copy=True)
        xs[circuit.get_node_index(refnode)] = 0.0
        for r, volts in (held or {}).items():
            xs[r] = volts
        devices = pcnr_devices(circuit)
        raw = v_lim_init(devices, xs, epar)
        if not np.array_equal(raw, v_lim_init(devices, xs, epar, limit=True)):
            got = solves(_rung(circuit, refnode, x0, raw, kw))
            if got is not None:
                return got + ('plain',)
    try:
        got = solves(solve_dc(circuit, refnode, x0=x0, v_lim0=v_lim0, **kw))
        if got is not None:
            return got + ('plain',)
        first = RuntimeError(
            'PCNR converged to a point that does not solve the DC equations '
            f'(residual {rejected[-1][0]:.3g} at row {rejected[-1][1]})')
    except _RUNG_FAILURES as exc:
        first = exc
    got = _gshunt_ladder(circuit, refnode, x0, kw)
    if got is not None:
        x, v_lim, total = got
        pure = solves(_rung(circuit, refnode, x, v_lim, kw))
        if pure is not None:
            return pure[0], pure[1], total + pure[2], 'gshunt'
    got = solves(_jgmin_ladder(circuit, refnode, x0, kw))
    if got is not None:
        return got + ('jgmin',)
    got = solves(_source_steps(circuit, refnode, x0, kw))
    if got is not None:
        return got + ('source',)
    got = solves(_ptc_ladder(circuit, refnode, x0, kw))
    if got is not None:
        return got + ('ptc',)
    raise first


#: what a failed rung raises: no convergence, a singular Schur matrix, an
#: arithmetic fault
_RUNG_FAILURES = (RuntimeError, ArithmeticError, np.linalg.LinAlgError)

#: THE DC ACCEPTANCE CHECK's bound (`dc_residual`): the largest residual a DC
#: answer may leave -- A on a node row, V on a branch row.  The Newton tests
#: are per row and RELATIVE, ``|F| < reltol (|J||x| + |F|) + abstol``, and at
#: a point far off -- a converter's supply side at 1.4e10 V -- ``|J||x|`` is
#: ~2.4e91, so a KCL residual of 5.6e79 A passed both PCNR's test and the
#: chain's (2026-10-21; `DC()` from zeros returned it, the chain from zeros
#: finds 1.655 V).  A true root at reltol 1e-6 exceeds 1 only on a row
#: carrying over 1e6 A.
DC_ACCEPT_RESIDUAL = 1.0


def dc_residual(circuit, refnode, x, epar=defaultepar, held=None, u=None):
    """``(r, row)``: the largest ``|i(x) + u_dc|`` -- the DC equations, the
    devices at `x`'s node voltages (a stateful limiter synced to `x` for the
    read and put back: a read, not a limiting) -- over the rows a solve
    tests (not the reference's, not a `held` node's), and its row; ``(inf,
    row)`` where an entry is not finite, ``(inf, -1)`` where the evaluation
    raises.  `u`: the DC sources as the caller reads them (None: `solve_dc`'s
    call).  `DC`'s acceptance check (`DC_ACCEPT_RESIDUAL`), on PCNR's answers
    and the chain's."""
    from pycircuit.circuit._limiting import (
        limit_sync,
        limiter_snapshot,
        state_restore,
        stateful_limiters,
    )
    x = np.asarray(x, dtype=float)
    ## (the sources first: they do not read a limiter's state, and read after
    ## the sync writes it they would re-arm their stamp, `_stamp_plan`)
    if u is None:
        u = circuit.u(0.0, epar, analysis='dc')
    lims = stateful_limiters(circuit)
    snap = limiter_snapshot(lims) if lims else None
    try:
        if lims:
            limit_sync(circuit, x, epar, lims)
        with np.errstate(all='ignore'):
            f = np.abs(np.asarray(circuit.i(x, epar), dtype=float) + np.asarray(u, dtype=float))
    except (ArithmeticError, ValueError):
        return float('inf'), -1
    finally:
        if snap is not None:
            state_restore(snap)
    f[circuit.get_node_index(refnode)] = 0.0
    for r in (held or ()):
        f[r] = 0.0
    bad = np.flatnonzero(~np.isfinite(f))
    if len(bad):
        return float('inf'), int(bad[0])
    k = int(np.argmax(f))
    return float(f[k]), k


def _rung(circuit, refnode, x, v_lim, kw, **deform):
    """One rung: `(x, v_lim, iterations)`, or None where it raises."""
    try:
        return solve_dc(circuit, refnode, x0=x, v_lim0=v_lim, **kw, **deform)
    except _RUNG_FAILURES:
        return None


def _gshunt_ladder(circuit, refnode, x0, kw):
    """gshunt from `GSHUNT_START` down by decades to `GSHUNT_END`, a failed
    rung retried at the geometric mean with the last that converged:
    `(x, v_lim, iterations)` at the last rung (the pure solve is the
    caller's), or None where a rung failed on the way down."""
    got = _rung(circuit, refnode, x0, None, kw, gshunt=GSHUNT_START)
    if got is None:
        return None
    x, v_lim, total = got
    g = GSHUNT_START
    ## (with a relative slack: ten divisions by 10 leave 1.0000000000000002e-12,
    ## which IS the last decade -- compared exactly, it bought one more rung at
    ## GSHUNT_END from a point converged two ulps away, and gm1 failed it five
    ## times at the iteration limit)
    while g > GSHUNT_END * (1.0 + 1e-9):
        g_try = max(g / 10.0, GSHUNT_END)
        for _retry in range(GSHUNT_RETRIES + 1):
            got = _rung(circuit, refnode, x, v_lim, kw, gshunt=g_try)
            if got is not None:
                break
            g_try = (g_try * g) ** 0.5
        if got is None:
            return None
        x, v_lim, its = got
        total += its
        g = g_try
    return x, v_lim, total


#: `_source_steps`: the first step, the largest, and the smallest before it
#: gives up -- three failing steps from the start (0.05, 0.025, 0.0125).  0.05,
#: not 0.1: on opampal, the one gold deck whose operating point steps its
#: sources, every try at 0.1 failed (~300 iterations); with the junction-gmin
#: ladder first in a failed step (`_source_step`) its operating point went
#: from 1.90 to 0.99 s, the answer within 3.7e-12 V.  A smallest step of 1e-3
#: spent seven 100-iteration failures on each deck it did not rescue
SOURCE_STEP, SOURCE_STEP_MAX, SOURCE_STEP_MIN = 0.05, 0.5, 0.0125



def _jgmin_ladder(circuit, refnode, x0, kw):
    """`jgmin` from `JGMIN_START` down by decades to `JGMIN_END`, a failed
    rung retried at the geometric mean with the last that converged, then
    the pure solve: `(x, v_lim, iterations)`, or None -- at once where the
    circuit has no junction pair."""
    if not pcnr_junction_pairs(circuit):
        return None
    got = _rung(circuit, refnode, x0, None, kw, jgmin=JGMIN_START)
    if got is None:
        return None
    x, v_lim, total = got
    g = JGMIN_START
    while g > JGMIN_END * (1 + 1e-9):
        g_try = max(g / 10.0, JGMIN_END)
        for _r in range(GSHUNT_RETRIES + 1):
            got = _rung(circuit, refnode, x, v_lim, kw, jgmin=g_try)
            if got is not None:
                break
            g_try = (g_try * g) ** 0.5
        if got is None:
            return None
        x, v_lim, its = got
        total += its
        g = g_try
    pure = _rung(circuit, refnode, x, v_lim, kw)
    if pure is None:
        return None
    return pure[0], pure[1], total + pure[2]


def _ptc_ladder(circuit, refnode, x0, kw):
    """PSEUDO-TRANSIENT CONTINUATION, PCNR's last route (2026-10-21): the
    chain's `PseudoTransientNewton` -- backward-Euler pseudo-time steps, each
    a gshunt rung ANCHORED at its seed, ``g (x - x_seed)`` on the node rows
    (`solve_dc`'s `anchor`), so the homotopy follows the pseudo-transient
    trajectory from the start rather than a deformation anchored at 0 V --
    on its schedule (`nrsolver._adaptive_conductance_ladder`: g from 1 down
    to 1e-12 by halving steps in exponent space, escalating to 1e6 before
    the first rung converges), the limited unknowns carried from rung to
    rung, then the pure solve: ``(x, v_lim, iterations)``, or None.

    Last because it is the most robust and the slowest: measured alone on
    the gold decks (scratch, 2026-10-21) it solved 24 of 27 -- not e1480 and
    opampal, which the ladders above solve -- and took slowlatch and todd3
    to the other state of their latches (it follows the trajectory from the
    zeros seed); last, no gold deck reaches it."""
    from pycircuit.circuit.analysis import NoConvergenceError
    from pycircuit.circuit.nrsolver import _adaptive_conductance_ladder
    state = {'v_lim': None, 'its': 0}

    def rung(x_seed, g):
        a = np.zeros(circuit.n) if x_seed is None else np.asarray(x_seed, dtype=float)
        got = _rung(circuit, refnode, a, state['v_lim'], kw, gshunt=g, anchor=a)
        if got is None:
            raise NoConvergenceError(f'pseudo-transient rung at g={g:g}')
        state['v_lim'] = got[1]
        state['its'] += got[2]
        return got[0]

    def pure(x_seed):
        got = _rung(circuit, refnode, x_seed, state['v_lim'], kw)
        if got is None:
            raise NoConvergenceError('pseudo-transient: the pure solve')
        return got[0], got[1], state['its'] + got[2]

    try:
        return _adaptive_conductance_ladder(rung, pure, x0, e_start=PTC_E_START,
                                            e_end=PTC_E_END, e_max=PTC_E_MAX,
                                            label='PCNR pseudo-transient continuation')
    except NoConvergenceError:
        return None


def _source_step(circuit, refnode, x, v_lim, kw, target):
    """One source step: the plain rung, else the junction-gmin ladder, else
    the gshunt ladder and its pure solve -- all at `target`'s source level:
    `(x, v_lim, iterations)`, or None.  The junction ladder FIRST here (the
    route's own order is gshunt first): on opampal, at its source 0.05, the
    gshunt ladder spent ~900 iterations in a dead end (its 1e-9 rung and eight
    retries at the iteration limit) before the junction ladder succeeded in
    ~150."""
    got = _rung(circuit, refnode, x, v_lim, kw, source=target)
    if got is not None:
        return got
    kws = dict(kw, source=target)
    got = _jgmin_ladder(circuit, refnode, x, kws)
    if got is not None:
        return got
    lad = _gshunt_ladder(circuit, refnode, x, kws)
    if lad is not None:
        pure = _rung(circuit, refnode, lad[0], lad[1], kws)
        if pure is not None:
            return pure[0], pure[1], lad[2] + pure[2]
    return None


def _source_steps(circuit, refnode, x0, kw):
    """The sources raised from 0 to 1, each step by :func:`_source_step`,
    the step doubled after one that converged and halved after one that did
    not: `(x, v_lim, iterations)` at full sources, or None."""
    n = circuit.n
    x, v_lim, lam, step, total = np.zeros(n), None, 0.0, SOURCE_STEP, 0
    while lam < 1.0:
        target = min(1.0, lam + step)
        got = _source_step(circuit, refnode, x, v_lim, kw, target)
        if got is None:
            step /= 2.0
            if step < SOURCE_STEP_MIN:
                return None
            continue
        x, v_lim, its = got
        total += its
        lam, step = target, min(step * 2.0, SOURCE_STEP_MAX)
    return x, v_lim, total

