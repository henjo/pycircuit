"""THE CONSTANT-STAMP PLAN: a SubCircuit's assembly that does not re-stamp
the elements whose stamps cannot change (2026-10-01; the speed plan's P2).

Assembly was a Python loop over every element on every pass -- `x[nodemap]`,
a method call that built or fetched a stamp, list appends, then one
`bincount` -- five passes per Newton iterate (C, q, u, i, G), 73-86 % of a
transient on an ordinary circuit, where most elements (R, C, L, sources) have
stamps that never change.  This plan keeps that arithmetic and drops the
loop for them.  BIT-IDENTICAL BY CONSTRUCTION:

* A MATRIX pass (G, C) is still one `bincount` over the same values in the
  same order: the constant elements' stamp entries sit pre-filled in a
  template buffer, the others are called in order and written into their
  slots.  Exact zeros of a constant stamp are left out: a bin starts at
  +0.0 and so is never -0.0, and adding +-0.0 to it changes no bit.
* A VECTOR pass (i, q) computes a constant element's ``dot(S, x[nodemap])``
  -- exactly what `Circuit.i` / `Circuit.q` do -- as ONE `np.matmul` over
  every constant element of one stamp size: it reaches the same BLAS gemv
  per element, so it is bit-identical (measured 0 of 20000 differ; a
  per-process self-check per stamp size falls back to per-element `np.dot`
  if a numpy or BLAS ever made them differ).  All-zero stamp rows give +-0
  entries and are left out -- exact for a FINITE `x` only (`0 * inf` is
  NaN), so a non-finite `x` takes the legacy loop.

WHICH ELEMENTS ARE CONSTANT, per method pair (G with i, C with q), decided by
METHOD IDENTITY: the partner (`i` / `q`) is exactly `Circuit`'s default and
the matrix method is `Circuit`'s zero one, or the method of the nearest class
declaring `_constant_stamps` naming it (a stamp built only in `update()`,
independent of x, t and epar).  An override breaks the identity, so a
subclass that computes its own `G` drops out without saying so.
`Circuit.linear` is NOT used: linear is not constant (a TLine is linear and
re-stamped from its history), and BSource and the switches said linear until
2026-10-02.  Per element too: no instance attribute shadowing the methods,
the circuit's own toolkit, a real `(k, k)` ndarray stamp.

WHEN IT RUNS: the numeric toolkit exactly (not sparse, JAX or symbolic), no
`params_tree`, methods G, C, i, q (never `CY`), a 1-D float64 ndarray `x` of
the circuit's size, no `dtype` on a vector pass.  Anything else, and any
non-real stamp a non-constant element returns, takes the legacy loop.

WHEN IT IS REBUILT: the parameter epoch moved (`ParameterDict._epoch`, a new
sentinel on every notification an element observes -- every parameter write
goes through it), the element node map is a new object (topology), or the
circuit's size changed.  `invalidate()` for what none of those see (a stamp
mutated in place).  `ENABLED` (env `PYCIRCUIT_STAMP_PLAN=0` turns it off) is
read on every pass, so a measurement can A/B it.

THE C-BOUND HDL CLASSES ARE ONE CALL EACH (`_hdl_batch`, 2026-10-03; speed
round 6): the non-constant elements of a chained hdl class whose kernels are
bound run through one C driver per class per pass -- the kernel's own
function pointer on each element's `x[nm]` and the element's own pack, the
outputs written into the elements' slots -- instead of a wrapped call per
element through the loop above (56 -> 7.5 us for a 20-MosLevel1 `G` pass).
Which elements, and what is checked on every pass, is that module's note;
an element a batch hands back is called here as before, in its order.

History: `doc/transient_history.md`, `_stamp_plan`.
"""
import os

import numpy as np

from pycircuit.circuit import _hdl_batch

ENABLED = os.environ.get('PYCIRCUIT_STAMP_PLAN', '1') != '0'

_PARTNER = {'G': 'i', 'C': 'q'}
_MATRIX = {'i': 'G', 'q': 'C'}

#: per stamp size, whether a batched `np.matmul` reproduced per-element
#: `np.dot` in this process (checked once on the first stamps of that size)
_MATMUL_OK = {}


def pair_kind(cls, m):
    """'zero', 'cached' or None for the matrix method `m` ('G' or 'C') of the
    element class `cls` -- see the module note."""
    from pycircuit.circuit.circuit import Circuit
    p = _PARTNER[m]
    if getattr(cls, p) is not getattr(Circuit, p):
        return None
    f = getattr(cls, m)
    if f is getattr(Circuit, m):
        return 'zero'
    for k in cls.__mro__:
        decl = k.__dict__.get('_constant_stamps')
        if decl is not None:
            return 'cached' if (m in decl and f is getattr(k, m)) else None
    return None


def _element_kind(cir, el, m, k):
    """`pair_kind` for this element, with the per-element checks; and its
    stamp when 'cached'.  Returns `(kind, stamp)`."""
    if k == 0 or el.toolkit is not cir.toolkit:
        return None, None
    if m in el.__dict__ or _PARTNER[m] in el.__dict__:
        return None, None
    kind = pair_kind(type(el), m)
    if kind != 'cached':
        return kind, None
    from pycircuit.circuit.circuit import defaultepar
    S = getattr(el, m)(np.zeros(k), defaultepar)
    if (not isinstance(S, np.ndarray) or S.shape != (k, k)
            or S.dtype.kind not in 'fiub'):
        return None, None
    return 'cached', S


class _MatrixPlan:
    """One matrix method's plan: the template buffer, the flat indices, the
    non-constant calls with their slots (see the module note)."""

    __slots__ = ('batches', 'calls', 'flat', 'template')

    def __init__(self, cir, m):
        n = cir.n
        idxmap = cir._map_indices_2d
        nodemaps = cir.elementnodemap
        flats, parts, calls = [], [], []
        pos = 0
        for inst, el in cir.elements.items():
            nm = nodemaps[inst]
            rc = idxmap.get(inst)
            kind, S = _element_kind(cir, el, m, len(nm))
            if kind == 'zero':
                continue
            if kind == 'cached':
                if rc is None:
                    continue
                v = np.asarray(S, dtype=np.float64).ravel()
                keep = v != 0.0
                rows, cols = rc
                flats.append((np.asarray(rows, dtype=np.intp) * n
                              + np.asarray(cols, dtype=np.intp))[keep])
                parts.append(v[keep])
                pos += int(keep.sum())
                continue
            if rc is None:
                calls.append((inst, el, nm, -1, -1))
                continue
            rows, cols = rc
            size = len(rows)
            flats.append(np.asarray(rows, dtype=np.intp) * n
                         + np.asarray(cols, dtype=np.intp))
            parts.append(np.zeros(size))
            calls.append((inst, el, nm, pos, pos + size))
            pos += size
        self.flat = (np.concatenate(flats) if flats
                     else np.zeros(0, dtype=np.intp))
        self.template = np.concatenate(parts) if parts else np.zeros(0)
        self.calls, self.batches = _hdl_batch.split(cir, m, calls)


class _VectorPlan:
    """One vector method's plan: the constant elements grouped by stamp size
    for the batched product, the non-constant calls, the flat indices."""

    __slots__ = ('batches', 'calls', 'groups', 'idx', 'length')

    def __init__(self, cir, v):
        m = _MATRIX[v]
        idxmap = cir._map_indices_1d
        nodemaps = cir.elementnodemap
        idxs, calls = [], []
        by_k = {}
        pos = 0
        for inst, el in cir.elements.items():
            nm = nodemaps[inst]
            ind = idxmap.get(inst)
            kind, S = _element_kind(cir, el, m, len(nm))
            if kind == 'zero':
                continue
            if kind == 'cached' and len(nm) >= 2 and ind is not None:
                rows = np.flatnonzero(np.any(np.asarray(S) != 0, axis=1))
                if rows.size:
                    g = by_k.setdefault(len(nm), [[], [], [], []])
                    g[0].append(np.asarray(S, dtype=np.float64))
                    g[1].append(np.asarray(nm, dtype=np.intp))
                    g[2].append(rows + len(g[2]) * len(nm))
                    g[3].append(np.arange(pos, pos + rows.size))
                    idxs.append(np.asarray(ind, dtype=np.intp)[rows])
                    pos += rows.size
                continue
            if kind == 'cached':
                ## (a 1x1 stamp: `np.dot` on it is a scalar product, not
                ## gemv -- called per element like any other)
                pass
            if ind is None:
                calls.append((inst, el, nm, -1, -1))
                continue
            k = len(ind)
            idxs.append(np.asarray(ind, dtype=np.intp))
            calls.append((inst, el, nm, pos, pos + k))
            pos += k
        groups = []
        for k, (Ss, IDX, sel, dst) in by_k.items():
            S = np.ascontiguousarray(np.stack(Ss))
            I = np.stack(IDX)
            sel = np.concatenate(sel)
            dst = np.concatenate(dst)
            groups.append((S, I, sel, dst, _matmul_ok(k, S)))
        self.idx = (np.concatenate(idxs) if idxs
                    else np.zeros(0, dtype=np.intp))
        self.length = pos
        self.groups = groups
        self.calls, self.batches = _hdl_batch.split(cir, v, calls)


def _matmul_ok(k, S):
    """Whether a batched `np.matmul` reproduces per-element `np.dot` for
    stamps of size `k` in this process -- checked once, on these stamps."""
    ok = _MATMUL_OK.get(k)
    if ok is None:
        rng = np.random.default_rng(0)
        X = rng.standard_normal((S.shape[0], k)) * 10.0 ** rng.integers(
            -8, 3, (S.shape[0], 1))
        with np.errstate(all='ignore'):
            batched = np.matmul(S, X[:, :, None])[:, :, 0]
        single = np.array([np.dot(S[e], X[e]) for e in range(S.shape[0])])
        ok = _MATMUL_OK[k] = bool(np.array_equal(batched, single))
    return ok


class _Plan:
    __slots__ = ('builds', 'elements', 'epoch', 'methods', 'n', 'nodemap')

    def __init__(self, cir, epoch, builds):
        self.epoch = epoch
        ## (the ELEMENTS dict too: a shallow copy may keep the node map object
        ## and hold other elements)
        self.elements = cir.elements
        self.nodemap = cir.elementnodemap
        self.n = cir.n
        self.methods = {}
        self.builds = builds


def _plan_for(cir):
    """The circuit's current plan, rebuilt when it went stale."""
    from pycircuit.utilities.param import ParameterDict
    epoch = ParameterDict._epoch
    plan = cir.__dict__.get('_stamp_plan')
    if (plan is None or plan.epoch is not epoch
            or plan.elements is not cir.elements
            or plan.nodemap is not cir.elementnodemap
            or plan.n != len(cir.nodes) + len(cir.branches)):
        plan = _Plan(cir, epoch, 1 if plan is None else plan.builds + 1)
        cir.__dict__['_stamp_plan'] = plan
    return plan


def invalidate(cir):
    """Drop `cir`'s plan: for a change the plan cannot see (a stamp array
    mutated in place, a method patched on a class)."""
    cir.__dict__.pop('_stamp_plan', None)


def _eligible(cir, x):
    from pycircuit.circuit.toolkit import NumericToolkit
    return (ENABLED and type(cir.toolkit) is NumericToolkit
            and type(x) is np.ndarray and x.dtype == np.float64
            and x.ndim == 1 and x.shape[0] == len(cir.nodes) + len(cir.branches))


def _call(el, m, subx, args):
    ## the legacy MATRIX loop's wrapping, so a failing element reads the same
    ## (the vector loop never wrapped, and neither does `assemble_vector`)
    try:
        return getattr(el, m)(subx, *args)
    except Exception as e:  # noqa: BLE001 -- re-raised, as the loop does
        raise e.__class__(str(e) + ' at element ' + str(el)
                          + ', args=' + str(args))


def _call_plain(el, v, subx, args):
    return getattr(el, v)(subx, *args)


_DEFAULT_EPAR = None


def _run_batches(x, args, p, buf, got, call, m):
    """The plan's batches (`_hdl_batch`) into `buf`: True, or False with
    `got` extended by everything computed so far, for the legacy loop.  An
    element a batch hands back is called here, in its order."""
    global _DEFAULT_EPAR
    if args:
        epar = args[0]
    else:
        if _DEFAULT_EPAR is None:
            from pycircuit.circuit.circuit import defaultepar
            _DEFAULT_EPAR = defaultepar
        epar = _DEFAULT_EPAR
    filled = []
    for bt in p.batches:
        res = bt.run(x, epar)
        served = ()
        if res is not None:
            OUT, pos = res
            if pos is None:
                buf[bt.dst] = OUT.reshape(-1)
                filled.extend(bt.slots)
                continue
            for j, e in enumerate(pos):
                a, b = bt.slots[e]
                buf[a:b] = OUT[j].reshape(-1)
                filled.append((a, b))
            served = set(pos)
        for e, (inst, el, nm, a, b) in enumerate(bt.entries):
            if e in served:
                continue
            rhs = call(el, m, x[nm], args)
            val = np.asarray(rhs).ravel()
            got.append((a, b, rhs))
            if val.dtype.kind not in 'fiub' or val.size != b - a:
                for a2, b2 in filled:
                    got.append((a2, b2, buf[a2:b2]))
                return False
            buf[a:b] = val
    return True


def _every_call(p):
    """The plan's non-constant calls, batched or not."""
    yield from p.calls
    for bt in p.batches:
        yield from bt.entries


def assemble_matrix(cir, m, x, args):
    """`m` ('G' or 'C') at `x` through the plan, or None for the legacy loop."""
    if not _eligible(cir, x):
        return None
    plan = _plan_for(cir)
    mp = plan.methods.get(m)
    if mp is None:
        mp = plan.methods[m] = _MatrixPlan(cir, m)
    got = []
    for inst, el, nm, a, b in mp.calls:
        rhs = _call(el, m, x[nm], args)
        if a >= 0:
            got.append((a, b, rhs))
    buf = mp.template.copy()
    for a, b, rhs in got:
        val = np.asarray(rhs).ravel()
        if val.dtype.kind not in 'fiub' or val.size != b - a:
            return _legacy_matrix(cir, m, x, args, mp, got)
        buf[a:b] = val
    if mp.batches and not _run_batches(x, args, mp, buf, got, _call, m):
        return _legacy_matrix(cir, m, x, args, mp, got)
    n = plan.n
    if not mp.flat.size:
        ## (an EMPTY bincount is int64 even with weights: the loop's zeros)
        return np.zeros((n, n))
    return np.bincount(mp.flat, weights=buf, minlength=n * n).reshape(n, n)


def assemble_vector(cir, v, x, args):
    """`v` ('i' or 'q') at `x` through the plan, or None for the legacy loop."""
    if not _eligible(cir, x):
        return None
    if not np.isfinite(x).all():
        return None
    plan = _plan_for(cir)
    vp = plan.methods.get(v)
    if vp is None:
        vp = plan.methods[v] = _VectorPlan(cir, v)
    got = []
    for inst, el, nm, a, b in vp.calls:
        rhs = getattr(el, v)(x[nm], *args)
        if a >= 0:
            got.append((a, b, rhs))
    buf = np.empty(vp.length)
    for S, I, sel, dst, ok in vp.groups:
        if ok:
            with np.errstate(all='ignore'):
                Y = np.matmul(S, x[I][:, :, None])
        else:
            Y = np.array([np.dot(S[e], x[I[e]]) for e in range(S.shape[0])])
        buf[dst] = Y.reshape(-1)[sel]
    for a, b, rhs in got:
        val = np.asarray(rhs).ravel()
        if val.dtype.kind not in 'fiub' or val.size != b - a:
            return _legacy_vector(cir, v, x, args, vp, got)
        buf[a:b] = val
    if vp.batches and not _run_batches(x, args, vp, buf, got, _call_plain, v):
        return _legacy_vector(cir, v, x, args, vp, got)
    if not vp.idx.size:
        return np.zeros(plan.n)
    return np.bincount(vp.idx, weights=buf, minlength=plan.n)


def _legacy_matrix(cir, m, x, args, mp, got):
    """A non-real stamp from a non-constant element: the legacy pending lists,
    the constant elements re-called (they are pure), the non-constant ones'
    values reused -- never called twice."""
    computed = {a: rhs for a, b, rhs in got}
    by_inst = {inst: a for inst, el, nm, a, b in _every_call(mp)}
    idxmap = cir._map_indices_2d
    pending_rc, pending_val = [], []
    for inst, el in cir.elements.items():
        rc = idxmap.get(inst)
        if rc is None:
            continue
        a = by_inst.get(inst)
        rhs = (computed[a] if a is not None and a in computed
               else getattr(el, m)(x[cir.elementnodemap[inst]], *args))
        pending_rc.append(rc)
        pending_val.append(np.asarray(rhs).ravel())
    n = cir.n
    return cir._scatter_2d(cir.toolkit.zeros((n, n)), pending_rc, pending_val,
                           n)


def _legacy_vector(cir, v, x, args, vp, got):
    """The vector pass's fallback, as `_legacy_matrix`."""
    computed = {a: rhs for a, b, rhs in got}
    by_inst = {inst: a for inst, el, nm, a, b in _every_call(vp)}
    idxmap = cir._map_indices_1d
    pending_idx, pending_val = [], []
    for inst, el in cir.elements.items():
        ind = idxmap.get(inst)
        if ind is None:
            continue
        a = by_inst.get(inst)
        rhs = (computed[a] if a is not None and a in computed
               else getattr(el, v)(x[cir.elementnodemap[inst]], *args))
        pending_idx.append(ind)
        pending_val.append(np.asarray(rhs).ravel())
    n = cir.n
    return cir._scatter_1d(cir.toolkit.zeros(n), pending_idx, pending_val, n)
