"""THE STAMP PLANS BUILT BY SIZE (`_stamp_plan._MatrixPlan`, `_VectorPlan`,
`_hdl_batch.split` and `Batch`; plan 261's L7d): one pass in the elements'
order calls each element's own stamp, then each stamp size's nonzeros, rows
and slots are read at once and placed at each element's offset -- where every
element took its own small numpy operations (24 000 elements a plan on an
8 000-stage chain: 0.09 s a plan, eight plans a run).  The plans are the
per-element build's, entry for entry: the reference below is that build."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _hdl_batch, _stamp_plan, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import VS, C, L, R, SubCircuit, VPulse, gnd


def _ref_matrix(cir, m):
    """The per-element build of a matrix plan: `(flat, template, calls)`."""
    n = cir.n
    idxmap, nodemaps = cir._map_indices_2d, cir.elementnodemap
    flats, parts, calls, pos = [], [], [], 0
    for inst, el in cir.elements.items():
        nm, rc = nodemaps[inst], idxmap.get(inst)
        kind, S = _stamp_plan._element_kind(cir, el, m, len(nm))
        if kind == 'zero':
            continue
        if kind == 'cached':
            if rc is None:
                continue
            v = np.asarray(S, dtype=np.float64).ravel()
            keep = v != 0.0
            flats.append((np.asarray(rc[0], dtype=np.intp) * n
                          + np.asarray(rc[1], dtype=np.intp))[keep])
            parts.append(v[keep])
            pos += int(keep.sum())
            continue
        if rc is None:
            calls.append((inst, el, nm, -1, -1))
            continue
        size = len(rc[0])
        flats.append(np.asarray(rc[0], dtype=np.intp) * n + np.asarray(rc[1], dtype=np.intp))
        parts.append(np.zeros(size))
        calls.append((inst, el, nm, pos, pos + size))
        pos += size
    flat = np.concatenate(flats) if flats else np.zeros(0, dtype=np.intp)
    return flat, (np.concatenate(parts) if parts else np.zeros(0)), calls


def _ref_vector(cir, v):
    """The per-element build of a vector plan: `(idx, length, groups, calls)`."""
    m = _stamp_plan._MATRIX[v]
    idxmap, nodemaps = cir._map_indices_1d, cir.elementnodemap
    idxs, calls, by_k, pos = [], [], {}, 0
    for inst, el in cir.elements.items():
        nm, ind = nodemaps[inst], idxmap.get(inst)
        kind, S = _stamp_plan._element_kind(cir, el, m, len(nm))
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
        if ind is None:
            calls.append((inst, el, nm, -1, -1))
            continue
        idxs.append(np.asarray(ind, dtype=np.intp))
        calls.append((inst, el, nm, pos, pos + len(ind)))
        pos += len(ind)
    groups = [(np.ascontiguousarray(np.stack(Ss)), np.stack(I), np.concatenate(s),
               np.concatenate(d)) for _k, (Ss, I, s, d) in by_k.items()]
    idx = np.concatenate(idxs) if idxs else np.zeros(0, dtype=np.intp)
    return idx, pos, groups, calls


def _same_array(a, b):
    a, b = np.asarray(a), np.asarray(b)
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def _same_calls(a, b):
    return (len(a) == len(b) and all(p[0] == q[0] and p[1] is q[1] and p[2] is q[2]
                                     and p[3:] == q[3:] for p, q in zip(a, b)))


def _circuits():
    """A pulsed MOS chain with coupling capacitors and junctions; and a small
    circuit where a stamp size's first constant stamp is all zero (a 0 F
    capacitor) before another size's (an inductor's branch), and an element
    whose method is shadowed on the instance."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=2e-8, per=5e-8)
    for k in range(60):
        g = 'g0' if k == 0 else f'd{k - 1}'
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', g, gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
        c[f'cm{k}'] = C(f'd{k}', g, c=2e-15)
        if k % 5 == 0:
            c[f'D{k}'] = eh.DiodeSpiceHdl(gnd, f'd{k}', cjo=2e-14, vj=0.7, m=0.5)
    s = SubCircuit()
    s['v1'] = VS('a', gnd, v=1.0)
    s['c0'] = C('a', 'b', c=0.0)
    s['l1'] = L('b', 'c', L=1e-9)
    s['r1'] = R('c', gnd, r=50.0)
    s['c1'] = C('c', gnd, c=1e-12)
    s['r2'] = R('a', 'b', r=10.0)
    s['r2'].G = R.G.__get__(s['r2'])
    for cir in (c, s):
        cir.toolkit = circuit.numeric
    return [('chain', c), ('mixed', s)]


@pytest.mark.parametrize('name, cir', _circuits(), ids=lambda v: v if isinstance(v, str) else '')
def test_the_plans_by_size_are_the_per_element_ones(name, cir):
    """Each matrix plan's flat indices, template and calls (the same objects,
    the same slots) and each vector plan's indices, length, groups (in their
    order) and calls are the per-element build's."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for m in ('G', 'C'):
            got = _stamp_plan._MatrixPlan(cir, m)
            flat, template, calls = _ref_matrix(cir, m)
            want_calls, _b = _hdl_batch.split(cir, m, calls)
            assert _same_array(got.flat, flat) and _same_array(got.template, template), (name, m)
            assert _same_calls(got.calls, want_calls), (name, m)
        for v in ('i', 'q'):
            got = _stamp_plan._VectorPlan(cir, v)
            idx, length, groups, calls = _ref_vector(cir, v)
            want_calls, _b = _hdl_batch.split(cir, v, calls)
            assert _same_array(got.idx, idx) and got.length == length, (name, v)
            assert len(got.groups) == len(groups), (name, v)
            for g, w in zip(got.groups, groups):
                assert all(_same_array(a, b) for a, b in zip(g[:4], w)), (name, v)
            assert _same_calls(got.calls, want_calls), (name, v)


def test_a_zero_stamp_first_does_not_order_the_groups():
    """The vector plan's groups come in the order their sizes' first kept
    stamp does: the 0 F capacitor (all zero, size 2) before the inductor
    (size 3) does not put the size-2 group first."""
    _n, cir = _circuits()[1]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        got = _stamp_plan._VectorPlan(cir, 'q')
    sizes = [g[0].shape[1] for g in got.groups]
    _idx, _len, groups, _calls = _ref_vector(cir, 'q')
    assert sizes == [g[0].shape[1] for g in groups] and len(sizes) >= 1


def test_a_batch_is_the_per_element_one():
    """A batch's node maps and destination slots are the per-element
    stacking's: one row an element, its slots' indices in order."""
    _n, cir = _circuits()[0]
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        plan = _stamp_plan._MatrixPlan(cir, 'G')
    assert plan.batches, 'the MOS chain batches its transistors'
    for b in plan.batches:
        NM = np.ascontiguousarray(np.stack([np.asarray(e[2], dtype=np.int64) for e in b.entries]))
        dst = np.concatenate([np.arange(a, z) for a, z in b.slots])
        assert _same_array(b.NM, NM) and _same_array(b.dst, dst)
