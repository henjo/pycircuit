"""Prototype: what does a hoisted + bincount assembly buy on the leapfrog?

Monkeypatches SubCircuit._add_element_subvectors / _add_element_submatrices with
a version that (a) hoists the per-element hasattr()/getattr() probes out of the
loop and (b) accumulates with np.bincount instead of np.add.at.  Correctness is
checked against the original on every call before timing.
"""
import sys, time
import numpy as np
sys.path.insert(0, '/home/andreas/sources/pycircuit')

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, gnd, numeric
from pycircuit.circuit.circuit import SubCircuit as SC
from pycircuit.circuit.benchmark_circuits import build_leapfrog_network
from pycircuit.circuit.elements import BSource
from pycircuit.circuit.transient import Transient

F1 = 8e3
G3 = 1e-3

def build(stages=5):
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        node, amp, _ = build_leapfrog_network(cir, stages=stages, tones=((1e-3, F1),))
        cir['nl'] = BSource(amp[0]['e1'], gnd, amp[0]['e1'], gnd,
                            i_func=lambda v: G3 * v ** 3)
        return cir
    finally:
        circuit_module.default_toolkit = saved

orig_vec = SC._add_element_subvectors
orig_mat = SC._add_element_submatrices


def fast_subvectors(self, methodname, x, args, dtype=None, params_tree=None):
    n = self.n
    groups = getattr(self, '_eval_groups', None)
    elementnodemap = self.elementnodemap
    idxmap = self._map_indices_1d
    vals, idxs = [], []
    for instance, element in self.elements.items():
        if groups and element.__class__ in groups and x is not None:
            continue
        if x is not None:
            rhs = getattr(element, methodname)(x[elementnodemap[instance]], *args)
        else:
            rhs = getattr(element, methodname)(*args)
        ind = idxmap.get(instance)
        if ind is None:
            continue
        vals.append(np.asarray(rhs).ravel())
        idxs.append(ind)
    if not vals:
        return np.zeros(n, dtype=dtype)
    allv = np.concatenate(vals)
    alli = np.concatenate(idxs)
    if allv.dtype == object or (dtype is not None and np.dtype(dtype).kind == 'c'):
        lhs = np.zeros(n, dtype=dtype if dtype is not None else object)
        np.add.at(lhs, alli, allv)
        return lhs
    return np.bincount(alli, weights=allv, minlength=n)


def fast_submatrices(self, methodname, x, args, params_tree=None):
    n = self.n
    groups = getattr(self, '_eval_groups', None)
    elementnodemap = self.elementnodemap
    idxmap = self._map_indices_2d
    vals, rows, cols = [], [], []
    for instance, element in self.elements.items():
        if groups and element.__class__ in groups:
            continue
        if x is not None:
            rhs = getattr(element, methodname)(x[elementnodemap[instance]], *args)
        else:
            rhs = getattr(element, methodname)(*((None,) + tuple(args)))
        rc = idxmap.get(instance)
        if rc is None:
            continue
        vals.append(np.asarray(rhs).ravel())
        rows.append(rc[0]); cols.append(rc[1])
    lhs = np.zeros((n, n))
    if not vals:
        return lhs
    allv = np.concatenate(vals)
    flat = np.concatenate(rows) * n + np.concatenate(cols)
    return np.bincount(flat, weights=allv, minlength=n * n).reshape(n, n)


cir = build()
x = np.random.rand(cir.n) * 0.01
# correctness
for name, o, f in (('i', orig_vec, fast_subvectors), ('q', orig_vec, fast_subvectors)):
    a = o(cir, name, x, (cir.iparv,))
print('checking correctness...')
from pycircuit.circuit.circuit import defaultepar
for name in ('i', 'q'):
    a = orig_vec(cir, name, x, (defaultepar,))
    b = fast_subvectors(cir, name, x, (defaultepar,))
    print('  %s max abs diff %.3e' % (name, np.max(np.abs(a - b))))
for name in ('G', 'C'):
    a = orig_mat(cir, name, x, (defaultepar,))
    b = fast_submatrices(cir, name, x, (defaultepar,))
    print('  %s max abs diff %.3e' % (name, np.max(np.abs(a - b))))
a = orig_vec(cir, 'u', None, (0.0, defaultepar, 'tran'))
b = fast_subvectors(cir, 'u', None, (0.0, defaultepar, 'tran'))
print('  u max abs diff %.3e' % np.max(np.abs(a - b)))

# timing single assemblies
def bench(f, *a, **k):
    N = 30
    t0 = time.perf_counter()
    for _ in range(N):
        f(*a, **k)
    return (time.perf_counter() - t0) / N

print('\nper-call assembly time (s):')
for name in ('i', 'q'):
    print('  %s: orig %.5f  fast %.5f  speedup %.1fx' % (
        name, bench(orig_vec, cir, name, x, (defaultepar,)),
        bench(fast_subvectors, cir, name, x, (defaultepar,)),
        bench(orig_vec, cir, name, x, (defaultepar,)) / bench(fast_subvectors, cir, name, x, (defaultepar,))))
for name in ('G', 'C'):
    to = bench(orig_mat, cir, name, x, (defaultepar,))
    tf = bench(fast_submatrices, cir, name, x, (defaultepar,))
    print('  %s: orig %.5f  fast %.5f  speedup %.1fx' % (name, to, tf, to / tf))

# full transient A/B
def run():
    c = build()
    t0 = time.perf_counter()
    r = Transient(c, toolkit=numeric).solve(refnode=gnd, tend=20e-6,
                                            timestep=1.0 / (F1 * 128))
    return time.perf_counter() - t0, r

t_orig, r_orig = run()
SC._add_element_subvectors = fast_subvectors
SC._add_element_submatrices = fast_submatrices
t_fast, r_fast = run()
SC._add_element_subvectors = orig_vec
SC._add_element_submatrices = orig_mat

print('\nfull transient: orig %.2f s, fast %.2f s -> %.2fx' % (t_orig, t_fast, t_orig / t_fast))
print('waveform max abs diff: %.3e' % np.max(np.abs(np.asarray(r_orig.x) - np.asarray(r_fast.x))))
