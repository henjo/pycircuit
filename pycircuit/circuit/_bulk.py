"""THE BULK LINEAR ELEMENTS (plan 261 L21): many two-terminal linear
elements of one kind held as arrays in one `Bank` -- 'R' (each a
conductance 1/r between its ends: `elements.R`'s G and i stamps), 'C'
(each a capacitance: `elements.C`'s C and q) or 'I' (each a pulse current
source: `elements.IPulse`'s u, dudt and next_event) -- where an element
object costs ~5 KB of Python objects (the IBM power grids hold millions:
ibmpg3t's import alone filled 11 GB).

A circuit keeps its banks BESIDE its elements (`SubCircuit.add_bank`), not
among them: the walks over a circuit's elements never meet one.  What
assembles a matrix or a vector takes a bank's stamps from its arrays --
the stamp plan's constant groups, after every element
(`_stamp_plan._MatrixPlan`, `_VectorPlan`), the legacy loops
(`SubCircuit._add_element_submatrices` / `_subvectors`,
`_stamp_plan._legacy_*`) and the Jacobian's structure
(`analysis.jacobian_structure`); a view of the circuit carries them
(`_flat_view_of`, `pcnr.pcnr_view`, a copy).  Each element's stamp is the
element's own, `[[g, -g], [-g, g]]` with `g = 1.0 / r` (C: `c`): the same
values as `R` / `C` make them, summed after the elements rather than among
them (a run's rounding moves where elements and a bank would interleave).

A bank is constant: no parameters to write, no initial condition, no
element of its own to read results by, and no noise -- `CY` refuses a
circuit holding a resistor bank (`elements.R` is noisy by default).

A CURRENT-SOURCE BANK ('I') holds `IPulse` sources: per source its DC
value `i` (`value`) and its waveform's seven values (`wave`: i1, i2, td,
tr, tf, pw, per; tr and tf floored at `Pulse.MIN_EDGE`, as `Pulse.__init__`
floors them), and no AC value or noise -- an AC or noise analysis sees
nothing of it.  In the time domain ('dc', 'tran') each source puts `i +
f(t)` on its first node and its negative on its second, after every
element, each source's two entries added in order (`source_entries`) --
the direct source pass (`_stamp_plan._SourceDirect`) and its C table
(`_tran_sources`, the bank's rows after the elements') add them so too;
`dudt` its `dfdt(t)`; `next_event` the earliest over its sources.  Each
evaluated on the arrays as `func.Pulse` makes it for one source, operation
for operation -- the same bits (numpy's `remainder` and `floor_divide` are
Python's float `%` and `//`).
"""
import numpy as np

#: per kind, the matrix and the vector method a bank stamps (a current
#: source bank: its source methods, without x)
METHODS = {'R': ('G', 'i'), 'C': ('C', 'q'), 'I': ('u', 'dudt')}
#: the analyses in which a source has a value (`circuit.timedomain_analyses`)
_TIME_DOMAIN = ('dc', 'tran')
#: `Pulse.next_event`'s period index past which the Python's `int` is not
#: one the arrays make exactly (sources.c `PULSE_K_MAX`)
_K_MAX = 2.0 ** 52


class Bank:
    """`len(a)` elements of `kind` ('R', 'C' or 'I'): element e between
    `nodes[a[e]]` and `nodes[b[e]]`, its resistance, capacitance or DC
    current `value[e]` (a resistor's conductance its reciprocal) -- a
    current source's pulse `wave[e]` (i1, i2, td, tr, tf, pw, per)."""

    __slots__ = ('a', 'b', 'kind', 'name', 'nodes', 'value', 'wave')

    def __init__(self, name, kind, nodes, a, b, value, wave=None):
        if kind not in METHODS:
            raise ValueError(f'a bank of {kind!r}: R, C or I')
        self.name, self.kind, self.nodes = name, kind, list(nodes)
        self.a = np.ascontiguousarray(a, dtype=np.int64)
        self.b = np.ascontiguousarray(b, dtype=np.int64)
        self.value = np.ascontiguousarray(value, dtype=np.float64)
        n = len(self.nodes)
        if not (self.a.shape == self.b.shape == self.value.shape and self.a.ndim == 1):
            raise ValueError('a bank: a, b and value of one length')
        if self.a.size and (min(int(self.a.min()), int(self.b.min())) < 0
                            or max(int(self.a.max()), int(self.b.max())) >= n):
            raise ValueError('a bank: an end outside its nodes')
        self.wave = None
        if kind == 'I':
            from pycircuit.circuit.func import Pulse
            w = np.array(wave, dtype=np.float64).reshape(-1, 7) if wave is not None else None
            if w is None or w.shape[0] != self.a.size:
                raise ValueError('a current source bank: a pulse (i1, i2, td, tr, tf, pw, '
                                 'per) for each source')
            w[:, 3] = np.maximum(w[:, 3], Pulse.MIN_EDGE)
            w[:, 4] = np.maximum(w[:, 4], Pulse.MIN_EDGE)
            self.wave = np.ascontiguousarray(w)
        elif wave is not None:
            raise ValueError(f'a bank of {kind!r} has no waveform')

    def __len__(self):
        return int(self.a.size)

    def __repr__(self):
        return f'Bank({self.name!r}, {self.kind!r}, {len(self)} elements)'

    def stamps(self):
        """Every element's 2 x 2 stamp, `(E, 2, 2)`: `[[g, -g], [-g, g]]`,
        `g` the conductance `1.0 / r` (`R.update`'s `1/r`) or `c`."""
        g = 1.0 / self.value if self.kind == 'R' else self.value
        S = np.empty((g.size, 2, 2))
        S[:, 0, 0] = g
        S[:, 0, 1] = -g
        S[:, 1, 0] = -g
        S[:, 1, 1] = g
        return S


def _pulse_fold(W, t):
    """Each source's time in its period (`Pulse.f`'s `t % per` where `per`
    is not 0) and its edges `td`, `td + tr`, `td + tr + pw`, `td + tr + pw +
    tf`, summed left to right as the Python sums them."""
    per = W[:, 6]
    on = per != 0
    with np.errstate(all='ignore'):
        tt = np.where(on, np.remainder(t, np.where(on, per, 1.0)), t)
    e1 = W[:, 2] + W[:, 3]
    e2 = e1 + W[:, 5]
    return tt, e1, e2, e2 + W[:, 4]


def pulse_values(W, t):
    """`Pulse.f(t)` of each source, its `where` chain over the arrays."""
    v1, v2, td, tr, tf = W[:, 0], W[:, 1], W[:, 2], W[:, 3], W[:, 4]
    tt, e1, e2, e3 = _pulse_fold(W, t)
    with np.errstate(all='ignore'):
        v = np.where(tt < e3, v2 + (v1 - v2) / tf * (tt - e2), v1)
        v = np.where(tt < e2, v2, v)
        v = np.where(tt < e1, v1 + ((v2 - v1) / tr) * (tt - td), v)
        return np.where(tt < td, v1, v)


def pulse_slopes(W, t):
    """`Pulse.dfdt(t)` of each source: the segment entered at `t`."""
    v1, v2, td, tr, tf = W[:, 0], W[:, 1], W[:, 2], W[:, 3], W[:, 4]
    tt, e1, e2, e3 = _pulse_fold(W, t)
    with np.errstate(all='ignore'):
        rise, fall, zero = (v2 - v1) / tr, (v1 - v2) / tf, 0.0 * tt
        d = np.where(tt < e3, fall, zero)
        d = np.where(tt < e2, zero, d)
        d = np.where(tt < e1, rise, d)
        return np.where(tt < td, zero, d)


def pulse_next(W, t):
    """The earliest of the sources' `Pulse.next_event(t)` -- each one's
    candidates `j * per + edge` (j its period index by Python's `//` and the
    periods on either side; the four edges, then `per`), the first strictly
    after `t`; a source without a period its four edges -- or inf."""
    td, per = W[:, 2], W[:, 6]
    _tt, e1, e2, e3 = _pulse_fold(W, t)
    edges = np.stack([td, e1, e2, e3], axis=1)
    on = per != 0
    best = np.inf
    if (~on).any():
        c = edges[~on]
        later = c[c > t]
        if later.size:
            best = float(later.min())
    if on.any():
        p, E = per[on], np.concatenate([edges[on], per[on][:, None]], axis=1)
        with np.errstate(all='ignore'):
            k = np.floor_divide(t, p)
        if not np.all(np.abs(k) < _K_MAX):
            ## (a period index the Python's `int` makes otherwise, or raises
            ## on: each source as the Python asks it)
            from pycircuit.circuit.func import Pulse
            for row in W[on]:
                best = min(best, Pulse(*map(float, row)).next_event(t))
            return best
        for j in (-1.0, 0.0, 1.0):
            c = ((k + j) * p)[:, None] + E
            later = c[c > t]
            if later.size:
                best = min(best, float(later.min()))
    return best


def source_entries(bk, ia, ib, t, analysis, slope=False):
    """A current source bank's entries in `u` (`slope`: `dudt`) at `t`, its
    sources in order, each one's first node then its second: `(rows,
    values)`, the value `i + f(t)` and its negative (`dudt`: `dfdt(t)`, no
    DC value) -- or None outside the time domain, where a source puts
    nothing."""
    if analysis not in _TIME_DOMAIN:
        return None
    v = pulse_slopes(bk.wave, t) if slope else bk.value + pulse_values(bk.wave, t)
    return np.stack([ia, ib], axis=1).reshape(-1), np.stack([v, -v], axis=1).reshape(-1)


def next_event(cir, t):
    """The earliest breakpoint after `t` of `cir`'s current source banks,
    or inf."""
    best = np.inf
    for bk in banks(cir):
        if bk.kind == 'I' and len(bk):
            best = min(best, pulse_next(bk.wave, t))
    return best


def current_scale(cir):
    """The largest swing `cir`'s current source banks declare,
    `IPulse.signal_scale`'s max(|i1|, |i2|, |i2 - i1|) over their sources."""
    out = 0.0
    for bk in banks(cir):
        if bk.kind == 'I' and len(bk):
            i1, i2 = bk.wave[:, 0], bk.wave[:, 1]
            out = max(out, float(np.max(np.maximum(np.maximum(np.abs(i1), np.abs(i2)),
                                                   np.abs(i2 - i1)))))
    return out


def banks(cir):
    """The banks `cir` holds, in their order (empty where none)."""
    return cir.__dict__.get('_banks') or ()


def resolve(cir, node_index):
    """Each bank's ends as rows of `cir`'s solution vector, by
    `node_index` (`update_node_map`'s: node -> first place) -- kept on the
    circuit as `_bank_index` {name: (ia, ib)} with the node map."""
    out = {}
    for bk in banks(cir):
        idx = np.fromiter((node_index[nd] for nd in bk.nodes), dtype=np.intp,
                          count=len(bk.nodes))
        out[bk.name] = (idx[bk.a], idx[bk.b])
    cir._bank_index = out


def of(cir, m):
    """The banks stamping method `m` (a matrix method 'G' / 'C' or a vector
    'i' / 'q') with their ends' rows: `[(bank, ia, ib)]`.  `CY` (noise)
    refuses a circuit holding a resistor bank; any other method finds none."""
    found = banks(cir)
    if not found:
        return ()
    if m == 'CY' and any(bk.kind == 'R' for bk in found):
        raise NotImplementedError('noise: a resistor bank has no noise model (bulk elements)')
    cir._ensure_node_map()
    ix = cir.__dict__.get('_bank_index') or {}
    return [(bk, *ix[bk.name]) for bk in found if m in METHODS[bk.kind]]


def matrix_entries(bk, ia, ib):
    """A bank's matrix entries as an element's are scattered: per element
    its block's rows `(a, a, b, b)` and columns `(a, b, a, b)` (`np.repeat`
    and `np.tile` of its node map) with its stamp's values, flattened in
    element order."""
    rows = np.stack([ia, ia, ib, ib], axis=1).reshape(-1)
    cols = np.stack([ia, ib, ia, ib], axis=1).reshape(-1)
    return rows, cols, bk.stamps().reshape(-1)


def vector_entries(bk, ia, ib, x):
    """A bank's vector entries at `x`: each element's `S @ x[(a, b)]` (the
    plan's batched product of a constant stamp), in element order with its
    rows `(a, b)`."""
    I = np.stack([ia, ib], axis=1)
    with np.errstate(all='ignore'):
        Y = np.matmul(bk.stamps(), np.asarray(x, dtype=float)[I][:, :, None])[:, :, 0]
    return I.reshape(-1), Y.reshape(-1)
