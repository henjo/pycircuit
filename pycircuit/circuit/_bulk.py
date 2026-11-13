"""THE BULK LINEAR ELEMENTS (plan 261 L21): many two-terminal linear
elements of one kind held as arrays in one `Bank` -- 'R' (each a
conductance 1/r between its ends: `elements.R`'s G and i stamps) or 'C'
(each a capacitance: `elements.C`'s C and q) -- where an element object
costs ~5 KB of Python objects (the IBM power grids hold millions:
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
"""
import numpy as np

#: per kind, the matrix and the vector method a bank stamps
METHODS = {'R': ('G', 'i'), 'C': ('C', 'q')}


class Bank:
    """`len(a)` elements of `kind` ('R' or 'C'): element e between
    `nodes[a[e]]` and `nodes[b[e]]`, its resistance or capacitance
    `value[e]` (a resistor's conductance its reciprocal)."""

    __slots__ = ('a', 'b', 'kind', 'name', 'nodes', 'value')

    def __init__(self, name, kind, nodes, a, b, value):
        if kind not in METHODS:
            raise ValueError(f'a bank of {kind!r}: R or C')
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
