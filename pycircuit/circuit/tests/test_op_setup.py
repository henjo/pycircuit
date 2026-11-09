"""THE OPERATING POINT'S SETUP AT SIZE (plan 261's L8): the node map's 2-D
indices as `np.repeat` and `np.tile` make them (`Circuit.update_node_map`,
where a `meshgrid` an element was most of the map), the parameters read in
bulk (`ParameterDict.values_of`: `pcnr`'s device dicts, the C kernels' packs),
and the PCNR's first reduced Schur matrix on pattern pages
(`_pcnr_c._reduced_J`) -- each the same values as before."""
import warnings

import numpy as np

from pycircuit.circuit import _pcnr_c, _tran_core, circuit
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.tests.test_radau_pattern_pages import _chain
from pycircuit.circuit.transient import Transient
from pycircuit.utilities.param import Parameter, ParameterDict


def test_the_node_maps_indices_are_meshgrids():
    """Each element's 2-D indices are `np.meshgrid(nm, nm, indexing='ij')`
    flattened: values, dtype, order -- on a chain of resistors, capacitors,
    MOSFETs, junctions and sources (branches among them)."""
    cir = _chain(30)
    cir.update_node_map()
    for inst, nm in cir.elementnodemap.items():
        r, c = np.meshgrid(nm, nm, indexing='ij')
        got = cir._map_indices_2d[inst]
        for g, w in zip(got, (r.flatten(), c.flatten())):
            assert g.dtype == w.dtype and g.shape == w.shape and g.tobytes() == w.tobytes()
            assert g.flags.c_contiguous
        assert not np.shares_memory(got[0], got[1])


def test_values_of_is_getattr_for_every_name():
    """`values_of` reads each name as `getattr` does: a parameter's value, a
    name shadowed by an instance attribute or a class attribute (`append`,
    a method), an unknown name's AttributeError."""
    pd = ParameterDict(Parameter('r', default=2.0), Parameter('append', default=3.0),
                       Parameter('w', default=4.0))
    pd.r = 5.0
    pd.__dict__['w'] = 'shadow'
    names = ['r', 'append', 'w']
    got = pd.values_of(names)
    want = [getattr(pd, n) for n in names]
    assert got[0] == want[0] == 5.0 and got[2] == want[2] == 'shadow'
    assert got[1] == want[1] and callable(got[1])
    try:
        pd.values_of(['nope'])
    except AttributeError:
        pass
    else:
        raise AssertionError('an unknown name reads as getattr does: AttributeError')


def test_the_reduced_schur_matrix_on_pattern_pages_is_the_dense_one(monkeypatch):
    """At size the PCNR's first reduced Schur matrix is on pattern pages,
    its bytes `J[np.ix_(keep, keep)]`'s (J zero off the core's pattern)."""
    monkeypatch.setattr(_tran_core, 'PATTERN_PAGES_FROM', 1)
    cir = _chain(300)
    seen = {}
    release = _pcnr_c.release

    def keep(c):
        rec = c.__dict__['_pcnr_c']
        k = rec.red >= 0
        seen['dense'] = _pcnr_c._dense_J(rec)[np.ix_(k, k)].copy()
        seen['paged'] = _pcnr_c._reduced_J(rec)
        return release(c)

    monkeypatch.setattr(_pcnr_c, 'release', keep)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        Transient(cir, toolkit=circuit.numeric, integrator=RadauIIA3Integrator()).solve(
            tend=2e-9, timestep=1e-9)
    if 'dense' not in seen:
        import pytest
        pytest.skip('the op is not the PCNR C record\'s')
    d, p = seen['dense'], seen['paged']
    assert np.count_nonzero(d) > 0
    assert d.shape == p.shape and d.tobytes() == np.asarray(p).tobytes()
    assert _tran_core.pattern_of(p) is not None


def _ref_structure(cir, iref):
    """The per-element build of `jacobian_structure`."""
    n = cir.n
    flats = [np.arange(n, dtype=np.intp) * (n + 1)]
    for rows, cols in cir._map_indices_2d.values():
        flats.append(np.asarray(rows, dtype=np.intp) * n + np.asarray(cols, dtype=np.intp))
    r, c = np.divmod(np.unique(np.concatenate(flats)), n)
    keep = (r != iref) & (c != iref)
    r, c = r[keep] - (r[keep] > iref), c[keep] - (c[keep] > iref)
    order = np.lexsort((r, c))
    return n - 1, r[order], c[order]


def test_the_jacobian_structure_is_the_per_element_one():
    """`jacobian_structure` -- every element's block concatenated at once --
    is the per-element build's: size, rows, columns, dtype, CSC order."""
    from pycircuit.circuit.analysis import jacobian_structure
    cir = _chain(40)
    cir.update_node_map()
    for iref in (0, 3, cir.n - 1):
        got, want = jacobian_structure(cir, iref), _ref_structure(cir, iref)
        assert got[0] == want[0]
        assert all(g.dtype == w.dtype and g.tobytes() == w.tobytes() for g, w in zip(got[1:], want[1:]))


def test_the_plans_nonzeros_follow_the_plan_they_read(monkeypatch):
    """`_plan_nonzeros` keeps its last plan's sorted positions: a C of
    another circuit, read in turn, gets its own plan's -- the result each
    time the unkept one's."""
    from pycircuit.circuit import _tran_branch
    from pycircuit.circuit.tests.test_radau_pattern_pages import _chain as chain
    out = []
    for N in (12, 30, 12):
        cir = chain(N)
        cir.toolkit = circuit.numeric
        x = np.linspace(-0.3, 0.9, cir.n)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            C = np.asarray(cir.C(x), dtype=float)
        got = _tran_branch._plan_nonzeros(cir, C)
        monkeypatch.setattr(_tran_branch, '_UPOS', [None])
        want = _tran_branch._plan_nonzeros(cir, C)
        assert got is not None and all(np.array_equal(g, w) for g, w in zip(got, want))
        out.append(len(got[2]))
    assert out[0] == out[2] != out[1]
