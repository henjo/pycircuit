"""KLU's order tried against the structural pattern's always-zero places
(`_klu_order`, `linearsolver.StructuralKLU._analysed`): two resistor meshes
joined at every node by a 0 V via -- the IBM power grids' form, where KLU's
own order pivots on the vias' zero branch diagonals and fills past its
estimate -- take the matched order, which the trial passes, and the run is
the same with the C paths as without; a circuit whose trial passes keeps
KLU's own analysis; the trial's values are a circuit matrix's dominance (a
uniform draw misreads a MOS chain); the stamped places are the plan's."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _dc_newton_c,
    _klu_order,
    _paths,
    _tran_branch_c,
    _tran_loop_c,
    _tran_lte_c,
    _tran_newton_c,
    circuit,
)
from pycircuit.circuit.analysis import jacobian_structure, stamped_places
from pycircuit.circuit.elements import IS, VS, C, R, SubCircuit, gnd
from pycircuit.circuit.linearsolver import StructuralKLU
from pycircuit.circuit.tests.test_first_attempt_c import _chain
from pycircuit.circuit.transient import Transient


def _offsets():
    off = _klu_order._offsets()
    if off is None:
        pytest.skip('no KLU header or compiler here: no trial is made')
    return off


def _layers(k=8):
    """Two k x k resistor meshes (a, b) joined at every node by a 0 V via,
    sinks and capacitors on a, pads to the supply on b."""
    c = SubCircuit()
    c['vdd'] = VS('s', gnd, v=1.8)
    for L, r in (('a', 0.5), ('b', 0.1)):
        for i in range(k):
            for j in range(k):
                if i + 1 < k:
                    c[f'r{L}h{i}_{j}'] = R(f'{L}{i}_{j}', f'{L}{i + 1}_{j}', r=r)
                if j + 1 < k:
                    c[f'r{L}v{i}_{j}'] = R(f'{L}{i}_{j}', f'{L}{i}_{j + 1}', r=r)
    for i in range(k):
        for j in range(k):
            c[f'v{i}_{j}'] = VS(f'a{i}_{j}', f'b{i}_{j}', v=0.0)
            c[f'c{i}_{j}'] = C(f'a{i}_{j}', gnd, c=1e-14 * (1 + (i + j) % 3))
            c[f'i{i}_{j}'] = IS(f'a{i}_{j}', gnd, i=1e-4 * (1 + (i * j) % 5))
            if i % 4 == 0 and j % 4 == 0:
                c[f'rp{i}_{j}'] = R('s', f'b{i}_{j}', r=0.25)
    return c


def _structural(cir):
    """`(the StructuralKLU on cir's pattern, KLU's own symbolic, the stamped
    places)` -- the solver made as the analyses make it, unanalysed."""
    iref = cir.get_node_index(circuit.gnd)
    m, rows, cols = jacobian_structure(cir, iref)
    sk = StructuralKLU(m, rows, cols, fallback=None)
    sym = sk._analysed()
    mask = stamped_places(cir, iref, sk.Ai.astype(np.intp), np.repeat(np.arange(m), np.diff(sk.Ap)))
    return sk, sym, mask


def _run(on, monkeypatch):
    for mod in (_tran_newton_c, _tran_loop_c, _tran_lte_c, _dc_newton_c, _tran_branch_c):
        monkeypatch.setattr(mod, 'ENABLED', on)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_layers(), toolkit=circuit.numeric)
        res = tr.solve(tend=4e-10, timestep=1e-11)
    st = tr.statistics
    run = (np.asarray(res.x, float).tobytes(), np.asarray(res.sweep_values, float).tobytes(),
           st.accepted_steps, st.rejected_steps, st.newton_iterations)
    return run, tr._get_linearsolver()._choice, _paths.since(before)


def test_meshes_joined_by_zero_volt_vias_take_the_matched_order(monkeypatch):
    """KLU's own order on two meshes joined by 0 V vias fails the trial (it
    pivots on the vias' zero diagonals: reallocation, or fill past twice its
    estimate); the run's solver keeps the matched order, which passes it --
    counted, the operating point's solver too -- and the run is the same
    with the C paths as without (they share the one analysis)."""
    off = _offsets()
    py, sk_py, _ = _run(False, monkeypatch)
    got, sk, counts = _run(True, monkeypatch)
    assert got == py
    assert counts.get('klu_order:reordered', 0) >= 2, counts
    assert counts.get('newton_c:served', 0) > 0, counts
    own, own_sym, mask = _structural(_layers())
    assert _klu_order._fails(own, own_sym, off, mask) is True
    fill_own = _klu_order._trial(own, own_sym, off, mask)[1]
    for s in (sk, sk_py):
        assert type(s) is StructuralKLU
        assert _klu_order._trial(s, s._symbolic, off, mask)[1] < fill_own


def test_a_circuit_whose_trial_passes_keeps_klus_own_order(monkeypatch):
    """A 400-stage MOS chain: its branch rows under `BRANCH_SHARE` of its
    unknowns, no trial is made (no plan read: `stamped_places` None); made
    all the same (the share 0), it passes, and the solver's symbolic is
    KLU's own -- its P and Q a plain analysis's, so nothing it factors
    moves."""
    off = _offsets()
    cir = _chain(400)
    iref = cir.get_node_index(circuit.gnd)
    m, rows, cols = jacobian_structure(cir, iref)
    assert stamped_places(cir, iref, rows, cols) is None
    monkeypatch.setattr(_klu_order, 'BRANCH_SHARE', 0.0)
    before = _paths.snapshot()
    sk = StructuralKLU(m, rows, cols, fallback=None,
                       stamped=lambda r, c: stamped_places(cir, iref, r, c))
    sym = sk._analysed()
    counts = _paths.since(before)
    plain = StructuralKLU(m, rows, cols, fallback=None)
    psym = plain._analysed()
    assert counts.get('klu_order:trial') == 1 and 'klu_order:reordered' not in counts, counts
    for key in ('s_P', 's_Q'):
        assert np.array_equal(_klu_order._perm(sym, off, key, m), _klu_order._perm(psym, off, key, m))


def test_the_trials_values_are_a_circuit_matrixs_dominance(monkeypatch):
    """The trial's diagonals dominate their rows where an element stamps
    them: a 400-stage MOS chain passes it, and a uniform draw at the same
    places -- no dominance -- misreads it as an order built on zeros."""
    off = _offsets()
    monkeypatch.setattr(_klu_order, 'BRANCH_SHARE', 0.0)
    sk, sym, mask = _structural(_chain(400))
    assert _klu_order._fails(sk, sym, off, mask) is False
    monkeypatch.setattr(_klu_order, 'trial_values', lambda m, rows, cols, st: np.where(
        st, np.random.default_rng(_klu_order.SEED).uniform(-1.0, 1.0, len(rows)), 0.0))
    assert _klu_order._fails(sk, sym, off, mask) is True


def test_the_stamped_places_are_the_plans():
    """A 0 V via's branch diagonal is a place of the pattern no element
    stamps (the ladders' diagonal), its four couplings and a resistor's
    block are stamped -- and every stamped place lies in the pattern."""
    cir = _layers(2)
    iref = cir.get_node_index(circuit.gnd)
    _m, rows, cols = jacobian_structure(cir, iref)
    mask = stamped_places(cir, iref, rows, cols)
    idx = {nm: cir.get_node_index(nm) for nm in ('a0_0', 'b0_0', 'a0_1')}
    red = lambda i: i - (i > iref)
    br = next(x for x in cir.branches if x.plus.name == 'a0_0' and x.minus.name == 'b0_0')
    b = red(cir.get_branch_index(br))
    at = {(r, c): k for k, (r, c) in enumerate(zip(rows.tolist(), cols.tolist()))}
    assert not mask[at[(b, b)]]
    assert mask[at[(b, red(idx['a0_0']))]] and mask[at[(red(idx['b0_0']), b)]]
    assert mask[at[(red(idx['a0_0']), red(idx['a0_1']))]]
    assert mask.sum() < len(mask)


def test_a_matched_order_that_fills_no_less_is_not_kept(monkeypatch):
    """The matched order is kept only where its own trial fills less than
    KLU's: one that fills the same (KLU's own analysis again, standing in
    for it) is freed and KLU's order stands, counted."""
    _offsets()
    own, own_sym, mask = _structural(_layers())

    def same(sk, off, stamped):
        ct = sk._ct
        return sk._lib.klu_analyze(sk.m, sk._ptr(sk.Ap, ct.c_int), sk._ptr(sk.Ai, ct.c_int),
                                   ct.byref(sk._common))
    monkeypatch.setattr(_klu_order, '_matched', same)
    before = _paths.snapshot()
    assert _klu_order.better(own, own_sym, mask) is None
    assert _paths.since(before).get('klu_order:kept_own') == 1


def _mesh(k):
    """A k x k resistor mesh with a capacitor and a sink at every node and a
    pad to the supply at each corner: a 2-D grid, its one branch row the
    supply's (under `BRANCH_SHARE`)."""
    c = SubCircuit()
    c['vdd'] = VS('s', gnd, v=1.8)
    for i in range(k):
        for j in range(k):
            if i + 1 < k:
                c[f'h{i}_{j}'] = R(f'n{i}_{j}', f'n{i + 1}_{j}', r=0.5)
            if j + 1 < k:
                c[f'v{i}_{j}'] = R(f'n{i}_{j}', f'n{i}_{j + 1}', r=0.5)
            c[f'c{i}_{j}'] = C(f'n{i}_{j}', gnd, c=1e-14)
            c[f'i{i}_{j}'] = IS(f'n{i}_{j}', gnd, i=1e-4 * (1 + (i * j) % 5))
    for i, j in ((0, 0), (0, k - 1), (k - 1, 0), (k - 1, k - 1)):
        c[f'p{i}_{j}'] = R('s', f'n{i}_{j}', r=0.25)
    return c


def _analysed(cir):
    """`(the solver analysed as the analyses make it -- the stamped places
    given --, its symbolic, its size, the stamped places, the counts)`."""
    iref = cir.get_node_index(circuit.gnd)
    m, rows, cols = jacobian_structure(cir, iref)
    before = _paths.snapshot()
    sk = StructuralKLU(m, rows, cols, fallback=None,
                       stamped=lambda r, c: stamped_places(cir, iref, r, c))
    sym = sk._analysed()
    counts = _paths.since(before)
    mask = stamped_places(cir, iref, sk.Ai.astype(np.intp), np.repeat(np.arange(m), np.diff(sk.Ap)))
    return sk, sym, m, mask, counts


@pytest.mark.parametrize('build, kept', [(lambda: _mesh(70), 'dissected'), (lambda: _mesh(40), 'own'),
                                         (lambda: _layers(30), 'reordered')],
                         ids=['mesh70', 'mesh40', 'layers30'])
def test_from_its_size_a_dissection_is_tried_and_the_least_fill_kept(build, kept, monkeypatch):
    """From `ND_FROM` unknowns (lowered here) a nested dissection of A + A' is
    tried -- the stamped places made whatever the branch share -- and the
    order kept is the one whose trial fills least: on a 70 x 70 mesh the
    dissection (AMD fills past it), on a 40 x 40 one KLU's own, on layered
    meshes joined by 0 V vias the matched order (it beats the dissection
    there, as on the IBM grid ibmpg1t); the same order every time."""
    off = _offsets()
    if _klu_order._nd_fn() is None:
        pytest.skip('no CHOLMOD header: no dissection tried')
    monkeypatch.setattr(_klu_order, 'ND_FROM', 1000)
    sk, sym, m, mask, counts = _analysed(build())
    assert mask is not None and counts.get('klu_order:dissection') == 1, counts
    plain = StructuralKLU(m, sk.Ai, np.repeat(np.arange(m), np.diff(sk.Ap)), fallback=None)
    psym = plain._analysed()
    nd, mt = _klu_order.dissected(plain), _klu_order._matched(plain, off, mask)
    fills = {'own': _klu_order._trial(plain, psym, off, mask)[1],
             'dissected': _klu_order._trial(plain, nd, off, mask)[1]}
    if mt is not None:
        fills['reordered'] = _klu_order._trial(plain, mt, off, mask)[1]
    for extra in (nd, mt):
        if extra is not None:
            _klu_order._free(plain, extra)
    assert _klu_order._trial(sk, sym, off, mask)[1] == fills[kept]
    assert (kept == 'dissected') is ('klu_order:dissected' in counts), counts
    assert (kept == 'reordered') is ('klu_order:reordered' in counts), counts
    ## (the second solver kept alive: its symbolic is freed with it)
    sk2, again = _analysed(build())[:2]
    for key in ('s_P', 's_Q'):
        assert np.array_equal(_klu_order._perm(sym, off, key, m), _klu_order._perm(again, off, key, m))
    del sk2


def test_below_its_size_no_dissection_is_tried():
    """Below `ND_FROM` the 70 x 70 mesh -- its branch rows under
    `BRANCH_SHARE` -- makes no trial at all (no stamped places read), and
    its symbolic is KLU's own."""
    off = _offsets()
    sk, sym, m, mask, counts = _analysed(_mesh(70))
    assert mask is None and not any(k.startswith('klu_order:') for k in counts), counts
    plain = StructuralKLU(m, sk.Ai, np.repeat(np.arange(m), np.diff(sk.Ap)), fallback=None)
    psym = plain._analysed()
    for key in ('s_P', 's_Q'):
        assert np.array_equal(_klu_order._perm(sym, off, key, m), _klu_order._perm(psym, off, key, m))


def test_a_run_on_the_dissection_is_the_same_with_the_c(monkeypatch):
    """The 70 x 70 mesh's transient on the dissection (`ND_FROM` lowered):
    the C paths and the Python share the one analysis -- the same run with
    them as without, bit for bit -- and the dissection kept."""
    if _klu_order._nd_fn() is None:
        pytest.skip('no CHOLMOD header: no dissection tried')
    monkeypatch.setattr(_klu_order, 'ND_FROM', 1000)

    def run(on):
        for mod in (_tran_newton_c, _tran_loop_c, _tran_lte_c, _dc_newton_c, _tran_branch_c):
            monkeypatch.setattr(mod, 'ENABLED', on)
        before = _paths.snapshot()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            tr = Transient(_mesh(70), toolkit=circuit.numeric)
            res = tr.solve(tend=1e-10, timestep=1e-11)
        st = tr.statistics
        return ((np.asarray(res.x, float).tobytes(), st.accepted_steps, st.rejected_steps,
                 st.newton_iterations), _paths.since(before))
    py, _cpy = run(False)
    got, counts = run(True)
    assert got == py
    assert counts.get('klu_order:dissected', 0) >= 1 and counts.get('newton_c:served', 0) > 0, counts
