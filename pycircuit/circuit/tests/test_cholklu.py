"""THE SPD REDUCTION ON CHOLMOD (`_cholklu`): a linear grid's transient
matrix solved through KLU's entries over CHOLMOD's supernodal Cholesky --
the voltage sources' nodes fixed, the inductors' currents eliminated, the
free nodes' symmetric positive definite system factored -- as a dense solve
makes it; the plan declining what is not a grid; a matrix CHOLMOD will not
factor KLU's; a transient on it its KLU run, its C the Python's; the kept
factors, the step ladder, and the glue's own solve on CHOLMOD's factor."""
import ctypes
import hashlib
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _cholklu, _paths, linearsolver, numeric
from pycircuit.circuit.analysis import jacobian_structure
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import IS, VS, C, Diode, IPulse, R
from pycircuit.circuit.elements import L as Ind
from pycircuit.circuit.transient import Transient


def _grid(nx=6, ny=5):
    """An RC mesh fed by a supply through an inductor (an L from a fixed node
    to a free one), a second supply on the far corner, an inductor across
    the mesh and one to ground, a current source."""
    c = SubCircuit()
    nm = lambda i, j: f'n{i}_{j}'
    c['vdd'] = VS('s', gnd, v=1.0)
    c['lp'] = Ind('s', nm(0, 0), L=1e-10)
    c['v2'] = VS(nm(nx - 1, ny - 1), gnd, v=0.9)
    c['is'] = IS(nm(2, 2), gnd, i=1e-3)
    k = 0
    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                c[f'r{k}'] = R(nm(i, j), nm(i + 1, j), r=1.0 + 0.1 * (k % 7))
                k += 1
            if j + 1 < ny:
                c[f'r{k}'] = R(nm(i, j), nm(i, j + 1), r=2.0 + 0.1 * (k % 5))
                k += 1
            c[f'c{i}_{j}'] = C(nm(i, j), gnd, c=1e-12 * (1 + (i + j) % 3))
    c['lx'] = Ind(nm(1, 1), nm(3, 3), L=2e-10)
    c['lg'] = Ind(nm(4, 4), gnd, L=3e-10)
    ## (a supply made last, its inductor from an early free node: the free end
    ## the lower row, as the first supply's inductor has it the higher)
    c['v3'] = VS('late', gnd, v=0.8)
    c['ll'] = Ind(nm(1, 0), 'late', L=1.5e-10)
    return c


def _matrix(c, h):
    iref = c.get_node_index(gnd)
    x = np.zeros(c.n)
    J = np.asarray(c.G(x)) + np.asarray(c.C(x)) / h
    return np.delete(np.delete(J, iref, 0), iref, 1), iref


def _solver(c, plan=True):
    iref = c.get_node_index(gnd)
    m, rows, cols = jacobian_structure(c, iref)
    return linearsolver.StructuralKLU(m, rows, cols, fallback=linearsolver.SuperLUSolver,
                                      cholmod=_cholklu.plan(c, iref) if plan else None)


@pytest.fixture(autouse=True)
def _cholmod_here():
    if _cholklu.library(linearsolver.KLUSolver()._lib) is None:
        pytest.skip(f'no CHOLMOD here: {_cholklu.STATUS}')


def test_the_reduction_solves_as_a_dense_solve():
    """On the reduction: a factor, a refactor (h doubled), several
    right-hand sides, the transposed solve and an owned factor -- each the
    dense solve's to rounding, its backward error at the arithmetic's."""
    c = _grid()
    sk = _solver(c)
    assert sk.cholmod
    rng = np.random.default_rng(0)
    for h in (1e-12, 2e-12):
        J, _iref = _matrix(c, h)
        b = rng.standard_normal(J.shape[0])
        x = np.asarray(sk.solve(J, b, None))
        assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-13 * np.max(np.abs(x))
        r = J @ x - b
        assert np.max(np.abs(r)) <= 1e-15 * np.max(np.abs(J) @ np.abs(x) + np.abs(b))
        B = rng.standard_normal((J.shape[0], 3))
        X = np.asarray(sk.solve(J, B, None))
        assert np.max(np.abs(X - np.linalg.solve(J, B))) <= 1e-13 * np.max(np.abs(X))
        XT = np.asarray(sk.solve_transposed(J, B, None))
        assert np.max(np.abs(XT - np.linalg.solve(J.T, B))) <= 1e-13 * np.max(np.abs(XT))
        f = sk.factor(J, None)
        assert np.max(np.abs(np.asarray(f.solve(b)) - x)) <= 1e-13 * np.max(np.abs(x))
    assert sk.factors == 1 and sk.refactors >= 1 and sk.cholmod


def test_a_refactor_check_failing_on_the_reduction_is_kept():
    """The refactor's residual check refused (a tolerance of 0): the
    reduction's growth entry reports a sound factor (a Cholesky of a positive
    definite matrix has no pivot order to go stale), so the refactor is
    kept (THE REFACTOR KEPT ON ITS GROWTH) -- no fresh factor, counted."""
    c = _grid()
    sk = _solver(c)
    sk.REFACTOR_RESIDUAL_TOL = 0.0
    b = np.random.default_rng(3).standard_normal(sk.m)
    for h in (1e-12, 2e-12):
        J, _iref = _matrix(c, h)
        x = np.asarray(sk.solve(J, b, None))
        assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-13 * np.max(np.abs(x))
    assert (sk.factors, sk.refactors, sk.residual_fallbacks, sk.residual_kept) == (1, 1, 0, 1)


def test_the_plan_declines_what_is_not_a_grid():
    """A device, a voltage source between two nodes, two on one node, a
    subcircuit instance, an instance's own stamp: each named; the grid's
    rows counted by kind."""
    c = _grid()
    got = _cholklu.plan(c, c.get_node_index(gnd))
    assert list(np.bincount(got)) == [c.n - 1 - 7, 3, 4]
    cases = {}
    d = _grid()
    d['d1'] = Diode('n1_2', gnd)
    cases['element'] = d
    f = _grid()
    f['vf'] = VS('n1_2', 'n2_2', v=0.1)
    cases['vsource'] = f
    t = _grid()
    t['vt'] = VS('s', gnd, v=0.1)
    cases['vsource '] = t

    class Stage(SubCircuit):
        terminals = ('a', 'b')

        def __init__(self, *args, **kw):
            super().__init__(*args, **kw)
            self['R1'] = R('a', 'b', r=1.0)
    s = _grid()
    s['X1'] = Stage('n0_1', 'n1_1')
    cases['hierarchy'] = s
    o = _grid()
    G0 = o['r0'].G
    o['r0'].G = lambda x, *a: G0(x, *a)
    cases['instance'] = o
    for why, cir in cases.items():
        assert _cholklu.plan(cir, cir.get_node_index(gnd)) == why.strip(), why


def test_a_matrix_cholmod_refuses_is_klus(monkeypatch):
    """The inductors' diagonals zero (h infinite: the DC's matrix) -- CHOLMOD
    refuses the factor -- and a pattern a voltage source's row couples two
    nodes in: each solved by KLU, counted, the dense solve's answer."""
    c = _grid()
    sk = _solver(c)
    J, _iref = _matrix(c, np.inf)
    b = np.random.default_rng(1).standard_normal(J.shape[0])
    before = _paths.snapshot()
    x = np.asarray(sk.solve(J, b, None))
    assert _paths.since(before).get('solver.cholmod:klu') == 1
    assert np.max(np.abs(x - np.linalg.solve(J, b))) <= 1e-12 * np.max(np.abs(x))
    assert sk.cholmod and sk._ck_twin is not None and sk._ck_twin.factors == 1
    ## (a pattern not a grid's: KLU from the analysis, nothing made on CHOLMOD)
    iref = c.get_node_index(gnd)
    m, rows, cols = jacobian_structure(c, iref)
    kind = _cholklu.plan(c, iref).copy()
    kind[np.flatnonzero(kind == 0)[3]] = 1
    before = _paths.snapshot()
    bad = linearsolver.StructuralKLU(m, rows, cols, fallback=linearsolver.SuperLUSolver, cholmod=kind)
    J1, _iref = _matrix(c, 1e-12)
    x = np.asarray(bad.solve(J1, b, None))
    assert _paths.since(before).get('solver.cholmod:pattern') == 1 and not bad.cholmod
    assert np.max(np.abs(x - np.linalg.solve(J1, b))) <= 1e-12 * np.max(np.abs(x))


def test_the_threshold_the_switch_and_the_dc_keep_klu(monkeypatch):
    """Below `CHOLMOD_FROM`, with the switch off, and in the DC (its
    inductors shorts) or any analysis but the transient (`_SPD_PLAN`): KLU;
    from the threshold, a transient's matrix the reduction's
    (`solver.auto:cholmod`)."""
    def choice(make):
        c = _grid(10, 10)
        c['ip'] = IPulse('n2_3', 'n0_0', i1=0.0, i2=1e-3, td=1e-11, tr=1e-11, tf=1e-11, pw=1e-10,
                         per=4e-10)
        a = make(c)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            if isinstance(a, Transient):
                a.solve(tend=5e-11, timestep=1e-11)
            else:
                a.solve()
        return a._get_linearsolver()

    def tran(c):
        return Transient(c, toolkit=numeric)

    def dc(c):
        return DC(c, toolkit=numeric)
    assert not choice(tran)._choice.cholmod
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0)
    before = _paths.snapshot()
    assert choice(tran)._choice.cholmod and _paths.since(before).get('solver.auto:cholmod') == 1
    assert choice(dc).spd_plan is None
    from pycircuit.circuit.analysis_ss import AC
    assert AC(_grid(), toolkit=numeric)._get_linearsolver().spd_plan is None
    monkeypatch.setattr(_cholklu, 'ENABLED', False)
    assert not choice(tran)._choice.cholmod


def _run(chol, c_on, monkeypatch, edge=1e-10, **make):
    from pycircuit.circuit import _tran_loop_c, _tran_lte_c, _tran_newton_c
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0 if chol else 1 << 40)
    for mod in (_tran_loop_c, _tran_lte_c, _tran_newton_c):
        monkeypatch.setattr(mod, 'ENABLED', c_on)
    c = _grid(16, 16)
    c['ip'] = IPulse('n8_8', 'n0_0', i1=0.0, i2=5e-3, td=1e-10, tr=edge, tf=edge, pw=4e-10,
                     per=2e-9)
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(c, toolkit=numeric, **make)
        res = tr.solve(tend=3e-9, timestep=1e-11)
    st = tr.statistics
    x = np.asarray(res.x, float)
    sk = tr._get_linearsolver()._choice
    return (x, hashlib.sha1(x.tobytes()).hexdigest(), (st.accepted_steps, st.rejected_steps,
            st.newton_iterations), (sk.cholmod, sk.factors, sk.refactors, sk.residual_fallbacks,
            sk.chord_factors, sk.chord_refactors), _paths.since(before), sk.kept_counts(),
            sk.solve_counts())


@pytest.mark.parametrize('chord', [True, False], ids=['chord', 'full'])
def test_a_grid_transient_on_cholmod_is_its_klu_run_and_its_c_the_python(chord, monkeypatch):
    """A linear grid's Gear-2 run, the threshold at 0: on the reduction the C
    Newton, error test and stepping loop serve (their KLU entries the
    reduction's), the run the Python's byte for byte; against KLU's run the
    same steps, the answers within rounding."""
    ## (today's stepping on both: the ladder is the reduction's alone -- its own tests below)
    make = {'chord_jacobian': chord, 'step_ladder': False}
    k = _run(False, True, monkeypatch, **make)
    c = _run(True, True, monkeypatch, **make)
    p = _run(True, False, monkeypatch, **make)
    assert c[1:4] == p[1:4] and c[3][0], (c[2:4], p[2:4])
    assert c[2] == k[2] and not k[3][0]
    assert float(np.max(np.abs(c[0] - k[0]))) <= 1e-12
    d = c[4]
    assert d.get('solver.auto:cholmod') == 1 and d.get('loop_c:steps', 0) > 50, d
    assert d.get('newton_c:served', 0) > 50 and d.get('lte_c:sparse_served', 0) > 50, d
    assert not any(key.startswith(('newton_c:bail', 'lte_c:bail', 'solver.cholmod:')) for key in d), d


## THE KEPT FACTORS (`_cholklu`'s note, plan 261 L24)

#: a step-size sequence that revisits its matrices, as a ladder's corners do
_HS = tuple(h * 1e-12 for h in (1, 2, 4, 1, 2, 8, 1, 4, 2, 8, 16, 1, 2, 4))


def _replay(hs, keep):
    """The factorizations `keep` held factors make on the sequence of step
    sizes (each a matrix; the last one's repeat none): the least recently
    used kept, the first factor counted."""
    held, made = [], 0
    for h in hs:
        if held and held[-1] == h:
            continue
        if h in held:
            held.remove(h)
        else:
            made += 1
            if len(held) >= keep:
                held.pop(0)
        held.append(h)
    return made


def _kept_solves(keep, max_bytes=None, c=None):
    c = _grid() if c is None else c
    sk = _solver(c)
    sk.keep_factors(keep, max_bytes)
    b = np.random.default_rng(1).standard_normal(sk.m)
    xs = []
    for h in _HS:
        J, _iref = _matrix(c, h)
        xs.append(np.asarray(sk.solve(J, b, None)).tobytes())
    return xs, sk


def test_kept_factors_give_a_refactors_bytes():
    """A sequence revisiting its matrices: with 2, 3 and 8 factors kept
    every solution is the bytes of the run keeping none; the factorizations
    those the least recently used kept make, the rest served by a kept
    factor (no factorization: CHOLMOD's factor of the same values on the
    same symbolic is the same bytes)."""
    want, sk1 = _kept_solves(1)
    ## (the refactor calls: a step size other than the last one's)
    calls = sum(1 for i, h in enumerate(_HS) if i == 0 or h != _HS[i - 1])
    assert sk1.kept_counts() == {'hits': 0, 'factorizations': calls, 'slots': 0, 'held': 0}
    for keep in (2, 3, 8):
        got, sk = _kept_solves(keep)
        assert got == want, keep
        k = sk.kept_counts()
        assert k['factorizations'] == _replay(_HS, keep), (keep, k)
        assert k['hits'] + k['factorizations'] == calls, (keep, k)
        assert k['slots'] == min(keep - 1, len(set(_HS)) - 1), (keep, k)
    assert _replay(_HS, 8) < _replay(_HS, 3) < calls


def test_kept_factors_stay_within_their_bytes():
    """Factors kept only within the bytes allowed: below one slot (L's values
    and J's) none is kept, every refactor a factorization, the same bytes."""
    want, _sk1 = _kept_solves(1)
    got, sk = _kept_solves(8, max_bytes=1.0)
    assert got == want
    assert sk.kept_counts()['slots'] == 0 and sk.kept_counts()['hits'] == 0


def test_a_hand_back_is_restored_from_a_kept_factor():
    """A checkpoint opened on J held, a refactor of J' made inside it, then
    handed back (`linearsolver.held_back`): the restore -- a refactor of J on
    the backup values -- is served by the factor kept for J, no
    factorization; J solved next is the bytes it was."""
    from pycircuit.circuit.linearsolver import held_back, held_open
    c = _grid()
    sk = _solver(c)
    sk.keep_factors(3)
    J, _iref = _matrix(c, 1e-12)
    J1, _iref = _matrix(c, 2e-12)
    b = np.random.default_rng(2).standard_normal(sk.m)
    sk.solve(J, b, None)
    xa = np.asarray(sk.solve(J, b, None)).tobytes()
    own = held_open(sk)
    sk.solve(J1, b, None)
    before = sk.kept_counts()
    held_back(sk, own)
    after = sk.kept_counts()
    assert after['hits'] == before['hits'] + 1, (before, after)
    assert after['factorizations'] == before['factorizations'], (before, after)
    assert np.asarray(sk.solve(J, b, None)).tobytes() == xa


@pytest.mark.parametrize('c_on', [True, False], ids=['c', 'python'])
def test_a_grid_transient_keeping_factors_is_its_run(c_on, monkeypatch):
    """The grid's Gear-2 run on the reduction with eight factors kept: every
    byte, step and count the run keeping none makes, its C loop's and its
    Python's alike; the step sizes it comes back to after each corner
    served by a kept factor."""
    ## (today's stepping: the ladder, which keeps factors of its own accord, off)
    a = _run(True, c_on, monkeypatch, step_ladder=False)
    monkeypatch.setattr(_cholklu, 'KEEP', 8)
    b = _run(True, c_on, monkeypatch, step_ladder=False)
    assert a[1:4] == b[1:4], (a[2:4], b[2:4])
    assert a[5]['hits'] == 0 and b[5]['hits'] > 0, (a[5], b[5])
    assert b[5]['factorizations'] < a[5]['factorizations'], (a[5], b[5])


def test_a_kept_factor_is_matched_on_every_value():
    """A kept factor is served only to the very values it was made on: J'
    and J'' each differ from J at one free node's diagonal, both places
    outside the entries the comparison reads first (`ck_keep`'s 64) -- J''
    meets J's kept factor and is factored anew, its bytes the run keeping
    none makes, no refactor check failing."""
    c = _grid()
    iref = c.get_node_index(gnd)

    def sequence(keep):
        sk = _solver(c)
        sk.keep_factors(keep)
        nnz = sk.nnz
        nsig = min(64, nnz)
        sig = {int((k + 0.5) * nnz / nsig) for k in range(nsig)}
        J, _iref = _matrix(c, 1e-12)
        places = []
        for name in ('n2_3', 'n3_2', 'n4_1', 'n1_4'):
            i = c.get_node_index(name)
            r = i - (i > iref)
            p = sk.Ap[r] + int(np.flatnonzero(sk.Ai[sk.Ap[r]:sk.Ap[r + 1]] == r)[0])
            if p not in sig:
                places.append(r)
        assert len(places) >= 2, places
        J1, J2 = J.copy(), J.copy()
        J1[places[0], places[0]] *= 1.5
        J2[places[1], places[1]] *= 1.5
        b = np.random.default_rng(4).standard_normal(sk.m)
        xs = [np.asarray(sk.solve(A, b, None)).tobytes() for A in (J, J1, J2)]
        return xs, sk
    want, _sk1 = sequence(1)
    got, sk = sequence(3)
    assert got == want
    assert sk.kept_counts()['hits'] == 0 and sk.kept_counts()['factorizations'] == 3, sk.kept_counts()
    assert (sk.residual_fallbacks, sk.residual_kept) == (0, 0)


## THE STEP LADDER (`Transient.step_ladder`, `transient._ladder_rung`, loop.c's twin; plan 261 L24 stage LB)

def test_a_ladder_rung_is_the_largest_power_of_two_below():
    """`_ladder_rung`: the largest anchor x 2^-k not over the size asked --
    within 1e-9 of a rung, that rung -- exact (each a halving); the size
    itself where it is not below the anchor, not positive, or the rung would
    fall under `minstep`."""
    from pycircuit.circuit.transient import _ladder_rung
    a = 1e-11
    for h, want in ((1e-11, 1e-11), (9.9e-12, 5e-12), (5e-12, 5e-12), (4.99e-12, 2.5e-12),
                    (1e-11 * (1 - 1e-12), 1e-11), (5e-12 * (1 - 1e-12), 5e-12), (1.3e-12, 1.25e-12)):
        assert _ladder_rung(h, a, 1e-18) == want, (h, want)
    assert _ladder_rung(2e-11, a, 1e-18) == 2e-11
    assert _ladder_rung(0.0, a, 1e-18) == 0.0 and np.isnan(_ladder_rung(np.nan, a, 1e-18))
    assert _ladder_rung(1.3e-12, a, 1.26e-12) == 1.3e-12
    r = _ladder_rung(3e-14, a, 1e-18)
    assert r == a * 2.0 ** -round(np.log2(a / r)) and r <= 3e-14 < 2 * r


def _ladder_run(c_on, monkeypatch, chol=True, edge=1e-10, **make):
    """The grid's run (`_run`'s, the threshold at 0 where `chol`): x, digest,
    statistics, counts as `_run`, the paths, the kept counts -- and the
    accepted steps' sizes, their end times, the pulse's corners, the anchor."""
    from pycircuit.circuit import transient as trmod
    got = {}
    real = trmod._SteppingLoop.__init__

    def init(self, *a, **k):
        real(self, *a, **k)
        got['loop'] = self
    monkeypatch.setattr(trmod._SteppingLoop, '__init__', init)
    out = _run(chol, c_on, monkeypatch, edge=edge, **make)
    loop = got['loop']
    t = np.asarray(loop.timelist, float)
    hs = np.diff(np.concatenate(([0.0], t)))
    ## (the run's own breakpoints -- the pulse's corners and its period's
    ## starts -- and the end)
    corners, tt = {3e-9}, 0.0
    while True:
        nb = float(loop.family.next_breakpoint(tt))
        if not nb > tt or nb > 3e-9:
            break
        corners.add(nb)
        tt = nb
    return out, hs, t, np.asarray(sorted(corners)), loop.ladder


@pytest.mark.parametrize('edge', [1e-10, 1e-11], ids=['edge', 'short-edge'])
@pytest.mark.parametrize('c_on', [True, False], ids=['c', 'python'])
def test_a_grid_on_the_ladder_steps_on_its_rungs(c_on, edge, monkeypatch):
    """The grid's Gear-2 run on the reduction, the ladder 'auto' (on: the
    matrix is the SPD reduction's): every accepted step a rung of the run's
    `max_step` (x 2^-k) or one landing on a corner or the end -- a pulse's
    edges shorter than the step too, where the corner cut's size is a
    quarter of the gap, rounded down; the run keeps `LADDER_KEEP` factors and
    its corners' matrices are served by them."""
    from pycircuit.circuit.transient import LADDER_KEEP
    out, hs, t, corners, anchor = _ladder_run(c_on, monkeypatch, edge=edge)
    assert anchor > 0.0
    k = np.round(np.log2(anchor / hs))
    ## (relative alone: numpy's default atol, 1e-8, would take every picosecond step for a rung)
    rung = np.isclose(hs, anchor * 2.0 ** -k, rtol=1e-9, atol=0.0)
    landing = np.array([np.min(np.abs(corners - ti)) <= 1e-12 * ti for ti in t])
    assert np.all(rung | landing), [(float(ti), float(h)) for ti, h, r, l in zip(t, hs, rung, landing) if not (r or l)][:5]
    assert rung.sum() > 0.8 * len(hs), (rung.sum(), len(hs))
    kept = out[5]
    assert kept['hits'] > 0 and kept['factorizations'] > 0 and kept['slots'] == LADDER_KEEP - 1, kept


def test_the_ladders_c_loop_is_its_python(monkeypatch):
    """On the ladder the C stepping loop's run is the Python's, byte for
    byte -- the same steps (loop.c's `ladder_rung`, `corner_size`), the same
    counts and kept factors; off the ladder another run."""
    c = _ladder_run(True, monkeypatch)
    p = _ladder_run(False, monkeypatch)
    assert c[0][1:4] == p[0][1:4] and c[0][5] == p[0][5], (c[0][2:4], p[0][2:4])
    assert np.array_equal(c[1], p[1])
    off = _ladder_run(True, monkeypatch, step_ladder=False)
    assert off[0][1] != c[0][1] and off[4] == 0.0


def test_the_ladder_is_off_where_the_matrix_is_not_the_reductions(monkeypatch):
    """'auto' off below `CHOLMOD_FROM` (KLU's run: the run with the ladder
    off, byte for byte, no anchor), and for a circuit the reduction declines
    (a diode in the grid); True puts it on a KLU run too -- its C loop the
    Python's."""
    klu = _ladder_run(True, monkeypatch, chol=False)
    off = _ladder_run(True, monkeypatch, chol=False, step_ladder=False)
    assert klu[4] == 0.0 and klu[0][1] == off[0][1] and not klu[0][3][0]
    on_c = _ladder_run(True, monkeypatch, chol=False, step_ladder=True)
    on_p = _ladder_run(False, monkeypatch, chol=False, step_ladder=True)
    assert on_c[4] > 0.0 and on_c[0][1] != off[0][1]
    assert on_c[0][1:4] == on_p[0][1:4]
    from pycircuit.circuit import transient as trmod
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0)
    c = _grid(6, 5)
    c['d'] = Diode('n2_2', gnd)
    got, real = [], trmod._SteppingLoop.__init__

    def init(self, *a, **k):
        real(self, *a, **k)
        got.append(self.ladder)
    monkeypatch.setattr(trmod._SteppingLoop, '__init__', init)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(c, toolkit=numeric)
        tr.solve(tend=3e-10, timestep=1e-11)
    assert got == [0.0] and not getattr(tr._get_linearsolver()._choice, 'cholmod', False)


def test_a_ladder_setting_not_known_is_refused():
    """`step_ladder` other than 'auto', True or False: refused, named."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(_grid(), toolkit=numeric, step_ladder='on')
    with pytest.raises(ValueError, match='step_ladder'):
        tr._step_ladder(None, type('R', (), {'max_step': 1e-11})())


## THE OWN SOLVE (`_cholklu`'s note, plan 261 L25)

def _own_solves(own, monkeypatch):
    """The reduction's solves on the grid with THE OWN SOLVE on or off: a
    factor, a refactor (h doubled), one right-hand side and three, the
    transposed solve and an owned factor's, each with its matrix and its
    right-hand side."""
    monkeypatch.setattr(_cholklu, 'OWN_SOLVE', own)
    c = _grid(16, 16)
    sk = _solver(c)
    rng = np.random.default_rng(5)
    got = []
    for h in (1e-12, 2e-12):
        J, _iref = _matrix(c, h)
        b = rng.standard_normal(J.shape[0])
        B = rng.standard_normal((J.shape[0], 3))
        got += [(J, b, np.asarray(sk.solve(J, b, None))), (J, B, np.asarray(sk.solve(J, B, None))),
                (J.T, B, np.asarray(sk.solve_transposed(J, B, None))),
                (J, b, np.asarray(sk.factor(J, None).solve(b)))]
    return got, sk


def _near(x, xr, A, b):
    """`x` the solve `xr` is to rounding, its backward error at the arithmetic's."""
    assert np.max(np.abs(x - xr)) <= 1e-13 * np.max(np.abs(xr))
    r = A @ x - b
    assert np.max(np.abs(r)) <= 1e-15 * np.max(abs(A) @ np.abs(x) + np.abs(b))


def test_the_own_solve_is_cholmod_solves_to_rounding(monkeypatch):
    """Every solve the reduction makes -- a factor's, a refactor's, several
    right-hand sides, the transposed solve, an owned factor's -- the one
    `cholmod_solve` makes to rounding, its backward error at the
    arithmetic's; counted as the own solve's and none `cholmod_solve`'s, and
    `OWN_SOLVE` off the other way round."""
    own, sko = _own_solves(True, monkeypatch)
    ref, skc = _own_solves(False, monkeypatch)
    for (A, b, x), (_A, _b, xr) in zip(own, ref):
        _near(x, xr, A, b)
    so, sc = sko.solve_counts(), skc.solve_counts()
    assert so['own'] > 0 and so['cholmod'] == 0, so
    assert sc == {'own': 0, 'cholmod': so['own']}, (so, sc)


def _laplacian(n, seed=7):
    """An n x n mesh's conductances (each between 1 and 2 S) with a small
    conductance from every node to ground: symmetric positive definite, CSC
    with both triangles and sorted rows -- a grid's node block."""
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    idx = np.arange(n * n).reshape(n, n)
    a = np.concatenate([idx[:, :-1].ravel(), idx[:-1, :].ravel()])
    b = np.concatenate([idx[:, 1:].ravel(), idx[1:, :].ravel()])
    g = 1.0 + rng.random(a.size)
    off = sp.coo_matrix((-g, (a, b)), shape=(n * n, n * n))
    off = off + off.T
    diag = -np.asarray(off.sum(axis=1)).ravel() + 1e-3 * (1.0 + rng.random(n * n))
    A = (off + sp.diags(diag)).tocsc()
    A.sort_indices()
    return A


def _entries(A, own, B):
    """`_cholklu`'s own entries on `A` (every row a node): the symbolic with
    THE OWN SOLVE on or off, a factor, the columns of `B` solved; the
    answers, the numeric's counts and the factor's statistics."""
    lib = _cholklu.library(linearsolver.KLUSolver()._lib)
    ip = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    m = A.shape[0]
    Ap, Ai = A.indptr.astype(np.int32), A.indices.astype(np.int32)
    Ax, kind = np.ascontiguousarray(A.data, dtype=np.float64), np.zeros(m, dtype=np.int32)
    sym = lib.analyze(m, ip(Ap), ip(Ai), ip(kind), -1)
    assert sym
    lib.own(ctypes.c_void_p(sym), int(own))
    num = lib.klu_factor(ip(Ap), ip(Ai), dp(Ax), ctypes.c_void_p(sym), None)
    assert num
    X = np.asfortranarray(B, dtype=np.float64).copy(order='F')
    assert lib.klu_solve(ctypes.c_void_p(sym), ctypes.c_void_p(num), m, X.shape[1], dp(X), None)
    counts, info = np.zeros(7), np.zeros(7)
    lib.num_info(ctypes.c_void_p(num), dp(counts))
    lib.info(ctypes.c_void_p(sym), dp(info))
    lib.klu_free_numeric(ctypes.byref(ctypes.c_void_p(num)), None)
    lib.klu_free_symbolic(ctypes.byref(ctypes.c_void_p(sym)), None)
    return X, counts, info


def test_the_own_solve_reads_narrow_and_wide_supernodes():
    """A 100 x 100 mesh's factor: thousands of supernodes below
    `CK_BLAS_FROM` columns (the plain loops) and wide ones with rows below
    their triangle (the BLAS calls, their gather and scatter) -- three
    right-hand sides each the one `cholmod_solve` makes to rounding, its
    backward error at the arithmetic's; one own solve a column, none
    `cholmod_solve`'s."""
    A = _laplacian(100)
    B = np.random.default_rng(11).standard_normal((A.shape[0], 3))
    X, counts, info = _entries(A, True, B)
    Xr, counts_r, _info = _entries(A, False, B)
    assert info[4] > 1000 and info[6] >= 2, info
    _near(X, Xr, A, B)
    assert tuple(counts[5:]) == (3, 0) and tuple(counts_r[5:]) == (0, 3), (counts, counts_r)


@pytest.mark.parametrize('c_on', [True, False], ids=['c', 'python'])
def test_a_grid_transient_on_the_own_solve_is_its_cholmod_solve_run(c_on, monkeypatch):
    """The grid's Gear-2 run on the reduction as a large grid runs (the
    ladder on, eight factors kept): every solve -- the Newton's, the error
    test's, the Python's -- the own solve's, the run the one on
    `cholmod_solve` step for step, the answers within rounding."""
    monkeypatch.setattr(_cholklu, 'OWN_SOLVE', True)
    a = _run(True, c_on, monkeypatch)
    monkeypatch.setattr(_cholklu, 'OWN_SOLVE', False)
    b = _run(True, c_on, monkeypatch)
    assert a[2:4] == b[2:4] and a[3][0], (a[2:4], b[2:4])
    assert float(np.max(np.abs(a[0] - b[0]))) <= 1e-12
    assert a[5]['hits'] > 0, a[5]
    assert a[6]['own'] > 50 and a[6]['cholmod'] == 0, a[6]
    assert b[6] == {'own': 0, 'cholmod': a[6]['own']}, (a[6], b[6])
