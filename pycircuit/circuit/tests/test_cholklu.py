"""THE SPD REDUCTION ON CHOLMOD (`_cholklu`): a linear grid's transient
matrix solved through KLU's entries over CHOLMOD's supernodal Cholesky --
the voltage sources' nodes fixed, the inductors' currents eliminated, the
free nodes' symmetric positive definite system factored -- as a dense solve
makes it; the plan declining what is not a grid; a matrix CHOLMOD will not
factor KLU's; a transient on it its KLU run, its C the Python's; the kept
factors, the step ladder, the glue's own solve on CHOLMOD's factor, the
first Newton iterate taken on an affine circuit, and the affine step."""
import ctypes
import hashlib
import types
import warnings

import numpy as np
import pytest

from pycircuit.circuit import (
    _cholklu,
    _paths,
    _stamp_plan,
    _tran_newton,
    linearsolver,
    numeric,
)
from pycircuit.circuit.analysis import jacobian_structure
from pycircuit.circuit.circuit import SubCircuit, defaultepar, gnd
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
    ## (today's stepping on both: the ladder and the first iterate taken are
    ## the reduction's alone -- their own tests below)
    make = {'chord_jacobian': chord, 'step_ladder': False, 'linear_newton': False}
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


def _kept_solves(keep, max_bytes=None, c=None, **budget):
    c = _grid() if c is None else c
    sk = _solver(c)
    sk.keep_factors(keep, max_bytes, **budget)
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


def _resident():
    """This process's resident bytes, as THE BUDGET reads them."""
    import os
    with open('/proc/self/statm') as f:
        return int(f.read().split()[1]) * os.sysconf('SC_PAGESIZE')


def test_kept_factors_stay_within_the_memory_left(monkeypatch):
    """THE BUDGET: a new slot only while the process's resident memory, the
    slot and the margin fit under the limit -- a limit below the memory
    already resident, or a margin over the room the limit leaves, keeps
    none; a limit with room keeps every slot `keep` allows, `KEEP_BYTES`
    then not binding (a `max_bytes` given still does); the answers the run
    keeping none gives each time."""
    want, _sk1 = _kept_solves(1)
    room = min(7, len(set(_HS)) - 1)
    for limit, margin, slots in ((_resident() - 64 * 2 ** 20, 0.0, 0),
                                 (_resident() + 2 ** 30, 2.0 * 2 ** 30, 0),
                                 (_resident() + 2 ** 30, 0.0, room)):
        got, sk = _kept_solves(8, limit=limit, margin=margin)
        assert got == want, (limit, margin)
        assert sk.kept_counts()['slots'] == slots, (limit, margin, slots, sk.kept_counts())
    monkeypatch.setattr(_cholklu, 'KEEP_BYTES', 1.0)
    for budget, slots in (({'limit': _resident() + 2 ** 30, 'margin': 0.0}, room), ({}, 0),
                          ({'max_bytes': 1.0, 'limit': _resident() + 2 ** 30, 'margin': 0.0}, 0)):
        got, sk = _kept_solves(8, **budget)
        assert got == want, budget
        assert sk.kept_counts()['slots'] == slots, (budget, slots, sk.kept_counts())


def test_the_memory_limit_is_the_least_of_the_cgroups_and_the_machines(tmp_path):
    """THE BUDGET's limit: the least of the `memory.max` up the process's
    cgroup tree ('max' none) and its resident memory plus MemAvailable;
    None where neither can be read."""
    import os
    proc, cg = tmp_path / 'proc', tmp_path / 'cg'
    (proc / 'self').mkdir(parents=True)
    (proc / 'self' / 'statm').write_text(f'900000 {2 ** 30 // os.sysconf("SC_PAGESIZE")} 1000 1 0 1 0\n')
    (proc / 'meminfo').write_text('MemTotal:       32000000 kB\nMemAvailable:    8388608 kB\n')
    (proc / 'self' / 'cgroup').write_text('0::/user.slice/run-1.scope\n')
    (cg / 'user.slice' / 'run-1.scope').mkdir(parents=True)
    (cg / 'user.slice' / 'memory.max').write_text('max\n')

    def limit():
        return _cholklu.memory_limit(str(proc), str(cg))
    assert limit() == 9 * 2 ** 30
    (cg / 'user.slice' / 'run-1.scope' / 'memory.max').write_text(f'{6 * 2 ** 30}\n')
    assert limit() == 6 * 2 ** 30
    (cg / 'user.slice' / 'memory.max').write_text(f'{5 * 2 ** 30}\n')
    assert limit() == 5 * 2 ** 30
    (proc / 'meminfo').write_text('MemAvailable:    1048576 kB\n')
    assert limit() == 2 * 2 ** 30
    (proc / 'self' / 'statm').unlink()
    assert limit() == 5 * 2 ** 30
    assert _cholklu.memory_limit(str(tmp_path / 'none'), str(tmp_path / 'none')) is None


def test_the_memory_limit_the_user_sets(tmp_path, monkeypatch):
    """THE BUDGET's limit the user sets (`KEEP_MEMORY`, env
    PYCIRCUIT_CHOLMOD_KEEP_MEMORY): bytes, or K, M, G, T (powers of 1024);
    the least of it, the cgroup's and the machine's binding, it alone where
    neither can be read; a value not bytes warned of and not used."""
    import os
    import subprocess
    import sys
    pb = _cholklu.parse_bytes
    assert (pb('123'), pb('8G'), pb(' 1.5k'), pb('2T'), pb('3M'), pb('0'), pb('')) == (
        123.0, 8.0 * 2 ** 30, 1536.0, 2.0 * 2 ** 40, 3.0 * 2 ** 20, 0.0, None)
    for bad in ('8GB', 'lots', '-1G', 'nan', 'inf', 'G'):
        with pytest.raises(ValueError):
            pb(bad)
    proc, cg = tmp_path / 'proc', tmp_path / 'cg'
    (proc / 'self').mkdir(parents=True)
    (proc / 'self' / 'statm').write_text(f'900000 {2 ** 30 // os.sysconf("SC_PAGESIZE")} 1000 1 0 1 0\n')
    (proc / 'meminfo').write_text('MemAvailable:    8388608 kB\n')
    (proc / 'self' / 'cgroup').write_text('0::/run-1.scope\n')
    (cg / 'run-1.scope').mkdir(parents=True)
    (cg / 'run-1.scope' / 'memory.max').write_text(f'{6 * 2 ** 30}\n')
    monkeypatch.setattr(_cholklu, 'KEEP_MEMORY', 4.0 * 2 ** 30)
    assert _cholklu.memory_limit(str(proc), str(cg)) == 4 * 2 ** 30
    monkeypatch.setattr(_cholklu, 'KEEP_MEMORY', 7.0 * 2 ** 30)
    assert _cholklu.memory_limit(str(proc), str(cg)) == 6 * 2 ** 30
    assert _cholklu.memory_limit(str(tmp_path / 'none'), str(tmp_path / 'none')) == 7 * 2 ** 30
    code = ('import warnings; warnings.simplefilter("always"); '
            'from pycircuit.circuit import _cholklu; print(repr(_cholklu.KEEP_MEMORY))')

    def read(value):
        env = dict(os.environ, PYCIRCUIT_CHOLMOD_KEEP_MEMORY=value, JAX_PLATFORMS='cpu')
        r = subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True, timeout=300,
                           check=False)
        assert r.returncode == 0, r.stderr[-2000:]
        return r.stdout.split()[-1], r.stderr
    got, err = read('2G')
    assert got == repr(2.0 * 2 ** 30) and 'PYCIRCUIT_CHOLMOD_KEEP_MEMORY' not in err, (got, err[-2000:])
    got, err = read('8GB')
    assert got == 'None' and "PYCIRCUIT_CHOLMOD_KEEP_MEMORY='8GB' is not a number of bytes" in err, (got, err[-2000:])


def test_cholmod_is_the_library_given_or_built_else_the_systems(tmp_path):
    """THE LIBRARY: `PYCIRCUIT_CHOLMOD_LIB` opened first, else the build in
    pycircuit's cache, else the system's; one not had -- missing, or a
    CHOLMOD of another version -- falls back to the system's, warned;
    `STATUS` names the one opened (each case a process of its own: CHOLMOD
    is opened once a process)."""
    import os
    import shutil
    import subprocess
    import sys

    from pycircuit.circuit import _hdl_cache
    code = ('import warnings; warnings.simplefilter("always"); '
            'from pycircuit.circuit import _cholklu, linearsolver; '
            '_cholklu.library(linearsolver.KLUSolver()._lib); print(_cholklu.STATUS); '
            'print(sorted({ln.split()[-1] for ln in open("/proc/self/maps") if "libcholmod" in ln}))')

    def run(lib=None, py=None):
        env = {k: v for k, v in os.environ.items() if k != 'PYCIRCUIT_CHOLMOD_LIB'}
        env.update(JAX_PLATFORMS='cpu', XDG_CACHE_HOME=str(tmp_path / 'xdg'),
                   PYCIRCUIT_HDL_CACHE_DIR=_hdl_cache.cache_dir())
        if lib is not None:
            env['PYCIRCUIT_CHOLMOD_LIB'] = lib
        r = subprocess.run([sys.executable, '-c', py or code], env=env, capture_output=True, text=True,
                           timeout=300, check=False)
        assert r.returncode == 0, r.stderr[-2000:]
        return r.stdout.splitlines()[-2:], r.stderr
    ## the system's, its path as the loader finds it
    system = run(py='import ctypes; ctypes.CDLL("libcholmod.so.5"); print("-"); '
                    'print([ln.split()[-1] for ln in open("/proc/self/maps") if "libcholmod" in ln][0])')[0][1]
    (out, err) = run()
    assert out == ["on (the system's libcholmod)", repr([system])] and 'CHOLMOD' not in err, (out, err[-2000:])
    given = tmp_path / 'given' / 'libcholmod.so.5'
    given.parent.mkdir()
    shutil.copy(system, given)
    (out, err) = run(str(given))
    assert out == [f'on ({given})', repr([str(given)])] and 'CHOLMOD' not in err, (out, err[-2000:])
    built = tmp_path / 'xdg' / 'pycircuit' / 'cholmod' / 'lib' / 'libcholmod.so.5'
    built.parent.mkdir(parents=True)
    shutil.copy(system, built)
    (out, err) = run()
    assert out == [f'on ({built})', repr([str(built)])] and 'CHOLMOD' not in err, (out, err[-2000:])
    ## a path given wins over the build; one missing, or another version's, falls back to the system's
    (out, err) = run(str(given))
    assert out == [f'on ({given})', repr([str(given)])], (out, err[-2000:])
    built.unlink()
    shim = tmp_path / 'shim'
    shim.mkdir()
    (shim / 'v.c').write_text('int cholmod_version(int *v) { v[0] = 4; v[1] = 0; v[2] = 0; return 4000; }\n')
    subprocess.run([shutil.which('cc') or 'gcc', '-shared', '-fPIC', '-o', str(shim / 'libshim.so'),
                    str(shim / 'v.c'), system], check=True, capture_output=True)
    for bad in (tmp_path / 'missing.so', shim / 'libshim.so'):
        (out, err) = run(str(bad))
        assert out[0] == "on (the system's libcholmod)" and system in out[1], (bad, out, err[-2000:])
        assert f"CHOLMOD '{bad}' could not be opened at the version of the header" in err, (bad, err[-2000:])


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


def test_the_ladder_keeps_its_factors_within_the_memory_left(monkeypatch):
    """The ladder's kept factors under THE BUDGET, its limit
    `memory_limit()`: one no slot fits under, or one leaving less than
    `KEEP_MARGIN`, or the user's `KEEP_MEMORY` at 0, keeps none, the run the
    bytes, steps and counts of the run with room; none read: `KEEP_BYTES`
    binds."""
    room = _ladder_run(True, monkeypatch)[0]
    ## (the user's limit through `memory_limit()` itself: 0 is no slot, not no limit)
    monkeypatch.setattr(_cholklu, 'KEEP_MEMORY', 0.0)
    held = _ladder_run(True, monkeypatch)[0]
    monkeypatch.setattr(_cholklu, 'KEEP_MEMORY', None)
    monkeypatch.setattr(_cholklu, 'memory_limit', lambda: 1.0)
    none = _ladder_run(True, monkeypatch)[0]
    ## (room for the slots, not for `KEEP_MARGIN` beside them)
    monkeypatch.setattr(_cholklu, 'memory_limit', lambda: _resident() + 0.5 * _cholklu.KEEP_MARGIN)
    tight = _ladder_run(True, monkeypatch)[0]
    monkeypatch.setattr(_cholklu, 'memory_limit', lambda: None)
    monkeypatch.setattr(_cholklu, 'KEEP_BYTES', 1.0)
    capped = _ladder_run(True, monkeypatch)[0]
    assert room[5]['slots'] > 0 and room[5]['hits'] > 0, room[5]
    assert none[5]['slots'] == 0 and none[5]['hits'] == 0, none[5]
    assert tight[5]['slots'] == 0 and capped[5]['slots'] == 0, (tight[5], capped[5])
    assert none[1:3] == room[1:3] and capped[1:3] == room[1:3], (none[2], room[2], capped[2])
    assert held[5]['slots'] == 0 and held[5]['hits'] == 0 and held[1:3] == room[1:3], held[5]


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
    assert sc == {'own': 0, 'cholmod': so['own'], 'threaded': 0}, (so, sc)


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
    counts, info = np.zeros(8), np.zeros(7)
    lib.num_info(ctypes.c_void_p(num), dp(counts))
    lib.info(ctypes.c_void_p(sym), dp(info))
    lib.klu_free_numeric(ctypes.byref(ctypes.c_void_p(num)), None)
    lib.klu_free_symbolic(ctypes.byref(ctypes.c_void_p(sym)), None)
    return X, counts, info


def test_the_budget_counts_the_slot_it_would_make():
    """THE BUDGET counts the slot it would make (L's values and J's): a 300
    x 300 mesh's factor, its slot tens of MB -- with the limit half a slot
    over the resident memory no slot is made (the least recently used one
    takes the factor), with it two slots over one is."""
    import math
    A = _laplacian(300)
    lib = _cholklu.library(linearsolver.KLUSolver()._lib)
    ip = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
    dp = lambda a: a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    m = A.shape[0]
    Ap, Ai, kind = A.indptr.astype(np.int32), A.indices.astype(np.int32), np.zeros(m, dtype=np.int32)
    vals = [np.ascontiguousarray(A.data * (1.0 + 0.125 * k), dtype=np.float64) for k in range(4)]
    sym = ctypes.c_void_p(lib.analyze(m, ip(Ap), ip(Ai), ip(kind), -1))
    assert sym
    lib.keep(sym, 3, math.inf, 0.0, 0.0)
    num = ctypes.c_void_p(lib.klu_factor(ip(Ap), ip(Ai), dp(vals[0]), sym, None))
    assert num

    def refactor(k):
        assert lib.klu_refactor(ip(Ap), ip(Ai), dp(vals[k]), sym, num, None)
        out = np.zeros(8)
        lib.num_info(num, dp(out))
        return out
    try:
        out = refactor(1)
        one = 8.0 * (out[4] + A.nnz)
        assert out[2] == 1 and one > 16 * 2 ** 20, out
        lib.keep(sym, 3, math.inf, _resident() + 0.5 * one, 0.0)
        out = refactor(2)
        assert out[2] == 1 and out[1] == 3, out
        lib.keep(sym, 3, math.inf, _resident() + 2.0 * one, 0.0)
        out = refactor(3)
        assert out[2] == 2 and out[1] == 4, out
    finally:
        lib.klu_free_numeric(ctypes.byref(num), None)
        lib.klu_free_symbolic(ctypes.byref(sym), None)


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
    assert tuple(counts[5:]) == (3, 0, 0) and tuple(counts_r[5:]) == (0, 3, 0), (counts, counts_r)


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
    assert b[6] == {'own': 0, 'cholmod': a[6]['own'], 'threaded': 0}, (a[6], b[6])



## THE THREADS (`_cholklu`'s note): the own solve's subtrees on the C units' pool

def _laplacian3(n, seed=5):
    """An n x n x n mesh's conductances (each between 1 and 2 S), a small
    conductance from every node to ground: a factor whose supernodal tree
    branches at every level -- many subtrees, their updates to the rows
    above them, a top of wide separators."""
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    idx = np.arange(n ** 3).reshape(n, n, n)
    pairs = [(idx[:-1].ravel(), idx[1:].ravel()), (idx[:, :-1].ravel(), idx[:, 1:].ravel()),
             (idx[:, :, :-1].ravel(), idx[:, :, 1:].ravel())]
    a = np.concatenate([q[0] for q in pairs])
    b = np.concatenate([q[1] for q in pairs])
    off = sp.coo_matrix((-(1.0 + rng.random(a.size)), (a, b)), shape=(n ** 3, n ** 3))
    off = off + off.T
    diag = -np.asarray(off.sum(axis=1)).ravel() + 1e-3 * (1.0 + rng.random(n ** 3))
    A = (off + sp.diags(diag)).tocsc()
    A.sort_indices()
    return A


@pytest.fixture
def _threads(monkeypatch):
    """`set_threads` for the test, a factor of any size taken on the pool
    (`THREADS_FROM` 0); the count put back after."""
    from pycircuit.circuit import _tran_core
    monkeypatch.setattr(_cholklu, 'THREADS_FROM', 0)
    old = _tran_core.get_threads()
    yield _tran_core.set_threads
    _tran_core.set_threads(old)


def _differ(x, y):
    """How many of `x`'s values are not `y`'s, bit for bit."""
    return int(np.sum(np.asarray(x).view(np.uint64) != np.asarray(y).view(np.uint64)))


@pytest.mark.parametrize('nthreads', [2, 3, 4, 8])
@pytest.mark.parametrize('mesh', ['plane', 'cube'])
def test_the_solve_on_threads_is_the_solve_on_one(mesh, nthreads, _threads):
    """A mesh's factor solved with its subtrees on the pool -- each
    subtree's forward updates to the rows above it kept, then replayed with
    the top's supernodes in the serial order; the backward top first, then
    the subtrees -- three right-hand sides the bytes of the solve on one
    thread; every own solve the pool's, none on one."""
    A = _laplacian(100) if mesh == 'plane' else _laplacian3(22)
    B = np.random.default_rng(11).standard_normal((A.shape[0], 3))
    _threads(1)
    X1, c1, _info = _entries(A, True, B)
    _threads(nthreads)
    X, c, info = _entries(A, True, B)
    assert info[5] >= 1, info
    assert _differ(X, X1) == 0, _differ(X, X1)
    assert tuple(c1[5:]) == (3, 0, 0) and tuple(c[5:]) == (3, 0, 3), (c1, c)


def test_the_threads_take_a_factor_from_their_size_and_their_switch(monkeypatch, _threads):
    """A factor of fewer than `THREADS_FROM` values, or the switch off, is
    solved on the caller alone -- the same bytes."""
    A = _laplacian(100)
    B = np.random.default_rng(11).standard_normal((A.shape[0], 3))
    _threads(4)
    X, c, _info = _entries(A, True, B)
    monkeypatch.setattr(_cholklu, 'THREADS_FROM', 1 << 40)
    _threads(4)
    Xs, cs, _info = _entries(A, True, B)
    monkeypatch.setattr(_cholklu, 'THREADS_FROM', 0)
    monkeypatch.setattr(_cholklu, 'THREADS', False)
    _threads(4)
    Xo, co, _info = _entries(A, True, B)
    assert c[7] == 3 and cs[7] == 0 and co[7] == 0, (c, cs, co)
    assert _differ(Xs, X) == 0 and _differ(Xo, X) == 0


#: the own solve's answers on `_laplacian(18)` (three right-hand sides), as
#: the loops made them before a one-column supernode was taken inline and the
#: backward sums eight columns at a time: their bytes, kept (CHOLMOD's
#: ordering of the mesh moves them too -- its supernodes in `info` first)
OWN_SOLVE_BYTES = '6028ff72719e8574'


@pytest.mark.parametrize('nthreads', [1, 4])
def test_the_own_solve_keeps_the_bytes_it_had(nthreads, _threads):
    """An 18 x 18 mesh's factor -- one-column supernodes, others of eight
    columns or more with rows below, none wide enough for the BLAS (so the
    bytes are the loops' alone) -- solved on one thread and on four: the
    answers the loops made before, bit for bit."""
    import hashlib
    A = _laplacian(18)
    B = np.random.default_rng(11).standard_normal((A.shape[0], 3))
    _threads(nthreads)
    X, c, info = _entries(A, True, B)
    assert info[5] == 0, info
    assert hashlib.sha1(X.tobytes()).hexdigest()[:16] == OWN_SOLVE_BYTES
    assert c[7] == (3 if nthreads > 1 else 0), c


@pytest.mark.parametrize('c_on', [True, False], ids=['c', 'python'])
def test_a_grid_transient_on_threads_is_its_run_on_one(c_on, monkeypatch, _threads):
    """The grid's Gear-2 run with the own solve on four threads: every
    solve the pool's, the run the one on one thread byte for byte -- its
    answers, steps and counts."""
    _threads(1)
    one = _run(True, c_on, monkeypatch)
    _threads(4)
    four = _run(True, c_on, monkeypatch)
    assert four[1] == one[1] and four[2:4] == one[2:4], (four[2:4], one[2:4])
    assert one[6]['threaded'] == 0 and four[6]['threaded'] == four[6]['own'] > 50, (one[6], four[6])


## THE FIRST ITERATE ACCEPTED (`Transient.linear_newton`, plan 261 L26)

class _VaractorC(C):
    """A capacitor whose C moves with its voltage: its own `C`."""

    def C(self, x, epar=defaultepar, params_tree=None):
        return super().C(x, epar) * (1.0 + (x[0] - x[1]) ** 2)


class _LeakyR(R):
    """A resistor with a cubic leak: its own `i`, the class's `G`."""

    def i(self, x, epar=defaultepar, params_tree=None):
        d = x[0] - x[1]
        return super().i(x, epar) + 1e-3 * d ** 3 * np.array([1.0, -1.0])


class _ClampR(R):
    """A resistor that limits its Newton step."""

    def limit(self, x, x0, epar=defaultepar):
        return x


def test_the_certificate_is_the_constant_stamps_with_nothing_limiting(monkeypatch):
    """`Transient._affine`: the grid -- resistors, capacitors, inductors,
    sources, a pulse -- is affine; with a capacitor of its own `C`, a
    resistor of its own `i`, one that limits, a diode, or a resistor's `G`
    shadowed on the instance it is not; nor any circuit with the stamp plan
    off."""
    def affine(c):
        return Transient(c, toolkit=numeric)._affine()
    c = _grid()
    c['ip'] = IPulse('n2_2', 'n0_0', i1=0.0, i2=1e-3, td=1e-10, tr=1e-10, tf=1e-10, pw=4e-10)
    assert affine(c)
    for name, el in (('cv', _VaractorC('n1_1', gnd, c=1e-12)),
                     ('rl', _LeakyR('n1_1', 'n2_2', r=5.0)),
                     ('rc', _ClampR('n1_1', 'n2_2', r=5.0)),
                     ('d', Diode('n1_1', gnd))):
        c = _grid()
        c[name] = el
        assert not affine(c), name
    c = _grid()
    r = c['r0']
    r.G = types.MethodType(type(r).G, r)
    assert not affine(c)
    monkeypatch.setattr(_stamp_plan, 'ENABLED', False)
    assert not affine(_grid())


def _ln_run(chol, c_on, monkeypatch, nx=16, diode=False, fixed=False, **make):
    """The grid's run (`_run`'s, on an `nx` x `nx` mesh; a diode across it
    where asked; on the fixed grid where asked): x, digest, (accepted,
    rejected, Newton iterations), the paths, the loop's `first_iterate`,
    whether the run's mark outlived it."""
    from pycircuit.circuit import _tran_loop_c, _tran_lte_c, _tran_newton_c
    from pycircuit.circuit import transient as trmod
    monkeypatch.setattr(_cholklu, 'CHOLMOD_FROM', 0 if chol else 1 << 40)
    for mod in (_tran_loop_c, _tran_lte_c, _tran_newton_c):
        monkeypatch.setattr(mod, 'ENABLED', c_on)
    got = {}
    real = trmod._SteppingLoop.__init__

    def init(self, *a, **k):
        real(self, *a, **k)
        got['loop'] = self
    monkeypatch.setattr(trmod._SteppingLoop, '__init__', init)
    c = _grid(nx, nx)
    c['ip'] = IPulse(f'n{nx // 2}_{nx // 2}', 'n0_0', i1=0.0, i2=5e-3, td=1e-10, tr=1e-10,
                     tf=1e-10, pw=4e-10, per=2e-9)
    if diode:
        c['d'] = Diode('n2_2', 'n0_0')
    before = _paths.snapshot()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tr = Transient(c, toolkit=numeric, **make)
        res = tr.solve(tend=3e-9, timestep=1e-11, **({'fixed_timestep': True} if fixed else {}))
    st = tr.statistics
    x = np.asarray(res.x, float)
    return (x, hashlib.sha1(x.tobytes()).hexdigest(),
            (st.accepted_steps, st.rejected_steps, st.newton_iterations), _paths.since(before),
            got['loop'].first_iterate, '_first_iterate' in tr.__dict__)


@pytest.mark.parametrize('c_on', [True, False], ids=['c', 'python'])
def test_the_first_iterate_is_taken_on_the_reduction_and_where_asked(c_on, monkeypatch):
    """'auto': on the reduction (CHOLMOD) one Newton iteration an attempt,
    the steps those of the run without it, the answers within 1e-12; below
    the threshold (KLU) not taken -- True takes it there, the switch off
    never; the run's mark gone after it."""
    off = _ln_run(True, c_on, monkeypatch, linear_newton=False)
    on = _ln_run(True, c_on, monkeypatch)
    assert on[4] and not off[4] and not on[5] and not off[5]
    acc, rej, its = on[2]
    assert its == acc + rej and on[2][:2] == off[2][:2] and off[2][2] > its + 20, (on[2], off[2])
    assert float(np.max(np.abs(on[0] - off[0]))) <= 1e-12
    klu = _ln_run(False, c_on, monkeypatch)
    forced = _ln_run(False, c_on, monkeypatch, linear_newton=True)
    assert not klu[4] and forced[4] and not forced[5]
    assert forced[2][2] == sum(forced[2][:2]) and klu[2][2] > forced[2][2] + 20, (klu[2], forced[2])
    assert float(np.max(np.abs(forced[0] - klu[0]))) <= 1e-12
    monkeypatch.setattr(_tran_newton, 'LINEAR_NEWTON', False)
    sw = _ln_run(False, c_on, monkeypatch, linear_newton=True)
    assert not sw[4] and sw[1] == klu[1]


@pytest.mark.parametrize('kind', ['cholmod', 'klu', 'dense'])
def test_the_first_iterate_s_c_is_the_python(kind, monkeypatch):
    """With the first iterate taken -- on the reduction ('auto'), on KLU and
    on a dense matrix (True) -- the C Newton and stepping loop serve the
    run, and it is the Python's byte for byte."""
    chol, nx = kind == 'cholmod', 6 if kind == 'dense' else 16
    make = {} if chol else {'linear_newton': True}
    c = _ln_run(chol, True, monkeypatch, nx=nx, **make)
    p = _ln_run(chol, False, monkeypatch, nx=nx, **make)
    assert c[4] and p[4]
    assert c[1:3] == p[1:3] and c[2][2] == sum(c[2][:2]), (c[2], p[2])
    d = c[3]
    assert d.get('loop_c:steps', 0) > 50 and d.get('newton_c:served', 0) > 50, d
    assert not any(k.startswith(('newton_c:bail', 'newton_c:no')) for k in d), d


def test_the_first_iterate_is_neither_the_chord_s_nor_a_non_affine_run_s(monkeypatch):
    """The chord on the reduction ('auto' marks the run) iterates as it did;
    True on the grid with a diode across it is the run without it; a
    stage method's run is never marked."""
    from pycircuit.circuit.integrator import TRBDF2Integrator
    a = _ln_run(True, True, monkeypatch, chord_jacobian=True)
    b = _ln_run(True, True, monkeypatch, chord_jacobian=True, linear_newton=False)
    assert a[4] and a[1:3] == b[1:3], (a[2], b[2])
    a = _ln_run(False, True, monkeypatch, nx=8, diode=True, linear_newton=True)
    b = _ln_run(False, True, monkeypatch, nx=8, diode=True, linear_newton=False)
    assert not a[4] and a[1:3] == b[1:3], (a[2], b[2])
    a = _ln_run(False, True, monkeypatch, nx=6, integrator=TRBDF2Integrator(),
                linear_newton=True)
    assert not a[4] and not a[5]


@pytest.mark.parametrize('chol, trap, fixed', [(True, False, False), (False, False, False),
                                               (True, True, False), (False, True, False),
                                               (False, True, True)],
                         ids=['cholmod', 'klu', 'cholmod-trap', 'klu-trap', 'klu-trap-fixed'])
def test_the_affine_step_is_the_run_s_bytes_and_makes_nothing_twice(chol, trap, fixed,
                                                                     monkeypatch):
    """THE AFFINE STEP (csrc/newton.h) on the first iterate's run -- on the
    reduction ('auto') and on KLU (True), Gear-2's and the trapezoid's (whose
    J reads h and the formula, not a0: on the fixed grid Euler's restart
    after a corner and the trapezoid's next step share a0 = 0 and h): the
    run without it byte for byte; every error test of the C loop's attempts
    (an adaptive run's) on the Newton's own numeric (none refactors, none
    counted as the seed's), and the evaluations whose J kAx held forming
    none -- some, not all (the step sizes or the formula move); the switch
    off, neither."""
    from pycircuit.circuit.integrator import TrapezoidalIntegrator
    make = {} if chol else {'linear_newton': True}
    if trap:
        make = dict(make, integrator=TrapezoidalIntegrator(), linear_newton=True)
    ## (off first: a process's first fixed-grid run declines the C loop --
    ## `_sources`' stamp one behind the epoch -- not this test's subject)
    monkeypatch.setattr(_tran_newton, 'AFFINE_STEP', False)
    off = _ln_run(chol, True, monkeypatch, fixed=fixed, **make)
    monkeypatch.setattr(_tran_newton, 'AFFINE_STEP', True)
    on = _ln_run(chol, True, monkeypatch, fixed=fixed, **make)
    assert on[4] and on[1:3] == off[1:3], (on[2], off[2])
    d, e = on[3], off[3]
    steps = d.get('loop_c:steps', 0)
    assert steps > 50 and 0 < d.get('loop_c:affine_j', 0) < steps, d
    assert not any('affine' in k for k in e), e
    if not fixed:
        assert d.get('loop_c:affine_lte', 0) >= steps, d
        assert not d.get('loop_c:lte_on_pattern') and not d.get('lte:seed'), d
        assert e.get('loop_c:lte_on_pattern', 0) >= steps, e
