"""The BLAS thread limit (`_blas._single_threaded_blas`): held by the transient, PSS, the DC's dense LU, the
small-signal analyses below their size and every PAC surface."""
import pytest

from pycircuit.circuit import _blas
from pycircuit.circuit import transient as tmod


def test_the_blas_limit_scans_the_libraries_once_and_still_limits(monkeypatch):
    """A threadpoolctl controller rescans every loaded library when it is
    built (3.4 M instructions, measured 2026-10-04 -- two steps of a
    20-MosLevel1 transient); the limit builds one at most once while no
    module is imported, and each use still holds BLAS at one thread."""
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    import threadpoolctl
    probe = threadpoolctl.ThreadpoolController()
    made = []
    orig = threadpoolctl.ThreadpoolController.__init__

    def counting(self, *args, **kwargs):
        made.append(1)
        return orig(self, *args, **kwargs)
    monkeypatch.setattr(threadpoolctl.ThreadpoolController, '__init__', counting)
    for _ in range(3):
        with tmod._single_threaded_blas():
            inside = [c.num_threads for c in probe.lib_controllers if c.user_api == 'blas']
        assert inside and all(t == 1 for t in inside)
    assert len(made) <= 1


def test_a_pss_solve_holds_blas_at_one_thread(monkeypatch):
    """`PSS.solve` holds the limit for its whole run, as `Transient.solve`
    does: its walks drive the inner transient's steps directly, and until
    2026-10-09 every dense solve ran on OpenBLAS's pool -- radau's stage
    system of a 49-unknown Gilbert cell (147, past the threading size) took
    221 s against 1.2 s on one thread.  Read inside the shooting itself."""
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    import threadpoolctl

    from pycircuit.circuit import circuit
    from pycircuit.circuit.elements import VSin, R, C, gnd
    from pycircuit.circuit.shooting.pss import PSS
    probe = threadpoolctl.ThreadpoolController()
    if not any(c.user_api == 'blas' for c in probe.lib_controllers):
        pytest.skip('no BLAS threadpoolctl can see')
    c = circuit.SubCircuit()
    c['V'] = VSin('in', gnd, va=1.0, freq=1e3)
    c['R'] = R('in', 'out', r=1e3)
    c['C'] = C('out', gnd, c=1e-7)
    pss = PSS(c, method='radau')
    seen = []
    orig = PSS._shoot

    def shoot(self, run):
        seen.append([lc.num_threads for lc in probe.lib_controllers if lc.user_api == 'blas'])
        return orig(self, run)
    monkeypatch.setattr(PSS, '_shoot', shoot)
    ## (the suite's workers may run pinned to one thread already: open the
    ## pool first, or the check is vacuous)
    with probe.limit(limits=4, user_api='blas'):
        assert all(lc.num_threads == 4 for lc in threadpoolctl.ThreadpoolController().lib_controllers
                   if lc.user_api == 'blas')
        pss.solve(period=1e-3, timestep=1e-5)
    assert seen and all(t == 1 for t in seen[0]), seen


@pytest.mark.parametrize('stages,held', [(120, True), (20, False), (300, True)])
def test_a_dc_solve_holds_blas_at_one_thread_where_its_dense_lu_is_threaded(
        stages, held, monkeypatch):
    """`DC.solve` holds the limit where the AutoSolver solves from
    `_blas.THREADED_FROM` unknowns -- a diode ladder of 122, factored
    dense, and of 302, on KLU -- and leaves the pool alone below it (22):
    the C Newton's ladders make hundreds of small dense solves (the threads
    cost up to 30x: voter25's operating point, 2026-10-13), and at size the
    solve's other dense calls are threaded too (an active chain's DC at
    4 000 stages 1.34x).  Read inside the solve, the pool opened to 4
    threads first."""
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    import threadpoolctl

    from pycircuit.circuit import circuit, gnd
    from pycircuit.circuit import elements_hdl as eh
    from pycircuit.circuit.dcanalysis import DC
    from pycircuit.circuit.elements import VS, R, SubCircuit
    probe = threadpoolctl.ThreadpoolController()
    if not any(c.user_api == 'blas' for c in probe.lib_controllers):
        pytest.skip('no BLAS threadpoolctl can see')
    c = SubCircuit()
    c['v'] = VS('n0', gnd, v=5.0)
    for k in range(1, stages + 1):
        c[f'r{k}'] = R(f'n{k - 1}', f'n{k}', r=100.0)
        c[f'd{k}'] = eh.DiodeSpiceHdl(f'n{k}', gnd, IS=1e-14 * (1 + k % 7))
    seen = []
    orig = DC._solve_dc

    def spy(self, x0=None):
        seen.append([lc.num_threads for lc in probe.lib_controllers if lc.user_api == 'blas'])
        return orig(self, x0)
    monkeypatch.setattr(DC, '_solve_dc', spy)
    with probe.limit(limits=4, user_api='blas'):
        DC(c, toolkit=circuit.numeric).solve()
    assert seen and all(t == (1 if held else 4) for t in seen[0]), (c.n, seen)


@pytest.mark.parametrize('kind', ['ac', 'noise', 'transimpedance'])
@pytest.mark.parametrize('stages', [1, 450], ids=['small', 'large'])
def test_the_small_signal_analyses_hold_blas_at_one_thread_below_its_size(
        kind, stages, monkeypatch):
    """`AC.solve`, `Noise.solve` and `TransimpedanceAnalysis.solve` hold the
    limit for their whole run below `_blas.AC_THREADS_FROM` unknowns
    (`blas_held_below`): a dense complex solve or more at each frequency,
    1.19-1.48x from 171 to 295 unknowns on the gold decks' operating
    points; from 440 the threads pay, and the pool is left alone (a ladder
    of 450 sections).  Read inside the analysis, the pool opened to 4
    threads first."""
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    import numpy as np
    import threadpoolctl

    from pycircuit.circuit import analysis_ss, circuit
    from pycircuit.circuit.circuit import Branch
    from pycircuit.circuit.elements import IS, VS, C, R, gnd
    probe = threadpoolctl.ThreadpoolController()
    if not any(c.user_api == 'blas' for c in probe.lib_controllers):
        pytest.skip('no BLAS threadpoolctl can see')
    c = circuit.SubCircuit()
    c['vs'] = VS('in', gnd, vac=1.0)
    c['r'] = R('in', 'n1', r=1e3)
    for k in range(1, stages):
        c[f'r{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
    c['rl'] = R(f'n{stages}', 'out', r=1e3)
    c['c'] = C('out', gnd, c=1e-9)
    c['is'] = IS('out', gnd, iac=0.0)
    seen = []
    orig = analysis_ss.remove_row_col

    def spy(*args, **kwargs):
        seen.append([lc.num_threads for lc in probe.lib_controllers if lc.user_api == 'blas'])
        return orig(*args, **kwargs)
    monkeypatch.setattr(analysis_ss, 'remove_row_col', spy)
    freqs = np.array([1e3, 1e6])
    with probe.limit(limits=4, user_api='blas'):
        if kind == 'ac':
            analysis_ss.AC(c, toolkit=circuit.numeric).solve(freqs)
        elif kind == 'noise':
            analysis_ss.Noise(c, inputsrc='vs', outputnodes=('out', gnd),
                              toolkit=circuit.numeric).solve(freqs)
        else:
            analysis_ss.TransimpedanceAnalysis(c, toolkit=circuit.numeric).solve(
                1e3, [Branch('out', gnd)])
        outside = [lc.num_threads for lc in threadpoolctl.ThreadpoolController().lib_controllers
                   if lc.user_api == 'blas']
    held = c.n < _blas.AC_THREADS_FROM
    assert held == (stages == 1), c.n
    assert seen and all(t == (1 if held else 4) for s in seen for t in s), (c.n, seen)
    assert outside and all(t == 4 for t in outside), outside


@pytest.mark.parametrize('surface', ['solve', 'pnoise'])
def test_the_periodic_small_signal_analysis_holds_blas_at_one_thread(surface, monkeypatch):
    """Every public method of `PAC` holds the limit for its call
    (`_blas.blas_held_methods`): many circuit-size operations through the
    period, where the threads never paid -- one thread 0.97-1.76x for
    `solve`, 1.11-1.57x for `pnoise`, from 53 to 803 unknowns (2026-10-24).
    Read inside the surface, the pool opened to 4 threads first."""
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    import threadpoolctl

    from pycircuit.circuit import circuit
    from pycircuit.circuit.elements import C, R, VSin, gnd
    from pycircuit.circuit.shooting.pac import PAC
    from pycircuit.circuit.shooting.pss import PSS
    from pycircuit.circuit.simwarnings import AccuracyWarning
    from pycircuit.circuit.tests._warnpolicy import quiet
    probe = threadpoolctl.ThreadpoolController()
    if not any(c.user_api == 'blas' for c in probe.lib_controllers):
        pytest.skip('no BLAS threadpoolctl can see')
    c = circuit.SubCircuit()
    c['vs'] = VSin('in', gnd, vac=1.0, va=1.0, freq=1e3)
    c['R'] = R('in', 'out', r=1e3, noisy=True)
    c['C'] = C('out', gnd, c=1e-7)
    pss = PSS(c, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        pss.solve(period=1e-3, timestep=1e-3 / 60, maxiterations=40)
    pac = PAC(c, toolkit=circuit.numeric)
    seen = []
    orig = PAC._check_circuit

    def spy(self, pss_):
        seen.append([lc.num_threads for lc in probe.lib_controllers if lc.user_api == 'blas'])
        return orig(self, pss_)
    monkeypatch.setattr(PAC, '_check_circuit', spy)
    with probe.limit(limits=4, user_api='blas'), quiet(AccuracyWarning):
        if surface == 'solve':
            pac.solve(pss, [300.0])
        else:
            pac.pnoise(pss, 300.0, 'out')
        outside = [lc.num_threads for lc in threadpoolctl.ThreadpoolController().lib_controllers
                   if lc.user_api == 'blas']
    assert seen and all(t == 1 for s in seen for t in s), seen
    assert outside and all(t == 4 for t in outside), outside


#: a process of its own: warm every import, keep the hold's controller (a held
#: DC), then load CHOLMOD, then read the pool inside a held DC
_DLOPEN_SCRIPT = """
import threadpoolctl
from pycircuit.circuit import _cholklu, circuit, gnd, linearsolver
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, R, SubCircuit
def ladder(stages):
    c = SubCircuit()
    c['v'] = VS('n0', gnd, v=5.0)
    for k in range(1, stages + 1):
        c[f'r{k}'] = R(f'n{k - 1}', f'n{k}', r=100.0)
        c[f'd{k}'] = eh.DiodeSpiceHdl(f'n{k}', gnd, IS=1e-14 * (1 + k % 7))
    return c
probe = threadpoolctl.ThreadpoolController
before = len([c for c in probe().lib_controllers if c.user_api == 'blas'])
klu = linearsolver.KLUSolver()._lib
DC(ladder(300), toolkit=circuit.numeric).solve()
DC(ladder(120), toolkit=circuit.numeric).solve()
if _cholklu.library(klu) is None:
    print('SKIP', _cholklu.STATUS)
    raise SystemExit
pool = probe()
after = len([c for c in pool.lib_controllers if c.user_api == 'blas'])
if after <= before:
    print('SKIP CHOLMOD brings no BLAS of its own here')
    raise SystemExit
seen, orig = [], DC._solve_dc
def spy(self, x0=None):
    seen.append([c.num_threads for c in pool.lib_controllers if c.user_api == 'blas'])
    return orig(self, x0)
DC._solve_dc = spy
with pool.limit(limits=4, user_api='blas'):
    DC(ladder(300), toolkit=circuit.numeric).solve()
print('SEEN', before, after, seen)
"""


def _in_a_process(script, **environ):
    """`script`'s last line, run by a Python of its own on this tree
    (`environ`: variables set in its environment, None removing one); a line
    starting SKIP skips."""
    import os
    import subprocess
    import sys

    import pycircuit
    if not tmod.blas_single_thread_available():
        pytest.skip('no threadpoolctl')
    root = os.path.dirname(os.path.dirname(os.path.abspath(pycircuit.__file__)))
    env = {**os.environ, 'PYTHONPATH': root}
    for name, value in environ.items():
        if value is None:
            env.pop(name, None)
        else:
            env[name] = str(value)
    r = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True, timeout=600, check=False,
                       cwd=root, env=env)
    out = r.stdout.strip().splitlines()
    assert r.returncode == 0 and out, r.stderr[-2000:]
    if out[-1].startswith('SKIP'):
        pytest.skip(out[-1][5:])
    return out[-1]


def test_a_blas_a_loader_opens_is_held_too():
    """CHOLMOD, opened by `_cholklu`'s C (`dlopen`, no import), brings the
    system's OpenBLAS: a hold whose controller was kept from before -- every
    module imported, a held DC solved -- holds that library too once CHOLMOD
    has loaded (`_blas.library_loaded`), read inside a DC at 302 unknowns
    with the pool opened to 4 threads first.  In a process of its own: a
    worker may have loaded CHOLMOD already (the suite failed on that order)."""
    import ast
    line = _in_a_process(_DLOPEN_SCRIPT)
    _, before, after, seen = line.split(' ', 3)
    seen = ast.literal_eval(seen)
    assert int(after) > int(before), line
    assert seen and len(seen[0]) == int(after) and all(t == 1 for t in seen[0]), line


#: a process of its own, every OpenBLAS at the threads the test asks: the
#: first transient on CHOLMOD (a grid, the threshold at 0), run inside a hold
#: of the caller's, loads it inside the transient's own hold; the BLAS
#: libraries' thread counts read as the C loop starts, in that run and the
#: next, after the first inside the caller's hold, then after both
_LOAD_IN_HOLD_SCRIPT = """
import contextlib, threadpoolctl
from pycircuit.circuit import _blas, _cholklu, _tran_loop_c, numeric
from pycircuit.circuit.elements import IPulse
from pycircuit.circuit.tests.test_cholklu import _grid
from pycircuit.circuit.transient import Transient
_cholklu.CHOLMOD_FROM = 0
def blas():
    return {c.filepath: c.num_threads for c in threadpoolctl.ThreadpoolController().lib_controllers
            if c.user_api == 'blas'}
before = blas()
seen, orig = [], _tran_loop_c.run
def spy(loop):
    seen.append(blas())
    return orig(loop)
_tran_loop_c.run = spy
runs = []
for k in range(2):
    c = _grid(16, 16)
    c['ip'] = IPulse('n8_8', 'n0_0', i1=0.0, i2=5e-3, td=1e-10, tr=1e-10, tf=1e-10, pw=4e-10, per=2e-9)
    tr = Transient(c, toolkit=numeric)
    with _blas._single_threaded_blas() if k == 0 else contextlib.nullcontext():
        tr.solve(tend=1e-9, timestep=1e-11)
        if k == 0:
            outer = blas()
    if not tr._get_linearsolver()._choice.cholmod or not seen:
        print('SKIP', 'CHOLMOD ' + _cholklu.STATUS + ' or no C loop')
        raise SystemExit
    runs.append(seen[-1])
new = sorted(set(runs[0]) - set(before))
if not new:
    print('SKIP CHOLMOD brings no BLAS of its own here')
    raise SystemExit
print('SEEN', repr([before, runs + [outer], blas(), new]).replace(' ', ''))
"""


def test_a_blas_loaded_inside_a_hold_is_held_to_its_end():
    """The first transient on CHOLMOD loads it inside the transient's own
    hold, after the hold's limit was set: the BLAS it brings is held at one
    thread from then to the outermost hold's end (`_blas.library_loaded`) --
    a hold of the caller's around the transient still holds it after the
    transient's own has closed -- and back at its own count after.  On ibmpg3t that library ran on 24 threads in the
    run that loaded it (1.30x the time of one thread, other bits than the runs
    after it).  Every OpenBLAS at 4 threads here (fewer cores: as many),
    read as the C loop starts in that run and in the next, and in the
    caller's hold after the first."""
    import ast
    import os
    k = min(4, os.cpu_count() or 1)
    if k < 2:
        pytest.skip('one core: no thread count to tell apart')
    line = _in_a_process(_LOAD_IN_HOLD_SCRIPT, OPENBLAS_NUM_THREADS=k)
    before, runs, after, new = ast.literal_eval(line.split(' ', 1)[1])
    for run in runs:
        assert set(run) == set(after) and all(t == 1 for t in run.values()), line
    assert all(after[f] == t for f, t in before.items()), line
    assert all(after[f] == k for f in new), line


#: a process of its own: the BLAS libraries' thread counts before and after
#: CHOLMOD loads (outside any hold), and whether OPENBLAS_NUM_THREADS is in
#: the environment after
_START_AT_ONE_SCRIPT = """
import os, threadpoolctl
from pycircuit.circuit import _cholklu, linearsolver
def blas():
    return {c.filepath: c.num_threads for c in threadpoolctl.ThreadpoolController().lib_controllers
            if c.user_api == 'blas'}
klu = linearsolver.KLUSolver()._lib
before = blas()
if _cholklu.library(klu) is None:
    print('SKIP', _cholklu.STATUS)
    raise SystemExit
after = blas()
new = sorted(set(after) - set(before))
if not new:
    print('SKIP CHOLMOD brings no BLAS of its own here')
    raise SystemExit
print('SEEN', repr([before, after, new, 'OPENBLAS_NUM_THREADS' in os.environ]).replace(' ', ''))
"""


def test_cholmods_blas_starts_at_one_thread_unless_the_environment_says():
    """CHOLMOD's BLAS (the system OpenBLAS, opened with it by `_cholklu`'s C)
    reads OPENBLAS_NUM_THREADS as it loads; started on every core and held at
    one, CHOLMOD's solve on ibmpg3t ran 6.5 % slower.  So `_cholklu.library`
    sets the variable to 1 for that load alone when the environment does not
    set it.  In a process of its own with OPENBLAS_NUM_THREADS unset and
    OMP_NUM_THREADS -- OpenBLAS's fallback -- at 4 (fewer cores: as many):
    the libraries loaded before at 4, which is where CHOLMOD's would start
    without the variable; after, those at 4 still, CHOLMOD's at one, the
    variable gone.  (Set, as in the test above, CHOLMOD's follows it.)"""
    import ast
    import os
    k = min(4, os.cpu_count() or 1)
    if k < 2:
        pytest.skip('one core: no thread count to tell apart')
    line = _in_a_process(_START_AT_ONE_SCRIPT, OPENBLAS_NUM_THREADS=None, GOTO_NUM_THREADS=None, OMP_NUM_THREADS=k)
    before, after, new, still_set = ast.literal_eval(line.split(' ', 1)[1])
    assert before and all(t == k for t in before.values()), line
    assert all(after[f] == t for f, t in before.items()), line
    assert all(after[f] == 1 for f in new), line
    assert not still_set, line
