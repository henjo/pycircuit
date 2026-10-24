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


@pytest.mark.parametrize('stages,held', [(120, True), (20, False), (300, False)])
def test_a_dc_solve_holds_blas_at_one_thread_where_its_dense_lu_is_threaded(
        stages, held, monkeypatch):
    """`DC.solve` holds the limit where the AutoSolver factors dense from
    `_blas.THREADED_FROM` unknowns up to its sparse floor -- a diode ladder
    of 122 -- and leaves the pool alone below it (22) and on KLU (302): the
    C Newton's ladders make hundreds of small solves there, and the threads
    cost up to 30x (voter25's operating point, 2026-10-13).  Read inside
    the solve, the pool opened to 4 threads first."""
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
