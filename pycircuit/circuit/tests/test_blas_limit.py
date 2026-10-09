"""The transient's BLAS thread limit (`transient._single_threaded_blas`)."""
import pytest

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
