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
