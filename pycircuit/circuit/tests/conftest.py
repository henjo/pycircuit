import pytest
import pycircuit.circuit.circuit

## the shooting tests' shared helpers assert, and live outside `test_*`
## files (split out of test_analysis_shooting.py on 2026-09-27): without
## this their asserts would be plain, unexplained AssertionErrors
pytest.register_assert_rewrite('pycircuit.circuit.tests._shooting_fixtures',
                               'pycircuit.circuit.tests._shooting_elements')

@pytest.fixture(autouse=True)
def reset_global_toolkit(request):
    """Ensure the default toolkit is reset to numeric after every test.
    This prevents tests that test SymbolicToolkit from leaking it into other tests.

    It runs BEFORE the leak detector looks (it is torn down first), so it
    reports what it repairs to the detector (2026-10-03) -- until the 28
    unrestored writes are fixed and it goes (robust testing, stage 3).
    """
    from pycircuit.circuit.toolkit import numeric
    yield
    old = pycircuit.circuit.circuit.default_toolkit
    if old is not numeric:
        from pycircuit._testing import leaks
        leaks.note_toolkit_reset(request.config, request.node.nodeid, old)
    pycircuit.circuit.circuit.default_toolkit = numeric
