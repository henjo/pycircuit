import pytest
import pycircuit.circuit.circuit

## the shooting tests' shared helpers assert, and live outside `test_*`
## files (split out of test_analysis_shooting.py on 2026-09-27): without
## this their asserts would be plain, unexplained AssertionErrors
pytest.register_assert_rewrite('pycircuit.circuit.tests._shooting_fixtures',
                               'pycircuit.circuit.tests._shooting_elements')

@pytest.fixture(autouse=True)
def reset_global_toolkit():
    """Ensure the default toolkit is reset to numeric after every test.
    This prevents tests that test SymbolicToolkit from leaking it into other tests.
    """
    from pycircuit.circuit.toolkit import numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = numeric
