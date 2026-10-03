import pytest
import pycircuit.circuit.circuit

## the shooting tests' shared helpers assert, and live outside `test_*`
## files (split out of test_analysis_shooting.py on 2026-09-27): without
## this their asserts would be plain, unexplained AssertionErrors
pytest.register_assert_rewrite('pycircuit.circuit.tests._shooting_fixtures',
                               'pycircuit.circuit.tests._shooting_elements')

@pytest.fixture
def restore_default_toolkit():
    """Put back the default toolkit a test switched (2026-10-03).  For the
    modules whose tests deliberately run on the symbolic toolkit by setting
    `circuit.default_toolkit` (`pytestmark = pytest.mark.usefixtures(
    'restore_default_toolkit')`).  Until then an autouse fixture reset the
    toolkit after EVERY test, which silently repaired the 17 writes the leak
    detector then reported -- and would have repaired a library leak too.
    Any other test that leaves the toolkit changed now fails (the detector,
    `pycircuit/_testing/leaks.py`)."""
    old = pycircuit.circuit.circuit.default_toolkit
    yield
    pycircuit.circuit.circuit.default_toolkit = old
