import pytest
import pycircuit.circuit.circuit

@pytest.fixture(autouse=True)
def reset_global_toolkit():
    """Ensure the default toolkit is reset to numeric after every test.
    This prevents tests that test SymbolicToolkit from leaking it into other tests.
    """
    from pycircuit.circuit.toolkit import numeric
    yield
    pycircuit.circuit.circuit.default_toolkit = numeric
