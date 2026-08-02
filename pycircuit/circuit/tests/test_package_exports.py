"""STAGE 10 -- every analysis must be reachable from the package root.

`Transient` and `PSS` were the only two that were not. `DC`, `DCSweep`, `AC` and
`Noise` all arrive through the star-imports in ``circuit/__init__.py``, so

    from pycircuit.circuit import Transient

-- the import every transient example and every doc page uses -- raised
``ImportError`` while its neighbours worked. It is a one-line defect that
survived the whole of stages 0-12, which is the argument for a test rather than
a fix: nothing else in the suite imports the package the way a user does.
"""
import importlib

import pytest


ANALYSES = ['DC', 'DCSweep', 'AC', 'Noise', 'Transient', 'PSS']


@pytest.mark.parametrize('name', ANALYSES)
def test_analysis_is_importable_from_the_package_root(name):
    mod = importlib.import_module('pycircuit.circuit')
    assert hasattr(mod, name), \
        '%s is not exported from pycircuit.circuit' % name


@pytest.mark.parametrize('name', ANALYSES)
def test_the_from_import_form_works(name):
    """`hasattr` and `from ... import` can disagree, so both are checked."""
    mod = __import__('pycircuit.circuit', fromlist=[name])
    assert getattr(mod, name) is not None


def test_the_exported_object_is_the_real_class():
    """Guards against a name being shadowed by something else star-imported.

    `transient.py` and `shooting.py` both do `from ...analysis import *`, so a
    star-import here would re-export their transitive imports; the named imports
    exist to stop that, and this asserts the result rather than the mechanism.
    """
    from pycircuit.circuit import Transient, PSS
    from pycircuit.circuit.transient import Transient as T_direct
    from pycircuit.circuit.shooting import PSS as P_direct
    assert Transient is T_direct
    assert PSS is P_direct
