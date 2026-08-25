"""`uic=True` must seed a generated state from the model, not from a NAME.

`_apply_element_ics` used to read `element.iparv.ic` and then, on the
state path, never use the value -- so what looked like "does this
element declare an initial condition?" was really "does it happen to
have a parameter SPELLED `ic`?".

A generated model whose seed is called anything else fell through it.
The DC pin was still perfectly correct, because that comes from the
compiled `state_ic()`, so the operating point looked right and only
the `uic=True` waveform was wrong -- starting at zero, with no error
and no warning.  Found by writing a memristor whose state parameter is
naturally called `ic`... and then renaming it.

The tests below are written so they FAIL against the old gate: each
declares its seed under a name that is not `ic`.
"""

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit import gnd
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.elements import SubCircuit, VS, R
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, idt
from pycircuit.circuit.transient import Transient
from pycircuit.utilities.param import Parameter


class StateNamedX0(Behavioural):
    """An integrator whose seed is `x0`, not `ic`.

    `I(p,m) <+ g * (idt(V(p,m), x0) - target)` is not physics; it is the
    smallest thing that puts a SEEDED state on a row `uic` has to find.
    """

    instparams = [Parameter(name='x0', desc='initial state', unit='V',
                            default=0.0),
                  Parameter(name='g', desc='gain', unit='S', default=1e-3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * idt(b.V, x0))            # noqa: F821


def _build(x0):
    cm.default_toolkit = numeric
    c = SubCircuit(toolkit=numeric)
    c['vs'] = VS(cm.Node('a'), gnd, v=0.0)
    c['r'] = R(cm.Node('a'), gnd, r=1e3)
    c['e'] = StateNamedX0(cm.Node('a'), gnd, x0=x0, g=1e-3)
    return c


def test_the_element_declares_a_state_seed_without_a_parameter_called_ic():
    """The precondition.  If this fails the rest proves nothing."""
    e = StateNamedX0(cm.Node('a'), gnd, x0=0.75)
    e.update_iparv()
    assert getattr(e.iparv, 'ic', None) is None, \
        'this element must NOT have a parameter named ic, or the test is vacuous'
    assert hasattr(e, 'state_ic'), 'no compiled state_ic to seed from'
    assert getattr(e, 'IC_KIND', 'current') == 'state'
    ## The compiled seed carries the value, which is the whole point:
    ## the information was always there, the gate just never asked for it.
    assert any(abs(v - 0.75) < 1e-12 for _row, v in e.state_ic())


@pytest.mark.parametrize('x0', [0.75, -1.25, 3.0])
def test_uic_seeds_the_state_from_the_models_own_seed(x0):
    """The regression.  Against the old NAME gate every case starts at 0."""
    c = _build(x0)
    tran = Transient(c, toolkit=numeric)
    n = len(c.nodes) + len(c.branches)
    x_init = tran._apply_element_ics(np.zeros(n), gnd)
    rows = c.elementnodemap['e']
    seeded = [x_init[rows[r]] for r, _v in c['e'].state_ic()]
    assert seeded, 'no state row was located'
    assert seeded[0] == pytest.approx(x0, rel=1e-12), \
        'uic did not seed the state from state_ic()'


def test_a_zero_seed_is_still_a_seed_and_is_not_confused_with_absence():
    """`x0 = 0` must take the same path as any other value.

    Asserted because the natural way to write the fix -- treat a falsy
    seed as "no seed" -- reintroduces the bug for the one value a state
    most often starts at.
    """
    c = _build(0.0)
    tran = Transient(c, toolkit=numeric)
    n = len(c.nodes) + len(c.branches)
    x_init = tran._apply_element_ics(np.zeros(n), gnd)
    rows = c.elementnodemap['e']
    seeded = [x_init[rows[r]] for r, _v in c['e'].state_ic()]
    assert seeded[0] == pytest.approx(0.0, abs=1e-15)


def test_the_nested_subcircuit_guard_sees_it_too():
    """`_descendant_has_ic` carried the same NAME check.

    Under-detecting there is the dangerous direction: the guard exists to
    REFUSE an initial condition the flat walk cannot place, so failing to
    see one drops it silently -- exactly what the guard is for.
    """
    cm.default_toolkit = numeric
    inner = SubCircuit(toolkit=numeric)
    inner['e'] = StateNamedX0(cm.Node('a'), gnd, x0=0.5)
    outer = SubCircuit(toolkit=numeric)
    outer['x'] = inner
    tran = Transient(outer, toolkit=numeric)
    assert tran._descendant_has_ic(inner) is True
