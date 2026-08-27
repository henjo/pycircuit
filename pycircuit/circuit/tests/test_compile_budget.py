# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""The compile-time budget warning (roadmap sec. 36, S4 phase 1).

**A forgotten `var()` announces itself in compile TIME and nothing
else.**  Measured on `PspMosLongChannel`: with every hold dropped it does
not finish compiling in ten minutes, against 62 s held -- an expression
that mentions a previous result twice doubles the tree, so a model nested
`n` deep has ``2**n`` occurrences.

Nothing else is a usable signal.  The emitted line count goes *down* when
a hold is dropped (the definition is inlined rather than named), the
individual held expressions are small either way (median 3 ops across
PSP's 342), and the finiteness reach did not move on the very holds whose
docstring blames evaluation order.  So the clock is the instrument.

It warns rather than raises: the threshold is a property of the machine
as much as of the model.
"""

import warnings

import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import hdl, gnd
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution
from pycircuit.utilities.param import Parameter

cm.default_toolkit = numeric


def _build(name):
    """Compile a trivial class under whatever threshold is set."""
    return hdl.BehaviouralMeta(name, (Behavioural,), dict(
        terminals=['p', 'm'],
        instparams=[Parameter('r', default=1e3)],
        analog=staticmethod(
            lambda p, m: (Contribution(Branch(p, m, 'b').I,
                                       Branch(p, m, 'b').V / r),)),  # noqa
        __module__=__name__))


def test_a_normal_compile_is_silent():
    """The library must not warn.  Every class in `elements_hdl` compiles
    far inside the threshold, and a check that cried wolf would be turned
    off within a day."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        _build('BudgetQuiet')
    assert not [x for x in w if 'to compile' in str(x.message)], \
        [str(x.message) for x in w]


def test_a_pathological_compile_warns_and_says_why():
    """Threshold dropped to zero so any compile is 'pathological'.

    What is asserted is not the timing -- that would be a flaky test --
    but that the warning FIRES on the class's own name and names the
    cause an author can act on.
    """
    was = hdl.COMPILE_WARN_SECONDS
    hdl.COMPILE_WARN_SECONDS = 0.0
    try:
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            _build('BudgetLoud')
    finally:
        hdl.COMPILE_WARN_SECONDS = was
    msgs = [str(x.message) for x in w if 'to compile' in str(x.message)]
    assert len(msgs) == 1, msgs
    m = msgs[0]
    assert 'BudgetLoud' in m
    ## it must name the cause, not merely report a number
    assert 'var()' in m and 'more than once' in m
    assert 'hold_value.py' in m


def test_the_warning_is_a_warning_and_not_an_error():
    """A build that refused to finish on a slow machine would be worse
    than the problem it reports."""
    was = hdl.COMPILE_WARN_SECONDS
    hdl.COMPILE_WARN_SECONDS = 0.0
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cls = _build('BudgetNotFatal')
        assert cls is not None
        el = cls('p', gnd)
        el.update_iparv()
    finally:
        hdl.COMPILE_WARN_SECONDS = was


def test_the_threshold_is_above_the_largest_model_in_the_tree():
    """`PspMosLongChannel` is ~62 s.  A threshold at or below that would
    warn on a correct model every time it is compiled cold."""
    assert hdl.COMPILE_WARN_SECONDS > 100.0
