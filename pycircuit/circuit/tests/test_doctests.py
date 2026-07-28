# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""circuit.py and elements.py carry many doctests documenting their public
API, but none of them ran under pytest: they were gated behind
``if __name__ == "__main__": doctest.testmod()``, which only fires when a
module is executed as a script. Enabling ``--doctest-modules`` project-wide
surfaces 15 collection errors in unrelated, largely legacy modules
(``post/cds``, ``post/jwdb``, ``sim/gnucap`` -- some still importing Python-2
``StringIO`` -- and a couple of broken example scripts), which is a separate,
larger cleanup (see architecture.md P13). This scopes the fix to the two
modules the coverage survey actually flagged as highest-risk: most
depended-upon, least directly tested.

Two real, previously-undetected bugs were hiding behind this exact gap:

* ``Quantity.__repr__`` unconditionally raised ``ValueError("APA")`` before
  its actual logic (introduced 2010, unreachable ever since -- the same
  shadowed-dead-code shape as ``func.Tanh.fprime``'s 2009-2026 bug, P10).
* ``Circuit.name_state_vector`` had two bugs: ``x[:len(self.nodes)][0]``
  double-indexed into a scalar instead of slicing, and its branch-naming
  loop unpacked ``enumerate(zip(...))`` (a 2-tuple) into three names. Both
  were invisible because nothing in the live codebase calls this method and
  its own doctest's example happened not to exercise the branch path.

This test makes running every doctest in both modules a permanent part of
the ordinary suite, so a shadowed ``raise`` or a bad unpacking can't go
undetected behind a dead example again.
"""
import doctest

from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit import elements as elements_module


def test_circuit_module_doctests():
    results = doctest.testmod(circuit_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in circuit.py failed' %
        (results.failed, results.attempted))


def test_elements_module_doctests():
    results = doctest.testmod(elements_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in elements.py failed' %
        (results.failed, results.attempted))
