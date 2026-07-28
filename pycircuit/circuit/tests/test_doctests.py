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

``volterra.py`` (P14, the legacy corners this same doctest sweep found)
joined this file the same way: its own doctests turned up a chain of
real bugs the same gap had hidden -- a stale absolute import, ``self.c``
where the base class sets ``self.cir``, an undefined ``symbolic_linsolve``,
``InternalResult`` where only ``InternalResultDict`` exists (which itself
had a matching bug: ``__len__`` read ``self.results``, ``__init__`` sets
``self.items``), a stray ``self`` parameter on a module-level function, and
``NLVCCS`` (a test-only nonlinear element defined in the same file) hitting
the same G-vs-i mismatch P9 fixed for the shipped elements, plus a
sympy-version change in how ``diff`` applies to a plain array. All fixed;
``Volterra.solve()``/``.run()`` remain a deliberately unfinished stub (the
orchestration calling ``K()`` per element was commented out in 2008 and
never completed) -- the doctest reflects that honestly rather than
inventing the missing algorithm.

Fixing ``volterra.py`` in turn required fixing ``symbolicapprox.py``, whose
own doctest turned up two independent bugs in one line: ``series(...,
point=0, ...)`` (modern sympy renamed the keyword to ``x0``), and
``.subs({'t': 1}).removeO()`` -- substituting a concrete value into an
expression that still carries a symbolic ``O(t**3)`` term collapses the
whole thing to ``0`` before ``removeO()`` ever runs, silently. Swapped to
``.removeO().subs({'t': 1})``; the docstring's own claimed example output
was itself wrong (never verified, same 18-year-invisible pattern) and is
now the actual, checked value.
"""
import doctest

from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit import elements as elements_module
from pycircuit.circuit import volterra as volterra_module
from pycircuit.circuit import symbolicapprox as symbolicapprox_module


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


def test_volterra_module_doctests():
    results = doctest.testmod(volterra_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in volterra.py failed' %
        (results.failed, results.attempted))


def test_symbolicapprox_module_doctests():
    results = doctest.testmod(symbolicapprox_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in symbolicapprox.py failed' %
        (results.failed, results.attempted))
