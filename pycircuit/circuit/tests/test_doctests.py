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
from pycircuit.circuit import mos as mos_module
from pycircuit.circuit import nportanalysis as nportanalysis_module


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


def test_mos_module_doctests():
    """``mos.py`` joined this file at stage 5+.5, and the route it took is the
    argument for the file existing.

    ``MOS_ACM`` sat in that module for years: unconstructable (``super(MOS, self)``
    from a class with ``MOS`` not in its MRO), a verbatim copy of ``MOS`` with a
    mis-described parameter and a noise PSD referring to a symbol nothing binds.
    The only thing that would ever have caught it was the doctest in its own
    docstring, gated behind ``if __name__ == "__main__"``. It was deleted at 5+.2.

    Adding the module here then surfaced two more defects of the same shape in
    ``MOS``'s own example, neither related to the deleted class:

    * ``c = SubCircuit()`` took the *numeric* default toolkit while the parameters
      were ``Symbol(...)``, so construction died in ``elements.update``. The
      example predates the toolkit split and never ran after it.
    * ``freqs=array([Symbol('s')])`` tripped ``nportanalysis``'s
      ``assert not isiterable(freqs)``. That assertion was right in intent -- more
      than one symbolic frequency is silently reduced to the last -- and wrong in
      form, testing the container rather than the count. A length-1 array is one
      frequency and gives the identical answer, so it is now accepted and a longer
      one raises a ``ValueError`` that says why.

    Three defects, one module, all invisible for the same reason.
    """
    results = doctest.testmod(mos_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in mos.py failed' %
        (results.failed, results.attempted))


def test_nportanalysis_module_doctests():
    """17 of 33 failing when this module was first measured, from four causes.

    * ``import symbolic`` -- a Python-2 implicit relative import. One line, and
      because every later example in that block referred to names it was meant to
      bind, it cascaded into ten ``NameError`` failures.
    * ``res['mu'].y[0]`` -- results are no longer wrapped in a ``Waveform``. In the
      *numeric* path they are plain arrays, so it is ``res['mu'][0]``; in the
      *symbolic* path they are scalars, so it is ``res['mu']``. That asymmetry is
      the current API, documented here rather than quietly changed.
    * ``Symbol('R1', real=True)`` -- real is not enough for sympy to discharge the
      ``Abs`` in the noise expression, so ``Svn`` could not reduce. Resistances are
      positive, and saying so lets it.
    * ``print(mu)`` against a pasted ``(0.1+0j)`` -- the value is 0.1 to within
      7e-15 and the literal was never going to survive a change of BLAS or
      platform. Formatted to 4 significant figures instead.

    **The physics documented here was correct throughout.** ``simplify(old - new)``
    is 0 for the noise expression: the old expected ``(4*R1*R2*kT + 4*kT*R1**2)/R2``
    and the current ``4*R1*T*k*(R1 + R2)/R2`` are the same quantity, differing only
    in how ``kT`` is now composed and in what sympy could reduce. Only the
    plumbing had drifted.
    """
    results = doctest.testmod(nportanalysis_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in nportanalysis.py failed' %
        (results.failed, results.attempted))


def test_symbolicapprox_module_doctests():
    results = doctest.testmod(symbolicapprox_module, verbose=False)
    assert results.failed == 0, (
        '%d of %d doctests in symbolicapprox.py failed' %
        (results.failed, results.attempted))
