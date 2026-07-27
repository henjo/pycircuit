# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""JAX is an optional accelerator, and this is what keeps it that way.

``toolkit.py`` guards its JAX import so that ``jaxtoolkit`` becomes ``None``
when JAX is absent.  That guard was silently defeated for a while by an
unconditional ``from . import _jaxtoolkit`` higher up the same file: since
``_jaxtoolkit`` does ``import jax`` at module level, importing pycircuit at all
required JAX -- including the numeric and symbolic paths, which never touch it.

The failure was invisible in an environment where JAX happens to be installed,
which is why it survived.  These tests make it visible by removing JAX.
"""

import builtins
import sys

import pytest


def _reimport_without(blocked):
    """Import ``pycircuit.circuit.toolkit`` with ``blocked`` unavailable."""
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == blocked or name.startswith(blocked + '.'):
            raise ImportError('no %s (blocked by test)' % blocked)
        return real_import(name, *args, **kwargs)

    saved = {k: v for k, v in sys.modules.items()
             if k.startswith('pycircuit') or k == blocked
             or k.startswith(blocked + '.')}
    for key in saved:
        del sys.modules[key]

    builtins.__import__ = fake_import
    try:
        import pycircuit.circuit.toolkit as tk
        return tk
    finally:
        builtins.__import__ = real_import
        for key in list(sys.modules):
            if key.startswith('pycircuit'):
                del sys.modules[key]
        sys.modules.update(saved)


def test_toolkit_imports_without_jax():
    """The package must import when JAX is not installed."""
    tk = _reimport_without('jax')
    assert tk.numeric is not None
    assert tk.symbolic is not None


def test_jaxtoolkit_is_none_without_jax():
    """The guard's actual contract: absent JAX yields a None singleton."""
    tk = _reimport_without('jax')
    assert tk.jaxtoolkit is None


def test_toolkit_module_does_not_import_jaxtoolkit_eagerly():
    """Pin the specific regression: no unconditional _jaxtoolkit import.

    Checked as source text rather than behaviour because the behavioural test
    above passes for *either* reason -- guard working, or JAX installed -- and
    this is the line that made the difference.
    """
    import pycircuit.circuit.toolkit as tk
    with open(tk.__file__) as fh:
        source = fh.read()
    for line in source.split('\n'):
        stripped = line.strip()
        if stripped.startswith('#') or stripped.startswith('##'):
            continue
        assert stripped != 'from . import _jaxtoolkit', (
            'unconditional _jaxtoolkit import is back; it defeats the guarded '
            'import at the bottom of toolkit.py and makes JAX mandatory')


def test_elements_do_not_mention_jax():
    """Elements reach differentiation through the toolkit, not through JAX.

    Twenty elements used to carry ``if hasattr(toolkit, 'jax') ...: import jax``
    inline.  They now call ``toolkit.jacobian(...)`` or guard on
    ``supports('autodiff')``, so adding another differentiating backend does not
    mean editing the element layer again.
    """
    from pycircuit.circuit import elements
    with open(elements.__file__) as fh:
        source = fh.read()
    assert 'import jax' not in source
    assert "'jax'" not in source
