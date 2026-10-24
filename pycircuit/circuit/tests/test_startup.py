"""What a numeric start imports, and what it leaves for first use.

A SPICE transient's process imports the package, the netlist importer and the
device library.  JAX (the JAX toolkit's, ~0.2 s) and scipy's special,
interpolate and optimize (the post-processing's and the table function's,
~0.2 s) are imported where they are used, never by that start.  JAX imported
later -- by the JAX toolkit's first use, or by any module directly -- runs in
64-bit mode, as every use of it in the package assumes, and the JAX toolkit
is one instance whichever thread asks first.  Each in a fresh interpreter
(this one imported everything long ago), importing this tree, without the
JAX switch this process's toolkit set."""
import os
import subprocess
import sys
import textwrap

import pytest

import pycircuit

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(pycircuit.__file__)))


def _child(script):
    """`script`'s printed words, from a fresh interpreter."""
    env = {k: v for k, v in os.environ.items() if k != 'JAX_ENABLE_X64'}
    out = subprocess.run([sys.executable, '-c', textwrap.dedent(script)], cwd=_ROOT, env=env,
                         capture_output=True, text=True, timeout=600)
    assert out.returncode == 0, out.stderr[-3000:]
    return out.stdout.split()


def test_a_numeric_start_imports_neither_jax_nor_scipys_extras():
    """The package, the netlist importer and the device library imported:
    no JAX, no scipy.special, scipy.interpolate or scipy.optimize."""
    got = _child('''
        import sys
        import pycircuit.circuit
        from pycircuit.circuit import elements_hdl, spice_import
        assert elements_hdl.DiodeHdl and spice_import.import_netlist
        print(*[m for m in ('jax', 'scipy.special', 'scipy.interpolate', 'scipy.optimize')
                if m in sys.modules] or ['none'])
    ''')
    assert got == ['none']


def test_jax_imported_directly_later_is_64_bit():
    """JAX imported by a module of its own after the package: 64-bit."""
    pytest.importorskip('jax')
    got = _child('''
        import sys
        import pycircuit.circuit
        assert 'jax' not in sys.modules
        import jax
        import jax.numpy as jnp
        print(jax.config.jax_enable_x64, jnp.asarray(1.0).dtype)
    ''')
    assert got == ['True', 'float64']


def test_the_jax_toolkit_asked_for_at_once_is_one():
    """Eight threads asking for the JAX toolkit at once, JAX not yet
    imported: one instance (each thread's kept, so no address is reused), in
    64-bit mode."""
    pytest.importorskip('jax')
    got = _child('''
        import sys, threading
        import pycircuit.circuit.toolkit as tk
        assert 'jax' not in sys.modules
        start, seen = threading.Barrier(8), []
        def ask():
            start.wait()
            seen.append(tk.jaxtoolkit)
        ts = [threading.Thread(target=ask) for _ in range(8)]
        for t in ts:
            t.start()
        for t in ts:
            t.join()
        import jax.numpy as jnp
        from pycircuit.circuit.toolkit import jaxtoolkit
        print(len({id(s) for s in seen}), len(seen), all(s is jaxtoolkit for s in seen),
              type(jaxtoolkit).__name__, jnp.asarray(1.0).dtype)
    ''')
    assert got == ['1', '8', 'True', 'JAXToolkit', 'float64']
