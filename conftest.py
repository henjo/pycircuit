"""Session-wide pytest configuration.

Kept free of ``jax`` imports at module level: the environment variable below
only takes effect if it is set before JAX initialises its backend, and this
file is imported before any test module.
"""

import os

# JAX preallocates ~75% of the GPU's memory in each process the first time a
# device is used.  Under pytest-xdist every worker is a separate process, so on
# a GPU with modest VRAM all but the first worker fails to obtain memory --
# and the resulting failures surface as ordinary assertion errors rather than
# an out-of-memory message, which makes them easy to misread as numerical bugs
# in the solver.  Allocate on demand instead.
#
# Measured on an RTX A1000 (4 GB) with ``-n 8``: 15 failures with
# preallocation, 0 without.  Harmless on CPU-only installations.
#
# ``setdefault`` so an explicit setting in the environment still wins.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")

# Some example tests call ``pylab.show()``.  On a developer machine with DISPLAY
# set, matplotlib picks an interactive backend (tkagg) and show() opens a window
# and blocks until something closes it -- so the suite appears to hang, for
# minutes at a time, and nondeterministically depending on the window manager.
# test_my_elements.py::test_transient_plot took 148 s this way against 0.96 s
# headless.  Tests should never open a window; force the non-interactive
# backend.  Again setdefault, so running with MPLBACKEND=tkagg to eyeball a
# plot still works.
os.environ.setdefault("MPLBACKEND", "Agg")
