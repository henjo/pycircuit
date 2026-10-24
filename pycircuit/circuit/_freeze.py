"""THE START-UP'S OBJECTS FROZEN FOR THE COLLECTOR -- OPT-IN.

`freeze_imports()` is `gc.freeze()`: every object the garbage collector
tracks now moves to its permanent generation, which no collection walks
again.  The interpreter's exit runs several full collections as it tears
the modules down, each walking the ~146 000 objects a SPICE run's imports
leave -- 0.15 s of a short deck's ~0.8 s process (2026-10-25); frozen, the
exit takes ~0.025 s.  No number a simulation computes depends on it: what a
run makes is made after the freeze and collected as before.

Off unless asked -- `freeze_imports()` called, or `PYCIRCUIT_GC_FREEZE=1`
(the package's and the device library's imports then end by calling it) --
because it is a setting of the whole process, not of this package:

* everything alive at that moment is frozen, objects a program made before
  the import included: one that later becomes garbage in a reference cycle
  is never reclaimed;
* a module reloaded (`importlib.reload`, a notebook's autoreload) leaves
  its old version frozen and never reclaimed: memory grows with each reload;
* at exit, cyclic garbage among the frozen gets no collection, so its
  `__del__` and weakref callbacks do not run (Python promises neither at
  exit; without the freeze they run);
* code managing the collector itself (`gc.unfreeze()`, a pre-fork server's
  own freeze, a memory profiler reading the generations) finds a state it
  did not set.  `gc.unfreeze()` undoes it.
"""
import gc
import os


def freeze_imports():
    """`gc.freeze()` -- read this module's note before calling it."""
    gc.freeze()


def freeze_if_asked():
    """`freeze_imports()` where the environment has `PYCIRCUIT_GC_FREEZE=1`."""
    if os.environ.get('PYCIRCUIT_GC_FREEZE', '').strip() == '1':
        freeze_imports()
