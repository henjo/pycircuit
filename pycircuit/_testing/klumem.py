"""KLU's own count of the memory its objects hold (`klu_common.memusage`), for tests and the gate.

A numeric made and never freed changes no result and no count a comparison reads, so a leak is invisible to every twin check.  KLU counts what its objects hold whichever caller made them -- the C units through the pointers their glue hands them, the Python
through ctypes -- in the common block both pass.  The binding never reads that block (`linearsolver.KLUSolver`: no layout it cannot verify), so the offset here
comes from the compiler on the installed header, and is trusted only after KLU's own behaviour confirms it: 0 after `klu_defaults`, more after a factor, back
after the frees.  None where it cannot be had (no header, no compiler, KLU absent, the check failing).
"""
import ctypes

_SRC = ('#include <stddef.h>\n#include <suitesparse/klu.h>\n'
        'long pcx_klumem(void) { return (long) offsetof(klu_common, memusage); }\n')
_OFFSET = []


def _offset():
    if _OFFSET:
        return _OFFSET[0]
    _OFFSET.append(None)
    try:
        from pycircuit.circuit import _cbuild
        _ffi, fn, _key, _cold, _s = _cbuild.load_unit(_SRC, 'long pcx_klumem(void);', 'pcx_klumem')
        off = int(fn())
    except Exception:       # noqa: BLE001 -- no header, compiler or KLU: unmeasured, not a failure
        return None
    if _confirmed(off):
        _OFFSET[0] = off
    return _OFFSET[0]


def _read(common, off):
    return ctypes.c_size_t.from_buffer(common, off).value


def _confirmed(off):
    """KLU's own behaviour at `off`: nothing held after `klu_defaults`, a factor's memory after a solve, nothing after its objects are freed."""
    import numpy as np
    try:
        from pycircuit.circuit.linearsolver import KLUSolver
        sk = KLUSolver()
    except ImportError:
        return False
    if off < 0 or off + ctypes.sizeof(ctypes.c_size_t) > len(sk._common) or _read(sk._common, off) != 0:
        return False
    a = np.array([[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]])
    sk.solve(a, np.ones(3), None)
    held = _read(sk._common, off)
    sk._free_numeric()
    sk._free_symbolic()
    return held > 0 and _read(sk._common, off) == 0


def memusage(solver):
    """The bytes KLU's objects of `solver` hold now (its `_common`), or None where unmeasured or `solver` holds no KLU."""
    common = getattr(getattr(solver, '_choice', solver), '_common', None)
    off = _offset()
    if common is None or off is None:
        return None
    return _read(common, off)
