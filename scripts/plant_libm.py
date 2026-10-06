"""A SIMULATED CPU for the C backend's tests (speed round 12, 2026-10-06): a
pytest plugin that makes numpy's float64 `exp` and `log` (or the functions
named in `$PYCIRCUIT_PLANT_LIBM`, comma-separated) differ from the C
library's -- one ulp up or down on about a quarter of the arguments, normal
results only -- as numpy's own AVX-512 loops do on some CPUs.

    PYTHONPATH=scripts pytest -p plant_libm pycircuit --tier fast

`_hdl_cbackend.libm_check` then lists them (exp ~73, log ~68 of 960
probes), and every test that compares C with numpy must take the rule for
listed functions (`_hdl_cbackend`, "Fidelity").  Results that are not
float64 (complex, integer) pass through untouched.  jax is imported first:
its `ml_dtypes` registers loops on the real ufuncs at import.  Loaded
before pycircuit's circuit modules, which bind numpy's functions by name.
"""
import os

import numpy as np

try:
    import jax  # noqa: F401
except ImportError:
    pass

_ORIG = {}


def _make(name):
    f = getattr(np, name)
    _ORIG[name] = f

    def planted(x, *a, **k):
        out = f(x, *a, **k)
        if a or k or np.asarray(out).dtype != np.float64:
            return out
        with np.errstate(all='ignore'):
            o = np.asarray(out, dtype=float)
            bits = np.asarray(np.asarray(x, dtype=float)).view(np.int64)
            ## about a quarter of the arguments, by their low bits; up or
            ## down by the next one
            sel = ((bits & 3) == 3) & np.isfinite(o) & (np.abs(o) >= 2.2250738585072014e-308)
            r = np.where(sel, np.nextafter(o, np.where((bits & 4) == 4, np.inf, -np.inf)), o)
        return r[()] if np.ndim(out) == 0 else r
    planted.__name__ = name
    return planted


for _n in os.environ.get('PYCIRCUIT_PLANT_LIBM', 'exp,log').split(','):
    if _n.strip():
        setattr(np, _n.strip(), _make(_n.strip()))
