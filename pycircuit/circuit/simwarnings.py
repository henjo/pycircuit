"""The warnings pycircuit's analyses give, by what they mean, and where they
are attributed.

Every one is a `SimulationWarning`, itself a `RuntimeWarning`, so a filter
or a `pytest.warns(RuntimeWarning)` written before these existed still sees
it:

* `AccuracyWarning` -- the answer is less accurate than was asked, or its
  accuracy is not known: truncation error over tolerance, an unresolved
  quadrature or fit, a grid that cannot resolve what is asked of it.
* `ConvergenceWarning` -- a solve did not converge, or fell back to another
  route (an iteration budget, a Newton fallback, PCNR's).
* `ModelWarning` -- the analysis represents the circuit or its noise only
  in part, or approximately (a sign-blind square root, white sources only,
  an index-2 row), or a check of the model could not run.
* `CostWarning` -- a costlier route was taken; nothing is less accurate.
* `UsageWarning` -- an argument or a setting has no effect, or a
  combination is unusual.
* `PlatformWarning` -- on this machine an answer agrees with its other path
  to an ulp, not bitwise (the C backend where the CPU's numpy brings its own
  `exp`, `log`, ...: `_hdl_cbackend.libm_check`).

⚠ ONE ATTRIBUTION FOR EVERY WARNING.  `warn` attributes the warning to the
first frame OUTSIDE the library (`warnings.warn`'s `skip_file_prefixes`),
however deep it is raised -- a fixed `stacklevel` is right for one call
chain only, and the same helper is reached from several (it landed inside
the library for ~15 of the 76 sites, and two had none; the review's X8,
2026-10-01).  The tests and benchmarks are outside the library.
"""
import contextlib
import os
import warnings

__all__ = ['AccuracyWarning', 'ConvergenceWarning', 'CostWarning',
           'ModelWarning', 'PlatformWarning', 'SimulationWarning', 'UsageWarning',
           'summarised', 'warn']


class SimulationWarning(RuntimeWarning):
    """A warning from one of pycircuit's analyses."""


class AccuracyWarning(SimulationWarning):
    """The answer is less accurate than was asked, or its accuracy is not
    known."""


class ConvergenceWarning(SimulationWarning):
    """A solve did not converge, or fell back to another route."""


class ModelWarning(SimulationWarning):
    """The analysis represents the circuit or its noise only in part, or
    approximately -- or a check of the model could not run."""


class CostWarning(SimulationWarning):
    """A costlier route was taken; nothing is less accurate."""


class UsageWarning(SimulationWarning):
    """An argument or a setting has no effect, or a combination is
    unusual."""


class PlatformWarning(SimulationWarning):
    """On this machine an answer agrees with its other path to an ulp, not
    bitwise: nothing is less accurate, but a contract measured elsewhere as
    bit-identity holds here only to the last bit."""


def _library_files():
    """Every module of the package but its tests, as absolute paths."""
    pkg = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    out = []
    for dirpath, dirnames, filenames in os.walk(pkg):
        dirnames[:] = [d for d in dirnames
                       if d not in ('tests', '__pycache__')]
        out += [os.path.join(dirpath, f) for f in filenames
                if f.endswith('.py')]
    return tuple(sorted(out))


_LIBRARY = _library_files()


def warn(message, category=SimulationWarning):
    """`warnings.warn` attributed to the first frame outside the library."""
    warnings.warn(message, category, skip_file_prefixes=_LIBRARY)



@contextlib.contextmanager
def summarised(what, categories=(ConvergenceWarning, AccuracyWarning)):
    """Run a SUB-solve the answer is built on, and give its warnings as ONE
    per category that bears on the answer -- how many, and the first --
    rather than none (a blanket filter hid the sub-solve's own accuracy from
    the answer it defines: the monodromy twin, `lte_grid`'s run,
    `warping_estimate`'s; the review's F13, 2026-10-01) or each one,
    attributed inside the library.  Anything else it warns (a cost, a usage
    note, numpy's) is dropped, as before.  On an exception nothing is
    given: the caller's handling of it says what happened."""
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        yield
    seen = {}
    for r in rec:
        for cat in categories:
            if issubclass(r.category, cat):
                seen.setdefault(cat, []).append(str(r.message))
                break
    for cat in categories:
        msgs = seen.get(cat)
        if msgs:
            first = msgs[0] if len(msgs[0]) <= 600 else msgs[0][:600] + ' ...'
            warn(f'{what} warned {len(msgs)} time(s) ({cat.__name__}); the '
                 f'first: {first}', cat)
