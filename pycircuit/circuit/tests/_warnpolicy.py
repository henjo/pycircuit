"""The suite's warnings policy (the review's test policy, 2026-10-01).

`pytest.ini` makes every `SimulationWarning` -- pycircuit's own warnings,
`pycircuit.circuit.simwarnings` -- an ERROR, so a library warning a test
does not expect fails it.  A test that expects some says which, by
category:

    with quiet(AccuracyWarning, ConvergenceWarning):
        pss.solve(...)

ignores those two and anything that is not pycircuit's (numpy's overflow
and the like, as the blanket ``simplefilter('ignore')`` it replaces did);
any OTHER pycircuit warning raised inside is an error.  Until 2026-10-01
~800 blocks ignored everything, so a new library warning was invisible to
the suite.

``PYCIRCUIT_WARN_SURVEY=<dir>``: instead of erroring, every block records
the pycircuit categories raised inside it (one JSON line per block per
category, a file per process) -- how the blocks' categories were filled in.
"""
import json
import os
import sys
import warnings

from pycircuit.circuit.simwarnings import SimulationWarning

_SURVEY = os.environ.get('PYCIRCUIT_WARN_SURVEY')


class quiet:
    """``with quiet(*categories):`` -- see the module note."""

    def __init__(self, *expected):
        self.expected = expected
        f = sys._getframe(1)
        self.site = (os.path.abspath(f.f_code.co_filename), f.f_lineno)

    def __enter__(self):
        self._cm = warnings.catch_warnings(record=bool(_SURVEY))
        rec = self._cm.__enter__()
        self._rec = rec
        warnings.simplefilter('ignore')
        if _SURVEY:
            warnings.simplefilter('always', SimulationWarning)
        else:
            warnings.simplefilter('error', SimulationWarning)
            for cat in self.expected:
                warnings.simplefilter('ignore', cat)
        return self

    def __exit__(self, *exc):
        rec = self._rec
        out = self._cm.__exit__(*exc)
        if _SURVEY and rec:
            cats = sorted({r.category.__name__ for r in rec
                           if issubclass(r.category, SimulationWarning)})
            if cats:
                os.makedirs(_SURVEY, exist_ok=True)
                with open(os.path.join(_SURVEY, f'survey-{os.getpid()}.jsonl'),
                          'a') as fh:
                    fh.write(json.dumps({'file': self.site[0],
                                         'line': self.site[1],
                                         'cats': cats}) + '\n')
        return out
