"""`$limit` and `$param_given` in one model used to be a run-time TypeError.

Every compiled function in `hdl.py` is called through `_args_of`, which
is ``params + [T] + givenness flags``.  `limit_spec`'s parameter
expressions alone were lambdified over ``paramsyms + [TEMP]``, so a
model declaring both features compiled clean and then died on the first
Newton iteration with

    TypeError: _lambdifygenerated() takes N positional arguments
               but N+1 were given

from inside the generated `limit()` -- a place no author would look for
a parameter-passing mistake.

Nothing had ever used the two together.  `$limit` arrived with the
junction devices, `$param_given` with the SPICE-compatibility work, and
the first model to want both was the Gummel-Poon BJT: SPICE's ``RBM``
defaults to ``RB``, which a default VALUE cannot express, so it needs
to ask whether the parameter was given at all.

The trap when writing this test: `0.0*param_given('x')` is simplified
to `0` by sympy before the compiler ever sees it, `given_names` comes
out empty, and the bug does NOT reproduce.  The flag has to be
load-bearing.  `test_the_flag_is_load_bearing` guards that, because a
test that silently stops exercising the bug is worse than no test.
"""

import numpy as np
import pytest

import pycircuit.circuit.circuit as cm
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   limit_pnj, limexp, param_given)
from pycircuit.utilities.param import Parameter


class LimitAndGiven(Behavioural):
    """A junction that also asks whether `scale` was given."""

    instparams = [Parameter(name='iss', desc='IS', unit='A', default=1e-14),
                  Parameter(name='vtt', desc='VT', unit='V', default=0.02585),
                  Parameter(name='scale', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        ## Load-bearing: selects between the given value and a default of 2.
        s = (param_given('scale') * scale                      # noqa: F821
             + (1 - param_given('scale')) * 2.0)               # noqa: F821
        return Contribution(
            b.I, s * iss * (limexp(limit_pnj(b.V, iss, vtt)    # noqa: F821
                                   / vtt) - 1))                # noqa: F821


@pytest.fixture(autouse=True)
def _numeric():
    cm.default_toolkit = numeric


def test_the_flag_is_load_bearing():
    """Without this the rest of the file can pass while testing nothing."""
    assert LimitAndGiven._hdl_info['given_names'] == ['scale'], \
        'param_given was optimised away; this file no longer tests the bug'
    assert len(LimitAndGiven._hdl_info['limit_spec']) == 1


def test_limit_runs_on_a_model_that_also_uses_param_given():
    """The regression.  Against the narrow signature this is a TypeError."""
    e = LimitAndGiven(cm.Node('p'), cm.Node('n'))
    e.update_iparv()
    out = e.limit(np.array([0.9, 0.0]), np.array([0.1, 0.0]))
    assert np.all(np.isfinite(out))
    ## It must actually limit: 0.9 V from a 0.1 V iterate is a wild step
    ## for a junction, so the result has to land well below 0.9.
    assert out[1] == 0.0, 'the cathode must not move'
    assert 0.1 < out[0] < 0.9, out


def test_the_limiter_still_agrees_with_the_bare_kernel():
    """Widening the signature must not change what the limiter COMPUTES."""
    from pycircuit.circuit._limiting import _pnjlim
    e = LimitAndGiven(cm.Node('p'), cm.Node('n'))
    e.update_iparv()
    for vnew, vold in ((0.9, 0.1), (0.75, 0.7), (-2.0, 0.4), (0.4, 0.4)):
        got = e.limit(np.array([vnew, 0.0]), np.array([vold, 0.0]))[0]
        want = _pnjlim(vnew, vold, 0.02585, 1e-14, e.toolkit)
        assert got == pytest.approx(want, rel=1e-15, abs=1e-18)


@pytest.mark.parametrize('given', [True, False])
def test_the_flag_still_reaches_the_current(given):
    """The other half: widening the signature must not break `param_given`.

    `scale` given -> the current scales by it; not given -> by 2.0.
    """
    kw = dict(scale=5.0) if given else {}
    e = LimitAndGiven(cm.Node('p'), cm.Node('n'), **kw)
    e.update_iparv()
    ref = LimitAndGiven(cm.Node('p'), cm.Node('n'), scale=1.0)
    ref.update_iparv()
    x = np.array([0.6, 0.0])
    ratio = float(np.asarray(e.i(x))[0]) / float(np.asarray(ref.i(x))[0])
    assert ratio == pytest.approx(5.0 if given else 2.0, rel=1e-12)
