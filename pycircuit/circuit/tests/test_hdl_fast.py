"""The exact scalar fast paths of the numpy twins (`_hdl_fast`, the speed plan
round 2, stage C).

Each helper must answer what its numpy counterpart answers -- the same type
(a numpy SCALAR where numpy returns one; `_fwhere` the scalar of the dtype
numpy's 0-d array holds), the same bits, the same raise-mode behaviour --
for every input, fast branch or not.  The two places a numpy scalar and a
0-d array compute differently are pinned here as live (so the guards are
known to be needed), and the guards are tested: a `where` value reaching
`**` keeps `numpy.where`; a non-finite chosen value keeps it too.
"""
import itertools
import math
import warnings

import numpy as np
import pytest
import sympy

from pycircuit.circuit import _hdl_cse as cs
from pycircuit.circuit import _hdl_fast as hf
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.hdl import Behavioural, Branch, Contribution, var
from pycircuit.circuit.toolkit import numeric
from pycircuit.utilities.param import Parameter

H = hf.helpers()
NAN_NEG = -float('nan')
VALUES = [0.0, -0.0, 0.3, -0.3, 1.0, 1e308, -1e308, 5e-324, -5e-324,
          float('inf'), -float('inf'), float('nan'), NAN_NEG, 2.0 ** 53,
          2.0 ** 53 + 2]
INTS = [0, 1, -1, 2 ** 53, 2 ** 53 + 1, -(2 ** 60)]


def _forms(v):
    """`v` as every type the generated code can hand a helper."""
    out = [v, np.float64(v), np.array(v), np.array([v])]
    if isinstance(v, float) and v.is_integer() and abs(v) < 2 ** 62:
        out.append(np.int64(int(v)))
    return out


def _outcome(fn, *a):
    """(kind, type, dtype, bytes) under raise mode, plus the warnings the
    same call gives under warn mode."""
    try:
        with np.errstate(all='raise'):
            r = fn(*a)
    except FloatingPointError as e:
        res = ('raise', str(e).split(' in ')[0])
    except (TypeError, ValueError) as e:
        res = ('error', type(e).__name__)
    else:
        arr = np.asarray(r)
        res = ('ok', type(r).__name__, arr.dtype.str, arr.tobytes())
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        with np.errstate(all='warn'):
            try:
                fn(*a)
            except (TypeError, ValueError):
                pass
    return res, len(w)


_PAIRS = [(a, b) for x, y in itertools.product(VALUES, VALUES)
          for a in _forms(x) for b in _forms(y)][::3]
_PAIRS += [(i, j) for i in INTS for j in INTS + VALUES[:6]]
_PAIRS += [(i, np.float64(j)) for i in INTS for j in VALUES[:6]]

_BINARY = [('_fmaxc', hdl._maxc_numpy), ('_fminc', hdl._minc_numpy),
           ('_fmax', np.maximum), ('_fmin', np.minimum),
           ('_fless', np.less), ('_fless_equal', np.less_equal),
           ('_fgreater', np.greater),
           ('_fgreater_equal', np.greater_equal),
           ('_fequal', np.equal), ('_fnot_equal', np.not_equal),
           ('_fstep', hdl._step_numpy)]


@pytest.mark.parametrize('name,ref', _BINARY, ids=[n for n, _ in _BINARY])
def test_a_binary_fast_path_answers_as_its_numpy_counterpart(name, ref):
    """Type, dtype, bytes, raise-mode outcome and warning count, over
    signed zeros, ties, infinities, NaNs of both signs, subnormals, ints
    beyond 2**53, in every type the generated code passes."""
    fast = H[name]
    for a, b in _PAIRS:
        assert _outcome(fast, a, b) == _outcome(ref, a, b), (name, a, b)


def test_the_fast_maxc_keeps_the_symbolic_guard():
    x = sympy.Symbol('x')
    with pytest.raises(TypeError, match='maxc has no symbolic evaluation'):
        H['_fmaxc'](x, 1.0)


def test_the_fast_where_answers_as_numpy_where():
    """The value (bytes and dtype) of `numpy.where`'s 0-d array, as a
    scalar where the chosen arm is a finite float; numpy's own answer
    otherwise -- int arms, a 0-d or array argument, a non-finite choice."""
    conds = [True, False, np.True_, np.False_, np.float64(0.0),
             np.float64(np.nan), np.array(True), np.array([True, False])]
    arms = [v for x in VALUES for v in _forms(x)[:3]] + [1, 0, np.int64(3)]
    for c in conds:
        for a, b in itertools.product(arms[::2], arms[1::2]):
            (kr, *rr), wr = _outcome(np.where, c, a, b)
            (kf, *rf), wf = _outcome(H['_fwhere'], c, a, b)
            assert (kr, wr) == (kf, wf), (c, a, b)
            if kr == 'ok':
                assert rr[1:] == rf[1:], (c, a, b)        # dtype, bytes
                fast = rf[0] != 'ndarray'
                chosen = np.asarray(np.where(c, a, b))
                if fast:
                    assert chosen.dtype == np.float64
                    assert math.isfinite(float(chosen))


def test_a_numpy_scalar_and_a_0d_array_part_ways_only_where_guarded():
    """The guards' reason, kept alive: `**` with exponent 2, 1/2 or -1 on a
    0-d array differs from the scalar's glibc `pow` at some points, and an
    operation meeting two NaNs returns another one; every other operation
    the twins apply to a FINITE where value agrees bit for bit."""
    rng = np.random.default_rng(1)
    xs = rng.uniform(0.1, 10.0, 4000)
    for e in (2, 0.5, -1):
        assert any((np.array(x) ** e).tobytes()
                   != (np.float64(x) ** e).tobytes() for x in xs), e
    a, b = np.float64(np.nan), np.float64(NAN_NEG)
    assert (np.array(a) + b).tobytes() != (a + b).tobytes() or \
        (b + np.array(a)).tobytes() != (b + a).tobytes()
    others = [v for x in VALUES for v in (x, np.float64(x))]
    ops = [lambda p, q: p + q, lambda p, q: p - q, lambda p, q: p * q,
           lambda p, q: p / q, lambda p, q: q + p, lambda p, q: q - p,
           lambda p, q: q * p, lambda p, q: q / p, lambda p, q: -p,
           lambda p, q: abs(p), lambda p, q: np.sqrt(p),
           lambda p, q: np.exp(p), lambda p, q: np.log(p),
           lambda p, q: np.sign(p), lambda p, q: p ** 3,
           lambda p, q: p ** 2.5]
    for v in xs[:40].tolist() + [0.0, -0.0, 1e-300, -2.5, 1e300]:
        for q in others:
            for op in ops:
                with np.errstate(all='ignore'):
                    r0 = np.asarray(op(np.array(v), q))
                    r1 = np.asarray(op(np.float64(v), q))
                assert r0.tobytes() == r1.tobytes(), (v, q)


class _SquaredWhere(Behavioural):
    instparams = [Parameter(name='gg', desc='g', unit='S', default=1.0)]

    @staticmethod
    def analog(p, m):
        b = Branch(p, m)
        u = var(b.V, 'u')
        w = var(sympy.Piecewise((u + 1.0, u > 0), (1.0 - u, True)), 'w')
        w2 = var(w, 'w2')
        return Contribution(b.I, gg * (w2 ** 2 + w))          # noqa: F821


def test_a_where_value_that_is_squared_keeps_numpy_where():
    """The `**` guard, through a plain rename: the squared value's `where`
    stays `numpy.where` (a 0-d array, squared exactly), the others go fast,
    and the twin returns the reference's bytes where `pow` and the square
    differ."""
    fn = _SquaredWhere._hdl_info['funcs']['i']
    text = fn._src_cse
    tree = cs.ast.parse(text).body[0]
    assert cs._kept_wheres(tree), text
    assert 'numpy.where(' in text
    x0 = next(float(x) for x in np.random.default_rng(2).uniform(
        0.1, 3.0, 5000)
        if (np.array(x + 1.0) ** 2).tobytes()
        != (np.float64(x + 1.0) ** 2).tobytes())
    args = list(hdl._args_of(_inst(_SquaredWhere), defaultepar))
    x = np.array([x0, 0.0])
    assert np.asarray(fn(x, *args)).tobytes() == \
        np.asarray(fn._hdl_ref(x, *args)).tobytes()


def _inst(cls):
    cm.default_toolkit = numeric
    e = cls(*[Node(f'n{k}') for k in range(len(cls.terminals))])
    e.update_iparv()
    return e


def test_the_rewrite_inverts_and_refuses_a_shadowed_name():
    src = ('def _f(x, _fwhere):\n    _x0 = x[0]\n'
           '    _v1 = numpy.where(numpy.less(_x0, 0.0), _x0, _fwhere)\n'
           '    return [_v1]')
    with pytest.raises(cs.Refused):
        cs._fast_rewrite(src)
    ok = src.replace('_fwhere', 'p0')
    out = cs._fast_rewrite(ok)
    assert '_fwhere(_fless(' in out
    back = out
    for old, new in hf.RENAMES.items():
        back = back.replace(new + '(', old + '(')
    assert back == ok


def test_the_library_twins_run_the_fast_paths():
    from pycircuit.circuit import elements_hdl as eh
    f = eh.MosLevel1Hdl._hdl_info['funcs']['G']
    assert '_fwhere(' in f._src_cse and '_fstep(' in f._src_cse
    assert f.__globals__['_fwhere'] is H['_fwhere']
    assert 'exact scalar fast paths' in hdl.explain(eh.MosLevel1Hdl)
