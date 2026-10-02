"""The bit-identical common-subexpression elimination of the generated
chain functions (`pycircuit/circuit/_hdl_cse.py`; 2026-10-02, the
fused-evaluation plan's F1a): the rules it hoists and refuses by, the
verification it rests on, every chained library class byte for byte, and
the switches around it."""
import ast
import os
import subprocess
import sys
import textwrap

import numpy as np
import pytest

from pycircuit.circuit import _hdl_cse as cs
from pycircuit.circuit import circuit as cm
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import hdl
from pycircuit.circuit.circuit import defaultepar
from pycircuit.circuit.toolkit import numeric


def _fn(text):
    ns = {'numpy': np}
    exec(compile(text, '<t>', 'exec'), ns)  # noqa: S102 -- a test chain
    return ns['_f']


def _calls(text):
    return sum(isinstance(n, ast.Call) for n in ast.walk(ast.parse(text)))


def test_a_repeated_subtree_is_computed_once_and_a_whole_value_is_reused():
    src = textwrap.dedent('''\
        def _f(x):
            _x0 = x[0]
            _v1 = numpy.exp(_x0)
            _v2 = numpy.exp(_x0)*numpy.sqrt(_x0 + 1.0)
            _d_1 = numpy.sqrt(_x0 + 1.0) - 2.0*numpy.exp(_x0)
            return [_v1, _v2, _d_1]''')
    out, st = cs.cse_source(src)
    ## exp(_x0) is `_v1`'s whole value: reused, not recomputed; the sqrt
    ## appears twice: one temporary
    assert out.count('numpy.exp') == 1 and out.count('numpy.sqrt') == 1
    assert st['temporaries'] == 1
    x = np.array([0.37])
    a, b = np.asarray(_fn(src)(x)), np.asarray(_fn(out)(x))
    assert a.tobytes() == b.tobytes()


def test_the_pass_refuses_what_it_cannot_keep_exact():
    for body in ('_v1 = _x0 if _x0 > 0 else 1.0',          # lazy arm
                 '_v1 = _x0 and 1.0',                        # short circuit
                 '_v1 = sum(t for t in (_x0, 1.0))'):        # generator
        src = f'def _f(x):\n    _x0 = x[0]\n    {body}\n    return [_v1]'
        with pytest.raises(cs.Refused):
            cs.cse_source(src)
    ## and a name assigned twice is not SSA
    with pytest.raises(cs.Refused):
        cs.cse_source('def _f(x):\n    _a = x[0]\n    _a = _a + 1.0\n'
                      '    return [_a]')


def test_constants_are_keyed_by_type_and_sign_and_never_hoisted_alone():
    src = textwrap.dedent('''\
        def _f(x):
            _x0 = x[0]
            _v1 = _x0*0
            _v2 = _x0*0.0
            _v3 = _x0*-0.0
            _v4 = (2.0 + 3.0)*_x0 + (2.0 + 3.0)*_x0
            return [_v1, _v2, _v3, _v4]''')
    out, st = cs.cse_source(src)
    ## `_x0*0`, `_x0*0.0` and `_x0*-0.0` are three different products (an
    ## integer zero and a negative zero round differently); the constant
    ## `(2.0 + 3.0)` is not a temporary, its product with `_x0` is
    assert st['temporaries'] == 1
    for v in (-1.5, 0.0, -0.0, 2.0):
        x = np.array([v])
        assert (np.asarray(_fn(src)(x)).tobytes()
                == np.asarray(_fn(out)(x)).tobytes())


def test_a_deep_chain_is_handled_without_recursion():
    ## PSP's longest statement nests ~300 operations; the walk is iterative
    ## and the text is spliced, never unparsed
    expr = ' + '.join(f'numpy.sin(_x0*{k % 7}.0)' for k in range(600))
    src = f'def _f(x):\n    _x0 = x[0]\n    _v1 = {expr}\n    return [_v1]'
    out, st = cs.cse_source(src)
    assert st['temporaries'] == 7
    x = np.array([0.3])
    assert (np.asarray(_fn(src)(x)).tobytes()
            == np.asarray(_fn(out)(x)).tobytes())


def test_the_verification_rejects_a_changed_expression():
    src = ('def _f(x):\n    _x0 = x[0]\n    _v1 = (_x0 + 1.0) + 2.0\n'
           '    return [_v1]')
    wrong = ('def _f(x):\n    _x0 = x[0]\n    _v1 = _x0 + (1.0 + 2.0)\n'
             '    return [_v1]')
    fn = ast.parse(src).body[0]
    cs._verify(fn, src)
    with pytest.raises(AssertionError):
        cs._verify(fn, wrong)


def _chained():
    import inspect
    return sorted(n for n, c in vars(eh).items()
                  if inspect.isclass(c) and issubclass(c, hdl.Behavioural)
                  and c.__module__ == eh.__name__
                  and c._hdl_info.get('chained'))


@pytest.mark.parametrize('name', _chained())
def test_every_chained_library_class_is_byte_identical_to_its_reference(name):
    cm.default_toolkit = numeric
    cls = getattr(eh, name)
    e = cls(*[cm.Node(f'n{k}') for k in range(len(cls.terminals))])
    e.update_iparv()
    args = list(hdl._args_of(e, defaultepar))
    ## (the instance's own class: a model whose parameters collapse a node
    ## is built as a variant with its own, smaller, functions)
    funcs = type(e)._hdl_info['funcs']
    rng = np.random.default_rng(0)
    pts = [np.ascontiguousarray(p) for p in rng.uniform(-2, 2, (60, e.n))]
    pts += [np.full(e.n, v) for v in (1e30, -1e30, 100.0, -100.0, 0.7,
                                      0.0, -0.0)]
    for k in cs.FUNCS:
        f = funcs.get(k)
        if f is None:
            continue
        assert '_hdl_ref' in f.__dict__, (name, k)
        ref = f._hdl_ref
        ## never more calls emitted than the reference
        assert _calls(f._src_cse) <= _calls(ref._src), (name, k)
        assert f._src == ref._src           # the reference text is kept
        for x in pts:
            with np.errstate(all='ignore'):
                a = np.asarray(ref(x, *args), float)
                b = np.asarray(f(x, *args), float)
            assert a.tobytes() == b.tobytes(), (name, k, x)
        ## raise-iff-raise: the same operations, so the same conditions
        for x in pts[-7:]:
            outcome = []
            for g in (ref, f):
                try:
                    with np.errstate(all='raise'):
                        g(x, *args)
                    outcome.append(None)
                except FloatingPointError:
                    outcome.append('raised')
            assert outcome[0] == outcome[1], (name, k, x)


def test_the_auto_options_read_the_reference_size():
    """`compiled_jacobian_size` is calibrated on the reference bytecode
    (`AUTO_JACOBIAN_CODE`); the optimised `G` is a third of PSP's, and would
    flip the 'auto' Newton options without the pin (`_hdl_codelen`)."""
    from pycircuit.circuit._tran_newton import compiled_jacobian_size
    from pycircuit.circuit.compact import PspMosLongChannel
    cm.default_toolkit = numeric
    c = cm.SubCircuit()
    c['M'] = PspMosLongChannel('d', 'g', cm.gnd, cm.gnd)
    funcs = type(c['M'])._hdl_info['funcs']
    ref = sum(len(funcs[k]._hdl_ref.__code__.co_code) for k in ('G', 'C'))
    assert compiled_jacobian_size(c) == ref
    assert len(funcs['G'].__code__.co_code) < 0.5 * len(
        funcs['G']._hdl_ref.__code__.co_code)


def test_the_store_serves_the_text_without_redoing_the_pass(tmp_path,
                                                             monkeypatch):
    monkeypatch.setenv('PYCIRCUIT_HDL_CACHE_DIR', str(tmp_path))
    src = eh.MosLevel1Hdl._hdl_info['funcs']['G']._src
    first = cs._optimised_text(src)
    files = os.listdir(tmp_path / 'cse')
    assert len(files) == 1 and files[0].endswith('.py')

    def boom(_src):
        raise AssertionError('the store was not read')
    monkeypatch.setattr(cs, 'cse_source', boom)
    assert cs._optimised_text(src) == first


def _run(code, env):
    out = subprocess.run([sys.executable, '-c', code], capture_output=True,
                         text=True, env=dict(os.environ, **env), timeout=600,
                         check=False)
    assert out.returncode == 0, out.stderr[-2000:]
    return out.stdout


def test_the_switch_turns_the_pass_off():
    code = ('from pycircuit.circuit import elements_hdl as eh\n'
            'print("_hdl_ref" in eh.MosLevel1Hdl._hdl_info["funcs"]["G"]'
            '.__dict__)')
    assert _run(code, {'PYCIRCUIT_HDL_CSE': '0'}).strip() == 'False'
    assert _run(code, {}).strip() == 'True'


def test_the_text_does_not_depend_on_the_hash_seed():
    code = ('from pycircuit.circuit import _hdl_cse as cs, elements_hdl as eh\n'
            'import hashlib\n'
            'src = eh.MosLevel3Hdl._hdl_info["funcs"]["G"]._src\n'
            'print(hashlib.sha256(cs.cse_source(src)[0].encode()).hexdigest())')
    a = _run(code, {'PYTHONHASHSEED': '1'})
    b = _run(code, {'PYTHONHASHSEED': '12345'})
    assert a == b


def test_a_transient_is_byte_identical_with_the_reference_functions(
        monkeypatch):
    """End to end: a MosLevel1 inverter's transient with the optimised
    functions and with the references swapped back in."""
    from pycircuit.circuit.elements import VS, C, R, VPulse
    from pycircuit.circuit.transient import Transient
    cm.default_toolkit = numeric

    def inverter():
        c = cm.SubCircuit()
        c['vdd'] = VS('vdd', cm.gnd, v=1.8)
        c['vin'] = VPulse('in', cm.gnd, v1=0.0, v2=1.8, td=1e-9, tr=1e-10,
                          tf=1e-10, pw=4e-9, per=1e-8)
        c['M'] = eh.MosLevel1Hdl('out', 'in', cm.gnd, cm.gnd)
        c['RL'] = R('vdd', 'out', r=5e3)
        c['CL'] = C('out', cm.gnd, c=1e-13)
        return c
    a = Transient(inverter()).solve(tend=2e-8, timestep=1e-10)
    funcs = eh.MosLevel1Hdl._hdl_info['funcs']
    for k in cs.FUNCS:
        f = funcs.get(k)
        if f is not None and '_hdl_ref' in f.__dict__:
            monkeypatch.setitem(funcs, k, f._hdl_ref)
    b = Transient(inverter()).solve(tend=2e-8, timestep=1e-10)
    assert (np.asarray(a.x, float).tobytes()
            == np.asarray(b.x, float).tobytes())
