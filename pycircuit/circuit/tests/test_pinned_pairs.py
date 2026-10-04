"""PINNED PAIRS (2026-10-03): a Python object that DEFINES a behaviour and
the C text, or the second Python path, that REPRODUCES it, recorded
together -- so that a change to either side fails here, naming the twin,
until the record is re-made.  A prompt at the moment of editing, which
the agreement tests cannot give.

What this proves and what it does not: a digest says that something
changed, never that the two sides still agree.  Agreement is the job of
the sweeps that run both paths on the same inputs and assert the bytes
equal (`test_hdl_batch`, `test_hdl_climit`, `test_hdl_limit_walk`,
`test_stamp_plan`), of the recorded gate over every solve in the suite,
and of the private comparison.  This test adds the one thing those lack:
the author is made to look at the twin BEFORE any arithmetic has had the
chance to drift.

THE RULE: whoever changes one side re-records the pair -- after checking
or changing the twin -- and the record line lands in the same commit, so
a reviewer sees the code, its twin and the record move together.  To
re-record:

    python pycircuit/circuit/tests/test_pinned_pairs.py

prints the `RECORD` table to paste below.  Comments are part of a digest
on purpose: a comment that explains the arithmetic is part of what the
twin mirrors.  Trailing whitespace is not.

Each side is a list of `module:qualname` names; a name resolving to a
string is C text, anything else is read with `inspect.getsource`.
"""
import hashlib
import importlib
import inspect
import sys
import textwrap
import types

import pytest

PAIRS = {
    'the kernel call and the pass driver': {
        'why': 'the driver writes the temperature into the element\'s own pack and '
            'calls the kernel exactly as `CKernel.__call__` does; `Batch.run` '
            'mirrors its pack and temperature rules per pass',
        'reference': ['pycircuit.circuit._hdl_cbackend:CKernel.__call__',
                   'pycircuit.circuit._hdl_cbackend:CKernel.pack'],
        'twin': ['pycircuit.circuit._hdl_batch:PASS_C',
              'pycircuit.circuit._hdl_batch:Batch'],
    },
    'the limiting loop and the walk': {
        'why': 'the walk is the loop in C: dict order, both gathers before the call, '
            'the write-back after, the last duplicate row winning, `x0` aliasing '
            '`x`; it stops where the loop\'s statement must run',
        'reference': ['pycircuit.circuit.circuit:SubCircuit.limit',
                   'pycircuit.circuit._hdl_climit:CLimitKernel.__call__'],
        'twin': ['pycircuit.circuit._hdl_climit:_WALK_C',
              'pycircuit.circuit._hdl_climit:_Walk',
              'pycircuit.circuit._hdl_climit:limit_walk'],
    },
    'the limiter laws and their C prelude': {
        'why': 'the prelude transliterates the laws and the device write-back in '
            'Python\'s own forms (max as `(b > a) ? b : a`, the tuple order on '
            'the keys, the stable sort); the renderer prints the closure\'s '
            'ranking and write-back from the spec',
        'reference': ['pycircuit.circuit._limiting:apply_limit',
                   'pycircuit.circuit._limiting:device_writeback'],
        'twin': ['pycircuit.circuit._hdl_climit:_LIMIT_C',
              'pycircuit.circuit._hdl_climit:render'],
    },
    'the transient evaluation and the core': {
        'why': 'the core evaluates the four passes, the companion in each '
               'integrator\'s own operation order and the residual in one C '
               'call, asks the same C lookup first, and leaves the state '
               '`_companion_at` and `get_diff` leave; `_solve_timestep`\'s '
               '`jacobian_only` and `residual_only` are its other two sites',
        'reference': ['pycircuit.circuit._tran_companion:_CompanionModel._residual_and_jacobian',
                      'pycircuit.circuit._tran_companion:_CompanionModel._companion_at',
                      'pycircuit.circuit._tran_companion:_CompanionModel._C_at_state',
                      'pycircuit.circuit._tran_companion:_CompanionModel.get_diff',
                      'pycircuit.circuit._lte_kernels:bdf2_companion',
                      'pycircuit.circuit._lte_kernels:euler_companion',
                      'pycircuit.circuit._lte_kernels:trapezoidal_companion',
                      'pycircuit.circuit._lte_kernels:theta_companion',
                      'pycircuit.circuit._stamp_plan:assemble_matrix',
                      'pycircuit.circuit._stamp_plan:assemble_vector',
                      'pycircuit.circuit.transient:Transient._solve_timestep'],
        'twin': ['pycircuit.circuit._tran_core:CORE_C',
                 'pycircuit.circuit._tran_core:_Core',
                 'pycircuit.circuit._tran_core:evaluate'],
    },
    'the assembly loops and the plan': {
        'why': 'the plan is the loop\'s bincount over the same values in the same '
            'order, with the constant elements pre-filled, the C-bound classes '
            'batched and the zero sources skipped; its fallbacks are the loop',
        'reference': ['pycircuit.circuit.circuit:SubCircuit._add_element_submatrices',
                   'pycircuit.circuit.circuit:SubCircuit._add_element_subvectors'],
        'twin': ['pycircuit.circuit._stamp_plan:assemble_matrix',
              'pycircuit.circuit._stamp_plan:assemble_vector',
              'pycircuit.circuit._stamp_plan:_run_batches',
              'pycircuit.circuit._stamp_plan:_legacy_matrix',
              'pycircuit.circuit._stamp_plan:_legacy_vector',
              'pycircuit.circuit._hdl_batch:split',
              'pycircuit.circuit._hdl_batch:zero_source'],
    },
    'the Newton solve and its C': {
        'why': 'the C solve is `nrsolver`\'s plain and chord Newton on the reduced '
               'system (numpy\'s and SciPy\'s own LAPACK, the walk, the test in its '
               'order, `dx` re-taken where a limiter exists), `jacobian_only`\'s '
               'converged-point evaluation, `_newton`\'s statistics and branch check',
        'reference': ['pycircuit.circuit.nrsolver:StandardNewton.solve_system',
                      'pycircuit.circuit.nrsolver:ChordNewton.solve_system',
                      'pycircuit.circuit.dcanalysis:refnode_removed',
                      'pycircuit.circuit.analysis:insert_row',
                      'pycircuit.circuit.analysis:remove_row_col',
                      'pycircuit.circuit._tran_newton:_StepNewton._newton',
                      'pycircuit.circuit._tran_newton:_StepNewton._newton_limiter',
                      'pycircuit.circuit.transient:Transient._solve_timestep'],
        'twin': ['pycircuit.circuit._tran_newton_c:NEWTON_C',
                 'pycircuit.circuit._tran_newton_c:solve',
                 'pycircuit.circuit._tran_newton_c:_Ctx',
                 'pycircuit.circuit._tran_newton_c:_walk_full'],
    },
    'the PSP limiter and its C twin': {
        'why': 'the twin is `PspMosLongChannel.limit` transliterated -- d, g, b moved '
               'about the source in Python\'s order, `(vs + vold) + (+-lim)` -- with '
               '`vlimit` baked in per value and any non-finite value declined',
        'reference': ['pycircuit.circuit.compact:PspMosLongChannel.limit'],
        'twin': ['pycircuit.circuit._hdl_climit:PSP_LIMIT_C',
                 'pycircuit.circuit._hdl_climit:_handwritten',
                 'pycircuit.circuit._hdl_climit:_Handwritten'],
    },
    'the error test and its C': {
        'why': 'the C error test is `_charge_lte` (the trapezoid\'s and Gear-2\'s '
               'charge-form `compute_lte`, numpy\'s own LAPACK on the reduced '
               '`J`), `_normalised` (the running reference by `relref`, '
               '`sigglobal_reference`, `tolerance`, `normalised_error`) and '
               '`np.max`, the running reference written only on success',
        'reference': ['pycircuit.circuit.stepcontroller:StepController._charge_lte',
                      'pycircuit.circuit.stepcontroller:StepController._normalised',
                      'pycircuit.circuit.stepcontroller:StepController._reference',
                      'pycircuit.circuit.stepcontroller:StepController.tolerance',
                      'pycircuit.circuit.stepcontroller:StepController._max_error',
                      'pycircuit.circuit.stepcontroller:normalised_error',
                      'pycircuit.circuit.stepcontroller:sigglobal_reference',
                      'pycircuit.circuit.integrator:Gear2Integrator.compute_lte',
                      'pycircuit.circuit.integrator:TrapezoidalIntegrator.compute_lte',
                      'pycircuit.circuit._lte_kernels:third_divided_difference',
                      'pycircuit.circuit.analysis:remove_row_col'],
        'twin': ['pycircuit.circuit._tran_lte_c:LTE_C',
                 'pycircuit.circuit._tran_lte_c:max_error',
                 'pycircuit.circuit._tran_lte_c:_globals',
                 'pycircuit.circuit._tran_lte_c:_Ctx'],
    },
    'the stage predictor and its multistep fast path': {
        'why': 'behind the target nearest-first is newest-first: the fast path '
               'keeps the general path\'s order, its 1e-13 dedupe, the fit\'s '
               'IEEE arithmetic and memo key bytes, and the same clip ufunc',
        'reference': ['pycircuit.circuit._tran_predictor:_StagePredictor._predict_state',
                      'pycircuit.circuit._tran_predictor:_StagePredictor._fit'],
        'twin': ['pycircuit.circuit._tran_predictor:_StagePredictor._predict_fast',
                 'pycircuit.circuit._tran_predictor:_StagePredictor._pred_fast_nodes',
                 'pycircuit.circuit._tran_predictor:_StagePredictor._fit_fast'],
    },
}

#: `pair name: (reference digest, twin digest)` -- re-made by running this
#: module as a script, after the twin was checked.
RECORD = {
    'the kernel call and the pass driver': ('09a9e274bbe8', 'c900976facbc'),
    'the limiting loop and the walk': ('d151206fce18', '411e5dd4e655'),
    'the limiter laws and their C prelude': ('9b7944f9c0dd', '770bc2e815cc'),
    'the transient evaluation and the core': ('b3fb392a7acd', '8433a00608ff'),
    'the assembly loops and the plan': ('ed84458db1a3', '9578c5eca9ad'),
    'the stage predictor and its multistep fast path': ('739b71ddbace', 'f2c284eb392d'),
    'the Newton solve and its C': ('7e6a12179239', '6904a806ad8d'),
    'the error test and its C': ('a24ba5659ca4', 'd6903a790f88'),
    'the PSP limiter and its C twin': ('8e14bbc903f1', '58969d2d34a6'),
}

#: The generated methods a batch and the walk tell from their doubles by
#: CODE identity (`_hdl_batch.generated_code`): each must be defined
#: exactly once in `BehaviouralMeta.__init__`, or the marker is ambiguous.
GENERATED = ('i', 'G', 'q', 'C', 'limit', 'u', 'dudt')


def _resolve(name):
    mod, _, qual = name.partition(':')
    obj = importlib.import_module(mod)
    for part in qual.split('.'):
        obj = getattr(obj, part)
    return obj


def _source(name):
    obj = _resolve(name)
    src = obj if isinstance(obj, str) else inspect.getsource(obj)
    lines = [ln.rstrip() for ln in textwrap.dedent(src).splitlines()]
    while lines and not lines[-1]:
        lines.pop()
    return '\n'.join(lines) + '\n'


def digest(names):
    h = hashlib.sha256()
    for name in names:
        h.update(name.encode())
        h.update(b'\0')
        h.update(_source(name).encode())
        h.update(b'\0')
    return h.hexdigest()[:12]


def current():
    return {name: (digest(p['reference']), digest(p['twin'])) for name, p in PAIRS.items()}


@pytest.mark.parametrize('name', list(PAIRS))
def test_every_side_resolves_to_source(name):
    for side in ('reference', 'twin'):
        for obj in PAIRS[name][side]:
            assert len(_source(obj)) > 20, obj


@pytest.mark.parametrize('name', list(PAIRS))
def test_a_pinned_pair_is_as_recorded(name):
    assert name in RECORD, f'{name!r}: not recorded yet -- run this module as a script'
    ref, twin = current()[name]
    ref0, twin0 = RECORD[name]
    moved = [s for s, a, b in (('reference', ref, ref0), ('twin', twin, twin0)) if a != b]
    if moved:
        p = PAIRS[name]
        other = {'reference': 'twin', 'twin': 'reference'}
        what = ' and '.join(f'the {s} ({", ".join(p[s])})' for s in moved)
        still = [other[s] for s in moved if other[s] not in moved]
        hint = (f'; {", ".join(n for s in still for n in p[s])} did not' if still else '')
        pytest.fail(f'pinned pair {name!r} moved on {what}{hint}.  Why they are a pair: '
                    f'{p["why"]}.  Check the twin, then re-record: '
                    f'python {__file__}')


def test_a_moved_side_is_named(monkeypatch):
    """The failure is the prompt: it names the side that moved, the one
    that did not, and the way to re-record."""
    name = next(iter(PAIRS))
    ref, _twin = current()[name]
    monkeypatch.setitem(RECORD, name, (ref, 'ffffffffffff'))
    with pytest.raises(pytest.fail.Exception) as e:
        test_a_pinned_pair_is_as_recorded(name)
    msg = str(e.value)
    assert 'moved on the twin' in msg and 'did not' in msg and 'Check the twin' in msg
    monkeypatch.setitem(RECORD, name, ('ffffffffffff', 'ffffffffffff'))
    with pytest.raises(pytest.fail.Exception) as e:
        test_a_pinned_pair_is_as_recorded(name)
    assert 'the reference' in str(e.value) and 'the twin' in str(e.value)


def test_the_generated_methods_are_defined_once():
    from pycircuit.circuit.hdl import BehaviouralMeta
    counts = {}
    for c in BehaviouralMeta.__init__.__code__.co_consts:
        if isinstance(c, types.CodeType):
            counts[c.co_name] = counts.get(c.co_name, 0) + 1
    assert {m: counts.get(m, 0) for m in GENERATED} == dict.fromkeys(GENERATED, 1)


def main():
    print('RECORD = {')
    for name, (ref, twin) in current().items():
        print(f'    {name!r}: ({ref!r}, {twin!r}),')
    print('}')


if __name__ == '__main__':
    sys.exit(main())
