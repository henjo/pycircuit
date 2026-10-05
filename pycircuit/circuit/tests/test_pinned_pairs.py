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
                 'pycircuit.circuit._tran_core:evaluate',
                 'pycircuit.circuit._tran_core:passes'],
    },
    'the assembly loops and the plan': {
        'why': 'the plan is the loop\'s bincount over the same values in the same '
            'order, with the constant elements pre-filled, the C-bound classes '
            'batched and the zero sources skipped -- the source pass the same '
            'elements called in the same order; its fallbacks are the loop',
        'reference': ['pycircuit.circuit.circuit:SubCircuit._add_element_submatrices',
                   'pycircuit.circuit.circuit:SubCircuit._add_element_subvectors'],
        'twin': ['pycircuit.circuit._stamp_plan:assemble_matrix',
              'pycircuit.circuit._stamp_plan:assemble_vector',
              'pycircuit.circuit._stamp_plan:assemble_source',
              'pycircuit.circuit._stamp_plan:_SourcePlan',
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
    'the stepping loop and its unread Jacobian': {
        'why': 'a fixed multistep run reads no step\'s Jacobian -- `judge` returns '
               'first, `event` does nothing, the other pieces carry it -- so '
               '`j_unread` marks its attempts (the converged point without G), '
               'every piece checked as defined at the run\'s start',
        'reference': ['pycircuit.circuit.transient:_SteppingLoop.where',
                      'pycircuit.circuit.transient:_SteppingLoop.take',
                      'pycircuit.circuit.transient:_SteppingLoop.judge',
                      'pycircuit.circuit.transient:_SteppingLoop.event',
                      'pycircuit.circuit.transient:_SteppingLoop.accept',
                      'pycircuit.circuit.transient:_StepFamily.attempt',
                      'pycircuit.circuit.transient:Transient._attempt_step'],
        'twin': ['pycircuit.circuit.transient:_SteppingLoop.__init__',
                 'pycircuit.circuit.transient:_SteppingLoop.execute'],
    },
    'the printed pass kernels and their fused kernel': {
        'why': 'the fused kernel is the passes\' printed statements unioned by name, '
               'each run under the bits of the passes that print it, the outputs per '
               'pass -- it parses exactly the lines `_render_c` prints and refuses '
               'any other; the core calls it once a call and scatters from its '
               'staging in the passes\' order',
        'reference': ['pycircuit.circuit.hdl:_render_c'],
        'twin': ['pycircuit.circuit._hdl_cbackend:fuse_csrc',
                 'pycircuit.circuit._hdl_cbackend:fuse_source',
                 'pycircuit.circuit._hdl_cbackend:bind_fused',
                 'pycircuit.circuit._tran_core:_fuse_groups'],
    },
    'the stage paths\' passes and the core\'s': {
        'why': 'where the core serves, a stage path takes the passes from one core '
               'call (formula 4) for the circuit\'s own `q`, `i`, `G`, `C` calls, '
               'which stay as its fallback -- the same values, fresh arrays',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._coupled_stage_context',
                      'pycircuit.circuit._tran_radau:_RadauStages._coupled_stage_solver',
                      'pycircuit.circuit._tran_radau:_RadauStages._radau_error_estimate',
                      'pycircuit.circuit._tran_radau:_RadauStages._rk_step_transformed',
                      'pycircuit.circuit._tran_radau:_RadauStages._transform_loop',
                      'pycircuit.circuit._tran_radau:_RadauStages._stage_end_passes',
                      'pycircuit.circuit._tran_stages:_SequentialStages._solve_implicit_stage',
                      'pycircuit.circuit._tran_stages:_SequentialStages._finish_stage_step'],
        'twin': ['pycircuit.circuit._tran_core:passes',
                 'pycircuit.circuit._tran_core:CORE_C'],
    },
    'the coupled stage Newton and its C': {
        'why': 'the C runs `_stage_newton` undamped and unshunted: the assembly\'s '
               'operations and order over the full width, numpy\'s solve, the '
               'reference row\'s `+ 0.0`, the walk against the previous stage, the '
               'convergence test after the assembly, every assembly\'s memo record '
               'and the source memo\'s calls and hits',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._coupled_stage_solver',
                      'pycircuit.circuit._tran_radau:_RadauStages._coupled_stage_system',
                      'pycircuit.circuit._tran_radau:_RadauStages._stages_converged',
                      'pycircuit.circuit._tran_companion:_CompanionModel._memo_put',
                      'pycircuit.circuit._tran_companion:_CompanionModel._source_at'],
        'twin': ['pycircuit.circuit._tran_radau_c:RADAU_C',
                 'pycircuit.circuit._tran_radau_c:solve'],
    },
    'the transform solve and its frozen form': {
        'why': 'the frozen solve makes once a step what the per-iteration solve '
               'makes every iteration -- `P`, both factors, the real LU, the '
               'complex factor\'s marshalling and refactor -- and calls the same '
               'right-hand side and update helpers; it keeps numpy\'s LU only where '
               'the analysis solver\'s choice is the dense one, each piece as defined',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._radau_transform_solve',
                      'pycircuit.circuit._tran_radau:_RadauStages._radau_complex_solve',
                      'pycircuit.circuit.analysis:Analysis._get_linearsolver',
                      'pycircuit.circuit.linearsolver:AutoSolver',
                      'pycircuit.circuit.linearsolver:DenseSolver.solve',
                      'pycircuit.circuit.linearsolver:ComplexKLUSolver.solve'],
        'twin': ['pycircuit.circuit._tran_radau:_FrozenTransform',
                 'pycircuit.circuit._tran_radau:_frozen_chain',
                 'pycircuit.circuit._tran_radau:_RadauStages._radau_frozen',
                 'pycircuit.circuit.linearsolver:ComplexKLUSolver.solve_prepared'],
    },
    'the transform Newton and its C': {
        'why': 'the C runs `_transform_loop` with the step\'s frozen factors: the stages\' '
               '`q` and `i` through the core, `_coupled_stage_system`\'s residual, '
               '`_transform_rhs` and `_transform_back` with numpy\'s complex product as '
               '`cmul_mode` reads it, the kept LU\'s `dgesv`-then-`dgetrs`, '
               '`solve_prepared`\'s refactor once a record and its residual check '
               '(decided with a margin), the walk, the test after the update',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._transform_loop',
                      'pycircuit.circuit._tran_radau:_RadauStages._coupled_stage_system',
                      'pycircuit.circuit._tran_radau:_RadauStages._stages_converged',
                      'pycircuit.circuit._tran_radau:_FrozenTransform.solve',
                      'pycircuit.circuit._tran_radau:_transform_rhs',
                      'pycircuit.circuit._tran_radau:_transform_back',
                      'pycircuit.circuit.linearsolver:NumpyLU.solve',
                      'pycircuit.circuit.linearsolver:ComplexKLUSolver.solve_prepared'],
        'twin': ['pycircuit.circuit._tran_radau_tc:RADAU_TC',
                 'pycircuit.circuit._tran_radau_tc:solve',
                 'pycircuit.circuit._tran_radau_tc:cmul_mode'],
    },
    'the stage block and its one pass': {
        'why': 'the one pass makes the block loop\'s elementwise products `(h A_ij) G_j` '
               'and diagonal sums `C_i + ...` in its operand order, placed by one '
               'transpose and reshape; it declines non-finite blocks (a NaN sum\'s '
               'payload follows numpy\'s loop, not the operand order) and anything '
               'that would raise, leaving them to the loop',
        'reference': ['pycircuit.circuit.shooting._pss_walks:_PeriodWalks._stage_step'],
        'twin': ['pycircuit.circuit.shooting._pss_walks:_stage_block'],
    },
    'numpy\'s solve and its kept LU': {
        'why': 'the kept LU is numpy\'s own call -- its OpenBLAS\'s `dgesv` on a '
               'Fortran copy, one right-hand side -- then `dgetrs` on its factors; '
               'not `dgetrf` first (threaded from 100 unknowns: other bits), and '
               'numpy\'s error at every solve of a singular matrix',
        'reference': ['pycircuit.circuit._numeric:linearsolver'],
        'twin': ['pycircuit.circuit.linearsolver:_numpy_lapack',
                 'pycircuit.circuit.linearsolver:NumpyLU'],
    },
    'the complex factor\'s CSC, SciPy\'s and numpy\'s': {
        'why': 'numpy builds the arrays `csc_matrix(A).astype(complex128)` holds -- '
               'the nonzeros (a NaN is one, a signed zero is not) column by column, '
               'rows ascending, copied -- and the product as SciPy\'s '
               '`_matmul_vector` runs it, a zero vector and `csc_matvec`',
        'reference': ['pycircuit.circuit.linearsolver:ComplexKLUSolver.prepare'],
        'twin': ['pycircuit.circuit.linearsolver:_csc_of_dense',
                 'pycircuit.circuit.linearsolver:_csc_dot',
                 'pycircuit.circuit.linearsolver:_csc_matvec'],
    },
    'the readiness reads and their stamps': {
        'why': 'a stamp stands for a full check while no dict it read has changed: '
               'every read of the plan lookup, the generated-pass test, the pack, a '
               'hand-written limiter\'s twin, the linear solver and the tolerances '
               'must be a watched dict or compared on every call',
        'reference': ['pycircuit.circuit._stamp_plan:_plan_for',
                      'pycircuit.circuit._hdl_batch:is_generated',
                      'pycircuit.circuit._hdl_cbackend:CKernel.pack',
                      'pycircuit.circuit._hdl_climit:_Handwritten',
                      'pycircuit.circuit.analysis:Analysis._get_linearsolver',
                      'pycircuit.circuit._tran_newton:_StepNewton._newton_tolerances'],
        'twin': ['pycircuit.circuit._watch:WATCH_C',
                 'pycircuit.circuit._watch:arm',
                 'pycircuit.circuit._tran_core:_Core._stamp',
                 'pycircuit.circuit._tran_core:_plan_stamp',
                 'pycircuit.circuit._tran_newton_c:_par_stamped',
                 'pycircuit.circuit._tran_newton_c:_par_stamp'],
    },
}

#: `pair name: (reference digest, twin digest)` -- re-made by running this
#: module as a script, after the twin was checked.
RECORD = {
    'the kernel call and the pass driver': ('09a9e274bbe8', 'c900976facbc'),
    'the limiting loop and the walk': ('d151206fce18', '411e5dd4e655'),
    'the limiter laws and their C prelude': ('9b7944f9c0dd', '770bc2e815cc'),
    'the transient evaluation and the core': ('ce82a66a3520', '9a913d20faa9'),
    'the assembly loops and the plan': ('4c77e2ad0e1e', '951ecc14d234'),
    'the stage predictor and its multistep fast path': ('739b71ddbace', 'f2c284eb392d'),
    'the Newton solve and its C': ('a4c57d1d7bee', '96a7d2fa7fe9'),
    'the error test and its C': ('a24ba5659ca4', 'd6903a790f88'),
    'the PSP limiter and its C twin': ('8e14bbc903f1', '58969d2d34a6'),
    'the stepping loop and its unread Jacobian': ('16a27211f0a2', '9be60e00fe1a'),
    'the printed pass kernels and their fused kernel': ('e61140203114', 'df928813e8cd'),
    'the readiness reads and their stamps': ('48b75bd9cace', 'e52c723f0eb8'),
    "the stage paths' passes and the core's": ('8e4b7cf1fa24', 'a3ae0f6c8d98'),
    'the coupled stage Newton and its C': ('5f7c56901f9a', 'e9156eef9c41'),
    'the transform solve and its frozen form': ('945834192f2a', 'f7146ba0f56e'),
    'the transform Newton and its C': ('d7c24c3c7e7d', 'b4413fe8ec29'),
    'the stage block and its one pass': ('a177e4f1decf', '24be0c26211b'),
    "numpy's solve and its kept LU": ('1f8908b895c7', 'f0bcd0f84207'),
    "the complex factor's CSC, SciPy's and numpy's": ('9de9a726e6a5', '9ce76450fa2a'),
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
