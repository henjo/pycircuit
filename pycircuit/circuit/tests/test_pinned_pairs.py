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
    'the controller law and its C': {
        'why': 'the C makes `IntegralController._evaluate` past its error test with the '
               'default band: `0.9 ** p`, the shrink with `_retry_exponent` (the same point '
               'an index where the Python compares by `is`), the grow, the cap, `_judged` -- '
               'Python\'s float operations, `pow` and `log`, its `min` and `max` in order; '
               'an early abort\'s rejection (`reject_early`) is the same law at its estimate, '
               'which loop.c\'s `early_reject` hands it',
        'reference': ['pycircuit.circuit.stepcontroller:IntegralController._evaluate',
                      'pycircuit.circuit.stepcontroller:IntegralController._shrink',
                      'pycircuit.circuit.stepcontroller:IntegralController.reject_early',
                      'pycircuit.circuit.stepcontroller:StepController._retry_exponent',
                      'pycircuit.circuit.stepcontroller:StepController._band_decision',
                      'pycircuit.circuit.stepcontroller:StepController._band_target',
                      'pycircuit.circuit.stepcontroller:StepController._judged'],
        'twin': ['pycircuit.circuit._tran_ctrl_c:CTRL_C',
                 'pycircuit.circuit._tran_ctrl_c:new_state'],
    },
    'the stepping loop and its C': {
        'why': 'the C makes `_SteppingLoop`\'s multistep attempts, fixed and adaptive -- '
               'where (the grid or the controller\'s step, the next breakpoint as '
               '`_next_breakpoint` makes it, the run\'s end merging a corner within minbreak, crossed or landed on, the order drop a landing '
               'or a force asks for), take (the predictor, `check_order_drop`, '
               '`bdf2_alphas`, the C Newton set up as its glue sets it up, the screen with the '
               'run\'s certificate where the Python made it), judge '
               '(the C error test as its glue sets it up, the controller\'s law, a rejection '
               'retried within the budget) and accept (the rings, `_pred_note`, the '
               'statistics) -- and writes back what the Python attempts leave; the early '
               'abort armed as `early_arm` arms it, an abandoned attempt rejected as `take` '
               'rejects it (the law, the point\'s counts, nothing recorded as made)',
        'reference': ['pycircuit.circuit.transient:_SteppingLoop.where',
                      'pycircuit.circuit.transient:_SteppingLoop.take',
                      'pycircuit.circuit.transient:_SteppingLoop.early_arm',
                      'pycircuit.circuit.transient:_SteppingLoop.judge',
                      'pycircuit.circuit.transient:_SteppingLoop.accept',
                      'pycircuit.circuit.transient:_LMMSteps.judge',
                      'pycircuit.circuit.transient:Transient._solve_timestep',
                      'pycircuit.circuit._tran_events:_StepEvents._next_breakpoint',
                      'pycircuit.circuit._tran_history:_RunHistory._roll_history',
                      'pycircuit.circuit._tran_history:_RunHistory._push_history',
                      'pycircuit.circuit._tran_history:_RunHistory._lte_inputs',
                      'pycircuit.circuit._tran_predictor:_StagePredictor._pred_note',
                      'pycircuit.circuit._tran_predictor:_StagePredictor._pred_promote',
                      'pycircuit.circuit.integrator:Gear2Integrator.check_order_drop',
                      'pycircuit.circuit.integrator:TrapezoidalIntegrator.check_order_drop',
                      'pycircuit.circuit._lte_kernels:bdf2_alphas',
                      'pycircuit.circuit.stepcontroller:IntegralController._evaluate',
                      'pycircuit.circuit._tran_newton_c:solve',
                      'pycircuit.circuit._tran_lte_c:max_error'],
        'twin': ['pycircuit.circuit._tran_loop_c:LOOP_C',
                 'pycircuit.circuit._tran_loop_c:run',
                 'pycircuit.circuit._tran_loop_c:_write_back',
                 'pycircuit.circuit._tran_loop_c:_stats_back',
                 'pycircuit.circuit._tran_loop_c:_count',
                 'pycircuit.circuit._tran_loop_c:_why_not_the_run',
                 'pycircuit.circuit._tran_loop_c:_sources',
                 'pycircuit.circuit._tran_loop_c:_read',
                 'pycircuit.circuit._tran_loop_c:_bound',
                 'pycircuit.circuit._tran_loop_c:_load',
                 'pycircuit.circuit._tran_loop_c:_call',
                 'pycircuit.circuit._tran_loop_c:_bind',
                 'pycircuit.circuit._tran_loop_c:_judge'],
    },
    'the predictor and its C': {
        'why': 'the C makes `_predict_fast`: the nodes newest first with the 1e-13 dedupe '
               'and the tie test, the fit\'s times in Python\'s float order, np.vander\'s '
               'running products, numpy\'s own dgesv and dgemv, the clamp in clip\'s '
               'order; the general path\'s cases and a reported flag are handed back',
        'reference': ['pycircuit.circuit._tran_predictor:_StagePredictor._predict_fast',
                      'pycircuit.circuit._tran_predictor:_StagePredictor._pred_fast_nodes',
                      'pycircuit.circuit._tran_predictor:_StagePredictor._fit_fast'],
        'twin': ['pycircuit.circuit._tran_pred_c:PRED_C',
                 'pycircuit.circuit._tran_pred_c:predict'],
    },
    'the branch screen and its C': {
        'why': 'the C makes `_branch_screen_C`\'s cheap level -- |C|\'s largest entry '
               '(on the core\'s pattern where the core wrote C), Python\'s `max` of the '
               'running maximum, the diagonal\'s largest and smallest positive entry, the '
               'two tests -- and where it fires the certificate (`_certify`: the entries '
               'between B\'s blocks, each block\'s numpy inv as dgesv on the identity, its '
               'sums of squares in cumsum\'s order, the bound against tol (1 + delta)) on '
               'the run\'s data as `_branch_cert_data` makes them, at the same screen -- the '
               'same comparisons on the same values; the SVD stays the Python\'s '
               '(`_branch_svd`); a non-finite entry is handed back',
        'reference': ['pycircuit.circuit._tran_branch:_BranchCheck._branch_screen_C',
                      'pycircuit.circuit._tran_branch:_BranchCheck._branch_screen',
                      'pycircuit.circuit._tran_branch:_BranchCheck._branch_certify',
                      'pycircuit.circuit._tran_branch:_BranchCheck._branch_cert_data',
                      'pycircuit.circuit._tran_branch:_BranchCheck._branch_cert_kept',
                      'pycircuit.circuit._tran_branch:_BranchCheck._branch_structural_rank',
                      'pycircuit.circuit._tran_branch:_certify',
                      'pycircuit.circuit._tran_branch:_cert_data'],
        'twin': ['pycircuit.circuit._tran_branch_c:BRANCH_C',
                 'pycircuit.circuit._tran_branch_c:screen',
                 'pycircuit.circuit._tran_branch_c:cert_bind',
                 'pycircuit.circuit._tran_branch_c:_cert'],
    },
    'the sources and their C': {
        'why': 'the C makes each waveform\'s `f` in its operations and order (Python\'s '
               'float `%` and `//` for the pulse\'s period, the `where` chain\'s last true '
               'branch, `_sin_scalar` and its guards), the direct pass\'s fast add and '
               'its scatter in element order from +0.0, and each `next_event`; it '
               'declines where the Python takes another route',
        'reference': ['pycircuit.circuit.func:_sin_scalar',
                      'pycircuit.circuit.func:Sin.f', 'pycircuit.circuit.func:Sin.next_event',
                      'pycircuit.circuit.func:Pulse.f', 'pycircuit.circuit.func:Pulse.next_event',
                      'pycircuit.circuit.func:PWL.f', 'pycircuit.circuit.func:PWL.next_event',
                      'pycircuit.circuit.func:Exp.f', 'pycircuit.circuit.func:Exp.next_event',
                      'pycircuit.circuit.func:AM.f', 'pycircuit.circuit.func:SFFM.f',
                      'pycircuit.circuit._stamp_plan:_SourceDirect.run'],
        'twin': ['pycircuit.circuit._tran_sources:SOURCES_C',
                 'pycircuit.circuit._tran_sources:table_for',
                 'pycircuit.circuit._tran_sources:_waveform',
                 'pycircuit.circuit._tran_sources:Table'],
    },
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
    'the sources and their direct pass': {
        'why': 'a pass of VS- and IS-family sources is each source\'s parameter plus '
               'its time function at `t` -- `u`\'s own add, its warnings re-emitted '
               'from `u`\'s line -- scattered as the loop\'s bincount (every bin '
               'from +0.0, element order), in Python floats; its stamp watches what '
               'its full check read, the rest compared on every call',
        'reference': ['pycircuit.circuit.elements:VS.u',
                      'pycircuit.circuit.elements:IS.u',
                      'pycircuit.circuit.circuit:SubCircuit.u',
                      'pycircuit.circuit.circuit:SubCircuit._add_element_subvectors',
                      'pycircuit.circuit.circuit:SubCircuit._scatter_1d'],
        'twin': ['pycircuit.circuit._stamp_plan:source_direct',
                 'pycircuit.circuit._stamp_plan:_SourceDirect',
                 'pycircuit.circuit._stamp_plan:_added',
                 'pycircuit.circuit._stamp_plan:_u_add_site',
                 'pycircuit.circuit._stamp_plan:_genuine'],
    },
    'the sine and its scalar form': {
        'why': '`Sin.f` at a scalar time is its expression in Python floats around '
               'the toolkit\'s own `exp` and `sin`: the same IEEE operations in the '
               'same order, served only where nothing the expression does could warn',
        'reference': ['pycircuit.circuit.func:Sin.f'],
        'twin': ['pycircuit.circuit.func:_sin_scalar'],
    },
    'the reduce and its small take': {
        'why': 'a square matrix of up to `_TAKE_2D` rows or a vector of up to '
               '`_TAKE_1D` entries, `n` an int within it, is one `take` of the kept '
               'entries\' flat indices: the slices\' entries bit for bit, a fresh '
               'C-contiguous array of the input\'s dtype',
        'reference': ['pycircuit.circuit.analysis:_reduce_ndarray'],
        'twin': ['pycircuit.circuit.analysis:_reduce_small'],
    },
    'the insert and its small copy': {
        'why': 'a float64 vector of up to `_TAKE_1D` entries, `n` an int within it '
               'or at its end, is zeros and one indexed copy of the vector: the '
               'slices\' values bit for bit, the inserted entry +0.0',
        'reference': ['pycircuit.circuit.analysis:insert_row'],
        'twin': ['pycircuit.circuit.analysis:_insert_small'],
    },
    'the Newton solve and its C': {
        'why': 'the C solve is `nrsolver`\'s plain and chord Newton on the reduced '
               'system (numpy\'s and SciPy\'s own LAPACK, the walk, the test in its '
               'order, `dx` re-taken where a limiter exists), `jacobian_only`\'s '
               'converged-point evaluation, `_newton`\'s statistics and branch check; the '
               'early abort\'s check at the first iterate (`_early_estimate` on the seed\'s '
               'factors -- numpy\'s solve of J0 is dgetrs on dgesv\'s LU, the chord\'s '
               '`fac`) and its stop (`_early`: the counts `_newton` makes up to it); '
               'not the noncon test (`_limiting.NONCON`, 2026-10-15): only a stateful '
               'limiter counts, and the C declines every circuit with one '
               '(`newton_c:stateful`); the chord\'s settled test (`_chord_settled`, '
               'newton.c\'s `chord_settled` and `chord_dominant`, roadmap 136); the full '
               'Newton\'s sparse solve (`StructuralKLU.solve` on its refactor path -- '
               'the run\'s numeric refactored, the residual\'s validation by numpy\'s '
               'matmul -- newton.c\'s `sparse_full`, every other course handed back; '
               '`_record` gives it the run\'s numeric, `_serve` counts its refactors, '
               '2026-10-16)',
        'reference': ['pycircuit.circuit.nrsolver:StandardNewton.solve_system',
                      'pycircuit.circuit.nrsolver:ChordNewton.solve_system',
                      'pycircuit.circuit.nrsolver:_chord_settled',
                      'pycircuit.circuit.dcanalysis:refnode_removed',
                      'pycircuit.circuit.analysis:insert_row',
                      'pycircuit.circuit.analysis:remove_row_col',
                      'pycircuit.circuit._tran_newton:_StepNewton._newton',
                      'pycircuit.circuit._tran_newton:_StepNewton._newton_limiter',
                      'pycircuit.circuit._tran_newton:_StepNewton._early_estimate',
                      'pycircuit.circuit.transient:Transient._solve_timestep',
                      ## (the sparse chord: StructuralKLU's factor, the
                      ## owned numeric's solve -- stage 2, step d)
                      'pycircuit.circuit.linearsolver:StructuralKLU.factor',
                      'pycircuit.circuit.linearsolver:StructuralKLU._factor',
                      'pycircuit.circuit.linearsolver:_KLUNumeric._run',
                      ## (the sparse full Newton: StructuralKLU's refactor path)
                      'pycircuit.circuit.linearsolver:StructuralKLU.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU._ours'],
        'twin': ['pycircuit.circuit._tran_newton_c:NEWTON_C',
                 'pycircuit.circuit._tran_newton_c:solve',
                 'pycircuit.circuit._tran_newton_c:_klu_bind',
                 'pycircuit.circuit._tran_newton_c:_record',
                 'pycircuit.circuit._tran_newton_c:_serve',
                 'pycircuit.circuit._tran_newton_c:_Ctx',
                 'pycircuit.circuit._tran_newton_c:_walk_full',
                 'pycircuit.circuit._tran_newton_c:_arm',
                 'pycircuit.circuit._tran_newton_c:_arm_set',
                 'pycircuit.circuit._tran_newton_c:_early'],
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
               '`np.max`, the running reference written only on success; an unknown '
               '`lte_mask` leaves out is +0.0 before the maximum; on the rows '
               '`lte_charge` marks the absolute tolerance is `charge_floor`\'s; the early '
               'abort\'s estimate (`_early_estimate`) is the same chain at the first iterate, '
               'J^-1 by factors in hand (lte.c\'s `factored_solve`, its context '
               '`early_context`)',
        'reference': ['pycircuit.circuit.stepcontroller:StepController._charge_lte',
                      'pycircuit.circuit._tran_newton:_StepNewton._early_estimate',
                      'pycircuit.circuit.stepcontroller:StepController._normalised',
                      'pycircuit.circuit.stepcontroller:StepController._reference',
                      'pycircuit.circuit.stepcontroller:StepController.tolerance',
                      'pycircuit.circuit.stepcontroller:StepController._max_error',
                      'pycircuit.circuit.stepcontroller:normalised_error',
                      'pycircuit.circuit.stepcontroller:sigglobal_reference',
                      'pycircuit.circuit.stepcontroller:charge_floor',
                      'pycircuit.circuit.integrator:Gear2Integrator.compute_lte',
                      'pycircuit.circuit.integrator:TrapezoidalIntegrator.compute_lte',
                      'pycircuit.circuit.integrator:EulerIntegrator.compute_lte',
                      ## (the sparse solve: StructuralKLU's, stage 3 step b)
                      'pycircuit.circuit.linearsolver:AutoSolver.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU._ours',
                      'pycircuit.circuit._lte_kernels:third_divided_difference',
                      'pycircuit.circuit.analysis:remove_row_col'],
        'twin': ['pycircuit.circuit._tran_lte_c:LTE_C',
                 'pycircuit.circuit._tran_lte_c:max_error',
                 'pycircuit.circuit._tran_lte_c:early_context',
                 'pycircuit.circuit._tran_lte_c:_load',
                 'pycircuit.circuit._tran_lte_c:_globals',
                 'pycircuit.circuit._tran_lte_c:_Ctx',
                 'pycircuit.circuit._tran_lte_c:_klu_bind'],
    },
    'the DC base solve and its C': {
        'why': 'the C base solve is `StandardNewton.solve_system` on the DC operating '
               'point\'s compositions: `dcanalysis`\'s residual (`func`, '
               '`source_callback`, `refnode_removed`, the reduced G), gshunt (`F + x '
               'g`, `_plus_diagonal`), junction gmin (`_junction_gmin`), source '
               'stepping\'s factor; the solve (numpy\'s dgesv, or StructuralKLU\'s '
               'refactor, fresh factor and residual check), the DC\'s limiter, the '
               'convergence test and the iteration limit\'s exception; not the noncon '
               'test (`_limiting.NONCON`, 2026-10-15): only a stateful limiter counts, '
               'and the C declines every circuit with one (`dcnewton_c:unservable`)',
        'reference': ['pycircuit.circuit.nrsolver:StandardNewton.solve_system',
                      'pycircuit.circuit.nrsolver:_standard_failure',
                      'pycircuit.circuit.nrsolver:_plus_diagonal',
                      'pycircuit.circuit.nrsolver:_junction_gmin',
                      'pycircuit.circuit.nrsolver:GminSteppingNewton.solve_system',
                      'pycircuit.circuit.nrsolver:JunctionGminSteppingNewton.solve_system',
                      'pycircuit.circuit.nrsolver:SourceSteppingNewton.solve_system',
                      'pycircuit.circuit.nrsolver:SourceSteppingNewton._adaptive',
                      'pycircuit.circuit.dcanalysis:DC._solve_dc',
                      'pycircuit.circuit.dcanalysis:DC._newton',
                      'pycircuit.circuit.dcanalysis:refnode_removed',
                      'pycircuit.circuit.dcanalysis:_ReducedG',
                      'pycircuit.circuit.linearsolver:AutoSolver.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU._ours'],
        'twin': ['pycircuit.circuit._dc_newton_c:DCNEWTON_C',
                 'pycircuit.circuit._dc_newton_c:solve',
                 'pycircuit.circuit._dc_newton_c:composition',
                 'pycircuit.circuit._dc_newton_c:_load',
                 'pycircuit.circuit._dc_newton_c:_record',
                 'pycircuit.circuit._dc_newton_c:_solver',
                 'pycircuit.circuit._dc_newton_c:_walk',
                 'pycircuit.circuit._dc_newton_c:_klu_bind',
                 'pycircuit.circuit._dc_newton_c:_write_back',
                 'pycircuit.circuit._dc_newton_c:_Rec'],
    },
    'the PCNR DC solve and its C': {
        'why': 'the C solve is `pcnr.solve_dc`\'s loop from its start: the ordinary '
               'part on the circuit\'s view (`pcnr_view`, the core on it), each '
               'participant at its own unknowns (`_kernel_map`, `_kernel_stamp`: the '
               'walk, the mirrors, the remainder, the block), `augmented_system`\'s '
               'order, gshunt\'s rows and the junction-gmin term\'s (`jgmin` across '
               '`pcnr_junction_pairs`, in pair order), `schur_reduce`\'s rank-one '
               'updates, the reduction (the reference and each held row out: '
               '`predict`\'s keep, the C\'s row map `red`), the solve '
               '(numpy\'s dgesv, or StructuralKLU\'s sequence after the AutoSolver\'s '
               'own first choice), `dx_lim_of`, `refine` -- each class\'s law '
               'parameters, its laws (`apply_limit`) and `limit_block`\'s Kruskal and '
               'potentials -- the three tests, and the iteration limit\'s, a '
               'non-finite iterate\'s and a singular solve\'s exceptions',
        'reference': ['pycircuit.circuit.pcnr:solve_dc',
                      'pycircuit.circuit.pcnr:augmented_system',
                      'pycircuit.circuit.pcnr:schur_reduce',
                      'pycircuit.circuit.pcnr:predict',
                      'pycircuit.circuit.pcnr:dx_lim_of',
                      'pycircuit.circuit.pcnr:refine',
                      'pycircuit.circuit.pcnr:limit_block',
                      'pycircuit.circuit.pcnr:lim_converged',
                      'pycircuit.circuit.pcnr:pcnr_view',
                      'pycircuit.circuit.pcnr:_remainder_of',
                      'pycircuit.circuit.pcnr:PcnrDevice.stamp',
                      'pycircuit.circuit.pcnr:PcnrDevice._kernel_map',
                      'pycircuit.circuit.pcnr:PcnrDevice._kernel_stamp',
                      'pycircuit.circuit.pcnr:PcnrDevice.limit',
                      'pycircuit.circuit._limiting:apply_limit',
                      'pycircuit.circuit.linearsolver:AutoSolver.solve',
                      'pycircuit.circuit.linearsolver:AutoSolver._select',
                      'pycircuit.circuit.linearsolver:StructuralKLU.solve',
                      'pycircuit.circuit.linearsolver:StructuralKLU._ours',
                      'pycircuit.circuit.linearsolver:StructuralKLU._factor'],
        'twin': ['pycircuit.circuit._pcnr_c:PCNRDC_C',
                 'pycircuit.circuit._pcnr_c:REFINE_C',
                 'pycircuit.circuit._pcnr_c:render_refine', 'pycircuit.circuit._pcnr_c:_chains',
                 'pycircuit.circuit._pcnr_c:_block_c',
                 'pycircuit.circuit._pcnr_c:refine',
                 'pycircuit.circuit._pcnr_c:solve',
                 'pycircuit.circuit._pcnr_c:_outcome',
                 'pycircuit.circuit._pcnr_c:_settings',
                 'pycircuit.circuit._pcnr_c:_solver',
                 'pycircuit.circuit._pcnr_c:_choose',
                 'pycircuit.circuit._pcnr_c:_record',
                 'pycircuit.circuit._pcnr_c:_class_ok',
                 'pycircuit.circuit._pcnr_c:_tables',
                 'pycircuit.circuit._pcnr_c:_Rec',
                 'pycircuit.circuit._hdl_climit:_parameter_chains'],
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
               'every piece checked as defined at the run\'s start; `take`\'s early abort '
               'is never armed on a fixed grid (`early_arm`)',
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
    "the step's end passes and the C's": {
        'why': 'once converged the C makes the calls `_stage_end_passes` makes next -- '
               '`passes`\'s core call, its bits, outputs and dummies, at the last stage '
               'with all four passes and at the first two with `C` and `G` where the '
               'shooting reads them -- and `_stage_end_passes` takes them where each call '
               'would find its readiness as stamped (`passes_held`: `core_for`\'s record, '
               'the probe), making the counters each call makes (`passes_count`)',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._stage_end_passes',
                      'pycircuit.circuit._tran_core:passes',
                      'pycircuit.circuit._tran_core:core_for',
                      'pycircuit.circuit._tran_core:_Core.ready',
                      'pycircuit.circuit._tran_core:_Core.probe'],
        'twin': ['pycircuit.circuit._tran_radau_tc:RADAU_TC',
                 'pycircuit.circuit._tran_radau_tc:solve',
                 'pycircuit.circuit._tran_core:passes_held',
                 'pycircuit.circuit._tran_core:passes_count'],
    },
    'the frozen factors and their C': {
        'why': 'the C makes `_radau_frozen`\'s two factors as numpy does -- the real '
               'one in two roundings, the complex one with `Cr` cast to complex and '
               'numpy\'s product in `cmul_mode`\'s form, then `Gr` cast to complex -- '
               'with the four flags its `errstate(all=\'raise\')` raises on, a '
               'non-finite input left to numpy, and the complex one\'s nonzeros where '
               'they fall in the last record\'s pattern, in `_csc_of_dense`\'s order; '
               '`prepare_values` makes `prepare`\'s record of them',
        'reference': ['pycircuit.circuit._tran_radau:_RadauStages._radau_frozen',
                      'pycircuit.circuit.linearsolver:_csc_of_dense',
                      'pycircuit.circuit.linearsolver:ComplexKLUSolver.prepare',
                      'pycircuit.circuit._tran_radau_tc:cmul_mode'],
        'twin': ['pycircuit.circuit._tran_radau_tc:FOLD_C',
                 'pycircuit.circuit._tran_radau_tc:fold',
                 'pycircuit.circuit._tran_radau_tc:_fold_bufs',
                 'pycircuit.circuit._tran_radau_tc:_fold_nz',
                 'pycircuit.circuit._tran_radau:_lam_folds',
                 'pycircuit.circuit.linearsolver:_csc_of_values',
                 'pycircuit.circuit.linearsolver:ComplexKLUSolver.prepare_values'],
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
    'the stage step map and its C': {
        'why': 'the plain dense walk\'s step in one call: the block column-major (as '
               '`lu_factor` copies it) from the stages\' full `C` and `G`, the reference '
               'row and column skipped as the readers reduce them, in `_stage_block`\'s '
               'operations (`h A_ij` first, `C_i` added on the diagonal), SciPy\'s own '
               '`getrf` and `getrs` on `C_n P` stacked per stage, the new `P` the last '
               'stage block in `getrs`\' Fortran order; it hands back wherever the Python '
               'step declines, warns or raises, and counts the readers\' sync skips',
        'reference': ['pycircuit.circuit.shooting._pss_walks:_PeriodWalks._stage_step',
                      'pycircuit.circuit.shooting._pss_walks:_stage_block',
                      'pycircuit.circuit.shooting._pss_walks:_stage_reads',
                      'pycircuit.circuit.shooting._steps:_StageStep.solve',
                      'pycircuit.circuit.shooting._numerics:_lu_solve_split',
                      'pycircuit.circuit.shooting._numerics:_complex_solve',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._C_at',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._G_at',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._factorise',
                      'pycircuit.circuit.linearsolver:DenseSolver.factor',
                      'pycircuit.circuit.linearsolver:AutoSolver.factor',
                      'pycircuit.circuit.linearsolver:lu_solve'],
        'twin': ['pycircuit.circuit.shooting._sens_c:SENS_C_SRC',
                 'pycircuit.circuit.shooting._sens_c:step',
                 'pycircuit.circuit.shooting._sens_c:_map',
                 'pycircuit.circuit.shooting._sens_c:_Ctx'],
    },
    'the stage readers and their one read': {
        'why': 'where the readers would read just that -- no junction, nothing keeping '
               'limiting state (the sync skips), each full stage\'s reference entry +0.0 '
               '(the state `_insert_refnode` rebuilds), the memo holding both matrices -- '
               'each stage\'s `C` and `G` are the step\'s memo record at the full stage, '
               'reduced by the readers\' own `remove_row_col`, the sync\'s skip counted '
               'before each `G` as `_G_at` counts it',
        'reference': ['pycircuit.circuit.shooting._pss_inner:_InnerTransient._C_at',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._G_at',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._sync_limit_at',
                      'pycircuit.circuit.shooting._pss_inner:_InnerTransient._insert_refnode',
                      'pycircuit.circuit._tran_companion:_CompanionModel._memo_get'],
        'twin': ['pycircuit.circuit.shooting._pss_walks:_stage_reads',
                 'pycircuit.circuit.shooting._pss_walks:_readers'],
    },
    'numpy\'s solve and its kept LU': {
        'why': 'the kept LU is numpy\'s own call -- its OpenBLAS\'s `dgesv` on a '
               'Fortran copy, one right-hand side -- then `dgetrs` on its factors; '
               'not `dgetrf` first (threaded from 100 unknowns: other bits), and '
               'numpy\'s error at every solve of a singular matrix; a refilled LU '
               'copies its matrix as a fresh one does',
        'reference': ['pycircuit.circuit._numeric:linearsolver'],
        'twin': ['pycircuit.circuit.linearsolver:_numpy_lapack',
                 'pycircuit.circuit.linearsolver:NumpyLU'],
    },
    'the complex factor\'s CSC, SciPy\'s and numpy\'s': {
        'why': 'numpy builds the arrays `csc_matrix(A).astype(complex128)` holds -- '
               'the nonzeros (a NaN is one, a signed zero is not) column by column, '
               'rows ascending, copied -- and the product as SciPy\'s '
               '`_matmul_vector` runs it, a zero vector and `csc_matvec`; the last '
               'record\'s `Ap`, `Ai` and key kept only where the nonzeros fall '
               'where they fell',
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
    'the behavioural source and its HDL form': {
        'why': '`BSourceHdl` makes `BSource`\'s stamps from the same callables -- the same '
               'terminals and orientation, the current and charge as the callables trace, '
               'the derivative exact where `BSource` takes the toolkit\'s central '
               'difference -- and its tests hold it to `BSource`\'s values',
        'reference': ['pycircuit.circuit.elements:BSource',
                      'pycircuit.circuit.toolkit:Toolkit.derivative'],
        'twin': ['pycircuit.circuit.elements_hdl:BSourceHdl',
                 'pycircuit.circuit.elements_hdl:_bsource_class',
                 'pycircuit.circuit.elements_hdl:_bsource_trace'],
    },
}

#: `pair name: (reference digest, twin digest)` -- re-made by running this
#: module as a script, after the twin was checked.
RECORD = {
    'the controller law and its C': ('1524666662a6', '716b6f29ec13'),
    'the stepping loop and its C': ('c7c3ea89fda0', 'ed148322c3cc'),
    'the predictor and its C': ('f2c284eb392d', '2f1a2f96e90f'),
    'the branch screen and its C': ('adb832d7aae0', '1da8be0e9b2e'),
    'the sources and their C': ('55cf2d3441f8', '41b0a7d59bbe'),
    'the kernel call and the pass driver': ('09a9e274bbe8', 'f44eb4ff41df'),
    'the limiting loop and the walk': ('58c4ce9c8f4d', '4ed13cca30aa'),
    'the limiter laws and their C prelude': ('9b7944f9c0dd', '7e8f23b73019'),
    'the transient evaluation and the core': ('ecdae2231dae', 'c1699314005b'),
    'the assembly loops and the plan': ('fa280abb9bf9', '14651de8716f'),
    'the sources and their direct pass': ('61c7a4bed013', 'fced12f39974'),
    'the sine and its scalar form': ('5096cada457f', '827ecf798b4d'),
    'the reduce and its small take': ('f4d40012ff79', '41b65c2ead50'),
    'the insert and its small copy': ('8b4df58556f0', '4827f496e6dd'),
    'the Newton solve and its C': ('07b811013bb8', 'b393c22f15df'),
    'the PSP limiter and its C twin': ('8e14bbc903f1', '58969d2d34a6'),
    'the error test and its C': ('1736570f3f0e', 'a1bf1ec171ad'),
    'the DC base solve and its C': ('0db2906ca846', '60e5fa192106'),
    'the PCNR DC solve and its C': ('818043efa34e', '64a5819e64c9'),
    'the stage predictor and its multistep fast path': ('739b71ddbace', 'f2c284eb392d'),
    'the stepping loop and its unread Jacobian': ('de8b4b288065', 'a0ac4588446c'),
    'the printed pass kernels and their fused kernel': ('e61140203114', '1cc7e58f8006'),
    "the stage paths' passes and the core's": ('f83268db1289', '131833be2877'),
    'the coupled stage Newton and its C': ('c9e9adf33989', '2895432d069c'),
    'the transform solve and its frozen form': ('fae3de6ce9eb', '9b7970f8fbb7'),
    'the transform Newton and its C': ('d7c24c3c7e7d', 'deb779100728'),
    "the step's end passes and the C's": ('5501fef13bd5', 'cc482038c73e'),
    'the frozen factors and their C': ('c27dee545b25', '00491721f3f9'),
    'the stage block and its one pass': ('a013cd007e0c', '24be0c26211b'),
    'the stage step map and its C': ('7507f30134bb', 'a2de0dd5aa86'),
    'the stage readers and their one read': ('d25552fe5ede', '08f3606179d1'),
    "numpy's solve and its kept LU": ('1f8908b895c7', '119066a97f50'),
    "the complex factor's CSC, SciPy's and numpy's": ('2e7218131272', '033c4ac51e63'),
    'the readiness reads and their stamps': ('61efd629dd90', '76e34cda0583'),
    'the behavioural source and its HDL form': ('eb811b119033', '6428eaf454a9'),
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
