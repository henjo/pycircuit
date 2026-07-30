Goal: Implement JAX batched evaluation via `vmap` and pure functions.

1. In `elements.py` and `semiconductors.py`, refactor `Diode`, `NonLinearVCCS`, `BSource`, `VSwitch`, `ISwitch`, `CCCS`, `BJT`, `JFET` to implement `eval_i_pure(x, params, epar, toolkit)` and `eval_q_pure(x, params, epar, toolkit)`.
   - Update `def i(self, x, epar=None)` to call `self.eval_i_pure(x, self.iparv_dict, epar, self.toolkit)`.
2. In `circuit.py`, rewrite `_add_element_submatrices` and `_add_element_subvectors` to use `self._eval_groups`.
   - If `hasattr(self, '_eval_groups') and self._eval_groups` and JAX is used:
     Iterate through `_eval_groups`.
     Slice `X_batch = x[group['nodemaps_2d']]` (or `1d` depending on vector/matrix).
     Evaluate using `batched_eval = self.toolkit.generate_batched_eval(cls, 'i' if methodname in ['i', 'G'] else 'q')`
     Add `jac_batch.flatten()` (or `val_batch.flatten()`) to the sparse arrays / lhs.
     Remove these batched instances from the sequential python loop.
3. Verify with `pytest pycircuit/circuit/tests/test_analysis_transient_stress.py`.
