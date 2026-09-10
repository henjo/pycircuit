# GLM construction scripts — the provenance of the shipped tableaux

These produced the coefficients in `pycircuit/circuit/integrator.py`'s `GLM2Integrator`,
`GLM3Integrator` and `GLM4Integrator`, and the negative results recorded in
`doc/pss_roadmap_260902.md`. They were written in a session scratchpad and are kept here because
the roadmap cites them: a citation that points at a temporary directory is not a record.

| script | what it is |
|---|---|
| `glm_tools.py` | Nordsieck order conditions (`U = C − ACK`, `V = E − BCK`), polynomial-exactness check, `M(z)`, `M_∞`, and a free numerical search |
| `wright_ctor.py` | Wright's §3.7 closed form: `δ`, `N_k`, `F = exp(K(I+λK)⁻¹)`, `Ψ`, and `B̃` from (3.7.8) |
| `wright_full.py` | the whole method — `B̃` back-transformed by (3.5.10) `B = ΨB̃`, `Ã = B̃⁻¹JB̃` from (3.5.9), `A = Ã + λI`. **Run it: it reproduces Wright's printed p = 2 tableau to 5.8e-16, which is what validates the route.** |
| `wright_sa.py` | the stiffly accurate sub-class (§3.9.2) — produced GLM4's first tableau |
| `wright_p4.py` | the general p = 3 / p = 4 driver: least-squares on `A~ = B~^-1 J B~` coming out strictly lower triangular, per lambda and per permutation. `python wright_p4.py 4` |
| `wright_cb.py` | the same with abscissae bounded to [0, 1] — produced the SHIPPED GLM4 tableau |
| `wright_all.py` | the full 24-permutation sweep with every constraint (3 h 32 min; found no feasible point — the recorded frontier) |
| `glm_irks.py`, `glm_irks_v.py` | the IRKS-conditioned searches; `_v` (nilpotency as a POLYNOMIAL rather than an eigenvalue condition) is what produced GLM3 |

The stage predictor's own measurement moved out of this folder when it stopped being a GLM
thing: `benchmarks/stage_predictor.py` covers every integrator family.

NOT kept: two earlier drafts of the stiffly-accurate search (`wright_sc`, `wright_sc1`) produced
ZERO output — the same `c`-with-`p-1`-free-entries defect that `wright_sa` had, where every
evaluation raised. `wright_sa.py` is the one that works; there is nothing in the drafts to keep.

⚠ Two traps these scripts encode, both measured the hard way:
the Nordsieck convention is `y_k = h^k y^(k)` **without** `1/k!`, and `F` is the exponential of the
whole product `K(I+λK)⁻¹`, not `exp(K)·(I+λK)⁻¹`.

## Adaptive stepping (2026-09-10)

`glm_adaptive.py` measures the cost and the restart count of the GLM under the step controller, and the
order on a deliberately non-uniform grid — which is the measurement Voigtmann Thm 9.5 does NOT cover, since
it is stated at constant stepsize.
