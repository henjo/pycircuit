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
| `wright_cb.py` | the same with abscissae bounded to [0, 1] — produced the SHIPPED GLM4 tableau |
| `wright_all.py` | the full 24-permutation sweep with every constraint (3 h 32 min; found no feasible point — the recorded frontier) |
| `glm_irks.py`, `glm_irks_v.py` | the IRKS-conditioned searches; `_v` (nilpotency as a POLYNOMIAL rather than an eigenvalue condition) is what produced GLM3 |
| `stage_predictor.py` | the stage predictor's cost measurement — a driven RC with a STATE-FREE exponential, device `i` evaluations and the WORST stage's Newton iterations, with `glm_predictor='none'` as the control. ⚠ `Diode` cannot measure this at all: its `G` linearises around a stored `_vlim`, so an all-zeros seed gives the same iteration histogram as the exact one |

⚠ Two traps these scripts encode, both measured the hard way:
the Nordsieck convention is `y_k = h^k y^(k)` **without** `1/k!`, and `F` is the exponential of the
whole product `K(I+λK)⁻¹`, not `exp(K)·(I+λK)⁻¹`.
