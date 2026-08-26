# Native PCNR in pycircuit — what it would take, and what it would touch

**Status:** design study. Nothing here is implemented. PCNR ships today as a
*layer* (`pycircuit/circuit/pcnr.py`, `pcnr=True`, off by default); this document
describes what it would mean to support it *natively* instead, what that costs,
and what else in pycircuit moves when it happens.

This is the reasoning document, not a plan. It is written to be argued with.

**Provenance.** Line references are as of `4a4ca33` and will rot; the surrounding
names are the durable part. Every measured number below names the script that
produces it, so it can be re-measured rather than trusted:
`benchmarks/transient_review/companion_consistency.py` (§2),
`benchmarks/transient_review/stage13_pcnr.py` (§4 bill 1), and the gate 13-6
instrumentation recorded in `doc/transient_work_plan.md` (§4 bill 2).

---

## 1. What PCNR is

Aadithya, Keiter & Mei (Sandia), *"Predictor/Corrector Newton-Raphson"*. The
problem it addresses is **Newton limiting** for exponential devices.

A diode's `exp(v/VT)` overflows or oscillates if Newton takes an undamped step,
so simulators limit the step. SPICE's `pnjlim` does this by choosing a
**linearisation point** `VD` near the new iterate but not beyond a safe
excursion, and the device then presents a linear companion model about `VD`.

PCNR's observation is that `VD` is a *hidden, shared* quantity. Two devices on
the same branch each adjust the branch voltage, and the second undoes the first;
the limiter's output depends on the order devices are visited. PCNR's fix is to
make each limited quantity an **explicit unknown of the system**:

```
x = [x_MNA ; x_lim]
```

with one residual row per limited quantity, `v_Dk − (e_a − e_b) = 0`, and each
device evaluating at **its own** `v_Dk`. A Newton iteration then splits into

- **predict** — a full Newton step on the enlarged system, and
- **correct** — each device limiting *only the variables it owns*.

Because no two devices share an unknown, the clash cannot arise by construction
rather than by ordering. The enlarged system is not more expensive to solve: the
`lim/lim` block is the identity, so the border collapses by Schur complement to a
rank-`k` update on a matrix of the original MNA size.

---

## 2. A correction that changes what "native" means

**The existing documentation states that pycircuit's limiting makes the Jacobian
and the residual inconsistent. It does not, and this was verified rather than
reasoned about.**

`doc/src/circuit/pcnr_limiting.rst` says `Diode.G` linearises around the limited
voltage while `Diode.i` evaluates at the node voltage, "so the Jacobian and the
residual are taken at different points". Reading `elements.py`, `Diode.i` ends:

```python
I_eff = I - g * (VD - v_nodes)
```

which is the **companion form**: `i(x) = I(VD) + g(VD)·(v(x) − VD)`, affine in
`x` with slope `g(VD)`. And `G` returns exactly `g(VD)`. So `G` *is* `di/dx`,
exactly. Measured by central difference at `_vlim = 0.55` V:

| node v | `G[0,0]` | `d(i)/dx[0]` | equal? | `g` at the node voltage |
|---|---|---|---|---|
| 0.50 | 6.769891e-04 | 6.769891e-04 | yes | 9.779492e-05 |
| 0.55 | 6.769891e-04 | 6.769891e-04 | yes | 6.769891e-04 |
| 0.70 | 6.769891e-04 | 6.769891e-04 | yes | 2.245835e-01 |
| −2.00 | 6.769891e-04 | 6.769891e-04 | yes | 9.485184e-47 |

**Why this matters for the design.** If the defect were "J ≠ df/dx", native PCNR
would be a *correctness* fix and would rank accordingly. It is not. The residual
and Jacobian are mutually consistent; what is questionable is the **choice of
linearisation point** — `VD` comes from a limiter that is a function of the
previous iterate and of the order devices are visited. So native PCNR is a
**robustness and clarity** change (convergence behaviour, no hidden shared
state), not a correctness fix, and the case for it has to be made on those terms.

**But the exactness is per-iteration, and that is the whole subtlety.** Within
one linear solve, `G` is exactly `di/dx` of the residual *as implemented*, so the
solve is consistent. Across iterations `VD` moves, so the iteration is **not**
Newton's method on the original nonlinear system — it is a modified Newton on a
sequence of tangent points chosen by the limiter. Reading "`G` is exactly `di/dx`"
as "so limiting is a true Newton method" would be wrong, and it is the reading the
table invites.

That is precisely the quantity PCNR makes explicit: it promotes the tangent point
from a hidden per-device variable, updated as a side effect between iterations,
into an unknown the solver carries and converges. Stated that way the method's
claim is legible — **PCNR does not change the equations, it changes what is
allowed to be implicit** — and its benefit is correspondingly about convergence
robustness, not about getting a different answer.

The last row is also worth keeping in view: at −2 V the companion presents a
conductance of 6.77e-04 S where the true diode has 9.49e-47 S. That is not a bug
— it is the tangent line at `VD`, and Newton is entitled to it — but it means a
stale `VD` is silently a *plausible* matrix, never an obviously broken one. This
is exactly how the gate 13-6 defect survived three gates.

---

## 3. How pycircuit assembles today

Every element implements stamp methods — `G`, `C`, `i`, `q`, `u`, `dudt`, `CY` —
returning its own small dense matrix or vector over its own terminals. The
circuit sums them (`circuit.py:1337`):

```python
def G(self, x, epar=defaultepar, params_tree=None):
    return self._add_element_submatrices('G', x, (epar,), params_tree=params_tree)
```

`_add_element_submatrices` (`circuit.py`) is not a simple loop. It carries:

- an **evaluation-group** mechanism, so a toolkit can evaluate many instances of
  one class in bulk (`toolkit.batched_contributions`);
- a **sparse** path (`toolkit.build_sparse`) accumulating `(data, rows, cols)`;
- a **dense** path with a single scatter at the end;
- hoisted `hasattr`/`getattr` probes, deliberately, for speed.

Anything that changes *which elements are stamped* has to be threaded through all
four. **This, and not the device interface, is where the work lives.**

### Two limiting conventions, now one

`Semiconductor.limit` (`semiconductors.py:133`) **returns** a limited vector and
is state-free. `elements.Diode.limit` used to mutate `self._vlim` instead; both
now route through `_limiting.py`. `Diode` still *keeps* `_vlim` because its
non-autodiff `G`/`i` need a linearisation point.

### The autodiff split, which halves the problem

`Diode.G` has two branches:

```python
if self.toolkit.supports('autodiff'):
    return self.toolkit.jacobian(self.eval_i_pure, x, ...)   # pure in x
... else: use self._vlim                                      # stateful
```

**On an autodiff toolkit there is no `_vlim` and no hidden state at all.** The
JAX transient (`jaxtransient.py`) never limits — it uses a damped Newton with a
0.5 V step cap. So everything below about stale linearisation points applies to
the **numeric toolkit only**.

---

## 4. How PCNR is bolted on today, and its two bills

`pcnr.augmented_system` builds the enlarged system **by subtraction**:

```python
J_mm = np.array(circuit.G(x, epar), dtype=float)     # every device, diodes included
...
J_mm[np.ix_(rows, rows)] -= element.G(sub, epar)     # take the diodes back out
```

then re-stamps each junction at its own `v_lim` via `didv`. This is sound —
element stamps are additive, so "everything except the PCNR devices" is
expressible as a subtraction — and it was a deliberate win: **PCNR needed no
change to MNA assembly at all.** It bills twice for that.

### Bill 1 — cost (gate 13-4)

Measured on 60 diodes: classic assembly (`circuit.i` + `circuit.G`) 1115 µs,
`augmented_system` 2965 µs, of which the subtract loop is only **14%**. The cost
is not the subtracting; it is **running the full ordinary assembly first and then
doing more work on top**. Per Newton iteration this came out **+60% to +80%**.

Two optimisations were tried and neither closed it: exploiting the rank-one
collapse instead of a dense `J_ml @ J_lm`, and skipping the `(n,k)`/`(k,n)`
allocations. Both are in the shipped code; the ratio did not move. The cost is
structural to the layer.

### Bill 2 — a value that is not what its type says (gate 13-6)

On the PCNR path `cir.G(x)` returns a matrix **with no diode conductance in it**,
because `_vlim` is frozen at whatever it was first set to — PCNR never calls
`limit`. Measured on a mains rectifier: 1 `limit` call against 2283 `G` calls,
`_vlim` stuck at 0.0 V while the junction swung −18.47 to +0.75 V.

It cancels inside `augmented_system` — added by `cir.G(x)`, subtracted by the
per-device loop, *the same wrong number twice* — so the solved system is right
and every converged voltage is right. **The correctness of the solve does not
depend on `_vlim` being meaningful; it depends only on the two calls agreeing.**
Nothing in the solve can therefore detect that it is garbage.

The consequence is that **`cir.G(x)` on the PCNR path is not a Jacobian. It is
half of an expression** — and it has the type, shape and plausible magnitude of a
Jacobian, which is what makes it dangerous. It escaped through the one consumer
that used it for something other than the solve (the transient step controller,
forming `J⁻¹·Eg`), and produced steps 6.6× too large.

**Both bills are the same design decision seen from two sides.** Native PCNR is
the removal of that decision.

---

## 5. What "native" means — three designs

### Design A — assembly takes an exclusion set

```python
circuit.G(x, epar, skip=frozenset_of_instances)
```

`_add_element_submatrices` omits those instances; PCNR then stamps `didv`
directly. Nothing is built and destroyed.

- **Device contract:** unchanged. No device is touched.
- **Where the work is:** threading `skip` through the batched path, the sparse
  path, the dense path, and the evaluation groups — including
  `toolkit.batched_contributions`, whose JAX implementation vmaps a whole class
  at once and would need to mask or partition the group.
- **Cost:** reaches the paper's cost. Assembly walks *fewer* elements than the
  classic path, so PCNR assembly should be marginally cheaper than classic, not
  60-80% dearer.
- **Does it fix bill 2?** **No.** `cir.G(x)` with no `skip` still returns the
  same misleading matrix. The trap survives, guarded by discipline.
- **The specific risk, which is not the loop.** `_map_indices_2d` and the
  evaluation groups are **precomputed from the element list**. A runtime `skip`
  either invalidates those caches or has to be baked into a second set of them.
  Getting this wrong is not a crash — it is a scatter into the wrong entries,
  which looks like a subtly wrong Jacobian. Anyone starting here should decide
  *first* whether `skip` is a construction-time partition (two cached index maps)
  or a per-call argument, and the answer is probably the former.
- **Granularity.** `skip` is per *instance*, but a device could in principle have
  some junctions PCNR-managed and some not. Nothing needs that today, and
  per-instance is the cheaper contract; note it so the choice is deliberate.

### Design B — devices declare linear and nonlinear parts separately

The device contract grows: a device with PCNR junctions provides its nonlinear
contribution only through `pcnr_i`/`pcnr_didv`, and its `G`/`i` cover only what
is left. Assembly asks for the right form and there is nothing to subtract.

- **Device contract:** changed, for every nonlinear device.
- **Cost:** the paper's, same as A.
- **Does it fix bill 2?** **Yes, structurally** — there is no node-voltage
  nonlinear stamp anywhere, so no stale one to read.
- **Price:** this is the full cross-cutting change. See §7.

### Design C — an assembly mode threaded through `epar`

Devices consult a mode flag and stamp zero for their nonlinear part under
`mode='pcnr'`.

- **Device contract:** changed, but only by adding a branch.
- **Assessment:** **not recommended.** It puts solver state into the environment
  parameters, which are meant to describe the physical environment (temperature),
  and it makes every device's stamp silently dependent on a global. It also
  reintroduces exactly the failure mode of §4 bill 2: a caller that forgets the
  mode gets a plausible wrong matrix. **Reconsider if** a future need makes
  several assembly variants necessary, at which point the mode belongs on an
  explicit assembly context object, not on `epar`.

### The narrow fix that is *not* native PCNR, and should be considered first

Bill 2 alone can be closed without any of the above: have `Diode.G`/`Diode.i`
consult `_vlim` **only when limiting is active**, and use the node voltage
otherwise. The cancellation in `augmented_system` is invariant to what `G`
returns provided both calls agree, so PCNR is unaffected; classic limiting is
unaffected because the flag is on there.

This is a change to one device's state dependence rather than to the assembly
contract, and it makes `cir.G(x)` an honest Jacobian for every caller.

**It is one to two orders of magnitude cheaper than native PCNR and closes the
hazard that has actually cost debugging time.** It does nothing for bill 1.
*(This argument is static reasoning about the cancellation; it has not been
implemented or measured.)*

---

## 6. The device contract under Design B

Today, `Diode` declares:

```python
pcnr_junctions = ((0, 1),)          # (anode, cathode) indices into terminals

@staticmethod
def pcnr_i(v, params, epar, toolkit): ...      # terminal currents from its OWN v
@staticmethod
def pcnr_didv(v, params, epar, toolkit): ...   # d(terminal currents)/dv
```

Native support needs this from **every** nonlinear device. The current inventory:

| device | file | has `pcnr_*` | junctions | charge |
|---|---|---|---|---|
| `Diode` | elements.py | **yes** | 1 | no |
| `BJT` | semiconductors.py | no | 2 (B-E, B-C) | yes |
| `JFET` | semiconductors.py | no | 2 | yes |
| `ZenerDiode` | semiconductors.py | no | 0 *(deliberate)* | no |
| `Varactor` | semiconductors.py | no | — | yes |
| `VSwitch`, `ISwitch` | elements.py | no | — | no |
| `Idtmod`, `BSource`, `VCVS_limited` | elements.py | no | — | no |

Three structural problems appear here that the single-diode implementation never
had to face:

**(a) Multi-junction devices.** A BJT has two junctions sharing the base. The
existing `junctions` tuple carries a `move` index precisely because limiting both
by adjusting the anode would have the second undo the first. Under PCNR the two
become two unknowns — fine in principle, but `didv` is currently a 2-vector for a
2-terminal branch, and a BJT junction's current depends on *both* junction
voltages. The `lim/lim` block **stays the identity regardless** — the residual row is
`v_Dk − (e_a − e_b) = 0`, whose derivative with respect to `v_Dj` is `δ_jk` no
matter how the device's currents couple. What changes is the sparsity either
side: column `k` of `J_ml` gains a nonzero per device terminal (three for a BJT,
not two), while row `k` of `J_lm` still has exactly two. The Schur correction
stays a **sum of `k` rank-one terms**, now touching six entries each instead of
four — still `O(k)`, so this is a code generalisation and **not** a performance
concern. `pcnr.schur_reduce` hardcodes the two-terminal shape today.

**(b) Charge storage is currently refused.** `augmented_system` raises
`NotImplementedError` if a participating device has non-zero `C` — because the
charge term would also have to move to the limited unknown, and does not. Every
device in the table that matters for real circuits (BJT, JFET, Varactor, and any
future MOSFET) **has charge**. So native PCNR for anything beyond a diode
requires solving the charge problem first. **This is the single largest piece of
work in this document and it is easy to miss**, because the diode-only prototype
never encounters it.

**(c) `ZenerDiode` declares no junctions deliberately.** A forward-junction
limiter would step through the breakdown knee, which is a second exponential
running the other way; it gets a direction-agnostic exponent clamp instead. Under
native PCNR "no junctions" must keep meaning "not a PCNR device", not "a PCNR
device with zero unknowns". Recorded because it has already been mistaken for an
omission once.

---

## 7. Blast radius

### Directly changed

| area | Design A | Design B |
|---|---|---|
| `circuit.py` `_add_element_submatrices` / `_add_element_subvectors` | **yes** — skip-set through 4 code paths | yes |
| `toolkit.py` `batched_contributions` | **yes** — group masking/partition | yes |
| JAX `batched_contributions` (vmap per class) | **yes** | yes |
| `pcnr.py` | simplifies — subtract loop deleted | simplifies |
| every nonlinear device | no | **yes** |
| `semiconductors.py` `Semiconductor` base | no | **yes** |

### Consumers of `cir.G` / `cir.i` that must be audited either way

Enumerated from the tree:

- `dcanalysis.py:122,126` — the DC Newton
- `transient.py:981` (`_residual_and_jacobian`, coupled path), `:1562`, `:1588`
  (`solve_timestep`, classic path), and the PCNR timestep
- `analysis_ss.py:247,495` — small-signal / steady-state
- `feedback.py:20,121,177,178` — loop-gain analysis, which builds a *modified*
  circuit and takes differences of `G`
- `shooting.py:136` — periodic steady state
- `jaxtransient.py:226,252,468`

`feedback.py` deserves specific attention: it computes `G` of a circuit and of a
loop-broken copy and subtracts. Under any design where `G` means something
different depending on solver mode, **that subtraction becomes a question**, and
it is not obvious the two calls would be in the same mode.

### Not affected

- **Symbolic and DDD paths.** They consume `G` symbolically where devices are
  linear or already differentiated; no limiting is involved.
- **AC analysis.** Linearisation happens at a converged operating point; there is
  no limiter in the loop.
- **The JAX transient**, on its own terms — it uses autodiff Jacobians and a
  damped Newton, never `_vlim`. But if Design B changes what `G` *contains*, the
  JAX assembly consumes `circuit.G` too and would need the PCNR contribution
  routed in.

### Tests

`test_pcnr.py`, `test_dc_pcnr.py`, and the gate 13-6/13-7 regression tests all
assert against the *layer's* structure — in particular that PCNR and classic
limiting agree step-for-step on the standard path, and that `augmented_system` is
called. Design B removes `augmented_system`'s subtract loop, so the call-counting
test in gate 13-7 needs a different probe.

---

## 8. Recommendation

**Do the narrow fix from §5 now; do not do native PCNR yet.**

The two bills are **independent**, and that is the useful structure — the linear
list below obscures it. Pick by which bill is actually being paid:

| what hurts | the fix | size |
|---|---|---|
| bill 2 only (wrong Jacobian to other callers) | narrow fix — `Diode.G` stops reading stale `_vlim` | one device, hours |
| bill 1 only (assembly cost) | Design A — `skip` set through assembly | assembly internals, days |
| both | A + narrow fix | as above, additive |
| both, plus a clean device contract | Design B | cross-cutting, weeks |

**The end state worth aiming at is "A + narrow fix", not B.** That was not obvious
before writing this: A reaches the paper's cost with no device changes, and the
narrow fix closes the hazard structurally, so B's remaining advantage is contract
tidiness alone — which is not worth a cross-cutting change to every nonlinear
device on an opt-in path.

1. **Close bill 2 (`Diode.G`'s state dependence).** Small, local, testable, and
   it removes a trap that has already cost real debugging time. It is the only
   item here with a demonstrated cost of *not* doing it.
2. **Leave bill 1 (cost) parked.** PCNR is off by default; a 60-80% per-iteration
   penalty on an opt-in path is a real number but not a bill anyone is paying.
3. **If native PCNR is wanted, Design A first, then B if needed.** A reaches the
   paper's cost with **no device changes at all**, which was not obvious and is
   the main finding of this study. B's extra benefit over A is that it closes
   bill 2 structurally — but the narrow fix already closes bill 2, which removes
   most of B's remaining case.
4. **Nothing beyond a diode is reachable until charge storage is solved** (§6b).
   Any plan that says "extend PCNR to BJTs" is really a plan to move charge terms
   onto limited unknowns, and should be costed as that.

### The option this document nearly failed to list

**Removing PCNR entirely is a legitimate fourth answer**, and a survey that only
offers ways to invest more is not a survey. The case for removal: it is off by
default; it costs 60-80% per iteration; it converges no more often and in the same
iteration counts; it produces identical answers and, on the standard transient
path, identical steps; its one demonstrated effect on this repository was a defect
it introduced (gate 13-6). The clash it prevents has been investigated here once
and **did not reproduce**.

The case against removal: it is implemented, tested and documented, so carrying
cost is now low; the clash is real in the literature and the circuits that show it
(tightly coupled multi-diode nets) are ones this repository does not yet have; and
`pcnr.py` is where the correct-Jacobian knowledge now lives, so deleting it would
re-scatter that.

**On balance: keep, do not extend.** But the decision should be taken knowingly
rather than by default, and if §6b (charge storage) is ever costed honestly, the
removal option deserves re-reading first.

### Reconsider-ifs

- **PCNR becomes the default** → the layer is the wrong shape; Design A becomes
  necessary and bill 1 becomes a bill someone pays.
- **A second consumer of `cir.G` appears on the PCNR path** → bill 2 goes from
  latent to active; the narrow fix becomes urgent rather than advisable.
- **A BJT or MOSFET circuit needs PCNR to converge** → §6b is the project, and
  this document underestimates it; re-survey before planning.
- **The numeric toolkit is retired in favour of autodiff everywhere** → `_vlim`
  disappears, bill 2 disappears with it, and PCNR's remaining case is only the
  device-clash argument, which has never been reproduced on a real circuit here
  (the clash hypothesis was investigated once and did not reproduce).

---

## 9. How this would have to be tested

Gate 13-6 is the reason this section exists. That defect made `cir.G(x)` return a
matrix with **no diode conductance at all**, and it passed every DC test, three
gates, and every assertion on a converged voltage — because the wrong value
cancelled inside the one expression that used it. **Tests that check answers
cannot see assembly defects that cancel.** A native-PCNR change is entirely a
change to how assembly is composed, so it is squarely in that blind spot.

What actually catches this class:

1. **Assert on the matrix, not on the solution.** Compare the assembled Jacobian
   against a finite difference of the assembled residual, per device and at the
   circuit level. This is the check that would have caught 13-6 immediately, and
   `companion_consistency.py` is a one-device instance of it.
2. **Assert on call counts.** `test_gate_13_7_...` counts calls precisely because
   the two paths agree to within truncation error, so waveform comparison cannot
   distinguish "ran" from "was silently skipped". Any `skip`-set design needs the
   same: that the skipped elements were *actually* skipped.
3. **Compare against an independent reference, never the integrator with itself.**
   Gate 13-7b failed as declared because it compared the two paths to each other;
   the question was only settled by measuring both against a tight reference. Any
   accuracy claim here must be against closed forms or a tight run of a *different*
   path.
4. **Verify by injection.** Every regression test above was confirmed to fail
   against the defect reinstated deliberately. For an assembly change, injection
   means: skip an element that should not be skipped, and confirm something goes
   red.

## 10. What this document does not establish

- **That PCNR is worth having.** On measurement it converges where limiting does,
  in the same iteration counts, to the same answers, and — on the standard
  transient path — takes *identical* steps. Its case rests on robustness
  arguments that have not been demonstrated to bite on a circuit in this
  repository.
- **That Design A's assembly threading is cheap.** It is asserted to be the
  smaller change because it leaves the device contract alone. The batched/JAX
  group path is the risk, and it has not been prototyped.
- **Any number for Design B's cost.** No implementation exists to measure.
