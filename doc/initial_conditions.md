# Initial conditions — what works, what does not, and why

**Status 2026-08-02.** Node voltages and inductor currents can be set for `uic=True`.
Capacitor initial voltages cannot. This page says precisely where the boundary is, because
the two halves are separated by a real mathematical difference rather than by how much time
was available.

---

## 1. What ships

```python
tran = Transient(circuit, uic=True, ic={'out': 1.2})   # node voltages
circuit['L1'] = L('a', 'b', L=1e-6, ic=0.5)            # inductor current, amperes
```

`uic=True` means "do not solve an operating point; start from these values". Anything not
named starts at zero, as before.

**Why this is not a convenience feature.** `uic=True` used to mean a vector of zeros, and
that leaves a class of circuit *unsimulable* rather than merely awkward:

* an **LC tank at zero is at an equilibrium** — it stays there forever, so an oscillator
  could not be started at all;
* a **latch at zero sits on its metastable point** — it has no reason to fall either way.

`tests/test_initial_conditions.py` asserts the tank stays dead *without* `ic`, so the feature
cannot be mistaken for optional, and then that it oscillates at `1/(2π√(LC))` with it — the
frequency measured from zero crossings against the analytic value, so the test cannot be
satisfied by a circuit that merely wobbles.

## 2. What raises, and why silence was the alternative

Every one of these could have been a no-op. Each is a request the solver cannot honour, and
accepting it quietly would leave the caller believing an initial condition had been applied.

| situation | behaviour |
|---|---|
| `ic` names a node that does not exist | `ValueError`, listing the circuit's nodes |
| `ic` names the reference node | `ValueError` — it is held at 0 V by construction |
| `ic` or an element `ic` given without `uic=True` | `ValueError` naming both forms |
| an element `ic` inside a nested subcircuit | `NotImplementedError` naming the instance |
| an element declares `ic` but owns ≠ 1 branch row | `ValueError` |

The `uic` guard is worth singling out. SPICE's `.ic` **without** UIC constrains the operating
point and then releases it — a genuinely different feature, and one that is *not* implemented
here. Raising says which of the two is missing; ignoring would do neither and report neither.

## 3. The boundary: why `L` works and `C` does not

This is the substance of the page.

**An inductor's initial condition is a branch current, and its unknown already exists.**
MNA carries a row for every inductor current, so setting one is an *assignment*: find the
row, write the number.

**A capacitor's initial condition is a branch voltage, and no unknown corresponds to it.**
`C1 ... IC=2.5` constrains `v(plus) − v(minus) = 2.5` — a **difference** of two node unknowns,
not either of them. A set of such constraints does not assign values; it defines a system.

Concretely, three capacitors in a floating chain `a—b—c—d` with `d` grounded and initial
voltages `v_ab`, `v_bc`, `v_cd` determine every node. The same three capacitors with **no**
node grounded determine the nodes only up to a constant — the chain floats, and any solution
plus an offset is equally valid. Deciding which offset to take is a choice about the circuit,
not arithmetic, and picking one silently would be exactly the class of quiet wrong answer
this project keeps finding.

The general form is a spanning-tree problem: take the graph whose edges are the capacitors
carrying an `IC`, root each connected component (at ground where a component touches it), and
propagate voltages along a tree. Components not touching ground need a reference chosen or an
error raised; cycles need consistency checked, because `v_ab + v_bc ≠ v_ac` is a contradiction
a user can easily write.

**That is the deferred work**, and it is a self-contained piece: a graph traversal plus a
consistency check, with a decision to make about ungrounded components.

## 4. How an element's branch row is found

Recorded when the branch is created, not reconstructed later:

```python
branch_start = len(self.branches)
self.append_branches(*newbranches)
self._instance_branch_span[instancename] = (branch_start, len(self.branches))
```

**Reconstruction is silently wrong for parallel elements**, which is why the map exists.
`Branch.__eq__` compares node pairs, so two inductors between the same two nodes produce
*equal* branches and `branches.index()` returns the first for both. Measured before the map
existed: `get_branch_index` gave row 2 for each of two parallel inductors, so an initial
current given to the second would have landed on the first one's unknown with nothing
whatever to indicate it.

The span is stored in **branch-list coordinates, not solution-vector rows**. A branch's row
is `len(nodes) + offset`, and nodes are still being added while elements are, so a stored row
would go stale. `instance_branch_indices(name)` adds the offset at lookup time.

## 5. Not implemented

* **`C ... IC=`** — §3. The spanning-tree solve, plus the ungrounded-component decision.
* **`.nodeset`** — a genuinely different feature: a *hint* that seeds the DC solve and is then
  released, where `.ic` under UIC is a *starting value* that never is. It belongs with stage
  5's convergence-aid ladder, not here.
* **`.ic` without UIC** — constraining the operating point and releasing it; see §2.
* **Element `ic` inside a subcircuit** — the span of a nested instance covers all its
  children's branches. Resolving `'X1.L1'` means recursing into the child's own map with the
  parent's offset. Mechanical, but a second piece of work.

**Reconsider the deferrals if** a circuit needs a floating capacitor's initial voltage, a
starting current on an inductor inside a subcircuit, or a DC-convergence hint. None of the
three is expressible by naming node voltages, and none is a workaround away.

## 6. A note on the tests

The feature code needed no correction during this work. **Four test circuits did**, and every
one of them failed in a way indistinguishable from the feature being broken:

* two ideal `VCVS` cross-coupled directly — node voltages become *algebraically* determined,
  so `v1 = 9·v1` forces `v1 = 0` and there is no state for an initial condition to set;
* `VCVS` written as `(outp, outn, inp, inn)` — it is `(inp, inn, outp, outn)`, so each source
  *drove* the node it was meant to sense and pushed 1 A into the capacitor it should have
  been reading;
* `assert abs(v[0]) < 1e-9` — the result excludes `x0`, so the first sample is the first
  *stepped* point, exactly `0.5 A · 1e-13 s / 1 nF = 5e-5 V`. The assertion was wrong, not the
  code;
* a bare `SubCircuit()` assigned as an element — it must be a subclass with `terminals`,
  instantiated with the nodes it connects to, or `_instance_nodes` raises `KeyError`.

The first two both produce "v(1) sits at zero", which is also exactly what a broken `ic`
produces. That ambiguity is the reason each is recorded in the test docstrings rather than
quietly fixed: the next person to see a latch that will not start needs to know these three
explanations exist before concluding the fourth.
