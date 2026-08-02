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
circuit['C1'] = C('a', 'b', c=1e-9, ic=2.5)            # capacitor voltage, volts
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
| capacitor `ic`s form a group touching neither ground nor a node `ic` | `ValueError` — they fix the nodes only up to a constant |
| capacitor `ic`s contradict each other, or a node `ic` | `ValueError` naming the element and both values |
| a capacitor with both terminals on one node and `ic ≠ 0` | `ValueError` — it constrains `0 == ic` |
| `ic` or an element `ic` given without `uic=True` | `ValueError` naming both forms |
| an element `ic` inside a nested subcircuit | `NotImplementedError` naming the instance |
| an element declares `ic` but owns ≠ 1 branch row | `ValueError` |

The `uic` guard is worth singling out. SPICE's `.ic` **without** UIC constrains the operating
point and then releases it — a genuinely different feature, and one that is *not* implemented
here. Raising says which of the two is missing; ignoring would do neither and report neither.

## 3. The boundary: why `L` is an assignment and `C` is a solve

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

**That is what §4a implements**: a graph traversal plus a consistency check, with the
ungrounded-component question answered by raising.

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

## 4a. `C ... IC=` — investigated 2026-08-02, then implemented as designed

### What the investigation established

**The graph is directly constructible.** `term_node_map[instance]['plus'|'minus']` gives the
parent node and `get_node_index` gives its row, so no new plumbing is needed to know which
two unknowns a capacitor's `ic` constrains.

**There are no graph helpers anywhere in the package.** No traversal utilities, no networkx
dependency. The component search and the tree walk are new code, with no existing pattern in
the codebase to follow.

**A capacitor has no independent state, and this is the real reason.** `C.q(x)` is computed
purely from node voltages — `C·(v_plus − v_minus)`. So the obstacle is not merely that an
`ic` constrains a *difference*: even if one wanted to store an initial charge on the device,
**there is nowhere to put it.** SPICE can implement `C ... IC=` as a direct assignment
because its transient state for a capacitor *is* the branch voltage. pycircuit's MNA carries
no such variable. That is a structural difference between the two simulators rather than a
gap in effort, and it is why this feature is a solve rather than an assignment.

**And it unblocks nothing.** Measured on a series-capacitor circuit whose midpoint has no DC
path to ground:

| | result |
|---|---|
| `uic=False` | fails — singular Jacobian, *"'m' appears in no equation"* |
| `uic=True` | runs, 230 points, `v(m)` starting at 0 |

Unlike the node initial conditions, which made oscillators and latches simulable *at all*,
`C.ic` chooses a different starting state for a circuit that already runs. The need is real
— pre-charged capacitors, bootstrap capacitors, switched-capacitor circuits — but it is
expressiveness rather than capability, and the priority argument should say so.

### The algorithm

1. Collect capacitors carrying an `ic`; each is an edge between two node rows with a required
   difference.
2. Seed the known voltages: the reference node at 0, plus anything the analysis-level `ic`
   dict names.
3. Find connected components of that edge set.
4. Root each component at a seeded node and breadth-first assign
   `v(neighbour) = v(current) ± ic`.
5. On reaching an already-assigned node by a second edge, check the implied value agrees.

### Decisions, recorded because they are choices rather than derivations

**A component containing no seeded node raises.** Its voltages are determined only up to a
constant — the differences the user asked for are satisfied by infinitely many assignments.
Two alternatives were considered and rejected: silently choosing a reference, which is the
quiet-wrong-answer shape this project keeps finding, since the absolute node voltages appear
in the output waveform; and choosing one with a warning, which is the same thing with a line
of text. SPICE effectively does choose silently, and can afford to because its per-capacitor
state makes the question moot. **Reconsider if** a real circuit turns up whose floating
component genuinely has no natural reference and the error is merely obstructive — the fix
would be an explicit opt-in, not a change of default.

**The `ic` flavour is declared, not inferred.** `L`'s `ic` is a current and `C`'s is a
voltage, and they are handled by completely different mechanisms. Inferring from "does the
element own a branch" would work today and silently mis-handle the first element that breaks
the correlation, so each element states its own: `IC_KIND = 'current'` or `'voltage'`.

**Contradictions raise rather than being resolved by ordering.** A cycle whose `ic` values do
not sum to zero, two parallel capacitors disagreeing, a capacitor `ic` disagreeing with a
node `ic`, and a capacitor whose terminals are the same node with a nonzero `ic` are all
statements the caller cannot have meant. Silently taking the last one written would make the
result depend on element insertion order.

**Node ICs are applied first** and act as seeds, so a capacitor chain anchored by an explicit
node voltage propagates from it.

## 5. Not implemented

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
