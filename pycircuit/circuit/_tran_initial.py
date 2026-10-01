"""The run's initial state: initial conditions and the operating point.  A
theme of `Transient` (see `transient.py`).
"""

import numpy as np

from pycircuit.circuit.analysis import (
    NoConvergenceError,
    SingularMatrix,
)


class _InitialState:
    """The run's initial state: initial conditions and the operating point.  A
    theme of `Transient` (see `transient.py`)."""

    def _initial_state(self, refnode):
        """The `uic=True` starting vector: zeros, plus whatever `ic` names.

        STAGE 10.3.  A vector of zeros alone makes a class of circuit
        unsimulable rather than merely inconvenient: **an LC tank at zero is at
        an equilibrium** and will sit there forever, and a latch at zero is on
        its metastable point.

        The analysis-level `ic` names node voltages only.  Element-level
        initial conditions -- SPICE's ``L ... IC=i`` and ``C ... IC=v``, and a
        state element's seed -- are applied after it, by `_apply_element_ics`
        (a ``C``'s constrains a DIFFERENCE of two node unknowns and is solved
        as a spanning tree in `_apply_voltage_ics`).

        The starting vector is deliberately NOT made consistent with the circuit
        equations. Under `uic` there is no operating point by definition; the
        first Newton solve at `t = h` sees these values as history and produces a
        consistent solution from them, which is what SPICE does too.

        History: `doc/transient_history.md`, `Transient._initial_state`.
        """
        n = self.cir.n
        ## numpy, not self.toolkit: this vector is built and MUTATED before
        ## the loop starts, and the JAX backend shares these methods verbatim
        ## (P12) -- a traced array cannot take item assignment, and pre-loop
        ## state is numpy on both backends anyway (the caller converts).
        x0 = np.zeros(n)

        ic = self.par.ic
        if not ic:
            return self._apply_element_ics(x0, refnode)

        irefnode = self.cir.get_node_index(refnode)
        for node, value in dict(ic).items():
            try:
                idx = self.cir.get_node_index(node)
            except ValueError:
                raise ValueError(
                    "ic names node %r, which is not in the circuit. Nodes are: "
                    "%s" % (node, ', '.join(str(nd) for nd in self.cir.nodes)))
            ## Naming the reference node is not a harmless no-op -- it is a
            ## statement the solver cannot honour, since that node is held at
            ## zero by construction, so silently dropping it would leave the
            ## caller believing an initial condition was applied.
            if idx == irefnode:
                raise ValueError(
                    "ic sets node %r, which is the reference node and is held "
                    "at 0 V by construction" % (node,))
            x0[idx] = value

        return self._apply_element_ics(x0, refnode)

    def _apply_element_ics(self, x0, refnode):
        """Write each element's ``ic`` instparam into its own branch rows.

        STAGE 10.3.  Only elements whose initial condition is a branch CURRENT
        can be handled this way -- `L` today. The row is found through the
        circuit's recorded instance-to-branch span, NOT by searching
        `self.branches` for a matching `Branch`: that search is ambiguous for
        parallel elements, whose branches compare equal, so an initial current
        given to the second of two parallel inductors would land on the first
        one's unknown with nothing to indicate it.
        """
        cir = self.cir
        elements = getattr(cir, 'elements', None)
        if not elements:
            return x0

        for name, element in elements.items():
            ## A nested subcircuit owns a span covering its children's branches,
            ## so an `ic` inside one cannot be placed by this flat walk. Detected
            ## and refused rather than skipped: skipping would accept the
            ## parameter and ignore it.
            if getattr(element, 'elements', None):
                if self._descendant_has_ic(element):
                    raise NotImplementedError(
                        "element initial conditions inside a subcircuit (%r) are "
                        "not supported: the branch rows of a nested instance are "
                        "not individually resolvable yet. Move the element to the "
                        "top level, or set the node voltage with the analysis's "
                        "`ic` instead." % name)
                continue

            ## State-flavoured ICs (Idt/Idtmod) live on the element's own
            ## private row -- neither a branch current nor a node-voltage
            ## difference -- and the element hands back `(local_row, value)`
            ## pairs it has already converted (e.g. wrapped into the modulus
            ## range).  `elementnodemap` maps its local x-indices, private
            ## nodes included, onto this circuit's rows.
            ##
            ## GATED ON `state_ic`, NOT ON A PARAMETER SPELLED `ic`: a
            ## generated model's state seed may be called anything else
            ## (`x0`, `phi0`), and a NAME check would let `uic=True` start that
            ## state at zero -- no error, no warning, a wrong waveform.
            ## `IC_KIND == 'state'` and `state_ic` are installed under the
            ## same condition (hdl.py, `state_meta['dc_pins']`), so this is
            ## the same question asked of the thing that answers it.
            ## History: `doc/transient_history.md`, `Transient._apply_element_ics`.
            if getattr(element, 'IC_KIND', 'current') == 'state':
                if not hasattr(element, 'state_ic'):
                    continue
                rows = cir.elementnodemap[name]
                for local_row, value in element.state_ic():
                    x0[rows[local_row]] = value
                continue

            ic = getattr(getattr(element, 'iparv', None), 'ic', None)
            if ic is None:
                continue

            ## Voltage-flavoured ICs constrain a difference of two node unknowns
            ## and are solved together, after this loop.
            if getattr(element, 'IC_KIND', 'current') != 'current':
                continue

            rows = cir.instance_branch_indices(name)
            if len(rows) != 1:
                raise ValueError(
                    "%r declares an ic but owns %d branch rows; an initial "
                    "condition is only meaningful for an element with exactly "
                    "one branch current" % (name, len(rows)))
            x0[rows[0]] = ic

        return self._apply_voltage_ics(x0, refnode)

    def _apply_voltage_ics(self, x0, refnode):
        """Solve the capacitor initial voltages as a spanning tree.

        STAGE 10.3.  A capacitor has no state variable of its own -- `q` is
        derived from the node voltages -- so `C ... IC=v` cannot be assigned
        anywhere. It constrains ``v(plus) - v(minus) = v``, a DIFFERENCE of two
        unknowns, and a set of such constraints is a system rather than a list of
        assignments. See `doc/initial_conditions.md` sec. 4a.

        Each constraint is an edge; the reference node and anything the
        analysis-level `ic` named are seeds; each connected component is walked
        breadth-first from a seed, and a node reached twice must agree both times.

        **A component with no seed raises.** Its voltages are determined only up
        to a constant, so infinitely many assignments satisfy what was asked for,
        and the absolute values reach the output waveform. Choosing one silently
        is the defect shape this stage exists to avoid.
        """
        cir = self.cir
        elements = getattr(cir, 'elements', None)
        if not elements:
            return x0

        edges = {}
        for name, element in elements.items():
            if getattr(element, 'elements', None):
                continue
            ic = getattr(getattr(element, 'iparv', None), 'ic', None)
            if ic is None or getattr(element, 'IC_KIND', 'current') != 'voltage':
                continue
            terms = cir.term_node_map[name]
            p = cir.get_node_index(terms['plus'])
            m = cir.get_node_index(terms['minus'])
            if p == m:
                if ic != 0.0:
                    raise ValueError(
                        "%r has both terminals on the same node but ic=%g; that "
                        "constrains 0 == %g" % (name, ic, ic))
                continue
            edges.setdefault(p, []).append((m, +float(ic), name))
            edges.setdefault(m, []).append((p, -float(ic), name))

        if not edges:
            return x0

        ## Seeds: the reference node, plus whatever the node-level `ic` fixed.
        ## `refnode` is passed rather than read from `self.irefnode`, which is
        ## set by the caller a few lines earlier -- a hidden ordering dependency
        ## that would break silently if this were ever called first.
        seeded = {cir.get_node_index(refnode)}
        for node in dict(self.par.ic or {}):
            seeded.add(self.cir.get_node_index(node))

        ## Tolerance for the consistency check: relative to the largest voltage
        ## involved, so a chain at kilovolts is not judged by the same absolute
        ## slack as one at millivolts.
        scale = max([abs(v) for row in edges.values() for _, v, _ in row] + [1.0])
        tol = 1e-9 * scale

        assigned = set(seeded)
        for start in sorted(edges):
            if start in assigned:
                continue
            ## Walk this component to see whether it contains any seed at all.
            comp, stack = set(), [start]
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                stack.extend(nb for nb, _, _ in edges.get(cur, ()))
            if not (comp & assigned):
                names = sorted({nm for nd in comp
                                for _, _, nm in edges.get(nd, ())})
                raise ValueError(
                    "the initial voltages on %s form a group with no connection "
                    "to ground or to any node given in `ic`, so they fix the node "
                    "voltages only up to a constant. Ground one node of the "
                    "group, or name one in the analysis's `ic`."
                    % ', '.join(repr(n) for n in names))

        ## Breadth-first assignment from every seed.
        from collections import deque
        queue = deque(sorted(assigned & set(edges)))
        while queue:
            cur = queue.popleft()
            for nb, dv, name in edges.get(cur, ()):
                value = x0[cur] - dv
                if nb in assigned:
                    if abs(x0[nb] - value) > tol:
                        raise ValueError(
                            "initial voltages are contradictory at node index %d: "
                            "%r implies %g V, but %g V was already established. "
                            "Check for a loop of capacitor ics that does not sum "
                            "to zero, or an ic that disagrees with the analysis's "
                            "`ic`." % (nb, name, value, x0[nb]))
                    continue
                x0[nb] = value
                assigned.add(nb)
                queue.append(nb)

        return x0

    def _descendant_has_ic(self, circuit, include_state=True):
        """Does anything below this instance carry a set ``ic``?

        ``include_state=False`` skips elements with ``IC_KIND='state'``
        (Idt/Idtmod): their ``ic`` pins the DC operating point per the LRM,
        so unlike ``L``/``C`` it is meaningful WITHOUT ``uic=True`` and must
        not trip the guard that rejects that combination.  The nested-
        subcircuit refusal in ``_apply_element_ics`` keeps the default, so a
        state ic buried where the flat uic walk cannot reach it still fails
        loudly instead of being dropped.
        """
        for element in getattr(circuit, 'elements', {}).values():
            ## As in `_apply_element_ics`: a state element declares its seed
            ## through `state_ic`, not through a parameter named `ic`.  Asking
            ## about `ic` here would make the guard UNDER-detect, which is the
            ## direction that silently drops an initial condition instead of
            ## refusing it.
            ## History: `doc/transient_history.md`, `Transient._descendant_has_ic`.
            kind = getattr(element, 'IC_KIND', 'current')
            if kind == 'state':
                carries = hasattr(element, 'state_ic')
            else:
                carries = getattr(getattr(element, 'iparv', None),
                                  'ic', None) is not None
            if carries:
                if include_state or kind != 'state':
                    return True
            if self._descendant_has_ic(element, include_state=include_state):
                return True
        return False

    def _solve_operating_point(self, refnode):
        """Solve the DC operating point that seeds the transient.

        A failure here **raises**.  A substituted vector of zeros would return a
        complete, plausible-looking waveform computed from a bias point that was
        never found, indistinguishable from a successful run: it does not fail, it
        lies.

        The inner `DC` is constructed from *this* analysis's configuration rather
        than from `DC`'s defaults -- the transient's toolkit, environment
        parameters, tolerances, solver and scaler -- so the operating point is
        solved at the same temperature, to the same accuracy, and with the same
        Newton strategy as every step that follows it.

        History: `doc/transient_history.md`, `Transient._solve_operating_point`.
        """
        from pycircuit.circuit.dcanalysis import DC

        ## Only forward what DC actually declares, so a parameter that exists on
        ## Transient but not on DC (e.g. `integrator`) does not raise a KeyError,
        ## and so this keeps working if either parameter list changes.
        dc_par_names = {p.name for p in DC.parameters}
        shared = {}
        for name in ('reltol', 'iabstol', 'vabstol', 'maxiter',
                     'bypass', 'bypasstol', 'epar', 'nrsolver', 'scaler'):
            if name not in dc_par_names:
                continue
            try:
                value = getattr(self.par, name)
            except (AttributeError, KeyError):
                continue
            if value is not None:
                shared[name] = value

        dc = DC(self.cir, toolkit=self.toolkit, refnode=refnode, **shared)
        try:
            return dc.solve().x
        except (NoConvergenceError, SingularMatrix) as exc:
            raise NoConvergenceError(
                "Transient could not find a DC operating point to start from: %s\n"
                "The transient has NOT been run.  Either fix the bias condition, or "
                "start deliberately from a known state:\n"
                "  * Transient(..., uic=True)  -- start from zeros (SPICE's 'use "
                "initial conditions'); note this is a Transient() argument, NOT a "
                "solve() one, or\n"
                "  * solve(x0=<vector>)        -- start from an operating point you "
                "supply.\n"
                "Both are explicit choices; substituting zeros silently is what this "
                "error replaced." % (exc,)) from exc
