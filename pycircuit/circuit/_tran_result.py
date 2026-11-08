"""A transient's result (`TransientResult`) and what `save` names
(`resolve_save`).  A theme of `Transient` (see `transient.py`).

A run keeps every unknown unless `save` names what to keep: a node's
voltage by its name, a current by `i(term)` with the string `res.i()`
takes.  A source's or an inductor's current is an unknown (a branch row);
a device's terminal current is made from the rows the device reads, at
each point's time and with its stateful limiters at each point.  A device
holding charge is refused: its capacitive current is the integrator's
derivative of its own charge, which the run does not keep.
"""

import re

import numpy as np

from pycircuit.circuit._limiting import devices_at, stateful_limiters
from pycircuit.circuit.analysis import CircuitResult
from pycircuit.circuit.circuit import Circuit, SubCircuit, defaultepar

_CURRENT = re.compile(r'^i\((.+)\)$')


def holds_charge(c):
    """Whether `c` (an element, or a subcircuit through its elements) has a
    charge: its class gives `q` or `C` other than `Circuit`'s zeros."""
    if isinstance(c, SubCircuit):
        return any(holds_charge(e) for e in c.elements.values())
    return type(c).q is not Circuit.q or type(c).C is not Circuit.C


def device_of(cir, term):
    """`(element, its terminal, its rows in the run's x)` for a terminal
    named through the hierarchy, as `SubCircuit.extract_i` walks it; rows
    None for the circuit itself."""
    c, rows = cir, None
    while isinstance(c, SubCircuit) and c.get_terminal_branch(term) is None and '.' in term:
        first, term = term.split('.', 1)
        nm = np.asarray(c.elementnodemap[first], dtype=np.intp)
        rows = nm if rows is None else rows[nm]
        c = c.elements[first]
    return c, term, rows


def _refusal(c, term, path):
    """Why `path`'s current cannot be made from the run's rows, or None."""
    if term not in c.nodenames:
        return f'i({path}): no terminal {term!r}'
    hidden = c.hidden_state_elements()
    if hidden:
        return (f'i({path}): {", ".join(hidden)} carries state the result does not keep '
                '(a delay line\'s history), so its current cannot be made from the run\'s rows')
    if holds_charge(c):
        return (f'i({path}): the device holds charge, and its capacitive current is the integrator\'s '
                'derivative of that charge, which the run does not keep; only currents of branches '
                '(sources, inductors) and of devices without charge can be had from a transient')
    return None


class Saved:
    """What `save` resolved to: the names in order and `rows`, the run's
    rows to keep (sorted, the reference node's always among them)."""

    def __init__(self, names, rows):
        self.names, self.rows = tuple(names), rows


def resolve_save(cir, names, irefnode):
    """`save`'s names against the circuit: a `Saved`, or ValueError naming
    the first that does not resolve -- before the run spends anything."""
    if isinstance(names, str):
        raise TypeError(f'save takes a list of names, not the string {names!r}')
    rows = {int(irefnode)}
    for name in names:
        m = _CURRENT.match(str(name).strip())
        if m is None:
            try:
                rows.add(int(cir.get_node_index(name)))
            except ValueError:
                raise ValueError(f'save: {name!r} is not a node of the circuit (a current is i(term))') from None
            continue
        term = m.group(1).strip()
        try:
            bs = cir.get_terminal_branch(term)
            if bs is not None:
                rows.add(int(cir.get_branch_index(bs[0])))
                continue
            c, local, sub = device_of(cir, term)
        except (KeyError, ValueError):
            raise ValueError(f'save: {name!r} names no terminal of the circuit') from None
        why = _refusal(c, local, term)
        if why is not None:
            raise ValueError('save: ' + why)
        rows.update(int(r) for r in (np.arange(cir.n) if sub is None else sub))
    return Saved(names, np.array(sorted(rows), dtype=np.intp))


class _Rows:
    """The run's x seen through the saved rows: `[k]` and `[ks]` by the
    run's row numbers, a row not saved refused naming `save`."""

    def __init__(self, x, rows):
        self.x = x
        self.at = {int(r): j for j, r in enumerate(rows)}

    def _j(self, k):
        try:
            return self.at[int(k)]
        except KeyError:
            raise ValueError(f'row {int(k)} was not saved: save=[...] names what a run keeps') from None

    def __getitem__(self, k):
        if np.ndim(k) == 0:
            return self.x[self._j(k)]
        return self.x[[self._j(i) for i in np.asarray(k).ravel()]]


class TransientResult(CircuitResult):
    """A transient's result: every unknown, or the rows `save` kept
    (`rows`, the run's row numbers of `x`'s rows; None: all).

    `i(term)` of a branch is its row.  Of a device's terminal it is made
    from the device's rows: its static current with the sources at each
    point's time, a device with stateful limiters read at each point
    (`devices_at`).  A device holding charge, or carrying state the result
    does not keep, is refused (`_refusal`)."""

    def __init__(self, circuit, x, sweep_values, rows=None, saved=None, epar=None):
        super().__init__(circuit, x, xdot=None, sweep_values=sweep_values, sweep_label='time', sweep_unit='s')
        self.rows = rows
        self.saved = saved
        self.epar = defaultepar if epar is None else epar

    def _x(self):
        return self.x if self.rows is None else _Rows(self.x, self.rows)

    def v(self, plus, minus=None):
        if self.rows is None:
            return super().v(plus, minus)
        result = self.circuit.extract_v(self._x(), plus, minus)
        ylabel = f'v({plus},{minus})' if minus is not None else f'v({plus})'
        return self.build_waveform(result, ylabel, 'V')

    def i(self, term):
        cir, x = self.circuit, self._x()
        if cir.get_terminal_branch(term) is not None:
            result = cir.extract_i(x, term)
        else:
            c, local, sub = device_of(cir, term)
            if c.get_terminal_branch(local) is not None:
                result = Circuit.extract_i(c, x[np.arange(cir.n) if sub is None else sub], local)
            else:
                why = _refusal(c, local, term)
                if why is not None:
                    raise ValueError(why)
                result = self._device_current(c, local, x[np.arange(cir.n) if sub is None else sub])
        return self.build_waveform(result, f'i({term})', 'A')

    def _device_current(self, c, term, subx):
        """`Circuit.extract_i`'s terminal current, `i(x) + u(t)`, at each
        point's own time, the device at each point where it has stateful
        limiters."""
        k = c.get_node_index(c.nodenames[term])
        t = np.asarray(self.sweep_values, dtype=float)
        epar = self.epar
        ## (`stateful_limiters` finds a circuit's elements, not the circuit)
        lims = [c] if getattr(c, 'reset_limit_state', None) is not None else stateful_limiters(c)
        if lims:
            out = np.empty(t.size)
            for j in range(t.size):
                xj = subx[:, j]
                with devices_at(c, xj, epar, lims):
                    out[j] = c.i(xj, epar)[k] + c.u(t[j], epar=epar, analysis='tran')[k]
            return out
        static = c.i(subx, epar)[k]
        if type(c).u is Circuit.u:
            ## (no sources: `Circuit.u`'s zeros, as `extract_i` adds them)
            return static + c.u(0.0, analysis='tran')[k]
        return static + np.array([c.u(tj, epar=epar, analysis='tran')[k] for tj in t])
