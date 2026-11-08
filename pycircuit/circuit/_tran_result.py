"""A transient's result (`TransientResult`) and what `save` names
(`resolve_save`).  A theme of `Transient` (see `transient.py`).

A run keeps every unknown unless `save` names what to keep: a node's
voltage by its name, a current by `i(term)` with the string `res.i()`
takes.  A source's or an inductor's current is an unknown (a branch row);
a device's terminal current is made from the rows the device reads, at
each point's time and with its stateful limiters at each point.  A device
holding charge has its current made during the run, where `save` names it
(`DeviceCurrents`): its capacitive current is the integrator's derivative
of its own charge, which a run that keeps every unknown does not keep, so
a result's `i()` refuses it otherwise.
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
                'derivative of that charge, which a run keeps only for the currents `save` names '
                f"(save=['i({path})'])")
    return None


class Saved:
    """What `save` resolved to: the names in order, `rows`, the run's rows
    to keep (sorted, the reference node's always among them), and
    `charged`, the currents of devices holding charge made during the run
    (`(term, device, its terminal, its rows)` each)."""

    def __init__(self, names, rows, charged=()):
        self.names, self.rows, self.charged = tuple(names), rows, tuple(charged)


def _charge_refusal(c, path, charge):
    """Why a device holding charge cannot have its current made in this
    run (`charge`: the run's reason, or None), or None."""
    if charge is not None:
        return f'save: i({path}): the device holds charge, and ' + charge
    if c.periodic_states():
        return (f'save: i({path}): the device\'s state is rewrapped as the run goes (a periodic state), '
                'so its charge has no derivative its rows give')
    return None


def resolve_save(cir, names, irefnode, charge=None):
    """`save`'s names against the circuit: a `Saved`, or an error naming
    the first that does not resolve -- before the run spends anything.
    `charge`: why this run cannot make a charge-holding device's current
    (its method), or None."""
    if isinstance(names, str):
        raise TypeError(f'save takes a list of names, not the string {names!r}')
    rows, charged = {int(irefnode)}, []
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
        sub = np.arange(cir.n) if sub is None else sub
        why = _refusal(c, local, term)
        if why is not None and holds_charge(c) and local in c.nodenames and not c.hidden_state_elements():
            why = _charge_refusal(c, term, charge)
            if why is None:
                charged.append((term, c, local, sub))
                continue
            raise ValueError(why)
        if why is not None:
            raise ValueError('save: ' + why)
        rows.update(int(r) for r in sub)
    return Saved(names, np.array(sorted(rows), dtype=np.intp), charged)


class DeviceCurrents:
    """`save`'s currents of devices holding charge, made at each accepted
    point as the run goes: the device's static current, with its sources
    at the point's time and its stateful limiters at the point, plus its
    charge's derivative by the step's own companion, `iq_n = sum_k
    alphas[k] q_{n-k} + b iq_{n-1}` (`companion_coefficients`) on the
    device's own charges.  The run's `iq` row is the same sum over a node's
    devices, so the currents at a node close as its KCL does.  The start
    point's derivative is zero, as the run's ring starts: right at an
    operating point; at a start the run imposed (`uic`, `ic`, an `x0`) the
    device's capacitive current there is not one the run made, so its
    current there is its static part and KCL does not close at that point.

    A stage step (Radau IIA, `w` its `A`'s inverse's last row) hands its
    step and its stages 1 and 2 (`('stage', h, Y1, Y2, on_rows)`: on
    `srows`, the saved devices' rows, where `on_rows`, else whole): the
    derivative is the step's collocation derivative at its end on the
    device's own stage charges, `sum_j w_j (Q_j - q_n) / h` (stage 3 the
    point itself).  Summed over a node's devices it is the step's
    collocation derivative of the node's charge, which the run's own end
    derivative (the equation's `-(i + u)`) equals only as well as the stage
    solve converged, divided by the step: on the gold decks up to about 2 %
    of the node's current at fast edges with the cost transform, a few
    tenths of a percent with the dense Newton.  The run says so
    (`transient._RADAU_CURRENTS`).

    `extra[k]` is point k's values, in `charged`' order (a store takes them
    beside the point's rows and drops them)."""

    def __init__(self, charged, epar, x0, w=None):
        self.epar, self.devs, self.slots, self.w = epar, [], [], w
        at = {}
        for j, (_term, c, local, sub) in enumerate(charged):
            if id(c) not in at:
                lims = [c] if getattr(c, 'reset_limit_state', None) is not None else stateful_limiters(c)
                at[id(c)] = len(self.devs)
                self.devs.append({'c': c, 'sub': np.asarray(sub, dtype=np.intp), 'lims': lims,
                                  'u': type(c).u is not Circuit.u, 'qs': None, 'iq': None})
            self.slots.append((at[id(c)], c.get_node_index(c.nodenames[local])))
        ## the saved devices' rows, each device's at its place among them
        self.srows = np.unique(np.concatenate([d['sub'] for d in self.devs]))
        for d in self.devs:
            d['pos'] = np.searchsorted(self.srows, d['sub'])
        self.extra = []
        self._point(x0, 0.0, None)

    def add(self, x, t, step):
        """Point k's values, `step` its accepted step's `(alphas, b)`."""
        if step is None:
            raise RuntimeError('a device current needs its step\'s companion coefficients')
        self._point(x, t, step)

    def _eval(self, d, xe, what):
        """The device's `what` (its methods' names) at its own `xe`, its
        stateful limiters at `xe`."""
        c, epar = d['c'], self.epar
        if d['lims']:
            with devices_at(c, xe, epar, d['lims']):
                return [getattr(c, f)(xe, epar) for f in what]
        return [getattr(c, f)(xe, epar) for f in what]

    def _point(self, x, t, step):
        epar, curs = self.epar, []
        stage = step is not None and step[0] == 'stage'
        if stage:
            _tag, h, Y1, Y2, on_rows = step
            if not on_rows:
                Y1, Y2 = Y1[self.srows], Y2[self.srows]
        for d in self.devs:
            c = d['c']
            q, i = self._eval(d, x[d['sub']], ('q', 'i'))
            if step is None:
                iq, d['qs'] = 0.0 * q, [q, q]
            elif stage:
                w, qn = self.w, d['qs'][0]
                (q1,), (q2,) = self._eval(d, Y1[d['pos']], ('q',)), self._eval(d, Y2[d['pos']], ('q',))
                iq = (w[0] * (q1 - qn) + w[1] * (q2 - qn) + w[2] * (q - qn)) / h
                d['qs'] = [q, qn]
            else:
                alphas, b = step
                iq = alphas[0] * q
                for k in range(1, len(alphas)):
                    iq = iq + alphas[k] * d['qs'][k - 1]
                if b:
                    iq = iq + b * d['iq']
                d['qs'] = [q, d['qs'][0]]
            d['iq'] = iq
            cur = i + iq
            if d['u']:
                cur = cur + c.u(t, epar=epar, analysis='tran')
            curs.append(cur)
        self.extra.append(np.array([curs[dev][k] for dev, k in self.slots], dtype=float))


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

    def __init__(self, circuit, x, sweep_values, rows=None, saved=None, epar=None, currents=None):
        super().__init__(circuit, x, xdot=None, sweep_values=sweep_values, sweep_label='time', sweep_unit='s')
        self.rows = rows
        self.saved = saved
        self.epar = defaultepar if epar is None else epar
        ## the currents of devices holding charge the run made (`DeviceCurrents`)
        self.currents = currents or {}

    def _x(self):
        return self.x if self.rows is None else _Rows(self.x, self.rows)

    def v(self, plus, minus=None):
        if self.rows is None:
            return super().v(plus, minus)
        result = self.circuit.extract_v(self._x(), plus, minus)
        ylabel = f'v({plus},{minus})' if minus is not None else f'v({plus})'
        return self.build_waveform(result, ylabel, 'V')

    def i(self, term):
        if term in self.currents:
            return self.build_waveform(self.currents[term], f'i({term})', 'A')
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
