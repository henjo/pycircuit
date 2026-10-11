"""THE SOURCES IN C: the direct source pass (`_stamp_plan.source_direct`)
-- `u(t)` of a circuit of VS- and IS-family sources -- and the earliest of
their `next_event(t)`, as `csrc/sources.c`.

The Python evaluates per source `v = a + f(t)` (`a` its DC value, `f` its
`function`: `func.TimeFunction`, `Sin`, `Pulse`, `PWL`, `Exp`, `AM`,
`SFFM`), scattered in element order into bins from +0.0.  A `Table` holds
the same sources as the C reads them: one call makes the vector `u(t)`,
another the earliest event after `t` -- bit for bit the Python's, or
declined.

OWNS nothing on others' objects: a `Table` holds what it read.  The C loop
keeps the table it uses on the stepping loop (`c_tab`, `_tran_loop_c`),
checks `holds` once an entry and calls the C on the table's struct itself,
every attempt.  `_paths` prefix `sources_c:`: `sources_c:built` (a table
made), `sources_c:u` and `sources_c:next` (the Python entries' serves), the
declines (`sources_c:<reason>`) and the hand-backs (`sources_c:bail:<why>`);
inside the C loop the pass counts as the Python's (`src.direct:served`,
`u:called`), re-created by `_tran_loop_c._write_back`.  Switch: `ENABLED`
(env `PYCIRCUIT_SOURCES_C=0`), read by `table_for` at every build.

DECLINE, NEVER DIFFER.  A table is built only where every source's value
comes the fast float way -- the waveform one of the seven exactly (no
subclass), the numeric toolkit's own `exp`, `sin`, `where` and `pi`, every
parameter a float (or a small int: its float exact), the DC value and the
waveform not both ints (the Python adds those as ints), an exponential's
`tau` nonzero (the Python raises), numpy's `exp` and `sin` the C library's
on this CPU (`_hdl_cbackend.libm_check`).  A call where the C would leave
the Python's route -- a non-finite value, one past the fast add's 1e300
guard, a sine past `_sin_scalar`'s guards, numpy's error state not
ignoring underflow -- returns None, and the Python makes it.  A table
holds its inputs as it read them and `holds` compares them all (each
source's `function`, its class, its attributes, a PWL's points, the DC
value): stale, it is not used.

`Table.u` and `Table.next_event` are the Python entries (the twin's tests);
the C loop checks what `_call` checks once an entry (`_tran_loop_c.run`).
History: `doc/transient_history.md`.
"""
import os

import numpy as np

from pycircuit.circuit import _cbuild, _paths

_PC = _paths.COUNTS
_no = _paths.no

ENABLED = os.environ.get('PYCIRCUIT_SOURCES_C', '1') != '0'
STATUS = 'not loaded'

SOURCES_C = _cbuild.csrc('sources.c')
SOURCES_CDEF = _cbuild.csrc('sources.h')
_CODES = _cbuild.constants(SOURCES_CDEF)
BAIL = _cbuild.bails(_CODES, 'SOURCES_BAIL_')

_driver = None
#: an int parameter's bound in magnitude (exclusive): below it its float is
#: exact, and so is the difference of two of them
_INT_MAX = 2 ** 52
_F64 = np.float64


def driver():
    """`(ffi, cfn, differ)` once loaded -- `differ` the C library functions
    whose numpy twin gives other bits on this CPU -- or None (`STATUS`)."""
    global _driver, STATUS
    if _driver is None:
        _driver = False
        from pycircuit.circuit import _hdl_cbackend as cb
        try:
            ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(SOURCES_C, SOURCES_CDEF, 'pcx_sources')
        except (_cbuild.CompileError, OSError) as e:
            STATUS = f'off ({e})'
            return None
        _driver = (ffi, cfn, cb.libm_check())
        STATUS = 'on'
    return _driver or None


def _num(v):
    """`v` as the C's double, or None: a float, numpy's, or an int whose
    float is exact (not a bool)."""
    t = type(v)
    if t is float or t is _F64:
        return float(v)
    if t is int and -_INT_MAX < v < _INT_MAX:
        return float(v)
    return None


def _is_int(v):
    return type(v) is int


class Table:
    """One direct source pass's sources as the C reads them, and what they
    were read from (`holds`)."""

    __slots__ = ('arrays', 'cfn', 'entries', 'ffi', 'has_sin', 'kinds', 'n', 'ns', 'out', 'pts',
                 's', 'tk', 'uses_libm')

    def __init__(self, ffi, cfn, n, rows, kinds, a, params, npts, entries, tk):
        self.ffi, self.cfn, self.n, self.ns = ffi, cfn, n, len(kinds)
        poff, off = [], 0
        for p in params:
            poff.append(off)
            off += len(p)
        i64 = np.int64
        self.arrays = arr = {
            'kind': np.array(kinds, dtype=i64),
            'isv': np.array([r[0] for r in rows], dtype=i64),
            'row0': np.array([r[1] for r in rows], dtype=i64),
            'row1': np.array([r[2] for r in rows], dtype=i64),
            'row2': np.array([r[3] for r in rows], dtype=i64),
            'poff': np.array(poff, dtype=i64),
            'npts': np.array(npts, dtype=i64),
            'a': np.array(a, dtype=float),
            'par': np.array([x for p in params for x in p] or [0.0], dtype=float),
        }
        self.out = np.empty(max(n, 1))
        self.s = s = ffi.new('sources_t *')
        s.n, s.ns = n, self.ns
        fb = ffi.from_buffer
        self.pts = [fb('int64_t *' if v.dtype == i64 else 'double *', v) for v in arr.values()]
        for (k, _v), h in zip(arr.items(), self.pts):
            setattr(s, k, h)
        self.pts.append(fb('double *', self.out))
        s.out = self.pts[-1]
        self.kinds = frozenset(kinds)
        self.uses_libm = bool(self.kinds & {_CODES['SOURCES_SIN'], _CODES['SOURCES_EXP'],
                                            _CODES['SOURCES_AM'], _CODES['SOURCES_SFFM']})
        self.has_sin = _CODES['SOURCES_SIN'] in self.kinds
        self.entries, self.tk = entries, tk

    def holds(self):
        """Whether every input the table was built from is as it was read:
        each source's `function` (the object, its class, its attributes,
        a PWL's points) and its DC value."""
        for el, key, fn, cls, attrs, pts, a in self.entries:
            d = el.__dict__
            if d.get('function', fn) is not fn or type(fn) is not cls:
                return False
            if d['iparv']._values.get(key) is not a:
                return False
            fd = fn.__dict__
            if len(fd) != len(attrs) or any(fd.get(k) is not v for k, v in attrs):
                return False
            if pts is not None and (tuple(fn.times), tuple(fn.values)) != pts:
                return False
        return True

    def _call(self, mode, t):
        ## (what the Python reads on every call: the toolkit's functions, the
        ## sine's scalar switch, numpy's error state for underflow)
        tk = self.tk
        if not (tk.exp is np.exp and tk.sin is np.sin and tk.where is np.where):
            return _no('sources_c:toolkit')
        if self.uses_libm:
            from pycircuit.circuit import func
            if (self.has_sin and not func.SIN_SCALAR) or not func._underflow_ignored():
                return _no('sources_c:errstate')
        if not (type(t) is float or type(t) is _F64):
            return _no('sources_c:t')
        s = self.s
        s.mode, s.t = mode, float(t)
        status = self.cfn(s)
        if status != _CODES['SOURCES_DONE']:
            return _no('sources_c:bail:' + BAIL.get(status, str(status)))
        return s

    def u(self, t):
        """`u(t)` as the direct pass makes it (a new array), or None."""
        s = self._call(_CODES['SOURCES_U'], t)
        if s is None:
            return None
        _PC['sources_c:u'] += 1
        return self.out[:self.n].copy()

    def next_event(self, t):
        """The earliest of the sources' `next_event(t)`, or None."""
        s = self._call(_CODES['SOURCES_NEXT'], t)
        if s is None:
            return None
        _PC['sources_c:next'] += 1
        return float(s.t_next)


def _waveform(fn, func, tk, differ, a):
    """`(kind, params, npts, pts)` of a source's `function` as the C reads
    it, or a decline reason (str)."""
    C = _CODES
    cls = type(fn)
    if cls is func.TimeFunction:
        ## (`f` is the int 0: an int DC value adds as ints in the Python)
        return 'int' if _is_int(a) else (C['SOURCES_DC'], (), 0, None)
    if cls is func.Sin:
        if not func.SIN_SCALAR or {'exp', 'sin'} & differ:
            return 'libm'
        names = ('td', 'theta', 'omega', 'phase', 'offset', 'amplitude')
        kind = C['SOURCES_SIN']
        ## (`_sin_scalar`'s guard on its parameters: past it, numpy's route)
        if not all(type(getattr(fn, k)) in func._NUMBER_T
                   and -func._BIG < getattr(fn, k) < func._BIG for k in names):
            return 'param'
    elif cls is func.Pulse:
        ## (a zero edge -- `__init__` floors both at `MIN_EDGE`, an assignment
        ## after it does not: the Python's eager `where` divides by it and
        ## raises)
        if not (fn.tr and fn.tf):
            return 'edge'
        names = ('v1', 'v2', 'td', 'tr', 'tf', 'pw', 'per')
        kind = C['SOURCES_PULSE']
    elif cls is func.Exp:
        if 'exp' in differ:
            return 'libm'
        if not (fn.tau1 and fn.tau2):
            return 'tau'
        ## (`f` is `v1` itself before `td1`: an int added as ints)
        if _is_int(a) and _is_int(fn.v1):
            return 'int'
        names = ('v1', 'v2', 'td1', 'tau1', 'td2', 'tau2')
        kind = C['SOURCES_EXP']
    elif cls is func.AM or cls is func.SFFM:
        if 'sin' in differ:
            return 'libm'
        names = (('vo', 'va', 'fc', 'fm', 'm') if cls is func.AM
                 else ('vo', 'va', 'fc', 'mdi', 'fm'))
        kind = C['SOURCES_AM'] if cls is func.AM else C['SOURCES_SFFM']
        p = [_num(getattr(fn, k)) for k in names] + [_num(tk.pi)]
        return 'param' if None in p else (kind, tuple(p), 0, None)
    elif cls is func.PWL:
        T, V = tuple(fn.times), tuple(fn.values)
        if not T or len(T) != len(V):
            return 'pwl'
        p = [_num(x) for x in T + V]
        if None in p or (_is_int(a) and any(map(_is_int, V))):
            return 'param' if None in p else 'int'
        return C['SOURCES_PWL'], tuple(p), len(T), (T, V)
    else:
        return 'class'
    p = [_num(getattr(fn, k)) for k in names]
    return 'param' if None in p else (kind, tuple(p), 0, None)


def table_for(sd):
    """The `Table` of a direct source pass `sd` (`_stamp_plan._SourceDirect`:
    its sources, rows and order), or None (declined, counted)."""
    if not ENABLED:
        return _no('sources_c:off')
    drv = driver()
    if drv is None:
        return _no('sources_c:driver')
    ffi, cfn, differ = drv
    from pycircuit.circuit import func
    from pycircuit.circuit.toolkit import NumericToolkit
    n = sd.n
    if not sd.ent:
        return _no('sources_c:empty')
    rows, kinds, avals, params, npts, entries = [], [], [], [], [], []
    tk0 = None
    for el, key, _desc, r in sd.ent:
        d = el.__dict__
        fn = d.get('function')
        if fn is None:
            fn = el.function
        tk = fn.__dict__.get('toolkit')
        if type(tk) is not NumericToolkit or tk.exp is not np.exp or tk.sin is not np.sin \
                or tk.where is not np.where:
            return _no('sources_c:toolkit')
        if tk0 is None:
            tk0 = tk
        a = d['iparv']._values.get(key)
        af = _num(a)
        if af is None:
            return _no('sources_c:param')
        w = _waveform(fn, func, tk, differ, a)
        if w.__class__ is str:
            return _no('sources_c:' + w)
        kind, p, npt, pts = w
        if key == 'v':
            if not 0 <= r[2] < n:
                return _no('sources_c:rows')
            rows.append((1, 0, 0, r[2]))
        else:
            if not (0 <= r[0] < n and 0 <= r[1] < n):
                return _no('sources_c:rows')
            rows.append((0, r[0], r[1], 0))
        kinds.append(kind)
        avals.append(af)
        params.append(p)
        npts.append(npt)
        entries.append((el, key, fn, type(fn), tuple(fn.__dict__.items()), pts, a))
    _PC['sources_c:built'] += 1
    return Table(ffi, cfn, n, rows, kinds, avals, params, npts, entries, tk0)
