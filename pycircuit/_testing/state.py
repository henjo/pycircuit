"""The hdl backend state of a process, read without changing it (2026-10-03).

Every behavioural class carrying its own `_hdl_info` -- library classes,
collapse variants, fold classes, classes a test defined -- found by walking
`hdl.Behavioural`'s subclasses, grouped by the `info` dict they share (a
plain subclass inherits its parent's `info`, so one `info` can be read
through several classes; its OWNER is the class with `_hdl_info` in its own
`__dict__`).

`hdl_state()` is what a gate dumps per worker at the end of a session
(`PYCIRCUIT_STATE_DUMP`), so two workers can be diffed after an anomaly;
`invariant_breaks()` is the consistency check the leak detector runs: a class
marked bound to the C backend must carry its kernels, and one marked unbound
must carry none.  The anomaly that motivated it (gates G97/G98, 2026-10-03):
a MosLevel1 collapse variant with `_c_bound` True and no `_hdl_c` on its
chain functions, so numpy ran in one worker of eight.

Nothing is imported that the process has not imported already: a session
that never touched the hdl compiler reports nothing.
"""
import sys

#: the chain functions that carry a C kernel (`_hdl_cbackend.C_FUNCS`)
C_FUNCS = ('i', 'G', 'q', 'C', 'i_dc', 'G_dc')


def hdl_classes():
    """Every live class with `_hdl_info` in its own `__dict__`, or [] when
    the hdl module was never imported."""
    hdl = sys.modules.get('pycircuit.circuit.hdl')
    if hdl is None:
        return []
    base = getattr(hdl, 'Behavioural', None)
    if base is None:
        return []
    seen, stack, out = set(), [base], []
    while stack:
        c = stack.pop()
        for s in type.__subclasses__(c):
            if s in seen:
                continue
            seen.add(s)
            stack.append(s)
            if '_hdl_info' in s.__dict__:
                out.append(s)
    return out


def _kernel_key(fn):
    k = fn.__dict__.get('_hdl_c') if fn is not None else None
    return None if k is None else getattr(k, 'key', repr(k))


def info_state(cls):
    """One owner class's backend state, as plain data: comparable within a
    process and across workers (no object identities)."""
    info = cls.__dict__['_hdl_info']
    funcs = info.get('funcs') or {}
    out = {
        'class': f'{cls.__module__}.{cls.__qualname__}',
        'status': cls.__dict__.get('_hdl_backend_status'),
        'pin': cls.__dict__.get('hdl_backend'),
        'chained': bool(info.get('chained')),
        'c_bound': info.get('_c_bound'),
        'pending': info.get('_backend_pending'),
        'seen': info.get('_backend_seen'),
        'limit_status': info.get('_c_limit_status'),
        'funcs': {},
    }
    for name in C_FUNCS:
        fn = funcs.get(name)
        if fn is None:
            continue
        out['funcs'][name] = {
            'csrc': getattr(fn, '_csrc', None) is not None,
            'kernel': _kernel_key(fn),
            'optimised': '_hdl_ref' in fn.__dict__,
            'qualname': getattr(fn, '__qualname__', None),
        }
    return out


def hdl_state():
    """`{owner class name: info_state}` for every live hdl class, plus the
    statuses of the subclasses that read an inherited `info`."""
    out = {}
    owners = hdl_classes()
    by_info = {id(c.__dict__['_hdl_info']): c for c in owners}
    for c in owners:
        out[f'{c.__module__}.{c.__qualname__}'] = info_state(c)
    hdl = sys.modules.get('pycircuit.circuit.hdl')
    base = getattr(hdl, 'Behavioural', None) if hdl is not None else None
    if base is not None:
        seen, stack = set(), [base]
        while stack:
            c = stack.pop()
            for s in type.__subclasses__(c):
                if s in seen:
                    continue
                seen.add(s)
                stack.append(s)
                if '_hdl_info' not in s.__dict__ and '_hdl_backend_status' in s.__dict__:
                    owner = by_info.get(id(getattr(s, '_hdl_info', None)))
                    if owner is not None:
                        key = f'{owner.__module__}.{owner.__qualname__}'
                        out[key].setdefault('readers', {})[
                            f'{s.__module__}.{s.__qualname__}'] = s.__dict__['_hdl_backend_status']
    return out


def invariant_breaks(state=None):
    """The owners whose C binding flag disagrees with their kernels:
    `{class name: reason}`.  Bound means at least one chain function
    carries C source and every such function carries its kernel; unbound
    means no chain function carries a kernel.  The limiter is not read (it
    may refuse itself on a bound class)."""
    state = hdl_state() if state is None else state
    bad = {}
    for name, st in state.items():
        fs = st['funcs']
        with_src = [f for f, d in fs.items() if d['csrc']]
        with_kern = [f for f, d in fs.items() if d['kernel'] is not None]
        if st['c_bound']:
            missing = [f for f in with_src if fs[f]['kernel'] is None]
            if not with_src:
                bad[name] = 'bound to C, but no chain function carries C source'
            elif missing:
                bad[name] = f'bound to C, but {missing} carry no kernel'
        elif with_kern:
            bad[name] = f'not bound to C, but {with_kern} carry a kernel'
    return bad
