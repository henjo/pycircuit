"""EXACT INVALIDATION: one counter that moves on any change to a watched
dict, so that a C path's readiness check that passed holds, unchecked,
while the counter stands.  No Python reference: without it every checker
(`_tran_core`, `_tran_newton_c`, `_tran_lte_c`, `_tran_loop_c`,
`_tran_radau`, `_tran_radau_tc`, `_stamp_plan`) runs its check in full on
every call.

OWNS.  Nothing on others' objects but CPython's watcher mark on each dict a
checker arms (one dict watcher id, registered by `load`, cleared at exit by
`_clear`); the stamps are the checkers' own.  No `_paths` counts.  Switch:
`PYCIRCUIT_WATCH=0` (`ENABLED`, read at import; `load` consults it once).

THE COUNTER.  CPython (3.12+) lets a C callback watch a dict: here one
callback (`csrc/watch.c`) bumps ONE counter on any change to a watched dict
(an item set, added or deleted, the dict cleared, cloned or freed).  The
checks it serves -- is every batch's kernel still its class's bound one,
every element without an instance shadow, every pack the one mirrored, the
limiter walk's tables as built, the plan the circuit's -- read mostly
dicts, and their answer does not change between steps.  A checker that
passed stamps the counter (`arm`); while the counter stands, no dict it
read has changed, and the check passes again.

`arm(dicts, before)` watches every dict the checker read (idempotent) and
returns the stamp: the counter's value, or -1 where a stamp would not be
exact.  The counter read BEFORE the full check (`before`) must still stand
after the arming (a change during the check -- the checker's own write
into a watched dict, a thread -- leaves the next call to check again).
What no dict watcher sees each checker compares on every call, cheaply: an
instance's `__dict__` REASSIGNED (the objects' identities), a class
attribute (a pass's code, `ParameterDict`'s epoch, a hand-written
limiter's `vlimit`), a list's length.  (CPython's TYPE watchers fire only
while the type holds a version tag, which no public API can confirm: none
are used.)

`EPOCH` is the counter (a `ctypes.c_int64` over the C variable), or None
where watching is off -- `ENABLED` false, the API missing (before 3.12),
the build failing or the watcher cleared at exit (`STATUS` says): every
check then runs in full.  History: `doc/transient_history.md`,
`_watch.py`.
"""
import atexit
import ctypes
import os

from pycircuit.circuit import _cbuild

ENABLED = os.environ.get('PYCIRCUIT_WATCH', '1') != '0'
STATUS = 'not loaded'
EPOCH = None
#: A CHECKER STOPS ARMING after this many stamps: a circuit whose watched
#: dicts change every step (a stateful limiter's state, written into its
#: element's dict each iteration) would otherwise pay the full check AND
#: the arming on every call -- past the cap it pays the full check alone
#: (its stamp -1)
MAX_ARMS = 16
#: ... COUNTING ONLY THE STAMPS THAT CHURN: a stamp that held `HOLD` checks
#: or more before it broke watched dicts that change between runs, not
#: within them -- a parameter set between the solves of a sweep, a Monte
#: Carlo draw -- and its re-arm starts the count again (`counted`), so a
#: circuit object solved many times keeps its stamps
HOLD = 8

WATCH_C = _cbuild.csrc('watch.c')
WATCH_CDEF = _cbuild.csrc('watch.h')

_API = {}


def load():
    """Build and register the counter once: `EPOCH`, or None (`STATUS`)."""
    global EPOCH, STATUS
    if STATUS != 'not loaded':
        return EPOCH
    STATUS = 'off'
    if not ENABLED:
        STATUS = 'off (PYCIRCUIT_WATCH=0)'
        return None
    api = ctypes.pythonapi
    try:
        add_d, watch_d = api.PyDict_AddWatcher, api.PyDict_Watch
    except AttributeError:
        STATUS = 'off (no dict watchers: CPython 3.12+)'
        return None
    try:
        ffi, cfn, _key, _cold, _secs = _cbuild.load_unit(WATCH_C, WATCH_CDEF)
    except (_cbuild.CompileError, OSError) as e:
        STATUS = f'off ({e})'
        return None
    add_d.restype, add_d.argtypes = ctypes.c_int, [ctypes.c_void_p]
    watch_d.restype, watch_d.argtypes = ctypes.c_int, [ctypes.c_int, ctypes.py_object]
    did = add_d(int(ffi.cast('uintptr_t', cfn(1))))
    if did < 0:
        STATUS = 'off (no dict watcher id left)'
        return None
    _API.update(ffi=ffi, cfn=cfn, did=did, watch_d=watch_d)
    EPOCH = ctypes.c_int64.from_address(int(ffi.cast('uintptr_t', cfn(0))))
    STATUS = 'on'
    ## (the watcher removed at exit, before the modules are torn down: a
    ## watched dict freed then must not call into an unloaded object)
    atexit.register(_clear)
    return EPOCH


def _clear():
    """The watcher removed (registered at exit by `load`): `EPOCH` None,
    every check in full from then."""
    global EPOCH, STATUS
    did = _API.get('did')
    if did is not None:
        clear = ctypes.pythonapi.PyDict_ClearWatcher
        clear.restype, clear.argtypes = ctypes.c_int, [ctypes.c_int]
        clear(did)
        _API.pop('did')
    EPOCH = None
    STATUS = 'off (cleared at exit)'


def counted(arms, held):
    """A checker's arm count after one more arm: `arms + 1`, or 1 where
    its stamp before held `HOLD` checks or more (`HOLD`)."""
    return 1 if held >= HOLD else arms + 1


def now():
    """The counter, or -1 (watching off: -1 is never a stamp)."""
    ep = EPOCH if STATUS == 'on' else load()
    return -1 if ep is None else ep.value


def arm(dicts, before):
    """Watch every dict of `dicts` (what a check that passed read); the
    stamp for it: the counter, where it still stands at `before` (`now()`
    before the check), else -1."""
    if before < 0 or EPOCH is None:
        return -1
    did, watch_d = _API['did'], _API['watch_d']
    for d in dicts:
        if not isinstance(d, dict) or watch_d(did, d) != 0:
            return -1
    ep = EPOCH.value
    return ep if ep == before else -1
