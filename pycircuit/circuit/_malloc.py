"""THE ALLOCATOR'S MODE, pinned on request (2026-10-21; glibc only).

glibc serves an allocation above `M_MMAP_THRESHOLD` (128 KB at start: a
dense n x n float64 array from n ~ 130) with fresh mmap'ed pages, and the
first time such a chunk is freed raises the threshold to its size, and
`M_TRIM_THRESHOLD` (when the heap's free top is handed back to the OS) to
twice that.  A process therefore settles into a MODE by its allocation
history: the same deck ran 26.6 or 28.4 s (arom, by whether the coupled
stage C had run) and 31.6 .. 43 s (opampal), and a change that allocates
differently can flip a deck's mode and read as a cost it does not have
(roadmap 156's +5.4 % on arom).

Pinned -- both thresholds set once, which also switches glibc's adjustment
off -- a large array lives in the heap and a freed one is reused, not
unmapped and faulted in again: `MMAP_THRESHOLD` 32 MB (glibc's own ceiling
for it on 64 bits), `TRIM_THRESHOLD` 256 MB.  Every answer is the same
bytes; only where the memory comes from changes.

NOT the default: measured on the 27 gold decks' default transient (three
rounds interleaved, best of 3) pinning moved the summed time +0.0 % and
Radau's on its 20 decks -1.0 %, against a 3 % line for pinning by default;
a pinned process keeps more at its end (mux8 459 -> 508 MB).  So the
timing harnesses pin it (`benchmarks/_bench.pin_malloc`, the gold census),
where a mode flip would read as a change's cost.

Env `PYCIRCUIT_MALLOC_PIN=1` pins when `pycircuit.circuit` is imported, or
call `pin()`.  A threshold the user set (`MALLOC_MMAP_THRESHOLD_`,
`MALLOC_TRIM_THRESHOLD_`, `MALLOC_TOP_PAD_`, `MALLOC_MMAP_MAX_` or a
`glibc.malloc.` tunable) wins: the process's allocator is then theirs, and
nothing is set.  `STATUS` says which.
"""
import ctypes
import os
import sys

#: mallopt's parameters (malloc.h)
M_TRIM_THRESHOLD = -1
M_MMAP_THRESHOLD = -3

MMAP_THRESHOLD = 32 * 1024 * 1024
TRIM_THRESHOLD = 256 * 1024 * 1024

_USER_ENV = ('MALLOC_MMAP_THRESHOLD_', 'MALLOC_TRIM_THRESHOLD_',
             'MALLOC_TOP_PAD_', 'MALLOC_MMAP_MAX_')

ENABLED = os.environ.get('PYCIRCUIT_MALLOC_PIN', '0') != '0'
#: 'not pinned', or what `pin()` found: 'pinned', 'user' (a threshold the
#: user set), 'not glibc', 'refused' (mallopt said no)
STATUS = 'not pinned'


def _glibc():
    """The C library when it is glibc, else None (musl's mallopt is a stub)."""
    if not sys.platform.startswith('linux'):
        return None
    try:
        libc = ctypes.CDLL(None)
    except OSError:
        return None
    if not hasattr(libc, 'gnu_get_libc_version') or not hasattr(libc, 'mallopt'):
        return None
    libc.mallopt.argtypes = (ctypes.c_int, ctypes.c_int)
    libc.mallopt.restype = ctypes.c_int
    return libc


def _user_set():
    return any(k in os.environ for k in _USER_ENV) \
        or 'glibc.malloc.' in os.environ.get('GLIBC_TUNABLES', '')


def pin():
    """Pin both thresholds, once a process; returns (and sets) `STATUS`."""
    global STATUS
    if STATUS != 'not pinned':
        return STATUS
    if _user_set():
        STATUS = 'user'
    else:
        libc = _glibc()
        if libc is None:
            STATUS = 'not glibc'
        elif libc.mallopt(M_MMAP_THRESHOLD, MMAP_THRESHOLD) \
                and libc.mallopt(M_TRIM_THRESHOLD, TRIM_THRESHOLD):
            STATUS = 'pinned'
        else:
            STATUS = 'refused'
    return STATUS


if ENABLED:
    pin()
