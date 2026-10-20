"""The allocator's mode pinned on request (`_malloc`, 2026-10-21).

Each case is a fresh process: the thresholds are a process's, set when
`pycircuit.circuit` is imported with `PYCIRCUIT_MALLOC_PIN=1` or by `pin()`.
What pinning does is read from glibc's own count of mmap'ed chunks
(`mallinfo2().hblks`) across a 24 MB allocation: pinned it comes from the
heap, under a threshold of 128 KB from mmap.
"""
import ctypes
import os
import subprocess
import sys

import pytest

from pycircuit.circuit import _malloc

pytestmark = pytest.mark.skipif(
    _malloc._glibc() is None or not hasattr(ctypes.CDLL(None), 'mallinfo2'),
    reason='glibc >= 2.33 only')

_PROBE = '''
import ctypes, sys
from pycircuit.circuit import _malloc
import numpy as np
class MI2(ctypes.Structure):
    _fields_ = [(n, ctypes.c_size_t) for n in (
        'arena', 'ordblks', 'smblks', 'hblks', 'hblkhd', 'usmblks', 'fsmblks',
        'uordblks', 'fordblks', 'keepcost')]
libc = ctypes.CDLL(None)
libc.mallinfo2.restype = MI2
if sys.argv[1:] == ['pin']:
    _malloc.pin()
h0 = libc.mallinfo2().hblks
a = np.empty(3 << 20)
print(_malloc.STATUS, libc.mallinfo2().hblks - h0, sep='|')
'''

_OURS = ('PYCIRCUIT_MALLOC_PIN', 'GLIBC_TUNABLES') + _malloc._USER_ENV


def _probe(*args, **env):
    base = {k: v for k, v in os.environ.items() if k not in _OURS}
    out = subprocess.run([sys.executable, '-c', _PROBE, *args],
                         capture_output=True, text=True,
                         env=dict(base, **env), timeout=600, check=False)
    assert out.returncode == 0, out.stderr[-2000:]
    status, mapped = out.stdout.strip().split('|')
    return status, int(mapped)


def test_glibcs_own_mode_by_default():
    assert _probe()[0] == 'not pinned'
    assert _probe(PYCIRCUIT_MALLOC_PIN='0')[0] == 'not pinned'


def test_a_24_mb_array_lives_in_the_heap_once_pinned():
    assert _probe(PYCIRCUIT_MALLOC_PIN='1') == ('pinned', 0)
    assert _probe('pin') == ('pinned', 0)


def test_a_threshold_the_user_set_wins():
    assert _probe(PYCIRCUIT_MALLOC_PIN='1', MALLOC_MMAP_THRESHOLD_='131072') \
        == ('user', 1)
    assert _probe('pin', GLIBC_TUNABLES='glibc.malloc.mmap_threshold=131072') \
        == ('user', 1)
