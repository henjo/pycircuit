"""The sanitizer build of all C (`PYCIRCUIT_C_SANITIZE=1`; testing for
development, stage 4, 2026-10-04).

Every C object the package builds can be compiled with AddressSanitizer
and UndefinedBehaviorSanitizer (`_hdl_cbackend.SANITIZE`), and
`scripts/sanitize_suite.sh` runs the suite that way.  Here, in the ordinary
suite: a sanitized process stops on a write one past a kernel's output and
names it (the instrument is alive), runs the same kernel in bounds, keys
its objects apart from the ordinary ones, and refuses to start without the
runtime preloaded (it would otherwise fall back to numpy everywhere and
test nothing).
"""
import os
import subprocess
import sys

import pytest

from pycircuit.circuit import _hdl_cbackend as cb
from pycircuit.circuit import hdl


def _gcc_lib(name):
    try:
        p = subprocess.run(['gcc', f'-print-file-name={name}'], capture_output=True,
                           text=True, timeout=60, check=False).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return p if os.path.isabs(p) and os.path.exists(p) else None


#: the runtime, and libstdc++ beside it: Python is not C++, and without
#: libstdc++ loaded when ASan starts its `__cxa_throw` interceptor has no
#: real function to call -- jaxlib throws at import and ASan aborts
ASAN = _gcc_lib('libasan.so')
STDCXX = _gcc_lib('libstdc++.so')
needs_asan = pytest.mark.skipif(ASAN is None or STDCXX is None or cb.find_compiler()[0] is None,
                                reason='no gcc AddressSanitizer runtime or no compiler here')


def _python(code, sanitize=True, preload=True):
    env = {k: v for k, v in os.environ.items()
           if k not in ('LD_PRELOAD', 'PYCIRCUIT_C_SANITIZE', 'ASAN_OPTIONS')}
    if sanitize:
        env['PYCIRCUIT_C_SANITIZE'] = '1'
    if preload:
        env['LD_PRELOAD'] = f'{ASAN} {STDCXX}'
        env['ASAN_OPTIONS'] = 'detect_leaks=0:halt_on_error=1'
        env['PYTHONMALLOC'] = 'malloc'
    return subprocess.run([sys.executable, '-c', code], env=env, capture_output=True,
                          text=True, timeout=600, check=False)


#: a chain-function kernel writing `out[0..__N__]` into a 256-double output
#: (2 KiB: past numpy's small-block cache, a block of its own with a redzone)
PLANT = r'''
import numpy as np
from pycircuit.circuit import _hdl_cbackend as cb, hdl
src = ("void %s(const double *x, const double *p, double *out) "
       "{ for (int k = 0; k <= __N__; k++) out[k] = x[0] + k; }" % hdl._C_ENTRY)
ffi, cfn, key, cold, secs = cb.load_kernel(src)
x = np.ones(1)
out = np.zeros(256)
cfn(ffi.from_buffer('double *', x), ffi.NULL, ffi.from_buffer('double *', out))
print('RETURNED', out[0], out[255], key)
'''


@needs_asan
def test_a_sanitized_process_stops_on_a_write_past_the_output():
    bad = _python(PLANT.replace('__N__', '256'))
    assert bad.returncode != 0, bad.stdout
    assert 'heap-buffer-overflow' in bad.stderr and 'RETURNED' not in bad.stdout
    assert hdl._C_ENTRY in bad.stderr           # the frame is named
    ok = _python(PLANT.replace('__N__', '255'))
    assert ok.returncode == 0, ok.stderr[-2000:]
    assert 'RETURNED 1.0 256.0' in ok.stdout


@needs_asan
def test_the_sanitized_objects_are_keyed_apart():
    code = ("from pycircuit.circuit import _hdl_cbackend as cb; "
            "print(cb.SANITIZE, cb.CFLAGS[-len(cb.SANITIZE_FLAGS):] == cb.SANITIZE_FLAGS, "
            "cb.source_key('int f;'))")
    on = _python(code)
    off = _python(code, sanitize=False, preload=False)
    assert on.returncode == 0 and off.returncode == 0, on.stderr[-2000:] + off.stderr[-2000:]
    s_on, flags_on, key_on = on.stdout.split()
    s_off, _flags, key_off = off.stdout.split()
    assert (s_on, flags_on, s_off) == ('True', 'True', 'False')
    assert key_on != key_off


def test_a_sanitized_process_without_the_runtime_refuses_to_start():
    r = _python('import pycircuit.circuit._hdl_cbackend', preload=False)
    assert r.returncode != 0
    assert 'needs the AddressSanitizer runtime preloaded' in r.stderr
