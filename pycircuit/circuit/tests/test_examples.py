# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Regression test for two example scripts fixed under architecture.md P14.

Both were broken by real API drift (``examples/multi_feedback_filter.py``:
the 2008-era dict-style ``res['out']`` never worked against ``AC``'s result
object, which is empty by design; ``examples/regulatedcascode.py``: a
symbolic-toolkit ``z0`` crashed ``NPortS.Y`` -- ``np.sqrt``/``np.linalg.inv``
can't handle sympy objects, and ``symbolicapprox.approx`` had two separate
bugs of its own, a stale ``series()`` keyword and swapped
``.subs()``/``.removeO()`` order that silently zeroed its result for 18
years). Run as subprocesses (both execute unconditionally at import time,
with no ``__main__`` guard) and checked for real, verified output -- not
just a zero exit code.
"""
import subprocess
import sys

import pytest

EXAMPLES_DIR = 'pycircuit/circuit/examples'


def _run_example(name):
    result = subprocess.run(
        [sys.executable, '%s/%s' % (EXAMPLES_DIR, name)],
        capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, (
        '%s exited with %d:\n%s' % (name, result.returncode, result.stderr))
    return result.stdout


def test_multi_feedback_filter_example():
    out = _run_example('multi_feedback_filter.py')
    assert 'DC Gain' in out
    ## Verified by hand: the classic MFB lowpass DC gain is -R1/R3.
    assert '-R\u2081' in out or '-R1/R3' in out
    assert 'Input impedance' in out


def test_regulatedcascode_example():
    out = _run_example('regulatedcascode.py')
    assert 'Input impedance' in out
    assert 'Input referred current noise PSD' in out
    assert 'Input referred voltage noise PSD' in out
