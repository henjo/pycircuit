# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Tests for the symbolic-backend benchmark fixtures.

These guard two things that the whole DDD comparison rests on:

1. the fixtures really are built from pycircuit circuits through the normal MNA
   path (not hand-written matrices), and they solve to the right answer;
2. the measurement code reproduces the operation counts already published for
   the sequence-of-expressions solver.  A harness that cannot reproduce our own
   numbers cannot referee anyone else's.
"""

import os
import sys

import numpy as np
import pytest
import sympy

from pycircuit.circuit import benchmark_circuits as bc

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                os.pardir, os.pardir, os.pardir, 'benchmarks'))


def _solve_numerically(system, freq=1e6):
    """Substitute every symbol and solve with numpy -- the reference answer."""
    env = dict(system.params)
    env[system.s] = 1j * 2 * np.pi * freq
    n = system.A.rows
    A = np.array([[complex(system.A[i, j].subs(env)) for j in range(n)]
                  for i in range(n)], dtype=complex)
    b = np.array([complex(system.b[i].subs(env)) for i in range(n)], dtype=complex)
    return np.linalg.solve(A, b)


@pytest.mark.parametrize('N', [3, 4, 6])
def test_rc_ladder_is_built_from_a_circuit(N):
    """The ladder comes from a real SubCircuit, not a hand-written matrix."""
    system = bc.rc_ladder(N)
    assert system.cir is not None
    ## MNA carries the voltage-source branch current, so the system is one
    ## larger than the node count.
    assert system.dim == N + 1
    assert system.s in system.A.free_symbols


@pytest.mark.parametrize('N', [3, 4, 6])
def test_rc_ladder_is_a_lowpass(N):
    """Sanity-check the physics: |H| falls with frequency and is <= 1."""
    lo = _solve_numerically(bc.rc_ladder(N), freq=1.0)
    hi = _solve_numerically(bc.rc_ladder(N), freq=1e9)
    system = bc.rc_ladder(N)
    glo = abs(lo[system.out_index] / lo[system.in_index])
    ghi = abs(hi[system.out_index] / hi[system.in_index])
    assert glo <= 1.0 + 1e-9
    assert ghi < glo


def test_mfb_filter_builds_and_solves():
    system = bc.mfb_filter()
    assert system.cir is not None
    x = _solve_numerically(system)
    assert np.isfinite(x).all()


@pytest.mark.parametrize('N', [4, 6])
def test_semi_symbolic_keeps_only_requested_symbols(N):
    system = bc.rc_ladder_semi_symbolic(N, n_symbolic=2)
    free = system.A.free_symbols - {system.s}
    ## Two resistors symbolic; every capacitor and the rest numeric.
    assert len(free) == 2
    assert all(str(sym).startswith('R') for sym in free)


@pytest.mark.parametrize('n', [3, 4, 5])
def test_dense_matrix_is_full_and_distinct(n):
    system = bc.dense_symbolic_matrix(n)
    assert system.dim == n
    assert len(system.A.free_symbols) == n * n      # every entry its own symbol
    assert system.cir is None                       # deliberately not a circuit


@pytest.mark.parametrize('N', [3, 4])
def test_legacy_ladder_matches_the_prototype_matrix(N):
    """The calibration fixture must stay identical to the published one."""
    system = bc.legacy_soe_ladder(N)
    s = system.s
    R = [sympy.Symbol('R%d' % i, positive=True) for i in range(N)]
    C = [sympy.Symbol('C%d' % i, positive=True) for i in range(N)]
    for i in range(N):
        expected = 1 / R[i] + s * C[i] + (1 / R[i - 1] if i > 0 else 0)
        assert sympy.simplify(system.A[i, i] - expected) == 0


def test_standard_suite_is_stable():
    """Guard the fixed comparison set -- changing it invalidates published runs."""
    names = [sysm.name for _, sysm in bc.standard_suite()]
    assert names == [
        'rc_ladder_N4', 'rc_ladder_N8', 'rc_ladder_N12', 'rc_ladder_N16',
        'rc_ladder_semi_N4_k2', 'rc_ladder_semi_N8_k2',
        'rc_ladder_semi_N12_k2', 'rc_ladder_semi_N16_k2',
        'mfb_filter',
        'dense_3', 'dense_4', 'dense_5', 'dense_6',
    ]


## -- the harness itself ---------------------------------------------------

@pytest.mark.parametrize('N,expected', [(4, 73), (8, 157), (12, 241), (16, 325)])
def test_harness_reproduces_published_soe_op_counts(N, expected):
    """Stage B's gate: reproduce the numbers in ``soe_symbolic.rst``."""
    from symbolic_bench import measure
    rec = measure('legacy_soe:%d' % N, 'soe')
    assert rec['status'] == 'ok', rec.get('error')
    assert rec['size'] == expected
    assert rec['correct'] is True


@pytest.mark.parametrize('backend', ['symbolic_poly', 'soe'])
def test_backends_agree_with_a_numeric_solve(backend):
    from symbolic_bench import measure
    rec = measure('rc_ladder:4', backend)
    assert rec['status'] == 'ok', rec.get('error')
    assert rec['correct'] is True


def test_harness_reports_rather_than_raises_on_a_bad_backend():
    """A failing backend must be recorded, not blow up the run."""
    from symbolic_bench import measure, BACKENDS

    def boom(system):
        raise RuntimeError('deliberate')

    BACKENDS['_boom'] = boom
    try:
        rec = measure('rc_ladder:3', '_boom')
    finally:
        del BACKENDS['_boom']
    assert rec['status'] == 'error'
    assert 'deliberate' in rec['error']
