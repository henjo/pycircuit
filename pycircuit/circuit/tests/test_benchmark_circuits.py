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


def test_leapfrog_is_a_fifth_order_lowpass():
    """The fixture must be the filter it claims to be.

    A leapfrog realises an LC ladder as integrators with feedback between
    adjacent stages, and its terminating integrators must be **damped** to
    simulate the ladder's source and load resistances.  The first version of
    this fixture omitted that: it integrated at DC, never developed a
    passband, and rolled off at roughly half the intended rate -- a circuit
    that would have been a perfectly plausible-looking test case for
    something that was not a fifth-order filter.

    So this checks the two things that would catch that: a flat passband at
    the doubly-terminated ladder's 0.5, and a stopband slope of 100 dB per
    decade, which is 5 x 20.
    """
    import numpy as np
    import sympy
    from pycircuit.circuit import benchmark_circuits as bc

    system = bc.leapfrog_5th_order()
    n = system.dim
    rows = sympy.Matrix(system.A)
    source = sympy.Matrix(system.b)
    vin = list(system.b.free_symbols)[0]

    def gain(frequency):
        env = {system.s: 1j*2*np.pi*frequency, vin: 1.0}
        M = np.array([[complex(rows[i, j].subs(env)) for j in range(n)]
                      for i in range(n)], dtype=complex)
        b = np.array([complex(source[i].subs(env)) for i in range(n)])
        return abs(np.linalg.solve(M, b)[system.out_index])

    ## Doubly-terminated ladders have a passband gain of one half.
    assert abs(gain(10.0) - 0.5) < 0.02
    assert abs(gain(1e3)/gain(10.0) - 1.0) < 0.05, 'passband must be flat'

    ## Fifth order: 100 dB per decade, measured across an octave well into
    ## the stopband.
    decade = 20*np.log10(gain(2e5)/gain(1e5))/np.log10(2.0)
    assert -110 < decade < -90, 'expected ~-100 dB/decade, got %.1f' % decade


def test_leapfrog_reuses_the_calibration_amplifier_unchanged():
    """Five µA741s, and the single-amplifier fixture must be untouched.

    ``ua741()`` is the circuit the published DDD sizes are calibrated
    against, so extracting its device stamping for reuse had to leave its
    matrix bit-identical.  This pins the shapes that would move if it did
    not: node count, nonzero count, and that the leapfrog is five amplifiers
    plus its passive network.
    """
    from pycircuit.circuit import benchmark_circuits as bc

    single = bc.ua741()
    assert single.dim == 26
    assert sum(1 for e in single.A if e != 0) == 103

    filt = bc.leapfrog_5th_order()
    ## 5 amplifiers of 25 nodes each, the input node, and one more unknown:
    ## MNA carries a branch current for the voltage source.
    assert filt.dim == 5*len(bc._UA741_NODE_NAMES) + 2
    assert filt.dim == 127
    assert single.dim == len(bc._UA741_NODE_NAMES) + 1, (
        'the single amplifier is its own nodes plus the source branch')


def test_fully_symbolic_ua741_substitutes_back_to_the_numeric_one():
    """``fully_symbolic=True`` must change the symbols and nothing else.

    The option exists for cancellation work: a device parameter left numeric
    merges into its matrix entry's arithmetic, so its cancelling partner at the
    mirrored position becomes undetectable.  With the default fixture every one of
    the matrix's additive contributions is a pure number, which is why a
    de-cancelling expansion has nothing to see there.

    The property that makes it trustworthy is that it is the *same circuit*:
    substituting the recorded nominal values must reproduce the numeric fixture
    entry for entry.
    """
    import numpy as np
    import sympy
    from pycircuit.circuit import benchmark_circuits as bc

    num = bc.ua741()
    sym = bc.ua741(fully_symbolic=True)

    assert sym.dim == num.dim
    assert sum(1 for e in sym.A if e != 0) == sum(1 for e in num.A if e != 0)
    ## Every transistor parameter plus every passive, against gm alone.
    assert len(sym.A.free_symbols) > 100
    assert len(num.A.free_symbols) == 1              # only s

    probe = {num.s: 2j * np.pi * 1e3}
    back = sym.A.subs({k: v for k, v in sym.params.items()})
    for i in range(num.dim):
        for j in range(num.dim):
            a, b = num.A[i, j], back[i, j]
            va = complex(a.subs(probe)) if getattr(a, 'free_symbols', None) else complex(a)
            vb = complex(b.subs({sym.s: probe[num.s]})) if getattr(b, 'free_symbols', None) else complex(b)
            scale = max(abs(va), abs(vb), 1e-30)
            assert abs(va - vb) / scale < 1e-12, (i, j, va, vb)


def test_fully_symbolic_ua741_leaves_only_structural_numbers():
    """Only the source's incidence entries may stay numeric.

    That is the measurable statement of "one symbol per device": after the
    change, the sole additive contributions without a device symbol should be the
    voltage source's +-1 incidence pair, which is topology rather than a device.
    """
    import sympy
    from pycircuit.circuit import benchmark_circuits as bc

    sym = bc.ua741(fully_symbolic=True)
    devices = {s for s in sym.A.free_symbols if s is not sym.s}
    numeric_only = 0
    for i in range(sym.dim):
        for j in range(sym.dim):
            if sym.A[i, j] == 0:
                continue
            for term in sympy.Add.make_args(sympy.expand(sym.A[i, j])):
                if not (term.free_symbols & devices):
                    numeric_only += 1
    assert numeric_only == 2, numeric_only
