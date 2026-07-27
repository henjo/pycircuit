# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Tests for the DDD family of analyses.

The point of these is as much structural as numerical.  The DDD analyses must
agree with the originals -- but they must also *not reimplement* them, because
two copies of the port bookkeeping or the return-ratio derivation would drift
apart and the originals would stop being a usable reference.  So there are tests
for the agreement and tests for the thinness.
"""

import numpy as np
import pytest
import sympy

import pycircuit.circuit.circuit as circuit_module
from pycircuit.circuit import SubCircuit, R, VS, VCVS, gnd
from pycircuit.circuit import benchmark_circuits as bc
from pycircuit.circuit.dddanalysis import (DDDAnalysisMixin, DDDLoopGain,
                                           DDDTransimpedance, DDDTwoPort)
from pycircuit.circuit.feedback import FeedbackLoopAnalysis, LoopProbe
from pycircuit.circuit.nportanalysis import TwoPortAnalysis
from pycircuit.circuit.toolkit import ddd_toolkit, symbolic


@pytest.fixture
def symbolic_default():
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = symbolic
    yield
    circuit_module.default_toolkit = saved


## -- the primitives the reuse rests on ------------------------------------

def test_det_and_cofactor_agree_with_sympy():
    """Overriding these two is the whole backend of this stage."""
    A = bc.rc_ladder(4).A
    assert sympy.simplify(ddd_toolkit.det(A) - symbolic.det(A)) == 0
    for i in range(A.rows):
        for j in range(A.cols):
            assert sympy.simplify(ddd_toolkit.cofactor(A, i, j)
                                  - symbolic.cofactor(A, i, j)) == 0


def test_one_family_serves_the_determinant_and_every_cofactor():
    """A two-port asks for a determinant then cofactors of the *same* matrix.

    Rebuilding per call would discard exactly the sharing that makes a diagram
    worth having, so the toolkit keeps the family.
    """
    A = bc.rc_ladder(5).A
    ddd_toolkit.det(A)
    first = ddd_toolkit._family_cache[1]
    for i in range(A.rows):
        ddd_toolkit.cofactor(A, i, 0)
    assert ddd_toolkit._family_cache[1] is first


def test_a_different_matrix_gets_a_different_family():
    ddd_toolkit.det(bc.rc_ladder(4).A)
    first = ddd_toolkit._family_cache[1]
    ddd_toolkit.det(bc.rc_ladder(5).A)
    assert ddd_toolkit._family_cache[1] is not first


## -- agreement with the originals -----------------------------------------

def test_two_port_matches_the_original(symbolic_default):
    R1, R2 = sympy.symbols('R1 R2', real=True)
    s = sympy.Symbol('s')

    def build():
        cir = SubCircuit(toolkit=symbolic)
        n1, n2 = cir.add_nodes('n1', 'n2')
        cir['R1'] = R(n1, n2, r=R1)
        cir['R2'] = R(n2, gnd, r=R2)
        return cir, n1, n2

    cir, n1, n2 = build()
    reference = TwoPortAnalysis(cir, n1, gnd, n2, gnd,
                                toolkit=symbolic).solve(s, complexfreq=True)
    cir, n1, n2 = build()
    got = DDDTwoPort(cir, n1, gnd, n2, gnd).solve(s, complexfreq=True)

    for name in ('mu', 'gamma', 'zeta', 'beta'):
        assert sympy.simplify(sympy.simplify(got[name])
                              - sympy.simplify(reference[name])) == 0


def test_loop_gain_matches_the_original(symbolic_default):
    A, R1, R2 = sympy.symbols('A R1 R2')
    s = sympy.Symbol('s')

    def build():
        cir = SubCircuit(toolkit=symbolic)
        cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g=A, toolkit=symbolic)
        cir['R1'] = R('in', 'int', r=R1)
        cir['probe'] = LoopProbe('out', gnd, 'out_R2', gnd, toolkit=symbolic)
        cir['R2'] = R('int', 'out_R2', r=R2)
        cir['VS'] = VS('in', gnd)
        return cir

    reference = FeedbackLoopAnalysis(build(),
                                     toolkit=symbolic).solve(s, complexfreq=True)
    got = DDDLoopGain(build()).solve(s, complexfreq=True)

    assert sympy.simplify(sympy.simplify(got['loopgain'])
                          - sympy.simplify(reference['loopgain'])) == 0
    ## And it is the textbook answer, not merely self-consistent.
    assert sympy.simplify(got['loopgain'] + A * R1 / (R1 + R2)) == 0


## -- thinness --------------------------------------------------------------

@pytest.mark.parametrize('cls', [DDDTwoPort, DDDLoopGain, DDDTransimpedance])
def test_the_frontends_do_not_reimplement_the_analysis(cls):
    """Nothing in the DDD family redefines ``solve``.

    If one ever does, the reuse route has been abandoned and there are two
    implementations to keep in step -- which is the thing this design exists to
    avoid.
    """
    assert 'solve' not in cls.__dict__
    assert issubclass(cls, DDDAnalysisMixin)


@pytest.mark.parametrize('cls', [DDDTwoPort, DDDLoopGain, DDDTransimpedance])
def test_the_frontends_bind_the_ddd_toolkit(cls):
    assert cls.__mro__.index(DDDAnalysisMixin) < len(cls.__mro__)


def test_an_explicit_toolkit_still_wins(symbolic_default):
    """Binding a default must not take the choice away."""
    cir = SubCircuit(toolkit=symbolic)
    n1, n2 = cir.add_nodes('n1', 'n2')
    cir['R1'] = R(n1, n2, r=sympy.Symbol('R1', real=True))
    cir['R2'] = R(n2, gnd, r=sympy.Symbol('R2', real=True))
    analysis = DDDTwoPort(cir, n1, gnd, n2, gnd, toolkit=symbolic)
    assert analysis.toolkit is symbolic


def test_diagram_accessors_expose_what_the_base_class_cannot(symbolic_default):
    """The one thing the frontends add: a window onto the representation."""
    cir = SubCircuit(toolkit=symbolic)
    n1, n2 = cir.add_nodes('n1', 'n2')
    cir['R1'] = R(n1, n2, r=sympy.Symbol('R1', real=True))
    cir['R2'] = R(n2, gnd, r=sympy.Symbol('R2', real=True))

    analysis = DDDTwoPort(cir, n1, gnd, n2, gnd)
    A = bc.rc_ladder(5).A
    assert analysis.family(A).denominator.size > 0
    assert analysis.determinant_diagram(A).term_count() > 0


def test_loop_sensitivity_comes_from_the_shared_family(symbolic_default):
    """Return ratio is built from a derivative the family already holds."""
    cir = SubCircuit(toolkit=symbolic)
    cir['A1'] = VCVS(gnd, 'int', 'out', gnd, g=sympy.Symbol('A'),
                     toolkit=symbolic)
    cir['R1'] = R('in', 'int', r=sympy.Symbol('R1'))
    cir['probe'] = LoopProbe('out', gnd, 'out_R2', gnd, toolkit=symbolic)
    cir['R2'] = R('int', 'out_R2', r=sympy.Symbol('R2'))
    cir['VS'] = VS('in', gnd)

    analysis = DDDLoopGain(cir)
    system = bc.rc_ladder(4)
    parameter = sorted(system.A.free_symbols - {system.s}, key=str)[0]
    combination = analysis.loop_sensitivity(system.A, parameter)
    assert sympy.simplify(combination.eval()
                          - sympy.diff(system.A.det(), parameter)) == 0
