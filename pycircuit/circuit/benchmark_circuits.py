# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Reference circuits for comparing symbolic solve/representation backends.

These are *fixtures*, deliberately fixed in advance so that a backend comparison
cannot be chosen after seeing results (see ``doc/ddd_plan.md``, Stage B).  Every
circuit here is built with the ordinary :class:`~pycircuit.circuit.SubCircuit`
API and stamped through the normal MNA path, so what a backend is measured on is
exactly what an analysis would hand it -- never a hand-written matrix.

The one exception is :func:`legacy_soe_ladder`, which reproduces the hand-built
matrix used by the published SoE table so the harness can be validated against
known numbers.  It is a calibration aid, not a benchmark circuit.

Each builder returns a :class:`BenchSystem`: the circuit, the linear system
``A x = b`` at complex frequency ``s`` (reference node removed), and the indices
needed to form a transfer function.
"""

import numpy as np
import sympy

from pycircuit.circuit import SubCircuit, R, C, VS, IS, VCCS, Nullor, gnd
from pycircuit.circuit.toolkit import symbolic
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.analysis_ss import dc_steady_state


class BenchSystem:
    """A linear system to hand to a symbolic backend, plus its provenance.

    Attributes:
        name: Short identifier used in result tables.
        cir: The originating circuit, or ``None`` for pure-matrix fixtures.
        A: ``s*C + G`` with the reference node removed (sympy ``Matrix``).
        b: Right-hand side ``-u`` (sympy ``Matrix``, column).
        s: The complex frequency symbol.
        params: Numeric substitution for every free symbol except ``s``,
            used for correctness checks and numeric evaluation.
        out_index, in_index: Row indices of the output and input unknowns, so a
            transfer function is ``x[out_index] / x[in_index]``.
    """

    def __init__(self, name, A, b, s, params, out_index, in_index, cir=None):
        self.name = name
        self.cir = cir
        self.A = A
        self.b = b
        self.s = s
        self.params = params
        self.out_index = out_index
        self.in_index = in_index

    @property
    def dim(self):
        return self.A.rows

    def __repr__(self):
        return '<BenchSystem %s dim=%d>' % (self.name, self.dim)


def system_from_circuit(cir, name, params, out_index, in_index, s=None):
    """Stamp ``cir`` through the normal MNA path and return a `BenchSystem`.

    This is the only route by which a circuit becomes a matrix here -- it uses
    :func:`~pycircuit.circuit.analysis_ss.dc_steady_state` exactly as
    :class:`~pycircuit.circuit.analysis_ss.AC` does, so the matrix under test is
    the analysis's matrix.
    """
    s = s if s is not None else sympy.Symbol('s')
    G, Cm, CY, u, x, ss = dc_steady_state(cir, s, gnd, symbolic, complexfreq=True)
    irefnode = cir.get_node_index(gnd)
    G, Cm, CY, u = remove_row_col((G, Cm, CY, u), irefnode, symbolic)
    A = sympy.Matrix(s * Cm + G)
    b = -sympy.Matrix(u)
    return BenchSystem(name, A, b, s, params, out_index, in_index, cir=cir)


def rc_ladder(N):
    """Fully symbolic ``N``-section RC ladder, driven by a voltage source.

    Section ``i`` is a series ``R{i}`` to the next node and a shunt ``C{i}`` to
    ground; every component is a distinct symbol.  This is the workhorse case:
    fully symbolic, sparse, and with a structure whose exploded ``N(s)/D(s)`` is
    exponentially large while a shared representation stays small.

    Note the system is ``(N+1)x(N+1)``, not ``NxN`` -- MNA carries the voltage
    source's branch current as an extra unknown.  It is therefore *not* the same
    system as :func:`legacy_soe_ladder`.
    """
    cir = SubCircuit(toolkit=symbolic)
    nodes = [cir.add_node('n%d' % i) for i in range(N)]
    Rs = [sympy.Symbol('R%d' % i, positive=True) for i in range(N - 1)]
    Cs = [sympy.Symbol('C%d' % i, positive=True) for i in range(N)]
    vin = sympy.Symbol('vin')

    cir['vs'] = VS(nodes[0], gnd, vac=vin)
    for i in range(N - 1):
        cir['R%d' % i] = R(nodes[i], nodes[i + 1], r=Rs[i])
    for i in range(N):
        cir['C%d' % i] = C(nodes[i], gnd, c=Cs[i])

    params = {Rs[i]: 100.0 * (i + 1) for i in range(N - 1)}
    params.update({Cs[i]: 1e-9 * (i + 1) for i in range(N)})
    params[vin] = 1.0

    # Node ordering follows the circuit's own node list; the reference node is
    # removed above, and gnd is appended last by add_node, so n0..n{N-1} keep
    # indices 0..N-1.
    return system_from_circuit(cir, 'rc_ladder_N%d' % N, params,
                               out_index=N - 1, in_index=0)


def mfb_filter():
    """The multiple-feedback filter of Example 10 -- exercises a `Nullor`.

    Included because a nullor stamps unit entries with no diagonal contribution,
    giving a structurally different sparsity pattern from the passive ladder and
    a classic source of cancelling terms.  Fully symbolic.
    """
    R1, R2, R3 = sympy.symbols('R1 R2 R3', positive=True)
    C1, C2 = sympy.symbols('C1 C2', positive=True)
    i_s = sympy.Symbol('i_s')

    cir = SubCircuit(toolkit=symbolic)
    n1, n2, n3 = (cir.add_node(n) for n in ('n1', 'n2', 'n3'))
    cir['R1'] = R(n1, n3, r=R1)
    cir['R2'] = R(n1, n2, r=R2)
    cir['R3'] = R(n1, gnd, r=R3)
    cir['C1'] = C(n1, gnd, c=C1)
    cir['C2'] = C(n2, n3, c=C2)
    cir['nullor'] = Nullor(n2, gnd, n3, gnd)
    cir['is'] = IS(n1, gnd, iac=i_s)

    params = {R1: 1e3, R2: 2e3, R3: 4e3, C1: 1e-9, C2: 2e-9, i_s: 1.0}
    return system_from_circuit(cir, 'mfb_filter', params,
                               out_index=2, in_index=0)


def rc_ladder_semi_symbolic(N, n_symbolic=2):
    """RC ladder with only ``n_symbolic`` components left symbolic.

    The semi-symbolic middle: most parameters numeric, a few of interest kept as
    symbols, plus ``s``.  This is the regime numeric terminals target, and the
    one where the GiNaC backend previously stalled on exact-rational growth.
    """
    cir = SubCircuit(toolkit=symbolic)
    nodes = [cir.add_node('n%d' % i) for i in range(N)]
    vin = sympy.Symbol('vin')
    params = {vin: 1.0}

    cir['vs'] = VS(nodes[0], gnd, vac=vin)
    for i in range(N - 1):
        if i < n_symbolic:
            sym = sympy.Symbol('R%d' % i, positive=True)
            params[sym] = 100.0 * (i + 1)
            rval = sym
        else:
            rval = 100.0 * (i + 1)
        cir['R%d' % i] = R(nodes[i], nodes[i + 1], r=rval)
    for i in range(N):
        cir['C%d' % i] = C(nodes[i], gnd, c=1e-9 * (i + 1))

    return system_from_circuit(cir, 'rc_ladder_semi_N%d_k%d' % (N, n_symbolic),
                               params, out_index=N - 1, in_index=0)


def dense_symbolic_matrix(n):
    """A dense ``n x n`` matrix of distinct symbols -- *not* a circuit.

    Used only for the analytic identity from Shi (TCAS-II 2010): under a
    row-wise expansion order the optimal DDD size for a full matrix is exactly
    ``n * 2**(n-1)``.  That makes this the one fixture with an exact expected
    answer, which is why it is the primary structural check that sharing works.
    """
    a = [[sympy.Symbol('a%d_%d' % (i, j)) for j in range(n)] for i in range(n)]
    A = sympy.Matrix(a)
    b = sympy.zeros(n, 1)
    b[0] = 1

    ## Values come from a seeded RNG rather than a closed-form pattern.  An
    ## arithmetic pattern (e.g. ``(i*n+j) % 7 + 1``) makes the matrix regular
    ## enough that an unpivoted elimination can land on an exactly-zero pivot --
    ## which is a real limitation of division-based solvers, but not one a
    ## fixture should trigger by accident.  Diagonally dominant here, so every
    ## backend sees a well-conditioned, non-degenerate system.
    rng = np.random.default_rng(20260727 + n)
    params = {}
    for i in range(n):
        for j in range(n):
            v = float(rng.uniform(1.0, 9.0))
            if i == j:
                v += n
            params[a[i][j]] = v
    return BenchSystem('dense_%d' % n, A, b, sympy.Symbol('s'), params,
                       out_index=n - 1, in_index=0)


def legacy_soe_ladder(N):
    """The hand-built matrix behind the published SoE op-count table.

    Reproduces ``doc/prototypes/soe_symbolic.py::_ladder`` so the harness can be
    checked against known numbers (73/157/241/325 ops at N=4/8/12/16).  It is a
    Norton-driven ``NxN`` system with no explicit source branch, so it is *not*
    equivalent to :func:`rc_ladder`; it exists purely to validate the
    measurement code, not to benchmark backends.
    """
    s = sympy.Symbol('s')
    Rs = [sympy.Symbol('R%d' % i, positive=True) for i in range(N)]
    Cs = [sympy.Symbol('C%d' % i, positive=True) for i in range(N)]
    A = sympy.zeros(N, N)
    b = sympy.zeros(N, 1)
    for i in range(N):
        A[i, i] = 1 / Rs[i] + s * Cs[i] + (1 / Rs[i - 1] if i > 0 else 0)
        if i + 1 < N:
            A[i, i + 1] = -1 / Rs[i]
            A[i + 1, i] = -1 / Rs[i]
    b[0] = 1 / Rs[0]
    params = {Rs[i]: 100.0 * (i + 1) for i in range(N)}
    params.update({Cs[i]: 1e-9 * (i + 1) for i in range(N)})
    return BenchSystem('legacy_soe_ladder_N%d' % N, A, b, s, params,
                       out_index=N - 1, in_index=0)


## The fixed benchmark set.  Adding to this list is fine; changing or removing
## an entry invalidates comparisons already published, so don't.
LADDER_SIZES = (4, 8, 12, 16)


def standard_suite():
    """The benchmark set fixed by Stage B, as ``(group, BenchSystem)`` pairs."""
    suite = [('ladder', rc_ladder(N)) for N in LADDER_SIZES]
    suite += [('semi', rc_ladder_semi_symbolic(N)) for N in LADDER_SIZES]
    suite += [('mfb', mfb_filter())]
    suite += [('dense', dense_symbolic_matrix(n)) for n in range(3, 7)]
    return suite


## ------------------------------------------------ paper calibration case --

def add_small_signal_bjt(cir, name, b, c, e, gm, rpi, ro, cpi=None, cmu=None):
    """Stamp a hybrid-pi small-signal BJT into ``cir``.

    Built entirely from primitives pycircuit already has -- a transconductance
    (`VCCS`), base and output resistances, and the two junction capacitances --
    so no new element type is needed to model an amplifier symbolically.

    Args:
        cir: Circuit to add to.
        name: Instance prefix for the stamped elements.
        b, c, e: Base, collector and emitter nodes.
        gm, rpi, ro: Transconductance, base-emitter and output resistance.
        cpi, cmu: Optional junction capacitances; omit for a DC-only model.
    """
    cir[name + '_rpi'] = R(b, e, r=rpi)
    cir[name + '_ro'] = R(c, e, r=ro)
    cir[name + '_gm'] = VCCS(b, e, c, e, gm=gm)
    if cpi is not None:
        cir[name + '_cpi'] = C(b, e, c=cpi)
    if cmu is not None:
        cir[name + '_cmu'] = C(b, c, c=cmu)


def opamp_741_like(symbolic_devices=(), miller=True):
    """A µA741-class operational amplifier, small-signal.

    **Not a bit-exact µA741 netlist.**  It follows the same signal path at
    comparable node count -- differential input pair into common-base stages, a
    degenerated current-mirror active load, an emitter-follower into a
    Miller-compensated common-emitter gain stage, and a push-pull output -- so
    that the diagram sizes it produces can be compared *in order of magnitude*
    with published figures for a real µA741.

    Those figures depend on the authors' own netlist, symbol convention and
    partitioning, so exact agreement is not the expectation.  This is
    calibration, not a target: landing within a small factor says the
    implementation is sound, being orders out says it is not.  See
    ``doc/ddd_conclusions.md``.

    Args:
        symbolic_devices: Names of devices whose ``gm`` stays a symbol.  Empty
            gives an entirely numeric small-signal circuit (plus ``s``).
        miller: Include the compensation capacitor.

    Returns:
        BenchSystem.
    """
    ## VCCS builds its stamp at construction time from ``default_toolkit``, not
    ## from the parent's toolkit, so a symbolic ``gm`` needs the global set while
    ## the circuit is assembled -- the same thing the symbolic analysis tests do.
    import pycircuit.circuit.circuit as _circuit_module
    _saved_toolkit = _circuit_module.default_toolkit
    _circuit_module.default_toolkit = symbolic
    try:
        return _build_opamp_741_like(symbolic_devices, miller)
    finally:
        _circuit_module.default_toolkit = _saved_toolkit


def _build_opamp_741_like(symbolic_devices, miller):
    cir = SubCircuit(toolkit=symbolic)
    names = ('inp', 'inn', 'e1', 'e2', 'c3', 'c4', 'e5', 'e6',
             'e16', 'c17', 'e17', 'e23', 'out')
    n = {name: cir.add_node(name) for name in names}

    vin = sympy.Symbol('vin')
    params = {vin: 1.0}
    cir['vs'] = VS(n['inp'], gnd, vac=vin)
    cir['rinn'] = R(n['inn'], gnd, r=1e3)

    def gm_of(dev, nominal):
        if dev in symbolic_devices:
            sym = sympy.Symbol('gm_%s' % dev, positive=True)
            params[sym] = nominal
            return sym
        return nominal

    ## Input pair: emitter followers into common-base stages (the 741's cascode
    ## input), sharing degenerated tails.
    add_small_signal_bjt(cir, 'q1', n['inp'], gnd, n['e1'],
                         gm_of('q1', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    add_small_signal_bjt(cir, 'q2', n['inn'], gnd, n['e2'],
                         gm_of('q2', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    add_small_signal_bjt(cir, 'q3', gnd, n['c3'], n['e1'],
                         gm_of('q3', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    add_small_signal_bjt(cir, 'q4', gnd, n['c4'], n['e2'],
                         gm_of('q4', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    cir['rt1'] = R(n['e1'], gnd, r=500e3)
    cir['rt2'] = R(n['e2'], gnd, r=500e3)

    ## Degenerated current-mirror load: q5 diode-connected, q6 mirrors it.
    add_small_signal_bjt(cir, 'q5', n['c3'], n['c3'], n['e5'],
                         gm_of('q5', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    add_small_signal_bjt(cir, 'q6', n['c3'], n['c4'], n['e6'],
                         gm_of('q6', 2e-3), 25e3, 1e6, 3e-12, 0.5e-12)
    cir['r5'] = R(n['e5'], gnd, r=1e3)
    cir['r6'] = R(n['e6'], gnd, r=1e3)

    ## Emitter follower into the Miller-compensated gain stage.
    add_small_signal_bjt(cir, 'q16', n['c4'], gnd, n['e16'],
                         gm_of('q16', 3e-3), 20e3, 500e3, 5e-12, 1e-12)
    cir['r16'] = R(n['e16'], gnd, r=50e3)
    add_small_signal_bjt(cir, 'q17', n['e16'], n['c17'], n['e17'],
                         gm_of('q17', 5e-3), 10e3, 200e3, 10e-12, 1e-12)
    cir['r17'] = R(n['e17'], gnd, r=100)
    cir['rl17'] = R(n['c17'], gnd, r=100e3)
    if miller:
        cir['cc'] = C(n['e16'], n['c17'], c=30e-12)

    ## Output: follower pair driving the load.
    add_small_signal_bjt(cir, 'q23', n['c17'], gnd, n['e23'],
                         gm_of('q23', 8e-3), 8e3, 100e3, 8e-12, 1e-12)
    cir['r23'] = R(n['e23'], gnd, r=50e3)
    add_small_signal_bjt(cir, 'q14', n['e23'], gnd, n['out'],
                         gm_of('q14', 20e-3), 5e3, 50e3, 20e-12, 2e-12)
    cir['rload'] = R(n['out'], gnd, r=2e3)

    return system_from_circuit(cir, 'opamp_741_like', params,
                               out_index=names.index('out'),
                               in_index=names.index('inp'))
