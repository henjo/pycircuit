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


## Nominal small-signal parameters by device role.  Values are representative of
## a bipolar op-amp biased at a few tens of microamps; the point of this fixture
## is the topology and matrix structure, not predicting a datasheet number.
_UA741_ROLES = {
    ## role          gm      rpi     ro      cpi     cmu
    'input':      (0.4e-3,  250e3,  2.0e6,  3e-12,  0.5e-12),
    'mirror':     (0.4e-3,  250e3,  2.0e6,  3e-12,  0.5e-12),
    'bias':       (0.6e-3,  180e3,  1.5e6,  3e-12,  0.5e-12),
    'gain':       (2.0e-3,   50e3,  400e3,  8e-12,  1.0e-12),
    'output':     (8.0e-3,   12e3,  100e3, 20e-12,  2.0e-12),
    ## Protection devices sit cut off in normal operation: negligible
    ## transconductance, very high impedances.  They still contribute nodes and
    ## a little capacitive coupling, which is why they are not simply omitted.
    'off':        (1e-9,      1e9,    1e9,  0.5e-12, 0.2e-12),
}


def ua741(symbolic_devices=(), miller=True, fully_symbolic=False):
    """The µA741 operational amplifier, small-signal -- the calibration circuit.

    Transcribed from Fig. 15 of Tan & Shi, *Hierarchical Symbolic Analysis of
    Analog Integrated Circuits via Determinant Decision Diagrams* (IEEE TCAD
    19(4), 2000), the circuit their published DDD sizes are measured on: 24
    transistors and 11 resistors, reported there as a 24x24 MNA matrix with 89
    nonzeros, a flat diagram of 6654 vertices standing for 119 011 product terms.

    Supply rails are AC grounds, and each transistor is the hybrid-pi model of
    :func:`add_small_signal_bjt`, so the matrix structure follows the schematic
    while the device values are representative rather than taken from a DC
    operating point.  Q13's split collector is modelled as a single collector.

    Args:
        symbolic_devices: Device names (``'q1'`` ... ``'q24'``) whose ``gm``
            stays symbolic.  Empty gives a numeric circuit plus ``s``.
        miller: Include the 30 pF compensation capacitor.
        fully_symbolic: Give **every** device parameter its own symbol -- each
            transistor's ``gm``, ``rpi``, ``ro``, ``cpi``, ``cmu`` and every
            resistor and capacitor -- rather than only the named
            transconductances.  Overrides ``symbolic_devices``.

            This is what a *cancellation-free* expansion needs, and the reason is
            not aesthetic: term cancellation in an MNA determinant is a device
            appearing at two matrix positions with opposite signs, so a parameter
            left numeric merges into its entry's arithmetic and its cancelling
            partner becomes undetectable.  With the default fixture 159 of the
            matrix's 215 additive contributions are pure numbers.  See
            ``doc/cancellation_ranking_conclusions.md`` §18.

            The cost is a far larger symbol set (roughly 130 symbols against 24),
            so a fully symbolic diagram is much bigger; this is for cancellation
            work, not for the calibration sizes.

    Returns:
        BenchSystem.
    """
    import pycircuit.circuit.circuit as _circuit_module
    saved = _circuit_module.default_toolkit
    _circuit_module.default_toolkit = symbolic
    try:
        return _build_ua741(symbolic_devices, miller, fully_symbolic)
    finally:
        _circuit_module.default_toolkit = saved


def _stamp_ua741_devices(cir, n, bjt, miller, pfx='', res=None, cap=None):
    """Stamp the µA741's devices, given nodes already created.

    Split out of :func:`_build_ua741` so the amplifier can be
    instantiated several times in one circuit -- the leapfrog filter
    needs five.  Node *creation* deliberately stays with the caller:
    node order fixes the matrix row order, and ``ua741()`` is the
    calibration circuit whose published diagram sizes depend on it.

    Args:
        cir: Circuit to stamp into.
        n: Mapping from the µA741 node names to nodes.
        bjt: ``bjt(device, role, b, c, e)``, supplied by the caller so
            it decides which transconductances stay symbolic.
        miller: Include the 30 pF compensation capacitor.
        pfx: Element-name prefix, so instances do not collide.
        res, cap: ``res(name, a, b, value)`` / ``cap(name, a, b, value)``, the
            same hook idea as ``bjt`` extended to the passives.  Default to
            stamping the numeric value.  A caller that wants *every* device
            parameter symbolic -- which is what de-cancellation needs, since a
            numeric contribution merges into its matrix entry and its
            cancellation becomes invisible -- supplies all three.
    """
    def _numeric_res(name, a, b, value):
        cir[name] = R(a, b, r=value)

    def _numeric_cap(name, a, b, value):
        cir[name] = C(a, b, c=value)

    res = res if res is not None else _numeric_res
    cap = cap if cap is not None else _numeric_cap
    ## -- input stage: emitter followers Q1/Q2 into common-base Q3/Q4 --------
    ## Q1/Q2 are followers whose emitters drive the common-base pair Q3/Q4 --
    ## two separate signal paths, not a shared tail node.
    ## Q1/Q2 are followers: their collectors carry bias, not signal, and tie
    ## together at the Q8 collector.  Splitting them into two nodes is both an
    ## extra unknown and electrically wrong.
    bjt('q1', 'input', n['inp'], n['c12'], n['e1'])
    bjt('q2', 'input', n['inn'], n['c12'], n['e2'])
    bjt('q3', 'input', n['nb34'], n['c3'], n['e1'])
    bjt('q4', 'input', n['nb34'], n['c4'], n['e2'])

    ## Q8/Q9 mirror loading the input pair collectors.
    bjt('q8', 'mirror', n['nq9'], n['c12'], gnd)
    bjt('q9', 'mirror', n['nq9'], n['nq9'], gnd)
    res(pfx + 'rq9', n['c12'], n['nq9'], 1e3)

    ## -- Q5/Q6 active load with Q7 buffering the mirror base ----------------
    bjt('q5', 'mirror', n['nb56'], n['c3'], n['e5'])
    bjt('q6', 'mirror', n['nb56'], n['c4'], n['e6'])
    bjt('q7', 'mirror', n['c3'], gnd, n['nb56'])
    res(pfx + 'R1', n['e5'], gnd, 1e3)
    res(pfx + 'R2', n['nb56'], gnd, 50e3)
    res(pfx + 'R3', n['e6'], gnd, 1e3)

    ## -- Widlar bias chain: Q10/Q11 with R4, mirrored by Q12 into Q13 -------
    bjt('q10', 'bias', n['nb1011'], n['nb34'], n['e10'])
    bjt('q11', 'bias', n['nb1011'], n['nb1011'], gnd)
    bjt('q12', 'bias', n['nr5'], n['nr5'], gnd)
    res(pfx + 'R4', n['e10'], gnd, 5e3)
    res(pfx + 'R5', n['nr5'], n['nb1011'], 39e3)

    ## -- second stage: Q16 follower into Q17, Miller compensated ------------
    bjt('q16', 'gain', n['c4'], gnd, n['e16'])
    bjt('q17', 'gain', n['e16'], n['c17'], n['e17'])
    res(pfx + 'R9', n['e16'], gnd, 50e3)
    res(pfx + 'R8', n['e17'], gnd, 100)
    bjt('q13', 'bias', n['nr5'], n['c17'], gnd)    # current-source load
    if miller:
        cap(pfx + 'cc', n['c4'], n['c17'], 30e-12)

    ## -- output stage: Q23 driver, Q18/Q19 Vbe multiplier, Q14/Q20 pair -----
    bjt('q23', 'gain', n['c17'], gnd, n['e23'])
    res(pfx + 'R11', n['e23'], gnd, 50e3)
    bjt('q18', 'bias', n['n19'], n['nb14'], n['nb20'])
    bjt('q19', 'bias', n['nb20'], n['n19'], gnd)
    res(pfx + 'R10', n['n19'], n['nb20'], 40e3)
    res(pfx + 'rdrv', n['e23'], n['nb14'], 1e3)

    bjt('q14', 'output', n['nb14'], gnd, n['nr6'])
    bjt('q20', 'output', n['nb20'], gnd, n['nr7'])
    res(pfx + 'R6', n['nr6'], n['out'], 27)
    res(pfx + 'R7', n['nr7'], n['out'], 22)

    ## -- short-circuit protection and remaining bias devices ----------------
    bjt('q15', 'off', n['nr6'], n['nb14'], gnd)
    bjt('q21', 'off', n['nr7'], n['nb20'], gnd)
    bjt('q22', 'off', n['e17'], n['c17'], gnd)
    bjt('q24', 'off', n['nb1011'], n['nb1011'], gnd)

    res(pfx + 'rload', n['out'], gnd, 2e3)


def _build_ua741(symbolic_devices, miller, fully_symbolic=False):
    cir = SubCircuit(toolkit=symbolic)
    names = (
        'inp', 'inn',            # inputs
        'e1', 'e2', 'c12',       # input pair emitters, common collector node
        'nb34', 'nq9',           # cascode bases, Q8/Q9 mirror diode node
        'c3', 'c4',              # first-stage outputs (c4 is the high-Z node)
        'nb56', 'e5', 'e6',      # load mirror base node and degeneration
        'nb1011', 'e10', 'nr5',  # Widlar bias chain, Q12/Q13 mirror node
        'e16', 'e17', 'c17',     # second stage
        'e23', 'nb14', 'nb20',   # output drive
        'n19', 'nr6', 'nr7',     # Vbe multiplier, output degeneration
        'out',
    )
    n = {name: cir.add_node(name) for name in names}

    vin = sympy.Symbol('vin')
    params = {vin: 1.0}
    cir['vs'] = VS(n['inp'], gnd, vac=vin)
    ## Stamped here rather than in the helper, so it needs the same treatment.
    cir['rinn'] = R(n['inn'], gnd,
                    r=sympy.Symbol('r_rinn', positive=True)
                    if fully_symbolic else 1e6)
    if fully_symbolic:
        params[sympy.Symbol('r_rinn', positive=True)] = 1e6

    def symbolise(value, kind, dev):
        """One symbol per device parameter, with its nominal value in params."""
        sym = sympy.Symbol('%s_%s' % (kind, dev), positive=True)
        params[sym] = value
        return sym

    def bjt(dev, role, b, c, e):
        gm, rpi, ro, cpi, cmu = _UA741_ROLES[role]
        if fully_symbolic:
            gm = symbolise(gm, 'gm', dev)
            rpi = symbolise(rpi, 'rpi', dev)
            ro = symbolise(ro, 'ro', dev)
            cpi = None if cpi is None else symbolise(cpi, 'cpi', dev)
            cmu = None if cmu is None else symbolise(cmu, 'cmu', dev)
        elif dev in symbolic_devices:
            gm = symbolise(gm, 'gm', dev)
        add_small_signal_bjt(cir, dev, b, c, e, gm, rpi, ro, cpi, cmu)

    def res(name, a, b, value):
        cir[name] = R(a, b, r=symbolise(value, 'r', name)
                      if fully_symbolic else value)

    def cap(name, a, b, value):
        cir[name] = C(a, b, c=symbolise(value, 'c', name)
                      if fully_symbolic else value)

    _stamp_ua741_devices(cir, n, bjt, miller, res=res, cap=cap)

    return system_from_circuit(cir, 'ua741', params,
                               out_index=names.index('out'),
                               in_index=names.index('inp'))


_UA741_NODE_NAMES = (
    'inp', 'inn', 'e1', 'e2', 'c12', 'nb34', 'nq9', 'c3', 'c4',
    'nb56', 'e5', 'e6', 'nb1011', 'e10', 'nr5', 'e16', 'e17', 'c17',
    'e23', 'nb14', 'nb20', 'n19', 'nr6', 'nr7', 'out',
)


def leapfrog_5th_order(symbolic_devices=(), miller=True):
    """A 5th-order leapfrog filter built from five µA741 amplifiers.

    Leapfrog (ladder-simulation) filters realise an LC ladder as a chain of
    integrators with **feedback between adjacent stages** -- that bidirectional
    coupling is the point of the topology, and it is what makes this a harder
    test than a cascade: the matrix does not decompose into independent blocks
    and every stage's response depends on its neighbours in both directions.

    Five inverting integrators simulate the ``L1 C2 L3 C4 L5`` prototype.  Each
    stage's summing junction takes a resistor from the previous stage, a
    resistor back from the *next* stage, and an integrating capacitor from its
    own output.

    Each amplifier is the full µA741 of :func:`ua741` -- 24 transistors -- so
    this is around 120 unknowns, an order of magnitude beyond the single
    amplifier and well beyond the cascade fixtures.

    Args:
        symbolic_devices: Device names whose ``gm`` stays symbolic.  Names are
            per-stage, ``'s0_q1'`` through ``'s4_q24'``.
        miller: Include each amplifier's compensation capacitor.

    Returns:
        BenchSystem.
    """
    import pycircuit.circuit.circuit as _circuit_module
    saved = _circuit_module.default_toolkit
    _circuit_module.default_toolkit = symbolic
    try:
        return _build_leapfrog(symbolic_devices, miller)
    finally:
        _circuit_module.default_toolkit = saved


def _build_leapfrog(symbolic_devices, miller):
    cir = SubCircuit(toolkit=symbolic)
    vin = sympy.Symbol('vin')
    params = {vin: 1.0}

    STAGES = 5
    names = ['in']
    node = {'in': cir.add_node('in')}
    cir['vs'] = VS(node['in'], gnd, vac=vin)

    ## One amplifier per reactive element of the LC prototype.  Nodes are
    ## created stage by stage in the µA741's own order, so each block of the
    ## matrix has the structure the single-amplifier fixture has.
    amp = []
    for k in range(STAGES):
        local = {}
        for base in _UA741_NODE_NAMES:
            name = 's%d_%s' % (k, base)
            local[base] = cir.add_node(name)
            names.append(name)
        amp.append(local)

    def bjt_for(stage):
        def bjt(dev, role, b, c, e):
            gm, rpi, ro, cpi, cmu = _UA741_ROLES[role]
            full = 's%d_%s' % (stage, dev)
            if full in symbolic_devices:
                sym = sympy.Symbol('gm_%s' % full, positive=True)
                params[sym] = gm
                gm = sym
            add_small_signal_bjt(cir, full, b, c, e, gm, rpi, ro, cpi, cmu)
        return bjt

    for k in range(STAGES):
        _stamp_ua741_devices(cir, amp[k], bjt_for(k), miller,
                             pfx='s%d_' % k)
        ## Non-inverting input grounded: these are inverting integrators.
        cir['sg%d' % k] = R(amp[k]['inp'], gnd, r=1.0)

    ## Integrating capacitor around each amplifier.
    for k in range(STAGES):
        cir['ci%d' % k] = C(amp[k]['out'], amp[k]['inn'], c=1e-9)

    ## The terminating stages are LOSSY integrators, damped by a resistor
    ## across the capacitor: they simulate the ladder's source and load
    ## resistances.  Without them the chain integrates at DC, the passband
    ## never flattens, and the result is not a filter at all -- which is
    ## exactly what the first version of this fixture did.
    cir['rd0'] = R(amp[0]['out'], amp[0]['inn'], r=10e3)
    cir['rd%d' % (STAGES - 1)] = R(amp[STAGES-1]['out'],
                                   amp[STAGES-1]['inn'], r=10e3)

    ## Forward path: input, then each stage into the next.
    cir['rin'] = R(node['in'], amp[0]['inn'], r=10e3)
    for k in range(1, STAGES):
        cir['rf%d' % k] = R(amp[k-1]['out'], amp[k]['inn'], r=10e3)

    ## THE LEAPFROG FEEDBACK: each stage but the last is also driven from the
    ## stage *after* it.  Without these the circuit is an integrator chain,
    ## not a ladder simulation.
    for k in range(STAGES - 1):
        cir['rb%d' % k] = R(amp[k+1]['out'], amp[k]['inn'], r=10e3)

    return system_from_circuit(cir, 'leapfrog_5th_order', params,
                               out_index=names.index('s%d_out'
                                                     % (STAGES - 1)),
                               in_index=names.index('in'))


def cascaded_opamps(blocks=1, symbolic_devices=()):
    """``blocks`` Miller-compensated gain stages in cascade.

    The structure Tan & Shi cascade for their comparison against SCAPP (TCAD
    2000, Fig. 13 and Table II), which is the one published measurement of a
    diagram against a *sequence of expressions* -- the same representation as
    :doc:`soe_symbolic`.  Their circuit grows by 16 unknowns per block; this one
    grows by a similar amount, so the interesting quantity, the ratio between the
    two representations, is comparable even though the blocks are not identical.

    Args:
        blocks: How many stages to cascade.
        symbolic_devices: Device names whose ``gm`` stays symbolic.

    Returns:
        BenchSystem.
    """
    import pycircuit.circuit.circuit as _circuit_module
    saved = _circuit_module.default_toolkit
    _circuit_module.default_toolkit = symbolic
    try:
        return _build_cascade(blocks, symbolic_devices)
    finally:
        _circuit_module.default_toolkit = saved


def _build_cascade(blocks, symbolic_devices):
    cir = SubCircuit(toolkit=symbolic)
    vin = sympy.Symbol('vin')
    params = {vin: 1.0}
    names = ['in']
    node = {'in': cir.add_node('in')}
    cir['vs'] = VS(node['in'], gnd, vac=vin)

    def gm_of(dev, nominal):
        if dev in symbolic_devices:
            sym = sympy.Symbol('gm_%s' % dev, positive=True)
            params[sym] = nominal
            return sym
        return nominal

    previous = node['in']
    for b in range(blocks):
        stage = ['e%d' % b, 'c%d' % b, 'd%d' % b, 'g%d' % b, 'o%d' % b]
        for name in stage:
            node[name] = cir.add_node(name)
            names.append(name)

        ## Differential pair against a mirror load, then a Miller-compensated
        ## common-emitter stage into a follower.
        add_small_signal_bjt(cir, 'q%da' % b, previous, node['c%d' % b],
                             node['e%d' % b], gm_of('q%da' % b, 2e-3),
                             25e3, 1e6, 3e-12, 0.5e-12)
        add_small_signal_bjt(cir, 'q%db' % b, gnd, node['d%d' % b],
                             node['e%d' % b], gm_of('q%db' % b, 2e-3),
                             25e3, 1e6, 3e-12, 0.5e-12)
        cir['rt%d' % b] = R(node['e%d' % b], gnd, r=500e3)
        add_small_signal_bjt(cir, 'q%dm' % b, node['c%d' % b], node['d%d' % b],
                             gnd, gm_of('q%dm' % b, 2e-3), 25e3, 1e6,
                             3e-12, 0.5e-12)
        cir['rc%d' % b] = R(node['c%d' % b], gnd, r=100e3)

        add_small_signal_bjt(cir, 'q%dg' % b, node['d%d' % b], node['g%d' % b],
                             gnd, gm_of('q%dg' % b, 5e-3), 10e3, 200e3,
                             10e-12, 1e-12)
        cir['rg%d' % b] = R(node['g%d' % b], gnd, r=100e3)
        cir['cc%d' % b] = C(node['d%d' % b], node['g%d' % b], c=5e-12)

        add_small_signal_bjt(cir, 'q%do' % b, node['g%d' % b], gnd,
                             node['o%d' % b], gm_of('q%do' % b, 20e-3),
                             5e3, 50e3, 20e-12, 2e-12)
        cir['ro%d' % b] = R(node['o%d' % b], gnd, r=10e3)
        previous = node['o%d' % b]

    return system_from_circuit(cir, 'cascaded_opamps_%d' % blocks, params,
                               out_index=names.index('o%d' % (blocks - 1)),
                               in_index=names.index('in'))


def cauer_lowpass(sections=2, symbolic=True):
    """An elliptic (Cauer) LC low-pass ladder -- the paper's worked example.

    Tan & Shi illustrate hierarchical analysis on a Cauer filter (TCAD 2000
    §VI).  It is also the one fixture here containing inductors, which matters:
    pycircuit gives an inductor a branch current rather than a ``1/(sL)``
    admittance, so the matrix entries stay *polynomial* in ``s`` and the
    s-expansion applies -- an admittance-form inductor would not expand.

    Each section is a parallel ``L``-``C`` trap in the series arm (which places
    the transmission zeros) and a shunt capacitor to ground.

    Args:
        sections: Number of ladder sections.
        symbolic: Keep component values symbolic.

    Returns:
        BenchSystem.
    """
    import pycircuit.circuit.circuit as _circuit_module
    saved = _circuit_module.default_toolkit
    _circuit_module.default_toolkit = _symbolic_toolkit_for(symbolic)
    try:
        return _build_cauer(sections, symbolic)
    finally:
        _circuit_module.default_toolkit = saved


def _symbolic_toolkit_for(is_symbolic):
    return symbolic


def _build_cauer(sections, use_symbols):
    from pycircuit.circuit import L as Inductor

    cir = SubCircuit(toolkit=symbolic)
    names = ['in'] + ['n%d' % i for i in range(sections)]
    node = {name: cir.add_node(name) for name in names}
    params = {}

    def value(prefix, i, nominal):
        if use_symbols:
            sym = sympy.Symbol('%s%d' % (prefix, i), positive=True)
            params[sym] = nominal
            return sym
        return nominal

    vin = sympy.Symbol('vin')
    params[vin] = 1.0
    cir['vs'] = VS(node['in'], gnd, vac=vin)
    cir['rs'] = R(node['in'], node['n0'], r=value('Rs', 0, 50.0))

    for i in range(sections):
        here = node['n%d' % i]
        there = node['n%d' % (i + 1)] if i + 1 < sections else gnd
        if i + 1 < sections:
            ## Parallel trap in the series arm places a transmission zero.
            cir['Ls%d' % i] = Inductor(here, there, L=value('Ls', i, 1e-3))
            cir['Cs%d' % i] = C(here, there, c=value('Cs', i, 1e-9))
        cir['Cp%d' % i] = C(here, gnd, c=value('Cp', i, 2e-9))
    cir['rl'] = R(node['n%d' % (sections - 1)], gnd, r=value('Rl', 0, 50.0))

    return system_from_circuit(cir, 'cauer_%d' % sections, params,
                               out_index=names.index('n%d' % (sections - 1)),
                               in_index=names.index('in'))
