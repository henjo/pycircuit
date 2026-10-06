"""The constant-stamp plan (`pycircuit/circuit/_stamp_plan.py`, the speed
plan's P2, 2026-10-01): the assembly that does not re-stamp elements whose
stamps cannot change must give the legacy loop's matrices and vectors BIT FOR
BIT, rebuild when a parameter or the topology changes, and stay out of the
paths it was not built for."""
import contextlib
import warnings

import numpy as np

from pycircuit.circuit import _stamp_plan, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.elements import (
    CCCS,
    CCVS,
    IS,
    VCCS,
    VCVS,
    BSource,
    C,
    Diode,
    Gyrator,
    L,
    R,
    SubCircuit,
    Transformer,
    VSin,
    gnd,
)


@contextlib.contextmanager
def plan(on):
    was = _stamp_plan.ENABLED
    _stamp_plan.ENABLED = on
    try:
        yield
    finally:
        _stamp_plan.ENABLED = was


def mixed():
    """Every kind the plan sorts: constant G (R, controlled sources,
    gyrator, transformer, the source's branch), constant C (C, L), zero
    stamps (IS), non-constant hand-written elements (Diode, BSource), hdl
    elements and a nested SubCircuit."""
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['r1'] = R('a', 'b', r=1e3)
    c['c1'] = C('b', gnd, c=1e-9)
    c['l1'] = L('b', 'c', L=1e-6)
    c['r2'] = R('c', gnd, r=50.0)
    c['is'] = IS('c', gnd, i=1e-3)
    c['g1'] = VCCS('a', gnd, 'd', gnd, gm=1e-3)
    c['r3'] = R('d', gnd, r=2e3)
    c['e1'] = VCVS('b', gnd, 'e', gnd, g=2.0)
    c['r4'] = R('e', 'f', r=1e3)
    c['h1'] = CCVS('f', gnd, 'g', gnd, r=10.0)
    c['f1'] = CCCS('a', gnd, 'h', gnd, F=0.5)
    c['r5'] = R('h', gnd, r=300.0)
    c['gy'] = Gyrator('b', gnd, 'k', gnd)
    c['tr'] = Transformer('k', gnd, 'm', gnd, n=2.0)
    c['r6'] = R('m', gnd, r=75.0)
    c['d1'] = Diode('c', gnd)
    c['bs'] = BSource('d', gnd, gnd, 'd', i_func=lambda u: 1e-3 * u ** 3)
    c['rh'] = eh.RHdl('e', gnd, r=1e3)
    c['dh'] = eh.DiodeHdl('h', gnd)
    c['sub'] = _Sub('f')
    return c


class _Sub(SubCircuit):
    """A nested SubCircuit: assembled as one element of its parent, with
    a plan of its own."""
    terminals = ('p',)

    def __init__(self, *args, **kvargs):
        super().__init__(*args, **kvargs)
        self['rs'] = R(self.nodenames['p'], gnd, r=1e3)
        self['cs'] = C(self.nodenames['p'], gnd, c=1e-12)


def states(n, k=200, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for j in range(k):
        x = rng.standard_normal(n) * 10.0 ** rng.integers(-12, 4, n)
        x[rng.random(n) < 0.15] = 0.0
        x[rng.random(n) < 0.05] = -0.0
        out.append(x)
    return out


def test_the_plan_assembles_bit_for_bit_what_the_element_loop_does():
    cir = mixed()
    n = cir.n
    for x in states(n):
        for m in ('C', 'q', 'i', 'G', 'G', 'i'):
            ## (the extreme states overflow the diodes on both paths alike)
            with plan(True), np.errstate(all='ignore'):
                a = getattr(cir, m)(x)
            with plan(False), np.errstate(all='ignore'):
                b = getattr(cir, m)(x)
            assert a.dtype == b.dtype and a.shape == b.shape, m
            assert a.tobytes() == b.tobytes(), (m, np.max(np.abs(a - b)))
    assert cir.__dict__['_stamp_plan'].builds == 1


def test_an_assembly_executes_no_import_statement(monkeypatch):
    """Speed round 12: the plan's checks and the assembly's passes ran a
    function-level import on every call (`_plan_for`, `_ineligible`, the
    submatrix and subvector passes: 0.8-1.7 k instructions each, ~2 % of a
    small circuit's PSS).  Assemblies through the plan now execute none --
    every element kind, the nested circuit's too (the parent: four or more
    each pass)."""
    import builtins
    cir = mixed()
    x = states(cir.n, k=1, seed=3)[0]
    with np.errstate(all='ignore'):
        for m in ('G', 'C', 'i', 'q'):
            getattr(cir, m)(x)
    seen = []
    real = builtins.__import__

    def counting(name, *a, **k):
        seen.append(name)
        return real(name, *a, **k)
    monkeypatch.setattr(builtins, '__import__', counting)
    with np.errstate(all='ignore'):
        for _ in range(5):
            for m in ('G', 'C', 'i', 'q'):
                getattr(cir, m)(x)
    monkeypatch.undo()
    assert seen == [], sorted(set(seen))


def test_a_pass_with_nothing_to_stamp_is_the_loops_float_zeros():
    """A capacitor's G, a resistor's C: every entry an exact zero, so the
    plan's bincount is EMPTY -- and an empty bincount is int64 even with
    weights (the gate's doctest caught it, 2026-10-01)."""
    circuit.default_toolkit = circuit.numeric
    for el, ms in ((C('a', gnd, c=1e-12), ('G', 'i')),
                   (R('a', gnd, r=1e3), ('C', 'q'))):
        cir = SubCircuit()
        cir['e'] = el
        for x in states(cir.n, 5):
            for m in ms:
                with plan(True):
                    a = getattr(cir, m)(x)
                with plan(False):
                    b = getattr(cir, m)(x)
                assert a.dtype == b.dtype == np.float64, (m, a.dtype)
                assert a.tobytes() == b.tobytes(), m


def test_a_non_finite_state_takes_the_element_loop():
    cir = mixed()
    x = states(cir.n, 1)[0]
    x[1] = np.inf
    x[2] = np.nan
    for m in ('G', 'C', 'i', 'q'):
        with plan(True), np.errstate(all='ignore'):
            a = getattr(cir, m)(x)
        with plan(False), np.errstate(all='ignore'):
            b = getattr(cir, m)(x)
        assert a.tobytes() == b.tobytes(), m


def test_a_parameter_or_topology_change_rebuilds_the_plan():
    cir = mixed()
    x = states(cir.n, 1)[0]
    with plan(True):
        g0 = cir.G(x)
    builds = cir.__dict__['_stamp_plan'].builds
    cir['r1'].ipar.r = 2e3
    cir.update_iparv()
    with plan(True):
        g1 = cir.G(x)
    with plan(False):
        g1_ref = cir.G(x)
    assert g1.tobytes() == g1_ref.tobytes()
    assert not np.array_equal(g0, g1), 'the new resistance was not seen'
    assert cir.__dict__['_stamp_plan'].builds > builds
    ## topology: an element added after a pass
    cir['r7'] = R('a', 'm', r=10.0)
    x = states(cir.n, 1, seed=3)[0]
    with plan(True):
        a = cir.i(x)
    with plan(False):
        b = cir.i(x)
    assert a.tobytes() == b.tobytes()


def test_which_elements_the_plan_treats_as_constant():
    from pycircuit.circuit.circuit import IProbe
    kinds = {cls.__name__: (_stamp_plan.pair_kind(cls, 'G'),
                            _stamp_plan.pair_kind(cls, 'C'))
             for cls in (R, C, L, VSin, IS, VCCS, VCVS, CCVS, CCCS, Gyrator,
                         Transformer, Diode, BSource, eh.RHdl, eh.DiodeHdl,
                         IProbe, SubCircuit)}
    assert kinds == {
        'R': ('cached', 'zero'), 'C': ('zero', 'cached'),
        'L': ('cached', 'cached'), 'VSin': ('cached', 'zero'),
        'IS': ('zero', 'zero'), 'VCCS': ('cached', 'zero'),
        'VCVS': ('cached', 'zero'), 'CCVS': ('cached', 'zero'),
        'CCCS': ('cached', 'zero'), 'Gyrator': ('cached', 'zero'),
        'Transformer': ('cached', 'zero'),
        ## computed per state, or generated i/q -- never constant
        'Diode': (None, 'zero'), 'BSource': (None, None),
        'RHdl': (None, None), 'DiodeHdl': (None, None),
        'IProbe': ('cached', 'zero'), 'SubCircuit': (None, None)}

    ## an override drops out by method identity, whatever it inherits
    class _R2(R):
        def G(self, x, epar=None):
            return R.G(self, x)
    assert _stamp_plan.pair_kind(_R2, 'G') is None


def test_the_plan_stays_out_of_the_paths_it_was_not_built_for():
    cir = mixed()
    x = states(cir.n, 1)[0]
    with plan(True):
        ## a complex state, a parameter tree, a dtype on a vector pass
        assert _stamp_plan.assemble_matrix(cir, 'G', x.astype(complex),
                                           (None,)) is None
        assert _stamp_plan.assemble_vector(cir, 'i', x[:-1], (None,)) is None
    ## (and the circuit's own answer through that path is the loop's -- whose
    ## real-pinned `matrix_from_entries` drops the imaginary part with a
    ## ComplexWarning, plan or no plan)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', np.exceptions.ComplexWarning)
        with plan(True):
            a = cir.G(x.astype(complex))
        with plan(False):
            b = cir.G(x.astype(complex))
    assert a.tobytes() == b.tobytes()


def test_a_transient_builds_the_plan_once_and_steps_bit_for_bit():
    from pycircuit.circuit.transient import Transient
    circuit.default_toolkit = circuit.numeric

    def ladder():
        c = SubCircuit()
        c['vs'] = VSin('n0', gnd, va=2.0, freq=1e3)
        for k in range(8):
            c[f'R{k}'] = R(f'n{k}', f'n{k + 1}', r=1e3)
            c[f'C{k}'] = C(f'n{k + 1}', gnd, c=1e-8)
        c['D'] = Diode('n8', gnd)
        return c
    got = {}
    for on in (True, False):
        cir = ladder()
        with plan(on):
            res = Transient(cir, reltol=1e-5).solve(tend=2e-3, timestep=1e-5)
        got[on] = np.asarray(res.x, dtype=float)
        if on:
            assert cir.__dict__['_stamp_plan'].builds == 1
    assert got[True].tobytes() == got[False].tobytes()
