"""The device bypass (`Transient(bypass=True)`; the evaluation core's,
core.h's bp_ fields): a compiled device whose first input, and each other
input's difference to it, moved within reltol * |v| + vabstol of its last
full evaluation is not evaluated but linearized there -- its currents and
charges from that evaluation's plus its G and C times the change, G and C
its own.  Off by default; armed for a Gear-2 run's stepping alone (the
method its accuracy was measured on).  SPICE's bypass, loaded as SPICE loads
a bypassed device (its companion at the reference)."""
import numpy as np
import pytest

from pycircuit.circuit import _paths, _tran_core, circuit
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit, gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, C, R, VPulse
from pycircuit.circuit.integrator import RadauIIA3Integrator, TrapezoidalIntegrator
from pycircuit.circuit.transient import Transient

#: a tolerance every move below sits well inside, and the move
LOOSE, STEP = (0.1, 1e-3), 1e-4


def _chain(n=6):
    """A pulsed chain of gate-charge MOS stages: idle between the edges."""
    c = SubCircuit()
    c.add_node('vdd')
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c.add_node('g0')
    c['vin'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-8, tr=1e-8, tf=1e-8, pw=2.4e-7, per=5e-7)
    for k in range(n):
        c.add_node(f'd{k}')
        c[f'rl{k}'] = R('vdd', f'd{k}', r=2e4)
        c[f'M{k}'] = eh.MosLevel1GateChargeHdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
        c[f'c{k}'] = C(f'd{k}', gnd, c=2e-14)
    return c


def _tr(integrator=None, **kw):
    if integrator is not None:
        kw['integrator'] = integrator
    return Transient(_chain(), toolkit=circuit.numeric, **kw)


def _run(tr, rt=None):
    if rt is not None:
        tr.par.reltol = rt
    return tr.solve(tend=1e-6, timestep=1e-8)


def _core():
    """A transient after a short run, its core, and a state of the run's."""
    tr = _tr()
    res = _run(tr)
    core = _tran_core.core_for(tr)
    if core is None:
        pytest.skip('the core does not serve this circuit here')
    x = np.array([np.asarray(res.v(nd).y, float)[-1] if nd is not gnd else 0.0
                  for nd in tr.cir.nodes] + [0.0] * len(tr.cir.branches))
    return tr, core, x


def _moved(tr, x, step=STEP, seed=0):
    """`x` with every node row but ground moved by up to `step`."""
    rng = np.random.default_rng(seed)
    y = x.copy()
    g = tr.cir.nodes.index(gnd)
    for r in range(len(tr.cir.nodes)):
        if r != g:
            y[r] += step * rng.uniform(-1.0, 1.0)
    return y


def test_off_by_default_and_disarmed_after_a_run():
    """`bypass` is off by default: nothing bypassed.  On, the run's stepping
    bypasses and the core is disarmed after it -- the operating point (the
    run's first point) the same as without -- and no `bypasstol` is set:
    the Python devices' cache needs it by name."""
    off = _tr()
    a = _run(off)
    assert off.par.bypass is False and off.statistics.bypassed_evaluations == 0
    on = _tr(bypass=True)
    b = _run(on)
    core = _tran_core.core_for(on)
    assert on.statistics.bypassed_evaluations > 0
    assert core.cs.bp_at == core.ffi.NULL
    assert on.epar.bypasstol == -1.0
    assert np.asarray(a.v('d5').y)[0] == np.asarray(b.v('d5').y)[0]
    assert Transient(_chain(), bypass=True, bypasstol=1e-5).epar.bypasstol == 1e-5
    assert DC(_chain(), bypass=True).epar.bypasstol == -1.0


def test_a_run_bypasses_on_its_c_path_within_its_tolerance():
    """Gear-2: devices bypassed, the run in C (the Newton) as without,
    every step count the same here, and each node's error against a
    reltol/100 reference at most 1.25x the error without bypass plus reltol
    (measured: 3.0e-2 both ways -- the edges' timing)."""
    runs, paths = {}, {}
    for name, on in (('off', False), ('on', True), ('ref', False)):
        tr = _tr(None, bypass=on)
        before = _paths.snapshot()
        runs[name] = (_run(tr, 1e-6 if name == 'ref' else None), tr.statistics)
        paths[name] = {k for k, v in _paths.since(before).items()
                       if k.endswith(':served') and k.split(':')[0] in ('newton_c', 'radau_c')}
    assert runs['on'][1].bypassed_evaluations > 0 and runs['off'][1].bypassed_evaluations == 0
    assert paths['on'] == paths['off'] and paths['on']
    assert runs['on'][1].accepted_steps == runs['off'][1].accepted_steps
    for k in range(6):
        r = runs['ref'][0].v(f'd{k}')
        tr_, yr = np.asarray(r.x[0], float), np.asarray(r.y, float)
        err = {}
        for name in ('off', 'on'):
            w = runs[name][0].v(f'd{k}')
            y = np.interp(tr_, np.asarray(w.x[0], float), np.asarray(w.y, float))
            err[name] = np.max(np.abs(y - yr)) / (yr.max() - yr.min())
        assert err['on'] <= 1.25 * err['off'] + 1e-4, (k, err)


@pytest.mark.parametrize('integrator', [RadauIIA3Integrator, TrapezoidalIntegrator],
                         ids=['radau', 'trap'])
def test_another_method_bypasses_nothing(integrator):
    """Another method than Gear-2 -- a stage method, another multistep
    formula -- evaluates every device: nothing bypassed, the run the same
    bytes as without."""
    out = []
    for on in (False, True):
        tr = _tr(integrator(), bypass=on)
        res = _run(tr)
        assert tr.statistics.bypassed_evaluations == 0
        out.append([np.asarray(res.v(f'd{k}').y).tobytes() for k in range(6)])
    assert out[0] == out[1]


def test_a_bypass_run_is_repeatable():
    """Two runs, and a second run on the same transient (every reference
    cleared at its arming), are the same bytes.  (A transient's first run
    makes two evaluations of its own -- bypass or not -- so the counts are
    compared between first runs.)"""
    tr = _tr(bypass=True)
    a = _run(tr)
    na = tr.statistics.bypassed_evaluations
    b = _run(tr)
    tc = _tr(bypass=True)
    c = _run(tc)
    for nd in ('d0', 'd3', 'd5'):
        ya = np.asarray(a.v(nd).y)
        assert ya.tobytes() == np.asarray(b.v(nd).y).tobytes() == np.asarray(c.v(nd).y).tobytes()
    assert na == tc.statistics.bypassed_evaluations > 0


def test_arming_clears_every_reference():
    """An arming starts with no reference: the first call after it
    evaluates every device, however close to the last arming's references,
    and the count starts at zero."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        _tran_core.passes(tr, x0, 'GCiq')
        _tran_core.passes(tr, _moved(tr, x0), 'GCiq')
        assert core.bypassed() == 6
        core.bypass(LOOSE)
        assert core.bypassed() == 0
        _tran_core.passes(tr, _moved(tr, x0, seed=4), 'GCiq')
        assert core.bypassed() == 0
    finally:
        core.bypass(None)


def test_a_bypassed_device_is_linear_at_its_reference():
    """Within the tolerance every device is bypassed: i and q are the
    reference's plus its G and C times the move -- not the reference's
    outputs (stale) nor a fresh evaluation -- and G and C the reference's."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        P0 = _tran_core.passes(tr, x0, 'GCiq')
        assert core.bypassed() == 0
        x1 = _moved(tr, x0)
        P1 = _tran_core.passes(tr, x1, 'GCiq')
        assert core.bypassed() == 6
    finally:
        core.bypass(None)
    exact = _tran_core.passes(tr, x1, 'GCiq')
    d = x1 - x0
    assert P1['G'].tobytes() == P0['G'].tobytes() and P1['C'].tobytes() == P0['C'].tobytes()
    np.testing.assert_allclose(P1['i'], P0['i'] + P0['G'] @ d, rtol=1e-12, atol=1e-20)
    np.testing.assert_allclose(P1['q'], P0['q'] + P0['C'] @ d, rtol=1e-12, atol=1e-30)
    assert not np.array_equal(P1['i'], P0['i']) and not np.array_equal(P1['i'], exact['i'])
    np.testing.assert_allclose(P1['i'], exact['i'], rtol=1e-3, atol=1e-9)


def test_beyond_the_tolerance_a_device_is_evaluated_and_is_the_reference():
    """A move past the tolerance on one gate: the devices it reaches are
    evaluated (their values a fresh evaluation's), the others bypassed with
    no move (their values the same); then the new point is their reference.
    A common-mode move of every input (the differences the same) is past the
    tolerance too: the first input's own value is tested."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        _tran_core.passes(tr, x0, 'GCiq')
        x2 = x0.copy()
        x2[tr.cir.get_node_index(tr.cir.get_node('d2'))] += 0.5
        P2 = _tran_core.passes(tr, x2, 'GCiq')
        n2 = core.bypassed()
        x3 = _moved(tr, x2, seed=1)
        P3 = _tran_core.passes(tr, x3, 'GCiq')
        n3 = core.bypassed()
        x4 = x3 + 0.5
        _tran_core.passes(tr, x4, 'GCiq')
        n4 = core.bypassed()
    finally:
        core.bypass(None)
    exact = _tran_core.passes(tr, x2, 'GCiq')
    ## (M2 and M3 read d2: four of the six bypassed, at no move)
    assert n2 == 4 and n3 == n2 + 6 and n4 == n3
    for m in 'GCiq':
        assert np.array_equal(P2[m], exact[m]), m
    np.testing.assert_allclose(P3['i'], P2['i'] + P2['G'] @ (x3 - x2), rtol=1e-12, atol=1e-20)


def test_a_partial_call_leaves_the_reference_whole():
    """A call that makes G and C elsewhere without being a reference (no i)
    leaves the staging another point's: a bypassed call after it has the
    reference's G and C, and its i linear by them."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        P0 = _tran_core.passes(tr, x0, 'GCiq')
        far = _moved(tr, x0, step=0.5, seed=2)
        Pf = _tran_core.passes(tr, far, 'GC')
        x1 = _moved(tr, x0, seed=3)
        P1 = _tran_core.passes(tr, x1, 'GCiq')
        assert core.bypassed() == 6
    finally:
        core.bypass(None)
    assert not np.array_equal(Pf['G'], P0['G'])
    assert P1['G'].tobytes() == P0['G'].tobytes() and P1['C'].tobytes() == P0['C'].tobytes()
    np.testing.assert_allclose(P1['i'], P0['i'] + P0['G'] @ (x1 - x0), rtol=1e-12, atol=1e-20)


def test_another_temperature_or_pack_evaluates():
    """Another temperature, or a device's pack another (its parameters
    packed again), and the device is evaluated, not bypassed by a reference
    of the old state."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        P0 = _tran_core.passes(tr, x0, 'GCiq')
        tr.epar.T = 310.0
        P1 = _tran_core.passes(tr, x0, 'GCiq')
        assert core.bypassed() == 0
        del tr.cir['M3'].__dict__['_hdl_cp']
        P2 = _tran_core.passes(tr, x0, 'GCiq')
        assert core.bypassed() == 5
    finally:
        core.bypass(None)
    exact = _tran_core.passes(tr, x0, 'GCiq')
    for m in 'GCiq':
        assert np.array_equal(P1[m], exact[m]) and np.array_equal(P2[m], exact[m]), m
    assert not np.array_equal(P0['i'], P1['i'])


def test_a_parameter_change_ends_the_bypass():
    """A device's parameter changed: the circuit's next evaluation is on
    another core, which the run did not arm -- every device evaluated, at
    the new parameters."""
    tr, core, x0 = _core()
    core.bypass(LOOSE)
    try:
        _tran_core.passes(tr, x0, 'GCiq')
        tr.cir['M3'].ipar.vto = 0.35
        P = _tran_core.passes(tr, _moved(tr, x0), 'GCiq')
        assert core.bypassed() == 0 and _tran_core.core_for(tr) is not core
    finally:
        core.bypass(None)
    exact = _tran_core.passes(tr, _moved(tr, x0), 'GCiq')
    for m in 'GCiq':
        assert P[m].tobytes() == exact[m].tobytes(), m


def test_a_group_with_another_input_never_bypasses():
    """A group whose inputs are not all node rows (a branch's current among
    them) never bypasses: its references' offset -1."""
    tr, core, x0 = _core()
    core.bp, core.nn = None, 3
    core.bypass(LOOSE)
    try:
        assert (core.bp[0][:core.nz] == -1).all()
        _tran_core.passes(tr, x0, 'GCiq')
        _tran_core.passes(tr, _moved(tr, x0), 'GCiq')
        assert core.bypassed() == 0
    finally:
        core.bypass(None)
        core.bp, core.nn = None, len(tr.cir.nodes)
    core.bypass(LOOSE)
    core.bypass(None)
    assert (core.bp[0][:core.nz] >= 0).all()
