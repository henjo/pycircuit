"""The stepping loop's state is declared (the C stepping loop plan's stage
R2a, 2026-10-10).

`transient._RUN_STATE` names every attribute a stepping loop writes on the
transient, by its role: CARRIED from one step to the next (what a C loop
marshals), decided once a RUN, or parked within one STEP.  These tests hold
it: every family's loop is observed -- `sys.monitoring` on
`_SteppingLoop.execute` and `accept`, nothing on the loop replaced, so the
paths are a plain run's -- and each attribute whose binding in the
transient's `__dict__` changes while it runs must be declared.  A declared
name must be one the package assigns.  The observer is proved alive on a
planted write.
"""
import glob
import os
import re
import sys

import pytest

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit import gnd, numeric
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import VS, C, R, VPulse
from pycircuit.circuit.integrator import (
    ESDIRK43Integrator,
    EulerIntegrator,
    Gear2Integrator,
    GLM2Integrator,
    GLM3Integrator,
    RadauIIA3Integrator,
    ThetaIntegrator,
    TrapezoidalIntegrator,
    TRBDF2Integrator,
)
from pycircuit.circuit.simwarnings import AccuracyWarning
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.transient import Transient

MON = sys.monitoring


def _pulsed_rc():
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=1e-7, tf=1e-7, pw=2e-5, per=5e-5)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c


def _mos_chain(n=4):
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['vg'] = VPulse('g0', gnd, v1=0.0, v2=1.8, td=2e-9, tr=1e-9, tf=1e-9, pw=5e-9, per=2e-8)
    for k in range(n):
        c[f'R{k}'] = R('vdd', f'd{k}', r=5e3)
        c[f'C{k}'] = C(f'd{k}', gnd, c=1e-14)
        c[f'M{k}'] = eh.MosLevel1Hdl(f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


def _rc(integ):
    return lambda: Transient(_pulsed_rc(), toolkit=numeric, integrator=integ(),
                             reltol=1e-4).solve(tend=1e-4, timestep=1e-6)


#: a run of every step family: each integrator adaptive on a pulsed RC (its
#: breakpoints), the coupled family, a fixed grid, and an HDL MOS chain (the
#: C Newton, the evaluation core, the stage solve in C)
RUNS = {
    **{f'rc-{i.__name__}': _rc(i) for i in (
        Gear2Integrator, EulerIntegrator, TrapezoidalIntegrator, ThetaIntegrator,
        RadauIIA3Integrator, TRBDF2Integrator, ESDIRK43Integrator, GLM2Integrator,
        GLM3Integrator)},
    'rc-coupled': lambda: Transient(_pulsed_rc(), toolkit=numeric, reltol=1e-4).solve(
        tend=1e-4, timestep=1e-6, coupled_lte=True),
    'rc-fixed': lambda: Transient(_pulsed_rc(), toolkit=numeric).solve(
        tend=1e-4, timestep=1e-6, fixed_timestep=True),
    'mos-fixed': lambda: Transient(_mos_chain(), toolkit=numeric).solve(
        tend=4e-8, timestep=2e-10, fixed_timestep=True),
    'mos-gear': lambda: Transient(_mos_chain(), toolkit=numeric, reltol=1e-4).solve(
        tend=4e-8, timestep=2e-10),
    'mos-radau': lambda: Transient(_mos_chain(), toolkit=numeric, reltol=1e-4,
                                   integrator=RadauIIA3Integrator()).solve(
        tend=4e-8, timestep=2e-10),
}


def _written(run):
    """The transient attributes whose binding changes while `run`'s stepping
    loops run (each `__dict__` entry's identity, compared at the loop's
    start and at every accepted step)."""
    tool = next((i for i in range(6) if MON.get_tool(i) is None), None)
    if tool is None:
        pytest.skip('no free sys.monitoring tool id')
    state = {'tr': None, 'base': None}
    written = set()
    execute, accept = trmod._SteppingLoop.execute.__code__, trmod._SteppingLoop.accept.__code__

    def start(code, offset):
        loop = sys._getframe(1).f_locals.get('self')
        if isinstance(loop, trmod._SteppingLoop) and state['tr'] is None:
            state['tr'] = loop.tr
            state['base'] = {k: id(v) for k, v in loop.tr.__dict__.items()}

    def on_return(code, offset, retval):
        loop = sys._getframe(1).f_locals.get('self')
        if state['tr'] is None or not isinstance(loop, trmod._SteppingLoop) \
                or loop.tr is not state['tr']:
            return
        now = {k: id(v) for k, v in state['tr'].__dict__.items()}
        base = state['base']
        written.update(k for k in set(now) | set(base) if now.get(k) != base.get(k))
        state['base'] = now
        if code is execute:
            state['tr'] = None

    MON.use_tool_id(tool, 'test_run_state')
    try:
        MON.register_callback(tool, MON.events.PY_START, start)
        MON.register_callback(tool, MON.events.PY_RETURN, on_return)
        MON.set_local_events(tool, execute, MON.events.PY_START | MON.events.PY_RETURN)
        MON.set_local_events(tool, accept, MON.events.PY_RETURN)
        ## (a force-accept on the coarse runs warns: what is written is
        ## observed here, not the accuracy)
        with quiet(AccuracyWarning):
            run()
    finally:
        MON.set_local_events(tool, execute, 0)
        MON.set_local_events(tool, accept, 0)
        MON.free_tool_id(tool)
    return written


@pytest.mark.parametrize('name', sorted(RUNS))
def test_every_attribute_a_loop_writes_is_declared(name):
    written = _written(RUNS[name])
    assert written, 'the observer saw no write: is it alive?'
    assert not written - set(trmod._RUN_STATE), sorted(written - set(trmod._RUN_STATE))


def test_each_declared_name_has_a_role_and_is_assigned_in_the_package():
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    texts = []
    for path in glob.glob(os.path.join(here, '*.py')) + glob.glob(os.path.join(here, 'shooting', '*.py')):
        with open(path, encoding='utf-8') as fh:
            texts.append(fh.read())
    src = '\n'.join(texts)
    for name, (role, meaning) in trmod._RUN_STATE.items():
        assert role in ('carried', 'run', 'step') and meaning, name
        assert re.search(r'\.' + name + r'\s*=[^=]|\[\'' + name + r'\'\]\s*=', src), name


def test_the_observer_sees_a_planted_write(monkeypatch):
    """Alive: a write the manifest does not declare, made inside a step, is
    reported."""
    roll = Transient._roll_history

    def planted(self, *a, **k):
        self._planted_state = object()
        return roll(self, *a, **k)
    monkeypatch.setattr(Transient, '_roll_history', planted)
    written = _written(RUNS['rc-Gear2Integrator'])
    assert '_planted_state' in written - set(trmod._RUN_STATE)
