"""A method set on an element's INSTANCE -- its `i`, `G`, `q`, `C`, `limit`,
`accept_step`, a source's `u` and `next_event` -- is what the Python calls,
so every C path and second path either calls it too or declines: each
analysis with every twin on is the run with every twin off, bit for bit,
the override planted either way.  (The PCNR operating point's C made a
participant's class kernels where its instance held its own `i` or `G`
until roadmap 203; the C loop's sources made a time function's class
waveform where it held its own `f` until 201.)"""
import importlib
import math
import warnings

import numpy as np
import pytest

from pycircuit.circuit import _paths, circuit, gnd
from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.elements import VS, C, R, VSin
from pycircuit.circuit.integrator import RadauIIA3Integrator
from pycircuit.circuit.shooting import PSS
from pycircuit.circuit.transient import Transient

#: every twin claimed bit-identical to the path it stands in for, by its
#: switch: off, the Python's generic path runs (the dict watcher, the
#: algorithmic choices -- the stage KLU, the predictor's fast path -- left as
#: they are)
OFF = (('_dc_newton_c', 'ENABLED'), ('_hdl_batch', 'ENABLED'), ('_hdl_batch', 'SKIP_ZERO_SOURCE'),
       ('_hdl_climit', 'ENABLED'), ('_hdl_climit', 'PSP_LIMIT'), ('_hdl_climit', 'WALK'),
       ('_pcnr_c', 'ENABLED'), ('shooting._lmmsens_c', 'LMMSENS_C'),
       ('shooting._sens_c', 'SENS_C'), ('shooting._sensklu_c', 'SENSKLU_C'),
       ('shooting._stagewalk_c', 'STAGEWALK_C'), ('_stamp_plan', 'ENABLED'),
       ('_stamp_plan', 'SOURCE_DIRECT'), ('_stamp_plan', 'SOURCE_PLAN'),
       ('_tran_branch_c', 'ENABLED'), ('_tran_core', 'CORE'), ('_tran_core', 'CORE_PASSES'),
       ('_tran_core', 'PATTERN'), ('_tran_core', 'SKIP_UNREAD_J'), ('_tran_ctrl_c', 'ENABLED'),
       ('_tran_loop_c', 'ENABLED'), ('_tran_lte_c', 'ENABLED'), ('_tran_newton_c', 'ENABLED'),
       ('_tran_pred_c', 'ENABLED'), ('_tran_radau_c', 'ENABLED'), ('_tran_radau', 'STAGE_FUSE'),
       ('_tran_radau_tc', 'ENABLED'), ('_tran_radau_tc', 'END'), ('_tran_radau_tc', 'FOLD'),
       ('_tran_sources', 'ENABLED'))

#: the period of the drives (1 MHz)
T = 1e-6


def _diode():
    """A charged SPICE diode (a PCNR participant, its limiter) on a 2 V sine."""
    c = SubCircuit()
    c['src'] = VSin('a', gnd, vo=0.5, va=2.0, freq=1e6)
    c['R'] = R('a', 'b', r=1e3)
    c['dev'] = eh.DiodeSpiceHdl('b', 'k', cjo=1e-12, tt=1e-9)
    c['RL'] = R('k', gnd, r=1e3)
    c['CL'] = C('k', gnd, c=1e-11)
    return c


def _mos():
    """Three gate-charge MOSFETs (a batch of a class) on a sine."""
    c = SubCircuit()
    c['vdd'] = VS('vdd', gnd, v=1.8)
    c['src'] = VSin('g0', gnd, v=0.9, va=0.8, freq=1e6)
    for k in range(3):
        c[f'rl{k}'] = R('vdd', f'd{k}', r=5e3)
        c['dev' if k == 1 else f'M{k}'] = eh.MosLevel1GateChargeHdl(
            f'd{k}', 'g0' if k == 0 else f'd{k - 1}', gnd, gnd)
    return c


CIRCUITS = {'diode': _diode, 'mos': _mos}
#: (the element, the method) planted: the device's evaluation and limiter,
#: its step hook; the source's value and events
PLANTS = (('dev', 'i'), ('dev', 'G'), ('dev', 'q'), ('dev', 'C'), ('dev', 'limit'),
          ('dev', 'accept_step'), ('src', 'u'), ('src', 'next_event'))


def _plant(c, name, meth, calls):
    """The method `meth` of element `name` set on its instance: the class's,
    its value moved (`i`, `G`, `q`, `C`, `u` by 1e-3), no limiting, every
    call counted, a breakpoint every period / 7 added."""
    el = c.elements[name]
    f = getattr(type(el), meth).__get__(el)
    if meth in ('i', 'G', 'q', 'C', 'u'):
        def w(*a, **k):
            return np.asarray(f(*a, **k), dtype=float) * 1.001
    elif meth == 'limit':
        def w(x, x0, *a, **k):
            return x
    elif meth == 'accept_step':
        def w(*a, **k):
            calls.append(1)
            return f(*a, **k)
    else:
        def w(t):
            nb = (math.floor(t / (T / 7)) + 1) * (T / 7)
            got = f(t)
            return nb if got is None else min(got, nb)
    setattr(el, meth, w)


def _answer(build, analysis, plant):
    """The analysis' answer (its bytes, or the error it raised) with `plant`
    -- `(element, method)` or None -- on a fresh circuit."""
    c, calls = build(), []
    if plant is not None:
        _plant(c, *plant, calls)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        try:
            if analysis == 'dc':
                out = [DC(c, toolkit=circuit.numeric).solve().x]
            elif analysis in ('gear', 'radau'):
                kw = {} if analysis == 'gear' else {'integrator': RadauIIA3Integrator()}
                tr = Transient(c, toolkit=circuit.numeric, **kw)
                r = tr.solve(tend=2 * T, timestep=T / 50)
                out = [r.x, r.sweep_values, [tr.statistics.accepted_steps]]
            else:
                p = PSS(c, method=analysis.split('_')[1])
                out = [p.solve(period=T, timestep=T / 40)['tpss'].x]
        except Exception as e:                                 # noqa: BLE001
            return 'raised ' + type(e).__name__ + ': ' + str(e)
    return [np.asarray(a, dtype=float).tobytes() for a in out] + [len(calls)]


@pytest.fixture(autouse=True)
def _numeric(monkeypatch):
    monkeypatch.setattr(circuit, 'default_toolkit', circuit.numeric)


def _differing(monkeypatch, build, analysis):
    """The plants (None: none) whose answer with every twin on is not the
    one with every twin off."""
    on = {p: _answer(build, analysis, p) for p in (None,) + PLANTS}
    with monkeypatch.context() as mp:
        for mod, attr in OFF:
            mp.setattr(importlib.import_module('pycircuit.circuit.' + mod), attr, False)
        off = {p: _answer(build, analysis, p) for p in (None,) + PLANTS}
    return [p for p in on if on[p] != off[p]]


@pytest.mark.parametrize('analysis', ['dc', 'gear', 'radau', 'pss_radau', 'pss_gear'])
@pytest.mark.parametrize('name', sorted(CIRCUITS))
def test_an_instance_override_is_every_paths(monkeypatch, name, analysis):
    assert _differing(monkeypatch, CIRCUITS[name], analysis) == []


def test_the_comparison_sees_a_twin_that_ignores_an_override(monkeypatch):
    """Alive: the PCNR operating point's C with its decline emptied makes the
    participant's class kernels, and the comparison sees its own `i` and `G`
    ignored -- and nothing else."""
    from pycircuit.circuit import _pcnr_c
    monkeypatch.setattr(_pcnr_c, '_OWN', ())
    assert _differing(monkeypatch, _diode, 'dc') == [('dev', 'i'), ('dev', 'G')]


def test_the_pcnr_operating_point_declines_a_participants_own_i():
    """A participant whose instance holds its own `i` is the Python PCNR's:
    the C counts the decline, and the operating point moves with it."""
    base = _answer(_diode, 'dc', None)
    before = _paths.snapshot()
    got = _answer(_diode, 'dc', ('dev', 'i'))
    assert _paths.since(before).get('pcnr_c:own', 0) > 0
    assert got != base
