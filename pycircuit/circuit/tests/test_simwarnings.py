"""The warnings policy (the review's X8 / F13, 2026-10-01): what the analyses
warn with (`pycircuit.circuit.simwarnings`), where it is attributed, and
what is said once instead of per step."""
import warnings

import numpy as np
import pytest

from pycircuit.circuit import SubCircuit, circuit, gnd
from pycircuit.circuit.analysis import NoConvergenceError
from pycircuit.circuit.elements import C, R, VSin
from pycircuit.circuit.simwarnings import (
    AccuracyWarning,
    ConvergenceWarning,
    CostWarning,
    ModelWarning,
    PlatformWarning,
    SimulationWarning,
    UsageWarning,
    summarised,
    warn,
)
from pycircuit.circuit.tests._warnpolicy import quiet
from pycircuit.circuit.transient import Transient


def test_every_category_is_a_simulation_warning_and_a_runtime_warning():
    """A filter or `pytest.warns(RuntimeWarning)` written before the
    categories existed still sees every one of them."""
    for cat in (AccuracyWarning, ConvergenceWarning, ModelWarning,
                CostWarning, UsageWarning, PlatformWarning):
        assert issubclass(cat, SimulationWarning)
        assert issubclass(cat, RuntimeWarning)


def _rc():
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e6, vac=1.0)
    c['R'] = R('a', 'b', r=1e3, noisy=True)
    c['C'] = C('b', gnd, c=1e-10)
    return c


def test_a_warning_raised_deep_in_the_library_lands_on_the_callers_line():
    """`warn` attributes to the first frame outside the library, however deep
    it is raised.  A fixed `stacklevel` is right for one call chain only:
    pnoise's sideband-cap warning (`stacklevel=2`) landed inside the
    library when `band_spread` called it, and ~15 others landed there from
    every caller (until 2026-10-01)."""
    from pycircuit.circuit.shooting import PAC, PSS
    c = _rc()
    p = PSS(c, method='gear', reltol=1e-10)
    with quiet(AccuracyWarning):
        p.solve(period=1e-6, timestep=1e-6 / 40, maxiterations=20)
    pac = PAC(c)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        pac.pnoise(p, 0.13e6, 'b', maxsidebands=1, sweeptype='absolute')
    hits = [r for r in rec if 'maxsidebands=1' in str(r.message)]
    ## (outside the library: this file, or a wrapper of the call -- the
    ## gate's recorder wraps PAC's methods from outside it)
    from pycircuit.circuit.simwarnings import _LIBRARY
    assert hits and all(r.filename not in _LIBRARY for r in hits), \
        [(r.filename, r.lineno) for r in hits]
    assert any(r.filename == __file__ for r in hits) or all(
        'tran_recorder' in r.filename for r in hits), \
        [(r.filename, r.lineno) for r in hits]
    assert issubclass(hits[0].category, AccuracyWarning)


def test_a_fixed_grid_run_says_once_how_many_steps_it_could_not_take():
    """The fixed grid's Newton fallbacks and the force-accepts were one
    warning PER STEP, each with its own `t` -- which Python's
    once-per-location filter cannot collapse -- and the PCNR fallbacks a
    `logging.warning` per step.  Now one warning a run with the count and
    the first time (three injected failures: three warnings until
    2026-10-01)."""
    c = _rc()
    tran = Transient(c, reltol=1e-6)
    real = tran.solve_timestep
    fails = {'n': 0}

    def flaky(*args, **kwargs):
        ## (the first attempt at three grid points fails; their retries at a
        ## quarter step succeed)
        if fails['n'] < 3 and getattr(tran, '_dt', 0.0) > 5e-9:
            fails['n'] += 1
            raise NoConvergenceError('injected')
        return real(*args, **kwargs)
    tran.solve_timestep = flaky
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        tran.solve(refnode=gnd, tend=1e-6, timestep=2e-8, fixed_timestep=True)
    hits = [r for r in rec if 'no longer uniform' in str(r.message)]
    assert fails['n'] == 3
    assert len(hits) == 1, [str(r.message)[:80] for r in rec]
    assert '3 step(s)' in str(hits[0].message)
    assert issubclass(hits[0].category, ConvergenceWarning)


def test_pcnr_fallbacks_are_one_warning_a_run_not_a_log_line_a_step():
    """A PCNR failure falls back to device limiting for that step; the run
    said so in a `logging.warning` per failed step, outside the warnings
    machinery (no filter or `pytest.warns` could see it).  Now one
    `ConvergenceWarning` a run with the count and the first failure."""
    import logging

    from pycircuit.circuit.elements import Diode
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()
    c['V'] = VSin('a', gnd, va=2.0, freq=1e3)
    c['D'] = Diode('a', 'b')
    c['R'] = R('b', gnd, r=1e3)
    tran = Transient(c, pcnr=True)

    def boom(*a, **kw):
        raise RuntimeError('injected PCNR failure')
    tran._solve_timestep_pcnr = boom
    logged = []
    handler = logging.Handler()
    handler.emit = lambda r: logged.append(r.getMessage())
    logging.getLogger().addHandler(handler)
    try:
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            tran.solve(refnode=gnd, tend=2e-3, timestep=1e-5)
    finally:
        logging.getLogger().removeHandler(handler)
    assert tran.pcnr_fallbacks > 1
    hits = [r for r in rec if 'PCNR failed' in str(r.message)]
    assert len(hits) == 1, [str(r.message)[:80] for r in hits]
    assert issubclass(hits[0].category, ConvergenceWarning)
    assert f'on {tran.pcnr_fallbacks} solve(s)' in str(hits[0].message)
    assert 'injected PCNR failure' in str(hits[0].message)
    assert not any('PCNR' in m for m in logged), logged


def test_a_sub_solve_s_warnings_are_summarised_not_silenced():
    """`summarised` (the review's F13): the monodromy twin, `lte_grid`'s run
    and `warping_estimate`'s were run under a blanket 'ignore', hiding the
    accuracy of the very solve the answer was built on.  One warning a
    category that bears on the answer -- the count and the first -- and
    nothing for the rest."""
    def sub():
        warn('first accuracy note', AccuracyWarning)
        warn('second accuracy note', AccuracyWarning)
        warn('a cost note', CostWarning)
        warnings.warn('a foreign one', RuntimeWarning)
    with warnings.catch_warnings(record=True) as rec:
        warnings.simplefilter('always')
        with summarised('the sub-solve'):
            sub()
    msgs = [(r.category, str(r.message)) for r in rec]
    assert len(msgs) == 1, msgs
    cat, text = msgs[0]
    assert cat is AccuracyWarning
    assert text.startswith('the sub-solve warned 2 time(s) (AccuracyWarning)')
    assert 'first accuracy note' in text
    with pytest.raises(ValueError), summarised('raising'):
        raise ValueError('the caller handles it')


def test_lte_grid_passes_on_its_adaptive_run_s_accuracy():
    """`lte_grid` builds the grid from an adaptive run whose force-accepts
    (an LTE over tolerance, accepted anyway) were silenced; they reach the
    caller now, summarised."""
    from pycircuit.circuit.shooting import PSS
    c = _rc()
    p = PSS(c, method='gear', reltol=1e-6)
    tr_cls = type(p._new_transient(p._integrator_for('gear')))
    real = tr_cls._warn_run_summary

    def forced(self, fixed_fallbacks, forced_, timestep):
        return real(self, fixed_fallbacks, forced_ or [(1e-7, 1e-9, 4)],
                    timestep)
    tr_cls._warn_run_summary = forced
    try:
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter('always')
            p.lte_grid(1e-6, tstab=2e-6)
    finally:
        tr_cls._warn_run_summary = real
    hits = [r for r in rec if str(r.message).startswith(
        'PSS.lte_grid: its adaptive run warned')]
    assert len(hits) == 1, [str(r.message)[:80] for r in rec]
    assert issubclass(hits[0].category, AccuracyWarning)
    assert 'still above tolerance' in str(hits[0].message)
    assert np.isfinite(float(p.lte_period))


def test_no_library_module_warns_outside_the_policy():
    """Every warning the library gives goes through `simwarnings.warn` with a
    `SimulationWarning` category: no plain `warnings.warn` -- under any
    alias, or as `from warnings import warn` -- anywhere in the package but
    `simwarnings` itself.  The review's X8 moved the analyses' sites; six
    modules outside its scope kept fourteen (`hdl`, `_hdl_cache`,
    `_hdl_cbackend`, `ddd`, `analysis_ss`, `jaxtransient`) until
    2026-10-01: bare `RuntimeWarning`s (two the default `UserWarning`) that
    the suite's `error::SimulationWarning` never saw, attributed by fixed
    stacklevels."""
    import ast
    import os

    import pycircuit
    root = os.path.dirname(os.path.abspath(pycircuit.__file__))
    policy = os.path.join(root, 'circuit', 'simwarnings.py')
    found = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in ('tests', '__pycache__')]
        for name in filenames:
            path = os.path.join(dirpath, name)
            if not name.endswith('.py') or path == policy:
                continue
            with open(path) as fh:
                tree = ast.parse(fh.read())
            names = {'warnings'}
            for node in ast.walk(tree):
                if isinstance(node, ast.Import):
                    names |= {a.asname for a in node.names
                              if a.name == 'warnings' and a.asname}
                elif (isinstance(node, ast.ImportFrom)
                      and node.module == 'warnings'
                      and any(a.name == 'warn' for a in node.names)):
                    found.append(f'{os.path.relpath(path, root)}:'
                                 f'{node.lineno} (from warnings import warn)')
            for node in ast.walk(tree):
                if (isinstance(node, ast.Call)
                        and isinstance(node.func, ast.Attribute)
                        and node.func.attr == 'warn'
                        and isinstance(node.func.value, ast.Name)
                        and node.func.value.id in names):
                    found.append(f'{os.path.relpath(path, root)}:'
                                 f'{node.lineno}')
    assert not found, found
