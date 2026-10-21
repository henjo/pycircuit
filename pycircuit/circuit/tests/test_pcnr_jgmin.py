"""The junction-gmin term in PCNR's equations, held nodes, and the continued
solve's route.

`solve_dc(jgmin=g)` puts `g` across every junction pair
(`pcnr_junction_pairs`) IN the equations -- SPICE's gmin, stepped: true gmin
stepping -- in the Python loop and in C (`pcnrdc.c`'s `assemble`), its twin.
`solve_dc_continued` tries the plain solve, the gshunt ladder and its pure
solve, the junction-gmin ladder, then source steps -- a step that fails
retried by both ladders at its own source level, the ordinary chain's shape
(opampal reaches PCNR's answer only so).
"""
from types import SimpleNamespace

import numpy as np
import pytest

from pycircuit.circuit import gnd, pcnr
from pycircuit.circuit.elements import G as Gel
from pycircuit.circuit.tests.test_dc_pcnr import _bjt_mirror_20v
from pycircuit.circuit.tests.test_pcnr_c import _twin


@pytest.mark.parametrize('g', [1e-2, 1e-6])
def test_the_junction_term_in_c_is_the_pythons(g, monkeypatch):
    """Bit for bit, the C serving: x, the limited unknowns, the iterations."""
    c, _x0 = _bjt_mirror_20v()
    got, _counts = _twin(c, monkeypatch, jgmin=g)
    assert got[0] == 'ok', got


@pytest.mark.parametrize('kw', [{}, {'source': 0.5}, {'jgmin': 1e-4}],
                         ids=['full', 'half-sources', 'with-jgmin'])
def test_a_held_solve_in_c_is_the_pythons(kw, monkeypatch):
    """`held` -- the mirror's output held at 2 V, at `source` times that --
    the reduction a row map in C (`red`), the step and the KCL test of the
    held row out: bit for bit the Python's, the C serving."""
    c, _x0 = _bjt_mirror_20v()
    held = {c.get_node_index('no'): 2.0}
    got, _counts = _twin(c, monkeypatch, held=held, **kw)
    assert got[0] == 'ok', got
    x = np.frombuffer(got[1])
    assert x[c.get_node_index('no')] == kw.get('source', 1.0) * 2.0


def test_the_junction_term_is_a_conductance_across_each_junction():
    """The same circuit with `g` stamped as conductances across the junction
    pairs: its solution is the term's -- from it the term's solve converges
    at once and moves nothing beyond rounding."""
    g = 1e-3
    d, _ = _bjt_mirror_20v()
    pairs = pcnr.pcnr_junction_pairs(d)
    assert pairs
    for j, (_inst, _dev, ra, rb) in enumerate(pairs):
        d[f'gj{j}'] = Gel(d.nodes[ra], d.nodes[rb], g=g)
    x_el, v_el, _its = pcnr.solve_dc(d, gnd, reltol=1e-9)
    c, _ = _bjt_mirror_20v()
    x, _v, its = pcnr.solve_dc(c, gnd, x0=x_el, v_lim0=v_el, jgmin=g, reltol=1e-9)
    assert its == 1, its
    assert np.max(np.abs(x - x_el)) <= 1e-12 * np.max(np.abs(x_el)), np.max(np.abs(x - x_el))


class _FakeSolve:
    """`pcnr.solve_dc` standing in for a circuit: `ok(source, kind, calls)`
    decides each rung, every call recorded as `(source, kind, value)` --
    kind 'p' plain, 'g' a gshunt rung, 'j' a junction-gmin rung -- and its
    outcome in `done`."""

    def __init__(self, ok):
        self.ok, self.calls, self.done = ok, [], []

    def __call__(self, circuit, refnode, x0=None, v_lim0=None, gshunt=0.0, jgmin=0.0, source=1.0,
                 **_kw):
        kind = 'g' if gshunt else ('j' if jgmin else 'p')
        self.calls.append((source, kind, gshunt or jgmin))
        good = self.ok(source, kind, self)
        self.done.append(good)
        if not good:
            raise RuntimeError('fake: no convergence')
        return np.zeros(circuit.n), np.zeros(1), 1


def _route(monkeypatch, ok):
    fake = _FakeSolve(ok)
    monkeypatch.setattr(pcnr, 'solve_dc', fake)
    monkeypatch.setattr(pcnr, 'pcnr_junction_pairs', lambda _c: [('pair',)])
    ## (the routes alone: the fake's answers have no circuit to be checked
    ## against -- the acceptance check, `test_dc_pcnr`)
    monkeypatch.setattr(pcnr, 'dc_residual', lambda *a, **k: (0.0, 0))
    out = pcnr.solve_dc_continued(SimpleNamespace(n=3), gnd)
    return out[3], fake.calls


def _after_a_junction_rung(fake):
    """The previous call a junction-gmin rung that converged, at this call's
    source level: the ladder's pure solve."""
    c = fake.calls
    return len(c) > 1 and c[-2][1] == 'j' and fake.done[-1] and c[-2][0] == c[-1][0]


def test_the_junction_gmin_ladder_comes_after_the_gshunt_ladder(monkeypatch):
    """Plain, gshunt from its first rung, then the junction-gmin ladder by
    decades from 1e-2 to 1e-11 and its pure solve."""
    path, calls = _route(monkeypatch, lambda s, kind, fake: kind == 'j' or _after_a_junction_rung(fake))
    assert path == 'jgmin'
    assert [k for _s, k, _v in calls] == ['p', 'g'] + ['j'] * 10 + ['p'], calls
    assert [v for _s, k, v in calls if k == 'j'] == pytest.approx(
        [10.0 ** -e for e in range(2, 12)], rel=1e-9)


def test_a_failed_source_step_is_retried_by_both_ladders_at_its_level(monkeypatch):
    """Everything fails at full sources until the stepping starts; the first
    step (0.05) is reached through the junction-gmin ladder, which a failed
    step tries first, after its own plain rung, at 0.05; and where that fails
    too the gshunt ladder follows at the same level."""
    def ok(source, kind, fake):
        if abs(source - 0.05) < 1e-12:
            return kind == 'j' or _after_a_junction_rung(fake)
        return any(c[0] < 1.0 for c in fake.calls)
    path, calls = _route(monkeypatch, ok)
    assert path == 'source'
    at = [k for s, k, _v in calls if abs(s - 0.05) < 1e-12]
    assert at == ['p'] + ['j'] * 10 + ['p'], at
    assert calls[-1][0] == 1.0 and calls[-1][1] == 'p', calls[-3:]

    def gshunt_only(source, kind, fake):
        if abs(source - 0.05) < 1e-12:
            return kind == 'g' or (kind == 'p' and len(fake.calls) > 1 and fake.calls[-2][1] == 'g'
                                   and fake.done[-1])
        return any(c[0] < 1.0 for c in fake.calls)
    path, calls = _route(monkeypatch, gshunt_only)
    at = [k for s, k, _v in calls if abs(s - 0.05) < 1e-12]
    assert path == 'source' and at[:2] == ['p', 'j'] and 'g' in at, at


def test_the_plain_solve_needs_no_ladder(monkeypatch):
    path, calls = _route(monkeypatch, lambda s, kind, fake: True)
    assert path == 'plain' and len(calls) == 1
