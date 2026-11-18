"""Steps around a breakpoint (`_SteppingLoop._corner_size`): LANDING AHEAD
(`Transient.land_ahead`) and THE CORNER'S START (`Transient.corner_cut`), in
the Python and in the C loop (loop.c `corner_size`)."""
import re
import types

import pytest

from pycircuit.circuit import _cbuild, _tran_loop_c, circuit
from pycircuit.circuit import transient as trmod
from pycircuit.circuit.elements import C, R, SubCircuit, VPulse, gnd
from pycircuit.circuit.tests.test_loop_c import _pulsed, _run, _same
from pycircuit.circuit.transient import Transient


def _size(t, h, t_break, on_corner=False, land_ahead=True, corner_cut=0.25, max_step=10.0,
          tend=100.0):
    lp = types.SimpleNamespace(
        tr=types.SimpleNamespace(par=types.SimpleNamespace(land_ahead=land_ahead,
                                                           corner_cut=corner_cut)),
        run=types.SimpleNamespace(max_step=max_step), tend=tend, landing_ahead=None)
    got = trmod._SteppingLoop._corner_size(lp, t, h, t_break, on_corner)
    return got, lp.landing_ahead


def test_landing_ahead_stretches_or_halves_the_last_stretch():
    """A breakpoint beyond the step but within two of it: within 1.1x and
    max_step a landing (the size left for `where` to set), else half the
    way; at or past two steps, at the step itself, past the run's end, or
    switched off: the size as asked."""
    assert _size(0.0, 1.0, 3.0) == (1.0, False)
    assert _size(0.0, 1.0, 2.0) == (1.0, False)
    assert _size(0.0, 1.0, 1.0) == (1.0, False)
    assert _size(0.0, 1.0, 1.05) == (1.0, True)
    assert _size(0.0, 1.0, 1.1) == (1.0, True)
    assert _size(0.0, 1.0, 1.05, max_step=1.0) == (0.525, False)
    assert _size(0.0, 1.0, 1.5) == (0.75, False)
    assert _size(0.0, 1.0, 1.9) == (0.95, False)
    assert _size(2.0, 1.0, 3.5, tend=3.5) == (1.0, False)
    assert _size(2.0, 1.0, 3.5, tend=3.6) == (0.75, False)
    assert _size(0.0, 1.0, 1.5, land_ahead=False) == (1.0, False)


def test_the_corner_cut_starts_small_from_a_landed_corner():
    """From the breakpoint just landed on: at most corner_cut x min(the size,
    the gap to the next); elsewhere, or with the cut 0, the size."""
    assert _size(0.0, 1.0, 3.0, on_corner=True) == (0.25, False)
    assert _size(0.0, 1.0, 0.4, on_corner=True) == (0.1, False)
    assert _size(0.0, 0.2, 3.0, on_corner=True, corner_cut=0.5) == (0.1, False)
    assert _size(0.0, 1.0, 3.0, on_corner=False) == (1.0, False)
    assert _size(0.0, 1.0, 3.0, on_corner=True, corner_cut=0.0) == (1.0, False)
    ## (the cut first, then the window on the cut size)
    assert _size(0.0, 1.0, 0.4, on_corner=True, corner_cut=0.3) == (0.12, False)


def test_the_stretch_is_one_value():
    """loop.c's `LAND_STRETCH` is `transient.LAND_STRETCH`."""
    m = re.search(r'^#define LAND_STRETCH (\S+)$', _cbuild.csrc('loop.c'), re.MULTILINE)
    assert m and float(m.group(1)) == trmod.LAND_STRETCH


def _trains():
    """An RC behind two pulse trains out of step: corners every few steps,
    landed on from every distance."""
    c = SubCircuit()
    c['v1'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1.3e-9, tr=2e-10, tf=3e-10, pw=1.7e-9, per=4.1e-9)
    c['r1'] = R('a', 'b', r=1e3)
    c['c1'] = C('b', gnd, c=1e-12)
    c['v2'] = VPulse('d', gnd, v1=0.0, v2=0.5, td=0.7e-9, tr=1e-10, tf=1e-10, pw=2.3e-9, per=5.3e-9)
    c['r2'] = R('d', 'b', r=2e3)
    return c


def _reached(build, monkeypatch, make=None, **kw):
    """The Python run's `_corner_size` outcomes: stretched, halved, cut (the
    C loops off: every attempt the Python's)."""
    from pycircuit.circuit import _tran_stageloop_c, _tran_tcloop_c
    for mod in (_tran_loop_c, _tran_stageloop_c, _tran_tcloop_c):
        monkeypatch.setattr(mod, 'ENABLED', False)
    seen = {'stretch': 0, 'split': 0, 'cut': 0}
    real = trmod._SteppingLoop._corner_size

    def counted(self, t, h, t_break, on_corner):
        got = real(self, t, h, t_break, on_corner)
        rem = float(t_break - t)
        cut = self.tr.par.corner_cut * min(h, rem)
        seen['stretch'] += self.landing_ahead
        seen['split'] += got == 0.5 * rem and not self.landing_ahead
        seen['cut'] += bool(on_corner and self.tr.par.corner_cut and cut < h)
        return got
    monkeypatch.setattr(trmod._SteppingLoop, '_corner_size', counted)
    tr = Transient(build(), toolkit=circuit.numeric, **(make or {}))
    tr.solve(fixed_timestep=False, **kw)
    monkeypatch.undo()
    return seen


@pytest.mark.parametrize('name, build, kw', [
    ('trains', _trains, {'tend': 4e-8, 'timestep': 1e-10}),
    ('pulse', _pulsed, {'tend': 6e-8, 'timestep': 1e-9}),
])
def test_the_c_loop_sizes_corners_as_the_python(name, build, kw, monkeypatch):
    """Each outcome reached in the Python run, and the C loop's run of the
    same circuit the Python's byte for byte (`test_loop_c._same`).  The
    linear trains: everything but the tallies -- the branch check's
    `branch_c:served` counts each C step there and the Python's once, with
    the corners' rules off as on."""
    seen = _reached(build, monkeypatch, **kw)
    assert seen['split'] > 0 and seen['cut'] > 0, seen
    if name != 'trains':
        _same(build, monkeypatch, steps_min=10, fixed_timestep=False, **kw)
        return
    assert seen['stretch'] > 0, seen
    a, ca = _run(build, True, monkeypatch, fixed_timestep=False, **kw)
    b, _cb = _run(build, False, monkeypatch, fixed_timestep=False, **kw)
    for k in ('x', 't', 'stats', 'warnings', 'state'):
        assert a[k] == b[k], k
    assert ca.get('loop_c:steps', 0) >= 100, ca


@pytest.mark.parametrize('loop', ['stage', 'transform'])
def test_the_radau_loops_size_corners_as_the_python(loop, monkeypatch):
    """Radau's C loops -- the stage loop, the transform loop -- at the
    corners' defaults (their own suites run with the cut off): the halving
    and the cut reached in the Python run, and each C loop's run the
    Python's (`test_stageloop_c._same`, `test_tcloop_c._same`)."""
    from pycircuit.circuit.integrator import RadauIIA3Integrator
    from pycircuit.circuit.tests import test_stageloop_c, test_tcloop_c
    mod = test_stageloop_c if loop == 'stage' else test_tcloop_c
    kw = {'tend': 4e-8, 'timestep': 2.5e-10}
    seen = _reached(_pulsed, monkeypatch, make={'integrator': RadauIIA3Integrator(),
                                                'radau_transform': loop == 'transform'}, **kw)
    assert seen['split'] > 0 and seen['cut'] > 0, seen
    mod._same(_pulsed, monkeypatch, steps_min=10, make={'corner_cut': 0.25}, **kw)


def test_no_step_before_a_corner_is_a_sliver(monkeypatch):
    """Landing ahead off, a step that lands on a corner is often under half
    the one before it (the remainder left as a step: 31 of the trains' 80
    landings); on, a tenth of that at most (one, after a rejected landing).
    The run's last step, onto its end, aside."""
    def slivers(**make):
        rec = []
        real = trmod._SteppingLoop.accept

        def accept(self):
            rec.append((self.h, bool(self.landing)))
            return real(self)
        monkeypatch.setattr(trmod._SteppingLoop, 'accept', accept)
        tr = Transient(_trains(), toolkit=circuit.numeric, corner_cut=0.0, **make)
        tr.solve(tend=4e-8, timestep=1e-10, fixed_timestep=False)
        monkeypatch.undo()
        landed = [i for i in range(1, len(rec) - 1) if rec[i][1]]
        return sum(1 for i in landed if rec[i][0] < 0.5 * rec[i - 1][0]), len(landed)
    on, n_on = slivers()
    off, n_off = slivers(land_ahead=False)
    assert n_on > 50 and n_off > 50, (n_on, n_off)
    assert off >= 20 and on <= off // 10, (on, off)
