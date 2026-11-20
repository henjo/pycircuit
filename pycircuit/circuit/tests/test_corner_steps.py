"""Steps around a breakpoint (`_SteppingLoop._corner_size`): LANDING AHEAD
(`Transient.land_ahead`) and THE CORNER'S START (`Transient.corner_cut`), in
the Python and in the C loop (loop.c `corner_size`)."""
import collections
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


def _entries(mp):
    """THE C LOOP'S ENTRIES over a run, for the message where it served
    none: `_tran_loop_c.attach` and `.run` wrapped -- what `attach` gave the
    stepping loop, and at each entry which of `run`'s uncounted returns
    held (switched off, the entries spaced out, watching off, a kept
    decline, the C Newton's or the C error test's last call not served
    since the watch counter moved) and what the entry returned.  A decline
    is a `loop_c:` count; these returns are not, so a run one of them kept
    out says nothing else.  (The labels restate `run`'s gate; the run
    itself is the real one's.)"""
    seen = collections.Counter()
    L = _tran_loop_c
    attach, run = L.attach, L.run

    def attached(loop):
        got = attach(loop)
        seen[f'attach {type(loop.family).__name__}: '
             f'{"a state" if got is not None else "None"}'] += 1
        return got

    def entered(loop):
        now, c, tr = L._watch.now(), loop.c, loop.tr
        rec = tr.__dict__.get('_newton_c')
        ctrl = tr.step_controller
        lte = ctrl.__dict__.get('_lte_c') if ctrl is not None else None
        held = [why for why, on in (
            ('off', not L.ENABLED), ('spaced', bool(c.skip)), ('watching off', now == -1),
            ('kept', c.kept == now), ('newton not served', rec is None or rec.served != now),
            ('lte not served', lte is None or lte.served != now)) if on]
        got = run(loop)
        seen[f'entry {" + ".join(held) or "open"}: {bool(got)}'] += 1
        return got
    mp.setattr(L, 'attach', attached)
    mp.setattr(L, 'run', entered)
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
    with monkeypatch.context() as mp:
        entries = _entries(mp)
        a, ca = _run(build, True, monkeypatch, fixed_timestep=False, **kw)
        watching = _tran_loop_c._watch.STATUS
    b, _cb = _run(build, False, monkeypatch, fixed_timestep=False, **kw)
    for k in ('x', 't', 'stats', 'warnings', 'state'):
        assert a[k] == b[k], k
    ## (the message a string, whole: a dict's repr pytest cuts short, and
    ## the C loop kept out by an uncounted return names it only here)
    assert ca.get('loop_c:steps', 0) >= 100, (
        f"the C loop served {ca.get('loop_c:steps', 0)} steps; watching {watching}; "
        f"entries {dict(entries)}; counts {dict(sorted(ca.items()))}")
    ## (and served through the wrapped entry: the message reads the run's path)
    assert entries['entry open: True'] >= 1, dict(entries)


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
