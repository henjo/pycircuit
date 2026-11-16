"""Does the "narrow fix" of `doc/pcnr_native_design.md` section 5 survive PCNR?

The narrow fix: outside classic limiting, `Diode.G`/`Diode.i` evaluate at the
NODE voltage instead of the stale `_vlim`, so `cir.G(x)` becomes an honest
Jacobian for every caller.  The design study argued PCNR is unaffected because
`augmented_system` adds the diode through `cir.i`/`cir.G` and subtracts it again
through `element.i`/`element.G`, and "the cancellation is invariant to what G
returns provided both calls agree".

That is exact-arithmetic reasoning.  PCNR does not limit NODE voltages -- only
each device's own unknown -- so the node voltage the diode would be evaluated at
during an iteration is unbounded.  This script measures what the cancelled term
actually is, today and under the narrow fix.

What it simulates, and what it does not: during every `augmented_system` call
the diode is switched to node-voltage evaluation (NARROW) or left as shipped
(TODAY); both modes go through the same recording wrappers.  It does NOT
implement the flag for the classic limiting path, because that path is not in
question here -- only PCNR's add-then-subtract is.

Recorded per circuit and mode: PCNR assemblies, the largest junction node
voltage the diode is evaluated at inside them, the largest diode current that is
added and subtracted again, assemblies with a non-finite residual or Jacobian,
and the answer against TODAY.

Result on 2026-11-16, tree at 8048e42 (see pcnr_narrow_fix.log):

  1 A into a diode, PCNR DC: the first predict puts 1 A through the diode's
  off-conductance IS/VT = 3.9e-12 S, so the node lands at 2.58e11 V.  TODAY the
  stale-`_vlim` companion is affine in the node voltage and cancels 1 A;
  NARROW overflows `exp`, cancels inf, and 199 of 200 assemblies are non-finite.
  PCNR stops converging on exactly the overshoot it exists to absorb.

  Elsewhere the narrow fix still converges, with the same step counts and
  answers within 7.3e-12, but the cancelled term grows by 10-15 orders of
  magnitude (fig1: 3.9e-8 A -> 6.4e7 A, against ~1 mA of real current).  The
  answers agree because the rounding happened to be small on these circuits;
  nothing bounds it.

Run:
    PYTHONPATH=<repo> MPLBACKEND=Agg python3 -u benchmarks/transient_review/pcnr_narrow_fix.py
"""
import warnings
import numpy as np

from pycircuit.circuit import gnd, numeric
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VS, VSin, IS, Diode
from pycircuit.circuit import pcnr as P
from pycircuit.circuit.transient import Transient

STATE = {'pure': False, 'in_aug': False, 'vmax': -np.inf, 'imax': 0.0, 'nonfinite': 0, 'calls': 0}

_orig_G, _orig_i = Diode.G, Diode.i


def _pure_G(self, x, epar=None):
    if not STATE['pure']:
        return _orig_G(self, x, epar) if epar is not None else _orig_G(self, x)
    v = float(x[0] - x[1])
    VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
    with np.errstate(over='ignore'):
        g = self.iparv.IS * np.exp(v / VT) / VT
    return np.array([[g, -g], [-g, g]])


def _pure_i(self, x, epar=None):
    if not STATE['pure']:
        out = _orig_i(self, x, epar) if epar is not None else _orig_i(self, x)
    else:
        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        with np.errstate(over='ignore'):
            i = self.iparv.IS * (np.exp(float(x[0] - x[1]) / VT) - 1.0)
        out = np.array([i, -i])
    if STATE['in_aug']:
        STATE['vmax'] = max(STATE['vmax'], float(x[0] - x[1]))
        STATE['imax'] = max(STATE['imax'], float(abs(out[0])))
    return out


_orig_aug = P.augmented_system


def _aug(*a, **kw):
    STATE['in_aug'] = True
    STATE['pure'] = STATE['mode_pure']
    try:
        out = _orig_aug(*a, **kw)
    finally:
        STATE['pure'] = False
        STATE['in_aug'] = False
    STATE['calls'] += 1
    g_mna, _gl, J_mm = out[0], out[1], out[2]
    if not (np.all(np.isfinite(g_mna)) and np.all(np.isfinite(J_mm))):
        STATE['nonfinite'] += 1
    return out


def patch(narrow):
    """Both modes run through the recording wrappers; only `narrow` switches
    the diode to node-voltage evaluation inside `augmented_system`."""
    Diode.G, Diode.i, P.augmented_system = _pure_G, _pure_i, _aug
    STATE['mode_pure'] = narrow
    for k, v in (('pure', False), ('in_aug', False), ('vmax', -np.inf),
                 ('imax', 0.0), ('nonfinite', 0), ('calls', 0)):
        STATE[k] = v


def stats():
    return STATE['calls'], STATE['vmax'], STATE['imax'], STATE['nonfinite']


def one_amp():
    c = SubCircuit()
    c['is'] = IS(gnd, 1, i=1.0)
    c['D1'] = Diode(1, gnd)
    return c


def fig1():
    c = SubCircuit()
    c['vs'] = VS('e1', gnd, v=1.0)
    c['R'] = R('e1', 'e2', r=1e3)
    c['D1'] = Diode('e2', gnd, IS=1e-15)
    c['D2'] = Diode('e2', gnd, IS=1e-9)
    return c


def rect():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=5.0, freq=1e3)
    c['D'] = Diode('a', 'b', IS=1e-15)
    c['R'] = R('b', gnd, r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c


def mains():
    c = SubCircuit()
    c['VS'] = VSin(1, gnd, va=10.0, freq=50.0)
    c['D1'] = Diode(1, 2)
    c['C1'] = C(2, gnd, c=100e-6)
    c['R1'] = R(2, gnd, r=1e3)
    return c


def run_dc(builder):
    c = builder()
    x, v_lim, its = P.solve_dc(c, gnd)
    return np.asarray(x, float), its


def run_tran(builder, node, tend, coupled=False):
    tran = Transient(builder(), toolkit=numeric, reltol=1e-5, pcnr=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = tran.solve(tend=tend, timestep=1e-5 if tend < 1e-2 else 1e-4,
                         coupled_lte=coupled)
    y = np.asarray(res.v(node).y, float).ravel()
    return y, tran.statistics.accepted_steps


CASES = [
    ('dc   1 A into a diode', lambda: run_dc(one_amp)),
    ('dc   fig1 two diodes ', lambda: run_dc(fig1)),
    ('tran 5 V rectifier   ', lambda: run_tran(rect, 'b', 2e-3)),
    ('tran mains rectifier ', lambda: run_tran(mains, 2, 40e-3)),
    ('coup 5 V rectifier   ', lambda: run_tran(rect, 'b', 2e-3, coupled=True)),
]

if __name__ == '__main__':
    for name, fn in CASES:
        patch(False)
        try:
            ref, ref_n = fn()
            ref_err = None
        except Exception as e:  # noqa: BLE001
            ref, ref_n, ref_err = None, None, '%s: %s' % (type(e).__name__, str(e)[:80])
        r_calls, r_vmax, r_imax, r_nonfin = stats()
        patch(True)
        try:
            got, got_n = fn()
            err = None
        except Exception as e:  # noqa: BLE001
            got, got_n, err = None, None, '%s: %s' % (type(e).__name__, str(e)[:80])
        calls, vmax, imax, nonfin = stats()
        Diode.G, Diode.i, P.augmented_system = _orig_G, _orig_i, _orig_aug
        print('%s | TODAY  aug %5d | max v_node %10.4g V | max |i_D| cancelled %10.4g A | non-finite %4d'
              % (name, r_calls, r_vmax, r_imax, r_nonfin))
        line = '%s | NARROW aug %5d | max v_node %10.4g V | max |i_D| cancelled %10.4g A | non-finite %4d |' % (
            name, calls, vmax, imax, nonfin)
        if err or ref_err:
            line += ' ref=%s  narrow=%s' % (ref_err or 'ok(%s)' % ref_n, err or 'ok(%s)' % got_n)
        else:
            if ref.shape == got.shape:
                d = float(np.max(np.abs(ref - got)) / max(np.max(np.abs(ref)), 1e-30))
                line += ' n %s -> %s | max rel diff %.3e' % (ref_n, got_n, d)
            else:
                line += ' n %s -> %s | shapes differ %s vs %s' % (ref_n, got_n, ref.shape, got.shape)
        print(line, flush=True)
