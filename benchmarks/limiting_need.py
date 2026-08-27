"""Which models need limiting, and therefore PCNR -- measured, not argued.

Limiting exists because Newton steps into a region where the
linearisation is worthless.  That happens when a device's CONDUCTANCE
grows super-polynomially with a branch voltage the solver controls: the
step is computed from a `G` that is wrong by tens of decades by the time
it lands.  It does NOT happen merely because a device is nonlinear.

Two criteria were tried first and are wrong.  They are recorded because
each looks obviously right until it is measured:

* **Does the current overflow?**  `DiodeHdl` leaves float range at
  **18.35 V** -- the classic `709*V_T` signature.  But `DiodeSpiceHdl`
  and `GummelPoonNpnHdl` never overflow at any voltage, because the
  DSL's `expl` already bounds them, and they plainly still need
  limiting.  Overflow is a symptom of the DSL lacking a regulariser,
  not of the device needing a limiter.

* **How many decades does `|G|` span?**  This gives a FALSE POSITIVE.
  `ChargePumpHdl` spans **8.1 decades** over a +/-1 V sweep, more than
  `EkvNmosHdl` (6.8) or `PhotodiodeHdl` (4.7), both of which are
  eligible.  Its nonlinearity is a `tanh` with a narrow switching
  width, so the range comes from the TAILS going tiny -- its max `|G|`
  actually FALLS, from 1e-3 to 8e-12.  Dynamic range conflates "G
  becomes huge" with "G becomes small", and only the first is a problem.

What works is the GROWTH LAW of `max|G|` against branch voltage.  Run
this and the two populations do not overlap by eighty decades.

⚠ Two limits of the method, both real:

* A low reading is NOT proof a device is safe.  `EkvNmosHdl` reads
  bounded here because a probe from zero bias never reaches its
  exponential region, and it genuinely needs limiting.  The criterion
  is sound in one direction only: **large growth proves need**.  For the
  refused models the conclusion rests additionally on their
  nonlinearities being structurally bounded -- `tanh`, `maxc`/`minc`,
  products -- rather than merely un-probed.
* The sweep drives one row at a time from zero.  A device whose
  exponential needs two terminals moved together is under-probed.  Every
  row is swept, which is what catches an exponential living on an
  INTERNAL junction: probing only row 0 makes `GummelPoonNpnHdl` and
  `MosLevel1Hdl` look flat, and that nearly went into a report as data.

Usage::

    .venv/bin/python benchmarks/limiting_need.py
"""

import warnings

import numpy as np

from pycircuit.circuit import elements_hdl as eh
from pycircuit.circuit.circuit import Node, defaultepar
from pycircuit.circuit.hdl import Behavioural

#: Branch volts to probe.  The top of the ladder is deliberately far
#: outside any operating point: the question is what Newton sees when it
#: OVERSHOOTS, not what the device does in normal use.
VOLTS = (0.5, 1.0, 2.0, 5.0, 10.0)

#: `max|G|(10 V) / max|G|(1 V)` above which a device is taken to need
#: limiting.  Measured, the two populations sit at ~1e87 and ~1e4, so
#: anything in six orders of magnitude of slack would do; 1e12 is chosen
#: to be far above a plausible high-order polynomial (a 10th-order term
#: gives 1e10 across that decade) and far below the exponentials.
NEEDS_LIMITING = 1e12


def library_classes():
    """Every compiled `Behavioural` in the library, by name."""
    return sorted(n for n, c in vars(eh).items()
                  if isinstance(c, type) and issubclass(c, Behavioural)
                  and c is not Behavioural and '_hdl_info' in c.__dict__)


def growth(name, volts=VOLTS):
    """``(max|G| per volt, ratio)`` for `name`, sweeping every row.

    Returns ``(None, None)`` when the device cannot be evaluated.
    """
    cls = getattr(eh, name)
    try:
        el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))])
        el.update_iparv()
    except Exception:
        return None, None
    out = []
    for v in volts:
        best = 0.0
        for row in range(el.n):
            x = np.zeros(el.n)
            x[row] = v
            try:
                g = np.asarray(el.G(x, defaultepar), dtype=float)
            except Exception:
                continue
            m = float(np.max(np.abs(g)))
            if np.isfinite(m):
                best = max(best, m)
        out.append(best)
    if len(out) < 2 or out[1] <= 0.0:
        return out, None
    return out, out[-1] / out[1]


def eligible(name):
    """Does `name` participate in PCNR at its default card?"""
    cls = getattr(eh, name)
    try:
        el = cls(*[Node('n%d' % i) for i in range(len(cls.terminals))])
        el.update_iparv()
    except Exception:
        return None
    return bool(getattr(el, 'pcnr_probes', None)
                or getattr(el, 'pcnr_junctions', ()))


def main():
    warnings.simplefilter('ignore')
    rows = []
    for name in library_classes():
        g, ratio = growth(name)
        rows.append((name, g, ratio, eligible(name)))

    print('max|G| over every row, at branch V = %s\n'
          % ', '.join('%g' % v for v in VOLTS))
    print('%-24s %-11s %-8s %s' % ('model', 'G(10)/G(1)', 'PCNR?', 'verdict'))
    print('-' * 74)
    disagree = []
    for name, g, ratio, ok in sorted(rows, key=lambda r: -(r[2] or 0.0)):
        if ratio is None:
            ## A device with no resistive current at all -- CHdl is a
            ## pure capacitor.  `max|G|` is identically zero, so there is
            ## no ratio to take, and "needs limiting" is not a question
            ## that applies to it.
            print('%-24s %-11s %-8s %s'
                  % (name, 'n/a', 'yes' if ok else 'no',
                     'no resistive current'))
            continue
        needs = ratio > NEEDS_LIMITING
        verdict = 'NEEDS LIMITING' if needs else 'bounded / polynomial'
        flag = ''
        if needs and not ok:
            flag = '   <-- GROWS AND IS REFUSED'
            disagree.append(name)
        print('%-24s %-11.2e %-8s %s%s'
              % (name, ratio, 'yes' if ok else 'no', verdict, flag))

    print()
    if disagree:
        print('⚠ %d model(s) grow super-polynomially and PCNR refuses them: %s'
              % (len(disagree), ', '.join(disagree)))
    else:
        print('No model grows super-polynomially without being PCNR-eligible.')
    print()
    print('⚠ The converse does NOT hold and this table cannot show it: a low')
    print('   ratio is not proof a device is safe.  EkvNmosHdl reads bounded')
    print('   because a probe from zero bias never reaches its exponential')
    print('   region, and it does need limiting.  See the module docstring.')


if __name__ == '__main__':
    main()
