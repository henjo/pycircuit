"""STAGE 12B -- what does flooring `Pulse.tr`/`tf` at MIN_EDGE cost?

SPICE's default rise and fall times are zero, which makes the edge a true
discontinuity with no derivative -- awkward for Fang's `p = df_ckt/dh`.  Flooring
them at the smallest step the analysis would ever take keeps the waveform a
function, so `du/dt` is finite (very large) and is the true slope of the ramp
being integrated.

THE CONCERN THIS MEASURES: flooring puts a breakpoint at `td + MIN_EDGE`, and a
breakpoint 1e-18 s after another one would force a 1e-18 s step -- twelve orders
below anything else in the run.

MEASURED: it does not, because `next_event`'s minimum-spacing guard absorbs a
breakpoint that close.  At MIN_EDGE=1e-18 the smallest accepted step is 1e-9,
the same as at MIN_EDGE=1e-9, and step count and wall time are within noise
(1036 steps / 1.10 s against 1047 / 1.03 s).  The waveform is identical to four
figures either way.

Run:
    MPLBACKEND=Agg PYTHONPATH=. python benchmarks/transient_review/stage12b_min_edge.py
"""
import warnings, time
import numpy as np
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.elements import R, C, VPulse
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.func import Pulse

def build():
    c = SubCircuit()
    c['vs'] = VPulse('a', gnd, v1=0.0, v2=1.0, td=1e-5, tr=0.0, tf=0.0,
                     pw=2e-5, per=5e-5)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-9)
    return c

print('%-12s %8s %8s %10s %12s %10s' %
      ('MIN_EDGE', 'steps', 'reject', 'min h', 'wall (s)', 'v(b) max'))
for me in (1e-18, 1e-15, 1e-12, 1e-9):
    Pulse.MIN_EDGE = me
    cir = build()
    tran = Transient(cir, toolkit=numeric, reltol=1e-5)
    hs = []
    from pycircuit.circuit import stepcontroller as sc
    orig = sc.IntegralController.evaluate_step
    def cap(self, *a, _o=orig, **k):
        acc, hn = _o(self, *a, **k)
        if acc: hs.append(k['h_curr'])
        return acc, hn
    sc.IntegralController.evaluate_step = cap
    t0 = time.perf_counter()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=1.2e-4, timestep=1e-6, coupled_lte=False)
    finally:
        sc.IntegralController.evaluate_step = orig
    wall = time.perf_counter() - t0
    st = tran.statistics
    vb = np.asarray(res.v('b').y, float).ravel()
    print('%-12.0e %8d %8d %10.2e %12.3f %10.4f' %
          (me, st.accepted_steps, st.rejected_steps,
           min(hs) if hs else float('nan'), wall, vb.max()))
Pulse.MIN_EDGE = 1e-18
