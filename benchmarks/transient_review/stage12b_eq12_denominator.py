"""Is eq (14)'s denominator recoverable analytically?  Three ways, one truth.

d(f_lte)/dh is the TOTAL derivative of the normalised LTE with respect to the
step size, including the circuit's response.  Eq (12) forms it as
`q^T dxh + d`; the claim under test is that this is a catastrophic subtraction
and that `err * w'(h)/w(h)` gives the same number without one.

Ground truth is a central finite difference of `err` itself, with the circuit
RE-SOLVED TO CONVERGENCE at each perturbed h -- that is what "total" means.
"""
import warnings
import numpy as np
from pycircuit.circuit import numeric, gnd
from pycircuit.circuit.circuit import SubCircuit
from pycircuit.circuit.elements import R, C, VSin
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.stepcontroller import SolutionLTEController
from pycircuit.circuit._lte_kernels import extrapolation_error_weight


def build():
    c = SubCircuit()
    c['vs'] = VSin('a', gnd, va=1.0, freq=1e3)
    c['R'] = R('a', 'b', r=1e3)
    c['C'] = C('b', gnd, c=1e-7)
    return c, dict(tend=6e-5, timestep=1e-6)


## Run a short transient, then stop and probe the state it reached.
cir, kw = build()
tran = Transient(cir, toolkit=numeric, reltol=1e-5)
snap = {}
real = Transient.fang_timestep
def spy(self, x_prev, t_prev, h, x_hist, **k):
    out = real(self, x_prev, t_prev, h, x_hist, **k)
    if len(snap) == 0 and t_prev > 2e-5 and not k.get('hold_h'):
        snap.update(x_prev=np.array(x_prev), t_prev=t_prev, h=out[1],
                    x_hist=[np.array(z) for z in x_hist],
                    dt_last=getattr(self, '_dt_last', None),
                    dt_last2=getattr(self, '_dt_last2', None),
                    qlast=np.array(self._qlast), iqlast=np.array(self._iqlast))
    return out
Transient.fang_timestep = spy
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    tran.solve(coupled_lte=True, **kw)
Transient.fang_timestep = real

h0 = snap['h']; t_prev = snap['t_prev']; x_hist = snap['x_hist']
h_hist = [z for z in (snap['dt_last'], snap['dt_last2']) if z is not None]
print('probe point: t=%.6e  h=%.4e  history %s' %
      (t_prev, h0, ['%.3e' % z for z in h_hist]))

ctrl = SolutionLTEController(); ctrl.set_relref(tran.par.relref)
tran._qlast = snap['qlast']; tran._iqlast = snap['iqlast']
tran._dt_last = snap['dt_last']; tran._dt_last2 = snap['dt_last2']
tran._is_first_step = False; tran._no_history = False
tran.irefnode = cir.get_node_index(gnd)
tran.base_integrator = tran._get_integrator()

def converged_err(h):
    """Solve the circuit to convergence at this h, return normalised err."""
    tran._dt = h
    x = snap['x_prev'].copy()
    for _ in range(40):
        f, J = tran._residual_and_jacobian(x, t_prev + h)
        fr = np.delete(f, tran.irefnode)
        Jr = np.delete(np.delete(J, tran.irefnode, 0), tran.irefnode, 1)
        dx = np.insert(np.linalg.solve(Jr, -fr), tran.irefnode, 0.0)
        xn = cir.limit(x + dx, x, tran.epar)
        if np.max(np.abs(xn - x)) < 1e-14 * max(1.0, np.max(np.abs(xn))):
            x = xn; break
        x = xn
    etol = tran._lte_tolerance(ctrl, x, snap['x_prev'], h_hist)
    ok, err = tran._lte_in_band(ctrl, x, x_hist, h_hist, h, etol, 0.7, 3.0)
    return err, x, etol

err0, x0, etol0 = converged_err(h0)
eps = 1e-4 * h0
ep, _, _ = converged_err(h0 + eps)
em, _, _ = converged_err(h0 - eps)
truth = (ep - em) / (2 * eps)

## (1) analytic:  err * w'(h)/w(h)
acc, wp_over_w = 0.0, 1.0 / h0
for hh in h_hist:
    acc += hh
    wp_over_w += 1.0 / (h0 + acc)
analytic = err0 * wp_over_w

## (2) the subtraction eq (12) forms
tran._dt = h0
f, J = tran._residual_and_jacobian(x0, t_prev + h0)
Jr = np.delete(np.delete(J, tran.irefnode, 0), tran.irefnode, 1)
i_ctrl, q_val, d = ctrl.lte_gradients(x0, x_hist, h_hist, h0, etol0)
p = tran.residual_dh(x0, t_prev + h0, h0)
pr = np.delete(p, tran.irefnode)
dxh = np.insert(np.linalg.solve(Jr, -pr), tran.irefnode, 0.0)
qt_dxh = q_val * dxh[i_ctrl]
subtraction = qt_dxh + d

print()
print('err at h0                       : %+.6e' % err0)
print('TRUTH  d(err)/dh (re-solved FD) : %+.6e' % truth)
print('  (1) analytic err*w\'/w         : %+.6e   ratio %.3f' % (analytic, analytic/truth))
print('  (2) eq (12) q^T dxh + d       : %+.6e   ratio %.3f' % (subtraction, subtraction/truth))
print()
print('the two terms being subtracted  : q^T dxh = %+.6e' % qt_dxh)
print('                                        d = %+.6e' % d)
print('  cancellation: result is %.2e of the larger term' % abs(subtraction/qt_dxh))
