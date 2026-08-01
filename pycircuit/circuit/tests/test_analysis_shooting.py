from pycircuit.circuit import *
from pycircuit.circuit.shooting import *
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest

class myC(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(np.zeros(2))
    array([[ 0.,  0.],
           [ 0.,  0.]])
    >>> c.C(np.zeros(2))
    array([[  1.0000e-12,  -1.0000e-12],
           [ -1.0000e-12,   1.0000e-12]])

    """

    terminals = ('plus', 'minus')
    instparams = [Parameter(name='c0', desc='Capacitance', 
                            unit='F', default=1e-12),
                  Parameter(name='c1', desc='Nonlinear capacitance', 
                            unit='F', default=0.5e-12),
                  Parameter(name='v0', desc='Voltage for nominal capacitance', 
                            unit='V', default=1),
                  Parameter(name='v1', desc='Slope voltage ...?', 
                            unit='V', default=1)]

    def C(self, x, epar=defaultepar): 
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        c = c0+c1*self.toolkit.tanh((v-v0)/v1)
        return self.toolkit.array([[c, -c],
                                  [-c, c]])

    def q(self, x, epar=defaultepar):
        v=x[0]-x[1]
        c0 = self.ipar.c0
        c1 = self.ipar.c1
        v0 = self.ipar.v0
        v1 = self.ipar.v1
        q = c0*v+c1*v1*self.toolkit.ln(self.toolkit.cosh((v-v0)/v1))
        return self.toolkit.array([q, -q])

def test_shooting():
    """PSS of a linear RC network must match the AC steady state.

    The shooting solver uses backward Euler, which is first order, so the
    error in the periodic steady state scales like O(1/N) in the number of
    timesteps per period.  N is chosen so the discretisation error is a few
    per-mille; the assertions below allow a comfortable margin on top of that
    rather than demanding near-exact equality (which no first-order method
    reaches at a practical N).
    """
    circuit.default_toolkit = circuit.numeric

    cir = SubCircuit()

    N = 500
    period = 1e-3

    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=1/period, phase=20)
    cir['R'] = R(1,2, r=1e4)
    cir['C'] = C(2,gnd, c=1e-8)

    resac = AC(cir).solve(1/period)

    pss = PSS(cir)

    res = pss.solve(period=period, timestep = period/N)

    v2ac = resac.v(2,gnd)
    v2pss = res['tpss'].v(2,gnd)

    t,dt = numeric.linspace(0,period,num=N,endpoint=True,
                            retstep=True)

    v2ref = numeric.imag(v2ac * numeric.exp(2j*numeric.pi*1/period*t))

    w2ref = Waveform(t,v2ref,ylabel='reference', yunit='V',
                     xunits=('s',), xlabels=('vref(2,gnd!)',))

    ## Check amplitude of the fundamental against the AC result
    v2rms_ac = np.abs(v2ac) / np.sqrt(2)
    v2rms_pss = np.abs(res['fpss'].v(2,gnd)).value(1/period)
    relerr = np.abs(v2rms_pss - v2rms_ac) / v2rms_ac
    assert relerr < 1e-2, 'amplitude rel. error=%g too high' % relerr

    ## Check error of the full waveform against the AC-reconstructed reference
    rmserror = np.sqrt(average((v2pss-w2ref)**2))
    assert rmserror < 1e-2, 'rmserror=%f too high' % rmserror

 
def test_PSS_nonlinear_C():
    """PSS with a tanh-nonlinear capacitor: the result must be PERIODIC.

    STAGE 10.2.  This test asserted NOTHING -- it called `pss.solve(...)` and
    checked no property of the result, so it passed whatever PSS returned,
    including nothing sensible.  It is the same class of defect as
    `test_sparse_toolkit` (passed while never exercising the sparse path) and
    `test_PAC` (skipped): a test that reads as coverage and is not.

    Periodicity is the right assertion because it is exactly what PSS claims to
    deliver -- x(0) = x(T) -- and it needs no external reference to check
    against, so it cannot be fitted. Measured at 2.2e-05 relative; asserted at
    1e-3, which is ~45x margin.

    The second assertion matters as much: without it the test would also pass on
    a degenerate solution that never leaves the capacitor's linear region, where
    `myC` is just a capacitor and the nonlinearity under test is not exercised.
    """
    circuit.default_toolkit = circuit.numeric
    c = SubCircuit()

    c['VSin'] = VSin(gnd, 1, va=10, freq=50e3)
    c['R1'] = R(1, 2, r=1e6)
    c['C'] = myC(2, gnd)
    #c['L'] = L(2,gnd, L=1e-3)
    pss = PSS(c)
    res = pss.solve(period=1/50e3,timestep=1/50e3/20)

    X = np.asarray(res['tpss'].x, dtype=float)
    assert np.isfinite(X).all(), 'PSS returned non-finite values'

    scale = np.abs(X).max()
    periodicity = np.abs(X[:, 0] - X[:, -1]).max() / scale
    assert periodicity < 1e-3, \
        'the steady state is not periodic: |x(0)-x(T)|/|x| = %.3e' % periodicity

    ## The nonlinearity must actually be traversed: myC's capacitance is
    ## c0 + c1*tanh((v - v0)/v1) with v0 = 1 V, so a solution confined near v0
    ## would exercise a plain capacitor and prove nothing about this circuit.
    v = X[c.get_node_index(2)]
    assert v.max() - v.min() > 2.0, \
        'v(C) spans only %.3f V -- the tanh knee is not being crossed' % (v.max() - v.min())


def test_PAC_is_withdrawn():
    """Stage 11: PAC says it is unimplemented instead of allocating 420 GiB.

    It was `@unittest.skip("Skip failing test")` -- advertised in the analysis
    inventory, never validated, and forming the whole (N*M)x(N*M) operator densely.
    An analysis that announces its absence is strictly better than one that fails
    somewhere deep in an allocation.
    """
    circuit.default_toolkit = circuit.numeric
    cir = SubCircuit()
    cir['vs'] = VSin(1, gnd, vac=2.0, va=2.0, freq=1e6, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2, gnd)
    cir['C'] = C(2, gnd, c=1e-12)
    pss = PSS(cir)
    res = pss.solve(period=1e-6, timestep=1e-6/10)
    with pytest.raises(NotImplementedError, match='withdrawn as unimplemented'):
        PAC(cir).solve(pss, np.array([1e6]))


@unittest.skip("Superseded by test_PAC_is_withdrawn; kept for the rewrite")
def test_PAC():
    circuit.default_toolkit = circuit.numeric
    N = 10
    fc = 1e6

    cir = SubCircuit()
    cir['vs'] = VSin(1,gnd, vac=2.0, va=2.0, freq=fc, phase=20)
    cir['R'] = R(1, 2, r=1e6)
    cir['D'] = Diode(2,gnd)
    cir['C'] = C(2,gnd, c=1e-12)
    
    pss = PSS(cir)

    res = pss.solve(period=1/fc, timestep = 1/(fc*N))

    pac = PAC(cir)
    res = pac.solve(pss, freqs = fc + np.array([1e3, 2e3, 4e3]))
    
    assert False, "Test should compare with spectre simulation"


## ---------------------------------------------------------------------------
## STAGE 11 -- PSS: `method` now selects something, and the inverse is a solve.
## ---------------------------------------------------------------------------

def _series_rlc(Lv=1e-3, Cv=1e-9, Rv=50.0, va=1.0):
    """Series RLC driven AT resonance, where |v(C)| = Q * va analytically."""
    import numpy as _np
    circuit.default_toolkit = circuit.numeric
    f0 = 1.0 / (2 * _np.pi * _np.sqrt(Lv * Cv))
    c = SubCircuit()
    c['vs'] = VSin(1, gnd, va=va, freq=f0)
    c['R'] = R(1, 2, r=Rv)
    c['L'] = L(2, 3, L=Lv)
    c['C'] = C(3, gnd, c=Cv)
    return c, f0, (1.0 / Rv) * _np.sqrt(Lv / Cv)


def _pss_peak(method, steps=20):
    import warnings
    cir, f0, Q = _series_rlc()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        res = PSS(cir, method=method).solve(period=1.0 / f0,
                                            timestep=1.0 / (f0 * steps))
    v = np.asarray(res['tpss'].v(3, gnd), dtype=float)
    return 0.5 * (v.max() - v.min()), Q


def test_pss_method_parameter_is_actually_read():
    """It was declared with default='euler' and never read anywhere in the file.

    A knob that advertises a choice it does not make is the "thin advertised
    feature" 0.1c warns about -- the same defect `lte_formula` was removed for in
    9(f).  Here the choice is worth having, so it is wired rather than deleted.
    """
    euler_peak, _Q = _pss_peak('euler')
    trap_peak, _Q = _pss_peak('trap')
    assert abs(trap_peak - euler_peak) > 0.1 * euler_peak, \
        'method= changes nothing: euler %.4f, trap %.4f' % (euler_peak, trap_peak)

    with pytest.raises(ValueError, match='method must be'):
        _pss_peak('bogus')


def test_backward_euler_damps_the_limit_cycle_and_trapezoidal_does_not():
    """The defect 0.1c names, measured against an ANALYTIC answer.

    Backward Euler's numerical damping attenuates exactly the limit cycle PSS
    exists to find.  On a Q=20 resonator driven at resonance, where the analytic
    peak is Q*va = 20 V:

        steps/period   euler          trapezoidal
                  20   2.63 V (13%)   19.32 V (97%)
                 200  12.20 V (61%)   19.23 V (96%)

    Euler is not merely less accurate, it is WRONG BY A FACTOR and gets worse as
    the step coarsens -- silently.
    """
    trap_peak, Q = _pss_peak('trap', steps=20)
    euler_peak, _ = _pss_peak('euler', steps=20)
    assert trap_peak > 0.9 * Q, \
        'trapezoidal recovers only %.1f%% of the analytic amplitude' % (100 * trap_peak / Q)
    assert euler_peak < 0.5 * Q, \
        'the euler damping this test documents has gone: %.1f%%' % (100 * euler_peak / Q)


def test_pss_uses_a_solve_not_an_explicit_inverse():
    """`inv(Jf) @ C @ Jshoot` formed a dense inverse per timestep per iteration.

    Asserted structurally: the quantity wanted is the solution of
    `Jf X = C @ Jshoot`, and at N=137/M=1000 with 20 shooting iterations the old
    form was 20,000 dense inversions.  A timing test would be a flake; the source
    check says exactly what changed.
    """
    import inspect
    from pycircuit.circuit import shooting
    src = inspect.getsource(shooting.PSS.solve)
    assert 'linalg.inv' not in src, 'the explicit inverse is back'
    assert 'linearsolver' in src


def test_pss_still_matches_the_ac_reference_with_a_fine_step():
    """Both methods must converge to the same, correct answer.

    At a coarse step neither is reliable and Euler's closeness is coincidence --
    measured, it is 0.9886 of the AC answer at dt = RC but 1.3283 at dt = RC/4.
    With a fine enough step both land on 1.0000, which is what makes the
    resonator comparison above a statement about damping rather than about luck.
    """
    import warnings
    from pycircuit.circuit.analysis_ss import AC
    circuit.default_toolkit = circuit.numeric

    def build():
        c = SubCircuit()
        c['VSin'] = VSin(gnd, 1, va=10, freq=50e3, vac=10)
        c['R1'] = R(1, 2, r=1e6)
        c['C'] = C(2, gnd, c=1e-12)
        c['L'] = L(2, gnd, L=1e-3)
        return c

    f = 50e3
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        ac = AC(build()).solve(freqs=np.array([f]))
    ref = abs(complex(np.asarray(ac.v(2, gnd)).ravel()[0]))

    for method in ('euler', 'trap'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = PSS(build(), method=method).solve(period=1 / f,
                                                    timestep=1 / f / 1280)
        v = np.asarray(res['tpss'].v(2, gnd), dtype=float)
        amp = 0.5 * (v.max() - v.min())
        assert amp == pytest.approx(ref, rel=0.02), \
            '%s gives %.6f against the AC reference %.6f' % (method, amp, ref)
