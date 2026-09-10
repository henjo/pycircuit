import os, sys, warnings, numpy as np; warnings.simplefilter('ignore')
sys.path.insert(0, '/home/andreas/source/pycircuit/pycircuit/circuit/tests')
src = open(os.path.dirname(os.path.abspath(__file__)) + '/cyclo_gate.py').read().split('cA = build')[0]
exec(src)
from pycircuit.circuit.hdl import flicker_noise
import test_analysis_shooting as TAS
class ModFlicker(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]
    @staticmethod
    def analog(p, outp, outn, b, bn):
        return Contribution(Branch(outp, outn).I, flicker_noise((p.k * Branch(b, bn).V) ** 2, 1))
def build_f(kind):
    c = SubCircuit()
    for n in ('lo', 'mid', 'out'):
        c.add_node(n)
    c['vlo'] = VSin('lo', gnd, va=1.0, vo=0.3, freq=f0)
    if kind == 'A':
        c.add_node('n'); c['xi'] = TAS._Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0); c['Rn'] = R('n', gnd, r=Rn)
        c['M1'] = Mult('mid', gnd, 'n', gnd, 'lo', gnd, k=k1)
    else:
        c['src'] = ModFlicker('mid', gnd, 'lo', gnd, k=k1 * Rn)
    c['Rm'] = R('mid', gnd, r=Rm); c['M2'] = Mult('out', gnd, 'mid', gnd, 'lo', gnd, k=k2)
    c['Ro'] = R('out', gnd, r=Ro); c['Co'] = C('out', gnd, c=Co)
    return c
cA = build('A'); pA, pacA = solve(cA); oA = [str(n) for n in cA.nodes].index('out')
cB = build('B'); pB, pacB = solve(cB); oB = [str(n) for n in cB.nodes].index('out')
for f in (0.13 * f0, 1.37 * f0):
    sA, _ = pacA.pnoise(pA, f, oA, maxsidebands=16); sB, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, cyclostationary=True)
    sR, _ = pacA.pnoise(pA, f, oA, maxsidebands=16, cyclostationary=True)
    print('WHITE   f/f0=%.2f  identity rel %.1e   reduction rel %.1e' % (f / f0, abs(sB / sA - 1), abs(sR / sA - 1)))
cA = build_f('A'); pA, pacA = solve(cA); oA = [str(n) for n in cA.nodes].index('out')
cB = build_f('B'); pB, pacB = solve(cB); oB = [str(n) for n in cB.nodes].index('out')
for f in (0.13 * f0, 1.37 * f0):
    sA, _ = pacA.pnoise(pA, f, oA, maxsidebands=16); sB, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, cyclostationary=True)
    sBm, _ = pacB.pnoise(pB, f, oB, maxsidebands=16, modulated=True)
    print('FLICKER f/f0=%.2f  A %.8e  B %.8e  identity rel %.1e   B/A %.6f   cycle-avg/A %.3f' % (f / f0, sA, sB, abs(sB / sA - 1), sB / sA, sBm / sA), flush=True)
## diagnostics
for f in (0.13 * f0, 1.37 * f0):
    sA, _ = pacA.pnoise(pA, f, oA, maxsidebands=16); sAc, _ = pacA.pnoise(pA, f, oA, maxsidebands=16, cyclostationary=True)
    print('FLICKER REDUCTION on A  f/f0=%.2f  stationary %.6e  cyclostationary %.6e  rel %.1e' % (f / f0, sA, sAc, abs(sAc / sA - 1)), flush=True)
VA_SMALL = 0.2
src2 = src.replace("va=1.0, vo=0.3", "va=%g, vo=0.3" % VA_SMALL)
def build_f2(kind):
    c = SubCircuit()
    for n in ('lo', 'mid', 'out'):
        c.add_node(n)
    c['vlo'] = VSin('lo', gnd, va=VA_SMALL, vo=0.3, freq=f0)
    if kind == 'A':
        c.add_node('n'); c['xi'] = TAS._Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0); c['Rn'] = R('n', gnd, r=Rn)
        c['M1'] = Mult('mid', gnd, 'n', gnd, 'lo', gnd, k=k1)
    else:
        c['src'] = ModFlicker('mid', gnd, 'lo', gnd, k=k1 * Rn)
    c['Rm'] = R('mid', gnd, r=Rm); c['M2'] = Mult('out', gnd, 'mid', gnd, 'lo', gnd, k=k2)
    c['Ro'] = R('out', gnd, r=Ro); c['Co'] = C('out', gnd, c=Co)
    return c
cA2 = build_f2('A'); pA2, pacA2 = solve(cA2); oA2 = [str(n) for n in cA2.nodes].index('out')
cB2 = build_f2('B'); pB2, pacB2 = solve(cB2); oB2 = [str(n) for n in cB2.nodes].index('out')
for f in (0.13 * f0, 1.37 * f0):
    sA, _ = pacA2.pnoise(pA2, f, oA2, maxsidebands=16); sB, _ = pacB2.pnoise(pB2, f, oB2, maxsidebands=16, cyclostationary=True)
    print('FLICKER va=0.2 (no crossing) f/f0=%.2f  B/A %.6f' % (f / f0, sB / sA), flush=True)
## sign-definite modulation: gain k V_lo^2 on both sides (va = 1, crossing in V_lo but not in the gain)
class Mult2(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='k', desc='gain', unit='A/V^3', default=1.0)]
    @staticmethod
    def analog(p, outp, outn, a, an, b, bn):
        return Contribution(Branch(outp, outn).I, p.k * Branch(a, an).V * Branch(b, bn).V ** 2)
class ModFlicker2(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]
    @staticmethod
    def analog(p, outp, outn, b, bn):
        return Contribution(Branch(outp, outn).I, flicker_noise((p.k * Branch(b, bn).V ** 2) ** 2, 1))
def build_f3(kind):
    c = SubCircuit()
    for n in ('lo', 'mid', 'out'):
        c.add_node(n)
    c['vlo'] = VSin('lo', gnd, va=1.0, vo=0.3, freq=f0)
    if kind == 'A':
        c.add_node('n'); c['xi'] = TAS._Flicker('n', gnd, i=0.0, noisePSD=1.0, fref=1.0); c['Rn'] = R('n', gnd, r=Rn)
        c['M1'] = Mult2('mid', gnd, 'n', gnd, 'lo', gnd, k=k1)
    else:
        c['src'] = ModFlicker2('mid', gnd, 'lo', gnd, k=k1 * Rn)
    c['Rm'] = R('mid', gnd, r=Rm); c['M2'] = Mult('out', gnd, 'mid', gnd, 'lo', gnd, k=k2)
    c['Ro'] = R('out', gnd, r=Ro); c['Co'] = C('out', gnd, c=Co)
    return c
cA3 = build_f3('A'); pA3, pacA3 = solve(cA3); oA3 = [str(n) for n in cA3.nodes].index('out')
cB3 = build_f3('B'); pB3, pacB3 = solve(cB3); oB3 = [str(n) for n in cB3.nodes].index('out')
for f in (0.13 * f0, 1.37 * f0):
    sA, _ = pacA3.pnoise(pA3, f, oA3, maxsidebands=16); sB, _ = pacB3.pnoise(pB3, f, oB3, maxsidebands=16, cyclostationary=True)
    print('FLICKER sign-definite gain (va=1) f/f0=%.2f  B/A %.9f' % (f / f0, sB / sA), flush=True)
