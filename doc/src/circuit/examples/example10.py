"""Example 10: MFB filter -- symbolic post-processing with GiNaC.

Solves the same Multiple-FeedBack (MFB) active filter as Example 2, but performs
the symbolic post-processing -- transfer function ``N(s)/D(s)``, DC gain, natural
frequency, quality factor and poles -- through the native GiNaC backend
(``GinacToolkit`` + :mod:`pycircuit.circuit._ginac`), replicating the classic
hand ("Symby") analysis but with the collection/cancellation done in GiNaC.

Requires the optional compiled ``_ginac_ext`` extension
(``pycircuit/circuit/build_ginac_ext.sh``).
"""
import sympy
from sympy import Symbol, symbols

from pycircuit.circuit import SubCircuit, R, C, IS, Nullor, gnd
from pycircuit.circuit.analysis_ss import AC
from pycircuit.circuit.toolkit import symbolic, ginac_toolkit
from pycircuit.circuit import _ginac

if ginac_toolkit is None:
    raise SystemExit("GiNaC extension not built; see build_ginac_ext.sh")

R1, R2, R3, C1, C2, i_s = symbols('R1 R2 R3 C1 C2 i_s', real=True, positive=True)
s = Symbol('s', complex=True)


def build_mfb():
    cir = SubCircuit(toolkit=symbolic)
    cir['R1'] = R(1, 3, r=R1)
    cir['R2'] = R(1, 2, r=R2)
    cir['R3'] = R(1, gnd, r=R3)
    cir['C1'] = C(1, gnd, c=C1)
    cir['C2'] = C(2, 3, c=C2)
    cir['Nullor'] = Nullor(2, gnd, 3, gnd)
    cir['ISource'] = IS(1, gnd, iac=i_s)
    return cir


cir = build_mfb()

## Solve the AC transfer function v(3)/i_s through the GiNaC backend.
res = AC(cir, toolkit=ginac_toolkit).solve(s, complexfreq=True)
H = res.v(3, gnd) / i_s

## Collect into N(s)/D(s) coefficient vectors (low order first) -- natively in
## GiNaC (normal() + collect), only the compact coefficients cross back to sympy.
num, den = _ginac.poly_coeffs(H, s)
print("N(s) coeffs (low->high):", num)
print("D(s) coeffs (low->high):", den)

## DC gain = N(0)/D(0)
dc_gain = sympy.simplify(num[0] / den[0])
print("DC gain:", dc_gain)

## Second-order lowpass: D(s) = d0 + d1 s + d2 s^2 -> omega0, Q
d0, d1, d2 = den
omega0 = sympy.simplify(sympy.sqrt(d0 / d2))
Q = sympy.simplify(sympy.sqrt(d0 * d2) / d1)
print("omega0:", omega0)
print("Q:", Q)

## Poles = roots of the (quadratic) denominator
D = sum(c * s**k for k, c in enumerate(den))
poles = [sympy.simplify(p) for p in sympy.solve(D, s)]
print("poles:", poles)

## Cross-check the whole thing against the pure-sympy ("Symby") reference.
H_ref = sympy.simplify(
    AC(cir, toolkit=symbolic).solve(s, complexfreq=True).v(3, gnd) / i_s)
N = sum(c * s**k for k, c in enumerate(num))
assert sympy.simplify(N / D - H_ref) == 0
print("GiNaC result matches the sympy reference.")
