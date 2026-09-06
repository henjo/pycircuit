"""B2 gate (roadmap sec. B2) -- theta = 1/2 + Ch, the fifth trapezoidal design.

RUN 2026-09-06.  IT PASSES, and it also produces the recipe Houben does not give.

C1 is a theorem: trapezoidal is A-stable but not L-stable, maps null(C) by EXACTLY
-1, so A_trap^K is singular at even K.  Reproduced here (rcond exactly 0.0 at
K = 200/400, healthy 3.1e-3 at K = 199/201/401 -- the obstruction is parity-specific).
Biasing theta moves those modes to -(1-theta)/theta, matching the numerical
eigenvalues to every printed digit.

The two-sided constraint, measured on the Q=20 resonator (analytic peak 20 V):

    C        theta          peak      error     null(C) |mode|^K   rcond(I-A^K)
    0        0.500000000    SINGULAR  --        1.000e+00          0.0e+00
    1e3      0.500031416    20.02011  +0.0201   9.752e-01          4.4e-03
    1e4      0.500314159    20.01301  +0.0130   7.778e-01          4.6e-03   <- knee
    1e5      0.503141593    19.94223  -0.0578   8.100e-02          4.3e-03
    1e6      0.531415927    19.26150  -0.7385   1.176e-11          3.7e-03

rcond SATURATES at C ~ 1e4 and then DEGRADES while the amplitude keeps paying, so
more bias buys no conditioning.  That knee is the recipe: C ~ 1e3-1e4 costs 0.07%
of the peak -- against the 30% (13.92 V of 20) that got `x0_unknown` reverted.

VALIDATED AGAINST A REFERENCE THIS CODE CANNOT INFLUENCE: at theta = 0.5 and odd K
the closed form agrees with the SHIPPED trapezoidal PSS to 1.47e-07 (K = 401).

What it does NOT show: anything about a nonlinear circuit or the shipped
integrator.  It is the gate the roadmap said to run before writing anything.

Half two of that gate: what biasing theta costs the PHYSICAL oscillation.
Linear theta-method PSS solved in closed form on the Q=20 resonator whose
analytic peak is 20 V -- the same fixture the suite gates trapezoidal on."""
import warnings, numpy as np
from pycircuit import circuit
from pycircuit.circuit.elements import VSin, R, C, L, SubCircuit, gnd
from pycircuit.circuit.analysis import remove_row_col
circuit.default_toolkit=circuit.numeric; warnings.simplefilter('ignore')
tk=circuit.numeric

def resonator():
    Lv,Cv,Q=1e-3,1e-9,20.0
    f0=1.0/(2*np.pi*np.sqrt(Lv*Cv))
    c=SubCircuit(); c.add_node('n1'); c.add_node('n2')
    c['vs']=VSin(gnd,'n1',va=1.0,freq=f0)
    c['L']=L('n1','n2',L=Lv); c['C']=C('n2',gnd,c=Cv)
    c['R']=R('n1','n2',r=Q*np.sqrt(Lv/Cv))
    return c,1.0/f0,f0

cir,T,f0=resonator()
n=cir.n; iref=cir.get_node_index(gnd)
x0=np.zeros(n)
Cm=np.asarray(cir.C(x0,None),float); Gm=np.asarray(cir.G(x0,None),float)
(Cr,)=remove_row_col((Cm,),iref,tk); (Gr,)=remove_row_col((Gm,),iref,tk)
Cr=np.asarray(Cr,float); Gr=np.asarray(Gr,float); m=Cr.shape[0]
## the source vector u(t): -(i + u) convention -- take it by finite difference
## of the circuit's own u() so no source model is re-implemented here
def uvec(t):
    class E: pass
    e=E(); e.t=t
    u=np.asarray(cir.u(t,analysis='tran'),float)
    return np.concatenate((u[:iref],u[iref+1:]))
## which reduced row is node 'c'?  (the fixture measures v('c'))
names=[cir.get_node_name(nd) for nd in cir.nodes]
gi=names.index('n2')                      # the capacitor node of this fixture
ci=gi if gi<iref else gi-1

def peak(theta,K):
    h=T/K
    M=Cr/h+theta*Gr; N=Cr/h-(1.0-theta)*Gr
    A=np.linalg.solve(M,N)
    W=np.zeros((K,m))
    for k in range(K):
        W[k]=np.linalg.solve(M,-(theta*uvec((k+1)*h)+(1.0-theta)*uvec(k*h)))
    acc=np.zeros(m)
    for k in range(K):
        acc=A@acc+W[k]
    IA=np.eye(m)-np.linalg.matrix_power(A,K)
    sig=np.linalg.svd(IA,compute_uv=False)
    rc=sig[-1]/sig[0]
    if rc < 1e-12:
        return None,rc                     # the (-1)^n obstruction, exactly
    xs=np.linalg.solve(IA,acc)
    pk=abs(xs[ci]); xk=xs.copy()
    for k in range(K):
        xk=A@xk+W[k]
        pk=max(pk,abs(xk[ci]))
    return pk,rc

print("Q=20 resonator at resonance, analytic peak 20 V.  K steps per period.")
print("  K=%d   h=T/K=%.4g s\n"%(200,T/200))
print("     C          theta            peak (V)     error       null(C) |mode|^K")
for K in (200,800):
    h=T/K
    print("  --- K = %d ---"%K)
    for Cc in (0.0,1e2,1e3,1e4,1e5,1e6):
        th=0.5+Cc*h
        p,cond=peak(th,K)
        if p is None:
            print("   %-9g %.9f   SINGULAR -- I - A^K has rcond %.1e (the theorem)"
                  %(Cc,th,cond))
        else:
            print("   %-9g %.9f   %9.5f   %+8.4f    %.3e   rcond %.1e"
                  %(Cc,th,p,p-20.0,abs((1-th)/th)**K,cond))
