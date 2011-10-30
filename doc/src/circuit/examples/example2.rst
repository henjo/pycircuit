Example 2
----------

Multi feedBack (MFB) filter 
``````````````````````````````````````````````
.. image:: mfb.*

Find symbolic expression of transfer function from Is to V(3,0)

.. sympy::
   :persistent:

   from pycircuit.circuit import *
   from sympy import symbols, simplify, ratsimp, sympify, factor, limit, solve, pprint, fraction, collect, powsimp, powdenest, Add, sqrtdenest, I, Rational
   from sympy.simplify.simplify import fraction_expand


Symbol definition
++++++++++++++++++

.. sympy::
   :persistent:

   R1,R2,R3,C1,C2,i_s,Z_in = symbols('R1 R2 R3 C1 C2 i_s Z_{in}', real=True, positive=True, bounded=True)
   s = Symbol('s', complex = True)   
   w = Symbol('omega', real = True)   

Circuit definition
+++++++++++++++++++

.. sympy::
   :persistent:

   ## Create circuit
   cir = SubCircuit(toolkit=symbolic)
   cir['R1'] = R(1, 3, r = R1)
   cir['R2'] = R(1, 2, r = R2)
   cir['R3'] = R(1, gnd, r = R3)
   cir['C1'] = C(1, gnd, c = C1)
   cir['C2'] = C(2, 3, c = C2)
   cir['Nullor'] = Nullor(2, gnd, 3, gnd)

   # Current source for AC stimuli
   cir['ISource'] = IS(1, gnd, iac=i_s)

AC analysis
++++++++++++

.. sympy::
   :persistent:

   ## Run symbolic AC analysis     
   ac = AC(cir)
   result = ac.solve(freqs=s, complexfreq=True)


Transfer function 
///////////////////

From the voltage source to net 2

.. sympy::
   :persistent:

   simplify(result.v(3, gnd) / i_s)

Denominator of transfer function:

.. sympy::
   :persistent:
   :include-source: False
   
   denominator = fraction(simplify(result.v(3, gnd) / i_s))[1]
   denominator


DC gain
///////////////////

.. sympy::
   :persistent:
   
   simplify(simplify(result.v(3, gnd) / i_s).subs({s:0}))

Poles
///////////////////

.. sympy::
   :persistent:
   :include-source: False

   a = use(powsimp(solve(denominator, s)[0], deep = True), powsimp, level=2) 
   a

.. sympy::
   :persistent:
   :include-source: False

   b = use(powsimp(solve(denominator, s)[1], deep = True), powdenest, level=1)
   b

Poles in omega

.. sympy::
   :persistent:
   :include-source: False

   omega11, omega12, omega13, omega22 = symbols('omega_11 omega_12 omega_13 omega_22', real=True, positive=True, bounded=True)

   aa = expand((fraction_expand(a.args[0])*I)**2)   
   aa = aa.subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   aaa = omega11*omega22+factor(collect(aa-omega11*omega22,[omega11,omega12]))
   aaa

.. sympy::
   :persistent:
   :include-source: False

   aa1 = a.args[1].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   aa2 = a.args[2].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   aa3 = a.args[3].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   aaaa = I*sqrt(aaa)+aa1+aa2+aa3
   aaaa
 
.. sympy::
   :persistent:
   :include-source: False

   ab = expand((fraction_expand(b.args[1])*I)**2)   
   ab = ab.subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   aba = omega11*omega22+factor(collect(ab-omega11*omega22,[omega11,omega12]))
   aba

.. sympy::
   :persistent:
   :include-source: False

   ab1 = b.args[0].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   ab2 = b.args[2].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   ab3 = b.args[3].subs({simplify(1/(C1*R3)):omega13,simplify(1/(C1*R2)):omega12,simplify(1/(C1*R1)):omega11,simplify(1/(C1*R3)):omega13,simplify(1/(C2*R2)):omega22})
   abaa = -I*sqrt(aba)+ab1+ab2+ab3
   abaa

a1

.. sympy::
   :persistent:
   :include-source: False

   a1_s,a0_s,Q_s,omega_0,zeta, G = symbols('a_1,a_0,Q,omega_0 zeta G', real=True, positive=True, bounded=True)

   a1 = -(aaaa + abaa)
   Eq(a1_s,a1)

a0

.. sympy::
   :persistent:
   :include-source: False

   a0 = -((aaaa-abaa)**2-a1**2)/4
   Eq(a0_s,a0)

omega0

.. sympy::
   :persistent:
   :include-source: False

   omega0 = sqrt(a0) 
   Eq(omega_0,omega0)

.. sympy::
   :persistent:
   :include-source: False

   omega0r = omega0.subs({omega13:simplify(1/(C1*R3)),omega12:simplify(1/(C1*R2)),omega11:simplify(1/(C1*R1)),omega13:simplify(1/(C1*R3)),omega22:simplify(1/(C2*R2))}) 
   Eq(omega_0,omega0r)

zeta

.. sympy::
   :persistent:
   :include-source: False

   z = a1/(2*omega0) 
   Eq(zeta,z)

.. sympy::
   :persistent:
   :include-source: False

   zr = z.subs({omega13:simplify(1/(C1*R3)),omega12:simplify(1/(C1*R2)),omega11:simplify(1/(C1*R1)),omega13:simplify(1/(C1*R3)),omega22:simplify(1/(C2*R2))})
   Eq(zeta,zr)

Q

.. sympy::
   :persistent:
   :include-source: False

   Q = omega0/a1
   Eq(Q_s,Q)

.. sympy::
   :persistent:
   :include-source: False

   Qr = Q.subs({omega13:simplify(1/(C1*R3)),omega12:simplify(1/(C1*R2)),omega11:simplify(1/(C1*R1)),omega13:simplify(1/(C1*R3)),omega22:simplify(1/(C2*R2))})
   Eq(Q_s,Qr)

.. sympy::
   :persistent:
   :include-source: False

   apa = solve([omega_0-omega0r,zeta-zr,G-R1],[C1,C2,R1]) 

.. sympy::
   :persistent:
   :include-source: False

   Eq(C1,apa[0][0])

.. sympy::
   :persistent:
   :include-source: False

   Eq(C2,apa[0][1])

.. sympy::
   :persistent:
   :include-source: False

   Eq(C2/C1,ratsimp(apa[0][1]/apa[0][0]))

.. sympy::
   :persistent:
   :include-source: False

   Eq(R1,apa[0][2])

.. sympy::
   :persistent:
   :include-source: False

   tf = ratsimp(simplify(result.v(3, gnd) / i_s).subs({R1:G,C1:apa[0][0],C2:apa[0][1]}))
   tf 


Input impedance
///////////////////

.. sympy::
   :persistent:
   :include-source: False

   zin = simplify(result.v(1, gnd) / i_s)
   Eq(Z_in,simplify(-zin))

Noise analysis
++++++++++++++

Input current noise

.. sympy::
   :persistent:
   :include-source: False

   in_s = symbols('i_{n}', real=True, positive=True, bounded=True)

   ## Run symbolic Noise analysis     
   noise = Noise(cir, inputsrc='ISource', outputnodes=('3', gnd), toolkit=symbolic)
   resultNoise = noise.solve(freqs=0, complexfreq=False)

   a = ratsimp(collect(resultNoise['Sininp'],[noise.toolkit.kboltzmann*noise.par.epar.T]).subs({R1:G}))
   expr = use(a,factor,level=3)
   expr 

.. sympy::
   :persistent:
   :include-source: False

   in_r = collect(expand(expr),[noise.toolkit.kboltzmann*noise.par.epar.T])
   Eq(in_s,in_r) 

.. sympy::
   :persistent:
   :include-source: False

   iin_s,iout_s = symbols('i_{in}, i_{out}', real=True, positive=True, bounded=True)
   inr = collect(expand(ratsimp(solve(iin_s-in_r,R2)[0]*(G+R3)**2/(G*R3))),[R3]) 
   inr = G*R3/((G+R3)**2)*collect(use(inr,factor,level=2),[R3,R3*R3])
   Eq(R2,inr)

.. sympy::
   :persistent:
   :include-source: False

   Eq(R2/Rational(10e3,1),ratsimp(inr.subs({R3:Rational(10e3,1),G:Rational(10e3,1)})/Rational(10e3,1)))
   
Output voltage noise

.. sympy::
   :persistent:
   :include-source: False

   a     = ratsimp(resultNoise['Svnout']).subs({R1:G})
   expr  = use(a,fraction_expand,level=2)
   in_rr = collect(expand(expr),[noise.toolkit.kboltzmann*noise.par.epar.T])
   Eq(iout_s,in_rr) 

.. sympy::
   :persistent:
   :include-source: False

   in_rr = use(ratsimp(in_rr.subs({R2:inr})),factor,level=1)
   Eq(iout_s,in_rr) 
