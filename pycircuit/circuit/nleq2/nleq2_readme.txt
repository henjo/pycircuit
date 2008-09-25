NLEQ2 package - release 2.3 at January 3, 1992

Extract from nleq2.f
====================================================================

C*  Title
C
C     Numerical solution of nonlinear (NL) equations (EQ)
C     especially designed for numerically sensitive problems.
C
C*  Written by        U. Nowak, L. Weimann 
C*  Purpose           Solution of systems of highly nonlinear equations
C*  Method            Damped affine invariant Newton method with rank-
C                     strategy (see references below)
C*  Category          F2a. - Systems of nonlinear equations
C*  Keywords          Nonlinear equations, Newton methods
C*  Version           2.3
C*  Revision          September 1991
C*  Latest Change     July 2000
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad-Zuse-Zentrum fuer
C                     Informationstechnik Berlin (ZIB)
C                     Takustrasse 7, D-14195 Berlin-Dahlem
C                     phone : + 49/30/84185-0
C                     fax   : + 49/30/84185-125
C*  Contact           Lutz Weimann
C                     ZIB, Division Scientific Computing, 
C                          Department Scientific Software
C                     phone : + 49/30/84185-185
C                     fax   : + 49/30/84185-107
C                     e-mail: weimann@zib.de
C
C*    References:
C
C     /1/ P. Deuflhard:
C         Newton Techniques for Highly Nonlinear Problems -
C         Theory and Algorithms.
C         Academic press Inc. (To be published)
C
C     /2/ U. Nowak, L. Weimann:
C         A Family of Newton Codes for Systems of Highly Nonlinear
C         Equations - Algorithm, Implementation, Application.
C         ZIB, Technical Report TR 90-10 (December 1990)
C
C  ---------------------------------------------------------------
C
C* Licence
C    You may use or modify this code for your own non commercial
C    purposes for an unlimited time. 
C    In any case you should not deliver this code without a special 
C    permission of ZIB.
C    In case you intend to use the code commercially, we oblige you
C    to sign an according licence agreement with ZIB.
C
C* Warranty 
C    This code has been tested up to a certain level. Defects and
C    weaknesses, which may be included in the code, do not establish
C    any warranties by ZIB. ZIB does not take over any liabilities
C    which may follow from acquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C     ------------------------------------------------------------
C
C*    Summary:
C     ========
C     Damped Newton-algorithm with rank strategy for systems of 
C     highly nonlinear equations - damping strategy due to Ref.(1).
C
C     (The iteration is done by subroutine N2INT currently. NLEQ2
C      itself does some house keeping and builds up workspace.)
C
C     Jacobian approximation by numerical differences or user
C     supplied subroutine JAC.
C
C     The numerical solution of the arising linear equations is
C     done by means of the subroutines DECCON and SOLCON (QR de-
C     composition with subcondition estimation, rank decision and
C     computation of the rank-deficient pseudoinverse) .
C     For special purposes these routines may be substituted.
C
C     This is a driver routine for the core solver N2INT.
C
C
