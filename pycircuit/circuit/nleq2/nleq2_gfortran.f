      SUBROUTINE NLEQ2(N,FCN,JAC,X,XSCAL,RTOL,IOPT,IERR,
     $LIWK,IWK,LRWK,RWK)
C*    Begin Prologue NLEQ2
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
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
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    External subroutines (to be supplied by the user)
C     =================================================
C 
C     (Caution: Arguments declared as (input) must not
C               be altered by the user subroutines ! )
C
C     FCN(N,X,F,IFAIL) Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       X(N)           Dble   Vector of unknowns (input)
C       F(N)           Dble   Vector of function values (output)
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a value <= 2.
C                             If <0: NLEQ2 will be terminated with 
C                                    error code = 82, and IFAIL stored
C                                    to IWK(23).
C                             If =1: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduct. factor, which 
C                                    must be output through F(1) by FCN,
C                                    and it's value must be >0 and < 1.
C                             Note, that if IFAIL = 1 or 2, additional
C                             conditions concerning the damping factor,
C                             e.g. the minimum damping factor or the
C                             bounded damping strategy may also influ-
C                             ence the value of the reduced damping 
C                             factor.
C
C
C     JAC(N,LDJAC,X,DFDX,IFAIL) 
C                       Ext    Jacobian matrix subroutine
C       N                 Int    Number of vector components (input)
C       LDJAC             Int    Leading dimension of Jacobian array
C                                (input)
C       X(N)              Dble   Vector of unknowns (input)
C       DFDX(LDJAC,N)     Dble   DFDX(i,k): partial derivative of
C                                I-th component of FCN with respect
C                                to X(k) (output)
C       IFAIL             Int    JAC evaluation-failure indicator. 
C                                (output)
C                                Has always value 0 (zero) on input.
C                                Indicates failure of JAC evaluation
C                                and causes termination of NLEQ2,
C                                if set to a negative value on output
C
C
C*    Input parameters of NLEQ2
C     =========================
C
C     N              Int    Number of unknowns
C   * X(N)           Dble   Initial estimate of the solution
C   * XSCAL(N)       Dble   User scaling (lower threshold) of the 
C                           iteration vector X(N)
C   * RTOL           Dble   Required relative precision of
C                           solution components -
C                           RTOL.GE.EPMACH*TEN*N
C   * IOPT(50)       Int    Array of run-time options. Set to zero
C                           to get default values (details see below)
C
C*    Output parameters of NLEQ2
C     ==========================
C
C   * X(N)           Dble   Solution values ( or final values,
C                           respectively )
C   * XSCAL(N)       Dble   After return with IERR.GE.0, it contains
C                           the latest internal scaling vector used
C                           After return with IERR.EQ.-1 in onestep-
C                           mode it contains a possibly adapted 
C                           (as described below) user scaling vector:
C                           If (XSCAL(I).LT. SMALL) XSCAL(I) = SMALL ,
C                           If (XSCAL(I).GT. GREAT) XSCAL(I) = GREAT .
C                           For SMALL and GREAT, see section machine
C                           constants below  and regard note 1.
C   * RTOL           Dble   Finally achieved (relative) accuracy
C                           The estimated absolute error of component i
C                           of x_out is approximately given by
C                             abs_err(i) = RTOL * XSCAL_out(i) ,
C                           where (approximately)
C                             XSCAL_out(i) = 
C                                max(abs(X_out(i)),XSCAL_in(i)).
C     IERR           Int    Return value parameter
C                           =-1 sucessfull completion of one iteration
C                               step, subsequent iterations are needed 
C                               to get a solution. (stepwise mode only)
C                           = 0 successfull completion of iteration
C                           > 0 see list of error messages below
C
C     Note 1.
C        The machine dependent values SMALL, GREAT and EPMACH are
C        gained from calls of the machine constants function D1MACH.
C        As delivered, this function is adapted to use constants 
C        suitable for all machines with IEEE arithmetic. If you use
C        another type of machine, you have to change the DATA state-
C        ments for IEEE arithmetic in D1MACH into comments and to 
C        uncomment the set of DATA statements suitable for your machine.
C
C*    Workspace parameters of NLEQ2
C     =============================
C
C     LIWK           Int    Declared dimension of integer
C                           workspace.
C                           Required minimum (for standard linear system
C                           solver) : N+52
C   * IWK(LIWK)      Int    Integer Workspace
C     LRWK           Int    Declared dimension of real workspace.
C                           Required minimum (for standard linear system
C                           solver and Jacobian computed by numerical
C                           approximation - if the Jacobian is computed
C                           by a user subroutine JAC, decrease the 
C                           expression noted below by N):
C                           (N+NBROY+15)*N+61
C                           NBROY = Maximum number of Broyden steps
C                           (Default: if Broyden steps are enabled, e.g.
C                                                IOPT(32)=1            -
C                                       NBROY=MAX(N,10)
C                                     else (if IOPT(32)=0) - 
C                                       NBROY=0 ;
C                            see equally named IWK-field below)
C   * RWK(LRWK)      Dble   Real Workspace
C
C     Note 2a.  A test on sufficient workspace is made. If this
C               test fails, IERR is set to 10 and an error-message
C               is issued from which the minimum of required
C               workspace size can be obtained.
C
C     Note 2b.  The first 50 elements of IWK and RWK are partially 
C               used as input for internal algorithm parameters (for
C               details, see below). In order to set the default values
C               of these parameters, the fields must be set to zero.
C               Therefore, it's recommanded always to initialize the
C               first 50 elements of both workspaces to zero.
C
C*   Options IOPT:
C    =============
C
C     Pos. Name   Default  Meaning
C
C       1  QSUCC  0        =0 (.FALSE.) initial call:
C                             NLEQ2 is not yet initialized, i.e. this is
C                             the first call for this nonlinear system.
C                             At successfull return with MODE=1,
C                             QSUCC is set to 1.
C                          =1 (.TRUE.) successive call:
C                             NLEQ2 is initialized already and is now
C                             called to perform one or more following
C                             Newton-iteration steps.
C                             ATTENTION:
C                                Don't destroy the contents of
C                                IOPT(i) for 1 <= i <= 50 ,
C                                IWK(j)  for 1 <= j < NIWKFR and
C                                RWK(k)  for 1 <= k < NRWKFR.
C                                (Nevertheless, some of the options, e.g.
C                                 FCMIN, SIGMA, MPR..., can be modified
C                                 before successive calls.)
C       2  MODE   0        =0 Standard mode initial call:
C                             Return when the required accuracy for the
C                             iteration vector is reached. User defined
C                             parameters are evaluated and checked.
C                             Standard mode successive call:
C                             If NLEQ2 was called previously with MODE=1,
C                             it performs all remaining iteration steps.
C                          =1 Stepwise mode:
C                             Return after one Newton iteration step.
C       3  JACGEN 0        Method of Jacobian generation
C                          =0 Standard method is JACGEN=2
C                          =1 User supplied subroutine JAC will be 
C                             called to generate Jacobian matrix
C                          =2 Jacobian approximation by numerical
C                             differentation (see subroutine N2JAC)
C                          =3 Jacobian approximation by numerical
C                             differentation with feedback control
C                             (see subroutine N2JCF)
C       4..8               Reserved
C       9  ISCAL  0        Determines how to scale the iterate-vector:
C                          =0 The user supplied scaling vector XSCAL is
C                             used as a (componentwise) lower threshold
C                             of the current scaling vector
C                          =1 The vector XSCAL is always used as the
C                             current scaling vector
C      10                  Reserved
C      11  MPRERR 0        Print error messages
C                          =0 No output
C                          =1 Error messages
C                          =2 Warnings additionally
C                          =3 Informal messages additionally
C      12  LUERR  6        Logical unit number for error messages
C      13  MPRMON 0        Print iteration Monitor
C                          =0 No output
C                          =1 Standard output
C                          =2 Summary iteration monitor additionally
C                          =3 Detailed iteration monitor additionally
C                          =4,5,6 Outputs with increasing level addi-
C                             tional increasing information for code
C                             testing purposes. Level 6 produces
C                             in general extremely large output!
C      14  LUMON  6        Logical unit number for iteration monitor
C      15  MPRSOL 0        Print solutions
C                          =0 No output
C                          =1 Initial values and solution values
C                          =2 Intermediate iterates additionally
C      16  LUSOL  6        Logical unit number for solutions
C      17..18              Reserved
C      19  MPRTIM 0        Output level for the time monitor
C                          = 0 : no time measurement and no output
C                          = 1 : time measurement will be done and
C                                summary output will be written -
C                                regard note 4a.
C      20  LUTIM  6        Logical output unit for time monitor
C      21..30              Reserved
C      31  NONLIN 3        Problem type specification
C                          =1 Linear problem
C                             Warning: If specified, no check will be
C                             done, if the problem is really linear, and
C                             NLEQ2 terminates unconditionally after one
C                             Newton-iteration step.
C                          =2 Mildly nonlinear problem
C                          =3 Highly nonlinear problem
C                          =4 Extremely nonlinear problem
C      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
C                             approximation are inhibited.
C                          =1 (.TRUE.) Rank-1 updates by Broyden-
C                             approximation are allowed.
C      33..34              Reserved
C      35  QNSCAL 0        Inhibit automatic row scaling: 
C                          =0 (.FALSE.) Automatic row scaling of
C                             the linear system is activ: 
C                             Rows i=1,...,N will be divided by
C                             max j=1,...,N (abs(a(i,j))) 
C                          =1 (.TRUE.) No row scaling of the linear
C                             system. Recommended only for well row-
C                             scaled nonlinear systems.
C      36..37              Reserved
C      38  IBDAMP          Bounded damping strategy switch:
C                          =0 The default switch takes place, dependent
C                             on the setting of NONLIN (=IOPT(31)):
C                             NONLIN = 0,1,2,3 -> IBDAMP = off ,
C                             NONLIN = 4 -> IBDAMP = on
C                          =1 means always IBDAMP = on 
C                          =2 means always IBDAMP = off 
C      39  IORMON          Convergence order monitor 
C                          =0 Standard option is IORMON=2 
C                          =1 Convergence order is not checked,
C                             the iteration will be always proceeded
C                             until the solution has the required 
C                             precision RTOL (or some error condition
C                             occured)
C                          =2 Use additional 'weak stop' criterion:
C                             Convergence order is monitored
C                             and termination due to slowdown of the
C                             convergence may occur.
C                          =3 Use additional 'hard stop' criterion:
C                             Convergence order is monitored
C                             and termination due to superlinear 
C                             convergence slowdown may occur. 
C                          In case of termination due to convergence
C                          slowdown, the warning code IERR=4 will be
C                          set.
C                          In cases, where the Newton iteration con-
C                          verges but superlinear convergence order has
C                          never been detected, the warning code IERR=5 
C                          is returned.
C      40..45              Reserved
C      46..50              User options (see note 4b)
C
C     Note 3:
C         If NLEQ2 terminates with IERR=2 (maximum iterations)
C         or  IERR=3 (small damping factor), you may try to continue
C         the iteration by increasing NITMAX or decreasing FCMIN
C         (see RWK) and setting QSUCC to 1.
C
C     Note 4a:
C        The integrated time monitor calls the machine dependent
C        subroutine SECOND to get the current time stamp in form
C        of a real number (Single precision). As delivered, this
C        subroutine always return 0.0 as time stamp value. Refer
C        to the compiler- or library manual of the FORTRAN compiler
C        which you currently use to find out how to get the current
C        time stamp on your machine.
C
C     Note 4b:
C         The user options may be interpreted by the user replacable
C         routines N2SOUT, N2FACT, N2SOLV - the distributed version
C         of N2SOUT currently uses IOPT(46) as follows:
C         0 = standard plotdata output (may be postprocessed by a user-
C             written graphical program)
C         1 = plotdata output is suitable as input to the graphical
C             package GRAZIL (based on GKS), which has been developed
C             at ZIB. 
C
C
C*   Optional INTEGER input/output in IWK:
C    =======================================
C
C     Pos. Name          Meaning
C
C      1   NITER  IN/OUT Number of Newton-iterations
C      2                 reserved
C      3   NCORR  IN/OUT Number of corrector steps
C      4   NFCN   IN/OUT Number of FCN-evaluations
C      5   NJAC   IN/OUT Number of Jacobian generations or
C                        JAC-calls
C      6                 reserved
C      7                 reserved
C      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
C                        approximation
C      9   NREJR1 IN/OUT Number of rejected Newton iteration steps
C                        done with a rank-1 approximated Jacobian
C     10..11             Reserved
C     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
C                        fication number ppppvvvv, consisting of the
C                        program code pppp and the version code vvvv.
C                        Input: If containing a negative number,
C                        it will only be overwritten by the identi-
C                        fication number, immediately followed by
C                        a return to the calling program.      
C     13..15             Reserved
C     16   NIWKFR OUT    First element of IWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solver: N+53
C     17   NRWKFR OUT    First element of RWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solver and numerically 
C                        approximated Jacobian computed by the 
C                        expression: (N+9+NBROY)*N+62 
C                        If the Jacobian is computed by a user routine
C                        JAC, subtract N in this expression.
C     18   LIWKA  OUT    Length of IWK currently required
C     19   LRWKA  OUT    Length of RWK currently required
C     20..22             Reserved
C     23   IFAIL  OUT    Set in case of failure of N2FACT (IERR=80),
C                        N2SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
C                        to the nonzero IFAIL value returned by the 
C                        routine indicating the failure .
C     24..30             Reserved
C     31   NITMAX IN     Maximum number of permitted iteration
C                        steps (default: 50)
C     32   IRANK  IN     Initial rank of the Jacobian 
C                        (at the iteration starting point)
C                        =0 : full rank N
C                        =k with min_rank <= k < N : 
C                           deficient rank assumed,
C                        where min_rank = max (1,N-max(N/10,10))
C     33   NEW    IN/OUT Count of consecutive rank-1 updates
C     34..35             Reserved
C     36   NBROY  IN     Maximum number of possible consecutive 
C                        iterative Broyden steps. The total real 
C                        workspace needed (RWK) depends on this value
C                        (see LRWK above).
C                        Default is N (see parameter N),
C                        if MSTOR=0 (=IOPT(4)), 
C                        and ML+MU+1 (=IOPT(6)+IOPT(7)+1), if MSTOR=1
C                        (but minimum is always 10) - 
C                        provided that Broyden is allowed. 
C                        If Broyden is inhibited, NBROY is always set to
C                        zero.
C     37..50             Reserved
C
C*   Optional REAL input/output in RWK:
C    ====================================
C
C     Pos. Name          Meaning
C
C      1..16             Reserved
C     17   CONV   OUT    The achieved relative accuracy after the  
C                        current step
C     18   SUMX   OUT    Natural level (not Normx of printouts)
C                        of the current iterate, i.e. Sum(DX(i)**2),
C                        where DX = scaled Newton correction.
C     19   DLEVF  OUT    Standard level (not Normf of printouts)
C                        of the current iterate, i.e. Norm2(F(X)),
C                        where F =  nonlinear problem function.
C     20   FCBND  IN     Bounded damping strategy restriction factor
C                        (Default is 10)
C     21   FCSTRT IN     Damping factor for first Newton iteration -
C                        overrides option NONLIN, if set (see note 5)
C     22   FCMIN  IN     Minimal allowed damping factor (see note 5)
C     23   SIGMA  IN     Broyden-approximation decision parameter
C                        Required choice: SIGMA.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C                        Rank-1 updates are inhibited, if 
C                        SIGMA.GT.1/FCMIN is set. (see note 5)
C     24   SIGMA2 IN     Decision parameter about increasing damping
C                        factor to corrector if predictor is small.
C                        Required choice: SIGMA2.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs rank-1 updates.
C     25   COND   IN     Maximum permitted subcondition for rank-
C                        decision by linear solver.
C                        (Default: 1/epmach, epmach: relative
C                         machine precision) 
C     26   AJDEL  IN     Jacobian approximation without feedback:
C                        Relative pertubation for components
C                        (Default: sqrt(epmach*10), epmach: relative
C                         machine precision) 
C     27   AJMIN  IN     Jacobian approximation without feedback:
C                        Threshold value (Default: 0.0d0)
C                          The absolute pertubation for component k is
C                          computed by 
C                          DELX := AJDEL*max(abs(Xk),AJMIN)
C     28  ETADIF  IN     Jacobian approximation with feedback:
C                        Target value for relative pertubation ETA of X
C                        (Default: 1.0d-6)
C     29  ETAINI  IN     Jacobian approximation with feedback:
C                        Initial value for denominator differences
C                        (Default: 1.0d-6)
C     30..50             Reserved
C
C     Note 5:
C       The default values of the internal parameters may be obtained
C       from the monitor output with at least IOPT field MPRMON set to 2
C       and by initializing the corresponding RWK-fields to zero. 
C
C*   Error and warning messages:
C    ===========================
C
C      1    Termination at stationary point (rank deficient Jacobian
C           and termination criterion fulfilled)
C      2    Termination after NITMAX iterations ( as indicated by
C           input parameter NITMAX=IWK(31) )
C      3    Termination, since damping factor became to small and
C           rank of Jacobian matrix is already zero
C      4    Warning: Superlinear or quadratic convergence slowed down
C           near the solution.
C           Iteration has been stopped therefore with an approximation
C           of the solution not such accurate as requested by RTOL,
C           because possibly the RTOL requirement may be too stringent
C           (i.e. the nonlinear problem is ill-conditioned)
C      5    Warning: Iteration stopped with termination criterion 
C           (using RTOL as requested precision) satisfied, but no 
C           superlinear or quadratic convergence has been indicated yet.
C           Therefore, possibly the error estimate for the solution may
C           not match good enough the really achieved accuracy.
C     10    Integer or real workspace too small
C     20    Bad input to dimensional parameter N
C     21    Nonpositive value for RTOL supplied
C     22    Negative scaling value via vector XSCAL supplied
C     30    One or more fields specified in IOPT are invalid
C           (for more information, see error-printout)
C     80    Error signalled by linear solver routine N2FACT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used by standard routine N2FACT)
C     81    Error signalled by linear solver routine N2SOLV,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used by standard routine N2SOLV)
C     82    Error signalled by user routine FCN (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     83    Error signalled by user routine JAC (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C
C     Note 6 : in case of failure:
C        -    use non-standard options
C        -    use another initial guess
C        -    or reformulate model
C        -    or apply continuation techniques
C
C*    Machine dependent constants used:
C     =================================
C
C     DOUBLE PRECISION EPMACH  in  NLEQ2, N2PCHK, N2INT
C     DOUBLE PRECISION GREAT   in  N2PCHK
C     DOUBLE PRECISION SMALL   in  N2PCHK, N2INT, N2SCAL
C
C*    Subroutines called: N2PCHK, N2INT, D1MACH
C
C     ------------------------------------------------------------
C*    End Prologue
C
C*    Summary of changes:
C     ===================
C      
C     2.2.1  91, June  3    Time monitor included
C     2.2.2  91, June  3    Bounded damping strategy implemented
C     2.2.3  91, July 26    AJDEL, AJMIN as RWK-options for JACGEN.EQ.2,
C                           ETADIF, ETAINI as RWK-opts. for JACGEN.EQ.3
C                           FCN-count changed for anal. Jacobian
C     2.2.4  91, August 16  Convergence order monitor included
C     2.2.5  91, August 19  Standard Broyden updates replaced by
C                           iterative Broyden
C            91, Sept.      Rank strategy modified
C                           DECCON with new fail exit, for the case that
C                           the square root of a negative number will
C                           appear during pseudo inverse computation.
C                           (Occured, although theoretical impossible!)
C     2.2.6  91, Sept.  17  Damping factor reduction by FCN-fail imple-
C                           mented
C     2.3    91, Dec.   20  New Release for CodeLib
C            00, July   12  RTOL output-value bug fixed
C   
C     ------------------------------------------------------------
C
C     PARAMETER (IRWKI=xx, LRWKI=yy)  
C     IRWKI: Start position of internally used RWK part
C     LRWKI: Length of internally used RWK part
C     (current values see parameter statement below)
C
C     INTEGER L4,L41,L5,L51,L6,L61,L62,L63,L7,L71,L9,L11,L12,L121,
C             L13,L14,L20
C     Starting positions in RWK of formal array parameters of internal
C     routine N1INT (dynamically determined in driver routine NLEQ1,
C     dependent on N and options setting)
C
C     Further RWK positions (only internally used)
C
C     Position  Name     Meaning
C
C     IRWKI     FCKEEP   Damping factor of previous successfull iter.
C     IRWKI+1   FCA      Previous damping factor
C     IRWKI+2   FCPRI    A priori estimate of damping factor
C     IRWKI+3   DMYCOR   Number My of latest corrector damping factor
C                        (kept for use in rank-1 decision criterium)
C     IRWKI+(4..LRWKI-1) Free
C
C     Internal arrays stored in RWK (see routine N2INT for descriptions)
C
C     Position  Array         Type   Remarks
C
C     L4        QA(N,N)       Perm   Arrays QA and DXSAVE need never to
C     L4        DXSAVE(N,NBROY)      be kept the same time and therefore
C                             Perm   may be stored to the same workspace
C                                    part
C     L41       A(N,N)        Perm
C     L5        DX(N)         Perm  
C     L51       DXQ(N)        Perm 
C     L6        XA(N)         Perm
C     L61       F(N)          Perm
C     L62       FW(N)         Perm
C     L63       XWA(N)        Perm
C     L7        FA(N)         Perm
C     L71       ETA(N)        Perm   Only used for JACGEN=IOPT(3)=3
C     L8                      Perm   Start position of array workspace 
C                                    needed for linear solver  
C     L9        XW(N)         Temp
C     L11       DXQA(N)       Temp
C     L111      QU(N)         Temp
C     L12       T1(N)         Temp
C     L121      T2(N)         Temp
C     L13       T3(N)         Temp   Not yet used or even reserved
C                                    (for future band mode implementat.)
C     
C
      EXTERNAL D1MACH,N2INT
      INTRINSIC DBLE
      INTEGER IRWKI, LRWKI
      PARAMETER (IRWKI=51, LRWKI=10)  
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER IRANK,NITMAX,LUERR,LUMON,LUSOL,MPRERR,MPRMON,
     $MPRSOL,M1,M2,NRWKFR,NRFRIN,NRW,NIWKFR,NIFRIN,NIW,NONLIN,JACGEN
      INTEGER L4,L41,L5,L51,L6,L61,L62,L63,L7,L71,L8,L9,L11,L12,L121,
     $        L13,L20
      DOUBLE PRECISION COND,D1MACH,FC,FCMIN,PERCI,PERCR
      DOUBLE PRECISION EPMACH
      LOGICAL QINIMO,QRANK1,QFCSTR,QSUCC,QBDAMP
      CHARACTER CHGDAT*20, PRODCT*8
C     Which version ?
      LOGICAL QVCHK
      INTEGER IVER
      PARAMETER( IVER=21122301 )
C
C     Version: 2.3               Latest change:
C     -----------------------------------------
C
      DATA      CHGDAT      /'July 12, 2000       '/
      DATA      PRODCT      /'NLEQ2   '/
C*    Begin
      EPMACH = D1MACH(3)
      IERR = 0
      QVCHK = IWK(12).LT.0
      IWK(12) = IVER
      IF (QVCHK) RETURN
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .EQ. 0) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C        Print iteration monitor?
      MPRMON = IOPT(13)
      LUMON = IOPT(14)
      IF (LUMON .LE. 0 .OR. LUMON .GT. 99) THEN
        LUMON = 6
        IOPT(14)=LUMON
      ENDIF
C        Print intermediate solutions?
      MPRSOL = IOPT(15)
      LUSOL = IOPT(16)
      IF (LUSOL .EQ. 0) THEN
        LUSOL = 6
        IOPT(16)=LUSOL
      ENDIF
C        Print time summary statistics?
      MPRTIM = IOPT(19)
      LUTIM = IOPT(20)
      IF (LUTIM .EQ. 0) THEN
        LUTIM = 6
        IOPT(20)=LUTIM
      ENDIF
      QSUCC = IOPT(1).EQ.1
      QINIMO = MPRMON.GE.1.AND..NOT.QSUCC
C     Print NLEQ2 heading lines
      IF(QINIMO)THEN
10000   FORMAT('    N L E Q 2  *****  V e r s i o n  ',
     $         '2 . 3 ***',//,1X,'Newton-Method ',
     $         'for the solution of nonlinear systems',//)
        WRITE(LUMON,10000)
      ENDIF
C     Check input parameters and options
      CALL N2PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C     Exit, if any parameter error was detected till here
      IF (IERR.NE.0) RETURN 
      M1=N
      M2=N
      JACGEN=IOPT(3)
      IF (JACGEN.EQ.0) JACGEN=2
      IOPT(3)=JACGEN
      QRANK1=IOPT(32).EQ.1
      IF (QRANK1) THEN
        NBROY=IWK(36)
        IF (NBROY.EQ.0) NBROY=MAX(M2,10)
        IWK(36)=NBROY
      ELSE
        NBROY=0
      ENDIF
C     WorkSpace: RWK
      L4=IRWKI+LRWKI
      L41=L4+NBROY*N
      L5=L41+M1*N
      L51=L5+N
      L6=L51+N
      L61=L6+N
      L62=L61+N
      L63=L62+N
      L7=L63+N
      L71=L7+N
      IF (JACGEN.NE.3) THEN
        L8=L71
      ELSE
        L8=L71+N
      ENDIF
      NRWKFR = L8
      L9=LRWK+1-N
      L11=L9-N
      L111=L11-N
      L12=L111-N
      L121=L12-N
C     L13 : Work array T3, for future band mode implementation
      L13=L121
      NRW=NRWKFR+LRWK-L13+1
C     End WorkSpace at NRW
C     WorkSpace: IWK
      L20=51
      NIWKFR = L20
      NIW = NIWKFR-1
      NIFRIN = NIWKFR
      NRFRIN = NRWKFR
      LIWL=N+2
      LRWL=2*N+1
      IF (QRANK1) THEN
        NRWKFR=NRWKFR+LRWL
        NIWKFR=NIWKFR+LIWL
      ENDIF
      NRW = NRW+LRWL
      NIW = NIW+LIWL
C     End WorkSpace at NIW
      IWK(16) = NIWKFR
      IWK(17) = NRWKFR
C
      IF(NRW.GT.LRWK.OR.NIW.GT.LIWK)THEN
        IERR=10
      ELSE
        IF(QINIMO)THEN
          PERCR = DBLE(NRW)/DBLE(LRWK)*100.0D0
          PERCI = DBLE(NIW)/DBLE(LIWK)*100.0D0
C         Print statistics concerning workspace usage
10050     FORMAT(' Real    Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//,
     $    ' Integer Workspace declared as ',I9,
     $    ' is used up to ',I9,' (',F5.1,' percent)',//)
          WRITE(LUMON,10050)LRWK,NRW,PERCR,LIWK,NIW,PERCI
        ENDIF
        IF(QINIMO)THEN
10051     FORMAT(/,' N =',I4,//,' Prescribed relative ',
     $    'precision',D10.2,/)
          WRITE(LUMON,10051)N,RTOL
10052     FORMAT(' The Jacobian is supplied by ',A)
          IF (JACGEN.EQ.1) THEN
            WRITE(LUMON,10052) 'a user subroutine'
          ELSE IF (JACGEN.EQ.2) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (without feedback strategy)'
          ELSE IF (JACGEN.EQ.3) THEN
             WRITE(LUMON,10052) 
     $        'numerical differentiation (feedback strategy included)'
          ENDIF
10057     FORMAT(' Automatic row scaling of the Jacobian is ',A,/)
          IF (IOPT(35).EQ.1) THEN
            WRITE(LUMON,10057) 'inhibited'
          ELSE
            WRITE(LUMON,10057) 'allowed'
          ENDIF
        ENDIF
        NONLIN=IOPT(31)
        IF (IOPT(38).EQ.0) QBDAMP = NONLIN.EQ.4
        IF (IOPT(38).EQ.1) QBDAMP = .TRUE.
        IF (IOPT(38).EQ.2) QBDAMP = .FALSE.
        IF (QBDAMP) THEN
          IF (RWK(20).LT.ONE) RWK(20) = TEN
        ENDIF
10064   FORMAT(' Rank-1 updates are ',A)
        IF (QINIMO) THEN
          IF (QRANK1) THEN
            WRITE(LUMON,10064) 'allowed'
          ELSE
            WRITE(LUMON,10064) 'inhibited'
          ENDIF
10065     FORMAT(' Problem is specified as being ',A)
          IF (NONLIN.EQ.1) THEN
            WRITE(LUMON,10065) 'linear'
          ELSE IF (NONLIN.EQ.2) THEN
            WRITE(LUMON,10065) 'mildly nonlinear'
          ELSE IF (NONLIN.EQ.3) THEN
            WRITE(LUMON,10065) 'highly nonlinear'
          ELSE IF (NONLIN.EQ.4) THEN
            WRITE(LUMON,10065) 'extremely nonlinear'
          ENDIF
10066     FORMAT(' Bounded damping strategy is ',A,:,/, 
     $           ' Bounding factor is ',D10.3)
          IF (QBDAMP) THEN
            WRITE(LUMON,10066) 'active', RWK(20)
          ELSE
            WRITE(LUMON,10066) 'off'
          ENDIF
        ENDIF
C       Maximum permitted number of iteration steps
        NITMAX=IWK(31)
        IF (NITMAX.LE.0) NITMAX=50
        IWK(31)=NITMAX
10068   FORMAT(' Maximum permitted number of iteration steps : ',
     $         I6)
        IF (QINIMO) WRITE(LUMON,10068) NITMAX
C       Initial damping factor for highly nonlinear problems
        QFCSTR=RWK(21).GT.ZERO
        IF (.NOT.QFCSTR) THEN
          RWK(21)=1.0D-2
          IF (NONLIN.EQ.4) RWK(21)=1.0D-4
        ENDIF
C       Minimal permitted damping factor
        IF (RWK(22).LE.ZERO) THEN
          RWK(22)=1.0D-4
          IF (NONLIN.EQ.4) RWK(22)=1.0D-8
        ENDIF
        FCMIN=RWK(22)
C       Rank1 decision parameter SIGMA
        IF (RWK(23).LT.ONE) RWK(23)=3.0D0
        IF (.NOT.QRANK1) RWK(23)=10.0D0/FCMIN
C       Decision parameter about increasing too small predictor
C       to greater corrector value
        IF (RWK(24).LT.ONE) RWK(24)=10.0D0/FCMIN       
C       Starting value of damping factor (FCMIN.LE.FC.LE.1.0)
        IF(NONLIN.LE.2.AND..NOT.QFCSTR)THEN
C         for linear or mildly nonlinear problems
          FC = ONE
        ELSE
C         for highly or extremely nonlinear problems
          FC = RWK(21)
        ENDIF
        RWK(21)=FC
C       Initial rank
        IRANK = IWK(32)
        IF (IRANK.LE.0.OR.IRANK.GT.N) IWK(32) = N
C       Maximum permitted sub condition number of matrix A
        COND = RWK(25)
        IF (COND.LT.ONE) COND = ONE/EPMACH
        RWK(25) = COND
        IF (MPRMON.GE.2.AND..NOT.QSUCC) THEN
10069     FORMAT(//,' Internal parameters:',//,
     $      ' Starting value for damping factor FCSTART = ',D9.2,/,
     $      ' Minimum allowed damping factor FCMIN = ',D9.2,/,
     $      ' Rank-1 updates decision parameter SIGMA = ',D9.2,/,
     $      ' Initial Jacobian pseudo-rank IRANK =',I6,/,
     $      ' Maximum permitted subcondition COND = ',D9.2)
          WRITE(LUMON,10069) RWK(21),FCMIN,RWK(23),IWK(32),COND
        ENDIF
C       Store lengths of currently required workspaces
        IWK(18) = NIFRIN-1
        IWK(19) = NRFRIN-1
C
C       Initialize and start time measurements monitor
C
        IF ( IOPT(1).EQ.0 .AND. MPRTIM.NE.0 ) THEN
          CALL MONINI (' NLEQ2',LUTIM)
          CALL MONDEF (0,'NLEQ2')
          CALL MONDEF (1,'FCN')
          CALL MONDEF (2,'Jacobi')
          CALL MONDEF (3,'Lin-Fact')
          CALL MONDEF (4,'Lin-Sol')
          CALL MONDEF (5,'Output')
          CALL MONSRT ()
        ENDIF
C
C
        IERR=-1
C       If IERR is unmodified on exit, successive steps are required
C       to complete the Newton iteration
        IF (NBROY.EQ.0) NBROY=1
        CALL N2INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IWK(32),IOPT,
     $  IERR,LRWK,RWK,NRFRIN,LRWL,LIWK,IWK,NIFRIN,LIWL,M1,M2,NBROY,
     $  RWK(L4),RWK(L41),RWK(L4),RWK(L5),RWK(L51),RWK(L6),RWK(L63),
     $  RWK(L61),
     $  RWK(L7),RWK(L71),RWK(L9),RWK(L62),RWK(L11),RWK(L111),RWK(L12),
     $  RWK(L121),RWK(L13),RWK(21),RWK(22),RWK(23),RWK(24),RWK(IRWKI+1),
     $  RWK(IRWKI),RWK(IRWKI+2),COND,RWK(IRWKI+3),RWK(17),RWK(18),
     $  RWK(19),MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,IWK(1),IWK(3),
     $  IWK(4),IWK(5),IWK(8),IWK(9),IWK(33),QBDAMP)
C
      IF (MPRTIM.NE.0.AND.IERR.NE.-1.AND.IERR.NE.10) CALL MONEND
C
C       Free workspaces, so far not used between steps
        IWK(16) = NIWKFR
        IWK(17) = NRWKFR
      ENDIF
C     Print statistics
      IF (MPRMON.GE.1.AND.IERR.NE.-1.AND.IERR.NE.10) THEN
10080   FORMAT(/, '   ******  Statistics * ', A8, ' *******', /,
     $            '   ***  Newton iterations : ', I7,'  ***', /,
     $            '   ***  Corrector steps   : ', I7,'  ***', /,
     $            '   ***  Rejected rk-1 st. : ', I7,'  ***', /,
     $            '   ***  Jacobian eval.    : ', I7,'  ***', /,
     $            '   ***  Function eval.    : ', I7,'  ***', /,
     $            '   ***  ...  for Jacobian : ', I7,'  ***', /,
     $            '   *************************************', /)
        WRITE (LUMON,10080) PRODCT,IWK(1),IWK(3),IWK(9),IWK(5),
     $  IWK(4),IWK(8)
      ENDIF
C     Print workspace requirements, if insufficient
      IF (IERR.EQ.10) THEN
10090   FORMAT(///,20('*'),'Workspace Error',20('*'))
        IF (MPRERR.GE.1) WRITE(LUERR,10090)
        IF(NRW.GT.LRWK)THEN
10091     FORMAT(/,' Real Workspace dimensioned as',1X,I9,
     $    1X,'must be enlarged at least up to ',
     $    I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10091)LRWK,NRFRIN-1
        ENDIF
        IF(NIW.GT.LIWK)THEN
10092     FORMAT(/,' Integer Workspace dimensioned as ',
     $    I9,' must be enlarged at least up ',
     $    'to ',I9,//)
          IF (MPRERR.GE.1) WRITE(LUERR,10092)LIWK,NIFRIN-1
        ENDIF
      ENDIF
C     End of subroutine NLEQ2
      RETURN
      END
C
      SUBROUTINE N2PCHK(N,X,XSCAL,RTOL,IOPT,IERR,LIWK,IWK,LRWK,RWK)
C*    Begin Prologue N2PCHK
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 2 P C H K : Checking of input parameters and options
C                   for NLEQ2.
C
C*    Parameters:
C     ===========
C
C     See parameter description in driver routine.
C
C*    Subroutines called: D1MACH
C
C*    Machine dependent constants used:
C     =================================
C
C     EPMACH = relative machine precision
C     GREAT = squareroot of maxreal divided by 10
C     SMALL = squareroot of "smallest positive machine number
C             divided by relative machine precision"
      DOUBLE PRECISION EPMACH,GREAT,SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
C
      EXTERNAL D1MACH
      INTRINSIC DBLE
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=1.0D1)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C
      PARAMETER (NUMOPT=50)
      INTEGER IOPTL(NUMOPT),IOPTU(NUMOPT)
      DOUBLE PRECISION D1MACH,TOLMIN,TOLMAX,DEFSCL
C
      DATA IOPTL /0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,
     $            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,0,0,0,
     $            -9999999,-9999999,-9999999,-9999999,-9999999/
      DATA IOPTU /1,1,3,0,0,0,0,0,1,0,3,99,6,99,3,99,0,0,1,99,
     $            0,0,0,0,0,0,0,0,0,0,4,1,0,0,1,
     $            0,0,2,3,0,0,0,0,0,0,
     $            9999999,9999999,9999999,9999999,9999999/
C
      EPMACH = D1MACH(3)
      GREAT  = DSQRT(D1MACH(2)/TEN)
      SMALL  = D1MACH(6)
      IERR = 0
C        Print error messages?
      MPRERR = IOPT(11)
      LUERR = IOPT(12)
      IF (LUERR .LE. 0 .OR. LUERR .GT. 99) THEN
        LUERR = 6
        IOPT(12)=LUERR
      ENDIF
C
C     Checking dimensional parameter N
      IF ( N.LE.0 ) THEN
        IF (MPRERR.GE.1)  WRITE(LUERR,10011) N
10011   FORMAT(/,' Error: Bad input to dimensional parameter N supplied'
     $         ,/,8X,'choose N positive, your input is: N = ',I5)
        IERR = 20
      ENDIF
C
C     Problem type specification by user
      NONLIN=IOPT(31)
      IF (NONLIN.EQ.0) NONLIN=3
      IOPT(31)=NONLIN
C
C     Checking and conditional adaption of the user-prescribed RTOL
      IF (RTOL.LE.ZERO) THEN
        IF (MPRERR.GE.1) 
     $      WRITE(LUERR,'(/,A)') ' Error: Nonpositive RTOL supplied'
        IERR = 21
      ELSE
        TOLMIN = EPMACH*TEN*DBLE(N)
        IF(RTOL.LT.TOLMIN) THEN
          RTOL = TOLMIN
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'increased ','smallest',RTOL
        ENDIF
        TOLMAX = 1.0D-1
        IF(RTOL.GT.TOLMAX) THEN
          RTOL = TOLMAX
          IF (MPRERR.GE.2) 
     $      WRITE(LUERR,10012) 'decreased ','largest',RTOL
        ENDIF
10012   FORMAT(/,' Warning: User prescribed RTOL ',A,'to ',
     $         'reasonable ',A,' value RTOL = ',D11.2)
      ENDIF
C     
C     Test user prescribed accuracy and scaling on proper values
      IF (N.LE.0) RETURN 
      IF (NONLIN.GE.3) THEN
        DEFSCL = RTOL
      ELSE
        DEFSCL = ONE
      ENDIF
      DO 10 I=1,N
        IF (XSCAL(I).LT.ZERO) THEN
          IF (MPRERR.GE.1) THEN 
            WRITE(LUERR,10013) I
10013       FORMAT(/,' Error: Negative value in XSCAL(',I5,') supplied')
          ENDIF
          IERR = 22
        ENDIF
        IF (XSCAL(I).EQ.ZERO) XSCAL(I) = DEFSCL
        IF ( XSCAL(I).GT.ZERO .AND. XSCAL(I).LT.SMALL ) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10014) I,XSCAL(I),SMALL
10014       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too small, ',
     $             'increased to',D9.2)
          ENDIF
          XSCAL(I) = SMALL
        ENDIF
        IF (XSCAL(I).GT.GREAT) THEN
          IF (MPRERR.GE.2) THEN
            WRITE(LUERR,10015) I,XSCAL(I),GREAT
10015       FORMAT(/,' Warning: XSCAL(',I5,') = ',D9.2,' too big, ',
     $             'decreased to',D9.2)
          ENDIF
          XSCAL(I) = GREAT
        ENDIF
10    CONTINUE
C     Checks options
      DO 20 I=1,30
        IF (IOPT(I).LT.IOPTL(I) .OR. IOPT(I).GT.IOPTU(I)) THEN
          IERR=30
          IF (MPRERR.GE.1) THEN
            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';',
     $             /,3X,'range of permitted values is ',I8,' to ',I8)
          ENDIF
        ENDIF
20    CONTINUE
C     End of subroutine N2PCHK
      RETURN
      END
C
      SUBROUTINE N2INT(N,FCN,JAC,X,XSCAL,RTOL,NITMAX,NONLIN,IRANK,IOPT,
     $IERR,LRWK,RWK,NRWKFR,LRWL,LIWK,IWK,NIWKFR,LIWL,M1,M2,NBROY,
     $QA,A,DXSAVE,DX,DXQ,XA,XWA,F,FA,ETA,XW,FW,DXQA,QU,T1,T2,T3,FC,
     $FCMIN,SIGMA,SIGMA2,FCA,FCKEEP,FCPRI,COND,DMYCOR,CONV,SUMX,DLEVF,
     $MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,NCORR,NFCN,NJAC,
     $NFCNJ,NREJR1,NEW,QBDAMP)
C*    Begin Prologue N2INT
      INTEGER N
      EXTERNAL FCN,JAC
      DOUBLE PRECISION X(N),XSCAL(N)
      DOUBLE PRECISION RTOL
      INTEGER NITMAX,NONLIN,IRANK
      INTEGER IOPT(50)
      INTEGER IERR
      INTEGER LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER NRWKFR,LRWL,LIWK
      INTEGER IWK(LIWK)
      INTEGER NIWKFR,LIWL,M1,M2,NBROY
      DOUBLE PRECISION QA(M2,N),A(M1,N),DXSAVE(N,NBROY)
      DOUBLE PRECISION DX(N),DXQ(N),XA(N),XWA(N),F(N),FA(N),ETA(N)
      DOUBLE PRECISION XW(N),FW(N),DXQA(N),QU(N),T1(N),T2(N),T3(N)
      DOUBLE PRECISION FC,FCMIN,SIGMA,SIGMA2,FCA,FCKEEP,COND,CONV,SUMX,
     $                 DLEVF,FCPRI,DMYCOR
      INTEGER MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL,NITER,
     $NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW
      LOGICAL QBDAMP
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 2 I N T : Core routine for NLEQ2 .
C     Damped Newton-algorithm with rank-strategy for systems of 
C     highly nonlinear equations especially designed for
C     numerically sensitive problems.
C
C*    Parameters:
C     ===========
C
C       N,FCN,JAC,X,XSCAL,RTOL   
C                         See parameter description in driver routine
C
C       NITMAX     Int    Maximum number of allowed iterations
C       NONLIN     Int    Problem type specification
C                         (see IOPT-field NONLIN)
C       IRANK      Int    Initially proposed (in) and final (out) rank
C                         of Jacobian
C       IOPT       Int    See parameter description in driver routine
C       IERR       Int    See parameter description in driver routine
C       LRWK       Int    Length of real workspace
C       RWK(LRWK)  Dble   Real workspace array
C       NRWKFR     Int    First free position of RWK on exit 
C       LIWK       Int    Length of integer workspace
C       IWK(LIWK)  Int    Integer workspace array
C       NIWKFR     Int    First free position of IWK on exit 
C       M1         Int    Leading dimension of Jacobian array A
C                         for full case Jacobian: N
C                         (other matrix types are not yet implemented)
C       M2         Int    Leading dimension of Jacobian array QA
C                         for full case Jacobian: N
C       NBROY      Int    Maximum number of possible consecutive
C                         iterative Broyden steps. (See IWK(36))
C       QA(M2,N)   Dble   Holds the originally computed Jacobian
C                         or the pseudo inverse in case of rank-
C                         deficiency
C       A(M1,N)    Dble   Holds the Jacobian matrix (decomposed form
C                         after call of linear decomposition
C                         routine)
C       DXSAVE(X,NBROY)
C                  Dble   Used to save the quasi Newton corrections of
C                         all previously done consecutive Broyden
C                         steps.
C       DX(N)      Dble   Current Newton correction
C       DXQ(N)     Dble   Simplified Newton correction J(k-1)*X(k)
C       XA(N)      Dble   Previous Newton iterate
C       XWA(N)     Dble   Scaling factors used for latest decomposed
C                         Jacobian for column scaling - may differ
C                         from XW, if Broyden updates are performed
C       F(N)       Dble   Function (FCN) value of current iterate
C       FA(N)      Dble   Function (FCN) value of previous iterate
C       ETA(N)     Dble   Jacobian approximation: updated scaled
C                         denominators
C       XW(N)      Dble   Scaling factors for iteration vector
C       FW(N)      Dble   Scaling factors for rows of the system
C       DXQA(N)    Dble   Previous Newton correction
C       QU(N)      Dble   Savespace for right hand side belonging
C                         to upper triangular linear system
C       T1(N)      Dble   Workspace for linear solvers and internal
C                         subroutines
C       T2(N)      Dble   Workspace array for internal subroutines
C       T3(N)      Dble   Workspace array for internal subroutines
C       FC         Dble   Current Newton iteration damping factor.
C       FCMIN      Dble   Minimum permitted damping factor. If
C                         FC becomes smaller than this value, one
C                         of the following may occur:
C                         a.    Recomputation of the Jacobian
C                               matrix by means of difference
C                               approximation (instead of Rank1
C                               update), if Rank1 - update
C                               previously was used
C                         b.    Rank reduction of Jacobi
C                               matrix ,  if difference
C                               approximation was used previously
C                               and Rank(A).NE.0
C                         c.    Fail exit otherwise
C       SIGMA      Dble   Decision parameter for rank1-updates.
C       SIGMA2     Dble   Decision parameter for damping factor
C                         increasing to corrector value
C       FCA        Dble   Previous Newton iteration damping factor.
C       FCKEEP     Dble   Keeps the damping factor as it is at start
C                          of iteration step.
C       COND       Dble   Maximum permitted subcondition for rank-
C                         decision by linear solver.
C       CONV       Dble   Scaled maximum norm of the Newton-
C                         correction. Passed to RWK-field on output.
C       SUMX       Dble   Square of the natural level (see equal-
C                         named IOPT-output field)
C       DLEVF      Dble   Square of the standard level (see equal-
C                         named IOPT-output field)
C       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL :
C                         See description of equal named IOPT-fields
C                         in the driver subroutine
C       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
C                         See description of equal named IWK-fields
C                         in the driver subroutine
C       QBDAMP     Logic  Flag, that indicates, whether bounded damping
C                         strategy is active:
C                         .true.  = bounded damping strategy is active
C                         .false. = normal damping strategy is active
C
C
C*    Internal double variables
C     =========================
C
C       AJDEL    See RWK(26) (num. diff. without feedback)
C       AJMIN    See RWK(27) (num. diff. without feedback)
C       COND1    Gets the subcondition of the linear system
C                as estimated by the linear solver (N2FACT)
C       CONVA    Holds the previous value of CONV .
C       DEL      Gets the projection defect in case of rank-
C                deficiency.
C       DMUE     Temporary value used during computation of damping 
C                factors predictor.
C       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
C       ETADIF   See description of RWK(28) (num. diff. with feedback)
C       ETAINI   Initial value for all ETA-components (num. diff. fb.)
C       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
C       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
C       FCDNM    Used to compute the denominator of the damping 
C                factor FC during computation of it's predictor,
C                corrector and aposteriori estimate (in the case of
C                performing a Rank1 update) .
C       FCK2     Aposteriori estimate of FC.
C       FCMIN2   FCMIN**2 . Used for FC-predictor computation.
C       FCMINH   DSQRT(FCMIN).
C                Used in rank decision logical expression.
C       FCNUMP   Gets the numerator of the predictor formula for FC.
C       FCNMP2   Temporary used for predictor numerator computation.
C       FCNUMK   Gets the numerator of the corrector computation 
C                of FC .
C       SENS1    Gets the sensitivity of the Jacobian as
C                estimated by the linear solver N2FACT.
C       SUMXA    Natural level of the previous iterate.
C       TH       Temporary variable used during corrector- and 
C                aposteriori computations of FC.
C
C*    Internal integer variables
C     ==========================
C
C     IFAIL      Gets the return value from subroutines called from
C                N2INT (N2FACT, N2SOLV, FCN, JAC)
C     ISCAL      Holds the scaling option from the IOPT-field ISCAL      
C     MODE       Matrix storage mode (see IOPT-field MODE) 
C     NRED       Count of successive corrector steps
C     NILUSE     Gets the amount of IWK used by the linear solver
C     NRLUSE     Gets the amount of RWK used by the linear solver
C     NIWLA      Index of first element of IWK provided to the
C                linear solver
C     NRWLA      Index of first element of RWK provided to the
C                linear solver
C     LIWL       Holds the maximum amount of integer workspace
C                available to the linear solver
C     LRWL       Holds the maximum amount of real workspace
C                available to the linear solver
C
C*    Internal logical variables
C     ==========================
C
C     QGENJ      Jacobian updating technique flag:
C                =.TRUE.  : Call of analytical subroutine JAC or
C                           numerical differentiation
C                =.FALSE. : rank1- (Broyden-) update
C     QINISC     Iterate initial-scaling flag:
C                =.TRUE.  : at first call of N2SCAL
C                =.FALSE. : at successive calls of N2SCAL
C     QSUCC      See description of IOPT-field QSUCC.
C     QJCRFR     Jacobian refresh flag:
C                set to .TRUE. if damping factor gets too small
C                and Jacobian was computed by rank1-update. 
C                Indicates, that the Jacobian needs to be recomputed
C                by subroutine JAC or numerical differentation.
C     QLINIT     Initialization state of linear solver workspace:
C                =.FALSE. : Not yet initialized
C                =.TRUE.  : Initialized - N2FACT has been called at
C                           least one time.
C     QREPET     Operation mode flag for linear solver:
C                =.FALSE. : Normal operation (full rank matrix)
C                =.TRUE.  : Special operation in rank deficient case:
C                           Compute rank-deficient pseudo-inverse,
C                           partial recomputation when solving the
C                           linear system again prescribing a lower
C                           rank as before.
C     QNEXT      Set to .TRUE. to indicate success of the current
C                Newton-step, i.e. : sucessfull monotonicity-test.
C     
C     QREDU      Set to .TRUE. to indicate that rank-reduction (or
C                refreshment of the Jacobian) is needed - if the
C                computed damping factor gets too small.
C     QSCALE     Holds the value of .NOT.QNSCAL. See description
C                of IOPT-field QNSCAL.
C
C*    Subroutines called:
C     ===================
C
C       N2FACT, N2SOLV, N2JAC,  N2JCF,  N2LVLS, N2PRJN,
C       N2SCRF, N2SOUT, N2PRV1, N2PRV2, N2SCAL,
C       MONON,  MONOFF
C
C*    Functions called:
C     =================
C
C       D1MACH, WNORM
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION EPMACH,SMALL
C 
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL N2FACT, N2SOLV, N2JAC,  N2JCF, N2LVLS, N2PRJN,
     $         N2SCRF, N2SOUT, N2PRV1, N2PRV2, N2SCAL,
     $         MONON,  MONOFF, D1MACH, WNORM
      INTRINSIC DSQRT,DMIN1
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TEN
      PARAMETER (TEN=10.0D0)
      INTEGER IFAIL,ILOOP,ISCAL,K,MODE,NRED,NILUSE,NRLUSE,NIWLA,
     $NRWLA,L1,JACGEN,MINRNK
      DOUBLE PRECISION AJDEL,AJMIN,ALFA1,ALFA2,ALFA,BETA,
     $COND1,CONVA,DLEVXA,DEL,DMYPRI,D1MACH,DXANRM,DXNRM,WNORM,EPDIFF,
     $ETAMIN,ETAMAX,ETAINI,ETADIF,FCDNM,FCBND,FCBH,FCK2,FCMIN2,FCMINH,
     $FCNUMP,FCCOR,FCNMP2,FCNUMK,FCREDU,DLEVFN,SENS1,SUMXA,SUM1,SUM2,
     $S1,TH,RSMALL,APREC
      LOGICAL QGENJ,QINISC,QSUCC,QJCRFR,QLINIT,QNEXT,QRANK1,QREPET,
     $        QREDU,QREP,QSCALE,QMIXIO
CWEI
      INTRINSIC DLOG
      DOUBLE PRECISION CLIN0,CLIN1,CALPHA,CALPHK,ALPHAE,ALPHAK,ALPHAA,
     $                 SUMXA0,SUMXA1,SUMXA2,SUMXTE,FCMON,DLOG
      INTEGER ICONV, IORMON
      LOGICAL QMSTOP
      SAVE CLIN0,CLIN1,CALPHA,ALPHAE,ALPHAK,ALPHAA,SUMXA0,SUMXA1,SUMXA2,
     $     ICONV,QMSTOP
C
      EPMACH = D1MACH(3)
      SMALL  = D1MACH(6)
C*    Begin
C       ----------------------------------------------------------
C       1 Initialization
C       ----------------------------------------------------------
C       1.1 Control-flags and -integers
        QSUCC = IOPT(1).EQ.1
        QSCALE = .NOT. IOPT(35).EQ.1
        QRANK1 = IOPT(32).EQ.1
        IORMON = IOPT(39)
        IF (IORMON.EQ.0) IORMON=2
        ISCAL = IOPT(9)
        MODE = IOPT(2)
        JACGEN = IOPT(3)
        QMIXIO = LUMON.EQ.LUSOL .AND. MPRMON.NE.0 .AND. MPRSOL.NE.0
        MPRTIM = IOPT(19)
C       ----------------------------------------------------------
C       1.2 Derivated dimensional parameters
        MINRNK = MAX0(1,N-MAX0(INT(FLOAT(N)/10.0),10))
C       ----------------------------------------------------------
C       1.3 Derivated internal parameters
        FCMIN2 = FCMIN*FCMIN
        FCMINH = DSQRT(FCMIN)
        TOLMIN = DSQRT(TEN*EPMACH)
        RSMALL = DSQRT(TEN*RTOL)
C       ----------------------------------------------------------
C       1.4 Adaption of input parameters, if necessary
        IF(FC.LT.FCMIN) FC = FCMIN
        IF(FC.GT.ONE) FC = ONE
C       ----------------------------------------------------------
C       1.5 Initial preparations
        QJCRFR = .FALSE.
        QLINIT = .FALSE.
        QREPET = .FALSE.
        IFAIL = 0
        FCBND = ZERO
        IF (QBDAMP) FCBND = RWK(20)
C       ----------------------------------------------------------
C       1.5.1 Numerical differentation related initializations
        IF (JACGEN.EQ.2) THEN
          AJDEL = RWK(26)
          IF (AJDEL.LE.SMALL) AJDEL = DSQRT(EPMACH*TEN)
          AJMIN = RWK(27)
        ELSE IF (JACGEN.EQ.3) THEN
          ETADIF = RWK(28)
          IF (ETADIF .LE. SMALL) ETADIF = 1.0D-6
          ETAINI = RWK(29)
          IF (ETAINI .LE. SMALL) ETAINI = 1.0D-6
          EPDIFF = DSQRT(EPMACH*TEN)
          ETAMAX = DSQRT(EPDIFF)
          ETAMIN = EPDIFF*ETAMAX
        ENDIF
C       ----------------------------------------------------------
C       1.5.2 Miscellaneous preparations of first iteration step
        IF (.NOT.QSUCC) THEN
          NITER = 0
          NCORR = 0
          NREJR1 = 0
          NFCN = 0
          NJAC = 0
          NFCNJ = 0
          QGENJ = .TRUE.
          QINISC = .TRUE.
          FCKEEP = FC
          FCA = FC
          FCPRI = FC
          FCK2 = FC
          CONV = ZERO
          IF (JACGEN.EQ.3) THEN
            DO 1520 L1=1,N
              ETA(L1)=ETAINI
1520        CONTINUE
          ENDIF
          DO 1521 L1=1,N
            XA(L1)=X(L1)
1521      CONTINUE
CWEI      
          ICONV = 0
          ALPHAE = ZERO
          SUMXA1 = ZERO
          SUMXA0 = ZERO
          CLIN0  = ZERO
          QMSTOP = .FALSE.
C         ------------------------------------------------------
C         1.6 Print monitor header
          IF(MPRMON.GE.2 .AND. .NOT.QMIXIO)THEN
16003       FORMAT(///,2X,66('*'))
            WRITE(LUMON,16003)
16004       FORMAT(/,8X,'It',7X,'Normf ',10X,'Normx ',8X,
     $             'Damp.Fct.',3X,'New',6X,'Rank',8X,'Cond')
            WRITE(LUMON,16004)
          ENDIF
C         --------------------------------------------------------
C         1.7 Startup step
C         --------------------------------------------------------
C         1.7.1 Computation of the residual vector
          IF (MPRTIM.NE.0) CALL MONON(1)
          CALL FCN(N,X,F,IFAIL)
          IF (MPRTIM.NE.0) CALL MONOFF(1)
          NFCN = NFCN+1
C     Exit, if ...
          IF (IFAIL.NE.0) THEN
            IERR = 82
            GOTO 4299
          ENDIF
        ELSE
          QINISC = .FALSE.
        ENDIF
C
C       Main iteration loop
C       ===================
C
C       Repeat
2       CONTINUE
C         --------------------------------------------------------
C         2 Startup of iteration step
          IF (.NOT.QJCRFR) THEN
C           ------------------------------------------------------
C           2.1 Scaling of variables X(N)
            CALL N2SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
            QINISC = .FALSE.
            IF(NITER.NE.0)THEN
C             Preliminary pseudo-rank
              IRANKA = IRANK
C             ----------------------------------------------------
C             2.2 Aposteriori estimate of damping factor
              DO 2200 L1=1,N
                DXQA(L1)=DXQ(L1)
2200          CONTINUE
              FCNUMP = ZERO
              DO 2201 L1=1,N
                FCNUMP=FCNUMP+(DX(L1)/XW(L1))**2
2201          CONTINUE
              TH = FC-ONE
              FCDNM = ZERO
              DO 2202 L1=1,N
                FCDNM=FCDNM+((DXQA(L1)+TH*DX(L1))/XW(L1))**2
2202          CONTINUE
C             ----------------------------------------------------
C             2.2.2 Decision criterion for Jacobian updating
C                   technique:
C                   QGENJ.EQ..TRUE. numerical differentation,
C                   QGENJ.EQ..FALSE. rank1 updating
              QGENJ = .TRUE.
              IF (FC.EQ.FCPRI) THEN
                QGENJ = FC.LT.ONE.OR.FCA.LT.ONE.OR.DMYCOR.LE.FC*SIGMA
     $                  .OR. .NOT.QRANK1 .OR. NEW+2.GT.NBROY 
                FCA = FC
              ELSE
                DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMP/FCDNM)
                IF (NONLIN.LE.3) THEN
                  FCCOR = DMIN1(ONE,DMYCOR)
                ELSE
                  FCCOR = DMIN1(ONE,HALF*DMYCOR)
                ENDIF
                FCA = DMAX1(DMIN1(FC,FCCOR),FCMIN)
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,22201) FCCOR, FC, DMYCOR, FCNUMP,
     $                               FCDNM
22201             FORMAT (/, ' +++ aposteriori estimate +++', /,
     $                    ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                    ' DMYCOR = ', D18.10, '  FCNUMP = ', D18.10, /,
     $                    ' FCDNM  = ', D18.10, /,
     $                       ' ++++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
              ENDIF
              FCK2 = FCA
C             ------------------------------------------------------
C             2.2.1 Computation of the numerator of damping
C                   factor predictor
              FCNMP2 = ZERO
              DO 221 L1=1,N
                FCNMP2=FCNMP2+(DXQA(L1)/XW(L1))**2
221           CONTINUE
              FCNUMP = FCNUMP*FCNMP2
            ENDIF
          ENDIF
          QJCRFR =.FALSE.
C         --------------------------------------------------------
C         2.3 Jacobian matrix (stored to array A(M1,N))
C         --------------------------------------------------------
C         2.3.1 Jacobian generation by routine JAC or
C               difference approximation (If QGENJ.EQ..TRUE.)
C               - or -
C               Rank-1 update of Jacobian (If QGENJ.EQ..FALSE.)
          IF(QGENJ)THEN
            NEW = 0
            IF (JACGEN.EQ.1) THEN
               IF (MPRTIM.NE.0) CALL MONON(2)
               CALL JAC(N,M1,X,A,IFAIL)
               IF (MPRTIM.NE.0) CALL MONOFF(2)
            ELSE
              IF (MPRTIM.NE.0) CALL MONON(2)
              IF (JACGEN.EQ.3) 
     $          CALL N2JCF(FCN,N,M1,X,F,A,XW,ETA,ETAMIN,ETAMAX,
     $                     ETADIF,CONV,NFCNJ,T1,IFAIL)
              IF (JACGEN.EQ.2) 
     $          CALL N2JAC(FCN, N, M1, X, F, A, XW,  AJDEL, AJMIN,
     $                     NFCNJ, T1, IFAIL)
              IF (MPRTIM.NE.0) CALL MONOFF(2)
             ENDIF
            NJAC = NJAC + 1
C     Exit, If ...
            IF (JACGEN.EQ.1 .AND. IFAIL.LT.0) THEN
              IERR = 83
              GOTO 4299
            ENDIF
            IF (JACGEN.NE.1 .AND. IFAIL.NE.0) THEN
              IERR = 82
              GOTO 4299
            ENDIF
          ELSE
            NEW = NEW+1
          ENDIF
          IF ( NEW.EQ.0 ) THEN
C           ------------------------------------------------------
C           2.3.2 Save scaling values
            DO 232 L1=1,N
              XWA(L1) = XW(L1)
232         CONTINUE
C           --------------------------------------------------------
C           2.4 Prepare solution of the linear system
C           --------------------------------------------------------
C           2.4.1 internal column scaling of matrix A
            DO 2410 K=1,N
              S1 =-XW(K)
              DO 2412 L1=1,N
                A(L1,K)=A(L1,K)*S1
2412          CONTINUE
2410        CONTINUE
C           ------------------------------------------------------
C           2.4.2 Row scaling of matrix A
            IF (QSCALE) THEN
              CALL N2SCRF(N,N,A,FW)
            ELSE
              DO 242 K=1,N
                FW(K)=ONE
242           CONTINUE
            ENDIF
          ENDIF
C         --------------------------------------------------------
C         2.4.3 Save and scale values of F(N)
          DO 243 L1=1,N
            FA(L1)=F(L1)
            T1(L1)=F(L1)*FW(L1)
243       CONTINUE
          IRANKA = IRANK
          IF (NITER.NE.0) IRANK = N
          QNEXT = .FALSE.
C         --------------------------------------------------------
C         3 Central part of iteration step
C
C         Pseudo-rank reduction loop
C         ==========================
C         DO (Until)
3         CONTINUE
C           --------------------------------------------------------
C           3.1 Solution of the linear system
C           --------------------------------------------------------
C           3.1.1 Decomposition of (N,N)-matrix A
            IF (.NOT.QLINIT) THEN
              NIWLA = NIWKFR
              NRWLA = NRWKFR
            ENDIF
            IF (NEW.EQ.0) THEN
              COND1 = COND
              IF (QREPET) THEN
                IWK(NIWLA) = 1
              ELSE
                IWK(NIWLA) = 0
              ENDIF
              IF (MPRTIM.NE.0) CALL MONON(3)
              CALL N2FACT(N,M1,N,ML,MU,A,QA,COND1,IRANK,IOPT,IFAIL,LIWL,
     $                    IWK(NIWLA),NILUSE,LRWL,RWK(NRWLA),NRLUSE)
              IF (MPRTIM.NE.0) CALL MONOFF(3)
C     Exit Repeat If ...
              IF(IFAIL.NE.0) THEN
                IERR = 80
                GOTO 4299
              ENDIF
              IF (.NOT.QLINIT) THEN
                NIWKFR = NIWKFR+NILUSE
                NRWKFR = NRWKFR+NRLUSE
C               Store lengths of currently required workspaces
                IWK(18) = NIWKFR-1
                IWK(19) = NRWKFR-1
              ENDIF
              SENS1 = RWK(NRWLA)
            ENDIF
            QLINIT = .TRUE.
C           --------------------------------------------------------
C           3.1.2 Solution of linear (N,N)-system
            IF(NEW.EQ.0) THEN 
              IF (MPRTIM.NE.0) CALL MONON(4)
              CALL N2SOLV(N,M1,N,ML,MU,A,QA,T1,T2,IRANK,IOPT,IFAIL,LIWL,
     $                   IWK(NIWLA),IDUMMY,LRWL,RWK(NRWLA),IDUMMY)
              IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
              IF(IFAIL.NE.0)  THEN
                IERR = 81
                GOTO 4299
              ENDIF
              IF(.NOT.QREPET.AND.IRANK.NE.0)THEN
                DO 312 L1=1,N
                  QU(L1)=T1(L1)
312             CONTINUE
              ENDIF
            ELSE  
              ALFA1=ZERO
              ALFA2=ZERO
              DO 3121 I=1,N
                ALFA1=ALFA1+DX(I)*DXQ(I)/XW(I)**2
                ALFA2=ALFA2+DX(I)**2/XW(I)**2
3121          CONTINUE
              ALFA=ALFA1/ALFA2
              BETA=ONE-ALFA
              DO 3122 I=1,N
                T2(I)=(DXQ(I)+(FCA-ONE)*ALFA*DX(I))/BETA
3122          CONTINUE
              IF(NEW.EQ.1) THEN
                DO 3123 I=1,N
                  DXSAVE(I,1)=DX(I)
3123            CONTINUE
              ENDIF
              DO 3124 I=1,N
                DXSAVE(I,NEW+1)=T2(I)
                DX(I)=T2(I)
                T2(I)=T2(I)/XW(I)
3124          CONTINUE
            ENDIF
C           --------------------------------------------------------
C           3.2 Evaluation of scaled natural level function SUMX
C               scaled maximum error norm CONV
C               evaluation of (scaled) standard level function
C               DLEVF ( DLEVF only, if MPRMON.GE.2 )
C               and computation of ordinary Newton corrections DX(
C               N)
            CALL N2LVLS(N,T2,XW,F,DX,CONV,SUMX,DLEVF,MPRMON,NEW.EQ.0)
            DO 32 L1=1,N
              XA(L1)=X(L1)
32          CONTINUE
            SUMXA = SUMX
            DLEVXA = DSQRT(SUMXA/DBLE(FLOAT(N)))
            CONVA = CONV
            DXANRM = WNORM(N,DX,XW)
C           --------------------------------------------------------
C           3.3 A - priori estimate of damping factor FC
            QREDU = .FALSE.
            IF(NITER.NE.0.AND.NONLIN.NE.1)THEN
CWei;              IF(NEW.EQ.0.AND.(IRANK.EQ.N.OR.IRANKA.EQ.N).OR.
CWei;     *           QREPET)THEN
              IF(NEW.EQ.0.OR.QREPET)THEN
C               ------------------------------------------------------
C               3.3.1 Computation of the denominator of a-priori
C                     estimate
                FCDNM = ZERO
                DO 331 L1=1,N
                  FCDNM=FCDNM+((DX(L1)-DXQ(L1))/XW(L1))**2
331             CONTINUE
                IF(IRANK.NE.N)THEN
C                 --------------------------------------------
C                 3.3.2 Rank-deficient case (if previous rank
C                           was full) computation of the projected
C                       denominator of a-priori estimate
                  DO 332 L1=1,N
                    T1(L1)=DXQ(L1)/XW(L1)
332               CONTINUE
C                 Norm of projection of reduced component T1(N)
                  CALL N2PRJN(N,IRANK,DEL,T1,RWK(NRWLA+1),T2,QA,
     $                       IWK(NIWLA+2))
                  FCDNM = FCDNM-DEL
                ENDIF
                FCDNM = FCDNM*SUMX
C               ------------------------------------------------------
C               3.3.3 New damping factor
                IF(FCDNM.GT.FCNUMP*FCMIN2 .OR.
     $            (NONLIN.EQ.4 .AND. FCA**2*FCNUMP .LT. 4.0D0*FCDNM)) 
     $          THEN
                  DMYPRI = FCA*DSQRT(FCNUMP/FCDNM)
                  FCPRI = DMIN1(DMYPRI,ONE)
                  IF (NONLIN.EQ.4) FCPRI = DMIN1(HALF*DMYPRI,ONE)
                ELSE
                  FCPRI = ONE
C$Test-begin
                  DMYPRI = -1.0D0
C$Test-end
                ENDIF
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,33201) FCPRI, FC, FCA, DMYPRI, FCNUMP,
     $                               FCDNM
33201             FORMAT (/, ' +++ apriori estimate +++', /,
     $                   ' FCPRI  = ', D18.10, '  FC     = ', D18.10, /,
     $                   ' FCA    = ', D18.10, '  DMYPRI = ', D18.10, /,
     $                   ' FCNUMP = ', D18.10, '  FCDNM  = ', D18.10, /,
     $                       ' ++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
                FC = FCPRI
                IF ( IRANK.EQ.N .OR. IRANK.LE.MINRNK )
     $            FC=DMAX1(FC,FCMIN)
                IF (QBDAMP) THEN
                  FCBH = FCA*FCBND
                  IF (FC.GT.FCBH) THEN
                    FC = FCBH
                    IF (MPRMON.GE.4)
     $                WRITE(LUMON,*)' *** incr. rest. act. (a prio) ***'
                  ENDIF
                  FCBH = FCA/FCBND
                  IF (FC.LT.FCBH) THEN
                    FC = FCBH
                    IF (MPRMON.GE.4)
     $                WRITE(LUMON,*)' *** decr. rest. act. (a prio) ***'
                  ENDIF
                ENDIF
              ENDIF
              QREDU = FC.LT.FCMIN
            ENDIF
            QREPET = .FALSE.
CWEI
            IF (IORMON.GE.2) THEN
              SUMXA2=SUMXA1
              SUMXA1=SUMXA0
              SUMXA0=DLEVXA
              IF (SUMXA0.EQ.ZERO) SUMXA0=SMALL
C             Check convergence rates (linear and superlinear)
C             ICONV : Convergence indicator
C                     =0: No convergence indicated yet
C                     =1: Damping factor is 1.0d0
C                     =2: Superlinear convergence detected (alpha >=1.2)
C                     =3: Quadratic convergence detected (alpha > 1.8)
              FCMON = DMIN1(FC,FCMON)
              IF (FCMON.LT.ONE) THEN
                ICONV = 0
                ALPHAE = ZERO
              ENDIF
              IF (FCMON.EQ.ONE .AND. ICONV.EQ.0) ICONV=1
              IF (NITER.GE.1) THEN
                CLIN1 = CLIN0
                CLIN0 = SUMXA0/SUMXA1
              ENDIF
              IF (ICONV.GE.1.AND.NITER.GE.2) THEN
                ALPHAK = ALPHAE
                ALPHAE = ZERO
                IF (CLIN1.LE.0.95D0) ALPHAE = DLOG(CLIN0)/DLOG(CLIN1)
                IF (ALPHAK.NE.ZERO) ALPHAK =0.5D0*(ALPHAE+ALPHAK)
                ALPHAA = DMIN1(ALPHAK,ALPHAE)
                CALPHK = CALPHA
                CALPHA = ZERO
                IF (ALPHAE.NE.ZERO) CALPHA = SUMXA1/SUMXA2**ALPHAE
                SUMXTE = DSQRT(CALPHA*CALPHK)*SUMXA1**ALPHAK-SUMXA0
                IF (ALPHAA.GE.1.2D0 .AND. ICONV.EQ.1) ICONV = 2
                IF (ALPHAA.GT.1.8D0) ICONV = 3
                IF (MPRMON.GE.4)  WRITE(LUMON,32001) ICONV, ALPHAE, 
     $                              CALPHA, CLIN0, ALPHAK, SUMXTE
32001           FORMAT(' ** ICONV: ',I1,'  ALPHA: ',D9.2,
     $                '  CONST-ALPHA: ',D9.2,'  CONST-LIN: ',D9.2,' **',
     $                /,' **',11X,'ALPHA-POST: ',D9.2,' CHECK: ',D9.2,
     $                25X,'**')
                IF ( ICONV.GE.2 .AND. ALPHAA.LT.0.9D0 ) THEN
                   IF (IORMON.EQ.3) THEN
                     IERR = 4
                     GOTO 4299
                   ELSE
                     QMSTOP = .TRUE.
                   ENDIF 
                ENDIF
              ENDIF
            ENDIF
            FCMON = FC
C
            IF (MPRMON.GE.2) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N2PRV1(DLEVF,DLEVXA,FCKEEP,NITER,NEW,IRANK,MPRMON,
     $                    LUMON,QMIXIO,COND1)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            IF(.NOT.QREDU)THEN
C             --------------------------------------------------------
C             3.4 Save natural level for later computations of
C                 corrector and print iterate
              FCNUMK = SUMX
              NRED = 0
              QREP  = .FALSE.   
C             QREP = ITER .GT. ITMAX   or  QREP = ITER .GT. 0
C             Damping-factor reduction loop
C             ================================
C             DO (Until)
34            CONTINUE
C               ------------------------------------------------------
C               3.5 Preliminary new iterate
                DO 35 L1=1,N
                  X(L1)=XA(L1)+DX(L1)*FC
35              CONTINUE
C               -----------------------------------------------------
C               3.5.2 Exit, if problem is specified as being linear
C     Exit Repeat If ...
                IF( NONLIN.EQ.1 )THEN
                  IERR = 0
                  GOTO 4299
                ENDIF
C               ------------------------------------------------------
C               3.6.1 Computation of the residual vector
                IF (MPRTIM.NE.0) CALL MONON(1)
                CALL FCN(N,X,F,IFAIL)
                IF (MPRTIM.NE.0) CALL MONOFF(1)
                NFCN = NFCN+1
C     Exit, if ...
                IF (IFAIL.LT.0) THEN
                  IERR = 82
                  GOTO 4299
                ENDIF
                IF(IFAIL.EQ.1 .OR. IFAIL.EQ.2) THEN
                  IF (IFAIL.EQ.1) THEN
                    FCREDU = HALF
                  ELSE
                    FCREDU = F(1)
C     Exit, If ...
                    IF (FCREDU.LE.0 .OR. FCREDU.GE.1) THEN
                      IERR = 83
                      GOTO 4299
                    ENDIF
                  ENDIF
                  IF (MPRMON.GE.2) THEN
36101               FORMAT(8X,I2,' FCN could not be evaluated  ',
     $                     8X,F7.5,4X,I2,6X,I4)
                    WRITE(LUMON,36101)NITER,FC,NEW,IRANK
                  ENDIF
                  FCH = FC
                  FC = FCREDU*FC
                  IF (FCH.GT.FCMIN) FC = DMAX1(FC,FCMIN)
                  IF (QBDAMP) THEN
                    FCBH = FCH/FCBND
                    IF (FC.LT.FCBH) THEN
                      FC = FCBH
                      IF (MPRMON.GE.4) WRITE(LUMON,*)
     $                   ' *** decr. rest. act. (FCN redu.) ***'
                    ENDIF
                  ENDIF
                  IF (FC.LT.FCMIN) THEN
                    IERR = 3
                    GOTO 4299
                  ENDIF  
C     Break DO (Until) ...
                  GOTO 3109
                ENDIF
                DO 361 L1=1,N
                  T1(L1)=F(L1)*FW(L1)
361             CONTINUE
C               ------------------------------------------------------
C               3.6.2 Solution of linear (N,N)-system
                IF (QREPET) THEN
                  IWK(NIWLA) = 1
                ELSE
                  IWK(NIWLA) = 0
                ENDIF
                IF (MPRTIM.NE.0) CALL MONON(4)
                CALL N2SOLV(N,M1,N,ML,MU,A,QA,T1,T2,IRANK,IOPT,IFAIL,
     $                      LIWL,IWK(NIWLA),IDUMMY,LRWL,RWK(NRWLA),
     $                      IDUMMY)
                IF (MPRTIM.NE.0) CALL MONOFF(4)
C     Exit Repeat If ...
                IF(IFAIL.NE.0)  THEN
                  IERR = 81
                  GOTO 4299
                ENDIF
                IF(NEW.GT.0) THEN 
                  DO 3630 I=1,N
                    DXQ(I) = T2(I)*XWA(I)
3630              CONTINUE                   
                  DO 363 ILOOP=1,NEW 
                    SUM1=ZERO
                    SUM2=ZERO
                    DO 3631 I=1,N
                      SUM1=SUM1+(DXQ(I)*DXSAVE(I,ILOOP))/ XW(I)**2
                      SUM2=SUM2+(DXSAVE(I,ILOOP)/XW(I))**2
3631                CONTINUE
                    BETA=SUM1/SUM2
                    DO 3632 I=1,N
                      DXQ(I)=DXQ(I)+BETA*DXSAVE(I,ILOOP+1)
                      T2(I) = DXQ(I)/XW(I)
3632                CONTINUE
363               CONTINUE
                ENDIF
C               ------------------------------------------------------
C               3.6.3 Evaluation of scaled natural level function
C                     SUMX
C                     scaled maximum error norm CONV and evaluation
C                     of (scaled) standard level function DLEVF
                CALL N2LVLS(N,T2,XW,F,DXQ,CONV,SUMX,DLEVFN,MPRMON,
     $                      NEW.EQ.0)
                DXNRM = WNORM(N,DXQ,XW)
C               ------------------------------------------------------
C               3.6.4 Convergence test
C     Exit Repeat If ...
                IF ( DXNRM.LE.RTOL .AND. DXANRM.LE.RSMALL .AND. 
     $              FC.EQ.ONE ) THEN
                  IERR = 0
                  GOTO 4299
                ENDIF
C           
                FCA = FC
C               ----------------------------------------------------
C               3.6.5 Evaluation of reduced damping factor
                TH = FCA-ONE
                FCDNM = ZERO
                DO 39 L1=1,N
                  FCDNM=FCDNM+((DXQ(L1)+TH*DX(L1))/XW(L1))**2
39              CONTINUE
                IF (FCDNM.NE.ZERO) THEN
                  DMYCOR = FCA*FCA*HALF*DSQRT(FCNUMK/FCDNM)
                ELSE
                  DMYCOR = 1.0D+35
                ENDIF
                IF (NONLIN.LE.3) THEN
                  FCCOR = DMIN1(ONE,DMYCOR)
                ELSE
                  FCCOR = DMIN1(ONE,HALF*DMYCOR)
                ENDIF
C$Test-begin
                IF (MPRMON.GE.5) THEN
                  WRITE(LUMON,39001) FCCOR, FC, DMYCOR, FCNUMK,
     $                               FCDNM, FCA
39001             FORMAT (/, ' +++ corrector computation +++', /,
     $                   ' FCCOR  = ', D18.10, '  FC     = ', D18.10, /,
     $                   ' DMYCOR = ', D18.10, '  FCNUMK = ', D18.10, /,
     $                   ' FCDNM  = ', D18.10, '  FCA    = ', D18.10, /,
     $                       ' +++++++++++++++++++++++++++++', /)
                ENDIF
C$Test-end 
C               ------------------------------------------------------
C               3.7 Natural monotonicity test
                IF(SUMX.GT.SUMXA)THEN
C                 ----------------------------------------------------
C                 3.8 Output of iterate
                  IF(MPRMON.GE.3) THEN
                    IF (MPRTIM.NE.0) CALL MONON(5)
                    CALL N2PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                          NITER,MPRMON,LUMON,QMIXIO,'*')
                    IF (MPRTIM.NE.0) CALL MONOFF(5)
                  ENDIF
                  IF (QMSTOP) THEN
                    IERR = 4
                    GOTO 4299
                  ENDIF
                  FC = DMIN1(FCCOR,HALF*FC)
                  IF ((IRANK.EQ.N .OR. IRANK.EQ.MINRNK) .AND.
     $               FCA.GT.FCMIN)
     $               FC=DMAX1(FC,FCMIN)
                  IF (QBDAMP) THEN
                    FCBH = FCA/FCBND
                    IF (FC.LT.FCBH) THEN
                      FC = FCBH
                      IF (MPRMON.GE.4) WRITE(LUMON,*) 
     $                    ' *** decr. rest. act. (a post) ***'
                    ENDIF
                  ENDIF
CWEI
                  FCMON = FC
C
C$Test-begin
                  IF (MPRMON.GE.5) THEN
                    WRITE(LUMON,39002) FC
39002               FORMAT (/, ' +++ corrector setting 1 +++', /,
     $                      ' FC     = ', D18.10, /,
     $                         ' +++++++++++++++++++++++++++', /)
                  ENDIF
C$Test-end 
                  QREP = .TRUE.
                  NCORR = NCORR+1
                  NRED = NRED+1
C                 ----------------------------------------------------
C                  3.10 If damping factor is too small:
C                       Refresh Jacobian,if current Jacobian was computed
C                       by a Rank1-update, else reduce pseudo-rank
                  QREDU  = FC.LT.FCMIN.OR.NEW.GT.0.AND.NRED.GT.1
                ELSE
                  IF (.NOT.QREP .AND. FCCOR.GT.SIGMA2*FC) THEN
                    IF(MPRMON.GE.3) THEN
                      IF (MPRTIM.NE.0) CALL MONON(5)
                      CALL N2PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                            NITER,MPRMON,LUMON,QMIXIO,'+')
                      IF (MPRTIM.NE.0) CALL MONOFF(5)
                    ENDIF
                    FC = FCCOR
C$Test-begin
                    IF (MPRMON.GE.5) THEN
                      WRITE(LUMON,39003) FC
39003                 FORMAT (/, ' +++ corrector setting 2 +++', /,
     $                        ' FC     = ', D18.10, /,
     $                           ' +++++++++++++++++++++++++++', /)
                    ENDIF
C$Test-end 
                    QREP = .TRUE.
                  ELSE
                    QNEXT = .TRUE.
                  ENDIF
                ENDIF
3109          CONTINUE
              IF(.NOT.(QNEXT.OR.QREDU)) GOTO  34
C             UNTIL ( expression - negated above)
C             End of damping-factor reduction loop
C           =======================================
            ENDIF
            IF(QREDU)THEN
C             ------------------------------------------------------
C             3.11 Restore former values for repeting step
C                  step
              NREJR1 = NREJR1+1
              DO 3111 L1=1,N
                X(L1)=XA(L1)
3111          CONTINUE
              DO 3112 L1=1,N
                F(L1)=FA(L1)
3112          CONTINUE
              DO 3113 L1=1,N
                DXQ(L1)=DXQA(L1)
3113          CONTINUE
              IF(MPRMON.GE.2)THEN
31130           FORMAT(8X,I2,' Not accepted damping ',
     $                 'factor',8X,F7.5,4X,I2,6X,I4)
                WRITE(LUMON,31130)NITER,FC,NEW,IRANK
              ENDIF
              FC = FCKEEP
              FCA = FCK2
              IF(NITER.EQ.0)THEN
                FC = FCMIN
              ENDIF
              IF(NEW.GT.0)THEN
                QGENJ = .TRUE.
                QJCRFR = .TRUE.
                QREDU = .FALSE.
              ELSE
C               ------------------------------------------------
C               3.12 Pseudo-rank reduction
                QREPET = .TRUE.
                DO 42 L1=1,N
                  T1(L1)=QU(L1)
42              CONTINUE
                IRANK = IRANK-1
                IF(IRANK.LT.MINRNK)THEN
                  IERR = 3
                  GOTO 4299
                ENDIF
              ENDIF
            ENDIF
          IF(.NOT.(.NOT.QREDU)) GOTO  3
C         UNTIL ( expression - negated above)
C
C         End of pseudo-rank reduction loop
C         =================================
          IF (QNEXT) THEN
C           ------------------------------------------------------
C           4 Preparations to start the following iteration step
C           ------------------------------------------------------
C           4.1 Print values
            IF(MPRMON.GE.3) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N2PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,NITER+1,
     $                    MPRMON,LUMON,QMIXIO,'*')
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
C           Print the natural level of the current iterate and return
C           it in one-step mode
            SUMX = SUMXA
            IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N2SOUT(N,XA,2,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
              IF (MPRTIM.NE.0) CALL MONON(5)
              CALL N2SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
              IF (MPRTIM.NE.0) CALL MONOFF(5)
            ENDIF
            NITER = NITER+1
            DLEVF = DLEVFN
C     Exit Repeat If ...
            IF(NITER.GE.NITMAX)THEN
              IERR = 2
              GOTO 4299
            ENDIF
            FCKEEP = FC
C           ------------------------------------------------------
C           4.2 Return, if in one-step mode
C Exit Subroutine If ...
            IF (MODE.EQ.1) THEN
              IWK(18)=NIWLA-1
              IWK(19)=NRWLA-1
              IOPT(1)=1
              RETURN
            ENDIF
          ENDIF
        GOTO 2
C       End Repeat
4299    CONTINUE
C
C       End of main iteration loop
C       ==========================
C       ----------------------------------------------------------
C       9 Exits
C       ----------------------------------------------------------
C       9.1 Solution exit
        APREC = -1.0D0
C
        IF(IERR.EQ.0 .OR. IERR.EQ.4)THEN
          IF (NONLIN.NE.1) THEN
            IF ( IERR.EQ.0 ) THEN
              APREC = DSQRT(SUMX/DBLE(FLOAT(N)))
              DO 91 L1=1,N
                X(L1)=X(L1)+DXQ(L1)
91            CONTINUE
            ELSE 
              APREC = DSQRT(SUMXA/DBLE(FLOAT(N)))
              IF (ALPHAA.GT.ZERO .AND. IORMON.EQ.3) THEN
                DO 92 L1=1,N
                  X(L1)=X(L1)+DX(L1)
92              CONTINUE
              ENDIF
            ENDIF
            IF(IRANK.LT.N)THEN
              IERR = 1
            ENDIF
C           Print final monitor output
            IF(MPRMON.GE.2) THEN
              IF (IERR.EQ.0) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL N2PRV2(DLEVFN,DSQRT(SUMX/DBLE(FLOAT(N))),FC,
     $                      NITER+1,MPRMON,LUMON,QMIXIO,'*')
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ELSE IF (IORMON.EQ.3) THEN
                IF (MPRTIM.NE.0) CALL MONON(5)
                CALL N2PRV1(DLEVFN,DSQRT(SUMXA/DBLE(FLOAT(N))),FC,
     $                      NITER,NEW,IRANK,MPRMON,LUMON,QMIXIO,COND1)
                IF (MPRTIM.NE.0) CALL MONOFF(5)
              ENDIF
            ENDIF
            IF (  IORMON.GE.2 ) THEN
              IF ( ICONV.LE.1 .AND. ALPHAE .NE. ZERO 
     $                        .AND. ALPHAK .NE. ZERO ) IERR = 5
            ENDIF
C
            IF(MPRMON.GE.1.AND. IERR.NE.1) THEN
91001         FORMAT(///' Solution of nonlinear system ',
     $        'of equations obtained within ',I3,
     $        ' iteration steps',//,' Achieved relative accuracy',D10.3)
              IF (IERR.EQ.4) THEN
                WRITE(LUMON,91001) NITER,APREC
              ELSE
                WRITE(LUMON,91001) NITER+1,APREC
              ENDIF 
            ENDIF
          ELSE
            IF(MPRMON.GE.1) THEN
91002         FORMAT(///' Solution of linear system ',
     $        'of equations obtained by NLEQ2',//,' No estimate ',
     $        'available for the achieved relative accuracy')
                WRITE(LUMON,91002)
            ENDIF
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2 Fail exit messages
C       ----------------------------------------------------------
C       9.2.1 Termination at stationary point
        IF(IERR.EQ.1.AND.MPRERR.GE.1)THEN
92101     FORMAT(/,' Iteration terminates at stationary point',/)
          WRITE(LUERR,92101)
        ENDIF
C       ----------------------------------------------------------
C       9.2.2 Termination after more than NITMAX iterations
        IF(IERR.EQ.2.AND.MPRERR.GE.1)THEN
92201     FORMAT(/,' Iteration terminates after NITMAX ',
     $    '=',I3,'  Iteration steps')
          WRITE(LUERR,92201)NITMAX
        ENDIF
C       ----------------------------------------------------------
C       9.2.3 Newton method fails to converge
        IF(IERR.EQ.3.AND.MPRERR.GE.1)THEN
92301     FORMAT(/,' Newton method fails to converge')
          WRITE(LUERR,92301)
        ENDIF
CWEI
C       ----------------------------------------------------------
C       9.2.4.1 Superlinear convergence slowed down
        IF(IERR.EQ.4.AND.MPRERR.GE.1)THEN
92401     FORMAT(/,' Warning: Monotonicity test failed after ',A,
     $           ' convergence was already checked;',/,
     $    ' RTOL requirement may be too stringent',/)
92402     FORMAT(/,' Warning: ',A,' convergence slowed down;',/,
     $    ' RTOL requirement may be too stringent',/)
          IF (QMSTOP) THEN
            IF (ICONV.EQ.2) WRITE(LUERR,92401) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92401) 'quadratic'
          ELSE
            IF (ICONV.EQ.2) WRITE(LUERR,92402) 'superlinear'
            IF (ICONV.EQ.3) WRITE(LUERR,92402) 'quadratic'
          ENDIF
        ENDIF
C       ----------------------------------------------------------
C       9.2.4.2 Convergence criterion satisfied before superlinear
C               convergence has been established
        IF(IERR.EQ.5.AND.MPRERR.GE.1)THEN
92410     FORMAT(/,' Warning: No quadratic or superlinear convergence ',
     $           'established yet',/,
     $           10X,'your solution may perhaps may be less accurate ',
     $           /,10X,'as indicated by the standard error estimate')
          WRITE(LUERR,92410)
        ENDIF
C       ----------------------------------------------------------
C       9.2.5 Error exit due to linear solver routine N2FACT
        IF(IERR.EQ.80.AND.MPRERR.GE.1)THEN
92501     FORMAT(/,' Error ',I5,' signalled by linear solver N2FACT')
          WRITE(LUERR,92501) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.6 Error exit due to linear solver routine N2SOLV
        IF(IERR.EQ.81.AND.MPRERR.GE.1)THEN
92601     FORMAT(/,' Error ',I5,' signalled by linear solver N2SOLV')
          WRITE(LUERR,92601) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.7 Error exit due to fail of user function FCN
        IF(IERR.EQ.82.AND.MPRERR.GE.1)THEN
92701     FORMAT(/,' Error ',I5,' signalled by user function FCN')
          WRITE(LUERR,92701) IFAIL
        ENDIF
C       ----------------------------------------------------------
C       9.2.8 Error exit due to fail of user function JAC
        IF(IERR.EQ.83.AND.MPRERR.GE.1)THEN
92801     FORMAT(/,' Error ',I5,' signalled by user function JAC')
          WRITE(LUERR,92801) IFAIL
        ENDIF
        IF(IERR.GE.80.AND.IERR.LE.83) IWK(23) = IFAIL
        IF ((IERR.EQ.82.OR.IERR.EQ.83).AND.NITER.LE.1.AND.MPRERR.GE.1)
     $  THEN
          WRITE (LUERR,92810)
92810     FORMAT(' Try to find a better initial guess for the solution')
        ENDIF
C       ----------------------------------------------------------
C       9.3 Common exit
        IF (MPRERR.GE.3.AND.IERR.NE.0.AND.IERR.NE.4.AND.NONLIN.NE.1)
     $    THEN
93100     FORMAT(/,'    Achieved relative accuracy',D10.3,2X)
          WRITE(LUERR,93100)CONVA
          APREC = CONVA
        ENDIF
        IF(MPRMON.GE.1)THEN
93001     FORMAT(/,3X,'Subcondition ( 1,',I4,')',1X,D10.3,2X,
     $    /,3X,'Sensitivity ( 1,',I4,')',1X,D10.3,2X,/)
          WRITE(LUMON,93001)IRANK,COND1,IRANK,SENS1
        ENDIF
        RTOL = APREC
        SUMX = SUMXA
        IF(MPRSOL.GE.2.AND.NITER.NE.0) THEN
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL N2SOUT(N,XA,2,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ELSE IF(MPRSOL.GE.1.AND.NITER.EQ.0)THEN
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL N2SOUT(N,XA,1,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
        IF (IERR.NE.4) NITER = NITER+1
        DLEVF = DLEVFN
        IF(MPRSOL.GE.1)THEN
C         Print Solution or final iteration vector
          IF(IERR.EQ.0)THEN
             MODEFI = 3
          ELSE
             MODEFI = 4
          ENDIF
          IF (MPRTIM.NE.0) CALL MONON(5)
          CALL N2SOUT(N,X,MODEFI,IOPT,RWK,LRWK,IWK,LIWK,MPRSOL,LUSOL)
          IF (MPRTIM.NE.0) CALL MONOFF(5)
        ENDIF
C       Return the latest internal scaling to XSCAL
        DO 93 I=1,N
          XSCAL(I)=XW(I)
93      CONTINUE
C       End of exits
C       End of subroutine N2INT
      RETURN
      END
C
      SUBROUTINE N2SCAL(N,X,XA,XSCAL,XW,ISCAL,QINISC,IOPT,LRWK,RWK)
C*    Begin Prologue SCALE
      INTEGER N
      DOUBLE PRECISION X(N),XSCAL(N),XA(N),XW(N)
      INTEGER ISCAL
      LOGICAL QINISC
      INTEGER IOPT(50),LRWK
      DOUBLE PRECISION RWK(LRWK)
C     ------------------------------------------------------------
C
C*    Summary :
C    
C     S C A L E : To be used in connection with NLEQ2 .
C       Computation of the internal scaling vector XW used for the
C       Jacobian matrix, the iterate vector and it's related
C       vectors - especially for the solution of the linear system
C       and the computations of norms to avoid numerical overflow.
C
C*    Input parameters
C     ================
C
C     N         Int     Number of unknowns
C     X(N)      Dble    Current iterate
C     XA(N)     Dble    Previous iterate
C     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
C                       of interface routine NLEQ2
C     ISCAL     Int     Option ISCAL passed from IOPT-field
C                       (for details see description of IOPT-fields)
C     QINISC    Logical = .TRUE.  : Initial scaling
C                       = .FALSE. : Subsequent scaling
C     IOPT(50)  Int     Options array passed from NLEQ2 parameter list
C     LRWK      Int     Length of real workspace
C     RWK(LRWK) Dble    Real workspace (see description above)
C
C*    Output parameters
C     =================
C
C     XW(N)     Dble   Scaling vector computed by this routine
C                      All components must be positive. The follow-
C                      ing relationship between the original vector
C                      X and the scaled vector XSCAL holds:
C                      XSCAL(I) = X(I)/XW(I) for I=1,...N
C
C*    Subroutines called: D1MACH
C
C*    Machine constants used
C     ======================
C
      DOUBLE PRECISION SMALL
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL D1MACH
      INTRINSIC DABS,DMAX1
      DOUBLE PRECISION D1MACH,HALF
      PARAMETER (HALF=0.5D0)
      INTEGER MPRMON,LUMON
      SMALL  = D1MACH(6)
C*    Begin
      DO 1 L1=1,N
        IF (ISCAL.EQ.1) THEN
          XW(L1) = XSCAL(L1)
        ELSE
          XW(L1)=DMAX1(XSCAL(L1),(DABS(X(L1))+DABS(XA(L1)))*HALF,SMALL)
        ENDIF
1     CONTINUE
C$Test-Begin
      MPRMON = IOPT(13)
      IF (MPRMON.GE.6) THEN
        LUMON = IOPT(14)
        WRITE(LUMON,*) ' '
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) '      X-components   Scaling-components    '
        WRITE(LUMON,10) (X(L1), XW(L1), L1=1,N)
10      FORMAT('  ',D18.10,'  ',D18.10)
        WRITE(LUMON,*) ' ++++++++++++++++++++++++++++++++++++++++++'
        WRITE(LUMON,*) ' '
      ENDIF
C$Test-End
C     End of subroutine N2SCAL
      RETURN
      END
C
      SUBROUTINE N2SCRF(M,N,A,FW)
C*    Begin Prologue SCROWF
      INTEGER M,N
      DOUBLE PRECISION A(M,N),FW(M)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
C                   mode
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       M           Int    Number of rows of the matrix
C       N           Int    Number of columns of the matrix
C     * A(M,N)      Dble   Matrix to be scaled
C
C*    Output parameters
C     =================
C
C       FW(M)       Dble   Row scaling factors - FW(i) contains
C                          the factor by which the i-th row of A
C                          has been multiplied
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER J,K
      DOUBLE PRECISION S1,S2
C*    Begin
      DO 1 K=1,M
        S1=ZERO
        DO 2 J=1,N
          S2=DABS(A(K,J))
          IF (S2.GT.S1) S1=S2
2       CONTINUE
        IF (S1.GT.ZERO) THEN
          S1=ONE/S1
          FW(K)=S1
          DO 3 J=1,N
            A(K,J)=A(K,J)*S1
3         CONTINUE
        ELSE
          FW(K)=ONE
        ENDIF
1     CONTINUE
C     End of subroutine N1SCRF
      RETURN
      END
C
      SUBROUTINE N2FACT(N,LDA,LDAINV,ML,MU,A,AINV,COND,IRANK,IOPT,
     $IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK)
C*    Begin Prologue FACT
      INTEGER N,LDA,LDAINV,ML,MU
      DOUBLE PRECISION A(LDA,N),AINV(LDAINV,N)
      DOUBLE PRECISION COND
      INTEGER IRANK
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LAIWK,LRWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     F A C T : Call linear algebra subprogram for factorization of
C               a (N,N)-matrix with rank decision and casual compu-
C               tation of the rank deficient pseudo-inverse matrix
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C     N             Int    Order of the linear system
C     LDA           Int    Leading dimension of the matrix array A
C     LDAINV        Int    Leading dimension of the matrix array AINV
C     ML            Int    Lower bandwidth of the matrix (only for
C                          banded systems)
C     MU            Int    Upper bandwidth of the matrix (only for
C                          banded systems)
C   * A(LDA,N)      Dble   Matrix storage.
C   * COND          Dble   On Input, COND holds the maximum permitted
C                          subcondition for the prescribed rank
C                          On Output, it holds the estimated subcon-
C                          dition of A
C     IOPT(50)      Int    Option vector passed from NLEQ2
C
C*    Output parameters
C     =================
C
C     AINV(LDAINV,N) Dble   If matrix A is rank deficient, this array
C                           holds the pseudo-inverse of A
C     IFAIL          Int    Error indicator returned by this routine:
C                           = 0 matrix decomposition successfull
C                           =10 supplied (integer) workspace too small
C
C*    Workspace parameters
C     ====================
C
C     LIWK           Int    Length of integer workspace passed to this
C                           routine (In)
C     IWK(LIWK)      Int    Integer Workspace supplied for this routine
C     LAIWK          Int    Length of integer Workspace used by this 
C                           routine (out)       
C     LRWK           Int    Length of real workspace passed to this
C                           routine (In)                  
C     RWK(LRWK)      Dble   Real Workspace supplied for this routine
C     LARWK          Int    Length of real Workspace used by this 
C                           routine (out)
C
C*    Subroutines called:  DECCON
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL DECCON
      INTRINSIC DABS
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER IREPET,MCON
C*    Begin
      MPRERR=IOPT(11)
      LUERR=IOPT(12)
      LAIWK = N+2
      LARWK = 2*N+1
      IF (LIWK.GE.LAIWK.AND.LRWK.GE.LARWK) THEN
        MCON = 0
        IREPET = -IWK(1)
        IF (IREPET.EQ.0)  IWK(2) = MCON
        CALL DECCON(A,LDA,N,MCON,N,N,IWK(2),IRANK,COND,RWK(2),IWK(3),
     $              IREPET,AINV,RWK(N+2),IFAIL)
        IF (IFAIL.EQ.-2 .AND. MPRERR.GT.0) WRITE(LUERR,10001)
10001   FORMAT(1X,
     $       'DECCON failed to compute rank-deficient QR-decomposition',
     $        /)
        IF(IRANK.NE.0)THEN
          COND = DABS(RWK(2)/RWK(IRANK+1))
          RWK(1) = DABS(RWK(2))
        ELSE
          COND = ONE
          RWK(1) = ZERO
        ENDIF
      ELSE
        IFAIL = 10
10002   FORMAT(/,' Insuffient workspace for linear solver,',
     $         ' at least needed more needed : ',/,
     $         ' ',A,' workspace : ',I4)
        IF (LIWK.LT.LAIWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10002) 'Integer',LAIWK-LIWK
        IF (LRWK.LT.LARWK.AND.MPRERR.GT.0) 
     $    WRITE(LUERR,10002) 'Double',LARWK-LRWK
      ENDIF
      RETURN
      END
C
      SUBROUTINE N2SOLV(N,LDA,LDAINV,ML,MU,A,AINV,B,Z,IRANK,IOPT,
     $IFAIL,LIWK,IWK,LAIWK,LRWK,RWK,LARWK)
C*    Begin Prologue SOLVE
      INTEGER N,LDA,LDAINV,ML,MU
      DOUBLE PRECISION A(LDA,N),AINV(LDAINV,N)
      DOUBLE PRECISION B(N),Z(N)
      INTEGER IRANK
      INTEGER IOPT(50)
      INTEGER IFAIL
      INTEGER LIWK
      INTEGER IWK(LIWK)
      INTEGER LRWK,LAIWK
      DOUBLE PRECISION RWK(LRWK)
      INTEGER LARWK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L V E : Call linear algebra subprogram for solution of
C                 the linear system A*Z = B
C
C*    Parameters
C     ==========
C
C     N,LDA,LDAINV,ML,MU,A,AINV,IRANK,IOPT,IFAIL,LIWK,IWK,LAIWK,LRWK,
C     RWK,LARWK :
C                        See description for subroutine N2FACT.          
C     B          Dble    In:  Right hand side of the linear system
C                        Out: Rhs. transformed to the upper trian-
C                             gular part of the linear system
C     Z          Dble    Out: Solution of the linear system
C
C     Subroutines called: SOLCON
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL SOLCON
      INTEGER IREPET
C*    Begin
      MCON = 0
      IREPET = -IWK(1)
      CALL SOLCON(A,LDA,N,MCON,N,N,Z,B,IWK(2),IRANK,RWK(2),IWK(3),
     $            IREPET,AINV,RWK(N+2))
      IFAIL = 0
      RETURN
      END
C
      SUBROUTINE N2LVLS(N,DX1,XW,F,DXQ,CONV,SUMX,DLEVF,MPRMON,QDSCAL)
C*    Begin Prologue LEVELS
      INTEGER N,MPRMON
      DOUBLE PRECISION DX1(N),XW(N),F(N),DXQ(N)
      DOUBLE PRECISION CONV,SUMX,DLEVF
      LOGICAL QDSCAL
C     ------------------------------------------------------------
C
C*    Summary :
C
C     L E V E L S : To be used in connection with NLEQ2 .
C     provides descaled solution, error norm and level functions
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       DX1(N)         Dble   array containing the scaled Newton
C                             correction
C       XW(N)          Dble   Array containing the scaling values
C       F(N)           Dble   Array containing the residuum
C
C*    Output parameters
C     =================
C
C       DXQ(N)         Dble   Array containing the descaled Newton
C                             correction
C       CONV           Dble   Scaled maximum norm of the Newton
C                             correction
C       SUMX           Dble   Scaled natural level function value
C       DLEVF          Dble   Standard level function value (only
C                             if needed for print)
C       MPRMON         Int    Print information parameter (see
C                             driver routine NLEQ2 )
C       QDSCAL         Logic  .TRUE., if descaling of DX1 required,
C                             else .FALSE.
C
C     ------------------------------------------------------------
C*    End Prologue
      INTRINSIC DABS
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1
      DOUBLE PRECISION S1
C*    Begin
      IF (QDSCAL) THEN
C       ------------------------------------------------------------
C       1.2 Descaling of solution DX1 ( stored to DXQ )
        DO 12 L1=1,N
          DXQ(L1)=DX1(L1)*XW(L1)
12      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     2 Evaluation of scaled natural level function SUMX and
C       scaled maximum error norm CONV
      CONV = ZERO
      DO 20 L1=1,N
        S1 = DABS(DX1(L1))
        IF(S1.GT.CONV) CONV=S1
20    CONTINUE
      SUMX = ZERO
      DO 21 L1=1,N
        SUMX = SUMX+DX1(L1)**2
21    CONTINUE
C     ------------------------------------------------------------
C     3 Evaluation of (scaled) standard level function DLEVF (only
C       if needed for print)
      IF(MPRMON.GE.2)THEN
        DLEVF = ZERO
        DO 3 L1=1,N
          DLEVF = DLEVF+F(L1)**2
3       CONTINUE
        DLEVF = DSQRT(DLEVF/DBLE(FLOAT(N)))
      ENDIF
C     End of subroutine N2LVLS
      RETURN
      END
C
      SUBROUTINE N2JAC (FCN, N, LDA, X, FX, A, YSCAL, AJDEL, AJMIN,
     $                  NFCN, FU, IFAIL)
C* Begin Prologue N2JAC
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), AJDEL, AJMIN
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Evaluation of a dense Jacobian matrix using finite difference
C  approximation adapted for use in nonlinear systems solver NLEQ2
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   January 1991
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N2JAC is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading Dimension of array A
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
C  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DSIGN
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K
      DOUBLE PRECISION U, W
C
C* Begin
C
      IFAIL = 0
      DO 1 K = 1,N
         W = X(K)
         U = DSIGN(DMAX1(DABS(X(K)),AJMIN,YSCAL(K))*AJDEL, X(K))
         X(K) = W + U
C
         CALL FCN (N, X, FU, IFAIL)
         NFCN = NFCN + 1
         IF (IFAIL .NE. 0) GOTO 99
C
         X(K) = W
         DO 11 I = 1,N
            A(I,K) = (FU(I) - FX(I)) / U  
 11      CONTINUE
 1    CONTINUE
C
99    CONTINUE
      RETURN
C
C
C* End of N2JAC
C
      END
      SUBROUTINE N2JCF (FCN, N, LDA, X, FX, A, YSCAL, ETA, ETAMIN,
     $     ETAMAX, ETADIF, CONV, NFCN, FU, IFAIL)
C* Begin Prologue N2JCF
      EXTERNAL FCN
      INTEGER N, LDA
      DOUBLE PRECISION X(N), FX(N), A(LDA,N), YSCAL(N), ETA(N),
     $     ETAMIN, ETAMAX, ETADIF, CONV
      INTEGER NFCN
      DOUBLE PRECISION FU(N)
      INTEGER IFAIL
C
C  ---------------------------------------------------------------------
C
C* Title
C
C  Approximation of dense Jacobian matrix for nonlinear systems solver
C  NLEQ2 with feed-back control of discretization and rounding errors
C
C* Environment       Fortran 77
C                    Double Precision
C                    Sun 3/60, Sun OS
C* Latest Revision   May 1990
C
C
C* Parameter list description
C  --------------------------
C
C* External subroutines (to be supplied by the user)
C  -------------------------------------------------
C
C  FCN        Ext     FCN (N, X, FX, IFAIL)
C                     Subroutine in order to provide right-hand
C                     side of first-order differential equations
C    N        Int     Number of rows and columns of the Jacobian
C    X(N)     Dble    The current scaled iterates
C    FX(N)    Dble    Array containing FCN(X)
C    IFAIL    Int     Return code
C                     Whenever a negative value is returned by FCN
C                     routine N2JCF is terminated immediately.
C
C
C* Input parameters (* marks inout parameters)
C  ----------------
C
C  N          Int     Number of rows and columns of the Jacobian
C  LDA        Int     Leading dimension of A (LDA .GE. N)
C  X(N)       Dble    Array containing the current scaled
C                     iterate
C  FX(N)      Dble    Array containing FCN(X)
C  YSCAL(N)   Dble    Array containing the scaling factors
C  ETA(N)     Dble *  Array containing the scaled denominator
C                     differences
C  ETAMIN     Dble    Minimum allowed scaled denominator
C  ETAMAX     Dble    Maximum allowed scaled denominator
C  ETADIF     Dble    DSQRT (1.1*EPMACH)
C                     EPMACH = machine precision
C  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
C  NFCN       Int  *  FCN - evaluation count
C
C* Output parameters (* marks inout parameters)
C  -----------------
C
C  A(LDA,N)   Dble    Array to contain the approximated
C                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
C  ETA(N)     Dble *  Scaled denominator differences adjusted
C  NFCN       Int  *  FCN - evaluation count adjusted
C  IFAIL      Int     Return code non-zero if Jacobian could not
C                     be computed.
C
C* Workspace parameters
C  --------------------
C
C  FU(N)      Dble    Array to contain FCN(x+dx) for evaluation of
C                     the numerator differences
C
C* Called
C  ------
C
      INTRINSIC DABS, DMAX1, DMIN1, DSIGN, DSQRT
C
C* Constants
C  ---------
C
      DOUBLE PRECISION SMALL2, ZERO
      PARAMETER (SMALL2 = 0.1D0,
     $           ZERO   = 0.0D0)
C
C  ---------------------------------------------------------------------
C
C* End Prologue
C
C* Local variables
C  ---------------
C
      INTEGER I, K, IS
      DOUBLE PRECISION FHI, HG, U, SUMD, W
      LOGICAL QFINE
C
C* Begin
C
      DO 1 K = 1,N
         IS = 0
C        DO (Until)
 11         CONTINUE
            W = X(K)
            U = DSIGN (ETA(K)*YSCAL(K), X(K))
            X(K) = W + U
            CALL FCN (N, X, FU, IFAIL)
            NFCN = NFCN + 1
C           Exit, If ...
            IF (IFAIL .NE. 0) GOTO 99
            X(K) = W
            SUMD = ZERO
            DO 111 I = 1,N
               HG = DMAX1 (DABS (FX(I)), DABS (FU(I)))
               FHI = FU(I) - FX(I)
               IF (HG .NE. ZERO) SUMD = SUMD + (FHI/HG)**2
               A(I,K) = FHI / U
 111        CONTINUE
            SUMD = DSQRT (SUMD / DBLE(N))
            QFINE = .TRUE.
            IF (SUMD .NE. ZERO .AND. IS .EQ. 0)THEN
               ETA(K) = DMIN1 (ETAMAX,
     $              DMAX1 (ETAMIN, DSQRT (ETADIF / SUMD)*ETA(K)))
               IS = 1
               QFINE = CONV .LT. SMALL2 .OR. SUMD .GE. ETAMIN
            ENDIF
            IF (.NOT.(QFINE)) GOTO  11
C        UNTIL ( expression - negated above)
 1    CONTINUE
C
C     Exit from DO-loop
 99   CONTINUE
C
      RETURN
C
C* End of subroutine N2JCF
C
      END
      SUBROUTINE N2PRJN(N,IRANK,DEL,U,D,V,QE,PIVOT)
C*    Begin Prologue PRJCTN
      INTEGER IRANK,N
      INTEGER PIVOT(N)
      DOUBLE PRECISION DEL
      DOUBLE PRECISION QE(N,N)
      DOUBLE PRECISION U(N),D(N),V(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     P R J C T N :
C     To be used in connection with either DECOMP/SOLVE or 
C     DECCON/SOLCON .
C     Provides the projection to the appropriate subspace in case
C     of rank - reduction
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       N              Int    Number of parameters to be estimated
C       IRANK                 Pseudo rank of decomposed Jacobian
C                             matrix
C       U(N)           Dble   Scaled Newton correction
C       D(N)           Dble   Diagonal elements of upper
C                             triangular matrix
C       QE(N,N)        Dble   Part of pseudoinverse Jacobian
C                             matrix ( see QA of DECCON )
C       PIVOT(N)       Dble   Pivot vector resulting from matrix
C                             decomposition (DECCON)
C       V(N)           Dble   Real work array
C
C*    Output parameters
C     =================
C
C       DEL            Dble   Defekt
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,I,IRK1
      DOUBLE PRECISION S,SH
C*    Begin
      DO 1 I=1,N
        V(I)=U(PIVOT(I))
1     CONTINUE
      IRK1 = IRANK+1
      DEL = ZERO
      DO 2 I=IRK1,N
        SH = 0.0
        DO 21 L1=1,I-1
          SH = SH+QE(L1,I)*V(L1)
21      CONTINUE
        S =(V(I)-SH)/D(I)
        DEL = S*S+DEL
        V(I)=S
2     CONTINUE
C     End of subroutine N2PRJN
      RETURN
      END
C
      SUBROUTINE N2PRV1(DLEVF,DLEVX,FC,NITER,NEW,IRANK,MPRMON,LUMON,
     $                  QMIXIO,COND1)
C*    Begin Prologue N2PRV1
      DOUBLE PRECISION DLEVF,DLEVX,FC,COND1
      INTEGER NITER,NEW,IRANK,MPRMON,LUMON
      LOGICAL QMIXIO
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 2 P R V 1 : Printing of intermediate values (Type 1 routine)
C
C     Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N2INT
C     FC,NITER,NEW,IRANK,MPRMON,LUMON,COND1
C                  See parameter descr. of subroutine N2INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',20X,'New',6X,'Rank',
     $         8X,'Cond')
        IF (MPRMON.GE.3) WRITE(LUMON,2)
3       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.',3X,'New',
     $         6X,'Rank',8X,'Cond')
        IF (MPRMON.EQ.2) WRITE(LUMON,3)
      ENDIF
4     FORMAT(6X,I4,5X,D10.3,2X,4X,D10.3,17X,I2,6X,I4,2X,D10.3)
      IF (MPRMON.GE.3.OR.NITER.EQ.0) 
     $  WRITE(LUMON,4) NITER,DLEVF,DLEVX,NEW,IRANK,COND1
5     FORMAT(6X,I4,5X,D10.3,6X,D10.3,6X,F7.5,4X,I2,6X,I4,2X,D10.3)
      IF (MPRMON.EQ.2.AND.NITER.NE.0) 
     $  WRITE(LUMON,5) NITER,DLEVF,DLEVX,FC,NEW,IRANK,COND1
      IF(QMIXIO)THEN
6       FORMAT(2X,66('*'))
        WRITE(LUMON,6)
      ENDIF
C     End of subroutine N2PRV1
      RETURN
      END
C
      SUBROUTINE N2PRV2(DLEVF,DLEVX,FC,NITER,MPRMON,LUMON,QMIXIO,
     $                  CMARK)
C*    Begin Prologue N2PRV2
      DOUBLE PRECISION DLEVF,DLEVX,FC
      INTEGER NITER,MPRMON,LUMON
      LOGICAL QMIXIO
      CHARACTER*1 CMARK
C     ------------------------------------------------------------
C
C*    Summary :
C
C     N 2 P R V 2 : Printing of intermediate values (Type 2 routine)
C
C*    Parameters
C     ==========
C
C     DLEVF, DLEVX   See descr. of internal double variables of N2INT
C     FC,NITER,MPRMON,LUMON
C                  See parameter descr. of subroutine N2INT
C     QMIXIO Logical  = .TRUE.  , if LUMON.EQ.LUSOL
C                     = .FALSE. , if LUMON.NE.LUSOL
C     CMARK Char*1    Marker character to be printed before DLEVX
C
C     ------------------------------------------------------------
C*    End Prologue
C     Print Standard - and natural level, and damping
C     factor
      IF(QMIXIO)THEN
1       FORMAT(2X,66('*'))
        WRITE(LUMON,1)
2       FORMAT(8X,'It',7X,'Normf ',10X,'Normx ',8X,'Damp.Fct.')
        WRITE(LUMON,2)
      ENDIF
3     FORMAT(6X,I4,5X,D10.3,4X,A1,1X,D10.3,6X,F7.5)
      WRITE(LUMON,3)NITER,DLEVF,CMARK,DLEVX,FC
      IF(QMIXIO)THEN
4       FORMAT(2X,66('*'))
        WRITE(LUMON,4)
      ENDIF
C     End of subroutine N2PRV2
      RETURN
      END
C
      SUBROUTINE N2SOUT(N,X,MODE,IOPT,RWK,NRW,IWK,NIW,MPRINT,LUOUT)
C*    Begin Prologue SOLOUT
      INTEGER N
      DOUBLE PRECISION X(N)
      INTEGER NRW
      INTEGER MODE
      INTEGER IOPT(50)
      DOUBLE PRECISION RWK(NRW)
      INTEGER NIW
      INTEGER IWK(NIW)
      INTEGER MPRINT,LUOUT
C     ------------------------------------------------------------
C
C*    Summary :
C
C     S O L O U T : Printing of iterate (user customizable routine)
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     X(N)   Dble   iterate vector
C     MODE          =1 This routine is called before the first
C                      Newton iteration step
C                   =2 This routine is called with an intermedi-
C                      ate iterate X(N)
C                   =3 This is the last call with the solution
C                      vector X(N)
C                   =4 This is the last call with the final, but
C                      not solution vector X(N)
C     IOPT(50)  Int The option array as passed to the driver
C                   routine (elements 46 to 50 may be used
C                   for user options)
C     MPRINT    Int Solution print level 
C                   (see description of IOPT-field MPRINT)
C     LUOUT     Int the solution print unit 
C                   (see description of see IOPT-field LUSOL)
C
C
C*    Workspace parameters
C     ====================
C
C     NRW, RWK, NIW, IWK    see description in driver routine
C
C*    Use of IOPT by this routine
C     ===========================
C
C     Field 46:       =0 Standard output
C                     =1 GRAZIL suitable output
C
C     ------------------------------------------------------------
C*    End Prologue
      LOGICAL QGRAZ,QNORM
C*    Begin
      QNORM = IOPT(46).EQ.0
      QGRAZ = IOPT(46).EQ.1
      IF(QNORM) THEN
1        FORMAT('  ',A,' data:',/)
         IF (MODE.EQ.1) THEN
101        FORMAT('  Start data:',/,'  N =',I5,//,
     $            '  Format: iteration-number, (x(i),i=1,...N), ',
     $            'Normf , Normx ',/)
           WRITE(LUOUT,101) N
           WRITE(LUOUT,1) 'Initial'
         ELSE IF (MODE.EQ.3) THEN
           WRITE(LUOUT,1) 'Solution'
         ELSE IF (MODE.EQ.4) THEN
           WRITE(LUOUT,1) 'Final'
         ENDIF
2        FORMAT(' ',I5)
C        WRITE          NITER
         WRITE(LUOUT,2) IWK(1)
3        FORMAT((12X,3(D18.10,1X)))
         WRITE(LUOUT,3)(X(L1),L1=1,N)
C        WRITE          DLEVF,  DLEVX
         WRITE(LUOUT,3) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
         IF(MODE.EQ.1.AND.MPRINT.GE.2) THEN
           WRITE(LUOUT,1) 'Intermediate'
         ELSE IF(MODE.GE.3) THEN
           WRITE(LUOUT,1) 'End'
         ENDIF
      ENDIF
      IF(QGRAZ) THEN
        IF(MODE.EQ.1) THEN
10        FORMAT('&name com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,10)(I,I=1,N+2)
15        FORMAT('&def  com',I3.3,:,255(7(', com',I3.3,:),/))
          WRITE(LUOUT,15)(I,I=1,N+2)
16        FORMAT(6X,': X=1, Y=',I3)
          WRITE(LUOUT,16) N+2
        ENDIF
20      FORMAT('&data ',I5)
C        WRITE          NITER
        WRITE(LUOUT,20) IWK(1) 
21      FORMAT((6X,4(D18.10)))
        WRITE(LUOUT,21)(X(L1),L1=1,N)
C        WRITE          DLEVF,  DLEVX
        WRITE(LUOUT,21) RWK(19),DSQRT(RWK(18)/DBLE(FLOAT(N)))
        IF(MODE.GE.3) THEN
30        FORMAT('&wktype 3111',/,'&atext x ''iter''')
          WRITE(LUOUT,30)
35        FORMAT('&vars = com',I3.3,/,'&atext y ''x',I3,'''',
     $           /,'&run')
          WRITE(LUOUT,35) (I,I,I=1,N)
36        FORMAT('&vars = com',I3.3,/,'&atext y ''',A,'''',
     $           /,'&run')
          WRITE(LUOUT,36) N+1,'Normf ',N+2,'Normx '
C39       FORMAT('&stop')
C         WRITE(LUOUT,39)
        ENDIF
      ENDIF
C     End of subroutine N2SOUT
      RETURN
      END
C
      DOUBLE PRECISION FUNCTION WNORM(N,Z,XW)
      INTEGER N
      DOUBLE PRECISION Z(N), XW(N)
C     ------------------------------------------------------------
C
C*    Summary :
C
C     E N O R M : Return the norm to be used in exit (termination)
C                 criteria
C
C*    Input parameters
C     ================
C
C     N         Int Number of equations/unknowns
C     Z(N)     Dble  The vector, of which the norm is to be computed
C     XW(N)    Dble  The scaling values of Z(N)
C
C*    Output
C     ======
C
C     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
C                          to the scaling values in XW(N):
C                          = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
C
C     ------------------------------------------------------------
C*    End Prologue
      INTEGER I
      DOUBLE PRECISION S
C*    Begin
      S = 0.0D0
      DO 10 I=1,N
        S = S + ( Z(I)/XW(I) ) ** 2
10    CONTINUE
      WNORM = DSQRT( S / DBLE(FLOAT(N)) )
C     End of function WNORM
      RETURN
      END
C
C
C*    Group  Linear Solver subroutines (Code DECCON/SOLCON)
C
      SUBROUTINE DECCON(A,NROW,NCOL,MCON,M,N,IRANKC,IRANK,COND,D,PIVOT,
     *KRED,AH,V,IERR)
C*    Begin Prologue DECCON
      INTEGER IRANKC,IRANK,MCON
      INTEGER M,N,NROW,NCOL,KRED
      INTEGER PIVOT(NCOL)
      DOUBLE PRECISION COND
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION D(NCOL),V(NCOL)
      INTEGER IERR
C     ------------------------------------------------------------
C
C*  Title
C
C*    Deccon - Constrained Least Squares QR-Decomposition
C
C*  Written by        P. Deuflhard, U. Nowak, L. Weimann 
C*  Purpose           Solution of least squares problems, optionally
C                     with equality constraints.
C*  Method            Constrained Least Squares QR-Decomposition
C                     (see references below)
C*  Category          D9b1. -  Singular, overdetermined or
C                              underdetermined systems of linear 
C                              equations, generalized inverses. 
C                              Constrained Least Squares solution
C*  Keywords          Linear Least Square Problems, constrained, 
C                     QR-decomposition, pseudo inverse.
C*  Version           1.1
C*  Revision          March 1989
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C*  Copyright     (c) Konrad Zuse Zentrum fuer
C                     Informationstechnik Berlin
C                     Heilbronner Str. 10, D-1000 Berlin 31
C                     phone 0049+30+89604-0, 
C                     telefax 0049+30+89604-125
C*  Contact           Lutz Weimann 
C                     ZIB, Numerical Software Development 
C                     phone: 0049+30+89604-185 ;
C                     e-mail:
C                     RFC822 notation: weimann@sc.zib-berlin.de
C                     X.400: C=de;A=dbp;P=zib-berlin;OU=sc;S=Weimann
C
C*    References:
C     ===========
C
C       /1/ P.Deuflhard, V.Apostolescu:
C           An underrelaxed Gauss-Newton method for equality
C           constrained nonlinear least squares problems.
C           Lecture Notes Control Inform. Sci. vol. 7, p.
C           22-32 (1978)
C       /2/ P.Deuflhard, W.Sautter:
C           On rank-deficient pseudoinverses.
C           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
C    
C*    Related Programs:     SOLCON
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
C     Constrained QR-decomposition of (M,N)-system  with
C     computation of pseudoinverse in case of rank-defeciency .
C     First MCON rows belong to equality constraints.
C
C     ------------------------------------------------------------
C
C*    Parameters list description (* marks inout parameters)
C     ======================================================
C
C*    Input parameters
C     ================
C
C       A(NROW,NCOL) Dble   Array holding the (M,N)-Matrix to be 
C                           decomposed
C       NROW         Int    Declared number of rows of array A
C       NCOL         Int    Declared number of columns of array A and 
C                           rows and columns of array AH
C       MCON         Int    Number of equality constraints (MCON.LE.N)
C                           Internally reduced if equality constraints
C                           are linearly dependent
C       M            Int    Current number of rows of matrix A
C       N            Int    Current number of columns of matrix A
C     * IRANKC       Int    Prescribed maximum pseudo-rank of 
C                           constrained part of matrix A (IRANKC.LE.MCON)
C     * IRANK        Int    Prescribed maximum pseudo-rank of matrix A
C                           (IRANK.LE.N)
C     * COND         Dble   Permitted upper bound of DABS(D(1)/D(IRANK))
C                           and of DABS(D(IRANK+1)/D(IRANK))
C                           (subcondition numbers of A)
C       KRED         Int    Type of operation
C                           >=0  Householder triangularization
C                                (build up pseudo-inverse,if IRANK.LT.N)
C                           < 0  Reduction of pseudo-rank of matrix A, 
C                                skipping Householder triangularization, of
C                                 build-up new pseudo-inverse
C
C*    Output parameters
C     =================
C
C       A(NROW,NCOL)  Dble   Array holding the (M,N)-output consisting
C                            of the transformed matrix in the upper 
C                            right triangle and the performed House-
C                            holder transf. in the lower left triangle.
C     * IRANKC        Int    New pseudo-rank of constrained part of
C                            matrix A, determined so that
C                            DABS(D(1)/D(IRANKC))<COND
C     * IRANK         Int    New pseudo-rank of matrix A, determined
C                            so that DABS(D(IRANKC+1)/D(IRANK)) < COND
C       D(IRANK)      Dble   Diagonal elements of upper triangular matr.
C       PIVOT(N)      Int    Index vector storing permutation of columns
C                            due to pivoting
C     * COND          Dble   The maximum of the two sub-condition
C                            numbers belonging to the constrained
C                            part and the remaining part of A.
C                            (in case of rank reduction:
C                             sub-condition number which led to
C                             rank reduction)
C       AH(NCOL,NCOL) Dble   In case of rank-defect used to compute the
C                            pseudo-inverse (currently used will be an
C                            (N,N)-part of this array)
C       IERR          Int    Error indicator:
C                            = 0 : DECCON computations are successfull.
C                            =-2 : Numerically negative diagonal element
C                                  encountered during computation of
C                                  pseudo inverse - due to extremely bad
C                                  conditioned Matrix A. DECCON is
C                                  unable to continue rank-reduction.
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C*    Subroutines called: D1MACH
C
C*    Machine constants used
C     ======================
C
C     EPMACH = relative machine precision
      DOUBLE PRECISION EPMACH
C
C     ------------------------------------------------------------
C*    End Prologue
      EXTERNAL D1MACH
      INTRINSIC DABS,DSQRT
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION REDUCE
      PARAMETER (REDUCE=0.05D0)
      INTEGER L1
      DOUBLE PRECISION S1
      INTEGER I,II,IRANKH,IRK1,I1,J,JD,JJ,K,K1,LEVEL,MH
      DOUBLE PRECISION DD,D1MACH,H,HMAX,S,SH,SMALLS,T
C*    Begin
C     --------------------------------------------------------------
C     1 Initialization
      EPMACH = D1MACH(3)
      SMALLS = DSQRT(EPMACH*10.0D0)
      IF(IRANK.GT.N) IRANK = N
      IF(IRANK.GT.M) IRANK = M
C     --------------------------------------------------------------
C     1.1 Special case M=1 and N=1
      IF(M.EQ.1.AND.N.EQ.1)THEN
        PIVOT(1)=1
        D(1)=A(1,1)
        COND = ONE
        RETURN
      ENDIF
      IF(KRED.GE.0)THEN
C       ------------------------------------------------------------
C       1.1 Initialize pivot-array
        DO 11 J=1,N
          PIVOT(J)=J
11      CONTINUE
C       ------------------------------------------------------------
C       2. Constrained Householder triangularization    
        JD = 1
        IRANC1 = IRANKC + 1
        MH = MCON
        IRANKH = IRANKC
        IF(MH.EQ.0) THEN
          IRANKH = IRANK
          MH = M
        ENDIF
        IRK1 = IRANK
        DO 2  K=1,IRK1
2000      LEVEL = 1
          IF(K.NE.N)THEN
            K1 = K+1
C           DO (Until)
20          CONTINUE
              IF(JD.NE.0)THEN
                DO 201 J=K,N
                  S = ZERO
                  DO 2011 L1=K,MH
                    S = S+A(L1,J)**2
2011              CONTINUE
                  D(J)=S
201             CONTINUE
              ENDIF
C             ------------------------------------------------------
C             2.1 Column pivoting
              S1 = D(K)
              JJ = K
              DO 21 L1=K,N
                IF(D(L1).GT.S1) THEN
                  S1=D(L1)
                  JJ = L1
                ENDIF
21            CONTINUE
              H = D(JJ)
              IF(JD.EQ.1) HMAX = H/(COND* REDUCE)
              JD = 0
              IF(H.LT.HMAX) JD = 1
            IF(.NOT.(H.GE.HMAX)) GOTO 20
C           UNTIL ( expression - negated above)
            IF(JJ.NE.K)THEN
C             ------------------------------------------------------
C             2.2 Column interchange
              I = PIVOT(K)
              PIVOT(K)=PIVOT(JJ)
              PIVOT(JJ)=I
              D(JJ)=D(K)
              DO 221 L1=1,M
                S1=A(L1,JJ)
                A(L1,JJ)=A(L1,K)
                A(L1,K)=S1
221           CONTINUE
            ENDIF
          ENDIF
          H = ZERO
          DO 222 L1=K,MH
            H = H+A(L1,K)**2
222       CONTINUE
          T = DSQRT(H)
C         ----------------------------------------------------------
C         2.3.0 A-priori test on pseudo-rank
          IF  ( K.EQ.1 .OR. K.EQ.IRANC1 )  DD = T/COND
          IF(T.LE.DD .OR. K.GT.IRANKH)THEN
C           ------------------------------------------------------
C           2.3.1 Rank reduction
            IRANKH = K-1
            IF  (MH.NE.MCON)  THEN
              IRANK = IRANKH
              IF (IRANKC.EQ.IRANK) THEN
                LEVEL = 4
              ELSE
                LEVEL = 3
              ENDIF
            ELSE
              IRANKC = IRANKH
              IF  (IRANKC.NE.MCON) THEN
                MH = M
                IRANKH = IRANK
                JD = 1
                GOTO 2000
              ENDIF
            ENDIF
          ENDIF
          IF (LEVEL.EQ.1) THEN
C           ------------------------------------------------------
C           2.4 Householder step
            S = A(K,K)
            T = -DSIGN(T,S)
            D(K)=T
            A(K,K)=S-T
            IF(K.NE.N)THEN
              T = ONE/(H-S*T)
              DO 24 J=K1,N
                S = ZERO
                DO 241 L1=K,MH
                  S = S+A(L1,K)*A(L1,J)
241             CONTINUE
                S = S*T
                S1 =-S
                DO 242 L1=K,M
                  A(L1,J) = A(L1,J)+A(L1,K)*S1
242             CONTINUE
                D(J)=D(J)-A(K,J)**2
24            CONTINUE
              IF(K.EQ.IRANKC)THEN
                MH = M
                JD = 1
                IRANKH=IRANK
                IF (K.EQ.IRK1) LEVEL=3
              ENDIF
            ELSE
              LEVEL = 4
            ENDIF
          ENDIF
C       Exit Do 2 If ... 
          IF(LEVEL.GT.1) GOTO  2999
2       CONTINUE
C       ENDDO
2999    CONTINUE
      ELSE
        K = -1
        LEVEL = 3
      ENDIF
C     --------------------------------------------------------------
C     3 Rank-deficient pseudo-inverse
      IF(LEVEL.EQ.3)THEN
        IRK1 = IRANK+1
        DO 3 J=IRK1,N
          DO 31 II=1,IRANK
            I = IRK1-II
            S = A(I,J)
            IF(II.NE.1)THEN
              SH = ZERO
              DO 311 L1=I1,IRANK
                SH=SH+A(I,L1)*V(L1)
311           CONTINUE
              S = S-SH
            ENDIF
            I1 = I
            V(I)=S/D(I)
            AH(I,J)=V(I)
31        CONTINUE
          DO 32 I=IRK1,J
            S = ZERO
            DO 321 L1=1,I-1
              S = S+AH(L1,I)*V(L1)
321         CONTINUE
            IF(I.NE.J)THEN
              V(I)=-S/D(I)
              AH(I,J)=-V(I)
            ENDIF
32        CONTINUE
          IF (S.GT.-ONE) THEN
            D(J)=DSQRT(S+ONE)
          ELSE 
            IERR=-2
            GOTO 999
          ENDIF
3       CONTINUE
      ENDIF
C    --------------------------------------------------------------
C     9 Exit
      IF (IRANKC.NE.0) THEN
        SH = D(IRANKC)
        IF(SH.NE.ZERO) SH = DABS(D(1)/SH)
      ELSE
        SH=ZERO
      ENDIF
      IF (K.EQ.IRANK) T = D(IRANK)
      IF (IRANKC+1.LE.IRANK .AND. T.NE.ZERO) THEN
        S = DABS(D(IRANKC+1)/T)
      ELSE
        S = ZERO
      ENDIF
      IF (S.GT.SH) THEN
        COND=S
      ELSE
        COND=SH
      ENDIF
      IERR=0
999   RETURN
      END
      SUBROUTINE SOLCON(A,NROW,NCOL,MCON,M,N,X,B,IRANKC,IRANK,D,PIVOT,
     *KRED,AH,V)
C*    Begin Prologue SOLCON
      DOUBLE PRECISION A(NROW,NCOL),AH(NCOL,NCOL)
      DOUBLE PRECISION X(NCOL),B(NROW),D(NCOL),V(NCOL)
      INTEGER NROW,NCOL,MCON,M,N,IRANKC,IRANK,KRED
      INTEGER PIVOT(NCOL)
C     ------------------------------------------------------------
C    
C*    Summary
C     =======
C
C     Best constrained linear least squares solution of (M,N)-
C     system . First MCON rows comprise MCON equality constraints.
C     To be used in connection with subroutine DECCON
C     References:       See DECCON
C     Related Programs: DECCON
C    
C*    Parameters:
C     ===========
C
C*    Input parameters (* marks inout parameters)
C     ===========================================
C
C       A(M,N), NROW, NCOL, M, N, MCON, IRANKC, IRANK,
C       D(N), PIVOT(N), AH(N,N), KRED
C                           See input- respective output-parameters
C                           description of subroutine DECCON
C     * B(M)         Dble   Right-hand side of linear system, if
C                           KRED.GE.0
C                           Right-hand side of upper linear system,
C                           if KRED.LT.0
C
C*    Output parameters
C     =================
C
C       X(N)         Dble   Best LSQ-solution of linear system
C       B(M)         Dble   Right-hand of upper trigular system
C                           (transformed right-hand side of linear
C                            system)
C
C*    Workspace parameters
C     ====================
C
C       V(N)         Dble   Workspace array
C
C     ------------------------------------------------------------
C*    End Prologue
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      INTEGER L1,L2
      INTEGER I,II,I1,IRANC1,IRK1,J,JJ,J1,MH
      DOUBLE PRECISION S,SH
C*    Begin
C     ------------------------------------------------------------
C     1 Solution for pseudo-rank zero
      IF(IRANK.EQ.0)THEN
        DO 11 L1=1,N
          X(L1)=ZERO
11      CONTINUE
        RETURN
      ENDIF
      IF  ((IRANK.LE.IRANKC).AND.(IRANK.NE.N)) THEN
        IRANC1 = IRANKC + 1
        DO 12 L1=IRANC1,N
          V(L1) = ZERO
12      CONTINUE
      ENDIF
      IF(KRED.GE.0.AND.(M.NE.1.OR.N.NE.1))THEN
C       ----------------------------------------------------------
C       2 Constrained householder transformations of right-hand side
        MH = MCON
        IF(IRANKC.EQ.0) MH = M
        DO 21 J=1,IRANK
          S = ZERO
          DO 211 L1=J,MH
            S = S+A(L1,J)*B(L1)
211       CONTINUE
          S = S/(D(J)*A(J,J))
          DO 212 L1=J,M
            B(L1)=B(L1)+A(L1,J)*S
212       CONTINUE
          IF(J.EQ.IRANKC) MH = M
21      CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     3 Solution of upper triangular system
      IRK1 = IRANK+1
      DO 31 II=1,IRANK
        I = IRK1-II
        I1 = I+1
        S = B(I)
        IF(II.NE.1)THEN 
          SH = ZERO
          DO 311 L1=I1,IRANK 
            SH=SH+A(I,L1)*V(L1)
311       CONTINUE
          S = S-SH
        ENDIF
        V(I)=S/D(I)
31    CONTINUE
      IF((IRANK.NE.N).AND.(IRANK.NE.IRANKC))THEN
C       ----------------------------------------------------------
C       3.2 Computation of the best constrained least squares-
C           solution
        DO 321 J=IRK1,N
          S = ZERO
          DO 3211 L1=1,J-1
            S = S+AH(L1,J)*V(L1)
3211      CONTINUE
          V(J)=-S/D(J)
321     CONTINUE
        DO 322 JJ=1,N
          J = N-JJ+1
          S = ZERO
          IF(JJ.NE.1)THEN
            DO 3221 L1=J1,N
              S=S+AH(J,L1)*V(L1)
3221        CONTINUE
          ENDIF
          IF(JJ.NE.1.AND.J.LE.IRANK)THEN
            V(J)=V(J)-S
          ELSE
            J1 = J
            V(J)=-(S+V(J))/D(J)
          ENDIF
322     CONTINUE
      ENDIF
C     ------------------------------------------------------------
C     4 Back-permutation of solution components
      DO 4 L1=1,N
        L2 = PIVOT(L1)
        X(L2) = V(L1)
4     CONTINUE
      RETURN
      END
C*    End package
C
C
C*    Group  Time monitor package
C
C*    Begin Prologue
C     ------------------------------------------------------------
C
C*  Title
C    
C     Monitor - A package for making multiple time measurements and
C               summary statistics
C
C*  Written by        U. Nowak, L. Weimann 
C*  Version           1.0
C*  Revision          January 1991
C*  Latest Change     January 1991
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C*  Copyright     (c) Konrad Zuse Zentrum fuer
C                     Informationstechnik Berlin
C                     Heilbronner Str. 10, D-1000 Berlin 31
C                     phone 0049+30+89604-0, 
C                     telefax 0049+30+89604-125
C*  Contact           Lutz Weimann 
C                     ZIB, Numerical Software Development 
C                     phone: 0049+30+89604-185 ;
C                     e-mail: 
C                     RFC822 notation: weimann@sc.zib-berlin.de
C                     X.400: C=de;A=dbp;P=zib-berlin;OU=sc;S=Weimann
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
C    which may follow from aquisition or application of this code.
C
C* Software status 
C    This code is under care of ZIB and belongs to ZIB software class 1.
C
C  ---------------------------------------------------------------
C
C*    Summary:
C
C     Monitor is a package for generating time and summary statistics
C     about the execution of multiple program parts of any program.
C     Nested measurements of program parts are possible.
C     ------------------------------------------------------------
C
C*    Usage:
C
C     The usage of Monitor is naturally divided into three phases:
C     1. the initialization and setup phase before the start of
C        the program or subroutines package to be measured;
C     2. the run phase of the program to be measured;
C     3. the final evaluation call.
C
C     The phase 1 must start with exactly one call of the subroutine
C     MONINI, which passes a title string and a logical unit for
C     later statistics output and possible error messages to the
C     package. This call follows a number of calls of the subroutine
C     MONDEF, where each call associates an identification string
C     to a positive integer number, called the measurement index
C     - up to maxtab, where maxtab is a package constant. Multiple
C     measurement indices may be used for measurements of multiple
C     program parts. The index 0 must also be associated with some
C     identification string, and corresponds to all parts of the
C     measured program from the measurement start call till the final
C     evaluation call, which are not associated with specific positive
C     measurement indices. After all necessary MONDEF calls are done,
C     the measurements are started at begin of the program to be
C     measured by a parameterless call of MONSRT.
C     In phase 2, each program part to be measured must be immediately
C     preceeded by a call of the subroutine MONON with the associated 
C     measurement index, and must be immediately followed by a call of
C     the subroutine MONOFF with the same measurement index. Measure-
C     ments of nested program parts are possible, and nesting is allowed
C     up to the number mnest, where mnest is a package constant.
C     Calling MONOFF without a preceeding MONON call with the same 
C     measurement index, or calling one of these subroutines with a
C     measurement index not previously defined by a MONDEF call causes
C     an error stop of the program. 
C     Finally at the end of the program to be measured, the parameter-
C     less call of the subroutine MONEND closes all measurements and
C     prints the summary statistics.
C     As delivered, maxtab has a value 20 and mnest a value 10, but
C     both constants may be increased, if needed, to any possible
C     integer value, by simply changing it's values in the first 
C     parameter statement of the subroutine MONTOR below.
C
C*    Subroutines and their parameters:
C     =================================
C
C     MONINI(CIDENT,LUMON)  : Initialize Monitor
C       CIDENT  char*20  Identification string for the total measurement
C                        ( printed in summary )
C       LUMON   int      The logical unit for printing out the summary
C
C     MONDEF(MESIND,CIDMES) : Define one measurement index
C       MESIND  int      >=1 : measurement index for a specific part
C                        = 0 : measurement index for all remaining parts
C                              (i.e. not belonging to parts with 
C                               index >=1)
C       CIDMES  char*15  Identification string for the part associated
C                        with MESIND ( printed in summary )
C
C     MONSRT                : Start measurements
C       (no parameters)
C
C     MONON(MESIND)         : Start measurement of a specific part
C       MESIND  int      >=1 : measurement index for a specific part
C
C     MONOFF(MESIND)        : Stop measurement of a specific part
C       MESIND  int      >=1 : measurement index for a specific part
C
C     MONEND                : Finish measurements and print summary
C       (no parameters)
C
C
C*    Example:
C       Calling sequence:
C
C       CALL MONINI (' Example',6)
C       CALL MONDEF (0,'Solver')
C       CALL MONDEF (1,'User function')
C       CALL MONDEF (2,'User matrix')
C       CALL MONSRT ()
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C 1     CONTINUE      
C       ...
C       CALL MONON (2)
C       ...  user matrix code ...
C       CALL MONOFF(2)
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C       CALL MONON (1)
C       ...  user function code ...
C       CALL MONOFF(1)
C       ...
C       program to be measured (part without specific measurement index)
C       ...
C       IF (no termination) GOTO 1
C       ...
C       CALL MONEND ()
C     ------------------------------------------------------------
C 
      SUBROUTINE MONTOR
      PARAMETER(MAXTAB=20,MNEST=10)
      CHARACTER*15 NAME(MAXTAB),NAME0
      CHARACTER*20 TEXT 
      CHARACTER*(*) TEXTH 
      CHARACTER*(*) NAMEH   
      REAL SEC(MAXTAB),ASEC(MAXTAB),PC1(MAXTAB),PC2(MAXTAB)
      INTEGER COUNT(MAXTAB),INDACT(MNEST)
      LOGICAL QON(MAXTAB)
      INTEGER IOUNIT
C
      SAVE SEC,COUNT,ASEC,PC1,PC2,INDXO,TIME1,TIME0,MAXIND,NAME
      SAVE SEC0,NAME0,TEXT,MONI,QON,IONCNT,INDACT
C
C
      DATA MONI/6/ , INFO/1/ , IGRAPH/1/
C
      RETURN
C
C     initialize monitor
C
      ENTRY MONINI (TEXTH,IOUNIT)
C
      MONI=IOUNIT
      MAXIND=0
      TEXT=TEXTH
      DO 100 I=1,MAXTAB
        SEC(I)=0.
        ASEC(I)=0.
        COUNT(I)=0
        QON(I)=.FALSE.
100   CONTINUE
      DO 105 I=1,MNEST
        INDACT(I)=0
105   CONTINUE
C
      SEC0=0.
      IONCNT=0
      RETURN
C
C     define one monitor entry
C
      ENTRY MONDEF(INDX,NAMEH)
      IF(INDX.LT.0 .OR. INDX.GT.MAXTAB) GOTO 1190
      IF (INDX.GT.MAXIND) MAXIND=INDX
      IF (INDX.GT.0) THEN
        IF (COUNT(INDX).GT.0) GOTO 1290
      ENDIF
      IF (INDX.EQ.0) THEN
        NAME0 = NAMEH
      ELSE
        NAME(INDX) = NAMEH
      ENDIF
      RETURN
C
C     start monitor measurements
C 
      ENTRY MONSRT()
      CALL SECOND (TIME1)
C
C      if(igraph.gt.0) call gmini(maxind,name0,name)
C
      RETURN
C
C     start one measurement
C
      ENTRY MONON (INDX)
      IF(INDX.GT.MAXIND.OR.INDX.LE.0) GOTO 1010
      IF (QON(INDX)) GOTO 1030
      CALL SECOND(ASEC(INDX))
      QON(INDX)=.TRUE.
      IF (IONCNT.EQ.0) THEN
        SEC0=SEC0+ASEC(INDX)-TIME1
      ELSE
        INDXO=INDACT(IONCNT)
        SEC(INDXO)=SEC(INDXO)+ASEC(INDX)-ASEC(INDXO)
      ENDIF
      IONCNT=IONCNT+1
      INDACT(IONCNT)=INDX
      IF(INFO.GT.1) WRITE(MONI,*) ' enter',NAME(INDX),ASEC(INDX)
C
C      if(igraph.gt.0) call gmon(indx,sec0)
C
      RETURN
C
C     stop one measurement
C
      ENTRY MONOFF (INDX)
      IF(INDX.GT.MAXIND.OR.INDX.LE.0) GOTO 1010
      IF (.NOT. QON(INDX)) GOTO 1040
      CALL SECOND(TIME2)
      QON(INDX)=.FALSE.
      SEC(INDX)=SEC(INDX)+TIME2-ASEC(INDX)
      COUNT(INDX)=COUNT(INDX)+1
      IONCNT=IONCNT-1
      IF (IONCNT.EQ.0) THEN
        TIME1=TIME2
      ELSE
        ASEC(INDACT(IONCNT))=TIME2
      ENDIF
      IF(INFO.GT.1) WRITE(MONI,*) ' exit ',NAME(INDX),TIME2
C
C      if(igraph.gt.0) call gmoff(indx,sec(indx))
C
      RETURN
C
C     terminate monitor and print statistics
C
      ENTRY MONEND
      CALL SECOND (TIME0)
      SEC0=SEC0+TIME0-TIME1
C
      SUM=1.E-10
      DO 200 I=1,MAXIND
      SUM=SUM+SEC(I)
      IF(COUNT(I).LE.0) GOTO 200
      ASEC(I)=SEC(I)/FLOAT(COUNT(I))
200   CONTINUE
      SUM0=SUM+SEC0
C
      DO 250 I=1,MAXIND
      PC1(I)=100.*SEC(I)/SUM0
      PC2(I)=100.*SEC(I)/SUM
250   CONTINUE
      PC10=100.*SEC0/SUM0
      PC20=100.*SEC0/SUM
C
      WRITE(MONI,9500)
      WRITE(MONI,9510)
      WRITE(MONI,9505)
9500  FORMAT(///)
9510  FORMAT(1X,75('#'))
9505  FORMAT(' #',73X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9512) TEXT
9512  FORMAT(' #   Results from time monitor program for: ',A29,2X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9514) SUM0,SUM
9514  FORMAT(' #   Total time:',F11.3,5X,'Sum of parts:',F11.3,19X,'#')
      WRITE(MONI,9505)
      WRITE(MONI,9520)
9520  FORMAT(' #   ',2X,'name',12X,'calls',7X,'time',4X,'av-time',
     1       4X,'% total',6X,'% sum   #')
C
      I0=1
      WRITE(MONI,9550) NAME0,I0,SEC0,SEC0,PC10,PC20
9550  FORMAT(' #   ',A15,I8,F11.3,F11.4,F11.2,F11.2,'   #')
C
      DO 300 I=1,MAXIND
      WRITE(MONI,9550) NAME(I),COUNT(I),SEC(I),ASEC(I),PC1(I),PC2(I)
300   CONTINUE
C
C
      WRITE(MONI,9505)
      WRITE(MONI,9510)
      WRITE(MONI,9500)
C
C
C      IF(IGRAPH.GT.0) CALL GMEND
C
      RETURN
C
C  error exits
C
1010  CONTINUE
      WRITE(MONI,9010) INDX
9010  FORMAT(/,' error in subroutine monon or monoff',/,
     $         '   indx out of range    indx=',I4)
      GOTO 1111
C
1020  CONTINUE
      WRITE(MONI,9020) INDX
9020  FORMAT(/,' error in subroutine monoff',/,'   indx out of range',/,
     1         '   indx=',I4)
      GOTO 1111
C
1030  CONTINUE
      WRITE(MONI,9030) INDX
9030  FORMAT(/,' error in subroutine monon',/,
     $         '   measurement is already running for ',
     1         '   indx=',I4)
      GOTO 1111
C
1040  CONTINUE
      WRITE(MONI,9040) INDX
9040  FORMAT(/,' error in subroutine monoff',/,
     $         '   measurement has never been activated for ',
     1         '   indx=',I4)
      GOTO 1111
C
1190  CONTINUE
      WRITE(MONI,9190) MAXTAB,INDX
9190  FORMAT(/,' error in subroutine mondef',/,'   indx gt ',I4,/,
     1         '   indx=',I4)
      GOTO 1111
C
1290  CONTINUE
      WRITE(MONI,9290) INDX
9290  FORMAT(/,' error in subroutine mondef',/,'   indx = ',I4,
     1         '   already in use' )
      GOTO 1111
C
1111  STOP
C
C  end subroutine monitor
C
      END
C
C*    Group  Machine dependent subroutines and functions
C
      SUBROUTINE SECOND(TIME)
      REAL TIME
      TIME=0.0E0
      RETURN
      END
      DOUBLE PRECISION FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES AND MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST.
C
C       DATA SMALL(1),SMALL(2) /    1048576,          0 /
C       DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
C       DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
C       DATA DIVER(1),DIVER(2) / 1018167296,          0 /
C       DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES AND 8087-BASED
C     MICROS, SUCH AS THE IBM PC AND AT&T 6300, IN WHICH THE LEAST
C     SIGNIFICANT BYTE IS STORED FIRST.
C
      DATA SMALL(1),SMALL(2) /          0,    1048576 /
      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / 00564000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C
C      DATA LARGE(1) / 37757777777777777777B /
C      DATA LARGE(2) / 37157777777777777774B /
C
C      DATA RIGHT(1) / 15624000000000000000B /
C      DATA RIGHT(2) / 00000000000000000000B /
C
C      DATA DIVER(1) / 15634000000000000000B /
C      DATA DIVER(2) / 00000000000000000000B /
C
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / O"00564000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C
C      DATA LARGE(1) / O"37757777777777777777" /
C      DATA LARGE(2) / O"37157777777777777774" /
C
C      DATA RIGHT(1) / O"15624000000000000000" /
C      DATA RIGHT(2) / O"00000000000000000000" /
C
C      DATA DIVER(1) / O"15634000000000000000" /
C      DATA DIVER(2) / O"00000000000000000000" /
C
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /
C
C     MACHINE CONSTANTS FOR CONVEX C-1
C
C      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
C      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
C      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
C      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
C      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
C      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
C      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
C      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
C
C      DATA SMALL(1),SMALL(2) /        128,           0 /
C      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
C      DATA DIVER(1),DIVER(2) /       9472,           0 /
C      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER
C
C      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
C
C      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
C      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
C      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
C      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /
C
C/6S
      IF (I .LT. 1  .OR.  I .GT. 6) GOTO 999
      IF (I .LE. 5 ) THEN
        D1MACH = DMACH(I)
      ELSE IF (I .EQ. 6) THEN
C       D1MACH = DSQRT(DMACH(1)/DMACH(3))
        D1MACH = 4.94D-32
      ENDIF
      RETURN
  999 WRITE(6,1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
