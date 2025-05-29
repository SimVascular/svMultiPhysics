c     MUSC 2 Immediate postop
c     Ethan Kung  keo@ucsd.edu

c     Created by Mahdi Esmaily Moghadam 12-01-2010
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c This subroutine initializes the parameters, you need to read the
c comments and specify them carefuly
c--------------------------------------------------------------------
c This is an example for RCR boundary condition with parameters
c Rd, R, and C which are distal, and proximal resistance and capacitor.

      SUBROUTINE INITIALIZE(nTimeStep)
      USE COM
      IMPLICIT NONE
      INTENT(OUT) nTimeStep

      LOGICAl ierr
      INTEGER i, nTimeStep
      REAL(KIND=8), ALLOCATABLE :: tZeroX(:)
c
c********************************************************************
c For instance if pressure in 3D solver is in cgs and here mmHg
c pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
c here is mL/s qConv=1000=1D3. In the case both solver are using same
c unites you can set these two conversion coefficients to 1D0
      pConv = 1D0
      qConv = 1D0

c Only when all the surfaces of you model are coupled with NeumannSrfs
c you may set this to .TRUE.
      pCorr = .FALSE.
      qCorr = .FALSE.

c********************************************************************
c Block for your inputs

c These two value should match to that of defined in solver.inp
      nDirichletSrfs = 0
      nNeumannSrfs   = 2
c Number of unknowns that you need inside your lumped parameter network
      nUnknowns      = 3
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 10

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 3

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      ALLOCATE (tZeroX(nUnknowns), srfToXdPtr(nDirichletSrfs))       !
      ALLOCATE (srfToXPtr(nNeumannSrfs))                             !
      tZeroX = 0D0
c--------------------------------------------------------------------

      INCLUDE "initial_values_final.f"

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      INQUIRE (FILE='InitialData', EXIST=ierr)                       !
      IF (.NOT.ierr) THEN                                            !
c         PRINT *, 'Initializing unknowns in LPM'                     !
         OPEN (1, FILE='InitialData',STATUS='NEW',FORM='UNFORMATTED')!
         WRITE (1) 0D0                                               !
         DO i=1, nUnknowns                                           !
            WRITE (1) tZeroX(i)                                      !
         END DO                                                      !
         CLOSE(1)                                                    !
      END IF                                                         !
c--------------------------------------------------------------------

c Surface to X pointer: this defines which Unknown corresponds to which
c suface in "List of Neumann Surfaces" inside solver.inp
c For example, if you have "List of Neumann Surfaces= 2 8 4 ...."
c and you set "srfToXPtr = (/5,3,9,.../)"
C this means X(5) corresponds to surface 2, X(3) <-> surface 8,
c and X(9) <-> surface 4
c Also, Q(1) corresponds to the first element in srfToXPtr, in this
c example, Q(1)<=>X(5), Q(2)<=>X(3)
      srfToXPtr = (/1,2/)
c      srfToXdPtr = (/1/)

      END SUBROUTINE INITIALIZE

c####################################################################
c Here you should find the f_i=dx_i/dt, based on the following parameters:
c  current x_i:                   x(i)
c  Current time:                  t
c  Flowrates from 3D code:        Q(i)
c  Pressure from Dirichlet faces: P(i)

      SUBROUTINE FINDF(t, x, f, Q, P)
      USE COM
      IMPLICIT NONE
      INTENT(IN) t, Q
      INTENT(OUT) f

      REAL(KIND=8) t, x(nUnknowns), f(nUnknowns), Q(nNeumannSrfs),
     2   P(nDirichletSrfs)

!     These are the dumy variables
      REAL(KIND=8) Pivc,Ppwc,Pcon,Qout, Phigh, Pmax, tramp, R

      INCLUDE "parameters_final.f"

!********************************************************************
!     This is the only case that x is directly manipulated, ugly, but the
!     only way!!!

!     Implement a resistance filling network so that the ventricle pressure
!     Pv = Phigh - Q*R
!     This fills the sphere at an approximately constant rate

!     x(1): Pressure applied to endo surface. Actually, we will assign the pressure
!           by setting the offset instead
!     x(2): Pressure applied to endo cap surface (fake)
!     x(3): Accumulated volume during inflation. Not necessary, but for insight

!     Define a pressure ramp for Phigh (1e7 dyne/cm^2 in 0.11 second, then constant)
      Pmax = 1.0D+7
      tramp = 0.1
      IF (t .LT. tramp) THEN
            Phigh = Pmax * t / tramp
      ELSE
            Phigh = Pmax
      ENDIF
      R = 1.0D+5
      
!     The main body of equations
      f(1)  = 0 ! The endo surface pressure (we will set it with offset)
      f(2)  = 0 ! The endo cap surface pressure (fake)
      f(3)  = -Q(1) ! Integrate the flow rate to compute accumulated volume

!     Define offset pressure.
      offset(1) = Phigh + Q(1) * R ! Endo pressure
      offset(2) = Phigh + Q(2) * R ! Endo cap pressure

c     Assign the additional parameters to be printed
      Xprint(1) = offset(1) ! This is the actual value of pressure for this sim
      Xprint(2) = Q(1)      ! The flowrate/dVdt due to endo + cap surface
      Xprint(3) = Q(2)      ! The flowrate/dVdt due to endo cap surface    
      
      RETURN
      END SUBROUTINE FINDF


