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
c     Pressure sent to 3D is Pa
      pConv = 1D0
c     Flow from 3D is m^3/s
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
      nUnknowns      = 6
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 1

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 9

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

!     These are the dummy variables
      REAL(KIND=8)Cp, Cd, Rp, Lp, Rav, Rsl, Rd, Rmax1, Rmin1,
     2 Pref, C1, R2, Ct, tact, RLV, At, Pat, Qtot, Rmax2, Rmin2,
     3 Tc, Tm, steep

      INCLUDE "parameters_final.f"

!********************************************************************
      Tc   = 1
      Tm   = MOD(t,Tc)

!     Parameters used (from Martin's paper). All units are SI
      Lp = 1.3D05       ! kg / m^4
      Cp = 7.7D-09      ! m^4 s^2 / kg
      Cd = 8.7D-09      ! m^4 s^2 / kg
      Rp = 7.3D06       ! kg / m^4 s
      Rd = 1.0D08       ! kg / m^4 s
      Pref = 0          ! Pa

! great results
!      Rmax1 = 100
!      Rmin1 = 0.01
!      Rmax2 = 100
!      Rmin2 = 0.01
!      steep = 0.01

      Rmax1 = 1D9      ! kg / m^4 s;     
      Rmin1 = 1D6       ! kg / m^4 s
      Rmax2 = 1D9     ! kg / m^4 s;      
      Rmin2 = 1D6       ! kg / m^4 s
      steep = 1D-3     ! Pa

!     Prescribed atrium pressure with maximum = 14 mmHg = 14(133.322) Pa
      At = 8
      IF (Tm .GE. 0.0 .AND. Tm .LE. 0.2) THEN
         Pat = At * 5D-1*(1D0 - COS(2D0*pi*(Tm-0.0)/0.2)) + 6
      ELSE
         Pat = 6D0
      END IF
      Pat = Pat*133.322 ! Pa

      
!     Determine the resistance of valve based on the pressure at two sides
!      Rav = Rmin1 + (Rmax1-Rmin1)* EXP(x(1)-Pat)/(EXP(x(1)-Pat)+1) 
!      Rsl = Rmin2 + (Rmax2-Rmin2)* EXP(x(3)-x(1))/(EXP(x(3)-x(1))+1) 
      Rav = Rmin1 + (Rmax1-Rmin1)* (1+tanh((x(1)-Pat)/steep)) / 2
      Rsl = Rmin2 + (Rmax2-Rmin2)* (1+tanh((x(3)-x(1))/steep)) / 2
      
!     The main body of equations
!     x(1) = p_v      (ventricle pressure)
!     x(2) = p_v_cap  (ventricle pressure on cap (not used, just need it to compute resistance for cap))
!     x(3) = p_p      (proximal pressure)
!     x(4) = q_p      (proximal flowrate)
!     x(5) = p_d      (distal pressure)
!     x(6) = V_v      (ventricle volume)
      f(1) = 0
      f(2) = 0
      f(3) = 1/Cp * ((x(1)-x(3))/Rsl - x(4))
      f(4) = Rp/Lp * (-x(4) + (x(3)-x(5))/Rp)
      f(5) = 1/Cd * (x(4) - (x(5)-Pref)/Rd)
      f(6) = -Q(1)

!     These sets p_v = x(1) + offset(1), and likewise for p_v_cap
      offset(1) =  (Q(1) + Pat/Rav + x(3)/Rsl)/(1/Rav+1/Rsl)
      offset(2) =  (-Q(2) + Pat/Rav + x(3)/Rsl)/(1/Rav+1/Rsl)
      

c     Assign the additional parameters to be printed
      Xprint(1) = t                       ! Time
      Xprint(2) = offset(1)               ! p_v
      Xprint(3) = offset(2)               ! p_v_cap
      Xprint(4) = (Pat-offset(1))/Rav     ! Flow through AV valve
      Xprint(5) = (offset(1)-x(3))/Rsl    ! Flow through SL valve
      Xprint(6) = Rav                     ! AV valve resistance
      Xprint(7) = Rsl                     ! SL valve resistance
      Xprint(8) = Pat                     ! Atrial pressure (prescribed)
      Xprint(9) = Q(1)                    ! Flow rate into ventricle

      RETURN
      END SUBROUTINE FINDF