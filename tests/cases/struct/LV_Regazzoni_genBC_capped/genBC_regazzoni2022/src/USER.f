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
c     Pressure: 0D = mmHg, 3D = Pa
      pConv = 133.322
c     Flow: 0D = mL/s, 3D = m^3/s
      qConv = 1.0d-6   

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
      nUnknowns      = 14
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

      INCLUDE "initial_values.f"

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

c####################################################################
c     Activation function based on the double Hill equation
      SUBROUTINE activation_double_Hill(t, t_C, tau1, tau2, 
     2                                  T_HB, m1, m2, phi)
      REAL(KIND=8), INTENT(IN) :: t, t_C, tau1, tau2, T_HB, m1, m2
      REAL(KIND=8), INTENT(OUT) :: phi
      REAL(KIND=8) :: g1, g2, phi_max, tol
      
c     Compute the activation function based on the double Hill equation
      g1 = (MODULO(t - t_C, T_HB) / tau1)**m1
      g2 = (MODULO(t - t_C, T_HB) / tau2)**m2
      phi = g1 / (1.D0 + g1) * 1.D0 / (1.D0 + g2)

c     Normalize phi
      tol = 1.D-6
c     Parameters for ventriclular activation
      IF ((ABS(t_C - 0.207) < tol) .AND.
     2    (ABS(tau1 - 0.18551724137931036) < tol) .AND.
     2    (ABS(tau2 - 0.31172413793103454) < tol) .AND.
     2    (ABS(m1 - 1.32) < tol) .AND.
     2    (ABS(m2 - 27.4) < tol) .AND.
     2    (ABS(T_HB - 0.6896551724137931) < tol)) THEN
      phi_max = 0.6094935065418275
      ELSE IF ((ABS(t_C - 0.025) < tol) .AND.
     2    (ABS(tau1 - 0.07586206896551724) < tol) .AND.
     2    (ABS(tau2 - 0.12413793103448276) < tol) .AND.
     2    (ABS(m1 - 1.32) < tol) .AND.
     2    (ABS(m2 - 13.1) < tol) .AND.
     2    (ABS(T_HB - 0.6896551724137931) < tol)) THEN
      phi_max = 0.5584798709826506
      ELSE
      PRINT *, 'ERROR: Unknown parameters for activation_double_Hill'
      PRINT *, 't_C = ', t_C
      PRINT *, 'tau1 = ', tau1
      PRINT *, 'tau2 = ', tau2
      PRINT *, 'm1 = ', m1
      PRINT *, 'm2 = ', m2
      PRINT *, 'T_HB = ', T_HB

      STOP
      END IF

      phi = phi / phi_max

      END SUBROUTINE activation_double_Hill
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

!     Declare local variables
      REAL(KIND=8) p_LA, p_RA, p_LV, p_RV, p_LV_cap, p_RV_cap
      REAL(KIND=8) p_AR_SYS, p_VEN_SYS, p_AR_PUL, p_VEN_PUL
      REAL(KIND=8) Q_LV, Q_RV, Q_LV_cap, Q_RV_cap 
      REAL(KIND=8) Q_AR_SYS, Q_VEN_SYS, Q_AR_PUL, Q_VEN_PUL
      REAL(KIND=8) Q_MV, Q_AV, Q_TV, Q_PV
      REAL(KIND=8) V_LA, V_RA, V_LV, V_RV
      REAL(KIND=8) A_LA, A_RA, A_RV
      REAL(KIND=8) R_MV, R_AV, R_TV, R_PV


      INCLUDE "parameters.f"

!********************************************************************

!     Variable definitions
!     x(1) = p_LV       (left ventricle pressure)
!     x(2) = p_LV_cap   (left ventricle pressure on cap (not used, just need it to compute resistance for cap))
!     x(3) = V_RV       (right ventricle volume)
!     x(4) = V_LA       (left atrium volume)
!     x(5) = V_RA       (right atrium volume)
!     x(6) = p_AR_SYS   (arterial systemic pressure)
!     x(7) = p_VEN_SYS  (venous systemic pressure)
!     x(8) = p_AR_PUL   (arterial pulmonic pressure)
!     x(9)= p_VEN_PUL  (venous pulmonic pressure)
!     x(10)= Q_AR_SYS   (arterial systemic flow rate)
!     x(11)= Q_VEN_SYS  (venous systemic flow rate)
!     x(12)= Q_AR_PUL   (arterial pulmonic flow rate)
!     x(13)= Q_VEN_PUL  (venous pulmonic flow rate)
!
!     Auxillary variables
!     x(14)= V_LV       (left ventricle volume)
!     
!     Note:
!     Q(1) = LV flow rate (out of LV)
!     Q(2) = LV cap flow rate
!    --------------- EQUATIONS --------------
      
!     Unpack variables from x for convenience
      p_LV = x(1)             ! mmHg; Left ventricle pressure
      p_LV_cap = x(2)         ! mmHg; Left ventricle pressure on cap (not used, just need it to compute resistance for cap)
      V_RV = x(3)             ! mL; Right ventricle volume
      V_LA = x(4)             ! mL; Left atrium volume
      V_RA = x(5)             ! mL; Right atrium volume
      p_AR_SYS = x(6)         ! mmHg; Arterial systemic pressure
      p_VEN_SYS = x(7)        ! mmHg; Venous systemic pressure
      p_AR_PUL = x(8)         ! mmHg; Arterial pulmonic pressure
      p_VEN_PUL = x(9)        ! mmHg; Venous pulmonic pressure
      Q_AR_SYS = x(10)        ! mL/s; Arterial systemic flow rate
      Q_VEN_SYS = x(11)       ! mL/s; Venous systemic flow rate
      Q_AR_PUL = x(12)        ! mL/s; Arterial pulmonic flow rate
      Q_VEN_PUL = x(13)       ! mL/s; Venous pulmonic flow rate
      V_LV = x(14)            ! mL; Left ventricle volume    

!     Unpack flowrates from Q for convenience
      Q_LV = Q(1)       ! LV flow rate (out of LV)
      Q_LV_cap = Q(2)   ! LV cap flow rate

!     Compute the atrial time-varying elastances (E_LA and E_RA)
!      CALL time_varying_elastance(t, E_LA_act, E_LA_pas, time_C_LA,  
!     2                              T_C_LA, T_R_LA, T_HB, E_LA)
!      CALL time_varying_elastance(t, E_RA_act, E_RA_pas, time_C_RA,
!     2                              T_C_RA, T_R_RA, T_HB, E_RA)
      CALL activation_double_Hill(t, time_C_RV, tau1_RV, tau2_RV, T_HB,
     2                               m1_RV, m2_RV, A_RV)
      CALL activation_double_Hill(t, time_C_LA, tau1_LA, tau2_LA, T_HB, 
     2                               m1_LA, m2_LA, A_LA)
      CALL activation_double_Hill(t, time_C_RA, tau1_RA, tau2_RA, T_HB,
     2                               m1_RA, m2_RA, A_RA)

!     Compute atrial pressures
!      p_LA = E_LA * (V_LA - V0_LA)
!      p_RA = E_RA * (V_RA - V0_RA)
      p_RV = (E_RV_pas + E_RV_act * A_RV) * (V_RV - V0_RV)
      p_LA = (E_LA_pas + E_LA_act * A_LA) * (V_LA - V0_LA)
      p_RA = (E_RA_pas + E_RA_act * A_RA) * (V_RA - V0_RA)

!     Determine the resistance of valve based on the pressure at two sides
!     of the valve. The resistance is a function of the pressure difference
      R_MV = merge(R_min, R_max, p_LA >= p_LV)
      R_AV = merge(R_min, R_max, p_LV >= p_AR_SYS)
      R_TV = merge(R_min, R_max, p_RA >= p_RV)
      R_PV = merge(R_min, R_max, p_RV >= p_AR_PUL)

!     Compute the flow rates through the valves
      Q_MV = (p_LA - p_LV) / R_MV
      Q_AV = (p_LV - p_AR_SYS) / R_AV
      Q_TV = (p_RA - p_RV) / R_TV
      Q_PV = (p_RV - p_AR_PUL) / R_PV

!     The main body of equations
      f(1) = 0    ! LV pressure is determined by algebraic equation offset(1) 
      f(2) = 0    ! LV cap pressure is determined by algebraic equation offset(2)
      f(3) = Q_TV - Q_PV        ! RV volume
      f(4) = Q_VEN_PUL - Q_MV   ! LA volume
      f(5) = Q_VEN_SYS - Q_TV   ! RA volume
      f(6) = 1/C_AR_SYS * (Q_AV - Q_AR_SYS)   ! Arterial systemic pressure
      f(7) = 1/C_VEN_SYS * (Q_AR_SYS - Q_VEN_SYS)   ! Venous systemic pressure
      f(8) = 1/C_AR_PUL * (Q_PV - Q_AR_PUL)   ! Arterial pulmonic pressure
      f(9)= 1/C_VEN_PUL * (Q_AR_PUL - Q_VEN_PUL)   ! Venous pulmonic pressure
      f(10) = R_AR_SYS/L_AR_SYS * 
     2       (-Q_AR_SYS - (p_VEN_SYS - p_AR_SYS) / R_AR_SYS)   ! Arterial systemic flow rate
      f(11) = R_VEN_SYS/L_VEN_SYS * 
     2       (-Q_VEN_SYS - (p_RA - p_VEN_SYS) / R_VEN_SYS)   ! Venous systemic flow rate
      f(12) = R_AR_PUL/L_AR_PUL * 
     2       (-Q_AR_PUL - (p_VEN_PUL - p_AR_PUL) / R_AR_PUL)   ! Arterial pulmonic flow rate
      f(13) = R_VEN_PUL/L_VEN_PUL * 
     2       (-Q_VEN_PUL - (p_LA - p_VEN_PUL) / R_VEN_PUL)   ! Venous pulmonic flow rate

      f(14) = -Q_LV   ! LV volume
      

!     These set p_lv = x(1) + offset(1), and likewise for p_lv_cap, p_rv, and p_rv_cap 
!     offset(3,4) only need to contain Q(3,4), so that the cap resistance is
!     non-zero. The value of the cap resistance is not used at all.
!     
!     NOTE: These equations are implicit in the ventricular pressures, since the
!     valve resistances depend on the ventricular pressures. However, we are solving
!     the equations explicitly by computing the valve resistances using the values 
!     of ventricular pressure from the previous 3D time step. There is some error
!     associated with this splitting, but this also prevents valves from switching
!     between 3D Newton iterations, causing cycling of Newton's method.
      offset(1) =  (Q_LV + p_LA/R_MV+ p_AR_SYS/R_AV)/(1/R_MV+1/R_AV) ! LV pressure
      offset(2) = (Q_LV_cap + p_LA/R_MV+ p_AR_SYS/R_AV)/(1/R_MV+1/R_AV) ! Fake LV cap pressure

c     Assign the additional parameters to be printed
      Xprint(1) = t                       ! Time
      Xprint(2) = R_MV                    ! Mitral valve resistance
      Xprint(3) = R_AV                    ! Aortic valve resistance
      Xprint(4) = R_TV                    ! Tricuspid valve resistance
      Xprint(5) = R_PV                    ! Pulmonary valve resistance
      Xprint(6) = Q_LV                    ! Flow rate out of LV
      Xprint(7) = A_RV                    ! RV activation
      Xprint(8) = A_LA                    ! LA activation
      Xprint(9) = A_RA                    ! RA activation
      RETURN
      END SUBROUTINE FINDF