c     Created by Mahdi Esmaily Moghadam 05-25-2011
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c Here data are received from Phasta and it setup the data required for
c integration of ODE's inside FINDX
      PROGRAM GenBC
      USE COM
      USE, INTRINSIC :: IEEE_ARITHMETIC
      IMPLICIT NONE

      INTEGER i, n, nTimeStep
      REAL(KIND=8) t, dt, tFinal, pMin, temp
      CHARACTER(LEN=2048) string
      CHARACTER(LEN=32) sTmp
      CHARACTER flag

      REAL(KIND=8), ALLOCATABLE :: Qi(:), Qf(:), Pi(:), Pf(:), X(:),
     2   Xo(:), f(:,:)

c     flag = I : Initializing
c     flag = T : Iteration loop
c     flag = L : Last iteration
c     flag = D : Derivative
c
c   Here is the period of time that you should integrate for
c   each of these flags:
c
c Flags              I                                            D&T&L
c                    ^                                              ^
c 3D code time step: N.............................................N+1
c 0D code time step: 1......................................nTimeStep+1
c Flowrates:         Qi............................ ................Qf
c Time, t:         tInitial..............................tInitial+tFinal


      CALL INITIALIZE(nTimeStep)

c********************************************************************
c Block for reading the data from phasta
      OPEN (1, FILE='GenBC.int', STATUS='OLD', FORM='UNFORMATTED');
      READ (1) flag
      READ (1) tFinal
      READ (1) i
      IF (i .NE. nDirichletSrfs) THEN
         PRINT *, 'Error: Number of Dirichlet Surfaces from Phasta is:',
     2      i
         PRINT *, 'While nDirichletSrfs is equal to:', nDirichletSrfs
         PRINT *
         PRINT *, 'Number of Dirichlet surfaces defined in solver.inp',
     2      ' should match with nDirichletSrfs defined in USER.f'
         STOP
      END IF
      READ (1) i
      IF (i .NE. nNeumannSrfs) THEN
         PRINT *, 'Error: Number of Neumann Surfaces from Phasta is:',
     2      i
         PRINT *, 'While nNeumannSrfs is equal to:', nNeumannSrfs
         PRINT *
         PRINT *, 'Number of Neumann surfaces defined in solver.inp',
     2      ' should match with nNeumannSrfs defined in USER.f'
         STOP
      END IF
      IF (nDirichletSrfs .GT. 0) THEN
         IF (qCorr .OR. pCorr) THEN
            PRINT *, 'You should only use P/Q correction when all',
     2         ' the surfaces are Neumann surfaces'
            STOP
         END IF
      END IF

      ALLOCATE (Pi(nDirichletSrfs), Pf(nDirichletSrfs),
     2   PDirichlet(nDirichletSrfs,4))
      ALLOCATE (Qi(nNeumannSrfs), Qf(nNeumannSrfs),
     2   QNeumann(nNeumannSrfs,4))

      DO i=1, nDirichletSrfs
         READ(1) Pi(i)
         READ(1) Pf(i)
      END DO
      Pi = Pi/pConv
      Pf = Pf/pConv

      DO i=1,nNeumannSrfs
         READ (1) Qi(i)
         READ (1) Qf(i)
c        Print to console and file pressure from genBC to svFSI for debugging
         PRINT*, 'BC i = ', i, 
     2      'Flowrates to genBC: ', Qi(i), Qf(i)
      END DO
      CLOSE(1)
      Qi = Qi/qConv
      Qf = Qf/qConv

c********************************************************************
c Block for initializing the unknowns

      IF (flag .EQ. 'I') THEN
         nTimeStep = 0
      ELSE
         dt = tFinal/REAL(nTimeStep,8)
      END IF
      ALLOCATE (X(nUnknowns), Xo(nUnknowns), f(nUnknowns,4),
     2   offset(nUnknowns), Xprint(nXprint))
      Xprint = 0D0

      OPEN (1, FILE='InitialData', STATUS='OLD', FORM='UNFORMATTED')
      READ (1) t
      DO i=1,nUnknowns
          READ (1) Xo(i)
      END DO
      CLOSE (1)

c********************************************************************
c Setting up the system of equations
      offset = 0D0
      DO n=1, nTimeStep
         DO i=1, 4
            temp = (REAL(n-1,8) + REAL(i-1,8)/3D0)/REAL(nTimeStep,8)

            QNeumann(:,i)   = Qi + (Qf-Qi)*temp
            PDirichlet(:,i) = Pi + (Pf-Pi)*temp

            IF (qCorr) THEN
               temp = SUM(QNeumann(:,i))/REAL(nNeumannSrfs,8)
               QNeumann(:,i) = QNeumann(:,i) - temp
            END IF
         END DO

         X = Xo
         CALL FINDF (t, X, f(:,1), QNeumann(:,1),
     2      PDirichlet(:,1))
         X = Xo + dt*f(:,1)/3D0

         CALL FINDF (t+dt/3D0, X, f(:,2), QNeumann(:,2),
     2      PDirichlet(:,2))
         X = Xo - dt*f(:,1)/3D0 + dt*f(:,2)

         CALL FINDF (t+dt*2D0/3D0, X, f(:,3), QNeumann(:,3),
     2      PDirichlet(:,3))
         X = Xo + dt*f(:,1) - dt*f(:,2) + dt*f(:,3)

         CALL FINDF (t+dt, X, f(:,4), QNeumann(:,4),
     2      PDirichlet(:,4))

         f(:,1) = (f(:,1) + 3D0*f(:,2) + 3D0*f(:,3) + f(:,4))/8D0
         Xo = Xo + dt*f(:,1)
         t = t + dt
      END DO

c********************************************************************
c Time to write the results
      X = Xo
      IF (pCorr .AND. flag.NE.'D') THEN
         pMin = X(srfToXPtr(1))
         DO i=2, nNeumannSrfs
            IF (X(srfToXPtr(i)) .LT. pMin) THEN
               pMin = X(srfToXPtr(i))
            END IF
         END DO
      ELSE
         pMin = 0D0
      END IF

c Writing nDirichlet flowrates here
      OPEN (1, FILE='GenBC.int', STATUS='OLD', FORM='UNFORMATTED')
      DO i=1, nDirichletSrfs
         IF(IEEE_IS_NAN(X(srfToXdPtr(i)))) THEN
            PRINT*, 'Error! NAN encountered..'
            STOP
         END IF
         WRITE (1) X(srfToXdPtr(i))*qConv
      END DO

c Writing nNeumannSrfs pressures here
      DO i=1, nNeumannSrfs
         IF(IEEE_IS_NAN(X(srfToXPtr(i)))) THEN
            PRINT*, 'Error! NAN encountered..'
            STOP
         END IF
c        Print out pressure from genBC to svFSI for debugging
         PRINT*, 'BC i = ', i, 
     2      'Pressure from genBC: ', offset(srfToXPtr(i))

c        Here, saving offset value to X itself. Necessary for keeping track of 
c        pressure from one time step to the next
         X(srfToXPtr(i))  = offset(srfToXPtr(i))
         
c         WRITE (1) (X(srfToXPtr(i)) - pMin
c     2      + offset(srfToXPtr(i)))*pConv
         WRITE (1) (X(srfToXPtr(i)))*pConv
      END DO
      CLOSE(1)

      IF (flag .EQ. 'L') THEN
         OPEN(1, FILE='InitialData', STATUS='OLD', FORM='UNFORMATTED');
         WRITE (1) t
         DO i=1,nUnknowns
            WRITE (1) X(i)
         END DO
         CLOSE(1)

c         PRINT *,'Before AllData'

         OPEN(1, FILE='AllData', STATUS='UNKNOWN', ACCESS='APPEND');
         string = ''
         DO i=1, nUnknowns
            WRITE (sTmp,"(ES14.6E2)") X(i)
            string = TRIM(string)//sTmp
         END DO
         DO i=1, nXprint
            WRITE (sTmp,"(ES14.6E2)") Xprint(i)
            string = TRIM(string)//sTmp
         END DO
         WRITE (1,"(A)") TRIM(string)
         CLOSE(1)
      END IF

      DEALLOCATE (Pi)
      DEALLOCATE (Pf)
      DEALLOCATE (PDirichlet)
      DEALLOCATE (Qi)
      DEALLOCATE (Qf)
      DEALLOCATE (QNeumann)
      DEALLOCATE (X)
      DEALLOCATE (Xo)
      DEALLOCATE (f)
      DEALLOCATE (srfToXdPtr)
      DEALLOCATE (srfToXPtr)
      DEALLOCATE (offset)
      DEALLOCATE (Xprint)

      END PROGRAM GenBC
