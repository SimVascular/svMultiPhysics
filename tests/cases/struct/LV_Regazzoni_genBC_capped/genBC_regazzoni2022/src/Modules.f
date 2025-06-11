c     Created by Mahdi Esmaily Moghadam 05-25-2011
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

      MODULE COM

      LOGICAL pCorr, qCorr

      INTEGER nUnknowns, nDirichletSrfs, nNeumannSrfs, nXprint

      REAL(KIND=8) pConv, qConv

      INTEGER, ALLOCATABLE :: srfToXPtr(:), srfToXdPtr(:)

      REAL(KIND=8), ALLOCATABLE :: QNeumann(:,:), PDirichlet(:,:),
     2   offset(:), Xprint(:)

      END MODULE COM
