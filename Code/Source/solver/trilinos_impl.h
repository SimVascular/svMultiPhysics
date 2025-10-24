// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef TRILINOS_LINEAR_SOLVER_H
#define TRILINOS_LINEAR_SOLVER_H
/*!
  \file    trilinos_linear_solver.h
  \brief   wrap Trilinos solver functions
*/

/**************************************************************/
/*                          Includes                          */
/**************************************************************/

#include <stdio.h>
#include <vector>
#include <iostream>
#include "mpi.h"
#include <time.h>
#include <numeric>

// Theuchos includes
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultComm.hpp"
#include <Teuchos_Time.hpp> 

#include "Kokkos_Core.hpp" 
#include "Tpetra_KokkosCompat_ClassicNodeAPI_Wrapper.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp" 
#include "MueLu_TpetraOperator.hpp"  
#include "Xpetra_Matrix.hpp"

// Kokkos includes
// #include "KokkosCompat_KokkosSerialWrapperNode.hpp"  // For Node type

// Tpetra includes
#include "Tpetra_Core.hpp"                           // For Tpetra::initialize, finalize
#include "Tpetra_Map.hpp"                            // For Tpetra::Map
#include "Tpetra_CrsMatrix.hpp"                      // If you use Tpetra::CrsMatrix
#include "Tpetra_MultiVector.hpp" 
#include "Tpetra_Map_decl.hpp"
#include "NOX_TpetraTypedefs.hpp"

//Belos includes
#include "BelosSolverFactory_Tpetra.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBiCGStabSolMgr.hpp"
#include "BelosPseudoBlockCGSolMgr.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosStatusTestResNorm.hpp"
#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include <BelosStatusTestImpResNorm.hpp>
#include "BelosStatusTestCombo.hpp"

//Ifpack2 includes
#include "Ifpack2_Factory.hpp"   

// MueLu includes
#include "MueLu_CreateTpetraPreconditioner.hpp"

/**************************************************************/
/*                      Types Definitions                     */
/**************************************************************/
/* Scalar types aliases */
using Scalar_d = double;
using Scalar_i = int;
using Scalar_c = std::complex<double>;

/* Ordinals and node aliases */
using LO = int;
using GO = int;
using Node = Tpetra::Map<>::node_type;

/* Tpetra type aliases */
using Tpetra_Map            = Tpetra::Map<LO, GO, Node>;
using Tpetra_CrsMatrix      = Tpetra::CrsMatrix<Scalar_d, LO, GO, Node>;
using Tpetra_BlockCrsMatrix = Tpetra::BlockCrsMatrix<Scalar_d, LO, GO, Node>;
using Tpetra_MultiVector    = Tpetra::MultiVector<Scalar_d, LO, GO, Node>; 
using Tpetra_Vector         = Tpetra::Vector<Scalar_d, LO, GO, Node>;
using Tpetra_Import         = Tpetra::Import<LO, GO, Node>;
using Tpetra_CrsGraph       = Tpetra::CrsGraph<LO, GO, Node>;
using Tpetra_Operator       = Tpetra::Operator<Scalar_d, LO, GO, Node>;

/* Belos aliases */
using Belos_LinearProblem = Belos::LinearProblem<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;
using Belos_SolverFactory =  Belos::TpetraSolverFactory<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;
using Belos_SolverManager = Belos::SolverManager<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;
using Belos_StatusTestGenResNorm = Belos::StatusTestGenResNorm<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;
using Belos_StatusTestCombo = Belos::StatusTestCombo<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;
using Belos_StatusTestMaxIters = Belos::StatusTestMaxIters<Scalar_d, Tpetra_MultiVector, Tpetra_Operator>;

/* IFPACK2 preconditioner aliases */  
using Ifpack2_Preconditioner = Ifpack2::Preconditioner<Scalar_d, LO, GO, Node>;

/* MueLu preconditioner aliases */
using MueLu_Preconditioner = Tpetra_Operator;

/**************************************************************/
/*                      Macro Definitions                     */
/**************************************************************/

// Define linear solver as following naming in FSILS_struct
#define TRILINOS_CG_SOLVER 798
#define TRILINOS_GMRES_SOLVER 797
#define TRILINOS_BICGSTAB_SOLVER 795

// Define preconditioners as following naming in FSILS_struct
#define NO_PRECONDITIONER 700
#define TRILINOS_DIAGONAL_PRECONDITIONER 702
#define TRILINOS_BLOCK_JACOBI_PRECONDITIONER 703
#define TRILINOS_ILU_PRECONDITIONER 704
#define TRILINOS_ILUT_PRECONDITIONER 705
#define TRILINOS_RILUK0_PRECONDITIONER 706
#define TRILINOS_RILUK1_PRECONDITIONER 707
#define TRILINOS_ML_PRECONDITIONER 708

/// @brief Initialize all Epetra types we need separate from Fortran
struct Trilinos
{
  Teuchos::RCP<const Tpetra_Map> Map;
  Teuchos::RCP<const Tpetra_Map> ghostMap;
  Teuchos::RCP<Tpetra_MultiVector> F;
  Teuchos::RCP<Tpetra_MultiVector> ghostF;
  Teuchos::RCP<Tpetra_CrsMatrix> K;
  Teuchos::RCP<Tpetra_Vector> X;
  Teuchos::RCP<Tpetra_Vector> ghostX;
  Teuchos::RCP<Tpetra_Import> Importer;
  std::vector<Teuchos::RCP<Tpetra_MultiVector>> bdryVec_list;
  Teuchos::RCP<const Teuchos::Comm<int>> comm;
  Teuchos::RCP<Tpetra_CrsGraph> K_graph;

  Teuchos::RCP<Tpetra_Operator> MueluPrec;
  Teuchos::RCP<Ifpack2_Preconditioner> ifpackPrec;
  Trilinos() : MueluPrec(nullptr), ifpackPrec(nullptr) {}
};

/**
 * \class TrilinosMatVec
 * \brief This class implements the pure virtual class Epetra_Operator for the
 *        AztecOO iterative solve which only uses the Apply() method to compute
 *        the matrix vector product
 */
class TrilinosMatVec: public Tpetra_Operator
{
public:

  /** Define matrix vector operation at each iteration of the linear solver
   *  adds on the coupled neuman boundary contribution to the matrix
   *
   *  \param x vector to be applied on the operator
   *  \param y result of sparse matrix vector multiplication
   */
  TrilinosMatVec(const Teuchos::RCP<Trilinos>& trilinos) : trilinos_(trilinos) {}

  /* Y = beta * Y + alpha * A^mode * X */
  void apply(const Tpetra_MultiVector& X, Tpetra_MultiVector& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar_d alpha = Teuchos::ScalarTraits< Scalar_d >::one(), 
           Scalar_d beta  = Teuchos::ScalarTraits<Scalar_d>::zero()) const override;

  /*  
    Returns the map describing the layout of the domain vector space.
    This map defines the distribution of the input vectors to the operator.
  */
  Teuchos::RCP<const Tpetra_Map> getDomainMap() const override
  {
    return trilinos_->K->getDomainMap();
  }

  /* 
    Returns the map describing the layout of the range vector space.
    This map defines the distribution of the output vectors from the operator.
  */
  Teuchos::RCP<const Tpetra_Map> getRangeMap() const override
  {
    return trilinos_->K->getRangeMap();
  }

  private:
    Teuchos::RCP<Trilinos> trilinos_;
};// class TrilinosMatVec

//  --- Functions to be called in fortran -------------------------------------

#ifdef __cplusplus
  extern "C"
  {
#endif
  /// Give function definitions which will be called through fortran
  void trilinos_lhs_create(const Teuchos::RCP<Trilinos> &trilinos_, const int numGlobalNodes, const int numLocalNodes,
          const int numGhostAndLocalNodes, const int nnz, const Vector<int>& ltgSorted,
          const Vector<int>& ltgUnsorted, const Vector<int>& rowPtr, const Vector<int>& colInd,
          const int dof, const int cpp_index, const int proc_id, const int numCoupledNeumannBC);

  /**
   * \param v           coeff in the scalar product
   * \param isCoupledBC determines if coupled resistance BC is turned on
   */
  void trilinos_bc_create_(const Teuchos::RCP<Trilinos> &trilinos_, const std::vector<Array<double>> &v_list, bool &isCoupledBC);

  void trilinos_doassem_(const Teuchos::RCP<Trilinos> &trilinos_, int &numNodesPerElement, const int *eqN,
          const double *lK, double *lR);

  void trilinos_global_solve_(const Teuchos::RCP<Trilinos> &trilinos_, const double *Val, const double *RHS,
          double *x, const double *dirW, double &resNorm, double &initNorm,
          int &numIters, double &solverTime, double &dB, bool &converged,
          int &lsType, double &relTol, int &maxIters, int &kspace,
          int &precondType);

  void trilinos_solve_(const Teuchos::RCP<Trilinos> &trilinos_, double *x, const double *dirW, double &resNorm,
          double &initNorm, int &numIters, double &solverTime,
          double &dB, bool &converged, int &lsType, double &relTol,
          int &maxIters, int &kspace, int &precondType, bool &isFassem);

#ifdef __cplusplus  /* this brace matches the one on the extern "C" line */
  }
#endif

// --- Define functions to only be called in C++ ------------------------------
void setPreconditioner(const Teuchos::RCP<Trilinos> &trilinos_, int precondType, 
  Teuchos::RCP<Belos_LinearProblem>& BelosProblem);

void setMueLuPreconditioner(Teuchos::RCP<MueLu_Preconditioner>& MueLuPrec, 
  const Teuchos::RCP<Tpetra_CrsMatrix>& A);

void checkDiagonalIsZero(const Teuchos::RCP<Trilinos> &trilinos_);

void constructJacobiScaling(const Teuchos::RCP<Trilinos> &trilinos_, const double *dirW,
              Tpetra_Vector &diagonal);

// --- Debugging functions ----------------------------------------------------
void printMatrixToFile(const Teuchos::RCP<Trilinos> &trilinos_);

void printRHSToFile(const Teuchos::RCP<Trilinos> &trilinos_);

void printSolutionToFile(const Teuchos::RCP<Trilinos> &trilinos_);

#endif //TRILINOS_LINEAR_SOLVER_H
