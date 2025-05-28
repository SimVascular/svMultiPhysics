/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//--------------------------------------------------------------------
//
// Interface to Trilinos linear solver library.
//
//--------------------------------------------------------------------

/*!
  \file    trilinos_linear_solver.cpp
  \brief   wrap Trilinos solver functions
*/

#include "trilinos_impl.h"
#include "ComMod.h"
#define NOOUTPUT

// --- Define global Trilinos variables to be used in below functions ---------

/// Unique block map consisting of nodes owned by each processor
Epetra_BlockMap *Trilinos::blockMap;

/// Global block force vector
Epetra_FEVector *Trilinos::F;

/// Global block stiffness matrix
Epetra_FEVbrMatrix *Trilinos::K;

/// Solution vector consisting of unique nodes owned by the processor
Epetra_Vector *Trilinos::X;

/// Solution vector with nodes owned by processor followed by its ghost nodes
Epetra_Vector *Trilinos::ghostX;

/// Import ghostMap into blockMap to create ghost map
Epetra_Import *Trilinos::Importer;

/// Contribution from coupled neumann boundary conditions. One vector for each coupled Neumann boundary
std::vector<Epetra_FEVector*> Trilinos::bdryVec_list;

Epetra_FECrsGraph *Trilinos::K_graph;

ML_Epetra::MultiLevelPreconditioner* MLPrec = NULL;

Ifpack_Preconditioner* ifpackPrec = NULL;

// --- Define global variables to be shared amongst the below functions -------

/// Nodal degrees of freedom
int dof;

/// Total number of nodes including the ghost nodes
int ghostAndLocalNodes;

/// Nodes owned by processor
int localNodes;

/// Converts local proc column indices to global indices to be inserted
std::vector<int> globalColInd;

/// Converts local indices to global indices in unsorted ghost node order
std::vector<int> localToGlobalUnsorted;

/// Stores number of nonzeros per row for the topology
std::vector<int> nnzPerRow;

std::vector<int> localToGlobalSorted;

int timecount = 0;

bool coupledBC;

// ----------------------------------------------------------------------------
/**
 * Define the matrix vector multiplication operation to do at each iteration of
 * an iterative linear solver. This function is called by the AztecOO
 * (K + v1*v1' + v2*v2' + ...) *x = K*x + v1*(v1'*x) + v2*(v2'*x), 
 * where [v1, v2, ...] = bdryVec_list
 *
 * For coupled Neumann boundary terms (v1*v1', v2*v2'), we use efficient an vectorized 
 * operation rather than explicitly forming the rank 1 outer product matrix and 
 * adding it to the global stiffness matrix K.
 *
 * \param x vector to be applied on the operator
 * \param y result of sprase matrix vector multiplication
 */
int TrilinosMatVec::Apply(const Epetra_MultiVector &x,
        Epetra_MultiVector &y) const
{
  // Store initial matrix vector product result in y (y = K*x)
  Trilinos::K->Apply(x, y); //K*x

  // Now add on the coupled Neumann boundary contribution y += v1*(v1'*x) + v2*(v2'*x) + ...
  if (coupledBC)
  {
    // Declare dot product v_i'*x
    double dot = 0.0;

    // Loop over all coupled Neumann boundary vectors
    for (auto bdryVec : Trilinos::bdryVec_list)
    {
      // Compute dot product dot = v_i'*x
      bdryVec->Dot(x, &dot);

      // y = 1*y + dot*v_i
      y.Update(dot,
              *bdryVec,
              1.0);
    }
  }
  return 0;
}

// ----------------------------------------------------------------------------
/**
 * make the graph global and optimize storage on it
 * this function creates the LHS structure topology and RHS vector based on the
 * block map
 * \param numGlobalNodes total/global number of nodes each node can have dof
 *                       for the spatial dim
 * \param numLocalNodes  number of nodes owned by this proc in block coordiantes
 * \param numGhostAndLocalNodes number of nodes owned and shared by this
 *                              processor includes ghost nodes
 * \param nnz            number of nonzeros in LHS matrix of calling proc
 * \param ltgSorted      integer pointer of size numLocalNodes returns global
 *                       node number of local node a to give blockrow
 * \param ltgUnsorted    unsorted/not-reordered version
 * \param rowPtr         CSR row ptr of size numLocalNodes + 1 to block rows
 * \param colInd         CSR column indices ptr (size nnz points) to block rows
 * \param Dof            size of each block element to give dim of each block
 * 
 * \param numCoupledNeumannBC number of coupled Neumann BCs

 */

void trilinos_lhs_create(const int numGlobalNodes, const int numLocalNodes,
        const int numGhostAndLocalNodes, const int nnz, const Vector<int>& ltgSorted,
        const Vector<int>& ltgUnsorted, const Vector<int>& rowPtr, const Vector<int>& colInd,
        const int Dof, const int cpp_index, const int proc_id, const int numCoupledNeumannBC)
{
  #ifdef debug_trilinos_lhs_create
  std::string msg_prefix;
  msg_prefix = std::string("[trilinos_lhs_create:") + std::to_string(proc_id) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== trilinos_lhs_create ==========" << std::endl;
  std::cout << msg_prefix << "Dof: " << Dof << std::endl;
  std::cout << msg_prefix << "cpp_index: " << cpp_index << std::endl;
  #endif

  int indexBase = 1; //0 based indexing for C/C++ / 1 for Fortran
  if (cpp_index == 1) {
    indexBase = 0; 
  }

  dof = Dof; //constant size dof blocks
  ghostAndLocalNodes = numGhostAndLocalNodes;
  localNodes = numLocalNodes;
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  #ifdef debug_trilinos_lhs_create
  std::cout <<  msg_prefix << "indexBase: " << indexBase << std::endl;
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "ghostAndLocalNodes: " << ghostAndLocalNodes << std::endl;
  std::cout << msg_prefix << "localNodes: " << localNodes << std::endl;
  std::cout << msg_prefix << "localToGlobalSorted.size(): " << localToGlobalSorted.size() << std::endl;
  #endif

  bool new_mapping_pattern = false; //could change with remeshing

  // check if mesh was changed and vectors need to be resized
  if (localToGlobalSorted.size() != numGhostAndLocalNodes) {
    localToGlobalSorted.clear();   // sets size to zero
    localToGlobalUnsorted.clear();
    nnzPerRow.clear();
    globalColInd.clear();
  }

  // allocate memory for vectors
  if (localToGlobalSorted.size() == 0) { // vectors not allocated or cleared
    localToGlobalSorted.reserve(numGhostAndLocalNodes);
    localToGlobalUnsorted.reserve(numGhostAndLocalNodes);
    nnzPerRow.reserve(numGhostAndLocalNodes);
    globalColInd.reserve(nnz);
    new_mapping_pattern = true;
  }

  // Define localtoglobal to be used for unqiue partition map
  //only take ltgSorted(1:numLocalNodes) since those are owned by the processor
  //
  if (new_mapping_pattern) {
    for (unsigned i = 0; i < numLocalNodes; ++i) {
      //any nodes following are ghost nodes so that subset
      localToGlobalSorted.emplace_back(ltgSorted[i]);
    }
  }

  // Create blockmap for arbitrary distributon of constant size, dof, element
  // amongst the procs. This is a non-overlapping map defining a unique partition
  // localToGlobalSorted[0] contains global index value of ith owned element on
  // the calling processor
  //
  Trilinos::blockMap = new Epetra_BlockMap(numGlobalNodes, numLocalNodes,
                       &localToGlobalSorted[0], dof, indexBase, comm);

  // Create blockmap which includes the ghost nodes to be imported into the
  // solution vector at the end
  //
  int inc_ghost = -1; // since including ghost nodes sum will not equal total
  Epetra_BlockMap ghostMap(inc_ghost, numGhostAndLocalNodes, ltgSorted.data(), dof,
          indexBase, comm);

  // Calculate nnzPerRow to pass into graph constructor
  //
  if (new_mapping_pattern) {
    nnzPerRow.clear (); //destroy vector
    for (unsigned i = 0; i < numLocalNodes; ++i) {
      nnzPerRow.emplace_back(rowPtr[i+1] - rowPtr[i]);
    }
  } 

  // Construct graph based on nnz per row
  //
  Trilinos::K_graph = new Epetra_FECrsGraph(Copy, *Trilinos::blockMap, &nnzPerRow[0]);

  unsigned nnzCount = 0; //cumulate count of block nnz per rows

  // Use unsortedltg to map local col ind to global indices
  //
  if (globalColInd.size() != nnz) { // only if nnz changed
    globalColInd.clear(); // destroy
    for (unsigned i = 0; i < nnz; ++i) {
      // Convert to global indexing subtract 1 for fortran vs C array indexing
      globalColInd.emplace_back(ltgUnsorted[colInd[i] - indexBase]);
    }
  } 

  //loop over block rows owned by current proc using localToGlobal index pointer
  for (unsigned i = 0; i < numGhostAndLocalNodes; ++i)
  {
    if (new_mapping_pattern && i >= numLocalNodes)
      nnzPerRow.emplace_back(rowPtr[i+1] - rowPtr[i]);
    int numEntries = nnzPerRow[i]; //nnz per row
    int error = 0;
    int num_rows_inserting = 1; //number of rows-inserting one row at a time
    error = Trilinos::K_graph->InsertGlobalIndices(num_rows_inserting,
            &ltgUnsorted[i], numEntries, &globalColInd[nnzCount]);

    //colInd is indexed by nnz per rows processed so far
    if (error != 0)
    {
      std::cout << "ERROR: Inserting Global Indices into Map" << std::endl;
      exit(1);
    }

    // final result should be nnz and so will not go out of bounds on colInd
    // up to nnz -1
    nnzCount += numEntries;
    //store into global vector
    if (new_mapping_pattern)
      localToGlobalUnsorted.emplace_back(ltgUnsorted[i]);
  } // for

  //by end of iterations nnzCount should equal nnz-otherwise there is an error
  if (nnzCount != nnz)
  {
    std::cout << "Number of Entries Entered in Graph does not equal nnz " << std::endl;
    exit(1);
  }

  //check if trilinos methods are successful
  if (Trilinos::K_graph->GlobalAssemble() != 0) //Calls FillComplete
  {
    std::cout << "ERROR: Calling Fill Complete on Graph" << std::endl;
    exit(1);
  }

  //Trilinos::K_graph->OptimizeStorage(); //TODO: that causes it to fail

  // --- Create block finite element matrix from graph with fillcomplete ------
  // construct matrix from filled graph
  Trilinos::K = new Epetra_FEVbrMatrix(Copy, *Trilinos::K_graph);
  //construct RHS force vector F topology
  Trilinos::F = new Epetra_FEVector(*Trilinos::blockMap);
  // Construct a boundary vector for each coupled Neumann boundary condition
  Trilinos::bdryVec_list.clear();
  for (int i = 0; i < numCoupledNeumannBC; ++i)
  {
    Trilinos::bdryVec_list.push_back(new Epetra_FEVector(*Trilinos::blockMap));
  }

  // Initialize solution vector which is unique and does not include the ghost
  // indices using the unique map
  Trilinos::X = new Epetra_Vector(*Trilinos::blockMap);
  //initialize vector which will import the ghost nodes using the ghost map
  Trilinos::ghostX = new Epetra_Vector(ghostMap);
  //Create importer of the two maps
  Trilinos::Importer = new Epetra_Import(ghostMap, *Trilinos::blockMap);
  //need to give v a map-same as F
  if (new_mapping_pattern)
    for (unsigned i = numLocalNodes; i < numGhostAndLocalNodes; ++i)
      localToGlobalSorted.emplace_back(ltgSorted[i]);
} // trilinos_lhs_create_

// ----------------------------------------------------------------------------
/**
 * This function assembles the dense element block stiffness matrix and element
 * block force vector into the global block stiffness K and global block force
 * vector F.
 * Happens on a single element and will be called repeatedly as we loop over
 * all the elements.
 *
 * \param numNodesPerElement number of nodes (d) for element e
 * \param eqN                converts local element node number to local proc
 *                           equation number-size is numNodesPerElement
 * \param lK                 element local stiffness of size
 *                           (dof*dof, numNodesPerElement, numNodesPerElement)
 * \param lR                 element local force vector of size
 *                           (dof, numNodesPerElement)
 */
void trilinos_doassem_(int &numNodesPerElement, const int *eqN, const double *lK, double *lR)
{
  #ifdef debug_trilinos_doassem
  std::cout << "[trilinos_doassem_] ========== trilinos_doassem_ ===========" << std::endl;
  std::cout << "[trilinos_doassem_] dof: " << dof << std::endl;
  std::cout << "[trilinos_doassem_] numNodesPerElement: " << numNodesPerElement << std::endl;
  #endif

  //dof values per global ID in the force vector
  int numValuesPerID = dof;

  //converts eqN in local proc values to global values
  std::vector<int> localToGlobal(numNodesPerElement);

  for (int i = 0; i < numNodesPerElement; ++i)
    localToGlobal[i] = localToGlobalUnsorted[eqN[i]];

  //loop over local nodes on the element
  for (int a = 0; a < numNodesPerElement; ++a)
  {
    // Sum into contributions from element node-global assemble will take those
    // from shared nodes on other processors since FE routine.
    int error = Trilinos::K->BeginSumIntoGlobalValues(localToGlobal[a],
            numNodesPerElement, &localToGlobal[0]);

    if (error != 0)
    {
      std::cout << "ERROR: Setting block row and block column of summing into "
                << "global values!"
                << std::endl;
      exit(1);
    }

    //submit global F values
    int num_block_rows = 1; //number of global block rows put in 1 at a time
    Trilinos::F->SumIntoGlobalValues (num_block_rows, &localToGlobal[a],
            &numValuesPerID, &lR[a*dof]);

    if (error != 0) //1 or more indices not associated with calling processor!
    {
      std::cout << error << std::endl;
      std::cout << "ERROR: Summing into global vector values!" << std::endl;
      exit(1);
    }

    //loop over local nodes for columns
     std::vector<double> values(dof*dof);
    for (int b = 0; b < numNodesPerElement; ++b)
    {
      //transpose block since Trilinos takes in SerialMAtrix in column major
      for (int i = 0; i < dof; ++i)
      {
        //premult by diagonal W so W(a) multiplies row a
        for (int j = 0; j < dof; ++j)
        {
          //taking transpose of block so flip i & j
          values[i*dof + j]
              = lK[b*dof*dof*numNodesPerElement + a*dof*dof + j*dof + i];
        }
      }

      error = Trilinos::K->SubmitBlockEntry(&values[0], dof, dof, dof);

      if (error != 0)
      {
        std::cout << "ERROR: Inputting values of summing into matrix "
                  << "global values!"
                  << std::endl;
        exit(1);
      }
    }

    //finish submitting the block entry for the current row
    error = Trilinos::K->EndSubmitEntries();
    if (error != 0)
    {
      std::cout << "[trilinos_doassem_] ERROR: End submitting block entries!" << std::endl;
      exit(1);
    }
  }
} // trilinos_doassem_

// ----------------------------------------------------------------------------
/**
 * Tthis function creates the LHS structure topology and RHS vector based on
 * the block map using the global stiffness matrix and global force vector
 * assuming the assembly part has been done in svFSI
 *
 * \param Val         nonzero values of assembled global stiffness
 * \param RHS         global force vector
 * \param x           solution vector
 * \param dirW        gives information about the Dirichlet BC
 * \param resNorm     norm of solution
 * \param initNorm    initial norm of RHS
 * \param numIters    linear solver number of iterations
 * \param solverTime  time in linear solver
 * \param dB          log ratio for output
 * \param converged   neither or not converged in the max number of iterations
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 */

void trilinos_global_solve_(const double *Val, const double *RHS, double *x,
        const double *dirW, double &resNorm, double &initNorm, int &numIters,
        double &solverTime, double &dB, bool &converged, int &lsType,
        double &relTol, int &maxIters, int &kspace, int &precondType)
{
  int nnzCount = 0; //cumulate count of block nnz per rows
  int count = 0;
  int numValuesPerID = dof; //dof values per id pointer to dof
  std::vector<double> values(dof*dof); // holds local matrix entries

  // loop over block rows owned by current proc using localToGlobal index pointer
  //
  for (int i = 0; i < ghostAndLocalNodes; ++i) {
    int numEntries = nnzPerRow[i]; //block per of entries per row
    // Copy global stiffness values
    int error = Trilinos::K->BeginReplaceGlobalValues(localToGlobalUnsorted[i],
            numEntries, &globalColInd[nnzCount]);

    // Need to check if globalColInd and localToGlobal equal for block diag
    //
    if (error != 0) {
      std::cout << error << std::endl;
      std::cout << "ERROR: Setting block row and block column of setting "
                << "global values!"
                << std::endl;
      exit(1);
    }

    // Copy F values-loop over dof for bool
    //
    int num_block_rows = 1; //number of global block rows put in 1 at a time

    error = Trilinos::F->ReplaceGlobalValues (num_block_rows,
            &localToGlobalUnsorted[i], &numValuesPerID, &RHS[i*dof]);

    //check is bool true or false whether to give 0 if on diagonal 1
    if (error != 0) {
      std::cout << "ERROR: Setting global vector values!" << std::endl;
      exit(1);
    }

    for (int j = 0; j < numEntries; ++j) {
      for (int l = 0; l < dof; ++l) { //loop over dof for bool to contruct
        for (int m = 0; m < dof; ++m) {
          values[l*dof + m] = Val[count*dof*dof + m*dof + l]; //transpose it
        }
      }

      // Submit square dof*dof blocks
      error = Trilinos::K->SubmitBlockEntry(&values[0], dof, dof, dof);

      if (error != 0) {
        std::cout << "ERROR: Inputting values of setting matrix global "
                  << "values!" << std::endl;
        exit(1);
      }
      count++;
    }

    error = Trilinos::K->EndSubmitEntries(); //for current block row

    if (error != 0) {
      std::cout << "ERROR: End submitting block entries!" << std::endl;
      exit(1);
    }
    nnzCount += numEntries;
  }

  // Call solver code which assembles K and F for shared processors
  //
  bool flagFassem = false;

  trilinos_solve_(x, dirW, resNorm, initNorm, numIters,
          solverTime, dB, converged, lsType,
          relTol, maxIters, kspace, precondType, flagFassem);

} // trilinos_global_solve_

// ----------------------------------------------------------------------------
/**
 * This function uses the established data structure to solve the block linear
 * system and passes the solution vector back to fortran with the ghost nodes
 * filled in
 *
 * \param x           solution vector
 * \param dirW        gives information about Dirichlet BC
 * \param resNorm     norm of the final residual
 * \param initNorm    norm of initial resiual x_init = 0 ||b||
 * \param numIters    number of iterations in the linear solver
 * \param solverTime  time in the linear solver
 * \param dB          log ratio for output
 * \param converged   can pass in solver and preconditioner type too
 * \param lsType      defines type of linear solver to use
 * \param relTol      default relative tolerance or as set in input file
 * \param maxIters    default max number of iterations for gmres per restart
 * \param kspace      specific for gmres dim of the stored Krylov space vectors
 * \param precondType defines type of preconditioner to use
 * \param isFassem    determines if F is already assembled at ghost nodes
 */
void trilinos_solve_(double *x, const double *dirW, double &resNorm,
        double &initNorm, int &numIters, double &solverTime, double &dB,
        bool &converged, int &lsType, double &relTol, int &maxIters,
        int &kspace, int &precondType, bool &isFassem)
{
  #define n_debug_trilinos_solve
  #ifdef debug_trilinos_solve
  std::cout << "[trilinos_solve] ========== trilinos_solve ==========" << std::endl;
  std::cout << "[trilinos_solve] resNorm: " << resNorm << std::endl;
  std::cout << "[trilinos_solve] initNorm: " << initNorm << std::endl;
  std::cout << "[trilinos_solve] lsType: " << lsType << std::endl;
  std::cout << "[trilinos_solve] precondType: " << precondType << std::endl;
  std::cout << "[trilinos_solve] isFassem: " << isFassem << std::endl;
  #endif
  bool flagFassem = isFassem;

  // Already filled from graph so does not need to call fillcomplete
  // routine will sum in contributions from elements on shared nodes amongst
  // processors
  //
  int error = Trilinos::K->GlobalAssemble(false);
  if (error != 0) {
    std::cout << "ERROR: Global Assembling stiffness matrix" << std::endl;
    exit(1);
  }

  //very important for performance-makes memory contiguous
  Trilinos::K->OptimizeStorage();

  if (flagFassem) {
    //sum in values from shared nodes amongst the processors
    error = Trilinos::F->GlobalAssemble();

    if (error != 0) {
      std::cout << "ERROR: Global Assembling force vector" << std::endl;
      exit(1);
    }
  }

  // Construct Jacobi scaling vector which uses dirW to take the Dirichlet BC
  // into account
  //
  Epetra_Vector diagonal(*Trilinos::blockMap);
  constructJacobiScaling(dirW, diagonal);

  // Compute norm of preconditioned multivector F that we will be solving
  // problem with
  Trilinos::F->Norm2(&initNorm); //pass preconditioned norm W*F

  // Define Epetra_Operator which is global stiffness with coupled Neumann BC
  // contributions included
  TrilinosMatVec K_bdry;

  // Define linear problem if v is 0 does standard matvec product with K
  Epetra_LinearProblem Problem(&K_bdry, Trilinos::X, Trilinos::F);

  AztecOO Solver(Problem);

  // Can set output solver parameter options below
#ifdef NOOUTPUT
  Solver.SetAztecOption(AZ_diagnostics, AZ_none);
  Solver.SetAztecOption(AZ_output, AZ_none);
#endif

  setPreconditioner(precondType, Solver);

  // Set convergence type as relative ||r|| <= relTol||b||
  Solver.SetAztecOption(AZ_conv, AZ_rhs);

  //Set solver options
  int numRestarts = 1; //also changes for gmres
  int maxItersPerRestart = maxIters;

  // Solver is GMRES by default
  if (lsType ==TRILINOS_GMRES_SOLVER)
  {
    //special parameters to set orthog and kspace
    numRestarts = maxIters; //different definition for gmres
    Solver.SetAztecOption(AZ_orthog, AZ_modified); //modified GS
    Solver.SetAztecOption(AZ_kspace, kspace);
    Solver.SetAztecOption(AZ_solver, AZ_gmres);
    maxItersPerRestart = kspace; //total maxIters is kspace * numRestarts
  }
  else if (lsType == TRILINOS_BICGSTAB_SOLVER)
    Solver.SetAztecOption(AZ_solver, AZ_bicgstab);

  else if (lsType == TRILINOS_CG_SOLVER)
    Solver.SetAztecOption(AZ_solver, AZ_cg);

  //checkStatus to calculate residual norm
  AztecOO_StatusTestResNorm restartResNorm(K_bdry, *Trilinos::X,
                      (Epetra_Vector&) Trilinos::F[0], relTol);
  restartResNorm.DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                 AztecOO_StatusTestResNorm::TwoNorm);
  Solver.SetStatusTest(&restartResNorm);

  int status;
  numIters = 0;
  solverTime = 0.0;
  for (int iter = 0; iter < numRestarts; ++iter)
  {
    dB = (iter > 0) ? resNorm : initNorm;
    if (numRestarts > 1) numIters += 1; //gmres case
    status = Solver.Iterate(maxItersPerRestart, relTol);
    numIters += Solver.NumIters();
    resNorm = restartResNorm.GetResNormValue();
    solverTime += Solver.SolveTime();
    if (resNorm < relTol * initNorm) break;
  }
  converged = (status == 0) ? true : false;
  dB = 10 * log(restartResNorm.GetResNormValue()/dB); //fits with gmres def

  //Right scaling so need to multiply x by diagonal
  Trilinos::X->Multiply(1.0, *Trilinos::X, diagonal, 0.0);

  //Fill ghost X with x communicating ghost nodes amongst processors
  error = Trilinos::ghostX->Import(*Trilinos::X, *Trilinos::Importer, Insert);
  //check imported correctly
  if (error != 0)
  {
    std::cout << "ERROR: Map ghost node importer!" << std::endl;
     exit(1);
  }

    error = Trilinos::ghostX->ExtractCopy(x);
   if (error != 0)
   {
     std::cout << "ERROR: Extracting copy of solution vector!" << std::endl;
     exit(1);
  }
  //set to 0 for the next time iteration
  Trilinos::K->PutScalar(0.0);
  Trilinos::F->PutScalar(0.0);
  if (coupledBC) {
    for (auto bdryVec : Trilinos::bdryVec_list)
      bdryVec->PutScalar(0.0);

  }
  //0 out initial guess for iteration
  Trilinos::X->PutScalar(0.0);
  // Free memory if MLPrec is invoked
  if (ifpackPrec) {
      delete ifpackPrec;
      ifpackPrec = NULL;
  }
  if (MLPrec) {
      delete MLPrec;
      MLPrec = NULL;
  }
} // trilinos_solve_

// ----------------------------------------------------------------------------
void setPreconditioner(int precondType, AztecOO &Solver)
{
  //initialize reordering for ILU/ILUT preconditioners
  Solver.SetAztecOption(AZ_reorder, 1);
  //solve precond structure into separate function
  if (precondType == TRILINOS_DIAGONAL_PRECONDITIONER ||
      precondType == NO_PRECONDITIONER )
    Solver.SetAztecOption(AZ_precond, AZ_none);
  else if (precondType == TRILINOS_BLOCK_JACOBI_PRECONDITIONER)
  {
    Solver.SetAztecOption(AZ_precond, AZ_Jacobi);
    checkDiagonalIsZero();
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if(precondType == TRILINOS_ILU_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilu);
    Solver.SetAztecOption(AZ_overlap,1);
    Solver.SetAztecOption(AZ_graph_fill,0);
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if (precondType == TRILINOS_ILUT_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    Solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
    Solver.SetAztecOption(AZ_overlap,1);
    Solver.SetAztecOption(AZ_graph_fill,4);
    Solver.SetAztecParam(AZ_ilut_fill, 2);
    Solver.SetAztecParam(AZ_drop, 1e-2);

    // Sets the global stiffness prior to rank 1 update as the matrix off of
    // which the preconditioner to calculated from to utilize the native
    // preconditioners from trilinos
    Solver.SetPrecMatrix(Trilinos::K);
  }
  else if (precondType == TRILINOS_IC_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    setIFPACKPrec(Solver);
  }
  else if (precondType == TRILINOS_ICT_PRECONDITIONER)
  {
    checkDiagonalIsZero();
    setIFPACKPrec(Solver); //add in parameter for string for different types
  }
  else if (precondType == TRILINOS_ML_PRECONDITIONER)
    setMLPrec(Solver);
  else
  {
    std::cout << "ERROR: Preconditioner Type is undefined" << std::endl;
    exit(1);
  }
} // setPreconditioner

// ----------------------------------------------------------------------------
/**
 * Tune parameters for htis and IFPACK
 * Ref: https://trilinos.org/oldsite/packages/ml/mlguide5.pdf
 */
void setMLPrec(AztecOO &Solver)
{
  //break up into initializer
  Teuchos::ParameterList MLList;
  int *options = new int[AZ_OPTIONS_SIZE];
  double *params = new double[AZ_PARAMS_SIZE];
  ML_Epetra::SetDefaults("SA",MLList, options, params);
  //MLList.set("XML input file", "/home/augustin/programs/MUPFES/trunk/examples/conf-files/ml_ParameterList_sGS.xml");

  // ML general options
  // output level, 0 being silent and 10 verbose
  bool verbose = false;
  if (verbose)
    MLList.set("ML output", 10);
  else
    MLList.set("ML output", 0);

  // ML Cycle options
  // maximum number of levels possible
  MLList.set("max levels",4);
  //If set to increasing, level 0 will correspond to the finest level.
  //If set to decreasing,max levels - 1 will correspond to the finest level.
  //Default: increasing.
  MLList.set("increasing or decreasing","increasing");

  // Aggregation & prolongator parameters
  // coarsening options:  Uncoupled, MIS, Uncoupled-MIS (uncoupled on the finer grids, then switch to MIS)
  //MLList.set("aggregation: type", "Uncoupled");
  MLList.set("aggregation: type", "MIS");
  MLList.set("aggregation: threshold", 0.0);

  // Smoother parameters
  //common smoother options: Chebyshev, Gauss-Seidel, symmetric Gauss-Seidel, Jacobi, ILU, IC
  MLList.set("smoother: ifpack type","ILU");
  MLList.set("smoother: ifpack overlap",1);
  MLList.sublist("smoother: ifpack list").set("fact: level-of-fill",1);
  MLList.sublist("smoother: ifpack list").set("schwarz: reordering type","rcm");
  // use both pre and post smoothing
  MLList.set("smoother: pre or post", "both");
  MLList.set("subsmoother: type", "symmetric Gauss-Seidel");

  //Coarse grid parameters
  MLList.set("coarse: max size", 15);
  // Coarse solver Default: AmesosKLU.
  MLList.set("coarse: type","symmetric Gauss-Seidel");

  // Load balancing options
  MLList.set("repartition: enable",1);
  MLList.set("repartition: max min ratio",1.3);
  MLList.set("repartition: min per proc",500);
  MLList.set("repartition: partitioner","Zoltan");
  MLList.set("repartition: Zoltan dimensions",2);

  // create the preconditioner object based on options in MLList and compute hierarchy
  if (MLPrec == NULL)
    MLPrec = new ML_Epetra::MultiLevelPreconditioner(*Trilinos::K, MLList, false);

  timecount = 0;
  if (timecount == 0)
    MLPrec->ComputePreconditioner();
  else {//switch to recompute
    MLPrec->ReComputePreconditioner();
  }
  Solver.SetPrecOperator(MLPrec);

  //solver to separte function only recompute at separate time iter
  timecount += 1;

  delete[] options;
  delete[] params;
}// setMLPrec

// ----------------------------------------------------------------------------
/**
 * pass in IC, ICT
 * pass in string for which to turn on right now set to IC
 */
void setIFPACKPrec(AztecOO &Solver)
{
  //Ifpack Factory;
  //std::string PrecType = "ILUT"; // exact solve on each subdomain
  //int OverlapLevel = 0; // one row of overlap among the processes
  //ifpackPrec = Factory.Create(PrecType, Trilinos::K, OverlapLevel);
  //Teuchos::ParameterList List;
  //List.set("fact: level-of-fill", 2);
  //List.set("fact: drop tolerance", 1e-2);
  //ifpackPrec->SetParameters(List);
  //ifpackPrec->Initialize();
  //ifpackPrec->Compute();
  //Solver.SetPrecOperator(&*ifpackPrec);

  Teuchos::ParameterList List;
  int OverlapLevel = 0;
  ifpackPrec = new Ifpack_AdditiveSchwarz<Ifpack_ILUT> (Trilinos::K, OverlapLevel);
  List.set("fact: ict level-of-fill", 2.0);
  List.set("fact: drop tolerance", 1e-2);
  ifpackPrec->SetParameters(List);
  ifpackPrec->Initialize();
  ifpackPrec->Compute();
  Solver.SetPrecOperator(&*ifpackPrec);

} // setIFPACKPrec

// ----------------------------------------------------------------------------
/**
 * This routine is to be used with preconditioners such as ILUT which require
 * 1s on the diagonal
 */
void checkDiagonalIsZero()
{
  Epetra_Vector diagonal(*Trilinos::blockMap);
  Trilinos::K->ExtractDiagonalCopy(diagonal);
  bool isZeroDiag = false; //initialize to false
  for (int i = 0; i < diagonal.MyLength(); ++i)
  {
    //if diagonal is 0 change it to be 1
    if (diagonal[i] == 0.0)
    {
      diagonal[i] = 1.0;
      isZeroDiag = true;
    }
  }
  if (isZeroDiag) Trilinos::K->ReplaceDiagonalValues(diagonal);
} // void checkDiagonalIsZero()

// ----------------------------------------------------------------------------
/**
 * To be called within solve-output diagonal from here
 *
 * \param dirW    pass in array with Dirichlet boundary face nodes marked
 * \paramdiagonal diagonal scaling vector need to output to multiply solution by
 */
void constructJacobiScaling(const double *dirW, Epetra_Vector &diagonal)
{
  // Loop over nodes owned by that processor
  //
  for (int i = 0; i < localNodes; ++i) {
    for (int j = 0; j < dof; ++j) {
      int error = diagonal.ReplaceGlobalValue(localToGlobalSorted[i],
                          j, //block offset
                          0, //multivector of 1 vector
                          dirW[i*dof+j]); //value to insert
      if (error != 0) {
        std::cout << "ERROR: Setting Dirichlet diagonal scaling values!" << std::endl;
        exit(1);
      }
    }
  }

  //Extract diagonal of K
  Epetra_Vector Kdiag(*Trilinos::blockMap);
  Trilinos::K->ExtractDiagonalCopy(Kdiag);
  for (int i = 0; i < Kdiag.MyLength(); ++i)
  {
    if (Kdiag[i] == 0.0) Kdiag[i] = 1.0;
    Kdiag[i] = 1 / sqrt(abs(Kdiag[i])); //Jacobi scaling 1 / sqrt|aii|
  }

  //multiply the two vectors together
  diagonal.Multiply(1.0, diagonal, Kdiag, 0.0);

  //Multiply K and Fon the left by diagonal
  //Let diag = W, solving W*K*W*y = W*F, where X = W*y
  Trilinos::K->LeftScale(diagonal);

  // Elem by elem multiplication to support diagonal matrix multiply of vector
  // this -> 0.0*this + 1.0*F*diag
  Trilinos::F->Multiply(1.0, *Trilinos::F, diagonal, 0.0);

  //Right scaling of K need to multiply solution vector as well
  Trilinos::K->RightScale(diagonal);

  //Need to multiply v in vv' by diagonal
  if (coupledBC) {
    for (auto bdryVec : Trilinos::bdryVec_list)
      bdryVec->Multiply(1.0, *bdryVec, diagonal, 0.0);
  }
} // void constructJacobiScaling()

// ----------------------------------------------------------------------------
/**
 * \param  v            coupled boundary vector
 * \param  isCoupledBC  determines if coupled resistance BC is turned on
 */
void trilinos_bc_create_(const std::vector<Array<double>> &v_list, bool &isCoupledBC)
{
  //store as global to determine which matvec multiply to use in solver
  coupledBC = isCoupledBC;

  int error = 0;

  if (isCoupledBC)
  {
    //loop over block rows owned by current proc using localToGlobal index ptr
    for (int i = 0; i < ghostAndLocalNodes; ++i)
    {
      //sum values into v fe case
      int num_global_rows = 1; //number of global block rows put in 1 at a time

      // Loop over all the coupled Neumann BCs
      for (int k = 0; k < v_list.size(); ++k)
      {
        auto v = v_list[k].data();
        auto bdryVec = Trilinos::bdryVec_list[k];
        //set the coupled boundary vector
        error = bdryVec->ReplaceGlobalValues(
                num_global_rows,
                &localToGlobalSorted[i], //global idx id-inputting sorted array
                &dof, //dof values per id pointer to dof
                &v[i*dof]); //values of size dof
        if (error != 0)
        {
          throw std::runtime_error("Setting boundary vector values failed");
        }
      }
    }
  }
} // void trilinos_bc_create_()

// ----------------------------------------------------------------------------
/**
 * free the memory of the global data structures
 */
void trilinos_lhs_free_()
{
  if (MLPrec) {
      MLPrec->DestroyPreconditioner();
      delete MLPrec;
      MLPrec = NULL;
  }
  if (Trilinos::blockMap) {
      delete Trilinos::blockMap;
      Trilinos::blockMap = NULL;
  }
  if (Trilinos::F) {
      delete Trilinos::F;
      Trilinos::F = NULL;
  }
  if (Trilinos::K) {
      delete Trilinos::K;
      Trilinos::K = NULL;
  }
  if (Trilinos::X) {
      delete Trilinos::X;
      Trilinos::X = NULL;
  }
  if (Trilinos::ghostX) {
      delete Trilinos::ghostX;
      Trilinos::ghostX = NULL;
  }
  if (Trilinos::Importer) {
      delete Trilinos::Importer;
      Trilinos::Importer = NULL;
  }
  if (Trilinos::bdryVec_list.size() > 0)
  {
    for (auto bdryVec : Trilinos::bdryVec_list)
      delete bdryVec;
      Trilinos::bdryVec_list.clear();
  }
  if (Trilinos::K_graph) {
      delete Trilinos::K_graph;
      Trilinos::K_graph = NULL;
  }

}

// ----------------------------------------------------------------------------
/**
 * for debugging purposes here are routines to print the matrix and RHS vector
 */
void printMatrixToFile()
{
  std::ofstream Kfile("K.txt");
  Kfile << std::scientific;
  Kfile.precision(17);

  //loop over block rows owned by current proc using localToGlobal index pointer
  for (int i = 0; i < ghostAndLocalNodes; ++i)
  {
    int numEntries = nnzPerRow[i]; //block per of entries per row
    //copy global stiffness values
    int rowDim, numBlockEntries;
    std::vector<int> blockIndices(numEntries);
    //int *blockIndices = new int[numEntries];
    std::vector<int> colDims(numEntries);
    //int * colDim = new int[numEntries];
    Trilinos::K->BeginExtractGlobalBlockRowCopy(localToGlobalUnsorted[i],
                        numEntries, rowDim, numBlockEntries, &blockIndices[0],
                        &colDims[0]);

    for (int j = 0; j < numEntries; ++j)
    {
      std::vector<double> values(dof*dof);
      int sizeofValues = dof*dof;
      int LDA = dof;
      Trilinos::K->ExtractEntryCopy(sizeofValues, &values[0], LDA, false);
      //print returned block
      for (int k = 0; k < dof; ++k)
      {
        for (int l = 0; l < dof; ++l)
        {
          Kfile << values[l*dof + k] << " ";
        }
      }
    }
  }
  Kfile.close();
} // printMatrixToFile()


void printRHSToFile()
{
  std::ofstream Ffile("F.txt");
  Ffile.precision(17);
  std::vector<double> F(Trilinos::F->MyLength());
  //extract copy to print values
  Trilinos::F->ExtractCopy(&F[0], 0);
  //print values to file
  for (int i = 0; i < Trilinos::F->MyLength(); ++i)
    Ffile << F[i] << std::endl; //Jacobi preconditioning on the lefts;
  Ffile.close();
}

// ----------------------------------------------------------------------------
/**
 */
void printSolutionToFile()
{
  std::ofstream Xfile("X.txt");
  Xfile.precision(17);
  std::vector<double> X(Trilinos::X->MyLength());
  //extract copy to print values
  Trilinos::X->ExtractCopy(&X[0], 0);
  //print values to file
  for (int i = 0; i < Trilinos::X->MyLength(); ++i)
    Xfile << X[i] << std::endl; //Jacobi preconditioning on the lefts;
  Xfile.close();
}

/////////////////////////////////////////////////////////////////
//                  T r i l i n o s I m p l                    //
/////////////////////////////////////////////////////////////////

//--------------
// TrilinosImpl 
//--------------
// The TrilinosImpl private class hides Trilinos data structures
// and functions.
//
class TrilinosLinearAlgebra::TrilinosImpl {
  public:
    TrilinosImpl();
    void alloc(ComMod& com_mod, eqType& lEq);
    void assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR);
    void initialize(ComMod& com_mod);
    void solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    void solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res);
    void init_dir_and_coup_neu(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res);
    void set_preconditioner(consts::PreconditionerType preconditioner);

    consts::PreconditionerType preconditioner_;

    /// @brief Local to global mapping
    Vector<int> ltg_;

    /// @brief Factor for Dirichlet BCs
    Array<double> W_;

    /// @brief Residual
    Array<double> R_;
};

TrilinosLinearAlgebra::TrilinosImpl::TrilinosImpl()
{
}

/// @brief Allocate Trilinos arrays.
void TrilinosLinearAlgebra::TrilinosImpl::alloc(ComMod& com_mod, eqType& lEq) 
{
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;
  #ifdef n_debug_alloc
  std::cout << "[TrilinosImpl.alloc] dof: " << dof << std::endl;
  std::cout << "[TrilinosImpl.alloc] tnNo: " << tnNo << std::endl;
  std::cout << "[TrilinosImpl.alloc] gtnNo: " << gtnNo << std::endl;
  std::cout << "[TrilinosImpl.alloc] ltg_.size(): " << ltg_.size() << std::endl;
  #endif

  if (W_.size() != 0) {
    W_.clear();
    R_.clear();
    trilinos_lhs_free_();
  }

  W_.resize(dof,tnNo); 
  R_.resize(dof,tnNo);

  int cpp_index = 1;
  int task_id = com_mod.cm.idcm();

  trilinos_lhs_create(gtnNo, lhs.mynNo, tnNo, lhs.nnz, ltg_, com_mod.ltg, com_mod.rowPtr, 
      com_mod.colPtr, dof, cpp_index, task_id, com_mod.lhs.nFaces);

}

/// @brief Assemble local element arrays.
void TrilinosLinearAlgebra::TrilinosImpl::assemble(ComMod& com_mod, const int num_elem_nodes, const Vector<int>& eqN,
        const Array3<double>& lK, const Array<double>& lR)
{
  trilinos_doassem_(const_cast<int&>(num_elem_nodes), eqN.data(), lK.data(), lR.data());
}

/// @brief Set data for Dirichlet and coupled Neumann boundary conditions.
void TrilinosLinearAlgebra::TrilinosImpl::init_dir_and_coup_neu(ComMod& com_mod, const Vector<int>& incL, const Vector<double>& res)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  int dof = com_mod.dof;
  int gtnNo = com_mod.gtnNo;
  int tnNo = com_mod.tnNo;
  auto& lhs = com_mod.lhs;

  if (lhs.nFaces != 0) {
    for (auto& face : lhs.face) {
      face.incFlag = true;
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      if (incL(faIn) == 0)  {
        lhs.face[faIn].incFlag = false;
      }
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      auto& face = lhs.face[faIn];
      face.coupledFlag = false;
      if (!face.incFlag) {
        continue;
      }

      bool flag = (face.bGrp == BcType::BC_TYPE_Neu);
      if (flag && res(faIn) != 0.0) {
        face.res = res(faIn);
        face.coupledFlag = true;
      }
    }
  }

  W_ = 1.0;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    if (!face.incFlag) {
      continue;
    }

    int faDof = std::min(face.dof,dof);

    if (face.bGrp == BcType::BC_TYPE_Dir) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          W_(i,Ac) = W_(i,Ac) * face.val(i,a);
        }
      }
    }
  }

  std::vector<Array<double>> v_list;
  bool isCoupledBC = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    // Extract face
    auto& face = lhs.face[faIn]; 

    // Create a new array for each face and add it to the list
    v_list.push_back(Array<double>(dof,tnNo));
    auto& v = v_list.back();
    
    if (face.coupledFlag) {
      isCoupledBC = true;
      int faDof = std::min(face.dof,dof);

      // Compute the coupled Neumann BC vector and store it in v
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          v(i,Ac) = v(i,Ac) + sqrt(fabs(res(faIn))) * face.val(i,a);
        }
      }
    }
  }

  // Add the v vectors to global bdryVec_list
  trilinos_bc_create_(v_list, isCoupledBC);

  
}

/// @brief Initialze an array used for something.
void TrilinosLinearAlgebra::TrilinosImpl::initialize(ComMod& com_mod)
{
  #ifdef n_debug_initialize
  std::cout << "[TrilinosImpl] ---------- initialize ---------- " << std::endl;
  std::cout << "[TrilinosImpl.initialize] com_mod.tnNo: " << com_mod.tnNo << std::endl;
  #endif
  ltg_.resize(com_mod.tnNo);

  for (int a = 0; a < com_mod.tnNo; a++) {
    ltg_(com_mod.lhs.map(a)) = com_mod.ltg(a);
  }
}

/// @brief Set the preconditioner.
void TrilinosLinearAlgebra::TrilinosImpl::set_preconditioner(consts::PreconditionerType prec_type)
{
  preconditioner_ = prec_type;
}

/// @brief Solve a system of linear equations assembled by fsils.
void TrilinosLinearAlgebra::TrilinosImpl::solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, 
    const Vector<double>& res)
{
  init_dir_and_coup_neu(com_mod, incL, res);

  auto& Val = com_mod.Val;
  auto& R = com_mod.R;
  int solver_type = static_cast<int>(lEq.ls.LS_type);
  int prec_type = static_cast<int>(preconditioner_);
  #define n_debug_solve
  #ifdef debug_solve
  std::cout << "[TrilinosImpl::solve] ---------- solve ---------- " << std::endl;
  auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
  std::cout << "[TrilinosImpl::solve] solver_type: " << solver_type << std::endl;
  std::cout << "[TrilinosImpl::solve] prec_type: " << prec_name << std::endl;
  std::cout << "[TrilinosImpl::solve] Val.size(): " << Val.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] R_.size(): " << R_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] W_.size(): " << W_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve] Call trilinos_global_solve_ " << std::endl;
  #endif

  if (consts::trilinos_preconditioners.count(preconditioner_) == 0) {
    auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
    throw std::runtime_error("[TrilinosLinearAlgebra::solve] ERROR: '" + prec_name + "' is not a valid Trilinos preconditioner.");
  }

  trilinos_global_solve_(Val.data(), R.data(), R_.data(), W_.data(), lEq.FSILS.RI.fNorm,
      lEq.FSILS.RI.iNorm, lEq.FSILS.RI.itr, lEq.FSILS.RI.callD, lEq.FSILS.RI.dB, lEq.FSILS.RI.suc,
      solver_type, lEq.FSILS.RI.relTol, lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD, prec_type);

  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < com_mod.R.nrows(); i++) {
      com_mod.R(i,a) = R_(i,com_mod.lhs.map(a));
    }
  } 
}

/// @brief Solve a system of linear equations assembled by Trilinos.
void TrilinosLinearAlgebra::TrilinosImpl::solve_assembled(ComMod& com_mod, eqType& lEq, const Vector<int>& incL, const Vector<double>& res)
{
  lEq.FSILS.RI.suc = false; 
  int solver_type = static_cast<int>(lEq.ls.LS_type);
  int prec_type = static_cast<int>(preconditioner_);
  bool assembled = true;
  #define n_debug_solve_assembled
  #ifdef debug_solve_assembled
  auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_); 
  std::cout << "[TrilinosImpl::solve_assembled] ---------- solve_assembled ---------- " << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] solver_type: " << solver_type << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] prec_type: " << prec_name << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] R_.size(): " << R_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] W_.size(): " << W_.size() << std::endl;
  std::cout << "[TrilinosImpl::solve_assembled] lEq.FSILS.RI.suc: " << lEq.FSILS.RI.suc << std::endl;
  #endif

  if (consts::trilinos_preconditioners.count(preconditioner_) == 0) {
    auto prec_name = consts::preconditioner_type_to_name.at(preconditioner_);
    throw std::runtime_error("[TrilinosLinearAlgebra::solve_assembled] ERROR: '" + prec_name + "' is not a valid Trilinos preconditioner.");
  }

  init_dir_and_coup_neu(com_mod, incL, res);

  trilinos_solve_(R_.data(), W_.data(), lEq.FSILS.RI.fNorm, lEq.FSILS.RI.iNorm, 
      lEq.FSILS.RI.itr, lEq.FSILS.RI.callD, lEq.FSILS.RI.dB, lEq.FSILS.RI.suc, 
      solver_type, lEq.FSILS.RI.relTol, lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD, 
      prec_type, assembled);

  for (int a = 0; a < com_mod.tnNo; a++) {
    for (int i = 0; i < com_mod.R.nrows(); i++) {
      com_mod.R(i,a) = R_(i,com_mod.lhs.map(a));
    }
  }

}