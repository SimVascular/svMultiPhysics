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

// The functions defined here replicate the Fortran functions defined in LHSA.f.

#include "lhsa.h"

#include "consts.h"
#include "utils.h"

namespace lhsa_ns {

void add_col(const int tnNo, const int row, const int col, int& mnnzeic, Array<int>& uInd)
{
  int i = -1;

  while (true) {
    i = i + 1;
    if (i == mnnzeic) {
      resiz(tnNo, mnnzeic, uInd);
    }

    // If current entry is zero, then  fill it up
    if (uInd(i,row) == -1) {
       uInd(i,row) = col;
       break;
    }

    // If current entry is still smaller, keep going
    if (col > uInd(i,row)) { 
      continue;
    }

    // If column entry already exists, exit
    if (col == uInd(i,row)) { 
      break;
    }

    // If we are this point, then then the current entry is bigger.
    // Shift all the entries from here to the end of the list. If
    // list is full, we request a larger list, otherwise we shift
    // and add the item at the current entry position.
    //
    if (uInd(mnnzeic-1,row) != -1) {
      resiz(tnNo, mnnzeic, uInd);
    }

    for (int j = mnnzeic - 1; j >= i+1; j--) {
       uInd(j, row) = uInd(j-1, row);
    }
    uInd(i,row) = col;
    break;
  } 
}

/// @brief This subroutine assembles the element stiffness matrix into the 
/// global stiffness matrix (Val sparse matrix formatted as a vector). Also
/// assembles the element residual into the global residual (R).
///
/// Parameters:
///   d - Number of nodes? (eNoN)
///
///
/// Modifies
///   com_mod.R - Residual
///   com_mod.Val - LHS matrix
///
/// Replicates 'SUBROUTINE DOASSEM (d, eqN, lK, lR)'.
//
void do_assem(ComMod& com_mod, const int d, const Vector<int>& eqN, const Array3<double>& lK, const Array<double>& lR)
{
  //std::cout << "[lhs::do_assem] ======================= do_assem ==============" << std::endl;
  auto& R = com_mod.R;
  auto& Val = com_mod.Val;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  //std::cout << "[lhs::do_assem] R.size(): " << R.size() << std::endl;
  //std::cout << "[lhs::do_assem] Val.size(): " << Val.size() << std::endl;

  for (int a = 0; a < d; a++) {
    int rowN = eqN(a);
    if (rowN == -1) {
      continue;
    }

    for (int i = 0; i < R.nrows(); i++) {
      R(i,rowN) = R(i,rowN) + lR(i,a);
    }

    for (int b = 0; b < d; b++) {
      int colN = eqN(b);
      if (colN == -1) {
        continue;
      }

      int left = rowPtr(rowN);
      int right = rowPtr(rowN+1);
      int ptr = (right + left) / 2;

      while (colN != colPtr(ptr)) {
        if (colN > colPtr(ptr)) { 
          left  = ptr;
        } else { 
          right = ptr;
        }
        ptr = (right + left) / 2;
      }

      for (int i = 0; i < Val.nrows(); i++) {
        Val(i,ptr) = Val(i,ptr) + lK(i,a,b);
      }
    }
  }
}

//------
// lhsa
//------
// Create data structure and assembling LHS sparse matrix.
//
// Modifies:
//   com_mod.idMap
//   com_mod.colPtr.resize(nnz); 
//   com_mod.rowPtr.resize(tnNo+1);
//
void lhsa(Simulation* simulation, int& nnz)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& idMap = com_mod.idMap;
  int tnNo = com_mod.tnNo;
  idMap.resize(tnNo);

  for (int a = 0; a < tnNo; a++) { 
    idMap[a] = a;
  }

  int max_enon = std::numeric_limits<int>::min();
  for (auto& mesh : com_mod.msh) { 
    max_enon = std::max(mesh.eNoN, max_enon);
  }

  int mnnzeic = 10 * max_enon; 

  // First fill uInd array depending on mesh connectivity as is.
  //
  // Set uInd = -1 for checking end of rows/columns.
  //
  Array<int> uInd(mnnzeic, tnNo);
  uInd = -1;
  int nMsh = com_mod.nMsh;
  bool risFlag = com_mod.risFlag;

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    if (com_mod.shlEq && msh.eType == ElementType::TRI3) {
      continue;
    }
    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < msh.eNoN; a++) { 
        int rowN = msh.IEN(a,e);
        for (int b = 0; b < msh.eNoN; b++) {
          int colN = msh.IEN(b,e);
          add_col(tnNo, rowN, colN, mnnzeic, uInd);
        }

        // Add extra connections for cooresponding nodes in case of RIS
        int jMRIS = 0;
        if (risFlag) {
          auto& RIS = com_mod.ris;
          // If rowN is in the list of ris nodes
          // Loop through all projections and find the one associated with the current mesh
          for (int iProj = 0; iProj < RIS.nbrRIS; iProj++) {
            if (RIS.lst(0,0,iProj) == iM) {
              jMRIS = 1;
            } else if (RIS.lst(1,0,iProj) == iM) {
              jMRIS = 0;
            } else {
              continue; // Skip to the next iteration
            }
            std::array<int, 2> mapIdx;
            utils::find_loc(com_mod.grisMapList[iProj].map, rowN, mapIdx);
            if (mapIdx[0] != -1) {
              int rowNR = com_mod.grisMapList[iProj].map(jMRIS,mapIdx[1]);
              if (rowNR == -1) continue;
              for (int b = 0; b < msh.eNoN; b++) {
                int colN = msh.IEN(b,e);
                add_col(tnNo, rowNR, colN, mnnzeic, uInd);
              }
            }
          }
        }
      }
    }
  }

  // Treat shells with triangular elements here
  //
  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];

    if (!com_mod.shlEq || !msh.lShl) { 
      continue;
    }

    if (msh.eType != ElementType::TRI3) {
      continue;
    }

    if (msh.eType == ElementType::NRB) {
      continue;
    }

    for (int e = 0; e < msh.nEl; e++) {
      for (int a = 0; a < 2*msh.eNoN; a++) { 
        int rowN;

        if (a < msh.eNoN) {
          rowN = msh.IEN(a,e);
        } else { 
          rowN = msh.eIEN(a-msh.eNoN,e);
        }

        if (rowN == -1) { 
          continue;
        }

        int colN;

        for (int b = 0; b < 2*msh.eNoN; b++) {
          if (b < msh.eNoN) {
            colN = msh.IEN(b,e);
          } else { 
            colN = msh.eIEN(b-msh.eNoN,e);
          }
          if (colN == -1) { 
            continue;
          }

          add_col(tnNo, rowN, colN, mnnzeic, uInd);
        }
      }
    }
  }

  // Now reset idMap for undeforming Neumann BC faces. Then insert
  // master node as a column entry in each row for all the slave nodes.
  // This step is performed even for ghost master nodes where the idMap
  // points to the ghost master node.
  //
  bool flag = false; 

  for (int i = 0; i < com_mod.nEq; i++) {
    auto& eq = com_mod.eq[i];

    for (int j = 0; j <  eq.nBc; j++) {
      auto& bc = eq.bc[j];
      int iM = bc.iM;
      int iFa = bc.iFa;

      if (utils::btest(bc.bType, enum_int(BoundaryConditionType::bType_undefNeu))) {
        int masN = bc.masN;
        if (masN == 0) {
          continue;
        }
        for (int a = 0; a < com_mod.msh[iM].fa[iFa].nNo; a++) {
          int rowN = com_mod.msh[iM].fa[iFa].gN[a];
          if (rowN == masN) {
            continue;
          }
          idMap(rowN) = masN;
          // Insert master to the row if not already present
          add_col(tnNo, rowN, masN, mnnzeic, uInd);
        }

        flag = true; 
      }
    }
  }

  // Change uInd if idMap has been changed
  //
  if (flag) {
    for (int a = 0; a < tnNo; a++) {
      int rowN = idMap(a);

      // If the mapping is not changed, examine the mapping of the
      // column entries of uInd. Don't do anything if the mapping is
      // unchanged. If the mapping is changed, then add to the column\
      // indices of uInd if the entry is not already present.
      //
      if (rowN == a) {
        int i = -1;

        while (true) { 
          i = i + 1;
          int b = uInd(i,rowN);
          // Terminate if end of column entries are reached
          if (b == -1) {
            break;
          }
          int colN = idMap(b);
          // Ignore if the column entry mapping is not changed.
          // This entry is already present and will be used to
          // assemble mass (D) matrix
          if (b == colN) {
            continue;
          }

          // As the column entry is now mapped to a new node,
          // search all column entries and insert the new node if
          // it is not present. This step is performed to assemble
          // the Divergence (C) matrix
          add_col(tnNo, rowN, colN, mnnzeic, uInd);
          if (i == mnnzeic-1) {
            break;
          }
        }

      // If the row mapping is changed, insert the mapped/unmapped
      // column entries of the old row into the new row if not
      // already present.

      } else { 
        int i = -1;
        while (true) { 
          i = i + 1;
          int b = uInd(i,a);
          // Terminate if end of column entries are reached
          if (b == -1) {
            break;
          }

          // Add unmapped column to assemble gradient matrix
          add_col(tnNo, rowN, b, mnnzeic, uInd);

          // If column is mapped, add the mapped column to assemble
          // stiffness matrix
          int colN = idMap(b);
          if (b != colN) {
            add_col(tnNo, rowN, colN, mnnzeic, uInd);
          }
          if (i == mnnzeic-1) {
            break;
          }
        }
      }
    }
  }

  // Finding number of non-zeros in colPtr vector
  nnz = 0;
  for (int rowN = 0; rowN < tnNo; rowN++) { 
    if (uInd(0,rowN) == -1) {
      throw std::runtime_error("An isolated node " + std::to_string(rowN) + " has been found during assambly");
    }
    for (int i = 0; i < mnnzeic; i++) {
      if (uInd(i,rowN) != -1) {
        nnz = nnz + 1;
      }
    }
  }

  // Now constructing compact form of rowPtr and colPtr
  //
  com_mod.colPtr.resize(nnz); 
  com_mod.rowPtr.resize(tnNo+1);
  int j = 0;
  com_mod.rowPtr(0) = 0;

  for (int rowN = 0; rowN < tnNo; rowN++) { 
    for (int i = 0; i < mnnzeic; i++) {
      if (uInd(i, rowN) != -1) {
        com_mod.colPtr(j) = uInd(i, rowN);
        j = j + 1;
      }
    }
    com_mod.rowPtr(rowN+1) = j;
  }
}

//-------
// resiz
//-------
//
void resiz(const int tnNo, int& mnnzeic, Array<int>& uInd)
{
  int n = mnnzeic;
  Array<int> tmp(n, tnNo);
  tmp = -1;

  for (int i = 0; i < uInd.nrows(); i++) {
    for (int j = 0; j < uInd.ncols(); j++) {
      tmp(i,j) = uInd(i,j);
    }
  }

  mnnzeic = n + std::max(5, n/5);
  uInd.resize(mnnzeic, tnNo);
  uInd = -1;

  for (int i = 0; i < tmp.nrows(); i++) {
    for (int j = 0; j < tmp.ncols(); j++) {
      uInd(i,j) = tmp(i,j);
    }
  }
}

};


