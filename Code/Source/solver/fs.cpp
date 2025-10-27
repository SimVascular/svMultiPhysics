// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// Replicates Fortran FS.f.

#include "fs.h"
#include "consts.h"
#include "nn.h"

namespace fs {


/// @brief Allocates arrays within the function space type. Assumes that 
/// fs%eNoN and fs%nG are already defined
///
/// Replicates 'SUBROUTINE ALLOCFS(fs, insd)'.
//
void alloc_fs(fsType& fs, const int nsd, const int insd)
{
  int nG = fs.nG;
  int eNoN = fs.eNoN;

  fs.w.resize(nG); 
  fs.xi.resize(insd,nG); 
  fs.xib.resize(2,nsd); 
  fs.N.resize(eNoN,nG); 
  fs.Nb.resize(2,eNoN); 
  fs.Nx.resize(insd,eNoN,nG);

  int ind2 = std::max(3*(insd-1), 1);
  fs.Nxx.resize(ind2,eNoN,nG);
}


void get_thood_fs(ComMod& com_mod, std::array<fsType,2>& fs, const mshType& lM, const bool lStab, const int iOpt)
{
  #define n_debug_get_thood_fs 
  #ifdef debug_get_thood_fs 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lStab: " << lStab;
  dmsg << "iOpt: " << iOpt;
  #endif

  int nsd = com_mod.nsd;

  if (lStab) {
    for (int i = 0; i < 2; i++) {
      fs[i].nG = lM.fs[0].nG;
      fs[i].eType = lM.fs[0].eType;
      fs[i].lShpF = lM.fs[0].lShpF;
      fs[i].eNoN  = lM.fs[0].eNoN;

      alloc_fs(fs[i], nsd, nsd);

      fs[i].w   = lM.fs[0].w;
      fs[i].xi  = lM.fs[0].xi;
      fs[i].N   = lM.fs[0].N;
      fs[i].Nx  = lM.fs[0].Nx;
      fs[i].xib = lM.fs[0].xib;
      fs[i].Nb  = lM.fs[0].Nb;

      // [NOTE] Not sure what this is all about. 
      /*
      if (fs[i].Nxx.size() != 0) {
        fs[i].Nxx = lM.fs[0].Nxx;
      }
      */
      fs[i].Nxx = lM.fs[0].Nxx;
    }

  } else { 
    if (iOpt == 1) {
      fs[0].nG = lM.fs[0].nG;
      fs[0].eType = lM.fs[0].eType;
      fs[0].lShpF = lM.fs[0].lShpF;
      fs[0].eNoN = lM.fs[0].eNoN;

      alloc_fs(fs[0], nsd, nsd);

      fs[0].w = lM.fs[0].w;
      fs[0].xi = lM.fs[0].xi;
      fs[0].N = lM.fs[0].N;
      fs[0].Nx = lM.fs[0].Nx;
      fs[0].xib = lM.fs[0].xib;
      fs[0].Nb = lM.fs[0].Nb;

      if (fs[0].Nxx.size() != 0) {
        fs[0].Nxx = lM.fs[0].Nxx;
      }

      fs[1].nG = lM.fs[0].nG;
      fs[1].eType = lM.fs[1].eType;
      fs[1].lShpF = lM.fs[1].lShpF;
      fs[1].eNoN = lM.fs[1].eNoN;

      alloc_fs(fs[1], nsd, nsd);

      fs[1].w = lM.fs[0].w;
      fs[1].xi = lM.fs[0].xi;

      for (int g = 0; g < fs[1].nG; g++) {
        nn::get_gnn(nsd, fs[1].eType, fs[1].eNoN, g, fs[1].xi, fs[1].N, fs[1].Nx);
      }
      nn::get_nn_bnds(nsd, fs[1].eType, fs[1].eNoN, fs[1].xib, fs[1].Nb);

    } else if (iOpt == 2) {
      fs[1].nG    = lM.fs[1].nG;
      fs[1].eType = lM.fs[1].eType;
      fs[1].lShpF = lM.fs[1].lShpF;
      fs[1].eNoN  = lM.fs[1].eNoN;

      alloc_fs(fs[1], nsd, nsd);

      fs[1].w   = lM.fs[1].w;
      fs[1].xi  = lM.fs[1].xi;
      fs[1].N   = lM.fs[1].N;
      fs[1].Nx  = lM.fs[1].Nx;
      fs[1].xib = lM.fs[1].xib;
      fs[1].Nb  = lM.fs[1].Nb;

      fs[0].nG = lM.fs[1].nG;
      fs[0].eType = lM.fs[0].eType;
      fs[0].lShpF = lM.fs[0].lShpF;
      fs[0].eNoN  = lM.fs[0].eNoN;

      alloc_fs(fs[0], nsd, nsd);

      fs[0].w = lM.fs[1].w;
      fs[0].xi = lM.fs[1].xi;

      for (int g = 0; g < fs[0].nG; g++) {
        nn::get_gnn(nsd, fs[0].eType, fs[0].eNoN, g, fs[0].xi, fs[0].N, fs[0].Nx);
      }
      nn::get_nn_bnds(nsd, fs[0].eType, fs[0].eNoN, fs[0].xib, fs[0].Nb);
    }
  }
}


void init_fs(fsType& fs,  const int nsd, const int insd)
{
  alloc_fs(fs, nsd, insd);

  // Get Gauss points and shape functions
  //
  // xi(insd,nG)
  //
  // w(nG)
  //
  nn::get_gip(insd, fs.eType, fs.nG, fs.w, fs.xi);

  for (int g = 0; g < fs.nG; g++) {
    nn::get_gnn(insd, fs.eType, fs.eNoN, g, fs.xi, fs.N, fs.Nx);
  }

  nn::get_nn_bnds(nsd, fs.eType, fs.eNoN, fs.xib, fs.Nb);
}


void init_fs_face(const ComMod& com_mod, mshType& lM, faceType& lFa)
{
  using namespace consts;

  #define n_debug_init_fs_face
  #ifdef debug_init_fs_face
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  int nsd = com_mod.nsd;
  int insd = nsd - 1;

  if (lM.lShl) {
    insd = nsd - 1;
  }

  if (lM.lFib) {
    insd = 0;
  }

  lFa.nFs = lM.nFs;
  lFa.fs.resize(lFa.nFs);

  // The first set of basis is inherited directly from face basis
  lFa.fs[0].lShpF = lM.lShpF;
  lFa.fs[0].eType = lFa.eType;
  lFa.fs[0].eNoN  = lFa.eNoN;
  lFa.fs[0].nG    = lFa.nG;

  alloc_fs(lFa.fs[0], nsd, insd);

  if (lFa.eType != ElementType::NRB) {
     lFa.fs[0].w   = lFa.w;
     lFa.fs[0].xi  = lFa.xi;
     lFa.fs[0].N   = lFa.N;
     lFa.fs[0].Nx  = lFa.Nx;
     nn::get_nn_bnds(nsd, lFa.fs[0].eType, lFa.fs[0].eNoN, lFa.fs[0].xib, lFa.fs[0].Nb);
  }

  // Sets Taylor-Hood basis if invoked by user (fluid, ustruct, FSI)
  if (lFa.nFs == 2) {
     // Select Taylor-Hood element
     set_thood_fs(lFa.fs[1], lFa.fs[0].eType);

     // Allocate arrays
     alloc_fs(lFa.fs[1], nsd, insd);

     // Get Gauss points and shape functions
     //
     nn::get_gip(insd, lFa.fs[1].eType, lFa.fs[1].nG, lFa.fs[1].w, lFa.fs[1].xi);

     for (int g = 0; g < lFa.fs[1].nG; g++) {
        nn::get_gnn(insd, lFa.fs[1].eType, lFa.fs[1].eNoN, g, lFa.fs[1].xi, lFa.fs[1].N, lFa.fs[1].Nx);
     }

     nn::get_nn_bnds(nsd, lFa.fs[1].eType, lFa.fs[1].eNoN, lFa.fs[1].xib, lFa.fs[1].Nb);
  }
}


/// @brief Modifies:
///  lM.fs.resize(lM.nFs)
///  lM.fs[0].lShpF = lM.lShpF;
///  lM.fs[0].eType = lM.eType;
///  lM.fs[0].eNoN  = lM.eNoN;
///  lM.fs[0].nG    = lM.nG;
///
/// Replicates 'SUBROUTINE INITFSMSH(lM)'.
//
void init_fs_msh(const ComMod& com_mod, mshType& lM)
{
  using namespace consts;

  #define n_debug_init_fs_msh
  #ifdef debug_init_fs_msh
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.nFs: " << lM.nFs;
  dmsg << "lM.eNoN: " << lM.eNoN;
  dmsg << "lM.nG: " << lM.nG;
  dmsg << "lM.lShpF: " << lM.lShpF;
  #endif

  int nsd = com_mod.nsd;
  int insd = nsd;

  if (lM.lShl) {
    insd = nsd - 1;
  }

  if (lM.lFib) {
    insd = 1;
  }
  #ifdef debug_init_fs_msh
  dmsg << "nsd: " << nsd;
  dmsg << "insd: " << insd;
  #endif

  lM.fs.resize(lM.nFs);

  // The first set of basis is inherited directly from mesh basis
  lM.fs[0].lShpF = lM.lShpF;
  lM.fs[0].eType = lM.eType;
  lM.fs[0].eNoN  = lM.eNoN;
  lM.fs[0].nG    = lM.nG;

  alloc_fs(lM.fs[0], nsd, insd);

  if (lM.eType != ElementType::NRB) {
    lM.fs[0].w = lM.w;
    lM.fs[0].xi  = lM.xi;
    lM.fs[0].xib = lM.xib;
    lM.fs[0].N   = lM.N;
    lM.fs[0].Nb  = lM.Nb;
    lM.fs[0].Nx  = lM.Nx;
  }
  // Second order derivatives for vector function space
  //
  if (!lM.fs[0].lShpF) {
    int ind2 = std::max(3*(insd-1), 1);
    for (int g = 0; g < lM.fs[0].nG; g++) {
      nn::get_gn_nxx(insd, ind2, lM.fs[0].eType, lM.fs[0].eNoN, g, lM.fs[0].xi, lM.fs[0].Nxx);
    }
  }

  // Sets Taylor-Hood basis [fluid, stokes, ustruct, FSI)
  if (lM.nFs == 2) {
    // Select Taylor-Hood element
    set_thood_fs(lM.fs[1], lM.fs[0].eType);

    // Initialize the function space
    init_fs(lM.fs[1], nsd, insd);
  }
}


/// @brief Sets Tayloor-Hood basis for a parent element type
///
/// Replicates 'SUBROUTINE SETTHOODFS(fs, eType)'.
//
void set_thood_fs(fsType& fs, consts::ElementType eType)
{
  using namespace consts;

  switch (eType) {

    case ElementType::HEX20:
    case ElementType::HEX27:
      fs.eType = ElementType::HEX8;
      fs.lShpF = false; 
      fs.eNoN  = 8;
      fs.nG    = 8;
    break;

    case ElementType::LIN2:
      fs.eType = ElementType::LIN1;
      fs.lShpF = true; 
      fs.eNoN  = 2;
      fs.nG    = 2;
    break;

    case ElementType::QUD8:
    case ElementType::QUD9:
      fs.eType = ElementType::QUD4;
      fs.lShpF = false; 
      fs.eNoN  = 4;
      fs.nG    = 4;
    break;

    case ElementType::TET10:
      fs.eType = ElementType::TET4;
      fs.lShpF = true; 
      fs.eNoN  = 4;
      fs.nG    = 4;
    break;

    case ElementType::TRI6:
      fs.eType = ElementType::TRI3;
      fs.lShpF = true; 
      fs.eNoN  = 3;
      fs.nG    = 3;
    break;

    default:
      throw std::runtime_error("Cannot choose Taylor-Hood basis");
    break;
  }
}

/// @brief Set the residual of the continuity equation and its tangent matrix
/// due to variation with pressure to 0 on all the edge nodes. This
/// step is done only for P2P1 type discretization for mixed saddle
/// point type problems such as fluid, stokes, ustruct, and fsi.
///
/// Modifies: com_mod.Val, com_mod.R
///
/// Reproduces 'SUBROUTINE THOOD_ValRC()'
//
void thood_val_rc(ComMod& com_mod)
{
  using namespace consts;

  #define n_debug_thood_val_rc 
  #ifdef debug_thood_val_rc 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  if (std::set<EquationType>{Equation_stokes, Equation_fluid, Equation_ustruct, Equation_FSI}.count(eq.phys) == 0) {
    return;
  }

  const int tnNo = com_mod.tnNo;
  const int nsd = com_mod.nsd;
  const auto& rowPtr = com_mod.rowPtr;
  const auto& colPtr = com_mod.colPtr;
  auto& Val = com_mod.Val;
  auto& R = com_mod.R;
     
  bool THflag = false; 

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    if (com_mod.msh[iM].nFs == 2) {
      THflag = true;
      break;
    }
  }

  #ifdef debug_thood_val_rc 
  dmsg << "THflag: " << THflag;
  #endif

  if (THflag) {
    Vector<int> eNds(tnNo);

    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];
      if (msh.nFs == 1) {
        continue;
      }

      int i = msh.fs[1].eNoN;

      for (int e = 0; e < msh.nEl; e++) {
        for (int a = i; a < msh.fs[0].eNoN; a++) {
          int Ac = msh.IEN(a,e);
          eNds(Ac) = 1;
        }
      }
    }

    for (int a = 0; a < tnNo; a++) {
      if (eNds(a) == 1) {
        R(nsd,a) = 0.0;
        int s = (nsd+1)*(nsd+1) - 1;

        for (int i = rowPtr(a); i <= rowPtr(a+1)-1; i++) {
          int c = colPtr(i);
          if (c == a) {
            Val(s,i) = 1.0;
          } else {
            Val(s,i) = 0.0;
          }
        }
      }
    }
  }
}

};
