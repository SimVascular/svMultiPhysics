// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// This routine computes body force for the current equation and assembles it to the residual

#include "bf.h"

#include "all_fun.h"
#include "cmm.h"
#include "consts.h"
#include "fft.h"
#include "lhsa.h"
#include "shells.h"
#include "utils.h"

namespace bf {

/// @brief This subroutine is reached only for shell follower pressre loads
/// or applying initialization pressure for CMM method. nsd must be
/// equal to 3.
///
/// Reproduces 'SUBROUTINE BFCONSTRUCT(lM, e, eNoN, idof, xl, dl, bfl, ptr)'.
//
void bf_construct(ComMod& com_mod, const mshType& lM, const int e, const int eNoN, const int idof, Array<double>& xl, 
    const Array<double>& dl, const Array<double>& bfl, const Vector<int>& ptr)
{
  using namespace consts;

  #define n_debug_bf_construct 
  auto& cm = com_mod.cm;
  #ifdef debug_bf_construct 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;

  auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  Vector<double> N(eNoN); 
  Array<double> Nx(2,eNoN);
  Array<double> lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);

  cDmn = all_fun::domain(com_mod, lM, cEq, e);
  auto cPhys = eq.dmn[cDmn].phys;
  #ifdef debug_bf_construct 
  dmsg << "cPhys: " << cPhys;
  #endif

  //  Updating the shape functions, if neccessary
  if (lM.eType == ElementType::NRB) {
    //CALL NRBNNX(lM, e)
  }

  // Setting intial values
  //
  for (int g = 0; g < lM.nG; g++) {
    double w = lM.w(g);
    Vector<double> N = lM.N.col(g);
    Array<double> Nx = lM.Nx.slice(g);

    switch (cPhys) {
      case EquationType::phys_shell:
        // [NOTE] passing Array 'bfl' to Vector 'tfl' arg in shell_fp does not work.
        shells::shell_fp(com_mod, eNoN, w, N, Nx, dl, xl, bfl, lR, lK);
        //CALL SHELLFP(eNoN, w, N, Nx, dl, xl, bfl, lR, lK)
        //throw std::runtime_error("[bf_construct] Shell follower pressure loads not implemented.");
      break;

      case EquationType::phys_CMM:
        cmm::bcmmi(com_mod, eNoN, idof, w, N, Nx, xl, bfl, lR);
      break;

      default: 
        throw std::runtime_error("[bf_construct] Undefined physics.");
    }
  }

  // Now doing the assembly part
  eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
}

/// @brief Modifes: com_mod.Bf, Dg
//
void set_bf(ComMod& com_mod, const Array<double>& Dg)
{
  #define n_debug_set_bf
  auto& cm = com_mod.cm;
  #ifdef debug_set_bf 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  auto& Bf = com_mod.Bf;
  Bf = 0.0;
  #ifdef debug_set_bf 
  dmsg << "eq.nBf: " << eq.nBf;
  #endif

  for (int iBf = 0; iBf < eq.nBf; iBf++) {
    int iM = eq.bf[iBf].iM;
    set_bf_l(com_mod, eq.bf[iBf], com_mod.msh[iM], Dg);
  }
}

/// @brief Modifies
/// \code {.cpp}
///  com_mod.Bf
/// \endcode
///
/// Reproduces 'SUBROUTINE SETBFL(lBf, lM, Dg)'.
//
void set_bf_l(ComMod& com_mod, bfType& lBf, mshType& lM, const Array<double>& Dg)
{
  #define n_debug_set_bf_l 
  auto& cm = com_mod.cm;
  #ifdef debug_set_bf_l 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lBf.file_name: " << lBf.file_name;
  #endif

  using namespace consts;

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const int tDof = com_mod.tDof;
  const int nNo = lM.nNo;
  const int idof = lBf.dof;
  const int eNoN = lM.eNoN;
  #ifdef debug_set_bf_l 
  dmsg << "nsd: " << nsd;
  dmsg << "idof: " << idof;
  dmsg << "eNoN: " << eNoN;
  dmsg << "lBf.bType: " << lBf.bType;
  #endif

  Vector<double> f(idof);
  Array<double> bfl; 

  if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_std))) {
    f = lBf.b;

  } else if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_ustd))) {
    Vector<double> rtmp(1);
    ifft(com_mod, lBf.bt, f, rtmp);

  } else if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_gen))) {
    bfl.resize(idof,nNo); 
    Array<double> xl(idof,nNo);
    igbc(com_mod, lBf.bm, bfl, xl);
  }

  Array<double> bfg(idof,tnNo);

  if (utils::btest(lBf.bType,enum_int(BodyForceType::bfType_gen))) {
    for (int a = 0; a < lM.nNo; a++) {
      int Ac = lM.gN(a);
      bfg.set_col(Ac, bfl.col(a));
    }

  } else if (utils::btest(lBf.bType,enum_int(BodyForceType::bfType_spl))) {
    for (int a = 0; a < lM.nNo; a++) {
      int Ac = lM.gN(a);
      bfg.set_col(Ac, lBf.bx.col(a));
    }

  } else {
    for (int a = 0; a < lM.nNo; a++) {
      int Ac = lM.gN(a);
      bfg.set_col(Ac, f);
    }
  }

  // Assemble pressure/traction load (shells/CMM initialization) to
  // residual. For general body force (vector), assemble later with
  // other volumetric forces
  //
  if (utils::btest(lBf.bType, enum_int(BodyForceType::bfType_vol))) {
    for (int a = 0; a < nNo; a++) {
      int Ac = lM.gN(a);
      com_mod.Bf.set_col(Ac, bfg.col(Ac));
     }

  } else {
    Array<double> bfl(idof,eNoN); 
    Array<double> xl(nsd,eNoN); 
    Array<double> dl(tDof,eNoN);
    Vector<int> ptr(eNoN);

    for (int e = 0; e < lM.nEl; e++) {
      for (int a = 0; a < eNoN; a++) {
        int Ac = lM.IEN(a,e);
        ptr(a) = Ac;
        xl.set_col(a, com_mod.x.col(Ac));
        dl.set_col(a, Dg.col(Ac));
        bfl.set_col(a, bfg.col(Ac));
      }
      bf_construct(com_mod, lM, e, eNoN, idof, xl, dl, bfl, ptr);
    }
  }
}

};


