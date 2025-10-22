// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "mesh.h"

#include "all_fun.h"
#include "consts.h"
#include "fluid.h"
#include "fs.h"
#include "l_elas.h"
#include "lhsa.h"
#include "nn.h"
#include "sv_struct.h"
#include "utils.h"

#include <array>
#include <iomanip>
#include <math.h>

namespace mesh {

void construct_mesh(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag, const Array<double>& Dg)
{
  #define n_debug_construct_mesh 
  #ifdef debug_construct_mesh
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  auto& cem = cep_mod.cem;
  const int nsd  = com_mod.nsd;
  const int tDof = com_mod.tDof;
  const int dof = com_mod.dof;
  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;
  const int nsymd = com_mod.nsymd;
  auto& Do = com_mod.Do;
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;
  auto& pSa = com_mod.pSa;
  bool pstEq = com_mod.pstEq;

  // Start and end DOF
  int is = nsd + 1;
  int ie = 2*nsd;
  int eNoN = lM.eNoN;
  #ifdef debug_construct_mesh
  dmsg << "cEq: " << cEq;
  dmsg << "is: " << is;
  dmsg << "ie: " << ie;
  dmsg << "cDmn: " << cDmn;
  #endif

  Vector<int> ptr(eNoN);
  Vector<double> pSl(nsymd), ya_l(eNoN), N(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), dl(tDof,eNoN),
                dol(nsd,eNoN), pS0l(nsymd,eNoN), Nx(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);
  Array<double> ksix(nsd,nsd), bfl(nsd,eNoN);

  for (int e = 0; e < lM.nEl; e++) {
    // Update domain and proceed if domain phys and eqn phys match
    cDmn = all_fun::domain(com_mod, lM, cEq, e);
    auto cPhys = eq.dmn[cDmn].phys;
    if (cPhys != EquationType::phys_mesh) {
      continue;
    }

    // Update shape functions for NURBS
    if (lM.eType == ElementType::NRB) {
      //CALL NRBNNX(lM, e)
    }

    // Create local copies
    bfl = 0.0;
    pS0l = 0.0;

    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.IEN(a,e);
      ptr(a) = Ac;

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
        dol(i,a) = Do(is+i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        al(i,a) = Ag(i,Ac);
        dl(i,a) = Dg(i,Ac);
      }
    }

    // For MESH, the reference configuration is the one at the
    // beginning of the time step. Update displacements accordingly
    for (int i = 0; i < nsd; i++) {
      for (int j = 0; j < eNoN; j++) {
        xl(i,j) = xl(i,j) + dol(i,j);
        dl(i+is,j) = dl(i+is,j) - dol(i,j);
      }
    }

    // Gauss integration
    //
    lR = 0.0;
    lK = 0.0;
    double Jac{0.0};

    for (int g = 0; g < lM.nG; g++) {
      if (g == 0 || !lM.lShpF) {
        auto Nx_g = lM.Nx.rslice(g);
        nn::gnn(eNoN, nsd, nsd, Nx_g, xl, Nx, Jac, ksix);
        if (utils::is_zero(Jac)) {
          throw std::runtime_error("[construct_mesh] Jacobian for element " + std::to_string(e) + " is < 0.");
        }
      }

      double w = lM.w(g);
      N = lM.N.col(g);
      pS0l = 0.0;

      if (nsd == 3) {
        l_elas::l_elas_3d(com_mod, eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR, lK);

      } else if (nsd == 2) {
        l_elas::l_elas_2d(com_mod, eNoN, w, N, Nx, al, dl, bfl, pS0l, pSl, lR, lK);
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

};
