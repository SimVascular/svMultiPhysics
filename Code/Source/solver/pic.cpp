// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The code here replicates the Fortran code in PIC.f.
//
// See the publications below, section 4.4 for theory and derivation:
//  1.  Bazilevs, et al. "Isogeometric fluid-structure interaction:
//      theory, algorithms, and computations.", Computational Mechanics,
//      43 (2008): 3-37. doi: 10.1007/s00466-008-0315-x
//  2. Bazilevs, et al. "Variational multiscale residual-based 
//      turbulence modeling for large eddy simulation of incompressible 
//      flows.", CMAME (2007)

#include "pic.h"

#include "Simulation.h"
#include "all_fun.h"
#include "cep_ion.h"
#include "nn.h"
#include "utils.h"

#include "mpi.h"

#include <iostream>
#include <math.h>

namespace pic {

/// @brief This is the corrector. Decision for next eqn is also made here (modifies cEq global).
///
/// Modifies:
/// \code {.cpp}
///   com_mod.Ad
///   com_mod.An
///   com_mod.Dn
///   com_mod.Yn
///   cep_mod.Xion
///   com_mod.pS0
///   com_mod.pSa
///   com_mod.pSn
///
///   com_mod.cEq
///   eq.FSILS.RI.iNorm 
///   eq.pNorm
/// \endcode
//
void picc(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cep_mod = simulation->get_cep_mod();

  #define n_debug_picc
  #ifdef debug_picc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const double dt = com_mod.dt;

  const auto& R = com_mod.R;
  const auto& Rd = com_mod.Rd;

  auto& cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];

  auto& An = com_mod.An;
  auto& Ad = com_mod.Ad;
  auto& Dn = com_mod.Dn;
  auto& Yn = com_mod.Yn;

  auto& pS0 = com_mod.pS0;
  auto& pSa = com_mod.pSa;
  auto& pSn = com_mod.pSn;
  auto& Xion = cep_mod.Xion;

  int s = eq.s;
  int e = eq.e;

  std::array<double,4> coef;
  coef[0] = eq.gam * dt;
  coef[1] = eq.beta*dt*dt;
  coef[2] = 1.0 / eq.am;
  coef[3] = eq.af*coef[0]*coef[2];

  #ifdef debug_picc
  dmsg << "cEq: " << cEq;
  dmsg << "s: " << s;
  dmsg << "e: " << e;
  dmsg << "coef: " << coef[0] << " " << coef[1] << " " << coef[2] << " " << coef[3];
  dmsg << "sstEq: " << com_mod.sstEq;
  dmsg << "An nrows: " << An.nrows_;
  dmsg << "   ncols: " << An.ncols_;
  #endif

  // ustruct, FSI (ustruct)
  //
  if (com_mod.sstEq) {
    if (eq.phys == EquationType::phys_ustruct || eq.phys == EquationType::phys_FSI) {
      Vector<double> dUl(nsd);

      for (int a = 0; a < tnNo; a++) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,a) = An(i+s,a) - R(i,a);
          Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];
        }

        for (int i = 0; i < e-s; i++) {
          dUl(i) = Rd(i,a)*coef[2] + R(i,a)*coef[3];
          Ad(i,a) = Ad(i,a) - dUl(i);
          Dn(i+s,a) = Dn(i+s,a) - dUl(i)*coef[0];
        }
      }

    } else if (eq.phys == EquationType::phys_mesh) {
      for (int a = 0; a < tnNo; a++) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,a) = An(i+s,a) - R(i,a);
          Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];
          Dn(i+s,a) = Dn(i+s,a) - R(i,a)*coef[1];
        }
      }
    }

  } else {
    for (int a = 0; a < tnNo; a++) {
      for (int i = 0; i < e-s+1; i++) {
        // eqn 94 of Bazilevs 2007 // here, -R contains the acceleration update (obtained from Newton solve))?
        An(i+s,a) = An(i+s,a) - R(i,a);
        
        // eqn 95 of Bazilevs 2007
        Yn(i+s,a) = Yn(i+s,a) - R(i,a)*coef[0];
        
        Dn(i+s,a) = Dn(i+s,a) - R(i,a)*coef[1];
      }
    }
  }

  if (std::set<EquationType>{Equation_stokes, Equation_fluid, Equation_ustruct, Equation_FSI}.count(eq.phys) != 0) {
    pic_eth(simulation);
  }

  // If explicit geometric coupling is not used, update the acceleration, 
  // velocity, and displacement for the FSI equation as usual
  if (eq.phys == Equation_FSI && !eq.expl_geom_cpl) {
    if (com_mod.eq[1].phys != Equation_mesh) {
      throw std::runtime_error("PICC: FSI simulation requires a mesh motion equation as the second equation.");
    }
    int s = com_mod.eq[1].s;
    int e = com_mod.eq[1].e;
    #ifdef debug_picc
    dmsg << "eq.phys == Equation_FSI ";
    dmsg << "com_mod.eq[1].sym: " << com_mod.eq[1].sym;
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    #endif

    for (int Ac = 0; Ac < tnNo; Ac++) {
      if (all_fun::is_domain(com_mod, eq, Ac, Equation_struct) || 
          all_fun::is_domain(com_mod, eq, Ac, Equation_ustruct) || 
          all_fun::is_domain(com_mod, eq, Ac, Equation_lElas)) {
        for (int i = 0; i < e-s+1; i++) {
          An(i+s,Ac) = An(i,Ac);
          Yn(i+s,Ac) = Yn(i,Ac);
          Dn(i+s,Ac) = Dn(i,Ac);
        }
      }
    }
  }

  // Update Xion for cardiac electrophysiology
  //  
  if (eq.phys == Equation_CEP) {
    int s = eq.s;
    for (int a = 0; a < tnNo; a++) {
      Xion(0,a) = Yn(s,a);
    }
  }

  // Update prestress at the nodes and re-initialize
  //
  if (com_mod.pstEq) {
    all_fun::commu(com_mod, pSn);
    all_fun::commu(com_mod, pSa);

    for (int a = 0; a < tnNo; a++) {
      if (!utils::is_zero(pSa(a))) {
        for (int i = 0; i < pSn.nrows(); i++) {
          pSn(i,a) = pSn(i,a) / pSa(a);
        }
      }
    }

    pSa = 0.0;
  }

  // Filter out the non-wall displacements for CMM equation
  //
  if (eq.phys == Equation_CMM && !com_mod.cmmInit) {
    for (int a = 0; a < tnNo; a++) {
      double r1 = static_cast<double>(com_mod.cmmBdry(a));
      for (int i = 0; i < e-s; i++) {
        Dn(i+s,a) = Dn(i+s,a)*r1;
      }
    }
  }

  // IB treatment
  //if (ibFlag) CALL IB_PICC()

  // Computes norms and check for convergence of Newton iterations
  double eps = std::numeric_limits<double>::epsilon();

  if (utils::is_zero(eq.FSILS.RI.iNorm)) {
    eq.FSILS.RI.iNorm = eps;
  }

  if (utils::is_zero(eq.iNorm)) {
    eq.iNorm = eq.FSILS.RI.iNorm;
    #ifdef debug_picc
    dmsg << "eq.iNorm: " << eq.iNorm;
    #endif
  }

  if (eq.itr == 1) {
     eq.pNorm = eq.FSILS.RI.iNorm / eq.iNorm;
    #ifdef debug_picc
    dmsg << "eq.itr: " << eq.itr;
    dmsg << "eq.pNorm: " << eq.pNorm;
    #endif
  }

  double r1 = eq.FSILS.RI.iNorm / eq.iNorm;
  bool l1 = (eq.itr >= eq.maxItr);
  bool l2 = (r1 <= eq.tol);
  bool l3 = (r1 <= eq.tol*eq.pNorm);
  bool l4 = (eq.itr >= eq.minItr);

  #ifdef debug_picc
  dmsg << "eq.itr: " << eq.itr;
  dmsg << "eq.minItr: " << eq.minItr;
  dmsg << "r1: " << r1;
  dmsg << "l1: " << l1;
  dmsg << "l2: " << l2;
  dmsg << "l3: " << l3;
  dmsg << "l4: " << l4;
  #endif

  if (l1 || ((l2 || l3) && l4)) {
    eq.ok = true;
    #ifdef debug_picc
    dmsg << "eq.ok: " << eq.ok;
    dmsg << "com_mod.eq[0].ok: " << com_mod.eq[0].ok;
    dmsg << "com_mod.eq[1].ok: " << com_mod.eq[1].ok;
    #endif
  }

  auto& eqs = com_mod.eq;
  if (std::count_if(eqs.begin(),eqs.end(),[](eqType& eq){return eq.ok;}) == eqs.size()) {
    #ifdef debug_picc
    dmsg << "all ok";
    #endif
    return;
  }
  //if (ALL(eq.ok)) RETURN

  if (eq.coupled) {
    if (!eq.expl_geom_cpl) {
      // For coupled equations, if explicit geometric coupling is not used, 
      // increment the equation counter after each linear solve
      cEq = cEq + 1;
      #ifdef debug_picc
      dmsg << "eq " << " coupled ";
      dmsg << "1st update cEq: " << cEq;
      #endif

      auto& eqs = com_mod.eq;
      if (std::count_if(eqs.begin(),eqs.end(),[](eqType& eq){return !eq.coupled || eq.ok;}) == eqs.size()) {
        while (cEq < com_mod.nEq) {
          if (!eqs[cEq].coupled) {
            break;
          }
          cEq = cEq + 1;
        }

      } else {
        if (cEq >= com_mod.nEq) {
          cEq = 0;
        }

        while (!eqs[cEq].coupled) {
          cEq = cEq + 1;
          if (cEq >= com_mod.nEq) {
            cEq = 0;
          }
        }
      }
    } else {
      // If explicit geometric coupling is used for coupled equations, 
      // only update the equation counter if the current equation is converged
      if (eq.ok) {
        cEq = cEq + 1;
        if (cEq >= com_mod.nEq) {
            cEq = 0;
        }
      }
      // Update the acceleration, velocity, and displacement when the FSI equation is ok
      if (eq.ok && eq.phys == Equation_FSI) {
        if (com_mod.eq[1].phys != Equation_mesh) {
          throw std::runtime_error("PICC: FSI simulation requires a mesh motion equation as the second equation.");
        }
        int s = com_mod.eq[1].s;
        int e = com_mod.eq[1].e;
        for (int Ac = 0; Ac < tnNo; Ac++) {
          if (all_fun::is_domain(com_mod, eq, Ac, Equation_struct) || 
              all_fun::is_domain(com_mod, eq, Ac, Equation_ustruct) || 
              all_fun::is_domain(com_mod, eq, Ac, Equation_lElas)) {
            for (int i = 0; i < e-s+1; i++) {
              An(i+s,Ac) = An(i,Ac);
              Yn(i+s,Ac) = Yn(i,Ac);
              Dn(i+s,Ac) = Dn(i,Ac);
            }
          }
        }
      }
    }

  } else {
    if (eq.ok) {
      cEq = cEq + 1;
    }
  }
 #ifdef debug_picc
 dmsg << "eq " << " coupled ";
 dmsg << "2nd update cEq: " << cEq;
 #endif
}

//---------
// pic_eth
//---------
// Pressure correction at edge nodes for Taylor-Hood type element.
// Here, we interpolate pressure at the edge nodes by interpolating
// using a reduced basis (such as P1) applied on element vertices
// (i.e., corner nodes). For e.g., for a P2 element, pressure is
// interpolated at the edge nodes using P1 vertices.
//
// Modifies: com_mod.Yn
//
void pic_eth(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cep_mod = simulation->get_cep_mod();

  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  const double dt = com_mod.dt;
  const auto& cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];

  auto& cDmn = com_mod.cDmn;
  auto& Yn = com_mod.Yn;

  // Check for something ...
  //
  bool THflag = false;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    if (com_mod.msh[iM].nFs == 2) {
      THflag = true;
      break;
     }
  }

  if (!THflag) {
    return;
  }

  Vector<double> sA(tnNo), sF(tnNo);
  int s = eq.s;

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    auto& msh = com_mod.msh[iM];
    if (msh.nFs == 1) {
      continue;
    }

    auto eType = msh.fs[1].eType;
    int eNoN = msh.fs[0].eNoN;
    int eNoNq = msh.fs[1].eNoN;

    Array<double> xl(nsd,eNoN), xql(nsd,eNoNq), Nqx(nsd,eNoNq);
    Vector<double> pl(eNoNq), Nq(eNoNq); 

    Vector<double> xp(nsd), xi0(nsd), xi(nsd); 
    Array<double> ksix(nsd,nsd);

    for (int g = 0; g < msh.fs[1].nG; g++) {
      for (int i = 0; i < nsd; i++) {
        xi0(i) = xi0(i) + msh.fs[1].xi(i,g);
      }
    }

    xi0 = xi0 / static_cast<double>(msh.fs[1].nG);

    for (int e = 0; e < msh.nEl; e++) {
      cDmn = all_fun::domain(com_mod, msh, cEq, e);       // setting global cDmn
      if ((eq.dmn[cDmn].phys != Equation_stokes) && 
          (eq.dmn[cDmn].phys != Equation_fluid)  && 
          (eq.dmn[cDmn].phys != Equation_ustruct)) {
        continue;
      }

      for (int a = 0; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xl(i,a) = com_mod.x(i,Ac);
        }
      }

      for (int a = 0; a < eNoNq; a++) {
        int Ac = msh.IEN(a,e);
        pl(a) = Yn(s+nsd,Ac);
        for (int i = 0; i < nsd; i++) {
          xql(i,a) = xl(i,a);
        }
      }

      double eVol = 0.0;
      double Jac = 0.0;

      for (int g = 0; g < msh.fs[1].nG; g++) {
        if (g == 0 || !msh.fs[1].lShpF) {
          auto Nx = msh.fs[1].Nx.slice(g);
          nn::gnn(eNoNq, nsd, nsd, Nx, xql, Nqx, Jac, ksix);

          if (utils::is_zero(Jac)) {
            throw std::runtime_error("[pic_eth] Jacobian for element " + std::to_string(e) + " is < 0.");
          }
        }

        eVol = eVol + msh.fs[1].w(g)*Jac;
      }

      for (int a = eNoNq; a < eNoN; a++) {
        int Ac = msh.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          xp(i) = xl(i,a);
        }
        xi = xi0;
        nn::get_nnx(nsd, eType, eNoNq, xql, msh.fs[1].xib, msh.fs[1].Nb, xp, xi, Nq, Nqx);

        double p = 0.0;
        for (int b = 0; b < eNoNq; b++) {
          p = p + pl(b)*Nq(b);
        }

        sF(Ac) = sF(Ac) + p*eVol;
        sA(Ac) = sA(Ac) + eVol;
      }
    } // e-loop
  } // iM-loop

  all_fun::commu(com_mod, sA);
  all_fun::commu(com_mod, sF);

  for (int a = 0; a < tnNo; a++) {
    if (!utils::is_zero(sA(a))) {
      Yn(s+nsd,a) = sF(a) / sA(a);
    }
  }
}

//------
// pici
//------
// This is the initiator.  
//
// Uses Generalized α− Method for time stepping.
//
// Modifes Ag from combination of An and Ao defined by coefs from eq.am, eq.af, 
//   Ag = (1 - eq.am) * Ao  +  eq.am * An
//   Yg = (1 - eq.af) * Yo  +  eq.af * Yn
//   Dg = (1 - eq.af) * Do  +  eq.af * Dn
//
// Modifies:
//   Ag - acceleration
//   Yg - velocity
//   Dg - displacement
//
void pici(Simulation* simulation, Array<double>& Ag, Array<double>& Yg, Array<double>& Dg)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;

  #define n_debug_pici
  #ifdef debug_pici
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  const int cEq = com_mod.cEq;
  const int tnNo = com_mod.tnNo;
  auto& eq = com_mod.eq[cEq];
  auto& dof = com_mod.dof;
  eq.itr = eq.itr + 1;

  // [NOTE] Setting gobal variable 'dof'.
  dof = eq.dof;
  #ifdef debug_pici
  dmsg << "cEq: " << cEq;
  dmsg << "eq.itr: " << eq.itr;
  dmsg << "dof: " << dof;
  dmsg << "tnNo: " << tnNo;
  dmsg << "com_mod.pstEq: " << com_mod.pstEq;
  #endif

  const auto& Ao = com_mod.Ao;
  const auto& An = com_mod.An;
  const auto& Do = com_mod.Do;
  const auto& Dn = com_mod.Dn;
  const auto& Yo = com_mod.Yo;
  const auto& Yn = com_mod.Yn;

  for (int i = 0; i < com_mod.nEq; i++) {
    auto& eq = com_mod.eq[i];
    int s = eq.s;
    int e = eq.e;
    Vector<double> coef(4);
    coef(0) = 1.0 - eq.am;
    coef(1) = eq.am;
    coef(2) = 1.0 - eq.af;
    coef(3) = eq.af;
    #ifdef debug_pici
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    dmsg << "coef: " << coef[0] << " " << coef[1] << " " << coef[2] << " " << coef[3];
    #endif

    if ((eq.phys == Equation_heatF) && (com_mod.usePrecomp)){
        for (int a = 0; a < tnNo; a++) {
            for (int j = 0; j < com_mod.nsd; j++) {
                //Ag(j, a) = An(j, a);
                //Yg(j, a) = Yn(j, a);
                //Dg(j, a) = Dn(j, a);
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);
                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
        for (int a = 0; a < tnNo; a++) {
            for (int j = s; j <= e; j++) {
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);
                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
    } else {
        for (int a = 0; a < tnNo; a++) {
            for (int j = s; j <= e; j++) {
                // eqn 89 of Bazilevs 2007
                Ag(j, a) = Ao(j, a) * coef(0) + An(j, a) * coef(1);
                
                // eqn 90 of Bazilevs 2007
                Yg(j, a) = Yo(j, a) * coef(2) + Yn(j, a) * coef(3);
                
                Dg(j, a) = Do(j, a) * coef(2) + Dn(j, a) * coef(3);
            }
        }
    }
  }
  
  // prestress
  if (com_mod.pstEq) {
    com_mod.pSn = 0.0;
    com_mod.pSa = 0.0;
  }
}

//------
// picp
//------
// This is the predictor.
//
// Modifies:
//   pS0 
//   Ad 
//   Ao 
//   Yo 
//   Do 
//   An 
//   Yn 
//   Dn 
//
void picp(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;

  #define n_debug_picp
  #ifdef debug_picp
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "pstEq: " << com_mod.pstEq;
  #endif

  // Variables for prestress calculations
  auto& pS0 = com_mod.pS0;
  auto& pSn = com_mod.pSn;

  // time derivative of displacement
  auto& Ad = com_mod.Ad;
  
  auto& Ao = com_mod.Ao;
  auto& An = com_mod.An;
  auto& Yo = com_mod.Yo;
  auto& Yn = com_mod.Yn;
  auto& Do = com_mod.Do;
  auto& Dn = com_mod.Dn;

  // Prestress initialization
  if (com_mod.pstEq) {
     pS0 = pS0 + pSn;
     Ao = 0.0;
     Yo = 0.0;
     Do = 0.0;
  }

  // IB treatment: Set dirichlet BC and update traces. For explicit
  // coupling, compute FSI forcing and freeze it for the time step.
  // For implicit coupling, project IB displacement on background
  // mesh and predict quantities at next time step
  //
  // [NOTE] not implemented.
  /*
  if (ibFlag) {
    // Set IB Dirichlet BCs
    CALL IB_SETBCDIR(ib.Yb, ib.Ubo)

    // Update IB location and tracers
    CALL IB_UPDATE(Do)

    if (ib.cpld == ibCpld_E) {
      // FSI forcing for immersed bodies (explicit coupling)
      CALL IB_CALCFFSI(Ao, Yo, Do, ib.Auo, ib.Ubo)

    } else { if (ib.cpld == ibCpld_I) {
      // Project IB displacement (Ubn) to background mesh
      CALL IB_PRJCTU(Do)

    //  Predictor step for implicit coupling
    CALL IB_PICP()
    }
  }
  */

  const auto& dt = com_mod.dt;
  #ifdef debug_picp
  dmsg << "dt: " << dt;
  dmsg << "dFlag: " << com_mod.dFlag;
  #endif

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    auto& eq = com_mod.eq[iEq];
    int s = eq.s; // start row
    int e = eq.e; // end row

    #ifdef debug_picp
    dmsg << "----- iEq " << iEq << " -----";
    dmsg << "s: " << s;
    dmsg << "e: " << e;
    dmsg << "eq.gam: " << eq.gam;
    dmsg << "coef: " << coef;
    #endif

    // [TODO:DaveP] careful here with s amd e.
    double coef = (eq.gam - 1.0) / eq.gam;
    for (int i = s; i <= e; i++) {
      for (int j = 0; j < Ao.ncols(); j++) {
        // eqn 87 of Bazilevs 2007
        An(i,j) = Ao(i,j) * coef;
      }
    }

    // electrophysiology
    if (eq.phys == Equation_CEP) {
      cep_ion::cep_integ(simulation, iEq, e, Do);
    }

    // eqn 86 of Bazilevs 2007
    Yn.set_rows(s,e, Yo.rows(s,e));

    if (com_mod.dFlag) {

      // struct, lElas, FSI (struct, mesh)
      if (!com_mod.sstEq) {
        double coef = dt*dt*(0.5*eq.gam - eq.beta) / (eq.gam - 1.0);
        Dn.set_rows(s,e, Do.rows(s,e) + Yn.rows(s,e)*dt + An.rows(s,e)*coef);

      // ustruct, FSI
      //
      } else {

        if (eq.phys == Equation_ustruct || eq.phys == Equation_FSI) {
          double coef = (eq.gam - 1.0) / eq.gam;
          Ad = Ad*coef;
          Dn.set_rows(s,e, Do.rows(s,e));

        } else if (eq.phys == Equation_mesh) {
          double coef = dt*dt*(0.5*eq.gam - eq.beta) / (eq.gam - 1.0);
          Dn.set_rows(s,e, Do.rows(s,e) + Yn.rows(s,e)*dt + An.rows(s,e)*coef);
        }
      }
    } else {
      Dn.set_rows(s,e, Do.rows(s,e));
    }
  }
}

};

