// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "eq_assem.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "nn.h"
#include "utils.h"

#include "electrophysiology/cep.h"
#include "cmm.h"
#include "porous/darcy.h"
#include "fluid/fluid.h"
#include "fsi.h"
#include "thermal/heatf.h"
#include "thermal/heats.h"
#include "solid/l_elas.h"
#include "mesh.h"
#include "solid/shells.h"
#include "fluid/stokes.h"
#include "solid/sv_struct.h"
#include "solid/ustruct.h"

#include <fsils_api.hpp>

#include <math.h>

namespace eq_assem {

void b_assem_neu_bc(ComMod& com_mod, const faceType& lFa, const Vector<double>& hg, const SolutionStates& solutions)
{
  // Local aliases for solution arrays
  const auto& Do = solutions.old.get_displacement();
  const auto& Yg = solutions.intermediate.get_velocity();

  #define n_debug_b_assem_neu_bc
  #ifdef debug_b_assem_neu_bc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  const int cEq = com_mod.cEq;
  const auto& eq = com_mod.eq[cEq];
  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;

  const int iM = lFa.iM;
  const int eNoN = lFa.eNoN;
  const auto& msh = com_mod.msh[iM];
  auto& cDmn = com_mod.cDmn;

  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);
    cDmn = all_fun::domain(com_mod, msh, cEq, Ec);
    auto cPhys = eq.dmn[cDmn].phys;

    Vector<int> ptr(eNoN); 
    Vector<double> N(eNoN), hl(eNoN); 
    Array<double> yl(tDof,eNoN), lR(dof,eNoN); 
    Array3<double> lK(dof*dof,eNoN,eNoN);

    for (int a = 0; a < eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      ptr(a) = Ac;
      hl(a) = hg(Ac);
      for (int i = 0; i < tDof; i++) {
        yl(i,a) = Yg(i,Ac);
      }
    }

    // Updating the shape functions, if neccessary
    if (lFa.eType == ElementType::NRB) {
      //CALL NRBNNXB(msh(iM), lFa, e)
    }

    for (int g = 0; g < lFa.nG; g++) {
      auto Nx = lFa.Nx.rslice(g);
      Vector<double> nV = nn::gnnb(com_mod, lFa, e, g, Nx, solutions);
      double Jac = utils::norm(nV);
      nV = nV / Jac;
      double w = lFa.w(g)*Jac;
      N  = lFa.N.col(g);

      double h = 0.0;
      Vector<double> y(tDof);

      for (int a = 0; a < eNoN; a++) {
        h = h + N(a)*hl(a);
        y = y + N(a)*yl.col(a);
      }

      switch ( cPhys) {
        case EquationType::phys_fluid:
          fluid::b_fluid(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_CMM:
          fluid::b_fluid(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_heatS:
          heats::b_heats(com_mod, eNoN, w, N, h, lR);
        break;

        case EquationType::phys_heatF:
          heatf::b_heatf(com_mod, eNoN, w, N, y, h, nV, lR, lK);
        break;

        case EquationType::phys_darcy:
          darcy::b_darcy(com_mod, eNoN, w, N, h, lR);
        break;

        case EquationType::phys_lElas:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_struct:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_ustruct:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_shell:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_mesh:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_stokes:
          l_elas::b_l_elas(com_mod, eNoN, w, N, h, nV, lR);
        break;

        case EquationType::phys_CEP:
          cep::b_cep(com_mod, eNoN, w, N, h, lR);
        break;

        default:
          throw std::runtime_error("[b_assem_neu_bc] Undefined physics selection for assembly");
      }
    }

    eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
  }
}

/// @brief  For struct/ustruct - construct follower pressure load contribution
/// to the residual vector and stiffness matrix.
/// We use Nanson's formula to take change in normal direction with
/// deformation into account. Additional calculations based on mesh
/// need to be performed.
///
/// Reproduces 'SUBROUTINE BNEUFOLWP(lFa, hg, Dg)'
/// @param com_mod 
/// @param lBc 
/// @param lFa 
/// @param hg Pressure magnitude
/// @param Dg 
void b_neu_folw_p(ComMod& com_mod, const bcType& lBc, const faceType& lFa, const Vector<double>& hg, const SolutionStates& solutions)
{
  using namespace consts;
  using namespace utils;
  const auto& Dg = solutions.intermediate.get_displacement();

  #define n_debug_b_neu_folw_p
  #ifdef debug_b_neu_folw_p 
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lFa.name: " << lFa.name;
  #endif

  const int cEq = com_mod.cEq;
  const int nsd = com_mod.nsd;
  const int dof = com_mod.dof;
  const int tDof = com_mod.tDof;

  const int iM = lFa.iM;
  const auto& msh = com_mod.msh[iM];
  const int eNoN = msh.eNoN;
  const int eNoNb = lFa.eNoN;

  auto& eq = com_mod.eq[cEq];
  auto& cDmn = com_mod.cDmn;

  #ifdef debug_b_neu_folw_p 
  dmsg << "nsd: " << nsd;
  dmsg << "eNoN: " << nsd;
  #endif

  for (int e = 0; e < lFa.nEl; e++) {
    int Ec = lFa.gE(e);
    cDmn = all_fun::domain(com_mod, msh, cEq, Ec);  // Changes global
    auto cPhys = eq.dmn[cDmn].phys;

    Vector<int> ptr(eNoN); 
    Vector<double> hl(eNoN); 
    Array<double> xl(nsd,eNoN); 
    Array<double> dl(tDof,eNoN);
    Vector<double> N(eNoN); 
    Array<double> Nxi(nsd,eNoN); 
    Array<double> Nx(nsd,eNoN); 
    Array<double> lR(dof,eNoN);
    Array3<double> lK(dof*dof,eNoN,eNoN);
    Array3<double> lKd;

    if (cPhys == EquationType::phys_ustruct) {
      lKd.resize(dof*nsd,eNoN,eNoN);
    }

    // Create local copies
    for (int a = 0; a < eNoN; a++) {
      int Ac = msh.IEN(a,Ec);
      ptr(a) = Ac;
      hl(a) = hg(Ac);

      for (int i = 0; i < nsd; i++) {
        xl(i,a) = com_mod.x(i,Ac);
      }

      for (int i = 0; i < tDof; i++) {
        dl(i,a) = Dg(i,Ac);
      }
    }

    // Initialize parameteric coordinate for Newton's iterations
    Vector<double> xi0(nsd);
    for (int g = 0; g < msh.nG; g++) {
      xi0 = xi0 + msh.xi.col(g);
    }
    xi0 = xi0 / static_cast<double>(msh.nG);

    for (int g = 0; g < lFa.nG; g++) {
      Vector<double> xp(nsd);
      double Jac;

      for (int a = 0; a < eNoNb; a++) {
        int Ac = lFa.IEN(a,e);
        xp = xp + com_mod.x.col(Ac) * lFa.N(a,g);
      }

      auto xi = xi0;
      nn::get_nnx(nsd, msh.eType, eNoN, xl, msh.xib, msh.Nb, xp, xi, N, Nxi);

      if (g == 0 || !msh.lShpF) {
        Array<double> ksix(nsd,nsd);
        nn::gnn(eNoN, nsd, nsd, Nxi, xl, Nx, Jac, ksix);
      }

      // Get surface normal vector
      auto Nx_g = lFa.Nx.rslice(g);
      Vector<double> nV = nn::gnnb(com_mod, lFa, e, g, Nx_g, solutions);
      Jac = utils::norm(nV);
      nV = nV / Jac;
      double w = lFa.w(g)*Jac;

      // Compute residual and tangent contributions
      if (cPhys == EquationType::phys_ustruct) {
        if (nsd == 3) {
          ustruct::b_ustruct_3d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd);
        } else {
          ustruct::b_ustruct_2d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK, lKd);
        }

      } else if (cPhys == EquationType::phys_struct) {
        if (nsd == 3) {
          struct_ns::b_struct_3d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK);
        } else {
          struct_ns::b_struct_2d(com_mod, eNoN, w, N, Nx, dl, hl, nV, lR, lK);
        }
      }
    }

    if (cPhys == EquationType::phys_ustruct) {
      ustruct::ustruct_do_assem(com_mod, eNoN, ptr, lKd, lK, lR);
    } else if (cPhys == EquationType::phys_struct) {
      eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
    }
  }
}

void fsi_ls_upd(ComMod& com_mod, const bcType& lBc, const faceType& lFa, const SolutionStates& solutions)
{
  using namespace consts;
  using namespace utils;
  using namespace fsi_linear_solver;

  #define n_debug_fsi_ls_upd
  #ifdef debug_fsi_ls_upd
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lFa.name: " << lFa.name;
  #endif

  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;
  int nNo = lFa.nNo;

  // Load that a unit pressure on the face puts on each node: the surface
  // integral of the shape function times the normal vector
  Array<double> unit_pressure_load(nsd, tnNo);
  Array<double> unit_pressure_load_on_face_nodes(nsd, nNo);

  const auto &eq = com_mod.eq[com_mod.cEq];
  const auto &Yn = solutions.current.get_velocity();

  // The flow rate of a struct equation coupled to svZeroD or svOneD is
  // integrated on the deformed face, so it also depends on the displacement.
  const bool flowrate_depends_on_displacement =
      (eq.phys == EquationType::phys_struct) &&
      utils::btest(lBc.bType,
                   static_cast<int>(BoundaryConditionType::bType_Coupled));

  // Derivative of the flow rate with respect to the nodal velocities, which
  // is equal to the unit pressure load. Where the flow rate depends on the
  // displacement, its derivative with respect to the nodal displacements is
  // added, with the velocity held fixed. It is multiplied by beta*dt/gamma,
  // the change of displacement per unit change of velocity that an increment
  // of the unknown produces.
  Array<double> flowrate_derivative(nsd, tnNo);
  Array<double> flowrate_derivative_on_face_nodes(nsd, nNo);
  const double displacement_per_velocity = eq.beta * com_mod.dt / eq.gam;

  // Updating the value of the surface integrals using the deformed
  // configuration ('n' = new = timestep n+1)
  unit_pressure_load = 0.0;
  for (int e = 0; e < lFa.nEl; e++) {
    for (int g = 0; g < lFa.nG; g++) {
      auto Nx = lFa.Nx.rslice(g);

      auto cfg = MechanicalConfigurationType::new_timestep;

      const Vector<double> n = nn::gnnb(com_mod, lFa, e, g, Nx, solutions, cfg,
                                        com_mod.eq[com_mod.cEq].s);

      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          unit_pressure_load(i, Ac) += lFa.N(a, g) * lFa.w(g) * n(i);
          flowrate_derivative(i, Ac) += lFa.N(a, g) * lFa.w(g) * n(i);
        }
      }

      if (flowrate_depends_on_displacement) {
        const auto derivative_of_normal =
            nn::derivative_of_outward_normal_wrt_position(com_mod, lFa, e, Nx,
                                                          solutions, cfg, eq.s);

        // Retrieve the velocity at the Gauss point
        Vector<double> velocity(nsd);
        for (int a = 0; a < lFa.eNoN; a++) {
          int Ac = lFa.IEN(a, e);
          for (int i = 0; i < nsd; i++) {
            velocity(i) += lFa.N(a, g) * Yn(eq.s + i, Ac);
          }
        }

        // The flow rate density is velocity . normal, so its derivative with
        // respect to the node positions is velocity . derivative_of_normal
        for (int a = 0; a < lFa.eNoN; a++) {
          int Ac = lFa.IEN(a, e);
          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < nsd; j++) {
              flowrate_derivative(i, Ac) += displacement_per_velocity *
                                            lFa.w(g) * velocity(j) *
                                            derivative_of_normal(j, i, a);
            }
          }
        }
      }
    }
  }

  for (int a = 0; a < lFa.nNo; a++) {
    int Ac = lFa.gN(a);
    unit_pressure_load_on_face_nodes.set_col(a, unit_pressure_load.col(Ac));
    flowrate_derivative_on_face_nodes.set_col(a, flowrate_derivative.col(Ac));
  }

  // Update the unit pressure load and the flow rate derivative of lhs.face(i)
  fsils_bc_update(com_mod.lhs, lBc.lsPtr, unit_pressure_load_on_face_nodes,
                  flowrate_derivative_on_face_nodes);

  // If this is a Coupled BC with a cap, append the cap contribution to the
  // flow rate derivative of the linear solver face
  if (utils::btest(lBc.bType, static_cast<int>(BoundaryConditionType::bType_Coupled))) {
    auto& face = com_mod.lhs.face[lBc.lsPtr];
    auto cfg = MechanicalConfigurationType::new_timestep;

    auto& cpl = lBc.coupled_bc;
    if (cpl.has_cap()) {
      cpl.append_cap_flowrate_derivative(com_mod, face, cfg);
    }
  }
}

/// @brief This routine assembles the equation on a given mesh.
///
/// Ag(tDof,tnNo), Yg(tDof,tnNo), Dg(tDof,tnNo) are taken from solutions.intermediate.
//
void global_eq_assem(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions)
{
  #define n_debug_global_eq_assem
  #ifdef debug_global_eq_assem
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lM.name: " << lM.name;
  com_mod.timer.set_time();
  #endif

  using namespace consts;

  int cEq = com_mod.cEq;
  auto& eq = com_mod.eq[cEq];
  #ifdef debug_global_eq_assem
  dmsg << "cEq: " << cEq;
  dmsg << "eq.sym: " << eq.sym;
  dmsg << "eq.phys: " << eq.phys;
  #endif

  switch (eq.phys) {

    case EquationType::phys_fluid:
      fluid::construct_fluid(com_mod, lM, solutions);
    break;

    case EquationType::phys_heatF:
      heatf::construct_heatf(com_mod, lM, solutions);
    break;

    case EquationType::phys_heatS:
      heats::construct_heats(com_mod, lM, solutions);
    break;

    case EquationType::phys_darcy:
      darcy::construct_darcy(com_mod, lM, solutions);
    break;

    case EquationType::phys_lElas:
      l_elas::construct_l_elas(com_mod, lM, solutions);
    break;

    case EquationType::phys_struct:
      struct_ns::construct_dsolid(com_mod, cep_mod, lM, solutions);
    break;

    case EquationType::phys_ustruct:
      ustruct::construct_usolid(com_mod, cep_mod, lM, solutions);
    break;

    case EquationType::phys_CMM:
      cmm::construct_cmm(com_mod, lM, solutions);
    break;

    case EquationType::phys_shell:
      shells::construct_shell(com_mod, lM, solutions);
    break;

    case EquationType::phys_FSI:
      fsi::construct_fsi(com_mod, cep_mod, lM, solutions);
    break;

    case EquationType::phys_mesh:
      mesh::construct_mesh(com_mod, cep_mod, lM, solutions);
    break;

    case EquationType::phys_CEP:
      cep::construct_cep(com_mod, cep_mod, lM, solutions);
    break;

    case EquationType::phys_stokes:
      stokes::construct_stokes(com_mod, lM, solutions);
    break;

    default:
      throw std::runtime_error("[global_eq_assem] Undefined physics selection for assembly");
  } 

  #ifdef debug_global_eq_assem
  double elapsed_time = com_mod.timer.get_elapsed_time();
  dmsg << "elapsed_time: " << elapsed_time;
  #endif
}

};

