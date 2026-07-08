// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "uris.h"

#include "all_fun.h"
#include "consts.h"
#include "lhsa.h"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"
#include "set_bc.h"

#include "load_msh.h"
#include "vtk_xml.h"
#include "read_msh.h"
#include "VtkData.h"

namespace uris { 

void find_closest_element_centroid(const mshType& mesh, const Vector<double>& xp,
                                   double& minS, int& element_index, Vector<double>& xb);

/// @brief This subroutine computes the mean pressure and flux on the 
/// immersed surface 
void uris_meanp(ComMod& com_mod, CmMod& cm_mod, const int iUris, const SolutionStates& solutions) {
  // Local aliases for solution arrays
  auto& Yn = solutions.current.get_velocity();
  auto& Do = solutions.old.get_displacement();
  #define n_debug_uris_meanp
  #ifdef debug_uris_meanp
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "urisFlag: " << com_mod.urisFlag;
  #endif

  using namespace consts;

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& uris = com_mod.uris;
  auto& msh = com_mod.msh;
  auto& uris_obj = uris[iUris];

  const int nsd = com_mod.nsd;

  // Compute the mean pressure in the upstream and downstream regions of the fluid mesh 

  // Dimensionless factor scaling sdf_deps_close to set the outer signed distance limit 
  // of the upstream/downstream fluid bands used for mean pressure computation.
  double sdf_region_factor = 5.0; 

  // Set the limit of the upstream and downstream regions to 5 times the closed 
  // valve thickness. This should give a reasonable range for the upstream and 
  // downstream regions.
  double meanp_sdf_outer_limit = uris_obj.sdf_deps_close * sdf_region_factor;
  double volU = 0.0;
  double volD = 0.0;

  # ifdef debug_uris_meanp
  if (cm.mas(cm_mod)) {
    dmsg << "Computing upstream region from SDF -" + std::to_string(meanp_sdf_outer_limit) + " to -" 
              + std::to_string(uris_obj.sdf_deps_close) + " for: " + uris_obj.name << std::endl;
    dmsg << "Computing downstream region from SDF " + std::to_string(uris_obj.sdf_deps_close) 
              + " to " + std::to_string(meanp_sdf_outer_limit) + " for: " + uris_obj.name << std::endl;
  }
  # endif

  // Compute the upstream region: negative sdf side (opposite to valve normal), 
  // outside resistance region
  Array<double> sUPS(1, com_mod.tnNo);
  sUPS = 0.0;
  for (int j = 0; j < com_mod.tnNo; j++) {
    double sdf_j = uris_obj.sdf(j);
    if (sdf_j >= -meanp_sdf_outer_limit && sdf_j <= -uris_obj.sdf_deps_close) {
        sUPS(0,j) = 1.0;
    }
  }

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    volU += all_fun::integ(com_mod, cm_mod, iM, sUPS, solutions);
  }

  // Compute the downstream region: positive sdf side (valve normal direction), 
  // outside resistance region
  Array<double> sDST(1, com_mod.tnNo);
  sDST = 0.0;
  for (size_t j = 0; j < com_mod.tnNo; j++) {
    double sdf_j = uris_obj.sdf(j);
    if (sdf_j >= uris_obj.sdf_deps_close && sdf_j <= meanp_sdf_outer_limit) {
      sDST(0,j) = 1.0;
    }
  } 

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    volD += all_fun::integ(com_mod, cm_mod, iM, sDST, solutions);
  }

  // Print volume messages.
  # ifdef debug_uris_meanp
  if (cm.mas(cm_mod)) {
    dmsg << "volume upstream " + std::to_string(volU) + " for: " + uris_obj.name << std::endl;
    dmsg << "volume downstream " + std::to_string(volD) + " for: " + uris_obj.name << std::endl;
  }
  # endif

  double meanPU = 0.0;
  double meanPD = 0.0;

  int iEq = 0;
  int m = 1;
  int s = eq[iEq].s + nsd;
  // int e = s + m - 1;

  Array<double> tmpV(maxNSD, com_mod.tnNo);
  for (int j = 0; j < Yn.ncols(); j++) {
    tmpV(0,j) = Yn(s,j)*sUPS(j);
  }
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    meanPU += all_fun::integ(com_mod, cm_mod, iM, tmpV, solutions);
  }
  meanPU = meanPU / volU;

  for (int j = 0; j < Yn.ncols(); j++) {
    tmpV(0,j) = Yn(s,j)*sDST(j);
  }
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    meanPD += all_fun::integ(com_mod, cm_mod,iM, tmpV, solutions);
  }
  meanPD = meanPD / volD;

  uris_obj.meanPU = uris_obj.relax_factor * meanPU +
                       (1.0 - uris_obj.relax_factor) * uris_obj.meanPU;
  uris_obj.meanPD = uris_obj.relax_factor * meanPD +
                       (1.0 - uris_obj.relax_factor) * uris_obj.meanPD;

  # ifdef debug_uris_meanp
  if (cm.mas(cm_mod)) {
    dmsg << "mean P upstream " + std::to_string(meanPU) + " " 
        + std::to_string(uris_obj.meanPU) + " for: " + uris_obj.name << std::endl;
    dmsg << "mean P downstream " + std::to_string(meanPD) + " " 
        + std::to_string(uris_obj.meanPD) + " for: " + uris_obj.name << std::endl;
  }
  #endif

  //  If the uris has passed the closing state
  if (uris_obj.cnt > uris_obj.DxClose.nslices()) {
    if (uris_obj.meanPU > uris_obj.meanPD) {
      uris_obj.cnt = 1;
      uris_obj.clsFlg = false;
      com_mod.urisActFlag = true;
      # ifdef debug_uris_meanp
      if (cm.mas(cm_mod)) {
        dmsg << "** Set urisCloseFlag to FALSE for: " + uris_obj.name << std::endl;
      }
      # endif
    }
  }
  # ifdef debug_uris_meanp
  if (cm.mas(cm_mod)) {
    dmsg << "urisCloseFlag is: " + std::to_string(uris_obj.clsFlg) + " for: " + uris_obj.name << std::endl;
  }
  # endif
}

/// @brief This subroutine computes the mean velocity in the fluid elements 
/// near the immersed surface  
void uris_meanv(ComMod& com_mod, CmMod& cm_mod, const int iUris, const SolutionStates& solutions) {
  // Local aliases for solution arrays
  auto& Yn = solutions.current.get_velocity();
  auto& Do = solutions.old.get_displacement();
  #define n_debug_uris_meanv
  #ifdef debug_uris_meanv
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "urisFlag: " << com_mod.urisFlag;
  #endif

  using namespace consts;

  auto& cm = com_mod.cm;
  auto& eq = com_mod.eq;
  auto& uris = com_mod.uris;
  auto& msh = com_mod.msh;
  auto& uris_obj = uris[iUris];

  const int nsd = com_mod.nsd;

  // Compute the neighboring region with negative sdf within the 
  // valve's bounding box. When the valve is open, this region 
  // should roughly be valve orifice.
  int iEq = 0;

  double Deps = uris_obj.sdf_deps;
  Array<double> sImm(1, com_mod.tnNo);
  sImm = 0.0;
  double volI = 0.0;

  for (int i = 0; i < com_mod.tnNo; i++) {
    if (uris_obj.sdf(i) <= -Deps) {
      sImm(0,i) = 1.0;
    }
  }

  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    volI += all_fun::integ(com_mod, cm_mod, iM, sImm, solutions);
  }
  # ifdef debug_uris_meanv
  if (cm.mas(cm_mod)) {
    dmsg << "volume inside " + std::to_string(volI) + " for: " + uris_obj.name << std::endl;
  }
  # endif
  
  int m = nsd;
  int s = eq[iEq].s;

  Array<double> tmpV(maxNSD, com_mod.tnNo);
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < Yn.ncols(); j++) {
      tmpV(i,j) = Yn(s+i,j)*sImm(0,j);
    }
  }

  Array<double> tmpVNrm(1, com_mod.tnNo);
  for (int i = 0; i < com_mod.tnNo; i++) {
    tmpVNrm(0,i) = tmpV(0,i) * uris_obj.nrm(0) 
                 + tmpV(1,i) * uris_obj.nrm(1) 
                 + tmpV(2,i) * uris_obj.nrm(2);
  }

  double meanV = 0.0;
  for (int iM = 0; iM < com_mod.nMsh; iM++) {
    meanV += all_fun::integ(com_mod, cm_mod, iM, tmpVNrm, solutions)/volI;
  }
  
  # ifdef debug_uris_meanv
  if (cm.mas(cm_mod)) {
    dmsg << "mean velocity: " + std::to_string(meanV) + " for: " + uris_obj.name << std::endl;
  }
  # endif

  // If the uris has passed the open state
  if (uris_obj.cnt > uris_obj.DxOpen.nslices()) {
    if (meanV < 0.0) {
      uris_obj.cnt = 1;
      uris_obj.clsFlg = true;
      com_mod.urisActFlag = true;
      # ifdef debug_uris_meanv
      if (cm.mas(cm_mod)) {
        dmsg << "** Set urisCloseFlag to TRUE for: " + uris_obj.name << std::endl;
      }
      # endif
    }
  }
  # ifdef debug_uris_meanv
  if (cm.mas(cm_mod)) {
    dmsg << "urisCloseFlag is: " + std::to_string(uris_obj.clsFlg) + " for: " + uris_obj.name << std::endl;
  }
  # endif
}

/// @brief  This subroutine computes the displacement of the immersed 
/// surface with fem projection
void uris_update_disp(ComMod& com_mod, CmMod& cm_mod, const SolutionStates& solutions) {
  #define n_debug_uris_update_disp 
  #ifdef debug_uris_update_disp
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  // Local alias for solution array
  const auto& Dn = solutions.current.get_displacement();

  auto& cm = com_mod.cm;
  auto& uris = com_mod.uris;
  auto& msh = com_mod.msh;
  int nUris = com_mod.nUris;

  const int nsd = com_mod.nsd;

  Array<double> localYd, xl, Nxi;
  Vector<double> N;
  Vector<double> xp(nsd), xi(nsd), d(nsd);
  bool fl;

  for (int iUris = 0; iUris < nUris; iUris++) {
    // Localize each immersed surface node into a fluid element and elect
    // one owner rank per node. Result is cached in elemId and localNode;
    // the search is skipped on subsequent steps if the valve state is unchanged.
    uris_find_tetra(com_mod, cm_mod, iUris);

    auto& uris_obj = uris[iUris];
    localYd.resize(nsd, uris_obj.tnNo);
    localYd = 0.0;
    int local_xi_fail_count = 0;
    int local_elem_miss_count = 0;

    for (int nd = 0; nd < uris_obj.tnNo; nd++) {
      // Skip nodes not owned by this rank to avoid double-counting
      // in the subsequent MPI_SUM gather.
      if (!uris_obj.localNode(nd)) { continue; }

      int jM = uris_obj.elemId(0, nd);
      if (jM == -1) {
        // Defensive check: localNode guarantees a containing element was found
        // across all fluid meshes during uris_find_tetra, so jM == -1 here
        // indicates an inconsistency.
        local_elem_miss_count += 1;
        continue;
      }
      auto& mesh = msh[jM];

      int iEln = uris_obj.elemId(1, nd);
      Vector<double> xp = uris_obj.x.col(nd);
      xl.resize(nsd, mesh.eNoN);
      N.resize(mesh.eNoN);
      Nxi.resize(nsd, mesh.eNoN);

      for (int a = 0; a < mesh.eNoN; a++) {
        int Ac = mesh.IEN(a, iEln);
        xl.rcol(a) = com_mod.x.rcol(Ac);
      }

      // Localize xp inside the parent element to get reference coords xi
      xi = 0.0;
      nn::get_xi(nsd, mesh.eType, mesh.eNoN, xl, xp, xi, fl);
      if (!fl) {
        local_xi_fail_count += 1;
        if (cm.mas(cm_mod)) {
          std::cout << "[WARNING] URIS get_xi not converging!" << std::endl;
        }
        continue;
      }

      // Evaluate shape functions N at xi
      nn::get_gnn(nsd, mesh.eType, mesh.eNoN, xi, N, Nxi);

      // Interpolate displacement at node xp using shape functions.
      d = 0.0;
      for (int a = 0; a < mesh.eNoN; a++) {
        int Ac = mesh.IEN(a, iEln);
        d = d + N(a) * Dn.rows(nsd+1, 2*nsd, Ac);
      }
      // update uris disp  
      localYd.set_col(nd, d);
    }

    // Each node is written by exactly one rank (localNode), so MPI_SUM
    // correctly assembles the global displacement array without double-counting.
    MPI_Allreduce(localYd.data(), uris_obj.Yd.data(), uris_obj.tnNo*nsd,
                  cm_mod::mpreal, MPI_SUM, cm.com());

    int xi_fail_count = 0;
    int elem_miss_count = 0;
    MPI_Allreduce(&local_xi_fail_count, &xi_fail_count, 1,
                  cm_mod::mpint, MPI_SUM, cm.com());
    MPI_Allreduce(&local_elem_miss_count, &elem_miss_count, 1,
                  cm_mod::mpint, MPI_SUM, cm.com());
    if (cm.mas(cm_mod) && (xi_fail_count > 0 || elem_miss_count > 0)) {
      std::cout << "[URIS WARNING] URIS disp update failed for " << uris_obj.name
                << ": elem_miss=" << elem_miss_count
                << ", xi_fail=" << xi_fail_count << std::endl;
    }
  }
}

/// @brief  This subroutine computes the tetrahedral elements
void uris_find_tetra(ComMod& com_mod, CmMod& cm_mod, const int iUris) {
  #define n_debug_uris_find_tetra 
  #ifdef debug_uris_find_tetra
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  auto& uris = com_mod.uris;
  auto& uris_obj = uris[iUris];
  const int nsd = com_mod.nsd;

  // We need to check if the valve needs to move
  int cnt;
  if (!uris_obj.clsFlg) {
    cnt = std::min(uris_obj.cnt, uris_obj.DxOpen.nslices());
  } else {
    cnt = std::min(uris_obj.cnt, uris_obj.DxClose.nslices());
  }

  if (uris_obj.elemId.allocated() && cnt < uris_obj.cnt) {
    return;
  }

  // For each point in the immersed surface we need to localize it 
  // and find the fluid element that contains the node
  // Since the fluid element could be on another processor, we need to
  // gather the displacement values at the end 

  bool include_bdry = true;
  if (!uris_obj.elemId.allocated()) {
    uris_obj.elemId.resize(2, uris_obj.tnNo);
  }
  // localNode(nd) = 1 if this rank is the elected owner of node nd, 0 otherwise.
  // Ownership is resolved here via centroid distance + MPI_MIN so that
  // uris_update_disp does not need a second ownership reduce.
  if (!uris_obj.localNode.allocated()) {
    uris_obj.localNode.resize(uris_obj.tnNo);
  }

  uris_obj.elemId = -1;
  uris_obj.localNode = 0;

  Vector<int> local_counter(uris_obj.tnNo);
  Vector<int> global_counter(uris_obj.tnNo);
  // local_metric(nd): squared distance from node nd to the centroid of its
  // found element on this rank. Used to elect a single owner when multiple
  // ranks find the same node.
  Vector<double> local_metric(uris_obj.tnNo);
  Vector<double> global_metric(uris_obj.tnNo);
  Array<double>  xl;

  local_counter = 0;
  local_metric  = std::numeric_limits<double>::max();

  for (int nd = 0; nd < uris_obj.tnNo; nd++) {
    Vector<double> xp = uris_obj.x.col(nd);
    bool found = false;
    for (int jM = 0; jM < com_mod.nMsh && !found; jM++) {
      auto& mesh = com_mod.msh[jM];
      xl.resize(nsd, mesh.eNoN);
      for (int iEln = 0; iEln < mesh.nEl && !found; iEln++) {
        for (int a = 0; a < mesh.eNoN; a++) {
          int Ac = mesh.IEN(a, iEln);
          xl.rcol(a) = com_mod.x.rcol(Ac);
        }
        if (inside_tet(com_mod, mesh.eNoN, xp, xl, include_bdry)) {
          uris_obj.elemId(0, nd) = jM;
          uris_obj.elemId(1, nd) = iEln;
          local_counter(nd) += 1;

          // Compute squared distance from node to element centroid.
          // This metric is used below to elect one owner per node via
          // MPI_MIN, avoiding double-counting in uris_update_disp.
          Vector<double> x_cent(nsd);
          x_cent = 0.0;
          for (int a = 0; a < mesh.eNoN; a++) {
            int Ac = mesh.IEN(a, iEln);
            x_cent = x_cent + com_mod.x.rcol(Ac);
          }
          x_cent = x_cent / mesh.eNoN;
          Vector<double> dx = xp - x_cent;
          local_metric(nd) = dx * dx;
          found = true;
        }
      }
    }
  }

  // Accumulate how many ranks found each node (used for diagnostics only)
  MPI_Allreduce(local_counter.data(), global_counter.data(), 
                uris_obj.tnNo, cm_mod::mpint, MPI_SUM, cm.com());

  // Select one owner per node: first find the global minimum centroid
  // distance, then break ties deterministically by choosing the smallest
  // rank among all ranks that attained that minimum.
  MPI_Allreduce(local_metric.data(), global_metric.data(),
                uris_obj.tnNo, cm_mod::mpreal, MPI_MIN, cm.com());

  const double owner_tol = 1.0e-14;
  const int my_rank = cm.idcm();
  Vector<int> local_owner_rank(uris_obj.tnNo);
  Vector<int> global_owner_rank(uris_obj.tnNo);

  for (int nd = 0; nd < uris_obj.tnNo; nd++) {
    bool is_min_metric = (local_metric(nd) != std::numeric_limits<double>::max()) &&
                         (std::fabs(local_metric(nd) - global_metric(nd)) <= owner_tol);
    local_owner_rank(nd) = is_min_metric ? my_rank : std::numeric_limits<int>::max();
  }
  MPI_Allreduce(local_owner_rank.data(), global_owner_rank.data(),
                uris_obj.tnNo, cm_mod::mpint, MPI_MIN, cm.com());

  for (int nd = 0; nd < uris_obj.tnNo; nd++) {
    bool is_owner = (local_owner_rank(nd) != std::numeric_limits<int>::max()) &&
                    (global_owner_rank(nd) == my_rank);
    uris_obj.localNode(nd) = is_owner ? 1 : 0;
  }

  // Diagnostic: nodes found on more than one rank indicate ghost/shared
  // element overlap. Ownership is still unique because ties are broken
  // deterministically by rank before the MPI_SUM gather.
  int multi_owner_count = 0;
  for (int nd = 0; nd < uris_obj.tnNo; nd++) {
      if (global_counter(nd) > 1) {
          multi_owner_count += 1;
      }
  }
  if (cm.mas(cm_mod) && multi_owner_count > 0) {
      std::cout << "[URIS WARNING] Multi-owner nodes for " << uris_obj.name
                << ": " << multi_owner_count << std::endl;
  }
}


/// @brief This subroutine check if a node is inside a tetrahedron
bool inside_tet(ComMod& com_mod, int& eNoN, Vector<double>& xp, 
                Array<double>& xl, bool include_bdry) {
  #define n_debug_inside_tet 
  #ifdef debug_inside_tet
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& uris = com_mod.uris;
  const int nsd = com_mod.nsd;
  Vector<double> minb(nsd);
  Vector<double> maxb(nsd);

  // Create a bounding box around of the current solid location 
  // Make the bbox tolerance scale-aware (relative to element size)
  for (int i = 0; i < nsd; i++) {
    double min_val = std::numeric_limits<double>::max();
    double max_val = std::numeric_limits<double>::lowest();
    for (int j = 0; j < eNoN; j++) {
      double val = xl(i,j);
      if (val < min_val) {min_val = val;}
      if (val > max_val) {max_val = val;}
    }
    double range = std::max(max_val - min_val, 0.0);
    // Relative tolerance with a tiny floor for near-degenerate spans
    double tol = std::max(range * 1.0e-3, 1.0e-12);
    minb(i) = min_val - tol;
    maxb(i) = max_val + tol;
  }

  // Is the node inside the BBox? 
  for (int i = 0; i < nsd; ++i) {
    if (xp(i) < minb(i) || xp(i) > maxb(i)) {
      return false;
    }
  }
  return in_poly(xp, xl, include_bdry);
}

/// @brief Precompute whether each node belongs to a fluid-related domain.
void uris_build_fluid_node_mask(ComMod& com_mod) {
  using namespace consts;

  com_mod.urisFluidNodeMask.resize(com_mod.tnNo, 0);
  for (int a = 0; a < com_mod.tnNo; ++a) {
    for (int iEq = 0; iEq < com_mod.nEq; ++iEq) {
      const auto& eq = com_mod.eq[iEq];
      if (all_fun::is_domain(com_mod, eq, a, Equation_fluid) ||
          all_fun::is_domain(com_mod, eq, a, Equation_CMM) ||
          all_fun::is_domain(com_mod, eq, a, Equation_stokes)) {
        com_mod.urisFluidNodeMask[a] = 1;
        break;
      }
    }
  }
}

/// @brief Build an expanded bounding box around coordinates.
void uris_compute_expanded_bbox(const Array<double>& x, const int nsd, const double expansion,
                                Vector<double>& minb, Vector<double>& maxb) {
  minb.resize(nsd);
  maxb.resize(nsd);
  Vector<double> min_val(nsd);
  Vector<double> max_val(nsd);

  for (int i = 0; i < nsd; i++) {
    min_val(i) = std::numeric_limits<double>::max();
    max_val(i) = std::numeric_limits<double>::lowest();
    for (int j = 0; j < x.ncols(); j++) {
      const double val = x(i,j);
      if (val < min_val(i)) { min_val(i) = val; }
      if (val > max_val(i)) { max_val(i) = val; }
    }
  }
  // Compute the diagonal length of the bounding box for use in degenerate cases where max_val == min_val
  double diag_length = utils::norm((max_val - min_val));

  for (int i = 0; i < nsd; i++) {
    double extra = 0.0;
    if (max_val(i) > min_val(i)) {
      extra = (max_val(i) - min_val(i)) * expansion;
    } else if (max_val(i) == min_val(i)) {
      // When points are all the same in this dimension, use the diagonal length of the bounding box
      extra = expansion * (diag_length > 0.0 ? diag_length : 1.0);
    } else {
      throw std::runtime_error("Invalid bounding box: max_val < min_val for dimension " + std::to_string(i));
    }
    minb(i) = min_val(i) - extra;
    maxb(i) = max_val(i) + extra;
  }
}

/// @brief Check if a point lies inside a bounding box.
bool uris_point_in_bbox(const Vector<double>& xp, const Vector<double>& minb,
                        const Vector<double>& maxb, const int nsd) {
  for (int i = 0; i < nsd; ++i) {
    if (xp(i) < minb(i) || xp(i) > maxb(i)) {
      return false;
    }
  }
  return true;
}

/// @brief Load a shell mesh from a VTU file and initialize element metadata.
void load_shell_mesh_from_file(Simulation* simulation, mshType& mesh,
                               const std::string& mesh_path,
                               const bool check_element_ordering) {
  mesh.lShl = true;
  vtk_xml::read_vtu(mesh_path, mesh);

  const int nsd = simulation->com_mod.nsd;
  if (nsd == 1) {
    throw std::runtime_error("The number of spatial dimensions (" + std::to_string(nsd) +
        ") is not consistent with the mesh '" + mesh.name + "' which contains shell elements.");
  }

  nn::select_ele(simulation->com_mod, mesh);

  if (check_element_ordering) {
    read_msh_ns::check_ien(simulation, mesh);
  }
}

/// @brief Load and initialize the optional URIS scaffold mesh.
///
/// @param uris_name Name used only to identify the URIS object in error messages.
mshType load_scaffold_from_file(Simulation* simulation,
                                const std::string& scaffold_file_path,
                                const std::string& uris_name,
                                const double scaffold_scale) {
  mshType scaffold_mesh;

  try {
    load_shell_mesh_from_file(simulation, scaffold_mesh, scaffold_file_path, true);
  } catch (const std::exception& e) {
    throw std::runtime_error("Failed to read URIS scaffold mesh for '" + uris_name + "': " + e.what());
  } catch (...) {
    throw std::runtime_error("Failed to read URIS scaffold mesh for '" + uris_name + "'.");
  }

  // Scale the scaffold mesh coordinates by scF to match the URIS mesh scale.
  for (int a = 0; a < scaffold_mesh.gnNo; a++) {
    scaffold_mesh.x.rcol(a) = scaffold_mesh.x.rcol(a) * scaffold_scale;
  }

  int b = 0;
  scaffold_mesh.nNo = scaffold_mesh.gnNo;
  scaffold_mesh.gN.resize(scaffold_mesh.nNo);
  scaffold_mesh.gN = 0;
  scaffold_mesh.lN.resize(scaffold_mesh.nNo);
  scaffold_mesh.lN = 0;
  for (int a = 0; a < scaffold_mesh.nNo; a++) {
    scaffold_mesh.gN(a) = b;
    scaffold_mesh.lN(b) = a;
    b++;
  }

  scaffold_mesh.nEl = scaffold_mesh.gnEl;
  scaffold_mesh.IEN.resize(scaffold_mesh.eNoN, scaffold_mesh.nEl);
  for (int e = 0; e < scaffold_mesh.nEl; e++) {
    for (int a = 0; a < scaffold_mesh.eNoN; a++) {
      int Ac = scaffold_mesh.gIEN(a,e);
      Ac = scaffold_mesh.gN(Ac);
      scaffold_mesh.IEN(a,e) = Ac;
    }
  }
  scaffold_mesh.gIEN.clear();

  return scaffold_mesh;
}

/// @brief Read the URIS mesh separately 
void uris_read_msh(Simulation* simulation) {
  #define n_debug_uris_read_msh 
  #ifdef debug_uris_read_msh
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& com_mod = simulation->get_com_mod();
  auto& uris = com_mod.uris;
  const int nsd = com_mod.nsd;

  com_mod.urisFlag = true;
  com_mod.urisActFlag = true;

  int nUris = simulation->parameters.URIS_mesh_parameters.size();
  com_mod.nUris = nUris;
  # ifdef debug_uris_read_msh
  dmsg << "Number of immersed surfaces for uris: " + std::to_string(nUris) << std::endl;
  # endif
  uris.resize(nUris);

  for (int iUris = 0; iUris < nUris; iUris++) {
    auto param = simulation->parameters.URIS_mesh_parameters[iUris];
    auto& uris_obj = uris[iUris];
    uris_obj.name = param->name();
    # ifdef debug_uris_read_msh
    dmsg << "** Reading URIS mesh: " + uris_obj.name << std::endl;
    # endif

    uris_obj.scF = param->mesh_scale_factor();
    uris_obj.nFa = param->URIS_face_parameters.size();
    uris_obj.msh.resize(uris_obj.nFa);
    uris_obj.nrm.resize(nsd);
    Array<double> gX(0,0);

    std::string positive_flow_normal_file_path = param->positive_flow_normal_file_path();
    // lPtr => lPM%get(fTmp, "Positive flow normal file")
    // fid = fTmp%open()
    // READ (fid,*) uris(iUris)%nrm(:)
    // CLOSE (fid)
    std::ifstream file_stream;
    file_stream.open(positive_flow_normal_file_path);
    if (!file_stream.is_open()) {
      throw std::runtime_error("Failed to open the open positive flow normal file '" 
                               + positive_flow_normal_file_path + "'.");
    }
    for (int i = 0; i < nsd; i++) {
      file_stream >> uris_obj.nrm(i);
      if (!file_stream) {
        throw std::runtime_error(
            "Failed to read positive flow normal (component=" + std::to_string(i)
            + ") from positive flow normal file '"
            + positive_flow_normal_file_path + "'.");
      }
    }
    file_stream.close();

    uris_obj.sdf_deps = param->thickness();
    uris_obj.sdf_deps_close = param->close_thickness();
    // Use large default value for the signed distance function to indicate that 
    // the fluid node is far away from the valve. 
    uris_obj.sdf_default = param->close_thickness() * 1e6;
    uris_obj.resistance = param->resistance();
    uris_obj.clsFlg = param->valve_starts_as_closed();
    uris_obj.invert_normal = param->invert_normal();
    uris_obj.include_uris_velocity = param->include_uris_velocity();
    std::string scaffold_file_path = param->scaffold_file_path();

    if (scaffold_file_path != "") {
      uris_obj.scaffold_flag = true;
      uris_obj.scaffold_msh = load_scaffold_from_file(simulation, scaffold_file_path,
          uris_obj.name, uris_obj.scF);
      
      # ifdef debug_uris_read_msh
      dmsg << "Scaffold mesh is included for: " + uris_obj.name << std::endl;
      dmsg << "Scaffold mesh nodes: " + std::to_string(uris_obj.scaffold_msh.gnNo) << std::endl;
      dmsg << "Scaffold mesh elements: " + std::to_string(uris_obj.scaffold_msh.gnEl) << std::endl;
      # endif
    }

    // uris_obj.tnNo = 0;
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      // Set as shell
      auto mesh_param = param->URIS_face_parameters[iM];
      auto& mesh = uris_obj.msh[iM];
      mesh.lShl = true;
      mesh.name = mesh_param->name();
      # ifdef debug_uris_read_msh
      dmsg << "-- Reading URIS face: " + mesh.name << std::endl;
      # endif

      // Read mesh nodal coordinates and element connectivity.
      uris_read_sv(simulation, mesh, mesh_param);

      // // [HZ] eType is not NA, neglecting this for now
      // IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
      // CALL READCCNE(lPN, uris(iUris)%msh(iM))
      // END IF
      // IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
      //     CALL READNRB(lPN, uris(iUris)%msh(iM))
      // END IF
      // IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
      //     CALL READGAMBIT(lPN, uris(iUris)%msh(iM))
      // END IF
      // IF (uris(iUris)%msh(iM)%eType .EQ. eType_NA) THEN
      //     err = " Failed to identify format of the uris mesh"
      // END IF

      # ifdef debug_uris_read_msh
      dmsg << "Number of uris nodes: " + std::to_string(mesh.gnNo) << std::endl;
      dmsg << "Number of uris elements: " + std::to_string(mesh.gnEl) << std::endl;
      # endif

      // Read valve motion: note that this motion is defined on the 
      // reference configuration 
      std::string open_motion_file_path = mesh_param->open_motion_file_path();
      // std::ifstream file_stream;
      file_stream.open(open_motion_file_path);
      if (!file_stream.is_open()) {
        throw std::runtime_error("Failed to open the open motion file '" 
                                 + open_motion_file_path + "'.");
      }
      int dispNtOpen, dispNnOpen;
      file_stream >> dispNtOpen >> dispNnOpen;
      if (!file_stream) {
        throw std::runtime_error(
            "Failed to read time steps and node numbers from open motion file '"
            + open_motion_file_path + "'.");
      }
      // std::cout << "dispNtOpen: " << dispNtOpen << std::endl;
      // std::cout << "dispNnOpen: " << dispNnOpen << std::endl;

      if (dispNnOpen != mesh.gnNo) {
        throw std::runtime_error("Mismatch in node numbers between URIS mesh and displacements.");
      }

      Array3<double> dispOpen(nsd, dispNnOpen, dispNtOpen);
      for (int k = 0; k < dispNtOpen; k++) {
        for (int j = 0; j < dispNnOpen; j++) {
          for (int i = 0; i < nsd; i++) {
            file_stream >> dispOpen(i,j,k);
            if (!file_stream) {
              throw std::runtime_error(
                  "Failed to read displacement (time=" + std::to_string(k)
                  + ", node=" + std::to_string(j)
                  + ", component=" + std::to_string(i)
                  + ") from open motion file '"
                  + open_motion_file_path + "'.");
            }
          }
        }
      }
      file_stream.close();

      std::string close_motion_file_path = mesh_param->close_motion_file_path();
      // std::ifstream file_stream;
      file_stream.open(close_motion_file_path);
      if (!file_stream.is_open()) {
        throw std::runtime_error("Failed to open the close motion file '" 
                                 + close_motion_file_path + "'.");
      }
      int dispNtClose, dispNnClose;
      file_stream >> dispNtClose >> dispNnClose;
      if (!file_stream) {
        throw std::runtime_error(
            "Failed to read time steps and node numbers from close motion file '"
            + close_motion_file_path + "'.");
      }
      // std::cout << "dispNtClose: " << dispNtClose << std::endl;
      // std::cout << "dispNnClose: " << dispNnClose << std::endl;

      if (dispNnClose != mesh.gnNo) {
        throw std::runtime_error("Mismatch in node numbers between URIS mesh and displacements.");
      }

      Array3<double> dispClose(nsd, dispNnClose, dispNtClose);
      for (int k = 0; k < dispNtClose; k++) {
        for (int j = 0; j < dispNnClose; j++) {
          for (int i = 0; i < nsd; i++) {
            file_stream >> dispClose(i,j,k);
            if (!file_stream) {
              throw std::runtime_error(
                  "Failed to read displacement (time=" + std::to_string(k)
                  + ", node=" + std::to_string(j)
                  + ", component=" + std::to_string(i)
                  + ") from open motion file '"
                  + close_motion_file_path + "'.");
            }
          }
        }
      }
      file_stream.close();

      if (iM > 0) {
        if (dispNtOpen != uris_obj.DxOpen.nslices()) {
          throw std::runtime_error(
              "Mismatch in open motion time steps for URIS mesh '" + uris_obj.name +
              "', face '" + mesh.name + "'. Expected " +
              std::to_string(uris_obj.DxOpen.nslices()) + ", got " +
              std::to_string(dispNtOpen) + ".");
        }
        if (dispNtClose != uris_obj.DxClose.nslices()) {
          throw std::runtime_error(
              "Mismatch in close motion time steps for URIS mesh '" + uris_obj.name +
              "', face '" + mesh.name + "'. Expected " +
              std::to_string(uris_obj.DxClose.nslices()) + ", got " +
              std::to_string(dispNtClose) + ".");
        }
      }

      // To scale the mesh, while attaching x to gX
      int a = uris_obj.tnNo + mesh.gnNo;

      if (iM == 0) {
        gX.resize(nsd, a);
        gX = 0.0;
        uris_obj.DxOpen.resize(nsd, a, dispNtOpen);
        uris_obj.DxClose.resize(nsd, a, dispNtClose);
      } else{
        Array<double> tmpX(nsd, uris_obj.tnNo);
        tmpX = gX;
        gX.resize(nsd, a);
        for (int j = 0; j < uris_obj.tnNo; j++) {
          gX.rcol(j) = tmpX.rcol(j);
        }

        // Move data for open
        auto tmpDxOpen = uris_obj.DxOpen;
        uris_obj.DxOpen.resize(nsd, a, dispNtOpen);
        for (int i = 0; i < dispNtOpen; i++) {
          for (int j = 0; j < uris_obj.tnNo; j++) {
            for (int k = 0; k < nsd; k++) {
              uris_obj.DxOpen(k,j,i) = tmpDxOpen(k,j,i);
            }
          }
        }
        // Move data for close
        auto tmpDxClose = uris_obj.DxClose;
        uris_obj.DxClose.resize(nsd, a, dispNtClose);
        for (int i = 0; i < dispNtClose; i++) {
          for (int j = 0; j < uris_obj.tnNo; j++) {
            for (int k = 0; k < nsd; k++) {
              uris_obj.DxClose(k,j,i) = tmpDxClose(k,j,i);
            }
          }
        }
      }

      for (int i = uris_obj.tnNo; i < a; i++) {
        gX.rcol(i) = mesh.x.rcol(i-uris_obj.tnNo) * uris_obj.scF;
      }

      for (int k = 0; k < dispNtOpen; k++) {
        for (int j = uris_obj.tnNo; j < a; j++) {
          for (int i = 0; i < nsd; i++) {
            uris_obj.DxOpen(i,j,k) = dispOpen(i,j-uris_obj.tnNo,k) * uris_obj.scF;
          }
        }
      }

      for (int k = 0; k < dispNtClose; k++) {
        for (int j = uris_obj.tnNo; j < a; j++) {
          for (int i = 0; i < nsd; i++) {
            uris_obj.DxClose(i,j,k) = dispClose(i,j-uris_obj.tnNo,k) * uris_obj.scF;
          }
        }
      }

      uris_obj.tnNo = a;
    }

    uris_obj.x.resize(nsd, uris_obj.tnNo);
    // Set the valve position in the initial position
    if (uris_obj.clsFlg) {
      int dispNtClose = uris_obj.DxClose.nslices();
      uris_obj.x = uris_obj.DxClose.rslice(dispNtClose-1);
    } else {
      int dispNtOpen = uris_obj.DxOpen.nslices();
      uris_obj.x = uris_obj.DxOpen.rslice(dispNtOpen-1);
    }
    uris_obj.Yd.resize(nsd, uris_obj.tnNo);
    uris_obj.Yd = 0.0;

    if (uris_obj.include_uris_velocity) {
      uris_obj.x_prev.resize(nsd, uris_obj.tnNo);
      uris_obj.x_prev = uris_obj.x;
      uris_obj.valve_velocity.resize(nsd, uris_obj.tnNo);
      uris_obj.valve_velocity = 0.0;
    }

    // Setting mesh.gN, mesh.lN parameter
    int b = 0;
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      auto& mesh = uris_obj.msh[iM];
      mesh.nNo = mesh.gnNo;
      mesh.gN.resize(mesh.nNo);
      mesh.gN = 0;
      mesh.lN.resize(uris_obj.tnNo);
      mesh.lN = 0;
      for (int a = 0; a < mesh.nNo; a++) {
        mesh.gN(a) = b;
        mesh.lN(b) = a;
        b++;
      }
    }
    if (b != uris_obj.tnNo) {
      throw std::runtime_error("Mismatch in uris.tnNo. Correction needed.");
    }

    // Remap mesh.gIEN to mesh.IEN
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      auto& mesh = uris_obj.msh[iM];
      mesh.nEl = mesh.gnEl;
      mesh.IEN.resize(mesh.eNoN, mesh.nEl);
      for (int e = 0; e < mesh.nEl; e++) {
        for (int a = 0; a < mesh.eNoN; a++) {
          int Ac = mesh.gIEN(a,e);
          Ac = mesh.gN(Ac);
          mesh.IEN(a,e) = Ac;
        }
      }
      mesh.gIEN.clear();
    }

    if (uris_obj.nFa > 0) {
      # ifdef debug_uris_read_msh
      dmsg << "Total number of uris nodes: " + std::to_string(uris_obj.tnNo) << std::endl;
      # endif
      int total_nel = 0;
      for (int iM = 0; iM < uris_obj.nFa; iM++) {
          total_nel += uris_obj.msh[iM].nEl;
      }
      # ifdef debug_uris_read_msh
      dmsg << "Total number of uris elements: " + std::to_string(total_nel) << std::endl;
      # endif
    }
  }
  # ifdef debug_uris_read_msh
  dmsg << "URIS mesh data imported successfully." << std::endl;
  # endif
}

/// @brief Write URIS solution to a vtu file
void uris_write_vtus(ComMod& com_mod) {
  #define n_debug_uris_write_vtus 
  #ifdef debug_uris_write_vtus
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  auto& uris = com_mod.uris;
  const int nsd = com_mod.nsd;
  const int nUris = com_mod.nUris;

  // we plot coord + displacement
  int nOut = 2;
  int outDof = nOut * nsd;
  std::vector<std::string> outNames(nOut);
  Vector<int> outS(nOut+1);

  // Prepare all solultions in to dataType d
  std::vector<dataType> d;
  for (int iUris = 0; iUris < nUris; iUris++) {
    auto& uris_obj = uris[iUris];
    d.resize(uris_obj.nFa);
    int nNo = 0;
    int nEl = 0;
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      auto& mesh = uris_obj.msh[iM];
      int cOut = 0;
      outS(cOut) = 0;
      outS(cOut+1) = nsd;
      // outNames[cOut] = "";

      // outS = [0, 3]

      if (mesh.eType == ElementType::NRB) {
        throw std::runtime_error("Outputs for NURBS data is under development.");
      }
      d[iM].nNo = mesh.nNo;
      d[iM].nEl = mesh.nEl;
      d[iM].eNoN = mesh.eNoN;
      d[iM].vtkType = mesh.vtkType;
      d[iM].x.resize(outDof, mesh.nNo);
      d[iM].IEN.resize(mesh.eNoN, mesh.nEl);

      for (int a = 0; a < mesh.nNo; a++) {
        int Ac = mesh.gN(a);
        for (int i = 0; i < nsd; i++) {
          // Scale the valve displacement to the unit of the mesh
          d[iM].x(i,a) = uris_obj.x(i,Ac) / uris_obj.scF;
        }
      }
      for (int e = 0; e < mesh.nEl; e++) {
        for (int i = 0; i < mesh.eNoN; i++) {
          d[iM].IEN(i,e) = mesh.IEN(i,e);
        }
      }

      int l = nsd; // l = 3
      int s = 0;
      int e = s + l; // e = 3

      cOut += 1; // cOut = 1
      int is = outS(cOut); // is = outS[1] = 3
      int ie = is + l; // ie = 3 + 3 = 6
      outS(cOut+1) = ie; // outS[2] = 7, outS = [0, 3, 6]
      outNames[0] = "coordinates";
      outNames[1] = "URIS_displacement";

      for (int a = 0; a < mesh.nNo; a++) {
        int Ac = mesh.gN(a);
        for (int i = 0; i < nsd; i++) {
          d[iM].x(is+i,a) = uris_obj.Yd(s+i,Ac) / uris_obj.scF;
        }
      }

      nNo += mesh.nNo;
      nEl += mesh.nEl;
    }

    // Writing to vtu file (master only)
    char fName_num[100];
    if (com_mod.cTS >= 1000) {
      sprintf(fName_num, "%d", com_mod.cTS);
    } else {
      sprintf(fName_num, "%03d", com_mod.cTS);
    }
    std::string fName = com_mod.saveName + "_uris_" + uris_obj.name + "_" + fName_num + ".vtu";

    auto vtk_writer = VtkData::create_writer(fName);
    // Writing the position data
    int iOut = 0;
    int s = outS(iOut);
    // int e = outS(iOut+1) - 1;
    int nSh = 0;
    Array<double> tmpV(maxNSD, nNo);
    tmpV = 0.0;
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      for (int a = 0; a < d[iM].nNo; a++) {
        for (int i = 0; i < nsd; i++) {
          tmpV(i,a+nSh) = d[iM].x(s+i,a);
        }
      }
      nSh += d[iM].nNo;
    }
    vtk_writer->set_points(tmpV);

    // Writing the connectivity data
    Array<int> tmpI;
    for (int iM = 0; iM < uris_obj.nFa; iM++) {
      tmpI.resize(d[iM].eNoN, d[iM].nEl);
      for (int e = 0; e < d[iM].nEl; e++) {
        for (int i = 0; i < d[iM].eNoN; i++) {
          tmpI(i,e) = d[iM].IEN(i,e);
        }
      }
      vtk_writer->set_connectivity(nsd, tmpI);
    }

    // Writing all solutions
    for (int iOut = 0; iOut < nOut; iOut++) {
      int s = outS(iOut);
      int e = outS(iOut+1);
      int l = e - s;
      tmpV.resize(l, nNo);
      tmpV = 0.0;
      int nSh = 0;
      for (int iM = 0; iM < uris_obj.nFa; iM++) {
        for (int a = 0; a < d[iM].nNo; a++) {
          for (int i = 0; i < nsd; i++) {
            tmpV(i,a+nSh) = d[iM].x(s+i,a);
          }
        }
        nSh += d[iM].nNo;
      }
      vtk_writer->set_point_data(outNames[iOut], tmpV);
    }

    vtk_writer->write();

    delete vtk_writer;

  }

}

/// @brief  Checks if a probe lies inside or outside an immersed boundary
void uris_calc_sdf(ComMod& com_mod) {
  #define n_debug_uris_calc_sdf 
  #ifdef debug_uris_calc_sdf
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  using namespace consts;

  auto& cm = com_mod.cm;
  auto& uris = com_mod.uris;
  const int nsd = com_mod.nsd;
  const int nUris = com_mod.nUris;

  if (com_mod.urisFluidNodeMask.size() != static_cast<size_t>(com_mod.tnNo)) {
    uris_build_fluid_node_mask(com_mod);
  }

  for (int iUris = 0; iUris < nUris; iUris++) {
    // First check if the valve needs to move 
    auto& uris_obj = uris[iUris];

    const Array3<double>& Dx = uris_obj.clsFlg ? uris_obj.DxClose : uris_obj.DxOpen;
    int cnt = std::min(uris_obj.cnt, Dx.nslices());
    uris_obj.x = Dx.rslice(cnt - 1);

    if (uris_obj.include_uris_velocity) {
      if (cnt < uris_obj.cnt) {
        // Frozen at last frame: valve not changing
        uris_obj.x_prev = uris_obj.x;
      } else if (cnt > 1) {
        // Mid-motion: previous prescribed frame
        uris_obj.x_prev = Dx.rslice(cnt - 2);
      } else {
        // cnt == 1: first frame of sequence
        uris_obj.x_prev = uris_obj.x;
      }
      for (int i = 0; i < uris_obj.tnNo; i++) {
        uris_obj.valve_velocity.rcol(i) = (uris_obj.x.rcol(i) - uris_obj.x_prev.rcol(i)) / com_mod.dt;
      }
      uris_obj.valve_velocity_fluid = 0.0;
    }

    const bool compute_valve_sdf = !(cnt < uris_obj.cnt && uris_obj.sdf_computed);
    const bool compute_scaffold_udf = uris_obj.scaffold_flag && !uris_obj.scaffold_udf_computed;
    if (!compute_valve_sdf && !compute_scaffold_udf) {
      continue;
    }

    const double bbox_expansion = 0.1;
    Vector<double> minb(nsd);
    Vector<double> maxb(nsd);
    if (compute_valve_sdf) {
      #ifdef debug_uris_calc_sdf
      dmsg << "Recomputing SDF for " << uris_obj.name;
      #endif
      uris_obj.sdf = uris_obj.sdf_default;
      // The valve BBox is 10% larger than the current valve coordinates.
      uris_compute_expanded_bbox(uris_obj.x, nsd, bbox_expansion, minb, maxb);
    }

    auto& scaffold_mesh = uris_obj.scaffold_msh;
    Vector<double> minb_scaf(nsd);
    Vector<double> maxb_scaf(nsd);
    if (compute_scaffold_udf) {
      uris_obj.scaffold_udf = uris_obj.sdf_default;
      uris_compute_expanded_bbox(scaffold_mesh.x, nsd, bbox_expansion, minb_scaf, maxb_scaf);
    }

    // The SDF is computed on the reference configuration, which
    // means that the valves will be morphed based on the fluid mesh
    // motion. If the fluid mesh stretches near the valve, the valve
    // leaflets will also be streched. Note that
    // this is a simplifying assumption. 
    Vector<double> xp(nsd);
    for (int ca = 0; ca < com_mod.tnNo; ca++) {
      xp = com_mod.x.rcol(ca);
      if (!com_mod.urisFluidNodeMask[ca]) {
        continue;
      }

      if (compute_valve_sdf && uris_point_in_bbox(xp, minb, maxb, nsd)) {
        double minS = std::numeric_limits<double>::max();
        int Ec = -1;
        int jM = -1;
        Vector<double> xb(nsd);
        Vector<double> unitNormal(nsd);
        uris_find_closest_face_centroid(uris_obj, xp, nsd, minS, Ec, jM, xb);
        uris_face_unit_normal(uris_obj, nsd, jM, Ec, unitNormal);
        const double dotp = (xp - xb) * unitNormal;
        const double sdf_sign = uris_compute_sdf_sign(uris_obj, xp, xb, dotp);
        uris_obj.sdf[ca] = sdf_sign * minS;

        if (uris_obj.include_uris_velocity) {
          Vector<double> interp_valve_vel(nsd);
          uris_interp_valve_velocity(uris_obj, xp, nsd, jM, Ec, dotp, unitNormal, interp_valve_vel);
          uris_obj.valve_velocity_fluid.rcol(ca) = interp_valve_vel;
        }
      }

      if (compute_scaffold_udf && uris_point_in_bbox(xp, minb_scaf, maxb_scaf, nsd)) {
        double minS_scaf = std::numeric_limits<double>::max();
        int Ec = -1;
        Vector<double> xb(nsd);
        find_closest_element_centroid(scaffold_mesh, xp, minS_scaf, Ec, xb);
        uris_obj.scaffold_udf[ca] = minS_scaf;
      }
    } // ca: loop

    if (compute_valve_sdf) {
      uris_obj.sdf_computed = true;
    }
    if (compute_scaffold_udf) {
      uris_obj.scaffold_udf_computed = true;
    }
  } // iUris: loop

}


/// @brief Create data for a mesh.
///
/// Replicates Fortran READSV subroutine defined in LOADMSH.f.
///
///   SUBROUTINE READSV(list, lM)
//

void uris_read_sv(Simulation* simulation, mshType& mesh, const URISFaceParameters* mesh_param) {
  auto mesh_path = mesh_param->face_file_path();
  auto mesh_name = mesh_param->name();

  #define n_dbg_read_sv
  #ifdef dbg_read_sv
    DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
    dmsg.banner();
    dmsg << "Mesh name: " << mesh_name;
    dmsg << "Mesh path: " << mesh_path;
    dmsg << "mesh.lShl: " << mesh.lShl;
  #endif

  load_shell_mesh_from_file(simulation, mesh, mesh_path, simulation->get_com_mod().ichckIEN);
}


/// @brief Check whether point P lies inside a convex polygon (2D) or
/// tetrahedron (3D) defined by vertices P1. Returns true if P is inside,
/// or if P lies on the boundary and include_bdry is true; otherwise
/// returns false.
bool in_poly(const Vector<double>& P, const Array<double>& P1, bool include_bdry) {
  #define n_dbg_in_poly
  #ifdef dbg_in_poly
    DebugMsg dmsg(__func__, 0);
    dmsg.banner();
    dmsg << "checking in_poly";
  #endif

  const int nd = P1.nrows();

  if (nd == 2) {
      Vector<double> N(nd);
      for (int i = 0; i <= nd; i++) {
          const int j = (i < nd) ? i + 1 : 0;  // wrap last edge back to 0
          N(0) = P1(1,i) - P1(1,j);
          N(1) = P1(0,j) - P1(0,i);
          const double dotP = N(0)*(P(0)-P1(0,i)) + N(1)*(P(1)-P1(1,i));
          if (dotP < 0.0) { return false; }
      }
      return true;
  } else if (nd == 3) {
    // 3D: all four faces must pass the same-side test
    const Vector<double> v1 = P1.col(0);
    const Vector<double> v2 = P1.col(1);
    const Vector<double> v3 = P1.col(2);
    const Vector<double> v4 = P1.col(3);

    return same_side(v1, v2, v3, v4, P, include_bdry)
        && same_side(v2, v3, v4, v1, P, include_bdry)
        && same_side(v3, v4, v1, v2, P, include_bdry)
        && same_side(v4, v1, v2, v3, P, include_bdry);
  } else {
    throw std::runtime_error("Invalid number of dimensions for in_poly");
  }
}

/// @brief Check if a point is on the same side of another point wrt a triangle in 3D
bool same_side(const Vector<double>& v1, const Vector<double>& v2,
               const Vector<double>& v3, const Vector<double>& v4,
               const Vector<double>& p,  bool include_bdry) {
  #define n_dbg_same_side
  #ifdef dbg_same_side
    DebugMsg dmsg(__func__, 0);
    dmsg.banner();
    dmsg << "checking same_side";
  #endif

  Vector<double> v21 = v2 - v1;
  Vector<double> v31 = v3 - v1;
  Array<double> V(3, 2);
  V.set_col(0, v21);
  V.set_col(1, v31);

  const Vector<double> N = utils::cross(V);
  const double dotV4 = N * (v4 - v1);
  const double dotP = N * (p - v1);

  // P and v4 are on the same side if their dot products share a sign
  const bool sameside = (dotP >= 0.0) == (dotV4 >= 0.0);

  // If not strictly same side, accept if P lies on the face within tolerance
  if (include_bdry && !sameside) {
    const double tol = 1.0e-9 * std::max(std::fabs(dotV4), 1.0e-30);
    return std::fabs(dotP) <= tol;
  }

  return sameside;
}

/// @brief Barycenter of URIS surface/shell element (jM, Ec) in current valve coordinates uris_obj.x.
void surface_element_barycenter(const urisType& uris_obj, int jM, int Ec, Vector<double>& xb) {
  const auto& mesh = uris_obj.msh[jM];
  xb = 0.0;
  for (int a = 0; a < mesh.eNoN; a++) {
    const int Ac = mesh.IEN(a, Ec);
    xb = xb + uris_obj.x.rcol(Ac);
  }
  xb = xb / mesh.eNoN;
}

/// @brief Barycenter of a fixed shell element using mesh coordinates.
Vector<double> mesh_element_barycenter(const mshType& mesh, const int element_index) {
  Vector<double> xb(mesh.x.nrows());
  xb = 0.0;
  for (int a = 0; a < mesh.eNoN; a++) {
    const int Ac = mesh.IEN(a, element_index);
    xb = xb + mesh.x.rcol(Ac);
  }
  xb = xb / mesh.eNoN;
  return xb;
}

/// @brief Find the closest URIS shell element centroid to a point.
///
/// @param[in] uris_obj URIS object containing the shell meshes and current
/// valve coordinates.
/// @param[in] xp Background mesh point used to search for the nearest
/// URIS shell element centroid.
/// @param[in] nsd Number of spatial dimensions.
/// @param[in,out] minS Current minimum centroid distance; updated when a
/// closer centroid is found.
/// @param[out] Ec Element index of the closest shell element centroid.
/// @param[out] jM Mesh index containing the closest shell element centroid.
/// @param[out] xb Coordinates of the closest shell element centroid.
void uris_find_closest_face_centroid(const urisType& uris_obj, const Vector<double>& xp,
                                     const int nsd, double& minS, int& Ec, int& jM,
                                     Vector<double>& xb) {
  #define n_dbg_uris_find_closest_face_centroid
  #ifdef dbg_uris_find_closest_face_centroid
  DebugMsg dmsg(__func__, 0);
  dmsg.banner();
  dmsg << "finding closest face centroid";
  #endif

  Vector<double> face_centroid(nsd);
  for (int iM = 0; iM < uris_obj.nFa; iM++) {
    const auto& mesh = uris_obj.msh[iM];
    for (int e = 0; e < mesh.nEl; e++) {
      surface_element_barycenter(uris_obj, iM, e, face_centroid);
      const double dS = utils::norm((xp - face_centroid));
      if (dS < minS) {
        minS = dS;
        Ec = e;
        jM = iM;
        xb = face_centroid;
      }
    }
  }
}

/// @brief Find the closest fixed mesh element centroid to a point.
///
/// @param[in] mesh Mesh containing the elements to search.
/// @param[in] xp Background mesh point used to search for the nearest
/// fixed mesh element centroid.
/// @param[in,out] minS Current minimum centroid distance; updated when a
/// closer centroid is found.
/// @param[out] element_index Index of the element with the closest centroid.
/// @param[out] xb Coordinates of the closest element centroid.
void find_closest_element_centroid(const mshType& mesh, const Vector<double>& xp,
                                   double& minS, int& element_index, Vector<double>& xb) {
  for (int e = 0; e < mesh.nEl; e++) {
    const Vector<double> elem_centroid = mesh_element_barycenter(mesh, e);
    const double dS = utils::norm((xp - elem_centroid));
    if (dS < minS) {
      minS = dS;
      element_index = e;
      xb = elem_centroid;
    }
  }
}

/// @brief Unit normal of URIS face element (jM, Ec) from parametric tangents 
/// (optionally flipped by uris_obj.invert_normal).
void uris_face_unit_normal(const urisType& uris_obj, const int nsd, const int jM, const int Ec,
                                 Vector<double>& unitNormal) {
  #define n_dbg_uris_face_unit_normal
  #ifdef dbg_uris_face_unit_normal
  DebugMsg dmsg(__func__, 0);
  dmsg.banner();
  dmsg << "computing URIS face unit normal";
  #endif

  const auto& mesh = uris_obj.msh[jM];
  Array<double> xXi(nsd, nsd - 1);
  Array<double> lX(nsd, mesh.eNoN);

  xXi = 0.0;
  lX = 0.0;

  for (int a = 0; a < mesh.eNoN; a++) {
    const int Ac = mesh.IEN(a, Ec);
    lX.rcol(a) = uris_obj.x.rcol(Ac);
  }

  for (int a = 0; a < mesh.eNoN; a++) {
    for (int i = 0; i < nsd - 1; i++) {
      xXi.rcol(i) = xXi.rcol(i) + lX.rcol(a) * mesh.Nx(i, a, 0);
    }
  }

  unitNormal = utils::cross(xXi);
  const auto Jac = utils::norm(unitNormal);
  if (uris_obj.invert_normal) {
    unitNormal = -unitNormal / Jac;
  } else {
    unitNormal = unitNormal / Jac;
  }
}

/// @brief Compute SDF sign for open/closed URIS states.
double uris_compute_sdf_sign(const urisType& uris_obj, const Vector<double>& xp,
  const Vector<double>& xb, const double dotP) {
  #define n_dbg_uris_compute_sdf_sign
  #ifdef dbg_uris_compute_sdf_sign
    DebugMsg dmsg(__func__, 0);
    dmsg.banner();
    dmsg << "computing SDF sign";
  #endif

  // Improved implementation for SDF sign. For closed state, sign is less 
  // sensitive to local face normal orientation inconsistency and more 
  // aligned with the intended physical flow direction
  if (uris_obj.clsFlg) {
    auto dot_nrm = (xp - xb) * uris_obj.nrm;
    return (dot_nrm < 0.0 && dotP < 0.0) ? -1.0 : 1.0;
  }
  return (dotP < 0.0) ? -1.0 : 1.0;
}


/// @brief Interpolate valve velocity to fluid node at the given point xp
void uris_interp_valve_velocity(const urisType& uris_obj, const Vector<double>& xp, const int nsd, 
                                const int jM, const int Ec, const double dotP, 
                                const Vector<double>& unitNormal, Vector<double>& interp_valve_vel) {
  #define n_dbg_uris_interp_valve_velocity
  #ifdef dbg_uris_interp_valve_velocity
    DebugMsg dmsg(__func__, 0);
    dmsg.banner();
    dmsg << "interpolating valve velocity";
  #endif

  const auto& mesh = uris_obj.msh[jM];                          
  Vector<double> xp_plane(nsd), E1(nsd), E2(nsd), v(nsd);
  Array<double> lX(nsd, mesh.eNoN);

  for (int a = 0; a < mesh.eNoN; a++) {
    const int Ac = mesh.IEN(a, Ec);
    lX.rcol(a) = uris_obj.x.rcol(Ac);
  }
  
  // project xp onto the triangle plane
  xp_plane = xp - dotP * unitNormal;
  
  // compute barycentric coordinates (xi, eta)
  E1 = lX.rcol(1) - lX.rcol(0);
  E2 = lX.rcol(2) - lX.rcol(0);
  v = xp_plane - lX.rcol(0);

  auto g11 = E1 * E1;
  auto g12 = E1 * E2; 
  auto g22 = E2 * E2;
  auto b1 = v * E1;
  auto b2 = v * E2;
  double det = g11 * g22 - g12 * g12;

  double xi = (g22 * b1 - g12 * b2) / det;
  double eta = (g11 * b2 - g12 * b1) / det;

  // shape functions:
  double N1 = 1.0 - xi - eta;
  double N2 = xi;
  double N3 = eta;

  // interpolate the valve velocity:
  interp_valve_vel = N1 * uris_obj.valve_velocity.rcol(mesh.IEN(0,Ec)) +
                     N2 * uris_obj.valve_velocity.rcol(mesh.IEN(1,Ec)) +
                     N3 * uris_obj.valve_velocity.rcol(mesh.IEN(2,Ec));
  
}

/// @brief Evaluate total Brinkman factor and weighted valve velocity at element quadrature points.
void eval_uris_ris_factors_quadrature(const ComMod& com_mod, const mshType& lM, const fsType& fs, 
  const int e, Vector<double>& uris_factor_total_el, Array<double>& uris_valve_vel_term_total_el) {
  #define n_dbg_eval_uris_ris_factors_quadrature
  #ifdef dbg_eval_uris_ris_factors_quadrature
    DebugMsg dmsg(__func__, 0);
    dmsg.banner();
    dmsg << "computing RIS factor";
  #endif

  const int nUris = com_mod.nUris;
  uris_factor_total_el.resize(fs.nG);
  uris_factor_total_el = 0.0;
  uris_valve_vel_term_total_el.resize(com_mod.nsd, fs.nG);
  uris_valve_vel_term_total_el = 0.0;

  if (!com_mod.urisActFlag) {
    return;
  }

  Vector<double> dist_srf(nUris);
  Vector<double> dist_scaffold(nUris);
  Array<double> valve_velocity(com_mod.nsd, nUris);

  for (int g = 0; g < fs.nG; g++) {
    dist_srf = 0.0;
    dist_scaffold = 0.0;
    valve_velocity = 0.0;
    for (int a = 0; a < fs.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      for (int iUris = 0; iUris < nUris; iUris++) {
        dist_srf(iUris) += fs.N(a,g) * std::fabs(com_mod.uris[iUris].sdf(Ac));
        if (com_mod.uris[iUris].scaffold_flag) {
          dist_scaffold(iUris) += fs.N(a,g) * std::fabs(com_mod.uris[iUris].scaffold_udf(Ac));
        }
        if (com_mod.uris[iUris].include_uris_velocity) {
          valve_velocity.rcol(iUris) = valve_velocity.rcol(iUris) +
              fs.N(a,g) * com_mod.uris[iUris].valve_velocity_fluid.rcol(Ac);
        }
      }
    }

    double sdf_deps;
    double delta_eps;
    double delta_eps_scaffold;
    for (int iUris = 0; iUris < nUris; iUris++) {
      sdf_deps = 0.0;
      delta_eps = 0.0;
      delta_eps_scaffold = 0.0;
      double start_deps, end_deps;
      int n_steps;
      if (com_mod.uris[iUris].clsFlg) {
        start_deps = com_mod.uris[iUris].sdf_deps;
        end_deps   = com_mod.uris[iUris].sdf_deps_close;
        n_steps  = com_mod.uris[iUris].DxClose.nslices();
      } else {
        start_deps = com_mod.uris[iUris].sdf_deps_close;
        end_deps   = com_mod.uris[iUris].sdf_deps;
        n_steps  = com_mod.uris[iUris].DxOpen.nslices();
      }

      if (n_steps <= 0) {
        sdf_deps = end_deps;
      } else if (com_mod.uris[iUris].cnt >= n_steps) {
        sdf_deps = end_deps;
      } else if (com_mod.uris[iUris].cnt <= 0) {
        sdf_deps = start_deps;
      } else {
        // Linear ramping: start_deps -> end_deps over n_steps.
        // This alleviates the sudden change in the resistance factor when the 
        // valve status changes from open to close and reduces oscillations
        double progress = static_cast<double>(com_mod.uris[iUris].cnt) / static_cast<double>(n_steps);
        sdf_deps = start_deps + progress * (end_deps - start_deps);
      }
      if (dist_srf(iUris) < sdf_deps && sdf_deps > 0.0) {
        delta_eps = (1 + cos(consts::pi*dist_srf(iUris)/sdf_deps))/(2*sdf_deps*sdf_deps);
      }
      if (com_mod.uris[iUris].scaffold_flag) {
        // Compute the scaffold resistance factor based on the unsigned distance function (UDF)
        // The scaffold surface uses the same thickness parameter as the closed valve surface
        const double scaffold_deps = com_mod.uris[iUris].sdf_deps_close;
        if (dist_scaffold(iUris) < scaffold_deps && scaffold_deps > 0.0) {
          delta_eps_scaffold = (1 + cos(consts::pi*dist_scaffold(iUris)/scaffold_deps))/(2*scaffold_deps*scaffold_deps);
        }
      }

      uris_factor_total_el(g) += com_mod.uris[iUris].resistance * (delta_eps + delta_eps_scaffold);

      if (com_mod.uris[iUris].include_uris_velocity) {
        uris_valve_vel_term_total_el.rcol(g) = uris_valve_vel_term_total_el.rcol(g) 
          + com_mod.uris[iUris].resistance * delta_eps*valve_velocity.rcol(iUris);
      }
    } // iUris: loop
  } // g: loop
}

}
