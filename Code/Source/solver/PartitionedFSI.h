// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PARTITIONED_FSI_H
#define PARTITIONED_FSI_H

#include "Simulation.h"
#include "Integrator.h"
#include "Array.h"
#include "CmMod.h"

#include <fstream>
#include <memory>
#include <vector>
#include <string>

/// @brief Coupling method for interface relaxation
enum class CouplingMethod { constant, aitken };

/// @brief Configuration for partitioned FSI coupling, read from XML input.
struct PartitionedFSIConfig {
  int max_coupling_iterations = 50;
  double coupling_tolerance = 1e-6;
  double initial_relaxation = 1.0;
  double omega_max = 1.0;
  CouplingMethod coupling_method = CouplingMethod::aitken;

  // Face names for the FSI interface
  std::string fluid_interface_face;
  std::string solid_interface_face;
};

/// @brief Partitioned FSI coupling with 3 independent sub-Simulations.
///
/// Each sub-field (fluid, struct, mesh) has its own Simulation object with
/// independent mesh, solution arrays, and linear system. No shared global
/// arrays, no DOF offsets (each eq.s=0), no regularization of inactive nodes.
///
/// Implements Dirichlet-Neumann coupling with Aitken relaxation:
///   1. Transfer solid displacement to mesh interface, solve mesh equation
///   2. Deform fluid mesh using mesh displacement, solve fluid equation
///   3. Extract fluid traction, apply to solid, solve solid equation
///   4. Extract solid displacement, apply Aitken relaxation
///   5. Check coupling convergence
///
/// Related to GitHub issue #431: Implement partitioned FSI in svMultiPhysics
class PartitionedFSI {
public:
  PartitionedFSI(Simulation* main_simulation, const PartitionedFSIConfig& config,
                 const std::string& xml_file_path);

  ~PartitionedFSI();

  void run();
  bool step();

private:
  Simulation* main_sim_;
  PartitionedFSIConfig config_;
  std::string xml_file_path_;

  // Sub-simulations (owned)
  std::unique_ptr<Simulation> fluid_sim_;
  std::unique_ptr<Simulation> solid_sim_;
  std::unique_ptr<Simulation> mesh_sim_;

  // Interface face pointers within each sub-sim
  const faceType* fluid_face_ = nullptr;
  const faceType* solid_face_ = nullptr;
  const faceType* mesh_face_ = nullptr;

  // Mesh pointers within each sub-sim
  const mshType* fluid_mesh_ = nullptr;
  const mshType* solid_mesh_ = nullptr;
  const mshType* mesh_mesh_ = nullptr;

  // Global face node counts and per-rank offsets (for MPI distribution)
  int solid_face_global_nNo_ = 0;
  int solid_face_local_offset_ = 0;
  int fluid_face_global_nNo_ = 0;
  int fluid_face_local_offset_ = 0;
  int mesh_face_global_nNo_ = 0;
  int mesh_face_local_offset_ = 0;

  // Node maps between interface faces: global_src_idx → global_tgt_idx
  std::vector<int> solid_to_fluid_map_;
  std::vector<int> fluid_to_solid_map_;
  std::vector<int> solid_to_mesh_map_;

  // Canonical slot maps: canonical_[i] = first global slot with same physical node as i
  std::vector<int> solid_face_canonical_;
  std::vector<int> fluid_face_canonical_;
  std::vector<int> mesh_face_canonical_;

  // Per-local-node owner flag for the solid interface face: 1 if THIS rank is
  // the min-rank owner of the node, else 0.  Used to apply the interface
  // traction exactly once per physical node (it is summed by all_fun::commu).
  std::vector<int> solid_face_owner_;

  // Coupling state — indexed by GLOBAL solid face nodes (replicated on all ranks)
  Array<double> disp_prev_;
  Array<double> vel_prev_;
  double omega_;
  double first_res_norm_ = 0.0;

  // Aitken state — sized for global solid face
  std::vector<double> r_prev_;

  // Output files for coupling convergence history
  std::ofstream coupling_log_;
  std::ofstream histor_log_;

  void resolve_faces();
  void build_node_maps();

  /// Solve fluid equation with current interface velocity and ALE mesh velocity
  bool solve_fluid(const Array<double>& mesh_vel_Yo, const Array<double>& mesh_vel_Yn);

  /// Extract fluid traction, transfer to solid, solve solid equation
  bool solve_solid();

  /// Solve mesh equation with relaxed displacement, deform fluid mesh
  bool solve_mesh(const Array<double>& x_ref, int mesh_s);

  /// Compute vel_prev_ (global) from disp_prev_ (global) using Newmark relationship
  void compute_interface_velocity();

  void relax_interface(int cp, int nsd, const Array<double>& disp_current);
  void relax_constant(int cp, int nsd, const Array<double>& disp_current);
  void relax_aitken(int cp, int nsd, const Array<double>& disp_current);

  /// Compute global_nNo and local_offset for a distributed face
  static void compute_face_global_info(const faceType& face, cmType& cm,
                                       const CmMod& cm_mod,
                                       int& global_nNo, int& local_offset);

  /// All-gather local face data (nrows, local_nNo) to global (nrows, global_nNo)
  static Array<double> gather_face_data(const Array<double>& local_data,
                                        int global_nNo, int local_offset,
                                        cmType& cm, const CmMod& cm_mod);

  /// All-gather local (local_src → global_tgt) map to global (global_src → global_tgt) map
  static void gather_global_map(const std::vector<int>& local_map,
                                int global_src_nNo,
                                cmType& cm, const CmMod& cm_mod,
                                std::vector<int>& global_map);

  /// Build local face_a → global face_b node map using pre-gathered global face_b coords
  static void build_face_node_map(const faceType& face_a, const ComMod& com_a,
                                  int global_b_nNo, const Array<double>& global_b_coords,
                                  std::vector<int>& a_to_global_b);

  /// Build a minimal sub-simulation XML for the given role by extracting
  /// the tagged equation and its meshes from the main XML. Returns temp file path.
  static std::string build_sub_xml(const std::string& main_xml_path,
                                   const std::string& role);

  /// Transfer data from global src face to global tgt face using global map
  static Array<double> transfer_data(const std::vector<int>& src_to_tgt_map,
                                     const Array<double>& src_data, int tgt_nNo);

  /// Build canonical[i] = index of first occurrence of global node ID i across ranks
  std::vector<int> build_face_dedup_map(const faceType& face, const ComMod& com,
                                        int global_nNo, cmType& cm, const CmMod& cm_mod);

  /// Propagate canonical slot values to duplicate slots (canonical[i] != i)
  static void apply_dedup(Array<double>& arr, const std::vector<int>& canonical);

  /// Build per-local-node flags marking whether THIS rank is the min-rank
  /// owner of each face node (1 = owner, 0 = shared but owned by a lower rank).
  static std::vector<int> build_face_owner_flags(const faceType& face,
                                                 const ComMod& com, cmType& cm);

  void save_results();
};

#endif // PARTITIONED_FSI_H
