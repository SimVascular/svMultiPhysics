// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fsils_api.hpp"
#include "fils_struct.hpp"

#include <vector>

namespace fsi_linear_solver {

namespace {

/// @brief Convert the mesh indices of nodes to linear solver indices.
///
/// @param[in] lhs The left-hand side of the linear system.
/// @param[in] gNodes The mesh indices of the nodes.
/// @param[in] num_nodes The number of nodes.
///
/// @return The linear solver indices of the nodes.
Vector<int> linear_solver_node_indices(const FSILS_lhsType &lhs,
                                       const Vector<int> &gNodes,
                                       int num_nodes) {
  Vector<int> node_indices(num_nodes);

  for (int a = 0; a < num_nodes; a++) {
    node_indices(a) = lhs.map(gNodes(a));
  }

  return node_indices;
}

/// @brief Tell whether more than one process has nodes of a face.
///
/// Every process has to make this call.
///
/// @param[in] lhs The left-hand side of the linear system.
/// @param[in] nNo The number of nodes of the face on this process.
///
/// @return True if more than one process has nodes of the face.
bool is_shared_between_processes(const FSILS_lhsType &lhs, int nNo) {
  if (lhs.commu.nTasks == 1) {
    return false;
  }

  // 1 if this process has nodes of the face, 0 otherwise. Its sum over the
  // processes is the number of processes that have nodes of the face.
  const int has_nodes = nNo != 0 ? 1 : 0;
  int num_processes_with_nodes;

  MPI_Allreduce(&has_nodes, &num_processes_with_nodes, 1, cm_mod::mpint,
                MPI_SUM, lhs.commu.comm);

  return num_processes_with_nodes > 1;
}

/// @brief Sum the values at each node of a face over the processes that have
/// the node.
///
/// Does nothing if face.sharedFlag is false. Otherwise every process has to
/// make this call.
///
/// @param[in] lhs The left-hand side of the linear system.
/// @param[in] face The face, whose sharedFlag, dof and glob are read.
/// @param[in,out] values The values at the nodes of the face (face.dof x the
///   number of nodes, nodes ordered as in face.glob), replaced by their sums.
void sum_nodal_values_over_processes(const FSILS_lhsType &lhs,
                                     const FSILS_faceType &face,
                                     Array<double> &values) {
  if (!face.sharedFlag) {
    return;
  }

  Array<double> v(face.dof, lhs.nNo);

  for (int a = 0; a < values.ncols(); a++) {
    int Ac = face.glob(a);
    for (int i = 0; i < face.dof; i++) {
      v(i, Ac) = values(i, a);
    }
  }

  fsils_commuv(lhs, face.dof, v);

  for (int a = 0; a < values.ncols(); a++) {
    int Ac = face.glob(a);
    for (int i = 0; i < face.dof; i++) {
      values(i, a) = v(i, Ac);
    }
  }
}

/// @brief Append to the flow rate derivative of a Neumann face the values at
/// the nodes that this process owns.
///
/// @param[in] lhs The left-hand side that the face belongs to.
/// @param[in,out] face The face, whose flowrate_derivative_node_indices,
///   flowrate_derivative and scaled_flowrate_derivative are extended.
/// @param[in] node_indices The linear solver index of the node of each column
///   of 'values'.
/// @param[in] values The values at those nodes (face.dof x the number of
///   nodes).
void append_owned_flowrate_derivative(const FSILS_lhsType &lhs,
                                      FSILS_faceType &face,
                                      const Vector<int> &node_indices,
                                      const Array<double> &values) {
  // Columns of 'values' whose node this process owns. Those are the nodes
  // that the linear solver numbers first.
  std::vector<int> owned_columns;

  for (int a = 0; a < values.ncols(); a++) {
    if (node_indices(a) < lhs.mynNo) {
      owned_columns.push_back(a);
    }
  }

  if (owned_columns.empty()) {
    return;
  }

  const int num_old = face.flowrate_derivative.ncols();
  const int num_new = owned_columns.size();

  Vector<int> new_node_indices(num_old + num_new);
  Array<double> new_values(face.dof, num_old + num_new);

  for (int b = 0; b < num_old; b++) {
    new_node_indices(b) = face.flowrate_derivative_node_indices(b);
    new_values.set_col(b, face.flowrate_derivative.col(b));
  }

  for (int b = 0; b < num_new; b++) {
    const int a = owned_columns[b];
    new_node_indices(num_old + b) = node_indices(a);
    new_values.set_col(num_old + b, values.col(a));
  }

  face.flowrate_derivative_node_indices = new_node_indices;
  face.flowrate_derivative = new_values;
  face.scaled_flowrate_derivative.resize(face.dof, num_old + num_new);
}

/// @brief Set the flow rate derivative of a Neumann face from values given at
/// the nodes of the face.
///
/// The flow rate derivative keeps the nodes of the face that this process
/// owns.
///
/// @param[in] lhs The left-hand side that the face belongs to.
/// @param[in,out] face The face, whose flowrate_derivative_node_indices,
///   flowrate_derivative and scaled_flowrate_derivative are replaced.
/// @param[in] values The values at the nodes of the face (face.dof x face.nNo,
///   nodes ordered as in face.glob).
void set_flowrate_derivative_on_face_nodes(const FSILS_lhsType &lhs,
                                           FSILS_faceType &face,
                                           const Array<double> &values) {
  face.flowrate_derivative_node_indices.clear();
  face.flowrate_derivative.clear();
  face.scaled_flowrate_derivative.clear();

  append_owned_flowrate_derivative(lhs, face, face.glob, values);
}

} // namespace

void fsils_bc_create(FSILS_lhsType &lhs, int faIn, BcType BC_type,
                     const Vector<int> &gNodes, const Array<double> &Val) {
  using namespace consts;

  // Number of components per node
  const int dof = Val.nrows();

  // Number of nodes of the face on this process
  const int nNo = Val.ncols();

#define n_debug_fsils_bc_create
#ifdef debug_fsils_bc_create
  DebugMsg dmsg(__func__, lhs.commu.task);
  dmsg.banner();
  dmsg << "lhs.nNo: " << lhs.nNo;
  dmsg << "lhs.nFaces: " << lhs.nFaces;
  dmsg << "faIn: " << faIn;
  dmsg << "nNo: " << nNo;
  dmsg << "dof: " << dof;
  dmsg << "BC_type: " << BC_type;
  dmsg << "Val.size(): " << Val.size();
  dmsg << "Val.nrows: " << Val.nrows_;
  dmsg << "Val.ncols: " << Val.ncols_;
#endif

  if (faIn >= lhs.nFaces) {
    throw std::runtime_error(
        "FSILS: faIn is exceeding lhs structure maximum number of faces (" +
        std::to_string(lhs.nFaces) + ") is less than " + std::to_string(faIn) +
        ".");
  }

  if (faIn <= -1) {
    throw std::runtime_error("FSILS: faIn is smaller than zero");
  }

  lhs.face[faIn].nNo = nNo;
  lhs.face[faIn].dof = dof;
  lhs.face[faIn].bGrp = BC_type;
  lhs.face[faIn].glob = linear_solver_node_indices(lhs, gNodes, nNo);

  lhs.face[faIn].sharedFlag = is_shared_between_processes(lhs, nNo);

  // Store the values and synchronize them on nodes on boundary between
  // processes. The resize sets the size on a process that has no node of the
  // face, where the assignment does nothing.
  if (BC_type == BcType::BC_TYPE_Dir) {
    lhs.face[faIn].free_component_mask.resize(dof, nNo);
    lhs.face[faIn].free_component_mask = Val;

    sum_nodal_values_over_processes(lhs, lhs.face[faIn],
                                    lhs.face[faIn].free_component_mask);
  } else {
    lhs.face[faIn].unit_pressure_load.resize(dof, nNo);
    lhs.face[faIn].unit_pressure_load = Val;

    sum_nodal_values_over_processes(lhs, lhs.face[faIn],
                                    lhs.face[faIn].unit_pressure_load);

    lhs.face[faIn].scaled_unit_pressure_load.resize(dof, nNo);

    set_flowrate_derivative_on_face_nodes(lhs, lhs.face[faIn],
                                          lhs.face[faIn].unit_pressure_load);
  }
}

void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes)
{
  Array<double> Val(dof, nNo);
  fsils_bc_create(lhs, faIn, BC_type, gNodes, Val);
}

/// @brief Frees FSILS boundary condition structure
/// @param lhs 
/// @param faIn 
void fsils_bc_free(FSILS_lhsType& lhs, int faIn)
{
  //IF (.NOT.lhs%face(faIn)%foC) THEN
  //       PRINT *, 'FSILS: Cannot free a face that is not created yet'
  //      STOP "FSILS: FATAL ERROR"
  //END IF

  lhs.face[faIn].foC        = false;
  lhs.face[faIn].nNo        = 0;
  lhs.face[faIn].bGrp       = BcType::BC_TYPE_Dir;
  lhs.face[faIn].res        = 0.0;
  lhs.face[faIn].sharedFlag = false;

  //DEALLOCATE(lhs.face(faIn).glob, lhs.face(faIn).val, lhs.face(faIn).valM)

}

void fsils_bc_update(FSILS_lhsType &lhs, int faIn,
                     const Array<double> &unit_pressure_load,
                     const Array<double> &flowrate_derivative) {
  lhs.face[faIn].unit_pressure_load = unit_pressure_load;
  sum_nodal_values_over_processes(lhs, lhs.face[faIn],
                                  lhs.face[faIn].unit_pressure_load);

  Array<double> summed_flowrate_derivative = flowrate_derivative;
  sum_nodal_values_over_processes(lhs, lhs.face[faIn],
                                  summed_flowrate_derivative);
  set_flowrate_derivative_on_face_nodes(lhs, lhs.face[faIn],
                                        summed_flowrate_derivative);
}

void fsils_bc_append_flowrate_derivative(const FSILS_lhsType &lhs,
                                         FSILS_faceType &face,
                                         const Vector<int> &gNodes,
                                         const Array<double> &values) {

  append_owned_flowrate_derivative(
      lhs, face, linear_solver_node_indices(lhs, gNodes, values.ncols()),
      values);
}
};