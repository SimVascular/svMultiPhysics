// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef PRESTRAIN_H
#define PRESTRAIN_H

#include "ActiveStress.h"
#include "Array.h"
#include "Array3.h"
#include "mat_fun.h"

#include <cmath>

class CepMod;
class CmMod;
class ComMod;
class dmnType;
class MeshParameters;
class cmType;
class mshType;

/// @brief Prestrain of a body imaged in a loaded configuration.
///
/// The mesh is the imaged, loaded geometry. The prestrain is a deformation
/// gradient F0 at every Gauss point that maps an unknown stress-free
/// configuration to the mesh. A displacement u from the mesh then gives the deformation
/// gradient from the stress-free configuration
///
///     F = (I + Grad u) F0,
///
/// with Grad taken on the mesh. The solid kernels integrate on the stress-free
/// configuration: where the assembly loops compute the shape function
/// gradients, they pull them back in place, F0^T Grad N, with the volume
/// dV / det F0, and pull the fibers back as material lines, F0^-1 f normalized,
/// each on its own, so that a basis orthonormal on the mesh is orthonormal in
/// the imaged configuration, not in the unloaded one. The kernels then form
/// F = F0 + Grad_0 u and keep their total Lagrangian form.
/// The formulation is objective: rotating the deformed body
/// rigidly rotates F and the internal forces with it.
///
/// F0 is found by Newton's method on the imaged geometry. Every Newton
/// iteration starts from rest under the full load, so its residual is that of
/// the current F0. The displacement u it reaches is composed into the
/// prestrain, F0 <- (I + Grad u) F0, and discarded. At the fixed point the mesh
/// carries the load in Cauchy stress equilibrium without moving. This is the
/// incremental scheme of Gee, Forster and Wall, Int. J. Numer. Meth. Biomed.
/// Engng. 26 (2010) 52-72, section 3, with the increments composed on the
/// imaged mesh, so that equilibrium holds on the imaged geometry rather than
/// on a virtually deformed one. Because each iteration is composed, the
/// tangent of the composed residual differs from the kernel's own, and the
/// solid kernels correct it (see correct_composed_tangent). Max_iterations and
/// Tolerance bound the iterations of a time step as in any other run. With
/// Pseudo_transient the time step is adapted to the residual between steps
/// (see Integrator.h).
///
/// Each mesh holds F0 in a Prestrain, one tensor per Gauss point,
/// written to and read from VTU cell arrays named Prestrain_F_g<g>, one per
/// Gauss point.
///
class Prestrain {
public:
  /**
   * @brief Reads the prestrain from the Prestrain_F_g<g> cell arrays of the
   * VTU file named by Prestrain_file_path, in the mesh's original element
   * order.
   */
  void read(const MeshParameters& mesh_parameters, mshType& mesh, int nsd);

  /**
   * @brief Parallel distribution of the prestrain according to the mesh's
   * element partition.
   */
  void distribute(const CmMod& cm_mod, const cmType& cm, const mshType& mesh);

  /**
   * @brief Returns true if the mesh carries a prestrain. The same on every
   * rank once distributed, also on a rank that holds none of the mesh's
   * elements.
   */
  bool has_prestrain() const { return F0.nrows() != 0; }

  /**
   * @brief Sets the prestrain to the identity at every Gauss point.
   */
  void set_identity(int nsd, int nG, int nEl);

  /**
   * @brief Returns the prestrain at Gauss point g of element e, nsd x nsd.
   * The mesh must carry a prestrain (see has_prestrain).
   */
  Eigen::MatrixXd get(int nsd, int g, int e) const;

  /**
   * @brief Returns the element average of the prestrain, for output whose
   * quadrature differs from the mesh's. The mesh must carry a prestrain.
   */
  Eigen::MatrixXd get_mean(int nsd, int e) const;

  /**
   * @brief Gathers the prestrain onto the master in the mesh's original
   * element order, for writing: one row per component and Gauss point, one
   * column per element. Collective: every rank must call it, also a rank that
   * holds none of the mesh's elements.
   */
  void gather(const ComMod& com_mod, const CmMod& cm_mod, const mshType& mesh, Array<double>& gF0) const;

  /**
   * @brief F0 <- (I + Grad u) F0 at every Gauss point of the mesh, with u a
   * nodal displacement field on the mesh, nsd x tnNo.
   */
  void compose(const ComMod& com_mod, const mshType& mesh, const Array<double>& u);

private:
  /**
   * @brief Prestrain at the quadrature points.
   *
   * Indexed as @c (component, gauss_point, element):
   * - @c component packs the nsd x nsd tensor in the column-major order of
   *   Matrix<nsd>, @f$ i + j \cdot n_{sd} @f$ for row @f$ i @f$ and column
   *   @f$ j @f$;
   * - @c element is the local element index after @ref distribute (or the
   *   global element index before it).
   */
  Array3<double> F0;

  /// The tensor of Gauss point g of element e.
  double* block(const int g, const int e) const
  {
    return F0.data() + (static_cast<size_t>(e) * F0.ncols() + g) * F0.nrows();
  }
};

namespace prestrain {

/// @brief Make a solid kernel's tangent that of a prestrain iteration.
///
/// A prestrain iteration composes the displacement u it reaches into F0 and
/// starts again from rest, so the residual the next iteration sees is that of
/// the Cauchy stress on the mesh, which does not move, with the mesh's fibers
/// pulled back anew. The kernel's tangent differs from that residual's in two
/// ways, both corrected here.
///
/// The current configuration: the kernel's tangent carries the change of the
/// current volume and of the current shape function gradients with u; this
/// subtracts that part,
///
///     w (tau Grad_x N_a (x) Grad_x N_b - tau Grad_x N_b (x) Grad_x N_a),
///
/// with tau = P F^T the Kirchhoff stress, Grad_x N = F^-T Grad_0 N, and w the
/// stress-free weight. Without it each iteration is an inexact Newton step
/// whose error grows with the stress, and the prestrain can stop converging
/// once the pseudo time step is large.
///
/// The fibers: composing the update pulls the mesh's fibers back anew (see
/// the overview above), so the fibers the material law sees change with u. To
/// first order, with G = F0^-1 Grad u F0, the gradient on the mesh carried to
/// the stress-free configuration, each fiber f0 changes by
///
///     df0 = -(I - f^ f^T) G f0,
///
/// f^ being its unit vector: it turns without stretching. The kernel's tangent
/// holds the fibers fixed; this adds
///
///     w F dS Grad_0 N_a,
///
/// taking the gradient of the material law's S with respect to each fiber by
/// one-sided finite differences along the axes and applying it to df0, so it
/// costs nsd calls of the material law per fiber, plus one for the baseline.
/// A mesh without fibers pays nothing beyond the first part.
///
/// @param active_stress_evaluator the kernel's evaluator, updated for the element.
/// @param N shape functions at the Gauss point, for the active tension.
/// @param P the first Piola-Kirchhoff stress the kernel assembled, viscous and
/// prestress parts included.
/// @param Nx0 shape function gradients on the stress-free configuration, nsd x eNoN.
/// @param fibers the pulled-back fibers the material law saw, nsd x nFn.
/// @param w_afu the weight times the stiffness scaling of lK.
/// @param lK the element tangent, lK(i*dof + j, a, b).
template <int nsd>
void correct_composed_tangent(const ComMod& com_mod, const CepMod& cep_mod, const dmnType& dmn,
    const ActiveStress::Evaluator& active_stress_evaluator, const Eigen::Ref<const Eigen::VectorXd>& N,
    const mat_fun::Matrix<nsd>& F0, const mat_fun::Matrix<nsd>& F, const mat_fun::Matrix<nsd>& P,
    const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>& Nx0,
    const Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>& fibers,   // mat_models::FiberRef<nsd>
    double w_afu, int dof, Array3<double>& lK);

/// @brief Validate the settings and the prestrains read from file, and seed
/// the identity where a prestrain run has none, so that every mesh of a
/// prestrain run carries a prestrain. Called once the meshes are partitioned.
/// Rejects a prestrain run with more than one equation, one that is not struct,
/// or a line search, and any prestrain on a ustruct solid or with remeshing.
void init(ComMod& com_mod);

};

#endif
