// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NN_H 
#define NN_H 

#include "Simulation.h"
#include "SolutionStates.h"
#include "ComMod.h"
#include "consts.h"

namespace nn {

  void get_gip(const int insd, consts::ElementType eType, const int nG, Vector<double>& w, Array<double>& xi);
  void get_gip(Simulation* simulation, faceType& face);
  void get_gip(mshType& mesh);

  void get_gnn(const int insd, consts::ElementType eType, const int eNoN, const int g, Array<double>& xi,
      Array<double>& N, Array3<double>& Nx);
  void get_gnn(const int insd, consts::ElementType eType, const int eNoN, Vector<double>& xi, Vector<double>& N, 
      Array<double>& Nx);

  void get_gnn(Simulation* simulation, int gaus_pt, faceType& face);
  void get_gnn(int gaus_pt, mshType& mesh);
  void get_gn_nxx(const int insd, const int ind2, consts::ElementType eType, const int eNoN, const int gaus_pt,
    const Array<double>& xi, Array3<double>& Nxx);

  void get_nn_bnds(const int nsd, consts::ElementType eType, const int eNoN, Array<double>& xib, Array<double>& Nb);
  void get_nn_bnds(const ComMod& com_mod, mshType& mesh);

  void get_nodal_xi(const int insd, consts::ElementType eType, const int eNoN, Array<double>& xi);

  void get_nnx(const int nsd, const consts::ElementType eType, const int eNoN, const Array<double>& xl, 
      const Array<double>& xib, const Array<double>& Nb, const Vector<double>& xp, Vector<double>& xi, 
      Vector<double>& N, Array<double>& Nx);

  void get_xi(const int nsd, consts::ElementType eType, const int eNoN, const Array<double>& xl, const Vector<double>& xp, 
    Vector<double>& xi, bool& flag);

  void gnn(const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& x, Array<double>& Nx, 
      double& Jac, Array<double>& ks);

  /**
   * @brief Return the area-weighted surface normal at a given element and Gauss
   * point.
   *
   * Returns the outward normal at element 'e', Gauss point 'g' of face 'lFa',
   * weighted by the surface Jacobian: Jac = norm(n) is the Jacobian of the
   * mapping from the parent surface element to the reference/old/new
   * configuration.
   *
   * For simulations involving structural displacement, this function allows
   * computing the normal vector in any of the following configurations:
   * - reference configuration (the mesh is not displaced);
   * - current configuration (the mesh is displaced by the current displacement
   *   field);
   * - old configuration (the mesh is displaced by the displacement field from
   *   previous time step).
   *
   * @todo[michelebucelli] This could be renamed to outward_normal, to match the
   *   function below.
   *
   * @param[in] com_mod The common module.
   * @param[in] lFa The boundary face for which the normal vector is computed.
   * @param[in] e The face-local index of the element for which the normal
   *   vector is computed.
   * @param[in] g The Gauss point index for which the normal vector is computed.
   * @param[in] Nx The shape function derivatives at the Gauss point. Its shape
   *   (insd x eNoNb) determines the surface's intrinsic dimension and the
   *   number of nodes per face element.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] cfg The configuration in which the normal vector is computed
   *   (reference, old, current or intermediate).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   *
   * @return The area-weighted outward normal vector.
   */
  Vector<double> gnnb(const ComMod &com_mod, const faceType &lFa, const int e,
                      const int g, const Array<double> &Nx,
                      const SolutionStates &solutions,
                      consts::MechanicalConfigurationType cfg,
                      const unsigned int displacement_index);

  /**
   * @brief Return the area-weighted surface normal at a given element and Gauss
   * point.
   *
   * Returns the outward normal at element 'e', Gauss point 'g' of face 'lFa',
   * weighted by the surface Jacobian: Jac = norm(n) is the Jacobian of the
   * mapping from the parent surface element to the mesh.
   *
   * This is the overload to use in the general case. The other overload adds
   * the ability to compute the normal in a displaced configuration, which is
   * only relevant for simulations involving structural displacement.
   *
   * @todo[michelebucelli] This could be renamed to outward_normal, to match the
   *   function below.
   *
   * @param[in] com_mod The common module.
   * @param[in] lFa The boundary face for which the normal vector is computed.
   * @param[in] e The face-local index of the element for which the normal
   *   vector is computed.
   * @param[in] g The Gauss point index for which the normal vector is computed.
   * @param[in] Nx The shape function derivatives at the Gauss point.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   *
   * @return The area-weighted outward normal vector.
   */
  Vector<double> gnnb(const ComMod &com_mod, const faceType &lFa, const int e,
                      const int g, const Array<double> &Nx,
                      const SolutionStates &solutions);

  /**
   * @brief Differentiate the cross product of the tangent vectors of a
   * surface element at a Gauss point with respect to the positions of the
   * element nodes.
   *
   * Tangent vector 'k' is the derivative of the position with respect to
   * parametric coordinate 'k': the sum over the element nodes of Nx(k,a) times
   * the position of node 'a'. The cross product is the one that utils::cross()
   * returns, an area-weighted normal of the element.
   *
   * @param[in] tangents The tangent vectors at the Gauss point, one per column
   *   (nsd x insd).
   * @param[in] Nx The shape function derivatives at the Gauss point
   *   (insd x eNoNb).
   *
   * @return The derivative (nsd x nsd x eNoNb): entry (j,i,a) is the derivative
   *   of component 'j' of the cross product with respect to coordinate 'i' of
   *   node 'a' of the element.
   */
  Array3<double> derivative_of_tangents_cross_product_wrt_position(
      const Array<double> &tangents, const Array<double> &Nx);

  /**
   * @brief Differentiate the area-weighted outward normal of a face element at
   * a Gauss point with respect to the positions of the element nodes.
   *
   * The normal is the one that gnnb() returns.
   *
   * Not available for the faces of shell meshes.
   *
   * @param[in] com_mod The common module.
   * @param[in] lFa The boundary face.
   * @param[in] e The face-local index of the element.
   * @param[in] Nx The shape function derivatives at the Gauss point. Its shape
   *   (insd x eNoNb) determines the surface's intrinsic dimension and the
   *   number of nodes per face element.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] cfg The configuration in which the normal is computed
   *   (reference, old, current or intermediate).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   *
   * @return The derivative (nsd x nsd x eNoNb): entry (j,i,a) is the derivative
   *   of component 'j' of the normal with respect to coordinate 'i' of node
   *   'a' of the face element.
   */
  Array3<double> derivative_of_outward_normal_wrt_position(
      const ComMod &com_mod, const faceType &lFa, const int e,
      const Array<double> &Nx, const SolutionStates &solutions,
      consts::MechanicalConfigurationType cfg,
      const unsigned int displacement_index);

  void gnns(const int nsd, const int eNoN, const Array<double>& Nxi, Array<double>& xl, Vector<double>& nV, 
      Array<double>& gCov, Array<double>& gCnv);

  void gn_nxx(const int l, const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& Nxi2, Array<double>& lx,
      Array<double>& Nx, Array<double>& Nxx);

  void select_ele(const ComMod& com_mod, mshType& mesh);

  void select_eleb(Simulation* simulation,  mshType& mesh, faceType& face);

};

#endif

