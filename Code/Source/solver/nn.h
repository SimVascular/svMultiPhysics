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

  void get_nnx(const int nsd, const consts::ElementType eType, const int eNoN, const Array<double>& xl, 
      const Array<double>& xib, const Array<double>& Nb, const Vector<double>& xp, Vector<double>& xi, 
      Vector<double>& N, Array<double>& Nx);

  void get_xi(const int nsd, consts::ElementType eType, const int eNoN, const Array<double>& xl, const Vector<double>& xp, 
    Vector<double>& xi, bool& flag);

  void gnn(const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& x, Array<double>& Nx, 
      double& Jac, Array<double>& ks);

  /**
   * @brief Return the surface normal vector at a given element and Gauss point.
   *
   * 'g' of face 'lFa' that is the normal weighted by Jac, i.e.
   * Jac = norm(n), the Jacobian of the mapping from parent surface element to
   * reference/old/new configuration.
   *
   * For simulations involving structural displacement, this function allows
   * computing the normal vector in any of the following configurations:
   * - reference configuration (the mesh is not displaced);
   * - current configuration (the mesh is displaced by the current displacement
   *   field);
   * - old configuration (the mesh is displaced by the displacement field from
   *   previous time step).
   *
   * @todo[michelebucelli] Many of the arguments passed to this function are
   *   redundant, as they can be inferred from the others. For example, insd and
   *   eNoNb can be inferred from the size of Nx. Those arguments should be
   *   removed to simplify the function signature and make it less error-prone.
   *
   * @param[in] com_mod The common module.
   * @param[in] lFa The boundary face for which the normal vector is computed.
   * @param[in] e The face-local index of the element for which the normal
   *   vector is computed.
   * @param[in] g The Gauss point index for which the normal vector is computed.
   * @param[in] nsd The number of spatial dimensions.
   * @param[in] insd The intrinsic dimension of the boundary surface (typically
   *   nsd - 1).
   * @param[in] eNoNb Number of nodes on a boundary face element.
   * @param[in] Nx The shape function derivatives at the Gauss point.
   * @param[out] n The computed normal vector.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] cfg The configuration in which the normal vector is computed
   *   (reference, old or current).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   */
  void gnnb(const ComMod &com_mod, const faceType &lFa, const int e,
            const int g, const int nsd, const int insd, const int eNoNb,
            const Array<double> &Nx, Vector<double> &n,
            const SolutionStates &solutions,
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

