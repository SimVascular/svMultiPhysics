// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ALL_FUN_H 
#define ALL_FUN_H 

#include "Array3.h"
#include "SolutionStates.h"
#include "Array.h"
#include "Vector.h"
#include "ComMod.h"

#include "consts.h"

#include <optional>
#include <string>

namespace all_fun {

  double aspect_ratio(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x);

  void commu(const ComMod& com_mod, Vector<double>& u);
  void commu(const ComMod& com_mod, Array<double>& u);

  int domain(const ComMod& com_mod, const mshType& lM, const int iEq, const int e);

  void find_face(const std::vector<mshType>& mesh_list, const std::string& faceName, int& iM, int& iFa);

  void find_msh(const std::vector<mshType>& mesh_list, const std::string& mesh_name, int& iM);

  Array<double> global(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, const Array<double>& U);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, int iM, const Array<double>& s, const SolutionStates& solutions);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, int dId, const Array<double>& s, int l, int u, 
      const SolutionStates& solutions, bool pFlag=false);

  /**
   * @brief Integrate a scalar field over a boundary face.
   *
   * Reproduces 'FUNCTION IntegS(lFa, s, pflag)'. The scalar field s is
   * integrated over the face, i.e. this computes
   * \f[
   *   \int_{\Gamma} s \, d\Gamma,
   * \f]
   * where \f$\Gamma\f$ is the boundary face.
   *
   * For simulations involving structural displacement, this function allows
   * computing the integral in any of the following configurations:
   * - reference configuration (the mesh is not displaced);
   * - current configuration (the mesh is displaced by the current displacement
   *   field);
   * - old configuration (the mesh is displaced by the displacement field from
   *   previous time step).
   *
   * @param[in] com_mod The common module.
   * @param[in] cm_mod The communication module containing MPI data.
   * @param[in] lFa The boundary face over which the integral is computed.
   * @param[in] s The scalar value at each node of the mesh.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] pFlag Whether to use the Taylor-Hood function space for the
   *   pressure field.
   * @param[in] cfg The configuration in which the integral is computed
   *   (reference, old or current).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Vector<double> &s, const SolutionStates &solutions,
               bool pFlag, consts::MechanicalConfigurationType cfg,
               const unsigned int displacement_index);

  /**
   * @brief Integrate a scalar field over a boundary face.
   *
   * This is the overload to use in the general case. The other overload adds
   * the ability to integrate over a displaced configuration, which is only
   * relevant for simulations involving structural displacement. See it for the
   * meaning of the arguments.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Vector<double> &s, const SolutionStates &solutions,
               bool pFlag);

  /**
   * @brief Integrate one or more components of a field over a boundary face.
   *
   * Reproduces 'FUNCTION IntegG(lFa, s, l, u, THflag)'. Rows l to u of s are
   * integrated over the face. When they span the spatial dimensions the field
   * is dotted with the outward surface normal (i.e. a flux is computed),
   * \f[
   *   \int_{\Gamma} \sum_{i=l}^{u} s_i \, n_i \, d\Gamma,
   * \f]
   * otherwise the single row l is integrated as a scalar,
   * \f[
   *   \int_{\Gamma} s_l \, d\Gamma,
   * \f]
   * where \f$\Gamma\f$ is the boundary face and \f$\mathbf{n}\f$ its outward
   * unit normal.
   *
   * For simulations involving structural displacement, this function allows
   * computing the integral in any of the following configurations:
   * - reference configuration (the mesh is not displaced);
   * - current configuration (the mesh is displaced by the current displacement
   *   field);
   * - old configuration (the mesh is displaced by the displacement field from
   *   previous time step).
   *
   * @param[in] com_mod The common module.
   * @param[in] cm_mod The communication module containing MPI data.
   * @param[in] lFa The boundary face over which the integral is computed.
   * @param[in] s The field value at each node of the mesh.
   * @param[in] l The first row of s to integrate.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] uo The last row of s to integrate. Defaults to l, i.e. a single
   *   component.
   * @param[in] THflag Whether to use the Taylor-Hood function space for the
   *   pressure field.
   * @param[in] cfg The configuration in which the integral is computed
   *   (reference, old or current).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Array<double> &s, const int l,
               const SolutionStates &solutions, std::optional<int> uo,
               bool THflag, consts::MechanicalConfigurationType cfg,
               const unsigned int displacement_index);

  /**
   * @brief Integrate one or more components of a field over a boundary face.
   *
   * This is the overload to use in the general case. The other overload adds
   * the ability to integrate over a displaced configuration, which is only
   * relevant for simulations involving structural displacement. See it for the
   * meaning of the arguments.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Array<double> &s, const int l,
               const SolutionStates &solutions, std::optional<int> uo,
               bool THflag);

  /**
   * @brief Integrate the flux of a vector field over a boundary face.
   *
   * Reproduces 'FUNCTION IntegV(lFa, s)'. The vector field s (one component per
   * spatial dimension at each node) is dotted with the outward surface normal
   * and integrated over the face, i.e. this computes
   * \f[
   *   \int_{\Gamma} \mathbf{s} \cdot \mathbf{n} \, d\Gamma,
   * \f]
   * where \f$\Gamma\f$ is the boundary face and \f$\mathbf{n}\f$ its outward
   * unit normal.
   *
   * For simulations involving structural displacement, this function allows
   * computing the integral in any of the following configurations:
   * - reference configuration (the mesh is not displaced);
   * - current configuration (the mesh is displaced by the current displacement
   *   field);
   * - old configuration (the mesh is displaced by the displacement field from
   *   previous time step).
   *
   * @param[in] com_mod The common module.
   * @param[in] cm_mod The communication module containing MPI data.
   * @param[in] lFa The boundary face over which the integral is computed.
   * @param[in] s The vector value at each node of the mesh.
   * @param[in] solutions The solution states that the displacement fields are
   *   extracted from.
   * @param[in] cfg The configuration in which the integral is computed
   *   (reference, old or current).
   * @param[in] displacement_index The index of the displacement field in the
   *   solution arrays. This should correspond to the start index of the
   *   equation that solves for the displacement.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Array<double> &s, const SolutionStates &solutions,
               consts::MechanicalConfigurationType cfg,
               const unsigned int displacement_index);

  /**
   * @brief Integrate the flux of a vector field over a boundary face.
   *
   * This is the overload to use in the general case. The other overload adds
   * the ability to integrate over a displaced configuration, which is only
   * relevant for simulations involving structural displacement. See it for the
   * meaning of the arguments.
   */
  double integ(const ComMod &com_mod, const CmMod &cm_mod, const faceType &lFa,
               const Array<double> &s, const SolutionStates &solutions);

  bool is_domain(const ComMod& com_mod, const eqType& eq, const int node, const consts::EquationType phys);

  double jacobian(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x, const Array<double>&Nxi);

  Vector<int> local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Vector<int>& u);

  Array<double> local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Array<double>& u);

  Array3<double> local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Array3<double>& u);

  Vector<double> mkc(const ComMod& com_mod, Vector<double>& U);
  Array<double> mkc(const ComMod& com_mod, Array<double>& U);

  void mkci(const ComMod& com_mod, Vector<double>& U);
  void mkci(const ComMod& com_mod, Array<double>& U);

  void set_dmn_id(mshType& mesh, const int iDmn, const int ifirst=consts::int_inf, const int ilast=consts::int_inf);

  double skewness(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x);

  void split_jobs(int tid, int m, int n, Array<double>& A, Vector<double>& b);

  /// @brief This routine is for calculating values by the inverse of general BC
  //
  void igbc(const ComMod &com_mod, const MBType &gm, Array<double> &Y,
            Array<double> &dY);
};

#endif

