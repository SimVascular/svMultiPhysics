// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef MAT_MODELS_H 
#define MAT_MODELS_H 

#include "ActiveStress.h"
#include "Array.h"
#include "CepMod.h"
#include "ComMod.h"
#include "Tensor4.h"

#include "mat_fun.h"

#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Dense"
#include "eigen3/unsupported/Eigen/CXX11/Tensor"

namespace mat_models {

using mat_fun::Matrix;
using mat_fun::Tensor;
using mat_fun::NodalMatrix;
using mat_fun::eigen_view;

void actv_strain(const ComMod& com_mod, const CepMod& cep_mod, const double gf,
    const int nfd, const Array<double>& fl, Array<double>& Fa);

void cc_to_voigt(const int nsd, const Tensor4<double>& CC, Array<double>& Dm);

void voigt_to_cc(const int nsd, const Array<double>& Dm, Tensor4<double>& CC);

/// @brief Fiber directions at a Gauss point, one direction per column
/// with row count fixed at compile time to the number of spatial dimensions.
template <int nsd>
using FiberRef = Eigen::Ref<const Eigen::Matrix<double, nsd, Eigen::Dynamic>>;

/**
 * @brief Compute 2nd Piola-Kirchhoff stress and material stiffness tensors
 * including both dilational and isochoric components.
 *
 * @param[in] com_mod Object containing global common variables.
 * @param[in] cep_mod Object containing electrophysiology-specific common
 *   variables.
 * @param[in] lDmn Domain object.
 * @param[in] F Deformation gradient tensor.
 * @param[in] nfd Number of fiber directions.
 * @param[in] fl Fiber directions.
 * @param[in] active_tension Active tension along the fiber, sheet and
 *   sheet-normal directions, and its derivatives with respect to the fiber
 *   stretch. The derivatives contribute the tangent of the active stress to
 *   @p Dm, which is therefore not symmetric unless the active stress acts along
 *   the fiber direction alone.
 * @param[out] S 2nd Piola-Kirchhoff stress tensor (modified in place).
 * @param[out] Dm Material stiffness tensor (modified in place).
 * @param[out] Ja Jacobian for active strain
 *
 * @return None, but modifies S, Dm, and Ja in place.
 *
 * @tparam nsd Number of spatial dimensions.
 */
template <int nsd>
void compute_pk2cc(const ComMod &com_mod, const CepMod &cep_mod,
                   const dmnType &lDmn, const Matrix<nsd> &F, const int nfd,
                   const FiberRef<nsd> &fl,
                   const ActiveStress::ActiveTension &active_tension,
                   Matrix<nsd> &S, Matrix<3 * (nsd - 1)> &Dm, double &Ja);

void compute_pk2cc_shlc(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0,
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml);

void compute_pk2cc_shli(const ComMod& com_mod, const dmnType& lDmn, const int nfd, const Array<double>& fNa0,
    const Array<double>& gg_0, const Array<double>& gg_x, double& g33, Vector<double>& Sml, Array<double>& Dml);

void compute_tau(const ComMod& com_mod, const dmnType& lDmn, const double detF, const double Je, double& tauM, double& tauC);

void compute_svol_p(const ComMod& com_mod, const CepMod& cep_mod, const stModelType& stM, const double J, 
    double& p, double& pl);

void g_vol_pen(const ComMod& com_mod, const dmnType& lDmn, const double p, 
    double& ro, double& bt, double& dro, double& dbt, const double Ja);


/// @brief Viscous 2nd Piola-Kirchhoff stress and tangent contributions at a
/// Gauss point, for the viscosity model configured for the domain.
///
/// Owns the buffers it computes into, so a kernel can hold one instance across
/// an element's Gauss points and skip the calls it does not need.
///
/// @tparam nsd Number of spatial dimensions.
///
template <int nsd>
class ViscousResponse {
  public:
    /// @brief Evaluate the domain's viscosity model at this Gauss point.
    ///
    /// Does nothing when recompute is false, preserving values from the
    /// previous call. Used for linear triangules and tetrahedra for which
    /// Nx, vx, and F are constant across the element's Gauss points.
    ///
    /// @param[in] lDmn Domain, supplying the viscosity model and its parameters.
    /// @param[in] eNoN Number of element nodes.
    /// @param[in] Nx Shape function spatial derivatives.
    /// @param[in] vx Velocity gradient.
    /// @param[in] F Deformation gradient.
    /// @param[in] recompute Whether the stored values must be recomputed.
    void update(const dmnType& lDmn, const int eNoN, const Array<double>& Nx,
                const Matrix<nsd>& vx, const Matrix<nsd>& F, const bool recompute);

    /// @brief Viscous 2nd Piola-Kirchhoff stress.
    const Matrix<nsd>& S() const { return Svis_; }

    /// @brief Tangent w.r.t. displacement. du(i*nsd + j, a, b) is the (i,j)
    /// entry of the block coupling nodes a and b.
    double du(const int ij, const int a, const int b) const { return Kvis_u_(ij, a, b); }

    /// @brief Tangent w.r.t. velocity.
    double dv(const int ij, const int a, const int b) const { return Kvis_v_(ij, a, b); }

  private:
    Matrix<nsd> Svis_;
    Array3<double> Kvis_u_, Kvis_v_;
};
};

#endif

