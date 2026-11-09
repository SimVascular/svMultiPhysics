// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef STRUCT_H 
#define STRUCT_H 

#include "ActiveStress.h"
#include "ComMod.h"
#include "SolutionStates.h"
#include "mat_fun.h"

namespace struct_ns {

void b_struct_2d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void b_struct_3d(const ComMod& com_mod, const int eNoN, const double w, const Vector<double>& N, 
    const Array<double>& Nx, const Array<double>& dl, const Vector<double>& hl, const Vector<double>& nV, 
    Array<double>& lR, Array3<double>& lK);

void construct_dsolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const SolutionStates& solutions);

/// @brief Momentum residual and tangent of the solid at one Gauss point.
///
/// w and Nx are the integration weight and shape function gradients on the
/// stress-free configuration, which F0 maps to the mesh, and fibers the fibers
/// pulled back to it. The density and damping, given per unit volume of the
/// mesh, are scaled by det F0. Without a prestrain F0 is the identity and the
/// stress-free configuration is the mesh (see prestrain.h).
void struct_2d(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const mat_fun::Matrix<2> &F0,
               const Array<double> &bfl,
               const Eigen::Ref<const Eigen::Matrix<double, 2, Eigen::Dynamic>> &fibers,
               const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc);

/// @brief Momentum residual and tangent of the solid at one Gauss point; see struct_2d.
void struct_3d(ComMod &com_mod, CepMod &cep_mod, const int eNoN,
               const double w, const Vector<double> &N, const Array<double> &Nx,
               const Array<double> &al, const Array<double> &yl,
               const Array<double> &dl, const mat_fun::Matrix<3> &F0,
               const Array<double> &bfl,
               const Eigen::Ref<const Eigen::Matrix<double, 3, Eigen::Dynamic>> &fibers,
               const Array<double> &pS0l,
               Vector<double> &pSl, const ActiveStress::Evaluator &active_stress_evaluator,
               Array<double> &lR, Array3<double> &lK, const bool recompute_visc);
};

#endif
