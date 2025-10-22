// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#include "fils_struct.hpp"

namespace omp_la {

using namespace fsi_linear_solver;

void omp_mul_s(const int nNo, const double r, Vector<double>& U);

void omp_mul_v(const int dof, const int nNo, const double r, Array<double>& U);

void omp_sum_s(const int nNo, const double r, Vector<double>& U, const Vector<double>& V);

void omp_sum_v(const int dof, const int nNo, const double r, Array<double>& U, const Array<double>& V);

};
