// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// A bunch of operation that benefits from OMP hyperthreading

#include "omp_la.h"

namespace omp_la {

/// @brief Reproduces 'SUBROUTINE OMPMULS (nNo, r, U)'.
//
void omp_mul_s(const int nNo, const double r, Vector<double>& U)
{
  for (int i = 0; i < nNo; i++) {
    U(i) = r * U(i);
  }
}

/// @brief Reproduces 'SUBROUTINE OMPMULV (dof, nNo, r, U)'.
//
void omp_mul_v(const int dof, const int nNo, const double r, Array<double>& U)
{
  switch (dof) {
    case 1:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r * U(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
        U(2,i) = r*U(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
        U(2,i) = r*U(2,i);
        U(3,i) = r*U(3,i);
      }
    break;

    default: 
      for (int i = 0; i < nNo; i++) {
        for (int j = 0; j < U.nrows(); j++) {
          U(j,i) = r*U(j,i);
        }
        //U(:,i) = r*U(:,i)
      }
  }
}

/// @brief Reproduces 'SUBROUTINE OMPSUMS (nNo, r, U, V)'.
//
void omp_sum_s(const int nNo, const double r, Vector<double>& U, const Vector<double>& V)
{
  for (int i = 0; i < nNo; i++) {
    U(i) = U(i) + r*V(i);
  }
}

/// @brief Reproduces 'SUBROUTINE OMPSUMV (dof, nNo, r, U, V)'.
//
void omp_sum_v(const int dof, const int nNo, const double r, Array<double>& U, const Array<double>& V)
{
  switch (dof) {

    case 1:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
        U(2,i) = U(2,i) + r*V(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
        U(2,i) = U(2,i) + r*V(2,i);
        U(3,i) = U(3,i) + r*V(3,i);
      }
    break;

    default:
      for (int j = 0; j < U.nrows(); j++) {
        for (int i = 0; i < nNo; i++) {
          U(j,i) = U(j,i) + r*V(j,i);
        }
        //U(:,i) = U(:,i) + r*V(:,i)
      } 
  }
}


};


