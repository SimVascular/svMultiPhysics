// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// Define prototypes for LAPACK.

extern "C" {
  extern int dgetrf_(int*, int*, double*, int*, int*, int*);
  extern int dgetri_(int*, double*, int* nd, int*, double*, int*, int*);
  extern int dgesv_(const int* N, const int* NRHS, double* A, const int* LDA, 
      int* IPIV, double* B, const int* LDB, int* INFO);
};

