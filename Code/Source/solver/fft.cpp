// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The code here replicates the code in the Fortran FFT.f file. 

#include "fft.h"
#include <math.h>

/// @brief This routine is for calculating values by the inverse of general BC
//
void igbc(const ComMod& com_mod, const MBType& gm, Array<double>& Y, Array<double>& dY)
{
  double t = fmod(com_mod.time, gm.period);
  int i = 0;

  for (int ii = 0; ii < gm.nTP - 1; ii++) {
    if (gm.t(ii+1) >= t) {
      Y = 0.0;
      dY = 0.0;
      i = ii;
      break; 
    }
  }

  double delT = gm.t(i+1) - gm.t(i);
  double tmp  = (t - gm.t(i)) / delT;

  for (int a = 0; a < gm.d.ncols(); a++) {
    for (int j = 0; j < gm.dof; j++) {
      Y(j,a) = tmp*gm.d(j,a,i+1) + gm.d(j,a,i)*(1.0-tmp);
      dY(j,a) = (gm.d(j,a,i+1) - gm.d(j,a,i)) / delT;
    }
  }
}