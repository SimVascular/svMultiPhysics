// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

// The code here replicates the code in the Fortran FFT.f file. 

#include "fft.h"
#include <math.h>

/// @brief Replicates Fortran 'SUBROUTINE FFT(fid, np, gt)'.
///
/// temporal_values contains all of the values read in from
/// the temporal values file.
//
void fft(const int np, const std::vector<std::vector<double>>& temporal_values, fcType& gt)
{
  using namespace consts;

  #define n_debug_fft
  #ifdef debug_fft
  DebugMsg dmsg(__func__, 0);
  dmsg.banner();
  dmsg << "np: " << np;
  #endif

  Vector<double> t(np);
  Array<double> q(gt.d, np);

  for (int i = 0; i < np; i++) {
    t[i] = temporal_values[i][0];
    for (int j = 0; j < gt.d; j++) {
      q(j,i) = temporal_values[i][j+1];
    }
  }


  gt.ti = t[0];
  gt.T = t[np-1] - t[0];
  #ifdef debug_fft
  dmsg << "pi: " << pi;
  dmsg << "gt.ti: " << gt.ti;
  dmsg << "gt.T: " << gt.T;
  dmsg << "gt.d: " << gt.d;
  dmsg << "gt.n: " << gt.n;
  dmsg << "np: " << np;
  dmsg << "gt.r.nrows(): " << gt.r.nrows();
  #endif

  for (int j = 0; j < gt.d; j++) {
    gt.qi[j] = q(j,0);
    gt.qs[j] = (q(j,np-1) - q(j,0)) / gt.T;
  }

  for (int i = 0; i < np; i++) {
    t[i] = t[i] - gt.ti;
    for (int j = 0; j < gt.d; j++) {
      q(j,i) = q(j,i) - gt.qi[j] - gt.qs[j]*t[i];
    }
  }

  int ns = 0;
  Vector<double> ns_array(512);

  for (int n = 0; n < gt.n; n++) {
    double tmp = static_cast<double>(n);
    gt.r.set_col(n, 0.0);
    gt.i.set_col(n, 0.0);

    for (int i = 0; i < np-1; i++) {
      double ko = 2.0 * pi * tmp * t[i] / gt.T;
      double kn = 2.0 * pi * tmp * t[i+1] / gt.T;

      for (int j = 0; j < gt.d; j++) {
        double s = (q(j,i+1) - q(j,i)) / (t[i+1] - t[i]);
        if (n == 0) {
          gt.r(j,n) = gt.r(j,n) + 0.5*(t[i+1] - t[i]) * (q(j,i+1) + q(j,i));
        } else {
          gt.r(j,n) = gt.r(j,n) + s*(cos(kn) - cos(ko));
          gt.i(j,n) = gt.i(j,n) - s*(sin(kn) - sin(ko));
        }
      }
    }

    if (n == 0) {
      for (int k = 0; k < gt.d; k++) {
        gt.r(k,n) = gt.r(k,n) / gt.T;
      }
    } else { 
      for (int k = 0; k < gt.d; k++) {
        gt.r(k,n) = 0.5 * gt.r(k,n) * gt.T / (pi*pi*tmp*tmp);
        gt.i(k,n) = 0.5 * gt.i(k,n) * gt.T / (pi*pi*tmp*tmp);
      }
    }
  }
}

/// @brief This is to calculate flow rate and flow acceleration (IFFT)
//
void ifft(const ComMod& com_mod, const fcType& gt, Vector<double>& Y, Vector<double>& dY) 
{
  using namespace consts;
  #define n_debug_ifft
  #ifdef debug_ifft
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  double time = com_mod.time;
  #ifdef debug_ifft
  dmsg << "time: " << time;
  dmsg << "gt.lrmp: " << gt.lrmp;
  dmsg << "gt.ti: " << gt.ti;
  dmsg << "gt.n: " << gt.n;
  dmsg << "pi: " << pi;
  #endif

  if (gt.lrmp) {
    double t = time - gt.ti;
    if (t <= 0.0) {
      t = fmax(t, 0.0);
    } else {
      t = fmin(t, gt.T);
    }
    #ifdef debug_ifft
    dmsg << "t: " <<  t;
    #endif

    for (int i = 0; i < gt.d; i++) { 
      #ifdef debug_ifft
      dmsg << " gt.qi(i): " <<  gt.qi(i);
      dmsg << " gt.qs(i): " <<  gt.qs(i);
      #endif
      Y(i) = gt.qi(i) + t*gt.qs(i);
      dY(i) = gt.qs(i);
    }

  } else {
    double t  = fmod(time - gt.ti, gt.T);
    double tmp = 2.0 * pi / gt.T;
    #ifdef debug_ifft
    dmsg << "t: " << t;
    dmsg << "tmp: " << tmp;
    dmsg << "gt.d: " << gt.d;
    #endif

    for (int j = 0; j < gt.d; j++) { 
      Y(j) = gt.qi(j) + t*gt.qs(j);
      dY(j) = gt.qs(j);
    }

    for (int i = 0; i < gt.n; i++) { 
      double dk = tmp * static_cast<double>(i);
      double K = t*dk;

      for (int j = 0; j < gt.d; j++) { 
        Y(j) = Y(j) +  gt.r(j,i)*cos(K) - gt.i(j,i)*sin(K);
        dY(j) = dY(j) - (gt.r(j,i)*sin(K) + gt.i(j,i)*cos(K))*dk;
      }
    }
  }
}

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


