// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef NN_H 
#define NN_H 

#include "Simulation.h"
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

  void gnnb(const ComMod& com_mod, const faceType& lFa, const int e, const int g, const int nsd, const int insd,
      const int eNoNb, const Array<double>& Nx, Vector<double>& n, consts::MechanicalConfigurationType cfg=consts::MechanicalConfigurationType::reference);

  void gnns(const int nsd, const int eNoN, const Array<double>& Nxi, Array<double>& xl, Vector<double>& nV, 
      Array<double>& gCov, Array<double>& gCnv);

  void gn_nxx(const int l, const int eNoN, const int nsd, const int insd, Array<double>& Nxi, Array<double>& Nxi2, Array<double>& lx,
      Array<double>& Nx, Array<double>& Nxx);

  void select_ele(const ComMod& com_mod, mshType& mesh);

  void select_eleb(Simulation* simulation,  mshType& mesh, faceType& face);

};

#endif

