// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef ALL_FUN_H 
#define ALL_FUN_H 

#include "Array3.h"
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

  double integ(const ComMod& com_mod, const CmMod& cm_mod, int iM, const Array<double>& s);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, int dId, const Array<double>& s, int l, int u, 
      bool pFlag=false);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Vector<double>& s, 
      bool pFlag=false, consts::MechanicalConfigurationType cfg=consts::MechanicalConfigurationType::reference);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& s, 
      const int l, std::optional<int> uo=std::nullopt, bool THflag=false, consts::MechanicalConfigurationType cfg=consts::MechanicalConfigurationType::reference);

  double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& s, consts::MechanicalConfigurationType cfg=consts::MechanicalConfigurationType::reference);

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

};

#endif

