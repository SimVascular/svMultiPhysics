// SPDX-FileCopyrightText: Copyright (c) Stanford University, The Regents of the University of California, and others.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef URIS_H 
#define URIS_H
 
#include "ComMod.h"
#include "SolutionStates.h"
#include "Simulation.h"

namespace uris {

void uris_meanp(ComMod& com_mod, CmMod& cm_mod, const int iUris, const SolutionStates& solutions);

void uris_meanv(ComMod& com_mod, CmMod& cm_mod, const int iUris, const SolutionStates& solutions);

void uris_update_disp(ComMod& com_mod, CmMod& cm_mod, const SolutionStates& solutions);

void uris_find_tetra(ComMod& com_mod, CmMod& cm_mod, const int iUris);

bool inside_tet(ComMod& com_mod, int& eNoN, Vector<double>& xp, Array<double>& xl, bool include_bdry);

void uris_read_msh(Simulation* simulation);

void uris_write_vtus(ComMod& com_mod);

void uris_calc_sdf(ComMod& com_mod);

void uris_read_sv(Simulation* simulation, mshType& mesh, const URISFaceParameters* mesh_param);

bool in_poly(const Vector<double>& P, const Array<double>& P1, bool include_bdry);

bool same_side(const Vector<double>& v1, const Vector<double>& v2, const Vector<double>& v3,
              const Vector<double>& v4, const Vector<double>& p, bool include_bdry);

void surface_element_barycenter(const urisType& uris_obj, int jM, int Ec, Vector<double>& xb);

void uris_find_closest_face_centroid(const urisType& uris_obj, const Vector<double>& xp,
  const int nsd, double& minS, int& Ec, int& jM, Vector<double>& xb);

void uris_face_unit_normal(const urisType& uris_obj, const int nsd, const int jM, const int Ec,
  Vector<double>& unitNormal);

double uris_compute_sdf_sign(const urisType& uris_obj, const Vector<double>& xp,
  const Vector<double>& xb, const double dotP);

void uris_build_fluid_node_mask(ComMod& com_mod);

}

#endif

