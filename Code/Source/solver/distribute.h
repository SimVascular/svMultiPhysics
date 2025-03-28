/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "Simulation.h"

#ifndef DISTRIBUTE_H
#define DISTRIBUTE_H

void distribute(Simulation* simulation);

void dist_bc(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bcType& lBc, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl);

void dist_bf(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, bfType& lBf);

void dist_ris(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm);

void dist_uris(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm);

void dist_uris_msh(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, mshType& lM, const int iUris);

void dist_eq(ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, const std::vector<mshType>& tMs,
             const Vector<int>& gmtl, CepMod& cep_mod, eqType& lEq);

void dist_mat_consts(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, stModelType& lStM);

void dist_fluid_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, fluidViscModelType& lVis);

void dist_solid_visc_model(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, solidViscModelType& lVis);

void part_face(Simulation* simulation, mshType& lM, faceType& lFa, faceType& gFa, Vector<int>& gmtl);

void part_msh(Simulation* simulation, int iM, mshType& lM, Vector<int>& mtl, int nP, Vector<float>& wgt);

#endif

